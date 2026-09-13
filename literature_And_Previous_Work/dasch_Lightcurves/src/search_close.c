// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 * Search for multiple images near the brightest stars.  This is an extension and adaptation of the filter_wedge.c program.
 *
 * search_close -v -s 1.786 -w 17412 -h 22026 -r i31090_00_01r180ww -i /home/scanner/Pipeline/match/i31090_00_01r180ww_tnx.tmp -o /home/scanner/Pipeline/match/i31090_00_01r180ww_close.db
 *
 * <plate>_close_gmt0.txt - 3d plot of (x, y, bin_count)
 * <plate>_close_gmt1.txt - convolution statistics
 * <plate>_close_gmt2.txt - mask table
 * <plate>_close_gmt3.txt - region file for masked points used to create the convolution table
 * <plate>_close_gmt4.txt - region file for all masked points
 * <plate>_close_plot.db  - maps of sextractor images
 *
 * Matching limiting_iso
 *    cat *_gmt1.txt | sort -u > los.db
 *    row 'selected > 0' < los.db | column Plate MAG_ISO_med MAG_ISO_rms | sort > los1.db
 *    select series, plateNumber, spatial_bin, limiting_iso from spatialbin<n> inner join photseries using(seriesId);
 *    perl ~/perl/zerofill2.perl
 *    ~/Pipeline/makestarbase.csh los1.tmp los1.db
 *    sorttable Plate < los.db > los2.db
 *    sorttable Plate < los1.db > los3.db
 *    jointable -j Plate los2.db los3.db > los4.db
 *
 * background map
 *
 *   column -i i31090_00_01r180ww_tnx.db NUMBER X_IMAGE Y_IMAGE | sorttable -n NUMBER > los1.db
 *   jointable -n -j NUMBER los1.db i31090_00_01r180ww_background.db > los2.db
 *   votable -i los2.db -o los.xml
 */

#include <math.h>
#include <time.h>

#include "table.h"

#include "mysql.h"

#include "libwcs/fitsfile.h"

#include "scandb.h"
#include "pipelineutils.h"

#define MAX_BUFFER 512
#define SEARCH_TABLE_LENGTH 2500
#define BORDER_EDGE 0.025
#define NUM_GMT_FILES 7
#define MATCH_ARRAY_ALLOC 1000000
#define MAG_ISO_RMS_FACTOR 4
#define SQUARE_BIN_OBJECTS 10
#define CLIP_RMS_FACTOR 1.0
#define EXPECTED_PEAKS 1
#define EXCESS_LIMIT 10.0
#define MAX_PIXEL_RADIUS 1000
#define MAX_CONVOLUTION_RADIUS 11
/* #define LOS_DEBUG 3   This is the desired bin number - use 3 for mc05077 */
#define USE_ANALYSIS_THRESHOLD 2.2
#define BACKGROUND_CLIP 3.0
#if 0
#define BACKGROUND_LIMIT  8.0 /* Limit without spatial bins */
#else
#define BACKGROUND_LIMIT  3.0 /* Test limit with spatial bins */
#endif
#define BACKGROUND_BIN_SIZE 100 /* Background bin size in pixels */
#define BACKGROUND_SMOOTH_RADIUS 1010 /* Background smoothing radius (Not a multiple of the bin size) */
#define MIN_BIN_COUNT 10 /* Minimum points to determine the bin count */
#define MAX_EXPANSION_BINS 10.0
/* #define BACKGROUND_PLOTS 1 */

#ifdef USE_ANALYSIS_THRESHOLD
/* Threshold = 2.0, use 8-15 */
#define COUNT_LIMIT_RMS_FACTOR 15
#else /* USE_ANALYSIS_THRESHOLD */
#define COUNT_LIMIT_RMS_FACTOR 30.0
#endif /* USE_ANALYSIS_THRESHOLD */

int maskBrightnessTable[SEARCH_TABLE_LENGTH][SEARCH_TABLE_LENGTH];

typedef struct _backgroundbin {
  double X_IMAGE;     /* X Center point */
  double Y_IMAGE;     /* Y Center point */
  double BACKGROUND;  /* Average background in bin */
  double FILTERED_BACKGROUND; /* Smoothed background */
  double expansionRadius;
  int count;         /* Count of objects in the bin */
} BACKGROUNDBIN, *PBACKGROUNDBIN;

typedef struct _inputimage {
  double BACKGROUND;
  double FILTERED_BACKGROUND;
  double FLUX_MAX;
  double THRESHOLD;
  double MAG_ISO;
  double FWHM_IMAGE;
  double X_IMAGE;
  double Y_IMAGE;
  int NUMBER; /* Sextractor reference number */
  int BFLAGS; /* flags word */
  int ISO4;
  int ISO5;
  int spatial_bin;
} INPUTIMAGE, *PINPUTIMAGE;

typedef struct _searchentry {
  int sextractorIndex;          /* Index into the sextractor table */
  double MAG_ISO;     /* Brightness of this object */
  double FLUX_MAX;
  double THRESHOLD;
  int NUMBER;
  int X_IMAGE;
  int Y_IMAGE;
} SEARCHENTRY, *PSEARCHENTRY;

typedef struct _starimage {
  double MAG_ISO;
  double FLUX_MAX;
  double THRESHOLD;
  int NUMBER;         /* Sextractor reference number */
  int MATCH_NUMBER;   /* If nonzero, matching bright star */
  int BFLAGS;          /* flags word */
  int FWHM_IMAGE;
  int X_IMAGE;
  int Y_IMAGE;
  int maskIndex; /* Matching mask index */
  int convolutionFlag; /* Image contributed to the convolution table if nozero.   */
  int searchStar;      /* Image is in the search table if nonzero */
  struct _starimage *pSextractor1;
} STARIMAGE, *PSTARIMAGE;

SEARCHENTRY searchTable[SEARCH_TABLE_LENGTH];

typedef struct _matchentry {
  double xValue;
  double yValue;
  double magsum;  /* dim + bright iso magnitude */
  double magdiff; /* dim - bright iso magnitude */
  double drad;    /* Distance from centroid */
  int brightIndex;
  int dimIndex;
  int selected;
  int binNumber;
} MATCHENTRY, *PMATCHENTRY;

typedef struct _magsum_index {
  double magsum;
  int index;
} MAGSUM_INDEX, *PMAGSUM_INDEX;

typedef struct maskorder {
  int maskIndex;
  double MAG_ISO;
} MASKORDER, *PMASKORDER;


/* Given an angle and a radius, compute the limits */
void
GetLimits(double binAngle, double radius, double *pMinX, double *pMaxX, double *pMinY, double *pMaxY)
{
  double xval = radius * cos(binAngle*DEGREES_TO_RAD);
  double yval = radius * sin(binAngle*DEGREES_TO_RAD);
  if (xval < *pMinX) {
    *pMinX = xval;
  }
  if (xval > *pMaxX) {
    *pMaxX = xval;
  }
  if (yval < *pMinY) {
    *pMinY = yval;
  }
  if (yval > *pMaxY) {
    *pMaxY = yval;
  }
#if 0
  printf("binAngle %f radius %f, xval %f yval %f\n", binAngle, radius, xval, yval);
#endif
  return;
}

int
yCompare(const void *first, const void *second)
{
  double yFirst = ((PINPUTIMAGE)first)->Y_IMAGE;
  double ySecond = ((PINPUTIMAGE)second)->Y_IMAGE;
  if (yFirst > ySecond) {
    return(1);
  } else if (yFirst < ySecond) {
    return(-1);
  } else {
    return(0);
  }
}

int
maskOrderCompare(const void *first, const void *second)
{
  double yFirst = ((PMASKORDER)first)->MAG_ISO;
  double ySecond = ((PMASKORDER)second)->MAG_ISO;
  if (yFirst > ySecond) {
    return(1);
  } else if (yFirst < ySecond) {
    return(-1);
  } else {
    return(0);
  }
}

/* This is a reverse sort with the highest counts first */
int
maskBrightnessCompare(const void *first, const void *second)
{
  int yFirst  = ((PMULTMASK)first)->maskIndex;
  int ySecond = ((PMULTMASK)second)->maskIndex;
  if (yFirst < ySecond) {
    return(1);
  } else if (yFirst > ySecond) {
    return(-1);
  } else {
    return(0);
  }
}

int
searchCompare(const void *first, const void *second)
{
  int yFirst = ((PSEARCHENTRY)first)->Y_IMAGE;
  int ySecond = ((PSEARCHENTRY)second)->Y_IMAGE;
  if (yFirst > ySecond) {
    return(1);
  } else if (yFirst < ySecond) {
    return(-1);
  } else {
    return(0);
  }
}

/* NOTE: The following is a reverse sort */
int
magsumCompare(const void *first, const void *second)
{
  double magsumFirst = ((PMAGSUM_INDEX)first)->magsum;
  double magsumSecond = ((PMAGSUM_INDEX)second)->magsum;
  if (magsumFirst < magsumSecond) {
    return(1);
  } else if (magsumFirst > magsumSecond) {
    return(-1);
  } else {
    return(0);
  }
}

/* Search the reassignment table for the lowest possible mask number.  Return this
   number and change all mask values in the reassign table to this number */
int
GetMask(int initialMask, int *reassignTable, int *reassignStack, int maskCount, char *fileroot)
{
  int minMask = initialMask;
  int nextMask;
  int stackDepth = 0;
  int index;
  reassignStack[stackDepth] = initialMask;
  stackDepth++;
  nextMask = reassignTable[initialMask];
  while (nextMask != 0) {
    if (nextMask < minMask) {
      minMask = nextMask;
    }
    reassignStack[stackDepth] = nextMask;
    stackDepth++;
    if (stackDepth > maskCount) {
      printf("ERROR: GetMask stack overflow for %s\n", fileroot);
      exit(-1);
    }
    nextMask = reassignTable[nextMask];
  }
  if (stackDepth > 1) {
    for (index = 0; index < stackDepth; index++) {
      if (reassignStack[index] == minMask) {
        reassignTable[minMask] = 0;
      } else {
        reassignTable[reassignStack[index]] = minMask;
      }
    }
  }

  return(minMask);
}

void
FindMasks(int count_limit, int *correlation_table2, int *correlation_table3, PMULTMASK *pMaskTable, int* pTotalMasks, char *fileroot)
{
  int xIndex;
  int yIndex;
  int xIndex2;
  int yIndex2;
  int correlationIndex; /* Index into correlation table (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1))) */
  int correlationIndex2; /* Index into correlation table (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1))) */
  int maskIndex = 0;
  int maskCount;
  int maskIndex1;
  int maskIndex2;
  int *reassignTable;
  int *reassignStack;
  int maxCorrelationIndex = 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1);
  int nextMask;
  int nextMask2;
  int workDone = 1;
  int finalMaskCount = 0;
  int whileIteration = 0;
  PMULTMASK maskTable;
  PMULTMASK pMask = NULL;

  /* First assign a mask to each pixel over the convolution limit */

  for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
    for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
      correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
      if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
        printf("ERROR: correlationIndex %d exceeds limits (1) %d\n", correlationIndex, maxCorrelationIndex);
        exit(-1);
      }
      if ((correlation_table2[correlationIndex]) > count_limit) {
        correlation_table3[correlationIndex] = ++maskIndex;
#if 0
        printf("Mask %d at xIndex %d yIndex %d\n", maskIndex, xIndex, yIndex);
#endif
      }

    }
  }
  maskCount = maskIndex;
#if 0
  printf("Number of masks %d for %s\n", maskIndex, fileroot);
#endif
  if (maskIndex <= 0) {
    return;
  }
  reassignTable = (int *)calloc(maskCount+1, sizeof(int));
  if (reassignTable == NULL) {
    printf("ERROR: failed to allocate reassignTable of size %d\n", maskCount+1);
  }
  reassignStack = (int *)calloc(maskCount+1, sizeof(int));
  if (reassignStack == NULL) {
    printf("ERROR: failed to allocate reassignStack of size %d\n", maskCount+1);
  }
  /* Now examine all the neighbors of a hot pixel and assign a mask number equal to the lowest hot neighbor */
  for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
    for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
      correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
      if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
        printf("ERROR: correlationIndex %d exceeds limits (2) %d\n", correlationIndex, maxCorrelationIndex);
        exit(-1);
      }
      if (correlation_table3[correlationIndex] == 0) {
        continue;
      }
      for (xIndex2 = (xIndex-1); xIndex2 <= (xIndex+1); xIndex2++) {
        if ((xIndex2 < -MAX_PIXEL_RADIUS) || (xIndex2 > MAX_PIXEL_RADIUS)) {
          continue;
        }
        for (yIndex2 = yIndex-1; yIndex2 <= yIndex+1; yIndex2++) {
          if ((yIndex2 < -MAX_PIXEL_RADIUS) || (yIndex2 > MAX_PIXEL_RADIUS)) {
            continue;
          }
          correlationIndex2 = (xIndex2+MAX_PIXEL_RADIUS) + ((yIndex2+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
          if ((correlationIndex2 < 0) || (correlationIndex2 >= maxCorrelationIndex)) {
            printf("ERROR: correlationIndex %d exceeds limits (3) %d\n", correlationIndex, maxCorrelationIndex);
            exit(-1);
          }
          if (correlationIndex == correlationIndex2) {
            continue;
          }
          maskIndex1 = GetMask(correlation_table3[correlationIndex], reassignTable, reassignStack, maskCount, fileroot);
          maskIndex2 = GetMask(correlation_table3[correlationIndex2], reassignTable, reassignStack, maskCount, fileroot);
          if (maskIndex2 == 0) {
            continue;
          }
          if (maskIndex1 == maskIndex2) {
            continue;
          } else if (maskIndex1 < maskIndex2) {
            nextMask = correlation_table3[correlationIndex2];
            correlation_table3[correlationIndex2] = maskIndex1;
            while (nextMask != 0) {
              nextMask2 = reassignTable[nextMask];
              if (nextMask2 == maskIndex1) {
                break;
              }
              reassignTable[nextMask] = maskIndex1;
              nextMask = reassignTable[nextMask2];
            }

          } else {

            nextMask = correlation_table3[correlationIndex];
            correlation_table3[correlationIndex] = maskIndex2;
            while (nextMask != 0) {
              nextMask2 = reassignTable[nextMask];
              if (nextMask2 == maskIndex2) {
                break;
              }
              reassignTable[nextMask] = maskIndex2;
              nextMask = reassignTable[nextMask2];
            }
          }
        }
      }
    }
  }

  maskTable = (PMULTMASK)calloc(maskCount, sizeof(MULTMASK));
  if (maskTable == NULL) {
    printf("ERROR: failed to allocate maskTable of size %d\n", maskCount);
  }

  /* Now make sure that the reassign table has a complete chain */
  while (workDone) {
    workDone = 0;
    whileIteration++;
    finalMaskCount = 0;
    for (maskIndex = 1; maskIndex <= maskCount; maskIndex++) {

      if ((nextMask2 = reassignTable[maskIndex]) == 0) {
        pMask = &maskTable[finalMaskCount];
        pMask->maskIndex = maskIndex;
        finalMaskCount++;
        continue;
      }
      nextMask = GetMask(maskIndex, reassignTable, reassignStack, maskCount, fileroot);

      if (nextMask != nextMask2) {
        workDone = 1;
        if (whileIteration > maskCount) {
          printf("ERROR: too many iterations in FindMasks for %s\n", fileroot);
          for (maskIndex = 1; maskIndex <= maskCount; maskIndex++) {
            printf("Iteration %d maskIndex %d reassignTable %d\n", whileIteration, maskIndex, reassignTable[maskIndex]);
          }
          exit(-1);
        }
        break;
      }

    }
  }

  /* Now get the mask characteristics */
  for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
    for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
      correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
      if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
        printf("ERROR: correlationIndex %d exceeds limits (4) %d\n", correlationIndex, maxCorrelationIndex);
        exit(-1);
      }
      if (correlation_table3[correlationIndex] == 0) {
        continue;
      }
      maskIndex1 = GetMask(correlation_table3[correlationIndex], reassignTable, reassignStack, maskCount, fileroot);
      for (maskIndex2 = 0; maskIndex2 < finalMaskCount; maskIndex2++) {
        pMask = &maskTable[maskIndex2];
        if (pMask->maskIndex == maskIndex1) {
          break;
        }
      }
      if (maskIndex2 == finalMaskCount) {
        printf("ERROR: maskIndex %d %d not found for %s\n", correlation_table3[correlationIndex], maskIndex1, fileroot);
        exit(-1);
      }
      if (pMask->maskArea == 0) {
        pMask->maskXMin = xIndex;
        pMask->maskXMax = xIndex+1;
        pMask->maskYMin = yIndex;
        pMask->maskYMax = yIndex+1;
      } else {
        if (pMask->maskXMin > xIndex) {
          pMask->maskXMin = xIndex;
        }
        if (pMask->maskXMax < xIndex+1) {
          pMask->maskXMax = xIndex+1;
        }
        if (pMask->maskYMin > yIndex) {
          pMask->maskYMin = yIndex;
        }
        if (pMask->maskYMax < yIndex+1) {
          pMask->maskYMax = yIndex+1;
        }
      }
      pMask->maskArea++;
      if (pMask->maxBinCount < correlation_table2[correlationIndex]) {
        pMask->maxBinCount = correlation_table2[correlationIndex];
      }
    }
  }

  for (maskIndex2 = 0; maskIndex2 < finalMaskCount; maskIndex2++) {
    double xValue;
    double yValue;
    pMask = &maskTable[maskIndex2];
    if (pMask->maskArea == 0) {
      printf("ERROR: Mask %d has no entries\n", pMask->maskArea);
      exit(-1);
    }
    /* Centroids */
    xValue = 0.5*(pMask->maskXMax+pMask->maskXMin);
    yValue = 0.5*(pMask->maskYMax+pMask->maskYMin);
    pMask->maskCenterDistance = sqrt(sqr(xValue)+sqr(yValue));
    if ((pMask->maskXMin <= 0) &&
        (pMask->maskXMax >= 0) &&
        (pMask->maskYMin <= 0) &&
        (pMask->maskYMax >= 0)) {
      pMask->maskCenterFlag = 1;
    }
    xValue = (pMask->maskXMax-pMask->maskXMin);
    yValue = (pMask->maskYMax-pMask->maskYMin);
    pMask->maskAreaRatio = (1.0*pMask->maskArea)/(xValue * yValue);
#if 0
    printf("Mask %d xIndex %d %d yIndex %d %d area %d ratio %.3f maxBinCount %d maskCenterFlag %d maskCenterDistance %.1f\n",
           pMask->maskIndex,
           pMask->maskXMin,
           pMask->maskXMax,
           pMask->maskYMin,
           pMask->maskYMax,
           pMask->maskArea,
           pMask->maskAreaRatio,
           pMask->maxBinCount,
           pMask->maskCenterFlag,
           pMask->maskCenterDistance);

#endif
  }
#if 0
  printf("FinalMaskCount is %d for %s\n", finalMaskCount, fileroot);
#endif

  free(reassignTable);
  free(reassignStack);
  *pTotalMasks = finalMaskCount;
  *pMaskTable = maskTable;
  return;
}

int
GetItbg(double X_IMAGE, double Y_IMAGE, int mosaicWidth, int mosaicHeight, int nxbg, int nybg)
{
  int ixbg;
  int iybg;
  int ntbg = nxbg*nybg;
  int itbg;

  ixbg = (X_IMAGE*nxbg)/(1.0*mosaicWidth);
  iybg = (Y_IMAGE*nybg)/(1.0*mosaicHeight);
  if (ixbg < 0) {
    ixbg = 0;
  }
  if (iybg < 0) {
    iybg = 0;
  }
  if (ixbg >= nxbg) {
    ixbg = nxbg-1;
  }
  if (iybg >= nybg) {
    iybg = nybg-1;
  }

  itbg =  ixbg + (nxbg*iybg);
  if (itbg >= ntbg) {
    printf("ERROR: itbg overflow\n");
    exit(-1);
  }
  return(itbg);
}

int
main(int argc, char *argv[])
{
  char *argstr;
  char cmdchar;
  int nvals;
  int errorFlag = 0;
  char fileroot[MAX_BUFFER];
  char fileroot2[MAX_BUFFER];
  char *charPtr;
  char outfile[MAX_BUFFER];
  FILE *outHandle = NULL;
  FILE *gmtHandle[NUM_GMT_FILES];
  char gmtname[NUM_GMT_FILES][MAX_BUFFER];
  char emulsionName[MAX_BUFFER];
  FILE *emulsionHandle = NULL;
  char tmpname[MAX_BUFFER];
  int gmtIndex;
  char *slashPtr;
  char *curPtr;
  int searchIndex;
  int insertIndex;
  int searchTableSize = 0;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int candidateCount = 0;
  int outputCount = 0;
  int populatedBinCount = 0;
  int psfsaturatedCount = 0;
  int maxBinCount[MAX_CONVOLUTION_RADIUS];
  double maxBinSNR[MAX_CONVOLUTION_RADIUS];
  double maxBinSNRVal = 0;
  int max_count_limit = 0;

  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  PINPUTIMAGE sextractor_table = NULL;
  PINPUTIMAGE pInputImage = NULL;
  size_t sextractor_nrecs = 0;
  size_t vectorSize;
  size_t vectorIndex;
  int sextractorIndex;
  int sextractorIndex1;
  int sextractorIndex2;
  int minSextractorIndex;
  PSTARIMAGE image_table;
  PSTARIMAGE pSextractor1 = NULL;
  PSTARIMAGE pSextractor2 = NULL;
  PSTARIMAGE pSextractor3 = NULL;
  PSTARIMAGE pSextractorBright = NULL;
  PSTARIMAGE pSextractorDim = NULL;
  int *correlation_table1 = NULL;
  int *correlation_table2 = NULL;
  int *correlation_table3 = NULL;

  int minX;
  int maxX;
  int minY;
  int maxY;

  int verbose = 0;
  int includeSaturated = 0;
  int doPlots = 0;
  int emulsionPlots = 0;
  time_t startTime;
  time_t curTime;

  char series[MAX_SERIES_STRING];
  int plateNumber;

  PMATCHENTRY pMatchArray = NULL;
  PMATCHENTRY pSelectedMatchArray = NULL;
  int matchArraySize = 0;
  int selectedCount = 0;

  PMAGSUM_INDEX pMagsumIndexList = NULL;
  int plotCount = 0;
  double arcsecPerPixel;

  double * vector = NULL;
#ifdef USE_ANALYSIS_THRESHOLD
  /* This is the cutoff for the ratio (FLUX_MAX/THRESHOLD) */
  double analysis_threshold = USE_ANALYSIS_THRESHOLD;
#else /* USE_ANALYSIS_THRESHOLD */
  double limiting_iso;
  int curClipCount;
  double MAG_ISO_med;
  double MAG_ISO_rms;
#endif /* USE_ANALYSIS_THRESHOLD */
  double count_med[MAX_CONVOLUTION_RADIUS];
  double count_rms[MAX_CONVOLUTION_RADIUS];
  int count_limit[MAX_CONVOLUTION_RADIUS];
  int quality;
  int xIndex;
  int yIndex;
  int xIndex2;
  int yIndex2;
  int correlationIndex; /* Index into correlation table (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1))) */
  int correlationIndex2; /* Index into correlation table (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1))) */
  int convolutionCount; /* Pixel count in the correlation area */
  int convolutionArea[MAX_CONVOLUTION_RADIUS+1]; /* Area in pixels of the convolution function */
  int maxCorrelationIndex = 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1);
  int convolutionRadius = 0;
  int convolutionRadius2;
  int convolutionRadius3;
  int NUMBER;
  PSEARCHENTRY pSearchEntry1;
  PMULTMASK maskTable = NULL;
  PMULTMASK pMask = NULL;
  int totalMasks = 0;
  int maskIndex;
  int maskIndex2;
  int selected = 0;
  int binning = 1;
  PMASKORDER maskOrderTable = NULL;
  PMASKORDER pMaskOrder;
  PMASKORDER pMaskOrder2;


  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  int readFlag = 0;
  int FitWCS;
  char queryString[MAX_QUERY_STRING];
  int res;
  int maskCenterIndex = 0;
  int mosaicNumber;
  int tempbinning;
  int rotation;
  int ix;
  int iy;
  int local_bin;
  int nx = 0;
  int ny = 0;
  int ixVector[X_DMAGBINS_NORMAL];
  int iyVector[Y_DMAGBINS_NORMAL];
  int localVector[X_DMAGBINS_NORMAL*Y_DMAGBINS_NORMAL];
  int backgroundStudy = 0;
  int spatial_bin;
  double *backgroundVector[MAX_SPATIAL_BINS+1];
  double *curVector;
  int sextractorBinIndex[MAX_SPATIAL_BINS+1];
  double background_med[MAX_SPATIAL_BINS+1];
  double background_rms[MAX_SPATIAL_BINS+1];
  double total_background_med = 0;
  double total_background_rms = 0;
  int background_clip[MAX_SPATIAL_BINS+1];
  int total_background_clip = 0;
  double background_ratio;
  char backgroundfile[MAX_BUFFER];
  FILE *backgroundHandle = NULL;
#ifdef BACKGROUND_PLOTS
  char regionfile[MAX_BUFFER];
  char backgroundmap[MAX_BUFFER];
  char imagemap[MAX_BUFFER];
  FILE *regionHandle = NULL;
  FILE *backgroundmapHandle = NULL;
  FILE *imagemapHandle = NULL;
#endif /* BACKGROUND_PLOTS */
  int nxbg; /* Total x background bins */
  int nybg; /* Total y background bins */
  int ntbg; /* Total backround bins */
  int ixbg; /* x background index */
  int iybg; /* y background index */
  int itbg; /* Background bin index itgb = ixbg + (nxbg*iybg) */
  int itbg2;
  PBACKGROUNDBIN pBinMap = NULL;
  PBACKGROUNDBIN pBackBin;
  PBACKGROUNDBIN pBackBin2;
  double curExpansionRadius;

  double baseExpansionRadius;
  double maxExpansionRadius;
  int badBinCount = 0;
  int curCount = 0;
  double curBackground;
  double edgeDist;

  outfile[0] = 0;
  backgroundfile[0] = 0;
#ifdef BACKGROUND_PLOTS
  regionfile[0] = 0;
#endif /* BACKGROUND_PLOTS */

  time(&startTime);

  memset(backgroundVector, 0, sizeof(backgroundVector));
  memset(sextractorBinIndex, 0, sizeof(sextractorBinIndex));
  memset(ixVector, 0, sizeof(ixVector));
  memset(iyVector, 0, sizeof(iyVector));
  memset(localVector, 0, sizeof(localVector));
  memset(background_med, 0, sizeof(background_med));
  memset(background_rms, 0, sizeof(background_rms));
  memset(background_clip, 0, sizeof(background_clip));
  memset(maxBinCount, 0, sizeof(maxBinCount));
  memset(count_med, 0, sizeof(count_med));
  memset(count_rms, 0, sizeof(count_rms));
  memset(count_limit, 0, sizeof(count_limit));
  memset(maxBinSNR, 0, sizeof(maxBinSNR));
  memset(convolutionArea, 0, sizeof(convolutionArea));
  memset(gmtHandle, 0, sizeof(gmtHandle));

  convolutionArea[0] = 1;

#ifdef LOS_DEBUG
  printf("ERROR: LOS_DEBUG is set \n");
#endif /* LOS_DEBUG */

  sextractor_name[0] = 0;
  fileroot[0] = 0;
  fileroot2[0] = 0;
  arcsecPerPixel = 0.0;
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      printf("ERROR: argument %s does not have a qualifier\n", argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile, *++argv, MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n", *argv);
            }
          }
          break;

        case 'i': /* sextractor file name */
        case 'I':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sextractor_name, *++argv, MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n", *argv);
            }
          }
          break;

        case 'r': /* file root */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(fileroot, *++argv, MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n", *argv);
            }
            strcpy(fileroot2, fileroot);
            charPtr = strstr(fileroot2, "_");
            if (charPtr != NULL) {
              *charPtr = 0;
            }
          }
          break;

        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv, "%d", &mosaicWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic width %s\n", *argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* mosaic height */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv, "%d", &mosaicHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic height %s\n", *argv);
              errorFlag = 1;
            }
          }
          break;

        case 'b': /* binning */
        case 'B':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv, "%d", &binning);
            if (nvals != 1) {
              printf("ERROR: Unable to decode binning %s\n", *argv);
              errorFlag = 1;
            }
          }
          break;


        case 's': /* arcsec/pixel */
        case 'S':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv, "%lf", &arcsecPerPixel);
            if (nvals != 1) {
              printf("ERROR: Can not decode number of arcsec per pixel\n");
              errorFlag = 1;
            }
          }
          break;
#ifdef USE_ANALYSIS_THRESHOLD
        case 't': /* analysis_threshold */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv, "%lf", &analysis_threshold);
            if (nvals != 1) {
              printf("ERROR: Can not decode analysis threshold\n");
              errorFlag = 1;
            }
          }
          break;
#endif /* USE_ANALYSIS_THRESHOLD */

        case 'u': /* include saturated PSF values */
        case 'U':
          includeSaturated = 1;
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'd': /* search for calibration squares */
        case 'D':
          backgroundStudy = 1;
          argc--;
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(backgroundfile, *++argv, MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n", *argv);
            }
#ifdef BACKGROUND_PLOTS
            strncpy(regionfile, *++argv, MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n", *argv);
            }
#else /* BACKGROUND_PLOTS */
            ++argv;
#endif /* BACKGROUND_PLOTS */

          }
          break;

          break;

        case 'p': /* do plotting */
        case 'P':
          doPlots = 1;
          break;

        case 'q': /* do emulsion plots */
        case 'Q':
          emulsionPlots = 1;
          break;

        default:
          printf("ERROR:  unknown command -%c\n", cmdchar);
          errorFlag = 1;

        }
      }
    }
  }

  if (mosaicWidth == 0) {
    printf("ERROR: No mosaic width specified\n");
    errorFlag = 1;
  }

  if (mosaicHeight == 0) {
    printf("ERROR: No mosaic height specified\n");
    errorFlag = 1;
  }

  if (arcsecPerPixel == 0.0) {
    printf("ERROR: No arcsec per pixel plate scale was specified\n");
    errorFlag = 1;
  }

  if (outfile[0] == 0) {
    printf("ERROR: No output filename was specified\n");
    errorFlag = 1;
  }

  if (sextractor_name[0] == 0) {
    printf("ERROR: No sextractor filename was specified\n");
    errorFlag = 1;
  }

  if (fileroot[0] == 0) {
    printf("ERROR: No file root was specified\n");
    errorFlag = 1;
  }

  if (ParseFilename(sextractor_name, series, &plateNumber, &mosaicNumber, &tempbinning, &rotation) == 0) {
    printf("ERROR: Failed to parse filename %s\n", sextractor_name);
    errorFlag = 1;
  }

  sextractor_handle = Open(sextractor_name, "r");
  if (sextractor_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the sextractor file %s\n", sextractor_name);
  } else {
    if (verbose) {
      printf("Found sextractor file %s\n", sextractor_name);
    }
  }

  if (doPlots) {
    for (gmtIndex = 0; gmtIndex < NUM_GMT_FILES; gmtIndex++) {
      strcpy(gmtname[gmtIndex], outfile);
      curPtr = gmtname[gmtIndex];
      while ((slashPtr = strstr(curPtr, "/")) != NULL) {
        curPtr = slashPtr+1;
      }
      if (curPtr != NULL) {
        *curPtr = 0;
      }
      strcat(gmtname[gmtIndex], fileroot);
      sprintf(tmpname, "_close_gmt%d.txt", gmtIndex);
      strcat(gmtname[gmtIndex], tmpname);

    }
  }

  if (backgroundStudy) {
    if (backgroundfile[0] == 0) {
      printf("ERROR: No background filename specified\n");
      errorFlag = 1;
    }
#ifdef BACKGROUND_PLOTS
    if (regionfile[0] == 0) {
      printf("ERROR: No background filename specified\n");
      errorFlag = 1;
    }
#endif /* BACKGROUND_PLOTS */
  }

  if (errorFlag) {
    printf("Usage: search_close\n");
    printf("                 -b  <binning> Set to 16 for regions on thumbnail mosaics\n");
    printf("                 -d  <background file> <region file>search for calibration squares\n");
    printf("                 -h  <mosaic height>\n");
    printf("                 -i <sextractor file> \n");
    printf("                 -o  <output file>\n");
    printf("                 -p  write GMT plotting files\n");
    printf("                 -q  write emulsion plotting files\n");
    printf("                 -r <file root>\n");
    printf("                 -s <arcsec per pixel>\n");
#ifdef USE_ANALYSIS_THRESHOLD
    printf("                 -t <analysis_threshold>, cutoff ratio of FLUX_MAX/THRESHOLD\n");
#endif /* USE_ANALYSIS_THRESHOLD */

    printf("                 -u  include saturated PSFs\n");
    printf("                 -v  verbose\n");
    printf("                 -w  <mosaic width>\n");
    return(-1);
  }

  nx = XDmagBins(mosaicWidth, mosaicHeight);
  ny = YDmagBins(mosaicWidth, mosaicHeight);

  if (verbose) {
    printf("search_close of %s %s \n Input Filename %s\n Output Filename %s mosaic width %d mosaic height %d scale %f arcsec/pixel",
           __DATE__, __TIME__, sextractor_name, outfile,
           mosaicWidth,
           mosaicHeight,
           arcsecPerPixel);
#ifdef USE_ANALYSIS_THRESHOLD
    printf(" threshold: %f", analysis_threshold);
#endif /* USE_ANALYSIS_THRESHOLD */
    printf("\n");
  }

  dasch_init_scandb(pConnection);

  /* Clear the multiple mask status */
  SetMosaicFitWCS(pConnection, "NoMultipleMask", fileroot, 0, 0, readFlag, &FitWCS);
  SetMosaicFitWCS(pConnection, "HaveMultipleMask", fileroot, 0, 0, readFlag, &FitWCS);

  sprintf(queryString, "DELETE FROM masks where series = '%s' and plateNumber = %d and mosaicNumber = %d;",
          series,
          plateNumber,
          mosaicNumber);
  if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
    printf("ERROR: search_close MAX_QUERY_STRING exceeded\n");
    exit(-1);
  }

  res = ExecuteQuery(pConnection, queryString);
  if (res) {
    exit(-1);
  }

  sextractor_header = table_header(sextractor_handle, TABLE_PARSE);
  if (sextractor_header == NULL) {
    printf("ERROR: Failed to read header for %s\n", sextractor_name);
    return(-1);
  }

  sextractor_table = table_loadva(sextractor_handle,
                                  &sextractor_header,
                                  NULL, /* hbase */
                                  NULL, /* rows */
                                  NULL,
                                  sizeof(INPUTIMAGE),
                                  &sextractor_nrecs,
                                  TblInt, "NUMBER", TblOff(PINPUTIMAGE, NUMBER),
                                  TblInt, "BFLAGS", TblOff(PINPUTIMAGE, BFLAGS),
                                  TblInt, "ISO5", TblOff(PINPUTIMAGE, ISO5),
                                  TblInt, "ISO4", TblOff(PINPUTIMAGE, ISO4),
                                  TblDbl, "MAG_ISO", TblOff(PINPUTIMAGE, MAG_ISO),
                                  TblDbl, "BACKGROUND", TblOff(PINPUTIMAGE, BACKGROUND),
                                  TblDbl, "FLUX_MAX", TblOff(PINPUTIMAGE, FLUX_MAX),
                                  TblDbl, "THRESHOLD", TblOff(PINPUTIMAGE, THRESHOLD),
                                  TblDbl, "FWHM_IMAGE", TblOff(PINPUTIMAGE, FWHM_IMAGE),
                                  TblDbl, "X_IMAGE", TblOff(PINPUTIMAGE, X_IMAGE),
                                  TblDbl, "Y_IMAGE", TblOff(PINPUTIMAGE, Y_IMAGE),
                                  0, "end", 0);

  if (sextractor_table == NULL) {
    printf("ERROR: Failed to read table for %s\n", sextractor_name);
    return(-1);
  }

  if (verbose) {
    printf("read %zu records for %s\n", sextractor_nrecs, sextractor_name);
  }

  /* Calculate the spatial bin for each record */
  for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
    pInputImage = &sextractor_table[sextractorIndex];
    pInputImage->spatial_bin =  CalculateBin(mosaicWidth, mosaicHeight, pInputImage->X_IMAGE, pInputImage->Y_IMAGE, &edgeDist);
  }

  if (backgroundStudy) {
    nxbg = (mosaicWidth+BACKGROUND_BIN_SIZE-1)/BACKGROUND_BIN_SIZE;
    nybg = (mosaicHeight+BACKGROUND_BIN_SIZE-1)/BACKGROUND_BIN_SIZE;
    ntbg = nxbg*nybg;
    pBinMap = (PBACKGROUNDBIN)calloc(ntbg, sizeof(BACKGROUNDBIN));
    if (pBinMap == NULL) {
      printf("ERROR: failed to allocate the background map\n");
      exit(-1);
    }
    for (itbg = 0; itbg < ntbg; itbg++) {
      ixbg = itbg % nxbg;
      iybg = itbg / nxbg;
      pBackBin = &pBinMap[itbg];
      pBackBin->X_IMAGE = ((1.0*mosaicWidth)*(0.5+(1.0*ixbg)))/(1.0*nxbg);
      pBackBin->Y_IMAGE = ((1.0*mosaicHeight)*(0.5+(1.0*iybg)))/(1.0*nybg);
      pBackBin->count = 0;
    }

    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      backgroundVector[spatial_bin] = (double *)calloc(sextractor_nrecs, sizeof(double));
      if (backgroundVector[spatial_bin] == NULL) {
        printf("ERROR Failed to allocate the background backgroundVector\n");
        exit(-1);
      }
    }
    /* Compute the average background in each bin */
#ifdef BACKGROUND_PLOTS
    strcpy(backgroundmap, regionfile);
    slashPtr = strrchr(backgroundmap, '/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = backgroundmap;
    }
    strcpy(slashPtr, "map.db");
    backgroundmapHandle = fopen(backgroundmap, "wt");
    if (backgroundmapHandle == NULL) {
      printf("ERROR: failed to open %s\n", backgroundmap);
      exit(-1);
    }
    fprintf(backgroundmapHandle, "X_IMAGE\tY_IMAGE\tBACKGROUND\tFILTERED_BACKGROUND\tcount\texpansionRadius\n");
    fprintf(backgroundmapHandle, "-------\t-------\t----------\t-----------------\t-----\t---------------\n");
#endif /* BACKGROUND_PLOTS */
    for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
      pInputImage = &sextractor_table[sextractorIndex];
      itbg = GetItbg(pInputImage->X_IMAGE, pInputImage->Y_IMAGE, mosaicWidth, mosaicHeight, nxbg, nybg);

      pBackBin = &pBinMap[itbg];
      pBackBin->BACKGROUND += pInputImage->BACKGROUND;
      pBackBin->count++;
    }
    if (mosaicWidth < mosaicHeight) {
      baseExpansionRadius = (1.0*mosaicWidth)/(1.0*nxbg);
    } else {
      baseExpansionRadius = (1.0*mosaicHeight)/(1.0*nybg);
    }
    maxExpansionRadius = MAX_EXPANSION_BINS * baseExpansionRadius;
    for (itbg = 0; itbg < ntbg; itbg++) {
      pBackBin = &pBinMap[itbg];
      curExpansionRadius = baseExpansionRadius/2.0;
      while (curExpansionRadius < maxExpansionRadius) {
        curCount = 0;
        curBackground = 0.0;
        if (curExpansionRadius < baseExpansionRadius) {
          curCount += pBackBin->count;
          curBackground = pBackBin->BACKGROUND;
        } else {
          for (ixbg = 0; ixbg < nxbg; ixbg++) {
            pBackBin2 = &pBinMap[ixbg];
            if (fabs(pBackBin->X_IMAGE-pBackBin2->X_IMAGE) >= curExpansionRadius) {
              continue;
            }
            for (iybg = 0; iybg < nybg; iybg++) {
              itbg2 =  ixbg + (nxbg*iybg);
              pBackBin2 = &pBinMap[itbg2];
              if ((pBackBin2->Y_IMAGE-pBackBin->Y_IMAGE) > curExpansionRadius) {
                break;
              }
              if (sqrt(sqr(pBackBin->X_IMAGE-pBackBin2->X_IMAGE)+sqr(pBackBin->Y_IMAGE-pBackBin2->Y_IMAGE)) < curExpansionRadius) {
                curCount += pBackBin2->count;
                curBackground += pBackBin2->BACKGROUND;
              }
            }

          }
        }

        /* Temporarily use the smoothed background as a placeholder */

        if (curCount >= MIN_BIN_COUNT) {
          pBackBin->FILTERED_BACKGROUND = curBackground/(1.0*curCount);
          pBackBin->expansionRadius = curExpansionRadius;
#ifdef BACKGROUND_PLOTS
          fprintf(backgroundmapHandle, "%.1f\t%.1f\t%.1f\t%.1f\t%d\t%.1f\n", pBackBin->X_IMAGE, pBackBin->Y_IMAGE, pBackBin->BACKGROUND, pBackBin->FILTERED_BACKGROUND, curCount, pBackBin->expansionRadius);
#endif /* BACKGROUND_PLOTS */
          break;
        }
        curExpansionRadius += baseExpansionRadius;
      }
      if (curCount < MIN_BIN_COUNT) {
        badBinCount++;
      }
    }

    for (itbg = 0; itbg < ntbg; itbg++) {
      pBackBin = &pBinMap[itbg];
      pBackBin->BACKGROUND = pBackBin->FILTERED_BACKGROUND;
      pBackBin->FILTERED_BACKGROUND = 0;
      pBackBin->count = 0;
    }

    /* Now compute the smoothed average background in each bin */
    for (itbg = 0; itbg < ntbg; itbg++) {
      pBackBin = &pBinMap[itbg];
      for (ixbg = 0; ixbg < nxbg; ixbg++) {
        pBackBin2 = &pBinMap[ixbg];
        if (fabs(pBackBin->X_IMAGE-pBackBin2->X_IMAGE) >= (1.0*BACKGROUND_SMOOTH_RADIUS)) {
          continue;
        }
        for (iybg = 0; iybg < nybg; iybg++) {
          itbg2 =  ixbg + (nxbg*iybg);
          if (itbg2 == itbg) {
            continue;
          }
          pBackBin2 = &pBinMap[itbg2];
          if ((pBackBin2->Y_IMAGE-pBackBin->Y_IMAGE) > (1.0*BACKGROUND_SMOOTH_RADIUS)) {
            break;
          }
          if ((sqr(pBackBin2->X_IMAGE-pBackBin->X_IMAGE) + sqr(pBackBin2->Y_IMAGE - pBackBin->Y_IMAGE)) < (1.0*BACKGROUND_SMOOTH_RADIUS*BACKGROUND_SMOOTH_RADIUS)) {
            pBackBin->count++;
            pBackBin->FILTERED_BACKGROUND += pBackBin2->BACKGROUND;
          }

        }
      }
    }

    for (itbg = 0; itbg < ntbg; itbg++) {
      pBackBin = &pBinMap[itbg];
      if (pBackBin->count > 0) {
        pBackBin->FILTERED_BACKGROUND = pBackBin->FILTERED_BACKGROUND/(1.0*pBackBin->count);
      }
#if 0
#ifdef BACKGROUND_PLOTS

      fprintf(backgroundmapHandle, "%.1f\t%.1f\t%.1f\t%.1f\t%d\t%.1f\n", pBackBin->X_IMAGE, pBackBin->Y_IMAGE, pBackBin->BACKGROUND, pBackBin->FILTERED_BACKGROUND, pBackBin->count, pBackBin->expansionRadius);
#endif /* BACKGROUND_PLOTS */
#endif
    }

#ifdef BACKGROUND_PLOTS
    strcpy(imagemap, regionfile);
    slashPtr = strrchr(imagemap, '/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = imagemap;
    }
    strcpy(slashPtr, "image.db");
    imagemapHandle = fopen(imagemap, "wt");
    if (imagemapHandle == NULL) {
      printf("ERROR: failed to open %s\n", imagemap);
      exit(-1);
    }
    fprintf(imagemapHandle, "X_IMAGE\tY_IMAGE\tBACKGROUND\tFILTERED_BACKGROUND\tcount\tbackground_ratio\n");
    fprintf(imagemapHandle, "-------\t-------\t----------\t-----------------\t-----\t----------\n");
#endif /* BACKGROUND_PLOTS */

    for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
      pInputImage = &sextractor_table[sextractorIndex];
      itbg = GetItbg(pInputImage->X_IMAGE, pInputImage->Y_IMAGE, mosaicWidth, mosaicHeight, nxbg, nybg);
      pBackBin = &pBinMap[itbg];
      pInputImage->FILTERED_BACKGROUND = pInputImage->BACKGROUND - pBackBin->FILTERED_BACKGROUND;
      curVector = backgroundVector[pInputImage->spatial_bin];
      curVector[sextractorBinIndex[pInputImage->spatial_bin]] = pInputImage->FILTERED_BACKGROUND;
      sextractorBinIndex[pInputImage->spatial_bin]++;
#if 0
#ifdef BACKGROUND_PLOTS
      fprintf(imagemapHandle, "%.1f\t%.1f\t%.1f\t%.1f\t%d\n", pInputImage->X_IMAGE, pInputImage->Y_IMAGE, pInputImage->BACKGROUND, pInputImage->FILTERED_BACKGROUND, pBackBin->count);
#endif /* BACKGROUND_PLOTS */
#endif
    }

    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      curVector = backgroundVector[pInputImage->spatial_bin];
      background_clip[spatial_bin] = CalcMedianAndRMS(sextractorBinIndex[spatial_bin], 10, curVector, &background_med[spatial_bin], &background_rms[spatial_bin], 1, BACKGROUND_CLIP, 0);
      total_background_clip += background_clip[spatial_bin];
      if (background_clip[spatial_bin] == 0) {
        printf("ERRORD: Insufficient background stars in spatial bin %d for %s\n", spatial_bin, sextractor_name);
        background_rms[spatial_bin] = 0;
      } else if (background_rms[spatial_bin] == 0) {
        printf("ERRORD: Zero background rms in spatial bin %d for %s\n", spatial_bin, sextractor_name);
      } else {
        total_background_med += background_med[spatial_bin];
        total_background_rms += background_rms[spatial_bin];
      }
    }

    int clip_count = 0;

    if (total_background_clip > 0) {
      backgroundHandle = fopen(backgroundfile, "wt");
      if (backgroundHandle == NULL) {
        fprintf(stderr, "ERROR: Failed to open the output file %s\n", backgroundfile);
        exit(-1);
      } else {
        fprintf(backgroundHandle, "NUMBER\tAFLAGS\tbackground_ratio\n");
        fprintf(backgroundHandle, "------\t------\t----------------\n");
      }
#ifdef BACKGROUND_PLOTS
      regionHandle = fopen(regionfile, "wt");
      if (regionHandle == NULL) {
        fprintf(stderr, "ERROR: Failed to open the output file %s\n", regionfile);
        exit(-1);
      }
#endif /* BACKGROUND_PLOTS */

      for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
        pInputImage = &sextractor_table[sextractorIndex];
        itbg = GetItbg(pInputImage->X_IMAGE, pInputImage->Y_IMAGE, mosaicWidth, mosaicHeight, nxbg, nybg);
        pBackBin = &pBinMap[itbg];
        if (background_rms[pInputImage->spatial_bin] > 0) {

          background_ratio = (pInputImage->FILTERED_BACKGROUND-background_med[pInputImage->spatial_bin])/background_rms[pInputImage->spatial_bin];
#ifdef BACKGROUND_PLOTS
          fprintf(imagemapHandle, "%.1f\t%.1f\t%.1f\t%.1f\t%d\t%.1f\n", pInputImage->X_IMAGE, pInputImage->Y_IMAGE, pInputImage->BACKGROUND, pInputImage->FILTERED_BACKGROUND, pBackBin->count, background_ratio);
#endif /* BACKGROUND_PLOTS */

          if (pInputImage->FILTERED_BACKGROUND > (background_med[pInputImage->spatial_bin] + (BACKGROUND_CLIP*background_rms[pInputImage->spatial_bin]))) {
            if (((pInputImage->FILTERED_BACKGROUND-background_med[pInputImage->spatial_bin])/background_rms[pInputImage->spatial_bin]) > BACKGROUND_LIMIT) {
              fprintf(backgroundHandle, "%d\t%d\t%.2f\n", pInputImage->NUMBER, (1<<FILTER_AFLAG_BACKGROUND), background_ratio);
#ifdef BACKGROUND_PLOTS
              fprintf(regionHandle, "CIRCLE(%f, %f, 5) # text = {%d}\n", pInputImage->X_IMAGE/binning, pInputImage->Y_IMAGE/binning, background_ratio);
#endif /* BACKGROUND_PLOTS */
              clip_count++;
            } else {
#ifdef BACKGROUND_PLOTS
              fprintf(regionHandle, "CIRCLE(%f, %f, 2) # color = blue\n", pInputImage->X_IMAGE/binning, pInputImage->Y_IMAGE/binning);
#endif /* BACKGROUND_PLOTS */
            }
          } else {
#ifdef BACKGROUND_PLOTS
            fprintf(regionHandle, "CIRCLE(%f, %f, 2) # color = red\n", pInputImage->X_IMAGE/binning, pInputImage->Y_IMAGE/binning);
#endif /* BACKGROUND_PLOTS */
          }
        }
      }

      printf("background median %8.1f, background rms %6.1f clipped objects %5zu plotted objects %6d bad bins %5d of %5d for %s%05d_%02d\n", total_background_med/MAX_SPATIAL_BINS, total_background_rms/MAX_SPATIAL_BINS, sextractor_nrecs-total_background_clip, clip_count, badBinCount, ntbg, series, plateNumber, mosaicNumber);
    }
  }

  /* Sort the table in increasing Y_IMAGE values */
  qsort(sextractor_table, sextractor_nrecs, sizeof(INPUTIMAGE), yCompare);
  /* Now copy the table over, converting X_IMAGE and Y_IMAGE to integer pixels and setting the FILTER_BFLAG_PSFSATURATED flag */
  image_table = (PSTARIMAGE)calloc(sextractor_nrecs, sizeof(STARIMAGE));
  if (image_table == NULL) {
    printf("ERROR: failed to allocate image_table of size %zu\n", sextractor_nrecs);
    exit(-1);
  }

  for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
    pInputImage = &sextractor_table[sextractorIndex];
    pSextractor3 = &image_table[sextractorIndex];

    if (emulsionPlots) {
      /* Reject saturated and bad images,
         images that are too large,
         and image below the analysis threshold
      */
      if (((pInputImage->BFLAGS & (~((1 << FILTER_BFLAG_NEIGHBORS) | (1 << FILTER_BFLAG_BLEND)))) == 0) &&
          (pInputImage->FWHM_IMAGE <= MAX_PIXEL_RADIUS) &&
          ((pInputImage->FLUX_MAX/pInputImage->THRESHOLD) >= analysis_threshold)) {
        ix = (pInputImage->X_IMAGE * nx) /(1.0 * mosaicWidth);
        iy = (pInputImage->Y_IMAGE * ny) /(1.0 * mosaicHeight);
        local_bin = ix + (nx*iy);
        if ((ix >= 0) && (ix < nx)) {
          ixVector[ix]++;
        }
        if ((iy >= 0) && (iy < ny)) {
          iyVector[iy]++;
        }
        if ((local_bin >= 0) && (local_bin < (nx*ny))) {
          localVector[local_bin]++;
        }
      } else {
#if 0
        printf("NUMBER %d BFLAGS %d %x, FWHM_IMAGE %f, FLUX_MAX/THRESHOLD %f\n",
               pInputImage->NUMBER,
               pInputImage->BFLAGS,
               pInputImage->BFLAGS,
               pInputImage->FWHM_IMAGE,
               pInputImage->FLUX_MAX/pInputImage->THRESHOLD);
#endif

      }

    }

    pSextractor3->NUMBER = pInputImage->NUMBER;
    pSextractor3->BFLAGS = pInputImage->BFLAGS & FILTER_BMASK_SEXTRACTOR;
    pSextractor3->MAG_ISO = pInputImage->MAG_ISO;
    pSextractor3->FLUX_MAX = pInputImage->FLUX_MAX;
    pSextractor3->THRESHOLD = pInputImage->THRESHOLD;
    pSextractor3->FWHM_IMAGE = pInputImage->FWHM_IMAGE;
    pSextractor3->X_IMAGE = pInputImage->X_IMAGE+0.5;
    pSextractor3->Y_IMAGE = pInputImage->Y_IMAGE+0.5;

    /* Set the FILTER_BFLAG_PSFSATURATED flag */
    if (includeSaturated == 0) {
      if ((pInputImage->ISO4 > 0) &&
          (((1.0*pInputImage->ISO5)/(1.0*pInputImage->ISO4)) > PSFSATURATED_ISO) &&
          ((pInputImage->BACKGROUND + pInputImage->FLUX_MAX) > PSFSATURATED_FLUX)) {
#if 0
        printf("%d %f %f\n", pInputImage->NUMBER, (1.0*pInputImage->ISO5)/(1.0*pInputImage->ISO4), pInputImage->BACKGROUND + pInputImage->FLUX_MAX);
#endif
        pSextractor3->BFLAGS |= (1 << FILTER_BFLAG_PSFSATURATED);
        psfsaturatedCount++;
      }
    }
  }

  if (emulsionPlots != 0) {
    strcpy(emulsionName, outfile);
    curPtr = emulsionName;
    while ((slashPtr = strstr(curPtr, "/")) != NULL) {
      curPtr = slashPtr+1;
    }
    if (curPtr != NULL) {
      *curPtr = 0;
    }
    strcat(emulsionName, fileroot);
    sprintf(tmpname, "_close_plot.db");
    strcat(emulsionName, tmpname);
    emulsionHandle = fopen(emulsionName, "wt");

    if (emulsionHandle == NULL) {
      printf("ERROR: failed to open the emulsion plot name %s\n", emulsionName);
    } else {
      fprintf(emulsionHandle, "ix\tiy\tcount\tseries\tplateNumber\tmosaicNumber\n");
      fprintf(emulsionHandle, "--\t--\t-----\t------\t-----------\t------------\n");
      for (ix = 0; ix < nx; ix++) {
        fprintf(emulsionHandle, "%d\t-1\t%d\t%s\t%d\t%d\n", ix, ixVector[ix], series, plateNumber, mosaicNumber);
        for (iy = 0; iy < ny; iy++) {
          if (ix == 0) {
            fprintf(emulsionHandle, "-1\t%d\t%d\t%s\t%d\t%d\n", iy, iyVector[iy], series, plateNumber, mosaicNumber);
          }
          local_bin = ix + (nx*iy);
          fprintf(emulsionHandle, "%d\t%d\t%d\t%s\t%d\t%d\n", ix, iy, localVector[local_bin], series, plateNumber, mosaicNumber);
        }
      }
      fclose(emulsionHandle);
    }
  }

  /* At this point, we are done with the sextractor_table.  Free memory */
  if (sextractor_table != NULL) {
    Free(sextractor_table);
    sextractor_table = NULL;
  }

  if (sextractor_header != NULL) {
    table_hdrfree(sextractor_header);
    sextractor_header = NULL;
  }

  if (sextractor_handle != NULL) {
    Close(sextractor_handle);
    sextractor_handle = NULL;
  }

#if 0
  matchArrayAlloc = MATCH_ARRAY_ALLOC;
  pMatchArray = (PMATCHENTRY)calloc(matchArrayAlloc, sizeof(MATCHENTRY));
  if (pMatchArray == NULL) {
    printf("ERROR Failed to allocate pMatchArray\n");
    exit(-1);
  }
#endif

  candidateCount = 0;
  matchArraySize = 0;
  searchTableSize = 0;
  selectedCount = 0;
  memset(searchTable, 0, sizeof(searchTable));
  for (searchIndex = 0; searchIndex < SEARCH_TABLE_LENGTH; searchIndex++) {
    pSearchEntry1 = &searchTable[searchIndex];
    pSearchEntry1->sextractorIndex = -1;
  }

  minX = mosaicWidth * BORDER_EDGE;
  minY = mosaicHeight * BORDER_EDGE;
  maxX = mosaicWidth - minX;
  maxY = mosaicHeight - minY;
#if 0
  printf("Min, Max X: %d %d, Min Max Y %d %d\n", minX, maxX, minY, maxY);
#endif

  /* We assume that most of the objects are dim.  Calculate the clipped median to
     get the limiting dim MAG_ISO */
  if (vector != NULL) {
    free(vector);
    vector = NULL;
  }

  vectorSize = 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1);
  if (sextractor_nrecs > vectorSize) {
    vectorSize = sextractor_nrecs;
  }

  vector = (double *)calloc(vectorSize, sizeof(double));

#ifndef USE_ANALYSIS_THRESHOLD
  for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
    pSextractorBright = &image_table[sextractorIndex];
    vector[sextractorIndex] = pSextractorBright->MAG_ISO;
  }

  curClipCount = sextractor_nrecs;
  curClipCount = CalcMedianAndRMS(curClipCount, 10, vector, &MAG_ISO_med, &MAG_ISO_rms, 1, CLIP_RMS_FACTOR, 0);
  limiting_iso = MAG_ISO_med - (MAG_ISO_RMS_FACTOR * MAG_ISO_rms);

  if (verbose) {
    printf("curClipCount %d MAG_ISO median %f, MAG_ISO rms %f\n", curClipCount, MAG_ISO_med, MAG_ISO_rms);
  }
#endif /* USE_ANALYSIS_THRESHOLD */

  /* At this point, we will produce a sorted list of the brightest objects on the plate
     with a reasonable FWHM */

  for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
    pSextractorBright = &image_table[sextractorIndex];
    /* Reject saturated and bad images */
    if ((pSextractorBright->BFLAGS & (~((1 << FILTER_BFLAG_NEIGHBORS) | (1 << FILTER_BFLAG_BLEND)))) != 0) {
#if 0
      printf("Rejecting BFLAGS %d\n", pSextractorBright->BFLAGS);
#endif
      continue;
    }

    /* Reject images that are too large */
    if (pSextractorBright->FWHM_IMAGE > MAX_PIXEL_RADIUS) {
#if 0
      printf("Rejecting Object %6d MAG_ISO %10f FWHM_WORLD %d\n",
             sextractorIndex,
             pSextractorBright->MAG_ISO,
             pSextractorBright->FWHM_IMAGE);
#endif
      continue;
    }

#ifdef USE_ANALYSIS_THRESHOLD
    if ((pSextractorBright->FLUX_MAX/pSextractorBright->THRESHOLD) < analysis_threshold) {
      continue;
    }
#else /* USE_ANALYSIS_THRESHOLD */
    /* Reject images dimmer than the estimated limiting magnitude */
    if (pSextractorBright->MAG_ISO > limiting_iso) {
      continue;
    }
#endif /* ANALYSIS_THRESHOLD */

    if ((pSextractorBright->X_IMAGE < minX) ||
        (pSextractorBright->X_IMAGE > maxX) ||
        (pSextractorBright->Y_IMAGE < minY) ||
        (pSextractorBright->Y_IMAGE > maxY)) {
#if 0
      printf("IMAGE ON EDGE %10d %10d\n",
             pSextractorBright->X_IMAGE,
             pSextractorBright->Y_IMAGE);
#endif
      continue;
    }

    /* Insert this object in the search table */
    if ((searchTableSize < SEARCH_TABLE_LENGTH) ||
        (pSextractorBright->MAG_ISO < searchTable[SEARCH_TABLE_LENGTH-1].MAG_ISO)) {
      for (searchIndex = 0; searchIndex < SEARCH_TABLE_LENGTH; searchIndex++) {
        pSearchEntry1 = &searchTable[searchIndex];
        if (pSearchEntry1->sextractorIndex < 0) {
          pSearchEntry1->sextractorIndex = sextractorIndex;
          pSearchEntry1->MAG_ISO = pSextractorBright->MAG_ISO;
          pSearchEntry1->FLUX_MAX = pSextractorBright->FLUX_MAX;
          pSearchEntry1->THRESHOLD = pSextractorBright->THRESHOLD;
          pSearchEntry1->NUMBER = pSextractorBright->NUMBER;
          pSearchEntry1->X_IMAGE     = pSextractorBright->X_IMAGE;
          pSearchEntry1->Y_IMAGE     = pSextractorBright->Y_IMAGE;
          searchTableSize++;
          break;
        } else {
          if (pSextractorBright->MAG_ISO < pSearchEntry1->MAG_ISO) {
            /* Our entry goes here - move everything else down */
            for (insertIndex = (SEARCH_TABLE_LENGTH-1); insertIndex > searchIndex; insertIndex--) {
              memcpy(&searchTable[insertIndex], &searchTable[insertIndex-1], sizeof(SEARCHENTRY));
            }
            pSearchEntry1->sextractorIndex = sextractorIndex;
            pSearchEntry1->MAG_ISO = pSextractorBright->MAG_ISO;
            pSearchEntry1->FLUX_MAX = pSextractorBright->FLUX_MAX;
            pSearchEntry1->THRESHOLD = pSextractorBright->THRESHOLD;
            pSearchEntry1->NUMBER = pSextractorBright->NUMBER;
            pSearchEntry1->X_IMAGE     = pSextractorBright->X_IMAGE;
            pSearchEntry1->Y_IMAGE     = pSextractorBright->Y_IMAGE;
            if (searchTableSize < SEARCH_TABLE_LENGTH) {
              searchTableSize++;
            }
            break;
          }
        }
      } /* Insertion loop */
#if 0
      for (insertIndex = 0; insertIndex < SEARCH_TABLE_LENGTH; insertIndex++) {
        pSearchEntry2 = &searchTable[insertIndex];
        if (pSearchEntry2->sextractorIndex >= 0) {
          pSextractor2 = &image_table[pSearchEntry2->sextractorIndex];
          printf("Entry %7d  Index %7d MAG_ISO %10f NUMBER %6d X_IMAGE %10d Y_IMAGE %10d\n",
                 insertIndex,
                 pSearchEntry2->sextractorIndex,
                 pSearchEntry2->MAG_ISO,
                 pSextractor2->NUMBER,
                 pSextractor2->X_IMAGE,
                 pSextractor2->Y_IMAGE);

        }
      }
#endif
    } /* Star bright enough for insertion */
#if 0
    if (sextractorIndex > (2 * SEARCH_TABLE_LENGTH)) {
      printf("ERROR: early abort\n");
      exit(-1);
    }
#endif
  } /* Sextractor table index */

  correlation_table1 = (int*)calloc(4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1), sizeof(int));
  if (correlation_table1 == NULL) {
    printf("ERROR: failed to allocate correlation table of size %d\n", 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1));
  }

  correlation_table2 = (int*)calloc(4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1), sizeof(int));
  if (correlation_table2 == NULL) {
    printf("ERROR: failed to allocate correlation table of size %d\n", 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1));
  }

  correlation_table3 = (int*)calloc(4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1), sizeof(int));
  if (correlation_table3 == NULL) {
    printf("ERROR: failed to allocate correlation table of size %d\n", 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1));
  }

  qsort(searchTable, searchTableSize, sizeof(SEARCHENTRY), searchCompare);

#if 0
  for (insertIndex = 0; insertIndex < searchTableSize; insertIndex++) {
    pSearchEntry2 = &searchTable[insertIndex];
    if (pSearchEntry2->sextractorIndex >= 0) {
      pSextractor2 = &image_table[pSearchEntry2->sextractorIndex];
      printf("Entry %6d  Index %7d MAG_ISO %10f NUMBER %6d X_IMAGE %10d Y_IMAGE %10d\n",
             insertIndex,
             pSearchEntry2->sextractorIndex,
             pSearchEntry2->MAG_ISO,
             pSextractor2->NUMBER,
             pSextractor2->X_IMAGE,
             pSextractor2->Y_IMAGE);

    }
  }
#endif

  for (searchIndex = 0; searchIndex < searchTableSize; searchIndex++) {
    pSearchEntry1 = &searchTable[searchIndex];
    for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
      pSextractorDim = &image_table[sextractorIndex];
      if (pSearchEntry1->NUMBER == pSextractorDim->NUMBER) {
        pSextractorDim->searchStar = 1;
        break;
      }
    }
  }

  /* Now go through the sextractor list for each of these objects and save likely candidates */
  minSextractorIndex = 0;

  for (searchIndex = 0; searchIndex < searchTableSize; searchIndex++) {
    pSearchEntry1 = &searchTable[searchIndex];

    NUMBER = image_table[pSearchEntry1->sextractorIndex].NUMBER;

    minX = pSearchEntry1->X_IMAGE  - MAX_PIXEL_RADIUS;
    maxX = pSearchEntry1->X_IMAGE  + MAX_PIXEL_RADIUS;
    minY = pSearchEntry1->Y_IMAGE  - MAX_PIXEL_RADIUS;
    maxY = pSearchEntry1->Y_IMAGE  + MAX_PIXEL_RADIUS;

    for (sextractorIndex = minSextractorIndex; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
      pSextractorDim = &image_table[sextractorIndex];
      if (pSextractorDim->NUMBER == NUMBER) {
        pSextractor2 = pSextractorDim;
        continue;
      }

      if ((sextractorIndex > minSextractorIndex) &&
          (pSextractorDim->Y_IMAGE < pSextractor2->Y_IMAGE)) {
        printf("ERROR: sextractor table is not sorted by Y_IMAGE %d %d %s\n",
               pSextractorDim->Y_IMAGE,
               pSextractor2->Y_IMAGE,
               sextractor_name);
        exit(-1);
      }

      pSextractor2 = pSextractorDim;

      if (pSextractorDim->Y_IMAGE < minY) {
        minSextractorIndex = sextractorIndex;
        continue;
      }
      if (pSextractorDim->Y_IMAGE > maxY) {
        break;
      }
      if ((pSextractorDim->X_IMAGE < minX) ||
          (pSextractorDim->X_IMAGE > maxX)) {
        continue;
      }
      if (pSextractorDim->MAG_ISO < pSearchEntry1->MAG_ISO) {
        /* The search entry must be brighter */
        continue;
      }
#ifdef USE_ANALYSIS_THRESHOLD
      if ((pSextractorDim->FLUX_MAX/pSextractorDim->THRESHOLD) < analysis_threshold) {
        continue;
      }

#else /* USE_ANALYSIS_THRESHOLD */
         /* Reject images dimmer than the estimated limiting magnitude */
      if (pSextractorDim->MAG_ISO > limiting_iso) {
        continue;
      }
#endif /* ANALYSIS_THRESHOLD */

      xIndex = pSextractorDim->X_IMAGE - pSearchEntry1->X_IMAGE;
      yIndex = pSextractorDim->Y_IMAGE - pSearchEntry1->Y_IMAGE;
      if ((xIndex < -MAX_PIXEL_RADIUS) ||
          (xIndex >  MAX_PIXEL_RADIUS) ||
          (yIndex < -MAX_PIXEL_RADIUS) ||
          (yIndex >  MAX_PIXEL_RADIUS)) {
        continue;
      }
      pSextractorDim->convolutionFlag = 1;
      correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
      if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
        printf("ERROR: correlationIndex %d exceeds limits (5) %d\n", correlationIndex, maxCorrelationIndex);
        exit(-1);
      }
#if 0
      if ((xIndex == -1) && (yIndex == -136)) {
        printf("At x %d y %d NUMBER %d %d\n",
               xIndex, yIndex, pSextractorDim->NUMBER, pSearchEntry1->NUMBER);
      }

#endif
      correlation_table1[correlationIndex]++;

      if (correlation_table1[correlationIndex] > maxBinCount[0]) {
        maxBinCount[0] = correlation_table1[correlationIndex];
      }
    }
  }

  vectorIndex = 0;

  for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
    for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
      correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
      if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
        printf("ERROR: correlationIndex %d exceeds limits (6) %d\n", correlationIndex, maxCorrelationIndex);
        exit(-1);
      }

      if (doPlots) {
        if (gmtHandle[5] == NULL) {
          gmtHandle[5] = fopen(gmtname[5], "wt");
          if (gmtHandle[5] == NULL) {
            printf("ERROR: Failed to open the gmt file %s\n", gmtname[5]);
            exit(-1);
          } else {
            if (verbose) {
              printf("GMT file %s\n", gmtname[5]);
            }
          }
          fprintf(gmtHandle[5], "xIndex\tyIndex\tcorrelationCount\n");
          fprintf(gmtHandle[5], "------\t------\t----------------\n");
        }
        fprintf(gmtHandle[5], "%d\t%d\t%d\t\n", xIndex, yIndex, correlation_table1[correlationIndex]);
      }

      vector[vectorIndex] = correlation_table1[correlationIndex];
      vectorIndex++;
      if ((correlation_table1[correlationIndex]) > 0) {
        populatedBinCount++;
      }

    }
  }

  CalcMedianAndRMS(vectorIndex, 10, vector, &count_med[0], &count_rms[0], 1, 3.0, 0);

  if (count_rms[0] > 0) {
    maxBinSNR[0] = ((1.0*maxBinCount[0])-count_med[0])/count_rms[0];
  }

  count_limit[0] = count_med[0] + (COUNT_LIMIT_RMS_FACTOR*count_rms[0]);

  if (verbose) {
    printf("count_med %f count_rms %f, count_limit %d maxBinSNR %f for %s\n", count_med[0], count_rms[0], count_limit[0], maxBinSNR[0], fileroot);
  }

  if (doPlots) {
    if (gmtHandle[1] == NULL) {
      gmtHandle[1] = fopen(gmtname[1], "wt");
      if (gmtHandle[1] == NULL) {
        printf("ERROR: Failed to open the gmt file %s\n", gmtname[1]);
        exit(-1);
      } else {
        if (verbose) {
          printf("GMT file %s\n", gmtname[1]);
        }
      }
      fprintf(gmtHandle[1], "convolutionRadius\tconvolutionArea\tcount_med\tcount_rms\tcount_limit\tmaxBinCount\tmaxBinSNR\tselected\tmaskCount");
#ifdef USE_ANALYSIS_THRESHOLD
      fprintf(gmtHandle[1], "\tANALYSIS_THRESHOLD");
#else /* USE_ANALYSIS_THRESHOLD */
      fprintf(gmtHandle[1], "\tMAG_ISO_med\tMAG_ISO_rms");
#endif /* USE_ANALYSIS_THRESHOLD */
      fprintf(gmtHandle[1], "\tPlate\n");

      fprintf(gmtHandle[1], "-----------------\t---------------\t---------\t---------\t-----------\t-----------\t---------\t--------\t---------");
#ifdef USE_ANALYSIS_THRESHOLD
      fprintf(gmtHandle[1], "\t------------------");
#else /* USE_ANALYSIS_THRESHOLD */
      fprintf(gmtHandle[1], "\t-----------\t-----------");
#endif /* USE_ANALYSIS_THRESHOLD */
      fprintf(gmtHandle[1], "\t-----\n");
    }

    fprintf(gmtHandle[1], "%d\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d", convolutionRadius, convolutionArea[0], count_med[0], count_rms[0], (1.0*count_limit[0]), maxBinCount[0], maxBinSNR[0], selected, totalMasks);
#ifdef USE_ANALYSIS_THRESHOLD
    fprintf(gmtHandle[1], "\t%f", analysis_threshold);

#else /* USE_ANALYSIS_THRESHOLD */
    fprintf(gmtHandle[1], "\t%f\t%f\n", MAG_ISO_med, MAG_ISO_rms);
#endif /* USE_ANALYSIS_THRESHOLD */
    fprintf(gmtHandle[1], "\t%s\n", fileroot2);
  }

  /* Now perform a correlation analysis */
  for (convolutionRadius2 = 1; convolutionRadius2 < MAX_CONVOLUTION_RADIUS; convolutionRadius2++) {

#ifdef LOS_DEBUG
    if (convolutionRadius2 == (LOS_DEBUG+1)) {
      convolutionRadius2 = MAX_CONVOLUTION_RADIUS-1;
    } else {
      convolutionRadius2 = LOS_DEBUG;
    }
#endif /* LOS_DEBUG */
    if (convolutionRadius2 < (MAX_CONVOLUTION_RADIUS-1)) {
      convolutionRadius = convolutionRadius2;
    } else {
      /* Go back and select the correlation radius with the higest value */
      convolutionRadius = 1;
      maxBinSNRVal = maxBinSNR[1];
      max_count_limit = count_limit[1];
      for (convolutionRadius3 = 2; convolutionRadius3 < (MAX_CONVOLUTION_RADIUS-1); convolutionRadius3++) {
        if (maxBinSNRVal < maxBinSNR[convolutionRadius3]) {
          maxBinSNRVal = maxBinSNR[convolutionRadius3];
          max_count_limit = count_limit[convolutionRadius3];
          convolutionRadius = convolutionRadius3;
        }
      }
      selected = 1;
#if 0
      printf("Maximum SNR %f at radius %d for %s\n", maxBinSNRVal, convolutionRadius, fileroot);
#endif
    }

    memset(correlation_table2, 0, (4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1))*(sizeof(int)));
    convolutionArea[convolutionRadius] = 0;
    for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
      for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
        correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
        if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
          printf("ERROR: correlationIndex %d exceeds limits (7) %d\n", correlationIndex, maxCorrelationIndex);
          exit(-1);
        }
        convolutionCount = 0;
        for (xIndex2 = (xIndex-convolutionRadius); xIndex2 <= (xIndex+convolutionRadius); xIndex2++) {
          if ((xIndex2 < -MAX_PIXEL_RADIUS) || (xIndex2 > MAX_PIXEL_RADIUS)) {
            continue;
          }
          for (yIndex2 = yIndex-convolutionRadius; yIndex2 <= yIndex+convolutionRadius; yIndex2++) {
            if ((yIndex2 < -MAX_PIXEL_RADIUS) || (yIndex2 > MAX_PIXEL_RADIUS)) {
              continue;
            }
            convolutionCount++;
            correlationIndex2 = (xIndex2+MAX_PIXEL_RADIUS) + ((yIndex2+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
            if ((correlationIndex2 < 0) || (correlationIndex2 >= maxCorrelationIndex)) {
              printf("ERROR: correlationIndex %d exceeds limits (8) %d\n", correlationIndex, maxCorrelationIndex);
              exit(-1);
            }
            correlation_table2[correlationIndex2] += correlation_table1[correlationIndex];
            if (correlation_table2[correlationIndex2] > maxBinCount[convolutionRadius]) {
              maxBinCount[convolutionRadius] = correlation_table2[correlationIndex2];
            }
          }
        }
        if (convolutionCount > convolutionArea[convolutionRadius]) {
          convolutionArea[convolutionRadius] = convolutionCount;
        }
      }
    }

    vectorIndex = 0;

    for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
      for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
        correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
        if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
          printf("ERROR: correlationIndex %d exceeds limits (9) %d\n", correlationIndex, maxCorrelationIndex);
          exit(-1);
        }
        vector[vectorIndex] = correlation_table2[correlationIndex];
        vectorIndex++;
      }
    }

    CalcMedianAndRMS(vectorIndex, 10, vector, &count_med[convolutionRadius], &count_rms[convolutionRadius], 1, 3.0, 0);

    if (count_rms[convolutionRadius] > 0) {
      maxBinSNR[convolutionRadius] = ((1.0*maxBinCount[convolutionRadius])-count_med[convolutionRadius])/count_rms[convolutionRadius];
    }
    count_limit[convolutionRadius] = count_med[convolutionRadius] + (COUNT_LIMIT_RMS_FACTOR*count_rms[convolutionRadius]);
    if (verbose) {
      printf("convolutionRadius %d convolutionArea %d count_med %f count_rms %f, count_limit %f maxBinCount %d maxBinSNR %f for %s\n", convolutionRadius, convolutionArea[convolutionRadius], count_med[convolutionRadius], count_rms[convolutionRadius], (1.0*count_limit[convolutionRadius]), maxBinCount[convolutionRadius], maxBinSNR[convolutionRadius], fileroot);
    }
    if (doPlots) {
      if (selected == 0) {
        fprintf(gmtHandle[1], "%d\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d", convolutionRadius, convolutionArea[convolutionRadius], count_med[convolutionRadius], count_rms[convolutionRadius], (1.0*count_limit[convolutionRadius]), maxBinCount[convolutionRadius], maxBinSNR[convolutionRadius], selected, totalMasks);
#ifdef USE_ANALYSIS_THRESHOLD
        fprintf(gmtHandle[1], "\t%f", analysis_threshold);

#else /* USE_ANALYSIS_THRESHOLD */
        fprintf(gmtHandle[1], "\t%f\t%f\n", MAG_ISO_med, MAG_ISO_rms);
#endif /* USE_ANALYSIS_THRESHOLD */
        fprintf(gmtHandle[1], "\t%s\n", fileroot2);
      }
    }
  }

  if (doPlots) {
    for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
      for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
        correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
        if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
          printf("ERROR: correlationIndex %d exceeds limits (10) %d\n", correlationIndex, maxCorrelationIndex);
          exit(-1);
        }
        if ((correlation_table2[correlationIndex]) > count_limit[convolutionRadius]) {
          if (gmtHandle[0] == NULL) {
            gmtHandle[0] = fopen(gmtname[0], "wt");
            if (gmtHandle[0] == NULL) {
              printf("ERROR: Failed to open the gmt file %s\n", gmtname[0]);
              exit(-1);
            } else {
              if (verbose) {
                printf("GMT file %s\n", gmtname[0]);
              }
            }
            fprintf(gmtHandle[0], "ix\tiy\tbin_count\n");
            fprintf(gmtHandle[0], "--\t--\t---------\n");
          }
          fprintf(gmtHandle[0], "%d\t%d\t%d\n", xIndex, yIndex, correlation_table2[correlationIndex]);
        }

      }
    }
  }

  /* Now we need to identify the peak correlation masks and print a table of mask characteristics */

  FindMasks(count_limit[convolutionRadius], correlation_table2, correlation_table3, &maskTable, &totalMasks, fileroot);
  if (totalMasks >= SEARCH_TABLE_LENGTH) {
    printf("ERROR search_close totalMasks %d greater then SEARCH_TABLE_LENGTH %d for %s\n", totalMasks, SEARCH_TABLE_LENGTH, fileroot);
    exit(-1);
  }

  if (doPlots) {
    fprintf(gmtHandle[1], "%d\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d", convolutionRadius, convolutionArea[convolutionRadius], count_med[convolutionRadius], count_rms[convolutionRadius], (1.0*count_limit[convolutionRadius]), maxBinCount[convolutionRadius], maxBinSNR[convolutionRadius], selected, totalMasks);
#ifdef USE_ANALYSIS_THRESHOLD
    fprintf(gmtHandle[1], "\t%f", analysis_threshold);

#else /* USE_ANALYSIS_THRESHOLD */
    fprintf(gmtHandle[1], "\t%f\t%f\n", MAG_ISO_med, MAG_ISO_rms);
#endif /* USE_ANALYSIS_THRESHOLD */
    fprintf(gmtHandle[1], "\t%s\n", fileroot2);
  }

  /* Expand all entries in the mask table by the correlation radius and find the limits.  Set the centerFlag if necessary */
  for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
    pMask = &maskTable[maskIndex];
    pMask->maskXMin -= convolutionRadius;
    pMask->maskXMax += convolutionRadius;
    pMask->maskYMin -= convolutionRadius;
    pMask->maskYMax += convolutionRadius;
    if ((pMask->maskXMin <= 0) &&
        (pMask->maskXMax >= 0) &&
        (pMask->maskYMin <= 0) &&
        (pMask->maskYMax >= 0)) {
      pMask->maskCenterFlag = 1;
      maskCenterIndex = 1;
    }
  }

  /* Here we need to know the relative brightness of the secondary images.  Because MAG_ISO is nonlinear, we can only compare objects for a given host star.
     1.  For each host star, find the brightest mask candidate that matches the star.
     2.  Sort the mask candidates to get a sequence of masks from brightest to dimmest
     3.  Using the sorted sequence, increment the mask brightness count table entry [i, j] if mask i is brighter than mask j
     4.  After all host stars are processed find the order of masks based on the highest mask brightness counts
  */
  if (totalMasks > 1) {
    maskOrderTable = (PMASKORDER)calloc(totalMasks, sizeof(MASKORDER));
    if (maskOrderTable == NULL) {
      printf("ERROR: failed to allocate maskOrderTable of size %d\n", totalMasks);
      exit(-1);
    }
    memset(maskBrightnessTable, 0, sizeof(maskBrightnessTable));
    for (searchIndex = 0; searchIndex < searchTableSize; searchIndex++) {
      pSearchEntry1 = &searchTable[searchIndex];
      for (maskIndex = 0; maskIndex < totalMasks; maskIndex++) {
        pMaskOrder = &maskOrderTable[maskIndex];
        pMaskOrder->maskIndex = -1;
        pMaskOrder->MAG_ISO = 99.00;
      }

      for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
        pSextractorDim = &image_table[sextractorIndex];
        if (pSextractorDim->NUMBER == pSearchEntry1->NUMBER) {
          continue;
        }
        if (pSextractorDim->MAG_ISO < pSearchEntry1->MAG_ISO) {
          continue;
        }
#ifdef USE_ANALYSIS_THRESHOLD
        if ((pSextractorDim->FLUX_MAX/pSextractorDim->THRESHOLD) < analysis_threshold) {
          continue;
        }
#else /* USE_ANALYSIS_THRESHOLD */
         /* Reject images dimmer than the estimated limiting magnitude */
        if (pSextractorDim->MAG_ISO > limiting_iso) {
          continue;
        }
#endif /* ANALYSIS_THRESHOLD */
        for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
          pMask = &maskTable[maskIndex];
          pMaskOrder = &maskOrderTable[maskIndex];
          if (((pMask->maskXMin+pSearchEntry1->X_IMAGE) <= (pSextractorDim->X_IMAGE)) &&
              ((pMask->maskXMax+pSearchEntry1->X_IMAGE) >= (pSextractorDim->X_IMAGE)) &&
              ((pMask->maskYMin+pSearchEntry1->Y_IMAGE) <= (pSextractorDim->Y_IMAGE)) &&
              ((pMask->maskYMax+pSearchEntry1->Y_IMAGE) >= (pSextractorDim->Y_IMAGE))) {
            /* We have a match! */
            if (pMaskOrder->MAG_ISO > pSextractorDim->MAG_ISO) {
              pMaskOrder->maskIndex = maskIndex;
              pMaskOrder->MAG_ISO = pSextractorDim->MAG_ISO;
            }
          }
        } /* Dim star loop */
      } /* mask loop */

      for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
        pMaskOrder = &maskOrderTable[maskIndex];
        if (pMaskOrder->maskIndex > 0) {
          break;
        }
      }

      if (maskIndex < totalMasks) {
        /* Sort the mask table in order of decreasing brightness */
        qsort(maskOrderTable, totalMasks, sizeof(MASKORDER), maskOrderCompare);

#if 0
        for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
          pMaskOrder = &maskOrderTable[maskIndex];
          printf("%6.2f ", pMaskOrder->MAG_ISO);
        }
        printf("\n");
#endif
        for (maskIndex  = 0; maskIndex < (totalMasks-1); maskIndex++) {
          pMaskOrder = &maskOrderTable[maskIndex];
          if (pMaskOrder->maskIndex < 0) {
            continue;
          }
          for (maskIndex2 = maskIndex+1; maskIndex2 < totalMasks; maskIndex2++) {
            pMaskOrder2 = &maskOrderTable[maskIndex2];
            if (pMaskOrder2->maskIndex < 0) {
              continue;
            }
            if ((pMaskOrder->maskIndex >= SEARCH_TABLE_LENGTH) ||
                (pMaskOrder2->maskIndex >= SEARCH_TABLE_LENGTH)) {
              printf("ERROR: search_close maskIndex %d %d greater than %d for %s\n", pMaskOrder->maskIndex, pMaskOrder2->maskIndex, SEARCH_TABLE_LENGTH, fileroot);
              exit(-1);
            }
            maskBrightnessTable[pMaskOrder->maskIndex][pMaskOrder2->maskIndex]++;
          }
        }
      }

    } /* Bright star loop */
#if 0
    for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
      pMask = &maskTable[maskIndex];

      printf("%2d %3d:  ", pMask->maskIndex, maskIndex);
      for (maskIndex2 = 0; maskIndex2 < totalMasks; maskIndex2++) {
        printf("%4d ", maskBrightnessTable[maskIndex][maskIndex2]);
      }
      printf("\n");
    }
#endif
    /* Now rate the masks in terms of the number of times a given mask is brighter than another */
    for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
      pMask = &maskTable[maskIndex];
      pMask->maskIndex = 0;
      for (maskIndex2 = 0; maskIndex2 < totalMasks; maskIndex2++) {
        if (maskIndex == maskIndex2) {
          continue;
        }
        if (maskBrightnessTable[maskIndex][maskIndex2] > maskBrightnessTable[maskIndex2][maskIndex]) {
          pMask->maskIndex++;
        }
      }
    }

#if 0
    for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
      pMask = &maskTable[maskIndex];

      printf("%2d %3d:  ", pMask->maskIndex, maskIndex);
      for (maskIndex2 = 0; maskIndex2 < totalMasks; maskIndex2++) {
        printf("%4d ", maskBrightnessTable[maskIndex][maskIndex2]);
      }
      printf("\n");
    }
#endif
    qsort(maskTable, totalMasks, sizeof(MULTMASK), maskBrightnessCompare);
#if 0
    for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
      pMask = &maskTable[maskIndex];

      printf("%2d %3d:  ", pMask->maskIndex, maskIndex);
      for (maskIndex2 = 0; maskIndex2 < totalMasks; maskIndex2++) {
        printf("%4d ", maskBrightnessTable[maskIndex][maskIndex2]);
      }
      printf("\n");
    }
#endif
  }

  /* Save the results in the database */
  if (totalMasks > 0) {
    SetMosaicFitWCS(pConnection, "HaveMultipleMask", fileroot, 0, 1, readFlag, &FitWCS);
    if ((totalMasks <= MAX_MASKS) &&
        (maskCenterIndex == 0)) {
      SetPlateQuality(pConnection, "multiple", fileroot, 1, 0, &quality);
    }
  } else {
    SetMosaicFitWCS(pConnection, "NoMultipleMask", fileroot, 0, 1, readFlag, &FitWCS);
  }

  for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
    pMask = &maskTable[maskIndex];
    pMask->maskIndex = maskIndex+1;
    /* Load the masks into the datbase, but only 1 mask if there is an unusual condition */
    if (((maskCenterIndex == 0) && (totalMasks <= MAX_MASKS)) ||
        ((maskCenterIndex != 0) && (pMask->maskCenterFlag != 0)) ||
        ((totalMasks > MAX_MASKS) && (maskIndex == 0))) {

      sprintf(queryString, "INSERT INTO masks (series, plateNumber, mosaicNumber, maskCount, maskIndex, maskXMin, maskXMax, maskYMin, maskYMax, maskArea, maxBinCount, maskCenterFlag, maskCenterDistance, maskAreaRatio, convolutionRadius) values ('%s', %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %f, %f, %d);",
              series,
              plateNumber,
              mosaicNumber,
              totalMasks,
              pMask->maskIndex,
              pMask->maskXMin,
              pMask->maskXMax,
              pMask->maskYMin,
              pMask->maskYMax,
              pMask->maskArea,
              pMask->maxBinCount,
              pMask->maskCenterFlag,
              pMask->maskCenterDistance,
              pMask->maskAreaRatio,
              convolutionRadius);
      if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
        printf("ERROR: search_close MAX_QUERY_STRING exceeded\n");
        exit(-1);
      }
      res = ExecuteQuery(pConnection, queryString);
      if (res) {
        exit(-1);
      }
    }

    if (doPlots) {
      if (gmtHandle[2] == NULL) {
        gmtHandle[2] = fopen(gmtname[2], "wt");
        if (gmtHandle[2] == NULL) {
          printf("ERROR: Failed to open the gmt file %s\n", gmtname[2]);
          exit(-1);
        } else {
          if (verbose) {
            printf("GMT file %s\n", gmtname[2]);
          }
        }

        fprintf(gmtHandle[2], "maskCount\tmaskIndex\tmaskXMin\tmaskXMax\tmaskYMin\tmaskYMax\tmaskArea\tmaxBinCount\tmaskCenterFlag\tmaskCenterDistance\tmaskAreaRatio\tconvolutionRadius\tPlate\n");
        fprintf(gmtHandle[2], "---------\t---------\t--------\t--------\t--------\t--------\t--------\t-----------\t--------------\t------------------\t-------------\t-----------------\t-----\n");
      }

      if (gmtHandle[6] == NULL) {
        gmtHandle[6] = fopen(gmtname[6], "wt");
        if (gmtHandle[6] == NULL) {
          printf("ERROR: Failed to open the gmt file %s\n", gmtname[6]);
          exit(-1);
        } else {
          if (verbose) {
            printf("GMT file %s\n", gmtname[6]);
          }
        }
        fprintf(gmtHandle[6], "xIndex\tyIndex\tmaskIndex\n");
        fprintf(gmtHandle[6], "------\t------\t---------\n");
      }

      fprintf(gmtHandle[2], "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.3f\t%d\t%s\n",
              totalMasks,
              pMask->maskIndex,
              pMask->maskXMin,
              pMask->maskXMax,
              pMask->maskYMin,
              pMask->maskYMax,
              pMask->maskArea,
              pMask->maxBinCount,
              pMask->maskCenterFlag,
              pMask->maskCenterDistance,
              pMask->maskAreaRatio,
              convolutionRadius,
              fileroot);
      for (xIndex = pMask->maskXMin; xIndex <= pMask->maskXMax; xIndex++) {
        for (yIndex = pMask->maskYMin; yIndex <= pMask->maskYMax; yIndex++) {
          fprintf(gmtHandle[6], "%d\t%d\t\%d\n", xIndex, yIndex, pMask->maskIndex);
        }
      }
    }
  }

  /* Update the quality bits for this plate */
  UpdateQuality(pConnection, series, plateNumber);

  /* Now we need to search through the Sextractor file and flag any objects that fall within the masks */

  if ((totalMasks > 0) && (totalMasks <= MAX_MASKS) && (maskCenterIndex == 0)) {
    for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
      pMask = &maskTable[maskIndex];
      if (maskIndex == 0) {
        minX = pMask->maskXMin;
        maxX = pMask->maskXMax;
        minY = pMask->maskYMin;
        maxY = pMask->maskYMax;
      } else {
        if (minX > pMask->maskXMin) {
          minX = pMask->maskXMin;
        }
        if (maxX < pMask->maskXMax) {
          maxX = pMask->maskXMax;
        }
        if (minY > pMask->maskYMin) {
          minY = pMask->maskYMin;
        }
        if (maxY < pMask->maskYMax) {
          maxY = pMask->maskYMax;
        }
      }
    }

    /* Now for every sextractor object pair, look for a match in the mask table */
    minSextractorIndex = 0;
    for (sextractorIndex1 = 0; sextractorIndex1 < sextractor_nrecs; sextractorIndex1++) {
      pSextractor1 = &image_table[sextractorIndex1];
#ifdef USE_ANALYSIS_THRESHOLD
      if ((pSextractor1->FLUX_MAX/pSextractor1->THRESHOLD) < analysis_threshold) {
        continue;
      }

#else /* USE_ANALYSIS_THRESHOLD */
         /* Reject images dimmer than the estimated limiting magnitude */
      if (pSextractor1->MAG_ISO > limiting_iso) {
        continue;
      }
#endif /* ANALYSIS_THRESHOLD */

      for (sextractorIndex2 = minSextractorIndex; sextractorIndex2 < sextractor_nrecs; sextractorIndex2++) {
        pSextractor2 = &image_table[sextractorIndex2];
        xIndex = pSextractor2->X_IMAGE - pSextractor1->X_IMAGE;
        yIndex = pSextractor2->Y_IMAGE - pSextractor1->Y_IMAGE;
#if 0
        /* pSextractor1  MAG_ISO = -14.4598, NUMBER = 179928, X_IMAGE = 13917, Y_IMAGE = 7607 */
        /* pSextractor2  MAG_ISO = -11.2915, NUMBER = 175947, X_IMAGE = 13916, Y_IMAGE = 7471 */
        /* xIndex = -1 yIndex = -136  for /mc05077_01_01r270ww */
        if ((pSextractor1->NUMBER == 179928) &&
            (pSextractor2->NUMBER == 175947)) {
          printf("At x %d y %d NUMBER %d %d\n",
                 xIndex, yIndex, pSextractor1->NUMBER, pSextractor2->NUMBER);
        }
#endif

        if (yIndex < minY) {
          minSextractorIndex = sextractorIndex2;
          continue;
        }
        if (yIndex > maxY) {
          break;
        }
        if (pSextractor2->MAG_ISO < pSextractor1->MAG_ISO) {
          continue;
        }
        if (pSextractor2->NUMBER == pSextractor1->NUMBER) {
          continue;
        }
#ifdef USE_ANALYSIS_THRESHOLD
        if ((pSextractor2->FLUX_MAX/pSextractor2->THRESHOLD) < analysis_threshold) {
          continue;
        }

#else /* USE_ANALYSIS_THRESHOLD */
         /* Reject images dimmer than the estimated limiting magnitude */
        if (pSextractor2->MAG_ISO > limiting_iso) {
          continue;
        }
#endif /* ANALYSIS_THRESHOLD */
        if ((pSextractor2->maskIndex == 0) &&
            (xIndex >= minX) &&
            (xIndex <= maxX)) {
          /* In the correct zone.  Look for a match with the table entries */
          for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
            pMask = &maskTable[maskIndex];
            if ((xIndex >= pMask->maskXMin) &&
                (xIndex <= pMask->maskXMax) &&
                (yIndex >= pMask->maskYMin) &&
                (yIndex <= pMask->maskYMax)) {
              pSextractor2->maskIndex = pMask->maskIndex;
              if (pSextractor1->searchStar != 0) {
                pSextractor2->pSextractor1 = pSextractor1;
              }
              outputCount++;
              break;
            }
          }
        }
      }
    }
  }

  if (outputCount > 0) {
    for (sextractorIndex1 = 0; sextractorIndex1 < sextractor_nrecs; sextractorIndex1++) {
      pSextractor1 = &image_table[sextractorIndex1];
      if (pSextractor1->maskIndex == 0) {
        continue;
      }
#ifdef USE_ANALYSIS_THRESHOLD
      if ((pSextractor1->FLUX_MAX/pSextractor1->THRESHOLD) < analysis_threshold) {
        continue;
      }
#else /* USE_ANALYSIS_THRESHOLD */
         /* Reject images dimmer than the estimated limiting magnitude */
      if (pSextractor1->MAG_ISO > limiting_iso) {
        continue;
      }
#endif /* ANALYSIS_THRESHOLD */

      if (outHandle == NULL) {
        outHandle = fopen(outfile, "wt");
        if (outHandle == NULL) {
          fprintf(stderr, "ERROR: Failed to open the output file %s\n", outfile);
          exit(-1);
        } else {
          if (verbose) {
            fprintf(stderr, "Output file %s\n", outfile);
          }
          fprintf(outHandle, "NUMBER\tmaskIndex\tX_IMAGE\tY_IMAGE\n");
          fprintf(outHandle, "------\t---------\t-------\t-------\n");

        }
      }

      fprintf(outHandle, "%d\t%d\t%d\t%d\n",
              pSextractor1->NUMBER,
              pSextractor1->maskIndex,
              pSextractor1->X_IMAGE,
              pSextractor1->Y_IMAGE);

      if (doPlots) {
        if (gmtHandle[4] == NULL) {
          gmtHandle[4] = fopen(gmtname[4], "wt");
          if (gmtHandle[4] == NULL) {
            printf("ERROR: Failed to open the gmt file %s\n", gmtname[4]);
            exit(-1);
          } else {
            if (verbose) {
              printf("GMT file %s\n", gmtname[4]);
            }
          }
        }

        fprintf(gmtHandle[4], "CIRCLE(%d, %d, 5) # text = {%d}\n", pSextractor1->X_IMAGE/binning, pSextractor1->Y_IMAGE/binning, pSextractor1->maskIndex);
      }

      if ((pSextractor1->convolutionFlag == 0) || (pSextractor1->pSextractor1 == NULL)) {
        continue;
      }

      if (doPlots) {
        if (gmtHandle[3] == NULL) {
          gmtHandle[3] = fopen(gmtname[3], "wt");
          if (gmtHandle[3] == NULL) {
            printf("ERROR: Failed to open the gmt file %s\n", gmtname[3]);
            exit(-1);
          } else {
            if (verbose) {
              printf("GMT file %s\n", gmtname[3]);
            }
          }
        }

        fprintf(gmtHandle[3], "LINE(%d, %d, %d, %d) # text = {%d} line = 0 1 \n", pSextractor1->pSextractor1->X_IMAGE/binning, pSextractor1->pSextractor1->Y_IMAGE/binning, pSextractor1->X_IMAGE/binning, pSextractor1->Y_IMAGE/binning, pSextractor1->maskIndex);
      }
    }
  }

  /* All done.  Clean up */

  if (outHandle != NULL) {
    fclose(outHandle);
  }

  if (doPlots) {
    for (gmtIndex = 0; gmtIndex < NUM_GMT_FILES; gmtIndex++) {
      if (gmtHandle[gmtIndex] != NULL) {
        fclose(gmtHandle[gmtIndex]);
      }
    }
  }

  if (pMatchArray != NULL) {
    free(pMatchArray);
  }

  if (pSelectedMatchArray != NULL) {
    free(pSelectedMatchArray);
  }

  if (pMagsumIndexList != NULL) {
    free(pMagsumIndexList);
  }

  if (vector != NULL) {
    free(vector);
  }

  if (image_table != NULL) {
    free(image_table);
  }

  if (correlation_table1 != NULL) {
    free(correlation_table1);
  }

  if (correlation_table2 != NULL) {
    free(correlation_table2);
  }

  if (maskTable != NULL) {
    free(maskTable);
  }

  if (maskOrderTable != NULL) {
    free(maskOrderTable);
  }

  if (vector == 0) {
    free(vector);
  }

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (backgroundVector[spatial_bin] == 0) {
      free(backgroundVector[spatial_bin]);
    }
  }

  if (backgroundHandle != NULL) {
    fclose(backgroundHandle);
  }

#ifdef BACKGROUND_PLOTS
  if (regionHandle != NULL) {
    fclose(regionHandle);
  }

  if (backgroundmapHandle != NULL) {
    fclose(backgroundmapHandle);
  }

  if (imagemapHandle != NULL) {
    fclose(imagemapHandle);
  }
#endif /* BACKGROUND_PLOTS */

  if (pBinMap != NULL) {
    free(pBinMap);
  }

  mysql_close(pConnection);

  time(&curTime);
  curTime -= startTime;

  if (verbose) {
    printf("Candidates: %d MatchArray: %d selected: %d plotted %d psfsaturated %d populatedBinCount %d maxBinCount %d outputCount %d for %s\n",
           candidateCount,
           matchArraySize,
           selectedCount,
           plotCount,
           psfsaturatedCount,
           populatedBinCount,
           maxBinCount[0],
           outputCount,
           fileroot);
  }

  printf("Maximum SNR %f at radius %d", maxBinSNRVal, convolutionRadius);
#ifdef USE_ANALYSIS_THRESHOLD
  printf(" threshold %f", analysis_threshold);
#endif /* USE_ANALYSIS_THRESHOLD */
  printf(" count_limit %d totalMasks %d seconds %ld for %s\n", max_count_limit, totalMasks, curTime, fileroot);

  return(0);
}

/*
 * Dec 20, 2011 Edward J. Los - Adapted from filter_wedge
 * Jan 24, 2012 Edward J. Los - Switch to ANALYSIS_THRESH approximately equal to FLUX_MAX/THRESHOLD
 * Feb 13, 2012 Edward J. Los - Consider only objects above the threashold when creating the masks.
 * Feb 22, 2012 Edward J. Los - If we have valid masks, set the "multiple" bit on the plate.
 * Feb 24, 2012 Edward J. Los - Add UpdateQuality
 * Mar  5, 2012 Edward J. Los - Do not output the table if the masks are invalid.
 * Apr 10, 2012 Edward J. Los - Add emulsion plot files
 * May 12, 2012 Edward J. Los - Revise emulsion plots to filter grain noise
 * Jun  1, 2012 Edward J. Los - Add search for calibration squares
 * Jun 25, 2012 Edward J. Los - Do background calculations on a spatial_bin basis.
 * Jul  3, 2012 Edward J. Los - Integrate high background object support with rest of pipeline
 * Apr 21, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
 * Mar  3, 2015 Edward J. Los - Correct segmentation fault when totalMasks exceeds SEARCH_TABLE_LENGTH
 * Jun 13, 2017 Edward J. Los - Correct usage
 *                              Add gmt5 which is a dump of the initial correlation table.
 *                              Add gmt6 which is a map of the final masks
 */
