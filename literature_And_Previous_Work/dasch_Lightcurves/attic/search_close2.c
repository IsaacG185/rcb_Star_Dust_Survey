// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* search_close2.c
 *  
 *  Search for multiple images near the brightest stars.  This is a revision of the search_close program
 *  On Mon, 8 May 2017 16:09:28 approximately 17% of 494811 plates are known to be multiple, including grating and wedge plates
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm -L/usr/lib64/mysql -L/usr/lib/mysql -lmysqlclient search_close2.c pipelineutils.a -ltable -lutil -lwcs -o search_close2 
 *`
 *  search_close2 -v -s 1.786 -w 17412 -h 22026 -r i3s1090_00_01r180ww -i /home/scanner/Pipeline/match/i31090_00_01r180ww_tnx.tmp -o /home/scanner/Pipeline/match/i31090_00_01r180ww_close.db
 *
 * <plate>_close_gmt0.txt - 3d plot of (x,y,bin_count)
 * <plate>_close_gmt1.txt - convolution statistics
 * <plate>_close_gmt2.txt - mask table
 * <plate>_close_gmt3.txt - region file for masked points used to create the convolution table
 * <plate>_close_gmt4.txt - region file for all masked points
 * <plate>_close_gmt5.txt - correlationCount for unbinned results
 * <plate>_close_gmt6.txt - dump of mask indices for debug
 * <plate>_close_gmt7.txt - all values of MAG_ISO vs FLUX_MAX vs analysis_threshold 
 * <plate>_close_gmt8.txt - new peak ratios for truncated limiting magnitudes
 * <plate>_close_gmt9.txt - CountMasks result for entire magnitude range
 * <plate>_close_plot.db  - maps of sextractor images
 *
 * Matching limiting_iso
   cat *_gmt1.txt | sort -u > los.db
   row 'selected > 0' < los.db | column Plate MAG_ISO_med MAG_ISO_rms | sort > los1.db
   select series,plateNumber,spatial_bin,limiting_iso from spatialbin<n> inner join photseries using(seriesId);
   perl ~/perl/zerofill2.perl
   ~/Pipeline/makestarbase.csh los1.tmp los1.db
   sorttable Plate < los.db > los2.db
   sorttable Plate < los1.db > los3.db
   jointable -j Plate los2.db los3.db > los4.db

 *
 * background map
 *
   column -i i31090_00_01r180ww_tnx.db NUMBER X_IMAGE Y_IMAGE | sorttable -n NUMBER > los1.db
   jointable -n -j NUMBER los1.db i31090_00_01r180ww_background.db > los2.db
   votable -i los2.db -o los.xml

 *
 * Posisson formula (Drake: Fundamentals of applied probability theory 1967 p. 136
 *  mean = lambda * t = mu
 *  pmf = (mu**k)exp(-mu)/(k!) where k is the expected count in a pixel.
 *  below mu is meanArrival
 *        exp(-k) is exponential
 *        (mu**k) is prefix
 *        
 * New database variables (Some of these may now be unused)
       maskSector
       nxs
       nys
       matchRatio2        
       searchTableSize
       countRMSFactor
       maskKernelType
       bin_sum
       xIndex
       yIndex
       correlationCount
       MAG_ISO_diff
       MAG_ISO_flag
       pixDist
       NUMBER_1
       NUMBER_2
       maxCorrelationIndex
       pmf2
       pmf2_valid
       
       MAG_ISO_dim
       MAG_ISO_bright
       analysis_threshold
       countmatch1
       countmatch2

       table1Peak
       table1MaxRatio
       table2Peak
       table2Ratio
       table1_dimlimit
 *

 * old search_close log:
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
 * Sep  5, 2017 Edward J. Los - Create new search_close2 revision.
 *                              Remove USE_ANALYSIS_THRESHOLD (make it default)
 *                              Add DEBUG_TRUNCATION for testing purposes to shorten execution
 *                              Implement Poisson statistics for estimating the significance
 *                              Bin the mosaic in approximately 1" squares
 *                              Do not
 * Dec 21, 2017 Edward J. Los   Implement new algorithm to observe correlation as a function of limiting magnitude
 * Jan  3, 2017 Edward J. Los   redefine the table1Ratio to be the number of peaks at 50% divided by the number of peaks at 100%
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"

#define MAX_SECTOR_COUNT 252 /* Current maximum sector count */
#define DEFAULT_CONVOLUTION_RADIUS 1
#define MAX_BUFFER 100
#define SEARCH_TABLE_LENGTH 2500
#define NUM_GMT_FILES 10
#define MATCH_ARRAY_ALLOC 1000000
#define MAG_ISO_RMS_FACTOR 4
#define SQUARE_BIN_OBJECTS 10
#define CLIP_RMS_FACTOR 1.0
#define EXPECTED_PEAKS 1
#define EXCESS_LIMIT 10.0
#define MAX_PIXEL_RADIUS 1000
/* #define LOS_DEBUG 3   This is the desired bin number - use 3 for mc05077 */
#define DEBUG_MAG_ISO 1
#define DEBUG_TRUNCATION 1 /* Halt processing early (ave time for 72 mosaics: 1305 seconds median 1134 seconds 5.67 hr total */
#define DEBUG_SINGLE_SECTOR 1 /* Put everything into one sector */
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
#define MIN_CORRELATION_PROBABILITY 10.0 /* Clamp the pmf2 minimum value to this level */
#define MAG_ISO_DELTA 0.1 /* Resolution of the MAG_ISO_count correlation table */
/* #define BACKGROUND_PLOTS 1 */
#ifdef DEBUG_SINGLE_SECTOR
#define SECTOR_PIXEL_WIDTH 211600
#else /* DEBUG_SINGLE_SECTOR */
#define SECTOR_PIXEL_WIDTH 2116 /* Minimum width of a sector = 23.27 mm or 0.91 inch */
#endif /* DEBUG_SINGLE_SECTOR */

/* Threshold = 2.0, use 8-15 */
#if 0
#define COUNT_LIMIT_RMS_FACTOR 15
#else
#define COUNT_LIMIT_RMS_FACTOR 8
#endif

int maskBrightnessTable[SEARCH_TABLE_LENGTH][SEARCH_TABLE_LENGTH];

typedef struct _backgroundbin {
  double X_IMAGE;     /* X Center point */
  double Y_IMAGE;     /* Y Center point */
  double BACKGROUND;  /* Average background in bin */
  double FILTERED_BACKGROUND; /* Smoothed background */
  double expansionRadius;
  int count;         /* Count of objects in the bin */
} BACKGROUNDBIN,*PBACKGROUNDBIN;

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
  int maskSector;
} INPUTIMAGE, *PINPUTIMAGE;

typedef struct _searchentry {
  int sextractorIndex;          /* Index into the sextractor table */
  double MAG_ISO;     /* Brightness of this object */
  double FLUX_MAX;
  double THRESHOLD;
  int NUMBER;
  int X_IMAGE;
  int Y_IMAGE;
} SEARCHENTRY,*PSEARCHENTRY;

typedef struct _starimage {
  double MAG_ISO;
  double FLUX_MAX;
  double THRESHOLD;
  double X_IMAGE;
  double Y_IMAGE;
  int NUMBER;         /* Sextractor reference number */
  int MATCH_NUMBER;   /* If nonzero, matching bright star */
  int BFLAGS;          /* flags word */
  int FWHM_IMAGE;
  int maskIndex; /* Matching mask index */
  int convolutionFlag; /* Image contributed to the convolution table if nozero.   */
  int searchStar;      /* Image is in the search table if nonzero */
  int maskSector;
  struct _starimage *pSextractor1;
} STARIMAGE,*PSTARIMAGE;

#define MAX_SECTOR_NEIGHBOR 9 /* Includes this current sector */
typedef struct _masksector {
  int maskSector;
  int neighborCount;
  int maskNeighbor[MAX_SECTOR_NEIGHBOR];
  int sextractorIndexLow;
  int sextractorIndexHigh;
  int maxBinCount0; /* For convolution radius = 0 */
  double maxBinSNR0;
  double count_med0;
  double count_rms0;
  int count_limit0;  
  int maxBinCount1; /* For convolution radius = 1 */
  double maxBinSNR1;
  double count_med1;
  double count_rms1;
  int count_limit1;  
  int convolutionArea; /* Area in pixels of the convolution function */
  double aveX;
  double aveY;
  int sectorCount;
  int unsaturatedCount;
} MASKSECTOR,*PMASKSECTOR;

SEARCHENTRY searchTable[SEARCH_TABLE_LENGTH];


typedef struct _magsum_index {
  double magsum;
  int index;
} MAGSUM_INDEX,*PMAGSUM_INDEX;


typedef struct maskorder {
  int maskIndex;
  double MAG_ISO;
} MASKORDER,*PMASKORDER;

typedef struct knownMask {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
} KNOWNMASK,*PKNOWNMASK;

KNOWNMASK knownMaskTable [] = {
  {"a",10404,0},
  {"b",16412,0},
  {"b",16412,1},
  {"b",31829,0},
  {"b",46256,0},
  {"br",620,1},
  {"i",25055,0},
  {"i",31013,0},
  {"i",31090,0},
  {"i",49523,1},
  {"ir",11831,0},
  {"mc",48,0},
  {"mc",308,0},
  {"mc",309,0},
  {"mc",5077,1},
  {"mc",8175,1},
  {"mc",8621,0},
  {"mc",35417,1},
  {"rb",2046,0},
  {"rb",6125,0},
};

int totalNumMask = sizeof(knownMaskTable)/sizeof(KNOWNMASK);

/* Given an angle and a radius, compute the limits */
void GetLimits(double binAngle,double radius,double *pMinX,double *pMaxX,double *pMinY,double *pMaxY)
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
  printf("binAngle %f radius %f, xval %f yval %f\n",binAngle,radius,xval,yval);
#endif
  return;
}
int MAG_ISOSectorCompare(const void *first, const void *second) 
{
  double MAG_ISOFirst = ((PINPUTIMAGE)first)->MAG_ISO;
  double ySecond = ((PINPUTIMAGE)second)->MAG_ISO;
  int maskSectorFirst = ((PINPUTIMAGE)first)->maskSector;
  int maskSectorSecond = ((PINPUTIMAGE)second)->maskSector;

  if (maskSectorFirst > maskSectorSecond) {
    return(1);
  } else if (maskSectorFirst < maskSectorSecond) {
    return(-1);
  } else {
    if (MAG_ISOFirst > ySecond) {
      return(1);
    } else if (MAG_ISOFirst < ySecond) {
      return(-1);
    } else {
      return(0);
    }
  }
}
int maskOrderCompare(const void *first, const void *second) 
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
int maskBrightnessCompare(const void *first, const void *second) 
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

int searchCompare(const void *first, const void *second) 
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
int magsumCompare(const void *first, const void *second) 
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
int GetMask(int initialMask,int *reassignTable,int *reassignStack,int maskCount,char *fileroot)
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
      printf("ERROR: GetMask stack overflow for %s\n",fileroot);
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
/* Returns the number of valid masks and the noise floor which is the point in the histogram where the correlations exceed the MAX_MASKS */
int CountMasks(int tablePeak, int correlationTablePixels,int * correlation_histogram,int *correlation_table, int* correlation_scratchA, int* correlation_scratchB,int* noiseFloor,char *fileroot) {
  int result = -1;
  int histogram_scratch[SEARCH_TABLE_LENGTH+5];
  int peak_histogram[SEARCH_TABLE_LENGTH+5];
  int stack_table[SEARCH_TABLE_LENGTH+5];
  int stackIndex = 0;
  int currentPeak;
  int peakIndex;
  int xIndex;
  int yIndex;
  int xIndex2;
  int yIndex2;
  int localPeak;
  int localIndex;
  int testPeak;
  int testIndex;
  int iterationCounter = 0;
  int iterationCounter2 = 0;
  int iterationCounter2max = 0;
  int totalPoints = 0;
  int totalHighestPeaks = 0;

  *noiseFloor = -1;
  memcpy(correlation_scratchA,correlation_table,correlationTablePixels*sizeof(int));
  memset(correlation_scratchB,0,correlationTablePixels*sizeof(int));
  memcpy(histogram_scratch,correlation_histogram,sizeof(histogram_scratch));
  memset(peak_histogram,0,sizeof(peak_histogram));
  while (1) {
    iterationCounter++;
    if (iterationCounter > 1000000) {
      printf("ERROR line %4d iterationCounter overflow\n",__LINE__);
      exit(-1);
    }
    /* Get the current peak */
    for (currentPeak = tablePeak; currentPeak > 0; currentPeak--) {
      if (histogram_scratch[currentPeak] != 0) {
        break;
      }
    }
    if (currentPeak <= 0) {
      *noiseFloor = currentPeak;
      break;
    }
    if ((peak_histogram[currentPeak] >= MAX_MASKS) || (totalHighestPeaks >= (MAX_MASKS-1))) {
      *noiseFloor = currentPeak;
      break;
    }


    /* We have the highest peak */
    peak_histogram[currentPeak]++;
    totalHighestPeaks++;
    /* Find this peak in the histogram table */
    for (peakIndex = 0; peakIndex < correlationTablePixels; peakIndex++) {
      if (correlation_scratchA[peakIndex] == currentPeak) {
        break;
      }
    }
    if (peakIndex >= correlationTablePixels) {
      printf("ERROR: line %4d inconsistent histogram table\n");
      exit(-1);
    }
    localPeak = currentPeak;
    localIndex = peakIndex;
    histogram_scratch[localPeak]--;
    correlation_scratchA[localIndex] = -localPeak;
    correlation_scratchB[localIndex] = iterationCounter;
#if 0
    xIndex = (localIndex % (2*(MAX_PIXEL_RADIUS+1))) - MAX_PIXEL_RADIUS;
    yIndex = localIndex /((2*(MAX_PIXEL_RADIUS+1))) - MAX_PIXEL_RADIUS;
    printf("line %4d currentPeak %4d at index %10d x %4d  y %4d totalPoints %10d totalHighestPeaks %3d\n",__LINE__,localPeak,localIndex,xIndex,yIndex,totalPoints,totalHighestPeaks);
#endif
    totalPoints++;
    iterationCounter2 = 0;
    while (1) {
      iterationCounter2++;
      if (iterationCounter2 > 1000000) {
        printf("ERROR line %4d iterationCounter2  overflow\n",__LINE__);
        exit(-1);
      }
      if (iterationCounter2 > iterationCounter2max) {
        iterationCounter2max = iterationCounter2;
      }

      xIndex = (localIndex % (2*(MAX_PIXEL_RADIUS+1))) - MAX_PIXEL_RADIUS;
      yIndex = localIndex /((2*(MAX_PIXEL_RADIUS+1))) - MAX_PIXEL_RADIUS;
      for (xIndex2 = xIndex-1; xIndex2 <= xIndex+1; xIndex2++) {
        if ((xIndex2 < -MAX_PIXEL_RADIUS) ||
            (xIndex2 >  MAX_PIXEL_RADIUS)) {
          continue;
        }
        for (yIndex2 = yIndex-1; yIndex2 <= yIndex+1; yIndex2++) {
          if ((yIndex2 < -MAX_PIXEL_RADIUS) ||
              (yIndex2 >  MAX_PIXEL_RADIUS)) {
            continue;
          }
          testIndex = (xIndex2+MAX_PIXEL_RADIUS) + ((yIndex2+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
          testPeak = correlation_scratchA[testIndex];
          if (testPeak <= 0) {
            continue;
          }
          if (testPeak > localPeak) {
            continue;
          }
          histogram_scratch[testPeak]--;
          totalPoints++;
#if 0
          printf("line %4d                 localPeak  %4d testPeak %4d at index %10d x %4d  y %4d totalPoints %10d\n",__LINE__,localPeak,testPeak,testIndex,xIndex2,yIndex2,totalPoints);
#endif
          correlation_scratchA[testIndex] = -testPeak;
          correlation_scratchB[testIndex] = iterationCounter;
          if (stackIndex >= SEARCH_TABLE_LENGTH) {
            printf("ERROR: line %4d stackPtr exceeded\n",__LINE__);
            exit(-1);
          }
          stack_table[stackIndex] = testIndex;
          stackIndex++;
        }
      }
      /* Pop the last item off the stack and process it */
      if (stackIndex <= 0) {
        localIndex = -1;
        break;
      } else {
        stackIndex--;
        localIndex = stack_table[stackIndex];
#if 0
        xIndex = (localIndex % (2*(MAX_PIXEL_RADIUS+1))) - MAX_PIXEL_RADIUS;
        yIndex = localIndex /((2*(MAX_PIXEL_RADIUS+1))) - MAX_PIXEL_RADIUS;
        localPeak = correlation_scratchA[localIndex];
        printf("line %4d currentPeak %4d localPeak %4d at index %10d x %4d  y %4d totalPoints %10d\n",__LINE__,currentPeak,localPeak,localIndex,xIndex,yIndex,totalPoints);
#endif
        localPeak = abs(correlation_scratchA[localIndex]);


      }
        
      
    }
#if 0
    printf("line %4d iterationCounter %d iterationCounter2max %d result %d totalHighestPeaks %2d\n",__LINE__,iterationCounter,iterationCounter2max,result,totalHighestPeaks);
#endif
  } /* Outermost loop */
#if 0
  /* Clear out all matched and unmatched objects below the noise floor */
  for (localIndex = 0; localIndex < correlationTablePixels; localIndex++) {
    if (correlation_scratchA[localIndex] >= *noiseFloor) {
      correlation_scratchB[localIndex] = 0;
    }
  }
#endif
  /* At this point, count the total number of peaks */
  result = 0;
  for (currentPeak = tablePeak; currentPeak > *noiseFloor; currentPeak--) {
    if (peak_histogram[currentPeak] >= MAX_MASKS) {
      break;
    }
    result += peak_histogram[currentPeak];
  }
#if 0
  printf("line %4d CountMasks currentPeak %4d  CountMasks %4d totalPoints %4d noiseFloor %4d\n",__LINE__,currentPeak,result,totalPoints,*noiseFloor);
#endif

  return(result);
}


void FindMasks(int count_limit,int *correlation_table2,int *correlation_table3,PMULTMASK *pMaskTable,int* pTotalMasks,char *fileroot) 
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
  PMULTMASK pMask;

  /* First assign a mask to each pixel over the convolution limit */

  for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
    for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
      correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
      if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
        printf("ERROR: correlationIndex %d exceeds limits (1) %d\n",correlationIndex,maxCorrelationIndex);
        exit(-1);
      }
      if ((correlation_table2[correlationIndex]) > count_limit) {
        correlation_table3[correlationIndex] = ++maskIndex;
#if 0
        printf("Mask %d at xIndex %d yIndex %d\n",maskIndex,xIndex,yIndex);
#endif
      }
      
    }
  }
  maskCount = maskIndex;
#if 0
  printf("Number of masks %d for %s\n",maskIndex,fileroot);
#endif
  if (maskIndex <= 0) {
    return;
  }
  reassignTable = (int *)calloc(maskCount+1,sizeof(int));
  if (reassignTable == NULL) {
    printf("ERROR: failed to allocate reassignTable of size %d\n",maskCount+1);
  }
  reassignStack = (int *)calloc(maskCount+1,sizeof(int));
  if (reassignStack == NULL) {
    printf("ERROR: failed to allocate reassignStack of size %d\n",maskCount+1);
  }
  /* Now examine all the neighbors of a hot pixel and assign a mask number equal to the lowest hot neighbor */
  for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
    for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
      correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
      if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
        printf("ERROR: correlationIndex %d exceeds limits (2) %d\n",correlationIndex,maxCorrelationIndex);
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
            printf("ERROR: correlationIndex %d exceeds limits (3) %d\n",correlationIndex,maxCorrelationIndex);
            exit(-1);
          }
          if (correlationIndex == correlationIndex2) {
            continue;
          }
          maskIndex1 = GetMask(correlation_table3[correlationIndex],reassignTable,reassignStack,maskCount,fileroot);
          maskIndex2 = GetMask(correlation_table3[correlationIndex2],reassignTable,reassignStack,maskCount,fileroot);
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
  maskTable = (PMULTMASK)calloc(maskCount,sizeof(MULTMASK));
  if (maskTable == NULL) {
    printf("ERROR: failed to allocate maskTable of size %d\n",maskCount);
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
      nextMask = GetMask(maskIndex,reassignTable,reassignStack,maskCount,fileroot);
      
      if (nextMask != nextMask2) {
        workDone = 1;
        if (whileIteration > maskCount) {
          printf("ERROR: too many iterations in FindMasks for %s\n",fileroot);
          for (maskIndex = 1; maskIndex <= maskCount; maskIndex++) {
            printf("Iteration %d maskIndex %d reassignTable %d\n",whileIteration,maskIndex,reassignTable[maskIndex]);
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
        printf("ERROR: correlationIndex %d exceeds limits (4) %d\n",correlationIndex,maxCorrelationIndex);
        exit(-1);
      }
      if (correlation_table3[correlationIndex] == 0) {
        continue;
      }
      maskIndex1 = GetMask(correlation_table3[correlationIndex],reassignTable,reassignStack,maskCount,fileroot);
      for (maskIndex2 = 0; maskIndex2 < finalMaskCount; maskIndex2++) {
        pMask = &maskTable[maskIndex2];
        if (pMask->maskIndex == maskIndex1) {
          break;
        }
      }
      if (maskIndex2 == finalMaskCount) {
        printf("ERROR: maskIndex %d %d not found for %s\n",correlation_table3[correlationIndex],maskIndex1,fileroot);
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
      printf("ERROR: Mask %d has no entries\n",pMask->maskArea);
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
  printf("FinalMaskCount is %d for %s\n",finalMaskCount,fileroot);
#endif

  
  free(reassignTable);
  free(reassignStack);
  *pTotalMasks = finalMaskCount;
  *pMaskTable = maskTable;
  return;
}

int GetItbg(double X_IMAGE,double Y_IMAGE,int mosaicWidth,int mosaicHeight,int nxbg,int nybg)
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
double Factorial(int iteration,int *pOverflowFlag) {
  double factorial = 1.;
  int index;
  int overflowFlag = *pOverflowFlag;
  if (iteration == 0) {
    return(factorial);
  }
  for (index = 1; index <= iteration;index++) {
    factorial = factorial*index;
    if ((overflowFlag == 0) && (isnormal(factorial) == 0)) {
      overflowFlag = 1;
    }
  }
  *pOverflowFlag = overflowFlag;
  return(factorial);
}
int MatchIndexCompare(const void *first, const void *second) 
{
  int matchIndexFirst = *((int *)first);
  int matchIndexSecond = *((int *)second);
  if (matchIndexFirst > matchIndexSecond) {
    return(1);
  } else if (matchIndexFirst < matchIndexSecond) {
    return(-1);
  } else {
    return(0);
  }
}

int main(int argc,char *argv[])
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
  int index;
  int searchIndex;
  int insertIndex;
  int seriesIndex = -1;
  int searchTableSize = 0;
  double matchRatio2;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int candidateCount = 0;
  int outputCount = 0;
  int populatedBinCount = 0;
  int psfsaturatedCount = 0;
  int neighborIndex;
  PMASKSECTOR maskSectorTable = NULL;
  PMASKSECTOR pMaskSector;
  PMASKSECTOR pMaskSector2;
  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  PINPUTIMAGE sextractor_table = NULL;
  PINPUTIMAGE pInputImage;
  size_t sextractor_nrecs = 0;
  size_t vectorSize;
  size_t vectorIndex;
  int sextractorIndex;
  int sextractorIndex1;
  int sextractorIndex2;
  int baseIndex;
  PSTARIMAGE image_table;
  PSTARIMAGE pSextractor1 = NULL;
  PSTARIMAGE pSextractor2 = NULL;
  PSTARIMAGE pSextractor3 = NULL;
  PSTARIMAGE pSextractorBright = NULL;
  PSTARIMAGE pSextractorDim = NULL;
  int maxBinCountOverall = 0;

  int correlation_sum0 = 0;
  int correlation_sum1 = 0;
  double pmf_sum0 = 0.0;
  double bin_sum0 = 0.0;
  double bin_ratio;
  int correlation_histogram1[SEARCH_TABLE_LENGTH+5];
  int correlation_histogram2[SEARCH_TABLE_LENGTH+5];
  int *correlation_table1 = NULL;
  int *correlation_table2 = NULL;
  int *correlation_table3 = NULL;
  int *correlation_table4 = NULL;
  int *correlation_scratch1A = NULL;
  int *correlation_scratch1B = NULL;
  int *correlation_scratch2A = NULL;
  int *correlation_scratch2B = NULL;
  int correlationTablePixels;

  int table1Points;
  int table1Total;
  int table1Peak;
  int table1OldTotal;  
  double table1Ratio;
  int table1NoiseFloor;
  int table1MaskCount;
  int histogramIndex;
  int table1MaxPoints;
  int table1MaxTotal;
  int table1MaxPeak;
  double table1MaxRatio;
  double table1_dimlimit;

  int table2Points;
  int table2Total;
  int table2Peak;
  double table2Ratio;
  int table2NoiseFloor;
  int table2MaskCount;
  


  int minX;
  int maxX;
  int minY;
  int maxY;
  double binAngleSize;
  int angleBinNumber;
  int averageBin;
  double sigmaSum = 0;
  double sigmaSqr = 0;
  double binMean;
  double binStd;
  double binValue;
  double binAngle;

  int verbose = 0;
  int includeSaturated = 0;
  int doPlots1 = 0;
  int doPlots2 = 0;
  int emulsionPlots = 0;
  time_t startTime;
  time_t curTime;

  char series[MAX_SERIES_STRING];
  int plateNumber;
  double excess;
  double maxExcess = -1000000.0;
  int maxIndex1 = 0;
  int maxIndex2 = 0;
  int index1;
  int index2;

  double minDec;
  double maxDec;
  double minRa;
  double maxRa;
  double drad;
  double xValue;
  double yValue;
  double angle;

  int totSqrBins; /* Total number of square bins */
  int numSqrBins; /* Number of sqare bins per axis */
  int xSqrBin;
  int ySqrBin;
  int iSqrBin;
  int matchIndex;
  int binIndex;
  int maxSqrBinValue;
  int maxSqrBinIndex;
  int aveBinValue;
  int startIndex;
  int endIndex;
  double brightestMagIso;
  int brightestMagIsoIndex;
  PMAGSUM_INDEX pMagsumIndexList = NULL;
  PMAGSUM_INDEX pMagsumIndex;
  int plotCount = 0;
  double arcsecPerPixel;
  double degreesPerPixel;
 
  double * vector = NULL;
  /* This is the cutoff for the ratio (FLUX_MAX/THRESHOLD) */
  double analysis_threshold = USE_ANALYSIS_THRESHOLD;
  double analysis_threshold2;  /* Computed for histogram */
  double MAG_ISO_crit;
  double MAG_ISO_bright;

  double MAG_ISO_dimmest;
  double MAG_ISO_brightest;
  int MAG_ISO_brightindex;
  int MAG_ISO_index;
  double MAG_ISO_brightlimit;
  double MAG_ISO_dimlimit;


  double MAG_ISO_match;
  int MAG_ISO_matchindex;
  int *MAG_ISO_countmatch1; 
  int *MAG_ISO_countmatch2; 
  int countmatch1;
  int countmatch2;

  int *MAG_ISO_count1;   /* Correlation between bright and dim */
  int *MAG_ISO_count2;  /* Simple histogram */
  int MAG_ISO_size;
  int MAG_ISO_alloc;
  int count;
  int numPeaks = 0;
  int abovePeak;
  int initAbovePeak;
  int finalAbovePeak;
  int solutionNumber = 0;
  int quality;
  int wedgeFlag = -1;
  int xIndex;
  int yIndex;
  int xIndex2;
  int yIndex2;
  int correlationIndex; /* Index into correlation table (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1))) */
  int correlationIndex2; /* Index into correlation table (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1))) */
  int correlationCount; /* Histogram index */
  int convolutionCount; /* Pixel count in the correlation area */
  int maxCorrelationIndex = 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1);
  int NUMBER;
  PSEARCHENTRY pSearchEntry1;
  PSEARCHENTRY pSearchEntry2;
  PMULTMASK maskTable = NULL;
  PMULTMASK pMask;
  int totalMasks = 0;
  int debugPrint = 0;
  int maskIndex;
  int maskIndex2;
  int selected = 0;
  int binning = 1;
  PMASKORDER maskOrderTable = NULL;
  PMASKORDER pMaskOrder;
  PMASKORDER pMaskOrder2;


  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  char *mysqlhost;
  char *username;
  char *password;
  int numExposures;
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
  int clip_factor;
  int clip_count;
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
  
  int nxs;
  int nys;
  int ixs;
  int iys;
  int ixs2;
  int iys2;
  int totalMaskSectors;
  int maskSector;
  int maskSector2;
  int oldMaskSector;
  double oldMAG_ISO;
  int minSectorCount;
  int maxSectorCount; 
  int minUnsaturatedCount;
  int maxUnsaturatedCount;
  int BFLAGSMASK = (~((1 << FILTER_BFLAG_NEIGHBORS) | (1 << FILTER_BFLAG_BLEND)));
  double pmf2;
  double meanArrival;
  int intMeanArrival;
  double factorial = 1.;
  double prefix;
  double exponential;
  double matchRatio;
  int overflowFlag;
  int pmf2_valid;
  int overflowPrint;
  int sequence;
  PKNOWNMASK pKnownMask;
  int knownMaskIndex;

  outfile[0] = 0;
  backgroundfile[0] = 0;
#ifdef BACKGROUND_PLOTS
  regionfile[0] = 0;
#endif /* BACKGROUND_PLOTS */

  time(&startTime);

  memset(backgroundVector,0,sizeof(backgroundVector));
  memset(sextractorBinIndex,0,sizeof(sextractorBinIndex));
  memset(ixVector,0,sizeof(ixVector));
  memset(iyVector,0,sizeof(iyVector));
  memset(localVector,0,sizeof(localVector));

  memset(background_med,0,sizeof(background_med));
  memset(background_rms,0,sizeof(background_rms));
  memset(background_clip,0,sizeof(background_clip));



  memset(gmtHandle,0,sizeof(gmtHandle));


#ifdef LOS_DEBUG 
  printf("ERROR: LOS_DEBUG is set \n");
#endif /* LOS_DEBUG */
#ifdef DEBUG_TRUNCATION
  printf("ERROR: DEBUG_TRUNCATION is set \n");
#endif /* DEBUG_TRUNCATION */
#ifdef DEBUG_MAG_ISO
  printf("ERROR: DEBUG_MAG_ISO is set \n");
#endif /* DEBUG_MAG_ISO */
#ifdef DEBUG_SINGLE_SECTOR
  printf("ERROR: DEBUG_SINGLE_SECTOR is set \n");
#endif /* DEBUG_SINGLE_SECTOR */

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
      printf("ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'i': /* sextractor file name */
        case 'I':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sextractor_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'r': /* file root */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(fileroot,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
            strcpy(fileroot2,fileroot);
            charPtr = strstr(fileroot2,"_");
            if (charPtr != NULL) {
              *charPtr = 0;
            }
          }
          break;

        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic width %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* mosaic height */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'b': /* binning */
        case 'B':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&binning);
            if (nvals != 1) {
              printf("ERROR: Unable to decode binning %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 's': /* arcsec/pixel */
        case 'S':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&arcsecPerPixel);
            if (nvals != 1) {
              printf("ERROR: Can not decode number of arcsec per pixel\n");
              errorFlag = 1;
            }
          }
          break;
        case 't': /* analysis_threshold */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&analysis_threshold);
            if (nvals != 1) {
              printf("ERROR: Can not decode analysis threshold\n");
              errorFlag = 1;
            }
          }
          break;
#if 0 
        case 'u': /* include saturated PSF values */
        case 'U':
          includeSaturated = 1;
          break;
#endif
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
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(backgroundfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
#ifdef BACKGROUND_PLOTS
            strncpy(regionfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
#else /* BACKGROUND_PLOTS */
            ++argv;
#endif /* BACKGROUND_PLOTS */

          }
          break;

          break;

        case 'p': /* do plotting */
        case 'P':
          doPlots1 = 1;
          break;

        case 'q': /* do emulsion plots */
        case 'Q':
          emulsionPlots = 1;
          break;

        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
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
  degreesPerPixel = arcsecPerPixel/3600.0;

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
  if (ParseFilename(sextractor_name,series,&plateNumber,&mosaicNumber,&tempbinning,&rotation) == 0) {
    printf("ERROR: Failed to parse filename %s\n",sextractor_name);
    errorFlag = 1;
  }


  

  sextractor_handle = Open(sextractor_name,"r");
  if (sextractor_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the sextractor file %s\n",sextractor_name);
  } else {
    if (verbose) {
      printf("Found sextractor file %s\n",sextractor_name);
    }
  }

   
  if (doPlots1) {
    for (gmtIndex = 0; gmtIndex < NUM_GMT_FILES; gmtIndex++) {
      
      strcpy(gmtname[gmtIndex],outfile);
      curPtr = gmtname[gmtIndex];
      while ((slashPtr = strstr(curPtr,"/")) != NULL) {
        curPtr = slashPtr+1;
      }
      if (curPtr != NULL) {
        *curPtr = 0;
      }
      strcat(gmtname[gmtIndex],fileroot);
      sprintf(tmpname,"_close_gmt%d.txt",gmtIndex);
      strcat(gmtname[gmtIndex],tmpname);
    
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
    printf("Usage: search_close2\n");
    printf("                 -b  <binning> Set to 16 for regions on thumbnail mosaics\n");
    printf("                 -d  <background file> <region file>search for calibration squares\n");
    printf("                 -h  <mosaic height>\n");
    printf("                 -i <sextractor file> \n");
    printf("                 -o  <output file>\n");
    printf("                 -p  write GMT plotting files\n");
    printf("                 -q  write emulsion plotting files\n");
    printf("                 -r <file root>\n");
    printf("                 -s <arcsec per pixel>\n");
    printf("                 -t <analysis_threshold>, cutoff ratio of FLUX_MAX/THRESHOLD\n");
#if 0
    printf("                 -u  include saturated PSFs\n");
#endif
    printf("                 -v  verbose\n");
    printf("                 -w  <mosaic width>\n");

    return(-1);
  }
  nx = XDmagBins(mosaicWidth,mosaicHeight);
  ny = YDmagBins(mosaicWidth,mosaicHeight);

  nxs = mosaicWidth/SECTOR_PIXEL_WIDTH;
  if (nxs < 1) {
    nxs = 1;
  }
  nys = mosaicHeight/SECTOR_PIXEL_WIDTH;
  if (nys < 1) {
    nys = 1;
  }
  totalMaskSectors = nxs*nys;
  maskSectorTable = (PMASKSECTOR)calloc(totalMaskSectors+1,sizeof(MASKSECTOR));
  if (maskSectorTable == NULL) {
    printf("ERROR: failed to allocate the maskSectorStable of size %d\n",totalMaskSectors);
    exit(-1);
  }

  if (verbose) {
    printf("search_close2 of %s %s \n Input Filename %s\n Output Filename %s mosaic width: %d mosaic height: %d scale: %f arcsec/pixel:",
           __DATE__,__TIME__,sextractor_name,outfile,
           mosaicWidth,
           mosaicHeight,
           arcsecPerPixel);
    printf(" threshold: %f",analysis_threshold);
    printf("\n");
    printf("search_close2 BFLAGSMASK is %d 0x%x ~BFLAGSMASK is %d 0x%x\n",BFLAGSMASK,BFLAGSMASK,~BFLAGSMASK,~BFLAGSMASK);

  }
  
  if ((totalMaskSectors) > MAX_SECTOR_COUNT) {
    printf("ERROR: MAX_SECTOR_COUNT %d is too large for %s\n",totalMaskSectors,fileroot);
    exit(-1);
  }


  /* Connect to the database */
  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    printf("ERROR: filter_multiple  DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    printf("ERROR: filter_multiple  DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    printf("ERROR: filter_multiple  DASCH_PASSWORD is not defined\n");
    return(-1);
  }
                   
  mysql_init(pConnection);


  if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pConnection)) {
      printf("ERROR: filter_multiple  MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
    }
    return(-1);
  }

#ifndef DEBUG_TRUNCATION
  /* Clear the multiple mask status */
  SetMosaicFitWCS(pConnection,"NoMultipleMask",fileroot,0,0,readFlag,&FitWCS);
  SetMosaicFitWCS(pConnection,"HaveMultipleMask",fileroot,0,0,readFlag,&FitWCS);
  
  sprintf(queryString,"DELETE FROM masks2 where series = '%s' and plateNumber = %d and mosaicNumber = %d;",
          series,
          plateNumber,
          mosaicNumber);
  if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
    printf("ERROR: search_close2 MAX_QUERY_STRING exceeded\n");
    exit(-1);
  }
  res = ExecuteQuery(pConnection,queryString);
  if (res) {
    exit(-1);
  }
#endif /* DEBUG_TRUNCATION */
  sprintf(queryString,"DELETE FROM masks3 where series = '%s' and plateNumber = %d and mosaicNumber = %d;",
          series,
          plateNumber,
          mosaicNumber);
  if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
    printf("ERROR: search_close2 MAX_QUERY_STRING exceeded\n");
    exit(-1);
  }
  res = ExecuteQuery(pConnection,queryString);
  if (res) {
    exit(-1);
  }


  sextractor_header = table_header(sextractor_handle,TABLE_PARSE);
  if (sextractor_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",sextractor_name);
    return(-1);
  }

  sextractor_table = table_loadva(sextractor_handle,
                                  &sextractor_header,
                                  NULL, /* hbase */
                                  NULL, /* rows */
                                  NULL,
                                  sizeof(INPUTIMAGE),
                                  &sextractor_nrecs,
                                  TblInt,"NUMBER",TblOff(PINPUTIMAGE,NUMBER),
                                  TblInt,"BFLAGS",TblOff(PINPUTIMAGE,BFLAGS),
                                  TblInt,"ISO5",TblOff(PINPUTIMAGE,ISO5),
                                  TblInt,"ISO4",TblOff(PINPUTIMAGE,ISO4),
                                  TblDbl,"MAG_ISO",TblOff(PINPUTIMAGE,MAG_ISO),
                                  TblDbl,"BACKGROUND",TblOff(PINPUTIMAGE,BACKGROUND),
                                  TblDbl,"FLUX_MAX",TblOff(PINPUTIMAGE,FLUX_MAX),
                                  TblDbl,"THRESHOLD",TblOff(PINPUTIMAGE,THRESHOLD),
                                  TblDbl,"FWHM_IMAGE",TblOff(PINPUTIMAGE,FWHM_IMAGE),
                                  TblDbl,"X_IMAGE",TblOff(PINPUTIMAGE,X_IMAGE),
                                  TblDbl,"Y_IMAGE",TblOff(PINPUTIMAGE,Y_IMAGE),
                                  0,"end",0);
  if (sextractor_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",sextractor_name);
    return(-1);
  }
  if (verbose) {
    printf("read %d records for %s\n",sextractor_nrecs,sextractor_name);
  }


  /* Calculate the spatial bin and maskSector for each record */
  for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
    pInputImage = &sextractor_table[sextractorIndex];
    pInputImage->spatial_bin =  CalculateBin(mosaicWidth,mosaicHeight,pInputImage->X_IMAGE,pInputImage->Y_IMAGE,&edgeDist);
    ixs = (pInputImage->X_IMAGE * nxs) /(1.0 * mosaicWidth);
    iys = (pInputImage->Y_IMAGE * nys) /(1.0 * mosaicHeight);
    if (ixs == nxs) {
      ixs--;
    }
    if (iys == nys) {
      iys--;
    }
    pInputImage->maskSector =  ixs + (nxs*iys);
    if ((pInputImage->maskSector < 0) || (pInputImage->maskSector >= totalMaskSectors)) {
      printf("ERROR: line %d illegal maskSector value %d\n",__LINE__,pInputImage->maskSector);
      exit(-1);
    }
    pMaskSector = &maskSectorTable[pInputImage->maskSector];
    pMaskSector->sectorCount++;
    pMaskSector = NULL;
  }

  if (backgroundStudy) {
    nxbg = (mosaicWidth+BACKGROUND_BIN_SIZE-1)/BACKGROUND_BIN_SIZE;
    nybg = (mosaicHeight+BACKGROUND_BIN_SIZE-1)/BACKGROUND_BIN_SIZE;
    ntbg = nxbg*nybg;
    pBinMap = (PBACKGROUNDBIN)calloc(ntbg,sizeof(BACKGROUNDBIN));
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
      backgroundVector[spatial_bin] = (double *)calloc(sextractor_nrecs,sizeof(double));
      if (backgroundVector[spatial_bin] == NULL) {
        printf("ERROR Failed to allocate the background backgroundVector\n");
        exit(-1);
      }
    }
    /* Compute the average background in each bin */
#ifdef BACKGROUND_PLOTS 
    strcpy(backgroundmap,regionfile);
    slashPtr = strrchr(backgroundmap,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = backgroundmap;
    }
    strcpy(slashPtr,"map.db");
    backgroundmapHandle = fopen(backgroundmap,"wt");
    if (backgroundmapHandle == NULL) {
      printf("ERROR: failed to open %s\n",backgroundmap);
      exit(-1);
    }
    fprintf(backgroundmapHandle,"X_IMAGE\tY_IMAGE\tBACKGROUND\tFILTERED_BACKGROUND\tcount\texpansionRadius\n");
    fprintf(backgroundmapHandle,"-------\t-------\t----------\t-----------------\t-----\t---------------\n");
#endif /* BACKGROUND_PLOTS */


    for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
      pInputImage = &sextractor_table[sextractorIndex];
      itbg = GetItbg(pInputImage->X_IMAGE,pInputImage->Y_IMAGE,mosaicWidth,mosaicHeight,nxbg,nybg);

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
#if 1
#ifdef BACKGROUND_PLOTS
          fprintf(backgroundmapHandle,"%.1f\t%.1f\t%.1f\t%.1f\t%d\t%.1f\n",pBackBin->X_IMAGE,pBackBin->Y_IMAGE,pBackBin->BACKGROUND,pBackBin->FILTERED_BACKGROUND,curCount,pBackBin->expansionRadius);
#endif /* BACKGROUND_PLOTS */
#endif
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

      fprintf(backgroundmapHandle,"%.1f\t%.1f\t%.1f\t%.1f\t%d\t%.1f\n",pBackBin->X_IMAGE,pBackBin->Y_IMAGE,pBackBin->BACKGROUND,pBackBin->FILTERED_BACKGROUND,pBackBin->count,pBackBin->expansionRadius);
#endif /* BACKGROUND_PLOTS */
#endif
    }

#ifdef BACKGROUND_PLOTS
    strcpy(imagemap,regionfile);
    slashPtr = strrchr(imagemap,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = imagemap;
    }
    strcpy(slashPtr,"image.db");
    imagemapHandle = fopen(imagemap,"wt");
    if (imagemapHandle == NULL) {
      printf("ERROR: failed to open %s\n",imagemap);
      exit(-1);
    }
    fprintf(imagemapHandle,"X_IMAGE\tY_IMAGE\tBACKGROUND\tFILTERED_BACKGROUND\tcount\tbackground_ratio\n");
    fprintf(imagemapHandle,"-------\t-------\t----------\t-----------------\t-----\t----------\n");
#endif /* BACKGROUND_PLOTS */
    for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
      pInputImage = &sextractor_table[sextractorIndex];
      itbg = GetItbg(pInputImage->X_IMAGE,pInputImage->Y_IMAGE,mosaicWidth,mosaicHeight,nxbg,nybg);
      pBackBin = &pBinMap[itbg];
      pInputImage->FILTERED_BACKGROUND = pInputImage->BACKGROUND - pBackBin->FILTERED_BACKGROUND;
      curVector = backgroundVector[pInputImage->spatial_bin];
      curVector[sextractorBinIndex[pInputImage->spatial_bin]] = pInputImage->FILTERED_BACKGROUND;
      sextractorBinIndex[pInputImage->spatial_bin]++;
#if 0
#ifdef BACKGROUND_PLOTS
      fprintf(imagemapHandle,"%.1f\t%.1f\t%.1f\t%.1f\t%d\n",pInputImage->X_IMAGE,pInputImage->Y_IMAGE,pInputImage->BACKGROUND,pInputImage->FILTERED_BACKGROUND,pBackBin->count);
#endif /* BACKGROUND_PLOTS */
#endif
    }
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      curVector = backgroundVector[pInputImage->spatial_bin];
      background_clip[spatial_bin] = CalcMedianAndRMS(sextractorBinIndex[spatial_bin],10,curVector,&background_med[spatial_bin],&background_rms[spatial_bin],1,BACKGROUND_CLIP,0);
      total_background_clip += background_clip[spatial_bin];
      if (background_clip[spatial_bin] == 0) {
        printf("ERRORD: Insufficient background stars in spatial bin %d for %s\n",spatial_bin,sextractor_name);
        background_rms[spatial_bin] = 0;
      } else if (background_rms[spatial_bin] == 0) {
        printf("ERRORD: Zero background rms in spatial bin %d for %s\n",spatial_bin,sextractor_name);
      } else {
        total_background_med += background_med[spatial_bin];
        total_background_rms += background_rms[spatial_bin];
      }
    }
    int clip_count = 0;
    if (total_background_clip > 0) {
      backgroundHandle = fopen(backgroundfile,"wt");
      if (backgroundHandle == NULL) {
        fprintf(stderr,"ERROR: Failed to open the output file %s\n",backgroundfile);
        exit(-1);
      } else {
        fprintf(backgroundHandle,"NUMBER\tAFLAGS\tbackground_ratio\n");
        fprintf(backgroundHandle,"------\t------\t----------------\n");
      }
#ifdef BACKGROUND_PLOTS
      regionHandle = fopen(regionfile,"wt");
      if (regionHandle == NULL) {
        fprintf(stderr,"ERROR: Failed to open the output file %s\n",regionfile);
        exit(-1);
      }

#endif /* BACKGROUND_PLOTS */
        
      
      for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
        pInputImage = &sextractor_table[sextractorIndex];
        itbg = GetItbg(pInputImage->X_IMAGE,pInputImage->Y_IMAGE,mosaicWidth,mosaicHeight,nxbg,nybg);  
        pBackBin = &pBinMap[itbg];
        if (background_rms[pInputImage->spatial_bin] > 0) {
          
          background_ratio = (pInputImage->FILTERED_BACKGROUND-background_med[pInputImage->spatial_bin])/background_rms[pInputImage->spatial_bin];
#ifdef BACKGROUND_PLOTS
          fprintf(imagemapHandle,"%.1f\t%.1f\t%.1f\t%.1f\t%d\t%.1f\n",pInputImage->X_IMAGE,pInputImage->Y_IMAGE,pInputImage->BACKGROUND,pInputImage->FILTERED_BACKGROUND,pBackBin->count,background_ratio);
#endif /* BACKGROUND_PLOTS */

          if (pInputImage->FILTERED_BACKGROUND > (background_med[pInputImage->spatial_bin] + (BACKGROUND_CLIP*background_rms[pInputImage->spatial_bin]))) {
            if (((pInputImage->FILTERED_BACKGROUND-background_med[pInputImage->spatial_bin])/background_rms[pInputImage->spatial_bin]) > BACKGROUND_LIMIT) {
              clip_factor = background_ratio;              
              fprintf(backgroundHandle,"%d\t%d\t%.2f\n",pInputImage->NUMBER,(1<<FILTER_AFLAG_BACKGROUND),background_ratio);
#ifdef BACKGROUND_PLOTS
              fprintf(regionHandle,"CIRCLE(%f,%f,5) # text = {%d}\n",pInputImage->X_IMAGE/binning,pInputImage->Y_IMAGE/binning,clip_factor);    
#endif /* BACKGROUND_PLOTS */
              clip_count++;
            } else {
#ifdef BACKGROUND_PLOTS
              fprintf(regionHandle,"CIRCLE(%f,%f,2) # color = blue\n",pInputImage->X_IMAGE/binning,pInputImage->Y_IMAGE/binning);    
#endif /* BACKGROUND_PLOTS */
            }
          } else {
#ifdef BACKGROUND_PLOTS
            fprintf(regionHandle,"CIRCLE(%f,%f,2) # color = red\n",pInputImage->X_IMAGE/binning,pInputImage->Y_IMAGE/binning);    
#endif /* BACKGROUND_PLOTS */
          }
        }
      }
        
      printf("background median %8.1f, background rms %6.1f clipped objects %5d plotted objects %6d bad bins %5d of %5d for %s%05d_%02d\n",total_background_med/MAX_SPATIAL_BINS,total_background_rms/MAX_SPATIAL_BINS,sextractor_nrecs-total_background_clip,clip_count,badBinCount,ntbg,series,plateNumber,mosaicNumber);

    }
#if 0
    exit(-1);
#endif
  }
  /* Sort the table in increasing maskSector and MAG_ISO values */
  oldMaskSector = -1;
  oldMAG_ISO = -9999999.;
  pMaskSector = NULL;
  qsort(sextractor_table,sextractor_nrecs,sizeof(INPUTIMAGE),MAG_ISOSectorCompare);
  /* Now copy the table over, setting the FILTER_BFLAG_PSFSATURATED flag */
  image_table = (PSTARIMAGE)calloc(sextractor_nrecs,sizeof(STARIMAGE));
  if (image_table == NULL) {
    printf("ERROR: failed to allocate image_table of size %d\n",sextractor_nrecs);
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
      if (((pInputImage->BFLAGS & BFLAGSMASK) == 0) &&
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
    pSextractor3->X_IMAGE = pInputImage->X_IMAGE;
    pSextractor3->Y_IMAGE = pInputImage->Y_IMAGE;
    if (sextractorIndex == 0) {
      MAG_ISO_dimmest = pInputImage->MAG_ISO;
      MAG_ISO_brightest = pInputImage->MAG_ISO;
    } else {
      if (MAG_ISO_dimmest < pInputImage->MAG_ISO) {
        MAG_ISO_dimmest = pInputImage->MAG_ISO;
      }
      if (MAG_ISO_brightest >  pInputImage->MAG_ISO) {
        MAG_ISO_dimmest = pInputImage->MAG_ISO;
      }
    }
    if (sextractorIndex == 0) {
      oldMAG_ISO = pSextractor3->MAG_ISO;
    }
    pSextractor3->maskSector = pInputImage->maskSector;
    if ((pSextractor3->maskSector < 0) || (pSextractor3->maskSector >= totalMaskSectors)) {
      printf("ERROR: line %d illegal maskSector value %d\n",__LINE__,pSextractor3->maskSector);
      exit(-1);
    }
    if (oldMaskSector != pSextractor3->maskSector) {
      if (oldMaskSector > pSextractor3->maskSector) {
        printf("ERROR maskSector sort error in line %d\n",__LINE__);
        exit(-1);
      }
      oldMaskSector = pSextractor3->maskSector;
      oldMAG_ISO = pSextractor3->MAG_ISO;
      pMaskSector = &maskSectorTable[oldMaskSector];
      pMaskSector->maskSector = oldMaskSector;
      pMaskSector->sextractorIndexLow = sextractorIndex;
    }
    pMaskSector->sextractorIndexHigh = sextractorIndex;
    if (pSextractor3->Y_IMAGE < oldMAG_ISO) {
      printf("ERROR Y_IMAGE sort error in line %d\n",__LINE__);
      exit(-1);
    }
    oldMAG_ISO = pSextractor3->MAG_ISO;    

    /* Set the FILTER_BFLAG_PSFSATURATED flag */
    if ((pInputImage->ISO4 > 0) &&
        (((1.0*pInputImage->ISO5)/(1.0*pInputImage->ISO4)) > PSFSATURATED_ISO) &&
        ((pInputImage->BACKGROUND + pInputImage->FLUX_MAX) > PSFSATURATED_FLUX)) {
#if 0
      printf("%d %f %f\n",pInputImage->NUMBER,(1.0*pInputImage->ISO5)/(1.0*pInputImage->ISO4),pInputImage->BACKGROUND + pInputImage->FLUX_MAX);
#endif
      pSextractor3->BFLAGS |= (1 << FILTER_BFLAG_PSFSATURATED);
      psfsaturatedCount++;
      pMaskSector2 = &maskSectorTable[pSextractor3->maskSector];
      pMaskSector2->unsaturatedCount--;
      pMaskSector2 = NULL;
    }
  }
  MAG_ISO_size = ((MAG_ISO_dimmest-MAG_ISO_brightest)/MAG_ISO_DELTA)+1.0;
  MAG_ISO_alloc = MAG_ISO_size*MAG_ISO_size;
  MAG_ISO_count1 = (int*)calloc(MAG_ISO_alloc,sizeof(int));
  
  if (MAG_ISO_count1 == NULL) {
    printf("ERROR: failed to allocate MAG_ISO_count1 of size %d\n",MAG_ISO_size);
    exit(-1);
  }

  MAG_ISO_count2 = (int*)calloc(MAG_ISO_size,sizeof(int));
  
  if (MAG_ISO_count2 == NULL) {
    printf("ERROR: failed to allocate MAG_ISO_count2 of size %d\n",MAG_ISO_size);
    exit(-1);
  }

  MAG_ISO_countmatch1 = (int*)calloc(MAG_ISO_size,sizeof(int));
  
  if (MAG_ISO_countmatch1 == NULL) {
    printf("ERROR: failed to allocate MAG_ISO_countmatch1 of size %d\n",MAG_ISO_size);
    exit(-1);
  }

  MAG_ISO_countmatch2 = (int*)calloc(MAG_ISO_size,sizeof(int));
  
  if (MAG_ISO_countmatch2 == NULL) {
    printf("ERROR: failed to allocate MAG_ISO_countmatch2 of size %d\n",MAG_ISO_size);
    exit(-1);
  }



  if (verbose) {
    printf("MAG_ISO_dimmest - MAG_ISO_brightest) = %f MAG_ISO_DELTA %f MAG_ISO_size %d for %s\n",MAG_ISO_dimmest-MAG_ISO_brightest,MAG_ISO_DELTA,MAG_ISO_size,fileroot);
  }

  maskSector = 0;
  pMaskSector = &maskSectorTable[maskSector];
  minSectorCount = pMaskSector->sectorCount;
  maxSectorCount = pMaskSector->sectorCount;
  minUnsaturatedCount = pMaskSector->unsaturatedCount+pMaskSector->sectorCount;
  maxUnsaturatedCount = pMaskSector->unsaturatedCount+pMaskSector->sectorCount;
#if 0
  printf("nxs %2d nys %2d\n",nxs,nys);
#endif
  for (maskSector = 0; maskSector < totalMaskSectors; maskSector++) {
    pMaskSector = &maskSectorTable[maskSector];
    pMaskSector->unsaturatedCount += pMaskSector->sectorCount;
    if (minSectorCount > pMaskSector->sectorCount) {
      minSectorCount = pMaskSector->sectorCount;
    }
    if (maxSectorCount < pMaskSector->sectorCount) {
      maxSectorCount = pMaskSector->sectorCount;
    }

    if (minUnsaturatedCount > pMaskSector->unsaturatedCount) {
      minUnsaturatedCount = pMaskSector->unsaturatedCount;
    }
    if (maxUnsaturatedCount < pMaskSector->unsaturatedCount) {
      maxUnsaturatedCount = pMaskSector->unsaturatedCount;
    }
    /* Now fill in the neighbors table */
    if (pMaskSector->neighborCount >= MAX_SECTOR_NEIGHBOR) {
      printf("ERROR MAX_SECTOR_NEIGHBOR exceeded in line %d\n",__LINE__);
      exit(-1);
    }
    pMaskSector->maskNeighbor[pMaskSector->neighborCount] = pMaskSector->maskSector;

    pMaskSector->neighborCount++;
    ixs2 = -1;
    iys2 = -1;
    ixs = pMaskSector->maskSector % nxs;
    iys = pMaskSector->maskSector / nxs;
#if 0
    printf("maskSector %2d ixs %2d iys %2d neigbor %d ixs2 %2d iys2 %2d maskSector2 %2d\n",maskSector,ixs,iys,pMaskSector->neighborCount,ixs2,iys2,pMaskSector->maskSector);
#endif
    for (ixs2 =  (ixs -1); ixs2 <= (ixs +1); ixs2++) {
#if 0
      printf("ixs2 %2d\n",ixs2);
#endif
      if ((ixs2 < 0) || ( ixs2 >= nxs)) {
        continue;
      }
      for (iys2 = (iys -1); iys2 <= (iys +1); iys2++) {
#if 0
        printf("       iys2 %2d\n",iys2);
#endif
        if ((iys2 < 0) || (iys2 >= nys)) {
          continue;
        }
        if ((ixs == ixs2) && (iys == iys2)) {
          continue;
        }
        maskSector2 = ixs2 + (nxs*iys2);
        if (pMaskSector->neighborCount >= MAX_SECTOR_NEIGHBOR) {
          printf("ERROR MAX_SECTOR_NEIGHBOR exceeded in line %d\n",__LINE__);
          exit(-1);
        }
        pMaskSector->maskNeighbor[pMaskSector->neighborCount] = maskSector2;
        pMaskSector->neighborCount++;
#if 0
        printf("maskSector %2d ixs %2d iys %2d neigbor %d ixs2 %2d iys2 %2d maskSector2 %2d\n",maskSector,ixs,iys,pMaskSector->neighborCount,ixs2,iys2,maskSector2);
#endif
      }
    } 
  }
  pMaskSector = NULL;
  if (verbose) {
    printf(" sections: %2d x %2d minSectorCount: %5d maxSectorCount %5d minUnsaturedCount %5d maxUnsaturedCount %5d for %s\n",nxs,nys,minSectorCount,maxSectorCount,minUnsaturatedCount,maxUnsaturatedCount,fileroot);
  }

  if (emulsionPlots != 0) {
    strcpy(emulsionName,outfile);
    curPtr = emulsionName;
    while ((slashPtr = strstr(curPtr,"/")) != NULL) {
      curPtr = slashPtr+1;
    }
    if (curPtr != NULL) {
      *curPtr = 0;
    }
    strcat(emulsionName,fileroot);
    sprintf(tmpname,"_close_plot.db",gmtIndex);
    strcat(emulsionName,tmpname);
    emulsionHandle = fopen(emulsionName,"wt");
    if (emulsionHandle == NULL) {
      printf("ERROR: failed to open the emulsion plot name %s\n",emulsionName);
    } else {
      fprintf(emulsionHandle,"ix\tiy\tcount\tseries\tplateNumber\tmosaicNumber\n");
      fprintf(emulsionHandle,"--\t--\t-----\t------\t-----------\t------------\n");
      for (ix = 0; ix < nx; ix++) {
        fprintf(emulsionHandle,"%d\t-1\t%d\t%s\t%d\t%d\n",ix,ixVector[ix],series,plateNumber,mosaicNumber);
        for (iy = 0; iy < ny; iy++) {
          if (ix == 0) {
            fprintf(emulsionHandle,"-1\t%d\t%d\t%s\t%d\t%d\n",iy,iyVector[iy],series,plateNumber,mosaicNumber);  
          }
          local_bin = ix + (nx*iy);
          fprintf(emulsionHandle,"%d\t%d\t%d\t%s\t%d\t%d\n",ix,iy,localVector[local_bin],series,plateNumber,mosaicNumber);
        }
      }         
      fclose(emulsionHandle);
    }
  }

  /* At this point, we are done with the sextractor_table.  Free memory */
  if (sextractor_table != NULL) {
    Free(sextractor_table);
    sextractor_table = NULL;  }
  if (sextractor_header != NULL) {
    table_hdrfree(sextractor_header);
    sextractor_header = NULL;
  }
  if (sextractor_handle != NULL) {
    Close(sextractor_handle);
    sextractor_handle = NULL;
  }
  correlationTablePixels = 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1);
  correlation_table1 = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_table1 == NULL) {
    printf("ERROR: failed to allocate correlation table 1 of size %d\n",correlationTablePixels);
  }
  correlation_table2 = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_table2 == NULL) {
    printf("ERROR: failed to allocate correlation table 2 of size %d\n",correlationTablePixels);
  }
  correlation_table3 = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_table3 == NULL) {
    printf("ERROR: failed to allocate correlation table 2 of size %d\n",correlationTablePixels);
  }
  correlation_table4 = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_table4 == NULL) {
    printf("ERROR: failed to allocate correlation table 4 of size %d\n",correlationTablePixels);
  }

  correlation_scratch1A = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_scratch1A == NULL) {
    printf("ERROR: failed to allocate correlation table scratch of size %d\n",correlationTablePixels);
  }
  correlation_scratch1B = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_scratch1B == NULL) {
    printf("ERROR: failed to allocate correlation table scratch of size %d\n",correlationTablePixels);
  }
  correlation_scratch2A = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_scratch2A == NULL) {
    printf("ERROR: failed to allocate correlation table scratch of size %d\n",correlationTablePixels);
  }
  correlation_scratch2B = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_scratch2B == NULL) {
    printf("ERROR: failed to allocate correlation table scratch of size %d\n",correlationTablePixels);
  }
  /* BEGINNING OF SECTOR LOOP (LOOP 1)*/
  for (maskSector = 0; maskSector < totalMaskSectors; maskSector++) {
#if DEBUG_TRUNCATION_X
    if (maskSector == 0) {
      maskSector = 14;
    } else {
      printf("ERROR: DEBUG_TRUNCATION abort in line %d\n",__LINE__);
      exit(-1);
    }
#endif /* DEBUG_TRUNCATION */
    pMaskSector = &maskSectorTable[maskSector];
    if (verbose) {
      printf("Processing sector %2d for %s\n",maskSector,fileroot);
    }

    pMaskSector->convolutionArea = 1;
    candidateCount = 0;
    searchTableSize = 0;
    maxIndex1 = 0;
    maxIndex2 = 0;
    maxExcess = -1000000.0;
    numPeaks = 0;
    sigmaSum = 0;
    sigmaSqr = 0;
    memset(searchTable,0,sizeof(searchTable));
    memset(correlation_histogram1,0,sizeof(correlation_histogram1));
    correlation_histogram1[0] = correlationTablePixels;
    memset(correlation_histogram2,0,sizeof(correlation_histogram2));
    correlation_histogram2[0] = correlationTablePixels;
    memset(correlation_table1,0,correlationTablePixels*sizeof(int));    
    memset(correlation_table2,0,correlationTablePixels*sizeof(int));
    memset(correlation_table3,0,correlationTablePixels*sizeof(int));
    memset(MAG_ISO_count1,0,MAG_ISO_size*MAG_ISO_size*sizeof(int));
    memset(MAG_ISO_count2,0,MAG_ISO_size*sizeof(int));
    memset(MAG_ISO_countmatch1,0,MAG_ISO_size*sizeof(int));
    memset(MAG_ISO_countmatch2,0,MAG_ISO_size*sizeof(int));

    table1Points = 0;
    table1Total = 0;
    table1Peak = 0;
    table1MaxPeak = 0;
    table1OldTotal = 0;
    table1MaxTotal = 0;
    table1MaxPoints = 0;
    table1MaxRatio = 0;
    table1_dimlimit = 0;
    table1Ratio = 0;
    table1MaxRatio = 0;
    
    vectorIndex = 0;
    correlation_sum0 = 0;
    correlation_sum1 = 0;
    pmf_sum0 = 0.0;
    bin_sum0 = 0.0;
    for (searchIndex = 0; searchIndex < SEARCH_TABLE_LENGTH; searchIndex++) {
      pSearchEntry1 = &searchTable[searchIndex];
      pSearchEntry1->sextractorIndex = -1;
    }
    ixs = maskSector % nxs;
    iys = maskSector / nxs;
    minX = (mosaicWidth * ixs)/nxs;
    minY = (mosaicHeight * iys)/nys;
    maxX = (mosaicWidth *(ixs+1))/nxs;
    maxY = (mosaicHeight * (iys+1))/nys;
    pMaskSector->aveX = (minX+maxX)/2.0;
    pMaskSector->aveY = (minY+maxY)/2.0;
    if (verbose) {
      printf("Sector %d Min, Max X: %d %d, Min Max Y %d %d\n",maskSector,minX,maxX,minY,maxY);
    }

 
    /* We assume that most of the objects are dim.  Calculate the clipped median to
       get the limiting dim MAG_ISO */
    if (vector != NULL) {
      free(vector);
      vector = NULL;
    }
    vectorSize = correlationTablePixels;
    if (sextractor_nrecs > vectorSize) {
      vectorSize = sextractor_nrecs;
    }
    vector = (double *)calloc(vectorSize,sizeof(double));
    /* At this point, we will produce a sorted list of the brightest objects on the plate
       with a reasonable FWHM */

    for (sextractorIndex = pMaskSector->sextractorIndexLow; sextractorIndex <= pMaskSector->sextractorIndexHigh; sextractorIndex++) {
      pSextractorBright = &image_table[sextractorIndex];
      if (pSextractorBright->maskSector != pMaskSector->maskSector) {
        printf("ERROR: line %d maskSector %d and %d do not agree\n",__LINE__,pSextractorBright->maskSector,pMaskSector->maskSector);
        exit(-1);
      }

      /* Reject saturated and bad images */
      if ((pSextractorBright->BFLAGS & BFLAGSMASK) != 0) {
#if 0
        printf("Rejecting BFLAGS %d\n",pSextractorBright->BFLAGS);
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
      MAG_ISO_brightindex = (pSextractorBright->MAG_ISO-MAG_ISO_brightest)/MAG_ISO_DELTA;
      if (MAG_ISO_brightindex < 0) {
        MAG_ISO_brightindex = 0;
      }
      if (MAG_ISO_brightindex >= MAG_ISO_size) {
        MAG_ISO_brightindex = MAG_ISO_size-1;
      }
      MAG_ISO_count2[MAG_ISO_brightindex]++;
      if (doPlots2) {
        if (gmtHandle[7] == NULL) {
          gmtHandle[7] = fopen(gmtname[7],"wt");
          if (gmtHandle[7] == NULL) {
            printf("ERROR: Failed to open the gmt file %s\n",gmtname[7]);
            exit(-1);
          } else {
            if (verbose) {
              printf("GMT file %s\n",gmtname[7]);
            }
          }
          fprintf(gmtHandle[7],"maskSector\tMAG_ISO\tFLUX_MAX\tanalysis_threshold\n");
          fprintf(gmtHandle[7],"----------\t-------\t--------\t-----------------\n");
        }
        if (pSextractorBright->THRESHOLD == 0) {
          analysis_threshold2 = 0;
        } else {
          analysis_threshold2 = pSextractorBright->FLUX_MAX/pSextractorBright->THRESHOLD;
        }
        fprintf(gmtHandle[7],"%d\t%f\t%f\t%f\n",
                maskSector,
                pSextractorBright->MAG_ISO,
                pSextractorBright->FLUX_MAX,
                analysis_threshold2);

      }
      if ((pSextractorBright->FLUX_MAX/pSextractorBright->THRESHOLD) < analysis_threshold) {
        continue;
      }

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
                memcpy(&searchTable[insertIndex],&searchTable[insertIndex-1],sizeof(SEARCHENTRY));
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
    if (verbose) {
      printf("Sector %d searchTableSize %d\n",maskSector,searchTableSize);
    }
#if 0  /* The table is already sorted according to MAG_ISO */     
    qsort(searchTable,searchTableSize,sizeof(SEARCHENTRY),searchCompare);
#endif
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
    
    /* BEGINNING OF THE BRIGHTNESS LOOP (LOOP 2) */
    for (MAG_ISO_index = 0; MAG_ISO_index < MAG_ISO_size; MAG_ISO_index++) {
      MAG_ISO_brightlimit = (MAG_ISO_index*MAG_ISO_DELTA) + MAG_ISO_brightest;
      MAG_ISO_dimlimit = ((MAG_ISO_index+1)*MAG_ISO_DELTA) + MAG_ISO_brightest;
      /* BEGINNING OF THE NEIGBOR SECTOR LOOP (LOOP 3) */
      for (neighborIndex = 0; neighborIndex < pMaskSector->neighborCount; neighborIndex++) {
        pMaskSector2 = &maskSectorTable[pMaskSector->maskNeighbor[neighborIndex]];
      

        /* Now go through the sextractor list for each of these objects and save likely candidates */
#if 0
        printf("line %d processing neighborIndex %2d maskSector %2d\n",__LINE__,neighborIndex,pMaskSector2->maskSector);
#endif
        /* BEGINNING OF THE MATCH DIM STAR LOOP (LOOP 4) */
        for (sextractorIndex =  pMaskSector2->sextractorIndexLow; sextractorIndex <= pMaskSector2->sextractorIndexHigh; sextractorIndex++) {
          pSextractorDim = &image_table[sextractorIndex];
          if (pSextractorDim->MAG_ISO < MAG_ISO_brightlimit) {
            continue; /* Too bright */
          }
          if (pSextractorDim->MAG_ISO > MAG_ISO_dimlimit) {
            break; /* Too dim */
          }
          if (pSextractorDim->maskSector != pMaskSector2->maskSector) {
            printf("ERROR: line %d maskSector %d and %d do not agree\n",__LINE__,pSextractorDim->maskSector,pMaskSector2->maskSector);
            exit(-1);
          }
          /* BEGINNING OF THE MATCH BRIGHT STAR LOOP (LOOP 5) */
          for (searchIndex = 0; searchIndex < searchTableSize; searchIndex++) {
            pSearchEntry1 = &searchTable[searchIndex];
            if (pSearchEntry1-> MAG_ISO >  MAG_ISO_dimlimit) {
              continue; /* Both stars are too dim */
            }
            NUMBER = image_table[pSearchEntry1->sextractorIndex].NUMBER;
            minX = pSearchEntry1->X_IMAGE  - MAX_PIXEL_RADIUS;
            maxX = pSearchEntry1->X_IMAGE  + MAX_PIXEL_RADIUS;
            minY = pSearchEntry1->Y_IMAGE  - MAX_PIXEL_RADIUS;
            maxY = pSearchEntry1->Y_IMAGE  + MAX_PIXEL_RADIUS;            
            if (pSextractorDim->NUMBER == NUMBER) {
              pSextractor2 = pSextractorDim;
              continue;
            }
#if 0
            if (neighborIndex > 0) {
              printf("line %d processing searchIndex %6d searchY_IMAGE %d sextractorIndex %6d maskSector %2d Y_IMAGE %f\n",__LINE__,searchIndex,pSearchEntry1->Y_IMAGE,sextractorIndex,pSextractorDim->maskSector,pSextractorDim->Y_IMAGE);
            }
#endif
            if (pSextractorDim->MAG_ISO < pSearchEntry1->MAG_ISO) {
              /* The search entry must be brighter */
              continue;
            }
#if 0 /* Excluded Oct 17, 2017 */
            if ((pSextractorDim->FLUX_MAX/pSextractorDim->THRESHOLD) < analysis_threshold) {
              continue;
            }
#endif


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
              printf("ERROR: correlationIndex %d exceeds limits (5) %d\n",correlationIndex,maxCorrelationIndex);
              exit(-1);
            }
#if 0
            if ((xIndex == -1) && (yIndex == -136)) {
              printf("At x %d y %d NUMBER %d %d\n",
                     xIndex,yIndex,pSextractorDim->NUMBER,pSearchEntry1->NUMBER);
            }
      
#endif
            
            if (correlation_table1[correlationIndex] == 0) {
              table1Points++;
            }
            correlation_histogram1[correlation_table1[correlationIndex]]--;
            correlation_table1[correlationIndex]++;
            correlation_histogram1[correlation_table1[correlationIndex]]++;
            
            
            table1Total++;
            if (table1Peak < correlation_table1[correlationIndex]) {
              table1Peak = correlation_table1[correlationIndex];
              if (table1Peak >= SEARCH_TABLE_LENGTH) {
                printf("ERROR: correlation_histogram1 index exceeded\n");
                exit(-1);
              }
            }

          } /* END OF SEARCH STAR LOOP (LOOP 5) */
        } /* END of DIM STAR LOOP (LOOP 4) */
      } /* END OF NEIGHBOR SECTOR SEARCH (LOOP 3) */


#ifdef COUNT_LOOP2
      if (table1Total > table1OldTotal) {
        if ((table1Points == 0) || (table1Peak < 2)) {
          table1Ratio = -1;
        } else {
#if 0
          if (table1MaskCount == 1) {
            printf("line %5d table1MaskCount %d\n",__LINE__,table1MaskCount);
          }
#endif
          table1MaskCount = CountMasks(table1Peak, correlationTablePixels,correlation_histogram1,correlation_table1,correlation_scratch1A,correlation_scratch1B,&table1NoiseFloor,fileroot);
          if (table1NoiseFloor == 0) {
            table1NoiseFloor = 1;
          }
          
          table1Ratio = (1.0*table1Peak)/(1.0*table1NoiseFloor);

        }      
        if (table1MaxRatio < table1Ratio) {
          table1MaxRatio = table1Ratio;
          table1_dimlimit = MAG_ISO_dimlimit;
#if 0
          memcpy(correlation_table4,correlation_table1,correlationTablePixels*sizeof(int));
#endif
          table1MaxTotal = table1Total;
          table1MaxPoints = table1Points;
          table1MaxPeak = table1Peak;
        }
#if 1
        printf("maskSector %2d MAG_ISO_dimlimit %5.1f table1Total %9d table1Points %8d table1Peak %5d table1NoiseFloor %d table1MaskCount %5d table1NoiseFloor %5d table1ratio %10.6f for %s\n",maskSector,MAG_ISO_dimlimit,table1Total,table1Points,table1Peak,table1NoiseFloor,table1MaskCount,table1NoiseFloor,table1Ratio,fileroot);
#endif
        table1OldTotal = table1Total;
      }
#endif /* COUNT_LOOP2 */
    } /* END OF MAG_ISO LOOP (LOOP 2) */
    printf("maskSector %2d table1Total %9d table1Points %8d table1Max %5d table1MaxRatio %10.6f table1_dimlimit %5.1f for %s\n",maskSector,table1MaxTotal,table1MaxPoints,table1MaxPeak,table1MaxRatio,table1_dimlimit,fileroot);
#if 0
    /* Copy back the best table */
    memcpy(correlation_table1,correlation_table4,correlationTablePixels*sizeof(int));
#endif

    for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
      for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
        correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
        if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
          printf("ERROR: correlationIndex %d exceeds limits (6) %d\n",correlationIndex,maxCorrelationIndex);
          exit(-1);
        }
        if (doPlots2) {
          if (gmtHandle[5] == NULL) {
            gmtHandle[5] = fopen(gmtname[5],"wt");
            if (gmtHandle[5] == NULL) {
              printf("ERROR: Failed to open the gmt file %s\n",gmtname[5]);
              exit(-1);
            } else {
              if (verbose) {
                printf("GMT file %s\n",gmtname[5]);
              }
            }
            fprintf(gmtHandle[5],"maskSector\txIndex\tyIndex\tcorrelationCount\n");
            fprintf(gmtHandle[5],"----------\t------\t------\t----------------\n");
          }
          if (correlation_table1[correlationIndex] > 0) {
            fprintf(gmtHandle[5],"%d\t%d\t%d\t%d\n",maskSector,xIndex,yIndex,correlation_table1[correlationIndex]);
          }
        }

        vector[vectorIndex] = correlation_table1[correlationIndex];
        vectorIndex++;
        if ((correlation_table1[correlationIndex]) > 0) {
          populatedBinCount++;
        }

      }
    }
    CalcMedianAndRMS(vectorIndex,10,vector,&pMaskSector->count_med0,&pMaskSector->count_rms0,1,3.0,0);



    if (pMaskSector->count_rms0 > 0) {
      pMaskSector->maxBinSNR0 = ((1.0*pMaskSector->maxBinCount0)-pMaskSector->count_med0)/pMaskSector->count_rms0;
    }
    if (searchTableSize == 0) {
      matchRatio2 = 0;
    } else {
      matchRatio2 = (1.0*pMaskSector->maxBinCount0)/(1.0*searchTableSize);
    }

    pMaskSector->count_limit0 = pMaskSector->count_med0 + (COUNT_LIMIT_RMS_FACTOR*pMaskSector->count_rms0);
    if (verbose) {
      printf("sector %2d initial count_med0 %f count_rms0 %f, count_limit0 %d maxBinSNR0 %f for %s\n",maskSector,pMaskSector->count_med0,pMaskSector->count_rms0,pMaskSector->count_limit0,pMaskSector->maxBinSNR0,fileroot);
    }
    if (doPlots1) {
      if (gmtHandle[1] == NULL) {
        gmtHandle[1] = fopen(gmtname[1],"wt");
        if (gmtHandle[1] == NULL) {
          printf("ERROR: Failed to open the gmt file %s\n",gmtname[1]);
          exit(-1);
        } else {
          if (verbose) {
            printf("GMT file %s\n",gmtname[1]);
          }
        }
        fprintf(gmtHandle[1],"maskSector\tX_IMAGE\tY_IMAGE\tconvolutionArea\tcount_med\tcount_rms\tcount_limit\tmaxBinCount\tmaxBinSNR\tselected\tmaskCount");
        fprintf(gmtHandle[1],"\tanalysis_threshold");
        fprintf(gmtHandle[1],"\tPlate\n");

        fprintf(gmtHandle[1],"----------\t--------\t------\t---------------\t---------\t---------\t-----------\t-----------\t---------\t--------\t---------");
        fprintf(gmtHandle[1],"\t------------------");
        fprintf(gmtHandle[1],"\t-----\n");


      }
#if 0
      fprintf(gmtHandle[1],"%d\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d",maskSector,pMaskSector->aveX,pMaskSector->aveY,pMaskSector->convolutionArea,pMaskSector->count_med0,pMaskSector->count_rms0,(1.0*pMaskSector->count_limit0),pMaskSector->maxBinCount0,pMaskSector->maxBinSNR0,selected,totalMasks);
      fprintf(gmtHandle[1],"\t%f",analysis_threshold);

      fprintf(gmtHandle[1],"\t%s\n",fileroot2);
#endif
    }


    table1MaskCount = CountMasks(table1Peak, correlationTablePixels,correlation_histogram1,correlation_table1,correlation_scratch1A,correlation_scratch1B,&table1NoiseFloor,fileroot);
    if (table1NoiseFloor == 0) {
      table1NoiseFloor = 1;
    }
          
    table1Ratio = (1.0*table1Peak)/(1.0*table1NoiseFloor);

    
    table1MaxRatio = table1Ratio;
    table1_dimlimit = MAG_ISO_dimlimit;
#if 0
    memcpy(correlation_table4,correlation_table1,correlationTablePixels*sizeof(int));
#endif
    table1MaxTotal = table1Total;
    table1MaxPoints = table1Points;
    table1MaxPeak = table1Peak;
#if 1
    printf("maskSector %2d MAG_ISO_dimlimit %5.1f table1Total %9d table1Points %8d table1Peak %5d table1MaskCount %5d table1NoiseFloor %5d table1ratio %10.6f for %s\n",maskSector,MAG_ISO_dimlimit,table1Total,table1Points,table1Peak,table1MaskCount,table1NoiseFloor,table1Ratio,fileroot);
#endif



    /* Now expand the convolution radius to the adjacent pixels */
    table2Points = 0;
    table2Total = 0;
    table2Peak = 0;
    table2Ratio = 0;


    memset(correlation_table2,0,(correlationTablePixels)*(sizeof(int)));
    memset(correlation_histogram2,0,sizeof(correlation_histogram2));
    correlation_histogram2[0] = correlationTablePixels;
    pMaskSector->convolutionArea = 0;
    for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
      for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
        correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
        if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
          printf("ERROR: correlationIndex %d exceeds limits (7) %d\n",correlationIndex,maxCorrelationIndex);
          exit(-1);
        }
        convolutionCount = 0;
        for (xIndex2 = (xIndex-DEFAULT_CONVOLUTION_RADIUS); xIndex2 <= (xIndex+DEFAULT_CONVOLUTION_RADIUS); xIndex2++) {
          if ((xIndex2 < -MAX_PIXEL_RADIUS) || (xIndex2 > MAX_PIXEL_RADIUS)) {
            continue;
          }
          for (yIndex2 = yIndex-DEFAULT_CONVOLUTION_RADIUS; yIndex2 <= yIndex+DEFAULT_CONVOLUTION_RADIUS; yIndex2++) {
            if ((yIndex2 < -MAX_PIXEL_RADIUS) || (yIndex2 > MAX_PIXEL_RADIUS)) {
              continue;
            }
            convolutionCount++;
            correlationIndex2 = (xIndex2+MAX_PIXEL_RADIUS) + ((yIndex2+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
            if ((correlationIndex2 < 0) || (correlationIndex2 >= maxCorrelationIndex)) {
              printf("ERROR: correlationIndex %d exceeds limits (8) %d\n",correlationIndex,maxCorrelationIndex);
              exit(-1);
            }
            if ((correlation_table2[correlationIndex2] == 0) &&  (correlation_table1[correlationIndex] != 0)) {
              table2Points++;
            }

            correlation_histogram2[correlation_table2[correlationIndex2]]--;
            correlation_table2[correlationIndex2] += correlation_table1[correlationIndex];
            correlation_histogram2[correlation_table2[correlationIndex2]]++;
            table2Total += correlation_table1[correlationIndex];
            if (table2Peak < correlation_table2[correlationIndex2]) {
              table2Peak = correlation_table2[correlationIndex2];
            }

            if (correlation_table2[correlationIndex2] > pMaskSector->maxBinCount1) {
              pMaskSector->maxBinCount1 = correlation_table2[correlationIndex2];
            }
          }
        }
        if (convolutionCount > pMaskSector->convolutionArea) {
          pMaskSector->convolutionArea = convolutionCount;
        }
      }
    }
    vectorIndex = 0;
    for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
      for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
        correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
        if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
          printf("ERROR: correlationIndex %d exceeds limits (9) %d\n",correlationIndex,maxCorrelationIndex);
          exit(-1);
        }
        vector[vectorIndex] = correlation_table2[correlationIndex];
        vectorIndex++;

      }
    }
    CalcMedianAndRMS(vectorIndex,10,vector,&pMaskSector->count_med1,&pMaskSector->count_rms1,1,3.0,0);
    if (pMaskSector->count_rms1 > 0) {
      pMaskSector->maxBinSNR1 = ((1.0*pMaskSector->maxBinCount1)-pMaskSector->count_med1)/pMaskSector->count_rms1;
    }
    table2MaskCount = CountMasks(table2Peak, correlationTablePixels,correlation_histogram2,correlation_table2,correlation_scratch2A,correlation_scratch2B,&table2NoiseFloor,fileroot);
    if (table2NoiseFloor == 0) {
      table2NoiseFloor = 2;
    }
          
    table2Ratio = (1.0*table2Peak)/(1.0*table2NoiseFloor);



    pMaskSector->count_limit1 = pMaskSector->count_med1 + (COUNT_LIMIT_RMS_FACTOR*pMaskSector->count_rms1);
    if (table2Points == 0) {
      table2Ratio = -1;
    } else {
      table2Ratio = (1.0*table2Peak)/((1.0*table2Total)/(1.0*table2Points));
    }      
    if (verbose) {
      printf("maskSector %2d table1Total %9d table1Points %8d table1Max %5d table1MaxRatio %10.6f table1_dimlimit %5.1f  table2Total %9d table2Points %8d table2Peak %5d table2Ratio %10.6f for %s\n",
             maskSector,
             table1MaxTotal,
             table1MaxPoints,
             table1MaxPeak,
             table1MaxRatio,
             table1_dimlimit,
             table2Total,
             table2Points,
             table2Peak,
             table2Ratio,fileroot);
    }
    if (doPlots1) {
      if (gmtHandle[9] == NULL) {
        gmtHandle[9] = fopen(gmtname[9],"wt");
        if (gmtHandle[9] == NULL) {
          printf("ERROR: Failed to open the gmt file %s\n",gmtname[9]);
          exit(-1);
        } else {
          if (verbose) {
            printf("GMT file %s\n",gmtname[9]);
          }
        }
        fprintf(gmtHandle[9],"table1MaskCount");
        fprintf(gmtHandle[9],"---------------");

      }

      fprintf(gmtHandle[9],"%d\t%f\t%f\t%d\t%f\t%d\t%f\t%f\t%s\n",maskSector,pMaskSector->aveX,pMaskSector->aveY,table1Peak,table1MaxRatio,table2Peak,table2Ratio,table1_dimlimit,fileroot2);
              
    }



    if (doPlots1) {
      if (gmtHandle[8] == NULL) {
        gmtHandle[8] = fopen(gmtname[8],"wt");
        if (gmtHandle[8] == NULL) {
          printf("ERROR: Failed to open the gmt file %s\n",gmtname[8]);
          exit(-1);
        } else {
          if (verbose) {
            printf("GMT file %s\n",gmtname[8]);
          }
        }
        fprintf(gmtHandle[8],"maskSector\tX_IMAGE\tY_IMAGE\ttable1Peak\ttable1MaxRatio\ttable2Peak\ttable2Ratio\ttable1_dimlimit\tPlate\n");
        fprintf(gmtHandle[8],"----------\t--------\t------\t----------\t--------------\t----------\t-----------\t---------------\t-----\n");

      }
      fprintf(gmtHandle[8],"%d\t%f\t%f\t%d\t%f\t%d\t%f\t%f\t%s\n",maskSector,pMaskSector->aveX,pMaskSector->aveY,table1Peak,table1MaxRatio,table2Peak,table2Ratio,table1_dimlimit,fileroot2);
              
    }

    if (verbose) {
      printf("sector %d convolutionArea %d count_med1 %f count_rms1 %f, count_limit1 %f maxBinCount1 %d maxBinSNR1 %f for %s\n",maskSector,pMaskSector->convolutionArea,pMaskSector->count_med1,pMaskSector->count_rms1,(1.0*pMaskSector->count_limit1),pMaskSector->maxBinCount1,pMaskSector->maxBinSNR1,fileroot);
    }


    if (doPlots1) {
      for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
        for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
          correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
          if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
            printf("ERROR: correlationIndex %d exceeds limits (10) %d\n",correlationIndex,maxCorrelationIndex);
            exit(-1);
          }
          /* Was ((correlation_table2[correlationIndex]) > pMaskSector->count_limit1) */
          if ((correlation_table1[correlationIndex]) > 0) {
            if (gmtHandle[0] == NULL) {
              gmtHandle[0] = fopen(gmtname[0],"wt");
              if (gmtHandle[0] == NULL) {
                printf("ERROR: Failed to open the gmt file %s\n",gmtname[0]);
                exit(-1);
              } else {
                if (verbose) {
                  printf("GMT file %s\n",gmtname[0]);
                }
              }
              fprintf(gmtHandle[0],"maskSector\tix\tiy\tbin_count\tmaskIndex\tsequence\n");
              fprintf(gmtHandle[0],"----------\t--\t--\t---------\t---------\t--------\n");
            }
#if 1
            fprintf(gmtHandle[0],"%d\t%d\t%d\t%d\t%d\t%d\n",maskSector,xIndex,yIndex,correlation_scratch2A[correlationIndex],correlation_scratch2B[correlationIndex],3);
#else
            fprintf(gmtHandle[0],"%d\t%d\t%d\t%d\t%d\t%d\n",maskSector,xIndex,yIndex,correlation_table1[correlationIndex],correlation_scratch1B[correlationIndex],0);
            fprintf(gmtHandle[0],"%d\t%d\t%d\t%d\t%d\t%d\n",maskSector,xIndex,yIndex,correlation_table2[correlationIndex],correlation_scratch1B[correlationIndex],1);
            fprintf(gmtHandle[0],"%d\t%d\t%d\t%d\t%d\t%d\n",maskSector,xIndex,yIndex,correlation_scratch1A[correlationIndex],correlation_scratch1B[correlationIndex],2);
#endif

          }

        }
      }
    }

    /* Now we need to identify the peak correlation masks and print a table of mask characteristics */
#ifdef DEBUG_TRUNCATION
    if (debugPrint == 0) {
      printf("ERROR: setting totalMasks to zero in line %d\n",__LINE__);
      debugPrint = 1;
#if 0
      exit(-1);
#endif
    }
    totalMasks = 0;
#else /* DEBUG_TRUNCATION */
    FindMasks(pMaskSector->count_limit1,correlation_table2,correlation_table3,&maskTable,&totalMasks,fileroot);
#endif /* DEBUG_TRUNCATION */

    if (totalMasks >= SEARCH_TABLE_LENGTH) {
      printf("ERROR search_close2 totalMasks %d greater then SEARCH_TABLE_LENGTH %d for %s\n",totalMasks,SEARCH_TABLE_LENGTH,fileroot);
      exit(-1);
    }

    if (doPlots1) {
      fprintf(gmtHandle[1],"%d\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d",maskSector,pMaskSector->aveX,pMaskSector->aveY,pMaskSector->convolutionArea,pMaskSector->count_med1,pMaskSector->count_rms1,(1.0*pMaskSector->count_limit1),pMaskSector->maxBinCount1,pMaskSector->maxBinSNR1,selected,totalMasks);
      fprintf(gmtHandle[1],"\t%f",analysis_threshold);
        
      fprintf(gmtHandle[1],"\t%s\n",fileroot2);
    }

    /* Expand all entries in the mask table by the correlation radius and find the limits.  Set the centerFlag if necessary */
    for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
      pMask = &maskTable[maskIndex];
      pMask->maskXMin -= DEFAULT_CONVOLUTION_RADIUS;
      pMask->maskXMax += DEFAULT_CONVOLUTION_RADIUS;
      pMask->maskYMin -= DEFAULT_CONVOLUTION_RADIUS;
      pMask->maskYMax += DEFAULT_CONVOLUTION_RADIUS;
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
       3.  Using the sorted sequence, increment the mask brightness count table entry [i,j] if mask i is brighter than mask j
       4.  After all host stars are processed find the order of masks based on the highest mask brightness counts
    */
    if (totalMasks > 1) {
      maskOrderTable = (PMASKORDER)calloc(totalMasks,sizeof(MASKORDER));
      if (maskOrderTable == NULL) {
        printf("ERROR: failed to allocate maskOrderTable of size %d\n",totalMasks);
        exit(-1);
      }
      memset(maskBrightnessTable,0,sizeof(maskBrightnessTable));
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
          if ((pSextractorDim->FLUX_MAX/pSextractorDim->THRESHOLD) < analysis_threshold) {
            continue;
          }
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
          qsort(maskOrderTable,totalMasks,sizeof(MASKORDER),maskOrderCompare);
     

#if 0
          for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
            pMaskOrder = &maskOrderTable[maskIndex];
            printf("%6.2f ",pMaskOrder->MAG_ISO);
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
                printf("ERROR: search_close2 maskIndex %d %d greater than %d for %s\n",pMaskOrder->maskIndex,pMaskOrder2->maskIndex,SEARCH_TABLE_LENGTH,fileroot);
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

        printf("%2d %3d:  ",pMask->maskIndex,maskIndex);
        for (maskIndex2 = 0; maskIndex2 < totalMasks; maskIndex2++) {
          printf("%4d ",maskBrightnessTable[maskIndex][maskIndex2]);
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

        printf("%2d %3d:  ",pMask->maskIndex,maskIndex);
        for (maskIndex2 = 0; maskIndex2 < totalMasks; maskIndex2++) {
          printf("%4d ",maskBrightnessTable[maskIndex][maskIndex2]);
        }
        printf("\n");    
      }
#endif
      qsort(maskTable,totalMasks,sizeof(MULTMASK),maskBrightnessCompare);
#if 0
      for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
        pMask = &maskTable[maskIndex];

        printf("%2d %3d:  ",pMask->maskIndex,maskIndex);
        for (maskIndex2 = 0; maskIndex2 < totalMasks; maskIndex2++) {
          printf("%4d ",maskBrightnessTable[maskIndex][maskIndex2]);
        }
        printf("\n");    
      }
#endif

    }
 
    /* Save the results in the database */
    if (totalMasks > 0) {
      SetMosaicFitWCS(pConnection,"HaveMultipleMask",fileroot,0,1,readFlag,&FitWCS);
      if ((totalMasks <= MAX_MASKS) &&
          (maskCenterIndex == 0)) {
        SetPlateQuality(pConnection,"multiple",fileroot,1,0,&quality);
      }


    } else {
      SetMosaicFitWCS(pConnection,"NoMultipleMask",fileroot,0,1,readFlag,&FitWCS);
    }


    for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
      pMask = &maskTable[maskIndex];
      pMask->maskIndex = maskIndex+1;
#ifndef DEBUG_TRUNCATION
      /* Load the masks into the datbase, but only 1 mask if there is an unusual condition */
      if (((maskCenterIndex == 0) && (totalMasks <= MAX_MASKS)) ||
          ((maskCenterIndex != 0) && (pMask->maskCenterFlag != 0)) ||
          ((totalMasks > MAX_MASKS) && (maskIndex == 0))) {
        sprintf(queryString,"INSERT INTO masks2 (series,plateNumber,mosaicNumber,nxs,nys,maskSector,maskCount,maskIndex,maskXMin,maskXMax,maskYMin,maskYMax,maskArea,maxBinCount,maskCenterFlag,maskCenterDistance,maskAreaRatio,convolutionRadius) values ('%s',%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%d);",
                series,
                plateNumber,
                mosaicNumber,
                nxs,
                nys,
                maskSector,
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
                DEFAULT_CONVOLUTION_RADIUS);
        if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
          printf("ERROR: search_close2 MAX_QUERY_STRING exceeded\n");
          exit(-1);
        }
#if 0
        printf("%s\n",queryString);
#endif

        res = ExecuteQuery(pConnection,queryString);
        if (res) {
          exit(-1);
        }
      }
#endif /* DEBUG_TRUNCATION */

      if (doPlots1) {
        if (gmtHandle[2] == NULL) {
          gmtHandle[2] = fopen(gmtname[2],"wt");
          if (gmtHandle[2] == NULL) {
            printf("ERROR: Failed to open the gmt file %s\n",gmtname[2]);
            exit(-1);
          } else {
            if (verbose) {
              printf("GMT file %s\n",gmtname[2]);
            }
          }

          fprintf(gmtHandle[2],"maskCount\tmaskIndex\tmaskXMin\tmaskXMax\tmaskYMin\tmaskYMax\tmaskArea\tmaxBinCount\tmaskCenterFlag\tmaskCenterDistance\tmaskAreaRatio\tPlate\n");
          fprintf(gmtHandle[2],"---------\t---------\t--------\t--------\t--------\t--------\t--------\t-----------\t--------------\t------------------\t-------------\t-----\n");
  
        }
        if (gmtHandle[6] == NULL) {
          gmtHandle[6] = fopen(gmtname[6],"wt");
          if (gmtHandle[6] == NULL) {
            printf("ERROR: Failed to open the gmt file %s\n",gmtname[6]);
            exit(-1);
          } else {
            if (verbose) {
              printf("GMT file %s\n",gmtname[6]);
            }
          }
          fprintf(gmtHandle[6],"xIndex\tyIndex\tmaskIndex\n");
          fprintf(gmtHandle[6],"------\t------\t---------\n");
  
        }
        fprintf(gmtHandle[2],"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.3f\t%s\n",
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
                fileroot);
        for (xIndex = pMask->maskXMin; xIndex <= pMask->maskXMax; xIndex++) {
          for (yIndex = pMask->maskYMin; yIndex <= pMask->maskYMax; yIndex++) {
            fprintf(gmtHandle[6],"%d\t%d\t\%d\n",xIndex,yIndex,pMask->maskIndex);
          }
        }
      }
    }

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
      for (sextractorIndex1 = 0; sextractorIndex1 < sextractor_nrecs; sextractorIndex1++) {
        pSextractor1 = &image_table[sextractorIndex1];
        if ((pSextractor1->FLUX_MAX/pSextractor1->THRESHOLD) < analysis_threshold) {
          continue;
        }
    



        for (sextractorIndex2 = 0; sextractorIndex2 < sextractor_nrecs; sextractorIndex2++) {
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
                   xIndex,yIndex,pSextractor1->NUMBER,pSextractor2->NUMBER);
          }
#endif


          if (pSextractor2->MAG_ISO < pSextractor1->MAG_ISO) {
            continue;
          }
          if (pSextractor2->NUMBER == pSextractor1->NUMBER) {
            continue;
          }
          if ((pSextractor2->FLUX_MAX/pSextractor2->THRESHOLD) < analysis_threshold) {
            continue;
          }
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
    if (maxBinCountOverall < pMaskSector->maxBinCount1) {
      maxBinCountOverall = pMaskSector->maxBinCount1;
    }

  } /* END OF SECTOR SEARCH LOOP (maskSector) (LOOP 1) */
  pMaskSector = NULL;
  if (outputCount > 0) {
    for (sextractorIndex1 = 0; sextractorIndex1 < sextractor_nrecs; sextractorIndex1++) {
      pSextractor1 = &image_table[sextractorIndex1];
      if (pSextractor1->maskIndex == 0) {
        continue;
      }
      if ((pSextractor1->FLUX_MAX/pSextractor1->THRESHOLD) < analysis_threshold) {
        continue;
      }    



      if (outHandle == NULL) {
        outHandle = fopen(outfile,"wt");
        if (outHandle == NULL) {
          fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
          exit(-1);
        } else {
          if (verbose) {
            fprintf(stderr,"Output file %s\n",outfile);
          }
          fprintf(outHandle,"NUMBER\tmaskIndex\tX_IMAGE\tY_IMAGE\n");
          fprintf(outHandle,"------\t---------\t-------\t-------\n");
          
        }
      }
      fprintf(outHandle,"%d\t%d\t%d\t%d\n",
              pSextractor1->NUMBER,
              pSextractor1->maskIndex,
              pSextractor1->X_IMAGE,
              pSextractor1->Y_IMAGE);


      if (doPlots1) {
        if (gmtHandle[4] == NULL) {
          gmtHandle[4] = fopen(gmtname[4],"wt");
          if (gmtHandle[4] == NULL) {
            printf("ERROR: Failed to open the gmt file %s\n",gmtname[4]);
            exit(-1);
          } else {
            if (verbose) {
              printf("GMT file %s\n",gmtname[4]);
            }
          }
        }
        fprintf(gmtHandle[4],"CIRCLE(%d,%d,5) # text = {%d}\n",pSextractor1->X_IMAGE/binning,pSextractor1->Y_IMAGE/binning,pSextractor1->maskIndex);    
      }
      if ((pSextractor1->convolutionFlag == 0) || (pSextractor1->pSextractor1 == NULL)) {
        continue;
      }
      if (doPlots1) {
        if (gmtHandle[3] == NULL) {
          gmtHandle[3] = fopen(gmtname[3],"wt");
          if (gmtHandle[3] == NULL) {
            printf("ERROR: Failed to open the gmt file %s\n",gmtname[3]);
            exit(-1);
          } else {
            if (verbose) {
              printf("GMT file %s\n",gmtname[3]);
            }
          }
        }
        fprintf(gmtHandle[3],"LINE(%d,%d,%d,%d) # text = {%d} line = 0 1 \n",pSextractor1->pSextractor1->X_IMAGE/binning,pSextractor1->pSextractor1->Y_IMAGE/binning,pSextractor1->X_IMAGE/binning,pSextractor1->Y_IMAGE/binning,pSextractor1->maskIndex);    
      }

    }

  }
  /* Update the quality bits for this plate */
  UpdateQuality(pConnection,series,plateNumber);

  /* All done.  Clean up */
 

  if (outHandle != NULL) {
    fclose(outHandle);
  }
  if (doPlots1) {
    for (gmtIndex = 0; gmtIndex < NUM_GMT_FILES; gmtIndex++) {
      if (gmtHandle[gmtIndex] != NULL) {
        fclose(gmtHandle[gmtIndex]);
      }
    }
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
  if (correlation_table3 != NULL) {
    free(correlation_table3);
  }
  if (correlation_table4 != NULL) {
    free(correlation_table4);
  }
  if (correlation_scratch1A != NULL) {
    free(correlation_scratch1A);
  }
  if (correlation_scratch1B != NULL) {
    free(correlation_scratch1B);
  }
  if (correlation_scratch2A != NULL) {
    free(correlation_scratch2A);
  }
  if (correlation_scratch2B != NULL) {
    free(correlation_scratch2B);
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
  if (maskSectorTable != NULL) {
    free(maskSectorTable);
  }
  if (MAG_ISO_count1 != NULL) {
    free(MAG_ISO_count1);
  }
  if (MAG_ISO_count2 != NULL) {
    free(MAG_ISO_count2);
  }
  if (pBinMap != NULL) {
    free(pBinMap);
  }
  if (MAG_ISO_countmatch1 != NULL) {
    free(MAG_ISO_countmatch1);
  }
  if (MAG_ISO_countmatch2 != NULL) {
    free(MAG_ISO_countmatch2);
  }


  mysql_close(pConnection);

  time(&curTime);
  curTime -= startTime;

  if (verbose) {
    printf("Candidates: %d plotted %d psfsaturated %d populatedBinCount %d maxBinCountOverall %d outputCount %d for %s\n",
           candidateCount,
           plotCount,
           psfsaturatedCount,
           populatedBinCount,
           maxBinCountOverall,
           outputCount,
           fileroot);

  }
  printf(" threshold %f",analysis_threshold);
  printf(" totalMasks %d seconds %d for %s\n",totalMasks,curTime,fileroot);


  return(0);
}

