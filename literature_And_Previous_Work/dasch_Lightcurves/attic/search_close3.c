// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* search_close3.c
 *  
 *  Search for multiple images near the brightest stars.  This is a revision of the search_close program
 *  On Mon, 8 May 2017 16:09:28 approximately 17% of 494811 plates are known to be multiple, including grating and wedge plates
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm -L/usr/lib64/mysql -L/usr/lib/mysql -lmysqlclient search_close3.c pipelineutils.a -ltable -lutil -lwcs -o search_close3 
 *`
 *  search_close3 -v -s 1.786 -w 17412 -h 22026 -r i31090_00_01r180ww -i /home/scanner/Pipeline/match/i31090_00_01r180ww_tnx.tmp -o /home/scanner/Pipeline/match/i31090_00_01r180ww_close.db
 *
 * <plate>_close_gmt0.txt - 3d plot of (x,y,bin_count)
 * <plate>_close_gmt1.txt - convolution statistics
 * <plate>_close_gmt2.txt - mask table
 * <plate>_close_gmt3.txt - region file for masked points used to create the convolution table
 * <plate>_close_gmt4.txt - region file for all masked points
 * <plate>_close_gmt5.txt - correlationCount for unbinned results
 * <plate>_close_gmt6.txt - dump of mask indices for debug
 * <plate>_close_gmt7.txt - Histogram for the 9x9 binned results
 * <plate>_close_plot.db  - maps of sextractor images
 *
 *
 * background map
 *
   column -i i31090_00_01r180ww_tnx.db NUMBER X_IMAGE Y_IMAGE | sorttable -n NUMBER > los1.db
   jointable -n -j NUMBER los1.db i31090_00_01r180ww_background.db > los2.db
   votable -i los2.db -o los.xml

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
 *
 * Mar 30, 2020 Edward J. Los   create search_close3 from oldest available version of search_close2 (Dec 28 2017 or Jan 3, 2018)
 * Apr  4, 2020 Edward J. Los   Add patternID, dimSeletedCount, searchTableSize to parameter values checked by parselogs.c
 *                              Add SEARCHTABLE_ONLY to restrict searches to the searchTable
 * Apr  6, 2020 Edward J. Los   Apply REJECT_GRAIN to apply the filter_wedge limit:  (pSextractorBright->MAG_ISO < (MAG_ISO_med - MAG_ISO_rms) to eliminate grain noise
 *                              Apply COUNT_LIMIT1 8 when building masks
 *                              Reject any zone where SEARCH_TABLE_LENGTH can not be obtained
 *                              Reduce SEARCH_TABLE_LENGTH to 1990 (to include br00777, which has the lowest good point count in the overlap zone)
 * Apr 10, 2020 Edward J. Los   Investigate using the area of each mask
 * Apr 12, 2020 Edward J. Los   Flag QUALITY_MULTIPLE condition from the old SEARCH_CLOSE3 algorithm.
 *                              Add "getfits" and rotation commands for visual examination of images.  Zone size  INSPECTION_WIDTH  x INSPECTION_HEIGHT
 *                              Add -g to specify the getfits target directory
 * Apr 14, 2020 Edward J. Los   If there are no suitable objects for getfits in the dual zones, use the single zone.
 *                              For getfits, do not include the bottom and top edges
 * Apr 16, 2020 Edward J. Los   Do not force searchTableSize to zero if there are insufficient stars in any maskSector
 * Oct 31, 2020 Edward J. Los   Add DEBUG_VERT4X5 for temperary buffix. 
 * Nov  3, 2020 Edward J. Los   Replace DEBUG_VERT4X5 with wideScanLocationPatternID 
 *                               for the widest scanLocation range representing the smallest mosaic range.  
 *                               NOTE: this fix has not been tested for all combinations.
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include <sys/stat.h>


#define SECTOR_PIXEL_WIDTH 2116; /* Minimum width of a sector = 23.27 mm or 0.91 inch */
#define MAX_SECTOR_COUNT 252 /* Current maximum sector count */
#define DEFAULT_CONVOLUTION_RADIUS 1
#define MAX_BUFFER 100
#define SEARCH_TABLE_LENGTH 1990 /* Changed from 2500 on Apr 7, 2020 */
#define NUM_GMT_FILES 8
#define MATCH_ARRAY_ALLOC 1000000
#define MAG_ISO_RMS_FACTOR 4
#define SQUARE_BIN_OBJECTS 10
#define CLIP_RMS_FACTOR 1.0
#define EXPECTED_PEAKS 1
#define EXCESS_LIMIT 10.0
#define MAX_PIXEL_RADIUS 1000
/* #define LOS_DEBUG 3   This is the desired bin number - use 3 for mc05077 */
/* #define DEBUG_TRUNCATION 1 *//* Halt processing early (ave time for 72 mosaics: 1305 seconds median 1134 seconds 5.67 hr total */
#define DEBUG_NO_MYSQL 1  /* If defined, make no changes to the MySQL databases */
#define DEBUG_DELETE_MASKS 1 /* Delete the masks after every mask sector */
#define SEARCHTABLE_ONLY /* Select dim candidates only from the search table */
#define REJECT_GRAIN  1   /* Reject grain noise */
#define COUNT_LIMIT1  8 /* Force FindMasks to use this limit */
#define SKIP_POISSON 1
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
#define INSPECTION_WIDTH 1920 /* Pixels along the original X axis */
#define INSPECTION_HEIGHT 764 /* Pixels along the original Y axis */

/* #define BACKGROUND_PLOTS 1 */

/* Threshold = 2.0, use 8-15 */
#define COUNT_LIMIT_RMS_FACTOR 15

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
  double X_IMAGE;
  double Y_IMAGE;
  int maskSector;
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
  int sextractorIndex; /* Index after the Y_IMAGE sort */
  int dimSelectedFlag; /* Image considered as the dim star of the pair if nonzero */
  int convolutionFlag; /* Image contributed to the convolution table if nozero.   */
  int searchStar;      /* Image is in the search table if nonzero */
  int maskSector;
  struct _starimage *pSextractor1;
} STARIMAGE,*PSTARIMAGE;
#define MASKSECTOR_FLAG 0x544345534b53414d /* MASKSECT */


#define MAX_SECTOR_NEIGHBOR 9 /* Includes this current sector */
typedef struct _masksector {
  long long sectorFlag;
  int maskSector;
  int neighborCount;
  int maskNeighbor[MAX_SECTOR_NEIGHBOR];
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
  double aveXY;
  int sectorCount;
  int unsaturatedCount;
  int dimSelectedCount;
  double MAG_ISO_med;
  double MAG_ISO_rms;
  int curClipCount;
  int maxMaskArea;
  int maxCenterMaskArea;
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

typedef struct _scanlocation {
  int patternID; /* scan pattern ID */
  int tile;      /* tile number */
  double xPos;      /* Aerotech X position for the tile */
  double yPos;      /* Aerotech Y position for the tile */
} SCANLOCATION, *PSCANLOCATION;
#define MAX_TILES 160 /* This should handle a 14" x 17" plate */
#define AXIS_FLAG_X 0 /* This zone maps to X_IMAGE */
#define AXIS_FLAG_Y 1 /* This zone maps to Y_IMAGE */
#define MAX_MIRROR_ERROR  2 /* If we invert the X or Y axis, we should get the same results to within this number of X or Y pixels.  Otherwise, we need to consider the flip, mirror, and rotate transformations */
typedef struct _overlapentry{
  int tileOverlapCount; /* Number of tiles overlapping withing this zone */
  int imageCount;       /* SExtractor images assigned to this zone */
  double startPos;
  double endPos;
  double startPixel;
  double endPixel;
  int startSextractorIndex;
  int endSextractorIndex;
  int* inspectionTileTableLower;  /* images on lower boundary of the entry */
  int* inspectionTileTableUpper;  /* images on upper boundary of the entry */
} OVERLAPENTRY,*POVERLAPENTRY;
  /* 
     Transforms:
     1.  flipFlag = 1    invertYArrayFlag = invertYArrayFlag % 1   
     2.  mirrorFlag = 1  invertXArrayFlag = invertXArrayFlag % 1
     3.  rotation     0
                     90  swap invertXArrayFlag and  invertXArrayFlag then invertXArrayFlag = invertXArrayFlag % 1
                    180  invertYArrayFlag = invertYArrayFlag % 1  and invertXArrayFlag = invertXArrayFlag % 1
                    270  swap invertXArrayFlag and  invertXArrayFlag then invertYArrayFlag = invertYArrayFlag % 1 

                                                     
     if xPosXImageFlag == 1, then  xImageArrayLength = naxis1;  yImageArrayLength = naxis2;
                                   X_IMAGE = (((xPos - xPosAve) / NOMINAL_MM_PER_PIXEL) + (naxis1/2))  or xPos = xPosAve + (NOMINAL_MM_PER_PIXEL * (X_IMAGE - (xImageArrayLength/2))); 
     and                           Y_IMAGE = (((yPos - yPosAve) / NOMINAL_MM_PER_PIXEL) + (naxis2/2))  or yPos = yPosAve + (NOMINAL_MM_PER_PIXEL * (Y_IMAGE - (yImageArrayLength/2))); 

and if  xPosXImageFlag == 0, then  xImageArrayLength = naxis1;  yImageArrayLength = naxis2;
                                   X_IMAGE = (((yPos - yPosAve) / NOMINAL_MM_PER_PIXEL) + (naxis1/2))  or yPos = yPosAve + (NOMINAL_MM_PER_PIXEL * (X_IMAGE - (xImageArrayLength/2) )); 
and                                Y_IMAGE = (((xPos - xPosAve) / NOMINAL_MM_PER_PIXEL) + (naxis1/2))  or xPos = xPosAve + (NOMINAL_MM_PER_PIXEL * (Y_IMAGE - (yImageArrayLength/2) )); 


  */
typedef struct _overlapzone {
  int origAxisFlag;  
  int curAxisFlag;
  int nZones; /* number of zones */
  int patternID;
  int numTiles;
  int flipFlag;
  int mirrorFlag;
  int rotation;
  double xPosAve;
  double yPosAve;
  int xPosXImageFlag;  /* If set, the X_IMAGE axis of SExtractor matches the xPos entries in the scan patterns.  Otherwise the Y_IMAGE axis matches */
  int maxInspectionTiles; /* number of inspection areas for getfits */
  int xImageArrayLength; 
  int invertXArrayFlag;
  int *pXImageArray;
  int yImageArrayLength;
  int invertYArrayFlag;
  int *pYImageArray;
  int unassignedImageCount; /* SExtractor images that did not fit any zone */
  int xZoneCount;
  int yZoneCount;
  OVERLAPENTRY XOverlapEntryTable[MAX_TILES];
  OVERLAPENTRY YOverlapEntryTable[MAX_TILES];
  PSCANLOCATION pScanLocationTable;
  
} OVERLAPZONE,*POVERLAPZONE;


int wideScanLocationPatternID[] = {50,51}; /* Bugfix of Nov 3, 2020 */
int numWideScanLocationPatternID = (sizeof(wideScanLocationPatternID)/sizeof(int));
typedef struct _xposximage {
  int rangeCondition1; /* patternID in wideScanLocationPatternID */
  int rangeCondition2; /* (xRange > yRange) */
  int rangeCondition3; /* (pMosaic->naxis1 > pMosaic->naxis2) */
  int xPosXImageFlag;
} XPOSXIMAGE,*PXPOSXIMAGE;

XPOSXIMAGE xPosXImageTable[] = {
  {0,1,1,1},
  {0,1,0,0},
  {0,0,1,0},
  {0,0,0,1},
  {1,1,1,0},
  {1,1,0,1},
  {1,0,1,1},
  {1,0,0,0},
};

int xPosXImageTableSize = sizeof(xPosXImageTable,sizeof(XPOSXIMAGE));

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
double fmin(double x, double y);
double fmax(double x, double y);
int yCompare(const void *first, const void *second) 
{
  double yFirst = ((PINPUTIMAGE)first)->Y_IMAGE;
  double ySecond = ((PINPUTIMAGE)second)->Y_IMAGE;
  int maskSectorFirst = ((PINPUTIMAGE)first)->maskSector;
  int maskSectorSecond = ((PINPUTIMAGE)second)->maskSector;


  if (maskSectorFirst > maskSectorSecond) {
    return(1);
  } else if (maskSectorFirst < maskSectorSecond) {
    return(-1);
  } else {
    if (yFirst > ySecond) {
      return(1);
    } else if (yFirst < ySecond) {
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
  double yFirst = ((PSEARCHENTRY)first)->Y_IMAGE;
  double ySecond = ((PSEARCHENTRY)second)->Y_IMAGE;
  int maskSectorFirst = ((PSEARCHENTRY)first)->maskSector;
  int maskSectorSecond = ((PSEARCHENTRY)second)->maskSector;


  if (maskSectorFirst > maskSectorSecond) {
    return(1);
  } else if (maskSectorFirst < maskSectorSecond) {
    return(-1);
  } else {
    if (yFirst > ySecond) {
      return(1);
    } else if (yFirst < ySecond) {
      return(-1);
    } else {
      return(0);
    }
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


void FindMasks(int maskSector,int tileOverlapCount,int count_limit,int *correlation_table2,int *correlation_table3,PMULTMASK *pMaskTable,int* pTotalMasks,char *fileroot) 
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
    printf("line %4d search_close3 Sector %2d overlap %2d Mask %4d xIndex %4d %4d yIndex %4d %4d maskArea %4d ratio %0.3f maxBinCount %4d maskCenterFlag %d maskCenterDistance %6.1f for %s\n",
           __LINE__,
           maskSector,
           tileOverlapCount,
           pMask->maskIndex,
           pMask->maskXMin,
           pMask->maskXMax,
           pMask->maskYMin,
           pMask->maskYMax,
           pMask->maskArea,
           pMask->maskAreaRatio,
           pMask->maxBinCount,
           pMask->maskCenterFlag,
           pMask->maskCenterDistance,
           fileroot);
          
#endif
  }
#if 0
  printf("line %4d search_close3 Sector %2d overlap %2d FinalMaskCount is %4d for %s\n",__LINE__,maskSector,tileOverlapCount,finalMaskCount,fileroot);
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
void FreeOverlapZone(POVERLAPZONE *ppOverlapZone) {
  if (*ppOverlapZone == NULL) {
    return;
  }
  POVERLAPZONE pOverlapZone = *ppOverlapZone;
  if (pOverlapZone->pScanLocationTable != NULL) {
    free(pOverlapZone->pScanLocationTable);
    pOverlapZone->pScanLocationTable = NULL;
  }

  if (pOverlapZone->pXImageArray != NULL) {
    free(pOverlapZone->pXImageArray);
    pOverlapZone->pXImageArray = NULL;
  }
  if (pOverlapZone->pYImageArray != NULL) {
    free(pOverlapZone->pYImageArray);
    pOverlapZone->pYImageArray = NULL;
  }

  free(pOverlapZone);
  *ppOverlapZone = NULL;
  return;
}
int GetOverlapZones(void *connection,char *series,int plateNumber,int mosaicNumber,int tempbinning,int rotation,PMOSAIC pMosaic,POVERLAPZONE *ppOverlapZone) {
  char queryString[MAX_QUERY_STRING]; 
  MYSQL *pConnection  = (MYSQL *)connection;
  int res = 0;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int nvals;
  int curTileCount = 0;
  int curTileIndex;
  POVERLAPZONE pOverlapZone = NULL;
  PSCANLOCATION pScanLocation = NULL;
  double xPosMin;
  double xPosMax;
  double yPosMin;
  double yPosMax;
  double xRange;
  double yRange;
  int imageMax;
  int imageIndex;
  int rangeCondition1;
  int rangeCondition2;
  int rangeCondition3;
  int xPosXImageFlag;
  int mappedImageIndex;
  double xPosValue;
  double yPosValue;
  double posValue;
  int XImageArrayFlag;
  int YImageArrayFlag;
  POVERLAPENTRY pOverlapEntry;
  int debugFlag1 = 0;
  int debugFlag2 = 0;
  int debugFlag3 = 0;
  int debugFlag4 = 0;
  int wideScanLocationPatternIndex;
  int xPosXImageTableIndex;
  PXPOSXIMAGE pXposXImage;


  if (*ppOverlapZone != NULL) {
    printf("ERROR: line %4d GetOverlapZones pOverlapZone is not null\n");
    return(-1);
  }
  

  if (GetMosaicInfo(pConnection,
                    series,
                    plateNumber,
                    mosaicNumber,
                    0,
                    pMosaic) != 1) {
    printf("ERROR: GetOverlapZones can not find mosaic %s%05d_%02d\n",series,plateNumber,mosaicNumber);
    return(-1);
  }
  if (pMosaic->patternID < 0) {
    printf("ERROR: GetOverlapZones can not find a patternID for mosaic %s%05d_%02d\n",series,plateNumber,mosaicNumber);
    return(-1);
  }
  pOverlapZone = (POVERLAPZONE)calloc(1,sizeof(OVERLAPZONE));
  if (pOverlapZone == NULL) {
    printf("ERROR: GetOverlapZones failed allocate the horizontal zone mosaic %s%05d_%02d\n",series,plateNumber,mosaicNumber);
    return(-1);
  }

  pOverlapZone->patternID = pMosaic->patternID;
  if ((pMosaic->transform & TRANSFORM_MIRROR) != 0) {
    pOverlapZone->mirrorFlag = 1;
  }
  if ((pMosaic->transform &  TRANSFORM_FLIP) != 0) {
    pOverlapZone->flipFlag = 1;
  }
  pOverlapZone->rotation = rotation;

  sprintf(queryString,"SELECT tile,xPos,yPos from scanLocations where patternID = %d\n",pMosaic->patternID);
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      pOverlapZone->numTiles =  (int)mysql_num_rows(res_ptr);
      if (pOverlapZone->numTiles == 0) {
        printf("ERROR: GetOverlapZones found no tiles for %s%05d_%02d\n",series,plateNumber,mosaicNumber);
        FreeOverlapZone(&pOverlapZone);
        return(-1);
      }
      pOverlapZone->pScanLocationTable = (PSCANLOCATION)calloc(pOverlapZone->numTiles,sizeof(SCANLOCATION));
      if (pOverlapZone->pScanLocationTable == NULL) {
        printf("ERROR: GetOverlapZones: failed to allocate %d MULTMASK entries\n",pOverlapZone->numTiles);
        FreeOverlapZone(&pOverlapZone);
        return(-1);
      }

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (curTileCount == pOverlapZone->numTiles) {
          printf("ERROR: GetOverlapZones: currentMaskCount %d pOverlapZone->numTiles %d for %s%5d\n",curTileCount,pOverlapZone->numTiles,series,plateNumber);
          FreeOverlapZone(&pOverlapZone);
          return(-1);
        }
        pScanLocation = &pOverlapZone->pScanLocationTable[curTileCount];
        memset(pScanLocation,0,sizeof(SCANLOCATION));
        pScanLocation->patternID = pMosaic->patternID;


        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pScanLocation->tile);
          if (nvals != 1) {
            printf("ERROR: GetOverlapZones: nvals is %d for mosaicNumber %s%5d tile %d\n",nvals,series,plateNumber,curTileCount);
            continue;
          }
        } else {
          printf("ERROR: GetOverlapZones: NULL mosaicNumber %s%5d tile %d \n",series,plateNumber,curTileCount);
          continue;
        }

        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%lf",&pScanLocation->xPos);
          if (nvals != 1) {
            printf("ERROR: GetOverlapZones: nvals is %d for mosaicNumber %s%5d xPos %d\n",nvals,series,plateNumber,curTileCount);
            continue;
          }
        } else {
          printf("ERROR: GetOverlapZones: NULL mosaicNumber %s%5d xPos %d \n",series,plateNumber,curTileCount);
          continue;
        }

        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%lf",&pScanLocation->yPos);
          if (nvals != 1) {
            printf("ERROR: GetOverlapZones: nvals is %d for mosaicNumber %s%5d yPos %d\n",nvals,series,plateNumber,curTileCount);
            continue;
          }
        } else {
          printf("ERROR: GetOverlapZones: NULL mosaicNumber %s%5d yPos %d \n",series,plateNumber,curTileCount);
          continue;
        }
        pOverlapZone->xPosAve += pScanLocation->xPos;
        pOverlapZone->yPosAve += pScanLocation->yPos;
        if (curTileCount == 0) {
          xPosMin = pScanLocation->xPos;
          xPosMax = pScanLocation->xPos;
          yPosMin = pScanLocation->yPos;
          yPosMax = pScanLocation->yPos;
        } else {
          xPosMin = fmin(xPosMin,pScanLocation->xPos);
          xPosMax = fmax(xPosMax,pScanLocation->xPos);
          yPosMin = fmin(yPosMin,pScanLocation->yPos);
          yPosMax = fmax(yPosMax,pScanLocation->yPos);
        }

        curTileCount++;


      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: GetOverlapZones mysql_store_result failed for %s%05d_%02d\n",series,plateNumber,mosaicNumber);
      FreeOverlapZone(&pOverlapZone);
      return(-1);
    }
  }
  if (curTileCount != pOverlapZone->numTiles) {
    printf("ERROR: GetOverlapZones mysql_store_result failed for %s%05d_%02d expected %d tiles, got only %d\n",series,plateNumber,mosaicNumber,pOverlapZone->numTiles,curTileCount);
    FreeOverlapZone(&pOverlapZone);
    return(-1);
  }

  pOverlapZone->xPosAve =  pOverlapZone->xPosAve/(1.*curTileCount);
  pOverlapZone->yPosAve =  pOverlapZone->yPosAve/(1.*curTileCount);
  xRange = xPosMax-xPosMin;
  yRange = yPosMax-yPosMin;
  

#if 1 /* Bugfix of Nov 3, 2020 */
  rangeCondition1 = 0;
  for (wideScanLocationPatternIndex = 0; wideScanLocationPatternIndex < numWideScanLocationPatternID ; wideScanLocationPatternIndex++) {
    if (pMosaic->patternID == wideScanLocationPatternID[wideScanLocationPatternIndex]) {
      rangeCondition1 = 1;
      break;
    }
  }
  if (xRange > yRange) {
    rangeCondition2 = 1;
  } else {
    rangeCondition2 = 0;
  }
  if (pMosaic->naxis1 > pMosaic->naxis2) {
    rangeCondition3 = 1;
    imageMax = pMosaic->naxis1;
  } else {
    rangeCondition3 = 0;
    imageMax = pMosaic->naxis2;
  }
  pOverlapZone->xPosXImageFlag = -1;
  for (xPosXImageTableIndex = 0; xPosXImageTableIndex < xPosXImageTableSize; xPosXImageTableIndex++) {
    pXposXImage = &xPosXImageTable[xPosXImageTableIndex];
    if ((rangeCondition1 ==  pXposXImage->rangeCondition1) &&
        (rangeCondition2 ==  pXposXImage->rangeCondition2) &&
        (rangeCondition3 ==  pXposXImage->rangeCondition3)) {
      pOverlapZone->xPosXImageFlag = pXposXImage->xPosXImageFlag;
      break;
    }
  } 
  if (pOverlapZone->xPosXImageFlag == -1) {
    printf("ERROR: GetOverlapZones line %4d illegal rangeCondition %d %d %d  in  %s%05d_%02d\n",__LINE__,rangeCondition1,rangeCondition2,rangeCondition3,series,plateNumber,mosaicNumber);
    return(-1);
  } else {
    printf("GetOverlapZones xPosImageFlag %d for rangeCondition %d %d %d  in  %s%05d_%02d\n",pXposXImage->xPosXImageFlag,rangeCondition1,rangeCondition2,rangeCondition3,series,plateNumber,mosaicNumber);
  }



#else /* bugfix of Nov 3, 2020  */
  if (xRange > yRange) {
    if (pMosaic->naxis1 > pMosaic->naxis2) {
      pOverlapZone->xPosXImageFlag = 1;
      imageMax = pMosaic->naxis1;
      if ((pOverlapZone->rotation != 0) && (pOverlapZone->rotation != 180)) {
        printf("ERROR: GetOverlapZones line %4d inconsistant rotation %3d for xPosXImageFlag %d in  %s%05d_%02d\n",__LINE__,pOverlapZone->rotation,pOverlapZone->xPosXImageFlag,series,plateNumber,mosaicNumber);
        return(-1);
      }
    } else {
      pOverlapZone->xPosXImageFlag = 0;
      imageMax = pMosaic->naxis2;
      if ((pOverlapZone->rotation != 90) && (pOverlapZone->rotation != 270)) {
        printf("ERROR: GetOverlapZones line %4d inconsistant rotation %3d for xPosXImageFlag %d in  %s%05d_%02d\n",__LINE__,pOverlapZone->rotation,pOverlapZone->xPosXImageFlag,series,plateNumber,mosaicNumber);
        return(-1);
      }
    }
  } else {
    if (pMosaic->naxis1 > pMosaic->naxis2) {
      imageMax = pMosaic->naxis1;
      pOverlapZone->xPosXImageFlag = 0;
      if ((pOverlapZone->rotation != 90) && (pOverlapZone->rotation != 270)) {
        printf("ERROR: GetOverlapZones line %4d inconsistant rotation %3d for xPosXImageFlag %d in  %s%05d_%02d\n",__LINE__,pOverlapZone->rotation,pOverlapZone->xPosXImageFlag,series,plateNumber,mosaicNumber);
        return(-1);
      }
    } else {
      pOverlapZone->xPosXImageFlag = 1;
      imageMax = pMosaic->naxis2;
      if ((pOverlapZone->rotation != 0) && (pOverlapZone->rotation != 180)) {
        printf("ERROR: GetOverlapZones line %4d inconsistant rotation %3d for xPosXImageFlag %d in  %s%05d_%02d\n",__LINE__,pOverlapZone->rotation,pOverlapZone->xPosXImageFlag,series,plateNumber,mosaicNumber);
        return(-1);
      }
    }
  }
#endif /* bugfix of Nov 3, 2020 */
  pOverlapZone->xImageArrayLength = pMosaic->naxis1;
  pOverlapZone->pXImageArray = (int *)calloc(imageMax+1,sizeof(int));
  pOverlapZone->yImageArrayLength = pMosaic->naxis2;
  pOverlapZone->pYImageArray = (int *)calloc(imageMax+1,sizeof(int));
  if (pOverlapZone->xPosXImageFlag != 0) {
    pOverlapZone->maxInspectionTiles = (pMosaic->naxis1 + (INSPECTION_HEIGHT -1))/INSPECTION_HEIGHT;
  } else {
    pOverlapZone->maxInspectionTiles = (pMosaic->naxis2 + (INSPECTION_HEIGHT -1))/INSPECTION_HEIGHT;
  }


  /* Do the transforms */
  if (pOverlapZone->flipFlag  != 0) {
    pOverlapZone->invertYArrayFlag = pOverlapZone->invertYArrayFlag % 1;
  } 
  if (pOverlapZone->mirrorFlag  != 0) {
    pOverlapZone->invertXArrayFlag = pOverlapZone->invertXArrayFlag % 1;
  } 
  {
    int tempXArrayFlag;
    switch (pOverlapZone->rotation) {
    case 0:
      break;

    case 90:
      tempXArrayFlag = pOverlapZone->invertYArrayFlag;
      pOverlapZone->invertYArrayFlag = pOverlapZone->invertXArrayFlag;
      pOverlapZone->invertXArrayFlag = tempXArrayFlag;
      pOverlapZone->invertXArrayFlag = pOverlapZone->invertXArrayFlag % 1;
      break;

    case 180:
      pOverlapZone->invertXArrayFlag = pOverlapZone->invertXArrayFlag % 1;
      pOverlapZone->invertYArrayFlag = pOverlapZone->invertYArrayFlag % 1;
      break;

    case 270:
      tempXArrayFlag = pOverlapZone->invertYArrayFlag;
      pOverlapZone->invertYArrayFlag = pOverlapZone->invertXArrayFlag;
      pOverlapZone->invertXArrayFlag = tempXArrayFlag;
      pOverlapZone->invertYArrayFlag = pOverlapZone->invertYArrayFlag % 1;
      break;

    default:
      printf("ERROR: GetOverlapZones line %4d illegal rotation %3d for xPosXImageFlag %d in  %s%05d_%02d\n",__LINE__,pOverlapZone->rotation,pOverlapZone->xPosXImageFlag,series,plateNumber,mosaicNumber);
      return(-1);
      break;

    }
  }

#if 0
  printf("line %4d xPosAve %f xImageArrayLength %d yPosAve %f yImageArrayLength %d imageMax %d for %s\n",__LINE__, pOverlapZone->xPosAve,pOverlapZone->xImageArrayLength,pOverlapZone->yPosAve,pOverlapZone->yImageArrayLength,imageMax,fileroot);
#endif


  if ((pOverlapZone->pXImageArray == NULL) || (pOverlapZone->pYImageArray == NULL)) {
    printf("ERROR: GetOverlapZones: failed to allocate %x %s IMAGEArray entries\n",pOverlapZone->pXImageArray,pOverlapZone->pYImageArray);
    FreeOverlapZone(&pOverlapZone);
    return(-1);
  }
  /* Now populate the arrays */
  for (imageIndex = 0; imageIndex < imageMax; imageIndex++) {
    for (curTileIndex = 0; curTileIndex < curTileCount; curTileIndex++) {
#if 0
      if (curTileIndex != 10) {
        continue;
      }
#endif
      pScanLocation = &pOverlapZone->pScanLocationTable[curTileIndex];
      XImageArrayFlag = 0;
      YImageArrayFlag = 0;
      if (pOverlapZone->xPosXImageFlag == 1) {
        xPosValue =  pOverlapZone->xPosAve + (NOMINAL_MM_PER_PIXEL * (imageIndex - (pOverlapZone->xImageArrayLength/2) )); 
        yPosValue =  pOverlapZone->yPosAve + (NOMINAL_MM_PER_PIXEL * (imageIndex - (pOverlapZone->yImageArrayLength/2) )); 
        if ((imageIndex < pOverlapZone->xImageArrayLength) &&
            (xPosValue < (pScanLocation->xPos + ((CCD_PIXEL_SIZE/2)*NOMINAL_MM_PER_PIXEL))) &&
            (xPosValue > (pScanLocation->xPos - ((CCD_PIXEL_SIZE/2)*NOMINAL_MM_PER_PIXEL)))) {
          pOverlapZone->pXImageArray[imageIndex]++;
          XImageArrayFlag = 1;
        }
        if ((imageIndex < pOverlapZone->yImageArrayLength) &&
            (yPosValue < (pScanLocation->yPos + ((CCD_PIXEL_SIZE/2)*NOMINAL_MM_PER_PIXEL))) &&
            (yPosValue > (pScanLocation->yPos - ((CCD_PIXEL_SIZE/2)*NOMINAL_MM_PER_PIXEL)))) {
          pOverlapZone->pYImageArray[imageIndex]++;
          YImageArrayFlag = 1;
        }

      } else {
        yPosValue = pOverlapZone->yPosAve + (NOMINAL_MM_PER_PIXEL * (imageIndex - (pOverlapZone->xImageArrayLength/2) )); 
        xPosValue = pOverlapZone->xPosAve + (NOMINAL_MM_PER_PIXEL * (imageIndex - (pOverlapZone->yImageArrayLength/2) )); 
        if ((imageIndex < pOverlapZone->xImageArrayLength) &&
            (yPosValue < (pScanLocation->yPos + ((CCD_PIXEL_SIZE/2)*NOMINAL_MM_PER_PIXEL))) &&
            (yPosValue > (pScanLocation->yPos - ((CCD_PIXEL_SIZE/2)*NOMINAL_MM_PER_PIXEL)))) {
          pOverlapZone->pXImageArray[imageIndex]++;
          XImageArrayFlag = 1;
        }
        if ((imageIndex < pOverlapZone->yImageArrayLength) &&
            (xPosValue < (pScanLocation->xPos + ((CCD_PIXEL_SIZE/2)*NOMINAL_MM_PER_PIXEL))) &&
            (xPosValue > (pScanLocation->xPos - ((CCD_PIXEL_SIZE/2)*NOMINAL_MM_PER_PIXEL)))) {
          pOverlapZone->pYImageArray[imageIndex]++;
          YImageArrayFlag = 1;
        }
      }
#if 0
      {
        int tmpXIndex = pOverlapZone->xImageArrayLength-imageIndex;
        int tmpYIndex = pOverlapZone->yImageArrayLength-imageIndex;
        char tmpChar = ' ';
        if (tmpXIndex < 0) {
          tmpXIndex = 0;
        }
        if (tmpYIndex < 0) {
          tmpYIndex = 0;
        }
        if ((imageIndex == 0) ||
            (imageIndex == (pOverlapZone->xImageArrayLength-1)) ||
            (imageIndex == (pOverlapZone->yImageArrayLength-1)) ||
            (imageIndex == (pOverlapZone->xImageArrayLength/2)) ||
            (imageIndex == (pOverlapZone->yImageArrayLength/2))) {
          tmpChar = '*';
        }
           
            
        

        printf("line %4d imageIndex %5d %5d %5d %c curTileIndex %2d, xPosAve %f xPosValue %f yPosAve %f yPosValue %f xFlag %d %d yFlag %d %d \n",
               __LINE__,
               imageIndex,
               tmpXIndex,
               tmpYIndex,
               tmpChar,
               curTileIndex,
               pOverlapZone->xPosAve,
               xPosValue,
               pOverlapZone->yPosAve,
               yPosValue,
               XImageArrayFlag,
               pOverlapZone->pXImageArray[imageIndex],
               YImageArrayFlag,
               pOverlapZone->pYImageArray[imageIndex]);
      }
#endif
    }
 
  }
#if 0 /* for the following, use: grep "line 1092" los.log | grep "pXImageArray" | awk  '{OFS="\t"}{print $4,$6,$8}' > los.txt */
  for (imageIndex = 0; imageIndex < imageMax; imageIndex++) {
    printf("line %4d imageIndex %5d pXImageArray %5d pYImageArray %5d\n",__LINE__,imageIndex,pOverlapZone->pXImageArray[imageIndex],pOverlapZone->pYImageArray[imageIndex]);
  }
#endif
    
  /* Here we are going to invert the pattern and look for the maximum # of mismatched rows/columns.  If it is large, then we will need to to the mirror, flip, and rotate transforms for X and Y */
  {
    int imageXInvIndex;
    int imageYInvIndex;
    int imageXInvValue;
    int imageYInvValue;
    int deltaX;
    int deltaY;
    int badXPixelCount = 0; 
    int badYPixelCount = 0;
    for (imageIndex = 0; imageIndex < imageMax; imageIndex++) {
      imageXInvIndex = pOverlapZone->xImageArrayLength - imageIndex -1;
      imageYInvIndex = pOverlapZone->yImageArrayLength - imageIndex -1;
      imageXInvValue = pOverlapZone->pXImageArray[imageXInvIndex];
      imageYInvValue = pOverlapZone->pYImageArray[imageYInvIndex];
      deltaX = pOverlapZone->pXImageArray[imageIndex] -pOverlapZone->pXImageArray[imageXInvIndex];
      deltaY = pOverlapZone->pYImageArray[imageIndex] -pOverlapZone->pYImageArray[imageYInvIndex];
      if (imageIndex >= pOverlapZone->xImageArrayLength) {
        imageXInvValue = 0;
        deltaX = -10;
        badXPixelCount = 0;
      }
      if (imageIndex >= pOverlapZone->yImageArrayLength) {
        deltaY = -10;
        imageYInvValue = 0;
        badYPixelCount = 0;
      }
      if (deltaX != 0) {
        badXPixelCount++;
      } else {
#if 0
        if (badXPixelCount >  MAX_MIRROR_ERROR) {
          printf("ERROR: line %4d GetOverlapZones: X-axis index %d badXPixelCount %d exceeds MAX_MIRROR_ERROR for mosaic  %s%05d_%02d\n",__LINE__,imageIndex,badXPixelCount,series,plateNumber,mosaicNumber);
          FreeOverlapZone(&pOverlapZone);
          return(-1);
        }
#endif
        badXPixelCount = 0;
      }
      if (deltaY != 0) {
        badYPixelCount++;
      } else {
#if 0
        if (badYPixelCount >  MAX_MIRROR_ERROR) {
          printf("ERROR: line %4d GetOverlapZones: Y-axis index %d badYPixelCount %d exceeds MAX_MIRROR_ERROR for mosaic  %s%05d_%02d\n",__LINE__,imageIndex,badYPixelCount,series,plateNumber,mosaicNumber);
          FreeOverlapZone(&pOverlapZone);
          return(-1);
        }
#endif
        badYPixelCount = 0;
      }
#if 0
      /* for below use grep "imageXInvIndex" los.log | awk '{print $12,$14,$22,$24}' | sort -u */
      printf("line %4d imageIndex %5d pXImageArray %5d imageXInvIndex %5d pXInvImageArray %5d deltaX %5d countX %d pYImageArray %5d imageYInvIndex %5d pYInvImageArray %5d deltaY %5d countY %d\n",
             __LINE__,
             imageIndex,
             pOverlapZone->pXImageArray[imageIndex],
             imageXInvIndex,
             imageXInvValue,
             deltaX,
             badXPixelCount,
             pOverlapZone->pYImageArray[imageIndex],
             imageYInvIndex,
             imageYInvValue,
             deltaY,
             badYPixelCount);
#endif
      /* Now fill in the zones  */
      /* Start with the horizontal zone */
      if (imageIndex < pOverlapZone->xImageArrayLength) {
        pOverlapEntry = &pOverlapZone->XOverlapEntryTable[pOverlapZone->xZoneCount]; 
        pOverlapEntry->startSextractorIndex = -1;
        pOverlapEntry->endSextractorIndex = -1;
        if (pOverlapZone->invertXArrayFlag == 0) {
          mappedImageIndex = imageIndex;
#if 0
          if (debugFlag1 == 0) {
            printf("line %4d DEBUG mappedImageIndex\n",__LINE__);
            debugFlag1 = 1;
          }
#endif
        } else {
          mappedImageIndex = pOverlapZone->xImageArrayLength - imageIndex -1;
#if 1
          if (debugFlag2 == 0) {
            printf("line %4d DEBUG mappedImageIndex\n",__LINE__);
            debugFlag2 = 1;
          }
#endif
        }
        if (pOverlapZone->xPosXImageFlag == 1) {
          posValue =  pOverlapZone->xPosAve + (NOMINAL_MM_PER_PIXEL * (mappedImageIndex - (pOverlapZone->xImageArrayLength/2) )); 
        } else {
          posValue = pOverlapZone->yPosAve + (NOMINAL_MM_PER_PIXEL * (mappedImageIndex - (pOverlapZone->xImageArrayLength/2) ));
        }
        pOverlapEntry->endPixel = imageIndex;
        pOverlapEntry->endPos = posValue;
        if (pOverlapZone->pXImageArray[imageIndex] != pOverlapEntry->tileOverlapCount) {
          if (pOverlapEntry->tileOverlapCount == 0) {
            pOverlapEntry->tileOverlapCount = pOverlapZone->pXImageArray[imageIndex];
            pOverlapEntry->startPixel = imageIndex;
            pOverlapEntry->startPos = posValue - NOMINAL_MM_PER_PIXEL;
          } else {
            pOverlapEntry->endPixel = imageIndex-1;
            pOverlapEntry->endPos = posValue;
            pOverlapZone->xZoneCount++;
            if (pOverlapZone->xZoneCount >= MAX_TILES) {
              printf("ERROR: line %4d GetOverlapZones: xZoneCount exceeds MAX_TILES for  mosaic  %s%05d_%02d\n",__LINE__,series,plateNumber,mosaicNumber);
              FreeOverlapZone(&pOverlapZone);
              return(-1);
            }
            pOverlapEntry = &pOverlapZone->XOverlapEntryTable[pOverlapZone->xZoneCount];
            pOverlapEntry->startSextractorIndex = -1;
            pOverlapEntry->endSextractorIndex = -1;
            pOverlapEntry->tileOverlapCount = pOverlapZone->pXImageArray[imageIndex];
            pOverlapEntry->startPixel = imageIndex;
            pOverlapEntry->startPos = posValue;
            pOverlapEntry->endPixel = imageIndex; 
            pOverlapEntry->endPos = posValue;

          }
        }
      }
      /* Now do the vertical zone */
      if (imageIndex < pOverlapZone->yImageArrayLength) {
        pOverlapEntry = &pOverlapZone->YOverlapEntryTable[pOverlapZone->yZoneCount]; 
        pOverlapEntry->startSextractorIndex = -1;
        pOverlapEntry->endSextractorIndex = -1;
        if (pOverlapZone->invertYArrayFlag == 0) {
          mappedImageIndex = imageIndex;
#if 0 
          if (debugFlag3 == 0) {
            printf("line %4d DEBUG mappedImageIndex\n",__LINE__);
            debugFlag3 = 1;
          }
#endif
        } else {
          mappedImageIndex = pOverlapZone->yImageArrayLength - imageIndex -1;
#if 1
          if (debugFlag4 == 0) {
            printf("line %4d DEBUG mappedImageIndex\n",__LINE__);
            debugFlag4 = 1;
          }
#endif
        }

        if (pOverlapZone->xPosXImageFlag == 1) {
          posValue =  pOverlapZone->yPosAve + (NOMINAL_MM_PER_PIXEL * (mappedImageIndex - (pOverlapZone->yImageArrayLength/2) )); 
        } else {
          posValue = pOverlapZone->xPosAve + (NOMINAL_MM_PER_PIXEL * (mappedImageIndex - (pOverlapZone->yImageArrayLength/2) ));
        }
        pOverlapEntry->endPixel = imageIndex;
        pOverlapEntry->endPos = posValue;
        if (pOverlapZone->pYImageArray[imageIndex] != pOverlapEntry->tileOverlapCount) {
          if (pOverlapEntry->tileOverlapCount == 0) {
            pOverlapEntry->tileOverlapCount = pOverlapZone->pYImageArray[imageIndex];
            pOverlapEntry->startPixel = imageIndex;
            pOverlapEntry->startPos = posValue - NOMINAL_MM_PER_PIXEL;
          } else {
            pOverlapEntry->endPixel = imageIndex-1;
            pOverlapEntry->endPos = posValue;
            pOverlapZone->yZoneCount++;
            if (pOverlapZone->yZoneCount >= MAX_TILES) {
              printf("ERROR: line %4d GetOverlapZones: yZoneCount exceeds MAX_TILES for  mosaic  %s%05d_%02d\n",__LINE__,series,plateNumber,mosaicNumber);
              FreeOverlapZone(&pOverlapZone);
              return(-1);
            }
            pOverlapEntry = &pOverlapZone->YOverlapEntryTable[pOverlapZone->yZoneCount];
            pOverlapEntry->startSextractorIndex = -1;
            pOverlapEntry->endSextractorIndex = -1;
            pOverlapEntry->tileOverlapCount = pOverlapZone->pYImageArray[imageIndex];
            pOverlapEntry->startPixel = imageIndex;
            pOverlapEntry->startPos = posValue;
            pOverlapEntry->endPixel = imageIndex; 
            pOverlapEntry->endPos = posValue;

          }
        }
      }

    }
    pOverlapZone->xZoneCount++;
    pOverlapZone->yZoneCount++;
#if 0
    if (badXPixelCount >  MAX_MIRROR_ERROR) {
      printf("ERROR: line %4d GetOverlapZones: X-axis index %d badXPixelCount %d exceeds MAX_MIRROR_ERROR for mosaic  %s%05d_%02d\n",__LINE__,imageIndex,badXPixelCount,series,plateNumber,mosaicNumber);
      FreeOverlapZone(&pOverlapZone);
      return(-1);
    }
#endif
#if 0
    if (badYPixelCount >  MAX_MIRROR_ERROR) {
      printf("ERROR: line %4d GetOverlapZones: Y-axis index %d badYPixelCount %d exceeds MAX_MIRROR_ERROR for mosaic  %s%05d_%02d\n",__LINE__,imageIndex,badYPixelCount,series,plateNumber,mosaicNumber);
      FreeOverlapZone(&pOverlapZone);
      return(-1);
    }
#endif
  }




  *ppOverlapZone = pOverlapZone;

  return(0);
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
  char getfitsdirectory[MAX_BUFFER];
  char mosaicFileName[MAX_BUFFER];
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
#ifdef SEARCHTABLE_ONLY
  int savedSearchTableSize = 0;
#endif /* SEARCHTABLE_ONLY */
  double matchRatio2;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int candidateCount = 0;
  int outputCount = 0;
  int populatedBinCount = 0;
  int psfsaturatedCount = 0;
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
  size_t vector2Size;
  size_t vectorIndex;
  int sextractorIndex;
  int sextractorIndex1;
  int sextractorIndex2;
  int minSextractorIndex;
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
  int *correlation_table1 = NULL;
  int *correlation_table2 = NULL;
  int *correlation_histogram2 = NULL;
  int *correlation_table3 = NULL;
  int correlationTablePixels;
#if 0
  int minX;
  int maxX;
  int minY;
  int maxY;
#endif
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
  int doPlots = 0;
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
  double * vector2 = NULL;
  /* This is the cutoff for the ratio (FLUX_MAX/THRESHOLD) */
  double analysis_threshold = USE_ANALYSIS_THRESHOLD;
  double MAG_ISO_crit;
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
  int totalMasksSector = 0;
  int debugPrint = 0;
  int debugCount = 0;
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
  int int_XY_IMAGE;
  int totalMaskSectors;
  int maskSector;
  int minSectorCount;
  int maxSectorCount; 
  int minUnsaturatedCount;
  int maxUnsaturatedCount;
  int BFLAGSMASK = (~((1 << FILTER_BFLAG_NEIGHBORS) | (1 << FILTER_BFLAG_BLEND)));
  double pmf2;
  double meanArrival;
  double factorial = 1.;
  double prefix;
  double exponential;
  double matchRatio;
  int overflowFlag;
  int overflowPrint;
  POVERLAPZONE pOverlapZone = NULL;
  int result;
  MOSAIC mosaicInfo;
  PMOSAIC pMosaic = &mosaicInfo;
  POVERLAPENTRY pOverlapEntryTable = NULL;
  POVERLAPENTRY pOverlapEntry = NULL;
  int maskResult = 0;
  int maskCount = 0;
  int maskFlag = 0;
  int maskCenterFlag = 0;
  PMULTMASK pMultipleMask = NULL;
  PMULTMASK pMultipleTable = NULL;

  /* The following is used for getfits */
  int inspectionTileIndex;
  int bestInspectionTileIndex;
  int bestInspectionTileLower;
  int bestInspectionTileCount;
  int maxTileOverlapCount = 0;
  int cutoffTileOverlapCount = 0; /* used to distinguish between single and dual overlap areas */
  int maxSingleMaskArea = 0;              /* highest object density when masks are generated */
  int maxSingleMaskAreaMaskSector = -1;
  int maxDualMaskArea = 0;              /* highest object density when masks are generated */
  int maxDualMaskAreaMaskSector = -1;
  int highestDualMaskInspectionCount = 0; /* highest object density when no masks are generated */
  int highestDualMaskMaskSector = -1;
  int directoryLength;

  outfile[0] = 0;
  getfitsdirectory[0] = 0;
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

#ifdef DEBUG_NO_MYSQL
  printf("ERROR: DEBUG_NO_MYSQL is set \n");
#endif /* DEBUG_NO_MYSQL */

#ifdef SKIP_POISSON
  printf("ERROR: SKIP_POISSON is set \n");
#endif /* SKIP_POISSON */

#ifdef DEBUG_DELETE_MASKS
  printf("ERROR: DEBUG_DELETE_MASKS is set \n");
#endif /* DEBUG_DELETE_MASKS */

#ifdef DEBUG_VERT4X5
  printf("ERROR: DEBUG_VERT4X5 is set \n");
#endif /* DEBUG_VERT4X5 */




#ifdef SEARCHTABLE_ONLY
  printf("ERROR: SEARCHTABLE_ONLY is set SEARCH_TABLE_LENGTH %d \n",SEARCH_TABLE_LENGTH);
#endif /* SEARCHTABLE_ONLY */

#ifdef REJECT_GRAIN
  printf("ERROR: REJECT_GRAIN is set \n");
#endif /* REJECT_GRAIN */

#ifdef COUNT_LIMIT1
  printf("ERROR: COUNT_LIMIT1 is set to %d\n",COUNT_LIMIT1);
#endif /* COUNT_LIMIT1 */

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

        case 'g': /* getfits directory */
        case 'G':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            directoryLength = strlen(*++argv);
            if (directoryLength >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
            strncpy(getfitsdirectory,*argv,MAX_BUFFER-2);
            if (getfitsdirectory[directoryLength-1] != '/') {
              strcat(getfitsdirectory,"/");
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
          doPlots = 1;
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

   
  if (doPlots) {
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
    printf("Usage: search_close3\n");
    printf("                 -b  <binning> Set to 16 for regions on thumbnail mosaics\n");
    printf("                 -d  <background file> <region file>search for calibration squares\n");
    printf("                 -h  <mosaic height>\n");
    printf("                 -i <sextractor file> \n");
    printf("                 -o  <output file>\n");
    printf("                 -o  <getfits directory>\n");
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


  if (verbose) {
    printf("search_close3 of %s %s \n Input Filename %s\n Output Filename %s mosaic width: %d mosaic height: %d scale: %f arcsec/pixel:",
           __DATE__,__TIME__,sextractor_name,outfile,
           mosaicWidth,
           mosaicHeight,
           arcsecPerPixel);
    printf(" threshold: %f",analysis_threshold);
    printf("\n");
    printf("search_close3 BFLAGSMASK is %d 0x%x ~BFLAGSMASK is %d 0x%x\n",BFLAGSMASK,BFLAGSMASK,~BFLAGSMASK,~BFLAGSMASK);

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
  /* see if this is already a close multiple exposure */
  maskResult = GetMultipleMask(pConnection,series,plateNumber,mosaicNumber,&maskCount,&pMultipleTable);
  if (maskCount > 0) {
    for (maskIndex = 0; maskIndex < maskCount; maskIndex++) {
      pMultipleMask = &pMultipleTable[maskIndex];
      if (pMultipleMask->maskCenterFlag != 0) {
        maskCenterFlag = 1;
      } else {
        maskFlag = 1;
      }
    }
    if ((maskCenterFlag != 0) || (maskFlag != 0)) {
      printf("line %4d search_close3 QUALITY_MULTIPLE maskCount %2d maskFlag %d maskCenterFlag %d for plate %s\n",__LINE__,maskCount,maskFlag,maskCenterFlag,fileroot);
    }

    if (pMultipleTable != NULL) {
      free(pMultipleTable);
      pMultipleTable= NULL;
    }
  }



#ifndef DEBUG_TRUNCATION
#ifndef DEBUG_NO_MYSQL
  /* Clear the multiple mask status */
  SetMosaicFitWCS(pConnection,"NoMultipleMask",fileroot,0,0,readFlag,&FitWCS);
  SetMosaicFitWCS(pConnection,"HaveMultipleMask",fileroot,0,0,readFlag,&FitWCS);

  
  sprintf(queryString,"DELETE FROM masks2 where series = '%s' and plateNumber = %d and mosaicNumber = %d;",
          series,
          plateNumber,
          mosaicNumber);
  if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
    printf("ERROR: search_close3 MAX_QUERY_STRING exceeded\n");
    exit(-1);
  }
  res = ExecuteQuery(pConnection,queryString);
  if (res) {
    exit(-1);
  }
#endif /* DEBUG_NO_MYSQL */
#endif /* DEBUG_TRUNCATION */
#ifndef DEBUG_NO_MYSQL
  sprintf(queryString,"DELETE FROM masks3 where series = '%s' and plateNumber = %d and mosaicNumber = %d;",
          series,
          plateNumber,
          mosaicNumber);
  if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
    printf("ERROR: search_close3 MAX_QUERY_STRING exceeded\n");
    exit(-1);
  }
  res = ExecuteQuery(pConnection,queryString);
  if (res) {
    exit(-1);
  }
#endif /* DEBUG_NO_MYSQL */

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
  result = GetOverlapZones(pConnection,series,plateNumber,mosaicNumber,tempbinning,rotation,pMosaic,&pOverlapZone);
  if (result != 0) {
    printf("ERROR: line %4d failed to get overlap zones for %s\n",__LINE__,fileroot);
    FreeOverlapZone(&pOverlapZone);
    exit(-1);
  }
  if ((pMosaic->naxis1 != mosaicWidth) || (pMosaic->naxis2 != mosaicHeight)) {
    printf("ERROR: line %4d mosaic width %d %d or mosaic height %d %d in error for  %s\n",__LINE__,pMosaic->naxis1,mosaicWidth,pMosaic->naxis2,mosaicHeight,fileroot);
    FreeOverlapZone(&pOverlapZone);
    exit(-1);
  }
  if (pOverlapZone->xPosXImageFlag == 1) {
    totalMaskSectors = pOverlapZone->yZoneCount;
    pOverlapEntryTable = &pOverlapZone->YOverlapEntryTable[0];
  } else {
    totalMaskSectors = pOverlapZone->xZoneCount;
    pOverlapEntryTable = &pOverlapZone->XOverlapEntryTable[0];
  } 
  if ((totalMaskSectors) > MAX_SECTOR_COUNT) {
    printf("ERROR: MAX_SECTOR_COUNT %d is too large for %s\n",totalMaskSectors,fileroot);
    exit(-1);
  }
 
  maskSectorTable = (PMASKSECTOR)calloc(totalMaskSectors+1,sizeof(MASKSECTOR));
  if (maskSectorTable == NULL) {
    printf("ERROR: failed to allocate the maskSectorStable of size %d\n",totalMaskSectors);
    exit(-1);
  }
  for (maskSector = 0; maskSector < totalMaskSectors; maskSector++) {
    pMaskSector = &maskSectorTable[maskSector];
    pMaskSector->sectorFlag = MASKSECTOR_FLAG;
    pMaskSector->MAG_ISO_med = -99.0;
    pMaskSector->MAG_ISO_rms = -99.0;     
    pOverlapEntry = &pOverlapEntryTable[maskSector];
    if (pOverlapEntry->tileOverlapCount >  maxTileOverlapCount) {
      maxTileOverlapCount = pOverlapEntry->tileOverlapCount;
    }

    if (pOverlapEntry->inspectionTileTableLower != NULL) {
      printf("ERROR: line %4d inspectionTileTableLower already allocated for %s\n",__LINE__,fileroot);
      exit(-1);
    }
    pOverlapEntry->inspectionTileTableLower = (int *)calloc(pOverlapZone->maxInspectionTiles,sizeof(int));
    if (pOverlapEntry->inspectionTileTableLower == NULL) {
      printf("ERROR: line %4d inspectionTileTableLower allocation failure for %s\n",__LINE__,fileroot);
      exit(-1);
    }

    if (pOverlapEntry->inspectionTileTableUpper != NULL) {
      printf("ERROR: line %4d inspectionTileTableUpper already allocated for %s\n",__LINE__,fileroot);
      exit(-1);
    }
    pOverlapEntry->inspectionTileTableUpper = (int *)calloc(pOverlapZone->maxInspectionTiles,sizeof(int));
    if (pOverlapEntry->inspectionTileTableUpper == NULL) {
      printf("ERROR: line %4d inspectionTileTableUpper allocation failure for %s\n",__LINE__,fileroot);
      exit(-1);
    }
  }
  /* Calculate the spatial bin and maskSector for each record */
  for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
    pInputImage = &sextractor_table[sextractorIndex];
    pInputImage->spatial_bin =  CalculateBin(mosaicWidth,mosaicHeight,pInputImage->X_IMAGE,pInputImage->Y_IMAGE,&edgeDist);
    if (pOverlapZone->xPosXImageFlag == 1) {
      int_XY_IMAGE = pInputImage->Y_IMAGE;
    } else {
      int_XY_IMAGE = pInputImage->X_IMAGE;
    }
    for (maskSector = 0; maskSector < totalMaskSectors; maskSector++) {
      pOverlapEntry = &pOverlapEntryTable[maskSector];
      if ((int_XY_IMAGE >= pOverlapEntry->startPixel) &&
          (int_XY_IMAGE <= pOverlapEntry->endPixel)) {
        break;
      }
    } 
    if (maskSector >= totalMaskSectors) {
      pInputImage->maskSector = -1;
      pOverlapZone->unassignedImageCount++;
    } else { 
      pOverlapEntry->imageCount++;
      pInputImage->maskSector = maskSector;
      pMaskSector = &maskSectorTable[pInputImage->maskSector];
      if (pMaskSector->sectorFlag != MASKSECTOR_FLAG) {
        printf("ERROR: line %4d Illegal mask sector %d for %s\n",__LINE__,pInputImage->maskSector,fileroot);
        exit(-1);
      }
      pMaskSector->sectorCount++;
      pMaskSector = NULL;
    }
  }


#if 0
  printf("line %4d DEBUG EXIT\n",__LINE__);
  exit(-1);
  int nxs;
  int nys;
  int ixs;
  int iys;
#endif

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
  /* Sort the table in increasing Y_IMAGE values */
  qsort(sextractor_table,sextractor_nrecs,sizeof(INPUTIMAGE),yCompare);
#if 0
  for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
    pInputImage = &sextractor_table[sextractorIndex];
    printf("line %4d Entry %7d NUMBER %6d X_IMAGE %12.6f Y_IMAGE %12.6f maskSector %3d\n",
           __LINE__,
           sextractorIndex,
           pInputImage->NUMBER,
           pInputImage->X_IMAGE,
           pInputImage->Y_IMAGE,
           pInputImage->maskSector);
  }
#endif


  /* Now copy the table over, setting the FILTER_BFLAG_PSFSATURATED flag */
  image_table = (PSTARIMAGE)calloc(sextractor_nrecs,sizeof(STARIMAGE));
  if (image_table == NULL) {
    printf("ERROR: failed to allocate image_table of size %d\n",sextractor_nrecs);
    exit(-1);
  }
  for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
    pInputImage = &sextractor_table[sextractorIndex];
    pSextractor3 = &image_table[sextractorIndex];
    pSextractor3->sextractorIndex = sextractorIndex;
    pSextractor3->maskSector = -1;

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
    pSextractor3->maskSector = pInputImage->maskSector;

    pOverlapEntry = &pOverlapEntryTable[pInputImage->maskSector];
    if (pOverlapEntry->startSextractorIndex < 0) {
      pOverlapEntry->startSextractorIndex = sextractorIndex;
    }
    pOverlapEntry->endSextractorIndex = sextractorIndex;

    /* Set the FILTER_BFLAG_PSFSATURATED flag */
    if ((pInputImage->ISO4 > 0) &&
        (((1.0*pInputImage->ISO5)/(1.0*pInputImage->ISO4)) > PSFSATURATED_ISO) &&
        ((pInputImage->BACKGROUND + pInputImage->FLUX_MAX) > PSFSATURATED_FLUX)) {
#if 0
      printf("%d %f %f\n",pInputImage->NUMBER,(1.0*pInputImage->ISO5)/(1.0*pInputImage->ISO4),pInputImage->BACKGROUND + pInputImage->FLUX_MAX);
#endif
      pSextractor3->BFLAGS |= (1 << FILTER_BFLAG_PSFSATURATED);
      psfsaturatedCount++;
      if ((pSextractor3->maskSector >= 0) &&
          (pSextractor3->maskSector < totalMaskSectors)) {
        pMaskSector2 = &maskSectorTable[pSextractor3->maskSector];
        if (pMaskSector2->sectorFlag != MASKSECTOR_FLAG) {
          printf("ERROR: line %4d Illegal mask sector %d for %s\n",__LINE__,pSextractor3->maskSector,fileroot);
          exit(-1);
        }

        pMaskSector2->unsaturatedCount--;
        pMaskSector2 = NULL;
      }
    }
  }
#if 0
  for (insertIndex = 0; insertIndex < sextractor_nrecs; insertIndex++) {
    pSextractor2 = &image_table[insertIndex];
    printf("line %4d Entry %7d %7d NUMBER %6d X_IMAGE %12.6f Y_IMAGE %12.6f maskSector %3d\n",
           __LINE__,
           insertIndex,
           pSextractor2->sextractorIndex,
           pSextractor2->NUMBER,
           pSextractor2->X_IMAGE,
           pSextractor2->Y_IMAGE,
           pSextractor2->maskSector);
  }
  exit(-1);
#endif
  


  maskSector = 0;
  pMaskSector = &maskSectorTable[maskSector];
  if (pMaskSector->sectorFlag != MASKSECTOR_FLAG) {
    printf("ERROR: line %4d Illegal mask sector %d for %s\n",__LINE__,maskSector,fileroot);
    exit(-1);
  }
 
  minSectorCount = pMaskSector->sectorCount;
  maxSectorCount = pMaskSector->sectorCount;
  minUnsaturatedCount = pMaskSector->unsaturatedCount+pMaskSector->sectorCount;
  maxUnsaturatedCount = pMaskSector->unsaturatedCount+pMaskSector->sectorCount;
  vector2 = (double *)calloc(sextractor_nrecs,sizeof(double));
  if (vector2 == NULL) {
    printf("ERROR: line %4d search_close3 failed to allocate vector2\n");
    exit(-1);
  }
  for (maskSector = 0; maskSector < totalMaskSectors; maskSector++) {
    pMaskSector = &maskSectorTable[maskSector];
    pOverlapEntry = &pOverlapEntryTable[maskSector]; 
    if (pMaskSector->sectorFlag != MASKSECTOR_FLAG) {
      printf("ERROR: line %4d Illegal mask sector %d for %s\n",__LINE__,maskSector,fileroot);
      exit(-1);
    }
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
    /* Here we use the filter_wedge.c algorithm to calculate MAG_ISO_med, MAG_ISO_rms, and curClipCount
       such that (pSextractorBright->MAG_ISO < (MAG_ISO_med - MAG_ISO_rms)) eliminates all of the grain noise 
    */
    vector2Size = 0;
    for (sextractorIndex = pOverlapEntry->startSextractorIndex; sextractorIndex < pOverlapEntry->endSextractorIndex; sextractorIndex++) {
      pInputImage = &sextractor_table[vector2Size];
      vector2[vector2Size] = pInputImage->MAG_ISO;
      vector2Size++;
    }
    pMaskSector->curClipCount = vector2Size;
    pMaskSector->curClipCount = CalcMedianAndRMS(pMaskSector->curClipCount,10,vector2,&pMaskSector->MAG_ISO_med,&pMaskSector->MAG_ISO_rms,1,CLIP_RMS_FACTOR,0);
    if (pMaskSector->curClipCount <= 0) {
      pMaskSector->MAG_ISO_med = -99.0;
      pMaskSector->MAG_ISO_rms = -99.0;
    }
    if (verbose) {
      printf("line %4d maskSector %2d startSextractorIndex %6d endSextractorIndex %6d vector2Size %6d curClipCount %6d MAG_ISO median %f, MAG_ISO rms %f\n",
             __LINE__,
             maskSector,
             pOverlapEntry->startSextractorIndex,
             pOverlapEntry->endSextractorIndex,
             vector2Size,
             pMaskSector->curClipCount,
             pMaskSector->MAG_ISO_med,
             pMaskSector->MAG_ISO_rms);
    }


  }
  if (vector2 = NULL) {
    free(vector2);
    vector2 = NULL;
  }

  pMaskSector = NULL;
  if (verbose) { 
    printf("line %4d search_close3 totalMaskSectors %2d  minSectorCount: %5d maxSectorCount %5d minUnsaturedCount %5d maxUnsaturedCount %5d flip %d mirror %d rotation %3d for %s\n",
           __LINE__,
           totalMaskSectors,
           minSectorCount,
           maxSectorCount,
           minUnsaturatedCount,
           maxUnsaturatedCount,
           pOverlapZone->flipFlag,
           pOverlapZone->mirrorFlag,
           pOverlapZone->rotation,
           fileroot);
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
  correlationTablePixels = 4*(MAX_PIXEL_RADIUS+1)*(MAX_PIXEL_RADIUS+1);
  correlation_table1 = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_table1 == NULL) {
    printf("ERROR: failed to allocate correlation table of size %d\n",correlationTablePixels);
  }
  correlation_table2 = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_table2 == NULL) {
    printf("ERROR: failed to allocate correlation table of size %d\n",correlationTablePixels);
  }
  correlation_histogram2 = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_histogram2 == NULL) {
    printf("ERROR: failed to allocate correlation histogram of size %d\n",correlationTablePixels);
  }
  correlation_table3 = (int*)calloc(correlationTablePixels,sizeof(int));
  if (correlation_table3 == NULL) {
    printf("ERROR: failed to allocate correlation table of size %d\n",correlationTablePixels);
  }
  /* BEGINNING OF SECTOR LOOP */
  for (maskSector = 0; maskSector < totalMaskSectors; maskSector++) {
    double minXRange;
    double maxXRange;
    double minYRange;
    double maxYRange;
    pMaskSector = &maskSectorTable[maskSector];
    if (pMaskSector->sectorFlag != MASKSECTOR_FLAG) {
      printf("ERROR: line %4d Illegal mask sector %d for %s\n",__LINE__,maskSector,fileroot);
      exit(-1);
    }

    pOverlapEntry = &pOverlapEntryTable[maskSector]; 
    cutoffTileOverlapCount = (maxTileOverlapCount/2)+1; /* Decide whether this is a single or dual overlap (See parselogs.c line 1803) */
#if DEBUG_TRUNCATION
    if (maskSector == 0) {
      maskSector = 3;
    } else {
      printf("ERROR: DEBUG_TRUNCATION abort in line %4d\n",__LINE__);
      exit(-1);
    }

#endif /* DEBUG_TRUNCATION */
    if (verbose) {
      printf("line %4d search_close3 Processing sector %2d for %s\n",__LINE__,maskSector,fileroot);
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
    memset(correlation_table1,0,correlationTablePixels*sizeof(int));
    memset(correlation_table2,0,correlationTablePixels*sizeof(int));
    memset(correlation_table3,0,correlationTablePixels*sizeof(int));
    correlation_sum0 = 0;
    correlation_sum1 = 0;
    pmf_sum0 = 0.0;
    bin_sum0 = 0.0;
    for (searchIndex = 0; searchIndex < SEARCH_TABLE_LENGTH; searchIndex++) {
      pSearchEntry1 = &searchTable[searchIndex];
      memset(pSearchEntry1,0,sizeof(SEARCHENTRY));
      pSearchEntry1->maskSector = -1;
      pSearchEntry1->sextractorIndex = -1;
    }
    pMaskSector->aveXY = (pOverlapEntry->startPixel+pOverlapEntry->endPixel)/2.0;
    if (verbose) {
      printf("line %4d search_close3 Sector %2d Min, Max XY: %f %f aveXY %f\n",__LINE__,maskSector,pOverlapEntry->startPixel,pOverlapEntry->endPixel,pMaskSector->aveXY);
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

    for (sextractorIndex = pOverlapEntry->startSextractorIndex; sextractorIndex < pOverlapEntry->endSextractorIndex; sextractorIndex++) {
      pSextractorBright = &image_table[sextractorIndex];

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

#ifdef REJECT_GRAIN
      if (pSextractorBright->MAG_ISO > (pMaskSector->MAG_ISO_med - pMaskSector->MAG_ISO_rms)) {
#if 0
        printf("line %4d Rejecting Bright Object %6d MAG_ISO %10f FWHM_WORLD %d limit %10f\n",
               __LINE__,
               sextractorIndex,
               pSextractorBright->MAG_ISO,
               pSextractorBright->FWHM_IMAGE,
               pMaskSector->MAG_ISO_med - pMaskSector->MAG_ISO_rms);
#endif
        continue;
      }
#endif /* REJECT_GRAIN */      

      if ((pSextractorBright->FLUX_MAX/pSextractorBright->THRESHOLD) < analysis_threshold) {
        continue;
      }

      if (pOverlapZone->xPosXImageFlag == 1) {
        int_XY_IMAGE = pSextractorBright->Y_IMAGE;
      } else {
        int_XY_IMAGE = pSextractorBright->X_IMAGE;
      }

      if ((int_XY_IMAGE < pOverlapEntry->startPixel) ||
          (int_XY_IMAGE > pOverlapEntry->endPixel)) {
#if 1
        printf("line %4d IMAGE ON EDGE int_XY_IMAGE: %10d for maskSector %d in %s \n",__LINE__,int_XY_IMAGE,maskSector,fileroot);
#endif
        continue;
      }

      /* Insert this object in the search table */
      if ((searchTableSize < SEARCH_TABLE_LENGTH) ||
          ((searchTable[SEARCH_TABLE_LENGTH-1].sextractorIndex >= 0) && (pSextractorBright->MAG_ISO < searchTable[SEARCH_TABLE_LENGTH-1].MAG_ISO))) {
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
            pSearchEntry1->maskSector  = pSextractorBright->maskSector;
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
              pSearchEntry1->maskSector  = pSextractorBright->maskSector;
              if (searchTableSize < SEARCH_TABLE_LENGTH) {
                searchTableSize++;
              }
              break;
            }
          }
        } /* Insertion loop */
      } /* Star bright enough for insertion */
#if 0
      if (sextractorIndex > (2 * SEARCH_TABLE_LENGTH)) {
        printf("ERROR: early abort\n");
        exit(-1);
      }
#endif
#if 0
      if (maskSector == 1) {
        printf("line %4d maskSector %2d sextractorIndex %5d searchTableSize %4d\n",__LINE__,maskSector,sextractorIndex,searchTableSize);
        for (insertIndex = 0; insertIndex < SEARCH_TABLE_LENGTH; insertIndex++) {
          pSearchEntry2 = &searchTable[insertIndex];
          if (pSearchEntry2->sextractorIndex >= 0) {
            pSextractor2 = &image_table[pSearchEntry2->sextractorIndex];
            printf("line %4d Entry %7d  Index %7d MAG_ISO %10f NUMBER %6d X_IMAGE %12.6f Y_IMAGE %12.6f searchStar %d maskSector %2d\n",
                   __LINE__,
                   insertIndex,
                   pSearchEntry2->sextractorIndex,
                   pSearchEntry2->MAG_ISO,
                   pSextractor2->NUMBER,
                   pSextractor2->X_IMAGE,
                   pSextractor2->Y_IMAGE,
                   pSextractor2->searchStar,
                   pSextractor2->maskSector);

          }
        }
      }
#endif

    } /* Sextractor table index */
    if (verbose) {
      printf("line %4d search_close3 Sector %2d searchTableSize %d\n",__LINE__,maskSector,searchTableSize);
    }
#ifdef SEARCHTABLE_ONLY
#if 0 /* Remove on Apr 16, 2020 - need get even low quality data */
    savedSearchTableSize = searchTableSize;
    if (searchTableSize < SEARCH_TABLE_LENGTH) {
      searchTableSize = 0;
    }
#endif
#endif /* SEARCHTABLE_ONLY */
      
    qsort(searchTable,searchTableSize,sizeof(SEARCHENTRY),searchCompare);


    for (searchIndex = 0; searchIndex < searchTableSize; searchIndex++) {
      pSearchEntry1 = &searchTable[searchIndex];
      /* See if this image is in an inspection tile */
      if (pOverlapZone->xPosXImageFlag == 1) {
        
        if ((pSearchEntry1->Y_IMAGE - pOverlapEntry->startPixel) < (INSPECTION_WIDTH/2)) {
          inspectionTileIndex = pSearchEntry1->X_IMAGE/INSPECTION_HEIGHT;
          if (inspectionTileIndex < pOverlapZone->maxInspectionTiles) {
            pOverlapEntry->inspectionTileTableLower[inspectionTileIndex]++;
            if ((pOverlapEntry->tileOverlapCount > cutoffTileOverlapCount) &&
                (pOverlapEntry->inspectionTileTableLower[inspectionTileIndex] > highestDualMaskInspectionCount)) {
              highestDualMaskInspectionCount = pOverlapEntry->inspectionTileTableLower[inspectionTileIndex];
              highestDualMaskMaskSector = maskSector;
            }
#if 0
            /* The following grep "line 2935 inspection " los.log | awk '{OFS="\t"}{print $5,$7,$9,$11}'  */
            printf("line %4d inspection X_IMAGE %f Y_IMAGE %f maskSector %2d sequence %2d for %s\n",
                   __LINE__,
                   pSearchEntry1->X_IMAGE,
                   pSearchEntry1->Y_IMAGE,
                   maskSector,
                   (2*maskSector),
                   fileroot);
#endif
          } else {
            printf("ERROR: line %4d inspectionTileTable index failure for %s\n",__LINE__,fileroot);
          }
        }
        if ((pOverlapEntry->endPixel - pSearchEntry1->Y_IMAGE) < (INSPECTION_WIDTH/2)) {
          inspectionTileIndex = pSearchEntry1->X_IMAGE/INSPECTION_HEIGHT;
          if (inspectionTileIndex < pOverlapZone->maxInspectionTiles) {
            pOverlapEntry->inspectionTileTableUpper[inspectionTileIndex]++;
            if ((pOverlapEntry->tileOverlapCount > cutoffTileOverlapCount) &&
                (pOverlapEntry->inspectionTileTableLower[inspectionTileIndex] > highestDualMaskInspectionCount)) {
              highestDualMaskInspectionCount = pOverlapEntry->inspectionTileTableLower[inspectionTileIndex];
              highestDualMaskMaskSector = maskSector;
            }
#if 0
            printf("line %4d inspection X_IMAGE %f Y_IMAGE %f maskSector %2d sequence %2d for %s\n",
                   __LINE__,
                   pSearchEntry1->X_IMAGE,
                   pSearchEntry1->Y_IMAGE,
                   maskSector,
                   (2*maskSector)+1,
                   fileroot);
#endif
          } else {
            printf("ERROR: line %4d inspectionTileTable index failure for %s\n",__LINE__,fileroot);
          }
        }
      } else {
        if ((pSearchEntry1->X_IMAGE - pOverlapEntry->startPixel) < (INSPECTION_WIDTH/2)) {
          inspectionTileIndex = pSearchEntry1->Y_IMAGE/INSPECTION_HEIGHT;
          if (inspectionTileIndex < pOverlapZone->maxInspectionTiles) {
            pOverlapEntry->inspectionTileTableLower[inspectionTileIndex]++;
            if ((pOverlapEntry->tileOverlapCount > cutoffTileOverlapCount) &&
                (pOverlapEntry->inspectionTileTableLower[inspectionTileIndex] > highestDualMaskInspectionCount)) {
              highestDualMaskInspectionCount = pOverlapEntry->inspectionTileTableLower[inspectionTileIndex];
              highestDualMaskMaskSector = maskSector;
            }
#if 0
            printf("line %4d inspection X_IMAGE %f Y_IMAGE %f maskSector %2d sequence %2d for %s\n",
                   __LINE__,
                   pSearchEntry1->X_IMAGE,
                   pSearchEntry1->Y_IMAGE,
                   maskSector,
                   (2*maskSector),
                   fileroot);
#endif
          } else {
            printf("ERROR: line %4d inspectionTileTable index failure for %s\n",__LINE__,fileroot);
          }
        }
        if ((pOverlapEntry->endPixel - pSearchEntry1->X_IMAGE) < (INSPECTION_WIDTH/2)) {
          inspectionTileIndex = pSearchEntry1->Y_IMAGE/INSPECTION_HEIGHT;
          if (inspectionTileIndex < pOverlapZone->maxInspectionTiles) {
            pOverlapEntry->inspectionTileTableUpper[inspectionTileIndex]++;
            if ((pOverlapEntry->tileOverlapCount > cutoffTileOverlapCount) &&
                (pOverlapEntry->inspectionTileTableLower[inspectionTileIndex] > highestDualMaskInspectionCount)) {
              highestDualMaskInspectionCount = pOverlapEntry->inspectionTileTableLower[inspectionTileIndex];
              highestDualMaskMaskSector = maskSector;
            }
#if 0
            printf("line %4d inspection X_IMAGE %f Y_IMAGE %f maskSector %2d sequence %2d for %s\n",
                   __LINE__,
                   pSearchEntry1->X_IMAGE,
                   pSearchEntry1->Y_IMAGE,
                   maskSector,
                   (2*maskSector)+1,
                   fileroot);
#endif
          } else {
            printf("ERROR: line %4d inspectionTileTable index failure for %s\n",__LINE__,fileroot);
          }
        }


      }

      for (sextractorIndex = pOverlapEntry->startSextractorIndex; sextractorIndex < pOverlapEntry->endSextractorIndex; sextractorIndex++) {
        pSextractorDim = &image_table[sextractorIndex];
        if (pSearchEntry1->NUMBER == pSextractorDim->NUMBER) {
          pSextractorDim->searchStar = 1;
          break;
        }
      }
    }



#if 0
    for (insertIndex = 0; insertIndex < searchTableSize; insertIndex++) {
      pSearchEntry2 = &searchTable[insertIndex];
      if (pSearchEntry2->sextractorIndex >= 0) {
        pSextractor2 = &image_table[pSearchEntry2->sextractorIndex];
        printf("line %4d Entry %6d  Index %7d %7d MAG_ISO %10f NUMBER %6d X_IMAGE %12.6f Y_IMAGE %12.6f searchStar %d maskSector %2d %2d\n",
               __LINE__,
               insertIndex,
               pSearchEntry2->sextractorIndex,
               pSextractor2->sextractorIndex,
               pSearchEntry2->MAG_ISO,
               pSextractor2->NUMBER,
               pSextractor2->X_IMAGE,
               pSextractor2->Y_IMAGE,
               pSextractor2->searchStar,
               pSearchEntry2->maskSector,
               pSextractor2->maskSector);

      }
    }
#endif

    /* Now go through the sextractor list for each of these objects and save likely candidates */
    minSextractorIndex = pOverlapEntry->startSextractorIndex;
    pSextractor2 = NULL;
    for (searchIndex = 0; searchIndex < searchTableSize; searchIndex++) {
      pSearchEntry1 = &searchTable[searchIndex];

      NUMBER = image_table[pSearchEntry1->sextractorIndex].NUMBER;

      minXRange = pSearchEntry1->X_IMAGE  - MAX_PIXEL_RADIUS;
      maxXRange = pSearchEntry1->X_IMAGE  + MAX_PIXEL_RADIUS;
      minYRange = pSearchEntry1->Y_IMAGE  - MAX_PIXEL_RADIUS; maxYRange = pSearchEntry1->Y_IMAGE  + MAX_PIXEL_RADIUS;

      for (sextractorIndex = minSextractorIndex; sextractorIndex < pOverlapEntry->endSextractorIndex; sextractorIndex++) {  
        pSextractorDim = &image_table[sextractorIndex];

        debugCount++;
#if 0
        if (debugCount > (147554747-1)) {
          printf("line %4d sextractorIndex min %7d loop %7d Dim %7d End %7d Y_IMAGE %12.6f maskSector %2d ",
                 __LINE__,
                 minSextractorIndex,
                 sextractorIndex,
                 pSextractorDim->sextractorIndex,
                 pOverlapEntry->endSextractorIndex,
                 pSextractorDim->Y_IMAGE,
                 pSextractorDim->maskSector);
          if (pSextractor2 != NULL) {
            printf("sextractorIndex2 %7d Y_IMAGE %12.6f maskSector %2d\n",
                   pSextractor2->sextractorIndex,
                   pSextractor2->Y_IMAGE,
                   pSextractor2->maskSector);
          }
        }
#endif


#if 1
        if (pSextractorDim->sextractorIndex > pOverlapEntry->endSextractorIndex) {
          printf("ERROR: line %4d endSextractorIndex %d exceeded %d debugCount %d\n",__LINE__,pOverlapEntry->endSextractorIndex,pSextractorDim->sextractorIndex,debugCount);
          exit(-1);
        }
#endif

        if (pSextractorDim->NUMBER == NUMBER) {
          pSextractor2 = pSextractorDim;
          continue;
        }

        if ((sextractorIndex > minSextractorIndex) &&
            (pSextractorDim->Y_IMAGE < pSextractor2->Y_IMAGE)) {
          printf("ERROR: line %4d sextractor table is not sorted by Y_IMAGE %f %f debugCount %d  %s\n",
                 __LINE__,
                 pSextractorDim->Y_IMAGE,
                 pSextractor2->Y_IMAGE,
                 debugCount,
                 sextractor_name);
          exit(-1);
        }


        pSextractor2 = pSextractorDim;

        if (pSextractorDim->Y_IMAGE < minYRange) {
          minSextractorIndex = sextractorIndex;
          continue;
        }
        if (pSextractorDim->Y_IMAGE > maxYRange) {
          break;
        }
        if ((pSextractorDim->X_IMAGE < minXRange) ||
            (pSextractorDim->X_IMAGE > maxXRange)) {
          continue;
        }
        if (pSextractorDim->MAG_ISO < pSearchEntry1->MAG_ISO) {
          /* The search entry must be brighter */
          continue;
        }
        if ((pSextractorDim->FLUX_MAX/pSextractorDim->THRESHOLD) < analysis_threshold) {
          continue;
        }
    

#ifdef REJECT_GRAIN
        if (pSextractorDim->MAG_ISO > (pMaskSector->MAG_ISO_med - pMaskSector->MAG_ISO_rms)) {
#if 0
          printf("line %4d Rejecting Dim Object %6d MAG_ISO %10f FWHM_WORLD %d limit %10f\n",
                 __LINE__,
                 sextractorIndex,
                 pSextractorDim->MAG_ISO,
                 pSextractorDim->FWHM_IMAGE,
                 pMaskSector->MAG_ISO_med - pMaskSector->MAG_ISO_rms);
#endif
          continue;
        }
#endif /* REJECT_GRAIN */      



#ifdef SEARCHTABLE_ONLY /* Select dim candidates only from the search table */
        if (pSextractorDim->searchStar == 0) {
          continue;
        }
#endif /* SEARCHTABLE_ONLY */



        if (pSextractorDim->dimSelectedFlag == 0) {
          pMaskSector->dimSelectedCount++;
          pSextractorDim->dimSelectedFlag = 1;
        }



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
        correlation_table1[correlationIndex]++;
      
        if (correlation_table1[correlationIndex] > pMaskSector->maxBinCount0) {
          pMaskSector->maxBinCount0 = correlation_table1[correlationIndex];
        }
        debugCount++;
#if 0
        if (debugCount > (147554747-1)) {
          printf("line %4d sextractorIndex min %7d loop %7d Dim %7d End %7d Y_IMAGE %12.6f maskSector %2d ",
                 __LINE__,
                 minSextractorIndex,
                 sextractorIndex,
                 pSextractorDim->sextractorIndex,
                 pOverlapEntry->endSextractorIndex,
                 pSextractorDim->Y_IMAGE,
                 pSextractorDim->maskSector);
          if (pSextractor2 != NULL) {
            printf("sextractorIndex2 %7d Y_IMAGE %12.6f maskSector %2d\n",
                   pSextractor2->sextractorIndex,
                   pSextractor2->Y_IMAGE,
                   pSextractor2->maskSector);
          }
        }
#endif

      }
    }
    debugCount++;
#if 0
    if (debugCount > (147554747-1)) {

      printf("line %4d sextractorIndex min %7d loop %7d Dim %7d End %7d Y_IMAGE %12.6f maskSector %2d ",
             __LINE__,
             minSextractorIndex,
             sextractorIndex,
             pSextractorDim->sextractorIndex,
             pOverlapEntry->endSextractorIndex,
             pSextractorDim->Y_IMAGE,
             pSextractorDim->maskSector);
      if (pSextractor2 != NULL) {
        printf("sextractorIndex2 %7d Y_IMAGE %12.6f maskSector %2d\n",
               pSextractor2->sextractorIndex,
               pSextractor2->Y_IMAGE,
               pSextractor2->maskSector);
      }
    }
#endif

    vectorIndex = 0;
    for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
      for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
        correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
        if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
          printf("ERROR: correlationIndex %d exceeds limits (6) %d\n",correlationIndex,maxCorrelationIndex);
          exit(-1);
        }
#if 0
        if (doPlots) {
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
          fprintf(gmtHandle[5],"%d\t%d\t%d\t%d\n",maskSector,xIndex,yIndex,correlation_table1[correlationIndex]);
        }
#endif
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
    printf("line %4d search_close3 sector %2d overlap %2d initial count_med0 %f count_rms0 %f, count_limit0 %d maxBinSNR0 %f for %s\n",__LINE__,maskSector,pOverlapEntry->tileOverlapCount,pMaskSector->count_med0,pMaskSector->count_rms0,pMaskSector->count_limit0,pMaskSector->maxBinSNR0,fileroot);
    if (doPlots) {
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
        fprintf(gmtHandle[1],"\tANALYSIS_THRESHOLD");
        fprintf(gmtHandle[1],"\tPlate\n");

        fprintf(gmtHandle[1],"----------\t--------\t------\t-----------------\t---------------\t---------\t---------\t-----------\t-----------\t---------\t--------\t---------");
        fprintf(gmtHandle[1],"\t------------------");
        fprintf(gmtHandle[1],"\t-----\n");


      }
      fprintf(gmtHandle[1],"%d\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d",maskSector,pMaskSector->aveXY,pMaskSector->convolutionArea,pMaskSector->count_med0,pMaskSector->count_rms0,(1.0*pMaskSector->count_limit0),pMaskSector->maxBinCount0,pMaskSector->maxBinSNR0,selected,totalMasks);
      fprintf(gmtHandle[1],"\t%f",analysis_threshold);

      fprintf(gmtHandle[1],"\t%s\n",fileroot2);
    }

    /* Now expand the convolution radius to the adjacent pixels */

    memset(correlation_table2,0,(correlationTablePixels)*(sizeof(int)));
    memset(correlation_histogram2,0,(correlationTablePixels)*(sizeof(int)));
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
            correlation_table2[correlationIndex2] += correlation_table1[correlationIndex];
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
    pMaskSector->count_limit1 = pMaskSector->count_med1 + (COUNT_LIMIT_RMS_FACTOR*pMaskSector->count_rms1);


    printf("line %4d search_close3 sector %2d overlap %2d convolutionArea %d count_med1 %f count_rms1 %f, count_limit1 %f maxBinCount1 %d maxBinSNR1 %f for %s\n",__LINE__,maskSector,pOverlapEntry->tileOverlapCount,pMaskSector->convolutionArea,pMaskSector->count_med1,pMaskSector->count_rms1,(1.0*pMaskSector->count_limit1),pMaskSector->maxBinCount1,pMaskSector->maxBinSNR1,fileroot);

#ifdef COUNT_LIMIT1
    if (searchTableSize > 0) {      
      pMaskSector->count_limit1 = COUNT_LIMIT1;
    }
#endif /* COUNT_LIMIT1 */
#if 0
    if (doPlots) {
      if (selected == 0) {
        fprintf(gmtHandle[1],"%d\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d",maskSector,pMaskSector->aveX,pMaskSector->aveY,pMaskSector->convolutionArea,pMaskSector->count_med1,pMaskSector->count_rms1,(1.0*pMaskSector->count_limit1),pMaskSector->maxBinCount1,pMaskSector->maxBinSNR1,selected,totalMasks);
        fprintf(gmtHandle[1],"\t%f",analysis_threshold);
        fprintf(gmtHandle[1],"\t%s\n",fileroot2);


      }
    }
#endif
    if (doPlots) {
      for (yIndex = -MAX_PIXEL_RADIUS; yIndex <= MAX_PIXEL_RADIUS; yIndex++) {
        for (xIndex = -MAX_PIXEL_RADIUS; xIndex <= MAX_PIXEL_RADIUS; xIndex++) {
          correlationIndex = (xIndex+MAX_PIXEL_RADIUS) + ((yIndex+MAX_PIXEL_RADIUS) * (2*(MAX_PIXEL_RADIUS+1)));
          if ((correlationIndex < 0) || (correlationIndex >= maxCorrelationIndex)) {
            printf("ERROR: correlationIndex %d exceeds limits (10) %d\n",correlationIndex,maxCorrelationIndex);
            exit(-1);
          }
          if ((correlation_table2[correlationIndex]) > pMaskSector->count_limit1) {
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
              fprintf(gmtHandle[0],"maskSector\tix\tiy\tbin_count\tsequence\n");
              fprintf(gmtHandle[0],"----------\t--\t--\t---------\t--------\n");
            }
            fprintf(gmtHandle[0],"%d\t%d\t%d\t%d\t%d\n",maskSector,xIndex,yIndex,correlation_table1[correlationIndex],0);
            fprintf(gmtHandle[0],"%d\t%d\t%d\t%d\t%d\n",maskSector,xIndex,yIndex,correlation_table2[correlationIndex],1);
          }

        }
      }
    }
    /* Find the histogram of correlation_table2 */
#ifndef SKIP_POISSON
    factorial = 1.;
    for (correlationIndex = 0; correlationIndex < maxCorrelationIndex; correlationIndex++) {
      correlation_sum0 += correlation_table1[correlationIndex];
      correlation_sum1 += correlation_table2[correlationIndex];
      correlationCount = correlation_table1[correlationIndex];
      if (correlationCount >= maxCorrelationIndex) {
        correlationCount = maxCorrelationIndex-1;
      }
      correlation_histogram2[correlationCount]++;
    }
    meanArrival = (1.0*correlation_sum0)/(1.0*maxCorrelationIndex);
    exponential = exp(-meanArrival);
    overflowPrint = 0;
    if ((isnormal(meanArrival) == 0) || (isnormal(exponential)) == 0) {
      overflowFlag = 1;
    } else {
      overflowFlag = 0;
    }
    for (correlationIndex = 0; correlationIndex < maxCorrelationIndex; correlationIndex++) {
      if (correlationIndex != 0) {
        factorial=factorial*(1.0*correlationIndex);
        if ((overflowFlag == 0) && (isnormal(factorial) == 0)) {
          overflowFlag = 1;
        }
      }
      prefix = pow(meanArrival,1.0*correlationIndex);
      pmf2 = prefix*exponential/factorial;
      if ((overflowFlag == 0) & ((isnormal(prefix) == 0) || (isnormal(pmf2) == 0))) {
        overflowFlag = 1;
      }
      if (overflowFlag == 1) {
        pmf2 = 0;
      }
      bin_ratio = (1.0*correlation_histogram2[correlationIndex]);
      if (pmf2 == 0) {
        matchRatio = 0;
      } else {
        matchRatio = bin_ratio/pmf2;
      }
#if 1
      if ((correlation_histogram2[correlationIndex] > 0) || (pmf2 > 1.0)) {
        printf("maskSector %d correlationIndex %d prefix %e exponential %e, factorial %e, pmf2 %e matchRatio %e overflowFlag %d line %4d\n",maskSector,correlationIndex,prefix,exponential,factorial,pmf2,matchRatio,overflowFlag,__LINE__);
      }
      
#endif
      if (pmf2 > 1.0) {
        printf("maskSector %d correlationIndex %d prefix %e exponential %e, factorial %e, pmf2 %e matchRatio %e overflowFlag %d line %4d\n",maskSector,correlationIndex,prefix,exponential,factorial,pmf2,matchRatio,overflowFlag,__LINE__);
        printf("ERROR: pmf is %f in maskSector %d correlationIndex %d for %s\n",pmf2,maskSector,correlationIndex,fileroot);
        exit(-1);
      }
      pmf_sum0 += pmf2;
      bin_sum0 += bin_ratio;
      if (correlation_histogram2[correlationIndex] > 0) {
        if (doPlots) {
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
            fprintf(gmtHandle[7],"maskSector\tcorrelationCount\tstd2\tbin_sum\tmatchRatio\tmaxCorrelationIndex\n");
            fprintf(gmtHandle[7],"----------\t----------------\t----\t-------\t----------\t-------------------\n");
          }
          if (bin_ratio > 0) {
            fprintf(gmtHandle[7],"%d\t%d\t%e\t%f\t%e\t%d\n",maskSector,correlationIndex,pmf2,bin_ratio,matchRatio,maxCorrelationIndex);
          }
        }
        if (bin_ratio > 0) {
#ifndef DEBUG_NO_MYSQL
          sprintf(queryString,"INSERT INTO masks3 (series,plateNumber,mosaicNumber,nxs,nys,maskSector,correlationCount,std2,bin_sum,maxCorrelationIndex,countRMSFactor,count_med,count_rms,searchTableSize) values ('%s',%d,%d,%d,%d,%d,%d,%e,%f,%d,%d,%f,%f,%d);",
                  series,
                  plateNumber,
                  mosaicNumber,
                  nxs,
                  nys,
                  maskSector,
                  correlationIndex,
                  pmf2,
                  bin_ratio,
                  maxCorrelationIndex,
                  COUNT_LIMIT_RMS_FACTOR,
                  pMaskSector->count_med0,
                  pMaskSector->count_rms0,
                  searchTableSize);
          if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
            printf("ERROR: search_close3 MAX_QUERY_STRING exceeded\n");
            exit(-1);
          }
          res = ExecuteQuery(pConnection,queryString);
          if (res) {
            exit(-1);
          }
#endif /* DEBUG_NO_MYSQL */
        }

      }
    }

    printf("Sector %2d radius 0 total %d  radius 1 total %d meanArrival %f pmf_sum0 %f bin_sum0 %f line %4d\n",maskSector,correlation_sum0,correlation_sum1,meanArrival,pmf_sum0,bin_sum0,__LINE__);
#endif /* SKIP_POISSON */
    

    /* Now we need to identify the peak correlation masks and print a table of mask characteristics */
#ifdef DEBUG_TRUNCATION_X
    if (debugPrint == 0) {
      printf("ERROR: setting totalMasks to zero in line %4d\n",__LINE__);
      debugPrint = 1;
#if 0
      exit(-1);
#endif
    }
    totalMasks = 0;
#else /* DEBUG_TRUNCATION */
    
    if (searchTableSize > 0) {
      FindMasks(maskSector,pOverlapEntry->tileOverlapCount,pMaskSector->count_limit1,correlation_table2,correlation_table3,&maskTable,&totalMasks,fileroot);
    }

    /* Expand all entries in the mask table by the correlation radius and find the limits.  Set the centerFlag if necessary */
    pMaskSector->maxMaskArea = 0;
    pMaskSector->maxCenterMaskArea = 0;

    for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
      pMask = &maskTable[maskIndex];
      if (pMask->maskCenterFlag != 0) {
        if (pMaskSector->maxCenterMaskArea < pMask->maskArea) {
          pMaskSector->maxCenterMaskArea = pMask->maskArea;
        }
      } else {
        if (pMaskSector->maxMaskArea < pMask->maskArea) {
          pMaskSector->maxMaskArea = pMask->maskArea;
          
          /* The following is used for the getfits algorithm */
          if (pOverlapEntry->tileOverlapCount > cutoffTileOverlapCount) {
            if (pMaskSector->maxMaskArea >  maxDualMaskArea) {
              maxDualMaskArea = pMaskSector->maxMaskArea;
              maxDualMaskAreaMaskSector = maskSector;
            }
          } else {
            if (pMaskSector->maxMaskArea >  maxSingleMaskArea) {
              maxSingleMaskArea = pMaskSector->maxMaskArea;
              maxSingleMaskAreaMaskSector = maskSector;
            }

          }
        }
      }
    }


#endif /* DEBUG_TRUNCATION */
#ifdef DEBUG_DELETE_MASKS

    totalMasksSector += totalMasks;
    /* for the following grep "line 3188" los.log | awk '{OFS="\t"}{print $5,$7,$9}' */
    printf("line %4d search_close3 Sector %2d overlap %2d  totalMasksSector %6d totalMasks %6d sectorCount %6d dimSelectedCount %6d searchTableSize %4d curClipCount %6d MAG_ISO_med %f MAG_ISO_RMS %f maxMaskArea %4d maxCenterMaskArea %4d for %s\n",
           __LINE__,
           maskSector,
           pOverlapEntry->tileOverlapCount,
           totalMasks,
           totalMasksSector,
           pMaskSector->sectorCount,
           pMaskSector->dimSelectedCount,
           savedSearchTableSize,
           pMaskSector->curClipCount,
           pMaskSector->MAG_ISO_med,
           pMaskSector->MAG_ISO_rms,
           pMaskSector->maxMaskArea,
           pMaskSector->maxCenterMaskArea,
           fileroot);
    if (maskTable != NULL) {
      free(maskTable);
      maskTable = NULL;
    }
    totalMasks = 0;

#else /* DEBUG_DELETE_MASKS */
    printf("line %4d search_close3 Sector %2d overlap %2d  totalMasksSector %6d totalMasks %6d dimSelectedCount %6d searchTableSize %4d searchTableSize %4d curClipCount %6d MAG_ISO_med %f MAG_ISO_RMS %f for %s\n",
           __LINE__,
           maskSector,
           pOverlapEntry->tileOverlapCount,
           totalMasksSector,
           totalMasks,
           pMaskSector->dimSelectedCount,
           searchTableSize,
           pMaskSector->curClipCount,
           pMaskSector->MAG_ISO_med,
           pMaskSector->MAG_ISO_rms,
           fileroot);
    totalMasksSector = 0;
#endif /* DEBUG_DELETE_MASKS  */

    if (totalMasks >= SEARCH_TABLE_LENGTH) {
      printf("ERROR: line %4d search_close3 totalMasks %d greater than SEARCH_TABLE_LENGTH %d for %s\n",__LINE__,totalMasks,SEARCH_TABLE_LENGTH,fileroot);
    } else {

      if (doPlots) {
        fprintf(gmtHandle[1],"%d\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d",maskSector,pMaskSector->aveXY,pMaskSector->convolutionArea,pMaskSector->count_med1,pMaskSector->count_rms1,(1.0*pMaskSector->count_limit1),pMaskSector->maxBinCount1,pMaskSector->maxBinSNR1,selected,totalMasks);
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

          for (sextractorIndex = pOverlapEntry->startSextractorIndex; sextractorIndex < pOverlapEntry->endSextractorIndex; sextractorIndex++) {
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
                  printf("ERROR: search_close3 maskIndex %d %d greater than %d for %s\n",pMaskOrder->maskIndex,pMaskOrder2->maskIndex,SEARCH_TABLE_LENGTH,fileroot);
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
#ifndef DEBUG_NO_MYSQL
        SetMosaicFitWCS(pConnection,"HaveMultipleMask",fileroot,0,1,readFlag,&FitWCS);
        if ((totalMasks <= MAX_MASKS) &&
            (maskCenterIndex == 0)) {
          SetPlateQuality(pConnection,"multiple",fileroot,1,0,&quality);
        }

#endif /* DEBUG_NO_MYSQL */
      } else {
#ifndef DEBUG_NO_MYSQL
        SetMosaicFitWCS(pConnection,"NoMultipleMask",fileroot,0,1,readFlag,&FitWCS);
#endif /* DEBUG_NO_MYSQL */
      }


      for (maskIndex  = 0; maskIndex < totalMasks; maskIndex++) {
        pMask = &maskTable[maskIndex];
        pMask->maskIndex = maskIndex+1;
#ifndef DEBUG_TRUNCATION
        /* Load the masks into the datbase, but only 1 mask if there is an unusual condition */
        if (((maskCenterIndex == 0) && (totalMasks <= MAX_MASKS)) ||
            ((maskCenterIndex != 0) && (pMask->maskCenterFlag != 0)) ||
            ((totalMasks > MAX_MASKS) && (maskIndex == 0))) {
#ifndef DEBUG_NO_MYSQL
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
            printf("ERROR: search_close3 MAX_QUERY_STRING exceeded\n");
            exit(-1);
          }
#if 0
          printf("%s\n",queryString);
#endif

          res = ExecuteQuery(pConnection,queryString);
          if (res) {
            exit(-1);
          }
#endif /* DEBUG_NO_MYSQL */

        }
#endif /* DEBUG_TRUNCATION */

        if (doPlots) {
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
            minXRange = pMask->maskXMin;
            maxXRange = pMask->maskXMax;
            minYRange = pMask->maskYMin;
            maxYRange = pMask->maskYMax;
          } else {
            if (minXRange > pMask->maskXMin) {
              minXRange = pMask->maskXMin;
            }
            if (maxXRange < pMask->maskXMax) {
              maxXRange = pMask->maskXMax;
            }
            if (minYRange > pMask->maskYMin) {
              minYRange = pMask->maskYMin;
            }
            if (maxYRange < pMask->maskYMax) {
              maxYRange = pMask->maskYMax;
            }
          }
        }
        /* Now for every sextractor object pair, look for a match in the mask table */
        minSextractorIndex = pOverlapEntry->startSextractorIndex;
        for (sextractorIndex1 = pOverlapEntry->startSextractorIndex; sextractorIndex1 < pOverlapEntry->endSextractorIndex; sextractorIndex1++) {
          pSextractor1 = &image_table[sextractorIndex1];
          if ((pSextractor1->FLUX_MAX/pSextractor1->THRESHOLD) < analysis_threshold) {
            continue;
          }
    



          for (sextractorIndex2 = minSextractorIndex; sextractorIndex2 < pOverlapEntry->endSextractorIndex; sextractorIndex2++) {  
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


            if (yIndex < minYRange) {
              minSextractorIndex = sextractorIndex2;
              continue;
            }
            if (yIndex > maxYRange) {
              break;
            }
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
                (xIndex >= minXRange) && 
                (xIndex <= maxXRange)) {
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
    } /* End of totalMasks processing */
  } /* End of maskSector loop */
  pMaskSector = NULL;
  
  if (maxDualMaskAreaMaskSector >= 0) {
    highestDualMaskMaskSector = maxDualMaskAreaMaskSector;
  }
  if (highestDualMaskMaskSector < 0) {
    highestDualMaskMaskSector = maxSingleMaskAreaMaskSector;
  }

  if (highestDualMaskMaskSector < 0) {
    printf("ERROR: line %4d bad highestDualMaskMaskSector for %s\n",__LINE__,fileroot);
  } else { 
    int centerXPixel;
    int startXPixel;
    int endXPixel;
    int centerYPixel;
    int startYPixel;
    int endYPixel;
    struct stat statbuf;
    int statResult;

    pMaskSector = &maskSectorTable[highestDualMaskMaskSector];
    pOverlapEntry = &pOverlapEntryTable[highestDualMaskMaskSector]; 
    bestInspectionTileIndex = -1;
    bestInspectionTileLower = -1;
    bestInspectionTileCount  = 0;
    if (pOverlapZone->maxInspectionTiles < 2) {
      printf("ERROR: line %4d bad maxInspectionTiles %de for %s\n",__LINE__,pOverlapZone->maxInspectionTiles,fileroot);
    }

    for(inspectionTileIndex = 1; inspectionTileIndex < (pOverlapZone->maxInspectionTiles-2); inspectionTileIndex++) {
      if (pOverlapEntry->inspectionTileTableLower[inspectionTileIndex] > bestInspectionTileCount) {
        bestInspectionTileCount = pOverlapEntry->inspectionTileTableLower[inspectionTileIndex];
        bestInspectionTileLower = 1;
        bestInspectionTileIndex = inspectionTileIndex;
      }
      if (pOverlapEntry->inspectionTileTableUpper[inspectionTileIndex] > bestInspectionTileCount) {
        bestInspectionTileCount = pOverlapEntry->inspectionTileTableLower[inspectionTileIndex];
        bestInspectionTileLower = 0;
        bestInspectionTileIndex = inspectionTileIndex;
      }
    }
    
    if (strstr(fileroot,"ww") == NULL) {
      sprintf(mosaicFileName,"/dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s.fit",
              pMosaic->diskLocation,
              pMosaic->series,
              pMosaic->plateNumber,
              pMosaic->mosaicNumber,
              fileroot);
      statResult = Stat(mosaicFileName,&statbuf);
    } else {


      sprintf(mosaicFileName,"/dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s_tnx.fit",
              pMosaic->diskLocation,
              pMosaic->series,
              pMosaic->plateNumber,
              pMosaic->mosaicNumber,
              fileroot);
      statResult = Stat(mosaicFileName,&statbuf);
      if (statResult != 0) {
        sprintf(mosaicFileName,"/dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s.fit",
                pMosaic->diskLocation,
                pMosaic->series,
                pMosaic->plateNumber,
                pMosaic->mosaicNumber,
                fileroot);
        statResult = Stat(mosaicFileName,&statbuf);
      }
    }
    if (statResult != 0) {
      printf("line %4d ERROR: could not find mosaic %s\n",__LINE__,mosaicFileName);
    }

    if (pOverlapZone->xPosXImageFlag == 1) {
      /* need to rotate the result 90 degrees */
      if (bestInspectionTileLower == 1) {
        centerYPixel = pOverlapEntry->startPixel;
      } else {
        centerYPixel =  pOverlapEntry->endPixel;
      }
      startYPixel = centerYPixel - (INSPECTION_WIDTH/2);
      if (startYPixel <= 0) {
        startYPixel = 1;
      }
      endYPixel = centerYPixel + (INSPECTION_WIDTH/2);
      if (endYPixel > mosaicHeight) {
        endYPixel = mosaicHeight;
      }
      centerXPixel = ((1.0*bestInspectionTileIndex)+0.5)*INSPECTION_HEIGHT;
      startXPixel = centerXPixel-INSPECTION_HEIGHT/2;
      if (startXPixel <= 0) {
        startXPixel = 1;
      }
      endXPixel = centerXPixel+INSPECTION_HEIGHT/2;
      if (endXPixel > mosaicWidth) {
        endXPixel = mosaicWidth;
      }
      printf("line %2d getfits -o %s%s_search_close3.tmp %s %d-%d %d-%d\n",
             __LINE__,
             getfitsdirectory,
             fileroot,
             mosaicFileName,
             startXPixel,
             endXPixel,
             startYPixel,
             endYPixel);
      printf("line %2d imrot -r 270 -o %s%s_search_close3.fit %s%s_search_close3.tmp \n",
             __LINE__,
             getfitsdirectory,
             fileroot,
             getfitsdirectory,
             fileroot);



    
    } else {
      if (bestInspectionTileLower == 1) {
        centerXPixel = pOverlapEntry->startPixel;
      } else {
        centerXPixel =  pOverlapEntry->endPixel;
      }
      startXPixel = centerXPixel - (INSPECTION_WIDTH/2);
      if (startXPixel <= 0) {
        startXPixel = 1;
      }
      endXPixel = centerXPixel + (INSPECTION_WIDTH/2);
      if (endXPixel > mosaicWidth) {
        endXPixel = mosaicWidth;
      }
      centerYPixel = ((1.0*bestInspectionTileIndex)+0.5)*INSPECTION_HEIGHT;
      startYPixel = centerYPixel-INSPECTION_HEIGHT/2;
      if (startYPixel <= 0) {
        startYPixel = 1;
      }
      endYPixel = centerYPixel+INSPECTION_HEIGHT/2;
      if (endYPixel > mosaicHeight) {
        endYPixel = mosaicHeight;
      }
      printf("line %2d getfits -o %s%s_search_close3.fit %s %d-%d %d-%d\n",
             __LINE__,
             getfitsdirectory,
             fileroot,
             mosaicFileName,
             startXPixel,
             endXPixel,
             startYPixel,
             endYPixel);

    }
    
  }
  pMaskSector = 0;
  pOverlapEntry = 0;
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


      if (doPlots) {
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
      if (doPlots) {
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
#ifndef DEBUG_NO_MYSQL
  /* Update the quality bits for this plate */
  UpdateQuality(pConnection,series,plateNumber);
#endif /* DEBUG_NO_MYSQL */
  
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
  if (maskSectorTable != NULL) {
    free(maskSectorTable);
  }

  if (pBinMap != NULL) {
    free(pBinMap);
  }

  mysql_close(pConnection);

  time(&curTime);
  curTime -= startTime;
#ifdef DEBUG_DELETE_MASKS
  if (totalMasksSector != 0) {
    totalMasks = totalMasksSector;
  }
#endif /* DEBUG_DELETE_MASKS  */


  printf("Candidates: %d plotted %d psfsaturated %d populatedBinCount %d maxBinCountOverall %d outputCount %d flip %d mirror %d rotation %3d for %s  line %4d\n",
         candidateCount,
         plotCount,
         psfsaturatedCount,
         populatedBinCount,
         maxBinCountOverall,
         outputCount,
         pOverlapZone->flipFlag,
         pOverlapZone->mirrorFlag,
         pOverlapZone->rotation,

         fileroot,
         __LINE__);

  printf(" threshold %f totalMasks %d patternID %2d seconds %d for %s\n",analysis_threshold,totalMasks,pMosaic->patternID,curTime,fileroot);

  return(0);
}

