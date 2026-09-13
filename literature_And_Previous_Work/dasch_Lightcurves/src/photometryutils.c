// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* Includes code from statstable.c,
 *                     Starbase Data Tables - An Ascii Database for UNIX
 *                     copyright 1991, 1993, 1995, 1999 John B. Roll jr.
 *
 * "Bevington" below refers to Philip R. Bevington "Data Reduction and Error Analysis for the Physical Sciences" McGraw-Hill 1969
 *
 * The Sky2KCompare,GetSky2kCatalog and SearchGSCBlends routines filter unmatched objects next to bright
 * GSC2.3.2 catalog stars using the Sky2K catalog.
 * The algorithm is based on the following measurements from plate bm00631 which includes Vega with one point
 * from plate rh03313 which includes gamma Cygni.  Below magnitude 12, the catalog appears to be complete.
 *
 *     Magnitude   Radius in     ln(Radius)
 *                  degrees
 *      y                              x
 *     0.55         0.324           -1.127
 *     2.65         0.0923          -2.383
 *     4.39         0.0504          -2.988
 *     5.60         0.070           -2.659
 *     5.7          0.065           -2.733
 *     6.78         0.038           -3.270
 *     8.0          0.015           -4.200
 *     9.63         0.008           -4.828
 *    10.45         0.0097          -4.636
 *    11.47         0.0044          -5.426
 *    12.54         0.0043          -5.449
 *
 *
 *     n = 11 points
 *     sum of x -39.699
 *     sum of x**2 162.685
 *     sum of y  77.76
 *     sum of y**2 691.17
 *     sum of xy -331.797
 *
 *     intercept = -2.44326
 *     slope = -2.6357
 *     mag = -2.44326 - ((2.6357)* ln(Radius))
 *
 *     if x = 0, y = -2.4432
 *     if y = 0  x = -0.9269
 *
 *     ln(Radius) = -0.9269 - ((0.3793)* (mag))
 */

#define PROPER_MOTION_EPOCH_SPAN 50
#define MIN_PROPER_MOTION_COUNT 30
#define MIN_FWHM_WORLD (3.0 / 3600.0) /* 3 arcsec */

#define _GNU_SOURCE /* for exp10() */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mysql.h"
#include <sys/types.h>
#include <grp.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>

#include <fcntl.h>
#include <limits.h>
#include <assert.h>
#include "pipelineutils.h"
#include "photometryutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "libwcs/wcscat.h"

#include "kdtree.h"
void *ptree = NULL;

#define DESIRED_UMASK 2 /* Desired value of umask */

#define MAX_FILENAME 512
#define MAX_PLOT_TEXT 32
#define MAX_GSC_FILE_BUFFER 512
#define GSC_NEIGHBORS 3

extern int parse_coo(char *coi, double *a, double *d);

extern GSCBIN gscBin01;
extern SERIESENTRY seriesTableX[MAX_SERIES];
extern int seriesTableXInit;

static int maxTreeIndex = 0;
static int allocTreeCount = 0;
static PSERIESTREE seriesTreeTable = NULL;

char *colortable[COLORFLAG_MAX] = {"", "linear", "metropolis", "copy"};
int damonSeriesId[NUM_DAMON_SERIES]; /* Series id for the damon series */
int excludeSeriesId[NUM_EXCLUDE_SERIES]; /* Series to exclude */

static PGSCIMAGE GetCatalogRecord(PFILECOMMON pFileCommon, long long REFNumber, int gsc_bin_index);

int
ProcessLongTransientCandidates(
  PFILECOMMON pFileCommon,
  PTRANCOMMON pTranCommon,
  PPARAMETERVECTORSTORE pVectorStore,
  PPHOTSTARIMAGE pMagnitudeTable,
  PGALAXYCOMMON pGalaxyCommon,
  PSTARENTRY pCurStarEntry,
  char* nearbyObjects,
  int AFLAGSMASK1,
  int QUALITYMASK
);


/* The return value comes from getenv(), and so should not be freed. */
char *
GetPhotFileBase(char *catalogString)
{
  char dasch_phot_magnitudes[40];
  char *photfilebase;
  int n;

  n = snprintf(dasch_phot_magnitudes, sizeof(dasch_phot_magnitudes), "DASCH_PHOT_MAGNITUDES%s", catalogString);
  photfilebase = getenv(dasch_phot_magnitudes);

  if (n < 0 || photfilebase == NULL) {
    printf("$DASCH_PHOT_MAGNITUDES%s is not defined\n", catalogString);
    exit(1);
  }

  return photfilebase;
}


void
InitSeriesTable(MYSQL *pConnection, MYSQL *pPhotConnection) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  int seriesId;
  int numSeries = 0;
  int maxSeriesId = 0;
  int numSeries2 = 0;
  int maxSeriesId2 = 0;
  int excludeSeriesIndex;
  int excludeSeriesCount = NUM_EXCLUDE_SERIES;
  unsigned char * charPtr;
  PSERIESTREE pCurTree;
  PSERIESTREE pNextTree;
  PSERIESTREE tmpSeriesTreeTable;
  int charCount;
  int stringLength;
  char charVal;
  int charIndex;

  if (seriesTableXInit != 0) {
    return;
  }

  memset(seriesTableX,0,sizeof(seriesTableX));
  /* First read everything in from the photseries table */
  sprintf(queryString,"SELECT seriesId,series from photseries;");
  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&seriesId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for seriesId in InitSeriesTable\n", nvals);
            exit(1);
          }
        } else {
          printf("ERROR: NULL seriesId in InitSeriesTable\n");
          exit(1);
        }
        if ((seriesId <= 0) || (seriesId >= MAX_SERIES)) {
          printf("ERROR: illegal seriesId %d in InitSeriesTable\n", seriesId);
          exit(1);
        }
        if (seriesId > maxSeriesId) {
          maxSeriesId = seriesId;
        }
        if (sqlrow[1]) {
          if (strlen(sqlrow[1]) > (MAX_SERIES_STRING-1)) {
            printf("ERROR: series %s has length %zu in InitSeriesTable\n",sqlrow[1],strlen(sqlrow[1]));
            exit(1);
          }
          strcpy(seriesTableX[seriesId].series,sqlrow[1]);
        } else {
          printf("ERROR: NULL series in InitSeriesTable\n");
          exit(1);
        }
        numSeries++;

      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }
  if (numSeries != maxSeriesId) {
    printf("ERROR: series table entries %d does not agree with maximum seriesId %d in InitSeriesTable\n",numSeries,maxSeriesId);
    exit(1);
  }
  if (pConnection != NULL) {
    /* Now read everything in from the series table */
    sprintf(queryString,"SELECT seriesId,series,nominalPlateScale,fittedPlateScale,sequestered+0 from series;");
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&seriesId);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for seriesId in InitSeriesTable\n", nvals);
              exit(1);
            }
          } else {
            printf("ERROR: NULL seriesId in InitSeriesTable\n");
            exit(1);
          }
          if ((seriesId <= 0) || (seriesId >= MAX_SERIES)) {
            printf("ERROR: illegal seriesId %d in InitSeriesTable\n", seriesId);
            exit(1);
          }
          if (seriesId > maxSeriesId2) {
            maxSeriesId2 = seriesId;
          }
          if (sqlrow[1]) {
            if (strlen(sqlrow[1]) > (MAX_SERIES_STRING-1)) {
              printf("ERROR: series %s has length %zu in InitSeriesTable\n",sqlrow[1],strlen(sqlrow[1]));
              exit(1);
            }
            if (strcmp(seriesTableX[seriesId].series,sqlrow[1]) != 0) {
              printf("ERROR: series entry %d %s and %s do not agree between series and photseries tables\n",seriesId,seriesTableX[seriesId].series,sqlrow[1]);
              exit(1);
            }
          } else {
            printf("ERROR: NULL series in InitSeriesTable\n");
            exit(1);
          }
          if (sqlrow[2]) {
            nvals = sscanf(sqlrow[2],"%lf",&seriesTableX[seriesId].nominalPlateScale);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for nominalPlateScale series %s in InitSeriesTable\n",nvals,seriesTableX[seriesId].series);
              exit(1);
            }
          }
          if (sqlrow[3]) {
            nvals = sscanf(sqlrow[3],"%lf",&seriesTableX[seriesId].fittedPlateScale);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for nominalPlateScale series %s in InitSeriesTable\n",nvals,seriesTableX[seriesId].series);
              exit(1);
            }
          }
          if (sqlrow[4]) {
            nvals = sscanf(sqlrow[4],"%d",&seriesTableX[seriesId].sequestered);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for sequestered series %s in InitSeriesTable\n",nvals,seriesTableX[seriesId].series);
              exit(1);
            }
            if ((seriesTableX[seriesId].sequestered != SEQUESTERED_NO) &&
                (seriesTableX[seriesId].sequestered != SEQUESTERED_YES)) {
              printf("ERROR: sequestered is %d in series %s in InitSeriesTable\n",seriesTableX[seriesId].sequestered,seriesTableX[seriesId].series);
              exit(1);
            }

          }

          numSeries2++;

        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
      }
    }
    if (numSeries2 != maxSeriesId2) {
      printf("ERROR: series table entries %d does not agree with maximum seriesId %d in InitSeriesTable\n",numSeries2,maxSeriesId2);
      exit(1);
    }
    if (maxSeriesId != maxSeriesId2) {
      printf("ERROR: series table entries %d %d do not agree in series and photseries table InitSeriesTable\n",maxSeriesId,maxSeriesId2);
      exit(1);
    }

  }
  if (pConnection == NULL) {
    seriesTableXInit = 1;
  } else {
    seriesTableXInit = 2;

  }
  allocTreeCount = SERIES_ALLOC_INCREMENT;
  seriesTreeTable = (PSERIESTREE)calloc(allocTreeCount,sizeof(SERIESTREE));
  if (seriesTreeTable == NULL) {
    printf("ERROR: failed to allocate the seriesTreeTable line %d\n",__LINE__);
    exit(1);
  }
  pCurTree = &seriesTreeTable[0];
  memset(pCurTree,0,sizeof(SERIESTREE));
  maxTreeIndex = 1;

  for (seriesId = 0; seriesId < MAX_SERIES; seriesId++) {
    if (seriesTableX[seriesId].series[0] == 0) {
      continue;
    }
    if (seriesId == 0) {
      printf("ERROR: Illegal seriesId of zero line %d\n",__LINE__);
      exit(1);
    }

    charCount = 0;
    charPtr = (unsigned char *) seriesTableX[seriesId].series;
    stringLength = strlen((char *) charPtr);
    charVal = *charPtr;
    while (1) {
      if (charVal == 0) {
        printf("ERROR: null character for series %s line %d\n",seriesTableX[seriesId].series,__LINE__);
      }
      charCount++;
      charIndex = charVal - 'a';
      if ((charIndex < 0) || (charIndex >= MAX_SERIES_CHAR)) {
        printf("ERROR: illegal series character  '%c; in %s line %d\n",charVal,seriesTableX[seriesId].series,__LINE__);
        exit(1);
      }
      if (charCount == 1) {
        pCurTree = &seriesTreeTable[0];
      }
      if (pCurTree->charIndex[charIndex] == 0) {
        /* Need to allocate another level */
        if (maxTreeIndex >= allocTreeCount) {
          /* Need to allocate another block */
          allocTreeCount += SERIES_ALLOC_INCREMENT;
          tmpSeriesTreeTable = (PSERIESTREE)realloc(seriesTreeTable,allocTreeCount*sizeof(SERIESTREE));
          if (tmpSeriesTreeTable == NULL) {
            printf("ERROR: failed to reallocate the seriesTreeTable of size %d line %d\n",allocTreeCount,__LINE__);
            exit(1);
          }
          seriesTreeTable = tmpSeriesTreeTable;
          tmpSeriesTreeTable = NULL;
        }
        pNextTree = &seriesTreeTable[maxTreeIndex];
        memset(pNextTree,0,sizeof(SERIESTREE));
        pNextTree->level=charCount;
        pNextTree->charVal = charVal;
        pNextTree->treeIndex = maxTreeIndex;
        pCurTree->charIndex[charIndex] = maxTreeIndex;
        maxTreeIndex++;
      } else {
        pNextTree = &seriesTreeTable[pCurTree->charIndex[charIndex]];
      }
#if 0
      printf("line %d level %d char %c\n",__LINE__,pNextTree->level,pNextTree->charVal);
#endif
      if (charCount == stringLength) {
        if (pNextTree->seriesId != 0) {
          printf("ERROR: duplicate seriesId %d for %s line %d\n",pNextTree->seriesId,seriesTableX[seriesId].series,__LINE__);
          exit(1);
        }
#if 0
        printf("line %d level %d seriesId %d series %s\n",__LINE__,pNextTree->level,seriesId,seriesTableX[seriesId].series);
#endif
        pNextTree->seriesId = seriesId;
        break;
      }
      /* go to the next level */
      pCurTree = pNextTree;
      charPtr++;
      charVal = *charPtr;
    }
  }
#if 0
    printf("line %4d maxTreeIndex is %d for allocTreeCount %d\n",__LINE__,maxTreeIndex,allocTreeCount);
#endif
  for (seriesId = 0; seriesId < MAX_SERIES; seriesId++) {
    if (seriesTableX[seriesId].series[0] == 0) {
      continue;
    }
    if (GetSeriesId(seriesTableX[seriesId].series,1) != seriesId) {
      printf("ERROR: incorrect series hash table\n");
      exit(1);
    }
  }

  damonSeriesId[0] = GetSeriesId("dnb",1);
  damonSeriesId[1] = GetSeriesId("dsb",1);
  damonSeriesId[2] = GetSeriesId("dnr",1);
  damonSeriesId[3] = GetSeriesId("dsr",1);
  damonSeriesId[4] = GetSeriesId("dny",1);
  damonSeriesId[5] = GetSeriesId("dsy",1);

  for (excludeSeriesIndex = 0; excludeSeriesIndex < excludeSeriesCount; excludeSeriesIndex++) {
    excludeSeriesId[excludeSeriesIndex] = GetSeriesId(excludeSeriesList[excludeSeriesIndex],1);
  }

  return;
}


int
GetSeriesId(char *series, int fatal) {
  int seriesId;
  PSERIESTREE pCurTree;
  char *charPtr;
  char charVal;
  int charIndex;
  int nextTreeIndex;

  if (seriesTableXInit == 0) {
    printf("ERROR: Series table not initialized with a call to InitSeriesTable\n");
    exit(1);
  }

  if (fatal) {
    pCurTree = &seriesTreeTable[0];
    charPtr = series;
    while (1) {
      charVal = *charPtr;
      if (charVal == 0) {
        if ((pCurTree->seriesId <= 0) && (pCurTree->seriesId >= MAX_SERIES)) {
          printf("ERROR: illegal seriesId %d; in %s photometryutils line %d\n",pCurTree->seriesId,series,__LINE__);
          exit(1);
        }
        return(pCurTree->seriesId);
      }
      charIndex = charVal - 'a';
      if ((charIndex < 0) || (charIndex >= MAX_SERIES_CHAR)) {
        printf("ERROR: illegal series character  '%c; in %s photometryutils line %d\n",charVal,series,__LINE__);
        exit(1);
      }
      nextTreeIndex =  pCurTree->charIndex[charIndex];
      if ((nextTreeIndex <= 0) || (nextTreeIndex > maxTreeIndex)) {
        printf("ERROR: illegal tree index %d; in %s photometryutils line %d\n",nextTreeIndex,series,__LINE__);
        exit(1);
      }
      pCurTree = &seriesTreeTable[nextTreeIndex];
      charPtr++;
    }


  } else {
    if (seriesTableXInit == 0) {
      printf("ERROR: Series table not initialized with a call to InitSeriesTable\n");
      exit(1);
    }
    for (seriesId = 0; seriesId < MAX_SERIES; seriesId++) {
      if (strcmp(series,seriesTableX[seriesId].series) == 0) {
        return(seriesId);
      }

    }
    if (fatal) {
      printf("ERROR: series %s is unrecognized in GetSeriesId\n",series);
      exit(1);
    } else {
      return(-1);
    }
  }
  return(-1);
}


/* Returns the plate scale in degrees/pixel */
double
GetFittedPlateScale(int seriesId, int plateNumber) {
  double plateScale;
  char *series;

  if (seriesTableXInit != 2) {
    printf("ERROR: Series table not initialized for plate scale with a call to InitSeriesTable\n");
    exit(1);
  }
  if ((seriesId <= 0) || (seriesId >= MAX_SERIES)) {
    printf("ERROR: Illegal series id %d in GetFittedPlateScale\n",seriesId);
    exit(1);
  }
  series = seriesTableX[seriesId].series;
  plateScale = seriesTableX[seriesId].fittedPlateScale;
  if (plateScale == 0.0) {
    plateScale = seriesTableX[seriesId].nominalPlateScale;
  }
  /* Now apply special knowledge that does not fit in the table (ugly hack) */
  /* This hack is repeated in the Catalog.cpp GetFittedScale routine */
  /* See /dasch/mysql/getscale.c */

  if (strcmp(series,"mc") == 0) {
    if  ((plateNumber < 4166) &&
         (plateNumber >= 3500)) {
      /* strcat(series,"2"); */
      plateScale = 97.96;
    } else if (plateNumber < 3499) {
      /* strcat(series,"1"); */
      plateScale = 93.21;
    } else {
      /* strcat(series,"3"); */ /* Vacuum "subsequent to 4171 */
      plateScale =  97.87;
    }
  }
  if (strcmp(series,"ac") == 0) {
    if ((plateNumber >=   1641) && (plateNumber <=   9169)) {
      /* strcat(series,"2"); */ /* Cooke #4665 */
      plateScale = 622.66;
    } else if (((plateNumber >=  9170) && (plateNumber <=  9207)) ||
               ((plateNumber >= 12361) && (plateNumber <= 12365))) {
      /* strcat(series,"3"); */ /* Cooke #832 */
      /* plateScale = 602.93; */  /* Only one entry! */
    } else if (((plateNumber >=  9208) && (plateNumber <= 12360)) ||
               ((plateNumber >= 12366) && (plateNumber <= 13050)) ||
               ((plateNumber >= 13051) && (plateNumber <= 13066))) {
      /* strcat(series,"4"); */ /* Cooke #17486 */
      plateScale = 599.67;
    } else if (plateNumber >= 13067) {
      /* strcat(series,"5"); */ /* Cooke #V30401 (or #17486 ??) */
      plateScale = 606.43;
    } else {
      /* strcat(series,"1"); */ /* Cooke no number */
      plateScale = 615.78;
    }

  }
  if (strcmp(series,"am") == 0) {
    if (plateNumber <= 16722) {
      /* strcat(series,"1"); */ /* 1 inch lens */
      plateScale = 610.66;
    } else {
      /* strcat(series,"2"); */ /* 1.5 inch lens */
      plateScale = 610.83;
    }
  }
  if (strcmp(series,"mb") == 0) {
    if (plateNumber <= 327) {
      /* strcat(series,"1"); */ /* 4 in Cooke lens  191"/mm */
      /* plateScale = 0;  */ /* No entries */
    } else if ((plateNumber >= 328) && (plateNumber <= 1776)) {
      /* strcat(series,"2"); */ /* 6 in             354"/mm */
      /* plateScale = 0;  */ /* No entries */
    } else {
      /* strcat(series,"3"); */ /* 3 in Ross Lundin 390"/mm */
      plateScale = 389.97;
    }
  }

  /* Covert from arcsec/mm to degrees/pixel */
  plateScale = (NOMINAL_MM_PER_PIXEL * plateScale) / 3600.0;
  return(plateScale);
}


int
GetPhotPlate(
  MYSQL *pPhotConnection,
  char *series,
  int plateNumber,
  PPHOTPLATES pPhotPlates,
  char *catalogString,
  int readOnly
) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  int gotAnswer = 0;
  int seriesId;
  memset(pPhotPlates,0,sizeof(PHOTPLATES));
  strcpy(pPhotPlates->series,series);
  pPhotPlates->plateNumber = plateNumber;

  seriesId = GetSeriesId(series,1);
  pPhotPlates->seriesId = seriesId;

#ifndef READONLY_PHOTOMETRY
  if (readOnly == 0) {
    sprintf(queryString,"INSERT IGNORE into photplates (seriesId,plateNumber) values (%d,%d);\n",
            seriesId,plateNumber);
    res = ExecuteQuery(pPhotConnection,queryString);
    if (res) {
      return(1);
    }
  }
#endif /* READONLY_PHOTOMETRY */

  sprintf(queryString,"SELECT versionId%s,THRESHOLD,scale,mosaicWidth,mosaicHeight,exposures,mosaicNumber,leftMargin,rightMargin,bottomMargin,topMargin,quality+0,unmatched from photplates where seriesId = %d and plateNumber = %d;\n",
          catalogString,seriesId,plateNumber);
  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (gotAnswer) {
          printf("ERROR: two photplates entries\n");
          memset(pPhotPlates,0,sizeof(PHOTPLATES));
          strcpy(pPhotPlates->series,series);
          pPhotPlates->plateNumber = plateNumber;
        }
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pPhotPlates->versionId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for versionId\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->versionId = 0;
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%lf",&pPhotPlates->THRESHOLD);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for THRESHOLD\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->THRESHOLD = 0.0;
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%lf",&pPhotPlates->scale);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for scale\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->scale = 0.0;
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pPhotPlates->mosaicWidth);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for mosaicWidth\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->mosaicWidth = 0;
        }
        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&pPhotPlates->mosaicHeight);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for mosaicHeight\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->mosaicHeight = 0;
        }
        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%d",&pPhotPlates->exposures);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for exposures\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->exposures = 0;
        }
        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%d",&pPhotPlates->mosaicNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for mosaicNumber\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->mosaicNumber = 0;
        }
        if (sqlrow[7]) {
          nvals = sscanf(sqlrow[7],"%d",&pPhotPlates->leftMargin);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for leftMargin\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->leftMargin = 0;
        }
        if (sqlrow[8]) {
          nvals = sscanf(sqlrow[8],"%d",&pPhotPlates->rightMargin);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for rightMargin\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->rightMargin = 0;
        }
        if (sqlrow[9]) {
          nvals = sscanf(sqlrow[9],"%d",&pPhotPlates->bottomMargin);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for bottomMargin\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->bottomMargin = 0;
        }
        if (sqlrow[10]) {
          nvals = sscanf(sqlrow[10],"%d",&pPhotPlates->topMargin);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for topMargin\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->topMargin = 0;
        }
        if (sqlrow[11]) {
          nvals = sscanf(sqlrow[11],"%d",&pPhotPlates->quality);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for quality\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->quality = 0;
        }
        if (sqlrow[12]) {
          nvals = sscanf(sqlrow[12],"%d",&pPhotPlates->unmatched);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for unmatched\n", nvals);
            continue;
          }
        } else {
          pPhotPlates->unmatched = 0;
        }
        gotAnswer++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }

  if (gotAnswer == 1) {
    return(0);
  } else {
    return(1);
  }
}


int
GetPhotometryGlobal(MYSQL *pPhotConnection, PPHOTGLOBAL pPhotGlobal)
{
  int nvals;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int gotAnswer = 0;
  char queryString[MAX_QUERY_STRING];
  char tempString[MAX_QUERY_STRING];
  int index;

  memset(pPhotGlobal,0,sizeof(PHOTGLOBAL));

  strcpy(queryString,"SELECT epoch,equinox,maxSpatialBins,keepalive+0,currentVersion,timeStamp,magnitudeFile+0,solutionNumber+0,minVersionId");
  for (index = 1; index < MAX_CATALOG_NUMBER; index++) {
    sprintf(tempString,",minVersionId%d",index);
    strcat(queryString,tempString);
  }
  strcat(queryString," from photglobal;\n");

  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (gotAnswer) {
          printf("ERROR: two photglobal entries\n");
        }
        gotAnswer = 1;
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%lf",&pPhotGlobal->epoch);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for epoch\n", nvals);
          }
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%lf",&pPhotGlobal->equinox);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for epoch\n", nvals);
          }
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pPhotGlobal->maxSpatialBins);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for maxSpatialBins\n", nvals);
          }
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pPhotGlobal->keepalive);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for keepalive\n", nvals);
          }
        }
        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&pPhotGlobal->currentVersion);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for currentVersion\n", nvals);
          }
        }
        if (sqlrow[5]) {
          strncpy(pPhotGlobal->timeStamp,sqlrow[5],MAX_TIMESTAMP_STRING);
          pPhotGlobal->timeStamp[MAX_TIMESTAMP_STRING] = 0;
        }

        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%d",&pPhotGlobal->magnitudeFile);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magnitudeFile\n", nvals);
            pPhotGlobal->magnitudeFile = 0;
          }
        } else {
          pPhotGlobal->magnitudeFile = 0;
        }
        if (sqlrow[7]) {
          nvals = sscanf(sqlrow[7],"%d",&pPhotGlobal->solutionNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for solutionNumber\n", nvals);
            pPhotGlobal->solutionNumber = 0;
          }
        } else {
          pPhotGlobal->solutionNumber = 0;
        }
        for (index = 0; index < MAX_CATALOG_NUMBER; index++) {
          if (sqlrow[8+index]) {
            nvals = sscanf(sqlrow[8+index],"%d",&pPhotGlobal->minVersionId[index]);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for minVersionId%d\n", nvals, index);
              pPhotGlobal->minVersionId[index] = 0;
            }
          } else {
            pPhotGlobal->minVersionId[index] = 0;
          }

        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }

  /* If we were successful, query the version files for the latest version information */
  if (gotAnswer != 1) {
    printf("ERROR: no entries found in photglobal\n");
    return(gotAnswer);
  }

  gotAnswer = 0;
  sprintf(queryString,"SELECT versionDate,versionName from photversions where versionId = '%d'\n;",pPhotGlobal->currentVersion);
  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (gotAnswer) {
          printf("ERROR: two entries for version %d in photversions\n",pPhotGlobal->currentVersion);
        }
        gotAnswer = 1;
        if (sqlrow[0]) {
          strncpy(pPhotGlobal->versionDate,sqlrow[0],MAX_TIMESTAMP_STRING);
          pPhotGlobal->versionDate[MAX_TIMESTAMP_STRING] = 0;
        }
        if (sqlrow[1]) {
          strncpy(pPhotGlobal->versionName,sqlrow[1],MAX_VERSION_NAME);
          pPhotGlobal->versionName[MAX_VERSION_NAME] = 0;
        }

      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }

  if (gotAnswer != 1) {
    printf("ERROR: %d entries found in photversions\n",gotAnswer);
  }

  return(gotAnswer);
}


int
UpdatePhotometryGlobal(
  MYSQL *pPhotConnection,
  PPHOTGLOBAL pPhotGlobal,
  int new_version,
  char *versionname
) {
  int res;
  char queryString[MAX_QUERY_STRING];
  int gotAnswer = 0;
  int nextVersion;
  time_t tp;
  struct tm *ptr;


  if (new_version) {
    /* Get our current version */
    gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
    if (gotAnswer != 1) {
      return(gotAnswer);
    }
    time (&tp);
    ptr = gmtime(&tp);
    strftime(pPhotGlobal->versionDate,MAX_TIMESTAMP_STRING, "%Y-%m-%dT%H-%M-%S", ptr);
    pPhotGlobal->versionDate[MAX_TIMESTAMP_STRING] = 0;
    /* New version, so bump the number */
    nextVersion = ++pPhotGlobal->currentVersion;
    strncpy(pPhotGlobal->versionName,versionname,MAX_VERSION_NAME);
    pPhotGlobal->versionName[MAX_VERSION_NAME] = 0;

    sprintf(queryString,"INSERT IGNORE into photversions (versionId,versionName,versionDate) values (%d,'%s','%s');\n",
            nextVersion,pPhotGlobal->versionName,pPhotGlobal->versionDate);
    res = ExecuteQuery(pPhotConnection,queryString);
    if (res) {
      gotAnswer = 0;
      return(gotAnswer);
    }
    sprintf(queryString,"UPDATE  photglobal set currentVersion = %d;\n",nextVersion);
    res = ExecuteQuery(pPhotConnection,queryString);
    if (res) {
      gotAnswer = 0;
      return(gotAnswer);
    }

  }
  gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
  return(gotAnswer);
}


int
TimeStampPhotometryGlobal(MYSQL *pPhotConnection)
{
  int res;
  char queryString[MAX_QUERY_STRING];
  time_t tp;
  struct tm *ptr;
  char timestr[MAX_TIMESTAMP_STRING+1];

  time (&tp);
  ptr = gmtime(&tp);
  strftime(timestr,MAX_TIMESTAMP_STRING, "%Y-%m-%dT%H-%M-%S", ptr);
  timestr[MAX_TIMESTAMP_STRING] = 0;
  sprintf(queryString,"UPDATE  photglobal set timeStamp = '%s';\n",timestr);
  res = ExecuteQuery(pPhotConnection,queryString);
  return(res);
}


int
GetPlateScale(
  MYSQL *pConnection,
  char * series,
  int plateNumber,
  int mosaicNumber,
  int solutionNumber,
  int *pMosaicWidth,
  int *pMosaicHeight,
  double *pScale
) {
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int nvals;
  double cd1_1;
  double cd2_2;
  double cd2_1;
  double cd1_2;
  int naxis1;
  int naxis2;
  char *ctype1;
  char *ctype2;
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  int WCSSource;
  double dRightAscension;
  double dDeclination;
  int rotation;
  int scanNumber;
  int transform;
  double cd[4];
  double xctr;
  double yctr;
  struct WorldCoor *wcs;
  double cra;
  double cdec;
  double era;
  double edec;
  double dist;
  double pixdist;
  double xmin;
  double ymin;
  double scale;

  char queryString[MAX_QUERY_STRING];
  int foundEntry = 0;

  *pMosaicWidth = 0;
  *pMosaicHeight = 0;
  *pScale = 0;
  /* Bugfix of Jun 13, 2011: change "exposureNumber" below to "solutionNumber" */
  sprintf(queryString,"SELECT cd1_1,cd1_2,cd2_1,cd2_2,naxis1,naxis2,ctype1,ctype2,crval1,crval2,crpix1,crpix2,WCSSource+0,rotation,dRightAscension,dDeclination,scanNumber,transform+0  FROM mosaics INNER JOIN exposures USING (series,plateNumber,exposureNumber) where series = '%s' and plateNumber = %d and mosaicNumber = %d and solutionNumber = %d;",series,plateNumber,mosaicNumber,solutionNumber);
  if (strlen(queryString) > MAX_QUERY_STRING) {
    printf("ERROR: MAX_QUERY_STRING (1) exceeded %zu\n", strlen(queryString));
    exit(1);
  }


  res = mysql_query(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);

    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (foundEntry > 0) {
          printf("ERROR: Two database entries found for %s%05d\n",series,plateNumber);
        }
        if (sqlrow[0] == NULL) {
          printf("ERROR: No WCS found for %s%05d\n",series,plateNumber);
          continue;
        }
        nvals = sscanf(sqlrow[0],"%lf",&cd1_1);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for cd1_1\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[1],"%lf",&cd1_2);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for cd1_2\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[2],"%lf",&cd2_1);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for cd2_1\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[3],"%lf",&cd2_2);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for cd2_2\n",nvals);
          continue;
        }

        nvals = sscanf(sqlrow[4],"%d",&naxis1);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for naxis1\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[5],"%d",&naxis2);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for naxis2\n",nvals);
          continue;
        }

        ctype1 = sqlrow[6];
        ctype2 = sqlrow[7];


        nvals = sscanf(sqlrow[8],"%lf",&crval1);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for crval1\n",nvals);
          continue;
        }

        nvals = sscanf(sqlrow[9],"%lf",&crval2);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for crval2\n",nvals);
          continue;
        }

        nvals = sscanf(sqlrow[10],"%lf",&crpix1);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for crpix1\n",nvals);
          continue;
        }

        nvals = sscanf(sqlrow[11],"%lf",&crpix2);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for crpix2\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[12],"%d",&WCSSource);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for WCSSource\n",nvals);
          continue;
        }

        if (sqlrow[13] != NULL) {
          nvals = sscanf(sqlrow[13],"%d",&rotation);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for rotation\n",nvals);
            continue;
          }
        } else {
          if (WCSSource < WCSSOURCE_IMWCS) {
            rotation = 0;
          } else {
            printf("ERROR: rotation field is null\n");
            continue;
          }
        }

        nvals = sscanf(sqlrow[14],"%lf",&dRightAscension);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for dRightAscension\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[15],"%lf",&dDeclination);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for dDeclination\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[16],"%d",&scanNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for scanNumber\n",nvals);
          continue;
        }
        if (sqlrow[17] != NULL) {
          nvals = sscanf(sqlrow[17],"%d",&transform);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for transform\n",nvals);
            continue;
          }
        } else {
          transform = 0;
        }


        foundEntry++;

      }

      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }

  } else {
    printf("ERROR: Select error %d: %s res %d\n",mysql_errno(pConnection),mysql_error(pConnection),res);
  }

  if (foundEntry != 1) {
    printf("ERROR: %d plate entry found for %s%05d_%02d s%d %s\n",foundEntry,series,plateNumber,mosaicNumber,solutionNumber,queryString);
    return(1);
  }
  if ((naxis1 == 0) ||
      (naxis2 == 0)) {
    printf("ERROR: naxis1 %d or naxis2 %d is zero for %s%05d\n", naxis1, naxis2, series, plateNumber);
    return(1);
  }

  if (strstr(ctype1,"DEC")) {
    char *tmpPtr;
    double dtmp;
    tmpPtr = ctype2;
    ctype2 = ctype1;
    ctype1 = tmpPtr;
    dtmp = crval1;
    crval1 = crval2;
    crval2 = dtmp;

    cd[0] = cd2_1;
    cd[1] = cd2_2;
    cd[2] = cd1_1;
    cd[3] = cd1_2;
  } else {
    cd[0] = cd1_1;
    cd[1] = cd1_2;
    cd[2] = cd2_1;
    cd[3] = cd2_2;
  }

  wcs = wcskinit(naxis1,
                 naxis2,
                 ctype1,
                 ctype2,
                 crpix1,
                 crpix2,
                 crval1,
                 crval2,
                 cd,
                 0,  /* cdelt1 */
                 0,  /* cdelt2 */
                 0,  /* crota */
                 2000, /* equinox */
                 0);   /* epoch */

  xctr = 0.5 + (0.5 *naxis1);
  yctr = 0.5 + (0.5 *naxis2);
  pix2wcs(wcs,xctr,yctr,&cra,&cdec);
  xmin = 0.5;
  ymin = 0.5;
  pix2wcs(wcs,xmin,ymin,&era,&edec);

  dist = wcsdist(cra,cdec,era,edec);

  pixdist = sqrt(((xctr-xmin)*(xctr-xmin)) + ((yctr-ymin)*(yctr-ymin)));
  if (pixdist > 0.0) {
    scale = dist/pixdist;
  } else {
    scale = 0;
    return(1);
  }

  *pMosaicWidth = naxis1;
  *pMosaicHeight = naxis2;
  *pScale = 3600.0 * scale;
  return(0);
}


/* Get a block of stars for processing.  This routine returns the number of stars found.  If maxstars is zero, then
   the routine will allocate the required block of stars
*/
int
GetStarEntry(
  MYSQL *pPhotConnection,
  PSTARENTRY *ppStarEntry,
  int maxstars,
  int updateFlag,
  char *REF,
  int gsc_bin_index,
  char *catalogString,
  int verbose,
  int nospace
) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  PSTARENTRY pStarEntry = *ppStarEntry;
  char queryString[MAX_QUERY_STRING];
  char tmpString[MAX_QUERY_STRING];
  int curIndex = 0;
  PSTARENTRY pCurStarEntry;
  long long REFNumber;
  int refType;
  static int gsc_bin_indexPrint = 0;

  if (maxstars > 0) {
    memset(pStarEntry,0,maxstars *sizeof(STARENTRY));
  }
  if ((REF != NULL) && (strlen(REF) > 0)) {
    GetREFNumber(REF,&REFNumber,&refType,nospace,verbose);
    sprintf(queryString,  "SELECT gsc_bin_index,versionId,Stdmag,color,ra,declination,REFNumber,updateFlag+0,MAGFlag,gscclass,VFlag,RaPM,DecPM from stars%s where REFNumber = %lld order by gsc_bin_index",catalogString,REFNumber);

  } else if (gsc_bin_index != 0) {
    if (gsc_bin_indexPrint == 0) {
      printf("ERROR: GetStarEntry is performing a gsc_bin_index search\n");
      gsc_bin_indexPrint = 1;
    }
    if (updateFlag == UPDATEFLAG_NEED_UPDATE) {
      sprintf(queryString,"SELECT gsc_bin_index,versionId,Stdmag,color,ra,declination,REFNumber,updateFlag+0,MAGFlag,gscclass,VFlag,RaPM,DecPM from stars%s where updateFlag = 'yes' and gsc_bin_index = %d",catalogString,gsc_bin_index);
    } else {
      sprintf(queryString,"SELECT gsc_bin_index,versionId,Stdmag,color,ra,declination,REFNumber,updateFlag+0,MAGFlag,gscclass,VFlag,RaPM,DecPM from stars%s where gsc_bin_index = %d",catalogString,gsc_bin_index);

    }
  } else {
    if (updateFlag == UPDATEFLAG_NEED_UPDATE) {
      sprintf(queryString,"SELECT gsc_bin_index,versionId,Stdmag,color,ra,declination,REFNumber,updateFlag+0,MAGFlag,gscclass,VFlag,RaPM,DecPM from stars%s where updateFlag = 'yes' order by gsc_bin_index",catalogString);
    } else {
      sprintf(queryString,"SELECT gsc_bin_index,versionId,Stdmag,color,ra,declination,REFNumber,updateFlag+0,MAGFlag,gscclass,VFlag,RaPM,DecPM from stars%s order by gsc_bin_index",catalogString);

    }
  }
  if (maxstars > 0) {
    sprintf(tmpString," limit %d;\n",maxstars);
  } else {
    sprintf(tmpString,";\n");

  }
  strcat(queryString,tmpString);
  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      if (maxstars == 0) {
        maxstars = mysql_affected_rows(pPhotConnection);
        if (maxstars > 0) {
          pStarEntry = (PSTARENTRY)calloc(maxstars,sizeof(STARENTRY));
          if (pStarEntry == NULL) {
            printf("ERROR: Failed to allocate pStarEntry\n");
            exit(1);
          }


        } else {
          pStarEntry = NULL;
        }
        *ppStarEntry = pStarEntry;
      }


      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (curIndex >= maxstars) {
          printf("ERROR: GetStarEntry requested %d stars and got %d stars\n",maxstars,curIndex);
          exit(1);
        }
        pCurStarEntry = &pStarEntry[curIndex];
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pCurStarEntry->gsc_bin_index);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for gsc_bin_index\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->gsc_bin_index = 0;
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%d",&pCurStarEntry->versionId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for versionId\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->versionId = 0;
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%lf",&pCurStarEntry->Stdmag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for Stdmag\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->Stdmag = 0.0;
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%lf",&pCurStarEntry->color);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for color\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->color = 0;
        }
        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%lf",&pCurStarEntry->ra);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ra\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->ra = 0;
        }
        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%lf",&pCurStarEntry->dec);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for dec\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->dec = 0;
        }


        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%lld",&pCurStarEntry->REFNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for REFNumber\n", nvals);
            continue;
          } else {
            GetREF(pCurStarEntry->REFNumber,pCurStarEntry->REF,nospace,1);
          }

        } else {
          printf("ERROR: REFNumber is blank in GetStarEntry");
          continue;
        }
        if (sqlrow[7]) {
          nvals = sscanf(sqlrow[7],"%d",&pCurStarEntry->updateflag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for updateflag\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->updateflag = UPDATEFLAG_IGNORE;
        }


        if (sqlrow[8]) {
          nvals = sscanf(sqlrow[8],"%d",&pCurStarEntry->MAGFlag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for MAGFlag\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->MAGFlag = 0;
        }
        if (sqlrow[9]) {
          nvals = sscanf(sqlrow[9],"%d",&pCurStarEntry->class);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for class\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->class = 0;
        }
        if (sqlrow[10]) {
          nvals = sscanf(sqlrow[10],"%d",&pCurStarEntry->VFlag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for VFlag\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->VFlag = 0;
        }

        if (sqlrow[11]) {
          nvals = sscanf(sqlrow[11],"%lf",&pCurStarEntry->RaPM);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for RaPM\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->RaPM = 0;
        }

        if (sqlrow[12]) {
          nvals = sscanf(sqlrow[12],"%lf",&pCurStarEntry->DecPM);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for DecPM\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->DecPM = 0;
        }

        curIndex++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return(-1);
  }
  return(curIndex);
}

/* Get the location for a single star that does not have the location already
 * embedded in the starID (such as DASCH and APASS catalog stars) This routine
 * bypasses the deprecated stars tables in favor of the new starcatalog table.
 * Do not return entries that have no valid positions.
 */
int
GetStarEntry2(
  MYSQL *pPhotConnection,
  PSTARENTRY pCurStarEntry,
  char *REF,
  int verbose,
  int nospace
) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  int curIndex = 0;
  int refType;
  int maxstars = 1;

  memset(pCurStarEntry,0,maxstars *sizeof(STARENTRY));
  strcpy(pCurStarEntry->REF,REF);

  if ((REF != NULL) && (strlen(REF) > 0)) {
    GetREFNumber(REF,&pCurStarEntry->REFNumber,&refType,nospace,verbose);
    sprintf(queryString,  "SELECT ra,declination,Stdmag,color,gscclass,VFlag,MAGFlag,RaPM,DecPM,gsc_bin_index from starcatalog where REFNumber = %lld;",pCurStarEntry->REFNumber);
  } else {
    printf("ERROR: Unable to decode REF in GetStarEntry2");
    return(-1);
  }
  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (curIndex >= maxstars) {
          printf("ERROR: GetStarEntry requested %d stars and got %d stars\n",maxstars,curIndex);
          exit(1);
        }
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%lf",&pCurStarEntry->ra);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ra\n", nvals);
            continue;
          }
          if ((pCurStarEntry->ra < 0) || (pCurStarEntry->ra >= 360.0)) {
            continue;
          }
        } else {
          pCurStarEntry->ra = 0;
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%lf",&pCurStarEntry->dec);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for dec\n", nvals);
            continue;
          }
          if ((pCurStarEntry->dec < -90.0) || (pCurStarEntry->dec > +90.0)) {
            continue;
          }
        } else {
          pCurStarEntry->dec = 0;
        }

        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%lf",&pCurStarEntry->Stdmag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for Stdmag\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->Stdmag = 0.0;
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%lf",&pCurStarEntry->color);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for color\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->color = 0;
        }

        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&pCurStarEntry->class);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for class\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->class = 0;
        }
        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%d",&pCurStarEntry->VFlag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for VFlag\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->VFlag = 0;
        }

        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%d",&pCurStarEntry->MAGFlag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for MAGFlag\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->MAGFlag = 0;
        }

        if (sqlrow[7]) {
          nvals = sscanf(sqlrow[7],"%lf",&pCurStarEntry->RaPM);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for RaPM\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->RaPM = 0;
        }

        if (sqlrow[8]) {
          nvals = sscanf(sqlrow[8],"%lf",&pCurStarEntry->DecPM);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for DecPM\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->DecPM = 0;
        }

        if (sqlrow[9]) {
          nvals = sscanf(sqlrow[9],"%d",&pCurStarEntry->gsc_bin_index);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for gsc_bin_index\n", nvals);
            continue;
          }
        } else {
          pCurStarEntry->gsc_bin_index = -1;
        }

        pCurStarEntry->updateflag =  UPDATEFLAG_NO_UPDATE;
        pCurStarEntry->versionId = -1;
        curIndex++;


      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return(-1);
  }
  return(curIndex);
}


/* Get the location from the 'identifier' table which is a local name resolver database */
int
GetStarEntry3(
  MYSQL *pPhotConnection,
  char *catalogname,
  double *pra,
  double *pdeclination
) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  int curIndex = 0;
  int maxstars = 1;

  *pra = 999;
  *pdeclination = 999;

  sprintf(queryString,"SELECT ra,declination from identifier where catalogname = '%s';",catalogname);
  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (curIndex >= maxstars) {
          printf("ERROR: GetStarEntry requested %d stars and got %d stars\n",maxstars,curIndex);
        }
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%lf",pra);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ra\n", nvals);
            continue;
          }
          if ((*pra < 0) || (*pra >= 360.0)) {
            continue;
          }
        } else {
          *pra = 999.;
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%lf",pdeclination);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for dec\n", nvals);
            continue;
          }
          if ((*pdeclination < -90.0) || (*pdeclination > +90.0)) {
            continue;
          }
        } else {
          *pdeclination = 0;
        }

        curIndex++;


      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return(-1);
  }
  return(curIndex);
}

int
ReadPhotSpatialBin(
  MYSQL *pPhotConnection,
  PPHOT_SPATIAL_BIN *ppSpatialBinTable,
  char *series,
  int plateNumber,
  int solutionNumber,
  int spatial_bin,
  char *catalogString,
  int solutionNumberFix
) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  char tempString[MAX_QUERY_STRING];
  PPHOT_SPATIAL_BIN pSpatialBinTable = *ppSpatialBinTable;
  PPHOT_SPATIAL_BIN pCurSpatialBin;
  int numRows = 0;
  int curRow = 0;
  int seriesId;
  int whereFlag = 0;
  char solutionNumberString[20];

  if (solutionNumberFix) {
    strcpy(solutionNumberString,"solutionNumber");
  } else {
    strcpy(solutionNumberString,"exposureNumber");
  }

  sprintf(queryString,"SELECT seriesId,plateNumber,%s,spatial_bin,rms,n1,n2,limiting_mag,limiting_iso,max_bright_mag,upper_limit,colorterm,errorcolor,colorflag+0,max_bright_iso,versionId from spatialbin%s ",solutionNumberString,catalogString);

  if (series != NULL) {
    seriesId = GetSeriesId(series,1);
    sprintf(tempString,"WHERE seriesId = %d ",seriesId);
    strcat(queryString,tempString);
    whereFlag = 1;
  }
  if (plateNumber != 0) {
    if (whereFlag) {
      strcat(queryString," AND ");
    } else {
      strcat(queryString," WHERE ");
      whereFlag = 1;
    }
    sprintf(tempString,"plateNumber  = %d ",plateNumber);
    strcat(queryString,tempString);
  }
  if (solutionNumber > 0) {
    if (whereFlag) {
      strcat(queryString," AND ");
    } else {
      strcat(queryString," WHERE ");
      whereFlag = 1;
    }
    sprintf(tempString,"%s  = %d ",solutionNumberString,solutionNumber);
    strcat(queryString,tempString);
  }
  if (spatial_bin > 0) {
    if (whereFlag) {
      strcat(queryString," AND ");
    } else {
      strcat(queryString," WHERE ");
      whereFlag = 1;
    }
    sprintf(tempString,"spatial_bin  = %d ",spatial_bin);
    strcat(queryString,tempString);
  }
  strcat(queryString,";");

  if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
    printf("ERROR: queryString exceeded in ReadPhotSpatialBin by %zu\n", strlen(queryString));
    exit(1);
  }

  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pPhotConnection);
      pSpatialBinTable = (PPHOT_SPATIAL_BIN)calloc(numRows,sizeof(PHOT_SPATIAL_BIN));
      if (pSpatialBinTable == NULL) {
        printf("ERROR: failed to allocate %d PHOT_SPATIAL_BIN rows in ReadPhotSpatialBin\n",numRows);
        exit(1);
      }
      *ppSpatialBinTable = pSpatialBinTable;


      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        pCurSpatialBin = &pSpatialBinTable[curRow];
        memset(pCurSpatialBin,0,sizeof(PHOT_SPATIAL_BIN));
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pCurSpatialBin->seriesId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for seriesId in ReadPhotSpatialBin\n", nvals);
            continue;
          } else {
            strcpy(pCurSpatialBin->series,GetSeriesString(pCurSpatialBin->seriesId,1));
          }
        } else {
          printf("ERROR: NULL seriesId in ReadPhotSpatialBin\n");
          continue;
        }


        if (sqlrow[1]) {

          nvals = sscanf(sqlrow[1],"%d",&pCurSpatialBin->plateNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for plateNumber\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->plateNumber = 0;
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pCurSpatialBin->solutionNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for solutionNumber\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->solutionNumber = 0;
        }

        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pCurSpatialBin->spatial_bin);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for spatial_bin\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->spatial_bin = 0;
        }

        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%lf",&pCurSpatialBin->rms);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for rms\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->rms = 0.0;
        }

        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%d",&pCurSpatialBin->n1);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for n1\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->n1 = 0;
        }


        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%d",&pCurSpatialBin->n2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for n2\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->n2 = 0;
        }


        if (sqlrow[7]) {
          nvals = sscanf(sqlrow[7],"%lf",&pCurSpatialBin->limiting_mag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for limiting_mag\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->limiting_mag = 0.0;
        }

        if (sqlrow[8]) {
          nvals = sscanf(sqlrow[8],"%lf",&pCurSpatialBin->limiting_iso);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for limiting_iso\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->limiting_iso = 0.0;
        }



        if (sqlrow[9]) {
          nvals = sscanf(sqlrow[9],"%lf",&pCurSpatialBin->max_bright_mag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for max_bright_mag\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->max_bright_mag = 0.0;
        }

        if (sqlrow[10]) {
          nvals = sscanf(sqlrow[10],"%lf",&pCurSpatialBin->upper_limit);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for upper_limit\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->upper_limit = 0.0;
        }

        if (sqlrow[11]) {
          nvals = sscanf(sqlrow[11],"%lf",&pCurSpatialBin->colorterm);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for colorterm\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->colorterm = 0.0;
        }

        if (sqlrow[12]) {
          nvals = sscanf(sqlrow[12],"%lf",&pCurSpatialBin->errorcolor);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for errorcolor\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->errorcolor = 0.0;
        }


        if (sqlrow[13]) {
          nvals = sscanf(sqlrow[13],"%d",&pCurSpatialBin->colorflag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for colorflag\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->colorflag = 0;
        }

        if (sqlrow[14]) {
          nvals = sscanf(sqlrow[14],"%lf",&pCurSpatialBin->max_bright_iso);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for max_bright_iso\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->max_bright_iso = 0.0;
        }
        if (sqlrow[15]) {
          nvals = sscanf(sqlrow[15],"%d",&pCurSpatialBin->versionId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for versionId\n", nvals);
            continue;
          }
        } else {
          pCurSpatialBin->versionId = 0;
        }

        curRow++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }

  return(curRow);
}


int
GetPhotLocalBin(
  MYSQL *pPhotConnection,
  LOCALBIN *pLocalBinTable,
  int localBinTableSize,
  char *series,
  int plateNumber,
  int solutionNumber,
  int local_bin_index,
  char *catalogString,
  int solutionNumberFix,
  int debugMode
) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  int numRows = 0;
  int curRow = 0;
  int index;
  LOCALBIN curLocalBin;
  PLOCALBIN pCurLocalBin = &curLocalBin;
  PLOCALBIN pTargetLocalBin;
  int seriesId;
  char solutionNumberString[20];

  if (solutionNumberFix) {
    strcpy(solutionNumberString,"solutionNumber");
  } else {
    strcpy(solutionNumberString,"exposureNumber");
  }

  seriesId = GetSeriesId(series,1);
  memset(pLocalBinTable,0,localBinTableSize * sizeof(LOCALBIN));
  for (index = 0; index < localBinTableSize; index++) {
    pTargetLocalBin = &pLocalBinTable[index];
    pTargetLocalBin->local_bin_index = -1;
  }
  if (local_bin_index >= 0) {
    if (localBinTableSize != 1) {
      printf("ERROR: localBinTableSize %d must be 1 if local_bin_index %d is specified\n",localBinTableSize,local_bin_index);
      exit(1);
    }
    sprintf(queryString,"SELECT seriesId,plateNumber,%s,local_bin_index,versionId,altitude,extinction,magcor_local,magcal_local_error,npoints_local,rejectFlag,drad_bin_count,drad_bin_size,drad_reject_count,draMedian,draRMS,ddecMedian,ddecRMS,drad_bin_count2,drad_bin_size2,drad_reject_count2,draMedian2,draRMS2,ddecMedian2,ddecRMS2,dradRMS2 from localbin%s WHERE seriesId = %d and plateNumber = %d and %s = %d and local_bin_index = %d;",solutionNumberString,catalogString,seriesId,plateNumber,solutionNumberString,solutionNumber,local_bin_index);

  } else {

    sprintf(queryString,"SELECT seriesId,plateNumber,%s,local_bin_index,versionId,altitude,extinction,magcor_local,magcal_local_error,npoints_local,rejectFlag,drad_bin_count,drad_bin_size,drad_reject_count,draMedian,draRMS,ddecMedian,ddecRMS,drad_bin_count2,drad_bin_size2,drad_reject_count2,draMedian2,draRMS2,ddecMedian2,ddecRMS2,dradRMS2 from localbin%s WHERE seriesId = %d and plateNumber = %d and %s = %d;",solutionNumberString,catalogString,seriesId,plateNumber,solutionNumberString,solutionNumber);
  }

  if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
    printf("ERROR: queryString exceeded in GetPhotLocalBin by %zu\n", strlen(queryString));
    exit(1);
  }

  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pPhotConnection);
      if (numRows != localBinTableSize) {
        if (localBinTableSize != 1) {
          if (debugMode == 0) {
            printf("ERROR: Plate %s%05d expects %d localbin rows and has %d localbin rows\n",series,plateNumber,localBinTableSize,numRows);
          }
        }
        mysql_free_result(res_ptr);
        return(1);
      }

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        memset(pCurLocalBin,0,sizeof(LOCALBIN));

        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pCurLocalBin->seriesId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for seriesId in GetPhotLocalBin\n", nvals);
            continue;
          } else {
            strcpy(pCurLocalBin->series,GetSeriesString(pCurLocalBin->seriesId,1));
          }
        } else {
          printf("ERROR: NULL seriesId in GetPhotLocalBin\n");
          continue;
        }



        if (sqlrow[1]) {

          nvals = sscanf(sqlrow[1],"%d",&pCurLocalBin->plateNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for plateNumber\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->plateNumber = 0;
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pCurLocalBin->solutionNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for solutionNumber\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->solutionNumber = 0;
        }

        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pCurLocalBin->local_bin_index);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for local_bin_index\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->local_bin_index = -1;
        }

        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&pCurLocalBin->versionId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for versionId\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->versionId = 0;
        }

        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%lf",&pCurLocalBin->altitude);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for altitude\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->altitude = 0;
        }


        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%lf",&pCurLocalBin->extinction);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for extinction\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->extinction = 0;
        }


        if (sqlrow[7]) {
          nvals = sscanf(sqlrow[7],"%lf",&pCurLocalBin->magcor_local);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magcor_local\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->magcor_local = 99.0;
        }

        if (sqlrow[8]) {
          nvals = sscanf(sqlrow[8],"%lf",&pCurLocalBin->magcal_local_error);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magcal_local_error\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->magcal_local_error = 99.0;
        }



        if (sqlrow[9]) {
          nvals = sscanf(sqlrow[9],"%d",&pCurLocalBin->npoints_local);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for npoints_local\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->npoints_local = 0;
        }

        if (sqlrow[10]) {
          nvals = sscanf(sqlrow[10],"%d",&pCurLocalBin->rejectFlag);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for rejectFlag\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->rejectFlag = 0;
        }

        if (sqlrow[11]) {
          nvals = sscanf(sqlrow[11],"%d",&pCurLocalBin->drad_bin_count);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for drad_bin_count\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->drad_bin_count = 0;
        }

        if (sqlrow[12]) {
          nvals = sscanf(sqlrow[12],"%d",&pCurLocalBin->drad_bin_size);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for drad_bin_size\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->drad_bin_size = 0;
        }

        if (sqlrow[13]) {
          nvals = sscanf(sqlrow[13],"%d",&pCurLocalBin->drad_reject_count);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for drad_reject_count\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->drad_reject_count = 0.0;
        }


        if (sqlrow[14]) {
          nvals = sscanf(sqlrow[14],"%lf",&pCurLocalBin->draMedian);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for draMedian\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->draMedian = 0.0;
        }

        if (sqlrow[15]) {
          nvals = sscanf(sqlrow[15],"%lf",&pCurLocalBin->draRMS);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for draRMS\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->draRMS = 0.0;
        }

        if (sqlrow[16]) {
          nvals = sscanf(sqlrow[16],"%lf",&pCurLocalBin->ddecMedian);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ddecMedian\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->ddecMedian = 0.0;
        }
        if (sqlrow[17]) {
          nvals = sscanf(sqlrow[17],"%lf",&pCurLocalBin->ddecRMS);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ddecRMS\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->ddecRMS = 0.0;
        }


        if (sqlrow[18]) {
          nvals = sscanf(sqlrow[18],"%d",&pCurLocalBin->drad_bin_count2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for drad_bin_count2\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->drad_bin_count2 = 0;
        }

        if (sqlrow[19]) {
          nvals = sscanf(sqlrow[19],"%d",&pCurLocalBin->drad_bin_size2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for drad_bin_size2\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->drad_bin_size2 = 0;
        }



        if (sqlrow[20]) {
          nvals = sscanf(sqlrow[20],"%d",&pCurLocalBin->drad_reject_count2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for drad_reject_count2\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->drad_reject_count2 = 0.0;
        }


        if (sqlrow[21]) {
          nvals = sscanf(sqlrow[21],"%lf",&pCurLocalBin->draMedian2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for draMedian2\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->draMedian2 = 0.0;
        }

        if (sqlrow[22]) {
          nvals = sscanf(sqlrow[22],"%lf",&pCurLocalBin->draRMS2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for draRMS2\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->draRMS2 = 0.0;
        }

        if (sqlrow[23]) {
          nvals = sscanf(sqlrow[23],"%lf",&pCurLocalBin->ddecMedian2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ddecMedian\2n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->ddecMedian2 = 0.0;
        }
        if (sqlrow[24]) {
          nvals = sscanf(sqlrow[24],"%lf",&pCurLocalBin->ddecRMS2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ddecRMS2\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->ddecRMS2 = 0.0;
        }
        if (sqlrow[25]) {
          nvals = sscanf(sqlrow[25],"%lf",&pCurLocalBin->dradRMS2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ddecRMS2\n", nvals);
            continue;
          }
        } else {
          pCurLocalBin->dradRMS2 = 0.0;
        }
        if (local_bin_index < 0) {
          if ((pCurLocalBin->local_bin_index < 0) ||
              (pCurLocalBin->local_bin_index >= localBinTableSize)) {
            printf("ERROR: plate %s%05d has an illegal local bin index of %d\n",series,plateNumber,pCurLocalBin->local_bin_index);
            return(1);
          }
          pTargetLocalBin = &pLocalBinTable[pCurLocalBin->local_bin_index];
          if (pTargetLocalBin->local_bin_index >= 0) {
            printf("ERROR: plate %s%05d has two local bin entries with index %d\n",series,plateNumber,pCurLocalBin->local_bin_index);
            return(1);
          }
        } else {
          if (pCurLocalBin->local_bin_index < local_bin_index) {
            printf("ERROR: plate %s%05d has an illegal local bin index of %d\n",series,plateNumber,pCurLocalBin->local_bin_index);
            return(1);
          }
          pTargetLocalBin = pLocalBinTable;
        }
        memcpy(pTargetLocalBin,pCurLocalBin,sizeof(LOCALBIN));
        curRow++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }
  if (curRow != localBinTableSize) {
    printf("ERROR: plate %s%05d has %d local bins instead of %d\n",series,plateNumber,curRow,localBinTableSize);
    return(1);
  }

  return(0);
}


void
InitMaxPlateNumber(MYSQL *pConnection, int *maxPlateNumber)
{
  PRANGES pCurRange;
  int rangeIndex;
  int seriesId;
  GetRangeTable(pConnection);
  for (rangeIndex = 0; rangeIndex < numRanges; rangeIndex++) {
    pCurRange = &pRangeTable[rangeIndex];
    seriesId = GetSeriesId(pCurRange->series,1);
    if (maxPlateNumber[seriesId] < pCurRange->lastPlate) {
      maxPlateNumber[seriesId] = pCurRange->lastPlate;
    }

  }
}


void
InitFileCommon(
  FILE *consoleHandle,
  PFILECOMMON pFileCommon,
  MYSQL *pConnection,
  MYSQL *pPhotConnection,
  int catalogNumber,
  int verbose
) {
  int cache_index;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;

  assert(MAX_ADJACENT_BINS < CACHE_MAGNITUDE_ENTRIES);
  memset(pFileCommon,0,sizeof(FILECOMMON));
  time(&pFileCommon->startTime);
  if (pPhotConnection != NULL) {
    if (GetPhotometryGlobal(pPhotConnection,pPhotGlobal) != 1) {
      fprintf(consoleHandle,"ERROR: failed to get the global photometry table\n");
      exit(1);
    }
  }
  pFileCommon->cache_current_sequence = 1;
  for (cache_index = 0; cache_index < CACHE_MAGNITUDE_ENTRIES; cache_index++) {
    pFileCommon->cache_gsc_bin_index[cache_index] = -1;
  }
  pFileCommon->pConnection = pConnection;
  pFileCommon->pPhotConnection = pPhotConnection;
  pFileCommon->catalogNumber = catalogNumber;
  pFileCommon->minVersionId = pPhotGlobal->minVersionId[catalogNumber];
  pFileCommon->currentMagnitudesBin = -1;
  pFileCommon->ra0 = 999.;
  pFileCommon->dec0 = 99.;
  pFileCommon->rad0 = -1.0;
  pFileCommon->catalogFD = -1;
  pFileCommon->indexFD = -1;
  pFileCommon->gsc_bin_index = -1;
  /* Initializations for search_none */
  pFileCommon->minimumGalacticLatitude = DEFAULT_GAL_LATITUDE;
  pFileCommon->minimumImageCount = -1;
  pFileCommon->startrank = 1;
  pFileCommon->maxproc = 1;
  time(&pFileCommon->beginTime);

  if (verbose) {
    fprintf(
      consoleHandle,
      "sizeof FILESTARHEADER %zu sizeof FILESTARIMAGE %zu\n",
      sizeof(FILESTARHEADER),
      sizeof(FILESTARIMAGE)
    );
  }

#ifdef USE_DEBUG_HANDLE
  {
    char debugName[] = "/home/scanner/Pipeline/losdebug.db";
    pFileCommon->debugHandle = fopen(debugName,"wt");
    if (pFileCommon->debugHandle == NULL) {
      fprintf(consoleHandle,"ERROR: Failed to open the debug handle %s\n",debugName);
    } else {
      fprintf(consoleHandle,"ERROR: Active debug handle %s\n",debugName);
    }
  }
#endif /* USE_DEBUG_HANDLE */
}


void
FreeFileCommon(PFILECOMMON pFileCommon, int verbose)
{
  int seriesId;
  int cache_index;
  PGALAXYCOMMON pGalaxyCommon = &pFileCommon->galaxyCommon;
  FreeGalaxyCommon(pGalaxyCommon,verbose);
  if (pFileCommon->magnitudeBuffer != NULL) {
    free(pFileCommon->magnitudeBuffer);
  }
  if (pFileCommon->pGscImageTable != NULL) {
    free(pFileCommon->pGscImageTable);
  }

  if (pFileCommon->gscBinIndexCount > 0) {
    if (pFileCommon->gscBinIndexList != NULL) {
      free(pFileCommon->gscBinIndexList);
      pFileCommon->gscBinIndexList = NULL;
    }
    pFileCommon->gscBinIndexCount = 0;
  }


  for (seriesId = 0; seriesId < MAX_SERIES; seriesId++) {
    float **seriesLimitingMagnitudeArray = pFileCommon->seriesLimitingMagnitudeArray;

    if (pFileCommon->seriesVersionTable[seriesId] != NULL) {
      free(pFileCommon->seriesVersionTable[seriesId]);
    }
    if (pFileCommon->seriesQualityMaskTable[seriesId] != NULL) {
      free(pFileCommon->seriesQualityMaskTable[seriesId]);
    }
    if (pFileCommon->seriesQualityValueTable[seriesId] != NULL) {
      free(pFileCommon->seriesQualityValueTable[seriesId]);
    }
    if (pFileCommon->seriesColortermTable[seriesId] != NULL) {
      free(pFileCommon->seriesColortermTable[seriesId]);
    }
    if (pFileCommon->seriesColorflagTable[seriesId] != NULL) {
      free(pFileCommon->seriesColorflagTable[seriesId]);
    }
    if (pFileCommon->seriesTimeAccuracyTable[seriesId] != NULL) {
      free(pFileCommon->seriesTimeAccuracyTable[seriesId]);
      pFileCommon->seriesTimeAccuracyCount = 0;
      pFileCommon->seriesTimeAccuracyAlloc = 0;
    }

    if (seriesLimitingMagnitudeArray != NULL) {
      if (seriesLimitingMagnitudeArray[seriesId]  != NULL) {
        free(seriesLimitingMagnitudeArray[seriesId]);
        seriesLimitingMagnitudeArray[seriesId] = NULL;
      }
    }
  }
  if (pFileCommon->seriesLimitingMagnitudeArray != NULL) {
    free(pFileCommon->seriesLimitingMagnitudeArray);
  }
  for (cache_index = 0; cache_index < CACHE_MAGNITUDE_ENTRIES; cache_index++) {
    pFileCommon->cache_gsc_bin_index[cache_index] = -1;
    if (pFileCommon->pCacheMagnitudeTable[cache_index] != NULL) {
      free(pFileCommon->pCacheMagnitudeTable[cache_index]);
      pFileCommon->pCacheMagnitudeTable[cache_index] = NULL;
    }
    if (pFileCommon->pCacheFileMagnitudeTable[cache_index] != NULL) {
      free(pFileCommon->pCacheFileMagnitudeTable[cache_index]);
      pFileCommon->pCacheFileMagnitudeTable[cache_index] = NULL;
    }
  }

  if (pFileCommon->none_stats_table != NULL) {
    free(pFileCommon->none_stats_table);
    pFileCommon->none_stats_table = NULL;
  }
  if (pFileCommon->pm_stats_table != NULL) {
    free(pFileCommon->pm_stats_table);
    pFileCommon->pm_stats_table = NULL;
  }
  if (pFileCommon->area_sort_table != NULL) {
    free(pFileCommon->area_sort_table);
    pFileCommon->area_sort_table = NULL;
  }

  if (pFileCommon->group_stats_table != NULL) {
    free(pFileCommon->group_stats_table);
    pFileCommon->group_stats_table = NULL;
  }

  if (pFileCommon->nearby_catalog_table != NULL) {
    free(pFileCommon->nearby_catalog_table);
    pFileCommon->nearby_catalog_table = NULL;
    pFileCommon->nearbyCatalogCount = 0;
  }

  if (pFileCommon->vector1 != NULL) {
    free(pFileCommon->vector1);
    pFileCommon->vector1 = NULL;
  }

  if (pFileCommon->vector2 != NULL) {
    free(pFileCommon->vector2);
    pFileCommon->vector2 = NULL;
  }

  if (pFileCommon->vector3 != NULL) {
    free(pFileCommon->vector3);
    pFileCommon->vector3 = NULL;
  }

  if (pFileCommon->catalogFD >= 0) {
    close(pFileCommon->catalogFD);
    pFileCommon->catalogFD = -1;
  }
  if (pFileCommon->indexFD >= 0) {
    close(pFileCommon->indexFD);
    pFileCommon->indexFD = -1;
  }
  if (pFileCommon->catalogTable != NULL) {
    free(pFileCommon->catalogTable);
    pFileCommon->catalogTable = NULL;
    pFileCommon->gsc_bin_index = -1;
    pFileCommon->catalogTableSize = 0;
    pFileCommon->catalogTableAlloc = 0;
  }
  if (pFileCommon->candidate_table != NULL ) {
    free(pFileCommon->candidate_table);
    pFileCommon->candidate_table = NULL;
  }

  if (verbose) {
    printf("FileCommon staleVersionCount %d, versionIdQueries %d, versionIdHits %d, excludeSeriesCount %d, magnitudeBufferAlloc %ld\n",
           pFileCommon->counterBlock.staleVersionCount,
           pFileCommon->counterBlock.versionIdQueries,
           pFileCommon->counterBlock.versionIdHits,
           pFileCommon->counterBlock.excludeSeriesCount,
           pFileCommon->counterBlock.magnitudeBufferAlloc);
    printf("Subroutine calls %d, staleREFCount %d readMagnitudes %lld\n",
           pFileCommon->counterBlock.subroutineCalls,
           pFileCommon->counterBlock.staleREFCount,
           pFileCommon->counterBlock.readMagnitudes);
    printf("FileCommon binRejectCount %d, nomatchRejectCount %d, flagsRejectCount %d,refRejectCount %d chipRejectCount %d, currentMagnitudes %lld\n",
           pFileCommon->counterBlock.binRejectCount,
           pFileCommon->counterBlock.nomatchRejectCount,
           pFileCommon->counterBlock.flagsRejectCount,
           pFileCommon->counterBlock.refRejectCount,
           pFileCommon->counterBlock.chipRejectCount,
           pFileCommon->counterBlock.currentMagnitudes);
    printf("FileCommon noHeaderCount %d, versionOneFile %d, corruptStarCount %d, newDataCount %d, noValidImages %d, qualityRejectCount %d, returnedMagnitudes %lld\n",
           pFileCommon->counterBlock.noHeaderCount,
           pFileCommon->counterBlock.versionOneFile,
           pFileCommon->counterBlock.corruptStarCount,
           pFileCommon->counterBlock.newDataCount,
           pFileCommon->counterBlock.noValidImages,
           pFileCommon->counterBlock.qualityRejectCount,
           pFileCommon->counterBlock.returnedMagnitudes);
    printf("FileCommon filteredVersionCount %d binZeroErr %d mixedCatalogCount %d sequesteredCount %d limitingmag1 %d limitingmag2 %d maxGroupDrad %f rejectedTCGroup %d\n",
           pFileCommon->counterBlock.filteredVersionCount,
           pFileCommon->counterBlock.binZeroErrors,
           pFileCommon->counterBlock.mixedCatalogCount,
           pFileCommon->counterBlock.sequesteredCount,
           pFileCommon->counterBlock.limitingmag1,
           pFileCommon->counterBlock.limitingmag2,
           pFileCommon->counterBlock.maxGroupDrad,
           pFileCommon->counterBlock.rejectedTCGroup);
    if ((pFileCommon->counterBlock.maxPMrarms > 0.0) || (pFileCommon->counterBlock.maxPMdecrms > 0.0)) {
      printf("FileCommon properMotionErrorCount %d maxPMrarms %f maxPMDecrms %f\n",pFileCommon->counterBlock.properMotionErrorCount,pFileCommon->counterBlock.maxPMrarms,pFileCommon->counterBlock.maxPMdecrms);
    }

    if (pFileCommon->counterBlock.cache_access_count > 0) {
      printf("FileCommon cacheAccessCount %d cacheHitCount %d cacheMissCount %d magnitudeReadCount %d magnitudeRejectCount %d\n",
             pFileCommon->counterBlock.cache_access_count,
             pFileCommon->counterBlock.cache_hit_count,
             pFileCommon->counterBlock.cache_miss_count,
             pFileCommon->counterBlock.cache_magnitude_read_count,
             pFileCommon->counterBlock.cache_magnitude_reject_count);

    }

  }

  if (pFileCommon->debugHandle != NULL) {
    fclose(pFileCommon->debugHandle);
  }
  if (pFileCommon->propermotionHandle != NULL) {
    fclose(pFileCommon->propermotionHandle);
  }

  memset(pFileCommon,0,sizeof(FILECOMMON));
}


// Sort routine for the photometry files.
static int
StarImageCompare(const void *first, const void *second)
{
  int gscFirst = ((PFILESTARIMAGE) first)->gsc_bin_index;
  int gscSecond = ((PFILESTARIMAGE) second)->gsc_bin_index;
  long long refFirst = ((PFILESTARIMAGE) first)->REFNumber;
  long long refSecond = ((PFILESTARIMAGE) second)->REFNumber;

  if (gscFirst > gscSecond) {
    return 1;
  }

  if (gscFirst < gscSecond) {
    return -1;
  }

  if (refFirst > refSecond) {
    return 1;
  }

  if (refFirst < refSecond) {
    return -1;
  }

  return 0;
}


static void
AddCatalogReference(PFILECOMMON pFileCommon, PFILESTARIMAGE pFileStarImage)
{
  PGSCIMAGE pCurGscImageX;

  if ((pFileStarImage->REFNumber == 0) ||
      ((pFileStarImage->RaPM < MAX_PROPER_MOTION) &&
       (pFileStarImage->DecPM < MAX_PROPER_MOTION) &&
       (pFileStarImage->ra_2 < 990.0) &&
       (pFileStarImage->dec_2 <= 90.0))) {
    return;
  }
  pCurGscImageX = GetCatalogRecord(pFileCommon,pFileStarImage->REFNumber,pFileStarImage->gsc_bin_index);
  if (pCurGscImageX != NULL) {
    double plateepoch;
    double ra;
    double dec;
    double factor;
    plateepoch = jd2ep(pFileStarImage->Date) - GSC_EQUINOX;
    dec = pCurGscImageX->dec + ((pCurGscImageX->DecPM * plateepoch)/(3600.0*1000.0));
    if ((dec < 90.0) && (dec > -90.0)) {
      factor = cos(dec * DEGREES_TO_RAD);
      ra = pCurGscImageX->ra + (((pCurGscImageX->RaPM *plateepoch)/(3600.0*1000.0))/factor);
    } else {
      ra = pCurGscImageX->ra;
      factor = 1.0;
    }
    pFileStarImage->RaPM = pCurGscImageX->RaPM;
    pFileStarImage->DecPM = pCurGscImageX->DecPM;
    pFileStarImage->ra_2 = ra;
    pFileStarImage->dec_2 = dec;
  }
}


/*
 *  NOTE: verbose = -1 means to just print the file name
 */
int
GetFileSummaryMagnitudes(
  PGSCBIN pGscBin,
  PFILECOMMON pFileCommon,
  int gsc_bin_index_low,
  int gsc_bin_index_high,
  int zeroAFLAGS,
  int REFonly,
  int magnitudeLimit,
  int dumpAllFlag,
  FILE *headerHandle,
  char *catalogString,
  int verbose,
  FILE *repairHandle,
  int readOnly,
  int nosequestered
) {
  FILESTARHEADER fileStarHeader;
  PFILESTARHEADER pFileStarHeader = &fileStarHeader;
  PFILESTARIMAGE pFileStarTable = NULL;
  PFILESTARIMAGE pFileStarImage = NULL;
  PFILESTARIMAGE pCurStarImage;
  int starImageCount;
  int readStarImageCount = 0;
  int writeStarImageCount = 0;
  int versionId;
  off_t curFileOffset;
  int gsc_base_index = gsc_bin_index_low - (gsc_bin_index_low % MAG_FILE_MODULUS);
  char filename[MAX_FILENAME];
  int fileNumber = -1;
  int statResult;
  struct stat filestats;
  off_t filesize;
  ssize_t readBytes;
  ssize_t readBytes2;
  ssize_t writeBytes;
  int result;
  long long header;
  long long headerFlag;
  long long headerVersion;
  char *photfilebase;
  int index;
  int bufferFullFlag = 0;
  off_t readOffset;
  int gsc_bin_index;
  off_t newOffset;
  PHOTPLATES basePhotPlates;
  PPHOTPLATES pPhotPlates = &basePhotPlates;
  char REF[MAX_REF];
  int tmp_gsc_bin_index;
  int decBin;
  int raBin;
  int catalogNumber = CATALOG_GSC232;
  int refType = REF_TYPE_GSC;
  int nvals;
  int writeHeader = 1;
  int magfiles_use_new_layout = dasch_photdb_magfiles_use_new_layout();

  if (catalogString[0] != 0) {
    nvals = sscanf(catalogString, "%d", &catalogNumber);
    if (nvals != 1) {
      printf("ERROR: failed to decode number catalog string %s\n", catalogString);
      exit(1);
    }
  }

  memset(pPhotPlates, 0, sizeof(PHOTPLATES));
  pPhotPlates->mosaicNumber = -1;
  pPhotPlates->versionId = -1;
  pPhotPlates->quality = QUALITY_UNINITIALIZED;

  if (gsc_bin_index_high < gsc_bin_index_low || gsc_bin_index_high >= (gsc_base_index + MAG_FILE_MODULUS)) {
    printf(
      "ERROR: gsc_bin_index_high %d should be between %d and %d\n",
      gsc_bin_index_high,
      gsc_bin_index_low,
      gsc_base_index + MAG_FILE_MODULUS
    );
    exit(1);
  }

  if (gsc_base_index == pFileCommon->currentMagnitudesBin) {
    // Second consecutive call with this index
    pFileCommon->currentMagnitudesBin++;
  } else if (gsc_base_index > pFileCommon->currentMagnitudesBin) {
    // First call with this index
    pFileCommon->currentMagnitudesBin = gsc_base_index;
  }

  photfilebase = GetPhotFileBase(catalogString);

  if (pFileCommon->startTime != 0) {
    time_t newTime;

    time(&newTime);
    newTime -= pFileCommon->startTime;

    if (newTime > WEB_TIMEOUT) {
      printf("ERROR: Request processing timed out in line %d after %ld seconds\n", __LINE__, newTime);
      exit(1);
    }
  }

  if (magfiles_use_new_layout) {
    // gsc_base_index is a multiple of 1024, so to hash effectively:
    int key = (gsc_base_index >> 10) % 400;
    sprintf(filename, "%s/%03d/%010d.dat", photfilebase, key, gsc_base_index);
  } else {
    int dir1 = gsc_base_index / (MAG_SUBDIR_MODULUS * MAG_SUBDIR_MODULUS);
    int dir2 = gsc_base_index / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
    sprintf(filename, "%s/mag%03d/mag%03d/mag%09d.dat", photfilebase, dir1, dir2, gsc_base_index);
  }

  statResult = stat(filename, &filestats);
  if (statResult != 0) {
    // Assume no such file
    pFileCommon->curMagnitudes = 0;
    return 0;
  }

  pFileCommon->counterBlock.subroutineCalls++;
  filesize = filestats.st_size;

  if (
    pFileCommon->counterBlock.magnitudeBufferAlloc < (filesize + sizeof(FILESTARHEADER))
  ) {
    pFileCommon->counterBlock.magnitudeBufferAlloc = filesize + sizeof(FILESTARHEADER);

    if (pFileCommon->magnitudeBuffer != NULL) {
      free(pFileCommon->magnitudeBuffer);
    }

    pFileCommon->magnitudeBuffer = (PFILESTARIMAGE) calloc(pFileCommon->counterBlock.magnitudeBufferAlloc, sizeof(char));

    if (pFileCommon->magnitudeBuffer == NULL) {
      printf(
        "ERROR: failed to allocate pFileCommon->magnitudeBuffer of size %ld\n",
        pFileCommon->counterBlock.magnitudeBufferAlloc
      );
      exit(1);
    }
  }

  pFileStarTable = (PFILESTARIMAGE) pFileCommon->magnitudeBuffer;
  fileNumber = open(filename, O_RDONLY, 0664);
  readBytes = read(fileNumber, pFileStarHeader, sizeof(FILESTARHEADER));

  if (readBytes < sizeof(long long)) {
    printf("ERROR: File %s has only %zd bytes\n", filename, readBytes);
    exit(1);
  }

  header = *((long long *) pFileStarHeader);
  headerFlag = header & 0xffffffff;
  headerVersion = header >> 32;
  writeStarImageCount = 0;

  if (headerFlag == MAG_HEADER_DATA && readBytes == sizeof(FILESTARHEADER)) {
    for (index = 0; index < MAG_FILE_MODULUS; index++) {
      writeStarImageCount += pFileStarHeader->starCount[index];
    }
  }

  if (
    headerFlag != MAG_HEADER_DATA ||
    headerVersion != MAG_HEADER_VERSION ||
    readBytes != sizeof(FILESTARHEADER) ||
    pFileStarHeader->headersize != sizeof(FILESTARHEADER) ||
    pFileStarHeader->starimagesize != sizeof(FILESTARIMAGE) ||
    pFileStarHeader->curFileSize != filesize ||
    repairHandle != NULL ||
    writeStarImageCount != (filesize - sizeof(FILESTARHEADER)) / sizeof(FILESTARIMAGE)
  ) {
    /* We have a stale file.  We need to reformat the file */
    /* First read the entire file in */
#ifdef READONLY_PHOTOMETRY
    printf("ERROR: stale file %s\n",filename);
    exit(1);
#endif /* READONLY_PHOTOMETRY */

    if (headerFlag == MAG_FORMAT_DATA) {
      pFileCommon->counterBlock.noHeaderCount++;
      memcpy(pFileCommon->magnitudeBuffer, pFileStarHeader, readBytes);

      if ((filesize - readBytes) > 0) {
        readBytes2 = read(fileNumber, &((char *)pFileCommon->magnitudeBuffer)[readBytes], filesize - readBytes);
        if (readBytes2 != (filesize - readBytes)) {
          printf(
            "ERROR: %d %s read %zd (1) instead of %ld bytes from %s\n",
            errno,
            strerror(errno),
            readBytes2,
            filesize - readBytes,
            filename
          );
          exit(1);
        }
      }

      /* file has no header.  Move this data to where it belongs in the magnitude buffer */
      if (((filesize % sizeof(FILESTARIMAGE)) != 0) || (headerVersion == MAG_HEADER_VERSION_FIVE)) {
        fprintf(stderr, "FATAL: removed ReformatV5Records code\n");
        exit(1);
      }

      starImageCount = filesize / sizeof(FILESTARIMAGE);
    } else if (headerFlag == MAG_HEADER_DATA) {
      if (headerVersion == MAG_HEADER_VERSION_ONE && readBytes == sizeof(FILESTARHEADER)) {
        /* Version one has a corrupt fileOffset table - we will write a new header */
        printf("ERROR: file format version 1 is no longer supported for %s\n",filename);
        exit(1);
      } else if (headerVersion == MAG_HEADER_VERSION_TWO)  {
        printf("ERROR: file format version 2 is no longer supported for %s\n",filename);
        exit(1);
      } else if (headerVersion == MAG_HEADER_VERSION_THREE)  {
        printf("ERROR: file format version 3 is no longer supported for %s\n",filename);
        exit(1);
      } else if (headerVersion == MAG_HEADER_VERSION_FOUR)  {
        printf("ERROR: file format version 4 is no longer supported for %s\n",filename);
        exit(1);
      } else if (headerVersion == MAG_HEADER_VERSION_FIVE) {
        printf("ERROR: file format version 5 is no longer supported for %s\n",filename);
        exit(1);
      } else if (
        headerVersion == MAG_HEADER_VERSION &&
        readBytes == sizeof(FILESTARHEADER) &&
        pFileStarHeader->headersize == sizeof(FILESTARHEADER) &&
        pFileStarHeader->starimagesize == sizeof(FILESTARIMAGE)
      ) {
        /* This is a good file, but data has been tacked onto the end of it */
        if ((filesize - readBytes) < (writeStarImageCount * sizeof(FILESTARIMAGE))) {
          pFileCommon->counterBlock.corruptStarCount++;
        } else {
          pFileCommon->counterBlock.newDataCount++;
        }

        if ((filesize - readBytes) > 0) {
          if (((filesize - readBytes) % sizeof(FILESTARIMAGE)) != 0) {
            printf("ERROR: file %s size %ld is not divisible by FILESTARIMAGE\n", filename, filesize - readBytes);
            if (repairHandle == NULL) {
              exit(1);
            }
          }

          readBytes2 = read(fileNumber,pFileCommon->magnitudeBuffer,filesize-readBytes);
          if (readBytes2 != (filesize-readBytes)) {
            printf("ERROR: %d %s read %zd (5) instead of %ld bytes from %s\n",errno,strerror(errno),readBytes2,filesize-readBytes,filename);
            exit(1);
          }
        }

        starImageCount = (filesize - readBytes) / sizeof(FILESTARIMAGE);
      } else {
        printf("ERROR: file %s has a corrupt header version %lld or header size %zd\n", filename, headerVersion, readBytes);
        exit(1);
      }
    } else {
      printf("ERROR: file %s has a corrupt header flag %08x \n", filename, (unsigned int) headerFlag);
      exit(1);
    }

    if (starImageCount <= 0) {
      if (readOnly == 0) {
        printf("ERROR: no star images found in %s\n", filename);
        exit(1);
      } else {
        pFileCommon->curMagnitudes = 0;
        return 0;
      }
    }

    /* Now we need to perform versionId checking */

    pCurStarImage = pFileStarTable;
    pFileStarImage = pFileStarTable;
    readStarImageCount = 0;
    writeStarImageCount = 0;
    pFileCommon->counterBlock.readMagnitudes += starImageCount;

    while (readStarImageCount < starImageCount) {
      long long tmp_REFNumber = pFileStarImage->REFNumber;
      double tmp_magcal_local = pFileStarImage->magcal_local;
      double tmp_ra = pFileStarImage->ra;
      double tmp_dec = pFileStarImage->dec;
      int tmp_versionId = pFileStarImage->versionId;
      int tmp_plateNumber = pFileStarImage->plateNumber;
      int tmp2_gsc_bin_index = pFileStarImage->gsc_bin_index;

      header = *((long long *) pFileStarImage);
      headerFlag = header & 0xffffffff;
      headerVersion = header >> 32;
      versionId = GetLatestVersionId(
        pFileCommon,
        pFileStarImage->seriesId,
        tmp_plateNumber,
        tmp_versionId,
        pFileStarImage->REFNumber,
        0
      );

      if (
        (headerFlag != MAG_FORMAT_DATA) ||
        (headerVersion < 0) ||
        (GetREF(tmp_REFNumber,REF,0,0) != 0) ||
        (tmp_plateNumber > 100000) ||
        (headerVersion > MAG_FORMAT_VERSION) ||
        (versionId < 0) ||
        (tmp_magcal_local > 1000.0)
      ) {
        printf(
          "ERROR: record %d after header has flag 0x%x 0x%x and version %lld"
          " magcal_local %f and plateNumber %d versionId %d for %s REFNumber %lld\n",
          readStarImageCount,
          (unsigned int) headerFlag,
          (unsigned int) headerFlag,
          headerVersion,
          tmp_magcal_local,
          tmp_plateNumber,
          versionId,
          filename,
          tmp_REFNumber
        );

        if (repairHandle == NULL) {
          exit(1);
        } else {
          WriteStarbaseRecord(pFileStarImage, pPhotPlates, repairHandle, NULL, &writeHeader, catalogNumber);
        }
      } else {
        if (versionId < pFileCommon->minVersionId) {
          versionId = pFileCommon->minVersionId;
        }

        if ((RELEASE_EXPERIMENTAL != 0) && (catalogNumber == CATALOG_EXPERIMENTAL)) {
          /* Accept this star */
          AddCatalogReference(pFileCommon,pFileStarImage);

          if (readStarImageCount != writeStarImageCount) {
            memcpy(pCurStarImage,pFileStarImage,sizeof(FILESTARIMAGE));
          }

          pCurStarImage++;
          writeStarImageCount++;
        } else  {
          if ((repairHandle == NULL) && (versionId > tmp_versionId)) {
            /* Stale version here */
            pFileCommon->counterBlock.staleVersionCount++;
          } else {
            if (tmp2_gsc_bin_index == 0) {
              tmp_gsc_bin_index = GetGSCBin(pGscBin,tmp_ra,tmp_dec,&decBin,&raBin,"binCheck");
              if (tmp_gsc_bin_index != tmp2_gsc_bin_index) {
                pFileCommon->counterBlock.binZeroErrors++;
              } else {
                /* Accept this star */
                AddCatalogReference(pFileCommon,pFileStarImage);
                if (readStarImageCount != writeStarImageCount) {
                  memcpy(pCurStarImage,pFileStarImage,sizeof(FILESTARIMAGE));
                }
                pCurStarImage++;
                writeStarImageCount++;

              }
            } else {
              /* Accept this star */
              AddCatalogReference(pFileCommon,pFileStarImage);
              if (readStarImageCount != writeStarImageCount) {
                memcpy(pCurStarImage,pFileStarImage,sizeof(FILESTARIMAGE));
              }
              pCurStarImage++;
              writeStarImageCount++;
            }
          }
        }
      }

      pFileStarImage++;
      readStarImageCount++;
    }

    /* When performing a search_none sweep, count only once ! */
    if (gsc_base_index == pFileCommon->currentMagnitudesBin) {
      pFileCommon->counterBlock.currentMagnitudes += writeStarImageCount;
    }

    if (writeStarImageCount == 0) {
      /* Here the file has no valid entries.  Delete it */
      if (readOnly == 0) {
        pFileCommon->counterBlock.noValidImages++;
        if (fileNumber >= 0) {
          close(fileNumber);
          fileNumber = -1;
        }
        remove(filename);
      }

      pFileCommon->curMagnitudes = 0;
      return 0;
    } else {
      /* Sort the images in the order of increasing gsc_bin_index */
      qsort((void*)pFileStarTable, writeStarImageCount, sizeof(FILESTARIMAGE), StarImageCompare);

      /* Format a new header */
      memset(pFileStarHeader, 0, sizeof(FILESTARHEADER));
      pFileStarHeader->versionTag = MAG_HEADER_VERSION;
      pFileStarHeader->versionTag = pFileStarHeader->versionTag << 32;
      pFileStarHeader->versionTag = pFileStarHeader->versionTag | MAG_HEADER_DATA;
      pFileStarHeader->curFileSize = sizeof(FILESTARHEADER) + (writeStarImageCount * sizeof(FILESTARIMAGE));
      pFileStarHeader->gsc_base_index = gsc_base_index;
      pFileStarHeader->headersize = sizeof(FILESTARHEADER);
      pFileStarHeader->starimagesize = sizeof(FILESTARIMAGE);
      pCurStarImage = pFileStarTable;
      readStarImageCount = 0;
      curFileOffset =  sizeof(FILESTARHEADER);

      for (index = 0; index < MAG_FILE_MODULUS; index++) {
        int tmp2_gsc_bin_index;

        pFileStarHeader->starCount[index] = 0;
        pFileStarHeader->fileOffset[index] = curFileOffset;
        tmp2_gsc_bin_index = pCurStarImage->gsc_bin_index;

        if (readStarImageCount < writeStarImageCount && index + gsc_base_index == tmp2_gsc_bin_index) {
          pFileStarHeader->starCount[index]++;
          curFileOffset += sizeof(FILESTARIMAGE);

          while(readStarImageCount < writeStarImageCount) {
            readStarImageCount++;
            pCurStarImage++;
            tmp2_gsc_bin_index = pCurStarImage->gsc_bin_index;

            if (readStarImageCount >= writeStarImageCount || index + gsc_base_index != tmp2_gsc_bin_index) {
              break;
            }

            pFileStarHeader->starCount[index]++;
            curFileOffset += sizeof(FILESTARIMAGE);
          }
        }

        if (readStarImageCount != writeStarImageCount && index + gsc_base_index + 1 > tmp2_gsc_bin_index) {
          printf(
            "ERROR: software consistency check at photometryutils.c line %d %d %d %d\n",
            __LINE__,
            index,
            gsc_base_index,
            tmp2_gsc_bin_index
          );
          exit(1);
        }
      }

      /* Recheck our math */
      readStarImageCount = 0;
      for (index = 0; index < MAG_FILE_MODULUS; index++) {
        readStarImageCount += pFileStarHeader->starCount[index];
      }

      if (
        readStarImageCount != writeStarImageCount ||
        curFileOffset != sizeof(FILESTARHEADER) + writeStarImageCount * sizeof(FILESTARIMAGE)
      ) {
        printf(
          "ERROR: software consistency check for %s failed at photometryutils.c line %d for count %d %d or size %lu %ld\n",
          filename,
          __LINE__,
          readStarImageCount,
          writeStarImageCount,
          sizeof(FILESTARHEADER) + writeStarImageCount * sizeof(FILESTARIMAGE),
          curFileOffset
        );
        exit(1);
      }

      bufferFullFlag = 1; /* Do not read in the data again! */

      /* At this point write out the file and get rid of the old file */
      if (readOnly == 0) {
        char newfile[MAX_FILENAME + 4];
        int newFileNumber = -1;

        snprintf(newfile, sizeof(newfile), "%s.tmp", filename);
        newFileNumber = open(newfile, O_CREAT|O_WRONLY|O_TRUNC, 0664);
        if (newFileNumber < 0) {
          printf("ERROR creating file %s result %d errno %d %s\n", newfile, newFileNumber, errno, strerror(errno));
          exit(1);
        }

        writeBytes = write(newFileNumber, pFileStarHeader, sizeof(FILESTARHEADER));
        if (writeBytes != sizeof(FILESTARHEADER)) {
          printf("ERROR: Wrote only %zd bytes to %s\n", writeBytes, newfile);
          exit(1);
        }

        writeBytes = write(newFileNumber, pFileStarTable, writeStarImageCount * sizeof(FILESTARIMAGE));
        if (writeBytes != (writeStarImageCount * sizeof(FILESTARIMAGE))) {
          printf("ERROR: Wrote only %zd bytes to %s\n",writeBytes,newfile);
          exit(1);
        }

        close(newFileNumber);
        newFileNumber = -1;

        if (fileNumber >= 0) {
          close(fileNumber);
          fileNumber = -1;
        }

        if (!magfiles_use_new_layout) {
          ChangeGroup(newfile);
        }

        result = rename(newfile, filename);
        if (result != 0) {
          printf("ERROR: Exiting with result %d rename %s %s\n", result, newfile, filename);
          exit(1);
        }
      }
    }
  } else {
    /* Format is good. Just count the contents */
    pFileCommon->counterBlock.readMagnitudes += writeStarImageCount;
    if (gsc_base_index == pFileCommon->currentMagnitudesBin) {
      pFileCommon->counterBlock.currentMagnitudes += writeStarImageCount;
    }
  }

  if (headerHandle != NULL) {
    fprintf(headerHandle, "index\tgsc_bin_index\tstarCount\tfileOffset\n");
    fprintf(headerHandle, "-----\t-------------\t---------\t----------\n");

    for (index = 0; index < MAG_FILE_MODULUS; index++) {
      if (pFileStarHeader->starCount[index] != 0) {
        fprintf(
          headerHandle,
          "%d\t%d\t%d\t%ld\n",
          index,
          index + gsc_base_index,
          pFileStarHeader->starCount[index],
          pFileStarHeader->fileOffset[index]
        );
      }
    }
  }

  if (bufferFullFlag == 0) {
    /* We need to fill in our buffer.  Decide on our read limits */
    if (dumpAllFlag) {
      writeStarImageCount = (filesize - sizeof(FILESTARHEADER)) / sizeof(FILESTARIMAGE);
      readOffset = sizeof(FILESTARHEADER);
    } else {
      writeStarImageCount = 0;
      readOffset = -1;

      for (gsc_bin_index = gsc_bin_index_low; gsc_bin_index <= gsc_bin_index_high; gsc_bin_index++) {
        index = gsc_bin_index - gsc_base_index;

        if (index < 0 || index >= MAG_FILE_MODULUS) {
          continue;
        }

        if (pFileStarHeader->starCount[index]) {
          if (readOffset < 0) {
            readOffset = pFileStarHeader->fileOffset[index];
          }

          writeStarImageCount += pFileStarHeader->starCount[index];
        }
      }
    }

    if (writeStarImageCount > 0) {
      newOffset = lseek(fileNumber, readOffset, SEEK_SET);

      if (newOffset != readOffset || newOffset < 0) {
        printf(
          "ERROR: Seek error %d %s for %s to %ld\n",
          errno,
          strerror(errno),
          filename,
          readOffset
        );
        exit(1);
      }

      readBytes = read(fileNumber, pFileCommon->magnitudeBuffer, writeStarImageCount * sizeof(FILESTARIMAGE));
      if (readBytes != writeStarImageCount * sizeof(FILESTARIMAGE)) {
        printf(
          "ERROR: %d %s read %zd (8) instead of %ld bytes from %s\n",
          errno,
          strerror(errno),
          readBytes,
          writeStarImageCount * sizeof(FILESTARIMAGE),
          filename
        );
        exit(1);
      }
    }
  }

  if (fileNumber >= 0) {
    close(fileNumber);
    fileNumber = -1;
  }

  if (dumpAllFlag) {
    pFileCommon->counterBlock.returnedMagnitudes += writeStarImageCount;

    for (readStarImageCount = 0; readStarImageCount < writeStarImageCount; readStarImageCount++) {
      pFileStarImage = &((PFILESTARIMAGE)pFileCommon->magnitudeBuffer)[readStarImageCount];

      /* Because MAX_LIMITING_MAG was changed on Feb 15, 2013, this code enforces the new setting */
      if ((pFileStarImage->limiting_mag_local - pFileStarImage->magcal_magdep) < MAX_LIMITING_MAG) {
        pFileStarImage->AFLAGS |= (1 << FILTER_AFLAG_LIMITING_MAG);
      } else {
        pFileStarImage->AFLAGS &= (~(1 << FILTER_AFLAG_LIMITING_MAG));
      }
    }

    pFileCommon->curMagnitudes = writeStarImageCount;
    return writeStarImageCount;
  }

  /* Now we need to perform bin and flags checking */

  pCurStarImage = (PFILESTARIMAGE) pFileCommon->magnitudeBuffer;
  pFileStarImage = pCurStarImage;
  starImageCount = writeStarImageCount;
  readStarImageCount = 0;
  writeStarImageCount = 0;

  while (readStarImageCount < starImageCount) {
    if (pFileStarImage->gsc_bin_index < gsc_bin_index_low || pFileStarImage->gsc_bin_index > gsc_bin_index_high) {
      pFileCommon->counterBlock.binRejectCount++;
      pFileStarImage++;
      readStarImageCount++;
      continue;
    }

    if (pFileStarImage->versionId < pFileCommon->minVersionId) {
      pFileStarImage++;
      readStarImageCount++;
      pFileCommon->counterBlock.filteredVersionCount++;

      if (gsc_base_index == pFileCommon->currentMagnitudesBin) {
        pFileCommon->counterBlock.currentMagnitudes--;
      }
      continue;
    }

    if (
      catalogNumber == CATALOG_KEPLER ||
      catalogNumber == CATALOG_APASS ||
      catalogNumber == CATALOG_GAIA ||
      catalogNumber == CATALOG_ATLAS ||
      catalogNumber == CATALOG_EXPERIMENTAL
    ) {
      refType = GetREFType(pFileStarImage->REFNumber);
    }

    if (catalogNumber == CATALOG_KEPLER && (refType == REF_TYPE_GSC || refType == REF_TYPE_APASS)) {
      pFileStarImage++;
      readStarImageCount++;
      pFileCommon->counterBlock.mixedCatalogCount++;
      continue;
    } else if (catalogNumber == CATALOG_APASS && refType == REF_TYPE_KEPLER) {
      pFileStarImage++;
      readStarImageCount++;
      pFileCommon->counterBlock.mixedCatalogCount++;
      continue;
    } else if (catalogNumber == CATALOG_GAIA && refType != REF_TYPE_GAIA2 && refType != REF_TYPE_NONE) {
      pFileStarImage++;
      readStarImageCount++;
      pFileCommon->counterBlock.mixedCatalogCount++;
      continue;
    } else if (catalogNumber == CATALOG_ATLAS && refType != REF_TYPE_ATLAS2 && refType != REF_TYPE_NONE) {
      pFileStarImage++;
      readStarImageCount++;
      pFileCommon->counterBlock.mixedCatalogCount++;
      continue;
    } else if (catalogNumber == CATALOG_EXPERIMENTAL && (refType == REF_TYPE_GSC || refType == REF_TYPE_KEPLER)) {
      if (RELEASE_EXPERIMENTAL == 0 || catalogNumber != CATALOG_EXPERIMENTAL) {
        pFileStarImage++;
        readStarImageCount++;
        pFileCommon->counterBlock.mixedCatalogCount++;
        continue;
      }
    }

    if (pFileStarImage->gsc_bin_index == 0) {
      tmp_gsc_bin_index = GetGSCBin(pGscBin, pFileStarImage->ra, pFileStarImage->dec, &decBin, &raBin, "binCheck");

      if (tmp_gsc_bin_index != pFileStarImage->gsc_bin_index) {
        pFileCommon->counterBlock.binZeroErrors++;
        pFileStarImage++;
        readStarImageCount++;
        continue;
      }
    }

    if ((pFileStarImage->AFLAGS & (1 << FILTER_AFLAG_BLEND_NOMATCH)) != 0) {
      pFileCommon->counterBlock.nomatchRejectCount++;
      pFileStarImage++;
      readStarImageCount++;
      continue;
    }

    if (zeroAFLAGS != 0 && pFileStarImage->AFLAGS != 0) {
      pFileCommon->counterBlock.flagsRejectCount++;
      pFileStarImage++;
      readStarImageCount++;
      continue;
    }

    if (REFonly > 0 && pFileStarImage->REFNumber == 0) {
      pFileCommon->counterBlock.refRejectCount++;
      pFileStarImage++;
      readStarImageCount++;
      continue;
    }

    if (REFonly < 0 && pFileStarImage->REFNumber != 0) {
      pFileCommon->counterBlock.refRejectCount++;
      pFileStarImage++;
      readStarImageCount++;
      continue;
    }

    if (nosequestered != 0 && GetSequesteredFlag(pFileStarImage->seriesId) == SEQUESTERED_YES) {
      pFileStarImage++;
      pFileCommon->counterBlock.sequesteredCount++;
      readStarImageCount++;
    }

    /* Because MAX_LIMITING_MAG was changed on Feb 15, 2013, this code enforces the new setting */
    if ((pFileStarImage->limiting_mag_local - pFileStarImage->magcal_magdep) < MAX_LIMITING_MAG) {
      pFileStarImage->AFLAGS |= (1 << FILTER_AFLAG_LIMITING_MAG);
    } else {
      pFileStarImage->AFLAGS &= (~(1 << FILTER_AFLAG_LIMITING_MAG));
    }

    /* Keep track of the deepest limiting magnitude on this plate */
    if (
      (pFileStarImage->AFLAGS & FILTER_AMASK_LIMITING) == 0 &&
      pFileCommon->seriesLimitingMagnitudeArray != NULL &&
      pFileStarImage->limiting_mag_local < 90.0
    ) {
      float **seriesLimitingMagnitudeArray = pFileCommon->seriesLimitingMagnitudeArray;
      float *seriesLimitingMagnitudeTable = seriesLimitingMagnitudeArray[pFileStarImage->seriesId];

      if (seriesLimitingMagnitudeTable == NULL) {
        seriesLimitingMagnitudeTable = (float *) calloc(
          pFileCommon->maxPlateNumber[pFileStarImage->seriesId] + 1,
          sizeof(float)
        );

        if (seriesLimitingMagnitudeTable == NULL) {
          printf(
            "ERROR: failed to allocate seriesLimitingMagnitudeTable of size %d\n",
            pFileCommon->maxPlateNumber[pFileStarImage->seriesId] + 1
          );
          exit(1);
        }

        seriesLimitingMagnitudeArray[pFileStarImage->seriesId] = seriesLimitingMagnitudeTable;
      }

      if (pFileStarImage->limiting_mag_local > seriesLimitingMagnitudeTable[pFileStarImage->plateNumber]) {
        seriesLimitingMagnitudeTable[pFileStarImage->plateNumber] = pFileStarImage->limiting_mag_local;
      }
    }

    if (readStarImageCount != writeStarImageCount) {
      memcpy(pCurStarImage, pFileStarImage, sizeof(FILESTARIMAGE));
    }

    pCurStarImage++;
    writeStarImageCount++;
    pFileStarImage++;
    readStarImageCount++;
  }

  pFileCommon->counterBlock.returnedMagnitudes += writeStarImageCount;
  pFileCommon->curMagnitudes = writeStarImageCount;
  return writeStarImageCount;
}


/* Search the photplates database for a version and return it */
int
GetLatestVersionId(
  PFILECOMMON pFileCommon,
  unsigned char seriesId,
  int plateNumber,
  int versionId,
  long long REFNumber,
  int fatal
) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  char catalogString[5];
  int plateVersionId = 0;
  int *seriesVersionArray;
  if (seriesId > MAX_SERIES) {
    printf("ERROR: illegal seriesId for seriesId %d, plateNumber %d, REFNumber %lld\n",seriesId,plateNumber,REFNumber);
    if (fatal) {
      exit(1);
    } else {
      return(-1);
    }
  }

  if ((plateNumber < 0) || (plateNumber > pFileCommon->maxPlateNumber[seriesId])) {
    printf("ERROR: illegal plateNumber for seriesId %d, plateNumber %d, NUMBER %lld\n",seriesId,plateNumber,REFNumber);
    if (fatal) {
      exit(1);
    } else {
      return(-1);
    }

  }

#if 0
  if ((plateNumber == LOS_PLATENUMBER) &&
      (seriesId == LOS_SERIESID)) {
    printf("Line %d GetLatestVersionID() At seriesId %d plateNumber %d\n",__LINE__,plateNumber,seriesId);
  }
#endif


  seriesVersionArray = pFileCommon->seriesVersionTable[seriesId];
  if (seriesVersionArray == NULL) {
    seriesVersionArray = (int *)calloc(pFileCommon->maxPlateNumber[seriesId]+1,sizeof(int));
    if (seriesVersionArray == NULL) {
      printf("ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[seriesId]+1);
      if (fatal) {
        exit(1);
      } else {
        return(-1);
      }

    }
    pFileCommon->seriesVersionTable[seriesId] = seriesVersionArray;
  }
  if (seriesVersionArray[plateNumber] == 0) {
    /* Need to query the database to find the versionId */
    if (pFileCommon->catalogNumber == 0) {
      catalogString[0] = 0;
    } else {
      sprintf(catalogString,"%d",pFileCommon->catalogNumber);
    }
    plateVersionId = -1;
    sprintf(queryString,"SELECT versionId%s from photplates where seriesId = %d and plateNumber = %d;",catalogString,seriesId,plateNumber);
    res = ExecuteQuery(pFileCommon->pPhotConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pFileCommon->pPhotConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          if (plateVersionId > 0) {
            printf("ERROR: Multiple plate version Ids for plate %s%05d\n",GetSeriesString(seriesId,0),plateNumber);
          }
          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&plateVersionId);
            if (nvals != 1) {
              printf(
                "ERROR: nvals is %d for versionId%s of plate %s%05d\n",
                nvals,
                catalogString,
                GetSeriesString(seriesId, 0),
                plateNumber
              );
              plateVersionId = -2;
            }
          }
        }

        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: GetLatestVersionId failed to fetch version for %s%05d\n",GetSeriesString(seriesId,0),plateNumber);
      }
    }
    if ((plateVersionId <= 0) && (versionId > 0)) {
      plateVersionId = versionId;
    }
    if (plateVersionId < pFileCommon->minVersionId) {
      plateVersionId = pFileCommon->minVersionId;
    }

    seriesVersionArray[plateNumber] = plateVersionId;

  } else {
    /* Already have the version in cache */
    pFileCommon->counterBlock.versionIdHits++;
  }
  return(seriesVersionArray[plateNumber]);
}


/* Search the photplates database for the plate quality and return it.  The plate quality is a bitmask
   where bit <n> is set if the spatial bin is bad
*/
static int
GetPlateQualityMask(
  PFILECOMMON pFileCommon,
  unsigned char seriesId,
  int plateNumber,
  int *pQuality,
  double **pColorterm,
  int **pColorflag,
  long long REFNumber
) {
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  int plateQuality = 0;
  int spatial_bin;
  double colorterm;
  int colorflag;
  int *seriesQualityMaskArray;
  int *seriesQualityValueArray;
  double *seriesColortermArray;
  int *seriesColorflagArray;
  char catalogString[4];

  if (seriesId > MAX_SERIES) {
    printf("ERROR: illegal seriesId for seriesId %d, plateNumber %d, REFNumber %lld\n",seriesId,plateNumber,REFNumber);
    exit(1);
  }

  if ((plateNumber < 0) || (plateNumber > pFileCommon->maxPlateNumber[seriesId])) {
    printf("ERROR: illegal plateNumber for seriesId %d, plateNumber %d, NUMBER %lld\n",seriesId,plateNumber,REFNumber);
    exit(1);
  }

  if (pFileCommon->catalogNumber == 0) {
    catalogString[0] = '\0';
  } else {
    sprintf(catalogString, "%d", pFileCommon->catalogNumber);
  }

  *pQuality = 0;

  seriesQualityMaskArray = pFileCommon->seriesQualityMaskTable[seriesId];
  if (seriesQualityMaskArray == NULL) {
    seriesQualityMaskArray = (int *)calloc(pFileCommon->maxPlateNumber[seriesId]+1,sizeof(int));
    if (seriesQualityMaskArray == NULL) {
      printf("ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[seriesId]+1);
      exit(1);
    }
    pFileCommon->seriesQualityMaskTable[seriesId] = seriesQualityMaskArray;

    seriesQualityValueArray = (int *)calloc(pFileCommon->maxPlateNumber[seriesId]+1,sizeof(int));
    if (seriesQualityValueArray == NULL) {
      printf("ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[seriesId]+1);
      exit(1);
    }
    pFileCommon->seriesQualityValueTable[seriesId] = seriesQualityValueArray;

    seriesColortermArray = (double *)calloc(MAX_SPATIAL_BINS*(pFileCommon->maxPlateNumber[seriesId]+1),sizeof(double));
    if (seriesColortermArray == NULL) {
      printf("ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[seriesId]+1);
      exit(1);
    }
    pFileCommon->seriesColortermTable[seriesId] = seriesColortermArray;

    seriesColorflagArray = (int *)calloc(MAX_SPATIAL_BINS*(pFileCommon->maxPlateNumber[seriesId]+1),sizeof(int));
    if (seriesColorflagArray == NULL) {
      printf("ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[seriesId]+1);
      exit(1);
    }
    pFileCommon->seriesColorflagTable[seriesId] = seriesColorflagArray;

  }
  seriesQualityValueArray = pFileCommon->seriesQualityValueTable[seriesId];
  seriesColortermArray = pFileCommon->seriesColortermTable[seriesId];
  seriesColorflagArray = pFileCommon->seriesColorflagTable[seriesId];

  if (seriesQualityMaskArray[plateNumber] == 0) {
    /* Need to query the database to find the plate quality */
    sprintf(queryString,"SELECT quality+0 from photplates where seriesId = %d and plateNumber = %d;",seriesId,plateNumber);
    res = ExecuteQuery(pFileCommon->pPhotConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pFileCommon->pPhotConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&plateQuality);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for quality of plate %s%05d\n",nvals,GetSeriesString(seriesId,0),plateNumber);
              exit(1);
            }
          } else {
            plateQuality = 0;
          }
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: GetLatestVersionId failed to fetch version for %s%05d\n",GetSeriesString(seriesId,0),plateNumber);
      }
    }
    seriesQualityValueArray[plateNumber] = plateQuality;
    /* Now query the database to find the colorterm and colorflag for each spatial bin */
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      seriesColortermArray[(MAX_SPATIAL_BINS*plateNumber)+spatial_bin-1] = 99.0;
      seriesColorflagArray[(MAX_SPATIAL_BINS*plateNumber)+spatial_bin-1] = 0;
    }

    sprintf(
      queryString,
      "SELECT spatial_bin,colorterm,colorflag+0 from spatialbin%s where seriesId = %d and plateNumber = %d;",
      catalogString,
      seriesId,
      plateNumber
    );

    res = ExecuteQuery(pFileCommon->pPhotConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pFileCommon->pPhotConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&spatial_bin);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for spatial_bin of plate %s%05d\n",nvals,GetSeriesString(seriesId,0),plateNumber);
              exit(1);
            }

            if ((spatial_bin < 1) || (spatial_bin > MAX_SPATIAL_BINS)) {
              printf(
                "ERROR: illegal spatial_bin %d of plate %s%05d\n",
                spatial_bin,
                GetSeriesString(seriesId,0),
                plateNumber
              );
              exit(1);
            }
          } else {
            continue;
          }

          if (sqlrow[1]) {
            nvals = sscanf(sqlrow[1],"%lf",&colorterm);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for colorterm of plate %s%05d\n",nvals,GetSeriesString(seriesId,0),plateNumber);
              exit(1);
            } else {
              seriesColortermArray[(MAX_SPATIAL_BINS*plateNumber)+spatial_bin-1] = colorterm;
            }
          } else {
            continue;
          }
          if (sqlrow[2] != NULL) {
            nvals = sscanf(sqlrow[2],"%d",&colorflag);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for colorflag of plate %s%05d\n",nvals,GetSeriesString(seriesId,0),plateNumber);
              exit(1);
            } else {
              seriesColorflagArray[(MAX_SPATIAL_BINS*plateNumber)+spatial_bin-1] = colorflag;
            }
          } else {
            continue;
          }
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: GetLatestVersionId failed to fetch version for %s%05d\n",GetSeriesString(seriesId,0),plateNumber);
      }
    }

    if (plateQuality != 0) {
      /* Here the entire plate has been flagged because of a quality problem.  Set all spatial bins bad */
      seriesQualityMaskArray[plateNumber] = 0x7fffffff;
    } else {
      /* Here the overall plate quality is good.  Check the colorterm for each spatial bin */
      /* First mark all spatial bins bad */
      seriesQualityMaskArray[plateNumber] = 0x7fffffff;
      for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
        colorterm = seriesColortermArray[(MAX_SPATIAL_BINS*plateNumber)+spatial_bin-1];
        colorflag = seriesColorflagArray[(MAX_SPATIAL_BINS*plateNumber)+spatial_bin-1];
        if ((colorterm < 90.0) && (colorflag == COLORFLAG_METROPOLIS)) {
          switch(pFileCommon->catalogNumber) {
          case CATALOG_GSC232:
            if ((colorterm >= MIN_BLUE_COLORTERM) &&
                (colorterm <= MAX_BLUE_COLORTERM)) {
              /* This is a good spatial bin. Clear the appropriate bit */
              seriesQualityMaskArray[plateNumber] &= ~(1 << spatial_bin);

            }
            break;
          case CATALOG_KEPLER:
            if ((colorterm >= MIN_KEPLER_BLUE_COLORTERM) &&
                (colorterm <= MAX_KEPLER_BLUE_COLORTERM)) {
              /* This is a good spatial bin. Clear the appropriate bit */
              seriesQualityMaskArray[plateNumber] &= ~(1 << spatial_bin);

            }
            break;
          case CATALOG_APASS:
            if ((colorterm >= MIN_APASS_BLUE_COLORTERM) &&
                (colorterm <= MAX_APASS_BLUE_COLORTERM)) {
              /* This is a good spatial bin. Clear the appropriate bit */
              seriesQualityMaskArray[plateNumber] &= ~(1 << spatial_bin);

            }
            break;
          case CATALOG_ATLAS:
            if ((colorterm >= MIN_ATLAS_BLUE_COLORTERM) &&
                (colorterm <= MAX_ATLAS_BLUE_COLORTERM)) {
              /* This is a good spatial bin. Clear the appropriate bit */
              seriesQualityMaskArray[plateNumber] &= ~(1 << spatial_bin);

            }
            break;
          case CATALOG_GAIA:
            if ((colorterm >= MIN_GAIA_BLUE_COLORTERM) &&
                (colorterm <= MAX_GAIA_BLUE_COLORTERM)) {
              /* This is a good spatial bin. Clear the appropriate bit */
              seriesQualityMaskArray[plateNumber] &= ~(1 << spatial_bin);

            }
            break;
          case CATALOG_EXPERIMENTAL:
            if ((colorterm >= MIN_EXPERIMENTAL_BLUE_COLORTERM) &&
                (colorterm <= MAX_EXPERIMENTAL_BLUE_COLORTERM)) {
              /* This is a good spatial bin. Clear the appropriate bit */
              seriesQualityMaskArray[plateNumber] &= ~(1 << spatial_bin);

            }
            break;
          default:
            printf("ERROR: unknown catalog %d\n",pFileCommon->catalogNumber);
            exit(1);
            break;
          }
        }
      }
    }
  } else {
    /* Already have the version in cache */
    plateQuality = seriesQualityValueArray[plateNumber];
    pFileCommon->counterBlock.versionIdHits++;
  }
  *pQuality = plateQuality;
  *pColorterm = &seriesColortermArray[(MAX_SPATIAL_BINS*plateNumber)];
  *pColorflag = &seriesColorflagArray[(MAX_SPATIAL_BINS*plateNumber)];

  return(seriesQualityMaskArray[plateNumber]);
}


static int
SearchGSCBlends(
  double ra,
  double dec,
  PSKY2KSTAR pSky2KTable,
  int sky2k_size,
  double maxGSCRadius
) {
  PSKY2KSTAR pSky2KStar;
  double x;
  double y;
  double z;
  double cosx;
  double pt[3];
  double pos[3];
  struct kdres *presults;
  double distance;

  if (sky2k_size <= 0) {
    return(0);
  }

  z = sin(DEGREES_TO_RAD*dec);
  cosx = cos(DEGREES_TO_RAD*dec);
  x = cosx*cos(DEGREES_TO_RAD*ra);
  y = cosx*sin(DEGREES_TO_RAD*ra);
  pt[0] = x;
  pt[1] = y;
  pt[2] = z;
  presults = kd_nearest_range( ptree, pt, maxGSCRadius*DEGREES_TO_RAD );
  while( !kd_res_end( presults ) ) {
    /* get the data and position of the current result item */
    pSky2KStar = (PSKY2KSTAR)kd_res_item( presults, pos );
    distance = sqrt(sqr(pt[0]-pos[0])+sqr(pt[1]-pos[1])+sqr(pt[2]-pos[2]));
    if (distance < pSky2KStar->radius*DEGREES_TO_RAD) {
      kd_res_free( presults );
      return(1);
    }
    /* go to the next entry */
    kd_res_next( presults );
  }
  kd_res_free( presults );
  return(0);
}


// This loads a tranche of magnitudes from a GSC-bin64 bin. The name is
// misleading; it doesn't limit itself to "NONE" detections.
static int
LoadNoneImages(
  PGSCBIN pGscBin,
  PFILECOMMON pFileCommon,
  int gsc_bin_index,
  PPHOTSTARIMAGE* ppMagnitudeTable,
  char *catalogString,
  int readOnly,
  int verbose
) {
  int curMagnitudes = 0;
  int fileMagnitudes;
  int cache_index;
  int cache_index2;
  int min_cache_sequence;
  int magnitudeIndex;
  int baseGscBinIndex = (gsc_bin_index / MAG_FILE_MODULUS) * MAG_FILE_MODULUS;
  PFILESTARIMAGE pFileStarImage;
  PFILESTARIMAGE pCacheFileStarImage;
  PPHOTSTARIMAGE pCurStarImage;
  int zeroAFLAGS = 0;
  int REFonly = 0;
  int magnitudeLimit = 0;
  int dumpAllFlag = 0;
  int nosequestered = 1;

  pFileCommon->counterBlock.cache_access_count++;

  for (cache_index = 0; cache_index < CACHE_MAGNITUDE_ENTRIES; cache_index++) {
    if (pFileCommon->cache_gsc_bin_index[cache_index] == baseGscBinIndex) {
      /* Found the entry in the cache */
      pFileCommon->counterBlock.cache_hit_count++;
      curMagnitudes = pFileCommon->cache_bin_magnitudes[cache_index];
      break;
    }
  }

  if (cache_index == CACHE_MAGNITUDE_ENTRIES) {
    min_cache_sequence = pFileCommon->cache_current_sequence;
    /* Entry not in cache, find a new cache entry */

    for (cache_index2 = 0; cache_index2 < CACHE_MAGNITUDE_ENTRIES; cache_index2++) {
      if (pFileCommon->cache_gsc_bin_index[cache_index2] < 0) {
        cache_index = cache_index2;
        break;
      }

      if (pFileCommon->cache_access_sequence[cache_index2] <=  min_cache_sequence) {
        min_cache_sequence = pFileCommon->cache_access_sequence[cache_index2];
        cache_index = cache_index2;
      }
    }

    if (cache_index == CACHE_MAGNITUDE_ENTRIES) {
      printf("ERROR: LoadNoneImages cache sanity check\n");
      exit(1);
    }

    /* Now read the requested bin */
    fileMagnitudes = GetFileSummaryMagnitudes(
      pGscBin,
      pFileCommon,
      baseGscBinIndex,
      baseGscBinIndex + MAG_FILE_MODULUS - 1,
      zeroAFLAGS,
      REFonly,
      magnitudeLimit,
      dumpAllFlag,
      NULL,
      catalogString,
      verbose,
      NULL,
      readOnly,
      nosequestered
    );
    pFileCommon->counterBlock.cache_magnitude_read_count += fileMagnitudes;

    if (pFileCommon->startTime != 0) {
      time_t newTime;
      time(&newTime);
      newTime -= pFileCommon->startTime;

      if (newTime > WEB_TIMEOUT) {
        printf(
          "ERROR: Request processing timed out in line %d after %ld seconds\n",
          __LINE__,
          newTime
        );
        exit(1);
      }
    }

    pFileCommon->cache_bin_magnitudes[cache_index] = 0;

    if (pFileCommon->pCacheMagnitudeTable[cache_index] != NULL) {
      free(pFileCommon->pCacheMagnitudeTable[cache_index]);
      pFileCommon->pCacheMagnitudeTable[cache_index] = NULL;
    }

    if (pFileCommon->pCacheFileMagnitudeTable[cache_index] != NULL) {
      free(pFileCommon->pCacheFileMagnitudeTable[cache_index]);
      pFileCommon->pCacheFileMagnitudeTable[cache_index] = NULL;
    }

    pFileCommon->cache_gsc_bin_index[cache_index] = baseGscBinIndex;

    if (fileMagnitudes != 0) {
      pFileCommon->counterBlock.cache_miss_count++;
      /* Allocate tables to hold all of the star images */
      pFileCommon->pCacheMagnitudeTable[cache_index] = (PPHOTSTARIMAGE)calloc(fileMagnitudes,sizeof(PHOTSTARIMAGE));

      if (pFileCommon->pCacheMagnitudeTable[cache_index] == NULL) {
        printf("ERROR: failed to allocate %d PHOTSTARIMAGES in LoadNoneImages\n",fileMagnitudes);
        exit(1);
      }

      pFileCommon->pCacheFileMagnitudeTable[cache_index] = (PFILESTARIMAGE)calloc(fileMagnitudes,sizeof(FILESTARIMAGE));

      if (pFileCommon->pCacheFileMagnitudeTable[cache_index] == NULL) {
        printf("ERROR: failed to allocate %d FILESTARIMAGES in LoadNoneImages\n",fileMagnitudes);
        exit(1);
      }

      curMagnitudes = 0;

      for (magnitudeIndex = 0; magnitudeIndex < fileMagnitudes; magnitudeIndex++) {
        pCurStarImage = &pFileCommon->pCacheMagnitudeTable[cache_index][curMagnitudes];
        pCacheFileStarImage = &pFileCommon->pCacheFileMagnitudeTable[cache_index][curMagnitudes];
        pFileStarImage = &pFileCommon->magnitudeBuffer[magnitudeIndex];
        memset(pCurStarImage,0,sizeof(PHOTSTARIMAGE));
        memcpy(pCacheFileStarImage,pFileStarImage,sizeof(FILESTARIMAGE));
        pCurStarImage->pFileStarImage = pCacheFileStarImage;
        strcpy(pCurStarImage->series,GetSeriesString(pFileStarImage->seriesId,1));
        GetREF(pFileStarImage->REFNumber,pCurStarImage->REF,0,1);

        if (pFileStarImage->solutionNumber == 0) {
          sprintf(pCurStarImage->Plate, "%s%05d", pCurStarImage->series, pFileStarImage->plateNumber);
        } else {
          sprintf(pCurStarImage->Plate, "%s%05d_s%d", pCurStarImage->series, pFileStarImage->plateNumber, pFileStarImage->solutionNumber);
        }

        curMagnitudes++;
      }

      if (curMagnitudes == 0) {
        /* No good stars present after all */
        pFileCommon->cache_bin_magnitudes[cache_index] = 0;

        if (pFileCommon->pCacheMagnitudeTable[cache_index] != NULL) {
          free(pFileCommon->pCacheMagnitudeTable[cache_index]);
          pFileCommon->pCacheMagnitudeTable[cache_index] = NULL;
        }

        if (pFileCommon->pCacheFileMagnitudeTable[cache_index] != NULL) {
          free(pFileCommon->pCacheFileMagnitudeTable[cache_index]);
          pFileCommon->pCacheFileMagnitudeTable[cache_index] = NULL;
        }
      }

      pFileCommon->cache_bin_magnitudes[cache_index] = curMagnitudes;
      pFileCommon->cache_gsc_bin_index[cache_index] = baseGscBinIndex;
    }
  }

  if (pFileCommon->cache_current_sequence > (INT_MAX-5)) {
    printf("ERROR: cache_gsc_bin_index exceeds maximum reference in LoadNoneImages\n");
    exit(1);
  }

  pFileCommon->cache_current_sequence++;
  pFileCommon->cache_access_sequence[cache_index] = pFileCommon->cache_current_sequence;
  *ppMagnitudeTable = pFileCommon->pCacheMagnitudeTable[cache_index];
  return curMagnitudes;
}


// Load a tranche of magnitudes from a GSC-bin64 bin, *and* magnitudes from
// nearby bins if they're sufficiently close. The name is misleading; it doesn't
// limit itself to "NONE" detections.
void
LocateNoneImages(
  PGSCBIN pGscBin,
  PFILECOMMON pFileCommon,
  int gsc_bin_index,
  PPHOTSTARIMAGE* ppMagnitudeTable,
  int *pAllocMagnitudes,
  int *pCurMagnitudes,
  char *catalogString,
  int readOnly,
  int verbose
) {
  PPHOTSTARIMAGE pMagnitudeTable = *ppMagnitudeTable;
  PPHOTSTARIMAGE pBinMagnitudeTable;
  PPHOTSTARIMAGE tempMagnitudeTable = NULL;
  PPHOTSTARIMAGE pSourceImage;
  PFILESTARIMAGE pFileSourceImage;
  PPHOTSTARIMAGE pTargetImage;

  int curMagnitudes;
  int tmpMagnitudes;
  int finalMagnitudes = 0;
  int allocMagnitudes = *pAllocMagnitudes;
  int gsc_bin_list[MAX_ADJACENT_BINS];
  int gsc_bin_count;
  int gsc_index;
  int star_index;

  PBININDEX pBinIndex;
  double deltaRa;
  double minRa;
  double maxRa;
  double minDec;
  double maxDec;
  double tempRa;
  int raBin;
  int decBin;

  // See if there's anything for our target bin. If not, exit immediately
  curMagnitudes = LoadNoneImages(
    pGscBin,
    pFileCommon,
    gsc_bin_index,
    &pBinMagnitudeTable,
    catalogString,
    readOnly,
    verbose
  );

  if (curMagnitudes == 0) {
    *pCurMagnitudes = 0;
    return;
  }

  // Allocate and copy data.

  if (curMagnitudes > allocMagnitudes) {
    allocMagnitudes = 2 * (curMagnitudes + 1000);

    if (pMagnitudeTable == NULL) {
      pMagnitudeTable = (PPHOTSTARIMAGE) calloc(allocMagnitudes, sizeof(PHOTSTARIMAGE));
      if (pMagnitudeTable == NULL) {
        printf("ERROR: failed to allocate pMagnitudeTable of size %d\n",allocMagnitudes);
        exit(1);
      }
    } else {
      tempMagnitudeTable = realloc(pMagnitudeTable, allocMagnitudes * sizeof(PHOTSTARIMAGE));
      if (tempMagnitudeTable == NULL) {
        printf("ERROR: failed to allocate tempMagnitudeTable of size %d\n",allocMagnitudes);
        exit(1);
      }

      pMagnitudeTable = tempMagnitudeTable;
      tempMagnitudeTable = NULL;
    }
  }

  for (star_index = 0; star_index < curMagnitudes; star_index++) {
    pSourceImage = &pBinMagnitudeTable[star_index];

    if (pSourceImage->pFileStarImage->gsc_bin_index == gsc_bin_index) {
      pTargetImage = &pMagnitudeTable[finalMagnitudes];
      memcpy(pTargetImage, pSourceImage, sizeof(PHOTSTARIMAGE));
      finalMagnitudes++;
    }
  }

  // Get current bin limits and expand; BIN_EXPANSION_RADIUS is currently 1/64 deg.

  pBinIndex = GetSubbins(pGscBin, gsc_bin_index, &raBin, &decBin, (char *) "LocateNoneImages");
  deltaRa = 360.00 / pBinIndex->numBins;
  minRa = (360.00 * (1.0 * raBin)) / pBinIndex->numBins;
  maxRa = minRa + deltaRa;
  minDec = pBinIndex->declination;
  maxDec = minDec + pGscBin->bin_size;

  minDec -= BIN_EXPANSION_RADIUS;
  maxDec += BIN_EXPANSION_RADIUS;
  minRa -= BIN_EXPANSION_RADIUS * deltaRa / pGscBin->bin_size;
  maxRa += BIN_EXPANSION_RADIUS * deltaRa / pGscBin->bin_size;

  // For each adjacent bin, load data ...

  FindAdjacentBins(pGscBin, gsc_bin_index, gsc_bin_list, &gsc_bin_count);

  for (gsc_index = 0; gsc_index < gsc_bin_count; gsc_index++) {
    tmpMagnitudes = LoadNoneImages(
      pGscBin,
      pFileCommon,
      gsc_bin_list[gsc_index],
      &pBinMagnitudeTable,
      catalogString,
      readOnly,
      verbose
    );

    for (star_index = 0; star_index < tmpMagnitudes; star_index++) {
      pSourceImage = &pBinMagnitudeTable[star_index];
      pFileSourceImage = pSourceImage->pFileStarImage;

      if (pSourceImage->pFileStarImage->gsc_bin_index != gsc_bin_list[gsc_index]) {
        continue;
      }

      // Only apply spatial filtering if it is a NONE entry -- ???

      if (pFileSourceImage->REFNumber == 0) {
        if (pFileSourceImage->dec < minDec || pFileSourceImage->dec > maxDec) {
          continue;
        }

        tempRa = pFileSourceImage->ra;

        if ((tempRa - minRa) > 180.0) {
          tempRa -= 360.0;
        } else if ((minRa - tempRa) > 180.0) {
          tempRa += 360.0;
        }

        if (tempRa < minRa || tempRa > maxRa) {
          continue;
        }
      }

      // Alloc and copy

      if (finalMagnitudes + 10 > allocMagnitudes) {
        allocMagnitudes = 2 * (finalMagnitudes + 1000);

        if (pMagnitudeTable == NULL) {
          pMagnitudeTable = (PPHOTSTARIMAGE)calloc(allocMagnitudes,sizeof(PHOTSTARIMAGE));
          if (pMagnitudeTable == NULL) {
            printf("ERROR: failed to allocate pMagnitudeTable of size %d\n",allocMagnitudes);
            exit(1);
          }
        } else {
          tempMagnitudeTable = realloc(pMagnitudeTable,allocMagnitudes*sizeof(PHOTSTARIMAGE));
          if (tempMagnitudeTable == NULL) {
            printf("ERROR: failed to allocate tempMagnitudeTable of size %d\n",allocMagnitudes);
            exit(1);
          }

          pMagnitudeTable = tempMagnitudeTable;
          tempMagnitudeTable = NULL;
        }
      }

      pTargetImage = &pMagnitudeTable[finalMagnitudes];
      memcpy(pTargetImage, pSourceImage, sizeof(PHOTSTARIMAGE));
      finalMagnitudes++;
    }
  }

  *ppMagnitudeTable = pMagnitudeTable;
  *pAllocMagnitudes = allocMagnitudes;
  *pCurMagnitudes = finalMagnitudes;
}


void
AllocFileCommonVectors(PFILECOMMON pFileCommon,int vectorSize)
{
  if (pFileCommon->vectorAllocCount < vectorSize) {
    if (pFileCommon->vector1 != NULL) {
      free(pFileCommon->vector1);
    }
    if (pFileCommon->vector2 != NULL) {
      free(pFileCommon->vector2);
    }
    if (pFileCommon->vector3 != NULL) {
      free(pFileCommon->vector3);
    }
    pFileCommon->vectorAllocCount = vectorSize+STATS_REALLOC_INCREMENT;
    pFileCommon->vector1 = (double *)calloc(pFileCommon->vectorAllocCount,sizeof(double));
    if (pFileCommon->vector1 == NULL) {
      printf("ERROR: Failed to allocate pFileCommon->vector1 of size %d\n",pFileCommon->vectorAllocCount);
      exit(1);
    }
    pFileCommon->vector2 = (double *)calloc(pFileCommon->vectorAllocCount,sizeof(double));
    if (pFileCommon->vector2 == NULL) {
      printf("ERROR: Failed to allocate pFileCommon->vector2 of size %d\n",pFileCommon->vectorAllocCount);
      exit(1);
    }
    pFileCommon->vector3 = (double *)calloc(pFileCommon->vectorAllocCount,sizeof(double));
    if (pFileCommon->vector3 == NULL) {
      printf("ERROR: Failed to allocate pFileCommon->vector3 of size %d\n",pFileCommon->vectorAllocCount);
      exit(1);
    }
  }
}


static void
ProcessProperMotions(PFILECOMMON pFileCommon,PNEARBYCATALOG pNearbyCatalog,int curGscBinIndex)
{
  PPMSTATS pPmStats;
  int properMotionIndex;
  int curCount = 0;
  PGSCIMAGE pCurGscImageX;
  double drad;
  int vectorIndex = 0;
  double ramedian;
  double rarms;
  double decmedian;
  double decrms;
  double rapmmedian;
  double rapmrms;
  double decpmmedian;
  double decpmrms;
  double raProperMotion = 0;
  double decProperMotion = 0;
  double minEpoch = 2.0*GSC_EQUINOX;
  double maxEpoch = 0.0;
  double factor;
  int count;

  /* The following change was made on July 17, 2015, but may not be fully correct */
  if (pFileCommon->enableTransientSearch != 0) {
    if (pFileCommon->catalogFD < 0) {
      InitCatalogAccess(pFileCommon);
    }
  }

  pCurGscImageX = GetCatalogRecord(pFileCommon,pNearbyCatalog->REFNumber,pNearbyCatalog->gsc_bin_index);
#if 0
  if (pNearbyCatalog->REFNumber == LOS_REF_NUMBER ) {
    printf("ProcessProperMotions count %d gsc_bin_index %d cur bin %d for REFNumber %lld",pNearbyCatalog->properMotionCount,pNearbyCatalog->gsc_bin_index,curGscBinIndex,pNearbyCatalog->REFNumber);
    if (pCurGscImageX != NULL) {
      printf(" catalog ra %.5f catalog dec %.5f RaPM %.2f DecPM %.2f",pCurGscImageX->ra,pCurGscImageX->dec,pCurGscImageX->RaPM,pCurGscImageX->DecPM);
    } else {
      printf(" ERROR: no catalog entry!\n");
    }
    printf("\n");
  }
#endif
  if (pNearbyCatalog->properMotionCount == 0) {
#if 1
    if (pCurGscImageX != NULL) {
      pNearbyCatalog->ra = pCurGscImageX->ra;
      pNearbyCatalog->dec = pCurGscImageX->dec;
    }
#endif
    return;
  }
  AllocFileCommonVectors(pFileCommon,pNearbyCatalog->properMotionCount);

  properMotionIndex = pNearbyCatalog->firstProperMotionIndex;

  while (properMotionIndex >= 0) {
    pPmStats = &pFileCommon->pm_stats_table[properMotionIndex];
    curCount++;
    if (pPmStats->epoch < minEpoch) {
      minEpoch = pPmStats->epoch;
    }
    if (pPmStats->epoch > maxEpoch) {
      maxEpoch = pPmStats->epoch;
    }
    pFileCommon->vector1[vectorIndex] = pPmStats->catra;
    pFileCommon->vector2[vectorIndex] = pPmStats->catdec;
    raProperMotion += pPmStats->RaPM;
    decProperMotion += pPmStats->DecPM;
    vectorIndex++;


    if (pCurGscImageX != NULL) {
      drad = fabs(pCurGscImageX->ra-pPmStats->catra);
      if (pFileCommon->enableTransientSearch == 0) {
        if (drad > 0.00001) {
          printf("ERROR in ra %.5f %.5f diff %f\n",pCurGscImageX->ra,pPmStats->catra,drad);
        }

        drad = fabs(pCurGscImageX->dec-pPmStats->catdec);
        if (drad > 0.00001) {
          printf("ERROR in dec %.5f %.5f diff %f\n",pCurGscImageX->dec,pPmStats->catdec,drad);
        }

        drad = fabs(pCurGscImageX->RaPM-pPmStats->RaPM);
        if (drad > 0.001) {
          printf("ERROR in RaPM %.5f %.5f diff %f\n",pCurGscImageX->RaPM,pPmStats->RaPM,drad);
        }

        drad = fabs(pCurGscImageX->DecPM-pPmStats->DecPM);
        if (drad > 0.001) {
          printf("ERROR in DecPM %.5f %.5f diff %f\n",pCurGscImageX->DecPM,pPmStats->DecPM,drad);
        }
      }
    }
#if 0
    if (pNearbyCatalog->REFNumber == LOS_REF_NUMBER) {
      printf("ra %.5f dec %.5f catra %.5f catdec %.5f ra_2 %.5f dec_2 %.5f RaPM %.2f DecPM %.2f epoch %.7f count %5d link %5d\n",
             pPmStats->ra,
             pPmStats->dec,
             pPmStats->catra,
             pPmStats->catdec,
             pPmStats->ra_2,
             pPmStats->dec_2,
             pPmStats->RaPM,
             pPmStats->DecPM,
             pPmStats->epoch,
             pPmStats->count,
             pPmStats->link);
    }
#endif
    properMotionIndex = pPmStats->link;
  }
  if (curCount != pNearbyCatalog->properMotionCount) {
    printf("ERROR: curCount %d expected %d\n",curCount,pNearbyCatalog->properMotionCount);
  } else {
#if 0
    printf("curCount %d expected %d\n",curCount,pNearbyCatalog->properMotionCount);
#endif
  }
  if (CalcMedianAndRMS(vectorIndex,0,pFileCommon->vector1,&ramedian,&rarms,0,3.0,0) <= 0) {
    rarms = 99.0;
  }
  if (CalcMedianAndRMS(vectorIndex,0,pFileCommon->vector2,&decmedian,&decrms,0,3.0,0) <= 0) {
    decrms = 99.0;
    decmedian = 0.0;
  }
  if (vectorIndex == 1) {
    decrms = 0.0;
    rarms = 0.0;
  }
#if 0
  printf("ramedian %.5f rarms %.5f decmedian %.5f decrms %.5f\n",ramedian,rarms,decmedian,decrms);
#endif
  /* If the rms values are good, substitute estimated catalog magnitudes with correct magnitudes */
  factor = cos(decmedian * DEGREES_TO_RAD);
  rarms = rarms*factor;
  if (rarms > pFileCommon->counterBlock.maxPMrarms) {
    pFileCommon->counterBlock.maxPMrarms = rarms;
  }
  if (decrms > pFileCommon->counterBlock.maxPMdecrms) {
    pFileCommon->counterBlock.maxPMdecrms = decrms;
  }

  if ((rarms < PMSTATS_MAX_RMS_ERROR) && (decrms < PMSTATS_MAX_RMS_ERROR)) {
    pNearbyCatalog->RaPM = raProperMotion/(1.0*vectorIndex);
    pNearbyCatalog->DecPM = decProperMotion/(1.0*vectorIndex);

    if (pNearbyCatalog->dradRMS2 >= 0) {
      pNearbyCatalog->ra = ramedian;
      pNearbyCatalog->dec = decmedian;

    }
    if (pNearbyCatalog->dradRMS2GoodMag >= 0) {
      pNearbyCatalog->raGoodMag = ramedian;
      pNearbyCatalog->decGoodMag = decmedian;
    }

  } else {
#if 0
    printf("ERROR: rarms %f or decrms %f too high gsc_bin_index %d for %lld\n",rarms,decrms,curGscBinIndex,pNearbyCatalog->REFNumber);
#endif
    pFileCommon->counterBlock.properMotionErrorCount++;
    pNearbyCatalog->dradRMS2 = -1;
    pNearbyCatalog->dradRMS2GoodMag = -1;
    pNearbyCatalog->ra = 999.0;
    pNearbyCatalog->dec = 99.0;
    pNearbyCatalog->raGoodMag = 999.0;
    pNearbyCatalog->decGoodMag = 99.0;

  }

  if ((pFileCommon->propermotionHandle != NULL) &&
      (pNearbyCatalog->gsc_bin_index == curGscBinIndex) &&
      (pNearbyCatalog->dec <= 90.0) &&
      ((maxEpoch-minEpoch) > PROPER_MOTION_EPOCH_SPAN) &&
      ((maxEpoch-minEpoch)  < 2000.0) &&
      (pNearbyCatalog->properMotionCount > MIN_PROPER_MOTION_COUNT)) {

    /* Here we want to calculate the observed RaPM and DecPM values */
    curCount = 0;
    vectorIndex = 0;
    properMotionIndex = pNearbyCatalog->firstProperMotionIndex;

    while (properMotionIndex >= 0) {
      pPmStats = &pFileCommon->pm_stats_table[properMotionIndex];

      if (((pPmStats->AFLAGS & FILTER_AMASK_LOWDRAD) == 0) &&
          ((pPmStats->BFLAGS & (1 << FILTER_BFLAG_PSFSATURATED)) == 0) &&
          (pPmStats->ELLIPTICITY < MAX_ELLIPTICITY)) {

        pFileCommon->vector1[vectorIndex] = pPmStats->epoch;
        pFileCommon->vector2[vectorIndex] = pPmStats->ra;
        pFileCommon->vector3[vectorIndex] = pPmStats->dec;
        if (vectorIndex == 0) {
          minEpoch = pPmStats->epoch;
          maxEpoch = pPmStats->epoch;
        } else {
          if (pPmStats->epoch < minEpoch) {
            minEpoch = pPmStats->epoch;
          }
          if (pPmStats->epoch > maxEpoch) {
            maxEpoch = pPmStats->epoch;
          }

        }

        curCount++;
        vectorIndex++;
      }

      properMotionIndex = pPmStats->link;
    }

    count = curCount;
    if ((curCount > MIN_PROPER_MOTION_COUNT) &&
        ((maxEpoch-minEpoch) > PROPER_MOTION_EPOCH_SPAN)) {
      double offset;
      double slope;
      double mean;
      double std;
      double residual;
      pPmStats = &pFileCommon->pm_stats_table[0];
      if ((pPmStats->catdec < 90.0) || (pPmStats->catdec > -90.0)) {
        factor = cos(pPmStats->catdec * DEGREES_TO_RAD);
      } else {
        factor = 1.0;
      }
      LinearFit(curCount,pFileCommon->vector1,pFileCommon->vector2,&slope,&offset,&mean,&std,&residual);
      rapmmedian = slope*3600.0*1000.0*factor;
      rapmrms = std*3600.0*1000.0;
      LinearFit(curCount,pFileCommon->vector1,pFileCommon->vector3,&slope,&offset,&mean,&std,&residual);
      decpmmedian = slope*3600.0*1000.0;
      decpmrms = std*3600.0*1000.0;

    } else {
      count = 0;
    }

    if (count > MIN_PROPER_MOTION_COUNT) {
      char REF[MAX_REF];
      GetREF(pNearbyCatalog->REFNumber,REF,1,0);

      fprintf(pFileCommon->propermotionHandle,"%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\n",
              REF,
              count,
              rapmmedian,
              rapmrms,
              decpmmedian,
              decpmrms,
              pNearbyCatalog->RaPM,
              pNearbyCatalog->DecPM);
    }
  }
}


static void
UpdateGroupRange(PGROUPSTATS pGroupStats,double ra, double dec)
{
  if (pGroupStats->maxRA < ra) {
    pGroupStats->maxRA = ra;
  }
  if (pGroupStats->maxDec < dec) {
    pGroupStats->maxDec = dec;
  }
  if (pGroupStats->minRA > ra) {
    pGroupStats->minRA = ra;
  }
  if (pGroupStats->minDec > dec) {
    pGroupStats->minDec = dec;
  }
  return;
}


static void
SaveGroupStatistics(
  PFILECOMMON pFileCommon,
  PGROUPSTATS pGroupStats,
  PFILESTARIMAGE pNoneMagnitudeTable,
  int magnitudeCount
) {
  double factor  = cos(DEGREES_TO_RAD*pGroupStats->dec);
  double dra;
  double ddec;
  double drad;
  double maxDrad = 0;
  int histogramIndex;

  int groupVectorIndex;
  int groupMemberIndex;
  PFILESTARIMAGE pGroupStarImage;
  PNONESTATS pGroupNoneStats;
  groupMemberIndex = pGroupStats->flink;
  for (groupVectorIndex = 0; groupVectorIndex < pGroupStats->groupCount; groupVectorIndex++) {
    pGroupStarImage =  &pNoneMagnitudeTable[groupMemberIndex];
    dra = (pGroupStarImage->ra - pGroupStats->ra)*factor;
    ddec = pGroupStarImage->dec - pGroupStats->dec;
    drad = 3600.0*sqrt(sqr(dra)+sqr(ddec)); /* Distance in arcsec */
    if (drad > maxDrad) {
      maxDrad = drad;
    }
    pGroupNoneStats = &pFileCommon->none_stats_table[groupMemberIndex];
    groupMemberIndex = pGroupNoneStats->flink;
    if (groupVectorIndex == pGroupStats->groupCount-1) {
      if (groupMemberIndex > 0) {
        printf("ERROR: Illegal end groupMemberIndex %d in line %d\n",groupMemberIndex,__LINE__);
        exit(1);
      }
    } else {
      if ((groupMemberIndex < 0) ||	(groupMemberIndex >= magnitudeCount)) {
        printf("ERROR: Illegal groupMemberIndex %d %d in line %d\n",groupMemberIndex,magnitudeCount,__LINE__);
        exit(1);
      }
    }

  }
  if (pFileCommon->counterBlock.maxGroupDrad < maxDrad) {
    pFileCommon->counterBlock.maxGroupDrad = maxDrad;
  }
  if (pGroupStats->groupCount > 1) {
    histogramIndex = maxDrad;
    if (histogramIndex > MAX_GROUP_DRAD_HISTOGRAM) {
      histogramIndex = MAX_GROUP_DRAD_HISTOGRAM;
    }
    pFileCommon->counterBlock.maxGroupHistogram[histogramIndex]++;
  }
  return;
}


/* Returns 1 if the specified plate from an unmatched object matches a matched object with the given REFNumber */
static int
CheckDuplicatePlate(
  PFILESTARIMAGE pNoneMagnitudeTable,
  int magnitudeCount,
  long long REFNumber,
  unsigned char seriesId,
  int plateNumber,
  int solutionNumber
) {
  int magnitudeIndex;
  PFILESTARIMAGE pFileStarImage;
  for (magnitudeIndex = 0; magnitudeIndex < magnitudeCount; magnitudeIndex++ ) {
    pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
    if ((pFileStarImage->REFNumber == REFNumber) &&
        (pFileStarImage->seriesId == seriesId) &&
        (pFileStarImage->plateNumber == plateNumber) &&
        (pFileStarImage->solutionNumber == solutionNumber)) {
      return(1);
    }
  }
  return(0);
}


static void
InsertBlendChain(
  PFILECOMMON pFileCommon,
  int magnitudeIndex1,
  int magnitudeIndex2
) {
  PNONESTATS pNoneStats1 = &pFileCommon->none_stats_table[magnitudeIndex1];
  PNONESTATS pNoneStats2 = &pFileCommon->none_stats_table[magnitudeIndex2];
  PGROUPSTATS pGroupStats1;
  PGROUPSTATS pGroupStats2;
  int primaryGroupNumber;
  int secondaryGroupNumber;
  int savedFlink;
  int nextFlink;
  int curFlink;
  int groupIndex;
  if ((pNoneStats1->blendGroupNumber >= 0) && (pNoneStats2->blendGroupNumber >= 0)) {
    /* Both of these objects are currently in a blend group.  If the same blend group, we are done.  Otherwise, we must
       join the two groups */
    if (pNoneStats1->blendGroupNumber == pNoneStats2->blendGroupNumber) {
      return;
    } else {
      if (pNoneStats1->blendGroupNumber < pNoneStats2->blendGroupNumber) {
        primaryGroupNumber = pNoneStats1->blendGroupNumber;
        secondaryGroupNumber = pNoneStats2->blendGroupNumber;
      } else {
        primaryGroupNumber = pNoneStats2->blendGroupNumber;
        secondaryGroupNumber = pNoneStats1->blendGroupNumber;
      }
      pGroupStats1 = &pFileCommon->group_stats_table[primaryGroupNumber];
      pGroupStats2 = &pFileCommon->group_stats_table[secondaryGroupNumber];

      curFlink = pGroupStats2->blendFlink;
      while (curFlink >= 0) {
        pNoneStats2 = &pFileCommon->none_stats_table[curFlink];
        pNoneStats2->blendGroupNumber = primaryGroupNumber;
        nextFlink = pNoneStats2->blendFlink;
        savedFlink = pGroupStats1->blendFlink;
        pNoneStats2->blendFlink = savedFlink;
        pGroupStats1->blendFlink = curFlink;
        pGroupStats2->blendGroupMemberCount--;
        pGroupStats1->blendGroupMemberCount++;
        if (pGroupStats2->blendGroupMemberCount < 0) {
          printf("ERROR: blendGroupMemberCount is %d in line %d\n",pGroupStats2->blendGroupMemberCount,__LINE__);
          exit(1);
        }
        curFlink = nextFlink;
      }
      if (pGroupStats2->blendGroupMemberCount != 0) {
        printf("ERROR: blendGroupMemberCount is %d in line %d\n",pGroupStats2->blendGroupMemberCount,__LINE__);
        exit(1);
      }
      pGroupStats2->blendGroupNumber = -1;
      pGroupStats2->blendFlink = -1;
      return;
    }

  } else if ((pNoneStats1->blendGroupNumber >= 0) || (pNoneStats2->blendGroupNumber >= 0)) {
    /* Here we are adding a member to a group */
    if (pNoneStats1->blendGroupNumber >= 0) {
      primaryGroupNumber = pNoneStats1->blendGroupNumber;
      curFlink = magnitudeIndex2;
    } else {
      primaryGroupNumber = pNoneStats2->blendGroupNumber;
      curFlink = magnitudeIndex1;
    }
    /* First be sure that we are not already a menber of the group */
    pGroupStats1 = &pFileCommon->group_stats_table[primaryGroupNumber];
    nextFlink = pGroupStats1->blendFlink;
    while (nextFlink >= 0) {
      if (nextFlink == curFlink) {
        return;
      }
      pNoneStats2 = &pFileCommon->none_stats_table[nextFlink];
      nextFlink = pNoneStats2->blendFlink;
    }
    pNoneStats2 = &pFileCommon->none_stats_table[curFlink];
    savedFlink = pGroupStats1->blendFlink;
    pGroupStats1->blendFlink = curFlink;
    pGroupStats1->blendGroupMemberCount++;
    pNoneStats2->blendFlink = savedFlink;
    pNoneStats2->blendGroupNumber = primaryGroupNumber;
    return;
  } else {
    /* Allocate a new  group */
    for (groupIndex = 0; groupIndex < pFileCommon->blendGroupCount; groupIndex++) {
      pGroupStats1 = &pFileCommon->group_stats_table[groupIndex];
      if (pGroupStats1->blendGroupNumber < 0) {
        break;
      }
    }
    if (groupIndex >= pFileCommon->blendGroupCount) {
      groupIndex = pFileCommon->blendGroupCount;
      pFileCommon->blendGroupCount++;
    }
    pGroupStats1 = &pFileCommon->group_stats_table[groupIndex];
    pGroupStats1->blendGroupNumber = groupIndex;
    pGroupStats1->blendFlink = magnitudeIndex1;
    pGroupStats1->blendGroupMemberCount = 2;
    pNoneStats1->blendGroupNumber = groupIndex;
    pNoneStats1->blendFlink = magnitudeIndex2;
    pNoneStats2->blendGroupNumber = groupIndex;
    return;
  }
}


static void
GetAFLAGSStats(PFILECOMMON pFileCommon,int AFLAGS)
{
  int bit;
  int bitmask;
  /* Ignore plate-wide errors and variables condition */
  int tempAFLAGS =  ~((1 << FILTER_AFLAG_QUALITY) | (1 << GSC_VARIABLE_BIT) | (1 << FILTER_AFLAG_UNCERTAIN_DATE));
  tempAFLAGS = AFLAGS & tempAFLAGS;
  if (tempAFLAGS == 0) {
    pFileCommon->counterBlock.AFLAGS_ZERO_COUNT++;
    return;
  }
  pFileCommon->counterBlock.AFLAGS_NONZERO_COUNT++;
  for (bit = 0; bit < 32; bit++) {
    bitmask = 1 << bit;
    if (tempAFLAGS == bitmask) {
      pFileCommon->counterBlock.AFLAGS_UNIQUE_COUNT[bit]++;
    }
    if ((tempAFLAGS & bitmask) != 0) {
      pFileCommon->counterBlock.AFLAGS_MULT_COUNT[bit]++;
    }
  }
}


static void
ProcessBlendGroups(
  PGSCBIN pGscBin,
  PFILECOMMON pFileCommon,
  PFILESTARIMAGE pNoneMagnitudeTable,
  int magnitudeCount,
  PSKY2KSTAR pSky2KTable,
  int sky2k_size,
  double maxGSCRadius,
  int curGscBinIndex,
  int verbose
) {
  int groupIndex;
  int groupIndex2;
  int dupFound;
  PGROUPSTATS pGroupStatsTC;
  PGROUPSTATS pGroupStatsTC2;
  PGROUPSTATS pGroupStatsBlend;
  PGROUPSTATS pGroupStats3;
  PNONESTATS pNoneStats;
  PNONESTATS pNoneStats2;
  PNONESTATS pNoneStats3;
  int nextFlink;
  int nextFlink2;
  int groupMemberCount;
  int groupMemberCount2;
  PFILESTARIMAGE pCurStarImage;
  PFILESTARIMAGE pCurStarImage2;
  double factor;
  double testra;
  double testdec;
  double testdra;
  double testddec;
  double testdrad;
  int maxRefcnt = 0;
  int maxRefcntIndex = -1;
  int maxRefcntIndex2 = -1;
  int maxRefcnt2 = 0;
  int bestFlag = 1;
  int groupVectorIndex;
  int groupMemberIndex;
  PFILESTARIMAGE pGroupStarImage;
  PNONESTATS pGroupNoneStats;
  int raBin;
  int decBin;
  int curFlink;
  int *groupFlinkPtr;
  double tempEpoch;
  double distance1;
  double distance2;
  int badFlink;
  double raAverage;
  double decAverage;
  double positionWeight;
  long long newREFNumber;
  int positionCount;
  int result;
  int skipRematchProcessing = 1;

  /* Start by selecting the object which has the most nearby objects within SEARCH_NONE_MAX_ARCSEC radius */

  for (groupIndex = 0; groupIndex <  pFileCommon->blendGroupCount; groupIndex++) {
    pGroupStatsBlend = &pFileCommon->group_stats_table[groupIndex];
    pGroupStats3 = NULL;

    if (pGroupStatsBlend->blendGroupNumber < 0) {
      continue;
    }

    skipRematchProcessing = 1;

    if (skipRematchProcessing != 0) {
      bestFlag = 1;
      maxRefcnt = 0;
      maxRefcnt2 = 0;
      maxRefcntIndex = -1;
      maxRefcntIndex2 = -1;
      nextFlink = pGroupStatsBlend->blendFlink;
      groupMemberCount = pGroupStatsBlend->blendGroupMemberCount;
      if (groupMemberCount == 1) {
        pCurStarImage = &pNoneMagnitudeTable[nextFlink];
        if ((pCurStarImage->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0) { /* change of Jul 13, 2015*/
          maxRefcnt = 1;
          maxRefcntIndex = nextFlink;
          maxRefcntIndex2 = nextFlink;
        }
      } else {
        maxRefcntIndex = -1;
        while (nextFlink >= 0) {
          pNoneStats = &pFileCommon->none_stats_table[nextFlink];
          pNoneStats->blendRefcnt = 0;
          pCurStarImage = &pNoneMagnitudeTable[nextFlink];

          if ((pCurStarImage->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0) { /* change of Jul 13, 2015*/


            nextFlink2 = pGroupStatsBlend->blendFlink;

            groupMemberCount2 = pGroupStatsBlend->blendGroupMemberCount;
            while (nextFlink2 >= 0) {
              pNoneStats2 = &pFileCommon->none_stats_table[nextFlink2];
              pCurStarImage2 = &pNoneMagnitudeTable[nextFlink2];

              if ((nextFlink2 != nextFlink) &&
                  ((pCurStarImage2->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0)) { /* change of Jul 13, 2015*/
                if ((pCurStarImage->seriesId != pCurStarImage2->seriesId) ||
                    (pCurStarImage->plateNumber != pCurStarImage2->plateNumber)) {
                  factor = cos(DEGREES_TO_RAD*((pCurStarImage->dec+pCurStarImage2->dec)/2.0));
                  testra = pCurStarImage->ra;
                  testdec = pCurStarImage->dec;
                  testdra = 3600.*(factor * (pCurStarImage2->ra - testra));
                  testddec = 3600.*(pCurStarImage2->dec - testdec);
                  testdrad =  sqrt(sqr(testdra)+sqr(testddec));
                  if (testdrad < SEARCH_NONE_MAX_ARCSEC) {
                    pNoneStats->blendRefcnt++;
                    if ((pCurStarImage->AFLAGS & (~FILTER_AMASK_PLOT)) == 0) {
                      if (pNoneStats->blendRefcnt > maxRefcnt) {
                        maxRefcnt = pNoneStats->blendRefcnt;
                      }
                    }
                    if (pNoneStats->blendRefcnt > maxRefcnt2) {
                      maxRefcnt2 = pNoneStats->blendRefcnt;
                    }
                  }
                }
              }

              nextFlink2 = pNoneStats2->blendFlink;
              if (--groupMemberCount2 < 0) {
                printf("ERROR: incorrect group member count in line %d\n",__LINE__);
                exit(1);
              }
            }
            if (groupMemberCount2 != 0) {
              printf("ERROR: incorrect group member count in line %d\n",__LINE__);
              exit(1);
            }

          }
          nextFlink = pNoneStats->blendFlink;
          if (--groupMemberCount < 0) {
            printf("ERROR: incorrect group member count in line %d\n",__LINE__);
            exit(1);
          }
        }
        if (groupMemberCount != 0) {
          printf("ERROR: incorrect group member count in line %d\n",__LINE__);
          exit(1);
        }
        if (maxRefcnt <= 0) {
          maxRefcnt = maxRefcnt2;
          bestFlag = 0;
        }

        /* We have the maximum reference count, now search for the best image with the smallest area */
        nextFlink = pGroupStatsBlend->blendFlink;
        while (nextFlink >= 0) {
          pNoneStats = &pFileCommon->none_stats_table[nextFlink];
          pCurStarImage = &pNoneMagnitudeTable[nextFlink];

          if ((pNoneStats->blendRefcnt == maxRefcnt) &&
              ((pCurStarImage->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0)) { /* change of Jul 13, 2015*/
            if ((bestFlag == 0) || ((pCurStarImage->AFLAGS & (~FILTER_AMASK_PLOT)) == 0)) {
              /* Second quality points first */
              if (maxRefcntIndex < 0) {
                maxRefcntIndex = nextFlink;
              } else {
                pNoneStats2 = &pFileCommon->none_stats_table[maxRefcntIndex];
                if (pNoneStats->imageArea  < pNoneStats2->imageArea) {
                  maxRefcntIndex = nextFlink;
                }
              }
              /* Now see if we have any points without Plate Defect or SExtractor blend set */
              if ((pCurStarImage->AFLAGS & ((1 << FILTER_AFLAG_DEFECT) | (1 << FILTER_AFLAG_BLEND))) == 0) {
                if (maxRefcntIndex2 < 0) {
                  maxRefcntIndex2 = nextFlink;
                } else {
                  pNoneStats2 = &pFileCommon->none_stats_table[maxRefcntIndex2];
                  if (pNoneStats->imageArea  < pNoneStats2->imageArea) {
                    maxRefcntIndex2 = nextFlink;
                  }
                }

              }
            }
          }
          nextFlink = pNoneStats->blendFlink;
        }

      }
      /* Select the best image if available */
      if (maxRefcntIndex2 >= 0) {
        maxRefcntIndex = maxRefcntIndex2;
      }

      /* Make this image the root of a new group */

      if (maxRefcntIndex < 0) {
        /* Change of Jul 13, 2015 - reject the entire group because it contains only defects */
        nextFlink = pGroupStatsBlend->blendFlink;
        while (nextFlink >= 0) {
          pNoneStats = &pFileCommon->none_stats_table[nextFlink];
          pCurStarImage = &pNoneMagnitudeTable[nextFlink];
          nextFlink = pNoneStats->blendFlink;
          pNoneStats->selected = 2*magnitudeCount;
        }
        pFileCommon->counterBlock.rejectedTCGroup++;
        pGroupStatsBlend->blendGroupNumber = -1;
        pGroupStatsBlend->groupNumber = -1;
        pGroupStatsBlend->blendGroupMemberCount = 0;
        continue;
      }
    }

    pNoneStats = &pFileCommon->none_stats_table[maxRefcntIndex];
    pCurStarImage = &pNoneMagnitudeTable[maxRefcntIndex];
    GetAFLAGSStats(pFileCommon,pCurStarImage->AFLAGS);
    pNoneStats->selected = maxRefcntIndex;
    pNoneStats->curdrad = 0.0;
    pGroupStatsTC = &pFileCommon->group_stats_table[pFileCommon->numGroups];
    pNoneStats->groupNumber = pFileCommon->numGroups;
    pGroupStatsTC->groupNumber = pFileCommon->numGroups;
    pFileCommon->numGroups++;
    pGroupStatsTC->groupIndex = maxRefcntIndex;
    pGroupStatsTC->flink = maxRefcntIndex;
    pGroupStatsTC->groupCount++;
    pGroupStatsTC->magcal_magdep_count++;
    pGroupStatsTC->magcal_magdep += pCurStarImage->magcal_magdep;
    while (pCurStarImage->ra >= 360.0) {  /* Fix introduced Sep 22, 2018 */
      pCurStarImage->ra -= 360.0;
    }
    while (pCurStarImage->ra < 0.0) {
      pCurStarImage->ra += 360.0;
    }
    pGroupStatsTC->ra = pCurStarImage->ra;
    pGroupStatsTC->dec = pCurStarImage->dec;
    pGroupStatsTC->maxRA = pCurStarImage->ra;
    pGroupStatsTC->maxDec = pCurStarImage->dec;
    pGroupStatsTC->minRA = pCurStarImage->ra;
    pGroupStatsTC->minDec = pCurStarImage->dec;

    result = GetDASCHNumber(pCurStarImage->ra,pCurStarImage->dec,&pGroupStatsTC->REFNumber,0,REF_TYPE_DASCH);
    if (result < 0) {
      printf("ERROR: line %d photometryutils GetDASCHNumber ra %f dec %f\n",__LINE__,pCurStarImage->ra,pCurStarImage->dec);
      exit(1);
    }
    /* At this point, make sure that the REFNumber does not agree with any previous member.  If so, tweak it by 1 arcsec */

    dupFound = 1;
    while (dupFound == 1) {
      dupFound = 0;
      for (groupIndex2 = 0; groupIndex2 < pFileCommon->numGroups; groupIndex2++) {
        if (groupIndex2 == pGroupStatsTC->groupNumber) {
          continue;
        }
        pGroupStatsTC2 = &pFileCommon->group_stats_table[groupIndex2];
        if (pGroupStatsTC2->REFNumber == pGroupStatsTC->REFNumber) {
          dupFound = 1;
          pGroupStatsTC->REFNumber++;
          break;
        }
      }
    }

    pNoneStats->REFNumber = pGroupStatsTC->REFNumber;
    pGroupStatsTC->factor = cos(DEGREES_TO_RAD*pCurStarImage->dec);
    pGroupStatsTC->drad2 = pNoneStats->drad2;
    if (pCurStarImage->Date == 0) {  /* Database corruption found Sep 24, 2018 */
      pGroupStatsTC->minEpoch = PIPELINE_MIN_DATE;
    } else {
      pGroupStatsTC->minEpoch =  jd2ep(pCurStarImage->Date);
    }
    pGroupStatsTC->maxEpoch = pGroupStatsTC->minEpoch;
    pGroupStatsTC->gsc_bin_index = GetGSCBin(pGscBin,pCurStarImage->ra,pCurStarImage->dec,&decBin,&raBin,"ProcessNoneImages");
    pGroupStatsTC->GSCBlendFlag = SearchGSCBlends(pCurStarImage->ra,pCurStarImage->dec,pSky2KTable,sky2k_size,maxGSCRadius);

    /* Now add all of the other group members to the group providing that there are no plate conflicts
       and the added members are all within SEARCH_NONE_MAX_ARCSEC of the root member */
    nextFlink = pGroupStatsBlend->blendFlink;
    while (nextFlink >= 0) {
      curFlink = nextFlink;
      pNoneStats2 = &pFileCommon->none_stats_table[curFlink];

      pGroupStatsBlend->blendFlink = pNoneStats2->blendFlink;
      pGroupStatsBlend->blendGroupMemberCount--;
      nextFlink = pNoneStats2->blendFlink;
      pNoneStats2->blendFlink = -1;
      pNoneStats2->blendGroupNumber = -1;
      if (curFlink != maxRefcntIndex) {
        pCurStarImage2 = &pNoneMagnitudeTable[curFlink];
        groupMemberIndex = pGroupStatsTC->flink;
        groupFlinkPtr = &pGroupStatsTC->flink;
        for (groupVectorIndex = 0; groupVectorIndex < pGroupStatsTC->groupCount; groupVectorIndex++) {
          pGroupStarImage =  &pNoneMagnitudeTable[groupMemberIndex];
          pGroupNoneStats = &pFileCommon->none_stats_table[groupMemberIndex];
          if ((pGroupStarImage->seriesId == pCurStarImage2->seriesId) &&
              /* Remove (pGroupStarImage->solutionNumber == pCurStarImage2->solutionNumber) Jul 13, 2015 */
              (pGroupStarImage->plateNumber == pCurStarImage2->plateNumber)) {
            /* Here we have a duplicate plate in the group. Find out which of the two is further from the group root member and move it into a new group */
            distance1 = sqrt(sqr(pGroupStatsTC->factor*(pGroupStarImage->ra-pCurStarImage->ra))+sqr(pGroupStarImage->dec-pCurStarImage->dec));
            distance2 = sqrt(sqr(pGroupStatsTC->factor*(pCurStarImage2->ra) -pCurStarImage->ra)+sqr(pCurStarImage2->dec -pCurStarImage->dec));
            if ((groupMemberIndex == pGroupStatsTC->groupIndex) || /* Never toss the group root! */
                (distance2 > SEARCH_NONE_MAX_ARCSEC) ||
                (distance1 < distance2)) {
              badFlink = curFlink;
            } else {
              /* Link the new one in the group */
              badFlink = groupMemberIndex;
              *groupFlinkPtr = curFlink;  /* Link the new one to the group */
              pNoneStats2->flink = pGroupNoneStats->flink;
              pNoneStats2->selected = pGroupStatsTC->groupIndex;
              pNoneStats2->REFNumber = pGroupStatsTC->REFNumber;
              pNoneStats2->groupNumber = pGroupStatsTC->groupNumber;
              pGroupNoneStats->flink = -1;
              pGroupNoneStats->blendFlink = -1;
              pGroupNoneStats->selected = -1;
              pGroupNoneStats->REFNumber = 0;
              pGroupNoneStats->groupNumber = -1;
              pCurStarImage2 = pGroupStarImage;
            }
            pNoneStats3 = &pFileCommon->none_stats_table[badFlink];
            if (pGroupStats3 == NULL) {
              /* Need to allocate a new group */
              pGroupStats3 = &pFileCommon->group_stats_table[pFileCommon->blendGroupCount];


              pGroupStats3->blendGroupNumber = pFileCommon->blendGroupCount;
              pGroupStats3->blendFlink = badFlink;
              pGroupStats3->blendGroupMemberCount = 1;
              pNoneStats3->blendGroupNumber = pGroupStats3->blendGroupNumber;
              pNoneStats3->blendFlink = -1;
              pNoneStats3->blendRefcnt = 0;
              pFileCommon->blendGroupCount++;

            } else {
              /* Adding this to the already created new group */
              pNoneStats3->blendFlink = pGroupStats3->blendFlink;
              pNoneStats3->blendGroupNumber =  pGroupStats3->blendGroupNumber;
              pNoneStats3->blendRefcnt = 0;
              pGroupStats3->blendFlink = badFlink;
              pGroupStats3->blendGroupMemberCount++;
            }
            break;
          }
          groupMemberIndex = pGroupNoneStats->flink;
          groupFlinkPtr = &pGroupNoneStats->flink;
          if (groupVectorIndex == pGroupStatsTC->groupCount-1) {
            if (groupMemberIndex > 0) {
              printf("ERROR: Illegal end groupMemberIndex %d in line %d\n",groupMemberIndex,__LINE__);
              exit(1);
            }
          } else {
            if ((groupMemberIndex < 0) ||	(groupMemberIndex >= magnitudeCount)) {
              printf("ERROR: Illegal groupMemberIndex %d %d in line %d\n",groupMemberIndex,magnitudeCount,__LINE__);
              exit(1);
            }
          }
        }
        if (groupVectorIndex == pGroupStatsTC->groupCount) {
          /* No match, so o.k to add this one */
          pNoneStats2->selected = pGroupStatsTC->groupIndex;
          pNoneStats2->REFNumber = pGroupStatsTC->REFNumber;
          pNoneStats2->groupNumber = pGroupStatsTC->groupNumber;
          pNoneStats2->flink = pGroupStatsTC->flink;
          pGroupStatsTC->flink = curFlink;
          pGroupStatsTC->groupCount++;
        }
      }


    }

  }

  /* Now that the membership of the blend groups is stable, update the statistics */
  /* AllocFileCommonVectors(pFileCommon,magnitudeCount); */
  for (groupIndex = 0; groupIndex <  pFileCommon->numGroups; groupIndex++) {

    pGroupStatsTC = &pFileCommon->group_stats_table[groupIndex];

    groupMemberIndex = pGroupStatsTC->flink;
    positionCount = 0;
    raAverage = 0;
    decAverage = 0;
    positionWeight = 0;
    for (groupVectorIndex = 0; groupVectorIndex < pGroupStatsTC->groupCount; groupVectorIndex++) {
      pGroupStarImage =  &pNoneMagnitudeTable[groupMemberIndex];
      pGroupNoneStats = &pFileCommon->none_stats_table[groupMemberIndex];
      if (((pGroupStarImage->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0) &&
          (pGroupStarImage->dradRMS2 > 0) &&
          (pGroupStarImage->dradRMS2 < 90.0)) {
        raAverage += (pGroupStarImage->ra+360.)/sqr(pGroupStarImage->dradRMS2);
        decAverage += (pGroupStarImage->dec+180.)/sqr(pGroupStarImage->dradRMS2);
        positionWeight += 1.0/sqr(pGroupStarImage->dradRMS2);
        positionCount++;
      }
      if (groupVectorIndex != 0) {

        UpdateGroupRange(pGroupStatsTC,pGroupStarImage->ra,pGroupStarImage->dec);
        pGroupStatsTC->magcal_magdep = pGroupStarImage->magcal_magdep + (pGroupStatsTC->magcal_magdep*pGroupStatsTC->magcal_magdep_count);
        pGroupStatsTC->magcal_magdep_count++;
        pGroupStatsTC->magcal_magdep = pGroupStatsTC->magcal_magdep/pGroupStatsTC->magcal_magdep_count;
        if (pGroupStarImage->Date != 0) {  /* Database corruption found Sep 24, 2018 */

          tempEpoch = jd2ep(pGroupStarImage->Date);
          if (pGroupStatsTC->minEpoch == PIPELINE_MIN_DATE) {
            pGroupStatsTC->minEpoch = tempEpoch;
            pGroupStatsTC->maxEpoch = tempEpoch;
          } else {
            if (tempEpoch < pGroupStatsTC->minEpoch) {
              pGroupStatsTC->minEpoch = tempEpoch;
            }
            if (tempEpoch > pGroupStatsTC->maxEpoch) {
              pGroupStatsTC->maxEpoch = tempEpoch;
            }
          }
        }
      }
      groupMemberIndex = pGroupNoneStats->flink;
      if (groupVectorIndex == pGroupStatsTC->groupCount-1) {
        if (groupMemberIndex > 0) {
          printf("ERROR: Illegal end groupMemberIndex %d in line %d\n",groupMemberIndex,__LINE__);
          exit(1);
        }
      } else {
        if ((groupMemberIndex < 0) ||	(groupMemberIndex >= magnitudeCount)) {
          printf("ERROR: Illegal groupMemberIndex %d %d in line %d\n",groupMemberIndex,magnitudeCount,__LINE__);
          exit(1);
        }
      }
    }
    if (positionCount > 1) {
      raAverage = raAverage/positionWeight;
      while (raAverage >= 360.0) {
        raAverage -= 360.0;
      }
      while (raAverage < 0.0) {
        raAverage += 360.0;
      }
      decAverage = decAverage/positionWeight;
      while (decAverage > 90.0) {
        decAverage -= 180.0;
      }
      result = GetDASCHNumber(raAverage,decAverage,&newREFNumber,0,REF_TYPE_DASCH);
      if (result < 0) {
        printf("ERROR: line %d photometryutils GetDASCHNumber ra %f dec %f\n",__LINE__,raAverage,decAverage);
        exit(1);
      }
      if (newREFNumber != pGroupStatsTC->REFNumber) {
        /* Again, check to see if the new REFNumber matches a previous REFNumber */
        dupFound = 1;
        while (dupFound == 1) {
          dupFound = 0;
          for (groupIndex2 = 0; groupIndex2 < pFileCommon->numGroups; groupIndex2++) {
            if (groupIndex2 == pGroupStatsTC->groupNumber) {
              continue;
            }
            pGroupStatsTC2 = &pFileCommon->group_stats_table[groupIndex2];
            if (pGroupStatsTC2->REFNumber == newREFNumber) {
              dupFound = 1;
              newREFNumber++;
              decAverage += 1.0/3600.;
              break;
            }
          }
        }

      }
      if (newREFNumber != pGroupStatsTC->REFNumber) {
        /* Now update all entries with the new REFNumber */
        pGroupStatsTC->REFNumber = newREFNumber;
        pGroupStatsTC->ra = raAverage;
        pGroupStatsTC->dec = decAverage;
        groupMemberIndex = pGroupStatsTC->flink;
        for (groupVectorIndex = 0; groupVectorIndex < pGroupStatsTC->groupCount; groupVectorIndex++) {
          pGroupStarImage =  &pNoneMagnitudeTable[groupMemberIndex];
          pGroupNoneStats = &pFileCommon->none_stats_table[groupMemberIndex];
          pGroupNoneStats->REFNumber = newREFNumber;
          groupMemberIndex = pGroupNoneStats->flink;
          if (groupVectorIndex == pGroupStatsTC->groupCount-1) {
            if (groupMemberIndex > 0) {
              printf("ERROR: Illegal end groupMemberIndex %d in line %d\n",groupMemberIndex,__LINE__);
              exit(1);
            }
          } else {
            if ((groupMemberIndex < 0) ||	(groupMemberIndex >= magnitudeCount)) {
              printf("ERROR: Illegal groupMemberIndex %d %d in line %d\n",groupMemberIndex,magnitudeCount,__LINE__);
              exit(1);
            }
          }
        }

      }

    }
  }
}


// Given a bunch of magnitudes loaded via LocateNoneImages ... do stuff?
//
// "This routine will modify the REFNumber and gsc_bin_index of the input array.
// Therefore, the input array should be a copy of the files read from the
// database."
void
ProcessNoneImagesX(
  PGSCBIN pGscBin,
  PFILECOMMON pFileCommon,
  PPHOTTARGET *pTarget_table,
  size_t *pTarget_nrecs,
  size_t *pTarget_alloc,
  PFILESTARIMAGE pNoneMagnitudeTable,
  int magnitudeCount,
  int verbose,
  int *histogramTable,
  int *pGroupNumber,
  int curGscBinIndex
) {
  char REF[MAX_REF];
  int initialGroupNumber = *pGroupNumber;
  PFILESTARIMAGE pCurStarImage;
  PFILESTARIMAGE pCurStarImage2;
  PPMSTATS pPmStats; /* Proper motion stats */
  PPMSTATS pPmStats2; /* Proper motion stats */
  PAREASORT pAreaSort; /* Object area sorting */
  PNONESTATS pNoneStats;
  PNONESTATS pNoneStats2;
  PGROUPSTATS pGroupStats;
  PNEARESTCATALOGSTAR pNearestCatalogStar;
  double testdrad;
  PPHOTTARGET target_table = NULL;
  PPHOTTARGET tmp_target_table = NULL;
  PPHOTTARGET pTarget = NULL;
  size_t target_nrecs = 0;
  size_t target_alloc = 0;
  int decCBin;
  int raCBin;
  PGSCBIN pGscBin01 = &gscBin01;
  char series[MAX_SERIES_STRING];
  int nearbyCatalogIndex;

  PNEARBYCATALOG pNearbyCatalog = NULL;
  double minNearestDist = 0;
  double minNearestDistGoodMag = 0;
  int minNearestIndex;
  int minNearestIndexGoodMag;
  double tmpdist;
  double tmpNearestDist;
  double tmpNearestDistGoodMag;
  int nearestIndex;
  int quality;
  int duplicatePlate;
  PGALAXYCOMMON pGalaxyCommon = &pFileCommon->galaxyCommon;
  int properMotionIndex = 0; /* Next free index in the proper motion table */
  int magnitudeIndex;
  int magnitudeIndex2;
  int groupIndex;
  double testra;
  double testdec;
  double testdra;
  double testddec;
  double factor;
  double testaLength;
  double testbLength;
  int blendCandidateCount = 0;
  double *pColorterm = NULL;
  int *pColorflag = NULL;

  if (magnitudeCount <= 0) {
    return;
  }

  if (pFileCommon->startTime != 0) {
    time_t newTime;
    time(&newTime);
    newTime -= pFileCommon->startTime;

    if (newTime > WEB_TIMEOUT) {
      printf("ERROR: Request processing timed out in line %d after %ld seconds\n", __LINE__, newTime);
      exit(1);
    }
  }

  pFileCommon->numGroups = 0;

  if (pTarget_table != NULL && pTarget_nrecs != NULL && pTarget_alloc != NULL) {
    target_table = *pTarget_table;
    target_nrecs = *pTarget_nrecs;
    target_alloc = *pTarget_alloc;
  }

  if (magnitudeCount > pFileCommon->statsAllocCount) {
    pFileCommon->statsAllocCount = magnitudeCount + STATS_REALLOC_INCREMENT;

    if (pFileCommon->none_stats_table != NULL) {
      free(pFileCommon->none_stats_table);
    }

    pFileCommon->none_stats_table = (PNONESTATS)calloc(pFileCommon->statsAllocCount+3,sizeof(NONESTATS));
    if (pFileCommon->none_stats_table == NULL) {
      printf("ERROR: failed to allocate %d magnitudes in none_stats_table\n", pFileCommon->statsAllocCount + 3);
      exit(1);
    }

    if (pFileCommon->pm_stats_table != NULL) {
      free(pFileCommon->pm_stats_table);
    }

    pFileCommon->pm_stats_table = (PPMSTATS)calloc(pFileCommon->statsAllocCount+3,sizeof(PMSTATS));
    if (pFileCommon->pm_stats_table == NULL) {
      printf("ERROR: failed to allocate %d magnitudes in pm_stats_table\n", pFileCommon->statsAllocCount + 3);
      exit(1);
    }

    if (pFileCommon->area_sort_table != NULL) {
      free(pFileCommon->area_sort_table);
    }

    pFileCommon->area_sort_table = (PAREASORT)calloc(pFileCommon->statsAllocCount+3,sizeof(AREASORT));
    if (pFileCommon->area_sort_table == NULL) {
      printf("ERROR: failed to allocate %d magnitudes in area_sort_table\n", pFileCommon->statsAllocCount + 3);
      exit(1);
    }

    if (pFileCommon->group_stats_table != NULL) {
      free(pFileCommon->group_stats_table);
    }

    pFileCommon->group_stats_table = (PGROUPSTATS)calloc(pFileCommon->statsAllocCount+3,sizeof(GROUPSTATS));
    if (pFileCommon->group_stats_table == NULL) {
      printf("ERROR: failed to allocate %d magnitudes in group_stats_table\n", pFileCommon->statsAllocCount + 3);
      exit(1);
    }

    if (pFileCommon->nearby_catalog_table != NULL) {
      free(pFileCommon->nearby_catalog_table);
    }

    pFileCommon->nearby_catalog_table = (PNEARBYCATALOG)calloc(pFileCommon->statsAllocCount+3,sizeof(NEARBYCATALOG));
    if (pFileCommon->nearby_catalog_table == NULL) {
      printf("ERROR: failed to allocate %d magnitudes in nearby_catalog_table\n", pFileCommon->statsAllocCount + 3);
      exit(1);
    }
  } else {
    memset(pFileCommon->none_stats_table, 0, pFileCommon->statsAllocCount * sizeof(NONESTATS));
    memset(pFileCommon->pm_stats_table, 0, pFileCommon->statsAllocCount * sizeof(PMSTATS));
    memset(pFileCommon->area_sort_table, 0, pFileCommon->statsAllocCount * sizeof(AREASORT));
    memset(pFileCommon->group_stats_table, 0, pFileCommon->statsAllocCount * sizeof(GROUPSTATS));
    memset(pFileCommon->nearby_catalog_table, 0, pFileCommon->statsAllocCount * sizeof(NEARBYCATALOG));
  }

  pFileCommon->nearbyCatalogCount = 0;
  pFileCommon->blendGroupCount = 0;

  // Pass 1: initialization

  for (magnitudeIndex = 0; magnitudeIndex < magnitudeCount; magnitudeIndex++) {
    pPmStats = &pFileCommon->pm_stats_table[magnitudeIndex];
    pNoneStats = &pFileCommon->none_stats_table[magnitudeIndex];
    pGroupStats = &pFileCommon->group_stats_table[magnitudeIndex];
    pAreaSort = &pFileCommon->area_sort_table[magnitudeIndex];
    pCurStarImage = &pNoneMagnitudeTable[magnitudeIndex];

    pNoneStats->selected = -1;
    pNoneStats->magnitudeIndex = magnitudeIndex;
    pNoneStats->groupNumber = -1;
    pNoneStats->flink = -1;
    pNoneStats->REFNumber = 0;

    pGroupStats->groupIndex = -1;
    pGroupStats->groupCount = 0;
    pGroupStats->GSCBlendFlag = 0;
    pGroupStats->flink = -1;
    pGroupStats->magcal_magdep_count = 0;
    pGroupStats->magcal_magdep = 0;
    memset(pPmStats, 0, sizeof(PMSTATS));
    pPmStats->link = -1;
    InitNearestCatalogStar(&pGroupStats->nearestCatalogStar);
    pAreaSort->magnitudeIndex = magnitudeIndex;
    pAreaSort->pNoneStats = pNoneStats;
    pNoneStats->blendGroupNumber = -1;
    pNoneStats->blendFlink = -1;
    pGroupStats->blendGroupNumber = -1;
    pGroupStats->blendFlink = -1;
  }

  // Pass 2

  for (magnitudeIndex = 0; magnitudeIndex < magnitudeCount; magnitudeIndex++ ) {
    pNoneStats = &pFileCommon->none_stats_table[magnitudeIndex];
    pCurStarImage = &pNoneMagnitudeTable[magnitudeIndex];
    pPmStats = &pFileCommon->pm_stats_table[properMotionIndex];

    if (pCurStarImage->REFNumber != 0) {
      // "This is a catalog object. All we need to save are the position and
      // reference number. The magnitudeBuffer array should be sorted by
      // REFNumber."

      // Search for a match in the "nearby catalog".
      for (nearbyCatalogIndex = 0; nearbyCatalogIndex < pFileCommon->nearbyCatalogCount; nearbyCatalogIndex++) {
        pNearbyCatalog = &pFileCommon->nearby_catalog_table[nearbyCatalogIndex];
        if (pNearbyCatalog->REFNumber == pCurStarImage->REFNumber) {
          break;
        }
      }

      if (nearbyCatalogIndex == pFileCommon->nearbyCatalogCount) {
        // No exact match in the nearby objects catalog. Allocate a new entry
        // for just this source.

        pNearbyCatalog = &pFileCommon->nearby_catalog_table[pFileCommon->nearbyCatalogCount];
        memset(pNearbyCatalog, 0, sizeof(NEARBYCATALOG));
        pNearbyCatalog->firstProperMotionIndex = -1;
        pNearbyCatalog->lastProperMotionIndex = -1;
        pNearbyCatalog->REFNumber = pCurStarImage->REFNumber;
        pNearbyCatalog->gsc_bin_index = pCurStarImage->gsc_bin_index;

        // Stats for all points:

        pNearbyCatalog->ra = 999.0;
        pNearbyCatalog->dec = 99.0;
        pNearbyCatalog->dradRMS2 = pCurStarImage->dradRMS2;

        if (pCurStarImage->magcal_magdep < 90.0) {
          pNearbyCatalog->magcal_magdep = pCurStarImage->magcal_magdep;
          pNearbyCatalog->magcal_magdep_count++;
        }

        if ((pCurStarImage->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) == 0 && pCurStarImage->magcal_magdep < 90.0) {
          // Stats for "good magnitudes":

          pNearbyCatalog->raGoodMag = 999.0;
          pNearbyCatalog->decGoodMag = 99.0;
          pNearbyCatalog->dradRMS2GoodMag = pCurStarImage->dradRMS2;

          if (pCurStarImage->magcal_magdep < 90.0) {
            pNearbyCatalog->magcal_magdepGoodMag = pCurStarImage->magcal_magdep;
            pNearbyCatalog->magcal_magdepGoodMag_count++;
          }
        } else {
          pNearbyCatalog->dradRMS2GoodMag = -1.0;
        }

        pFileCommon->nearbyCatalogCount++;
      } else {
        // We found an exact match in the nearby objects catalog

        if (pCurStarImage->magcal_magdep < 90.0) {
          pNearbyCatalog->magcal_magdep += pCurStarImage->magcal_magdep;
          pNearbyCatalog->magcal_magdep_count++;
        }

        if (pCurStarImage->dradRMS2 < pNearbyCatalog->dradRMS2 || pNearbyCatalog->dradRMS2 < 0.0) {
          pNearbyCatalog->ra = 999.0;
          pNearbyCatalog->dec = 99.0;
          pNearbyCatalog->dradRMS2 = pCurStarImage->dradRMS2;
        }

        if ((pCurStarImage->AFLAGS & (1<<FILTER_AFLAG_LIMITING_MAG)) == 0 && pCurStarImage->magcal_magdep < 90.0) {
          // Augment stats for "good magnitudes":

          pNearbyCatalog->magcal_magdepGoodMag += pCurStarImage->magcal_magdep;
          pNearbyCatalog->magcal_magdepGoodMag_count++;

          if ((pCurStarImage->dradRMS2 < pNearbyCatalog->dradRMS2GoodMag) || (pNearbyCatalog->dradRMS2GoodMag < 0.0)) {
            pNearbyCatalog->raGoodMag = 999.0;
            pNearbyCatalog->decGoodMag = 99.0;
            pNearbyCatalog->dradRMS2GoodMag = pCurStarImage->dradRMS2;
          }
        }
      }

      if (
        (pCurStarImage->RaPM < MAX_PROPER_MOTION) &&
        (pCurStarImage->DecPM < MAX_PROPER_MOTION) &&
        (GetFittedPlateScale(pCurStarImage->seriesId,pCurStarImage->plateNumber) < PATROL_PLATE_SCALE)
      ) {
        // "This star has a valid proper motion. Calculate its catalog
        // position." Accumulate proper motion stats.
        pPmStats->REFNumber = pCurStarImage->REFNumber;
        pPmStats->ra = pCurStarImage->ra;
        pPmStats->dec = pCurStarImage->dec;
        pPmStats->ra_2 = pCurStarImage->ra_2;
        pPmStats->dec_2 = pCurStarImage->dec_2;
        pPmStats->RaPM = pCurStarImage->RaPM;
        pPmStats->DecPM = pCurStarImage->DecPM;
        pPmStats->AFLAGS = pCurStarImage->AFLAGS;
        pPmStats->BFLAGS = pCurStarImage->BFLAGS;
        pPmStats->ELLIPTICITY = pCurStarImage->ELLIPTICITY;
        pPmStats->epoch = jd2ep(pCurStarImage->Date) - GSC_EQUINOX;
        pPmStats->catdec = pCurStarImage->dec_2 - (pCurStarImage->DecPM * pPmStats->epoch) / (3600.0 * 1000.0);

        if (pPmStats->catdec < 90.0 && pPmStats->catdec > -90.0) {
          factor = cos(pPmStats->catdec * DEGREES_TO_RAD);
          pPmStats->catra = pCurStarImage->ra_2 - ((pCurStarImage->RaPM * pPmStats->epoch) / (3600.0 * 1000.0)) / factor;
        } else {
          pPmStats->catra = pCurStarImage->ra_2;
          factor = 1.0;
        }

        if (pNearbyCatalog->properMotionCount == 0) {
          pNearbyCatalog->firstProperMotionIndex = properMotionIndex;
          pNearbyCatalog->lastProperMotionIndex = properMotionIndex;
        } else {
          pPmStats2 = &pFileCommon->pm_stats_table[pNearbyCatalog->lastProperMotionIndex];
          pPmStats2->link = properMotionIndex;
          pNearbyCatalog->lastProperMotionIndex = properMotionIndex;
        }

        pNearbyCatalog->properMotionCount++;
        pPmStats->count = pNearbyCatalog->properMotionCount;
        pPmStats->link = -1;
        properMotionIndex++;
      }
    }

    // That was all processing for records with REFnumber != 0. For *all* records ....
    // "Select the better unmatched objects"

    GetPlateQualityMask(
      pFileCommon,
      pCurStarImage->seriesId,
      pCurStarImage->plateNumber,
      &quality,
      &pColorterm,
      &pColorflag,
      pCurStarImage->REFNumber
    );

    pNoneStats->colorterm = pColorterm[pCurStarImage->spatial_bin - 1];
    pNoneStats->colorflag = pColorflag[pCurStarImage->spatial_bin - 1];

    if (pCurStarImage->REFNumber == 0 || pFileCommon->enableRematch != 0) {
      // For unmatched objects with REFnumber = 0, accumulate "none stats".
      pNoneStats->scale = GetFittedPlateScale(pCurStarImage->seriesId, pCurStarImage->plateNumber);
      pNoneStats->MAG_ISO = pCurStarImage->MAG_ISO;
      pNoneStats->drad = pNoneStats->scale * MaxDradPixels("None",1);

      if (pCurStarImage->dradRMS2 > 0 && pCurStarImage->dradRMS2 / 3600.0 < pNoneStats->drad) {
        pNoneStats->drad = pCurStarImage->dradRMS2 / 3600.0;
      }

      pNoneStats->drad2 = pNoneStats->drad;
      if (pNoneStats->drad2 < (SEARCH_NONE_MIN_DRAD * pNoneStats->scale)) {
        pNoneStats->drad2 = (SEARCH_NONE_MIN_DRAD * pNoneStats->scale);
      }

      if (pCurStarImage->ELLIPTICITY < 1.0) {
        pNoneStats->aLength = pNoneStats->scale * sqrt((0.5 * (pCurStarImage->ISO3 + pCurStarImage->ISO4)) / (PI_VALUE * (1.0 - pCurStarImage->ELLIPTICITY)));
      } else {
        pNoneStats->aLength = pNoneStats->scale * sqrt((0.5 * (pCurStarImage->ISO3 + pCurStarImage->ISO4)) / (PI_VALUE));
      }

      pNoneStats->bLength = pNoneStats->aLength * (1.0 - pCurStarImage->ELLIPTICITY);

      // "We will be sorting on this area to select our candidates."
      pNoneStats->imageArea = 4.0 * (pNoneStats->aLength + pNoneStats->drad2) * (pNoneStats->bLength + pNoneStats->drad2);
      pNoneStats->epoch = jd2ep(pCurStarImage->Date) - GSC_EQUINOX;
      pNoneStats->blendAngle =  pCurStarImage->THETA_J2000 + 90.;

      if (pNoneStats->blendAngle > 90) {
        pNoneStats->blendAngle -= 180;
      }
    }
  }

  // "Find the average magnitude of all nearby catalog objects."

  for (nearbyCatalogIndex = 0; nearbyCatalogIndex < pFileCommon->nearbyCatalogCount; nearbyCatalogIndex++) {
    pNearbyCatalog = &pFileCommon->nearby_catalog_table[nearbyCatalogIndex];

    if (pNearbyCatalog->magcal_magdep_count > 0) {
      pNearbyCatalog->magcal_magdep = pNearbyCatalog->magcal_magdep / (1.0 * pNearbyCatalog->magcal_magdep_count);
    } else {
      pNearbyCatalog->magcal_magdep = 99.0;
    }

    if (pNearbyCatalog->magcal_magdepGoodMag_count > 0) {
      pNearbyCatalog->magcal_magdepGoodMag = pNearbyCatalog->magcal_magdepGoodMag / (1.0 * pNearbyCatalog->magcal_magdepGoodMag_count);
    } else {
      pNearbyCatalog->magcal_magdepGoodMag = 99.0;
    }

    // ... and compute our own proper motions, it looks like?
    ProcessProperMotions(pFileCommon, pNearbyCatalog, curGscBinIndex);
  }

  // "Workaround to handle bug of Jul 3, 2013 without having to reload the
  // entire database. For each transient candidate, check to be sure that it is
  // not a blend with a matched catalog object. If it is, assign the REFNumber
  // to the matched object"

  {
    PFILESTARIMAGE pFileStarImage;
    double testfactor;
    double tempmagcal_magdep;
    double tempmagcal_magdep2;
    double tmpFlux;
    double blendedMag;

    // Pass 3

    for (magnitudeIndex = 0; magnitudeIndex < magnitudeCount; magnitudeIndex++ ) {
      pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
      pNoneStats = &pFileCommon->none_stats_table[magnitudeIndex];

      // Only work on "NONE" objects
      if (pFileStarImage->REFNumber != 0) {
        continue;
      }

      // For each catalog object ...
      for (nearbyCatalogIndex = 0; nearbyCatalogIndex < pFileCommon->nearbyCatalogCount; nearbyCatalogIndex++) {
        // Check for overlap, accounting for proper motion.
        pNearbyCatalog = &pFileCommon->nearby_catalog_table[nearbyCatalogIndex];
        testra = pNearbyCatalog->ra;
        testdec = pNearbyCatalog->dec;
        testfactor = cos(pNearbyCatalog->dec * DEGREES_TO_RAD);
        testdec += ((pNearbyCatalog->DecPM * pNoneStats->epoch) / (3600.0 * 1000.0));
        testra += (((pNearbyCatalog->RaPM *pNoneStats->epoch) / (3600.0 * 1000.0)) / testfactor);
        testdra = 3600. * testfactor * (pFileStarImage->ra - testra);
        testddec = 3600. * (pFileStarImage->dec - testdec);

        if (pNoneStats->epoch == 0) {
          printf("ERROR: pNoneStats epoch not initialized in line %d\n", __LINE__);
        }

        if (
          CheckBlend(
            -testdra, /* x axis is positive */
            testddec,
            pNoneStats->blendAngle,
            3600.0 * pNoneStats->aLength,
            3600.0 * pNoneStats->bLength,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          )
        ) {
          // This NONE appears to overlap the catalog object. Check whether we
          // seem to have a duplicate record for the same source from the same
          // mosaic+solution.
          duplicatePlate = 0;

          if (pNearbyCatalog->gsc_bin_index == curGscBinIndex || pFileStarImage->gsc_bin_index == curGscBinIndex) {
            duplicatePlate = CheckDuplicatePlate(
              pNoneMagnitudeTable,
              magnitudeCount,
              pNearbyCatalog->REFNumber,
              pFileStarImage->seriesId,
              pFileStarImage->plateNumber,
              pFileStarImage->solutionNumber
            );
          }

          if (duplicatePlate == 0) {
            // No, it's not a duplicate.

            if (pFileStarImage->REFNumber == 0) {
              // This is always true, from the precondition above. Copy over
              // data with the "late match" flag.
              pFileStarImage->REFNumber = pNearbyCatalog->REFNumber;
              pFileStarImage->ra_2 = testra;
              pFileStarImage->dec_2 = testdec;
              pFileStarImage->RaPM = pNearbyCatalog->RaPM;
              pFileStarImage->DecPM = pNearbyCatalog->DecPM;
              pFileStarImage->BFLAGS |= (1 << FILTER_BFLAG_LATEMATCH);

              pNoneStats->selected = 2 * magnitudeCount;

              if (pNearbyCatalog->magcal_magdepGoodMag < 90.0) {
                pNoneStats->magcal_magdep = pNearbyCatalog->magcal_magdepGoodMag;
              } else {
                pNoneStats->magcal_magdep = pNearbyCatalog->magcal_magdep;
              }
            } else {
              // This stanza should be unreachable based on the earlier precondition.

              if (pNearbyCatalog->magcal_magdepGoodMag < 90.0) {
                tempmagcal_magdep = pNearbyCatalog->magcal_magdepGoodMag;
              } else {
                tempmagcal_magdep = pNearbyCatalog->magcal_magdep;
              }

              tempmagcal_magdep2 = pNoneStats->magcal_magdep;

              if (tempmagcal_magdep < pNoneStats->magcal_magdep) {
                pFileStarImage->REFNumber = pNearbyCatalog->REFNumber;
                pFileStarImage->ra_2 = testra;
                pFileStarImage->dec_2 = testdec;
                pFileStarImage->RaPM = pNearbyCatalog->RaPM;
                pFileStarImage->DecPM = pNearbyCatalog->DecPM;
                pFileStarImage->BFLAGS |= (1 << FILTER_BFLAG_LATEMATCH);
                pNoneStats->magcal_magdep = tempmagcal_magdep;
              }

              pNoneStats->selected = 2*magnitudeCount;
              pFileStarImage->AFLAGS |= (1<<FILTER_AFLAG_BLEND);

              /* See if we meet the blend condition */
              tmpFlux = exp10(-tempmagcal_magdep/2.5) + exp10(-tempmagcal_magdep2/2.5);
              blendedMag =  -2.5*log10(tmpFlux);

              if (
                (fabs(blendedMag-tempmagcal_magdep) > MIN_BLEND_MAGNITUDE) &&
                (fabs(blendedMag-tempmagcal_magdep2) > MIN_BLEND_MAGNITUDE)
              ) {
                pFileStarImage->AFLAGS |= (1 << FILTER_AFLAG_CASEB);
              }
            }
          }
        }
      }
    }
  }

  if (pFileCommon->enableRematch == 0) {
    pFileCommon->blendGroupCount = 0;
  }

  // Pass 4

  for (magnitudeIndex = 0; magnitudeIndex < magnitudeCount; magnitudeIndex++ ) {
    pNoneStats = &pFileCommon->none_stats_table[magnitudeIndex];
    pCurStarImage = &pNoneMagnitudeTable[magnitudeIndex];
    pPmStats = &pFileCommon->pm_stats_table[properMotionIndex];

    // Only work on (remaining) "NONE" objects
    if (pCurStarImage->REFNumber != 0) {
      pNoneStats->selected = 2 * magnitudeCount;
      continue;
    }

    // Skip if: MULTIPLE_BLEND, DRADBIN, CASEB, CASEC, CASED, WEDGE, DRAD,
    // LIMITING_MAG. "High drad object, reject it."
    if ((pCurStarImage->AFLAGS & FILTER_AMASK_NONE_CANDIDATE2) != 0) {
      pNoneStats->selected = 2 * magnitudeCount;
      continue;
    }

    // "Mark this one. If we don't match it, then it will become a group of 1"

    if (pNoneStats->blendCandidate == 0) {
      blendCandidateCount++;
    }

    pNoneStats->blendCandidate = 1;

    if (pFileCommon->enableRematch == 0) {
      // N^2 inner pass if *not* in "rematch" mode

      for (magnitudeIndex2 = 0; magnitudeIndex2 < magnitudeCount; magnitudeIndex2++) {
        if (magnitudeIndex2 == magnitudeIndex) {
          continue;
        }

        pNoneStats2 = &pFileCommon->none_stats_table[magnitudeIndex2];
        pCurStarImage2 = &pNoneMagnitudeTable[magnitudeIndex2];
        pPmStats2 = &pFileCommon->pm_stats_table[magnitudeIndex2];

        // Only compare against other, remaining "NONE" objects
        if (pCurStarImage2->REFNumber != 0) {
          pNoneStats->selected = 2 * magnitudeCount;
          continue;
        }

        // Same skip filter as in the outer loop.
        if ((pCurStarImage2->AFLAGS & FILTER_AMASK_NONE_CANDIDATE2) != 0) {
          pNoneStats2->selected = 2 * magnitudeCount;
          continue;
        }

        // "Mark this one. If we don't match it, then it will become a group of 1.
        if (pNoneStats2->blendCandidate == 0) {
          blendCandidateCount++;
        }

        pNoneStats2->blendCandidate = 1;

        factor = cos(DEGREES_TO_RAD * ((pCurStarImage->dec + pCurStarImage2->dec) / 2.0));
        testra = pCurStarImage->ra;
        testdec = pCurStarImage->dec;
        testdra = 3600. * (factor * (pCurStarImage2->ra - testra));
        testddec = 3600. * (pCurStarImage2->dec - testdec);
        testdrad =  SEARCH_NONE_DRAD_FACTOR * sqrt(sqr(pNoneStats->drad2) + sqr(pNoneStats2->drad2));
        testaLength = 3600.0 * pNoneStats->aLength;

        if (testaLength < testdrad) {
          testaLength = testdrad;
        }

        if (testaLength < SEARCH_NONE_MIN_ARCSEC) {
          testaLength = SEARCH_NONE_MIN_ARCSEC;
        }

        if (testaLength > SEARCH_NONE_MAX_ARCSEC) {
          testaLength = SEARCH_NONE_MAX_ARCSEC;
        }

        testbLength = 3600.0 * pNoneStats->bLength;

        if (testbLength < testdrad) {
          testbLength = testdrad;
        }

        if (testbLength < SEARCH_NONE_MIN_ARCSEC) {
          testbLength = SEARCH_NONE_MIN_ARCSEC;
        }

        if (testbLength > SEARCH_NONE_MAX_ARCSEC) {
          testbLength = SEARCH_NONE_MAX_ARCSEC;
        }

        // "This test is less inclusive because it does not look for an overlap
        // of the images, but whether the center of each image is in the area of
        // the other image."

        if (
          CheckBlend(
            -testdra, /* x axis is positive */
            testddec,
            pNoneStats->blendAngle,
            testaLength,
            testbLength,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          )
        ) {
          // These two NONE objects seem to overlap. Mark them as part of a
          // blend group.
          InsertBlendChain(pFileCommon, magnitudeIndex, magnitudeIndex2);
        }
      }
    }
  }

  // Pass 5

  for (magnitudeIndex = 0; magnitudeIndex < magnitudeCount; magnitudeIndex++ ) {
    pNoneStats = &pFileCommon->none_stats_table[magnitudeIndex];
    pCurStarImage = &pNoneMagnitudeTable[magnitudeIndex];

    // If the object was considered in the above processing, but never blended with
    // anything else, give it its own blend "group" of size one.

    if (pNoneStats->blendCandidate != 0 && pNoneStats->blendGroupNumber < 0) {
      for (groupIndex = 0; groupIndex < pFileCommon->blendGroupCount; groupIndex++) {
        pGroupStats = &pFileCommon->group_stats_table[groupIndex];
        if (pGroupStats->blendGroupNumber < 0) {
          break;
        }
      }

      if (groupIndex >= pFileCommon->blendGroupCount) {
        groupIndex = pFileCommon->blendGroupCount;
        pGroupStats = &pFileCommon->group_stats_table[groupIndex];
        pFileCommon->blendGroupCount++;
      }

      pGroupStats = &pFileCommon->group_stats_table[groupIndex];
      pGroupStats->blendGroupNumber = groupIndex;
      pGroupStats->blendFlink = magnitudeIndex;
      pGroupStats->blendGroupMemberCount = 1;
      pNoneStats->blendGroupNumber = groupIndex;
      pNoneStats->blendFlink = -1;
      pNoneStats->selected = groupIndex;
    }
  }

  // For each blend group, calculate sensible merged properties.

  ProcessBlendGroups(
    pGscBin,
    pFileCommon,
    pNoneMagnitudeTable,
    magnitudeCount,
    NULL,
    0, // sky2k_size
    0.0, // maxGSCRadius
    curGscBinIndex,
    verbose
  );

  // "Now that we know the magnitude we can find the nearest object with a mean
  // magnitude closest to this magnitude."

  for (groupIndex = 0; groupIndex < pFileCommon->numGroups; groupIndex++) {
    minNearestIndex = -1;
    minNearestDist = 0.0;
    minNearestIndexGoodMag = -1;
    minNearestDistGoodMag = 0.0;

    pGroupStats = &pFileCommon->group_stats_table[groupIndex];
    if (pGroupStats->groupNumber < 0) {
      continue;
    }

    pNearestCatalogStar = &pGroupStats->nearestCatalogStar;
    InitNearestCatalogStar(pNearestCatalogStar);

    // Find the nearest catalog star(s)

    for (nearestIndex = 0; nearestIndex < pFileCommon->nearbyCatalogCount; nearestIndex++) {
      pNearbyCatalog = &pFileCommon->nearby_catalog_table[nearestIndex];

      // Do one search over everything, without requiring "good mag" data.
      if (pNearbyCatalog->dradRMS2 >= 0.0) {
        double nearby_magcal_magdep = 99.;

        factor = cos(DEGREES_TO_RAD * ((pGroupStats->dec + pNearbyCatalog->dec) / 2.0));
        nearby_magcal_magdep = pNearbyCatalog->magcal_magdepGoodMag;
        if (nearby_magcal_magdep > 90) {
          nearby_magcal_magdep = pNearbyCatalog->magcal_magdep;
        }

        {
          // Move the catalog star to the middle of the epoch range of the group
          double testepoch = ((pGroupStats->maxEpoch + pGroupStats->minEpoch) / 2.0) - GSC_EQUINOX;
          double testFactor;
          double testNearbydec = pNearbyCatalog->dec + ((pNearbyCatalog->DecPM * testepoch) / (3600.0 * 1000.0));
          double testNearbyra;

          if (pNearbyCatalog->dec < 90.0 && pNearbyCatalog->dec > -90.0) {
            testFactor = cos(pNearbyCatalog->dec * DEGREES_TO_RAD);
            testNearbyra = pNearbyCatalog->ra + (((pNearbyCatalog->RaPM * testepoch)/ (3600.0 * 1000.0)) / testFactor);
          } else {
            testNearbyra = pNearbyCatalog->ra;
            testFactor = 1.0;
          }

          tmpNearestDist = sqr(pGroupStats->dec - testNearbydec) + sqr(testFactor * (pGroupStats->ra - testNearbyra));
        }

        if (
          pFileCommon->enableTransientSearch > 0 &&
          nearby_magcal_magdep < 90.0 &&
          pGroupStats->gsc_bin_index == curGscBinIndex
        ) {
          int tempdrad = 3600. * sqrt(tmpNearestDist);

          if (tempdrad >= 0 && tempdrad < MAX_TC_DIST) {
            pFileCommon->catalogDistTotal[tempdrad]++;
            if (fabs(nearby_magcal_magdep - pGroupStats->magcal_magdep) < MAX_TC_NEARBY_MAG) {
              pFileCommon->catalogDistSimilarMag[tempdrad]++;
            } else {
              pFileCommon->catalogDistDifferentMag[tempdrad]++;
            }
          }
        }

        if (minNearestIndex < 0) {
          minNearestIndex = nearestIndex;
          minNearestDist = tmpNearestDist;
        } else if (tmpNearestDist < minNearestDist) {
          minNearestIndex = nearestIndex;
          minNearestDist = tmpNearestDist;
        }
      }

      // Do another search over only "good mag" stars

      if (
        pNearbyCatalog->dradRMS2GoodMag >= 0.0 &&
        fabs(pGroupStats->magcal_magdep - pNearbyCatalog->magcal_magdepGoodMag) < NEAREST_MAG_DIFFERENCE
      ) {
        factor = cos(DEGREES_TO_RAD * ((pGroupStats->dec + pNearbyCatalog->decGoodMag) / 2.0));

        {
          // Move the catalog star to the middle of the epoch range of the group
          double testepoch = ((pGroupStats->maxEpoch + pGroupStats->minEpoch) / 2.0) - GSC_EQUINOX;
          double testFactor;
          double testNearbydec = pNearbyCatalog->decGoodMag + ((pNearbyCatalog->DecPM * testepoch) / (3600.0 * 1000.0));
          double testNearbyra;

          if (pNearbyCatalog->decGoodMag < 90.0 && pNearbyCatalog->decGoodMag > -90.0) {
            testFactor = cos(pNearbyCatalog->decGoodMag * DEGREES_TO_RAD);
            testNearbyra = pNearbyCatalog->raGoodMag + (((pNearbyCatalog->RaPM * testepoch) / (3600.0 * 1000.0)) / testFactor);
          } else {
            testNearbyra = pNearbyCatalog->raGoodMag;
            testFactor = 1.0;
          }

          tmpNearestDistGoodMag = sqr(pGroupStats->dec - testNearbydec) + sqr(testFactor * (pGroupStats->ra - testNearbyra));
        }

        if (minNearestIndexGoodMag < 0) {
          minNearestIndexGoodMag = nearestIndex;
          minNearestDistGoodMag = tmpNearestDistGoodMag;
        } else if (tmpNearestDistGoodMag < minNearestDistGoodMag) {
          minNearestIndexGoodMag = nearestIndex;
          minNearestDistGoodMag = tmpNearestDistGoodMag;
        }
      }
    }

    // Log results from the all-stars search

    pNearestCatalogStar = &pGroupStats->nearestCatalogStar;
    tmpdist = sqrt(minNearestDist);

    if (minNearestIndex >= 0 && tmpdist <= BIN_EXPANSION_RADIUS) {
      if (minNearestIndex >= pFileCommon->statsAllocCount + 3) {
        printf("ERROR: line %d photometryutils curGscBinIndex %d\n", __LINE__, curGscBinIndex);
        exit(1);
      }

      pNearbyCatalog = &pFileCommon->nearby_catalog_table[minNearestIndex];
      GetREF(pNearbyCatalog->REFNumber, pNearestCatalogStar->nearestREF, 1, 1);
      pNearestCatalogStar->nearestREFarcsec = tmpdist * 3600.;
      pNearestCatalogStar->nearestREFmag = pNearbyCatalog->magcal_magdep;
      pNearestCatalogStar->nearestREFflag = 0;
    } else {
      pNearestCatalogStar->nearestREFflag = -1;
    }

    // Log results from the "good mag" search: "stars with magnitudes within 1
    // mag of the object of interest"

    tmpdist = sqrt(minNearestDistGoodMag);

    if (minNearestIndexGoodMag >= 0 && tmpdist <= BIN_EXPANSION_RADIUS) {
      pNearbyCatalog = &pFileCommon->nearby_catalog_table[minNearestIndexGoodMag];

      if (minNearestIndexGoodMag >= pFileCommon->statsAllocCount + 3) {
        printf("ERROR: line %d photometryutils curGscBinIndex %d\n", __LINE__, curGscBinIndex);
        exit(1);
      }

      GetREF(pNearbyCatalog->REFNumber, pNearestCatalogStar->nearestREFGoodMag, 1, 1);
      pNearestCatalogStar->nearestREFarcsecGoodMag = tmpdist * 3600.;
      pNearestCatalogStar->nearestREFmagGoodMag = pNearbyCatalog->magcal_magdepGoodMag;
      pNearestCatalogStar->nearestREFflagGoodMag = 1;
    } else {
      pNearestCatalogStar->nearestREFflagGoodMag = -1;
    }

    // Find ... the nearest galaxy from the galaxy database?

    FindNearestGalaxy(
      pGalaxyCommon,
      pGroupStats->ra,
      pGroupStats->dec,
      pNearestCatalogStar,
      pGroupStats->nearbyObjects,
      NULL
    );
  }

  // Fill out histogram of blend group sizes

  for (groupIndex = 0; groupIndex < pFileCommon->numGroups; groupIndex++) {
    pGroupStats = &pFileCommon->group_stats_table[groupIndex];

    if (pGroupStats->groupNumber < 0) {
      continue;
    }

    if (pGroupStats->gsc_bin_index == curGscBinIndex) {
      if (pGroupStats->groupCount < MAX_NONE_HISTOGRAM) {
        histogramTable[pGroupStats->groupCount]++;
      } else {
        histogramTable[MAX_NONE_HISTOGRAM - 1]++;
      }
    }
  }

  // Compute other statistics about the blend groups

  for (groupIndex = 0; groupIndex < pFileCommon->numGroups; groupIndex++) {
    pGroupStats = &pFileCommon->group_stats_table[groupIndex];

    if (pGroupStats->groupNumber < 0) {
      continue;
    }

    if (pGroupStats->groupCount > 0 && pGroupStats->gsc_bin_index == curGscBinIndex) {
      SaveGroupStatistics(pFileCommon, pGroupStats, pNoneMagnitudeTable, magnitudeCount);
    }
  }

  // Pass 6
  //
  // Fill in identifying information for stars newly matched with catalog sources.

  for (magnitudeIndex = 0; magnitudeIndex < magnitudeCount; magnitudeIndex++ ) {
    pCurStarImage = &pNoneMagnitudeTable[magnitudeIndex];
    pNoneStats = &pFileCommon->none_stats_table[magnitudeIndex];

    strcpy(series, GetSeriesString(pCurStarImage->seriesId, 1));

    if (pNoneStats->selected < magnitudeCount) {
      if (pNoneStats->selected < 0) {
        printf("ERROR: Bad pNoneStats->selected flag in line %d for index %d curGscBinIndex %d\n",__LINE__,magnitudeIndex,curGscBinIndex);
        exit(1);
      }

      if (pNoneStats->groupNumber < 0 || pNoneStats->groupNumber >= pFileCommon->numGroups) {
        if (pFileCommon->enableRematch == 0) {
          printf("ERROR: Bad pNoneStats->groupNumber line %d\n",__LINE__);
          exit(1);
        }
      }

      pGroupStats = &pFileCommon->group_stats_table[pNoneStats->groupNumber];

      if (pGroupStats->groupCount > 0 && pGroupStats->gsc_bin_index == curGscBinIndex) {
        GetREF(pNoneStats->REFNumber, REF, 1, 1);
        pCurStarImage->REFNumber = pNoneStats->REFNumber;
        pCurStarImage->gsc_bin_index = curGscBinIndex;
      }
    }
  }

  // If we're asked to assemble a "target table", do so.

  if (pTarget_table != NULL && pTarget_nrecs != NULL && pTarget_alloc != NULL) {
    for (groupIndex = 0; groupIndex < pFileCommon->numGroups; groupIndex++) {
      pGroupStats = &pFileCommon->group_stats_table[groupIndex];

      if (pGroupStats->groupNumber < 0) {
        continue;
      }

      if (pGroupStats->groupCount > 0 && pGroupStats->gsc_bin_index == curGscBinIndex) {
        magnitudeIndex = pGroupStats->groupIndex;
        pCurStarImage = &pNoneMagnitudeTable[magnitudeIndex];
        pNoneStats = &pFileCommon->none_stats_table[magnitudeIndex];

        if (pCurStarImage->gsc_bin_index != pGroupStats->gsc_bin_index) {
          printf("ERROR: line %d gsc_bin_index %d %d mismatch\n", __LINE__, pCurStarImage->gsc_bin_index,pGroupStats->gsc_bin_index);
          exit(1);
        }

        GetREF(pNoneStats->REFNumber, REF, 1, 1);

        /* Add this one to the target table */

        if (target_nrecs >= target_alloc) {
          target_alloc += 100;
          tmp_target_table = realloc(target_table,target_alloc*sizeof(PHOTTARGET));

          if (tmp_target_table == NULL) {
            printf("ERROR: failed to reallocate target_table of size %zu\n", target_alloc * sizeof(PHOTTARGET));
          }

          target_table = tmp_target_table;
          tmp_target_table = NULL;
        }

        pTarget = &target_table[target_nrecs];
        memset(pTarget,0,sizeof(PHOTTARGET));

        pTarget->groupCount = pGroupStats->groupCount;
        memcpy(&pTarget->nearestCatalogStar, &pGroupStats->nearestCatalogStar, sizeof(NEARESTCATALOGSTAR));
        strcpy(pTarget->nearbyObjects, pGroupStats->nearbyObjects);
        pTarget->magcal_magdep = pGroupStats->magcal_magdep;
        pTarget->Stdmag = 99.0;
        pTarget->color  = 99.0;

        pTarget->updateflag = UPDATEFLAG_NO_UPDATE;
        pTarget->ra = pGroupStats->ra;
        pTarget->dec = pGroupStats->dec;
        pTarget->haveLocation = 1;
        pTarget->gsc_bin_index = pCurStarImage->gsc_bin_index;
        pTarget->coverage_bin_index = GetGSCBin(pGscBin01, pTarget->ra, pTarget->dec, &decCBin, &raCBin, REF);
        strcpy(pTarget->REF, REF);
        strcpy(pTarget->src_name, REF);
        target_nrecs++;
      }
    }
  }

  // All done! Store outputs.

  *pGroupNumber = initialGroupNumber + pFileCommon->numGroups;

  if (pTarget_table != NULL && pTarget_nrecs != NULL && pTarget_alloc != NULL) {
    *pTarget_table = target_table;
    *pTarget_nrecs = target_nrecs;
    *pTarget_alloc = target_alloc;
  }
}



/* Read in the sky2k catalog */
void
GetSky2kCatalog(
  FILE *consoleHandle,
  PSKY2KSTAR *pSky2KTable,
  int *pSky2k_size,
  double *pMaxGSCRadius
) {
  struct StarCat *sc;	/* Star catalog data structure */
  struct Star *star;
  int istar;
  int istar1;
  int istar2;
  double magb;
  double magv;
  double magph;
  double magpv;
#if 0
  FILE *fileHandle;
#endif
  PSKY2KSTAR Sky2KTable;
  PSKY2KSTAR pSky2KStar;
  int sky2k_size = 0;
  int noMagCount = 0;
  int dimStarCount = 0;
  double lnRadius;
  double maxGSCRadius = 0.0;
  *pSky2k_size = 0;
  *pMaxGSCRadius = 0.0;

  sc = binopen("sky2kra");
#if 0
  fileHandle = fopen("/dasch/junk/los.db","wt");
  fprintf(fileHandle,"index\tra\tdec\tmagb\tmagv\tmagph\tmagpv\n");
  fprintf(fileHandle,"-----\t--\t---\t----\t----\t----\t----\n");
#endif
  istar1 = sc->star1;
  istar2 = sc->star0 + sc->nstars;
  Sky2KTable = (PSKY2KSTAR)calloc(sc->nstars,sizeof(SKY2KSTAR));
  if (Sky2KTable == NULL) {
    fprintf(consoleHandle,"ERROR: Failed to allocate Sky2KTable of size %d\n",sc->nstars);
    exit(1);
  }


  star = (struct Star *) calloc (1, sizeof (struct Star));
  star->num = 0.0;

  /* create a k-d tree for 3-dimensional points */
  ptree = kd_create( 3 );

  for (istar = istar1; istar <= istar2; istar++) {
    pSky2KStar = &Sky2KTable[sky2k_size];
    memset(pSky2KStar,0,sizeof(SKY2KSTAR));
    if (binstar (sc, star, istar)) {
      fprintf (stderr,"BINREAD: Cannot read star %d\n", istar);
      break;
    }

    magb = 0.01 * star->mag[0];
    magv = 0.01 * star->mag[1];
    magph = 0.01 * star->mag[2];
    magpv = 0.01 * star->mag[3];
    if (magb != 0.0) {
      pSky2KStar->magnitude = magb;
    } else if (magph != 0.0) {
      pSky2KStar->magnitude = magph;
    } else if (magv != 0.0) {
      pSky2KStar->magnitude = magv;
    } else if (magpv != 0.0) {
      pSky2KStar->magnitude = magpv;
    } else {
      noMagCount++;
      continue;
    }
    if (pSky2KStar->magnitude > 12) {
      dimStarCount++;
      continue;
    }
    pSky2KStar->ra = star->ra;
    pSky2KStar->dec = star->dec;
    lnRadius = -0.9269 - (0.3793 * pSky2KStar->magnitude);
    pSky2KStar->radius = exp(lnRadius);
    if ( pSky2KStar->radius > maxGSCRadius) {
      maxGSCRadius = pSky2KStar->radius;
    }

    /* Insert this node in the kdtree */
    {
      double x;
      double y;
      double z;
      double cosx;
      z = sin(DEGREES_TO_RAD*pSky2KStar->dec);
      cosx = cos(DEGREES_TO_RAD*pSky2KStar->dec);
      x = cosx*cos(DEGREES_TO_RAD*pSky2KStar->ra);
      y = cosx*sin(DEGREES_TO_RAD*pSky2KStar->ra);
      if (kd_insert3(ptree,x,y,z,pSky2KStar) != 0) {
        fprintf(consoleHandle,"ERROR: fatal return from kd_insert3\n");
        exit(1);
      }

    }

    sky2k_size++;
  }
  fprintf(consoleHandle,"Read %d stars from the Sky2K catalog, wrote %d stars, noMagCount %d dimStarCount %d maxGSCRadius %f\n",
          sc->nstars,sky2k_size,noMagCount,dimStarCount,maxGSCRadius);
  free(star);
  binclose(sc);
  *pSky2KTable = Sky2KTable;
  *pSky2k_size = sky2k_size;
  *pMaxGSCRadius = maxGSCRadius;
  return;
}


/* Starting with the quality bits from the photometry database, add what was formerly only for web display */
void
GetFullQuality(PFILESTARIMAGE pFileStarImage, int oldquality, int *pQuality)
{
  int quality = oldquality;
  if (oldquality != QUALITY_UNINITIALIZED) {
    if ((pFileStarImage->BFLAGS & (1<<FILTER_BFLAG_COLORTERM)) != 0) {
      quality |= QUALITY_COLORTERM;
    }
    if (GetFittedPlateScale(pFileStarImage->seriesId,pFileStarImage->plateNumber) >= PATROL_PLATE_SCALE) {
      quality |= QUALITY_PATROL;
    } else {
      quality |= QUALITY_NONPATROL;
    }
    if (pFileStarImage->ELLIPTICITY > MAX_ELLIPTICITY) {
      quality |= QUALITY_TRAILED;
    }
    if ((pFileStarImage->BFLAGS & (1 << FILTER_BFLAG_PSFSATURATED)) != 0) {
      quality |=  QUALITY_SATURATED;
    }
    if (((pFileStarImage->BFLAGS & (1 << FILTER_BFLAG_MAGDEP_MAGCOR)) == 0) &&
        ((quality & QUALITY_UNDETECTED) == 0)) {
      quality |=  QUALITY_NOMAGDEP;
    }
  }
  *pQuality = quality;
  return;
}


void
WriteStarbaseRecord(
  PFILESTARIMAGE pFileStarImage,
  PPHOTPLATES pPhotPlates,
  FILE *dbhandle,
  char *title,
  int *pWriteHeader,
  int catalogNumber
) {
  int quality;
  int tmpCatalogNumber;

  if (*pWriteHeader != 0) {
    if (title != NULL) {
      fprintf(dbhandle,"title %s\n",title);
    }
    fprintf(dbhandle,"REFNumber\tX_IMAGE\tY_IMAGE\tMAG_ISO\tra\tdec\tmagcal_iso\tmagcal_iso_rms\tmagcal_local\tmagcal_local_rms\tDate\tFLUX_ISO\tMAG_APER\tMAG_AUTO\tKRON_RADIUS\tBACKGROUND\tFLUX_MAX\tTHETA_J2000\tELLIPTICITY\tISOAREA_WORLD\tFWHM_IMAGE\tFWHM_WORLD\tplate_dist\tBlendedmag\tlimiting_mag_local\tmagcal_local_error\tdradRMS2\tmagcor_local\textinction\tgsc_bin_index\tseries\tplateNumber\tNUMBER\tversionId\tAFLAGS\tBFLAGS\tISO0\tISO1\tISO2\tISO3\tISO4\tISO5\tISO6\tISO7\tnpoints_local\trejectFlag\tlocal_bin_index\tseriesId\texposureNumber\tsolutionNumber\tspatial_bin\tmosaicNumber\tquality\tplateVersionId\tmagdep_bin\tmagcal_magdep\tmagcal_magdep_rms\tra_2\tdec_2\tRaPM\tDecPM\tA2FLAGS\tB2FLAGS\ttimeAccuracy\tmaskIndex\tcatalogNumber\n");
    fprintf(dbhandle,"---------\t-------\t-------\t-------\t--\t---\t----------\t--------------\t------------\t----------------\t----\t--------\t--------\t--------\t-----------\t----------\t--------\t-----------\t-----------\t-------------\t----------\t----------\t----------\t----------\t------------------\t------------------\t--------\t------------\t----------\t-------------\t------\t-----------\t------\t---------\t------\t------\t----\t----\t----\t----\t----\t----\t----\t----\t-------------\t----------\t---------------\t--------\t--------------\t--------------\t-----------\t------------\t-------\t--------------\t----------\t-------------\t-----------------\t----\t-----\t----\t-----\t-------\t-------\t-----------\t----------\t-------------\n");


    *pWriteHeader = 0;
  }
  GetFullQuality(pFileStarImage,pPhotPlates->quality,&quality);
#if 0
  if (quality != QUALITY_UNINITIALIZED) {
    printf("ERROR: quality is %d in WriteStarbaseRecord()\n",quality);
  }
#endif
  tmpCatalogNumber = catalogNumber;
  if ((RELEASE_EXPERIMENTAL != 0) && (catalogNumber == CATALOG_EXPERIMENTAL)) {
    tmpCatalogNumber = pFileStarImage->catalogNumber;
  } else {
    tmpCatalogNumber = catalogNumber;
  }

  fprintf(dbhandle,"%lld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%d\t%d\n",
          pFileStarImage->REFNumber,
          pFileStarImage->X_IMAGE,
          pFileStarImage->Y_IMAGE,
          pFileStarImage->MAG_ISO,
          pFileStarImage->ra,
          pFileStarImage->dec,
          pFileStarImage->magcal_iso,
          pFileStarImage->magcal_iso_rms,
          pFileStarImage->magcal_local,
          pFileStarImage->magcal_local_rms,
          pFileStarImage->Date,
          pFileStarImage->FLUX_ISO,
          pFileStarImage->MAG_APER,
          pFileStarImage->MAG_AUTO,
          pFileStarImage->KRON_RADIUS,
          pFileStarImage->BACKGROUND,
          pFileStarImage->FLUX_MAX,
          pFileStarImage->THETA_J2000,
          pFileStarImage->ELLIPTICITY,
          pFileStarImage->ISOAREA_WORLD,
          pFileStarImage->FWHM_IMAGE,
          pFileStarImage->FWHM_WORLD,
          pFileStarImage->plate_dist,
          pFileStarImage->Blendedmag,
          pFileStarImage->limiting_mag_local,
          pFileStarImage->magcal_local_error,
          pFileStarImage->dradRMS2,
          pFileStarImage->magcor_local,
          pFileStarImage->extinction,
          pFileStarImage->gsc_bin_index,
          GetSeriesString(pFileStarImage->seriesId,0),
          pFileStarImage->plateNumber,
          pFileStarImage->NUMBER,
          pFileStarImage->versionId,
          pFileStarImage->AFLAGS,
          pFileStarImage->BFLAGS,
          pFileStarImage->ISO0,
          pFileStarImage->ISO1,
          pFileStarImage->ISO2,
          pFileStarImage->ISO3,
          pFileStarImage->ISO4,
          pFileStarImage->ISO5,
          pFileStarImage->ISO6,
          pFileStarImage->ISO7,
          pFileStarImage->npoints_local,
          pFileStarImage->rejectFlag,
          pFileStarImage->local_bin_index,
          pFileStarImage->seriesId,
          pFileStarImage->exposureNumber,
          pFileStarImage->solutionNumber,
          pFileStarImage->spatial_bin,
          pPhotPlates->mosaicNumber,
          quality,
          pPhotPlates->versionId,
          pFileStarImage->magdep_bin,
          pFileStarImage->magcal_magdep,
          pFileStarImage->magcal_magdep_rms,
          pFileStarImage->ra_2,
          pFileStarImage->dec_2,
          pFileStarImage->RaPM,
          pFileStarImage->DecPM,
          pFileStarImage->A2FLAGS,
          pFileStarImage->B2FLAGS,
          pFileStarImage->timeAccuracy,
          pFileStarImage->maskIndex,
          tmpCatalogNumber);
}


void
ChangeGroup(char *filename)
{
  char cmdstr[MAX_FILENAME];
  int result;
  char *dasch_gidname;
  struct stat filestats;
  static gid_t gid = 0;
  struct group *group;

  if (gid == 0) {

    dasch_gidname = getenv("DASCH_GID");
    if (dasch_gidname == NULL) {
      printf("ERROR: DASCH_GID is not defined\n");
      exit(1);
    }
    group = getgrnam(dasch_gidname);
    if (group == NULL) {
      printf("ERROR: A group does not exist for %s\n",dasch_gidname);
      exit(1);
    }
    if (group->gr_gid == 0) {
      printf("ERROR: gid is 0 for group %s\n",dasch_gidname);
      exit(1);
    }
    gid = group->gr_gid;
  }
  result = stat(filename,&filestats);
  if (result != 0) {
    printf("ERROR: ChangeGroup file %s does not exist\n",filename);
    exit(1);
  }
  if (filestats.st_gid != gid) {
    sprintf(cmdstr,"chgrp scanner  %s",filename);
    result = system(cmdstr);
    if (result != 0) {
      printf("ERROR: failed to change group %s rom gid %d\n",cmdstr,filestats.st_gid);
      exit(1);
    }
  }
}


int
GetSequesteredFlag(int seriesId)
{
  if (seriesTableXInit != 2) {
    printf("ERROR: Series table not initialized for plate scale with a call to InitSeriesTable\n");
    exit(1);
  }

  if ((seriesId <= 0) || (seriesId >= MAX_SERIES)) {
    printf("ERROR: Illegal series id %d in GetSequesteredFlag\n",seriesId);
    exit(1);
  }

  return(seriesTableX[seriesId].sequestered);
}


PSERIESENTRY
GetSeriesEntry(int seriesId)
{
  if (seriesTableXInit != 2) {
    printf("ERROR: Series table not initialized for plate scale with a call to InitSeriesTable\n");
    exit(1);
  }

  if ((seriesId <= 0) || (seriesId >= MAX_SERIES)) {
    printf("ERROR: Illegal series id %d in GetSequesteredFlag\n",seriesId);
    exit(1);
  }

  return(&seriesTableX[seriesId]);
}


static int
LightcurveImageCompare(const void *first, const void *second)
{
  double numberFirst = ((PMALMQUIST)first)->limiting_mag_local;
  double numberSecond = ((PMALMQUIST)second)->limiting_mag_local;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}


static int
LightcurveDateCompare(const void *first, const void *second)
{
  double dateFirst = ((PMALMQUIST)first)->Date;
  double dateSecond = ((PMALMQUIST)second)->Date;
  if (dateFirst > dateSecond) {
    return(1);
  } else if (dateFirst < dateSecond) {
    return(-1);
  } else {
    return(0);
  }

}


static int
LightcurveMagnitudeCompare(const void *first, const void *second)
{
  double magnitudeFirst = ((PMALMQUIST)first)->magcal_magdep;
  double magnitudeSecond = ((PMALMQUIST)second)->magcal_magdep;
  if (magnitudeFirst > magnitudeSecond) {
    return(1);
  } else if (magnitudeFirst < magnitudeSecond) {
    return(-1);
  } else {
    return(0);
  }

}


static int
rcmp(double x, double y)
{
  if (x < y)		return -1;
  if (x > y)		return 1;
  return 0;
}


static void
rPsort2(double *x, int lo, int hi, int k)
{
  double v, w;
  int L, R, i, j;

  for (L = lo, R = hi; L < R; ) {
    v = x[k];
    for(i = L, j = R; i <= j;) {
      while (rcmp(x[i], v) < 0) i++;
      while (rcmp(v, x[j]) < 0) j--;
      if (i <= j) { w = x[i]; x[i++] = x[j]; x[j--] = w; }
    }
    if (j < k) L = i;
    if (k < i) R = j;
  }
}


static void
rPsort(double *x, int n, int k)
{
  rPsort2(x, 0, n-1, k);
}


static double
fcube(double x)
{
  return x * x * x;
}


static double
fsquare(double x)
{
  return x * x;
}


static void
calclimits(
  double *x,
  int span,
  int n,
  int loopIndex,
  double mx,
  int *pMinidx,
  int* pMaxidx,
  double *pDmax
) {
  int minidx;
  int maxidx;
  double dmax;
  *pMinidx = 0;
  *pMaxidx = 0;
  *pDmax = 0.0;

  minidx = loopIndex;
  maxidx = loopIndex;

  while((maxidx - minidx) < (span-1)) {
    if (minidx <= 0) {
      /* Can only bump upper loopIndex */
      if (maxidx >= (n-1)) {
        break; /* can't go any further */
      }
      maxidx++;
    } else if (maxidx >= (n-1)) {
      minidx--;
    } else {
      /* Can go either way.  Pick the smallest increment */
      if ((x[maxidx+1] - mx) < (mx - x[minidx-1])) {
        maxidx++;
      } else {
        minidx--;
      }
    }
  }
  dmax = x[maxidx] - mx;
  if ((mx - x[minidx]) > dmax) {
    dmax = (mx - x[minidx]);
  }
  while (minidx > 0) {
    if ((mx - x[minidx-1]) <= dmax) {
      minidx--;
    } else {
      break;
    }
  }
  while (maxidx < (n - 1)) {
    if ((x[maxidx+1] - mx) <= dmax) {
      maxidx++;
    } else {
      break;
    }
  }
  *pMinidx = minidx;
  *pMaxidx = maxidx;
  *pDmax = dmax;

  return;
}


/*
 *  Least mean square fit to minimize the sum of
 *    (A*x) + (B*y) - C;
 *  solve for x, the weighted "1" coefficient
 *
 */
static int
calcleastmeansquare1(
  double* A,
  double * B,
  double * C,
  double *x,
  double * y,
  int minidx,
  int maxidx,
  int loopIndex
) {
  int index;
  /* int N = maxidx-minidx+1; */
  double sqrA = 0.0;
  double sqrB = 0.0;
  double AB = 0.0;
  double AC = 0.0;
  double BC = 0.0;
  double D,E,F,G,H,I;
  double numerator;
  double denominator;
  double XVAL;
  double YVAL;
  *x = 0.0;
  *y = 0.0;

  for (index = minidx; index <= maxidx; index++) {
    sqrA += fsquare(A[index]);
    sqrB +=  fsquare(B[index]);
    AB += A[index]*B[index];
    AC += A[index]*C[index];
    BC += B[index]*C[index];
  }
  /*  We have (A*A*x) + (A*B*y) = AC;   = Dx+Ey=F;
      (A*B*x) + (B*B*y) = BC;   = Gx+Hy=I;
  */

  D = sqrA;
  E = AB;
  F = AC;
  G = AB;
  H = sqrB;
  I = BC;


  numerator = (F*H) - (I*E);
  denominator = (H*D) - (E*G);

  if (denominator == 0.0) {
    if (sqrA == 0.0) {
      return(1); /* Octave also returns a zero for this matrix inversion */
    } else {
      XVAL = AC/sqrA;
      *x = XVAL;
#if 0
      printf("XVAL is %f\n",XVAL);
      error("Return 0 test");
#endif
      return(1);
    }

  }

  XVAL = numerator/denominator;

  YVAL = (F - (D*XVAL))/E;

  *x = XVAL;
  *y = YVAL;
  return(1);

}


/*
 *  Least mean square fit to minimize the sum of
 *    (A*x) + (B*y) + (C*z) - D;
 *  solve for x, the weighted "1" coefficient
 *
 */
static int
calcleastmeansquare2(
  double* A,
  double * B,
  double * C,
  double *D,
  double *x,
  double * y,
  double *z,
  int minidx,
  int maxidx,
  int loopIndex
) {
  int index;
  /* int N = maxidx-minidx+1; */
  double sqrA = 0.0;
  double sqrB = 0.0;
  double sqrC = 0.0;
  double AB = 0.0;
  double AC = 0.0;
  double AD = 0.0;
  double BC = 0.0;
  double BD = 0.0;
  double CD = 0.0;
  double E,F,G,H,I;
  double J,K,L,M,P,Q,R,S,T,U,V,W;
  double D2;
  double numerator;
  double denominator;
  double XVAL;
  double YVAL;
  double ZVAL;
  *x = 0.0;
  *y = 0.0;
  *z = 0.0;

  for (index = minidx; index <= maxidx; index++) {

    sqrA += fsquare(A[index]);
    sqrB += fsquare(B[index]);
    sqrC += fsquare(C[index]);
    AB += A[index]*B[index];
    AC += A[index]*C[index];
    AD += A[index]*D[index];
    BC += B[index]*C[index];
    BD += B[index]*D[index];
    CD += C[index]*D[index];

  }
  /*  We have (A*A*x) + (A*B*y) + (A*C*z) = AD;   = Jx+Ky+Lz=M;
      (A*B*x) + (B*B*y) + (B*C*z) = BD;   = Px+Qy+Rz=S;
      (A*C*x) + (B*C*y) + (C*C*z) = CD;   = Tx+Uy+Vz=W;
  */
  J = sqrA;
  K = AB;
  L = AC;
  M = AD;
  P = AB;
  Q = sqrB;
  R = BC;
  S = BD;
  T = AC;
  U = BC;
  V = sqrC;
  W = CD;


  /* Eliminate y
     (Q*J-K*P)x + (Q*L-K*R)z = (Q*M-K*S) = D2x + Ez = F;
     (U*P-Q*T)x + (U*R-Q*V)z = (U*S-Q*W) =  Gx + Hz = I;
  */
  D2 = (Q*J) - (K*P);
  E  = (Q*L) - (K*R);
  F  = (Q*M) - (K*S);
  G  = (U*P) - (Q*T);
  H  = (U*R) - (Q*V);
  I  = (U*S) - (Q*W);


  numerator = (F*H) - (I*E);
  denominator = (H*D2) - (E*G);

  if (denominator == 0.0) {
    if (sqrA == 0.0) {
      return(1); /* Octave also returns a zero for this matrix inversion */
    } else {
      XVAL = AC/sqrA;
      *x = XVAL;
#if 0
      printf("XVAL is %f\n",XVAL);
      error("Return 0 test");
#endif
      return(1);
    }

  }

  XVAL = numerator/denominator;

  ZVAL = (F - (D2*XVAL))/E;

  YVAL = (S - (P*XVAL) - (R*ZVAL))/Q;




  *x = XVAL;
  *y = YVAL;
  *z = ZVAL;
  return(1);

}


static double
calcmedian(double *x, double * y, int n, int nmax)
{
#if 0 /* Bubble sort */
  int index;
  int index1;
  double temp;
  double returnVal;

  /* Sort the values with a bubble sort */
  for (index = 0; index < n; index++) {
    y[index] = x[index];
  }
  for (index = 0; index < n; index++) {
    for (index1 = 0; index1 < (n-1); index1++)
      if (y[index1] > y[index1+1]) {
        temp = y[index1];
        y[index1] = y[index1+1];
        y[index1+1] = temp;
      }
  }
  if ((n & 1) != 0) {
    returnVal = y[n/2];
  } else {
    returnVal = (y[n/2] + y[(n/2)-1])/2.0;
  }

#else
  int index;
  double returnVal;
  int m1;
  int m2;

  m1 = n/2;
  for (index = 0; index < n; index++) {
    y[index] = x[index];
  }
  rPsort(y,n,m1);
  returnVal = y[m1];
  if(n % 2 == 0) {
    m2 = n-m1-1;
    rPsort(y, n, m2);
    returnVal = (returnVal+y[m2])/2.0;
  }
#endif
  return(returnVal);
}


/*
 *   From Philip R. Bevington, "Data Reduction and Error Analysis for the Physical Sciences", McGraw Hill, New York, 1969.
 *   Find least mean square fit of Y = A + B*X (Equation 6.2, p 99)
 *
 */
static int
BevingtonFit(PMALMQUIST pPoints, int nPoints, double *pB, double *pSigmab, double *pA)
{
  int index;
  double sumx = 0.0;
  double sumx2 = 0.0;
  double sumy = 0.0;
  double sumy2 = 0.0;
  double sumxy = 0.0;
  double determinant;
  double A;
  double B;
  double sigma2 = 0.0;
  double sigmab2 = 0.0;

  PMALMQUIST pLinearFit;
  *pB = 99.0;
  *pA = 99.0;
  *pSigmab = 99.0;
  for (index = 0 ;index < nPoints; index++) {
    pLinearFit = &pPoints[index];
#if 0
    printf("Bevington Date %f mag %f\n",pLinearFit->Date,pLinearFit->magcal_magdep);
#endif

    sumx += pLinearFit->Date;
    sumx2 += (pLinearFit->Date) * (pLinearFit->Date);
    sumy += pLinearFit->magcal_magdep;
    sumy2 += (pLinearFit->magcal_magdep) * (pLinearFit->magcal_magdep);
    sumxy += (pLinearFit->Date) * (pLinearFit->magcal_magdep);

  }
  if (nPoints < 2) {
#if 0
    printf("ERROR: only %d points\n",nPoints);
#endif
    return(-1);
  }
  determinant = sumx2*nPoints - sqr(sumx);
  if (determinant == 0.0) {
#if 0
    printf("ERROR: determinant is zero\n");
#endif
    return(-1);
  }
  A = ((sumx2*sumy) - (sumx * sumxy))/determinant;
  B = ((sumxy*nPoints) - (sumx*sumy))/determinant;
#if 0
  printf("A: %g, B: %g",A,B);
#endif
  *pB = 1000.0 * B;
  *pA = 1000.0 * A;

  /* The following is equation 6-23 on p. 115 */


  sigma2 = (sumy2 + (sqr(A)*nPoints) + (sqr(B)*sumx2) - (2.0*A*sumy) - (2.0*B*sumxy) + (2.0*A*B*sumx))/(nPoints-2.0);
  sigmab2 = (sigma2 * nPoints)/determinant;
#if 0
  printf(" determinant %f 1000*sigma2 %f 1000*sigmab2: %f, 100*(sqrt sigmab2) %f\n",determinant,1000.0*sigma2,1000.0*sigmab2,1000.0*sqrt(sigmab2));
#endif
  *pSigmab = 1000.0*sqrt(sigmab2);

  return(0);
}


static int
countAdjacent(
  int nvals, /* Size of the index array */
  int *vector, /* Array of indices */
  int minAdjacent /* Minimum run of adjacent indices */
) {
  int result = 0;
  int curCount = 0;
  int index;
  for (index = 0; index < (nvals-1); index++) {
    if ((vector[index+1] - vector[index]) == 1) {
      curCount++;
    } else {
      if (curCount >= (minAdjacent-1)) {
        result += curCount - minAdjacent + 2;
      }
      curCount = 0;
    }

  }
  if (curCount >= (minAdjacent-1)) {
    result += curCount - minAdjacent + 2;
  }

  return(result);
}


static double
findmean(int nvals,double *vector) {
  double mean = 0;
  double result = 0;
  int index;

  if (nvals <= 1) {
    return(0);
  }
  for (index = 0; index < nvals; index++) {
    mean += vector[index];
  }
  result = mean/nvals;
  return(result);
}


static double
findstd(int nvals,double *vector) {
  double result = 0;
  int index;
  double meanValue = 0;
  if (nvals <= 1) {
    return(0);
  }
  meanValue = findmean(nvals,vector);

  for (index = 0; index < nvals; index++) {
    result += sqr(vector[index] - meanValue);
  }
  result = result * 1.0/((1.0*nvals) - 1.0);
  result = sqrt(result);
  return(result);
}

static double
corr2(int nvals,double *vectorA,double *vectorB) {
  double meanA = 0;
  double meanB = 0;
  double prodAB = 0;
  double sqrA = 0;
  double sqrB = 0;
  double result;
  int index;
  if (nvals <= 0) {
    return(0);
  }
  for (index = 0; index < nvals; index++) {
    meanA += vectorA[index];
    meanB += vectorB[index];

  }
  meanA = meanA/nvals;
  meanB = meanB/nvals;
  for (index = 0; index < nvals; index++) {
    prodAB += (vectorA[index] - meanA) * (vectorB[index] - meanB);
    sqrA += sqr(vectorA[index] - meanA);
    sqrB += sqr(vectorB[index] - meanB);
  }
  if ((sqrA == 0)  || (sqrB == 0)) {
    return(0);
  }
  result = prodAB/sqrt(sqrA*sqrB);
  return(result);

}

#define EPS 2.2204e-16

static void
clowess(
  double *x,
  double *y,
  double *c,
  int n,
  int robust,
  int span,
  int iter,
  int loessFlag
) {
  double seps = sqrt(EPS);
  int k;
  double mx;
  int loopIndex;
  int index;
  double dmax;
  int minidx;
  int maxidx;
  int weightFlag;
  double * weight = NULL;
  double * v1 = NULL;
  double * v2 = NULL;
  double * v3 = NULL;
  double * y1 = NULL;
  double * r = NULL;
  double * r1 = NULL;
  double * r2 = NULL;
  double mad;
  double minmad;
  double * rweight = NULL;
  int result = 1;
  double xval;
  double yval;
  double zval = 0;
#ifdef DEBUG_TIME
  time_t         startTime;
  time_t         endTime;
#endif /* DEBUG_TIME */


#ifdef DEBUG_TIME
  time(&startTime);
#endif /* DEBUG_TIME */

  weight = (double *)calloc(n,sizeof(double));
  v1 = (double *)calloc(n,sizeof(double));
  v2 = (double *)calloc(n,sizeof(double));
  v3 = (double *)calloc(n,sizeof(double));
  y1 = (double *)calloc(n,sizeof(double));
  if (robust) {
    r = (double *)calloc(n,sizeof(double));
    r1 = (double *)calloc(n,sizeof(double));
    r2 = (double *)calloc(n,sizeof(double));
    rweight = (double *)calloc(n,sizeof(double));
  }

  for (loopIndex = 0; loopIndex < n; loopIndex++) {
    if (loopIndex >= 1) {
      if (x[loopIndex] == x[loopIndex-1]) {
        c[loopIndex] = c[loopIndex-1];
#ifdef SAVE_BOUNDS
        if (robust) {
          lbound[loopIndex] = lbound[loopIndex-1];
          rbound[loopIndex] = rbound[loopIndex-1];
          dmaxv[loopIndex] = dmaxv[loopIndex-1];

        }
#endif /* SAVE_BOUNDS */
        continue;

      }
    }

    mx = x[loopIndex];
    calclimits(x,span,n,loopIndex,mx,&minidx,&maxidx,&dmax);
#ifdef SAVE_BOUNDS
    if (robust) {
      lbound[loopIndex] = minidx;
      rbound[loopIndex] = maxidx;
      dmaxv[loopIndex] = dmax;


    }
#endif /* SAVE_BOUNDS */

    /* Compute weights */
    weightFlag = 1;


    for (index = minidx; index <= maxidx; index++) {
      weight[index] = sqrt(fcube(1.0 - fcube((fabs(x[index]-mx)/dmax))));
      if (weight[index] > seps) {
        weightFlag = 0;
      }
    }
    if (weightFlag == 1) {
      /* All weights are zero, skip weighting */
      for (index = minidx; index <= maxidx; index++) {
        weight[index] = 1.0;
      }
    }

    for (index = minidx; index <= maxidx; index++) {
      v1[index] = 1.0 * weight[index];
      v2[index] = (x[index]-mx) * weight[index];
      if (loessFlag) {
        v3[index] =  (x[index]-mx) * (x[index]-mx) * weight[index];
      }
      y1[index] = y[index] * weight[index];
    }
    if (loessFlag == 0) {
      result = calcleastmeansquare1(v1,v2,y1,&xval,&yval,minidx,maxidx,loopIndex);
    } else {
      result = calcleastmeansquare2(v1,v2,v3,y1,&xval,&yval,&zval,minidx,maxidx,loopIndex);
    }
    if (result == 0) {
      printf("ERROR: result is 0 for loopIndex %d,minidx %d maxidx %d\n",loopIndex,minidx,maxidx);
    }
    c[loopIndex] = xval;



  } /* End of primary loopIndex loop */

  if (robust) {


    minmad = 0.0;
    for (index = 0; index <= n; index++) {
      if (fabs(y[index]) > minmad) {
        minmad = fabs(y[index]);
      }
    }
    minmad = minmad * EPS;

    for (k = 1; k <= iter; k++) {
#ifdef DEBUG_TIME_X
      time(&endTime);
      printf("Starting iteration (2) %3d at: %5d seconds\n",k,(endTime-startTime));
#endif /* DEBUG_TIME_X */


      for (index = 0; index < n; index++) {
        r[index] = y[index]-c[index];
      }
      for (loopIndex = 0; loopIndex < n; loopIndex++) {

        if (loopIndex >= 1) {
          if (x[loopIndex] == x[loopIndex-1]) {
            c[loopIndex] = c[loopIndex-1];
            continue;

          }
        }
        mx = x[loopIndex];
#ifdef SAVE_BOUNDS
        minidx = lbound[loopIndex];
        maxidx = rbound[loopIndex];
        dmax = dmaxv[loopIndex];
#else /* SAVE_BOUNDS */
        calclimits(x,span,n,loopIndex,mx,&minidx,&maxidx,&dmax);
#endif /* SAVE_BOUNDS */
        for (index = minidx; index <= maxidx; index++) {
          y1[index] = y[index];
          r1[index] = r[index];
        }

        weightFlag = 1;
        for (index = minidx; index <= maxidx; index++) {
          weight[index] = sqrt(fcube(1.0 - fcube((fabs(x[index]-mx)/dmax))));
          if (weight[index] > seps) {
            weightFlag = 0;
          }

        }
        if (weightFlag == 1) {
          /* All weights are zero, skip weighting */
          for (index = minidx; index <= maxidx; index++) {
            weight[index] = 1.0;
          }
        }



        mad = calcmedian(&r1[minidx],r2,maxidx-minidx+1,n);
        for (index = minidx; index <= maxidx; index++) {
          r1[index] = fabs(r1[index]-mad);
        }

        mad = calcmedian(&r1[minidx],r2,maxidx-minidx+1,n);
        if (mad > minmad) {
          for (index = minidx; index <= maxidx; index++) {
            rweight[index] = r1[index]/(6.0 * mad);
            if (rweight[index] <= 1.0) {
              rweight[index] = 1.0 - fsquare(rweight[index]);
            } else {
              rweight[index] = 0.0;
            }
            weight[index] = weight[index] * rweight[index];

          }

        }


        for (index = minidx; index <= maxidx; index++) {
          v1[index] = 1.0 * weight[index];
          v2[index] = (x[index]-mx) * weight[index];
          if (loessFlag) {
            v3[index] = (x[index]-mx) * (x[index]-mx) * weight[index];
          }
          y1[index] = y[index] * weight[index];
        }
        if (loessFlag == 0) {
          result = calcleastmeansquare1(v1,v2,y1,&xval,&yval,minidx,maxidx,loopIndex);
        } else {
          result = calcleastmeansquare2(v1,v2,v3,y1,&xval,&yval,&zval,minidx,maxidx,loopIndex);
        }

        if (result == 0) {
          printf("ERROR: result (2) is 0 for loopIndex %d iter %d\n",loopIndex,iter);
        }
        c[loopIndex] = xval;
      } /* End of loopIndex */

    } /* End of iteration liip */

  } /* End of robust iterations */




  free(weight);
  free(v1);
  free(v2);
  free(v3);
  free(y1);
  if (robust) {
#ifdef SAVE_BOUNDS
    free(lbound);
    free(rbound);
    free(dmaxv);
#endif /* SAVE_BOUNDS */
    free(r1);
    free(r2);
    free(rweight);
    free(r);
  }

#ifdef DEBUG_TIME
  time(&endTime);
  printf("Total Time: %5d seconds n: %d robust %d span %d iter %d loessFlag %d\n",(endTime-startTime),n,robust,span,iter,loessFlag);
#endif /* DEBUG_TIME */

}


static void
sgolay(double *x, double *y, double *c, int n, int f, int k)
{
  int HF;
  int index1;
  int index2;
  int tidx;
  int L;
  int R;
  int nspan = 0;
  double * v1 = NULL;
  double * v2 = NULL;
  double * v3 = NULL;
  double * q  = NULL;
  double * y1 = NULL;
  int k2;
  double xval;
  double yval;
  double zval = 0;

  v1 = (double *)calloc(n,sizeof(double));
  v2 = (double *)calloc(n,sizeof(double));
  v3 = (double *)calloc(n,sizeof(double));
  y1 = (double *)calloc(n,sizeof(double));
  q = (double *)calloc(n,sizeof(double));


  if (k != 2) {
    printf("ERROR: sgolay implemented only for degree of 2\n");
    exit(1);
  }
  if (f > n) {
    f = n;
  }
  if ((f & 1) == 0) {
    f -= 1; /* Span must be odd */
  }
  if (f <= k) {
    printf("ERROR: sgolay span %d is less than or equal to degree %d\n",f,k);
    return;
  }

  /* Make sure x is monotonically increasing */
  for (index1 = 0; index1 < (n-1); index1++) {
    if (x[index1+1] < x[index1]) {
      printf("ERROR: sgolay x is not monotonically increasing\n");
      exit(1);
    }

  }
  for (index1 = 0; index1 < n; index1++) {
    c[index1] = y[index1];

  }
  for (index1 = 0; index1 < n; index1++) {

    if (index1 > 0) {
      if (x[index1] == x[index1-1]) {
        c[index1] = c[index1-1];
        continue;
      }
    }
    L = index1;
    R = index1;
    while (1) {
      if (R >= (n-1)) {
        break;
      }
      if (x[R] != x[R+1]) {
        break;
      }
      R = R + 1;
    }
    while (1) {
      if (L <= 0) {
        break;
      }
      if (x[L-1] != x[L]) {
        break;
      }
      L = L - 1;
    }
    HF = ((f - ((R-L+1))) + 1)/2;
    if (HF < 0) {
      HF = 0;
    }
    L = L - HF;
    if (L < 0) {
      L = 0;
    }
    if ((n-f) < L) {
      L = n-f;
    }
    while (1) {
      if (L <= 0) {
        break;
      }
      if (x[L-1] != x[index1]) {  /* BUG? x(index1) should be x(L) */
        break;
      }
      L = L -1;
    }
    R = R + HF;
    if ((L+f-1) > R) {
      R = L+f-1;
    }
    if ((n-1) < R) {
      R = n-1;
    }
    while (1) {
      if (R >= (n-1)) {
        break;
      }
      if (x[R] != x[R+1]) {
        break;
      }
      R = R + 1;
    }
    nspan = R - L + 1;
    for (index2 = 0; index2 < nspan; index2++) {
      tidx = L + index2;
      q[index2] = x[tidx] - x[index1];
      v1[index2] = 1.0;
      v2[index2] = q[index2];
      v3[index2] = sqr(q[index2]);
      y1[index2] = y[tidx];
    }
    k2 = nspan-1;
    if (k2 < 1) {
      k2 = 1;
    }
    if (k2 > k) {
      k2 = k;
    }
    switch (k2) {
    case 1:
      calcleastmeansquare1(v1,v2,y1,&xval,&yval,0,nspan-1,0); /* result ignored */
      break;
    case 2:
      calcleastmeansquare2(v1,v2,v3,y1,&xval,&yval,&zval,0,nspan-1,0); /* result ignored */
      break;
    default:
      printf("ERROR: k2 is %d in sgolay\n",k2);
      xval = y[index1];
    }

    c[index1] = xval;
  }

  free(v1);
  free(v2);
  free(v3);
  free(y1);
  free(q);
  return;
}


static void
smooth(int full_ngood, double *X, double *Y, int span, char *method, double *C)
{
  int robust = 0;
  int iter = 5;
  int loessFlag = 0;
  int index;
  int degree = 2;
  for (index = 0; index < full_ngood; index++) {
    C[index] = Y[index];
  }
  if (method[1] == 'r') {
    robust = 1;
  }
  if (strstr(method,"loess") != NULL) {
    loessFlag = 1;
  }
  if (strstr(method,"sgolay") != NULL) {
    sgolay(X,Y,C,full_ngood,span,degree);
  } else {
    clowess(X,Y,C,full_ngood,robust,span,iter,loessFlag);
  }
}


static void
ReallocParameterVectorStore(PPARAMETERVECTORSTORE pVectorStore)
{
  int index;
  double *tmpVector;
  int *tmpivector;
  PMALMQUIST tmpmvector;

  for (index = 0; index < MAX_PARAMETERVECTOR; index++) {
    tmpVector = realloc(pVectorStore->vector[index], pVectorStore->vectorAlloc * sizeof(double));
    if (tmpVector == NULL) {
      fprintf(
        stderr,
        "ERROR: Failed to reallocate vector %d of size %zu\n",
        index,
        pVectorStore->vectorAlloc * sizeof(double)
      );
      exit(1);
    }

    pVectorStore->vector[index] = tmpVector;
  }

  for (index = 0; index < MAX_MALMQUIST_PARAMETER; index++) {
    tmpmvector = realloc(pVectorStore->mvector[index], pVectorStore->vectorAlloc * sizeof(MALMQUIST));
    if (tmpmvector == NULL) {
      fprintf(
        stderr,
        "ERROR: Failed to reallocate mvector %d of size %zu\n",
        index,
        pVectorStore->vectorAlloc * sizeof(double)
      );
      exit(1);
    }

    pVectorStore->mvector[index] = tmpmvector;
  }

  for (index = 0; index < MAX_IPARAMETERVECTOR; index++) {
    tmpivector = realloc(pVectorStore->ivector[index], pVectorStore->vectorAlloc * sizeof(int));
    if (tmpivector == NULL) {
      fprintf(
        stderr,
        "ERROR: Failed to reallocate ivector %d of size %zu\n",
        index,
        pVectorStore->vectorAlloc * sizeof(double)
      );
      exit(1);
    }

    pVectorStore->ivector[index] = tmpivector;
  }
}


void
FreeParameterVectorStore(PPARAMETERVECTORSTORE pVectorStore)
{
  int index;
  double *vector;
  int * ivector;
  PMALMQUIST mvector;
  for (index = 0; index < MAX_PARAMETERVECTOR; index++) {
    vector = pVectorStore->vector[index];
    if (vector != NULL) {
      free(vector);
    }
  }
  for (index = 0; index < MAX_MALMQUIST_PARAMETER; index++) {
    mvector = pVectorStore->mvector[index];
    if (mvector != NULL) {
      free(mvector);
    }
  }
  for (index = 0; index < MAX_IPARAMETERVECTOR; index++) {
    ivector = pVectorStore->ivector[index];
    if (ivector != NULL) {
      free(ivector);
    }
  }
  memset(pVectorStore,0,sizeof(PARAMETERVECTORSTORE));
}


#define INCLUDE_TRANSIENT_DEFECTS 1
#define INCLUDE_TRANSIENT_MULTIPLE 1
#define INCLUDE_TRANSIENT_ONE_MULTIPLE 1
/* #define  INCLUDE_TRANSIENT_DRADBIN 1 */

/* Search for flares, requested by Josh Wed 4/01/15 11:50 PM and following discussions
 *
 *          Flare amplitude > 1 mag
 *          ObsDur > 3 days < 90 days at least 3 points These measured from the peak and including the peak.
 *          Flare Dur > 0.6 ObsDur
 *          delta-drad < 9 (~max drad=18/2)
 *
 * Work only with Population 6 stars that have been filtered for Pickering Wedge and Multiple Plates
 * 1.  Find every peak
 * 2.  For each peak in order of brightness, take a range of 90 days forward and 90 days back.
 *     If peak is not centered reject (Tail after peak needs to be 60% of the lightcurve)
 *     If < 3 points reject
 *     If < 3 days, reject
 *     Recalculate drad and reject RMS mostion > 9 arcsec.
 *     Check peak vs minimum magnitude to 90 days, must be at least 0.5 mag
 *     Find slope from peak to end - must be a decline of at least 0.25 mag
 *     If good, output duration in peakDays peakYear peakMag, or 0.0, 0.0, and 0.0 and STOP here
 *
 * 3.  On rejection, eliminate considered data points on both sides down to the minimum magnitude.
 *
 */
static void
ProcessTransientCandidates(
  PFILECOMMON pFileCommon,
  PTRANCOMMON pTranCommon,
  PPARAMETERVECTORSTORE pVectorStore,
  PPHOTSTARIMAGE pMagnitudeTable,
  PGALAXYCOMMON pGalaxyCommon,
  PSTARENTRY pCurStarEntry,
  char *nearbyObjects,
  int AFLAGSMASK1,
  int QUALITYMASK
) {
  int index_ngood;
  PMALMQUIST pMalmquist1 = NULL;                    /* Preceeding star image */
  PPHOTSTARIMAGE pSumStarImage1 = NULL;
  PFILESTARIMAGE pFileStarImage1 = NULL;
  PMALMQUIST pMalmquist2 = NULL;              /* Star image of interest */
  PPHOTSTARIMAGE pSumStarImage2 = NULL;
  PFILESTARIMAGE pFileStarImage2 = NULL;
  PMALMQUIST pMalmquist3 = NULL;              /* Following star image */
  PPHOTSTARIMAGE pSumStarImage3 = NULL;
  PFILESTARIMAGE pFileStarImage3 = NULL;
  int peak_count = 0;
  PMALMQUIST pMalmquist4;                    /* temporary placeholder */
  PPHOTSTARIMAGE pSumStarImage4 = NULL;
  PFILESTARIMAGE pFileStarImage4 = NULL;
  PMALMQUIST pMalmquist5 = NULL;              /* temporary file */
  PPHOTSTARIMAGE pSumStarImage5 = NULL;
  PFILESTARIMAGE pFileStarImage5 = NULL;
  int peak_index;
  int index_start;
  int index_peak;
  int index_end;
  int temp_index;
  int rejectFlag;
  double dateRange;
  double min_magcal_magdep = 99.0;
  double min_nodefect_magcal_magdep = 99.0;
  double max_magcal_magdep = -99.0;
  double ave_magcal_magdep = 99.0;
  double ave_magcal_magdep2;
  double count_magcal_magdep2;
  double variancesum;
  int variancecount;
  double fitmag;
  double fitdiff;
  double slope;
  double Sigmab;
  char REF[MAX_REF];
  double dradmed = 99.0;
  double dradrms = 99.0;
  double dradRMS2med = 99.0;
  double dradRMS2rms = 99.0;
  double dradRMS3med = 99.0;
  double dradRMS3rms = 99.0;
  double dradRMS3;
  double minimum_magcal_magdep = 99.0;
  double peakOffset;
  int dradIndex;
  int nearbyCatalogIndex;
  PNEARBYCATALOG pNearbyCatalog;
  double factor;
  double curDistance;
  long long tempREFNumber;
  int refType;
  double nearby_magcal_magdep = 99.;
  double patrolScale;
  double raAverage = 0;
  double decAverage = 0;
  double positionWeight = 0;
  int positionCount = 0;
  double dateyear;
  int intyear;
  PPLATELIMITINGREC pPlateLimitingRec;
  int plateIndex;
  int plateLimitingYears[PIPELINE_MAX_DATE-PIPELINE_MIN_DATE+1];

  pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_ENTRY);

  if (pTranCommon->full_ngood < MIN_TC_POINTS) {
    pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_INSUFFICIENT_POINTS);
    return; /* not enough points in the curve */
  }


  /* This code is effectively disabled when MIN_TC_CLIP_NGOOD is greater than MAX_TC_POINTS as is currently the case */
  if (pTranCommon->clip_ngood > MIN_TC_CLIP_NGOOD) {
    /* We have a real lightcurve here.  No transient can be dimmer then three sigma above the median */
    minimum_magcal_magdep = pTranCommon->clip_med - (3.0*pTranCommon->clip_rms);
  }


  /* Resort the MAMQUIST selections according to date */
  qsort((void*)pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST],pTranCommon->full_ngood,sizeof(MALMQUIST),LightcurveDateCompare);

  for (index_ngood = 0; index_ngood < pTranCommon->full_ngood; index_ngood++) {
    pMalmquist2 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_ngood];
    pSumStarImage2 = &pMagnitudeTable[pMalmquist2->imageIndex];
    pFileStarImage2 = pSumStarImage2->pFileStarImage;
    if ((pMalmquist2->selectionFlags & SELECT_GOOD) == 0) {
      printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
      exit(1);
    }
    if (pMalmquist2->Date != pFileStarImage2->Date) {
      printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
      exit(1);
    }
    if (index_ngood == 0) {
      pMalmquist1 = NULL;    /* Preceeding star image */
      pSumStarImage1 = NULL;
      pFileStarImage1 = NULL;
    } else {
      pMalmquist1 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_ngood-1];
      pSumStarImage1 = &pMagnitudeTable[pMalmquist1->imageIndex];
      pFileStarImage1 = pSumStarImage1->pFileStarImage;
      if ((pFileStarImage2->Date-pFileStarImage1->Date) > MAX_TC_DAYS) {
        pMalmquist1 = NULL;    /* Preceeding star image is too early */
        pSumStarImage1 = NULL;
        pFileStarImage1 = NULL;
      }
    }

    if (index_ngood == (pTranCommon->full_ngood -1)) {
      pMalmquist3 = NULL;              /* Following star image */
      pSumStarImage3 = NULL;
      pFileStarImage3 = NULL;
    } else {
      pMalmquist3 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_ngood+1];
      pSumStarImage3 = &pMagnitudeTable[pMalmquist3->imageIndex];
      pFileStarImage3 = pSumStarImage3->pFileStarImage;
      if (pSumStarImage3 == pSumStarImage2) {
        printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
        exit(1);
      }
      if (pFileStarImage3 == pFileStarImage2) {
        printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
        exit(1);
      }

      if (pFileStarImage3->Date < pFileStarImage2->Date) {
        printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
        exit(1);
      }
      if ((pFileStarImage3->Date-pFileStarImage2->Date) > MAX_TC_DAYS) {
        pMalmquist3 = NULL;    /* Following star image is too late*/
        pSumStarImage3 = NULL;
        pFileStarImage3 = NULL;
      }
    }
    if (pFileStarImage2->magcal_magdep > 90) {
      printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
      exit(1);
    }
    if ((pMalmquist1 == NULL) && (pMalmquist3 == NULL)) {
      /* Only one peak in this section of the curve */
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_PEAKVALUES][peak_count];
      memcpy(pMalmquist4,pMalmquist2,sizeof(MALMQUIST));
      pMalmquist4->secondaryIndex = index_ngood;
      peak_count++;
    } else if (pMalmquist3 == NULL) {
      if (pMalmquist2->magcal_magdep < pMalmquist1->magcal_magdep) {
        /* Peak at the end of the curve */
        pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_PEAKVALUES][peak_count];
        memcpy(pMalmquist4,pMalmquist2,sizeof(MALMQUIST));
        pMalmquist4->secondaryIndex = index_ngood;
        peak_count++;
      }
    } else if (pMalmquist1 == NULL) {
      if (pMalmquist2->magcal_magdep < pMalmquist3->magcal_magdep) {
        /* Peak at the beginning of the curve */
        pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_PEAKVALUES][peak_count];
        memcpy(pMalmquist4,pMalmquist2,sizeof(MALMQUIST));
        pMalmquist4->secondaryIndex = index_ngood;
        peak_count++;
      }

    } else {
      if ((pMalmquist2->magcal_magdep < pMalmquist3->magcal_magdep)  &&
          (pMalmquist2->magcal_magdep < pMalmquist1->magcal_magdep)) {
        pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_PEAKVALUES][peak_count];
        memcpy(pMalmquist4,pMalmquist2,sizeof(MALMQUIST));
        pMalmquist4->secondaryIndex = index_ngood;
        peak_count++;
      }
    }
  }
  if (peak_count == 0) {
    pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_NOPEAKS);
    /* No peaks found here */
    return;
  }

  pTranCommon->peakEventCount3 = peak_count;

  qsort((void*)pVectorStore->mvector[PARAMETERVECTOR_PEAKVALUES],peak_count,sizeof(MALMQUIST),LightcurveMagnitudeCompare);
  for (peak_index = 0; peak_index < peak_count; peak_index++) {
    /* Start with the brightest peak */
    pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_PEAKVALUES][peak_index];

    if (pMalmquist4 == NULL) {
      printf("ERROR: Null pointer in line %d for lightcurve REFNumber %lld\n",__LINE__,pFileStarImage2->REFNumber);
      fflush(stdout);
      exit(1);
    }
    if (pMalmquist4->magcal_magdep > minimum_magcal_magdep) {
      /* We are in the main body of the lightcurve */
      continue;
    }
#if 0
    printf("peak_index %d magnitude %f\n",peak_index,pMalmquist4->magcal_magdep);
#endif
    /* Recover our date-ordered reference */
    pMalmquist2 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pMalmquist4->secondaryIndex];
    if (pMalmquist2 == NULL) {
      printf("ERROR: Null pointer in line %d for lightcurve REFNumber %lld\n",__LINE__,pFileStarImage2->REFNumber);
      fflush(stdout);
      exit(1);
    }
    if ((pMalmquist4->secondaryIndex < 0) ||
        (pMalmquist4->secondaryIndex >= pVectorStore->vectorAlloc)) {
      printf("ERROR: Bad pointer in line %d for lightcurve REFNumber %lld\n",__LINE__,pFileStarImage2->REFNumber);
      fflush(stdout);
      exit(1);
    }
    if ((pMalmquist2->selectionFlags & SELECT_TRANSIENT) != 0) {
      /* We already looked at this part of the lightcurve */
      continue;
    }
    /* Search for the beginning of the transient */
    index_peak = pMalmquist4->secondaryIndex;
    index_start = pMalmquist4->secondaryIndex-1;
    while (index_start >= 0) {
      pMalmquist1 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_start];
      if ((pMalmquist1->selectionFlags & SELECT_TRANSIENT) != 0) {
        /* We already looked at this part of the lightcurve */
        break;
      }
      if ((pMalmquist2->Date - pMalmquist1->Date) > (2*MAX_TC_DAYS)) {
        /* This point is too far in the past */
        break;
      }
      if (pMalmquist2->magcal_magdep > minimum_magcal_magdep) {
        /* We are in the main body of the lightcurve */
        break;
      }
      index_start--;
    }
    index_start++;
    pMalmquist1 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_start];
    /* Search for the end of the transient */
    index_end =  pMalmquist4->secondaryIndex + 1;
    while (index_end < pTranCommon->full_ngood) {
      pMalmquist3 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_end];
      if ((pMalmquist3->selectionFlags & SELECT_TRANSIENT) != 0) {
        /* We already looked at this part of the lightcurve */
        break;
      }
      if ((pMalmquist3->Date - pMalmquist2->Date) > (2*MAX_TC_DAYS)) {
        /* This point is too far in the future */
        break;
      }
      if (pMalmquist3->magcal_magdep > minimum_magcal_magdep) {
        /* We are in the main body of the lightcurve */
        break;
      }

      index_end++;
    }
    index_end--;
    pTranCommon->peakEventCount2++;

    pMalmquist3 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_end];
    rejectFlag = 0; /* After this point, we need to set the SELECT_TRANSIENT flag */
    /* What we have needs to be at least 90 days */
    dateRange = pMalmquist3->Date - pMalmquist1->Date;

    if (((index_end - index_peak + 1) < MIN_TC_POINTS) ||
        (dateRange < MIN_TC_DAYS) ||
        (dateRange > MAX_TC_DAYS)) {
      if (dateRange < MIN_TC_DAYS) {
        pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_TOOSHORTDATERANGE);
      } else if (dateRange > MAX_TC_DAYS) {
        pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_TOOLONGDATERANGE);
      } else {
        pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_INSUFFICIENT_POINTS3);
      }
      rejectFlag = 1;
    }
    if (rejectFlag == 0) {
      dradIndex = 0;
      ave_magcal_magdep2 = 0;
      count_magcal_magdep2 = 0;
      pTranCommon->peakMultipleCount = 0;
      for (temp_index = index_start; temp_index <= index_end; temp_index++) {
        pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
        pSumStarImage4 = &pMagnitudeTable[pMalmquist4->imageIndex];
        if ((pSumStarImage4->quality & QUALITY_MULTIPLE) != 0) {
          pTranCommon->peakMultipleCount++; /* Count multiple exposure plates within the lightcurve */
        }
        ave_magcal_magdep2 += pMalmquist4->magcal_magdep;
        count_magcal_magdep2++;
        if (temp_index == index_start) {
          if ((pMalmquist4->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) {
            min_nodefect_magcal_magdep = pMalmquist4->magcal_magdep;
          } else {
            min_nodefect_magcal_magdep = 99.9;
          }
          min_magcal_magdep = pMalmquist4->magcal_magdep;
          max_magcal_magdep = pMalmquist4->magcal_magdep;
        } else {
          if (min_magcal_magdep > pMalmquist4->magcal_magdep) {
            min_magcal_magdep = pMalmquist4->magcal_magdep;
          }
          if ((pMalmquist4->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) {
            if (min_nodefect_magcal_magdep > pMalmquist4->magcal_magdep) {
              min_nodefect_magcal_magdep = pMalmquist4->magcal_magdep;
            }
          }
          if (max_magcal_magdep < pMalmquist4->magcal_magdep) {
            max_magcal_magdep = pMalmquist4->magcal_magdep;
          }
        }


#if defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
        pVectorStore->vector[PARAMETERVECTOR_DRAD][dradIndex] = pMalmquist4->drad;
        pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = pMalmquist4->dradRMS2;
        dradIndex++;
#endif
#if defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
        if ((pMalmquist4->AFLAGS & (1 << FILTER_AFLAG_BACKGROUND)) == 0) {
          pVectorStore->vector[PARAMETERVECTOR_DRAD][dradIndex] = pMalmquist4->drad;
          pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = pMalmquist4->dradRMS2;
          dradIndex++;
        }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
        if (pMalmquist4->AFLAGS & (1 << FILTER_AFLAG_DEFECT) == 0) { /* change of Jul 13, 2015*/
          pVectorStore->vector[PARAMETERVECTOR_DRAD][dradIndex] = pMalmquist4->drad;
          pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = pMalmquist4->dradRMS2;
          dradIndex++;
        }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
        if ((pMalmquist4->AFLAGS & ((1 << FILTER_AFLAG_DEFECT) | (1 << FILTER_AFLAG_BACKGROUND)) == 0)) { /* change of Jul 13, 2015*/
          pVectorStore->vector[PARAMETERVECTOR_DRAD][dradIndex] = pMalmquist4->drad;
          pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = pMalmquist4->dradRMS2;
          dradIndex++;
        }
#endif



      }
#ifdef INCLUDE_TRANSIENT_ONE_MULTIPLE
      if (pTranCommon->peakMultipleCount > MAX_TC_MULTIPLE_POINTS) {
        pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_MULTIPLE1);

        rejectFlag = 1;
      }
#endif /* INCLUDE_TRANSIENT_ONE_MULTIPLE */
      if (count_magcal_magdep2 > 0) {
        ave_magcal_magdep2 = ave_magcal_magdep2/count_magcal_magdep2;
      } else {
        pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_SOFTWARE);
        rejectFlag = 1; /* should never happen! */
      }
      if (CalcMedianAndRMS(dradIndex,pTranCommon->minGoodStars,pVectorStore->vector[PARAMETERVECTOR_DRAD],&dradmed,&dradrms,0,3.0,1) == 0) {
        dradmed = 99.0;
        dradrms = 99.0;
      }
      if (CalcMedianAndRMS(dradIndex,pTranCommon->minGoodStars,pVectorStore->vector[PARAMETERVECTOR_DRADRMS2],&dradRMS2med,&dradRMS2rms,0,3.0,1) == 0) {
        dradRMS2med = 99.0;
        dradRMS2rms = 99.0;
      }

      if (((max_magcal_magdep - min_magcal_magdep) < MIN_TC_MAGNITUDE)
#if 0
          || (dradrms > MIN_TC_DRAD)
#endif
          ) {
        pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_LOW_AMPLITUDE);
        rejectFlag = 1;
      }
    }

    if (rejectFlag == 0) {
      /* Search the nearby_catalog_table for anything within 90 arcsec */
      factor = cos(DEGREES_TO_RAD*pCurStarEntry->dec);
      for (nearbyCatalogIndex = 0; nearbyCatalogIndex < pFileCommon->nearbyCatalogCount; nearbyCatalogIndex++) {
        pNearbyCatalog = &pFileCommon->nearby_catalog_table[nearbyCatalogIndex];
        if (pNearbyCatalog->REFNumber == pCurStarEntry->REFNumber) {
          continue;
        }
        if ((pNearbyCatalog->dec > 90.0) ||
            (pNearbyCatalog->ra > 990.0)) {
          continue;
        }
        nearby_magcal_magdep = pNearbyCatalog->magcal_magdepGoodMag;
        if (nearby_magcal_magdep > 90) {
          nearby_magcal_magdep = pNearbyCatalog->magcal_magdep;
        }
        if (nearby_magcal_magdep > 90.0) {
          continue;
        }
        curDistance = 3600.*sqrt(sqr(pNearbyCatalog->dec-pCurStarEntry->dec) + sqr(factor*(pNearbyCatalog->ra-pCurStarEntry->ra)));
        if (curDistance < MAX_TC_NEARBY_ARCSEC) {
          if (fabs(nearby_magcal_magdep-min_magcal_magdep) < MAX_TC_NEARBY_MAG) {
            pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_NEARBY_DIMMAG);
            rejectFlag = 1;
            break;
          } else if (fabs(nearby_magcal_magdep-ave_magcal_magdep) < MAX_TC_NEARBY_MAG) {
            pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_NEARBY_AVEMAG);
            rejectFlag = 1;
            break;
          }
        } else {
          if (fabs(nearby_magcal_magdep-min_magcal_magdep) < MAX_TC_NEARBY_MAG) {
            if (pTranCommon->peakNearbyDistance == 0) {
              pTranCommon->peakNearbyDistance = curDistance;
            } else if (curDistance <  pTranCommon->peakNearbyDistance) {
              pTranCommon->peakNearbyDistance = curDistance;
            }
          } else if (fabs(nearby_magcal_magdep-ave_magcal_magdep) < MAX_TC_NEARBY_MAG) {
            if (pTranCommon->peakNearbyDistance == 0) {
              pTranCommon->peakNearbyDistance = curDistance;
            } else if (curDistance <  pTranCommon->peakNearbyDistance) {
              pTranCommon->peakNearbyDistance = curDistance;
            }
          }

        }
      }
    }

    if (rejectFlag == 0) {
      /* This is a catalog star.  Make sure the flare is 2 mag brighter than the catalog magnitude */
      if ((pCurStarEntry->Stdmag < 90.0) &&
          (GetREFNumber(pCurStarEntry->REF,&tempREFNumber,&refType,0,0) == 0) &&
          (refType != REF_TYPE_DASCH)) {
        if ((pCurStarEntry->Stdmag-min_magcal_magdep) < MAX_TC_NEARBY_MAG) { /* fabs removed Jan 6, 2017 */
          pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_NEARBY_DIMMAG2);
          rejectFlag = 1;
        } else if ((pCurStarEntry->Stdmag-ave_magcal_magdep) < MAX_TC_NEARBY_MAG) { /* fabs removed Jan 6, 2017 */
          pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_NEARBY_AVEMAG2);
          rejectFlag = 1;
        }
      }

    }

    if (rejectFlag == 0) {
      /* Now look for nearby objects */
      if (FindNearestTycho2Star(pGalaxyCommon,pCurStarEntry->ra,pCurStarEntry->dec,MAX_TC_BRIGHT_MAG,MAX_TC_BRIGHT_ARCSEC) != 0) {
        pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_NEARTYCHO2);
        rejectFlag = 1;
      }

    }

    if (rejectFlag == 0) {
      if (pTranCommon->peakEventCount == 0) {
        pTranCommon->peakDays = dateRange;
        pTranCommon->peakYear = jd2ep(pMalmquist2->Date);
        pTranCommon->peakMag = max_magcal_magdep - min_magcal_magdep;
        if (BevingtonFit(pMalmquist2,index_end-index_peak+1,&slope,&Sigmab,&peakOffset) == 0) {
          pTranCommon->peakSlope = slope/1000.; /* The factor of 1000 is added by BevingtonFit */
          pTranCommon->peakSlopeRMS = Sigmab/1000.; /* The factor of 1000 is added by BevingtonFit */
        }
        pTranCommon->peakCount = index_end-index_start+1;
        pTranCommon->peakOutside = pTranCommon->full_ngood-pTranCommon->peakCount;
        pTranCommon->peakMaxMag = min_magcal_magdep;
        pTranCommon->peakNoDefectMag = min_nodefect_magcal_magdep;
        pTranCommon->peakMinMag = max_magcal_magdep;
        pTranCommon->peakNumber = peak_index+1;
        pTranCommon->peakEventCount = 1;
        pTranCommon->peakDradRMS = dradrms;
        pTranCommon->peakDradRMS2 = dradRMS2rms;
        /* Count wide and narrow field points for the lightcurve and find the transient weighted position */
        for (index_ngood = 0; index_ngood < pTranCommon->full_ngood; index_ngood++) {
          pMalmquist5 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_ngood];
          pSumStarImage5 = &pMagnitudeTable[pMalmquist5->imageIndex];
          pFileStarImage5 = pSumStarImage5->pFileStarImage;
          patrolScale = GetFittedPlateScale(pFileStarImage5->seriesId,pFileStarImage5->plateNumber);
#if 0
          qualityMask = GetPlateQualityMask(pFileCommon,pFileStarImage5->seriesId,pFileStarImage5->plateNumber,&quality,&pColorterm,&pColorflag,pFileStarImage5->REFNumber);
          if ((quality == 0) &&
              ((qualityMask & (1 << pFileStarImage5->spatial_bin)) != 0)) {
            printf("line %d Bad colorterm for %s plate %s spatial bin %d gsc_bin_index %d\n",__LINE__,pCurStarEntry->REF,pSumStarImage5->Plate,pFileStarImage5->spatial_bin,pFileStarImage5->gsc_bin_index);
          }
#endif
          if ((index_ngood < index_start) || (index_ngood > index_end)) {
            if (patrolScale > PATROL_PLATE_SCALE) {
              pTranCommon->peakOutsideWF++;
            } else {
              pTranCommon->peakOutsideNF++;
            }
          } else {
#if defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if ((pFileStarImage5->dradRMS2 > 0) &&
                (pFileStarImage5->dradRMS2 < 90.0)) {
              raAverage += (pFileStarImage5->ra+360.)/sqr(pFileStarImage5->dradRMS2);
              decAverage += (pFileStarImage5->dec+180.)/sqr(pFileStarImage5->dradRMS2);
              positionWeight += 1.0/sqr(pFileStarImage5->dradRMS2);
              positionCount++;
            }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if (((pFileStarImage5->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0) &&
                (pFileStarImage5->dradRMS2 > 0) &&
                (pFileStarImage5->dradRMS2 < 90.0)) {
              raAverage += (pFileStarImage5->ra+360.)/sqr(pFileStarImage5->dradRMS2);
              decAverage += (pFileStarImage5->dec+180.)/sqr(pFileStarImage5->dradRMS2);
              positionWeight += 1.0/sqr(pFileStarImage5->dradRMS2);
              positionCount++;
            }
#endif
#if defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if (((pFileStarImage5->AFLAGS & (1 << FILTER_AFLAG_BACKGROUND)) == 0) &&
                (pFileStarImage5->dradRMS2 > 0) &&
                (pFileStarImage5->dradRMS2 < 90.0)) {
              raAverage += (pFileStarImage5->ra+360.)/sqr(pFileStarImage5->dradRMS2);
              decAverage += (pFileStarImage5->dec+180.)/sqr(pFileStarImage5->dradRMS2);
              positionWeight += 1.0/sqr(pFileStarImage5->dradRMS2);
              positionCount++;
            }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if (((pFileStarImage5->AFLAGS & ((1 << FILTER_AFLAG_DEFECT) | (1 << FILTER_AFLAG_BACKGROUND))) == 0) &&
                (pFileStarImage5->dradRMS2 > 0) &&
                (pFileStarImage5->dradRMS2 < 90.0)) {
              raAverage += (pFileStarImage5->ra+360.)/sqr(pFileStarImage5->dradRMS2);
              decAverage += (pFileStarImage5->dec+180.)/sqr(pFileStarImage5->dradRMS2);
              positionWeight += 1.0/sqr(pFileStarImage5->dradRMS2);
              positionCount++;
            }
#endif

            if ((pFileStarImage5->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) != 0) {
              pTranCommon->peakDefectCount++; /* Count defects within the lightcurve */
            }
            if (patrolScale > PATROL_PLATE_SCALE) {
              pTranCommon->peakCountWF++;
            } else {
              pTranCommon->peakCountNF++;
            }

            /* Here we count points with uncertain or bad colorterm correction */
            refType = GetREFType(pCurStarEntry->REFNumber);
            if ((refType != REF_TYPE_DASCH) &&
                (pCurStarEntry->Stdmag < 90.0)) {
              /* This is a catalog star.  Flag bad color if the color is bad or the colorterm is bad */
              if ((pCurStarEntry->color > 90.0) ||
                  (pSumStarImage5->colorflag != COLORFLAG_METROPOLIS) ||
                  (pSumStarImage5->colorterm > 90.0)) {
#if 0
                printf("ERROR line %d badcolorflag for REF %s gsc_bin_index %d and plate %s\n",__LINE__,pSumStarImage5->REF,pFileStarImage5->gsc_bin_index,pSumStarImage5->Plate);
#endif
                pSumStarImage5->badcolorflag = 1;
                pTranCommon->peakBadColorCount++;
              }
            } else {
              /* This is an unmatched star or matched to the wrong catalog.  Flag bad color if the colorterm is bad or if
                 the colorterm is not blue */
              if ((pSumStarImage5->colorflag != COLORFLAG_METROPOLIS) ||
                  (pSumStarImage5->colorterm > 90.0) ||
                  (BlueColorterm(pFileCommon->catalogNumber,pSumStarImage5->colorterm) == 0)) {
#if 0
                printf("ERROR line %d badcolorflag for REF %s gsc_bin_index %d and plate %s\n",__LINE__,pSumStarImage5->REF,pFileStarImage5->gsc_bin_index,pSumStarImage5->Plate);
#endif
                pSumStarImage5->badcolorflag = 1;
                pTranCommon->peakBadColorCount++;
              }


            }

          }
        }
        if (positionCount > 1) {
          raAverage = raAverage/positionWeight;
          while (raAverage >= 360.0) {
            raAverage -= 360.0;
          }
          decAverage = decAverage/positionWeight;
          while (decAverage > 90.0) {
            decAverage -= 180.0;
          }
          pTranCommon->peakRA = raAverage;
          pTranCommon->peakDec = decAverage;


          dradIndex = 0;
          for (temp_index = index_start; temp_index <= index_end; temp_index++) {
            pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
            pSumStarImage4 = &pMagnitudeTable[pMalmquist4->imageIndex];
            pFileStarImage4 = pSumStarImage4->pFileStarImage;
            factor =  cos(DEGREES_TO_RAD*(pTranCommon->peakDec));
            dradRMS3 = sqrt(sqr(3600.0*(pFileStarImage4->dec-pTranCommon->peakDec))+sqr(3600*factor*(pFileStarImage4->ra-pTranCommon->peakRA)));
#if defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = dradRMS3;
            dradIndex++;
#endif
#if defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if ((pMalmquist4->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0) {
              pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = dradRMS3;
              dradIndex++;
            }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if ((pMalmquist4->AFLAGS & (1 << FILTER_AFLAG_BACKGROUND)) == 0) {
              pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = dradRMS3;
              dradIndex++;
            }

#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if (pMalmquist4->AFLAGS & ((1 << FILTER_AFLAG_DEFECT) | (1 << FILTER_AFLAG_BACKGROUND)) == 0) {
              pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = dradRMS3;
              dradIndex++;
            }
#endif
          }
          /* Now compute peakDradRMS3 with respect to this transient center */
          if (CalcMedianAndRMS(dradIndex,pTranCommon->minGoodStars,pVectorStore->vector[PARAMETERVECTOR_DRADRMS2],&dradRMS3med,&dradRMS3rms,0,3.0,1) == 0) {
            dradRMS3med = 99.0;
            dradRMS3rms = 99.0;
          }
        }
        pTranCommon->peakDradRMS3 = dradRMS3rms;
        /* Count how many points are outside the flare */
        pTranCommon->peakExtra = 0;
        for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
          pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
          if ((pMalmquist4->Date < pMalmquist1->Date) ||
              (pMalmquist4->Date > pMalmquist3->Date)) {
            pTranCommon->peakExtra++;
          }
        }
        ave_magcal_magdep =  (max_magcal_magdep+min_magcal_magdep)/2.0;
        GetREF(pFileStarImage2->REFNumber,REF,0,0);
        variancesum = 0.0;
        variancecount = 0;
        for (temp_index = index_start; temp_index <= index_end; temp_index++) {
          pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
          if  (pMalmquist4->magcal_magdep <= ave_magcal_magdep) {
            pTranCommon->peakUpperCount++;
          }
          fitmag = (peakOffset/1000.) + ( pTranCommon->peakSlope *pMalmquist4->Date);
          fitdiff = pMalmquist4->magcal_magdep - fitmag;
          if (temp_index >= index_peak) {
            variancesum += fitdiff*fitdiff;
            variancecount++;
          }
        }
        if (variancecount > 0) {
          pTranCommon->peakRMS = sqrt(variancesum/(1.0*variancecount));
        } else {
          pTranCommon->peakRMS = 99;
        }
#ifdef PLOT_LIMITING_MAGNITUDES
        LoadGalaxyTable(pGalaxyCommon,pTranCommon->peakRA,pTranCommon->peakDec);
#if 0
        DumpLimitingTable(pGalaxyCommon);
#endif
        memset(plateLimitingYears,0,sizeof(plateLimitingYears));
        /* Convert the Julian Date to a year, and eliminate all points with known magnitudes */
        for (plateIndex = 0; plateIndex < pGalaxyCommon->plateCount; plateIndex++) {
          pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
          /* Skip plates that don't go as deep as the average brightness of the flare */
          if (pPlateLimitingRec->limiting_mag_local < ((pTranCommon->peakMaxMag+pTranCommon->peakMinMag)/2.0)) {
            continue;
          }
          pTranCommon->peakLimitingPoints++;
          dateyear = jd2ep(pPlateLimitingRec->geoJulianDate);
          intyear = dateyear-PIPELINE_MIN_DATE;
          if (intyear > (PIPELINE_MAX_DATE-PIPELINE_MIN_DATE)) {
            printf("ERROR: intyear %d is too large\n",intyear+PIPELINE_MIN_DATE);
            intyear = PIPELINE_MAX_DATE-PIPELINE_MIN_DATE;
          }
          if (plateLimitingYears[intyear] == 0) {
            pTranCommon->peakLimitingYears++;
          }
          plateLimitingYears[intyear]++;
        }
#endif /* PLOT_LIMITING_MAGNITUDES */
        pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_ACCEPTED);

        PlotTransientCandidates(pFileCommon,pTranCommon,pVectorStore,pMagnitudeTable,pGalaxyCommon,pCurStarEntry,refType,nearbyObjects,AFLAGSMASK1,QUALITYMASK,pMalmquist1,pMalmquist3);
      } else {
        pTranCommon->peakEventCount++;
      }
    }
    for (temp_index = index_start; temp_index <= index_end; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
      pMalmquist4->selectionFlags |= SELECT_TRANSIENT;
    }


  }
  pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_FINISHED);

}


/* Process Long Nova-like transients as requested by Josh memo of Tue, 20 Mar
 * 2018 15:47:56 and Sun, 1 Apr 2018 20:30:06. Prototype is
 * APASS_J075731.1+201735. Need:  >=10pts with >=1mag rise from pre-outburst.
 * This routine was originally ProcessStar in update_summary and brought over on
 * September 28, 2012 to process unmatched stars. The
 * pFileCommon->extraColumnFlag outputs columns not in the MySQL database
 * summary tables */
int
ProcessLightcurveParameters(
  FILE *consoleHandle,
  PFILECOMMON pFileCommon,
  MYSQL *pPhotConnection,
  PPHOTSTARIMAGE pMagnitudeTable,
  int numMagnitudes,
  PSTARENTRY pCurStarEntry,
  PPARAMETERVECTORSTORE pVectorStore,
  int minGoodStars,
  int versionId,
  FILE *summary_handle,
  int *summaryCount,
  int debugMode,
  int *staleVersionIdCount,
  int *damonSeriesId,
  int *excludeSeriesId,
  int excludeSeriesCount,
  PPHOTTARGET pTarget,
  int gsc_bin_index
) {
  TRANCOMMON transientCandidateCommon;
  PTRANCOMMON pTranCommon = &transientCandidateCommon;
  PPHOTSTARIMAGE pCurSumstarimage;
  PFILESTARIMAGE pCurFilestarimage;
#if 0
  double *vector1 = *pVector1;
  double *vector2 = *pVector2;
  double *vector3 = *pVector3;
  double *vector4 = *pVector4;
  double *vector5 = *pVector5;
  double *vector5B = *pVector5B;
  PMALMQUIST vector6 = *pVector6;
  PMALMQUIST vector6B = *pVector6B;
  double *vector7 = *pVector7;
  double *vector8 = *pVector8;
  PMALMQUIST tmpVector6;
  PMALMQUIST tmpVector6B;
  double *tmpVector;
  int vectorAlloc = *pVectorAlloc;
#endif
  double malmquistVector[MALMQUIST_STARS];
  int malmquistIndex;
  int imageIndex;
  int temp_ngood = 0; /* Count of AFLAGS = 0 points for further processing */
  int temp_selected = 0;
  int temp_index;
  int full_ngoodB = 0; /* Count of second quality points (Population 2) */
  int Sextractor_Blend = 0; /* Number of good stars (full_ngood) with the Sextractor BLEND or NEIGHBORS mask set */
  double full_med;
  double full_rms;
  double min_local;
  double min_local2;
  double max_local2;
  double max_local;
  int span;

  double pmagmin_local;
  double pmagmin_local2;
  double pmagmax_local2;
  double pmagmax_local;
  double mmagmin_local;
  double mmagmin_local2;
  double mmagmax_local2;
  double mmagmax_local;
  double range_local;
  double range_local2;


  double min_date;
  double min_date2;
  double max_date2;
  double max_date;
  double median_iso;
  double median_iso_rms;
  double min_iso;
  double min_iso2;
  double max_iso2;
  double max_iso;
  double rawrms;
  double rawerrmed;
  double rawerrrms;
  double error_bar_factor;
  int outStars = 0;
  double drad = 99.0;
  double max_drad = 0.0;
  double max_dradB = 0.0;
  int doClip = 1;
  double factor;
  double dra;
  double ddec;
  double dradmed;
  double dradmedB;
  double dradrms;
  double dradrmsB = 99.0;
  double Malmquist_factor = 99.0;
  double Malmquist_factorB = 99.0;
  double Damon_factor = 99.0;
  int excludeSeriesIndex;
  int excludeSeriesFlag;
  double precessDec;
  double precessRa;
  double plateepoch;
  int AFLAGSMASK1 = FILTER_AMASK_PLOT;
  int QUALITYMASK = 0;
  /* The following new variables were requested by Sumin Tang on 7/03/10 10:53 */
  int nblend = 0;                  /*  number of blended lightcurve points (Population 7) */
  int nNonDamonBlue = 0;           /*  number of good lightcurve points excluding dnb and dsb plates (Population 9)*/
  double medNonDamonBlue = 99;    /*  median of good lightcuve magnitudes excluding dnb and dsb plates*/
  double nonDamonBluerms = 99;    /*  rms of good lightcurve magnitudes excluding dnb and dsb plates*/
  int nDamonBlue = 0;              /*  number of good lightcurve points in dnb and dsb plates (Population 8) */
  double medDamonBlue = 99;       /*  median of good lightcurve magnitudes in dnb and dsb plates*/
  double damonBlueRms = 99;       /*  rms of good lightcurve magnitudes in dnb and dsb plates*/
  double magvslimitingcorr = 99;  /*  correlation coef. between lightcurve magnitude and limiting magnitude*/
  double magvsracorr = 99;        /*  correlation coef. between lightcurve magnitude and ra*/
  double magvsdeccorr = 99;       /*  correlation coef. between lightcurve magnitude and dec*/
  double rmsdradrms2 = 99;        /*  standard deviation of lightcurve dradRMS2*/
  double rarms = 99;              /*  standard deviation of lightcurve ra*/
  double decrms = 99;             /*  standard deviation of lightcurve dec*/
  int nburst = -1;                  /*  number of good lightcurve points 0.8 mag above the median value*/
  int nburst2 = -1;                 /*  number of good lightcurve points 0.5 mag above the median value*/
  int nburst3 = -1;                 /*  number of good lightcurve points 0.4 mag above the median value*/
  int nburst4 = -1;                 /*  number of good lightcurve points 3-sigma above the median value*/
  int ndip = -1;                    /*  number of good lightcurve points 0.8 mag below the median value*/
  int ndip2 = -1;                   /*  number of good lightcurve points 0.5 mag below the median value*/
  int ndip3 = -1;                   /*  number of good lightcurve points 0.4 mag below the median value*/
  int ndip4 = -1;                   /*  number of good lightcurve points 3-sigma below the median value*/
  int ndev3 = -1;                   /*  number of good lightcurve points 3-sigma from the median value*/
  int ndev2 = -1;                   /*  number of good lightcurve points 2-sigma from the median value*/
  int adjacentburstdip = -1;        /*  a measure of the number of adjacent nburst3/4 and ndip3/4 points in a lightcurve*/
  int adjacentburstdip2 = -1;       /*  a measure of number of >5 adjacent nburst3/4 and ndip3/4 points in a lightcurve*/
  int adjacentburstdip3 = -1;       /*  a measure of number of >7 adjacent points 2-sigma above or below the median value of a lightcurve*/
  double slope_all = 99;            /*  Slope (millimag/yr) of all good lightcurve points */
  double slope_all_err = 99;        /*  Uncertainty in slope_all */
  int nslope_all = 0;              /*  Number of points used to determine slope_all */
  double slope_60 = 99;            /*  Slope (millimag/yr) of good lightcurve points before 1960 */
  double slope_60_err = 99;        /*  Uncertainty in slope_60 */
  int nslope_60 = 0;              /*  Number of points used to determine slope_60 */

  double lightcurverms1 = 99;     /*  5-sigma clipped rms of good lightcurve points*/
  double lightcurverms2 = 99;     /*  5-sigma clipped rms of good lightcurve points after detrending (smooth(x,y,0.4,'lowess'))*/
  double lightcurverms3 = 99;     /*  5-sigma clipped rms of good lightcurve points after detrending (smooth(y,0.8,'lowess'))*/
  double lightcurverms4 = 99;     /*  5-sigma clipped rms of good lightcurve points after detrending (smooth(y,10,'sgolay'))*/
  double lightcurverms5 = 99;     /*  5-sigma clipped rms of good lightcurve points after detrending (smooth(y,15,'loess'))*/
  double rms_factor = 99;         /*  rms of the lightcurve divided by the median error; real variables have larger rmsf */
  double dmagcatalog = 99;        /*  the differences between catalog magnitude and median magnitude in the plates*/
  /* End of new variables were requested by Sumin Tang on 7/03/10 10:53 */

  PMALMQUIST pMalmquist;
  double tempramed;
  double temprarms;
  double tempdecmed;
  double tempdecrms;
  int errorFlag = 0;
  double mag_median = 0;
  int num_iburst3 = 0;
  int num_idip3 = 0;
  int num_iburst4 = 0;
  int num_idip4 = 0;
  int num_iburst5 = 0;
  int num_idip5 = 0;
  int iadj;
  int idip_index;
  int iburst_index;
  int index;
  int kadj2a;
  int kadj2b;
  int kadj2c;
  int kadj2d;
  int kadj2;
  double meanvalue;
  double stdvalue;
  int trendCount;
  double errreal_med;
  double errreal_rms;
  int keplerField = KEPLERFIELD_NO;
  const char *keplerFieldString[3] = {"no","yes","no"};
  int chn_n;
  double p_row;
  double p_coln;
#ifdef NEW_LIMITING_MAG
  double temp_median;
  double temp_rms;
#endif /* NEW_LIMITING_MAG */
  int latestVersionId;
  int quality;
  char nearbyObjects[MAX_NEARBY_OBJECTS_STRING + 2];
  char nearbyObjects2[MAX_NEARBY_OBJECTS_STRING + 4];
  PGALAXYCOMMON pGalaxyCommon = &pFileCommon->galaxyCommon;
  NEARESTCATALOGSTAR nullCatalogStar;
  PNEARESTCATALOGSTAR pNullCatalogStar = &nullCatalogStar;
  PNEARESTCATALOGSTAR pNearestCatalogStar;
  double lat;
  double lon;
  double min_FWHM_WORLD;
  int releaseField;
  double offset;
  double closestDistanceArcsec;
  static int printFlags = 0;
  double *pColorterm = NULL;
  int *pColorflag = NULL;

  memset(pTranCommon,0,sizeof(TRANCOMMON));
  pTranCommon->peakRA = 999.0;
  pTranCommon->peakDec = 99.0;
  pTranCommon->gsc_bin_index = gsc_bin_index;

  if (pFileCommon->enableTransientSearch != 0) {
#if 0
    /* These magnitudes are garbage */
    AFLAGSMASK1 |=  (1 << FILTER_AFLAG_TOO_BRIGHT);
#endif
    AFLAGSMASK1 |=  (1 << FILTER_AFLAG_BIN9) | (1 << FILTER_AFLAG_MULTIPLE_NONE); /* BIN9 and MULTIPLE_NONE added on August 1 */
    AFLAGSMASK1 &=  ~((1 << FILTER_AFLAG_WEDGE));  /* Defect removed on Jul 21, 2015 and added again on Aug 30, 2015  */
    AFLAGSMASK1 |=  (GSC_CLASS_NONSTAR << GSC_CLASS_BIT) ;  /* Search for Sumin's nova in Baade's window Added on Nov 11, 2015 */
#ifndef INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND
    AFLAGSMASK1 &=  ~((1 << FILTER_AFLAG_BACKGROUND));  /* Added Feb 12, 2017 */
#endif /* INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND */
#ifndef INCLUDE_TRANSIENT_DEFECTS
    AFLAGSMASK1 &=  ~((1 << FILTER_AFLAG_DEFECT));  /* Defect removed on Jul 21, 2015 and added again on Aug 30, 2015  */
#endif /* INCLUDE_TRANSIENT_DEFECTS */
#ifdef  INCLUDE_TRANSIENT_DRADBIN
    AFLAGSMASK1 |=  ((1 << FILTER_AFLAG_DRADBIN) | (GSC_CLASS_NONSTAR << GSC_CLASS_BIT)) ;  /* Search for Sumin's nova in Baade's window  */
#endif /* INCLUDE_TRANSIENT_DRADBIN */
#ifdef INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND
    AFLAGSMASK1 |=  (1 << FILTER_AFLAG_LOCAL_RMS);  /* LIMITING_MAG Added Feb 12, 2017 and removed Mar 15, 2017; LOCAL_RMS added Mar 13, 2017 */
#endif /* INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND */



    QUALITYMASK = QUALITY_GRATING | QUALITY_SPECTRA | QUALITY_COLOR;  /* On July 21, 2015, add QUALITY_SPECTRA   */
    /* On Apr 28, 2017, remove QUALITY_WEDGE */
#ifndef INCLUDE_TRANSIENT_MULTIPLE
    QUALITYMASK |= QUALITY_MULTIPLE;
#endif /* INCLUDE_TRANSIENT_MULTIPLE */
#ifdef INCLUDE_TRANSIENT_DRADBIN
    QUALITYMASK &= (~(QUALITY_WEDGE));
#endif /* INCLUDE_TRANSIENT_DRADBIN */

  }
  if (printFlags == 0) {
    fprintf(consoleHandle,"AFLAGSMASK1 is %d complement is %d\n",AFLAGSMASK1,((~AFLAGSMASK1) & 0x7FFFFFFF));
    fprintf(consoleHandle,"FILTER_AMASK_NONE_CANDIDATE2 is %d complement is %d\n",FILTER_AMASK_NONE_CANDIDATE2,((~FILTER_AMASK_NONE_CANDIDATE2) & 0x7FFFFFFF));
    if (strcmp(GetSeriesString(MF_TC_SERIESID,1),"mf") != 0) {
      fprintf(consoleHandle,"ERROR: Series id %d is %s instead of mf\n",MF_TC_SERIESID,GetSeriesString(MF_TC_SERIESID,1));
      exit(1);
    }

    if (pFileCommon->enableTransientSearch != 0) {
#ifdef INCLUDE_TRANSIENT_DEFECTS
      fprintf(consoleHandle,"ERROR: INCLUDE_TRANSIENT_DEFECTS is enabled\n");
#endif  /* INCLUDE_TRANSIENT_DEFECTS */
#ifdef INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND
      fprintf(consoleHandle,"ERROR: INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND is enabled\n");
#endif  /* INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND */

#ifdef INCLUDE_TRANSIENT_MULTIPLE
      fprintf(consoleHandle,"ERROR: INCLUDE_TRANSIENT_MULTIPLE is enabled\n");
#endif /* INCLUDE_TRANSIENT_MULTIPLE */

#ifdef INCLUDE_TRANSIENT_ONE_MULTIPLE
      fprintf(consoleHandle,"ERROR: INCLUDE_TRANSIENT_ONE_MULTIPLE is enabled\n");
#endif /* INCLUDE_TRANSIENT_ONE_MULTIPLE */

#ifdef INCLUDE_TRANSIENT_DRADBIN
      fprintf(consoleHandle,"ERROR: INCLUDE_TRANSIENT_DRADBIN is enabled\n");
#endif  /* INCLUDE_TRANSIENT_DRADBIN */
      fprintf(consoleHandle,"candidate image suffix is %s.png\n",pFileCommon->timestr);
    }
    printFlags = 1;
  }
#if 0
  fprintf(consoleHandle,"FILTER_AMASK_PLOT is %d complement is %d\n",FILTER_AMASK_PLOT,((~FILTER_AMASK_PLOT) & 0x7FFFFFFF));
  fprintf(consoleHandle,"((1 << FILTER_AFLAG_WEDGE) | (1 << FILTER_AFLAG_DEFECT)) is %d complement is %d\n",((1 << FILTER_AFLAG_WEDGE) | (1 << FILTER_AFLAG_DEFECT)),((~((1 << FILTER_AFLAG_WEDGE) | (1 << FILTER_AFLAG_DEFECT))) & 0x7FFFFFFF));


  fprintf(consoleHandle,"FILTER_AMASK_BLEND is %d\n",FILTER_AMASK_BLEND);
  fprintf(consoleHandle,"FILTER_BLEND_MASK_SEXTRACTOR is %d\n",FILTER_BLEND_MASK_SEXTRACTOR);
  fprintf(consoleHandle,"AFLAGSMASK3 is %d\n",AFLAGSMASK3);
  fprintf(consoleHandle,"Processing %s\n",pCurStarEntry->REF);
#endif

#if 0
  if (strcmp(pCurStarEntry->REF,"K5292705") == 0) {
    fprintf(consoleHandle,"Processing %s\n",pCurStarEntry->REF);
  }

#endif

  if (RA2fpix(pCurStarEntry->ra,pCurStarEntry->dec,&chn_n,&p_row,&p_coln) != 0) {
    keplerField = KEPLERFIELD_YES;
  }
  InitNearestCatalogStar(pNullCatalogStar);

  for (imageIndex = 0; imageIndex < numMagnitudes; imageIndex++) {
    pCurSumstarimage = &pMagnitudeTable[imageIndex];
    pCurFilestarimage = pCurSumstarimage->pFileStarImage;

    if (pCurFilestarimage->gsc_bin_index > pCurStarEntry->gsc_bin_index) {
      break;
    }

    if (strcmp(pCurSumstarimage->REF,pCurStarEntry->REF) == 0) {
      latestVersionId = GetLatestVersionId(pFileCommon,pCurFilestarimage->seriesId,pCurFilestarimage->plateNumber,0,pCurFilestarimage->REFNumber,1);

      if ((latestVersionId < 0) || (pCurFilestarimage->versionId < latestVersionId )) {
        (*staleVersionIdCount)++;
        continue;
      }
      if (pTranCommon->npoints >= pVectorStore->vectorAlloc) {
        pVectorStore->vectorAlloc += PARAMETERVECTOR_ALLOC_INCREMENT;
        ReallocParameterVectorStore(pVectorStore);
      }
      pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][pTranCommon->npoints].imageIndex = imageIndex;
      pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][pTranCommon->npoints].selectionFlags = 0;
      pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][pTranCommon->npoints].Date = pCurFilestarimage->Date;
      pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][pTranCommon->npoints].magcal_magdep = pCurFilestarimage->magcal_magdep;
      pTranCommon->npoints++; /* Population 0 (Population numbers are referenced in the photometry paper */

      if ((pCurFilestarimage->AFLAGS & (~AFLAGSMASK1)) != 0) {
        continue;  /* Reject Population 1 */
      }
      if (pFileCommon->enableTransientSearch) {
        /* Reject any image from mf plates at Bloemfontain (mf10104+) where the image is within 1.9 degrees of the top or bottom. See memo of Tue 9/29/15 11:02 AM */
        if ((pCurFilestarimage->seriesId == MF_TC_SERIESID) && (pCurFilestarimage->plateNumber >=  MF_TC_EXCLUSION_PIXELS)) {
          if ((pCurFilestarimage->Y_IMAGE < MF_TC_EXCLUSION_PIXELS) || ((MF_TC_PLATE_HEIGHT-pCurFilestarimage->Y_IMAGE) < MF_TC_EXCLUSION_PIXELS)) {
            pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_EXCLUSIONZONE);
#if 0
            fprintf(consoleHandle,".");
#endif
            continue;
          }
        }
      }
      if (full_ngoodB >= pVectorStore->vectorAlloc) {
        pVectorStore->vectorAlloc += PARAMETERVECTOR_ALLOC_INCREMENT;
        ReallocParameterVectorStore(pVectorStore);
      }
      plateepoch = jd2ep(pCurFilestarimage->Date) - GSC_EQUINOX;
      if ((plateepoch+GSC_EQUINOX > PIPELINE_MIN_DATE) &&
          (plateepoch+GSC_EQUINOX < PIPELINE_MAX_DATE)) {

        precessDec =  pCurStarEntry->dec + (pCurStarEntry->DecPM * plateepoch)/(3600.0*1000.0);
        factor = cos(DEGREES_TO_RAD*((pCurFilestarimage->dec + precessDec)/2.0));
        if (factor != 0) {
          precessRa = pCurStarEntry->ra + ((pCurStarEntry->RaPM *plateepoch)/(3600.0*1000.0))/factor;
        } else {
          precessRa = pCurStarEntry->ra;
        }
      } else {
        fprintf(
          consoleHandle,
          "ERROR: image date %f not between %f and %f\n",
          plateepoch + GSC_EQUINOX,
          (float) PIPELINE_MIN_DATE,
          (float) PIPELINE_MAX_DATE
        );
        exit(1);
      }

      dra = 3600*factor*(pCurFilestarimage->ra - precessRa);
      ddec = 3600*(pCurFilestarimage->dec - precessDec);
      drad = sqrt((dra*dra) + (ddec*ddec));

      pVectorStore->vector[PARAMETERVECTOR_BDRAD][full_ngoodB] = drad;
      if (drad > max_dradB) {
        max_dradB = drad;
      }
#ifdef DUMP_VALUES
      if (pCurStarEntry->REFNumber == LOS_REF_NUMBER) {
        fprintf(consoleHandle,"IndexM %d PARAMETERVECTOR_BMALMQUIST %f\n",full_ngoodB,pCurFilestarimage->magcal_magdep);
      }
#endif /* DUMP_VALUES */
      pVectorStore->mvector[PARAMETERVECTOR_BMALMQUIST][full_ngoodB].limiting_mag_local = pCurFilestarimage->limiting_mag_local;
      pVectorStore->mvector[PARAMETERVECTOR_BMALMQUIST][full_ngoodB].magcal_magdep = pCurFilestarimage->magcal_magdep;
      full_ngoodB++; /* Population 2 */

#if 0
      if ((strcmp(pCurFilestarimage->series,"mb") == 0) &&
          (pCurFilestarimage->plateNumber == 2450)) {
        fprintf(consoleHandle,"At plate %s%05d\n",pCurFilestarimage->series,pCurFilestarimage->plateNumber);
      }

#endif
      if (pFileCommon->enableTransientSearch != 0) {
        /* Reject multiple exposure plates (Change of Mar 24, 2015) */
        GetPlateQualityMask(
          pFileCommon,
          pCurFilestarimage->seriesId,
          pCurFilestarimage->plateNumber,
          &quality,
          &pColorterm,
          &pColorflag,
          pCurFilestarimage->REFNumber
        );
        pCurSumstarimage->quality = quality;
        pCurSumstarimage->colorterm = pColorterm[pCurFilestarimage->spatial_bin - 1];
        pCurSumstarimage->colorflag = pColorflag[pCurFilestarimage->spatial_bin - 1];
#if 0
        if (quality != 0) {
          fprintf(consoleHandle,"quality 0x%x mask 0x%x for seriesId %d plateNumber %d NUMBER %d\n",quality,(QUALITY_MULTIPLE | QUALITY_GRATING |  QUALITY_WEDGE),pCurFilestarimage->seriesId,pCurFilestarimage->plateNumber,pCurFilestarimage->NUMBER);
        }
#endif
#if 0 /* Test of Apr 28, 2017 */
        if ((quality & QUALITY_WEDGE) != 0) {
          printf("line %d quality is %d 0x%x for plate %s gsc_bin_index %d NUMBER %d WEDGE %d\n",__LINE__,quality,quality,pCurSumstarimage->Plate,gsc_bin_index,pCurSumstarimage->pFileStarImage->NUMBER,(pCurSumstarimage->pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_WEDGE)) != 0 );

        }
#endif
        if ((quality & QUALITYMASK) != 0) { /* On July 21, 2015, add QUALITY_SPECTRA   */
#if 0 /* Test of Apr 28, 2017 */
          if ((quality & QUALITY_WEDGE) != 0) {
            printf("line %d quality is %d 0x%x for plate %s gsc_bin_index %d NUMBER %d WEDGE %d\n",__LINE__,quality,quality,pCurSumstarimage->Plate,gsc_bin_index,pCurSumstarimage->pFileStarImage->NUMBER,(pCurSumstarimage->pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_WEDGE)) != 0 );
          }
#endif
          pFileCommon->counterBlock.qualityRejectCount++;
#if 0
          fprintf(consoleHandle,"qualityRejectCount %d seriesId %d plateNumber %d NUMBER %d\n",pFileCommon->counterBlock.qualityRejectCount,pCurFilestarimage->seriesId,pCurFilestarimage->plateNumber,pCurFilestarimage->NUMBER);
#endif
          continue;
        }
#if 0 /* Test of Apr 28, 2017 */
        if ((quality & QUALITY_WEDGE) != 0) {
          printf("line %d quality is %d 0x%x for plate %s gsc_bin_index %d NUMBER %d WEDGE %d\n",__LINE__,quality,quality,pCurSumstarimage->Plate,gsc_bin_index,pCurSumstarimage->pFileStarImage->NUMBER,(pCurSumstarimage->pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_WEDGE)) != 0 );

        }
#endif

      }

#ifndef PERMISSIVE_ID_TABLE
      /* Reject second quality colorterm values (Change of Jun 13, 2011) */
      qualityMask = GetPlateQualityMask(pFileCommon,pCurFilestarimage->seriesId,pCurFilestarimage->plateNumber,&quality,&pColorterm,&pColorflag,pCurFilestarimage->REFNumber);
      if ((qualityMask & (1 << pCurFilestarimage->spatial_bin)) != 0) {
        pFileCommon->counterBlock.qualityRejectCount++;
        continue;
      }

      tmpAFLAGS = pCurFilestarimage->AFLAGS & (~AFLAGSMASK3);
      /* Select this star if it is blended */
      if (((tmpAFLAGS & FILTER_AMASK_BLEND) != 0) &&
          ((tmpAFLAGS & (~(FILTER_AMASK_BLEND))) == 0)) {
        /* Save this one because it may be tallied against the blend count */
        pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_selected].imageIndex = imageIndex;
        pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_selected].selectionFlags = SELECT_TEMP_BLEND;
        pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_selected].Date = pCurFilestarimage->Date;
        temp_selected++; /* Population 3 */

        continue;
      }
      if ((pCurFilestarimage->AFLAGS & (~AFLAGSMASK2)) != 0) {
        continue;
      }
#endif /* PERMISSIVE_ID_TABLE */


#ifdef NEW_LIMITING_MAG
#ifndef PERMISSIVE_ID_TABLE
      /* Reject stars that do not have magnitude-dependent correction */
      if ((pCurFilestarimage->BFLAGS &  (1 << FILTER_BFLAG_MAGDEP_MAGCOR)) == 0) {
        continue;
      }
#endif /* PERMISSIVE_ID_TABLE */
      /* (((pCurFilestarimage->limiting_mag_local-pCurStarEntry->Stdmag) <= 0.5) ||
         ((pCurFilestarimage->limiting_mag_local-pCurFilestarimage->magcal_magdep) <= 0.75) ||
         (((2*pCurFilestarimage->limiting_mag_local)-pCurStarEntry->Stdmag-pCurFilestarimage->magcal_magdep) <= 1.5)) */

#ifndef INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND
      /* Reject this star if it does not satisfy Sumin's criteria of Thu 12/15/11 1:52 PM */
      if ((pCurFilestarimage->limiting_mag_local-pCurFilestarimage->magcal_magdep) <= 0.75)  {
        pFileCommon->counterBlock.limitingmag1++;
#if 0
        fprintf(consoleHandle,"Rejecting %s\n",pCurStarEntry->REF);
#endif
        continue;
      }
#endif /* INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND */
#endif /* NEW_LIMITING_MAG */


      /* The quality mask check was here before Jun 13, 2011 */


      excludeSeriesFlag = 0;
      for (excludeSeriesIndex = 0; excludeSeriesIndex < excludeSeriesCount; excludeSeriesIndex++) {

        if (pCurFilestarimage->seriesId == excludeSeriesId[excludeSeriesIndex]) {
          pFileCommon->counterBlock.excludeSeriesCount++;
          excludeSeriesFlag = 1;
          break;
        }
      }
      if (excludeSeriesFlag) {
        continue;
      }
#ifdef PERMISSIVE_ID_TABLE
      pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_selected].selectionFlags = SELECT_TEMP_GOOD;  /* Population 5 */
#else /* PERMISSIVE_ID_TABLE */
      if ((pCurFilestarimage->BFLAGS & FILTER_BLEND_MASK_SEXTRACTOR) != 0) {
        /* Reject Sextractor blends */
        pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_selected].selectionFlags = SELECT_TEMP_BLEND; /* Population 4 */
      } else {
        /* This point is tentatively good.  Save */
        pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_selected].selectionFlags = SELECT_TEMP_GOOD;  /* Population 5 */
      }
#endif /* PERMISSIVE_ID_TABLE */
      pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_selected].imageIndex = imageIndex;
      pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_selected].Date = pCurFilestarimage->Date;
      temp_selected++;
    }
  }
  /* Exit if nothing more to do */
  if (temp_selected < minGoodStars) {
    return(outStars);
  }
  /* Aug 3, 2015 - use George Miller's criterion to look for flare tops only */
  if (pFileCommon->enableTransientSearch != 0) {
    if (pTranCommon->npoints < MIN_TC_POINTS) {
      pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_INSUFFICIENT_POINTS2);
      return(outStars);
    }
  }

  /* Sort the selection vector for increasing date */
  qsort((void*)pVectorStore->mvector[PARAMETERVECTOR_SELECTION],temp_selected,sizeof(MALMQUIST),LightcurveDateCompare);

  /* Now go through our selected stars and find the median, rms dec and ra */
  for (temp_index = 0; temp_index < temp_selected; temp_index++) {
    pMalmquist = &pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_index];
    if ((pMalmquist->selectionFlags & SELECT_TEMP_GOOD) != 0) {
      pCurSumstarimage = &pMagnitudeTable[pMalmquist->imageIndex];
      pCurFilestarimage = pCurSumstarimage->pFileStarImage;
      pVectorStore->vector[PARAMETERVECTOR_RA][temp_ngood] = pCurFilestarimage->ra;
      pVectorStore->vector[PARAMETERVECTOR_DEC][temp_ngood] = pCurFilestarimage->dec;
#ifdef NEW_LIMITING_MAG
      pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP5][temp_ngood] = pCurFilestarimage->magcal_magdep;
#endif /* NEW_LIMITING_MAG */
      temp_ngood++;
    }


  }
  /* Exit if nothing more to do */
  if (temp_ngood < minGoodStars) {
    return(outStars);
  }
  /* Find the median and rms right ascension and declination */
  if (CalcMedianAndRMS(temp_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_RA],&tempramed,&temprarms,0,3.0,0) == 0) {
    return(outStars);
  }
  if (CalcMedianAndRMS(temp_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_DEC],&tempdecmed,&tempdecrms,0,3.0,0) == 0) {
    return(outStars);
  }
#ifdef NEW_LIMITING_MAG
  if (CalcMedianAndRMS(temp_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP5],&temp_median,&temp_rms,0,3.0,0) == 0) {
    return(outStars);
  }
#endif /* NEW_LIMITING_MAG */

  /* Avoid divide by zero.  Use 0.1 arcsec minimum */
  if (temprarms == 0.0) {
    temprarms = 0.000028;
  }
  if (tempdecrms == 0.0) {
    tempdecrms = 0.000028;
  }

  /* Now make the final selections based on ra and dec */
  min_FWHM_WORLD = 3600.; /* Start with a very large value */
  for (temp_index = 0; temp_index < temp_selected; temp_index++) {
    pMalmquist = &pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_index];
    pCurSumstarimage = &pMagnitudeTable[pMalmquist->imageIndex];
    pCurFilestarimage = pCurSumstarimage->pFileStarImage;
#if 0
    if ((strcmp(pCurFilestarimage->series,"dnb") == 0) &&
        (pCurFilestarimage->plateNumber == 3741)) {
      fprintf(consoleHandle,"At plate %s%05d\n",pCurFilestarimage->series,pCurFilestarimage->plateNumber);
      fprintf(consoleHandle,"ra %f tempramed %f temprarms %f, error %f\n",
              pCurFilestarimage->ra,tempramed,temprarms,((fabs(pCurFilestarimage->ra - tempramed)/temprarms) - 4.0));
      fprintf(consoleHandle,"dec %f tempdecmed %f tempdecrms %f, error %f\n",
              pCurFilestarimage->dec,tempdecmed,tempdecrms,((fabs(pCurFilestarimage->dec - tempdecmed)/tempdecrms) - 4.0));
    }

#endif


    if (((fabs(pCurFilestarimage->ra - tempramed)/temprarms) > 4.0)  ||
        ((fabs(pCurFilestarimage->dec - tempdecmed)/tempdecrms) > 4.0)) {
      continue;
    }
#ifdef NEW_LIMITING_MAG
    /* Reject this star if it does not satisfy Sumin's criteria of Thu 12/15/11 1:52 PM */
    if (((pCurFilestarimage->limiting_mag_local-temp_median) <= 0.5) ||
        (((2*pCurFilestarimage->limiting_mag_local)-temp_median-pCurFilestarimage->magcal_magdep) <= 1.5)) {
#ifndef INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND
      pFileCommon->counterBlock.limitingmag2++;
#if 0
      fprintf(consoleHandle,"Rejecting %s\n",pCurStarEntry->REF);
#endif
      continue;
#endif /* INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND */
    }

#endif /* NEW_LIMITING_MAG */

    if ((pMalmquist->selectionFlags & SELECT_TEMP_BLEND) != 0) {
      pMalmquist->selectionFlags |= SELECT_BLEND;

      nblend++; /* Population 7 */
    }

    if ((pMalmquist->selectionFlags & SELECT_TEMP_GOOD) != 0) {
      double errreal;
      pMalmquist->selectionFlags |= SELECT_GOOD;  /* Population 6 */

      pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][pTranCommon->full_ngood] = pCurFilestarimage->magcal_magdep;
      pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP2][pTranCommon->full_ngood] = pCurFilestarimage->magcal_magdep;
      pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][pTranCommon->full_ngood] = pCurFilestarimage->magcal_magdep;
      pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP4][pTranCommon->full_ngood] = pCurFilestarimage->magcal_magdep;
      pVectorStore->vector[PARAMETERVECTOR_XVALS][pTranCommon->full_ngood] = pTranCommon->full_ngood+1;

      pVectorStore->vector[PARAMETERVECTOR_JULIAN_DAY][pTranCommon->full_ngood] = pCurFilestarimage->Date;
      pVectorStore->vector[PARAMETERVECTOR_DATEREAL][pTranCommon->full_ngood] = pCurFilestarimage->Date - 2400000.;
      pVectorStore->vector[PARAMETERVECTOR_MAGCAL_ISO][pTranCommon->full_ngood] = pCurFilestarimage->magcal_iso;
      pVectorStore->vector[PARAMETERVECTOR_MAGCAL_LOCAL_RMS][pTranCommon->full_ngood] = pCurFilestarimage->magcal_local_rms;
      errreal = sqrt(sqr(pCurFilestarimage->magcal_local_error) + sqr(pCurFilestarimage->magcal_iso_rms))/2.0;
      pVectorStore->vector[PARAMETERVECTOR_ERRREAL][pTranCommon->full_ngood] = errreal;
      pVectorStore->vector[PARAMETERVECTOR_PMAG][pTranCommon->full_ngood] = pCurFilestarimage->magcal_magdep + errreal;
      pVectorStore->vector[PARAMETERVECTOR_MMAG][pTranCommon->full_ngood] = pCurFilestarimage->magcal_magdep - errreal;



      dra = 3600*factor*(pCurFilestarimage->ra - precessRa);
      ddec = 3600*(pCurFilestarimage->dec - precessDec);
      drad = sqrt((dra*dra) + (ddec*ddec));
#ifdef DUMP_VALUES
      if (pCurStarEntry->REFNumber == LOS_REF_NUMBER) {
        fprintf(consoleHandle,"Index %d PARAMETERVECTOR_DRAD %f\n",pTranCommon->full_ngood,drad);
      }
#endif /* DUMP_VALUES */

      pVectorStore->vector[PARAMETERVECTOR_DRAD][pTranCommon->full_ngood] = drad;
      if (drad > max_drad) {
        max_drad = drad;
      }
      pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood].limiting_mag_local = pCurFilestarimage->limiting_mag_local;
      pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood].magcal_magdep = pCurFilestarimage->magcal_magdep;
      pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood].imageIndex = pMalmquist->imageIndex;
      pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood].selectionFlags = pMalmquist->selectionFlags;
      pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood].Date = pCurFilestarimage->Date;
      pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood].drad = drad;
      pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood].dradRMS2 = pCurFilestarimage->dradRMS2;
      pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood].AFLAGS = pCurFilestarimage->AFLAGS;


      /* check systematic differences between dnb plates and non-dn plates */

      if ((pCurFilestarimage->seriesId == damonSeriesId[0]) ||
          (pCurFilestarimage->seriesId == damonSeriesId[1])) {
        pVectorStore->vector[PARAMETERVECTOR_DAMON][nDamonBlue] = pCurFilestarimage->magcal_magdep;
        nDamonBlue++; /* Population 9 */
      } else if ((pCurFilestarimage->seriesId != damonSeriesId[2]) &&
                 (pCurFilestarimage->seriesId != damonSeriesId[3]) &&
                 (pCurFilestarimage->seriesId != damonSeriesId[4]) &&
                 (pCurFilestarimage->seriesId != damonSeriesId[5])) {
        pVectorStore->vector[PARAMETERVECTOR_NONDAMON][nNonDamonBlue] = pCurFilestarimage->magcal_magdep;
        nNonDamonBlue++; /* Population 8 */

      }
#if 0

      fprintf(consoleHandle,"%s%05d\n",pCurFilestarimage->series,pCurFilestarimage->plateNumber);
#endif
      /* min_FWHM_WORLD is used to check against nearby catalog objects */
      if (min_FWHM_WORLD > pCurFilestarimage->FWHM_WORLD) {
        min_FWHM_WORLD = pCurFilestarimage->FWHM_WORLD;
      }

      pTranCommon->full_ngood++; /* Population 6 */
      if ((pCurFilestarimage->BFLAGS & FILTER_BLEND_MASK_SEXTRACTOR) != 0) {
        Sextractor_Blend++;
      }
    } /* Good stars */
  } /* selection loop */
  if (min_FWHM_WORLD < MIN_FWHM_WORLD) {
    min_FWHM_WORLD = MIN_FWHM_WORLD;
  }


  if (pTranCommon->full_ngood < minGoodStars) {
    return(outStars);
  }
  /* Here we have enough stars to write an output record */


  if (CalcMedianAndRMS(pTranCommon->full_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP2],&full_med,&full_rms,0,3.0,0) == 0) {
    errorFlag = 1;
  }

  if (errorFlag == 0) {
    mag_median = full_med;
    ComputeRange(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP2],&min_local,&min_local2,&max_local2,&max_local);
    ComputeRange(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_JULIAN_DAY],&min_date,&min_date2,&max_date2,&max_date);
    ComputeRange(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_PMAG],&pmagmin_local,&pmagmin_local2,&pmagmax_local2,&pmagmax_local);
    ComputeRange(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_MMAG],&mmagmin_local,&mmagmin_local2,&mmagmax_local2,&mmagmax_local);
    range_local = mmagmax_local - pmagmin_local;
    range_local2 = mmagmax_local2 - pmagmin_local2;

    CalcMedianAndRMS(pTranCommon->full_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_ISO],&median_iso,&median_iso_rms,0,3.0,0);
    ComputeRange(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_ISO],&min_iso,&min_iso2,&max_iso2,&max_iso);

    /* Find the zero-based clipped rms of our error bars */
    rawrms = full_rms;
    if (CalcMedianAndRMS(pTranCommon->full_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_LOCAL_RMS],&rawerrmed,&rawerrrms,0,3.0,1) == 0) {
      rawerrmed = 0.0;
      rawerrrms = 99.0;
    }
    if ((rawrms == 0.0) ||
        (rawrms == 99.0) ||
        (rawerrrms == 0.0) ||
        (rawerrrms == 99.0)) {
      error_bar_factor = 1.0;
    } else {
      error_bar_factor = rawrms/rawerrrms;
    }
    if (error_bar_factor > 1.0) {
      error_bar_factor = 1.0;
    }

#ifdef DUMP_VALUES
    if (pCurStarEntry->REFNumber == LOS_REF_NUMBER) {
      DumpVector(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_DRAD]);
    }
#endif /* DUMP_VALUES */


    if (CalcMedianAndRMS(pTranCommon->full_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_DRAD],&dradmed,&dradrms,0,3.0,1) == 0) {
      dradmed = 99.0;
      dradrms = 99.0;
    }
    if (CalcMedianAndRMS(full_ngoodB,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_BDRAD],&dradmedB,&dradrmsB,0,3.0,1) == 0) {
      dradmedB = 99.0;
      dradrmsB = 99.0;
    }

    if (pTranCommon->full_ngood > (2*MALMQUIST_STARS)) {
      double lowermed;
      double lowerrms;
      double uppermed;
      double upperrms;
      /* We can create a Malmquist parameter.
         1. Sort the table in limiting magnitude
         2. Find the median of MALMQUIST_STARS with the lowest limiting magnitude.
         3. Find the median of MALMQUIST_STARS with the highest limiting magnitude.
         4. Take the difference between steps 3 and 4 above.
      */
      qsort((void*)pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST],pTranCommon->full_ngood,sizeof(MALMQUIST),LightcurveImageCompare);
      for (malmquistIndex = 0; malmquistIndex < MALMQUIST_STARS; malmquistIndex++) {
        malmquistVector[malmquistIndex] = pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][malmquistIndex].magcal_magdep;
      }
      if (CalcMedianAndRMS(MALMQUIST_STARS,0,malmquistVector,&lowermed,&lowerrms,1,3.0,0) != 0) {
        for (malmquistIndex = 0; malmquistIndex < MALMQUIST_STARS; malmquistIndex++) {
          malmquistVector[malmquistIndex] = pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][pTranCommon->full_ngood+malmquistIndex-MALMQUIST_STARS].magcal_magdep;
        }


        if (CalcMedianAndRMS(MALMQUIST_STARS,0,malmquistVector,&uppermed,&upperrms,1,3.0,0) != 0) {

          Malmquist_factor = uppermed - lowermed;
        }
      }
    }
    /* Repeat the above calculation for lower-quality points */
    if (full_ngoodB > (2*MALMQUIST_STARS)) {
      double lowermed;
      double lowerrms;
      double uppermed;
      double upperrms;
      /* We can create a Malmquist parameter.
         1. Sort the table in limiting magnitude
         2. Find the median of MALMQUIST_STARS with the lowest limiting magnitude.
         3. Find the median of MALMQUIST_STARS with the highest limiting magnitude.
         4. Take the difference between steps 3 and 4 above.
      */
#ifdef DUMP_VALUES
      if (pCurStarEntry->REFNumber == LOS_REF_NUMBER) {
        DumpMVector(full_ngoodB,pVectorStore->mvector[PARAMETERVECTOR_BMALMQUIST]);
      }

#endif /* DUMP_VALUES */
      qsort((void*)pVectorStore->mvector[PARAMETERVECTOR_BMALMQUIST],full_ngoodB,sizeof(MALMQUIST),LightcurveImageCompare);
      for (malmquistIndex = 0; malmquistIndex < MALMQUIST_STARS; malmquistIndex++) {
        malmquistVector[malmquistIndex] = pVectorStore->mvector[PARAMETERVECTOR_BMALMQUIST][malmquistIndex].magcal_magdep;
      }
      if (CalcMedianAndRMS(MALMQUIST_STARS,0,malmquistVector,&lowermed,&lowerrms,1,3.0,0) != 0) {
        for (malmquistIndex = 0; malmquistIndex < MALMQUIST_STARS; malmquistIndex++) {
          malmquistVector[malmquistIndex] = pVectorStore->mvector[PARAMETERVECTOR_BMALMQUIST][full_ngoodB+malmquistIndex-MALMQUIST_STARS].magcal_magdep;
        }


        if (CalcMedianAndRMS(MALMQUIST_STARS,0,malmquistVector,&uppermed,&upperrms,1,3.0,0) != 0) {

          Malmquist_factorB = uppermed - lowermed;
        }
      }
    }
    if (nDamonBlue > 2) {
      if (CalcMedianAndRMS(nDamonBlue,0,pVectorStore->vector[PARAMETERVECTOR_DAMON],&medDamonBlue,&damonBlueRms,0,3.0,0) > 0) {
      } else {
        medDamonBlue = 99;
        damonBlueRms = 99;
      }
    } else {
      medDamonBlue = 99;
      damonBlueRms = 99;
    }
    if (nNonDamonBlue > 2) {
      if (CalcMedianAndRMS(nNonDamonBlue,0,pVectorStore->vector[PARAMETERVECTOR_NONDAMON],&medNonDamonBlue,&nonDamonBluerms,0,3.0,0) > 0) {
      } else {
        medNonDamonBlue = 99;
        nonDamonBluerms = 99;
      }
    } else {
      medNonDamonBlue = 99;
      nonDamonBluerms = 99;
    }




    if ((nNonDamonBlue >= DAMON_STARS) &&
        (nDamonBlue >= DAMON_STARS)) {
      double nonDamonmed;
      double nonDamonrms;
      double Damonmed;
      double Damonrms;
      /* Bugfix of Jun 13, 2011 */
      int nNonDamonBlueTmp;
      int nDamonBlueTmp;
      nNonDamonBlueTmp = CalcMedianAndRMS(nNonDamonBlue,0,pVectorStore->vector[PARAMETERVECTOR_NONDAMON],&nonDamonmed,&nonDamonrms,1,3.0,0);
      nDamonBlueTmp    = CalcMedianAndRMS(nDamonBlue   ,0,pVectorStore->vector[PARAMETERVECTOR_DAMON],&Damonmed,&Damonrms,1,3.0,0);
      if ((nNonDamonBlueTmp >= DAMON_STARS) &&
          (nDamonBlueTmp >= DAMON_STARS)) {
        Damon_factor = Damonmed - nonDamonmed;
      }
    }

    if ((pTranCommon->clip_ngood = CalcMedianAndRMS(pTranCommon->full_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP2],&pTranCommon->clip_med,&pTranCommon->clip_rms,doClip,3.0,0)) == 0) {
      errorFlag = 1;
    }
    if ((errorFlag == 0) && (pFileCommon->provideWebSummary == 0) && (pTranCommon->full_ngood >= MIN_GOODSTARS)) {
      /* Compute the remaining parameters */
      temp_ngood = 0;
      nburst = 0;
      nburst2 = 0;
      nburst3 = 0;
      nburst4 = 0;
      ndip = 0;
      ndip2 = 0;
      ndip3 = 0;
      ndip4 = 0;
      ndev2 = 0;
      ndev3 = 0;
      for (temp_index = 0; temp_index < temp_selected; temp_index++) {
        pMalmquist = &pVectorStore->mvector[PARAMETERVECTOR_SELECTION][temp_index];
        pCurSumstarimage = &pMagnitudeTable[pMalmquist->imageIndex];
        pCurFilestarimage = pCurSumstarimage->pFileStarImage;
        if ((pMalmquist->selectionFlags & SELECT_GOOD) == 0) {
          continue;
        }
        /* need new ra and dec medians */
        pVectorStore->vector[PARAMETERVECTOR_RA][temp_ngood] = pCurFilestarimage->ra;
        pVectorStore->vector[PARAMETERVECTOR_DEC][temp_ngood] = pCurFilestarimage->dec;
        pVectorStore->vector[PARAMETERVECTOR_LIMITING_MAG][temp_ngood] = pCurFilestarimage->limiting_mag_local;
        pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][temp_ngood] = pCurFilestarimage->dradRMS2;

        if ((mag_median - pVectorStore->vector[PARAMETERVECTOR_PMAG][temp_ngood]) > 0.8) {
          nburst++;
        }
        if ((mag_median - pVectorStore->vector[PARAMETERVECTOR_PMAG][temp_ngood]) > 0.5) {
          nburst2++;
        }
        if ((mag_median - pVectorStore->vector[PARAMETERVECTOR_PMAG][temp_ngood]) > 0.4) {
          nburst3++;
          pVectorStore->ivector[PARAMETERVECTOR_IBURST3][num_iburst3] = temp_ngood;
          num_iburst3++;
        }
        /* 3-sigma burst */
        if ((-pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][temp_ngood]+mag_median) > (3 * pVectorStore->vector[PARAMETERVECTOR_ERRREAL][temp_ngood])) {
          nburst4++;
          pVectorStore->ivector[PARAMETERVECTOR_IBURST4][num_iburst4] = temp_ngood;
          num_iburst4++;
        }

        /* 2-sigma burst */
        if ((-pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][temp_ngood]+mag_median) > (2 * pVectorStore->vector[PARAMETERVECTOR_ERRREAL][temp_ngood])) {
          pVectorStore->ivector[PARAMETERVECTOR_IBURST5][num_iburst5] = temp_ngood;
          num_iburst5++;
        }


        if ((pVectorStore->vector[PARAMETERVECTOR_MMAG][temp_ngood] - mag_median) > 0.8) {
          ndip++;
        }
        if ((pVectorStore->vector[PARAMETERVECTOR_MMAG][temp_ngood] - mag_median) > 0.5) {
          ndip2++;
        }
        if ((pVectorStore->vector[PARAMETERVECTOR_MMAG][temp_ngood] - mag_median) > 0.4) {
          ndip3++;
          pVectorStore->ivector[PARAMETERVECTOR_IDIP3][num_idip3] = temp_ngood;
          num_idip3++;
        }

        /*  3-sigma dip */
        if ((pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][temp_ngood]-mag_median) > (3 * pVectorStore->vector[PARAMETERVECTOR_ERRREAL][temp_ngood])) {
          ndip4++;
          pVectorStore->ivector[PARAMETERVECTOR_IDIP4][num_idip4] = temp_ngood;
          num_idip4++;
        }

        /* 2-sigma dip */
        if ((pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][temp_ngood]-mag_median) > (2 * pVectorStore->vector[PARAMETERVECTOR_ERRREAL][temp_ngood])) {
          pVectorStore->ivector[PARAMETERVECTOR_IDIP5][num_idip5] = temp_ngood;
          num_idip5++;
        }

        if (fabs(pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][temp_ngood]-mag_median) > (3 * pVectorStore->vector[PARAMETERVECTOR_ERRREAL][temp_ngood])) {
          ndev3++;
        }

        /* Slope calculations */

        if (fabs(pCurFilestarimage->magcal_magdep - full_med) < (3.0 * full_rms)) {
          /* Population 10 */
          pVectorStore->mvector[PARAMETERVECTOR_SLOPE_ALL][nslope_all].magcal_magdep = pCurFilestarimage->magcal_magdep;
          pVectorStore->mvector[PARAMETERVECTOR_SLOPE_ALL][nslope_all].Date = jd2ep(pCurFilestarimage->Date);
          nslope_all++;
          if ( pCurFilestarimage->Date < 2436934.50) {  /* Jan 1, 1960 */
            pVectorStore->mvector[PARAMETERVECTOR_SLOPE_60][nslope_60].magcal_magdep = pCurFilestarimage->magcal_magdep;
            pVectorStore->mvector[PARAMETERVECTOR_SLOPE_60][nslope_60].Date = jd2ep(pCurFilestarimage->Date);
            nslope_60++;
          }
        }


#if 0
        fprintf(consoleHandle,"plate %s%05d magcal_magdep %f mag_median %f errreal %f difference %f",
                pCurSumstarimage->series,
                pCurFilestarimage->plateNumber,
                pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][temp_ngood],
                mag_median,
                pVectorStore->vector[PARAMETERVECTOR_ERRREAL][temp_ngood],
                fabs(pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][temp_ngood]-mag_median) - (2 * pVectorStore->vector[PARAMETERVECTOR_ERRREAL][temp_ngood]));
#endif
        if (fabs(pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP][temp_ngood]-mag_median) > (2 * pVectorStore->vector[PARAMETERVECTOR_ERRREAL][temp_ngood])) {
#if 0
          fprintf(consoleHandle," SELECTED");
#endif
          ndev2++;
        }
#if 0
        fprintf(consoleHandle,"\n");
#endif



        temp_ngood++;
      }

    }

    if ((errorFlag == 0) && (pFileCommon->provideWebSummary == 0) && (pTranCommon->full_ngood >= MIN_GOODSTARS)) {
      if (pTranCommon->full_ngood != temp_ngood) {
        fprintf(consoleHandle,"ERROR: pTranCommon->full_ngood %d not temp_ngood %d in line %d\n",pTranCommon->full_ngood,temp_ngood,__LINE__);
        exit(1);
      }
      /* remove galaxies, noises, and mis-matches by doing the correlations */
      magvslimitingcorr = corr2(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP],pVectorStore->vector[PARAMETERVECTOR_LIMITING_MAG]);
      magvsracorr       = corr2(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP],pVectorStore->vector[PARAMETERVECTOR_RA]);
      magvsdeccorr      = corr2(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP],pVectorStore->vector[PARAMETERVECTOR_DEC]);

      /* rms of coord */
      rmsdradrms2       = findstd(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_DRADRMS2]);
      rarms             = findstd(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_RA]);
      decrms            = findstd(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_DEC]);

      /* add a parameter to indicate if more than 1 adjacent burst3/dip3/dev3 points */
      adjacentburstdip = 0;
      if (num_iburst3 > 1) {
        iadj = 0;
        for (iburst_index = 0; iburst_index < num_iburst3;iburst_index++) {
          if ((pVectorStore->ivector[PARAMETERVECTOR_IBURST3][iburst_index+1] - pVectorStore->ivector[PARAMETERVECTOR_IBURST3][iburst_index]) == 1) {
            iadj++;
          }
        }
        adjacentburstdip  += iadj;
      }
      if (num_idip3 > 1) {
        iadj = 0;
        for (idip_index = 0; idip_index < num_idip3;idip_index++) {
          if ((pVectorStore->ivector[PARAMETERVECTOR_IDIP3][idip_index+1] - pVectorStore->ivector[PARAMETERVECTOR_IDIP3][idip_index]) == 1) {
            iadj++;
          }
        }
        adjacentburstdip  += iadj;
      }
      if (num_iburst4 > 1) {
        iadj = 0;
        for (iburst_index = 0; iburst_index < num_iburst4;iburst_index++) {
          if ((pVectorStore->ivector[PARAMETERVECTOR_IBURST4][iburst_index+1] - pVectorStore->ivector[PARAMETERVECTOR_IBURST4][iburst_index]) == 1) {
            iadj++;
          }
        }
        adjacentburstdip  += iadj;
      }
      if (num_idip4 > 1) {
        iadj = 0;
        for (idip_index = 0; idip_index < num_idip4;idip_index++) {
          if ((pVectorStore->ivector[PARAMETERVECTOR_IDIP4][idip_index+1] - pVectorStore->ivector[PARAMETERVECTOR_IDIP4][idip_index]) == 1) {
            iadj++;
          }
        }
        adjacentburstdip  += iadj;
      }
      /*  add another parameter to indicate if more than 1 adjacent burst3-4/dip3-4 points, only for either burst or dip, and
       * only continuous >5 points count
       */
      kadj2a = 0;
      kadj2b = 0;
      kadj2c = 0;
      kadj2d = 0;
      if (num_iburst3 > 1) {
        kadj2a = countAdjacent(num_iburst3,pVectorStore->ivector[PARAMETERVECTOR_IBURST3],6);
      }
      if (num_idip3 > 1) {
        kadj2b = countAdjacent(num_idip3,pVectorStore->ivector[PARAMETERVECTOR_IDIP3],6);
      }
      if (num_iburst4 > 1) {
        kadj2c = countAdjacent(num_iburst4,pVectorStore->ivector[PARAMETERVECTOR_IBURST4],6);
      }
      if (num_idip4 > 1) {
        kadj2d = countAdjacent(num_idip4,pVectorStore->ivector[PARAMETERVECTOR_IDIP4],6);
      }
      adjacentburstdip2 = kadj2a;
      if (adjacentburstdip2 < kadj2b) {
        adjacentburstdip2 = kadj2b;
      }
      if (adjacentburstdip2 < kadj2c) {
        adjacentburstdip2 = kadj2c;
      }
      if (adjacentburstdip2 < kadj2d) {
        adjacentburstdip2 = kadj2d;
      }

      kadj2a = 0;
      kadj2b = 0;
      /*  for more than 8 points in iburst5/idip5 */


      if (num_iburst5 > 1) {
        kadj2c = countAdjacent(num_iburst5,pVectorStore->ivector[PARAMETERVECTOR_IBURST5],8);
      }
      if (num_idip5 > 1) {
        kadj2 = countAdjacent(num_idip5,pVectorStore->ivector[PARAMETERVECTOR_IDIP5],8);
        if (kadj2 > kadj2c) {
          kadj2c = kadj2;
        }
      }
      adjacentburstdip3 = kadj2c;
      /* check systematic differences between dnb plates and non-dn plates */


      for (index = 0; index < pTranCommon->full_ngood;index++) {
        if ((pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index] == min_local) ||
            (pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index] == max_local)) {
          pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index] = mag_median;
        }
      }
      span = 4000 * (pTranCommon->full_ngood+3999)/10000;
      smooth(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_DATEREAL],pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3],span,"lowess",pVectorStore->vector[PARAMETERVECTOR_MAGS]);
      span = 8000 * (pTranCommon->full_ngood+7999)/10000;
      smooth(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_XVALS],pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3],span,"lowess",pVectorStore->vector[PARAMETERVECTOR_MAGSB]);
      smooth(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_XVALS],pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3],10,"sgolay",pVectorStore->vector[PARAMETERVECTOR_MAGSC]);
      smooth(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_XVALS],pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3],15,"loess",pVectorStore->vector[PARAMETERVECTOR_MAGSD]);


      /* Calculate lightcurverms1 */
      trendCount = 0;
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3]);
      stdvalue = findstd(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3]);

      for (index = 0; index < pTranCommon->full_ngood;index++) {
        if (fabs(pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index] - meanvalue) < (5 *stdvalue)) {
          pVectorStore->vector[PARAMETERVECTOR_IDTREND][trendCount] = pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index];
          trendCount++;
        }
      }
      lightcurverms1 = findstd(trendCount,pVectorStore->vector[PARAMETERVECTOR_IDTREND]);

      /* Now calculate lightcurverms2 */
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] = pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index] - pVectorStore->vector[PARAMETERVECTOR_MAGS][index];
      }
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]);
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] = pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] - meanvalue;
      }
      trendCount = 0;
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]); /* This should be zero! */
      stdvalue = findstd(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]);
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        if (fabs(pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] - meanvalue) < (5 *stdvalue)) {
          pVectorStore->vector[PARAMETERVECTOR_IDTREND][trendCount] = pVectorStore->vector[PARAMETERVECTOR_SCATTER][index];
          trendCount++;
        }
      }
      lightcurverms2 = findstd(trendCount,pVectorStore->vector[PARAMETERVECTOR_IDTREND]);

      /* Now calculate lightcurverms3 */
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] = pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index] - pVectorStore->vector[PARAMETERVECTOR_MAGSB][index];
      }
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]);
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] = pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] - meanvalue;
      }
      trendCount = 0;
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]); /* This should be zero! */
      stdvalue = findstd(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]);
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        if (fabs(pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] - meanvalue) < (5 *stdvalue)) {
          pVectorStore->vector[PARAMETERVECTOR_IDTREND][trendCount] = pVectorStore->vector[PARAMETERVECTOR_SCATTER][index];
          trendCount++;
        }
      }
      lightcurverms3 = findstd(trendCount,pVectorStore->vector[PARAMETERVECTOR_IDTREND]);

      /* Now calculate lightcurverms4 */
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] = pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index] - pVectorStore->vector[PARAMETERVECTOR_MAGSC][index];
      }
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]);
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] = pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] - meanvalue;
      }
      trendCount = 0;
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]); /* This should be zero! */
      stdvalue = findstd(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]);
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        if (fabs(pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] - meanvalue) < (5 *stdvalue)) {
          pVectorStore->vector[PARAMETERVECTOR_IDTREND][trendCount] = pVectorStore->vector[PARAMETERVECTOR_SCATTER][index];
          trendCount++;
        }
      }
      lightcurverms4 = findstd(trendCount,pVectorStore->vector[PARAMETERVECTOR_IDTREND]);

      /* Now calculate lightcurverms5 */
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] = pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP3][index] - pVectorStore->vector[PARAMETERVECTOR_MAGSD][index];
      }
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]);
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] = pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] - meanvalue;
      }
      trendCount = 0;
      meanvalue = findmean(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]); /* This should be zero! */
      stdvalue = findstd(pTranCommon->full_ngood,pVectorStore->vector[PARAMETERVECTOR_SCATTER]);
      for (index = 0; index < pTranCommon->full_ngood;index++) {
        if (fabs(pVectorStore->vector[PARAMETERVECTOR_SCATTER][index] - meanvalue) < (5 *stdvalue)) {
          pVectorStore->vector[PARAMETERVECTOR_IDTREND][trendCount] = pVectorStore->vector[PARAMETERVECTOR_SCATTER][index];
          trendCount++;
        }
      }
      lightcurverms5 = findstd(trendCount,pVectorStore->vector[PARAMETERVECTOR_IDTREND]);


      qsort((void*)pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP4],pTranCommon->full_ngood,sizeof(double),dcmp);

      if (CalcMedianAndRMS(pTranCommon->full_ngood,minGoodStars,pVectorStore->vector[PARAMETERVECTOR_ERRREAL],&errreal_med,&errreal_rms,0,3.0,0) >= 0) {
        if (errreal_med != 0.0) {
          rms_factor = findstd(pTranCommon->full_ngood-2,&pVectorStore->vector[PARAMETERVECTOR_MAGCAL_MAGDEP4][1]);
          rms_factor = rms_factor/errreal_med;
        }
      }
      dmagcatalog = pCurStarEntry->Stdmag  - full_med;


    }
    if (nslope_all >= MIN_SLOPESTARS) {
      BevingtonFit(pVectorStore->mvector[PARAMETERVECTOR_SLOPE_ALL],nslope_all,&slope_all,&slope_all_err,&offset);
    }
    if (nslope_60 >= MIN_SLOPESTARS) {
      BevingtonFit(pVectorStore->mvector[PARAMETERVECTOR_SLOPE_60],nslope_60,&slope_60,&slope_60_err,&offset);
    }



    if (errorFlag == 0) {
      outStars++;
      (*summaryCount)++;
      pTranCommon->nearbyREFflag = 0;
      if (pFileCommon->provideWebSummary == 0) {
        if (pTarget != NULL) {
          if ((pFileCommon->enableTransientSearch != 0) && (GetREFType(pCurStarEntry->REFNumber) == REF_TYPE_DASCH)) {
            pNearestCatalogStar = &pTarget->nearestCatalogStar;
            closestDistanceArcsec = (1.5*(min_FWHM_WORLD*3600.0/2.0))+dradrms;
            if (closestDistanceArcsec < 18) {
              closestDistanceArcsec = 18;
            }
            if ((pNearestCatalogStar->nearestREFflagGoodMag >= 0) &&
                (strlen(pNearestCatalogStar->nearestREFGoodMag) > 0) &&
                (pNearestCatalogStar->nearestREFarcsecGoodMag < 9990.0) &&
                (closestDistanceArcsec > pNearestCatalogStar->nearestREFarcsecGoodMag) &&
                (pNearestCatalogStar->nearestREFmag < 16.0)) {
              pTranCommon->nearbyREFflag = 1;
            } else if ((pNearestCatalogStar->nearestREFflag >= 0) &&
                       (strlen(pNearestCatalogStar->nearestREF) > 0) &&
                       (pNearestCatalogStar->nearestREFarcsec < 9990.0) &&
                       (closestDistanceArcsec > pNearestCatalogStar->nearestREFarcsec) &&
                       (pNearestCatalogStar->nearestREFmag < 16.0)) {
              pTranCommon->nearbyREFflag = 1;
            }
          }


          FindNearestGalaxy(pGalaxyCommon,pCurStarEntry->ra,pCurStarEntry->dec,&pTarget->nearestCatalogStar,nearbyObjects,NULL);
        } else {
          FindNearestGalaxy(pGalaxyCommon,pCurStarEntry->ra,pCurStarEntry->dec,pNullCatalogStar,nearbyObjects,NULL);
        }

        if (pFileCommon->enableTransientSearch != 0) {
          pTranCommon->minGoodStars = minGoodStars;
          if (ProcessLongTransientCandidates(pFileCommon,pTranCommon,pVectorStore,pMagnitudeTable,pGalaxyCommon,pCurStarEntry,nearbyObjects,AFLAGSMASK1,QUALITYMASK) == 0) {

            if (pTranCommon->npoints > MAX_TC_POINTS) {
              pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_TOOMANYPOINTS);
            } else {
              ProcessTransientCandidates(pFileCommon,pTranCommon,pVectorStore,pMagnitudeTable,pGalaxyCommon,pCurStarEntry,nearbyObjects,AFLAGSMASK1,QUALITYMASK);
            }
          }
        }
        if ((pFileCommon->enableTransientSearch == 0) ||
            (pTranCommon->peakDays > 0) ||
            (pTranCommon->peakEvaluation != 0)) {

          if (pTarget == NULL) {
            if (strlen(nearbyObjects) > 0) {
              sprintf(nearbyObjects2, "\"%s\"", nearbyObjects);
            } else {
              if (pFileCommon->extraColumnFlag == 0) {
                sprintf(nearbyObjects2,"\\N");
              } else {
                nearbyObjects2[0] = 0;
              }
            }

            fprintf(summary_handle,"%lld\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.7f\t%d\t%d\t%.4f\t%.7f\t%d\t%.4f\t%.4f\t%4f\t%f\t%f\t%f\t%d\t%d\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%s\t%s",
                    pCurStarEntry->REFNumber,
                    PASSBIT_GSC232,
                    versionId,
                    dradrms,
                    dradrmsB,
                    max_drad,
                    max_dradB,
                    jd2ep(min_date),
                    jd2ep(max_date),
                    pTranCommon->npoints,
                    min_local,
                    max_local,
                    range_local,
                    min_local2,
                    max_local2,
                    range_local2,
                    full_med,
                    full_rms,
                    pTranCommon->full_ngood,
                    full_ngoodB,
                    pTranCommon->clip_med,
                    pTranCommon->clip_rms,
                    pTranCommon->clip_ngood,
                    median_iso,
                    max_iso-min_iso,
                    error_bar_factor,
                    Malmquist_factor,
                    Malmquist_factorB,
                    Damon_factor,
                    Sextractor_Blend,
                    nblend,
                    nNonDamonBlue,
                    medNonDamonBlue,
                    nonDamonBluerms,
                    nDamonBlue,
                    medDamonBlue,
                    damonBlueRms,
                    magvslimitingcorr,
                    magvsracorr,
                    magvsdeccorr,
                    rmsdradrms2,
                    rarms,
                    decrms,
                    nburst,
                    nburst2,
                    nburst3,
                    nburst4,
                    ndip,
                    ndip2,
                    ndip3,
                    ndip4,
                    ndev3,
                    ndev2,
                    adjacentburstdip,
                    adjacentburstdip2,
                    adjacentburstdip3,
                    lightcurverms1,
                    lightcurverms2,
                    lightcurverms3,
                    lightcurverms4,
                    lightcurverms5,
                    slope_all,
                    slope_all_err,
                    nslope_all,
                    slope_60,
                    slope_60_err,
                    nslope_60,
                    rms_factor,
                    dmagcatalog,
                    keplerFieldString[keplerField],
                    nearbyObjects2);

          } else {

            fprintf(summary_handle,"%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.7f\t%d\t%d\t%.4f\t%.7f\t%d\t%.4f\t%.4f\t%4f\t%f\t%f\t%f\t%d\t%d\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.2f\t%.2f\t%.6f\t%.7f\t%d\t%d\t%d\t%.5f\t%.5f\t%s\t%d\t%s",
                    pCurStarEntry->REF,  /* REF */
                    dradrms,             /* drad */
                    dradrmsB,            /* dradB */
                    max_drad,            /* max_drad */
                    max_dradB,           /* max_dradB */
                    jd2ep(min_date),     /* yrbegin */
                    jd2ep(max_date),     /* yrend */
                    pTranCommon->npoints,             /* npoints */
                    min_local,           /* min_local */
                    max_local,           /* max_local */
                    range_local,         /* range_local */
                    min_local2,          /* min_local2 */
                    max_local2,          /* max_local2 */
                    range_local2,        /* range_local2 */
                    full_med,            /* median_local */
                    full_rms,            /* rms_local */
                    pTranCommon->full_ngood,          /* ngood */
                    full_ngoodB,         /* ngoodB */
                    pTranCommon->clip_med,            /* clip_median_local */
                    pTranCommon->clip_rms,            /* clip_rms_local */
                    pTranCommon->clip_ngood,          /* clip_ngood */
                    median_iso,          /* median_iso */
                    max_iso-min_iso,     /* range_iso */
                    error_bar_factor,    /* error_bar_factor */
                    Malmquist_factor,    /* Malmquist_factor */
                    Malmquist_factorB,   /* Malmquist_factorB */
                    Damon_factor,        /* Damon_factor */
                    Sextractor_Blend,    /* Sextractor_Blend */
                    nblend,              /* nblend */
                    nNonDamonBlue,       /* nNonDamonBlue */
                    medNonDamonBlue,     /* medNonDamonBlue */
                    nonDamonBluerms,     /* nonDamonBluerms */
                    nDamonBlue,          /* nDamonBlue */
                    medDamonBlue,        /* medDamonBlue */
                    damonBlueRms,        /* damonBlueRms */
                    magvslimitingcorr,   /* magvslimitingcorr */
                    magvsracorr,         /* magvsracorr */
                    magvsdeccorr,        /* magvsdeccorr */
                    rmsdradrms2,         /* rmsdradrms2 */
                    rarms,               /* rarms */
                    decrms,              /* decrms */
                    nburst,              /* nburst */
                    nburst2,             /* nburst2 */
                    nburst3,             /* nburst3 */
                    nburst4,             /* nburst4 */
                    ndip,                /* ndip */
                    ndip2,               /* ndip2 */
                    ndip3,               /* ndip3 */
                    ndip4,               /* ndip4 */
                    ndev3,               /* ndev3 */
                    ndev2,               /* ndev2 */
                    adjacentburstdip,    /* adjacentburstdip */
                    adjacentburstdip2,   /* adjacentburstdip2 */
                    adjacentburstdip3,   /* adjacentburstdip3 */
                    lightcurverms1,      /* lightcurverms1 */
                    lightcurverms2,      /* lightcurverms2 */
                    lightcurverms3,      /* lightcurverms3 */
                    lightcurverms4,      /* lightcurverms4 */
                    lightcurverms5,      /* lightcurverms5 */
                    slope_all,           /* slope_all */
                    slope_all_err,       /* slope_all_err */
                    nslope_all,          /* nslope_all */
                    slope_60,            /* slope_60 */
                    slope_60_err,        /* slope_60_err */
                    nslope_60,           /* nslope_60 */
                    rms_factor,          /* rms_factor */
                    dmagcatalog,         /* dmagcatalog */
                    pTarget->Stdmag,     /* Stdmag */
                    pTarget->color,      /* color */
                    pTarget->ra,         /* ra */
                    pTarget->dec,        /* declination */
                    pTarget->MAGFlag,    /* MAGFlag */
                    pTarget->class,      /* gscclass */
                    pTarget->VFlag,      /* VFlag */
                    pTarget->RaPM,       /* RaPM */
                    pTarget->DecPM,      /* DecPM */
                    keplerFieldString[keplerField],  /* keplerField */
                    versionId,           /* starsummary2.versionId */
                    nearbyObjects);

          }
          if (pFileCommon->extraColumnFlag != 0) {
            lon = pCurStarEntry->ra;
            lat = pCurStarEntry->dec;
            CheckAuthorization("none", pCurStarEntry->ra,pCurStarEntry->dec,&releaseField);
            wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
            fprintf(summary_handle,"\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t%f\t%d\n",
                    lat,
                    pTranCommon->nearbyREFflag,
                    releaseField,
                    pTranCommon->peakDays,
                    pTranCommon->peakYear,
                    pTranCommon->peakMag,
                    pTranCommon->peakSlope,
                    pTranCommon->peakRMS,
                    pTranCommon->peakCount,
                    pTranCommon->peakUpperCount,
                    pTranCommon->peakNumber,
                    pTranCommon->peakMaxMag,
                    pTranCommon->peakMinMag,
                    pTranCommon->peakEventCount,
                    pTranCommon->peakDradRMS,
                    pTranCommon->peakEventCount2,
                    pTranCommon->peakEventCount3,
                    pTranCommon->peakExtra,
                    pTranCommon->peakCountNF,
                    pTranCommon->peakCountWF,
                    pTranCommon->peakOutside,
                    pTranCommon->peakOutsideNF,
                    pTranCommon->peakOutsideWF,
                    pTranCommon->peakDefectCount,
                    pTranCommon->peakRA,
                    pTranCommon->peakDec,
                    pTranCommon->peakNearbyDistance,
                    pTranCommon->peakMultipleCount,
                    pTranCommon->gsc_bin_index,
                    pTranCommon->peakDradRMS2,
                    pTranCommon->peakDradRMS3,
                    pTranCommon->peakBadColorCount,
                    pTranCommon->peakLimitingYears,
                    pTranCommon->peakLimitingPoints,
                    pTranCommon->peakEvaluation,
                    pTranCommon->peakNoDefectMag,
                    pTranCommon->peakLongOutburst);
          }

          if (debugMode) {
            fflush(summary_handle);
          }
#if 0
          fprintf(consoleHandle,"index\tvariable\tvalue\n");
          fprintf(consoleHandle,"-----\t--------\t-----\n");
          fprintf(consoleHandle,"1\tngood\t%d\n",pTranCommon->full_ngood);
          fprintf(consoleHandle,"2\trange_local\t%6.3f\n",range_local);
          fprintf(consoleHandle,"3\trange_local2\t%6.3f\n",range_local2);
          fprintf(consoleHandle,"4\tnblend\t%d\n",nblend);
          fprintf(consoleHandle,"5\tnNonDamonBlue\t%d\n",nNonDamonBlue);
          fprintf(consoleHandle,"6\tmedNonDamonBlue\t%6.3f\n",medNonDamonBlue);
          fprintf(consoleHandle,"7\tnonDamonBluerms\t%6.3f\n",nonDamonBluerms);
          fprintf(consoleHandle,"8\tnDamonBlue\t%d\n",nDamonBlue);
          fprintf(consoleHandle,"9\tmedDamonBlue\t%6.3f\n",medDamonBlue);
          fprintf(consoleHandle,"10\tdamonBlueRms\t%6.3f\n",damonBlueRms);
          fprintf(consoleHandle,"11\tmagvslimitingcorr\t%6.3f\n",magvslimitingcorr);
          fprintf(consoleHandle,"12\tmagvsracorr\t%6.3f\n",magvsracorr);
          fprintf(consoleHandle,"13\tmagvsdeccorr\t%6.3f\n",magvsdeccorr);
          fprintf(consoleHandle,"14\trmsdradrms2\t%6.3f\n",rmsdradrms2);
          fprintf(consoleHandle,"15\trarms\t%7.4f\n",rarms);
          fprintf(consoleHandle,"16\tdecrms\t%7.4f\n",decrms);
          fprintf(consoleHandle,"17\tnburst\t%d\n",nburst);
          fprintf(consoleHandle,"18\tnburst2\t%d\n",nburst2);
          fprintf(consoleHandle,"19\tnburst3\t%d\n",nburst3);
          fprintf(consoleHandle,"20\tnburst4\t%d\n",nburst4);
          fprintf(consoleHandle,"21\tndip\t%d\n",ndip);
          fprintf(consoleHandle,"22\tndip2\t%d\n",ndip2);
          fprintf(consoleHandle,"23\tndip3\t%d\n",ndip3);
          fprintf(consoleHandle,"24\tndip4\t%d\n",ndip4);
          fprintf(consoleHandle,"25\tndev3\t%d\n",ndev3);
          fprintf(consoleHandle,"26\tndev2\t%d\n",ndev2);
          fprintf(consoleHandle,"27\tadjacentburstdip\t%d\n",adjacentburstdip);
          fprintf(consoleHandle,"28\tadjacentburstdip2\t%d\n",adjacentburstdip2);
          fprintf(consoleHandle,"29\tadjacentburstdip3\t%d\n",adjacentburstdip3);
          fprintf(consoleHandle,"30\tlightcurverms1\t%6.3f\n",lightcurverms1);
          fprintf(consoleHandle,"31\tlightcurverms2\t%6.3f\n",lightcurverms2);
          fprintf(consoleHandle,"32\tlightcurverms3\t%6.3f\n",lightcurverms3);
          fprintf(consoleHandle,"33\tlightcurverms4\t%6.3f\n",lightcurverms4);
          fprintf(consoleHandle,"34\tlightcurverms5\t%6.3f\n",lightcurverms5);
          fprintf(consoleHandle,"35\trms_factor\t%6.3f\n",rms_factor);
          fprintf(consoleHandle,"36\tdmagcatalog\t%6.3f\n",dmagcatalog);
#endif
        }
      } /* provideWebSummary */
    } /* errorFlag == 0 */
  }

  return(outStars);
}


void
InitCatalogAccess(PFILECOMMON pFileCommon)
{
  char catalogname[MAX_FILENAME];
  char indexname[MAX_FILENAME];
  char* catalogdir;
  char* slashPtr;
  char *dotPtr;

  if ((pFileCommon->catalogFD >= 0) ||
      (pFileCommon->indexFD >= 0)) {
    printf("ERROR: InitCatalogAccess has already been called\n");
    return;
  }
  catalogdir = getenv("DASCH_CATALOG");
  if (catalogdir == NULL) {
    printf("ERROR: DASCH_CATALOG is not defined\n");
    return;
  }
  strcpy(catalogname,catalogdir);

  if (pFileCommon->catalogNumber == CATALOG_KEPLER) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"kepler.dat");
  }

  if (pFileCommon->catalogNumber == CATALOG_APASS) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"apass.dat");
  }

  if (pFileCommon->catalogNumber == CATALOG_ATLAS) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"atlas.dat");
  }

  if (pFileCommon->catalogNumber == CATALOG_GAIA) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"gaiadr2.dat");
  }

  if (pFileCommon->catalogNumber == CATALOG_EXPERIMENTAL) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"experimental.dat");
  }

  strcpy(indexname,catalogname);
  dotPtr = strrchr(indexname,'.');
  if (dotPtr != NULL) {
    *dotPtr = 0;
  }
  strcat(indexname,".idx");



  pFileCommon->catalogFD = open(catalogname,O_RDONLY,S_IRWXU|S_IRGRP);
  if (pFileCommon->catalogFD < 0) {
    printf("ERROR Could not open catalog file %s\n",catalogname);
    return;
  }
  pFileCommon->indexFD = open(indexname,O_RDONLY,S_IRWXU|S_IRGRP);
  if (pFileCommon->indexFD < 0) {
    printf("ERROR Could not open catalog index file %s\n",indexname);
    close(pFileCommon->catalogFD);
    pFileCommon->catalogFD = -1;

    return;
  }

}


static int
GscImageXCompare(const void *first, const void *second)
{
  long long coverage_bin_indexFirst = ((PGSCIMAGE)first)->REFNumber;
  long long coverage_bin_indexSecond = ((PGSCIMAGE)second)->REFNumber;
  if (coverage_bin_indexFirst > coverage_bin_indexSecond) {
    return(1);
  } else if (coverage_bin_indexFirst < coverage_bin_indexSecond) {
    return(-1);
  } else {
    return(0);
  }
}


static int
ReadCatalogBin(PFILECOMMON pFileCommon,int gsc_bin_index)
{
  STARINDEX curStarIndex;
  PSTARINDEX pCurStarIndex = &curStarIndex;
  int readItems;
  int readBytes;
  int itemIndex;
  PGSCIMAGE pCurGscImageX;
  PGSCIMAGE pCurGscImage;

  if (gsc_bin_index == pFileCommon->gsc_bin_index) {
    return(0);
  }
  if ((pFileCommon->catalogFD < 0) ||
      (pFileCommon->indexFD < 0)) {
    return(-1);
  }
  lseek(pFileCommon->indexFD,gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
  readBytes = read(pFileCommon->indexFD,pCurStarIndex,sizeof(STARINDEX));
  readItems = readBytes/sizeof(STARINDEX);
  if (readItems != 1) {
    printf("ERROR reading star index file errno: %d %s gsc_bin_index %d\n",
           errno,
           strerror(errno),
           gsc_bin_index);
    return(-1);
  } else {

    if (pCurStarIndex->numStars > pFileCommon->catalogTableAlloc) {
      pFileCommon->catalogTableAlloc = pCurStarIndex->numStars+1000;
      pCurGscImageX = realloc(pFileCommon->catalogTable,pFileCommon->catalogTableAlloc*sizeof(GSCIMAGE));
      if (pCurGscImageX == NULL) {
        fprintf(
          stderr,
          "ERROR: failed to realloc pFileCommon->catalogTable of size %zu\n",
          pFileCommon->catalogTableAlloc * sizeof(GSCIMAGE)
        );
        return(-1);
      }
      pFileCommon->catalogTable = pCurGscImageX;
      pCurGscImageX = NULL;
    }
    lseek(pFileCommon->catalogFD,pCurStarIndex->offset,SEEK_SET);
    readBytes = read(pFileCommon->catalogFD,pFileCommon->catalogTable,sizeof(GSCIMAGE)*pCurStarIndex->numStars);
    readItems = readBytes/sizeof(GSCIMAGE);
    if (readItems != pCurStarIndex->numStars) {
      fprintf(stderr,"ERROR reading gsc catalog file at index %d\n",gsc_bin_index);
      return(-1);
    }
  }
  /* At this point, we are going to convert the ASCII REF into REFNumber and sort the array in increasing REFNumber (Not necessary after Aug 12, 2013) */
  for (itemIndex = 0; itemIndex < readItems; itemIndex++) {
    pCurGscImageX = &pFileCommon->catalogTable[itemIndex];
    pCurGscImage = (PGSCIMAGE)(&pFileCommon->catalogTable[itemIndex]);
    pCurGscImageX->REFNumber = pCurGscImage->REFNumber;
  }
  qsort((void *)pFileCommon->catalogTable,readItems,sizeof(GSCIMAGE),GscImageXCompare);
  pFileCommon->gsc_bin_index = gsc_bin_index;
  pFileCommon->catalogTableSize = readItems;
  return(0);
}


static PGSCIMAGE
GetCatalogRecord(PFILECOMMON pFileCommon, long long REFNumber, int gsc_bin_index)
{
  PGSCIMAGE pCurGscImageX;
  int itemIndex;

  if (REFNumber == 0) {
    return(NULL);
  }
  if (ReadCatalogBin(pFileCommon,gsc_bin_index) != 0) {
    return(NULL);
  }
  for (itemIndex = 0; itemIndex < pFileCommon->catalogTableSize; itemIndex++) {
    pCurGscImageX = &pFileCommon->catalogTable[itemIndex];
    if (pCurGscImageX->REFNumber == REFNumber) {
      return(pCurGscImageX);
    }
  }

  pFileCommon->counterBlock.staleREFCount++;
#if 0
  printf("ERROR: Failed to find in gsc_bin_index %d REFNumber %lld\n",gsc_bin_index,REFNumber);
#endif
  return(NULL);
}


/* The gsc bin list is used for debugging or quick repeat runs in a particular part of the sky
   For every bin index derived from the ra and dec, use the one ahead of it and the one behind it
*/
int
InitGscBinList(PFILECOMMON pFileCommon, char *gscbinfile) {
  FILE *gscbinhandle;
  char *inBuffer;
  char inLine[MAX_GSC_FILE_BUFFER];
  int lineLen;
  int lineCounter = 0;
  int* gscBinIndexList = NULL;
  int gscBinIndexCount = 0;
  int nvals;
  int gsc_bin_index;
  int index2;

  if (pFileCommon->gscBinIndexCount > 0) {
    if (pFileCommon->gscBinIndexList != NULL) {
      free(pFileCommon->gscBinIndexList);
      pFileCommon->gscBinIndexList = NULL;
    }
    pFileCommon->gscBinIndexCount = 0;
  }
  if (gscbinfile[0] == 0) {
    return(0);
  }

  gscbinhandle = fopen(gscbinfile,"rt");
  if (gscbinhandle == NULL) {
    printf("ERROR InitGscBinList: failed to open the gsc bin file %s\n",gscbinfile);
    return(-1);
  }

  /* First find out how many thumbnails we have */
  while (1) {
    inBuffer = fgets(inLine,MAX_GSC_FILE_BUFFER,gscbinhandle);
    if (inBuffer == NULL) {
      break;
    }
    lineLen = strlen(inBuffer);
    /* Trim off the carriage return */
    if (inBuffer[lineLen-1] == 10) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }

    lineCounter++;
  }

  fclose(gscbinhandle);

  if (lineCounter == 0) {
    printf("ERROR InitGscBinList: no lines are in %s\n", gscbinfile);
    return(-1);
  }

  gscBinIndexCount = 0;
  gscBinIndexList = (int *)calloc(GSC_NEIGHBORS*(lineCounter+2),sizeof(int));
  if (gscBinIndexList == NULL) {
    printf("ERROR InitGscBinList: failed to allocate the gscBinIndexList of size %d\n",lineCounter+1);
    return(-1);
  }
  /* Add in two dummy entries to simplify logic checks */
  gscBinIndexList[gscBinIndexCount] = -1;
  gscBinIndexCount++;
  gscBinIndexList[gscBinIndexCount] = -1;
  gscBinIndexCount++;


  gscbinhandle = fopen(gscbinfile,"rt");
  if (gscbinhandle == NULL) {
    free(gscBinIndexList);
    printf("Could not open list file %s\n",gscbinfile);
    return(-1);
  }
  while (1) {
    inBuffer = fgets(inLine,MAX_GSC_FILE_BUFFER,gscbinhandle);
    if (inBuffer == NULL) {
      break;
    }
    lineLen = strlen(inBuffer);
    /* Trim off the carriage return */
    if (inBuffer[lineLen-1] == 10) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    /* Trim trailing spaces */
    while ((lineLen >= 0) && (inBuffer[lineLen-1] == ' ')) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }

    nvals = sscanf(inBuffer, "%d", &gsc_bin_index);

    if (nvals == 1) {
      if (((gsc_bin_index -1) != gscBinIndexList[gscBinIndexCount-1]) &&
          ((gsc_bin_index -1) != gscBinIndexList[gscBinIndexCount-2])) {
        gscBinIndexList[gscBinIndexCount] = gsc_bin_index-1;
        gscBinIndexCount++;
      }
      if (gsc_bin_index != gscBinIndexList[gscBinIndexCount-1]) {
        gscBinIndexList[gscBinIndexCount] = gsc_bin_index;
        gscBinIndexCount++;
      }
      gscBinIndexList[gscBinIndexCount] = gsc_bin_index+1;
      gscBinIndexCount++;

      if (gscBinIndexCount > GSC_NEIGHBORS*(lineCounter+2)) {
        printf("ERROR InitGscBinList: index %d exceeds %d in file %s\n",gscBinIndexCount,lineCounter,gscbinfile);
      }
    }
  }

  fclose(gscbinhandle);

  /* Now move every entry down by two */
  for (index2 = 0; index2 < (gscBinIndexCount-2); index2++) {
    gscBinIndexList[index2] = gscBinIndexList[index2+2];
  }

  gscBinIndexCount -= 2;

  /* Perform some checks */
  for (index2 = 1; index2 < gscBinIndexCount; index2++) {
    if (gscBinIndexList[index2] == gscBinIndexList[index2-1]) {
      printf("ERROR InitGscBinList: duplicate entry %d\n",gscBinIndexList[index2]);
      exit(1);
    } else if (gscBinIndexList[index2] < gscBinIndexList[index2-1]) {
      printf("ERROR InitGscBinList: decreasing entry %d and %d\n",gscBinIndexList[index2], gscBinIndexList[index2-1]);
      exit(1);
    }
  }

#if 0
  for (index2 = 0; index2 < gscBinIndexCount; index2++) {
    printf("%d\n",gscBinIndexList[index2]);
  }
#endif

  if (gscBinIndexCount == 0) {
    printf("ERROR InitGscBinList: no valid lines are in %s\n", gscbinfile);
    free(gscBinIndexList);
    return(-1);
  }

  pFileCommon->gscBinIndexCount = gscBinIndexCount;
  pFileCommon->gscBinIndexList = gscBinIndexList;
  printf("ERROR InitGscBinList: using reduced search from %s of length %d\n",gscbinfile,gscBinIndexCount);
  return(0);
}


/* Now sort by seriesId,plateNumber,NUMBER and catalogNumber */
static int
CatalogImageCompare(const void *first, const void *second)
{
  PPHOTSTARIMAGE pPhotFirst = (PPHOTSTARIMAGE)first;
  PPHOTSTARIMAGE pPhotSecond = (PPHOTSTARIMAGE)second;
  PFILESTARIMAGE pFirst = pPhotFirst->pFileStarImage;
  PFILESTARIMAGE pSecond = pPhotSecond->pFileStarImage;

  if (pFirst->seriesId > pSecond->seriesId) {
    return(1);
  } else if (pFirst->seriesId < pSecond->seriesId) {
    return(-1);
  } else {

    if (pFirst->plateNumber > pSecond->plateNumber) {
      return(1);
    } else if (pFirst->plateNumber < pSecond->plateNumber) {
      return(-1);
    } else {

      if (pFirst->NUMBER > pSecond->NUMBER) {
        return(1);
      } else if (pFirst->NUMBER < pSecond->NUMBER) {
        return(-1);
      } else {
        if (pFirst->solutionNumber > pSecond->solutionNumber) {
          return(1);
        } else if (pFirst->solutionNumber < pSecond->solutionNumber) {
          return(-1);
        } else {
#if 0
          if (((pFirst->REFNumber == 120132012267L) && (pSecond->REFNumber == 405144112674112L)) ||
              ((pSecond->REFNumber == 120132012267L) && (pFirst->REFNumber == 405144112674112L))) {
            printf("DUMP line %4d gsc_bin_index %d seriesId %2d %3s plateNumber %d NUMBER %d catalogNumber %d REFNumber %lld\n",
                   __LINE__,
                   pFirst->gsc_bin_index,
                   pFirst->seriesId,
                   GetSeriesString(pFirst->seriesId,0),
                   pFirst->plateNumber,
                   pFirst->NUMBER,
                   pFirst->catalogNumber,
                   pFirst->REFNumber);
          }
#endif
          /* Favor the ATLAS catalog */
          if ((pFirst->catalogNumber == CATALOG_ATLAS) && pSecond->catalogNumber != CATALOG_ATLAS) {
            return(1);
          } else if ((pFirst->catalogNumber != CATALOG_ATLAS) && pSecond->catalogNumber == CATALOG_ATLAS) {
            return(-1);
          } else {
            /* Favor the GAIA catalog */
            if ((pFirst->catalogNumber == CATALOG_GAIA) && pSecond->catalogNumber != CATALOG_GAIA) {
              return(1);
            } else if ((pFirst->catalogNumber != CATALOG_GAIA) && pSecond->catalogNumber == CATALOG_GAIA) {
              return(-1);
            } else {
              /* Now favor the APASS catalog */
              if ((pFirst->catalogNumber == CATALOG_APASS) && pSecond->catalogNumber != CATALOG_APASS) {
                return(1);
              } else if ((pFirst->catalogNumber != CATALOG_APASS) && pSecond->catalogNumber == CATALOG_APASS) {
                return(-1);
              } else {
                if (pFirst->catalogNumber > pSecond->catalogNumber) {
                  return(1);
                } else if (pFirst->catalogNumber < pSecond->catalogNumber) {
                  return(-1);
                } else {
                  return(0);
                }
              }
            }
          }
        }
      }
    }
  }
}


/* Reverse sort by count */
static int
SubstitutionCompare(const void *first, const void *second)
{
  PSUBSTITUTION pSubstitutionFirst = (PSUBSTITUTION)first;
  PSUBSTITUTION pSubstitutionSecond = (PSUBSTITUTION)second;
  if (pSubstitutionFirst->count < pSubstitutionSecond->count) {
    return(1);
  } else if (pSubstitutionFirst->count > pSubstitutionSecond->count) {
    return(-1);
  } else {
    return(0);
  }
}


void
CombineMultipleTables(
  PPHOTGLOBAL pPhotGlobal,
  PPHOTSTARIMAGE *ppMagnitudeTable,
  PFILESTARIMAGE *ppFileStarImageTable,
  PPHOTSTARIMAGE pMagnitudeTable1,
  PPHOTSTARIMAGE pMagnitudeTable2,
  PSUBSTITUTION* pSubstitutionTable,
  int numMagnitudes1,
  int numMagnitudes2,
  int catalogNumber1,
  int catalogNumber2,
  int* pAllocMagnitudes,
  int* pNumMagnitudes,
  int combineAlgorithm,
  int gsc_bin_index
) {
  int allocMagnitudes = *pAllocMagnitudes;
  int numMagnitudes = numMagnitudes1+numMagnitudes2;
  int curMagnitudes = 0;
  int magnitudeIndex0;
  int magnitudeIndex1;
  int magnitudeIndex2;
  PPHOTSTARIMAGE pSourceImage;
  PPHOTSTARIMAGE pTargetImage;
  PPHOTSTARIMAGE pStarImage0;
  PPHOTSTARIMAGE pStarImage1;
  PPHOTSTARIMAGE pStarImage2;

  PFILESTARIMAGE pFileStarImage0;
  PFILESTARIMAGE pFileStarImage1;
  PFILESTARIMAGE pFileStarImage2;
  PPHOTSTARIMAGE pMagnitudeTable = *ppMagnitudeTable;
  PFILESTARIMAGE pFileStarImageTable = *ppFileStarImageTable;
  PSUBSTITUTION substitutionTable = *pSubstitutionTable;
  PSUBSTITUTION pSubstitution;
  PSUBSTITUTION pSubstitution2;
  int substitutionCount = 0;
  int substitutionIndex;
  int substitutionIndex2;
  int bestChoice;

  *pNumMagnitudes = 0;
  if ((numMagnitudes1 == 0) && (numMagnitudes2 == 0)) {
    return;
  }
#if 0
  printf("numMagnitudes1 %d numMagnitudes2 %d, pMagnitudeTable %x allocMagnitudes %d\n",
         numMagnitudes1,numMagnitudes2,pMagnitudeTable,allocMagnitudes);
#endif

  if (allocMagnitudes < (numMagnitudes1+numMagnitudes2)) {
    if (pMagnitudeTable != NULL) {
      free(pMagnitudeTable);
    }
    allocMagnitudes = 2*(numMagnitudes1+numMagnitudes2+1000);
    pMagnitudeTable =  (PPHOTSTARIMAGE)calloc(allocMagnitudes,sizeof(PHOTSTARIMAGE));
    if (pMagnitudeTable == NULL) {
      printf("ERROR: failed to allocate pMagnitudeTable of size %d\n",allocMagnitudes);
      exit(1);
    }

    if (pFileStarImageTable != NULL) {
      free(pFileStarImageTable);
    }
    pFileStarImageTable =  (PFILESTARIMAGE)calloc(allocMagnitudes,sizeof(FILESTARIMAGE));
    if (pFileStarImageTable == NULL) {
      printf("ERROR: failed to allocate pFileStarImageTable of size %d\n",allocMagnitudes);
      exit(1);
    }

    if (substitutionTable != NULL) {
      free(substitutionTable);
    }

    substitutionTable =  (PSUBSTITUTION)calloc(allocMagnitudes,sizeof(SUBSTITUTION));
    if (substitutionTable == NULL) {
      printf("ERROR: failed to allocate substitutionTable of size %d\n",allocMagnitudes);
      exit(1);
    }
    *pSubstitutionTable = substitutionTable;
    *pAllocMagnitudes = allocMagnitudes;
    *ppMagnitudeTable = pMagnitudeTable;
    *ppFileStarImageTable = pFileStarImageTable;

  } else {
    memset(substitutionTable,0,allocMagnitudes*sizeof(SUBSTITUTION));
  }


  if (numMagnitudes2 == 0) {
    for (magnitudeIndex1 = 0; magnitudeIndex1 < numMagnitudes1; magnitudeIndex1++) {
      pTargetImage = &pMagnitudeTable[magnitudeIndex1];
      pFileStarImage0 = &pFileStarImageTable[magnitudeIndex1];
      pSourceImage = &pMagnitudeTable1[magnitudeIndex1];
      memcpy(pTargetImage,pSourceImage,sizeof(PHOTSTARIMAGE));
      memcpy(pFileStarImage0,pSourceImage->pFileStarImage,sizeof(FILESTARIMAGE));
      pTargetImage->pFileStarImage = pFileStarImage0;
      pTargetImage->pFileStarImage->catalogNumber = catalogNumber1;
      pTargetImage->pFileStarImage->versionId = pPhotGlobal->currentVersion;
    }
    *pNumMagnitudes = numMagnitudes;
    return;
  }
  if (numMagnitudes1 == 0) {
    for (magnitudeIndex2 = 0; magnitudeIndex2 < numMagnitudes2; magnitudeIndex2++) {
      pTargetImage = &pMagnitudeTable[magnitudeIndex2];
      pFileStarImage0 = &pFileStarImageTable[magnitudeIndex2];
      pSourceImage = &pMagnitudeTable2[magnitudeIndex2];
      memcpy(pTargetImage,pSourceImage,sizeof(PHOTSTARIMAGE));
      memcpy(pFileStarImage0,pSourceImage->pFileStarImage,sizeof(FILESTARIMAGE));
      pTargetImage->pFileStarImage = pFileStarImage0;
      pTargetImage->pFileStarImage->catalogNumber = catalogNumber2;
      pTargetImage->pFileStarImage->versionId = pPhotGlobal->currentVersion;
    }
    *pNumMagnitudes = numMagnitudes;
    return;
  }
  for (magnitudeIndex1 = 0; magnitudeIndex1 < numMagnitudes1; magnitudeIndex1++) {
    pTargetImage = &pMagnitudeTable[magnitudeIndex1];
    pFileStarImage0 = &pFileStarImageTable[magnitudeIndex1];
    pSourceImage = &pMagnitudeTable1[magnitudeIndex1];
    memcpy(pTargetImage,pSourceImage,sizeof(PHOTSTARIMAGE));
    memcpy(pFileStarImage0,pSourceImage->pFileStarImage,sizeof(FILESTARIMAGE));

    pTargetImage->pFileStarImage = pFileStarImage0;
    pTargetImage->pFileStarImage->catalogNumber = catalogNumber1;
    pTargetImage->pFileStarImage->versionId = pPhotGlobal->currentVersion;

  }
  for (magnitudeIndex2 = 0; magnitudeIndex2 < numMagnitudes2; magnitudeIndex2++) {
    pTargetImage = &pMagnitudeTable[magnitudeIndex2+numMagnitudes1];
    pFileStarImage0 = &pFileStarImageTable[magnitudeIndex2+numMagnitudes1];
    pSourceImage = &pMagnitudeTable2[magnitudeIndex2];
    memcpy(pTargetImage,pSourceImage,sizeof(PHOTSTARIMAGE));
    memcpy(pFileStarImage0,pSourceImage->pFileStarImage,sizeof(FILESTARIMAGE));

    pTargetImage->pFileStarImage = pFileStarImage0;
    pTargetImage->pFileStarImage->catalogNumber = catalogNumber2;
    pTargetImage->pFileStarImage->versionId = pPhotGlobal->currentVersion;
  }
  /* Now sort by gsc_bin_number, seriesId,plateNumber,NUMBER, and catalogNumber */
  qsort((void*)pMagnitudeTable,numMagnitudes,sizeof(PHOTSTARIMAGE),CatalogImageCompare);

  substitutionCount = 0;

  /* For all cases, substitute the APASS id for any non-APASS id */
  for (magnitudeIndex0 = 0; magnitudeIndex0 < numMagnitudes; magnitudeIndex0++) {
    pStarImage1 = &pMagnitudeTable[magnitudeIndex0];
    pFileStarImage1 = pStarImage1->pFileStarImage;
    if ((magnitudeIndex0+1) < numMagnitudes) {
      pStarImage2 =  &pMagnitudeTable[magnitudeIndex0+1];
      pFileStarImage2 = pStarImage2->pFileStarImage;
      if ((pFileStarImage1->solutionNumber == pFileStarImage2->solutionNumber) &&
          (pFileStarImage1->seriesId == pFileStarImage2->seriesId) &&
          (pFileStarImage1->plateNumber == pFileStarImage2->plateNumber) &&
          (pFileStarImage1->NUMBER == pFileStarImage2->NUMBER) &&
          (pFileStarImage1->catalogNumber != CATALOG_APASS) &&
          (pFileStarImage2->catalogNumber == CATALOG_APASS)) {
        if ((pFileStarImage1->REFNumber != pFileStarImage2->REFNumber) &&
            (GetREFType(pFileStarImage1->REFNumber) != GetREFType(pFileStarImage2->REFNumber)) &&
            (pFileStarImage1->REFNumber != REF_TYPE_NONE) &&
            (pFileStarImage2->REFNumber != REF_TYPE_NONE)) {
          /* Here we want to replace pFileStarImage1->catalogNumber with pFileStarImage2->catalogNumber */
          for (substitutionIndex = 0; substitutionIndex < substitutionCount; substitutionIndex++) {
            pSubstitution = &substitutionTable[substitutionIndex];
            if ((pSubstitution->badREFNumber == pFileStarImage1->REFNumber) &&
                (pSubstitution->goodREFNumber == pFileStarImage2->REFNumber)) {
              pSubstitution->count++;
              pSubstitution->valid = 1;
              break; /* Already in the table */
            }
          }
          if (substitutionIndex >= substitutionCount) {
            if ((substitutionCount+1) >= allocMagnitudes) {
              printf("ERROR: substitutionTable overrun\n");
              exit(1);
            }
            pSubstitution = &substitutionTable[substitutionCount];
            pSubstitution->badREFNumber = pFileStarImage1->REFNumber;
            pSubstitution->goodREFNumber = pFileStarImage2->REFNumber;
            pSubstitution->valid = 1;
            pSubstitution->count++;
            substitutionCount++;
          }
        }

      }

    }

  }
  qsort((void*)substitutionTable,substitutionCount,sizeof(SUBSTITUTION),SubstitutionCompare);
  /* Now remove duplicates */
  for (substitutionIndex = 0; substitutionIndex < substitutionCount; substitutionIndex++) {
    pSubstitution = &substitutionTable[substitutionIndex];
    if (pSubstitution->valid == 0) {
      continue;
    }
    for (substitutionIndex2 = (substitutionIndex+1); substitutionIndex2 < substitutionCount; substitutionIndex2++) {
      pSubstitution2 = &substitutionTable[substitutionIndex2];
      if (pSubstitution2->valid == 0) {
        continue;
      }
      if ((pSubstitution->badREFNumber == pSubstitution2->badREFNumber) ||
          (pSubstitution->goodREFNumber == pSubstitution2->goodREFNumber)) {
        pSubstitution2->valid = 0;
      }
    }
  }

  /* Now perform the substitutions */
  for (magnitudeIndex0 = 0; magnitudeIndex0 < numMagnitudes; magnitudeIndex0++) {
    pStarImage1 = &pMagnitudeTable[magnitudeIndex0];
    pFileStarImage1 = pStarImage1->pFileStarImage;
    for (substitutionIndex = 0; substitutionIndex < substitutionCount; substitutionIndex++) {
      pSubstitution = &substitutionTable[substitutionIndex];
      if (pSubstitution->badREFNumber == pFileStarImage1->REFNumber)  {
        pFileStarImage1->REFNumber = pSubstitution->goodREFNumber;
        break;
      }
    }
  }
  if (combineAlgorithm == COMBINE_ALGORITHM_BOTH) {
    *pNumMagnitudes = numMagnitudes;
  } else if (combineAlgorithm == COMBINE_ALGORITHM_BEST) {
    curMagnitudes = 0;
    for (magnitudeIndex0 = 0; magnitudeIndex0 < numMagnitudes; magnitudeIndex0++) {
      pStarImage1 = &pMagnitudeTable[magnitudeIndex0];
      pFileStarImage1 = pStarImage1->pFileStarImage;
      if ((magnitudeIndex0+1) < numMagnitudes) {
        pStarImage2 =  &pMagnitudeTable[magnitudeIndex0+1];
        pFileStarImage2 = pStarImage2->pFileStarImage;
        if ((pFileStarImage1->solutionNumber == pFileStarImage2->solutionNumber) &&
            (pFileStarImage1->seriesId == pFileStarImage2->seriesId) &&
            (pFileStarImage1->REFNumber == pFileStarImage2->REFNumber) &&
            (pFileStarImage1->plateNumber == pFileStarImage2->plateNumber) &&
            (pFileStarImage1->NUMBER == pFileStarImage2->NUMBER) &&
            (pFileStarImage1->catalogNumber != CATALOG_APASS) &&
            (pFileStarImage2->catalogNumber == CATALOG_APASS)) {
          /* Here we select the best of both */
          /*
            Astrometry rejection
            *     15      32768       8000 Object rejected because bin has unknown drad.
            *     23    8388608     800000 Object has a drad three times the bin drad or is in a bad spatial bin or local bin

            Photometry rejection
            *     29  536870912   20000000 Object is too bright for accurate magnitudes
            *     13       8192       2000         Star too close to the limiting magnitude
            *     11       2048        800         Maximum isophotonic rms exceeded
            *     28  268435456   10000000 Rejected because of high smoothing correction
            *     12       4096       1000         Maximum locally smoothed rms exceeded

            Irrelevant
            *     24   16777216    1000000 Object is a Pickering Wedge object
            *     25   33554432    2000000 Object fails the defect filter
            *     27  134217728    8000000 Rejected blended object - avoid plotting a limiting magnitude.
            *      6         64         40         Object has high BACKGROUND level
            *      8        256        100         Unmatched object on a multiple-exposure plate
            *      9        512        200         Object has an uncertain time for extinction calculation
            *     14      16384       4000         Star is in bin 9
            *     10       1024        400         Object is a multiple exposure blend
            *     19     524288      80000         Variable or Uncertain Catalog magnitude
            *     20    1048576     100000         Case B - blended
            *     21    2097152     200000         Case C - overlapping Sextractor hits for one GSC star
            *     22    4194304     400000         Case D - mixture of cases B and C
            *     26   67108864    4000000 Copy of the Sextractor blended flag
            *     30 1073741824   40000000 Object is within 23.5 degrees of horizon
            */
          bestChoice = 0;
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)) != 0)) {
            bestChoice = 2;
          }

          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRAD)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_DRAD)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRAD)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRAD)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_DRAD)) != 0)) {
            bestChoice = 2;
          }



          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)) != 0)) {
            bestChoice = 2;
          }


          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) != 0)) {
            bestChoice = 2;
          }


          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)) != 0)) {
            bestChoice = 2;
          }


          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)) != 0)) {
            bestChoice = 2;
          }


          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)) != 0)) {
            bestChoice = 2;
          }


          /* now compare in the following priorities:
             magcal_magdep_rms
             magcal_local_rms
             magcal_iso_rms */
          if ((bestChoice == 0) &&
              (pFileStarImage1->magcal_magdep_rms < 90.0) &&
              (pFileStarImage2->magcal_magdep_rms < 90.0) &&
              (pFileStarImage1->magcal_magdep_rms != pFileStarImage2->magcal_magdep_rms)) {
            if (pFileStarImage1->magcal_magdep_rms < pFileStarImage2->magcal_magdep_rms) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              (pFileStarImage1->magcal_local_rms < 90.0) &&
              (pFileStarImage2->magcal_local_rms < 90.0) &&
              (pFileStarImage1->magcal_local_rms != pFileStarImage2->magcal_local_rms)) {
            if (pFileStarImage1->magcal_local_rms < pFileStarImage2->magcal_local_rms) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              (pFileStarImage1->magcal_iso_rms < 90.0) &&
              (pFileStarImage2->magcal_iso_rms < 90.0) &&
              (pFileStarImage1->magcal_iso_rms != pFileStarImage2->magcal_iso_rms)) {
            if (pFileStarImage1->magcal_iso_rms < pFileStarImage2->magcal_iso_rms) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if (bestChoice == 0) {
            bestChoice = 2;
          }
          if (bestChoice == 1) {
            if (curMagnitudes < magnitudeIndex0) {
              pStarImage0 = &pMagnitudeTable[curMagnitudes];
              pFileStarImage0 = pStarImage0->pFileStarImage;
              memcpy(pStarImage0,pStarImage1,sizeof(PHOTSTARIMAGE));
              memcpy(pFileStarImage0,pStarImage1->pFileStarImage,sizeof(FILESTARIMAGE));
              pStarImage0->pFileStarImage = pFileStarImage0;
            }
            curMagnitudes++;
            magnitudeIndex0++; /* Skip next */
          } else {
            if (curMagnitudes < magnitudeIndex0) {
              pStarImage0 = &pMagnitudeTable[curMagnitudes];
              pFileStarImage0 = pStarImage0->pFileStarImage;
              memcpy(pStarImage0,pStarImage2,sizeof(PHOTSTARIMAGE));
              memcpy(pFileStarImage0,pStarImage2->pFileStarImage,sizeof(FILESTARIMAGE));
              pStarImage0->pFileStarImage = pFileStarImage0;
            }
            curMagnitudes++;
            magnitudeIndex0++; /* Skip next */
          }


        } else {
          /* No match, just copy over */
          if (curMagnitudes < magnitudeIndex0) {
            pStarImage0 = &pMagnitudeTable[curMagnitudes];
            pFileStarImage0 = pStarImage0->pFileStarImage;
            memcpy(pStarImage0,pStarImage1,sizeof(PHOTSTARIMAGE));
            memcpy(pFileStarImage0,pStarImage1->pFileStarImage,sizeof(FILESTARIMAGE));
            pStarImage0->pFileStarImage = pFileStarImage0;
          }
          curMagnitudes++;
        }
      } else {
        /* End of the table, just copy over */
        if (curMagnitudes < magnitudeIndex0) {
          pStarImage0 = &pMagnitudeTable[curMagnitudes];
          pFileStarImage0 = pStarImage0->pFileStarImage;
          memcpy(pStarImage0,pStarImage1,sizeof(PHOTSTARIMAGE));
          memcpy(pFileStarImage0,pStarImage1->pFileStarImage,sizeof(FILESTARIMAGE));
          pStarImage0->pFileStarImage = pFileStarImage0;
        }
        curMagnitudes++;
      }
    }
#if 0
    printf("curMagnitudes %d numMagnitudes %d\n",curMagnitudes,numMagnitudes);
#endif
    numMagnitudes = curMagnitudes;
    *pNumMagnitudes = numMagnitudes;
    return;
  } else if (combineAlgorithm == COMBINE_ALGORITHM_AVERAGE) {
    curMagnitudes = 0;
    for (magnitudeIndex0 = 0; magnitudeIndex0 < numMagnitudes; magnitudeIndex0++) {
      pStarImage1 = &pMagnitudeTable[magnitudeIndex0];
      pFileStarImage1 = pStarImage1->pFileStarImage;


      if ((magnitudeIndex0+1) < numMagnitudes) {
        pStarImage2 =  &pMagnitudeTable[magnitudeIndex0+1];
        pFileStarImage2 = pStarImage2->pFileStarImage;
        if ((pFileStarImage1->solutionNumber == pFileStarImage2->solutionNumber) &&
            (pFileStarImage1->seriesId == pFileStarImage2->seriesId) &&
            (pFileStarImage1->REFNumber == pFileStarImage2->REFNumber) &&
            (pFileStarImage1->plateNumber == pFileStarImage2->plateNumber) &&
            (pFileStarImage1->NUMBER == pFileStarImage2->NUMBER) &&
            (pFileStarImage1->catalogNumber != CATALOG_APASS) &&
            (pFileStarImage2->catalogNumber == CATALOG_APASS)) {
          /* Use the same AFLAGS criteria as in COMBINE_ALGORITHM_BOTH */
#if 0
          if (pFileStarImage1->REFNumber != pFileStarImage2->REFNumber) {
            printf("ERROR line %d REFNumber %lld %lld does not agree\n",__LINE__,pFileStarImage1->REFNumber,pFileStarImage2->REFNumber);
            exit(1);
          }
#endif

          bestChoice = 0;
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)) == 0) {
              bestChoice = 1;
            } else {


              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_DRADBIN)) != 0)) {
            bestChoice = 2;
          }

          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRAD)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_DRAD)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRAD)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_DRAD)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_DRAD)) != 0)) {
            bestChoice = 2;
          }


          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)) == 0) {

              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_TOO_BRIGHT)) != 0)) {
            bestChoice = 2;
          }

          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) != 0)) {
            bestChoice = 2;
          }

          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_ISO_RMS)) != 0)) {
            bestChoice = 2;
          }

          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_HIZOUT)) != 0)) {
            bestChoice = 2;
          }

          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)) != (pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)))) {
            if ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)) == 0) {
              bestChoice = 1;
            } else {
              bestChoice = 2;
            }
          }
          if ((bestChoice == 0) &&
              ((pFileStarImage1->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)) != 0) &&
              ((pFileStarImage2->AFLAGS & (1 << FILTER_AFLAG_LOCAL_RMS)) != 0)) {
            bestChoice = 2;
          }

          if (bestChoice == 0) {
            double mag1;
            double mag2;
            double rms1;
            double rms2;
            bestChoice = 1;
            pFileStarImage1->AFLAGS |= pFileStarImage2->AFLAGS; /* Combine all problems */
            pFileStarImage1->catalogNumber = CATALOG_EXPERIMENTAL;
            mag1 =  pFileStarImage1->magcal_iso;
            rms1 =  sqr(pFileStarImage1->magcal_iso_rms);
            mag2 =  pFileStarImage2->magcal_iso;
            rms2 =  sqr(pFileStarImage2->magcal_iso_rms);
            pFileStarImage1->magcal_iso = ((mag1/rms1)+(mag2/rms2))/((1/rms1)+(1/rms2)); /* Bevington Equation 5-6 p 70. */
            pFileStarImage1->magcal_iso_rms = sqrt(1/((1/rms1)+(1/rms2))); /* Bevington Equation 5-10 p 71. */

            mag1 =  pFileStarImage1->magcal_local;
            rms1 =  sqr(pFileStarImage1->magcal_local_rms);
            mag2 =  pFileStarImage2->magcal_local;
            rms2 =  sqr(pFileStarImage2->magcal_local_rms);
            pFileStarImage1->magcal_local = ((mag1/rms1)+(mag2/rms2))/((1/rms1)+(1/rms2)); /* Bevington Equation 5-6 p 70. */
            pFileStarImage1->magcal_local_rms = sqrt(1/((1/rms1)+(1/rms2))); /* Bevington Equation 5-10 p 71. */

            if ((pFileStarImage1->magcal_magdep_rms < 90.0) &&
                (pFileStarImage2->magcal_magdep_rms < 90.0)) {
              mag1 =  pFileStarImage1->magcal_magdep;
              rms1 =  sqr(pFileStarImage1->magcal_magdep_rms);
              mag2 =  pFileStarImage2->magcal_magdep;
              rms2 =  sqr(pFileStarImage2->magcal_magdep_rms);
              pFileStarImage1->magcal_magdep = ((mag1/rms1)+(mag2/rms2))/((1/rms1)+(1/rms2)); /* Bevington Equation 5-6 p 70. */
              pFileStarImage1->magcal_magdep_rms = sqrt(1/((1/rms1)+(1/rms2))); /* Bevington Equation 5-10 p 71. */
            } else {
              pFileStarImage1->magcal_magdep = pFileStarImage1->magcal_local;
              pFileStarImage1->magcal_magdep_rms = 99.0;
            }
          }
          if (bestChoice == 1) {
            if (curMagnitudes < magnitudeIndex0) {
              pStarImage0 = &pMagnitudeTable[curMagnitudes];
              pFileStarImage0 = pStarImage0->pFileStarImage;
              memcpy(pStarImage0,pStarImage1,sizeof(PHOTSTARIMAGE));
              memcpy(pFileStarImage0,pStarImage1->pFileStarImage,sizeof(FILESTARIMAGE));

              pStarImage0->pFileStarImage = pFileStarImage0;
            }
            curMagnitudes++;
            magnitudeIndex0++; /* Skip next */
          } else {
            if (curMagnitudes < magnitudeIndex0) {
              pStarImage0 = &pMagnitudeTable[curMagnitudes];
              pFileStarImage0 = pStarImage0->pFileStarImage;
              memcpy(pStarImage0,pStarImage2,sizeof(PHOTSTARIMAGE));
              memcpy(pFileStarImage0,pStarImage2->pFileStarImage,sizeof(FILESTARIMAGE));
              pStarImage0->pFileStarImage = pFileStarImage0;
            }
            curMagnitudes++;
            magnitudeIndex0++; /* Skip next */
          }


        } else {
          /* No match, just copy over */
          if (curMagnitudes < magnitudeIndex0) {
            pStarImage0 = &pMagnitudeTable[curMagnitudes];
            pFileStarImage0 = pStarImage0->pFileStarImage;
            memcpy(pStarImage0,pStarImage1,sizeof(PHOTSTARIMAGE));
            memcpy(pFileStarImage0,pStarImage1->pFileStarImage,sizeof(FILESTARIMAGE));
            pStarImage0->pFileStarImage = pFileStarImage0;
          }
          curMagnitudes++;
        }
      } else {
        /* End of the table, just copy over */
        if (curMagnitudes < magnitudeIndex0) {
          pStarImage0 = &pMagnitudeTable[curMagnitudes];
          pFileStarImage0 = pStarImage0->pFileStarImage;
          memcpy(pStarImage0,pStarImage1,sizeof(PHOTSTARIMAGE));
          memcpy(pFileStarImage0,pStarImage1->pFileStarImage,sizeof(FILESTARIMAGE));
          pStarImage0->pFileStarImage = pFileStarImage0;
        }
        curMagnitudes++;
      }
    }
#if 0
    printf("curMagnitudes %d numMagnitudes %d\n",curMagnitudes,numMagnitudes);
#endif
    numMagnitudes = curMagnitudes;
    *pNumMagnitudes = numMagnitudes;
    return;
  } else {
    /* Unsupported is COMBINE_ALGORITHM_AVERAGE */
    printf("ERROR: unsupported combineAlgorithm %d\n",combineAlgorithm);
    exit(1);
  }
  return;
}


int
ProcessLongTransientCandidates(
  PFILECOMMON pFileCommon,
  PTRANCOMMON pTranCommon,
  PPARAMETERVECTORSTORE pVectorStore,
  PPHOTSTARIMAGE pMagnitudeTable,
  PGALAXYCOMMON pGalaxyCommon,
  PSTARENTRY pCurStarEntry,
  char *nearbyObjects,
  int AFLAGSMASK1,
  int QUALITYMASK
) {
  int index_ngood;
  PMALMQUIST pMalmquist1 = NULL;                    /* Preceeding star image */
  PMALMQUIST pMalmquist1X = NULL;
  PMALMQUIST pMalmquist2 = NULL;              /* Middle image of interest */
  PPHOTSTARIMAGE pSumStarImage2 = NULL;
  PFILESTARIMAGE pFileStarImage2 = NULL;
  PMALMQUIST pMalmquist3 = NULL;              /* Following star image */
  PMALMQUIST pMalmquist3X = NULL;
  PPHOTSTARIMAGE pSumStarImage3 = NULL;
  PFILESTARIMAGE pFileStarImage3 = NULL;
  PMALMQUIST pMalmquist4;                    /* temporary placeholder */
  PPHOTSTARIMAGE pSumStarImage4 = NULL;
  PFILESTARIMAGE pFileStarImage4 = NULL;
  PMALMQUIST pMalmquist5 = NULL;              /* temporary file */
  PPHOTSTARIMAGE pSumStarImage5 = NULL;
  PFILESTARIMAGE pFileStarImage5 = NULL;
  int peak_index = 0;
  int temp_index;
  int rejectFlag = 1;
  double min_magcal_magdep;
  double min_nodefect_magcal_magdep;
  double max_magcal_magdep;
  double ave_magcal_magdep = 99.0;
  double ave_magcal_magdep2;
  double count_magcal_magdep2;
  char REF[MAX_REF];
  double dradmed = 99.0;
  double dradrms = 99.0;
  double dradRMS2med = 99.0;
  double dradRMS2rms = 99.0;
  double dradRMS3med = 99.0;
  double dradRMS3rms = 99.0;
  double dradRMS3;
  double debugPrintFlag = 0;
  int dradIndex;
  int nearbyCatalogIndex;
  PNEARBYCATALOG pNearbyCatalog;
  double factor;
  double curDistance;
  int refType;
  double nearby_magcal_magdep = 99.;
  double patrolScale;
  double raAverage = 0;
  double decAverage = 0;
  double positionWeight = 0;
  int positionCount = 0;
  double dateyear;
  int intyear;
  PPLATELIMITINGREC pPlateLimitingRec;
  int plateIndex;
  int plateLimitingYears[PIPELINE_MAX_DATE-PIPELINE_MIN_DATE+1];

  if (pTranCommon->full_ngood < (2*LONG_TC_POINTS)) {
    pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_LONG_INSUFFICIENT_POINTS);
    return(0); /* not enough stars to be of interest */
  }


  /* Resort the MAMQUIST selections according to date */
  qsort((void*)pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST],pTranCommon->full_ngood,sizeof(MALMQUIST),LightcurveDateCompare);

  int index_first = 0; /* First year  */
  int count_first = 0;
  double date_first = 0;

  double sum_first = 0;
  int index_second = 0; /* Second year */
  int count_second = 0;
  int index_third = 0;
  double date_second = 0;
  double sum_second = 0;
  double date_third = 0; /* Current year */
  double date_skip = -LONG_TC_SKIP_DAYS;

  for (index_third = 0; index_third < pTranCommon->full_ngood; index_third++) {
    pMalmquist3 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_third];
    pSumStarImage3 = &pMagnitudeTable[pMalmquist3->imageIndex];
    pFileStarImage3 = pSumStarImage3->pFileStarImage;
    if (pFileStarImage3->magcal_magdep > 90) {
      printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
      exit(1);
    }
    if ((pMalmquist3->selectionFlags & SELECT_GOOD) == 0) {
      printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
      exit(1);
    }
    if (pMalmquist3->Date != pFileStarImage3->Date) {
      printf("ERROR: ProcessTransientCandidates selection error in line %d\n",__LINE__);
      exit(1);
    }


    pMalmquist2 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_second];
    pSumStarImage2 = &pMagnitudeTable[pMalmquist2->imageIndex];
    pFileStarImage2 = pSumStarImage2->pFileStarImage;
    pMalmquist1 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_first];

    count_second++;
    sum_second += pMalmquist3->magcal_magdep;
    date_third = pMalmquist3->Date;

    if (index_third == 0) {
      index_first = 0;
      date_first = date_third;
      index_second = 0;
      date_second = date_third;
    } else {
      while ((date_third - date_second) > LONG_TC_FLARE_DAYS) {
        /* Reduce the count by one */
        sum_second -= pMalmquist2->magcal_magdep;
        sum_first += pMalmquist2->magcal_magdep;
        count_second--;
        count_first++;
        /* Advance the pointer */
        index_second++;
        pMalmquist2 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_second];
        pSumStarImage2 = &pMagnitudeTable[pMalmquist2->imageIndex];
        pFileStarImage2 = pSumStarImage2->pFileStarImage;
        date_second = pMalmquist2->Date;
      }

      while ((date_second - date_first) > LONG_TC_PRE_DAYS) {
        /* Reduce the count by one */
        sum_first -= pMalmquist1->magcal_magdep;
        count_first--;
        /* Advance the pointer */
        index_first++;
        pMalmquist1 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_first];
        date_first = pMalmquist1->Date;
      }
    }

    if ((count_first >= LONG_TC_POINTS) && (count_second >= LONG_TC_POINTS) && ((date_third-date_first) > LONG_TC_MIN_DAYS) && (date_third > date_skip)) {
      if (((sum_first/(1.0*count_first))-(sum_second/(1.0*count_second))) >  LONG_TC_MAGNITUDE) {
        if (debugPrintFlag == 1) {
          printf("EVENT!! %2d year: %7.2f count_first %d average_first %.2f count_second %d average_second %.2f\n",pTranCommon->peakLongOutburst+1,jd2ep(pMalmquist2->Date),count_first,(sum_first/(1.0*count_first)),count_second,(sum_second/(1.0*count_second)));
        }
        pTranCommon->peakLongOutburst++;
        if (pTranCommon->peakLongOutburst == 1) {
          /* Fill in information about this outburst */
          min_magcal_magdep = (sum_first/(1.0*count_first));
          max_magcal_magdep = (sum_second/(1.0*count_second));
          rejectFlag = 0;
          pTranCommon->peakYear = jd2ep(pMalmquist2->Date);
          pTranCommon->peakDays = date_third-date_second;
          pTranCommon->peakMag = max_magcal_magdep - min_magcal_magdep;
          pTranCommon->peakCount = count_second;
          pTranCommon->peakOutside = pTranCommon->full_ngood-pTranCommon->peakCount;
          pTranCommon->peakNoDefectMag = 99.0;
          pTranCommon->peakMaxMag = min_magcal_magdep;
          pTranCommon->peakMinMag = max_magcal_magdep;
          pTranCommon->peakNumber = peak_index+1;
          pTranCommon->peakEventCount = 1;
          pMalmquist1X = pMalmquist1;
          pMalmquist3X = pMalmquist3;
          dradIndex = 0;
          ave_magcal_magdep2 = 0;
          count_magcal_magdep2 = 0;
          pTranCommon->peakMultipleCount = 0;
          for (temp_index = index_second; temp_index <= index_third; temp_index++) {
            pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
            pSumStarImage4 = &pMagnitudeTable[pMalmquist4->imageIndex];
            if ((pSumStarImage4->quality & QUALITY_MULTIPLE) != 0) {
              pTranCommon->peakMultipleCount++; /* Count multiple exposure plates within the lightcurve */
            }
            ave_magcal_magdep2 += pMalmquist4->magcal_magdep;
            count_magcal_magdep2++;
            if (temp_index == index_second) {
              if ((pMalmquist4->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) {
                min_nodefect_magcal_magdep = pMalmquist4->magcal_magdep;
              } else {
                min_nodefect_magcal_magdep = 99.9;
              }
              min_magcal_magdep = pMalmquist4->magcal_magdep;
              max_magcal_magdep = pMalmquist4->magcal_magdep;
            } else {
              if (min_magcal_magdep > pMalmquist4->magcal_magdep) {
                min_magcal_magdep = pMalmquist4->magcal_magdep;
              }
              if ((pMalmquist4->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) {
                if (min_nodefect_magcal_magdep > pMalmquist4->magcal_magdep) {
                  min_nodefect_magcal_magdep = pMalmquist4->magcal_magdep;
                }
              }
              if (max_magcal_magdep < pMalmquist4->magcal_magdep) {
                max_magcal_magdep = pMalmquist4->magcal_magdep;
              }
            }


#if defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            pVectorStore->vector[PARAMETERVECTOR_DRAD][dradIndex] = pMalmquist4->drad;
            pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = pMalmquist4->dradRMS2;
            dradIndex++;
#endif
#if defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if ((pMalmquist4->AFLAGS & (1 << FILTER_AFLAG_BACKGROUND)) == 0) {
              pVectorStore->vector[PARAMETERVECTOR_DRAD][dradIndex] = pMalmquist4->drad;
              pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = pMalmquist4->dradRMS2;
              dradIndex++;
            }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if (pMalmquist4->AFLAGS & (1 << FILTER_AFLAG_DEFECT) == 0) { /* change of Jul 13, 2015*/
              pVectorStore->vector[PARAMETERVECTOR_DRAD][dradIndex] = pMalmquist4->drad;
              pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = pMalmquist4->dradRMS2;
              dradIndex++;
            }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
            if ((pMalmquist4->AFLAGS & ((1 << FILTER_AFLAG_DEFECT) | (1 << FILTER_AFLAG_BACKGROUND)) == 0)) { /* change of Jul 13, 2015*/
              pVectorStore->vector[PARAMETERVECTOR_DRAD][dradIndex] = pMalmquist4->drad;
              pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = pMalmquist4->dradRMS2;
              dradIndex++;
            }
#endif
          }
          if (CalcMedianAndRMS(dradIndex,pTranCommon->minGoodStars,pVectorStore->vector[PARAMETERVECTOR_DRAD],&dradmed,&dradrms,0,3.0,1) == 0) {
            dradmed = 99.0;
            dradrms = 99.0;
          }
          if (CalcMedianAndRMS(dradIndex,pTranCommon->minGoodStars,pVectorStore->vector[PARAMETERVECTOR_DRADRMS2],&dradRMS2med,&dradRMS2rms,0,3.0,1) == 0) {
            dradRMS2med = 99.0;
            dradRMS2rms = 99.0;
          }
          pTranCommon->peakDradRMS = dradrms;
          pTranCommon->peakDradRMS2 = dradRMS2rms;
          /* Count wide and narrow field points for the lightcurve and find the transient weighted position */
          for (index_ngood = 0; index_ngood < pTranCommon->full_ngood; index_ngood++) {
            pMalmquist5 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][index_ngood];
            pSumStarImage5 = &pMagnitudeTable[pMalmquist5->imageIndex];
            pFileStarImage5 = pSumStarImage5->pFileStarImage;
            patrolScale = GetFittedPlateScale(pFileStarImage5->seriesId,pFileStarImage5->plateNumber);
#if 0
            qualityMask = GetPlateQualityMask(pFileCommon,pFileStarImage5->seriesId,pFileStarImage5->plateNumber,&quality,&pColorterm,&pColorflag,pFileStarImage5->REFNumber);
            if ((quality == 0) &&
                ((qualityMask & (1 << pFileStarImage5->spatial_bin)) != 0)) {
              printf("line %d Bad colorterm for %s plate %s spatial bin %d gsc_bin_index %d\n",__LINE__,pCurStarEntry->REF,pSumStarImage5->Plate,pFileStarImage5->spatial_bin,pFileStarImage5->gsc_bin_index);
            }
#endif
            if ((index_ngood < index_second) || (index_ngood > index_third)) {
              if (patrolScale > PATROL_PLATE_SCALE) {
                pTranCommon->peakOutsideWF++;
              } else {
                pTranCommon->peakOutsideNF++;
              }
            } else {
#if defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
              if ((pFileStarImage5->dradRMS2 > 0) &&
                  (pFileStarImage5->dradRMS2 < 90.0)) {
                raAverage += (pFileStarImage5->ra+360.)/sqr(pFileStarImage5->dradRMS2);
                decAverage += (pFileStarImage5->dec+180.)/sqr(pFileStarImage5->dradRMS2);
                positionWeight += 1.0/sqr(pFileStarImage5->dradRMS2);
                positionCount++;
              }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
              if (((pFileStarImage5->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0) &&
                  (pFileStarImage5->dradRMS2 > 0) &&
                  (pFileStarImage5->dradRMS2 < 90.0)) {
                raAverage += (pFileStarImage5->ra+360.)/sqr(pFileStarImage5->dradRMS2);
                decAverage += (pFileStarImage5->dec+180.)/sqr(pFileStarImage5->dradRMS2);
                positionWeight += 1.0/sqr(pFileStarImage5->dradRMS2);
                positionCount++;
              }
#endif
#if defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
              if (((pFileStarImage5->AFLAGS & (1 << FILTER_AFLAG_BACKGROUND)) == 0) &&
                  (pFileStarImage5->dradRMS2 > 0) &&
                  (pFileStarImage5->dradRMS2 < 90.0)) {
                raAverage += (pFileStarImage5->ra+360.)/sqr(pFileStarImage5->dradRMS2);
                decAverage += (pFileStarImage5->dec+180.)/sqr(pFileStarImage5->dradRMS2);
                positionWeight += 1.0/sqr(pFileStarImage5->dradRMS2);
                positionCount++;
              }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
              if (((pFileStarImage5->AFLAGS & ((1 << FILTER_AFLAG_DEFECT) | (1 << FILTER_AFLAG_BACKGROUND))) == 0) &&
                  (pFileStarImage5->dradRMS2 > 0) &&
                  (pFileStarImage5->dradRMS2 < 90.0)) {
                raAverage += (pFileStarImage5->ra+360.)/sqr(pFileStarImage5->dradRMS2);
                decAverage += (pFileStarImage5->dec+180.)/sqr(pFileStarImage5->dradRMS2);
                positionWeight += 1.0/sqr(pFileStarImage5->dradRMS2);
                positionCount++;
              }
#endif

              if ((pFileStarImage5->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) != 0) {
                pTranCommon->peakDefectCount++; /* Count defects within the lightcurve */
              }
              if (patrolScale > PATROL_PLATE_SCALE) {
                pTranCommon->peakCountWF++;
              } else {
                pTranCommon->peakCountNF++;
              }
              /* Here we count points with uncertain or bad colorterm correction */
              refType = GetREFType(pCurStarEntry->REFNumber);
              if ((refType != REF_TYPE_DASCH) &&
                  (pCurStarEntry->Stdmag < 90.0)) {
                /* This is a catalog star.  Flag bad color if the color is bad or the colorterm is bad */
                if ((pCurStarEntry->color > 90.0) ||
                    (pSumStarImage5->colorflag != COLORFLAG_METROPOLIS) ||
                    (pSumStarImage5->colorterm > 90.0)) {
#if 0
                  printf("ERROR line %d badcolorflag for REF %s gsc_bin_index %d and plate %s\n",__LINE__,pSumStarImage5->REF,pFileStarImage5->gsc_bin_index,pSumStarImage5->Plate);
#endif
                  pSumStarImage5->badcolorflag = 1;
                  pTranCommon->peakBadColorCount++;
                }
              } else {
                /* This is an unmatched star or matched to the wrong catalog.  Flag bad color if the colorterm is bad or if
                   the colorterm is not blue */
                if ((pSumStarImage5->colorflag != COLORFLAG_METROPOLIS) ||
                    (pSumStarImage5->colorterm > 90.0) ||
                    (BlueColorterm(pFileCommon->catalogNumber,pSumStarImage5->colorterm) == 0)) {
#if 0
                  printf("ERROR line %d badcolorflag for REF %s gsc_bin_index %d and plate %s\n",__LINE__,pSumStarImage5->REF,pFileStarImage5->gsc_bin_index,pSumStarImage5->Plate);
#endif
                  pSumStarImage5->badcolorflag = 1;
                  pTranCommon->peakBadColorCount++;
                }
              }
            }
          }

          if (positionCount > 1) {
            raAverage = raAverage/positionWeight;
            while (raAverage >= 360.0) {
              raAverage -= 360.0;
            }
            decAverage = decAverage/positionWeight;
            while (decAverage > 90.0) {
              decAverage -= 180.0;
            }
            pTranCommon->peakRA = raAverage;
            pTranCommon->peakDec = decAverage;


            dradIndex = 0;
            for (temp_index = index_second; temp_index <= index_third; temp_index++) {
              pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
              pSumStarImage4 = &pMagnitudeTable[pMalmquist4->imageIndex];
              pFileStarImage4 = pSumStarImage4->pFileStarImage;
              factor =  cos(DEGREES_TO_RAD*(pTranCommon->peakDec));
              dradRMS3 = sqrt(sqr(3600.0*(pFileStarImage4->dec-pTranCommon->peakDec))+sqr(3600*factor*(pFileStarImage4->ra-pTranCommon->peakRA)));
#if defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
              pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = dradRMS3;
              dradIndex++;
#endif
#if defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
              if ((pMalmquist4->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) == 0) {
                pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = dradRMS3;
                dradIndex++;
              }
#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
              if ((pMalmquist4->AFLAGS & (1 << FILTER_AFLAG_BACKGROUND)) == 0) {
                pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = dradRMS3;
                dradIndex++;
              }

#endif
#if !defined(INCLUDE_TRANSIENT_DEFECTS) && !defined(INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND)
              if (pMalmquist4->AFLAGS & ((1 << FILTER_AFLAG_DEFECT) | (1 << FILTER_AFLAG_BACKGROUND)) == 0) {
                pVectorStore->vector[PARAMETERVECTOR_DRADRMS2][dradIndex] = dradRMS3;
                dradIndex++;
              }
#endif
            }
            /* Now compute peakDradRMS3 with respect to this transient center */
            if (CalcMedianAndRMS(dradIndex,pTranCommon->minGoodStars,pVectorStore->vector[PARAMETERVECTOR_DRADRMS2],&dradRMS3med,&dradRMS3rms,0,3.0,1) == 0) {
              dradRMS3med = 99.0;
              dradRMS3rms = 99.0;
            }
          }
          pTranCommon->peakDradRMS3 = dradRMS3rms;
          /* Count how many points are outside the flare */
          pTranCommon->peakExtra = 0;
          for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
            pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
            if ((pMalmquist4->Date < pMalmquist1->Date) ||
                (pMalmquist4->Date > pMalmquist3->Date)) {
              pTranCommon->peakExtra++;
            }
          }
          ave_magcal_magdep =  (max_magcal_magdep+min_magcal_magdep)/2.0;
          GetREF(pFileStarImage2->REFNumber,REF,0,0);
          pTranCommon->peakRMS = 99;
#ifdef PLOT_LIMITING_MAGNITUDES
          LoadGalaxyTable(pGalaxyCommon,pTranCommon->peakRA,pTranCommon->peakDec);
          memset(plateLimitingYears,0,sizeof(plateLimitingYears));

          /* Convert the Julian Date to a year, and eliminate all points with known magnitudes */
          for (plateIndex = 0; plateIndex < pGalaxyCommon->plateCount; plateIndex++) {
            pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
            /* Skip plates that don't go as deep as the average brightness of the flare */
            if (pPlateLimitingRec->limiting_mag_local < ((pTranCommon->peakMaxMag+pTranCommon->peakMinMag)/2.0)) {
              continue;
            }
            pTranCommon->peakLimitingPoints++;
            dateyear = jd2ep(pPlateLimitingRec->geoJulianDate);
            intyear = dateyear-PIPELINE_MIN_DATE;
            if (intyear > (PIPELINE_MAX_DATE-PIPELINE_MIN_DATE)) {
              printf("ERROR: intyear %d is too large\n",intyear+PIPELINE_MIN_DATE);
              intyear = PIPELINE_MAX_DATE-PIPELINE_MIN_DATE;
            }
            if (plateLimitingYears[intyear] == 0) {
              pTranCommon->peakLimitingYears++;
            }
            plateLimitingYears[intyear]++;
          }
#endif /* PLOT_LIMITING_MAGNITUDES */
        }
#if 1
        /* Change of May 2, 2018: continue with averaging, but do not
           report anything until LONG_TC_SKIP_DAYS after */
        date_skip = date_third+LONG_TC_SKIP_DAYS;
#else
        /* Now reset all of the counts and pointers */
        index_first = index_third;
        index_second = index_third;
        date_third = date_third;
        date_first = date_third;
        date_second = date_third;
        count_first = 0;
        count_second = 0;
        pMalmquist2 = pMalmquist3;
        pSumStarImage2 = pSumStarImage3;
        pFileStarImage2 = pFileStarImage3;
        pMalmquist1 = pMalmquist3;
        pSumStarImage1 = pSumStarImage3;
        pFileStarImage1 = pFileStarImage3;
#endif

      }
    }
  }
  if (rejectFlag == 0) {
    if (pTranCommon->peakLongOutburst <= LONG_TC_MAX_EVENTS) {
      /* Search the nearby_catalog_table for anything within 90 arcsec */
      factor = cos(DEGREES_TO_RAD*pCurStarEntry->dec);
      for (nearbyCatalogIndex = 0; nearbyCatalogIndex < pFileCommon->nearbyCatalogCount; nearbyCatalogIndex++) {
        pNearbyCatalog = &pFileCommon->nearby_catalog_table[nearbyCatalogIndex];
        if (pNearbyCatalog->REFNumber == pCurStarEntry->REFNumber) {
          continue;
        }
        if ((pNearbyCatalog->dec > 90.0) ||
            (pNearbyCatalog->ra > 990.0)) {
          continue;
        }
        nearby_magcal_magdep = pNearbyCatalog->magcal_magdepGoodMag;
        if (nearby_magcal_magdep > 90) {
          nearby_magcal_magdep = pNearbyCatalog->magcal_magdep;
        }
        if (nearby_magcal_magdep > 90.0) {
          continue;
        }
        curDistance = 3600.*sqrt(sqr(pNearbyCatalog->dec-pCurStarEntry->dec) + sqr(factor*(pNearbyCatalog->ra-pCurStarEntry->ra)));
        if (fabs(nearby_magcal_magdep-min_magcal_magdep) < MAX_TC_NEARBY_MAG) {
          if (pTranCommon->peakNearbyDistance == 0) {
            pTranCommon->peakNearbyDistance = curDistance;
          } else if (curDistance <  pTranCommon->peakNearbyDistance) {
            pTranCommon->peakNearbyDistance = curDistance;
          }
        } else if (fabs(nearby_magcal_magdep-ave_magcal_magdep) < MAX_TC_NEARBY_MAG) {
          if (pTranCommon->peakNearbyDistance == 0) {
            pTranCommon->peakNearbyDistance = curDistance;
          } else if (curDistance <  pTranCommon->peakNearbyDistance) {
            pTranCommon->peakNearbyDistance = curDistance;
          }
        }
      }

      PlotTransientCandidates(pFileCommon,pTranCommon,pVectorStore,pMagnitudeTable,pGalaxyCommon,pCurStarEntry,refType,nearbyObjects,AFLAGSMASK1,QUALITYMASK,pMalmquist1X,pMalmquist3X);
      pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_LONG_ACCEPTED);

      return(1);
    } else {
      /* Too many events!!! flush everything except peakLongOutburst*/
      pTranCommon->peakBadColorCount = 0;
      pTranCommon->peakCount = 0;
      pTranCommon->peakCountNF = 0;
      pTranCommon->peakCountWF = 0;
      pTranCommon->peakDays = 0;
      pTranCommon->peakDec = 99.;
      pTranCommon->peakDefectCount = 0;
      pTranCommon->peakDradRMS = 0;
      pTranCommon->peakDradRMS2 = 0;
      pTranCommon->peakDradRMS3 = 0;
      pTranCommon->peakEventCount = 0;
      pTranCommon->peakExtra = 0;
      pTranCommon->peakLimitingPoints = 0;
      pTranCommon->peakLimitingYears = 0;
      pTranCommon->peakMag = 0;
      pTranCommon->peakMaxMag = 0;
      pTranCommon->peakMinMag = 0;
      pTranCommon->peakMultipleCount = 0;
      pTranCommon->peakNoDefectMag = 0;
      pTranCommon->peakNumber = 0;
      pTranCommon->peakOutside = 0;
      pTranCommon->peakOutsideNF = 0;
      pTranCommon->peakOutsideWF = 0;
      pTranCommon->peakRA = 999.;
      pTranCommon->peakRMS = 0;
      pTranCommon->peakYear = 0;
      pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_LONG_TOO_MANY_EVENTS);
      return(0);
    }
  } else {
    pTranCommon->peakEvaluation |= (1 << PEAKEVALUATION_LONG_REJECTED);
    return(0);
  }
}


/* Mar  3, 2009 Edward J. Los - Initial Version
 * May  6, 2009 Edward J. Los - add mosaicNumber to photplates
 * May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 * Jun 23, 2009 Edward J. Los - Add plate edges and quality field to photplates.
 * Jun 26, 2009 Edward J. Los - Add passBits to the magnitudes and starsummary tables
 * Jun 29, 2009 Edward J. Los - Add InsertNewMagnitude
 *                            - Add PurgeMagnitudes to exclude stale entries
 * Jul 19, 2009 Edward J. Los - Reject FILTER_AFLAG_BLEND_NOMATCH stars in GetMagnitudes
 * Aug 20, 2009 Edward J. Los - Add Kepler Input Catalog support
 * Aug 31, 2009 Edward J. Los - Move photometry conversion routines to pipelineutils.h
 * Sep 29, 2009 Edward J. Los - Add GetFittedPlateScale
 * Nov 19, 2009 Edward J. Los - Change STARIMAGE to PHOTSTARIMAGE
 * Jan 27, 2010 Edward J. Los - Support rejectFlag in the file-based magnitude table
 * Feb 12, 2010 Edward J. Los - Make photometry files group writable
 * Feb 23, 2010 Edward J. Los - Move ProcessNoneImages here from search_none
 *                              Include printout of spatial bin and sextractor flags
 *                              Search radius is the sum of both dradRMS2 added in
 *                               quadrature times 2.5
 * Mar  2, 2010 Edward J. Los - Filter unmatched objects near bright GSC2.3.2 stars
 * Mar 17, 2010 Edward J. Los - Add subroutineCalls, readMagnitudes, currentMagnitudes, and returnedMagnitudes counters
 * Apr  2, 2010 Edward J. Los - Support kepler magnitude files
 * Apr 16, 2010 Edward J. Los - Add excludeSeriesCount
 * May  4, 2010 Edward J. Los - Add file semaphore support
 * Sep  6, 2010 Edward J. Los - Correct Solaris compiler warnings
 *                              Add READONLY_PHOTOMETRY to avoid modifying photometry database files
 *                              Add SOLARIS_BUILD to swap bytes of the database files
 * Nov  1, 2010 Edward J. Los - Correct long long references in a print statement.
 * Nov  2, 2010 Edward J. Los - Add photometry file repair capability
 * Nov  5, 2010 Edward J. Los - Change group ownership of files to the scanner group
 * Dec 18, 2010 Edward J. Los - Skip change group error check - it will fail if the owner is scanner and the
 *                              user is not.
 * Dec 31, 2010 Edward J. Los - Added min versionID and filteredVersionCount
 * Feb 14, 2011 Edward J. Los - Correct GetMagnitudeFileName for Kepler
 * Mar 23, 2011 Edward J. Los - Add apass catalog support
 * Apr  8, 2011 Edward J. Los - Add mixedCatalogCount
 * Apr 26, 2011 Edward J. Los - Add series to the photometry database
 * Jun  1, 2011 Edward J. Los - Ensure that the value of umask is 2 for group rw access
 * Jun  3, 2011 Edward J. Los - Make all web reads readonly
 * Jun 13, 2011 Edward J. Los - Correct bug in GetPlateScale which used "exposureNumber" instead of "solutionNumber"
 * Aug 31, 2011 Edward J. Los - Temporarily remove the REFNumber from the header sanity check (TEMP_REFNUMBER_FIX)
 * Nov  9, 2011 Edward J. Los - Add magnitude-dependent calibration
 * Nov 14, 2011 Edward J. Los - Make minVersionId dependent on the catalog
 * Nov 22, 2011 Edward J. Los - Correct starbase second line in WriteStarbaseRecord
 * Nov 29, 2011 Edward J. Los - Sequester non-Harvard series.
 * Dec 16, 2011 Edward J. Los - Test new limiting magnitude criteria of Sumin's memo of Thu 12/15/11 1:52 PM
 * Feb 21, 2012 Edward J. Los - When reading magnitudes in read-only mode, return what would have been written.
 * Jul 18, 2012 Edward J. Los - Make sure that all calls to GetRef do not return spaces in object names
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support.
 * Aug 27, 2012 Edward J. Los - Add SEARCH_NONE_MIN_DRAD
 * Sep 10, 2012 Edward J. Los - Reorganize code to optimize magnitude file support.
 * Sep 18, 2012 Edward J. Los - Expand search for unmatched objects beyong gsc bin borders.
 * Oct  8, 2012 Edward J. Los - Improve SearchGSCBlends performance and accuracy by using a k-d tree
 * Oct 12, 2012 Edward J. Los - Flag an error for GetStarEntry gsc_bin_index searches
 * Oct 22, 2012 Edward J. Los - ProcessNoneImages will create DASCH images only to good images, but will
 *                              assign a DASCH number to an image which is bad, but would have been included
 *                              in the group
 * Oct 30, 2012 Edward J. Los - Implement SINGLE_BIN_CACHE to reduce memory consumption of the magnitude cache
 *                              at the cost of more disk reads.
 *                              Expand ProcessNoneImages search by aLength, the highest ellipse semi-radius.
 * Nov  8, 2012 Edward J. Los - Add AFLAGS, BFLAGS (replacing SFLAGS), lat, lon, ELLIPTICITY, THETA_J2000, and limiting_mag_local to the search_none.db file
 *                              Allow all AFLAGS but FILTER_AMASK_NONE_CANDIDATE in search_none.db and quality except for multiple, grating, and spectra.
 *                              Correct colorterm code in GetPlateQualityMask
 * Nov 20, 2012 Edward J. Los - Eliminate duplicate plate reference when searching for unmatched objects.
 * Dec  1, 2012 Edward J. Los - Add a star to a group if the positions decode to the same REFNumber
 * Dec  7, 2012 Edward J. Los - Correct a bug which rejected GSC objects in the APASS catalog.
 * Dec 10, 2012 Edward J. Los - Add PERMISSIVE_ID_TABLE to include more points in the id and summary tables
 *                              Add ra and dec to WriteStarbaseRecord, short version.  Add a title to WriteStarbaseRecord headers
 * Jan 15, 2012 Edward J. Los - Correct inclusion of two stars from the same plate for second-quality unmatched stars.
 * Jan 25, 2012 Edward J. Los - Limit unmatched object search to BIN_EXPANSION_RADIUS
 * Feb 15, 2013 Edward J. Los - Enforce the new MAX_LIMITING_MAG limit when reading data files
 *                              In ProcessNoneImagesX, if the nearest catalog object does not have a good position, use it anyway with a "?".
 * Feb 22, 2013 Edward J. Los - In ProcessNoneImagesX, allow all catalog objects for the first nearby object
 *                                                     allow only objects 0.5 mag above the limiting magnitude and +/- 1 mag different from host star.
 * Mar 11, 2013 Edward J. Los - Add RaPM, DecPM, ra_2 and dec_2 to the magnitude file (Version 5);
 * Mar 25, 2013 Edward J. Los - Use a different semaphores for each catalog database.
 * Jul  8, 2013 Edward J. Los - Modify LoadNoneImages to be sure that unmatched images affected by a bug found in filterblended.c on July 3, 2013 (CHECKBLEND_BUGFIX)
 * Jul 18, 2013 Edward J. Los - For proper motion studies, do not consider saturated and trailed images
 *                              When calculating aLength and bLength use (ISO3+ISO4)/2 to estimate the 50% flux size.
 * Aug 12, 2013 Edward J. Los   Add proper motion standard deviations to GSCIMAGE
 *                              Reformat GSCIMAGE, eliminating GSCIMAGEX
 * Aug 27, 2013 Edward J. Los   Correct problem with stale pNoneStats structures.
 * Sep  3, 2013 Edward J. Los   Collect statistics on the range on Transient Candidates in RA and declination.
 * Sep 10, 2013 Edward J. Los   Correct the bug in the FILTER_BFLAG_LATEMATCH which allowed duplicate plates.
 * Sep 17, 2013 Edward J. Los - Add SEARCH_NONE_MIN_ARCSEC
 * Oct 28, 2013 Edward J. Los - Add BLEND_TC_SEARCH for the new search for transient candidates.
 * Nov  8, 2013 Edward J. Los - Make sure that no two blend groups have the same REFNumber by adding 1 arcsec declination to one of them.
 *                              Eliminate all group members that are greater than SEARCH_NONE_MAX_ARCSEC from the root.
 * Dec 10, 2013 Edward J. Los - Get statistics of AFLAGS settings for root Transient candidate objects
 * Dec 24, 2013 Edward J. Los - Avoid making Plate Defects and SExtractor blends the root Transient Candidate object.
 * Jan 22, 2014 Edward J. Los - Add a website timeout
 * Mar 31, 2014 Edward J. Los - Correct new proper motion algorithm for short span cases.
 * Apr  2, 2014 Edward J. Los - Implement a series hash table
 * Apr 10, 2014 Edward J. Los - Record limiting magnitudes only from good local bins.
 * Nov 10, 2014 Edward J. Los - Replace FILTER_AMASK_LOWDRAD with FILTER_AMASK_NONE_CANDIDATE2 which adds the FILTER_AFLAG_LIMITING_MAG to rejected transient candidates.
 * Dec  2, 2014 Edward J. Los - Add temporary code to filter out GSC and Tycho objects from LoadNoneImages (USE_ONLY_APASS defined)
 * Dec 23, 2014 Edward J. Los - V6 data format: add A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex
 *                              Remove SOLARIS_BUILD support
 * Jan 27, 2015 Edward J. Los - Add GetTimeAccuracy
 * Mar 24, 2015 Edward J. Los - Add REJECT_MULTIPLE_FOR_ID to reject multiple exposure plates for the ID tables
 * Apr 13, 2015 Edward J. Los - Add extraColumnFlag to support the galactic latitude column
 * May  4, 2015 Edward J. Los - Move above to changes into filecommon flags: extraColumnFlag and enableTransientSearch flag
 *                              Accept FILTER_AFLAG_TOO_BRIGHT when enableTransientSearch flag is set
 * May  6, 2015 Edward J. Los - Add nearbyREFflag if a TC has a catalog star within FWHM/2 and the TC is dimmer than the catalog object or the TC is within 0.5 mag of the catalog object (modified Josh request of Tue 5/05/15 9:15 AM)
 * Jun 15, 2015 Edward J. Los - Reject transient candidates with a dateRange > 100 days.  Reduce MIN_TC_CLIP_NGOOD to 10.  Increase MAX_TC_DAYS to 100
 *                              Change peakSlopeRMS to peakRMS, add peakUpperCount, peakNumber, peakMaxMag, and peakMinMag
 * Jun 19, 2015 Edward J. Los - Increase MIN_TC_CLIP_NGOOD to 100
 *                              Turn off the requirment for MAX_PEAK_RATIO
 * Jun 20, 2015 Edward J. Los - Correct divide by zero error for variancecount == 0
 * Jun 23, 2015 Edward J. Los - Add Stdmag and color for the the transient search
 * Jun 24, 2015 Edward J. Los - Set nearbyREFflag  if the TC has a non-SDSS catalog star within min(3 arcsec,FWHM/2) and brighter than 16th magnitude.  Add peakEventCount
 * Jun 25, 2015 Edward J. Los - Remove MIN_TC_DRAD astrometry check for potential asteroids
 * Jul 13, 2015 Edward J. Los - Remove BLEND_TC_SEARCH conditional (Now the default).  Do not allow FILTER_AFLAG_DEFECT to become a root object
 *                              Redefine peakDays as the full range of the transient (dateRange).  The new peakDays must be greater then MIN_TC_DAYS, instead of the time after the peak.
 *                              Add peakEventCount2, the number of valid peaks with at least one point
 *                              Recalculate a transient position as being the dradRMS2 weighted mean of non-defects in the group.
 *                              Redefine nearbyREFflag to be set at min(18 arcsec,(1.5*FWHM)+dradrms)
 * Jul 16, 2015 Edward J. Los - For Transient Candidate searches, split the lightcurve if there is a MAX_TC_DAYS gap
 * Jul 17, 2015 Edward J. Los - Correct peakEventCount2 and add peakEventCount3
 *                              In ProcessProperMotions, use the catalog position if no estimate is available
 * Jul 21, 2015 Edward J. Los - Implement Josh Grindlay's changes of 7/19/15 10:03 PM:
 *                                Reject DEFECT, QUALITY_SPECTRA, and QUALITY_COLOR
 * Aug  1, 2015 Edward J. Los - Include FILTER_AFLAG_BIN9 and FILTER_AFLAG_MULTIPLE_NONE in transient searches
 *                              Stop considering stars with MIN_TC_POINTS = 3  <= npoints <= MAX_TC_POINTS = 50
 *                              Redefine nearbyREFflag to be set if (1) A 10th magnitude (MAX_TC_BRIGHT_MAG) or brighter star exists within 360 arcsec (MAX_TC_BRIGHT_ARCSEC) of the transient or
 *                                 (2) A catalog star within 90 arcsec (MAX_TC_NEARBY_ARCSEC) is within 2 mags (MAX_TC_NEARBY_MAG) of the average or peak brightness of the good points in the lightcurve (George used all points)
 * Aug 10, 2015 Edward J. Los - Add PARAMETERVECTOR_NPOINTS to determine peakExcess = npoints - (points outside flare candidate)
 *                              Add "catalogDistance" histogram of distances of DASCH objects from catalog objects
 *                              Correct flare search for matched catalog objects.
 * Aug 30, 2015 Edward J. Los - Restore DEFECT to transient candidate search
 * Sep 15, 2015 Edward J. Los - Rename "peakExcess" to "peakExtra"
 *                              Add peakCountNF, peakCountWF, peakOutside, peakOutsideNF, and peakOutsideWF
 *                              Change weighting of DASCH position calculations from dradRMS2 to sqr(dradRMS2)
 *                              Add peakRA and peakDec with 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window
 *                              Add peakDefectCount for the number of lightcurve defects in the selected transient window.
 * Sep 18, 2015 Edward J. Los   Do not set nearbyREFflag = 1 for peaks < 0.5 mag
 * Sep 21, 2015 Edward J. Los   Back out nearbyREFflag changes made after Jul 13, 2015
 *                              Correct computation of peakCount and peakDefectCount
 *                              Add peakNearbyDistance for the distance to nearby stars that are within 2 magnitudes of the brightest or average magnitude of the transient candidate
 * Sep 22, 2015 Edward J. Los   Correct a bug in finding peakNearbyDistance - look for the smallest nearby star that meets the criteria.
 * Oct  2, 2015 Edward J. Los   Add INCLUDE_TRANSIENT_DEFECTS to control whether defects are considered in transients
 *                              Add INCLUDE_TRANSIENT_MULTIPLE to control whether multiple expousres are considered in transients
 * Oct  4, 2015 Edward J. Los - Add peakMultipleCount for the number of multiple exposure plates in the selected transient window
 * Oct 13, 2015 Edward J. Los - Add INCLUDE_TRANSIENT_DRADBIN to test Sumin's nova in Baade's window: S301233173413
 * Oct 26, 2015 Edward J. Los - Add gsc_bin_index for the current gsc bin of the transient
 *                              Add peakDradRMS2 for the dradRMS2 values added in quadrature
 *                              Reject any image from mf plates at Bloemfontain (mf10104+) where the image is within 0.86
 *                                degrees of the top or bottom. See memo of Tue 9/29/15 11:02 AM
 *                              When defects are enabled, include the defects in the calculation of peakDradRMS and peakRA and peakDec
 * Oct 29, 2015 Edward J. Los - Add peakDradRMS3 for the drad rms of points in the selected 90 day window relative to peakRA and peakDec
 *                              Do not generate the timestamp from the current time - have the user pick the string.
 * Dec 24, 2015 Edward J. Los - Correct ProcessLightcurveParameters() console output for multiprocessor operation
 * Jan  5, 2015 Edward J. Los - Add pGscImageTable to pFileCommon for ProcessMatchedImages
 * Feb 26, 2016 Edward J. Los - Add GetFullQuality() and add it to WriteStarbaseRecord() to make all of the web-based quality bits established.
 * Dec 23, 2016 Edward J. Los - Add seriesColortermTable and seriesColorflagTable.  Update GetPlateQualityMask to return these values
 *                              If a selected measurement is for a catalog star with a bad color or magnitude, mark it with badcolorflag and count with peakBadColorCount
 *                              If a selected measurement is for a catalog star with a good color and Stdmag but is in a plate bin with a bad colorterm, mark it with badcolorflag and count with peakBadColorCount
 *                              If a selected unmatched measurement is from a plate with a bad or non-blue colorterm, mark it with badcolorflag and count with peakBadColorCount
 *                              Correct a bug in GetPlateQuality mask that used only the GSC2.3.2 database entries.
 *                              Set the Stdmag and color to 99.0 if a named star is not in the binary catalog database.
 * Jan  6, 2017 Edward J. Los   Correct the catalog star Stdmag comparison with the peak and average brightness of the transient
 *                              Add INCLUDE_TRANSIENT_ONE_MULTIPLE to accept a single multiple exposure plate
 * Feb 22, 2017 Edward J. Los   Add INCLUDE_TRANSIENT_LIMITING_MAG_BACKGROUND to include limiting magnitudes and high background
 *                               use this conditional to back out Sumin's request of Thu 12/15/11 1:52 PM
 * Mar  8, 2017 Edward J. Los   Move the readiing of flare candidate limiting magnitudes to this module.
 *                              Calculate peakLimitingYears and peakLimitingPoints
 *                              For transient searches, limit the id_table output to flare candidates only
 * Mar 13, 2017 Edward J. Los   Add FILTER_AFLAG_LOCAL_RMS to INCLUDE_TRANSIENT_LIMITING_MAG_BACKGROUND group
 *                              Remove FILTER_AFLAG_LIMITING_MAG from INCLUDE_TRANSIENT_LIMITING_MAG_BACKGROUND Rename latter to INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND
 * Mar 24, 2017 Edward J. Los   Add GetStarEntry2() to access the new starcatalog table instead of the old stars<n> tables.
 * Mar 28, 2017 Edward J. Los   Add peakEvaluation bitmap to indicate reasons for rejection of transient candidate flares
 * Apr  4, 2017 Edward J. Los   Split PEAKEVALUATION_DATERANGE into PEAKEVALUATION_TOOSHORTDATERANGE, PEAKEVALUATION_TOOLONGDATERANGE, and PEAKEVALUATION_PEAKEVALUATION_INSUFFICIENT_POINTS3
 *                              Add PEAKEVALUATION_EXCLUSIONZONE, PEAKEVALUATION_INSUFFICIENT_POINTS2, and PEAKEVALUATION_TOOMANYPOINTS
 * Apr 24, 2017 Edward J. Los   Add fatal flag to GetSeriesId().  Revert to a slow search if the fatal flag is clear.
 * Apr 28, 2017 Edward J. Los   Remove QUALITY_WEDGE for flare candidate searches
 * Jul  1, 2017 Edward J. Los   Support pFileCommon->provideWebSummary for website statistics only
 * Jul 24, 2017 Edward J. Los   Remove xbins and ybins from photglobal because of variable numbers of smoothing bins depending on plate size
 * Aug  3, 2017 Edward J. Los   Add debugMode for GetPhotLocalBin();
 * Dec 20, 2017 Edward J. Los   Add GetStarEntry3 to access the identifier table
 * Jan 23, 2018 Edward J. Los   Support the merged experimental table: add catalogNumber to WriteStarbaseRecord.
 * Mar  5, 2018 Edward J. Los   Add peakNoDefectMag
 * Mar 30, 2018 Edward J. Los   Add combined catalog support.
 * Apr  2, 2018 Edward J. Los   Add peakLongOutburst (Josh memo of Sun, 1 Apr 2018 20:30:06)
 *                              Add LONG_TC_DAYS = 365.0, LONG_TC_POINTS = 10, LONG_TC_MAGNITUDE = 1.0.
 *                              Declare a peakLongOutburst for peakYear of the average of LONG_TC_POINTS or more since (peakYear - 365.0)
 *                               is less than LONG_TC_MAGNITUDE of the average of LONG_TC_POINTS or more until (peakYear + 365.0).
 *                               Add ProcessLongTransientCandidates to find these long outbursts.
 * Apr 17, 2018 Edward J. Los   Correct logic in CombineMultiple for
 *                               (1) both AFLAGS bits set,
 *                               (2) duplicates in the substitution table, and
 *                               (3) clear substitiution table between invocations.
 * May  2, 2018 Edward J. Los   Split LONG_TC_DAYS into
 *                                    LONG_TC_FLARE_DAYS 365.0 days the flare is active
 *                                    LONG_TC_PRE_DAYS (5*365.0) days integration before the transient
 *                                    LONG_TC_SKIP_DAYS (2*365.0) days to skip after finding a transient
 * May 29, 2018 Edward J. Los   Add support for the GAIA catalog
 * Jul 20, 2018 Edward J. Los   Replace the series hash table with a search tree
 * Aug 13, 2018 Edward J. Los   Define enableRematch to optimize location of transients
 *                              Add ProcessNoneImagesY to handle rematch
 * Sep 22, 2018 Edward J. Los   Fix illegal RA for GetDASCHNumber
 * Sep 24, 2018 Edward J. Los   Handle zero corrupted Date from database
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 * Jan  4, 2019 Edward J. Los   Add verbose to LocateNoneImages and LoadNoneImages to display the files searched.
 * Dec  9, 2020 Edward J. Los   Add InitFlareCandidateList
 * Jan 20, 2021 Edward J. Los   Add gsc_bin_index to FLARECANDIDATE
 */
