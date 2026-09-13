// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* parseyield.c
 *
 *  This program extracts the yield per series from the pipeline run logs
 *
 *  Usage:  parseyield -l  <list of logs>  test with run_932.log for i24330
 *
 *           parseyield  -v  -l /home/scanner/Pipeline/los.list -o /home/scanner/Pipeline/los.db
 *
            ./losyield.csh
            parseyield -v -l losyield.txt

            parseyield -v -l /home/scanner/Pipeline/aibifi_run.list
           
 *
 *     To find runs
       cd ~/backup
       find . -name "run*.log" -exec egrep "/ai|/bi0|/bi01|/fa" {} \; -print > ~/junk/runyield.log
       grep "\./" ~/junk/runyield.log | sort -u 
 *
 * cc -ggdb -O2  -I/usr/include/mysql  -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm parseyield.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib/mysql  -l mysqlclient -o parseyield
 * 
 *
 * Oct 20, 2017 Edward J. Los - Initial version  
 * Nov 15, 2017 Edward J. Los - Correct for partial runs
 * Jan 23, 2018 Edward J. Los - Add quotes filenames to look for stray spaces
 * Apr 13, 2018 Edward J. Los - Add additional group 3 series
 * May 28, 2018 Edward J. Los - Add gaia support
 * Jun  8, 2018 Edward J. Los - Stop the seriesId from being overwritten
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 * Dec  9, 2019 Edward J. Los - for FitsMosaic, get the mosaic number from the output file, not the input
 * Dec 10, 2019 Edward J. Los - Add up the output column in all "Illegal" lines to estimate the total photometry space
 *                              Remove gaia support
 */   


#include <math.h>
#include <string.h>
#include "table.h"
#include "time.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER     1024
#define MAX_TIME_LENGTH 28
#define MAX_SOLUTION_NUMBER 10
#define MAX_PLATE_COUNT 500000
/* #define DEBUG_PARSE 1 */

#if 0
#define DEBUG_SERIES "i"
#define DEBUG_PLATENUMBER 46253
#endif

/* Table of all mosaics */
typedef struct _platestats {
	char series[MAX_SERIES_STRING];
  int seriesId;
	int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  long long outputCount;
  int totalMosaics;
  int totalAstrometry;  /* Astrometry.net succeeded */
  int totalWCS16tools;       /* imwcs succeeded 16 bin */
  int totalWCS01tools;       /* imwcs succeeded  1 bin */
  int totalGSCIllegal;        /* pipeline succeeded for gsc2.3.2 */
  int totalAPASSIllegal;   /* piperine succeeded for apass */
  int totalATLASIllegal;
  long long totalGSCoutputCount;
  long long totalAPASSoutputCount;
  long long totalATLASoutputCount;
} PLATESTATS,*PPLATESTATS;

/* Table of all series */
typedef struct _seriesstats {
  char series[MAX_SERIES_STRING];
  int seriesId;
  int totalMosaics;
  int totalAstrometry;
  int totalWCS16tools;
  int totalWCS01tools;
  int totalAPASSIllegal;
  int totalATLASIllegal;
  int totalGSCIllegal;
  long long totalOutputCount;
  long long totalGSCoutputCount;
  long long totalAPASSoutputCount;
  long long totalATLASoutputCount;

} SERIESSTATS,*PSERIESSTATS;


typedef struct _pparsecommon {
	FILE *out_handle;
	time_t startTime;
	time_t totalStepTime;
	time_t totalTime;
	time_t totalChargeTime; 
	time_t totalTimestampTime;
	int verbose;
  int partialFlag;
	int maxLineLen;
  int totalPlates;
	PPLATESTATS plateStatsTable;
  SERIESSTATS seriesTable[MAX_SERIES];
  int maxSeriesId;
} PARSECOMMON,*PPARSECOMMON;

int GetSeriesStatsId(PPARSECOMMON pParseCommon,char* series)
{
  int seriesId;
  PSERIESSTATS pSeriesStats;

  if (strcmp(series,"mask") == 0) {
    return(-1);
  }
  for (seriesId = 0; seriesId < pParseCommon->maxSeriesId; seriesId++) {
    pSeriesStats = &pParseCommon->seriesTable[seriesId];
    if (strcmp(pSeriesStats->series,series) == 0) {
      return(seriesId);
    }
  }
#ifdef DEBUG_SERIES
      if (strcmp(series,DEBUG_SERIES) == 0) {
        printf("line %d at %s\n",__LINE__,series);
      }
#endif /* DEBUG_SERIES */

  pSeriesStats = &pParseCommon->seriesTable[pParseCommon->maxSeriesId];
  strcpy(pSeriesStats->series,series);
  seriesId = pParseCommon->maxSeriesId;
  pSeriesStats->seriesId = seriesId;
  if (pParseCommon->maxSeriesId >= (MAX_SERIES-1)) {
    printf("ERROR: MAX_SERIES exceeded in line %d\n",__LINE__);
    exit(-1);
  }
  pParseCommon->maxSeriesId++;
  return(seriesId);
}
PPLATESTATS GetPlateStats(PPARSECOMMON  pParseCommon,char *series,int plateNumber,int mosaicNumber)
{
  int plateIndex;
  PPLATESTATS pPlateStats;
  for (plateIndex = 0; plateIndex < pParseCommon->totalPlates; plateIndex++) {
    pPlateStats = &pParseCommon->plateStatsTable[plateIndex];
    if ((pPlateStats->plateNumber == plateNumber) &&
        (pPlateStats->mosaicNumber == mosaicNumber) &&
        (strcmp(pPlateStats->series,series) == 0))
      return(pPlateStats);
  }
  if (pParseCommon->totalPlates >= (MAX_PLATE_COUNT-1)) {
    printf("ERROR: MAX_PLATE_COUNT exceeded in line %d\n",__LINE__);
    exit(-1);
  }
  pPlateStats = &pParseCommon->plateStatsTable[pParseCommon->totalPlates];
  strcpy(pPlateStats->series,series);
  pPlateStats->plateNumber = plateNumber;
  pPlateStats->mosaicNumber = mosaicNumber;
  pParseCommon->totalPlates++;
  return(pPlateStats);
}


int ProcessLogBuffer(PPARSECOMMON pParseCommon, char* log_name) {
	char *inBuffer;
	char copyLine[MAX_BUFFER];
	char inLine[MAX_BUFFER];
	char parseLine[MAX_BUFFER];
	char *parsePtr1;
	char *parsePtr2;
	size_t tmpBytesCopied;
	int parseCount;
	int solutionIndex;
	char token[MAX_BUFFER];
	char startToken[MAX_BUFFER];
	int nvals;
  long long outputCount;
	FILE *log_handle = NULL;
	int lineLen;
	int lineCount = 0;
	int recStatCount = 0;
	int parseState = 1; /* 1 = idle; 2 = FitsMosaic; 3 astrometry.csh; 4 WCS16; 5 WCS01 */
	char *charPtr;
	char *charPtr2;
  char *charPtr3;
  char *charPtr4;
  char *charPtr5;
	char series[MAX_SERIES_STRING];
	int plateNumber;
	int processCount = 0;
	time_t lastTime = 0;
	time_t prevTime = 0;
	PPLATESTATS pPlateStatsTable = pParseCommon->plateStatsTable;
	double realSeconds;
	int intSeconds;
	int forwardFlag;
	int saveStep;
	int tmpSolutionNumber;
	int catalogNumber;
  int doneFlag = 0; /* set when we encounter a "_s" in "Illegal" */
  PLATESTATS plateStats;
  PPLATESTATS pTempPlateStats = &plateStats; /* Note: access to this structure is controlled by parseState */
  PPLATESTATS pPlateStats;

	printf("Processing '%s'\n",log_name);
	log_handle = fopen(log_name,"rt");
	startToken[0] = 0;
	if (log_handle == NULL) {
		printf("ERROR: Could not open list file '%s'\n",log_name);
		return(-1);
	}
	while (1) {
		inBuffer = fgets(inLine,MAX_BUFFER,log_handle);
		if (inBuffer == NULL) {
			break;
			
		}
		lineLen = strlen(inBuffer);
		if (lineLen > pParseCommon->maxLineLen) {
			pParseCommon->maxLineLen = lineLen;
		}


		/* Trim off the carriage return */
		if (inBuffer[lineLen-1] == 10) {
			inBuffer[lineLen-1] = 0;
			lineLen--;
		}
		/* Trim off the line feed */
		if (inBuffer[lineLen-1] == 13) {
			inBuffer[lineLen-1] = 0;
			lineLen--;
		}
		lineCount++;
#if 0
		if (lineCount == 984) {
			printf("line %d At lineCount %d %s\n",__LINE__,lineCount,inLine);
		}
#endif
    strcpy(copyLine,inBuffer);
		if ((doneFlag == 0) && (strstr(inBuffer,"CHILD")) != 0) {
			charPtr = inBuffer+strlen("CHILD:");
			while ((*charPtr == ' ') && (*charPtr != 0)) {
				charPtr++;
			}
			charPtr2 = charPtr;
			while ((*charPtr2 != ' ') && (*charPtr2 != 0)) {
				charPtr2++;
			}
			*charPtr2 = 0;
			strcpy(token,charPtr);
      if (strstr(copyLine,"FitsMosaic") != NULL) {
        memset(pTempPlateStats,0,sizeof(PLATESTATS));
        charPtr2++;
        charPtr3 = strstr(charPtr2,"/Mosaics/");
        if (charPtr3 == NULL) {
          printf("ERROR: line %d Failed to parse lineCount %d %s of %s\n",__LINE__,lineCount,copyLine,log_name);
          fclose(log_handle);
          return(-1);
        }
        charPtr3 += strlen("/Mosaics/");
        charPtr4 = strstr(charPtr3,"/");
        if (charPtr4 == NULL) {
          printf("ERROR: line %d Failed to parse lineCount %d %s of %s\n",__LINE__,lineCount,copyLine,log_name);
          fclose(log_handle);
          return(-1);
        }
        charPtr4++;
        charPtr3 = strstr(charPtr4,"/");
        if (charPtr3 == NULL) {
          printf("ERROR: line %d Failed to parse lineCount %d %s of %s\n",__LINE__,lineCount,copyLine,log_name);
          fclose(log_handle);
          return(-1);
        }
        charPtr3++;
        charPtr4 = strstr(charPtr3,".");
        if (charPtr4 == NULL) {
          printf("ERROR: line %d Failed to parse lineCount %d %s of %s\n",__LINE__,lineCount,copyLine,log_name);
          fclose(log_handle);
          return(-1);
        }
        *charPtr4 = 0;
        strcpy(token,charPtr3);


        strcpy(startToken,token);
        if (pParseCommon->verbose) {
          printf("Have token %s in line %d of %s\n",startToken,lineCount,log_name);
        }
        parseState = 2;
#ifdef DEBUG_PARSE
        printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        if (ParseFilename(startToken,pTempPlateStats->series,&pTempPlateStats->plateNumber,&pTempPlateStats->mosaicNumber,&pTempPlateStats->binning,&pTempPlateStats->rotation) == 0) {
          if ((strstr(copyLine,"FitsLinearity") != 0) ||
              (strstr(copyLine,"FitsMedian") != 0)) {
#ifdef DEBUG_PARSE
            parseState = 1; /* Flush rest */
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
          } else {
            printf("ERROR: Failed to parse startToken %s in line %d %s of %s\n",startToken,lineCount,copyLine,log_name);
            fclose(log_handle);
            return(-1);
          }
        }
        pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);

#ifdef DEBUG_SINGLE
        if (strstr(startToken,DEBUG_SINGLE) != NULL) {
          printf("%s\n",copyLine);
        }
#endif /* DEBUG_SINGLE */
        if (parseState == 2) {
          pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
          if ((pTempPlateStats->seriesId < 0) || (pTempPlateStats->seriesId >= MAX_SERIES)) {
            parseState = 1;
#ifdef DEBUG_PARSE
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
          }
        }
      } if (strstr(copyLine,"find_astrometry.csh") != NULL) {
        memset(pTempPlateStats,0,sizeof(PLATESTATS));
        strcpy(startToken,token);
        if (pParseCommon->verbose) {
          printf("Have token %s in line %d of %s\n",startToken,lineCount,log_name);
        }
        parseState = 3;
#ifdef DEBUG_PARSE
        printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        if (ParseFilename(startToken,pTempPlateStats->series,&pTempPlateStats->plateNumber,&pTempPlateStats->mosaicNumber,&pTempPlateStats->binning,&pTempPlateStats->rotation) == 0) {
          if ((strstr(copyLine,"FitsLinearity") != 0) ||
              (strstr(copyLine,"FitsMedian") != 0)) {
#ifdef DEBUG_PARSE
            parseState = 1; /* Flush rest */
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
          } else {
            printf("ERROR: Failed to parse startToken %s in line %d %s of %s\n",startToken,lineCount,copyLine,log_name);
            fclose(log_handle);
            return(-1);
          }
        }
        pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
#ifdef DEBUG_SINGLE
        if (strstr(startToken,DEBUG_SINGLE) != NULL) {
          printf("%s\n",copyLine);
        }
#endif /* DEBUG_SINGLE */
        if (parseState == 3) {
          pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
          if ((pTempPlateStats->seriesId < 0) || (pTempPlateStats->seriesId >= MAX_SERIES)) {
            parseState = 1;
#ifdef DEBUG_PARSE
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
          }
        }
      } if (strstr(copyLine,"Starting AstrometryWCS") != NULL) {
        memset(pTempPlateStats,0,sizeof(PLATESTATS));
        strcpy(startToken,token);
        if (pParseCommon->verbose) {
          printf("Have token %s in line %d of %s\n",startToken,lineCount,log_name);
        }
        parseState = 4;
#ifdef DEBUG_PARSE
        printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        if (ParseFilename(startToken,pTempPlateStats->series,&pTempPlateStats->plateNumber,&pTempPlateStats->mosaicNumber,&pTempPlateStats->binning,&pTempPlateStats->rotation) == 0) {
          if ((strstr(copyLine,"FitsLinearity") != 0) ||
              (strstr(copyLine,"FitsMedian") != 0)) {
#ifdef DEBUG_PARSE
            parseState = 1; /* Flush rest */
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
          } else {
            printf("ERROR: Failed to parse startToken %s in line %d %s of %s\n",startToken,lineCount,copyLine,log_name);
            fclose(log_handle);
            return(-1);
          }
        }
        pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
#ifdef DEBUG_SINGLE
        if (strstr(startToken,DEBUG_SINGLE) != NULL) {
          printf("%s\n",copyLine);
        }
#endif /* DEBUG_SINGLE */
        if (parseState == 4) {
          pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
          if ((pTempPlateStats->seriesId < 0) || (pTempPlateStats->seriesId >= MAX_SERIES)) {
            parseState = 1;
#ifdef DEBUG_PARSE
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
          }
        }
      } if ((strstr(copyLine,"fulldb.csh") != NULL) ||
            ((pParseCommon->partialFlag != 0) && (strstr(copyLine,"full.csh") != NULL))) {
        memset(pTempPlateStats,0,sizeof(PLATESTATS));
        strcpy(startToken,token);
        if (pParseCommon->verbose) {
          printf("Have token %s in line %d of %s\n",startToken,lineCount,log_name);
        }
        parseState = 5;
#ifdef DEBUG_PARSE
        printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        if (ParseFilename(startToken,pTempPlateStats->series,&pTempPlateStats->plateNumber,&pTempPlateStats->mosaicNumber,&pTempPlateStats->binning,&pTempPlateStats->rotation) == 0) {
          if ((strstr(copyLine,"FitsLinearity") != 0) ||
              (strstr(copyLine,"FitsMedian") != 0)) {
#ifdef DEBUG_PARSE
            parseState = 1; /* Flush rest */
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
          } else {
            printf("ERROR: Failed to parse startToken %s in line %d %s of %s\n",startToken,lineCount,copyLine,log_name);
            fclose(log_handle);
            return(-1);
          }
        }
        pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);

#ifdef DEBUG_SINGLE
        if (strstr(startToken,DEBUG_SINGLE) != NULL) {
          printf("%s\n",copyLine);
        }
#endif /* DEBUG_SINGLE */
#if 0
        if (strcmp(pTempPlateStats->series,"bi") == 0) {
          printf("At series %s\n",pTempPlateStats->series);
        }

#endif

        if (parseState == 5) {
          pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
          if ((pTempPlateStats->seriesId < 0) || (pTempPlateStats->seriesId >= MAX_SERIES)) {
            parseState = 1;
#ifdef DEBUG_PARSE
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
          } else {
            parseState = 1;
#ifdef DEBUG_PARSE
            printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */   
            pPlateStats = GetPlateStats(pParseCommon,pTempPlateStats->series,pTempPlateStats->plateNumber,pTempPlateStats->mosaicNumber);
            pPlateStats->totalWCS01tools++;
            if ((pPlateStats->seriesId != 0) && 
                (pTempPlateStats->seriesId != pPlateStats->seriesId)) {
              printf("ERROR: line %d seriesId overwrite %d %d for %s%05d\n",__LINE__,pTempPlateStats->seriesId,pPlateStats->seriesId,pPlateStats->series,pPlateStats->seriesId);
              exit(-1);
            }
            pPlateStats->seriesId = pTempPlateStats->seriesId;
          }
        }
      }
    } else if ((doneFlag == 0) && (parseState == 2)) {
#if 0
      printf("%s\n",copyLine);
#endif
      if (strstr(copyLine,"***L") != 0) {
        pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
        if ((pTempPlateStats->seriesId < 0) || (pTempPlateStats->seriesId >= MAX_SERIES)) {
          parseState = 1;
#ifdef DEBUG_PARSE
          printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        } else {
          parseState = 1;
          pPlateStats = GetPlateStats(pParseCommon,pTempPlateStats->series,pTempPlateStats->plateNumber,pTempPlateStats->mosaicNumber);
          pPlateStats->totalMosaics++;
#ifdef DEBUG_SERIES
          if ((pPlateStats->plateNumber == DEBUG_PLATENUMBER) &&
              (strcmp(pPlateStats->series,DEBUG_SERIES) == 0)) {
            printf("line %d at %s%05d_%02d Mosaics %d Astrometry %d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,pPlateStats->mosaicNumber,pPlateStats->totalMosaics,pPlateStats->totalAstrometry);
          }
#endif /* DEBUG_SERIES */

          if ((pPlateStats->seriesId != 0) && 
              (pTempPlateStats->seriesId != pPlateStats->seriesId)) {
            printf("ERROR: line %d seriesId overwrite %d %d for %s%05d\n",__LINE__,pTempPlateStats->seriesId,pPlateStats->seriesId,pPlateStats->series,pPlateStats->seriesId);
            exit(-1);
          }
          pPlateStats->seriesId = pTempPlateStats->seriesId;
#ifdef DEBUG_PARSE
          printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        }
      }
    } else if ((doneFlag == 0) && (parseState == 3)) {
#if 0
      printf("%s\n",copyLine);
#endif
      if (strstr(copyLine,"Solution found at star") != 0) {
        pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
        if ((pTempPlateStats->seriesId < 0) || (pTempPlateStats->seriesId >= MAX_SERIES)) {
          parseState = 1;
#ifdef DEBUG_PARSE
          printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        } else {
          parseState = 1;
          pPlateStats = GetPlateStats(pParseCommon,pTempPlateStats->series,pTempPlateStats->plateNumber,pTempPlateStats->mosaicNumber);
          pPlateStats->totalAstrometry++;
#ifdef DEBUG_SERIES
          if ((pTempPlateStats->plateNumber == DEBUG_PLATENUMBER) &&
              (strcmp(pTempPlateStats->series,DEBUG_SERIES) == 0)) {
            printf("line %d at %s%05d_%02d Mosaics %d Astrometry %d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,pPlateStats->mosaicNumber,pPlateStats->totalMosaics,pPlateStats->totalAstrometry);
          }
#endif /* DEBUG_SERIES */
          pPlateStats->seriesId = pTempPlateStats->seriesId;
#ifdef DEBUG_PARSE
          printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        }
      }
    } else if ((doneFlag == 0) && (parseState == 4)) {
#if 0
      printf("%s\n",copyLine);
#endif
      if (strstr(copyLine,"_16.fit: written successfully.") != 0) {
        pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
        if ((pTempPlateStats->seriesId < 0) || (pTempPlateStats->seriesId >= MAX_SERIES)) {
          parseState = 1;
#ifdef DEBUG_PARSE
          printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        } else {
          parseState = 1;
          pPlateStats = GetPlateStats(pParseCommon,pTempPlateStats->series,pTempPlateStats->plateNumber,pTempPlateStats->mosaicNumber);
          pPlateStats->totalWCS16tools++;
          if ((pPlateStats->seriesId != 0) && 
              (pTempPlateStats->seriesId != pPlateStats->seriesId)) {
            printf("ERROR: line %d seriesId overwrite %d %d for %s%05d\n",__LINE__,pTempPlateStats->seriesId,pPlateStats->seriesId,pPlateStats->series,pPlateStats->seriesId);
            exit(-1);
          }
          pPlateStats->seriesId = pTempPlateStats->seriesId;
#ifdef DEBUG_PARSE
          printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
        }
      }
    }
    if (strstr(copyLine,"Illegal") != 0) {
      parseState = 1;
#ifdef DEBUG_PARSE
      printf("line %4d parseState is %d\n",__LINE__,parseState);
#endif /* DEBUG_PARSE */
      charPtr = strstr(copyLine," for ");
      if (charPtr != NULL) {
        charPtr += 5;
        if (strstr(charPtr,"_s") != 0) {
          /* We are done */
          doneFlag = 1;
        }
        if (ParseFilename(charPtr,pTempPlateStats->series,&pTempPlateStats->plateNumber,&pTempPlateStats->mosaicNumber,&pTempPlateStats->binning,&pTempPlateStats->rotation) != 0) {
          pTempPlateStats->seriesId = GetSeriesStatsId(pParseCommon,pTempPlateStats->series);
#ifdef DEBUG_SERIES
          if ((pTempPlateStats->plateNumber == DEBUG_PLATENUMBER) &&
              (strcmp(pTempPlateStats->series,DEBUG_SERIES) == 0)) {
            printf("line %d at %s%05d\n",__LINE__,pTempPlateStats->series,pTempPlateStats->plateNumber);
          }
#endif /* DEBUG_SERIES */
          pPlateStats = GetPlateStats(pParseCommon,pTempPlateStats->series,pTempPlateStats->plateNumber,pTempPlateStats->mosaicNumber);
          outputCount = -1;
          charPtr5 = strstr(copyLine,"output "); 
          if (charPtr5 != NULL) {
            charPtr5 += 7;
            nvals = sscanf(charPtr5,"%lld",&outputCount);
            if (nvals != 1) {
              outputCount = -1;
            }
          }
          

          if (outputCount > 0) {
            if (strstr(charPtr,"_apass") != 0) {
                pPlateStats->totalAPASSoutputCount += outputCount;
            } else if (strstr(charPtr,"_atlas") != 0) {
                pPlateStats->totalATLASoutputCount += outputCount;
            } else {
                pPlateStats->totalGSCoutputCount += outputCount;
            }
          }
          if (doneFlag == 0) {
            if (strstr(charPtr,"_apass") != 0) {
              pPlateStats->totalAPASSIllegal++;
            } else if (strstr(charPtr,"_atlas") != 0) {
              pPlateStats->totalATLASIllegal++;
            } else {
              pPlateStats->totalGSCIllegal++;
              if ((outputCount >= 0) && (pPlateStats->outputCount < outputCount)) {
                pPlateStats->outputCount = outputCount;
              }
            }
          }
        }
        if ((pPlateStats->seriesId != 0) && 
            (pTempPlateStats->seriesId != pPlateStats->seriesId)) {
          printf("ERROR: line %d seriesId overwrite %d %d for %s%05d\n",__LINE__,pTempPlateStats->seriesId,pPlateStats->seriesId,pPlateStats->series,pPlateStats->seriesId);
          exit(-1);
        }
        pPlateStats->seriesId = pTempPlateStats->seriesId;
      }
        
    }
  

  }

	printf("Finished processing %s lines %d recStatCount %d\n",log_name,lineCount,recStatCount);
	fclose(log_handle);
	return(0);
}

int main(int argc,char *argv[])
{
	char *argstr;
	char *inBuffer;
	char inLine[MAX_BUFFER];
	char list_name[MAX_INPUT_NAME];
	char out_name[MAX_INPUT_NAME];
	FILE *list_handle = NULL;
	int errorFlag = 0;
	char cmdchar;
	int lineLen;
	time_t curTime;
	PARSECOMMON parsecommon;
	PPARSECOMMON pParseCommon = &parsecommon;
	int stepIndex;

	time_t pTotalStepTime = 0;
	int stepCounter;
	double percentage;
	double percentage2;
	time_t averageTime;
	time_t totalAverageTime = 0;
	int solutionIndex;
  int plateIndex;
  int seriesId;
  int iteration;
  PPLATESTATS pPlateStats;
  PSERIESSTATS pSeriesStats;
  int totalMosaics;
  int totalGSCIllegal;
  int totalAPASSIllegal;
  int totalATLASIllegal;
  int totalAstrometry;
  int totalWCS01tools;
  int totalWCS16tools;
  long long totalOutputCount;
  long long tempOutputCount;
  long long totalGSCoutputCount;
  long long totalAPASSoutputCount;
  long long totalATLASoutputCount;
  long long allOutputCount;
  double totalBytes;
  double totalMegaBytes;
  double totalTeraBytes;

	list_name[0] = 0;
	out_name[0] = 0;
	memset(pParseCommon,0,sizeof(PARSECOMMON));
  pParseCommon->maxSeriesId = 1;
  pParseCommon->plateStatsTable = (PPLATESTATS)calloc(MAX_PLATE_COUNT,sizeof(PLATESTATS));
  if (pParseCommon->plateStatsTable == NULL) {
    printf("ERROR: failed to allocate plateStatsTable\n");
    exit(-1);
  }

	/* Loop through the arguments */
	for (argv++; --argc > 0; argv++) {
		argstr = *argv;
		/* Decode arguments */
		if (argstr[0] != '-') {
			errorFlag = 1;
			printf("ERROR: stray argument %s\n",argstr);
		} else {
			while ((cmdchar = *++argstr) != 0) {
				switch(cmdchar) {
				case 'v':
				case 'V':
					pParseCommon->verbose = 1;
					break;

				case 'p':
				case 'P':
					pParseCommon->partialFlag = 1;
					break;

				case 'l': /* list file name */
				case 'L':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(list_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;

				case 'o': /* output file name */
				case 'O':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(out_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;


				default:
					printf("* illegal command -%c-",cmdchar);
					errorFlag = 1;
					break;
				}
        
			}
		}
	}



	if (list_name[0] == 0) {
		printf("ERROR: List file name not specified\n");
		errorFlag = 1;
	} else {
		list_handle = fopen(list_name,"rt");
		if (list_handle == NULL) {
			errorFlag = 1;
			printf("Could not open list file %s\n",list_name);
		}
	}
#if 0
	if (out_name[0] == 0) {
		printf("ERROR: Out file name not specified\n");
		errorFlag = 1;
	} else {
		pParseCommon->out_handle = fopen(out_name,"wt");
		if (pParseCommon->out_handle == NULL) {
			errorFlag = 1;
			printf("Could not open out file %s\n",out_name);
		}
	}
#endif


	if (errorFlag) {
		printf("Usage: parseyield [-v] [-p] -l <list of log files> -o <outfile>\n");
		printf("       where -v is the verbose flag\n");
    printf("       where -p is for partial runs not containing a mosaic build\n");
    printf("       NOTE: outfile is currently not supported\n");

		return(-1);
	}

	printf("parseyield of %s %s List Filename %s Output Filename %s partialFlag %d\n",
				 __DATE__,__TIME__,list_name,out_name,pParseCommon->partialFlag);
 

	time(&pParseCommon->startTime);

	while (1) {
		inBuffer = fgets(inLine,MAX_BUFFER,list_handle);
		if (inBuffer == NULL) {
			break;
			
		}
		lineLen = strlen(inBuffer);
		/* Trim off the carriage return */
		if (inBuffer[lineLen-1] == 10) {
			inBuffer[lineLen-1] = 0;
			lineLen--;
		}
		/* Trim off the line feed */
		if (inBuffer[lineLen-1] == 13) {
			inBuffer[lineLen-1] = 0;
			lineLen--;
		}
		if (inBuffer[0] != '#') {
			ProcessLogBuffer(pParseCommon,inBuffer);
		}

#ifdef DEBUG_SERIES_X
    {
      
        if ((pPlateStats->plateNumber == DEBUG_PLATENUMBER) &&
            (strcmp(pPlateStats->series,DEBUG_SERIES) == 0)) {
          printf("line %d at %s%05d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber);
        }
    }
#endif /* DEBUG_SERIES */

	}
  for (plateIndex = 0; plateIndex < pParseCommon->totalPlates; plateIndex++) {
    pPlateStats = &pParseCommon->plateStatsTable[plateIndex];
    pSeriesStats = &pParseCommon->seriesTable[pPlateStats->seriesId];
    if (strcmp(pSeriesStats->series,pPlateStats->series) != 0) {
      printf("ERROR line %d series %s %s does not match\n",__LINE__,pSeriesStats->series,pPlateStats->series);
#if 0
      exit(-1);
#endif
    }

#ifdef DEBUG_SERIES
        if ((pPlateStats->plateNumber == DEBUG_PLATENUMBER) &&
            (strcmp(pPlateStats->series,DEBUG_SERIES) == 0)) {
          printf("line %d at %s%05d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber);
        }
#endif /* DEBUG_SERIES */

    if ((pParseCommon->partialFlag != 0) && 
        (pPlateStats->totalMosaics == 0)) {
      if ((pPlateStats->totalAstrometry != 0) ||
          (pPlateStats->totalWCS16tools != 0) ||
          (pPlateStats->totalWCS01tools != 0) ||
          (pPlateStats->totalGSCIllegal != 0) ||
          (pPlateStats->totalATLASIllegal != 0) ||
          (pPlateStats->totalAPASSIllegal != 0)) {
        pPlateStats->totalMosaics++;
      }
    }


    if (pPlateStats->totalMosaics != 0) {
      pSeriesStats->totalMosaics++;
      if (pPlateStats->totalAstrometry != 0) {
        pSeriesStats->totalAstrometry++;
      }
      if (pPlateStats->totalWCS16tools != 0) {
        pSeriesStats->totalWCS16tools++;
      }
      if (pPlateStats->totalWCS01tools != 0) {
        pSeriesStats->totalWCS01tools++;
      }
      if (pPlateStats->totalGSCIllegal != 0) {
        pSeriesStats->totalGSCIllegal++;
        pSeriesStats->totalOutputCount += pPlateStats->outputCount;
        pSeriesStats->totalGSCoutputCount += pPlateStats->totalGSCoutputCount;
      }
      
      if (pPlateStats->totalAPASSIllegal != 0) {
        pSeriesStats->totalAPASSIllegal++;
        pSeriesStats->totalAPASSoutputCount += pPlateStats->totalAPASSoutputCount;
      }
      if (pPlateStats->totalATLASIllegal != 0) {
        pSeriesStats->totalATLASIllegal++;
        pSeriesStats->totalATLASoutputCount += pPlateStats->totalATLASoutputCount;
      }
    }
  }
  totalGSCoutputCount = 0;
  totalAPASSoutputCount = 0;
  totalATLASoutputCount = 0;
  for (iteration = 0; iteration < 3; iteration++) {
    totalMosaics = 0;
    totalAstrometry = 0;
    totalWCS01tools = 0;
    totalWCS16tools = 0;
    totalGSCIllegal = 0;
    totalAPASSIllegal = 0;
    totalATLASIllegal = 0;
    totalOutputCount = 0;
 

    for (seriesId = 0; seriesId < pParseCommon->maxSeriesId; seriesId++) {
      pSeriesStats = &pParseCommon->seriesTable[seriesId];

#ifdef DEBUG_SERIES
      if (strcmp(pSeriesStats->series,DEBUG_SERIES) == 0) {
        printf("line %d at %s\n",__LINE__,pSeriesStats->series);
      }
#endif /* DEBUG_SERIES */


      if (pSeriesStats->totalMosaics != 0) {
        if ((iteration == 1) && 
            ((strcmp(pSeriesStats->series,"bi") == 0) ||
             (strcmp(pSeriesStats->series,"ai") == 0) ||
             (strcmp(pSeriesStats->series,"al") == 0) ||
             (strcmp(pSeriesStats->series,"ka") == 0) ||
             (strcmp(pSeriesStats->series,"kb") == 0) ||
             (strcmp(pSeriesStats->series,"ke") == 0) ||
             (strcmp(pSeriesStats->series,"kf") == 0) ||
             (strcmp(pSeriesStats->series,"kg") == 0) ||
             (strcmp(pSeriesStats->series,"kge") == 0) ||
             (strcmp(pSeriesStats->series,"kh") == 0) ||
             (strcmp(pSeriesStats->series,"meteor") == 0) ||
             (strcmp(pSeriesStats->series,"pz") == 0) ||
             (strcmp(pSeriesStats->series,"fa") == 0))) {
#if 0
          printf("rejecting iteration %d series %s\n",iteration,pSeriesStats->series);
#endif
          continue;
        }
        if ((iteration == 2) && 
            ((strcmp(pSeriesStats->series,"bi") != 0) &&
             (strcmp(pSeriesStats->series,"ai") != 0) &&
             (strcmp(pSeriesStats->series,"al") != 0) &&
             (strcmp(pSeriesStats->series,"ka") != 0) &&
             (strcmp(pSeriesStats->series,"kb") != 0) &&
             (strcmp(pSeriesStats->series,"ke") != 0) &&
             (strcmp(pSeriesStats->series,"kf") != 0) &&
             (strcmp(pSeriesStats->series,"kg") != 0) &&
             (strcmp(pSeriesStats->series,"kge") != 0) &&
             (strcmp(pSeriesStats->series,"kh") != 0) &&
             (strcmp(pSeriesStats->series,"meteor") != 0) &&
             (strcmp(pSeriesStats->series,"pz") != 0) &&
             (strcmp(pSeriesStats->series,"fa") != 0))) {
#if 0
          printf("rejecting iteration %d series %s\n",iteration,pSeriesStats->series);
#endif
          continue;
        }
#if 0    
        printf("iteration %d series %s\n",iteration,pSeriesStats->series);
#endif
        if (iteration == 0) {
          totalGSCoutputCount += pSeriesStats->totalGSCoutputCount;
          totalATLASoutputCount += pSeriesStats->totalATLASoutputCount;
          totalAPASSoutputCount += pSeriesStats->totalAPASSoutputCount;
        }
        if (pSeriesStats->totalGSCIllegal != 0) {
          tempOutputCount = pSeriesStats->totalOutputCount/pSeriesStats->totalGSCIllegal;
        } else {
          tempOutputCount = 0;
        }
        printf("Series %6s Mosaics %5d Astrometry %5d (%5.1f\%) WCS16 %5d (%5.1f\%) WCS01 %5d (%5.1f\%) GSC %5d (%5.1f\%) APASS %5d (%5.1f\%) ATLAS %5d (%5.1f\%) GSC output %6lld\n",
               pSeriesStats->series,
               pSeriesStats->totalMosaics,
               pSeriesStats->totalAstrometry,
               ((100.*pSeriesStats->totalAstrometry)/pSeriesStats->totalMosaics),
               pSeriesStats->totalWCS16tools,
               ((100.*pSeriesStats->totalWCS16tools)/pSeriesStats->totalMosaics),
               pSeriesStats->totalWCS01tools,
               ((100.*pSeriesStats->totalWCS01tools)/pSeriesStats->totalMosaics),
               pSeriesStats->totalGSCIllegal,
               ((100.*pSeriesStats->totalGSCIllegal)/pSeriesStats->totalMosaics),
               pSeriesStats->totalAPASSIllegal,
               ((100.*pSeriesStats->totalAPASSIllegal)/pSeriesStats->totalMosaics),
               pSeriesStats->totalATLASIllegal,
               ((100.*pSeriesStats->totalATLASIllegal)/pSeriesStats->totalMosaics),
               tempOutputCount);
        totalMosaics += pSeriesStats->totalMosaics;
        totalAstrometry += pSeriesStats->totalAstrometry;
        totalWCS16tools += pSeriesStats->totalWCS16tools;
        totalWCS01tools += pSeriesStats->totalWCS01tools;
        totalGSCIllegal += pSeriesStats->totalGSCIllegal;
        totalAPASSIllegal += pSeriesStats->totalAPASSIllegal;
        totalATLASIllegal += pSeriesStats->totalATLASIllegal;
        totalOutputCount += pSeriesStats->totalOutputCount;
      }
    }
    if (totalMosaics != 0) {
      if (totalGSCIllegal != 0) {
        tempOutputCount = totalOutputCount/totalGSCIllegal;
      } else {
        tempOutputCount = 0;
      }
      printf("Total         Mosaics %5d Astrometry %5d (%5.1f\%) WCS16 %5d (%5.1f\%) WCS01 %5d (%5.1f\%) GSC %5d (%5.1f\%) APASS %5d (%5.1f\%) ATLAS %5d (%5.1f\%) GSC output %6lld\n",
             totalMosaics,
             totalAstrometry,
             ((100.*totalAstrometry)/totalMosaics),
             totalWCS16tools,
             ((100.*totalWCS16tools)/totalMosaics),
             totalWCS01tools,
             ((100.*totalWCS01tools)/totalMosaics),
             totalGSCIllegal,
             ((100.*totalGSCIllegal)/totalMosaics),
             totalAPASSIllegal,
             ((100.*totalAPASSIllegal)/totalMosaics),
             totalATLASIllegal,
             ((100.*totalATLASIllegal)/totalMosaics),
             tempOutputCount);
    }
  }

	if (list_handle != NULL) {
		fclose(list_handle);
	}
	if (pParseCommon->out_handle != NULL) {
		fclose(pParseCommon->out_handle);
	}
	time(&curTime);
	curTime -= pParseCommon->startTime;

  if (pParseCommon->verbose) {
    for (plateIndex = 0; plateIndex < pParseCommon->totalPlates; plateIndex++) {
      pPlateStats = &pParseCommon->plateStatsTable[plateIndex];
      if ((pPlateStats->totalWCS16tools > 0) && (pPlateStats->totalWCS01tools == 0)) {
        printf("plate %s%05d_%02d failed at WCS01\n",pPlateStats->series,pPlateStats->plateNumber,pPlateStats->mosaicNumber);
      }
    }
  }
  allOutputCount = totalGSCoutputCount+totalAPASSoutputCount+totalATLASoutputCount;
  totalBytes = (1.0*allOutputCount) * sizeof( FILESTARIMAGE);
  totalMegaBytes = totalBytes/1000000.;
  totalTeraBytes = totalBytes/(1024.*1024.*1024.*1024.);
  printf("totalOutput: GSC %lld records, APASS %lld records, ATLAS %lld records, all %lld records, %.0f MegaBytes, %.6f Terabytes\n",totalGSCoutputCount,totalAPASSoutputCount,totalATLASoutputCount,allOutputCount,totalMegaBytes,totalTeraBytes);
	printf("Execution Time: %lld seconds. maxLineLen %d total series %d total plates %d\n",
				 curTime,
				 pParseCommon->maxLineLen,
         pParseCommon->maxSeriesId,
         pParseCommon->totalPlates);

  if (pParseCommon->plateStatsTable != NULL) {
    free(pParseCommon->plateStatsTable);
  }


	return(EXIT_SUCCESS);
}
