// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* parsesearch.c
 *
 *  This program reads the log of a multiprocessor search_none job and estimates the total completion
 *
 *  Usage:  parsesearch -l  <run_phot log>  
 *
    parsesearch  -v  -l /home/scanner/junk/phot_605.log -o /home/scanner/junk/los.db
 
   cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm parsesearch.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  plottransientstub.o -o parsesearch
 * 
 *
 * Nov 23, 2018 Edward J. Los - Initial version  
 * Nov 26, 2018 Edward J. Los - check completion of each thread
 */   


#include <math.h>
#include <string.h>
#include "table.h"
#include "time.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER     2048
#define MAX_TIME_LENGTH 28
#define MAX_SOLUTION_NUMBER 10
#define MAX_PLATE_COUNT 500000
#define MAX_RUN 1000
#define TAB_GSC_BIN_INDEX 13
/* #define DEBUG_PARSE 1 */

/* Table of all runs */
typedef struct _runstats {
  int numstarted;
  int rank;
  int doneFlag;
  int subprocessTime;
  int startGscBinIndex;
  int endGscBinIndex;
  int curGscBinIndex;
} RUNSTATS,*PRUNSTATS;




typedef struct _pparsecommon {
	FILE *out_handle;
	time_t startTime;
  int totalRuns;
	int verbose;
	int maxLineLen;
  RUNSTATS runTable[MAX_RUN];
#if 0
	time_t totalStepTime;
	time_t totalTime;
	time_t totalChargeTime; 
	time_t totalTimestampTime;
#endif
} PARSECOMMON,*PPARSECOMMON;

int main(int argc,char *argv[])
{
	char log_name[MAX_INPUT_NAME];
	char out_name[MAX_INPUT_NAME];
  char db_name[MAX_INPUT_NAME];
  char temp_name[MAX_INPUT_NAME];
	PARSECOMMON parsecommon;
	PPARSECOMMON pParseCommon = &parsecommon;
	char *argstr;
	int errorFlag = 0;
	char cmdchar;
	FILE *log_handle = NULL;
	char *inBuffer;
	char inLine[MAX_BUFFER];
	int lineLen;
	time_t curTime;
  int runIndex;
  RUNSTATS templateStats;
  PRUNSTATS pTemplateStats = &templateStats;
  PRUNSTATS pRunStats;
  int statResult;
  struct stat statbuf;
  int entry;
  int doneCount = 0;
#if 0
	int stepIndex;

	time_t pTotalStepTime = 0;
	int stepCounter;
	double percentage;
	double percentage2;
	time_t averageTime;
	time_t totalAverageTime = 0;
	int solutionIndex;
  int seriesId;
  int iteration;
  int totalMosaics;
  int totalGSCIllegal;
  int totalAPASSIllegal;
  int totalGAIAIllegal;
  int totalATLASIllegal;
  int totalAstrometry;
  int totalWCS01tools;
  int totalWCS16tools;
  long long totalOutputCount;
  long long tempOutputCount;
#endif
  char* charPtr;
  int fd;
  off_t seekResult;
  off_t readBytes;
  off_t seekBytes;
  off_t bufferBytes;
  char* recordStart;
  char* recordEnd;
  off_t maxRecordLength = 0;
  char charVal;
  int tabCount;
  int nvals;
  double pctdone;
  int finishedBins = 0;
  int totalBins = 0;


	log_name[0] = 0;
	out_name[0] = 0;
  db_name[0] = 0;
	memset(pParseCommon,0,sizeof(PARSECOMMON));

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

				case 'l': /* list file name */
				case 'L':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(log_name,*++argv,MAX_INPUT_NAME-2);
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



	if (log_name[0] == 0) {
		printf("ERROR: List file name not specified\n");
		errorFlag = 1;
	} else {
		log_handle = fopen(log_name,"rt");
		if (log_handle == NULL) {
			errorFlag = 1;
			printf("Could not open list file %s\n",log_name);
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
		printf("Usage: parsesearch [-v] [-p] -l <list of log files> -o <outfile>\n");
		printf("       where -v is the verbose flag\n");
    printf("       NOTE: outfile is currently not supported\n");

		return(-1);
	}

	printf("parsesearch of %s %s List Filename %s Output Filename %s\n",
				 __DATE__,__TIME__,log_name,out_name);
 

	time(&pParseCommon->startTime);

	while (1) {
		inBuffer = fgets(inLine,MAX_BUFFER,log_handle);
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
    charPtr = strstr(inBuffer,"Output file");
    if ((db_name[0] == 0) && (charPtr == inBuffer)) {
      charPtr += (strlen("Output file")+1);
      if (strlen(charPtr) >= MAX_INPUT_NAME-2) {
        printf("ERROR: db_name %s is too large\n",charPtr);
        exit(-1);
      }
      strcpy(db_name,charPtr);
      charPtr = strstr(db_name,".log");
      if (charPtr == NULL) {
        printf("ERROR: db_name format error\n");
        exit(-1);
      }
      *charPtr = 0;
    }
    charPtr = strstr(inBuffer,"numpending");
    if ((pParseCommon->totalRuns == 0) && (charPtr != NULL)) {
      charPtr += (strlen("numpending")+1);
      nvals = sscanf(charPtr,"%d",&pParseCommon->totalRuns);
      if (nvals != 1) {
        printf("ERROR: failed to decode totalRuns in %s\n",inBuffer);
        exit(-1);
      }
    }

    charPtr = strstr(inBuffer,"PARENT: Created subprocess");
    if (charPtr == inBuffer) {
      charPtr += (strlen("PARENT: Created subprocess")+1);
      memset(pTemplateStats,0,sizeof(RUNSTATS));

      nvals = sscanf(charPtr,"%3d on %3d at %d seconds for bins %9d to %9d",
                     &pTemplateStats->numstarted,
                     &pTemplateStats->rank,
                     &pTemplateStats->subprocessTime,
                     &pTemplateStats->startGscBinIndex,
                     &pTemplateStats->endGscBinIndex);
      if (nvals == 5) {
        if ((pTemplateStats->rank <= 0) || (pTemplateStats->rank > pParseCommon->totalRuns)) {
          printf("ERROR: illegal rank %d totalRuns %d\n",pTemplateStats->rank,pParseCommon->totalRuns);
          exit(-1);
        }
        pRunStats = &pParseCommon->runTable[pTemplateStats->rank-1];
        if (pRunStats->rank != 0) {
          printf("ERROR: duplicate rank %d\n",pTemplateStats->rank);
          exit(-1);
        }
        memcpy(pRunStats,pTemplateStats,sizeof(RUNSTATS));
        
      } else {
        printf("ERROR: failed to parse %s\n",inBuffer);
        exit(-1);
      }

    }
    charPtr = strstr(inBuffer,"Response result:");
    if ((pParseCommon->totalRuns != 0) &&(charPtr == inBuffer)) {
      charPtr = strstr(inBuffer,"entry");
      if (charPtr == NULL) {
        printf("ERROR: can not find entry in  %s\n",inBuffer);
        exit(-1);
      }
      charPtr += strlen("entry")+1;
      nvals = sscanf(charPtr,"%d",&entry);
      if (nvals != 1) {
        printf("ERROR: can not decode entry in  %s\n",inBuffer);
        exit(-1);
      }
      if (entry <= 0 && (entry > pParseCommon->totalRuns)) {
        printf("ERROR: illegal rank %d totalRuns %d\n",pTemplateStats->rank,pParseCommon->totalRuns);
        exit(-1);
      }
      pRunStats = &pParseCommon->runTable[entry-1];
      if (pRunStats->doneFlag == 1) {
        printf("ERROR: doneFlag set twice for %s\n",inBuffer);
        exit(-1);
      }
      pRunStats->doneFlag = 1;
      doneCount++;
    }

	}
  for (runIndex = 0; runIndex < pParseCommon->totalRuns; runIndex++) {
    pRunStats = &pParseCommon->runTable[runIndex];
    sprintf(temp_name,"%s.db_%02d.log",db_name,pRunStats->rank);
    statResult = stat(temp_name,&statbuf);
    if (statResult != 0) {
      printf("ERROR: failed to find file of rank %d: %s\n",pRunStats->rank,temp_name);
    } else {
      printf("Found %s of size %lld\n",temp_name,statbuf.st_size);
      fd = open(temp_name,O_RDONLY,S_IRWXU|S_IRGRP);
      if (fd < 0) {
        printf("ERROR: failed to open %s\n",temp_name);
      } else {
        bufferBytes = MAX_BUFFER;
        seekBytes = statbuf.st_size-MAX_BUFFER-1;
        if (seekBytes < 0) {
          bufferBytes-= seekBytes;
          seekBytes = 0;
        }

        seekResult = lseek(fd,seekBytes,SEEK_SET);
        if (seekResult < 0) {
          printf("ERROR: failed to seek to %lld in %s\n",seekBytes,temp_name);
        } else {
          readBytes = read(fd,inLine,bufferBytes);
          if (readBytes != bufferBytes) {
            printf("ERROR: read %lld of %lld bytes from %s status %d %s\n",readBytes,bufferBytes,temp_name,errno,strerror(errno));
          }
          charPtr = &inLine[0];
          recordStart = NULL;
          recordEnd = NULL;
          while(bufferBytes > 0) {
            charVal = *charPtr;
            if (charVal == '\n') {
              recordStart = recordEnd;
              recordEnd = charPtr;
              if (recordStart != NULL) {
                if (maxRecordLength < (recordEnd-recordStart)) {
                  maxRecordLength =  (recordEnd-recordStart);
                }
              }

            }
            charPtr++;
            bufferBytes--;
          }
          /* we have the last record, now step forward to TAB_GSC_BIN_INDEX */
          if ((recordStart != NULL) && (recordEnd != NULL)) {
            *recordEnd = 0;
            charPtr = recordStart;
            charVal = *charPtr;
            tabCount = 0;
            while(charVal != 0) {
              if (charVal == '\t') {
                tabCount++;
                if (tabCount == TAB_GSC_BIN_INDEX) {
                  charPtr++;
                  nvals = sscanf(charPtr,"%d",&pRunStats->curGscBinIndex);
                  if (nvals != 1) {
                    printf("ERROR: failed to decode curGscBinIndex for rank %d %s\n",pRunStats->rank,temp_name);
                  }
                  break;
                }
              }
              charPtr++;
              charVal = *charPtr;
            }

          }
            
        

        }
        close(fd);
      }

    }

  }
	if (log_handle != NULL) {
		fclose(log_handle);
	}
	if (pParseCommon->out_handle != NULL) {
		fclose(pParseCommon->out_handle);
	}
	time(&curTime);
	curTime -= pParseCommon->startTime;

  if (pParseCommon->verbose) {
    for (runIndex = 0; runIndex < pParseCommon->totalRuns; runIndex++) {
      pRunStats = &pParseCommon->runTable[runIndex];
      if (( pRunStats->endGscBinIndex-pRunStats->startGscBinIndex) == 0) {
        pctdone = 0;
      } else {
        pctdone = (100.0*(pRunStats->curGscBinIndex-pRunStats->startGscBinIndex))/(1.0*(pRunStats->endGscBinIndex-pRunStats->startGscBinIndex));
      }
      finishedBins += pRunStats->curGscBinIndex-pRunStats->startGscBinIndex;
      totalBins += pRunStats->endGscBinIndex-pRunStats->startGscBinIndex;
      printf("rank %2d startGscBinIndex %9d curGscBinIndex %9d endGscBinIndex %9d total bins %9d (%5.1f%) done %d\n",
             pRunStats->rank,
             pRunStats->startGscBinIndex,
             pRunStats->curGscBinIndex,
             pRunStats->endGscBinIndex,
             pRunStats->endGscBinIndex-pRunStats->startGscBinIndex,
             pctdone,
             pRunStats->doneFlag);
    }
    if (totalBins == 0) {
      pctdone = 0;
    } else {
       pctdone = (100.0*finishedBins)/(1.0*totalBins);
    }
    printf("Total bins %d finished bins %d (%5.1f%)\n",totalBins,finishedBins,pctdone);
  }
	printf("Execution Time: %d seconds. maxLineLen %d total runs %d finished: %d maxRecordLength %d\n",
				 curTime,
				 pParseCommon->maxLineLen,
         pParseCommon->totalRuns,
         doneCount,
         maxRecordLength);



	return(EXIT_SUCCESS);
}
