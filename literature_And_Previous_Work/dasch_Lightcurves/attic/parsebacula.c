// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* parsebacula.c
 *
 *  This program examines all of the bacula logs and obtains the usage history of all the tapes.  It is used to select tapes for verification.
 *
 *  Usage:  parsebacula -l  <list of logs>  test with run_932.log for i24330
 *
 *           find /home/scanner/backup -name "bacula.log" | sort -u > ~/junk/bacula.list
 *           parsebacula  -v  -l /home/scanner/junk/bacula.list -o /home/scanner/junk/los.db
 *
 * example 000119L5

./2014_12_11/bacula.log  label command  JobId 318 write 320 aborted 321 verify
./2014_12_18/bacula.log  job 330 write  JobId 336 verify 
./2017_01_26/bacula.log  job 693 restore

 parsebacula  -v  -l /home/scanner/junk/bacula.test -o /home/scanner/junk/los.db

 *
 *
 * cc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm parsebacula.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib/mysql  -l mysqlclient -o parsebacula
 * 
 *
 * May  3, 2019 Edward J. Los - Initial version  
 */   


#include <math.h>
#include <string.h>
#include "table.h"
#include "time.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER     1024
#define MAX_TIME_LENGTH 28
#define MAX_SOLUTION_NUMBER 10
#define MAX_PLATE_COUNT 500000
#define MAX_JOBS 2000
/* #define DEBUG_PARSE 1 */

#if 0
#define DEBUG_SERIES "ai"
#define DEBUG_PLATENUMBER 18190
#endif

typedef struct _jobtype {
  const char* jobname;
  const char* briefname;
} JOBTYPE,*PJOBTYPE;

/* Table of all jobs */
typedef struct _jobs {
  int jobtype;
  char date[MAX_DATE_STRING];
  int jobid;
  int verifyjobid;
} JOB,*PJOB;

#define MAX_VOLUME_NAME 10
#define MAX_TAPES 1000
#define MAX_JOBS_PER_TAPE 20
/* Table of all tapes */
typedef struct _tape {
  char volume[MAX_VOLUME_NAME];
  int numJobs;
  int jobid[MAX_JOBS_PER_TAPE];
} TAPE,*PTAPE;


typedef struct _pparsecommon {
	FILE *out_handle;
	time_t startTime;
	time_t totalStepTime;
	time_t totalTime;
	time_t totalChargeTime; 
	time_t totalTimestampTime;
	int verbose1; /* high level errors */
  int verbose2;
  int partialFlag;
	int maxLineLen;
  int errorCount;
  int totalJobs;
	PJOB jobTable;
  int totalTapes;
  PTAPE tapeTable;
  int maxJobsPerTape;
} PARSECOMMON,*PPARSECOMMON;

JOBTYPE jobtypetable[] = {
  {"NULL","NL"},
  {"BackupClientPrimary","BP"},
  {"BackupClientSecondary","BS"},
  {"VerifyPrimary","VP"},
  {"VerifySecondary","VS"},
  {"BackupMosaics","BM"},
  {"VerifyMosaics","VM"},
  {"BackupDatabaseA","BA"},
  {"BackupDatabaseB","BB"},
  {"RestoreFiles","RS"},
  {"BackupJpeg","BJ"}
};
int numJobTypes = sizeof(jobtypetable)/sizeof(JOBTYPE);


void GetJobDate(PPARSECOMMON pParseCommon,char* inBuffer,int line,int lineCount,char * log_name,void *date)
{
  char *charPtr;
  char *datePtr = (char *)date;
  charPtr = strstr(inBuffer,".20");
  if (charPtr == NULL) {
    printf("ERROR: no expected date stamp %s in line %d of %s\n",inBuffer,lineCount,log_name);
    pParseCommon->errorCount++;
    return;
  }
  charPtr++;
  strncpy(datePtr,charPtr,MAX_DATE_STRING-1);
  datePtr[10] = ' ';
  datePtr[13] = ':';
  datePtr[16] = ':';
  datePtr[19] = 0;
  datePtr[MAX_DATE_STRING-1] = 0;
 
#if 1
  printf("line %d GetJobDate '%s' in %s\n",line,datePtr,inBuffer);
#endif


}


int GetJobType(PPARSECOMMON pParseCommon,char* inBuffer,int line,int lineCount,char* log_name) {
  int jobTypeIndex;
  
  for (jobTypeIndex = 1; jobTypeIndex < numJobTypes; jobTypeIndex++) {
    if (strstr(inBuffer,jobtypetable[jobTypeIndex].jobname) != NULL) {
#if 1
      printf("line %d Found JobName %s\n",line,inBuffer);
#endif
      break;
    }
  }
  if (jobTypeIndex >= numJobTypes) {
    if (pParseCommon->verbose1) {
      printf("ERROR: unrecognized jobType %s in line %d of %s\n",inBuffer,lineCount,log_name);
    }
    pParseCommon->errorCount++;
        
  } else {
    return(jobTypeIndex);
  }
  return(-1);
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
  int outputCount;
	FILE *log_handle = NULL;
	int lineLen;
	int lineCount = 0;
	int recStatCount = 0;
	int parseState = 1; /* 1 = idle; 2 = FitsMosaic; 3 astrometry.csh; 4 WCS16; 5 WCS01 */
	char *charPtr;
	char *charPtr2;
  char *charPtr3;
  char *charPtr4;
	char series[MAX_SERIES_STRING];
	int plateNumber;
	int processCount = 0;
	time_t lastTime = 0;
	time_t prevTime = 0;
	double realSeconds;
	int intSeconds;
	int forwardFlag;
	int saveStep;
	int tmpSolutionNumber;
	int catalogNumber;
  JOB tmpJob;
  PJOB pTmpJob = &tmpJob;
  PJOB pCurJob = NULL;
  int tmpjobid;
  PJOB pJob;
  PJOBTYPE pJobType;
  int jobIndex;
  int jobTypeIndex;
  int tapeIndex;
  PTAPE pTape;
  PTAPE pLastTape = NULL;
  int curJobValid = 0;
  int lastTapeValid = 0;

  if (pParseCommon->verbose1) {
    printf("Processing '%s'\n",log_name);
  }
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
		if (lineCount == 138) {
			printf("At line %d %s\n",lineCount,inLine);
		}
#endif
    strcpy(copyLine,inBuffer);
    if (strstr(inBuffer,"JobName:") == inBuffer) {
      /* Got a new job */
      memset(pTmpJob,0,sizeof(JOB));
      jobTypeIndex = GetJobType(pParseCommon,inBuffer,__LINE__,lineCount,log_name);
      if (jobTypeIndex > 0) {
        pTmpJob->jobtype = jobTypeIndex;
      }
    }
#if 1
    if (lineCount == 157) {
      printf("line %d at lineCount %d\n",__LINE__,lineCount);
    }
#endif
    charPtr = strstr(inBuffer,"When:");
    if (charPtr == inBuffer) {
      /* Got the job date */
      if (pTmpJob->jobtype == 0) {
        if (pParseCommon->verbose1) {
          printf("ERROR: dangling date string  %s in line %d of %s\n",inBuffer,lineCount,log_name);
        }
        pParseCommon->errorCount++;
      } else {
        charPtr2 = strstr(charPtr,"20");
        if (charPtr2 == NULL) {
          printf("ERROR: unrecognized date string %s  in line %d of %s\n",inBuffer,lineCount,log_name);
          pParseCommon->errorCount++;
        } else {
          if (strlen(charPtr2) > (MAX_DATE_STRING-1)) {
            printf("ERROR: too long date string %s in line %d of %s\n",inBuffer,strlen(charPtr2),lineCount,log_name);
            pParseCommon->errorCount++;
          } else {
            strcpy(pTmpJob->date,charPtr2);
#if 1
            printf("line %d Found When: %s\n",__LINE__,inBuffer);
#endif
          }
        }
      }
    }
    curJobValid = 0;
    charPtr = strstr(inBuffer,"JobId ");
    if (charPtr != NULL) {
      charPtr += strlen("JobId ");
      nvals = sscanf(charPtr,"%d",&tmpjobid);
      if (nvals != 1) {
        if (pParseCommon->verbose1) {
          printf("ERROR: tmpjobid string %s  in line %d of %s\n",inBuffer,lineCount,log_name);
        }
        pParseCommon->errorCount++;
      } else {
        if ((pCurJob == NULL) || (pCurJob->jobid != tmpjobid)) {
          /* Switched to a new job */
          for (jobIndex = 0; jobIndex < pParseCommon->totalJobs; jobIndex++) {
            pJob = &pParseCommon->jobTable[jobIndex];
            if (pJob->jobid == tmpjobid) {
              break;
            }
          }
          if (jobIndex >= pParseCommon->totalJobs) {
            /* Allocate a new job here */
            pJob = &pParseCommon->jobTable[pParseCommon->totalJobs];
            memset(pJob,0,sizeof(JOB));
            pJob->jobid = tmpjobid;
            pJobType = &jobtypetable[pJob->jobtype];
            if ((charPtr3 = strstr(inBuffer,"Start Backup")) != NULL) {
              jobTypeIndex = GetJobType(pParseCommon,inBuffer,__LINE__,lineCount,log_name);
              if (jobTypeIndex > 0) {
                pJob->jobtype = jobTypeIndex;
              }
              GetJobDate(pParseCommon,charPtr3,__LINE__,lineCount,log_name,(void *)pJob->date);
            } else if ((charPtr3 = strstr(inBuffer,"Start Restore")) != NULL) {
              jobTypeIndex = GetJobType(pParseCommon,inBuffer,__LINE__,lineCount,log_name);
              if (jobTypeIndex > 0) {
                pJob->jobtype = jobTypeIndex;
              }
              GetJobDate(pParseCommon,charPtr3,__LINE__,lineCount,log_name,(void *)pJob->date);

            } else if ((charPtr3 = strstr(inBuffer,"Verifying against")) != NULL) {
              jobTypeIndex = GetJobType(pParseCommon,inBuffer,__LINE__,lineCount,log_name);
              if (jobTypeIndex > 0) {
                pJob->jobtype = jobTypeIndex;
              }
              GetJobDate(pParseCommon,charPtr3,__LINE__,lineCount,log_name,(void *)pJob->date);
            }
            pJobType = &jobtypetable[pJob->jobtype];        
#if 1
            printf("line %d New Job index %d type %s date %s jobId %4d\n",__LINE__,pParseCommon->totalJobs,pJobType->briefname,pJob->date,pJob->jobid);            
#endif
            pParseCommon->totalJobs++;
            pCurJob = pJob;
            curJobValid = 1;
          } else {
            pCurJob = pJob;
            curJobValid = 1;
#if 1 
            printf("line %d switching to jobid %d in %s\n",__LINE__,pCurJob->jobid,inBuffer);
#endif
          }
        } else {
          /* pCurJob is valid */
          if ((pCurJob != NULL) && (pCurJob->jobid == tmpjobid)) {
            curJobValid = 1;
          }
        }
      }
    }

    charPtr = strstr(inBuffer,"JobId=");
    if (charPtr != NULL) {
      /* Got the job id */
      charPtr += strlen("JobId=");
      nvals = sscanf(charPtr,"%d",&tmpjobid);
      if (nvals != 1) {
        printf("ERROR: jobID string %s  in line %d of %s\n",inBuffer,lineCount,log_name);
        pParseCommon->errorCount++;
      } else {
        if (strstr(inBuffer,"Job queued.") == inBuffer) {
          pTmpJob->jobid = tmpjobid;
          if ((pTmpJob->jobtype == 0) || (pTmpJob->date[0] == 0)) {
            if (pParseCommon->verbose1) {
              printf("ERROR: dangling date string  %s in line %d of %s\n",inBuffer,lineCount,log_name);
            }
            pParseCommon->errorCount++;
          } else {
            for (jobIndex = 0; jobIndex < pParseCommon->totalJobs; jobIndex++) {
              pJob = &pParseCommon->jobTable[jobIndex];
              if (pJob->jobid == pTmpJob->jobid) {
#if 1
                printf("line %d Found JobId: %s\n",__LINE__,inBuffer);
#endif
                break;
              }
            }
            if (jobIndex >= pParseCommon->totalJobs) {
              if (pParseCommon->totalJobs >= (MAX_JOBS-1)) {
                printf("ERROR MAX_JOBS exceeded\n");
                exit(-1);
              }
              pJob = &pParseCommon->jobTable[pParseCommon->totalJobs];
              memcpy(pJob,pTmpJob,sizeof(JOB));
              memset(pTmpJob,0,sizeof(JOB));
              pJobType = &jobtypetable[pJob->jobtype];
#if 1
              printf("line %d New Job index %d type %s date %s jobId %4d\n",__LINE__,pParseCommon->totalJobs,pJobType->briefname,pJob->date,pJob->jobid);
            
#endif
              pParseCommon->totalJobs++;
            }
          }
        } else if (strstr(inBuffer,"Verifying against") != NULL) {
          if (pCurJob == NULL) {
            printf("ERROR: null current job  %s in line %d of %s\n",inBuffer,lineCount,log_name);            
            pParseCommon->errorCount++;
          } else {     
            if (pCurJob->verifyjobid == 0) {
              pCurJob->verifyjobid = tmpjobid;
#if 1
              printf("line %d start verify %s  jobId %4d\n",__LINE__,inBuffer,pCurJob->jobid);
            
#endif
            } else {
              if (pParseCommon->verbose1) {
                printf("ERROR: double verify start for string  %s in line %d of %s\n",inBuffer,lineCount,log_name);            
              }
              pParseCommon->errorCount++;

            }
          

          }

        } else if (strstr(inBuffer,"Start Verify") != NULL) {
          pTmpJob->jobid = tmpjobid;
          if (pCurJob == NULL) {
            printf("ERROR: null current job  %s in line %d of %s\n",inBuffer,lineCount,log_name);            
            pParseCommon->errorCount++;
          } else {
            if (pCurJob->jobid != pTmpJob->jobid) {
              printf("ERROR: inconsistent JobId string  %s in line %d of %s\n",inBuffer,lineCount,log_name);            
              pParseCommon->errorCount++;
            }
          }
        } else {
          if (pParseCommon->verbose1) {
            printf("ERROR: dangling JobId string  %s in line %d of %s\n",inBuffer,lineCount,log_name);  
          }          
          pParseCommon->errorCount++;
        }
      }
      
    }
    lastTapeValid = 0;
    if ((curJobValid != 0) && (pCurJob != NULL)) {
      if ((pLastTape != NULL) && (strstr(inBuffer,pLastTape->volume) != NULL)) {
        lastTapeValid = 1;
      } else {
        charPtr = strstr(inBuffer,"L5\"");
        if (charPtr == NULL) {
          charPtr = strstr(inBuffer,"L6\"");
        }
        if (charPtr != NULL) {
          charPtr[3] = 0;
          charPtr -= 6;
          if ((charPtr >= inBuffer) & (strlen(charPtr) < MAX_VOLUME_NAME)) {
            
            for (tapeIndex = 0; tapeIndex < pParseCommon->totalTapes; tapeIndex++) {
              pTape = &pParseCommon->tapeTable[tapeIndex];
              if (strcmp(charPtr,pTape->volume) == 0) {
                break;
              }
            }
            if (tapeIndex >=  pParseCommon->totalTapes) {
              if (pParseCommon->totalTapes >= MAX_TAPES) {
                printf("ERROR: MAX_TAPES exceeded\n");
                exit(-1);
              }
              pTape = &pParseCommon->tapeTable[pParseCommon->totalTapes];
              pParseCommon->totalTapes++;
              memset(pTape,0,sizeof(TAPE));
              strcpy(pTape->volume,charPtr);
            }
            pLastTape = pTape;
            lastTapeValid = 1;
          }
        } 

      }  
      if (lastTapeValid) {
        for (jobIndex = 0; jobIndex < pLastTape->numJobs; jobIndex++) {
          if (pCurJob->jobid == pLastTape->jobid[jobIndex]) {
            break;
          }
        }
        if (jobIndex >=  pLastTape->numJobs) {
          if (pLastTape->numJobs >= MAX_JOBS_PER_TAPE) {
            printf("ERROR: MAX_JOBS_PER_TAPE exceeded\n");
            exit(-1);
          }
          pLastTape->jobid[pLastTape->numJobs] = pCurJob->jobid;
          pLastTape->numJobs++;
          if (pLastTape->numJobs > pParseCommon->maxJobsPerTape) {
            pParseCommon->maxJobsPerTape = pLastTape->numJobs;
          }
        }

      }

    }
  

  }
  if (pParseCommon->verbose1) {

    printf("Finished processing %s lines %d recStatCount %d\n",log_name,lineCount,recStatCount);
  }
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
  PJOBTYPE pJobType;
  PJOB pJob;
  int jobIndex;
  int tapeIndex;
  int tapeJobIndex;
  PTAPE pTape;
  char tmpDate[MAX_DATE_STRING];

	list_name[0] = 0;
	out_name[0] = 0;
	memset(pParseCommon,0,sizeof(PARSECOMMON));
  pParseCommon->jobTable = (PJOB)calloc(MAX_JOBS,sizeof(JOB));
  if (pParseCommon->jobTable == NULL) {
    printf("ERROR: failed to allocate jobTable\n");
    exit(-1);
  }
  pParseCommon->tapeTable = (PTAPE)calloc(MAX_TAPES,sizeof(TAPE));
  if (pParseCommon->tapeTable == NULL) {
    printf("ERROR: failed to allocate tapeTable\n");
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
					pParseCommon->verbose1 = 1;
					pParseCommon->verbose2 = 1;
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
		printf("Usage: parsebacula [-v] [-p] -l <list of log files> -o <outfile>\n");
		printf("       where -v is the verbose flag\n");
    printf("       where -p is for partial runs not containing a mosaic build\n");
    printf("       NOTE: outfile is currently not supported\n");

		return(-1);
	}

	printf("parsebacula of %s %s List Filename %s Output Filename %s partialFlag %d\n",
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
	}
  if (pParseCommon->verbose1) {

    for (jobIndex = 0; jobIndex < pParseCommon->totalJobs; jobIndex++) {
      pJob = &pParseCommon->jobTable[jobIndex];
      pJobType = &jobtypetable[pJob->jobtype];
#if 1
      printf("Job index %d type %s date %s jobId %4d\n",jobIndex,pJobType->briefname,pJob->date,pJob->jobid);
#endif
    }
  }
  for (tapeIndex = 0; tapeIndex < pParseCommon->totalTapes; tapeIndex++) {
    pTape = &pParseCommon->tapeTable[tapeIndex];
    printf("tape '%8s' with %2d jobs",pTape->volume,pTape->numJobs);
    for (tapeJobIndex = 0; tapeJobIndex < pTape->numJobs; tapeJobIndex++) {
#if 0
      printf(" %4d",pTape->jobid[tapeJobIndex]);
#else
      for (jobIndex = 0; jobIndex < pParseCommon->totalJobs; jobIndex++) {
        pJob = &pParseCommon->jobTable[jobIndex];
        if (pJob->jobid == pTape->jobid[tapeJobIndex]) {
          pJobType = &jobtypetable[pJob->jobtype];
          strcpy(tmpDate,pJob->date);
          tmpDate[10] = 0;
          printf(" %4d:%2s %10s",pJob->jobid,pJobType->briefname,tmpDate);
          break;
        }
      }
      if (jobIndex >= pParseCommon->totalJobs) {
        printf(" %4d:NOT IN TABLE ",pJob->jobid,pJobType->briefname,pJob->date);
      }
#endif
    }
    printf("\n");
  }
 



  printf("totalJobs %d totalTapes %d maxJobsPerTape %d\n",pParseCommon->totalJobs,pParseCommon->totalTapes,pParseCommon->maxJobsPerTape);



	printf("Execution Time: %lld seconds. maxLineLen %d errors %d\n",
				 curTime,
				 pParseCommon->maxLineLen,
         pParseCommon->errorCount);

  if (pParseCommon->jobTable != NULL) {
    free(pParseCommon->jobTable);
  }
  if (pParseCommon->tapeTable != NULL) {
    free(pParseCommon->tapeTable);
  }

	return(EXIT_SUCCESS);
}
