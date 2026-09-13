// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* parsescratchlog.c
 *
 *  This program extracts scratch statistics from the filter_scratch log output by parsing lines of the following format:
 *
 *   Sample FalsePositive scratchActual 0 defect 0 scratchDeclared 1 chainLength 1 lengthRatio 0 for rb12499_00 NUMBER 343764 X 8288 Y 10364 m:11: individual stars
 *
 *  Usage:  parsescratch -l <logfile>
 *
 *    parsescratchlog -l los3.log -o los.db
 *
 *
 * cc -ggdb -O2 -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I /dasch/install/include -L /dasch/install/lib -lm parsescratchlog.c pipelineutils.a -ltable -lutil -lwcs -L/usr/lib${lib64}/mysql -l mysqlclient -o parsescratchlog
 * 
 *
 * Dec 26, 2013 Edward J. Los - Initial version  
 * Jan 21, 2013 Edward J. Los - Find the mosaic processing time
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

typedef struct _resultTable {
  int countTotal;
  int countGood;
  double ratio;
  int scratchesFound;
  char *resultDescription;
} RESULTTABLE,*PRESULTTABLE;

RESULTTABLE resultTable[] = {
  {0,0,0.0,0,"Chain Length Points"},
  {0,0,0.0,0,"Length Ratio Points"},
  {0,0,0.0,0,"All points"},
};
#define CHAIN_LENGTH 0
#define LENGTH_RATIO 1
#define ALL_POINTS   2

int resultTableLength = sizeof(resultTable)/sizeof(RESULTTABLE);


int main(int argc,char *argv[])
{
	char *argstr;
	char *inBuffer;
	char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
	char log_name[MAX_INPUT_NAME];
	char out_name[MAX_INPUT_NAME];
	FILE *log_handle = NULL;
	int errorFlag = 0;
	char cmdchar;
	int lineLen;
  int lineNumber = 0;
	time_t curTime;
  time_t startTime;
	int stepIndex;
	time_t pTotalStepTime = 0;
	int stepCounter;
	double percentage;
	double percentage2;
	time_t averageTime;
  time_t totalTime;
	time_t totalAverageTime = 0;
	int solutionIndex;
  int verbose = 0;
  int defectOnly = 0;
  int excludeDefect = 0;
  int reprocessFlag = 0;
  int maxLineLen = 0;
  FILE* out_handle = NULL;
  int scratchActual;
  char *pScratchActual;
  int defect;
  char *pDefect;
  int scratchDeclared;
  char *pScratchDeclared;
  int chainLength;
  char *pChainLength;
  int lengthRatio;
  char *pLengthRatio;
  char *plateString;
  int NUMBER;
  char *pNumber;
  char *charPtr;
  int nvals;
  int goodLines = 0;
  PRESULTTABLE pResult;
  int resultIndex;


	/* Loop through the arguments */
	log_name[0] = 0;
	out_name[0] = 0;


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
					verbose = 1;
					break;

				case 'd':
				case 'D':
					defectOnly = 1;
					break;

				case 'n':
				case 'N':
					excludeDefect = 1;
					break;

        


				case 'r':
				case 'R':
					reprocessFlag = 1;
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
		printf("ERROR: log file name not specified\n");
		errorFlag = 1;
	} else {
		log_handle = fopen(log_name,"rt");
		if (log_handle == NULL) {
			errorFlag = 1;
			printf("Could not open list file %s\n",log_name);
		}
	}

	if (out_name[0] == 0) {
		printf("ERROR: Out file name not specified\n");
		errorFlag = 1;
	} else {
		out_handle = fopen(out_name,"wt");
		if (out_handle == NULL) {
			errorFlag = 1;
			printf("Could not open out file %s\n",out_name);
		}
	}



	if (errorFlag) {
		printf("Usage: parsescratchlog [-v] -l <list of log files> -o <outfile>\n");
		printf("       where -v is the verbose flag\n");
		printf("             -r prints summary statistics only for the reprocessing states\n");
    printf("             -d means include filter_defect defects only \n");
    printf("             -n means exclude filter_defect defects\n");
		return(-1);
	}

	printf("parsescratchlog of %s %s List Filename %s Output Filename %s defectOnly %d excludeDefect %d\n",
				 __DATE__,__TIME__,log_name,out_name,defectOnly,excludeDefect);
 

	time(&startTime);

	while (1) {
		inBuffer = fgets(inLine,MAX_BUFFER,log_handle);
		if (inBuffer == NULL) {
			break;
			
		}
    lineNumber++;
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
    if (lineLen > maxLineLen) {
      maxLineLen = lineLen;
    }
    strcpy(copyLine,inBuffer);
    if (strstr(inBuffer,"Sample ") != NULL) {
      pScratchActual = strstr(inBuffer,"scratchActual");
      pDefect = strstr(inBuffer,"defect");
      pScratchDeclared = strstr(inBuffer,"scratchDeclared");
      pChainLength = strstr(inBuffer,"chainLength");
      pLengthRatio = strstr(inBuffer,"lengthRatio");
      plateString = strstr(inBuffer,"for");
      pNumber = strstr(inBuffer,"NUMBER");

      if ((pScratchActual != NULL) && 
          (pDefect != NULL) && 
          (pScratchDeclared != NULL) && 
          (pChainLength != NULL) && 
          (pLengthRatio != NULL) && 
          (plateString != NULL) && 
          (pNumber != NULL)) {
        pScratchActual += 14;
        *(pScratchActual+1) = 0;
        pDefect += 7;
        *(pDefect+1) = 0;
        pScratchDeclared += 16;
        *(pScratchDeclared+1) = 0;
        pChainLength += 12;
        *(pChainLength+1) = 0;
        pLengthRatio += 12;
        *(pLengthRatio+1) = 0;
        plateString += 4;
        charPtr = strstr(plateString," ");
        if (charPtr != NULL) {
          *charPtr = 0;
        }
        pNumber += 7;
        nvals = sscanf(pNumber,"%d",&NUMBER);
        if (nvals == 1) {
          goodLines++;
          if (*pScratchActual == '0') {
            scratchActual = 0;
          } else {
            scratchActual = 1;
          }

          if (*pDefect == '0') {
            defect = 0;
          } else {
            defect = 1;
          }

          if (*pScratchDeclared == '0') {
            scratchDeclared = 0;
          } else {
            scratchDeclared = 1;
          }

          if (*pChainLength == '0') {
            chainLength = 0;
          } else {
            chainLength = 1;
          }

          if (*pLengthRatio == '0') {
            lengthRatio = 0;
          } else {
            lengthRatio = 1;
          }


#if 0
          printf("\n%s\n",copyLine);
          printf("\nline %d, scratchActual %s defect %s scratchDeclared %s chainLength %s lengthRatio %s plate %s NUMBER %6d\n",
                 lineNumber,
                 pScratchActual,
                 pDefect,
                 pScratchDeclared,
                 pChainLength,
                 pLengthRatio,
                 plateString,
                 NUMBER);
#endif
         if ((excludeDefect == 1) && (defect == 1)) {
           continue;
         }
         if ((defectOnly == 1) && (defect == 0)) {
           continue;
         }
#if 0
          printf("line %d, scratchActual %d defect %d scratchDeclared %d chainLength %d lengthRatio %d plate %s NUMBER %6d\n",
                 lineNumber,
                 scratchActual,
                 defect,
                 scratchDeclared,
                 chainLength,
                 lengthRatio,
                 plateString,
                 NUMBER);
#endif        
          pResult = &resultTable[CHAIN_LENGTH];
          pResult->countTotal++;
          if ((scratchActual ^ chainLength) == 0) {
            pResult->countGood++;
            if (scratchActual != 0) {
              pResult->scratchesFound++;
            }
          }
          pResult = &resultTable[LENGTH_RATIO];
          pResult->countTotal++;
          if ((scratchActual ^ lengthRatio) == 0) {
            pResult->countGood++;
            if (scratchActual != 0) {
              pResult->scratchesFound++;
            }
          }
          pResult = &resultTable[ALL_POINTS];
          pResult->countTotal++;
          if (scratchActual != 0) {
            pResult->scratchesFound++;
          }
          if ((scratchActual ^ scratchDeclared) == 0) {
            pResult->countGood++;
          }
        }
      }
    }
	}
  for (resultIndex = 0; resultIndex < resultTableLength; resultIndex++) {
    pResult = &resultTable[resultIndex];
    if (pResult->countTotal != 0) {
      pResult->ratio = (1.0*pResult->countGood)/(1.0*pResult->countTotal);
    } 
    printf("Count %3d  Good %3d Ratio %.2f scratchesFound %d %s\n",
           pResult->countTotal,
           pResult->countGood,
           pResult->ratio,
           pResult->scratchesFound,
           pResult->resultDescription);
  }

	if (log_handle != NULL) {
		fclose(log_handle);
	}
	if (out_handle != NULL) {
		fclose(out_handle);
	}
	time(&curTime);
	curTime -= startTime;



	printf("Execution Time: %lld seconds. maxLineLen %d goodLines %d totalTime %lld totalAverageTime %lld\n",
				 curTime,
				 maxLineLen,
         goodLines,
				 totalTime,
				 totalAverageTime);





	return(EXIT_SUCCESS);
}
