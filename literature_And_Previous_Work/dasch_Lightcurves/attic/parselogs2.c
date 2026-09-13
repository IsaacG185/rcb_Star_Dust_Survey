// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* parselogs2.c
 *
 *  This program checks execution from a single log file
 *
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm parselogs2.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib/mysql  -l mysqlclient -o parselogs2
 * 
      parselogs2 -l /home/scanner/junk/run_1095.log -o /home/scanner/Pipeline/losparse.log -g /home/scanner/Pipeline/losgoodlines.log -b /home/scanner/Pipeline/losbadlines.log -i /home/scanner/Pipeline/losidlelines.log  -r /home/scanner/Pipeline/losretry.list -d /home/scanner/Pipeline/losdone.list -a /home/scanner/Pipeline/losbadapass.list

 *
 * Feb  1, 2017 Edward J. Los - Initial version
 * Feb  9, 2017 Edward J. Los - Add 'a' option to print a list of plates showing no APASS photometry.
 * Mar 17, 2017 Edward J  Los - Skip parsing problems for calibration files.
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
#define MAX_PLATE_COUNT 50000
#define STATE_IDLE 0
#define STATE_ACTIVE 1
#define STATE_COMPLETE 2
#define MAX_PLATE_NAME 25
#define MAX_PHASE 40
#define MAX_ERROR_STRING 80

typedef struct _errortype {
  int index;
  int activecount[MAX_PHASE];
  int idlecount[MAX_PHASE];
  char errorString[MAX_ERROR_STRING];

} ERRORTYPE,*PERRORTYPE;

typedef struct _errortypex {
  int index;
  const char *errorString;
} ERRORTYPEX,*PERRORTYPEX;

#define ERROR_TYPE_PROCESSES  0
#define ERROR_TYPE_RESOURCE   1
#define ERROR_TYPE_SORT       2
#define ERROR_TYPE_CANTOPEN   3
#define ERROR_TYPE_CANTREAD   4
#define ERROR_TYPE_NODEFECT   5
#define ERROR_TYPE_TAWKWRITE  6
#define ERROR_TYPE_NOSEXTRAC  7
#define ERROR_TYPE_NOBACKGRN  8
#define ERROR_TYPE_WRONGSIZE  9
#define ERROR_TYPE_TERMINATE 10
#define ERROR_TYPE_NOTSORTED 11
#define MAX_ERROR_TYPE 12

ERRORTYPEX errorTypeTableX[MAX_ERROR_TYPE] = {
  {ERROR_TYPE_PROCESSES,"No more processes"},
  {ERROR_TYPE_RESOURCE ,"Resource temporarily unavailable"},
  {ERROR_TYPE_SORT     ,"sort exits abnormally"},
  {ERROR_TYPE_CANTOPEN ,"column: can't open input file"},
  {ERROR_TYPE_CANTREAD ,"can't read table header from input"},
  {ERROR_TYPE_NODEFECT ,"No defect file found"},
  {ERROR_TYPE_TAWKWRITE,"tawk: write failure (Broken pipe)"},
  {ERROR_TYPE_NOSEXTRAC,"search_close found no sextractor file for"},
  {ERROR_TYPE_NOBACKGRN,"No background file found"},
  {ERROR_TYPE_WRONGSIZE,"ERROR: wrong file sizes for"},
  {ERROR_TYPE_TERMINATE,"Terminated"},
  {ERROR_TYPE_NOTSORTED,"ERROR: Input sextractor file is not sorted by declination"},
#if 0
  {ERROR_TYPE_HEADCANT ,"head: cannot open"},               /* Usually allobjects file */
  {ERROR_TYPE_DOESNOTEX," does not exist"},                 /* local calibration file   
  {ERROR_TYPE_ERRORBACK,"ERROR: Failed to read table for"}, /* Should be the background file */
#endif
};



int numErrorTypes = sizeof(errorTypeTableX)/sizeof(ERRORTYPEX);

typedef struct _platestats {
  char plateName[MAX_PLATE_NAME];
	char series[MAX_SERIES_STRING];
	int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int index;
  int state;
  int calibrationjob;
  int astrometrywcsjob;
  int apassvalid;
  int anyapassvalid;
  int gscvalid; 
  int errorTypeTable[MAX_ERROR_TYPE+1];
  int errorTypeTotal[MAX_PHASE];
  int startLine[MAX_PHASE];
  int endLine[MAX_PHASE];
} PLATESTATS,*PPLATESTATS;

PPLATESTATS GetNextLine(PPLATESTATS pPlateTable,int numPlates,int lineCounter) {
  int plateIndex;
  int phaseCounter;
  PPLATESTATS pPlateStats;
  PPLATESTATS pHighestStats = pPlateTable;
  PPLATESTATS pNextStats = NULL;
  int highestLineCounter = 0;
  int nextLineCounter = -1;

  for (plateIndex = 0; plateIndex < numPlates; plateIndex++) {
    pPlateStats = &pPlateTable[plateIndex];
    for (phaseCounter = 0; phaseCounter < MAX_PHASE; phaseCounter++) {
      if (pPlateStats->endLine[phaseCounter] > highestLineCounter) {
        highestLineCounter = pPlateStats->endLine[phaseCounter];
        pHighestStats = pPlateStats;
      }
      if ((pPlateStats->startLine[phaseCounter] > 0) &&
          (pPlateStats->startLine[phaseCounter] > lineCounter)) {
        if (nextLineCounter < 0) {
          nextLineCounter = pPlateStats->startLine[phaseCounter];
          pNextStats = pPlateStats;
        } else {
          if (pPlateStats->startLine[phaseCounter] < nextLineCounter) {
            nextLineCounter = pPlateStats->startLine[phaseCounter];
            pNextStats = pPlateStats;
          }
        }
      }
    }
  }
  if (nextLineCounter > 0) {
    return(pNextStats);
  } else {
    return(pHighestStats);
  }

}

int main(int argc,char *argv[])
{
  char *argstr;
  char *inBuffer;
  FILE *out_handle = NULL;
  FILE *good_handle = NULL;
  FILE *idle_handle = NULL;
  FILE *bad_handle = NULL;
  FILE *retry_handle = NULL;
  FILE *done_handle = NULL;
  FILE *apass_handle = NULL;

  char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
  char tempLine[MAX_BUFFER];
  char log_name[MAX_INPUT_NAME];
  char out_name[MAX_INPUT_NAME];
  char good_name[MAX_INPUT_NAME];
  char bad_name[MAX_INPUT_NAME];
  char idle_name[MAX_INPUT_NAME];
  char retry_name[MAX_INPUT_NAME];
  char done_name[MAX_INPUT_NAME];
  char apass_name[MAX_INPUT_NAME];
  FILE *log_handle = NULL;
  int errorFlag = 0;
  char cmdchar;
  int lineLen;
  time_t curTime;
  int stepIndex;
  time_t pTotalStepTime = 0;
  int stepCounter;
  double percentage;
  double percentage2;
  time_t averageTime;
  time_t totalAverageTime = 0;
  int solutionIndex;
  int verbose = 0;
  int maxLineLen = 0;
  time_t startTime;
  int lineCounter = 0;
  int goodCounter = 0;
  int badCounter = 0;
  int idleCounter = 0;
  int retryCounter = 0;
  int doneCounter = 0;
  int badapassCounter = 0;
  char *charPtr2;
  PPLATESTATS pPlateTable;
  PPLATESTATS pPlateStats;
  PPLATESTATS pPlateStats2;
  PPLATESTATS pCurStats = NULL;
  PLATESTATS plateBuffer;
  PPLATESTATS pPlateBuffer = &plateBuffer;
  int numPlates = 0;
  int plateIndex;
  int totalapassvalid = 0;
  int totalgscvalid = 0;
  int totalinvalidapass = 0;
  int errorIndex;
  PERRORTYPE pErrorType;
  PERRORTYPEX pErrorTypeX;
  int outputState = STATE_IDLE;
  int phaseCounter = 0;
  int phaseNumber;

  ERRORTYPE errorTypeTable[MAX_ERROR_TYPE];


  /* Loop through the arguments */
  log_name[0] = 0;
  out_name[0] = 0;
  good_name[0] = 0;
  bad_name[0] = 0;
  idle_name[0] = 0;
  retry_name[0] = 0;
  done_name[0] = 0;
  apass_name[0] = 0;




  if (MAX_ERROR_TYPE != numErrorTypes) {
    printf("ERROR in numErrorTypes %d not %d\n",numErrorTypes,MAX_ERROR_TYPE);
    exit(-1);
  }

  memset(errorTypeTable,0,sizeof(errorTypeTable));
  for (errorIndex = 0; errorIndex < numErrorTypes; errorIndex++) {
    pErrorTypeX = &errorTypeTableX[errorIndex];
    pErrorType = &errorTypeTable[errorIndex];
    pErrorType->index = pErrorTypeX->index;
    if (strlen(pErrorTypeX->errorString) > (MAX_ERROR_STRING-2)) {
      printf("ERROR: MAX_ERROR_STRING exceeded for %d\n",strlen(pErrorTypeX->errorString));
      exit(-1);
    }
    strcpy(pErrorType->errorString,pErrorTypeX->errorString);
  }
  for (errorIndex = 0; errorIndex < numErrorTypes; errorIndex++) {
    pErrorType = &errorTypeTable[errorIndex];
    if (pErrorType->index != errorIndex) {
      printf("ERROR: errorTypeTable format error for index %d\n",errorIndex);
      exit(-1);
    }
  }

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

        case 'l': /* log file name */
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

        case 'i': /* idle file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(idle_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 'g': /* good file name */
        case 'G':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(good_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 'b': /* bad file name */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(bad_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 'r': /* retry file name */
        case 'R':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(retry_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 'd': /* done file name */
        case 'D':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(done_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 'a': /* good APASS photometry */
        case 'A':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(apass_name,*++argv,MAX_INPUT_NAME-2);
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

  if (bad_name[0] != 0) {
    bad_handle = fopen(bad_name,"wt");
    if (bad_handle == NULL) {
      errorFlag = 1;
      printf("Could not open bad file %s\n",bad_name);
    } else {
      printf("Opened bad file %s\n",bad_name);
    }
  }

  if (idle_name[0] != 0) {
    idle_handle = fopen(idle_name,"wt");
    if (idle_handle == NULL) {
      errorFlag = 1;
      printf("Could not open idle file %s\n",idle_name);
    } else {
      printf("Opened idle file %s\n",idle_name);
    }

  }

  if (good_name[0] != 0) {
    good_handle = fopen(good_name,"wt");
    if (good_handle == NULL) {
      errorFlag = 1;
      printf("Could not open good file %s\n",good_name);
    } else {
      printf("Opened good file %s\n",good_name);
    }

  }


  if (retry_name[0] != 0) {
    retry_handle = fopen(retry_name,"wt");
    if (retry_handle == NULL) {
      errorFlag = 1;
      printf("Could not open retry file %s\n",retry_name);
    } else {
      printf("Opened retry file %s\n",retry_name);
    }

  }

  if (done_name[0] != 0) {
    done_handle = fopen(done_name,"wt");
    if (done_handle == NULL) {
      errorFlag = 1;
      printf("Could not open done file %s\n",done_name);
    } else {
      printf("Opened done file %s\n",done_name);
    }

  }

  if (apass_name[0] != 0) {
    apass_handle = fopen(apass_name,"wt");
    if (apass_handle == NULL) {
      errorFlag = 1;
      printf("Could not open apass file %s\n",apass_name);
    } else {
      printf("Opened apass file %s\n",apass_name);
    }

  }




  if (errorFlag) {
    printf("Usage: parselogs2 [-v] -l <log file> -o <outfile>\n");
    printf("       where -v is the verbose flag\n");

    return(-1);
  }

  printf("parselogs2 of %s %s Log Filename %s Output Filename %s\n",
         __DATE__,__TIME__,log_name,out_name);
  pPlateTable = (PPLATESTATS)calloc(MAX_PLATE_COUNT,sizeof(PLATESTATS));
  if (pPlateTable == NULL) {
    printf("allocation error in line %d\n",__LINE__);
    exit(-1);
  }
    


  time(&startTime);

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
    if (lineLen > maxLineLen) {
      maxLineLen = lineLen;
    }
    strcpy(copyLine,inBuffer);
    if (strstr(inBuffer,"Phase") == inBuffer) {
      printf("Phase %d lineCounter %d %s\n",phaseCounter,lineCounter,inBuffer);
      phaseCounter++;
      if (phaseCounter >= MAX_PHASE) {
        printf("MAX_PHASE exceeded in lineCounter %d line %d\n",lineCounter,__LINE__);
        exit(-1);
      }
      for (plateIndex = 0; plateIndex < numPlates; plateIndex++) {
        pPlateStats = &pPlateTable[plateIndex];
        pPlateStats->calibrationjob = 0;
        pPlateStats->astrometrywcsjob = 0;
        pPlateStats->state = STATE_IDLE;
        pPlateStats->gscvalid = 0;
        pPlateStats->apassvalid = 0;
        for (errorIndex = 0; errorIndex < MAX_ERROR_TYPE; errorIndex++) {
          pPlateStats->errorTypeTable[errorIndex] = 0;
        }
      }
    }

    if (strstr(inBuffer,"CHILD") == inBuffer) {
      strcpy(tempLine,inBuffer);
      charPtr2 = strtok(tempLine," ");
      if (charPtr2 == NULL) {
        printf("Null pointer in lineCounter %d line %d\n",lineCounter,__LINE__);
        exit(-1);
      }
      charPtr2 = strtok(NULL," ");
      if (charPtr2 == NULL) {
        printf("Null pointer in lineCounter %d line %d\n",lineCounter,__LINE__);
        exit(-1);
      }
      if (strlen(charPtr2) > (MAX_PLATE_NAME-1)) {
        printf("plate name length exceeded in lineCounter %d line %d\n",lineCounter,__LINE__);
        exit(-1);
      }
      memset(pPlateBuffer,0,sizeof(PLATESTATS));
      strcpy(pPlateBuffer->plateName,charPtr2);
  
      if (strstr(inBuffer,"AstrometryWCS") != NULL) {
        pPlateBuffer->astrometrywcsjob = 1;
      }
      if ((strstr(inBuffer,"FitsLinearity") != NULL) ||
          (strstr(inBuffer,"FitsMedian") != NULL)) {
        pPlateBuffer->calibrationjob = 1;
      }
      if (ParseFilename(pPlateBuffer->plateName,pPlateBuffer->series,&pPlateBuffer->plateNumber,&pPlateBuffer->mosaicNumber,&pPlateBuffer->binning,&pPlateBuffer->rotation) == 0) {
        if (pCurStats == NULL) {
          printf("WARNING parse error in lineCounter %d line %d for %s\n",lineCounter,__LINE__,pPlateBuffer->plateName);
#if 0
          exit(-1);
#endif
        } else if (strcmp(pCurStats->plateName,pPlateBuffer->plateName) != 0) {
          if ((strcmp(pPlateBuffer->plateName,"Exited") == 0) &&
              (pCurStats->astrometrywcsjob != 0)) {
            pPlateBuffer->astrometrywcsjob = 1;
            strcpy(pPlateBuffer->plateName,pCurStats->plateName);
            inBuffer[7] = 'X'; /* Change "Exited" to "XXited" */
          } else {
            printf("parse error in lineCounter %d line %d\n",lineCounter,__LINE__);
            exit(-1);
          }
        } else if (pCurStats->calibrationjob != 0) {
          pPlateBuffer->calibrationjob = 1;
        } else if (pCurStats->astrometrywcsjob != 0) {
          pPlateBuffer->astrometrywcsjob = 1;
          inBuffer[7] = 'X'; /* Change "Exited" to "XXited" */
        } else {
          printf("parse error in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);
        }
      }
      
      for (plateIndex = 0; plateIndex < numPlates; plateIndex++) {
        pPlateStats = &pPlateTable[plateIndex];
        if (strcmp(pPlateStats->plateName,pPlateBuffer->plateName) == 0) {
          break;
        }
      }
      if (plateIndex >= numPlates) {
        if (numPlates > (MAX_PLATE_COUNT-2)) {
          printf("allocation error in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);
        }        
        pPlateStats = &pPlateTable[numPlates];
        memcpy(pPlateStats,pPlateBuffer,sizeof(PLATESTATS));
        pPlateStats->index = numPlates;
        plateIndex = numPlates;
        numPlates++;
        
      } else {
        if (pPlateBuffer->calibrationjob != 0) {
          pPlateStats->calibrationjob = pPlateBuffer->calibrationjob;
        }
        if (pPlateBuffer->astrometrywcsjob != 0) {
          pPlateStats->astrometrywcsjob = pPlateBuffer->astrometrywcsjob;
        }
      }
        
      if (strstr(inBuffer,"Starting") != 0) {
        if (pPlateStats->state != STATE_IDLE) {
          printf("start in idle state in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
        if (pCurStats != NULL) {
          printf("two plates active in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
        pCurStats = pPlateStats;
        pPlateStats->state = STATE_ACTIVE;
        pPlateStats->startLine[phaseCounter] = lineCounter;
#if 0
        printf("ACTIVE: index %d %s\n",plateIndex,pPlateStats->plateName);
#endif
      } else if (strstr(inBuffer,"Exited") != 0) {
        if (pCurStats == NULL) {
          printf("no plate active in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
        if (pCurStats != pPlateStats) {
          printf("wrong plate active in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
        if (pCurStats->state != STATE_ACTIVE) {
          printf("state not active in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
        if (pCurStats->apassvalid != 0) {
          totalapassvalid++;
        }
        if (pCurStats->gscvalid != 0) {
          totalgscvalid++;
        }
        for (errorIndex = 0; errorIndex < MAX_ERROR_TYPE; errorIndex++) {
          pErrorType = &errorTypeTable[errorIndex];
          if (pCurStats->errorTypeTable[errorIndex] != 0) {
            pCurStats->errorTypeTable[MAX_ERROR_TYPE]++;
            pCurStats->errorTypeTotal[phaseCounter]++;
            pErrorType->activecount[phaseCounter]++;
          }
        }

        if (pCurStats->apassvalid == 0) {
          totalinvalidapass++;
        }

        pPlateStats->endLine[phaseCounter] = lineCounter;
        pCurStats->state = STATE_COMPLETE;
        pCurStats = NULL;
#if 0
        printf("COMPLETE: index %d %s\n",plateIndex,pPlateStats->plateName);
#endif
      }  else if (strstr(inBuffer,"Xxited") != 0) {
        if (pCurStats == NULL) {
          printf("no plate active in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
        if (pCurStats != pPlateStats) {
          printf("wrong plate active in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
        if (pCurStats->astrometrywcsjob == 0) {
          printf("not astrometrywcsjob in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
      } else {
        printf("unrecognized CHILD in lineCounter %d line %d\n",lineCounter,__LINE__);
        exit(-1);
      }        

    } else {
      if (strstr(inBuffer,"Illegal") != NULL) {
        if (pCurStats == NULL) {
          printf("Illegal without plate in lineCounter %d line %d\n",lineCounter,__LINE__);
          exit(-1);          
        }
        if (strstr(inBuffer,"_apass") != NULL) {
          pCurStats->apassvalid++;
          pCurStats->anyapassvalid++;
          if (pCurStats->apassvalid > 1) {
            printf("multiple apass valid in lineCounter %d line %d\n",lineCounter,__LINE__);
            exit(-1);          
          }
        } else {
          pCurStats->gscvalid++;
          if (pCurStats->gscvalid > 1) {
            printf("multiple gsc valid in lineCounter %d line %d\n",lineCounter,__LINE__);
            exit(-1);
          }          
        }
      }
      for (errorIndex = 0; errorIndex < MAX_ERROR_TYPE; errorIndex++) {
        pErrorType = &errorTypeTable[errorIndex];
        if (strstr(inBuffer,pErrorType->errorString) != NULL) {
          if (pCurStats == NULL) {
            pErrorType->idlecount[phaseCounter]++;
          } else {
            pCurStats->errorTypeTable[errorIndex]++;
          }
        }
      }

    }
    lineCounter++;
  }

  for (plateIndex = 0; plateIndex < numPlates; plateIndex++) {
    for (phaseNumber = 0; phaseNumber < phaseCounter; phaseNumber++) {
      pPlateStats = &pPlateTable[plateIndex];
      if ((pPlateStats->startLine[phaseNumber] != 0) &&
          (pPlateStats->endLine[phaseNumber] != 0) &&
          (pPlateStats->startLine[phaseNumber] >= pPlateStats->endLine[phaseNumber])) {
        printf("start and end line error in lineCounter %d line %d\n",lineCounter,__LINE__);
        exit(-1);
      }
    }
  }


  if (log_handle != NULL) {
    fclose(log_handle);
  }
  log_handle = fopen(log_name,"rt");
  if (log_handle == NULL) {
    printf("ERROR: line %d Could not open list file %s\n",__LINE__,log_name);
  }
  lineCounter = 0;
  plateIndex = 0;
  phaseNumber = 0;
  pPlateStats = GetNextLine(pPlateTable,numPlates,lineCounter);
  if (pPlateStats == NULL) {
     printf("ERROR: line %d null next plate\n",__LINE__);
  }
   
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
    if (strstr(inBuffer,"Phase") == inBuffer) {
      printf("Phase %d lineCounter %d %s\n",phaseNumber,lineCounter,inBuffer);
      fprintf(good_handle,"Phase %d lineCounter %d %s\n",phaseNumber,lineCounter,inBuffer);
      fprintf(bad_handle,"Phase %d lineCounter %d %s\n",phaseNumber,lineCounter,inBuffer);
      
      phaseNumber++;
      if (phaseNumber >= MAX_PHASE) {
        printf("MAX_PHASE exceeded in lineCounter %d line %d\n",lineCounter,__LINE__);
        exit(-1);
      }
    }
    if (outputState == STATE_IDLE) {
      fprintf(idle_handle,"%s\n",inBuffer);
      idleCounter++;
    } else {
      if ((pPlateStats->calibrationjob != 0) ||
          (pPlateStats->errorTypeTotal[phaseNumber] == 0)) {
        fprintf(good_handle,"%s\n",inBuffer);
        goodCounter++;
      } else {
        fprintf(bad_handle,"%s\n",inBuffer);
        badCounter++;
      }

    }

    lineCounter++;
    if ((lineCounter >= pPlateStats->startLine[phaseNumber]) && 
        (lineCounter <= pPlateStats->endLine[phaseNumber])) {
      if (outputState == STATE_IDLE) {
        outputState = STATE_ACTIVE;
      }
    } else {
      if (outputState == STATE_ACTIVE) {
        outputState = STATE_IDLE;
        pPlateStats = GetNextLine(pPlateTable,numPlates,lineCounter);
        
      }
    }
  }
        
  
  for (plateIndex = 0; plateIndex < numPlates; plateIndex++) {
    pPlateStats = &pPlateTable[plateIndex];
    if (pPlateStats->calibrationjob != 0) {
      continue;
    }
    fprintf(out_handle,"SELECT series,plateNumber,solutionNumber FROM localbin2 INNER JOIN photseries USING (seriesId) where local_bin_index = 0 and series = '%s' and plateNumber = %d;\n",pPlateStats->series,pPlateStats->plateNumber);
    if (pPlateStats->errorTypeTable[MAX_ERROR_TYPE] != 0) {
      retryCounter++;
      fprintf(retry_handle,"%s\n",pPlateStats->plateName);
    } else {
      doneCounter++;
      fprintf(done_handle,"%s\n",pPlateStats->plateName);
      if (pPlateStats->anyapassvalid == 0) {
        badapassCounter++;
        fprintf(apass_handle,"%s\n",pPlateStats->plateName);
      }
    }
  }


  if (log_handle != NULL) {
    fclose(log_handle);
  }
  if (out_handle != NULL) {
    fclose(out_handle);
  }
  if (good_handle != NULL) {
    fclose(good_handle);
  }
  if (bad_handle != NULL) {
    fclose(bad_handle);
  }
  if (idle_handle != NULL) {
    fclose(idle_handle);
  }
  if (retry_handle != NULL) {
    fclose(retry_handle);
  }
  if (done_handle != NULL) {
    fclose(done_handle);
  }
  if (apass_handle != NULL) {
    fclose(apass_handle);
  }
  time(&curTime);
  curTime -= startTime;

  printf("Execution Time: max line length %d total lines %d good lines %d bad lines %d idle lines %d total plates %d %lld seconds. \n",
         maxLineLen,
         lineCounter,
         goodCounter,
         badCounter,
         idleCounter,
         numPlates,
         curTime);
  printf("phaseCounter %d\n",phaseCounter);
  printf("doneCounter %d retryCounter %d badapassCounter %d\n",doneCounter,retryCounter,badapassCounter);
  printf("Total apass valid %d apass invalid %d\n",totalapassvalid,totalinvalidapass);
  printf("Total gsc valid %d\n",totalgscvalid);
  for (phaseNumber = 0; phaseNumber < phaseCounter; phaseNumber++) { 
    for (errorIndex = 0; errorIndex < MAX_ERROR_TYPE; errorIndex++) {
      pErrorType = &errorTypeTable[errorIndex];
      if ((pErrorType->activecount[phaseNumber] != 0) ||
          (pErrorType->idlecount[phaseNumber] != 0)) {
        printf("Phase %2d Active: %5d Idle %5d for %s\n",phaseNumber,pErrorType->activecount[phaseNumber],pErrorType->idlecount[phaseNumber],pErrorType->errorString);
      }
    }
  }




  return(EXIT_SUCCESS);
}
