// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* web_queue.c
 *
 *  This program implements a web server daemon which periodically looks for queue requests 
 *
 *
 *  
 *  gcc -ggdb -O0 -I/dasch/install/include -D_FILE_OFFSET_BITS=64      web_queue.c -lm -pthread -o web_queue 
 *
 *
 *  Usage web_queue -d <root name>
 *
 * Feb 21, 2014 Edward J. Los - Initial version
 * Feb 28, 2014 Edward J. Los - Add a status printout for all pending requests.
 * Mar  3, 2014 Edward J. Los - For daschprivate, allow multiple active operations
 * Mar 21, 2014 Edward J. Los - Correct the -b qualifier to web_query to add the temporary directory
 * Apr 28, 2014 Edward J. Los - Add a retry on the syncronization lock
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 *
 *  /dasch/Pipeline/web_queue -d /home/scanner/web/dasch/tmp/Yxk9zE -v
 *
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <dirent.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/wait.h>
#include <errno.h>
#include "pipelineutils.h"
#include "daschunistd.h"
#define MAX_BUFFER 500

/* The  MAX_QUEUE_ENTRIES * (sizeof(REQUESTQUEUE)+ MAX_REQUEST) is the amount of memory used.  
 * The websites have 1GB of memory and at least 1GB of swap space 
 * In a denial of service attack, the web_queue program should not use more than 2 MB ot 0.2% of memory or about 300 entries
 */
#define MAX_QUEUE_ENTRIES 300 /* Maximum queue entries. */
#define MAX_REQUEST 50000
#define MAX_PATH   512
#define CURRENT_REQUEST_VERSION 4 /* This must match lightcurve_queue.php request version */
#define LIST_REALLOC_INCREMENT 10
/* #define CHILD_DEBUG 1 */

typedef struct _requestqueue {
  long long whenqueued;  /* time the request was queued */
  pid_t childPid;        /* Child PID */
  int requestDirectoryFound; /* Request directory still exists */
  int requestHandled;    /* Request already handled */
  int requestInProgress; /* Request is in progress */
  int requestError;      /* Error decoding the request */
  int haveRequest;       /* request.txt present */
  int haveResult;        /* result.txt present */
  int requestVersion;    /* Request version number */
  int nmin;              /* Minumum number of stars */
  int box;               /* Arcsecond radius */
  int objectCount;       /* Number of objects */
  time_t requestSeconds; /* Time of the request.txt file */
  long requestNsec;      /* Time of the request.txt file */
  time_t resultSeconds; /* Time of the result.txt file */
  long resultNsec;      /* Time of the result.txt file */
  
  char subdir[MAX_PATH]; /* Request subdirectory name */
  char subDirectoryName[MAX_PATH]; /* Full subdirectory path name */
  char source[MAX_BUFFER]; /* catalog name */
  char researcher[MAX_BUFFER]; /* researcher name */
  char frameformat[MAX_BUFFER]; /* Frame format */
  char binaries[MAX_BUFFER]; /* Directory of the pipeline web routines */
  char lsbin[MAX_BUFFER];    /* Directory of the system 'ls' command */
  char tsocks[MAX_BUFFER];  /* Contains 'tsocks' if indirect web access is necessary */
  char tmpsubdir[MAX_BUFFER+MAX_PATH]; /* Contains the subdirectory (should be the same as subdir) */
  char resultFileName[MAX_PATH];
  char remoteaddr[MAX_BUFFER];
  char *coo;                /* Pointer to the coordinate list */
} REQUESTQUEUE,*PREQUESTQUEUE;

typedef struct _common {
  int verbose;
  int requestCount;
  int busyFlag;
  PREQUESTQUEUE requestTable;
  char errorBuffer[MAX_BUFFER];
} COMMON, *PCOMMON;

/* Deposit a single error message in the result file */
void WriteResultFile(PCOMMON pCommon,PREQUESTQUEUE pRequest) {
  FILE *tmpFileHandle;
  char tmpFileName[MAX_PATH];
  char cmdStr[2*MAX_PATH];
  int result;
  pRequest->requestHandled = 1;
  strcpy(tmpFileName,pRequest->subDirectoryName);
  strcat(tmpFileName,"tmpresult.txt");
  tmpFileHandle = fopen(tmpFileName,"wt");
  if (tmpFileHandle == NULL) {
    if (pCommon->verbose) {
      printf("ERROR: the temporary result file %s could not be opened\n",tmpFileName);
    }
    return;
  }
  fprintf(tmpFileHandle,"%s",pCommon->errorBuffer);
  fclose(tmpFileHandle);
  sprintf(cmdStr,"mv %s %s\n",tmpFileName,pRequest->resultFileName);
  result = system(cmdStr);
  if (pCommon->verbose) {
    printf("Result is %d for %s, errorBuffer %s\n",result,cmdStr,pCommon->errorBuffer);
  }
  return;
  
}
int RequestCompare(const void *first, const void *second) 
{
  long long whenqueuedFirst = *((long long *)first);
  long long whenqueuedSecond = *((long long *)second);
  if (whenqueuedFirst > whenqueuedSecond) {
    return(1);
  } else if (whenqueuedFirst < whenqueuedSecond) {
    return(-1);
  } else {
    return(0);
  }
}

void WriteStatusFiles(PCOMMON pCommon,PREQUESTQUEUE pActiveRequest,PREQUESTQUEUE pDeletionCandidate) {
  int requestIndexOwner; /* This is the directory that we are going to write into */
  PREQUESTQUEUE pRequestOwner;
  int requestIndexQueue; /* This is the queue entry that we are going to write into the directory */
  PREQUESTQUEUE pRequestQueue;
  FILE *statusFile = NULL;
  char tmpStatusFileName[MAX_PATH];
  char statusFileName[MAX_PATH];
  char cmdStr[3*MAX_PATH];
  char thisRequest[MAX_BUFFER];
  char requestStatus[MAX_BUFFER];
  int result;
  int queueCount = 0;
  double requestTime;
  double resultTime;
  for (requestIndexOwner = pCommon->requestCount-1; requestIndexOwner >= -1; requestIndexOwner --) {
    if (requestIndexOwner == -1) {
      if (pDeletionCandidate != NULL) {
        pRequestOwner = pDeletionCandidate; 
#if 0
        if ((pRequestOwner->requestError == 0) && (pRequestOwner->haveResult == 0)) {
          printf("Good deletion candidate\n");
        }
#endif
      } else {
#if 0
        printf("No deletion candidate\n");
#endif
        continue;
      }
    } else {
      pRequestOwner = &pCommon->requestTable[requestIndexOwner];
    }
#if 1
    if ((pRequestOwner->requestError != 0) ||
        (pRequestOwner->haveResult != 0)) {
      continue;
    }
#else /* Here for debugging */
    if (pRequestOwner->requestError != 0) {
      continue;
    }
#endif
    strcpy(statusFileName,pRequestOwner->subDirectoryName);
    strcat(statusFileName,"status.txt");
    strcpy(tmpStatusFileName,pRequestOwner->subDirectoryName);
    strcat(tmpStatusFileName,"tmpstatus.txt");
    statusFile = fopen(tmpStatusFileName,"wt");
    if (statusFile == NULL) {
      continue;
    }
    fprintf(statusFile,"<h2>Current Queue Status</h2>\n");
    if (pCommon->busyFlag) {
      fprintf(statusFile,"<br /><br />WARNING: The system is currently too busy to list all queued requests.<br />");
    }
    fprintf(statusFile,"<table class=\"bordered\" summary=\"Table of queue entries\">");
    fprintf(statusFile,"<tr><th>Entry Number</th><th>Time Queued</th><th>Coordinate Pairs</th><th>Search Radius (arcsec)</th><th>Owner</th><th>Status</th></tr>");
    queueCount = 0;
    for (requestIndexQueue = pCommon->requestCount-1; requestIndexQueue >= 0; requestIndexQueue --) {
      pRequestQueue = &pCommon->requestTable[requestIndexQueue];
      if (pRequestQueue->requestError != 0) {
        continue;
      }
      queueCount++;
      if (strcmp(pRequestOwner->subdir,pRequestQueue->subdir) == 0) {
        strcpy(thisRequest,"This request.");
      } else {
        if (strcmp(pRequestOwner->remoteaddr,pRequestQueue->remoteaddr) == 0) {
          strcpy(thisRequest,"Request from this client.");
        } else {
          strcpy(thisRequest,"Request from other client.");
        }
      }
      if ((pActiveRequest != NULL) && (strcmp(pActiveRequest->subdir,pRequestQueue->subdir) == 0)) {
        strcpy(requestStatus,"In Progress");
      } else {
        if (pRequestQueue->requestInProgress != 0) {
          strcpy(requestStatus,"In Progress");
        } else if (pRequestQueue->haveResult == 0) {
          strcpy(requestStatus,"Pending");
        } else {
          resultTime = (1.0 * (pRequestQueue->resultSeconds)) + ((1.0*pRequestQueue->resultNsec)/1000000000.0);
          requestTime = (1.0 * (pRequestQueue->requestSeconds)) + ((1.0*pRequestQueue->requestNsec)/1000000000.0);
          sprintf(requestStatus,"Completed in %.1f seconds\n",resultTime-requestTime);
        }
      }
      fprintf(statusFile,"<tr><td>%d</td><td>%s</td><td>%d</td><td>%d</td><td>%s</td><td>%s</td></tr>",queueCount,ctime(&pRequestQueue->requestSeconds),pRequestQueue->objectCount,pRequestQueue->box,thisRequest,requestStatus);
      if (pCommon->verbose) {
        printf("Entry %d Time Queued %s Whose: %s Status %s\n",queueCount,ctime(&pRequestQueue->requestSeconds),thisRequest,requestStatus);
      }
      


    }
    fprintf(statusFile,"</table>");
    if (queueCount == 0) {
      fprintf(statusFile,"<br />No pending queue entries.");
    } 

    fclose(statusFile);
    sprintf(cmdStr,"mv %s %s\n",tmpStatusFileName,statusFileName);
    result = system(cmdStr);
    if (pCommon->verbose) {
      printf("Result is %d for %s\n",result,cmdStr);
    }
    
  }
}


void ChildProcess(PCOMMON pCommon,PREQUESTQUEUE pBestRequest,char *cmdStr)
{
  int result;
  int statResult;
  struct stat statbuf;
  char requestFileName[MAX_PATH];
  char tmpResultFileName[MAX_PATH];
  char tmpdirectory[MAX_PATH];  
  char *slashPtr;
  char tmprootdir[MAX_BUFFER];



  strcpy(requestFileName,pBestRequest->subDirectoryName);
  strcat(requestFileName,"request.txt");
  strcpy(tmpResultFileName,pBestRequest->subDirectoryName);
  strcat(tmpResultFileName,"tmpresult.txt");


  strcpy(tmpdirectory,pBestRequest->subDirectoryName);
          
  slashPtr = strrchr(tmpdirectory,'/');
  if (*(slashPtr+1) == 0) {
    *slashPtr = 0;
    slashPtr = strrchr(tmpdirectory,'/');
  }
  if (slashPtr == NULL) {
    strcpy(tmprootdir,tmpdirectory);
    strcpy(pBestRequest->tmpsubdir,tmpdirectory);
  } else {
    *slashPtr = 0;
    strcpy(tmprootdir,tmpdirectory);
    slashPtr++;
    strcpy(pBestRequest->tmpsubdir,slashPtr);
    /* Now put back the subdirectory */
    strcat(tmpdirectory,"/");
    strcat(tmpdirectory,pBestRequest->tmpsubdir);

  }
  

  sprintf(cmdStr,"chmod 666 %s/web_queue%s.log\n",tmprootdir,pBestRequest->tmpsubdir);
  printf("cmdStr %s",cmdStr);
  system(cmdStr);
  printf("Done executing the command\n");
  sprintf(cmdStr,"%s %s/web_query -f %s -n %d -q %s -r %d -b %s -d %s -l %s -e %s  '%s' >& %s\n",
          pBestRequest->tsocks,pBestRequest->binaries,"guest",pBestRequest->nmin,pBestRequest->source,pBestRequest->box,tmpdirectory,pBestRequest->tmpsubdir,pBestRequest->lsbin,pBestRequest->binaries,pBestRequest->coo,tmpResultFileName);
  printf("cmdStr %s",cmdStr);
  WriteStatusFiles(pCommon,pBestRequest,NULL);
  system(cmdStr);
  WriteStatusFiles(pCommon,NULL,NULL);
  printf("Done executing the command\n");
  statResult = stat(tmpResultFileName,&statbuf);
  if (statResult == 0) {
    /* We have a result */
    printf("Found %s\n",tmpResultFileName);
    sprintf(cmdStr,"mv %s %s\n",tmpResultFileName,pBestRequest->resultFileName);
  } else {
    printf("Did not find %s\n",tmpResultFileName);
    /* We are in trouble here, just copy the request file to the result file */
    sprintf(cmdStr,"mv %s %s\n",requestFileName,pBestRequest->resultFileName);
  }
  result = system(cmdStr);
  if (pCommon->verbose) {
    printf("Result is %d for %s\n",result,cmdStr);
  }          
#ifndef CHILD_DEBUG
  free(cmdStr);
  exit(0);
#endif /* CHILD_DEBUG */


}

int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  char charVal;
  char *slashPos;
  char *lastSlashPos;
  char *charPtr;
  time_t startTime;
  time_t curTime;
  char timestr[100];
  struct tm *ptr;
  int lineLen;
  int nlines = 0;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  int errorFlag = 0;
  char *tmpdir = NULL;
  char masterLog[MAX_PATH];
  char keepaliveFile[MAX_PATH];
  char tmpDirectoryRoot[MAX_PATH];
  char requestFileName[MAX_PATH];
  char tmpResultFileName[MAX_PATH];
  char tmpdirectory[MAX_PATH];  
  char *slashPtr;
  FILE *masterLogHandle = NULL;
  FILE *keepaliveHandle = NULL;
  int iteration = 0;
  int result;
  int workDone;
  struct dirent **directoryList = NULL;
  int directoryCount;
  int directoryIndex;
  int statResult;
  struct stat statbuf;
  FILE *requestHandle = NULL;
  int lineCounter;
  char lsbin[MAX_BUFFER];
  char tmprootdir[MAX_BUFFER];
  char *tmpcoo;
  char *cmdStr = NULL;
  COMMON commonBlock;
  PCOMMON pCommon = &commonBlock;
  REQUESTQUEUE newRequest;
  PREQUESTQUEUE pNewRequest = &newRequest;
  PREQUESTQUEUE tmpRequestTable = NULL;
  PREQUESTQUEUE pBestRequest = NULL;
  PREQUESTQUEUE pRequest;
  PREQUESTQUEUE pRequest2;
  PREQUESTQUEUE pRequest3;
  int requestAlloc = 0;
  int requestIndex;
  int errorRequestIndex;  /* Pointer to a request which has an error and can be deleted */
  int finishedRequestIndex; /* Pointer to a finished request which can be deleted */
  int latestRequestIndex; /* Pointer to the last issued valid request */
  int requestSkip = 0;
  int activeJobs = 0;
  pid_t pid;
  int stat_val;
  int requestSize;
  int lockretrycount = 3;
  int lockretrydelay;
  memset(pNewRequest,0,sizeof(REQUESTQUEUE));

  memset(pCommon,0,sizeof (COMMON));
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      fprintf(stderr,"ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'd': /* temporary subdirectory name */
        case 'D':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            tmpdir = *++argv;
          }
          break;

        case 'v': /* Verbose */
        case 'V':
          pCommon->verbose = 1;
          break;

        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  if (tmpdir == NULL) {
    printf("ERROR: no subdirectory specified for web_queue\n");
    errorFlag = 1;
  }
  if (errorFlag == 1) {
    printf("Usage: web_queue -d <temporary directory>\n");
    exit(-1);
  }
  time(&startTime);
  ptr = localtime(&startTime);
  strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
  printf("web_queue of %s %s starting at %s with directory %s sizeof(REQUESTQUEUE) is %d\n",__TIME__,__DATE__,timestr,tmpdir,sizeof(REQUESTQUEUE));
  
  /* The first order of business is to see if any other instances of web_queue are running by attempting to get write access to the log file */

  strcpy(masterLog,tmpdir);
  charPtr = strstr(masterLog,"/tmp/");
  if (charPtr == NULL) {
    printf("ERROR: failed to find 'tmp' in %s\n",masterLog);
    exit(-1);
  }
  charPtr += strlen("/tmp/");
  *charPtr = 0;
  strcpy(tmpDirectoryRoot,masterLog);
  strcat(masterLog,"web_queue.log");
  masterLogHandle = fopen(masterLog,"a+t");
  if (masterLogHandle == NULL) {
    printf("ERROR: failed to open %s\n",masterLog);
    exit(-1);
  }
  strcpy(keepaliveFile,tmpDirectoryRoot);
  strcat(keepaliveFile,"keepalive.txt");
 

  /* Now attempt to place an exclusive lock on this file */
  while (1) {
    time(&curTime);
    result = flock(masterLogHandle->_fileno,LOCK_EX | LOCK_NB);
    if (result != 0) {
      lockretrydelay = (rand()/(RAND_MAX/10)) + 1;
      printf("ERROR: Failed to lock %s errno: %d %s at time %d\n",masterLog,errno,strerror(errno),curTime);
      printf("At line %d RAND_MAX %d lockretrycount %d delay %d pid %d getppid %d\n",__LINE__,RAND_MAX,lockretrycount,lockretrydelay,getpid(),getppid());
      if (--lockretrycount <= 0) {
        exit(-1);
      }
      sleep(lockretrydelay);
      

    } else {
      break;
    }
  }
  result = chmod(masterLog,S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH);
  if (result != 0) {
    printf("ERROR: Failed to change protection on %s errno: %d %s\n",masterLog,errno,strerror(errno));
    fprintf(masterLogHandle,"ERROR: Failed to change protection on %s errno: %d %s\n",masterLog,errno,strerror(errno));
  }

  fprintf(masterLogHandle,"New instance of web_queue at %s from %s pid %d ppid %d time %d\n",timestr,tmpdir,getpid(),getppid(),curTime);
  keepaliveHandle = fopen(keepaliveFile,"wt");
  fprintf(keepaliveHandle,"New instance of web_queue at %s from %s pid %d ppid %d time %d\n",timestr,tmpdir,getpid(),getppid(),curTime);
  fclose(keepaliveHandle);
  keepaliveHandle = NULL;
  result = chmod(keepaliveFile,S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH);
  if (result != 0) {
    printf("ERROR: Failed to change protection on %s errno: %d %s\n",keepaliveFile,errno,strerror(errno));
  }

  while (1) {
    iteration++;
    time(&curTime);
    ptr = localtime(&curTime);
    strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
#if 0
    fprintf(masterLogHandle,"Iteration %d Timestamp %s\n",iteration,timestr);
#endif
    /* Now go look for work */
    workDone = 1;
    while (workDone == 1) {
      workDone = 0;
      keepaliveHandle = fopen(keepaliveFile,"rt");
      if (keepaliveHandle == NULL) {
        fprintf(masterLogHandle,"Exiting on Iteration %d Timestamp %s\n",iteration,timestr);
        exit(-1);
      } else {
        fclose(keepaliveHandle);
      }

      directoryCount = scandir(tmpDirectoryRoot,&directoryList,0,alphasort);
      if (directoryCount >= 0) {
        if (directoryCount > requestAlloc) {
          /* Here we need to reallocate a request table */
          if (requestAlloc == 0) {
            requestAlloc = directoryCount + LIST_REALLOC_INCREMENT;
            pCommon->requestTable = (PREQUESTQUEUE)calloc(requestAlloc,sizeof(REQUESTQUEUE));
            if (pCommon->requestTable == NULL) {
              printf("ERROR: failed allocate a request table of size %d\n",requestAlloc);
              exit(-1);
            }
          } else {
            requestAlloc = directoryCount + LIST_REALLOC_INCREMENT;
            tmpRequestTable = (PREQUESTQUEUE)realloc(pCommon->requestTable,requestAlloc*sizeof(REQUESTQUEUE));
            if (tmpRequestTable == NULL) {
              printf("ERROR: failed to realloc request queue of size %d\n",requestAlloc);
              exit(-1);
            }
            pCommon->requestTable = tmpRequestTable;
            tmpRequestTable = NULL;
          }
          
        }

        pBestRequest = NULL;
        for (directoryIndex = 0; directoryIndex < directoryCount; directoryIndex++) {
          if (pNewRequest->coo != 0) {
            printf("ERROR memory leak from stale object list in line %d\n",__LINE__);
          }
          memset(pNewRequest,0,sizeof(REQUESTQUEUE));
          strcpy(pNewRequest->subdir,directoryList[directoryIndex]->d_name);
          /* Go to the next if we have already handled this request */
          requestSkip = 0;
          for (requestIndex = 0; requestIndex < pCommon->requestCount; requestIndex++) {
            pRequest = &pCommon->requestTable[requestIndex];
            if (strcmp(pRequest->subdir,pNewRequest->subdir) == 0) {
              if ((pRequest->haveResult != 0) ||
                  (pRequest->requestError != 0)) {
                requestSkip = 1;
                break;
               
              }
            }
          }
          if (requestSkip == 0) {

            if (strstr(directoryList[directoryIndex]->d_name,".") != NULL) {
              pNewRequest->requestError = 1;
            }
            if (strlen(directoryList[directoryIndex]->d_name) > (MAX_PATH-2)) {
              pNewRequest->requestError = 1;
            }
            if (pNewRequest->requestError == 0) {

              sprintf(pNewRequest->subDirectoryName,"%s%s/",tmpDirectoryRoot,pNewRequest->subdir);
              strcpy(pNewRequest->resultFileName,pNewRequest->subDirectoryName);
              strcat(pNewRequest->resultFileName,"result.txt");
              strcpy(tmpResultFileName,pNewRequest->subDirectoryName);
              strcat(tmpResultFileName,"tmpresult.txt");
#if 0
              if (strstr(pNewRequest->subDirectoryName,"1dBfjq") != 0) {
                printf("DEBUG at %s\n",requestFileName);
              }
#endif
              if (pCommon->verbose) {
                printf("In subdirectory %5d %5d %s\n",iteration,directoryIndex,pNewRequest->subDirectoryName);
              }
              strcpy(requestFileName,pNewRequest->subDirectoryName);
              strcat(requestFileName,"request.txt");
              statResult = stat(requestFileName,&statbuf);
              if (statResult != 0) {
                pNewRequest->requestError = 1;             
              } else {
                pNewRequest->haveRequest = 1;
                pNewRequest->requestSeconds = statbuf.st_mtim.tv_sec;
                pNewRequest->requestNsec = statbuf.st_mtim.tv_nsec;
              }
            }
            if (pNewRequest->requestError == 0) {
              statResult = stat(pNewRequest->resultFileName,&statbuf);
              if (statResult == 0) {
                if (pCommon->verbose) {
                  printf("Result file already exists %s\n",pNewRequest->resultFileName);
                }            
                pNewRequest->haveResult = 1;
                pNewRequest->resultSeconds = statbuf.st_mtim.tv_sec;
                pNewRequest->resultNsec = statbuf.st_mtim.tv_nsec;
              }
            }

            if (pNewRequest->requestError == 0) {
              if (pCommon->verbose) {
                printf("Found %s\n",requestFileName);

              }
              requestHandle = fopen(requestFileName,"rt");
              if (requestHandle == NULL)  {
                pNewRequest->requestError == 1;
              } 
            }
            lineCounter = 0;
            while (pNewRequest->requestError == 0) {
              inBuffer = fgets(inLine,MAX_BUFFER,requestHandle);
              if (inBuffer == NULL) {
                break;
              }
              lineCounter++;
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
              switch(lineCounter) {
              case 1: /* version line */
                nvals = sscanf(inBuffer,"version %d",&pNewRequest->requestVersion);
                if (nvals == 0) {
                  pNewRequest->requestVersion = 0;
                }
                if (pNewRequest->requestVersion != CURRENT_REQUEST_VERSION) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Stale request: file %s has the wrong version have %d expected %d<br />",requestFileName,pNewRequest->requestVersion,CURRENT_REQUEST_VERSION);
                  WriteResultFile(pCommon,pNewRequest);                
                  break;
                }
                if (pCommon->verbose) {
                  printf("Found %s\n",requestFileName);
                }
                break;
              case 2: /* nmin line */
                nvals = sscanf(inBuffer,"%d",&pNewRequest->nmin);
                if (nvals != 1) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse nmin in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                }
                break;
              case 3: /* box line */
                nvals = sscanf(inBuffer,"%d",&pNewRequest->box);
                if (nvals != 1) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse box in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                }
                break;
              case 4: /* whenqueued line */
                nvals = sscanf(inBuffer,"%lld",&pNewRequest->whenqueued);
                if (nvals != 1) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse whenqueued in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                }
                break;
              case 5: /* source line */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse source in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  strcpy(pNewRequest->source,inBuffer);
                }
                break;
              case 6: /* researcher line */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse researcher in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  strcpy(pNewRequest->researcher,inBuffer);
                }
                break;
              case 7: /* frameformat line */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse frameformat in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  strcpy(pNewRequest->frameformat,inBuffer);
                }
                break;
              case 8: /* binaries line */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse binaries in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  strcpy(pNewRequest->binaries,inBuffer);
                }
                break;
              case 9: /* lsbin line */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse lsbin in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  strcpy(pNewRequest->lsbin,inBuffer);
                }
                break;
              case 10: /* tsocks line */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse tsocks in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  strcpy(pNewRequest->tsocks,inBuffer);
                }
                break;
              case 11: /* tmpsubdir line */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse tmpsubdir in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  strcpy(pNewRequest->tmpsubdir,inBuffer);
                  if (strcmp(pNewRequest->subdir,pNewRequest->tmpsubdir) != 0) {
                    pNewRequest->requestError = 1;
                    sprintf(pCommon->errorBuffer,"<br />ERROR: Discrepancy of subdir %s and tmpsubdir %s in %s<br />",pNewRequest->subdir,pNewRequest->tmpsubdir,requestFileName);
                    WriteResultFile(pCommon,pNewRequest);                

                  }
                }
                break;
              case 12: /* remoteaddr line */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to parse remoteaddr in %s<br />",requestFileName);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  strcpy(pNewRequest->remoteaddr,inBuffer);
                }
                break;
              default: /* coo (coordinates) */
                if (strlen(inBuffer) > (MAX_BUFFER-2)) {
                  pNewRequest->requestError = 1;
                  sprintf(pCommon->errorBuffer,"<br />ERROR: coordinates in  %s exceed line length %d<br />",requestFileName,MAX_BUFFER);
                  WriteResultFile(pCommon,pNewRequest);                
                } else {
                  if (pNewRequest->coo == NULL) {
                    pNewRequest->coo = calloc(strlen(inBuffer)+2,1);
                    if (pNewRequest->coo == NULL) {
                      sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to allocate coordinates in %s<br />",requestFileName);
                      pNewRequest->requestError = 1;
                      WriteResultFile(pCommon,pNewRequest);                
                    } else {
                      strcpy(pNewRequest->coo,inBuffer);
                      pNewRequest->objectCount++;
                    }
                  } else {
                    requestSize = strlen(pNewRequest->coo)+strlen(inBuffer)+10;
                    if (requestSize > MAX_REQUEST) {
                      pNewRequest->requestError = 1;
                      sprintf(pCommon->errorBuffer,"<br />ERROR: coordinates in  %s exceed total length %d<br />",requestFileName,MAX_REQUEST);
                      WriteResultFile(pCommon,pNewRequest);                

                    } else {
                      tmpcoo = calloc(requestSize,1);
                      if (tmpcoo == NULL) {
                        sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to allocate coordinates in %s<br />",requestFileName);
                        pNewRequest->requestError = 1;
                        WriteResultFile(pCommon,pNewRequest);                
                      } else {
                        strcpy(tmpcoo,pNewRequest->coo);
                        strcat(tmpcoo,"\r\n");
                        strcat(tmpcoo,inBuffer);
                        pNewRequest->objectCount++;
                        free(pNewRequest->coo);
                        pNewRequest->coo = tmpcoo;
                        tmpcoo = NULL;
                      }
                    }
                  }


                }
                break;
              }
            }
            if (requestHandle != NULL) {
              fclose(requestHandle);
              requestHandle = NULL;
            }
            if ((pNewRequest->requestError == 0) && (pNewRequest->coo == NULL)) {
              /* No request was entered */
              sprintf(pCommon->errorBuffer,"<br />ERROR: Requested coordinates field is blank %s<br />",requestFileName);
              pNewRequest->requestError = 1;
              WriteResultFile(pCommon,pNewRequest);                
            
            } 
            /* Strip out any HTML insertion attacks */
            
            charPtr = pNewRequest->coo;
            if (charPtr != NULL) {
              charVal = *charPtr;
              while (charVal != 0) {
                if ((charVal == '<') || (charVal == '>')) {
                  *charPtr = '_';
                }
                charPtr++;
                charVal = *charPtr;              
              }
            }
            /* At this point, see if we already have this request */
            errorRequestIndex = -1;
            finishedRequestIndex = -1;
            latestRequestIndex = -1;
            for (requestIndex = 0; requestIndex < pCommon->requestCount; requestIndex++) {
              pRequest = &pCommon->requestTable[requestIndex];
              if (pNewRequest->requestError == 0) {

                if (pRequest->requestError != 0) {
                  errorRequestIndex = requestIndex;
                } else if (((pRequest->haveResult != 0) || 
                            (pRequest->requestHandled != 0)) &&
                           (pRequest->requestInProgress == 0) &&
                           (pNewRequest->requestInProgress != 0) &&
                           (pNewRequest->haveResult == 0) &&
                           (pNewRequest->requestHandled != 0)) {
                  if (finishedRequestIndex < 0) {
                    finishedRequestIndex = requestIndex;
                  } else {
                    pRequest3 = &pCommon->requestTable[finishedRequestIndex];
                    if (pRequest->whenqueued > pRequest3->whenqueued) {
                      finishedRequestIndex = requestIndex;
                    }
                  }
                } else if ((pRequest->requestInProgress == 0) &&
                           (pNewRequest->haveResult == 0)) {
                  if (latestRequestIndex < 0) {
                    latestRequestIndex = requestIndex;
                  } else {
                    pRequest3 = &pCommon->requestTable[latestRequestIndex];
                    if (pRequest->whenqueued > pRequest3->whenqueued) {
                      latestRequestIndex = requestIndex;
                    }
                  }
                }

              }
            
          
              if (strcmp(pRequest->subdir,pNewRequest->subdir) == 0) {
                break;
              }
            }
            if (requestIndex == pCommon->requestCount) {
              /* New request appeared */
              if (pCommon->requestCount >= MAX_QUEUE_ENTRIES) {
                pCommon->busyFlag = 1; /* We have to discard something */
                pRequest3 = NULL;
                if (errorRequestIndex > 0) {
                  pRequest3 = &pCommon->requestTable[errorRequestIndex];
                } else if (finishedRequestIndex > 0) {
                  pRequest3 = &pCommon->requestTable[finishedRequestIndex];
                } else if (latestRequestIndex > 0) {
                  pRequest3 = &pCommon->requestTable[latestRequestIndex];
                }
                if (pRequest3 != NULL) {
                  /* We have a candidate, so substitute this request */
                  WriteStatusFiles(pCommon,NULL,pRequest3);

                  if (pRequest3->coo != NULL) {
                    free(pRequest3->coo);
                    pRequest3->coo = NULL;
                  }
                  memcpy(pRequest3,pNewRequest,sizeof(REQUESTQUEUE));
                  pNewRequest->coo = NULL;
                } else {
                  WriteStatusFiles(pCommon,NULL,pNewRequest);
                }
              } else {
                pCommon->requestCount++;
                pRequest = &pCommon->requestTable[requestIndex];
                if (pRequest->coo != NULL) {
                  free(pRequest->coo);
                  pRequest->coo = NULL;
                }
                memcpy(pRequest,pNewRequest,sizeof(REQUESTQUEUE));
                pNewRequest->coo = NULL;
              }
            } else {
              /* Old request, ignore it if it has already been handled.  Otherwise, update it */
              if ((pRequest->requestHandled != 0) &&
                  (pRequest->requestInProgress != 0)) {
                pRequest->haveResult = pNewRequest->haveResult;
                pRequest->resultSeconds = pNewRequest->resultSeconds;
                pRequest->resultNsec = pNewRequest->resultNsec;
              } else {
                if (pRequest->coo != NULL) {
                  free(pRequest->coo);
                  pRequest->coo = NULL;
                }
                memcpy(pRequest,pNewRequest,sizeof(REQUESTQUEUE));
                pNewRequest->coo = NULL;
              }
            }
            pRequest->requestDirectoryFound = 1;
            if (pRequest->requestError) {
              pRequest->requestHandled = 1;
            }
            
          } /* Not skipping the request */
          if (pNewRequest->coo != NULL) {
            free(pNewRequest->coo);
            pNewRequest->coo = NULL;
          }
          free(directoryList[directoryIndex]);
        }
        if (pNewRequest->coo != NULL) {
          printf("WARNING: memory leak at line %d\n",__LINE__);
        }
        free(directoryList);
        directoryList = NULL;
      }
      /* Go through the list and delete any directories that have since disappeared */

      for (requestIndex = 0; requestIndex < pCommon->requestCount; requestIndex++) {
        pRequest = &pCommon->requestTable[requestIndex];
        if (pRequest->requestDirectoryFound == 0) {
          if (requestIndex != (pCommon->requestCount -1)) {
            pRequest2 = &pCommon->requestTable[pCommon->requestCount -1];
            if (pRequest->coo != NULL) {
              free(pRequest->coo);
              pRequest->coo = NULL;
            }
            memcpy(pRequest,pRequest2,sizeof(REQUESTQUEUE));
            pRequest2->coo = NULL;
          }
          pCommon->requestCount--;
          if (pCommon->requestCount < (MAX_QUEUE_ENTRIES-1)) {
            pCommon->busyFlag = 0;
          }
        }
      }
      /* Now sort the queue in increasing time */
      qsort((void*)pCommon->requestTable,pCommon->requestCount,sizeof(REQUESTQUEUE),RequestCompare);
      /* Pick the first good entry */
      pBestRequest = NULL;
      for (requestIndex = 0; requestIndex < pCommon->requestCount; requestIndex++) {
        pRequest = &pCommon->requestTable[requestIndex];
        if ((pRequest->requestHandled == 0) &&
            (pRequest->requestInProgress == 0) &&
            (pRequest->requestError == 0) &&
            (pRequest->haveRequest == 1) &&
            (pRequest->haveResult == 0)) {
          pBestRequest = pRequest;
          break;
        }
      }
      if ((pBestRequest != NULL) && 
          ((RELEASE_LEVEL == (RELEASE_LEVEL_ALL+0)) || (activeJobs == 0))) {
        workDone = 1;
        pBestRequest->requestHandled = 1;
        printf("Earliest request is in %s\n",pBestRequest->subDirectoryName);
        printf("nmin %d, box %d, whenqueued %lld, source %s,researcher %s,frameformat %s,coo %s, binaries %s, lsbin %s, tsocks %s tmpsubdir %s\n",pBestRequest->nmin,pBestRequest->box,pBestRequest->whenqueued,pBestRequest->source,pBestRequest->researcher,pBestRequest->frameformat,pBestRequest->coo,pBestRequest->binaries,pBestRequest->lsbin,pBestRequest->tsocks,pBestRequest->tmpsubdir); 
        cmdStr = calloc((2*MAX_PATH)+MAX_BUFFER+strlen(pBestRequest->coo),1);
        if (cmdStr == NULL) {
          sprintf(pCommon->errorBuffer,"<br />ERROR: Failed to allocate a command string of length %d for %s<br />",MAX_BUFFER+strlen(pBestRequest->coo),requestFileName);
          pBestRequest->requestError = 1;
          WriteResultFile(pCommon,pBestRequest);                
        } else {
          pBestRequest->requestInProgress = 1;
          activeJobs++;
          pid = fork();
          switch(pid) {
          case -1:
            printf("Fork failed with status %d %s\n",errno,strerror(errno));
            return(-1);
            break;
          case 0:
#ifdef CHILD_DEBUG
            exit(-1);
#else /* CHILD_DEBUG */
            ChildProcess(pCommon,pBestRequest,cmdStr);
#endif /* CHILD_DEBUG */
            break;
          default:
#ifdef CHILD_DEBUG
            ChildProcess(pCommon,pBestRequest,cmdStr);
#endif /* CHILD_DEBUG */
            pBestRequest->childPid = pid;
          } /* End of case statement */
          free(cmdStr);
        }
        
      }
    }
    if (activeJobs > 0) {
      pid = waitpid(-1,&stat_val,WNOHANG);
      if (pid > 0) {
        for (requestIndex = 0; requestIndex < pCommon->requestCount; requestIndex++) {
          pRequest = &pCommon->requestTable[requestIndex];
          if (pRequest->childPid == pid) {
            break;
          }
        }
        if (requestIndex != pCommon->requestCount) {
          if (pCommon->verbose) {
            printf("PARENT: Child pid %d exited for %s",pid,pRequest->subdir);
            if (WIFEXITED(stat_val)) {
              printf(" with code %d\n",WEXITSTATUS(stat_val));
            } else {
              printf(" abnormally\n");
            }
          }
          if (pRequest->requestInProgress == 0) {
            printf("ERROR: Child pid %d exited for %s has no request in progress\n");
          } else {
            activeJobs--;
            pRequest->requestInProgress = 0;
          }

        } else {
          printf("PARENT: unknown hild pid %d exited",pid);
        }
      }
    }
    if (pCommon->verbose) {
      for (requestIndex = 0; requestIndex < pCommon->requestCount; requestIndex++) {
        pRequest = &pCommon->requestTable[requestIndex];
        printf("index %d whenqueued %lld, directory found %d handled %d inProgress %d error %d haveRequest %d haveResult %d subdir %s\n",
               requestIndex,
               pRequest->whenqueued,
               pRequest->requestDirectoryFound,
               pRequest->requestHandled,
               pRequest->requestInProgress,
               pRequest->requestError,
               pRequest->haveRequest,
               pRequest->haveResult,
               pRequest->subdir);
        
      }
    }
    WriteStatusFiles(pCommon,NULL,NULL);
    
    fflush(masterLogHandle);
    sleep(5);

  }



  return(0);
}

