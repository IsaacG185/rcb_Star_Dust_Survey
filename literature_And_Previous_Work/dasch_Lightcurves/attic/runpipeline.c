// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 * Given a list of files, run with the requested number of fork processes until the list has
 * been exhausted.  Each fork process will process one mosaic until completion.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <mpi.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <dirent.h>
#include "mysql.h"
#include "pipelineutils.h"

/* #define LOS_DEBUG 1 */
/* #define SHOW_DISK_STATS 1 */
/* #define DEBUG_MOSAIC_SEARCH 1 */

#define MAX_FILENAME 512
#define MAX_BUFFER 700
#define MAX_ALTLOGFILES 3
#define MAX_LIST_STRING 50
#define WCSSOURCE_LOGBOOK 1
#define WCSSOURCE_IMWCS 2
#define NETWORK_MULTIPLIER 1000L
#define MAX_CPUS 512

/* Message tags */

#define WORKTAG 1
#define RESPONSETAG 2
#define STATUSTAG 3
#define DIETAG 4

/* Types of work */
#define WORK_NONE              0
#define WORK_FLATFIELD         1
#define WORK_LINEARITY         2
#define WORK_MOSAIC            3
#define WORK_FIND_ASTROMETRY   4
#define WORK_RUN_ASTROMETRYWCS 5
#define WORK_RUN_FULL_CSH      6
#define WORK_RUN_FULLDB_CSH    7
#define WORK_RUN_SYSTEM_CMD    8
#define MAX_WORK_TYPE          9

static char* workTypeString[MAX_WORK_TYPE] = {
  (char *) "NONE     ",
  (char *) "Flat     ",
  (char *) "Linearity",
  (char *) "Mosaic   ",
  (char *) "FindAstr ",
  (char *) "AstWCS   ",
  (char *) "Full     ",
  (char *) "FullDB   ",
  (char *) "SystemCmd",
};

#define REQUEUE_FLAG_NOTRUN 0
#define REQUEUE_FLAG_STARTED 1 /* Transient state, reset to NOTRUN */
#define REQUEUE_FLAG_COMPLETE 2

/* Make GCC happy about potential overflows */
#define BUFPRINTF(BUF, FMT, ARGS...) \
  do { \
    int rv = snprintf((BUF), sizeof(BUF), (FMT), ARGS); \
    if (rv < 0 || rv >= sizeof(BUF)) \
      abort(); \
  } while (0)

/* Work queue */
typedef struct _workqueue {
  int queueIndex;      /* Index into the work queue */
  int workType;        /* Type of work requested */
  int skipFlag;        /* Skip this item if nonzero */
  int workResult;      /* Status of the work */
  int diskLocation;
  int plateNumber;     /* plate number */
  int scanNumber;      /* scan or mosaic number */
  int mosaicNumber;    /* New mosaic number for FitsMosaic step */
  int rotation;        /* mosaic rotation */
  char series[MAX_SERIES_STRING]; /* Series */
  char plateName[MAX_LIST_STRING]; /* Plate name */
  int solutionNumber; /* current WCS solution */
  int delay;          /* Linearity test delay in microseconds */
  int fileCheckActive; /* CheckFullProgress is looking for a file */
  int fileCheckIteration; /* CheckFullProgress file iteration */
  int fileCheck1OK;         /* CheckFullProgress allobjects file found */
  int fileCheck2OK;        /* CheckFullProgress next sextractor file found */
  int FitWCS;              /* Status of full.csh */
  int requeueFlag;         /* requeue status */

  /* The response follows */
  int errorFlag;
  off_t priFilesize;
  off_t altFilesize[MAX_ALTLOGFILES];
  int seconds;
  char priLogname[MAX_BUFFER];
  char systemCommandBuffer[MAX_BUFFER];
  char altLogname[MAX_ALTLOGFILES][MAX_BUFFER];
} WORKQUEUE, *PWORKQUEUE;

typedef struct _nodestats {
  int workType;   /* Last completed work type */
  char nodename[MPI_MAX_PROCESSOR_NAME];
  /* Overall Cumulative stats */
  int taskscompleted;
  int totalSeconds;
  int elapsedSeconds;
  /* Stats since the last STATUSTAG message */
  int newTaskscompleted;
  int newTotalSeconds;
  int newElapsedSeconds;
  long long rxrate; /* mB/sec received */
  long long txrate; /* mB/sec sent */
  long long rxdiskrate; /* mB/sec received */
  long long txdiskrate; /* mB/sec sent */

} NODESTATS, *PNODESTATS;

#define MAX_MESSAGE_SIZE sizeof(WORKQUEUE)

/* Global common for all processes */

typedef struct _runcommon {
  /* MPI specific fields */
  int myrank;             /* Process rank (0 = master) */
  char processorName[MPI_MAX_PROCESSOR_NAME];  /* CPU name */
  int namelength;                     /* Length of CPU name */
  int maxproc;                        /* Number of processes */
  int verbose;                        /* Verbose print mode */
  int pauseMode;                      /* Wait until start.txt gets written at the end of the file copy */
  int noGSC;                          /* Not processing GSC data (apass or kepler only) */
  int startrank;                      /* Start of the loop */
  int statusRequestNumber; /* Number of status request cycles */
  int queueIteration;  /* Work queue iteration */
  int mosaicLocation;
  int tileLocation; /* Tile and calibration location (suffix required) */
  int skipFileCheck;
  int enableRequeue;  /* Enable Odyssey requeue option if 1 */

  /* Old runpipline fields */
  time_t beginTime;  /* Overall execution time */
  char timestr[100]; /* timestamp */
  char runtag[MAX_FILENAME];  /* run tag for multiple instances of runpipeline */
  char listname[MAX_FILENAME]; /* List of plates */
  int solutionNumber; /* current WCS solution */
  char startDateString[MAX_DATE_STRING+2];
  char endDateString[MAX_DATE_STRING+2];
  char mosaicDateString[MAX_DATE_STRING+2]; /* If mosaics are built, then advance the startDate for subsequent phases */
  char suffix[MAX_FLATFIELDS_SUFFIX+1];
  char password[MAX_BUFFER]; /* command line override of DASCH_PASSWORD */
  char mysqlhost[MAX_BUFFER]; /* command line override of DASCH_MYSQLHOST */
  char username[MAX_BUFFER];  /* command line override of DASCH_USERNAME */
  int runFindAstrometry; /* Run astrometry.net, first pass */
  int runAstrometryWCS1; /* Run AstrometryWCS, first pass */
  int runAstrometryWCS2; /* Run AstrometryWCS, additional passes */
  int runCalibration;    /* Run flats and linearity test */
  int runSystemCommand;  /* Run a system command */
  int buildMosaics;      /* Create mosaics */
  int reprocessTiles;    /* Reprocess tiles */
  int runFull1;          /* Run full.csh, first pass */
  int runFull2;          /* Run full.csh, additional passes */
  int multipleMode;      /* Perform multiple pass processing */
  int databaseFlag1;     /* Run fulldb.csh instead of full.csh, first pass */
  int databaseFlag2;     /* Run fulldb.csh instead of full.csh, additional passes */
  int canceledFlag;      /* If nonzero, the user is requesting an abort */
  int completeCount;     /* Number of successfully competed entries in an iteration */
  int completeCountMax;  /* Maximum number of successfully completed entries for all iterations */
  char keepalivename[MAX_FILENAME];
  char systemCommandBuffer[MAX_BUFFER];
  char qualifier[MAX_QUALIFIER];
  PWORKQUEUE pWorkQueue; /* Work queue */
  int workQueueSize;
  FILE *logFileHandle;
  int remoteCPUTotal;
  char remoteProcessorName[MAX_CPUS][MPI_MAX_PROCESSOR_NAME];  /* CPU name */
  int remotePlatesCompleted[MAX_CPUS];
} RUNCOMMON, *PRUNCOMMON;

void ReprocessCleanup(PRUNCOMMON pCommon, PWORKQUEUE preceive_result);
int master(PRUNCOMMON pCommon);
int slave(PRUNCOMMON pCommon);
void ChildProcess(PRUNCOMMON pCommon, PWORKQUEUE pwork, MPI_Status *status, int *pSeconds, int *pNewSeconds);

#ifdef SHOW_DISK_STATS
void
ShowDiskStats(time_t beginTime, int stamp)
{
  time_t curTime;
  struct tm *ptr;
  char timestr[100];

  time(&curTime);
  ptr = localtime(&curTime);
  strftime(timestr, 25, "%Y-%m-%dT%H-%M-%S", ptr);
  printf("Timestamp number %d is %d seconds for %s\n", stamp, (int)(curTime-beginTime), timestr);
  system("cat /proc/diskstats");
}
#endif /* SHOW_DISK_STATS */

void
CleanWorkQueue(PRUNCOMMON pCommon)
{
  if (pCommon->pWorkQueue != NULL) {
    free(pCommon->pWorkQueue);
    pCommon->pWorkQueue = NULL;
  }

  pCommon->workQueueSize = 0;
}

void
ReadLogFile(PWORKQUEUE preceive_result)
{
  int statResult;
  struct stat filestats;
  int iteration = 180;
  FILE *logFile;
  FILE *altFile = NULL;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  char altfilename[MAX_BUFFER];
  int altLogIndex;

  altfilename[0] = 0;
  statResult = -1;

  if (preceive_result->errorFlag) {
    printf("ERROR: errorFlag is set for %s\n", preceive_result->priLogname);
  } else {
    while (statResult != 0) {
      statResult = stat(preceive_result->priLogname, &filestats);

      if (statResult != 0) {
        printf("WARNING: master did not find %s\n", preceive_result->priLogname);
      } else if (filestats.st_size != preceive_result->priFilesize) {
        printf("WARNING: master expected size %ld got size %ld for %s\n",
               preceive_result->priFilesize,
               filestats.st_size,
               preceive_result->priLogname);
        statResult = -1;
      } else {
        break;
      }

      sleep(1);

      if (iteration-- <= 0) {
        break;
      }
    }

    if (iteration <= 0) {
      printf(
        "ERROR: Master could not read %s of size %ld\n",
        preceive_result->priLogname,
        preceive_result->priFilesize
      );
      return;
    }

    logFile = fopen(preceive_result->priLogname, "rt");
    if (logFile == NULL) {
      printf("ERROR: Master could not open %s\n", preceive_result->priLogname);
      return;
    }

    while (1) {
      inBuffer = fgets(inLine, MAX_BUFFER, logFile);
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

      printf("%s\n", inBuffer);
    }

    fclose(logFile);
    unlink(preceive_result->priLogname);

    for (altLogIndex = 0; altLogIndex < MAX_ALTLOGFILES; altLogIndex++) {
      altfilename[0] = 0;

      if (preceive_result->altLogname[altLogIndex][0] != 0) {
        if (preceive_result->altFilesize[altLogIndex] <= 0) {
          printf(
            "ERROR: alternate file size is %ld for %s\n",
            preceive_result->altFilesize[altLogIndex],
            preceive_result->altLogname[altLogIndex]
          );
        } else {
          /* Here we have an alternate log file */
          iteration = 180;
          statResult = -1;

          while (statResult != 0) {
            statResult = stat(preceive_result->altLogname[altLogIndex], &filestats);

            if (statResult != 0) {
              printf("WARNING: master did not find %s\n", preceive_result->altLogname[altLogIndex]);
            } else if (filestats.st_size != preceive_result->altFilesize[altLogIndex]) {
              printf("WARNING: master expected size %ld got size %ld for %s\n",
                     preceive_result->altFilesize[altLogIndex],
                     filestats.st_size,
                     preceive_result->altLogname[altLogIndex]);
              statResult = -1;
            } else {
              break;
            }

            sleep(1);
            if (iteration-- <= 0) {
              break;
            }
          }

          if (iteration <= 0) {
            printf(
              "ERROR: Master could not read %s of size %ld\n",
              preceive_result->altLogname[altLogIndex],
              preceive_result->altFilesize[altLogIndex]
            );
            return;
          }

          logFile = fopen(preceive_result->altLogname[altLogIndex], "rt");
          if (logFile == NULL) {
            printf("ERROR: Master could not open %s\n", preceive_result->altLogname[altLogIndex]);
            return;
          }

          while (1) {
            inBuffer = fgets(inLine, MAX_BUFFER, logFile);
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

            /* The first line is the true file name */
            if (altfilename[0] == 0) {
              strcpy(altfilename, inBuffer);
              altFile = fopen(altfilename, "a+t");
              if (altFile == NULL) {
                printf("ERROR: failed to open %s - dumping to stdout\n", altfilename);
              }
            } else {
              if (altFile == NULL) {
                printf("%s\n", inBuffer);
              } else {
                fprintf(altFile, "%s\n", inBuffer);
              }
            }
          }

          if (altFile != NULL) {
            fclose(altFile);
          }

          fclose(logFile);
          unlink(preceive_result->altLogname[altLogIndex]);
        }
      }
    }
  }
}


int
ProcessWorkQueue(PRUNCOMMON pCommon, int maxproc)
{
  int rank;
  PWORKQUEUE pWorkEntry;
  PWORKQUEUE pReceiveWorkEntry;
  int totalSeconds = 0;
  int newTotalSeconds = 0;
  PWORKQUEUE preceive_result;
  char receive_buffer[MAX_MESSAGE_SIZE];
  MPI_Status status;
  int mpitestflag;
  MPI_Request request;
  int mpiresult;
  int numpending = 0;
  int numstarted = 0;
  int numprocessed = 0;
  int workIndex;
  struct stat statbuf;
  int result;
  time_t startTime;
  time_t curTime;
  int sendNewWorkRequest = 0;
  PNODESTATS pstatusmessage;
  double averageTime;
  int workType = 0;
  char rxstring[MAX_BUFFER];
  char txstring[MAX_BUFFER];
  char rxdiskstring[MAX_BUFFER];
  char txdiskstring[MAX_BUFFER];
  int count;
  int cpuIndex;

  if (pCommon->verbose) {
    printf("Have maxproc %d\n", pCommon->maxproc);
  }

  if (pCommon->workQueueSize <= 0 || pCommon->canceledFlag != 0) {
    /* Nothing to do! */
    return 0;
  }

  pCommon->queueIteration++;
  pCommon->completeCount = 0;
  time(&startTime);

  /* Add these fields for sanity */
  for (workIndex = 0; workIndex < pCommon->workQueueSize; workIndex++) {
    pWorkEntry = &pCommon->pWorkQueue[workIndex];
    pWorkEntry->solutionNumber = pCommon->solutionNumber;
  }

  /* Seed the slaves by asking them to report their status */
  for (cpuIndex = 0; cpuIndex < pCommon->remoteCPUTotal; cpuIndex++) {
    printf("CPU %8s Completed %5d\n", pCommon->remoteProcessorName[cpuIndex], pCommon->remotePlatesCompleted[cpuIndex]);
    pCommon->remotePlatesCompleted[cpuIndex] = 0;
  }

  pCommon->statusRequestNumber++;

  for (rank = pCommon->startrank; rank < pCommon->maxproc; ++rank) {
    mpiresult = MPI_Send(0, 0, MPI_BYTE, rank, STATUSTAG, MPI_COMM_WORLD);
    numpending++;
    if (mpiresult != MPI_SUCCESS) {
      printf("ERROR: Master sending status request to %d with result %d\n", pCommon->maxproc, mpiresult);
    }
  }

  /* Loop over getting new work requests until there is no more work
     to be done */
  pWorkEntry = pCommon->pWorkQueue;

  while (numpending > 0) {
    /* Receive results from a slave */
    sleep(1);

    mpiresult = MPI_Irecv(
      &receive_buffer,  /* message buffer */
      MAX_MESSAGE_SIZE, /* one data item */
      MPI_BYTE,         /* of type double real */
      MPI_ANY_SOURCE,   /* receive from any sender */
      MPI_ANY_TAG,      /* any type of message */
      MPI_COMM_WORLD,   /* default communicator */
      &request
    );

    if (mpiresult == MPI_SUCCESS) {
      while (1) {
        mpiresult = MPI_Test(&request, &mpitestflag, &status);
        if (mpiresult != MPI_SUCCESS) {
          break;
        }

        if (mpitestflag) {
          break;
        }

        sleep(1);
      }
    }

    MPI_Get_count(&status, MPI_BYTE, &count);

    if (mpiresult != MPI_SUCCESS) {
      printf("ERROR: Master receiving tag %d, err %d, count %3d, cancelled %d mpiresult %d from source %d",
             status.MPI_TAG, status.MPI_ERROR, count, status._cancelled, mpiresult, status.MPI_SOURCE);
    } else {
      if (pCommon->verbose) {
        printf("Master receiving tag %d, err %d, count %3d, cancelled %d mpiresult %d from source %d\n",
               status.MPI_TAG, status.MPI_ERROR, count, status._cancelled, mpiresult, status.MPI_SOURCE);
      }
    }

    sendNewWorkRequest = 0;

    switch (status.MPI_TAG) {
    case WORKTAG:
      /* Process the work immediately */
      if (sizeof(WORKQUEUE) == count) {
        pReceiveWorkEntry = (PWORKQUEUE) &receive_buffer;
        if (pCommon->verbose) {
          printf(" work item %d from itself\n", pReceiveWorkEntry->queueIndex);
        }

        ChildProcess(pCommon, pReceiveWorkEntry, &status, &totalSeconds, &newTotalSeconds);
      } else {
        printf("ERROR: master expected %lu bytes, got %d bytes\n", sizeof(WORKQUEUE), count);
      }
      break;

    case RESPONSETAG:
      preceive_result = (PWORKQUEUE) &receive_buffer;
      if (pCommon->verbose) {
        printf(" value %3d work item %3d\n", preceive_result->workResult, preceive_result->queueIndex);
      }

      /* Open and print out the results buffer */
      ReadLogFile(preceive_result);
      if (pCommon->reprocessTiles != 0) {
        ReprocessCleanup(pCommon, preceive_result);
      }

      numpending--;
      numprocessed++;
      pCommon->completeCount++;

      if (pCommon->runSystemCommand  == 0) {
        sendNewWorkRequest = 1;
      }
      break;

    case STATUSTAG:
      numpending--;

      if (count == 0) {
        printf("Master receiving STATUSTAG\n");
      } else if (sizeof(NODESTATS) == count) {
        pstatusmessage = (PNODESTATS) &receive_buffer;
        sprintf(rxstring, "%lld", pstatusmessage->rxrate);
        sprintf(txstring, "%lld", pstatusmessage->txrate);
        sprintf(rxdiskstring, "%lld", pstatusmessage->rxdiskrate);
        sprintf(txdiskstring, "%lld", pstatusmessage->txdiskrate);

        for (cpuIndex = 0; cpuIndex < pCommon->remoteCPUTotal; cpuIndex++) {
          if (strcmp(pstatusmessage->nodename, pCommon->remoteProcessorName[cpuIndex]) == 0) {
            break;
          }
        }

        if (cpuIndex < MAX_CPUS) {
          if (cpuIndex == pCommon->remoteCPUTotal) {
            strcpy(pCommon->remoteProcessorName[cpuIndex], pstatusmessage->nodename);
            pCommon->remoteCPUTotal++;
          }
          pCommon->remotePlatesCompleted[cpuIndex] += pstatusmessage->newTaskscompleted;
        } else {
          printf("ERROR: cpuIndex %d greater than MAX_CPUS\n", cpuIndex);
          exit(-1);
        }

        printf(
          "Recstat %d node %3d %s %-9s proc: %4d %4d Totalsec %6d %6d Elapsec %6d %6d nrx %6s ntx %6s kB drx %6s dtx %6s kB/sec\n",
          pCommon->statusRequestNumber,
          status.MPI_SOURCE,
          pstatusmessage->nodename,
          workTypeString[pstatusmessage->workType],
          pstatusmessage->taskscompleted,
          pstatusmessage->newTaskscompleted,
          pstatusmessage->totalSeconds,
          pstatusmessage->newTotalSeconds,
          pstatusmessage->elapsedSeconds,
          pstatusmessage->newElapsedSeconds,
          rxstring,
          txstring,
          rxdiskstring,
          txdiskstring
        );
      } else {
        printf("ERROR: STATUSTAG with size %d received\n", count);
      }

      sendNewWorkRequest = 1;
      break;

    default:
      printf("ERROR: received unknown tag %d\n", status.MPI_TAG);
      break;
    }

    if (sendNewWorkRequest) {
      /* Find the next item of work to do */
      while (pWorkEntry->skipFlag != 0 || pWorkEntry->requeueFlag != REQUEUE_FLAG_NOTRUN) {
        pWorkEntry++;
        numstarted++;
      }

      if (pCommon->canceledFlag == 0 && numstarted < pCommon->workQueueSize) {
        /* Send it to each rank */
        if (pWorkEntry->workType != WORK_NONE) {
          workType = pWorkEntry->workType;
        }

        mpiresult = MPI_Send(
          pWorkEntry,        /* message buffer */
          sizeof(WORKQUEUE), /* one data item */
          MPI_BYTE,          /* data item is an integer */
          status.MPI_SOURCE, /* destination process rank */
          WORKTAG,           /* user chosen message tag */
          MPI_COMM_WORLD     /* default communicator */
        );

        if (mpiresult != MPI_SUCCESS) {
          printf("ERROR: Master send to %d with result %d\n", rank, mpiresult);
        } else {
          if (pCommon->logFileHandle) {
            fprintf(
              pCommon->logFileHandle,
              "%s s%d %s\n",
              pWorkEntry->plateName,
              pCommon->solutionNumber,
              workTypeString[pWorkEntry->workType]
            );
            fflush(pCommon->logFileHandle);
          }

          time(&curTime);
          curTime -= pCommon->beginTime;
          printf(
            "PARENT: %17s Created subprocess %3d on %3d at %ld seconds\n",
            pWorkEntry->plateName,
            numstarted,
            rank,
            curTime
          );

          pWorkEntry++;
          numstarted++;
          numpending++;
        }
      }

      if (pCommon->verbose) {
        printf("numpending %d, numstarted %d, numprocessed %d (1)\n", numpending, numstarted, numprocessed);
      }
    }

    if (pCommon->verbose) {
      printf("numpending %d, numstarted %d, numprocessed %d (3)\n", numpending, numstarted, numprocessed);
    }

    result = stat(pCommon->keepalivename, &statbuf);

    if (result != 0) {
      printf(
        "Terminating at line %d because keepalive file %s is missing\n",
        numstarted + 1,
        pCommon->keepalivename
      );
      pCommon->canceledFlag = 1;
    } else if (pCommon->verbose) {
      printf("Stat on %s is %d\n", pCommon->keepalivename, result);
    }
  }

  if (pCommon->verbose) {
    printf("numpending %d, numstarted %d, numprocessed %d (4)\n", numpending, numstarted, numprocessed);
  }

  pCommon->queueIteration++;
  time(&curTime);
  curTime -= startTime;

  if (pCommon->completeCount == 0) {
    averageTime = 0.0;
  } else {
    averageTime = (1.0 * curTime) / pCommon->completeCount;
  }

  printf(
    "Phase %d Type %-9s Time: %ld seconds (%.1f sec/plate) for %d mosaics completed out of %d mosaics %d processors\n",
    pCommon->queueIteration,
    workTypeString[workType],
    curTime,
    averageTime,
    pCommon->completeCount,
    pCommon->workQueueSize,
    pCommon->maxproc
  );

  return 0;
}


/* Read the fitted.list file left by AstrometryWCS, return the number of entries
 * found */
int
CreateFittedList(PRUNCOMMON pCommon)
{
  char fittedFile[MAX_BUFFER];
  char cmdStr[MAX_BUFFER];
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  int curline = 0;
  FILE *fittedHandle;
  int nlines = 0;
  char *scriptDir;
  PWORKQUEUE pWorkEntry;
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  char series[MAX_SERIES_STRING];

  CleanWorkQueue(pCommon);

  scriptDir = getenv("DASCH_SCRIPTS");
  if (scriptDir == NULL) {
    printf("ERROR: DASCH_SCRIPTS is not defined\n");
    return 0;
  }

  strcpy(fittedFile, scriptDir);
  strcat(fittedFile, "/fitted.list");

  fittedHandle = fopen(fittedFile, "rt");
  if (fittedHandle == NULL) {
    printf("runpipeline failed to open fitted list %s line %d\n", fittedFile, __LINE__);
    return 0;
  }

  /* Make a backup of the fitted.list file */
  sprintf(cmdStr, "cp %s %s/fitted%d.tmp\n", fittedFile, scriptDir, pCommon->solutionNumber);
  system(cmdStr);

  /* Count the number of records in the file */
  while (1) {
    inBuffer = fgets(inLine, MAX_BUFFER, fittedHandle);
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

    if (lineLen > 5) {
      nlines++;
    }
  }

  fclose(fittedHandle);

  if (nlines == 0) {
    printf("Fitted file %s has no entries\n", fittedFile);
    return 0;
  }

  fittedHandle = fopen(fittedFile, "rt");
  if (fittedHandle == NULL) {
    printf("runpipeline failed to open fitted list %s line %d\n", fittedFile, __LINE__);
    return 0;
  }

  pCommon->pWorkQueue = (PWORKQUEUE) calloc(nlines + 2, sizeof(WORKQUEUE));
  if (pCommon->pWorkQueue == NULL) {
    fprintf(stderr, "ERROR: failed to allocate pCommon->pWorkQueue within CreateFittedList\n");
    return 0;
  }

  while (1) {
    inBuffer = fgets(inLine, MAX_BUFFER, fittedHandle);
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

    if (lineLen > 5) {
      pWorkEntry = &pCommon->pWorkQueue[curline];
      memset(pWorkEntry, 0, sizeof(WORKQUEUE));

      if (strlen(inBuffer) >= MAX_LIST_STRING - 1) {
        printf("ERROR: %s exceeds MAX_LIST_STRING\n", inBuffer);
      } else {
        if (ParseFilename(inBuffer, series, &plateNumber, &mosaicNumber, &binning, &rotation)) {
          strcpy(pWorkEntry->plateName, inBuffer);
          pWorkEntry->queueIndex = curline;
          strcpy(pWorkEntry->series, series);
          pWorkEntry->plateNumber = plateNumber;
          pWorkEntry->scanNumber = mosaicNumber;
          curline++;
        } else {
          printf("ERROR: failed to parse %s from fitted.list\n", inBuffer);
        }
      }
    }
  }

  fclose(fittedHandle);

  pCommon->workQueueSize = curline;
  return curline;
}

int
fitsfilter(const struct dirent *direntry)
{
  if (direntry->d_type == DT_REG || direntry->d_type == DT_UNKNOWN) {
    if (strstr(direntry->d_name, ".fit") != NULL) {
      return 1;
    }
  }

  return 0;
}

int
restartfilter(const struct dirent *direntry)
{
  if (direntry->d_type == DT_REG || direntry->d_type == DT_UNKNOWN) {
    if (strstr(direntry->d_name, ".xxx") != NULL) {
      return 1;
    }
  }

  return 0;
}

/* Build a list of flatframes and linearity tests to process */
int
CreateCalibrationList(PRUNCOMMON pCommon)
{
  MYSQL my_connection;
  char *mysqlhost;
  char *username;
  char *password;
  char filename[MAX_BUFFER];
  PFLATFIELDS FlatFieldTable = NULL;
  PFLATFIELDS pFlatfields;
  PFLATFIELDS pCurFlatfields;
  int flatIndex;
  int totalCount = 0;
  int foundCount = 0;
  int allocCount;
  PWORKQUEUE pWorkEntry;
  struct stat statbuf;
  int result;
  struct dirent **filelist;
  int fileIndex;
  int reprocessCount = 0;
  int numFiles;

  CleanWorkQueue(pCommon);

  if (pCommon->mysqlhost[0] != 0) {
    mysqlhost = &pCommon->mysqlhost[0];
  } else {
    mysqlhost = getenv("DASCH_MYSQLHOST");
  }

  if (mysqlhost == NULL) {
    fprintf(stderr, "DASCH_MYSQLHOST is not defined\n");
    return 1;
  }

  if (pCommon->username[0] != 0) {
    username = &pCommon->username[0];
  } else {
    username = getenv("DASCH_USERNAME");
  }

  if (username == NULL) {
    fprintf(stderr, "DASCH_USERNAME is not defined\n");
    return 1;
  }

  if (pCommon->password[0] != 0) {
    password = &pCommon->password[0];
  } else {
    password = getenv("DASCH_PASSWORD");
  }

  if (password == NULL) {
    fprintf(stderr, "DASCH_PASSWORD is not defined\n");
    return 1;
  }

  if (pCommon->pWorkQueue != NULL) {
    fprintf(stderr, "ERROR: pCommon->pWorkQueue already allocated\n");
    return 1;
  }

  mysql_init(&my_connection);

  if (mysql_real_connect(&my_connection, mysqlhost, username, password, "scanner", 0, NULL, CLIENT_FOUND_ROWS)) {
    if (pCommon->reprocessTiles) {
      allocCount = GetFlatfieldRecords(
        &my_connection,
        pCommon->startDateString,
        pCommon->endDateString,
        -1,
        -1,
        -1,
        NULL,
        NULL,
        NULL,
        &FlatFieldTable
      );
    } else {
      allocCount = GetFlatfieldRecords(
        &my_connection,
        pCommon->startDateString,
        pCommon->endDateString,
        -1,
        CALIBRATION_STATUS_ACQUIRED,
        -1,
        NULL,
        NULL,
        NULL,
        &FlatFieldTable
      );
    }

    if (allocCount > 0) {
      pCommon->pWorkQueue = (PWORKQUEUE) calloc(allocCount + 2, sizeof(WORKQUEUE));
      if (pCommon->pWorkQueue == NULL) {
        fprintf(stderr, "ERROR: failed to allocate pCommon->pWorkQueue within CreateListFile\n");
        return 1;
      }

      if (pCommon->reprocessTiles) {
        /* Reprocessing tiles -- Take only one entry for each type */
        pCurFlatfields = &FlatFieldTable[0];
        pCurFlatfields->alternateDate[0] = 0; /* Ignore alternateDate */
        pCurFlatfields->diskLocation = pCommon->tileLocation; /* Ignore disk */
        pCurFlatfields->delay = 0;   /* Ignore delay */
        pCurFlatfields->calibrationStatus = CALIBRATION_STATUS_ACQUIRED; /* Ignore calibration status */
        pCurFlatfields->JobId = -1; /* Ignore JobId */

        if (strcmp(pCurFlatfields->suffix, pCommon->suffix) != 0) {
          pCurFlatfields->suffix[0] = 0; /* Save only the suffix of interest */
        }

        reprocessCount = 1;

        for (flatIndex = 1; flatIndex < allocCount; flatIndex++) {
          pFlatfields = &FlatFieldTable[flatIndex];

          if (
            strcmp(pFlatfields->calibrationDate, pCurFlatfields->calibrationDate) == 0 &&
            pFlatfields->calibrationType == pCurFlatfields->calibrationType &&
            strcmp(pFlatfields->suffix, pCurFlatfields->suffix) == 0
          ) {
            /* This is a duplicate, ignore it */
            continue;
          } else {
            if (
              strcmp(pFlatfields->calibrationDate, pCurFlatfields->calibrationDate) == 0 &&
              pFlatfields->calibrationType == pCurFlatfields->calibrationType
            ) {
              /* Here only the suffix is different.  Save the suffix of interest */
              if (strcmp(pFlatfields->suffix, pCommon->suffix) == 0) {
                strcpy(pCurFlatfields->suffix, pFlatfields->suffix);
              }

              continue;
            }

            pCurFlatfields = &FlatFieldTable[reprocessCount];
            if (reprocessCount != flatIndex) {
              memcpy(pCurFlatfields, pFlatfields, sizeof(FLATFIELDS));
            }

            pCurFlatfields->alternateDate[0] = 0; /* Ignore alternateDate */
            pCurFlatfields->diskLocation = pCommon->tileLocation; /* Ignore disk */
            pCurFlatfields->delay = 0;   /* Ignore delay */
            pCurFlatfields->calibrationStatus = CALIBRATION_STATUS_ACQUIRED; /* Ignore calibration status */
            pCurFlatfields->JobId = -1; /* Ignore JobId */
            reprocessCount++;
          }
        }

        /* Reduce the size of the allocation */
        allocCount = reprocessCount;
      }

      for (flatIndex = 0; flatIndex < allocCount; flatIndex++) {
        pFlatfields = &FlatFieldTable[flatIndex];

        if (
          (pFlatfields->calibrationType != CALIBRATION_TYPE_LINEARITY && pFlatfields->calibrationType != CALIBRATION_TYPE_FLAT) ||
          pFlatfields->calibrationStatus != CALIBRATION_STATUS_ACQUIRED ||
          pFlatfields->diskLocation <= 0) {
          continue;
        }

        totalCount++;

        if (pFlatfields->calibrationType == CALIBRATION_TYPE_FLAT) {
          /* See if we have already done this work */
          sprintf(filename, "/dasch/raid%03d/ExposureData/FlatFrames/%s", pFlatfields->diskLocation, pFlatfields->calibrationDate);
          if (pFlatfields->suffix[0] != 0) {
            strcat(filename, pFlatfields->suffix);
          }

          result = stat(filename, &statbuf);
          if (result == 0) {
            /* Directory exists, now see if it has FITS files */
            numFiles = scandir(filename, &filelist, fitsfilter, alphasort);

            if (numFiles > 0) {
              /* Found files, just delete the list and ignore this entry */
              for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
                free(filelist[fileIndex]);
              }

              free(filelist);
              continue;
            }
          }

          /* Now see if the requested tile is present */
          sprintf(filename, "/dasch/raid%03d/ExposureData/FlatFrames/Flat_%s.fit", pFlatfields->diskLocation, pFlatfields->calibrationDate);
          result = stat(filename, &statbuf);
          if (result != 0) {
            /* File does not exist */
            continue;
          }

          pWorkEntry = &pCommon->pWorkQueue[foundCount];
          memset(pWorkEntry, 0, sizeof(WORKQUEUE));
          strcpy(pWorkEntry->plateName, pFlatfields->calibrationDate);
          pWorkEntry->workType = WORK_FLATFIELD;
          pWorkEntry->diskLocation = pFlatfields->diskLocation;
          pWorkEntry->queueIndex = foundCount;
          pWorkEntry->skipFlag = 0;
          foundCount++;

          if (foundCount > allocCount) {
            fprintf(stderr, "ERROR: foundCount  %d exceeds allocCount\n", foundCount);
            exit(-1);
          }
        } else if (pFlatfields->calibrationType == CALIBRATION_TYPE_LINEARITY) {
          /* See if we have already done this work */
          if (pFlatfields->delay == 0 || pFlatfields->delay == LINEARITY_DELAY_USEC) {
            sprintf(filename, "/dasch/raid%03d/ExposureData/FlatFrames/FlatStats_%s", pFlatfields->diskLocation, pFlatfields->calibrationDate);
          } else {
            sprintf(filename, "/dasch/raid%03d/ExposureData/FlatFrames/FlatStats_%sD%d", pFlatfields->diskLocation, pFlatfields->calibrationDate, pFlatfields->delay);
          }

          if (pFlatfields->suffix[0] != 0) {
            strcat(filename, pFlatfields->suffix);
          }

          strcat(filename, ".bin");

          result = stat(filename, &statbuf);
          if (result == 0) {
            /* Already did this one */
            continue;
          }

          /* See if the input directory exists */
          sprintf(filename, "/dasch/raid%03d/ExposureData/FlatStats/%s", pFlatfields->diskLocation, pFlatfields->calibrationDate);
          result = stat(filename, &statbuf);
          if (result == 0) {
            /* Directory exists, now see if it has FITS files */
            numFiles = scandir(filename, &filelist, fitsfilter, alphasort);

            if (numFiles > 0) {
              /* Found files, delete the list then process the data */
              for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
                free(filelist[fileIndex]);
              }

              free(filelist);
            } else {
              /* No fits files, just continue */
              continue;
            }
          } else {
            /* No raw data directory */
            continue;
          }

          pWorkEntry = &pCommon->pWorkQueue[foundCount];
          memset(pWorkEntry, 0, sizeof(WORKQUEUE));
          strcpy(pWorkEntry->plateName, pFlatfields->calibrationDate);
          pWorkEntry->workType = WORK_LINEARITY;
          pWorkEntry->diskLocation = pFlatfields->diskLocation;
          pWorkEntry->queueIndex = foundCount;
          pWorkEntry->skipFlag = 0;
          pWorkEntry->delay = pFlatfields->delay;
          foundCount++;

          if (foundCount > allocCount) {
            fprintf(stderr, "ERROR: foundCount  %d exceeds allocCount\n", foundCount);
            exit(-1);
          }
        } else {
          printf("ERROR: CreateCalibrationList has unknown type %d\n", pFlatfields->calibrationType);
          continue;
        }
      }

      if (foundCount == 0) {
        if (pCommon->pWorkQueue != NULL) {
          free(pCommon->pWorkQueue);
          pCommon->pWorkQueue = NULL;
        }
      }
    }

    pCommon->workQueueSize = foundCount;

    if (FlatFieldTable != NULL) {
      free(FlatFieldTable);
    }

    mysql_close(&my_connection);
  } else {
    fprintf(stderr, "Connection failed\n");

    if (mysql_errno(&my_connection)) {
      fprintf(stderr, "Connection error %d: %s\n", mysql_errno(&my_connection), mysql_error(&my_connection));
    }
  }

  printf(
    "CreateCalibrationList found %d database entries found for date %s and %s\n",
    foundCount,
    pCommon->startDateString,
    pCommon->endDateString
  );
  return 0;
}


/* Build a list from the MySQL database for files/mosaics created after a specific date
 *
 *  wcsFlag =  0  - the list should contain non-fitted scans from the mosaics table
 *  wcsFlag =  1  - the list should contain fitted mosaics from the mosaics table
 */
int
CreateListFile(PRUNCOMMON pCommon, int wcsFlag)
{
  int res;
  MYSQL my_connection;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char *series;
  int plateNumber;
  int mosaicNumber;
  int rotation;
  char *mysqlhost;
  char *username;
  char *password;
  char *scriptdirectory;
  char queryString[MAX_QUERY_STRING];
  int totalCount = 0;
  int foundCount = 0;
  int allocCount;
  int nvals;
  int WCSSource;
  PWORKQUEUE pWorkEntry;

  CleanWorkQueue(pCommon);

  if (pCommon->mysqlhost[0] != 0) {
    mysqlhost = &pCommon->mysqlhost[0];
  } else {
    mysqlhost = getenv("DASCH_MYSQLHOST");
  }

  if (mysqlhost == NULL) {
    fprintf(stderr, "DASCH_MYSQLHOST is not defined\n");
    return 1;
  }

  if (pCommon->username[0] != 0) {
    username = &pCommon->username[0];
  } else {
    username = getenv("DASCH_USERNAME");
  }

  if (username == NULL) {
    fprintf(stderr, "DASCH_USERNAME is not defined\n");
    return 1;
  }

  if (pCommon->password[0] != 0) {
    password = &pCommon->password[0];
  } else {
    password = getenv("DASCH_PASSWORD");
  }

  if (password == NULL) {
    fprintf(stderr, "DASCH_PASSWORD is not defined\n");
    return 1;
  }

  scriptdirectory = getenv("DASCH_SCRIPTS");

  if (pCommon->pWorkQueue != NULL) {
    fprintf(stderr, "ERROR: pCommon->pWorkQueue already allocated\n");
    return 1;
  }

  mysql_init(&my_connection);

  if (mysql_real_connect(&my_connection, mysqlhost, username, password, "scanner", 0, NULL, CLIENT_FOUND_ROWS)) {
    if (pCommon->reprocessTiles != 0 && pCommon->mosaicDateString[0] != 0) {
      sprintf(
        queryString,
        "SELECT series, plateNumber, mosaicNumber, rotation, WCSSource+0 "
        "from mosaics where mosaicDate > '%s' and naxis1 IS NOT NULL "
        "and series != 'mask' and ((mosaicComment IS NULL) or (not(mosaicComment regexp('deleted')))) ;",
        pCommon->mosaicDateString
      );
    } else {
      sprintf(
        queryString,
        "SELECT series, plateNumber, mosaicNumber, rotation, WCSSource+0 "
        "from mosaics where mosaicDate > '%s' and mosaicDate < '%s' and naxis1 IS NOT NULL "
        "and series != 'mask' and ((mosaicComment IS NULL) or (not(mosaicComment regexp('deleted')))) ;",
        pCommon->startDateString,
        pCommon->endDateString
      );
    }

    if (strlen(queryString) > MAX_QUERY_STRING) {
      fprintf(stderr, "ERROR: MAX_QUERY_STRING exceeded %lu\n", strlen(queryString));
      return 1;
    }

    res = mysql_query(&my_connection, queryString);

    if (!res) {
      res_ptr = mysql_store_result(&my_connection);

      if (res_ptr) {
        allocCount = (int) mysql_num_rows(res_ptr);
        if (allocCount > 0) {
          pCommon->pWorkQueue = (PWORKQUEUE) calloc(allocCount + 2, sizeof(WORKQUEUE));

          if (pCommon->pWorkQueue == NULL) {
            fprintf(stderr, "ERROR: failed to allocate pCommon->pWorkQueue within CreateListFile\n");
            return 1;
          }
        }

        while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
          totalCount++;

          if (sqlrow[0] == NULL) {
            fprintf(stderr, "No series found \n");
            continue;
          }

          series = sqlrow[0];
          nvals = sscanf(sqlrow[1], "%d", &plateNumber);
          if (nvals != 1) {
            fprintf(stderr, "ERROR: nvals is %d for plateNumber\n", nvals);
            continue;
          }

          nvals = sscanf(sqlrow[2], "%d", &mosaicNumber);
          if (nvals != 1) {
            fprintf(stderr, "ERROR: nvals is %d for mosaicNumber\n", nvals);
            continue;
          }

          if (sqlrow[3] != NULL) {
            nvals = sscanf(sqlrow[3], "%d", &rotation);
            if (nvals != 1) {
              fprintf(stderr, "ERROR: nvals is %d for rotation\n", nvals);
              continue;
            }
          } else {
            rotation = 0;
          }

          if (sqlrow[4] != NULL) {
            nvals = sscanf(sqlrow[4], "%d", &WCSSource);
            if (nvals != 1) {
              fprintf(stderr, "ERROR: nvals is %d for WCSSource\n", nvals);
              continue;
            }
          } else {
            WCSSource = 0;
          }

          if (allocCount > 0) {
            if (
              (wcsFlag == 0 && (WCSSource & WCSSOURCE_IMWCS) == 0) ||
              (wcsFlag != 0 && (WCSSource & WCSSOURCE_IMWCS) != 0)
            ) {
              pWorkEntry = &pCommon->pWorkQueue[foundCount];
              memset(pWorkEntry, 0, sizeof(WORKQUEUE));

              if ((WCSSource & WCSSOURCE_IMWCS) == 0) {
                sprintf(pWorkEntry->plateName, "%s%05d_%02d_01", series, plateNumber, mosaicNumber);
              } else if (rotation == 0) {
                sprintf(pWorkEntry->plateName, "%s%05d_%02d_01ww", series, plateNumber, mosaicNumber);
              } else {
                sprintf(pWorkEntry->plateName, "%s%05d_%02d_01r%dww", series, plateNumber, mosaicNumber, rotation);
              }

              pWorkEntry->workType = WORK_FIND_ASTROMETRY;
              pWorkEntry->queueIndex = foundCount;
              pWorkEntry->skipFlag = 0;
              strcpy(pWorkEntry->series, series);
              pWorkEntry->plateNumber = plateNumber;
              pWorkEntry->scanNumber = mosaicNumber;
              pWorkEntry->rotation = rotation;
              foundCount++;

              if (foundCount > allocCount) {
                fprintf(stderr, "ERROR: foundCount  %d exceeds allocCount\n", foundCount);
                exit(-1);
              }
            }
          }
        }

        pCommon->workQueueSize = foundCount;
        mysql_free_result(res_ptr);
      }
    } else {
      printf("Select error %d: %s res %d\n", mysql_errno(&my_connection), mysql_error(&my_connection), res);
    }

    mysql_close(&my_connection);
  } else {
    fprintf(stderr, "Connection failed\n");
    if (mysql_errno(&my_connection)) {
      fprintf(stderr, "Connection error %d: %s\n", mysql_errno(&my_connection), mysql_error(&my_connection));
    }
  }

  /* Write out runpipeline.tmp so we can check progress with monitor.csh */
  if (foundCount > 0) {
    char tmpFileName[MAX_FILENAME];
    FILE* tmpFileHandle = NULL;
    int index;

    BUFPRINTF(tmpFileName, "%s/runpipeline_%s.tmp", scriptdirectory, pCommon->runtag);
    tmpFileHandle = fopen(tmpFileName, "wt");

    if (tmpFileHandle == 0) {
      fprintf(stderr, "ERROR: failed to open temporary list file %s\n", tmpFileName);
    } else {
      for (index = 0; index < foundCount; index++) {
        pWorkEntry = &pCommon->pWorkQueue[index];
        fprintf(tmpFileHandle, "%s\n", pWorkEntry->plateName);
      }

      fclose(tmpFileHandle);
    }
  }

  printf(
    "CreateListFile found %d database entries found for date %s, WCS flag %d\n",
    foundCount,
    pCommon->startDateString,
    wcsFlag
  );
  return 0;
}


int
CreateMosaicList(PRUNCOMMON pCommon)
{
  int res;
  MYSQL my_connection;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char *series;
  int plateNumber;
  int scanNumber;
  char *mysqlhost;
  char *username;
  char *password;
  char queryString[MAX_QUERY_STRING];
  int totalCount = 0;
  int foundCount = 0;
  int allocCount;
  int nvals;
  int diskLocation;
  PWORKQUEUE pWorkEntry;
  char filename[MAX_BUFFER];
  struct dirent **filelist;
  int fileIndex;
  int numFiles;
  struct stat statbuf;
  int result;
  MOSAIC mosaicInfo;
  int mosaicNumber;
  int foundMosaic;
  int solutionNumber = 0;
  char lastScanDate[MAX_DATE_STRING];
  char curScanDate[MAX_DATE_STRING];
  int foundRowNoiseFix;
  int maxMosaicNumber;
  int verbose = 0;
  int pendingFlag = 0;
  int mosaicInfoFlag = 1;
  int selectAllFlag = 1;
  int numMosaicRecords = 0;
  PMOSAICLIST pMosaicTable = NULL;
  PMOSAIC pMosaic2 = &mosaicInfo;
  int mosaicListIndex;
  int totalSolution0Count = 0;
  int totalSolutionsCount = 0;

  CleanWorkQueue(pCommon);

  if (pCommon->mysqlhost[0] != 0) {
    mysqlhost = &pCommon->mysqlhost[0];
  } else {
    mysqlhost = getenv("DASCH_MYSQLHOST");
  }

  if (mysqlhost == NULL) {
    fprintf(stderr, "DASCH_MYSQLHOST is not defined\n");
    return 1;
  }

  if (pCommon->username[0] != 0) {
    username = &pCommon->username[0];
  } else {
    username = getenv("DASCH_USERNAME");
  }

  if (username == NULL) {
    fprintf(stderr, "DASCH_USERNAME is not defined\n");
    return 1;
  }

  if (pCommon->password[0] != 0) {
    password = &pCommon->password[0];
  } else {
    password = getenv("DASCH_PASSWORD");
  }

  if (password == NULL) {
    fprintf(stderr, "DASCH_PASSWORD is not defined\n");
    return 1;
  }

  if (pCommon->pWorkQueue != NULL) {
    fprintf(stderr, "ERROR: pCommon->pWorkQueue already allocated\n");
    return 1;
  }

  mysql_init(&my_connection);

  if (mysql_real_connect(&my_connection, mysqlhost, username, password, "scanner", 0, NULL, CLIENT_FOUND_ROWS)) {
    sprintf(
      queryString,
      "SELECT series, plateNumber, scanNumber, diskLocation from scans "
      "where scanDate > '%s' and scanDate < '%s'  and "
      "((scanComment IS NULL) or (not(scanComment regexp('deleted')))) ;",
      pCommon->startDateString,
      pCommon->endDateString
    );

    if (strlen(queryString) > MAX_QUERY_STRING) {
      fprintf(stderr, "ERROR: MAX_QUERY_STRING exceeded %lu\n", strlen(queryString));
      return 1;
    }

    res = mysql_query(&my_connection, queryString);

    if (!res) {
      res_ptr = mysql_store_result(&my_connection);

      if (res_ptr) {
        allocCount = (int) mysql_num_rows(res_ptr);

        if (allocCount > 0) {
          pCommon->pWorkQueue = (PWORKQUEUE) calloc(allocCount + 2, sizeof(WORKQUEUE));

          if (pCommon->pWorkQueue == NULL) {
            fprintf(stderr, "ERROR: failed to allocate pCommon->pWorkQueue within CreateListFile\n");
            return 1;
          }
        }

        while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
          totalCount++;

          if (sqlrow[0] == NULL) {
            fprintf(stderr, "No series found \n");
            continue;
          }

          series = sqlrow[0];
          nvals = sscanf(sqlrow[1], "%d", &plateNumber);
          if (nvals != 1) {
            fprintf(stderr, "ERROR: nvals is %d for plateNumber\n", nvals);
            continue;
          }

          nvals = sscanf(sqlrow[2], "%d", &scanNumber);
          if (nvals != 1) {
            fprintf(stderr, "ERROR: nvals is %d for scanNumber\n", nvals);
            continue;
          }

          if (scanNumber > MAX_SCAN_MOSAIC_MOSAIC_NUMBER) {
            printf("ERROR: MAX_SCAN_MOSAIC_MOSAIC_NUMBER exceeded by %d for %s%05d in runpipeline line %d\n", scanNumber, series, plateNumber, __LINE__);
          }

          if (pCommon->reprocessTiles == 0) {
            if (sqlrow[3] != NULL) {
              nvals = sscanf(sqlrow[3], "%d", &diskLocation);
              if (nvals != 1) {
                fprintf(stderr, "ERROR: nvals is %d for rotation\n", nvals);
                continue;
              }
            } else {
              continue;
            }
          } else {
            diskLocation = pCommon->tileLocation;
          }

          if (allocCount > 0) {
            /* Look for the tile directory */
            sprintf(filename, "/dasch/raid%03d/ExposureData/Plates/%s/%05d_%02d", diskLocation, series, plateNumber, scanNumber);
            result = stat(filename, &statbuf);

            if (result == 0) {
              /* Directory exists, now see if it has FITS files */
              numFiles = scandir(filename, &filelist, fitsfilter, alphasort);

              if (numFiles > 0) {
                /* Found files, just delete the list and ignore this entry */
                for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
                  free(filelist[fileIndex]);
                }

                free(filelist);
              } else {
                /* No fits files here. Skip */
                continue;
              }
            } else {
              /* No directory found. Skip */
              continue;
            }

            /* Check that this mosaic has not already been built */
            mosaicNumber = 0;
            foundMosaic = 0;
            foundRowNoiseFix = 0;
            maxMosaicNumber = -1;
            memset(lastScanDate, 0, sizeof(lastScanDate));
            memset(curScanDate, 0, sizeof(curScanDate));

            mosaicNumber = -1;
            BuildMosaicList(
              &my_connection,
              series,     /* if NULL or zero length, consider all series */
              plateNumber,  /* if -1, consider all plates */
              mosaicNumber,  /* if -1, consider all mosaics */
              solutionNumber, /* if -1, consider all solutions */
              pendingFlag, /* if 1, then return unscanned plates */
              mosaicInfoFlag, /* if 1, then being called from GetMosaicInfo: no exposure or plates table info */
              selectAllFlag,   /* if 1, then return all mosaic records: good, deleted, and stale */
              verbose,         /* if 1, display warnings */
              &totalSolutionsCount,
              &totalSolution0Count,
              &numMosaicRecords,
              &pMosaicTable
            );

            for (mosaicListIndex = 0; mosaicListIndex < numMosaicRecords; mosaicListIndex++) {
              pMosaic2 = &pMosaicTable[mosaicListIndex];
              maxMosaicNumber = pMosaic2->maxMosaicNumber;

              if (pMosaic2->scanNumber == scanNumber) {
                strcpy(curScanDate, pMosaic2->scanDate);
                foundMosaic = 1;
              }

              if (strstr(pMosaic2->mosaicComment, "Deleted") == NULL) {
                if (strstr(pMosaic2->mosaicComment, "RowNoiseFix") != NULL) {
                  foundRowNoiseFix = 1;
                }

                if (lastScanDate[0] == 0) {
                  strcpy(lastScanDate, pMosaic2->scanDate);
                } else {
                  result = strcmp(pMosaic2->scanDate, lastScanDate);
                  if (result > 0) {
                    strcpy(lastScanDate, pMosaic2->scanDate);
                  }
                }
              }
            }

            if (pMosaicTable != NULL) {
              free(pMosaicTable);
              pMosaicTable = NULL;
            }

            if (pCommon->reprocessTiles == 0) {
              if (foundMosaic) {
                continue;
              }
            } else {
              /* Row noise already corrected */
              if (foundRowNoiseFix == 1) {
                printf("WARNING: skipping row noise fix for %s%05d scanNumber %02d because it has already been fixed\n", series, plateNumber, scanNumber);
                continue;
              }

              /* Do not fix the row noise on masks */
              if (strcmp(series, "mask") == 0) {
                printf("WARNING: skipping row noise fix for %s%05d scanNumber %02d because of the 'mask' series\n", series, plateNumber, scanNumber);
                continue;
              }

              if (foundMosaic == 0) {
                printf("WARNING: skipping row noise fix for %s%05d scanNumber %02d because no valid mosaic built\n", series, plateNumber, scanNumber);
              } else {
                result = strcmp(lastScanDate, curScanDate);
                if (result > 0) {
                  printf("WARNING: skipping row noise fix for %s%05d scanNumber %02d because curScanDate %s is earlier than lastScanDate %s \n", series, plateNumber, scanNumber, curScanDate, lastScanDate);
                  continue; /* Found a later scan of this plate */
                }
              }
            }

            mosaicNumber = maxMosaicNumber+1;
            pWorkEntry = &pCommon->pWorkQueue[foundCount];
            memset(pWorkEntry, 0, sizeof(WORKQUEUE));
            pWorkEntry->workType = WORK_MOSAIC;
            pWorkEntry->queueIndex = foundCount;
            pWorkEntry->skipFlag = 0;
            strcpy(pWorkEntry->series, series);
            pWorkEntry->plateNumber = plateNumber;
            pWorkEntry->scanNumber = scanNumber;

            if (mosaicNumber > scanNumber) { /* Bugfix of Apr 29, 2015 */
              pWorkEntry->mosaicNumber = mosaicNumber;
            } else {
              pWorkEntry->mosaicNumber = scanNumber;
            }

            if (mosaicNumber > MAX_SCAN_MOSAIC_MOSAIC_NUMBER) {
              printf("ERROR: MAX_SCAN_MOSAIC_MOSAIC_NUMBER exceeded by %d for %s%05d in runpipeline line %d\n", mosaicNumber, series, plateNumber, __LINE__);
            }

            pWorkEntry->diskLocation = diskLocation;
            sprintf(pWorkEntry->plateName, "%s%05d_%02d_01", series, plateNumber, scanNumber);

            foundCount++;

            if (foundCount > allocCount) {
              fprintf(stderr, "ERROR: foundCount  %d exceeds allocCount\n", foundCount);
              exit(-1);
            }
          }
        }

        pCommon->workQueueSize = foundCount;
        mysql_free_result(res_ptr);
      }
    } else {
      printf("Select error %d: %s res %d\n", mysql_errno(&my_connection), mysql_error(&my_connection), res);
    }

    mysql_close(&my_connection);
  } else {
    fprintf(stderr, "Connection failed\n");

    if (mysql_errno(&my_connection)) {
      fprintf(stderr, "Connection error %d: %s\n", mysql_errno(&my_connection), mysql_error(&my_connection));
    }
  }

  printf(
    "CreateMosaicList found %d database entries found for date %s and %s\n",
    foundCount,
    pCommon->startDateString,
    pCommon->endDateString
  );
  return 0;
}


int
BuildSystemCommandList(PRUNCOMMON pCommon)
{
  int foundCount;
  PWORKQUEUE pWorkEntry;

  CleanWorkQueue(pCommon);

  pCommon->pWorkQueue = (PWORKQUEUE) calloc(pCommon->maxproc, sizeof(WORKQUEUE));

  if (pCommon->pWorkQueue == NULL) {
    fprintf(stderr, "ERROR: failed to allocate pCommon->pWorkQueue within CreateListFile\n");
    return 1;
  }

  for (foundCount = 0; foundCount < pCommon->maxproc; foundCount++) {
    pWorkEntry = &pCommon->pWorkQueue[foundCount];
    memset(pWorkEntry, 0, sizeof(WORKQUEUE));
    pWorkEntry->workType = WORK_RUN_SYSTEM_CMD;
    pWorkEntry->queueIndex = foundCount;
    pWorkEntry->skipFlag = 0;
    strcpy(pWorkEntry->systemCommandBuffer, pCommon->systemCommandBuffer);
  }

  pCommon->workQueueSize = foundCount;

  printf(
    "BuildSystemCommandList found %d database entries\n",
    foundCount
  );
  return 0;
}


/* This routine looks for *_allobjects.db file names, makes sure that they are
 * at least 1000 bytes long and returns 0 if it finds at least one and 1 if it
 * does not find any. If the "-n" qualifier is specified, we are performing
 * apass or kepler only processing. Look for a Sextractor file for the next
 * iteration.
 */
int
CheckFullProgress(PRUNCOMMON pCommon)
{
  int curline;
  char *ingestDirectory;
  char *astrometryDirectory;
  char allobjects_name[MAX_FILENAME];
  char astrometry_name[MAX_FILENAME];
  char tmpFileName[MAX_FILENAME];
  FILE *tmpFileHandle = NULL;
  struct stat statbuf;
  int result;
  int fileCount = 0;
  char *scriptdirectory;
  PWORKQUEUE pWorkEntry;
  char *matchDirectory;
  MYSQL my_connection;
  char *mysqlhost;
  char *username;
  char *password;
  int haveProgress;

  if (pCommon->mysqlhost[0] != 0) {
    mysqlhost = &pCommon->mysqlhost[0];
  } else {
    mysqlhost = getenv("DASCH_MYSQLHOST");
  }

  if (mysqlhost == NULL) {
    fprintf(stderr, "DASCH_MYSQLHOST is not defined\n");
    return 1;
  }

  if (pCommon->username[0] != 0) {
    username = &pCommon->username[0];
  } else {
    username = getenv("DASCH_USERNAME");
  }

  if (username == NULL) {
    fprintf(stderr, "DASCH_USERNAME is not defined\n");
    return 1;
  }

  if (pCommon->password[0] != 0) {
    password = &pCommon->password[0];
  } else {
    password = getenv("DASCH_PASSWORD");
  }

  if (password == NULL) {
    fprintf(stderr, "DASCH_PASSWORD is not defined\n");
    return 1;
  }

  scriptdirectory = getenv("DASCH_SCRIPTS");
  if (scriptdirectory == NULL) {
    fprintf(stderr, "ERROR: DASCH_SCRIPTS is not defined\n");
    return -1;
  }

  astrometryDirectory = getenv("DASCH_ASTROMETRY");
  if (astrometryDirectory == NULL) {
    fprintf(stderr, "ERROR: DASCH_ASTROMETRY is not defined\n");
    return -1;
  }

  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    fprintf(stderr, "ERROR: DASCH_INGEST is not defined\n");
    return -1;
  }

  matchDirectory = getenv("DASCH_MATCH");
  if (matchDirectory == NULL) {
    fprintf(stderr, "ERROR: DASCH_MATCH is not defined\n");
    return -1;
  }

  mysql_init(&my_connection);

  if (mysql_real_connect(&my_connection, mysqlhost, username, password, "scanner", 0, NULL, CLIENT_FOUND_ROWS)) {
    BUFPRINTF(tmpFileName, "%s/runpipeline%d_%s.tmp", scriptdirectory, pCommon->solutionNumber + 1, pCommon->runtag);
    tmpFileHandle = fopen(tmpFileName, "wt");
    if (tmpFileHandle == 0) {
      fprintf(stderr, "ERROR: failed to open temporary list file %s\n", tmpFileName);
    }

    /* First loop, initialize */
    for (curline = 0; curline < pCommon->workQueueSize; curline++) {
      pWorkEntry = &pCommon->pWorkQueue[curline];
      pWorkEntry->fileCheckActive = 1;
      pWorkEntry->fileCheckIteration = 180;
      pWorkEntry->fileCheck1OK = 0;
      pWorkEntry->fileCheck2OK = 0;

      if (pCommon->noGSC == 0) {
        if (SetMosaicFitWCS(&my_connection, "NoAllobjects", pWorkEntry->plateName, pCommon->solutionNumber, 0, 1, &pWorkEntry->FitWCS) != 0) {
          /* Nothing to look for */
          pWorkEntry->fileCheckActive = 0;
          pWorkEntry->fileCheckIteration = 0;
          pWorkEntry->skipFlag = 1;
          allobjects_name[0] = 0;
        } else {
          if ((pWorkEntry->FitWCS & FITWCS_NOALLOBJECTS) == 0) {
            pWorkEntry->fileCheckIteration = 180;
          } else {
            /* Nothing to look for */
            pWorkEntry->fileCheckActive = 0;
            pWorkEntry->fileCheckIteration = 0;
            pWorkEntry->skipFlag = 1;
            allobjects_name[0] = 0;
          }
        }
      }
    }

    /* Second loop, look for allobjects */
    while (1) {
      haveProgress = 0;

      for (curline = 0; curline < pCommon->workQueueSize; curline++) {
        pWorkEntry = &pCommon->pWorkQueue[curline];
        if (pWorkEntry->skipFlag != 0) {
          continue;
        }

        if (pWorkEntry->fileCheckActive == 0) {
          continue;
        }

        if (pWorkEntry->fileCheck1OK != 0) {
          continue;
        }

        if (pCommon->skipFileCheck) {
          pWorkEntry->fileCheck1OK = 1;
        } else {
          if (pCommon->noGSC == 0) {
            if (pCommon->solutionNumber > 0) {
              sprintf(allobjects_name, "%s/%s_s%d_allobjects.db", ingestDirectory, pWorkEntry->plateName, pCommon->solutionNumber);
            } else {
              sprintf(allobjects_name, "%s/%s_allobjects.db", ingestDirectory, pWorkEntry->plateName);
            }
          } else if (pCommon->solutionNumber > 0) {
            sprintf(allobjects_name, "%s/%s_s%d_%s_allobjects.db", ingestDirectory, pWorkEntry->plateName, pCommon->solutionNumber, pCommon->qualifier);
          } else {
            sprintf(allobjects_name, "%s/%s_%s_allobjects.db", ingestDirectory, pWorkEntry->plateName, pCommon->qualifier);
          }

          result = -1;
          result = stat(allobjects_name, &statbuf);

          if (result == 0 && statbuf.st_size > 1000) {
            pWorkEntry->fileCheck1OK = 1;
          } else if (pWorkEntry->fileCheckIteration-- <= 0) {
            pWorkEntry->fileCheckActive = 0;
          } else {
            haveProgress = 1;
          }
        }
      }

      if (haveProgress) {
        sleep(1);
      } else {
        break;
      }
    }

    /* Third loop, initialize for an astrometry result */
    for (curline = 0; curline < pCommon->workQueueSize; curline++) {
      pWorkEntry = &pCommon->pWorkQueue[curline];

      if (pWorkEntry->skipFlag != 0) {
        continue;
      }

      if (pWorkEntry->fileCheckActive == 0) {
        continue;
      }

      if (pWorkEntry->fileCheck1OK == 0) {
        continue;
      }

      if (pWorkEntry->fileCheck2OK != 0) {
        continue;
      }

      if (pCommon->skipFileCheck) {
        pWorkEntry->fileCheck2OK = 1;
        pWorkEntry->fileCheckActive = 0;
      } else if (pCommon->runAstrometryWCS2 != 0) {
        /* See if find_astrometry2 produced a fit (Astrometry2Succeeded) */
        if ((pWorkEntry->FitWCS & FITWCS_ASTROMETRY2) != 0) {
          pWorkEntry->fileCheckIteration = 180;
        } else {
          pWorkEntry->fileCheckIteration = 1; /* Try at least once if the MySQL entry is old and stale */
        }
      } else {
        /* find_astrometry2 never ran */
        pWorkEntry->fileCheck2OK = 1;
        pWorkEntry->fileCheckActive = 0;
      }
    }

    /* Fourth loop, look for for an astromtry result */
    while (1) {
      haveProgress = 0;

      for (curline = 0; curline < pCommon->workQueueSize; curline++) {
        pWorkEntry = &pCommon->pWorkQueue[curline];

        if (pWorkEntry->skipFlag != 0) {
          continue;
        }

        if (pWorkEntry->fileCheckActive == 0) {
          continue;
        }

        if (pWorkEntry->fileCheck1OK == 0) {
          continue;
        }

        if (pWorkEntry->fileCheck2OK != 0) {
          continue;
        }

        sprintf(astrometry_name, "%s/%s_s%d.db", astrometryDirectory, pWorkEntry->plateName, pCommon->solutionNumber+1);
        result = -1;
        result = stat(astrometry_name, &statbuf);

        if (result == 0) {
          pWorkEntry->fileCheck2OK = 1;
        } else if (pWorkEntry->fileCheckIteration-- <= 0) {
          pWorkEntry->fileCheckActive = 0;
        } else {
          haveProgress = 1;
        }
      }

      if (haveProgress) {
        sleep(1);
      } else {
        break;
      }
    }

    /* Fifth Loop, store results */
    for (curline = 0; curline < pCommon->workQueueSize; curline++) {
      pWorkEntry = &pCommon->pWorkQueue[curline];
      if (pWorkEntry->skipFlag) {
        continue;
      }

      if (pWorkEntry->fileCheck1OK == 0) {
        printf("CheckFullProgress failed to find an allobjects result for %s_s%d\n", pWorkEntry->plateName, pCommon->solutionNumber);
        pWorkEntry->skipFlag = 1;
      } else if (pWorkEntry->fileCheck2OK == 0) {
        printf("CheckFullProgress failed to find an astrometry result for %s_s%d\n", pWorkEntry->plateName, pCommon->solutionNumber+1);
        pWorkEntry->skipFlag = 1;
      } else {
        if (tmpFileHandle != NULL) {
          fprintf(tmpFileHandle, "%s\n", pWorkEntry->plateName);
        }
        fileCount++;
      }
    }

    printf(
      "CheckFullProgress found %d of %d files for solution %d\n",
      fileCount,
      pCommon->workQueueSize,
      pCommon->solutionNumber
    );

    if (tmpFileHandle != NULL) {
      fclose(tmpFileHandle);
    }

    mysql_close(&my_connection);
  } else {
    fprintf(stderr, "Connection failed\n");

    if (mysql_errno(&my_connection)) {
      fprintf(stderr, "Connection error %d: %s\n", mysql_errno(&my_connection), mysql_error(&my_connection));
    }
  }

  if (fileCount == 0) {
    return 1;
  } else {
    return 0;
  }
}


void
InitializeRestartState(PRUNCOMMON pCommon)
{
  struct dirent **filelist;
  int fileIndex;
  int numFiles;
  char *scriptDir;
  PWORKQUEUE pWorkEntry;
  int queueIndex;
  int lowerIndex = -1;
  int upperIndex = -1;
  int requeueNotrunCount = 0;
  int requeueStartedCount = 0;
  int requeueCompleteCount = 0;

  if (pCommon->enableRequeue == 0 || pCommon->workQueueSize == 0) {
    return;
  }

  scriptDir = getenv("DASCH_SCRIPTS");
  if (scriptDir == NULL) {
    printf("ERROR: line %d DASCH_SCRIPTS is not defined\n", __LINE__);
    exit(-1);
  }

  numFiles = scandir(scriptDir, &filelist, restartfilter, alphasort);

  if (numFiles > 0) {
    /* Found files, just delete the list and ignore this entry */
    for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
      for (queueIndex = 0; queueIndex < pCommon->workQueueSize; queueIndex++) {
        pWorkEntry = &pCommon->pWorkQueue[queueIndex];

        if (strstr(filelist[fileIndex]->d_name, pWorkEntry->plateName) != NULL) {
          printf("Found %s in work entry %d\n", filelist[fileIndex]->d_name, queueIndex);
          requeueStartedCount++;
          pWorkEntry->requeueFlag = REQUEUE_FLAG_STARTED;

          if (lowerIndex < 0) {
            lowerIndex = queueIndex;
            upperIndex = queueIndex;
          } else {
            if (queueIndex < lowerIndex) {
              lowerIndex = queueIndex;
            }

            if (queueIndex > upperIndex) {
              upperIndex = queueIndex;
            }
          }

          break;
        }
      }

      if (queueIndex >= pCommon->workQueueSize) {
        printf("Found %s, but not a work entry\n", filelist[fileIndex]->d_name);
      }

      free(filelist[fileIndex]);
    }

    free(filelist);

    /* Now set all of the flags */
    for (queueIndex = 0; queueIndex < pCommon->workQueueSize; queueIndex++) {
      pWorkEntry = &pCommon->pWorkQueue[queueIndex];

      if (queueIndex <= upperIndex) {
        if (pWorkEntry->requeueFlag != REQUEUE_FLAG_STARTED) {
          pWorkEntry->requeueFlag = REQUEUE_FLAG_COMPLETE;
          requeueCompleteCount++;
        } else {
          pWorkEntry->requeueFlag = REQUEUE_FLAG_NOTRUN;
        }
      } else {
        pWorkEntry->requeueFlag = REQUEUE_FLAG_NOTRUN;
        requeueNotrunCount++;
      }
    }

    printf(
      "found requeueCompleteCount %d requeueStartedCount %d requeueNotrunCount %d  lowerIndex %d upperIndex %d\n",
      requeueCompleteCount,
      requeueStartedCount,
      requeueNotrunCount,
      lowerIndex,
      upperIndex
    );
  }
}


int
main(int argc, char *argv[])
{
  RUNCOMMON pipelinecommon;
  PRUNCOMMON pCommon = &pipelinecommon;
  struct tm *ptr;
  int errorFlag = 0;
  char *argstr;
  char cmdchar;
  int nvals;
  char *argstr2;
  char cmdchar2;
  char pauseFile[MAX_BUFFER];
  char *scriptDir;
  struct stat statbuf;
  int result;
  time_t startTime = 0;
  time_t curTime;
  char *charPtr;
  char *charPtr2;
  int startDisk;
  char *atPtr;
  char tempusername[MAX_BUFFER];

  assert(sizeof(WORKQUEUE) <= MAX_MESSAGE_SIZE);
  assert(sizeof(NODESTATS) <= MAX_MESSAGE_SIZE);

  memset(pCommon, 0, sizeof(RUNCOMMON));
  pCommon->skipFileCheck = 1;
  pCommon->runFull1 = 1;
  pCommon->startrank = 1;

#ifdef SHOW_DISK_STATS
  ShowDiskStats(pCommon->beginTime, 1);
#endif /* SHOW_DISK_STATS */

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      if (strlen(pCommon->listname) > 0) {
        errorFlag = 1;
        printf("ERROR: list %s is being overwritten by %s\n", pCommon->listname, argstr);
      } else {
        strncpy(pCommon->listname, argstr, MAX_FILENAME - 1);
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 's': /* Solution Number */
        case 'S':
          argc--;
          nvals = sscanf(*++argv, "%d", &pCommon->solutionNumber);
          if (nvals != 1) {
            printf("ERROR: Can not decode solutionNumber\n");
            errorFlag = 1;
          }
          break;

        case 'd': /* Build lists from the MySQL database */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(pCommon->startDateString, *++argv, MAX_DATE_STRING);
            pCommon->startDateString[MAX_DATE_STRING] = 0;
          }
          break;

        case 'e': /* End date for the MySQL search */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(pCommon->endDateString, *++argv, MAX_DATE_STRING);
            pCommon->endDateString[MAX_DATE_STRING] = 0;
          }
          break;

        case 'E':
          pCommon->reprocessTiles = 1;
          argc--;

          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(pCommon->suffix, *++argv, MAX_FLATFIELDS_SUFFIX);
            pCommon->suffix[MAX_FLATFIELDS_SUFFIX] = 0;

            if (strlen(pCommon->suffix) == 0) {
              printf("ERROR: Zero length suffix specfied\n");
              errorFlag = 1;
            } else {
              charPtr = &pCommon->suffix[0];
              while (*charPtr != 0) {
                if ((*charPtr >= 'a' && *charPtr <= 'z') || (*charPtr >= 'A' && *charPtr <= 'Z')) {
                  charPtr++;
                } else {
                  printf("ERROR: Illegal suffix '%s' specified\n", pCommon->suffix);
                  errorFlag = 1;
                  break;
                }
              }
            }
          }
          break;

        case 'a': /* Run find_astrometry.csh */
        case 'A':
          pCommon->runFindAstrometry = 1;
          break;

        case 'c': /* Run calibration frames */
        case 'C':
          pCommon->runCalibration = 1;
          break;

        case 'g': /* Execute a system command */
        case 'G':
          pCommon->runSystemCommand = 1;
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(pCommon->systemCommandBuffer, *++argv, MAX_BUFFER);
            pCommon->systemCommandBuffer[MAX_BUFFER-1] = 0;
          }
          break;

        case 'm': /* Create mosaics */
        case 'M':
          pCommon->buildMosaics = 1;
          argc--;
          nvals = sscanf(*++argv, "%d", &pCommon->mosaicLocation);
          if (nvals != 1) {
            printf("ERROR: Can not decode mosaicLocation\n");
            errorFlag = 1;
          }
          break;

        case 'i': /* Disk with reprocessed tiles */
        case 'I':
          pCommon->reprocessTiles = 1;
          argc--;
          nvals = sscanf(*++argv, "%d", &pCommon->tileLocation);
          if (nvals != 1) {
            printf("ERROR: Can not decode tileLocation\n");
            errorFlag = 1;
          }
          break;

        case 'w': /* Run AstrometryWCS  */
          pCommon->runAstrometryWCS1 = 1;
          break;

        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            if (strlen(*++argv) < ((MAX_BUFFER / 2) - 2)) {
              strcpy(pCommon->password, *argv);
            } else {
              printf("ERROR argument is too long %s\n", *argv);
              errorFlag = 1;
            }
          }
          break;

        case 'x': /* Do not run full.csh */
        case 'X':
          pCommon->runFull1 = 0;
          break;

        case '2': /* Multiple pass processing */
          pCommon->multipleMode = 1;
          if (argc > 1 && (*(argv+1))[0] != '-') {
            argc--;
            ++argv;
            argstr2 = *argv;

            while ((cmdchar2 = *(argstr2++)) != 0) {
              switch(cmdchar2) {
              case 'w':
                pCommon->runAstrometryWCS2 = 1;
                break;

              case 'f':
                pCommon->runFull2 = 1;
                pCommon->databaseFlag2 = 0;
                break;

              case 'd':
                pCommon->runFull2 = 1;
                pCommon->databaseFlag2 = 1;
                break;

              default:
                printf("ERROR: Unknown multiple pass flag %c\n", cmdchar2);
                errorFlag = 1;
              }
            }
          }
          break;

        case 'v': /* Verbose mode */
        case 'V':
          pCommon->verbose = 1;
          break;

        case 'q': /* file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(pCommon->qualifier, *++argv, MAX_QUALIFIER-2);
            if (strlen(*argv) >= MAX_QUALIFIER-2) {
              printf("ERROR: MAX_QUALIFIER too small for %s\n", *argv);
            }
          }
          pCommon->skipFileCheck = 0;
          pCommon->noGSC = 1;
          break;

        case 'h': /* start disk */
          pCommon->pauseMode = 1;
          argc--;
          nvals = sscanf(*++argv, "%d", &startDisk);
          if (nvals != 1) {
            printf("ERROR: Can not decode startDisk\n");
            errorFlag = 1;
          }
          break;

        case 'H': /* SQL database user/host */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            if (strlen(*++argv) < ((MAX_BUFFER/2)-2)) {
              strcpy(tempusername, *argv);
              atPtr = strchr(tempusername, '@');
              if (atPtr == NULL) {
                strcpy(pCommon->mysqlhost, *argv);
              } else {
                *atPtr = 0;
                atPtr++;
                strcpy(pCommon->mysqlhost, atPtr);
                strcpy(pCommon->username, tempusername);
              }
            } else {
              printf("ERROR argument is too long %s\n", *argv);
              errorFlag = 1;
            }
          }
          break;


        case 'r': /* enable requeueing */
        case 'R':
          pCommon->enableRequeue = 1;
          break;

        default:
          printf("ERROR: * illegal command -%c-", cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  /* Validate arguments */

  if (pCommon->runSystemCommand) {
    pCommon->buildMosaics = 0;
    pCommon->reprocessTiles = 0;
    pCommon->runCalibration = 0;
    pCommon->runFindAstrometry = 0;
    pCommon->runAstrometryWCS1 = 0;
    pCommon->runFull1 = 0;
    pCommon->runFull2 = 0;
    pCommon->multipleMode = 0;
  } else {
    if (pCommon->listname[0] == 0 && pCommon->enableRequeue != 0) {
      printf("ERROR the requeue option '-r' requires a plate list\n");
      errorFlag = 1;
    }

    if (pCommon->startDateString[0] == 0 && pCommon->listname[0] == 0) {
      printf("ERROR: Neither -d <date> or plate list specified\n");
      errorFlag = 1;
    }

    if (pCommon->startDateString[0] != 0 && pCommon->listname[0] != 0) {
      printf("ERROR: Both -d %s or plate list %s can not be specified\n", pCommon->startDateString, pCommon->listname);
      errorFlag = 1;
    }

    if (pCommon->endDateString[0] != 0 && pCommon->startDateString[0] == 0) {
      printf("ERROR: End date %s specified without start date\n", pCommon->endDateString);
      errorFlag = 1;
    }

    if (pCommon->startDateString[0] != 0) {
      if (pCommon->endDateString[0] == 0) {
        strcpy(pCommon->endDateString, "2200");
      }

      if (strcmp(pCommon->startDateString, pCommon->endDateString) >= 0) {
        printf("ERROR: end time %s is less than or equal to start time %s\n", pCommon->endDateString, pCommon->startDateString);
        errorFlag = 1;
      }
    }

    if (pCommon->buildMosaics == 0 &&
        pCommon->runCalibration == 0 &&
        pCommon->runFindAstrometry == 0 &&
        pCommon->runAstrometryWCS1 == 0 &&
        pCommon->runFull1 == 0 &&
        pCommon->runAstrometryWCS2 == 0 &&
        pCommon->runFull2 == 0) {
      printf("ERROR: Nothing to do buildMosaics %d runCalibration %d find_astrometry: %d AstrometryWCS: %d %d runFull %d %d \n",
             pCommon->buildMosaics,
             pCommon->runCalibration,
             pCommon->runFindAstrometry,
             pCommon->runAstrometryWCS1,
             pCommon->runAstrometryWCS2,
             pCommon->runFull1,
             pCommon->runFull2);
      errorFlag = 1;
    }

    if (pCommon->multipleMode != 0 && pCommon->runAstrometryWCS2 == 0 && pCommon->runFull2 == 0) {
      printf("ERROR: Nothing to do on multiple pass AstrometryWCS2: %d pCommon->runFull2 %d\n", pCommon->runAstrometryWCS2, pCommon->runFull2);
      errorFlag = 1;
    }

    if (pCommon->runCalibration != 0 && pCommon->listname[0] != 0) {
      printf("ERROR: calibration mode does not use a list, but requires a starting date\n");
      errorFlag = 1;
    }

    if (pCommon->buildMosaics != 0 && pCommon->listname[0] != 0) {
      printf("ERROR: mosaic creation does not use a list, but requires a starting date\n");
      errorFlag = 1;
    }

    if (pCommon->reprocessTiles != 0 && (pCommon->suffix[0] == 0 || pCommon->tileLocation <= 0)) {
      printf("ERROR: reprocessing requires both -E and -i\n");
      errorFlag = 1;
    }
  }

  if (errorFlag) {
    printf("Usage: runpipeline [<plate list>| -d <startDate>] [-e <endDate>] [-c] [-a] [-x]  -p <subprocess count>\n");
    printf("  where -d creates the plate list from the MySql database (first pass only)\n");
    printf("        -e end date for the list\n");
    printf("        -E <suffix> append this suffix to the linearity file\n");
    printf("        -a run find_astrometry.csh (first pass only)\n");
    printf("        -w run AstrometryWCS (first pass only)\n");
    printf("        -W  <Database password>\n");
    printf("        -x Do not run full.csh (first pass only\n");
    printf("        -c Process calibration frames (flats and linearity)\n");
    printf("        -m <mosaicLocation> Build mosaics on specified disk\n");
    printf("        -r Enable requeing of a job for Odyssey\n");
    printf("        -s <solutionNumber\n");
    printf("        -q <file name qualifier> Add the catalog qualifier to output file names, and check only output from this calibration\n");
    printf("        -v verbose mode\n");
    printf("        -h <startDisk> Halt until /dasch/raid<startDisk>/ExposureData/start.txt is written into the scripts directory\n");
    printf("        -H <Database host>\n");
    printf("        -i <tileLocation> Find tiles and calibration on specified disk\n");
    printf("        -2 w|[f|d] run iteratively, performing the following actions on each successive pass:\n");
    printf("           w run AstrometryWCS\n");
    printf("           f run full.csh\n");
    printf("           d run fulldb.csh\n");
    printf("        -g <command string> Execute a system command on all processors\n");
    exit(-1);
  }

  if (pCommon->listname[0] == 0) {
    time (&pCommon->beginTime);
    ptr = localtime(&pCommon->beginTime);
    strftime(pCommon->runtag, 25, "%Y-%m-%dT%H-%M-%S", ptr);
  } else {
    charPtr2 = strrchr(pCommon->listname, '/');
    if (charPtr2 == NULL) {
      strcpy(pCommon->runtag, pCommon->listname);
    } else {
      charPtr2++;
      strcpy(pCommon->runtag, charPtr2);
    }

    charPtr2 = strchr(pCommon->runtag, '.');
    if (charPtr2 != NULL) {
      *charPtr2 = 0;
    }
  }

  /* Initialize MPI */

  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &pCommon->myrank);
  MPI_Get_processor_name(pCommon->processorName, &pCommon->namelength);
  if ((charPtr = strstr(pCommon->processorName, ".rc.fas.harvard.edu")) != NULL) {
    *charPtr = 0;
    pCommon->namelength -= 19;
  }

  /* At this point, we will delay until the appearance of the "start.txt" file indicates that the file transfer is complete */
  if (pCommon->pauseMode) {
    scriptDir = getenv("DASCH_SCRIPTS");
    if (scriptDir == NULL) {
      printf("ERROR: DASCH_SCRIPTS is not defined\n");
      return 0;
    }

    sprintf(pauseFile, "/dasch/raid%03d/ExposureData/start.txt", startDisk);

    while (1) {
      result = stat(pauseFile, &statbuf);
      if (result == 0) {
        break;
      }

      if (startTime == 0 && pCommon->myrank == 0) {
        printf("runpipeline pausing for file transfer, looking for %s\n", pauseFile);
        time(&startTime);
      }

      sleep(10);
    }

    if (startTime != 0 && pCommon->myrank == 0) {
      time(&curTime);
      curTime = curTime - startTime;
      printf("runpipeline resuming after %ld seconds\n", curTime);
    }
  }

  time (&pCommon->beginTime);
  ptr = localtime(&pCommon->beginTime);
  strftime(pCommon->timestr, 25, "%Y-%m-%dT%H-%M-%S", ptr);

  if (pCommon->myrank == 0) {
    master(pCommon);
  } else {
    slave(pCommon);
  }

  /* Shut down MPI */
  MPI_Finalize();
  return 0;
}


int
master(PRUNCOMMON pCommon)
{
  FILE * listFile;
  int errorFlag = 0;
  char *inBuffer;
  char fittedFile[MAX_BUFFER];
  char inLine[MAX_BUFFER];
  int lineLen;
  int nlinesmax = 0;
  char *scriptdirectory;
  FILE * keepalivefile;
  char logfilename[MAX_BUFFER];
  int curline = 0;
  time_t curTime;
  time_t astrometryTime;
  time_t runFullTime;
  double averageTime;
  int haveFittedList = 0;
  int queueIndex;
  PWORKQUEUE pWorkEntry;
  int numpending = 0;
  int mpiresult;
  char receive_buffer[MAX_MESSAGE_SIZE];
  MPI_Status status;
  int mpitestflag;
  MPI_Request request;
  PNODESTATS pdiemessage;
  int numstarted = 0;
  int numprocessed = 0;
  int rank;
  char rxstring[MAX_BUFFER];
  char txstring[MAX_BUFFER];
  char rxdiskstring[MAX_BUFFER];
  char txdiskstring[MAX_BUFFER];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  char series[MAX_SERIES_STRING];
  int count;
  int cpuIndex;
  struct tm *ptr;
  char* logsDir = NULL;

  /* Find out how many processes there are in the default communicator */

  MPI_Comm_size(MPI_COMM_WORLD, &pCommon->maxproc);
  printf("MPI size %d run tag %s\n", pCommon->maxproc, pCommon->runtag);
  if (pCommon->maxproc == 1) {
    pCommon->startrank = 0;
  }

  CleanWorkQueue(pCommon);

  /* Attempt to open the list */

  if (pCommon->listname[0] != 0) {
    listFile = fopen(pCommon->listname, "rt");

    if (listFile == NULL) {
      printf("Could not open file %s\n", pCommon->listname);
      errorFlag = 1;
    } else {
      /* Count the number of records in the file */
      while (1) {
        inBuffer = fgets(inLine, MAX_BUFFER, listFile);
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

        if (lineLen > 5) {
          pCommon->workQueueSize++;
        }
      }

      fclose(listFile);

      if (pCommon->workQueueSize == 0) {
        printf("Error: no lines in the list file\n");
        errorFlag = 1;
      }
    }
  }

  nlinesmax = pCommon->workQueueSize;
  scriptdirectory = getenv("DASCH_SCRIPTS");
  BUFPRINTF(pCommon->keepalivename, "%s/keepalive%s.txt", scriptdirectory, pCommon->runtag);
  keepalivefile = fopen(pCommon->keepalivename, "wt");

  if (keepalivefile == NULL) {
    errorFlag = 1;
    printf("Failed to open keepalive file %s\n", pCommon->keepalivename);
  } else {
    fprintf(keepalivefile, "Delete this file to stop pipeline processing\n");
    fclose(keepalivefile);
    keepalivefile = NULL;
  }

  printf(
    "runpipeline of %s %s, list file %s or scan date %s, entries %d, "
    "calibration %d buildMosaics on %d find_astrometry: %d AstrometryWCS: %d "
    "runFull1 %d runSystemCommand %d pCommon->solutionNumber %d "
    "suffix '%s' host '%s' tileDisk %d\n",
    __DATE__,
    __TIME__,
    pCommon->listname,
    pCommon->startDateString,
    pCommon->workQueueSize,
    pCommon->runCalibration,
    pCommon->mosaicLocation,
    pCommon->runFindAstrometry,
    pCommon->runAstrometryWCS1,
    pCommon->runFull1,
    pCommon->runSystemCommand,
    pCommon->solutionNumber,
    pCommon->suffix,
    pCommon->mysqlhost,
    pCommon->tileLocation
  );
  printf("Size of WORKQUEUE: %lu, NODESTATS: %lu\n", sizeof(WORKQUEUE), sizeof(NODESTATS));

  logsDir = getenv("DASCH_LOGS_ODYSSEY");
  if (logsDir == NULL) {
    sprintf(logfilename, "/dasch/mosaic/ExposureData/Logs/runpipeline_%s.log", pCommon->runtag);
  } else {
    sprintf(logfilename, "%s/runpipeline_%s.log", logsDir, pCommon->runtag);
  }

  if (pCommon->logFileHandle == NULL) {
    pCommon->logFileHandle = fopen(logfilename, "a+t");
    if (pCommon->logFileHandle) {
      fprintf(pCommon->logFileHandle, "# %s\n", pCommon->timestr);
      fflush(pCommon->logFileHandle);
    } else {
      printf("ERROR: Failed to open %s\n", logfilename);
      errorFlag = 1;
    }
  }

  if (pCommon->runCalibration != 0) {
    errorFlag = CreateCalibrationList(pCommon);

    if (errorFlag == 0 && pCommon->workQueueSize > 0) {
      errorFlag = ProcessWorkQueue(pCommon, pCommon->maxproc);

#ifdef SHOW_DISK_STATS
      ShowDiskStats(pCommon->beginTime, 2);
#endif /* SHOW_DISK_STATS */
    } else {
      printf("No calibration entries found\n");
    }
  }

  if (pCommon->buildMosaics != 0) {
    errorFlag = CreateMosaicList(pCommon);

    if (errorFlag == 0 && pCommon->workQueueSize > 0) {
      ptr = localtime(&pCommon->beginTime);
      strftime(pCommon->mosaicDateString, 25, "%Y-%m-%dT%H-%M-%S", ptr);
      errorFlag = ProcessWorkQueue(pCommon, pCommon->maxproc);

#ifdef SHOW_DISK_STATS
      ShowDiskStats(pCommon->beginTime, 2);
#endif /* SHOW_DISK_STATS */
    } else {
      printf("No mosaic entries found\n");
    }
  }

  /* Now allocate space for all of the records */
  if (pCommon->listname[0] != 0) {
    pCommon->pWorkQueue = (PWORKQUEUE) calloc(pCommon->workQueueSize + 1, sizeof(WORKQUEUE));

    if (pCommon->pWorkQueue == NULL) {
      printf("Failed to allocate pCommon->pWorkQueue %llx\n", (unsigned long long) pCommon->pWorkQueue);
      errorFlag = 1;
    } else {
      for (queueIndex = 0; queueIndex < pCommon->workQueueSize; queueIndex++) {
        pWorkEntry = &pCommon->pWorkQueue[queueIndex];
        strcpy(pWorkEntry->plateName, "UNKNOWN");
        pWorkEntry->queueIndex = queueIndex;
        pWorkEntry->skipFlag = 0;
      }
    }

    /* Attempt to open the list */
    listFile = fopen(pCommon->listname, "rt");

    if (listFile == NULL) {
      printf("Could not open file %s\n", pCommon->listname);
      return -1;
    } else {
      /* Read the file records */
      while (1) {
        pWorkEntry = &pCommon->pWorkQueue[curline];
        memset(pWorkEntry, 0, sizeof(WORKQUEUE));

        inBuffer = fgets(inLine, MAX_BUFFER, listFile);
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

        if (lineLen > 5) {
          if (strlen(inBuffer) >= MAX_LIST_STRING - 1) {
            printf("ERROR: MAX_LIST_STRING exceeded by %s\n", inBuffer);
            errorFlag = 1;
          } else {
            if (ParseFilename(inBuffer, series, &plateNumber, &mosaicNumber, &binning, &rotation)) {
              strcpy(pWorkEntry->plateName, inBuffer);
              pWorkEntry->queueIndex = curline;
              pWorkEntry->skipFlag = 0;
              strcpy(pWorkEntry->series, series);
              pWorkEntry->plateNumber = plateNumber;
              pWorkEntry->scanNumber = mosaicNumber;
              curline++;
            } else {
              printf("ERROR: failed to parse %s from fitted.list\n", inBuffer);
            }
          }
        }
      }

      fclose(listFile);
      pCommon->workQueueSize = curline;
    }
  }

  if (pCommon->startDateString[0] != 0) {
    pCommon->databaseFlag1 = 1;
  }

  InitializeRestartState(pCommon);

  if (errorFlag == 0 && pCommon->runFindAstrometry) {
    if (pCommon->startDateString[0] != 0) {
      /* Build our list from the MySQL database */
      if (errorFlag == 0) {
        errorFlag = CreateListFile(pCommon, 0);
      }

      nlinesmax = pCommon->workQueueSize;
    }

    if (errorFlag == 0 && pCommon->workQueueSize > 0) {
      errorFlag = ProcessWorkQueue(pCommon, pCommon->maxproc);

#ifdef SHOW_DISK_STATS
      ShowDiskStats(pCommon->beginTime, 2);
#endif /* SHOW_DISK_STATS */
    } else {
      printf("No entries found for find_astrometry.csh\n");
    }
  }

  if (errorFlag == 0 && pCommon->runSystemCommand) {
    errorFlag = BuildSystemCommandList(pCommon);
    nlinesmax = pCommon->workQueueSize;

    if (errorFlag == 0 && pCommon->workQueueSize > 0) {
      errorFlag = ProcessWorkQueue(pCommon, pCommon->maxproc);
    }
  }

  /* Begin iteration over pCommon->solutionNumber */
  if (errorFlag == 0) {
    while (1) {
      if (errorFlag == 0 && pCommon->runAstrometryWCS1) {
        int tmpmaxproc = (pCommon->maxproc + 1) / 2;

        if (tmpmaxproc == 0) {
          tmpmaxproc = 1;
        }

        time(&astrometryTime);

        for (queueIndex = 0; queueIndex < pCommon->workQueueSize; queueIndex++) {
          pWorkEntry = &pCommon->pWorkQueue[queueIndex];
          pWorkEntry->workType = WORK_RUN_ASTROMETRYWCS;
        }

        strcpy(fittedFile, scriptdirectory);
        strcat(fittedFile, "/fitted.list");
        remove(fittedFile);

        errorFlag = ProcessWorkQueue(pCommon, tmpmaxproc);

        if (errorFlag == 0 && pCommon->canceledFlag == 0) {
          if (CreateFittedList(pCommon) == 0) {
            pCommon->canceledFlag = 1;
          } else {
            haveFittedList = 1;
          }
        }

#ifdef SHOW_DISK_STATS
        ShowDiskStats(pCommon->beginTime, 3);
#endif /* SHOW_DISK_STATS */
      }

      if (errorFlag == 0 && pCommon->canceledFlag == 0 && pCommon->runFull1 != 0) {
        if (pCommon->startDateString[0] != 0) {
          if (haveFittedList == 0) {
            /* Build our list from the MySQL database */
            if (errorFlag == 0) {
              errorFlag = CreateListFile(pCommon, 1);
            }

            pCommon->databaseFlag1 = 1;
          }
        }

        if (errorFlag == 0 && pCommon->workQueueSize > 0) {
          time(&runFullTime);

          for (queueIndex = 0; queueIndex < pCommon->workQueueSize; queueIndex++) {
            pWorkEntry = &pCommon->pWorkQueue[queueIndex];
            if (pCommon->databaseFlag1 == 0) {
              pWorkEntry->workType =  WORK_RUN_FULL_CSH;
            } else {
              pWorkEntry->workType =  WORK_RUN_FULLDB_CSH;
            }
          }
          errorFlag = ProcessWorkQueue(pCommon, pCommon->maxproc);
        } else {
          printf("No entries found for full.csh\n");
        }

        if (errorFlag == 0 && pCommon->canceledFlag == 0) {
          if (CheckFullProgress(pCommon)) {
            pCommon->canceledFlag = 1;
          }
        }

        time(&curTime);
        curTime = curTime - runFullTime;

        if (pCommon->completeCount == 0) {
          averageTime = 0.0;
        } else {
          averageTime = (1.0 * curTime) / pCommon->completeCount;
        }

        if (pCommon->completeCountMax < pCommon->completeCount) {
          pCommon->completeCountMax = pCommon->completeCount;
        }

#ifdef SHOW_DISK_STATS
        ShowDiskStats(pCommon->beginTime, 4);
#endif /* SHOW_DISK_STATS */

        printf(
          "full.csh time: %ld seconds (%.1f sec/plate) for %d mosaics "
          "completed out of %d mosaics %d processors\n",
          curTime,
          averageTime,
          pCommon->completeCount,
          pCommon->workQueueSize,
          pCommon->maxproc
        );
      }

      /* End of the loop */
      if (pCommon->multipleMode == 0 || pCommon->canceledFlag == 1) {
        printf(
          "Multiple Loop mode exiting with mode %d canceled %d at solution %d\n",
          pCommon->multipleMode,
          pCommon->canceledFlag,
          pCommon->solutionNumber
        );
        break;
      }

      /* Go on to the next solution */
      pCommon->solutionNumber++;
      pCommon->runAstrometryWCS1 = pCommon->runAstrometryWCS2;
      pCommon->runFull1 = pCommon->runFull2;
      pCommon->databaseFlag1 = pCommon->databaseFlag2;
      printf(
        "runpipeline solution %d with AstrometryWCS %d, runFull %d, "
        "pCommon->databaseFlag %d pCommon->workQueueSize %d\n",
        pCommon->solutionNumber,
        pCommon->runAstrometryWCS1,
        pCommon->runFull1,
        pCommon->databaseFlag1,
        pCommon->workQueueSize
      );
    }
  }

  /* All tasks are complete, send out the request to shut down */
  /* Tell all the slaves to exit by sending an empty message with the DIETAG. */

  for (cpuIndex = 0; cpuIndex < pCommon->remoteCPUTotal; cpuIndex++) {
    printf("CPU %8s Completed %5d\n", pCommon->remoteProcessorName[cpuIndex], pCommon->remotePlatesCompleted[cpuIndex]);
    pCommon->remotePlatesCompleted[cpuIndex] = 0;
  }

  pCommon->statusRequestNumber++;

  for (rank = pCommon->startrank; rank < pCommon->maxproc; ++rank) {
    mpiresult = MPI_Send(0, 0, MPI_BYTE, rank, DIETAG, MPI_COMM_WORLD);
    numpending++;

    if (mpiresult != MPI_SUCCESS) {
      printf("ERROR: Master sending DIETAG to %d with result %d\n", pCommon->maxproc, mpiresult);
    }
  }

  while (numpending > 0) {
    /* Receive results from a slave */
    sleep(1);
    mpiresult = MPI_Irecv(
      &receive_buffer,  /* message buffer */
      MAX_MESSAGE_SIZE, /* one data item */
      MPI_BYTE,         /* of type double real */
      MPI_ANY_SOURCE,   /* receive from any sender */
      MPI_ANY_TAG,      /* any type of message */
      MPI_COMM_WORLD,   /* default communicator */
      &request
    );

    if (mpiresult == MPI_SUCCESS) {
      while (1) {
        mpiresult = MPI_Test(&request, &mpitestflag, &status);
        if (mpiresult != MPI_SUCCESS) {
          break;
        }

        if (mpitestflag) {
          break;
        }

        sleep(1);
      }
    }

    MPI_Get_count(&status, MPI_BYTE, &count);

    if (mpiresult != MPI_SUCCESS) {
      printf("ERROR: Master receiving tag %d, err %d, count %3d, cancelled %d mpiresult %d from source %d",
             status.MPI_TAG, status.MPI_ERROR, count, status._cancelled, mpiresult, status.MPI_SOURCE);
    } else {
      if (pCommon->verbose) {
        printf("Master receiving tag %d, err %d, count %3d, cancelled %d mpiresult %d from source %d\n",
               status.MPI_TAG, status.MPI_ERROR, count, status._cancelled, mpiresult, status.MPI_SOURCE);
      }
    }

    switch (status.MPI_TAG) {
    case DIETAG:
      numpending--;

      if (count == 0) {
        printf("Master receiving DIETAG\n");
      } else if (sizeof(NODESTATS) == count) {
        pdiemessage = (PNODESTATS) &receive_buffer;
        sprintf(rxstring, "%lld", pdiemessage->rxrate);
        sprintf(txstring, "%lld", pdiemessage->txrate);
        sprintf(rxdiskstring, "%lld", pdiemessage->rxdiskrate);
        sprintf(txdiskstring, "%lld", pdiemessage->txdiskrate);

        for (cpuIndex = 0; cpuIndex < pCommon->remoteCPUTotal; cpuIndex++) {
          if (strcmp(pdiemessage->nodename, pCommon->remoteProcessorName[cpuIndex]) == 0) {
            break;
          }
        }

        if (cpuIndex < MAX_CPUS) {
          if (cpuIndex == pCommon->remoteCPUTotal) {
            strcpy(pCommon->remoteProcessorName[cpuIndex], pdiemessage->nodename);
            pCommon->remoteCPUTotal++;
          }
          pCommon->remotePlatesCompleted[cpuIndex] += pdiemessage->newTaskscompleted;
        } else {
          printf("ERROR: cpuIndex %d greater than MAX_CPUS\n", cpuIndex);
          exit(-1);
        }

        printf(
          "Recstat %d node %3d %s %-9s proc: %4d %4d Totalsec %6d %6d "
          "Elapsec %6d %6d nrx %6s ntx %6s kB drx %6s dtx %6s kB/sec(DIETAG)\n",
          pCommon->statusRequestNumber,
          status.MPI_SOURCE,
          pdiemessage->nodename,
          workTypeString[pdiemessage->workType],
          pdiemessage->taskscompleted,
          pdiemessage->newTaskscompleted,
          pdiemessage->totalSeconds,
          pdiemessage->newTotalSeconds,
          pdiemessage->elapsedSeconds,
          pdiemessage->newElapsedSeconds,
          rxstring,
          txstring,
          rxdiskstring,
          txdiskstring
        );
      } else {
        printf("ERROR: DIETAG with size %d received\n", count);
      }

      /* We are done, having received our own DIETAG */
      break;

    default:
      printf("ERROR: received unknown tag %d\n", status.MPI_TAG);
      break;
    }

    if (pCommon->verbose) {
      printf("numpending %d, numstarted %d, numprocessed %d (2)\n", numpending, numstarted, numprocessed);
    }
  }

  if (pCommon->verbose) {
    printf("numpending %d, numstarted %d, numprocessed %d (3)\n", numpending, numstarted, numprocessed);
  }

  for (cpuIndex = 0; cpuIndex < pCommon->remoteCPUTotal; cpuIndex++) {
    printf("CPU %8s Completed %5d\n", pCommon->remoteProcessorName[cpuIndex], pCommon->remotePlatesCompleted[cpuIndex]);
    pCommon->remotePlatesCompleted[cpuIndex] = 0;
  }

  time(&curTime);
  curTime -= pCommon->beginTime;

  if (pCommon->completeCountMax == 0) {
    averageTime = 0.0;
  } else {
    averageTime = (1.0 * curTime) / pCommon->completeCountMax;
  }

#ifdef SHOW_DISK_STATS
  ShowDiskStats(pCommon->beginTime, 5);
#endif /* SHOW_DISK_STATS */

  printf(
    "Execution Time: %ld seconds (%.1f sec/plate) for %d mosaics "
    "completed out of %d mosaics %d processors\n",
    curTime,
    averageTime,
    pCommon->completeCountMax,
    nlinesmax,
    pCommon->maxproc
  );

  if (pCommon->pWorkQueue != NULL) {
    free(pCommon->pWorkQueue);
  }

  if (pCommon->logFileHandle) {
    fclose(pCommon->logFileHandle);
  }

  return EXIT_SUCCESS;
}


void
GetNetworkCounts(long long *prxbytes, long long *ptxbytes)
{
  char filename[] = "/proc/net/dev";
  FILE *filehandle = NULL;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  int lineNumber;
  long long rxbytes;
  long long txbytes;
  long long eth_rxbytes = 0;
  long long eth_txbytes = 0;
  long long bond_rxbytes = 0;
  long long bond_txbytes = 0;
  long long ib_rxbytes = 0;
  long long ib_txbytes = 0;
  long long skipval;
  char *charPtr;
  int nvals;
  int verbose = 0;

  *prxbytes = 0;
  *ptxbytes = 0;

  filehandle = fopen(filename, "rt");
  if (filehandle == NULL) {
    printf("ERROR opening %s\n", filename);
    return;
  }

  lineNumber = 0;

  while (1) {
    inBuffer = fgets(inLine, MAX_BUFFER, filehandle);
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

    lineNumber++;

    if ((charPtr = strstr(inBuffer, ":")) != NULL) {
      if (verbose) {
        printf("%s\n", inBuffer);
      }

      *charPtr = 0;
      charPtr++;

      nvals = sscanf(
        charPtr,
        "%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld",
        &rxbytes, &skipval, &skipval, &skipval, &skipval, &skipval, &skipval, &skipval,
        &txbytes, &skipval, &skipval, &skipval, &skipval, &skipval, &skipval, &skipval
      );

      if (nvals != 16) {
        printf("ERROR: nvals is %d in line %s\n", nvals, charPtr);
      } else {
        if (verbose) {
          printf("Link %s rx %6lld tx %6lld\n", inBuffer, rxbytes, txbytes);
        }

        if (strstr(inBuffer, "eth") != 0 || verbose != 0) {
          eth_rxbytes += rxbytes;
          eth_txbytes += txbytes;
        } else if (strstr(inBuffer, "bond") != 0 || verbose != 0) {
          bond_rxbytes += rxbytes;
          bond_txbytes += txbytes;
        } else if (strstr(inBuffer, "ib") != 0 || verbose != 0) {
          ib_rxbytes += rxbytes;
          ib_txbytes += txbytes;
        }
      }
    }
  }

  if (bond_rxbytes + bond_txbytes != 0) {
    *prxbytes = bond_rxbytes + ib_rxbytes;
    *ptxbytes = bond_txbytes + ib_txbytes;
  } else {
    *prxbytes = eth_rxbytes + ib_rxbytes;
    *ptxbytes = eth_txbytes + ib_txbytes;
  }

  fclose(filehandle);
}


void
GetDiskCounts(long long *prxbytes, long long *ptxbytes)
{
  char filename[] = "/proc/diskstats";
  FILE *filehandle = NULL;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  int lineNumber;
  long long rxbytes;
  long long txbytes;
  long long skipval;
  char *charPtr;
  int nvals;
  int verbose = 0;

  *prxbytes = 0;
  *ptxbytes = 0;

  filehandle = fopen(filename, "rt");
  if (filehandle == NULL) {
    printf("ERROR opening %s\n", filename);
    return;
  }

  lineNumber = 0;

  while (1) {
    inBuffer = fgets(inLine, MAX_BUFFER, filehandle);
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

    lineNumber++;

    if ((charPtr = strstr(inBuffer, "sdc1")) != NULL || (charPtr = strstr(inBuffer, "sdb")) != NULL) {
      if (verbose) {
        printf("%s\n", inBuffer);
      }

      charPtr += 4;
      nvals = sscanf(
        charPtr,
        "%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld",
        &skipval, &skipval, &rxbytes, &skipval,
        &skipval, &skipval, &txbytes, &skipval,
        &skipval, &skipval, &skipval
      );

      if (nvals != 11) {
        printf("ERROR: nvals is %d in line %s\n", nvals, charPtr);
      } else {
        rxbytes = rxbytes * 1024L;
        txbytes = txbytes * 1024L;

        if (verbose) {
          printf("Disk rx %6lld tx %6lld\n", rxbytes, txbytes);
        }

        *prxbytes += rxbytes;
        *ptxbytes += txbytes;

        if (verbose) {
          printf("\n");
        }
      }
    }
  }

  fclose(filehandle);
}


int
slave(PRUNCOMMON pCommon)
{
  WORKQUEUE workItem;
  PWORKQUEUE pWorkEntry = &workItem;
  MPI_Status status;
  int mpitestflag;
  MPI_Request request;
  int mpiresult;
  NODESTATS diemessage;
  time_t startTime;
  time_t newStartTime;
  time_t curTime;
  int totalSeconds = 0;
  int newTotalSeconds = 0;
  long long oldrxbytes;
  long long oldtxbytes;
  long long rxbytes;
  long long txbytes;
  long long denominator;
  long long oldrxdiskbytes;
  long long oldtxdiskbytes;
  long long rxdiskbytes;
  long long txdiskbytes;

  time(&startTime);
  time(&newStartTime);

  GetNetworkCounts(&oldrxbytes, &oldtxbytes);
  GetDiskCounts(&oldrxdiskbytes, &oldtxdiskbytes);

  memset(&diemessage, 0, sizeof(NODESTATS));
  strcpy(diemessage.nodename, pCommon->processorName);

  if (pCommon->verbose) {
    printf("Slave %d running on %s\n", pCommon->myrank, diemessage.nodename);
  }

  while (1) {
    /* Receive a message from the master */
    sleep(1);

    mpiresult = MPI_Irecv(
      pWorkEntry,        /* message buffer */
      sizeof(WORKQUEUE), /* one data item */
      MPI_BYTE,          /* of type double real */
      0,                 /* receive from any sender */
      MPI_ANY_TAG,       /* any type of message */
      MPI_COMM_WORLD,    /* default communicator */
      &request
    );

    if (mpiresult == MPI_SUCCESS) {
      while (1) {
        mpiresult = MPI_Test(&request, &mpitestflag, &status);

        if (mpiresult != MPI_SUCCESS) {
          break;
        }

        if (mpitestflag) {
          break;
        }

        sleep(1);
      }
    }

    if (status.MPI_TAG == DIETAG) {
      if (pCommon->verbose) {
        printf("received DIETAG: Slave %d\n", pCommon->myrank);
      }

      time(&curTime);

      diemessage.elapsedSeconds = curTime - startTime;
      diemessage.totalSeconds = totalSeconds;
      diemessage.newElapsedSeconds = curTime - newStartTime;
      diemessage.newTotalSeconds = newTotalSeconds;

      GetNetworkCounts(&rxbytes, &txbytes);
      GetDiskCounts(&rxdiskbytes, &txdiskbytes);

      if (diemessage.newElapsedSeconds == 0) {
        diemessage.rxrate = 0;
        diemessage.txrate = 0;
        diemessage.rxdiskrate = 0;
        diemessage.txdiskrate = 0;
      } else {
        denominator = NETWORK_MULTIPLIER * diemessage.newElapsedSeconds;
        diemessage.rxrate = (rxbytes - oldrxbytes) / denominator;
        diemessage.txrate = (txbytes - oldtxbytes) / denominator;
        diemessage.rxdiskrate = (rxdiskbytes - oldrxdiskbytes) / denominator;
        diemessage.txdiskrate = (txdiskbytes - oldtxdiskbytes) / denominator;
      }

      oldrxbytes = rxbytes;
      oldtxbytes = txbytes;
      oldrxdiskbytes = rxdiskbytes;
      oldtxdiskbytes = txdiskbytes;

      mpiresult = MPI_Send(&diemessage, sizeof(diemessage), MPI_BYTE, 0, DIETAG, MPI_COMM_WORLD);
      if (mpiresult != MPI_SUCCESS) {
        printf("ERROR: Slave returning DIETAG from %d with result %d\n", pCommon->myrank, mpiresult);
      }

      return 0;
    } else if (status.MPI_TAG == STATUSTAG) {
      if (pCommon->verbose) {
        printf("received STATUSTAG: Slave %d\n", pCommon->myrank);
      }

      time(&curTime);

      diemessage.elapsedSeconds = curTime - startTime;
      diemessage.newElapsedSeconds = curTime - newStartTime;
      diemessage.totalSeconds = totalSeconds;
      diemessage.newTotalSeconds = newTotalSeconds;

      GetNetworkCounts(&rxbytes, &txbytes);
      GetDiskCounts(&rxdiskbytes, &txdiskbytes);

      if (diemessage.newElapsedSeconds == 0) {
        diemessage.rxrate = 0;
        diemessage.txrate = 0;
        diemessage.rxdiskrate = 0;
        diemessage.txdiskrate = 0;
      } else {
        denominator = NETWORK_MULTIPLIER * diemessage.newElapsedSeconds;
        diemessage.rxrate = (rxbytes - oldrxbytes) / denominator;
        diemessage.txrate = (txbytes - oldtxbytes) / denominator;
        diemessage.rxdiskrate = (rxdiskbytes - oldrxdiskbytes) / denominator;
        diemessage.txdiskrate = (txdiskbytes - oldtxdiskbytes) / denominator;
      }

      oldrxbytes = rxbytes;
      oldtxbytes = txbytes;
      oldrxdiskbytes = rxdiskbytes;
      oldtxdiskbytes = txdiskbytes;

      mpiresult = MPI_Send(&diemessage, sizeof(diemessage), MPI_BYTE, 0, STATUSTAG, MPI_COMM_WORLD);
      if (mpiresult != MPI_SUCCESS) {
        printf("ERROR: Slave returning STATUSTAG from %d with result %d\n", pCommon->myrank, mpiresult);
      }

      time(&newStartTime);

      diemessage.newTaskscompleted = 0;
      diemessage.workType = WORK_NONE;
      newTotalSeconds = 0;
    }  else if (status.MPI_TAG == WORKTAG) {
      diemessage.taskscompleted++;
      diemessage.newTaskscompleted++;

      if (pWorkEntry->workType != WORK_NONE) {
        diemessage.workType = pWorkEntry->workType;
      }

      ChildProcess(pCommon, pWorkEntry, &status, &totalSeconds, &newTotalSeconds);
    } else {
      printf("ERROR: Slave received unknown tag %d\n", status.MPI_TAG);
      exit(-1);
    }
  }
}


void
ChildProcess(PRUNCOMMON pCommon, PWORKQUEUE pWorkEntry, MPI_Status *status, int *pSeconds, int *pNewSeconds)
{
  int mpiresult;
  char *scriptdirectory;
  char priFilename[MAX_BUFFER];
  char appendFilename[MAX_BUFFER];
  char altFilename[MAX_ALTLOGFILES][MAX_BUFFER];
  int altFileIndex;
  static int iteration = 1;
  int fileNumber = -1;
  FILE *logFile = NULL;
  FILE *appendHandle = NULL;
  time_t startTime;
  time_t curTime;
  int statResult;
  struct stat prifilestats;
  struct stat altfilestats[MAX_ALTLOGFILES];
  off_t filesize;
  int errorFlag = 0;
  int result;
  char cmdStr[MAX_BUFFER];
  FILE *tempFile;
  char tempName[MAX_FILENAME];
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  int count;
  char suffixString[2*MAX_FLATFIELDS_SUFFIX];
  char hostString[MAX_BUFFER];
  char passwordString[MAX_BUFFER];

  memset(altFilename, 0, sizeof(altFilename));
  memset(altfilestats, 0, sizeof(altfilestats));
  memset(suffixString, 0, sizeof(suffixString));
  memset(hostString, 0, sizeof(suffixString));
  memset(passwordString, 0, sizeof(suffixString));

  if (pCommon->suffix[0] != 0) {
    /* Note: FitsMosaic requires -E */
    sprintf(suffixString, " -E %s ", pCommon->suffix);
  }

  if (pCommon->mysqlhost[0] != 0) {
    /* Note: FitsMosaic requires -H */
    if (pCommon->username[0] == 0) {
      BUFPRINTF(hostString, " -H %s ", pCommon->mysqlhost);
    } else {
      BUFPRINTF(hostString, " -H %s@%s ", pCommon->username, pCommon->mysqlhost);
    }
  }

  if (pCommon->password[0] != 0) {
    /* Note: FitsMosaic requires -W */
    BUFPRINTF(passwordString, " -W %s ", pCommon->password);
  }

  scriptdirectory = getenv("DASCH_SCRIPTS");
  if (scriptdirectory == NULL) {
    sprintf(pWorkEntry->priLogname, "ERROR in getenv(DASCH_SCRIPTS)\n");
    errorFlag = 1;
  }

  if (errorFlag == 0) {
    BUFPRINTF(
      priFilename,
      "%s/mpi_%s_%03d_%03d_%03d_%s.log",
      scriptdirectory,
      pCommon->processorName,
      pCommon->myrank,
      iteration,
      pWorkEntry->queueIndex,
      pCommon->runtag
    );

    if (strlen(priFilename) > MAX_BUFFER - 1) {
      sprintf(pWorkEntry->priLogname, "ERROR, priFilename exceeds MAX_BUFFER\n");
      errorFlag = 1;
    }
  }

  if (errorFlag == 0) {
    BUFPRINTF(
      appendFilename,
      "%s/mpi_%s_%03d_%03d_%03dbash_%s.log",
      scriptdirectory,
      pCommon->processorName,
      pCommon->myrank,
      iteration,
      pWorkEntry->queueIndex,
      pCommon->runtag
    );

    if (strlen(priFilename) > MAX_BUFFER - 1) {
      sprintf(pWorkEntry->priLogname, "ERROR, appendFilename exceeds MAX_BUFFER\n");
      errorFlag = 1;
    }
  }

  if (errorFlag == 0) {
    fileNumber = open(priFilename, O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU|S_IRWXG);
    if (fileNumber < 0) {
      BUFPRINTF(pWorkEntry->priLogname, "ERROR, Failed to open %s\n", priFilename);
      errorFlag = 1;
    }
  }

  if (errorFlag == 0) {
    logFile = fdopen(fileNumber, "wt");
    if (logFile == NULL) {
      BUFPRINTF(pWorkEntry->priLogname, "ERROR, Failed to open stream for %s\n", priFilename);
    }
  }

  time(&startTime);

  if (errorFlag == 0) {
    fprintf(logFile, "Start transaction for %s\n", priFilename);

    MPI_Get_count(status, MPI_BYTE, &count);

    if (pCommon->verbose) {
      fprintf(
        logFile,
        "receiving tag %d, err %d, count %d, cancelled %d work %3d Slave %d from source %d\n",
        status->MPI_TAG,
        status->MPI_ERROR,
        count,
        status->_cancelled,
        pWorkEntry->queueIndex,
        status->MPI_SOURCE,
        pCommon->myrank
      );
    }
  }

  /* Do the work */
  if (errorFlag == 0) {
    strcpy(tempName, scriptdirectory);
    strcat(tempName, "/");
    strcat(tempName, pWorkEntry->plateName);
    strcat(tempName, ".xxx");

    tempFile = fopen(tempName, "wt");

    if (tempFile == NULL) {
      fprintf(logFile, "ERROR CHILD:  %17s %3d Failed to open file %s\n", pWorkEntry->plateName, pWorkEntry->queueIndex, tempName);
      errorFlag = 1;
    }
  }

  if (errorFlag == 0) {
    fprintf(tempFile, "%s", pWorkEntry->plateName);
    fclose(tempFile);

    switch(pWorkEntry->workType) {

    case WORK_FLATFIELD:
      BUFPRINTF(
        altFilename[0],
        "%s/mpi_%s_%03d_%03d_%03dflatfield_%s.log",
        scriptdirectory,
        pCommon->processorName,
        pCommon->myrank,
        iteration,
        pWorkEntry->queueIndex,
        pCommon->runtag
      );
      BUFPRINTF(
        cmdStr,
        "FitsMedian -m %s %s -n 4 /dasch/raid%03d/ExposureData/FlatFrames/Flat_%s.fit &> %s\n",
        altFilename[0],
        suffixString,
        pWorkEntry->diskLocation,
        pWorkEntry->plateName,
        appendFilename
      );
      break;

    case WORK_LINEARITY:
      BUFPRINTF(
        altFilename[0],
        "%s/mpi_%s_%03d_%03d_%03dlinearity_%s.log",
        scriptdirectory,
        pCommon->processorName,
        pCommon->myrank,
        iteration,
        pWorkEntry->queueIndex,
        pCommon->runtag
      );

      if (pWorkEntry->delay == 0 || pWorkEntry->delay == LINEARITY_DELAY_USEC) {
        BUFPRINTF(
          cmdStr,
          "FitsLinearity -n %s %s %s %s /dasch/raid%03d/ExposureData/FlatStats/%s &> %s\n",
          altFilename[0],
          suffixString,
          hostString,
          passwordString,
          pWorkEntry->diskLocation,
          pWorkEntry->plateName,
          appendFilename
        );
      } else {
        BUFPRINTF(
          cmdStr,
          "FitsLinearity -n %s -d %d %s %s %s /dasch/raid%03d/ExposureData/FlatStats/%s &> %s\n",
          altFilename[0],
          pWorkEntry->delay,
          suffixString,
          hostString,
          passwordString,
          pWorkEntry->diskLocation,
          pWorkEntry->plateName,
          appendFilename
        );
      }
      break;

    case WORK_MOSAIC:
      BUFPRINTF(
        altFilename[0],
        "%s/mpi_%s_%03d_%03d_%03dmosaic_%s.log",
        scriptdirectory,
        pCommon->processorName,
        pCommon->myrank,
        iteration,
        pWorkEntry->queueIndex,
        pCommon->runtag
      );

      BUFPRINTF(
        cmdStr,
        "FitsMosaic -h %s /dasch/raid%03d/ExposureData/Plates/%s/%05d_%02d "
        "-v %d -u -b 1 -o /dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d_01.fit "
        "-n 16 -q /dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d_16.fit "
        "-c -g -l /dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d.log "
        "%s %s %s &> %s\n",
        altFilename[0],
        pWorkEntry->diskLocation, pWorkEntry->series, pWorkEntry->plateNumber, pWorkEntry->scanNumber,
        pCommon->mosaicLocation,
        pCommon->mosaicLocation,
        pWorkEntry->series, pWorkEntry->plateNumber, pWorkEntry->mosaicNumber,
        pWorkEntry->series, pWorkEntry->plateNumber, pWorkEntry->mosaicNumber,
        pCommon->mosaicLocation,
        pWorkEntry->series, pWorkEntry->plateNumber, pWorkEntry->mosaicNumber,
        pWorkEntry->series, pWorkEntry->plateNumber, pWorkEntry->mosaicNumber,
        pWorkEntry->diskLocation,
        pWorkEntry->series, pWorkEntry->plateNumber, pWorkEntry->mosaicNumber,
        pWorkEntry->series, pWorkEntry->plateNumber, pWorkEntry->mosaicNumber,
        suffixString,
        hostString,
        passwordString,
        appendFilename
      );
      break;

    case WORK_FIND_ASTROMETRY:
      BUFPRINTF(
        cmdStr,
        "%s/find_astrometry.csh %s &> %s\n",
        scriptdirectory,
        tempName,
        appendFilename
      );
      break;

    case WORK_RUN_SYSTEM_CMD:
      BUFPRINTF(
        cmdStr,
        "%s >& %s\n",
        pWorkEntry->systemCommandBuffer,
        appendFilename
      );
      break;

    case WORK_RUN_ASTROMETRYWCS:
      BUFPRINTF(
        altFilename[0],
        "%s/mpi_%s_%03d_%03d_%03dwcssql_%s.log",
        scriptdirectory,
        pCommon->processorName,
        pCommon->myrank,
        iteration,
        pWorkEntry->queueIndex,
        pCommon->runtag
      );

      BUFPRINTF(
        altFilename[1],
        "%s/mpi_%s_%03d_%03d_%03dwcslog_%s.log",
        scriptdirectory,
        pCommon->processorName,
        pCommon->myrank,
        iteration,
        pWorkEntry->queueIndex,
        pCommon->runtag
      );

      BUFPRINTF(
        altFilename[2],
        "%s/mpi_%s_%03d_%03d_%03dfitted_%s.log",
        scriptdirectory,
        pCommon->processorName,
        pCommon->myrank,
        iteration,
        pWorkEntry->queueIndex,
        pCommon->runtag
      );

      BUFPRINTF(
        cmdStr,
        "AstrometryWCS -a %s -b %s -c %s -t %s -p 1 -s %s -n %d -m %d -e %d &> %s\n",
        altFilename[0],
        altFilename[1],
        altFilename[2],
        pCommon->timestr,
        pWorkEntry->series,
        pWorkEntry->plateNumber,
        pWorkEntry->scanNumber,
        pWorkEntry->solutionNumber,
        appendFilename
      );
      break;

    case WORK_RUN_FULL_CSH:
      BUFPRINTF(
        cmdStr,
        "%s/full.csh %d %s &> %s",
        scriptdirectory,
        pWorkEntry->solutionNumber,
        tempName,
        appendFilename
      );
      break;

    case WORK_RUN_FULLDB_CSH:
      BUFPRINTF(
        cmdStr,
        "%s/fulldb.csh %d %s &> %s\n",
        scriptdirectory,
        pWorkEntry->solutionNumber,
        tempName,
        appendFilename
      );
      break;
    }

    fprintf(
      logFile,
      "CHILD:  %17s %3d %s Starting %s",
      pWorkEntry->plateName,
      pWorkEntry->queueIndex,
      pCommon->processorName,
      cmdStr
    );

    result = system(cmdStr);
    if (result == -1) {
      printf("ERROR: system command failed with result %d %s for %s\n", errno, strerror(errno), pWorkEntry->plateName);
    }

    /* Now open the appended file name and write the results to the log file */
    appendHandle = fopen(appendFilename, "rt");

    if (appendHandle == NULL) {
      fprintf(logFile, "\nERROR: failed to open append file %s\n", appendFilename);
    } else {
      while (1) {
        inBuffer = fgets(inLine, MAX_BUFFER, appendHandle);
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

        fprintf(logFile, "%s\n", inBuffer);
      }

      fclose(appendHandle);
      unlink(appendFilename);
    }

    time(&curTime);
    curTime = curTime - startTime;

    pWorkEntry->workResult = result;

    if (logFile != NULL) {
      if (WIFSIGNALED(result) && (WTERMSIG(result) == SIGINT || WTERMSIG(result) == SIGQUIT)) {
        fprintf(
          logFile,
          "CHILD:  %17s %3d %s Exited with signal %d at %d seconds\n",
          pWorkEntry->plateName,
          pWorkEntry->queueIndex,
          pCommon->processorName,
          WTERMSIG(result),
          (int) curTime
        );
      }

      if (WIFEXITED(result)) {
        fprintf(
          logFile,
          "CHILD:  %17s %3d %s Exited with status %d at %d seconds\n",
          pWorkEntry->plateName,
          pWorkEntry->queueIndex,
          pCommon->processorName,
          WEXITSTATUS(result),
          (int) curTime
        );
      } else {
        fprintf(
          logFile,
          "CHILD:  %17s %3d %s Exited with unknown result %d at %d seconds\n",
          pWorkEntry->plateName,
          pWorkEntry->queueIndex,
          pCommon->processorName,
          result,
          (int) curTime
        );
      }
    }

    remove(tempName);
  }

  /* Record statistics */
  *pSeconds += curTime;
  *pNewSeconds += curTime;

  if (errorFlag == 0) {
    fflush(logFile);
    fsync(fileNumber);
    fclose(logFile);
    statResult = stat(priFilename, &prifilestats);

    if (statResult != 0) {
      BUFPRINTF(pWorkEntry->priLogname, "ERROR: failed to find file statistics on %s\n", priFilename);
      errorFlag = 1;
    }
  }

  for (altFileIndex = 0; altFileIndex < MAX_ALTLOGFILES; altFileIndex++) {
    if (errorFlag == 0 && altFilename[altFileIndex][0] != 0) {
      statResult = stat(altFilename[altFileIndex], &altfilestats[altFileIndex]);

      if (statResult != 0) {
        /* File was never created, just ignore it */
        altFilename[altFileIndex][0] = 0;
      }
    }
  }

  if (errorFlag == 0) {
    filesize = prifilestats.st_size;
    pWorkEntry->priFilesize = filesize;
    strcpy(pWorkEntry->priLogname, priFilename);

    for (altFileIndex = 0; altFileIndex < MAX_ALTLOGFILES; altFileIndex++) {
      if (altFilename[altFileIndex][0] != 0) {
        if (altfilestats[altFileIndex].st_size >= 0) {
          filesize = altfilestats[altFileIndex].st_size;
          pWorkEntry->altFilesize[altFileIndex] = filesize;
          strcpy(pWorkEntry->altLogname[altFileIndex], altFilename[altFileIndex]);
        } else {
          pWorkEntry->altFilesize[altFileIndex] = 0;
        }
      }
    }
  }

  pWorkEntry->errorFlag = errorFlag;
  pWorkEntry->seconds = curTime;
  iteration++;

  /* Send the result back */

  mpiresult = MPI_Send(
    pWorkEntry,
    sizeof(WORKQUEUE),
    MPI_BYTE,
    0,
    RESPONSETAG,
    MPI_COMM_WORLD
  );

  if (mpiresult != MPI_SUCCESS) {
    printf(
      "ERROR: sending work item %3d with result %d workResult %d Slave %d \n",
      pWorkEntry->queueIndex,
      mpiresult,
      pWorkEntry->workResult,
      pCommon->myrank
    );
  }
}


void
ReprocessCleanup(PRUNCOMMON pCommon, PWORKQUEUE pWorkEntry)
{
  int res;
  MYSQL my_connection;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char *mysqlhost;
  char *username;
  char *password;
  char queryString[MAX_QUERY_STRING];
  PMOSAIC MosaicList = NULL;
  PMOSAIC pMosaic;
  PMOSAIC pGoodMosaic;
  int mosaicIndex;
  int allocCount = 0;
  int mosaicCount = 0;
  int nvals;
  int goodMosaicEntry = -1;
  char mosaicComment[2*MAX_COMMENT_STRING];

  if (pWorkEntry->workType != WORK_MOSAIC) {
    return;
  }

  if (strcmp(pCommon->suffix, "a") != 0) {
    printf("ERROR: ReprocessCleanup can not handle suffix %s\n", pCommon->suffix);
    return;
  }

  if (pCommon->mysqlhost[0] != 0) {
    mysqlhost = &pCommon->mysqlhost[0];
  } else {
    mysqlhost = getenv("DASCH_MYSQLHOST");
  }

  if (mysqlhost == NULL) {
    fprintf(stderr, "DASCH_MYSQLHOST is not defined\n");
    return;
  }

  if (pCommon->username[0] != 0) {
    username = &pCommon->username[0];
  } else {
    username = getenv("DASCH_USERNAME");
  }

  if (username == NULL) {
    fprintf(stderr, "DASCH_USERNAME is not defined\n");
    return;
  }

  if (pCommon->password[0] != 0) {
    password = &pCommon->password[0];
  } else {
    password = getenv("DASCH_PASSWORD");
  }

  if (password == NULL) {
    fprintf(stderr, "DASCH_PASSWORD is not defined\n");
    return;
  }

  mysql_init(&my_connection);

  if (mysql_real_connect(&my_connection, mysqlhost, username, password, "scanner", 0, NULL, CLIENT_FOUND_ROWS)) {
    sprintf(
      queryString,
      "SELECT mosaicNumber, diskLocation, mosaicComment, rotation, WCSSource+0, FitWCS+0, scanNumber "
      "from mosaics where solutionNumber = 0 and series = '%s' and plateNumber = %d "
      "and ((mosaicComment IS NULL) or (not(mosaicComment regexp('deleted')))) ;",
      pWorkEntry->series,
      pWorkEntry->plateNumber
    );

    if (strlen(queryString) > MAX_QUERY_STRING) {
      fprintf(stderr, "ERROR: MAX_QUERY_STRING exceeded %lu\n", strlen(queryString));
      return;
    }

    res = mysql_query(&my_connection, queryString);

    if (!res) {
      res_ptr = mysql_store_result(&my_connection);

      if (res_ptr) {
        allocCount = (int) mysql_num_rows(res_ptr);

        if (allocCount > 0) {
          MosaicList = (PMOSAIC) calloc(allocCount + 2, sizeof(MOSAIC));

          if (MosaicList == NULL) {
            fprintf(stderr, "ERROR: failed to allocate MosaicList within ReprocessCleanup\n");
            return;
          }
        }

        mosaicCount = 0;

        while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
          pMosaic = &MosaicList[mosaicCount];
          memset(pMosaic, 0, sizeof(MOSAIC));

          if (sqlrow[0] == 0) {
            continue;
          } else {
            nvals = sscanf(sqlrow[0], "%d", &pMosaic->mosaicNumber);
            if (nvals != 1) {
              continue;
            }
          }

          if (sqlrow[1] == 0) {
            continue;
          } else {
            nvals = sscanf(sqlrow[1], "%d", &pMosaic->diskLocation);
            if (nvals != 1) {
            }
          }

          if (sqlrow[2] != NULL) {
            strncpy(pMosaic->mosaicComment, sqlrow[2], MAX_COMMENT_STRING-1);
            pMosaic->mosaicComment[MAX_COMMENT_STRING-1] = 0;
            if (strstr(pMosaic->mosaicComment, "RowNoiseFix") != NULL) {
              goodMosaicEntry = mosaicCount;
            }
          }

          if (sqlrow[3] == 0) {
            pMosaic->rotation = 0;
          } else {
            nvals = sscanf(sqlrow[3], "%d", &pMosaic->rotation);
            if (nvals != 1) {
              if (goodMosaicEntry == mosaicCount) {
                goodMosaicEntry = -1;
              }
              continue;
            }
          }

          if (sqlrow[4] == 0) {
            if (goodMosaicEntry == mosaicCount) {
              goodMosaicEntry = -1;
            }
            continue;
          } else {
            nvals = sscanf(sqlrow[4], "%d", &pMosaic->WCSSource);
            if (nvals != 1) {
              if (goodMosaicEntry == mosaicCount) {
                goodMosaicEntry = -1;
              }
              continue;
            }
          }

          if (sqlrow[5] == 0) {
            if (goodMosaicEntry == mosaicCount) {
              goodMosaicEntry = -1;
            }
            continue;
          } else {
            nvals = sscanf(sqlrow[5], "%d", &pMosaic->FitWCS);
            if (nvals != 1) {
              if (goodMosaicEntry == mosaicCount) {
                goodMosaicEntry = -1;
              }
              continue;
            }
          }

          if (sqlrow[6] == 0) {
            if (goodMosaicEntry == mosaicCount) {
              goodMosaicEntry = -1;
            }
            continue;
          } else {
            nvals = sscanf(sqlrow[6], "%d", &pMosaic->scanNumber);
            if (nvals != 1) {
              if (goodMosaicEntry == mosaicCount) {
                goodMosaicEntry = -1;
              }
              continue;
            }
          }

          mosaicCount++;
        }

        mysql_free_result(res_ptr);
      }
    } else {
      printf("Select error %d: %s res %d\n", mysql_errno(&my_connection), mysql_error(&my_connection), res);
    }

    mysql_close(&my_connection);
  } else {
    fprintf(stderr, "Connection failed\n");

    if (mysql_errno(&my_connection)) {
      fprintf(stderr, "Connection error %d: %s\n", mysql_errno(&my_connection), mysql_error(&my_connection));
    }
  }

  if (goodMosaicEntry > 0) {
    pGoodMosaic = &MosaicList[goodMosaicEntry];

    for (mosaicIndex = 0; mosaicIndex < mosaicCount; mosaicIndex++) {
      pMosaic = &MosaicList[mosaicIndex];

      if (mosaicIndex != goodMosaicEntry && pMosaic->scanNumber == pGoodMosaic->scanNumber) {
        if (pMosaic->mosaicComment[0] == 0) {
          strcpy(mosaicComment, "Deleted for RowNoiseFix");
        } else {
          strcpy(mosaicComment, pMosaic->mosaicComment);
          strcat(mosaicComment, ", Deleted for RowNoiseFix");
        }

        if (strlen(mosaicComment) > MAX_COMMENT_STRING - 1) {
          printf(
            "ERROR: mosaicComment is too long for %s%05d_%02d\n",
            pWorkEntry->series,
            pWorkEntry->plateNumber,
            pMosaic->mosaicNumber
          );
        }

        printf(
          "UPDATE mosaics set mosaicComment = '%s' "
          "WHERE series = '%s' and plateNumber = %d and mosaicNumber = %d;\n",
          mosaicComment,
          pWorkEntry->series,
          pWorkEntry->plateNumber,
          pMosaic->mosaicNumber
        );

        if ((pMosaic->FitWCS & 1) == 1) {
          if (pMosaic->WCSSource < 2) {
            printf(
              "ReprocessCleanup delete /dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d_01.fit\n",
              pMosaic->diskLocation,
              pWorkEntry->series,
              pWorkEntry->plateNumber,
              pMosaic->mosaicNumber,
              pWorkEntry->series,
              pWorkEntry->plateNumber,
              pMosaic->mosaicNumber
            );
          } else if (pMosaic->rotation == 0) {
            printf(
              "ReprocessCleanup delete /dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d_01ww_tnx.fit\n",
              pMosaic->diskLocation,
              pWorkEntry->series,
              pWorkEntry->plateNumber,
              pMosaic->mosaicNumber,
              pWorkEntry->series,
              pWorkEntry->plateNumber,
              pMosaic->mosaicNumber
            );
          } else {
            printf(
              "ReprocessCleanup delete /dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d_01r%dww_tnx.fit\n",
              pMosaic->diskLocation,
              pWorkEntry->series,
              pWorkEntry->plateNumber,
              pMosaic->mosaicNumber,
              pWorkEntry->series,
              pWorkEntry->plateNumber,
              pMosaic->mosaicNumber,
              pMosaic->rotation
            );
          }
        }
      }
    }
  }

  if (MosaicList != NULL) {
    free(MosaicList);
  }
}

/*
 *   gcc -ggdb -O0  -I/usr/include/mysql  -I/dasch/install/include -I /usr/include/openmpi-i386 -L /dasch/install/lib  -L /usr/lib/openmpi/lib runpipeline.c  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread -lm -lmpi -o runpipeline
 *
 * Typical execution after scanning plates:
 *
 *   mpirun --host dasch3, dasch4, dasch5 -np 37 runpipeline -d 2011-08-10 -e 2011-08-12 -c -m 17 -a -w -2 wd
 *
 * Typical execution from a list
 *
 *   runpipeline good_scans.debug
 *   runpipeline -d 2008-08-11  -a -w -x
 *
 *  Multithreaded execution:
 *
 *   mpirun -np 2 runpipeline /dasch/Pipeline/los8.tmp
 *   mpirun --host dasch3, dasch4, dasch5 -np 28 runpipeline /dasch/Pipeline/los8.tmp
 *   mpirun -np 3 ./runpipeline -d 2008-08-11 -a -w -x -v > & los1.log
 *   mpirun -np 3 runpipeline -d 2011-08-11 -c -a -w  -2 wd >& /dasch/Pipeline/los.log
 *
 * please refer to the "orte-clean" tool for assistance.
 *
 * Multiple exposure execution:
 *
 *   echo "runpipeline -d 2009-11-01 -s 1  -x -w -2 wd" | at now
 *
 * Run a system command
 *
 *   mpirun -np 3 runpipeline -g 'df -B 100000'
 *   mpirun --hostfile /dasch/Pipeline/mpihostalldasch -np 6 runpipeline -g '/dasch/raid021/Pipeline/Pipeline/setenvironment.csh'
 *   echo "mpirun --hostfile /dasch/Pipeline/mpihostalldasch -np 6 runpipeline -g '/dasch/raid021/Pipeline/Pipeline/setenvironment.csh' >& /dasch/raid021/junk/diskutilization.log" | at now
 *
 * Mar  5, 2008  Edward J. Los - Initial version
 * Jul  8, 2008  Edward J. Los - Integrate with astrometry.net
 * Jul 14, 2008  Edward J. Los - Correct allocCount limit test
 * Aug 25, 2008  Edward J. Los - Improve Execution statistics
 * Nov 25, 2008  Edward J. Los - Avoid a problem with duplicate pids
 * Feb 20, 2009  Edward J. Los - Correct summary calculations when WCS fitting is involved.
 * Apr  1, 2009  Edward J. Los - Write list of mosaics being processed to runpipeline.tmp
 * Apr  7, 2009  Edward J. Los - Add timeout support (Found to be not useful since the system disk filled
 *                               up with undeleted files in the /tmp directory
 * Oct  5, 2009  Edward J. Los - Add multiple exposure support
 * Dec 21, 2009  Edward J. Los - For multiple exposures, continue processing only when the allobjects file is present
 * Dec 15, 2010  Edward J. Los - Add disk performance statistics.
 * May 26, 2011  Edward J. Los - Ask AstrometryWCS to run on half the requested processes.
 * Aug  8, 2011  Edward J. Los - Convert to openmpi for odyssey use
 *                               Remove timeout support
 * Aug 19, 2011  Edward J. Los - Create a separate append file to fix out-of-date bash syntax on Odyssey
 *                               Do not overwrite a good work item with a null work item
 *                               Set in the queue index with reading fitted.list
 * Aug 29, 2011  Edward J. Los - Add DT_UNKNOWN to FitsFilter because of problem on Odyssey.
 * Sep 23, 2011  Edward J. Los - Do not make a mosaic that has already been made
 * Dec  7, 2011  Edward J. Los - Do not use the endDate when searching for mosaics since they could just have been created
 * Jan 13, 2012  Edward J. Los - Replace _count with MPI_Get_count
 * Feb 27, 2012  Edward J. Los - Add the "-n" qualifier for future non-GSC multiple exposure processing.
 * Apr  9, 2012  Edward J. Los - Correct network totals
 * Jan 10, 2013  Edward J. Los - Apply end date to mosaics as well.
 * May 10, 2014  Edward J. Los - Correct a bug where FitsMosaic could write into an existing mosaic directory.
 * Nov 18, 2014  Edward J. Los - Correct FitsMosaic "v" parameter to be the target mosaic disk.
 * Dec  1, 2014  Edward J. Los - In CheckFullProgress, check the NoAllobject flag before searching for the file
 *                               and add a 2 minute timeout for this search.  Finally, look also for the find_astrometry2 result
 *                               if we are going to run AstrometryWCS again.
 * Dec  5, 2014  Edward J. Los - Add SKIP_FILE_CHECK in CheckFullProgress  to skip the file checks which do not seem to work
 * Apr 29, 2015  Edward J. Los - Correct FitsMosaic when a plate is scanned twice
 * Sep 17, 2016  Edward J. Los - move start.txt from /dasch/Pipeline to /dasch/raid<startDisk>/ExposureData/Pipeline
 * Feb 28, 2017  Edward J. Los - Make SKIP_FILE_CHECK a variable skipFileCheck, normally '1' if -n is not specified
 *                               Replace -n with -q which takes a qualifer of the file to look for and look for the allobjects file, not the next SExtractor file
 * Jun  1, 2018  Edward J. Los - Add a carriage to "failed to open append file" for easier log searching, and decode a -1 system() command status
 * Jan 29, 2018  Edward J. Los - Increase MAX_CPUS from 10 to 512 for Odyssey shared queues.
 * Mar 12, 2018  Edward J. Los - Support Odyssey job requeueing with the '-r' option
 * Nov 29, 2019  Edward J. Los - Support NoiseRepair effort:
 *                                redefine -E to mean <suffix> which adds a suffix calibration files
 *                                redefine -H to mean <Database host> or <Database host@Database username>
 *                                redefine -W to mean <Database password>
 *                                correct GetMosaicInfo loop to search through MAX_SCAN_MOSAIC_MOSAIC_NUMBER
 * Dec 10, 2019  Edward J. Los - Put the NoiseRepair qualifier to FitsMosaic at the end to avoid problems with post-processing scripts
 *                               If the Mosaic phase is executed, then adjust the startDateString to the start of the job.
 * Dec 14, 2019  Edward J. Los - Add mosaic number to "skipping row noise fix"
 * Jan 10, 2020  Edward J. Los - Add -H and -W to FitsLinearity FitsMosaic
 * May  5, 2020  Edward J. Los - Remove "no valid mosaic build" as a valid reason for skipping the plate.  Apparently, such a case has not happened to date
 * May 21, 2020  Edward J. Los - Correct row noise sql command to set "Deleted" for all solutions of the mosaic
 * Sep  2, 2020  Edward J. Los - Support DASCH_LOGS_ODYSSEY for runpipeline*.log
 */
