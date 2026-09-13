// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* runfpack.c
 *
 * Searches a disk for uncompressed calibration files and compresses them.
 * 
 *  gcc -g -O0  -D_FILE_OFFSET_BITS=64 runfpack.c -o runfpack
 * 
 *  runfpack -d <disk number> [-r]
 *
 *  2011-11-21 Edward J. Los Initial Implementation
 *  2011-11-25 Edward J. Los Add a statistics phase
 *  2020-01-17 Edward J. Los Add -t for a different target disk
 *  2020-08-17 Edward J. Los Correct problem with unitialized "targetDisk" variable.
 *                           Add a check to ensure that "ExposureData" exists on the target disk
 *                           Add checks for zero size files
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <dirent.h>
#define MAX_PATH 512
#define KEEPALIVE_CHECK 20  /* Check keepalive every 20 files */
/* #define DEBUG_FLAG 1 */
typedef struct _common {
  int getFileCounts;
  int abortFlag;
  int unpackFlag;
  int directoryCount;
  int allfileCount;
  int fileCount;
  int plateCount;
  int calibrationCount;
  int rawFlatCount;
  int targetDiskLocation;
  off_t plateInitialSize;
  off_t plateFinalSize;
  off_t calibrationInitialSize;
  off_t calibrationFinalSize;
  char keepalivename[MAX_PATH];
} COMMON, *PCOMMON;
void ProcessDirectory(PCOMMON pCommon,char *rootDirectory)
{
  struct dirent **filelist;
  int numFiles;
  int fileIndex;
  char filename[MAX_PATH];
  struct stat statbuf;
  char cmdBuf[MAX_PATH];
  char cmdBuf2[MAX_PATH];
  char targetFile[MAX_PATH];
  char targetDirectory[MAX_PATH];
  int result;
  int plateFlag;
  off_t initialSize;
  off_t finalSize;
  char* charPtr;
  char* charPtr2;
  int targetLength;
  int firstFileFlag = 1;
  
  if (pCommon->abortFlag != 0) {
    return;
  }
  if (strstr(rootDirectory,"Mosaics")) {
    /* Never process any mosaics */
    return;
  }
  pCommon->directoryCount++;
  numFiles = scandir(rootDirectory,&filelist,0,alphasort);
  if (numFiles > 0) {
    for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
      if (filelist[fileIndex]->d_name[0] == '.') {
        /* skip "." and ".." */
        continue;
      }
      pCommon->allfileCount++;
      strcpy(filename,rootDirectory);
      strcat(filename,"/");
      strcat(filename,filelist[fileIndex]->d_name);
      if (strstr(filelist[fileIndex]->d_name,".")) {
        if (strstr(filename,".fit") == NULL) {
          continue;
        }
        if (strstr(filename,".head.fit")) {
          if (pCommon->getFileCounts == 0) {
            printf("WARNING: Header file in wrong location %s\n",filename);
          }
          continue;
        }
        if (pCommon->unpackFlag) {
           
          if ((charPtr = strstr(filename,".fit.fz")) == NULL) {
            /* not compressed */
            continue;
          }

        } else {
          if (strstr(filename,".fit.fz")) {
            /* already compressed */
            continue;
          }
          if (strcmp(&filename[strlen(filename)-4],".fit") != 0) {
            continue;
          }

        }

        if (strstr(filename,"FlatFrames/Flat_")) {
          /* Do not compress raw flat frames -- these need to be deleted */
          pCommon->rawFlatCount++;
          if (pCommon->getFileCounts == 0) {
            printf("WARNING: Raw flat frame %s\n",filename);
          }
          continue;
        }

        /* We need to compress/uncompress this file */
          
        pCommon->fileCount++;
        if ((strstr(filename,"Plates") != NULL) &&
            (strstr(filename,"mask") == NULL)) {
          plateFlag = 1;
          pCommon->plateCount++;
        } else {
          plateFlag = 0;
          pCommon->calibrationCount++;
        }
        if (pCommon->getFileCounts == 0) {

          result = stat(filename,&statbuf);
          if (result != 0) {
            printf("ERROR: stat failed on %s\n",filename);
            pCommon->abortFlag = 1;
            return;
          }
          initialSize = statbuf.st_size;
          if (initialSize == 0) {
            printf("ERROR: %s has zero bytes\n",filename);
            pCommon->abortFlag = 1;
            return;
          }
            
          if (pCommon->targetDiskLocation > 0) {
            charPtr = strstr(filename,"ExposureData");
            if (charPtr == NULL) {
              printf("ERROR: could not find ExposureData in %s\n",filename);
              pCommon->abortFlag = 1;
              return;
            }
            sprintf(targetFile,"/dasch/raid%03d/%s",pCommon->targetDiskLocation,charPtr);
            charPtr2 = strstr(targetFile,"/ExposureData"); 
            if (charPtr2 == NULL) {
              printf("ERROR: could not find ExposureData in %s\n",targetFile);
              pCommon->abortFlag = 1;
              return;
            }
            charPtr2 += strlen("/ExposureData");
            *charPtr2 = 0;
            result = stat(targetFile,&statbuf);
            if (result != 0) {
              printf("ERROR: stat failed on %s\n",targetFile);
              pCommon->abortFlag = 1;
              return;
            }
            sprintf(targetFile,"/dasch/raid%03d/%s",pCommon->targetDiskLocation,charPtr);
 
            if (pCommon->unpackFlag) {
              targetLength = strlen(targetFile);
              if (strcmp(&targetFile[targetLength -3],".fz") == 0) {
                targetFile[targetLength -3] = 0;
              } else {
                printf("ERROR: could not find '.fz' in %s\n",targetFile);
                pCommon->abortFlag = 1;
                return;
              }
              sprintf(cmdBuf,"funpack -v -D -O %s %s\n",targetFile,filename);
            } else {
              strcat(targetFile,".fz");
              /* NOTE: fpack does not support -v and -D when -S is specified */
              sprintf(cmdBuf,"fpack -S  %s > %s\n",filename,targetFile);
            }
            if (firstFileFlag != 0) {
              firstFileFlag = 0;
              strcpy(targetDirectory,targetFile);
              charPtr = strrchr(targetDirectory,'/');
              if (charPtr == NULL) {
                printf("ERROR: could not find '/' in %s\n",targetDirectory);
                pCommon->abortFlag = 1;
                return;
              }
              *charPtr = 0;
              sprintf(cmdBuf2,"mkdir -p %s\n",targetDirectory);
              system(cmdBuf2);
            }
          } else {
            if (pCommon->unpackFlag) {
              sprintf(cmdBuf,"funpack -v -D %s\n",filename);
            } else {
              sprintf(cmdBuf,"fpack -v -D %s\n",filename);
            }
          }
#ifdef DEBUG_FLAG
          printf("%s",cmdBuf);
#else /* DEBUG_FLAG */
          system(cmdBuf);
          if (pCommon->targetDiskLocation <= 0) {
            if (pCommon->unpackFlag) {
              charPtr = strstr(filename,".fit.fz");
              charPtr += 4;
              *charPtr = 0;
            } else {
              strcat(filename,".fz");
            }
          }
#endif /* DEBUG_FLAG */

          if (pCommon->targetDiskLocation > 0) {
            result = stat(targetFile,&statbuf);
            if (result != 0) {
              printf("ERROR: stat failed on %s\n",targetFile);
              pCommon->abortFlag = 1;
              return;
            }
            if (statbuf.st_size == 0) {
              printf("ERROR: %s has zero bytes\n",targetFile);
              pCommon->abortFlag = 1;
              return;
            }
            if (pCommon->unpackFlag == 0) {
              /* Because we specified -S, we must make our verbose report and then delete the input */
              printf("%s -> %s\n",filename,targetFile);
              result = unlink(filename);
              if (result != 0) {
                printf("ERROR: failed to remove %s\n",filename);
              }
            }
          } else {
            result = stat(filename,&statbuf);
            if (result != 0) {
              printf("ERROR: stat failed on %s\n",filename);
              pCommon->abortFlag = 1;
              return;
            }
            if (statbuf.st_size == 0) {
              printf("ERROR: %s has zero bytes\n",filename);
              pCommon->abortFlag = 1;
              return;
            }
          }
          finalSize = statbuf.st_size;
          if (plateFlag) {
            pCommon->plateInitialSize += initialSize;
            pCommon->plateFinalSize += finalSize;
          } else {
            pCommon->calibrationInitialSize += initialSize;
            pCommon->calibrationFinalSize += finalSize;
          }
        }


        if ((pCommon->fileCount % KEEPALIVE_CHECK) == 0) {
#ifdef DEBUG_FLAG
          printf("Keepalive check at count %d\n",pCommon->fileCount);
#endif /* DEBUG_FLAG */
          result = stat(pCommon->keepalivename,&statbuf);

          if (result != 0) {
            printf("Terminating at file %d because keepalive file %s is missing\n",pCommon->fileCount,pCommon->keepalivename);
            pCommon->abortFlag = 1;
            return;
          } 
            
        }
        

      } else {
        /* Assume this is a directory and recurse into it */
        ProcessDirectory(pCommon,filename);
        if (pCommon->abortFlag) {
          return;
        }
      }
    }

    for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
      free(filelist[fileIndex]);
    }
    free(filelist);

  }
  return;

}


int main(int argc,char *argv[])
{
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  int diskLocation = -1;
  int nvals;
  time_t curTime;
  time_t startTime;
  char directoryName[MAX_PATH];
  COMMON commonTable;
  PCOMMON pCommon = &commonTable;
  char *scriptdirectory;
  FILE * keepalivefile;
  int initialPlateMB;
  int finalPlateMB;
  int initialCalibrationMB;
  int finalCalibrationMB;
  int totalInitialMB;
  int totalFinalMB;
  double plateRatio;
  double calibrationRatio;
  double totalRatio;
  int iteration;
  double rate;
  int unpackFlag = 0;
  int targetDiskLocation = -1;

  memset(pCommon,0,sizeof(COMMON));
  pCommon->targetDiskLocation = -1;

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      printf("ERROR: unknown argument %s \n",argstr);
      errorFlag = 1;
 
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
	
        case 'd': /* disk drive */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&diskLocation);
            if (nvals != 1) {
              printf("ERROR: Unable to decode diskLocation %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 't': /* target disk drive */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&pCommon->targetDiskLocation);
            if (nvals != 1) {
              printf("ERROR: Unable to decode diskLocation %s\n",*argv);
              errorFlag = 1;
            } else {
              if (pCommon->targetDiskLocation <= 0) {
                printf("ERROR: Illegal target diskLocation %s\n",*argv);
                errorFlag = -1;
              }
            }
            targetDiskLocation = pCommon->targetDiskLocation;
          }
          break;

        case 'u': /* unpack */
        case 'U':
          pCommon->unpackFlag = 1;
          unpackFlag = 1;
          break;


        default:
          printf("ERROR: illegal command -%c-",cmdchar);
          errorFlag = 1;
          break;
        }
        
      }
    }
  }
  if (diskLocation < 0) {
    printf("ERROR: Disk location is not specified\n");
    errorFlag = 1;
  }
  if (errorFlag) {
    printf("Usage: runfpack  -d <source disk location>\n");
    printf("                [-t <target disk location>]\n");
    printf("                [-u run funpack rather than fpack]\n");
    exit(-1);
  }
  time(&startTime);
  printf("runfpack of %s %s, diskLocation %d\n",__DATE__,__TIME__,diskLocation);

  scriptdirectory = getenv("DASCH_SCRIPTS");
  strcpy(pCommon->keepalivename,scriptdirectory);
  strcat(pCommon->keepalivename,"/fpack_keepalive.txt");
  keepalivefile = fopen(pCommon->keepalivename,"wt");
  if (keepalivefile == NULL) {
    errorFlag = 1;
    printf("Failed to open keepalive file %s\n",pCommon->keepalivename);
  } else {
    fprintf(keepalivefile,"Delete this file to stop pipeline processing\n");
    fclose(keepalivefile);
    keepalivefile = NULL;
  }



  sprintf(directoryName,"/dasch/raid%03d/ExposureData",diskLocation);

  pCommon->getFileCounts = 1; /* First iteration provides file counts only */
  for (iteration = 0; iteration < 2; iteration++) {
    ProcessDirectory(pCommon,directoryName);
    if (iteration == 0) {
      printf("Disk: %3d AllFiles %d Directories %d rawFlat %d compressed %d plates %d calibration %d\n",
             diskLocation,
             pCommon->allfileCount,
             pCommon->directoryCount,
             pCommon->rawFlatCount,
             pCommon->fileCount,
             pCommon->plateCount,
             pCommon->calibrationCount);
      fflush(stdout);
      memset(pCommon,0,sizeof(COMMON));
      pCommon->unpackFlag = unpackFlag;
      pCommon->targetDiskLocation = targetDiskLocation;

      strcpy(pCommon->keepalivename,scriptdirectory);
      strcat(pCommon->keepalivename,"/fpack_keepalive.txt");
    }

  }

  time(&curTime);
  curTime -= startTime;
  if (curTime > 0) {
    rate = (1.0*pCommon->fileCount)/(1.0*curTime);
  } else {
    rate = 0.0;
  }
  printf("Execution time: %d sec AllFiles %d Directories %d rawFlat %d compressed %d plates %d calibration %d rate %f files/second\n",
         curTime,
         pCommon->allfileCount,
         pCommon->directoryCount,
         pCommon->rawFlatCount,
         pCommon->fileCount,
         pCommon->plateCount,
         pCommon->calibrationCount,
         rate);

  initialPlateMB = pCommon->plateInitialSize/1000000;
  finalPlateMB = pCommon->plateFinalSize/1000000;
  initialCalibrationMB = pCommon->calibrationInitialSize/1000000;
  finalCalibrationMB = pCommon->calibrationFinalSize/1000000;
  totalInitialMB = initialPlateMB + initialCalibrationMB;
  totalFinalMB = finalPlateMB + finalCalibrationMB;
  if (finalPlateMB > 0) {
    plateRatio = (1.0*initialPlateMB)/(1.0*finalPlateMB);
  } else {
    plateRatio = 0;
  }
  if (finalCalibrationMB > 0) {
    calibrationRatio = (1.0*initialCalibrationMB)/(1.0*finalCalibrationMB);
  } else {
    calibrationRatio = 0;
  }
  if (totalFinalMB > 0) {
    totalRatio = (1.0*totalInitialMB)/(1.0*totalFinalMB);
  } else {
    totalRatio = 0;
  }
  printf("Disk: %3d Plates: %d MB -> %d MB %.3fX Calibration: %d MB -> %d MB  %.3fX Total: %d MB -> %d MB %.3fX\n",diskLocation,initialPlateMB,finalPlateMB,plateRatio,initialCalibrationMB,finalCalibrationMB,calibrationRatio,totalInitialMB,totalFinalMB,totalRatio);

  return(0);
}
