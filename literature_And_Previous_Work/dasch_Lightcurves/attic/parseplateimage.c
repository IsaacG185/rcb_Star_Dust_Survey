// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* parseplateimage.c
 *
 *  This program checks the PlateImage*.log files for a variety of errors
 *
 *  Usage:  parseplateimage -l  <list of logs|single log|log directory>  test with run_932.log for i24330
 *
 *           parseplateimage  -v  -l /home/scanner/junk/2020_02_28/2020-02-24/PlateImage2020-02-24T07-24-32.log  -o /home/scanner/2020_02_28/parseplateimage1.db
 *           parseplateimage  -v  -l /home/scanner/junk/2020_02_28/parseplateimage.list  -o /home/scanner/junk/2020_02_28/parseplateimage2.db
 *           parseplateimage  -v  -l /home/scanner/backup  -o /home/scanner/junk/2020_02_28/parseplateimage3.db
 *                                    
 *           the output file, los.db will contain a list of possible saturated tiles (0 MEDIAN)
 *
 *
            ./losyield.csh
            parseplateimage -v -l losyield.txt


The zero median sequence count poduces data like this: for PlateImage2020-02-19T08-04-42.log
0 MEDIAN found in line 464635 at fileNameCount 18300

fileNameCount 20242 20242 lineCount 519154 preceeding count 18300 lineNumber 464635 following count    -1 lineNumber     -1 for 'File Name: Y:\ExposureData\Plates\ai\33833_01\Exposure_000103'
fileNameCount 20243 20243 lineCount 519181 preceeding count 18300 lineNumber 464635 following count    -1 lineNumber     -1 for 'File Name: Y:\ExposureData\Plates\ai\33833_01\Exposure_000104'
fileNameCount 20244 20244 lineCount 519208 preceeding count 18300 lineNumber 464635 following count    -1 lineNumber     -1 for 'File Name: Y:\ExposureData\Plates\ai\33833_01\Exposure_000105'
fileNameCount 20245 20245 lineCount 519235 preceeding count 20247 lineNumber 519286 following count    -1 lineNumber     -1 for 'File Name: Y:\ExposureData\Plates\ai\33833_01\Exposure_000106'  Save this! new preceeding lineNumber
fileNameCount 20246 20246 lineCount 519262 preceeding count 20247 lineNumber 519286 following count 20247 lineNumber 519286 for 'File Name: Y:\ExposureData\Plates\ai\33833_01\Exposure_000107'  Save this!

0 MEDIAN found in line 519286 at fileNameCount 20247

fileNameCount 20247 20247 lineCount 519289 preceeding count 20247 lineNumber 519286 following count    -1 lineNumber     -1 for 'File Name: Y:\ExposureData\Plates\ai\33833_01\Exposure_000108'  Save this 1st following line number
fileNameCount 20248 20248 lineCount 519316 preceeding count 20247 lineNumber 519286 following count    -1 lineNumber     -1 for 'File Name: Y:\ExposureData\Plates\ai\33833_01\Exposure_000109'  Save this 1st following line number
fileNameCount 20249 20249 lineCount 519343 preceeding count 20247 lineNumber 519286 following count    -1 lineNumber     -1 for 'File Name: Y:\ExposureData\Plates\ai\36398_01\Exposure_000110'
fileNameCount 20250 20250 lineCount 519370 preceeding count 20247 lineNumber 519286 following count    -1 lineNumber     -1 for 'File Name: Y:\ExposureData\Plates\ai\36398_01\Exposure_000111'

           
 *
 * 
cc -ggdb -O0  -I/usr/include/mysql  -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm parseplateimage.c pipelineutils.a -L/usr/lib${lib64}/mysql  -lmysqlclient -ltable -lutil -lwcs  -L/usr/lib/mysql  -l mysqlclient -o parseplateimage
 * 
 *
 * Mar  6, 2020 Edward J. Los - Initial version, adopted from parseyield
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
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>

#define MAX_INPUT_NAME 512
#define MAX_BUFFER     1024
#define MAX_TIME_LENGTH 28
#define MAX_SOLUTION_NUMBER 10
#define MAX_PLATE_COUNT 500000
#define MAX_TAIL_LENGTH 36
#define POST_FLOOD "PlateImage2016-07-21"
/* #define DEBUG_PARSE 1 */
#define XY_TOLERANCE 0.0003 /* Must match AerotechQlib.h */
#define Z_TOLERANCE  0.0055 /* Must match AerotechQlib.h */

typedef struct _zeromedian {
  int fileNameCount;  
  int lineCount; /* Line number of this buffer */
  char fileBuffer[MAX_BUFFER];
  int preceedingZeroMedianFileNameCount; /* The last zero median file preceeding this buffer */
  int followingZeroMedianFileNameCount;  /* The first zero median file following this buffer */
  int preceedingZeroMedianLineNumber; /* The last zero median lineNumber preceeding this buffer */
  int followingZeroMedianLineNumber;  /* The first zero median lineNumber following this buffer */
} ZEROMEDIAN,*PZEROMEDIAN;


typedef struct _pparsecommon {
	time_t startTime;
	FILE *out_handle;
	int verbose;
  int abortFlag;
	int maxLineLen;
  int wrongPositionCount;
  int xDeltaCount;
  int yDeltaCount;
  int zDeltaCount;
  int uDeltaCount;
  int axisFaultCount;
  int markerFaultCount;
  int unrecognizedCount;
  int timeoutCount;
  int zeroMedianCount;
  int fileCount;
  int directoryCount;
  int allFileCount;
} PARSECOMMON,*PPARSECOMMON;

int ProcessLogBuffer(PPARSECOMMON pParseCommon, char* log_name) {
	char *inBuffer;
	char copyLine[MAX_BUFFER];
  char prevCopyLine[MAX_BUFFER];
	char inLine[MAX_BUFFER];
	char parseLine[MAX_BUFFER];
	FILE *log_handle = NULL;
	int lineLen;
	int lineCount = 0;
	char *charPtr;
  int wrongPositionCount = 0;
  int axisFaultCount = 0;
  int markerFaultCount = 0;
  int unrecognizedCount = 0;
  int timeoutCount = 0;
  int zeroMedianCount = 0;
  int tailLength;
  char tailName[MAX_TAIL_LENGTH];
  char barcodeSeries[MAX_SERIES_STRING];
  int barcodePlateNumber = 0;
  char noneBarcodeSeries[MAX_SERIES_STRING];
  int noneBarcodePlateNumber = 0;
  char leftBarcodeSeries[MAX_SERIES_STRING];
  int leftBarcodePlateNumber = 0;
  char rightBarcodeSeries[MAX_SERIES_STRING];
  int rightBarcodePlateNumber = 0;
  int lastScanPattern = SCAN_PATTERN_NONE;
  int fileNameTotal = 0;
  int fileNameCount = 0;
  int last0MedianFileNameCount = -1;
  int last0MedianFileLineNumber = -1;
  int last1MedianFileNameCount = -1;
  int last1MedianFileLineNumber = -1;
  int last2MedianFileNameCount = -1;
  int last2MedianFileLineNumber = -1;
  PZEROMEDIAN zeroMedianTable = NULL;
  PZEROMEDIAN pZeroMedianCurrent = NULL;
  PZEROMEDIAN pZeroMedianLast1 = NULL;
  PZEROMEDIAN pZeroMedianLast2 = NULL;
  PZEROMEDIAN pZeroMedian = NULL;
  double xDelta = -1.0;
  double yDelta = -1.0;
  double zDelta = -1.0;
  double uDelta = -1.0;
  int xDeltaCount = 0;
  int yDeltaCount = 0;
  int zDeltaCount = 0;
  int uDeltaCount = 0;
  int nvals;


  barcodeSeries[0] = 0;
  noneBarcodeSeries[0] = 0;
  leftBarcodeSeries[0] = 0;
  rightBarcodeSeries[0] = 0;
  prevCopyLine[0] = 0;
 

  charPtr = strrchr(log_name,'/');
  if (charPtr == NULL) {
    charPtr = log_name;
  } else {
    charPtr++;
  }
  tailLength = strlen(charPtr);
  if (tailLength >=  (MAX_TAIL_LENGTH-1)) {
    printf("ERROR: File Tail length %d is too large in %s\n",tailLength,log_name);
    return(-1);
  }
  strcpy(tailName,charPtr);
  
  if (pParseCommon->verbose != 0) {
    printf("Processing '%s'\n",log_name);
  }
	log_handle = fopen(log_name,"rt");
	if (log_handle == NULL) {
		printf("ERROR: Could not open list file '%s'\n",log_name);
		return(-1);
	}
  pParseCommon->fileCount++;
  /* First count the number of File: entries in this file */
	while (1) {
		inBuffer = fgets(inLine,MAX_BUFFER,log_handle);
		if (inBuffer == NULL) {
			break;
		}
		lineLen = strlen(inBuffer);
		if (lineLen > pParseCommon->maxLineLen) {
			pParseCommon->maxLineLen = lineLen;
		}
    if (strstr(inBuffer,"File Name") != NULL) {
      fileNameTotal++;
    }
  }  

  fclose(log_handle);
  log_handle = NULL;

  zeroMedianTable = (PZEROMEDIAN)calloc(fileNameTotal+10,sizeof(ZEROMEDIAN));
  if (zeroMedianTable == NULL) {
    printf("ERROR: failed to allocate zeroMedianTable of size %d\n",fileNameTotal);
    exit(-1);
  }


  log_handle = fopen(log_name,"rt");

  if (log_handle == NULL) {
    printf("ERROR: Could not open list file (2) '%s'\n",log_name);
    if (zeroMedianTable != NULL) {
      free(zeroMedianTable);
    }
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
    strcpy(prevCopyLine,copyLine);
    strcpy(copyLine,inBuffer);
    if (strstr(inBuffer,"Wrong Position") != NULL) {
      wrongPositionCount++;
      charPtr = strstr(prevCopyLine,"At Desired Position: delta x:");
      if (charPtr == NULL) {
        printf("ERROR: no 'delta' found line %d At lineCount %d %s\n",__LINE__,lineCount-1,prevCopyLine);
      } else {
        charPtr += strlen("At Desired Position: delta x:");
        nvals = sscanf(charPtr,"%lg y: %lg z: %lg, u: %lg",&xDelta,&yDelta,&zDelta,&uDelta);
        if (nvals != 4) {
          printf("ERROR: nvals is %d in '%s' '%s'\n",nvals,copyLine,tailName);
        } else {
          if (xDelta > (XY_TOLERANCE)) {
            xDeltaCount++;
          }
          if (yDelta > (XY_TOLERANCE)) {
            yDeltaCount++;
          }
          if (zDelta > Z_TOLERANCE) {
            zDeltaCount++;
          }

          if (uDelta > (XY_TOLERANCE)) {
            uDeltaCount++;
          }
        }
#if 0
        printf("%s x %d y %d z %d u %d in lineCount %d of %s\n",
               prevCopyLine,
               xDeltaCount,
               yDeltaCount,
               zDeltaCount,
               uDeltaCount,
               lineCount -1,
               log_name);
#endif
      }
    }
    if (strstr(inBuffer,"An axis fault occurred") != NULL) {
      axisFaultCount++;
    }
    

    if (strstr(inBuffer,"MarkerInput") != NULL) { 
      if (strstr(inBuffer,"MarkerInput x: 0 y: 0 z: 0 u: 0") == NULL) {
        printf("ERROR: '%s' in '%s'\n",inBuffer,tailName);
        markerFaultCount++;
      }
    }
    charPtr = strstr(inBuffer,"Barcode Scanned:");
    if (charPtr != NULL) {
      if (strstr(inBuffer,"LEFT") != NULL) {
        lastScanPattern = SCAN_PATTERN_LEFT;
      } else if (strstr(inBuffer,"RIGHT") != NULL) {
        lastScanPattern = SCAN_PATTERN_RIGHT;
      } else {
        charPtr += strlen("Barcode Scanned:")+1;
        local_strlwr(charPtr);
        if (ParseFilename2(charPtr,barcodeSeries,&barcodePlateNumber) != 0) {
          switch(lastScanPattern) {
          case SCAN_PATTERN_NONE:
            strcpy(noneBarcodeSeries,barcodeSeries);
            noneBarcodePlateNumber = barcodePlateNumber;
            break;
          case SCAN_PATTERN_LEFT:
            strcpy(leftBarcodeSeries,barcodeSeries);
            leftBarcodePlateNumber = barcodePlateNumber;
            break;
          case SCAN_PATTERN_RIGHT:
            strcpy(rightBarcodeSeries,barcodeSeries);
            rightBarcodePlateNumber = barcodePlateNumber;
            break;
          default:
            printf("ERROR: illegal pattern in line %d lineCount %d for %s\n",__LINE__,lineCount,tailName);
            break;
          }
        }
      }
    }


    if (strstr(inBuffer,"Unrecognized") != NULL) {
      printf("%s ",copyLine);
      unrecognizedCount++;
      if (lastScanPattern != SCAN_PATTERN_NONE) {
        if (strstr(inBuffer,"Left") != NULL) {
          printf("%05s%05d\n",leftBarcodeSeries,leftBarcodePlateNumber);
        } else if (strstr(inBuffer,"Right") != NULL) {
          printf("%05s%05d\n",rightBarcodeSeries,rightBarcodePlateNumber);
        } else {
          printf("ERROR: unknown plate in line %d lineCount %d for %s\n",__LINE__,lineCount,tailName);
        }
      } else {
        charPtr = strrchr(inBuffer,':');
        if (charPtr == NULL) {
          printf("unsupported\n");
        } else {
          charPtr++;
          charPtr++;
          local_strlwr(charPtr);
          if (ParseFilename2(charPtr,barcodeSeries,&barcodePlateNumber) != 0) {
            printf("%05s%05d\n",barcodeSeries,barcodePlateNumber);
          } else {
            printf("unsupported\n");          
          }
        }
      }
    }
    if (strstr(inBuffer,"timeout") != NULL) {
      timeoutCount++;
    }
    if (strstr(inBuffer,"File Name") != NULL) {
#if 0
      printf("File: lineCount %d fileNameCount %d last0MedianFileLineNumber %d last0MedianFileNameCount %d last1MedianFileLineNumber %d last1MedianFileNameCount %d last2MedianFileLineNumber %d last2MedianFileNameCount %d \n",
             lineCount,
             fileNameCount,
             last0MedianFileNameCount,
             last0MedianFileLineNumber,
             last1MedianFileNameCount,
             last1MedianFileLineNumber,
             last2MedianFileNameCount,
             last2MedianFileLineNumber);
#endif
      pZeroMedianLast2 = pZeroMedianLast1;
      pZeroMedianLast1 = pZeroMedianCurrent;
      pZeroMedianCurrent = &zeroMedianTable[fileNameCount];
      if (strlen(copyLine) >= (MAX_BUFFER-5)) {
        printf("ERROR: MAX_BUFFER exceeded in line %d lineCount %d %s\n",__LINE__,lineCount,copyLine);
      }
      strcpy(pZeroMedianCurrent->fileBuffer,copyLine);
      pZeroMedianCurrent->fileNameCount = fileNameCount;
      pZeroMedianCurrent->lineCount = lineCount;
      pZeroMedianCurrent->preceedingZeroMedianFileNameCount = -1;
      pZeroMedianCurrent->preceedingZeroMedianLineNumber = -1;
      pZeroMedianCurrent->followingZeroMedianFileNameCount = last1MedianFileNameCount;
      pZeroMedianCurrent->followingZeroMedianLineNumber = last1MedianFileLineNumber;
      if ((pZeroMedianLast1 != NULL) && (last0MedianFileNameCount >= 0)) {
        pZeroMedianLast1->preceedingZeroMedianFileNameCount = last0MedianFileNameCount; /* Save the first before the actual overexposure */
        pZeroMedianLast1->preceedingZeroMedianLineNumber = last0MedianFileLineNumber; /* Save the first before the actual overexposure */
      }
      if ((pZeroMedianLast2 != NULL) && (last0MedianFileNameCount >= 0)) {
        pZeroMedianLast2->preceedingZeroMedianFileNameCount = last0MedianFileNameCount; /* Save the second before the actual overexposure */
        pZeroMedianLast2->preceedingZeroMedianLineNumber = last0MedianFileLineNumber; /* Save the second before the actual overexposure */
      }
      if ((pZeroMedianLast2 != NULL) && (last2MedianFileNameCount >= 0)) {
        pZeroMedianLast2->followingZeroMedianFileNameCount = last2MedianFileNameCount; /* Save the first before the actual overexposure */
        pZeroMedianLast2->followingZeroMedianLineNumber = last2MedianFileLineNumber; /* Save the first before the actual overexposure */
      }

#if 0
      printf("LINE %d last2MedianFileLineNumber %d last1MedianFileLineNumber %d last0MedianFileLineNumber %d\n",__LINE__,last2MedianFileLineNumber,last1MedianFileLineNumber,last0MedianFileLineNumber);
      if (pZeroMedianLast2 != NULL) {
        printf("fileNameCount %5d %5d lineCount %6d preceeding count %5d lineNumber %6d following count %5d lineNumber %6d for '%s'\n",
               fileNameCount-2,
               pZeroMedianLast2->fileNameCount,
               pZeroMedianLast2->lineCount,
               pZeroMedianLast2->preceedingZeroMedianFileNameCount,
               pZeroMedianLast2->preceedingZeroMedianLineNumber,
               pZeroMedianLast2->followingZeroMedianFileNameCount,
               pZeroMedianLast2->followingZeroMedianLineNumber,
               pZeroMedianLast2->fileBuffer);
      }
      if (pZeroMedianLast1 != NULL) {
        printf("fileNameCount %5d %5d lineCount %6d preceeding count %5d lineNumber %6d following count %5d lineNumber %6d for '%s'\n",
               fileNameCount-1,
               pZeroMedianLast1->fileNameCount,
               pZeroMedianLast1->lineCount,
               pZeroMedianLast1->preceedingZeroMedianFileNameCount,
               pZeroMedianLast1->preceedingZeroMedianLineNumber,
               pZeroMedianLast1->followingZeroMedianFileNameCount,
               pZeroMedianLast1->followingZeroMedianLineNumber,
               pZeroMedianLast1->fileBuffer);
      }
      printf("fileNameCount %5d %5d lineCount %6d preceeding count %5d lineNumber %6d following count %5d lineNumber %6d for '%s'\n",
             fileNameCount,
             pZeroMedianCurrent->fileNameCount,
             pZeroMedianCurrent->lineCount,
             pZeroMedianCurrent->preceedingZeroMedianFileNameCount,
             pZeroMedianCurrent->preceedingZeroMedianLineNumber,
             pZeroMedianCurrent->followingZeroMedianFileNameCount,
             pZeroMedianCurrent->followingZeroMedianLineNumber,
             pZeroMedianCurrent->fileBuffer);
#endif
      last2MedianFileNameCount = last1MedianFileNameCount;
      last2MedianFileLineNumber = last1MedianFileLineNumber;
      last1MedianFileNameCount = last0MedianFileNameCount;
      last1MedianFileLineNumber = last0MedianFileLineNumber;
      last0MedianFileNameCount = -1;
      last0MedianFileLineNumber = -1;

      fileNameCount++;
    }
    
    if (strstr(inBuffer," 0 MEDIAN") != NULL) {
#if 0
      printf("0 MEDIAN found in line %d at fileNameCount %d\n",lineCount,fileNameCount);
#endif
      last0MedianFileNameCount = fileNameCount;
      last0MedianFileLineNumber = lineCount;
        
      zeroMedianCount++;
    }
  }
#if 0
  printf("fileNameCount %d fileNameTotal %d\n",fileNameCount,fileNameTotal);
#endif
  fileNameTotal = fileNameCount;
  for (fileNameCount = 0; fileNameCount < fileNameTotal; fileNameCount++) {
    pZeroMedian = &zeroMedianTable[fileNameCount];
#if 0
    printf("fileNameCount %5d %5d lineCount %6d preceeding count %5d lineNumber %6d following count %5d lineNumber %6d for '%s'\n",
           fileNameCount,
           pZeroMedian->fileNameCount,
           pZeroMedian->lineCount,
           pZeroMedian->preceedingZeroMedianFileNameCount,
           pZeroMedian->preceedingZeroMedianLineNumber,
           pZeroMedian->followingZeroMedianFileNameCount,
           pZeroMedian->followingZeroMedianLineNumber,
           pZeroMedian->fileBuffer);
#endif
    if (strstr(pZeroMedian->fileBuffer,"Plates") != NULL) {
      if ((pZeroMedian->preceedingZeroMedianFileNameCount >= 0) || (pZeroMedian->followingZeroMedianFileNameCount >= 0)) {
        fprintf(pParseCommon->out_handle,"%s\n",pZeroMedian->fileBuffer);
      }
    } else {
      if (pZeroMedian->followingZeroMedianFileNameCount >= 0) {
        fprintf(pParseCommon->out_handle,"%s\n",pZeroMedian->fileBuffer);
      }
    }
  }


  printf("File: %s wrongPositionCount %d x %d y %d z %d u %d axisFaultCount %d markerFaultCount %d unrecognizedCount %d timeoutCount %d zeroMedianCount %d\n",
         tailName,
         wrongPositionCount,
         xDeltaCount,
         yDeltaCount,
         zDeltaCount,
         uDeltaCount,
         axisFaultCount,
         markerFaultCount,
         unrecognizedCount,
         timeoutCount,
         zeroMedianCount);
  pParseCommon->wrongPositionCount  += wrongPositionCount;
  pParseCommon->xDeltaCount += xDeltaCount;
  pParseCommon->yDeltaCount += yDeltaCount;
  pParseCommon->zDeltaCount += zDeltaCount;
  pParseCommon->uDeltaCount += uDeltaCount;
  pParseCommon->axisFaultCount  += axisFaultCount;
  pParseCommon->markerFaultCount  += markerFaultCount;
  pParseCommon->unrecognizedCount  += unrecognizedCount;
  pParseCommon->timeoutCount  += timeoutCount;
  pParseCommon->zeroMedianCount  += zeroMedianCount;

  if (pParseCommon->verbose != 0) {
    printf("Finished processing %s lines %d\n",log_name,lineCount);
  }
  fclose(log_handle);
  if (zeroMedianTable != NULL) {
    free(zeroMedianTable);
  }
  return(0);
}

void ProcessDirectory(PPARSECOMMON pParseCommon,char *rootDirectory)
{
  struct dirent **filelist;
  int numFiles;
  int fileIndex;
  char filename[MAX_BUFFER];
  struct stat statbuf;
  char cmdBuf[MAX_BUFFER];
  char cmdBuf2[MAX_BUFFER];
  char targetFile[MAX_BUFFER];
  char targetDirectory[MAX_BUFFER];
  int result;
  int plateFlag;
  off_t initialSize;
  off_t finalSize;
  char* charPtr;
  int targetLength;
  int firstFileFlag = 1;
  if (pParseCommon->verbose != 0) {
    printf("Processing Directory %s \n",rootDirectory);
  }

  if (pParseCommon->abortFlag != 0) {
    return;
  }
  if (strlen(rootDirectory) > MAX_BUFFER-10) {
    printf("ERROR: rootDirectory %s too long in line %d\n",rootDirectory,__LINE__);
    return;
  }

  pParseCommon->directoryCount++;
  numFiles = scandir(rootDirectory,&filelist,0,alphasort);
  if (numFiles > 0) {
    for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
      if (filelist[fileIndex]->d_name[0] == '.') {
        /* skip "." and ".." */
        continue;
      }
      /* If it has a dot, assume it is not a directory */
      pParseCommon->allFileCount++;
      strcpy(filename,rootDirectory);
      strcat(filename,"/");
      strcat(filename,filelist[fileIndex]->d_name);
      if ((strlen(rootDirectory)+strlen(filelist[fileIndex]->d_name)) > MAX_BUFFER-10) {
        printf("ERROR: filename %s/%s too long in line %d\n",rootDirectory,filelist[fileIndex]->d_name,__LINE__);
        continue;
      }
      if (strstr(filelist[fileIndex]->d_name,".") != NULL) {
        if (strstr(filename,".log") == NULL) {
          continue;
        }
        charPtr = strstr(filename,"PlateImage");
        if (charPtr == NULL) {
          continue;
        }
#ifdef POST_FLOOD
        if (strcmp(charPtr,POST_FLOOD) < 0) {
          continue;
        }
#endif /* POST_FLOOD */
        ProcessLogBuffer(pParseCommon,filename);
      }
#if 0
      if (strcmp(filelist[fileIndex]->d_name,"Logs") == 0) {
        printf("At line %d filename: %s\n",__LINE__,filename);
      }
#endif


      result = stat(filename,&statbuf);
      if (result != 0) {
        printf("ERROR: stat failed on %s\n",filename);
        pParseCommon->abortFlag = 1;
        return;
      }


      if  (S_ISLNK(statbuf.st_mode) != 0) {
        /* Skip links */
#if 0
        printf("Skip link at line %d\n",__LINE__);
#endif
        continue;
      }      

      /* This is a directory.  Recurse into it */
      ProcessDirectory(pParseCommon,filename);
          
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
  char out_name[MAX_INPUT_NAME];
  char *argstr;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  char list_name[MAX_INPUT_NAME];
  FILE *list_handle = NULL;
  int errorFlag = 0;
  char cmdchar;
  int lineLen;
  time_t curTime;


  PARSECOMMON parsecommon;
  PPARSECOMMON pParseCommon = &parsecommon;


  int statResult;
  struct stat filestats;
 
  list_name[0] = 0;
  out_name[0] = 0;
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
    statResult = stat(list_name,&filestats);
    if (statResult != 0) {
      errorFlag = 1;
      printf("Could not find list file %s\n",list_name);
    }
  }
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


  if (errorFlag) {
    printf("Usage: parseplateimage [-v]   -o <outfile>\n");
    printf("       where -l is <list of log files ending in .list>\n");
    printf("                or <single file with name PlateImage*.log>\n");
    printf("                or <directory tree to search\n");
    printf("       where -v is the verbose flag\n");
    printf("             -o <outfile> contains candiates with '0 MEDIAN'\n");

    return(-1);
  }

  printf("parseplateimage of %s %s List Filename %s Output Filename %s\n",
         __DATE__,__TIME__,list_name,out_name);
 

  time(&pParseCommon->startTime);
  if (strstr(list_name,".list") != NULL) {
    /* This must be a list of files */
    list_handle = fopen(list_name,"rt");
    if (list_handle == NULL) {
      printf("ERROR: Could not open list file %s\n",list_name);
    }
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
  } else if ((strstr(list_name,"PlateImage") != NULL) && (strstr(list_name,".log") != NULL)) {
    /* This is a single file */
    ProcessLogBuffer(pParseCommon,list_name);       
  } else {
    /* This is must be a directory.  Recurse into it */
    ProcessDirectory(pParseCommon,list_name);
  }


  if (list_handle != NULL) {
    fclose(list_handle);
  }
  if (pParseCommon->out_handle != NULL) {
    fclose(pParseCommon->out_handle);
  }
  time(&curTime);
  curTime -= pParseCommon->startTime;

  if (pParseCommon->fileCount > 1) {
    printf("Total: wrongPositionCount %d x %d y %d z %d u %d  axisFaultCount %d markerFaultCount %d unrecognizedCount %d timeoutCount %d zeroMedianCount %d\n",
           pParseCommon->wrongPositionCount,
           pParseCommon->xDeltaCount,
           pParseCommon->yDeltaCount,
           pParseCommon->zDeltaCount,
           pParseCommon->uDeltaCount,
           pParseCommon->axisFaultCount,
           pParseCommon->markerFaultCount,
           pParseCommon->unrecognizedCount,
           pParseCommon->timeoutCount,
           pParseCommon->zeroMedianCount);
  }

  printf("Execution Time: %lld seconds. maxLineLen %d\n",
         curTime,
         pParseCommon->maxLineLen);



  return(EXIT_SUCCESS);
}
