// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* checkdeletedplates.c
 *
 *  This program searches /dasch/Pipeline/ingest for stale SExtactor files not in /dasch/Pipeline/total.list
 *
 *
 *
 * cc -ggdb -O2 -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I /dasch/install/include -L /dasch/install/lib -lm checkdeletedplates.c pipelineutils.a -ltable -lutil -lwcs -L/usr/lib${lib64}/mysql -l mysqlclient -o checkdeletedplates
 *
 *
 * On April 28, 2014, 12199 files from 11498 scans were removed, freeing 1625205 MB, or 141 MB/scan.
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
#include <dirent.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER     1024
#define MAX_TIME_LENGTH 28
#define MAX_SOLUTION_NUMBER 10
#define MAX_PLATES 5000000
#define MAX_PATH 512

typedef struct _seriesEntry {
  char series[MAX_SERIES_STRING];
  int minIndex;
  int maxIndex;
}  SERIESENTRYX,*PSERIESENTRYX;

typedef struct _plateEntry {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int foundGood;
  int foundBad;
} PLATEENTRY,*PPLATEENTRY;

int PlateCompare(const void *first, const void *second) 
{
  
  PPLATEENTRY pFirst = (PPLATEENTRY)first;
  PPLATEENTRY pSecond = (PPLATEENTRY)second;
  int result;

  result = strcmp(pFirst->series,pSecond->series);
  if (result != 0) {
    return(result);
  }
  if (pFirst->plateNumber > pSecond->plateNumber) {
    return(1);
  } else if (pFirst->plateNumber < pSecond->plateNumber) {
    return(-1);
  } else {
    printf("ERROR: two plates for %s%05d\n",pFirst->series,pFirst->plateNumber);
    exit(-1);
  }

}


int main(int argc,char *argv[])
{
	char *argstr;
	char *inBuffer;
	char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
	char total_name[] = "/dasch/Pipeline/total.list";
  char rootDirectory[] = "/dasch/Pipeline/match";
  char filename[MAX_PATH];
	FILE *total_handle = NULL;
	int errorFlag = 0;
	char cmdchar;
	int lineLen;
  int lineNumber = 0;
	time_t curTime;
  time_t startTime;
  int verbose = 0;
  int goodLines = 0;
  time_t totalTime = 0;
	time_t totalAverageTime = 0;
  int maxLineLen = 0;
  PPLATEENTRY plate_table = NULL;
  PPLATEENTRY pPlate;
  int total_plates = 0;
  int plateIndex;
  SERIESENTRYX series_table[MAX_SERIES];
  PSERIESENTRYX pSeries;
  int seriesCount = 0;
  int seriesIndex = 0;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  struct dirent **filelist;
  int numFiles;
  int fileIndex;
  char *charPtr;
  int goodPlateCount = 0;
  int badPlateCount = 0;
  int result;
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
					verbose = 1;
					break;



				default:
					printf("* illegal command -%c-",cmdchar);
					errorFlag = 1;
					break;
				}
        
			}
		}
	}



  total_handle = fopen(total_name,"rt");
  if (total_handle == NULL) {
    errorFlag = 1;
    printf("Could not open list file %s\n",total_name);
  }




	if (errorFlag) {
		printf("Usage: checkdeletedplates [-v] \n");
		printf("       where -v is the verbose flag\n");
		return(-1);
	}

	printf("checkdeletedplates of %s %s Total Filename %s\n",
				 __DATE__,__TIME__,total_name);
 

	time(&startTime);
  plate_table = (PPLATEENTRY)calloc(MAX_PLATES,sizeof(PLATEENTRY));



	while (1) {
		inBuffer = fgets(inLine,MAX_BUFFER,total_handle);
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
    pPlate = &plate_table[total_plates];
    memset(pPlate,0,sizeof(PLATEENTRY));
    if (ParseFilename(inBuffer,pPlate->series,&pPlate->plateNumber,&pPlate->mosaicNumber,&pPlate->binning,&pPlate->rotation) == 0) {
      printf("ERROR: failed to parse %s in total.list\n",inBuffer);
    } else {
      if (total_plates < MAX_PLATES) {
        total_plates++;
        goodLines++;
      } else {
        printf("ERROR: total_plates %d exceeds MAX_PLATES\n",total_plates);
      }
        

    }

	}

	if (total_handle != NULL) {
		fclose(total_handle);
	}
  qsort((void*)plate_table,total_plates,sizeof(PLATEENTRY),PlateCompare);
  series[0] = 0;
  memset(series_table,0,sizeof(series_table));
  pSeries = &series_table[0];
  memset(pSeries,0,sizeof(SERIESENTRYX));
  /* Now form an index into each series */
  for (plateIndex = 0; plateIndex < total_plates; plateIndex++) { 
    pPlate = &plate_table[plateIndex];
    if (strcmp(pPlate->series,pSeries->series) != 0) {
      if (seriesCount == 0) {
        strcpy(pSeries->series,pPlate->series);
        pSeries->minIndex = plateIndex;
        seriesCount++;        
      } else {
        pSeries->maxIndex = plateIndex;
        pSeries = &series_table[seriesCount];
        memset(pSeries,0,sizeof(SERIESENTRYX));
        strcpy(pSeries->series,pPlate->series);
        pSeries->minIndex = plateIndex;
        seriesCount++;
      }
    }
#if 0
    printf("Entry %6d %5s%05d\n",plateIndex,pPlate->series,pPlate->plateNumber);
#endif
  }
  pSeries->maxIndex = plateIndex;
  

	time(&curTime);
	curTime -= startTime;
  numFiles = scandir(rootDirectory,&filelist,0,alphasort);
  printf("Found %d files in %s\n",numFiles,rootDirectory);
  if (numFiles > 0) {
    for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
      if (filelist[fileIndex]->d_name[0] == '.') {
        /* skip "." and ".." */
        continue;
      }
      strcpy(filename,rootDirectory);
      strcat(filename,"/");
      strcat(filename,filelist[fileIndex]->d_name);
      charPtr = strstr(filelist[fileIndex]->d_name,"match_");
      if (charPtr == NULL) {
        charPtr = filename;
      } else if (charPtr != filelist[fileIndex]->d_name) {
        charPtr = filename;
      } else {
        charPtr = charPtr + 6;
      }
      if (ParseFilename(charPtr,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
        printf("ERROR: failed to parse %s\n",filename);
      } else {
        /* Look up the plate in the table */
        for (seriesIndex = 0; seriesIndex < seriesCount; seriesIndex++) {
          pSeries = &series_table[seriesIndex];
          if (strcmp(series,pSeries->series) == 0) {
            break;
          }
        }
        if (seriesIndex >= seriesCount) {
          printf("ERROR: could not find series in %s\n",filename);
        } else {
          for (plateIndex = pSeries->minIndex; plateIndex < pSeries->maxIndex; plateIndex++) { 
            pPlate = &plate_table[plateIndex];
            if (pPlate->plateNumber == plateNumber) {
              if (pPlate->mosaicNumber ==mosaicNumber) {
                /* Got a good entry */
                pPlate->foundGood++;
                goodPlateCount++;
              } else {
                pPlate->foundBad--;
                printf("ERROR entry %s does not agree with plate %s%05d_%02d\n",filename,pPlate->series,pPlate->plateNumber,pPlate->mosaicNumber);
                badPlateCount++;
              }
                
            }
          }
        }
      }
   
    }
    for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
      free(filelist[fileIndex]);
    }
    free(filelist);
  
  }


	printf("Execution Time: %lld seconds. lines %d maxLineLen %d goodLines %d goodPlateCount %d badPlateCount %d totalTime %lld totalAverageTime %lld\n",
				 curTime,
         lineNumber,
				 maxLineLen,
         goodLines,
         goodPlateCount,
         badPlateCount,
				 totalTime,
				 totalAverageTime);





	return(EXIT_SUCCESS);
}
