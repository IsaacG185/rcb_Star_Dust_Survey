// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* sortplatelist.c
 *
 *  Takes a list of plates and sorts them into long focus, patrol, and meteor scales 
 *
 * 
gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I/usr/include/mysql  -I /dasch/install/include  -L /dasch/install/lib -lm sortplatelist.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread plottransientstub.o -o sortplatelist

sortplatelist -i /home/scanner/junk/gaia/gaia_dasch8.list -o /home/scanner/junk/gaia/los8.tmp

creates four files: los0.tmp through los4.tmp
sortplatelist -n 4 -i /home/scanner/junk/atlas/atlas_pending.list -o /home/scanner/junk/atlas/losX.tmp

sortplatelist -s -i /home/scanner/backup/2018_07_12/getseries.log -o /home/scanner/junk/los2.log

The following looks for saturated plates in the scanner log files:

cd ~/backup/2018_11_15/
find . -name "PlateImage*.log" -exec egrep -i "File Name| 0 MEDIAN" {} \; > /home/scanner/junk/plates/los.tmp
cd /home/scanner/junk/plates/
grep -A 1 -B 1 " 0 MEDIAN" los.tmp | sort -u > los1.tmp


sortplatelist -m -i /home/scanner/junk/plates/los.bat -o /home/scanner/junk/plates/los2.bat -i /home/scanner/junk/plates/los1.tmp


 *
 * Jul 11, 2018 Edward J. Los - Initial version  
 * Jul 13, 2018 Edward J. Los - Add '-s' qualifier to parse only the series
 * Jul 27, 2018 Edward J. Los - sort by series instead of id
 *                              Add '-p' for a plate type prefix
 * Nov 20, 2018 Edward J. Los - Add '-m' to handle a "MosaicCheck" list
 * Nov 22, 2018 Edward J. Los - Change /dasch/raid020 -> /n/dasch14 to /n/dasch15
 *                            - Add '-e" for zero median error list to look for saturated plates
 *                            - Reverse sort '-m' plate numbers to get the last plate in each series
 * Dec 11, 2018 Edward J. Los - correct initial value of mosaicCheckFlag
 * Mar  6, 2019 Edward J. Los - use -n <numfiles> to create a numfiles list with round-robin placement
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
#define MAX_LINE_LENGTH 110
#define MAX_BUFFER     1024
#define PLATE_TYPE_LONGFOCUS 1
#define PLATE_TYPE_PATROL    2
#define PLATE_TYPE_METEOR    3
#define PLATE_TYPE_MAX       4
#define MASK_SERIES_ID  (MAX_SERIES+1)
#define MIN_PLATES 2
#define MAX_PLATES 3
#define PLATES_CUTOFF 50 /* above this count, output MAX_PLATES */
#define MAX_FILES 8

char *plateTypeString[PLATE_TYPE_MAX] = {
  "UNKNOWN  ",
  "LONGFOCUS",
  "PATROL   ",
  "METEOR   "
  };

typedef struct _platestats {
  char fullname[MAX_LINE_LENGTH];
	char series[MAX_SERIES_STRING];
	int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int plateType;
  int seriesId;
  int medianErrorFlag; /* Set if the PlateImage<date>.log file has a "MEDIAN 0" error near when this file was created */
  double fittedPlateScale; /* This is degrees per pixel */
} PLATESTATS,*PPLATESTATS;

typedef struct _seriesstats {
  char series[MAX_SERIES_STRING];
  int plateCount;
  int outputCount1; /* all plates output */
  int outputCount2; /* non-special plates output */
  int allOutputFlag;
} SERIESSTATS,*PSERIESSTATS;



int PlateCompare(const void *first, const void *second) 
{
  
  PPLATESTATS pFirst = (PPLATESTATS)first;
  PPLATESTATS pSecond = (PPLATESTATS)second;
  int result;

  if (pFirst->plateType < pSecond->plateType) {
    return(-1);
  } else if (pFirst->plateType > pSecond->plateType) {
    return(1);
  } else {
    result = strcmp(pFirst->series,pSecond->series);
    if (result < 1) {
      return(-1);
    } else if (result > 1) {
      return(1);
    } else {
      if (pFirst->plateNumber < pSecond->plateNumber) {
        return(-1);
      } else if (pFirst->plateNumber > pSecond->plateNumber) {
        return(1);
      } else {
        return(0);
      }
    }
  }
}
int ReversePlateCompare(const void *first, const void *second) 
{
  
  PPLATESTATS pFirst = (PPLATESTATS)first;
  PPLATESTATS pSecond = (PPLATESTATS)second;
  int result;

  if (pFirst->plateNumber < pSecond->plateNumber) {
    return(1);
  } else if (pFirst->plateNumber > pSecond->plateNumber) {
    return(-1);
  } else {
    return(0);
  }
}

int main(int argc,char *argv[])
{
	char *argstr;
	char *inBuffer;
	char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
	char input_name[MAX_INPUT_NAME];
	FILE *input_handle = NULL;
  char out_name_template[MAX_INPUT_NAME];
	char out_name[MAX_FILES][MAX_INPUT_NAME];
  FILE *out_handle[MAX_FILES] ;
	char median_name[MAX_INPUT_NAME];
  FILE *median_handle = NULL;


	int errorFlag = 0;
	char cmdchar;
	int lineLen;
  time_t startTime;
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
  char *mysqlhost;
  char *username;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  char *password;
  int maxPlateCount = 0;
  int curPlateCount = 0;
  int plateIndex;
  int maxPlateString = 0;
  int lineLength;
  int longFocusPlateCount = 0;
  int patrolPlateCount = 0;
  int meteorPlateCount = 0;
  int lineCounter = 0;

  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  PPLATESTATS plateTable = NULL;
  PPLATESTATS pPlate;
  SERIESSTATS seriesStatsTable[MAX_SERIES+2];
  PSERIESSTATS pSeriesStats;
  int seriesIndex;
  int totalOutputCount = 0;
  int parseSeriesOnlyFlag = 0;
  int prefixFlag = 0;
  int mosaicCheckFlag = 0;
  int mosaicPrintFlag = 0;
  int mosaicSpecialFlag = 0;
  char *charPtr;
  char *charPtr2;
  char *charPtr3;
  char *charPtr4;
  char charVal;
  int charCount;
  int fileIndex = 0;
  int fileModulus = -1;
  int numFiles = 1;
  int nvals;


	/* Loop through the arguments */
	input_name[0] = 0;
	out_name_template[0] = 0;
  memset(out_name,0,sizeof(out_name));
  memset(out_handle,0,sizeof(out_handle));
  median_name[0] = 0;
  memset(seriesStatsTable,0,sizeof(seriesStatsTable));
  pSeriesStats = &seriesStatsTable[MASK_SERIES_ID];
  pSeriesStats->allOutputFlag = 1;

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

				case 's':
				case 'S':
					parseSeriesOnlyFlag = 1;
					break;

				case 'p':
				case 'P':
					prefixFlag = 1;
					break;

				case 'm':
				case 'M':
					mosaicCheckFlag = 1;
					break;

				case 'i': /* list file name */
				case 'I':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(input_name,*++argv,MAX_INPUT_NAME-2);
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
						strncpy(out_name_template,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;

				case 'n': /* number of files */
				case 'N':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
            nvals = sscanf(*++argv,"%d",&numFiles);
            if (nvals != 1) {
              printf("ERROR: failed to decode numFiles\n");
              errorFlag =1;
            } else {
              if ((numFiles < 1) || (numFiles > MAX_FILES)) {
                printf("ERROR: illegal numFiles. Must be 1 to %d\n",MAX_FILES);
                errorFlag = 1;
              }
            }
					}
					break;

				case 'e': /* median error file name */
				case 'E':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(median_name,*++argv,MAX_INPUT_NAME-2);
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



	if (input_name[0] == 0) {
		printf("ERROR: List file name not specified\n");
		errorFlag = 1;
	} else {
		input_handle = fopen(input_name,"rt");
		if (input_handle == NULL) {
			errorFlag = 1;
			printf("Could not open list file %s\n",input_name);
		}
	}
  



	if (out_name_template[0] == 0) {
		printf("ERROR: Out file name not specified\n");
		errorFlag = 1;
	} else {
    if (numFiles == 1) {
      strcpy(out_name[0],out_name_template);
      out_handle[0] = fopen(out_name[0],"wt");
      if (out_handle[0] == NULL) {
        errorFlag = 1;
        printf("Could not open out file %s\n",out_name[0]);
      }
    } else {
      if (strchr(out_name_template,'X') == NULL) {
        printf("ERROR: %s must contain an 'X'\n",out_name_template);
        errorFlag =1;
      }
      for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
        strcpy(out_name[fileIndex],out_name_template);
        charPtr4 = strchr(out_name[fileIndex],'X');
        if (charPtr4 == NULL) {
          printf("ERROR: bad charPtr4\n");
          exit(-1);
        }
        *charPtr4 = '0' + fileIndex;
        out_handle[fileIndex] = fopen(out_name[fileIndex],"wt");
        if (out_handle[fileIndex] == NULL) {
          errorFlag = 1;
          printf("Could not open out file %s\n",out_name[fileIndex]);
        } else {
          printf("Opened output file %d %s\n",fileIndex,out_name[fileIndex]);
        }
      }
    }
  }

  if (mosaicCheckFlag == 1) {
    if (median_name[0] == 0) {
      printf("ERROR: Median file name not specified\n");
      errorFlag = 1;
    } else {
      median_handle = fopen(median_name,"rt");
      if (median_handle == NULL) {
        errorFlag = 1;
        printf("Could not open out file %s\n",median_name);
      }
    }
  }


  if (errorFlag) {
    printf("Usage: sortplatelist  -i <input list> -o <outfile>\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s parses the series only\n");
    printf("       where -m handles the production mosaic check list\n");
    printf("       where -e is the 'MEDIAN 0' error file name\n");
    printf("       where -n <numfiles> is the number of files to split the list into\n");
    printf("                the output file becomes a template which must have one 'X'\n");
    return(-1);
  }

  printf("sortplatelist of %s %s List Filename %s Output Filename %s 'MEDIAN 0' Filename %s numFiles %d\n",
         __DATE__,__TIME__,input_name,out_name,median_name,numFiles);
 

  time(&startTime);


  /* Connect to the database */
  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    printf("ERROR: DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    printf("ERROR: DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    printf("ERROR: DASCH_PASSWORD is not defined\n");
    return(-1);
  }
                   
  mysql_init(pConnection);

  if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
    }
    return(-1);
  }


  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    printf("DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }

  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    printf("DASCH_PHOT_USERNAME is not defined\n");
    return(-1);
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    printf("DASCH_PHOT_PASSWORD is not defined\n");
    return(-1);
  }
  mysql_init(pPhotConnection);


  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
    }
    return(-1);
  }



  InitSeriesTable(pConnection,pPhotConnection);

  mysql_close(pConnection);
  mysql_close(pPhotConnection);


  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
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
      maxPlateCount++;
      lineLength = strlen(inBuffer);
      if (lineLength > maxPlateString) {
        maxPlateString = lineLength;
      }
    }
  }
  if (input_handle != NULL) {
    fclose(input_handle);
  }
  printf("Found %d plates with %d maxPlateString in %s\n",maxPlateCount,maxPlateString,input_name);
  
  plateTable = (PPLATESTATS)calloc(maxPlateCount,sizeof(PLATESTATS));
  if (plateTable == NULL) {
    printf("ERROR: failed to allocate plateTable\n");
    exit(-1);
  }

  input_handle = fopen(input_name,"rt");
  if (input_handle == NULL) {
    errorFlag = 1;
    printf("Could not open list file (2) %s\n",input_name);
  }
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
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
      pPlate = &plateTable[curPlateCount];
      memset(pPlate,0,sizeof(PLATESTATS));
      lineLength = strlen(inBuffer);
      if (lineLength > (MAX_LINE_LENGTH-2)) {
        printf("ERROR: lineLength %d exceeded for %s\n",lineLength,inBuffer);
        continue;
      }
      strcpy(pPlate->fullname,inBuffer);
      if (parseSeriesOnlyFlag == 0) {

        if (ParseFilename(pPlate->fullname,pPlate->series,&pPlate->plateNumber,&pPlate->mosaicNumber,&pPlate->binning,&pPlate->rotation) == 0) {
          printf("ERROR: failed to parse %s\n",pPlate->fullname);
          continue;
        }
      } else { /* parseSeriesOnlyFlag != 0 */
        /* Find the series only */
        charVal = -1;
        charPtr = inBuffer;
        charPtr2 = pPlate->series;
        charCount = 0;
        while (charVal != 0) {
          charVal = *charPtr;
          charPtr++;
          if (charVal == ' ') {
            continue;
          }
          if ((charVal >= 'a') && (charVal <= 'z')) {
            *charPtr2 = charVal;
            *charPtr2++;
            charCount++;
            if (charCount >= MAX_SERIES_STRING) {
              printf("ERROR: series length exceeded in %s\n",inBuffer);
              break;
            } else {
              continue;
            }
          }
          *charPtr2 = 0;
          break;
        }
        if (charCount >= MAX_SERIES_STRING) {
          continue;
        }

      }
      pPlate->seriesId = GetSeriesId(pPlate->series,0);
      if (pPlate->seriesId < 0) {
        if ((mosaicCheckFlag == 1) && (strcmp(pPlate->series,"mask") == 0)) {
          pPlate->seriesId = MASK_SERIES_ID;
        } else {
          printf("ERROR: failed to get seriesId for %s\n",pPlate->fullname);
          continue;
        }
      }
      if (mosaicCheckFlag == 1) {
        pSeriesStats = &seriesStatsTable[pPlate->seriesId];
        if (pSeriesStats->plateCount == 0) {
          strcpy(pSeriesStats->series,pPlate->series);
          if (strcmp(pSeriesStats->series,"a") == 0) {
            pSeriesStats->allOutputFlag = 1;
          }
        }
        pSeriesStats->plateCount++;

      }

#if 0
      if (pPlate->seriesId == 5) {
        printf("Got seriesId %d\n",pPlate->seriesId);
      }
#endif
      if (mosaicCheckFlag == 0) {
        pPlate->fittedPlateScale = GetFittedPlateScale(pPlate->seriesId,pPlate->plateNumber);
        if (pPlate->fittedPlateScale > (MIN_METEOR_ASTROMETRY_SCALE/3600.)) {
          pPlate->plateType = PLATE_TYPE_METEOR;
          meteorPlateCount++;
        } else if (pPlate->fittedPlateScale > (MIN_PATROL_ASTROMETRY_SCALE/3600.)) {
          pPlate->plateType = PLATE_TYPE_PATROL;
          patrolPlateCount++;
        } else {
          pPlate->plateType = PLATE_TYPE_LONGFOCUS;
          longFocusPlateCount++;
        }
      }
      curPlateCount++;

    }
  }
  printf("Parsed %d plates: %d long focus, %d patrol, %d meteor\n",curPlateCount,longFocusPlateCount,patrolPlateCount,meteorPlateCount);
  if (mosaicCheckFlag == 0) {
    qsort((void*)plateTable,curPlateCount,sizeof(PLATESTATS),PlateCompare);
  } else {
    qsort((void*)plateTable,curPlateCount,sizeof(PLATESTATS),ReversePlateCompare);
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,median_handle);
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
      lineCounter++;
      charPtr = strstr(inBuffer,"File Name");
      if (charPtr != inBuffer) {
        continue;
      }
      strcpy(copyLine,inLine);
      charPtr = strstr(inLine,"Exposure_");
      if (charPtr != NULL) {
        charPtr = strstr(inLine,"Plates");
        if (charPtr == NULL) {
          printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
          continue;
        }
        charPtr2 = strtok(charPtr,"\\");
        if (charPtr2 == NULL) {
          printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
          continue;
        }
        charPtr2 = strtok(NULL,"\\");
        if (charPtr2 == NULL) {
          printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
          continue;
        }
        charPtr3 = strtok(NULL,"\\");
        if (charPtr3 == NULL) {
          printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
          continue;
        }
#if 0
        printf("line %d found %s %s\n",__LINE__,charPtr2,charPtr3);
#endif
        for (plateIndex = 0; plateIndex < curPlateCount; plateIndex++) {
          pPlate = &plateTable[plateIndex];
          if ((strcmp(charPtr2,pPlate->series) == 0) &&
              (strstr(pPlate->fullname,charPtr3) != 0)) {
            if (pPlate->medianErrorFlag == 0) {
              pPlate->medianErrorFlag = 1;
              printf("ERROR: possible saturated plate: %s%05d_%02d\n",pPlate->series,pPlate->plateNumber,pPlate->mosaicNumber);
              break;
            }            
          }
        }
        
      } else {
        charPtr = strstr(inLine,"DarkFrames");
        if (charPtr == NULL) {
          charPtr = strstr(inLine,"FlatFrames");
        }
        if (charPtr != NULL) {
          charPtr2 = strtok(charPtr,"\\");
          if (charPtr2 == NULL) {
            printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
            continue;
          }
          charPtr2 = strtok(NULL,"\\");
          if (charPtr2 == NULL) {
            printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
            continue;
          }
          printf("ERROR saturated file: find . -name \"PlateImage*.log\" -exec egrep \"%s|MEDIAN\" {} \\; | grep -B 2 -A 2 \"%s\"\n",charPtr2,charPtr2);

        } else {
          charPtr = strstr(inLine,"FlatStats");
          if (charPtr != NULL) {
            charPtr2 = strtok(charPtr,"\\");
            if (charPtr2 == NULL) {
              printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
              continue;
            }
            charPtr2 = strtok(NULL,"\\");
            if (charPtr2 == NULL) {
              printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
              continue;
            }
            charPtr3 = strtok(NULL,"\\");
            if (charPtr3 == NULL) {
              printf("ERROR: line %d illegal plate reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
              continue;
            }
            printf("ERROR saturated tile: find . -name \"PlateImage*.log\" -exec egrep \"%s\\\\%s|MEDIAN\" {} \\; | grep -B 2 -A 2 \"%s\\\\%s\"\n",charPtr2,charPtr3,charPtr2,charPtr3);
          } else {
            printf("ERROR: line %d illegal exposure reference in %s line %d %s\n",__LINE__,median_name,lineCounter,copyLine);
            continue;
          }

        }

      }
    }
  }
  for (plateIndex = 0; plateIndex < curPlateCount; plateIndex++) {
    if (numFiles == 1) {
      fileModulus = 0;
    } else {
      fileModulus = (fileModulus+1) % numFiles;
    }
    pPlate = &plateTable[plateIndex];
    if (mosaicCheckFlag == 0) {


#if 1
      if ((pPlate->plateType >= PLATE_TYPE_MAX) ||
          (pPlate->plateType <= 0)) {
        printf("ERROR: illegal plate type for %s\n",pPlate->fullname);
        exit(-1);
      }
      if (prefixFlag != 0) {
        fprintf(out_handle[fileModulus],"%s %s\n",plateTypeString[pPlate->plateType],pPlate->fullname);
      } else {
        fprintf(out_handle[fileModulus],"%s\n",pPlate->fullname);
      }
#else 
      fprintf(out_handle[fileModulus],"plateType %d seriesId %2d plateNumber %2d '%s'\n",pPlate->plateType,pPlate->seriesId,pPlate->plateNumber,pPlate->fullname);
#endif
    } else {
      /* Mosaic check mode */
      pSeriesStats = &seriesStatsTable[pPlate->seriesId];
      mosaicPrintFlag = 0;
      mosaicSpecialFlag = 0;
      if ((pSeriesStats->allOutputFlag != 0) ||
          (pPlate->medianErrorFlag != 0)) {
        mosaicPrintFlag = 1;
        mosaicSpecialFlag = 1;
      } else {
        if (pSeriesStats->plateCount > PLATES_CUTOFF) {
          if (pSeriesStats->outputCount2 <  MAX_PLATES) {
            mosaicPrintFlag = 1;
          }
        } else {
          if (pSeriesStats->outputCount2 <  MIN_PLATES) {
            mosaicPrintFlag = 1;
          }
        }
      }
      if (mosaicPrintFlag == 1) {
        totalOutputCount++;
        if (mosaicSpecialFlag == 0) {
          pSeriesStats->outputCount2++;
        } 
        pSeriesStats->outputCount1++;
            
        charPtr3 = strstr(pPlate->fullname,"/dasch/raid028");
        if (charPtr3 != NULL) {
          charVal = charPtr3[14];
          strcpy(charPtr3,"    /n/dasch15");
          charPtr3[14] = charVal;
        }

        fprintf(out_handle[fileModulus],"%s\n",pPlate->fullname);
      }
    }

  }
  if (mosaicCheckFlag != 0) {
    for (seriesIndex = 0; seriesIndex <= MASK_SERIES_ID; seriesIndex++) {
      pSeriesStats = &seriesStatsTable[seriesIndex];
      if (pSeriesStats->plateCount > 0) {
        printf("seriesId %3d series %5s plateCount %4d outputCount %4d mask,a,saturated: %4d allOutputFlag %d\n",
               seriesIndex,
               pSeriesStats->series,
               pSeriesStats->plateCount,
               pSeriesStats->outputCount1,
               pSeriesStats->outputCount1 - pSeriesStats->outputCount2,
               pSeriesStats->allOutputFlag);
      }

    }
    printf("totalOutputCount %3d\n",totalOutputCount);
  }
  if (input_handle != NULL) {
    fclose(input_handle);
  }

  for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
    if (out_handle[fileIndex] != NULL) {
      fclose(out_handle[fileIndex]);
    }
  }
  if (median_handle != NULL) {
    fclose(median_handle);
  }
  if (plateTable != NULL) {
    free(plateTable);
  }

  time(&curTime);
  curTime -= startTime;

  return(EXIT_SUCCESS);
}
