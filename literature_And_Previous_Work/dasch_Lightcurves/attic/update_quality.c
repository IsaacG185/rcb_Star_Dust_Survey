// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* update_quality.c
 *
 *  This program reconciles the "quality" field in the plates table of the scanner database and the "quality" field in the 
 *  photplates table of the photometry database with all of the sources of this field from other tables, includeing the "stacks" database.
 *
 * cc -ggdb -O0  -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include update_quality.c pipelineutils.a -L /dasch/install/lib -lm -lcfitsio  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread    -ltable -lutil  -lwcs plottransientstub.o -o update_quality
 *
 *  update_quality -v -o losquality.txt >& losquality.log
 *
 *     Results as of Jan 31, 2019:
 * total plates 555059 processed photplates. 6957 unknown plates at 6 seconds
 * Totals: 49784 logbookMultiple, 123300 allMultiple, 39746 wedge for 555059 plates
 * Wedge plates accepted 16456 plates update 9482 photometry update 6070
 * Total plates: 555059 Multiple exposure plates: 431759 Single exposure plates: 123300 22.2%
 * Condition  0 totalCount 10410 uniqueCount  2132 'Grating Plate'
 * Condition  1 totalCount 39683 uniqueCount  6562 'Pickering Wedge Plate'
 * Condition  2 totalCount 31843 uniqueCount  9453 'Class indicated multiple'
 * Condition  3 totalCount 33730 uniqueCount  9948 'Multiple logbook entries'
 * Condition  4 totalCount 36905 uniqueCount  5048 'Noted during scanning or inspection'
 * Condition  5 totalCount 32047 uniqueCount  3949 'Multiple astrometry.net solutions'
 * Condition  6 totalCount 71126 uniqueCount 12063 'Close Correlation'
 * Condition  7 totalCount   331 uniqueCount     0 'Wedge Plate noted during inspection'
 * Plates with a unique multiple condition: 4915
 *
 * Apr 18, 2017 Edward J. Los - Initial version, adapted from /dasch/scanner/linux/apps/temp/UpdateQuality.cpp
 * May  6, 2017 Edward J. Los - Add totals 
 * Feb 12, 2017 Edward J. Los - provide detailed counts
 *
 *  TODO: Reconcile with UpdateQuality in photometryutils.c
 *        search of the mosaics, scans, and masks table must recognize the mosaic number and deleted mosaics. 
 */


#include <math.h>
#include "table.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>
#include <errno.h>
#include "mysql.h"
#include "errmsg.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"

#include "pipelineutils.h"

#include "photometryutils.h"
#include "magdeputil.h"
#include "daschunistd.h"

/* #define LOS_DEBUG 1 */
#define LOS_SERIES "mc"
#define LOS_PLATENUMBER 4551


#define MAX_LIST_STRING 25
#define MAX_FILENAME 512
#define MAX_BUFFER 256
#define ALLOC_INCREMENT 5000


#define DATABASE_EXPOSURES      0
#define DATABASE_PLATES         1
#define DATABASE_PLATECONDITION 2
#define DATABASE_PLATEEVENTS    3
#define DATABASE_SCANS          4
#define DATABASE_MOSAICS        5
#define DATABASE_MASKS          6
#define DATABASE_PHOTPLATES     7
#define DATABASE_MAX            8

char *databaseText[DATABASE_MAX] = 
  {
  "exposures table ???",
  "CameraPC",
  "stacks plate condition",
  "stacks plate event",
  "scans ???",
  "mosaics ???",
  "masks table ???",
  "photplates table ???"};

typedef struct _plateentry {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int seriesId;
#if 0
  int locationId;
#endif
  int galacticbin;
  char plateClass[MAX_CLASS_STRING];
  int jacketJpeg;
  int jacketJpeg2;
  int plateJpeg;
  int plateJpeg2;
  int source;
  int curplatesquality;
  int platesqualityvalid;
  int curphotquality;
  int photqualityvalid;
  int newquality;
  int updateNeeded;
  int database;   
  int rejected;  /* Not to be scanned */
  int unknownPlateFlag; /* Not in the plates database ! */
  int unknownSeriesFlag; /* Not in the series table */
  int noTranscriptionFlag; /* No transcription */
  int sequesteredFlag; /* Sequestered */
  int maskCount;
  int maskXMin;
  int maskXMax;
  int maskYMin;
  int maskYMax;
  char plateComment[MAX_COMMENT_STRING];
} PLATEENTRY,*PPLATEENTRY;

typedef struct _MULTIPLEBIT
{
  int qualityMask;
  int totalCount;
  int uniqueCount;
  const char *qualityDescr;
} MULTIPLEBIT,*PMULTIPLEBIT;

MULTIPLEBIT multipleMasks[] = {
  {QUALITY_GRATING,0,0,         "Grating Plate"},   /* Grating exposure */
  {QUALITY_WEDGE,0,0,           "Pickering Wedge Plate"},   /* Pickering Wedge plate */
  {QUALITY_MULTIPLE_CLASS,0,0,  "Class indicated multiple"}, /* bit 13 0x00002000 Class designation is multiple */
  {QUALITY_MULTIPLE_LOGBOOK,0,0,"Multiple logbook entries"},  /* bit 14 0x00004000 More than one logbook entry */
  {QUALITY_MULTIPLE_COMMENT,0,0,"Noted during scanning or inspection"}, /* bit 15 0x00008000 a comment by stacks personal */
  {QUALITY_MULTIPLE_FITWCS,0,0, "Multiple astrometry.net solutions"},  /* bit 16 0x00010000 multiple astronomy.net solutions */
  {QUALITY_MULTIPLE_MASK,0,0,   "Close Correlation"},   /* bit 17 0x00020000 is close correlation with search_close.c */
  {QUALITY_WEDGE_COMMENT,0,0,   "Wedge Plate noted during inspection"} /* bit 28 0x00040000 is noted as a wedge plate */
};

int numConditions = sizeof(multipleMasks)/sizeof(MULTIPLEBIT);

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

void NewPlate(PPLATEENTRY pPlateEntry,int plateCountNew, int plateAlloc, FILE *outHandle,int lineno) 
{
#if 0
  printf("ERROR: line %d new plate plateIndex %5d plateAlloc %d %s%05d\n",lineno,plateCountNew,plateAlloc,pPlateEntry->series,pPlateEntry->plateNumber);
#endif
  int seriesId;
  seriesId = GetSeriesId(pPlateEntry->series,0);
  if (seriesId < 0) {
#if 0    
    printf("ERROR: line %d unknown series  %s%05d\n",lineno,pPlateEntry->series,pPlateEntry->plateNumber);
#endif
    fprintf(outHandle,"%s%05d %s unknown series\n",pPlateEntry->series,pPlateEntry->plateNumber,databaseText[pPlateEntry->database]);

    pPlateEntry->unknownSeriesFlag = 1;
  } else {
    if (GetSequesteredFlag(seriesId) == SEQUESTERED_YES) {
      pPlateEntry->sequesteredFlag = 1;
    } else {
#if 0
#if 1
      fprintf(outHandle,"%s%05d %s\n",pPlateEntry->series,pPlateEntry->plateNumber,databaseText[pPlateEntry->database]);
#else 
      fprintf(outHandle,"ERROR line %d new plate plateIndex %5d plateAlloc %d %s%05d\n",lineno,plateCountNew,plateAlloc,pPlateEntry->series,pPlateEntry->plateNumber);
#endif
#endif
    }
  }


  return;

}

int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  char cmdchar;
  char outfile[MAX_BUFFER];
  FILE *outHandle;
  int verbose = 0;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  MYSQL my_stacks_connection;
  MYSQL *pStacksConnection = &my_stacks_connection;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  char *mysqlhost;
  char *username;
  char *password;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char queryString[MAX_QUERY_STRING];
  int numRows;
  int unknownPlateCount = 0;
  int plateCount = 0;
  int plateAlloc = 0;
  int plateCountNew = 0;
  int logbookMultipleCount = 0;
  int allMultipleCount = 0;
  int wedgeCount = 0;
  int plateIndex;
  int plateTestIndex;
  PPLATEENTRY pPlatesTable = NULL;
  PPLATEENTRY pPlateEntry;
#if 0
  PPLATEENTRY pTestEntry;
#endif
  PPLATEENTRY pTestEntry2;
  int result;
  int nvals;
  EXPOSURE exposureTable;
  PEXPOSURE pExposure = &exposureTable;
  time_t startTime;
  time_t curTime;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int FitWCS;
  int noTranscriptionCount = 0;
  int noLogbookEntryCount = 0;
  int wedgeAcceptedCount = 0;
  int updatePlatesCount = 0;
  int updatePhotometryCount = 0;
  int wedgeAcceptedFlag;
  PWEDGEENTRY pWedgeEntry = NULL;
  int wedgeIndex;
  double plateScale;
  double platePixels;
  double maskMaxPixels;
  double maskMinPixels;
  double maskPixels;
  int allmask = 0;       /* Bits considered */
  int qualitybits;
  double pctmask;
  int multipleIndex;
  int allPlateCount = 0;
  int nonMultipleCount = 0;
  int multipleCount = 0;
  int totalUniqueCount = 0;
  PMULTIPLEBIT pMultipleBit;


  outfile[0] = 0;
  time(&startTime);
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d\n",argstr,argc);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(outfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;



        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;



        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }



  /* Connect to the database */
  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    printf("DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }


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


  mysql_init(pStacksConnection);


  if (!mysql_real_connect(pStacksConnection,mysqlhost,username,password,"stacks",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pStacksConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pStacksConnection),mysql_error(pStacksConnection));
    }
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

  if (outfile[0] == 0) {
    printf("ERROR: no output file specified\n");
    errorFlag = 1;
  } else {
    outHandle = fopen(outfile,"wt");
    if (outHandle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to open the output file %s\n",outfile);
    } else {
      if (verbose) {
        printf("Output file %s\n",outfile);
      }
    }
  }

  if (errorFlag) {
    printf("Usage: update_quality options\n");
    printf("  options: -v verbose\n");
    printf("           -o <output file>\n");
    return(-1);
  }

  printf("update_quality of %s %s\n",
         __DATE__,__TIME__);
  /* Start by building a row of exposures, noting which ones are multiple */
  printf("line %d Searching the exposures table \n",__LINE__);
  res = ExecuteQuery(pConnection,"SELECT series,plateNumber,exposureNumber from exposures where exposureNumber < 2 order by series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pConnection);
      printf("There are currently %d + %d rows allocated to the plates table\n",numRows,ALLOC_INCREMENT);
      plateAlloc = numRows+ALLOC_INCREMENT;
      pPlatesTable = (PPLATEENTRY)calloc(plateAlloc,sizeof(PLATEENTRY));
      if (pPlatesTable == NULL) {
        printf("ERROR: failed to allocate %d PLATEENTRY rows\n",plateAlloc);
        exit(-1);
      }

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        pPlateEntry = &pPlatesTable[plateCount];
        memset(pPlateEntry,0,sizeof(PLATEENTRY));
        strcpy(pPlateEntry->series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&pPlateEntry->plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber\n");
          continue;
        }
        /* Skip illegal series and sequestered series */
        pPlateEntry->seriesId = GetSeriesId(pPlateEntry->series,0);
        if (pPlateEntry->seriesId < 0) {
          printf("ERROR: line %d unknown series  %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          continue;
        } else {
          if (GetSequesteredFlag(pPlateEntry->seriesId) == SEQUESTERED_YES) {
            continue;
          }
        }
        
#ifdef LOS_DEBUG
        if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
            (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
        }
#endif /* LOS_DEBUG */

        if (plateCount > 0) {
          pTestEntry2 = &pPlatesTable[plateCount-1];
#ifdef LOS_DEBUG
          if ((pTestEntry2->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pTestEntry2->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pTestEntry2->series,pTestEntry2->plateNumber);
          }
#endif /* LOS_DEBUG */
          if ((pTestEntry2->plateNumber == pPlateEntry->plateNumber) &&
              (strcmp(pTestEntry2->series,pPlateEntry->series) == 0)) {
            /* This is a multiple exposure */
            pTestEntry2->newquality |= QUALITY_MULTIPLE_LOGBOOK; 
          } else {
            plateCount++;
          }
        } else {
          plateCount++;
        }
        if (plateCount >= plateAlloc) {
          printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
          exit(-1);
        }

      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
#if 0
  plateTestIndex = 1; /* Temporarily a plate number */
  for (plateIndex = 0; plateIndex < plateCount; plateIndex++) {
    pPlateEntry = &pPlatesTable[plateIndex];
    if (strcmp(pPlateEntry->series,"a") != 0) {
      break;
    }
    if (plateTestIndex == pPlateEntry->plateNumber) {
      plateTestIndex++;
    } else {
      while (plateTestIndex < pPlateEntry->plateNumber) {
        fprintf(outHandle,"ERROR: No transcription for plate a%05d\n",plateTestIndex);
        plateTestIndex++;
      }
    }
  }
#endif
  time(&curTime);
  curTime -= startTime;
  printf("processed %d plates from the exposures table at %d seconds\n",plateCount,curTime);

  /* Now check the plates table */
  plateCountNew = plateCount;
  plateTestIndex = 0;
 
  printf("line %d Searching the plates table\n",__LINE__);
  res = ExecuteQuery(pConnection,"SELECT series,plateNumber,class,quality+0 from plates order by series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((plateNumber == LOS_PLATENUMBER) &&
            (strcmp(series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,series,plateNumber);
        }
#endif /* LOS_DEBUG */
        while (plateTestIndex < plateCount) {
          pPlateEntry = &pPlatesTable[plateTestIndex];
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
          result = strcmp(pPlateEntry->series,series);
          if (result == 0) {
            if (pPlateEntry->plateNumber < plateNumber) {
              result = -1;
            } else if (pPlateEntry->plateNumber > plateNumber) {
              result = +1;
            } else {
              result = 0;
            }
          }
          if (result == 0) {
            /* We have a match.  Stay on this entry in case we have duplicates */
            break;
          }
          if (result < 0) {
            /* Test entry is larger, check the next entry */
            plateTestIndex++;
          } else {
            /* Here we have an unknown plate. Create a new one */
            for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
              pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
              if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                  (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
                printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
              }
#endif /* LOS_DEBUG */
#if 0
              printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
              if ((pPlateEntry->plateNumber == plateNumber) &&
                  (strcmp(pPlateEntry->series,series) == 0)) {
                result = 0;
                break;
              }
            }
            if (plateIndex >= plateCountNew) {
              unknownPlateCount++;
              pPlateEntry = &pPlatesTable[plateCountNew];
              memset(pPlateEntry,0,sizeof(PLATEENTRY));
              strcpy(pPlateEntry->series,series);
              pPlateEntry->plateNumber = plateNumber;
              pPlateEntry->unknownPlateFlag = 1;
              pPlateEntry->database =  DATABASE_PLATES;
              NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);

              plateCountNew++;
              if (plateCountNew >= plateAlloc) {
                printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
                exit(-1);
              }
              result = 0; /* We now have a match! */          
              break;
            } else {
              break;
            }
          }            
        }
        if (result != 0) {
          for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
            pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
#if 0
            printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
            if ((pPlateEntry->plateNumber == plateNumber) &&
                (strcmp(pPlateEntry->series,series) == 0)) {
              result = 0;
              break;
            }
          }
          if (plateIndex >= plateCountNew) {
            /* Ran off the end of the list: an unknown plate */
            unknownPlateCount++;
            pPlateEntry = &pPlatesTable[plateCountNew];
            memset(pPlateEntry,0,sizeof(PLATEENTRY));
            strcpy(pPlateEntry->series,series);
            pPlateEntry->plateNumber = plateNumber;
            pPlateEntry->unknownPlateFlag = 1;
            pPlateEntry->database =  DATABASE_PLATES;
            NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
            plateCountNew++;
            if (plateCountNew >= plateAlloc) {
              printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
              exit(-1);
            }
          }
        }  
      

        if (sqlrow[2]) {
          strcpy(pPlateEntry->plateClass,sqlrow[2]);
        }
        pPlateEntry->platesqualityvalid = 1; /* O.K. here even if null */
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pPlateEntry->curplatesquality);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for quality\n");
            continue;
          }
        } 
#ifdef LOS_DEBUG
        if ((plateNumber == LOS_PLATENUMBER) &&
            (strcmp(series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,series,plateNumber);
        }
#endif /* LOS_DEBUG */
#if 0
        if ((strcmp(pPlateEntry->series,"ma") == 0) &&
            (pPlateEntry->plateNumber == 4979)) {
          printf("At plate %s%05d\n",pPlateEntry->series,pPlateEntry->plateNumber);
        }

#endif
        

        if (IsMultipleExposure(pPlateEntry->series,pPlateEntry->plateClass) != 0)  {
          pPlateEntry->newquality |= QUALITY_MULTIPLE_CLASS;
        } 
        if (IsGratingExposure(pPlateEntry->series,pPlateEntry->plateClass) != 0) {
          pPlateEntry->newquality |= QUALITY_GRATING;
        }
        if (IsColorFilterPlate(pPlateEntry->series,pPlateEntry->plateClass) != 0) {
          pPlateEntry->newquality |= QUALITY_COLOR;
        }
        if (IsSpectraPlate(pPlateEntry->series,pPlateEntry->plateClass) != 0)  {
          pPlateEntry->newquality |= QUALITY_SPECTRA;
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
        }

 
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
  time(&curTime);
  curTime -= startTime;
  printf("processed %d plates from the plates table at %d seconds\n",plateCount,curTime);
#ifdef LOS_DEBUG
  printf("line %d qsort\n",__LINE__);
#endif /* LOS_DEBUG */
  qsort((void*)pPlatesTable,plateCountNew,sizeof(PLATEENTRY),PlateCompare);
  plateCount = plateCountNew;
  plateCountNew = plateCount;
  plateTestIndex = 0;
  /* Now look for relevant keywords in the stacks platecondition database */
  printf("line %d Searching the platecondition table\n",__LINE__);
  res = ExecuteQuery(pStacksConnection,"SELECT series,plateNumber,pickedComment,pickednotes from platecondition where valid = 'yes' order by series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pStacksConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pStacksConnection);
      printf("There are currently %d rows in the platecondition table\n",numRows);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber in plateconditon\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((plateNumber == LOS_PLATENUMBER) &&
            (strcmp(series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,series,plateNumber);
        }
#endif /* LOS_DEBUG */
        while (plateTestIndex < plateCount) {
          pPlateEntry = &pPlatesTable[plateTestIndex];
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
          result = strcmp(pPlateEntry->series,series);
          if (result == 0) {
            if (pPlateEntry->plateNumber < plateNumber) {
              result = -1;
            } else if (pPlateEntry->plateNumber > plateNumber) {
              result = +1;
            } else {
              result = 0;
            }
          }
          if (result == 0) {
            /* We have a match.  Stay on this entry in case we have duplicates */
            break;
          }
          if (result < 0) {
            /* Test entry is larger, check the next entry */
            plateTestIndex++;
          } else {
            /* Here we have an unknown plate. Create a new one */
            for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
              pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
              if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                  (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
                printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
              }
#endif /* LOS_DEBUG */
#if 0
              printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
              if ((pPlateEntry->plateNumber == plateNumber) &&
                  (strcmp(pPlateEntry->series,series) == 0)) {
                result = 0;
                break;
              }
            }
            if (plateIndex >= plateCountNew) {
              unknownPlateCount++;
              pPlateEntry = &pPlatesTable[plateCountNew];
              memset(pPlateEntry,0,sizeof(PLATEENTRY));
              strcpy(pPlateEntry->series,series);
              pPlateEntry->plateNumber = plateNumber;
              pPlateEntry->unknownPlateFlag = 1;
              pPlateEntry->database =  DATABASE_PLATECONDITION;
              NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
              plateCountNew++;
              if (plateCountNew >= plateAlloc) {
                printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
                exit(-1);
              }
              result = 0; /* We now have a match! */          
              break;
            } else {
              break;
            }
          }            
        }
        if (result != 0) {
          for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
            pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
#if 0
            printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
            if ((pPlateEntry->plateNumber == plateNumber) &&
                (strcmp(pPlateEntry->series,series) == 0)) {
              result = 0;
              break;
            }
          }
          if (plateIndex >= plateCountNew) {
            /* Ran off the end of the list: an unknown plate */
            unknownPlateCount++;
            pPlateEntry = &pPlatesTable[plateCountNew];
            memset(pPlateEntry,0,sizeof(PLATEENTRY));
            strcpy(pPlateEntry->series,series);
            pPlateEntry->plateNumber = plateNumber;
            pPlateEntry->unknownPlateFlag = 1;
            pPlateEntry->database =  DATABASE_PLATECONDITION;
            NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
            plateCountNew++;
            if (plateCountNew >= plateAlloc) {
              printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
              exit(-1);
            }
          }
        }  
#ifdef LOS_DEBUG
        if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
            (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
        }
#endif /* LOS_DEBUG */

        if (sqlrow[2] != NULL) {
          if (strstr(sqlrow[2],"Multiple")) {
            pPlateEntry->newquality |= QUALITY_MULTIPLE_COMMENT;
          }
          if (strstr(sqlrow[2],"Yellow Filter")) {
            pPlateEntry->newquality |= QUALITY_COLOR;
          }
          if (strstr(sqlrow[2],"Red Filter")) {
            pPlateEntry->newquality |= QUALITY_COLOR;
          }
          if (strstr(sqlrow[2],"Fine Grating") != NULL) {
            pPlateEntry->rejected = 1;
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
          }

        }
        if (sqlrow[3] != NULL) {
          local_strlwr(sqlrow[3]);
          if (strstr(sqlrow[3],"wedge")) {
            pPlateEntry->newquality |= QUALITY_WEDGE_COMMENT;
          }
          if (strstr(sqlrow[3],"spectra")) {
            pPlateEntry->newquality |= QUALITY_SPECTRA;
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
          }
              
        }

      } 
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
  plateCount = plateCountNew;

  time(&curTime);
  curTime -= startTime;
  printf("total plates %d processed platecondition. %d unknown plates at %d seconds\n",plateCount,unknownPlateCount,curTime);
#ifdef LOS_DEBUG
  printf("line %d qsort\n",__LINE__);
#endif /* LOS_DEBUG */
  qsort((void*)pPlatesTable,plateCountNew,sizeof(PLATEENTRY),PlateCompare);
  plateCount = plateCountNew;
  plateCountNew = plateCount;
  plateTestIndex = 0;
  /* Now look for relevant keywords in the stacks plateevents database */
  printf("line %d Searching the plateevents table\n",__LINE__);
  res = ExecuteQuery(pStacksConnection,"SELECT series,plateNumber,pickedStatus from plateevents where valid = 'yes' order by series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pStacksConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pStacksConnection);
      printf("There are currently %d rows in the plateevents table\n",numRows);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber in plateconditon\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((plateNumber == LOS_PLATENUMBER) &&
            (strcmp(series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,series,plateNumber);
        }
#endif /* LOS_DEBUG */
        while (plateTestIndex < plateCount) {
          pPlateEntry = &pPlatesTable[plateTestIndex];
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
          result = strcmp(pPlateEntry->series,series);
          if (result == 0) {
            if (pPlateEntry->plateNumber < plateNumber) {
              result = -1;
            } else if (pPlateEntry->plateNumber > plateNumber) {
              result = +1;
            } else {
              result = 0;
            }
          }
          if (result == 0) {
            /* We have a match.  Stay on this entry in case we have duplicates */
            break;
          }
          if (result < 0) {
            /* Test entry is larger, check the next entry */
            plateTestIndex++;
          } else {
            /* Here we have an unknown plate. Create a new one */
            for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
              pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
              if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                  (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
                printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
              }
#endif /* LOS_DEBUG */
#if 0
              printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
              if ((pPlateEntry->plateNumber == plateNumber) &&
                  (strcmp(pPlateEntry->series,series) == 0)) {
                result = 0;
                break;
              }
            }
            if (plateIndex >= plateCountNew) {
              unknownPlateCount++;
              pPlateEntry = &pPlatesTable[plateCountNew];
              memset(pPlateEntry,0,sizeof(PLATEENTRY));
              strcpy(pPlateEntry->series,series);
              pPlateEntry->plateNumber = plateNumber;
              pPlateEntry->unknownPlateFlag = 1;
              pPlateEntry->database =  DATABASE_PLATEEVENTS;
              NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
              plateCountNew++;
              if (plateCountNew >= plateAlloc) {
                printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
                exit(-1);
              }
              result = 0; /* We now have a match! */          
              break;
            } else {
              break;
            }
          }            
        }
        if (result != 0) {
          for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
            pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
#if 0
            printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
            if ((pPlateEntry->plateNumber == plateNumber) &&
                (strcmp(pPlateEntry->series,series) == 0)) {
              result = 0;
              break;
            }
          }
          if (plateIndex >= plateCountNew) {
            /* Ran off the end of the list: an unknown plate */
            unknownPlateCount++;
            pPlateEntry = &pPlatesTable[plateCountNew];
            memset(pPlateEntry,0,sizeof(PLATEENTRY));
            strcpy(pPlateEntry->series,series);
            pPlateEntry->plateNumber = plateNumber;
            pPlateEntry->unknownPlateFlag = 1;
            pPlateEntry->database =  DATABASE_PLATEEVENTS;
            NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
            plateCountNew++;
            if (plateCountNew >= plateAlloc) {
              printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
              exit(-1);
            }
          }
        }  
#ifdef LOS_DEBUG
        if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
            (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
        }
#endif /* LOS_DEBUG */
        if (strstr(sqlrow[2],"No Transcription") != NULL) {
          pPlateEntry->noTranscriptionFlag = 1;
        }
        if ((strstr(sqlrow[2],"Discarded") != NULL) ||
            (strstr(sqlrow[2],"Rejected") != NULL) ||
            (strstr(sqlrow[2],"Resseau Grid") != NULL) ||
            (strstr(sqlrow[2],"Reseau Grid") != NULL) ||
            (strstr(sqlrow[2],"Severe Delamination") != NULL) ||
            (strstr(sqlrow[2],"Fine Grating") != NULL) ||
            (strstr(sqlrow[2],"Long Trails") != NULL) ||
            (strstr(sqlrow[2],"Many Exposures") != NULL) ||
            (strstr(sqlrow[2],"Multi Direction Trailing") != NULL) ||
            (strstr(sqlrow[2],"Too Thick") != NULL) ||
            (strstr(sqlrow[2],"Too Dark") != NULL)) {
          pPlateEntry->rejected = 1;
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
        }

        if (sqlrow[2] != NULL) {
          if (strstr(sqlrow[2],"Spectra")) {
            pPlateEntry->newquality |= QUALITY_SPECTRA;
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
          }
          if (strstr(sqlrow[2],"Grating")) {
            pPlateEntry->newquality |= QUALITY_GRATING;
          }
          if (strstr(sqlrow[2],"Many Exposures")) {
            pPlateEntry->newquality |= QUALITY_MULTIPLE_COMMENT;
          }

        }

      } 
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
  plateCount = plateCountNew;

  time(&curTime);
  curTime -= startTime;
  printf("total plates %d processed plateevents. %d unknown plates at %d seconds\n",plateCount,unknownPlateCount,curTime);
#ifdef LOS_DEBUG
  printf("line %d qsort\n",__LINE__);
#endif /* LOS_DEBUG */
  qsort((void*)pPlatesTable,plateCountNew,sizeof(PLATEENTRY),PlateCompare);
  plateCount = plateCountNew;
  plateCountNew = plateCount;
  plateTestIndex = 0;
  /* Now look for notes in the scans database */
  printf("line %d Searching the scans table\n",__LINE__);
  res = ExecuteQuery(pConnection,"SELECT series,plateNumber,scanComment from scans where series != 'mask' and series != '' and series != 'hale' and ((scanComment is not null) and (scanComment != '')) order by series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pConnection);
      printf("There are currently %d rows in the scans table\n",numRows);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber in plateconditon\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((plateNumber == LOS_PLATENUMBER) &&
            (strcmp(series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,series,plateNumber);
        }
#endif /* LOS_DEBUG */
        while (plateTestIndex < plateCount) {
          pPlateEntry = &pPlatesTable[plateTestIndex];
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
          result = strcmp(pPlateEntry->series,series);
          if (result == 0) {
            if (pPlateEntry->plateNumber < plateNumber) {
              result = -1;
            } else if (pPlateEntry->plateNumber > plateNumber) {
              result = +1;
            } else {
              result = 0;
            }
          }
          if (result == 0) {
            /* We have a match.  Stay on this entry in case we have duplicates */
            break;
          }
          if (result < 0) {
            /* Test entry is larger, check the next entry */
            plateTestIndex++;
          } else {
            /* Here we have an unknown plate. Create a new one */
            for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
              pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
              if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                  (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
                printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
              }
#endif /* LOS_DEBUG */
#if 0
              printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
              if ((pPlateEntry->plateNumber == plateNumber) &&
                  (strcmp(pPlateEntry->series,series) == 0)) {
                result = 0;
                break;
              }
            }
            if (plateIndex >= plateCountNew) {
              unknownPlateCount++;
              pPlateEntry = &pPlatesTable[plateCountNew];
              memset(pPlateEntry,0,sizeof(PLATEENTRY));
              strcpy(pPlateEntry->series,series);
              pPlateEntry->plateNumber = plateNumber;
              pPlateEntry->unknownPlateFlag = 1;
              pPlateEntry->database =  DATABASE_SCANS;
              NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
              plateCountNew++;
              if (plateCountNew >= plateAlloc) {
                printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
                exit(-1);
              }
              result = 0; /* We now have a match! */          
              break;
            } else {
              break;
            }
          }            
        }
        if (result != 0) {
          for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
            pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
#if 0
            printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
            if ((pPlateEntry->plateNumber == plateNumber) &&
                (strcmp(pPlateEntry->series,series) == 0)) {
              result = 0;
              break;
            }
          }
          if (plateIndex >= plateCountNew) {
            /* Ran off the end of the list: an unknown plate */
            unknownPlateCount++;
            pPlateEntry = &pPlatesTable[plateCountNew];
            memset(pPlateEntry,0,sizeof(PLATEENTRY));
            strcpy(pPlateEntry->series,series);
            pPlateEntry->plateNumber = plateNumber;
            pPlateEntry->unknownPlateFlag = 1;
            pPlateEntry->database =  DATABASE_SCANS;
            NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
            plateCountNew++;
            if (plateCountNew >= plateAlloc) {
              printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
              exit(-1);
            }
          }
        }  
#ifdef LOS_DEBUG
        if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
            (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
        }
#endif /* LOS_DEBUG */
        if (sqlrow[2] != NULL) {
          local_strlwr(sqlrow[2]);
          if (strstr(sqlrow[2],"multiple")) {
            pPlateEntry->newquality |= QUALITY_MULTIPLE_COMMENT;
          }
          if (strstr(sqlrow[2],"grating")) {
            pPlateEntry->newquality |= QUALITY_GRATING;
          }
          if (strstr(sqlrow[2],"specta")) {
            pPlateEntry->newquality |= QUALITY_SPECTRA;
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
          }
          if (strstr(sqlrow[2],"red filter")) {
            pPlateEntry->newquality |= QUALITY_COLOR;
          }
          if (strstr(sqlrow[2],"yellow filter")) {
            pPlateEntry->newquality |= QUALITY_COLOR;
          }


        }
    

      } 
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
  plateCount = plateCountNew;

  time(&curTime);
  curTime -= startTime;
  printf("total plates %d processed scans. %d unknown plates at %d seconds\n",plateCount,unknownPlateCount,curTime);
#ifdef LOS_DEBUG
  printf("line %d qsort\n",__LINE__);
#endif /* LOS_DEBUG */
  qsort((void*)pPlatesTable,plateCountNew,sizeof(PLATEENTRY),PlateCompare);
  plateCount = plateCountNew;
  plateCountNew = plateCount;
  plateTestIndex = 0;
  /* Now look for state in the mosaics database */
  printf("line %d Searching the mosaics table\n",__LINE__);
  res = ExecuteQuery(pConnection,"SELECT series,plateNumber,FitWCS+0 from mosaics where series != 'mask' and series != 'hale' and solutionNumber = 0 order by series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pConnection);
      printf("There are currently %d rows in the mosaics table\n",numRows);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber in plateconditon\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((plateNumber == LOS_PLATENUMBER) &&
            (strcmp(series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,series,plateNumber);
        }
#endif /* LOS_DEBUG */
        while (plateTestIndex < plateCount) {
          pPlateEntry = &pPlatesTable[plateTestIndex];
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
          result = strcmp(pPlateEntry->series,series);
          if (result == 0) {
            if (pPlateEntry->plateNumber < plateNumber) {
              result = -1;
            } else if (pPlateEntry->plateNumber > plateNumber) {
              result = +1;
            } else {
              result = 0;
            }
          }
          if (result == 0) {
            /* We have a match.  Stay on this entry in case we have duplicates */
            break;
          }
          if (result < 0) {
            /* Test entry is larger, check the next entry */
            plateTestIndex++;
          } else {
            /* Here we have an unknown plate. Create a new one */
            for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
              pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
              if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                  (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
                printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
              }
#endif /* LOS_DEBUG */
#if 0
              printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
              if ((pPlateEntry->plateNumber == plateNumber) &&
                  (strcmp(pPlateEntry->series,series) == 0)) {
                result = 0;
                break;
              }
            }
            if (plateIndex >= plateCountNew) {
              unknownPlateCount++;
              pPlateEntry = &pPlatesTable[plateCountNew];
              memset(pPlateEntry,0,sizeof(PLATEENTRY));
              strcpy(pPlateEntry->series,series);
              pPlateEntry->plateNumber = plateNumber;
              pPlateEntry->unknownPlateFlag = 1;
              pPlateEntry->database =  DATABASE_MOSAICS;
              NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
              plateCountNew++;
              if (plateCountNew >= plateAlloc) {
                printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
                exit(-1);
              }
              result = 0; /* We now have a match! */          
              break;
            } else {
              break;
            }
          }            
        }
        if (result != 0) {
          for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
            pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
#if 0
            printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
            if ((pPlateEntry->plateNumber == plateNumber) &&
                (strcmp(pPlateEntry->series,series) == 0)) {
              result = 0;
              break;
            }
          }
          if (plateIndex >= plateCountNew) {
            /* Ran off the end of the list: an unknown plate */
            unknownPlateCount++;
            pPlateEntry = &pPlatesTable[plateCountNew];
            memset(pPlateEntry,0,sizeof(PLATEENTRY));
            strcpy(pPlateEntry->series,series);
            pPlateEntry->plateNumber = plateNumber;
            pPlateEntry->unknownPlateFlag = 1;
            pPlateEntry->database =  DATABASE_MOSAICS;
            NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
            plateCountNew++;
            if (plateCountNew >= plateAlloc) {
              printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
              exit(-1);
            }
          }
        }  

        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&FitWCS);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for FitWCS\n");
            continue;
          }
        } else {
          FitWCS = 0;
        }

        if ((FitWCS & FITWCS_MULTSUCCEEDED) != 0) {
          pPlateEntry->newquality |=  QUALITY_MULTIPLE_FITWCS;
        }
        if ((FitWCS & FITWCS_WEDGE) != 0) {
          pPlateEntry->newquality |= QUALITY_WEDGE;
        }
      
 

      } 
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
  plateCount = plateCountNew;

  time(&curTime);
  curTime -= startTime;
  printf("total plates %d processed mosaics. %d unknown plates at %d seconds\n",plateCount,unknownPlateCount,curTime);
#ifdef LOS_DEBUG
  printf("line %d qsort\n",__LINE__);
#endif /* LOS_DEBUG */
  qsort((void*)pPlatesTable,plateCountNew,sizeof(PLATEENTRY),PlateCompare);
  plateCount = plateCountNew;
  plateCountNew = plateCount;
  plateTestIndex = 0;
  /* Now look for state in the masks database */
  printf("line %d Searching the masks table\n",__LINE__);
  res = ExecuteQuery(pConnection,"SELECT series,plateNumber,maskCount,maskXMin,maskXMax,maskYMin,maskYMax from masks where series != 'mask' and series != 'hale' and maskIndex = 1 order by series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pConnection);
      printf("There are currently %d rows in the masks table\n",numRows);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber in plateconditon\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((plateNumber == LOS_PLATENUMBER) &&
            (strcmp(series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,series,plateNumber);
        }
#endif /* LOS_DEBUG */
        while (plateTestIndex < plateCount) {
          pPlateEntry = &pPlatesTable[plateTestIndex];
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
          result = strcmp(pPlateEntry->series,series);
          if (result == 0) {
            if (pPlateEntry->plateNumber < plateNumber) {
              result = -1;
            } else if (pPlateEntry->plateNumber > plateNumber) {
              result = +1;
            } else {
              result = 0;
            }
          }
          if (result == 0) {
            /* We have a match.  Stay on this entry in case we have duplicates */
            break;
          }
          if (result < 0) {
            /* Test entry is larger, check the next entry */
            plateTestIndex++;
          } else {
            /* Here we have an unknown plate. Create a new one */
            for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
              pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
              if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                  (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
                printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
              }
#endif /* LOS_DEBUG */
#if 0
              printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
              if ((pPlateEntry->plateNumber == plateNumber) &&
                  (strcmp(pPlateEntry->series,series) == 0)) {
                result = 0;
                break;
              }
            }
            if (plateIndex >= plateCountNew) {
              unknownPlateCount++;
              pPlateEntry = &pPlatesTable[plateCountNew];
              memset(pPlateEntry,0,sizeof(PLATEENTRY));
              strcpy(pPlateEntry->series,series);
              pPlateEntry->plateNumber = plateNumber;
              pPlateEntry->unknownPlateFlag = 1;
              pPlateEntry->database =  DATABASE_MASKS;
              NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
              plateCountNew++;
              if (plateCountNew >= plateAlloc) {
                printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
                exit(-1);
              }
              result = 0; /* We now have a match! */          
              break;
            } else {
              break;
            }
          }            
        }
        if (result != 0) {
          for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
            pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
#if 0
            printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
            if ((pPlateEntry->plateNumber == plateNumber) &&
                (strcmp(pPlateEntry->series,series) == 0)) {
              result = 0;
              break;
            }
          }
          if (plateIndex >= plateCountNew) {
            /* Ran off the end of the list: an unknown plate */
            unknownPlateCount++;
            pPlateEntry = &pPlatesTable[plateCountNew];
            memset(pPlateEntry,0,sizeof(PLATEENTRY));
            strcpy(pPlateEntry->series,series);
            pPlateEntry->plateNumber = plateNumber;
            pPlateEntry->unknownPlateFlag = 1;
            pPlateEntry->database =  DATABASE_MASKS;
            NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
            plateCountNew++;
            if (plateCountNew >= plateAlloc) {
              printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
              exit(-1);
            }
          }
        }  
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pPlateEntry->maskCount);
          if (pPlateEntry->maskCount > 0) {
            pPlateEntry->newquality |=  QUALITY_MULTIPLE_MASK;
            if (nvals != 1) {
              printf("ERROR: nvals is %d for maskCount\n");
              continue;
            }

            if (sqlrow[3]) {
              nvals = sscanf(sqlrow[3],"%d",&pPlateEntry->maskXMin);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for maskXMin\n");
                continue;
              }
            }
            if (sqlrow[4]) {
              nvals = sscanf(sqlrow[4],"%d",&pPlateEntry->maskXMax);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for maskXMax\n");
                continue;
              }
            }
            if (sqlrow[5]) {
              nvals = sscanf(sqlrow[5],"%d",&pPlateEntry->maskYMin);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for maskYmin\n");
                continue;
              }
            }
            if (sqlrow[6]) {
              nvals = sscanf(sqlrow[6],"%d",&pPlateEntry->maskYMax);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for maskYMax\n");
                continue;
              }
            }

          }
        }

        
 

      } 
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
  plateCount = plateCountNew;

  time(&curTime);
  curTime -= startTime;
  printf("total plates %d processed masks. %d unknown plates at %d seconds\n",plateCount,unknownPlateCount,curTime);
#ifdef LOS_DEBUG
  printf("line %d qsort\n",__LINE__);
#endif /* LOS_DEBUG */
  qsort((void*)pPlatesTable,plateCountNew,sizeof(PLATEENTRY),PlateCompare);
  plateCount = plateCountNew;
  plateCountNew = plateCount;
  plateTestIndex = 0;
  /* Now look for the quality in the photplates database  */
  printf("line %d Searching the photplates table\n",__LINE__);
  res = ExecuteQuery(pPhotConnection,"SELECT series,plateNumber,quality+0 from photplates inner join photseries using(seriesId) order by series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pPhotConnection);
      printf("There are currently %d rows in the photplates table\n",numRows);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber in photplates\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((plateNumber == LOS_PLATENUMBER) &&
            (strcmp(series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,series,plateNumber);
        }
#endif /* LOS_DEBUG */
        while (plateTestIndex < plateCount) {
          pPlateEntry = &pPlatesTable[plateTestIndex];
#ifdef LOS_DEBUG
          if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
          result = strcmp(pPlateEntry->series,series);
          if (result == 0) {
            if (pPlateEntry->plateNumber < plateNumber) {
              result = -1;
            } else if (pPlateEntry->plateNumber > plateNumber) {
              result = +1;
            } else {
              result = 0;
            }
          }
          if (result == 0) {
            /* We have a match.  Stay on this entry in case we have duplicates */
            break;
          }
          if (result < 0) {
            /* Test entry is larger, check the next entry */
            plateTestIndex++;
          } else {
            /* Here we have an unknown plate. Create a new one */
            for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
              pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
              if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                  (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
                printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
              }
#endif /* LOS_DEBUG */
#if 0
              printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
              if ((pPlateEntry->plateNumber == plateNumber) &&
                  (strcmp(pPlateEntry->series,series) == 0)) {
                result = 0;
                break;
              }
            }
            if (plateIndex >= plateCountNew) {
              unknownPlateCount++;
              pPlateEntry = &pPlatesTable[plateCountNew];
              memset(pPlateEntry,0,sizeof(PLATEENTRY));
              strcpy(pPlateEntry->series,series);
              pPlateEntry->plateNumber = plateNumber;
              pPlateEntry->unknownPlateFlag = 1;
              pPlateEntry->database =  DATABASE_PHOTPLATES;
              NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
              plateCountNew++;
              if (plateCountNew >= plateAlloc) {
                printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
                exit(-1);
              }
              result = 0; /* We now have a match! */          
              break;
            } else {
              break;
            }
          }            
        }
        if (result != 0) {
          for (plateIndex = plateCount; plateIndex < plateCountNew; plateIndex++) {
            pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
            if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
                (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
              printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
            }
#endif /* LOS_DEBUG */
#if 0
            printf("line %d testing %s%05d and %s%05d at plateCountNew %d\n",__LINE__,series,plateNumber,pPlateEntry->series,pPlateEntry->plateNumber,plateCountNew);
#endif
            if ((pPlateEntry->plateNumber == plateNumber) &&
                (strcmp(pPlateEntry->series,series) == 0)) {
              result = 0;
              break;
            }
          }
          if (plateIndex >= plateCountNew) {
            /* Ran off the end of the list: an unknown plate */
            unknownPlateCount++;
            pPlateEntry = &pPlatesTable[plateCountNew];
            memset(pPlateEntry,0,sizeof(PLATEENTRY));
            strcpy(pPlateEntry->series,series);
            pPlateEntry->plateNumber = plateNumber;
            pPlateEntry->unknownPlateFlag = 1;
            pPlateEntry->database =  DATABASE_PHOTPLATES;
            NewPlate(pPlateEntry,plateCountNew,plateAlloc,outHandle,__LINE__);
            plateCountNew++;
            if (plateCountNew >= plateAlloc) {
              printf("ERROR: plateAlloc  %d is too small\n",plateAlloc);
              exit(-1);
            }
          }
        }  
        pPlateEntry->photqualityvalid = 1; /* OK here even if null */
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pPlateEntry->curphotquality);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for quality\n");
            continue;
          }
        }

        
 

      } 
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
  plateCount = plateCountNew;

  time(&curTime);
  curTime -= startTime;
  printf("total plates %d processed photplates. %d unknown plates at %d seconds\n",plateCount,unknownPlateCount,curTime);
#ifdef LOS_DEBUG
  printf("line %d qsort\n",__LINE__);
#endif /* LOS_DEBUG */
  qsort((void*)pPlatesTable,plateCountNew,sizeof(PLATEENTRY),PlateCompare);
  plateCountNew = plateCount;


  /* get counts of multiples */
  for (plateIndex = 0; plateIndex < plateCount; plateIndex++) {
    pPlateEntry = &pPlatesTable[plateIndex];
    if ((pPlateEntry->newquality & (QUALITY_MULTIPLE_CLASS|QUALITY_MULTIPLE_LOGBOOK)) != 0) {
      logbookMultipleCount++;
    }
    if ((pPlateEntry->newquality & (QUALITY_MULTIPLE_CLASS|QUALITY_MULTIPLE_LOGBOOK|QUALITY_MULTIPLE_COMMENT|QUALITY_MULTIPLE_FITWCS|QUALITY_MULTIPLE_MASK|QUALITY_MULTIPLE|QUALITY_GRATING|QUALITY_WEDGE|QUALITY_WEDGE_COMMENT)) != 0) {
      allMultipleCount++;
    }
    if ((pPlateEntry->newquality & (QUALITY_WEDGE|QUALITY_WEDGE_COMMENT)) != 0) {
      wedgeCount++;
    }
  }
  printf("Totals: %d logbookMultiple, %d allMultiple, %d wedge for %d plates\n",logbookMultipleCount,allMultipleCount,wedgeCount,plateCount);

#if 1  /* Uncomment this block to show necessary updates */
  for (plateIndex = 0; plateIndex < plateCount; plateIndex++) {
    pPlateEntry = &pPlatesTable[plateIndex];
    if (((pPlateEntry->newquality & QUALITY_WEDGE) != 0) &&
        (pPlateEntry->maskCount <= 1) && 
        ((pPlateEntry->newquality & (QUALITY_GRATING | QUALITY_MULTIPLE_LOGBOOK | QUALITY_SPECTRA | QUALITY_MULTIPLE_CLASS | QUALITY_MULTIPLE_FITWCS)) == 0)) {
      for (wedgeIndex = 0; wedgeIndex < wedgeTableSize; wedgeIndex++) {
        pWedgeEntry = &wedgeTable[wedgeIndex];
        if (strlen(pWedgeEntry->series) == 0) {
          break;
        }
        if (strcmp(pWedgeEntry->series,pPlateEntry->series) == 0) {
          wedgeAcceptedFlag = 0;
          maskMinPixels = 0;
          maskMaxPixels = 0;
          maskPixels = 0;
          if (pPlateEntry->maskCount == 0) {
            wedgeAcceptedFlag = 1;
          } else {

            plateScale = GetFittedPlateScale(pPlateEntry->seriesId,pPlateEntry->plateNumber); /* Plate scale is in degrees per pixel */
            platePixels = ((pWedgeEntry->separation/3600.0)/plateScale);
            maskMaxPixels = sqrt(sqr(pPlateEntry->maskXMax) + sqr(pPlateEntry->maskYMax));
            maskMinPixels = maskMaxPixels;
            maskPixels = sqrt(sqr(pPlateEntry->maskXMin) + sqr(pPlateEntry->maskYMax));
            if (maskPixels > maskMaxPixels) {
              maskMaxPixels = maskPixels; 
            }
            if (maskPixels < maskMinPixels) {
              maskMinPixels = maskPixels;
            }
            maskPixels = sqrt(sqr(pPlateEntry->maskXMax) + sqr(pPlateEntry->maskYMin));
            if (maskPixels > maskMaxPixels) {
              maskMaxPixels = maskPixels;
            }
            if (maskPixels < maskMinPixels) {
              maskMinPixels = maskPixels;
            }
            maskPixels = sqrt(sqr(pPlateEntry->maskXMin) + sqr(pPlateEntry->maskYMin));
            if (maskPixels > maskMaxPixels) {
              maskMaxPixels = maskPixels;
            }
            if (maskPixels < maskMinPixels) {
              maskMinPixels = maskPixels;
            }
            if ((platePixels < maskMaxPixels) &&
                (platePixels > maskMinPixels)) {
              wedgeAcceptedFlag = 1;
            }
          }
          if (wedgeAcceptedFlag) {
            wedgeAcceptedCount++;
            fprintf(outHandle,"%s%05d wedge plate newquality %6d 0x%06x maskCount %d Min %.1f Act %.1f Max %.1f\n",
                    pPlateEntry->series,
                    pPlateEntry->plateNumber,
                    pPlateEntry->newquality,
                    pPlateEntry->newquality,
                    pPlateEntry->maskCount,
                    maskMinPixels,
                    platePixels,
                    maskMaxPixels);
            if ((pPlateEntry->platesqualityvalid != 0) &&
                ((pPlateEntry->newquality &  DATABASE_QUALITY_MASK) != pPlateEntry->curplatesquality)) {
              fprintf(outHandle,"ERROR: plates quality is 0x%x and should be 0x%x in %s%05d\n",
                      pPlateEntry->curplatesquality,
                      (pPlateEntry->newquality &  DATABASE_QUALITY_MASK),
                      pPlateEntry->series,
                      pPlateEntry->plateNumber);
              updatePlatesCount++;
              fprintf(outHandle,"UPDATE plates set quality = %d where series = '%s' and plateNumber = %d;\n",(pPlateEntry->newquality &  DATABASE_QUALITY_MASK),pPlateEntry->series,pPlateEntry->plateNumber);
            }

            if ((pPlateEntry->photqualityvalid != 0) && 
                ((pPlateEntry->newquality &  DATABASE_QUALITY_MASK) != pPlateEntry->curphotquality)) {
              fprintf(outHandle,"ERROR: photplates quality is 0x%x and should be 0x%x in %s%05d\n",
                      pPlateEntry->curphotquality,
                      (pPlateEntry->newquality &  DATABASE_QUALITY_MASK),
                      pPlateEntry->series,
                      pPlateEntry->plateNumber
                      );
              updatePhotometryCount++;
              fprintf(outHandle,"UPDATE photplates set quality = %d where seriesId = %d and plateNumber = %d;\n",(pPlateEntry->newquality &  DATABASE_QUALITY_MASK),pPlateEntry->seriesId,pPlateEntry->plateNumber);
            }


            break;
          }
        }
      }
    }
  }
  printf("Wedge plates accepted %d plates update %d photometry update %d\n",wedgeAcceptedCount,updatePlatesCount,updatePhotometryCount);
#endif
#if 1  /* Uncomment this block to display multiple condition totals */
  allmask = 0;
  for (multipleIndex = 0; multipleIndex < numConditions; multipleIndex++) {
    pMultipleBit = &multipleMasks[multipleIndex];
    allmask |= pMultipleBit->qualityMask;
  }
  for (plateIndex = 0; plateIndex < plateCount; plateIndex++) {
    pPlateEntry = &pPlatesTable[plateIndex];
    allPlateCount++;
    if ((pPlateEntry->newquality & allmask) == 0) {
      nonMultipleCount++;
    } else {
      multipleCount++;
      qualitybits = (pPlateEntry->newquality & allmask);
      for (multipleIndex = 0; multipleIndex < numConditions; multipleIndex++) {
        pMultipleBit = &multipleMasks[multipleIndex];
        if ((qualitybits & pMultipleBit->qualityMask) != 0) {
          pMultipleBit->totalCount++;
          if (pMultipleBit->qualityMask == qualitybits) {
            pMultipleBit->uniqueCount++;
            totalUniqueCount++;
          }

        }

      }

    }
  }
  pctmask = (100.*multipleCount)/(1.0*allPlateCount);
  printf("Total plates: %d Single exposure plates: %d Multiple exposure plates: %d %.1f%%\n",allPlateCount,nonMultipleCount,multipleCount,pctmask); 
  for (multipleIndex = 0; multipleIndex < numConditions; multipleIndex++) {
    pMultipleBit = &multipleMasks[multipleIndex];
    printf("Condition %2d totalCount %5d uniqueCount %5d '%s'\n",multipleIndex,pMultipleBit->totalCount,pMultipleBit->uniqueCount,pMultipleBit->qualityDescr);
  }
  printf("Plates with a unique multiple condition: %d\n",totalUniqueCount);

#endif
#if 0  /* Uncomment this block to display a table of non-transcribed plates */
  for (plateIndex = 0; plateIndex < plateCount; plateIndex++) {
    pPlateEntry = &pPlatesTable[plateIndex];
#ifdef LOS_DEBUG
    if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
        (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
      printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
    }
#endif /* LOS_DEBUG */
    if (pPlateEntry->sequesteredFlag != 0) {
      continue;
    }
    if (pPlateEntry->unknownPlateFlag == 0) {
      continue;
    }
    if (pPlateEntry->unknownSeriesFlag != 0) {
      continue;
    }
    if ((pPlateEntry->newquality & QUALITY_SPECTRA) != 0) {
      continue;
    }
    if (pPlateEntry->rejected != 0) {
      continue;
    }
    noLogbookEntryCount++;
    if (pPlateEntry->noTranscriptionFlag != 0) {
      noTranscriptionCount++;
      fprintf(outHandle,"%s%05d %s\n",pPlateEntry->series,pPlateEntry->plateNumber,databaseText[pPlateEntry->database]);
    } else {
      fprintf(outHandle,"%s%05d %s pickedStatus does not include 'No Transcription' \n",pPlateEntry->series,pPlateEntry->plateNumber,databaseText[pPlateEntry->database]);      
    }
  }
  printf("noLogbookEntryCount %d noTranscriptionCount %d\n",noLogbookEntryCount,noTranscriptionCount);
#endif
  mysql_close(pConnection);
  mysql_close(pPhotConnection);
  mysql_close(pStacksConnection);
}
