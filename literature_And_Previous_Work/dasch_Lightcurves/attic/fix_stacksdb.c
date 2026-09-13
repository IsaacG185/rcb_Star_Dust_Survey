// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* fix_stacksdb.c
 *
 *  This program checks the eventnotes field for evindence of "sticky" keyboard keys and corrects any records showing this problem.
 *  
 cc -ggdb -O0  -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include fix_stacksdb.c pipelineutils.a -L /dasch/install/lib -lm -lcfitsio  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread    -ltable -lutil  -lwcs plottransientstub.o -o fix_stacksdb

 *
 *  ./fix_stacksdb -v -o losstacks.txt  > losstacks.log
 *
 * Oct  9, 2018 Edward J. Los - Initial version, adapted from update_quality.c
 *
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
#define MAX_BUFFER 10000
#define MAX_PATH 1024
#define ALLOC_INCREMENT 5000
#define MAX_RUN 5
#define CHECK_PLATECONDITION 1

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
    fprintf(outHandle,"%s%05d  unknown series\n",pPlateEntry->series,pPlateEntry->plateNumber);

    pPlateEntry->unknownSeriesFlag = 1;
  } else {
    if (GetSequesteredFlag(seriesId) == SEQUESTERED_YES) {
      pPlateEntry->sequesteredFlag = 1;
    } else {
#if 0
#if 1
      fprintf(outHandle,"%s%05d\n",pPlateEntry->series,pPlateEntry->plateNumber);
#else 
      fprintf(outHandle,"ERROR line %d new plate plateIndex %5d plateAlloc %d %s%05d\n",lineno,plateCountNew,plateAlloc,pPlateEntry->series,pPlateEntry->plateNumber);
#endif
#endif
    }
  }


  return;

}

char scratchBuffer[MAX_BUFFER];
char outputBuffer[MAX_BUFFER];

int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  char cmdchar;
  char outfile[MAX_PATH];
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
  int noteslength;
  int maxnoteslength = 0;
  char *inCharPtr;
  char inCharVal;
  char startCharVal;
  int startCharCounter;
  char *outCharPtr;
  char *outStartCharPtr;
  char outCharVal;
  int deleteflag;
  char oldCharVal;
  int totalCharsDeleted = 0;
  char* pickedStatus;
  char nostatus[] = {"No_pickedStatus"};

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
            strncpy(outfile,*++argv,MAX_PATH-2);
            if (strlen(outfile) >= MAX_PATH-3) {
              printf("ERROR: MAX_PATH too small for %s\n",*argv);
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
    printf("Usage: fix_stacksdb options\n");
    printf("  options: -v verbose\n");
    printf("           -o <output file>\n");
    return(-1);
  }

  printf("fix_stacksdb of %s %s\n",
         __DATE__,__TIME__);
  /* Now look for relevant keywords in the stacks plateevents database */
#ifdef  CHECK_PLATECONDITION
  printf("line %d Searching the platecondition table\n",__LINE__);
  res = ExecuteQuery(pStacksConnection,"SELECT series,plateNumber,date,versionId,pickedComment,boxnumber,pickednotes,valid from platecondition where pickednotes != '' order by series,plateNumber;");

#else /* CHECK_PLATECONDITION */
  printf("line %d Searching the plateevents table\n",__LINE__);
#if 0
  res = ExecuteQuery(pStacksConnection,"SELECT series,plateNumber,date,versionId,pickedStatus,stackLocation,eventnotes,valid from plateevents where eventnotes != '' and series = 'ee' and plateNumber = 1185 order by series,plateNumber;");
#else 
  res = ExecuteQuery(pStacksConnection,"SELECT series,plateNumber,date,versionId,pickedStatus,stackLocation,eventnotes,valid from plateevents where eventnotes != '' order by series,plateNumber;");
#endif
#endif /* CHECK_PLATECONDITION */
  if (!res) {
    res_ptr = mysql_store_result(pStacksConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pStacksConnection);
#ifdef CHECK_PLATECONDITION
      printf("There are currently %d rows in the platecondition table\n",numRows);
#else /* CHECK_PLATECONDITION */
      printf("There are currently %d rows in the plateevents table\n",numRows);
#endif
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber\n");
          continue;
        }
        if (sqlrow[4] != NULL) {
          pickedStatus = sqlrow[4];
        } else {
          pickedStatus = nostatus;
        }

        if (sqlrow[6] != NULL) {
          noteslength = strlen(sqlrow[6]);
          if (noteslength > maxnoteslength) {
            printf("length %d for %s%05d\n",noteslength,series,plateNumber);
            maxnoteslength = noteslength;
          }
          if (noteslength > (MAX_BUFFER-2)) {
            printf("ERROR: noteslength %d too long\n",noteslength);
            continue;
          }
          strcpy(scratchBuffer,sqlrow[6]);
          inCharPtr = scratchBuffer;
          inCharVal = *inCharPtr;
          startCharVal = 0;
          startCharCounter = 0;
          outCharPtr = outputBuffer;
          outStartCharPtr = outputBuffer;
          oldCharVal = 0;
          deleteflag = 0;
          while (inCharVal != 0) {
            if (inCharVal != oldCharVal) {
              if (startCharCounter >= MAX_RUN) {
                /* We have a run of characters.  Delete them all */
                printf("Deleting a run of %d characters '%c' from %s%05d\n",startCharCounter,oldCharVal,series,plateNumber);
                outCharPtr = outStartCharPtr;
                deleteflag = 1;
              }
              oldCharVal = inCharVal;
              outStartCharPtr= outCharPtr;
              startCharVal = inCharVal;
              startCharCounter = 1;
              *outCharPtr = *inCharPtr;
              outCharPtr++;
              inCharPtr++;
              inCharVal = *inCharPtr;
            } else {
              startCharCounter++;
              *outCharPtr = *inCharPtr;
              outCharPtr++;
              inCharPtr++;
              inCharVal = *inCharPtr;
            }
          }
          if (startCharCounter >= MAX_RUN) {
            /* We have a run of characters.  Delete them all */
            printf("Deleting a run of %d characters from %s%05d\n",startCharCounter,series,plateNumber);
            *outStartCharPtr = 0;
            deleteflag = 1;
          }
          if (strcmp(outputBuffer,"c xxc s+d") == 0) {
            outputBuffer[0] = 0;
          }
          if (deleteflag) {
#if 1
            printf("XXXXXXX '%s'\n",scratchBuffer);
#endif
#ifdef  CHECK_PLATECONDITION
            printf("Changed to '%s' input length %d output length %d for %s%05d pickedComment %s\n",outputBuffer,strlen(scratchBuffer),strlen(outputBuffer),series,plateNumber,pickedStatus);
            totalCharsDeleted += strlen(scratchBuffer)-strlen(outputBuffer);
            fprintf(outHandle,"UPDATE platecondition set pickednotes = '%s' where series = '%s' and plateNumber = %d;\n",outputBuffer,series,plateNumber);
#else /* CHECK_PLATECONDITION */
            printf("Changed to '%s' input length %d output length %d for %s%05d pickedStatus %s\n",outputBuffer,strlen(scratchBuffer),strlen(outputBuffer),series,plateNumber,pickedStatus);
            totalCharsDeleted += strlen(scratchBuffer)-strlen(outputBuffer);
            fprintf(outHandle,"UPDATE plateevents set eventnotes = '%s' where series = '%s' and plateNumber = %d;\n",outputBuffer,series,plateNumber);
#endif /* CHECK_PLATECONDITION */
          }
        }
        plateCount++;
      } 
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  } else {
    return;
  }
  printf("maxnoteslength %d plateCount %d totalCharsDeleted %d\n",maxnoteslength,plateCount,totalCharsDeleted);
  time(&curTime);
  curTime -= startTime;
  mysql_close(pConnection);
  mysql_close(pPhotConnection);
  mysql_close(pStacksConnection);
}
