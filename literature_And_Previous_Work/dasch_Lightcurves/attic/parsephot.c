// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* parsephot.c
 
 *  This program searches old phot_*.log files and populates the mosaicNumber<n> fields of photplates
 *
 *  Usage:  parsephot -l  <list of logs>  test with run_932.log for i24330
 *
 *           parsephot  -v  -l /home/scanner/junk/photlogs.list -o /home/scanner/Pipeline/los.sql -p /home/scanner/Pipeline/lostotal.list  -t /home/scanner/Pipeline/total.list
 *
 *           parsephot -s -t photometry.tmp to check for duplicates
 *           parsephot  -j -s -t /home/scanner/junk/JPEG/los22.tmp to check for JPEG duplicates
 *           parsephot -k -s -t  /home/scanner/Pipeline/lossearch_close3.list -m  /home/scanner/Pipeline/lossearch_close3.out

 * -p file format: <plate from logs|NOPHOT> <plate from total.list|NOTOTAL> <AGREE|DIFFER|UNK> [<catalog>:<mosaicNumber>:<versionNumber>[X for failed photometry][,REPROCESS for deleted mosaic source]
 *


cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I /dasch/install/include  -L /dasch/install/lib -lm parsephot.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread plottransientstub.o -o parsephot


awk '{print " " $1}' lostotal.list | sort > los.tmp
awk '{print " " $1}' total.list | sort > los1.tmp
diff los.tmp los1.tmp | grep ">" | wc
diff los.tmp los1.tmp | grep "<" | wc

TODO:
ERROR: Failed to find spatial bins file /dasch/Pipeline/ingest/x16422_00_01r90ww_apass.out.spatial_bins.db

 *
 * Oct  8, 2019 Edward J. Los - Initial version  
 * Oct 14, 2019 Edward J. Los - Handle 'table_load: missing column: rejectFlag' lines
 *                              Add total.list produced by dasch_scat
 * Oct 15, 2019 Edward J. Los - Do not output sequestered series
 * Oct 18, 2019 Edward J. Los - Allow check of total.list only
 * Oct 20, 2019 Edward J. Los - Add -s to process total.list entries of with no binning nor rotation
 * Mar  3, 2020 Edward J. Los - Add -j to process jpegissue.txt
 * Apr 15, 2020 Edward J. Los - Add -m to produce a merged file and -k to parse records with only plate and mosaic numbers.
 * May 29, 2020 Edward J. Los - Integrate with MySQL to correct any discrepancies in the photplates database.
 *                            - Integrate with BuildMosaicList() to get the current best mosaics. (deprecate -t option)
 * Jun  1, 2020 Edward J. Los - Correct jpegissue plate and jacket duplicate entry handling
 * Jun  5, 2020 Edward J. Los - Incorporate the plate lists from resort_magfiles
 * Aug 11, 2020 Edward J. Los - Skip further output in the -s mode.
 * Sep 16, 2020 Edward J. Los - Change list_name to loglist_name and list_handle to loglist_handle.  
 *                              Add -s to look at the log file only.  Remove enough of DEBUG_LOGS_ONLY to make this work.
 */   


#include <math.h>
#include <string.h>
#include "mysql.h"
#include "table.h"
#include "time.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER     512
#define MAX_TIME_LENGTH 28
#define JPEG_TYPE_NONE 0
#define JPEG_TYPE_PLATE 1
#define JPEG_TYPE_JACKET 2

#define SKIP_GAIA 1
/* #define DEBUG_PHOTINSERT 1 */
#ifdef DEBUG_PHOTINSERT
#define DEBUG_SERIES "ac"
#define DEBUG_SERIESID 2
#endif /* DEBUG_PHOTINSERT */
/* #define DEBUG_NO_LOGS 1 */

/* Define the following to look at all of the phot*.log files to date and populate the photinsert table */
/* Undefine the following to read the photinsert table and decide what needs to be flushed from the photometry database and what needs to be reprocessed */
/* #define DEBUG_LOGS_ONLY 1  */


#define ERROR_01  1 /*   */
#define ERROR_02  2 /*   */
#define ERROR_03  3 /*   */
#define ERROR_04  4 /*   */
#define ERROR_05  5 /*   */
#define ERROR_06  6 /*   */
#define ERROR_07  7 /*   */
#define ERROR_08  8 /*   */
#define ERROR_09  9 /*   */
#define ERROR_10 10 /*   */
#define ERROR_11 11 /*   */
#define ERROR_12 12 /*   */
#define ERROR_13 13 /*   */
#define ERROR_14 14 /*   */
#define ERROR_15 15 /*   */
#define ERROR_16 16 /*   */
#define ERROR_17 17 /*   */
#define ERROR_18 18 /*   */
#define ERROR_19 19 /*   */
#define ERROR_20 20 /*   */
#define ERROR_21 21 /*   */
#define ERROR_22 22 /*   */
#define ERROR_23 23 /*   */
#define ERROR_24 24 /*   */
#define ERROR_25 25 /*   */
#define ERROR_26 26 /*   */
#define ERROR_27 27 /*   */
#define ERROR_28 28 /*   */
#define ERROR_29 29 /*   */
#define ERROR_30 30 /*   */
#define ERROR_31 31 /*   */
#define ERROR_32 32 /*   */
#define ERROR_33 33 /*   */
#define ERROR_34 34 /*   */
#define ERROR_35 35 /*   */
#define ERROR_36 36 /*   */
#define ERROR_37 37 /*   */
#define ERROR_38 38 /*   */
#define ERROR_39 39 /*   */
#define ERROR_40 40 /*   */
#define ERROR_41 41 /*   */
#define ERROR_42 42 /*   */
#define ERROR_43 43 /*   */
#define ERROR_44 44 /*   */
#define ERROR_45 45 /*   */
#define ERROR_46 46 /*   */
#define ERROR_47 47 /*   */
#define ERROR_48 48 /*   */
#define ERROR_49 49 /*   */

typedef struct _platestats {
	char series[MAX_SERIES_STRING];
	int plateNumber;
  /* These entries are from BuildMosaicList() */
  int bestFlag; /* Found in BuildMosaicList() */
  int bestMosaicNumber;
  int bestRotation;
  int bestBinning;
  int bestWCSSource;
  char bestMosaicComment[MAX_COMMENT_STRING];
  
  /* These entries are from the photplates database */
  int photNullMosaicNumberFlag;
  int photBestFlag;
  int photBestMosaicNumber;
  int photBestHighestVersionId;
  int photBestVersionId[MAX_CATALOG_NUMBER];
  int photLatestFlag;
  int photHighestVersionId;
  int photLatestMosaicNumber;
  int photLatestVersionId[MAX_CATALOG_NUMBER];
  int logFileNumber[MAX_CATALOG_NUMBER];   /* <nnn> in phot_<nnn><v>.log */
  int logFileVersion[MAX_CATALOG_NUMBER];   /* <v> in phot_<nnn><v>.log  where 0 = '', 1 = 'a' ... */
  /* End of photplates entry */

  int resortListVersionId[MAX_CATALOG_NUMBER]; /* version ID read from the resort_magfiles plate lists */

  int versionId[MAX_CATALOG_NUMBER];    /* This versionId comes from the phot* logs  and eventually from the photinsert table */
  int mosaicNumber[MAX_CATALOG_NUMBER]; /* This mosaicNumber comes from the phot* logs and eventually from the photinsert table */
  int rotation[MAX_CATALOG_NUMBER];
  int binning[MAX_CATALOG_NUMBER];
  int plateJpegFlag;
  int jacketJpegFlag;
  int totallist_mosaicNumber;
  int totallist_rotation;
  int totallist_binning;

  int photFailFlag[MAX_CATALOG_NUMBER]; /* This flag comes from the phot* logs and eventually from the photinsert table */
  time_t timestamp[MAX_CATALOG_NUMBER];
  char plateName[MAX_PLATE_NAME];
  char jacketComments[MAX_BUFFER]; /* Also generic comments for  outputMergeFlag = 1 and  pParseCommon->jpegissueForm = 1 */
  char plateComments[MAX_BUFFER]; 
} PLATESTATS,*PPLATESTATS;

typedef struct _photversions {
  int versionId;
  char versionName[MAX_VERSION_NAME+1];
  char versionDate[MAX_TIMESTAMP_STRING+1];
} PHOTVERSION,*PPHOTVERSION;


typedef struct _pparsecommon {
	FILE *out_handle;
  FILE *merge_handle;
  FILE *platelist_handle;
  FILE *totallist_handle;
  int maxLineLen;
	time_t startTime;
	int verbose;
  int shortTotalForm;
  int shortNameForm;
  int jpegissueForm;
  int parseLogsOnly;
  int debugLogsOnly;
  int outputMergeFlag;
  int maxVersionId; /* Maximum versionId found in the logs */
	PPLATESTATS plateStats[MAX_SERIES];
  PPHOTVERSION photVersionTable;
  PHOTGLOBAL photGlobal;
  int maxPhotVersion;
  MYSQL my_connection;
  MYSQL *pConnection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection;
} PARSECOMMON,*PPARSECOMMON;

char *weekdays[] = {"Sun ","Mon ","Tue ","Wed ","Thu ","Fri ","Sat "};
char numWeekdays = sizeof(weekdays)/sizeof(char *);
char *monthstr[] = {"Jan ","Feb ","Mar ","Apr ","May ","Jun ","Jul ","Aug ","Sep ","Oct ","Nov ","Dec "};
char numMonthstr = sizeof(monthstr)/sizeof(char *);
/*
 *  Formats time of the form: "Fri Dec 20 13:43:57 EST 2013"
 *  Script started or Script done
 */
int ParseTime(PPARSECOMMON pParseCommon,char *inLine,time_t* timeResult,int *pVersionId,int lineCount) {
	int day;
	int month;
	char *charPtr = NULL;
	char timebuf[MAX_TIME_LENGTH+1];
	time_t tp;
	struct tm tmstruct;
	struct tm *ptr;
	char timestr[100];
	struct tm *gmptr;
	char gmtimestr[100];
	int nvals;
  int versionId;
	*timeResult = 0;
  *pVersionId = 0;
	for (day = 0; day < numWeekdays; day++) {
		charPtr = strstr(inLine,weekdays[day]);
		if (charPtr != NULL) {
			break;
		}
	}
	if (charPtr == NULL) {
		return(-1);
	}
	strncpy(timebuf,charPtr,MAX_TIME_LENGTH);
  timebuf[MAX_TIME_LENGTH] = 0;
	for (month = 0; month < numMonthstr; month++) {
		charPtr = strstr(timebuf,monthstr[month]);
		if (charPtr != NULL) {
			break;
		}
	}
	if (charPtr == NULL) {
		printf("ERROR: failed to find month in %s line %d\n",timebuf,lineCount);
		return(-1);
	}

	timebuf[MAX_TIME_LENGTH] = 0;
	memset(&tmstruct,0,sizeof(struct tm));
	if (strstr(timebuf,"EST") != NULL) {
		tmstruct.tm_isdst = 0;
	} else if (strstr(timebuf,"EDT") != NULL) {
		tmstruct.tm_isdst = 1;
	} else {
		if ((strstr(inLine,"Script started") != NULL) ||
				(strstr(inLine,"Script done") != NULL)) {
			return(-1);
		} else {
      if (pParseCommon->verbose) {
        printf("ERROR: unrecognized time type in %s line %d\n",inLine,lineCount);
      }
			return(-1);
		}
	}
	if (tmstruct.tm_isdst == 0) {
		nvals = sscanf(&timebuf[8],"%d %d:%d:%d EST %d",&tmstruct.tm_mday,&tmstruct.tm_hour,&tmstruct.tm_min,&tmstruct.tm_sec,&tmstruct.tm_year);
	} else {
		nvals = sscanf(&timebuf[8],"%d %d:%d:%d EDT %d",&tmstruct.tm_mday,&tmstruct.tm_hour,&tmstruct.tm_min,&tmstruct.tm_sec,&tmstruct.tm_year);
	}
	if (nvals != 5) {
		printf("ERROR: failed to decode numbers in %s line %d nvals %d %d %d:%d:%d yr %d\n",inLine,lineCount,nvals,tmstruct.tm_mday,tmstruct.tm_hour,tmstruct.tm_min,tmstruct.tm_sec,tmstruct.tm_year);
		return(-1);
	}
	
	tmstruct.tm_year -= 1900;
	tmstruct.tm_mon = month;
	
	tp = mktime(&tmstruct);
	ptr = localtime(&tp);
	strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
	gmptr = gmtime(&tp);
	strftime(gmtimestr,25,"%Y-%m-%dT%H-%M-%S", gmptr);
	if (day != tmstruct.tm_wday) {
		printf("ERROR: day disagreement %d %d for %s lineCount %d timestr %s\n",day,tmstruct.tm_wday,inLine,lineCount,timestr);
		return(-1);
	}
  /* The search needs to use GMT */
  versionId = GetPhotVersionId(pParseCommon,(char*)gmtimestr);
	if (pParseCommon->verbose) {
		printf("Found date %s  versionId %d timestr %s gmtimestr %s lineCount %d\n",inLine,versionId,timestr,gmtimestr,lineCount);
	}
  
  *pVersionId = versionId;
	*timeResult = tp;
	return(0);
}
void MakePlateName(char* plateName,char *series,int plateNumber,int mosaicNumber,int rotation,int binning,int WCSSource) {
  char nameString[MAX_PLATE_NAME];
  char tempString[MAX_PLATE_NAME];
  sprintf(nameString,"%s%05d_%02d",series,plateNumber,mosaicNumber);
  if ((binning & 1)  != 0) {
    strcat(nameString,"_01");
  } else {
    strcat(nameString,"_XX");
  }
  if (rotation != 0) {
    sprintf(tempString,"r%d",rotation);
    strcat(nameString,tempString);
  }
  if (WCSSource == WCSSOURCE_IMWCS) {
    strcat(nameString,"ww");
  }
  if (strlen(nameString) < MAX_PLATE_NAME) {
    strcpy(plateName,nameString);
  } else {
    strcpy(plateName,"ERR");
  }
  return;
}


PPLATESTATS GetPlateStats(PPARSECOMMON pParseCommon,char *series,int plateNumber) {
  int seriesId;
  PPLATESTATS pPlateStats;
  int plateNumberIndex;
  int catalogIndex;

#ifdef DEBUG_SERIES
  if (strcmp(series,DEBUG_SERIES) != 0) {
    printf("ERROR: line %4d GetPlateStats has incorrect series %s\n",__LINE__,series);
    exit(-1);
  }
#endif /* DEBUG_SERIES */


  seriesId = GetSeriesId(series,0);
  if (seriesId < 0) {
    printf("ERROR line %4d GetPlateStats bad series %s\n",__LINE__,series);
    return(NULL);
  }
  if (pParseCommon->plateStats[seriesId] == NULL) {
    pParseCommon->plateStats[seriesId] = (PPLATESTATS)calloc(MAX_PLATE_NUMBER,sizeof(PLATESTATS));
    if (pParseCommon->plateStats[seriesId] == NULL) {
      printf("ERROR line %d allocating platestats for series %s\n",__LINE__,series);
      return(NULL);
    }
    for (plateNumberIndex = 0; plateNumberIndex < MAX_PLATE_NUMBER; plateNumberIndex++) {
      pPlateStats = &(pParseCommon->plateStats[seriesId])[plateNumberIndex];
      pPlateStats->plateNumber = -1;
      pPlateStats->bestFlag = 0;
      pPlateStats->bestMosaicNumber = -1;
      pPlateStats->bestRotation = -1;
      pPlateStats->bestBinning = -1;
      pPlateStats->bestWCSSource = WCSSOURCE_NONE;
      pPlateStats->photBestFlag = -1;
      pPlateStats->photBestMosaicNumber = -1;
      pPlateStats->photLatestFlag = -1;
      pPlateStats->photHighestVersionId = -1;
      pPlateStats->photBestHighestVersionId = -1;
      pPlateStats->photLatestMosaicNumber = -1;
      pPlateStats->totallist_mosaicNumber = -1;
      pPlateStats->totallist_rotation = -1;
      pPlateStats->totallist_binning = -1;
                 
      for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
        pPlateStats->versionId[catalogIndex] = -1;
        pPlateStats->mosaicNumber[catalogIndex] = -1;
        pPlateStats->rotation[catalogIndex] = -1;
        pPlateStats->binning[catalogIndex] = -1;
        pPlateStats->photBestVersionId[catalogIndex] = -1;
        pPlateStats->photLatestVersionId[catalogIndex] = -1;
        pPlateStats->resortListVersionId[catalogIndex] = -1;
        pPlateStats->logFileNumber[catalogIndex] = -1;
        pPlateStats->logFileVersion[catalogIndex] = -1;
      }
    }
  }
            
  pPlateStats = &(pParseCommon->plateStats[seriesId])[plateNumber];

  return(pPlateStats);

}


void InitVersionTable(PPARSECOMMON pParseCommon) {
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  PPHOTVERSION pPhotVersion;
  int nvals;
  unsigned long numRows;
  unsigned long curRow = 0;
  char queryString[MAX_QUERY_STRING];
  if ((pParseCommon->photVersionTable != NULL) ||
      (pParseCommon->maxPhotVersion != 0)) {
    printf("ERROR: InitVersionTable line %d version table already initialized\n",__LINE__);
    exit(-1);
  }
  sprintf(queryString,"SELECT versionId,versionName,versionDate from photversions order by versionDate\n;");
  res = ExecuteQuery(pParseCommon->pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pParseCommon->pPhotConnection);
    numRows = mysql_num_rows(res_ptr);
    pParseCommon->photVersionTable = (PPHOTVERSION)calloc(numRows,sizeof(PHOTVERSION));
    if (pParseCommon->photVersionTable == NULL) {
      printf("ERROR: InitVersionTable line %d failed to allocate version table of size %d\n",__LINE__,numRows);
      exit(-1);
      
    }
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if ((curRow >= numRows) || (curRow < 0)) {
          printf("ERROR: InitVersionTable line %d exceeded allocated size %d\n",__LINE__,numRows);
          exit(-1);
        }
        pPhotVersion = &pParseCommon->photVersionTable[curRow];
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pPhotVersion->versionId);
          if (nvals != 1) {
            printf("ERROR: InitVersionTable line %d failed to parse versionId\n",__LINE__);
            exit(-1);
          }
        } else {
          printf("ERROR: InitVersionTable line %d missing versionId\n",__LINE__);
          exit(-1);    
        }

        if (sqlrow[1]) {
          strncpy(pPhotVersion->versionName,sqlrow[1],MAX_VERSION_NAME);
          pPhotVersion->versionName[MAX_VERSION_NAME] = 0;
        } else {
          printf("ERROR: InitVersionTable line %d missing versionName\n",__LINE__);
          exit(-1);    
        }
        if (sqlrow[2]) {
          strncpy(pPhotVersion->versionDate,sqlrow[2],MAX_TIMESTAMP_STRING);
          pPhotVersion->versionDate[MAX_TIMESTAMP_STRING] = 0;
        } else {
          printf("ERROR: InitVersionTable line %d missing versionDate\n",__LINE__);
          exit(-1);    
        }
        curRow++;
        if (curRow != pPhotVersion->versionId) {
          printf("ERROR: InitVersionTable line %d row %d versionId %d mismatch\n",__LINE__);
          exit(-1);    
        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }
  pParseCommon->maxPhotVersion = curRow;


  return;
} /* End of InitVersionTable() */
/* NOTE: timestr must be in GMT */
int GetPhotVersionId(PPARSECOMMON pParseCommon,char* timestr) {
  int versionIdIndex;
  PPHOTVERSION pPhotVersion;
  int result;
  int lastVersionId = -1;

  for (versionIdIndex = 0; versionIdIndex < pParseCommon->maxPhotVersion; versionIdIndex++) {
    pPhotVersion = &pParseCommon->photVersionTable[versionIdIndex];
    result = strcmp(timestr,pPhotVersion->versionDate);
    lastVersionId = pPhotVersion->versionId;
#if 0
    printf("line %4d result %2d timestr %s versionIdIndex %d versionId %d versionName %s versionDate %s\n",
           __LINE__,
           result,
           timestr,
           versionIdIndex,
           pPhotVersion->versionId,
           pPhotVersion->versionName,
           pPhotVersion->versionDate);
#endif
    if (result < 0) {
      break;
    }
  }
  return(lastVersionId);
}
void ReadPhotPlatesTable(PPARSECOMMON pParseCommon) {
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  PPHOTVERSION pPhotVersion;
  int nvals;
  unsigned long numRows;
  unsigned long curRow = 0;
  int seriesId;
  int plateNumber;
  char series[MAX_SERIES_STRING];
  char queryString[MAX_QUERY_STRING];
  PPLATESTATS pPlateStats;
  int versionId[MAX_CATALOG_NUMBER+1];
  int mosaicNumber;
  int catalogNumber;
  int highestVersionId;


#ifdef DEBUG_SERIESID
  sprintf(queryString,"SELECT seriesId,plateNumber,mosaicNumber,versionID,versionID0,versionID1,versionID2,versionID3,versionID4 from photplates where seriesId = %d\n",DEBUG_SERIESID);
#else
  /*                 catalogNumber                                 GSC 0         1      KIC 2    APASS 3     GAIA 4    ATLAS 5                   */
  sprintf(queryString,"SELECT seriesId,plateNumber,mosaicNumber,versionID,versionID0,versionID1,versionID2,versionID3,versionID4 from photplates\n;");
#endif
  res = ExecuteQuery(pParseCommon->pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pParseCommon->pPhotConnection);
    numRows = mysql_num_rows(res_ptr);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&seriesId);
          if (nvals != 1) {
            printf("ERROR: ReadPhotPlatesTable line %d failed to parse seriesId\n",__LINE__);
            exit(-1);
          }
        } else {
          printf("ERROR: ReadPhotPlatesTable line %d missing seriesId\n",__LINE__);
          exit(-1);    
        }
        if (GetSequesteredFlag(seriesId) == SEQUESTERED_YES) {
#if 0
          printf("line %d seriesId %d id sequestered\n",__LINE__,seriesId);
#endif
          continue;
        }
        strcpy(series,GetSeriesString(seriesId,1));



        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%d",&plateNumber);
          if (nvals != 1) {
            printf("ERROR: ReadPhotPlatesTable line %d failed to parse plateNumber\n",__LINE__);
            exit(-1);
          }
        } else {
          printf("ERROR: ReadPhotPlatesTable line %d missing plateNumber\n",__LINE__);
          exit(-1);    
        }
        pPlateStats = GetPlateStats(pParseCommon,series,plateNumber);
        if (pPlateStats == NULL) {
          printf("ERROR line %d GetPLateStats failed for mosaic %s%05d\n",__LINE__,series,plateNumber);
          exit(-1);
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&mosaicNumber);
          if (nvals != 1) {
            printf("ERROR: ReadPhotPlatesTable line %d failed to parse mosaicNumber\n",__LINE__);
            exit(-1);
          }
          int highestVersionId  = -1;
          memset(versionId,0,sizeof(versionId));
          for (catalogNumber = 0; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
            
            if (sqlrow[3+catalogNumber]) {
              nvals = sscanf(sqlrow[3+catalogNumber],"%d",&versionId[catalogNumber]);
              if (nvals != 1) {
                printf("ERROR: ReadPhotPlatesTable line %d failed to parse versionId%d\n",__LINE__,catalogNumber);
                exit(-1);
              }
              if (versionId[catalogNumber] > highestVersionId) {
                highestVersionId = versionId[catalogNumber];
              }
            } else {
              versionId[catalogNumber] = -1;
            }
          }
          if (versionId[1] >= 0) {
            printf("ERROR: ReadPhotPlatesTable line %d non-null versionId0\n",__LINE__);
            exit(-1);
          }
          if (mosaicNumber == pPlateStats->bestMosaicNumber) {
            pPlateStats->photBestFlag = 1;
            pPlateStats->photBestMosaicNumber = mosaicNumber;
            pPlateStats->photBestHighestVersionId = highestVersionId;
            for (catalogNumber = 0; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
              if (catalogNumber == 0) {
                pPlateStats->photBestVersionId[0] = versionId[0];
              } else {
                pPlateStats->photBestVersionId[catalogNumber] = versionId[catalogNumber+1];
              }
            }
          }
          if (highestVersionId > pPlateStats->photHighestVersionId) {
            pPlateStats->photLatestFlag = 1;
            pPlateStats-> photLatestMosaicNumber = mosaicNumber;
            pPlateStats->photHighestVersionId = highestVersionId;
            for (catalogNumber = 0; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
              if (catalogNumber == 0) {
                pPlateStats->photLatestVersionId[0] = versionId[0];
              } else {
                pPlateStats->photLatestVersionId[catalogNumber] = versionId[catalogNumber+1];
              }
            }
          }
        } else {
          pPlateStats->photNullMosaicNumberFlag = 1;
        }
        curRow++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }
  printf("ReadPhotPlatesTable found %d entries in line %d\n",curRow,__LINE__);

  return;
} /* End of ReadPhotPlatesTable() */


int ProcessTotalList(PPARSECOMMON pParseCommon) {
  char *inBuffer;
  char inLine[MAX_BUFFER];
  int lineLen;
  int lineCount = 0;
  char copyLine[MAX_BUFFER];
  char plateName[MAX_PLATE_NAME];
  char comments[MAX_BUFFER];  /* Contains everything following the first token */
  char series[MAX_SERIES_STRING];
  int seriesId;
  int plateNumber;
  int mosaicNumber;
  int rotation;
  int binning;
  PPLATESTATS pPlateStats;
  int errorCount = 0;
  char * charPtr;
  int jpegType;
  int tokenLength;
  char *spacePtr;
  int commentsLength;
  int prevLength;
  int plateNameLength;
  int curRow = 0;

  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,pParseCommon->totallist_handle);
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
    strcpy(copyLine,inLine);
    lineCount++;
    jpegType = JPEG_TYPE_NONE;
    spacePtr = strchr(inLine,' ');
    if (spacePtr == NULL) {
      comments[0] = 0;
    } else {
      *spacePtr = 0;
      *spacePtr++;
      strcpy(comments,spacePtr);
    }
    plateNameLength = strlen(inLine);
    if (plateNameLength > (MAX_PLATE_NAME-1)) {
      printf("ERROR: line %d MAX_PLATE_NAME exceeded with %d in %s\n",__LINE__,plateNameLength,copyLine);
      exit(-1);
    }
    strcpy(plateName,inLine);
    if (pParseCommon->jpegissueForm == 0) {
      if (ParseFilename(inLine,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
        if (pParseCommon->shortNameForm == 1) {
          mosaicNumber = 0;
          binning = 0;
          rotation = 0;
          if (ParseFilename2(inLine,series,&plateNumber) == 0) {
            printf("ERROR line %d ProcessTotalList failed to parse lineCount %d %s\n",__LINE__,lineCount,copyLine);
            exit(-1);
          }

        } else {
          if (pParseCommon->shortTotalForm != 0) {
            printf("ERROR line %d ProcessTotalList failed to parse lineCount %d %s\n",__LINE__,lineCount,copyLine);
            exit(-1);
          } else {
            strcat(inLine,"_01");
            if (ParseFilename(inLine,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
              printf("ERROR line %d ProcessTotalList failed to parse lineCount %d %s\n",__LINE__,lineCount,copyLine);
              exit(-1);
            }
          }
        }
      }
    } else {

      charPtr = strstr(inLine,"_j");
      if (charPtr == NULL) {
        charPtr = strstr(inLine,"_p");
        if (charPtr == NULL) {
          printf("ERROR line %d ProcessTotalList failed find '_p' nor '_j' in lineCount %d %s\n",__LINE__,lineCount,copyLine);
          exit(-1);
        } else {
          jpegType = JPEG_TYPE_PLATE;
        }
      } else {
        jpegType = JPEG_TYPE_JACKET;
      }
      *charPtr = 0;
      mosaicNumber = 0;
      binning = 0;
      rotation = 0;
      if (ParseFilename2(inLine,series,&plateNumber) == 0) {
        printf("ERROR line %d ProcessTotalList failed to parse lineCount %d %s\n",__LINE__,lineCount,copyLine);
        exit(-1);
      }
    }
#if 0
    printf("line %d ProcessTotalList processing %s\n",__LINE__,inLine);
#endif
#ifdef DEBUG_SERIES
    if (strcmp(series,DEBUG_SERIES) != 0) {
      continue;
    }
#endif /* DEBUG_SERIES */
  

   pPlateStats = GetPlateStats(pParseCommon,series,plateNumber);
    if (pPlateStats == NULL) {
      printf("ERROR line %d GetPLateStats failed with lineCount %d %s\n",__LINE__,lineCount,copyLine);
      exit(-1);
    }
    if (pPlateStats->plateNumber != plateNumber) {
      strcpy(pPlateStats->series,series);
      pPlateStats->plateNumber = plateNumber;
    }
    if (jpegType == JPEG_TYPE_JACKET) {
      if (pPlateStats->jacketJpegFlag != 0) {
        printf("ERROR line %d duplicate plate jacket lineCount %d %s\n",__LINE__,lineCount,copyLine); 
        errorCount++;
      }
      pPlateStats->jacketJpegFlag = 1;
    }
    if (jpegType == JPEG_TYPE_PLATE) {
      if (pPlateStats->plateJpegFlag != 0) {
        printf("ERROR line %d duplicate plate jacket lineCount %d %s\n",__LINE__,lineCount,copyLine); 
        errorCount++;
      }
      pPlateStats->plateJpegFlag = 1;
    }
    /* Copy in the plateName if it is longer */
    if (pPlateStats->plateName[0] == 0) {
      strcpy(pPlateStats->plateName,plateName);
    } else {
      if (strlen(pPlateStats->plateName) < plateNameLength) {
        strcpy(pPlateStats->plateName,plateName);
      }
    }
    /* Here we append the comments to the previously gathered comments */
    
    if ((pParseCommon->outputMergeFlag == 1) ||
        (pParseCommon->jpegissueForm == 0) ||  
        (jpegType == JPEG_TYPE_JACKET))  {
      if (pPlateStats->jacketComments[0] == 0) {
        strcpy(pPlateStats->jacketComments,comments);
      } else {
        if (comments[0] != 0) {
          commentsLength = strlen(comments);
  
          prevLength = strlen(pPlateStats->jacketComments);
          if ((commentsLength+prevLength+2) >= MAX_BUFFER) {
            printf("ERROR: line %d comments length exceeded by %s\n",__LINE__,copyLine);
          } else {
            strcat(comments," ");
            strcat(comments,pPlateStats->jacketComments);
            strcpy(pPlateStats->jacketComments,comments);
          }
        }
      }
    } else {
      if (pPlateStats->plateComments[0] == 0) {
        strcpy(pPlateStats->plateComments,comments);
      } else {
        if (comments[0] != 0) {
          commentsLength = strlen(comments);
  
          prevLength = strlen(pPlateStats->plateComments);
          if ((commentsLength+prevLength+2) >= MAX_BUFFER) {
            printf("ERROR: line %d comments length exceeded by %s\n",__LINE__,copyLine);
          } else {
            strcat(comments," ");
            strcat(comments,pPlateStats->plateComments);
            strcpy(pPlateStats->plateComments,comments);
          }
        }
      }
    }
    if ((pPlateStats->totallist_mosaicNumber > 0) ||
        (pPlateStats->totallist_rotation > 0) ||
        (pPlateStats->totallist_binning > 0)) {
      printf("ERROR line %d ProcessTotalList found duplicate plate %s%05d\n",__LINE__,series,plateNumber);
      errorCount++;
    } else {
      pPlateStats->totallist_mosaicNumber = mosaicNumber;
      pPlateStats->totallist_rotation = rotation;
      pPlateStats->totallist_binning = binning;
    }
    curRow++;
  }
  printf("ProcessTotalList processed %d of %d plates errors: %d unique plates %d\n",curRow,lineCount,errorCount,lineCount-errorCount);


} /* End of ProcessTotalList() */

int PrintPlateList(PPARSECOMMON pParseCommon) {
  int seriesIdIndex;
  int plateNumberIndex;
  int catalogIndex;
  int mosaicNumber = -1;
  int bestCatalogIndex = -1;
  char curPlate[2*MAX_PLATE_NAME];
  PPLATESTATS pPlateStats;
  char* photFailFlagString[2] = {"","X"};
  char* reprocessString[2] = {"",",REPROCESS"};
  char totallistPlate[2*MAX_PLATE_NAME];
  char agreeString[2*MAX_PLATE_NAME];
  int reprocessFlag;
  char *charPtr;
  char bestPlateName[MAX_PLATE_NAME];
  char logPlateName[MAX_PLATE_NAME];
  int outputFlag[MAX_CATALOG_NUMBER];
  char prefix[MAX_BUFFER];
  char testprefix[MAX_BUFFER];
  int errorFlag = 0;
  int testOutputFlag = 0;
  int seriesId = -1;
  PPHOTGLOBAL pPhotGlobal = &pParseCommon->photGlobal;

  for (seriesIdIndex = 0; seriesIdIndex < MAX_SERIES; seriesIdIndex++) {
    if (pParseCommon->plateStats[seriesIdIndex] == NULL) {
      continue;
    }    
    if (GetSequesteredFlag(seriesIdIndex) == SEQUESTERED_YES) {
#if 0
      printf("line %d seriesId %d id sequestered\n",__LINE__,seriesIdIndex);
#endif
      continue;
    }
    for (plateNumberIndex = 0; plateNumberIndex < MAX_PLATE_NUMBER; plateNumberIndex++) {            
      pPlateStats = &(pParseCommon->plateStats[seriesIdIndex])[plateNumberIndex];
      if (pPlateStats->plateNumber > 0) {
        if (pParseCommon->jpegissueForm != 0) {
          if (pPlateStats->jacketJpegFlag != 0) {
            if (pPlateStats->jacketComments[0] != 0) {
              fprintf(pParseCommon->out_handle,"%s%05d_j %s\n",pPlateStats->series,pPlateStats->plateNumber,pPlateStats->jacketComments);
            } else {
              fprintf(pParseCommon->out_handle,"%s%05d_j\n",pPlateStats->series,pPlateStats->plateNumber);
            }
          }
          if (pPlateStats->plateJpegFlag != 0) {
            if (pPlateStats->plateComments[0] != 0) {
              fprintf(pParseCommon->out_handle,"%s%05d_p %s\n",pPlateStats->series,pPlateStats->plateNumber,pPlateStats->plateComments);
            } else {
              fprintf(pParseCommon->out_handle,"%s%05d_p\n",pPlateStats->series,pPlateStats->plateNumber);
            }
          }
          continue;
        }
        if (pParseCommon->shortTotalForm == 1) {
          continue;
        }
          
        if (pParseCommon->outputMergeFlag) {
          if (pPlateStats->jacketComments[0] != 0) {
            fprintf(pParseCommon->merge_handle,"%s %s\n",pPlateStats->plateName,pPlateStats->jacketComments);
          } else {
            fprintf(pParseCommon->merge_handle,"%s\n",pPlateStats->plateName);
          }

          continue;
        }

        /* New processing */
        errorFlag = 0;

        if (pParseCommon->debugLogsOnly == 1) { /* ifdef debugLogsOnly */
          for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
            if (pPlateStats->versionId[catalogIndex] > 0) {
              if ((pPlateStats->logFileNumber[catalogIndex] >= 0) &&
                  (pPlateStats->logFileVersion[catalogIndex] >= 0)) {
#if 1
                seriesId = GetSeriesId(pPlateStats->series,0);
                printf("INSERT IGNORE INTO  photinsert (seriesId,plateNumber,catalogNumber) VALUES (%d,%d,%d);\n",
                       seriesId,
                       pPlateStats->plateNumber,
                       catalogIndex);
                printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                       pPlateStats->mosaicNumber[catalogIndex],
                       pPlateStats->versionId[catalogIndex],
                       pPlateStats->photFailFlag[catalogIndex],
                       seriesId,
                       pPlateStats->plateNumber,
                       catalogIndex);
#endif
              }
            }
          }
        } /* ifdef debugLogsOnly */
        if (pPlateStats->bestFlag != 0) {
          MakePlateName(&bestPlateName[0],
                        pPlateStats->series,
                        pPlateStats->plateNumber,
                        pPlateStats->bestMosaicNumber,
                        pPlateStats->bestRotation,
                        pPlateStats->bestBinning,
                        pPlateStats->bestWCSSource);
        } else {
          memset(bestPlateName,0,sizeof(bestPlateName));
        }
        sprintf(prefix,"\nNEW %s FULL AGREEMENT",bestPlateName);
        memset(outputFlag,0,sizeof(outputFlag));
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }

#ifdef DEBUG_PHOTINSERT
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if (outputFlag[catalogIndex] == 0) {
#if 1
            if ((strcmp(pPlateStats->series,"ai") == 0) &&
                (pPlateStats->plateNumber == 39438)) {
              printf("\nline %4d catalogNumber %d DEBUG %s for %s%05d\n",__LINE__,catalogIndex,bestPlateName,pPlateStats->series,pPlateStats->plateNumber);
            }
#endif

            printf("\nline %4d %s catalogIndex %d %d %d %d %d %d %d %d %d %d %d\n",
                   __LINE__,
                   bestPlateName,
                   catalogIndex,
                   (pPlateStats->resortListVersionId[catalogIndex] > 0),
                   (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]),
                   (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]),
                   (pPlateStats->bestFlag == 1),
                   (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS),
                   (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]),
                   (pPlateStats->photFailFlag[catalogIndex] == 0),
                   (strcmp(pPlateStats->plateName,bestPlateName) == 0),
                   (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") != NULL),
                   (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)));
            if (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)) {
              printf("\nline %4d ERROR DEBUG bestPlateName %s\n",
                     __LINE__,
                     bestPlateName);
            }
          }
        }
#endif /* DEBUG_PHOTINSERT */


        sprintf(prefix,"\nNEW %s ERROR%02d Fix PHOTINSERT (FixJPEG deletion)",bestPlateName,ERROR_01);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex] -1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              ((pPlateStats->bestMosaicNumber) == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->resortListVersionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);
            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE WCSSOURCE_IMWCS Rejected for Markings",bestPlateName,ERROR_02);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] == -1) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->mosaicNumber[catalogIndex] == -1) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") != NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed (FixJPEG deletion)",bestPlateName,ERROR_03);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex] -1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              ((pPlateStats->bestMosaicNumber) != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;            
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed (FixJPEG deletion)",bestPlateName,ERROR_04);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] < pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("INSERT IGNORE INTO  photinsert (seriesId,plateNumber,catalogNumber) VALUES (%d,%d,%d);\n",
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }




        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed (0 entries)",bestPlateName,ERROR_05);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->bestMosaicNumber,
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }



        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed",bestPlateName,ERROR_06);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex]+1,
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d RETRY PURGE Fix PHOTINSERT purge needed retry with best mosaic",bestPlateName,ERROR_07);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex]+1,
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d RETRY PURGE Fix PHOTINSERT purge needed retry with best mosaic",bestPlateName,ERROR_08);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d REPROCESS (check run and phot logs  atlas_completed.list)",bestPlateName,ERROR_09);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */

          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] > pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {

            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d REPROCESS (check run and phot logs atlas_nofirstdata.list)",bestPlateName,ERROR_10);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->mosaicNumber[catalogIndex] == -1) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE ",bestPlateName,ERROR_11);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE WCSSOURCE_LOGBOOK (FixJPEG deletion)",bestPlateName,ERROR_12);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] == -1) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->mosaicNumber[catalogIndex] == -1) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0]  == 0)) {

            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE WCSSOURCE_LOGBOOK check photBestVersionId",bestPlateName,ERROR_13);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] <= pPlateStats->photLatestVersionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE WCSSOURCE_LOGBOOK for new mosaic",bestPlateName,ERROR_14);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE WCSSOURCE_LOGBOOK for new mosaic (0 entries)",bestPlateName,ERROR_15);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE WCSSOURCE_LOGBOOK for mosaic",bestPlateName,ERROR_16);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0]  == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d RETRY rebuild reinsert old mosaic",bestPlateName,ERROR_17);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */

          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d RETRY rebuild reinsert",bestPlateName,ERROR_18);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s %s%05d ERROR%02d IGNORE",bestPlateName,pPlateStats->series,pPlateStats->plateNumber,ERROR_19);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag != 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }



        sprintf(prefix,"\nNEW %s %s%05d ERROR%02d PURGE Fix PHOTINSERT purge needed",bestPlateName,pPlateStats->series,pPlateStats->plateNumber,ERROR_20);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)) &&
              (pPlateStats->bestFlag != 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }




        sprintf(prefix,"\nNEW %s %s%05d ERROR%02d PURGE Fix PHOTINSERT purge needed (0 entries)",bestPlateName,pPlateStats->series,pPlateStats->plateNumber,ERROR_21);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)) &&
              (pPlateStats->bestFlag != 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->mosaicNumber[catalogIndex]  >= 0) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }

        sprintf(prefix,"\nNEW %s %s%05d ERROR%02d IGNORE ",bestPlateName,pPlateStats->series,pPlateStats->plateNumber,ERROR_22);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] == -1) &&
              (pPlateStats->bestFlag != 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->mosaicNumber[catalogIndex] == -1) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0))  {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }
        }

        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE ",bestPlateName,ERROR_23);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }
        }

        sprintf(prefix,"\nNEW %s ERROR%02d Fix PHOTINSERT, no purge",bestPlateName,ERROR_24);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->versionId[catalogIndex] == (pPlateStats->photBestVersionId[catalogIndex]+1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("INSERT IGNORE INTO  photinsert (seriesId,plateNumber,catalogNumber) VALUES (%d,%d,%d);\n",
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed (0 entries)",bestPlateName,ERROR_25);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] <  pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") != NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE WCSSOURCE_IMWCS Rejected for Markings",bestPlateName,ERROR_26);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] == -1) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->mosaicNumber[catalogIndex] == -1) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") != NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed",bestPlateName,ERROR_27);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] < pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d Fix PHOTINSERT no purge",bestPlateName,ERROR_28);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->versionId[catalogIndex] == -1) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->mosaicNumber[catalogIndex] == -1) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            printf("INSERT IGNORE INTO  photinsert (seriesId,plateNumber,catalogNumber) VALUES (%d,%d,%d);\n",
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->bestMosaicNumber,
                   pPlateStats->resortListVersionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed",bestPlateName,ERROR_29);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] < pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }



        sprintf(prefix,"\nNEW %s %s%05d ERROR%02d PURGE Fix PHOTINSERT purge needed (FixJPEG deletion)",bestPlateName,pPlateStats->series,pPlateStats->plateNumber,ERROR_30);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag != 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->bestMosaicNumber == -1) &&
              (pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex]+1,
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d PURGE purge needed (0 entries)",bestPlateName,ERROR_31);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] < pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d Fix PHOTINSERT no purge",bestPlateName,ERROR_32);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] == 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE ",bestPlateName,ERROR_33);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] == -1) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->bestMosaicNumber >= 0) && 
              (pPlateStats->mosaicNumber[catalogIndex] == -1) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }
        }

        sprintf(prefix,"\nNEW %s ERROR%02d RETRY PURGE Fix PHOTINSERT purge needed retry reinsert",bestPlateName,ERROR_34);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;            
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex]+1,
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d RETRY build and insert ",bestPlateName,ERROR_35);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] < pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }
        }

        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed",bestPlateName,ERROR_36);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex] -1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              ((pPlateStats->bestMosaicNumber) == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;            
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);
            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d RETRY PURGE Fix PHOTINSERT purge needed retry with best mosaic (0 entries)",bestPlateName,ERROR_37);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) != 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPhotGlobal->minVersionId[catalogIndex],
                   pPlateStats->versionId[catalogIndex]+1,
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE (FixJPEG deletion)",bestPlateName,ERROR_38);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) != 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }
        }

        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed (0 entries)",bestPlateName,ERROR_39);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex] -1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              ((pPlateStats->bestMosaicNumber) == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) != 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;            
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);
            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE (plate d11023)",bestPlateName,ERROR_40);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] == -1) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->mosaicNumber[catalogIndex] == -1) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }
        }

        sprintf(prefix,"\nNEW %s %s%05d ERROR%02d IGNORE (FixJPEG deletion)",bestPlateName,pPlateStats->series,pPlateStats->plateNumber,ERROR_41);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] == -1) &&
              (pPlateStats->versionId[catalogIndex] > 0) &&
              (pPlateStats->bestFlag != 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0] != 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d RETRY PURGE Fix PHOTINSERT purge needed retry with best mosaic (0 entries)",bestPlateName,ERROR_42);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

          }          
        }

        sprintf(prefix,"\nNEW %s %s%05d ERROR%02d PURGE Fix PHOTINSERT purge needed",bestPlateName,pPlateStats->series,pPlateStats->plateNumber,ERROR_43);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex] -1)) &&
              (pPlateStats->bestFlag != 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] != 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed (FixJPEG deletion)",bestPlateName,ERROR_44);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex] -1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (strcmp(pPlateStats->plateName,bestPlateName) != 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            pPlateStats->photFailFlag[catalogIndex] = 1;
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d Fix PHOTINSERT no purge (FixJPEG deletion)",bestPlateName,ERROR_45);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex] -1)) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->resortListVersionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed (FixJPEG deletion)",bestPlateName,ERROR_46);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] < pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag != 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_NONE) &&
              (pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) != 0)) {
            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }

        sprintf(prefix,"\nNEW %s ERROR%02d PURGE Fix PHOTINSERT purge needed (FixJPEG deletion)",bestPlateName,ERROR_47);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] < pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_LOGBOOK) &&
              (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 1) &&
              (strcmp(pPlateStats->plateName,bestPlateName) != 0)) {

            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }


        sprintf(prefix,"\nNEW %s ERROR%02d IGNORE",bestPlateName,ERROR_48);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if ((outputFlag[catalogIndex] == 0) && 
              (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
              (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
              (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
              (pPlateStats->bestFlag == 1) &&
              (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
              (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
              (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
              (pPlateStats->photFailFlag[catalogIndex] == 0) &&
              (pPlateStats->plateName[0] == 0)) {

            outputFlag[catalogIndex] = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            seriesId = GetSeriesId(pPlateStats->series,0);
            printf("UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                   pPlateStats->mosaicNumber[catalogIndex],
                   pPlateStats->versionId[catalogIndex],
                   pPlateStats->photFailFlag[catalogIndex],
                   seriesId,
                   pPlateStats->plateNumber,
                   catalogIndex);

            prefix[0] = 0;
          }          
        }


#if 0
        if ((outputFlag[catalogIndex] == 0) && 
            (pPlateStats->resortListVersionId[catalogIndex] > 0) &&
            (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]) &&
            (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]) &&
            (pPlateStats->bestFlag == 1) &&
            (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS) &&
            (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]) &&
            (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL) &&
            (pPlateStats->photFailFlag[catalogIndex] == 0) &&
            (strcmp(pPlateStats->plateName,bestPlateName) == 0)) {

        }
#endif

        testOutputFlag = 0;
        sprintf(prefix,"\nNEW %s %s%05d ERROR%02d NOT HANDLED for new mosaic (0 entries)",bestPlateName,pPlateStats->series,pPlateStats->plateNumber,ERROR_49);
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
#if 1
          if (outputFlag[catalogIndex] == 0) { 
            printf("\n%d %d %d %d %d %d %d %d %d %d",
                   (pPlateStats->resortListVersionId[catalogIndex] > 0),
                   (pPlateStats->resortListVersionId[catalogIndex] >= pPhotGlobal->minVersionId[catalogIndex]),
                   (pPlateStats->resortListVersionId[catalogIndex] == pPlateStats->versionId[catalogIndex]),
                   (pPlateStats->bestFlag == 1),
                   (pPlateStats->bestWCSSource ==  WCSSOURCE_IMWCS),
                   (pPlateStats->bestMosaicNumber == pPlateStats->mosaicNumber[catalogIndex]),
                   (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") == NULL),
                   (pPlateStats->photFailFlag[catalogIndex] == 0),
                   (strcmp(pPlateStats->plateName,bestPlateName) == 0),
                   (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)));
          }
#endif
          if (outputFlag[catalogIndex] == 0) {
            strcpy(testprefix,prefix);
            if (pPlateStats->resortListVersionId[catalogIndex] <= 0) {
              printf("%s TEST 1 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }
            if (pPlateStats->resortListVersionId[catalogIndex] < pPhotGlobal->minVersionId[catalogIndex])  {
              printf("%s TEST 2 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }
            if (pPlateStats->resortListVersionId[catalogIndex] != pPlateStats->versionId[catalogIndex])  {
              printf("%s TEST 3 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }
            if (pPlateStats->bestFlag != 1)  {
              printf("%s TEST 4 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }
            if (pPlateStats->bestWCSSource !=  WCSSOURCE_IMWCS)  {
              printf("%s TEST 5 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }
            if (pPlateStats->bestMosaicNumber != pPlateStats->mosaicNumber[catalogIndex])  {
              printf("%s TEST 6 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }
            if (pPlateStats->photFailFlag[catalogIndex] != 0)  {
              printf("%s TEST 7 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }
            if (strcmp(pPlateStats->plateName,bestPlateName) != 0)  {
              printf("%s TEST 8 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }
            if (strstr(pPlateStats->bestMosaicComment,"Rejected for Markings") != NULL) {
              printf("%s TEST 9 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }

            if (pPlateStats->resortListVersionId[catalogIndex] == (pPlateStats->versionId[catalogIndex]-1)) {
              printf("%s TEST 10 catalogIndex %d ",testprefix,catalogIndex);
              testprefix[0] = 0;
            }

            outputFlag[catalogIndex] = 1;
            testOutputFlag = 1;
            fprintf(pParseCommon->platelist_handle,"%s %s",prefix,catalogText[catalogIndex]);
            prefix[0] = 0;
          }          
        }
        if (testOutputFlag == 1) {
          printf("\n");
          testOutputFlag = 0;
        }

        fprintf(pParseCommon->platelist_handle,"\n");

        totallistPlate[0] = 0;
        strcpy(agreeString,"UNK");
        if (pPlateStats->totallist_mosaicNumber >= 0) {
          if (pPlateStats->rotation[bestCatalogIndex] == 0) {
            sprintf(totallistPlate,"%s%05d_%02d_%02dww",
                    pPlateStats->series,
                    pPlateStats->plateNumber,
                    pPlateStats->totallist_mosaicNumber,
                    pPlateStats->totallist_binning);
          } else {
            sprintf(totallistPlate,"%s%05d_%02d_%02dr%dww",
                    pPlateStats->series,
                    pPlateStats->plateNumber,
                    pPlateStats->totallist_mosaicNumber,
                    pPlateStats->totallist_binning,
                    pPlateStats->totallist_rotation);
          }
          

        } else {
          strcpy(totallistPlate,"NOTOTAL");
        }


        /* We want the catalog with the highest mosaic number */
        mosaicNumber = -1;
        bestCatalogIndex = -1;
        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          if (pPlateStats->mosaicNumber[catalogIndex] >  mosaicNumber) {
            mosaicNumber = pPlateStats->mosaicNumber[catalogIndex];
            bestCatalogIndex = catalogIndex;
          }
        }
        if (bestCatalogIndex >= 0) {
          if (pPlateStats->rotation[bestCatalogIndex] == 0) {
            sprintf(curPlate,"%s%05d_%02d_%02dww",
                    pPlateStats->series,
                    pPlateStats->plateNumber,
                    pPlateStats->mosaicNumber[bestCatalogIndex],
                    pPlateStats->binning[bestCatalogIndex]);
          } else {
            sprintf(curPlate,"%s%05d_%02d_%02dr%dww",
                    pPlateStats->series,
                    pPlateStats->plateNumber,
                    pPlateStats->mosaicNumber[bestCatalogIndex],
                    pPlateStats->binning[bestCatalogIndex],
                    pPlateStats->rotation[bestCatalogIndex]);
          }
          if (pPlateStats->totallist_mosaicNumber >= 0) {
            if (pPlateStats->totallist_mosaicNumber ==  pPlateStats->mosaicNumber[bestCatalogIndex]) {
              strcpy(agreeString,"AGREE");
            } else {
              strcpy(agreeString,"DIFFER");
            }
          }
          fprintf(pParseCommon->platelist_handle,"OLD %s %s %s",curPlate,totallistPlate,agreeString);

          for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
            if ((catalogIndex == CATALOG_GAIA) ||
                (catalogIndex == CATALOG_KEPLER) ||
                (catalogIndex == CATALOG_EXPERIMENTAL)) {
              continue;
            }
#endif /* SKIP_GAIA */
            if (pPlateStats->versionId[catalogIndex] >= 0) {
              reprocessFlag = 0;
              if ((pPlateStats->totallist_mosaicNumber >= 0) && 
                  (pPlateStats->mosaicNumber[catalogIndex] > pPlateStats->totallist_mosaicNumber)) {
                reprocessFlag = 1;
              }

              fprintf(pParseCommon->platelist_handle," %s:%02d:%3d%s%s",
                      catalogText[catalogIndex],
                      pPlateStats->mosaicNumber[catalogIndex],
                      pPlateStats->versionId[catalogIndex],
                      photFailFlagString[pPlateStats->photFailFlag[catalogIndex]],
                      reprocessString[reprocessFlag]);
            }
          }
          fprintf(pParseCommon->platelist_handle,"\n");
        
        } else {
          fprintf(pParseCommon->platelist_handle,"NOPHOT %s %s\n",totallistPlate,agreeString);
        }
      }
    }

  }
} /* End of PrintPlateList() */

void PrintUpdateList(PPARSECOMMON pParseCommon) {
  PPLATESTATS pPlateStats;
  int plateNumberIndex;
  int catalogIndex;
  int seriesIdIndex;
  int mosaicNumber = -1;
  int versionId = -1;
  int seriesId = -1;
  char catalogString[MAX_CATALOG_NUMBER];

  for (seriesIdIndex = 0; seriesIdIndex < MAX_SERIES; seriesIdIndex++) {
#ifdef DEBUG_SERIESID
    if (seriesIdIndex != DEBUG_SERIESID) {
      continue;
    }
#endif /* DEBUG_SERIESID */
    if (pParseCommon->plateStats[seriesIdIndex] == NULL) {
      continue;
    }    
    if (GetSequesteredFlag(seriesIdIndex) == SEQUESTERED_YES) {
#if 0
      printf("line %d seriesId %d id sequestered\n",__LINE__,seriesIdIndex);
#endif
      continue;
    }
    for (plateNumberIndex = 0; plateNumberIndex < MAX_PLATE_NUMBER; plateNumberIndex++) {            
      pPlateStats = &(pParseCommon->plateStats[seriesIdIndex])[plateNumberIndex];
      if (pPlateStats->plateNumber > 0) {


#ifdef DEBUG_PHOTINSERT
        if ((strcmp(pPlateStats->series,"mc") == 0) &&
            (pPlateStats->plateNumber == 23515)) {
          printf("\nline %4d for %s%05d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber);
        }
#endif /* DEBUG_PHOTINSERT */


        for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
          if ((catalogIndex == CATALOG_GAIA) ||
              (catalogIndex == CATALOG_KEPLER) ||
              (catalogIndex == CATALOG_EXPERIMENTAL)) {
            continue;
          }
#endif /* SKIP_GAIA */
          mosaicNumber = -1;
          versionId = -1;
          seriesId = -1;
        
          if ((pPlateStats->mosaicNumber[catalogIndex] >= 0) &&
              (pPlateStats->versionId[catalogIndex] > 0)) {
            if ((pPlateStats->photBestFlag == 1) &&
                (pPlateStats->photLatestFlag != 1)) {
              printf("ERROR: line %5d illegal condition for %s%05d catalogIndex %d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,catalogIndex);
            } else if ((pPlateStats->photBestFlag != 1) &&
                       (pPlateStats->photLatestFlag != 1)) {
              /* Nothing to do here */
            } else if ((pPlateStats->photBestFlag != 1) &&
                       (pPlateStats->photLatestFlag == 1)) {
              if (pPlateStats->versionId[catalogIndex] > pPlateStats->photLatestVersionId[catalogIndex]) {
                mosaicNumber = pPlateStats->mosaicNumber[catalogIndex];
                versionId = pPlateStats->versionId[catalogIndex];
              }
            } else { 
              if ((pPlateStats->photBestFlag != 1) ||
                  (pPlateStats->photLatestFlag != 1)) {
                printf("ERROR: line %5d illegal condition for %s%05d catalogIndex %d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,catalogIndex);
              }
              if (pPlateStats->photBestMosaicNumber == pPlateStats->photLatestMosaicNumber) {
                if (pPlateStats->mosaicNumber[catalogIndex] >= 0) {
                  if (pPlateStats->photBestMosaicNumber ==  pPlateStats->mosaicNumber[catalogIndex]) {
                    if (pPlateStats->versionId[catalogIndex] > pPlateStats->photLatestVersionId[catalogIndex]) {
                      mosaicNumber = pPlateStats->mosaicNumber[catalogIndex];
                      versionId = pPlateStats->versionId[catalogIndex];
                    } /* else: nothing needed to do */
                  } else {   
                    if (pPlateStats->photLatestMosaicNumber ==  pPlateStats->mosaicNumber[catalogIndex]) {
                      printf("ERROR: line %5d illegal condition for %s%05d catalogIndex %d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,catalogIndex);
                    } else {
                      if (pPlateStats->versionId[catalogIndex] > pPlateStats->photLatestVersionId[catalogIndex]) {
                        mosaicNumber = pPlateStats->photLatestMosaicNumber;
                        versionId = pPlateStats->versionId[catalogIndex];
                      } else {
                        if (pPlateStats->versionId[catalogIndex] <= pPlateStats->resortListVersionId[catalogIndex]) {
                          printf("ERROR: line %5d illegal condition for %s%05d catalogIndex %d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,catalogIndex);                        
                        } /* else: nothing needed to do */
                      }
                    }
                  }
                } else {
                  printf("ERROR: line %5d illegal condition for %s%05d catalogIndex %d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,catalogIndex);
                }
              }  else {
                printf("ERROR: line %5d illegal condition for %s%05d catalogIndex %d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,catalogIndex);
              }
            }
            
          }
          if ((mosaicNumber >= 0) &&
              (versionId > 0)) {
            if (catalogIndex == 0) {
              catalogString[0] = 0;
            } else {
              sprintf(catalogString,"%d",catalogIndex);
            }
            seriesId = GetSeriesId(pPlateStats->series,0);
            if (seriesId < 0) {
              printf("ERROR line %4d GetPlateStats bad series %s\n",__LINE__,pPlateStats->series);
              exit(-1);
            }


            fprintf(pParseCommon->out_handle,"UPDATE photplates set mosaicNumber = %d,versionId%s = %d  WHERE plateNumber = %d and seriesId = %d;\n",
                    pPlateStats->mosaicNumber[catalogIndex],
                    catalogString,
                    versionId,
                    pPlateStats->plateNumber,
                    seriesId);
          


          }
        }
      }
    }
  }
} /* End of PrintUpdateList() */

void ProcessResortList(PPARSECOMMON pParseCommon) {
  char *inBuffer;
  char inLine[MAX_BUFFER];
  int lineLen;
  int lineCount = 0;
  char *spacePtr;
  char copyLine[MAX_BUFFER];
  int plateNameLength;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  PPLATESTATS pPlateStats;
  int errorCount = 0;
  int catalogIndex;
  char *ingestDirectory;
  char plateListName[MAX_BUFFER];
  FILE *plateListHandle = NULL;
  char plateName[MAX_PLATE_NAME];
  int nvals;
  int curRow;
  
  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_INGEST is not defined\n");
    exit(-1);
  }

  for (catalogIndex = 0; catalogIndex < MAX_CATALOG_NUMBER; catalogIndex++) {
#ifdef SKIP_GAIA
    if ((catalogIndex == CATALOG_GAIA) ||
        (catalogIndex == CATALOG_EXPERIMENTAL)) {
      continue;
    }
#endif /* SKIP_GAIA */
    curRow = 0;
    lineCount = 0;
    strcpy(plateListName,ingestDirectory);
    strcat(plateListName,"/plates_");
    if (catalogIndex == CATALOG_GSC232) {
      strcat(plateListName,"gsc");
    } else {
      strcat(plateListName,catalogText[catalogIndex]);
    }
    strcat(plateListName,".txt");
    plateListHandle = fopen(plateListName,"rt");
    if (plateListHandle == NULL) {
      printf("ERROR: line %4d failed to open plate list name %s\n",__LINE__,plateListName);
      exit(-1);
    }
    
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,plateListHandle);
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
      strcpy(copyLine,inLine);
      lineCount++;
      if (lineCount <= 2) {
        continue;
      }
      spacePtr = strchr(inLine,'\t');
      if (spacePtr == NULL) {
        printf("ERROR:  line %d missing versionId with %d in %s for %s\n",__LINE__,copyLine,plateListName);
        exit(-1);
      }
      *spacePtr = 0;
      spacePtr++;
      plateNameLength = strlen(inLine);
      if (plateNameLength > (MAX_PLATE_NAME-1)) {
        printf("ERROR: line %d MAX_PLATE_NAME exceeded with %d in %s for %s\n",__LINE__,plateNameLength,copyLine,plateListName);
        exit(-1);
      }
      strcpy(plateName,inLine);
      if (ParseFilename2(inLine,series,&plateNumber) == 0) {
        printf("ERROR line %d ProcessResortList failed to parse lineCount %d %s in %s\n",__LINE__,lineCount,copyLine,plateListName);
        exit(-1);
      }
#ifdef DEBUG_SERIES
      if (strcmp(series,DEBUG_SERIES) != 0) {
        continue;
      }
#endif /* DEBUG_SERIES */

      pPlateStats = GetPlateStats(pParseCommon,series,plateNumber);
      if (pPlateStats == NULL) {
        printf("ERROR line %d GetPLateStats failed with lineCount %d %s in %s \n",__LINE__,lineCount,copyLine,plateListName);
        exit(-1);
      }
      if (pPlateStats->plateNumber != plateNumber) {
        strcpy(pPlateStats->series,series);
        pPlateStats->plateNumber = plateNumber;
      }
      nvals = sscanf(spacePtr,"%d",&pPlateStats->resortListVersionId[catalogIndex]);
      if (nvals != 1) {
        if (pPlateStats == NULL) {
          printf("ERROR line %d GetPLateStats failed with lineCount %d %s\n",__LINE__,lineCount,copyLine,plateListName);
          exit(-1);
        }
      }
      curRow++;
    }
    fclose(plateListHandle);
    plateListHandle = NULL;
    printf("ProcessResortList processed %d of %d lines for catalog %s\n",curRow,lineCount,catalogText[catalogIndex]);
  }
  return;
} /* End of ProcessResortList() */
void ReadPhotInsertTable(PPARSECOMMON pParseCommon) {
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  PPHOTVERSION pPhotVersion;
  int nvals;
  unsigned long numRows;
  unsigned long curRow = 0;
  int seriesId;
  int plateNumber;
  char series[MAX_SERIES_STRING];
  char queryString[MAX_QUERY_STRING];
  PPLATESTATS pPlateStats;
  int mosaicNumber;
  int catalogNumber;
  int versionId;
  int photFailFlag;
  int highestVersionId = -1;
  PPHOTGLOBAL pPhotGlobal = &pParseCommon->photGlobal;

#ifdef DEBUG_SERIESID
  sprintf(queryString,"SELECT seriesId,plateNumber,catalogNumber,mosaicNumber,versionId,photFailFlag from photinsert where seriesId = %d\n",DEBUG_SERIESID);
#else
  sprintf(queryString,"SELECT seriesId,plateNumber,catalogNumber,mosaicNumber,versionId,photFailFlag from photinsert\n;");
#endif
  res = ExecuteQuery(pParseCommon->pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pParseCommon->pPhotConnection);
    numRows = mysql_num_rows(res_ptr);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&seriesId);
          if (nvals != 1) {
            printf("ERROR: ReadPhotInsertTable line %d failed to parse seriesId\n",__LINE__);
            exit(-1);
          }
        } else {
          printf("ERROR: ReadPhotInsertTable line %d missing seriesId\n",__LINE__);
          exit(-1);    
        }
        if (GetSequesteredFlag(seriesId) == SEQUESTERED_YES) {
#if 0
          printf("line %d seriesId %d id sequestered\n",__LINE__,seriesId);
#endif
          continue;
        }
        strcpy(series,GetSeriesString(seriesId,1));

        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%d",&plateNumber);
          if (nvals != 1) {
            printf("ERROR: ReadPhotInsertTable line %d failed to parse plateNumber\n",__LINE__);
            exit(-1);
          }
        } else {
          printf("ERROR: ReadPhotInsertTable line %d missing plateNumber\n",__LINE__);
          exit(-1);    
        }
        pPlateStats = GetPlateStats(pParseCommon,series,plateNumber);
        if (pPlateStats == NULL) {
          printf("ERROR line %d GetPLateStats failed for mosaic %s%05d\n",__LINE__,series,plateNumber);
          exit(-1);
        }
        if (pPlateStats->plateNumber != plateNumber) {
          strcpy(pPlateStats->series,series);
          pPlateStats->plateNumber = plateNumber;
        }

        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&catalogNumber);
          if (nvals != 1) {
            printf("ERROR: ReadPhotInsertTable line %d failed to parse catalogNumber for %s%05d\n",__LINE__,series,plateNumber);
            exit(-1);
          }
        }
        if ((catalogNumber < 0) || (catalogNumber >= MAX_CATALOG_NUMBER)) {
          printf("ERROR: ReadPhotInsertTable line %d found incorrect catalog number %d for %s%05d\n",__LINE__,catalogNumber,series,plateNumber);
          exit(-1);
        }
      

        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&mosaicNumber);
          if (nvals != 1) {
            printf("ERROR: ReadPhotInsertTable line %d failed to parse mosaicNumber for %s%05d\n",__LINE__,series,plateNumber);
            exit(-1);
          }
        }
        if ((mosaicNumber < 0) || (mosaicNumber >= MAX_SCAN_MOSAIC_MOSAIC_NUMBER)) {
          printf("ERROR: ReadPhotInsertTable line %d illegal mosaicNumber %d for %s%05d\n",__LINE__,mosaicNumber,series,plateNumber);
          exit(-1);
        }

        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&versionId);
          if (nvals != 1) {
            printf("ERROR: ReadPhotInsertTable line %d failed to parse versionId\n",__LINE__);
            exit(-1);
          }
        }
        if (versionId < 0) {
          printf("ERROR: ReadPhotInsertTable line %d illegal versionId %d for %s%05d\n",__LINE__,versionId,series,plateNumber);
          exit(-1);
        }
        if (versionId > pPhotGlobal->currentVersion) {
          printf("ERROR: ReadPhotInsertTable line %d illegal high versionId %d currentVersion %d for %s%05d\n",__LINE__,versionId,pPhotGlobal->currentVersion,series,plateNumber);
          exit(-1);
        }
        if (versionId > highestVersionId) {
          highestVersionId = versionId;
        }


        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%d",&photFailFlag);
          if (nvals != 1) {
            printf("ERROR: ReadPhotInsertTable line %d failed to parse photFailFlag for %s%05d\n",__LINE__);
            exit(-1);
          }
        }
        if ((photFailFlag != 0) & (photFailFlag != 1)) {
          printf("ERROR: ReadPhotInsertTable line %d illegal photFailFlag %d for %s%05d\n",__LINE__,photFailFlag,series,plateNumber);
          exit(-1);
        }
        if (pPlateStats->versionId[catalogNumber] >= 0) {
          printf("ERROR: ReadPhotInsertTable line %d catalog %d already initialized for %s%05d\n",__LINE__,catalogNumber,series,plateNumber);
          exit(-1);
        }
        pPlateStats->versionId[catalogNumber] = versionId;
        pPlateStats->mosaicNumber[catalogNumber] = mosaicNumber;
        pPlateStats->photFailFlag[catalogNumber] = photFailFlag;
        if (pParseCommon->maxVersionId < versionId) {
          pParseCommon->maxVersionId = versionId;
        }
       

        curRow++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
    }
  }                                                                                                                        
  printf("ReadPhotInsertTable found %d entries. highestVersionId %d currentVersion %d in line %d\n",curRow,highestVersionId,pPhotGlobal->currentVersion,__LINE__);

  return;
} /* End of ReadPhotInsertTable() */

int main(int argc,char *argv[])
{
  char *argstr;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  char loglist_name[MAX_INPUT_NAME];
  char out_name[MAX_INPUT_NAME];
  char merge_name[MAX_INPUT_NAME];
  char platelist_name[MAX_INPUT_NAME];
  char totallist_name[MAX_INPUT_NAME];
  FILE *loglist_handle = NULL;
  int errorFlag = 0;
  char cmdchar;
  int lineLen;
  time_t curTime;
  PARSECOMMON parsecommon;
  PPARSECOMMON pParseCommon = &parsecommon;
  int stepIndex;
  time_t pTotalStepTime = 0;
  int stepCounter;
  double percentage;
  double percentage2;
  time_t averageTime;
  time_t totalAverageTime = 0;
  int solutionIndex;
  char *mysqlhost;
  char *username;
  char *password;
  char *photometry;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;

  char series[MAX_SERIES_STRING];
  int mosaicListResult = 0;
  int plateNumber = -1;
  int mosaicNumber = -1;
  int solutionNumber = 0;
  int mosaicInfoFlag = 1;
  int verbose = 0; /* for BuildMosaicList only */
  int selectAllFlag = 0;
  int pendingFlag = 0;
  int totalSolution0Count = 0;
  int totalSolutionsCount = 0;
  int numMosaicRecords = 0;
  int mosaicRecordIndex = 0;
  int numScannedMosaics = 0;
  PMOSAICLIST pMosaicTableEntry;
  PMOSAICLIST pMosaicTable = NULL;
  PPLATESTATS pPlateStats;
#ifdef DEBUG_SERIES
  printf("ERROR: DEBUG_SERIES is %s\n",DEBUG_SERIES);
#endif /* DEBUG_SERIES */
#ifdef DEBUG_SERIESID
  printf("ERROR: DEBUG_SERIESID is %d\n",DEBUG_SERIESID);
#endif /* DEBUG_SERIESID */
#ifdef DEBUG_NO_LOGS
  printf("ERROR: DEBUG_NO_LOGS is %d\n",DEBUG_NO_LOGS);
#endif /* DEBUG_NO_LOGS  */
#ifdef DEBUG_PHOTINSERT
  printf("ERROR: DEBUG_PHOTINSERT is %d\n",DEBUG_PHOTINSERT);
#endif /* DEBUG_PHOTINSERT  */

  series[0] = 0;


  /* Loop through the arguments */
  loglist_name[0] = 0;
  out_name[0] = 0;
  platelist_name[0] = 0;
  totallist_name[0] = 0;
  merge_name[0] = 0;
  memset(pParseCommon,0,sizeof(PARSECOMMON));
  pParseCommon->pConnection = &pParseCommon->my_connection;
  pParseCommon->pPhotConnection = &pParseCommon->my_phot_connection;


#ifdef DEBUG_LOGS_ONLY
  printf("ERROR: DEBUG_LOGS_ONLY is %d\n",DEBUG_LOGS_ONLY);
  pParseCommon->debugLogsOnly = 1;
#endif /* DEBUG_LOGS_ONLY  */

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

        case 's':
        case 'S':
          pParseCommon->shortTotalForm = 1;
          break;

        case 'n':
        case 'N':
          pParseCommon->parseLogsOnly = 1;
          pParseCommon->debugLogsOnly = 1;
          break;

        case 'j':
        case 'J':
          pParseCommon->jpegissueForm = 1;
          break;

        case 'k':
        case 'K':
          pParseCommon->shortNameForm = 1;
          break;


        case 'l': /* list file name */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(loglist_name,*++argv,MAX_INPUT_NAME-2);
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

        case 'm': /* merge file name */
        case 'M':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(merge_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
              exit(-1);
            } else {
              pParseCommon->outputMergeFlag = 1;
            }

          }
          break;




        case 'p': /* plate list name */
        case 'P':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(platelist_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 't': /* total list name */
        case 'T':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(totallist_name,*++argv,MAX_INPUT_NAME-2);
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
  
  if (totallist_name[0] == 0) {
    printf("ERROR: Totallist file name not specified\n");
    errorFlag = 1;
  } else {
    pParseCommon->totallist_handle = fopen(totallist_name,"rt");
    if (pParseCommon->totallist_handle == NULL) {
      errorFlag = 1;
      printf("Could not open totallist file %s\n",totallist_name);
    }
    if ((loglist_name[0] != 0) ||
        (platelist_name[0] != 0)) {

      if (loglist_name[0] == 0) {
        printf("ERROR: List file name not specified\n");
        errorFlag = 1;
      } else {
        loglist_handle = fopen(loglist_name,"rt");
        if (loglist_handle == NULL) {
          errorFlag = 1;
          printf("Could not open list file %s\n",loglist_name);
        }
      }

      if (platelist_name[0] == 0) {
        printf("ERROR: Platelist file name not specified\n");
        errorFlag = 1;
      } else {
        pParseCommon->platelist_handle = fopen(platelist_name,"wt");
        if (pParseCommon->platelist_handle == NULL) {
          errorFlag = 1;
          printf("Could not open platelist file %s\n",platelist_name);
        }
      }
    }
  }
  
  if (out_name[0] != 0) {
    pParseCommon->out_handle = fopen(out_name,"wt");
    if (pParseCommon->out_handle == NULL) {
      errorFlag = 1;
      printf("Could not open out file %s\n",out_name);
    }
  }

  if (pParseCommon->outputMergeFlag) {
    if (merge_name[0] == 0) {
      printf("ERROR: Merge file name not specified\n");
      errorFlag = 1;
    } else {
      pParseCommon->merge_handle = fopen(merge_name,"wt");
      if (pParseCommon->merge_handle == NULL) {
        errorFlag = 1;
        printf("Could not open merge file %s\n",merge_name);
      }
    }
  }

  if (errorFlag) {
    printf("Usage: parsephot [-v] -l <list of log files> -o <outfile> -p <platelist> -t <total.list> -s \n");
    printf("       where <total.list> is produced by dasch_scat\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s is a short form total.list, ignoring rotation and binning\n");
    printf("       where -j is used for jpegissue.txt, ignoring _p and _j\n");
    printf("       where -m <merge file> produces a merger of the two files, including excess information\n");
    printf("       where -k checks only series and plateNumber\n");
    printf("       where -n lists only the plates parsed in the logfile\n");
    return(-1);
  }

  printf("parsephot of %s %s List Filename %s Output Filename %s\n",
         __DATE__,__TIME__,loglist_name,out_name);
 

  time(&pParseCommon->startTime);

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
  mysql_init(pParseCommon->pConnection);
  if (!mysql_real_connect(pParseCommon->pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pParseCommon->pConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pParseCommon->pConnection),mysql_error(pParseCommon->pConnection));
    }
    return(-1);
  }
  photometry = getenv("DASCH_PHOTOMETRY");
  if (password == NULL) {
    printf("ERROR: DASCH_PHOTOMETRY is not defined\n");
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
  mysql_init(pParseCommon->pPhotConnection);
  if (!mysql_real_connect(pParseCommon->pPhotConnection,mysqlphothost,photusername,photpassword,photometry,0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pParseCommon->pPhotConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pParseCommon->pPhotConnection),mysql_error(pParseCommon->pPhotConnection));
    }
    return(-1);
  }

  InitSeriesTable(pParseCommon->pConnection,pParseCommon->pPhotConnection);
  InitVersionTable(pParseCommon);
  if (GetPhotometryGlobal(pParseCommon->pPhotConnection,&pParseCommon->photGlobal) != 1) {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }



  if ((loglist_name[0] == 0) &&
      (pParseCommon->parseLogsOnly == 0) &&
      (platelist_name[0] == 0)) {
    
    ProcessTotalList(pParseCommon);
    PrintPlateList(pParseCommon);
    if (pParseCommon->merge_handle != NULL) {
      fclose(pParseCommon->merge_handle);
      pParseCommon->merge_handle = NULL;
    }
    if (pParseCommon->totallist_handle != NULL) {
      fclose(pParseCommon->totallist_handle);
      pParseCommon->totallist_handle = NULL;
    }

    return(EXIT_SUCCESS);
  }

  /* Execute BuildMosaicLists() to get a list of all of the preferred scans */
#ifdef DEBUG_SERIES 
  strcpy(series,DEBUG_SERIES);
#else /* DEBUG_SERIES */
  series[0] = 0;
#endif /* DEBUG_SERIES */
  if (pParseCommon->debugLogsOnly == 0) { /* ifndef debugLogsOnly */
    if (pParseCommon->parseLogsOnly == 0) {
      mosaicListResult = BuildMosaicList(pParseCommon->pConnection,
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
                                         &pMosaicTable);
    } 
    for (mosaicRecordIndex = 0; mosaicRecordIndex < numMosaicRecords; mosaicRecordIndex++) {
      pMosaicTableEntry = &pMosaicTable[mosaicRecordIndex];
#ifdef DEBUG_SERIES
      if (strcmp(pMosaicTableEntry->series,DEBUG_SERIES) != 0) {
        continue;
      }
#endif /* DEBUG_SERIES */

      if (strcmp(pMosaicTableEntry->series,"mask") == 0) {
        continue;
      }
      pPlateStats = GetPlateStats(pParseCommon,pMosaicTableEntry->series,pMosaicTableEntry->plateNumber);
      if (pPlateStats == NULL) {
        printf("ERROR line %4d GetPLateStats failed\n",__LINE__);
        exit(-1);
      }
      if (pPlateStats->bestFlag != 0) {
        printf("ERROR line %4d GetPLateStats bestFlag already initialized\n",__LINE__);
        exit(-1);
      }
      if (pPlateStats->series[0] != 0) {
        printf("ERROR line %4d GetPLateStats series already initialized\n",__LINE__);
        exit(-1);
      }
      if (pPlateStats->plateNumber > 0) {
        printf("ERROR line %4d GetPLateStats plateNumber already initialized\n",__LINE__);
        exit(-1);
      }
      strcpy(pPlateStats->series,pMosaicTableEntry->series);
      pPlateStats->plateNumber = pMosaicTableEntry->plateNumber;
      pPlateStats->bestFlag = 1;
      pPlateStats->bestMosaicNumber = pMosaicTableEntry->mosaicNumber;
      pPlateStats->bestRotation =  pMosaicTableEntry->rotation;
      pPlateStats->bestBinning = pMosaicTableEntry->binning;
      pPlateStats->bestWCSSource = pMosaicTableEntry->WCSSource;
      strcpy(pPlateStats->bestMosaicComment,pMosaicTableEntry->mosaicComment);
      numScannedMosaics++;

    }

    if (pMosaicTable != NULL) {
      free(pMosaicTable);
      pMosaicTable = NULL;
    }
    
    printf("Found %6d scanned mosaics at line %d\n",numScannedMosaics,__LINE__);
    if (pParseCommon->parseLogsOnly == 0) {
      ProcessResortList(pParseCommon);
      ReadPhotPlatesTable(pParseCommon);
    }
  } /* ifndef debugLogsOnly */
  if (pParseCommon->debugLogsOnly == 0) { /* ifndef debugLogsOnly */
    if (pParseCommon->parseLogsOnly == 0) {
      ReadPhotInsertTable(pParseCommon);
    }
  } else { /* ifndef debugLogsOnly */
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,loglist_handle);
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
#ifndef DEBUG_NO_LOGS
      if (inBuffer[0] != '#') {
        ProcessLogBuffer(pParseCommon,inBuffer);
      }
#endif /* DEBUG_NO_LOGS */
    }
  } /* ifndef debugLogsOnly */


  if (pParseCommon->parseLogsOnly == 0) {
    if (pParseCommon->totallist_handle != NULL) {
      ProcessTotalList(pParseCommon);
    }
    if (pParseCommon->platelist_handle != NULL) {
      PrintPlateList(pParseCommon);
    }

    PrintUpdateList(pParseCommon);
  }

  if (loglist_handle != NULL) {
    fclose(loglist_handle);
  }
  if (pParseCommon->out_handle != NULL) {
    fclose(pParseCommon->out_handle);
    pParseCommon->out_handle = NULL;
  }
  if (pParseCommon->merge_handle != NULL) {
    fclose(pParseCommon->merge_handle);
    pParseCommon->merge_handle = NULL;
  }
  if (pParseCommon->platelist_handle != NULL) {
    fclose(pParseCommon->platelist_handle);
    pParseCommon->platelist_handle = NULL;
  }
  if (pParseCommon->totallist_handle != NULL) {
    fclose(pParseCommon->totallist_handle);
    pParseCommon->totallist_handle = NULL;
  }
  time(&curTime);
  curTime -= pParseCommon->startTime;



  printf("Execution Time: %lld seconds. maxLineLen %d maxVersionId %d\n",
         curTime,
         pParseCommon->maxLineLen,
         pParseCommon->maxVersionId);





  return(EXIT_SUCCESS);
} /* end of main() */

int ProcessLogBuffer(PPARSECOMMON pParseCommon, char* log_name) {
  FILE *log_handle = NULL;
  char *inBuffer;
  int nvals;
  int lineLen;
  int lineCount = 0;
  char copyLine[MAX_BUFFER];
  char tableLoadLine1[MAX_BUFFER];
  char tableLoadLine2[MAX_BUFFER];
  int tableLoadFlag = 0;
  char inLine[MAX_BUFFER];
  char parseLine[MAX_BUFFER];
  char *parsePtr1;
  char *parsePtr2;
  size_t tmpBytesCopied;
  int parseCount;
  int solutionIndex;
  char token[MAX_BUFFER];
  char startToken[MAX_BUFFER];
  int recStatCount = 0;
  int parseState = 1; /* 1 = idle; 2 = in progress ; 3 = flush */
  char *charPtr;
  char *charPtr2;
  char *charPtr3;
  int photFailFlag = 0;
  char series[MAX_SERIES_STRING];
  int seriesId;
  int plateNumber;
  int mosaicNumber;
  int rotation;
  int binning;
  int processCount = 0;
  time_t lastTime = 0;
  time_t tempTime = 0;
  int lastVersionId = -1;
  int tempVersionId = -1;
  int catalogNumber = -1;
  int parseErrorCount = 0;
  char prevPlate[2*MAX_PLATE_NAME];
  char curPlate[2*MAX_PLATE_NAME];
  PPLATESTATS pPlateStats;
  int plateNumberIndex;
  int catalogIndex;
  int logFileNumber = -1;   /* <nnn> in phot_<nnn><v>.log */
  int logFileVersion = -1;  /* <v> in phot_<nnn><v>.log  where 0 = '', 1 = 'a' ... */
  char charVal;
  prevPlate[0] = 0;
  curPlate[0] = 0;

  if (pParseCommon->verbose) {  
    printf("Processing %s\n",log_name);
  }

  strcpy(copyLine,log_name);
  charPtr = strstr(copyLine,"phot_");
  if (charPtr == NULL) {
    printf("ERROR: line %4d could not parse file name %s\n",__LINE__,log_name);
    return(-1);
  }
  charPtr += strlen("phot_");
  nvals = sscanf(charPtr,"%d",&logFileNumber);
  if (nvals == 0) {
    printf("ERROR: line %4d could not parse file name %s\n",__LINE__,log_name);
    return(-1);
  }
  charVal = *charPtr;
  while ((charVal >= '0') && (charVal <= '9')) {
    charPtr++;
    charVal = *charPtr;
  }
  if (charVal == '.') {
    logFileVersion = 0;
  } else {
    logFileVersion = charVal-'a'+1;
    if ((logFileVersion < 1) ||
        (logFileVersion > 26)) {
      printf("ERROR: line %4d could not parse file name %s\n",__LINE__,log_name);
      return(-1);
    }
    charPtr++;
    charVal = *charPtr;
    if (charVal != '.') {
      printf("ERROR: line %4d could not parse file name %s\n",__LINE__,log_name);
      return(-1);
    }
  }




  log_handle = fopen(log_name,"rt");
  startToken[0] = 0;
  if (log_handle == NULL) {
    printf("ERROR: line %4d Could not open list file %s\n",__LINE__,log_name);
    return(-1);
  }
  tableLoadFlag = 0;
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
    charPtr = strstr(inBuffer,"table_load: missing column: rejectFlag");
    if (charPtr != NULL) {
      if (tableLoadFlag != 0) {
        *charPtr = 0;
        charPtr += strlen("table_load: missing column: rejectFlag");
        charPtr2 = strstr(charPtr,"table_load: missing column: rejectFlag");
        if (charPtr2 != NULL) {
          printf("ERROR line %d lineCounter %d double table_load (1) %s in %s\n",__LINE__,lineCount,inBuffer,log_name);
          continue;
        }
        strcpy(tableLoadLine2,inBuffer);
        lineLen = strlen(tableLoadLine1)+strlen(tableLoadLine2);
        if (lineLen >= (MAX_BUFFER-1)) {
          printf("ERROR line %d MAX_BUFFER exceeded: length %d %s\n",__LINE__,lineLen,copyLine);
          exit(-1);
        }
        strcat(tableLoadLine1,tableLoadLine2);
        continue;
      } else {
        *charPtr = 0;
        charPtr += strlen("table_load: missing column: rejectFlag");
        charPtr2 = strstr(charPtr,"table_load: missing column: rejectFlag");
        if (charPtr2 != NULL) {
          printf("ERROR line %d lineCounter %d double table_load (2) %s in %s\n",__LINE__,lineCount,inBuffer,log_name);
          continue;
        }
        strcpy(tableLoadLine1,inBuffer);
        tableLoadFlag = 1;
        continue;
      }
    }
    if (tableLoadFlag != 0) {
      tableLoadFlag = 0;
      strcpy(tableLoadLine2,inBuffer);
      lineLen = strlen(tableLoadLine1)+strlen(tableLoadLine2);
      if (lineLen >= (MAX_BUFFER-1)) {
        printf("ERROR line %d MAX_BUFFER exceeded: length %d %s\n",__LINE__,lineLen,copyLine);
        exit(-1);
      }
      strcat(tableLoadLine1,tableLoadLine2);
      strcpy(inLine,tableLoadLine1);
    }
#if 0
    if (lineCount == 138) {
      printf("At line %d %s\n",lineCount,inLine);
    }
#endif
    strcpy(copyLine,inLine);
    strcpy(parseLine,inLine);
    if (ParseTime(pParseCommon,parseLine,&tempTime,&tempVersionId,lineCount) == 0) {
      /* We have a valid time for this file */
      lastTime = tempTime;
      lastVersionId = tempVersionId;
    }
    if ((pParseCommon->verbose != 0) && 
        (strstr(inLine,"New version ") == inLine)) {
      printf("%s in line %d\n",inLine,lineCount);
    }
    if (strstr(inLine,"update_photometry of ") == inLine) {
      charPtr = strstr(parseLine,"catalog ");
      if (charPtr != NULL) {
        charPtr += strlen("catalog ");
      
        for (catalogNumber = 0; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
          charPtr2 = strstr(charPtr,catalogText[catalogNumber]);
          if (charPtr2 != NULL) {
            break;
          }
        }
        if (catalogNumber >= MAX_CATALOG_NUMBER) {
          printf("ERROR: unrecognized catalog in line %d %s\n",lineCount,copyLine);
          catalogNumber = -1;
        } else {
          if (pParseCommon->verbose != 0) {
            printf("Found catalog in line %d catalog number %d %s %s\n",__LINE__,catalogNumber,catalogText[catalogNumber],copyLine);
          }
        }
      }

    } else if ((strstr(inLine,"resort_magfiles of ") == inLine) ||
               (strstr(inLine,"update_limiting of ") == inLine) ||
               (strstr(inLine,"search_none of ") == inLine) ||
               (strstr(inLine,"update_summary2 of ") == inLine) ||
               (strstr(inLine,"Execution Time: ") == inLine)) {
      if (catalogNumber >= 0) {
        if (pParseCommon->verbose != 0) {
          printf("Finished in line %d for catalogNumber %d %s %s\n",lineCount,catalogNumber,catalogText[catalogNumber],copyLine);
        }
      }
      catalogNumber = -1;
      prevPlate[0] = 0;
      curPlate[0] = 0;
    }
    photFailFlag = 0;
    if ((catalogNumber >= 0) &&
        ((strstr(inLine,"ERROR: Failed to find spatial bins file") != NULL) ||
         (strstr(inLine,"scale ") == inLine) ||
         (strstr(inLine,"rate ") == inLine))) {
      charPtr = NULL;
      charPtr3 = strstr(inLine,"ERROR: Failed to find spatial bins file");
      if (charPtr3 != NULL) {
        charPtr = strrchr(parseLine,'/');
        if (charPtr != NULL) {
          charPtr++;
          photFailFlag = 1;

        }
      } else {
        charPtr = strstr(parseLine,"for ");
        if (charPtr != NULL) {
          charPtr += strlen("for ");
        }
      }
      if (charPtr != NULL) {
        if (ParseFilename(charPtr,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 1) {
          if (rotation == 0) {
            sprintf(curPlate,"%05s%05d_%02d_%02dww",series,plateNumber,mosaicNumber,binning);
          } else {
            sprintf(curPlate,"%05s%05d_%02d_%02dr%dww",series,plateNumber,mosaicNumber,binning,rotation);
          }
          if (strlen(curPlate) >= MAX_PLATE_NAME) {
            printf("ERROR line %d MAX_PLATE_NAME exceeded in  %s\n",__LINE__,copyLine);
            exit(-1);
          }
          if (pParseCommon->parseLogsOnly) {
            fprintf(pParseCommon->out_handle,"%s\n",curPlate);
            continue;
          }

#ifdef DEBUG_SERIES
          if (strcmp(series,DEBUG_SERIES) != 0) {
            continue;
          }
#endif /* DEBUG_SERIES */

          if (pParseCommon->verbose != 0) {
            printf("line %d found %s photFailFlag %d catalogNumber %d lastVersionId %d phot_%d_%d in %s\n",__LINE__,curPlate,photFailFlag,catalogNumber,lastVersionId,logFileNumber,logFileVersion,copyLine);
          }
          if (strcmp(curPlate,prevPlate) != 0) {
            strcpy(prevPlate,curPlate);
            seriesId = GetSeriesId(series,0);
            if (seriesId > 0) {
              pPlateStats = GetPlateStats(pParseCommon,series,plateNumber);
              if (pPlateStats->plateNumber != plateNumber) {
                strcpy(pPlateStats->series,series);
                pPlateStats->plateNumber = plateNumber;
              }
#if 0
              printf("PHOTINSERT logFileNumber %4d logFileVersion %2d catalogNumber %d series %5s plateNumber %5d mosaicNumber %2d versionId %3d photFailFlag %d timestamp %7d\n",
                     logFileNumber,
                     logFileVersion,
                     catalogNumber,
                     series,
                     plateNumber,
                     mosaicNumber,
                     lastVersionId,
                     photFailFlag,
                     lastTime);
#endif
              if (pPlateStats->timestamp[catalogNumber] < lastTime) {
                pPlateStats->versionId[catalogNumber] = lastVersionId;
                pPlateStats->mosaicNumber[catalogNumber] = mosaicNumber;
                pPlateStats->rotation[catalogNumber] = rotation;
                pPlateStats->binning[catalogNumber] = binning;
                pPlateStats->timestamp[catalogNumber] = lastTime;
                pPlateStats->photFailFlag[catalogNumber] = photFailFlag;
                pPlateStats->logFileNumber[catalogNumber] = logFileNumber;
                pPlateStats->logFileVersion[catalogNumber] = logFileVersion;
                if (pParseCommon->maxVersionId < lastVersionId) {
                  pParseCommon->maxVersionId = lastVersionId;
                }
              }
            } else {
#if 0
              printf("ERROR line %d bad series in %s\n",__LINE__,copyLine);
#endif
            }

          }
          parseErrorCount = 0;
        } else {
          parseErrorCount++;
          if (parseErrorCount < 2) {
#if 0
            printf("WARNING line %d failed to parse %s\n",__LINE__,copyLine);
#endif
          } else {
            printf("ERROR line %d failed to parse %s\n",__LINE__,copyLine);
          }
        }
      }
    }
  }
  if (pParseCommon->verbose) {
    printf("Finished processing %s lines %d recStatCount %d\n",log_name,lineCount,recStatCount);
  }
  fclose(log_handle);
  return(0);
} /* end of ProcessLogBuffer() */

