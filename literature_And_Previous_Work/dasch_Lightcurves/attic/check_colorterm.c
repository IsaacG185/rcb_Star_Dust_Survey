// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* check_colorterm.c
 *
 *  This program examines all of the colorterm values for each plate under the expectation that they should all be the same unless a systematic error occurs
 *  SELECT series,seriesId,plateNumber,spatial_bin,solutionNumber,colorterm,errorcolor from spatialbin2 inner join photseries using(seriesId) where colorflag = "metropolis";
 *   gives 1093001 spatial bins on April 24, 2017 for the apass calibration
 *                  mean          SD        Minimum       Maximum
 *     colorterm  -0.007095     0.255834     -2.234        1.659 
 *     errorcolor  0.002974     0.004814       0           0.1
 *
 *     for 1078402 entries, RMS_0 of errorcolor = 0.00564043
 *
 *     For DNR series
 *     colorterm  -0.118233     0.085625     -1.338        0.408
 *     errorcolor  0.001953     0.003154      0.0          0.099
 *     for 40150 entries, RMS_0 of errorcolor =    0.0037095
 *     (errorcolor < 0.01) && (colorterm > -0.294) && (colorterm < 0.0575)
 *
 *     See memo of Sat, 6 May 2017 16:03:41 for results on Damons.  Graphs are damoncolortermvsdeclination.png  damoncolortermvsra.png  damonplatecenters.png in the /dasch/scanner/backup/2017_05_04 directory.
 * 
 * cc -ggdb -O0  -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include check_colorterm.c pipelineutils.a -L /dasch/install/lib -lm -lcfitsio  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread    -ltable -lutil  -lwcs plottransientstub.o -o check_colorterm

 *    
 *   check_colorterm -v -o /home/scanner/Pipeline/loscolorterm.log
 *
 * Apr 18, 2017 Edward J. Los - Initial version
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

#define LOS_DEBUG 1
#define LOS_SERIES "x"
#define LOS_PLATENUMBER 18967


#define MAX_LIST_STRING 25
#define MAX_FILENAME 512
#define MAX_BUFFER 256
#define NUM_OUT_FILES 2

typedef struct _plateentry {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int seriesId;
  int numSpatialBins;
  int spatialBinMask;
  int mosaicNumber;
  double colorterm[MAX_SPATIAL_BINS+1];
  double errorcolor[MAX_SPATIAL_BINS+1];
  double median9bin; /* Median for 9 spatial bins */
  double rms9bin;    /* rms for 9 spatial bins */
  double dRightAscension;
  double dDeclination;
} PLATEENTRY,*PPLATEENTRY;



int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char outfiles[NUM_OUT_FILES][MAX_BUFFER];
  int fileIndex;
  FILE *outHandles[NUM_OUT_FILES];
  int verbose = 0;
  char *pChar;
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
  int plateCount = 0;
  int plateCount2 = 0;
  int plateCountIncrement = 0;
  int plateAlloc = 0;
  int plateAlloc2 = 0;
  int plateCountNew = 0;
  int plateIndex;
  int plateIndex2;
  int plateTestIndex;
  int nineBinRMSCount = 0;
  int matchCount = 0;
  PPLATEENTRY pPlatesTable = NULL;
  PPLATEENTRY pPlateEntry;
  PPLATEENTRY pPrevEntry;
  PPLATEENTRY pPlatesTable2 = NULL;
  PPLATEENTRY pPlateEntry2;
  int result;
  int nvals;
  EXPOSURE exposureTable;
  PEXPOSURE pExposure = &exposureTable;
  time_t startTime;
  time_t curTime;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int FitWCS;
  int seriesId;
  int spatial_bin;
  int numSpatialBinCount[MAX_SPATIAL_BINS+1];
  double vector[MAX_SPATIAL_BINS+1];
  int vectorCount;
#ifdef DEFER_EXPOSURE_CHECK
  printf("ERROR: DEFER_EXPOSURE_CHECK is set\n");
#endif /* DEFER_EXPOSURE_CHECK */

  memset(numSpatialBinCount,0,sizeof(numSpatialBinCount));
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
    pChar = strchr(outfile,'.');
    if (pChar != NULL) {
      *pChar = 0;
    }
    for (fileIndex = 0;fileIndex < NUM_OUT_FILES; fileIndex++) {
      sprintf(outfiles[fileIndex],"%s_%d.db",outfile,fileIndex);
      


      outHandles[fileIndex] = fopen(outfiles[fileIndex],"wt");
      if (outHandles[fileIndex] == NULL) {
        errorFlag = 1;
        printf("ERROR: Failed to open the output file %s\n",outfiles[fileIndex]);
      } else {
        if (verbose) {
          printf("Output file %s\n",outfiles[fileIndex]);
        }
      }
    }
  }
  if (errorFlag) {
    printf("Usage: check_colorterm options\n");
    printf("  options: -v verbose\n");
    printf("           -o <output file>\n");
    return(-1);
  }

  printf("check_colorterm of %s %s\n",
         __DATE__,__TIME__);
  printf("line %d Loading the plates table \n",__LINE__);
  fprintf(outHandles[0],"Plate\tspatial_bin\tcolorterm\terrorcolor\n");
  fprintf(outHandles[0],"-----\t-----------\t---------\t----------\n");


  res = ExecuteQuery(pPhotConnection,"SELECT series,plateNumber,seriesId,spatial_bin,colorterm,errorcolor,mosaicNumber FROM spatialbin2 INNER JOIN photseries USING(seriesId) INNER JOIN photplates USING(seriesId,plateNumber) WHERE (series = 'dnb' or series = 'dsb') and colorflag = 'metropolis' AND solutionNumber = 0 ORDER BY series,plateNumber,spatial_bin;");
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pPhotConnection);
      printf("There are currently %d rows allocated to the plates table\n",numRows);
      plateAlloc = numRows;
      pPlatesTable = (PPLATEENTRY)calloc(plateAlloc,sizeof(PLATEENTRY));
      if (pPlatesTable == NULL) {
        printf("ERROR: failed to allocate %d PLATEENTRY rows\n",plateAlloc);
        exit(-1);
      }

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        plateCountIncrement = 0;
        pPlateEntry = &pPlatesTable[plateCount];
        memset(pPlateEntry,0,sizeof(PLATEENTRY));
        strcpy(pPlateEntry->series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&pPlateEntry->plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((pPlateEntry->plateNumber == LOS_PLATENUMBER) &&
            (strcmp(pPlateEntry->series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,pPlateEntry->series,pPlateEntry->plateNumber);
        }
#endif /* LOS_DEBUG */

        if (plateCount > 0) {
          pPrevEntry = &pPlatesTable[plateCount-1];
#ifdef LOS_DEBUG
          if ((pPrevEntry->plateNumber == LOS_PLATENUMBER) &&
              (strcmp(pPrevEntry->series,LOS_SERIES) == 0)) {
            printf("line %d plate %s%05d\n",__LINE__,pPrevEntry->series,pPrevEntry->plateNumber);
          }
#endif /* LOS_DEBUG */
          if ((pPrevEntry->plateNumber == pPlateEntry->plateNumber) &&
              (strcmp(pPrevEntry->series,pPlateEntry->series) == 0)) {
            pPlateEntry = pPrevEntry; /* Continuation of the previous plate */
          } else {
            plateCountIncrement++;
          } 
        } else {
          plateCountIncrement++;
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pPlateEntry->seriesId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for seriesId\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[6] != NULL) {
          nvals = sscanf(sqlrow[6],"%d",&pPlateEntry->mosaicNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for mosaicNumber\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&spatial_bin);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for spatial_bin\n");
            continue;
          }
          if ((spatial_bin < 1) || (spatial_bin > MAX_SPATIAL_BINS)) {
            printf("ERROR illegal spatial_bin %d\n",spatial_bin);
            continue;
          }
          if ((pPlateEntry->spatialBinMask & (1 << spatial_bin)) != 0) {
            printf("ERROR: duplicate spatial_bin %d\n",spatial_bin);
            continue;
          }
          nvals = sscanf(sqlrow[4],"%lf",&pPlateEntry->colorterm[spatial_bin]);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for colorterm\n");
            continue;
          }
          nvals = sscanf(sqlrow[5],"%lf",&pPlateEntry->errorcolor[spatial_bin]);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for errorcolor\n");
            continue;
          }
          pPlateEntry->numSpatialBins++;
          pPlateEntry->spatialBinMask |= (1 << spatial_bin);
          fprintf(outHandles[0],"%s%05d\t%d\t%f\t%f\n",
                  pPlateEntry->series,
                  pPlateEntry->plateNumber,
                  spatial_bin,
                  pPlateEntry->colorterm[spatial_bin],
                  pPlateEntry->errorcolor[spatial_bin]);


        } else {
          continue;
        }
        pPlateEntry->dRightAscension = 999.0;
        pPlateEntry->dDeclination = 99.0;

        plateCount += plateCountIncrement;
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
  time(&curTime);
  curTime -= startTime;
  printf("processed %d plates from the spatialbin table at %d seconds\n",plateCount,curTime);


  
  res = ExecuteQuery(pConnection,"SELECT series,plateNumber,dRightAscension,dDeclination FROM mosaics INNER JOIN exposures USING (series,plateNumber,exposureNumber)  WHERE  (series = 'dnb' or series = 'dsb')  AND solutionNumber = 0 AND centerSource = 'imWCS' ORDER BY series,plateNumber;");
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pConnection);
      printf("There are currently %d rows allocated to the plates2 table\n",numRows);
      plateAlloc2 = numRows + 1;
      pPlatesTable2 = (PPLATEENTRY)calloc(plateAlloc2,sizeof(PLATEENTRY));
      if (pPlatesTable2 == NULL) {
        printf("ERROR: failed to allocate %d PLATEENTRY rows\n",plateAlloc2);
        exit(-1);
      }

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        pPlateEntry2 = &pPlatesTable2[plateCount2];
        memset(pPlateEntry2,0,sizeof(PLATEENTRY));
        strcpy(pPlateEntry2->series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&pPlateEntry2->plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber\n");
          continue;
        }
#ifdef LOS_DEBUG
        if ((pPlateEntry2->plateNumber == LOS_PLATENUMBER) &&
            (strcmp(pPlateEntry2->series,LOS_SERIES) == 0)) {
          printf("line %d plate %s%05d\n",__LINE__,pPlateEntry2->series,pPlateEntry2->plateNumber);
        }
#endif /* LOS_DEBUG */

        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%lf",&pPlateEntry2->dRightAscension);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for seriesId\n");
            continue;
          }
        } else {
          continue;
        }

        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%lf",&pPlateEntry2->dDeclination);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for seriesId\n");
            continue;
          }
        } else {
          continue;
        }

        plateCount2++;
        if (plateCount2 >= plateAlloc2) {
          printf("ERROR: plateAlloc2  %d is too small\n",plateAlloc2);
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
  time(&curTime);
  curTime -= startTime;
  printf("processed %d plates from the exposures table at %d seconds\n",plateCount2,curTime);
  
  /* Now combine the two tables */
  plateIndex2 = 0;
  pPlateEntry2 = &pPlatesTable2[plateIndex2];
  for (plateIndex = 0; plateIndex < plateCount; plateIndex++) {
    pPlateEntry = &pPlatesTable[plateIndex];
    while (1) {
      if (plateIndex2 >= plateCount2) {
        break;
      }
      pPlateEntry2 = &pPlatesTable2[plateIndex2];
      result = strcmp(pPlateEntry->series,pPlateEntry2->series);
      if (result == 0) {
        if (pPlateEntry->plateNumber == pPlateEntry2->plateNumber) {
          result = 0;
        } else if (pPlateEntry->plateNumber > pPlateEntry2->plateNumber) {
          result = 1;
        } else {
          result = 0;
        }
      }
#if 0
      printf("result %d for %s %s and %d %d\n",result,pPlateEntry->series,pPlateEntry2->series,pPlateEntry->plateNumber,pPlateEntry2->plateNumber);
#endif
      if (result == 0) {
        matchCount++;
        pPlateEntry->dRightAscension = pPlateEntry2->dRightAscension;
        pPlateEntry->dDeclination = pPlateEntry2->dDeclination;
        plateIndex2++;
        break;
      } else if (result < 0) {
        break;
      } else {
        plateIndex2++;
      }
    }
  }
  printf("Total matches %d\n",matchCount);

#if 0
  fprintf(outHandles[1],"Plate\tmedian\trms\n");
  fprintf(outHandles[1],"-----\t------\t---\n");
#else
  fprintf(outHandles[1],"colorterm1\tcolorterm2\tcolorterm3\tcolorterm4\tcolorterm5\tcolorterm6\tcolorterm7\tcolorterm8\tcolorterm9\tPlate\tmedian\trms\tdRightAscension\tdDeclination\tseriesId\n");
  fprintf(outHandles[1],"----------\t----------\t----------\t----------\t----------\t----------\t----------\t----------\t----------\t-----\t------\t---\t---------------\t------------\t--------\n");
#endif
  for (plateIndex = 0; plateIndex < plateCount; plateIndex++) {
    pPlateEntry = &pPlatesTable[plateIndex];
    numSpatialBinCount[pPlateEntry->numSpatialBins]++;
    if (pPlateEntry->numSpatialBins == MAX_SPATIAL_BINS) {
      vectorCount = 0;
      memset(vector,0,sizeof(vector));
      for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
        vector[vectorCount] = pPlateEntry->colorterm[spatial_bin];
        fprintf(outHandles[1],"%0.3f\t",vector[vectorCount]);
        vectorCount++;
      }
      if (CalcMedianAndRMS(MAX_SPATIAL_BINS,0,vector,&pPlateEntry->median9bin,&pPlateEntry->rms9bin,0,3.0,0) == MAX_SPATIAL_BINS) {
        fprintf(outHandles[1],"%s%05d\t%f\t%f\t%f\t%f\t%d\n",
                pPlateEntry->series,
                pPlateEntry->plateNumber,
                pPlateEntry->median9bin,
                pPlateEntry->rms9bin,
                pPlateEntry->dRightAscension,
                pPlateEntry->dDeclination,
                pPlateEntry->seriesId);
        nineBinRMSCount++;
      }

    }
  }
  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    printf("%5d plates for %d spatial bins\n",numSpatialBinCount[spatial_bin],spatial_bin);
  }
  printf("Wrote %d records to %s\n",nineBinRMSCount,outfiles[1]);

  for (fileIndex = 0; fileIndex < NUM_OUT_FILES; fileIndex++) {
    fclose(outHandles[fileIndex]);
  }

  mysql_close(pConnection);
  mysql_close(pPhotConnection);
  mysql_close(pStacksConnection);
}
