// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* find_drad5.c
 *
 * Get the median and maximum dradRMS2 as a function of spatial bin and series using the photometry "localbin" table
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   find_drad5.c pipelineutils.a -ltable -lutil  -lwcs -o find_drad5 -L/usr/lib64/mysql -L/usr/lib/mysql -lmysqlclient
 *
 * May  3, 2010 Edward J. Los - Original Version
 * May  7, 2010 Edward J. Los - Add dradMedian
 * Dec 20, 2010 Edward J. Los - Add a plate count and limit results to solution 0
 * Mar 28, 2012 Edward J. Los - Rewrite to include plates that did not generate any photometry data
 *
 * NOTE: the following bin counts add up to 2500 if all bins are accepted
 *    smoothingBinCount  1      2      3      4      5      6      7      8      9  plates percent
 *                      224    232    216    232    230    280    282    345    459      1  100.0
 * 
 * 
 */


#include <math.h>
#include "table.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <assert.h>
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#define MAX_FILENAME 256
#define MAX_BUFFER 256
/* #define ALL_BINS 1 */

typedef struct _localbinEntry {
  char series[MAX_SERIES_STRING];
  int seriesId;
  int plateNumber;
  int local_bin_index;
  int localbin_versionId;
  int photplates_versionId;
  int mosaicWidth;
  int mosaicHeight;
  int nx;
  int ny;
  int spatial_bin;
  double draMedian;
  double ddecMedian;
  double draMedian2;
  double ddecMedian2;
  double dradRMS2;
  double scale;

} LOCALBINENTRY,*PLOCALBINENTRY;



typedef struct _plateEntry {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int naxis1;  /* Mosaic width in pixels */
  int naxis2; /* Mosaic height in pixels */
  int WCSSource;
  int FitWCS;
  int seriesId;
} PLATEENTRY,*PPLATEENTRY;


int main(int argc,char *argv[])
{
  char *argstr;
  char qualifier[MAX_BUFFER];
  char catalogString[MAX_BUFFER];
  char outputname[MAX_BUFFER];
  FILE *outputHandle = NULL;
  int errorFlag = 0;
  char cmdchar;
  int verbose = 0;
  int catalogNumber = 0;
  char *mysqlhost;
  char *username;
  char *password;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;

  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  char *photfilebase;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  int magnitudeFileFlag = 0;
  int solutionNumberFix = 0;
  char solutionNumberString[20];
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  char tempString[MAX_QUERY_STRING];
  int numRows = 0;
  int maximumPlates = 0;
  int totalPlates = 0;
  int goodRows = 0;
  int curRow = 0;
  PLOCALBINENTRY local_bin_table = NULL;
  PLOCALBINENTRY tmp_local_bin_table;
  int allocRows = 0;
  int totalRows = 0;
  int totalRowsRead = 0;
  PLOCALBINENTRY pCurEntry = NULL;
  int staleVersionIdCount = 0;
  int highDradCount = 0;
  int missingFieldCount = 0;
  int badSpatialBinCount = 0;
  int smoothingBinCount[MAX_SERIES][MAX_SPATIAL_BINS+1];
  int seriesCount[MAX_SERIES];       /* Number of entries in a series */
  int plateCount[MAX_SERIES];        /* Number of plates in a series with valid data */
  char *seriesPlateList[MAX_SERIES]; /* Bitmap of valid plate numbers */
  char *pPlateList;
  int curBins;

  int seriesIndex;
  int spatial_bin;
  int maxVectorCount = 0;
  int ix;
  int iy;
  int X_IMAGE;
  int Y_IMAGE;
  double edgeDist;
  double *vector;
  int vectorIndex;
  double median;
  double rms;
  time_t startTime;
  time_t curTime;
  char selectedSeries[MAX_SERIES_STRING];
  int selectedPlateNumber;
  PPLATEENTRY pPlateTable = NULL;
  PPLATEENTRY pPlateEntry;
  PLATEENTRY tempPlateEntry;
  PPLATEENTRY pNewPlateEntry = &tempPlateEntry;
  PLATEENTRY bestPlateEntry;
  PPLATEENTRY pBestPlateEntry = &bestPlateEntry;
  int plateIndex;
  int plotVector[X_DMAGBINS_NORMAL];
  int plotIndex;
  int plotCount = 0;

#ifdef ALL_BINS
  printf("ERROR: ALL_BINS is set!\n");
#endif /* ALL_BINS */
  outputname[0] = 0;
  time(&startTime);
  memset(smoothingBinCount,0,sizeof(smoothingBinCount));
  memset(seriesCount,0,sizeof(seriesCount));
  memset(plateCount,0,sizeof(plateCount));
  memset(seriesPlateList,0,sizeof(seriesPlateList));

  SetQueryCount(0);
  qualifier[0] = 0;
  catalogString[0] = 0;
  selectedSeries[0] = 0;
  selectedPlateNumber = -1;
  
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
        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;
          
        case 'p': /* plot good bins as a function of longer axis */
        case 'P':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outputname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;
          
        case 'n': /* selected plate number */
        case 'N':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&selectedPlateNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the selectedPlateNumber %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 's': /* Selected series */
        case 'S':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(selectedSeries,*++argv,MAX_SERIES_STRING);
            selectedSeries[MAX_SERIES_STRING-1] = 0;
          }
          break;
            


        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            
            if (catalogNumber < 0) {
              fprintf(stderr,"ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              sprintf(catalogString,"%d",catalogNumber);
              sprintf(qualifier,"_%s",catalogText[catalogNumber]);
            }
          }
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
    fprintf(stderr,"DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }


  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    fprintf(stderr,"ERROR: DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    fprintf(stderr,"ERROR: DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    fprintf(stderr,"ERROR: DASCH_PASSWORD is not defined\n");
    return(-1);
  }
                   
  mysql_init(pConnection);


  if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pConnection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
    }
    return(-1);
  }



  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    fprintf(stderr,"DASCH_PHOT_USERNAME is not defined\n");
    return(-1);
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    fprintf(stderr,"DASCH_PHOT_PASSWORD is not defined\n");
    return(-1);
  }
  mysql_init(pPhotConnection);
  
  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
    }
    return(-1);
  }
  InitSeriesTable(pConnection,pPhotConnection);

  if (GetPhotometryGlobal(pPhotConnection,pPhotGlobal) != 1) {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }
  if (pPhotGlobal->magnitudeFile == PHOT_MAGNITUDEFILE_YES) {
    magnitudeFileFlag = 1;
  }
  if (pPhotGlobal->solutionNumber == PHOT_SOLUTIONNUMBER_YES) {
    solutionNumberFix = 1;
    strcpy(solutionNumberString,"solutionNumber");
  } else {
    strcpy(solutionNumberString,"exposureNumber");
  }

  if (outputname[0] != 0) {
    outputHandle = fopen(outputname,"wt");
    if (outputHandle == NULL) {
      printf("ERROR: failed to open output file %s\n",outputname);
      errorFlag = 1;
    } else {
      fprintf(outputHandle,"iy\tcount\tseriesId\tplateNumber\tseries\n");
      fprintf(outputHandle,"--\t-----\t--------\t-----------\t------\n");
    }
  }

  if (errorFlag) {
    printf("Usage: find_drad5 options\n");
    printf("  options: -v verbose\n");
    printf("           -s series\n");
    printf("           -n plateNumber\n");
    printf("           -p <output summary>\n");
    printf("           -q <catalog>\n");
    printf(" The plotting feature totals the good smoothing bins along the\n");
    printf(" short axis (count) and plots them vs the long axis (iy)\n");
    return(-1);
  }

  printf("find_drad5 of %s %s for series '%s' and plateNumber %d and catalog %s\n",__DATE__,__TIME__,selectedSeries,selectedPlateNumber,catalogText[catalogNumber]);



  /* Start by building a list of plates to scan */
  sprintf(queryString,"SELECT series,plateNumber,mosaicNumber,naxis1,naxis2,WCSSource+0,FitWCS+0 FROM mosaics INNER JOIN scans using (series,plateNumber,scanNumber) where solutionNumber = 0 and aveADU is not NULL and FitWCS regexp('Selected') and ((mosaicComment IS NULL) or (not(mosaicComment regexp('deleted'))))");
  if (selectedSeries[0] != 0) {
    sprintf(tempString," and series = '%s'",selectedSeries);
    strcat(queryString,tempString);
  }
  if (selectedPlateNumber > 0) {
    sprintf(tempString," and plateNumber = %d",selectedPlateNumber);
    strcat(queryString,tempString);
  }
  strcat(queryString," order by series,plateNumber,scanNumber,mosaicNumber");
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      maximumPlates = mysql_affected_rows(pConnection);
      pPlateTable = (PPLATEENTRY)calloc(maximumPlates,sizeof(PLATEENTRY));
      if (pPlateTable == NULL) {
        printf("ERROR: failed to allocate pPlateTable\n");
        exit(-1);
      }
      memset(pBestPlateEntry,0,sizeof(PLATEENTRY));
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        memset(pNewPlateEntry,0,sizeof(PLATEENTRY));
        if (sqlrow[0]) {
          strcpy(pNewPlateEntry->series,sqlrow[0]);
        } else {
          continue;
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%d",&pNewPlateEntry->plateNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for plateNumber",nvals);
            continue;
          }
        } else {
          continue;
        }
        pNewPlateEntry->seriesId = GetSeriesId(pNewPlateEntry->series,1);
        if (pNewPlateEntry->seriesId < 0) {
          printf("ERROR: no seriesId for %s%05d\n",pNewPlateEntry->series,pNewPlateEntry->plateNumber);
          continue;
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pNewPlateEntry->mosaicNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for mosaicNumber",nvals);
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pNewPlateEntry->naxis1);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for naxis1",nvals);
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&pNewPlateEntry->naxis2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for naxis2",nvals);
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%d",&pNewPlateEntry->WCSSource);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for WCSSource",nvals);
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%d",&pNewPlateEntry->FitWCS);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for FitWCS",nvals);
            continue;
          }
        } else {
          continue;
        }
        if (pBestPlateEntry->series[0] == 0) {
          memcpy(pBestPlateEntry,pNewPlateEntry,sizeof(PLATEENTRY));
        } else {
          if ((pBestPlateEntry->plateNumber != pNewPlateEntry->plateNumber) ||
              (strcmp(pBestPlateEntry->series,pNewPlateEntry->series) != 0)) {
            pPlateEntry = &pPlateTable[totalPlates];
            memcpy(pPlateEntry,pBestPlateEntry,sizeof(PLATEENTRY));
            totalPlates++;
            memcpy(pBestPlateEntry,pNewPlateEntry,sizeof(PLATEENTRY));
          } else {
            if ((pNewPlateEntry->FitWCS & (FITWCS_SELECTED|FITWCS_COMPLETED)) >= (pBestPlateEntry->FitWCS & (FITWCS_SELECTED|FITWCS_COMPLETED))) {
              memcpy(pBestPlateEntry,pNewPlateEntry,sizeof(PLATEENTRY));
            }
          }

        }
      }
      pPlateEntry = &pPlateTable[totalPlates];
      memcpy(pPlateEntry,pBestPlateEntry,sizeof(PLATEENTRY));
      totalPlates++;

      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in find_drad5.c line %d\n",__LINE__);
    }
  }
  if (verbose) {
    printf("Found %d plates from a maximum of %d entries\n",totalPlates,maximumPlates);
  }

  for (plateIndex = 0; plateIndex < totalPlates; plateIndex++) {
    pPlateEntry = &pPlateTable[plateIndex];
    memset(plotVector,0,sizeof(plotVector));
    plotCount = 0;
#ifdef ALL_BINS 
    sprintf(queryString,"SELECT seriesId,plateNumber,local_bin_index,localbin%s.versionId,photplates.versionId%s,mosaicWidth,mosaicHeight,draMedian2,ddecMedian2,dradRMS2,scale,draMedian,ddecMedian FROM localbin%s INNER JOIN photplates USING (seriesId,plateNumber) INNER JOIN photseries USING (seriesId) where solutionNumber = 0 and seriesId = %d and plateNumber = %d;",catalogString,catalogString,catalogString,pPlateEntry->seriesId,pPlateEntry->plateNumber);
#else /* ALL_BINS */
    sprintf(queryString,"SELECT seriesId,plateNumber,local_bin_index,localbin%s.versionId,photplates.versionId%s,mosaicWidth,mosaicHeight,draMedian2,ddecMedian2,dradRMS2,scale,draMedian,ddecMedian FROM localbin%s INNER JOIN photplates USING (seriesId,plateNumber) INNER JOIN photseries USING (seriesId) where rejectFlag = 0 and solutionNumber = 0 and seriesId = %d and plateNumber = %d;",catalogString,catalogString,catalogString,pPlateEntry->seriesId,pPlateEntry->plateNumber);
#endif /* ALL_BINS */
    if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
      printf("ERROR: queryString exceeded in GetPhotLocalBin by %d\n",strlen(queryString));
      exit(-1);
    }
    res = ExecuteQuery(pPhotConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pPhotConnection);
      if (res_ptr) {
        numRows = mysql_affected_rows(pPhotConnection);
        if ((numRows+totalRows+10) > allocRows) {
          allocRows += (TOTAL_DMAGBINS_NORMAL * 1000);
          tmp_local_bin_table = (PLOCALBINENTRY)realloc(local_bin_table,allocRows * sizeof(LOCALBINENTRY));
          if (tmp_local_bin_table == NULL) {
            printf("ERROR: failed to reallocate local_bin_table of size %d\n",allocRows);
            exit(-1);
          }
          local_bin_table = tmp_local_bin_table;
          tmp_local_bin_table = NULL;        
        }
#if 0
        if (verbose && (numRows > 0)) {
          printf("Found %d rows in the localbin table\n",numRows);
        }
#endif
        totalRowsRead += numRows;
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          pCurEntry = &local_bin_table[totalRows];
          memset(pCurEntry,0,sizeof(LOCALBINENTRY));

          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&pCurEntry->seriesId);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for seriesId in GetPhotLocalBin\n");
              continue;
            } else {
              strcpy(pCurEntry->series,GetSeriesString(pCurEntry->seriesId,1));
            }
          } else {
            printf("ERROR: NULL seriesId in GetPhotLocalBin\n");
            continue;
          }



          if (sqlrow[1]) {
        
            nvals = sscanf(sqlrow[1],"%d",&pCurEntry->plateNumber);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for plateNumber\n");
              continue;
            }
          } else {
            continue;
          }
          if ((pCurEntry->plateNumber <= 0) ||
              (pCurEntry->plateNumber >= MAX_PLATE_NUMBER)) {
            printf("ERROR: plateNumber is %d for series %s\n",pCurEntry->plateNumber,pCurEntry->series);
            exit(-1);
          }

          if (sqlrow[2]) {
            nvals = sscanf(sqlrow[2],"%d",&pCurEntry->local_bin_index);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for local_bin_index\n");
              continue;
            }
          } else {
            continue;
          }

          if (sqlrow[3]) {
            nvals = sscanf(sqlrow[3],"%d",&pCurEntry->localbin_versionId);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for localbin_versionId\n");
              continue;
            }
          } else {
            continue;
          }

          if (sqlrow[4]) {
            nvals = sscanf(sqlrow[4],"%d",&pCurEntry->photplates_versionId);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for photplates_versionId\n");
              continue;
            }
          } else {
            continue;
          }

          if (sqlrow[5]) {
            nvals = sscanf(sqlrow[5],"%d",&pCurEntry->mosaicWidth);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for mosaicWidth\n");
              continue;
            }
          } else {
            continue;
          }


          if (sqlrow[6]) {
            nvals = sscanf(sqlrow[6],"%d",&pCurEntry->mosaicHeight);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for mosaicHeight\n");
              continue;
            }
          } else {
            continue;
          }

          pCurEntry->nx = XDmagBins(pCurEntry->mosaicWidth,pCurEntry->mosaicHeight);
          pCurEntry->ny = YDmagBins(pCurEntry->mosaicWidth,pCurEntry->mosaicHeight);

          if (sqlrow[7]) {
            nvals = sscanf(sqlrow[7],"%lf",&pCurEntry->draMedian2);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for draMedian2\n");
              continue;
            }
          } else {
            continue;
          }

          if (sqlrow[8]) {
            nvals = sscanf(sqlrow[8],"%lf",&pCurEntry->ddecMedian2);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for ddecMedian2\n");
              continue;
            }
          } else {
            continue;
          }



          if (sqlrow[9]) {
            nvals = sscanf(sqlrow[9],"%lf",&pCurEntry->dradRMS2);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for dradRMS2\n");
              continue;
            }
          } else {
            continue;
          }
          if (sqlrow[10]) {
            nvals = sscanf(sqlrow[10],"%lf",&pCurEntry->scale);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for scale\n");
              continue;
            }
          } else {
            continue;
          }



          if (sqlrow[11]) {
            nvals = sscanf(sqlrow[11],"%lf",&pCurEntry->draMedian);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for draMedian\n");
              continue;
            }
          } else {
            continue;
          }

          if (sqlrow[12]) {
            nvals = sscanf(sqlrow[12],"%lf",&pCurEntry->ddecMedian);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for ddecMedian\n");
              continue;
            }
          } else {
            continue;
          }




          if (pCurEntry->localbin_versionId != pCurEntry->photplates_versionId) {
            staleVersionIdCount++;
            continue;
          }
#if 0
          if (pCurEntry->seriesId == 31) {
            printf("Plate %5s%05d\n",pCurEntry->series,pCurEntry->plateNumber);
          }
#endif
#ifndef ALL_BINS
          if ((pCurEntry->draMedian > 90.0) ||
              (pCurEntry->ddecMedian > 90.0) ||
              (pCurEntry->draMedian2 > 90.0) ||
              (pCurEntry->ddecMedian2 > 90.0) ||
              (pCurEntry->dradRMS2 > 90.0)) {
            highDradCount++;
            continue;
          }
#endif /* ALL_BINS */
          if (pCurEntry->seriesId >= MAX_SERIES) {
            printf("ERROR: seriesId %d exceeds maximum\n",pCurEntry->seriesId);
            exit(-1);
          }
          /* Calculate the spatial bin of the center of the local bin */
        
          ix = pCurEntry->local_bin_index % pCurEntry->nx;
          iy = pCurEntry->local_bin_index / pCurEntry->ny;
          X_IMAGE = 1.0 * pCurEntry->mosaicWidth  * (0.5 + ix)/(1.0*pCurEntry->nx);
          Y_IMAGE = 1.0 * pCurEntry->mosaicHeight * (0.5 + iy)/(1.0*pCurEntry->ny);
          if (pCurEntry->mosaicHeight > pCurEntry->mosaicWidth) {
            plotVector[iy]++;
            plotCount++;
          } else {
            plotVector[ix]++;
            plotCount++;
          }
          pCurEntry->spatial_bin = CalculateBin(pCurEntry->mosaicWidth,pCurEntry->mosaicHeight,X_IMAGE,Y_IMAGE,&edgeDist);
          totalRows++;
          seriesCount[pCurEntry->seriesId]++;
          if (seriesPlateList[pCurEntry->seriesId] == NULL) {
            seriesPlateList[pCurEntry->seriesId] = (char *)calloc(MAX_PLATE_NUMBER,sizeof(char));
            if (seriesPlateList[pCurEntry->seriesId] == NULL) {
              printf("ERROR: Failed to allocate seriesPlateList for series %s\n",pCurEntry->series);
              exit(-1);
            }
          }
          pPlateList = seriesPlateList[pCurEntry->seriesId];
          if (pPlateList[pCurEntry->plateNumber] == 0) {
            pPlateList[pCurEntry->plateNumber] = 1;
            plateCount[pCurEntry->seriesId]++;

#if 0
            if (pCurEntry->seriesId == 31) {
              printf("Plate Count %5s%05d\n",pCurEntry->series,pCurEntry->plateNumber);
            }
#endif

          }
	  
          smoothingBinCount[pCurEntry->seriesId][pCurEntry->spatial_bin]++;
          if (smoothingBinCount[pCurEntry->seriesId][pCurEntry->spatial_bin] > maxVectorCount) {
            maxVectorCount = smoothingBinCount[pCurEntry->seriesId][pCurEntry->spatial_bin];
          }
          if (pCurEntry->spatial_bin == 0) {
            badSpatialBinCount++;
          }
        }
        if ((outputHandle != NULL) && (plotCount > 0)) {
          for (plotIndex = 0; plotIndex < pCurEntry->ny; plotIndex++) {
            fprintf(outputHandle,"%d\t%d\t%d\t%d\t%s\n",plotIndex,plotVector[plotIndex],pCurEntry->seriesId,pCurEntry->plateNumber,pCurEntry->series);
          }
          
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: mysql_store_result failed in find_drad5.c line %d\n",__LINE__);
      }
    }
  }
  goodRows = totalRows;
  missingFieldCount = totalRowsRead-goodRows-staleVersionIdCount -highDradCount;
  vector = (double*)calloc(maxVectorCount+1,sizeof(double));
  if (vector == NULL) {
    printf("ERROR: Failed to allocate vector of size %d\n",maxVectorCount);
    exit(-1);
  }
  for (seriesIndex = 0; seriesIndex < MAX_SERIES; seriesIndex++) {
    if (seriesCount[seriesIndex] == 0) {
      continue;
    }
    printf("Series %-6s Plate Count %5d\n",GetSeriesString(seriesIndex,1),plateCount[seriesIndex]);
  }

  printf("\nsmoothingBinCount   1       2       3       4       5       6       7       8       9  plates percent\n");
  for (seriesIndex = 0; seriesIndex < MAX_SERIES; seriesIndex++) {
    curBins = 0;
    if (seriesCount[seriesIndex] == 0) {
      continue;
    }
    printf("Series %-6s",GetSeriesString(seriesIndex,1));
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      curBins += smoothingBinCount[seriesIndex][spatial_bin];
      printf(" %7d",smoothingBinCount[seriesIndex][spatial_bin]);
    }
    printf("  %5d    %4.1f\n",plateCount[seriesIndex],(100.0*curBins)/(plateCount[seriesIndex]*TOTAL_DMAGBINS_NORMAL));
	
  }
  /* Print the initial spatial bin adjustment */

  printf("\ndradMedian  (pix) 1     2     3     4     5     6     7     8     9\n");
  for (seriesIndex = 0; seriesIndex < MAX_SERIES; seriesIndex++) {
    if (seriesCount[seriesIndex] == 0) {
      continue;
    }
    printf("Series %-6s",GetSeriesString(seriesIndex,1));
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      vectorIndex = 0;
      for (curRow = 0; curRow < numRows; curRow++) {
        pCurEntry = &local_bin_table[curRow];
        if ((pCurEntry->seriesId != seriesIndex) ||
            (pCurEntry->spatial_bin != spatial_bin)) {
          continue;
        }
        vector[vectorIndex] = sqrt(sqr(pCurEntry->draMedian)+sqr(pCurEntry->ddecMedian))/pCurEntry->scale;
        vectorIndex++;
        if (vectorIndex > maxVectorCount) {
          printf("ERROR: vectorIndex %d exceeds maxVectorCount %d\n",vectorIndex,maxVectorCount);
          exit(-1);
        }
      }
      CalcMedianAndRMS(vectorIndex,0,vector,&median,&rms,0,3.0,0);
      printf(" %5.1f",median);
    }
    printf("\n");
	
  }

  /* Print the final spatial bin adjustment */

  printf("\ndradMedian2 (pix) 1     2     3     4     5     6     7     8     9\n");
  for (seriesIndex = 0; seriesIndex < MAX_SERIES; seriesIndex++) {
    if (seriesCount[seriesIndex] == 0) {
      continue;
    }
    printf("Series %-6s",GetSeriesString(seriesIndex,1));
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      vectorIndex = 0;
      for (curRow = 0; curRow < totalRows; curRow++) {
        pCurEntry = &local_bin_table[curRow];
        if ((pCurEntry->seriesId != seriesIndex) ||
            (pCurEntry->spatial_bin != spatial_bin)) {
          continue;
        }
        vector[vectorIndex] = sqrt(sqr(pCurEntry->draMedian2)+sqr(pCurEntry->ddecMedian2))/pCurEntry->scale;
        vectorIndex++;
        if (vectorIndex > maxVectorCount) {
          printf("ERROR: vectorIndex %d exceeds maxVectorCount %d\n",vectorIndex,maxVectorCount);
          exit(-1);
        }
      }
      CalcMedianAndRMS(vectorIndex,0,vector,&median,&rms,0,3.0,0);
      printf(" %5.1f",median);
    }
    printf("\n");
	
  }

  /* Print the maximum image error compared with the catalog error */
  printf("\ndradRMS2 (pixels) 1     2     3     4     5     6     7     8     9\n");
  for (seriesIndex = 0; seriesIndex < MAX_SERIES; seriesIndex++) {
    if (seriesCount[seriesIndex] == 0) {
      continue;
    }
    printf("Series %-6s",GetSeriesString(seriesIndex,1));
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      vectorIndex = 0;
      for (curRow = 0; curRow < totalRows; curRow++) {
        pCurEntry = &local_bin_table[curRow];
        if ((pCurEntry->seriesId != seriesIndex) ||
            (pCurEntry->spatial_bin != spatial_bin)) {
          continue;
        }
        vector[vectorIndex] = pCurEntry->dradRMS2/pCurEntry->scale;
        vectorIndex++;
        if (vectorIndex > maxVectorCount) {
          printf("ERROR: vectorIndex %d exceeds maxVectorCount %d\n",vectorIndex,maxVectorCount);
          exit(-1);
        }
      }
      CalcMedianAndRMS(vectorIndex,0,vector,&median,&rms,0,3.0,0);
      printf(" %5.1f",median);
    }
    printf("\n");
	
  }


  if (local_bin_table != NULL) {
    free(local_bin_table);
  }
  if (vector != NULL) {
    free(vector);
  }

  printf("\nExecution complete numRows %d goodRows %d missingFieldCount %d staleVersionIdCount %d highDradCount %d badSpatialBinCount %d\n",totalRowsRead,goodRows,missingFieldCount,staleVersionIdCount,highDradCount,badSpatialBinCount);
  time(&curTime);
  curTime -= startTime;


  printf("maxVectorCount %d seconds %d\n",maxVectorCount,curTime);
  mysql_close(pPhotConnection);
  mysql_close(pConnection);

  for (seriesIndex = 0; seriesIndex < MAX_SERIES; seriesIndex++) {
    if (seriesPlateList[seriesIndex] != NULL) {
      free(seriesPlateList[seriesIndex]);
    }
  }
  if (outputHandle != NULL) {
    free(outputHandle);
  }

  return(EXIT_SUCCESS);
}
