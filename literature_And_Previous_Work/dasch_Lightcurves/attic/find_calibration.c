// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* find_calibration.c
 *
 *  This program analyzes the calibration table of the photometry database
 *  looking for the maximum differences between adjacent spatial bins as 
 *  a function of magnitude 
 *
 * cc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio  -L/usr/lib/mysql  -l mysqlclient  find_calibration.c pipelineutils.a -ltable -lutil  -lwcs -pthread -o find_calibration
 *
 *  find_calibration -q kepler -o /dasch/backup/2011_02_03/find_calibration_kepler.db
 *
 * Feb  2, 2011 Edward J. Los - Original Version
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
#define MAX_GRID_ENTRY 4000
#define MAX_BUFFER 256
#define MIN_MAG 20

typedef struct _pcalpoint {
  float maggrid;
  float isogrid;
  float griderr;
  int flaggrid;
} CALPOINT,*PCALPOINT;

typedef struct _pcalbin {
  double limiting_mag;
  double max_bright_mag;
} CALBIN,*PCALBIN;

typedef struct _pcalibration {
  char series[MAX_SERIES_STRING]; /* Series */
  int seriesId; /* Series ID */
  int plateNumber;
  int solutionNumber;
  int errorFlag;
  int maxmagflag[MAX_SPATIAL_BINS+1][MIN_MAG];
  float maxmagdiff[MAX_SPATIAL_BINS+1][MIN_MAG];
  CALBIN calbins[MAX_SPATIAL_BINS+1];

} CALIBRATION,*PCALIBRATION; 

/* Here find the mangnitude point for the next calibration curve and the same isogrid value */
int InterpolateIso(PCALPOINT pCalResult,PCALPOINT pCalPoint,int numRows, PCALPOINT calpointTable)
{
  PCALPOINT pNextCalPoint = NULL;
  PCALPOINT pCurCalPoint;
  int curRow;
  memset(pCalResult,0,sizeof(CALPOINT));


  for (curRow = 0; curRow < numRows; curRow++) {
    pCurCalPoint = &calpointTable[curRow];
    pNextCalPoint = &calpointTable[curRow+1];
    if ((pCalPoint->isogrid > pCurCalPoint->isogrid) &&
        (pCalPoint->isogrid <= pNextCalPoint->isogrid)) {
      if ((pCurCalPoint->flaggrid != 1) ||
          (pNextCalPoint->flaggrid != 1) ||
          (pCurCalPoint->griderr > MAX_ISO_RMS) ||
          (pNextCalPoint->griderr > MAX_ISO_RMS)) {
        return(-1);
      }
      memcpy(pCalResult,pCurCalPoint,sizeof(CALPOINT));
      if (pNextCalPoint->griderr > pCurCalPoint->griderr) {
        pCalResult->griderr = pNextCalPoint->griderr;
      }
      pCalResult->isogrid = pCalPoint->isogrid;
      if ((pNextCalPoint->isogrid-pCurCalPoint->isogrid) == 0) {
        return(-1);
      } else { 
        pCalResult->maggrid = pCurCalPoint->maggrid + ((pNextCalPoint->maggrid-pCurCalPoint->maggrid)*(pCalPoint->isogrid-pCurCalPoint->isogrid)/(pNextCalPoint->isogrid-pCurCalPoint->isogrid));
      }

#if 0
      {
        int tmpRow;
        PCALPOINT pTmpCalPoint;
  
        for (tmpRow = 0; tmpRow < numRows; tmpRow++) {
          pTmpCalPoint = &calpointTable[tmpRow];
          printf("Row %5d maggrid %9.5f isogrid %9.5f griderr %9.5f flaggrid %3d\n",tmpRow,pTmpCalPoint->maggrid,pTmpCalPoint->isogrid,pTmpCalPoint->griderr,pTmpCalPoint->flaggrid);
          

        }
        exit(-1);
      }
#endif

      return(0);
    }

  }
  return(-1);
}

int main(int argc,char *argv[])
{
  char *mysqlphothost;
  char *catalogallDirectory;
  char *scratchDirectory;
  char *mysqlhost;
  char *username;
  char *password;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  char *photusername;
  char *photpassword;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  char queryString[MAX_QUERY_STRING];
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  int errorFlag = 0;
  char cmdchar;
  char qualifier[MAX_BUFFER];
  int verbose = 0;
  char *argstr;
  int res;
  int numRows = 0;
  int curRow = 0;
  int curPlate = 0;
  int numPlates = 0;
  PCALIBRATION pCalibrationTable = NULL;
  PCALIBRATION pCalibration;
  PCALBIN pCalBin;
  PCALBIN pNextCalBin;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int nvals;
  int solutionNumberFix = 0;
  PPHOT_SPATIAL_BIN spatial_bin_table = NULL;
  PPHOT_SPATIAL_BIN pSpatialBinEntry;
  int spatial_bin_nrecs = 0;
  int spatial_bin_index;
  int spatial_bin;
  
  int calibrationRows[MAX_SPATIAL_BINS+1];
  
  CALPOINT calpointTable[MAX_SPATIAL_BINS+1][MAX_GRID_ENTRY];
  PCALPOINT pCalPoint;
  CALPOINT calresult;
  PCALPOINT pCalResult = &calresult;
  int next_spatial_bin;
  int result;
  float curmagdiff;
  int magintval;
  char outfile[MAX_BUFFER];
  FILE *outHandle = NULL;

  qualifier[0] = 0;
  outfile[0] = 0;
  catalogString[0] = 0;

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
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
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
  if (outfile[0] == 0) {
    printf("ERROR: No output filename was specified\n");
    errorFlag = 1;
  }
  outHandle = fopen(outfile,"wt");
  if (outHandle == NULL) {
    fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
    errorFlag = 1;
  } else {
    if (verbose) {
      fprintf(stderr,"Output file %s\n",outfile);
    }
  }

  fprintf(outHandle,"seriesId\tplateNumber\tsolutionNumber\tspatial_bin\tmagcal_iso\tmaxmagdiff\n");
  fprintf(outHandle,"--------\t-----------\t--------------\t-----------\t----------\t----------\n");

  if (errorFlag) {
    printf("Usage: update_photometry options\n");
    printf("  options: -v verbose\n");
    printf("           -o <output file name>\n");
    printf("           -q input catalog and filename qualifer\n");

    return(-1);
  }
  printf("find_calibration of %s %s\n",
         __DATE__,__TIME__);



  /* Connect to the database */
  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    fprintf(stderr,"DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }
  catalogallDirectory = getenv("DASCH_CATALOGALL");
  if (catalogallDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_CATALOGALL is not defined\n");
    return(-1);
  }

  scratchDirectory = getenv("DASCH_SCRATCH");
  if (scratchDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_SCRATCH is not defined\n");
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
  if (pPhotGlobal->solutionNumber == PHOT_SOLUTIONNUMBER_YES) {
    solutionNumberFix = 1;
  } 
  /* Start by getting a list of plates */

#if 0
  sprintf(queryString,"SELECT seriesId,plateNumber FROM photplates WHERE versionId%s = %d;",catalogString,pPhotGlobal->currentVersion);
#else
  sprintf(queryString,"SELECT seriesId,plateNumber FROM photplates; ");
#endif
  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pPhotConnection);
      pCalibrationTable = (PCALIBRATION)calloc(numRows,sizeof(CALIBRATION));
      if (pCalibrationTable == NULL) {
        printf("ERROR: failed to allocate %d CALIBRATION rows\n",numRows);
        exit(-1);
      }
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        pCalibration = &pCalibrationTable[curRow];
        memset(pCalibration,0,sizeof(CALIBRATION));
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pCalibration->seriesId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for seriesId in GetMagnitudes\n");
            continue;
          } else {
            strcpy(pCalibration->series,GetSeriesString(pCalibration->seriesId,1));
          }
        } else {
          printf("ERROR: NULL seriesId in GetMagnitudes\n");
          continue;
        }




        if (sqlrow[1]) {
        
          nvals = sscanf(sqlrow[1],"%d",&pCalibration->plateNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for plateNumber\n");
            continue;
          }
        } else {
          pCalibration->plateNumber = 0;
        }
        curRow++;
      }

      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed line %d\n",__LINE__);
    }
  } else {
    return(-1);
  }
  if (verbose) {
    printf("Fetched %d valid plates from a total of %d\n",curRow,numRows);
  }
  numPlates = curRow;
 
  /* Process each plate in turn */
  for (curPlate = 0; curPlate < numPlates; curPlate++) {
    pCalibration = &pCalibrationTable[curPlate];
#if 0
    if ((pCalibration->plateNumber == 502) &&
        (strcmp(pCalibration->series,"ac") == 0)) {
      printf("At plate %s%05d\n",pCalibration->series,pCalibration->plateNumber);
    }

#endif
    
    /* Now read in the spatial bin table for this plate to get the limiting_mag and max_bright_mag */
    spatial_bin_nrecs = ReadPhotSpatialBin(pPhotConnection,&spatial_bin_table,pCalibration->series,pCalibration->plateNumber,pCalibration->solutionNumber,0,catalogString,solutionNumberFix);
    if (spatial_bin_nrecs == 0) {
      pCalibration->errorFlag = 1;
      continue;
    }
    for (spatial_bin_index = 0; spatial_bin_index < spatial_bin_nrecs; spatial_bin_index++) {
      pSpatialBinEntry = &spatial_bin_table[spatial_bin_index];
      spatial_bin = pSpatialBinEntry->spatial_bin;
      if ((spatial_bin > 0) && (spatial_bin <= MAX_SPATIAL_BINS)) {
        pCalBin = &pCalibration->calbins[spatial_bin];
        pCalBin->limiting_mag = pSpatialBinEntry->limiting_mag;
        pCalBin->max_bright_mag = pSpatialBinEntry->max_bright_mag;
      } else {
        printf("ERROR: Illegal spatial bin %d\n",spatial_bin);
        exit(-1);
      }
    }
    if (spatial_bin_table != NULL) {
      free(spatial_bin_table);
      spatial_bin_table = NULL;
    }

  }
  /* Now read in the calibration data */
  for (curPlate = 0; curPlate < numPlates; curPlate++) {
    pCalibration = &pCalibrationTable[curPlate];
    if (pCalibration->errorFlag != 0) {
      continue;
    }
#if 0
    if ((pCalibration->plateNumber == 502) &&
        (strcmp(pCalibration->series,"ac") == 0)) {
      printf("At plate %s%05d\n",pCalibration->series,pCalibration->plateNumber);
    }

#endif

    memset(calibrationRows,0,sizeof(calibrationRows));
    memset(calpointTable,0,sizeof(calpointTable));
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      curRow = 0;

      sprintf(queryString,"SELECT maggrid,isogrid,griderr,flaggrid FROM calibration%s where seriesID = %d and plateNumber = %d and spatial_bin = %d; ",catalogString,pCalibration->seriesId,pCalibration->plateNumber,spatial_bin);
      res = ExecuteQuery(pPhotConnection,queryString);
      if (!res) {
        res_ptr = mysql_store_result(pPhotConnection);
        if (res_ptr) {
          numRows = mysql_affected_rows(pPhotConnection);
          if (numRows > MAX_GRID_ENTRY) {
            printf("ERROR: numRows %d exceeds MAX_GRID_ENTRY %d for %s%05d %d\n",numRows,MAX_GRID_ENTRY,pCalibration->series,pCalibration->plateNumber,spatial_bin);
            exit(-1);
          }
          calibrationRows[spatial_bin] = numRows;


          while ((sqlrow = mysql_fetch_row(res_ptr))) {
            pCalPoint = &calpointTable[spatial_bin][curRow];
            memset(pCalPoint,0,sizeof(CALPOINT));
            if (sqlrow[0]) {
              nvals = sscanf(sqlrow[0],"%f",&pCalPoint->maggrid);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for for maggrid\n");
                continue;
              }
            } else {
              printf("ERROR: NULL maggrid\n");
              continue;
            }
            if (sqlrow[1]) {
              nvals = sscanf(sqlrow[1],"%f",&pCalPoint->isogrid);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for for isogrid\n");
                continue;
              }
            } else {
              printf("ERROR: NULL isogrid\n");
              continue;
            }
            if (sqlrow[2]) {
              nvals = sscanf(sqlrow[2],"%f",&pCalPoint->griderr);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for for griderr\n");
                continue;
              }
            } else {
              printf("ERROR: NULL griderr\n");
              continue;
            }
            if (sqlrow[3]) {
              nvals = sscanf(sqlrow[3],"%d",&pCalPoint->flaggrid);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for for flaggrid\n");
                continue;
              }
            } else {
              printf("ERROR: NULL flaggrid\n");
              continue;
            }
            curRow++;
            if (curRow > numRows) {
              printf("ERROR: curRow %d exceeds numRows %d\n",curRow,numRows);
              exit(-1);
            }
          }
          calibrationRows[spatial_bin] = curRow;
          if (verbose) {
            if (curRow != 0) {
              printf("Found %d rows for spatial_bin %d in plate %s%05d\n",curRow,spatial_bin,pCalibration->series,pCalibration->plateNumber);
            }
#if 0
            for (curRow = 0; curRow < calibrationRows[spatial_bin]; curRow++) {
              pCalPoint = &calpointTable[spatial_bin][curRow];
              printf("Row %5d maggrid %9.5f isogrid %9.5f griderr %9.5f flaggrid %3d\n",curRow,pCalPoint->maggrid,pCalPoint->isogrid,pCalPoint->griderr,pCalPoint->flaggrid);
            }
#endif
          }
          mysql_free_result(res_ptr);
        } else {
          printf("ERROR: mysql_store_result failed line %d\n",__LINE__);
        }
      } else {
        return(-1);
      }
      

    }
    for (spatial_bin = 1; spatial_bin < MAX_SPATIAL_BINS; spatial_bin++) {
      pCalBin = &pCalibration->calbins[spatial_bin];
#if 0
      printf("spatial_bin %d limiting_mag %9.5f max_bright_mag %9.5f\n",spatial_bin,pCalBin->limiting_mag,pCalBin->max_bright_mag);
#endif
      next_spatial_bin = spatial_bin+1;
      if ((calibrationRows[spatial_bin] == 0) ||
          (calibrationRows[next_spatial_bin] == 0)) {
        continue;
      }
      for (curRow = 0; curRow < calibrationRows[spatial_bin]; curRow++) {
        pCalPoint = &calpointTable[spatial_bin][curRow];
        /* Pick only good points in the first spatial bin */
        if ((pCalPoint->flaggrid != 1) ||
            (pCalPoint->griderr > MAX_ISO_RMS) ||
            (pCalPoint->maggrid < pCalBin->max_bright_mag) ||
            (pCalPoint->maggrid > (pCalBin->limiting_mag-MAX_LIMITING_MAG))) {
          continue;
        }
        result = InterpolateIso(pCalResult,pCalPoint,calibrationRows[next_spatial_bin],calpointTable[next_spatial_bin]);
        if (result == 0) {
          pNextCalBin = &pCalibration->calbins[next_spatial_bin];
          if ((pCalResult->maggrid < pNextCalBin->max_bright_mag) ||
              (pCalResult->maggrid > (pNextCalBin->limiting_mag-MAX_LIMITING_MAG))) {
            continue;
          }
          curmagdiff = fabsf(pCalResult->maggrid-pCalPoint->maggrid);

          magintval = pCalPoint->maggrid;
          if ((magintval < 0) || (magintval >= MIN_MAG)) {
            continue;
          }
#if 0
          printf("maggrid %9.5f %9.5f isogrid %9.5f %9.5f griderr %9.5f %9.5f flaggrid %3d %3d curmagdiff %9.5f magintval %2d\n",
                 pCalPoint->maggrid,
                 pCalResult->maggrid,
                 pCalPoint->isogrid,
                 pCalResult->isogrid,
                 pCalPoint->griderr,
                 pCalResult->griderr,
                 pCalPoint->flaggrid,
                 pCalResult->flaggrid,
                 curmagdiff,
                 magintval);
#endif       

          if (curmagdiff > pCalibration->maxmagdiff[spatial_bin][magintval]) {
            pCalibration->maxmagdiff[spatial_bin][magintval] = curmagdiff;
            pCalibration->maxmagflag[spatial_bin][magintval] = 1;
          }
        }


      }

    }
    for (spatial_bin = 1; spatial_bin < MAX_SPATIAL_BINS; spatial_bin++) {
      for (magintval = 0; magintval < MIN_MAG; magintval++) {
        if (pCalibration->maxmagflag[spatial_bin][magintval] != 0) {
          fprintf(outHandle,"%d\t%d\t%d\t%d\t%d\t%4.2f\n",pCalibration->seriesId,pCalibration->plateNumber,pCalibration->solutionNumber,spatial_bin,magintval,pCalibration->maxmagdiff[spatial_bin][magintval]);
        }
      }
    }
  }
  

  if (outHandle != NULL) {
    fclose(outHandle);
  }


  mysql_close(pPhotConnection);
  mysql_close(pConnection);
  if (pCalibrationTable != NULL) {
    free(pCalibrationTable);
  }

  return(EXIT_SUCCESS);
} 
