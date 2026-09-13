// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* find_drad3.c
 * This routine takes a concatenation of all of the match_<plate>_drad.db files and generates a table of median_drad as a function of series and 
 * spatial bin.
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I/usr/include/mysql  -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio -lmysqlclient find_drad3.c pipelineutils.a -ltable -lutil  -lwcs -o find_drad3
 *
 * /dasch/Pipeline/find_drad3 -i /dasch/backup/2008_06_26/drad.db
 *
 * Jul 1, 2008 Edward J. Los - Initial version
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 *
 */


#include <math.h>
#include <time.h>
#include "mysql.h"
#include "table.h"
#include "pipelineutils.h"
#include "fitsio.h"
#include "longnam.h"
#define MAX_BUFFER 100
#define BIN_FACTOR 64
#define MAX_PLATE_SIZE 25
#define MAX_PLATES 1000000

typedef struct _image {
  char plate[MAX_PLATE_SIZE];
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int spatial_bin;
  double median_drad;
} PLATEDRAD,*PPLATEDRAD;

typedef struct _series {
  char series[MAX_SERIES_STRING];
  double ave_median_drad[MAX_SPATIAL_BINS+1];
  double median_median_drad[MAX_SPATIAL_BINS+1];
  int median_drad_count[MAX_SPATIAL_BINS+1];
  double plateScale;
  double printFlag;
} SERIES,*PSERIES;

void GetScales(PSERIES pSeriesTable,int numSeries) {
  int res;
  MYSQL my_connection;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char *mysqlhost;
  char *username;
  char *password;
  char *series;
  char queryString[MAX_QUERY_STRING];
  double nominalPlateScale = 0.0;
  double fittedPlateScale = 0.0;
  int nvals;
  int index;
  PSERIES pSeries;


  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    fprintf(stderr,"DASCH_MYSQLHOST is not defined\n");
    exit(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    fprintf(stderr,"DASCH_USERNAME is not defined\n");
    exit(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    fprintf(stderr,"DASCH_PASSWORD is not defined\n");
    exit(-1);
  }
                   

  mysql_init(&my_connection);


  if (mysql_real_connect(&my_connection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
#if 0
    fprintf(stderr,"Connection success\n");
#endif
    sprintf(queryString,"SELECT series,nominalPlateScale,fittedPlateScale FROM series;");
    if (strlen(queryString) > MAX_QUERY_STRING) {
      fprintf(stderr,"ERROR: MAX_QUERY_STRING exceeded %d\n",strlen(queryString));
      exit(-1);
    }


    res = mysql_query(&my_connection,queryString);
    if (!res) {
#if 0
      fprintf(stderr,"insert ID: %lu rows res %d\n",(unsigned long)mysql_affected_rows(&my_connection),res);
#endif
      res_ptr = mysql_store_result(&my_connection);
      if (res_ptr) {
        MYSQL_FIELD *field_ptr;
#if 0
        fprintf(stderr,"Retrieved %lu rows\n",(unsigned long)mysql_num_rows(res_ptr));
#endif
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          nominalPlateScale = 0.0;
          fittedPlateScale = 0.0;
          
          if (sqlrow[0] == NULL) {
            continue;
          }
          series = sqlrow[0];
          if (sqlrow[1] != NULL) {
            nvals = sscanf(sqlrow[1],"%lf",&nominalPlateScale);
            if (nvals != 1) {
              nominalPlateScale = 0.0;
            }
          }
          if (sqlrow[2] != NULL) {
            nvals = sscanf(sqlrow[2],"%lf",&fittedPlateScale);
            if (nvals != 1) {
              fittedPlateScale = 0.0;
            }
          }
          if (fittedPlateScale == 0.0) {
            fittedPlateScale = nominalPlateScale;
          }
          if (fittedPlateScale == 0.0) {
            continue;
          }
          for (index = 0; index < numSeries; index++) {
            pSeries = &pSeriesTable[index];
            if (strcmp(series,pSeries->series) == 0) {
              pSeries->plateScale = fittedPlateScale * NOMINAL_MM_PER_PIXEL;
              printf("Series %5s, arcsec/pixel %7.2f\n",series,pSeries->plateScale);
              break;
            }
          


          }

        }

        mysql_free_result(res_ptr); 
      }
    } else {
      fprintf(stderr,"ERROR: Select error %d: %s res %d\n",mysql_errno(&my_connection),mysql_error(&my_connection),res);
      exit(-1);
    }



    mysql_close(&my_connection);
  } else {
    fprintf(stderr,"ERROR: Connection failed\n");
    if (mysql_errno(&my_connection)) {
      fprintf(stderr,"ERROR: Connection error %d: %s\n",mysql_errno(&my_connection),mysql_error(&my_connection));
    }
    exit(-1);
  }


}


int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char infile[MAX_BUFFER];
  char cmdchar;
  int errorFlag = 0;
  char *charPtr;
  char *basePtr;
  char charVal;

  File input_handle = NULL;
  char input_name[MAX_BUFFER];
  TableHead input_header = NULL;
  PPLATEDRAD input_table = NULL;
  size_t input_nrecs = 0;
  int input_index;
  PPLATEDRAD pInput;
  PSERIES pSeriesTable = NULL;
  PSERIES pSeries;
  int numSeries = 0;
  int seriesIndex;
  double *vector;
  int vector_count;
  int verbose = 0;
  double rms;
  int index;
  time_t startTime;
  time_t curTime;

  time(&startTime);

  input_name[0] = 0;

  pSeriesTable = (PSERIES)calloc(MAX_SERIES,sizeof(SERIES));
  vector = (double *)calloc(MAX_PLATES,sizeof(double));
  if (pSeriesTable == NULL) {
    printf("ERROR: Failed to allocate the series table\n");
  }
  if (vector == NULL) {
    printf("ERROR: Failed to allocate vector\n");
  }
  


  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      fprintf(stderr,"ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'i': /* input file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(input_name,*++argv,MAX_BUFFER-2);
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
          fprintf(stderr,"ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */

  if (input_name[0] == 0) {
    fprintf(stderr,"ERROR: No input filename was specified\n");
    errorFlag = 1;
  }




  /* Open the input file */
  input_handle = Open(input_name,"r");
  if (input_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the input file %s\n",input_name);
  } else {
    fprintf(stderr,"Found input file %s\n",input_name);
  }



  if (errorFlag) {
    fprintf(stderr," Usage: find_drad3 -i <input file> [-v]\n");

    return(-1);
  }

  /* Now read in the input file */

  input_header = table_header(input_handle,TABLE_PARSE);
  if (input_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",input_name);
    return(-1);
  }

  input_table = table_loadva(input_handle,
                             &input_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(PLATEDRAD),
                             &input_nrecs,
                             TblInt,"spatial_bin",TblOff(PPLATEDRAD,spatial_bin),
                             TblDbl,"median_drad",TblOff(PPLATEDRAD,median_drad),
                             TblBuf,"plate"    ,TblOff(PPLATEDRAD,plate),MAX_PLATE_SIZE,
                             0,"end",0);
  if (input_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",input_name);
    return(-1);
  }
  fprintf(stderr,"read %d records for %s\n",input_nrecs,input_name);

  /* Now go through the input data and populate the arrays */
  for (index = 0; index < input_nrecs; index++) {
    pInput = &input_table[index];
    if (ParseFilename2(pInput->plate,pInput->series,&pInput->plateNumber) != 1) {
      printf("ERROR: Failed to parse %s\n",pInput->plate);
      continue;
    }
    if ((pInput->spatial_bin < 0) ||
        (pInput->spatial_bin > MAX_SPATIAL_BINS)) {
      printf("ERROR: Spatial bin %d is illegal for %s\n",pInput->spatial_bin,pInput->plate);
      continue;
    }

    for (seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
      pSeries = &pSeriesTable[seriesIndex];
      if (strcmp(pInput->series,pSeries->series) == 0) {
        break;
      }

    }
    if (seriesIndex == numSeries) {
      if (numSeries >= MAX_SERIES) {
        printf("ERROR: MAX_SERIES %d is exceeded\n",MAX_SERIES);
        return(-1);
      }
      pSeries = &pSeriesTable[numSeries];
      strcpy(pSeries->series,pInput->series);
      numSeries++;

    }
    pSeries->ave_median_drad[pInput->spatial_bin] += pInput->median_drad;
    pSeries->median_drad_count[pInput->spatial_bin]++;
  }
  GetScales(pSeriesTable,numSeries);


  /* Compute clipped medians */
  for (seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
    pSeries = &pSeriesTable[seriesIndex];
    for (index = 1; index <= MAX_SPATIAL_BINS; index++) {
      if (pSeries->median_drad_count[index] != pSeries->median_drad_count[1]) {
        printf("ERROR: median_drad count does not agree for series %s\n",pSeries->series);
        exit(-1);
      }
      vector_count = 0;
      for (input_index = 0; input_index < input_nrecs; input_index++) {
        pInput = &input_table[input_index];
        if ((pInput->spatial_bin == index) &&
            (strcmp(pInput->series,pSeries->series) == 0) &&
            (pInput->median_drad != 0)) {
          vector[vector_count] = pInput->median_drad;
          vector_count++;

        }
      }
      if ((vector_count != pSeries->median_drad_count[index]) &&
          (pSeries->printFlag == 0)) {
        printf("Vector count %d disagrees with median_drad_count %d for series %s\n",vector_count,pSeries->median_drad_count[index],pSeries->series);
        pSeries->printFlag = 1;
      }
      CalcMedianAndRMS(vector_count,0,vector,&pSeries->median_median_drad[index],&rms,1,3.0,0);
      
    }
  }

  /* Display our results */
  if (verbose) {
    printf("ave, median Drad for each bin in arcsec\n");
  } else {
    printf("median Drad for each bin in arcsec\n");
  }
  for (seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
    pSeries = &pSeriesTable[seriesIndex];
    if (verbose) {
      printf("Series %3d %6s Plates: %3d |",seriesIndex,pSeries->series,pSeries->median_drad_count[1]);
    } else {
      printf("Series %6s Plates: %3d |",pSeries->series,pSeries->median_drad_count[1]);      
    }
    for (index = 1; index <= MAX_SPATIAL_BINS; index++) {
      if (pSeries->median_drad_count[index] != 0) {
        pSeries->ave_median_drad[index] = pSeries->ave_median_drad[index]/(1.0 * pSeries->median_drad_count[index]);
      }
      if (verbose) {
        printf(" %4.1f %4.1f |",pSeries->ave_median_drad[index],pSeries->median_median_drad[index]);
      } else {
        printf(" %4.1f |",pSeries->median_median_drad[index]);
      }

    }
    printf("\n");
  }

  if (verbose) {
    printf("ave, median Drad for each bin in pixels\n");
  } else {
    printf("median Drad for each bin in pixels\n");
  }
  for (seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
    pSeries = &pSeriesTable[seriesIndex];
    if (pSeries->plateScale == 0.0) {
      continue;
    }
    if (verbose) {
      printf("Series %3d %6s Plates: %3d |",seriesIndex,pSeries->series,pSeries->median_drad_count[1]);
    } else {
      printf("Series %6s Plates: %3d |",pSeries->series,pSeries->median_drad_count[1]);      
    }
    for (index = 1; index <= MAX_SPATIAL_BINS; index++) {
      if (verbose) {
        printf(" %4.2f %4.2f |",pSeries->ave_median_drad[index]/pSeries->plateScale,
               pSeries->median_median_drad[index]/pSeries->plateScale);
      } else {
        printf(" %4.2f |",pSeries->median_median_drad[index]/pSeries->plateScale);

      }

    }
    printf("\n");
  }
  /* All done.  Clean up */


  if (input_table != NULL) {
    Free(input_table);
  }
  if (input_header != NULL) {
    table_hdrfree(input_header);
  }
  if (input_handle != NULL) {
    Close(input_handle);
  }


  
  if (pSeriesTable != NULL) {
    free(pSeriesTable);
  }

  if (vector != NULL) {
    free(vector);
  }


  time(&curTime);
  curTime -= startTime;
  fprintf(stdout,"maxseconds %d for %s\n",
          curTime,
          input_name);
    


  return(0);
}

