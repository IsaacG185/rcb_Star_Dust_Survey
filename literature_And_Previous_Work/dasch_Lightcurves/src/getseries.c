// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* getlocation.c
 *
 *  gcc -ggdb -O2  -I/usr/include/mysql -I /dasch/install/include -L /dasch/install/lib getseries.c pipelineutils.a -lmysqlclient -lwcs -o getseries
 *
 * This file takes the plate name the minimum mean and maximum plate scale
 *
 *  Other inputs: DASCH_USERNAME - MySQL user name
 *                DASCH_PASSWORD - MySQL password
 *
 *
 * May  2, 2007  Edward J. Los - Initial version
 * Oct  6, 2020  For the "na" dummy series, force the minimum scale to 10 arcsec/mm and the maxium to 1600 arcsec/mm to cover
 *                the entire range for these unknown plates.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <mysql.h>

#include <fitsio.h>
#include <longnam.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"
#include "pipelineutils.h"


#define MAX_COMMENT_STRING 132
#define MAX_CENTERSOURCE_STRING 10
#define MAX_RIGHTASCENSION_STRING 13
#define MAX_DECLINATION_STRING 16
#define MAX_DATE_STRING 25
#define MAX_FILENAME_LEN 256
#define A1_FACTOR 0.23666  /* From annular9.m */

int main(int argc,char *argv[])
{
  int res;
  MYSQL my_connection;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int binning = 1;
  int nvals;
  char filename[MAX_FILENAME_LEN];
  char queryString[MAX_QUERY_STRING];
  double nominalPlateScale = 0.0;
  double fittedPlateScale = 0.0;
  double factor = 1.1;
  double factor2;
  char *argstr;
  char cmdchar;
  int verboseFlag = 0;
  int totalCount = 0;
  int errorFlag = 0;
  int foundEntry = 0;

  /* Loop through the arguments */
  filename[0] = 0;
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the plate designation */
      if (strlen(filename) > 0) {
        errorFlag = 1;
        printf("ERROR: list %s is being overwritten by %s\n",filename,argstr);
      } else {
        strncpy(filename,argstr,MAX_FILENAME_LEN-1);
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'v': /* verbose flag */
        case 'V':
          verboseFlag = 1;
          break;
        case 'b':
        case 'B': /* Binning */
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&binning);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the number of bins %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;
        case 'f':
        case 'F': /* Ratio factor */
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&factor);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the factor %s\n",*argv);
              errorFlag = 1;
            } else {
              if (factor <= 1.0) {
                fprintf(stderr,"ERROR: factor %d must be greater than 1.0\n",factor);
                errorFlag = 1;
              }
            }
          }
          break;


        default:
          printf("* illegal command -%c-",cmdchar);
          errorFlag = 1;

        }

      }

    }
  }
  if (filename[0] == 0) {
    fprintf(stderr,"No filename was specified\n");
    errorFlag = 1;
  }

  if (errorFlag) {
    fprintf(stderr,"Usage: getseries <filename> [-v] [-b <binning>] [-f <factor>] \n");
    fprintf(stderr," Returns: min mean max arcsec per pixel\n");
    fprintf(stderr,"          -v is the verbosity flag\n");
    fprintf(stderr,"          -b is the mosaic binning (default: 1)\n");
    fprintf(stderr,"          -f is the ratio of mean/min and max/mean (default: 1.1)\n");
    return(-1);
  }


  if (ParseFilename2(filename,series,&plateNumber) == 0) {
    fprintf(stderr,"ERROR: Failed to parse filename %s as %s%05d\n",
            filename,series,plateNumber);
    return(-1);
  }

  dasch_init_scandb(&my_connection);

  sprintf(queryString,"SELECT nominalPlateScale,fittedPlateScale FROM series where series = '%s';",series);
  if (strlen(queryString) > MAX_QUERY_STRING) {
    fprintf(stderr,"ERROR: MAX_QUERY_STRING exceeded %d\n",strlen(queryString));
    return(-1);
  }

  res = mysql_query(&my_connection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(&my_connection);

    if (res_ptr) {
      MYSQL_FIELD *field_ptr;

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        totalCount++;
        if (sqlrow[0] == NULL) {
          fprintf(stderr,"No nominal plate scale found for %s\n",filename);
          continue;
        }
        nvals = sscanf(sqlrow[0],"%lf",&nominalPlateScale);
        if (nvals != 1) {
          fprintf(stderr,"ERROR: nvals is %d for nominalPlateScale\n");
          continue;
        }
        if (sqlrow[1] != NULL) {
          nvals = sscanf(sqlrow[1],"%lf",&fittedPlateScale);
          if (nvals != 1) {
            fprintf(stderr,"ERROR: nvals is %d for fittedPlateScale\n");
          }
        }
        foundEntry++;

      }

      mysql_free_result(res_ptr);
    }
  } else {
    fprintf(stderr,"Select error %d: %s res %d\n",mysql_errno(&my_connection),mysql_error(&my_connection),res);
  }

  mysql_close(&my_connection);

  if (foundEntry != 1) {
    fprintf(stderr,"ERROR: %d database entries found\n",foundEntry);
    return(-1);
  }
  if (fittedPlateScale == 0.0) {
    fittedPlateScale = nominalPlateScale;
  }
  if (fittedPlateScale == 0.0) {
    fprintf(stderr,"ERROR: plate scale is %f for %s\n",fittedPlateScale,filename);
    return(-1);
  }
  fittedPlateScale = fittedPlateScale * binning * NOMINAL_MM_PER_PIXEL ;
  if (strcmp(series,"na") == 0) {
    fittedPlateScale = 126.491;
    factor2 = 12.6491;
    printf("%f %f %f\n",fittedPlateScale/factor2,fittedPlateScale,fittedPlateScale*factor2);
  } else {
    printf("%f %f %f\n",fittedPlateScale/factor,fittedPlateScale,fittedPlateScale*factor);
  }
  return(EXIT_SUCCESS);
}
