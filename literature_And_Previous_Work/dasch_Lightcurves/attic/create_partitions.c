// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* create_partitions.c
 *
 * This program creates the MAGNITUDES_MODULUS magnitudes tables labeled magnitudesmnnn where nnn is 0 to 
 * MAGNITUDE_MODULUS-1 and m, currently 0, is the pass indicator
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   create_partitions.c pipelineutils.a -ltable -lutil  -lwcs -o create_partitions -L/usr/lib64/mysql -lmysqlclient
 *
 * Jun 23, 2008 Edward J. Los - Original Version
 *
 *
 */


#include <math.h>
#include <errno.h>
#include "table.h"
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include "pipelineutils.h"
#include "photometryutils.h"


#define MAX_BUFFER 256
#define MAX_QUERY_STRING_LARGE 5000

int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  char cmdchar;
  char *mysqlphothost;
  char username[MAX_BUFFER];
  char password[MAX_BUFFER];
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  time_t startTime;
  time_t curTime;
  char queryString[MAX_QUERY_STRING_LARGE];
  int res;
  int querySize;
  int partitionIndex;
  int catalogNumber = 0;
  int nvals;
  username[0] = 0;
  password[0] = 0;

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
        case 'u': /* user name */
        case 'U':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for list file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(username,*++argv,MAX_BUFFER-2);
            if (strlen(username) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'p': /* password */
        case 'P':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(password,*++argv,MAX_BUFFER-2);
            if (strlen(password) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'c': /* Catalog number */
        case 'C':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&catalogNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode catalogNumber %s\n",*argv);
              errorFlag = 1;
            } else {
              if ((catalogNumber < 0) || (catalogNumber > 9)) {
                fprintf(stderr,"ERROR: Illegal catalog number %d\n",catalogNumber);
              }
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
    if (username[0] == 0) {
      fprintf(stderr,"ERROR: MySQL user name was specified\n");
      errorFlag = 1;
    }
  if (errorFlag) {
    printf("Usage: create_partitions options\n");
    printf("  options: -u <username>\n");
    printf("           -p <password>\n");
    printf("           -c <catalog number (0 through 9)>\n");

    return(-1);
  }


  SetQueryCount(0);
  printf("create_partitions of %s %s, MAGNITUDES_MODULUS %d\n",
         __DATE__,__TIME__,MAGNITUDES_MODULUS);
    

  /* Connect to the database */
  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    fprintf(stderr,"DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }



  mysql_init(pPhotConnection);
  
  if (!mysql_real_connect(pPhotConnection,mysqlphothost,username,password,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
    }
    return(-1);
  }
  time(&startTime);

  for (partitionIndex = 0; partitionIndex < MAGNITUDES_MODULUS; partitionIndex++) {
#if 0 /* Highly dangerous! Backup the database first */
    sprintf(queryString,"DROP TABLE IF EXISTS magnitudes%d%03d;",catalogNumber,partitionIndex);
    printf("%s\n",queryString);
    res = ExecuteQuery(pPhotConnection,queryString);
#endif
    /* The following string is copied directly from photometry.sql */
    sprintf(queryString,"CREATE TABLE magnitudes%d%03d (   gsc_bin_index INT,                               /* GSC bin index */   seriesId TINYINT UNSIGNED,                       /* Plate series */   plateNumber MEDIUMINT UNSIGNED,                                 /* Plate number */   SEXTRNUMBER INT,                                 /* Sextractor object number */   passBits SMALLINT UNSIGNED,                      /* Pass identifier, currently only bit 1: GSC2.3.2 pass1 is recognized */   versionId INT,                                   /* Pipeline version Id */   exposureNumber TINYINT UNSIGNED,                              /* Exposure number (zero based) */   spatial_bin TINYINT UNSIGNED,                                 /* Spatial bin */   local_bin_index SMALLINT UNSIGNED,                             /* Local bin index: ix+(nx*iy) where ix is the horizontal zero-based index */   REFNumber BIGINT,                                /* Star reference designation - first digit is catalog */                                                    /*    For GSC2.3.2, replace a 'N' with 11 and a 'S' with 12 */                                                    /*    NONE - 0 */   X_IMAGE DOUBLE,                                  /* X pixel location */   Y_IMAGE DOUBLE,                                  /* Y pixel location */   AFLAGS INT,                                      /* Sextractor and pipeline flags */   BFLAGS INT,                                      /* Sextractor and pipeline flags */   MAG_ISO DOUBLE,                                  /* Sextractor isophotal estimate */   ra DOUBLE,                                       /* Object Right Ascension in degrees */   declination DOUBLE,                              /* Object Declination in degrees */   magcal_iso DOUBLE,                               /* Lowess magnitude estimate */   magcal_iso_rms DOUBLE,                           /* Lowess magnitude estimate error */   magcal_local DOUBLE,                             /* Local magnitude estimate */   magcal_local_rms DOUBLE,                         /* Local magnitude estimate error */   Date DOUBLE,                                     /* Julian date of exposure */   FLUX_ISO DOUBLE,                                 /* Peak Sextractor flux */   MAG_APER DOUBLE,                                 /* Sextractor aperture magnitude */   MAG_AUTO DOUBLE,                                 /* Sextractor auto magnitude */   KRON_RADIUS DOUBLE,                              /* Sextractor Kron radius */   BACKGROUND DOUBLE,                               /* Sextractor background estimate */   FLUX_MAX DOUBLE,                                 /* Sextractor maximum flux in ADU */   THETA_J2000 DOUBLE,                              /* Sextractor major axis angle */   ELLIPTICITY DOUBLE,                              /* Sextractor ellipticity */   ISOAREA_WORLD DOUBLE,                            /* Sextractor area in square degrees */   FWHM_IMAGE DOUBLE,                               /* Width of the image in pixels */   FWHM_WORLD DOUBLE,                               /* Width of the image in degrees */   ISO0 INT,                                        /* Isophotal level 0 area in sq pixels */   ISO1 INT,                                        /* Isophotal level 1 area in sq pixels */   ISO2 INT,                                        /* Isophotal level 2 area in sq pixels */   ISO3 INT,                                        /* Isophotal level 3 area in sq pixels */   ISO4 INT,                                        /* Isophotal level 4 area in sq pixels */   ISO5 INT,                                        /* Isophotal level 5 area in sq pixels */   ISO6 INT,                                        /* Isophotal level 6 area in sq pixels */   ISO7 INT,                                        /* Isophotal level 7 area in sq pixels */   plate_dist DOUBLE,                               /* Plate distance from center in degrees */   Blendedmag DOUBLE,                               /* Blended magnitude */   limiting_mag_local DOUBLE,                       /* Locally corrected limiting magnitude */   PRIMARY KEY (gsc_bin_index,seriesId,plateNumber,SEXTRNUMBER,passBits)      ) COMMENT 'Table of all magnitude estimates';",catalogNumber,partitionIndex);

    querySize = strlen(queryString);
    printf("querySize is %d for partition %d%03d\n",querySize,catalogNumber,partitionIndex);
    if (querySize > (MAX_QUERY_STRING_LARGE)) {
      fprintf(stderr,"ERROR: MAX_QUERY_STRING_LARGE exceeded %d\n",querySize);
      exit(-1);
    }


    res = ExecuteQuery(pPhotConnection,queryString);

  }


  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds %d queries\n",curTime,GetQueryCount());



  mysql_close(pPhotConnection);

  return(EXIT_SUCCESS);
}

