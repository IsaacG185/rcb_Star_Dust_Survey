// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* dumplimitingbin.c
 *
 * This routine dumps the contents of limiting magnitude bins
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   dumplimitingbin.c pipelineutils.a -ltable -lutil  -lwcs -pthread -o dumplimitingbin -L/usr/lib/mysql -L/usr/lib64/mysql -lmysqlclient
 *
 * May 11, 2013 Edward J. Los - Initial version, adapted from dumpbin.c
 * Apr  4, 2018 Edward J. Los - Search the index file for the most populated bins
 * Jun  1, 2018 Edward J. Los - Add plotting option (-g)
 *
 * 
    ./dumplimitingbin -v -n losnearby1.db -l loslimiting1.db -p v1 -q apass -r 130.092492 -d 19.672065
    ./dumplimitingbin -v -n losnearby2.db -l loslimiting2.db  -q apass -r 130.092492 -d 19.672065



dumplimitingbin -q gaia -g /home/scanner/scanner/linux/gmt/gaia_limiting.txt
plothammer -c -l 0 -u 163 -t 'GAIA calibrated plates' -i /home/scanner/scanner/linux/gmt/gaia_limiting.txt -o /home/scanner/scanner/linux/gmt/gaia_coverage.png
eog /home/scanner/scanner/linux/gmt/gaia_coverage.png

dumplimitingbin -q atlas -g /home/scanner/scanner/linux/gmt/atlas_limiting.txt
plothammer -c -l 0 -u 163 -t 'ATLAS calibrated plates' -i /home/scanner/scanner/linux/gmt/atlas_limiting.txt -o /home/scanner/scanner/linux/gmt/atlas_coverage.png
eog /home/scanner/scanner/linux/gmt/atlas_coverage.png


 */


#include <math.h>
#include <errno.h>
#include "table.h"
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include "pipelineutils.h"
#include "photometryutils.h"
#include "galaxyutils.h"

#define MAX_BUFFER 256

/* #define TEMP_REFNUMBER_FIX 1 */

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;


int main(int argc,char *argv[])
{
  int decBin;
  int raBin;
  char *argstr;
  int errorFlag = 0;
  int verbose = 0;
  int repair = 0;
  char cmdchar;
	char prefix[MAX_BUFFER];
  char nearbyfile[MAX_BUFFER];
	char limitingmagfile[MAX_BUFFER];
  char coveragefile[MAX_BUFFER];
	char qualifier[MAX_BUFFER];
  int nvals;
  FILE *nearbyHandle = NULL;
  FILE *limitingmagHandle = NULL;
  int curMagnitudes;
  PFILESTARIMAGE pFileStarImage = NULL;

  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon = &fileCommon;
  int magnitudeIndex;

  int tmp_gsc_bin_index = -1;
	int limiting_gsc_bin_index = -1;
	int nearby_gsc_bin_index = -1;

  char *mysqlhost;
  char *username;
  char *password;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  int catalogNumber = 0;
  int dumpAllFlag = 1;
  char catalogString[MAX_BUFFER];
  double rightAscension = -99.0;
  double declination = -99.0;
  double limitingRightAscension = -99.0;
  double limitingDeclination = -99.0;
	GALAXYCOMMON galaxycommon;
  PGALAXYCOMMON pGalaxyCommon = &galaxycommon;
	PGALAXYREC pGalaxyRec;
	int objectIndex;
	PPLATELIMITINGREC pPlateRec;
  int searchFlag = 0;
  int printLimit = 100; /* print 100 entries */

  catalogString[0] = 0;
	nearbyfile[0] = 0;
	limitingmagfile[0] = 0;
  coveragefile[0] = 0;
	qualifier[0] = 0;
	prefix[0] = 0;
  memset(pGalaxyCommon,0,sizeof(GALAXYCOMMON));

 
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


        case 'p': /* file prefix */
        case 'P':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(prefix,*++argv,MAX_BUFFER-2);
            if (strlen(prefix) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'n': /* nearby object file name */
        case 'N':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(nearbyfile,*++argv,MAX_BUFFER-2);
            if (strlen(nearbyfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'g': /* limiting magnitude coverge plot */
        case 'G':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(coveragefile,*++argv,MAX_BUFFER-2);
            if (strlen(coveragefile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'l': /* limiting magnitude file */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(limitingmagfile,*++argv,MAX_BUFFER-2);
            if (strlen(limitingmagfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        
        

        case 'b': /* BinNumber */
        case 'B':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&limiting_gsc_bin_index);
            if (nvals != 1) {
              printf("ERROR: Unable to decode limiting_gsc_bin_index %s\n",*argv);
              errorFlag = 1;
            } else {
							nearby_gsc_bin_index = limiting_gsc_bin_index;
						}
					}
					break;
          

        case 's': /* search with limit */
        case 'S':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&printLimit);
            if (nvals != 1) {
              printf("ERROR: Unable to decode printLimit %s\n",*argv);
              errorFlag = 1;
            } else {
              searchFlag = 1;
						}
					}
					break;
          

				case 'q': /* Catalog and file name qualifier */
				case 'Q':
					argc--;
					if (argc < 1) {
						printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						catalogNumber = GetCatalogNumber(*++argv);
            
						if (catalogNumber < 0) {
							printf("ERROR: Illegal catalog name %s\n",*argv);
							errorFlag = 1;
						} else {
							sprintf(catalogString,"%d",catalogNumber);
							strcpy(qualifier,*argv);
						}
					}
					break;

				case 'r': /* Right Ascension */
				case 'R':
					argc--;
					if (argc < 1) {
						printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						nvals = sscanf(*++argv,"%lf",&rightAscension);
						if (nvals != 1) {
							printf("ERROR: Unable to decode right ascension  %s\n",*argv);
							errorFlag = 1;
						}
					}
					break;
          
				case 'd': /* Declination */
				case 'D':
					argc--;
					if (argc < 1) {
						printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						nvals = sscanf(*++argv,"%lf",&declination);
						if (nvals != 1) {
							printf("ERROR: Unable to decode declination %s\n",*argv);
							errorFlag = 1;
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
	/* Validate arguments */

	if (OpenGalaxyFiles(pGalaxyCommon,qualifier,prefix) != 0) {
		printf("ERROR: failed to open the galaxy files for %s\n",qualifier);
		exit(-1);
	}
	if (verbose) {
		printf("Limiting magnitude file has %d gsc_bins with a bin size of %f degrees\n",pGalaxyCommon->pLimitingGscBin->total_gsc_bins,pGalaxyCommon->pLimitingGscBin->bin_size);
		if (pGalaxyCommon->nearbyDataHandle != NULL) {
			printf("Nearby Object file has %d gsc_bins with a bin size of %f degrees\n",pGalaxyCommon->pNearbyGscBin->total_gsc_bins,pGalaxyCommon->pNearbyGscBin->bin_size);
		}
	}

  if (coveragefile[0] != 0) {
    SearchLimiting(pGalaxyCommon,printLimit,coveragefile);
    exit(0);
  }

  if (searchFlag) {
    SearchLimiting(pGalaxyCommon,printLimit,NULL);
    exit(0);
  }


	if ((pGalaxyCommon->nearbyDataHandle != NULL) && 
			(pGalaxyCommon->pNearbyGscBin != pGalaxyCommon->pLimitingGscBin)) {
		nearby_gsc_bin_index = -1;
	}


	if ((declination > -98.0) &&
			(rightAscension > -98.0)) {
		tmp_gsc_bin_index = GetGSCBin(pGalaxyCommon->pLimitingGscBin,rightAscension,declination,&decBin,&raBin,"dumplimitingbin");
		if (limiting_gsc_bin_index < 0) {
			limiting_gsc_bin_index = tmp_gsc_bin_index;
		}
		if (pGalaxyCommon->nearbyDataHandle != NULL) {
			tmp_gsc_bin_index = GetGSCBin(pGalaxyCommon->pNearbyGscBin,rightAscension,declination,&decBin,&raBin,"dumplimitingbin");
			if (nearby_gsc_bin_index < 0) {
				nearby_gsc_bin_index = tmp_gsc_bin_index;
			}
		} else {
			nearby_gsc_bin_index = limiting_gsc_bin_index;
		}
	}
	if (verbose) {
		printf("Using limiting gsc_bin_index %d and nearby gsc_bin_index %d\n",limiting_gsc_bin_index,nearby_gsc_bin_index);
	}


	if (limiting_gsc_bin_index < 0) {    
		printf("ERROR: No gsc bin index was specified\n");
		errorFlag = 1;
	}

	if ((limitingmagfile[0] == 0) &&
			(nearbyfile[0] == 0)) {
		printf("ERROR either -f or -l must be specified\n");
		errorFlag == 1;
	}
	if (nearbyfile[0] != 0) {
		nearbyHandle = fopen(nearbyfile,"wt");
		if (nearbyHandle == NULL) {
			errorFlag = 1;
			printf("ERROR: Failed to open the nearby file %s\n",nearbyfile);
		} else {
			if (verbose > 0) {
				printf("Nearby file %s\n",nearbyfile);
			}
		}
	}

	if (limitingmagfile[0] != 0) {
		limitingmagHandle = fopen(limitingmagfile,"wt");
		if (limitingmagHandle == NULL) {
			errorFlag = 1;
			printf("ERROR: Failed to open the limitingmagput file %s\n",limitingmagfile);
		} else {
			if (verbose > 0) {
				printf("Limitingmagput file %s\n",limitingmagfile);
			}
		}
	}

  
	if (errorFlag) {
		printf("Usage: dumplimitingbin options\n");
		printf("  options: -v verbose\n");
		printf("  options: -n <nearby object file>\n");
		printf("           -l <limiting magnitude table>\n");
		printf("           -p <file prefix>\n");
		printf("           -b <gsc_bin_index>\n");
		printf("           -q input catalog and filename qualifer\n");
		printf("           -r Right Ascension\n");
		printf("           -d Declination\n");
    printf("           -s (search index for most populated bins\n");
    printf("           -g <coverage plot> provide a plot of limiting magnitude coverage\n");


		return(-1);
	}
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

#ifndef SKIP_MYSQL
	if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
		if (mysql_errno(pConnection)) {
			printf("ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
		}
		return(-1);
	}
#endif /* SKIP_MYSQL */

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

#ifndef SKIP_MYSQL
	if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
		if (mysql_errno(pPhotConnection)) {
			printf("ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
		}
		return(-1);
	}
#endif /* SKIP_MYSQL */

#ifndef SKIP_MYSQL
	InitSeriesTable(pConnection,pPhotConnection);
#endif /* SKIP_MYSQL */
#if 0
	InitFileCommon(stdout,pFileCommon,pConnection,pPhotConnection,catalogNumber,1);
#endif
#ifndef SKIP_MYSQL
#if 0
	InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);
#endif
	mysql_close(pConnection);
#endif /* SKIP_MYSQL */
#if 0 
	InitCatalogAccess(pFileCommon);
#endif
	if ((declination > -98.0) &&
			(rightAscension > -98.0)) {

		printf("Right Ascension %f, Declination %f, limiting gsc_bin_index %d nearby gsc_bin_index %d\n",rightAscension,declination,limiting_gsc_bin_index,nearby_gsc_bin_index);

	}
	if ((declination > -98.0) &&
			(rightAscension > -98.0)) {
		if (GetBinCenter(pGalaxyCommon->pLimitingGscBin,limiting_gsc_bin_index,&limitingDeclination,&limitingRightAscension,"dumplimitingbin") != 0) {
			printf("ERROR: GetBinCenter failed\n");
			exit(-1);
		}
	} else {
		limitingDeclination = declination;
		limitingRightAscension = rightAscension;
	}

	LoadGalaxyTable(pGalaxyCommon,limitingDeclination,limitingRightAscension);
	if (nearbyfile[0] != 0) {

		if (pGalaxyCommon->galaxyCount == 0) {
			printf("No nearby objects found in the limiting magnitude file\n");
		} else {
			if (verbose) {
				printf("Nearby objects in the table: %d\n",pGalaxyCommon->galaxyCount);
			}
			fprintf(nearbyHandle,"galaxyversion\tvariableFlag\tra\tdec\tcatalogmag\tradius\tcatalogname\tgalaxytype\tgalaxyflag\n");
			fprintf(nearbyHandle,"-------------\t------------\t--\t---\t----------\t------\t-----------\t----------\t----------\n");

			for (objectIndex = 0; objectIndex < pGalaxyCommon->galaxyCount; objectIndex++) {
				pGalaxyRec = &pGalaxyCommon->galaxyInputBuffer[objectIndex];
				fprintf(nearbyHandle,"%d\t%d\t%f\t%f\t%f\t%f\t%s\t%s\t%d\n",
								pGalaxyRec->galaxyversion,
								pGalaxyRec->variableFlag,
								pGalaxyRec->ra,
								pGalaxyRec->dec,
								pGalaxyRec->catalogmag,
								pGalaxyRec->radius,
								pGalaxyRec->catalogname,
								pGalaxyRec->galaxytype,							 
								pGalaxyRec->galaxyflag);
			}
		}

	}
	if (limitingmagfile[0] != 0) {

		if (pGalaxyCommon->plateCount == 0) {
			printf("No limiting magnitudes found\n");
		} else {
			fprintf(limitingmagHandle,"galaxyversion\tseriesId\tseries\tplateNumber\tmosaicNumber\tsolutionNumber\tversionId\tlimiting_mag_local\tgeoJulianDate\tgalaxyflag\n");
			fprintf(limitingmagHandle,"-------------\t--------\t------\t-----------\t------------\t--------------\t---------\t------------------\t-------------\t----------\n");
			if (verbose) {
				printf("Limiting magnitudes in the table: %d\n",pGalaxyCommon->plateCount);
			}

			for (objectIndex = 0; objectIndex < pGalaxyCommon->plateCount; objectIndex++) {
				pPlateRec = &pGalaxyCommon->plateLimitingBuffer[objectIndex];
				fprintf(limitingmagHandle,"%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\n",
								pPlateRec->galaxyversion,
								pPlateRec->seriesId,
								GetSeriesString(pPlateRec->seriesId,1),
								pPlateRec->plateNumber,
								pPlateRec->mosaicNumber,
								pPlateRec->solutionNumber,
								pPlateRec->versionId,
								pPlateRec->limiting_mag_local,
								pPlateRec->geoJulianDate,
								pPlateRec->galaxyflag);
			}
		}

	}
		



  
	mysql_close(pPhotConnection);

#if 0
	FreeFileCommon(pFileCommon,1);
#endif
	if (nearbyHandle != NULL) {
		fclose(nearbyHandle);
	}
	if (limitingmagHandle != NULL) {
		fclose(limitingmagHandle);
	}
	return(0);
}

