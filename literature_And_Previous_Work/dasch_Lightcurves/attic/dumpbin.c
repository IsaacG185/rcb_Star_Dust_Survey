// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* dumpbin.c
 *
 * This routine dumps the contents of a GSC bin.  It works only with the file version of the magnitudes file
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   dumpbin.c pipelineutils.a -ltable -lutil  -lwcs -pthread -o dumpbin -L/usr/lib/mysql -L/usr/lib64/mysql -lmysqlclient
 *
 * Nov 24, 2009 Edward J. Los - Initial version
 * Apr  2, 2010 Edward J. Los - Add a "q" qualifier for kepler catalog support 
 * Sep  7. 2010 Edward J. Los - Add SKIP_MYSQL for testing on the ice cluster
 * Jan  3, 2011 Edward J. Los - Add versionId
 * Mar 11, 2011 Edward J. Los - Add local_bin_index
 * Jun  3, 2011 Edward J. Los - Allow specification of the bin via ra and declination
 * Aug 31, 2011 Edward J. Los - Temporarily remove the REFNumber from the header sanity check (TEMP_REFNUMBER_FIX)
 * Nov  8, 2011 Edward J. Los - Add readonly mode
 *                              Add full dump of all fields
 *                              Support magnitude-dependent magnitude correction
 * Sep 21, 2012 Edward J. Los - Correct the sense of readonly (-u)
 *                              Add (-g) to create a script to check all files 
 * Mar 11, 2013 Edward J. Los - Add RaPM, DecPM, ra_2 and dec_2 to the magnitude file (Version 5);
 * Aug 19, 2013 Edward J. Los - Add support to dump catalog bin contents
 * Dec 23, 2014 Edward J. Los - V6 data format: add A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex
 * May 28, 2018 Edward J. Los - Add ucac5 support and gaia support
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 *
 *   dumpbin -b 111014150 -o bin1.tmp -h hdr1.tmp
 *   dumpbin -b 113039725 -o bin2.tmp
 *   dumpbin -r 335.31328527 -d -88.37324361 -n S310000131038 -o los.tmp -v
 *
 *   for test:    dumpbin -s -o los.dmp -h los.hdr -b 414720  -u
 *   for repair:  dumpbin -v -o los.dmp -h los.hdr -b 414720 -a  -f los.repair
 *
 *   to check all files: 
 *   dumpbin -o los.tmp -g -b 0 [-q apass]
 *   echo "source los.tmp >& dumpbin.log" | at now
 *   Took 45.83 hours on  Sep 25 04:29:27 EDT 2012
 *
 *   od -A d -t x4 -N 312 /dasch/filetest/mag111/mag013/mag111013888.dat
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

#define MAX_BUFFER 256

/* #define TEMP_REFNUMBER_FIX 1 */
/* #define OLD_GSC_CATALOG_FORMAT 1 */

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;


int main(int argc,char *argv[])
{
  int decBin;
  int raBin;
  char *argstr;
  int errorFlag = 0;
  int verbose = 0;
  int silent = 0;
  int repair = 0;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char headerfile[MAX_BUFFER];
  char repairfile[MAX_BUFFER];
  FILE *repairHandle = NULL;
  int nvals;
  FILE *outHandle = NULL;
  FILE *headerHandle = NULL;
  int curMagnitudes;
  PFILESTARIMAGE pFileStarImage = NULL;

  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon = &fileCommon;
  int magnitudeIndex;

  int tmp_gsc_bin_index = -1;
  int gsc_bin_index = -1;
  int gsc_bin_index_low;
  int gsc_bin_index_high;

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
  char InputREF[MAX_REF];
  int readOnly = 0;
  int listAllColumns = 0;
  int generateScript = 0;
	int starCatalogFlag = 0;

  char* catalogdir;
  char catalogname[MAX_BUFFER];
	char qualifier[MAX_BUFFER];
  char* slashPtr;
  char * dotPtr;
  char indexname[MAX_BUFFER];
  File catalogHandle;
  File indexHandle;
  int readItems;
  STARINDEX curStarIndex;
  PSTARINDEX pCurStarIndex = &curStarIndex;
	PGSCIMAGE pGscImageTable = NULL;
	PGSCIMAGE pGscImage;
	int gscIndex;
	char REF[MAX_REF];

  catalogString[0] = 0;
  outfile[0] = 0;
  headerfile[0] = 0;
  repairfile[0] = 0;
  InputREF[0] = 0;
  qualifier[0] = 0;
 
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d argv '%s'\n",argstr,argc,argstr);
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

        case 'h': /* header file name */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(headerfile,*++argv,MAX_BUFFER-2);
            if (strlen(headerfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;
        
        case 'N': /* object name */
        case 'n':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(InputREF,*++argv,MAX_REF-1);
            if (strlen(InputREF) >= (MAX_REF-1)) {
              printf("ERROR: InputREF name is too large for %s\n",*argv);
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
            nvals = sscanf(*++argv,"%d",&gsc_bin_index);
            if (nvals != 1) {
              printf("ERROR: Unable to decode gsc_bin_index %s\n",*argv);
              errorFlag = 1;
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
						strcpy(qualifier,*++argv);
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

        case 's': /* print only the file name, no output */
        case 'S':
          silent = 1;
          verbose = -1;
          break;

        case 'u': /*  */
        case 'U':
          readOnly = 1;
          break;

        case 'g': /*  */
        case 'G':
          generateScript = 1;
          break;
          
        case 'a': /*  */
        case 'A':
          listAllColumns = 1;
          break;

        case 'c': /*  */
        case 'C':
          starCatalogFlag = 1;
          break;

        case 'f': /* repair file name */
        case 'F':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(repairfile,*++argv,MAX_BUFFER-2);
            if (strlen(repairfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
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
  /* Validate arguments */

  if (outfile[0] == 0) {
    printf("ERROR: No output filename was specified\n");
    errorFlag = 1;
  }

  if ((declination > -98.0) &&
      (rightAscension > -98.0)) {
    tmp_gsc_bin_index = GetGSCBin(pGscBin,rightAscension,declination,&decBin,&raBin,"dumpbin");
  }
  if (gsc_bin_index < 0) {
    gsc_bin_index = tmp_gsc_bin_index;
  }


  if (gsc_bin_index < 0) {    
    printf("ERROR: No gsc bin index was specified\n");
    errorFlag = 1;
  }

  if ((silent != 0) && (generateScript != 0)) {
    printf("ERROR: -s and -g can not be used together\n");
    errorFlag = 1;
  }
  if (silent == 0) {

    outHandle = fopen(outfile,"wt");
    if (outHandle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to open the output file %s\n",outfile);
    } else {
      if (verbose > 0) {
        printf("Output file %s\n",outfile);
      }
    }
    if (generateScript) {
      printf("Generating integrity check commands for all files\n");
      fprintf(outHandle,"date\n");
      for (gsc_bin_index = 0; gsc_bin_index < pGscBin->total_gsc_bins; gsc_bin_index += MAG_FILE_MODULUS) {
        if (catalogString[0] == 0) {
          fprintf(outHandle,"dumpbin -s -o los.dmp -b %9d -u\n",gsc_bin_index); 
        } else {
          fprintf(outHandle,"dumpbin -s -o los.dmp -b %9d -q %s -u\n",gsc_bin_index,catalogString); 
        }
      }
      fprintf(outHandle,"date\n");

      fclose(outHandle);
      return(0);
    }


    if (headerfile[0] != 0) {
      headerHandle = fopen(headerfile,"wt");
      if (headerHandle == NULL) {
        errorFlag = 1;
        printf("ERROR: Failed to open the header file %s\n",headerfile);
      } else {
        if (verbose > 0) {
          printf("Header file %s\n",headerfile);
        }
      }
    }

    if (repairfile[0] != 0) {
      repairHandle = fopen(repairfile,"wt");
      if (repairHandle == NULL) {
        errorFlag = 1;
        printf("ERROR: Failed to open the repair file %s\n",repairfile);
      } else {
        if (verbose > 0) {
          printf("Repair file %s\n",repairfile);
        }
      }
    }
  } /* silent check */
  
  if (errorFlag) {
    printf("Usage: dumpbin options\n");
    printf("  options: -v verbose\n");
    printf("  options: -s silent print only the file name, no output\n");
    printf("           -o <output file>\n");
    printf("           -h <header file>\n");
    printf("           -b <gsc_bin_index>\n");
    printf("           -f <repair file name>\n");
    printf("           -q input catalog and filename qualifer\n");
    printf("           -r Right Ascension\n");
    printf("           -d Declination\n");
    printf("           -n Object name\n");
    printf("           -u Read Only: Do not attempt to write database files\n");
    printf("           -g Generate a script to check all files\n");
    printf("           -a List all columns\n");
    printf("           -c Read a star catalog file instead of the output data file\n");


    return(-1);
  }

	if (starCatalogFlag) {
		if (gsc_bin_index < 0) {
			printf("ERROR: gsc_bin_index %d is illegal for a star catalog dump\n");
			exit(-1);
		}
		if ((declination > -98.0) &&
				(rightAscension > -98.0)) {

			printf("Right Ascension %f, Declination %f, gsc_bin_index %d for REF %s\n",rightAscension,declination,gsc_bin_index,InputREF);
		}
		/* We are dumping from a star catalog */
	  catalogdir = getenv("DASCH_CATALOG");
		if (catalogdir == NULL) {
			printf("DASCH_CATALOG is not defined\n");
			printf("0\n");
			return(-1);
		}
		strcpy(catalogname,catalogdir);
		if (strstr(qualifier,"kepler")) {
			slashPtr = strrchr(catalogname,'/');
			if (slashPtr != NULL) {
				slashPtr++;
			} else {
				slashPtr = catalogname;
			}
			*slashPtr = 0;
			strcat(catalogname,"kepler.dat");
		}
		if (strstr(qualifier,"apass")) {
			slashPtr = strrchr(catalogname,'/');
			if (slashPtr != NULL) {
				slashPtr++;
			} else {
				slashPtr = catalogname;
			}
			*slashPtr = 0;
			strcat(catalogname,"apass.dat");
		}
		if (strstr(qualifier,"atlas")) {
			slashPtr = strrchr(catalogname,'/');
			if (slashPtr != NULL) {
				slashPtr++;
			} else {
				slashPtr = catalogname;
			}
			*slashPtr = 0;
			strcat(catalogname,"atlas.dat");
		}
		if (strstr(qualifier,"gaia")) {
			slashPtr = strrchr(catalogname,'/');
			if (slashPtr != NULL) {
				slashPtr++;
			} else {
				slashPtr = catalogname;
			}
			*slashPtr = 0;
			strcat(catalogname,"gaiadr2.dat");
		}
		if (strstr(qualifier,"experimental")) {
			slashPtr = strrchr(catalogname,'/');
			if (slashPtr != NULL) {
				slashPtr++;
			} else {
				slashPtr = catalogname;
			}
			*slashPtr = 0;
			strcat(catalogname,"experimental.dat");
		}
		if (strstr(qualifier,"ucac4")) {
			slashPtr = strrchr(catalogname,'/');
			if (slashPtr != NULL) {
				slashPtr++;
			} else {
				slashPtr = catalogname;
			}
			*slashPtr = 0;
			strcat(catalogname,"ucac4.dat");
		}
		if (strstr(qualifier,"ucac5")) {
			slashPtr = strrchr(catalogname,'/');
			if (slashPtr != NULL) {
				slashPtr++;
			} else {
				slashPtr = catalogname;
			}
			*slashPtr = 0;
			strcat(catalogname,"ucac5.dat");
		}
		strcpy(indexname,catalogname);
		dotPtr = strstr(indexname,".");
		if (dotPtr != NULL) {
			*dotPtr = 0;
		}
		strcat(indexname,".idx");
  
		catalogHandle = Open(catalogname,"r");
		if (catalogHandle == NULL) {
			printf("ERROR Could not open catalog file %s\n",catalogname);
			errorFlag = 1;
		} 
		indexHandle = Open(indexname,"r");
		if (indexHandle == NULL) {
			printf("ERROR Could not open catalog index file %s\n",indexname);
			errorFlag = 1;
		} 
		Seek(indexHandle,gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
		readItems = Read(indexHandle,pCurStarIndex,sizeof(STARINDEX),1);
		if (readItems != 1) {
			printf("ERROR: failed to read index of %s at gsc_bin_index %d\n",indexname,gsc_bin_index);
			exit(-1);
		}
		if (pCurStarIndex->numStars == 0) {
			printf("WARNING: No stars found in gsc_bin_index %d for %s\n",gsc_bin_index,catalogname);
		} else {
			pGscImageTable = (PGSCIMAGE)calloc(pCurStarIndex->numStars,sizeof(GSCIMAGE));
			if (verbose) {
				printf("Reading %d from %s\n",gsc_bin_index,catalogname);
			}
			if (pGscImageTable == NULL) {
				printf("ERROR: failed to allocate gsc image table of size %d\n",pCurStarIndex->numStars);
				exit(-1);
			}
			Seek(catalogHandle,pCurStarIndex->offset,SEEK_SET);
			readItems = Read(catalogHandle,pGscImageTable,sizeof(GSCIMAGE),pCurStarIndex->numStars);
			
			if (readItems != pCurStarIndex->numStars) {
				printf("ERROR: read only %d stars from a bin sized %d stars\n",readItems,pCurStarIndex->numStars);
				exit(-1);
			}
#ifdef OLD_GSC_CATALOG_FORMAT
      fprintf(outHandle,"REF\tra\tdec\tStdmag\tcolor\tRaPM\tDecPM\tclass\tVFlag\tMAGFlag\tflag\n");
      fprintf(outHandle,"---\t--\t---\t------\t-----\t----\t-----\t-----\t-----\t-------\t----\n");
#else /* OLD_GSC_CATALOG_FORMAT */
      fprintf(outHandle,"REF\tra\tdec\tStdmag\tcolor\tRaPM\tDecPM\tRaSigmaPM\tDecSigmaPM\tclass\tVFlag\tMAGFlag\tflag\n");
      fprintf(outHandle,"---\t--\t---\t------\t-----\t----\t-----\t---------\t----------\t-----\t-----\t-------\t----\n");
#endif /* OLD_GSC_CATALOG_FORMAT */

			for (gscIndex = 0; gscIndex < pCurStarIndex->numStars; gscIndex++) {
				pGscImage = &pGscImageTable[gscIndex];
#ifdef OLD_GSC_CATALOG_FORMAT
				fprintf(outHandle,"%s\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%d\n",
								pGscImage->REF,
								pGscImage->ra,
								pGscImage->dec,
								pGscImage->Stdmag,
								pGscImage->color,
								pGscImage->RaPM,
								pGscImage->DecPM,
								pGscImage->class,
								pGscImage->VFlag,
								pGscImage->MAGFlag,
								pGscImage->flag);

#else /* OLD_GSC_CATALOG_FORMAT */
				GetREF(pGscImage->REFNumber,REF,1,1);
				fprintf(outHandle,"%s\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%d\n",
								REF,
								pGscImage->ra,
								pGscImage->dec,
								pGscImage->Stdmag,
								pGscImage->color,
								pGscImage->RaPM,
								pGscImage->DecPM,
								pGscImage->RaSigmaPM,
								pGscImage->DecSigmaPM,
								pGscImage->class,
								pGscImage->VFlag,
								pGscImage->MAGFlag,
								pGscImage->flag);
#endif /* OLD_GSC_CATALOG_FORMAT */
			}
		}
		Close(indexHandle);
		Close(catalogHandle);
	} else {
		
		if (qualifier[0] != 0) {
			catalogNumber = GetCatalogNumber(qualifier);
    
			if (catalogNumber < 0) {
				printf("ERROR: Illegal catalog name %s\n",*argv);
				exit(-1);
			} else {
				sprintf(catalogString,"%d",catalogNumber);
			}
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
		InitFileCommon(stdout,pFileCommon,pConnection,pPhotConnection,catalogNumber,silent == 0);
#ifndef SKIP_MYSQL
		InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);
		mysql_close(pConnection);
#endif /* SKIP_MYSQL */
 
		InitCatalogAccess(pFileCommon);

		if (silent == 0) {
			printf("FILESTARHEADER size %d FILESTARIMAGE size %d readOnly %d allcolumns %d\n",sizeof(FILESTARHEADER),sizeof(FILESTARIMAGE),readOnly,listAllColumns);
		}
		if ((declination > -98.0) &&
				(rightAscension > -98.0)) {

			printf("Right Ascension %f, Declination %f, gsc_bin_index %d for REF %s\n",rightAscension,declination,gsc_bin_index,InputREF);
#if 1
			printf("UPDATE stars set gsc_bin_index = %d, updateFlag = 'yes' where REFNumber = %s;\n",gsc_bin_index,InputREF);
#endif

		}


		gsc_bin_index_low = gsc_bin_index;
		gsc_bin_index_high = gsc_bin_index;

		curMagnitudes = GetFileSummaryMagnitudes(pGscBin,pFileCommon,gsc_bin_index_low,gsc_bin_index_high,0,1,0,dumpAllFlag,headerHandle,catalogString,verbose,repairHandle,readOnly,0);

		if (silent == 0) {
			if (listAllColumns) {
				fprintf(outHandle,"REF\tseries\tX_IMAGE\tY_IMAGE\tMAG_ISO\tra\tdec\tmagcal_iso\tmagcal_iso_rms\tmagcal_local\tmagcal_local_rms\tDate\tFLUX_ISO\tMAG_APER\tMAG_AUTO\tKRON_RADIUS\tBACKGROUND\tFLUX_MAX\tTHETA_J2000\tELLIPTICITY\tISOAREA_WORLD\tFWHM_IMAGE\tFWHM_WORLD\tplate_dist\tBlendedmag\tlimiting_mag_local\tmagcal_local_error\tdradRMS2\tmagcor_local\textinction\tmagcal_magdep\tmagcal_magdep_rms\tgsc_bin_index\tplateNumber\tNUMBER\tversionId\tAFLAGS\tBFLAGS\tISO0\tISO1\tISO2\tISO3\tISO4\tISO5\tISO6\tISO7\tnpoints_local\trejectFlag\tmagdep_bin\tpassBits\tlocal_bin_index\tseriesId\texposureNumber\tsolutionNumber\tspatial_bin\tra_2\tdec_2\tRaPM\tDecPM\tA2FLAGS\tB2FLAGS\ttimeAccuracy\tmaskIndex\n");
				fprintf(outHandle,"---\t------\t-------\t-------\t-------\t--\t---\t----------\t--------------\t------------\t----------------\t----\t--------\t--------\t--------\t-----------\t----------\t--------\t-----------\t-----------\t-------------\t----------\t----------\t----------\t----------\t------------------\t------------------\t--------\t------------\t----------\t-------------\t-----------------\t-------------\t-----------\t------\t---------\t------\t------\t----\t----\t----\t----\t----\t----\t----\t----\t-------------\t----------\t----------\t--------\t---------------\t--------\t--------------\t--------------\t-----------\t----\t-----\t----\t-----\t-------\t-------\t-----------\t----------\n");
			} else {

				fprintf(outHandle,"REF\tseries\tplateNumber\tsolutionNumber\tgsc_bin_index\trejectFlag\tmagcal_local\tmagcal_local_error\tmagcal_iso_rms\tversionId\tlocal_bin_index\n");
				fprintf(outHandle,"---\t------\t-----------\t--------------\t-------------\t----------\t------------\t------------------\t--------------\t---------\t---------------\n");
			}
		}
		for (magnitudeIndex = 0; magnitudeIndex < curMagnitudes;magnitudeIndex++) {
			char REF[MAX_REF]; 
			pFileStarImage = &pFileCommon->magnitudeBuffer[magnitudeIndex];
#ifndef TEMP_REFNUMBER_FIX  
			if (GetREF(pFileStarImage->REFNumber,REF,1,0) == 0) {

#else /* TEMP_REFNUMBER_FIX */
				GetREF(pFileStarImage->REFNumber,REF,1,0);
#endif /* TEMP_REFNUMBER_FIX */

				if (silent == 0) {
					if (listAllColumns) {
						fprintf(outHandle,"%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%d\n",
										REF,
#if 0
										GetSeriesString(pFileStarImage->seriesId),
#else
										GetSeriesString(pFileStarImage->seriesId,0),
#endif
										pFileStarImage->X_IMAGE,
										pFileStarImage->Y_IMAGE,
										pFileStarImage->MAG_ISO,
										pFileStarImage->ra,
										pFileStarImage->dec,
										pFileStarImage->magcal_iso,
										pFileStarImage->magcal_iso_rms,
										pFileStarImage->magcal_local,
										pFileStarImage->magcal_local_rms,
										pFileStarImage->Date,
										pFileStarImage->FLUX_ISO,
										pFileStarImage->MAG_APER,
										pFileStarImage->MAG_AUTO,
										pFileStarImage->KRON_RADIUS,
										pFileStarImage->BACKGROUND,
										pFileStarImage->FLUX_MAX,
										pFileStarImage->THETA_J2000,
										pFileStarImage->ELLIPTICITY,
										pFileStarImage->ISOAREA_WORLD,
										pFileStarImage->FWHM_IMAGE,
										pFileStarImage->FWHM_WORLD,
										pFileStarImage->plate_dist,
										pFileStarImage->Blendedmag,
										pFileStarImage->limiting_mag_local,
										pFileStarImage->magcal_local_error,
										pFileStarImage->dradRMS2,
										pFileStarImage->magcor_local,
										pFileStarImage->extinction,
										pFileStarImage->magcal_magdep,
										pFileStarImage->magcal_magdep_rms,
										pFileStarImage->gsc_bin_index,
										pFileStarImage->plateNumber,
										pFileStarImage->NUMBER,
										pFileStarImage->versionId,
										pFileStarImage->AFLAGS,
										pFileStarImage->BFLAGS,
										pFileStarImage->ISO0,
										pFileStarImage->ISO1,
										pFileStarImage->ISO2,
										pFileStarImage->ISO3,
										pFileStarImage->ISO4,
										pFileStarImage->ISO5,
										pFileStarImage->ISO6,
										pFileStarImage->ISO7,
										pFileStarImage->npoints_local,
										pFileStarImage->rejectFlag,
										pFileStarImage->magdep_bin,
										pFileStarImage->passBits,
										pFileStarImage->local_bin_index,
										pFileStarImage->seriesId,
										pFileStarImage->exposureNumber,
										pFileStarImage->solutionNumber,
										pFileStarImage->spatial_bin,
										pFileStarImage->ra_2,
										pFileStarImage->dec_2,
										pFileStarImage->RaPM,
										pFileStarImage->DecPM,
                    pFileStarImage->A2FLAGS,
                    pFileStarImage->B2FLAGS,
                    pFileStarImage->timeAccuracy,
                    pFileStarImage->maskIndex);

					} else {
						fprintf(outHandle,"%s\t%s\t%05d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
										REF,
#if 0
										GetSeriesString(pFileStarImage->seriesId),
#else
										GetSeriesString(pFileStarImage->seriesId,0),
#endif
										pFileStarImage->plateNumber,
										pFileStarImage->solutionNumber,
										pFileStarImage->gsc_bin_index,
										pFileStarImage->rejectFlag,
										pFileStarImage->magcal_local,
										pFileStarImage->magcal_local_error,
										pFileStarImage->magcal_iso_rms,
										pFileStarImage->versionId,
										pFileStarImage->local_bin_index);
					}
				}
#ifndef TEMP_REFNUMBER_FIX 
			}
#endif /* TEMP_REFNUMBER_FIX */
		}
  
  
		mysql_close(pPhotConnection);
		FreeFileCommon(pFileCommon,silent == 0);

	} 

  if (silent == 0) {
    fclose(outHandle);
    if (headerfile[0] != 0) {
      fclose(headerHandle);
    }
    if (repairfile[0] != 0) {
      fclose(repairHandle);
    }
  }
  return(0);
}

