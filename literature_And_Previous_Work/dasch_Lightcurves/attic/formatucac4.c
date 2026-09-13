// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatucac4.c
 *
 * See "The fourth U.S. Naval Observatory CCD Astrograph Catalog (UCAC4)",
 *     Zacharias N., Finch C.T., Girard T.M., Henden A., Bartlett J.L., Monet D.G., Zacharias M.I.
 *     <Astron. J. (to be published)>=2012yCat.1322....0Z
 *
 * Read the ucac4 input catalog and put it into the gsc232bin.db format
 *
 * NOTE: maga (UCAC aperture magnitude) is given for Stdmag, no color is given.
 *       the obj_type is copied into the gsc class field.  Only zero obj_types have good astrometry.
 *     
 * NOTE: the highest proper motions are RaPM -4.41 to +6.77 arcsec/year and DecPM -5.82 to +10.33 arcsec/yr.
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm formatucac4.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -o formatucac4
 * 
 *  formatucac4 -v -o /dasch/Pipeline/catalogs/ucac4.db
 *
 * Jan 26, 2013 Edward J. Los - Initial version
 * Aug 12, 3013 Edward J. Los - Add proper motion standard deviations, Reformat GSCIMAGE
 */   


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#define MAX_INPUT_NAME 512

#define MAGBINS  21   /* Number of magnitude bins */
#define MAXRECORDS 287000 /* Largest file is zone 306 with 286833 records */


extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  
/* #define HIGH_PROPER_MOTION_CUTOFF 200.0 */ /* Proper motion cutoff in mas/yr */
/* #define DUMP_PROPER_MOTION 1 */   /* Print proper motion values */
#define MAXSIGMAPM 650  /* Units of 0.1 mas/yr */

#define OBJ_TYPE_CLEAN        1
#define OBJ_TYPE_OVEREXPOSED  2
#define OBJ_TYPE_STREAK       3
#define OBJ_TYPE_HPM          4
#define OBJ_TYPE_HPM_EXTERNAL 5
#define OBJ_TYPE_POORMOTION   6
#define OBJ_TYPE_SUPPLEMENT   7
#define OBJ_TYPE_HPMNOPPMXL   8
#define OBJ_TYPE_HPMBADPPMXL  9
/* The following is from ucac4.h */
#pragma pack( 1)
typedef struct ucac4_star {
   int32_t ra, spd;         /* RA/dec at J2000.0,  ICRS,  in milliarcsec */
   uint16_t mag1, mag2;     /* UCAC fit model & aperture mags, .001 mag */
   uint8_t mag_sigma;
   uint8_t obj_type, double_star_flag;
   int8_t ra_sigma, dec_sigma;    /* sigmas in RA and dec at central epoch */
   uint8_t n_ucac_total;      /* Number of UCAC observations of this star */
   uint8_t n_ucac_used;      /* # UCAC observations _used_ for this star */
   uint8_t n_cats_used;      /* # catalogs (epochs) used for prop motion */
   uint16_t epoch_ra;        /* Central epoch for mean RA, minus 1900, .01y */
   uint16_t epoch_dec;       /* Central epoch for mean DE, minus 1900, .01y */
   int16_t pm_ra;            /* prop motion, .1 mas/yr = .01 arcsec/cy */
   int16_t pm_dec;           /* prop motion, .1 mas/yr = .01 arcsec/cy */
   int8_t pm_ra_sigma;       /* sigma in same units */
   int8_t pm_dec_sigma;
   uint32_t twomass_id;        /* 2MASS pts_key star identifier */
   uint16_t mag_j, mag_h, mag_k;  /* 2MASS J, H, K_s mags,  in millimags */
   uint8_t icq_flag[3];
   uint8_t e2mpho[3];          /* 2MASS error photometry (in centimags) */
   uint16_t apass_mag[5];      /* in millimags */
   uint8_t apass_mag_sigma[5]; /* also in millimags */
   uint8_t yale_gc_flags;      /* Yale SPM g-flag * 10 + c-flag */
   uint32_t catalog_flags;
   uint8_t leda_flag;          /* LEDA galaxy match flag */
   uint8_t twomass_ext_flag;   /* 2MASS extended source flag */
   uint32_t id_number;
   uint16_t ucac2_zone;
   uint32_t ucac2_number;
} UCAC4STAR,*PUCAC4STAR;
#pragma pack( )


#define UCAC4_ASCII_SIZE 280

         /* This flag suppresses stars that were matched with Tycho       */
         /* stars.  In some of my software,  stars are drawn from both    */
         /* UCAC-4 and Tycho.  By setting this flag in the output_format  */
         /* field,  I keep those stars from being plotted twice.          */
#define UCAC4_OMIT_TYCHO_STARS            0x1

         /* By default,  zero magnitudes and proper motions are written    */
         /* out as zeroes.  Setting this 'output_format' flag causes them  */
         /* to be written out as spaces.                                   */
#define UCAC4_WRITE_SPACES                0x2
#define UCAC4_FORTRAN_STYLE               0x4

         /* "Raw binary" means the data is written out in the same 78 byte */
         /* per star format in which it was read: no reformatting is done. */
#define UCAC4_RAW_BINARY                  0x8
/* Fields in the u4hpm.dat file */
#define MAX_U4HPM_RECORDS                 50
typedef struct u4hpm {
  int rnm;   /* unique star identifier (col. 51 main data) */
  int zn;    /* UCAC4 zone number for this star */
  int rnz;   /* unning record number along that zone for this star */
  int pmrc;  /* actual proper motion in RA*cos(dec) [0.1 mas/yr] */
  int pmd;   /* actual proper motion in Dec         [0.1 mas/yr] */
  int ra;    /* col. 1 of main data file (redundant) */
  int dec;   /* col. 2 of main data file (redundant) */
  int maga;  /* col. 4 of main data file (redundant) */
} U4HPM,*PU4HPM;

int UCAC4Compare(const void *first, const void *second) 
{
  int numberFirst = ((PUCAC4STAR)first)->spd;
  int numberSecond = ((PUCAC4STAR)second)->spd;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}


int main(int argc,char *argv[])
{
  char *argstr;
  char ucacpath[] = "/dasch/ucac4/u4b";
  char zonepath[MAX_INPUT_NAME];
  char hpmpath[] = "/dasch/ucac4/u4i/u4hpm.dat";
  FILE *fcat;
  FILE *hpmHandle;
  char inLine[MAX_INPUT_NAME];
  char *inBuffer;
  U4HPM u4hpm_table[MAX_U4HPM_RECORDS];
  PU4HPM pU4hpm = NULL;
  int hpmLines = 0;
  int hpmIndex;
  int hpmCount = 0;
	int hpmFlag;
  int lineLen;
  int nvals;
  char output_name[MAX_INPUT_NAME];
  File output_handle = NULL;
  int outputCount = 0;
  int zeromagCount = 0;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int linecounter = 0;
  char *charPtr;
  char cmdchar;
  int verbose = 0;
  char *dotPtr = NULL;
  int compareResult;
  int binNumber;
  int oldBinNumber = -1;
  int dec_bin;
  int dec_bin2;
  int oldDecBin = 0;
  int savedDecBin = -2;
  int ra_bin;
  int maxDecBin = 0;
  double maxRaError = 0;
  double maxDecError = 0;
  int maxWriteCount = 0;
  int totalBinCount = 0;

  File input_handle = NULL;
  int input_nrecs = 0;
  int zone_nrecs = 0;
  int input_index;
  size_t nbr;
  
  PUCAC4STAR input_table;
  PUCAC4STAR pInput;
  GSCIMAGE gsc_record;
  PGSCIMAGE pGscImage = &gsc_record;

 
  double tmpRaPM;
  double tmpDecPM;
  double tmpPM;
  double maxRaPM = 0;
  double maxDecPM = 0;
  double maxPM = 0;
  int countPM = 0;
  int tempmag;

  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int magcount[MAGBINS+1];
  int zone;
  int statResult;
  struct stat filestats;
  off_t filesize;
  int filerecords;
  int min_id_number = -1;
  int max_id_number = -1;
	int properMotionCount = 0;
	int pm_sigma;
	double PMSigma;
	FILE *properMotionHandle = NULL;

	int sigmaPMHist[MAXSIGMAPM+1];
	double totSigma;
	int sigmaPMIndex;
	int sigmaPMCumulative = 0;
	char REF[MAX_REF];
	int refType;

	memset(sigmaPMHist,0,sizeof(sigmaPMHist));
#ifdef HIGH_PROPER_MOTION_CUTOFF
	printf("ERROR: HIGH_PROPER_MOTION_CUTOFF is set to %f\n",HIGH_PROPER_MOTION_CUTOFF);
#endif /* HIGH_PROPER_MOTION_CUTOFF */
#ifdef DUMP_PROPER_MOTION
	printf("ERROR: DUMP_PROPER_MOTION is set, writing to lospm.db\n");
	properMotionHandle = fopen("lospm.db","wt");
	if (properMotionHandle == NULL) {
		printf("ERROR writing lospm.db\n");
		exit(-1);
	}
	fprintf(properMotionHandle,"Stdmag\tclass\tRaSigmaPM\tDecSigmaPM\n");
	fprintf(properMotionHandle,"------\t-----\t---------\t----------\n");
#endif /* DUMP_PROPER_MOTION */

  /* Loop through the arguments */
  output_name[0] = 0;


  memset(magcount,0,sizeof(magcount));
  

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
          verbose = 1;
          break;

        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(output_name,*++argv,MAX_INPUT_NAME-2);
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



  if (output_name[0] == 0) {
    printf("ERROR: Output file name not specified\n");
    errorFlag = 1;
  } else {
    output_handle = fopen(output_name,"wt");
    if (output_handle == NULL) {
      errorFlag = 1;
      printf("Could not open output file %s\n",output_name);
    }
  }



  if (errorFlag) {
    printf("Usage: formatucac4 [-v][-s]\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s provides statistics only and skips writing output files\n");

    return(-1);
  }

  printf("formatucac4 of %s %s Output Filename %s\n",
         __DATE__,__TIME__,output_name);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of UCAC4STAR is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(UCAC4STAR),sizeof(GSCIMAGE));
 

  time(&startTime);
  input_table = (PUCAC4STAR)calloc(MAXRECORDS,sizeof(UCAC4STAR));
  if (input_table == NULL) {
    fprintf(stderr,"ERROR: failed to allocate %d bytes for %d records\n",MAXRECORDS*sizeof(UCAC4STAR),MAXRECORDS);
  }

  
  
  fprintf(output_handle,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tRaSigmaPM\tDecSigmaPM\n");
  fprintf(output_handle,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t---------\t----------\n");
  /* Read in the file of high proper motion stars */
  
  hpmHandle = fopen(hpmpath,"rt");
  if (hpmHandle == NULL) {
    printf("ERROR: failed to open %s\n",hpmpath);
    exit(-1);
  }
  hpmLines = 0;
  while(1) {
    if (hpmLines >= MAX_U4HPM_RECORDS) {
      printf("ERROR: MAX_U4HPM_RECORDS exceeded\n");
      exit(-1);
    }
    pU4hpm = &u4hpm_table[hpmLines];
    memset(pU4hpm,0,sizeof(U4HPM));
    inBuffer = fgets(inLine,MAX_INPUT_NAME,hpmHandle);
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
    nvals = sscanf(inBuffer,"%d %d %d %d %d %d %d %d",
                   &pU4hpm->rnm,
                   &pU4hpm->zn,
                   &pU4hpm->rnz,
                   &pU4hpm->pmrc,
                   &pU4hpm->pmd,
                   &pU4hpm->ra,
                   &pU4hpm->dec,
                   &pU4hpm->maga);
    if (nvals != 8) {
      printf("Error reading line %d of u4hpm %s\n",hpmLines,inLine);
    } else {            
      hpmLines++;
    }
  }
  fclose(hpmHandle);
  if (verbose) {
    printf("Read %d lines from %s\n",hpmLines,hpmpath);
  }
  /* Now read in the input file */


  for (zone = 1; zone <= 900; zone++) {
    sprintf (zonepath, "%s/z%03d", ucacpath, zone);
    statResult = stat(zonepath,&filestats);
    if (statResult != 0) {
      /* No such file exists */
      fprintf (stderr,"UCACOPEN: failed to obtain file stats for UCAC file %s\n",zonepath);
      exit(-1);
    }
    filesize = filestats.st_size;
    filerecords = (int)(filesize/sizeof(UCAC4STAR));
    if (filerecords > MAXRECORDS) {
      fprintf (stderr,"UCACOPEN: UCAC file %s has %d records while MAXRECORDS is %d\n",zonepath,filerecords,MAXRECORDS);
      exit(-1);
    }
    if (!(fcat = fopen (zonepath, "r"))) {
      fprintf (stderr,"UCACOPEN: UCAC file %s cannot be read\n",zonepath);
      exit(-1);
    }
    nbr = fread (input_table, sizeof(UCAC4STAR), filerecords, fcat);
    if (nbr != filerecords) {
      fprintf (stderr,"UCACOPEN: read %d of %d records for %s\n",nbr,filerecords,zonepath,filerecords);
      exit(-1);      
    }
    fclose(fcat);
    /* Now sort the records in increasing declination */
    qsort((void*)input_table,filerecords,sizeof(UCAC4STAR),UCAC4Compare);
 
    for (zone_nrecs = 0; zone_nrecs < filerecords; zone_nrecs++) {
			hpmFlag = 0;
      pInput = &input_table[zone_nrecs];
      memset(pGscImage,0,sizeof(GSCIMAGE));
      input_nrecs++;
      if (verbose) {
        if (((input_nrecs+1) % 1000000) == 0) {
          time(&curTime);
          curTime -= startTime;
          printf("At record %7d, time %5d sec\n",input_nrecs,curTime);
        }
      }
      if (pInput->mag2 == 20000) {
        tempmag = MAGBINS;
      } else {
        tempmag = (1.0*pInput->mag2)/1000.0;
      }
      if (tempmag < 0) {
        tempmag = 0;
      }
      if (tempmag >= MAGBINS) {
        tempmag = MAGBINS;
      }
      magcount[tempmag]++;
      if (min_id_number == -1) {
        min_id_number = pInput->id_number;
        max_id_number = pInput->id_number;
      } else {
        if (pInput->id_number < min_id_number) {
          min_id_number = pInput->id_number;
        }
        if (pInput->id_number > max_id_number) {
          max_id_number = pInput->id_number;
        }
      }
        

      /* Now fill in our GSC record */
      sprintf(REF,"U4%09d",pInput->id_number); /* A preceeding 'U4' flags ucac4 ids */
			if (GetREFNumber(REF,&pGscImage->REFNumber,&refType,0,0) != 0) {
				printf("ERROR: formatucac4: invalid REF %s\n",REF);
				exit(-1);
			}
			pGscImage->ra =  (1.0 * pInput->ra)  / 3600000.0;	/* mas */
      pGscImage->dec = ((1.0 * pInput->spd) / 3600000.0);
      pGscImage->dec = pGscImage->dec - 90.0;	/* referenced from the South Pole */
      pGscImage->class = pInput->obj_type; /* Temporary place for the object type.  only object types of zero are used to update astrometry */
      if ((pInput->apass_mag[0] != 20000) && (pInput->apass_mag[1] != 20000)) {
        /* Use APASS magnitudes */
        pGscImage->Stdmag = (1.0 * pInput->apass_mag[0]) / 1000.0;
        pGscImage->color = (1.0 * (pInput->apass_mag[0]-pInput->apass_mag[1])) / 1000.0;
      } else {
        /* Use the aperture magnitude */
        if (pInput->mag2 == 20000) {
          pGscImage->Stdmag = 99.0;
        } else { 
          pGscImage->Stdmag = (1.0 * pInput->mag2) / 1000.0;
        }
        pGscImage->color = 99.0;
      }
      if ((pInput->pm_ra == 32767) || (pInput->pm_dec == 32767)) {
        for (hpmIndex = 0; hpmIndex < hpmLines; hpmIndex++) {
          pU4hpm = &u4hpm_table[hpmIndex];
          if ((zone == pU4hpm->zn) &&
              (pInput->id_number == pU4hpm->rnm) &&
              (pInput->ra == pU4hpm->ra) &&
              (pInput->spd == pU4hpm->dec) &&
              (pInput->mag2 == pU4hpm->maga)) {
            hpmCount++;
						hpmFlag = 1;
            /* We have a match */
            if (verbose) {
              printf("Found HPM entry %2d\n",hpmIndex);
            }
            pGscImage->RaPM = (1.0 * pU4hpm->pmrc) / 10.0; /* 0.1 mas/yr */
            pGscImage->DecPM = (1.0 * pU4hpm->pmd) / 10.0; /* 0.1 mas/yr */
            break;
          }
        }
        if (hpmIndex == hpmLines) {
          printf("ERROR: no proper motion found for zone %d entry %d\n",zone,pInput->id_number);
        }
      } else {

        pGscImage->RaPM = (1.0 * pInput->pm_ra) / 10.0; /* 0.1 mas/yr */
        pGscImage->DecPM = (1.0 * pInput->pm_dec) / 10.0; /* 0.1 mas/yr */
      }
#if 0
			if ((pGscImage->ra < 187.27) ||
					(pGscImage->ra > 187.28) ||
					(pGscImage->dec < 2.05) ||
					(pGscImage->dec > 2.06)) {
				continue;
			}
#endif
			if ((pGscImage->class == OBJ_TYPE_HPMNOPPMXL) ||
					(pGscImage->class == OBJ_TYPE_HPMBADPPMXL)) {
				pGscImage->RaSigmaPM = 0;
				pGscImage->RaPM = 0;
				pGscImage->DecSigmaPM = 0;
				pGscImage->DecPM = 0;
			} else {

				pm_sigma =  pInput->pm_ra_sigma + 128;
				switch (pm_sigma) {
				case 251:
					PMSigma = 27.5;
					break;
				case 252:
					PMSigma = 32.5;
					break;
				case 253:
					PMSigma = 37.5;
					break;
				case 254:
					PMSigma = 45.5;
					break;
				case 255:
					PMSigma = 50.0;
					break;
				default:
					if (pm_sigma >= 251) {
						printf("ERROR: pm_sigma %d\n",pm_sigma);
						exit(-1);
					}
					PMSigma = (1.0 * pm_sigma)/ 10.0; /* 0.1 mas/yr */
					break;
				}
				if (PMSigma < 50.0) {
					pGscImage->RaSigmaPM = PMSigma;
				} else {
					pGscImage->RaSigmaPM = 0;
					pGscImage->RaPM = 0;
				}

				pm_sigma =  pInput->pm_dec_sigma + 128;
				switch (pm_sigma) {
				case 251:
					PMSigma = 27.5;
					break;
				case 252:
					PMSigma = 32.5;
					break;
				case 253:
					PMSigma = 37.5;
					break;
				case 254:
					PMSigma = 45.5;
					break;
				case 255:
					PMSigma = 50.0; /* No data case */
					break;
				default:
					if (pm_sigma >= 251) {
						printf("ERROR: pm_sigma %d\n",pm_sigma);
						exit(-1);
					}
					PMSigma = (1.0 * pm_sigma)/ 10.0; /* 0.1 mas/yr */
					break;
				}
				if (PMSigma < 50.0) {
					pGscImage->DecSigmaPM = PMSigma;
				} else {
					pGscImage->DecSigmaPM = 0;
					pGscImage->DecPM = 0;
				}
#ifdef DUMP_PROPER_MOTION
				fprintf(properMotionHandle,"%.2f\t%d\t%.2f\t%.2f\n",pGscImage->Stdmag,pGscImage->class,pGscImage->RaSigmaPM,pGscImage->DecSigmaPM);
				if (++properMotionCount > 1000000) {
					printf("ERROR: properMotionCount exceeds %d\n",properMotionCount);
					exit(-1);
				}

#endif /* DUMP_PROPER_MOTION */

			}
			totSigma = sqrt(sqr(pGscImage->DecSigmaPM)+sqr(pGscImage->RaSigmaPM));
			sigmaPMIndex = totSigma*10.0;
			if (sigmaPMIndex >= MAXSIGMAPM) {
				sigmaPMHist[MAXSIGMAPM]++;
			} else {
				sigmaPMHist[sigmaPMIndex]++;
			}

#ifdef HIGH_PROPER_MOTION_CUTOFF
#if 0
			if (hpmFlag == 0) {
				continue;
			}
#else
			if (((pGscImage->RaPM*pGscImage->RaPM) + (pGscImage->DecPM*pGscImage->DecPM)) < 
					((HIGH_PROPER_MOTION_CUTOFF*HIGH_PROPER_MOTION_CUTOFF))) {
				continue;
			}
#endif
			
#endif /* HIGH_PROPER_MOTION_CUTOFF */
#if 1
			if (strcmp(REF,"U4000200137") == 0) {
				printf("%s\n",REF);
			}
#endif


			outputCount++;
			fprintf(output_handle,"%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",
							pGscImage->REFNumber,
							pGscImage->ra,
							pGscImage->dec,
							pGscImage->Stdmag,
							pGscImage->color,
							pGscImage->class,
							pGscImage->VFlag,
							pGscImage->MAGFlag,
							pGscImage->RaPM,
							pGscImage->DecPM,
							pGscImage->RaSigmaPM,
							pGscImage->DecSigmaPM);
		}
		if (verbose) {
			printf("Zone %3d has %8d records\n",zone,zone_nrecs);
		}
	}
	if (input_table != NULL) {
		free(input_table);
	}

	Close(output_handle);
 
	time(&curTime);
	curTime -= startTime;

	printf("min_id_number %d, max_id_number %d hpmStars %d\n",min_id_number,max_id_number,hpmCount);
	printf("Execution Time: %d seconds for %d records; output %d \n",curTime,input_nrecs,outputCount);

	for (tempmag = 0; tempmag <= MAGBINS; tempmag++) {
		printf("Mag %3d     ",tempmag);
		printf("%8d",magcount[tempmag]);
		printf("\n");
	}
	printf("Proper Motion Error (mas/yr),count,cumulative count:\n");
	for (sigmaPMIndex = 0; sigmaPMIndex <= MAXSIGMAPM; sigmaPMIndex++) {
		if (sigmaPMHist[sigmaPMIndex] > 0) {
			totSigma = (1.0*sigmaPMIndex)/10.0;
			sigmaPMCumulative += sigmaPMHist[sigmaPMIndex];
			printf("%4.1f %9d %9d\n",totSigma,sigmaPMHist[sigmaPMIndex],sigmaPMCumulative);
		}

	}

	if (properMotionHandle != NULL) {
		fclose(properMotionHandle);
	}

	return(EXIT_SUCCESS);
}
