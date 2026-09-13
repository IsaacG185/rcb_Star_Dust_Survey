// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatucac3.c
 *
 * Read the ucac3 input catalog and put it into the gsc232bin.db format
 *  1.  Use this program to reformat the ucac3 catalog
 *      NOTE: the environment variable UCAC3_PATH must be defined
 *       setenv UCAC3_PATH "/dasch/ucac3"
 *      formatucac3 -v -o /dasch/Pipeline/catalogs/ucac3.db
 *  2.  Sort the results into 0.01 degree bins:
 *      
 *      index -mb -n ucac3.db dec
 *
 *  4.  Run formatgsc to put this file in binary format and to index it
 *      
 *      formatgsc  /dasch/Pipeline/catalogs/ucac3.db
 * 
 *      For statistics only:
 *      echo "formatgsc -s /dasch/Pipeline/catalogs/ucac3.db" | at now
 *
 *  scat -dh -c ucac3 -n 400000 -r 1  19:21:48.513 +37:43:30.15  J2000 
 *  NOTE: Stdmag is mmag (UCAC fit model magnitude)  Color is bmag-rmag (SC Bmag - SC R2mag)      
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm formatucac3.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib/mysql  -l mysqlclient -o formatucac3
 * 
 *
 * Dec 11, 2010 Edward J. Los - Initial version  adapted from Doug Mink's ucacread.c in WCSTools
 */   


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/stat.h>
#define MAX_INPUT_NAME 512

#define MAGBINS  19   /* Number of magnitude bins */
#define MAXRECORDS 637350 /* Largest file is zone 123 with 637320 records */


extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  

typedef struct {
    int rasec, decsec;
    short mmag, amag, sigmag;
    char objt, dsf;
    short sigra, sigdec;
    char na1, nu1, us1, cn1;
    short cepra, cepdec;
    int rapm, decpm;
    short sigpmr, sigpmd;
    int id2m;
    short jmag, hmag, kmag;
    char jq, hq, kq;
    char e2mpho[3];
    short bmag, rmag, imag;
    char clbl;
    char bq, rq, iq;
    char catflg[10];
    char g1, c1, leda, x2m;
    int rn;
} UCAC3STAR, *PUCAC3STAR;



int UCAC3Compare(const void *first, const void *second) 
{
  int numberFirst = ((PUCAC3STAR)first)->decsec;
  int numberSecond = ((PUCAC3STAR)second)->decsec;
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
  char *ucacpath;
  char zonepath[MAX_INPUT_NAME];
  FILE *fcat;
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
  
  PUCAC3STAR input_table;
  PUCAC3STAR pInput;
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
	char REF[MAX_REF];
	int refType;

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

  ucacpath = getenv("UCAC3_PATH");
  if (ucacpath == NULL) {
    printf("ERROR: UCAC3_PATH is not defined\n");
    errorFlag = 1;
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
    printf("Usage: formatucac3 [-v][-s]\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s provides statistics only and skips writing output files\n");

    return(-1);
  }

  printf("formatucac3 of %s %s Output Filename %s\n",
         __DATE__,__TIME__,output_name);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of UCAC3STAR is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(UCAC3STAR),sizeof(GSCIMAGE));
 

  time(&startTime);
  input_table = (PUCAC3STAR)calloc(MAXRECORDS,sizeof(UCAC3STAR));
  if (input_table == NULL) {
    fprintf(stderr,"ERROR: failed to allocate %d bytes for %d records\n",MAXRECORDS*sizeof(UCAC3STAR),MAXRECORDS);
  }

  
  
  fprintf(output_handle,"REF\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\n");
  fprintf(output_handle,"---\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\n");

  /* Now read in the input file */


  for (zone = 1; zone <= 360; zone++) {
    sprintf (zonepath, "%s/z%03d", ucacpath, zone);
    statResult = stat(zonepath,&filestats);
    if (statResult != 0) {
      /* No such file exists */
      fprintf (stderr,"UCACOPEN: failed to obtain file stats for UCAC file %s\n",zonepath);
      exit(-1);
    }
    filesize = filestats.st_size;
    filerecords = (int)(filesize/sizeof(UCAC3STAR));
    if (filerecords > MAXRECORDS) {
      fprintf (stderr,"UCACOPEN: UCAC file %s has %d records while MAXRECORDS is %d\n",zonepath,filerecords,MAXRECORDS);
      exit(-1);
    }
    if (!(fcat = fopen (zonepath, "r"))) {
      fprintf (stderr,"UCACOPEN: UCAC file %s cannot be read\n",zonepath);
      exit(-1);
    }
    nbr = fread (input_table, sizeof(UCAC3STAR), filerecords, fcat);
    if (nbr != filerecords) {
      fprintf (stderr,"UCACOPEN: read %d of %d records for %s\n",nbr,filerecords,zonepath,filerecords);
      exit(-1);      
    }
    fclose(fcat);
    /* Now sort the records in increasing declination */
    qsort((void*)input_table,filerecords,sizeof(UCAC3STAR),UCAC3Compare);
 
    for (zone_nrecs = 0; zone_nrecs < filerecords; zone_nrecs++) {
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
      if (pInput->mmag == 0) {
				tempmag = MAGBINS;
      } else {
				tempmag = ((double)pInput->mmag)/1000.0;
      }
      if (tempmag < 0) {
				tempmag = 0;
      }
      if (tempmag >= MAGBINS) {
				tempmag = MAGBINS;
      }
      magcount[tempmag]++;


      /* Now fill in our GSC record */

      sprintf(REF,"U%09d",pInput->rn); /* A preceeding 'U' flags ucac3 ids */
			if (GetREFNumber(REF,&pGscImage->REFNumber,&refType,0,0) != 0) {
				printf("ERROR: formatucac3: invalid REF %s\n",REF);
				exit(-1);
			}

      pGscImage->ra =  (double) pInput->rasec  / 3600000.0;	/* mas */
      pGscImage->dec = (double) pInput->decsec / 3600000.0;	/* mas */
      pGscImage->dec = pGscImage->dec - 90.0;
      if (pInput->mmag == 0) {
				pGscImage->Stdmag = 99.0;
      } else { 
				pGscImage->Stdmag = ((double) pInput->mmag) / 1000.0;
      }
      if ((pInput->bmag == 0) || (pInput->rmag == 0)) {
				pGscImage->color = 99.0;
      } else {
				pGscImage->color = (((double) pInput->bmag) - ((double) pInput->rmag))/1000.0;
      }
      pGscImage->RaPM = ((double)pInput->rapm) / 10.0; /* 0.1 mas/yr */
      pGscImage->DecPM = ((double)pInput->decpm) / 10.0; /* 0.1 mas/yr */
      
      
      outputCount++;
      fprintf(output_handle,"%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\n",
							pGscImage->REFNumber,
							pGscImage->ra,
							pGscImage->dec,
							pGscImage->Stdmag,
							pGscImage->color,
							pGscImage->class,
							pGscImage->VFlag,
							pGscImage->MAGFlag,
							pGscImage->RaPM,
							pGscImage->DecPM);
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

  printf("Execution Time: %d seconds for %d records; output %d \n",curTime,input_nrecs,outputCount);

  for (tempmag = 0; tempmag <= MAGBINS; tempmag++) {
    printf("Mag %3d     ",tempmag);
    printf("%8d",magcount[tempmag]);
    printf("\n");
  }    




  return(EXIT_SUCCESS);
}
