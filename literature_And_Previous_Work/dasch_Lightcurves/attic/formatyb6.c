// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatyb6.c
 *
 * Read the yb6 input catalog and put it into the gsc232bin.db format
 *  1.  Use this program to reformat the yb6 catalog
 *      NOTE: the environment variable YB6_PATH must be defined
 *       setenv YB6_PATH "/dasch/yb6"
 *      formatyb6 -v -o /dasch/Pipeline/catalogs/yb6.db
 *  2.  Sort the results into 0.01 degree bins:
 *      
 *      index -mb -n yb6.db dec
 *
 *  4.  Run formatgsc to put this file in binary format and to index it
 *      
 *      formatgsc  /dasch/Pipeline/catalogs/yb6.db
 * 
 *      For statistics only:
 *      echo "formatgsc -s /dasch/Pipeline/catalogs/yb6.db" | at now
 *
 *    class is the proper motion quality   0 = Tycho -> 10
 *              reject proper motions below pmqual = 2        
 * 
 *  /data/astrocat/usnoyb6/
 *  scat -dh -c yb6 -n 400000 -r 1  19:51:23.12 17:05:57.894  J2000 
    USNO_YB6_number    RA2000     Dec2000   MagB   MagR   MagJ   MagH   MagK      Peak  Arcsec
    1070.0000002 297.8463333  17.0994167   16.53  17.12  15.20  14.58  14.38  121000   0.01
    REF            ra         dec       Stdmag     color    class VFlag MAGFlag RaPM      DecPM    
   Y1070.0000002   297.846333 17.099417 16.530000 -0.590000 10    0     0       -8.000000 -24.000000
    print rapm0 = -2.3249963409498298e-06 arcsec/year    -2.2222222222221645e-06 deg/yr -0.007999999999999792 arcsec/yr
    decpm0      = -6.666666666666609e-06                  -6.666666666666609e-06 deg/yr -0.023999999999999792 arcsec/yr

   star = {rasec = 107224680, decsec = 38555790, pm = 49884996, pmerr = 100000000, poserr = 0, mag = {16530, 17120, 15198, 14576, 14381}, magerr = {0, 0, 0, 0, 0}, index = {0, 0, 0, 0, 0}}

    grep '19:51:23.1' yb6s.txt | grep '17:05:57'
               19:51:23.12 +17:05:57.9 16.53 17.12 -0.0074 -0.0230  15.198  14.576  14.381
    055-04.all:19:51:23.12 +17:05:57.9 16.53 17.12 -0.0074 -0.0230  15.198  14.576  14.381
    yb6.cat:   19:51:23.12 +17:05:57.9 -0.0074 -0.0230 16.53 17.12  15.198  14.576  14.381
    yb6.txt:   19:51:23.12 +17:05:57.9 16.53 17.12 -0.0074 -0.0230  15.198  14.576  14.381

 *  scat -dh -c yb6 -n 400000 -r 1  19:21:48.513 +37:43:30.15  J2000 
          18:09:46.33 +28:22:04.3 17.71 17.37  0.0071 -0.0180  16.417  16.010  15.701
          19:29:38.31 +20:00:55.9 18.35 17.65  0.0053 -0.0168  14.782  14.319  14.140
 *  
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm formatyb6.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib/mysql  -l mysqlclient -o formatyb6
 * 
 *
 * Dec 11, 2010 Edward J. Los - Initial version  adapted from Doug Mink's ubcread.c in WCSTools
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
    int rasec, decsec, pm, pmerr, poserr, mag[5], magerr[5], index[5];
} YB6STAR1, *PYB6STAR1;
typedef struct {
  int rasec, decsec, pm, pmerr, poserr, mag[5], magerr[5], index[5], recnum;
} YB6STAR2, *PYB6STAR2;

/* UBCSWAP -- Reverse bytes of UB Catalog entry */

void ubcswap (char *string)	/* Start of vector of 4-byte ints */
{
char *sbyte, *slast;
char temp0, temp1, temp2, temp3;
 int nbytes = sizeof(YB6STAR1); /* Number of bytes to reverse */

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
	temp3 = sbyte[0];
	temp2 = sbyte[1];
	temp1 = sbyte[2];
	temp0 = sbyte[3];
	sbyte[0] = temp0;
	sbyte[1] = temp1;
	sbyte[2] = temp2;
	sbyte[3] = temp3;
	sbyte = sbyte + 4;
	}
    return;
}

/* UBCRA -- returns right ascension in degrees from the UB star structure */

double ubcra (int rasec)	/* RA in 100ths of arcseconds from UB catalog entry */
{
    return ((double) (rasec) / 360000.0);
}


/* UBCDEC -- returns the declination in degrees from the UB star structure */

double ubcdec (int decsec)	/* Declination in 100ths of arcseconds from UB catalog entry */
{
    return ((double) (decsec - 32400000) / 360000.0);
}
/* UBCMAG -- returns a magnitude from the UBC star structure */

double ubcmag (int magetc)	/* Magnitude 4 bytes from UB catalog entry */
{
    double xmag;

    xmag =  0.001 * (double) magetc;
    if (xmag == 0.00)
	xmag = 99.999;
    return (xmag);
}

/* UBCPRA -- returns RA proper motion in mas/year from UBC star structure */

double ubcpra (int magetc)	/* Proper motion field from UB catalog entry */
{
    double pm;

    if (magetc < 0)
	pm = (double) (-magetc % 10000);
    else
	pm = (double) (magetc % 10000);
    pm = ((pm * 0.002) - 10.0);
    pm = pm *1000.0;
    return (pm);
}


/* UBCPDEC -- returns Dec proper motion in mas/year from UBC star structure */

double ubcpdec (int magetc)	/* Proper motion field from UB catalog entry */
{
    double pm;
    pm = (double) (magetc / 10000);
    pm = ((pm * 0.002) - 10.0);
    pm = pm *1000.0;
    return (pm);
}


/* UBCPMQ -- returns proper motion probability (1-9) */

int ubcpmq (int magetc)	/* Quality, plate, and magnitude from UB catalog entry */
{
    if (magetc < 0)
	return (-magetc / 100000000);
    else
	return (magetc / 100000000);
}




int YB6Compare(const void *first, const void *second) 
{
  int numberFirst = ((PYB6STAR2)first)->decsec;
  int numberSecond = ((PYB6STAR2)second)->decsec;
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
  
  PYB6STAR1 input_table1 = NULL;
  PYB6STAR1 pInput1;
  PYB6STAR2 input_table2 = NULL;
  PYB6STAR2 pInput2;
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

  ucacpath = getenv("YB6_PATH");
  if (ucacpath == NULL) {
    printf("ERROR: YB6_PATH is not defined\n");
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
    printf("Usage: formatyb6 [-v][-s]\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s provides statistics only and skips writing output files\n");

    return(-1);
  }

  printf("formatyb6 of %s %s Output Filename %s\n",
         __DATE__,__TIME__,output_name);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of YB6STAR is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(YB6STAR1),sizeof(GSCIMAGE));
 

  time(&startTime);
  input_table1 = (PYB6STAR1)calloc(MAXRECORDS,sizeof(YB6STAR1));
  if (input_table1 == NULL) {
    fprintf(stderr,"ERROR: failed to allocate %d bytes for %d records\n",MAXRECORDS*sizeof(YB6STAR1),MAXRECORDS);
  }
  input_table2 = (PYB6STAR2)calloc(MAXRECORDS,sizeof(YB6STAR2));
  if (input_table2 == NULL) {
    fprintf(stderr,"ERROR: failed to allocate %d bytes for %d records\n",MAXRECORDS*sizeof(YB6STAR2),MAXRECORDS);
  }

  
  
  fprintf(output_handle,"REF\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\n");
  fprintf(output_handle,"---\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\n");

  /* Now read in the input file */


  for (zone = 0; zone <= 1799; zone++) {
#if 0
    zone = 1070;
#endif
    sprintf (zonepath, "%s/%03d/b%04d.cat", ucacpath,zone/10, zone);
    statResult = stat(zonepath,&filestats);
    if (statResult != 0) {
      /* No such file exists */
      /* fprintf (stderr,"UCACOPEN: failed to obtain file stats for UCAC file %s\n",zonepath); */
      continue;
    }
    filesize = filestats.st_size;
    filerecords = (int)(filesize/sizeof(YB6STAR1));
    if (filerecords > MAXRECORDS) {
      fprintf (stderr,"UCACOPEN: UCAC file %s has %d records while MAXRECORDS is %d\n",zonepath,filerecords,MAXRECORDS);
      exit(-1);
    }
    if (!(fcat = fopen (zonepath, "r"))) {
      fprintf (stderr,"UCACOPEN: UCAC file %s cannot be read\n",zonepath);
      exit(-1);
    }
    nbr = fread (input_table1, sizeof(YB6STAR1), filerecords, fcat);
    if (nbr != filerecords) {
      fprintf (stderr,"UCACOPEN: read %d of %d records for %s\n",nbr,filerecords,zonepath,filerecords);
      exit(-1);      
    }
    fclose(fcat);
    /* Swap the bytes */
    for (zone_nrecs = 0; zone_nrecs < filerecords; zone_nrecs++) {
      pInput1 = &input_table1[zone_nrecs];
      pInput2 = &input_table2[zone_nrecs];
      ubcswap ((char *)pInput1);
      memcpy(pInput2,pInput1,sizeof(YB6STAR1));
      pInput2->recnum = zone_nrecs+1;
    }
    /* Now sort the records in increasing declination */
    qsort((void*)input_table2,filerecords,sizeof(YB6STAR2),YB6Compare);
 
    for (zone_nrecs = 0; zone_nrecs < filerecords; zone_nrecs++) {
      pInput2 = &input_table2[zone_nrecs];
#if 0
      if (pInput2->recnum == 2) {
				printf("At record %d\n",pInput2->recnum);
      }
#endif

      memset(pGscImage,0,sizeof(GSCIMAGE));
      input_nrecs++;
      if (verbose) {
				if (((input_nrecs+1) % 1000000) == 0) {
					time(&curTime);
					curTime -= startTime;
					printf("At record %7d, time %5d sec\n",input_nrecs,curTime);
				}
      }
      if (pInput2->mag[0] == 0) {
				tempmag = MAGBINS;
      } else {
				tempmag = ubcmag(pInput2->mag[0]);
      }
      if (tempmag < 0) {
				tempmag = 0;
      }
      if (tempmag >= MAGBINS) {
				tempmag = MAGBINS;
      }
      magcount[tempmag]++;


      /* Now fill in our GSC record */
      sprintf(REF,"Y%04d.%07d",zone,pInput2->recnum); /* A preceeding 'Y' flags yb6 ids */
			if (GetREFNumber(REF,&pGscImage->REFNumber,&refType,0,0) != 0) {
				printf("ERROR: formatyb6: invalid REF %s\n",REF);
				exit(-1);
			}

      pGscImage->ra =  ubcra(pInput2->rasec);
      pGscImage->dec = ubcdec(pInput2->decsec);
      if (pInput2->mag[0] == 0) {
				pGscImage->Stdmag = 99.0;
      } else { 
				pGscImage->Stdmag = ubcmag(pInput2->mag[0]);
      }
      if ((pInput2->mag[0] == 0) || (pInput2->mag[1] == 0)) {
				pGscImage->color = 99.0;
      } else {
				pGscImage->color = ubcmag(pInput2->mag[0]) - ubcmag(pInput2->mag[1]);
      }
      pGscImage->RaPM = ubcpra(pInput2->pm);
      pGscImage->DecPM = ubcpdec(pInput2->pm);
      pGscImage->class = ubcpmq (pInput2->pm);
      if ((pGscImage->class) == 0) {
				pGscImage->class = 10;
      }
      
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
  if (input_table1 != NULL) {
    free(input_table1);
  }
  if (input_table2 != NULL) {
    free(input_table2);
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
