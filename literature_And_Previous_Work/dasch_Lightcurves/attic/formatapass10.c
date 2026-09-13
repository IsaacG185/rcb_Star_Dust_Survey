// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatapass10.c
 *
 * Read the apass input catalog and put it into the gsc232bin.db format
 *        Accept everything:               formatapass10 -v -n -o /dasch/Pipeline/catalogs/apass.db
 *
 *  1.  Use this program to reformat the apass catalog
 *        Accept nobs > 2, rms cutoff 0.1: formatapass10 -v -c 0.1 -o /dasch/Pipeline/catalogs/apass_temp.db
 *        Execution Time: 529 seconds; rejected because of color 1390 stars written 7173325
 *  
 *
 *  2.  Sort the results into 0.01 degree bins:
 *      
 *      index -mb -n apass_temp.db dec
 *      (takes 2m24s on dell)
 *
 *  4.  Run formatgsc to put this file in binary format and to index it
 *      
 *      formatgsc -r 1.0 -q apass /dasch/Pipeline/catalogs/apass_temp.db
 * 
 *  5.  build matchcatalogs.c with APASS_CATALOG, INCLUDE_KEPLER_COLOR, WRITE_MATCH_STARBASE and  WRITE_MERGE_CATALOG defined (no longer necessary)
 *      
 *      /dasch/Pipeline/matchcatalogs -c a -t 2.0 -k -m -w -f /dasch/Pipeline/catalogs/ucac3.dat -s /dasch/Pipeline/catalogs/apass_temp.dat
 *      cd /dasch/Pipeline/catalogs
 *      mv merge.idx apass.idx
 *      mv merge.dat apass.dat
 *      rm apass_temp.dat
 *      rm apass_temp.idx
 *      rm apass_temp_zero.db
 *
 *

cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  formatapass10.c pipelineutils.a -L /dasch/install/lib -lm   -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -o formatapass10


 *     matchcatalogs of Mar 23 2011 17:19:47, 
 *       matching /dasch/Pipeline/catalogs/ucac3.dat and /dasch/Pipeline/catalogs/apass_temp.dat 
 *       output /dasch/Pipeline/catalogs/matchgsckepler.db total bins 168966386 TOLERANCE 0.000556 degrees
 *      Done  populatedBinCount 6121155 matches 6292508 for 7173325 kepler stars at 699 seconds
 *            median degDrad is 0.000085 deg or 0.3 arcsec max is 0.000556 or 2.0 arcsec.
 *      mv merge.dat apass.dat
 *      mv merge.idx apass.idx
 *      rm apass_temp*
 *
 * Mar 18, 2019 Edward J. Los - Support DR10, adopted from formatapass.c
 * Apr  8, 2019 Edward J. Los - Reduce memory consumption, eliminate general statistics section
 * Apr  9, 2019 Edward J. Los - Add write status checking. Result:  "output_handle write failed -1 27 File too large"
 *                              convert to binary writes
 */ 
  


#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER 1024
#define MAGBINS  24   /* Number of magnitude bins */
#define APASS_VMAG   0
#define APASS_BMAG   1
#define APASS_GMAG   2
#define APASS_RMAG   3
#define APASS_IMAG   4
#define APASS_MAXMAG 5
/* #define PROPERMOTIONSTUDY 1 */


extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  
typedef struct _apassimage {
  int zone; /* I3.3   APASS south-polar-distance (SPD) zone */
  char hyphen; /* c      "-" */
  int starnum; /* I7.7   zone sequential star number */
  double radeg; /* F11.6  RA (decimal degrees) */
  double raerr; /* F7.3   RA error (arcsec) */
  double decdeg; /* F11.6  DEC (decimal degrees) */
  double decerr; /* F7.3   DEC error (arcsec) */
  unsigned short nobsB; /* I5     number of observed nights for B  */
  unsigned short nobsV; /* I5     number of observed nights for V */
  unsigned short nobsSU; /* I5     number of observed nights for u  (SU) */
  unsigned short nobsSG; /* I5     number of observed nights for g  (SG) */
  unsigned short nobsSR; /* I5     number of observed nights for r  (SR) */
  unsigned short nobsSI; /* I5     number of observed nights for i  (SI) */
  unsigned short nobsSZ; /* I5     number of observed nights for z_s  (SZ) */
  unsigned short nobsY; /* I5     number of observed nights for Y */
  double magB; /* F7.3   B */
  double magV; /* F7.3   V */
  double magSU; /* F7.3   u  (SU) */
  double magSG; /* F7.3   g  (SG) */
  double magSR; /* F7.3   r  (SR) */
  double magSI; /* F7.3   i  (SI) */
  double magSZ; /* F7.3   z_s (SZ) */
  double magY; /* F7.3   Y */
  double emagB; /* F7.3   err(B) */
  double emagV; /* F7.3   err(V) */
  double emagSU; /* F7.3   err(u)    (SU) */
  double emagSG; /* F7.3   err(g)    (SG) */
  double emagSR; /* F7.3   err(r)    (SR) */
  double emagSI; /* F7.3   err(i)    (SI) */
  double emagSZ; /* F7.3   err(i)    (SI) */
  double emagY; /* F7.3   err(Y) */
} APASSIMAGE,*PAPASSIMAGE;



int main(int argc,char *argv[])
{
  char *argstr;
  char output_name[MAX_INPUT_NAME];
  char raw_name[MAX_INPUT_NAME];
  int output_fd = 0;
  int raw_fd = 0;
  int outputCount = 0;
  int colorReject = 0;
  int skipCount = 0;
  int rawSkipCount = 0;
  int zeromagCount = 0;
  int rawSkipFlag;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int linecounter = 0;
  int index1;
  int index2;
  char *charPtr;
  char *charPtr2;
  char cmdchar;
  int verbose = 0;
  int acceptNegative = 0;
  char *dotPtr = NULL;
  int compareResult;
  int result;
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
  FILE *input_handle;
  char *apassMeasurementText[APASS_MAXMAG] = {"magV","magB","magSG","magSR","magSI"};
  double newrms;
  double scalerms = 0;
  double maxrms2 = 0;
  int intrms;
  int printf_result;
#if 0
#define MAX_INPUT_LINES   128632700
  /* dr10 128632651  input lines */
#define INPUT_FILES 1
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/apassdr10/zm90.txt",
  };

#endif
#if 0
#define MAX_INPUT_LINES   128632700
  /* dr10 128632651  input lines */
#define INPUT_FILES 3
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/apassdr10/zm35.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm30.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm25.txt",
  };

#endif
#ifdef PROPERMOTIONSTUDY

#define MAX_INPUT_LINES   128632700
  /* dr10 128632651  input lines */
#define INPUT_FILES 1
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/apassdr10/zm05.txt",    
	};
#endif /* PROPERMOTIONSTUDY */


#if 1                     
#define MAX_INPUT_LINES   128632700
  /* dr10 128632651  input lines */
#define INPUT_FILES 36
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/apassdr10/zm90.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm85.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm80.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm75.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm70.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm65.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm60.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm55.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm50.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm45.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm40.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm35.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm30.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm25.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm20.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm15.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm10.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zm05.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp00.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp05.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp10.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp15.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp20.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp25.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp30.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp35.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp40.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp45.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp50.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp55.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp60.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp65.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp70.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp75.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp80.txt",
    "/dasch/Pipeline/catalogs/apassdr10/zp85.txt"
	};
#endif


  int inputFile;


  int input_nrecs = 0;
  PAPASSIMAGE pInput = NULL;
  APASSIMAGE input_record;
  GSCAPASSIMAGE gsc_record;
  PGSCAPASSIMAGE pGscApassImage = &gsc_record;

  int tempmag;


  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int apassCount[APASS_MAXMAG];
  int apassDual[APASS_MAXMAG][APASS_MAXMAG];
  int magcount[APASS_MAXMAG][MAGBINS+1];
  int apassNoCQCount = 0;
  int apassHeaderCount = 0;
  int maxcqLength = 0;
  char maxcq[MAX_REF];
  double deltacolor;
  double maxdeltacolor = 0.0;
  int maxapasscolorid;
  char inLine[MAX_BUFFER];
  char output_buffer[MAX_BUFFER];
  int lineLen;
  char *inBuffer;
  int numLines;
  int nvals;
  double rmsCutoff = 9.0;
  double value1;
  double rms1;
  double value2;
  double rms2;
  long long REFNumber;
  int skipFlag; /* Skip this flag because of neg rms or too high rms */
  int skipFlagBase;
  int zeroMagBCount = 0;
  int zeroMagVCount = 0;
  int negZeroCount = 0;
  int posZeroCount = 0;
  ssize_t writeBytes;


  memset(apassCount,0,sizeof(apassCount));
  memset(apassDual,0,sizeof(apassDual));
  memset(magcount,0,sizeof(magcount));


  /* Loop through the arguments */
  output_name[0] = 0;
	raw_name[0] = 0;

  memset(magcount,0,sizeof(magcount));
  
#ifdef PROPERMOTIONSTUDY
  printf("ERROR: PROPERMOTIONSTUDY is defined\n");
#endif /* PROPERMOTIONSTUDY */


  maxcq[0] = 0;
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

        case 'n':
        case 'N':
          acceptNegative = 1;
          break;

        case 'c': /* RMS cutoff */
        case 'C':
          argc--;
          nvals = sscanf(*++argv,"%lf",&rmsCutoff);
          if (nvals != 1) {
            fprintf(stderr,"ERROR: Can not decode the rmsCutoff\n");
            errorFlag = 1;
          }
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

        case 'r': /* raw file name */
        case 'R':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(raw_name,*++argv,MAX_INPUT_NAME-2);
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

  /* Attempt to open the list */
  if (output_name[0] == 0) {
    printf("ERROR: Output file name not specified\n");
    errorFlag = 1;
  } else {
    output_fd = open(output_name,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRGRP);
    if (output_fd < 0) {
      errorFlag = 1;
      printf("Could not open output file %s\n",output_name);
    }
  }

  if (raw_name[0] != 0) {
		raw_fd = open(raw_name,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRGRP);
    if (raw_fd < 0) {
      errorFlag = 1;
      printf("Could not open raw file %s\n",raw_name);
    } else {
			sprintf(output_buffer,"ra\tdec\tStdmag\trms\n");
      if (strlen(output_buffer) > (MAX_BUFFER-2)) {
        printf("ERROR: line %d MAX_BUFFER exceeded\n",__LINE__);
        exit(-1);
      }
      writeBytes = write(raw_fd,output_buffer,strlen(output_buffer));
      if (writeBytes != strlen(output_buffer) ) {
        printf("ERROR line %d write error %d, %d %s\n",__LINE__,writeBytes,errno,strerror(errno));
        exit(-1);
      }


			sprintf(output_buffer,"--\t---\t------\t---\n");
      if (strlen(output_buffer) > (MAX_BUFFER-2)) {
        printf("ERROR: line %d MAX_BUFFER exceeded\n",__LINE__);
        exit(-1);
      }
      writeBytes = write(raw_fd,output_buffer,strlen(output_buffer));
      if (writeBytes != strlen(output_buffer) ) {
        printf("ERROR line %d write error %d, %d %s\n",__LINE__,writeBytes,errno,strerror(errno));
        exit(-1);
      }

		}
  }

  scalerms = rmsCutoff * sqrt(2.0);
  maxrms2 = sqr(scalerms);

  if (errorFlag) {
    printf("Usage: formatapass10 -o {output name} [-r <raw file>] [-v][-n][-c {cutoff}]\n");
		printf("       where -r <filename> is a raw table of RA and DEC for plotting \n");
    printf("       where -v is the verbose flag\n");
    printf("       where -c {cutoff} gives the RMS cutoff values\n");
    printf("       where -n accepts negative RMS values indicating only one measurement\n");

    return(-1);
  }

  printf("formatapass10 of %s %s Output Filename %s Raw filename %s acceptNegative %d,rmsCutoff %f\n",
         __DATE__,__TIME__,output_name,raw_name,acceptNegative,rmsCutoff);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of APASSIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(APASSIMAGE));
 

  time(&startTime);
  
  sprintf(output_buffer,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tgmag\trmag\timag\tvmagerr\tbmagerr\tgmagerr\trmagerr\timagerr\n");
  if (strlen(output_buffer) > (MAX_BUFFER-2)) {
    printf("ERROR: line %d MAX_BUFFER exceeded\n",__LINE__);
    exit(-1);
  }
  writeBytes = write(output_fd,output_buffer,strlen(output_buffer));
  if (writeBytes != strlen(output_buffer) ) {
    printf("ERROR line %d write error %d, %d %s\n",__LINE__,writeBytes,errno,strerror(errno));
    exit(-1);
  }


  sprintf(output_buffer,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t----\t----\t----\t-------\t-------\t-------\t-------\t-------\n");
  if (strlen(output_buffer) > (MAX_BUFFER-2)) {
    printf("ERROR: line %d MAX_BUFFER exceeded\n",__LINE__);
    exit(-1);
  }
  writeBytes = write(output_fd,output_buffer,strlen(output_buffer));
  if (writeBytes != strlen(output_buffer) ) {
    printf("ERROR line %d write error %d, %d %s\n",__LINE__,writeBytes,errno,strerror(errno));
    exit(-1);
  }




  input_nrecs = 0;
  for (inputFile = 0; inputFile < INPUT_FILES; inputFile++) {
    printf("Opening %s\n",inputnames[inputFile]);
    input_handle = fopen(inputnames[inputFile],"rt");
    if (input_handle == NULL) {
      printf("ERROR: failed to open %s\n",inputnames[inputFile]);
      exit(-1);
    }
    numLines = 0;
    while(1) {
      pInput = &input_record;
      memset(pInput,0,sizeof(APASSIMAGE));
      inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
      if (inBuffer == NULL) {
        break;
      }
      numLines++;
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
      charPtr2 = inBuffer;
      while ((charPtr = strstr(charPtr2,"-0.000")) != NULL) {
        charPtr += 5;
        charPtr2 = charPtr;
        negZeroCount++;
#if 0 /* Removed for APASS DR10 */
        *charPtr = '1';
#endif
      }
      charPtr2 = inBuffer;
      while ((charPtr = strstr(charPtr2,"0.000")) != NULL) {
        charPtr += 4;
        charPtr2 = charPtr;
        posZeroCount++;
#if 0  /* Removed for APASS DR10 */
        *charPtr = '1';
#endif
      }


      nvals = sscanf(inBuffer,"%3d%c%7d %lf %lf %lf %lf %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",

                     &pInput->zone,
                     &pInput->hyphen,
                     &pInput->starnum,
                     &pInput->radeg,
                     &pInput->raerr,
                     &pInput->decdeg,
                     &pInput->decerr,
                     &pInput->nobsB,
                     &pInput->nobsV,
                     &pInput->nobsSU,
                     &pInput->nobsSG,
                     &pInput->nobsSR,
                     &pInput->nobsSI,
                     &pInput->nobsSZ,
                     &pInput->nobsY,
                     &pInput->magB,
                     &pInput->magV,
                     &pInput->magSU,
                     &pInput->magSG,
                     &pInput->magSR,
                     &pInput->magSI,
                     &pInput->magSZ,
                     &pInput->magY,
                     &pInput->emagB,
                     &pInput->emagV,
                     &pInput->emagSU,
                     &pInput->emagSG,
                     &pInput->emagSR,
                     &pInput->emagSI,
                     &pInput->emagSZ,
                     &pInput->emagY);

      if (nvals != 31) {
        printf("ERROR: nvals is %d in line %d\n",nvals,numLines);
      } else {
#if 0
        printf("XXX nobs %d  mobs %d\n",pInput->nobs,pInput->mobs);
#endif
        rawSkipFlag = 0;
        input_nrecs++;
				if (raw_fd > 0) {
          if ((pInput->emagV > 0) && 
              (pInput->emagB > 0) &&
              (pInput->emagV < rmsCutoff) &&
              (pInput->emagB < rmsCutoff)) {
            newrms = sqr(pInput->emagV) + sqr(pInput->emagB);
            if (newrms < maxrms2) {
              newrms = sqrt(newrms);
            } else {
              newrms = 99.0;
              rawSkipCount++;
              rawSkipFlag = 1;
            }
          } else {
            newrms = 99.0;
            rawSkipCount++;
            rawSkipFlag = 1;
          }
          sprintf(output_buffer,"%f\t%f\t%.2f\t%2f\n",pInput->radeg,pInput->decdeg,pInput->magB,newrms);
          if (strlen(output_buffer) > (MAX_BUFFER-2)) {
            printf("ERROR: line %d MAX_BUFFER exceeded\n",__LINE__);
            exit(-1);
          }
          writeBytes = write(raw_fd,output_buffer,strlen(output_buffer));
          if (writeBytes != strlen(output_buffer) ) {
            printf("ERROR line %d raw_fd write error %d, %d %s\n",__LINE__,writeBytes,errno,strerror(errno));
            exit(-1);
          }


          
				}
				if (pInput->emagV == 0) {
					zeroMagVCount++;
				}
				if (pInput->emagB == 0) {
					zeroMagBCount++;
				}

				tempmag = pInput->magV-0.001;
				if (tempmag >= MAGBINS) {
					tempmag = MAGBINS;
				}

				for (index1 = 0; index1 < APASS_MAXMAG; index1++) {
					skipFlag = 0;
					skipFlagBase = 0;
					switch(index1) {
					case APASS_VMAG:
						value1 = pInput->magV;
						rms1   = pInput->emagV;
						break;
					case APASS_BMAG:
						value1 = pInput->magB;
						rms1   = pInput->emagB;
						break;
					case APASS_GMAG:
						value1 = pInput->magSG;
						rms1   = pInput->emagSG;
						break;
					case APASS_RMAG:
						value1 = pInput->magSR;
						rms1   = pInput->emagSR;
						break;
					case APASS_IMAG:
						value1 = pInput->magSI;
						rms1   = pInput->emagSI;
						break;
					default:
						printf("ERROR: line %d\n",__LINE__);
						exit(-1);
					}
					if (value1 > 90.00) {
						skipFlagBase = 1;
					}
					if ((rms1 < 0) && (acceptNegative == 0)) {
						skipFlagBase = 1;
						rms1 = - rms1;
					} else {
						if (rms1 < 0) {
							rms1 = - rms1;
						}
					}
					if (rms1 > rmsCutoff) {
						skipFlagBase = 1;
					}
					if (skipFlagBase == 0) {
						magcount[index1][tempmag]++;
						apassCount[index1]++;
					}
					for (index2 = 0; index2 < APASS_MAXMAG; index2++) {
						skipFlag = skipFlagBase;
						switch(index2) {
						case APASS_VMAG:
							value2 = pInput->magV;
							rms2   = pInput->emagV;
							break;
						case APASS_BMAG:
							value2 = pInput->magB;
							rms2   = pInput->emagB;
							break;
						case APASS_GMAG:
							value2 = pInput->magSG;
							rms2   = pInput->emagSG;
							break;
						case APASS_RMAG:
							value2 = pInput->magSR;
							rms2   = pInput->emagSR;
							break;
						case APASS_IMAG:
							value2 = pInput->magSI;
							rms2   = pInput->emagSI;
							break;
						default:
							printf("ERROR: line %d\n",__LINE__);
							exit(-2);
						}
						if (value2 > 90.00) {
							skipFlag = 1;
						}
						if ((rms2 < 0) && (acceptNegative == 0)) {
							skipFlag = 1;
							rms2 = -rms2;
						} else {
							if (rms2 < 0) {
								rms2 = - rms2;
							}
						}
						if (rms2 > rmsCutoff) {
							skipFlag = 1;
						}
						if ((index1 == APASS_BMAG) &&
                (index2 == APASS_VMAG)) {

              memset(pGscApassImage,0,sizeof(GSCAPASSIMAGE));
              result = GetDASCHNumber(pInput->radeg,pInput->decdeg,&REFNumber,1,REF_TYPE_APASS);
              pGscApassImage->REFNumber = REFNumber;
              pGscApassImage->ra = pInput->radeg;
              pGscApassImage->dec = pInput->decdeg;
              pGscApassImage->Stdmag =  pInput->magB;
              pGscApassImage->color = pInput->magB - pInput->magV;
              pGscApassImage->gmag = pInput->magSG;
              pGscApassImage->rmag = pInput->magSR;
              pGscApassImage->imag = pInput->magSI;
              pGscApassImage->vmagerr = pInput->emagB;
              pGscApassImage->bmagerr = pInput->emagV;
              pGscApassImage->gmagerr = pInput->emagSG;
              pGscApassImage->rmagerr= pInput->emagSR;
              pGscApassImage->imagerr = pInput->emagSI;
              if ((pGscApassImage->color < MIN_APASS_COLOR) ||
                  (pGscApassImage->color > MAX_APASS_COLOR)) {
                colorReject++;
                skipFlag = 1;
              }
              outputCount++;
              /* Covert the rms added in quadrature to an integer */
              newrms = (127.*(sqrt(sqr(rms1) + sqr(rms2))))/scalerms;
              if (newrms > 127) {
                newrms = 127;
#if 0
                if ( pGscApassImage->MAGFlag < 0) {
                  printf("MAGFLAG (2) is %d for outputCount %d\n", pGscApassImage->MAGFlag,outputCount);
                }

#endif
                
              }
              
              intrms = newrms;

              pGscApassImage->MAGFlag = newrms;
#if 0
              if ( pGscApassImage->MAGFlag < 0) {
                printf("MAGFLAG (1) is %d for outputCount %d\n", pGscApassImage->MAGFlag,outputCount);
              }

#endif
              if ((skipFlag) || (sqrt(sqr(rms1) + sqr(rms2)) > rmsCutoff))  {
                pGscApassImage->VFlag = 1;
                skipCount++;
              }
#if 0
              if ((pGscApassImage->ra < 185.5) ||
                  (pGscApassImage->ra > 186.0) ||
                  (pGscApassImage->dec < 63.5) ||
                  (pGscApassImage->dec > 64.5)) {
                continue;
              }
#endif
								
#if 0
              if ((rawSkipFlag == 1)  && (pGscApassImage->VFlag == 0)) {
                printf("X");
              }

#endif

#ifdef PROPERMOTIONSTUDY
              /* "Poster child" in notes13.txt */
              if ((fabs(pGscApassImage->ra - 323.81167) > 0.05) ||
                  (fabs(pGscApassImage->dec - -3.31201) > 0.05)) {
                continue;
              }
              if (pGscApassImage->Stdmag > 13.0) {
                continue;
              }
#endif /* PROPERMOTIONSTUDY */
              
              sprintf(output_buffer,"%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                      pGscApassImage->REFNumber,
                      pGscApassImage->ra,
                      pGscApassImage->dec,
                      pGscApassImage->Stdmag,
                      pGscApassImage->color,
                      pGscApassImage->class,
                      pGscApassImage->VFlag,
                      pGscApassImage->MAGFlag,
                      pGscApassImage->RaPM,
                      pGscApassImage->DecPM,
                      pGscApassImage->gmag,
                      pGscApassImage->rmag,
                      pGscApassImage->imag,
                      pGscApassImage->vmagerr,
                      pGscApassImage->bmagerr,
                      pGscApassImage->gmagerr,
                      pGscApassImage->rmagerr,
                      pGscApassImage->imagerr);
              if (strlen(output_buffer) > (MAX_BUFFER-2)) {
                printf("ERROR: line %d MAX_BUFFER exceeded\n",__LINE__);
                exit(-1);
              }
              writeBytes = write(output_fd,output_buffer,strlen(output_buffer));
              if (writeBytes != strlen(output_buffer) ) {
                printf("ERROR line %d output_fd write error %d, %d %s\n",__LINE__,writeBytes,errno,strerror(errno));
                exit(-1);
              }
           


            }
						if (skipFlag == 0) {
							apassDual[index1][index2]++;
						}

					}

				}


        

				if (input_nrecs > MAX_INPUT_LINES-10) {
					printf("ERROR: table too small at %d records\n",input_nrecs);
					break;
				}
			}
			if (verbose && ((input_nrecs % 100000) == 0)) {
				printf("Reading record %d\n",input_nrecs);
        time(&curTime);
        curTime -= startTime;

        printf("zeroMagBCount %d zeroMagVCount %d negZeroCount %d posZeroCount %d\n",zeroMagBCount,zeroMagVCount,negZeroCount,posZeroCount);
        printf("Execution Time: %d seconds; rejected %d rejected because of color %d stars written %d\n",curTime,skipCount,colorReject,outputCount);
        printf("rawSkipCount %d\n",rawSkipCount);

			}
        
		}
		printf("Input lines %d, input records %d\n",numLines,input_nrecs);
		fclose(input_handle);
    if (input_nrecs > MAX_INPUT_LINES-10) {
      break;
    }
	}
	if (raw_fd > 0) {
		close(raw_fd);
	}
	close(output_fd);
 
	time(&curTime);
	curTime -= startTime;

	printf("Measurement:");
	for (index1 = 0; index1 < APASS_MAXMAG; index1++) {
		printf("%10s",apassMeasurementText[index1]);
	}
	printf("\n");

	printf("Count:      ");
	for (index1 = 0; index1 < APASS_MAXMAG; index1++) {
		printf("%10d",apassCount[index1]);
	}
	printf("\n");

 
	for (index1 = 0; index1 < APASS_MAXMAG; index1++) {
		printf("%10s  ",apassMeasurementText[index1]);
		for (index2 = 0; index2 < APASS_MAXMAG; index2++) {
			printf("%10d",apassDual[index1][index2]);
		}
		printf("\n");
	}





	for (tempmag = 0; tempmag <= MAGBINS; tempmag++) {
		printf("MagB %2d     ",tempmag);
		for (index1 = 0; index1 < APASS_MAXMAG; index1++) {
			printf("%10d",magcount[index1][tempmag]);
		}
		printf("\n");
	}    



	time(&curTime);
	curTime -= startTime;

	printf("zeroMagBCount %d zeroMagVCount %d negZeroCount %d posZeroCount %d\n",zeroMagBCount,zeroMagVCount,negZeroCount,posZeroCount);
	printf("Execution Time: %d seconds; rejected %d rejected because of color %d stars written %d\n",curTime,skipCount,colorReject,outputCount);
  printf("rawSkipCount %d\n",rawSkipCount);

	return(EXIT_SUCCESS);
}
