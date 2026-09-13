// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formattycho2.c
 *
 * Read the Tycho2  catalog and put it into the gsc232bin.db format
 *        Accept everything:               formattycho2 -s -v -n -o /dasch/Pipeline/catalogs/tycho2.db
 *
 *  1.  Use this program to reformat the Tycho2 catalog
 *        Accept nobs > 2, rms cutoff 0.1: formattycho2 -v -c 0.1 -o /dasch/Pipeline/catalogs/tycho_temp.db
 *        Execution Time: 529 seconds; rejected because of color 1390 stars written 7173325
 *  2.  Combine this catalog with tycho_temp.db and continue with the instructions in formatapass.c 
 *      Be sure that the -c <cutoff> value is identical for both catalogs.
 *
 *      NOTE: The transformation used here comes from the following URL:
 *       http://www.cfa.harvard.edu/kepler/kic/format/format.html
 *      TYBV          The parent catalog is Tycho-2, and this gives B and
 *                    V magnitudes.  The KEPMAG value was computed from
 *                           sdssg = 0.54*b + 0.46*v - 0.07;
 *                           sdssr = -0.44*b + 1.44*v + 0.12;
 *                    and the same transformation as above to compute KEPMAG.
 *                    The errors in KEPMAG follow those for the Tycho-2
 *                    catalog, being a few hundredths for brighter stars
 *                    and perhaps as large as a tenth for fainter ones,
 *                    with a few hundredths added in quadrature for the
 *                    uncertainty of the transformation.
 * 
 * gcc -ggdb -O0  -I/n/sw/plplot-5.9.9/include/plplot -I/usr/include/plplot -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include plothammer.c pipelineutils.a -L /dasch/install/lib -L/n/sw/plplot-5.9.9/lib -lplplotd  -lm  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread    -ltable -lutil  -lwcs -lgd  tclstub.o -o plothammer
 * 
 *  formattycho2 -v -s -n -o /dasch/Pipeline/catalogs/tycho2_temp.db
 * 
 * Aug 14, 2011 Edward J. Los - Initial version adapted from formatapass.c
 * Sep 16, 2011 Edward J. Los - Use B and V magnitudes; correct rms cut.
 * Aug  8, 2013 Edward J. Los - Output REFNumber instead of REF.
 * Aug  1, 2015 Edward J. Los - Because of the need for the transient candidate flare search, store the catalog on disk in the same format as other catalogs        
 *                              We will create tycho2.dat sorted by gsc bin index and tycho2.idx.  The original binary will be Unlike
 *                              the other catalogs, each entry will have a copy in an adjacent bin for the purposes of sorting.
 * Mar 14, 2015 Edward J. Los - correct missing REFNumber in the binary catalog
*/ 
  

#include "table.h"
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "libwcs/wcscat.h"
#include "kdtree.h"

extern GSCBIN gscBin08;  
PGSCBIN pGscBin = &gscBin08; /* WARNING: any change here must be reflected in FindNearestTycho2Star() of galaxyutils.c */

int ty2regn(int	region,		/* Region to find */
            int	*star1,		/* First star number in region (returned)*/
            int	*star2,		/* Last star number in region (returned)*/
            int	verbose);	/* 1 for diagnostics */

struct StarCat *ty2open (int nstar, int nread);
void ty2close (struct StarCat *sc);
int ty2star (
             struct StarCat *sc,	/* Star catalog descriptor */
             struct Star *st,	/* Current star entry */
             int istar);	/* Star sequence number in Tycho 2 catalog region file */

#define MAX_INPUT_NAME 512
#define MAX_BUFFER 256
#define MAGBINS  19   /* Number of magnitude bins */
#define TYCHO2_BTMAG   0
#define TYCHO2_VTMAG   1
#define TYCHO2_MAXMAG  2
/* #define SKIP_TRANSFORM 1 */
#define USE_BV 1


/* #define LOS_DEBUG 1 */  

typedef struct _apassimage {
  long long starid;
  double ra;
  double dec;
  double BTmag;
  double VTmag;
  double gmag;
  double rmag;
  double gminusrmag;
  double gmagerr;
  double rmagerr;
  double gminusrmagerr;
  
} TYCHO2IMAGE,*PTYCHO2IMAGE;

typedef struct _sortindex {
  int gsc_bin_index;
  int record_number;
} SORTINDEX,*PSORTINDEX;

int IndexCompare(const void *first, const void *second)  {
  PSORTINDEX indexFirst = ((PSORTINDEX)first);
  PSORTINDEX indexSecond = ((PSORTINDEX)second);
  if (indexFirst->gsc_bin_index > indexSecond->gsc_bin_index) {
    return(1);
  } else if (indexFirst->gsc_bin_index < indexSecond->gsc_bin_index) {
    return(-1);
  } else {
    return(0);
  }
}


void FindStats(PTYCHO2IMAGE input_table,double *vector,int input_nrecs,char *label,int offset,int intFlag,int zeroFlag,double cutoff,int acceptNegative,int rmsoffset,double rmsCutoff)
{
  PTYCHO2IMAGE pInput;
  int index;
  int count;
  int result;
  int ignoreNegative = 0;
  char *intType[2] = {"double","int"};
  double med = 0;
  double rms = 0;
  double rmsValue;
  count = 0;
  
  if ((strcmp(label,"dec") == 0) || (strcmp(label,"gminusrmag") == 0)) {
    acceptNegative = 1;
    ignoreNegative = 1;
  }

  for (index = 0; index < input_nrecs; index++) {
    pInput = &input_table[index];

    if (intFlag == 0) {
      vector[count] = *((double *) ((char *)pInput + offset));
    } else {
      vector[count] = *((int *) ((char *)pInput + offset));
    }

    if (rmsoffset >= 0) {
      rmsValue = *((double *) ((char *)pInput + rmsoffset));
    } else {
      rmsValue = 0;
    }
    if ((acceptNegative == 0) &&
        (rmsValue < 0.0)) {
      continue;
    } else {
      
      if ((ignoreNegative == 0) && (rmsValue < 0)) {
        rmsValue = - rmsValue;
      }
    }
    
    if (rmsValue > rmsCutoff) {
      continue;
    }
    
    


    if (acceptNegative == 0) {
      if (vector[count] < 0.0) {
        continue;
      }
    } else {
      if ((ignoreNegative == 0) && (vector[count] < 0.0)){
        vector[count] = -vector[count];
      }
    }


    if (cutoff < 0) {
      count++;
    } else {
      if (vector[count] < cutoff) {        
        count++;
      }
    }
  }
  result = CalcMedianAndRMS(count,0,vector,&med,&rms,0,3.0,zeroFlag);
  if (zeroFlag) {
    /* Recover the median */
    if (result % 2) {
      med =  vector[(result/2)-1];
    } else {
      med = (vector[(result/2)-1] + vector[result/2])/2;
    }
  }

  printf("%13s %6s %10d %15.6f %15.6f %15.6f %15.6f\n",label,intType[intFlag],count,vector[0],med,vector[count-1],rms);

}


int main(int argc,char *argv[])
{
  char *argstr;
  char output_name[MAX_INPUT_NAME];
  char binary_name[MAX_INPUT_NAME];
  char index_name[MAX_INPUT_NAME];
  char* slashPtr;
  FILE* output_handle = NULL;
  File binary_handle = NULL;
  File index_handle = NULL;
  PSORTINDEX sort_index_table = NULL;
  PSORTINDEX pSortIndex;
  int sortIndexSize = 0;
  int curSortIndex;
  int curRecord = 0;

  int outputCount = 0;
  int colorReject = 0;
  int skipCount = 0;
  int skipCount2 = 0;
  int skipCount3 = 0;
  int zeromagCount = 0;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int linecounter = 0;
  int index1;
  int index2;
  char *charPtr;
  char cmdchar;
  int verbose = 0;
  int acceptNegative = 0;
  int doStats = 0;
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
  char *tycho2MeasurementText[TYCHO2_MAXMAG] = {"gmag","rmag"};
  double newrms;
  double scalerms;
  int intrms;
  int max_input_lines = 0;
  struct StarCat *starcat = NULL;
  struct Star starent;
  struct Star *star = &starent;
 

  int inputFile;


  int input_nrecs = 0;
  PTYCHO2IMAGE pInput;
  PTYCHO2IMAGE input_table = NULL;
  TYCHO2IMAGE input_record;
  GSCIMAGE gsc_record;
  PGSCIMAGE pGscImage = &gsc_record;

  int tempmag;


  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int tycho2Count[TYCHO2_MAXMAG];
  int tycho2Dual[TYCHO2_MAXMAG][TYCHO2_MAXMAG];
  int magcount[TYCHO2_MAXMAG][MAGBINS+1];
  int tycho2NoCQCount = 0;
  int tycho2HeaderCount = 0;
  int maxcqLength = 0;
  char maxcq[MAX_REF];
  double deltacolor;
  double maxdeltacolor = 0.0;
  int maxtycho2colorid;
  char inLine[MAX_BUFFER];
  int lineLen;
  char *inBuffer;
  int numLines;
  int nvals;
  double *vector = NULL;
  double rmsCutoff = 9.0;
  double value1;
  double rms1;
  double value2;
  double rms2;
  long long REFNumber;
  int skipFlag; /* Skip this flag because of neg rms or too high rms */
  int skipFlagBase;
  int noPositionCount = 0;
  int istar1;
  int istar2;
	char REF[MAX_REF];
	int refType;

  int gsc_bin_index;
  int cur_bin_index;
  int raBin;
  int decBin;
  int bin_count;
  int bin_index;
  int bin_list[MAX_ADJACENT_BINS+1];
  int writeItems;
  int readItems;


	STARINDEX tycho2Index;
	PSTARINDEX pTycho2Index = &tycho2Index;
  off_t curOffset;

  memset(tycho2Count,0,sizeof(tycho2Count));
  memset(tycho2Dual,0,sizeof(tycho2Dual));
  memset(magcount,0,sizeof(magcount));


  /* Loop through the arguments */
  output_name[0] = 0;


  memset(magcount,0,sizeof(magcount));
  
#ifdef SKIP_TRANSFORM
  printf("ERROR: SKIP_TRANSFORM IS SET\n");
#endif /* SKIP_TRANSFORM */
#ifdef USE_BV
  printf("WARNING: USE_BV is set\n");
#endif /* USE_BV */

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

        case 's':
        case 'S':
          doStats = 1;
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
    output_handle = fopen(output_name,"wt");
    if (output_handle == NULL) {
      errorFlag = 1;
      printf("Could not open output file %s\n",output_name);
    }
  }

  scalerms = rmsCutoff * sqrt(2.0);

  if (errorFlag) {
    printf("Usage: formattycho2 -o {output name} [-v][-n][-s][-c {cutoff}]\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s computes general statistics\n");
    printf("       where -c {cutoff} gives the RMS cutoff values\n");
    printf("       where -n accepts negative RMS values indicating only one measurement\n");

    return(-1);
  }

  printf("formattycho2 of %s %s Output Filename %s acceptNegative %d,rmsCutoff %f\n",
         __DATE__,__TIME__,output_name,acceptNegative,rmsCutoff);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of TYCHO2IMAGE is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(TYCHO2IMAGE),sizeof(GSCIMAGE));
 

  time(&startTime);
  
  strcpy(binary_name,output_name);
  slashPtr = strrchr(binary_name,'/');
  if (slashPtr == NULL) {
    binary_name[0] = 0;
  } else {
    slashPtr++;
    *slashPtr = 0;
  }
  strcpy(index_name,binary_name);
  strcat(binary_name,"tycho2.dat");
  strcat(index_name,"tycho2.idx");

  binary_handle = Open(binary_name,"w");
  if (binary_handle == NULL) {
    printf("Could not open sorted binary file %s\n",binary_name);
    exit(-1);
  }    
  index_handle = Open(index_name,"w");
  if (index_handle == NULL) {
    printf("Could not open index file %s\n",index_name);
    exit(-1);
  }

  
  fprintf(output_handle,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\n");
  fprintf(output_handle,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\n");

  /* The following line is a hack to get the catalog count correct */
  ty2regn(1,&istar1,&istar2,0);
  starcat = ty2open (1,100);
  if (starcat == NULL) {
		printf ("ERROR: failed to open the Tycho2 catalog\n");
		return (0);
  } else {
    printf("Tycho2 catalog opened with %d stars\n",starcat->nstars);
  }
  max_input_lines = starcat->nstars;
  ty2close(starcat);
  starcat = ty2open (1,max_input_lines);
  if (starcat == NULL) {
		printf ("ERROR: failed to open the Tycho2 catalog\n");
		return (0);
  }

  if (doStats) {
    input_table = (PTYCHO2IMAGE)calloc(max_input_lines,sizeof(TYCHO2IMAGE));
    if (input_table == NULL) {
      printf("ERROR: failed to allocate input table of size %d\n",(max_input_lines*sizeof(TYCHO2IMAGE)));
      exit(-1);
    }
    sort_index_table = (PSORTINDEX)calloc(max_input_lines*MAX_ADJACENT_BINS,sizeof(SORTINDEX));
    if (sort_index_table == NULL) {
      printf("ERROR: failed to allocate the sort index table of size %d\n",max_input_lines*MAX_ADJACENT_BINS);
      exit(-1);
    }
   

  }
  vector = (double *)calloc(max_input_lines,sizeof(double));
  if (vector == NULL) {
    printf("ERROR: failed to allocate vector of size %d\n",(max_input_lines*sizeof(TYCHO2IMAGE)));
    exit(-1);
  }

  input_nrecs = 0;
  numLines = 0;
  for (input_nrecs = 1; input_nrecs <= max_input_lines; input_nrecs++) {
    skipFlag = 0;
    skipFlagBase = 0;
    if (doStats) {
      pInput = &input_table[numLines];
    } else {
      pInput = &input_record;
    }
    memset(pInput,0,sizeof(TYCHO2IMAGE));
    result = ty2star(starcat,star,input_nrecs);
    if (result != 0) {
      if (result == 6) {
        noPositionCount++;
        continue;
      } else {
        printf("ERROR: ty2star returned %d for star %d\n",result,numLines);
        exit(-1);
      }
    }
    pInput->starid = (star->num * 100000.)+0.1;
#if 0
    printf("line %d starid %lld num %f\n",__LINE__,pInput->starid,star->num);
#endif
#if 0
    if (pInput->starid == 103541) {
      printf("At T%09lld at record %d\n",pInput->starid,input_nrecs);
    }
#endif
    pInput->ra = star->ra;
    pInput->dec = star->dec;
    pInput->BTmag =  star->xmag[0];
    if (pInput->BTmag == 0) {
      pInput->BTmag = 99.0;
      skipFlagBase = 1;
    }
    pInput->VTmag = star->xmag[1];
    if (pInput->VTmag == 0) {
      pInput->VTmag = 99.0;
      skipFlagBase = 1;
    }
    if (skipFlagBase == 0) {
#ifdef USE_BV
      pInput->gmag    = pInput->BTmag;
      pInput->rmag    = pInput->VTmag;
#else /* USE_BV */
      pInput->gmag    = 0.54*pInput->BTmag + 0.46*pInput->VTmag - 0.07;
      pInput->rmag    = -0.44*pInput->BTmag + 1.44*pInput->VTmag + 0.12;
#endif /* USE_BV */
      pInput->gmagerr = star->xmag[2]; /* Approximate! */
      pInput->rmagerr = star->xmag[3]; /* Approximate! */
      pInput->gminusrmag = pInput->gmag - pInput->rmag;
      pInput->gminusrmagerr = sqrt(sqr(pInput->gmagerr) + sqr(pInput->rmagerr));

    } else {
      pInput->gmag    = 99.0;
      pInput->gmagerr = 99.0;
      pInput->rmag    = 99.0;
      pInput->rmagerr = 99.0;
      pInput->gminusrmag = 99.0;
      pInput->gminusrmagerr = 99.0;

    }
    numLines++;

    tempmag = pInput->gmag-0.001;
    if (tempmag >= MAGBINS) {
      tempmag = MAGBINS;
    }

    for (index1 = 0; index1 < TYCHO2_MAXMAG; index1++) {
      switch(index1) {
      case TYCHO2_BTMAG:
        value1 = pInput->gmag;
        rms1   = pInput->gmagerr;
        break;
      case TYCHO2_VTMAG:
        value1 = pInput->rmag;
        rms1   = pInput->rmagerr;
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
        tycho2Count[index1]++;
      }
      for (index2 = 0; index2 < TYCHO2_MAXMAG; index2++) {
        skipFlag = skipFlagBase;
        switch(index2) {
        case TYCHO2_BTMAG:
          value2 = pInput->gmag;
          rms2   = pInput->gmagerr;
          break;
        case TYCHO2_VTMAG:
          value2 = pInput->rmag;
          rms2   = pInput->rmagerr;
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
        if ((index1 == TYCHO2_BTMAG) &&
            (index2 == TYCHO2_VTMAG)) {

          memset(pGscImage,0,sizeof(GSCIMAGE));
					
          sprintf(REF,"T%09lld",pInput->starid);
					if (GetREFNumber(REF,&pGscImage->REFNumber,&refType,0,0) != 0) {
						printf("ERROR: formattycho2 invalid REF %s\n",REF);
						exit(-1);
					}
#if 1					
          if (pGscImage->REFNumber == 5) {
            printf("ERROR: line %d invalid REFNumber %lld\n",__LINE__,pGscImage->REFNumber);
          }
#endif				
          pGscImage->ra = pInput->ra;
          pGscImage->dec = pInput->dec;
          pGscImage->Stdmag =  pInput->gmag;
          pGscImage->color = pInput->gmag - pInput->rmag;
          if ((pGscImage->color < MIN_APASS_COLOR) ||
              (pGscImage->color > MAX_APASS_COLOR)) {
            colorReject++;
            skipFlag = 1;
          }
          outputCount++;
          /* Covert the rms added in quadrature to an integer */
          newrms = (127.*(sqrt(sqr(pInput->gmagerr) + sqr(pInput->rmagerr))))/scalerms;
          if (newrms > 127) {
            newrms = 127;
#if 0
            if ( pGscImage->MAGFlag < 0) {
              printf("MAGFLAG (2) is %d for outputCount %d\n", pGscImage->MAGFlag,outputCount);
            }

#endif
                
          }
             
 
          intrms = newrms;

          pGscImage->MAGFlag = newrms;
#if 0
          if ( pGscImage->MAGFlag < 0) {
            printf("MAGFLAG (1) is %d for outputCount %d\n", pGscImage->MAGFlag,outputCount);
          }

#endif
          if ((skipFlag) ||(sqrt(sqr(pInput->gmagerr) + sqr(pInput->rmagerr)) > rmsCutoff))  {
#if 1
            pGscImage->VFlag = 1;
            skipCount++;

#else 
            /* Check the effects of not adding magnitudes in quadrature */
            if ((pInput->gmagerr > rmsCutoff) || (pInput->rmagerr > rmsCutoff)) {
              pGscImage->VFlag = 1;
              skipCount++;
            } else if (sqrt(sqr(pInput->gmagerr) + sqr(pInput->rmagerr)) > rmsCutoff) {
              pGscImage->VFlag = 2;
              skipCount2++;
            } else {
              skipCount3++;
            }
#endif
          }

#if 0
          if ((pGscImage->MAGFlag > 90)  && (pGscImage->MAGFlag == 0)) {
            printf("at ref %s\n",pGscImage->REF);
          }
#endif
#ifdef SKIP_TRANSFORM
          fprintf(output_handle,"%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\n",
                  pGscImage->REFNumber,
                  pGscImage->ra,
                  pGscImage->dec,
                  pInput->BTmag,
                  pInput->VTmag,
                  pGscImage->class,
                  pGscImage->VFlag,
                  pGscImage->MAGFlag,
                  pGscImage->RaPM,
                  pGscImage->DecPM);
#else /* SKIP_TRANSFORM */
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
#endif /* SKIP_TRANSFORM */
          
        }
        if (skipFlag == 0) {
          tycho2Dual[index1][index2]++;
        }
      } 
    }

    if (numLines > max_input_lines) {
      printf("ERROR: table too small at %d records\n",input_nrecs);
      break;
    }
    if (verbose && ((input_nrecs % 100000) == 0)) {
      printf("Reading record %d\n",input_nrecs);
    }
  }
        
  printf("Input lines %d, input records %d at %d\n",numLines,input_nrecs-1);
  ty2close(starcat);
  if (doStats) {
    time(&curTime);
    curTime -= startTime;
    printf("Preparing sort index at %d seconds\n",curTime);
    for (curRecord = 0; curRecord < numLines; curRecord++) {
      pInput = &input_table[curRecord];
      gsc_bin_index = GetGSCBin(pGscBin,pInput->ra,pInput->dec,&decBin,&raBin,"formatsdss");
      FindAdjacentBins(pGscBin,gsc_bin_index,bin_list,&bin_count);
      bin_list[bin_count] = gsc_bin_index;
      bin_count++;
      for (bin_index = 0; bin_index < bin_count; bin_index++) {
        pSortIndex = &sort_index_table[sortIndexSize];
        pSortIndex->gsc_bin_index = bin_list[bin_index];
        pSortIndex->record_number = curRecord;
        sortIndexSize++;
        if (sortIndexSize >= numLines*MAX_ADJACENT_BINS) {
          printf("ERROR: sortIndexSize %d exceeds %d\n",sortIndexSize,numLines*MAX_ADJACENT_BINS);
          exit(-1);
        }
      }
    }
    time(&curTime);
    curTime -= startTime;
    printf("Sorting %d items at %d seconds\n",sortIndexSize,curTime);
    qsort((void *)sort_index_table,sortIndexSize,sizeof(SORTINDEX),IndexCompare);
    time(&curTime);
    curTime -= startTime;
    printf("Writing %d items at %d seconds\n",sortIndexSize,curTime);
  
    curOffset = 0;
    curSortIndex = 0;
    pSortIndex = &sort_index_table[curSortIndex];
    gsc_bin_index = pSortIndex->gsc_bin_index;
    for (cur_bin_index = 0; cur_bin_index < pGscBin->total_gsc_bins; cur_bin_index++) {
      pTycho2Index->offset = curOffset;
      pTycho2Index->binNumber = cur_bin_index;
      pTycho2Index->numStars = 0;
      while (1) {
        if ((cur_bin_index < gsc_bin_index) ||
            (cur_bin_index >= pGscBin->total_gsc_bins)) {
          break;
        }
        if (cur_bin_index == gsc_bin_index) {
          pInput = &input_table[pSortIndex->record_number];
          memset(pGscImage,0,sizeof(GSCIMAGE));
          sprintf(REF,"T%09lld",pInput->starid);
					if (GetREFNumber(REF,&pGscImage->REFNumber,&refType,0,0) != 0) {
						printf("ERROR: formattycho2 (2) invalid REF %s\n",REF);
						exit(-1);
					}
#if 0					
          if (pGscImage->REFNumber == 5) {
            printf("ERROR: line %d invalid REFNumber %lld\n",__LINE__,pGscImage->REFNumber);
          }
#endif

#if 0
        if (pGscImage->REFNumber == 59502000411L) {
          printf("At REF %lld\n",pGscImage->REFNumber);
        }
#endif

          pGscImage->ra = pInput->ra;
          pGscImage->dec = pInput->dec;
          pGscImage->Stdmag = pInput->BTmag;
          pGscImage->color = pInput->BTmag - pInput->VTmag;
          


          writeItems = Write(binary_handle,pGscImage,sizeof(GSCIMAGE),1);
          if (writeItems != 1) {
            printf("ERROR: writing binary file\n");
            exit(-1);
          }
          curOffset += sizeof(GSCIMAGE);
          pTycho2Index->numStars++;
          curSortIndex++;
          if (curSortIndex < sortIndexSize) {
            pSortIndex = &sort_index_table[curSortIndex];
            gsc_bin_index =  pSortIndex->gsc_bin_index;
          } else {
            gsc_bin_index =  pGscBin->total_gsc_bins;
            curRecord = numLines;
            break;
          }
        }
      }
      writeItems = Write(index_handle,pTycho2Index,sizeof(STARINDEX),1);
      if (writeItems != 1) {
        printf("ERROR writing index file\n");
        exit(-1);
      }
    }
  }


  time(&curTime);
  curTime -= startTime;
  printf("Done writing index at %d seconds\n",curTime);
  


  if (doStats) {
    printf("     variable  type      count           min             med             max             rms\n");
    FindStats(input_table,vector,numLines,"starid"       ,offsetof(TYCHO2IMAGE,starid)       ,1,0,-1.0     ,acceptNegative,-1,9.0);
    FindStats(input_table,vector,numLines,"ra"           ,offsetof(TYCHO2IMAGE,ra)           ,0,0,-1.0     ,acceptNegative,-1,9.0);
    FindStats(input_table,vector,numLines,"dec"          ,offsetof(TYCHO2IMAGE,dec)          ,0,0,-1.0     ,acceptNegative,-1,9.0);
    FindStats(input_table,vector,numLines,"BTmag"        ,offsetof(TYCHO2IMAGE,BTmag)        ,0,0,90.0     ,acceptNegative,offsetof(TYCHO2IMAGE,gmagerr),rmsCutoff);
    FindStats(input_table,vector,numLines,"VTmag"        ,offsetof(TYCHO2IMAGE,VTmag)        ,0,0,90.0     ,acceptNegative,offsetof(TYCHO2IMAGE,rmagerr),rmsCutoff);
    FindStats(input_table,vector,numLines,"gmag"         ,offsetof(TYCHO2IMAGE,gmag)         ,0,0,90.0     ,acceptNegative,offsetof(TYCHO2IMAGE,gmagerr),rmsCutoff);
    FindStats(input_table,vector,numLines,"rmag"         ,offsetof(TYCHO2IMAGE,rmag)         ,0,0,90.0     ,acceptNegative,offsetof(TYCHO2IMAGE,rmagerr),rmsCutoff);
    FindStats(input_table,vector,numLines,"gmagerr"      ,offsetof(TYCHO2IMAGE,gmagerr)      ,0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);
    FindStats(input_table,vector,numLines,"rmagerr"      ,offsetof(TYCHO2IMAGE,rmagerr)      ,0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);

    FindStats(input_table,vector,numLines,"gminusrmag"   ,offsetof(TYCHO2IMAGE,gminusrmag)   ,0,0,90.0     ,acceptNegative,offsetof(TYCHO2IMAGE,gminusrmagerr),rmsCutoff);
    FindStats(input_table,vector,numLines,"gminusrmagerr",offsetof(TYCHO2IMAGE,gminusrmagerr),0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);
  }
  fclose(output_handle);
  if (doStats) {
    Close(binary_handle);
    Close(index_handle);

  }


  time(&curTime);
  curTime -= startTime;

  printf("Measurement:");
  for (index1 = 0; index1 < TYCHO2_MAXMAG; index1++) {
    printf("%10s",tycho2MeasurementText[index1]);
  }
  printf("\n");

  printf("Count:      ");
  for (index1 = 0; index1 < TYCHO2_MAXMAG; index1++) {
    printf("%10d",tycho2Count[index1]);
  }
  printf("\n");

 
  for (index1 = 0; index1 < TYCHO2_MAXMAG; index1++) {
    printf("%10s  ",tycho2MeasurementText[index1]);
    for (index2 = 0; index2 < TYCHO2_MAXMAG; index2++) {
      printf("%10d",tycho2Dual[index1][index2]);
    }
    printf("\n");
  }





  for (tempmag = 0; tempmag <= MAGBINS; tempmag++) {
    printf("rmag %2d     ",tempmag);
    for (index1 = 0; index1 < TYCHO2_MAXMAG; index1++) {
      printf("%10d",magcount[index1][tempmag]);
    }
    printf("\n");
  }    
  if (sort_index_table != NULL) {
    free(sort_index_table);
  }

  if (input_table != NULL) {
    free(input_table);
  }

  if (vector != NULL) {
    free(vector);
  }

  time(&curTime);
  curTime -= startTime;

  printf("Execution Time: %d seconds; noPosition %d rejected %d %d %d rejected because of color %d stars written %d\n",curTime,noPositionCount,skipCount,skipCount2,skipCount3,colorReject,outputCount);

  return(EXIT_SUCCESS);
}
