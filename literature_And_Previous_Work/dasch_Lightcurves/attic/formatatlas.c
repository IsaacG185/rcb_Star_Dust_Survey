// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatatlas.c
 *
 *   Catalog comes from https://archive.stsci.edu/prepds/atlas-refcat2/
 *   See /dasch/data/scanner/backup/2018_09_27/atlas/README.txt for catalog format
 *   See /dasch/data/scanner/backup/2018_10_04/ATLAS1809.09157.pdf
 *   http://mastweb.stsci.edu/ps1casjobs/  Pan-STARRS DR1 at the STScI MAST archive
 *   http://www.fallingstar.com
 *   min_pmdec -5818 mas or 1/41 degrees.  Need at least two GSC bins resident
 *
 * Read the atlas input catalog and put it into the gsc232bin.db format
 *        Accept everything:               formatatlas -s -v -n -o /dasch/Pipeline/catalogs/apass.db
 *
 echo "formatatlas -v -c 0.1 -o /dasch/Pipeline/catalogs/atlas_temp.db >& /dasch/raid020/junk/atlas.log" | at now

 formatgsc  -q atlas /dasch/Pipeline/catalogs/atlas_temp.db

 *
 *
 cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  formatatlas.c pipelineutils.a -L /dasch/install/lib -lm   -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -o formatatlas

 *  1.  Use this program to reformat the apass catalog
 *        Accept nobs > 2, rms cutoff 0.1: formatatlas -v -c 0.1 -o /dasch/Pipeline/catalogs/apass_temp.db
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
 *     matchcatalogs of Mar 23 2011 17:19:47, 
 *       matching /dasch/Pipeline/catalogs/ucac3.dat and /dasch/Pipeline/catalogs/apass_temp.dat 
 *       output /dasch/Pipeline/catalogs/matchgsckepler.db total bins 168966386 TOLERANCE 0.000556 degrees
 *      Done  populatedBinCount 6121155 matches 6292508 for 7173325 kepler stars at 699 seconds
 *            median degDrad is 0.000085 deg or 0.3 arcsec max is 0.000556 or 2.0 arcsec.
 *      mv merge.dat apass.dat
 *      mv merge.idx apass.idx
 *      rm apass_temp*
 *
 * Oct 15, 2018 Edward J. Los - Initial versoin, adapted from formatatlas.c
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
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER 1024
#define MAGBINS  30  /* Number of magnitude bins */
#define ATLAS_BP     0
#define ATLAS_RP     1
#define ATLAS_G      2
#define ATLAS_R      3
#define ATLAS_I      4
#define ATLAS_Z      5
#define ATLAS_MAXMAG 6
#define MAX_MODULUS  1 /* Base 10 modulus */
#define MODULUS_BASE 10
#define MAX_SKIP_BIN 1000 /* dec bins skipped between entries */
#define SKIP_BUFFER 4 /* extra lines to read before the skip point */
#define ATLAS_MAXRMS 1.0
#define RMSBINS 100
#define RESIDENT_DEC_BINS 128
#define MAX_DEC_BIN 232860 /* dec_bin 3903 at Declination  -29.016 */
#define MAX_INPUT_LINES 992637834
#define GAIADR2_EPOCH 2015.5
/* #define PLOT_RMS_GMAG */ /* starbase count of gmag vs rms */
extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  
#define VARIABLE 1      /* see dupvar below */
#define NOT_AVAILABLE 2 /* see dupvar below */
#define DUPLICATE 4     /* see dupvar below */
typedef struct _atlasimage {
  long long objid; /* (none) Object ID           */
  double RA; /* (degrees) Right ascension from Gaia DR2, J2000, epoch 2015.5     */
  double Dec; /* (degrees) Declination from Gaia DR2, J2000, epoch 2015.5      */
  double plx; /* (mas) Parallax from Gaia DR2         */
  double dplx; /* (mas) Parallax uncertainty from Gaia DR2        */
  double pmra; /* (mas/yr) Proper motion in right ascension from Gaia DR2     */
  double dpmra; /* (mas/yr) Proper motion uncertainty in right ascension       */
  double pmdec; /* (mas/yr) Proper motion in declination from Gaia DR2      */
  double dpmdec; /* (mas/yr) Proper motion uncertainty in declination        */
  double Gaia; /* (mag) Gaia G magnitude          */
  double dGaia; /* (mag) Gaia G magnitude uncertainty         */
  double BP; /* (mag) Gaia G_bp magnitude          */
  double dBP; /* (mag) Gaia G_bp magnitude uncertainty         */
  double RP; /* (mag) Gaia G_rp magnitude          */
  double dRP; /* (mag) Gaia G_rp magnitude uncertainty         */
  int Teff; /* (K) Gaia stellar effective temperature         */
  double AGaia; /* (mag) Gaia estimate of G-band extinction for this star     */
  int dupvar; /* (none) Gaia variability and duplicate flags, 0/1/2 for "CONSTANT"/"VARIABLE"/"NOT AVAILABLE" + 4*DUPLICATE  */
  double Ag; /* (mag) SFD estimate of total g-band extinction       */
  double rp1; /* (arcsec) Radius where cummulative G flux exceeds 0.1 x this star   */
  double r1; /* (arcsec) Radius where cummulative G flux exceeds 1.0 x this star   */
  double r10; /* (arcsec) Radius where cummulative G flux exceeds 10.0 x this star   */
  double g; /* (mag) PanSTARRS g magnitude          */
  double dg; /* (mag) PanSTARRS g magnitude uncertainty         */
  double gchi; /* (none) chi^2 / DOF for contributors        */
  int gcontrib; /* (none) Bitmap of conributing catalogs to g       */
  double r; /* (mag) PanSTARRS r magnitude          */
  double dr; /* (mag) PanSTARRS r magnitude uncertainty         */
  double rchi; /* (none) chi^2 / DOF for contributors        */
  int rcontrib; /* (none) Bitmap of conributing catalogs to r       */
  double i; /* (mag) PanSTARRS i magnitude          */
  double di; /* (mag) PanSTARRS i magnitude uncertainty         */
  double ichi; /* (none) chi^2 / DOF for contributors        */
  int icontrib; /* (none) Bitmap of conributing catalogs to i       */
  double z; /* (mag) PanSTARRS z magnitude          */
  double dz; /* (mag) PanSTARRS z magnitude uncertainty         */
  double zchi; /* (none) chi^2 / DOF for contributors        */
  int zcontrib; /* (none) Bitmap of conributing catalogs to z       */
  int nstat; /* (none) Count of griz outliers rejected        */
  double J; /* (mag) 2MASS J magnitude          */
  double dJ; /* (mag) 2MASS J magnitude uncertainty         */
  double H; /* (mag) 2MASS H magnitude          */
  double dH; /* (mag) 2MASS H magnitude uncertainty         */
  double K; /* (mag) 2MASS K magnitude          */
  double dK; /* (mag) 2MASS K magnitude uncertainty         */  
} ATLASIMAGE,*PATLASIMAGE;

typedef struct _gscatlasimage {
	long long REFNumber; /* Encoded catalog reference number */
  long long objid; /* (none) Object ID           */
  int atlasNumber;     /* This is just the place in the .csv file set */
  double ra;           /* Right Ascension in degrees */
  double dec;          /* Declination in degrees */
  float Stdmag;        /* Blue magnitude */
  float color;         /* color (bmag - vmag) */
  float RaPM;          /* Right ascension proper motion in mas/yr */
  float DecPM;         /* Declination proper motion in mas/yr */
	float RaSigmaPM;     /* Error of the Right Ascension proper motion */
	float DecSigmaPM;    /* Error of the Declination proper motion */
  float rmag;          /* SDSS r color */
  float imag;          /* SDSS i color */
  float zmag;          /* SDSS z color */
  float gmagerr;       /* magnitude errors */
  float rmagerr;
  float imagerr;
  float zmagerr;
  char class;          /* class */
  char VFlag;          /* variable flag */
  char MAGFlag;        /* magnitude flag */
  char flag;           /* padding - used by check_photometry and formatgsc */
} GSCATLASIMAGE,*PGSCATLASIMAGE;

typedef struct _fileTable {
  const char *fileName;
  int inputLines;
  int skipLines;
} FILETABLE,*PFILETABLE;

FILETABLE fileTable[] = {
  {"/dasch/raid015/atlas_refcat2/hlsp_atlas-refcat2_atlas_ccd_m90-m53_multi_v1_cat.csv",191518323,191493646}, /*  3 skips < -2000;  -321 to +325; 3 skips > 2000 */
  {"/dasch/raid015/atlas_refcat2/hlsp_atlas-refcat2_atlas_ccd_m53-m33_multi_v1_cat.csv",190200973,190104958}, /*  2 skips < -1200; -18 to +18; 2 skips > 1200 */
  {"/dasch/raid015/atlas_refcat2/hlsp_atlas-refcat2_atlas_ccd_m33-m15_multi_v1_cat.csv",209919540,209796490}, /*  2 skips < -1100; -15 to +16; 3 skips > 1100 */
  {"/dasch/raid015/atlas_refcat2/hlsp_atlas-refcat2_atlas_ccd_m15-p19_multi_v1_cat.csv",200250764,200181514}, /*  2 skips < -2100; -25 to +25; 2 skips > 2100 */ 
  {"/dasch/raid015/atlas_refcat2/hlsp_atlas-refcat2_atlas_ccd_p19-p90_multi_v1_cat.csv",200748234,200680156}, /*  1 skip < -4300; -341 to +354; 11 skips > 1100 */  
};
int fileTableSize = sizeof(fileTable)/sizeof(FILETABLE);

/* Sort the records in decreasing declination */
int ATLASCompare(const void *first, const void *second) 
{
  double numberFirst = ((PGSCATLASIMAGE )first)->dec;
  double numberSecond = ((PGSCATLASIMAGE )second)->dec;
  if (numberFirst > numberSecond) {
    return(-1);
  } else if (numberFirst < numberSecond) {
    return(1);
  } else {
    return(0);
  }
}

void FlushOutputTable(FILE* output_handle,PGSCATLASIMAGE pGscAtlasTable,int maxAllocOutput,int max_read_dec_bin,int* pCurOutputCount,int *pOutputWrittenCount,int *poutput_dec_bin,int *decBinCountsOutput,int finalFlag)
{
  int curOutputCount = *pCurOutputCount;
  int outputWrittenCount = *pOutputWrittenCount;
  int output_dec_bin = *poutput_dec_bin;
  int initial_dec_bin = output_dec_bin;
  int curIndex;
  int newOutputCount = 0;
  int dec_bin;
  
  PGSCATLASIMAGE pGscAtlasImage = NULL;  
  if (curOutputCount == 0) {
    return;
  }
  qsort((void*)pGscAtlasTable,curOutputCount,sizeof(GSCATLASIMAGE),ATLASCompare);
  for (curIndex = (curOutputCount-1); curIndex >= 0; curIndex--) {
    pGscAtlasImage = &pGscAtlasTable[curIndex];
    dec_bin = GetDecBin(pGscBin,pGscAtlasImage->dec,"formatatlas");
    if (dec_bin <= initial_dec_bin) {
      printf("ERROR: line %d sort error dec_bin %d output_dec_bin %d, max_read_dec_bin %d\n",__LINE__,dec_bin,output_dec_bin,max_read_dec_bin);
    }
    if ((finalFlag == 0) && 
        (dec_bin > (max_read_dec_bin - RESIDENT_DEC_BINS))) { 
      break;
    }
    output_dec_bin = dec_bin;
    newOutputCount++;
    decBinCountsOutput[dec_bin]++;
    fprintf(output_handle,"%lld\t%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
            pGscAtlasImage->REFNumber,
            pGscAtlasImage->objid,
            pGscAtlasImage->ra,
            pGscAtlasImage->dec,
            pGscAtlasImage->Stdmag,
            pGscAtlasImage->color,
            pGscAtlasImage->class,
            pGscAtlasImage->VFlag,
            pGscAtlasImage->MAGFlag,
            pGscAtlasImage->RaPM,
            pGscAtlasImage->DecPM,
            pGscAtlasImage->RaSigmaPM,
            pGscAtlasImage->DecSigmaPM,
            pGscAtlasImage->rmag,
            pGscAtlasImage->imag,
            pGscAtlasImage->zmag,
            pGscAtlasImage->gmagerr,
            pGscAtlasImage->rmagerr,
            pGscAtlasImage->imagerr,
            pGscAtlasImage->zmagerr);
  }
  if (newOutputCount == 0) {
    printf("ERROR: line %d failed to flush the buffer\n",__LINE__);
    exit(-1);
  }
  printf("FlushOutputTable wrote %d records from dec bin %d to %d curIndex %d\n",newOutputCount,initial_dec_bin,dec_bin-1,curIndex);
  *poutput_dec_bin = dec_bin -1;
  *pOutputWrittenCount = outputWrittenCount+newOutputCount;
  *pCurOutputCount = curIndex+1;
}

int main(int argc,char *argv[])
{
  char *argstr;
  char output_name[MAX_INPUT_NAME];
  char raw_name[MAX_INPUT_NAME];
  FILE* output_handle = NULL;
  FILE* raw_handle = NULL;
  int outputCount = 0;
  int colorReject = 0;
  int skipCount = 0;
  int rawSkipCount = 0; /* Count of raw file output skipped */
  int zeromagCount = 0;
  int rawSkipFlag; /* entry skipped for raw file output */
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
  char *dotPtr = NULL;
  int compareResult;
  int result;
  int binNumber;
  int oldBinNumber = -1;
  int dec_bin;
  int old_dec_bin;
  double oldDec = -999;
  int skip_bin_index;
  int skip_bin_size;
  int oldDecBin = 0;
  int savedDecBin = -2;
  int ra_bin;
  int maxDecBin = 0;
  double maxRaError = 0;
  double maxDecError = 0;
  int maxWriteCount = 0;
  int totalBinCount = 0;
  FILE *input_handle;
  char *apassMeasurementText[ATLAS_MAXMAG] = {"bp","rp","g","r","i","z"};
  double newrms;
  double scalerms = 0;
  double maxrms2 = 0;
  int intrms;
  int intrms2;
  long long max_objid;
  long long min_objid;
  unsigned long long objidbits[MAX_MODULUS]; 
  unsigned long long modobjidbits[MAX_MODULUS];
  int modindex;
  unsigned long long temp_objid;
  unsigned long long modulus;
  int old_bin_index = -1;
  PFILETABLE pFileTable;
  double factor;
  double dec2;
  double ra2;
  int flipFlag; /* Require 12 hour telescope flip when crossing the pole */
  int dec_bin2;

  memset(objidbits,0,sizeof(objidbits)); 
  memset(modobjidbits,0,sizeof(modobjidbits));




  int inputFile;

  int input_nrecs = 0;
  ATLASIMAGE input_record;
  PATLASIMAGE pInput = &input_record;
  PGSCATLASIMAGE pGscAtlasImage = NULL;
  PGSCATLASIMAGE pGscAtlasTable = NULL;
  int passCounter;

  int tempmag;


  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int apassCount[ATLAS_MAXMAG];
  int apassDual[ATLAS_MAXMAG][ATLAS_MAXMAG];
  int magcount[ATLAS_MAXMAG][MAGBINS+5];
  int rmscount[RMSBINS+5];
  int rmscountTotal[RMSBINS+5][MAGBINS+5];
  
  int apassNoCQCount = 0;
  int apassHeaderCount = 0;
  int maxcqLength = 0;
  char maxcq[MAX_REF];
  double deltacolor;
  double maxdeltacolor = 0.0;
  int maxapasscolorid;
  char inLine[MAX_BUFFER];
  int lineLen;
  int maxLineLen = 0;
  char *inBuffer;
  int numLines;
  int numFileLines;
  int nvals;
  double rmsCutoff = ATLAS_MAXRMS;
  double value1;
  double rms1;
  double value2;
  double rms2;
  long long REFNumber;
  int skipFlag; /* Skip this entry because of zero mag or too high rms */
  int skipFlagBase;  /* Applies to the first of a magnitude pair */
  int zeroBmagCount = 0;
  int zeroVmagCount = 0;
  int negZeroCount = 0;
  int posZeroCount = 0;
  double maxpmdec = 0;
  double minpmdec = 0;
  double mindec = 0;
  double maxdec = 0;
  double minra = 0;
  double maxra = 0;
  int *decBinCounts = NULL;
  int *decBinCountsInitial = NULL;
  int *decBinCountsFinal = NULL;
  int *decBinCountsOutput = NULL;
  int *skipBinCounts = NULL;
  int *skipBinCountsTotal = NULL;
  int maxDecBinCount = 0;
  int maxDecBinValue = -1;
  int expectedRecords = 0;
  int outputIndex;
  int curOutputCount = 0;
  int maxAllocOutput = 0;
  int outputWrittenCount = 0;
  long long maxAllocBytes;
  int maxSkipLines = 0;
  int skipLines;
  int allocOutputLines;
  int output_dec_bin = -1; /* Last bin to be written */
  int max_read_dec_bin = -1;  /* Highest bin read in pass 1 */ 

  memset(apassCount,0,sizeof(apassCount));
  memset(apassDual,0,sizeof(apassDual));
  memset(magcount,0,sizeof(magcount));
  memset(rmscount,0,sizeof(rmscount));
  memset(rmscountTotal,0,sizeof(rmscountTotal));
  decBinCounts = (int *)calloc(pGscBin->dec_bins,sizeof(int));
  decBinCountsInitial = (int *)calloc(pGscBin->dec_bins,sizeof(int));
  decBinCountsFinal = (int *)calloc(pGscBin->dec_bins,sizeof(int));
  decBinCountsOutput = (int *)calloc(pGscBin->dec_bins,sizeof(int));
  skipBinCounts = (int *)calloc((2*pGscBin->dec_bins)+2,sizeof(int));
  skipBinCountsTotal = (int *)calloc((2*pGscBin->dec_bins)+2,sizeof(int));

  /* Loop through the arguments */
  output_name[0] = 0;
	raw_name[0] = 0;

  memset(magcount,0,sizeof(magcount));
  

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
    output_handle = fopen(output_name,"wt");
    if (output_handle == NULL) {
      errorFlag = 1;
      printf("Could not open output file %s\n",output_name);
    }
  }

  if (raw_name[0] != 0) {
		raw_handle = fopen(raw_name,"wt");
    if (raw_handle == NULL) {
      errorFlag = 1;
      printf("Could not open raw file %s\n",raw_name);
    } else {
			fprintf(raw_handle,"ra\tdec\tStdmag\trms\n");
			fprintf(raw_handle,"--\t---\t------\t---\n");

		}
  }

  scalerms = rmsCutoff * sqrt(2.0);
  maxrms2 = sqr(scalerms);

  if (errorFlag) {
    printf("Usage: formatatlas -o {output name} [-r <raw file>] [-v][-n][-s][-c {cutoff}]\n");
		printf("       where -r <filename> is a raw table of RA and DEC for plotting \n");
    printf("       where -v is the verbose flag\n");
    printf("       where -c {cutoff} gives the RMS cutoff values\n");

    return(-1);
  }

  printf("formatatlas of %s %s Output Filename %s Raw filename %s,rmsCutoff %f\n",
         __DATE__,__TIME__,output_name,raw_name,rmsCutoff);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of ATLASIMAGE is %d Size of GSCATLASIMAGE is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(ATLASIMAGE),sizeof(GSCATLASIMAGE),sizeof(GSCIMAGE));
 

  time(&startTime);
  
  
  fprintf(output_handle,"REFNumber\tobjid\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tRaSigmaPM\tDecSigmaPM\trmag\timag\tzmag\tgmagerr\trmagerr\timagerr\tzmagerr\n");
  fprintf(output_handle,"---------\t-----\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t---------\t----------\t----\t----\t----\t-------\t-------\t-------\t-------\n");


  for (inputFile = 0; inputFile < fileTableSize; inputFile++) {
    pFileTable = &fileTable[inputFile];
    expectedRecords += pFileTable->inputLines;
    skipLines = pFileTable->inputLines-pFileTable->skipLines-SKIP_BUFFER;
    if (skipLines > maxSkipLines) {
      maxSkipLines = skipLines;
    }
    allocOutputLines = (RESIDENT_DEC_BINS*MAX_DEC_BIN) + skipLines;
    if (allocOutputLines > maxAllocOutput) {
      maxAllocOutput = allocOutputLines;
    }
  }
  printf("expectedRecords %d maxSkipLines %d maxAllocOutput %lld\n",expectedRecords,maxSkipLines,maxAllocOutput);
  maxAllocOutput += maxSkipLines;
  maxAllocBytes = maxAllocOutput;
  maxAllocBytes *= sizeof(GSCATLASIMAGE);
  /* The following line shows that an attempt to sort the file internally is not feasable */
  printf("Need to allocate %d records, %lld bytes\n",maxAllocOutput,maxAllocBytes);
  pGscAtlasTable = (PGSCATLASIMAGE)calloc(maxAllocOutput,sizeof(GSCATLASIMAGE));
  if (pGscAtlasTable == NULL) {
    printf("ERROR: failed to allocate pGscAtlasTable\n");
    exit(-1);
  }


  input_nrecs = 0;
  for (inputFile = 0; inputFile < fileTableSize; inputFile++) {
  
    pFileTable = &fileTable[inputFile];
    
    printf("Opening %s\n",pFileTable->fileName);
    numFileLines = 0;
    for (passCounter = 0; passCounter < 2; passCounter++) {
      input_handle = fopen(pFileTable->fileName,"rt");
      if (input_handle == NULL) {
        printf("ERROR: failed to open %s\n",pFileTable->fileName);
        exit(-1);
      }
      numLines = 0;
      while(1) {
        memset(pInput,0,sizeof(ATLASIMAGE));
        inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
        if (inBuffer == NULL) {
          break;
        }
        numLines++;
        if (passCounter == 0) {
          /* In the first pass, read the skip lines */
          if (numLines < (pFileTable->skipLines-SKIP_BUFFER)) {
            continue;
          }
        } else {
          /* In the second pass, read up to the skip lines */
          if (numLines >= (pFileTable->skipLines-SKIP_BUFFER)) {
            continue;
          }
        }
        numFileLines++;
        lineLen = strlen(inBuffer);
        if (lineLen > maxLineLen) {
          maxLineLen = lineLen;
        }
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

        nvals = sscanf(inBuffer,"%lld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%lf,%d,%lf,%lf,%lf,%d,%lf,%lf,%lf,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf",
                       &pInput->objid, /* (none) Object ID           */
                       &pInput->RA, /* (degrees) Right ascension from Gaia DR2, J2000, epoch 2015.5     */
                       &pInput->Dec, /* (degrees) Declination from Gaia DR2, J2000, epoch 2015.5      */
                       &pInput->plx, /* (mas) Parallax from Gaia DR2         */
                       &pInput->dplx, /* (mas) Parallax uncertainty from Gaia DR2        */
                       &pInput->pmra, /* (mas/yr) Proper motion in right ascension from Gaia DR2     */
                       &pInput->dpmra, /* (mas/yr) Proper motion uncertainty in right ascension       */
                       &pInput->pmdec, /* (mas/yr) Proper motion in declination from Gaia DR2      */
                       &pInput->dpmdec, /* (mas/yr) Proper motion uncertainty in declination        */
                       &pInput->Gaia, /* (mag) Gaia G magnitude          */
                       &pInput->dGaia, /* (mag) Gaia G magnitude uncertainty         */
                       &pInput->BP, /* (mag) Gaia G_bp magnitude          */
                       &pInput->dBP, /* (mag) Gaia G_bp magnitude uncertainty         */
                       &pInput->RP, /* (mag) Gaia G_rp magnitude          */
                       &pInput->dRP, /* (mag) Gaia G_rp magnitude uncertainty         */
                       &pInput->Teff, /* (K) Gaia stellar effective temperature         */
                       &pInput->AGaia, /* (mag) Gaia estimate of G-band extinction for this star     */
                       &pInput->dupvar, /* (none) Gaia variability and duplicate flags, 0/1/2 for "CONSTANT"/"VARIABLE"/"NOT AVAILABLE" + 4*DUPLICATE  */
                       &pInput->Ag, /* (mag) SFD estimate of total g-band extinction       */
                       &pInput->rp1, /* (arcsec) Radius where cummulative G flux exceeds 0.1 x this star   */
                       &pInput->r1, /* (arcsec) Radius where cummulative G flux exceeds 1.0 x this star   */
                       &pInput->r10, /* (arcsec) Radius where cummulative G flux exceeds 10.0 x this star   */
                       &pInput->g, /* (mag) PanSTARRS g magnitude          */
                       &pInput->dg, /* (mag) PanSTARRS g magnitude uncertainty         */
                       &pInput->gchi, /* (none) chi^2 / DOF for contributors        */
                       &pInput->gcontrib, /* (none) Bitmap of conributing catalogs to g       */
                       &pInput->r, /* (mag) PanSTARRS r magnitude          */
                       &pInput->dr, /* (mag) PanSTARRS r magnitude uncertainty         */
                       &pInput->rchi, /* (none) chi^2 / DOF for contributors        */
                       &pInput->rcontrib, /* (none) Bitmap of conributing catalogs to r       */
                       &pInput->i, /* (mag) PanSTARRS i magnitude          */
                       &pInput->di, /* (mag) PanSTARRS i magnitude uncertainty         */
                       &pInput->ichi, /* (none) chi^2 / DOF for contributors        */
                       &pInput->icontrib, /* (none) Bitmap of conributing catalogs to i       */
                       &pInput->z, /* (mag) PanSTARRS z magnitude          */
                       &pInput->dz, /* (mag) PanSTARRS z magnitude uncertainty         */
                       &pInput->zchi, /* (none) chi^2 / DOF for contributors        */
                       &pInput->zcontrib, /* (none) Bitmap of conributing catalogs to z       */
                       &pInput->nstat, /* (none) Count of griz outliers rejected        */
                       &pInput->J, /* (mag) 2MASS J magnitude          */
                       &pInput->dJ, /* (mag) 2MASS J magnitude uncertainty         */
                       &pInput->H, /* (mag) 2MASS H magnitude          */
                       &pInput->dH, /* (mag) 2MASS H magnitude uncertainty         */
                       &pInput->K, /* (mag) 2MASS K magnitude          */
                       &pInput->dK); /* (mag) 2MASS K magnitude uncertainty         */  
        if (nvals != 45) {
          printf("ERROR: nvals is %d in line %d\n",nvals,numLines);
        } else {
          if (numFileLines == 1) {
            mindec = pInput->Dec;
            maxdec = pInput->Dec;
            minra = pInput->RA;
            maxra = pInput->RA;
          } else {
            if (pInput->Dec > maxdec) {
              maxdec = pInput->Dec;
            }
            if (pInput->Dec < mindec) {
              mindec = pInput->Dec;
            }
         
            if (pInput->RA > maxra) {
              maxra = pInput->RA;
            }
            if (pInput->RA < minra) {
              minra = pInput->RA;
            }
          }
          dec_bin = GetDecBin(pGscBin,pInput->Dec,"formatatlas2");
          if (passCounter != 0) {
            if (max_read_dec_bin < dec_bin) {
              max_read_dec_bin = dec_bin;
            }
          }
          decBinCounts[dec_bin]++;
          decBinCountsInitial[dec_bin]++;
          skip_bin_size = dec_bin - old_dec_bin;
          if ((skip_bin_size > MAX_SKIP_BIN) ||
              (skip_bin_size < -MAX_SKIP_BIN)) {
            printf("Skip of %5d dec bins %5d (%6.2f) -> %5d (%6.2f) at line %10d\n",skip_bin_size,old_dec_bin,oldDec,dec_bin,pInput->Dec,numLines);
          }

          skip_bin_index = skip_bin_size + pGscBin->dec_bins;
          old_dec_bin = dec_bin;
          oldDec = pInput->Dec;
          if (skip_bin_index < 0) {
            skip_bin_index = 0;
          }
          if (skip_bin_index > (2*pGscBin->dec_bins)+1) {
            skip_bin_index > (2*pGscBin->dec_bins)+1;
          }
          skipBinCounts[skip_bin_index]++;
          skipBinCountsTotal[skip_bin_index]++;


          if ((input_nrecs+numFileLines) == 0) {
            max_objid = pInput->objid;
            min_objid = pInput->objid;
          } else {
            if (pInput->objid > max_objid) {
              max_objid = pInput->objid;
            }
            if (pInput->objid < min_objid) {
              min_objid = pInput->objid;
            }

          }

          if (pInput->pmdec > maxpmdec) {
            maxpmdec = pInput->pmdec;
          } 
          if (pInput->pmdec < minpmdec) {
            minpmdec = pInput->pmdec;
          }


  
          temp_objid = pInput->objid;
          modulus = 1;
          for (modindex = 0; modindex < MAX_MODULUS; modindex++) {
            objidbits[modindex] |= temp_objid/modulus;
            modobjidbits[modindex] |= temp_objid % modulus;
            modulus = modulus * MODULUS_BASE;
          }


        

          rawSkipFlag = 0;
          if (raw_handle != NULL) {
            if ((pInput->dg > 0) && 
                (pInput->dr > 0) &&
                (pInput->dg < rmsCutoff) &&
                (pInput->dr < rmsCutoff)) {
              newrms = sqr(pInput->dg) + sqr(pInput->dr);
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
            fprintf(raw_handle,"%f\t%f\t%.2f\t%2f\n",pInput->RA,pInput->Dec,pInput->g,newrms);
          
          }

          tempmag = pInput->g-0.001;
          if (tempmag >= MAGBINS) {
            tempmag = MAGBINS+1;
          }
          if (tempmag < 0) {
            tempmag = MAGBINS+2;
          }

          for (index1 = 0; index1 < ATLAS_MAXMAG; index1++) {
            skipFlag = 0;
            skipFlagBase = 0;
            switch(index1) {
            case ATLAS_BP:
              value1 = pInput->BP;
              rms1   = pInput->dBP;
              break;
            case ATLAS_RP:
              value1 = pInput->RP;
              rms1   = pInput->dRP;
              break;
            case ATLAS_G:
              value1 = pInput->g;
              rms1   = pInput->dg;
              break;
            case ATLAS_R:
              value1 = pInput->r;
              rms1   = pInput->dr;
              break;
            case ATLAS_I:
              value1 = pInput->i;
              rms1   = pInput->di;
              break;
            case ATLAS_Z:
              value1 = pInput->z;
              rms1   = pInput->dz;
              break;
            default:
              printf("ERROR: line %d\n",__LINE__);
              exit(-1);
            }
            if ((rms1 > rmsCutoff) || (rms1 == 0)) {
              skipFlagBase = 1;
            }
            if (skipFlagBase == 0) {
              magcount[index1][tempmag]++;
              apassCount[index1]++;
            }
            for (index2 = 0; index2 < ATLAS_MAXMAG; index2++) {
              skipFlag = skipFlagBase;
              switch(index2) {
              case ATLAS_BP:
                value2 = pInput->BP;
                rms2   = pInput->dBP;
                break;
              case ATLAS_RP:
                value2 = pInput->RP;
                rms2   = pInput->dRP;
                break;
              case ATLAS_G:
                value2 = pInput->g;
                rms2   = pInput->dg;
                break;
              case ATLAS_R:
                value2 = pInput->r;
                rms2   = pInput->dr;
                break;
              case ATLAS_I:
                value2 = pInput->i;
                rms2   = pInput->di;
                break;
              case ATLAS_Z:
                value2 = pInput->z;
                rms2   = pInput->dz;
                break;
              default:
                printf("ERROR: line %d\n",__LINE__);
                exit(-2);
              }
              if ((rms2 > rmsCutoff) && (rms2 == 0)) {
                skipFlag = 1;
              }
              if ((index1 == ATLAS_G) &&
                  (index2 == ATLAS_R)) { 
                newrms = sqrt(sqr(rms1) + sqr(rms2));

#if 0
                printf("line %d rms1 %f rms2 %f sum of squares %f sqrt %f\n",__LINE__,rms1,rms2,sqr(rms1) + sqr(rms2),newrms);
#endif
                intrms = (newrms*RMSBINS)/ATLAS_MAXRMS;
                if (intrms >= RMSBINS) {
                  intrms = RMSBINS+1;
                }
                if (intrms < 0) {
                  intrms = RMSBINS+2;
                }
                rmscount[intrms]++;
                rmscountTotal[intrms][tempmag]++;


              }
              if ((index1 == ATLAS_G) &&
                  (index2 == ATLAS_R))  {
                if (curOutputCount < maxAllocOutput) {
                } else {
                  if (passCounter != 0) {
                    FlushOutputTable(output_handle,pGscAtlasTable,maxAllocOutput,max_read_dec_bin,&curOutputCount,&outputWrittenCount,&output_dec_bin,decBinCountsOutput,0);
                  } else {
                    printf("ERROR: unable to flush the output table in pass %d curOutputCount %d maxAllocOutput %d\n",passCounter,curOutputCount,maxAllocOutput);
                    exit(-1);
                  }
                }

                pGscAtlasImage = &pGscAtlasTable[curOutputCount];
                curOutputCount++;
                memset(pGscAtlasImage,0,sizeof(GSCATLASIMAGE));
                pGscAtlasImage->atlasNumber = input_nrecs+numLines;
                pGscAtlasImage->REFNumber = pGscAtlasImage->atlasNumber + 9000000000L;
                pGscAtlasImage->objid = pInput->objid;
                pGscAtlasImage->ra = pInput->RA;
                pGscAtlasImage->dec = pInput->Dec;
                pGscAtlasImage->Stdmag =  pInput->g;
                pGscAtlasImage->color = pInput->g - pInput->r;
                pGscAtlasImage->rmag = pInput->r;
                pGscAtlasImage->imag = pInput->i;
                pGscAtlasImage->zmag = pInput->z;
                pGscAtlasImage->gmagerr = pInput->dg;
                pGscAtlasImage->rmagerr = pInput->dr;
                pGscAtlasImage->imagerr = pInput->di;
                pGscAtlasImage->zmagerr = pInput->dz;
                pGscAtlasImage->RaPM = pInput->pmra;
                pGscAtlasImage->DecPM = pInput->pmdec;
                pGscAtlasImage->RaSigmaPM = pInput->dpmra;
                pGscAtlasImage->DecSigmaPM = pInput->dpmdec;
                  
                /* Shift the positions to epoch J2000 (see formatgaiadr2.c) */
                dec2 = pGscAtlasImage->dec + (pGscAtlasImage->DecPM * (GSC_EQUINOX-GAIADR2_EPOCH))/(3600.0*1000.0);
                factor = cos(DEGREES_TO_RAD*(dec2));
                if (factor == 0) {
                  ra2 = pGscAtlasImage->ra;
                } else {
                  ra2  = pGscAtlasImage->ra + ((pGscAtlasImage->RaPM * (GSC_EQUINOX-GAIADR2_EPOCH))/(3600.0*1000.0*factor));
                }
                while (ra2 < 0) {
                  ra2 += 360.0;
                } 
                while (ra2 >= 360.0) {
                  ra2 -= 360.0;
                }
                flipFlag = 0;
                if (dec2 == 90.0) {
                  dec2 -= 0.0001;
                } else if (dec2 > 90.0) {
                  dec2 = 180.0-dec2;
                  flipFlag = 1;
                }
                if (dec2 == -90.0) {
                  dec2 += 0.0001;
                } else if (dec2 < -90) {
                  dec2 = -180-dec2;
                  flipFlag = 1;
                }
                if (flipFlag) {
                  ra2 += 180.0;
                  while (ra2 >= 360.0) {
                    ra2 -= 360.0;
                  }
                }
                dec_bin2 = GetDecBin(pGscBin,dec2,"formatatlas3");
                decBinCountsFinal[dec_bin2]++;
                pGscAtlasImage->dec = dec2;
                pGscAtlasImage->ra = ra2;
                

                if ((pGscAtlasImage->color < MIN_ATLAS_COLOR) ||
                    (pGscAtlasImage->color > MAX_ATLAS_COLOR)) {
                  colorReject++;
                  skipFlag = 1;
                }
                outputCount++;
                /* Covert the rms added in quadrature to an integer */
                newrms = (127.*(sqrt(sqr(rms1) + sqr(rms2))))/scalerms;
#if 0
                printf("line %d rms1 %f rms2 %f sum of squares %f sqrt %f scalerms %f \n",__LINE__,rms1,rms2,sqrt(sqr(rms1) + sqr(rms2)),newrms,scalerms);
#endif
                if (newrms > 127) {
                  newrms = 127;
                
                }
              
                intrms = newrms;

                if ((skipFlag) || ((pInput->dupvar & VARIABLE) != 0) || (sqrt(sqr(rms1) + sqr(rms2)) > rmsCutoff))  {
                  pGscAtlasImage->VFlag = 1;
                  skipCount++;
                }
#if 0
                if ((pGscAtlasImage->ra < 185.5) ||
                    (pGscAtlasImage->ra > 186.0) ||
                    (pGscAtlasImage->dec < 63.5) ||
                    (pGscAtlasImage->dec > 64.5)) {
                  continue;
                }
#endif
								
#if 0
                if ((rawSkipFlag == 1)  && (pGscAtlasImage->VFlag == 0)) {
                  printf("X");
                }

#endif

              }
              if (skipFlag == 0) {
                apassDual[index1][index2]++;
              }
            }

          }

          if (verbose && (((input_nrecs+numFileLines) % 10000000) == 0)) {
            time(&curTime);
            curTime -= startTime;
            printf("Reading record %d at mindec %f maxdec %f minra %f maxra %f %d seconds\n",(input_nrecs+numFileLines),mindec,maxdec,minra,maxra,curTime);
            mindec = pInput->Dec;
            maxdec = pInput->Dec;
            minra = pInput->RA;
            maxra = pInput->RA;
     
          }

        
        }
        
      }
      time(&curTime);
      curTime -= startTime;
      printf("Input lines %d, input records %d maximum line length %d\n",numFileLines,(input_nrecs+numFileLines),maxLineLen);
      printf("Reading record %d at mindec %f maxdec %f minra %f maxra %f %d seconds\n",(input_nrecs+numFileLines),mindec,maxdec,minra,maxra,curTime);
      mindec = pInput->Dec;
      maxdec = pInput->Dec;
      minra = pInput->RA;
      maxra = pInput->RA;
#if 0
      for (dec_bin = 0; dec_bin < pGscBin->dec_bins; dec_bin++) {
        if (decBinCounts[dec_bin] > 0) {
          printf("dec_bin %d Declination %8.3f count %d\n",dec_bin,(pGscBin->bin_size*dec_bin)-90.0,decBinCounts[dec_bin]);
        }
      }
#endif
      memset(decBinCounts,0,pGscBin->dec_bins*sizeof(int));
#if 0
      for (skip_bin_index = 0; skip_bin_index < (2*pGscBin->dec_bins)+2; skip_bin_index++) {
        if (skipBinCounts[skip_bin_index] > 0) {
          printf("skip_bin_index %d count %d\n",skip_bin_index-pGscBin->dec_bins,skipBinCounts[skip_bin_index]);
        }
      }
#endif
      memset(skipBinCounts,0,((2*pGscBin->dec_bins)+2)*sizeof(int));
      fclose(input_handle);
    } /* passCounter loop */
    input_nrecs += numFileLines;
    FlushOutputTable(output_handle,pGscAtlasTable,maxAllocOutput,max_read_dec_bin,&curOutputCount,&outputWrittenCount,&output_dec_bin,decBinCountsOutput,0);
  }
  FlushOutputTable(output_handle,pGscAtlasTable,maxAllocOutput,max_read_dec_bin,&curOutputCount,&outputWrittenCount,&output_dec_bin,decBinCountsOutput,1);

  printf("min_pmdec %f max_pmdec %f mindec %f maxdec %f minra %f maxra %f\n",minpmdec,maxpmdec,mindec,maxdec,minra,maxra);
  
  printf("min_objid %lld, max_objid %lld\n",min_objid,max_objid);
  modulus = 1;
  for (modindex = 0; modindex < MAX_MODULUS; modindex++) {
    objidbits[modindex] |= temp_objid/modulus;
    modobjidbits[modindex] |= temp_objid % modulus;
    printf("modindex %2d modulus %20lld objidbits 0x%016llx modobjidbits 0x%016llx\n",modindex,modulus,objidbits[modindex],modobjidbits[modindex]);
    modulus = modulus * MODULUS_BASE;
  }



  if (raw_handle != NULL) {
    fclose(raw_handle);
  }
  fclose(output_handle);
 
  time(&curTime);
  curTime -= startTime;

  printf("Measurement:");
  for (index1 = 0; index1 < ATLAS_MAXMAG; index1++) {
    printf("%10s",apassMeasurementText[index1]);
  }
  printf("\n");

  printf("Count:      ");
  for (index1 = 0; index1 < ATLAS_MAXMAG; index1++) {
    printf("%10d",apassCount[index1]);
  }
  printf("\n");

 
  for (index1 = 0; index1 < ATLAS_MAXMAG; index1++) {
    printf("%10s  ",apassMeasurementText[index1]);
    for (index2 = 0; index2 < ATLAS_MAXMAG; index2++) {
      printf("%10d",apassDual[index1][index2]);
    }
    printf("\n");
  }





  for (tempmag = 0; tempmag <= MAGBINS+2; tempmag++) {
    printf("gmag %2d     ",tempmag);
    for (index1 = 0; index1 < ATLAS_MAXMAG; index1++) {
      printf("%10d",magcount[index1][tempmag]);
    }
    printf("\n");
  }    

  for (intrms = 0; intrms <= RMSBINS+2; intrms++) {
    if (rmscount[intrms] > 0) {
      printf("g-r rms  (mag) %.2f, count %6d\n",(ATLAS_MAXRMS*intrms)/(1.0*RMSBINS),rmscount[intrms]);
    }
  }
#ifdef PLOT_RMS_GMAG
  printf("gmag\trms\tcount\n");
  printf("----\t---\t-----\n");
  

  for (intrms = 0; intrms <= (RMSBINS+2); intrms++) {
    for (tempmag = 0; tempmag <= MAGBINS+2; tempmag++) {
      if (rmscountTotal[intrms][tempmag] != 0) {
        printf("%f\t%f\t%d\n",1.0*tempmag,(ATLAS_MAXRMS*intrms)/(1.0*RMSBINS),rmscountTotal[intrms][tempmag]);
      }

    }
  }
#endif /* PLOT_RMS_GMAG */
  for (dec_bin = 0; dec_bin < pGscBin->dec_bins; dec_bin++) {
    if ((decBinCountsInitial[dec_bin] > 0) && (decBinCountsFinal[dec_bin] > 0) && (decBinCountsOutput[dec_bin] > 0)) {
      if (decBinCountsInitial[dec_bin] > maxDecBinCount) {
        maxDecBinCount = decBinCountsInitial[dec_bin];
        maxDecBinValue = dec_bin;
      }
#if 1
      printf("dec_bin %5d Declination %8.3f count initial: %6d final %6d output %6d\n",dec_bin,(pGscBin->bin_size*dec_bin)-90.0,decBinCountsInitial[dec_bin],decBinCountsFinal[dec_bin],decBinCountsOutput[dec_bin]);
#endif
    }
  }
  printf("maxDecBinValue %d maxDecBinCount %d total\n",maxDecBinValue,maxDecBinCount);
#if 0
  for (skip_bin_index = 0; skip_bin_index < (2*pGscBin->dec_bins)+2; skip_bin_index++) {
    if (skipBinCountsTotal[skip_bin_index] > 0) {
      printf("skip_bin_index %d count %d total\n",skip_bin_index-pGscBin->dec_bins,skipBinCountsTotal[skip_bin_index]);
    }
  }
#endif
  if (pGscAtlasTable != NULL) {
    free(pGscAtlasTable);
  }

  time(&curTime);
  curTime -= startTime;

  printf("zeroBmagCount %d zeroVmagCount %d negZeroCount %d posZeroCount %d\n",zeroBmagCount,zeroVmagCount,negZeroCount,posZeroCount);
  printf("Execution Time: %d seconds; rejected %d rejected because of color %d stars written %d\n",curTime,skipCount,colorReject,outputCount);
  printf("outputWrittenCount %d\n",outputWrittenCount);
  printf("rawSkipCount %d\n",rawSkipCount);

  return(EXIT_SUCCESS);

}
