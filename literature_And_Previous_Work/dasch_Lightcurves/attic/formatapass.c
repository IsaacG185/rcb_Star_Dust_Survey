// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatapass.c
 *
 * Read the apass input catalog and put it into the gsc232bin.db format
 *        Accept everything:               formatapass -s -v -n -o /dasch/Pipeline/catalogs/apass.db
 *
 *  1.  Use this program to reformat the apass catalog
 *        Accept nobs > 2, rms cutoff 0.1: formatapass -v -c 0.1 -o /dasch/Pipeline/catalogs/apass_temp.db
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
 * gcc -g -O0   -D_FILE_OFFSET_BITS=64 -I/usr/include/mysql formatapass.c pipelineutils.a -L /dasch/install/lib  -L/usr/lib/mysql  -l mysqlclient -lwcs -lm -o formatapass -ldl -pthread 
 *     matchcatalogs of Mar 23 2011 17:19:47, 
 *       matching /dasch/Pipeline/catalogs/ucac3.dat and /dasch/Pipeline/catalogs/apass_temp.dat 
 *       output /dasch/Pipeline/catalogs/matchgsckepler.db total bins 168966386 TOLERANCE 0.000556 degrees
 *      Done  populatedBinCount 6121155 matches 6292508 for 7173325 kepler stars at 699 seconds
 *            median degDrad is 0.000085 deg or 0.3 arcsec max is 0.000556 or 2.0 arcsec.
 *      mv merge.dat apass.dat
 *      mv merge.idx apass.idx
 *      rm apass_temp*
 *
 * Mar 16, 2011 Edward J. Los - Initial version
 * Apr  8, 2011 Edward J. Los - Put the rms in MAGFlag for use in eliminating duplicates.
 * Aug 14, 2011 Edward J. Los - Do not allocate input table unless finding statistics
 * Sep 16, 2011 Edward J. Los - Covert to B and V and correct rms checks
 * Feb 13, 2012 Edward J. Los - Support dr4, use B and V magnitudes.
 *                              Distinguish between "-0.000" and "+0.000"
 * Feb 21, 2012 Edward J. Los - Support dr5, using B and V magnitudes
 * Aug  8, 2013 Edward J. Los - Output REFNumber instead of REF.
 * Apr 15, 2014 Edward J. Los - Add magnitudes and magnitude errors for all colors.
 * Dec 30, 2014 Edward J. Los - Support DR8; add bmag to the raw text output
 * Jan 19, 2015 Edward J. Los - Add rms to the raw text output and print only rms < rmsCutoff
 * Jan 29, 2015 Edward J. Los - Support DR9;
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
#define MAX_BUFFER 256
#define MAGBINS  19   /* Number of magnitude bins */
#define APASS_VMAG   0
#define APASS_BMAG   1
#define APASS_GMAG   2
#define APASS_RMAG   3
#define APASS_IMAG   4
#define APASS_MAXMAG 5



extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  

typedef struct _apassimage {
  int version;
  double ra;
  double raerr;
  double dec;
  double decerr;
  int nobs;
  int mobs;
  double vmag;
  double bminusvmag;
  double bmag;
  double gmag;
  double rmag;
  double imag;
  double vmagerr;
  double bminusvmagerr;
  double bmagerr;
  double gmagerr;
  double rmagerr;
  double imagerr;
  
} APASSIMAGE,*PAPASSIMAGE;

void FindStats(PAPASSIMAGE input_table,double *vector,int input_nrecs,char *label,int offset,int intFlag,int zeroFlag,double cutoff,int acceptNegative,int rmsoffset,double rmsCutoff)
{
  PAPASSIMAGE pInput;
  int index;
  int count;
  int result;
  int ignoreNegative = 0;
  char *intType[2] = {"double","int"};
  double med = 0;
  double rms = 0;
  double rmsValue;
  count = 0;
  
  if ((strcmp(label,"dec") == 0) || (strcmp(label,"bminusvmag") == 0)) {
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
  char raw_name[MAX_INPUT_NAME];
  FILE* output_handle = NULL;
  FILE* raw_handle = NULL;
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
  FILE *input_handle;
  char *apassMeasurementText[APASS_MAXMAG] = {"vmag","bmag","gmag","rmag","imag"};
  double newrms;
  double scalerms = 0;
  double maxrms2 = 0;
  int intrms;
#if 0
#define MAX_INPUT_LINES 100000
#define INPUT_FILES 1
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/master_apass_south_dr2.txt"};

#endif
#if 0
#define MAX_INPUT_LINES 9000000
#define INPUT_FILES 2
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/master_apass_dr1.txt",
    "/dasch/Pipeline/catalogs/master_apass_south_dr2.txt"};
#endif
#if 0
#define MAX_INPUT_LINES 11000000
#define INPUT_FILES 1
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/master_apass_dr2_1.txt"};

#endif
#if 0
#define MAX_INPUT_LINES 19000000
  /* dr3 has 18987754 input lines */
#define INPUT_FILES 1
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/master_apass_dr3.txt"};

#endif
#if 0
#define MAX_INPUT_LINES 100000
  /* test of dr4 */
#define INPUT_FILES 1
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/apassdr4/zm90.sum"
	};

#endif

#if 0
#define MAX_INPUT_LINES 23841000
  /* dr4 23840787 input lines */
#define INPUT_FILES 24
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/apassdr4/zm5.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm10.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm15.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm20.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm25.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm30.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm35.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm40.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm45.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm50.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm55.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm60.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm65.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm70.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm75.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm80.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm85.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zm90.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zp0.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zp5.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zp10.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zp15.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zp20.sum",
    "/dasch/Pipeline/catalogs/apassdr4/zp60.sum"
	};

#endif

#if 0
#define MAX_INPUT_LINES 30121000
  /* dr5 30120683 input lines */
#define INPUT_FILES 36
  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/apassdr5/zm5.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm10.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm15.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm20.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm25.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm30.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm35.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm40.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm45.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm50.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm55.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm60.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm65.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm70.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm75.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm80.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm85.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zm90.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp0.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp5.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp10.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp15.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp20.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp25.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp30.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp35.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp40.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp45.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp50.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp55.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp60.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp65.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp70.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp75.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp80.sum",
    "/dasch/Pipeline/catalogs/apassdr5/zp85.sum"
	};

#endif

#if 0
#define MAX_INPUT_LINES 61200000
  /* dr6_1 61079526 input lines */
#define INPUT_FILES 36
  char *inputnames[INPUT_FILES] = {
		"/dasch/Pipeline/catalogs/apassdr6/zm05_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm10_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm15_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm20_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm25_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm30_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm35_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm40_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm45_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm50_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm55_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm60_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm65_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm70_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm75_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm80_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm85_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zm90_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp00_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp05_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp10_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp15_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp20_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp25_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp30_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp35_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp40_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp45_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp50_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp55_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp60_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp65_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp70_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp75_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp80_6.sum",
		"/dasch/Pipeline/catalogs/apassdr6/zp85_6.sum"

	};

#endif

#if 0
#define MAX_INPUT_LINES  50660000
  /* dr7 50659257 input lines */
#define INPUT_FILES 36
  char *inputnames[INPUT_FILES] = {
		"/dasch/Pipeline/catalogs/apassdr7/zm05_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm10_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm15_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm20_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm25_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm30_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm35_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm40_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm45_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm50_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm55_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm60_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm65_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm70_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm75_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm80_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm85_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zm90_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp00_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp05_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp10_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp15_7.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp20_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp25_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp30_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp35_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp40_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp45_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp50_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp55_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp60_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp65_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp70_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp75_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp80_6.sum",
		"/dasch/Pipeline/catalogs/apassdr7/zp85_6.sum"

	};

#endif

#if 0
#define MAX_INPUT_LINES  57000000 
  /* dr8 55393945  input lines */
#define INPUT_FILES 36
  char *inputnames[INPUT_FILES] = {
		"/dasch/Pipeline/catalogs/apassdr8/zm05_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm10_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm15_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm20_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm25_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm30_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm35_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm40_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm45_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm50_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm55_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm60_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm65_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm70_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm75_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm80_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm85_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zm90_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp00_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp05_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp10_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp15_7.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp20_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp25_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp30_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp35_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp40_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp45_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp50_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp55_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp60_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp65_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp70_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp75_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp80_8.sum",
		"/dasch/Pipeline/catalogs/apassdr8/zp85_8.sum"

	};

#endif

#if 1
#define MAX_INPUT_LINES  62000000 
  /* dr9  61174922 input lines */
#define INPUT_FILES 36
  char *inputnames[INPUT_FILES] = {
		"/dasch/Pipeline/catalogs/apassdr9/zm05_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm10_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm15_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm20_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm25_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm30_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm35_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm40_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm45_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm50_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm55_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm60_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm65_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm70_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm75_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm80_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm85_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zm90_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp00_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp05_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp10_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp15_9.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp20_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp25_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp30_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp35_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp40_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp45_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp50_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp55_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp60_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp65_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp70_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp75_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp80_8.sum",
		"/dasch/Pipeline/catalogs/apassdr9/zp85_8.sum"

	};

#endif



  int inputFile;


  int input_nrecs = 0;
  PAPASSIMAGE pInput;
  PAPASSIMAGE input_table = NULL;
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
  int zeroBmagCount = 0;
  int zeroVmagCount = 0;
  int negZeroCount = 0;
  int posZeroCount = 0;


  memset(apassCount,0,sizeof(apassCount));
  memset(apassDual,0,sizeof(apassDual));
  memset(magcount,0,sizeof(magcount));


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
    printf("Usage: formatapass -o {output name} [-r <raw file>] [-v][-n][-s][-c {cutoff}]\n");
		printf("       where -r <filename> is a raw table of RA and DEC for plotting \n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s computes general statistics\n");
    printf("       where -c {cutoff} gives the RMS cutoff values\n");
    printf("       where -n accepts negative RMS values indicating only one measurement\n");

    return(-1);
  }

  printf("formatapass of %s %s Output Filename %s Raw filename %s acceptNegative %d,rmsCutoff %f\n",
         __DATE__,__TIME__,output_name,raw_name,acceptNegative,rmsCutoff);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of APASSIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(APASSIMAGE));
 

  time(&startTime);
  
  
  fprintf(output_handle,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tgmag\trmag\timag\tvmagerr\tbmagerr\tgmagerr\trmagerr\timagerr\n");
  fprintf(output_handle,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t----\t----\t----\t-------\t-------\t-------\t-------\t-------\n");


  if (doStats) {
    input_table = (PAPASSIMAGE)calloc(MAX_INPUT_LINES,sizeof(APASSIMAGE));
    if (input_table == NULL) {
      printf("ERROR: failed to allocate input table of size %d\n",(MAX_INPUT_LINES*sizeof(APASSIMAGE)));
      exit(-1);
    }
  }

  vector = (double *)calloc(MAX_INPUT_LINES,sizeof(double));
  if (vector == NULL) {
    printf("ERROR: failed to allocate vector of size %d\n",(MAX_INPUT_LINES*sizeof(APASSIMAGE)));
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
      if (doStats) {
        pInput = &input_table[input_nrecs];
      } else {
        pInput = &input_record;
        memset(pInput,0,sizeof(APASSIMAGE));
      }
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
      while ((charPtr = strstr(inBuffer,"-0.000")) != NULL) {
        charPtr += 5;
        negZeroCount++;
        *charPtr = '1';
      }
      while ((charPtr = strstr(inBuffer,"0.000")) != NULL) {
        charPtr += 4;
        posZeroCount++;
        *charPtr = '1';
      }


      nvals = sscanf(inBuffer,"%d %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                     &pInput->version,
                     &pInput->ra,
                     &pInput->raerr,
                     &pInput->dec,
                     &pInput->decerr,
                     &pInput->nobs,
                     &pInput->mobs,
                     &pInput->vmag,
                     &pInput->bminusvmag,
                     &pInput->bmag,
                     &pInput->gmag,
                     &pInput->rmag,
                     &pInput->imag,
                     &pInput->vmagerr,
                     &pInput->bminusvmagerr,
                     &pInput->bmagerr,
                     &pInput->gmagerr,
                     &pInput->rmagerr,
                     &pInput->imagerr);
      if (nvals != 19) {
        printf("ERROR: nvals is %d in line %d\n",nvals,numLines);
      } else {
#if 0
        printf("XXX nobs %d  mobs %d\n",pInput->nobs,pInput->mobs);
#endif
        rawSkipFlag = 0;
        input_nrecs++;
				if (raw_handle != NULL) {
          if ((pInput->vmagerr > 0) && 
              (pInput->bmagerr > 0) &&
              (pInput->vmagerr < rmsCutoff) &&
              (pInput->bmagerr < rmsCutoff)) {
            newrms = sqr(pInput->vmagerr) + sqr(pInput->bmagerr);
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
          fprintf(raw_handle,"%f\t%f\t%.2f\t%2f\n",pInput->ra,pInput->dec,pInput->bmag,newrms);
          
				}
				if (pInput->vmagerr == 0) {
					zeroVmagCount++;
				}
				if (pInput->bmagerr == 0) {
					zeroBmagCount++;
				}

				tempmag = pInput->vmag-0.001;
				if (tempmag >= MAGBINS) {
					tempmag = MAGBINS;
				}

				for (index1 = 0; index1 < APASS_MAXMAG; index1++) {
					skipFlag = 0;
					skipFlagBase = 0;
					switch(index1) {
					case APASS_VMAG:
						value1 = pInput->vmag;
						rms1   = pInput->vmagerr;
						break;
					case APASS_BMAG:
						value1 = pInput->bmag;
						rms1   = pInput->bmagerr;
						break;
					case APASS_GMAG:
						value1 = pInput->gmag;
						rms1   = pInput->gmagerr;
						break;
					case APASS_RMAG:
						value1 = pInput->rmag;
						rms1   = pInput->rmagerr;
						break;
					case APASS_IMAG:
						value1 = pInput->imag;
						rms1   = pInput->imagerr;
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
							value2 = pInput->vmag;
							rms2   = pInput->vmagerr;
							break;
						case APASS_BMAG:
							value2 = pInput->bmag;
							rms2   = pInput->bmagerr;
							break;
						case APASS_GMAG:
							value2 = pInput->gmag;
							rms2   = pInput->gmagerr;
							break;
						case APASS_RMAG:
							value2 = pInput->rmag;
							rms2   = pInput->rmagerr;
							break;
						case APASS_IMAG:
							value2 = pInput->imag;
							rms2   = pInput->imagerr;
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
						if
							((index1 == APASS_BMAG) &&
							 (index2 == APASS_VMAG)) 
								{

									memset(pGscApassImage,0,sizeof(GSCAPASSIMAGE));
									result = GetDASCHNumber(pInput->ra,pInput->dec,&REFNumber,1,REF_TYPE_APASS);
									pGscApassImage->REFNumber = REFNumber;
									pGscApassImage->ra = pInput->ra;
									pGscApassImage->dec = pInput->dec;
									pGscApassImage->Stdmag =  pInput->bmag;
									pGscApassImage->color = pInput->bmag - pInput->vmag;
                  pGscApassImage->gmag = pInput->gmag;
                  pGscApassImage->rmag = pInput->rmag;
                  pGscApassImage->imag = pInput->imag;
                  pGscApassImage->bmagerr = pInput->bmagerr;
                  pGscApassImage->vmagerr = pInput->vmagerr;
                  pGscApassImage->gmagerr = pInput->gmagerr;
                  pGscApassImage->rmagerr = pInput->rmagerr;
                  pGscApassImage->imagerr = pInput->imagerr;
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
								
#if 1
                  if ((rawSkipFlag == 1)  && (pGscApassImage->VFlag == 0)) {
                    printf("X");
                  }

#endif

									fprintf(output_handle,"%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
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
			}
        
		}
		printf("Input lines %d, input records %d\n",numLines,input_nrecs);
		fclose(input_handle);
	}
	if (raw_handle != NULL) {
		fclose(raw_handle);
	}
	if (doStats) {
		printf("     variable  type      count           min             med             max             rms\n");
		FindStats(input_table,vector,input_nrecs,"version"      ,offsetof(APASSIMAGE,version)      ,1,0,-1.0     ,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"ra"           ,offsetof(APASSIMAGE,ra)           ,0,0,-1.0     ,acceptNegative,offsetof(APASSIMAGE,raerr),9.0);
		FindStats(input_table,vector,input_nrecs,"raerr"        ,offsetof(APASSIMAGE,raerr)        ,0,1,-1.0     ,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"dec"          ,offsetof(APASSIMAGE,dec)          ,0,0,-1.0     ,acceptNegative,offsetof(APASSIMAGE,decerr),9.0);
		FindStats(input_table,vector,input_nrecs,"decerr"       ,offsetof(APASSIMAGE,decerr)       ,0,1,-1.0     ,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"nobs"         ,offsetof(APASSIMAGE,nobs)         ,1,0,-1.0     ,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"mobs"         ,offsetof(APASSIMAGE,mobs)         ,1,0,-1.0     ,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"vmag"         ,offsetof(APASSIMAGE,vmag)         ,0,0,90.0     ,acceptNegative,offsetof(APASSIMAGE,vmagerr),rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"bmag"         ,offsetof(APASSIMAGE,bmag)         ,0,0,90.0     ,acceptNegative,offsetof(APASSIMAGE,bmagerr),rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"gmag"         ,offsetof(APASSIMAGE,gmag)         ,0,0,90.0     ,acceptNegative,offsetof(APASSIMAGE,gmagerr),rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"rmag"         ,offsetof(APASSIMAGE,rmag)         ,0,0,90.0     ,acceptNegative,offsetof(APASSIMAGE,rmagerr),rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"imag"         ,offsetof(APASSIMAGE,imag)         ,0,0,90.0     ,acceptNegative,offsetof(APASSIMAGE,imagerr),rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"vmagerr"      ,offsetof(APASSIMAGE,vmagerr)      ,0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"bmagerr"      ,offsetof(APASSIMAGE,bmagerr)      ,0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"gmagerr"      ,offsetof(APASSIMAGE,gmagerr)      ,0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"rmagerr"      ,offsetof(APASSIMAGE,rmagerr)      ,0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"imagerr"      ,offsetof(APASSIMAGE,imagerr)      ,0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);

		FindStats(input_table,vector,input_nrecs,"bminusvmag"   ,offsetof(APASSIMAGE,bminusvmag)   ,0,0,90.0     ,acceptNegative,offsetof(APASSIMAGE,bminusvmagerr),rmsCutoff);
		FindStats(input_table,vector,input_nrecs,"bminusvmagerr",offsetof(APASSIMAGE,bminusvmagerr),0,1,rmsCutoff,acceptNegative,-1,rmsCutoff);
	}
	fclose(output_handle);
 
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
		printf("BMag %2d     ",tempmag);
		for (index1 = 0; index1 < APASS_MAXMAG; index1++) {
			printf("%10d",magcount[index1][tempmag]);
		}
		printf("\n");
	}    


	if (input_table != NULL) {
		free(input_table);
	}

	if (vector != NULL) {
		free(vector);
	}

	time(&curTime);
	curTime -= startTime;

	printf("zeroBmagCount %d zeroVmagCount %d negZeroCount %d posZeroCount %d\n",zeroBmagCount,zeroVmagCount,negZeroCount,posZeroCount);
	printf("Execution Time: %d seconds; rejected %d rejected because of color %d stars written %d\n",curTime,skipCount,colorReject,outputCount);
  printf("rawSkipCount %d\n",rawSkipCount);

	return(EXIT_SUCCESS);
}
