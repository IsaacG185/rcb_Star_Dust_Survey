// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* testmagiso.c
 *
 *  This temporary routine tests MAG_ISO vs Stdmag filtering in matchstars.c and filterblended.c
 * 
 *  gcc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm  testmagiso.c pipelineutils.a -ltable -lutil -lwcs  -L/usr/lib64/mysql -lmysqlclient -o testmagiso 
 * 
 *  Jun 15, 2009 Edward J. Los - Initial version
 *
 */   

#include <math.h>
#include <errno.h>
#include "table.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

/* Sorting routine borrowed from statstable.c by John Roll */
/* zeroFlag is 0 for normal median and RMS 
 *             1 to force the median to zero for a zero-based RMS calculation 
 */
               

int CalcMedianAndRMSBin(int vectorCount,int minGoodStars,int *vectorCount,double *vectorValue,double *med,double *rms,int doClip,double clipFactor,int zeroFlag) 
{
  int index;
  int tmpObsCount;
  int pass;
  double magcal_local;
  double curMedian;
  double curRMS;
  double magcal_diff;
  *med = 0;
  *rms = 0;
  for (pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      /* Clip the existing data set at three sigma */
      tmpObsCount = 0;
      for (index = 0; index < curObsCount; index++) {
        magcal_local = vector[index];
        magcal_diff = curMedian - magcal_local;
        if (magcal_diff < 0.0) {
          magcal_diff = - magcal_diff;
        }
        if (magcal_diff <= (clipFactor*curRMS)) {
          vector[tmpObsCount] = magcal_local;
          tmpObsCount++;
        }
      }
#if 0
      if (curObsCount != tmpObsCount) {
        printf("curObsCount %d tmpObsCount %d\n",curObsCount,tmpObsCount);
      }
#endif
      curObsCount = tmpObsCount;
    }



    if ((curObsCount >= minGoodStars) && (curObsCount > 0)) {
      double aveMagnitude = 0.0;
      double sqrMagnitude = 0.0;
      for (index = 0; index < curObsCount; index++) {
        magcal_local = vector[index];
        aveMagnitude += magcal_local;
        sqrMagnitude += (magcal_local * magcal_local);
      }

      /* Median and Standard deviation borrowed from John Roll's statstable.c */
      double squ;
      double mean;
      double sum;
      double nval;
      int nred;
      nval = curObsCount;
      nred = nval;
      sum = aveMagnitude;
      squ = sqrMagnitude;
      mean = sum/nval;
      qsort(vector,nred,sizeof(double),dcmp);
      if ( nred == 1 ) {
        curMedian = vector[0];
      } else {
        if ( nred % 2 ) {
          curMedian =  vector[(nred/2)-1];
        } else {
          curMedian = (vector[(nred/2)-1] + vector[nred/2])/2;
        }
      }
      if (zeroFlag) {
        curMedian = 0.0;
        mean = 0.0;
      }
      if (nval == 1) {
        curRMS = 99.0;
      } else {
        curRMS = sqrt((squ - 2 * mean * sum + nval * mean * mean)/(nval-1));
      }
    } else {
      return(0);
    }
    if (doClip == 0) {
      break;
    }
  }
  *med = curMedian;
  *rms = curRMS;
  return(curObsCount);
}



#define MAX_BUFFER 512
	
int main(int argc,char *argv[])
{
  int *magnitudeTable = NULL;
  char estimatename[MAX_BUFFER] = "/dasch/Pipeline/match/match_dsb00263_00_01ww_tnx_e.db";
  int iStdmag;
  int iMAG_ISO;
  int binCount;
  double MAG_ISO;
  double Stdmag;
  int spatial_bin = 1;
  char binaryName[MAX_BUFFER];
  char *suffixPtr;
  FILE *binaryHandle;
  size_t inBytes;
  int vectorCount[MAG_ISO_BINS];
  double vectorValue[MAG_ISO_BINS];
  int result;
  double median;
  double rms;

  magnitudeTable = (int *)calloc(MAGNITUDE_BINS*MAX_SPATIAL_BINS,sizeof(int));

  if (magnitudeTable == NULL) {
    fprintf(stderr,"ERROR: Failed to allocate the magnitudesTable\n");
    exit(-1);
  }


  strcpy(binaryName,estimatename);
  suffixPtr = strstr(binaryName,"_e.db");
  if (suffixPtr == NULL) {
    fprintf(stderr,"ERROR: estimatename %s has unknown suffix\n",estimatename);
  } else {
    *suffixPtr = 0;
    strcat(binaryName,"_e.bin");
    binaryHandle = fopen(binaryName,"rb");
    if (binaryHandle == NULL) {
      fprintf(stderr,"ERROR: failed to open %s\n",binaryName);
      exit(-1);
    } else {
      inBytes = fread(magnitudeTable,sizeof(int),MAGNITUDE_BINS*MAX_SPATIAL_BINS,binaryHandle);
      if (inBytes != (MAGNITUDE_BINS*MAX_SPATIAL_BINS)) {
        fprintf(stderr,"ERROR: wrote only %d of %d to %s\n",inBytes,MAGNITUDE_BINS*MAX_SPATIAL_BINS,binaryName);
        exit(-1);
      } else {
        printf("Read %d of %d entries to %s\n",inBytes,MAGNITUDE_BINS*MAX_SPATIAL_BINS,binaryName);
      }
      fclose(binaryHandle);
    }
  }

  for (iStdmag = 0; iStdmag < STDMAG_BINS; iStdmag++) {
    Stdmag = MIN_STDMAG + ((1.0 * iStdmag)*((MAX_STDMAG-MIN_STDMAG)/STDMAG_BINS));
      for (iMAG_ISO = 0; iMAG_ISO < MAG_ISO_BINS; iMAG_ISO++) {
        



        MAG_ISO = MIN_MAG_ISO + ((1.0 * iMAG_ISO) * ((MAX_MAG_ISO-MIN_MAG_ISO)/MAG_ISO_BINS));
        vectorValue[iMAG_ISO] = MAG_ISO;
        binCount = magnitudeTable[iMAG_ISO + (MAG_ISO_BINS * iStdmag)+ ((spatial_bin-1) * MAGNITUDE_BINS)];
        vectorCount[iMAG_ISO] = binCount;
        
      }
      result = CalcMedianAndRMSBin(MAG_ISO_BINS,0,vectorCount,vectorValue,&median,&rms,1,3.0,0);
  }
  return;
}
