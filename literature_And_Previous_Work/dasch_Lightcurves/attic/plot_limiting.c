// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* plot_limiting.c
 *
 *  Plot the limiting magnitudes for a series.  This program was written for aavso5.pdf graphs.
 *      Original results are in /dasch/backup/2009_11_26/aavso5
 *                          and /dasch/junk/limiting4773/mag14
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   plot_limiting.c pipelineutils.a -ltable -lutil  -lwcs -o plot_limiting -L/usr/lib64/mysql -lmysqlclient
 *
 * Apr 10, 2009 Edward J. Los - Original Version
 * Apr 28, 2009 Edward J. Los - Add limiting magnitude printout
 * Sep 25, 2009 Edward J. Los - Make deepLimitingMagnitude an input parameter
 * Nov 23, 2009 Edward J. Los - Correct for multiple exposures.  Each solution of a plate will be counted as a separate plate.
 * May 25, 2012 Edward J. Los - Add the summary table for solution 0 limiting magnitudes
 * 
 *
 * plot_limiting 
 * 
 *  
 *  Create per-magnitude plots.
    plot_limiting -v -m 13 -o /home/scanner/junk/plot_limiting/mag13.txt
    plot_limiting -v -m 14 -o /home/scanner/junk/plot_limiting/mag14.txt
    plot_limiting -v -m 15 -o /home/scanner/junk/plot_limiting/mag15.txt
    plot_limiting -v -m 16 -o /home/scanner/junk/plot_limiting/mag16.txt
    plot_limiting -v -m 17 -o /home/scanner/junk/plot_limiting/mag17.txt
    plot_limiting -v -m 18 -o /home/scanner/junk/plot_limiting/mag18.txt


 *  plot_limiting -v -c -q kepler -o /home/scanner/junk/plot_limiting/los.tmp
 *  
 *  The root file (los.tmp) plots all bin 1 limiting magnitudes as a function of year.
 *
 *  psxy los.tmp -R1880/1990/9/19 -Y4i -JX6i/5.5i -P -B20:"Year of Scanned Plate":/1.0:"Bin 1 Limiting Magnitude":WSne -Sc0.05 > los.ps
 *
 *   The limiting file (loslimiting.tmp) contains a listing of all spatial bin entries
 *   that exceed the deepLimitingMagntiude (-m) parameter.
 * 
 *  The "nnn" file contains the bin 1 limiting magnitude percentile as a function of the
 *  year of the scanned plate.
 *
    psxy los000.tmp -R1880/1990/9/19 -Y4i -JX6i/5.5i -P -B20:"Year of Scanned Plate":/1.0:"Bin 1 Limiting Magnitude":WSne -W0.1 -K >! los1.ps
    psxy los025.tmp -R  -JX -P -B -K -O -W5.0 >> los1.ps
    psxy los050.tmp -R  -JX -P -B -K -O -W20.0 >> los1.ps
    psxy los075.tmp -R  -JX -P -B -K -O -W5.0 >> los1.ps
    psxy los100.tmp -R  -JX -P -B -K -O -W0.1 >> los1.ps
    pstext plot_limiting.txt  -R  -JX -P -B -G0/0/0 -O  >> los1.ps
    gv los1.ps &
 *
 * The Bin Output "nnn" file contains the limiting magntude percentile as a function
 * of the spatial bin.
 *

    psxy losbin000.tmp -R0/10/8.0/19.0 -Y4i -JX6i/5.5i -P -B1:"spatial_bin":/1.0:"Limiting Magnitude":WSne -W0.1 -K >! los2.ps
    psxy losbin025.tmp -R  -JX -P -B -K -O -W5.0 >> los2.ps
    psxy losbin050.tmp -R  -JX -P -B -K -O -W20.0 >> los2.ps
    psxy losbin075.tmp -R  -JX -P -B -K -O -W5.0 >> los2.ps
    psxy losbin100.tmp -R  -JX -P -B -K -O -W0.1 >> los2.ps
    pstext plot_limiting2.txt  -R  -JX -P -B -G0/0/0 -O  >> los2.ps
    gv los2.ps &

*
*   The *Shallow* plates plot the year and exposures of those plates which fail to reach
*   the selected deepLimitingMagnitude
*
*   The *Deep* plates plot the year and exposures of those plates which reach the 
*   selected deepLimitingMagnitude
*


    minmax losDeep_a.tmp
    minmax losShallow_a.tmp
    psxy losDeep_a.tmp  -R1880/1960/0/130 -Y4i -JX6i/5.5i -P -B10:"A plates - Year":/10.0:"Exposure":WSne -Sc0.1 -K >! los_a.ps
    psxy losShallow_a.tmp -R  -JX -P -B  -O -S-0.2 >> los_a.ps
    gv los_a.ps &

    minmax losDeep_ac.tmp
    minmax losShallow_ac.tmp
    psxy losDeep_ac.tmp  -R1880/1960/0/130 -Y4i -JX6i/5.5i -P -B10:"AC plates - Year":/10.0:"Exposure":WSne -Sc0.1 -K >! los_ac.ps
    psxy losShallow_ac.tmp -R  -JX -P -B  -O -S-0.2 >> los_ac.ps
    gv los_ac.ps &

    minmax losDeep_ax.tmp
    minmax losShallow_ax.tmp
    psxy losDeep_ax.tmp  -R1880/1960/0/130 -Y4i -JX6i/5.5i -P -B10:"AX plates - Year":/10.0:"Exposure":WSne -Sc0.1 -K >! los_ax.ps
    psxy losShallow_ax.tmp -R  -JX -P -B  -O -S-0.2 >> los_ax.ps
    gv los_ax.ps &

    minmax losDeep_ay.tmp
    minmax losShallow_ay.tmp
    psxy losDeep_ay.tmp  -R1880/1960/0/130 -Y4i -JX6i/5.5i -P -B10:"AY plates - Year":/10.0:"Exposure":WSne -Sc0.1 -K >! los_ay.ps
    psxy losShallow_ay.tmp -R  -JX -P -B  -O -S-0.2 >> los_ay.ps
    gv los_ay.ps &

    minmax losDeep_b.tmp
    minmax losShallow_b.tmp
    psxy losDeep_b.tmp  -R1880/1960/0/200 -Y4i -JX6i/5.5i -P -B10:"B plates - Year":/20.0:"Exposure":WSne -Sc0.1 -K >! los_b.ps
    psxy losShallow_b.tmp -R  -JX -P -B -K -O -S-0.2 >> los_b.ps
    pstext plot_limiting3.txt  -R  -JX -P -B -G0/0/0 -O  >> los_b.ps
    gv los_b.ps &



    minmax losDeep_bm.tmp
    minmax losShallow_bm.tmp
    psxy losDeep_bm.tmp  -R1920/1970/0/200 -Y4i -JX6i/5.5i -P -B10:"BM plates - Year":/20.0:"Exposure":WSne -Sc0.1 -K >! los_bm.ps
    psxy losShallow_bm.tmp -R  -JX -P -B -O -S-0.2 >> los_bm.ps
    gv los_bm.ps &

    minmax losDeep_dnb.tmp
    minmax losShallow_dnb.tmp
    psxy losDeep_dnb.tmp  -R1910/1990/0/200 -Y4i -JX6i/5.5i -P -B10:"DNB plates - Year":/10.0:"Exposure":WSne -Sc0.1 -K >! los_dnb.ps
    psxy losShallow_dnb.tmp -R  -JX -P -B  -O -S-0.2 >> los_dnb.ps
    gv los_dnb.ps &

    
    minmax losDeep_ma.tmp
    minmax losShallow_ma.tmp
    psxy losDeep_ma.tmp  -R1900/1980/0/140 -Y4i -JX6i/5.5i -P -B10:"MA plates - Year":/10.0:"Exposure":WSne -Sc0.1 -K >! los_ma.ps
    psxy losShallow_ma.tmp -R  -JX -P -B  -O -S-0.2 >> los_ma.ps
    gv los_ma.ps &

    minmax losDeep_mc.tmp
    minmax losShallow_mc.tmp
    psxy losDeep_mc.tmp  -R1910/1990/0/200 -Y4i -JX6i/5.5i -P -B10:"MC plates - Year":/10.0:"Exposure":WSne -Sc0.1 -K >! los_mc.ps
    psxy losShallow_mc.tmp -R  -JX -P -B  -O -S-0.2 >> los_mc.ps
    gv los_mc.ps &

    minmax losDeep_mf.tmp
    minmax losShallow_mf.tmp
    psxy losDeep_mf.tmp  -R1910/1960/0/70 -Y4i -JX6i/5.5i -P -B10:"MF plates - Year":/20.0:"Exposure":WSne -Sc0.1 -K >! los_mf.ps
    psxy losShallow_mf.tmp -R  -JX -P -B  -O -S-0.2 >> los_mf.ps
    gv los_mf.ps &



    minmax losDeep_rb.tmp
    minmax losShallow_rb.tmp
    psxy losDeep_rb.tmp  -R1920/1970/0/200 -Y4i -JX6i/5.5i -P -B10:"RB plates - Year":/20.0:"Exposure":WSne -Sc0.1 -K >! los_rb.ps
    psxy losShallow_rb.tmp -R  -JX -P -B -O -S-0.2 >> los_rb.ps
    gv los_rb.ps &
    

    minmax losDeep_rh.tmp
    minmax losShallow_rh.tmp
    psxy losDeep_rh.tmp  -R1920/1970/0/200 -Y4i -JX6i/5.5i -P -B10:"RH plates - Year":/20.0:"Exposure":WSne -Sc0.1 -K >! los_rh.ps
    psxy losShallow_rh.tmp -R  -JX -P -B -K -O -S-0.2 >> los_rh.ps
    pstext plot_limiting3.txt  -R  -JX -P -B -G0/0/0 -O  >> los_rh.ps
    gv los_rh.ps &


 */


#include <math.h>
#include "table.h"
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#define MAX_FILENAME 256
#define MAX_BUFFER 256
#define MAX_LIST_STRING 25
#define SAMPLE_POINTS 49
#define ALLOC_INCREMENT 1000
#define MINIMUM_COUNT 20
#define MAX_QUERY_COUNT 200
#define DRAD_ERROR_FACTOR 3.0
#define MAX_YEARS 110
#define YEAR_INCREMENT 2
#define DEEPEST_MAGNITUDE 20

extern GSCBIN gscBin01;
PGSCBIN pGscBin01 = &gscBin01;
int *gridarray[DEEPEST_MAGNITUDE];

typedef struct _plateentry {
  char series[MAX_SERIES_STRING]; /* Series */
  int plateNumber;
  int mosaicNumber;
  int solutionNumber;
  int exposureNumber;
  double julianDate;   /* Julian Date */
  double year;         /* plate year */
  double limiting_mag[MAX_SPATIAL_BINS+1];
  double max_limiting_mag;
  double exposure;
  
  int naxis1;
  int naxis2;
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double cd1_1;
  double cd1_2;
  double cd2_1;
  double cd2_2;
  int rotation;
  int WCSSource;

  char ctype1[MAX_CTYPE_STRING];
  char ctype2[MAX_CTYPE_STRING];


} PLATEENTRY,*PPLATEENTRY;

typedef struct _series {
  char series[MAX_SERIES_STRING];
  double max_limiting_mag;
} SERIES,*PSERIES;

void IncrementGridArray(int cur_magnitude,int binNumber)
{
  int mag_index;
  int *pGrid;
  /* Increment the array for this magnitude and all brighter magnitudes */
  for (mag_index = cur_magnitude; mag_index >= 0; mag_index--) {
    pGrid = gridarray[mag_index];
    if (pGrid != NULL) {
      pGrid[binNumber]++;
    }

  }

}

int CheckBin(char *coverageArray,int cur_magnitude,struct WorldCoor *wcs,int decBinIndex,int raBinIndex)
{
  PBININDEX pBinIndex;
  int binNumber;
  double ramin;
  double ractr;
  double ramax;
  double decmin;
  double decctr;
  double decmax;
  double decpixel;
  double rapixel;
  int offscl;
  pBinIndex = &pGscBin01->pBinMasterIndex[decBinIndex];
  if (raBinIndex >= pBinIndex->numBins) {
    printf("ERROR: CheckBin has bad raBinIndex\n");
    exit(-1);
  }
  binNumber = raBinIndex+pBinIndex->startBin;
  if (coverageArray[binNumber] != 0) {
    /* This bin has already been studied */
    return(coverageArray[binNumber]-1);
  }
  decmin = (decBinIndex * pGscBin01->bin_size) - 90.0;
  decctr = decmin +  (pGscBin01->bin_size/2.0);
  decmax = decmin + pGscBin01->bin_size;
  ramin =  (360.0 * raBinIndex)/(pBinIndex->numBins);
  ractr = ramin + (180.0/(pBinIndex->numBins));
  ramax = ramin + (360.0/(pBinIndex->numBins));

                 
  /* Check bin center and all edges */
  wcs2pix(wcs,ractr,decctr,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;
    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  wcs2pix(wcs,ramin,decmin,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;

    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  wcs2pix(wcs,ramin,decctr,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;
    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  wcs2pix(wcs,ramin,decmax,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;
    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  wcs2pix(wcs,ractr,decmin,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;
    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  wcs2pix(wcs,ractr,decmax,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;
    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  wcs2pix(wcs,ramax,decmin,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;
    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  wcs2pix(wcs,ramax,decctr,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;
    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  wcs2pix(wcs,ramax,decmax,&decpixel,&rapixel,&offscl);
  if (offscl == 0) {
    coverageArray[binNumber] = 2;
    IncrementGridArray(cur_magnitude,binNumber);
    return(1);
  }
  coverageArray[binNumber] = 1;
  return(0);
}

void StoreResult(PPLATEENTRY pPlateEntry,int verbose,int skipFlag,int *mosaicCount,int *solutionCount,int *fittedCount,int *fittedSolutionCount,int *zeroCrossingCount,char *coverageArray)
{
  double cd[4];
  struct WorldCoor *wcs;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  /* The following array gives RA and DEC for the following plate coordinates 
   *          left    center   right
   *  top      0        1        2
   *  center   3        4        5
   *  bottom   6        7        8
   */
  

  double sra[9]; 
  double sdec[9];
  int ssize = 9;
  int sindex;
  
  double cdec;
  double cra;
  double xctr;
  double yctr;
  int RAindex;
  int DECindex;
  int RAcounter;
  int DECcounter;
  double decval;

  double raval;
  double decpixel;
  double rapixel;
  int offscl;
  int activeFlag;
  int southPoleFlag;
  int northPoleFlag;
  int ariesFlag;
  double minRa;
  double maxRa;
  double minDec;
  double maxDec;
  int minDecBin;
  int maxDecBin;
  int minRaBin;
  int maxRaBin;
  int ctrRaBin = 0;
  int decBinIndex;
  int raBinIndex;
  PBININDEX pBinIndex;
  char plateName[10+MAX_SERIES_STRING];
  int cur_magnitude;

  if (pPlateEntry->plateNumber == 0) {
    return;
  }
  cur_magnitude = pPlateEntry->max_limiting_mag;

  sprintf(plateName,"%6s%05d",pPlateEntry->series,pPlateEntry->plateNumber);
  if (verbose) {
    if (pPlateEntry->WCSSource !=  WCSSOURCE_IMWCS) {
      printf("PLATE: %6s%05d_%02d_01\n",pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->mosaicNumber);


    } else {
      if (pPlateEntry->rotation == 0) {
        printf("PLATE: %6s%05d_%02d_01ww\n",pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->mosaicNumber);
      } else {
        printf("PLATE: %6s%05d_%02d_01r%dww\n",pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->mosaicNumber,pPlateEntry->rotation);
      }
    }
  }
  (*solutionCount)++;
  if (pPlateEntry->solutionNumber == 0) {
    (*mosaicCount)++;
  }
  
  if (pPlateEntry->WCSSource == WCSSOURCE_IMWCS) {
    (*fittedSolutionCount)++;
    if (pPlateEntry->solutionNumber == 0) {
      (*fittedCount)++;
    }

          
    if (strstr(pPlateEntry->ctype1,"DEC")) {
      char tmpPtr[MAX_CTYPE_STRING];
      double dtmp;
      int itmp;
      strcpy(tmpPtr,pPlateEntry->ctype2);
      strcpy(pPlateEntry->ctype2,pPlateEntry->ctype1);
      strcpy(pPlateEntry->ctype1,tmpPtr);
      dtmp = pPlateEntry->crval1;
      pPlateEntry->crval1 = pPlateEntry->crval2;
      pPlateEntry->crval2 = dtmp;
            


      cd[0] = pPlateEntry->cd2_1;
      cd[1] = pPlateEntry->cd2_2;
      cd[2] = pPlateEntry->cd1_1;
      cd[3] = pPlateEntry->cd1_2;

    } else {

      cd[0] = pPlateEntry->cd1_1;
      cd[1] = pPlateEntry->cd1_2;
      cd[2] = pPlateEntry->cd2_1;
      cd[3] = pPlateEntry->cd2_2;

    }



    wcs = wcskinit(pPlateEntry->naxis1,
                   pPlateEntry->naxis2,
                   pPlateEntry->ctype1,
                   pPlateEntry->ctype2,
                   pPlateEntry->crpix1,
                   pPlateEntry->crpix2,
                   pPlateEntry->crval1,
                   pPlateEntry->crval2,
                   cd,
                   0,  /* cdelt1 */
                   0,  /* cdelt2 */
                   0,  /* crota */
                   2000, /* equinox */
                   0);   /* epoch */

    if (nowcs(wcs)) {
      wcsfree(wcs);
      printf("ERROR no wcs\n");
      exit(-1);
    }
    xmin = 0.5;
    ymin = 0.5;
    xmax = 0.5 + (1.0* pPlateEntry->naxis1);
    ymax = 0.5 + (1.0* pPlateEntry->naxis2);
    xctr = 0.5 + (0.5 *pPlateEntry->naxis1);
    yctr = 0.5 + (0.5 *pPlateEntry->naxis2);

    pix2wcs(wcs,xmin,ymax,&sra[0],&sdec[0]);
    pix2wcs(wcs,xctr,ymax,&sra[1],&sdec[1]);
    pix2wcs(wcs,xmax,ymax,&sra[2],&sdec[2]);

    pix2wcs(wcs,xmin,yctr,&sra[3],&sdec[3]);
    pix2wcs(wcs,xctr,yctr,&sra[4],&sdec[4]);
    pix2wcs(wcs,xmax,yctr,&sra[5],&sdec[5]);

    pix2wcs(wcs,xmin,ymin,&sra[6],&sdec[6]);
    pix2wcs(wcs,xctr,ymin,&sra[7],&sdec[7]);
    pix2wcs(wcs,xmax,ymin,&sra[8],&sdec[8]);

    cra = sra[4];
    cdec = sdec[4];
    


    if (skipFlag == 0) {

      

      memset(coverageArray,0,pGscBin01->total_gsc_bins);
      southPoleFlag = 0;
      northPoleFlag = 0;
      ariesFlag = 0;

      /* First see if we are at the pole */
      wcs2pix(wcs,0,-90.0,&decpixel,&rapixel,&offscl);
      if (offscl == 0) {
        southPoleFlag = 1;
#if 0
        printf("South Pole Plate %s\n",plateName);
#endif
      }
      wcs2pix(wcs,0,90.0,&decpixel,&rapixel,&offscl);
      if (offscl == 0) {
        northPoleFlag = 1;
#if 0
        printf("North Pole Plate %s\n",plateName);
#endif
      }
      
      minRa = sra[0];
      maxRa = minRa;
      minDec = sdec[0];
      maxDec = minDec;
      for (sindex = 1; sindex < ssize; sindex++) {
        if (minRa > sra[sindex]) {
          minRa = sra[sindex];
        }
        if (maxRa < sra[sindex]) {
          maxRa = sra[sindex];
        }
        if (minDec > sdec[sindex]) {
          minDec = sdec[sindex];
        }
        if (maxDec < sdec[sindex]) {
          maxDec = sdec[sindex];
        }
       
      }
      if (southPoleFlag == 1) {
        minDec = -90.0;
      }
      if (northPoleFlag == 1) {
        maxDec = 90.0;
      }
      if (maxDec >= 90.0) {
        maxDec = 90.0 -  (pGscBin01->bin_size/2.0);
      }
      if (minDec <= -90.0) {
        minDec = -90.0;
      }
      minDecBin = GetDecBin(pGscBin01,minDec,plateName);
      maxDecBin = GetDecBin(pGscBin01,maxDec,plateName);
      ariesFlag = 0;
      for (decBinIndex = minDecBin; decBinIndex <= maxDecBin; decBinIndex++) {
        pBinIndex = &pGscBin01->pBinMasterIndex[decBinIndex];
        if ((southPoleFlag) || (northPoleFlag)) {
          minRaBin = 0;
          maxRaBin = pBinIndex->numBins - 1;
        } else {
          minRaBin =((minRa * pBinIndex->numBins)/360.0);
          maxRaBin = ((maxRa * pBinIndex->numBins)/360.0);
          if (!CheckBin(coverageArray,cur_magnitude,wcs,decBinIndex,(minRaBin+maxRaBin)/2)) {
            ariesFlag = 1;
#if 0
            printf("Aries Plate %s\n",plateName);
            CheckBin(coverageArray,cur_magnitude,wcs,decBinIndex,(minRaBin+maxRaBin)/2);
#endif
            minRaBin = 0;
            maxRaBin = pBinIndex->numBins - 1;
          } 
        }
        for (raBinIndex = minRaBin; raBinIndex <= maxRaBin; raBinIndex++) {
          CheckBin(coverageArray,cur_magnitude,wcs,decBinIndex,raBinIndex);
        }
          
      }
    

    }
        

    wcsfree(wcs);


  }
  return;

}




int GetGridElement(int *gridarray,int binIndex)
{
  if (gridarray == NULL) {
    printf("ERROR: gridarray is not initialized\n");
    exit(-1);
  }
  return(gridarray[binIndex]);
}


int main(int argc,char *argv[])
{
  char *dotPtr;
  char outfile[MAX_BUFFER];
  char prefix[MAX_BUFFER];
  char suffix[MAX_BUFFER];
  char outlimitingfile[MAX_BUFFER];

  char out000file[MAX_BUFFER];
  char out025file[MAX_BUFFER];
  char out050file[MAX_BUFFER];
  char out075file[MAX_BUFFER];
  char out100file[MAX_BUFFER];

  char outbin000file[MAX_BUFFER];
  char outbin025file[MAX_BUFFER];
  char outbin050file[MAX_BUFFER];
  char outbin075file[MAX_BUFFER];
  char outbin100file[MAX_BUFFER];
  
  char outsummaryfile[MAX_BUFFER];


  FILE *outHandle = NULL;
  FILE *outlimitinghandle = NULL;
  FILE *outsummaryhandle = NULL;
 
  FILE *out000Handle = NULL;
  FILE *out025Handle = NULL;
  FILE *out050Handle = NULL;
  FILE *out075Handle = NULL;
  FILE *out100Handle = NULL;

  FILE *outbin000Handle = NULL;
  FILE *outbin025Handle = NULL;
  FILE *outbin050Handle = NULL;
  FILE *outbin075Handle = NULL;
  FILE *outbin100Handle = NULL;

  int errorFlag = 0;
  int verbose = 0;
  char cmdchar;

  char *argstr;
  char *mysqlhost;
  char *username;
  char *password;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;

  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  time_t startTime;
  time_t curTime;
  EXPOSURE exposureTable;
  PEXPOSURE pExposure = &exposureTable;
  PSERIES pSeriesTable = NULL;
  PSERIES pSeries;

  SERIES allSeries;
  PSERIES pAllSeries = &allSeries;

  int series_nrecs = 0;
  int series_index;

  PPHOT_SPATIAL_BIN spatial_bin_table = NULL;
  PPHOT_SPATIAL_BIN pSpatialBinEntry;

  int spatial_bin_nrecs = 0;
  int spatial_bin_index;

  PPLATEENTRY plate_table = NULL;
  PPLATEENTRY tmp_plate_table;
  PPLATEENTRY pPlateEntry;
  int plate_alloc = 0;
  int plate_nrecs = 0;
  int plate_index;
  double minYear = 9999;
  double maxYear = 0;
  int baseYear;
  int topYear;
  int yearIndex;
  int yearCount[MAX_YEARS];
  double *yearVector[MAX_YEARS];
  double *binVector[MAX_SPATIAL_BINS+1];
  int binCount[MAX_SPATIAL_BINS+1];
  int maxYearCount = 0;
  double *vector;
  int spatial_bin;
  double binmedian;
  double binRMS;
  int nvals;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  int deepLimitingMagnitude = DEEPEST_MAGNITUDE;
  MOSAIC mosaicInfo;
  PMOSAIC pMosaic = &mosaicInfo;
  PHOTPLATES basePhotPlates;
  PPHOTPLATES pPhotPlates = &basePhotPlates;
  int gotAnswer;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  int solutionNumberFix = 0;
  double *seriesVector = NULL;
  int seriesVectorCount = 0;
  int oldSeriesVectorCount;
  double limiting_mag_median;
  double limiting_mag_rms;
  int series_plate_count;
  double max_limiting_mag = 0.00;
  double min_limiting_mag = 99.00;
  int plotCoverage = 0;
  int plotGMT = 0;
  int coverage_start_mag = 0;
  int coverage_end_mag = 0;
  int coverage_index;
  char equcoveragefile[DEEPEST_MAGNITUDE][MAX_BUFFER];
  FILE *equcoveragehandle[DEEPEST_MAGNITUDE];
  char galcoveragefile[DEEPEST_MAGNITUDE][MAX_BUFFER];
  FILE *galcoveragehandle[DEEPEST_MAGNITUDE];
  char *coverageArray = NULL;
 

  int binIndex;
  int decBinIndex;
  int raBinIndex;
  PBININDEX pBinIndex; 
  double decval;
  double raval;
  double decout;
  double raout;
  int mosaicCount = 0;
  int solutionCount = 0;
  int fittedCount = 0;
  int fittedSolutionCount = 0;
  int zeroCrossingCount = 0;

  

  time(&startTime);
  InitBinIndex(pGscBin01);
  memset(pAllSeries,0,sizeof(SERIES));
  memset(equcoveragefile,0,sizeof(equcoveragefile));
  memset(equcoveragehandle,0,sizeof(equcoveragehandle));
  memset(gridarray,0,sizeof(gridarray));
  memset(galcoveragefile,0,sizeof(galcoveragefile));
  memset(galcoveragehandle,0,sizeof(galcoveragehandle));
  strcpy(pAllSeries->series,"all");

  pSeriesTable = (PSERIES)calloc(MAX_SERIES,sizeof(SERIES));
  


  catalogString[0] = 0;
  outfile[0] = 0;

  for (yearIndex = 0; yearIndex < MAX_YEARS; yearIndex++) {
    yearCount[yearIndex] = 0;
    yearVector[yearIndex] = NULL;
  }
  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    binCount[spatial_bin] = 0;
    binVector[spatial_bin] = NULL;
  }

 
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

        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(outfile) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;
        case 'm': /* deepLimitingMagnitude */
        case 'M':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&deepLimitingMagnitude);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode deep limiting magnitude %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              fprintf(stderr,"ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              sprintf(catalogString,"%d",catalogNumber);
            }
          }
          break;






        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'c': /* plot coverage tables */
        case 'C':
          plotCoverage = 1;
          break;

        case 'g': /* plot GMT tables */
        case 'G':
          plotGMT = 1;
          break;


        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Validate arguments */
 


  if (outfile[0] != 0) {
    strcpy(prefix,outfile);
    dotPtr = strstr(prefix,".");
    if (dotPtr == NULL) {
      fprintf(stderr,"ERROR: No suffix in %s\n",outfile);
      exit(-1);
    }
    strcpy(suffix,dotPtr);
    *dotPtr = 0;

    sprintf(outlimitingfile,"%slimiting%s",prefix,suffix);
    sprintf(outsummaryfile,"%ssummary%s",prefix,suffix);

    sprintf(out000file,"%s%03d%s",prefix,0,suffix);
    sprintf(out025file,"%s%03d%s",prefix,25,suffix);
    sprintf(out050file,"%s%03d%s",prefix,50,suffix);
    sprintf(out075file,"%s%03d%s",prefix,75,suffix);
    sprintf(out100file,"%s%03d%s",prefix,100,suffix);


    sprintf(outbin000file,"%sbin%03d%s",prefix,0,suffix);
    sprintf(outbin025file,"%sbin%03d%s",prefix,25,suffix);
    sprintf(outbin050file,"%sbin%03d%s",prefix,50,suffix);
    sprintf(outbin075file,"%sbin%03d%s",prefix,75,suffix);
    sprintf(outbin100file,"%sbin%03d%s",prefix,100,suffix);

    


    outHandle = fopen(outfile,"wt");
    if (outHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output file %s\n",outfile);
      }
    }


    outlimitinghandle = fopen(outlimitingfile,"wt");
    if (outlimitinghandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the Limiting output file %s\n",outlimitingfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Limiting file %s\n",outlimitingfile);
      }
    }

    outsummaryhandle = fopen(outsummaryfile,"wt");
    if (outsummaryhandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the Summary output file %s\n",outsummaryfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Summary file %s\n",outsummaryfile);
      }
      fprintf(outsummaryhandle,"Bin statistics solution 0 of plates that produce photometry data.\n");

      fprintf(outsummaryhandle,"series plates spatial_bins ave_bins_per_plate limiting_mag_median limiting_mag_rms  limiting_mag_range\n");
    }


    out000Handle = fopen(out000file,"wt");
    if (out000Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 000 output file %s\n",out000file);
    } else {
      if (verbose) {
        fprintf(stderr,"Output 000 file %s\n",out000file);
      }
    }
    out025Handle = fopen(out025file,"wt");
    if (out025Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 025 output file %s\n",out025file);
    } else {
      if (verbose) {
        fprintf(stderr,"Output 025 file %s\n",out025file);
      }
    }
    out050Handle = fopen(out050file,"wt");
    if (out050Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 050 output file %s\n",out050file);
    } else {
      if (verbose) {
        fprintf(stderr,"Output 050 file %s\n",out050file);
      }
    }
    out075Handle = fopen(out075file,"wt");
    if (out075Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 075 output file %s\n",out075file);
    } else {
      if (verbose) {
        fprintf(stderr,"Output 075 file %s\n",out075file);
      }
    }
    out100Handle = fopen(out100file,"wt");
    if (out100Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 100 output file %s\n",out100file);
    } else {
      if (verbose) {
        fprintf(stderr,"Output 100 file %s\n",out100file);
      }
    }



    outbin000Handle = fopen(outbin000file,"wt");
    if (outbin000Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 000 bin output file %s\n",outbin000file);
    } else {
      if (verbose) {
        fprintf(stderr,"Bin Output 000 file %s\n",outbin000file);
      }
    }
    outbin025Handle = fopen(outbin025file,"wt");
    if (outbin025Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 025 bin output file %s\n",outbin025file);
    } else {
      if (verbose) {
        fprintf(stderr,"Bin Output 025 file %s\n",outbin025file);
      }
    }
    outbin050Handle = fopen(outbin050file,"wt");
    if (outbin050Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 050 bin output file %s\n",outbin050file);
    } else {
      if (verbose) {
        fprintf(stderr,"Bin Output 050 file %s\n",outbin050file);
      }
    }
    outbin075Handle = fopen(outbin075file,"wt");
    if (outbin075Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 075 bin output file %s\n",outbin075file);
    } else {
      if (verbose) {
        fprintf(stderr,"Bin Output 075 file %s\n",outbin075file);
      }
    }
    outbin100Handle = fopen(outbin100file,"wt");
    if (outbin100Handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 100 bin output file %s\n",outbin100file);
    } else {
      if (verbose) {
        fprintf(stderr,"Bin Output 100 file %s\n",outbin100file);
      }
    }




  } else {
    printf("ERROR: No output file specified\n");
    errorFlag = 1;
  }

  if (errorFlag) {
    printf("Usage: plot_limiting options\n");
    printf("  options: -v verbose\n");
    printf("           -o <output file>\n");
    printf("           -q <catalog>\n");
    printf("           -c plot coverage tables\n");
    printf("           -g invert ra for gmt plotting\n");
    printf("           -m <deep limiting magnitude\n");
    return(-1);
  }

  printf("plot_limiting of %s %s, out file %s deepLimitingMagnitude %d\n",
         __DATE__,__TIME__,outfile,deepLimitingMagnitude);
    

  /* Connect to the database */
  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    fprintf(stderr,"DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }

  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    fprintf(stderr,"ERROR: DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    fprintf(stderr,"ERROR: DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    fprintf(stderr,"ERROR: DASCH_PASSWORD is not defined\n");
    return(-1);
  }
                   
  mysql_init(pConnection);


  if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pConnection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
    }
    return(-1);
  }



  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    fprintf(stderr,"DASCH_PHOT_USERNAME is not defined\n");
    return(-1);
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    fprintf(stderr,"DASCH_PHOT_PASSWORD is not defined\n");
    return(-1);
  }
  mysql_init(pPhotConnection);


  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
    }
    return(-1);
  }

  InitSeriesTable(pConnection,pPhotConnection);

  gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
  if (gotAnswer != 1)  {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }
  if (pPhotGlobal->solutionNumber == PHOT_SOLUTIONNUMBER_YES) {
    solutionNumberFix = 1;
  }


  spatial_bin_nrecs = ReadPhotSpatialBin(pPhotConnection,&spatial_bin_table,NULL,0,0,0,catalogString,solutionNumberFix);
  if (verbose) {
    printf("Read %d records from the spatial bin table\n",spatial_bin_nrecs);
  }
  for (spatial_bin_index = 0; spatial_bin_index < spatial_bin_nrecs; spatial_bin_index++) {
    pSpatialBinEntry = &spatial_bin_table[spatial_bin_index];
    if (pSpatialBinEntry->limiting_mag >= deepLimitingMagnitude) {
      fprintf(outlimitinghandle,"Limiting magnitude %f for plate %s %5d bin %d\n",
              pSpatialBinEntry->limiting_mag,
              pSpatialBinEntry->series,
              pSpatialBinEntry->plateNumber,
              pSpatialBinEntry->spatial_bin);
    }
  }

  /* Now go through the table and extract the information that we are interested in */
  for (spatial_bin_index = 0; spatial_bin_index < spatial_bin_nrecs; spatial_bin_index++) {
    pSpatialBinEntry = &spatial_bin_table[spatial_bin_index];

    for (series_index = 0; series_index < series_nrecs; series_index++) {
      pSeries = &pSeriesTable[series_index];
      if (strcmp(pSeries->series,pSpatialBinEntry->series) == 0) {
        break;
      }
    }

    if (series_index >= series_nrecs) {
      pSeries = &pSeriesTable[series_nrecs];
      series_nrecs++;
      if (series_nrecs >= MAX_SERIES) {
        fprintf(stderr,"ERROR: MAX_SERIES exceeded\n");
        exit(-1);
      }
      strcpy(pSeries->series,pSpatialBinEntry->series);
      pSeries->max_limiting_mag = pSpatialBinEntry->limiting_mag;

    }


    for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
      pPlateEntry = &plate_table[plate_index];
      if ((pPlateEntry->plateNumber == pSpatialBinEntry->plateNumber) &&
          (pPlateEntry->solutionNumber == pSpatialBinEntry->solutionNumber) &&
          (strcmp(pPlateEntry->series,pSpatialBinEntry->series) == 0)) {
        pPlateEntry->limiting_mag[pSpatialBinEntry->spatial_bin] = pSpatialBinEntry->limiting_mag;        
        if (pSpatialBinEntry->limiting_mag > pPlateEntry->max_limiting_mag) {
          pPlateEntry->max_limiting_mag = pSpatialBinEntry->limiting_mag;
        }


        if (pSeries->max_limiting_mag < pSpatialBinEntry->limiting_mag) {
          pSeries->max_limiting_mag = pSpatialBinEntry->limiting_mag;
        }


        break;
      }
      
    }
    if (plate_index >= plate_nrecs) {
      /* Need to allocate a new plate entry */
      if (plate_nrecs <= plate_alloc) {
        plate_alloc += 1000;
        tmp_plate_table = realloc(plate_table,plate_alloc * sizeof(PLATEENTRY));
        if (tmp_plate_table == NULL) {
          fprintf(stderr,"ERROR: failed to reallocate the plate table with size %d\n",plate_alloc * sizeof(PLATEENTRY));
          exit(-1);
        }
        plate_table = tmp_plate_table;
        tmp_plate_table = NULL;
      }
      pPlateEntry = &plate_table[plate_nrecs];
      memset(pPlateEntry,0,sizeof(PLATEENTRY));
      strcpy(pPlateEntry->series,pSpatialBinEntry->series);
      pPlateEntry->plateNumber = pSpatialBinEntry->plateNumber;
      pPlateEntry->solutionNumber = pSpatialBinEntry->solutionNumber;
      pPlateEntry->limiting_mag[pSpatialBinEntry->spatial_bin] = pSpatialBinEntry->limiting_mag;  
      pPlateEntry->max_limiting_mag = pSpatialBinEntry->limiting_mag; 
#if 0 
      if ((strcmp(pPlateEntry->series,"b") == 0) &&
          (pPlateEntry->plateNumber == 45200)) {
        printf("At plate %s%05d\n",pPlateEntry->series,pPlateEntry->plateNumber);
      }
#endif
      if (GetPhotPlate(pPhotConnection,
                       pPlateEntry->series,
                       pPlateEntry->plateNumber,
                       pPhotPlates,
                       "",0) != 0) {
        fprintf(stderr,"ERROR: Failed to get photplates record for %s%05d s%d\n",pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->solutionNumber);
      } else {
        pPlateEntry->mosaicNumber = pPhotPlates->mosaicNumber;
        if(GetMosaicInfo(pConnection,
                         pPlateEntry->series,
                         pPlateEntry->plateNumber,
                         pPlateEntry->mosaicNumber,
                         pPlateEntry->solutionNumber,
                         pMosaic) != 1) {
          fprintf(stderr,"ERROR: Failed to get mosaic record for %s%05d s%d\n",pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->solutionNumber);
	
        } else {
          pPlateEntry->exposureNumber = pMosaic->exposureNumber;
          pPlateEntry->naxis1 = pMosaic->naxis1;
          pPlateEntry->naxis2 = pMosaic->naxis2;
          pPlateEntry->crval1 = pMosaic->crval1;
          pPlateEntry->crval2 = pMosaic->crval2;
          pPlateEntry->crpix1 = pMosaic->crpix1;
          pPlateEntry->crpix2 = pMosaic->crpix2;
          pPlateEntry->cd1_1 = pMosaic->cd1_1;
          pPlateEntry->cd1_2 = pMosaic->cd1_2;
          pPlateEntry->cd2_1 = pMosaic->cd2_1;
          pPlateEntry->cd2_2 = pMosaic->cd2_2;

          pPlateEntry->rotation = pMosaic->rotation;
          pPlateEntry->WCSSource = pMosaic->WCSSource;
          
          if (pPlateEntry->exposureNumber >= 0) {
            if(GetExposureInfo(pConnection,
                               pPlateEntry->series,
                               pPlateEntry->plateNumber,
                               pPlateEntry->exposureNumber,
                               pExposure,
                               0)) {
              double ephemerisDate;
              pPlateEntry->julianDate = fd2jd(pExposure->date);
              ephemerisDate = jd2ep(pPlateEntry->julianDate);
              pPlateEntry->year = ephemerisDate;
              if (pPlateEntry->year > maxYear) {
                maxYear = pPlateEntry->year;
              }
              if (pPlateEntry->year < minYear) {
                minYear = pPlateEntry->year;
              }
              if (strlen(pExposure->exposure) > 0) {
                nvals = sscanf(pExposure->exposure,"%lf",&pPlateEntry->exposure);
                if (nvals != 1) {
                  fprintf(stderr,"ERROR: failed to decode exposure %10s for plate %6s%05d s%d\n",pExposure->exposure,pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->solutionNumber);
                  pPlateEntry->exposure = 0;
                } else {
                  if (pPlateEntry->exposure == 0) {
                    fprintf(stderr,"ERROR: zero exposure for plate %6s%05d s%d\n",pExposure->exposure,pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->solutionNumber);

                  }
                }
              } else {
                fprintf(stderr,"ERROR: null exposure for plate %6s%05d s%d\n",pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->solutionNumber);
                pPlateEntry->exposure = 0;
              }
            

        
        
            } else {
              fprintf(stderr,"ERROR: Failed to get exposure record for %s%05d s%d\n",pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->solutionNumber);
              printf("0\n");
            }
          } else {
            fprintf(stderr,"ERROR: ExposureNumber is %d for %s%05d s%d\n",pPlateEntry->exposureNumber,pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->solutionNumber);

          }
        }
      }
      plate_nrecs++;

    }

  }
  seriesVector = (double *)calloc(MAX_SPATIAL_BINS * plate_nrecs,sizeof(double));
  if (seriesVector == NULL) {
    printf("ERROR: Failed to allocate the series vector\n");
    exit(-1);
  }

  /* Find out how many plates we have for each year */
  baseYear = minYear;
  topYear = maxYear+1;
  if (((topYear-baseYear)/YEAR_INCREMENT) > MAX_YEARS) {
    fprintf(stderr,"ERROR: MAX_YEARS is %d and needs to be %d\n",MAX_YEARS,(topYear-baseYear)/YEAR_INCREMENT);
    exit(-1);
  }
  

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    binVector[spatial_bin] = (double *)calloc(plate_nrecs,sizeof(double));
    if (binVector[spatial_bin] == NULL) {
      fprintf(stderr,"ERROR: failed too allocate binVector index %d of size %d\n",spatial_bin,plate_nrecs);
      exit(-1);
    }
  }


  for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
    pPlateEntry = &plate_table[plate_index];
    yearIndex = (pPlateEntry->year-baseYear)/YEAR_INCREMENT;
    yearCount[yearIndex]++;
  }
  for (yearIndex = 0; yearIndex < ((MAX_YEARS+1)/YEAR_INCREMENT); yearIndex++) {
    if (yearCount[yearIndex] > maxYearCount) {
      maxYearCount = yearCount[yearIndex];
    }
  }
  if (verbose) {
    printf("Maxiumum number of plates in a year is %d\n",maxYearCount);
  }
  for (yearIndex = 0; yearIndex < ((MAX_YEARS+1)/YEAR_INCREMENT); yearIndex++) {
    yearVector[yearIndex] = (double *)calloc(maxYearCount,sizeof(double));
    if (yearVector[yearIndex] == NULL) {
      fprintf(stderr,"ERROR: failed too allocate yearVector index %d of size %d\n",yearIndex,maxYearCount);
      exit(-1);
    }
  }
  memset(yearCount,0,sizeof(yearCount));
  for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
    pPlateEntry = &plate_table[plate_index];
    if ((pPlateEntry->limiting_mag[1] != 0) &&
        (pPlateEntry->year >= baseYear)) {
      yearIndex = (pPlateEntry->year-baseYear)/YEAR_INCREMENT;
      vector = yearVector[yearIndex];
      vector[yearCount[yearIndex]] = pPlateEntry->limiting_mag[1];
      yearCount[yearIndex]++;
    }

    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      if (pPlateEntry->limiting_mag[spatial_bin] != 0) {
        vector=binVector[spatial_bin];
        vector[binCount[spatial_bin]] = pPlateEntry->limiting_mag[spatial_bin];
        binCount[spatial_bin]++;
      }
    }


  }
  for (yearIndex = 0; yearIndex < ((MAX_YEARS+1)/YEAR_INCREMENT); yearIndex++) {
    if (yearCount[yearIndex] > 0) {
      vector = yearVector[yearIndex];
      qsort(vector,yearCount[yearIndex],sizeof(double),dcmp);
      fprintf(out000Handle,"%f %f\n",1.0*((yearIndex * YEAR_INCREMENT)+baseYear-1),vector[0]);
      fprintf(out025Handle,"%f %f\n",1.0*((yearIndex * YEAR_INCREMENT)+baseYear-1),vector[ (yearCount[yearIndex]-1)/4]);
      fprintf(out050Handle,"%f %f\n",1.0*((yearIndex * YEAR_INCREMENT)+baseYear-1),vector[ (yearCount[yearIndex]-1)/2]);
      fprintf(out075Handle,"%f %f\n",1.0*((yearIndex * YEAR_INCREMENT)+baseYear-1),vector[((yearCount[yearIndex]-1)/4)*3]);
      fprintf(out100Handle,"%f %f\n",1.0*((yearIndex * YEAR_INCREMENT)+baseYear-1),vector[ (yearCount[yearIndex]-1)]);
      

    }
  }

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (binCount[spatial_bin] > 0) {
      vector = binVector[spatial_bin];
      qsort(vector,binCount[spatial_bin],sizeof(double),dcmp);
      fprintf(outbin000Handle,"%d %f\n",spatial_bin,vector[0]);
      fprintf(outbin025Handle,"%d %f\n",spatial_bin,vector[ (binCount[spatial_bin]-1)/4]);
      fprintf(outbin050Handle,"%d %f\n",spatial_bin,vector[ (binCount[spatial_bin]-1)/2]);
      fprintf(outbin075Handle,"%d %f\n",spatial_bin,vector[((binCount[spatial_bin]-1)/4)*3]);
      fprintf(outbin100Handle,"%d %f\n",spatial_bin,vector[ (binCount[spatial_bin]-1)]);

    }

  }





  for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
    pPlateEntry = &plate_table[plate_index];
    if (pPlateEntry->limiting_mag[1] != 0) {
      fprintf(outHandle,"%f %f\n",pPlateEntry->year,pPlateEntry->limiting_mag[1]);
    }
  }
  /* Now dump the individual series graphs - one for below the limiting magnitude and the other for above */
  
  for (series_index = 0; series_index < series_nrecs; series_index++) {
    pSeries = &pSeriesTable[series_index];
    FILE *seriesHandleDeep;
    FILE *seriesHandleShallow;
    char seriesDeepFile[MAX_BUFFER];
    char seriesShallowFile[MAX_BUFFER];
    int shallowPoints = 0;
    int deepPoints = 0;
    int zeroPoints = 0;

    if (pSeries->max_limiting_mag < deepLimitingMagnitude) {
      /* Not of interest if we are below the maximum deepest limiting magnitude */
      continue;
    }

    sprintf(seriesDeepFile,"%sDeep_%s%s",prefix,pSeries->series,suffix);
    sprintf(seriesShallowFile,"%sShallow_%s%s",prefix,pSeries->series,suffix);

    seriesHandleDeep = fopen(seriesDeepFile,"wt");
    if (seriesHandleDeep == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 000 output file %s\n",seriesDeepFile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Deep file %s\n",seriesDeepFile);
      }
    }
   

    seriesHandleShallow = fopen(seriesShallowFile,"wt");
    if (seriesHandleShallow == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the 000 output file %s\n",seriesShallowFile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Shallow file %s\n",seriesShallowFile);
      }
    }
   

    for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
      pPlateEntry = &plate_table[plate_index];
      if (strcmp(pPlateEntry->series,pSeries->series) == 0) {
        if (pPlateEntry->exposure == 0.0) {
          zeroPoints++;
        } else {
          if (pPlateEntry->max_limiting_mag < deepLimitingMagnitude) {
            fprintf(seriesHandleShallow,"%f %f\n",pPlateEntry->year,pPlateEntry->exposure);
            shallowPoints++;
          } else {

            fprintf(seriesHandleDeep,"%f %f\n",pPlateEntry->year,pPlateEntry->exposure);
            deepPoints++;
          }
        }

      }
    }
    printf("series %6s totalPoints %4d shallowPoints %4d deepPoints %4d zeroPoints %4d\n",pSeries->series,shallowPoints+deepPoints,shallowPoints,deepPoints,zeroPoints);
    fclose(seriesHandleDeep);
    fclose(seriesHandleShallow);

  }
  /* Here we create the summary file, counting plates and bins in each series for solution 0 only and calculate the median and std of the limiting magnitude */
  for (series_index = 0; series_index <= series_nrecs; series_index++) {
    seriesVectorCount = 0;
    series_plate_count = 0;
    max_limiting_mag = 0.00;
    min_limiting_mag = 99.00;

    if (series_index == series_nrecs) {
      pSeries = pAllSeries;
    } else {
      pSeries = &pSeriesTable[series_index];
    }
    for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
      pPlateEntry = &plate_table[plate_index];
      if (pPlateEntry->solutionNumber != 0) {
        continue;
      }
      if ((series_index != series_nrecs) &&
          (strcmp(pSeries->series,pPlateEntry->series) != 0)) {
        continue;
      }
      oldSeriesVectorCount = seriesVectorCount;
      if (pPlateEntry->max_limiting_mag > max_limiting_mag) {
        max_limiting_mag = pPlateEntry->max_limiting_mag;
      }
      if (pPlateEntry->max_limiting_mag < min_limiting_mag) {
        min_limiting_mag =  pPlateEntry->max_limiting_mag;
      }


      for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
        if (pPlateEntry->limiting_mag[spatial_bin] != 0) {
          seriesVector[seriesVectorCount] =  pPlateEntry->limiting_mag[spatial_bin];
          seriesVectorCount++;

          if (seriesVectorCount >= (MAX_SPATIAL_BINS * plate_nrecs)) {
            printf("ERROR: seriesVectorCount limit exceeded\n");
            exit(-1);
          }
        }
      }
      if (oldSeriesVectorCount != seriesVectorCount) {
        series_plate_count++;
      }

    }
    /* Get series statistics */
    if (seriesVectorCount > 0) {
      CalcMedianAndRMS(seriesVectorCount,0,seriesVector,&limiting_mag_median,&limiting_mag_rms,0,3.0,0);
      fprintf(outsummaryhandle,"%5s %6d      %6d          %3.1f               %5.2f             %5.2f           %5.2f to %5.2f\n",pSeries->series,series_plate_count,seriesVectorCount,(1.0*seriesVectorCount)/(1.0*series_plate_count),limiting_mag_median,limiting_mag_rms,min_limiting_mag,max_limiting_mag);
    }
  }


  /* Plot coverage maps as a function of magnitude */
  if (plotCoverage) {
    coverageArray = (char *)calloc(pGscBin01->total_gsc_bins,sizeof(char));
    coverage_start_mag = min_limiting_mag;
    coverage_end_mag = max_limiting_mag;
    if (coverage_end_mag >= DEEPEST_MAGNITUDE) {
      coverage_end_mag = DEEPEST_MAGNITUDE-1;
    }
    for (coverage_index = coverage_start_mag; coverage_index <= coverage_end_mag; coverage_index++) {
      gridarray[coverage_index] = (int *)calloc(pGscBin01->total_gsc_bins,sizeof(int));
      if (gridarray[coverage_index] == NULL) {
        printf("ERROR: failed to allocate a grid array of size %d\n",pGscBin01->total_gsc_bins);
        exit(-1);
      }

      sprintf(equcoveragefile[coverage_index],"%sequcoverage%02d%s",prefix,coverage_index,suffix);
      equcoveragehandle[coverage_index] = fopen(equcoveragefile[coverage_index],"wt");
      if (equcoveragehandle[coverage_index] == NULL) {
        errorFlag = 1;
        fprintf(stderr,"ERROR: Failed to open the equatorial coverage file %s\n",equcoveragefile[coverage_index]);
      } else {
        if (verbose) {
          fprintf(stderr,"Equatorial Coverage file %s\n",equcoveragefile[coverage_index]);
        }
      }
      sprintf(galcoveragefile[coverage_index],"%sgalcoverage%02d%s",prefix,coverage_index,suffix);
      galcoveragehandle[coverage_index] = fopen(galcoveragefile[coverage_index],"wt");
      if (galcoveragehandle[coverage_index] == NULL) {
        errorFlag = 1;
        fprintf(stderr,"ERROR: Failed to open the galactic coverage file %s\n",galcoveragefile[coverage_index]);
      } else {
        if (verbose) {
          fprintf(stderr,"Galactic Coverage file %s\n",galcoveragefile[coverage_index]);
        }
      }
      
    }

    for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
      pPlateEntry = &plate_table[plate_index];
      if (pPlateEntry->solutionNumber != 0) {
        continue;
      }
      StoreResult(pPlateEntry,0,0,&mosaicCount,&solutionCount,&fittedCount,&fittedSolutionCount,&zeroCrossingCount,coverageArray);
    }


    /* Now write all of the files */
    for (binIndex = 0; binIndex < pGscBin01->total_gsc_bins; binIndex++) {
      pBinIndex = GetSubbins(pGscBin01,binIndex,&raBinIndex,&decBinIndex,(char *)"plot_limiting");
      raval = (360.0 * raBinIndex)/(pBinIndex->numBins) + (180.0/(pBinIndex->numBins));
      decval = (decBinIndex * pGscBin01->bin_size) - 90.0 +  (pGscBin01->bin_size/2.0);
      if ((raval < 0) || (raval > 360.0)) {
        printf("raval error\n");
      }
      if ((decval < -90.0) || (decval > 90.0)) {
        printf("decval error\n");
      }
      raout = raval;
      decout = decval;
      for (coverage_index = coverage_start_mag; coverage_index <= coverage_end_mag; coverage_index++) {
        if (plotGMT) {
          fprintf(equcoveragehandle[coverage_index],"%f %f %d\n",360.0 - raout,decout,GetGridElement(gridarray[coverage_index],binIndex));
        } else {
          fprintf(equcoveragehandle[coverage_index],"%f %f %d\n",raout,decout,GetGridElement(gridarray[coverage_index],binIndex));
        }
      }

      wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&raout,&decout,2000.0);
      if (plotGMT) {
        raout = raout+180.0;
        if (raout >= 360.0) {
          raout = raout-360.0;
        }
      }
      for (coverage_index = coverage_start_mag; coverage_index <= coverage_end_mag; coverage_index++) {
        if (plotGMT) {
          fprintf(galcoveragehandle[coverage_index],"%f %f %d\n",360.0 - raout,decout,GetGridElement(gridarray[coverage_index],binIndex));
        } else {
          fprintf(galcoveragehandle[coverage_index],"%f %f %d\n",raout,decout,GetGridElement(gridarray[coverage_index],binIndex));
        }
      }

    }
   
    for (coverage_index = coverage_start_mag; coverage_index <= coverage_end_mag; coverage_index++) {
      
    }
  }



  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds %d plates  %d series   year: %f-%f\n",curTime,plate_nrecs,series_nrecs,minYear,maxYear);

  /* Done with the MySQL connection */



  mysql_close(pPhotConnection);
  mysql_close(pConnection);

  for (yearIndex = 0; yearIndex < MAX_YEARS; yearIndex++) {
    if (yearVector[yearIndex] != NULL) {
      free(yearVector[yearIndex]);
    }
  }

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (binVector[spatial_bin] != NULL) {
      free(binVector[spatial_bin]);
    }
  }


  if (spatial_bin_table != NULL) {
    free(spatial_bin_table);
  }
  if (plate_table != NULL) {
    free(plate_table);
  }
  if (pSeriesTable != NULL) {
    free(pSeriesTable);
  }
  if (seriesVector != NULL) {
    free(seriesVector);
  }
  if (coverageArray != NULL) {
    free(coverageArray);
  }


  if (plotCoverage) {
    for (coverage_index = coverage_start_mag; coverage_index <= coverage_end_mag; coverage_index++) {
      if (equcoveragehandle[coverage_index] != NULL) {
        fclose(equcoveragehandle[coverage_index]);
      }
      if (galcoveragehandle[coverage_index] != NULL) {
        fclose(galcoveragehandle[coverage_index]);
      }
    }
  }
  if (outHandle != NULL) {
    fclose(outHandle);
  }

  if (outsummaryhandle != NULL) {
    fclose(outsummaryhandle);
  }
  if (outlimitinghandle != NULL) {
    fclose(outlimitinghandle);
  }
  if (out000Handle != NULL) {
    fclose(out000Handle);
  }
  if (out025Handle != NULL) {
    fclose(out025Handle);
  }
  if (out050Handle != NULL) {
    fclose(out050Handle);
  }
  if (out075Handle != NULL) {
    fclose(out075Handle);
  }
  if (out100Handle != NULL) {
    fclose(out100Handle);
  }

  if (outbin000Handle != NULL) {
    fclose(outbin000Handle);
  }
  if (outbin025Handle != NULL) {
    fclose(outbin025Handle);
  }
  if (outbin050Handle != NULL) {
    fclose(outbin050Handle);
  }
  if (outbin075Handle != NULL) {
    fclose(outbin075Handle);
  }
  if (outbin100Handle != NULL) {
    fclose(outbin100Handle);
  }

  return(EXIT_SUCCESS);
}

