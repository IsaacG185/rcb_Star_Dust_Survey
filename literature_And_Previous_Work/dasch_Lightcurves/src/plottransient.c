// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* plottransient.c
 *
 *  This module contains PlotTransientCandidates to perform plotting from the ProcessTransientCandidates() routine called by search_none and update_summary2
 *
 *  For search_none and update_summary2:
 *
 *   gcc -fPIC -ggdb -c -O2   -I/n/sw/plplot-5.9.9/include/plplot -I/usr/include/mysql -I/dasch/install/include -I/usr/include/plplot plottransient.c  -D_FILE_OFFSET_BITS=64 -o plottransient.o
 *
 *  For all other users of photometryutils:
 *
 *    gcc -fPIC -ggdb -c -O2  -DPLOTTRANSIENTSTUB -I/n/sw/plplot-5.9.9/include/plplot -I/usr/include/mysql -I/dasch/install/include -I/usr/include/plplot plottransient.c  -D_FILE_OFFSET_BITS=64 -o plottransientstub.o
 *  
 *
 * Aug 25, 2015 Edward J. Los - Split off from photometryutils because of a problem with libraries on the Odyssey websites.
 * Aug 29, 2015 Edward J. Los - Make the transient points red in the environment plot
 *                              Correct value of "near" count of images in the environment plot
 *                              Reverse the RA sign on the environment plot
 * Aug 30, 2015 Edward J. Los - Add QUALITYFLAGS printout
 * Sep 15, 2015 Edward J. Los - Rename "peakExcess" to "peakExtra"
 *                              Rename "near" to "nearCount".  Do not plot Transient Candidates not assocated with the selected candidate that have FILTER_AMASK_NONE_CANDIDATE2 bits set.
 *                              Plot matched catalog objects in the environment plot with the cyan color
 *                               Add peakRA and peakDec with 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window
 * Sep 28, 2015 Edward J. Los - Add EXCLUSION_STUDY to produce Starbase tables of transients in the exclusion zones defined by /dasch/data/scanner/backup/2015_09_24/candidatedistribution.png
 * Oct  2, 2015 Edward J. Los - Add plotting of multiple exposure plates
 * Oct 28, 2015 Edward J. Los - Add peakDradRMS2, rename peakRA -> RA_TC, peakDec -> Dec_TC, peakCount -> Npts_TC, peakExtra -> Npts_ext and peakDradRMS2 -> Sig_TC
 * Oct 29, 2015 Edward J. Los - Add peakDradRMS3 for the drad rms of points in the selected 90 day window relative to peakRA and peakDec
 *                              Make Sig_TC = max(peakDradRMS2,peakDradRMS3);
 *                              Do not generate the timestamp from the current time - have the user pick the string.
 * Dec 24, 2016 Edward J. Los   If an observation is marked with the badcolorflag for a bad colorterm correction, plot it with red color in the "Outburst Detail" and "Full Lightcurve" plots 
 *                              Add peakBadColorCount for the number of images of the selected transient with bad colorterm correction.  Include this measurement in the id table.
 *                              Add peakDradRMS3 for the positional scatter of the images of the selected transient.
 *                              Create /dasch/Pipeline/candidates/README.png which is the key to all of the plotting symbols.
 * Jan 30, 2017 Edward J. Los - Add limiting magnitude support   
 * Feb 17, 2017 Edward J. Los   Add PLOT_TC_MAG_MARGIN and PLOT_TC_MAG_RANGE to ensure that limiting magnitudes to 3.0 dimmer than the brightest point are plotted.
 *                              Correct date range of second plot.
 * Mar  8, 2017 Edward J. Los   Move the reading of the limiting magnitudes into photometryutils.c
 *                              Move PLOT_LIMITING_MAGNITUDES to photometryutils.h
 * Mar 14, 2017 Edward J. Los   Remove peakDradRMS3, add peakLimitingYears, peakLimitingPoints
 * Apr  4, 2017 Edward J. Los   Use plmtex to show the location of the README files.
 * Dec 30, 2017 Edward J. Los   Add LARGE_TC_MAGNITUDE = 1.0 for selection of largest flares
 * Mar  5, 2018 Edward J. Los   Change the defect symbol from a filled star to an open star for better readability.
 *                              Add peakNoDefectMag to plot
 * Apr  3, 2018 Edward J. Los   Add peakLongOutburst
 *                              Correct bug where there are no limiting magnitudes.
 * May  2, 2018 Edward J. Los   Split LONG_TC_DAYS into
 *                                    LONG_TC_FLARE_DAYS 365.0 days the flare is active
 *                                    LONG_TC_PRE_DAYS (5*365.0) days integration before the transient
 *                                    LONG_TC_SKIP_DAYS (2*365.0) days to skip after finding a transient
 * Aug 13, 2018 Edward J. Los   Define enableRematch to optimize location of transients
 */                             

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mysql.h"
#include <sys/types.h>
#include <grp.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <assert.h>

#include <fcntl.h>
#include <limits.h>
#include <assert.h>
#include "pipelineutils.h"
#include "photometryutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "libwcs/wcscat.h"
#include "plplot.h"
#include "plplotP.h"


#define PLOT_ENVIRONMENT_WIDTH (56.25) /* Size of environment in arcsec */
#define MAX_FILENAME 512
#define HEIGHT_SCALE 0.8   
#define MAX_PLOT_TEXT 40 /* 32 when scale = 1.0; 40 when scale = 0.8; 54 when scale = 0.6 */
#define MAX_LINE_LENGTH 62 /* Line length for the search parameters */
/* #define LOS_PLOT_DEBUG 1 */ /* Debug mode: remove date suffix (no longer needed) */
/* #define LOS_PLOT_HERSHEY 1 *//* Plot Hershey symbols instead of the full lightcurve  */
/* #define LOS_PLOT_HERSHEY2 1 */ /* Plot Hershey symbols instead of the environment plot */
/* #define EXCLUSION_STUDY 1 */ /* Study lightcurves of exclusion zones */
/* Exclusion zones come from /dasch/data/scanner/backup/2015_09_24/candidatedistribution.png */
extern EXCLUSION exclusionTable[];
extern int exclusionTableSize;

double fmax(double x, double y);

/* Generate a plot of interesting transient candidates for placement in the $DASCH_CANDIDATES directory */
void PlotTransientCandidates(PFILECOMMON pFileCommon,PTRANCOMMON pTranCommon,PPARAMETERVECTORSTORE pVectorStore,PPHOTSTARIMAGE pMagnitudeTable,PGALAXYCOMMON pGalaxyCommon,PSTARENTRY pCurStarEntry,int refType,char *nearbyObjects,int AFLAGSMASK1,int QUALITYMASK,PMALMQUIST pMalmquist1,PMALMQUIST pMalmquist3) 
{
#ifdef PLOTTRANSIENTSTUB
  printf("ERROR: Illegal call to PlotTransientCandidates\n");
  exit(-1);
#else /* PLOTTRANSIENTSTUB */
  char *candidatesDirectory;
  char outfile[MAX_FILENAME];
  char buffer[MAX_FILENAME];
  PLFLT height;
  char nearbyObjects2[MAX_NEARBY_OBJECTS_STRING+2];
  int stringlength;
  int stringindex;
  PLFLT *x = NULL;
  PLFLT *y = NULL;
  int* plottedFlag = NULL;
  double xval;
  double yval;
  PLINT npoints;
  double dmin;
  double dmax;
  double yearmin;
  double yearmax;
  double limitingyearmin;
  double limitingyearmax;
  double margin;
  double year;
  int temp_index;
  PMALMQUIST pMalmquist4;
	PPHOTSTARIMAGE pSumStarImage;
	PFILESTARIMAGE pFileStarImage;
  double factor;
  double Sig_TC = 99.0;
  int nearCount = 0;
  static int printDebugError = 1;
  static int printExclusionError = 1;
  int redIndex;
  int haveGalaxyTable;
	int plateIndex;
	PPLATELIMITINGREC pPlateLimitingRec;
  double dateyear;
  int vectorSize;

#ifdef LOS_PLOT_DEBUG
  if (printDebugError) {
    printf("ERROR: Processing %s without a date stamp\n",pCurStarEntry->REF);
    printDebugError = 0;
  }

#endif /* LOS_PLOT_DEBUG */
  assert(sizeof(double) >= sizeof(PLFLT));
#ifdef PLOT_LIMITING_MAGNITUDES
  /* Convert the Julian Date to a year, and eliminate all points with known magnitudes */
  for (plateIndex = 0; plateIndex < pGalaxyCommon->plateCount; plateIndex++) {
    pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
    dateyear = jd2ep(pPlateLimitingRec->geoJulianDate);
    if (plateIndex == 0) {
      limitingyearmin = dateyear;
      limitingyearmax = dateyear;
    } else {
      if (dateyear < limitingyearmin) {
        limitingyearmin = dateyear;
      }
      if (dateyear > limitingyearmax) {
        limitingyearmax = dateyear;
      }
    }
    pPlateLimitingRec->unused = 0;  /* Used to filter out plates with known magnitudes */
    for (temp_index = 0; temp_index < pTranCommon->full_ngood; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      pFileStarImage = pSumStarImage->pFileStarImage;
      if ((pPlateLimitingRec->seriesId == pFileStarImage->seriesId) &&
          (pPlateLimitingRec->plateNumber == pFileStarImage->plateNumber) &&
          (pPlateLimitingRec->solutionNumber == pFileStarImage->solutionNumber)) {
        pPlateLimitingRec->unused = 1; /* This plate has a valid observation -- do not plot a limiting magnitude */
      }
    }
  }
#endif /* PLOT_LIMITING_MAGNITUDES */  

  vectorSize = pFileCommon->allMagnitudeCount;
  if (pTranCommon->npoints > vectorSize) {
    vectorSize = pTranCommon->npoints;
  }

#ifdef PLOT_LIMITING_MAGNITUDES
  if (pGalaxyCommon->plateCount > vectorSize) {
    vectorSize = pGalaxyCommon->plateCount;
  }
#endif /* PLOT_LIMITING_MAGNITUDES */  
  AllocFileCommonVectors(pFileCommon,vectorSize);



  x = (PLFLT *)pFileCommon->vector1;
  y = (PLFLT *)pFileCommon->vector2;
  plottedFlag = (int *)pFileCommon->vector3;

  candidatesDirectory = getenv("DASCH_CANDIDATES");
  if (candidatesDirectory == NULL) {
    printf("ERROR: DASCH_CANDIDATES is not defined\n");
    exit(-1);
  }
  strcpy(outfile,candidatesDirectory);
  strcat(outfile,"/");
  strcat(outfile,pCurStarEntry->REF);
#ifndef LOS_PLOT_DEBUG
  if (strlen(pFileCommon->timestr) > 0) {
    strcat(outfile,"_");
    strcat(outfile,pFileCommon->timestr);
  }
#endif /* LOS_PLOT_DEBUG */
#ifdef EXCLUSION_STUDY
  FILE *exclusionHandle;
  int exclusionIndex;
  PEXCLUSION pExclusion;

  if (printExclusionError) {
    printf("ERROR: Processing %s with exclusion zones enabled\n",pCurStarEntry->REF);
    printExclusionError = 0;
  }
  for (exclusionIndex = 0; exclusionIndex < exclusionTableSize; exclusionIndex++) {
    pExclusion = &exclusionTable[exclusionIndex];
    if ((pCurStarEntry->ra > pExclusion->ra1) &&
        (pCurStarEntry->ra < pExclusion->ra2) &&
        (pCurStarEntry->dec > pExclusion->dec1) &&
        (pCurStarEntry->dec < pExclusion->dec2)) {
      break;
    }  
  }
  if (exclusionIndex >= exclusionTableSize) {
    return;
  }

  strcat(outfile,".db");
  exclusionHandle = fopen(outfile,"wt");
  if (exclusionHandle == NULL) {
    printf("ERROR: plottransient failed to open %s\n",outfile);
    exit(-1);
  }
  fprintf(exclusionHandle,"REF\tDate\tmagcal_magdep\tX_IMAGE\tY_IMAGE\tseries\tplateNumber\tsolutionNumber\n");
  fprintf(exclusionHandle,"---\t----\t-------------\t-------\t-------\t------\t-----------\t--------------\n");
  
  /* Plot all of the selected points in the transient */
  npoints = 0;
  for (temp_index = 0; temp_index < pTranCommon->full_ngood; temp_index++) {
    pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
    pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
    pFileStarImage = pSumStarImage->pFileStarImage;
    if ((pMalmquist4->Date >= pMalmquist1->Date) &&
        (pMalmquist4->Date <= pMalmquist3->Date)) {
      year = jd2ep(pMalmquist4->Date);
      fprintf(exclusionHandle,"%s\t%f\t%.2f\t%f\t%f\t%s\t%d\t%d\n",
              pCurStarEntry->REF,
              year,
              pMalmquist4->magcal_magdep,
              pFileStarImage->X_IMAGE,
              pFileStarImage->Y_IMAGE,
              GetSeriesString(pFileStarImage->seriesId,1),
              pFileStarImage->plateNumber,
              pFileStarImage->solutionNumber);
    }
  }

  fclose(exclusionHandle);

#else /* EXCLUSION_STUDY */
  strcat(outfile,".png");

  plsdev("png");
  plsetopt("-o",outfile);
  plsetopt("-geometry","1440x1080");
  plscolbg(255,255,255);	/* Force the background colour to white */  
  plscol0(1, 0,0,0);		/* Force the foreground colour to black */
  plscol0(15,255,0,0);		/* Move red to 15 */
  plssub(2,2); /* Call for 4 panes */
  plinit(); 



  /* The first plot is the peak detail (upper left) */
  npoints = 0;
  dmin = 0.0;
  dmax = 1.0;
  memset(plottedFlag,0,pFileCommon->vectorAllocCount*sizeof(int));
  for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
    pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
    

    if ((pMalmquist4->Date >= pMalmquist1->Date) &&
        (pMalmquist4->Date <= pMalmquist3->Date)) {
#if 0
      printf("npoints %d x %f, y %f\n",npoints,x[npoints],y[npoints]);
#endif
      if (npoints == 0) {
        dmin = pMalmquist4->magcal_magdep;
        dmax = pMalmquist4->magcal_magdep;
      } else {
        if (dmin > pMalmquist4->magcal_magdep) {
          dmin = pMalmquist4->magcal_magdep;
        }
        if (dmax < pMalmquist4->magcal_magdep) {
          dmax = pMalmquist4->magcal_magdep;
        }
      }
      npoints++;
    }
  }
#ifdef PLOT_LIMITING_MAGNITUDES
#if 0
  printf("dmin %f dmax %f dateyear1 %f dateyear3 %f \n",dmin,dmax,pMalmquist1->Date,pMalmquist3->Date);
  printf("ERROR: test code is enabled in plottransient.c line %d\n",__LINE__);
  pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[0];
  pPlateLimitingRec->geoJulianDate = (pMalmquist1->Date+pMalmquist3->Date)/2.0;
  pPlateLimitingRec->limiting_mag_local = (dmax+dmin)/2.0;
#endif
  dmax += PLOT_TC_MAG_MARGIN; /* Add a bottom margin */
  if ((dmax-dmin) < PLOT_TC_MAG_RANGE) {
    dmax = dmin + PLOT_TC_MAG_RANGE;  /* Make sure the bottom margin is at least 3 magnitudes for the limiting magnitude */
  }
#endif /* PLOT_LIMITING_MAGNITUDES */

  plenv(-0.5,pMalmquist3->Date-pMalmquist1->Date+0.5, dmax+0.1,dmin-0.1, 0, 1);
  pllab("Days","Magnitude","Outburst Detail");
#ifdef PLOT_LIMITING_MAGNITUDES
  /* Plot limiting magnitudes first */
  npoints = 0;
  plcol0(7); /* Set color to gray */
  for (plateIndex = 0; plateIndex < pGalaxyCommon->plateCount; plateIndex++) {
    pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
    if (pPlateLimitingRec->unused != 0) {
      continue;
    }
    
    if ((pPlateLimitingRec->geoJulianDate < pMalmquist1->Date ) ||
        (pPlateLimitingRec->geoJulianDate > pMalmquist3->Date)) {
      continue;
    }
    x[npoints] = pPlateLimitingRec->geoJulianDate - pMalmquist1->Date;
    if (pPlateLimitingRec->limiting_mag_local > (dmax+0.05)) {
      y[npoints] = dmax+0.05;
    } else if (pPlateLimitingRec->limiting_mag_local < (dmin-0.05)) {
      y[npoints] = dmin-0.05;
    } else {
      y[npoints] = pPlateLimitingRec->limiting_mag_local;
    }
#if 0
    printf("npoints %d dateyear %f %f mag %f plot mag %f\n",npoints,pPlateLimitingRec->geoJulianDate,x[npoints],pPlateLimitingRec->limiting_mag_local, y[npoints] );
#endif
    npoints++;
  }
  plssym(0.0,0.5); /* Set to half size */
  plpoin(npoints,x,y,31); /* 31 is a downward pointing arrow */
  plssym(0.0,1.0); /* Set back to default */
#endif /* PLOT_LIMITING_MAGNITUDES */


  for (redIndex = 0; redIndex <= 1; redIndex++) {
    if (redIndex == 0) {
      plcol0(1); /* Set color to black */
    } else {
      plcol0(15); /* set color to red */
    }
    /* Plot only good points */
    npoints = 0;
    
    for (temp_index = 0; temp_index < pTranCommon->full_ngood; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      pFileStarImage = pSumStarImage->pFileStarImage;
      if (((redIndex == 0) && (pSumStarImage->badcolorflag != 0)) ||
          ((redIndex != 0) && (pSumStarImage->badcolorflag == 0))) {
        continue;
      }
      if (((pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) &&
          ((pSumStarImage->quality & QUALITY_MULTIPLE) == 0)) {
        if ((pMalmquist4->Date >= pMalmquist1->Date) &&
            (pMalmquist4->Date <= pMalmquist3->Date)) {
    
          plottedFlag[pMalmquist4->imageIndex] = 1;
    

          x[npoints] = pMalmquist4->Date - pMalmquist1->Date;
          y[npoints] = pMalmquist4->magcal_magdep;
#if 0
          printf("npoints %d x %f, y %f\n",npoints,x[npoints],y[npoints]);
#endif
          npoints++;
        }
      }
    }
    plpoin(npoints,x,y,17); /* 17 is a small black filled circle */
    /* Plot defects */
    npoints = 0;
    for (temp_index = 0; temp_index < pTranCommon->full_ngood; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
      if (plottedFlag[pMalmquist4->imageIndex] != 0) {
        continue;
      }
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      if (((redIndex == 0) && (pSumStarImage->badcolorflag != 0)) ||
          ((redIndex != 0) && (pSumStarImage->badcolorflag == 0))) {
        continue;
      }
      pFileStarImage = pSumStarImage->pFileStarImage;
      if ((pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) {
        continue;
      }
      if ((pMalmquist4->Date >= pMalmquist1->Date) &&
          (pMalmquist4->Date <= pMalmquist3->Date)) {    
        plottedFlag[pMalmquist4->imageIndex] = 1;
        x[npoints] = pMalmquist4->Date - pMalmquist1->Date;
        y[npoints] = pMalmquist4->magcal_magdep;
#if 0
        printf("npoints %d x %f, y %f\n",npoints,x[npoints],y[npoints]);
#endif
        npoints++;
      }
    }
    plpoin(npoints,x,y,12); /* 12 is a small open star */
    /* Now plot multiple exposure plates */
    npoints = 0;
    for (temp_index = 0; temp_index < pTranCommon->full_ngood; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
      if (plottedFlag[pMalmquist4->imageIndex] != 0) {
        continue;
      }
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      if (((redIndex == 0) && (pSumStarImage->badcolorflag != 0)) ||
          ((redIndex != 0) && (pSumStarImage->badcolorflag == 0))) {
        continue;
      }
      pFileStarImage = pSumStarImage->pFileStarImage;
      if ((pMalmquist4->Date >= pMalmquist1->Date) &&
          (pMalmquist4->Date <= pMalmquist3->Date)) {    
        plottedFlag[pMalmquist4->imageIndex] = 1;
        x[npoints] = pMalmquist4->Date - pMalmquist1->Date;
        y[npoints] = pMalmquist4->magcal_magdep;
#if 0
        printf("line %d npoints %d x %f, y %f\n",__LINE__,npoints,x[npoints],y[npoints]);
#endif
        npoints++;
      }
    }
    plpoin(npoints,x,y,14); /* 14 is an open cross */



    /* Now plot everything else */
    npoints = 0;
    for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
      if (plottedFlag[pMalmquist4->imageIndex] != 0) {
        continue;
      }
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      if (((redIndex == 0) && (pSumStarImage->badcolorflag != 0)) ||
          ((redIndex != 0) && (pSumStarImage->badcolorflag == 0))) {
        continue;
      }
      if ((pMalmquist4->Date >= pMalmquist1->Date) &&
          (pMalmquist4->Date <= pMalmquist3->Date)) {
        x[npoints] = pMalmquist4->Date - pMalmquist1->Date;
        y[npoints] = pMalmquist4->magcal_magdep;
#if 0
        printf("npoints %d x %f, y %f\n",npoints,x[npoints],y[npoints]);
#endif
        npoints++;
      }
    }
    plpoin(npoints,x,y,21); /* 21 is a small open circle */
  }
  plcol0(1); /* Set color to black */

#ifndef LOS_PLOT_HERSHEY
  /* The second plot is the full lightcurve (upper right) */
  npoints = 0;
  dmin = 0.0;
  dmax = 1.0;
  yearmax = 0.0;
  yearmin = 0.0;
  memset(plottedFlag,0,pFileCommon->vectorAllocCount*sizeof(int));
  for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
    pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
    year = jd2ep(pMalmquist4->Date);
    if (npoints == 0) {
      dmin = pMalmquist4->magcal_magdep;
      dmax = pMalmquist4->magcal_magdep;
      yearmax = year;
      yearmin = yearmax;
    } else {
      if (yearmin > year) {
        yearmin = year;
      }
      if (yearmax < year) {
        yearmax = year;
      }
      if (dmin > pMalmquist4->magcal_magdep) {
        dmin = pMalmquist4->magcal_magdep;
      }
      if (dmax < pMalmquist4->magcal_magdep) {
        dmax = pMalmquist4->magcal_magdep;
      }
       
    }
#if 0
    printf("npoints %d,x %f y %f,dmin %f dmax %f yearmin %f yearmax %f\n",npoints,x[npoints],y[npoints],dmin,dmax,yearmin,yearmax);
#endif
    npoints++;
    
  }
#ifdef PLOT_LIMITING_MAGNITUDES
  if (pGalaxyCommon->plateCount > 0) {
    if (yearmax < limitingyearmax) {
      yearmax = limitingyearmax;
    }
    if (yearmin > limitingyearmin) {
      yearmin = limitingyearmin; 
    }
  }
  dmax += PLOT_TC_MAG_MARGIN; /* Add a bottom margin */
  if ((dmax-dmin) < PLOT_TC_MAG_RANGE) {
    dmax = dmin + PLOT_TC_MAG_RANGE;  /* Make sure the bottom margin is at least 3 magnitudes for the limiting magnitude */
  }

#endif /* PLOT_LIMITING_MAGNITUDES */
  margin = 0.1*(yearmax-yearmin);
  plenv(yearmin-margin,yearmax+margin, dmax+0.1,dmin-0.1, 0, 1);
  pllab("Year","Magnitude","Full Lightcurve");
#ifdef PLOT_LIMITING_MAGNITUDES
  /* Plot limiting magnitudes first */
  npoints = 0;
  plcol0(7); /* Set color to gray */
  for (plateIndex = 0; plateIndex < pGalaxyCommon->plateCount; plateIndex++) {
    pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
    if (pPlateLimitingRec->unused != 0) {
      continue;
    }
    dateyear = jd2ep(pPlateLimitingRec->geoJulianDate);
    x[npoints] = dateyear;
    if (pPlateLimitingRec->limiting_mag_local > (dmax+0.05)) {
      y[npoints] = dmax+0.05;
    } else if (pPlateLimitingRec->limiting_mag_local < (dmin-0.05)) {
      y[npoints] = dmin-0.05;
    } else {
      y[npoints] = pPlateLimitingRec->limiting_mag_local;
    }
    npoints++;
  }
  plssym(0.0,0.5); /* Set to half size */
  plpoin(npoints,x,y,31); /* 31 is a downward pointing arrow */
  plssym(0.0,1.0); /* Set back to default */

#endif /* PLOT_LIMITING_MAGNITUDES */
  /* Plot only good points (no defects or multiple) */
  for (redIndex = 0; redIndex <= 1; redIndex++) {
    if (redIndex == 0) {
      plcol0(1); /* Set color to black */
    } else {
      plcol0(15); /* set color to red */
    }
    npoints = 0;
    for (temp_index = 0; temp_index < pTranCommon->full_ngood; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      if (((redIndex == 0) && (pSumStarImage->badcolorflag != 0)) ||
          ((redIndex != 0) && (pSumStarImage->badcolorflag == 0))) {
        continue;
      }
      pFileStarImage = pSumStarImage->pFileStarImage;
      if (((pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) &&
          ((pSumStarImage->quality & QUALITY_MULTIPLE) == 0)) {
    
        plottedFlag[pMalmquist4->imageIndex] = 1;
        year = jd2ep(pMalmquist4->Date);
        x[npoints] = year;
        y[npoints] = pMalmquist4->magcal_magdep;
        npoints++;
      }
    }
    plpoin(npoints,x,y,17); /* 17 is a small black filled circle */
    /* Now plot defects */
    npoints = 0;
    for (temp_index = 0; temp_index < pTranCommon->full_ngood; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
      if (plottedFlag[pMalmquist4->imageIndex] != 0) {
        continue;
      }
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      if (((redIndex == 0) && (pSumStarImage->badcolorflag != 0)) ||
          ((redIndex != 0) && (pSumStarImage->badcolorflag == 0))) {
        continue;
      }
      pFileStarImage = pSumStarImage->pFileStarImage;
      if ((pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) {
        continue;
      }
      plottedFlag[pMalmquist4->imageIndex] = 1;
      year = jd2ep(pMalmquist4->Date);
      x[npoints] = year;
      y[npoints] = pMalmquist4->magcal_magdep;
      npoints++;
    }
    plpoin(npoints,x,y,12); /* 12 is a small open star */
    /* Now plot multiple exposure plates */
    npoints = 0;
    for (temp_index = 0; temp_index < pTranCommon->full_ngood; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_MALMQUIST][temp_index];
      if (plottedFlag[pMalmquist4->imageIndex] != 0) {
        continue;
      }
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      if (((redIndex == 0) && (pSumStarImage->badcolorflag != 0)) ||
          ((redIndex != 0) && (pSumStarImage->badcolorflag == 0))) {
        continue;
      }
      pFileStarImage = pSumStarImage->pFileStarImage;
      if ((pSumStarImage->quality & QUALITY_MULTIPLE) == 0) {
        continue;
      }
      plottedFlag[pMalmquist4->imageIndex] = 1;
      year = jd2ep(pMalmquist4->Date);
      x[npoints] = year;
      y[npoints] = pMalmquist4->magcal_magdep;
#if 0
      printf("line %d npoints %d x %f, y %f\n",__LINE__,npoints,x[npoints],y[npoints]);
#endif
      npoints++;
    }
    plpoin(npoints,x,y,14); /* 14 is an open cross */
    npoints = 0;
    /* Now plot everything else */
    for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
      pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
      if (plottedFlag[pMalmquist4->imageIndex] != 0) {
        continue;
      }
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      if (((redIndex == 0) && (pSumStarImage->badcolorflag != 0)) ||
          ((redIndex != 0) && (pSumStarImage->badcolorflag == 0))) {
        continue;
      }
      year = jd2ep(pMalmquist4->Date);
      x[npoints] = year;
      y[npoints] = pMalmquist4->magcal_magdep;
#if 0
      printf("line %d npoints %d,x %f y %f,dmin %f dmax %f yearmin %f yearmax %f\n",__LINE__,npoints,x[npoints],y[npoints],dmin,dmax,yearmin,yearmax);
#endif
      npoints++;
    
    }
    plpoin(npoints,x,y,21); /* 21 is a small open circle */
  }
  plcol0(1); /* Set color to black */

  /* Now highlight the transient */
  npoints = 2;
  x[0] = jd2ep(pMalmquist1->Date);
  y[0] = dmin-0.1;
  x[1] = jd2ep(pMalmquist1->Date);
  y[1] = dmax+0.1;
  pllsty(2);
  plline(npoints,x,y);
  npoints = 2;
  x[0] = jd2ep(pMalmquist3->Date);
  y[0] = dmin-0.1;
  x[1] = jd2ep(pMalmquist3->Date);
  y[1] = dmax+0.1;
  plline(npoints,x,y);
#else /* LOS_PLOT_HERSHEY */
  plenv(-1,32.0,-1.0,32.0, 0, 1 );
  pllab("number","number","Hershey Symbols");
  for (temp_index = 0; temp_index < 32; temp_index++) {
    npoints = 1;
    x[0] = temp_index;
    y[0] = temp_index;
    plpoin(npoints,x,y,temp_index);
  }
#endif /* LOS_PLOT_HERSHEY */


#ifndef LOS_PLOT_HERSHEY2
  /* The third plot is for the surrounding environment (lower left) */
  plenv(-PLOT_ENVIRONMENT_WIDTH,+PLOT_ENVIRONMENT_WIDTH,-PLOT_ENVIRONMENT_WIDTH,+PLOT_ENVIRONMENT_WIDTH,0,1);
  pllab("ra arsec","dec arcsec","Environment");
  npoints = 0;
  memset(plottedFlag,0,pFileCommon->vectorAllocCount*sizeof(int));
  for (temp_index = 0; temp_index < pFileCommon->allMagnitudeCount; temp_index++) {
    pSumStarImage = &pFileCommon->pAllMagnitudeTable[temp_index];
    pFileStarImage = pSumStarImage->pFileStarImage;
    if ((pFileStarImage->REFNumber == 0) &&
        ((pFileStarImage->AFLAGS & FILTER_AMASK_NONE_CANDIDATE2) != 0)) {
      /* High drad unmatched object, reject it */
      continue;
    }
    factor = cos(DEGREES_TO_RAD*pFileStarImage->dec);
    xval = 3600.*factor*(pFileStarImage->ra-pCurStarEntry->ra);
    yval =3600.*(pFileStarImage->dec-pCurStarEntry->dec);
#if 0
    if (yval > PLOT_ENVIRONMENT_WIDTH) {
      printf("At yval %f fabs(yval) %f \n",yval,fabs(yval));
    }
#endif
    if ((fabs(xval) <= PLOT_ENVIRONMENT_WIDTH) && (fabs(yval) <= PLOT_ENVIRONMENT_WIDTH)) {
      x[npoints] = -xval;
      y[npoints] = yval;
#if 0
      printf("line %d npoints %d,x %f y %f, AFLAGS %d REFNumber %lld REF %s\n",__LINE__,npoints,x[npoints],y[npoints],pFileStarImage->AFLAGS,pFileStarImage->REFNumber,pSumStarImage->REF);
#endif
      npoints++;
      
    }
  }
  plpoin(npoints,x,y,1); /* 1 is a dot: "." */
  nearCount += npoints;

  /* Now overplot matched catalog environment points in a different color*/
  plcol0(11); /* set color to cyan */
  npoints = 0;
  for (temp_index = 0; temp_index < pFileCommon->allMagnitudeCount; temp_index++) {
    pSumStarImage = &pFileCommon->pAllMagnitudeTable[temp_index];
    pFileStarImage = pSumStarImage->pFileStarImage;
    if (pFileStarImage->REFNumber == 0) {
      continue;
    }
    factor = cos(DEGREES_TO_RAD*pFileStarImage->dec);
    xval = 3600.*factor*(pFileStarImage->ra-pCurStarEntry->ra);
    yval =3600.*(pFileStarImage->dec-pCurStarEntry->dec);
#if 0
    if (yval > PLOT_ENVIRONMENT_WIDTH) {
      printf("At yval %f fabs(yval) %f \n",yval,fabs(yval));
    }
#endif
    if ((fabs(xval) <= PLOT_ENVIRONMENT_WIDTH) && (fabs(yval) <= PLOT_ENVIRONMENT_WIDTH)) {
      x[npoints] = -xval;
      y[npoints] = yval;
#if 0
      printf("npoints %d,x %f y %f, AFLAGS %d REFNumber %lld REF %s\n",npoints,x[npoints],y[npoints],pFileStarImage->AFLAGS,pFileStarImage->REFNumber,pSumStarImage->REF);
#endif
      npoints++;
      
    }
  }
  plpoin(npoints,x,y,1);  /* 1 is a dot: "." */

  /* Now overplot the transient points as "x" values, but do not do this for multiple exposures */
  plcol0(15); /* set color to red */
  npoints = 0;
  for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
    pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
    pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
    pFileStarImage = pSumStarImage->pFileStarImage;
    if ((pSumStarImage->quality & QUALITY_MULTIPLE) != 0) {
      continue;
    }
    factor = cos(DEGREES_TO_RAD*pFileStarImage->dec);
    xval = 3600.*factor*(pFileStarImage->ra-pCurStarEntry->ra);
    yval = 3600.*(pFileStarImage->dec-pCurStarEntry->dec);
    x[npoints] = xval;
    y[npoints] = yval;


    if ((fabs(xval) <= PLOT_ENVIRONMENT_WIDTH) && (fabs(yval) <= PLOT_ENVIRONMENT_WIDTH)) {
#if 0
      printf("npoints %d,x %f y %f, AFLAGS %d REFNumber %lld REF %s\n",npoints,x[npoints],y[npoints],pFileStarImage->AFLAGS,pFileStarImage->REFNumber,pSumStarImage->REF);
#endif
      if ((refType != REF_TYPE_DASCH) ||
          ((pFileStarImage->AFLAGS & FILTER_AMASK_NONE_CANDIDATE2) == 0)) {
        nearCount--; /* This was double-counted */
      }
    }
#if 0
    printf("npoints %d,x %f y %f, REFNumber %lld\n",npoints,x[npoints],y[npoints],pFileStarImage->REFNumber);
#endif
    npoints++;
  }

  plpoin(npoints,x,y,5); /* 5 is an "x" */

  /* Now overplot the selected transient as "o" values for non-defects */
  npoints = 0;
  for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
    pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
    if ((pMalmquist4->Date >= pMalmquist1->Date) &&
        (pMalmquist4->Date <= pMalmquist3->Date)) {
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      pFileStarImage = pSumStarImage->pFileStarImage;
      if (((pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) &&
          ((pSumStarImage->quality & QUALITY_MULTIPLE) == 0)) {
   
        plottedFlag[pMalmquist4->imageIndex] = 1;
     
        pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
        pFileStarImage = pSumStarImage->pFileStarImage;
        factor = cos(DEGREES_TO_RAD*pFileStarImage->dec);
        x[npoints] = 3600.*factor*(pFileStarImage->ra-pCurStarEntry->ra);
        y[npoints] = 3600.*(pFileStarImage->dec-pCurStarEntry->dec);
#if 0
        printf("npoints %d,x %f y %f, REFNumber %lld\n",npoints,x[npoints],y[npoints],pFileStarImage->REFNumber);
#endif
        npoints++;
      }
    }
  }
  plpoin(npoints,x,y,4); /* 4 is an open circle */
  /* Now overplot the selected transient as "box" values for defects */
  npoints = 0;
  for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
    pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
    if ((pMalmquist4->Date >= pMalmquist1->Date) &&
        (pMalmquist4->Date <= pMalmquist3->Date)) {
      if (plottedFlag[pMalmquist4->imageIndex] != 0) {
        continue;
      }
    
    
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      pFileStarImage = pSumStarImage->pFileStarImage;
      if ((pFileStarImage->AFLAGS & (1<<FILTER_AFLAG_DEFECT)) == 0) {
        continue;
      }     
      plottedFlag[pMalmquist4->imageIndex] = 1;
      factor = cos(DEGREES_TO_RAD*pFileStarImage->dec);
      x[npoints] = 3600.*factor*(pFileStarImage->ra-pCurStarEntry->ra);
      y[npoints] = 3600.*(pFileStarImage->dec-pCurStarEntry->dec);
#if 0
      printf("npoints %d,x %f y %f, REFNumber %lld\n",npoints,x[npoints],y[npoints],pFileStarImage->REFNumber);
#endif
      npoints++;
    
    }
  }
  plpoin(npoints,x,y,6); /* 6 is an open box */

  /* Now overplot the selected transient as an "asterisk" values for multiple exposures */
  npoints = 0;
  for (temp_index = 0; temp_index < pTranCommon->npoints; temp_index++) {
    pMalmquist4 = &pVectorStore->mvector[PARAMETERVECTOR_NPOINTS][temp_index];
    if ((pMalmquist4->Date >= pMalmquist1->Date) &&
        (pMalmquist4->Date <= pMalmquist3->Date)) {
      if (plottedFlag[pMalmquist4->imageIndex] != 0) {
        continue;
      }
    
    
      pSumStarImage = &pMagnitudeTable[pMalmquist4->imageIndex];
      pFileStarImage = pSumStarImage->pFileStarImage;
      if ((pSumStarImage->quality & QUALITY_MULTIPLE) == 0) {
        continue;
      }
      plottedFlag[pMalmquist4->imageIndex] = 1;
      factor = cos(DEGREES_TO_RAD*pFileStarImage->dec);
      x[npoints] = 3600.*factor*(pFileStarImage->ra-pCurStarEntry->ra);
      y[npoints] = 3600.*(pFileStarImage->dec-pCurStarEntry->dec);
#if 0
      printf("line %d npoints %d,x %f y %f, REFNumber %lld\n",__LINE__,npoints,x[npoints],y[npoints],pFileStarImage->REFNumber);
#endif
      npoints++;
    
    }
  }
  plpoin(npoints,x,y,3); /* 3 is an asterisk */
#else /* LOS_PLOT_HERSHEY2 */
  plenv(-1,32.0,-1.0,32.0, 0, 1 );
  pllab("number","number","Hershey Symbols");
  for (temp_index = 0; temp_index < 32; temp_index++) {
    npoints = 1;
    x[0] = temp_index;
    y[0] = temp_index;
    plpoin(npoints,x,y,temp_index);
  }
#endif /* LOS_PLOT_HERSHEY2 */

  plcol0(1); /* Set color to black */
  /* The lower right box contains descriptive information we can fit 32 characters on a line */
  plenv( 0.0, 100.0, 0.0, 100.0, 0, -2);
  plschr(0.0,0.785*HEIGHT_SCALE);
  sprintf(buffer,"See /n/dasch8/Pipeline/Pipeline/candidates/%s/README.png and README2.png",pFileCommon->timestr);
  plmtex("b",2.0,-0.29,0.0,buffer);
  plschr( 0.0, HEIGHT_SCALE);
  height = 100. - (4.0*HEIGHT_SCALE);
  if (pFileCommon->enableRematch == 0) {
    sprintf(buffer,"REF: %s",pCurStarEntry->REF);
  } else {
    sprintf(buffer,"REF: %s TC Optimize",pCurStarEntry->REF);
  }
  plptex(0.0,height,1.0,0.0,0,buffer);
  height -= (8.0*HEIGHT_SCALE);
  sprintf(buffer,"ra: %.4f, dec: %.4f",
          pCurStarEntry->ra,
          pCurStarEntry->dec);
  plptex(0.0,height,1.0,0.0,0,buffer);
  if ((pTranCommon->peakRA < 999.0) && (pTranCommon->peakDec < 99.0)) {
    height -= (8.0*HEIGHT_SCALE);
    sprintf(buffer,"RA_TC: %.4f, Dec_TC: %.4f",
            pTranCommon->peakRA,
            pTranCommon->peakDec);
    plptex(0.0,height,1.0,0.0,0,buffer);
    height -= (8.0*HEIGHT_SCALE);
    if (pTranCommon->peakDradRMS2 < 90.0) {
      if (pTranCommon->peakDradRMS3 < 90.0) {
        Sig_TC = fmax(pTranCommon->peakDradRMS2,pTranCommon->peakDradRMS3);
      } else {
        Sig_TC = pTranCommon->peakDradRMS2;
      }
    } else {
      if (pTranCommon->peakDradRMS3 < 90.0) {
        Sig_TC = pTranCommon->peakDradRMS3;
      }
    }
#if 0
    sprintf(buffer,"Sig_TC %.1f = max(%.1f,%.1f)",
            Sig_TC,
            pTranCommon->peakDradRMS2,
            pTranCommon->peakDradRMS3);
#else
    if (pTranCommon->peakLongOutburst > 0) {
      sprintf(buffer,"peakDradRMS3 %.1f peakLongOutburst %d",
              pTranCommon->peakDradRMS3,pTranCommon->peakLongOutburst);
    } else {
      sprintf(buffer,"peakDradRMS3 %.1f",
              pTranCommon->peakDradRMS3);
    }
#endif
    plptex(0.0,height,1.0,0.0,0,buffer);
  }

  sprintf(buffer,"npoints: %d, ngood: %d nearCount: %d",
          pTranCommon->npoints,
          pTranCommon->full_ngood,
          nearCount);
  height -= (8.0*HEIGHT_SCALE);
  plptex(0.0,height,1.0,0.0,0,buffer);
  sprintf(buffer,"Npts_TC: %d Npts_ext: %d peakBadColorCount: %d",
          pTranCommon->peakCount,
          pTranCommon->peakExtra,
          pTranCommon->peakBadColorCount);
  height -= (8.0*HEIGHT_SCALE);
  plptex(0.0,height,1.0,0.0,0,buffer);
  sprintf(buffer,"peakYear: %.3f; limiting points %2d & years %2d",
          pTranCommon->peakYear,
          pTranCommon->peakLimitingPoints,
          pTranCommon->peakLimitingYears);
  height -= (8.0*HEIGHT_SCALE);
  plptex(0.0,height,1.0,0.0,0,buffer);
  if ((refType != REF_TYPE_DASCH) && (pCurStarEntry->Stdmag < 90.0)) {
    sprintf(buffer,"peakMaxMag: %.1f w/o defect: %.1f Stdmag: %.1f",
            pTranCommon->peakMaxMag,
            pTranCommon->peakNoDefectMag,
            pCurStarEntry->Stdmag);
  } else {
    sprintf(buffer,"peakMaxMag: %.1f w/o defect: %.1f",
            pTranCommon->peakMaxMag,
            pTranCommon->peakNoDefectMag);
  }
  height -= (8.0*HEIGHT_SCALE);
  plptex(0.0,height,1.0,0.0,0,buffer);
  /* #define MAX_PLOT_TEXT 32 */


  stringlength = strlen(nearbyObjects);
  stringindex = 0;
  height -= (8.0*HEIGHT_SCALE);
  while (stringlength > 0) {
    height -= (8.0*HEIGHT_SCALE);
    strcpy(nearbyObjects2,&nearbyObjects[stringindex]);
    nearbyObjects2[MAX_PLOT_TEXT] = 0;
    plptex(0.0,height,1.0,0.0,0,nearbyObjects2);
    stringindex += MAX_PLOT_TEXT;
    stringlength -= MAX_PLOT_TEXT;
  }


#if 0
  height -= (8.0*HEIGHT_SCALE);
  sprintf(buffer,"000000000111111111122222222223333333333444444444455555555556666666666");
  plptex(0.0,height,1.0,0.0,0,buffer);
  height -= (8.0*HEIGHT_SCALE);
  sprintf(buffer,"1234567890123456789012345678901234567890123456789012345678901234567890");
  plptex(0.0,height,1.0,0.0,0,buffer);
  height -= (8.0*HEIGHT_SCALE);
#endif
  height = 16.0-(4.0*HEIGHT_SCALE);
  sprintf(buffer,"AFLAGS: %d QUALITY %d\n",AFLAGSMASK1,(~QUALITYMASK) & DATABASE_QUALITY_MASK);
  plptex(0.0,height,1.0,0.0,0,buffer);
  height -= (8.0*HEIGHT_SCALE);
  sprintf(buffer,"Plot Date: %s\n",pFileCommon->timestr);
  plptex(0.0,height,1.0,0.0,0,buffer);



  plend();
#endif /* EXCLUSION_STUDY */
#endif /* PLOTTRANSIENTSTUB */
} /* End of PlotTransientCandidates( */

void PlotSymbolKey(char *timestr,int enableRematch) 
{
#ifdef PLOTTRANSIENTSTUB
  printf("ERROR: Illegal call to PlotSymbolKey\n");
  exit(-1);
#else /* PLOTTRANSIENTSTUB */
  char *candidatesDirectory;
  char outfile[MAX_FILENAME];
  char buffer[2*MAX_LINE_LENGTH];
  PLFLT height;
  char nearbyObjects2[MAX_NEARBY_OBJECTS_STRING+2];
  int stringlength;
  int stringindex;
  int* plottedFlag = NULL;
  PLFLT x;
  PLFLT y;
  double xval;
  double yval;
  double pointoffset = 2.0;
  double textoffset = 4.0;
#ifdef PLOT_LIMITING_MAGNITUDES
  double lineoffset = 6.05; /* 16 lines */
#else /* PLOT_LIMITING_MAGNITUDES */
  double lineoffset = 6.46; /* 15 lines */
#endif /* PLOT_LIMITING_MAGNITUDES */
  PLINT npoints;
  double dmin;
  double dmax;
  double yearmin;
  double yearmax;
  double margin;
  double year;
  int temp_index;
  PMALMQUIST pMalmquist4;
  PPHOTSTARIMAGE pSumStarImage;
  PFILESTARIMAGE pFileStarImage;
  double factor;
  double Sig_TC = 99.0;
  int nearCount = 0;
  int linelength;
  static int printDebugError = 1;
  static int printExclusionError = 1;

  assert(sizeof(double) >= sizeof(PLFLT));

  candidatesDirectory = getenv("DASCH_CANDIDATES");
  if (candidatesDirectory == NULL) {
    printf("ERROR: DASCH_CANDIDATES is not defined\n");
    exit(-1);
  }
  strcpy(outfile,candidatesDirectory);
  strcat(outfile,"/");
  strcat(outfile,"README");
  strcat(outfile,".png");

  plsdev("png");
  plsetopt("-o",outfile);
  plsetopt("-geometry","1440x1080");
  plscolbg(255,255,255);	/* Force the background colour to white */  
  plscol0(1, 0,0,0);		/* Force the foreground colour to black */
  plscol0(15,255,0,0);		/* Move red to 15 */
  plinit(); 


  plenv(0,100, 0,100, 0, -1);
  if (enableRematch == 0) {
    sprintf(buffer,"Key to Symbols for Run of %s",timestr);
  } else {
    sprintf(buffer,"Key to Symbols for TC Optimized Run of %s",timestr);
  }
  pllab("","",buffer);
  x = 50.0;
  y = 100.0-lineoffset;  
  plptex(x,y,1.0,0.0,0.5,"Points for the Outburst Detail and the Full Lightcurve");
  x = pointoffset;
  y -= lineoffset;
  plpoin(1,&x,&y,17); /* 17 is a small black filled circle */
  plptex(textoffset,y,1.0,0.0,0,"Good point");
  y -= lineoffset;
  plpoin(1,&x,&y,12); /* 12 is a small open star */
  plptex(textoffset,y,1.0,0.0,0,"Defect");
  y -= lineoffset;
  plpoin(1,&x,&y,14); /* 14 is an open cross */
  plptex(textoffset,y,1.0,0.0,0,"Multiple Exposure");
  y -= lineoffset;
  plpoin(1,&x,&y,21); /* 21 is a small open circle */
  plptex(textoffset,y,1.0,0.0,0,"Not a transient point");
#ifdef PLOT_LIMITING_MAGNITUDES
  y -= lineoffset;
  
  plcol0(7); /* set color to grey */
  plssym(0.0,0.5); /* Set to half size */
  plpoin(1,&x,&y,31); /* 31 is an arrow pointing down */
  plssym(0.0,1.0); /* Set back to default */
  plptex(textoffset,y,1.0,0.0,0,"Limiting Magnitude");
#endif /* PLOT_LIMITING_MAGNITUDES */
  y -= lineoffset;
  plcol0(15); /* set color to red */
  plptex(textoffset,y,1.0,0.0,0,"Points with uncertain or missing color calibration are in red.");
  plcol0(1); /* sete color to black */

  y -= lineoffset;
  y -= lineoffset;

  x = 50.0;
  plptex(x,y,1.0,0.0,0.5,"Points for the Environment Plot");
  x = pointoffset;
  y -= lineoffset;
  plpoin(1,&x,&y,1); /* 1 is a dot: "." */
  plptex(textoffset,y,1.0,0.0,0,"Unmatched Object");
  plcol0(11); /* set color to cyan */
  y -= lineoffset;
  plpoin(1,&x,&y,1); /* 1 is a dot: "." */
  plptex(textoffset,y,1.0,0.0,0,"Matched Catalog Object");

  plcol0(15); /* set color to red */

  y -= lineoffset;
  plpoin(1,&x,&y,5); /* 5 is an "x" */
  plptex(textoffset,y,1.0,0.0,0,"Selected lightcurve points outside the transient window.");

  y -= lineoffset;
  plpoin(1,&x,&y,5); /* 5 is an "x" */
  plpoin(1,&x,&y,4); /* 4 is an open circle */
  plptex(textoffset,y,1.0,0.0,0,"Non-multiple and non-defect points in the transient window.");

  y -= lineoffset;
  plpoin(1,&x,&y,5); /* 5 is an "x" */
  plpoin(1,&x,&y,6); /* 6 is an open box */
  plptex(textoffset,y,1.0,0.0,0,"Defect points in the transient window.");

  y -= lineoffset;
  plpoin(1,&x,&y,5); /* 5 is an "x" */
  plpoin(1,&x,&y,3); /* 3 is an "*" */
  plptex(textoffset,y,1.0,0.0,0,"Multiple exposure plate points in the transient window.");




  plend();
  /* Now print the search parameters */


  strcpy(outfile,candidatesDirectory);
  strcat(outfile,"/");
  strcat(outfile,"README2");
  strcat(outfile,".png");

  plsdev("png");
  plsetopt("-o",outfile);
  plsetopt("-geometry","1440x1080");
  plscolbg(255,255,255);	/* Force the background colour to white */  
  plscol0(1, 0,0,0);		/* Force the foreground colour to black */
  plscol0(15,255,0,0);		/* Move red to 15 */
  plinit(); 


  plenv(0,100, 0,100, 0, -1);
  if (enableRematch == 0) {
    sprintf(buffer,"Search Parameters for Run of %s",timestr);
  } else {
    sprintf(buffer,"Search Parameters for TC Optimized Run of %s",timestr);
  }
  pllab("","",buffer);
  lineoffset = 6.05; /* 16 lines */
  x = 50.0;
  y = 100.0-lineoffset;
  textoffset = 1.0;
  x = pointoffset;
#if 0
  sprintf(buffer,"1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890");
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;
#endif
  sprintf(buffer,"Minimum points in the transient: %d\n",
          MIN_TC_POINTS);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"Maximum points in the lightcurve: %d\n",
          MAX_TC_POINTS);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  y -= (lineoffset/2);
  sprintf(buffer,"Minimum days duration of the flare: %.0f\n",
          MIN_TC_DAYS);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"Maximum days duration of the flare: %.0f\n",
          MAX_TC_DAYS);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  y -= (lineoffset/2);
  sprintf(buffer,"Minimum years with limiting magnitudes: %d\n",
          MIN_TC_LIMITING_YEARS);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"Minimum plates with limiting magnitudes: %d\n",
          MIN_TC_LIMITING_POINTS);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  y -= (lineoffset/2);
  sprintf(buffer,"Minimum amplitude of the flare peak: %.2f mag\n",
          MIN_TC_MAGNITUDE);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"Amplitude of largest flares: %.2f mag\n",
          LARGE_TC_MAGNITUDE);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"Avoid stars brighter than %d mag and closer than %d arcsec\n",
          MAX_TC_BRIGHT_MAG,MAX_TC_BRIGHT_ARCSEC);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"Avoid all stars within %d arcsec and \n",
          MAX_TC_NEARBY_ARCSEC);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"within %d mag of flare dimmest or average mag\n",
          MAX_TC_NEARBY_MAG);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset+4.0,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"Any change of > %.1f mag for > %d points over %.f days",
          LONG_TC_MAGNITUDE,LONG_TC_POINTS,LONG_TC_FLARE_DAYS);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;
  sprintf(buffer,"after %.0f days average; then skip %.0f days.",
          LONG_TC_PRE_DAYS,LONG_TC_SKIP_DAYS);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset+4.0,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;

  sprintf(buffer,"Search limits: ra %3.0f to %3.0f deg (%2.0f to %2.0f h) dec %3.0f to %3.0f\n",
          MIN_TC_RIGHTASCENSION,
          MAX_TC_RIGHTASCENSION,
          MIN_TC_RIGHTASCENSION/15.0,
          MAX_TC_RIGHTASCENSION/15.0,
          MIN_TC_DECLINATION,
          MAX_TC_DECLINATION);
  linelength = strlen(buffer);
  if (linelength > MAX_LINE_LENGTH) {
    buffer[MAX_LINE_LENGTH] = 0;
    printf("ERROR: length %d of '%s' is greater than %d\n",linelength,buffer,MAX_LINE_LENGTH);
  }
  plptex(textoffset,y,1.0,0.0,0.0,buffer);
  y -= lineoffset;


  plend();



#endif /* PLOTTRANSIENTSTUB */
} /* End of PlotSymbolKey( */

