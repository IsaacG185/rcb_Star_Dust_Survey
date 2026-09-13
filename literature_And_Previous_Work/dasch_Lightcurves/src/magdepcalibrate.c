// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* magdepcalibrate.c
 *
 *  gcc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I/usr/include/mysql  -I /dasch/install/include  -I /dasch/install/pgplot   magdepcalibrate.c -L /dasch/install/lib -lm pipelineutils.a -L /dasch/install/pgplot  -L/usr/lib${lib64}/mysql  -lmysqlclient  -lcpgplot -lpgplot -L/usr/X11R6/lib  -lX11 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 /usr/lib${lib64}/libg2c.so.0  -ldl -pthread   -ltable -lutil -lwcs -o magdepcalibrate
 *
 *
 * Oct  4, 2011 Edward J. Los -  modelled after Sumin Tang's bright_dubious_variales.m  of Sep 24, 2011 (Matlab)
 *                               ported to Octave as calibratebright.m
 *                               NOTE line references apply to ~/backup/2011_09_22/bright_dubious_variales.m 
 * Oct 27, 2011 Edward J. Los -  Remove debugging statments
 *                               Make plotting permanent
 *                               Increase the line 99 factor to 1.5
 *                               Change the following variables:
 *                                   nstar0 -> nstar_magdep
 *                                   xcoord -> xcoord_magdep
 *                                   ycoord -> ycoord_magdep
 *                                   bmag0  -> magdep_bin_edge
 *                                   M0     -> magdep_bin_median
 *                                   Zrms2  -> magcal_magdep_rms
 *                                   Zfinal -> magdep_bin_magcor
 *                               Do not attempt to  introduce a new magcal_magdep_rms because of problems with crossing bins.
 *                               Introduce magdep_bin_quality = (magdep_bin_magcor * sqrt(nstar_magdep))/magcal_magdep_rms which must
 *                                   be greater than 1 for the final correction to be applied.  When interpolating, assume that adjacent
 *                                   magnitude bins have the same quality.
 *                               Introduce bin size;
 *                               If there are more than 30 stars, clip by 3 sigma; otherwise clip by 4 sigma
 *                               Add two new panels and add statistics of bright stars with good magdep_bin_quality.
 * Nov 25, 2011 Edward J. Los - Add MAGDEP_V2 changes described in Sumin Tang memorandum of Fri 11/18/11 5:27 PM
 *                               Increase MAGCAL_MAGDEP_QUALITY_LIMIT from 1.0 to 2.5
 *                               Decrease MAX_RMS2 from 0.7 to 0.5
 *                               Force the correction to zero for the faintest bin
 *                               Increase MIN_EXPAND from 10 to 20
 *                               Do not expand the bin size by more than a factor of two (2x2)
 *                               Increase BIN_THRESHOLD from 3000 to 5000
 *                               Increase MIN_SATURATED from 1000 to 5000
 *                               Decrease MIN_BIN from 20 to 15
 *                               Remove the smoothing over the magnitude bins
 *                               Replace magcal_magdep_rms with the clipped RMS value when magdep_bin_magcor is calculated
 *                               V2: MAGCAL_MAGDEP_QUALITY_LIMIT 2.500000 MAX_RMS2 0.500000 MIN_EXPAND 20 BIN_THRESHOLD 5000 MIN_SATURATED 5000 MIN_BIN 15 
 * Nov 25, 2011 Edward J. Los - Add MAGDEP_V3 changes described in Sumin Tang memorandum of Mon 11/28/11 1:47 PM
 *                               Decrease MAGCAL_MAGDEP_QUALITY_LIMIT from 2.5 to 2.0
 *                               Decrease BIN_THRESHOLD from 5000 to 3000
 *                               Increase MIN_BIN from 15 to 20
 *                               Decrease MIN_SATURATED from 5000 to 2000
 *                               Allow bin expansion to 4x4
 *                               V3: MAGCAL_MAGDEP_QUALITY_LIMIT 2.000000 MAX_RMS2 0.500000 MIN_EXPAND 20 BIN_THRESHOLD 3000 MIN_SATURATED 2000 MIN_BIN 20 
 * Dec  9, 2011 Edward J. Los - Correct segmantation fault when nfaint2 is zero or less
 * Dec 16, 2011 Edward J. Los - Increase the maximum bin expansion to 6x6 in accordance with Sumin's memorandum of 2011-12-15 13:52:06
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
 * May 28, 2018 Edward J. Los - Add gaia 
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "cpgplot.h"

#include "pipelineutils.h"

double round(double x);        /* Needed for Odyssey build */
float fminf(float x, float y); /* Needed for Odyssey build */

#define MAX_BUFFER 512
#define MIN_MAG_THRESHOLD 0.5 /* mthre: minimum threshold spacing of magnitude bins */
#define MIN_STARS 20          /* nreq1: minimum average number of stars per magnitude and location bin */
#define MIN_NSTARS 10
#if (defined(MAGDEP_V2) || defined(MAGDEP_V3))
#define MAX_RMS2 0.5
#define MIN_EXPAND 20         /* nreq2: threshold at which the bin must be expanded spatially */
#else /* MAGDEP_V2 || MAGDEP_V3 */
#define MAX_RMS2 0.7
#define MIN_EXPAND 10         /* nreq2: threshold at which the bin must be expanded spatially */
#endif /* MAGDEP_V2 || MAGDEP_V3 */

#ifdef MAGDEP_V2
#define BIN_THRESHOLD 5000 /* nthre: threshold average number of stars in each magnitude bin */
#define MIN_SATURATED 5000   /* Need at least 1000 saturated stars */
#define MIN_BIN 15            /* nbin:  minimum number of spatial bins */
#else /* MAGDEP_V2 */
#define BIN_THRESHOLD 3000 /* nthre: threshold average number of stars in each magnitude bin */
#define MIN_BIN 20            /* nbin:  minimum number of spatial bins */
#ifdef MAGDEP_V3
#define MIN_SATURATED 2000   /* Need at least 2000 saturated stars */
#else /* MAGDEP_V3 */
#define MIN_SATURATED 1000   /* Need at least 1000 saturated stars */
#endif /* MAGDEP_V3 */
#endif /* MAGDEP_V2 */
#define GRID_STEP 0.01
#define CLIP_LIMIT 30    /* Thirty stars for 3-sigma clipping, else use 4-sigma clipping */
#define PLOT_RATIO 1     /* Plot ratio of final dmag and original dmag; else plot non-smoothed correction vs original correction */

#define SELECTED_SATURATED   1
#define SELECTED_UNSATURATED 2
#define SELECTED_BRIGHT      4

typedef struct _starimage {
  int NUMBER;         /* Sextractor reference number NOTE: Both _starimage and _blend must have NUMBER as the first entry */
  int BFLAGS;          /* flags word */
  int cal_local;       /* This entry has been selected for output (begins with cal_local flag) */
  int npoints_local;  /* Number of points used for local calibration */
  int spatial_bin;    /* Lowess magnitude spatial bin  number */
  int local_bin;      /* Local calibration bin number */
  int selected;          /* Selection flags */
  double X_IMAGE;        /* Sextractor X location in pixels */
  double Y_IMAGE;        /* Sextractor Y location in pixels */
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
  double Stdmag;       /* GSC2.3.2 magnitude */
  double color;        /* GSC2.3.2 color */
  double MAG_ISO;     /* Sextractor isophotonic magnitude */
  double magcal_iso;  /* Lowess magnitude calibration */
  double magcal_iso_rms; /* Lowess magnitude calibration error */
  double magcal_local; /* Local magnitude calibration */
  double magcal_local_error; /* Local magnitude calibration error */
  double magcor_local; /* zout correction */
  double extinction;    /* Extinction correction */
  double limiting_mag; /* iso limiting magnitude */
  double limiting_mag_local; /* locally corrected limiting magnitude */
  double dmag;
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
} STARIMAGE,*PSTARIMAGE;




extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
double prctile(double *vector,int vectorcount,double pct)
{
  int index1;
  int index2;
  double pct1;
  double pct2;
  double pct3;
  double result;
  if (vectorcount <= 0) {
    return(0.0);
  }
  qsort(vector,vectorcount,sizeof(double),dcmp);
  index1 = ((1.0*vectorcount*pct)/100.0) -1.5;
  index2 = index1+1;
  if (index2 >= vectorcount) {
    return(vector[vectorcount-1]);
  } else if (index1 < 0) {
    return(vector[0]);
  } else {
    /* Here we interpolate */
    pct1 = (100.0 *(0.5+index1))/(1.0*vectorcount);
    pct2 = (100.0 *(0.5+index2))/(1.0*vectorcount);
    pct3 = 1.0*pct;
    result = vector[index1] + (((vector[index2]-vector[index1])*(pct3-pct1))/(pct2-pct1));
    return(result);

  }
  

}
void moving(double* invec,double* outvec,int veclen,int span)
{
  int index;
  int index2;
  int tmpspan;
  if ((span & 1) == 0) { /* Span must be odd */
    span--;
  }
  if (span != 3) {
    printf("ERROR: tested only for span = 3\n");
    exit(-1);
  }
  if (span == 0) {
    printf("ERROR: span is 0 in moving\n");
    exit(-1);
  }
  for (index = 0; index < veclen; index++) {
    outvec[index] = 0;
    tmpspan = span/2;
    if (index < tmpspan) {
      tmpspan = index;
    }
    if ((veclen-1-index) < tmpspan) {
      tmpspan = (veclen-1-index);
    }
    for (index2 = (index-tmpspan); index2 <= (index+tmpspan); index2++) {
      outvec[index] += invec[index2];
    }
    outvec[index] = outvec[index]/(2*tmpspan+1);  

  }



}
void interp1(double *xsample,double *ysample,int samplelength,double *xgrid,double *ygrid,int gridlength)
{
  int gridindex;
  int sampleindex;
  double xvalue;
  for (gridindex = 0; gridindex < gridlength; gridindex++) {
    xvalue = xgrid[gridindex];
    if (xvalue <= xsample[0]) {
      ygrid[gridindex] = ysample[0];
    } else if (xvalue >= xsample[samplelength-1]) {
      ygrid[gridindex] = ysample[samplelength-1];
    } else {
      for (sampleindex = 0; sampleindex < (samplelength-1); sampleindex++) {
        if ((xvalue >= xsample[sampleindex]) &&
            (xvalue < xsample[sampleindex+1])) {
          break;
        }
      }
      if (sampleindex == (samplelength-1)) {
        printf("ERROR: interp1 max sampleindex error\n");
        exit(-1);
      }
      ygrid[gridindex] = ysample[sampleindex]+ ((xvalue-xsample[sampleindex])*(ysample[sampleindex+1] - ysample[sampleindex])/(xsample[sampleindex+1] - xsample[sampleindex]));
    }
  }
}

int main(int argc,char *argv[])
{
  time_t startTime;
  time_t curTime;
  char title[MAX_BUFFER];
  char title2[MAX_BUFFER];
  char calfile[MAX_BUFFER];
  char fileroot[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  char gridfile[MAX_BUFFER];
  char plotdev[MAX_BUFFER];
  FILE *outHandle = NULL;
  FILE *calfile_handle = NULL;
  FILE *gridfile_handle = NULL;
  char *ingestDirectory;
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  int catalogNumber =  CATALOG_GSC232;
  char *plotValue;
  int plotFlag = 0;
  int verbose = 0;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int solutionNumber = -1;
  int binning;
  int rotation;
  PSTARIMAGE pSextractor;
  PSTARIMAGE pSextractor2;
  int index;
  int index2;

  File local_handle = NULL;
  char local_name[MAX_BUFFER];
  TableHead local_header = NULL;
  PSTARIMAGE local_table = NULL;
  size_t local_nrecs = 0;
  int local_index;
  PSTARIMAGE pLocal;
  int nbright; /* preliminary number of mag bins in the bright part */
  int nbright2;
  int ibright;
  double increment;
  double pct;
  double *Ybright = NULL;
  double *dybright = NULL;
  double *Ybright2 = NULL;
  double *dybright2 = NULL;
  double bmag1 = 0;
  double bmag2;
  int ntest;
  int ntest2;
  int itest;
  int nfaint;
  int ifaint;
  int nfaint2;
  int nfaint2x;
  double *Yfaint = NULL;
  double *Yfaint2 = NULL;
  int imag;
  int *nmag = NULL;
  double nmagmedian;
  int rejectedCount = 0;
  int StdmagRejected = 0;
  double X_IMAGE_MIN;
  double X_IMAGE_MAX;
  double Y_IMAGE_MIN;
  double Y_IMAGE_MAX;
  int kxb1;
  int kyb1;
  int ix;
  int iy;
  double *xrange;
  double *yrange;
  double *Zrms1;
  double *Z0;
  double *Z1;
  double *Z2;
  double *M1;
  double *X1;
  double *Zmag;
  double *Zfinalgrid;
  double *maggrid;
  float *xplot;
  float *yplot;
  float x_max;
  float x_min;
  float y_max;
  float y_min;
  double *M2;
  double clipFactor;

  int nthre = BIN_THRESHOLD;         /* nthre: threshold average number of stars in each magnitude bin */
  double mthre =  MIN_MAG_THRESHOLD; /* mthre: minimum threshold spacing of magnitude bins */
  int nbin = MIN_BIN;                /* nbin:  minimum number of spatial bins */
  int nreq1 =  MIN_STARS;            /* nreq1: minimum average number of stars per magnitude and location bin */
  int nreq2 = MIN_EXPAND;            /* nreq2: threshold at which the bin must be expanded spatially */

  double *vector1 = NULL;
  double *vector2 = NULL;
  double *vector3 = NULL;
  int vector1count = 0;
  int vector2count = 0;
  int vector3count = 0;

  double Y1;
  double Y2;
  double Y;
  
  double thedmag_median;
  double thedmag_rms;
  double final_median = 0;
  double final_rms = 0;
  double orig_median = 0;
  double orig_rms = 0;
  double const_median = 0;
  double const_rms = 0;
  double Zrms1_sum;
  double Zrms2_sum;
  double thedmag_median2;
  double thedmag_rms2;
  double Stdmag_median2;
  double Stdmag_rms2;
  double min_magcal_local;
  double max_magcal_local;
  int grid_length;

  int magdep_bin;  /* magdep_bin = ix + (kxb * iy) + (kxb * kyb * imag) 
                      where ix is the x-bin (0 to kxb-1) 
                      iy ix the y-bin (0 to ky-1)
                      imag is the magnitude bin (0 to kmag-1)
                   */
                      

  int *intvec1 = NULL;
  int *intvec2 = NULL;
  int *intvec3 = NULL;
  int *intvec4 = NULL;
  int *intvec5 = NULL;
  int intvec1count = 0;
  int intvec2count = 0;
  int intvec3count = 0;
  int intvec4count = 0;
  int intvec5count = 0;
  int nstar;

  PMAGDEPCORRECTION pMagdepTable = NULL;
  PMAGDEPCORRECTION pMagdep = NULL;
  MAGDEPLIMITS magdepLimits;
  PMAGDEPLIMITS pMagdepLimits = &magdepLimits;
  int goodBinCount = 0;
  int goodStarCount = 0;
  double debug_value;

  int temp_bin;
  double temp_magdep_rms;
  double temp_bin_magcor;
  int goodPlotCount;
  int totalPlotCount;
  int calfileopen = 0;


  memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
  time(&startTime);
  plotValue = getenv("DASCH_PLOT");
  if (plotValue == NULL) {
    printf("ERROR: DASCH_PLOT is not defined\n");
    exit(-1);
  }
  if (strstr(plotValue,"YES") != NULL) {
    plotFlag = 1;
  }
#if 0
  printf("ERROR: plotFlag is 1 in magdepcalibrate\n");
  plotFlag = 1;
#endif


  local_name[0] = 0;
  calfile[0] = 0;
  fileroot[0] = 0;
  qualifier[0] = 0;

  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    printf("ERROR: DASCH_INGEST is not defined\n");
    return(-1);
  }

 
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      printf("ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {


        case 'c': /* output calibration file name */
        case 'C':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(calfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'q': /* file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(qualifier,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
              errorFlag = 1;
            }
            if (strstr(qualifier,"kepler") != NULL) {
              catalogNumber = CATALOG_KEPLER;
            }
            if (strstr(qualifier,"apass") != NULL) {
              catalogNumber = CATALOG_APASS;
            }            
            if (strstr(qualifier,"atlas") != NULL) {
              catalogNumber = CATALOG_ATLAS;
            }            
            if (strstr(qualifier,"gaia") != NULL) {
              catalogNumber = CATALOG_GAIA;
            }
            if (strstr(qualifier,"experimental") != NULL) {
              catalogNumber = CATALOG_EXPERIMENTAL;
            }
           
          }
          break;

        case 'i': /* local file name */
        case 'I':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(local_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'r': /* file root */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(fileroot,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;





        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */
  if (calfile[0] == 0) {
    printf("ERROR: No output calibration filename was specified\n");
    errorFlag = 1;
  }
  if (local_name[0] == 0) {
    printf("ERROR: No local filename was specified\n");
    errorFlag = 1;
  }
  if (fileroot[0] == 0) {
    printf("ERROR: No file root was specified\n");
    errorFlag = 1;
  }


  /* Open the local file */
  local_handle = Open(local_name,"r");
  if (local_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the local file %s\n",local_name);
  } else {
    if (verbose) {
      printf("Found local file %s\n",local_name);
    }
  }



  
  if (ParseFilename(fileroot,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
    printf("ERROR: Failed to parse %s\n",fileroot);
    errorFlag = 1;
  }


  if (errorFlag) {
    printf("Usage: magdepcalibrate -r <file root> \n");
    printf("                      -i <local calibration file> \n");
    printf("                      -c  <output calibration file>\n");
    printf("                      -o  <output file>\n");
    printf("                      -v  verbose\n");

    return(-1);
  }

  if (verbose) {
    printf("magdepcalibrate of %s %s MAGCAL_MAGDEP_QUALITY_LIMIT %f MAX_RMS2 %f MIN_EXPAND %d BIN_THRESHOLD %d MIN_SATURATED %d MIN_BIN %d using %s\n",
           __DATE__,__TIME__,
           MAGCAL_MAGDEP_QUALITY_LIMIT,
           MAX_RMS2,
           MIN_EXPAND,
           BIN_THRESHOLD,
           MIN_SATURATED,
           MIN_BIN,
           local_name);
  }

  /* Now read in the local file */

  local_header = table_header(local_handle,TABLE_PARSE);
  if (local_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",local_name);
    return(-1);
  }

  local_table = table_loadva(local_handle,
                             &local_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(STARIMAGE),
                             &local_nrecs,
                             TblInt,"NUMBER" ,TblOff(PSTARIMAGE,NUMBER),
                             TblInt,"cal_local"  ,TblOff(PSTARIMAGE,cal_local),
                             TblDbl,"refmag",TblOff(PSTARIMAGE,Stdmag),  
                             TblDbl,"MAG_ISO",TblOff(PSTARIMAGE,MAG_ISO),
                             TblDbl,"magcal_iso"  ,TblOff(PSTARIMAGE,magcal_iso),
                             TblDbl,"magcal_iso_rms"  ,TblOff(PSTARIMAGE,magcal_iso_rms),
                             TblDbl,"magcal_local_error"  ,TblOff(PSTARIMAGE,magcal_local_error),
                             TblDbl,"extinction"  ,TblOff(PSTARIMAGE,extinction),
                             TblDbl,"magcor_local"  ,TblOff(PSTARIMAGE,magcor_local),
                             TblDbl,"limiting_mag"  ,TblOff(PSTARIMAGE,limiting_mag),
                             TblDbl,"X_IMAGE",TblOff(PSTARIMAGE,X_IMAGE),
                             TblDbl,"Y_IMAGE",TblOff(PSTARIMAGE,Y_IMAGE),
                             TblInt,"BFLAGS"  ,TblOff(PSTARIMAGE,BFLAGS),
                             TblInt,"npoints_local"  ,TblOff(PSTARIMAGE,npoints_local),
                             TblInt,"spatial_bin"  ,TblOff(PSTARIMAGE,spatial_bin),
                             TblInt,"local_bin"  ,TblOff(PSTARIMAGE,local_bin),
                             TblDbl,"color",TblOff(PSTARIMAGE,color),  
                             TblDbl,"ra"   ,TblOff(PSTARIMAGE,ra),
                             TblDbl,"dec"  ,TblOff(PSTARIMAGE,dec),
                             TblDbl,"magcal_local"  ,TblOff(PSTARIMAGE,magcal_local),
                             TblBuf,"REF"    ,TblOff(PSTARIMAGE,REF),MAX_REF,
                             0,"end",0);
  if (local_table == NULL) {
    printf("ERRORD: Failed to read table for %s\n",local_name);
    return(-1);
  }
  if (verbose) {
    printf("read %d records for %s\n",local_nrecs,local_name);
  }
  vector1 = (double *)calloc(local_nrecs,sizeof(double));
  vector2 = (double *)calloc(local_nrecs,sizeof(double));
  vector3 = (double *)calloc(local_nrecs,sizeof(double));
  if ((vector1 == NULL) ||
      (vector2 == NULL) ||
      (vector3 == NULL)) {
    printf("ERROR: failed to allocate vector1 %x or vector2 %x or vector3 %x\n",vector1,vector2,vector3);
  }





  /* Cast out values where Stdmag is less than zero */
  index = 0;
  while (index < local_nrecs) {
    pSextractor = &local_table[index];
#if 0
    if (strcmp(pSextractor->REF,"N0133233560") == 0) {
      printf("At ref %s\n",pSextractor->REF);
    }
#endif
    if (pSextractor->Stdmag <= 0) {
      StdmagRejected++;
    }

    if (
        (pSextractor->Stdmag <= 0)) {
      if (index < (local_nrecs-1)) {
        pSextractor2 = &local_table[local_nrecs-1];
        memcpy(pSextractor,pSextractor2,sizeof(STARIMAGE));
      }
      local_nrecs--;
      rejectedCount++;
    } else {
      index++;
    }
  }


  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
#if 0
    if (pSextractor->NUMBER == 9) {
      printf("At NUMBER %d\n",pSextractor->NUMBER);
    }
#endif

    pSextractor->limiting_mag_local = pSextractor->limiting_mag - pSextractor->magcor_local;
    pSextractor->dmag = pSextractor->magcal_local-pSextractor->Stdmag;
    if (index == 0) {
      min_magcal_local = pSextractor->magcal_local;
      max_magcal_local = pSextractor->magcal_local;
    } else {
      if (pSextractor->magcal_local < min_magcal_local) {
        min_magcal_local = pSextractor->magcal_local;
      }
      if (pSextractor->magcal_local > max_magcal_local) {
        max_magcal_local = pSextractor->magcal_local;
      }
    }

    /* Line 58 */
    if ((pSextractor->Stdmag < (pSextractor->limiting_mag_local - MAX_LIMITING_MAG)) &&
        ((pSextractor->BFLAGS & (1 << FILTER_BFLAG_PSFSATURATED)) != 0)) {
      pSextractor->selected |= SELECTED_SATURATED;
      vector1[vector1count] = pSextractor->Stdmag;
      vector1count++;
    }
    /* Line 59 */
    if ((pSextractor->Stdmag < (pSextractor->limiting_mag_local + MAX_LIMITING_MAG)) &&
        ((pSextractor->BFLAGS & (1 << FILTER_BFLAG_PSFSATURATED)) == 0)) {
      pSextractor->selected |= SELECTED_UNSATURATED;
      vector2[vector2count] = pSextractor->Stdmag;
      vector2count++;
    }
    /* V3.6.5 - discard unreasonable colors */
    switch (catalogNumber) {
    case CATALOG_GSC232:
    case CATALOG_EXPERIMENTAL:
      if ((pSextractor->color < MIN_GSC_COLOR) ||
          (pSextractor->color > MAX_GSC_COLOR)) {
        printf("ERROR: Stale GSC color value %f\n",pSextractor->color);
        pSextractor->color = 99.0;
        exit(-1);

      }
      break;
    case CATALOG_KEPLER:
      if ((pSextractor->color < MIN_KEPLER_COLOR) ||
          (pSextractor->color > MAX_KEPLER_COLOR)) {
        printf("ERROR: Stale KIC color value %f\n",pSextractor->color);
        pSextractor->color = 99.0;
        exit(-1);
      }
      break;
    case CATALOG_APASS:
      if ((pSextractor->color < MIN_APASS_COLOR) ||
          (pSextractor->color > MAX_APASS_COLOR)) {
        printf("ERROR: Stale APASS color value %f\n",pSextractor->color);
        pSextractor->color = 99.0;
        exit(-1);
      }
      break;
    case CATALOG_ATLAS:
      if ((pSextractor->color < MIN_ATLAS_COLOR) ||
          (pSextractor->color > MAX_ATLAS_COLOR)) {
        printf("ERROR: Stale ATLAS color value %f\n",pSextractor->color);
        pSextractor->color = 99.0;
        exit(-1);
      }
      break;
    case CATALOG_GAIA:
      if ((pSextractor->color < MIN_GAIA_COLOR) ||
          (pSextractor->color > MAX_GAIA_COLOR)) {
        printf("ERROR: Stale GAIA color value %f\n",pSextractor->color);
        pSextractor->color = 99.0;
        exit(-1);
      }
      break;
    default:
      printf("ERROR: Unknown catalog in prepare_octave\n");
      exit(-1);
      break;
    }

    

  }
  Y1 = prctile(vector1,vector1count,95.); /* 95% percentile */
  Y2 = prctile(vector2,vector2count,5.);  /*  5% percentile */
  Y = (Y1+Y2)/2.0;
#if 0
  if (verbose) {
    printf("vector1count %d, vector2count %d, Y1 %f, Y2 %f\n",vector1count,vector2count,Y1,Y2);
  }
#endif
  vector1count = 0; /* line 64 */
  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
    if ((pSextractor->Stdmag < Y) && 
        (pSextractor->Stdmag < pSextractor->limiting_mag_local)) {
      vector1[vector1count] = pSextractor->Stdmag;
      vector1count++;
    }
  }
  if (vector1count < MIN_SATURATED) {
    vector1count = 0; /* Line 67 */
    for (index = 0; index < local_nrecs; index++) {
      pSextractor = &local_table[index];
      if (pSextractor->Stdmag < pSextractor->limiting_mag_local) {
        vector1[vector1count] = pSextractor->Stdmag;
        vector1count++;
      }
    }
    if (vector1count < MIN_SATURATED) {
      printf("ERRORD: magdepcalibrate has only %d stars in %s%s\n",vector1count,fileroot,qualifier);
      exit(-1);
    }
    qsort(vector1,vector1count,sizeof(double),dcmp);
    Y = vector1[MIN_SATURATED-1];
     

  }
  vector1count = 0; /* Line 71 */
  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
    if (pSextractor->Stdmag < Y) {
      pSextractor->selected |= SELECTED_BRIGHT;
      vector1[vector1count] = pSextractor->Stdmag;
      vector1count++;
    }
  }
  nbright = round(((1.0*vector1count)/(1.0*nthre)) + 1.0);
  Ybright = (double *)calloc(nbright,sizeof(double));
  Ybright2 = (double *)calloc(nbright,sizeof(double));
  dybright = (double *)calloc(nbright,sizeof(double));
  dybright2 = (double *)calloc(nbright,sizeof(double));
  increment = 100.0/nbright;
  ntest = -1;
  for (ibright = 0; ibright < nbright; ibright++) {
    pct = increment + (ibright*increment);
    Ybright[ibright] = prctile(vector1,vector1count,pct);
    Ybright2[ibright] = round(Ybright[ibright]*10)/10;
    if (ibright > 0) {
      dybright[ibright-1] = Ybright[ibright] - Ybright[ibright-1];
      dybright2[ibright-1] = round(dybright[ibright-1]*10)/10;
      if (dybright[ibright-1] < mthre) {
        dybright2[ibright-1] = mthre;
        if (ntest < 0) {
          ntest = ibright-1;
        }
        ntest2 = ibright-1;
      }
#if 0
      printf("ibright %d Ybright %f dybright %f dybright2 %f\n",ibright,Ybright[ibright],dybright[ibright-1],dybright2[ibright-1]);
#endif

    }
  }
  if (ntest >= 0) {
    ntest = round((Y-Ybright[ntest])/mthre) + ntest -2;
    if (ntest >= 0) {
      if (ntest > ntest2) {
        ntest = ntest2;        
      }
      for (itest = 0; itest < ntest+1; itest++) {
        Ybright2[itest+1] = Ybright2[itest]+dybright2[itest];
        
      }
      nbright2 = ntest+2;
    } else {
      nbright2 = 1;
    }
  } else {
    nbright2 = nbright; /* Line 93 */
  }
  for (ibright = 0; ibright < nbright2; ibright++) {
    if (bmag1 < Ybright2[ibright]) {
      bmag1 = Ybright2[ibright];
    }
  }
  vector2count = 0; /* Line 96 */
  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
    vector2[vector2count] = pSextractor->limiting_mag_local;
    vector2count++;
    
  }
  bmag2 = prctile(vector2,vector2count,50.0);
  vector2count = 0; /* Line 97 */
  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
    if ((pSextractor->Stdmag >= bmag1) && 
        (pSextractor->Stdmag < bmag2)) {
      vector2[vector2count] = pSextractor->Stdmag;
      vector2count++;
    }
  }
  /* Line 99 - a factor of 2 added because of a bug in that line 99 should have filtered by myflag. 
   *           However, the mean ratio of total/good is 3 for the examples given, so a factor of 
   *           1.5 is more appropriate the actual range is more like 2.1 to 23
   */
  nfaint = ceil((1.5*vector2count)/nthre);
  nfaint2 = nfaint;
  Yfaint = (double *)calloc(nfaint,sizeof(double));
  Yfaint2 = (double *)calloc(nfaint,sizeof(double));
  increment = 100.0/nfaint;
  for (ifaint = 0; ifaint < nfaint; ifaint++) {
    pct = increment + (ifaint*increment);
    Yfaint[ifaint] = (round(prctile(vector2,vector2count,pct)*10.0))/10.0;
    Yfaint2[ifaint] = Yfaint[ifaint];
  }
  if (nfaint > round((bmag2-bmag1)/mthre)) {
    free(Yfaint2);
    nfaint2 = ceil((bmag2-bmag1)/0.5)+ 2;
    nfaint2x = nfaint2;
    if (nfaint2 < 0) {
      nfaint2 = 0;
    }
    Yfaint2= (double *)calloc(nfaint2+1,sizeof(double));
    nfaint2 = 0;
    Yfaint2[0] = bmag1+0.5;
    while(1) {
      if (Yfaint2[nfaint2]+0.5 > bmag2) {
        break;
      }
      if ((nfaint2+1) >= nfaint2x) {
        printf("ERROR: nfaint2x allocation %d is too small for %s%s\n",nfaint2x,fileroot,qualifier);
        exit(-1);
      }
      Yfaint2[nfaint2+1] = Yfaint2[nfaint2]+0.5;
      nfaint2++;      
    }
    nfaint2++;
  }
  /* Line 109 */
  pMagdepLimits->kmagb = 1 + nfaint2 + nbright2;

  pMagdepLimits->magdep_bin_edge  = (double *)calloc(pMagdepLimits->kmagb,sizeof(double));


#if 0
  printf("nbright2 %d, nfaint2 %d, pMagdepLimits->kmagb %d\n",nbright2,nfaint2,pMagdepLimits->kmagb);
#endif
  for (imag = 0; imag < pMagdepLimits->kmagb; imag++) {
    if (imag == 0) {
      pMagdepLimits->magdep_bin_edge[imag] = 0;
    } else if ((imag > 0) && (imag < nbright2+1)) {
      pMagdepLimits->magdep_bin_edge[imag] = Ybright2[imag-1];
    } else {
      pMagdepLimits->magdep_bin_edge[imag] = Yfaint2[imag-nbright2-1];
    }
#if 0
    printf("index %d magdep_bin_edge %f",imag,pMagdepLimits->magdep_bin_edge[imag]);
    if (imag < nfaint2) {
      printf(" Yfaint2 %f",Yfaint2[imag]);
    }
    if (imag < nbright2) {
      printf(" Ybright2 %f",Ybright2[imag]);
    }
    printf("\n");
#endif    

  }
  /* Line 112 */
  nmag = (int *)calloc(pMagdepLimits->kmagb-1,sizeof(int));
  for (imag = 0; imag < pMagdepLimits->kmagb-1; imag++) {
    for (index = 0; index < local_nrecs; index++) {
      pSextractor = &local_table[index];
      if ((pSextractor->Stdmag >= pMagdepLimits->magdep_bin_edge[imag]) &&
          (pSextractor->Stdmag < pMagdepLimits->magdep_bin_edge[imag+1])) {
        nmag[imag]++;
      }
    }
    vector1[imag] = 1.0*nmag[imag];
  }

  /* Line 120 */
  nmagmedian = prctile(vector1,nbright2,50.0);
  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
#if 0
    if (strcmp(pSextractor->REF,"N0133233560") == 0) {
      printf("At ref %s\n",pSextractor->REF);
    }
#endif
    if (index == 0) {
      X_IMAGE_MIN = pSextractor->X_IMAGE;
      X_IMAGE_MAX = pSextractor->X_IMAGE;
      Y_IMAGE_MIN = pSextractor->Y_IMAGE;
      Y_IMAGE_MAX = pSextractor->Y_IMAGE;
    } else {
      if (X_IMAGE_MIN >  pSextractor->X_IMAGE) {
        X_IMAGE_MIN = pSextractor->X_IMAGE;
      }
      if (X_IMAGE_MAX <  pSextractor->X_IMAGE) {
        X_IMAGE_MAX = pSextractor->X_IMAGE;
      }
      if (Y_IMAGE_MIN >  pSextractor->Y_IMAGE) {
        Y_IMAGE_MIN = pSextractor->Y_IMAGE;
      }
      if (Y_IMAGE_MAX <  pSextractor->Y_IMAGE) {
        Y_IMAGE_MAX = pSextractor->Y_IMAGE;
      }
    }
  }
  if (Y_IMAGE_MAX == Y_IMAGE_MIN) {
    printf("ERROR: No range in Y_IMAGE %f for %s%s\n",Y_IMAGE_MIN,fileroot,qualifier);
    exit(-1);
  }
  if (X_IMAGE_MAX == X_IMAGE_MIN) {
    printf("ERROR: No range in X_IMAGE %f for %s%s\n",X_IMAGE_MIN,fileroot,qualifier);
    exit(-1);
  }
  pMagdepLimits->kxb = round(sqrt(((nmagmedian/nreq1)*(X_IMAGE_MAX-X_IMAGE_MIN))/(Y_IMAGE_MAX-Y_IMAGE_MIN)));
  kxb1 = round(1.0*nbin*sqrt((X_IMAGE_MAX-X_IMAGE_MIN)/(Y_IMAGE_MAX-Y_IMAGE_MIN))); 
  if (kxb1 > pMagdepLimits->kxb) {
    pMagdepLimits->kxb = kxb1;
  }
  pMagdepLimits->kyb=round(sqrt(((nmagmedian/nreq1)*(Y_IMAGE_MAX-Y_IMAGE_MIN))/(X_IMAGE_MAX-X_IMAGE_MIN)));
  kyb1 = round(1.0*nbin*sqrt((Y_IMAGE_MAX-Y_IMAGE_MIN)/(X_IMAGE_MAX-X_IMAGE_MIN)));
  if (kyb1 > pMagdepLimits->kyb) {
    pMagdepLimits->kyb = kyb1;
  }
  pMagdepLimits->xcoord_magdep = (double *)calloc(pMagdepLimits->kxb,sizeof(double));
  pMagdepLimits->ycoord_magdep = (double *)calloc(pMagdepLimits->kyb,sizeof(double));

  /* line 125  for different magnitude stars */
  pMagdepLimits->dx = (X_IMAGE_MAX-X_IMAGE_MIN)/pMagdepLimits->kxb;
  pMagdepLimits->dy = (Y_IMAGE_MAX-Y_IMAGE_MIN)/pMagdepLimits->kyb;
  xrange = (double *)calloc(pMagdepLimits->kxb,sizeof(double));
  xrange[0] = X_IMAGE_MIN;
  for (ix = 1; ix < pMagdepLimits->kxb; ix++) {
    xrange[ix] = xrange[ix-1] + pMagdepLimits->dx;
  }
  yrange = (double *)calloc(pMagdepLimits->kyb,sizeof(double));
  yrange[0] = Y_IMAGE_MIN;
  for (iy = 1; iy < pMagdepLimits->kyb; iy++) {
    yrange[iy] = yrange[iy-1] + pMagdepLimits->dy;
  }

  /* Line 131: clear the following variables */
  pMagdepLimits->nrecs = pMagdepLimits->kxb*pMagdepLimits->kyb*(pMagdepLimits->kmagb-1);
  pMagdepTable = (PMAGDEPCORRECTION)calloc(pMagdepLimits->nrecs,sizeof(MAGDEPCORRECTION));

  Zrms1 = (double*)calloc(pMagdepLimits->nrecs,sizeof(double));

  Z0 = (double *)calloc((pMagdepLimits->kmagb-1),sizeof(double));
  Z1 = (double *)calloc((pMagdepLimits->kmagb-1),sizeof(double));
  Z2 = (double *)calloc((pMagdepLimits->kmagb-1),sizeof(double));
  M1 = (double *)calloc((pMagdepLimits->kmagb-1),sizeof(double));
  X1 = (double *)calloc((pMagdepLimits->kmagb-1),sizeof(double));

  intvec1 = (int *)calloc(local_nrecs,sizeof(int));
  intvec2 = (int *)calloc(local_nrecs,sizeof(int));
  intvec3 = (int *)calloc(local_nrecs,sizeof(int));
  intvec4 = (int *)calloc(local_nrecs,sizeof(int));
  intvec5 = (int *)calloc(local_nrecs,sizeof(int));

  Zmag       = (double*)calloc(pMagdepLimits->kmagb+1,sizeof(double));
  M2         = (double*)calloc(pMagdepLimits->kmagb+1,sizeof(double));


  if (plotFlag) {
    sprintf(plotdev,"%s/%s%s_magdep1.ps/ps",ingestDirectory,fileroot,qualifier);
    cpgopen(plotdev);
    cpgslw (1);
    cpgsch (1.5);
    cpgsubp(6,6);
  }
  


  for (ix = 0; ix < pMagdepLimits->kxb; ix++) { /* spatial bin in the X direction Line 139 */
    pMagdepLimits->xcoord_magdep[ix] = xrange[ix] + (pMagdepLimits->dx/2);
    for (iy = 0; iy < pMagdepLimits->kyb; iy++) { /*  spatial bin in the Y direction */
      pMagdepLimits->ycoord_magdep[iy] = yrange[iy] + (pMagdepLimits->dy/2);
      for (imag = 0; imag < (pMagdepLimits->kmagb-1); imag++) {
        Z0[imag] = 0;
        Z1[imag] = 0;
        Z2[imag] = 0;
        M1[imag] = 0;
        X1[imag] = 0;
      }
      for (imag = 0; imag < (pMagdepLimits->kmagb-1); imag++) {
        magdep_bin = ix + (pMagdepLimits->kxb*iy) + (pMagdepLimits->kxb*pMagdepLimits->kyb*imag);
        pMagdep = &pMagdepTable[magdep_bin];
        pMagdep->xcoord_magdep = pMagdepLimits->xcoord_magdep[ix];
        pMagdep->ycoord_magdep = pMagdepLimits->ycoord_magdep[iy];

        intvec1count = 0;
        intvec2count = 0;
        intvec3count = 0;
        intvec4count = 0;
        intvec5count = 0;
        for (index = 0; index < local_nrecs; index++) {
          pSextractor = &local_table[index];
          if ((pSextractor->Stdmag < pMagdepLimits->magdep_bin_edge[imag]) ||
              (pSextractor->Stdmag >= pMagdepLimits->magdep_bin_edge[imag+1])) {
            continue;
          }
          intvec5count++;
          if ((fabs(pSextractor->X_IMAGE-pMagdepLimits->xcoord_magdep[ix]) < (3*pMagdepLimits->dx)) &&
              (fabs(pSextractor->Y_IMAGE-pMagdepLimits->ycoord_magdep[iy]) < (3*pMagdepLimits->dy))) {
            intvec4[intvec4count] = index;
            intvec4count++;
          } 
          if ((fabs(pSextractor->X_IMAGE-pMagdepLimits->xcoord_magdep[ix]) < (2*pMagdepLimits->dx)) &&
              (fabs(pSextractor->Y_IMAGE-pMagdepLimits->ycoord_magdep[iy]) < (2*pMagdepLimits->dy))) {
            intvec3[intvec3count] = index;
            intvec3count++;
          } 
          if ((fabs(pSextractor->X_IMAGE-pMagdepLimits->xcoord_magdep[ix]) < (pMagdepLimits->dx)) &&
              (fabs(pSextractor->Y_IMAGE-pMagdepLimits->ycoord_magdep[iy]) < (pMagdepLimits->dy))) {
            intvec2[intvec2count] = index;
            intvec2count++;
          } 
          if ((fabs(pSextractor->X_IMAGE-pMagdepLimits->xcoord_magdep[ix]) < (pMagdepLimits->dx/2)) &&
              (fabs(pSextractor->Y_IMAGE-pMagdepLimits->ycoord_magdep[iy]) < (pMagdepLimits->dy/2))) {
            intvec1[intvec1count] = index;
            intvec1count++;
#if 0
            if ((ix == 17) && (iy == 7) && (imag == 6)) {
              printf("count %d X_IMAGE %f Y_IMAGE %f Stdmag %f BFLAGS 0x%08x %d\n",
                     intvec1count,
                     pSextractor->X_IMAGE,
                     pSextractor->Y_IMAGE,
                     pSextractor->Stdmag,
                     pSextractor->BFLAGS,
                     pSextractor->BFLAGS);
            }
#endif
          }
        }
#if 0
        if ((ix == 17) && (iy == 7) && (imag == 6)) {
          printf("ix %d iy %d imag  %d\n",ix,iy,imag);
        }
#endif
        if (intvec1count >= nreq1) {
          nstar = intvec1count;
          pMagdep->magdep_bin_size = 1;
        } else if (intvec2count >= nreq2) {
          nstar = intvec2count;
          pMagdep->magdep_bin_size = 2;
          for (index = 0; index < nstar; index++) {
            intvec1[index] = intvec2[index];
          }
          intvec1count = intvec2count;
#ifndef MAGDEP_V2
        } else if (intvec3count >= nreq2) {
          nstar = intvec3count;
          pMagdep->magdep_bin_size = 4;
          for (index = 0; index < nstar; index++) {
            intvec1[index] = intvec3[index];
          }
          intvec1count = intvec3count;
#if ((!defined(MAGDEP_V2)) && (!defined(MAGDEP_V3b)))
        } else if (intvec4count >= nreq2) {
          nstar = intvec4count;
          pMagdep->magdep_bin_size = 6;
          for (index = 0; index < nstar; index++) {
            intvec1[index] = intvec4[index];
          }
          intvec1count = intvec4count;
#endif /* !MAGDEP_V2 && !MAGDEP_V3b */
#endif /* MAGDEP_V2 */
        } else {
          /* Not enough stars.  Line 184 */
#ifdef MAGDEP_V2
          nstar = intvec2count;
#else /* MAGDEP_V2 */
#ifdef MAGDEP_V3b
          nstar = intvec3count;
#else /* MAGDEP_V3b */
          nstar = intvec4count;
#endif /* MAGDEP_V3b */
#endif /* MAGDEP_V2 */
          intvec1count = 0;
          Zrms1[magdep_bin] = 99.0;
          pMagdep->magcal_magdep_rms = 99.0;
          Z0[imag] = 0.0;
          pMagdep->magdep_bin_size = 99;
          pMagdep->magdep_bin_median = (pMagdepLimits->magdep_bin_edge[imag]+pMagdepLimits->magdep_bin_edge[imag+1])/2.0;
        }
        if (nstar >= nreq2) {
          /* line 170: find the median of thedmag */
          vector1count = 0;
          for (index = 0; index < nstar; index++) {
            pSextractor = &local_table[intvec1[index]];
            vector1[vector1count] = pSextractor->dmag; /* thedmag */
            vector1count++;
          }
          CalcMedianAndRMS(vector1count,0,vector1,&thedmag_median,&thedmag_rms,0,0,0);
          /* intvec5 becomes iy2 in line 173 */
          intvec5count = 0;
          clipFactor = 3.0;
          if (nstar < CLIP_LIMIT) {
            clipFactor = 4.0;
          }

          for (index = 0; index < nstar; index++) {
            pSextractor = &local_table[intvec1[index]];
            if (fabs(pSextractor->dmag - thedmag_median) < (clipFactor*thedmag_rms)) {
              intvec5[intvec5count] = intvec1[index];
              intvec5count++;
            }
          }
          pMagdep->nstar_magdep = intvec5count;
          if (intvec5count > 1) {
            /* Line 177 */
            Zrms1_sum = 0.0;
            Zrms2_sum = 0.0;
            vector1count = 0;
            vector2count = 0;
            for (index = 0; index < intvec5count; index++) {
              pSextractor = &local_table[intvec5[index]];
              Zrms1_sum += sqr(pSextractor->dmag);
              Zrms2_sum += sqr(pSextractor->dmag - thedmag_median);
              vector1[vector1count] = pSextractor->dmag;
              vector1count++;
              vector2[vector2count] = pSextractor->Stdmag;
              vector2count++;
            }
            Zrms1[magdep_bin] = sqrt(Zrms1_sum/(intvec5count-1));
            CalcMedianAndRMS(vector1count,0,vector1,&thedmag_median2,&thedmag_rms2,0,0,0);
            Z0[imag] = thedmag_median2;
#if  (defined(MAGDEP_V2) || defined(MAGDEP_V3))
            pMagdep->magcal_magdep_rms = thedmag_rms2;
#else /* MAGDEP_V2 || MAGDEP_V3 */
            pMagdep->magcal_magdep_rms = sqrt(Zrms2_sum/(intvec5count-1));
#endif /* MAGDEP_V2 || MAGDEP_V3 */

            CalcMedianAndRMS(vector2count,0,vector2,&Stdmag_median2,&Stdmag_rms2,0,0,0);
            pMagdep->magdep_bin_median = Stdmag_median2;
#if  (defined(MAGDEP_V2) || defined(MAGDEP_V3))
            if (imag == (pMagdepLimits->kmagb-2)) {
              Zrms1[magdep_bin] = 99.0;
              pMagdep->magcal_magdep_rms = 99.0;
              Z0[imag] = 0.0;
              pMagdep->magdep_bin_median = (pMagdepLimits->magdep_bin_edge[imag]+pMagdepLimits->magdep_bin_edge[imag+1])/2.0;
            }
#endif /*  MAGDEP_V2 || MAGDEP_V3 */
          } else {
            Zrms1[magdep_bin] = 99.0;
            pMagdep->magcal_magdep_rms = 99.0;
            Z0[imag] = 0;
            pMagdep->magdep_bin_median = (pMagdepLimits->magdep_bin_edge[imag]+pMagdepLimits->magdep_bin_edge[imag+1])/2.0;
          }
          /* Line 191 */
          Z1[imag] = Z0[imag];
          if ((pMagdep->nstar_magdep < MIN_NSTARS) || 
              (pMagdep->magcal_magdep_rms > MAX_RMS2)) {
            Z1[imag] = 0;
          }

        } else { /* nstar > nreq  line 184*/
          Zrms1[magdep_bin] = 99.0;
          pMagdep->magcal_magdep_rms = 99.0;
          Z0[imag] = 0;
          pMagdep->magdep_bin_median = (pMagdepLimits->magdep_bin_edge[imag]+pMagdepLimits->magdep_bin_edge[imag+1])/2.0;
          pMagdep->nstar_magdep = nstar;
          Z1[imag] = 0;
          
        }
#if 0
        printf("ix %d iy %d imag %d, magdep_bin %d, nstar_magdep %d\n",ix,iy,imag,magdep_bin,pMagdep->nstar_magdep);
#endif        


      } /* imag index */
      /* Line 198: do a smoothing in magnitude bins  */
#if 0
      for (imag = 0; imag < (pMagdepLimits->kmagb-1); imag++) {
        printf("ix %2d iy %2d imag %2d Z0 %f Z1 %f\n",ix,iy,imag,Z0[imag],Z1[imag]);
      }
#endif
#if (defined(MAGDEP_V2) || defined(MAGDEP_V3))
      for (imag = 0; imag < (pMagdepLimits->kmagb-1); imag++) {
        Z2[imag] = Z1[imag];
      }
#else /* MAGDEP_V2 || MAGDEP_V3 */
      moving(Z1,Z2,pMagdepLimits->kmagb-1,3);
#endif /* MAGDEP_V2 || MAGDEP_V3 */
      Zmag[0] = Z2[0];
      Zmag[pMagdepLimits->kmagb] = Z2[pMagdepLimits->kmagb-2];
      M2[0] = min_magcal_local;
      M2[pMagdepLimits->kmagb] = max_magcal_local;
    
      for (imag = 0; imag < (pMagdepLimits->kmagb-1); imag++) {
        magdep_bin = ix + (pMagdepLimits->kxb*iy) + (pMagdepLimits->kxb*pMagdepLimits->kyb*imag);
        pMagdep = &pMagdepTable[magdep_bin];
        pMagdep->magdep_bin_magcor = Z2[imag];
        M1[imag] = pMagdep->magdep_bin_median;
        Zmag[imag+1] = Z2[imag];
        M2[imag+1] = M1[imag];
        if (pMagdep->magcal_magdep_rms == 0) {
          pMagdep->magdep_bin_quality = 0;
        } else {
          pMagdep->magdep_bin_quality = ((pMagdep->magdep_bin_magcor*sqrt(pMagdep->nstar_magdep))/(pMagdep->magcal_magdep_rms));
          if (fabs(pMagdep->magdep_bin_quality) >= MAGCAL_MAGDEP_QUALITY_LIMIT) {
            goodBinCount++;
          }
        }

      }
#if 0
      {
        int kkk;
        for (kkk = 0; kkk < pMagdepLimits->kmagb+1;kkk++) {
          printf(" %d M2 %f Zmag %f\n",kkk,M2[kkk],Zmag[kkk]);
        }
      }
#endif
      /* Line 202: linear interp */
      grid_length = (max_magcal_local-min_magcal_local)/GRID_STEP;
      Zfinalgrid = (double*)calloc(grid_length,sizeof(double));
      maggrid    = (double*)calloc(grid_length,sizeof(double));
      xplot    =   (float*)calloc(grid_length,sizeof(float));
      yplot    =   (float*)calloc(grid_length,sizeof(float));


      for (index = 0; index < grid_length; index++) {
        if (index == 0) {
          maggrid[index] = min_magcal_local;
        } else {
          maggrid[index] = maggrid[index-1]+GRID_STEP;
        }
      }
      interp1(M2,Zmag,pMagdepLimits->kmagb+1,maggrid,Zfinalgrid,grid_length);
      for (imag = 0; imag < pMagdepLimits->kmagb-1; imag++) {
        magdep_bin = ix + (pMagdepLimits->kxb*iy) + (pMagdepLimits->kxb*pMagdepLimits->kyb*imag);
        pMagdep = &pMagdepTable[magdep_bin];

        if (calfileopen == 0) {

          /* Open the calibration file */
          calfile_handle = Open(calfile,"wt");
          if (calfile_handle == NULL) {
            printf("ERROR: Failed to find the calibration file %s\n",calfile);
            exit(-1);
          } else {
            if (verbose) {
              printf("Found calibration file %s\n",calfile);
            }
            fprintf(calfile_handle,"ixb\tiyb\timagb\tnstar_magdep\tmagdep_bin_size\txcoord_magdep\tycoord_magdep\tmagdep_bin_edge\tmagdep_bin_median\tmagdep_bin_magcor\tmagcal_magdep_rms\tmagdep_bin_quality\n");
            fprintf(calfile_handle,"---\t---\t-----\t------------\t---------------\t-------------\t-------------\t---------------\t-----------------\t--------------\t-----------------\t------------------\n");
    
          }
          calfileopen = 1;
        }

        fprintf(calfile_handle,"%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                ix,iy,imag,pMagdep->nstar_magdep,pMagdep->magdep_bin_size,pMagdepLimits->xcoord_magdep[ix],pMagdepLimits->ycoord_magdep[iy],pMagdepLimits->magdep_bin_edge[imag+1],pMagdep->magdep_bin_median, pMagdep->magdep_bin_magcor,pMagdep->magcal_magdep_rms,pMagdep->magdep_bin_quality);
      }        
      if (plotFlag) {
        if (((ix+1) == round((1.0*pMagdepLimits->kxb)/4)) && ((iy+1) == round((pMagdepLimits->kyb*3.0)/4))) {
          cpgpanl(1,1);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title2,"X=1/4,y=3/4 (%d,%d)",ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title2);

          cpgline(grid_length,xplot,yplot);
        }
        if (((ix+1) == round((1.0*pMagdepLimits->kxb)/2)) && ((iy+1) == round((pMagdepLimits->kyb*3.)/4))) {
          cpgpanl(2,1);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
        
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title,"%s%s (%d,%d)",fileroot,qualifier,ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title);

          cpgline(grid_length,xplot,yplot);
        }
        if (((ix+1) == round((pMagdepLimits->kxb*3.)/4)) && ((iy+1) == round((pMagdepLimits->kyb*3.)/4))) {
          cpgpanl(3,1);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
        
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title2,"X=3/4,y=3/4 (%d,%d)",ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title2);

          cpgline(grid_length,xplot,yplot);
        }
        if (((ix+1) == round((pMagdepLimits->kxb*1.)/4)) && ((iy+1) == round((pMagdepLimits->kyb*2.)/4))) {
          cpgpanl(1,2);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
        
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title2,"X=1/4,y=2/4 (%d,%d)",ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title2);

          cpgline(grid_length,xplot,yplot);
        }
        if (((ix+1) == round((pMagdepLimits->kxb*2.)/4)) && ((iy+1) == round((pMagdepLimits->kyb*2.)/4))) {
          cpgpanl(2,2);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
        
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title2,"X=2/4,y=2/4 (%d,%d)",ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title2);

          cpgline(grid_length,xplot,yplot);
        }
        if (((ix+1) == round((pMagdepLimits->kxb*3.)/4)) && ((iy+1) == round((pMagdepLimits->kyb*2.)/4))) {
          cpgpanl(3,2);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
        
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title2,"X=3/4,y=2/4 (%d,%d)",ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title2);

          cpgline(grid_length,xplot,yplot);
        }
        if (((ix+1) == round((pMagdepLimits->kxb*1.)/4)) && ((iy+1) == round((pMagdepLimits->kyb*1.)/4))) {
          cpgpanl(1,3);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
        
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title2,"X=1/4,y=1/4 (%d,%d)",ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title2);

          cpgline(grid_length,xplot,yplot);
        }
        if (((ix+1) == round((pMagdepLimits->kxb*2.)/4)) && ((iy+1) == round((pMagdepLimits->kyb*1.)/4))) {
          cpgpanl(2,3);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
        
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title2,"X=2/4,y=1/4 (%d,%d)",ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title2);

          cpgline(grid_length,xplot,yplot);
        }
        if (((ix+1) == round((pMagdepLimits->kxb*3.)/4)) && ((iy+1) == round((pMagdepLimits->kyb*1.)/4))) {
          cpgpanl(3,3);
          if (Zmag[pMagdepLimits->kmagb] == Zmag[0]) {
            Zmag[pMagdepLimits->kmagb] = Zmag[0]+1;
            Zmag[0] = Zmag[0]-1;
          }
          for (index = 0; index < grid_length; index++) {
            xplot[index] = maggrid[index];
            yplot[index] = Zfinalgrid[index];
            if (index == 0) {
              x_max = xplot[index];
              x_min = xplot[index];
              y_max = yplot[index];
              y_min = yplot[index];
            } else {
              if (x_max < xplot[index]) {
                x_max = xplot[index];
              }
              if (x_min > xplot[index]) {
                x_min = xplot[index];
              }
              if (y_max < yplot[index]) {
                y_max = yplot[index];
              }
              if (y_min > yplot[index]) {
                y_min = yplot[index];
              }
            }
          }
          cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);
        
          cpgmtxt ("B",2.5,0.5,0.5,"mag");
          cpgmtxt ("L",2.0,0.5,0.5,"dmag-correction");
	  sprintf(title2,"X=3/4,y=1/4 (%d,%d)",ix,iy);
          cpgmtxt ("T",1.0,0.5,0.5,title2);

          cpgline(grid_length,xplot,yplot);
        }
      }
  
      if (maggrid != NULL) {
        free(maggrid);
        maggrid = NULL;
      }
      if (Zfinalgrid != NULL) {
        free(Zfinalgrid);
        Zfinalgrid = NULL;
      }
      if (xplot != NULL) {
        free(xplot);
        xplot = NULL;
      }
      if (yplot != NULL) {
        free(yplot);
        yplot = NULL;
      }
      

    } /* iy index */
  } /* ix index */
  /* Now get three statistics for comparison:
   *  1.  The original rms of dmag
   *  2.  The rms of dmag-magdep_bin_magcor, no interpolation
   *  3.  The rms of dmag-magdep_bin_magcor, with interpolation
   *  Do this only for bright stars coming from good quality bins
   */
  vector1count = 0;
  vector2count = 0;
  vector3count = 0;
  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
    magdep_bin = GetMagdepBin(pMagdepLimits,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&ix,&iy,&imag);
    pMagdep = &pMagdepTable[magdep_bin];
#if 0
    if ((ix == 11) && (iy == 13)) {
      printf("%f %f\n",pSextractor->Stdmag,GetMagdepBinMagcor(pMagdepLimits,pMagdepTable,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&temp_bin,&temp_magdep_rms,&temp_bin_magcor,0));
    }
#endif
    if (fabs(pMagdep->magdep_bin_quality) < MAGCAL_MAGDEP_QUALITY_LIMIT) {
      continue;
    }
    goodStarCount++;
    if (pSextractor->Stdmag > Y) {
      continue;
    }
    vector1[vector1count] = pSextractor->dmag;
    vector2[vector2count] = pSextractor->magcal_local - pSextractor->Stdmag - pMagdep->magdep_bin_magcor;
    vector3[vector3count] = pSextractor->magcal_local - pSextractor->Stdmag - GetMagdepBinMagcor(pMagdepLimits,pMagdepTable,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&temp_bin,&temp_magdep_rms,&temp_bin_magcor,1);
#if 0
    printf("index %d dmag %f const magcor %f interp magcor %f  vector %f %f %f\n",
	   vector1count,
	   pSextractor->dmag,
	   pMagdep->magdep_bin_magcor,
	   GetMagdepBinMagcor(pMagdepLimits,pMagdepTable,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&temp_bin,&temp_magdep_rms,&temp_bin_magcor,1),
	   vector1[vector1count],
	   vector2[vector2count],
	   vector3[vector3count]);
#endif 
    vector1count++;
    vector2count++;
    vector3count++;
    CalcMedianAndRMS(vector1count,0,vector1,&orig_median,&orig_rms,1,3.0,0);
    CalcMedianAndRMS(vector2count,0,vector2,&const_median,&const_rms,1,3.0,0);
    CalcMedianAndRMS(vector3count,0,vector3,&final_median,&final_rms,1,3.0,0);
  }
  

  if (plotFlag) {
    cpgiden();
    cpgend();

    sprintf(plotdev,"%s/%s%s_magdep2.ps/ps",ingestDirectory,fileroot,qualifier);
    cpgopen(plotdev);
    cpgslw (2);
    cpgsch (1.5);
    cpgsubp(3,3);
    if (local_nrecs > pMagdepLimits->nrecs) {
      xplot    =   (float*)calloc(local_nrecs,sizeof(float));
      yplot    =   (float*)calloc(local_nrecs,sizeof(float));
    } else {
      xplot    =   (float*)calloc(pMagdepLimits->nrecs,sizeof(float));
      yplot    =   (float*)calloc(pMagdepLimits->nrecs,sizeof(float));
    }
    cpgpanl(1,1);


    for (magdep_bin = 0; magdep_bin < (pMagdepLimits->nrecs); magdep_bin++) {
      pMagdep = &pMagdepTable[magdep_bin];
      xplot[magdep_bin] = pMagdep->magdep_bin_magcor;
      yplot[magdep_bin] = pMagdep->magcal_magdep_rms;
      if (yplot[magdep_bin] > 2.0) {
        yplot[magdep_bin] = 2.0;
      }
      if (magdep_bin == 0) {
        x_max = xplot[magdep_bin];
        x_min = xplot[magdep_bin];
        y_max = yplot[magdep_bin];
        y_min = yplot[magdep_bin];
      } else {
        if (x_max < xplot[magdep_bin]) {
          x_max = xplot[magdep_bin];
        }
        if (x_min > xplot[magdep_bin]) {
          x_min = xplot[magdep_bin];
        }
        if (y_max < yplot[magdep_bin]) {
          y_max = yplot[magdep_bin];
        }
        if (y_min > yplot[magdep_bin]) {
          y_min = yplot[magdep_bin];
        }
      }
    }
    cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);

    cpgmtxt ("B",2.5,0.5,0.5,"magdep_bin_magcor (dmag correction)");
    cpgmtxt ("L",2.0,0.5,0.5,"magcal_magdep_rms");
    sprintf(title,"%s%s",fileroot,qualifier);
    cpgmtxt ("T",1.0,0.5,0.5,title);

    cpgpt(grid_length,xplot,yplot,-1);

    cpgpanl(2,1);  


    for (magdep_bin = 0; magdep_bin < (pMagdepLimits->nrecs); magdep_bin++) {
      pMagdep = &pMagdepTable[magdep_bin];
      xplot[magdep_bin] = Zrms1[magdep_bin];
      if (xplot[magdep_bin] > 2.0) {
        xplot[magdep_bin] = 2.0;
      }
      yplot[magdep_bin] = pMagdep->magcal_magdep_rms;
      if (yplot[magdep_bin] > 2.0) {
        yplot[magdep_bin] = 2.0;
      }
#if 0
      printf("%f %f\n",xplot[magdep_bin],yplot[magdep_bin]);
#endif
      if (magdep_bin == 0) {
        x_max = xplot[magdep_bin];
        x_min = xplot[magdep_bin];
        y_max = yplot[magdep_bin];
        y_min = yplot[magdep_bin];
      } else {
        if (x_max < xplot[magdep_bin]) {
          x_max = xplot[magdep_bin];
        }
        if (x_min > xplot[magdep_bin]) {
          x_min = xplot[magdep_bin];
        }
        if (y_max < yplot[magdep_bin]) {
          y_max = yplot[magdep_bin];
        }
        if (y_min > yplot[magdep_bin]) {
          y_min = yplot[magdep_bin];
        }
      }
    }
#if 0
    printf("x_min %f x_max %f y_min %f y_max %f\n",x_min,x_max,y_min,y_max);
#endif
    cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);

    cpgmtxt ("B",2.5,0.5,0.5,"Zrms1 (correction)");
    cpgmtxt ("L",2.0,0.5,0.5,"magcal_magdep_rms");
    sprintf(title,"%s%s",fileroot,qualifier);
    cpgmtxt ("T",1.0,0.5,0.5,title);

    cpgpt(grid_length,xplot,yplot,-1);
  
    xplot[0] = 0;
    yplot[0] = 0;
    xplot[1] = fminf(x_max,y_max);
    yplot[1] = fminf(x_max,y_max);

    cpgline(2,xplot,yplot);


    cpgpanl(1,2);
    /* Histogram of nstar_magdep from 0 to 200 in increments of 2 */
    for (index = 0; index < 100; index++) {
      xplot[index] = (index*2.0)+1.0; /* bin centers */
      yplot[index] = 0;
    }
    x_min = 0;
    x_max = 200;
    y_min = 0;
    y_max = 0;

    for (magdep_bin = 0; magdep_bin < (pMagdepLimits->nrecs); magdep_bin++) {
      pMagdep = &pMagdepTable[magdep_bin];
      index = pMagdep->nstar_magdep/2;
      if (index >= 200) {
        continue;
      }
      yplot[index]++;
      if (yplot[index] > y_max) {
        y_max = yplot[index];
      }
    
    }
    cpgenv(x_min,x_max,y_min,y_max+0.1,0,0);

    cpgmtxt ("B",2.5,0.5,0.5,"Num of stars in the mag&spatial bins");
    cpgmtxt ("L",2.0,0.5,0.5,"Num of bins");
    sprintf(title,"%s%s",fileroot,qualifier);
    cpgmtxt ("T",1.0,0.5,0.5,title);
    cpgbin(100,xplot,yplot,1);


    cpgpanl(2,2);
    /* Histogram of magdep_bin_magcor./(magcal_magdep_rms./nstar_magdep.^0.5) from -20 to 20 in increments of 0.5 (80 bins) */
    for (index = 0; index < 100; index++) {
      xplot[index] = (index*0.5)-20.0+0.25; /* bin centers */
      yplot[index] = 0;
    }
    x_min = -20.;
    x_max =  20.;
    y_min = 0;
    y_max = 0;
    

    for (magdep_bin = 0; magdep_bin < (pMagdepLimits->nrecs); magdep_bin++) {
      pMagdep = &pMagdepTable[magdep_bin];
      index = ((pMagdep->magdep_bin_quality+20.0)/0.5);
#if 0
      printf("index %d  Value %f\n",index,pMagdep->magdep_bin_quality);
#endif
      if ((index < 0) ||  
          (index >= 200)) {
        continue;
      }
      yplot[index]++;
      if (yplot[index] > y_max) {
        y_max = yplot[index];
      }
    
    }
    cpgenv(x_min,x_max,y_min,y_max+0.1,0,0);

    cpgmtxt ("B",2.5,0.5,0.5,"(magdep_bin_magcor*(nstar_magdep^{0.5}))/magcal_magdep_rms");
    cpgmtxt ("L",2.0,0.5,0.5,"Num of bins");
    sprintf(title,"%s%s",fileroot,qualifier);
    cpgmtxt ("T",1.0,0.5,0.5,title);
    cpgbin(100,xplot,yplot,1);



#if 0  /* Dump the linear arrays */
    { 
      int indexFlag;
      index = 0;
      indexFlag = 1;
      printf("dx %f dy %f\n",pMagdepLimits->dx,pMagdepLimits->dy);
      while (indexFlag == 1) {
        indexFlag = 0;
        printf("index %d ",index);
        if (index < pMagdepLimits->kxb) {
          indexFlag = 1;
          printf("xcoord_magdep %f ",pMagdepLimits->xcoord_magdep[index]);
        }
        if (index < pMagdepLimits->kyb) {
          indexFlag = 1;
          printf("ycoord_magdep %f ",pMagdepLimits->ycoord_magdep[index]);
        }
        if (index < (pMagdepLimits->kmagb)) {
          indexFlag = 1;
          printf("magdep_bin_edge %f ",pMagdepLimits->magdep_bin_edge[index]);
        }
        printf("\n");
        if (indexFlag == 0) {
          break;
        }
        index++;
      }
    }
#endif
    debug_value = pMagdepLimits->xcoord_magdep[0];
    cpgpanl(1,3);
#if PLOT_RATIO
   /* Plot the corrected dmag against the original dmag */
    goodPlotCount = 0;
    totalPlotCount = 0;
    index2 = 0;
    for (index = 0; index < local_nrecs; index++) {
      pSextractor = &local_table[index];
      xplot[index2] = pSextractor->dmag;
      magdep_bin = GetMagdepBin(pMagdepLimits,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&ix,&iy,&imag);
      pMagdep = &pMagdepTable[magdep_bin];
      if (fabs(pMagdep->magdep_bin_quality) < MAGCAL_MAGDEP_QUALITY_LIMIT) {
        continue;
      }
      if (pSextractor->dmag == 0) {
        continue;
      }
      if (pSextractor->Stdmag > Y) {
        continue;
      }

#if 1
      yplot[index2] = fabs(pSextractor->magcal_local-pSextractor->Stdmag - GetMagdepBinMagcor(pMagdepLimits,pMagdepTable,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&temp_bin,&temp_magdep_rms,&temp_bin_magcor,1) )/fabs(pSextractor->dmag);
#else 
      yplot[index2] = fabs(pSextractor->magcal_local-pSextractor->Stdmag - pMagdep->magdep_bin_magcor)/fabs(pSextractor->dmag);
#endif
      totalPlotCount++;
      if (yplot[index2] <= 1.0) {
	goodPlotCount++;
      }
#if 0
      printf("orig %f %f, new %f %f result %f\n",
             pSextractor->dmag,
             fabs(pSextractor->dmag),
             pSextractor->magcal_local-pSextractor->Stdmag - pMagdep->magdep_bin_magcor,
             fabs(pSextractor->magcal_local-pSextractor->Stdmag - pMagdep->magdep_bin_magcor),
             yplot[index2]);
#endif
      if (xplot[index2] > 1.0) {
        xplot[index2] = 1.0;
      }
      if (xplot[index2] < -1.0) {
        xplot[index2] = -1.0;
      }
#if 1
      if (yplot[index2] > 5.0) {
        yplot[index2] = 5.0;
      }
      if (yplot[index2] < -1.0) {
        yplot[index2] = -1.0;
      }
#endif
      if (index2 == 0) {
        x_max = xplot[index2];
        x_min = xplot[index2];
        y_max = yplot[index2];
        y_min = yplot[index2];
      } else {
        if (x_max < xplot[index2]) {
          x_max = xplot[index2];
        }
        if (x_min > xplot[index2]) {
          x_min = xplot[index2];
        }
        if (y_max < yplot[index2]) {
          y_max = yplot[index2];
        }
        if (y_min > yplot[index2]) {
          y_min = yplot[index2];
        }
      }
      index2++;
    }
    cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);

    sprintf(title2,"%s%s good %d total %d ratio %.2f",fileroot,qualifier,goodPlotCount,totalPlotCount,(1.0*goodPlotCount)/(1.0*totalPlotCount));

    cpgmtxt ("B",2.5,0.5,0.5,"brightest original dmag");
    cpgmtxt ("L",2.0,0.5,0.5,"abs(corrected dmag)/abs(original dmag)");
    cpgmtxt ("T",1.0,0.5,0.5,title2);

    cpgpt(index2,xplot,yplot,-1);
    xplot[0] = fminf(x_min,y_min);
    yplot[0] = 1.0;
    xplot[1] = fminf(x_max,y_max);
    yplot[1] = 1.0;

    cpgline(2,xplot,yplot);

#else /* PLOT_RATIO */

    /* Plot the corrected dmag against the original dmag */
    index2 = 0;
    for (index = 0; index < local_nrecs; index++) {
      pSextractor = &local_table[index];
      xplot[index2] = pSextractor->dmag;
      magdep_bin = GetMagdepBin(pMagdepLimits,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&ix,&iy,&imag);
      pMagdep = &pMagdepTable[magdep_bin];
      yplot[index2] = pSextractor->magcal_local-pSextractor->Stdmag - pMagdep->magdep_bin_magcor;
      if (pSextractor->Stdmag > Y) {
        continue;
      }
      if (fabs(pMagdep->magdep_bin_quality) < MAGCAL_MAGDEP_QUALITY_LIMIT) {
	continue;
      }

#if 0
      printf("%f %f %d %f\n",xplot[index2],yplot[index2],magdep_bin,pMagdep->magdep_bin_magcor);
#endif
      if (xplot[index2] > 1.0) {
        xplot[index2] = 1.0;
      }
      if (xplot[index2] < -1.0) {
        xplot[index2] = -1.0;
      }
      if (yplot[index2] > 1.0) {
        yplot[index2] = 1.0;
      }
      if (yplot[index2] < -1.0) {
        yplot[index2] = -1.0;
      }
      if (index2 == 0) {
        x_max = xplot[index2];
        x_min = xplot[index2];
        y_max = yplot[index2];
        y_min = yplot[index2];
      } else {
        if (x_max < xplot[index2]) {
          x_max = xplot[index2];
        }
        if (x_min > xplot[index2]) {
          x_min = xplot[index2];
        }
        if (y_max < yplot[index2]) {
          y_max = yplot[index2];
        }
        if (y_min > yplot[index2]) {
          y_min = yplot[index2];
        }
      }
      index2++;
    }
    cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);

    cpgmtxt ("B",2.5,0.5,0.5,"brightest original dmag");
    cpgmtxt ("L",2.0,0.5,0.5,"brightest corrected constant dmag");
    sprintf(title,"%s%s",fileroot,qualifier);
    cpgmtxt ("T",1.0,0.5,0.5,title);

    cpgpt(index2,xplot,yplot,-1);
    xplot[0] = fminf(x_min,y_min);
    yplot[0] = fminf(x_min,y_min);
    xplot[1] = fminf(x_max,y_max);
    yplot[1] = fminf(x_max,y_max);

    cpgline(2,xplot,yplot);
#endif /* PLOT_RATIO */
    cpgpanl(2,3);
    /* Plot the interpolated corrected dmag against the original dmag */
    index2 = 0;
    for (index = 0; index < local_nrecs; index++) {
      pSextractor = &local_table[index];
      xplot[index2] = pSextractor->dmag;
      magdep_bin = GetMagdepBin(pMagdepLimits,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&ix,&iy,&imag);
      pMagdep = &pMagdepTable[magdep_bin];
      yplot[index2] = pSextractor->magcal_local-pSextractor->Stdmag - GetMagdepBinMagcor(pMagdepLimits,pMagdepTable,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,pSextractor->Stdmag,&temp_bin,&temp_magdep_rms,&temp_bin_magcor,1);
      if (pSextractor->Stdmag > Y) {
        continue;
      }
      if (fabs(pMagdep->magdep_bin_quality) < MAGCAL_MAGDEP_QUALITY_LIMIT) {
	continue;
      }

#if 0
      printf("%f %f %d %f\n",xplot[index2],yplot[index2],magdep_bin,pMagdep->magdep_bin_magcor);
#endif
      if (xplot[index2] > 1.0) {
        xplot[index2] = 1.0;
      }
      if (xplot[index2] < -1.0) {
        xplot[index2] = -1.0;
      }
      if (yplot[index2] > 1.0) {
        yplot[index2] = 1.0;
      }
      if (yplot[index2] < -1.0) {
        yplot[index2] = -1.0;
      }
      if (index2 == 0) {
        x_max = xplot[index2];
        x_min = xplot[index2];
        y_max = yplot[index2];
        y_min = yplot[index2];
      } else {
        if (x_max < xplot[index2]) {
          x_max = xplot[index2];
        }
        if (x_min > xplot[index2]) {
          x_min = xplot[index2];
        }
        if (y_max < yplot[index2]) {
          y_max = yplot[index2];
        }
        if (y_min > yplot[index2]) {
          y_min = yplot[index2];
        }
      }
      index2++;
    }
    cpgenv(x_min-0.1,x_max+0.1,y_min-0.1,y_max+0.1,0,0);

    cpgmtxt ("B",2.5,0.5,0.5,"brightest original dmag");
    cpgmtxt ("L",2.0,0.5,0.5,"brightest corrected interpolated dmag");
    sprintf(title,"%s%s",fileroot,qualifier);
    cpgmtxt ("T",1.0,0.5,0.5,title);

    cpgpt(index2,xplot,yplot,-1);
    xplot[0] = fminf(x_min,y_min);
    yplot[0] = fminf(x_min,y_min);
    xplot[1] = fminf(x_max,y_max);
    yplot[1] = fminf(x_max,y_max);

    cpgline(2,xplot,yplot);



    if (xplot != NULL) {
      free(xplot);
      xplot = NULL;
    }
    if (yplot != NULL) {
      free(yplot);
      yplot = NULL;
    }

    cpgiden();
    cpgend();
  }



  if (calfile_handle == NULL) {
    fclose(calfile_handle);
  }


  if (outHandle != NULL) {
    fclose(outHandle);
  }
  if (local_handle != NULL) {
    Close(local_handle);
    local_handle = NULL;
  }
  if (local_header != NULL) {
    table_hdrfree(local_header);
    local_header = NULL;
  }
  if (local_table != NULL) {
    Free(local_table);
    local_table = NULL;
  }
  if (vector1 != NULL) {
    free(vector1);
    vector1 = NULL;
  }
  if (vector2 != NULL) {
    free(vector2);
    vector2 = NULL;
  }
  if (vector3 != NULL) {
    free(vector3);
    vector3 = NULL;
  }

  if (Ybright != NULL) {
    free(Ybright);
    Ybright = NULL;
  }
  if (Ybright2 != NULL) {
    free(Ybright2);
    Ybright2 = NULL;
  }
  if (dybright != NULL) {
    free(dybright);
    dybright = NULL;
  }
  if (dybright2 != NULL) {
    free(dybright2);
    dybright2 = NULL;
  }
  if (Yfaint != NULL) {
    free(Yfaint);
    Yfaint = NULL;
  }
  if (nmag != NULL) {
    free(nmag);
    nmag = NULL;
  }
  if (xrange != NULL) {
    free(xrange);
    xrange = NULL;
  }
  if (yrange != NULL) {
    free(yrange);
    yrange = NULL;
  }
  if (Zrms1 != NULL) {
    free(Zrms1);
    Zrms1 = NULL;
  }
  if (Z0 != NULL) {
    free(Z0);
    Z0 = NULL;
  }
  if (Z1 != NULL) {
    free(Z1);
    Z1 = NULL;
  }
  if (Z2 != NULL) {
    free(Z2);
    Z2 = NULL;
  }
  if (M1 != NULL) {
    free(M1);
    M1 = NULL;
  }
  if (X1 != NULL) {
    free(X1);
    X1 = NULL;
  }
  if (intvec1 != NULL) {
    free(intvec1);
    intvec1 = NULL;
  }
  if (intvec2 != NULL) {
    free(intvec2);
    intvec2 = NULL;
  }
  if (intvec3 != NULL) {
    free(intvec3);
    intvec3 = NULL;
  }
  if (intvec4 != NULL) {
    free(intvec4);
    intvec4 = NULL;
  }
  if (intvec5 != NULL) {
    free(intvec5);
    intvec5 = NULL;
  }
  if (Zmag != NULL) {
    free(Zmag);
    Zmag = NULL;
  }
  if (M2 != NULL) {
    free(M2);
    M2 = NULL;
  }

  if (pMagdepTable != NULL) {
    free(pMagdepTable);
    pMagdepTable = NULL;
  }


  time(&curTime);
  curTime -= startTime;
  printf("magdepcalibrate seconds %d length %d goodStar %d %.2f rejected %d Stdmagrejected %d bins %d goodbins %d %f for %s%s\n",
         curTime,
         local_nrecs,
         goodStarCount,
         (1.0*goodStarCount)/local_nrecs,
	 rejectedCount,
         StdmagRejected,
         pMagdepLimits->nrecs,
         goodBinCount,
         (1.0*goodBinCount)/pMagdepLimits->nrecs,
         fileroot,qualifier);

  printf("magdepcalibrate orig median %.2f rms %.2f const median %.2f rms %.2f, final median %.2f rms %.2f for %s%s\n",
	 orig_median,orig_rms,
	 const_median,const_rms,
	 final_median,final_rms,
	 fileroot,qualifier);


  return(0);
}


