// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* contour_plot.c
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64     -I /dasch/install/include -I /dasch/install/pgplot   contour_plot.c pipelineutils.a -L /dasch/install/lib  -lwcs -lm  -L /dasch/install/pgplot -I/usr/include/mysql  -lcpgplot -lpgplot  -L/usr/X11R6/lib  -lX11 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 /usr/lib${lib64}/libg2c.so.0  -ltable -lutil -L/usr/lib${lib64}/mysql  -l mysqlclient   -o contour_plot 
 *
 * This is a port of contour_plot_clip_smooth_correct2.f written by Silas Laycock
 * from FORTRAN to C.
 *
 * Make 2D contour/greyscale plot (mean[Z]) from irregularly spaced X,Y,Z data points
 * 1. Bin the data into a grid (bin_value=average of member values)
 * 2. vizualize the array
 * This version contour_plot_clip.f uses an iterative sigma-clipping routine to 
 * derive mean, median & sigma. Default settings are 3 iterations, removing 3sigma outliers
 * Parameter "clip" specifies number of clipping iterations. If non integer (ie. M.N), 
 * the part following the decimal point is treated as a stopping condition, and 
 * the clipping-code will keep iterating until fractional improvement in sigma is less than 0.N
 * In contour_plot_clip_smooth.f, a weighting function is used to sample data from an
 * arbitrarily sized "boxcar" onto a finely spaced grid. -thus smoothing out fluctuations
 * and interpolating-over sparse areas. The smoothing scale is set by parameter "scale"
 * which is specified in the same units as X & Y.
 * 16Jul2007: Correct the input Z values using the smooth clipped-median error model.
 *
 *
 * Jan 22, 2008 Edward J. Los - Ported from contour_plot_clip_smooth_correct2.f
 *                              Corrected bad cell replacement algorithm at plate edges
 *                              Used center of grid square rather than lower left corner
 *                               for the boxcar algorithm
 *                              Corrected ERROR: checksum = 0.00000 for datapoint errors
 *                              Added a grid file to calibrate all sextractor objects
 *                              Added a grid bin column for cross-checking
 *                              Make the generation of plots dependent on an environment
 *                               variable
 * Feb 22, 2008 Edward J. Los - Add a verbose flag
 * Feb 25, 2008 Edward J. Los - Use 1.5 bins for a smoothing scale
 *                              Perform the local calibration twice:  first for every single bin, then
 *                              for blocks of 3x3 bins.  The smoothing radius will be 1.5 bins for the
 *                              first calculation and 4.5 bins for the second calculation.
 *                              If any single bin has at least 20 stars, the result will be accepted.
 *                              Otherwise, the result will be replaced with the result for the 3x3 bin 
 *                              block.  If the 3x3 bin block still has less than 20 stars, the standard
 *                              deviation will be set to 99 as a flag that local calibration on these
 *                              stars should be skipped. When local calibration is skipped, the final result
 *                              will use the lowess magnitude and RMS and will be appropriately flagged.
 * Feb 23, 2011 Edward J. Los   Increase the minimum star count per bin from 20 to 30
 *                              Restrict adjacent bins to the same spatial_bin
 *                              Expand to a 5x5 surrounding area if a 3x3 does not have sufficient stars                             
 * Feb 26, 2011 Edward J. Los   Restructure arrays for recursive processing
 *                              Perform calculation twice: the second time for bright stars only within two magnitudes of the brightest good star.
 * Mar 11, 2011 Edward J. Los   Move functionality previously in the run_local_calibration.csh script to this routine.
 *                              Convert all float to double
 *                              Filter out  too bright, too dim, and high rms stars
 * Mar 14, 2011 Edward J. Los   Get spatial_bin and local_bin from the input file
 * Jul 25, 2017 Edward J. Los   Remove the nx and ny parameters and make the number of smoothing bins plate size dependent.
 */

#include <math.h>
#include "table.h"
#define DO_PLOT 1
#define MIN_STARS_PER_BIN 30
#define BRIGHT_MAGNITUDE_RANGE 2  /* Add this to the brightest magcal_iso in the bin to get the dimmest star to consider when calculating the "bright star" parameters */
#ifdef DO_PLOT
#include "cpgplot.h"
#include "pipelineutils.h"
#endif /* DO_PLOT */
#define MAX_BUFFER 512
#define lim2 X_DMAGBINS_NORMAL
#define lim3 10
#define BOUNDS_CHECK 1
/* These are global because they take up too much stack */

static int printFlag = 0;


typedef struct _inputentry {
  char REF[MAX_REF+1];
  double X_IMAGE; /* Formerly X */
  double Y_IMAGE; /* Formerly Y */
  double ra;
  double dec;
  double refmag;
  double color;
  double MAG_ISO;
  double extinction;
  double magcal_iso_rms;
  double dmag_iso; /* Formerly Z */
  double magcal_iso;
  double limiting_mag;
  double max_bright_mag;
  double magcor_local;  /* Formerly ZOUT */
  double magcal_local_error; /* Formerly ERROUT */
  double bright_magcor_local;  /* magcal_local for brightest two magnitudes */
  double bright_magcor_error; /* magcal_local_error for brightest two magnitudes */
  double sortval; /* sorting value.  first dmag_iso; later magcal_local; */
  int NUMBER;
  int BFLAGS;
  int npoints_local; /* Formerly NPOUT */
  int bright_npoints_local; /* npoints local for brightest two magnitudes*/
  int local_bin; /* Formerly LOCALBIN */
  int spatial_bin;
} INPUTENTRY,*PINPUTENTRY;

typedef struct _common {
  int nlines; /* Size of the input table */
  int np; /* Number of valid points in the input table */
  int spatial_bin2[lim2][lim2];
  int nvalues; /* Number of valid points in the select table */
  PINPUTENTRY pInputTable;
  PINPUTENTRY pSelectTable1;
  PINPUTENTRY pSelectTable2;
} COMMON,*PCOMMON;

/* We need an array for each possible combination, unless we get the Fortran vs Linux array indexing correct */
#define PLOT_TYPE_N 1
#define PLOT_TYPE_SL 2
#define PLOT_TYPE_LS 3
float PLOTOUT_N[X_DMAGBINS_NORMAL][Y_DMAGBINS_NORMAL];            /* PLOT_TYPE_N */
float PLOTOUT_SL[SHORT_DMAGBINS_SMALL][LONG_DMAGBINS_SMALL];     /* PLOT_TYPE_SL */
float PLOTOUT_LS[LONG_DMAGBINS_SMALL][SHORT_DMAGBINS_SMALL];     /* PLOT_TYPE_LS */
double MED1;
double MED2;
double MED3;
double MED4;

double MED[lim2][lim2];
double POINT[lim2][lim2];
double SIGMA[lim2][lim2];
double bright_magcor_local[lim2][lim2];
double bright_magcor_error[lim2][lim2];
double bright_npoints_local[lim2][lim2];


double SIGMA1;
double SIGMA2;
double SIGMA3;
double SIGMA4;
float TR[6];
float CON[lim3];	
double MAP[lim2][lim2];
double doublesqr(double input) {
  return(input*input);
}
double cvtreal(int input) {
  return(1.0 * input);
}
int INT(double input) {
  int retval = input;
  return(retval);
}

int sortval_compare(const void *first, const void *second) 
{
  double sortvalFirst = ((PINPUTENTRY)first)->sortval;
  double sortvalSecond = ((PINPUTENTRY)second)->sortval;
  if (sortvalFirst > sortvalSecond) {
    return(1);
  } else if (sortvalFirst < sortvalSecond) {
    return(-1);
  } else {
    return(0);
  }

}
/*
 *  sigmaclip
 *  compute the median and RMS of a set of numbers
 *  using an iterative sigma clipping, arguments are
 *  the Nsigma and Niterations for clipping
 *  If Niter is non integer, then the decimal 
 *  part is treated as as a stopping condition
 *  and the code will iterate until the fractional 
 *  improvement in sigma is less that 0.N
 *  The routine returns the clipped statistical values
 *  and the array X modified to contain only the "good" points.
 */
void sigmaclip (PCOMMON pCommon,PINPUTENTRY pSelectTable1,PINPUTENTRY pSelectTable2,int npoints,int *ngood,double nsigma,
                double iter,double *mean,double *median,double *sigma) 
{
  double sigma2;
  double mean2;
  double median2;
  int ngood2;
  double sigma1;
  double mean1;
  double median1;
  double sum;
  double sumsq;
  double  min;
  double  max;
  double  stop = 0.0;
  double *temp1 = NULL;
  double *temp2 = NULL;
  int niter;
  int i;
  int j;
  PINPUTENTRY pSelectInput;
  PINPUTENTRY pSelectInput2;
  

  if (npoints <= 1) {
    *mean=0.0;
    *sigma=0.0;
    *median=0.0;
    *ngood = 0;
    return;

  }


 
  /* initialize values of certain variables */
  ngood2=npoints;
  if ((iter-cvtreal(INT(iter))) > 0.0) {
    niter = 1000;
    stop = iter-cvtreal(INT(iter));
  } else {
    niter = INT(iter);
  }

  /* First pass to get raw mean */
  sum = 0.0;
  for (i = 0; i < npoints; i++) {
    pSelectInput = &pSelectTable1[i];
#ifdef BOUNDS_CHECK
    if ((i >= pCommon->nlines) ||
        (i < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",i,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */


    sum = sum + pSelectInput->sortval;

  }
  /*        npoints = i-1; */
  mean1 = sum/cvtreal(npoints);
	
  /* Get the raw Median */
#ifdef BOUNDS_CHECK
    if ((npoints > pCommon->nlines) ||
        (npoints < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",npoints,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */


  qsort((void*)pSelectTable1,npoints,sizeof(INPUTENTRY),sortval_compare);

/* Sort routine based on the sortval  */
  if (npoints >= 2) {
    pSelectInput = &pSelectTable1[(npoints/2)-1];
#ifdef BOUNDS_CHECK
    if ((((npoints/2)-1) >= pCommon->nlines) ||
        (((npoints/2)-1) < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",(npoints/2)-1,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */

    median1 = pSelectInput->sortval;
  } else {
    if (npoints > 0) {
      pSelectInput = &pSelectTable1[0];
      median1 = pSelectInput->sortval;
    } else {
      median1 = 0;
    }
  }
      

  /* Second pass to get raw RMS or standard eviation */
  sumsq = 0.0;
  for (i = 0; i < npoints; i++) {
    pSelectInput = &pSelectTable1[i];
#ifdef BOUNDS_CHECK
    if ((i >= pCommon->nlines) ||
        (i < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",i,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */

    sumsq = sumsq + doublesqr(pSelectInput->sortval-mean1);    
  }
  sigma1 = sqrt(sumsq/cvtreal(npoints-1));

  /*	print*,'iteration 0',' ngood2',ngood2,' mean',mean1,' sigma',sigma1 */
	
  /* Now iterate the requested number of times, to obtain a cleaned up "clipped" mean and sigma */
  sigma2 = sigma1;
  mean2 = mean1;
  median2 = median1;

  for (i = 0; i < niter; i++) {



    mean1 = mean2;
    sigma1 = sigma2;
    min = mean2 - cvtreal(nsigma)*(sigma2);
    max = mean2 + cvtreal(nsigma)*(sigma2);
    ngood2 = 0 ;
    sum = 0.0;
    sumsq = 0.0;
    for (j = 0; j < npoints; j++) {
      pSelectInput = &pSelectTable1[j];
      pSelectInput2 = &pSelectTable2[ngood2];
#ifdef BOUNDS_CHECK
    if ((j >= pCommon->nlines) ||
        (j < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",j,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
    if ((ngood2 >= pCommon->nlines) ||
        (ngood2 < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",ngood2,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */


      if ((pSelectInput->sortval > min) && (pSelectInput->sortval < max)) {
        sum = sum + pSelectInput->sortval;
        ngood2 = ngood2 + 1;
#ifdef BOUNDS_CHECK
        if (((ngood2-1) < 0) || ((ngood2-1) >= npoints)) {
          printf("ERROR: %s %d Array bounds for temp1 exceeded %d %d\n",__FUNCTION__,__LINE__,ngood2,npoints);
          exit(-1);
        }
#endif /* BOUNDS_CHECK */
        memcpy(pSelectInput2,pSelectInput,sizeof(INPUTENTRY));
      }				
    }
    mean2 = sum/cvtreal(ngood2);
    qsort((void*)pSelectTable2,ngood2,sizeof(INPUTENTRY),sortval_compare);
    if (ngood2 >= 2) {
      pSelectInput2 = &pSelectTable2[(ngood2/2)-1];
#ifdef BOUNDS_CHECK
    if ((((ngood2/2)-1) >= pCommon->nlines) ||
        (((ngood2/2)-1) < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",((ngood2/2)-1),pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */


      median2 = pSelectInput2->sortval;
    } else if (ngood2 == 1) {
      pSelectInput2 = &pSelectTable2[0];
      median2 = pSelectInput2->sortval;      
    } else {
      median2 = 0;
    }
    for( j = 0; j < ngood2; j++) {
      pSelectInput2 = &pSelectTable2[j];
#ifdef BOUNDS_CHECK
    if ((j >= pCommon->nlines) ||
        (j < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",j,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */


      sumsq = sumsq + doublesqr(pSelectInput2->sortval-mean2);
    }
    sigma2 = sqrt(sumsq/cvtreal(ngood2-1));
		
    /* c		if ( (mean1-mean2).lt.stop ) go to 100  */
    if ( (sigma1-sigma2) < stop ) {
      break;
    }
  }
  /* Make sure that weird numbers dont get spat out. */
  if (ngood2 == 0) {
    mean2=0.0;
    sigma2=0.0;
    median2=0.0;
  }
  *mean=mean2;
  *sigma=sigma2;
  *median=median2;
  *ngood = ngood2;

  /*
   * print helpful de-bugging information
   */
#if 0
  printf("npoints=%d, ngoodpoints = %d,niter=%d\n",npoints,ngood2,niter);
  printf("mean1=%f, mean=%f\n",mean1,mean2);
  printf("median1=%f, median=%f\n",median1,median2);
  printf("sigma1=%f, sigma=%f\n",sigma1,sigma2);
  for (i=0; i < ngood2; i++) {
    printf("%f ",x[i]);
  }
  printf("\n");
#endif
  return;
}
void CalcParameters(PCOMMON pCommon,int i,int j, double x1, double y1,double scale,double clip,int iteration,int brightFlag,double *POINT,double *MED,double *SIGMA)
{
  int nvalues;
  int k;
  double dist;
  PINPUTENTRY pInputEntry;
  PINPUTENTRY pSelectEntry;
  PINPUTENTRY pSelectEntry2;
  double ave;
  double med;
  double rms;
  int goodvals;
  double expandedscale;
  double min_magcal_iso = 99.0;
  switch(iteration) {
  case -1:
    /* Last time for plotting */
    expandedscale = scale;
    break;
  case 0:
    expandedscale = scale;
    break;
  case 1:
    expandedscale = 3.0 * scale;
    break;
  case 2:
    expandedscale = 5.0 * scale;
    break;
  default:
    printf("ERROR: Illegal iteration %d\n",iteration);
    exit(-1);

    return;
  }


  pCommon->nvalues=0;
  for ( k=0; k< pCommon->np; k++) {
    pInputEntry = &pCommon->pInputTable[k];
    pSelectEntry = &pCommon->pSelectTable1[pCommon->nvalues];
#ifdef BOUNDS_CHECK
    if ((k >= pCommon->nlines) ||
        (k < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",k,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */


    dist= sqrt(doublesqr(x1-pInputEntry->X_IMAGE) + doublesqr(y1-pInputEntry->Y_IMAGE));

    if ((dist < expandedscale) && 
        ((iteration < 0) || (pInputEntry->spatial_bin == pCommon->spatial_bin2[i][j]))) {
#if 0
      printf("i=%d j=%d k=%d dist=%f X Y Z = %f %f %f spatial_bin %d %d \n",
             i+1,
             j+1,
             k+1,
             dist,
             pInputEntry->X_IMAGE,
             pInputEntry->Y_IMAGE,
             pInputEntry->sortval,
             pInputEntry->spatial_bin,
             pCommon->spatial_bin2[i][j]);
#endif
      pCommon->nvalues=pCommon->nvalues+1;
#ifdef BOUNDS_CHECK
      if (((pCommon->nvalues-1) < 0) || ((pCommon->nvalues-1) >= pCommon->nlines)) {
        printf("ERROR: %s %d Array bounds for TMP1 exceeded %d %d\n",
               __FUNCTION__,__LINE__,pCommon->nvalues,pCommon->nlines);
        exit(-1);
      }
#endif /* BOUNDS_CHECK */
      memcpy(pSelectEntry,pInputEntry,sizeof(INPUTENTRY));
      if (pSelectEntry->magcal_iso < min_magcal_iso) {
        min_magcal_iso = pSelectEntry->magcal_iso;
      }
    }
  }
  if (printFlag) {
    printf("Have (1) %d nvalues brightFlag %d\n",pCommon->nvalues,brightFlag);
    for (k = 0; k < pCommon->nvalues; k++) {
      pSelectEntry = &pCommon->pSelectTable1[k];
      printf("%.6f\n",pSelectEntry->magcal_iso);
    }
  }


  /* When the bright flag is set, consider only magnitudes from min_magcal_iso to min_magcal_iso+BRIGHT_MAGNITUDE_RANGE */
  if ((brightFlag != 0) && 
      (pCommon->nvalues > 0)) {
    k = 0;
    while (k < pCommon->nvalues) {
      pSelectEntry = &pCommon->pSelectTable1[k];
      if (pSelectEntry->magcal_iso > (min_magcal_iso+BRIGHT_MAGNITUDE_RANGE)) {
        pSelectEntry2 = &pCommon->pSelectTable1[pCommon->nvalues-1];
        memcpy(pSelectEntry,pSelectEntry2,sizeof(INPUTENTRY));
        pCommon->nvalues--;
      } else {
        k++;
      }
    }

    if (printFlag) {
      printf("Have (2) %d nvalues brightFlag %d\n",pCommon->nvalues,brightFlag);
      for (k = 0; k < pCommon->nvalues; k++) {
        pSelectEntry = &pCommon->pSelectTable1[k];
        printf("%.6f\n",pSelectEntry->magcal_iso);
      }
    }

  }

  ave=0.0;
  med=0.0;
  rms=0.0;
  sigmaclip (pCommon,pCommon->pSelectTable1,pCommon->pSelectTable2,pCommon->nvalues,&goodvals,3.0,clip,&ave,&med,&rms)  ;
  if (iteration == -1) {
    /* Single interation mode for plotting */
    *POINT=cvtreal(goodvals);
    *MED=med;
    *SIGMA=rms;
    return;
  }

  if (goodvals >= MIN_STARS_PER_BIN) {
    *POINT=cvtreal(goodvals);
    *MED=med;
    *SIGMA=rms;
  } else {
    *POINT=cvtreal(goodvals);
    *MED=med;
    /* No success, increase size of the bin */
    if (iteration < 2) {
      CalcParameters(pCommon,i,j,x1,y1,scale,clip,iteration+1,brightFlag,POINT,MED,SIGMA);
    } else {
      *SIGMA=99.0;      
    }

  }

#if 0
  printf("%d %d %d %d %f %f %f %d\n",
         i+1,j+1,pCommon->nvalues,goodvals,med,rms,ave,pCommon->spatial_bin2[i][j]);
#endif
}







int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  char xtext[MAX_BUFFER];
  char ytext[MAX_BUFFER];
  char plotdev[MAX_BUFFER];
  char outfile[MAX_BUFFER];
  char gridfile[MAX_BUFFER];
  FILE *outHandle = NULL;
  FILE *gridHandle = NULL;
  int plot_type = -1;
  int nx;
  int ny;
  int mx;
  int my;
  int ncon;
  int goodvals;
  int goodvaltmp;
  double goodval2;
  double x1;
  double x2;
  double y1;
  double y2;
  double dx;
  double dy;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double temp;
  double medmin;
  double medmax;
  double sigmax;
  double sigmin;
  double npmin;
  double npmax;
  double dz;
  double ave;
  double med;
  double rms;
  double avebright;
  double medbright;
  double rmsbright;
  double avedim;
  double meddim;
  double rmsdim;
  double clip;
  double scale;
  double smoothDegrees;
  double arcsecPerPixel;
  double dist;	
  double mapmin;
  double mapmax;
  int i;
  int j;
  int k;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  char *tabPtr;
  int lineLen;
  char *tabLoc;
  int iteration;
  int brightFlag;
  int tmp_spatial_bin;

  COMMON globalcommon;
  PCOMMON pCommon = &globalcommon;
  PINPUTENTRY pInputEntry;
  PINPUTENTRY pSelectEntry;

  char *plotValue;
  int plotFlag = 0;
  int verbose = 0;
  double edgeDist;
  int printIndex;

  File input_handle = NULL;
  char input_name[MAX_BUFFER];
  TableHead input_header = NULL;
  PINPUTENTRY input_table = NULL;
  size_t input_nrecs = 0;
  int input_index;
  PINPUTENTRY pInput;
  double magcal_local;
  double dmag_local;
  int cal_local = 1;
  int bright_reject_count = 0;
  int max_rms_count = 0;

  memset(pCommon,0,sizeof(COMMON));


  plotValue = getenv("DASCH_PLOT");
  if (strstr(plotValue,"YES") != NULL) {
    plotFlag = 1;
  }

  input_name[0] = 0;
  outfile[0] = 0;
  gridfile[0] = 0;
  plotdev[0] = 0;
  xtext[0] = 0;
  ytext[0] = 0;
  nx = -1;
  ny = -1;
  mx = -1;
  my = -1;
  ncon = -1;
  clip = 0.0;
  smoothDegrees = 0.0;
  arcsecPerPixel = 0.0;
  scale = 0.0;


  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      fprintf(stderr,"ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'i': /* Input file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(input_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;
        case 'g': /* grid file name */
        case 'G':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(gridfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;
        case 'p': /* plot file name */
        case 'P':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(plotdev,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;
        case 'l': /* width or height label */
        case 'L':
          cmdchar = cmdchar = *++argstr;
          switch(cmdchar) {
          case 'x': /* width label */
          case 'X':
            argc--;
            if (argc < 1) {
              fprintf(stderr,"ERROR: Insufficient arguments for -l%c\n",cmdchar);
              errorFlag = 1;
            } else {
              strncpy(xtext,*++argv,MAX_BUFFER-2);
              if (strlen(*argv) >= MAX_BUFFER-2) {
                fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
              }
            }
            break;
          case 'y': /* height label */
          case 'Y':
            argc--;
            if (argc < 1) {
              fprintf(stderr,"ERROR: Insufficient arguments for -l%c\n",cmdchar);
              errorFlag = 1;
            } else {
              strncpy(ytext,*++argv,MAX_BUFFER-2);
              if (strlen(*argv) >= MAX_BUFFER-2) {
                fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
              }
            }
            break;
          default:
            fprintf(stderr,"ERROR:  command -l%c-",cmdchar);
            errorFlag = 1;

          }
          break;
          

        case 'm': /* width or height */
        case 'M':
          cmdchar = cmdchar = *++argstr;
          switch(cmdchar) {
          case 'x': /* number of x bins */
          case 'X':
            argc--;
            if (argc < 1) {
              fprintf(stderr,"ERROR: Insufficient arguments for -m%c\n",cmdchar);
              errorFlag = 1;
            } else {
              nvals = sscanf(*++argv,"%d",&mx);
              if (nvals != 1) {
                fprintf(stderr,"ERROR: Can not decode width\n");
                errorFlag = 1;
              }
            }
            break;
          case 'y': /* height label */
          case 'Y':
            argc--;
            if (argc < 1) {
              fprintf(stderr,"ERROR: Insufficient arguments for -m%c\n",cmdchar);
              errorFlag = 1;
            } else {
              nvals = sscanf(*++argv,"%d",&my);
              if (nvals != 1) {
                fprintf(stderr,"ERROR: Can not decode height\n");
                errorFlag = 1;
              }
            }
            break;
          default:
            fprintf(stderr,"ERROR:  command -m%c-",cmdchar);
            errorFlag = 1;

          }
          break;


        case 'c': /* contour lines */
        case 'C':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&ncon);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Can not decode number of contour lines\n");
              errorFlag = 1;
            }
            if (ncon>lim3) {
              fprintf(stderr,"ERROR: ncon %d greater than lim3 %d\n",ncon,lim3);
              errorFlag = 1;
            }
          }
          break;

        case 'r': /* clipping ratio */
        case 'R':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%f",&clip);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Can not clipping ratio\n");
              errorFlag = 1;
            }
          }
          break;

        case 's': /* smoothing degrees */
        case 'S':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%f",&smoothDegrees);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Can not decode number of smoothing degrees\n");
              errorFlag = 1;
            }
          }
          break;


        case 'a': /* arcsec/pixel */
        case 'A':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%f",&arcsecPerPixel);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Can not decode number of arcsec per pixel\n");
              errorFlag = 1;
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        default:
          fprintf(stderr,"ERROR:  command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */
  if (input_name[0] == 0) {
    fprintf(stderr,"ERROR: No input filename was specified\n");
    errorFlag = 1;
  }
  if (outfile[0] == 0) {
    fprintf(stderr,"ERROR: No output filename was specified\n");
    errorFlag = 1;
  }
  if (gridfile[0] == 0) {
    fprintf(stderr,"ERROR: No grid filename was specified\n");
    errorFlag = 1;
  }
  if (plotdev[0] == 0) {
    fprintf(stderr,"ERROR: No plot filename was specified\n");
    errorFlag = 1;
  }
  if (xtext[0] == 0) {
    fprintf(stderr,"ERROR: No x label was specified\n");
    errorFlag = 1;
  }
  if (ytext[0] == 0) {
    fprintf(stderr,"ERROR: No y label was specified\n");
    errorFlag = 1;
  }


  if (mx == -1) {
    fprintf(stderr,"ERROR: No x width was specified\n");
    errorFlag = 1;
  }
  if (my == -1) {
    fprintf(stderr,"ERROR: No y height was specified\n");
    errorFlag = 1;
  }
  nx = XDmagBins(mx,my);
  ny = YDmagBins(mx,my);
  if ((nx == X_DMAGBINS_NORMAL) && (ny == Y_DMAGBINS_NORMAL)) {
    plot_type = PLOT_TYPE_N;
  } else if ((nx == SHORT_DMAGBINS_SMALL) && (ny == LONG_DMAGBINS_SMALL)) {
    plot_type = PLOT_TYPE_SL;
  } else if ((nx == LONG_DMAGBINS_SMALL) && (ny == SHORT_DMAGBINS_SMALL)) {
    plot_type = PLOT_TYPE_LS;
  } else {
    printf("ERROR: Unrecognized plot_type for nx %d and ny %d for %s\n",nx,ny,input_name);
    errorFlag = 1;
  }

  if (ncon == -1) {
    fprintf(stderr,"ERROR: No number of contour lines was specified\n");
    errorFlag = 1;
  }

  if (clip == 0.0) {
    fprintf(stderr,"ERROR: No clip ratio was specified\n");
    errorFlag = 1;
  }
  if (smoothDegrees == 0.0) {
    fprintf(stderr,"ERROR: No smoothing degrees was specified\n");
    errorFlag = 1;
  }

  if (arcsecPerPixel == 0.0) {
    fprintf(stderr,"ERROR: No arcsec per pixel plate scale was specified\n");
    errorFlag = 1;
  }


  input_handle = Open(input_name,"r");
  if (input_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the input file %s\n",input_name);
  } else {
    if (verbose) {
      fprintf(stderr,"Found input file %s\n",input_name);
    }
  }



  outHandle = fopen(outfile,"wt");
  if (outHandle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to open outfile %s\n",outfile);
  }
  fprintf(outHandle,"REF\tra\tdec\trefmag\tcolor\tMAG_ISO\tNUMBER\tBFLAGS\tX_IMAGE\tY_IMAGE\textinction\tmagcal_iso\tdmag_iso\tmagcal_iso_rms\tspatial_bin\tlimiting_mag\tmagcor_local\tmagcal_local_error\tnpoints_local\tlocal_bin\tmagcal_local\tdmag_local\tcal_local\n");
  fprintf(outHandle,"---\t--\t---\t------\t-----\t-------\t------\t------\t-------\t-------\t----------\t----------\t--------\t--------------\t-----------\t------------\t------------\t------------------\t-------------\t---------\t------------\t----------\t---------\n");

  gridHandle = fopen(gridfile,"wt");
  if (gridHandle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to open gridfile %s\n",gridfile);
  }



  if (errorFlag) {
    fprintf(stderr,"Usage: contour_plot -i <input file> \n");
    fprintf(stderr,"                    -lx <width label> \n");
    fprintf(stderr,"                    -ly <height label> \n");
    fprintf(stderr,"                    -mx <width in pixels> \n");
    fprintf(stderr,"                    -my <height in pixels> \n");
    fprintf(stderr,"                    -c  <Number of contour lines>\n");
    fprintf(stderr,"                    -r  <clipping ratio>\n");
    fprintf(stderr,"                    -a  <plate scale arcsec/pixel>\n");
    fprintf(stderr,"                    -s  <smoothing degrees>\n");
    fprintf(stderr,"                    -o  <output file>\n");
    fprintf(stderr,"                    -g  <output grid file>\n");
    fprintf(stderr,"                    -v  verbose mode\n");
    return(-1);
  }

  /* Convert the smoothing degrees to a radius */

  /* Read in the input file */

  input_header = table_header(input_handle,TABLE_PARSE);
  if (input_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",input_name);
    return(-1);
  }

  input_table = table_loadva(input_handle,
                             &input_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(INPUTENTRY),
                             &input_nrecs,
                             TblInt,"NUMBER",TblOff(PINPUTENTRY,NUMBER),
                             TblInt,"BFLAGS",TblOff(PINPUTENTRY,BFLAGS),
                             TblInt,"spatial_bin",TblOff(PINPUTENTRY,spatial_bin),
                             TblInt,"local_bin",TblOff(PINPUTENTRY,local_bin),
                             TblDbl,"ra",TblOff(PINPUTENTRY,ra),
                             TblDbl,"dec",TblOff(PINPUTENTRY,dec),
                             TblDbl,"refmag",TblOff(PINPUTENTRY,refmag),
                             TblDbl,"color",TblOff(PINPUTENTRY,color),
                             TblDbl,"extinction",TblOff(PINPUTENTRY,extinction),
                             TblDbl,"MAG_ISO",TblOff(PINPUTENTRY,MAG_ISO),
                             TblDbl,"magcal_iso_rms",TblOff(PINPUTENTRY,magcal_iso_rms),
                             TblDbl,"limiting_mag",TblOff(PINPUTENTRY,limiting_mag),
                             TblDbl,"max_bright_mag",TblOff(PINPUTENTRY,max_bright_mag),
                             TblDbl,"magcal_iso",TblOff(PINPUTENTRY,magcal_iso),
                             TblDbl,"dmag_iso",TblOff(PINPUTENTRY,dmag_iso),
                             TblDbl,"X_IMAGE",TblOff(PINPUTENTRY,X_IMAGE),
                             TblDbl,"Y_IMAGE",TblOff(PINPUTENTRY,Y_IMAGE),
                             TblBuf,"REF"    ,TblOff(PINPUTENTRY,REF),MAX_REF,
                             0,"end",0);
  if (input_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",input_name);
    return(-1);
  }
  
  pCommon->nlines = input_nrecs;
  pCommon->pInputTable = input_table;
  pCommon->np = input_nrecs;
  if (pCommon->nlines == 0) {
    fprintf(stderr,"ERROR: no lines in the input file %s\n",input_name);
    return(-1);
  }
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pSelectEntry = &pCommon->pInputTable[input_index];
#if 0
    if (pSelectEntry->NUMBER == 10942) {
      printf("At NUMBER %d\n",pSelectEntry->NUMBER);
    }
#endif
    tmp_spatial_bin = CalculateBin(mx,my,pSelectEntry->X_IMAGE,pSelectEntry->Y_IMAGE,&edgeDist);
    if (tmp_spatial_bin != pSelectEntry->spatial_bin) {
      printf("ERROR: spatial_bin mismatch %d %d for NUMBER %d in %s\n",pSelectEntry->spatial_bin,tmp_spatial_bin,pSelectEntry->NUMBER,input_name);
    }
    pSelectEntry->sortval = pSelectEntry->dmag_iso;
    
  }


  /* At this point, filter the input table, removing all of the stars within 1/2 magnitude of the limiting magnitude */
  input_index = 0;
  while (input_index < pCommon->np) {
    pSelectEntry = &pCommon->pInputTable[input_index];
    pInputEntry = &pCommon->pInputTable[pCommon->np-1];
#if 0
    if (pSelectEntry->NUMBER == 10942) {
      printf("At %d\n",pSelectEntry->NUMBER);
    }
#endif
    if ((pSelectEntry->magcal_iso < (pSelectEntry->limiting_mag-MAX_LIMITING_MAG)) &&
        (pSelectEntry->magcal_iso > pSelectEntry->max_bright_mag) &&
        (pSelectEntry->magcal_iso_rms < MAX_ISO_RMS)) {
      /* Good point, select it */
      input_index++;
    } else {
      /* Bad point, replace with the the item at the end of the array */
      memcpy(pSelectEntry,pInputEntry,sizeof(INPUTENTRY));
      if (pSelectEntry->magcal_iso < pSelectEntry->max_bright_mag+MAX_BRIGHT_ADJUST) {
        bright_reject_count++;
#if 0
        printf("Rejecting (1) %s magcal_iso %f magcal_iso_rms %f limiting_mag %f max_bright_mag %f\n",
               pSelectEntry->REF,
               pSelectEntry->magcal_iso,
               pSelectEntry->magcal_iso_rms,
               pSelectEntry->limiting_mag,
               pSelectEntry->max_bright_mag);
#endif
      } else if (pSelectEntry->magcal_iso_rms > MAX_ISO_RMS) {
#if 0
        printf("Rejecting (2) %s magcal_iso %f magcal_iso_rms %f limiting_mag %f max_bright_mag %f\n",
               pSelectEntry->REF,
               pSelectEntry->magcal_iso,
               pSelectEntry->magcal_iso_rms,
               pSelectEntry->limiting_mag,
               pSelectEntry->max_bright_mag);
#endif
        max_rms_count++;
      }
      pCommon->np--;
    }
  }
  printf("Read points %d, bright reject %d high rms %d output %d plot_type %d for %s\n",pCommon->nlines,bright_reject_count,max_rms_count,pCommon->np,plot_type,input_name);
  pCommon->nlines = pCommon->np;




  /* Allocate our arrays */
  pCommon->pSelectTable1 = (PINPUTENTRY)calloc(pCommon->nlines+1,sizeof(INPUTENTRY));
  pCommon->pSelectTable2 = (PINPUTENTRY)calloc(pCommon->nlines+1,sizeof(INPUTENTRY));

  if ((pCommon->pInputTable == NULL) ||
      (pCommon->pSelectTable1 == NULL) ||
      (pCommon->pSelectTable2 == NULL)) {
    fprintf(stderr,"ERROR: Memory allocation failure for %d points\n",pCommon->nlines);
  }



  /* The X and Y limits are now given as inputs */
  xmin = 0.0;
  xmax = 1.0 * mx;
  ymin = 0.0;
  ymax = 1.0 * my;






  dx=(xmax-xmin)/(1.0 * nx);
  dy=(ymax-ymin)/(1.0 *ny);

  if (dx > dy) {
    scale = dx * 1.5;
  } else {
    scale = dy * 1.5;
  }



  /* Report the input values for checking */
  if (verbose) {
    fprintf(stdout,"input_name= %s\n",input_name);
    fprintf(stdout,"xmin= %f  xmax= %f\n",xmin,xmax);
    fprintf(stdout,"ymin= %f ymax= %f\n",ymin,ymax);
    fprintf(stdout,"nx= %d ny= %d\n",nx,ny);
    fprintf(stdout,"grid-spacing: X= %f Y= %f\n",dx,dy);
    fprintf(stdout,"plate scale (arcsec/pixel) = %f\n",arcsecPerPixel);
    fprintf(stdout,"smoothing scale (pixels)= %f\n",scale);
    fprintf(stdout,"xtext= %s ytext= %s ncon= %d\n",xtext,ytext,ncon);
    fprintf(stdout,"plotdev= %s\n",plotdev);
    fprintf(stdout,"clip= %f\n",clip);	
  }

  /* Tranformation Vector for plotting. */
  TR[0]=xmin-(dx/2);
  TR[1]=dx;
  TR[2]=0;
  TR[3]=ymin-(dy/2);
  TR[4]=0;
  TR[5]=dy;


  /*Initialize the 2D Arrays */
  for (i=0; i < nx; i++) {
    x1=xmin+(cvtreal(i)*dx)+ (0.5*dx);
    for (j=0; j<ny; j++) {
      y1=ymin+(cvtreal(j)*dy) + (0.5*dy);
      MED[i][j]=0.0;
      SIGMA[i][j]=0.0;
      POINT[i][j]=0.0;
      bright_magcor_local[i][j] = 0.0;
      bright_magcor_error[i][j] = 0.0;
      bright_npoints_local[i][j] = 0;
      MAP[i][j]=0.0;
      pCommon->spatial_bin2[i][j] = CalculateBin(mx,my,x1,y1,&edgeDist);
    }
  }


  /*  Sample (bin) the points, using a boxcar of radius "scale".*/
  if (verbose) {
    fprintf(stdout,"Creating map of input data...\n");
  }
  for (i=0; i < nx; i++) {
    x1=xmin+(cvtreal(i)*dx)+ (0.5*dx);
    for( j=0; j < ny; j++) {
      y1=ymin+(cvtreal(j)*dy) + (0.5*dy);

#if 0
      if ((i == 24) && (j == 29)) {
        printf("At ix %d iy %d\n",i,j);
        printFlag = 1;
      } else {
        if (printFlag == 1) {
          printFlag = 0;
        }
      }
#endif

      iteration = 0;
      brightFlag = 0;
      CalcParameters(pCommon,i,j,x1,y1,scale,clip,iteration,brightFlag,&POINT[i][j],&MED[i][j],&SIGMA[i][j]);
      brightFlag = 1;
      CalcParameters(pCommon,i,j,x1,y1,scale,clip,iteration,brightFlag,&bright_npoints_local[i][j],&bright_magcor_local[i][j],&bright_magcor_error[i][j]);



    }
  }



  if (verbose) {
    fprintf(stdout,"Correcting input data...\n");
  }


  for (k=0; k < pCommon->np; k++) {
    pInputEntry = &pCommon->pInputTable[k];
#if 0
    if (pInputEntry->NUMBER == 10942) {
      printf("At %d\n",pInputEntry->NUMBER);
    }
#endif
#ifdef BOUNDS_CHECK
    if ((k >= pCommon->nlines) ||
        (k < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",k,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */

    i = (pInputEntry->X_IMAGE-xmin)/dx;
    j = (pInputEntry->Y_IMAGE-ymin)/dy;
    if (i < 0) {
      i = 0;
    }
    if (i >= nx) {
      i = nx-1;
    }
    if (j < 0) {
      j = 0;
    }
    if (j >= ny) {
      j = ny - 1;
    }
    if (pInputEntry->local_bin != i + (nx*j)) {
      printf("WARNING: local bin mismatch %d %d NUMBER %d for %s\n",pInputEntry->local_bin,i + (nx*j),pInputEntry->NUMBER,input_name);
      i = pInputEntry->local_bin % nx;
      j = pInputEntry->local_bin / nx;
    }
    

    pInputEntry->magcor_local=MED[i][j];
    pInputEntry->magcal_local_error=SIGMA[i][j];
    pInputEntry->npoints_local=POINT[i][j];
  }

#if 0

  for (i=0; i < pCommon->np; i++) {
    pInputEntry = &pCommon->pInputTable[i];
#ifdef BOUNDS_CHECK
    if ((i >= pCommon->nlines) ||
        (i < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",i,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */

    printf("%d %f %f %f %d %d\n",i+1,
           pInputEntry->magcor_local,
           pInputEntry->magcor_local_error,
           pInputEntry->npoints_local,
           pInputEntry->local_bin);

  }
#endif
  


  /* Create a map of the corrected data, for comparison with the input data. Plot them with identical scaling */
  if (verbose) {
    fprintf(stdout,"Creating map of corrected data...\n");
  }
  /* Change our sort value to magcor_local */
  for (i=0; i < pCommon->nlines; i++) {
    pInputEntry = &pCommon->pInputTable[i];
#ifdef BOUNDS_CHECK
    if ((i >= pCommon->nlines) ||
        (i < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",i,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */
    pInputEntry->sortval = pInputEntry->magcor_local;
    
  }
  

  for (i=0; i < nx; i++) {
    x1=xmin+(cvtreal(i)*dx)+(0.5*dx);
    for (j=0; j < ny; j++) {
      y1=ymin+(cvtreal(j)*dy + (0.5*dx));
      iteration = -1;
      brightFlag = 0;
      CalcParameters(pCommon,i,j,x1,y1,scale,clip,iteration,brightFlag,&goodval2,&med,&rms);
#if 0
      printf("%d %d %f\n",
             i+1,j+1,med);
#endif
      MAP[i][j]=med;
    }
  }

  /* Write out the corrected data to a file */
  if (verbose) {
    fprintf(stdout,"Writing corrected data to file...\n");
  }


  for (i=0; i < pCommon->np; i++) {
    pInputEntry = &pCommon->pInputTable[i];
#if 0
    if (pInputEntry->NUMBER == 10942) {
      printf("At %d\n",pInputEntry->NUMBER);
    }
#endif
#ifdef BOUNDS_CHECK
    if ((i >= pCommon->nlines) ||
        (i < 0)) {
      printf("ERROR: Bounds Check Failure %d %d at %s %d\n",i,pCommon->nlines,__FUNCTION__,__LINE__);
      exit(-1);
    }
#endif /* BOUNDS_CHECK */
#if 0
    fprintf(outHandle,"REF\tra\tdec\trefmag\tcolor\tMAG_ISO\tNUMBER\tBFLAGS\tX_IMAGE\tY_IMAGE\textinction\tmagcal_iso\tdmag_iso\tmagcal_iso_rms\tspatial_bin\tlimiting_mag\tmagcor_local\tmagcal_local_error\tnpoints_local\tlocal_bin\tmagcal_local\tdmag_local\tcal_local\n");
#endif
    magcal_local = pInputEntry->magcal_iso - pInputEntry->magcor_local;
    dmag_local = pInputEntry->dmag_iso - pInputEntry->magcor_local;
    fprintf(outHandle,"%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\t%d\t%d\t%.6f\t%.6f\t%d\n",
            pInputEntry->REF,
            pInputEntry->ra,
            pInputEntry->dec,
            pInputEntry->refmag,
            pInputEntry->color,
            pInputEntry->MAG_ISO,
            pInputEntry->NUMBER,
            pInputEntry->BFLAGS,
            pInputEntry->X_IMAGE,
            pInputEntry->Y_IMAGE,
            pInputEntry->extinction,
            pInputEntry->magcal_iso,
            pInputEntry->dmag_iso,
            pInputEntry->magcal_iso_rms,
            pInputEntry->spatial_bin,
            pInputEntry->limiting_mag,
            pInputEntry->magcor_local,
            pInputEntry->magcal_local_error,
            pInputEntry->npoints_local,
            pInputEntry->local_bin,
            magcal_local,
            dmag_local,
            cal_local);
  }
  /* Now write out the binning information to the grid file 
   * Make it starbase compatible 
   */


  fprintf(gridHandle,"nx\tny\tix\tiy\tzout\terrout\tnpout\tbright_magcor_local\tbright_magcor_error\tbright_npoints_local\n");
  fprintf(gridHandle,"--\t--\t--\t--\t----\t------\t-----\t-------------------\t-------------------\t--------------------\n");

  for (i=0; i < nx; i++) {
    for (j=0; j < ny; j++) {
      fprintf(gridHandle,"%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.0f\t%.6f\t%.6f\t%.0f\n",
              nx,ny,i,j,MED[i][j],SIGMA[i][j],POINT[i][j],bright_magcor_local[i][j],bright_magcor_error[i][j],bright_npoints_local[i][j]);
    }
  }

  fclose(outHandle);
  fclose(gridHandle);


  /* Find high and low values for plotting. */
  medmin=1.0e6;
  medmax=-1.0e6;
  sigmin=1.0e6;
  sigmax=-1.0e6;
  npmin=1.0e6;
  npmax=-1.0e6;
  mapmin=1.0e6;
  mapmax=-1.0e6;
  for (i=0; i < nx; i++) {
    for (j=0; j < ny; j++) {
      if(MED[i][j] > medmax) medmax=MED[i][j];
      if(MED[i][j] < medmin) medmin=MED[i][j] ;
      if((SIGMA[i][j] < 90.0) && (SIGMA[i][j] > sigmax)) sigmax=SIGMA[i][j];
      if((SIGMA[i][j] < 90.0) && (SIGMA[i][j] < sigmin)) sigmin=SIGMA[i][j] ;
      if(POINT[i][j] > npmax) npmax=POINT[i][j];
      if(POINT[i][j] < npmin) npmin=POINT[i][j] ;
      if(MAP[i][j] > mapmax) mapmax=MAP[i][j];
      if(MAP[i][j] < mapmin) mapmin=MAP[i][j]			 ;
    }
  } 
  if (verbose) {
    fprintf(stdout,"mapmin %f mapmax %f\n",mapmin,mapmax);
    fprintf(stdout,"medmin %f medmax %f\n",medmin,medmax);
  }

  /*	mapmin=medmin; */
  /*	mapmax=medmax; */

  /* set contour levels */
  if(ncon > 1)  { 
    dz=(medmax-medmin)/cvtreal(ncon);
    if (verbose) {
      fprintf(stdout,"Contour Levels for input data-map\n");
    }
    for (i=0; i < ncon; i++) {
      CON[i]=medmin+(cvtreal(i)*dz);
      if (verbose) {
        fprintf(stdout,"%f\n",CON[i]);
      }
    }
    if (verbose) {
      fprintf(stdout,"------------------------\n");
    }
  }


  /* PLOT THE RESULTS! */
  xmin= xmin-dx/2.;
  ymin= ymin-dy/2.;
  xmax= xmax+dx/2.;
  ymax= ymax+dy/2.;
  if (verbose) {
    fprintf(stdout,"%f %f %f %f\n",xmin,ymin,xmax,ymax);
    fprintf(stdout,"%f %f\n",medmax, medmin);
  }
#ifdef DO_PLOT
  if (plotFlag == 1) {
    /* open graphics device: */
    cpgopen(plotdev);
    cpgslw (4);
    cpgsch (1.5);

    /* subdivide the page into 4 panels;*/
    cpgsubp(2,2);

    /* Panel(1) : Clipped Median */
    cpgpanl(2,2);
    cpgenv(xmin, xmax, ymin, ymax,  0,  0);
    cpgmtxt ("B",2.5,0.5,0.5,xtext);
    cpgmtxt ("L",2.0,0.5,0.5,ytext);
    cpgmtxt ("T",3.0,0.0,0.0,plotdev);
    cpgmtxt ("T",1.0,0.5,0.5,"Clipped Median Map of Dataset");
    cpgmtxt ("T",1.0,1.0,0.0,"Median");

    /* Invert axes in coverting from C to Fortran */
    for (i=0; i < nx; i++) {
      for (j=0; j < ny; j++) {
        switch (plot_type) {
        case PLOT_TYPE_N:
          PLOTOUT_N[i][j] = MED[j][i];
          break;
        case PLOT_TYPE_SL:
          PLOTOUT_SL[i][j] = MED[j][i];
          break;
        case PLOT_TYPE_LS:
          PLOTOUT_LS[i][j] = MED[j][i];
          break;
        default:
          printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
          break;
        }
      }
    }
    switch (plot_type) {
    case PLOT_TYPE_N:
      cpggray ((float *)&PLOTOUT_N[0][0],nx,ny,1,nx,1,ny,medmax,medmin,TR);
      break;
    case PLOT_TYPE_SL:
      cpggray ((float *)&PLOTOUT_SL[0][0],nx,ny,1,nx,1,ny,medmax,medmin,TR);    
      break;
    case PLOT_TYPE_LS:
      cpggray ((float *)&PLOTOUT_LS[0][0],nx,ny,1,nx,1,ny,medmax,medmin,TR);          
      break;
    default:
      printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
      break;
    }
 
    cpgwedg("RG", 1.0, 4.0, medmax, medmin, "");
    cpgsls(1);
    if(ncon > 1)  { 
      ncon = ncon * -1;
      if (verbose) {
        fprintf(stdout,"ncon=%d\n",ncon);
      }
      switch (plot_type) {
      case PLOT_TYPE_N:
        cpgcont ((float *)&PLOTOUT_N[0][0],nx,ny,1,nx,1,ny,CON,ncon,TR);
        break;
      case PLOT_TYPE_SL:
        cpgcont ((float *)&PLOTOUT_SL[0][0],nx,ny,1,nx,1,ny,CON,ncon,TR);
        break;
      case PLOT_TYPE_LS:
        cpgcont ((float *)&PLOTOUT_LS[0][0],nx,ny,1,nx,1,ny,CON,ncon,TR);
        break;
      default:
        printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
        break;
      }

    }

    /* Panel(2) : Clipped SIGMA */
    cpgpanl(1,2);
    cpgenv(xmin, xmax, ymin, ymax,  0,  0);
    cpgmtxt ("B",2.5,0.5,0.5,xtext);
    cpgmtxt ("L",2.0,0.5,0.5,ytext);
    cpgmtxt ("T",1.0,0.5,0.5,"Scatter of Calibration Points");
    cpgmtxt ("T",1.0,1.0,0.0,"Sigma");
    for (i=0; i < nx; i++) {
      for (j=0; j < ny; j++) {
        switch (plot_type) {
        case PLOT_TYPE_N:
          if (SIGMA[i][j] < 90.0) {
            PLOTOUT_N[i][j] = SIGMA[j][i];
          } else {
            PLOTOUT_N[i][j] =  sigmax; 
          }

          break;
        case PLOT_TYPE_SL:
          if (SIGMA[i][j] < 90.0) {
            PLOTOUT_SL[i][j] = SIGMA[j][i];
          } else {
            PLOTOUT_SL[i][j] =  sigmax; 
          }

          break;
        case PLOT_TYPE_LS:
          if (SIGMA[i][j] < 90.0) {
            PLOTOUT_LS[i][j] = SIGMA[j][i];
          } else {
            PLOTOUT_LS[i][j] =  sigmax; 
          }

          break;
        default:
          printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
          break;
        }
      }
    }
    switch (plot_type) {
    case PLOT_TYPE_N:
      cpggray ((float *)&PLOTOUT_N[0][0],nx,ny,1,nx,1,ny,sigmax,sigmin,TR);
      break;
    case PLOT_TYPE_SL:
      cpggray ((float *)&PLOTOUT_SL[0][0],nx,ny,1,nx,1,ny,sigmax,sigmin,TR);
      break;
    case PLOT_TYPE_LS:
      cpggray ((float *)&PLOTOUT_LS[0][0],nx,ny,1,nx,1,ny,sigmax,sigmin,TR);
      break;
    default:
      printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
      break;
    }

    cpgwedg("RG", 1.0, 4.0, sigmax, sigmin, "");

    /* Panel(3) : Number of good points contributing to each grid-cell */
    cpgpanl(2,1);
    cpgenv(xmin, xmax, ymin, ymax,  0,  0);
    cpgmtxt ("B",2.5,0.5,0.5,xtext);
    cpgmtxt ("L",2.0,0.5,0.5,ytext);
    cpgmtxt ("T",1.0,0.5,0.5,"Number of Points Contributing to each Grid-Cell")         ;
    cpgmtxt ("T",1.0,1.0,-1.5,"N");
    for (i=0; i < nx; i++) {
      for (j=0; j < ny; j++) {
        switch (plot_type) {
        case PLOT_TYPE_N:
          PLOTOUT_N[i][j] = POINT[j][i];
          break;
        case PLOT_TYPE_SL:
          PLOTOUT_SL[i][j] = POINT[j][i];
          break;
        case PLOT_TYPE_LS:
          PLOTOUT_LS[i][j] = POINT[j][i];
          break;
        default:
          printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
          break;
        }

      }
    }
    switch (plot_type) {
    case PLOT_TYPE_N:
      cpggray ((float *)&PLOTOUT_N[0][0],nx,ny,1,nx,1,ny,npmax,npmin,TR);
      break;
    case PLOT_TYPE_SL:
      cpggray ((float *)&PLOTOUT_SL[0][0],nx,ny,1,nx,1,ny,npmax,npmin,TR);
      break;
    case PLOT_TYPE_LS:
      cpggray ((float *)&PLOTOUT_LS[0][0],nx,ny,1,nx,1,ny,npmax,npmin,TR);
      break;
    default:
      printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
      break;
    }

    cpgwedg("RG", 1.0, 4.0, npmax, npmin, "");

    /* Panel(4) : Map of corrected datapoints, Median */
    cpgpanl(1,1);
    cpgenv(xmin, xmax, ymin, ymax,  0,  0);
    cpgmtxt ("B",2.5,0.5,0.5,xtext);
    cpgmtxt ("L",2.0,0.5,0.5,ytext);
    cpgmtxt ("T",1.0,0.5,0.5,"Corrected Data Map");
    cpgmtxt ("T",1.0,1.0,0.0,"Median");
    for (i=0; i < nx; i++) {
      for (j=0; j < ny; j++) {
        switch (plot_type) {
        case PLOT_TYPE_N:
          PLOTOUT_N[i][j] = MAP[j][i];
          break;
        case PLOT_TYPE_SL:
          PLOTOUT_SL[i][j] = MAP[j][i];
          break;
        case PLOT_TYPE_LS:
          PLOTOUT_LS[i][j] = MAP[j][i];
          break;
        default:
          printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
          break;
        }

      }
    }
    switch (plot_type) {
    case PLOT_TYPE_N:
      cpggray ((float *)&PLOTOUT_N[0][0],nx,ny,1,nx,1,ny,mapmax,mapmin,TR);
      break;
    case PLOT_TYPE_SL:
      cpggray ((float *)&PLOTOUT_SL[0][0],nx,ny,1,nx,1,ny,mapmax,mapmin,TR);
      break;
    case PLOT_TYPE_LS:
      cpggray ((float *)&PLOTOUT_LS[0][0],nx,ny,1,nx,1,ny,mapmax,mapmin,TR);
      break;
    default:
      printf("ERROR: contour_plot.c software bug in line %d\n",__LINE__);
      break;
    }

    cpgwedg("RG", 1.0, 4.0, mapmax, mapmin, "");


    cpgiden();

    cpgend();
  }
#endif /* DO_PLOT */


  if (pCommon->pSelectTable1 != NULL) {
    free(pCommon->pSelectTable1);
  }
  if (pCommon->pSelectTable2 != NULL) {
    free(pCommon->pSelectTable2);
  }
  if (input_table != NULL) {
    Free(input_table);
    input_table = NULL;
    pCommon->pInputTable = NULL;
  }
  
  return(0);
}

