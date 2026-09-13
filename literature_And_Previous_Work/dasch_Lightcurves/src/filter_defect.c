// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* filter_defect.c
 *
 *  gcc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include -I /dasch/install/pgplot  filter_defect.c pipelineutils.a -L /dasch/install/lib -lm  -L/usr/lib${lib64}/mysql -L /dasch/install/pgplot -lmysqlclient -ldl -pthread  -ltable -lutil -lwcs  -lcpgplot -lpgplot  -L/usr/X11R6/lib  -lX11 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 /usr/lib${lib64}/libg2c.so.0  -o filter_defect 

 *
 * Jun  3, 2008 Edward J. Los - Initial version adapted from Sumin Tang's my_dasch_filter_2sigma.m script of May 27, 2008.
 * Jul  7, 2008 Edward J. Los - Apply the defect filter to all objects for which FLUX_MAX <= 1.75 * THRESHOLD)
 * Aug  6, 2008 Edward J. Los - Change 1.75 * THRESHOLD filtering to 2.3 * THRESHOLD
 * Aug 22, 2008 Edward J. Los - Move THRESHOLD_FACTOR to pipelineutils.h
 * Sep 26, 2008 Edward J. Los - change ixbind to include all stars so that high drad and blend stars are also checked by the defect filter
 * Jan 27, 2009 Edward J. Los - split FLAGS into AFLAGS and BFLAGS
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 * Feb 27, 2009 Edward J. Los - Check FILTER_AFLAG_DRADBIN when checking FILTER_AFLAG_DRAD
 * Jul  9, 2012 Edward J. Los - Modify for Sumin Tang's dasch_defectfilter_20120628.m (Documented in memoranda of Sun 6/24/12 6:13 PM 
 *                              with changes from memorandum of Fri 6/29/12 9:43 AM)
 *                                Compared with the 2008 version, the following improvements are made:
 *                                1. 2-sigma threshold revised to be 3-sigma (2.5sigma for theta and ellipticity)
 *                                2. theta vs ellipticity: corrected for discontinuity in theta (0 deg = 180 deg)
 *                                3. revised defination of real objects (position drad within 2-sigma)
 *                                4. revised binning in both spatial bins and parameter bins (more dynamic now)
 *                                5. add cut on FLUX_MAX vs ISO0 (both x and y-axis)
 *                              The THRESHOLD_FACTOR is now removed.
 * Sep  3, 2013 Edward J. Los - Do not use high proper motion error objects to compute drad
 * 
 *
 * Note: The following subarrays parallel those in the script of May 27, 2008
 *
 *    ixbin[]  of size ixbinCount    All but high drad and low FLUX_MAX
 *    ixbin2[] of size ixbin2Count   Members of ixbin[] not near edges and corners
 *    ixbinr[] of size ixbinrCount   Members of ixbin2[] with low drad (Now equal to ixbin2 with the use of filterblended DRAD flags)
 *    ixbinf[] of size ixfinfCount   Members of ixbin2[] with high drad (Now null - was unused anyways)
 * 
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include "cpgplot.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

double round(double x);
/* #define LOS_DEBUG 1 */
/* #define LOS_DEBUG_NUMBER 174458 */
#define MAX_BUFFER 512
#define MAX_FLUX_LEVELS 20
#define MAX_MAG_LEVELS 30
#define MAX_LEVEL_VECTOR 200

#if 0 /* Original octave flags */
/* #define FLAGS_OUTPUT 1 Output the new Sextractor flags value - clear this to emulate the Octave script */
/* #define NO_DRAD_FILTER 1 Do not perform drad filtering - clear this to emulate the Octave script */
#define USE_XY_RANGE 1  /* Set this to emulate the Octave script */
#else /* New flags */
#define FLAGS_OUTPUT 1 /* Output the new Sextractor flags value - clear this to emulate the Octave script */
#define NO_DRAD_FILTER 1 /* Do not perform drad filtering - clear this to emulate the Octave script */
/* #define NO_DRAD_COLUMN 1  */ /* No drad or AFLAGS column */
/* #define USE_XY_RANGE 1  */ /* Set this to emulate the Octave script */

#endif

#define INTERNAL_FLAG_IXF1    1   /* If set, the object passes test #1: fwhm is near the average */
#define INTERNAL_FLAG_IXF2    2   /* If set, the object passes test #2 */
#define INTERNAL_FLAG_IXF3    4   /* If set, the object passes test #3 */
#define INTERNAL_FLAG_IXF4    8   /* If set, the object passes test #4 */
#define INTERNAL_FLAG_IXF5   16   /* If set, the object passes test #5 */
#if 1
#define INTERNAL_FLAG_IXFALL 31   /* All above are set - accept this object */
#else
#define INTERNAL_FLAG_IXFALL  8   /* All above are set - accept this object */
#endif

#define INTERNAL_FLAG_FLUX   32   /* If set, the object is below the flux threshold */
#define INTERNAL_FLAG_REJECT 64   /* If set, the object had either failed drad testing in filterblended.c
                                     or is considered a blended star */



typedef struct _starimage {
  int NUMBER;         /* Sextractor reference number */
  double X_IMAGE;        /* Sextractor X location in pixels */
  double Y_IMAGE;        /* Sextractor Y location in pixels */
  double drad;        /* Matching error */
  double FLUX_MAX; 
  double FLUX_ISO;
  double MAG_ISO;     /* Sextractor isophotonic magnitude */
  double FWHM_WORLD;
  double THETA_J2000;
  double ELLIPTICITY;
  double plate_dist;
  int ISO0;
  int AFLAGS;          /* flags word */
  int BFLAGS;          /* flags word */
  /* Aliases for the match dataset that appear in the script */
  double x; /* Alias for X_IMAGE */
  double y; /* Alias for Y_IMAGE */
  double THETA; /* Alias for THETA_J2000; */
  double magr;
  double fwhmr;
  double thetar;
  double er;
  double fluxmaxr;
  double iso0r;
  double thetadistr;

  double magf1;
  double fwhmf1;
  double thetaf1;
  double ef1;
  double fluxmaxf1;
  double iso0f1;

  /* Aliases for the sextractor dataset that appear in the script */
  int NUMBER3;
  double x3; /* Alias for X_IMAGE */
  double y3; /* Alias for Y_IMAGE */
  double FLUX_MAX3;
  double MAG_ISO3;
  double FWHM_WORLD3;
  double THETA_J20003;
  double ELLIPTICITY3;
  int ISO03;
  double plate_dist3;
  double THETA3; /* Alias for THETA_J20003; */
  double mag;
  double fwhm;
  double theta;
  double e;
  double fluxmax;
  double iso0;
  double num;
  double thetadist;

  /* Internal flags */
  int internal;        /* Internal flags */
  int flag_match;      /* Object appears in the match table if nonzero */

} STARIMAGE,*PSTARIMAGE;

/* Sort routine based on the Sextractor reference number */
int ImageCompare(const void *first, const void *second) 
{
  int numberFirst = ((PSTARIMAGE)first)->NUMBER;
  int numberSecond = ((PSTARIMAGE)second)->NUMBER;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}
/* Find the mean of the vector */
double FindMean(double *vector,int count) 
{
  double sum = 0;
  int index;
  double result;
  if (count == 0) {
    return(0.0);
  }
  for (index = 0; index < count; index++) {
    sum += vector[index];
  }
  result = sum/(1.0*count);
  return(result);

}
/* Find the standard deviation of the vector */
double FindStd(double *vector,int count) 
{
  double mean;
  double sum = 0;
  int index;
  double result;
  if (count <= 1) {
    return(0.0);
  }
  mean = FindMean(vector,count);

  for (index = 0; index < count; index++) {
    sum += sqr(vector[index] - mean);
  }
  result = sqrt(sum/(1.0*(count-1)));
  return(result);

}
/*
 * Create a percentvector of of percentages of the form [0:100*15/length(magr):80/nmag, 100/nmag:100/nmag:100]
 * then find a magvector such that the percentiles of vector correspond to the percentages of precentvector
 *
 */
void prctile(double *vector,int vectorCount,int nmag,double *vectorVector,int *magvectorsize,int invertFlag,int printFlag) 
{
  double percentVector[MAX_LEVEL_VECTOR];

  double *vectorpercentVector;
  int sizepercentVector = nmag;
  int increment;
  int index;
  int vectorIndex;
  memset(percentVector,0,sizeof(percentVector));
  increment = ((80.0*(1.0*vectorCount))/(100.0*15.0*(1.0*nmag)))+1.0;
  sizepercentVector = nmag+increment;
  if (sizepercentVector > MAX_LEVEL_VECTOR) {
    printf("ERROR: prctile sizepercentVector %d too large for %d\n",sizepercentVector,MAX_LEVEL_VECTOR);
    exit(-1);
  }
  for (index = 0; index < sizepercentVector; index++) {
    if (index < increment) {
      percentVector[index] = (1.0*index)*(100.0*15.0)/(1.0*vectorCount);
    } else {
      percentVector[index] = (1.0*(index-increment+1))*(100.0)/(1.0*nmag);
    }
    if(invertFlag != 0) {
      percentVector[index] = 100.0-percentVector[index];
    }
  }
  qsort((void*)percentVector,sizepercentVector,sizeof(double),dcmp); 
  qsort((void*)vector,vectorCount,sizeof(double),dcmp);
  vectorpercentVector = (double *)calloc(vectorCount,sizeof(double));
  if (vectorpercentVector == NULL) {
    printf("ERROR: failed to allocate vectorpercentVector of size %d\n",vectorCount);
    exit(-1);
  }
  for (vectorIndex = 0; vectorIndex < vectorCount; vectorIndex++) {
    vectorpercentVector[vectorIndex] = 100*(vectorIndex+0.5)/(1.0*vectorCount);
  }
  /* Now fill in the final vector, interpolating percentages */
  for (index = 0; index < sizepercentVector; index++) {
    if (index == 0) {
      vectorVector[index] = vector[0];
    } else if (index == (sizepercentVector-1)) {
      vectorVector[index] = vector[vectorCount-1];
    } else {
      for (vectorIndex = 0; vectorIndex < vectorCount-1; vectorIndex++) {
        if ((vectorpercentVector[vectorIndex] <= percentVector[index]) &&
            (vectorpercentVector[vectorIndex+1] > percentVector[index])) {
          /* Add 0.0001 to the result to avoid having everything fall on a boundary */
          vectorVector[index] = vector[vectorIndex]+ 0.0001 + (percentVector[index]-vectorpercentVector[vectorIndex])*(vector[vectorIndex+1]-vector[vectorIndex])/(vectorpercentVector[vectorIndex+1]-vectorpercentVector[vectorIndex]);
          break;
        }

      }
    }
    if (printFlag) {
      printf("index %d, percent %f, vector %f\n",index,percentVector[index],vectorVector[index]);
    }
  }




  *magvectorsize = sizepercentVector;

  if (vectorpercentVector != NULL) {
    free(vectorpercentVector);
  }
  return;

}
double checklog10(double val)
{
  double temp = log10(val);
  if (isnan(temp)) {
    temp = 0;
  }
  return(temp);
}

int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char infile[MAX_BUFFER];
  char cmdchar;
  int errorFlag = 0;
  int numBins;
  char outfile[MAX_BUFFER];
  char workingDirectory[MAX_BUFFER];
  char *charPtr;
  char plotdev[MAX_BUFFER];
  char summaryfile[MAX_BUFFER];
  char fileroot[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  FILE *outHandle = NULL;
  FILE *summaryHandle = NULL;
  char title[MAX_BUFFER];

#ifdef LOS_DEBUG
  char debug_name[MAX_BUFFER];
  FILE *debugHandle;
#endif /* LOS_DEBUG */


  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  PSTARIMAGE sextractor_table = NULL;
  size_t sextractor_nrecs = 0;
  int sextractor_index;
  PSTARIMAGE pSextractor;

  File match_handle = NULL;
  char match_name[MAX_BUFFER];
  TableHead match_header = NULL;
  PSTARIMAGE match_table = NULL;
  size_t match_nrecs = 0;
  int match_index;
  PSTARIMAGE pMatch;

  int mosaicWidth = 0;
  int mosaicHeight = 0;

  double minimumX;
  double maximumX;
  double minimumY;
  double maximumY;
  double xrange;
  double yrange;
  double xmintest;
  double xmaxtest;
  double ymintest;
  double ymaxtest;

  int numpoints = 0;
  int irec = 0;

  int index;
  int indexMag;
  int indexFluxMax;
  int verbose = 0;
  int doPlots = 0;
  double *vector = NULL;
  int vectorCount;
  double maxmagiso = 0;
  double meddrad1;
  double stddrad1;
  double meddrad2;
  double stddrad2;
  double meddrad3;
  double stddrad3;

  double maggrid[MAX_LEVEL_VECTOR];
  double iso0grid[MAX_LEVEL_VECTOR];
  double fluxmaxgrid[MAX_LEVEL_VECTOR];
  double fluxmaxc[MAX_LEVEL_VECTOR];
  double rmsfluxmaxc[MAX_LEVEL_VECTOR];
  double iso0c[MAX_LEVEL_VECTOR];
  double rmsiso0c[MAX_LEVEL_VECTOR];
  double fwhmc[MAX_LEVEL_VECTOR];
  double rmsfwhmc[MAX_LEVEL_VECTOR];
  double ec[MAX_LEVEL_VECTOR];
  double rmsec[MAX_LEVEL_VECTOR];
  double x_bin[MAX_LEVEL_VECTOR];
  double y_bin[MAX_LEVEL_VECTOR];

  double fluxmaxc0;
  double rmsfluxmaxc0;
  int maggridsize;
  int maggridOK[MAX_LEVEL_VECTOR];
  int iso0gridOK[MAX_LEVEL_VECTOR];
  int iso0cOK[MAX_LEVEL_VECTOR];
  int iso0gridsize;
  int fluxmaxgridsize;

  int lowFluxCount = 0;
  int dradFilterCount = 0;
  int ixf1temp;
  int ixf2temp;
  int ixf4temp;

  int ixf1Count = 0;
  int ixf2Count = 0;
  int ixf3Count = 0;
  int ixf4Count = 0;
  int ixf5Count = 0;
  int ixfAllCount = 0;
  int ixfFlaggedCount = 0;
  int kxr = 0;
  int kall;
  int ixfCount = 0;
  int nbin;
  int iii;
  int jjj;
  int ibin;
  int matchIndex;
  int *matchVector;
  int matchVectorCount = 0;
  int *ixbinf1;
  int ixbinf1Count = 0;
  int *ixa;
  int ixaCount = 0;
  int *ixb;
  int ixbCount = 0;
  int *ixbin;
  int ixbinCount;
  int ixbinTotal = 0;
  int *ixbin2;
  int ixbin2Count;
  int *ixbinf;
  int ixbinfCount;
  int *ixbinr;
  int ixbinrCount;
  int niso0 = 0;
  double max_plate_dist = 0;
  double minmag;
  double maxmag;
  double minmagr;
  double maxmagr;
  double minmagd;
  double maxmagd;
  double dmag;
  int nmag;
  int *ix;
  int ixCount;
  int *ix0;
  int ix0Count;

  time_t startTime;
  time_t curTime;
  double edgeDist;
  double fwhmc0;
  double rmsfwhmc0;
  double ec0;
  double stdec0;
  double rmsec0;
  double thetac0;
  double stdthetac0;
  double rmsthetac0;
  double rmsthetac;
  double rmsthetac0Sum;
  double iso0c0;
  double rmsiso0c0;

  double minfluxmax;
  double maxfluxmax;
  double minfluxmaxr;
  double maxfluxmaxr;
  double minfluxmaxd;
  double maxfluxmaxd;
  int nfluxmax;
  double dfluxmax;
  int ixsrCount = 0;
  int ixsfCount = 0;
  int flagsrCount = 0;
  int flagsfCount = 0;
  int ixf1old = 0;
  int ixf2old = 0;
  int ixf3old = 0;
  int ixf4old = 0;
  int ixf5old = 0;
  int ixf1old2 = 0;
  int ixf2old2 = 0;
  int ixf3old2 = 0;
  int ixf4old2 = 0;
  int ixf5old2 = 0;

  float xmin;
  float xmax;
  float ymin;
  float ymax;
  float valx;
  float valy;
  float minthetar;
  float maxthetar;
  int plot_nrecs1;
  int plot_index1;
  float *xval1 = NULL;
  float *yval1 = NULL;

  int plot_nrecs2;
  int plot_index2;
  float *xval2 = NULL;
  float *yval2 = NULL;
#ifdef LOS_DEBUG_NUMBER
  PSTARIMAGE pDebugSextractor = NULL;
  float xval1Debug = 0;
  float yval1Debug = 0;
#endif /* LOS_DEBUG_NUMBER */
#ifdef NO_DRAD_COLUMN
  printf("ERROR: NO_DRAD_COLUMN is set\n");
#endif /* NO_DRAD_COLUMN */


  time(&startTime);

  match_name[0] = 0;
  sextractor_name[0] = 0;
  outfile[0] = 0;
  fileroot[0] = 0;
  qualifier[0] = 0;
  summaryfile[0] = 0;

 
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
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'r': /* summary file name */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(summaryfile,*++argv,MAX_BUFFER-2);
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
            }
          }
          break;



        case 'm': /* match file name */
        case 'M':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(match_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 's': /* sextractor file name */
        case 'S':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sextractor_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'p': /* file root */
        case 'P':
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

        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic width %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* mosaic height */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;




        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'g': /* verbose */
        case 'G':
          doPlots = 1;
          break;

        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */
  if (outfile[0] == 0) {
    printf("ERROR: No output filename was specified\n");
    errorFlag = 1;
  }
  if (summaryfile[0] == 0) {
    printf("ERROR: No summary filename was specified\n");
    errorFlag = 1;
  }
  if (match_name[0] == 0) {
    printf("ERROR: No match filename was specified\n");
    errorFlag = 1;
  }
  if (sextractor_name[0] == 0) {
    printf("ERROR: No sextractor filename was specified\n");
    errorFlag = 1;
  }
  if (fileroot[0] == 0) {
    printf("ERROR: No file root was specified\n");
    errorFlag = 1;
  }
  if (mosaicWidth == 0) {
    printf("ERROR: No mosaic width specified\n");
    errorFlag = 1;
  }
  if (mosaicHeight == 0) {
    printf("ERROR: No mosaic height specified\n");
    errorFlag = 1;
  }



  /* Open the sextractor file */
  sextractor_handle = Open(sextractor_name,"r");
  if (sextractor_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the sextractor file %s\n",sextractor_name);
  } else {
    if (verbose) {
      printf("Found sextractor file %s\n",sextractor_name);
    }
  }

  /* Open the match file */
  match_handle = Open(match_name,"r");
  if (match_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the match file %s\n",match_name);
  } else {
    if (verbose) {
      printf("Found match file %s\n",match_name);
    }
  }




  if (errorFlag) {
    printf("Usage: filter_defect -p <file root> \n");
    printf("                      -m <match file> \n");
    printf("                      -s <sextractor file> \n");
    printf("                      -o  <output file>\n");
    printf("                      -r  <summary file>\n");
    printf("                      -v  verbose\n");
    printf("                      -g  output graphs\n");

    return(-1);
  }



  /* Now read in the sextractor file */

  sextractor_header = table_header(sextractor_handle,TABLE_PARSE);
  if (sextractor_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",sextractor_name);
    return(-1);
  }

  sextractor_table = table_loadva(sextractor_handle,
                                  &sextractor_header,
                                  NULL, /* hbase */
                                  NULL, /* rows */
                                  NULL,
                                  sizeof(STARIMAGE),
                                  &sextractor_nrecs,
                                  TblInt,"NUMBER",TblOff(PSTARIMAGE,NUMBER),
                                  TblDbl,"X_IMAGE",TblOff(PSTARIMAGE,X_IMAGE),
                                  TblDbl,"Y_IMAGE",TblOff(PSTARIMAGE,Y_IMAGE),
                                  TblDbl,"FLUX_ISO",TblOff(PSTARIMAGE,FLUX_ISO),
                                  TblDbl,"FLUX_MAX",TblOff(PSTARIMAGE,FLUX_MAX),
                                  TblDbl,"MAG_ISO",TblOff(PSTARIMAGE,MAG_ISO),
                                  TblDbl,"FWHM_WORLD",TblOff(PSTARIMAGE,FWHM_WORLD),
                                  TblDbl,"THETA_J2000",TblOff(PSTARIMAGE,THETA_J2000),
                                  TblDbl,"ELLIPTICITY",TblOff(PSTARIMAGE,ELLIPTICITY),
                                  TblDbl,"plate_dist",TblOff(PSTARIMAGE,plate_dist),
                                  TblInt,"ISO0",TblOff(PSTARIMAGE,ISO0),

                                  TblInt,"AFLAGS",TblOff(PSTARIMAGE,AFLAGS),
                                  TblInt,"BFLAGS",TblOff(PSTARIMAGE,BFLAGS),

                                  0,"end",0);
  if (sextractor_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",sextractor_name);
    return(-1);
  }
  if (verbose) {
    printf("read %d records for %s\n",sextractor_nrecs,sextractor_name);
  }
  /* Sort the table according to sextractor NUMBER */
  qsort((void*)sextractor_table,sextractor_nrecs,sizeof(STARIMAGE),ImageCompare);
  for (index = 0; index < sextractor_nrecs; index++) {

    /* Use aliases to match the script */
    pSextractor = &sextractor_table[index];
    if (pSextractor->plate_dist > max_plate_dist) {
      max_plate_dist = pSextractor->plate_dist;
    }
    if (pSextractor->ISO0 <= 0) {
      pSextractor->ISO0 = 1;
    }
    if (pSextractor->FWHM_WORLD <= 0) {
      pSextractor->FWHM_WORLD = 1.0/3600.;
    }

    pSextractor->drad = 999.;        /* Matching error */
    pSextractor->internal = 0;        /* Internal flags */
    pSextractor->mag = pSextractor->MAG_ISO; /* MAG_ISO */
    
    pSextractor->fwhm = checklog10(3600.0 * pSextractor->FWHM_WORLD);  /* log10(3600.0*FWHM_WORLD) */
    pSextractor->THETA = pSextractor->THETA_J2000;
    pSextractor->theta = pSextractor->THETA;  /* THETA_J2000 */
    pSextractor->e = pSextractor->ELLIPTICITY;     /* Ellipticity */
#ifdef LOS_DEBUG_NUMBER
    if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
      printf("At NUMBER %d in line %d\n",pSextractor->NUMBER,__LINE__);
      pDebugSextractor = pSextractor;
      pDebugSextractor->magr = pDebugSextractor->MAG_ISO;
      pDebugSextractor->fwhmr = checklog10(3600.0 * pDebugSextractor->FWHM_WORLD);
      pDebugSextractor->er = pDebugSextractor->ELLIPTICITY;
      pDebugSextractor->THETA = pDebugSextractor->THETA_J2000;
      pDebugSextractor->thetar = pDebugSextractor->THETA;
      pDebugSextractor->fluxmaxr = checklog10(pDebugSextractor->FLUX_MAX);
      pDebugSextractor->iso0r = checklog10(1.0*pDebugSextractor->ISO0);
    }

#endif /* LOS_DEBUG_NUMBER */
    pSextractor->fluxmax = checklog10(pSextractor->FLUX_MAX); /* log10(FLUX_MAX) */
    pSextractor->iso0 = checklog10(1.0*pSextractor->ISO0); /* log10 ISO0 */
    pSextractor->flag_match = 0;
    pSextractor->x = pSextractor->X_IMAGE;
    pSextractor->y = pSextractor->Y_IMAGE;
    pSextractor->NUMBER3 = pSextractor->NUMBER;
    pSextractor->x3 = pSextractor->X_IMAGE;
    pSextractor->y3 = pSextractor->Y_IMAGE;
    pSextractor->FLUX_MAX3 = pSextractor->FLUX_MAX;
    pSextractor->MAG_ISO3 = pSextractor->MAG_ISO;
    pSextractor->FWHM_WORLD3 = pSextractor->FWHM_WORLD;
    pSextractor->THETA_J20003 = pSextractor->THETA_J2000;
    pSextractor->ELLIPTICITY3 = pSextractor->ELLIPTICITY;
    pSextractor->ISO03 = pSextractor->ISO0;
    pSextractor->plate_dist3 = pSextractor->plate_dist;
    pSextractor->THETA3 = pSextractor->THETA_J2000;
    pSextractor->magf1 = pSextractor->MAG_ISO;
    pSextractor->magf1 = pSextractor->MAG_ISO;
    pSextractor->fwhmf1 = checklog10(3600.0 * pSextractor->FWHM_WORLD);
    pSextractor->thetaf1 = pSextractor->THETA;
    pSextractor->ef1 = pSextractor->ELLIPTICITY;
    pSextractor->fluxmaxf1 = checklog10(pSextractor->FLUX_MAX);
    pSextractor->iso0f1 = checklog10(1.0*pSextractor->ISO0);
  }

  matchVector = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (matchVector == NULL) {
    printf("ERROR: failed to allocate matchVector for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  ixbinf1 = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ixbinf1 == NULL) {
    printf("ERROR: failed to allocate ixbinf1 for %d images\n",sextractor_nrecs);
    exit(-1);
  }


  xval1 = (float *)calloc(sextractor_nrecs,sizeof(float));
  if (xval1 == NULL) {
    printf("ERROR: failed to allocate xval1 for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  yval1 = (float *)calloc(sextractor_nrecs,sizeof(float));
  if (yval1 == NULL) {
    printf("ERROR: failed to allocate yval1 for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  xval2 = (float *)calloc(sextractor_nrecs,sizeof(float));
  if (xval2 == NULL) {
    printf("ERROR: failed to allocate xval2 for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  yval2 = (float *)calloc(sextractor_nrecs,sizeof(float));
  if (yval2 == NULL) {
    printf("ERROR: failed to allocate yval2 for %d images\n",sextractor_nrecs);
    exit(-1);
  }

  ixa = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ixa == NULL) {
    printf("ERROR: failed to allocate ixa for %d images\n",sextractor_nrecs);
    exit(-1);
  }

  ixb = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ixb == NULL) {
    printf("ERROR: failed to allocate ixb for %d images\n",sextractor_nrecs);
    exit(-1);
  }

  /* Vector for computing median values */
  vector = (double *)calloc(sextractor_nrecs,sizeof(double));
  if (vector == NULL) {
    printf("ERROR: failed to allocate vector for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  ixbin = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ixbin == NULL) {
    printf("ERROR: failed to allocate ixbin for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  ixbin2 = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ixbin2 == NULL) {
    printf("ERROR: failed to allocate ixbin2 for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  ixbinr = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ixbinr == NULL) {
    printf("ERROR: failed to allocate ixbinr for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  ixbinf = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ixbinf == NULL) {
    printf("ERROR: failed to allocate ixbinf for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  ix = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ix == NULL) {
    printf("ERROR: failed to allocate ix for %d images\n",sextractor_nrecs);
    exit(-1);
  }
  ix0 = (int *)calloc(sextractor_nrecs,sizeof(int));
  if (ix0 == NULL) {
    printf("ERROR: failed to allocate ix0 for %d images\n",sextractor_nrecs);
    exit(-1);
  }
 


  /* Now read in the match file */

  match_header = table_header(match_handle,TABLE_PARSE);
  if (match_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",match_name);
    return(-1);
  }

  match_table = table_loadva(match_handle,
                             &match_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(STARIMAGE),
                             &match_nrecs,
                             TblInt,"NUMBER",TblOff(PSTARIMAGE,NUMBER),
#ifndef NO_DRAD_COLUMN
                             TblDbl,"drad",TblOff(PSTARIMAGE,drad),
                             TblInt,"AFLAGS",TblOff(PSTARIMAGE,AFLAGS),
#endif /* NO_DRAD_COLUMN */
                             TblInt,"BFLAGS",TblOff(PSTARIMAGE,BFLAGS),
                             0,"end",0);
  if (match_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",match_name);
    return(-1);
  }
  if (verbose) {
    printf("read %d records for %s\n",match_nrecs,match_name);
  }
  /* Sort the table according to sextractor NUMBER */
  qsort((void*)match_table,match_nrecs,sizeof(STARIMAGE),ImageCompare);

  /* Now we step through the match table and copy matches to SExtractor */
  match_index = 0;
  pMatch = &match_table[match_index];
  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];
    while (pMatch->NUMBER < pSextractor->NUMBER) {
      match_index++;
      if (match_index >= match_nrecs) {
        break;
      }
      pMatch = &match_table[match_index];
    }
    if (pSextractor->NUMBER == pMatch->NUMBER) {
      /* Copy match data to SExtractor table */
#ifdef NO_DRAD_COLUMN
      pSextractor->drad = 0.0;
      pSextractor->AFLAGS = 0;
#else /* NO_DRAD_COLUMN */
      pSextractor->drad = pMatch->drad;
      pSextractor->AFLAGS = pMatch->AFLAGS;
#endif /* NO_DRAD_COLUMN */
      pSextractor->BFLAGS = pMatch->BFLAGS;
      pSextractor->flag_match = 1;
      pSextractor->x = pSextractor->X_IMAGE;
      pSextractor->y = pSextractor->Y_IMAGE;
      matchVector[matchVectorCount] = sextractor_index;
      matchVectorCount++;
    }
  }



  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];
    pSextractor->internal = 0; 
    /* Set all alases now */

    /* Clear any stale defect filter information */
    pSextractor->AFLAGS &= ~(1<<FILTER_AFLAG_DEFECT);
#ifdef NO_DRAD_FILTER


    if (((pSextractor->AFLAGS & ((1<<FILTER_AFLAG_DRAD)|(1<<FILTER_AFLAG_DRADBIN))) != 0) ||
        ((pSextractor->AFLAGS & FILTER_AMASK_BLEND) != 0) || 
        ((pSextractor->BFLAGS & FILTER_BMASK_BLEND) != 0) ||
				((pSextractor->BFLAGS & (1 << FILTER_BFLAG_PMERROR)) != 0)) {
      pSextractor->internal |= INTERNAL_FLAG_REJECT;
      dradFilterCount++;
    }    
#endif /* NO_DRAD_FILTER */
  }
  /* Script line 54: changed: make the peak of the distribution near the center */
  ixaCount = 0;
  ixbCount = 0;
  for (matchIndex = 0; matchIndex < matchVectorCount; matchIndex++) {
    pSextractor = &sextractor_table[matchVector[matchIndex]];
    if (fabs(pSextractor->THETA_J2000) < 10) {
      ixa[ixaCount] = matchVector[matchIndex];
      ixaCount++;
    }
    if (fabs(pSextractor->THETA_J2000) > 80) {
      ixb[ixbCount] = matchVector[matchIndex];
      ixbCount++;
    }
    pSextractor->THETA = pSextractor->THETA_J2000;
  }
#if 0
  printf("ixaCount %d ixbCount %d\n",ixaCount,ixbCount);
#endif
  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];
    pSextractor->THETA3 = pSextractor->THETA_J20003;
#ifdef USE_XY_RANGE
    if (numpoints == 0) {
      minimumX = pSextractor->X_IMAGE;
      maximumX = pSextractor->X_IMAGE;
      minimumY = pSextractor->Y_IMAGE;
      maximumY = pSextractor->Y_IMAGE;
    } else {
      if (minimumX > pSextractor->X_IMAGE) {
        minimumX = pSextractor->X_IMAGE;
      }
      if (maximumX < pSextractor->X_IMAGE) {
        maximumX = pSextractor->X_IMAGE;
      }
      if (minimumY > pSextractor->Y_IMAGE) {
        minimumY = pSextractor->Y_IMAGE;
      }
      if (maximumY < pSextractor->Y_IMAGE) {
        maximumY = pSextractor->Y_IMAGE;
      }

    }

#endif /* USE_XY_RANGE */
    numpoints++;
  }
  if (ixaCount < ixbCount) {
    for (matchIndex = 0; matchIndex < matchVectorCount; matchIndex++) {
      pSextractor = &sextractor_table[matchVector[matchIndex]];
      if (pSextractor->THETA_J2000 < 0) {
        pSextractor->THETA = pSextractor->THETA_J2000+180.0;
      }
    }
    for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
      pSextractor = &sextractor_table[sextractor_index];
      if (pSextractor->THETA_J2000 < 0) {
        pSextractor->THETA = pSextractor->THETA_J2000+180.0;
      }
      if (pSextractor->THETA_J20003 < 0) {
        pSextractor->THETA3 = pSextractor->THETA_J20003 + 180.0;
      }
      pSextractor->theta = pSextractor->THETA;  /* THETA_J2000 */
      pSextractor->thetaf1 = pSextractor->THETA;

    }
   
  }


  /* Script line 70: % to avoid crash for some problematic objects (already done earlier!) */

  for (matchIndex = 0; matchIndex < matchVectorCount; matchIndex++) {
    pSextractor = &sextractor_table[matchVector[matchIndex]];
    if ((pSextractor->internal & INTERNAL_FLAG_REJECT) != 0) {
      continue;
    }
      
  }
  /* Line 78: Now look for objects of high and low drad */
  
  index = 0;
  for (matchIndex = 0; matchIndex < matchVectorCount; matchIndex++) {
    pSextractor = &sextractor_table[matchVector[matchIndex]];
    if ((pSextractor->internal & (INTERNAL_FLAG_FLUX|INTERNAL_FLAG_REJECT)) != 0) {
      continue;
    }
    if (pSextractor->drad < 90.0) {
      vector[index] = pSextractor->drad;
      index++;
    }
  }

#ifndef USE_XY_RANGE
  minimumX = 0;
  maximumX = 1.0 * mosaicWidth;
  minimumY = 0;
  maximumY = 1.0 * mosaicHeight;
#endif /* USE_XY_RANGE */
  xrange = maximumX-minimumX;
  yrange = maximumY-minimumY;
  if (verbose) {
    printf("%d %d X: Min %f Max %f range %f  Y: Min %f Max %f range %f from numpoints %d maximum plate distance %f\n",minimumX,maximumX,xrange,minimumY,maximumY,yrange,numpoints,max_plate_dist);
  }
#ifdef NO_DRAD_FILTER
  meddrad1 = 99.0;
#else /* NO_DRAD_FILTER */
  CalcMedianAndRMS(index,0,vector,&meddrad1,&stddrad1,0,1.0,0);
#endif /* NO_DRAD_FILTER */
#if 0
  if (verbose) {
    printf("drad median is %f and rms is %f for %d stars\n",meddrad1,stddrad1,index);
  }
#endif


  for (matchIndex = 0; matchIndex < matchVectorCount; matchIndex++) {
    pSextractor = &sextractor_table[matchVector[matchIndex]];
    if ((pSextractor->internal & (INTERNAL_FLAG_FLUX|INTERNAL_FLAG_REJECT)) != 0) {
      continue;
    }
    if (pSextractor->drad <= (meddrad1+(2.0*stddrad1))) {
      kxr++;
    }
  

  }
  /* Line 85: % devided into 5x5 bins (ensure at least 200 points per bin), then do the local filtering  */
  if (verbose) {
    printf("kxr is %d\n",kxr);
  }

  if (kxr > 10) { /* This conditional omitted in version of 6/28/2012 */
    nbin = round(sqrt((1.0*kxr)/1000.0));
    if (nbin > 10) {
      nbin = 10;
    }
    if ((nbin+1) > MAX_LEVEL_VECTOR) {
      printf("ERROR: nbin %d greater than MAX_LEVEL_VECTOR %d\n",nbin,MAX_LEVEL_VECTOR);
      exit(-1);
    }

    for (iii = 0; iii < nbin; iii++) { 
      for (jjj = 0; jjj < nbin; jjj++) { 
        /* for (iii = 7; iii < 8; iii++) {   
           for (jjj = 9; jjj < 10; jjj++) {  */


        memset(maggrid,0,sizeof(maggrid));
        memset(iso0grid,0,sizeof(iso0grid));
        memset(fluxmaxgrid,0,sizeof(fluxmaxgrid));
        memset(fluxmaxc,0,sizeof(fluxmaxc));
        memset(rmsfluxmaxc,0,sizeof(rmsfluxmaxc));
        memset(iso0c,0,sizeof(iso0c));
        memset(rmsiso0c,0,sizeof(rmsiso0c));
        memset(fwhmc,0,sizeof(fwhmc));
        memset(rmsfwhmc,0,sizeof(rmsfwhmc));
        memset(ec,0,sizeof(ec));
        memset(rmsec,0,sizeof(rmsec));

        ixbinCount = 0;
        ixbin2Count = 0;
        ixbinfCount = 0;
        ixbinrCount = 0;
        ixbinf1Count = 0;
        ixCount = 0;
        ix0Count = 0;
        if (doPlots) {
          strcpy(workingDirectory,outfile);
          charPtr = strrchr(workingDirectory,'/');
          if (charPtr != NULL) {
            charPtr++;
            *charPtr = 0;
          } else {
            workingDirectory[0] = 0;
          }
          sprintf(plotdev,"%s%s_%d_%d_plot.eps/cps",workingDirectory,fileroot,iii+1,jjj+1);
          cpgopen(plotdev);
          cpgslw(2); /* line width is 2/200 inch */
          cpgsch(2.0); /* Character size is 2.0/40 the height of the view surface */
          /* subdivide the page into 4 panels;*/
          cpgsubp(2,2);
        }

        xmintest = (minimumX+(0.05*xrange)) + (((1.0 * iii)    *(0.9*xrange)/(1.0*nbin)));
        xmaxtest = (minimumX+(0.05*xrange)) + (((1.0 * (iii+1))*(0.9*xrange)/(1.0*nbin)));
        ymintest = (minimumY+(0.05*yrange)) + (((1.0 * jjj)    *(0.9*yrange)/(1.0*nbin)));
        ymaxtest = (minimumY+(0.05*yrange)) + (((1.0 * (jjj+1))*(0.9*yrange)/(1.0*nbin)));
        if (iii == 0) {
          xmintest = minimumX;
        }
        if (iii == (nbin-1)) {
          xmaxtest = maximumX+1;
        }
        if (jjj == 0) {
          ymintest = minimumY;
        }
        if (jjj == (nbin-1)) {
          ymaxtest = maximumY+1; 
        }
#if 0
        if (verbose) {
          printf("Bin %2d %2d, Xrange %f  %f  xrange %f Yrange %f %f  yrange %f\n",iii+1,jjj+1,
                 xmintest,
                 xmaxtest,
                 xrange,
                 ymintest,
                 ymaxtest,
                 yrange);
        }
#endif
        /* Line 106 % all objects */




        for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
          pSextractor = &sextractor_table[sextractor_index];

          if ((pSextractor->x3 >= xmintest) &&
              (pSextractor->x3 <  xmaxtest) &&
              (pSextractor->y3 >= ymintest) &&
              (pSextractor->y3 < ymaxtest)) {
#ifdef LOS_DEBUG_NUMBER
            if (pSextractor == pDebugSextractor) {
              if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
                printf("At NUMBER %d in line %d iii %d jjj %d\n",pSextractor->NUMBER,__LINE__,iii+1,jjj+1);
              }
             
            }
#endif /* LOS_DEBUG_NUMBER */
            ixbin[ixbinCount] = sextractor_index;
            ixbinCount++;
            if (ixbinCount == 0) {
              minmag = pSextractor->mag;
              maxmag = pSextractor->mag;
              minfluxmax = pSextractor->fluxmax;
              maxfluxmax = pSextractor->fluxmax;
            } else {
              if (pSextractor->mag < minmag) {
                minmag = pSextractor->mag;
              }
              if (pSextractor->mag > maxmag) {
                maxmag = pSextractor->mag;
              }
              if (pSextractor->fluxmax < minfluxmax) {
                minfluxmax = pSextractor->fluxmax;
              }
              if (pSextractor->fluxmax > maxfluxmax) {
                maxfluxmax = pSextractor->fluxmax;
              }
             
            }


          }
        }



        for (matchIndex = 0; matchIndex < matchVectorCount; matchIndex++) {


          pSextractor = &sextractor_table[matchVector[matchIndex]];
          /* Line 95: larger bins for the edge bins; equal effective area for each bin */

          /* INTERNAL_FLAG_REJECT test removed from here on Sep 26, 2008 and moved down seven lines */

          if ((pSextractor->x >= xmintest) &&
              (pSextractor->x <  xmaxtest) &&
              (pSextractor->y >= ymintest) &&
              (pSextractor->y < ymaxtest)) {
            if ((pSextractor->internal & INTERNAL_FLAG_REJECT) != 0) {
              continue;
            }

            if ((pSextractor->internal & INTERNAL_FLAG_FLUX) != 0) {
              continue;
            }
            if ((pSextractor->drad < 90.0) &&
                (pSextractor->x > (minimumX + (0.05*(1.0*xrange)))) &&
                (pSextractor->x < (maximumX - (0.05*(1.0*xrange)))) &&
                (pSextractor->y > (minimumY + (0.05*(1.0*yrange)))) &&
                (pSextractor->y < (maximumY - (0.05*(1.0*yrange)))) &&
                (pSextractor->plate_dist < (0.85*max_plate_dist))) {
              ixbin2[ixbin2Count] =  matchVector[matchIndex];
              vector[ixbin2Count] = pSextractor->drad;
              ixbin2Count++;
            }
          }
        }
        if (ixbin2Count <= 0) { 
          cpgend();
          continue;
        }
#ifdef NO_DRAD_FILTER
        meddrad2 = 99.0;
#else /* NO_DRAD_FILTER */
        CalcMedianAndRMS(ixbin2Count,0,vector,&meddrad2,&stddrad2,0,1.0,0);
#endif /* NO_DRAD_FILTER */
        /* Line 108: cut off 10% edge for real stars */
        ixbinf1Count = 0;
        ixbinrCount = 0;
        for (index = 0; index < ixbin2Count; index++) {
          pSextractor = &sextractor_table[ixbin2[index]];
          if (pSextractor->drad <= (2*meddrad2)) {
            ixbinr[ixbinrCount] = ixbin2[index];
            if (ixbinrCount == 0) {
              minmagr = pSextractor->mag;
              maxmagr = pSextractor->mag;
              minfluxmaxr = pSextractor->fluxmax;
              maxfluxmaxr = pSextractor->fluxmax;
            } else {
              if (pSextractor->mag < minmagr) {
                minmagr = pSextractor->mag;
              }
              if (pSextractor->mag > maxmagr) {
                maxmagr = pSextractor->mag;
              }
              if (pSextractor->fluxmax < minfluxmaxr) {
                minfluxmaxr = pSextractor->fluxmax;
              }
              if (pSextractor->fluxmax > maxfluxmaxr) {
                maxfluxmaxr = pSextractor->fluxmax;
              }
             

            }
            ixbinrCount++;
          }           
          /* Line 111:  cut off 20% edge for dubious stars */
          if (pSextractor->drad > (4*meddrad2)) {
            if ((pSextractor->x > (minimumX + (0.1*(1.0*xrange)))) &&
                (pSextractor->x < (maximumX - (0.1*(1.0*xrange)))) &&
                (pSextractor->y > (minimumY + (0.1*(1.0*yrange)))) &&
                (pSextractor->y < (maximumY - (0.1*(1.0*yrange)))) &&
                (pSextractor->plate_dist < (0.70*max_plate_dist))) {
              ixbinf1[ixbinf1Count] = ixbin2[index];
              ixbinf1Count++;
            }


          }
          
        }
        if (ixbinrCount > 50) {
          /* line 115: % real objects */
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            pSextractor->magr = pSextractor->MAG_ISO;
            pSextractor->fwhmr = checklog10(3600.0 * pSextractor->FWHM_WORLD);
            pSextractor->thetar = pSextractor->THETA;
            pSextractor->er = pSextractor->ELLIPTICITY;
            pSextractor->fluxmaxr = checklog10(pSextractor->FLUX_MAX);
            pSextractor->iso0r = checklog10(1.0*pSextractor->ISO0);
            vector[index] = pSextractor->magr;
          }
          /* line 123: %fake objects as 4-sigma object */
          for (index = 0; index < ixbinf1Count; index++) {
            pSextractor = &sextractor_table[ixbinf1[index]];
            pSextractor->magf1 = pSextractor->MAG_ISO;
            pSextractor->fwhmf1 = checklog10(3600.0 * pSextractor->FWHM_WORLD);
            pSextractor->thetaf1 = pSextractor->THETA;
            pSextractor->ef1 = pSextractor->ELLIPTICITY;
            pSextractor->fluxmaxf1 = checklog10(pSextractor->FLUX_MAX);
            pSextractor->iso0f1 = checklog10(1.0*pSextractor->ISO0);
         
          }
          /* Line 131: % all objects */
          for (index = 0; index < ixbinCount; index++) {
            pSextractor = &sextractor_table[ixbin[index]];
            pSextractor->num = pSextractor->NUMBER3;
            pSextractor->mag = pSextractor->MAG_ISO3;
            pSextractor->fwhm = checklog10(3600.0 * pSextractor->FWHM_WORLD3);
            pSextractor->theta = pSextractor->THETA3;
            pSextractor->e = pSextractor->ELLIPTICITY3;
#ifdef LOS_DEBUG_NUMBER
            if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
              printf("At NUMBER %d in line %d\n",pSextractor->NUMBER,__LINE__);
            }

#endif /* LOS_DEBUG_NUMBER */
            pSextractor->fluxmax = checklog10(pSextractor->FLUX_MAX3);
            pSextractor->iso0 = checklog10(1.0*pSextractor->ISO03);
          }

          minmagr += 1.0;


          /* Line 147 % change the binning to prctile; each bin contains similar number of stars */

          nmag = round((1.0 * ixbinrCount)/30.0);
          if (nmag > MAX_MAG_LEVELS) {
            nmag = MAX_MAG_LEVELS;
          }
#if 0
          if (verbose) {
            printf("Bin %d %d, ixbincount %5d, ixbin2count %5d, ixbinrcount %5d ixbinfcount %5d, nmag %3d\n",
                   iii+1,jjj+1,ixbinCount,ixbin2Count,ixbinrCount,ixbinfCount,nmag);
          
          }
#endif
          if (nmag > 0) {
            prctile(vector,ixbinrCount,nmag,maggrid,&maggridsize,0,0);
            if (maggridsize > MAX_LEVEL_VECTOR) {
              printf("ERROR: maggridsize %d is too large %d\n",maggridsize,MAX_LEVEL_VECTOR);
            }

#if 0 
            if ((verbose) && (iii == 2) && (jjj == 2)) {
              printf("At bin %d %d\n",iii+1,jjj+1);
            }
#endif
            for (indexMag = 0; indexMag < (maggridsize-1); indexMag++) {  

              /* for (indexMag = 5; indexMag < 6; indexMag++) {  */
              ixf1temp = 0;
              ixf2temp = 0;
              ixCount = 0;

              for (index = 0; index < ixbinrCount; index++) {
                pSextractor = &sextractor_table[ixbinr[index]];
                if ((pSextractor->magr >= maggrid[indexMag]) &&
                    (pSextractor->magr <  maggrid[indexMag+1])) {
                  ix[ixCount] = ixbinr[index];
                  vector[ixCount] = pSextractor->fwhm;
#if 0
                  printf("ixCount %d NUMBER %d fwhm %f\n",ixCount,pSextractor->NUMBER,pSextractor->fwhm);
#endif
                  ixCount++;
                }
              }
              maggridOK[indexMag] = ixCount;
              if (maggridOK[indexMag] > 0) {
                /* line 158: Test 1: the full width half maximum must be near the mean */
                fwhmc0 = FindMean(vector,ixCount);
                rmsfwhmc0 = FindStd(vector,ixCount);
                if (rmsfwhmc0 < 0.2) {
                  rmsfwhmc0 = 0.2;
                }
                /* Now three sigma clip */
                ix0Count = 0;
                for (index = 0; index < ixCount; index++) {
                  pSextractor = &sextractor_table[ix[index]];
                  if (fabs(pSextractor->fwhm - fwhmc0) <= (3.0*rmsfwhmc0)) {
                    vector[ix0Count] = pSextractor->fwhm;
                    ix0[ix0Count] = ix[index];
                    ix0Count++;
                  }
                }
                maggridOK[indexMag] = ix0Count;
                if (maggridOK[indexMag] > 0) {
                  rmsfwhmc0 = FindStd(vector,ix0Count);
                  ix0Count = 0;
                  for (index = 0; index < ixCount; index++) {
                    pSextractor = &sextractor_table[ix[index]];
                    if (fabs(pSextractor->fwhm - fwhmc0) <= (3.0*rmsfwhmc0)) {
                      vector[ix0Count] = pSextractor->fwhm;
                      ix0[ix0Count] = ix[index];
                      ix0Count++;
                    }
                  }
                }
                maggridOK[indexMag] = ix0Count;
                if (maggridOK[indexMag] > 0) {
                  rmsfwhmc0 = FindStd(vector,ix0Count);
                  ix0Count = 0;
                  for (index = 0; index < ixCount; index++) {
                    pSextractor = &sextractor_table[ix[index]];
                    if (fabs(pSextractor->fwhm - fwhmc0) <= (3.0*rmsfwhmc0)) {
                      vector[ix0Count] = pSextractor->fwhm;
                      ix0[ix0Count] = ix[index];
                      ix0Count++;
                    }
                  }
                }
                maggridOK[indexMag] = ix0Count;
                if (maggridOK[indexMag] > 0) {
                  fwhmc[indexMag] = FindMean(vector,ix0Count);
                  rmsfwhmc[indexMag] = FindStd(vector,ix0Count);
#if 0
                  printf("Bin %d %d %d fwhmc %f rmsfwhmc %f ix0 %d\n",iii+1,jjj+1,indexMag+1,fwhmc[indexMag],rmsfwhmc[indexMag],ix0Count);
#endif

                  for (index = 0; index < ixbinCount; index++) {
                    pSextractor = &sextractor_table[ixbin[index]];
                    if ((pSextractor->mag >= maggrid[indexMag]) &&
                        (pSextractor->mag <  maggrid[indexMag+1]) &&
                        (fabs(pSextractor->fwhm - fwhmc[indexMag]) < (3.0 * rmsfwhmc[indexMag]))) {
                      pSextractor->internal |= INTERNAL_FLAG_IXF1;
                      ixf1Count++;
                      ixf1temp++;
                    }
                  }
                }
#if 0
                if ((verbose) && (iii == 2) && (jjj == 2)) {
                  printf("Bin %d %d, i %2d fwhm %f rmsfwhmc %f, length(fwhm) %d length ixf1 %d\n",iii+1,jjj+1,indexMag,fwhmc[indexMag],rmsfwhmc[indexMag],ixbinCount,ixf1Count-ixf1old2);
                }
#endif
                ixf1old2 = ixf1Count;

                /* line 170 Test 2: The eccentricity must be near the mean */
                for (index = 0; index < ixCount; index++) {
                  pSextractor = &sextractor_table[ix[index]];
                  vector[index] = pSextractor->e;
                }
                ec0 = FindMean(vector,ixCount);
                rmsec0 = FindStd(vector,ixCount);
                /* Now three sigma clip */
                ix0Count = 0;
                for (index = 0; index < ixCount; index++) {
                  pSextractor = &sextractor_table[ix[index]];
                  if (fabs(pSextractor->e - ec0) <= (3.0*rmsec0)) {
                    vector[ix0Count] = pSextractor->e;
                    ix0[ix0Count] = ix[index];
                    ix0Count++;
                  }
                }
                rmsec0 = FindStd(vector,ix0Count);
                ix0Count = 0;
                for (index = 0; index < ixCount; index++) {
                  pSextractor = &sextractor_table[ix[index]];
                  if (fabs(pSextractor->e - ec0) <= (3.0*rmsec0)) {
                    vector[ix0Count] = pSextractor->e;
                    ix0[ix0Count] = ix[index];
                    ix0Count++;
                  }
                }
                rmsec0 = FindStd(vector,ix0Count);
                ix0Count = 0;
                for (index = 0; index < ixCount; index++) {
                  pSextractor = &sextractor_table[ix[index]];
                  if (fabs(pSextractor->e - ec0) <= (3.0*rmsec0)) {
                    vector[ix0Count] = pSextractor->e;
                    ix0[ix0Count] = ix[index];
                    ix0Count++;
                  }
                }
                ec[indexMag] = FindMean(vector,ix0Count);
                rmsec[indexMag] = FindStd(vector,ix0Count);
#if 0
                printf("Bin %d %d %d ec %f rmsec %f ix0 %d\n",iii+1,jjj+1,indexMag+1,ec[indexMag],rmsec[indexMag],ix0Count);
#endif
                for (index = 0; index < ixbinCount; index++) {
                  pSextractor = &sextractor_table[ixbin[index]];
                  if ((pSextractor->mag >= maggrid[indexMag]) &&
                      (pSextractor->mag <  maggrid[indexMag+1]) &&
                      (fabs(pSextractor->e - ec[indexMag]) < (3.0 * rmsec[indexMag]))) {
                    pSextractor->internal |= INTERNAL_FLAG_IXF2;
                    ixf2Count++;
                    ixf2temp++;
                  }
                }
#if 0
                if ((verbose) && (iii == 2) && (jjj == 2)) {
                  printf("Bin %d %d, i %2d ec[indexMag] %f rmsec[indexMag] %f, length(fwhm) %d length ixf2 %d\n",iii+1,jjj+1,indexMag,ec[indexMag],rmsec[indexMag],ixbinCount,ixf2Count-ixf2old2);
                }
#endif
                ixf2old2 = ixf2Count;
              }
#if 0
              printf("bin %d ixf1temp %d ixf2temp %d\n",indexMag+1,ixf1temp,ixf2temp);
#endif
            } /* Loop on magnitudes */

            if (doPlots) {
              /* cpgpanl(1,1); */
              plot_nrecs1 = ixbinrCount;
#ifdef LOS_DEBUG_NUMBER
                xval1Debug = pDebugSextractor->magr;
                yval1Debug = pDebugSextractor->fwhmr;
#endif  /* LOS_DEBUG_NUMBER */
              for (plot_index1 = 0; plot_index1 < plot_nrecs1; plot_index1++) {
                pSextractor = &sextractor_table[ixbinr[plot_index1]];
                xval1[plot_index1] = pSextractor->magr;
                yval1[plot_index1] = pSextractor->fwhmr;
                if (plot_index1 == 0) {
                  xmin = xval1[0];
                  xmax = xval1[0];
                  ymin = yval1[0];
                  ymax = yval1[0];
                } else {
                  if (xmin > xval1[plot_index1]) {
                    xmin = xval1[plot_index1];
                  }
                  if (xmax < xval1[plot_index1]) {
                    xmax = xval1[plot_index1];
                  }
                  if (ymin > yval1[plot_index1]) {
                    ymin = yval1[plot_index1];
                  }
                  if (ymax < yval1[plot_index1]) {
                    ymax = yval1[plot_index1];
                  }
 
                }
              }
              plot_nrecs2 = ixbinf1Count;
              for (plot_index2 = 0; plot_index2 < plot_nrecs2; plot_index2++) {
                pSextractor = &sextractor_table[ixbinf1[plot_index2]];
                xval2[plot_index2] = pSextractor->magf1;
                yval2[plot_index2] = pSextractor->fwhmf1;
                if (xmin > xval2[plot_index2]) {
                  xmin = xval2[plot_index2];
                }
                if (xmax < xval2[plot_index2]) {
                  xmax = xval2[plot_index2];
                }
                if (ymin > yval2[plot_index2]) {
                  ymin = yval2[plot_index2];
                }
                if (ymax < yval2[plot_index2]) {
                  ymax = yval2[plot_index2];
                }
              }
              for (indexMag = 0; indexMag < (maggridsize-1); indexMag++) {
                if (maggridOK[indexMag] > 0) {
                  valx = maggrid[indexMag];
                  valy = 3*rmsfwhmc[indexMag]+fwhmc[indexMag];
                  if (xmin > valx) {
                    xmin = valx;
                  }
                  if (xmax < valx) {
                    xmax = valx;
                  }
                  if (ymin > valy) {
                    ymin = valy;
                  }
                  if (ymax < valy) {
                    ymax = valy;
                  }
                  valy = -3*rmsfwhmc[indexMag]+fwhmc[indexMag];
                  if (ymin > valy) {
                    ymin = valy;
                  }
                  if (ymax < valy) {
                    ymax = valy;
                  }
                }
              }
              cpgenv(xmin-0.1,xmax+0.1,ymin-0.1,ymax+0.1,0,0);   
              sprintf(title,"%s [%d,%d]",fileroot,iii+1,jjj+1);
              cpglab("MAG_ISO","log(FWHM) (arcsec)",title);
              cpgsch(1.0); /* Character size is 1.0/40 the height of the view surface */
              cpgsci(4);  /* Set blue */
              cpgpt(plot_nrecs1,xval1,yval1,16);
              cpgsci(2); /* Set red */
#ifdef LOS_DEBUG_NUMBER
              if ((xval1Debug != 0) && (yval1Debug != 0)) {
                cpgpt(1,&xval1Debug,&yval1Debug,19);
                xval1Debug = 0.0;
                yval1Debug = 0.0;
              }
#endif /* LOS_DEBUG_NUMBER */

              cpgpt(plot_nrecs2,xval2,yval2,4);
              cpgsch(2.0); /* Character size is 2.0/40 the height of the view surface */

              cpgsci(1); /* Set black */
              for (indexMag = 0; indexMag < maggridsize-1; indexMag++) {
                if (maggridOK[indexMag] > 0) {
                  xval1[0] = maggrid[indexMag];
                  xval1[1] = maggrid[indexMag+1];
                  yval1[0] = 3*rmsfwhmc[indexMag]+fwhmc[indexMag];
                  yval1[1] = 3*rmsfwhmc[indexMag]+fwhmc[indexMag];
                  cpgline(2,xval1,yval1);
                  yval1[0] = -3*rmsfwhmc[indexMag]+fwhmc[indexMag];
                  yval1[1] = -3*rmsfwhmc[indexMag]+fwhmc[indexMag];
                  cpgline(2,xval1,yval1);
                }

              }

              cpgpanl(1,1);
              plot_nrecs1 = ixbinrCount;
#ifdef LOS_DEBUG_NUMBER
                xval1Debug = pDebugSextractor->magr;
                yval1Debug = pDebugSextractor->er;
#endif  /* LOS_DEBUG_NUMBER */

              for (plot_index1 = 0; plot_index1 < plot_nrecs1; plot_index1++) {
                pSextractor = &sextractor_table[ixbinr[plot_index1]];
                xval1[plot_index1] = pSextractor->magr;
                yval1[plot_index1] = pSextractor->er;
                if (plot_index1 == 0) {
                  xmin = xval1[0];
                  xmax = xval1[0];
                  ymin = yval1[0];
                  ymax = yval1[0];
                } else {
                  if (xmin > xval1[plot_index1]) {
                    xmin = xval1[plot_index1];
                  }
                  if (xmax < xval1[plot_index1]) {
                    xmax = xval1[plot_index1];
                  }
                  if (ymin > yval1[plot_index1]) {
                    ymin = yval1[plot_index1];
                  }
                  if (ymax < yval1[plot_index1]) {
                    ymax = yval1[plot_index1];
                  }
 
                }
              }
              plot_nrecs2 = ixbinf1Count;
              for (plot_index2 = 0; plot_index2 < plot_nrecs2; plot_index2++) {
                pSextractor = &sextractor_table[ixbinf1[plot_index2]];
                xval2[plot_index2] = pSextractor->magf1;
                yval2[plot_index2] = pSextractor->ef1;
                if (xmin > xval2[plot_index2]) {
                  xmin = xval2[plot_index2];
                }
                if (xmax < xval2[plot_index2]) {
                  xmax = xval2[plot_index2];
                }
                if (ymin > yval2[plot_index2]) {
                  ymin = yval2[plot_index2];
                }
                if (ymax < yval2[plot_index2]) {
                  ymax = yval2[plot_index2];
                }
              }
              for (indexMag = 0; indexMag < maggridsize; indexMag++) {
                valx = maggrid[indexMag];
                valy = 3*rmsec[indexMag]+ec[indexMag];
                if (xmin > valx) {
                  xmin = valx;
                }
                if (xmax < valx) {
                  xmax = valx;
                }
                if (ymin > valy) {
                  ymin = valy;
                }
                if (ymax < valy) {
                  ymax = valy;
                }
                valy = -3*rmsec[indexMag]+ec[indexMag];
                if (ymin > valy) {
                  ymin = valy;
                }
                if (ymax < valy) {
                  ymax = valy;
                }

              }
              cpgenv(xmin-0.1,xmax+0.1,ymin-0.1,ymax+0.1,0,0);   
              cpglab("MAG_ISO","ELLIPTICITY","");
              cpgsch(1.0); /* Character size is 1.0/40 the height of the view surface */
              cpgsci(4);  /* Set blue */
              cpgpt(plot_nrecs1,xval1,yval1,16);
              cpgsci(2); /* Set red */
#ifdef LOS_DEBUG_NUMBER
              if ((xval1Debug != 0) && (yval1Debug != 0)) {
                cpgpt(1,&xval1Debug,&yval1Debug,19);
                xval1Debug = 0.0;
                yval1Debug = 0.0;
              }
#endif /* LOS_DEBUG_NUMBER */
              cpgpt(plot_nrecs2,xval2,yval2,4);
              cpgsch(2.0); /* Character size is 2.0/40 the height of the view surface */

              cpgsci(1); /* Set black */
              for (indexMag = 0; indexMag < maggridsize-1; indexMag++) {
                xval1[0] = maggrid[indexMag];
                xval1[1] = maggrid[indexMag+1];
                yval1[0] = 3*rmsec[indexMag]+ec[indexMag];
                yval1[1] = 3*rmsec[indexMag]+ec[indexMag];
                cpgline(2,xval1,yval1);
                yval1[0] = -3*rmsec[indexMag]+ec[indexMag];
                yval1[1] = -3*rmsec[indexMag]+ec[indexMag];
                cpgline(2,xval1,yval1);
                

              }


            } /* End of plot */


          } /* nmag check */
            /* Line 217: filter 3: 3-sigma clipping on theta, 2.5-sigma on ellipticity */
        
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
 

            vector[index] = pSextractor->thetar;
          }
          CalcMedianAndRMS(ixbinrCount,0,vector,&thetac0,&stdthetac0,0,1.0,0);
          for (index = 0; index < ixbinCount; index++) {
            pSextractor = &sextractor_table[ixbin[index]];
            pSextractor->thetadist = fabs(pSextractor->theta-thetac0);
            if (pSextractor->thetadist > 90) {
              pSextractor->thetadist = 180-pSextractor->thetadist;
            }
          }
          rmsthetac0 = 0;
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            pSextractor->thetadistr = fabs(pSextractor->thetar-thetac0);
            if ( pSextractor->thetadistr > 90) {
              pSextractor->thetadistr = 180.0-pSextractor->thetadistr;
            }
            rmsthetac0 += sqr(pSextractor->thetadistr);
          }
          /* line 226 % revised std, to reflect the continuous distribution of theta, i.e. 0==180 */
          rmsthetac0 = sqrt(rmsthetac0/(1.0*ixbinrCount));
          ix0Count = 0;
          rmsthetac0Sum = 0;
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            if (pSextractor->thetadistr <= (2*rmsthetac0)) {
              ix0Count++;
              rmsthetac0Sum += sqr(pSextractor->thetadistr);
            }
          }
          rmsthetac0 = sqrt(rmsthetac0Sum/(1.0*ix0Count));
          ix0Count = 0;
          rmsthetac0Sum = 0;
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            if (pSextractor->thetadistr <= (2*rmsthetac0)) {
              ix0Count++;
              rmsthetac0Sum += sqr(pSextractor->thetadistr);
            }
          }
          rmsthetac0 = sqrt(rmsthetac0Sum/(1.0*ix0Count));
          ix0Count = 0;
          rmsthetac0Sum = 0;
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            if (pSextractor->thetadistr <= (2*rmsthetac0)) {
              ix0Count++;
              rmsthetac0Sum += sqr(pSextractor->thetadistr);
            }
          }
          rmsthetac = sqrt(rmsthetac0Sum/(1.0*ix0Count));
#if 0
          printf("Bin %d %d thetac0 %f rmsthetac0 %f, length(ix0) %d\n",iii+1,jjj+1,thetac0,rmsthetac,ix0Count);
#endif


          /* Now work with the eccentricity */
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            vector[index] = pSextractor->er;
          } 
          CalcMedianAndRMS(ixbinrCount,0,vector,&ec0,&stdec0,0,1.0,0);
          rmsec0 = FindStd(vector,ixbinrCount);
          if ((verbose) && (iii == 2) && (jjj == 2)) {
#if 0
            for (index = 0; index < ixbinrCount; index++) {
              printf("%f ",vector[index]);
            }
            printf("\n");
#endif
          }


          ix0Count = 0;
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            if (fabs(pSextractor->er - ec0) <= (2.0*rmsec0)) {
              vector[ix0Count] = pSextractor->er;
              ix0[ix0Count] = ixbinr[index];
              ix0Count++;
            }
          }
#if 0
          if ((verbose) && (iii == 2) && (jjj == 2)) {
            for (index = 0; index < ix0Count; index++) {
              printf("%f ",vector[index]);
            }
            printf("\n");
            printf("Bin %d %d ec0 %f rmsec0 %f, ix0count %d\n",iii+1,jjj+1,ec0,rmsec0,ix0Count);
          }
#endif

          rmsec0 = FindStd(vector,ix0Count);
          ix0Count = 0;
  


          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            if (fabs(pSextractor->er - ec0) <= (2.0*rmsec0)) {
              vector[ix0Count] = pSextractor->er;
              ix0[ix0Count] = ixbinr[index];
              ix0Count++;
            }
          }

          ec[0] = FindMean(vector,ix0Count);
          rmsec[0] = FindStd(vector,ix0Count);
#if 0
          printf("Bin %d %d ec %f rmsec %f, length(ix0) %d\n",iii+1,jjj+1,ec[0],rmsec[0],ix0Count);
#endif

          for (index = 0; index < ixbinCount; index++) {
            pSextractor = &sextractor_table[ixbin[index]];
            if ((pSextractor->thetadist < (2.5*rmsthetac)) &&
                (fabs(pSextractor->e - ec[0]) < (2.5*rmsec[0]))) {
              pSextractor->internal |= INTERNAL_FLAG_IXF3;
              ixf3Count++;
            }
          }
          /* Line 247 */
          if (doPlots) {
            cpgpanl(2,1);
            plot_nrecs1 = ixbinrCount;
#ifdef LOS_DEBUG_NUMBER
                xval1Debug = pDebugSextractor->thetar;
                yval1Debug = pDebugSextractor->er;
#endif  /* LOS_DEBUG_NUMBER */

            for (plot_index1 = 0; plot_index1 < plot_nrecs1; plot_index1++) {
              pSextractor = &sextractor_table[ixbinr[plot_index1]];
              xval1[plot_index1] = pSextractor->thetar;
              yval1[plot_index1] = pSextractor->er;
                
              if (plot_index1 == 0) {
                xmin = xval1[0];
                xmax = xval1[0];
                ymin = yval1[0];
                ymax = yval1[0];
              } else {
                if (xmin > xval1[plot_index1]) {
                  xmin = xval1[plot_index1];
                }
                if (xmax < xval1[plot_index1]) {
                  xmax = xval1[plot_index1];
                }
                if (ymin > yval1[plot_index1]) {
                  ymin = yval1[plot_index1];
                }
                if (ymax < yval1[plot_index1]) {
                  ymax = yval1[plot_index1];
                }
 
              }
            }
            minthetar = xmin;
            maxthetar = xmax;
            plot_nrecs2 = ixbinf1Count;
              
            for (plot_index2 = 0; plot_index2 < plot_nrecs2; plot_index2++) {
              pSextractor = &sextractor_table[ixbinf1[plot_index2]];
              xval2[plot_index2] = pSextractor->thetaf1;
              yval2[plot_index2] = pSextractor->ef1;
              if (xmin > xval2[plot_index2]) {
                xmin = xval2[plot_index2];
              }
              if (xmax < xval2[plot_index2]) {
                xmax = xval2[plot_index2];
              }
              if (ymin > yval2[plot_index2]) {
                ymin = yval2[plot_index2];
              }
              if (ymax < yval2[plot_index2]) {
                ymax = yval2[plot_index2];
              }
            }
            valx = thetac0+2.5*rmsthetac;
            if (xmin > valx) {
              xmin = valx;
            }
            if (xmax < valx) {
              xmax = valx;
            }
            valx = thetac0-2.5*rmsthetac;
            if (xmin > valx) {
              xmin = valx;
            }
            if (xmax < valx) {
              xmax = valx;
            }
            valy = rmsec[0]*2.5+ec[0];
            if (ymin > valy) {
              ymin = valy;
            }
            if (ymax < valy) {
              ymax = valy;
            }
            valy = rmsec[0]*2.5-ec[0];
            if (ymin > valy) {
              ymin = valy;
            }
            if (ymax < valy) {
              ymax = valy;
            }
            if (ymin > 0) {
              ymin = 0;
            }
            if (ymax < 1.0) {
              ymax = 1.0;
            }
            
            cpgenv(minthetar-1.0,maxthetar+1.0,ymin,ymax,0,0);   
            cpglab("\\gh","Ellipticity","");
            cpgsch(1.0); /* Character size is 1.0/40 the height of the view surface */
            cpgsci(4);  /* Set blue */
            cpgpt(plot_nrecs1,xval1,yval1,16);
            cpgsci(2); /* Set red */
#ifdef LOS_DEBUG_NUMBER
              if ((xval1Debug != 0) && (yval1Debug != 0)) {
                cpgpt(1,&xval1Debug,&yval1Debug,19);
                xval1Debug = 0.0;
                yval1Debug = 0.0;
              }
#endif /* LOS_DEBUG_NUMBER */

            cpgpt(plot_nrecs2,xval2,yval2,4);
            cpgsch(2.0); /* Character size is 2.0/40 the height of the view surface */
            cpgsls(2); /* set dotted */
            cpgsci(1); /* Set black */
            xval1[0] = minthetar;
            xval1[1] = maxthetar;
            yval1[0] = rmsec[0]*2.5 + ec[0];
            yval1[1] = rmsec[0]*2.5 + ec[0];
            cpgline(2,xval1,yval1);
            yval1[0] = -rmsec[0]*2.5 + ec[0];
            yval1[1] = -rmsec[0]*2.5 + ec[0];
            cpgline(2,xval1,yval1);

            xval1[0] = thetac0+2.5*rmsthetac;
            xval1[1] = thetac0+2.5*rmsthetac;
            yval1[0] = 0;
            yval1[1] = 1;
            cpgline(2,xval1,yval1);
            xval1[0] = thetac0-2.5*rmsthetac;
            xval1[1] = thetac0-2.5*rmsthetac;
            yval1[0] = 0;
            yval1[1] = 1;
            cpgline(2,xval1,yval1);
                
            cpgsls(1); /* set full line */

             
          } /* End of plot */

#if 0
          if ((verbose) && (iii == 2) && (jjj == 2)) {
            printf("Bin %d %d, thetac %f rmsthetac %f, ec %f, rmsec %f, theta length %d e length %d ixf3 length %d\n",iii+1,jjj+1,thetac0,rmsthetac0,ec0,rmsec0,ixbinCount,ixbinCount,ixf3Count-ixf3old2);
          }
#endif
          ixf3old2 = ixf3Count;

          /* Line 260: % filter 4: fluxmax vs iso0, 15 bins in mag */
          /* % binning changed to be that each bin contains similar number of stars */
          niso0 = ((1.0 * ixbinrCount)/30.0) + 0.5;
          if (niso0 > MAX_MAG_LEVELS) {
            niso0 = MAX_MAG_LEVELS;
          }
#if 0
          if (verbose) {
            printf("Bin %d %d, ixbincount %5d, ixbin2count %5d, ixbinrcount %5d ixbinfcount %5d, niso0 %3d\n",
                   iii+1,jjj+1,ixbinCount,ixbin2Count,ixbinrCount,ixbinfCount,niso0);
          
          }
#endif
          for (index = 0; index < ixbinrCount; index++) {
            pSextractor = &sextractor_table[ixbinr[index]];
            vector[index] = pSextractor->iso0r;
          }

          if (niso0 > 0) {
            prctile(vector,ixbinrCount,niso0,iso0grid,&iso0gridsize,1,0);
            if (iso0gridsize > MAX_LEVEL_VECTOR) {
              printf("ERROR: iso0gridsize %d is too large %d\n",iso0gridsize,MAX_LEVEL_VECTOR);
            }



            for (indexFluxMax = 0; indexFluxMax < (iso0gridsize-1); indexFluxMax++) {
              ixf4temp = 0;
              ixCount = 0;
#ifdef LOS_DEBUG
              sprintf(debug_name,"/home/scanner/Pipeline/match/losixf4_%d%d_%dc.tmp",iii+1,jjj+1,indexFluxMax+1);
              debugHandle = fopen(debug_name,"wt");
              if (debugHandle == NULL) {
                printf("ERROR: failed to open %s\n",debug_name);
              }
              fprintf(debugHandle,"NUMBER\tX_IMAGE\tY_IMAGE\tISO0\tcount\n");
              fprintf(debugHandle,"------\t-------\t-------\t----\t-----\n");
#endif /* LOS_DEBUG */
              for (index = 0; index < ixbinrCount; index++) {
                pSextractor = &sextractor_table[ixbinr[index]];
                if ((pSextractor->iso0r >= iso0grid[indexFluxMax]) &&
                    (pSextractor->iso0r <  iso0grid[indexFluxMax+1])) {
                  ix[ixCount] = ixbinr[index];
                  vector[ixCount] = pSextractor->fluxmaxr;
#ifdef LOS_DEBUG
                  fprintf(debugHandle,"%d\t%f\t%f\t%f\t%d\n",
                          pSextractor->NUMBER,
                          pSextractor->X_IMAGE,
                          pSextractor->Y_IMAGE,
                          pSextractor->iso0r,
                          (10000*(iii+1))+(100*(jjj+1))+indexFluxMax+1);
#endif /* LOS_DEBUG */


                  ixCount++;
                }
              }
#ifdef LOS_DEBUG
              fclose(debugHandle);
#endif /* LOS_DEBUG */

              iso0gridOK[indexFluxMax] = ixCount;
              if (iso0gridOK[indexFluxMax] > 0) {
                fluxmaxc0 = FindMean(vector,ixCount);
                rmsfluxmaxc0 = FindStd(vector,ixCount);
                if (rmsfluxmaxc0 < 0.3) {
                  rmsfluxmaxc0 = 0.3;
                }
                /* Now sigma clip */
                ix0Count = 0;
                for (index = 0; index < ixCount; index++) {
                  pSextractor = &sextractor_table[ix[index]];
                  if (fabs(pSextractor->fluxmaxr - fluxmaxc0) <= (2.0*rmsfluxmaxc0)) {
                    vector[ix0Count] = pSextractor->fluxmaxr;
                    ix0[ix0Count] = ix[index];
                    ix0Count++;
                  }
                }
                iso0gridOK[indexFluxMax] = ix0Count;
                if (iso0gridOK[indexFluxMax] > 0) {
              
                  rmsfluxmaxc0 = FindStd(vector,ix0Count);
                  ix0Count = 0;
                  for (index = 0; index < ixCount; index++) {
                    pSextractor = &sextractor_table[ix[index]];
                    if (fabs(pSextractor->fluxmaxr - fluxmaxc0) <= (2.0*rmsfluxmaxc0)) {
                      vector[ix0Count] = pSextractor->fluxmaxr;
                      ix0[ix0Count] = ix[index];
                      ix0Count++;
                    }
                  }
                }
                iso0gridOK[indexFluxMax] = ix0Count;
                if (iso0gridOK[indexFluxMax] > 0) {
                  rmsfluxmaxc0 = FindStd(vector,ix0Count);
                  ix0Count = 0;
                  for (index = 0; index < ixCount; index++) {
                    pSextractor = &sextractor_table[ix[index]];
                    if (fabs(pSextractor->fluxmaxr - fluxmaxc0) <= (3.0*rmsfluxmaxc0)) {
                      vector[ix0Count] = pSextractor->fluxmaxr;
                      ix0[ix0Count] = ix[index];
                      ix0Count++;
                    }
                  }
                }
                iso0gridOK[indexFluxMax] = ix0Count;
                if (iso0gridOK[indexFluxMax] > 0) {
                  fluxmaxc[indexFluxMax] = FindMean(vector,ix0Count);
                  rmsfluxmaxc[indexFluxMax] = FindStd(vector,ix0Count);
                  if (indexFluxMax > 0) {
                    if (rmsfluxmaxc[indexFluxMax] > fabs(fluxmaxc[indexFluxMax]-fluxmaxc[indexFluxMax-1])+rmsfluxmaxc[indexFluxMax-1]) {
                      rmsfluxmaxc[indexFluxMax] = (fabs(fluxmaxc[indexFluxMax]-fluxmaxc[indexFluxMax-1])+rmsfluxmaxc[indexFluxMax-1]);
                    }
                  }
                  
#if 0
                  printf("Bin %d %d %d fluxmaxc %f rmsfluxmaxc %f ix0 %d\n",iii+1,jjj+1,indexFluxMax+1,fluxmaxc[indexFluxMax],rmsfluxmaxc[indexFluxMax],ix0Count);
#endif

                  for (index = 0; index < ixbinCount; index++) {
                    pSextractor = &sextractor_table[ixbin[index]];
                    if ((pSextractor->iso0 >= iso0grid[indexFluxMax]) &&
                        (pSextractor->iso0 <  iso0grid[indexFluxMax+1])) {
                      if (fabs(pSextractor->fluxmax - fluxmaxc[indexFluxMax]) < (3.0 * rmsfluxmaxc[indexFluxMax])) {
                        pSextractor->internal |= INTERNAL_FLAG_IXF4;
                        ixf4Count++;
                        ixf4temp++;
                      }
                    }
                  }
                }
              }
#if 0
              if ((verbose) && (iii == 2) && (jjj == 2)) {
                printf("Bin %d %d, i %2d iso0 %f rmsiso0c %f, length(iso0) %d length ixf4 %d\n",iii+1,jjj+1,indexFluxMax,fluxmaxc[indexFluxMax],rmsfluxmaxc[indexFluxMax],ixbinCount,ixf4Count-ixf4old2);
              }
#endif
              ixf4old2 = ixf4Count;
#if 0
              printf("Index %d ixf4temp %d\n",indexFluxMax,ixf4temp);
#endif
            }
          } /* niso0 test */
          
            /* Line 289 bin in FLUX_MAX */

          niso0 = ((1.0 * ixbinrCount)/30.0) + 0.5;
          if (niso0 > MAX_MAG_LEVELS) {
            niso0 = MAX_MAG_LEVELS;
          }
#if 0
          if (verbose) {
            printf("Bin %d %d, ixbincount %5d, ixbin2count %5d, ixbinrcount %5d ixbinfcount %5d, niso0 %3d\n",
                   iii+1,jjj+1,ixbinCount,ixbin2Count,ixbinrCount,ixbinfCount,niso0);
          
          }
#endif
          if (niso0 > 0) {
            for (index = 0; index < ixbinrCount; index++) {
              pSextractor = &sextractor_table[ixbinr[index]];
              vector[index] = pSextractor->fluxmaxr;
            }
            prctile(vector,ixbinrCount,niso0,fluxmaxgrid,&fluxmaxgridsize,1,0);
            if (fluxmaxgridsize > MAX_LEVEL_VECTOR) {
              printf("ERROR: fluxmaxgridsize %d is too large %d\n",fluxmaxgridsize,MAX_LEVEL_VECTOR);
            }



            for (indexFluxMax = 0; indexFluxMax < (fluxmaxgridsize-1); indexFluxMax++) {



              ixCount = 0;
              for (index = 0; index < ixbinrCount; index++) {
                pSextractor = &sextractor_table[ixbinr[index]];
                if ((pSextractor->fluxmaxr >= fluxmaxgrid[indexFluxMax]) &&
                    (pSextractor->fluxmaxr <  fluxmaxgrid[indexFluxMax+1])) {
                  ix[ixCount] = ixbinr[index];
                  vector[ixCount] = pSextractor->iso0r;


                  ixCount++;
                }
              }

              iso0cOK[indexFluxMax] = ixCount;
              if (iso0cOK[indexFluxMax] > 0) {
                iso0c0 = FindMean(vector,ixCount);
                rmsiso0c0 = FindStd(vector,ixCount);
                if (rmsiso0c0 < 0.3) {
                  rmsiso0c0 = 0.3;
                }
                ix0Count = 0;
                for (index = 0; index < ixCount; index++) {
                  pSextractor = &sextractor_table[ix[index]];
                  if (fabs(pSextractor->iso0r - iso0c0) <= (2.0*rmsiso0c0)) {
                    vector[ix0Count] = pSextractor->iso0r;
                    ix0[ix0Count] = ix[index];
                    ix0Count++;
                  }
                }
                iso0cOK[indexFluxMax] = ix0Count;
                if (iso0cOK[indexFluxMax] > 0) {
                  rmsiso0c0 = FindStd(vector,ix0Count);
                  ix0Count = 0;
                  for (index = 0; index < ixCount; index++) {
                    pSextractor = &sextractor_table[ix[index]];
                    if (fabs(pSextractor->iso0r - iso0c0) <= (2.0*rmsiso0c0)) {
                      vector[ix0Count] = pSextractor->iso0r;
                      ix0[ix0Count] = ix[index];
                      ix0Count++;
                    }
                  }
                }
                iso0cOK[indexFluxMax] = ix0Count;
                if (iso0cOK[indexFluxMax] > 0) {
                  rmsiso0c0 = FindStd(vector,ix0Count);
                  ix0Count = 0;
                  for (index = 0; index < ixCount; index++) {
                    pSextractor = &sextractor_table[ix[index]];
                    if (fabs(pSextractor->iso0r - iso0c0) <= (3.0*rmsiso0c0)) {
                      vector[ix0Count] = pSextractor->iso0r;
                      ix0[ix0Count] = ix[index];
                      ix0Count++;
                    }
                  }
                }
                iso0cOK[indexFluxMax] = ix0Count;
                if (iso0cOK[indexFluxMax] > 0) {

                  iso0c[indexFluxMax] = FindMean(vector,ix0Count);
                  rmsiso0c[indexFluxMax] = FindStd(vector,ix0Count);
                  if (rmsiso0c[indexFluxMax] > 1.0) {
                    rmsiso0c[indexFluxMax] = 1.0;
                  }

#if 0
                  printf("Bin %d %d %d iso0c %f rmsiso0c %f ix0 %d\n",iii+1,jjj+1,indexFluxMax+1,iso0c[indexFluxMax],rmsiso0c[indexFluxMax],ix0Count);
#endif

                  for (index = 0; index < ixbinCount; index++) {
                    pSextractor = &sextractor_table[ixbin[index]];
                    if ((pSextractor->fluxmax >= fluxmaxgrid[indexFluxMax]) &&
                        (pSextractor->fluxmax <  fluxmaxgrid[indexFluxMax+1]) &&
                        (fabs(pSextractor->iso0 - iso0c[indexFluxMax]) < (3.0 * rmsiso0c[indexFluxMax]))) {
                      pSextractor->internal |= INTERNAL_FLAG_IXF5;
                      ixf5Count++;
                    }
                  }
                }
              }
#if 0
              if ((verbose) && (iii == 2) && (jjj == 2)) {
                printf("Bin %d %d, i %2d iso0c %f rmsiso0c %f, length(iso0) %d length ixf4 %d\n",iii+1,jjj+1,indexFluxMax,iso0c[indexFluxMax],rmsiso0c[indexFluxMax],ixbinCount,ixf5Count-ixf5old2);
              }
#endif
              ixf5old2 = ixf5Count;

            }
#if 0
            for (indexMag = 0; indexMag < (maggridsize-1); indexMag++) {
              printf("index %2d iso0grid %f plotOK1 %d fluxmaxc %f rmsfluxmaxc %f fluxmaxgrid %f plotOK2 %d iso0c %f rmsiso0c %f\n",indexMag,iso0grid[indexMag],iso0gridOK[indexMag],fluxmaxc[indexMag],rmsfluxmaxc[indexMag],fluxmaxgrid[indexMag],iso0cOK[indexMag],iso0c[indexMag],rmsiso0c[indexMag]);
            }
#endif
            if (doPlots) {
              cpgpanl(1,2);
              plot_nrecs1 = ixbinrCount;
#ifdef LOS_DEBUG_NUMBER
                xval1Debug = pDebugSextractor->iso0r;
                yval1Debug = pDebugSextractor->fluxmaxr;
#endif  /* LOS_DEBUG_NUMBER */
              for (plot_index1 = 0; plot_index1 < plot_nrecs1; plot_index1++) {
                pSextractor = &sextractor_table[ixbinr[plot_index1]];
                xval1[plot_index1] = pSextractor->iso0r;
                yval1[plot_index1] = pSextractor->fluxmaxr;
                if (plot_index1 == 0) {
                  xmin = xval1[0];
                  xmax = xval1[0];
                  ymin = yval1[0];
                  ymax = yval1[0];
                } else {
                  if (xmin > xval1[plot_index1]) {
                    xmin = xval1[plot_index1];
                  }
                  if (xmax < xval1[plot_index1]) {
                    xmax = xval1[plot_index1];
                  }
                  if (ymin > yval1[plot_index1]) {
                    ymin = yval1[plot_index1];
                  }
                  if (ymax < yval1[plot_index1]) {
                    ymax = yval1[plot_index1];
                  }
 
                }
              }
              plot_nrecs2 = ixbinf1Count;
              for (plot_index2 = 0; plot_index2 < plot_nrecs2; plot_index2++) {
                pSextractor = &sextractor_table[ixbinf1[plot_index2]];
#ifdef LOS_DEBUG_NUMBER
                if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
                  printf("At NUMBER %d in line %d\n",pSextractor->NUMBER,__LINE__);
                }
#endif /* LOS_DEBUG_NUMBER */
                xval2[plot_index2] = pSextractor->iso0f1;
                yval2[plot_index2] = pSextractor->fluxmaxf1;
                if (xmin > xval2[plot_index2]) {
                  xmin = xval2[plot_index2];
                }
                if (xmax < xval2[plot_index2]) {
                  xmax = xval2[plot_index2];
                }
                if (ymin > yval2[plot_index2]) {
                  ymin = yval2[plot_index2];
                }
                if (ymax < yval2[plot_index2]) {
                  ymax = yval2[plot_index2];
                }
              }
              for (indexMag = 0; indexMag < (maggridsize-1); indexMag++) {
                if (iso0gridOK[indexMag]) {
                  valx = iso0grid[indexMag];
                  valy = 3*rmsfluxmaxc[indexMag]+fluxmaxc[indexMag];
                  if (xmin > valx) {
                    xmin = valx;
                  }
                  if (xmax < valx) {
                    xmax = valx;
                  }
                  if (ymin > valy) {
                    ymin = valy;
                  }
                  if (ymax < valy) {
                    ymax = valy;
                  }
                  valy = -3*rmsfluxmaxc[indexMag]+fluxmaxc[indexMag];
                  if (ymin > valy) {
                    ymin = valy;
                  }
                  if (ymax < valy) {
                    ymax = valy;
                  }
                }
              }
              for (indexMag = 0; indexMag < (maggridsize-1); indexMag++) {
                if (iso0cOK[indexMag] > 0) {
                  valx = 3*rmsiso0c[indexMag]+iso0c[indexMag];
                  valy = fluxmaxgrid[indexMag];
                  if (xmin > valx) {
                    xmin = valx;
                  }
                  if (xmax < valx) {
                    xmax = valx;
                  }
                  if (ymin > valy) {
                    ymin = valy;
                  }
                  if (ymax < valy) {
                    ymax = valy;
                  }
                  valx = -3*rmsiso0c[indexMag]+iso0c[indexMag];
                  valy = fluxmaxgrid[indexMag+1];
                  if (xmin > valx) {
                    xmin = valx;
                  }
                  if (xmax < valx) {
                    xmax = valx;
                  }
                  if (ymin > valy) {
                    ymin = valy;
                  }
                  if (ymax < valy) {
                    ymax = valy;
                  }
                }
              }
              cpgenv(xmin-0.1,xmax+0.1,ymin-0.1,ymax+0.1,0,0);   
              sprintf(title,"%s [%d,%d]",fileroot,iii+1,jjj+1);
              cpglab("log(ISO0) (pixel**2)","FLUX_MAX",title);
              cpgsch(1.0); /* Character size is 1.0/40 the height of the view surface */
              cpgsci(4);  /* Set blue */
              cpgpt(plot_nrecs1,xval1,yval1,16);
              cpgsci(2); /* Set red */
#ifdef LOS_DEBUG_NUMBER
              if ((xval1Debug != 0) && (yval1Debug != 0)) {
                cpgpt(1,&xval1Debug,&yval1Debug,19);
                xval1Debug = 0.0;
                yval1Debug = 0.0;
              }
#endif /* LOS_DEBUG_NUMBER */
              cpgpt(plot_nrecs2,xval2,yval2,4);
              cpgsch(2.0); /* Character size is 2.0/40 the height of the view surface */

              cpgsci(1); /* Set black */
              for (indexMag = 0; indexMag < maggridsize-1; indexMag++) {
                xval1[0] = iso0grid[indexMag];
                xval1[1] = iso0grid[indexMag+1];
                yval1[0] = 3*rmsfluxmaxc[indexMag]+fluxmaxc[indexMag];
                yval1[1] = 3*rmsfluxmaxc[indexMag]+fluxmaxc[indexMag];
                cpgline(2,xval1,yval1);
                yval1[0] = -3*rmsfluxmaxc[indexMag]+fluxmaxc[indexMag];
                yval1[1] = -3*rmsfluxmaxc[indexMag]+fluxmaxc[indexMag];
                cpgline(2,xval1,yval1);
                xval1[0] = 3*rmsiso0c[indexMag]+iso0c[indexMag];
                xval1[1] = 3*rmsiso0c[indexMag]+iso0c[indexMag];
                yval1[0] = fluxmaxgrid[indexMag];
                yval1[1] = fluxmaxgrid[indexMag+1];
                cpgline(2,xval1,yval1);
                xval1[0] = -3*rmsiso0c[indexMag]+iso0c[indexMag];
                xval1[1] = -3*rmsiso0c[indexMag]+iso0c[indexMag];

                cpgline(2,xval1,yval1);
                

              }


            } /* End of plot */

         

          } /* niso0 test */

        } /* ixbinrCount test */
#if 1
        if (verbose) {
          printf("Bin %d %d, ixf1 %5d ixf2 %5d ixf3 %5d ixf4 %5d ixf5 %5d\n",iii+1,jjj+1,ixf1Count-ixf1old,ixf2Count-ixf2old,ixf3Count-ixf3old,ixf4Count-ixf4old,ixf5Count -ixf5old);
        }
#endif
        ixf1old = ixf1Count;
        ixf2old = ixf2Count;
        ixf3old = ixf3Count;
        ixf4old = ixf4Count;
        ixf5old = ixf5Count;
        cpgend();

      } /* Loop jjj on Y bins */
    } /* Loop iii on X bins */

    outHandle = fopen(outfile,"wt");
    if (outHandle == NULL) {
      printf("ERROR: Failed to open the output file %s\n",outfile);
      exit(-1);
    } else {
      if (verbose) {
        printf("Output file %s\n",outfile);
      }
    }
#ifdef FLAGS_OUTPUT
    fprintf(outHandle,"NUMBER\tAFLAGS\n");
    fprintf(outHandle,"------\t------\n");
#else /* FLAGS_OUTPUT */
    fprintf(outHandle,"NUMBER\tflag\n");
    fprintf(outHandle,"------\t----\n");
#endif /* FLAGS_OUTPUT */
    for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
      int flags = 0;
      pSextractor = &sextractor_table[sextractor_index];
#ifdef LOS_DEBUG_NUMBER
      if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
        printf("At NUMBER %d in line %d\n",pSextractor->NUMBER,__LINE__);
      }

#endif /* LOS_DEBUG_NUMBER */

      if ((pSextractor->internal & INTERNAL_FLAG_IXFALL) == INTERNAL_FLAG_IXFALL) {
        ixfAllCount++;
        flags = 1;
      } else {
        /* Mark this one as a defect.  After Sep 9, 2008 include high drad and blended stars  */
        pSextractor->AFLAGS |= (1<<FILTER_AFLAG_DEFECT);
        ixfFlaggedCount++;
      }
#ifdef FLAGS_OUTPUT
      fprintf(outHandle,"%d\t%d\n",pSextractor->NUMBER,pSextractor->AFLAGS); 
#else /* FLAGS_OUTPUT */
      fprintf(outHandle,"%d\t%d\n",pSextractor->NUMBER,flags); 
#endif /* FLAGS_OUTPUT */
    }

    /* Now calculate the summary parameters */
    summaryHandle = fopen(summaryfile,"wt");
    if (summaryHandle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to open the summary file %s\n",summaryfile);
    } else {
      if (verbose) {
        printf("Output file %s\n",summaryfile);
      }
    }
    fprintf(summaryHandle,"ratio1\tratio2\tratio3\tnumreal\tnumfake\tPlate\n");
    fprintf(summaryHandle,"------\t------\t------\t-------\t-------\t-----\n");
    index = 0;
    indexMag = 0;
    for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
      pSextractor = &sextractor_table[sextractor_index];
      if ((pSextractor->internal & (INTERNAL_FLAG_FLUX|INTERNAL_FLAG_REJECT)) != 0) {
        continue;
      }
      if (indexMag == 0) {
        maxmagiso = pSextractor->MAG_ISO;
      } else {
        if (pSextractor->MAG_ISO > maxmagiso) {
          maxmagiso = pSextractor->MAG_ISO;
        }
      }
      indexMag++;
      if (pSextractor->drad < 100.0) {
        vector[index] = pSextractor->drad;
        index++;
      }
    }
    CalcMedianAndRMS(index,0,vector,&meddrad3,&stddrad3,0,1.0,0);
    ixsrCount = 0;
    ixsfCount = 0;
    flagsrCount = 0;
    flagsfCount = 0;


    for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
      pSextractor = &sextractor_table[sextractor_index];
      if ((pSextractor->internal & INTERNAL_FLAG_FLUX) != 0) {
        continue;
      }
#ifdef NO_DRAD_FILTER
      if ((pSextractor->drad < 100.0) &&
          ((pSextractor->internal & INTERNAL_FLAG_REJECT) == 0) &&
          (pSextractor->MAG_ISO < (maxmagiso - 2.0)) &&
          (pSextractor->X_IMAGE > (minimumX + (0.05*(1.0*xrange)))) &&
          (pSextractor->X_IMAGE < (minimumX + (0.95*(1.0*xrange)))) &&
          (pSextractor->Y_IMAGE > (minimumY + (0.05*(1.0*yrange)))) &&
          (pSextractor->Y_IMAGE < (minimumY + (0.95*(1.0*yrange)))) &&
          (pSextractor->plate_dist < (0.85*max_plate_dist))) {
        ixsrCount++;
        if ((pSextractor->internal & INTERNAL_FLAG_IXFALL) == INTERNAL_FLAG_IXFALL) {
          flagsrCount++;
        }
        
      }
      if (((pSextractor->drad > 100.0) || ((pSextractor->internal & INTERNAL_FLAG_REJECT) != 0)) &&
          (pSextractor->MAG_ISO < (maxmagiso - 2.0)) &&
          (pSextractor->X_IMAGE > (minimumX + (0.05*(1.0*xrange)))) &&
          (pSextractor->X_IMAGE < (minimumX + (0.95*(1.0*xrange)))) &&
          (pSextractor->Y_IMAGE > (minimumY + (0.05*(1.0*yrange)))) &&
          (pSextractor->Y_IMAGE < (minimumY + (0.95*(1.0*yrange)))) &&
          (pSextractor->plate_dist < (0.85*max_plate_dist))) {
        ixsfCount++;
        if ((pSextractor->internal & INTERNAL_FLAG_IXFALL) != INTERNAL_FLAG_IXFALL) {
          flagsfCount++;
        }
      }


#else /* NO_DRAD_FILTER */

      if ((pSextractor->drad <= (2*meddrad3)) &&
          (pSextractor->MAG_ISO < (maxmagiso - 2.0)) &&
          (pSextractor->X_IMAGE > (minimumX + (0.05*(1.0*xrange)))) &&
          (pSextractor->X_IMAGE < (minimumX + (0.95*(1.0*xrange)))) &&
          (pSextractor->Y_IMAGE > (minimumY + (0.05*(1.0*yrange)))) &&
          (pSextractor->Y_IMAGE < (minimumY + (0.95*(1.0*yrange)))) &&
          (pSextractor->plate_dist < (0.85*max_plate_dist))) {
        ixsrCount++;
        if ((pSextractor->internal & INTERNAL_FLAG_IXFALL) == INTERNAL_FLAG_IXFALL) {
          flagsrCount++;
        }
        
      }
      if ((pSextractor->drad > (10*meddrad3)) &&
          (pSextractor->MAG_ISO < (maxmagiso - 2.0)) &&
          (pSextractor->X_IMAGE > (minimumX + (0.05*(1.0*xrange)))) &&
          (pSextractor->X_IMAGE < (minimumX + (0.95*(1.0*xrange)))) &&
          (pSextractor->Y_IMAGE > (minimumY + (0.05*(1.0*yrange)))) &&
          (pSextractor->Y_IMAGE < (minimumY + (0.95*(1.0*yrange)))) &&
          (pSextractor->plate_dist < (0.85*max_plate_dist))) {
        ixsfCount++;
        if ((pSextractor->internal & INTERNAL_FLAG_IXFALL) != INTERNAL_FLAG_IXFALL) {
          flagsfCount++;
        }
      }
#endif /* NO_DRAD_FILTER */

    }
    fprintf(summaryHandle,"%f\t%f\t%f\t%d\t%d\t%s\n",
            (1.0 * flagsrCount)/(1.0 * ixsrCount),
            (1.0 * flagsfCount)/(1.0 * ixsfCount),
            (1.0 * ixsfCount)/(1.0 *ixsrCount),
            ixsrCount,
            ixsfCount,
            fileroot);
          
    

  } else {
    printf("ERROR: Insufficient low drad stars  %d for plate %s\n",kxr,fileroot);

  }


  /* All done.  Clean up */
 

  if (sextractor_table != NULL) {
    Free(sextractor_table);
  }
  if (sextractor_header != NULL) {
    table_hdrfree(sextractor_header);
  }
  if (sextractor_handle != NULL) {
    Close(sextractor_handle);
  }

  if (match_table != NULL) {
    Free(match_table);
  }
  if (match_header != NULL) {
    table_hdrfree(match_header);
  }
  if (match_handle != NULL) {
    Close(match_handle);
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }

  if (summaryHandle != NULL) {
    fclose(summaryHandle);
  }
  if (vector != NULL) {
    free(vector);
  }
  if (ixbin != NULL) {
    free(ixbin);
  }
  if (ixbin2 != NULL) {
    free(ixbin2);
  }
  if (ixbinf != NULL) {
    free(ixbinr);
  }
  if (ix != NULL) {
    free(ix);
  }
  if (ix0 != NULL) {
    free(ix0);
  }


  if (matchVector != NULL) {
    free(matchVector);
  }
  if (xval1 != NULL) {
    free(xval1);
  }
  if (xval2 != NULL) {
    free(xval2);
  }
  if (yval1 != NULL) {
    free(yval1);
  }
  if (yval2 != NULL) {
    free(yval2);
  }
  if (ixa != NULL) {
    free(ixa);
  }
  if (ixb != NULL) {
    free(ixb);
  }
  if (ixbinf1 != NULL) {
    free(ixbinf1);
  }


  time(&curTime);
  curTime -= startTime;
  fprintf(stdout,"Stars %d, matched %d, dradFilter %d lowFlux %d filters %d %d %d %d %d all %d flagged %d seconds %d for %s\n",
          sextractor_nrecs,
          match_nrecs,
          dradFilterCount,
          lowFluxCount,
          ixf1Count,
          ixf2Count,
          ixf3Count,
          ixf4Count,
          ixf5Count,
          ixfAllCount,
          ixfFlaggedCount,
          curTime,
          fileroot);
    


  return(0);
}

