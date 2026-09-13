// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* filter_ccmap.c
 *
 *  Select matched points for iraf ccmap coordinate refinement
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm -lmysqlclient filter_ccmap.c pipelineutils.a -ltable -lutil -lwcs -o filter_ccmap 
 *
 * /dasch/Pipeline/filter_ccmap -v -i /dasch/Pipeline/match/match_ac42227_00_01ww_u.db -o /dasch/data/ExposureData/Mosaics/ac/42227_00/ac42227_00_01ww.coo -l 2 -u 8 -t -w 17493 -h 21953 -g /dasch/Pipeline/match/ac42227_00_01ww_ccmapgmt.txt
 *
 * /dasch/Pipeline/filter_ccmap -v -i /dasch/Pipeline/match/match_ac42227_00_01ww_u.db -o /dasch/data/ExposureData/Mosaics/ac/42227_00/ac42227_00_01ww.coo -l 2 -u 9 -t -w 17493 -h 21953 -g /dasch/Pipeline/match/ac42227_00_01ww_ccmapgmt.txt
 *
 *
 * Aug 26, 2008 Edward J. Los - Initial version
 * Oct 11, 2008 Edward J. Los - Eliminate thresholding for now
 * Feb  6, 2009 Edward J. Los - Replace FILTER_BFLAG_BLEND with FILTER_AFLAG_BLEND
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

/* #define LOS_DEBUG 1  */

#define MAX_BUFFER 100
#define STARS_PER_FLUX_BIN 20
#define NUMBER_FLUX_BINS 100
#define FLUX_BIN_ALLOC 1000
#define MIN_VECTOR 10
#define SKIP_SATURATION 1 /* The threshold idea did not pan out.  These headers are not reflected in ./mysql/variables.txt */

/* Define only one of the following */
/* #define GMT_PLOT_A 1 */ /* Original Input Y_IMAGE vs X_IMAGE */
/* #define GMT_PLOT_B 1 */ /* Initial Filter FWHM_IMAGE vs FLUX_MAX */
/* #define GMT_PLOT_C 1 */ /* Initial filter Y_IMAGE vs X_IMAGE */
/* #define GMT_PLOT_D 1 */  /* Sorted FLUX_MAX vs FLUX_IMAGE */
/* #define GMT_PLOT_E 1 */      /* Output FLUX_MAX vs FLUX_IMAGE */
#define GMT_PLOT_F 1  /* Output X_IMAGE and Y_IMAGE */

typedef struct _starimage {
  double X_IMAGE;        /* Sextractor X location in pixels */
  double Y_IMAGE;        /* Sextractor Y location in pixels */
  double ra_2;           /* GSC Right Ascension in degrees */
  double dec_2;          /* GSC Declination in degrees */
  double FLUX_MAX; 
  double THRESHOLD;
  double FWHM_IMAGE;
  int AFLAGS;
  int BFLAGS;
  int class;
  int spatial_bin;
  int selected;
} STARIMAGE,*PSTARIMAGE;

int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char infile[MAX_BUFFER];
  char cmdchar;
  int errorFlag = 0;
  char *matchDirectory;
  int numBins = 0;
  char outfile[MAX_BUFFER];
  char sumfile[MAX_BUFFER];
  char gmtfile[MAX_BUFFER];
  char rootname[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  char extinctionfile[MAX_BUFFER];
  FILE *outHandle = NULL;
  FILE *sumHandle = NULL;
  FILE *gmtHandle = NULL;

  File match_handle = NULL;
  char match_name[MAX_BUFFER];
  TableHead match_header = NULL;
  PSTARIMAGE match_table = NULL;
  size_t match_nrecs = 0;
  int match_index;
  PSTARIMAGE pMatch;

  int spatial_bin;
  int verbose = 0;

  int mosaicWidth = 0;
  int mosaicHeight = 0;

  time_t startTime;
  time_t curTime;
  double edgeDist;
  int highBinNumber = 0;
  char *dotPtr;

  

  int enableThreshold = 0;
  int minimumBin = 1;
  int maximumBin = MAX_SPATIAL_BINS;
  double maximumFlux = 0;
  double baseFlux = 65536.;
  int outputCount = 0;
  int selectCount = 0;
  int rejectCount = 0;
  int binNumber;
  double saturationFlux = 0;
  int fluxBinCount[NUMBER_FLUX_BINS];
  int fluxBinClipCount[NUMBER_FLUX_BINS];
  int fluxBinAlloc[NUMBER_FLUX_BINS];
  double * fluxBinVector[NUMBER_FLUX_BINS];
  double * tempBinVector = NULL;
  double fluxBinMedian[NUMBER_FLUX_BINS];
  double fluxBinRMS[NUMBER_FLUX_BINS];
  double vector[NUMBER_FLUX_BINS];
  int vectorCount = 0;
  int clipVectorCount = 0;
  int index;
  double vectorMedian;
  double vectorRMS;
  double maxMedian = 0;



  time(&startTime);

  for (index = 0; index < NUMBER_FLUX_BINS; index++) {
    fluxBinCount[index] = 0;
    fluxBinAlloc[index] = FLUX_BIN_ALLOC;
    fluxBinClipCount[index] = 0;
    fluxBinRMS[index] = 0.0;
    fluxBinMedian[index] = 0.0;
    fluxBinVector[index] = calloc(FLUX_BIN_ALLOC,sizeof(double));
    if (fluxBinVector[index] == NULL) {
      fprintf(stderr,"ERROR: Failed to allocate fluxBinVector for entry %d\n",index);
      exit(-1);
    }
  }
    

  match_name[0] = 0;
  outfile[0] = 0;
  gmtfile[0] = 0;
  sumfile[0] = 0;


  matchDirectory = getenv("DASCH_MATCH");
  if (matchDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_MATCH is not defined\n");
    return(-1);
  }

 
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

        case 's': /* summary file name */
        case 'S':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sumfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'g': /* plot file name */
        case 'G':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(gmtfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'i': /* match file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(match_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicWidth);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode mosaic width %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* mosaic height */
        case 'H':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicHeight);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'l': /* minimum bin */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&minimumBin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode minimum bin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'u': /* maximum bin */
        case 'U':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&maximumBin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode maximum bin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 't': /* Perform threshold filtering */
        case 'T':
          enableThreshold = 1;
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;


        default:
          fprintf(stderr,"ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */
  if (mosaicWidth == 0) {
    fprintf(stderr,"ERROR: No mosaic width specified\n");
    errorFlag = 1;
  }
  if (mosaicHeight == 0) {
    fprintf(stderr,"ERROR: No mosaic height specified\n");
    errorFlag = 1;
  }

  if (outfile[0] == 0) {
    fprintf(stderr,"ERROR: No output filename was specified\n");
    errorFlag = 1;
  }
#ifndef SKIP_SATURATION
  if (sumfile[0] == 0) {
    fprintf(stderr,"ERROR: No summary filename was specified\n");
    errorFlag = 1;
  }
#endif /* SKIP_SATURATION */
  if (match_name[0] == 0) {
    fprintf(stderr,"ERROR: No match filename was specified\n");
    errorFlag = 1;
  }

  if ((minimumBin < 1) ||
      (maximumBin < 1) ||
      (minimumBin > MAX_SPATIAL_BINS) ||
      (maximumBin > MAX_SPATIAL_BINS) ||
      (minimumBin > maximumBin)) {
    fprintf(stderr,"ERROR: Illegal minimum bin %d or maximum bin %d\n",minimumBin,maximumBin);
    errorFlag = 1;
  }
  strcpy(rootname,outfile);
  dotPtr = strstr(rootname,".");
  if (dotPtr != NULL) {
    *dotPtr = 0;
  }



  /* Open the match file */
  match_handle = Open(match_name,"r");
  if (match_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the match file %s\n",match_name);
  } else {
    if (verbose) {
      fprintf(stderr,"Found match file %s\n",match_name);
    }
  }

  outHandle = fopen(outfile,"wt");
  if (outHandle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
  } else {
    if (verbose) {
      fprintf(stderr,"Output file %s\n",outfile);
    }
  }

#ifndef SKIP_SATURATION
  sumHandle = fopen(sumfile,"wt");
  if (sumHandle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to open the sumput file %s\n",sumfile);
  } else {
    if (verbose) {
      fprintf(stderr,"Summary file %s\n",sumfile);
    }
  }
#endif /* SKIP_SATURATION */   

  if (gmtfile[0] != 0) {
    gmtHandle = fopen(gmtfile,"wt");
    if (gmtHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the plot file %s\n",gmtfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output file %s\n",gmtfile);
      }
    }
  
  }


  if (errorFlag) {
    fprintf(stderr,"Usage: filter_ccmap -i <match file> \n");
    fprintf(stderr,"                    -o <output file> \n");
    fprintf(stderr,"                    -s <summary file> \n");
    fprintf(stderr,"                    -w <mosaic width in pixels> \n");
    fprintf(stderr,"                    -h <mosaic height in pixels> \n");
    fprintf(stderr,"                    -l <smallest bin for FLUX_MAX filtering> \n");
    fprintf(stderr,"                    -u <largest bin for FLUX_MAX filtering> \n");
    fprintf(stderr,"                    -t : enable threshold filtering \n");
    fprintf(stderr,"                   [-g <gmt plotting file>]\n");
    fprintf(stderr,"                    -v : verbose\n");

    return(-1);
  }

  printf("filter_ccmap of %s %s lower bin %d upper bin %d threshold %d\n",
         __DATE__,__TIME__,minimumBin,maximumBin,enableThreshold);
#ifdef LOS_DEBUG
  printf("ERROR: LOS_DEBUG is set\n");
#endif /* LOS_DEBUG */

#ifdef SKIP_SATURATION_X
  printf("ERROR: SKIP_SATURATION is set\n");
#endif /* SKIP_SATURATION */
#ifndef SKIP_SATURATION
  fprintf(sumHandle,"baseFlux\tmaxFlux\tsaturationFlux\thighBin\tnrecs\toutputCount\tPlate\n");
  fprintf(sumHandle,"--------\t-------\t--------------\t-------\t-----\t-----------\t-----\n");
#endif /* SKIP_SATURATION */
  /* Now read in the match file */

  match_header = table_header(match_handle,TABLE_PARSE);
  if (match_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",match_name);
    return(-1);
  }

  match_table = table_loadva(match_handle,
                             &match_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(STARIMAGE),
                             &match_nrecs,
                             TblDbl,"X_IMAGE",TblOff(PSTARIMAGE,X_IMAGE),
                             TblDbl,"Y_IMAGE",TblOff(PSTARIMAGE,Y_IMAGE),
                             TblInt,"AFLAGS"  ,TblOff(PSTARIMAGE,AFLAGS),
                             TblInt,"BFLAGS"  ,TblOff(PSTARIMAGE,BFLAGS),
#ifndef LOS_DEBUG
                             TblInt,"class"  ,TblOff(PSTARIMAGE,class),
#endif /* LOS_DEBUG */
                             TblDbl,"ra_2"   ,TblOff(PSTARIMAGE,ra_2),
                             TblDbl,"dec_2"  ,TblOff(PSTARIMAGE,dec_2),
#ifndef LOS_DEBUG
                             TblDbl,"FLUX_MAX"  ,TblOff(PSTARIMAGE,FLUX_MAX),
                             TblDbl,"THRESHOLD"  ,TblOff(PSTARIMAGE,THRESHOLD),
#endif /* LOS_DEBUG */
                             TblDbl,"FWHM_IMAGE"  ,TblOff(PSTARIMAGE,FWHM_IMAGE),
                             0,"end",0);
  if (match_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",match_name);
    return(-1);
  }
  if (verbose) {
    fprintf(stderr,"read %d records for %s\n",match_nrecs,match_name);
  }
#ifdef GMT_PLOT_A
  if (gmtHandle != NULL) {
    for (match_index = 0; match_index < match_nrecs; match_index++) {
      pMatch = &match_table[match_index];
      fprintf(gmtHandle,"%.6f %.6f\n",
              pMatch->X_IMAGE,pMatch->Y_IMAGE);

    }
  }
#endif /* GMT_PLOT_A */

  
  for (match_index = 0; match_index < match_nrecs; match_index++) {
    pMatch = &match_table[match_index];
    pMatch->selected = 0;
#ifdef LOS_DEBUG
    pMatch->class = 0;
    pMatch->FLUX_MAX = 2.0 * CCMAP_THRESHOLD;
    pMatch->THRESHOLD = 1.0;
#endif /* LOS_DEBUG */
    pMatch->spatial_bin = CalculateBin(mosaicWidth,mosaicHeight,pMatch->X_IMAGE,pMatch->Y_IMAGE,&edgeDist);

    if (enableThreshold) {
      /* Do not pick blended objects, non-stars, and low flux objects */
      if (((pMatch->AFLAGS & (1 << FILTER_AFLAG_BLEND)) == 0) &&
          (pMatch->FLUX_MAX > CCMAP_THRESHOLD * pMatch->THRESHOLD) &&
          (pMatch->class == 0)) {
        pMatch->selected = 1;
        selectCount++;
      }
    } else {
      /* Select everything */
      pMatch->selected = 1;
      selectCount++;
    }
    if (pMatch->spatial_bin > maximumBin) {
      if (pMatch->selected != 0) {
        selectCount--;
      }
      pMatch->selected = 0;
    }
    
  }
#if (defined(GMT_PLOT_C) || defined(GMT_PLOT_B))
  if (gmtHandle != NULL) {
    for (match_index = 0; match_index < match_nrecs; match_index++) {
      pMatch = &match_table[match_index];
      if (pMatch->selected != 0) {
#ifdef GMT_PLOT_B
        fprintf(gmtHandle,"%.6f %.6f\n",
                pMatch->FLUX_MAX,
                pMatch->FWHM_IMAGE);
#else /* GMT_PLOT_B */
        fprintf(gmtHandle,"%.6f %.6f\n",
                pMatch->X_IMAGE,pMatch->Y_IMAGE);
#endif /* GMT_PLOT_B */
      }

    }
  }
#endif /* GMT_PLOT_C & GMT_PLOT_B */
  if (selectCount == 0) {
    fprintf(stderr,"ERROR: initial selection count is 0 for %s\n",outfile);
    exit(-1);
  }


                                         
  /* Find the maximum and minimum flux */
  for (match_index = 0; match_index < match_nrecs; match_index++) {
    pMatch = &match_table[match_index];
    if (pMatch->selected != 0) {
      if (pMatch->FLUX_MAX > maximumFlux) {
        maximumFlux = pMatch->FLUX_MAX;
      }
      if (pMatch->FLUX_MAX < baseFlux) {
        baseFlux = pMatch->FLUX_MAX;
      }
    }
  }
  if (maximumFlux == baseFlux) {
    fprintf(stderr,"ERROR: minimum %f and maximum %f fluxes are equal\n",baseFlux,maximumFlux);
    exit(-1);
  }
  maximumFlux += 0.001; /* Make sure we never have a round-up problem */
  /* Now populate the bins */
  for (match_index = 0; match_index < match_nrecs; match_index++) {
    pMatch = &match_table[match_index];
    if (pMatch->selected != 0) {
      binNumber = ((pMatch->FLUX_MAX-baseFlux) * NUMBER_FLUX_BINS) /(maximumFlux-baseFlux);
      if ((binNumber < 0) || (binNumber >= NUMBER_FLUX_BINS)) {
        fprintf(stderr,"ERROR in calculating binNumber %f %f %f\n",baseFlux,maximumFlux,pMatch->FLUX_MAX);
      }
      if ((fluxBinCount[binNumber]+1) >= fluxBinAlloc[binNumber]) {
        fluxBinAlloc[binNumber] += FLUX_BIN_ALLOC;
        /* Need to get more bin space */
        tempBinVector = realloc(fluxBinVector[binNumber],fluxBinAlloc[binNumber]*sizeof(double));
        if (tempBinVector == NULL) {
          fprintf(stderr,"ERROR: failed to reallocate tempBinVector of size %d for %s\n",fluxBinAlloc[binNumber],rootname);
        }
        fluxBinVector[binNumber] = tempBinVector;
        tempBinVector = NULL;
      }
      fluxBinVector[binNumber][fluxBinCount[binNumber]] = pMatch->FWHM_IMAGE;
      fluxBinCount[binNumber]++;
    }
  }
  /* Next calculate the median and RMS for each bin */
  for (binNumber = 0; binNumber < NUMBER_FLUX_BINS; binNumber++) {
    
    fluxBinClipCount[binNumber] = CalcMedianAndRMS(fluxBinCount[binNumber],STARS_PER_FLUX_BIN,fluxBinVector[binNumber],&fluxBinMedian[binNumber],&fluxBinRMS[binNumber],0,3.0,0);
    if (maxMedian < fluxBinMedian[binNumber]) {
      maxMedian = fluxBinMedian[binNumber];
    }

#ifdef GMT_PLOT_D
    if ((gmtHandle != NULL) && (fluxBinClipCount[binNumber] > STARS_PER_FLUX_BIN)) {
      fprintf(gmtHandle,"%.6f %.6f\n",((gmtHandle,binNumber * (maximumFlux - baseFlux))/NUMBER_FLUX_BINS)+baseFlux,fluxBinMedian[binNumber]);
      fprintf(gmtHandle,"%.6f %.6f\n",((gmtHandle,binNumber * (maximumFlux - baseFlux))/NUMBER_FLUX_BINS)+baseFlux,50.0+fluxBinRMS[binNumber]);
    }

#endif /* GMT_PLOT_D */
    if (fluxBinClipCount[binNumber] >= STARS_PER_FLUX_BIN) {
      vector[vectorCount] = fluxBinMedian[binNumber];
      vectorCount++;
    }
    
  }
  if (vectorCount < MIN_VECTOR) {
#ifdef SKIP_SATURATION
    fprintf(stderr,"WARNING: filter_ccmap vectorCount is %d for %s\n",vectorCount,rootname);
#else /* SKIP_SATURATION */
    fprintf(stderr,"ERROR: filter_ccmap vectorCount is %d for %s\n",vectorCount,rootname);
    exit(-1);
#endif /* SKIP_SATURATION */
  } else {
    clipVectorCount = CalcMedianAndRMS(vectorCount,0,vector,&vectorMedian,&vectorRMS,1,1.0,0);
    if (clipVectorCount <= 0) {
      fprintf(stderr,"ERROR: clipVectorCount is %d for %s\n",clipVectorCount,rootname);
      exit(-1);
    }
    /* Now find the highest bin for which we are below the one-sigma point */
    for (binNumber = 0; binNumber < NUMBER_FLUX_BINS; binNumber++) {
#if defined(GMT_PLOT_D) || defined(GMT_PLOT_E)
      if (gmtHandle != NULL) {
        fprintf(gmtHandle,"%.6f %.6f\n",((gmtHandle,binNumber * (maximumFlux - baseFlux))/NUMBER_FLUX_BINS)+baseFlux,vectorMedian);
        fprintf(gmtHandle,"%.6f %.6f\n",((gmtHandle,binNumber * (maximumFlux - baseFlux))/NUMBER_FLUX_BINS)+baseFlux,(vectorMedian+(3.0*vectorRMS)));
        fprintf(gmtHandle,"%.6f %.6f\n",((gmtHandle,binNumber * (maximumFlux - baseFlux))/NUMBER_FLUX_BINS)+baseFlux,(vectorMedian-(3.0*vectorRMS)));
      }
#endif /* GMT_PLOT_D or GMT_PLOT_E */
      if ((fluxBinMedian[binNumber] <= (vectorMedian+vectorRMS)) &&
          (fluxBinClipCount[binNumber] > STARS_PER_FLUX_BIN)) {
        highBinNumber = binNumber;
      }
    }
  }
#ifdef SKIP_SATURATION
  saturationFlux = 99999.0;
#else /* SKIP_SATURATION */
  saturationFlux = ((highBinNumber * (maximumFlux - baseFlux))/NUMBER_FLUX_BINS)+baseFlux;
#endif /* SKIP_SATURATION */
#if defined(GMT_PLOT_D) || defined(GMT_PLOT_E)
  if (gmtHandle != NULL) {
    for (index = 0; index < 50; index++) {
      fprintf(gmtHandle,"%.6f %.6f\n",saturationFlux,(maxMedian*index)/(50.0));
    }
  }
#endif /* GMT_PLOT_D or GMT_PLOT_E  */
  
  if (saturationFlux > baseFlux) {
    for (match_index = 0; match_index < match_nrecs; match_index++) {
      pMatch = &match_table[match_index];
#if 0
      /* The following is a bad idea: does not get corners */
      if ((pMatch->selected != 0) &&
          ((pMatch->FLUX_MAX > saturationFlux) ||
           (pMatch->FWHM_IMAGE > (vectorMedian+(3.0*vectorRMS))) ||
           (pMatch->FWHM_IMAGE < (vectorMedian-(3.0*vectorRMS))))) {
        pMatch->selected = 0;
        rejectCount++;
      } 

#endif

      if ((pMatch->selected != 0) &&
          (pMatch->FLUX_MAX > saturationFlux)) {
        pMatch->selected = 0;
        rejectCount++;
      } 

      if (pMatch->selected != 0) {
#ifdef GMT_PLOT_E
        if (gmtHandle != NULL) {
          fprintf(gmtHandle,"%.6f %.6f\n",pMatch->FLUX_MAX,pMatch->FWHM_IMAGE);
        }
#endif /* GMT_PLOT_E */
#ifdef GMT_PLOT_F
        if (gmtHandle != NULL) {
          fprintf(gmtHandle,"CIRCLE(%.6f,%.6f,5)\n",pMatch->X_IMAGE,pMatch->Y_IMAGE);
        }
#endif /* GMT_PLOT_F */
        fprintf(outHandle,"%.6f\t%.6f\t%.6f\t%.6f\n",
                pMatch->X_IMAGE,
                pMatch->Y_IMAGE,
                pMatch->ra_2,
                pMatch->dec_2);
        outputCount++;


      }


    }


  } else {
    fprintf(stderr,"WARNING: Insufficient flux %f for %s\n",
            saturationFlux,
            rootname);
  }

#ifndef SKIP_SATURATION
  fprintf(sumHandle,"%.1f\t%.1f\t%.1f\t%d\t%d\t%d\t%s\n",
          baseFlux,maximumFlux,saturationFlux,highBinNumber,match_nrecs,outputCount,rootname);
#endif /* SKIP_SATURATION */


  /* We are done with the match table.  Clear memory */
  if (match_table != NULL) {
    Free(match_table);
    match_table = NULL;
  }
  if (match_header != NULL) {
    table_hdrfree(match_header);
    match_header = NULL;
  }
  if (match_handle != NULL) {
    Close(match_handle);
    match_handle = NULL;
  }
  pMatch = NULL;

  if (outHandle != NULL) {
    fclose(outHandle);
  }

  if (gmtHandle != NULL) {
    fclose(gmtHandle);
  }
#ifndef SKIP_SATURATION
  if (sumHandle != NULL) {
    fclose(sumHandle);
  }
#endif /* SKIP_SATURATION */


  for (index = 0; index < NUMBER_FLUX_BINS; index++) {
    if (fluxBinVector[index] != NULL) {
      free(fluxBinVector[index]);
    }
  }

  time(&curTime);
  curTime -= startTime;
  fprintf(stderr,"filter_ccmap baseFlux %f maximumFlux %f saturationFlux %f highBin %d in %d selected %d rejected %d out %d seconds %d for %s\n",
          baseFlux,
          maximumFlux,
          saturationFlux,
          highBinNumber,
          match_nrecs,
          selectCount,
          rejectCount,
          outputCount,
          curTime,
          outfile);
    


  return(0);
}

