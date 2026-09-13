// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* estimate_lightcurves.c
 *
 * Perform photometry by comparing an unknown star with MAG_ISO values of nearby comparison stars.
 *
 *  gcc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include -I /dasch/install/pgplot -I/usr/include/mysql estimate_lightcurves.c pipelineutils.a -L /dasch/install/lib -lm  -L/usr/lib${lib64}/mysql -L /dasch/install/pgplot -lmysqlclient -lcpgplot -lpgplot  -L/usr/X11R6/lib  -lX11 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 /usr/lib${lib64}/libg2c.so.0 -ldl -pthread  -ltable -lutil -lwcs -o estimate_lightcurves 
 *
 *   Procedure: 
  
    1.  Select a variable of interest

    2.  Use the web site to create a tarball of fits files.  The extraction radius should be large enough to include sufficient comparison stars.

    3.  Store the tarball in a working directory

    4.  Run the program

    For N233012397 9:53:10 +33:53:53 normally 12.2 mag, similar to epsilon Aurigae

    /home/scanner/Pipeline/estimate_lightcurves -v -s  -d /home/scanner/junk/epsilondasch/ 
    evince lc_N233012397.ps &
    evince sequence_plot.ps &
    evince limiting_plot.ps &
 
   /home/scanner/Pipeline/estimate_lightcurves -v -s -r 20 -q kepler -d /home/scanner/junk/CHCyg/
    extraction width .498 degrees 900 arcsec radius
    N0303103149 for GSC
    K11913210 for Kepler  7-10.6 dimmer in recent
    K11913802 - mag 12.5 in sequence, but not in plot.
 
    DASCH_J192433.1+501430 for APASS (No match!)
 * 
 * Apr  6, 2009  Edward J. Los - Initial version
 * Apr 25, 2012  Edward J. Los - extensive revisions for V1033 Sco
 * May 14, 2012  Edward J. Los - Integrate with estimate_lightcurves.csh and automate the process as much as possible
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 */


#include <math.h>
#include <time.h>
#include <dirent.h>
#include <sys/stat.h>
#include "table.h"
#include "mysql.h"
#include "cpgplot.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "daschunistd.h"
#define MAX_PLATE_NAME 32
#define MAX_BUFFER 500
#define IGNORE_LIMITING  1   /* Debug only! */ 
#define MAX_LIST_INDEX 3000
/* #define TOPCAT_OUTPUT 1 */
#define DEFAULT_REFCNT 10
#define DEFAULT_FRACTION 0.9
#define DEFAULT_THRESHOLD 10
#define DIMMEST_MAGNITUDE 20
#define DEFAULT_ANALYSIS_THRESHOLD 2.2
#define BRIGHTEST_MAGNITUDE 3

typedef struct _listentry {
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  char src_name[MAX_SRC_LENGTH];   /* GSC2.3.2 reference number */
  double ra; 
  double dec;
  double Stdmag;
  int list_index;
  int queued;
  int brightnessIndex;
  int refcnt;
  int gsc_bin_index;
  double rawmed;
  double rawrms;
  struct _listentry *flink;
} LISTENTRY,*PLISTENTRY;

typedef struct _magisoentry {
	double MAG_ISO_med;
	double MAG_ISO_rms;
  char Plate[MAX_PLATE_NAME];
} MAGISOENTRY,*PMAGISOENTRY;

typedef struct _master {
  double Date;
  double FLUX_ISO;
  int spatial_bin;
  int gsc_bin_index;
  int AFLAGS;
  int BFLAGS;
  int list_index;
  double ra;
  double dec;
  double Stdmag;
  double FLUX_MAX;
  double THRESHOLD;
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  char Plate[MAX_PLATE_NAME];
} MASTER,*PMASTER;

typedef struct _magentry {
  int master_index;
  int list_index;
  int spatial_bin;
  int AFLAGS;
  double FLUX_ISO;
  double FLUX_MAX;
  double THRESHOLD;
  double Date;
  
} MAGENTRY,*PMAGENTRY;

typedef struct _plate {
  double Date;
  char fitsFileName[MAX_PLATE_NAME];
  char Plate[MAX_PLATE_NAME];
  MAGENTRY magentry[MAX_LIST_INDEX];
  MAGENTRY magentrySort[MAX_LIST_INDEX];
  double magEstimate[MAX_LIST_INDEX];
  double flag[MAX_LIST_INDEX];
  /* Dummy entries */
  double magcal_local_rms;
  int npoints_local;
  double magcal_iso_rms;
  double magcal_local_error;
  double magcal_iso;
  double limiting_mag_local1; /* estimate from MAG_ISO_med+MAG_ISO_rms */
  double limiting_mag_local2; /* estimate from dimmest matched star in sequence */
  double reject_reason1;
  double reject_reason2;
  double dradRMS2;
	double FLUX_ISO_med;
} PLATE,*PPLATE;

/* This will crash the stack if inside a routine */
int lowestCount[MAX_LIST_INDEX][MAX_LIST_INDEX];
int brightnessOrder[MAX_LIST_INDEX];


double ComputeMag(double testFlux,
                  double lowFlux,
                  double highFlux,
                  double lowStdmag,
                  double highStdmag)
{
  double logLowFlux = log(lowFlux);
  double logHighFlux = log(highFlux);
  double logTestFlux = log(testFlux);
  double testStdmag;
  if (logLowFlux != logHighFlux) {
    testStdmag = lowStdmag+ ((highStdmag - lowStdmag) * (logTestFlux - logLowFlux)/(logHighFlux - logLowFlux));
  } else {
    testStdmag = lowStdmag;
  }

  return(testStdmag);


}


int ReadMaster(PMASTER pMaster,
               File master_handle,
               TableHead master_header,
               TblDescriptor master_descriptor,
               TableRow* master_row)
{
  *master_row = table_rowget(master_handle,master_header,*master_row,NULL,NULL,0);
  if (*master_row == NULL) {
    return(0);
  }
  if (!table_loadrow(master_handle,master_header,*master_row,master_descriptor,(char *)pMaster)) {
    printf("ERROR: Read Master table_loadrow failed\n");
    return(0);
  }
  return(1);
}
/* Sort routine based on the Julian Day */
int PlateCompare(const void *first, const void *second) 
{
  double numberFirst = ((PPLATE)first)->Date;
  double numberSecond = ((PPLATE)second)->Date;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}

/* Sort routine based on Stdmag */
int ListCompare(const void *first, const void *second) 
{
  double numberFirst = ((PLISTENTRY)first)->Stdmag;
  double numberSecond = ((PLISTENTRY)second)->Stdmag;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}


/* Sort routine based on MAG_ISO (reverse sense) */
int MagEntryCompare(const void *first, const void *second) 
{
  double numberFirst = ((PMAGENTRY)first)->FLUX_ISO;
  double numberSecond = ((PMAGENTRY)second)->FLUX_ISO;
  if (numberFirst < numberSecond) {
    return(1);
  } else if (numberFirst > numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}

int fitsfilter(const struct dirent *direntry) 
{
  char* charPtr;
  if ((direntry->d_type == DT_REG) ||
      (direntry->d_type == DT_UNKNOWN)) {
    if (((charPtr = strstr((char *)direntry->d_name,".fit")) != NULL) &&
        (charPtr[4] == 0)) {
      return(1);
    }
  }
  return(0);

}

int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  char workingDirectory[MAX_BUFFER];
  int length;
  char plateListFile[MAX_BUFFER];
  char plotdev[MAX_BUFFER];
  FILE *plateListHandle = NULL;
  char rmsName[MAX_BUFFER];
  char label[MAX_BUFFER];
  FILE *rmsHandle = NULL;
  int min_refcnt = DEFAULT_REFCNT;
  double fraction = DEFAULT_FRACTION;
  double fraction1;
  double fraction2;
  double analysis_threshold = DEFAULT_ANALYSIS_THRESHOLD;
  int base_min_threshold = -1;
  int max_min_threshold;
  int cur_min_threshold;
  
  int min_threshold;
  int best_selection_count = 0;
  int best_min_threshold = 0;

  char *dotPos;
  char *slashPos;
  char *lastSlashPos;
  double *vector = NULL;

  time_t startTime;
  time_t curTime;
  int verbose = 0;
  int printReject = 0;
  int debugFlag = 0;
  int outputPlotCount = 0;
  int outputResult;
  char *ingestDirectory;
  char *catalogallDirectory;
  
  char inLine[MAX_BUFFER];
  int lineLen;
  int nlines = 0;
  char *inBuffer;
  double *vector1 = NULL;
  double *vector2 = NULL;
  int plateListLength = 0;
  int compareFlag;
  int match_count = 0;
  int plateIndex;
  char output_name[MAX_BUFFER];
  char output_full_name[MAX_BUFFER];
  char summary_name[MAX_BUFFER];
  char summary_full_name[MAX_BUFFER];

  PPLATE plate_table = NULL;
  PPLATE tmp_plate_table;
  int plate_nrecs = 0;
  int plate_alloc = 0;
  int plate_index;
  PPLATE pPlateEntry;
  char *charPtr;

  PMASTER master_table = NULL;
  PMASTER tmp_master_table;

  File master_handle = NULL;
  char master_name[MAX_BUFFER];
  TableHead master_header = NULL;
  int master_nrecs = 0;
  int master_alloc = 0;
  int master_index;
  PMASTER pMasterEntry;
  MASTER master_record;
  PMASTER pMaster = &master_record;
  TblDescriptor master_descriptor = NULL;
  TableRow master_row = NULL;

  PMAGENTRY pMagEntry1;
  PMAGENTRY pMagEntry2;

  int errorFlag = 0;

  File list_handle = NULL;
  char list_name[MAX_BUFFER];
  TableHead list_header = NULL;
  PLISTENTRY list_table = NULL;
  PLISTENTRY tmp_list_table = NULL;
  size_t list_nrecs = 0;
  size_t list_alloc = 0;

  File magiso_handle = NULL;
  char magiso_name[MAX_BUFFER];
  TableHead magiso_header = NULL;
  PMAGISOENTRY magiso_table = NULL;
	PMAGISOENTRY pMagIsoEntry;
  size_t magiso_nrecs = 0;
	size_t magiso_index;

  char okFlag;
  double last_Stdmag;
  int list_index;
  int list_index1;
  int list_index2;
  int list_index3;
  int brightness_index;
  PLISTENTRY pListEntry;
  PLISTENTRY pListEntry1;
  PLISTENTRY pListEntry2;
  PLISTENTRY pListEntry3;
  int compare1good;
  int compare2good;
  int compare3good;
  int workDone;
  LISTENTRY brightness_queue;
  int goodpoints;
  int validPoints;
  int foundPoints;
  int errorPoints;
  int minmag;
  int maxmag;
  int badPoints = 0;
  int wedgeCount = 0;
  int highDradCount = 0;
  int defectCount = 0;
  int hiZoutCount = 0;
  int tooBrightCount = 0;
  int lowAltitudeCount = 0;
  int blendCount = 0;
  int limitingPoints1 = 0;
  int limitingPoints2 = 0;
  int bin9Points = 0;
  int bothPoints = 0;

  double rawmed = 0;
  double rawrms = 0;
  double clipmed = 0;
  double cliprms = 0;
  char lightcurveName[MAX_BUFFER];
  FILE *lightcurveHandle;
  char lightcurveplot[MAX_BUFFER];
  double startJD = 0.0;
  double endJD;
  int sequence;
  double error_bar_factor = 1.0;
  double startYear = 0.0;
  double endYear;
  FILE *outputHandle = NULL;
  FILE *summaryHandle = NULL;
  int fitFileCount = 0;
  int fitFileIndex;;
  struct dirent **fitFileList;
  char cmdStr[MAX_BUFFER];
  int result;

  float *xval = NULL;
  float *yval = NULL;

  float *xmin_threshold = NULL;
  float *ysequencesize = NULL;
  float xthreshold_max;
  int threshold_nrecs = 0;
 
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  int plot_nrecs;
  int plot_index;
  double yearval;
  float *timeval = NULL;
  float *sequenceval = NULL;
  float *magval = NULL;
  int *symbolval = NULL;
  float sequencemin;
  float sequencemax;
  float timemin;
  float timemax;
  float magmin;
  float magmax;
  int curve_nrecs;
  int curve_index;
  int skipMaster = 0;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int entries_printed;


#if 1
  printf("WARNING: ra and dec are not corrected for precession\n");
#endif


  time(&startTime);
  list_name[0] = 0;
  master_name[0] = 0;
  workingDirectory[0] = 0;
  output_name[0] = 0;
  summary_name[0] = 0;
	magiso_name[0] = 0;
  catalogString[0] = 0;
  qualifier[0] = 0;

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


        case 'l': /* list file name */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(list_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'p': /* magiso file name */
        case 'P':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(magiso_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;



        case 'm': /* master file name */
        case 'M':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(master_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'r': /* minimum refcnt */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&min_refcnt);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the min_refcnt %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 't': /* minimum threshold */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&base_min_threshold);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the min_refcnt %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'f': /* fraction of plates required */
        case 'F':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&fraction);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the fraction %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'a': /* analysis threshold */
        case 'A':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&analysis_threshold);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the analysis_threshold %s\n",*argv);
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
              strcpy(qualifier,*argv);
            }
          }
          break;



        case 'd': /* lightcurve output directory */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(workingDirectory,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
            length = strlen(workingDirectory);
            if (workingDirectory[length-1] == '/') {
              workingDirectory[length-1] = 0;
            }
          }
          break;


        case 'v': /* verbose */
        case 'V': /* verbose */
          verbose += 1;
          break;

        case 's': /* skip master */
        case 'S': /* skip master */
          skipMaster = 1;
          break;



        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  if (base_min_threshold >= 0) {
    min_threshold = base_min_threshold;
  }

  /* Verify that we have everything */
  if (workingDirectory[0] == 0) {
    printf("ERROR: No working directory was specified\n");
    errorFlag = 1;
  }

  /* Open the output file */
  strcpy(output_name,"id_estimate.db");
  if (strstr(output_name,"/") != NULL) {
    strcpy(output_full_name,output_name);
  } else {
    strcpy(output_full_name,workingDirectory);
    strcat(output_full_name,"/");
    strcat(output_full_name,output_name);
  }

  outputHandle = fopen(output_full_name,"wt");
  if (outputHandle == NULL) {
    errorFlag = 1;
    printf("ERROR: failed to open the output file %s\n",output_full_name);
  } else {
    if (verbose) {
      printf("Output file is %s\n",output_full_name);
    }
    fprintf(outputHandle,"median_local\tStdmag\trms_local\tclip_median_local\tclip_rms_local\tnpoints\tngood\tyrbegin\tyrend\tREF\n");
    fprintf(outputHandle,"------------\t------\t---------\t-----------------\t--------------\t-------\t-----\t-------\t-----\t---\n");

  }
  /* Open the summary file */
  strcpy(summary_name,"summary_estimate.db");
  if (strstr(summary_name,"/") != NULL) {
    strcpy(summary_full_name,summary_name);
  } else {
    strcpy(summary_full_name,workingDirectory);
    strcat(summary_full_name,"/");
    strcat(summary_full_name,summary_name);
  }

  summaryHandle = fopen(summary_full_name,"wt");
  if (summaryHandle == NULL) {
    errorFlag = 1;
    printf("ERROR: failed to open the summary file %s\n",summary_full_name);
  } else {
    if (verbose) {
      printf("Summary file is %s\n",summary_full_name);
    }
    fprintf(summaryHandle,"minmag\tmaxmag\tcount\tmedian\n");
    fprintf(summaryHandle,"------\t------\t-----\t------\n");

  }
  /* Open the list file */
  if (list_name[0] != 0) {
    list_handle = Open(list_name,"r");
    if (list_handle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to find the list file %s\n",list_name);
    } else {
      if (verbose) {
        printf("Found list file %s\n",list_name);
      }
    }
  }



  /* Now search the directory for fits files and build a list */
  fitFileCount = scandir(workingDirectory,&fitFileList,fitsfilter,alphasort);
  if (fitFileCount <= 0) {
    printf("ERROR %d accessing fits files in %s\n",fitFileCount,workingDirectory);
    errorFlag = 1;
  } else {
    plate_table = (PPLATE)calloc(fitFileCount,sizeof(PLATE));
    plate_nrecs = fitFileCount;

    if (plate_table == NULL) {
      printf("ERROR allocating plate_table of size %d\n",fitFileCount);
      errorFlag = 1;
    } else {
      for (fitFileIndex = 0; fitFileIndex < fitFileCount; fitFileIndex++) {
        if (strlen(fitFileList[fitFileIndex]->d_name) > (MAX_PLATE_NAME-1)) {
          printf("ERROR: %s is greater than MAX_PLATE_NAME\n",fitFileList[fitFileIndex]->d_name);
          errorFlag = 1;
        } else {
          pPlateEntry = &plate_table[fitFileIndex];
          strcpy(pPlateEntry->fitsFileName,fitFileList[fitFileIndex]->d_name);
          charPtr = strstr(pPlateEntry->fitsFileName,".fit");
          if (charPtr != NULL) {
            *charPtr = 0;
          }
          /* Search forward to the first a-z character which will be the series name */
          charPtr = pPlateEntry->fitsFileName;
          while (*charPtr != 0) {
            if ((*charPtr >= 'a') && (*charPtr <= 'z')) {
              break;
            }
            charPtr++;
          }
          if (*charPtr == 0) {
            printf("ERROR obtaining the root from %s\n",pPlateEntry->fitsFileName);
            errorFlag = 1;
          } else {
            strcpy(pPlateEntry->Plate,charPtr);
          }
          for (list_index = 0; list_index < MAX_LIST_INDEX; list_index++) {
            pPlateEntry->magentry[list_index].master_index = -1;
            pPlateEntry->magentry[list_index].list_index = list_index;
          }

        }
        free(fitFileList[fitFileIndex]); 
      }
      free(fitFileList);
    }
  }
  if (errorFlag) {
    printf("Usage: estimate_lightcurves -l <comparison star list (deprecated)> \n");
    printf("                            -m <master_all file (deprecated)> \n");
    printf("                            -p <limiting MAG_ISO_med file name (deprecated)> \n");
    printf("                            -d <working directory > \n");
    printf("                            -q input catalog qualifer\n");
    printf("                            -r <minimum refcnt    (default = %d)>\n",DEFAULT_REFCNT);
    printf("                            -t <minimum threshold (default = %d)>\n",DEFAULT_THRESHOLD);
    printf("                            -f <fraction of plates brighter than (default = %.2f)>\n",DEFAULT_FRACTION);
    printf("                            -a <analysis_threshold - minimum FLUX_MAX/THRESHOLD (default = %.2f)>\n",DEFAULT_ANALYSIS_THRESHOLD);
    printf("                           [-v]  verbose\n");
 
    return(-1);
  }
  
  printf("estimate_lightcurves of %s %s for %d fits files min_refcnt: %d min_threshold: %d fraction %.2f analysis_threshold %.2f\n",__DATE__,__TIME__,fitFileCount,min_refcnt,min_threshold,fraction,analysis_threshold);

  /* Open our plot */
  strcpy(plotdev,workingDirectory);
  strcat(plotdev,"/sequence_plot.ps/ps");
  cpgopen(plotdev);
  cpgslw(4); /* line width is 4/200 inch */
  cpgsch(1.5); /* Character size is 1.5/40 the height of the view surface */
  /* subdivide the page into 4 panels;*/
  cpgsubp(2,2);



  /* run estimate_lightcurves.csh */
  strcpy(plateListFile,workingDirectory);
  strcat(plateListFile,"/platelist.xxx");
  plateListHandle = fopen(plateListFile,"wt");
  if (plateListHandle == NULL) {
    printf("ERROR: failed to open %s\n",plateListFile);
    exit(-1);
  }
  for (fitFileIndex = 0; fitFileIndex < fitFileCount; fitFileIndex++) {
    pPlateEntry = &plate_table[fitFileIndex];
    fprintf(plateListHandle,"%s%%%s\n",pPlateEntry->fitsFileName,pPlateEntry->Plate);
  }
  fclose(plateListHandle);
  if (skipMaster == 0) {
    /* Now create master.db, sequence.db and magiso.db */
    sprintf(cmdStr,"$DASCH_SCRIPTS/estimate_lightcurves.csh %s %s %s >& %s/estimate.log\n",plateListFile,workingDirectory,qualifier,workingDirectory);
    printf("Executing %s\n",cmdStr);
    result = system(cmdStr);
    if (result != 0) {
      printf("ERROR: result %d for %s\n",result,cmdStr);
    }
  } else {
    printf("WARNING: skipMaster is set\n");
  }
  /* Open the magiso file */
  if (magiso_name[0] == 0) {
    sprintf(magiso_name,"%s/magiso.db",workingDirectory);
  }



  magiso_handle = Open(magiso_name,"r");
  if (magiso_handle == NULL) {
    printf("ERROR: Failed to find the magiso file %s\n",magiso_name);
    exit(-1);
  } else {
    if (verbose) {
      printf("Found magiso file %s\n",magiso_name);
    }
  }



	/* Now read in the magiso file */
  magiso_header = table_header(magiso_handle,TABLE_PARSE);
  if (magiso_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",magiso_name);
    return(-1);
  }

  magiso_table = table_loadva(magiso_handle,
															&magiso_header,
															NULL, /* hbase */
															NULL, /* rows */
															NULL,
															sizeof(MAGISOENTRY),
															&magiso_nrecs,
															TblDbl,"MAG_ISO_med",TblOff(PMAGISOENTRY,MAG_ISO_med),
															TblDbl,"MAG_ISO_rms",TblOff(PMAGISOENTRY,MAG_ISO_rms),
															TblBuf,"Plate"    ,TblOff(PMAGISOENTRY,Plate),MAX_PLATE_NAME, 
															0,"end",0);
  if (magiso_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",magiso_name);
    return(-1);
  }


  /* Now read in the list file */
  if (list_handle != NULL) {
    list_header = table_header(list_handle,TABLE_PARSE);
    if (list_header == NULL) {
      printf("ERROR: Failed to read header for %s\n",list_name);
      return(-1);
    }

    list_table = table_loadva(list_handle,
                              &list_header,
                              NULL, /* hbase */
                              NULL, /* rows */
                              NULL,
                              sizeof(LISTENTRY),
                              &list_nrecs,
                              TblDbl,"ra",TblOff(PLISTENTRY,ra),
                              TblDbl,"dec",TblOff(PLISTENTRY,dec),
                              TblDbl,"Stdmag",TblOff(PLISTENTRY,Stdmag),
                              TblBuf,"src_name"    ,TblOff(PLISTENTRY,src_name),MAX_SRC_LENGTH,
                              TblBuf,"REF"    ,TblOff(PLISTENTRY,REF),MAX_REF, 
                              0,"end",0);
    if (list_table == NULL) {
      printf("ERROR: Failed to read table for %s\n",list_name);
      return(-1);
    }
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      pListEntry = &list_table[list_index];
      pListEntry->list_index = list_index;
      pListEntry->flink = NULL;
      pListEntry->brightnessIndex = -1;
      pListEntry->refcnt = 0;
    }

    if (verbose) {
      printf("read %d records for %s\n",list_nrecs,list_name);
    }
    if (list_nrecs > MAX_LIST_INDEX) {
      printf("ERROR: MAX_LIST_INDEX should be at least %d\n",list_nrecs);
      exit(-1);
    }

  }

  /* Open the master photometry file */
  if (master_name[0] == 0) {
    sprintf(master_name,"%s/master.db",workingDirectory);    
  }

  master_handle = Open(master_name,"r");
  if (master_handle == NULL) {
    printf("ERROR: Failed to find the master catalog file %s\n",master_name);
    exit(-1);
  } else {
    if (verbose) {
      printf("Found master catalog file %s\n",master_name);
    }
  }


  /* Read the header of the master table */
  master_header = table_header(master_handle,TABLE_PARSE);
  if (master_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",master_name);
    return(-1);
  }


  master_descriptor = table_create_descrip(&master_nrecs,
                                           TblDbl,"FLUX_ISO",TblOff(PMASTER,FLUX_ISO),
                                           TblDbl,"FLUX_MAX",TblOff(PMASTER,FLUX_MAX),
                                           TblDbl,"THRESHOLD",TblOff(PMASTER,THRESHOLD),
                                           TblDbl,"Date",TblOff(PMASTER,Date),
                                           TblDbl,"ra",TblOff(PMASTER,ra),
                                           TblDbl,"dec",TblOff(PMASTER,dec),
                                           TblDbl,"Stdmag",TblOff(PMASTER,Stdmag),
                                           TblInt,"spatial_bin",TblOff(PMASTER,spatial_bin),
                                           TblInt,"gsc_bin_index",TblOff(PMASTER,gsc_bin_index),
                                           TblInt,"AFLAGS",TblOff(PMASTER,AFLAGS),
                                           TblInt,"BFLAGS",TblOff(PMASTER,BFLAGS),
                                           TblBuf,"REF"    ,TblOff(PMASTER,REF),MAX_REF,
                                           TblBuf,"Plate"    ,TblOff(PMASTER,Plate),MAX_PLATE_NAME,
                                           0,"end",0);

  if (master_descriptor == NULL) {
    printf("ERROR: Failed to allocate descriptor for %s\n",master_name);
    return(-1);
  }
  table_loadmap(master_header,master_descriptor);


  if (verbose) {
    printf("Reading %s\n",master_name);
  }
  memset(pMaster,0,sizeof(MASTER));
  if (verbose > 1) {
    debugFlag = 1;
  }
  /* Now scan through the master table and extract entries for the stars that we are interested in */
  while (1) {

    if(ReadMaster(pMaster,master_handle,master_header,master_descriptor,&master_row) == 0) {
      break;
    }
    if ((pMaster->AFLAGS & (1 << FILTER_AFLAG_BLEND_NOMATCH)) != 0) {
      continue;
    }

    if (list_handle == NULL) {
      /* No list provided, so save everything */
      if (master_nrecs >= master_alloc) {
        master_alloc += 1000;
        tmp_master_table = realloc(master_table,master_alloc * sizeof(MASTER));
        if (tmp_master_table == NULL) {
          printf("ERROR: failed to reallocate the master table with size %d\n",master_alloc * sizeof(MASTER));
          exit(-1);
        }
        master_table = tmp_master_table;
        tmp_master_table = NULL;
      }
      pMaster->list_index = -1;
      memcpy(&master_table[master_nrecs],pMaster,sizeof(MASTER));
      master_nrecs++;
      
    } else {
      for (list_index = 0; list_index < list_nrecs; list_index++) {
        pListEntry = &list_table[list_index];
        if (strcmp(pListEntry->REF,pMaster->REF) == 0) {
          /* We want to save this entry */
          if (master_nrecs >= master_alloc) {
            master_alloc += 1000;
            tmp_master_table = realloc(master_table,master_alloc * sizeof(MASTER));
            if (tmp_master_table == NULL) {
              printf("ERROR: failed to reallocate the master table with size %d\n",master_alloc * sizeof(MASTER));
              exit(-1);
            }
            master_table = tmp_master_table;
            tmp_master_table = NULL;
          }
          pMaster->list_index = list_index;
          memcpy(&master_table[master_nrecs],pMaster,sizeof(MASTER));
          pListEntry->refcnt++;
          master_nrecs++;
          break;
        }
      }
    }

  }
  if (verbose) {
    printf("found %d records in %s\n",master_nrecs,master_name);
  }
  if (list_handle == NULL) {
    /* Allocate the list table using every reference found in the master table */
    for (master_index = 0; master_index < master_nrecs; master_index++) {
      pMasterEntry = &master_table[master_index];
      for (list_index = 0; list_index < list_nrecs; list_index++) {
        pListEntry = &list_table[list_index];
        if (strcmp(pListEntry->REF,pMasterEntry->REF) == 0) {
          /* This gives us the most recent ra and dec */
          pListEntry->ra = pMasterEntry->ra;
          pListEntry->dec = pMasterEntry->dec;
          break;
        }
      }
      if (list_index >= list_nrecs) {
        if (list_nrecs >= list_alloc) {
          list_alloc += 1000;
          tmp_list_table = realloc(list_table,list_alloc * sizeof(LISTENTRY));
          if (tmp_list_table == NULL) {
            printf("ERROR: failed to reallocate the list table with size %d\n",list_alloc * sizeof(LISTENTRY));
            exit(-1);
          }
          list_table = tmp_list_table;
          tmp_list_table = NULL;
        }
        pListEntry = &list_table[list_index];
        strcpy(pListEntry->REF,pMasterEntry->REF);
        strcpy(pListEntry->src_name,pMasterEntry->REF);
        pListEntry->ra = pMasterEntry->ra;
        pListEntry->dec = pMasterEntry->dec;
        pListEntry->Stdmag = pMasterEntry->Stdmag;
        pListEntry->gsc_bin_index = pMasterEntry->gsc_bin_index;
        list_nrecs++;
        if (list_nrecs >= (MAX_LIST_INDEX-1)) {
          printf("ERROR: MAX_LIST_INDEX exceeded\n");
          exit(-1);
        }
      }
    }
  }
  if (list_nrecs <= 0) {
    printf("ERROR: No star designations found in the master file\n");
    exit(-1);
  }
  if (verbose) {
    printf("found %d star designations in the master file\n",list_nrecs);
  }
  vector1 = (double *)calloc(list_nrecs,sizeof(double));
  /* Sort the list in terms of decreasing brightness */
  if (list_handle == NULL) {
    qsort((void *)list_table,list_nrecs,sizeof(LISTENTRY),ListCompare);
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      pListEntry = &list_table[list_index];
      pListEntry->list_index = list_index;
      pListEntry->flink = NULL;
      pListEntry->brightnessIndex = -1;
      pListEntry->refcnt = 0;
      pListEntry->rawmed = 99.0;
      pListEntry->rawrms = 99.0;
    }
  }

#if 0
  printf("AFTER SORT\n");
  for (list_index = 0; list_index < list_nrecs; list_index++) {
    pListEntry = &list_table[list_index];
    printf("List index %d REF %s\n",list_index,pListEntry->REF);
  }
#endif

  /* Now search through the master records and allocate the plates table */
  for (master_index = 0; master_index < master_nrecs; master_index++) {
    pMasterEntry = &master_table[master_index];
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      pListEntry = &list_table[list_index];
#if 0
      if (strcmp(pListEntry->REF,"N23301237043") == 0) {
        printf("at %s\n",pListEntry->REF);
      }
#endif
      if (strcmp(pListEntry->REF,pMasterEntry->REF) == 0) {
        pMasterEntry->list_index = list_index;
        pListEntry->refcnt++;
        break;
      }
    }
    if (list_index == list_nrecs) {
      printf("ERROR: list consistency check\n");
      exit(-1);
    }

    for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
      pPlateEntry = &plate_table[plate_index];
      if (strcmp(pMasterEntry->Plate,pPlateEntry->Plate) == 0) {
        double dateerr = pMasterEntry->Date - pPlateEntry->Date;
        if (dateerr < 0) {
          dateerr = - dateerr;
        }
        if (fitFileCount != 0) {
          pPlateEntry->Date = pMasterEntry->Date;
        } else {

          /* Allow for 18 minutes of solar system travel time */
          if (dateerr > 0.000208) {
            printf("ERROR: Julian dates %f and %f do not agree for plate %s\n",pMasterEntry->Date,pPlateEntry->Date,pPlateEntry->Plate);
            exit(-1);
          }
        }
        if (pMasterEntry->list_index < 0) {
          printf("ERROR: master entry list index consistency check\n");
          exit(-1);
        }
        pPlateEntry->magentry[pMasterEntry->list_index].master_index = master_index;
        pPlateEntry->magentry[pMasterEntry->list_index].FLUX_ISO = pMasterEntry->FLUX_ISO;
        pPlateEntry->magentry[pMasterEntry->list_index].FLUX_MAX = pMasterEntry->FLUX_MAX;
        pPlateEntry->magentry[pMasterEntry->list_index].THRESHOLD = pMasterEntry->THRESHOLD;
        pPlateEntry->magentry[pMasterEntry->list_index].list_index = pMasterEntry->list_index;
        pPlateEntry->magentry[pMasterEntry->list_index].Date = pMasterEntry->Date;
        pPlateEntry->magentry[pMasterEntry->list_index].spatial_bin = pMasterEntry->spatial_bin;
        pPlateEntry->magentry[pMasterEntry->list_index].AFLAGS = pMasterEntry->AFLAGS;
        break;
      }
    }
    if (plate_index >= plate_nrecs) {
      if (fitFileCount != 0) {
        printf("ERROR: plate allocation consistency check for %s\n",pMasterEntry->Plate);

        for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
          pPlateEntry = &plate_table[plate_index];
          printf("Index %d Plate %s\n",plate_index,pPlateEntry->Plate);
        }
        exit(-1);
      }
      if (plate_nrecs >= plate_alloc) {
        plate_alloc += 100;
        tmp_plate_table = realloc(plate_table,plate_alloc * sizeof(PLATE));
        if (tmp_plate_table == NULL) {
          printf("ERROR: failed to reallocate the plate table with entries %d size %d\n",plate_alloc,plate_alloc * sizeof(PLATE));
          exit(-1);
        }
        plate_table = tmp_plate_table;
        tmp_plate_table = NULL;
      }
      pPlateEntry = &plate_table[plate_nrecs];
      plate_nrecs++;
      memset(pPlateEntry,0,sizeof(PLATE));
      for (list_index = 0; list_index < MAX_LIST_INDEX; list_index++) {
        pPlateEntry->magentry[list_index].master_index = -1;
        pPlateEntry->magentry[list_index].list_index = list_index;
        
      }
      pPlateEntry->magentry[pMasterEntry->list_index].master_index = master_index;
      pPlateEntry->magentry[pMasterEntry->list_index].FLUX_ISO = pMasterEntry->FLUX_ISO;
      pPlateEntry->magentry[pMasterEntry->list_index].FLUX_MAX = pMasterEntry->FLUX_MAX;
      pPlateEntry->magentry[pMasterEntry->list_index].THRESHOLD = pMasterEntry->THRESHOLD;
      pPlateEntry->magentry[pMasterEntry->list_index].list_index = pMasterEntry->list_index;
      pPlateEntry->magentry[pMasterEntry->list_index].Date = pMasterEntry->Date;
      pPlateEntry->magentry[pMasterEntry->list_index].spatial_bin = pMasterEntry->spatial_bin;
      pPlateEntry->magentry[pMasterEntry->list_index].AFLAGS = pMasterEntry->AFLAGS;
      strcpy(pPlateEntry->Plate,pMasterEntry->Plate);
      pPlateEntry->Date = pMasterEntry->Date;

    }
  }
  if (verbose) {
    printf("found %d plates in %s\n",plate_nrecs,master_name);
  }
  /* Sort the plate table in order of increasing Julian Day */
  qsort((void*)plate_table,plate_nrecs,sizeof(PLATE),PlateCompare);
  vector = (double *)calloc(plate_nrecs,sizeof(double));
  if (vector == NULL) {
    printf("ERROR: failed to allocate the plate vector\n");
    exit(-1);
  }

  

  /* Next we need to establish the order of brightness.  Sort all of the magnitude entries
     in increasing FLUX_ISO */
  
  for (list_index = 0; list_index < MAX_LIST_INDEX; list_index++) {
    brightnessOrder[list_index] = -1;
    for (list_index2 = 0; list_index2 < MAX_LIST_INDEX; list_index2++) {
      
      lowestCount[list_index][list_index2] = 0;
    }
  }

  for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
    pPlateEntry = &plate_table[plate_index];
    if (pPlateEntry->Date == 0) {
      continue;
    }
#if 0
    printf("index %d JulianDate %f plate %s\n",plate_index,pPlateEntry->Date,pPlateEntry->Plate);
#endif
    memcpy(&pPlateEntry->magentrySort,&pPlateEntry->magentry,sizeof(pPlateEntry->magentry));
    qsort((void*)&pPlateEntry->magentrySort[0],list_nrecs,sizeof(MAGENTRY),MagEntryCompare);
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      pMagEntry1 = &pPlateEntry->magentrySort[list_index];
#if 0
      if (pMagEntry1->master_index >= 0) {
        printf("list_index %d %d master index %d FLUX_ISO %f Date %f\n",
               list_index,pMagEntry1->list_index,pMagEntry1->master_index,pMagEntry1->FLUX_ISO,pMagEntry1->Date);
      }
    
#endif
      for (list_index2 = list_index+1; list_index2 < list_nrecs; list_index2++) { 
        pMagEntry2 = &pPlateEntry->magentrySort[list_index2];
        if ((pMagEntry1->master_index >= 0) &&
            ((pMagEntry1->FLUX_MAX/pMagEntry1->THRESHOLD) >= analysis_threshold) &&
            (pMagEntry2->master_index >= 0) &&
            ((pMagEntry2->FLUX_MAX/pMagEntry2->THRESHOLD) >= analysis_threshold)) {
          /* pMagEntry1 is the brighter object and pMagEntry2 is the dimmer object */
          lowestCount[pMagEntry1->list_index][pMagEntry2->list_index]++;
        }
      }

    }

  }

#if 0


  for (list_index2 = 0; list_index2 < list_nrecs; list_index2++) {
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      if (lowestCount[list_index][list_index2] > 0) {
        printf("%d %d %d %d\n",list_index2,list_index,lowestCount[list_index][list_index2],lowestCount[list_index][list_index2]+lowestCount[list_index2][list_index]);
      }
    }
  }
  printf("\n");


#endif

  xval = (float *)calloc(list_nrecs,sizeof(float));
  yval = (float *)calloc(list_nrecs,sizeof(float));
  xmin_threshold = (float *)calloc(list_nrecs,sizeof(float));
  ysequencesize = (float *)calloc(list_nrecs,sizeof(float));
  

  /* Now establish a queue of brighness from the brightest to the dimmest based on the
     frequencies found above.  If we do not have a bas_min_threshold, we will loop until
     we find the queue with the longest length */
  if (base_min_threshold <= 0) {
    /* Iterate over all threshold values */
    cur_min_threshold = list_nrecs;
  } else {
    /* Only one iteration if we are using a fixed threshold */
    cur_min_threshold = 0;
  }
  while (cur_min_threshold >= 0) {
    if (cur_min_threshold > 0) {
      min_threshold = cur_min_threshold;
    } else {
      if (base_min_threshold > 0) {
        min_threshold = base_min_threshold;
      } else {
        min_threshold = best_min_threshold;
      }
    }
    memset(&brightness_queue,0,sizeof(LISTENTRY));

    for (list_index1 = 0; list_index1 < list_nrecs; list_index1++) {
      pListEntry = &list_table[list_index1];
      pListEntry->queued = 0;
      pListEntry->flink = NULL;
    }



    for (list_index1 = 0; list_index1 < list_nrecs; list_index1++) {
      pListEntry = &list_table[list_index1];
      if ((pListEntry->refcnt < min_refcnt) ||
          (pListEntry->Stdmag > 90.0)) {
        pListEntry->queued = -1;
        continue;
      }
      pListEntry1 = &brightness_queue;
      if (pListEntry1->flink == NULL) {
        pListEntry1->flink = pListEntry;
      } else {
        pListEntry2 = pListEntry1->flink;
        while (pListEntry2 != NULL) {
          if ((lowestCount[pListEntry->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry->list_index]) >= min_threshold) {
            fraction1 = (1.0*lowestCount[pListEntry->list_index][pListEntry2->list_index])/(1.0*(lowestCount[pListEntry->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry->list_index]));
          } else {
            fraction1 = 0.0;
          }
          if (fraction1 >= fraction) {
            /* The next entry is usually dimmer than the current.  Insert our  entry here */
            pListEntry->queued = 1;
            pListEntry->flink = pListEntry2;
            pListEntry1->flink = pListEntry;
            pListEntry = NULL;
            break;
          }
          pListEntry1 = pListEntry2;
          pListEntry2 = pListEntry1->flink;
        }
        if (pListEntry != NULL) {
          /* Not yet inserted, so it must be the dimmest */
          pListEntry1->flink = pListEntry;
          pListEntry->queued = 1;
        }
      }
    }
    /* Plot results */
    plot_nrecs = 0;
    pListEntry = brightness_queue.flink;
    while (pListEntry != NULL) {
      xval[plot_nrecs] = plot_nrecs;
      yval[plot_nrecs] = pListEntry->Stdmag;
      plot_nrecs++;
      pListEntry = pListEntry->flink;
    }
    if (plot_nrecs >= 2) {
      /* Plot this chain */

      xmin = xval[0];
      xmax = xval[0];
      ymin = yval[0];
      ymax = yval[0];
      for (plot_index = 1; plot_index < plot_nrecs; plot_index++) {
        if (xmin > xval[plot_index]) {
          xmin = xval[plot_index];
        }
        if (xmax < xval[plot_index]) {
          xmax = xval[plot_index];
        }
        if (ymin > yval[plot_index]) {
          ymin = yval[plot_index];
        }
        if (ymax < yval[plot_index]) {
          ymax = yval[plot_index];
        }
      }
      if (cur_min_threshold == 0) {
        /* cpgpanl(1,1); */
        cpgenv(xmin,xmax,ymax,ymin,0,0);
        cpglab("Sequence","Stdmag","Comparison Star Sequence 1");
        cpgpt(plot_nrecs,xval,yval,4);
      }
    }
    /* Now chain through the list and verify the order */
    pListEntry = brightness_queue.flink;
    while (pListEntry != NULL) {
      pListEntry2 = pListEntry->flink;
      if (pListEntry2 != NULL) {
        if ((lowestCount[pListEntry->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry->list_index]) < min_threshold) {
#if 0
          printf("ERROR: uncertain comparison for %d %s and %d %s\n",pListEntry->list_index,pListEntry->REF,pListEntry2->list_index,pListEntry2->REF);
#endif
          /* Unhook this one */
          pListEntry2->queued = -2;
          pListEntry->flink = pListEntry2->flink;
          /* And retry the test */
          pListEntry2 = pListEntry;
        } else {
          if ((lowestCount[pListEntry->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry->list_index]) >= min_threshold) {
            fraction1 = (1.0*lowestCount[pListEntry->list_index][pListEntry2->list_index])/(1.0*(lowestCount[pListEntry->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry->list_index]));
          } else {
            fraction1 = 9000.0;
          }
          if (fraction1 < fraction) {
#if 0
            printf("ERROR: uncertain comparison for %d %s and %d %s\n",pListEntry->list_index,pListEntry->REF,pListEntry2->list_index,pListEntry2->REF);
#endif
          }

        }

      }
      pListEntry = pListEntry2;
    }
    /* Plot results */
    plot_nrecs = 0;
    pListEntry = brightness_queue.flink;
    while (pListEntry != NULL) {
      xval[plot_nrecs] = plot_nrecs;
      yval[plot_nrecs] = pListEntry->Stdmag;
      plot_nrecs++;
      pListEntry = pListEntry->flink;
    }
    if (plot_nrecs >= 2) {
      /* Plot this chain */

      xmin = xval[0];
      xmax = xval[0];
      ymin = yval[0];
      ymax = yval[0];
      for (plot_index = 1; plot_index < plot_nrecs; plot_index++) {
        if (xmin > xval[plot_index]) {
          xmin = xval[plot_index];
        }
        if (xmax < xval[plot_index]) {
          xmax = xval[plot_index];
        }
        if (ymin > yval[plot_index]) {
          ymin = yval[plot_index];
        }
        if (ymax < yval[plot_index]) {
          ymax = yval[plot_index];
        }
      }
      if (cur_min_threshold == 0) {
        cpgpanl(1,1);
        cpgenv(xmin,xmax,ymax,ymin,0,0);
        cpglab("Sequence","Stdmag","Comparison Star Sequence 2");
        cpgpt(plot_nrecs,xval,yval,4);
      }
    }
    /* Now go throught the list and eliminate anything that is out of order */
    pListEntry = NULL;
    pListEntry1 = brightness_queue.flink;
    while (pListEntry1 != NULL) {
      compare1good = 0;
      compare2good = 0;
      compare3good = 0;
      pListEntry2 = pListEntry1->flink;
      if (pListEntry2 == NULL) {
        break;
      }
      pListEntry3 = pListEntry2->flink;
      if ((lowestCount[pListEntry1->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry1->list_index]) >= min_threshold) {
        fraction1 = (1.0*lowestCount[pListEntry1->list_index][pListEntry2->list_index])/(1.0*(lowestCount[pListEntry1->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry1->list_index]));
      } else {
        fraction1 = 0.0;
      }
      if ((fraction1 >= fraction) &&
          (pListEntry1->Stdmag < pListEntry2->Stdmag)) {
        compare1good = 1; /* The first two stars compare fine */
      }
      if (pListEntry3 != NULL) {
        if ((lowestCount[pListEntry2->list_index][pListEntry3->list_index]+lowestCount[pListEntry3->list_index][pListEntry2->list_index]) >= min_threshold) {
          fraction1 = (1.0*lowestCount[pListEntry2->list_index][pListEntry3->list_index])/(1.0*(lowestCount[pListEntry2->list_index][pListEntry3->list_index]+lowestCount[pListEntry3->list_index][pListEntry2->list_index]));
        } else {
          fraction1 = 0.0;
        }
        if ((pListEntry3 != NULL) &&
            (fraction1 >= fraction) &&
            (pListEntry2->Stdmag < pListEntry3->Stdmag)) {
          compare2good = 1; /* The next two stars compare fine */
        }
      }
      if (pListEntry3 != NULL) {
        if ((lowestCount[pListEntry1->list_index][pListEntry3->list_index]+lowestCount[pListEntry3->list_index][pListEntry1->list_index]) >= min_threshold) {
          fraction1 = (1.0*lowestCount[pListEntry1->list_index][pListEntry3->list_index])/(1.0*(lowestCount[pListEntry1->list_index][pListEntry3->list_index]+lowestCount[pListEntry3->list_index][pListEntry1->list_index]));
        } else {
          fraction1 = 0.0;
        }

        if ((pListEntry3 != NULL) &&
            (fraction1 >= fraction) &&
            (pListEntry1->Stdmag < pListEntry3->Stdmag)) {
          compare3good = 1; /* The first and last compare fine */
        }   
      }
      /*  truth table  1   2    3  
          0   0    0  Remove second
          0   0    1  Remove second
          0   1    0  Remove second
          1   0    0  Remove third
          1   0    1  Remove second
          1   1    0  Not possible
          1   1    1  O.K
      */
      if ((compare1good == 1) && (compare2good == 0) && (compare3good == 0)) {
        if (pListEntry3 != NULL) {
          pListEntry3->queued = -2;
          pListEntry3 = pListEntry3->flink;
          pListEntry2->flink = pListEntry3;
        } 
      } else if ((compare1good != 1) || (compare2good) != 1) {
        pListEntry2->queued = -2;
        pListEntry2 = pListEntry2->flink;
        pListEntry1->flink = pListEntry2;
        pListEntry2 = pListEntry1;
      }
      pListEntry1 = pListEntry2;
    }
    /* Plot results */
    plot_nrecs = 0;
    pListEntry = brightness_queue.flink;
    while (pListEntry != NULL) {
      xval[plot_nrecs] = plot_nrecs;
      yval[plot_nrecs] = pListEntry->Stdmag;
      plot_nrecs++;
      pListEntry = pListEntry->flink;
    }
    if (plot_nrecs >= 2) {
      /* Plot this chain */

      xmin = xval[0];
      xmax = xval[0];
      ymin = yval[0];
      ymax = yval[0];
      for (plot_index = 1; plot_index < plot_nrecs; plot_index++) {
        if (xmin > xval[plot_index]) {
          xmin = xval[plot_index];
        }
        if (xmax < xval[plot_index]) {
          xmax = xval[plot_index];
        }
        if (ymin > yval[plot_index]) {
          ymin = yval[plot_index];
        }
        if (ymax < yval[plot_index]) {
          ymax = yval[plot_index];
        }
      }
      if (cur_min_threshold == 0) {

        cpgpanl(2,1);
        cpgenv(xmin,xmax,ymax,ymin,0,0);
        cpglab("Sequence","Stdmag","Comparison Star Sequence 3");
        cpgpt(plot_nrecs,xval,yval,4);
      }
    }

    /* Attempt a second insert.  We are particularly interested in tacking entries onto the end */
    workDone = 1;
    while (workDone) {
      workDone = 0;
      for (list_index = 0; list_index < list_nrecs; list_index++) {
        pListEntry2 = &list_table[list_index];
        if ((pListEntry2->refcnt < min_refcnt) ||
            (pListEntry2->Stdmag > 90.0)) {
          continue;
        }
        if (pListEntry2->queued != -2) {
          continue;
        }
        pListEntry1 = brightness_queue.flink;
        while (pListEntry1 != NULL) {
          compare1good = 0;
          compare2good = 0;
          pListEntry3 = pListEntry1->flink;
          if (pListEntry3 == NULL) {
            break;
          }
          if ((pListEntry3 == pListEntry2) ||
              (pListEntry1 == pListEntry2)) {
            fprintf(stderr,"ERROR: consistency check!\n");
            exit(-1);
          }
          if ((lowestCount[pListEntry1->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry1->list_index]) >= min_threshold) {
            fraction1 = (1.0*lowestCount[pListEntry1->list_index][pListEntry2->list_index])/(1.0*(lowestCount[pListEntry1->list_index][pListEntry2->list_index]+lowestCount[pListEntry2->list_index][pListEntry1->list_index]));
          } else {
            fraction1 = 0.0;
          }
          if ((fraction1 >= fraction) &&
              (pListEntry1->Stdmag < pListEntry2->Stdmag)) {
            compare1good = 1; /* The first two stars compare fine */
          }
          if (pListEntry3 != NULL) {
            if ((lowestCount[pListEntry2->list_index][pListEntry3->list_index]+lowestCount[pListEntry3->list_index][pListEntry2->list_index]) >= min_threshold) {
              fraction1 = (1.0*lowestCount[pListEntry2->list_index][pListEntry3->list_index])/(1.0*(lowestCount[pListEntry2->list_index][pListEntry3->list_index]+lowestCount[pListEntry3->list_index][pListEntry2->list_index]));
            } else {
              fraction1 = 0.0;
            }
            if ((pListEntry3 != NULL) &&
                (fraction1 >= fraction) &&
                (lowestCount[pListEntry2->list_index][pListEntry3->list_index] >= lowestCount[pListEntry3->list_index][pListEntry2->list_index]) &&
                (pListEntry2->Stdmag < pListEntry3->Stdmag)) {
              compare2good = 1; /* The next two stars compare fine */
            }   
          }
          if ((compare1good == 1) && (compare2good == 1) ||
              (compare1good == 1) && (pListEntry3 == NULL)) {
            /* We can insert this one back into the queue */
            pListEntry2->queued = 1;
            pListEntry1->flink = pListEntry2;
            pListEntry2->flink = pListEntry3;
            workDone = 1;
            break;
          }
          pListEntry1 = pListEntry1->flink;
        }
        if (workDone) {
          break;
        }
      }
     
    }
    /* Plot results */
    plot_nrecs = 0;
    pListEntry = brightness_queue.flink;
    while (pListEntry != NULL) {
      xval[plot_nrecs] = plot_nrecs;
      yval[plot_nrecs] = pListEntry->Stdmag;
      plot_nrecs++;
      pListEntry = pListEntry->flink;
    }
    if (plot_nrecs >= 2) {
      /* Plot this chain */

      xmin = xval[0];
      xmax = xval[0];
      ymin = yval[0];
      ymax = yval[0];
      for (plot_index = 1; plot_index < plot_nrecs; plot_index++) {
        if (xmin > xval[plot_index]) {
          xmin = xval[plot_index];
        }
        if (xmax < xval[plot_index]) {
          xmax = xval[plot_index];
        }
        if (ymin > yval[plot_index]) {
          ymin = yval[plot_index];
        }
        if (ymax < yval[plot_index]) {
          ymax = yval[plot_index];
        }
      }
      if (cur_min_threshold == 0) {

        cpgpanl(1,2);
        cpgenv(xmin,xmax,ymax,ymin,0,0);
        cpglab("Sequence","Stdmag","Comparison Star Sequence 4");
        cpgpt(plot_nrecs,xval,yval,4);
      }
    }
#if 0
    if (verbose) {
      printf("min_threshold %d, plot_nrecs %d, best_min_threshold %d, best_selection_count %d\n",min_threshold,plot_nrecs,best_min_threshold,best_selection_count);
    }
#endif
    xmin_threshold[threshold_nrecs] = min_threshold;
    ysequencesize[threshold_nrecs] = plot_nrecs;
    threshold_nrecs++;
    if (cur_min_threshold > 0) {
      if (plot_nrecs > best_selection_count) {
        best_selection_count = plot_nrecs;
        best_min_threshold = min_threshold;
      }
      if (cur_min_threshold <= min_refcnt) {
        cur_min_threshold = 1;
      }
    } 
    cur_min_threshold--;
    
    
  } /* Threshold loop */
    /* Close the plot */
  cpgiden();
  cpgend();
  
  /* Size the threshold plot */
  if ((base_min_threshold <= 0) &&
      (threshold_nrecs > 1)) {
    xthreshold_max = 0.0;
    for (plot_index = 0; plot_index < threshold_nrecs-1; plot_index++) {
#if 0
      printf("plot_index %d, xmin_threshold %f ysequencesize %f xthreshold_max %f\n",plot_index,xmin_threshold[plot_index],ysequencesize[plot_index],xthreshold_max);
#endif
      if (ysequencesize[plot_index] > 1.0) {
        xthreshold_max  = xmin_threshold[plot_index];
        break;
      }
    }
  }



  /* Now chain through the list and set our brightness order */
  pListEntry = brightness_queue.flink;
  list_index = 0;
  list_index2 = 0;
  while (pListEntry != NULL) {
    pListEntry1 = pListEntry->flink;
    list_index = pListEntry->list_index;
    brightnessOrder[list_index2] = list_index;
    pListEntry->brightnessIndex = list_index2;
    list_index2++;
    pListEntry = pListEntry1;
  }
#if 1
  printf("Valid brightness order %d for %d entries\n",list_index2,list_nrecs);
  last_Stdmag = -1000;
  entries_printed = 0;
  for (list_index = 0; list_index < (list_nrecs); list_index++) {
    list_index2 = brightnessOrder[list_index];
    if (list_index2 >= 0) {
      printf("Order %3d index %3d",list_index,list_index2);
      pListEntry = &list_table[list_index2];
      if (pListEntry->Stdmag > last_Stdmag) {
        okFlag = '*';
        last_Stdmag = pListEntry->Stdmag;
      } else {
        okFlag = 'x';
      }
      entries_printed++;
      printf(" Stdmag %5.2f %c  refcnt %3d REF %s ",pListEntry->Stdmag,okFlag,pListEntry->refcnt,pListEntry->REF);
      printf("\n");
    }
  }
  if (entries_printed == 0) {
    printf("NO ENTRIES\n");
  }
  printf("Entries dequeued because of uncertain comparison or Stdmag out of order\n");
  entries_printed = 0;

  for (list_index = 0; list_index < (list_nrecs); list_index++) {
    pListEntry = &list_table[list_index];
    if (pListEntry->queued == -2) {
      entries_printed++;
      printf("index %3d Stdmag %5.2f refcnt %3d REF %s\n",list_index,pListEntry->Stdmag,pListEntry->refcnt,pListEntry->REF);
    }
  }
  if (entries_printed == 0) {
    printf("NO ENTRIES\n");
  }
  printf("Entries rejected because of insufficient plates or bad magnitude\n");
  entries_printed = 0;
  for (list_index = 0; list_index < (list_nrecs); list_index++) {
    pListEntry = &list_table[list_index];
    if (pListEntry->queued == -1) {
      entries_printed++;
      printf("index %3d Stdmag %5.2f refcnt %3d REF %s\n",list_index,pListEntry->Stdmag,pListEntry->refcnt,pListEntry->REF);
    }
  }
  if (entries_printed == 0) {
    printf("NO ENTRIES\n");
  }

#endif
  /* Step through the plate list and find our limiting magnitude, which is determined by the clipped MAG_ISO_med found by update_sextractor from the complete source list */
  for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
    pPlateEntry = &plate_table[plate_index];
    if (pPlateEntry->Date == 0) {
      continue;
    }
    for (magiso_index = 0; magiso_index < magiso_nrecs; magiso_index++) {
      pMagIsoEntry = &magiso_table[magiso_index];
      if (strcmp(pMagIsoEntry->Plate,pPlateEntry->Plate) == 0) {
        pPlateEntry->FLUX_ISO_med = exp10((pMagIsoEntry->MAG_ISO_med+pMagIsoEntry->MAG_ISO_rms)/-2.5);
        break;
      }
			
    }
    if (magiso_index == magiso_nrecs) {
      pPlateEntry->FLUX_ISO_med = 99.0;
      pPlateEntry->limiting_mag_local1 = 99.0;
    }
		
  }
  sequenceval = (float *)calloc(plate_nrecs,sizeof(float));

  timeval = (float *)calloc(plate_nrecs,sizeof(float));
  magval = (float *)calloc(plate_nrecs,sizeof(float));
  symbolval = (int *)calloc(plate_nrecs,sizeof(int));

  /* We can now make our magnitude estimates.  For each star, step through the plate list and do a linear
     interpolation of magnitudes. */
  for (list_index = 0; list_index < list_nrecs; list_index++) {
    goodpoints = 0;
    sequence = 0;
    validPoints = 0;
    errorPoints = 0;
    foundPoints = 0;
    pListEntry = &list_table[list_index];
#if 0
    if (strcmp(pListEntry->src_name,"K11812880") == 0) {
      printf("At  %s\n",pListEntry->src_name);
    }
#endif
    for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
      pPlateEntry = &plate_table[plate_index];
      if (pPlateEntry->Date == 0) {
        continue;
      }
#if 0
      if ((list_index == 0) && (strcmp(pPlateEntry->Plate,"mf35153_00") == 0)) {
        printf("At plate %s\n",pPlateEntry->Plate);
        for (list_index1 = 0; list_index1 < list_nrecs; list_index1++) {
          if (pPlateEntry->magentrySort[list_index1].master_index < 0) {
            /* Not a valid magnitude measurement */
            continue;
          }
          list_index2 = pPlateEntry->magentrySort[list_index1].list_index;
          pListEntry2 = &list_table[list_index2]; 
          printf("%3d %3d %6.2f %s\n",list_index2,list_index1,pListEntry2->Stdmag,pListEntry2->REF);
        }
      }
#endif
      pMagEntry1 = &pPlateEntry->magentry[list_index];
      if (pMagEntry1->master_index < 0) {
        /* Our star was not found on this plate */
        continue;
      }
      foundPoints++;
      /* Build a queue of magnitudes available on this plate, skipping our star of
         interest and skipping any magnitude estimate that is out of order */
      for (list_index1 = 0; list_index1 < list_nrecs; list_index1++) {
        pListEntry1 = &list_table[list_index1];
        pListEntry1->flink = NULL;
        pListEntry1->queued = 0;
      }
      
      pListEntry1 = &brightness_queue; /* Entry1 is our queue */
      pListEntry1->flink = 0;
      brightness_index = 0;
      /* Now step through the official sequence */
      for (list_index1 = 0; list_index1 < list_nrecs; list_index1++) {
        list_index2 = brightnessOrder[list_index1];
        if (list_index2 < 0) {
          continue; /* Not in the sequence */
        }
        if (list_index2 == list_index) {
          continue; /* This is the object of current intrest */
        }
        if (pPlateEntry->magentry[list_index2].master_index < 0) {
          /* Not a valid magnitude measurement on this plate */
          continue;
        }
        /* Put this one in the queue */
        pListEntry2 = &list_table[list_index2]; /* Entry2 is our candidate */
        pListEntry1->flink = pListEntry2;
        pListEntry1 = pListEntry2;
        brightness_index++;
        pListEntry2->queued = 1;
      }
#if 0
      if (strcmp(pListEntry->src_name,"K11812880") == 0) {
        printf("At %s for plate %s\n",pListEntry->src_name,pPlateEntry->Plate);
        pListEntry1 = brightness_queue.flink; 
        while (pListEntry1 != NULL) {
          printf("Order %3d MAG_ISO %f\n",pListEntry1->brightnessIndex,pPlateEntry->magentry[pListEntry1->list_index].FLUX_ISO);
          pListEntry1 = pListEntry1->flink;
        }
      }
#endif

      /* Now step through the queue and eliminate any entries with an out-of-order mag_iso */
      pListEntry1 = NULL;
      pListEntry2 = brightness_queue.flink;      
      while (pListEntry2 != NULL) {
        pListEntry3 = pListEntry2->flink;
        if (pListEntry3 == NULL) {
          break;
        }
        if (pPlateEntry->magentry[pListEntry3->list_index].FLUX_ISO >
            pPlateEntry->magentry[pListEntry2->list_index].FLUX_ISO ) {
          /* Dequeue the dimmest one */
          pListEntry2->queued = -2;
          if (pListEntry1 != NULL) {
            pListEntry1->flink = pListEntry2->flink;
          } else {
            brightness_queue.flink = pListEntry2->flink;
          }
          brightness_index--;
          /* Start over */
          pListEntry1 = NULL;
          pListEntry2 = brightness_queue.flink;      
        } else {
          pListEntry1 = pListEntry2;
          pListEntry2 = pListEntry3;
        }
      }

#if 0
      if (strcmp(pListEntry->src_name,"K11812880") == 0) {
        printf("At %s for plate %s\n",pListEntry->src_name,pPlateEntry->Plate);
        pListEntry1 = brightness_queue.flink; 
        while (pListEntry1 != NULL) {
          printf("Order %3d MAG_ISO %f\n",pListEntry1->brightnessIndex,pPlateEntry->magentry[pListEntry1->list_index].FLUX_ISO);
          pListEntry1 = pListEntry1->flink;
        }
      }
#endif
      /* Attempt a second insert   We are particularly interested in tacking entries onto the end */
      workDone = 1;
      while (workDone) {
        workDone = 0;
        for (list_index2 = 0; list_index2 < list_nrecs; list_index2++) {
          pListEntry2 = &list_table[list_index2];
          if (pListEntry2->queued != -2) {
            continue;
          }
          pListEntry1 = brightness_queue.flink;
          while (pListEntry1 != NULL) {
            compare1good = 0;
            compare2good = 0;
            pListEntry3 = pListEntry1->flink;
            if (pListEntry3 == NULL) {
              break;
            }
            if ((pListEntry3 == pListEntry2) ||
                (pListEntry1 == pListEntry2)) {
              fprintf(stderr,"ERROR: consistency check!\n");
              exit(-1);
            }
            if ((pListEntry1->brightnessIndex < pListEntry2->brightnessIndex) && 
                (pPlateEntry->magentry[pListEntry2->list_index].FLUX_ISO < pPlateEntry->magentry[pListEntry1->list_index].FLUX_ISO )) {
              compare1good = 1; /* The first two stars compare fine */
            }
            if ((pListEntry3 != NULL) &&
                (pListEntry2->brightnessIndex < pListEntry3->brightnessIndex) && 
                (pPlateEntry->magentry[pListEntry2->list_index].FLUX_ISO > pPlateEntry->magentry[pListEntry3->list_index].FLUX_ISO )) {
              compare2good = 1; /* The next two stars compare fine */
            }   
            if ((compare1good == 1) && (compare2good == 1) ||
                (compare1good == 1) && (pListEntry3 == NULL)) {
              /* We can insert this one back into the queue */
              pListEntry2->queued = 1;
              pListEntry1->flink = pListEntry2;
              pListEntry2->flink = pListEntry3;
              brightness_index++;
              workDone = 1;
              break;
            }
            pListEntry1 = pListEntry1->flink;
          }
          if (workDone) {
            break;
          }
        }
     
      }


#if 0
      if (strcmp(pListEntry->src_name,"K11812880") == 0) {
        printf("At %s for plate %s\n",pListEntry->src_name,pPlateEntry->Plate);
        pListEntry1 = brightness_queue.flink; 
        while (pListEntry1 != NULL) {
          printf("Order %3d MAG_ISO %f\n",pListEntry1->brightnessIndex,pPlateEntry->magentry[pListEntry1->list_index].FLUX_ISO);
          pListEntry1 = pListEntry1->flink;
        }
      }
#endif



#if 0
      if (strcmp(pPlateEntry->Plate,"dnb02276") == 0) {
        printf("At plate %s\n",pPlateEntry->Plate);
      }
#endif
      /* Now step through the queue and get our magnitude */
      pListEntry1 = brightness_queue.flink; 
      if (pListEntry1 == NULL) {
        errorPoints++;
        /* Nothing on the queue for this star */
        continue;
      }
      if (pPlateEntry->magentry[list_index].FLUX_ISO > pPlateEntry->magentry[pListEntry1->list_index].FLUX_ISO) {
        /* Our star is brighter than the brighest star in the sequence */
        pPlateEntry->magEstimate[list_index] = pListEntry1->Stdmag;
        pPlateEntry->flag[list_index] = SYMBOL_TOO_BRIGHT;
        validPoints++;

      } else {
        pListEntry2 = pListEntry1->flink;
        while (pListEntry2 != NULL) {
          if (pPlateEntry->magentry[list_index].FLUX_ISO > pPlateEntry->magentry[pListEntry2->list_index].FLUX_ISO) {
            /* We can interpolate */
            pPlateEntry->magEstimate[list_index] = ComputeMag(pPlateEntry->magentry[list_index].FLUX_ISO,
                                                              pPlateEntry->magentry[pListEntry2->list_index].FLUX_ISO,
                                                              pPlateEntry->magentry[pListEntry1->list_index].FLUX_ISO,
                                                              pListEntry2->Stdmag,
                                                              pListEntry1->Stdmag);
            if ((pMagEntry1->FLUX_MAX/pMagEntry1->THRESHOLD) >= analysis_threshold) {

              pPlateEntry->flag[list_index] = SYMBOL_GOOD;
              vector[goodpoints] = pPlateEntry->magEstimate[list_index];
              goodpoints++;
            } else {
              pPlateEntry->flag[list_index] = SYMBOL_LIMITING1;
            }
            validPoints++;
            break;
          }
          pListEntry1 = pListEntry2;
          pListEntry2 = pListEntry1->flink;
        }
        if (pListEntry2 == NULL) {
          pPlateEntry->magEstimate[list_index] = pListEntry1->Stdmag;
          pPlateEntry->flag[list_index] = SYMBOL_NOT_FOUND;
          validPoints++;
        }
      }
      /* If we do not yet have a limiting magnitude from MAG_ISO_med, step through the list again and compute it */
      if (pPlateEntry->limiting_mag_local1 == 0.0) {
        pListEntry1 = brightness_queue.flink; 
        if (pPlateEntry->FLUX_ISO_med > pPlateEntry->magentry[pListEntry1->list_index].FLUX_ISO) {
          /* Our star is brighter than the brighest star in the sequence */
          pPlateEntry->limiting_mag_local1 = pListEntry1->Stdmag;
        } else {
          pListEntry2 = pListEntry1->flink;
          while (pListEntry2 != NULL) {
            if (pPlateEntry->FLUX_ISO_med > pPlateEntry->magentry[pListEntry2->list_index].FLUX_ISO) {
              /* We can interpolate */
              pPlateEntry->limiting_mag_local1 = ComputeMag(pPlateEntry->FLUX_ISO_med,
                                                            pPlateEntry->magentry[pListEntry2->list_index].FLUX_ISO,
                                                            pPlateEntry->magentry[pListEntry1->list_index].FLUX_ISO,
                                                            pListEntry2->Stdmag,
                                                            pListEntry1->Stdmag);
              break;
            }
            pListEntry1 = pListEntry2;
            pListEntry2 = pListEntry1->flink;
          }
          if (pListEntry2 == NULL) {
            pPlateEntry->limiting_mag_local1 = pListEntry1->Stdmag;
          }
        }
      }
      /* Calculate a second limiting magnitude from the last entry */
      if (pPlateEntry->limiting_mag_local2 == 0.0) {
        pListEntry1 = brightness_queue.flink; 
        pListEntry2 = brightness_queue.flink; 
        while (pListEntry2 != NULL) {
          pListEntry1 = pListEntry2;
          pListEntry2 = pListEntry1->flink;
        }
        if (pListEntry1 != NULL) {
          pPlateEntry->limiting_mag_local2 = pListEntry1->Stdmag;
          
        } else {
          pPlateEntry->limiting_mag_local2 = 99.0;
        }
      }
    
  


    } /* End of plate loop */
#if 0
    /* Print the relative number of times we found this object, but not for the object in question */
    if (list_index != 0) {
      fprintf(outputHandle,"%f %f\n",pListEntry->Stdmag,(1.0*foundPoints)/(1.0*plate_nrecs));
    }
#endif
    if (goodpoints > 0) {
      if (CalcMedianAndRMS(goodpoints,2,vector,&rawmed,&rawrms,0,3.0,0) == 0) {
        rawmed = 0.0;
        rawrms = 99.0;
      }
      if (CalcMedianAndRMS(goodpoints,2,vector,&clipmed,&cliprms,1,3.0,0) == 0) {
        clipmed = 0.0;
        cliprms = 99.0;
      }
    


      strcpy(lightcurveName,workingDirectory);
      strcat(lightcurveName,"/lc_");
      strcat(lightcurveName,pListEntry->src_name);
      strcat(lightcurveName,"_all");
      strcat(lightcurveName,".db");
      lightcurveHandle = fopen(lightcurveName,"wt");
      if (lightcurveHandle == NULL) {
        printf("ERROR: Failed to open %s\n",lightcurveName);
        return(-1);
      }
      sprintf(lightcurveplot,"%s/lc_%s.ps/ps",workingDirectory,pListEntry->src_name);
      cpgopen(lightcurveplot);
      /* cpgeras(); */
      cpgslw(4); /* line width is 4/200 inch */
      cpgsch(1.5); /* Character size is 1.5/40 the height of the view surface */
      /* subdivide the page into 2 panels;*/
      cpgsubp(2,2);


#ifdef TOPCAT_OUTPUT
      fprintf(lightcurveHandle,"sequence\thJD\tflag\tmagcal_local\tPlate\n");
      fprintf(lightcurveHandle,"--------\t---\t----\t------------\t-----\n");

#else /* TOPCAT_OUTPUT */
      fprintf(lightcurveHandle,"sequence\thJD\tspatial_bin\tflag\tmagcal_local\tmagcal_local_rms\tnlocal\tPlate\tmagcal_iso_rms\tmagcal_local_error\tmagcal_iso\tlimiting_mag_local\tAFLAGS\treject_reason1\treject_reason2\tdradRMS2");
      fprintf(lightcurveHandle,"\n");

      fprintf(lightcurveHandle,"--------\t---\t-----------\t----\t------------\t----------------\t------\t-----\t--------------\t-----------------\t----------\t------------------\t-----\t--------------\t--------------\t--------");
      fprintf(lightcurveHandle,"\n");
#endif /* TOPCAT_OUTPUT */
      curve_nrecs = 0;
      for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
        pPlateEntry = &plate_table[plate_index];
        if (pPlateEntry->Date == 0) {
          continue;
        }
        if (pPlateEntry->magEstimate[list_index] != 0) {
          endJD = pPlateEntry->Date;
          if (startJD == 0.0) {
            startJD = endJD;
          }
          sequence++;
          sequenceval[curve_nrecs] = sequence;
          magval[curve_nrecs] = pPlateEntry->magEstimate[list_index];
          yearval = jd2ep(pPlateEntry->magentry[list_index].Date);
          timeval[curve_nrecs] = yearval;
          if (pPlateEntry->flag[list_index] == SYMBOL_TOO_BRIGHT) {
            symbolval[curve_nrecs] = 30;
          } else if (pPlateEntry->flag[list_index] == SYMBOL_GOOD) {
            symbolval[curve_nrecs] = 17;
          } else if (pPlateEntry->flag[list_index] == SYMBOL_LIMITING1) {
            symbolval[curve_nrecs] = 4;
          } else if (pPlateEntry->flag[list_index] == SYMBOL_NOT_FOUND) {
            symbolval[curve_nrecs] = 31;
          } else {
            fprintf(stderr,"ERROR: unrecognized symbol\n");
            exit(-1);
          }
#if 0
          printf("%d %f %f %d\n",curve_nrecs,sequenceval[curve_nrecs],magval[curve_nrecs],symbolval[curve_nrecs]);
#endif
          curve_nrecs++;
#ifdef TOPCAT_OUTPUT
          fprintf(lightcurveHandle,"%d\t%f\t%.1f\t%.2f\t%s\n",
                  sequence,
                  pPlateEntry->magentry[list_index].Date,
                  pPlateEntry->flag[list_index],
                  pPlateEntry->magEstimate[list_index],
                  pPlateEntry->Plate);
          if ((pPlateEntry->limiting_mag_local1 != 0.0) &&
              (pPlateEntry->limiting_mag_local1 != 99.0)) {
            fprintf(lightcurveHandle,"%d\t%f\t%.1f\t%.2f\t%s\n",
                    sequence,
                    pPlateEntry->magentry[list_index].Date,
                    SYMBOL_LIMITING1,
                    pPlateEntry->limiting_mag_local1,
                    pPlateEntry->Plate);

          }
          if ((pPlateEntry->limiting_mag_local2 != 0.0) &&
              (pPlateEntry->limiting_mag_local2 != 99.0)) {
            fprintf(lightcurveHandle,"%d\t%f\t%.1f\t%.2f\t%s\n",
                    sequence,
                    pPlateEntry->magentry[list_index].Date,
                    SYMBOL_LIMITING2,
                    pPlateEntry->limiting_mag_local2,
                    pPlateEntry->Plate);

          }
          
#else /* TOPCAT_OUTPUT */
          fprintf(lightcurveHandle,"%d\t%f\t%d\t%.1f\t%.2f\t%.2f\t%d\t%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\n",
                  sequence,
                  pPlateEntry->magentry[list_index].Date,
                  pPlateEntry->magentry[list_index].spatial_bin,
                  pPlateEntry->flag[list_index],
                  pPlateEntry->magEstimate[list_index],
                  pPlateEntry->magcal_local_rms*error_bar_factor,
                  pPlateEntry->npoints_local,
                  pPlateEntry->Plate,
                  pPlateEntry->magcal_iso_rms*error_bar_factor,
                  pPlateEntry->magcal_local_error*error_bar_factor,
                  pPlateEntry->magcal_iso,
                  pPlateEntry->limiting_mag_local1,
                  pPlateEntry->magentry[list_index].AFLAGS,
                  pPlateEntry->reject_reason1,
                  pPlateEntry->reject_reason2,
                  pPlateEntry->dradRMS2);
#endif /* TOPCAT_OUTPUT */
        }

      }
      if (curve_nrecs > 1) {
        timemin = timeval[0];
        timemax = timeval[0];
        magmin = magval[0];
        magmax = magval[0];
        for (curve_index = 1; curve_index < curve_nrecs; curve_index++) {
#if 0
          printf("%d %f %f %d\n",curve_index,timeval[curve_index],magval[curve_index],symbolval[curve_index]);
#endif
          if (timemin > timeval[curve_index]) {
            timemin = timeval[curve_index];
          }
          if (timemax < timeval[curve_index]) {
            timemax = timeval[curve_index];
          }
          if (magmin > magval[curve_index]) {
            magmin = magval[curve_index];
          }
          if (magmax < magval[curve_index]) {
            magmax = magval[curve_index];
          }
        }
        /* cpgpanl(1,1); */
        cpgenv(timemin,timemax,magmax,magmin,0,0);
        sprintf(label,"median_local %.2f rms_local %.2f clip_median_local %.2f clip_rms_local %.2f npoints %d ngood %d REF %s",
                rawmed,rawrms,clipmed,cliprms,validPoints,goodpoints,pListEntry->src_name);
        cpgmtxt("T",2.0,0.0,0.0,label);
        cpglab("Year","magcal_local","");
        cpgpnts(curve_nrecs,timeval,magval,symbolval,curve_nrecs);

        sequencemin = sequenceval[0];
        sequencemax = sequenceval[0];
        magmin = magval[0];
        magmax = magval[0];
        for (curve_index = 1; curve_index < curve_nrecs; curve_index++) {
#if 0
          printf("%d %f %f %d\n",curve_index,sequenceval[curve_index],magval[curve_index],symbolval[curve_index]);
#endif
          if (sequencemin > sequenceval[curve_index]) {
            sequencemin = sequenceval[curve_index];
          }
          if (sequencemax < sequenceval[curve_index]) {
            sequencemax = sequenceval[curve_index];
          }
          if (magmin > magval[curve_index]) {
            magmin = magval[curve_index];
          }
          if (magmax < magval[curve_index]) {
            magmax = magval[curve_index];
          }
        }
        cpgpanl(1,1);
        cpgenv(sequencemin,sequencemax,magmax,magmin,0,0);
        cpglab("Sequence","magcal_local","");
        cpgpnts(curve_nrecs,sequenceval,magval,symbolval,curve_nrecs);


      }

      if ((base_min_threshold <= 0) &&
          (threshold_nrecs > 1)) {
        cpgpanl(2,1);

        cpgenv(0,xthreshold_max,0,1.0*best_selection_count,0,0);
        /* Plot our iteration results */
        cpglab("min_threshold","Number of sequenc points","Search for the best selection threshold");
        cpgpt(threshold_nrecs,xmin_threshold,ysequencesize,4);
      }
      
      cpgpanl(1,2);
      cpgenv(xmin,xmax,ymax,ymin,0,0);
      /* Plot our comparison sequence */
      cpglab("Sequence","Stdmag","Comparison Star Sequence");
      cpgpt(plot_nrecs,xval,yval,4);



      fclose(lightcurveHandle);
      cpgend();
      /* Now write out the summary text */
      startYear = jd2ep(startJD);
      endYear = jd2ep(endJD);
      strcpy(lightcurveName,workingDirectory);
      strcat(lightcurveName,"/lc_");
      strcat(lightcurveName,pListEntry->src_name);
      strcat(lightcurveName,"_all");
      strcat(lightcurveName,".txt");
      lightcurveHandle = fopen(lightcurveName,"wt");
      fprintf(lightcurveHandle,"#raw median mag, raw median rms, clipped median mag, clipped median rms, total points, good points, start year, end year,error_bar_factor\n");
      fprintf(lightcurveHandle,"%.2f %.2f %.2f %.2f %d %d %f %f %f\n",
              rawmed,rawrms,clipmed,cliprms,validPoints,goodpoints,startYear,endYear,error_bar_factor);
      pListEntry->rawmed = rawmed;
      pListEntry->rawrms = rawrms;
      fprintf(outputHandle,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%f\t%f\t%s\n",rawmed,pListEntry->Stdmag,rawrms,clipmed,cliprms,validPoints,goodpoints,startYear,endYear,pListEntry->src_name);

      if (lightcurveHandle == NULL) {
        printf("ERROR: Failed to open %s\n",lightcurveName);
        return(-1);
      }
      fclose(lightcurveHandle);
      printf("Lightcurve written good: %4d err: %4d bad: %4d bin9 %4d for %s %s\n",goodpoints,errorPoints,badPoints,bin9Points,pListEntry->src_name,pListEntry->REF);
      outputPlotCount++;




    } else {
      printf("ERROR: No good points available for entry index %d %s\n",pListEntry->list_index,pListEntry->src_name);
    }

  } /* End of candidate loop */

    /* Have all the lightcurves, now prepare our summary statistics */
  minmag = BRIGHTEST_MAGNITUDE-1;
  maxmag = BRIGHTEST_MAGNITUDE;
  printf(" mag    count  median\n");
  printf("                 rms\n");
  while ((minmag != BRIGHTEST_MAGNITUDE) || (maxmag != DIMMEST_MAGNITUDE)) {
    minmag++;
    maxmag++;
    if (minmag >= DIMMEST_MAGNITUDE) {
      minmag = BRIGHTEST_MAGNITUDE;
      maxmag = DIMMEST_MAGNITUDE;
    }

    goodpoints = 0;
      
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      pListEntry = &list_table[list_index];
      if ((pListEntry->rawrms < 90.0) &&
          (pListEntry->rawrms > 0.0) &&
          (pListEntry->rawmed >= (1.0*minmag)) &&
          (pListEntry->rawmed <  (1.0*maxmag))) {
        vector1[goodpoints] = pListEntry->rawrms;
        goodpoints++;
      }
    }
    if (goodpoints > 1) {
      if (CalcMedianAndRMS(goodpoints,0,vector1,&rawmed,&rawrms,0,3.0,0) > 0) {
        fprintf(summaryHandle,"%d\t%d\t\%5d\t%6.3f\n",
                minmag,
                maxmag,
                goodpoints,
                rawmed);
        printf("%2d %2d %7d %8.6f\n",minmag,maxmag,goodpoints,rawmed);
      }
    }
  }
#if 0
  for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
    pPlateEntry = &plate_table[plate_index];
    if (pPlateEntry->Date == 0) {
      continue;
    }
    printf("Plate %8s |",pPlateEntry->Plate);
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      printf(" %7.0f %d |",pPlateEntry->magentrySort[list_index].FLUX_ISO,pPlateEntry->magentrySort[list_index].list_index);
    }
    printf("\n");
  }
#endif


#if 0
  for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
    pPlateEntry = &plate_table[plate_index];
    printf("Plate %8s |",pPlateEntry->Plate);
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      printf(" %7.2f %4.1f |",pPlateEntry->magEstimate[list_index],pPlateEntry->flag[list_index]);
    }
    printf("\n");
  }
#endif
  
#if 1
  strcpy(plotdev,workingDirectory);
  strcat(plotdev,"/limiting_plot.ps/ps");
  cpgopen(plotdev);
  cpgslw(4); /* line width is 4/200 inch */
  cpgsch(1.5); /* Character size is 1.5/40 the height of the view surface */
  cpgenv(ymax,ymin,ymax,ymin,0,0);
  cpglab("Stdmag","magcal_local","");
  /*  cpgsubp(2,2); */

  /* Here we print out plate-specific files which have magcal_local as a function of Stdmag */
  for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
    char outfile[MAX_BUFFER];
    FILE *outHandle;
    pPlateEntry = &plate_table[plate_index];
    plot_nrecs = 0;
#if 0
    sprintf(outfile,"%s/%s_calibration.txt",workingDirectory,pPlateEntry->Plate);
    outHandle = fopen(outfile,"wt");
    if (outHandle != NULL) {
      for (list_index = 0; list_index < list_nrecs; list_index++) {
        pListEntry = &list_table[list_index];
        if ((pPlateEntry->magEstimate[list_index] > 0.0) &&
            (pPlateEntry->magEstimate[list_index] < 90.0) &&
            (pPlateEntry->flag[list_index] == SYMBOL_GOOD)) {
          fprintf(outHandle,"%.2f %.2f\n",pListEntry->Stdmag,pPlateEntry->magEstimate[list_index]);
        }
      }

      fclose(outHandle);
    } else {
      printf("ERROR: failed to open %s\n",outfile);
      exit(-1);
    }
#endif
    for (list_index = 0; list_index < list_nrecs; list_index++) {
      pListEntry = &list_table[list_index];
      if ((pPlateEntry->magEstimate[list_index] > 0.0) &&
          (pPlateEntry->magEstimate[list_index] < 90.0) &&
          (pPlateEntry->flag[list_index] == SYMBOL_GOOD)) {
        if (pListEntry->Stdmag < ymax) {
          xval[plot_nrecs] = pListEntry->Stdmag;
          yval[plot_nrecs] = pPlateEntry->magEstimate[list_index];
          plot_nrecs++;
        }
      }
    }
    if (plot_nrecs >= 2) {
      /* cpgeras(); */
      cpgpt(plot_nrecs,xval,yval,4);
      cpgline(plot_nrecs,xval,yval);
      /* cpgpage();  */


    }


  }
  cpgend();
#endif





  /* All done.  Clean up */
  
  if (master_header != NULL) {
    table_hdrfree(master_header);
  }
  if (master_handle != NULL) {
    Close(master_handle);
  }
  if (master_descriptor != NULL) {
    Free(master_descriptor);
  }
  if (master_row != NULL) {
    table_rowfree(master_row);
  }
  if (vector != NULL) {
    free(vector);
  }
  if (vector1 != NULL) {
    free(vector1);
  }

  if (list_table != NULL) {
    free(list_table);
  }
  if (list_header != NULL) {
    table_hdrfree(list_header);
  }
  if (list_handle != NULL) {
    Close(list_handle);
  }


  if (magiso_table != NULL) {
    free(magiso_table);
  }
  if (magiso_header != NULL) {
    table_hdrfree(magiso_header);
  }
  if (magiso_handle != NULL) {
    Close(magiso_handle);
  }


  if (master_table != NULL) {
    free(master_table);
  }

  if (plate_table != NULL) {
    free(plate_table);
  }
  if (outputHandle != NULL) {
    fclose(outputHandle);
  }
  if (summaryHandle != NULL) {
    fclose(summaryHandle);
  }
  if (xval != NULL) {
    free(xval);
  }
  if (yval != NULL) {
    free(yval);
  }
  if (timeval != NULL) {
    free(timeval);
  }
  if (magval != NULL) {
    free(magval);
  }
  if (symbolval != NULL) {
    free(symbolval);
  }

  time(&curTime);
  curTime -= startTime;
  fprintf(stdout,"Output stars %d seconds %d\n",
          outputPlotCount,curTime);
    


  return(0);
}

