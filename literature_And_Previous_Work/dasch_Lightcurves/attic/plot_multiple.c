// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* plot_multiple.c
 *
 * Obtain statistics on the number of multiple exposure plates in the stacks as a function of year and plate scale 
 *
 * cc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -I /dasch/install/pgplot  -L /dasch/install/lib -lm -lcfitsio   plot_multiple.c pipelineutils.a -ltable -lutil  -lwcs -o plot_multiple -L/usr/lib${lib64}/mysql -lmysqlclient  -L /dasch/install/pgplot  -lcpgplot -lpgplot  -L/usr/X11R6/lib  -lX11 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 /usr/lib${lib64}/libg2c.so.0  -ldl -pthread 
 *
 * Sep  7, 2010 Edward J. Los - Original Version
 * Sep 24, 2010 Edward J. Los - Add grating plates
 * Feb 27, 2012 Edward J. Los - Add optional cal to CheckQuality
 * Jun 19, 2012 Edward J. Los - Add a pgplot file showing connecting multiple exposures for the
 *                              purposes of finding calibration fields.
 * Jun 29, 2012 Edward J. Los - Use the logbook dRightAscension and dDeclination vs any that were entered by imWCS.
 *
 *
 For plotting only, no mosaics:
 ./plot_multiple -o los.tmp -p

 ./plot_multiple -o los.tmp -m
 Execution Time: 2 seconds 217058 plates  5842 multiple exposure plates year: 1885.723364-1992.667571
 maxExposureNumber 28 
 minmax los.tmp
 psxy los.tmp -R1880/2000/0/30 -Y4i -JX6i/5.5i -P -B10:"Year":/5:"Percentage of Multiple Exposure Plates":WSne -W0.5 -K > ! los.ps
 gv los.ps &
 minmax lostotal.tmp
 psxy lostotal.tmp -R1880/2000/0/7000 -Y4i -JX6i/5.5i -P -B10:"Year":/1000:"Number of Plates":WSne -W0.5 -K > ! lostotal.ps
 gv lostotal.ps &
 minmax losscale.tmp
 psxy losscale.tmp -R0/8/0/20 -Y4i -JX6i/5.5i -P -B1:"arcsec/pixel":/5:"Percentage of Multiple Exposure Plates":WSne -Sc0.1 -K > ! losscale.ps
 gv losscale.ps


 minmax losseries.tmp
 psxy losseries.tmp -R0/8/0/50000  -JX6i/5.5i -Y4i -P -B1:"arcsec/pixel":/10000:"Number of Plates":WSne -Sc0.1 -K > ! losseries.ps
 gv losseries.ps

 minmax losgrating.tmp
 psxy losgrating.tmp -R1880/2000/0/0.4 -Y4i -JX6i/5.5i -P -B10:"Year":/0.1:"Percentage of Grating Plates":WSne -W0.5 -K > ! losgrating.ps
 gv losgrating.ps &


*/


#include <math.h>
#include "table.h"
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include "cpgplot.h"
#define MAX_FILENAME 256
#define MAX_BUFFER 256
#define MAX_LIST_STRING 25
#define SAMPLE_POINTS 49
#define ALLOC_INCREMENT 1000
#define MINIMUM_COUNT 20
#define MAX_QUERY_COUNT 200
#define DRAD_ERROR_FACTOR 3.0
#define MAX_YEARS (PIPELINE_MAX_DATE-PIPELINE_MIN_DATE)
#define YEAR_INCREMENT 2
#define MAX_SOLUTIONS 10
#if 1
#define PLOT_XMIN -10.0
#define PLOT_XMAX 370.0
#define PLOT_YMIN -95.0
#define PLOT_YMAX  95.0
#endif

#if 0
/* North pole */
#define PLOT_XMIN -10.0
#define PLOT_XMAX 370.0
#define PLOT_YMIN  85.0
#define PLOT_YMAX  91.0
#endif

#if 0
/* South pole */
#define PLOT_XMIN -10.0
#define PLOT_XMAX 370.0
#define PLOT_YMIN  -91.0
#define PLOT_YMAX  -85.0
#endif

#if 0
/* +20 dec */
#define PLOT_XMIN -10.0
#define PLOT_XMAX 370.0
#define PLOT_YMIN  15.0
#define PLOT_YMAX  25.0
#endif

#define MINIMUM_DISTANCE 5.0 /* The exposures in a double exposure must be at least this far apart */
#define HOUR_CONVERT 15.0
#define POLE_BIN 85.0

typedef struct _mosaicEntry {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int scanNumber;
  int mosaicNumber;
} MOSAICENTRY,*PMOSAICENTRY;

typedef struct _multipleEntry {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int exposureCount;
  EXPOSURE exposure[MAX_CATALOG_EXPOSURES];
} MULTIPLEENTRY,*PMULTIPLEENTRY;

/* Regions are centered on integral RA and DEC values +/- 0.5 degrees */
#define DEC_BINS 181
#define RA_BINS  361
#define REGION_BINS ((DEC_BINS) * (RA_BINS))
/* decIndex = 0-180 representing -90 to +90 degrees
   raIndex  = 0-360
   regionIndex = raIndex + ((RA_BINS)*decIndex);
*/

typedef struct _regionTable {
  double minRA;
  double maxRA;
  double minDec;
  double maxDec;
  double centerRA;
  double centerDec;
  int endpointcount;
  int selected;
} REGION,*PREGION;
float xval[REGION_BINS];
float yval[REGION_BINS];



#if 1
float TR[6];
REGION regionTable[REGION_BINS];  /* The original table */
#if 0
float PLOTOUT[DEC_BINS][RA_BINS];
#endif
#else


REGION regionTable[] = {
  { -1.0,361.0, 89.0, 91.0}, /* North Pole */
  { -1.0,361.0,-91.0,-89.0}, /* South Pole */
  { -1.0,361.0, 19.0, 21.0}, /* +20 dec */
  { -1.0,  1.0, -1.0,  1.0}, /* 0h 0d */
  {360.0,361.0, -1.0,  1.0}, /* 0h 0d */
  { -1.0,361.0, 14.0, 16.0}, /* +15 dec */
  { -1.0,361.0,-16.0,-14.0}, /* -15 dec */
  { -1.0,361.0,-46.0,-44.0}, /* -45 dec */


};

#endif
int regionTableSize = sizeof(regionTable)/sizeof(REGION);
int RegionCompare(const void *first, const void *second) 
{
  int result;
  PREGION pFirst  = (PREGION)first;
  PREGION pSecond = (PREGION)second;
  if ((pFirst->selected == 1) && (pSecond->selected == 0)) {
      return(1);
  } else if ((pFirst->selected == 0) && (pSecond->selected == 1)) {
      return(-1);
  } else if ((pFirst->selected == 0) && (pSecond->selected == 0)) {
    return(0);
  } else {
    /* This is a REVERSE sort */
    if (pFirst->endpointcount > pSecond->endpointcount) {
      return(-1);
    } else if (pFirst->endpointcount < pSecond->endpointcount) {
      return(1);
    } else {
      return(0);
    }
  }
}

int main(int argc,char *argv[])
{
  char *dotPtr;
  char outfile[MAX_BUFFER];
  char prefix[MAX_BUFFER];
  char suffix[MAX_BUFFER];
  FILE *outHandle = NULL;
  char outscalefile[MAX_BUFFER];
  FILE *outscalehandle = NULL;
  char outregionfile[MAX_BUFFER];
  FILE *outregionhandle = NULL;
  char outtotalfile[MAX_BUFFER];
  FILE *outtotalhandle = NULL;
  char outseriesfile[MAX_BUFFER];
  FILE *outserieshandle = NULL;
  char outgratingfile[MAX_BUFFER];
  FILE *outgratinghandle = NULL;
  char outvectorfile[MAX_BUFFER];

  int errorFlag = 0;
  int verbose = 0;
  int doPlots = 0;
  int checkMosaics = 0;
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
  int series_nrecs = 0;
  int series_index;


  int spatial_bin_nrecs = 0;
  int spatial_bin_index;

  int plate_alloc = 0;
  int plateCount = 0;
  int multipleCount = 0;
  int plate_index;
  double minYear = 9999;
  double maxYear = 0;
  int maxExposureNumber = 0;
  int baseYear;
  int topYear;
 
  int yearIndex;
  int yearCount0[MAX_YEARS]; /* All plates per year */
  int yearCount1[MAX_YEARS]; /* Multiple exposure paltes per year */

  int seriesIndex;
  int seriesCount0[MAX_SERIES];
  int seriesCount1[MAX_SERIES];
  int seriesId;

  int gratingCount0[MAX_YEARS];
  int gratingCount1[MAX_YEARS];
  int classCount = 0;
  int gratingCount = 0;

  double yearRatio;
  double gratingRatio;
  double scaleRatio;

  int maxYearCount = 0;
  double *vector;
  int spatial_bin;
  double binmedian;
  double binRMS;
  int nvals;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  MOSAIC mosaicInfo;
  PHOTPLATES basePhotPlates;
  PPHOTPLATES pPhotPlates = &basePhotPlates;
  int gotAnswer;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  int solutionNumberFix = 0;

  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char queryString[MAX_QUERY_STRING];
  char series[MAX_SERIES_STRING];
  int exposureNumber;
  char dateString[MAX_DATE_STRING];
  double epoch;
  int year;
  char class[MAX_CLASS_STRING];

  PMOSAICENTRY mosaicListTable = NULL;
  PMOSAICENTRY pMosaicEntryCur;
  PMOSAICENTRY pMosaicEntryNext;
  int numMosaicEntries;
  int curMosaicCount = 0;
  int duplicateMosaicCount = 0;
  int nonDuplicateMosaicCount = 0;
  int mosaicIndex;

  MOSAIC mosaicTable[MAX_SOLUTIONS];
  PMOSAIC pMosaic;
  PMOSAIC pMosaic2;
  int wedgeFlag;
  int solutionCount;
  EXPOSURE exposureTable[MAX_CATALOG_EXPOSURES];
  PEXPOSURE pExposure;
  int exposureCount;
  int totalSolutions = 0;
  int totalExposures = 0;
  int missedExposureCount = 0;
  int goodExposureCount = 0;
  int missedExposureTable[MAX_CATALOG_EXPOSURES];
  int totalExposureTable[MAX_CATALOG_EXPOSURES];
  int goodExposureTable[MAX_CATALOG_EXPOSURES];
  int exposureIndex;
  int plateNumber;

  int unmatchedSolutionCount = 0;
  int unmatchedSolutionTable[MAX_SOLUTIONS];
  int goodUnmatchedCount = 0;
  int goodUnmatchedTable[MAX_SOLUTIONS];
  int wedgeSolutionCount = 0;
  int wedgeSolutionTable[MAX_SOLUTIONS];
  int totalSolutionTable[MAX_SOLUTIONS];
  int solutionIndex;
  int separationTable[MAX_SOLUTIONS];
  int staleTable[MAX_SOLUTIONS];
  int noAllobjects[MAX_SOLUTIONS];
  int duplicateRef[MAX_SOLUTIONS];
  int checkQuality = 0;
  double oldRightAscension;
  double oldDeclination;
  double factor;
  double distance;

  PMULTIPLEENTRY pMultipleTable = NULL;
  PMULTIPLEENTRY pMultipleEntry;
  int maxMultipleExposures = 0;
  int curMultipleExposures = 0;
  int multIndex;
  int pointCount;
  int plotCount;
  PREGION pRegion;
  int regionIndex;
  int plotFlag;
  int pointIndex;
  char plotString[MAX_BUFFER];
  int rejectFlag;

  int decIndex;
  int raIndex;
  int iteration;
  int workDone;
  int maxEndpointCount;
  int topEndpointCount = 0;
  char locationText[20];

  memset(regionTable,0,sizeof(regionTable));
  for (decIndex = 0; decIndex < DEC_BINS; decIndex++) {
    for (raIndex = 0; raIndex < RA_BINS; raIndex++) {
      regionIndex = raIndex + (RA_BINS * decIndex);
      pRegion = &regionTable[regionIndex];
      pRegion->centerRA = 1.0 * raIndex;
      pRegion->centerDec = (1.0* decIndex) -90.0;
      pRegion->minRA = pRegion->centerRA-0.5;
      pRegion->maxRA = pRegion->centerRA+0.5;
      pRegion->minDec = pRegion->centerDec-0.5;
      pRegion->maxDec = pRegion->centerDec+0.5;
#if 0
      printf("decIndex %5d raIndex %5d regionIndex %6d centerRA %8.2f centerDec %8.2f\n",
             decIndex,raIndex,regionIndex,pRegion->centerRA,pRegion->centerDec);
#endif
    }
  }


  memset(yearCount0,0,sizeof(yearCount0));
  memset(yearCount1,0,sizeof(yearCount1));
  memset(seriesCount0,0,sizeof(seriesCount0));
  memset(seriesCount1,0,sizeof(seriesCount1));
  memset(gratingCount0,0,sizeof(gratingCount0));
  memset(gratingCount1,0,sizeof(gratingCount1));
  memset(missedExposureTable,0,sizeof(missedExposureTable));
  memset(totalExposureTable,0,sizeof(totalExposureTable));
  memset(goodExposureTable,0,sizeof(goodExposureTable));
  memset(unmatchedSolutionTable,0,sizeof(unmatchedSolutionTable));
  memset(goodUnmatchedTable,0,sizeof(goodUnmatchedTable));
  memset(wedgeSolutionTable,0,sizeof(wedgeSolutionTable));
  memset(totalSolutionTable,0,sizeof(totalSolutionTable));

  memset(separationTable,0,sizeof(separationTable));
  memset(noAllobjects,0,sizeof(noAllobjects));
  memset(staleTable,0,sizeof(staleTable));
  memset(duplicateRef,0,sizeof(duplicateRef));


  time(&startTime);





  catalogString[0] = 0;
  outfile[0] = 0;

 
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


        case 'c': /* Update the quality field */
        case 'C':
          checkQuality = 1;
          break;




        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'p': /* verbose */
        case 'P':
          doPlots = 1;
          break;

        case 'm': /* check mosaics */
        case 'M':
          checkMosaics = 1;
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

    sprintf(outscalefile,"%sscale%s",prefix,suffix);
    sprintf(outregionfile,"%sregion%s",prefix,suffix);
    sprintf(outtotalfile,"%stotal%s",prefix,suffix);
    sprintf(outseriesfile,"%sseries%s",prefix,suffix);
    sprintf(outgratingfile,"%sgrating%s",prefix,suffix);



    


    outHandle = fopen(outfile,"wt");
    if (outHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output file %s\n",outfile);
      }
    }


    outscalehandle = fopen(outscalefile,"wt");
    if (outscalehandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the Scale output file %s\n",outscalefile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Scale file %s\n",outscalefile);
      }
    }

    outregionhandle = fopen(outregionfile,"wt");
    if (outregionhandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the Region output file %s\n",outregionfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Region file %s\n",outregionfile);
      }
    }

    outtotalhandle = fopen(outtotalfile,"wt");
    if (outtotalhandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the Total output file %s\n",outtotalfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Total file %s\n",outtotalfile);
      }
    }

    outserieshandle = fopen(outseriesfile,"wt");
    if (outserieshandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the Series output file %s\n",outseriesfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Series file %s\n",outseriesfile);
      }
    }

    outgratinghandle = fopen(outgratingfile,"wt");
    if (outgratinghandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the Grating output file %s\n",outgratingfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output Grating file %s\n",outgratingfile);
      }
    }




  } else {
    printf("ERROR: No output file specified\n");
    errorFlag = 1;
  }

  if (errorFlag) {
    printf("Usage: plot_multiple options\n");
    printf("  options: -v verbose\n");    
    printf("           -p do plots\n");    
    printf("           -m check mosaics\n");
    printf("           -c Invoke CheckQuality to update quality fields\n");
    printf("           -o <output file>\n");
    return(-1);
  }

  printf("plot_multiple of %s %s, out file %s \n",
         __DATE__,__TIME__,outfile);
    

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
  
  if (checkQuality) {
    printf("Invoking CheckQuality\n");
    CheckQuality(pConnection);
  }



  if (doPlots) {
    sprintf(queryString,"SELECT series,exposureNumber,date,class,plateNumber from exposures INNER JOIN plates using (series,plateNumber)",series);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        maxMultipleExposures = mysql_affected_rows(pConnection);
        pMultipleTable = (PMULTIPLEENTRY)calloc(maxMultipleExposures,sizeof(MULTIPLEENTRY)); 
        if (pMultipleTable == NULL) {
          printf("ERROR: failed to allocate multiple exposure table\n");
          exit(-1);
        }
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          pMultipleEntry = &pMultipleTable[curMultipleExposures];
          memset(pMultipleEntry,0,sizeof(MULTIPLEENTRY));
          if (sqlrow[0]) {
            strcpy(series,sqlrow[0]);
          } else {
            continue;
          }
          if (sqlrow[1]) {
            nvals = sscanf(sqlrow[1],"%d",&exposureNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: nvals is %d for exposureNumber\n",nvals);
              continue;
            } 
          } else {
            continue;
          }

          if (sqlrow[2]) {
            strcpy(dateString,sqlrow[2]);
          } else {
            continue;
          }

          if (sqlrow[3]) {
            strcpy(class,sqlrow[3]);
          } else {
            class[0] = 0;
          }

          if (sqlrow[4]) {
            nvals = sscanf(sqlrow[4],"%d",&plateNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: nvals is %d for plateNumber\n",nvals);
              continue;
            } 
          } else {
            continue;
          }


          epoch = fd2ep(dateString);
          year = epoch;
          if ((year < PIPELINE_MIN_DATE) || (year >= PIPELINE_MAX_DATE)) {
            printf("ERROR: invalid year %d for %s%05d\n",year,series,plateNumber);
            continue;
          }
          if (epoch < minYear) {
            minYear = epoch;
          }
          if (epoch > maxYear) {
            maxYear = epoch;
          }
          yearIndex = year-PIPELINE_MIN_DATE;
          if (exposureNumber > maxExposureNumber) {
            maxExposureNumber = exposureNumber;
          }
          seriesId = GetSeriesId(series,1);
          if (exposureNumber == 0) {
            plateCount++;
            yearCount0[yearIndex]++;
            seriesCount0[seriesId]++;
            if (strlen(class) > 0) {
              classCount++;
              gratingCount0[yearIndex]++;
              if (strstr(class,"G") != NULL) {
                gratingCount++;
                gratingCount1[yearIndex]++;

              }
            }

          }



          if (exposureNumber == 1) {
            if (curMultipleExposures >= maxMultipleExposures) {
              printf("ERROR: too many multiple exposures %d %d\n",curMultipleExposures,maxMultipleExposures);
              exit(-1);
            }
            strcpy(pMultipleEntry->series,series);
            pMultipleEntry->plateNumber = plateNumber;
            curMultipleExposures++;
            multipleCount++;
            yearCount1[yearIndex]++;
            seriesCount1[seriesId]++;
          }
	

        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: mysql_store_result failed in plot_multiple.c line %d\n",__LINE__);
      }
    }
    for (yearIndex = 0; yearIndex < MAX_YEARS; yearIndex++) {
      if (yearCount0[yearIndex] == 0) {
        yearRatio = 0.0;
      } else {
        yearRatio = (100.0 * yearCount1[yearIndex])/(1.0 * yearCount0[yearIndex]);
      }
      fprintf(outHandle,"%d %f\n",yearIndex+PIPELINE_MIN_DATE,yearRatio);
      fprintf(outtotalhandle,"%d %d\n",yearIndex+PIPELINE_MIN_DATE,yearCount0[yearIndex]);

      if (gratingCount0[yearIndex] == 0) {
        gratingRatio = 0.0;
      } else {
        gratingRatio = (100.0 * gratingCount1[yearIndex])/(1.0 * gratingCount0[yearIndex]);
      }
      fprintf(outgratinghandle,"%d %f\n",yearIndex+PIPELINE_MIN_DATE,gratingRatio);


    }
  
    for (seriesIndex = 0; seriesIndex < MAX_SERIES; seriesIndex++) {
      if (seriesCount0[seriesIndex] > 0) {
        double plateScale = 3600*GetFittedPlateScale(seriesIndex,99999);
        scaleRatio = (100.0 * seriesCount1[seriesIndex])/(1.0 * seriesCount0[seriesIndex]);
        if (seriesCount0[seriesIndex] >= 1000) {
          fprintf(outscalehandle,"%f %f\n",plateScale,scaleRatio);
        }
        fprintf(outserieshandle,"%f %d\n",plateScale,seriesCount0[seriesIndex]);
        printf("%5s %f %05d %3.0f\n",GetSeriesString(seriesIndex,1),plateScale,seriesCount0[seriesIndex],scaleRatio);
      }
    }
    maxMultipleExposures = curMultipleExposures;
  } /* doPlots */

  /* Now we need to estimate the number of solutions found compared with the number registered in the logbooks */
  if (checkMosaics) {
    sprintf(queryString,"SELECT series,plateNumber,scanNumber,mosaicNumber from mosaics where solutionNumber = 0 and WCSSource = 'imWCS' and ((mosaicComment IS NULL) OR (mosaicComment NOT regexp('Deleted'))) order by series,plateNumber,scanNumber;");
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        numMosaicEntries = (int)mysql_num_rows(res_ptr);
        mosaicListTable = (PMOSAICENTRY)calloc(numMosaicEntries,sizeof(MOSAICENTRY));
        if (mosaicListTable == NULL) {
          printf("ERROR: failed to allocate mosaic table of size %d\n",numMosaicEntries);
          exit(-1);
        }

        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          pMosaicEntryCur = &mosaicListTable[curMosaicCount];
          memset(pMosaicEntryCur,0,sizeof(MOSAICENTRY));
          if (sqlrow[0]) {
            strcpy(pMosaicEntryCur->series,sqlrow[0]);
          } else {
            continue;
          }
          if (sqlrow[1]) {
            nvals = sscanf(sqlrow[1],"%d",&pMosaicEntryCur->plateNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: nvals is %d for plateNumber\n",nvals);
              continue;
            } 
          } else {
            continue;
          }
          if (sqlrow[2]) {
            nvals = sscanf(sqlrow[2],"%d",&pMosaicEntryCur->scanNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: nvals is %d for scanNumber\n",nvals);
              continue;
            } 
          } else {
            continue;
          }
          if (sqlrow[3]) {
            nvals = sscanf(sqlrow[3],"%d",&pMosaicEntryCur->mosaicNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: nvals is %d for mosaicNumber\n",nvals);
              continue;
            } 
          } else {
            continue;
          }

          curMosaicCount++;
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: mysql_store_result failed in plot_multiple.c line %d\n",__LINE__);
      }
    }
    if (curMosaicCount != numMosaicEntries) {
      printf("ERROR: Current Mosaic Count %d not equal to numMosaicEntries %d\n",curMosaicCount,numMosaicEntries);
    }
    for (mosaicIndex = 0; mosaicIndex < numMosaicEntries; mosaicIndex++) {
      pMosaicEntryCur = &mosaicListTable[mosaicIndex];
      pMosaicEntryNext = &mosaicListTable[mosaicIndex+1];
      if ((mosaicIndex+1 < numMosaicEntries) &&
          (pMosaicEntryCur->plateNumber == pMosaicEntryNext->plateNumber) &&
          (strcmp(pMosaicEntryCur->series,pMosaicEntryNext->series) == 0)) {
        duplicateMosaicCount++;
        continue;
      }
      nonDuplicateMosaicCount++;
      for (solutionCount = 0; solutionCount < MAX_SOLUTIONS; solutionCount++) {
        pMosaic = &mosaicTable[solutionCount];
        if (GetMosaicInfo(pConnection,pMosaicEntryCur->series,pMosaicEntryCur->plateNumber,pMosaicEntryCur->mosaicNumber,solutionCount,pMosaic) == 0) {
          break;
        }
        if ((pMosaic->FitWCS & FITWCS_WEDGE) != 0) {
          wedgeSolutionCount++;
          wedgeSolutionTable[solutionCount]++;
        }

        if ((pMosaic->FitWCS & FITWCS_MULTFAILEDSEPARATION) != 0) {
          separationTable[solutionCount]++;
        }
        if ((pMosaic->FitWCS & FITWCS_STALEEXPOSURENUMBER) != 0) {
          staleTable[solutionCount]++;
        }
        if ((pMosaic->FitWCS & FITWCS_NOALLOBJECTS) != 0) {
          noAllobjects[solutionCount]++;
        }
        if ((pMosaic->FitWCS & FITWCS_DUPLICATEREFCOUNT) != 0) {
          duplicateRef[solutionCount]++;
        }

        if (pMosaic->exposureNumber < 0) {
          unmatchedSolutionCount++;
          unmatchedSolutionTable[solutionCount]++;
          /* The previous mosaic must not have been a wedge */
          wedgeFlag = 0;
          if (solutionCount > 0) {
            pMosaic2 = &mosaicTable[solutionCount-1];
            if ((pMosaic2->FitWCS & FITWCS_WEDGE) != 0) {
              wedgeFlag = 1;
            }
          }
#if 1
          if (wedgeFlag == 0) { 
            goodUnmatchedCount++;
            goodUnmatchedTable[solutionCount]++;
          }
#else
          if ((wedgeFlag == 0) && 
              ((pMosaic->FitWCS & FITWCS_NOALLOBJECTS) == 0) &&
              ((pMosaic->FitWCS & FITWCS_DUPLICATEREFCOUNT) != 0)) {
            goodUnmatchedCount++;
            goodUnmatchedTable[solutionCount]++;
          }
#endif
        }
        totalSolutions++;
      }
      if ((solutionCount <= 0) || (solutionCount >= MAX_SOLUTIONS)) {
        printf("ERROR: MAX_SOLUTIONS %d exceeded\n",solutionCount);
        exit(-1);
      }
      totalSolutionTable[solutionCount-1]++;
      for (exposureCount = 0; exposureCount < MAX_CATALOG_EXPOSURES; exposureCount++) {
        pExposure = &exposureTable[exposureCount];
        if (GetExposureInfo(pConnection,pMosaicEntryCur->series,pMosaicEntryCur->plateNumber,exposureCount,pExposure,1) == 0) {
          break;
        }
        totalExposures++;
      }

      if ((exposureCount <= 0) || (exposureCount >= MAX_CATALOG_EXPOSURES)) {
        printf("ERROR: MAX_CATALOG_EXPOSURES %d %d exceeded\n",MAX_CATALOG_EXPOSURES,exposureCount);
        exit(-1);
      }
      totalExposureTable[exposureCount]++;
      if (exposureCount > solutionCount) {
        missedExposureCount++;
        if (verbose ) {
          printf("Missed: %s%05d\n",pMosaicEntryCur->series,pMosaicEntryCur->plateNumber);
        }
        missedExposureTable[exposureCount] += (exposureCount-solutionCount);
      }
      if (exposureCount == solutionCount) {
        goodExposureCount++;
        if (exposureCount > 1) {
          if (verbose) {
            printf("Good:   %s%05d\n",pMosaicEntryCur->series,pMosaicEntryCur->plateNumber);
          }
        }
        goodExposureTable[exposureCount]++;
      }


    }
  } /* Check Mosaics */

  if (doPlots) {
    sprintf(outvectorfile,"%svector.ps/cps",prefix);
    cpgopen(outvectorfile);
    cpgslw(4); /* line width is 4/200 inch */
    cpgsch(1.5); /* Character size is 1.5/40 the height of the view surface */

    cpgenv(PLOT_XMAX/HOUR_CONVERT,PLOT_XMIN/HOUR_CONVERT,PLOT_YMIN,PLOT_YMAX,0,0);
    cpglab("R.A. (Hours)","Dec (Degrees)","Multiple Exposure Vector Plot");
    sprintf(plotString,"Total Plates: %d",maxMultipleExposures);
    cpgmtxt("T",1.0,0.0,0.0,plotString);

    cpgslw(1); /* line width is 1/200 inch */
    cpgsci(4); /* Set blue */
    plotCount = 0;
    for (multIndex = 0; multIndex < maxMultipleExposures; multIndex++) {
      pMultipleEntry = &pMultipleTable[multIndex];
      pointCount = 0;
      plotFlag = 1;
      for (exposureCount = 0; exposureCount < MAX_CATALOG_EXPOSURES; exposureCount++) {
        pExposure = &exposureTable[exposureCount];
        if (GetExposureInfo(pConnection,pMultipleEntry->series,pMultipleEntry->plateNumber,exposureCount,pExposure,1) == 0) {
          break;
        }
        memcpy(&pMultipleEntry->exposure[pMultipleEntry->exposureCount],pExposure,sizeof(EXPOSURE));
        pMultipleEntry->exposureCount++;
        if ((pExposure->dRightAscension < 990.0) &&
            (pExposure->dDeclination < 98.0)) {
          if ((pExposure->dDeclination > POLE_BIN) ||
              (pExposure->dDeclination < -POLE_BIN)) {
            pExposure->dRightAscension = 180.0;
          }
#if 0
          for (regionIndex = 0; regionIndex < regionTableSize; regionIndex++) {
            pRegion = &regionTable[regionIndex];
            if ((pExposure->dRightAscension >= pRegion->minRA) &&
                (pExposure->dRightAscension <= pRegion->maxRA) &&
                (pExposure->dDeclination >= pRegion->minDec) &&
                (pExposure->dDeclination <= pRegion->maxDec)) {
              /* one of the points lies in an excluded region */
              plotFlag = 0;
            }
          }
#endif
          xval[pointCount] = pExposure->dRightAscension/HOUR_CONVERT;
          yval[pointCount] = pExposure->dDeclination;
          if (pointCount > 0) {
            factor = cos(DEGREES_TO_RAD*((pExposure->dDeclination + oldDeclination)/2.0));
            distance = sqrt(sqr(factor*(pExposure->dRightAscension-oldRightAscension))-sqr(pExposure->dDeclination-oldDeclination));
            if (distance < MINIMUM_DISTANCE) {
              /* Distance is too short. Discard this point. */
              pExposure->dRightAscension = oldRightAscension;
              pExposure->dDeclination = oldDeclination;
              pointCount--;
            }
          }

          pointCount++;
          oldRightAscension = pExposure->dRightAscension;
          oldDeclination = pExposure->dDeclination;
        }
      }
      if ((plotFlag != 0) && (pointCount > 1)) {
        plotCount++;
        cpgline(pointCount,xval,yval);

      }
     
      
    }
    cpgslw(4); /* line width is 4/200 inch */
    cpgsci(1); /* Set black */
    sprintf(plotString,"Plotted Plates: %d",plotCount);
    cpgmtxt("T",1.0,1.0,1.0,plotString);

    cpgiden();
    cpgend();

    /* Now repeat this process iteratively, */
    iteration = 0;
    workDone = 1;
    while (workDone == 1) {
      workDone = 0;
      iteration++;
    
      plotCount = 0;
      for (regionIndex = 0; regionIndex <= REGION_BINS; regionIndex++) {
        pRegion = &regionTable[regionIndex];
        if (pRegion->selected == 0) {
          pRegion->endpointcount = 0;
        }
      }
     
      for (multIndex = 0; multIndex < maxMultipleExposures; multIndex++) {
        pMultipleEntry = &pMultipleTable[multIndex];
        pointCount = 0;
        plotFlag = 1;
        for (exposureCount = 0; exposureCount < pMultipleEntry->exposureCount; exposureCount++) {
          pExposure = &pMultipleEntry->exposure[exposureCount];
          if ((pExposure->dRightAscension < 990.0) &&
              (pExposure->dDeclination < 98.0)) {
            if ((pExposure->dDeclination > POLE_BIN) ||
                (pExposure->dDeclination < -POLE_BIN)) {
              pExposure->dRightAscension = 180.0;
            }
            rejectFlag = 0;
            for (regionIndex = 0; regionIndex < REGION_BINS; regionIndex++) {
              pRegion = &regionTable[regionIndex];
              if (pRegion->selected == 0) {
                continue;
              }
              if ((pExposure->dRightAscension >= pRegion->minRA) &&
                  (pExposure->dRightAscension < pRegion->maxRA) &&
                  (pExposure->dDeclination >= pRegion->minDec) &&
                  (pExposure->dDeclination < pRegion->maxDec)) {
                /* one of the points lies in an excluded region */
#if 1
                rejectFlag = 1;
#else
                plotFlag = 0;
#endif
              }
            }
            if (rejectFlag == 0) {
              xval[pointCount] = pExposure->dRightAscension/HOUR_CONVERT;
              yval[pointCount] = pExposure->dDeclination;
              if (pointCount > 0) {
                factor = cos(DEGREES_TO_RAD*((pExposure->dDeclination + oldDeclination)/2.0));
                distance = sqrt(sqr(factor*(pExposure->dRightAscension-oldRightAscension))+sqr(pExposure->dDeclination-oldDeclination));
                if (distance < MINIMUM_DISTANCE) {
                  /* Distance is too short. Discard this point. */
                  pExposure->dRightAscension = oldRightAscension;
                  pExposure->dDeclination = oldDeclination;
                  pointCount--;
                }
              }

              pointCount++;
              oldRightAscension = pExposure->dRightAscension;
              oldDeclination = pExposure->dDeclination;
            }
          }
        }
        if ((plotFlag != 0) && (pointCount > 1)) {
          plotCount++;
          for (pointIndex = 0; pointIndex < pointCount; pointIndex++) {
            decIndex = (DEC_BINS-1)*(yval[pointIndex]+ 90.5)/(1.0*(DEC_BINS-1));
            raIndex = (RA_BINS-1)*((HOUR_CONVERT*xval[pointIndex])+0.5)/(1.0*(RA_BINS-1));
            if (decIndex < 0) {
              decIndex = 0;
            }
            if (decIndex > DEC_BINS) {
              decIndex = DEC_BINS;
            }
            if (raIndex < 0) {
              raIndex = 0;
            }
            if (raIndex > RA_BINS) {
              raIndex = RA_BINS;
            }
            regionIndex = raIndex+(RA_BINS*decIndex);
            pRegion = &regionTable[regionIndex];
            if (pRegion->selected == 0) {
              pRegion->endpointcount++;
            } else {
              printf("ERROR: selection problem\n");
              exit(-1);
            }
          }

        }


      } /* Multiple Exposure Plate Loop */
      /* Now select the region with the highest endpoint count */
      maxEndpointCount = 0;
      for (regionIndex = 0; regionIndex <= REGION_BINS; regionIndex++) {
        pRegion = &regionTable[regionIndex];
        if ((pRegion->selected == 0) &&
            (pRegion->endpointcount > maxEndpointCount)) {
          maxEndpointCount = pRegion->endpointcount;
          if (topEndpointCount < maxEndpointCount) {
            topEndpointCount = maxEndpointCount;
          }
        }
        
      }
      
      for (regionIndex = 0; regionIndex <= REGION_BINS; regionIndex++) {
        pRegion = &regionTable[regionIndex];
        if ((pRegion->selected == 0) &&
            (pRegion->endpointcount == maxEndpointCount)) {
          pRegion->selected = 1;
          workDone = 1;
        }
        
      }
      time(&curTime);
      curTime -= startTime;
      printf("iteration %5d, maxEndpointCount %4d, plotCount %4d seconds %3d\n",iteration,maxEndpointCount,plotCount,curTime);
      if (maxEndpointCount <= 2) {
        break;
      }
      
    } /* While loop */
    


    sprintf(outvectorfile,"%svector3.ps/cps",prefix);
    cpgopen(outvectorfile);
    cpgslw(4); /* line width is 4/200 inch */
    cpgsch(1.5); /* Character size is 1.5/40 the height of the view surface */
#if 0
    memset(PLOTOUT,0,sizeof(PLOTOUT));
#endif
    pointCount = 0;
    cpgenv(PLOT_XMAX/HOUR_CONVERT,PLOT_XMIN/HOUR_CONVERT,PLOT_YMIN,PLOT_YMAX,0,0);
    cpglab("R.A. (Hours)","Dec (Degrees)","Distribution of Multiple Exposure Fields");
    for (decIndex = 0; decIndex < DEC_BINS; decIndex++) {
      for (raIndex = 0; raIndex < RA_BINS; raIndex++) {
        regionIndex = raIndex + (RA_BINS * decIndex);
        if ((regionIndex < 0) || (regionIndex >= REGION_BINS)) {
          printf("ERROR: regionIndex %d exceeds %d\n",regionIndex,REGION_BINS);
          exit(-1);
        }
        pRegion = &regionTable[regionIndex];
#if 1
#if 0
        xval[pointCount] = pRegion->centerRA/HOUR_CONVERT;
        yval[pointCount] = pRegion->centerDec;
        pointCount++;
#else
        if ((pRegion->selected != 0) &&
            (pRegion->endpointcount >= 2)) {
          xval[pointCount] = pRegion->centerRA/HOUR_CONVERT;
          yval[pointCount] = pRegion->centerDec;
          pointCount++;
          sprintf(locationText,"%d",pRegion->endpointcount);
        }
#endif
#endif
#if 0

        if ((pRegion->selected != 0) &&
            (pRegion->endpointcount > 0)) {
          sprintf(locationText,"%d",pRegion->endpointcount);
          cpgptxt(pRegion->centerRA/HOUR_CONVERT,pRegion->centerDec,0.0,0.5,locationText);
        }
#endif
#if 0
        if (pRegion->selected) {
          PLOTOUT[decIndex][raIndex] = pRegion->endpointcount;
        } else {
          PLOTOUT[decIndex][raIndex] = 0;
        }
#endif
        if ((pRegion->selected != 0) &&
            (pRegion->endpointcount > 0)) {
          printf("ra %8.2f dec %8.2f endpointcount %5d\n",
                 pRegion->centerRA,
                 pRegion->centerDec,
                 pRegion->endpointcount);
      
        }
      }
    }
    cpgpt(pointCount,xval,yval,1);
    TR[0]=-0.5;
    TR[1]=1.0/HOUR_CONVERT;
    TR[2]=0;
    TR[3]=-90.5;
    TR[4]=0;
    TR[5]=1.0;

#if 0
    cpggray ((float *)&PLOTOUT[0][0],RA_BINS,DEC_BINS,1,RA_BINS-1,1,DEC_BINS-1,topEndpointCount/4,0,TR);
    cpgwedg("RG", 1.0, 4.0,topEndpointCount/2,0, "");
#endif
    cpgiden();
    cpgend();

    sprintf(outvectorfile,"%svector2.ps/cps",prefix);
    cpgopen(outvectorfile);
    cpgslw(4); /* line width is 4/200 inch */
    cpgsch(1.5); /* Character size is 1.5/40 the height of the view surface */

    cpgenv(PLOT_XMAX/HOUR_CONVERT,PLOT_XMIN/HOUR_CONVERT,PLOT_YMIN,PLOT_YMAX,0,0);
    cpglab("R.A. (Hours)","Dec (Degrees)","Purged Multiple Exposure Vector Plot");
    sprintf(plotString,"Total Plates: %d",maxMultipleExposures);
    cpgmtxt("T",1.0,0.0,0.0,plotString);

    cpgslw(1); /* line width is 1/200 inch */
    cpgsci(4); /* Set blue */
    for (multIndex = 0; multIndex < maxMultipleExposures; multIndex++) {
      pMultipleEntry = &pMultipleTable[multIndex];
      pointCount = 0;
      plotFlag = 1;
      for (exposureCount = 0; exposureCount < pMultipleEntry->exposureCount; exposureCount++) {
        pExposure = &pMultipleEntry->exposure[exposureCount];
        if ((pExposure->dRightAscension < 990.0) &&
            (pExposure->dDeclination < 98.0)) {
          if ((pExposure->dDeclination > POLE_BIN) ||
              (pExposure->dDeclination < -POLE_BIN)) {
            pExposure->dRightAscension = 180.0;
          }
          rejectFlag = 0;
          for (regionIndex = 0; regionIndex < REGION_BINS; regionIndex++) {
            pRegion = &regionTable[regionIndex];
            if (pRegion->selected == 0) {
              continue;
            }
            if ((pExposure->dRightAscension >= pRegion->minRA) &&
                (pExposure->dRightAscension < pRegion->maxRA) &&
                (pExposure->dDeclination >= pRegion->minDec) &&
                (pExposure->dDeclination < pRegion->maxDec)) {
              /* one of the points lies in an excluded region */
#if 1
              rejectFlag = 1;
#else
              plotFlag = 0;
#endif
            }
          }
          if (rejectFlag == 0) {
            xval[pointCount] = pExposure->dRightAscension/HOUR_CONVERT;
            yval[pointCount] = pExposure->dDeclination;
            if (pointCount > 0) {
              factor = cos(DEGREES_TO_RAD*((pExposure->dDeclination + oldDeclination)/2.0));
              distance = sqrt(sqr(factor*(pExposure->dRightAscension-oldRightAscension))+sqr(pExposure->dDeclination-oldDeclination));
              if (distance < MINIMUM_DISTANCE) {
                /* Distance is too short. Discard this point. */
                pExposure->dRightAscension = oldRightAscension;
                pExposure->dDeclination = oldDeclination;
                pointCount--;
              }
            }

            pointCount++;
            oldRightAscension = pExposure->dRightAscension;
            oldDeclination = pExposure->dDeclination;
          }
        }
      }
      if ((plotFlag != 0) && (pointCount > 1)) {
         
        plotCount++;
        cpgline(pointCount,xval,yval);

      }


    } /* Multiple Exposure Plate Loop */
    cpgslw(4); /* line width is 4/200 inch */
    cpgsci(1); /* Set black */
    sprintf(plotString,"Plotted Plates: %d",plotCount);
    cpgmtxt("T",1.0,1.0,1.0,plotString);

    cpgiden();
    cpgend();

    /* Now dump the region table */
    qsort((void *)regionTable,regionTableSize,sizeof(REGION),RegionCompare);
    fprintf(outregionhandle,"REGION regionTable[] = {\n");

    for (regionIndex = 0; regionIndex < REGION_BINS; regionIndex++) {
      pRegion = &regionTable[regionIndex];
      if (pRegion->selected == 0) {
        continue;
      }
      fprintf(outregionhandle,"  { %5.1f , %5.1f , %5.1f , %5.1f , %5.1f , %5.1f , %6d , %d },\n",
              pRegion->minRA,
              pRegion->maxRA,
              pRegion->minDec,
              pRegion->maxDec,
              pRegion->centerRA,
              pRegion->centerDec,
              pRegion->endpointcount,
              pRegion->selected);
    }
    fprintf(outregionhandle,"};\n");


    
  } /* doPlots */


  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds\n",curTime);
  if (doPlots) {
    printf("Plot results for %d plates  %d multiple exposure plates year: %f-%f\n",plateCount,multipleCount,minYear,maxYear);
    printf("Plates with known class %d, grating plates %d\n",classCount,gratingCount);
    printf("maxExposureNumber %d\n",maxExposureNumber);
  } /* doPlots */
  if (checkMosaics) {
    printf("Total mosaics: %d duplicate mosaics: %d non-duplicate mosaics %d \n",numMosaicEntries,duplicateMosaicCount,nonDuplicateMosaicCount);
    printf("Total solutions %d total exposures %d\n",totalSolutions,totalExposures);

    printf("missedExposureCount %d goodExposureCount %d\n",missedExposureCount,goodExposureCount);
    for (exposureIndex = 0; exposureIndex < MAX_CATALOG_EXPOSURES; exposureIndex++) {
      if ((missedExposureTable[exposureIndex] != 0) || 
          (totalExposureTable[exposureIndex] != 0) ||
          (goodExposureTable[exposureIndex] != 0)) {
        printf("exposureIndex %d, missedExposureTable %3d goodExposureTable %5d totalExposureTable %5d\n",exposureIndex,missedExposureTable[exposureIndex],goodExposureTable[exposureIndex],totalExposureTable[exposureIndex]);
      }
    }

    printf("unmatchedSolutionCount %d goodUnmatchedCount %d wedgeSolutionCount %d\n",unmatchedSolutionCount,goodUnmatchedCount,wedgeSolutionCount);
    for (solutionIndex = 0; solutionIndex < MAX_SOLUTIONS; solutionIndex++) {
      printf("solutionIndex %d, totalSolution %5d unmatchedSolution %3d goodUnmatched %3d wedgeSolution %3d separation %3d stale %3d noAllobjects %3d duplicateRef %3d\n",solutionIndex,totalSolutionTable[solutionIndex],unmatchedSolutionTable[solutionIndex],goodUnmatchedTable[solutionIndex],wedgeSolutionTable[solutionIndex],separationTable[solutionIndex],staleTable[solutionIndex],noAllobjects[solutionIndex],duplicateRef[solutionIndex]);
    }
  } /* checkMosaics */




    /* Done with the MySQL connection */



  mysql_close(pPhotConnection);
  mysql_close(pConnection);


  if (mosaicListTable) {
    free(mosaicListTable);
  }
  if (pMultipleTable) {
    free(pMultipleTable);
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }
  

  if (outscalehandle != NULL) {
    fclose(outscalehandle);
  }
  if (outregionhandle != NULL) {
    fclose(outregionhandle);
  }
  if (outtotalhandle != NULL) {
    fclose(outtotalhandle);
  }
  if (outserieshandle != NULL) {
    fclose(outserieshandle);
  }
  if (outgratinghandle != NULL) {
    fclose(outgratinghandle);
  }

  return(EXIT_SUCCESS);
}

