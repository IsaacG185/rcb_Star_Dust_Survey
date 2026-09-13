// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 * plotmisc.c - Do miscellaneous plots for the web site.
 * 
 * gcc -ggdb -O2  -I/usr/include/plplot -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include plotmisc.c pipelineutils.a -L /dasch/install/lib -lplplotd  -lm  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread    -ltable -lutil  -lwcs -lgd  -o plotmisc
 * 
 *  NOTE: to use plotmisc on Odyssey, switch the loaded plplot version:
 *  module unload hpc/plplot-5.9.4
 *  module load module load  hpc/plplot-5.9.9
 *
 * Extract limiting magnitudes and create a limiting magnitude histogram:
 *         plotmisc -t 'Limiting Magnitude Distribution' -l -i /home/scanner/Pipeline/los1.db  -o /home/scanner/scanner/linux/gmt/limiting_histogram.png  
 *         plotmisc -t 'Limiting Magnitude Distribution' -l -i /home/scanner/backup/2014_04_10/los.db -o   /home/scanner/backup/2014_04_10/limiting_histogram.png
 *
 * Plot the monthly scan graph and scanning rate graph
 *         plotmisc -t 'DASCH Scanning Progress' -s -i /dasch/scanner/linux/gmt/plotpngscanned.txt -o /home/scanner/web/dasch/scangraph.png
 *         plotmisc -t 'DASCH Scanning Rates' -m -i /dasch/scanner/linux/gmt/pngminimum.txt -o /home/scanner/web/dasch/plotminimum.png
 *         plotmisc -t 'DASCH Smoothed Scanning Rate' -r 90 -i /dasch/scanner/linux/gmt/plotpngscanned.txt -o /dasch/scanner/linux/gmt/scanrate.png
 *
 * Plot dasch public site use:
 *
 *

/home/scanner/backup/2013_04_29/dasch_old.log         12/07/10 15:57:43 to 11/02/12 1:31:39
/home/scanner/backup/2013_05_02/dasch2013_05_02.log   05/01/13 14:15:31 to 05/02/13 13:59:39
/home/scanner/backup/2013_06_06/dasch.log             05/02/13 15:52:35 to 06/06/13 13:03:45
/home/scanner/backup/2013_10_03/dasch2013_10_06.log   06/07/13 10:49:12 to 10/06/13 1:53:35 
/home/scanner/backup/2014_08_05/dasch2014_08_06.log   10/06/13 20:22:51 to 08/06/14 11:00:18
/home/scanner/backup/2015_03_19/dasch2015_03_25.log   08/06/14 14:04:26 to 12/08/14 15:32:24 
/home/scanner/backup/2015_07_30/dasch2015_08_03.log   03/25/15 10:48:27 to 08/03/15 17:26:20
/home/scanner/backup/2015_10_22/dasch2015_10_22.log   08/04/15 3:47:58 to  10/22/15 11:36:59
/home/scanner/backup/2016_03_31/dasch2016_04_01.log   11/02/15 18:57:17 to 04/01/16 10:17:08
/home/scanner/backup/2016_08_18/dasch2016_08_19.log   04/01/16 10:50:09 to 08/19/16  9:43:14
/home/scanner/backup/2017_03_09/dasch2017_03_11.log   11/29/16 13:04:36 to 03/10/17 19:03:07
/home/scanner/backup/2017_07_03/dasch2017_07_06.log   03/11/17 14:20:03 to 07/06/17  7:14:23
/home/scanner/backup/2017_09_14/dasch2017_09_15.log   07/06/17 12:33:02 to 09/15/17 17:13:21
/home/scanner/backup/2018_01_25/dasch2018_01_26.log   07/06/17 12:33:02 to 12/09/17  6:10:23  error: need to add extra files
/home/scanner/backup/2018_09_06/dasch2018_09_07.log   01/26/18 15:20:55 to 09/07/18  8:43:09 
/home/scanner/backup/2019_03_14/dasch2019_03_20.log   09/07/18  9:58:58 to 03/20/19 16:26:07
/home/scanner/backup/2019_07_18/dasch2019_07_23.log   03/20/19 18:53:58 to 07/23/19  6:21:39 
/home/scanner/backup/2019_09_05/dasch2019_09_09.log   07/23/19 11:04:37 to 09/09/19 12:27:44
/home/scanner/backup/2020_01_09/dasch2020_01_10.log   09/09/19 16:28:04 to 12/16/19  5:43:01 error: data lost from 12/16/19 through 01/10/20 
/home/scanner/backup/2020_02_06/dasch2020_02_06.log   01/10/20 17:53:02 to 02/07/19  9:48:42  exclude  pc140.astro.lu.se,84.217.57.67 84-217-57-67.customers.ownit.se,eduroam42.nat.wireless.lu.se, eduroam-employee-10.wireless.lu.se
/home/scanner/backup/2020_05_27/dasch2020_05_28.log   02/07/20 11:43:09 to 05/28/20 11:14:20
/home/scanner/backup/2020_09_26/dasch2020_09_30.log   05/28/20 13:34:14 to 09/30/20 19:20:59
/home/scanner/backup/2021_01_30/dasch2021_01_11.log   09/30/20 19:54:09 to 01/11/21 20:03:36 

cd /home/scanner/backup
cat  ./2013_05_02/dasch2013_05_02.log  ./2013_06_06/dasch.log ./2013_10_03/dasch2013_10_06.log ./2014_08_05/dasch2014_08_06.log ./2015_03_19/dasch2015_03_25.log ./2015_07_30/dasch2015_08_03.log ./2015_10_22/dasch2015_10_22.log ./2016_03_31/dasch2016_04_01.log ./2016_08_18/dasch2016_08_19.log ./2017_03_09/dasch2017_03_11.log  ./2017_07_03/dasch2017_07_06.log ./2017_09_14/dasch2017_09_15.log ./2018_01_25/dasch2018_01_26.log ./2018_09_06/dasch2018_09_07.log ./2019_03_14/dasch2019_03_20.log ./2019_07_18/dasch2019_07_23.log ./2019_09_05/dasch2019_09_09.log ./2020_01_09/dasch2020_01_10.log ./2020_02_06/dasch2020_02_06.log ./2020_05_27/dasch2020_05_28.log /home/scanner/backup/2020_09_26/dasch2020_09_30.log /home/scanner/backup/2021_01_30/dasch2021_01_11.log  > /home/scanner/junk/daschcum2021_01_11.log

-w /home/scanner/junk/losweb2020_09_30.log formate
Initial records: fullDate address webType 1 = search, 2 = cat_input
Final records:   pointIndex address totalCount acceptedFlag goodFlag


1,239 references by a single website 151.188.17.247 which consisted of downloading a large number of Kepler Input Catalog stars.  
1,178 references by a single website 149.154.65.151  2017-09-19T7:40:35 to 2017-09-19T23:21:18 (performance test? or denial of service attack?)
NOTE: Public Website was down from Thu Dec 15 14:49:31 2016 to Wed Jan 04 09:22:53 2017 because of an attack  61534 references
LUNDS_UNIVERSITY: 75774 requests from a bot

plotmisc -t 'DASCH Monthly Web Use' -w /home/scanner/junk/losweb2021_01_11.log -i /home/scanner/junk/daschcum2021_01_11.log -o /home/scanner/junk/webuse2021_01_11.png
eog /home/scanner/junk/webuse2021_01_11.png



 *
 * For the median lightcurve vs RMS plot
 *         Use the output from idstats.c for each of the id tables

           cat /home/scanner/junk/DR4/id_keplerDR4.db  /home/scanner/junk/DR4/id_gscDR4.db  /home/scanner/junk/DR4/id_apassDR4.db  > /home/scanner/junk/DR4/id_all.db
           plotmisc -t 'DR4 Median Lightcurve RMS for each Calibration Catalog' -a  -i  /home/scanner/junk/DR4/id_all.db -o /home/scanner/junk/DR4/median_rms.png

           cat /home/scanner/junk/DR4/id_keplerDR4a.db  /home/scanner/junk/DR4/id_gscDR4a.db  /home/scanner/junk/DR4/id_apassDR4a.db  > /home/scanner/junk/DR4/id_alla.db
           plotmisc -t 'DR4 Median Lightcurve RMS for each Calibration Catalog' -a  -i  /home/scanner/junk/DR4/id_alla.db -o /home/scanner/junk/DR4/median_rmsa.png
 
           cat /home/scanner/junk/DR4/id_keplerDR4b.db  /home/scanner/junk/DR4/id_gscDR4b.db  /home/scanner/junk/DR4/id_apassDR4b.db  > /home/scanner/junk/DR4/id_allb.db
           plotmisc -t 'DR4 Median Lightcurve RMS for each Calibration Catalog' -a  -i  /home/scanner/junk/DR4/id_allb.db -o /home/scanner/junk/DR4/median_rmsb.png
 
 *
 * For the Median Points per Lightcurve vs Magnitude plot
 *

           cat /home/scanner/junk/DR4/id_keplerDR4.db  /home/scanner/junk/DR4/id_gscDR4.db  /home/scanner/junk/DR4/id_apassDR4.db  > /home/scanner/junk/DR4/id_all.db
           plotmisc -t 'DR4 Median Points per Magnitude for each Calibration Catalog' -d  -i  /home/scanner/junk/DR4/id_all.db -o /home/scanner/junk/DR4/median_points.png

           cat /home/scanner/junk/DR4/id_keplerDR4a.db  /home/scanner/junk/DR4/id_gscDR4a.db  /home/scanner/junk/DR4/id_apassDR4a.db  > /home/scanner/junk/DR4/id_alla.db
           plotmisc -t 'DR4 Median Points per Magnitude for each Calibration Catalog' -d  -i  /home/scanner/junk/DR4/id_alla.db -o /home/scanner/junk/DR4/median_pointsa.png

           cat /home/scanner/junk/DR4/id_keplerDR4b.db  /home/scanner/junk/DR4/id_gscDR4b.db  /home/scanner/junk/DR4/id_apassDR4b.db  > /home/scanner/junk/DR4/id_allb.db
           plotmisc -t 'DR4 Median Points per Magnitude for each Calibration Catalog' -d  -i  /home/scanner/junk/DR4/id_allb.db -o /home/scanner/junk/DR4/median_pointsb.png


 *
 * For the clipped RMS histograms
 *

     plotmisc -q apass -t 'DR4 clip_rms_local Histogram for the APASS Calibration' -b  -i  /home/scanner/junk/DR4/rms_apassDR4.db -o /home/scanner/junk/DR4/rms_apass_histogram.png
     plotmisc -q gsc2.3.2 -t 'DR4 clip_rms_local Histogram for the GSC2.3.2 Calibration' -b  -i  /home/scanner/junk/DR4/rms_gscDR4.db -o /home/scanner/junk/DR4/rms_gsc_histogram.png
     plotmisc -q kepler -t 'DR4 clip_rms_local Histogram for the KIC Calibration' -b  -i  /home/scanner/junk/DR4/rms_keplerDR4.db -o /home/scanner/junk/DR4/rms_kepler_histogram.png


     plotmisc -q apass -t 'DR4 clip_rms_local Histogram for the APASS Calibration' -b  -i  /home/scanner/junk/DR4/rms_apassDR4a.db -o /home/scanner/junk/DR4/rms_apass_histograma.png
     plotmisc -q gsc2.3.2 -t 'DR4 clip_rms_local Histogram for the GSC2.3.2 Calibration' -b  -i  /home/scanner/junk/DR4/rms_gscDR4a.db -o /home/scanner/junk/DR4/rms_gsc_histograma.png
     plotmisc -q kepler -t 'DR4 clip_rms_local Histogram for the KIC Calibration' -b  -i  /home/scanner/junk/DR4/rms_keplerDR4a.db -o /home/scanner/junk/DR4/rms_kepler_histograma.png

     plotmisc -q apass -t 'DR4 clip_rms_local Histogram for the APASS Calibration' -b  -i  /home/scanner/junk/DR4/rms_apassDR4b.db -o /home/scanner/junk/DR4/rms_apass_histogramb.png
     plotmisc -q gsc2.3.2 -t 'DR4 clip_rms_local Histogram for the GSC2.3.2 Calibration' -b  -i  /home/scanner/junk/DR4/rms_gscDR4b.db -o /home/scanner/junk/DR4/rms_gsc_histogramb.png
     plotmisc -q kepler -t 'DR4 clip_rms_local Histogram for the KIC Calibration' -b  -i  /home/scanner/junk/DR4/rms_keplerDR4b.db -o /home/scanner/junk/DR4/rms_kepler_histogramb.png

 *
 * For the Magnitude measurements per magnitude plot
 *
   cat /home/scanner/junk/DR4/countvsmag_apass_DR4.db  /home/scanner/junk/DR4/countvsmag_gsc_DR4.db /home/scanner/junk/DR4/countvsmag_kepler_DR4.db > /home/scanner/junk/DR4/count_all.db
   plotmisc -t 'DR4 Measurements per Magnitude for Telescope Scale ' -c  -i  /home/scanner/junk/DR4/count_all.db -o /home/scanner/junk/DR4/measurementspermag.png

 *
 * For the Average Limiting Magnitude per Year plot
 *

   cat /home/scanner/junk/DR4/limvsyear_apass_DR4.db  /home/scanner/junk/DR4/limvsyear_gsc_DR4.db  > /home/scanner/junk/DR4/limvsyear_all.db
   plotmisc -t 'DR4 Average Limiting Magnitude per Year ' -e  -i  /home/scanner/junk/DR4/limvsyear_all.db -o /home/scanner/junk/DR4/limperyear.png


*
* For the Magnitude Measurements per Year plot
*

   plotmisc -t 'DR4 Magnitude Measurements per Year ' -f  -i  /home/scanner/junk/DR4/countvsyear_gsc_DR4.db -o /home/scanner/junk/DR4/countperyear.png

* 
*  For linearity tests
* 
   plotmisc -h -i /dasch/data/ExposureData/FlatFrames/FlatStats_2018-05-11T18-39-38.bin -o /dasch/data/ExposureData/FlatFrames/FlatStats_2018-05-11T18-39-38.png

   ls *.bin | awk '{print "plotmisc -h -i",$1,"-o",$1 "XXX"}' | awk '{sub(/binXXX/,"png",$0);print $0}' > los.bat

* 
*  For TC-720 linearity logs
*
   egrep "Apr 04|Apr 05" /home/scanner/junk/TC-720_2019_03_19.txt > /home/scanner/junk/TC-720_2019_03_19A.txt
   plotmisc -n -t 'Scanner CCD Thermal Control' -i /home/scanner/junk/TC-720_2019_03_19A.txt -o /home/scanner/junk/tc720.png

*
*  For plate cleaning analysis
*

echo "select date,pickedStatus,series,plateNumber from plateevents where (pickedStatus = 'Machine Cleaned' or pickedStatus = 'Hand Cleaned') and valid = 'yes' and date > '2017-09-05' order by date;" > /home/scanner/Pipeline/los.sql
mysql stacks < /home/scanner/Pipeline/los.sql > /home/scanner/Pipeline/los.txt
plotmisc  -t 'Plate Cleaning' -j  -i /home/scanner/Pipeline/los.txt -o /home/scanner/Pipeline/los.png
eog /home/scanner/Pipeline/los.png

Note: 2017-06-08T16:30:51 to 2017-06-08T17:00:37 345 plates!!!  CLEANING_MAX_RATE = 100
line 2158 sessionCount 345 hours 0.496451 rate 694.932751
2017-09-05 start of stable machine cleaning
totalHandCleanCount 82939 totalHandCleanHours 2454.556719 aveHandCleaningRate 33.789808
totalMachineCleanCount 30512 totalMachineCleanHours 1179.877582 aveMachineCleaningRate 25.860310


 * Mar 25, 2014 Edward J. Los - Initial version
 * May  6, 2014 Edward J. Los - Support plotting of the monthly scan graph and the scanning rate graph
 * Aug  6, 2014 Edward J. Los - Support ploting of dasch public website use
 * Dec 18, 2014 Edward J. Los - Maximum number of plates in a day is now 464.
 * Apr  6, 2015 Edward J. Los - Correct axes of the DASCH Scanning Rates plot
 * Jun  6, 2015 Edward J. Los - Add -a through -f options
 * Mar 24, 2017 Edward J. Los - Add -g <timestr> option to test the search_none README*.png files
 * May  1, 2017 Edward J. Los - Make the lowest date float for the scangraph (-s)
 *                              Add BASE_YEAR for post-flood graph
 * Nov 22, 2017 Edward J. Los - Correct countperyear pointers
 * May 14, 2018 Edward J. Los - Add graphing for linearity files
 * Jul 22, 2018 Edward J. Los - Add GAIA support
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support 
 * Mar  6, 2019 Edward J. Los - Adjust limits for DR7 support
 * Mar 29, 2019 Edward J. Los - Add a plot option for the TC-720 log file
 * Apr  1, 2019 Edward J. Los - Plot T1 only with different colors
 * Apr  5, 2019 Edward J. Los - Optionally plot the power level and the setpoint
 * Apr 16, 2019 Edward J. Los - Support ATLAS in the lightcurve rms plot
 * Apr 17, 2019 Edward J. Los - Add 'j' option to check performance of the plate cleaning machine
 * May 11, 2019 Edward J. Los - Change label positions for ATLAS
 * May 21, 2019 Edward J. Los - Use a 30 day average for the plate cleaning machine
 * Jul 24, 2019 Edward J. Los - Add code to identify the most frequently used addresses
 *                              Eliminate identical requests in less than 1 month
 */
#define _GNU_SOURCE 
#define __USE_XOPEN2K8 1
#include "table.h"
#include <string.h>
#include <time.h>
#include <math.h>
#include "plplot.h"
#include "plplotP.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
/* #define BASE_YEAR 2017.0 */
#define NUMDEC 181
#define NUMRA  361
#if 0
#define NSHADES 5
#define NCOLORBAR 6
#else
#define NSHADES 100
#define NCOLORBAR 11
#endif
#define NPOINTDEF 100
#define MAX_BUFFER 512
#define MAX_OBJECT_LENGTH 80
#define MIN_EPOCH  (31./365.25) /* One month interval between successive requests */
#define STEP 1
#define MINDEC -90.0
#define MAXDEC +90.0
#define MINRA    0.0
#define MAXRA  360.0
#define DEGREES_TO_RAD  (3.141592654/180.0)
#define RAD_TO_DEGREES  (180.0/3.141592654)
#define MAGNITUDE_MAX 25
#define MAGNITUDE_MULTIPLIER 10
#define HISTOGRAM_SIZE (MAGNITUDE_MAX*MAGNITUDE_MULTIPLIER)
#define YEARS_PER_DAY (1./365.25) 
#define YEARS_PER_WEEK (7./365.25)
#define MIN_POINTS_PER_CURVE_MAG 8.4
#define USE_CENTIGRADE 1
#define LUNDS_UNIVERSITY 1 /* Define this to eliminate the Lunds entries */
/* #define PLOT_DAY  1551330000L */ /* Mar 0, 2019 */
/* #define PLOT_DAY  1554004800L */ /* Apr 0, 2019 */
/* #define PLOT_DAY  1556596800L */ /* May 1, 2019 */
/* #define PLOT_DAY  1559275200L */ /* Jun 1, 2019 */
/* #define PLOT_DAY  1561953600L */ /* Jul 1, 2019 */ 
#define PLOT_DAY     1583035200L    /*  Mar 3, 2020 */

/* #define PLOT_T1_ONLY 1 */
#if !defined(PLOT_T1_ONLY) 
#define PLOT_SETPOINT 1
#define PLOT_POWER 1
#endif /*  !defined(PLOT_T1_ONLY) */


#if 0
#define LOW_CUTOFF 130.0
#define HIGH_CUTOFF 230.0
#else
#define LOW_CUTOFF 145.0
#define HIGH_CUTOFF 210.0
#endif
#define MAX_TIMESTR_SIZE 26

extern GSCBIN gscBin01;
PGSCBIN pGscBin01 = &gscBin01;
extern char *releaseFieldText[RELEASE_FIELD_MAX+1];

#define M_SQRT2 1.41421356237309504880
#define GMT_CONV_LIMIT	1.0e-8	/* Fairly tight convergence limit or "close to zero" limit */
#if 1
#define EQ_RAD (1.0/M_SQRT2)
#else
#define EQ_RAD 6371008.7714
#endif
#define M_PI 3.14159265358979323846

#define MAX_PLATE_STRING 30

typedef struct _limitingstats {
  int count;   /* count of plates */
  int histogram[HISTOGRAM_SIZE+1];
  double dim_limiting_mag_local;
  double bright_limiting_mag_local;
} LIMITINGSTATS,*PLIMITINGSTATS;
#define MAX_PLOT_LABEL 25
typedef struct _labeltable {
  char plotLabel[MAX_PLOT_LABEL];
  double xLabelPos;
  double yLabelPos;
  double xArrowPos;
  int color;
  int style;
  PLIMITINGSTATS statsTable;
} LABELTABLE,*PLABELTABLE;


typedef struct _limitingtable {
  char Plate[MAX_PLATE_STRING+2];  /* Plate identification */
  int versionId;
  int catalogNumber;
  int patrolPlate;
  double limiting_mag_local;
} LIMITINGTABLE,*PLIMITINGTABLE;

typedef struct _scangraph {
  double epoch;
  int count;
  double rate;
  double midepoch;
} SCANGRAPH,*PSCANGRAPH;

typedef struct _minimumgraph {
  double minimum;
  double average;
  double median;
  double pct600;
  double pct700;
  double pct800;
  double pct900;
  double pct950;
  double pct990;
  int count;
} SCANRATE,*PSCANRATE;
#define WEB_TYPE_PLATE_SEARCH  1
#define WEB_TYPE_OBJECT_SEARCH 2
#define MAX_URL_STRING 20
#define RELEASE_DATE 2013.32877 /* May 1 2013  for DR1 */
#define BASE_DATE 2013.000
#define MAX_MONTH (12*12)
typedef struct _webusetable {
  int sequence;
  int webtype;  /* WEB_TYPE_PLATE_SEARCH = 1 and WEB_TYPE_OBJECT_SEARCH = 2 */
  int month;
  double epoch;
  char date[MAX_DATE_STRING];
  char time[MAX_DATE_STRING];
  char url[MAX_URL_STRING];
  char object[MAX_OBJECT_LENGTH];
} WEBUSE,*PWEBUSE;

typedef struct _webstatstable {
  char url[MAX_URL_STRING];
  int totalCount;
  int acceptedFlag; /* 0 = rejected because of known bad address; 1 = accepted */
  int goodFlag;     /* acceptedFlag = 1, goodFlag = 1 if already vetted */
} WEBSTATS,*PWEBSTATS;

typedef struct _lightcurvermstable {
  int minmag;
  int maxmag;
  int count;
  double median;
  double normnburst3;
  int catalogNumber;
  int npoints;
} LIGHTCURVERMS,*PLIGHTCURVERMS;

typedef struct _measpermag {
  int sequence; /* PATROL = 1, NONPATROL = 2, TOTAL = 3 */
  int magcal_magdep;
  int catalogNumber;
  long long count;
} MEASPERMAG,*PMEASPERMAG;

typedef struct _limmagperyear {
  int sequence; /* PATROL = 1, NONPATROL = 2, TOTAL = 3 */
  int year; /* calendar year */
  int catalogNumber;
  double limiting_mag_local;
} LIMMAGPERYEAR,*PLIMMAGPERYEAR;

typedef struct _countperyear {
  int sequence; /* PATROL = 1, NONPATROL = 2, TOTAL = 3 */
  int year; /* calendar year */
  int catalogNumber;
  long long count;
} COUNTPERYEAR,*PCOUNTPERYEAR;
#define HAND_CLEAN 1
#define MACHINE_CLEAN 2
#define CLEAN_IDLE_INTERVAL (0.5/(24.0*365.25)) /* One half hour idle time */
#define CLEAN_AVERAGE_DAYS 30.0 /* Running average */
#define MIN_SESSION_COUNT 4
#define CLEANING_MAX_RATE 100.0
typedef struct _cleaning {
  int sequence; /* hand = 1; machine = 2 */
  int count;
  double epoch;
  int sessionCount;  /* session count */
  double hours; /* session time */
  double rate; /* Plates per hour */
  char date[MAX_DATE_STRING];
  char day[MAX_DATE_STRING];
} CLEANING,*PCLEANING;


typedef struct _tc720table {
  struct tm tmstruct; /* parsed time */
  time_t tp;
  char timestr[100];
  double T1; /* control temperature */
  double T2; /* auxiliary temperature */
  double setpoint; /* setpoint temperature */
  double output;  /* output (-100 to +100) */
  double T3; /* Unknown temperature */
  double hour;
	       /* int tm_sec;	   /* seconds */
	       /* int tm_min;	   /* minutes */
	       /* int tm_hour;	   /* hours */
	       /* int tm_mday;	   /* day of the month */
	       /* int tm_mon;	   /* month */
	       /* int tm_year;	   /* year since 1900 */
	       /* int tm_wday;	   /* day of the week */
	       /* int tm_yday;	   /* day in the year */
	       /* int tm_isdst;	   /* daylight saving time */
} TC720,*PTC720;

char *weekdays[] = {"Sun,","Mon,","Tue,","Wed,","Thu,","Fri,","Sat,"};
char numWeekdays = sizeof(weekdays)/sizeof(char *);
char *monthstr[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
char numMonthstr = sizeof(monthstr)/sizeof(char *);


typedef struct _cliprms {
  int minmag;
  int maxmag;
  double clip_rms_local;
  int count;
  int catalogNumber;
} CLIPRMS,*PCLIPRMS;


char *scanText[] = 
  {"Maximum",
   "99 Percentile",
   "95 Percentile",
   "90 Percentile",
   "80 Percentile",
   "70 Percentile",
   "60 Percentile",
   "Median",
   "Average"};

char *catalogTextUppercase[MAX_CATALOG_NUMBER] =
  {"GSC2.3.2",
   "KIC",
   "APASS",
   "GAIA",
   "ATLAS",
   "experimental"
  };



int LightcurveRMSCompare(const void *first, const void *second) 
{
  int minmagFirst = ((PLIGHTCURVERMS)first)->minmag;
  int minmagSecond = ((PLIGHTCURVERMS)second)->minmag;
  if (minmagFirst > minmagSecond) {
    return(1);
  } else if (minmagFirst < minmagSecond) {
    return(-1);
  } else {
    return(0);
  }

}

int WebStatsDuplicate(const void *first, const void *second) 
{
  PWEBUSE pWebUseFirst = (PWEBUSE)first;
  PWEBUSE pWebUseSecond = (PWEBUSE)second;
  int result;
  char urlfirst[MAX_URL_STRING];
  char urlsecond[MAX_URL_STRING];

  if (pWebUseFirst->webtype > pWebUseSecond->webtype) {
    result = 1;
  } else if (pWebUseFirst->webtype < pWebUseSecond->webtype) {
    result = -1;
  } else {
    strcpy(urlfirst,pWebUseFirst->url);
    strcpy(urlsecond,pWebUseSecond->url);
#if 0
    if ((strstr(urlfirst,"210.72.91.") != NULL) &&
        (strstr(urlsecond,"210.72.91.") != NULL)) {
      urlfirst[10] = 0;
      urlsecond[10] = 0;
    }
#endif
    result = strcmp(urlfirst,urlsecond);
    if (result == 0) {
      result = strcmp(pWebUseFirst->object,pWebUseSecond->object);
      if (result == 0) {
        if (pWebUseFirst->epoch > pWebUseSecond->epoch) {
          result = 1;
        } else if (pWebUseFirst->epoch < pWebUseSecond->epoch) {
          result = -1;
        } else {
          result = 0;
        }
      }
    }
  }
  return(result);

}
int WebStatsSequence(const void *first, const void *second) 
{
  PWEBUSE pWebUseFirst = (PWEBUSE)first;
  PWEBUSE pWebUseSecond = (PWEBUSE)second;
  int result;

  if (pWebUseFirst->sequence > pWebUseSecond->sequence) {
    result = 1;
  } else if (pWebUseFirst->sequence < pWebUseSecond->sequence) {
    result = -1;
  } else {
    result = 0;
  }
  return(result);

}
int WebStatsCompare(const void *first, const void *second) 
{
  PWEBSTATS webstatsFirst = (PWEBSTATS)first;
  PWEBSTATS webstatsSecond = (PWEBSTATS)second;
  int result;
  result = strcmp(webstatsFirst->url,webstatsSecond->url);
  return(result);

}
int WebStatsTotal(const void *first, const void *second) 
{
  PWEBSTATS webstatsFirst = (PWEBSTATS)first;
  PWEBSTATS webstatsSecond = (PWEBSTATS)second;
  /* Reverse sort: largest first */
  if (webstatsFirst->totalCount < webstatsSecond->totalCount) {
    return(1);
  } else if (webstatsFirst->totalCount > webstatsSecond->totalCount) {
    return(-1);
  }
  return(0);

}


int main(int argc,char *argv[])
{
  PLIMITINGTABLE limitingTable = NULL;
  PLIMITINGTABLE pLimiting;
  PSCANGRAPH scangraphTable = NULL;
  PSCANGRAPH pScangraph;
  PSCANGRAPH pScangraph1;
  PSCANGRAPH pScangraph2;
  PSCANRATE scanrateTable = NULL;
  PSCANRATE pScanrate;
  PLIGHTCURVERMS lightcurveRMSTable = NULL;
  PLIGHTCURVERMS pLightcurveRMS;
  PCLIPRMS clipRMSTable = NULL;
  PCLIPRMS pClipRMS;
  PMEASPERMAG measPerMagTable = NULL;
  PMEASPERMAG pMeasPerMag;
  PLIMMAGPERYEAR limMagPerYearTable = NULL;
  PLIMMAGPERYEAR pLimMagPerYear;
  PCOUNTPERYEAR countPerYearTable = NULL;
  PCOUNTPERYEAR pCountPerYear;
  PCLEANING cleaningTable = NULL;
  PCLEANING pCleaning;
  PCLEANING pLastCleaning;
  int lastPointIndex;
  int cumulativeSessionCount;
  double cumulativeHoursCount;
  
  int machineCleanCount = 0;
  int handCleanCount = 0;
  int maxCleaningCount = 0;
  double minHandCleanYear = 3000;
  double maxHandCleanYear = 0;



  double startHandCleanYear = -1;
  double lastHandCleanYear = -1;
  int tempHandCleanCount = 0;
  PCLEANING pLastHandCleaning;
  int totalHandCleanCount = 0;
  double totalHandCleanHours = 0;

  double startMachineCleanYear = -1;
  double lastMachineCleanYear = -1;
  int tempMachineCleanCount = 0;
  PCLEANING pLastMachineCleaning;
  int totalMachineCleanCount = 0;
  double totalMachineCleanHours = 0;


  int minmag = 99;
  int maxmag = 0;
  int minyear = 2000;
  int maxyear = 0;
  double minlimitingmag = 99.0;
  double maxlimitingmag = 0.0;
  long long maxcountperyear;
  int curmag;
  int maxnpoints = 0;
  double minmedian = 0;
  double maxmedian = 0;
  long long maxCount;
  int sequence;
  double justify;  /* Justify 0 = reference at left; 1 = reference at right */
  int labelIndex;
  PWEBUSE webuseTable = NULL;
  PWEBUSE pWebUse;
  PWEBUSE pWebUse2;
  int webuseCount = 0;
  PWEBSTATS webstatsTable = NULL;
  PWEBSTATS pWebStats;
  PWEBSTATS pWebStats2;
  int acceptedFlag;
  int nondupCount = 0;
  PTC720 tc720Table = NULL;
  PTC720 ptc720;
  PTC720 nxt720;
  TC720 mintc720;
  TC720 maxtc720;
  char* tabPtr;
	double ra;
	double dec;
	double count;
	int xindex;
	int yindex;
	int tindex = 0;
	double xval;
	double yval;
  char *argstr;
	char cmdchar;
  int verbose = 0;
	int galacticFlag = 0;
	int publicFlag = 0;
	int	releaseLabelFlag = 0;
  int lightcurveRMSFlag = 0;
  int pointsPerCurveFlag = 0;
  int clipRMSFlag = 0;
  int measPerMagFlag = 0;
  int limMagPerYearFlag = 0;
  int countPerYearFlag = 0;
  int plotREADMEFlag = 0;
  FILE *input_handle = NULL;
  char input_name[MAX_BUFFER];
  char outfile[MAX_BUFFER]; 
  char webfile[MAX_BUFFER];
  FILE *webhandle = NULL;
  char magnitudeString[MAX_BUFFER];
	char magnitudeFile[MAX_BUFFER];
	char tempBuffer[MAX_BUFFER];
	char *pDot;
	int errorFlag = 0;
	int complementFlag = 0;
  char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
	int lineCount = 0;
  int lineCount2 = 0;
	int curCount = 0;
	int authorizedCount = 0;
	int nvals;
	int index;
	int gsc_bin_index;
	int decBin;
	int raBin;
	int unitializedCount = 0;
	int authorized;
  time_t baseTime;
	char title[MAX_BUFFER];
  char title2[MAX_BUFFER];
  char xlabel[MAX_BUFFER];
	char colormap_scale[MAX_BUFFER];
	double upperLimit = -1.0;
	double lowerLimit = -1.0;
	int catalogNumber = -1;
	char source[MAX_BUFFER];
	char catalogString[MAX_BUFFER];
	char qualifier[MAX_BUFFER];
  char* catalogdir;
  GALAXYCOMMON galaxycommon;
  PGALAXYCOMMON pGalaxyCommon = &galaxycommon;
	int releaseField;
	double *pLimitingVector;
	double *limitingVector[RELEASE_FIELD_MAX+1][MAGNITUDE_MAX];
	int plateCountLimiting[RELEASE_FIELD_MAX+1][MAGNITUDE_MAX];
	int plateCountMaximum[RELEASE_FIELD_MAX+1];
	int plateCount;
	int releaseFieldIndex;
	int magnitudeIndex;
	PPLATELIMITINGREC pPlateRec;
	int limitingMagLocal;
	int releaseFieldCount[RELEASE_FIELD_MAX+1];
	int magnitudeCount[MAGNITUDE_MAX];
	int brightMagnitudeLimit;
	double lon;
	double lat;
  int scangraphFlag = 0;
  int limitingMagnitudeFlag = 0;
  int scanrateFlag = 0;
  int scanSmoothedRateFlag = 0;
  int smoothingDays = 0;
  int webuseFlag = 0;
  int switchCounter;
  PLABELTABLE pLabelEntry;
  PLIMITINGSTATS pLimitingStats;

  LIMITINGSTATS allPlates;
  PLIMITINGSTATS pAllPlates = &allPlates;

  LIMITINGSTATS gscNonpatrolPlates;
  PLIMITINGSTATS pGscNonpatrolPlates = &gscNonpatrolPlates;

  LIMITINGSTATS gscPatrolPlates;
  PLIMITINGSTATS pGscPatrolPlates = &gscPatrolPlates;

  LIMITINGSTATS apassNonpatrolPlates;
  PLIMITINGSTATS pApassNonpatrolPlates = &apassNonpatrolPlates;

  LIMITINGSTATS apassPatrolPlates;
  PLIMITINGSTATS pApassPatrolPlates = &apassPatrolPlates;

  LIMITINGSTATS atlasNonpatrolPlates;
  PLIMITINGSTATS pAtlasNonpatrolPlates = &atlasNonpatrolPlates;

  LIMITINGSTATS atlasPatrolPlates;
  PLIMITINGSTATS pAtlasPatrolPlates = &atlasPatrolPlates;

  int limitingMagIndex;
  int limitingMagCount;
  int statsIndex;
  int pointIndex;
  int pointIndex1;
  int pointIndex2;
  double currentDate;
  double laterDate;
  double peakCurrentDate;
  double peakLaterDate;
  int goodPoints;
  int goodCurves;
  int peakIndex;
  char *charPtr1;
  char *charPtr2;
  char *charPtr3;
  char charVal;
  char copyDate[MAX_DATE_STRING];
  char fullDate[MAX_DATE_STRING];
  double julianDate;
  double minEpoch = 0.0;
  double maxEpoch = 0.0;
  int month;
  int year;
  int plateSearchCount = 0;
  int lightcurveSearchCount = 0;
  int monthIndex;
  int maxMonthIndex = 0;
  int monthArray[MAX_MONTH];
  int nonharvardArray[MAX_MONTH];
  int nonharvardCount = 0;
  int nonharvardMaxMonth = 0;
  char timestr2[MAX_TIMESTR_SIZE];
  int linearityFlag = 0;
  int tc720Flag = 0;
  int cleaningFlag = 0;
  int quadrant;
  int actualValue;
  int result;
  char timestr[100];
  struct tm *ptr;
  int color;
  int objectLength;

  LINCOMMON linearityCommon;
  PLINCOMMON pLinearityCommon = &linearityCommon;

#if 0
  LABELTABLE statsLabels[] = 
    {
      {"GSC Patrol Plates"      ,16.0,92.0,11.2,14,2,pGscPatrolPlates},
      {"APASS Patrol Plates"    ,16.0,82.0,11.2,13,3,pApassPatrolPlates},
      {"ATLAS PatrolPlates"     ,16.0,72.0,11.2,11,4,pAtlasPatrolPlates},
      {"All Plates"             ,16.0,62.0,11.2,15,1,pAllPlates},
      {"GSC Non-Patrol Plates"  ,16.0,52.0,11.2, 3,2,pGscNonpatrolPlates},
      {"APASS Non-Patrol Plates",16.0,42.0,11.2, 9,3,pApassNonpatrolPlates},
      {"ATLAS Non-Patrol Plates",16.0,32.0,11.2, 4,4,pAtlasNonpatrolPlates},
    };
#else
  LABELTABLE statsLabels[] = 
    {
      {"GSC Patrol Plates"      ,16.0,92.0,13.0,14,2,pGscPatrolPlates},
      {"APASS Patrol Plates"    ,16.0,82.0,13.0,13,3,pApassPatrolPlates},
      {"ATLAS Patrol Plates"    ,11.0,42.0,14.0,11,4,pAtlasPatrolPlates},
      {"All Plates"             ,16.0,62.0,14.0,15,1,pAllPlates},
      {"GSC Non-Patrol Plates"  ,16.0,52.0,15.5, 3,2,pGscNonpatrolPlates},
      {"APASS Non-Patrol Plates",16.0,42.0,15.5, 9,3,pApassNonpatrolPlates},
      {"ATLAS Non-Patrol Plates",11.0,10.0,18.0, 4,4,pAtlasNonpatrolPlates},
    };
#endif

  int labelSize = sizeof(statsLabels)/sizeof(LABELTABLE);

  double xmin;
  double xmax;
  double ymin;
  double ymax;
	double xminmm;
	double xmaxmm;
	double yminmm;
	double ymaxmm;
  double maxcount = 0;
  PLFLT *x = NULL;
	PLFLT *y = NULL;
  PLFLT *T1 = NULL;
  PLFLT *T2 = NULL;
  PLFLT *T3 = NULL;
  PLFLT *OUTPUT = NULL;
  PLFLT *SETPOINT = NULL;
  PLFLT x1[2] = {0,0};
  PLFLT y1[2] = {0,0};
  int just = 0.0;
  int allocCount;

  strcpy(timestr2,"No date");
  
  memset(pAllPlates,0,sizeof(LIMITINGSTATS));
  memset(pGscNonpatrolPlates,0,sizeof(LIMITINGSTATS));
  memset(pGscPatrolPlates,0,sizeof(LIMITINGSTATS));
  memset(pApassNonpatrolPlates,0,sizeof(LIMITINGSTATS));
  memset(pApassPatrolPlates,0,sizeof(LIMITINGSTATS));
  memset(pAtlasNonpatrolPlates,0,sizeof(LIMITINGSTATS));
  memset(pAtlasPatrolPlates,0,sizeof(LIMITINGSTATS));
  memset(monthArray,0,sizeof(monthArray));
  memset(nonharvardArray,0,sizeof(nonharvardArray));
  memset(pLinearityCommon,0,sizeof(LINCOMMON));


	memset(limitingVector,0,sizeof(limitingVector));
	memset(plateCountLimiting,0,sizeof(plateCountLimiting));
	memset(plateCountMaximum,0,sizeof(plateCountMaximum));
	memset(pGalaxyCommon,0,sizeof(GALAXYCOMMON));
	memset(releaseFieldCount,0,sizeof(releaseFieldCount));
	memset(magnitudeCount,0,sizeof(magnitudeCount));

	title[0] = 0;
	colormap_scale[0] = 0;
	input_name[0] = 0;
	outfile[0] = 0;
  webfile[0] = 0;

	qualifier[0] = 0;
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d\n",argstr,argc);
    } else {
      while ((cmdchar = *++argstr) != 0) {

        /* Note: values j,k,p,u,x,y, and z are free */

				switch(cmdchar) {
				case 'i': /* input file */
				case 'I':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for input file -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(input_name,*++argv,MAX_BUFFER-2);
						if (strlen(input_name) >= MAX_BUFFER-3) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
							exit(-1);
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
						if (strlen(outfile) >= MAX_BUFFER-3) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
							exit(-1);
						}
					}
					break;


				case 'w': /* output file name */
				case 'W':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(webfile,*++argv,MAX_BUFFER-2);
						if (strlen(webfile) >= MAX_BUFFER-3) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
							exit(-1);
						}
					}
          webuseFlag = 1;
					break;



        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              strcpy(source,*argv);
              if (catalogNumber > 0) {
                sprintf(catalogString,"%d",catalogNumber);
                strcpy(qualifier,*argv);
              }
            }
          }
          break;

				case 't': /* title */
				case 'T':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for title -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(title,*++argv,MAX_BUFFER-2);
						if (strlen(title) >= MAX_BUFFER-3) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
							exit(-1);
						}
					}
					break;




				case 'v': /* verbose */
				case 'V':
					verbose = 1;
					break;

				case 'l': /* Input starbase for a limiting magnitude plot */
				case 'L':
          limitingMagnitudeFlag = 1;
					break;

				case 'a': /* median lightcurve RMS vs magnitude*/
				case 'A':
          lightcurveRMSFlag = 1;
					break;

				case 'b': /* limiting magnitude histogram */
				case 'B':
          clipRMSFlag = 1;
					break;

				case 'c': /* limiting magnitude histogram */
				case 'C':
          measPerMagFlag = 1;
					break;

				case 'd': /* median lightcurve RMS vs magnitude*/
				case 'D':
          pointsPerCurveFlag = 1;
					break;

				case 'e': /* average limiting magnitude per year*/
				case 'E':
          limMagPerYearFlag = 1;
					break;

				case 'f': /* lightcurve count per year*/
				case 'F':
          countPerYearFlag = 1;
					break;

				case 'g': /* README files for search_none */
				case 'G':
          plotREADMEFlag = 1;
          if ((argc > 1) && (argv[1][0] != '-')) {
            argc--;
            strncpy(timestr2,*++argv,MAX_TIMESTR_SIZE-2);
            timestr2[MAX_TIMESTR_SIZE-1] = 0;
            if (strlen(timestr2) >= MAX_TIMESTR_SIZE-3) {
              printf("ERROR: MAX_TIMESTR_SIZE too small for %s\n",*argv);
              exit(-1);
            }
          }


					break;

				case 's': /* scanning progress flag */
				case 'S':
          scangraphFlag = 1;
					break;

				case 'h': /* plot linearity data */
				case 'H':
          linearityFlag = 1;
					break;

				case 'n': /* TC-720 temperature controller */
				case 'N':
          tc720Flag = 1;
					break;

				case 'j': /* plate cleanint */
				case 'J':
          cleaningFlag = 1;
					break;

				case 'm': /* scanning daily rate flag */
				case 'M':
          scanrateFlag = 1;
					break;

        case 'r': /* plot smoothed rate flag */
 					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for input file -%c\n",cmdchar);
						errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&smoothingDays);
            if (nvals != 1) {
              printf("ERROR: Unable to decode smoothingDays %s\n",*argv);
              errorFlag = 1;
            } else {
              if (smoothingDays < 1) {
                printf("ERROR: smoothingDays %d must be at least 1\n",smoothingDays);
                errorFlag = 1;
              }
            }
					}
          scanSmoothedRateFlag = 1;
          break;


				default:
					printf("ERROR: * illegal command -%c-",cmdchar);
					errorFlag = 1;

				}
        
      }

    }
  }

  if (plotREADMEFlag != 0) {
    printf("Plotting the search_none README files\n");
    PlotSymbolKey(timestr2,0);
    exit(0);
  }

  if (input_name[0] == 0) {
		printf("ERROR: no output file specified\n");
		errorFlag = 1;

  }

	if (outfile[0] == 0) {
		printf("ERROR: no output file specified\n");
		errorFlag = 1;
	}
	if ((title[0] == 0) && (linearityFlag == 0)) {
		printf("WARNING: no title specified\n");
	}

  if ((limitingMagnitudeFlag | scangraphFlag | scanrateFlag | webuseFlag | lightcurveRMSFlag | clipRMSFlag | measPerMagFlag | pointsPerCurveFlag | limMagPerYearFlag | countPerYearFlag | scanSmoothedRateFlag | linearityFlag | tc720Flag | cleaningFlag) == 0) {
    printf("ERROR: neither -l nor -s nor -m nor -a nor -b nor -r nor -n nor -j is specified\n");
    errorFlag = 1;
  }
  if ((limitingMagnitudeFlag + scangraphFlag + scanrateFlag + webuseFlag + lightcurveRMSFlag + clipRMSFlag + measPerMagFlag + pointsPerCurveFlag + limMagPerYearFlag + countPerYearFlag + scanSmoothedRateFlag+linearityFlag+tc720Flag+cleaningFlag) > 1) {
    printf("ERROR: only one of -l (%d) or -s (%d) or -m (%d) or -w (%d) or -a (%d) or -b (%d) or -c (%d) or -d (%d) or -e (%d) or -f (%d) or -r (%d) or -h (%d) or -n (%d) or -j (%d) may be specified\n",limitingMagnitudeFlag,scangraphFlag,scanrateFlag,webuseFlag,lightcurveRMSFlag,clipRMSFlag,measPerMagFlag,pointsPerCurveFlag,limMagPerYearFlag,countPerYearFlag,scanSmoothedRateFlag,linearityFlag,tc720Flag,cleaningFlag);
    errorFlag = 1;
  }

	if (errorFlag) {
		printf("Usage: plotmisc [-q <catalog>|-i<input file>] options\n");
    printf("  options: -v verbose\n");
    printf("           -i <input file>\n");
		printf("           -q <calibration catalog> to plot limiting magnitudes\n");
    printf("           -o <output file>\n");
    printf("           -l Plot the limiting magnitude graph\n");
    printf("           -s Plot the scanning progress graph\n");
    printf("           -m Plot the scanning daily rate graph\n");
    printf("           -r <days> Plot the smoothed scanning rate graph with <days> average\n");
    printf("           -w <filename> Plot the website graph from dasch.log put url values in <filename>\n");
    printf("           -a Plot the Median Lightcurve RMS vs Magnitude graph\n");
    printf("           -b Plot the clip_rms_local histogram\n");
    printf("           -g <timestr> Plot the search_none README files\n");
    printf("           -h plot input as a linearity file\n");
    printf("           -n plot input as a TC-720 temperature controller log file\n");
		printf("           -t title\n");
    return(-1);
	}		
	printf("plotmisc of %s %s complementFlag %d\n",__DATE__,__TIME__,complementFlag);

	if (input_name[0] != 0) {
    if (linearityFlag != 0) {
      result = ReadLinearityFile(pLinearityCommon,input_name);
      if (result < 0) {
        printf("ERROR: ReadLinearityFile failed with status %d for %s\n",result,input_name);
        return(-1);
      }
    } else {
      input_handle = fopen(input_name,"rt");
      if (input_handle == NULL) {
        printf("ERROR: could not open input file %s\n",input_name);
        return(-1);
      }
      while(1) {
	
        /* Read in the input file and see how many lines it has */
        inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
        if (inBuffer == NULL) {
          break;
        }
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
        if (webuseFlag) {
          if ((strstr(inBuffer,"cat_input") == NULL) &&
              (strstr(inBuffer,"search") == NULL)) {
            continue;
          }
        }
        lineCount++;
      }

      fclose(input_handle);
      input_handle = NULL;
    }
		if (lineCount < pGscBin01->total_gsc_bins) {
			lineCount = pGscBin01->total_gsc_bins;
		}
    if (limitingMagnitudeFlag) {
      allocCount = HISTOGRAM_SIZE+1;
    } else {
      allocCount = lineCount;
    }
    if (linearityFlag) {
      allocCount = pLinearityCommon->maxValue;

    }
    x = (PLFLT*)calloc(allocCount,sizeof(PLFLT));
    y = (PLFLT*)calloc(allocCount,sizeof(PLFLT));
    if (tc720Flag != 0) {
      T1 = (PLFLT*)calloc(allocCount+1,sizeof(PLFLT));
      T2 = (PLFLT*)calloc(allocCount+1,sizeof(PLFLT));
      T3 = (PLFLT*)calloc(allocCount+1,sizeof(PLFLT));
      OUTPUT = (PLFLT*)calloc(allocCount+1,sizeof(PLFLT));
      SETPOINT = (PLFLT*)calloc(allocCount+1,sizeof(PLFLT));
    }



    if ((x == NULL) || (y == NULL)) {
      printf("ERROR: failed to allocate x,y of size %d\n",allocCount);
      exit(-1);
    }

		/* Now allocate our arrays */
    if (limitingMagnitudeFlag) {
      limitingTable = (PLIMITINGTABLE)calloc(lineCount,sizeof(LIMITINGTABLE));
    
      if (limitingTable == NULL){
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
    }
    if ((scangraphFlag != 0) ||
        (scanSmoothedRateFlag != 0)) {
      scangraphTable = (PSCANGRAPH)calloc(lineCount,sizeof(SCANGRAPH)); 

      if (scangraphTable == NULL){
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }

    }
    if (scanrateFlag) {
      scanrateTable = (PSCANRATE)calloc(lineCount,sizeof(SCANRATE)); 

      if (scanrateTable == NULL){
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }

    }
    if (webuseFlag) {
      webuseTable = (PWEBUSE)calloc(lineCount,sizeof(WEBUSE)); 

      if (webuseTable == NULL){
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
      webstatsTable = (PWEBSTATS)calloc(lineCount,sizeof(WEBSTATS)); 

      if (webstatsTable == NULL){
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }



      webhandle = fopen(webfile,"wt");
      if (webhandle == NULL) {
        printf("ERROR: failed to open the web log file %s\n",webfile);
        exit(-1);
      }

    }
    if (tc720Flag) {
      tc720Table = (PTC720)calloc(lineCount,sizeof(TC720)); 
      if (tc720Table == NULL){
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
    }
    if ((lightcurveRMSFlag) || (pointsPerCurveFlag)) {
      lightcurveRMSTable = (PLIGHTCURVERMS)calloc(lineCount,sizeof(LIGHTCURVERMS));
      if (lightcurveRMSTable == NULL) {
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
    }

    if (clipRMSFlag) {
      clipRMSTable = (PCLIPRMS)calloc(lineCount,sizeof(CLIPRMS));
      if (clipRMSTable == NULL) {
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
    }

    if (measPerMagFlag) {
      measPerMagTable = (PMEASPERMAG)calloc(lineCount,sizeof(MEASPERMAG));
      if (measPerMagTable == NULL) {
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
    }

    if (limMagPerYearFlag) {
      limMagPerYearTable = (PLIMMAGPERYEAR)calloc(lineCount,sizeof(LIMMAGPERYEAR));
      if (limMagPerYearTable == NULL) {
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
    }

    if (countPerYearFlag) {
      countPerYearTable = (PCOUNTPERYEAR)calloc(lineCount,sizeof(COUNTPERYEAR));
      if (countPerYearTable == NULL) {
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
    }
    if (cleaningFlag) {
      cleaningTable = (PCLEANING)calloc(lineCount,sizeof(CLEANING));
      if (cleaningTable == NULL) {
        printf("ERROR: failed to allocate tables of size %d\n",lineCount);
        exit(-1);
      }
    }

    /* Now read in the data */
    if (linearityFlag == 0) {
      input_handle = fopen(input_name,"rt");
      if (input_handle == NULL) {
        printf("ERROR: could not open input file %s\n",input_name);
        return(-1);
      }
      while(1) {
	
        /* Read in the input file and see how many lines it has */
        inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
        if (inBuffer == NULL) {
          break;
        }
        lineCount2++;
        if ((lineCount2 <= 2) && (webuseFlag == 0)) {
          continue;
        }
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
        if (limitingMagnitudeFlag) {
          tabPtr = strchr(inBuffer,'\t');
          if (tabPtr != NULL) {
            pLimiting = &limitingTable[curCount];
            *tabPtr = 0;
            tabPtr++;
            strncpy(pLimiting->Plate,inBuffer,MAX_PLATE_STRING);
            pLimiting->Plate[MAX_PLATE_STRING] = 0;
            nvals = sscanf(tabPtr,"%d\t%d\t%d\t%lf",&pLimiting->versionId,&pLimiting->catalogNumber,&pLimiting->patrolPlate,&pLimiting->limiting_mag_local);
            if (nvals == 4) {
              curCount++;
              limitingMagIndex = pLimiting->limiting_mag_local * MAGNITUDE_MULTIPLIER;
              if ((limitingMagIndex > 0) && (limitingMagIndex < HISTOGRAM_SIZE)) {
                pAllPlates->count++;
                pAllPlates->histogram[limitingMagIndex]++;
                if (pAllPlates->count == 1) {
                  pAllPlates->dim_limiting_mag_local = pLimiting->limiting_mag_local;
                  pAllPlates->bright_limiting_mag_local = pLimiting->limiting_mag_local;
                } else {
                  if (pAllPlates->dim_limiting_mag_local < pLimiting->limiting_mag_local) {
                    pAllPlates->dim_limiting_mag_local = pLimiting->limiting_mag_local;
                  }
                  if (pAllPlates->bright_limiting_mag_local > pLimiting->limiting_mag_local) {
                    pAllPlates->bright_limiting_mag_local = pLimiting->limiting_mag_local;
                  }
                }
              
                if (pLimiting->catalogNumber == 0) {
                  if (pLimiting->patrolPlate == 0) {
                    pLimitingStats = pGscNonpatrolPlates;
                  } else if (pLimiting->patrolPlate == 1) {
                    pLimitingStats = pGscPatrolPlates;
                  }
                } else if (pLimiting->catalogNumber == 2) {
                  if (pLimiting->patrolPlate == 0) {
                    pLimitingStats =pApassNonpatrolPlates;
                  } else if (pLimiting->patrolPlate == 1) {
                    pLimitingStats = pApassPatrolPlates; 
                  }
                } else if (pLimiting->catalogNumber == 4) {
                  if (pLimiting->patrolPlate == 0) {
                    pLimitingStats =pAtlasNonpatrolPlates;
                  } else if (pLimiting->patrolPlate == 1) {
                    pLimitingStats = pAtlasPatrolPlates; 
                  }
                } 
                pLimitingStats->count++;
                pLimitingStats->histogram[limitingMagIndex]++;
                if (pLimitingStats->count == 1) {
                  pLimitingStats->dim_limiting_mag_local = pLimiting->limiting_mag_local;
                  pLimitingStats->bright_limiting_mag_local = pLimiting->limiting_mag_local;
                } else {
                  if (pLimitingStats->dim_limiting_mag_local < pLimiting->limiting_mag_local) {
                    pLimitingStats->dim_limiting_mag_local = pLimiting->limiting_mag_local;
                  }
                  if (pLimitingStats->bright_limiting_mag_local > pLimiting->limiting_mag_local) {
                    pLimitingStats->bright_limiting_mag_local = pLimiting->limiting_mag_local;
                  }
                }

              } /* Good limiting magnitude */


            }
          }
        }
        if ((scangraphFlag != 0) ||
            (scanSmoothedRateFlag != 0)) {
          pScangraph = &scangraphTable[curCount];
          nvals = sscanf(inBuffer,"%lf %d",&pScangraph->epoch,&pScangraph->count);
          if (nvals != 2) {
            printf("ERROR in line %d %s\n",curCount+1,inBuffer);
            continue;
          }
          curCount++;
        }
        if (scanrateFlag) {
          pScanrate = &scanrateTable[curCount];
          nvals = sscanf(inBuffer,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",&pScanrate->count,&pScanrate->minimum,&pScanrate->average,&pScanrate->median,&pScanrate->pct600,&pScanrate->pct700,&pScanrate->pct800,&pScanrate->pct900,&pScanrate->pct950,&pScanrate->pct990);
          if (nvals != 10) {
            printf("ERROR in line %d %s\n",curCount+1,inBuffer);
            continue;
          }
          curCount++;
        }
        if ((lightcurveRMSFlag) || (pointsPerCurveFlag)) {
          pLightcurveRMS = &lightcurveRMSTable[curCount];
          nvals = sscanf(inBuffer,"%d\t%d\t%d\t%lf\t%lf\t%d\t%d",
                         &pLightcurveRMS->minmag,
                         &pLightcurveRMS->maxmag,
                         &pLightcurveRMS->count,
                         &pLightcurveRMS->median,
                         &pLightcurveRMS->normnburst3,
                         &pLightcurveRMS->catalogNumber,
                         &pLightcurveRMS->npoints);
          if (nvals != 7) {
            printf("WARNING in line %d %s\n",curCount+1,inBuffer);
            continue;
          }
          curCount++;
          if (pLightcurveRMS->minmag < minmag) {
            minmag = pLightcurveRMS->minmag;
          }
          if (pLightcurveRMS->maxmag > maxmag) {
            maxmag = pLightcurveRMS->maxmag;
          }      
          if ((pointsPerCurveFlag != 0) && (pLightcurveRMS->minmag > MIN_POINTS_PER_CURVE_MAG) && (pLightcurveRMS->npoints > maxnpoints)) {
            maxnpoints = pLightcurveRMS->npoints;
          }
          if ((minmedian == 0.0) || (pLightcurveRMS->median < minmedian)) {
            minmedian = pLightcurveRMS->median;
          }
          if ((maxmedian == 0.0) || (pLightcurveRMS->median > maxmedian)) {
            maxmedian = pLightcurveRMS->median;
          }

        
        }
        if (clipRMSFlag) {
          pClipRMS = &clipRMSTable[curCount];
          nvals = sscanf(inBuffer,"%d\t%d\t%lf\t%d\t%d",
                         &pClipRMS->minmag,
                         &pClipRMS->maxmag,
                         &pClipRMS->clip_rms_local,
                         &pClipRMS->count,
                         &pClipRMS->catalogNumber);
          if (nvals != 5) {
            printf("ERROR in line %d %s\n",curCount+1,inBuffer);
            continue;
          }
          curCount++;
          if (pClipRMS->minmag < minmag) {
            minmag = pClipRMS->minmag;
          }
          if (pClipRMS->maxmag > maxmag) {
            maxmag = pClipRMS->maxmag;
          }
          if (pClipRMS->count > maxCount) {
            maxCount = pClipRMS->count;
          }
        
        }

        if (measPerMagFlag) {
          pMeasPerMag = &measPerMagTable[curCount];
          nvals = sscanf(inBuffer,"%d\t%d\t%d\t%lld",
                         &pMeasPerMag->sequence,
                         &pMeasPerMag->magcal_magdep,
                         &pMeasPerMag->catalogNumber,
                         &pMeasPerMag->count);
          if (nvals != 4) {
            printf("WARNING in line %d %s\n",curCount+1,inBuffer);
            continue;
          }
          curCount++;
          if ((pMeasPerMag->sequence == 1) || (pMeasPerMag->sequence == 2)) {
            if (pMeasPerMag->magcal_magdep < minmag) {
              minmag = pMeasPerMag->magcal_magdep;
            }
            if (pMeasPerMag->magcal_magdep > maxmag) {
              maxmag = pMeasPerMag->magcal_magdep;
            }
            if (pMeasPerMag->count > maxCount) {
              maxCount = pMeasPerMag->count;
            }
          }
        }

        if (limMagPerYearFlag) {
          pLimMagPerYear = &limMagPerYearTable[curCount];
          nvals = sscanf(inBuffer,"%d\t%d\t%d\t%lf",
                         &pLimMagPerYear->sequence,
                         &pLimMagPerYear->year,
                         &pLimMagPerYear->catalogNumber,
                         &pLimMagPerYear->limiting_mag_local);
          if (nvals != 4) {
            printf("WARNING in line %d %s\n",curCount+1,inBuffer);
            continue;
          }
          curCount++;
          if ((pLimMagPerYear->sequence == 1) || (pLimMagPerYear->sequence == 2)) {
            if (pLimMagPerYear->year < minyear) {
              minyear = pLimMagPerYear->year;
            }
            if (pLimMagPerYear->year > maxyear) {
              maxyear = pLimMagPerYear->year;
            }

            if (pLimMagPerYear->limiting_mag_local < minlimitingmag) {
              minlimitingmag = pLimMagPerYear->limiting_mag_local;
            }
            if (pLimMagPerYear->limiting_mag_local > maxlimitingmag) {
              maxlimitingmag = pLimMagPerYear->limiting_mag_local;
            }



          }
        }
        if (cleaningFlag) {
          pCleaning = &cleaningTable[maxCleaningCount];
          memset(pCleaning,0,sizeof(CLEANING));
          if (strstr(inBuffer,"Machine Cleaned") != 0) {
            pCleaning->sequence = MACHINE_CLEAN;
            machineCleanCount++;
            pCleaning->count = machineCleanCount;
          } else if (strstr(inBuffer,"Hand Cleaned") != 0) {
            pCleaning->sequence = HAND_CLEAN;
            handCleanCount++;
            pCleaning->count = handCleanCount;
          } else {
            continue;
          }
          charPtr3 = strchr(inBuffer,'\t');
          if (charPtr3 != NULL) {
            *charPtr3 = 0;
          } else {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          if (strlen(inBuffer) >= (MAX_DATE_STRING-1)) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;

          }
          strcpy(pCleaning->date,inBuffer);
          charPtr3 = strchr(inBuffer,'T');
          if (charPtr3 != NULL) {
            *charPtr3 = 0;
          } else {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          strcpy(pCleaning->day,inBuffer);
          julianDate = fd2jd(pCleaning->date);
          pCleaning->epoch = jd2ep(julianDate);
          if (pCleaning->epoch < minHandCleanYear) {
            minHandCleanYear = pCleaning->epoch;
          }
          if (pCleaning->epoch > maxHandCleanYear) {
            maxHandCleanYear = pCleaning->epoch;
          }
          maxCleaningCount++;
        }


        if (countPerYearFlag) {
          pCountPerYear = &countPerYearTable[curCount];
          nvals = sscanf(inBuffer,"%d\t%d\t%d\t%lld",
                         &pCountPerYear->sequence,
                         &pCountPerYear->year,
                         &pCountPerYear->catalogNumber,
                         &pCountPerYear->count);
          if (nvals != 4) {
            printf("WARNING in line %d %s\n",curCount+1,inBuffer);
            continue;
          }
          curCount++;
          if ((pCountPerYear->sequence == 1) || (pCountPerYear->sequence == 2)) {
            if (pCountPerYear->year < minyear) {
              minyear = pCountPerYear->year;
            }
            if (pCountPerYear->year > maxyear) {
              maxyear = pCountPerYear->year;
            }

            if (pCountPerYear->count > maxcountperyear) {
              maxcountperyear = pCountPerYear->count;
            }



          }
        }

        if (webuseFlag) {
          pWebUse = &webuseTable[curCount];
          memset(pWebUse,0,sizeof(WEBUSE));
          pWebUse->sequence = curCount;
          if (strstr(inBuffer,"elos") != NULL) {
            continue;
          }
          if (strstr(inBuffer,"search") != NULL) {
            pWebUse->webtype = WEB_TYPE_PLATE_SEARCH;
          } else if (strstr(inBuffer,"cat_input") != NULL) {
            pWebUse->webtype = WEB_TYPE_OBJECT_SEARCH;
          } else {
            continue;
          }
          strcpy(copyLine,inBuffer);
          charPtr1 = strtok(copyLine," ");
          if (charPtr1 == NULL) {
            continue;
          }
          if (strstr(charPtr1,"search") != NULL) {
            pWebUse->webtype = WEB_TYPE_PLATE_SEARCH;
          } else if (strstr(charPtr1,"cat_input") != NULL) {
            pWebUse->webtype = WEB_TYPE_OBJECT_SEARCH;
          } else {
            printf("ERROR no type in line %s\n",inBuffer);
            continue;
          }
          charPtr1 = strtok(NULL," ");
          if (charPtr1 == NULL) {
            continue;
          }
          if (strlen(charPtr1) > MAX_DATE_STRING) {
            printf("ERROR no date in line %s\n",inBuffer);
            continue;
          }
          strcpy(pWebUse->date,charPtr1);

          charPtr1 = strtok(NULL," ");
          if (charPtr1 == NULL) {
            continue;
          }
          if (strlen(charPtr1) > MAX_DATE_STRING) {
            printf("ERROR no time in line %s\n",inBuffer);
            continue;
          }
          strcpy(pWebUse->time,charPtr1);

          charPtr1 = strtok(NULL," ");
          if (charPtr1 == NULL) {
            continue;
          }
          if (strlen(charPtr1) > MAX_URL_STRING) {
            printf("ERROR no URL in line %s\n",inBuffer);
            continue;
          }
          strcpy(pWebUse->url,charPtr1);
          objectLength = 0;
          while (1) {
            charPtr1 = strtok(NULL," ");
            if (charPtr1 == NULL) {
              break;
            }
            if ((strcmp(charPtr1,"wcsfit") == 0) || 
                (strcmp(charPtr1,"R:") == 0)) {
              break;
            }
            objectLength = strlen(pWebUse->object);
            if ((objectLength + strlen(charPtr1) + 3) > MAX_OBJECT_LENGTH) {
              break;
            }
            strcat(pWebUse->object,charPtr1);
            strcat(pWebUse->object," ");
          }
          objectLength = strlen(pWebUse->object);
          if (pWebUse->object[objectLength-1] == ' ') {
            pWebUse->object[objectLength-1] = 0;
          }

          /* Now find the ephemeris date */
          if ((strlen(pWebUse->date) + strlen(pWebUse->time) + 4) >= MAX_DATE_STRING) {
            printf("ERROR full date too long in line %s\n",inBuffer);
            continue;
          }
          strcpy(copyDate,pWebUse->date);
          strcpy(fullDate,"20");
          strcat(fullDate,&copyDate[6]);
          copyDate[2] = 0;
          strcat(fullDate,"-");
          strcat(fullDate,copyDate);
          copyDate[5] = 0;
          strcat(fullDate,"-");
          strcat(fullDate,&copyDate[3]);
          strcat(fullDate,"T");
          strcat(fullDate,pWebUse->time);
          julianDate = fd2jd(fullDate);
          pWebUse->epoch = jd2ep(julianDate);
          if (pWebUse->epoch < RELEASE_DATE) {
            if (curCount > 0) {
              printf("ERROR pre-release date %.2f in line %s\n",pWebUse->epoch,inBuffer);
            }
            continue;
          }
          if (minEpoch == 0.0) {
            minEpoch = pWebUse->epoch;
            maxEpoch = pWebUse->epoch;
          } else {
            if (minEpoch > pWebUse->epoch) {
              minEpoch = pWebUse->epoch;
            }
            if (maxEpoch < pWebUse->epoch) {
              maxEpoch = pWebUse->epoch;
            }
          }
          pWebUse->month = 12.0*(pWebUse->epoch-BASE_DATE);
          if (pWebUse->month >= MAX_MONTH) {
            printf("ERROR date too high %.2f in line %s\n",pWebUse->epoch,inBuffer);          
            continue;
          }
        
          curCount++;
#if 1
          fprintf(webhandle,"%s %s %d\n",fullDate,pWebUse->url,pWebUse->webtype);
#endif
          if (pWebUse->webtype == WEB_TYPE_PLATE_SEARCH) {
            plateSearchCount++;
          } else {
            lightcurveSearchCount++;
          }

        }
        if (tc720Flag) {
          ptc720 = &tc720Table[curCount];
          memset(ptc720,0,sizeof(TC720));
          strcpy(copyLine,inBuffer);

          /* Day of week */
          charPtr1 = strtok(copyLine," ");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          for (ptc720->tmstruct.tm_wday = 0; ptc720->tmstruct.tm_wday < numWeekdays; ptc720->tmstruct.tm_wday++) {
            charPtr2 = strstr(charPtr1,weekdays[ptc720->tmstruct.tm_wday]);
            if (charPtr2 != NULL) {
              break;
            }
          }
          if (charPtr2 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }

          /* Month */
          charPtr1 = strtok(NULL," ");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          for (ptc720->tmstruct.tm_mon = 0; ptc720->tmstruct.tm_mon < numMonthstr; ptc720->tmstruct.tm_mon++) {
            charPtr2 = strstr(charPtr1,monthstr[ptc720->tmstruct.tm_mon]);
            if (charPtr2 != NULL) {
              break;
            }
          }
          if (charPtr2 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }

          /* Day of month */
          charPtr1 = strtok(NULL," ");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%d",&ptc720->tmstruct.tm_mday);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }

          /* Year */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%d",&ptc720->tmstruct.tm_year);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          ptc720->tmstruct.tm_year -= 1900;

          /* Hour */
          charPtr1 = strtok(NULL," :\t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%d",&ptc720->tmstruct.tm_hour);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }

          /* minute */
          charPtr1 = strtok(NULL," :\t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%d",&ptc720->tmstruct.tm_min);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }

          /* second */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%d",&ptc720->tmstruct.tm_sec);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          

          /* AM or PM */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }

          if ((strstr(charPtr1,"PM") != NULL)) {
            if (ptc720->tmstruct.tm_hour < 12) {
              ptc720->tmstruct.tm_hour += 12;
            }
          } else {
            if (strstr(charPtr1,"AM") == NULL) {
              printf("ERROR line %d for %s\n",__LINE__,inBuffer);
              continue;
            }
          }
          ptc720->tmstruct.tm_isdst = -1;
#if 0
          printf("ERROR: PLOT_DAY base time test\n");
          ptc720->tmstruct.tm_sec = 0;
          ptc720->tmstruct.tm_min = 0;
          ptc720->tmstruct.tm_hour = 0;
          ptc720->tmstruct.tm_mon = 2;
          ptc720->tmstruct.tm_mday = 1;
          ptc720->tmstruct.tm_year = 119;
          ptc720->tmstruct.tm_isdst = -1;

          ptc720->tmstruct.tm_wday = 0;
          ptc720->tmstruct.tm_yday = 0;
          
#endif

          ptc720->tp = mktime(&ptc720->tmstruct);
          ptr = localtime(&ptc720->tp);
          strftime(ptc720->timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
          ptc720->hour = (1.0*ptc720->tmstruct.tm_hour)+ (ptc720->tmstruct.tm_min/60.0) + (ptc720->tmstruct.tm_sec/3600.);

          

          /* unknown column */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          if (strstr(charPtr1,"0") == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          
          /* T1 */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%lf",&ptc720->T1);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          

          /* T2 */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%lf",&ptc720->T2);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }

          /* setpoint */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%lf",&ptc720->setpoint);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }


          /* output */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%lf",&ptc720->output);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }

          /* T3 */
          charPtr1 = strtok(NULL," \t");
          if (charPtr1 == NULL) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
          nvals = sscanf(charPtr1,"%lf",&ptc720->T3);
          if (nvals != 1) {
            printf("ERROR line %d for %s\n",__LINE__,inBuffer);
            continue;
          }
#ifndef USE_CENTIGRADE
          /* convert to Farenheight */
          ptc720->T1 = ((9.0*ptc720->T1)/5.0)+ 32.;
          ptc720->T2 = ((9.0*ptc720->T2)/5.0)+ 32.;
          ptc720->T3 = ((9.0*ptc720->T3)/5.0)+ 32.;
          ptc720->setpoint = ((9.0*ptc720->setpoint)/5.0)+ 32.;
#endif /* USE_CENTIGRADE */

          if (curCount == 0) {
            memcpy(&mintc720,ptc720,sizeof(TC720));
            memcpy(&maxtc720,ptc720,sizeof(TC720));
          } else {
            if (ptc720->tp < mintc720.tp) {
              mintc720.tp = ptc720->tp;
            }
            if (ptc720->T1 < mintc720.T1) {
              mintc720.T1 = ptc720->T1;
            }
            if (ptc720->T2 < mintc720.T2) {
              mintc720.T2 = ptc720->T2;
            }
            if (ptc720->setpoint < mintc720.setpoint) {
              mintc720.setpoint = ptc720->setpoint;
            }
            if (ptc720->output < mintc720.output) {
              mintc720.output = ptc720->output;
            }
            if (ptc720->T3 < mintc720.T3) {
              mintc720.T3 = ptc720->T3;
            }
            if (ptc720->hour < mintc720.hour) {
              mintc720.hour = ptc720->hour;
            }


            if (ptc720->tp > maxtc720.tp) {
              maxtc720.tp = ptc720->tp;
            }
            if (ptc720->T1 > maxtc720.T1) {
              maxtc720.T1 = ptc720->T1;
            }
            if (ptc720->T2 > maxtc720.T2) {
              maxtc720.T2 = ptc720->T2;
            }
            if (ptc720->setpoint > maxtc720.setpoint) {
              maxtc720.setpoint = ptc720->setpoint;
            }
            if (ptc720->output > maxtc720.output) {
              maxtc720.output = ptc720->output;
            }
            if (ptc720->T3 > maxtc720.T3) {
              maxtc720.T3 = ptc720->T3;
            }
            if (ptc720->hour > maxtc720.hour) {
              maxtc720.hour = ptc720->hour;
            }

          }

#if 0
          printf("line %d %s %s %lld\n",__LINE__,inBuffer,ptc720->timestr,ptc720->tp);

#endif

          curCount++;


        }

      }
      printf("Records read %d, goodRecords %d\n",lineCount2,curCount);
          
      if (tc720Flag) {
        printf("T1 range %4.1f to %4.1f  T2 range %4.1f to %4.1f  T3 range %4.1f to %4.1f  setpoint range %4.1f to %4.1f  output range %4.1f to %4.1f hour range %5.2f to %5.2f seconds %lld hours %f\n",
               mintc720.T1,maxtc720.T1,  
               mintc720.T2,maxtc720.T2,  
               mintc720.T3,maxtc720.T3,  
               mintc720.setpoint,maxtc720.setpoint,  
               mintc720.output,maxtc720.output,  
               mintc720.hour,maxtc720.hour,
               maxtc720.tp-mintc720.tp,
               ((1.0*(maxtc720.tp-mintc720.tp))/3600.));
      }


      if (webuseFlag != 0) {
        /* Eliminate duplicate object requests from the same URL */
        qsort((void*)webuseTable,curCount,sizeof(WEBUSE),WebStatsDuplicate);
        webuseCount = 0;
        pWebUse2 = &webuseTable[webuseCount];
        for (pointIndex = 1; pointIndex < curCount; pointIndex++) {
          pWebUse = &webuseTable[pointIndex];
          if ((pWebUse->webtype == pWebUse2->webtype) &&
              (strcmp(pWebUse->url,pWebUse2->url) == 0) &&
              (strcmp(pWebUse->object,pWebUse2->object) == 0)) {
            if ((pWebUse->epoch-pWebUse2->epoch) < MIN_EPOCH) {
              continue;
            }
          }
          webuseCount++;
          pWebUse2 = &webuseTable[webuseCount];
          if (pointIndex != webuseCount) {
            memcpy(pWebUse2,pWebUse,sizeof(WEBUSE));
          }
#if 0
          printf("OBJECT: %s\n",pWebUse->object);
#endif
        }
        printf("Initial count %d after duplicates %d difference %d\n",curCount,webuseCount,curCount-webuseCount);
        curCount = webuseCount;
        qsort((void*)webuseTable,curCount,sizeof(WEBUSE),WebStatsSequence);
#if 0
        for (pointIndex = 1; pointIndex < curCount; pointIndex++) {
          pWebUse = &webuseTable[pointIndex];
          printf("pointIndex %5d sequence %5d type %d url %20s object %80s epoch %f\n",
                 pointIndex,
                 pWebUse->sequence,
                 pWebUse->webtype,
                 pWebUse->url,
                 pWebUse->object,
                 pWebUse->epoch);
        }
#endif


        for (pointIndex = 1; pointIndex < curCount; pointIndex++) {
          pWebUse = &webuseTable[pointIndex];
          monthArray[pWebUse->month]++;
          /* .cfa.harvard.edu terminates in .142.131
           *  dhcp-128-103-233-245.harvard.edu is 245.233.103.128
           *  .wrls.harvard.edu is .247.140 and .134.67
           *  .fas.harvard.edu is .247.140
           * cncdnh.east.myfairpoint.net terminates in .109.70
           * vpn-10-11-3-133.vpnclient.harvard.edu 49.149.101.46 
           * 
           */
#if 1
          if ((strstr(pWebUse->url,"77.88.5.121") != NULL)) {
            printf("line %d at %s\n",__LINE__,pWebUse->url);
          }
#endif
          acceptedFlag = 0;
          if ((strstr(pWebUse->url,"131.142.") == NULL) &&
              (strstr(pWebUse->url,"128.103.233.245") == NULL) &&
              (strstr(pWebUse->url,"140.247.") == NULL) &&
              (strstr(pWebUse->url,"67.134.") == NULL) &&
              (strstr(pWebUse->url,"70.109.") == NULL) &&
              (strstr(pWebUse->url,"49.149.101.46") == NULL) && /* vpn-10-11-3-133.vpnclient.harvard.edu */
              (strstr(pWebUse->url,"107.23.37.51") == NULL) && /* ec2-107-23-37-51.compute-1.amazonaws.com 2372 in 2015-07 */
              (strstr(pWebUse->url,"78.163.106.223") == NULL) && /* Attack of Dec 9, 2017 */
              (strstr(pWebUse->url,"149.154.65.151") == NULL) && /* Performance test of  2017-09-19 or denial of service attack? */

              (strstr(pWebUse->url,"146.6.15.11") == NULL) && /* Attack of 05/07/18 */
              (strstr(pWebUse->url,"52.56.185.44") == NULL) && /* Attack of 04/28/18 */

#if 0
              (strstr(pWebUse->url,"153.9.254.113") == NULL) && /* Many retries of same object */
              (strstr(pWebUse->url,"71.226.97.61") == NULL) && /* Many retries of same object */
              (strstr(pWebUse->url,"146.6.15.11") == NULL) && /* Many retries of same object garbage retries */
              (strstr(pWebUse->url,"45.32.100.122") == NULL) && /* Many retries of same object */
              (strstr(pWebUse->url,"37.17.224.72") == NULL) && /* Many retries of same object */
#endif         
#ifdef LUNDS_UNIVERSITY /* 75774 Legitimate */
              (strstr(pWebUse->url,"130.235.102.60") == NULL) && /* Legitimate, from  pc140.astro.lu.se 51,335 automated requests */
              (strstr(pWebUse->url,"84.217.57.67") == NULL) &&   /* Legitimate, from  84-217-57-67.customers.ownit.se: 19,869 automated requests */
              (strstr(pWebUse->url,"130.235.136.42") == NULL) &&   /* Legitimate, from eduroam42.nat.wireless.lu.se 3,777 automated requests */
              (strstr(pWebUse->url,"130.235.240.10") == NULL) &&   /* Legitimate, from eduroam-employee-10.wireless.lu.se  793 automated requests */
#endif
              
              (strstr(pWebUse->url,"10.11.3.131") == NULL) &&   /* White hat testing 3/25/20 ? */
              (strstr(pWebUse->url,"10.11.3.140") == NULL) &&   /* White hat testing 3/25/20 ? */
              (strstr(pWebUse->url,"10.11.3.133") == NULL) &&   /* White hat testing 3/25/20 ? */

              (strstr(pWebUse->url,"10.11.3.131") == NULL) &&   /* White hat testing 3/18/20 */
              (strstr(pWebUse->url,"10.11.3.140") == NULL) &&   /* White hat testing 4/09/20 */
              (strstr(pWebUse->url,"10.11.3.133") == NULL) &&   /* White hat testing 3/25/20 ? */
              (strstr(pWebUse->url,"71.168.78.209") == NULL) &&   /* cncdnh.east.myfairpoint.net */
              (strstr(pWebUse->url,"78.163.106.223") == NULL) &&   /* White hat testing 12/09/17 ? */
              (strstr(pWebUse->url,"116.193.140.123") == NULL) &&   /* White hat testing 3/25/20 ? */
              (strstr(pWebUse->url,"185.244.215.9") == NULL) &&   /* White hat testing 3/16/20 and 3/17/20 no-mans-land.m247.com max_hyppolite@harvard.edu  */
              (strstr(pWebUse->url,"178.175.132.72") == NULL) &&   /* apparent White hat testing 07/05/20  */
              (strstr(pWebUse->url,"209.97.153.75") == NULL) &&   /* apparent White hat testing 09/13/20  */
              

              (strstr(pWebUse->url,"185.173.205.143") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"185.173.205.137") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"52.247.192.191") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"104.237.156.251") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"73.60.100.95") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"10.11.3.22") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"95.217.250.108") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"118.71.137.183") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"123.16.48.171") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"51.223.78.26") == NULL) && /* Found in 2021_01_11 */
              (strstr(pWebUse->url,"159.0.78.80") == NULL) && /* Found in 2021_01_11 */


              /* White Hat security scan (10/29/2015 3:26) */
              (strstr(pWebUse->url,"63.128.163.") == NULL) &&
              (strstr(pWebUse->url,"63.128.163.8") == NULL) &&
              (strstr(pWebUse->url,"12.248.108.") == NULL) &&
              (strstr(pWebUse->url,"67.207.113.") == NULL) &&
              (strstr(pWebUse->url,"64.244.165.") == NULL) &&
              (strstr(pWebUse->url,"38.122.74.16") == NULL)
              ) {
            acceptedFlag = 1;
            if ((strcmp(pWebUse->url,"151.188.17.247") == 0) || /* Dump of 1269 KIC objects */
                (strcmp(pWebUse->url,"128.118.147.210") == 0) || /* periodic check of TU UMa */
                (strcmp(pWebUse->url,"210.72.91.90") == 0) || /* legitimate use */
                (strcmp(pWebUse->url,"98.209.10.84") == 0) || /* Peter.Winter c-98-209-10-84.hsd1.mi.comcast.net. does show possible javascript attacks */
                
                (strcmp(pWebUse->url,"210.72.91.69") == 0))  /* March 9, legitimate use  */
              { 
                pWebStats->goodFlag = 1;
              }
            nonharvardArray[pWebUse->month]++;
            if (nonharvardMaxMonth < nonharvardArray[pWebUse->month]) {
              nonharvardMaxMonth = nonharvardArray[pWebUse->month];
            }
            nonharvardCount++;
          }
          pWebStats = &webstatsTable[pointIndex];
          memset(pWebStats,0,sizeof(WEBSTATS));
          strcpy(pWebStats->url,pWebUse->url);
          pWebStats->acceptedFlag = acceptedFlag;
          pWebStats->totalCount = 1;
        }
        /* Eliminate duplicates */
        qsort((void*)webstatsTable,curCount,sizeof(WEBSTATS),WebStatsCompare);
        nondupCount = 0;
        pWebStats2 = &webstatsTable[nondupCount];
        for (pointIndex = 1; pointIndex < curCount; pointIndex++) {
          pWebStats = &webstatsTable[pointIndex];
          if (strcmp(pWebStats->url,pWebStats2->url) == 0) {
            pWebStats2->totalCount++;
            continue;
          }
          nondupCount++;
          pWebStats2 = &webstatsTable[nondupCount];
          if (pointIndex != nondupCount) {
            memcpy(pWebStats2,pWebStats,sizeof(WEBSTATS));
          }
        }
        /* Now sort by totalCount */
        qsort((void*)webstatsTable,nondupCount,sizeof(WEBSTATS),WebStatsTotal);
        for (pointIndex = 1; pointIndex < nondupCount; pointIndex++) {
          pWebStats = &webstatsTable[pointIndex];

#if 0
          if ((strstr(pWebStats->url,"84.217.57.67") != NULL)) {
            printf("line %d at %s\n",__LINE__,pWebUse->url);
          }
#endif

          fprintf(webhandle,"pointIndex %5d address %20s count %5d accepted %d good %d\n",pointIndex,pWebStats->url,pWebStats->totalCount,pWebStats->acceptedFlag,pWebStats->goodFlag);
        }


        printf("Minimum date: %.3f Maximum date: %.3f plate searches: %d lightcurve searches %d non-Harvard searches %d\n",minEpoch,maxEpoch,plateSearchCount,lightcurveSearchCount,nonharvardCount);
        month = 1;
        year = BASE_DATE;
        for (monthIndex = 0; monthIndex < MAX_MONTH; monthIndex++) {
          if (monthArray[monthIndex] != 0) {
#if 1
            printf("Year: %4d  Month %2d Searches: %5d non-Harvard: %5d\n",year,month,monthArray[monthIndex],nonharvardArray[monthIndex]);
#endif
            maxMonthIndex = monthIndex;
          }
          month++;
          if (month >= 13) {
            month = 1;
            year++;
          }
        }
      }
    
      if (limitingMagnitudeFlag) {
        for (statsIndex = 0; statsIndex < labelSize; statsIndex++) {
          pLabelEntry = &statsLabels[statsIndex];
          pLimitingStats = pLabelEntry->statsTable;
          printf("Index %d  Type %25s  Plates %6d  limiting_mag_local %5.2f to %5.2f\n",statsIndex,pLabelEntry->plotLabel,pLimitingStats->count,pLimitingStats->bright_limiting_mag_local,pLimitingStats->dim_limiting_mag_local);
          limitingMagCount = 0;
          /* Now convert the histogram to a cumulative histogram with 100% at the brightest magnitude */
          for (limitingMagIndex = HISTOGRAM_SIZE -1; limitingMagIndex >= 0; limitingMagIndex--) {
            limitingMagCount += pLimitingStats->histogram[limitingMagIndex];
            pLimitingStats->histogram[limitingMagIndex] = limitingMagCount;
          }
          if (limitingMagCount != pLimitingStats->count) {
            printf("ERROR: limitingMagCount %d is not count %d\n",limitingMagCount,pLimitingStats->count);
          }
        }
      }
      fclose(input_handle);
      input_handle = NULL;
    }
  }



  if (limitingTable != NULL) {
    free(limitingTable);
  }

  if (limitingMagnitudeFlag) {
    ymin = 0.0;
    ymax = 100.0;

    for (statsIndex = 0; statsIndex < labelSize; statsIndex++) {
      pLabelEntry = &statsLabels[statsIndex];
      pLimitingStats = pLabelEntry->statsTable;
      if (pLimitingStats->count > maxcount) {
        maxcount = pLimitingStats->count;
      }
      if (statsIndex == 0) {
        xmin = pLimitingStats->bright_limiting_mag_local;
        xmax = pLimitingStats->dim_limiting_mag_local;
      } else {
        if (xmin > pLimitingStats->bright_limiting_mag_local) {
          xmin = pLimitingStats->bright_limiting_mag_local;
        }
        if (xmax < pLimitingStats->dim_limiting_mag_local) {
          xmax = pLimitingStats->dim_limiting_mag_local;
        }
      }

    }
    xmin = 10.0;
    xmax = 22.0;
  }
  if (scangraphFlag) {
    ymin = 0.0;
    ymax = 0.0;
    for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
      pScangraph = &scangraphTable[pointIndex];
      if (pointIndex == 0) {
        xmax = pScangraph->epoch;
        xmin = pScangraph->epoch;
      } else {
        if (pScangraph->epoch > xmax) {
          xmax = pScangraph->epoch;
        }
        if (pScangraph->count > ymax) {
          ymax = pScangraph->count;
        }
      }
    }
    xmax += 0.1;
    ymax += 1000;
#ifdef BASE_YEAR
    xmax -= BASE_YEAR;
    xmin -= BASE_YEAR;
#endif /* BASE_YEAR */
  }
  if (scanSmoothedRateFlag) {
    ymin = 0.0;
    ymax = 0.0;
    pointIndex1 = 0; /* earlier date */
    pointIndex2 = 0; /* later date */
    while (pointIndex2 < curCount) {
      pScangraph1 = &scangraphTable[pointIndex1];
      pScangraph2 = &scangraphTable[pointIndex2];
      currentDate = pScangraph1->epoch;
      laterDate = pScangraph2->epoch;
      if ((laterDate-currentDate) < (smoothingDays*(YEARS_PER_DAY))) {
        pointIndex2++;
      } else {
        pScangraph1->rate = ((1.0*(pointIndex2-pointIndex1))/(laterDate-currentDate))*YEARS_PER_WEEK;
        pScangraph1->midepoch = (currentDate+((laterDate-currentDate)/2.0));
        if (pointIndex1 == 0) {
          xmin = currentDate;
          ymax = pScangraph1->rate;
          ymin = pScangraph1->rate;
          peakCurrentDate = currentDate;
          peakLaterDate = laterDate;
        } else {
          if (pScangraph1->rate > ymax) {
            ymax = pScangraph1->rate;
            peakCurrentDate = currentDate;
            peakLaterDate = laterDate;
          }
          if (pScangraph1->rate < ymin) {
            ymin = pScangraph1->rate;
          }
        }
#if 0
        printf("scanSmoothedRate pointIndex1 %d currentDate %f pointIndex2 laterDate %f midepoch %f rate %f %d xmin %f xmax %f ymin %f ymax %f\n",pointIndex1,currentDate,pointIndex2,laterDate,pScangraph1->midepoch,pScangraph1->rate,xmin,xmax,ymin,ymax);
#endif
        pointIndex1++;
      }
    }
    curCount = pointIndex1;
    xmax = currentDate;
    printf("scanSmoothedRate ymax %f plates per week peakCurrentDate %f peakLaterDate %f\n",ymax,peakCurrentDate,peakLaterDate);
#if 0
    printf("scanSmoothedRate xmin %f xmax %f ymin %f ymax %f\n",xmin,xmax,ymin,ymax);
#endif
    xmax += 0.1;
    ymax += 0.1;
#ifdef BASE_YEAR
    xmax -= BASE_YEAR;
    xmin  = 0;
#endif /* BASE_YEAR */
  }
  if (scanrateFlag) {
#if 1 /* Normal */
    xmin = 0;
    xmax = 510;
    ymin = 30;
    ymax = 70;
#endif
#if 0 /* Large Plates */
    xmin = 0;
    xmax = 100;
    ymin = 15;
    ymax = 22;

#endif 
#if 0 /* Small Plates */
    xmin = 0;
    xmax = 510;
    ymin = 50;
    ymax = 130;
#endif
  }
  if (lightcurveRMSFlag) {
    xmin = 1.0*minmag;
    xmax = 1.0*maxmag;
#if 1
    ymin = 0.0;
#else
    ymin = minmedian;
#endif
    ymax = maxmedian+0.025;
  }
  if (pointsPerCurveFlag) {
    xmin = MIN_POINTS_PER_CURVE_MAG;
    xmax = 1.0*maxmag;
    ymin = 0.0;
    ymax = 1.1*maxnpoints;
  }
  if (clipRMSFlag) {
    xmin = 0.0;
    xmax = 0.3;
    ymin = 0.0;
    ymax = 1.025*maxCount;
  }
  if (measPerMagFlag) {
    xmin = 8.0;
    xmax = 1.0*maxmag;
    ymin = 0.0;
    ymax = 1.025*maxCount;
  }
  if (limMagPerYearFlag) {
    xmin = 1.0*(minyear-1);
    xmax = 1.0*(maxyear+1);
    ymax = 0.975*minlimitingmag;
    ymin = 1.025*maxlimitingmag;
  }

  if (countPerYearFlag) {
    xmin = 1.0*(minyear-1);
    xmax = 1990.;
    ymin = 0;
    ymax = 1.025*maxcountperyear;
  }

  if (cleaningFlag) {
    xmin = minHandCleanYear;
    xmax = maxHandCleanYear;
    ymin = 0.0;
    ymax = 0.0;
#ifdef  CLEAN_IDLE_INTERVAL
    for (pointIndex = 0; pointIndex < maxCleaningCount; pointIndex++) {
      pCleaning = &cleaningTable[pointIndex];
      if (pCleaning->sequence == HAND_CLEAN) {
        if (startHandCleanYear < 0) {
          startHandCleanYear = pCleaning->epoch;
          lastHandCleanYear = pCleaning->epoch;
          tempHandCleanCount = 1;
          pLastHandCleaning = pCleaning;
          
        } else {
          if ((pCleaning->epoch - lastHandCleanYear) < 0) {
            printf("ERROR: line %d year not in order at %s\n",__LINE__,pCleaning->date);
          }

          if ((pCleaning->epoch - lastHandCleanYear) < CLEAN_IDLE_INTERVAL) {
            tempHandCleanCount++;
            lastHandCleanYear = pCleaning->epoch;
            pLastHandCleaning = pCleaning;
          } else {
            if ((pLastHandCleaning->epoch-startHandCleanYear) <= 0) {
              printf("ERROR: line %d zero interval %s\n",__LINE__,pLastHandCleaning->date);
            } else {
#if 1
              if (tempHandCleanCount == 345) {
                printf("line %d tempHandCleanCount %d\n",__LINE__,tempHandCleanCount);
              }
#endif
              if (tempHandCleanCount >= MIN_SESSION_COUNT) {
                pLastHandCleaning->sessionCount = tempHandCleanCount;
                pLastHandCleaning->hours = (pLastHandCleaning->epoch-startHandCleanYear)*365.25*24.0;
                pLastHandCleaning->rate = (1.0*pLastHandCleaning->sessionCount)/pLastHandCleaning->hours;
                if (pLastHandCleaning->rate < CLEANING_MAX_RATE) {
                  totalHandCleanCount += pLastHandCleaning->sessionCount;
                  totalHandCleanHours += pLastHandCleaning->hours;
#if 1
                  printf("line %d sessionCount %d hours %f rate %f epoch %f ymax %f machine clean\n",__LINE__,pLastHandCleaning->sessionCount,pLastHandCleaning->hours,pLastHandCleaning->rate,pLastHandCleaning->epoch,ymax);
#endif
                  if (pLastHandCleaning->rate > ymax) {
                    ymax = pLastHandCleaning->rate;
                  }
                } else { 
                  pLastHandCleaning->sessionCount = 0;
                  pLastHandCleaning->hours = 0;
                  pLastHandCleaning->rate = 0;

                }
              }
            }
            startHandCleanYear = pCleaning->epoch;
            lastHandCleanYear = pCleaning->epoch;
            tempHandCleanCount = 1;
            pLastHandCleaning = pCleaning;
          }
        }
      } else {
#if 0
        if (pLastMachineCleaning == NULL) {
          printf("tempMachineCleanCount %d epoch %f startMachineCleanYear %f, lastMachineCleanYear %f\n",tempMachineCleanCount,pCleaning->epoch,startMachineCleanYear,lastMachineCleanYear);
        } else {
          printf("tempMachineCleanCount %d epoch %f startMachineCleanYear %f, lastMachineCleanYear %f %f\n",tempMachineCleanCount,pCleaning->epoch,startMachineCleanYear,lastMachineCleanYear,pLastMachineCleaning->epoch);
        }
#endif
        if (startMachineCleanYear < 0) {
          startMachineCleanYear = pCleaning->epoch;
          lastMachineCleanYear = pCleaning->epoch;
          tempMachineCleanCount = 1;
          pLastMachineCleaning = pCleaning;
          
        } else {
          if ((pCleaning->epoch - lastMachineCleanYear) < 0) {
            printf("ERROR: line %d year not in order at %s\n",__LINE__,pCleaning->date);
          }

          if ((pCleaning->epoch - lastMachineCleanYear) < CLEAN_IDLE_INTERVAL) {
            tempMachineCleanCount++;
            lastMachineCleanYear = pCleaning->epoch;
            pLastMachineCleaning = pCleaning;
          } else {
            if ((pLastMachineCleaning->epoch-startMachineCleanYear) <= 0) {
              printf("ERROR: line %d zero interval %s\n",__LINE__,pLastMachineCleaning->date);
            } else {
#if 1
              if (tempMachineCleanCount == 345) {
                printf("line %d tempMachineCleanCount %d\n",__LINE__,tempMachineCleanCount);
              }
#endif
              if (tempMachineCleanCount >= MIN_SESSION_COUNT) {
                pLastMachineCleaning->sessionCount = tempMachineCleanCount;
                pLastMachineCleaning->hours = (pLastMachineCleaning->epoch-startMachineCleanYear)*365.25*24.0;
                pLastMachineCleaning->rate = (1.0*pLastMachineCleaning->sessionCount)/pLastMachineCleaning->hours;
                if (pLastMachineCleaning->rate < CLEANING_MAX_RATE) {
                  totalMachineCleanCount += pLastMachineCleaning->sessionCount;
                  totalMachineCleanHours += pLastMachineCleaning->hours;
#if 1
                  printf("line %d sessionCount %d hours %f rate %f epoch %f ymax %f machine clean\n",__LINE__,pLastMachineCleaning->sessionCount,pLastMachineCleaning->hours,pLastMachineCleaning->rate,pLastMachineCleaning->epoch,ymax);
#endif
                  if (pLastMachineCleaning->rate > ymax) {
                    ymax = pLastMachineCleaning->rate;
                  }
                } else {
                  pLastMachineCleaning->sessionCount = 0;
                  pLastMachineCleaning->hours = 0;
                  pLastMachineCleaning->rate = 0;
                }
              }
            }
            startMachineCleanYear = pCleaning->epoch;
            lastMachineCleanYear = pCleaning->epoch;
            tempMachineCleanCount = 1;
            pLastMachineCleaning = pCleaning;
          }
        }

      }
    }


    if ((startHandCleanYear > 0) && (tempHandCleanCount > MIN_SESSION_COUNT)) {
      if ((pLastHandCleaning->epoch-startHandCleanYear) <= 0) {
        printf("ERROR: line %d zero interval %s\n",__LINE__,pCleaning->date);
      } else {
        
        pLastHandCleaning->sessionCount = tempHandCleanCount-1;
        pLastHandCleaning->hours = (pLastHandCleaning->epoch-startHandCleanYear)*365.25*24.0;
        pLastHandCleaning->rate = (1.0*pLastHandCleaning->sessionCount)/pLastHandCleaning->hours;
        if ( pLastHandCleaning->rate < CLEANING_MAX_RATE) {
          totalHandCleanCount += pLastHandCleaning->sessionCount;
          totalHandCleanHours += pLastHandCleaning->hours;
#if 1
          printf("line %d sessionCount %d hours %f rate %f epoch %f ymax %f machine clean\n",__LINE__,pLastHandCleaning->sessionCount,pLastHandCleaning->hours,pLastHandCleaning->rate,pLastHandCleaning->epoch,ymax);
#endif

          if (pLastHandCleaning->rate > ymax) {
            ymax = pLastHandCleaning->rate;
          }
        } else {
          pLastHandCleaning->sessionCount = 0;
          pLastHandCleaning->hours = 0;
          pLastHandCleaning->rate = 0;

        }
      }
    }
    if ((startMachineCleanYear > 0) && (tempMachineCleanCount > MIN_SESSION_COUNT)) {
      if ((pLastMachineCleaning->epoch-startMachineCleanYear) <= 0) {
        printf("ERROR: line %d zero interval %s\n",__LINE__,pCleaning->date);
      } else {
        pLastMachineCleaning->sessionCount = tempMachineCleanCount-1;
        pLastMachineCleaning->hours = (pLastMachineCleaning->epoch-startMachineCleanYear)*365.25*24.0;
        pLastMachineCleaning->rate = (1.0*pLastMachineCleaning->sessionCount)/pLastMachineCleaning->hours;
        if ( pLastMachineCleaning->rate < CLEANING_MAX_RATE) {

          totalMachineCleanCount += pLastMachineCleaning->sessionCount;
          totalMachineCleanHours += pLastMachineCleaning->hours;
#if 1
          printf("line %d sessionCount %d hours %f rate %f epoch %f ymax %f  machine clean\n",__LINE__,pLastMachineCleaning->sessionCount,pLastMachineCleaning->hours,pLastMachineCleaning->rate,pLastMachineCleaning->epoch,ymax);
#endif

          if (pLastMachineCleaning->rate > ymax) {
            ymax = pLastMachineCleaning->rate;
          }
        } else {
          pLastMachineCleaning->sessionCount = 0;
          pLastMachineCleaning->hours = 0;
          pLastMachineCleaning->rate = 0;

        }
      }
    }


  
#else /* CLEAN_IDLE_INTERVAL */

    if ( machineCleanCount > handCleanCount) {
      ymax = 1.0*machineCleanCount;
    } else {
      ymax = 1.0* handCleanCount;
    }
#endif /* CLEAN_IDLE_INTERVAL */

  }
  
  if (webuseFlag) {
    xmin =  minEpoch - 0.1;
    xmax = maxEpoch + 0.1;
    ymin = 0;
    ymax = nonharvardMaxMonth + 1.0;

  }
  if (linearityFlag) {
    xmin = 0;
    xmax = pLinearityCommon->maxValue;
    ymin = 0;
    ymax = 0;
    for (quadrant = 0; quadrant < pLinearityCommon->maxQuadrant; quadrant++) {
      for (actualValue = 0; actualValue < pLinearityCommon->maxValue; actualValue++) {
        result = pLinearityCommon->reverseBuffer[actualValue + (pLinearityCommon->maxValue * quadrant)] - actualValue;
        if (ymin > result) {
          ymin = result;
        }
        if (ymax < result) {
          ymax = result;
        }
      }
    }
  }

  if (tc720Flag) {
#ifdef PLOT_DAY
    xmin = (1.0*(mintc720.tp - PLOT_DAY))/(24.*3600.)-0.1;
    xmax = (1.0*(maxtc720.tp - PLOT_DAY))/(24.*3600.)+0.1;
#else 
    xmin = mintc720.hour;
    xmax = maxtc720.hour;
#endif
    ymin =  mintc720.T1;
    ymax =  maxtc720.T1;

#ifndef PLOT_T1_ONLY

    if (mintc720.T2 < ymin) {
      ymin = mintc720.T2;
    }
    if (maxtc720.T2 > ymax) {
      ymax = maxtc720.T2;
    }
#if 0
    if (mintc720.T3 < ymin) {
      ymin = mintc720.T3;
    }
    if (maxtc720.T3 > ymax) {
      ymax = maxtc720.T3;
    }
#endif
#ifdef PLOT_SETPOINT
    if (mintc720.setpoint < ymin) {
      ymin = mintc720.setpoint;
    }
    if (maxtc720.setpoint > ymax) {
      ymax = maxtc720.setpoint;
    }

#endif /* PLOT_SETPOINT */

    
#endif /* PLOT_T1_ONLY */
    
  }

  plsdev("png");
  plsetopt("-o",outfile);
  if (linearityFlag) {
    plsetopt("-geometry","1440x1080");
  }
  plscolbg(255,255,255);	/* Force the background colour to white */  
  plscol0(1, 0,0,0);		/* Force the foreground colour to black */
  plscol0(15,255,0,0);		/* Move red to 15 */
  plinit();
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax,just, 2); 
  plwind(xmin,xmax,ymin,ymax);
  plgspa(&xminmm,&xmaxmm,&yminmm,&ymaxmm);
  printf("xminmm %f xmaxmm %f yminmm %f ymaxmm %f\n",xminmm,xmaxmm,yminmm,ymaxmm);
  printf("xmin %f xmax %f ymin %f ymax %f\n",xmin,xmax,ymin,ymax);

  plcol0(1);
  if (limitingMagnitudeFlag) {
    pllab("Deepest Limiting Magnitude","Percent of Plates",title);
    for (statsIndex = 0; statsIndex < labelSize; statsIndex++) {
      pLabelEntry = &statsLabels[statsIndex];
      pLimitingStats = pLabelEntry->statsTable;
      if (pLimitingStats->count == 0) {
        continue;
      }
      /* Now convert the histogram to a cumulative histogram with 100% at the brightest magnitude */
      for (limitingMagIndex = 0; limitingMagIndex <= HISTOGRAM_SIZE; limitingMagIndex++) {
        x[limitingMagIndex] = (1.0*limitingMagIndex)/(1.0*MAGNITUDE_MULTIPLIER);
        y[limitingMagIndex] = (100.0* pLimitingStats->histogram[limitingMagIndex]) / (1.0*pLimitingStats->count);
#if 0
        printf("statsIndex %d limitingMagIndex %3d x %4.1f y %6.2f\n",statsIndex,limitingMagIndex,x[limitingMagIndex],y[limitingMagIndex]);
#endif
      }
      plcol0(pLabelEntry->color);
      pllsty(pLabelEntry->style);
      plptex(pLabelEntry->xLabelPos,pLabelEntry->yLabelPos,1.0,0.0,0.0,pLabelEntry->plotLabel);
      plline(HISTOGRAM_SIZE,x,y);
      x[0] = pLabelEntry->xLabelPos;
      y[0] = pLabelEntry->yLabelPos;
      x[1] = pLabelEntry->xArrowPos;
      limitingMagIndex = (pLabelEntry->xArrowPos * MAGNITUDE_MULTIPLIER) + 0.5;
      y[1] = y[limitingMagIndex];
      plline(2,x,y);
    }
  }
  if (scangraphFlag) {
#ifdef BASE_YEAR
    sprintf(xlabel,"Year - %.1f",BASE_YEAR);
#else /* BASE_YEAR */
    strcpy(xlabel,"Year");
#endif /* BASE_YEAR */
    pllab(xlabel,"Plate Scans",title);
    plcol0(15);
    for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
      pScangraph = &scangraphTable[pointIndex];
#ifdef BASE_YEAR
      x[pointIndex] = pScangraph->epoch - BASE_YEAR;
#else /* BASE_YEAR */
      x[pointIndex] = pScangraph->epoch;
#endif /* BASE_YEAR */
      y[pointIndex] = pScangraph->count;
    }
    plline(curCount-1,x,y);

  }

  if (scanSmoothedRateFlag) {
#ifdef BASE_YEAR
    sprintf(xlabel,"Year - %.1f",BASE_YEAR);
#else /* BASE_YEAR */
    strcpy(xlabel,"Year");
#endif /* BASE_YEAR */
    pllab(xlabel,"Plates per Week",title);
    plcol0(15);
    for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
      pScangraph = &scangraphTable[pointIndex];
#ifdef BASE_YEAR
      x[pointIndex] = pScangraph->midepoch - BASE_YEAR;
#else /* BASE_YEAR */
      x[pointIndex] = pScangraph->midepoch;
#endif /* BASE_YEAR */
      y[pointIndex] = pScangraph->rate;
#if 0
      printf("point %d x %f y %f\n",pointIndex,x[pointIndex],y[pointIndex]);
#endif
    }
    plline(curCount-1,x,y);
    plcol0(9);
    sprintf(tempBuffer,"Smoothing period is %d days.",smoothingDays);
    plptex(xmin+0.5,(ymax+ymin)/2.0,0.0,0.0,0,tempBuffer);
  }


  if (scanrateFlag) {
    pllab("Plates Scanned","Plates per Hour",title);
    for (switchCounter = 0; switchCounter < 9; switchCounter++) {
      for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
        pScanrate = &scanrateTable[pointIndex];
        x[pointIndex] = pScanrate->count;
        switch(switchCounter) {
        case 0:
          y[pointIndex] = pScanrate->minimum;
          break;
        case 8:
          y[pointIndex] = pScanrate->average;
          break;
        case 7:
          y[pointIndex] = pScanrate->median;
          break;
        case 6:
          y[pointIndex] = pScanrate->pct600;
          break;
        case 5:
          y[pointIndex] = pScanrate->pct700;
          break;
        case 4:
          y[pointIndex] = pScanrate->pct800;
          break;
        case 3:
          y[pointIndex] = pScanrate->pct900;
          break;
        case 2:
          y[pointIndex] = pScanrate->pct950;
          break;
        case 1:
          y[pointIndex] = pScanrate->pct990;
          break;
        default:
          printf("ERROR bad switch %d\n",switchCounter);
          exit(-1);
        }
      }
      plcol0(15);
      plline(curCount-1,x,y);
      plcol0(1);
#if 1
      plptex(210,69-(2*switchCounter),1.0,0.0,0.0,scanText[switchCounter]);
#endif
#if 0 /* Large Plates */
      plptex(80.,21.5-(0.5*switchCounter),1.0,0.0,0.0,scanText[switchCounter]);
#endif
#if 0 /* Small Plates */
      plptex(300.,125.-(4.*switchCounter),1.0,0.0,0.0,scanText[switchCounter]);
#endif
      x[0] = x[30];
      y[0] = y[30];
#if 1
      x[1] = 210.;
      y[1] = 69-(2*switchCounter);
#endif
#if 0 /* Large Plates */
      x[1] = 80.;
      y[1] = 21.5-(0.5*switchCounter);
#endif
#if 0 /* Small Plates */
      x[1] = 300.;
      y[1] = 125.0-(4.*switchCounter);
#endif
      plline(2,x,y);

    }
  }


  if (lightcurveRMSFlag) {
    qsort((void*)lightcurveRMSTable,curCount,sizeof(LIGHTCURVERMS),LightcurveRMSCompare);
    pllab("Magnitude","Median Lightcurve RMS",title);
    for (catalogNumber = 0; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
      goodPoints = 0;
      switch(catalogNumber) {
      case 0:
        plcol0(15); /* red */
        break;
      case 1:
        plcol0(3);  /* green */
        break;
      case 2:
        plcol0(9);  /* blue */
        break;
      case 3:
        plcol0(11); /* cyan */
        break;
      case 4:
        plcol0(1); /* black */
        break;
      default:
        plcol0(0); 
        break;
      }
      for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
        pLightcurveRMS = &lightcurveRMSTable[pointIndex];
        if (pLightcurveRMS->catalogNumber != catalogNumber) {
          continue;
        }
        x[goodPoints] = (1.0*(pLightcurveRMS->minmag+pLightcurveRMS->maxmag))/2.0;
        y[goodPoints] = pLightcurveRMS->median;
        if (goodPoints == 1) {
          plptex(x[goodPoints],y[goodPoints],0.0,0.0,0,catalogTextUppercase[catalogNumber]);
        }
        goodPoints++;
      }
      if (goodPoints == 0) {
        continue;
      }
      pllsty(catalogNumber+1);
      plline(goodPoints,x,y);
                        

      
    }
  }
  if (pointsPerCurveFlag) {
    qsort((void*)lightcurveRMSTable,curCount,sizeof(LIGHTCURVERMS),LightcurveRMSCompare);
    pllab("Magnitude","Median Points per Lightcurve",title);
    for (catalogNumber = 0; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
      goodPoints = 0;
      switch(catalogNumber) {
      case 0:
        plcol0(15); /* red */
        break;
      case 1:
        plcol0(3); /* green */
        break;
      case 2:
        plcol0(9); /* blue */
        break;
      case 3:
        plcol0(11); /* cyan */
        break;
      case 4:
        plcol0(1); /* black */
        break;

      default:
        plcol0(0);
        break;
      }
      for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
        pLightcurveRMS = &lightcurveRMSTable[pointIndex];
        if (pLightcurveRMS->catalogNumber != catalogNumber) {
          continue;
        }
        x[goodPoints] = (1.0*(pLightcurveRMS->minmag+pLightcurveRMS->maxmag))/2.0;
        y[goodPoints] = 1.0*pLightcurveRMS->npoints;
        if (x[goodPoints] == 9.5) {
          plptex(x[goodPoints],y[goodPoints],0.0,0.0,0,catalogTextUppercase[catalogNumber]);
        }
        goodPoints++;
      }
      if (goodPoints == 0) {
        continue;
      }
      pllsty(catalogNumber+1);
      plline(goodPoints,x,y);
                        

      
    }
  }
  if (clipRMSFlag) {
    pllab("Magnitude","Lightcurve Count",title);
    goodCurves = 0;
    for (curmag = minmag; curmag < (maxmag-1); curmag++) {
      peakIndex = 0;
      goodPoints = 0;
      for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
        pClipRMS = &clipRMSTable[pointIndex];
        if (pClipRMS->catalogNumber != catalogNumber) {
          continue;
        }
        if (pClipRMS->minmag != curmag) {
          continue;
        }
        x[goodPoints] = pClipRMS->clip_rms_local;
        y[goodPoints] = 1.0 * pClipRMS->count;
        if (y[goodPoints] > y[peakIndex]) {
          peakIndex = goodPoints;
        }
        goodPoints++;
        
      }
      if (goodPoints == 0) {
        continue;
      }
      
      if (y[peakIndex] > (0.02*maxCount)) {
        switch(goodCurves) {
        case 0:
          plcol0(4);
          break;
        case 1:
          plcol0(5);
          break;
        case 2:
          plcol0(6);
          break;
        case 3:
          plcol0(15);
          break;
        case 4:
          plcol0(9);
          break;
        case 5:
          plcol0(8);
          break;
        case 6:
          plcol0(3);
          break;
        case 7:
          plcol0(11);
          break;
        case 8:
          plcol0(13);
          break;

        default:
          plcol0(1);
          break;
        }
        pllsty((goodCurves % 8)+1);
        plline(goodPoints,x,y);
        sprintf(magnitudeString,"Mag %d-%d",curmag,curmag+1);
        printf("%s x %f y %f goodCurves %d\n",magnitudeString,x[peakIndex],y[peakIndex],goodCurves);
        plptex(x[peakIndex],y[peakIndex],0.0,0.0,0,magnitudeString);
        goodCurves++;
      }
                        

      
    }

  }
  if (measPerMagFlag) {
    pllab("Magnitude","Magnitude Measurement Count",title);
    for (sequence = 1; sequence <= 2; sequence++) {
      for (catalogNumber = CATALOG_GSC232; catalogNumber <= CATALOG_ATLAS; catalogNumber++) {
        if (catalogNumber == CATALOG_KEPLER) {
          continue;
        }
        goodPoints = 0;
        for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
          pMeasPerMag = &measPerMagTable[pointIndex];
          if (pMeasPerMag->catalogNumber != catalogNumber) {
            continue;
          }
          if (pMeasPerMag->sequence != sequence) {
            continue;
          }
          x[goodPoints] = 1.0 * pMeasPerMag->magcal_magdep;
          y[goodPoints] = 1.0 * pMeasPerMag->count;

          if ((sequence == 1) && (catalogNumber == CATALOG_GSC232) && (pMeasPerMag->magcal_magdep == 11)) {
            labelIndex = goodPoints;
          } else if ((sequence == 2) && (catalogNumber == CATALOG_GSC232) && (pMeasPerMag->magcal_magdep == 16)) {
            labelIndex = goodPoints;
          } else if ((sequence == 1) && (catalogNumber == CATALOG_APASS) && (pMeasPerMag->magcal_magdep == 14)) {
            labelIndex = goodPoints;
          } else if ((sequence == 2) && (catalogNumber == CATALOG_APASS) && (pMeasPerMag->magcal_magdep == 16)) {
            labelIndex = goodPoints;
          } else if ((sequence == 1) && (catalogNumber == CATALOG_GAIA) && (pMeasPerMag->magcal_magdep == 14)) {
            labelIndex = goodPoints;
          } else if ((sequence == 2) && (catalogNumber == CATALOG_GAIA) && (pMeasPerMag->magcal_magdep == 16)) {
            labelIndex = goodPoints;
          } else if ((sequence == 1) && (catalogNumber == CATALOG_ATLAS) && (pMeasPerMag->magcal_magdep == 13)) {
            labelIndex = goodPoints;
          } else if ((sequence == 2) && (catalogNumber == CATALOG_ATLAS) && (pMeasPerMag->magcal_magdep == 14)) {
            labelIndex = goodPoints;
          }
          goodPoints++;
        
        }
        if (goodPoints == 0) {
          continue;
        }
        printf("sequence %d catalogNumber %d\n",sequence,catalogNumber);
        if ((sequence == 1) && (catalogNumber == CATALOG_GSC232)) {
          plcol0(15);
          pllsty(1);
          justify = 1.0;
          strcpy(tempBuffer,"GSC2.3.2 Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_GSC232)) {
          plcol0(3);
          pllsty(2);
          justify = 0.0;
          strcpy(tempBuffer,"GSC2.3.2 Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_APASS)) {
          plcol0(9);
          pllsty(3);
          justify = 0.0;
          strcpy(tempBuffer,"APASS Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_APASS)) {
          plcol0(13);
          pllsty(4);
          justify = 0.0;
          strcpy(tempBuffer,"APASS Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_GAIA)) {
          plcol0(9);
          pllsty(3);
          justify = 0.0;
          strcpy(tempBuffer,"GAIA Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_GAIA)) {
          plcol0(13);
          pllsty(4);
          justify = 0.0;
          strcpy(tempBuffer,"GAIA Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_ATLAS)) {
          plcol0(1);
          pllsty(3);
          justify = 0.0;
          strcpy(tempBuffer,"ATLAS Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_ATLAS)) {
          plcol0(13);
          pllsty(4);
          justify = 0.0;
          strcpy(tempBuffer,"ATLAS Non-Patrol");
        } else {
          printf("ERROR: illegal sequence %d or catalog %d\n",sequence,catalogNumber);
          exit(-1);
        }
        printf("sequence %d catalogNumber %d (2)\n",sequence,catalogNumber);
        plline(goodPoints,x,y);
        plptex(x[labelIndex],y[labelIndex],0.0,0.0,justify,tempBuffer);
        goodCurves++;
      }
                        
    }
      
  }



  if (limMagPerYearFlag) {
    pllab("Year","Average Limiting Magnitude",title);
    for (sequence = 1; sequence <= 2; sequence++) {
      for (catalogNumber = CATALOG_GSC232; catalogNumber <= CATALOG_ATLAS; catalogNumber++) {
        if (catalogNumber == CATALOG_KEPLER) {
          continue;
        }
        goodPoints = 0;
        for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
          pLimMagPerYear = &limMagPerYearTable[pointIndex];
          if (pLimMagPerYear->catalogNumber != catalogNumber) {
            continue;
          }
          if (pLimMagPerYear->sequence != sequence) {
            continue;
          }
#if 0
          if ((pLimMagPerYear->year >= 1934) && (pLimMagPerYear->year <= 1952)) {
            printf("WARNING: skipping year %d\n",pLimMagPerYear->year);
            continue;
          }
#endif
          x[goodPoints] = 1.0 * pLimMagPerYear->year;
          y[goodPoints] = 1.0 * pLimMagPerYear->limiting_mag_local;

          if ((sequence == 1) && (catalogNumber == CATALOG_GSC232) && (pLimMagPerYear->year == 1958)) {
            labelIndex = goodPoints; /* GSC patrol */
          } else if ((sequence == 2) && (catalogNumber == CATALOG_GSC232) && (pLimMagPerYear->year == 1955)) {
            labelIndex = goodPoints; /* GSC non-patrol */
          } else if ((sequence == 1) && (catalogNumber == CATALOG_APASS) && (pLimMagPerYear->year == 1958)) {
            labelIndex = goodPoints; /* APASS patrol */
          } else if ((sequence == 2) && (catalogNumber == CATALOG_APASS) && (pLimMagPerYear->year == 1938)) {
            labelIndex = goodPoints; /* APASS non-patrol */
          } else if ((sequence == 1) && (catalogNumber == CATALOG_GAIA) && (pLimMagPerYear->year == 1958)) {
            labelIndex = goodPoints; /* GAIA patrol */
          } else if ((sequence == 2) && (catalogNumber == CATALOG_GAIA) && (pLimMagPerYear->year == 1928)) {
            labelIndex = goodPoints; /* GAIA non-patrol */
          } else if ((sequence == 1) && (catalogNumber == CATALOG_ATLAS) && (pLimMagPerYear->year == 1980)) {
            labelIndex = goodPoints; /* ATLAS patrol */
          } else if ((sequence == 2) && (catalogNumber == CATALOG_ATLAS) && (pLimMagPerYear->year == 1928)) {
            labelIndex = goodPoints; /* ATLAS non-patrol */
          }
          goodPoints++;
        
        }
        if (goodPoints == 0) {
          continue;
        }
        printf("sequence %d catalogNumber %d\n",sequence,catalogNumber);
        if ((sequence == 1) && (catalogNumber == CATALOG_GSC232)) {
          plcol0(15);
          pllsty(1);
          justify = 1.0;
          x1[0] = 1950.;
          y1[0] = 12.5;
          strcpy(tempBuffer,"GSC2.3.2 Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_GSC232)) {
          plcol0(3);
          pllsty(2);
          justify = 1.0;
          x1[0] = 1925.;
          y1[0] = 17.5;
          strcpy(tempBuffer,"GSC2.3.2 Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_APASS)) {
          plcol0(9);
          pllsty(3);
          justify = 1.0;
          x1[0] = 1950.;
          y1[0] = 13.0;
          strcpy(tempBuffer,"APASS Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_APASS)) {
          plcol0(13);
          pllsty(4);
          justify = 1.0;
          x1[0] = 1925.;
          y1[0] = 17.0;
          strcpy(tempBuffer,"APASS Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_GAIA)) {
          plcol0(9);
          pllsty(3);
          justify = 1.0;
          x1[0] = 1950.;
          y1[0] = 13.0;
          strcpy(tempBuffer,"GAIA Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_GAIA)) {
          plcol0(13);
          pllsty(4);
          justify = 1.0;
          x1[0] = 1925.;
          y1[0] = 15.0;
          strcpy(tempBuffer,"GAIA Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_ATLAS)) {
          plcol0(1);
          pllsty(3);
          justify = 1.0;
          x1[0] = 1990.;
          y1[0] = 13.0;
          strcpy(tempBuffer,"ATLAS Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_ATLAS)) {
          plcol0(1);
          pllsty(4);
          justify = 1.0;
          x1[0] = 1925.;
          y1[0] = 15.0;
          strcpy(tempBuffer,"ATLAS Non-Patrol");
        } else {
          printf("ERROR: illegal sequence %d or catalog %d\n",sequence,catalogNumber);
          exit(-1);
        }
        printf("sequence %d catalogNumber %d (2)\n",sequence,catalogNumber);
        x1[1] = x[labelIndex];
        y1[1] = y[labelIndex];
        plline(2,x1,y1);
        plptex(x1[0],y1[0],0.0,0.0,justify,tempBuffer);
#if 0
        for (pointIndex = 0; pointIndex < goodPoints; pointIndex++) {
          printf("pointIndex %d, x %f y %f\n",pointIndex,x[pointIndex],y[pointIndex]);
        }
#endif
        plline(goodPoints,x,y); 
        goodCurves++;
      }
                        
    }
      
  }



  if (countPerYearFlag) {
    pllab("Year","Number of Magnitude Measurements",title);
    for (sequence = 1; sequence <= 2; sequence++) {
      for (catalogNumber = CATALOG_GSC232; catalogNumber <= CATALOG_ATLAS; catalogNumber++) {
        if (catalogNumber == CATALOG_KEPLER) {
          continue;
        }
        goodPoints = 0;
        for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
          pCountPerYear = &countPerYearTable[pointIndex];
          if (pCountPerYear->catalogNumber != catalogNumber) {
            continue;
          }
          if (pCountPerYear->sequence != sequence) {
            continue;
          }

#if 0
          if ((pCountPerYear->year >= 1934) && (pCountPerYear->year <= 1952)) {
            printf("WARNING: skipping year %d\n",pCountPerYear->year);
            continue;
          }
#endif
          x[goodPoints] = 1.0 * pCountPerYear->year;
          y[goodPoints] = 1.0 * pCountPerYear->count;

          if ((sequence == 1) && (catalogNumber == CATALOG_GSC232) && (pCountPerYear->year == 1940)) {
            labelIndex = goodPoints; /* GSC patrol */
          } else if ((sequence == 2) && (catalogNumber == CATALOG_GSC232) && (pCountPerYear->year == 1896)) {
            labelIndex = goodPoints; /* GSC non-patrol */
          } else if ((sequence == 1) && (catalogNumber == CATALOG_APASS) && (pCountPerYear->year == 1984)) {
            labelIndex = goodPoints; /* APASS patrol */
          } else if ((sequence == 2) && (catalogNumber == CATALOG_APASS) && (pCountPerYear->year == 1934)) {
            labelIndex = goodPoints; /* APASS non-patrol */
          } else if ((sequence == 1) && (catalogNumber == CATALOG_GAIA) && (pCountPerYear->year == 1921)) {
            labelIndex = goodPoints; /* GAIA patrol */
          } else if ((sequence == 2) && (catalogNumber == CATALOG_GAIA) && (pCountPerYear->year == 1935)) {
            labelIndex = goodPoints; /* GAIA non-patrol */
          } else if ((sequence == 1) && (catalogNumber == CATALOG_ATLAS) && (pCountPerYear->year == 1934)) {
            labelIndex = goodPoints; /* ATLAS patrol */
          } else if ((sequence == 2) && (catalogNumber == CATALOG_ATLAS) && (pCountPerYear->year == 1935)) {
            labelIndex = goodPoints; /* ATLAS non-patrol */
          }
          goodPoints++;
        
        }
        if (goodPoints == 0) {
          continue;
        }
        printf("sequence %d catalogNumber %d\n",sequence,catalogNumber);
        if ((sequence == 1) && (catalogNumber == CATALOG_GSC232)) {
          plcol0(15);
          pllsty(1);
          justify = 0.0;
          x1[0] = 1950.;
          y1[0] = ymax*0.9;
          strcpy(tempBuffer,"GSC2.3.2 Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_GSC232)) {
          plcol0(3);
          pllsty(2);
          justify = 0.0;
          x1[0] = 1890.;
          y1[0] = ymax*0.5;
          strcpy(tempBuffer,"GSC2.3.2 Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_APASS)) {
          plcol0(9);
          pllsty(3);
          justify = 1.0;
          x1[0] = 1980.;
          y1[0] = 10.5;
          strcpy(tempBuffer,"APASS Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_APASS)) {
          plcol0(13);
          pllsty(4);
          justify = 1.0;
          x1[0] = 1925.;
          y1[0] = 17.0;
          strcpy(tempBuffer,"APASS Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_GAIA)) {
          plcol0(9);
          pllsty(3);
          justify = 1.0;
          x1[0] = 1920.;
          y1[0] = ymax*0.5;
          strcpy(tempBuffer,"GAIA Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_GAIA)) {
          plcol0(13);
          pllsty(4);
          justify = 1.0;
          x1[0] = 1925.;
          y1[0] = ymax*0.6;
          strcpy(tempBuffer,"GAIA Non-Patrol");
        } else if ((sequence == 1) && (catalogNumber == CATALOG_ATLAS)) {
          plcol0(1);
          pllsty(3);
          justify = 1.0;
          x1[0] = 1920.;
          y1[0] = ymax*0.55;
          strcpy(tempBuffer,"ATLAS Patrol");
        } else if ((sequence == 2) && (catalogNumber == CATALOG_ATLAS)) {
          plcol0(13);
          pllsty(4);
          justify = 1.0;
          x1[0] = 1925.;
          y1[0] = ymax*0.6;
          strcpy(tempBuffer,"ATLAS Non-Patrol");
        } else {
          printf("ERROR: illegal sequence %d or catalog %d\n",sequence,catalogNumber);
          exit(-1);
        }
        printf("sequence %d catalogNumber %d (2)\n",sequence,catalogNumber);
        x1[1] = x[labelIndex];
        y1[1] = y[labelIndex];
        plline(2,x1,y1);
        plptex(x1[0],y1[0],0.0,0.0,justify,tempBuffer);
#if 0
        for (pointIndex = 0; pointIndex < goodPoints; pointIndex++) {
          printf("pointIndex %d, x %f y %f\n",pointIndex,x[pointIndex],y[pointIndex]);
        }
#endif
        plline(goodPoints,x,y); 
        goodCurves++;
      }
                        
    }
      
  }
  if (cleaningFlag) {
#ifdef  CLEAN_IDLE_INTERVAL
#ifdef CLEAN_AVERAGE_DAYS 
    pllab("Year","Average Plates/Hour",title);
#else /* CLEAN_AVERAGE_DAYS */
    pllab("Year","Plates/Hour",title);
#endif /* CLEAN_AVERAGE_DAYS */
    goodPoints = 0;
    lastPointIndex = 0;
    pLastCleaning =  pCleaning = &cleaningTable[lastPointIndex];
    cumulativeSessionCount = 0;
    cumulativeHoursCount = 0;
    for (pointIndex = 0; pointIndex < maxCleaningCount; pointIndex++) {
      pCleaning = &cleaningTable[pointIndex];
      if ((pCleaning->sequence == HAND_CLEAN) && (pCleaning->rate > 0)) {
#ifdef CLEAN_AVERAGE_DAYS
        cumulativeSessionCount += pCleaning->sessionCount;
        cumulativeHoursCount += pCleaning->hours;
        while ((pCleaning->epoch-pLastCleaning->epoch) >  (CLEAN_AVERAGE_DAYS/365.25)) {
          if (pLastCleaning->sequence == HAND_CLEAN) {
            cumulativeSessionCount -=  pLastCleaning->sessionCount;
            cumulativeHoursCount -= pLastCleaning->hours;
          }
          lastPointIndex++;
          if (lastPointIndex >= maxCleaningCount) {
            printf("ERROR: line %d index overflow\n",__LINE__);
            exit(-1);
          }
          pLastCleaning = &cleaningTable[lastPointIndex];          
        }
        if (cumulativeHoursCount < 0) {
          printf("ERROR line %d negative cumulativeHoursCount\n",__LINE__);
          exit(-1);
        }
        if (cumulativeHoursCount == 0) {
          continue;
        }
        pCleaning->rate = (1.0*cumulativeSessionCount)/cumulativeHoursCount;
#endif /* CLEAN_AVERAGE_DAYS */        

        x[goodPoints] = pCleaning->epoch;
        y[goodPoints] = pCleaning->rate;
#if 0
        printf("line %d x %f y %f\n",__LINE__,x[goodPoints],y[goodPoints]);
#endif
        goodPoints++;
      }
    }
    plcol0(15); /* red */
#if 1
    plpoin(goodPoints,x,y,2); 
#else
    plline(goodPoints,x,y); 
#endif
    goodPoints = 0;
    lastPointIndex = 0;
    pLastCleaning =  pCleaning = &cleaningTable[lastPointIndex];
    cumulativeSessionCount = 0;
    cumulativeHoursCount = 0;
    for (pointIndex = 0; pointIndex < maxCleaningCount; pointIndex++) {
#if 0
      if (pointIndex == 112905) {
        printf("line %d\n",__LINE__);
      }
#endif
      pCleaning = &cleaningTable[pointIndex];
      if ((pCleaning->sequence == MACHINE_CLEAN) && (pCleaning->rate > 0)) {
#ifdef CLEAN_AVERAGE_DAYS
        cumulativeSessionCount += pCleaning->sessionCount;
        cumulativeHoursCount += pCleaning->hours;
        while ((pCleaning->epoch-pLastCleaning->epoch) >  (CLEAN_AVERAGE_DAYS/365.25)) {
          if (pLastCleaning->sequence == MACHINE_CLEAN) {
            cumulativeSessionCount -=  pLastCleaning->sessionCount;
            cumulativeHoursCount -= pLastCleaning->hours;
          }
          lastPointIndex++;
          if (lastPointIndex >= maxCleaningCount) {
            printf("ERROR: line %d index overflow\n",__LINE__);
            exit(-1);
          }
          pLastCleaning = &cleaningTable[lastPointIndex];          
        }
        if (cumulativeHoursCount < 0) {
          printf("ERROR line %d negative cumulativeHoursCount\n",__LINE__);
          exit(-1);
        }
#if 0
        printf("line %d lastPointIndex %d pointIndex %d cumulativeSessionCount %d cumulativeHoursCount %f epoch %f maxCleaningCount %d goodPoints %d\n",
               __LINE__,
               lastPointIndex,
               pointIndex,
               cumulativeSessionCount,
               cumulativeHoursCount,
               pCleaning->epoch,
               maxCleaningCount,
               goodPoints);
#endif
        if (cumulativeHoursCount == 0) {
          continue;
        }
        pCleaning->rate = (1.0*cumulativeSessionCount)/cumulativeHoursCount;
#endif /* CLEAN_AVERAGE_DAYS */        


        x[goodPoints] = pCleaning->epoch;
        y[goodPoints] = pCleaning->rate;
#if 0
        printf("line %d x %f y %f\n",__LINE__,x[goodPoints],y[goodPoints]);
#endif
        goodPoints++;
      }
    }
    plcol0(9);  /* blue */
#if 1
    plpoin(goodPoints,x,y,5); 
#else
    plline(goodPoints,x,y); 
#endif
    
    printf("totalHandCleanCount %d totalHandCleanHours %f aveHandCleaningRate %f\n",totalHandCleanCount,totalHandCleanHours,(1.0*totalHandCleanCount)/totalHandCleanHours);
    printf("totalMachineCleanCount %d totalMachineCleanHours %f aveMachineCleaningRate %f\n",totalMachineCleanCount,totalMachineCleanHours,(1.0*totalMachineCleanCount)/totalMachineCleanHours);


#else /* CLEAN_IDLE_INTERVAL */
    pllab("Year","Count",title);
    goodPoints = 0;
    for (pointIndex = 0; pointIndex < maxCleaningCount; pointIndex++) {
      pCleaning = &cleaningTable[pointIndex];
      if (pCleaning->sequence == HAND_CLEAN) {
        x[goodPoints] = pCleaning->epoch;
        y[goodPoints] = pCleaning->count;
        printf("line %d x %f y %f\n",__LINE__,x[goodPoints],y[goodPoints]);
        goodPoints++;
      }
    }
    plcol0(15); /* red */
    plline(goodPoints,x,y); 
    goodPoints = 0;
    for (pointIndex = 0; pointIndex < maxCleaningCount; pointIndex++) {
      pCleaning = &cleaningTable[pointIndex];
      if (pCleaning->sequence == MACHINE_CLEAN) {
        x[goodPoints] = pCleaning->epoch;
        y[goodPoints] = pCleaning->count;
        printf("line %d x %f y %f\n",__LINE__,x[goodPoints],y[goodPoints]);
        goodPoints++;
      }
    }
    plcol0(9);  /* blue */
    plline(goodPoints,x,y); 
#endif /* CLEAN_IDLE_INTERVAL */

  }

  if (webuseFlag) {
    pllab("Year","Monthly searches",title);
    goodPoints = 0;
    for (monthIndex = 0; monthIndex <= maxMonthIndex; monthIndex++) {
      x[goodPoints] = ((1.0*monthIndex) / 12.0) + BASE_DATE;
      y[goodPoints] = nonharvardArray[monthIndex];
      goodPoints++;
    }
    plline(goodPoints,x,y); 
  }



  if (tc720Flag) {
#ifdef PLOT_DAY
    baseTime = PLOT_DAY;
    ptr = localtime(&baseTime);
    strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
    sprintf(title2,"Day since %s",timestr);
#else /* PLOT_DAY */
    strcpy(title2,"Eastern Time Hour");
#endif /* PLOT_DAY */

#ifndef USE_CENTIGRADE
    pllab(title2,"Degrees F",title);
#else /* USE_CENTIGRADE */
    pllab(title2,"Degrees C",title);
#endif /* USE_CENTIGRADE */

    goodPoints = 0;
    color = 3;
    for (pointIndex = 0; pointIndex < curCount; pointIndex++) {
      ptc720 = &tc720Table[pointIndex];
#ifdef PLOT_DAY
      x[goodPoints] = (1.0*(ptc720->tp-PLOT_DAY))/(24.0*3600.);
#else /* PLOT_DAY */
      x[goodPoints] = ptc720->hour;
#endif /* PLOT_DAY */
      T1[goodPoints] = ptc720->T1;
      T2[goodPoints] = ptc720->T2;
      T3[goodPoints] = ptc720->T3;
      SETPOINT[goodPoints] = ptc720->setpoint;
      OUTPUT[goodPoints] = ymax - 0.1 + ((ymax-ymin-0.2) * ptc720->output)/100;
      goodPoints++;
      nxt720 = &tc720Table[pointIndex+1];
      if ((pointIndex == curCount-1) ||
          ((nxt720->tp - ptc720->tp) > 10)) { /* no more than a 10 second gap */
        printf("Gap of %f hours at ptc720->timestr\n",(1.0*(nxt720->tp - ptc720->tp))/3600.,ptc720->timestr);
#ifdef PLOT_T1_ONLY
        plcol0(color);/* Set color to red */
        plline(goodPoints,x,T1);
        color++;
        if (color > 15) {
          color = 3;
        }
#else /* PLOT_T1_ONLY */
        plcol0(15);/* Set color to red */
        plline(goodPoints,x,T1); 
        plcol0(3); /* green */
        plline(goodPoints,x,T2); 
#if 0
        plcol0(9); /* blue */
        plline(goodPoints,x,T3);
#endif
#ifdef PLOT_SETPOINT
        plcol0(2); /* yellow */
        plline(goodPoints,x,SETPOINT);
#endif /* PLOT_SETPOINT */
#ifdef PLOT_POWER
        plcol0(5); /* pink */
        plline(goodPoints,x,OUTPUT);

#endif /* PLOT_POWER */



#endif /* PLOT_T1_ONLY */
        goodPoints = 0;
      }


    }
  }


  if (linearityFlag) {
#define LABEL_X_VALUE 1500
    ptr = localtime(&pLinearityCommon->scanTime);
    strftime(timestr, 25, "%Y-%m-%dT%H-%M-%S", ptr);
    pllab("Measured ADU","Actual ADU - Measured ADU",timestr);
    for (quadrant = 0; quadrant < pLinearityCommon->maxQuadrant; quadrant++) {
      for (actualValue = 0; actualValue < pLinearityCommon->maxValue; actualValue++) {
        result = pLinearityCommon->reverseBuffer[actualValue + (pLinearityCommon->maxValue * quadrant)] - actualValue;
        x[actualValue] = actualValue;
        y[actualValue] = result;
      }
      result = pLinearityCommon->reverseBuffer[LABEL_X_VALUE + (pLinearityCommon->maxValue * quadrant)] - LABEL_X_VALUE;
      switch (quadrant) {
      case 0:
        plcol0(1); /* black */
        plptex(LABEL_X_VALUE,result,0.0,0.0,0,"Q0");
        break;
      case 1:
        plcol0(15); /* red */
        plptex(LABEL_X_VALUE,result,0.0,0.0,0,"Q1");
        break;
      case 2:
        plcol0(9); /* blue */
        plptex(LABEL_X_VALUE,result,0.0,0.0,0,"Q2");
        break;
      case 3:
        plcol0(3); /* green */
        plptex(LABEL_X_VALUE,result,0.0,0.0,0,"Q3");
        break;
      default:
        plcol0(7); /* grey */
      }
      pllsty(quadrant+1);
      plline(pLinearityCommon->maxValue,x,y); 
    }

  }


  plend();

  if (scangraphTable != NULL) {
    free(scangraphTable);
  }
  if (scanrateTable != NULL) {
    free(scanrateTable);
  }
  if (webuseTable != NULL) {
    free(webuseTable);
  }
  if (webstatsTable != NULL) {
    free(webstatsTable);
  }
  if (webhandle != NULL) {
    fclose(webhandle);
  }
  if (lightcurveRMSTable != NULL) {
    free(lightcurveRMSTable);
  }
  if (linearityFlag) {
    InitLinearityCommon(pLinearityCommon);
  }
  return(0);
}
