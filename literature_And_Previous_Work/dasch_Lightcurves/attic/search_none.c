// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* search_none.c
 *
 * This routine steps through the entire photometry magnitudes table searching for objects not in the GSC2.3.2 catalog
 *
 * cc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include -I /usr/include/openmpi-i386 -I /usr/include/openmpi-x86_64 -L /n/sw/openmpi-1.4.2_gnu-4.1.2/lib  -L /usr/lib${lib64}/openmpi/lib search_none.c pipelineutils.a -L /dasch/install/lib -lm -lcfitsio    -ltable -lutil  -lwcs -pthread -o search_none -L/usr/lib${lib64}/mysql -lmysqlclient -ldl -l mpi
 *
 *  search_none -q apass -v -l -90 -c 1 -i   -o /dasch/Pipeline/ingest/search_none_apass_2014_03_05.log -d /dasch/Pipeline/ingest/search_none_apass_2014_03_05.db  -p /dasch/Pipeline/ingest/plate_outlier_apass_2014_03_05.db -u /dasch/Pipeline/ingest/id_unmatched_apass_2014_03_05.db -m /dasch/Pipeline/ingest/propermotion_apass_2014_03_05.log -a /dasch/Pipeline/ingest/dradhistogram_apass_2014_03_05.log
 *  mpirun -np 4 search_none -q apass -v -l -90 -c 1 -i   -o /dasch/Pipeline/ingest/search_none_apass_2014_03_06.log -d /dasch/Pipeline/ingest/search_none_apass_2014_03_06.db  -p /dasch/Pipeline/ingest/plate_outlier_apass_2014_03_06.db -u /dasch/Pipeline/ingest/id_unmatched_apass_2014_03_06.db -m /dasch/Pipeline/ingest/propermotion_apass_2014_03_06.log -a /dasch/Pipeline/ingest/dradhistogram_apass_2014_03_06.log -b /dasch/Pipeline/ingest/binhistogram_apass.txt
 *
 * Comparing databaseHandle (different group numbers)
 * awk '{OFS="\t"}{FS="\t"}{print $1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' search_none_apass_2014_03_05.db > los.tmp
 * awk '{OFS="\t"}{FS="\t"}{print $1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' search_none_apass_2014_03_06.db > los1.tmp

 * awk '{OFS="\t"}{FS="\t"}{print $2}' search_none_apass_2014_03_05.db | sort -u | wc   1292 entries
 * awk '{OFS="\t"}{FS="\t"}{print $2}' search_none_apass_2014_03_06.db | sort -u | wc   1292 entries

 * 
 *  Files written:
 *      
 *   -o /dasch/Pipeline/ingest/search_none_apass_2014_03_05.log         outHandle                                               CHECK total magnitudes incorrect, queries incorrect Images/Location incorrect
 *   -d /dasch/Pipeline/ingest/search_none_apass_2014_03_05.db          databaseHandle                                          OK
 *   -m /dasch/Pipeline/ingest/propermotion_apass_2014_03_05.log        InitProperMotions()   pFileCommon->propermotionHandle   OK
 *   -p /dasch/Pipeline/ingest/plate_outlier_apass_2014_03_05.db        plateHandle                                             OK
 *   -u /dasch/Pipeline/ingest/id_unmatched_apass_2014_03_05.db         idHandle                                                OK
 *   -a /dasch/Pipeline/ingest/dradhistogram_apass.log                  astrometryHandle                                        OK
 *      /dasch/Pipeline/ingest/search_none_apass_2014_03_05_console.log logHandle - for multiprocessing.                        CHECK  
 *   
 * Aug  4, 2009 Edward J. Los - Original Version
 * Aug 24, 2009 Edward J. Los - Add kepler catalog support
 * Dec 18, 2009 Edward J. Los - Add magnitude file support
 * Jan 13, 2010 Edward J. Los - Add Heliocentric Julian Date to the output file
 *                              Correct an indexing bug which missed a lot of stars.
 * Jan 14, 2010 Edward J. Los - Make the initial GSC bin index variable
 * Jan 15, 2010 Edward J. Los - Correct abort flag loop
 * Jan 18, 2010 Edward J. Los - Correct behavior for restart
 * Feb  5, 2010 Edward J. Los - Accept only points from plates that have no quality bits set
 *                              and that have a colorterm within range
 *                              Add kepler catalog support - use RA2fpix to display only entries on chips
 * Feb  9, 2010 Edward J. Los - Correct labelling of the histogram table.
 * Feb 13, 2010 Edward J. Los - Add spatial_bin to the output
 * Mar  1, 2010 Edward J. Los - Flag objects near bright GSC2.3.2 objects where the GSC catalog is likely to be incomplete.
 * Dec 20, 2010 Edward J. Los - Remove support for magnitudes stored in MySQL tables.
 * Dec 31, 2010 Edward J. Los - display minVersionId
 * Feb 14, 2011 Edward J. Los - Search for plates with many outliers.  From Sumin Tang's memo of 2/03/11 4:04 PM choose lightcurves with rms_local < 0.2 and
 *                              use nburst3  (number of good lightcurve points 0.4 mag above the median value) and ndip3 (number of good lightcurve points 0.4 mag 
 &                              below the median value
 * Mar 23, 2011 Edward J. Los   Add apass catalog support
 * Aug 31, 2011 Edward J. Los - Temporarily remove the REFNumber from the header sanity check (TEMP_REFNUMBER_FIX)
 * Nov  9, 2011 Edward J. Los - Add magnitude-dependent parameters (magdep_bin, magcal_magdep, and magcal_magdep_rms).
 * Nov 14, 2011 Edward J. Los - make minVersionId dependent on the catalog
 * Jul  3, 2011 Edward J. Los - Intercept sanity check early -- possible file corruption tripped it
 * Jul 23, 2012 Edward J. Los - Remove unused directory references
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Sep 10, 2012 Edward J. Los - Reorganize code to optimize magnitude file support.
 *                              Correct erroneously zero SFLAGS field
 * Sep 17, 2012 Edward J. Los - Expand search to adjacent bins.              Test case takes (466?) 796 seconds, originally 44 seconds
 *                              Generate an id table for unmatched objects.  Test case now takes 797 seconds -> 240 sec with optimization  -> 48 sec with additional optimization.
 * Dec  7, 2012 Edward J. Los - Add measurementCount, matchCount, and matchRatio to the plate outlier table.
 * Dec 10, 2012 Edward J. Los - Add PERMISSIVE_ID_TABLE to include more points in the id and summary tables
 * Jan  4, 2013 Edward J. Los - Replace nearestREF and nearestREFarcsec with nearbyObjects.
 * Feb 11, 2013 Edward J. Los - Add dradRMS2 and FWHM_WORLD to the output database file.
 * Mar 29, 2013 Edward J. Los - Add proper motion summary file
 * Apr  2, 2013 Edward J. Los - Add galactic limit qualifier
 * Nov  6, 2013 Edward J. Los - Create histograms dra and ddec for low drad stars
 * Dec 10, 2013 Edward J. Los - Disable OutputCounts
 * Dec 24, 2013 Edward J. Los - Add rootFlag to the group table
 * Mar  3, 2014 Edward J. Los - Convert to openmpi to improve performance
 * May  5, 2014 Edward J. Los - Use the bin histogram file for even partitioning of the work on the basis of total points 
 * Oct 25, 2014 Edward J. Los - Move all filtering to ProcessNoneImagesX so that the website choice of DASCH objects agrees with the id tables
 * Apr 13, 2015 Edward J. Los - Add galactic latitude to the id table
 * May  6, 2015 Edward J. Los - Add nearbyREFflag if a TC has a catalog star within FWHM/2 and the catalog object is at least 0.5 mag dimmer than the TC (modified Josh request of Tue 5/05/15 9:15 AM)
 *                              Also add releaseField, peakDays, peakYear, peakMag, peakSlope, peakSlopeRMS, and peakCount
 * Jun 15, 2015 Edward J. Los - Change peakSlopeRMS to peakRMS, add peakUpperCount, peakNumber, peakMaxMag, and peakMinMag
 * Jun 24, 2015 Edward J. Los - Add peakEventCount
 * Jun 29, 2015 Edward J. Los - Add peakDradRMS
 * Jul 13, 2015 Edward J. Los - Add peakEventCount2
 * Jul 17, 2015 Edward J. Los - Add peakEventCount3
 * Aug  1, 2015 Edward J. Los - Add -t parameter to enable transient flare searches
 * Aug 18, 2015 Edward J. Los - Add peakExcess = npoints outside of flare
 *                              Add "catalogDistance" histogram of distances of DASCH objects from catalog objects
 * Sep 15, 2015 Edward J. Los - Rename "peakExcess" to "peakExtra"
 *                              Add peakCountNF, peakCountWF, peakOutside, peakOutsideNF, and peakOutsideWF
 *                              Add daschunistd.h for table.h conflicts
 *                              Add peakRA and peakDec with 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window
 *                              Add peakDefectCount for the number of lightcurve defects in the selected transient window.
 * Sep 24, 2015 Edward J. Los - Add peakNearbyDistance for the distance to nearby stars that are within 2 magnitudes of the brightest or average magnitude of the transient candidate
 * Oct  4, 2015 Edward J. Los - Add peakMultipleCount for the number of multiple exposure plates in the selected transient window
 * Oct 26, 2015 Edward J. Los - Add gsc_bin_index for convience of preparing abbreviated sorts
 *                              Add peakDradRMS2 for the dradRMS2 values added in quadrature
 *                              Add 'g' qualifier for a table of limited gsc bin indices to search
 * Oct 29, 2015 Edward J. Los - Add peakDradRMS3 for the drad rms of points in the selected 90 day window relative to peakRA and peakDec
 *                              Do not generate the timestamp from the current time - have the user pick the string.
 * Dec 22, 2015 Edward J. Los - Merge update_summary2 search with this search for transient candidates.
 * Dec 28, 2015 Edward J. Los - Fix memory leak
 * Dec  2, 2017 Edward J. Los - Implement MIN_TC_DECLINATION, MAX_TC_DECLINATION, MIN_TC_RIGHTASCENSION, MAX_TC_RIGHTASCENSION to limit the search
 *                              Correct a bug for stars that are not in the catalog by giving them Stdmag = 99.0 and color = 99.0
 * Jan 30, 2017 Edward J. Los - Add limiting magnitude support   
 * Mar  8, 2017 Edward J. Los   Add peakLimitingYears and peakLimitingPoints
 * Mar 28, 2017 Edward J. Los   Add peakEvaluation bitmap to indicate reasons for rejection of stars
 * Mar  5, 2018 Edward J. Los   Add peakNoDefectMag, the brightest point in a transient that is hot a defect
 * Apr  2, 2018 Edward J. Los - Add peakLongOutburst to the id table
 * May 29, 2018 Edward J. Los - Support gaia
 * Aug 13, 2018 Edward J. Los   Define enableRematch to optimize location of transients (Redefine -O) to support this function
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 * Jan  4, 2019 Edward J. Los   Add verbose to LocateNoneImages and LoadNoneImages to display the files searched.
 * Aug 12, 2019 Edward J. Los   Split the last processor load into two to account for the very high number of plates at the North Equatorial Pole.  
 * Dec  9, 2020 Edward J. Los   Add InitFlareCandidateList add '-w' to look for flare candidates only
 * Jan 20, 2021 Edward J. Los   Add gsc_bin_index to FLARECANDIDATE 
 *
 *   search_none  -c 1 -i   -o /dasch/Pipeline/ingest/search_none.log -d /dasch/Pipeline/ingest/search_none.db -p /dasch/Pipeline/ingest/plate_outlier.db -u /dasch/Pipeline/ingest/id_unmatched_gsc.db
 *   echo "search_none -c 1 -i  -o /dasch/Pipeline/ingest/search_none.log -d /dasch/Pipeline/ingest/search_none.db  -p /dasch/Pipeline/ingest/plate_outlier.db -u /dasch/Pipeline/ingest/id_unmatched_gsc.db" | at now
 *
 *    
 *   Kepler field starts at GSC bin about 141505545
 *   search_none  -q kepler -c 1 -i   -o /dasch/Pipeline/ingest/search_none_kepler.log -d /dasch/Pipeline/ingest/search_none_kepler.db  -p /dasch/Pipeline/ingest/plate_outlier_kepler.db  -u /dasch/Pipeline/ingest/id_unmatched_kepler.db
 *   echo "search_none -q kepler -c 1 -i  -o /dasch/Pipeline/ingest/search_none_kepler.log -d /dasch/Pipeline/ingest/search_none_kepler.db  -p /dasch/Pipeline/ingest/plate_outlier_kepler.db  -u /dasch/Pipeline/ingest/id_unmatched_kepler.db" | at now
 *
 *  apass
 *   search_none  -q apass -c 1 -i   -o /dasch/Pipeline/ingest/search_none_apass.log -d /dasch/Pipeline/ingest/search_none_apass.db  -p /dasch/Pipeline/ingest/plate_outlier_apass.db  -u /dasch/Pipeline/ingest/id_unmatched_apass.db
 *   echo "search_none -q apass -c 1 -i  -o /dasch/Pipeline/ingest/search_none_apass.log -d /dasch/Pipeline/ingest/search_none_apass.db  -p /dasch/Pipeline/ingest/plate_outlier_apass.db  -u /dasch/Pipeline/ingest/id_unmatched_apass.db" | at now

 *
 *  experimental
 *   search_none  -q experimental -c 1 -i   -o /dasch/Pipeline/ingest/search_none_experimental.log -d /dasch/Pipeline/ingest/search_none_experimental.db  -p /dasch/Pipeline/ingest/plate_outlier_experimental.db  -u /dasch/Pipeline/ingest/id_unmatched_experimental.db
 *   echo "search_none -q experimental -c 1 -i  -o /dasch/Pipeline/ingest/search_none_experimental.log -d /dasch/Pipeline/ingest/search_none_experimental.db  -p /dasch/Pipeline/ingest/plate_outlier_experimental.db  -u /dasch/Pipeline/ingest/id_unmatched_experimental.db" | at now

 *
 * Performance on dasch3 on Feb 21, 2014 while rewriting binary files: network read rate 3.3 MB/sec, write rate 2.4 MB/sec.  Peak rate of 95.8 MB/sec on read
 */


#include <math.h>
#include <errno.h>
#include <mpi.h>
#include <assert.h>
#include "table.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "libwcs/wcscat.h"
#include "time.h"
#include <sys/time.h>
#include "pipelineutils.h"
#include "mysql.h"
#include "photometryutils.h"
#include "galaxyutils.h"
#include "searchgsc.h"
#include <sys/stat.h>
#include "daschunistd.h"

#define MAX_BUFFER 1024

#define VECTOR_ALLOC_INCREMENT 500
#if 0
#define MAGNITUDE_LIMIT 2000
#else
#define MAGNITUDE_LIMIT 50000
#endif
#define RESET_DELTA_ENTRY -1
#define RESET_ALL_ENTRY -2
/* #define OUTPUT_COUNTS 1 */
#define LOS_DEBUG 1
/* #define DEBUG_PARENT 1 */

#if 1 /* V404 Cyg */
#define LOS_GSC_BIN 131570315
#define LOS_REF_NUMBER 11033222260707L
#define LOS_PLATENUMBER 63778
#define LOS_SERIESID 6
#define LOS_NUMBER  328580
#else
/* S220000212462 plate j03056 */
#define LOS_GSC_BIN 81206654
#define LOS_REF_NUMBER 12220000212462L
#define LOS_PLATENUMBER 3056
#define LOS_SERIESID 21 
#define LOS_DEBUG_NUMBER 1  
#endif

#if 0
#define STAR_LIMIT     10000
#define SUMMARY_LIMIT  100
#else
#define STAR_LIMIT     50000
#define SUMMARY_LIMIT  50000
#endif


/* The following for N03102131317 at 311.22600006 7.72723982 raPM = 4.6 decPM = -9.7
                     median measured 311.229      7.72543
                     nearby object   311.234      7.71949
   
   It shows the malmquist effect 
   unknown object                    311.228      7.7252
   RA -.001 degrees in 70 years raPM = -51  mas/yr  (unmatched  slightly positive when matched)
    DEC +0.0004 degrees in 70 years  = +20  mas/yr  (unmatched, slightly negative when matched)
   unknown object         
*/
#if 0
#define START_GSC_BIN_INDEX  (131570315)        /* V404 Cyg N033222260707   bin 14682532  */
#define END_GSC_BIN_INDEX    (131570315 -1+2)
#endif

#if 0
#define START_GSC_BIN_INDEX  (81206654)        /*  S220000212462_2020_09_14.png gsc2.3.2 */
#define END_GSC_BIN_INDEX    (81206654 +1)
#endif



#if 0                         
#define START_GSC_BIN_INDEX  (123673965)        /* APASS_J121318.0+862324 Stdmag check */
#define END_GSC_BIN_INDEX    (123673965+1)
#endif

#if 0                         
#define START_GSC_BIN_INDEX  (42714112)        /* V840 Oph  */
#define END_GSC_BIN_INDEX    (42714112+1024)
#endif


#if 0
#define START_GSC_BIN_INDEX  (6309888)        /* exception Oct 15, 2018  */
#define END_GSC_BIN_INDEX    (6309888+40960)
#endif

#if 0
#define START_GSC_BIN_INDEX  ( 95048188 )        /*MAXI_J1820+07 */
#define END_GSC_BIN_INDEX    ( 95048188 +40960)
#endif


#if 0
#define START_GSC_BIN_INDEX  (57387008)        /* DASCH_J182400.3-184213   */
#define END_GSC_BIN_INDEX    (57408512+1024)
#endif

#if 0
#define START_GSC_BIN_INDEX  (113774592)        /* APASS_J075731.1+201735 Josh memo of  1 Apr 2018 20:30:06  'REFNumber == 407573111201735 gsc_bin_index = 113774987 */
#define END_GSC_BIN_INDEX    (113774592+1024)
#endif

#if 0
#define START_GSC_BIN_INDEX  (73902080)        /* S Vir REFNumber = 413330012071141, REF = "APASS_J133300.1-071141" gsc_bin_index = 73902853 */
#define END_GSC_BIN_INDEX    (73902080+1024)
#endif

#if 0
#define START_GSC_BIN_INDEX  (113125376)
#define END_GSC_BIN_INDEX    (113125376+1024)
#endif


#if 0
#define START_GSC_BIN_INDEX  (147304448)  /* KV UMa */
#define END_GSC_BIN_INDEX    (147304448+1024)
#endif


#if 0
#define START_GSC_BIN_INDEX  (65062912)  /* DASCH_J122716.4-13174 */
#define END_GSC_BIN_INDEX    (65062912+1024)
#endif


#if 0
#define START_GSC_BIN_INDEX  (99138560)  /* APASS_J085000.9+095947 SV Cnc */
#define END_GSC_BIN_INDEX    (99138560+1024)
#endif


#if 0
#define START_GSC_BIN_INDEX  (105669632)  /*  DASCH_J115014.1+143141, the infamous torus */
#define END_GSC_BIN_INDEX    (105699632+1024)
#endif



#if 0
#define START_GSC_BIN_INDEX  (6326272)  /*  DASCH_J053216.3-674028 - best candidate available on development machine */
#define END_GSC_BIN_INDEX    (6326272+1024)
#endif

#if 0
#define START_GSC_BIN_INDEX  (17603584)  /* S130212057005 bin 17603922 */
#define END_GSC_BIN_INDEX    (17603584+1024)
#endif

#if 0
#define START_GSC_BIN_INDEX  (111233024)  /* APASS_J090439.0+182751 bin	111233110 */
#define END_GSC_BIN_INDEX    (111233024+1024)
#endif


#if 0
#define START_GSC_BIN_INDEX  (140803072)  /* DASCH_J174033.5+414755  */
#define END_GSC_BIN_INDEX    (140803072+1024)
#endif


#define KEEPALIVE_CHECK_MODULUS (50*MAG_FILE_MODULUS)
#define FILE_PRINT_MODULUS (5000*MAG_FILE_MODULUS)
#define MIN_OUTLIER_SEARCH 20 /* Minimum number of stars for an outlier search */

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
extern FLAGSENTRY AflagsTable[];
extern int AflagsTableSize; 
extern FLAGSENTRY BflagsTable[];
extern int BflagsTableSize;
extern QUALITYBIT qualityMasks[];
extern int qualityTableSize;
char catalogname[MAX_BUFFER];
char indexname[MAX_BUFFER];


/* 
 *  Search the spatial bin for good matched stars with good points.  If one is found with more than 10 good points, then look for outliers
 *
 */
#define VECTOR_MAGCAL_MAGDEP1 0
#define VECTOR_MAGCAL_MAGDEP2 1  /* Necessary because the median calculation sorts the vectors */
#define MAX_VECTOR          2


#define VECTOR_PROCESSED 0
#define VECTOR_INDEX     1
#define MAX_IVECTOR      2
typedef struct _plateCount {
  int measurementCount;
  int matchCount;
  int lightcurveCount;
  int outlierCount;
} PLATECOUNT,*PPLATECOUNT;
typedef struct _bitcounts {
  long long count;
  long long outliers;
  long long plot2;  /* Images left over after filtering with FILTER_AMASK_PLOT2  */
  long long none_candidate; /* Images left over after filtering with FILTER_AMASK_NONE_CANDIDATE */
} BITCOUNTS,*PBITCOUNTS;


typedef struct _plateoutlier {
  char Plate[MAX_PLATE_NAME];
  int seriesId;
  int plateNumber;
  int measurementCount;
  int matchCount;
  double matchRatio;
  int lightcurveCount;
  int outlierCount;
  double outlierRatio;
} PLATEOUTLIER,*PPLATEOUTLIER;
typedef struct _vectorcounterblock {
  int vectorAlloc;   /* Size of the vectors */
  int outlierCount;
  int lightcurvePointCount;
  int lightcurveCount;
  int qualityRejectCount; /* Why is this here? duplicated in pFileCommon */
  int pointsRead2;
  int pointsRead;
  int pointsRejected;
  int insufficientPoints;
  int highRMSPoints;
  long long totalStars;
  long long totalOutliers;
  long long totalPlot2;
  long long totalNoneCandidates;
} VECTORCOUNTERBLOCK, *PVECTORCOUNTERBLOCK;

typedef struct _processmatchedcommon {
  time_t oldTime;
  int versionId;
} MATCHEDCOMMON,*PMATCHEDCOMMON;

typedef struct _vectorstorage {
  VECTORCOUNTERBLOCK vectorCounterBlock;
  BITCOUNTS AFLAGS_COUNTS[32];
  BITCOUNTS BFLAGS_COUNTS[32];
  BITCOUNTS QUALITY_COUNTS[32];
  double     *vector[MAX_VECTOR];
  int        *ivector[MAX_IVECTOR];
  PPLATECOUNT plateCountTable[MAX_SERIES+1];

} VECTORSTORE,*PVECTORSTORE;
#define STATE_INIT                  0
#define STATE_STATUS_SENT           1
#define STATE_STATUS_RECEIVED       2
#define STATE_WORK_SENT             3
#define STATE_WORK_RECEIVED         4
#define STATE_DIE_SENT              5
#define STATE_DIE_RECEIVED          6
#define STATE_ERROR                 7
#define STATE_MAX                   8

static char *stateString[STATE_MAX] = {
  (char *)"INIT           ",
  (char *)"STATUS_SENT    ",
  (char *)"STATUS_RECEIVED",
  (char *)"WORK_SENT      ",
  (char *)"WORK_RECEIVED  ",
  (char *)"DIE_SENT       ",
  (char *)"DIE_RECEIVED   ",
  (char *)"ERROR          ",
};


/* Work queue */
typedef struct _workqueue {
  int queueIndex;       /* Index into the work queue */
  int entryState;       /* State of this entry */
  int skipFlag;         /* Skip this item if nonzero */
  int workResult;       /* Status of the work */
  int startGscBinIndex; /* start GSC index */
  int endGscBinIndex;   /* Ending GSC index */
  char consoleLogname[MAX_BUFFER];
  char outLogname[MAX_BUFFER];
  char databaseLogname[MAX_BUFFER];
  char propermotionLogname[MAX_BUFFER];
  char plateLogname[MAX_BUFFER];
  char idLogname[MAX_BUFFER];
  char astrometryLogname[MAX_BUFFER];
  /* The response follows */
  int errorFlag;
  off_t consoleFilesize;   /* Console log */
  off_t outFilesize;   /* output log */
  off_t databaseFilesize; /* image database */
  off_t propermotionFilesize;  /* proper motion file */
  off_t plateFilesize;   /* plate quality database */
  off_t idFilesize;  /* id file database */
  off_t astrometryFilesize;  /* drad database */
  int seconds;       /* execution time */
  /* Additional statistics returned by children */
  int histogramTable[MAX_NONE_HISTOGRAM];
  int totalStars;
  int recoveredStars;
  int badModulusBin;
  int binCount;
  int groupNumber;
  int noGscEntryCount;
  int staleVersionIdCount;
  int curGscBinIndex;
  int totalProcessedMagnitudes;
  int totalNoneMagnitudes;
  int totalBlendNomatchMagnitudes;
  int totalMagnitudes;
  int queryCount;
  int unmatchedStars;
  int goodUnmatchedStars;
  FILECOUNTERBLOCK counterBlock; /* pFileCommon counters */
  VECTORCOUNTERBLOCK vectorCounterBlock; /* pVectorBlock counters */
  MATCHEDCOMMON matchedCommon; /* Variables used by ProcessMatchedImages() */
  int catalogDistTotal[MAX_TC_DIST];  /* Histogram of the distances of DASCH objects from nearby catalog objects */
  int catalogDistSimilarMag[MAX_TC_DIST]; /* Histogram of the distances of DASCH objects from nearby catalog objects within MAX_TC_NEARBY_MAG magnitudes */
  int catalogDistDifferentMag[MAX_TC_DIST]; /* Histogram of the distances of DASCH objects from nearby catalog objects differing by more than MAX_TC_NEARBY_MAG magnitudes */
} WORKQUEUE,*PWORKQUEUE;
#define MAX_MESSAGE_SIZE sizeof(WORKQUEUE)

typedef struct _nodestats {
  char nodename[MPI_MAX_PROCESSOR_NAME];
  /* Overall Cumulative stats */
  int taskscompleted;
  int totalSeconds;
  int elapsedSeconds;
  /* Stats since the last STATUSTAG message */
  int newTaskscompleted;
  int newTotalSeconds;
  int newElapsedSeconds;
  long long rxrate; /* mB/sec received */
  long long txrate; /* mB/sec sent */
  long long rxdiskrate; /* mB/sec received */
  long long txdiskrate; /* mB/sec sent */
  
} NODESTATS,*PNODESTATS;

void *kd_create(int k);

int master(FILE *consoleHandle,PFILECOMMON pFileCommon,PWORKQUEUE pWorkQueue,FILE *outHandle,FILE *databaseHandle,FILE *plateHandle,FILE *idHandle,FILE *astrometryHandle,char *propermotionfile,char *qualifier);
int slave(FILE *consoleHandle,PFILECOMMON pFileCommon,char *qualifier);
int SearchProcess(FILE *consoleHandle,PFILECOMMON pFileCommon,PWORKQUEUE pWorkEntry,char *qualifier);
void ProcessMatchedImages(FILE *consoleHandle,PFILECOMMON pFileCommon,PWORKQUEUE pWorkEntry,PPARAMETERVECTORSTORE pVectorStore,PSTARENTRY pStarTable,FILE *outHandle,FILE *idHandle,File catalogHandle,File indexHandle,int gsc_bin_index);

#define WORKTAG 1
#define RESPONSETAG 2
#define STATUSTAG 3
#define DIETAG 4
#define MAXTAG 5

static char *tagString[STATE_MAX] = {
  (char *)"Unknown    ",
  (char *)"WORKTAG    ",
  (char *)"RESPONSETAG",
  (char *)"STATUSTAG  ",
  (char *)"DIETAG     "
};

int binGscIndex[BIN_HISTOGRAM_SIZE];
long long binCount[BIN_HISTOGRAM_SIZE];


void ReallocVectorStore(FILE *consoleHandle,PVECTORSTORE pVectorStore)
{
  int index;
  double *tmpVector;
  int *tmpivector;
  for (index = 0; index < MAX_VECTOR; index++) {
    tmpVector = realloc(pVectorStore->vector[index],pVectorStore->vectorCounterBlock.vectorAlloc * sizeof(double));
    if (tmpVector == NULL) {
      fprintf(consoleHandle,"ERROR: Failed to reallocate vector %d of size %d\n",index,pVectorStore->vectorCounterBlock.vectorAlloc*sizeof(double));
      exit(-1);
    }
    pVectorStore->vector[index] = tmpVector;

  }

  for (index = 0; index < MAX_IVECTOR; index++) {
    tmpivector = realloc(pVectorStore->ivector[index],pVectorStore->vectorCounterBlock.vectorAlloc * sizeof(int));
    if (tmpivector == NULL) {
      fprintf(consoleHandle,"ERROR: Failed to reallocate ivector %d of size %d\n",index,pVectorStore->vectorCounterBlock.vectorAlloc*sizeof(double));
      exit(-1);
    }
    pVectorStore->ivector[index] = tmpivector;

  }

}

void FreeVectorStore(PVECTORSTORE pVectorStore)
{
  int index;
  double *vector;
  int * ivector;
  PPLATECOUNT  pPlateCount;
 
  for (index = 0; index < MAX_VECTOR; index++) {
    vector = pVectorStore->vector[index];
    if (vector != NULL) {
      free(vector);
    }
  }
  for (index = 0; index < MAX_IVECTOR; index++) {
    ivector = pVectorStore->ivector[index];
    if (ivector != NULL) {
      free(ivector);
    }
  }
  for (index = 0; index <= MAX_SERIES; index++) {
    pPlateCount = pVectorStore->plateCountTable[index];
    if (pPlateCount != NULL) {
      free(pPlateCount);
    }
  }

  memset(pVectorStore,0,sizeof(VECTORSTORE));
}

/* Search for plates with many outliers.  From Sumin Tang's memo of 2/03/11 4:04 PM choose lightcurves with rms_local < 0.2 and use nburst3  (number of good lightcurve points 0.4 mag above the median value) and ndip3 (number of good lightcurve points 0.4 mag below the median value.  The result is written in the form of files named "plate_outlier*.db"
 *
 * On November 6, 2013, we also search the DRAD range for catalog stars that have good astrometry.
 */
void ProcessStars(FILE *consoleHandle,PFILECOMMON pFileCommon,PVECTORSTORE pVectorStore,PPHOTSTARIMAGE pMagnitudeTable,int curMagnitudes,int gsc_bin_index, FILE *astrometryHandle)
{
  PPHOTSTARIMAGE pCurStarImage;
  PPHOTSTARIMAGE pNewStarImage;
  PFILESTARIMAGE pCurFileStarImage;
  PFILESTARIMAGE pNewFileStarImage;
  
  int index1;
  int index2;
  int index3;
  double full_med;
  double full_rms;
  double offset;
  double derrreal = 0.0;
  int vectorSize;
  int qualityMask;
  int quality;
  PPLATECOUNT plateCountTable;
  int bit;
  int *pColorflag;
  double *pColorterm;

  pVectorStore->vectorCounterBlock.pointsRead += curMagnitudes;

  memset(pVectorStore->ivector[VECTOR_PROCESSED],0,curMagnitudes*sizeof(int));
  memset(pVectorStore->ivector[VECTOR_INDEX],0,curMagnitudes*sizeof(int));
  for (index1 = 0; index1 < curMagnitudes; index1++) { 
    pCurStarImage = &pMagnitudeTable[index1];
    pCurFileStarImage = pCurStarImage->pFileStarImage;


		/* Handle the drad histogram here */
		if (astrometryHandle != NULL) {
			if ((pCurFileStarImage->REFNumber != 0) &&
					((pCurFileStarImage->AFLAGS & FILTER_AMASK_LOWDRAD) == 0)) {
				if ((pCurFileStarImage->ra_2 < 900) && (pCurFileStarImage->dec_2 < 91)) {
					int magindex = pCurFileStarImage->magcal_magdep;
					if ((magindex > 0) && (magindex < MAX_DRAD_HISTOGRAM_MAG)) {
						double factor = cos(DEGREES_TO_RAD*pCurFileStarImage->dec);
						int tempdra = (3600.0*(factor*(pCurFileStarImage->ra - pCurFileStarImage->ra_2)))+ MAX_DRAD_RANGE;
						if (tempdra < 0) {
							tempdra = 0;
						}
						if (tempdra >= (2*(MAX_DRAD_RANGE)+1)) {
							tempdra = 2*MAX_DRAD_RANGE;
						}
						pFileCommon->draHistogram[magindex][tempdra]++;
						pFileCommon->draHistogram[0][tempdra]++;
						int tempddec = (3600.0*(pCurFileStarImage->dec - pCurFileStarImage->dec_2)) + MAX_DRAD_RANGE;
						if (tempddec < 0) {
							tempddec = 0;
						}
						if (tempddec >= (2*(MAX_DRAD_RANGE)+1)) {
							tempddec = 2*MAX_DRAD_RANGE;
						}
						pFileCommon->ddecHistogram[magindex][tempddec]++;
						pFileCommon->ddecHistogram[0][tempddec]++;
						
					}
				}

			}


		}


#if 0
    if ((pCurFileStarImage->REFNumber == 412140552130929L) &&
        (pCurFileStarImage->plateNumber == 4907) && 
        (pCurFileStarImage->seriesId == 30)) {
      fprintf(consoleHandle,"At %lld\n",pCurFileStarImage->REFNumber);
    }
#endif


    if (pVectorStore->ivector[VECTOR_PROCESSED][index1] != 0) {
      continue;
    }
    if ((pCurFileStarImage->AFLAGS != 0) ||
        (pCurFileStarImage->gsc_bin_index != gsc_bin_index) ||
        (pCurFileStarImage->REFNumber == 0)) {
      pVectorStore->ivector[VECTOR_PROCESSED][index1] = 1;
      pVectorStore->vectorCounterBlock.pointsRejected++;
      continue;
    }
#if 0
    if ((pCurFileStarImage->magcal_local_error > 90.0) ||
        (pCurFileStarImage->magcal_iso_rms > 90.0)) {
#if 0
      if (pCurFileStarImage->REFNumber == 413131971115824L) {
        fprintf(consoleHandle,"At %lld\n",pCurFileStarImage->REFNumber);
      }
#endif
      fprintf(consoleHandle,"ERROR: sanity check failed (1) for AFLAGS %d, magcal_local_error %f, or magcal_iso_rms %f gsc_bin_index %d %s%05d REFNumber %lld\n",
             pCurFileStarImage->AFLAGS,
             pCurFileStarImage->magcal_local_error,
             pCurFileStarImage->magcal_iso_rms,
             pCurFileStarImage->gsc_bin_index,
             pCurStarImage->series,
             pCurFileStarImage->plateNumber,
             pCurFileStarImage->REFNumber);
      continue;
    }
#endif
    
    


    /* A good star.  Search for other good stars */
    pVectorStore->ivector[VECTOR_PROCESSED][index1] = 1;
    vectorSize = 0;
    pVectorStore->vector[VECTOR_MAGCAL_MAGDEP1][vectorSize] = pCurFileStarImage->magcal_magdep;
    pVectorStore->vector[VECTOR_MAGCAL_MAGDEP2][vectorSize] = pCurFileStarImage->magcal_magdep;
    pVectorStore->ivector[VECTOR_INDEX][vectorSize] = index1;
    vectorSize++;
    for (index2 = index1+1; index2 < curMagnitudes; index2++) {
      pNewStarImage = &pMagnitudeTable[index2];
      pNewFileStarImage = pNewStarImage->pFileStarImage;
      if (pNewFileStarImage->gsc_bin_index != pCurFileStarImage->gsc_bin_index) {
        break;
      }
      if ((pNewFileStarImage->AFLAGS != 0) ||
          (pNewFileStarImage->gsc_bin_index != gsc_bin_index) ||
          (pNewFileStarImage->REFNumber == 0)) {
        pVectorStore->ivector[VECTOR_PROCESSED][index1] = 1;
        pVectorStore->vectorCounterBlock.pointsRejected++;
        continue;
      }

      if (pNewFileStarImage->REFNumber == pCurFileStarImage->REFNumber) {
        pVectorStore->ivector[VECTOR_PROCESSED][index2] = 1;
#if 1
        if ((pNewFileStarImage->magcal_local_error > 90.0) ||
            (pNewFileStarImage->magcal_iso_rms > 90.0)) {
          fprintf(consoleHandle,"ERROR: sanity check failed (2) for AFLAGS %d, magcal_local_error %f, or magcal_iso_rms %f gsc_bin_index %d %s%05d REFNumber %lld\n",
                 pNewFileStarImage->AFLAGS,
                 pNewFileStarImage->magcal_local_error,
                 pNewFileStarImage->magcal_iso_rms,
                 pNewFileStarImage->gsc_bin_index,
                 pNewStarImage->series,
                 pNewFileStarImage->plateNumber,
                 pNewFileStarImage->REFNumber);
          continue;
        }
#endif


        /* We have a good one */
        
        
        pVectorStore->vector[VECTOR_MAGCAL_MAGDEP1][vectorSize] = pNewFileStarImage->magcal_magdep;
        pVectorStore->vector[VECTOR_MAGCAL_MAGDEP2][vectorSize] = pNewFileStarImage->magcal_magdep;
        pVectorStore->ivector[VECTOR_INDEX][vectorSize] = index2;
        vectorSize++;
        
      }
    }
    if (vectorSize < MIN_OUTLIER_SEARCH) {
      pVectorStore->vectorCounterBlock.insufficientPoints += vectorSize;
      continue;
    }
    if (CalcMedianAndRMS(vectorSize,0,pVectorStore->vector[VECTOR_MAGCAL_MAGDEP2],&full_med,&full_rms,0,3.0,0) == 0) {
      fprintf(consoleHandle,"ERROR: CalcMedianAndRMS failed\n");
      exit(-1);
      continue;        
    }
    if (full_rms > 0.14) { /* Start with 0.2, change to 0.14 on 02/17/2011) */
      pVectorStore->vectorCounterBlock.highRMSPoints += vectorSize;
      continue;
    }
    pVectorStore->vectorCounterBlock.lightcurveCount++;

    for (index3 = 0; index3 < vectorSize; index3++) {
      index2 = pVectorStore->ivector[VECTOR_INDEX][index3];
      pNewStarImage = &pMagnitudeTable[index2];
      pNewFileStarImage = pNewStarImage->pFileStarImage;
      if ((pNewFileStarImage->AFLAGS != 0) ||
          (pNewFileStarImage->REFNumber != pCurFileStarImage->REFNumber) ||
          (pNewFileStarImage->magcal_local_error > 90.0) ||
          (pNewFileStarImage->magcal_iso_rms > 90.0)) {
        fprintf(consoleHandle,"ERROR: sanity check failed (3) for AFLAGS %d, magcal_local_error %f, or magcal_iso_rms %f gsc_bin_index %d REFNumber %d %lld\n",
               pNewFileStarImage->AFLAGS,
               pNewFileStarImage->magcal_local_error,
               pNewFileStarImage->magcal_iso_rms,
               pNewFileStarImage->gsc_bin_index,
               (pNewFileStarImage->REFNumber != pCurFileStarImage->REFNumber),pNewFileStarImage->REFNumber);
        exit(-1);
      }

      pVectorStore->vectorCounterBlock.lightcurvePointCount++;
      if ((pNewFileStarImage->seriesId <=0) ||
          (pNewFileStarImage->seriesId > MAX_SERIES)) {
        fprintf(consoleHandle,"ERROR: invalid seriesId %d\n",pNewFileStarImage->seriesId);
        exit(-1);
      }
      plateCountTable = pVectorStore->plateCountTable[pNewFileStarImage->seriesId];
      if (plateCountTable == NULL) {
        plateCountTable = (PPLATECOUNT)calloc(pFileCommon->maxPlateNumber[pNewFileStarImage->seriesId]+1,sizeof(PLATECOUNT));
        if (plateCountTable == NULL) {
          fprintf(consoleHandle,"ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[pNewFileStarImage->seriesId]+1);
          exit(-1);

        }
        pVectorStore->plateCountTable[pNewFileStarImage->seriesId] = plateCountTable;
      }
      if ((pNewFileStarImage->plateNumber <= 0) ||
          (pNewFileStarImage->plateNumber > pFileCommon->maxPlateNumber[pNewFileStarImage->seriesId])) {
        fprintf(consoleHandle,"ERROR: illegal plate number %d\n",pNewFileStarImage->plateNumber);
        exit(-1);
      }
      plateCountTable[pNewFileStarImage->plateNumber].lightcurveCount++;

#if 0
#if 0
      derrreal = full_rms;
#else
      derrreal =  sqrt(sqr(pNewFileStarImage->magcal_local_error) + sqr(pNewFileStarImage->magcal_iso_rms))/2.0;
#endif
#endif
      offset = fabs(pVectorStore->vector[VECTOR_MAGCAL_MAGDEP1][index3]-full_med) - derrreal;
      
      if (offset > 0.4) {  
        pVectorStore->vectorCounterBlock.outlierCount++;
        if ((pNewFileStarImage->plateNumber <= 0) ||
            (pNewFileStarImage->plateNumber > pFileCommon->maxPlateNumber[pNewFileStarImage->seriesId])) {
          fprintf(consoleHandle,"ERROR: illegal plate number %d\n",pNewFileStarImage->plateNumber);
          exit(-1);
        }
        plateCountTable[pNewFileStarImage->plateNumber].outlierCount++;

      }

    }
    /* Here we consider all points that match our object of interest and tally up AFLAGS and BFLAGS counts for outliers */
    for (index2 = index1+1; index2 < curMagnitudes; index2++) {
      pNewStarImage = &pMagnitudeTable[index2];
      pNewFileStarImage = pNewStarImage->pFileStarImage;
      if (pNewFileStarImage->gsc_bin_index != pCurFileStarImage->gsc_bin_index) {
        break;
      }
      qualityMask = GetPlateQualityMask(pFileCommon,pNewFileStarImage->seriesId,pNewFileStarImage->plateNumber,&quality,&pColorterm,&pColorflag,pNewFileStarImage->REFNumber);
      if ((qualityMask & (1 << pNewFileStarImage->spatial_bin)) != 0) {
        quality |= QUALITY_COLORTERM;
      }
      pVectorStore->vectorCounterBlock.totalStars++;
      /* Keep track of the AFLAG and BFLAG counts */
      for (bit = 0; bit < 32; bit++) {
        if (((1 << bit) & pNewFileStarImage->AFLAGS) != 0) {
          pVectorStore->AFLAGS_COUNTS[bit].count++;
        }
        if (((1 << bit) & pNewFileStarImage->BFLAGS) != 0) {
          pVectorStore->BFLAGS_COUNTS[bit].count++;
        }
        if (((1 << bit) & quality) != 0) {
          pVectorStore->QUALITY_COUNTS[bit].count++;
        }
      }
      if ((pNewFileStarImage->REFNumber == pCurFileStarImage->REFNumber) &&
          ((fabs(pNewFileStarImage->magcal_magdep-full_med) - derrreal) > 0.4)) {
        /* This is an outlier.  Keep track of the AFLAG and BFLAG counts */
        pVectorStore->vectorCounterBlock.totalOutliers++;
        for (bit = 0; bit < 32; bit++) {
          if (((1 << bit) & pNewFileStarImage->AFLAGS) != 0) {
            pVectorStore->AFLAGS_COUNTS[bit].outliers++;
          }
          if (((1 << bit) & pNewFileStarImage->BFLAGS) != 0) {
            pVectorStore->BFLAGS_COUNTS[bit].outliers++;
          }
          if (((1 << bit) & quality) != 0) {
            pVectorStore->QUALITY_COUNTS[bit].outliers++;
          }
        }
        if ((pNewFileStarImage->AFLAGS & (~FILTER_AMASK_PLOT2)) == 0) {
          pVectorStore->vectorCounterBlock.totalPlot2++;
          for (bit = 0; bit < 32; bit++) {
            if (((1 << bit) & pNewFileStarImage->AFLAGS) != 0) {
              pVectorStore->AFLAGS_COUNTS[bit].plot2++;
            }
            if (((1 << bit) & pNewFileStarImage->BFLAGS) != 0) {
              pVectorStore->BFLAGS_COUNTS[bit].plot2++;
            }
            if (((1 << bit) & quality) != 0) {
              pVectorStore->QUALITY_COUNTS[bit].plot2++;
            }
          }

        }
        if ((pNewFileStarImage->AFLAGS & FILTER_AMASK_NONE_CANDIDATE) == 0) {
          pVectorStore->vectorCounterBlock.totalNoneCandidates++;
          for (bit = 0; bit < 32; bit++) {
            if (((1 << bit) & pNewFileStarImage->AFLAGS) != 0) {
              pVectorStore->AFLAGS_COUNTS[bit].none_candidate++;
            }
            if (((1 << bit) & pNewFileStarImage->BFLAGS) != 0) {
              pVectorStore->BFLAGS_COUNTS[bit].none_candidate++;
            }
            if (((1 << bit) & quality) != 0) {
              pVectorStore->QUALITY_COUNTS[bit].none_candidate++;
            }
          }

        }
      }
      
    }
  }
} /* end of ProcessStars() */
#ifdef OUTPUT_COUNTS
void OutputCounts(char *header,PFLAGSENTRY flagsTable,PQUALITYBIT qualityTable,int tableSize,PBITCOUNTS counts) {
  double outRatio1;
  double outRatio2;
  double outRatio3;
  int bit;
  int index;
  PFLAGSENTRY pFlagEntry;
  PQUALITYBIT pQualityEntry;
  for (bit = 0; bit < 32; bit++) {
    if (counts[bit].count != 0) {
      outRatio1 = (1.0*counts[bit].outliers)/(1.0*counts[bit].count);
      outRatio2 = (1.0*counts[bit].plot2)/(1.0*counts[bit].count);
      outRatio3 = (1.0*counts[bit].none_candidate)/(1.0*counts[bit].count);
      fprintf(consoleHandle,"%s bit %2d total %10lld outliers %10lld ratio %.6f; plot2 %10lld ratio %.6f; none %10lld ratio %.6f ",header,bit,counts[bit].count,
             counts[bit].outliers,outRatio1,
             counts[bit].plot2,outRatio2,
             counts[bit].none_candidate,outRatio3);
             
      if (flagsTable != NULL) {
        for (index = 0; index < tableSize; index++) {
          pFlagEntry = &flagsTable[index];
          if (bit == pFlagEntry->flagbit) {
            fprintf(consoleHandle," %s",pFlagEntry->webname);
            break;
          }
        }
      }
      if (qualityTable != NULL) {
        for (index = 0; index < tableSize; index++) {
          pQualityEntry = &qualityTable[index];
          if ((1 << bit) == pQualityEntry->qualityMask) {
            fprintf(consoleHandle," %s",pQualityEntry->qualityDescr);
            break;
          }
        }


      }
      fprintf(consoleHandle,"\n");
    }
  }
}
#endif /* OUTPUT_COUNTS */

void	OutputRootCounts(FILE *consoleHandle,PFILECOMMON pFileCommon)
{
	int index;
	PFLAGSENTRY pFlagEntry;
	fprintf(consoleHandle,"AFLAGS zero count %d nonzero count %d\n",pFileCommon->counterBlock.AFLAGS_ZERO_COUNT,pFileCommon->counterBlock.AFLAGS_NONZERO_COUNT);
	for (index = 0; index < AflagsTableSize; index++) {
		pFlagEntry = &AflagsTable[index];
		if (pFileCommon->counterBlock.AFLAGS_MULT_COUNT[pFlagEntry->flagbit] != 0) {
			fprintf(consoleHandle,"AFLAGS bit %2d multiple count %10d unique count %10d %s\n",
						 pFlagEntry->flagbit,
						 pFileCommon->counterBlock.AFLAGS_MULT_COUNT[pFlagEntry->flagbit],
						 pFileCommon->counterBlock.AFLAGS_UNIQUE_COUNT[pFlagEntry->flagbit],
						 pFlagEntry->webname);
		}
	}
}

void ReadGscBinHistogram(PFILECOMMON pFileCommon,FILE *consoleHandle,int startGscBinIndex,int endGscBinIndex,char *binhistogramfile,int *binGscIndex,long long *binCount,int *pBinArraySize,long long *pBinTotal)
{
  FILE *binhistogramHandle = NULL;
  int lineNumber;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  int lineLen;
  int nvals;
  int binIndex;
  int binArraySize = 0;
  double lon;
  double lat;
  int binIncrement;
  *pBinArraySize = 0;
  /* Here we use the more accurate bin histogram to divide up the work */
  binhistogramHandle = fopen(binhistogramfile,"rt");
  if (binhistogramHandle == NULL) {
    fprintf(consoleHandle,"ERROR opening the bin histogram file %s\n",binhistogramfile);
    exit(-1);
  }
  lineNumber = 0;
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,binhistogramHandle);
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
    lineNumber++;
    if (lineNumber < 3) {
      continue;
    }
    nvals = sscanf(inBuffer,"%d\t%lld",&binGscIndex[binArraySize],&binCount[binArraySize]);
    if (nvals != 2) {
      fprintf(consoleHandle,"ERROR: bin histogram file nvals is %d for line %d %s\n",nvals,lineNumber,inBuffer);
      exit(-1);
    }
    binArraySize++;
    if (binArraySize > BIN_HISTOGRAM_SIZE) {
      fprintf(consoleHandle,"ERROR: bin histogram file size %d exceeded\n",binArraySize);
      exit(-1);
    } 
  }
  if (binArraySize < 2) {
    fprintf(consoleHandle,"ERROR insufficient bins in the bin histogram file\n");
    exit(-1);
  }
  binIncrement = binGscIndex[1] - binGscIndex[0];

  /* Now add everything up */
  for (binIndex = 0; binIndex < binArraySize; binIndex++) {
    if (((binGscIndex[binIndex]+binIncrement+MAG_FILE_MODULUS) < startGscBinIndex) ||
        (binGscIndex[binIndex] > endGscBinIndex)) {
      binCount[binIndex] = 0;
    } else {
      GetBinCenter(pGscBin,binGscIndex[binIndex],&lon,&lat,"search_none");
      if ((lat <  MIN_TC_DECLINATION) ||
          (lat > MAX_TC_DECLINATION) ||
          (lon < MIN_TC_RIGHTASCENSION) ||
          (lon > MAX_TC_RIGHTASCENSION)) {
        binCount[binIndex] = 0;
      }
      if (pFileCommon->minimumGalacticLatitude > -90.0) {
        wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
        if (lat < pFileCommon->minimumGalacticLatitude) {
          binCount[binIndex] = 0;
        }
      }
    }
    (*pBinTotal) += binCount[binIndex];
  }


  fclose(binhistogramHandle);
  *pBinArraySize = binArraySize;
}
int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  time_t deltaTime;
  time_t oldTime;
  int ignoreVersionId = 0;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char gscbinfile[MAX_BUFFER];
  char consolefile[MAX_BUFFER];
	char astrometryfile[MAX_BUFFER];
  char databasefile[MAX_BUFFER];
  char platefile[MAX_BUFFER];
	char propermotionfile[MAX_BUFFER];
	char binhistogramfile[MAX_BUFFER];
  char idfile[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int nvals;
  FILE *outHandle = NULL;
  FILE *databaseHandle = NULL;
	FILE *astrometryHandle = NULL;
  FILE *plateHandle = NULL;
  FILE *idHandle = NULL;
  char *dotPtr;
  int processedCount = 0;

  char *mysqlhost;
  char *username;
  char *password;
  char *mysqlphothost;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  char *photusername;
  char *photpassword;

  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  char timestr[100];
  struct tm *ptr;




  int goodMagnitudeCount;

 
  int curStars;
  int starIndex;
  int starIndex2;

  int res;
  char queryString[MAX_QUERY_STRING];
  int querySize;
  int oldQueryCount = 0;
  int oldTotalStars = 0;
  int oldGoodStars = 0;
  int oldBinCount = 0;
  int oldGscBinIndex;
  int newGscBinIndex;
#ifdef START_GSC_BIN_INDEX
  int curGscBinIndex = START_GSC_BIN_INDEX;
  int startGscBinIndex = START_GSC_BIN_INDEX;
  int endGscBinIndex = END_GSC_BIN_INDEX;
#else /* START_GSC_BIN_INDEX */
  int curGscBinIndex = 0;
  int startGscBinIndex = 0;
  int endGscBinIndex =  pGscBin->total_gsc_bins;
#endif /* START_GSC_BIN_INDEX */
  char* catalogdir;
  char* slashPtr;
  char galaxy_name[MAX_BUFFER];
  PGSCIMAGE pCurGscImage;
  int readItems;
  int gscIndex;

  int usecIndex;
  SetQueryCount(0);
  double borderDist;
  double maxBorderDist = 0;
  int workToDo;
  int catalogNumber = 0;
  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon = &fileCommon;
  int fileMagnitudeIndex;
  int qualityMask;
  int quality;
  PSKY2KSTAR pSky2KStar = NULL;

  int bit;
  PFLAGSENTRY pFlagEntry;
  PQUALITYBIT pQualityEntry;
  PGALAXYCOMMON pGalaxyCommon = &pFileCommon->galaxyCommon;
  char *charPtr;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;

  PWORKQUEUE pWorkEntry;
  PWORKQUEUE pWorkQueue;
  int workIndex;
  int gscIncrement;
  char *charPtr2;
  int verbose = 0;
  int enableTransientSearch = 0;
  int  minimumImageCount = -1;
  double minimumGalacticLatitude = DEFAULT_GAL_LATITUDE;
  int myrank = 0;
  FILE *consoleHandle = stdout;
  char *scriptdir;
  char consoleLogname[MAX_BUFFER];
  char timestr2[MAX_TIMESTR_SIZE];
  long long binTotal = 0;
  int binArraySize = 0;
  int beginGscIndex = 0;
  int endGscIndex = 0;
  int curBinIndex = 0;
  long long binCurrentCount = 0;
  long long binTargetCount = 0;
  int enableRematch = 0;

  void *candidate_ptree; /* kdtree for the previous candidates */
  int preallocFlag = 0;
  char interesting_name[] = "/dasch/Pipeline/flarecandidates.txt";
  int flareCandidatesOnly = 0;
  PFLARECANDIDATE pFlareCandidate = NULL;

  timestr2[0] = 0;

  memset(binGscIndex,0,sizeof(binGscIndex));
  memset(binCount,0,sizeof(binCount));

  assert(DASCH_MAX_PROCESSOR_NAME >= MPI_MAX_PROCESSOR_NAME);


  /* Initialize MPI */

  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  scriptdir = getenv("DASCH_SCRIPTS");
  if (scriptdir == NULL) {
    printf("ERROR: DASCH_SCRIPTS is not defined\n");
    exit(-1);
  }
 

  if (myrank > 0) {
    /* We need to redirect the console to a log file */ 
    sprintf(consoleLogname,"%s/search_none_console_%02d.log",scriptdir,myrank);
    consoleHandle = fopen(consoleLogname,"wt");
    if (consoleHandle == NULL) {
      printf("ERROR: failed to open console log %s\n",consoleLogname);
      exit(-1);
    }
    if (flareCandidatesOnly != 0) {
      printf("ERROR: line %5d flareCandidatesOnly not supported in multiprocessor mode\n",__LINE__);
      exit(-1);      
    }

    fprintf(consoleHandle,"CHILD Process starting with rank %d console handle 0x%x\n",myrank,consoleHandle);
  } else {
    consoleLogname[0] = 0;
    fprintf(consoleHandle,"MASTER Process starting with rank %d console handle 0x%x\n",myrank,consoleHandle);
  }

  if ((MIN_TC_DECLINATION != -90.0) ||
      (MAX_TC_DECLINATION != 90.0) ||
      (MIN_TC_RIGHTASCENSION != 0) ||
      (MAX_TC_RIGHTASCENSION != 360.0)) {
    printf("ERROR: limited search for dec %f to %f degrees and ra %f to %f degrees\n",MIN_TC_DECLINATION,MAX_TC_DECLINATION,MIN_TC_RIGHTASCENSION,MAX_TC_RIGHTASCENSION);
  }

#ifdef START_GSC_BIN_INDEX
  fprintf(consoleHandle,"ERROR: START_GSC_BIN_INDEX is set\n");
#endif /* START_GSC_BIN_INDEX */
#ifdef END_GSC_BIN_INDEX
	fprintf(consoleHandle,"ERROR: END_GSC_BIN_INDEX is set\n");
#endif /* END_GSC_BIN_INDEX */
  fflush(consoleHandle);

#if 0 /* Testing of DASCH reference designation support */
  {
    double dec;
    double ra;
    long long REFNumber;
    int result;
    int refType;
    char REF[MAX_REF];
    dec = 14.120293;
    ra = 127.660377;
    result = GetDASCHNumber(ra,dec,&REFNumber,1,REF_TYPE_DASCH);
    fprintf(consoleHandle,"result is %d for ra %f dec %f, REFNumber %lld\n",result,ra,dec,REFNumber);
    GetREF(REFNumber,REF,0,1);
    fprintf(consoleHandle,"REF is %s for REFNumber %lld\n",REF,REFNumber);
    result = GetREFNumber(REF,&REFNumber,&refType,1);
    fprintf(consoleHandle,"result %d for %s is REFNumber %lld\n",result,REF,REFNumber);
    exit(-1);
                            
  }


#endif
#ifdef LOS_DEBUG
  printf("DEBUG: line %4d LOS_GSC_BIN %d LOS_REF_NUMBER %lld\n",__LINE__,LOS_GSC_BIN,LOS_REF_NUMBER);
#endif
  candidate_ptree = (void*)kd_create(3);
  if (InitFlareCandidateList(pGscBin,interesting_name,&pFileCommon->interesting_count,&pFileCommon->candidate_table,candidate_ptree,&startTime,preallocFlag) != 0) {
#ifdef LOS_DEBUG
    printf("DEBUG: line %4d LOS_GSC_BIN %d LOS_REF_NUMBER %lld\n",__LINE__,LOS_GSC_BIN,LOS_REF_NUMBER);
#endif
    exit(-1);
  }
#ifdef LOS_DEBUG
  printf("DEBUG: line %4d LOS_GSC_BIN %d LOS_REF_NUMBER %lld\n",__LINE__,LOS_GSC_BIN,LOS_REF_NUMBER);
#endif


#if 0 /* Testing of FindAdjacentBins. Takes 1054 seconds with -O2 */
    {
    int gsc_bin;
    int bin_count;
    int bin_index;
    int max_bin_count = 0;
    int min_bin_count;
    int ave_bin_count;
    int bin_list[MAX_ADJACENT_BINS];
    time(&startTime);
    for (gsc_bin = 0; gsc_bin < pGscBin->total_gsc_bins; gsc_bin++) {
    FindAdjacentBins(pGscBin,gsc_bin,bin_list,&bin_count);
    if (bin_count > max_bin_count) {
    max_bin_count = bin_count;
  }
    if (gsc_bin == 0) {
    min_bin_count = bin_count;
  } else if (bin_count < min_bin_count) {
    min_bin_count = bin_count;
  }
    ave_bin_count += bin_count;
    if ((gsc_bin % 10000000) == 0) {
    time(&curTime);
    curTime -= startTime;
    fprintf(consoleHandle,"%10d of %10d bins min_bin_count %d ave_bin_count %f max_bin_count is %d in %d seconds\n",gsc_bin,pGscBin->total_gsc_bins,min_bin_count,(1.0*ave_bin_count)/(gsc_bin+1),max_bin_count,curTime);
    for (bin_index = 0; bin_index < bin_count; bin_index++) {
    fprintf(consoleHandle,"%d ",bin_list[bin_index]);
  }
    fprintf(consoleHandle,"\n");
  }
  }
    time(&curTime);
    curTime -= startTime;
    fprintf(consoleHandle,"min_bin_count %d ave_bin_count %f max_bin_count is %d in %d seconds\n",min_bin_count,(1.0*ave_bin_count)/gsc_bin,max_bin_count,curTime);
  }
    exit(-1);
#endif


    memset(pFileCommon->histogramTable,0,sizeof(pFileCommon->histogramTable));


    outfile[0] = 0;
    databasefile[0] = 0;
    qualifier[0] = 0;
    platefile[0] = 0;
    propermotionfile[0] = 0;
    astrometryfile[0] = 0;
    binhistogramfile[0] = 0;
    gscbinfile[0] = 0;

    /* Loop through the arguments */
    for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
    /* This must be the list of plates */
    errorFlag = 1;
    fprintf(consoleHandle,"ERROR: unqualified argument %s argc: %d\n",argstr,argc);
  } else {
    while ((cmdchar = *++argstr) != 0) {
    switch(cmdchar) {

  case 'o': /* output file name */
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    strncpy(outfile,*++argv,MAX_BUFFER-2);
    if (strlen(outfile) >= MAX_BUFFER-3) {
    fprintf(consoleHandle,"ERROR: MAX_BUFFER too small for %s\n",*argv);
  }
  }
    break;

  case 'O':
    enableRematch = 1;
    break;


  case 'g': /* limited bin list to search */
  case 'G':
    argc--;
    if (argc < 1) {
    printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    strncpy(gscbinfile,*++argv,MAX_BUFFER-2);
    if (strlen(gscbinfile) >= MAX_BUFFER-3) {
    printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
  }
  }
    break;
        
  case 'd': /* database file name */
  case 'D':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    strncpy(databasefile,*++argv,MAX_BUFFER-2);
    if (strlen(databasefile) >= MAX_BUFFER-3) {
    fprintf(consoleHandle,"ERROR: MAX_BUFFER too small for %s\n",*argv);
  }
  }
    break;

  case 'a': /* astrometry file name */
  case 'A':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    strncpy(astrometryfile,*++argv,MAX_BUFFER-2);
    if (strlen(astrometryfile) >= MAX_BUFFER-3) {
    fprintf(consoleHandle,"ERROR: MAX_BUFFER too small for %s\n",*argv);
  }
  }
    break;

  case 'p': /* plate summary file name */
  case 'P':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    strncpy(platefile,*++argv,MAX_BUFFER-2);
    if (strlen(platefile) >= MAX_BUFFER-3) {
    fprintf(consoleHandle,"ERROR: MAX_BUFFER too small for %s\n",*argv);
  }
  }
    break;

  case 'm': /* proper motion file name */
  case 'M':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    strncpy(propermotionfile,*++argv,MAX_BUFFER-2);
    if (strlen(propermotionfile) >= MAX_BUFFER-3) {
    fprintf(consoleHandle,"ERROR: MAX_BUFFER too small for %s\n",*argv);
  }
  }
    break;

  case 'b': /* Bin histogram file */
  case 'B':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    strncpy(binhistogramfile,*++argv,MAX_BUFFER-2);
    if (strlen(binhistogramfile) >= MAX_BUFFER-3) {
    fprintf(consoleHandle,"ERROR: MAX_BUFFER too small for %s\n",*argv);
  }
  }
    break;

  case 'u': /* id file name */
  case 'U':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    strncpy(idfile,*++argv,MAX_BUFFER-2);
    if (strlen(idfile) >= MAX_BUFFER-3) {
    fprintf(consoleHandle,"ERROR: MAX_BUFFER too small for %s\n",*argv);
  }
  }
    break;


  case 'i':  /* Ignore the versionId */
  case 'I':
    ignoreVersionId = 1;
    break;
 
  case 'c': /* Minimum image count */
  case 'C':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    nvals = sscanf(*++argv,"%d",&minimumImageCount);
    if (nvals != 1) {
    fprintf(consoleHandle,"ERROR: Unable to decode minimum image count %s\n",*argv);
    errorFlag = 1;
  }
  }
   

    break;


  case 'l': /* Minimum galactic latitude */
  case 'L':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    nvals = sscanf(*++argv,"%lf",&minimumGalacticLatitude);
    if (nvals != 1) {
    fprintf(consoleHandle,"ERROR: Unable to decode minimum galactic latitude %s\n",*argv);
    errorFlag = 1;
  }
  }
   

    break;



  case 's': /* Starting GSC bin index */
  case 'S':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    nvals = sscanf(*++argv,"%d",&curGscBinIndex);
    if (nvals != 1) {
    fprintf(consoleHandle,"ERROR: Unable to decode initial gsc bin index %s\n",*argv);
    errorFlag = 1;
  }
  }
   

    break;


  case 'v': /* verbose */
  case 'V':
    verbose = 1;
    break;

  case 'w': /* Handle flare candidates only */
  case 'W':
    flareCandidatesOnly = 1;
    break;

  case 'T': /* Enable transient flare search */
  case 't': 
    enableTransientSearch = 1;
    if ((argc > 1) && (argv[1][0] != '-')) {
    argc--;
    strncpy(timestr2,*++argv,MAX_TIMESTR_SIZE-2);
    timestr2[MAX_TIMESTR_SIZE-1] = 0;
    if (strlen(timestr2) >= MAX_TIMESTR_SIZE-3) {
    printf("ERROR: MAX_TIMESTR_SIZE too small for %s\n",*argv);
    errorFlag = 1;
  }
  }
    break;

  case 'q': /* Catalog and file name qualifier */
  case 'Q':
    argc--;
    if (argc < 1) {
    fprintf(consoleHandle,"ERROR: Insufficient arguments for -%c\n",cmdchar);
    errorFlag = 1;
  } else {
    catalogNumber = GetCatalogNumber(*++argv);
    if (catalogNumber < 0) {
    fprintf(consoleHandle,"ERROR: Illegal catalog name %s\n",*argv);
    errorFlag = 1;
  } else {
    strcpy(qualifier,*argv);
  }
  }
    break;






  default:
    fprintf(consoleHandle,"ERROR: * illegal command -%c-",cmdchar);
    errorFlag = 1;

  }
        
  }
  }
  }
    /* Validate arguments */

    if (outfile[0] == 0) {
      fprintf(consoleHandle,"ERROR: No output filename was specified\n");
      errorFlag = 1;
    }
    if (databasefile[0] == 0) {
      fprintf(consoleHandle,"ERROR: No database filename was specified\n");
      errorFlag = 1;
    }
    if (platefile[0] == 0) {
      fprintf(consoleHandle,"ERROR: No plate filename was specified\n");
      errorFlag = 1;
    }

    if (minimumImageCount <= 0) {
      fprintf(consoleHandle,"ERROR: no minimum image count specified\n");
      errorFlag = 1;
    }

    outHandle = fopen(outfile,"a+t");
    if (outHandle == NULL) {
      errorFlag = 1;
      fprintf(consoleHandle,"ERROR: Failed to open the output file %s\n",outfile);
    } else {
      if (verbose) {
        fprintf(consoleHandle,"Output file %s\n",outfile);
      }
    }



    databaseHandle = fopen(databasefile,"wt");
    if (databaseHandle == NULL) {
      errorFlag = 1;
      fprintf(consoleHandle,"ERROR: Failed to open the output file %s\n",databasefile);
    } else {
      if (verbose) {
        fprintf(consoleHandle,"Database file %s\n",databasefile);
      }
      fprintf(databaseHandle,"REF\tgroup\tgroupCount\tDate\tyear\tra\tdec\tmagcal_magdep\tdrad\tPlate\tspatial_bin\tBFLAGS\tGSCBlendFlag\tgsc_bin_index\tnearbyObjects\tlat\tlon\tELLIPTICITY\tTHETA_J2000\tlimiting_mag_local\tAFLAGS\tNUMBER\tdradRMS2\tFWHM_WORLD\trootFlag\tcolorterm\tcolorflag\n");
      fprintf(databaseHandle,"---\t-----\t----------\t----\t----\t--\t---\t-------------\t----\t-----\t-----------\t------\t------------\t-------------\t-------------\t---\t---\t-----------\t-----------\t------------------\t------\t------\t--------\t----------\t--------\t---------\t---------\n");
    }



    astrometryHandle = fopen(astrometryfile,"wt");
    if (astrometryHandle == NULL) {
      errorFlag = 1;
      fprintf(consoleHandle,"ERROR: Failed to open the output file %s\n",astrometryfile);
    } else {
      if (verbose) {
        fprintf(consoleHandle,"Astrometry file %s\n",astrometryfile);
      }
      fprintf(astrometryHandle,"minmag\tdrad\tdra\tddec\n");
      fprintf(astrometryHandle,"------\t----\t---\t----\n");
    }


    plateHandle = fopen(platefile,"wt");
    if (plateHandle == NULL) {
      errorFlag = 1;
      fprintf(consoleHandle,"ERROR: Failed to open the output file %s line %d\n",platefile,__LINE__);
    } else {
      if (verbose) {
        fprintf(consoleHandle,"Plate file %s\n",platefile);
      }
      /* WARNING: Modifications here require modifications to the sscanf file in MergePlateOutlierFile */
      fprintf(plateHandle,"Plate\tseriesId\tplateNumber\tmeasurementCount\tmatchCount\tmatchRatio\tlightcurveCount\toutlierCount\toutlierRatio\n");
      fprintf(plateHandle,"-----\t--------\t-----------\t----------------\t----------\t----------\t---------------\t------------\t------------\n");
    }


    if (idfile[0] != 0) {
      idHandle = fopen(idfile,"wt");
      if (idHandle == NULL) {
        errorFlag = 1;
        fprintf(consoleHandle,"ERROR: Failed to open the output file %s\n",idfile);
      } else {
        if (verbose) {
          fprintf(consoleHandle,"Id file %s\n",idfile);
        }
        /* This header implies that pFileCommon->extraColumnFlag == 1 */
        fprintf(idHandle,"REF\tdrad\tdradB\tmax_drad\tmax_dradB\tyrbegin\tyrend\tnpoints\tmin_local\tmax_local\trange_local\tmin_local2\tmax_local2\trange_local2\tmedian_local\trms_local\tngood\tngoodB\tclip_median_local\tclip_rms_local\tclip_ngood\tmedian_iso\trange_iso\terror_bar_factor\tMalmquist_factor\tMalmquist_factorB\tDamon_factor\tSextractor_Blend\tnblend\tnNonDamonBlue\tmedNonDamonBlue\tnonDamonBluerms\tnDamonBlue\tmedDamonBlue\tdamonBlueRms\tmagvslimitingcorr\tmagvsracorr\tmagvsdeccorr\trmsdradrms2\trarms\tdecrms\tnburst\tnburst2\tnburst3\tnburst4\tndip\tndip2\tndip3\tndip4\tndev3\tndev2\tadjacentburstdip\tadjacentburstdip2\tadjacentburstdip3\tlightcurverms1\tlightcurverms2\tlightcurverms3\tlightcurverms4\tlightcurverms5\tslope_all\tslope_all_err\tnslope_all\tslope_60\tslope_60_err\tnslope_60\trms_factor\tdmagcatalog\tStdmag\tcolor\tra\tdeclination\tMAGFlag\tgscclass\tVFlag\tRaPM\tDecPM\tkeplerField\tversionId\tnearbyObjects\tlat\tnearbyREFflag\treleaseField\tpeakDays\tpeakYear\tpeakMag\tpeakSlope\tpeakRMS\tpeakCount\tpeakUpperCount\tpeakNumber\tpeakMaxMag\tpeakMinMag\tpeakEventCount\tpeakDradRMS\tpeakEventCount2\tpeakEventCount3\tpeakExtra\tpeakCountNF\tpeakCountWF\tpeakOutside\tpeakOutsideNF\tpeakOutsideWF\tpeakDefectCount\tpeakRA\tpeakDec\tpeakNearbyDistance\tpeakMultipleCount\tgsc_bin_index\tpeakDradRMS2\tpeakDradRMS3\tpeakBadColorCount\tpeakLimitingYears\tpeakLimitingPoints\tpeakEvaluation\tpeakNoDefectMag\tpeakLongOutburst\n");
        fprintf(idHandle,"---\t----\t-----\t--------\t---------\t-------\t-----\t-------\t---------\t---------\t-----------\t----------\t----------\t------------\t------------\t---------\t-----\t------\t-----------------\t--------------\t----------\t----------\t---------\t----------------\t----------------\t-----------------\t------------\t----------------\t------\t-------------\t---------------\t---------------\t----------\t------------\t------------\t-----------------\t-----------\t------------\t-----------\t-----\t------\t------\t-------\t-------\t-------\t----\t-----\t-----\t-----\t-----\t-----\t----------------\t-----------------\t-----------------\t--------------\t--------------\t--------------\t--------------\t--------------\t---------\t-------------\t----------\t--------\t------------\t---------\t----------\t-----------\t------\t-----\t--\t-----------\t-------\t--------\t-----\t----\t-----\t-----------\t---------\t-------------\t---\t-------------\t------------\t--------\t--------\t-------\t---------\t-------\t---------\t--------------\t----------\t----------\t----------\t--------------\t-----------\t---------------\t---------------\t----------\t-----------\t-----------\t-----------\t-------------\t-------------\t--------------\t------\t-------\t------------------\t-----------------\t-------------\t------------\t------------\t-----------------\t-----------------\t------------------\t--------------\t---------------\t----------------\n");


      }

    }

    catalogdir = getenv("DASCH_CATALOG");
    if (catalogdir == NULL) {
      fprintf(consoleHandle,"DASCH_CATALOG is not defined\n");
      fprintf(consoleHandle,"0\n");
      return(-1);
    }
    strcpy(catalogname,catalogdir);
    if (strstr(qualifier,"kepler")) {
      slashPtr = strrchr(catalogname,'/');
      if (slashPtr != NULL) {
        slashPtr++;
      } else {
        slashPtr = catalogname;
      }
      *slashPtr = 0;
      strcat(catalogname,"kepler.dat");
    }
    if (strstr(qualifier,"apass")) {
      slashPtr = strrchr(catalogname,'/');
      if (slashPtr != NULL) {
        slashPtr++;
      } else {
        slashPtr = catalogname;
      }
      *slashPtr = 0;
      strcat(catalogname,"apass.dat");
    }
    if (strstr(qualifier,"gaia")) {
      slashPtr = strrchr(catalogname,'/');
      if (slashPtr != NULL) {
        slashPtr++;
      } else {
        slashPtr = catalogname;
      }
      *slashPtr = 0;
      strcat(catalogname,"gaiadr2.dat");
    }
    if (strstr(qualifier,"atlas")) {
      slashPtr = strrchr(catalogname,'/');
      if (slashPtr != NULL) {
        slashPtr++;
      } else {
        slashPtr = catalogname;
      }
      *slashPtr = 0;
      strcat(catalogname,"atlas.dat");
    }
    if (strstr(qualifier,"experimental")) {
      slashPtr = strrchr(catalogname,'/');
      if (slashPtr != NULL) {
        slashPtr++;
      } else {
        slashPtr = catalogname;
      }
      *slashPtr = 0;
      strcat(catalogname,"experimental.dat");
    }
  



    strcpy(indexname,catalogname);
    dotPtr = strstr(indexname,".");
    if (dotPtr != NULL) {
      *dotPtr = 0;
    }
    strcat(indexname,".idx");
  
    strcpy(galaxy_name,catalogdir);
    slashPtr = strrchr(galaxy_name,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = galaxy_name;
    }
    *slashPtr = 0;
    strcat(galaxy_name,"galaxy.dat");

    if (errorFlag) {
      fprintf(consoleHandle,"Usage: search_none options\n");
      fprintf(consoleHandle,"  options: -v verbose\n");    
      fprintf(consoleHandle,"           -i ignore version id\n");
      fprintf(consoleHandle,"           -q <catalog>\n");
      fprintf(consoleHandle,"           -c <minimum image count for magnitude file>\n");
      fprintf(consoleHandle,"           -s <starting gsc bin index>\n");
      fprintf(consoleHandle,"           -o <output file>\n");
      fprintf(consoleHandle,"           -O optimize location of transients\n");
      fprintf(consoleHandle,"           -d <database file>\n");
      fprintf(consoleHandle,"           -u <unmatched summary file>\n");
      fprintf(consoleHandle,"           -p <plate file>\n");
      fprintf(consoleHandle,"           -m <proper motion file>\n");
      fprintf(consoleHandle,"           -a <astrometry>\n");
      fprintf(consoleHandle,"           -l <minimum galactic latitude>, default %f\n",DEFAULT_GAL_LATITUDE);
      fprintf(consoleHandle,"           -t Enable Transient Flare Search\n");


      return(-1);
    }




    /* Connect to the database */

    mysqlhost = getenv("DASCH_MYSQLHOST");
    if (mysqlhost == NULL) {
      fprintf(consoleHandle,"ERROR: DASCH_MYSQLHOST is not defined\n");
      return(-1);
    }
    username = getenv("DASCH_USERNAME");
    if (username == NULL) {
      fprintf(consoleHandle,"ERROR: DASCH_USERNAME is not defined\n");
      return(-1);
    }
    password = getenv("DASCH_PASSWORD");
    if (password == NULL) {
      fprintf(consoleHandle,"ERROR: DASCH_PASSWORD is not defined\n");
      return(-1);
    }

    mysql_init(pConnection);


    if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
      if (mysql_errno(pConnection)) {
        fprintf(consoleHandle,"ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
      }
      return(-1);
    }


    mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
    if (mysqlphothost == NULL) {
      fprintf(consoleHandle,"DASCH_PHOT_MYSQLHOST is not defined\n");
      return(-1);
    }

    photusername = getenv("DASCH_PHOT_USERNAME");
    if (photusername == NULL) {
      fprintf(consoleHandle,"DASCH_PHOT_USERNAME is not defined\n");
      return(-1);
    }
    photpassword = getenv("DASCH_PHOT_PASSWORD");
    if (photpassword == NULL) {
      fprintf(consoleHandle,"DASCH_PHOT_PASSWORD is not defined\n");
      return(-1);
    }
    mysql_init(pPhotConnection);


    if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
      if (mysql_errno(pPhotConnection)) {
        fprintf(consoleHandle,"ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
      }
      return(-1);
    }

    InitSeriesTable(pConnection,pPhotConnection);

    if (GetPhotometryGlobal(pPhotConnection,pPhotGlobal) != 1) {
      fprintf(consoleHandle,"ERROR: failed to get the global photometry table\n");
      exit(-1);
    }
    if (pPhotGlobal->magnitudeFile == PHOT_MAGNITUDEFILE_NO ) {
      fprintf(consoleHandle,"ERROR: magnitudes stored in a MySQL table is no longer supported\n");
      exit(-1);
    }

    InitFileCommon(consoleHandle,pFileCommon,pConnection,pPhotConnection,catalogNumber,1);
    pFileCommon->minimumImageCount = minimumImageCount;
    pFileCommon->minimumGalacticLatitude = minimumGalacticLatitude;
    pFileCommon->verbose = verbose;
    pFileCommon->myrank = myrank;

#ifdef REMATCH_BUILD
    if (enableRematch) {
      printf("ERROR: enableRematch = 1 REMATCH_RESOLUTION %f arcsec\n",REMATCH_RESOLUTION);
      pFileCommon->enableRematch = 1;
    }
#endif /* REMATCH_BUILD */

    if (InitGscBinList(pFileCommon,gscbinfile) != 0) {
      exit(-1);
    }

    pFileCommon->extraColumnFlag = 1;
    if (enableTransientSearch) {
      pFileCommon->enableTransientSearch = 1;
      strcpy(pFileCommon->timestr,timestr2);
      printf("ERROR: enableTransientSearch = 1\n");
    }

  if (flareCandidatesOnly) {
#if !defined(START_GSC_BIN_INDEX) && !defined(END_GSC_BIN_INDEX)
    curGscBinIndex = endGscBinIndex;
    for (pFileCommon->interesting_index = 0; pFileCommon->interesting_index < pFileCommon->interesting_count; pFileCommon->interesting_index) {
      pFlareCandidate = &pFileCommon->candidate_table[pFileCommon->interesting_index];
      if (pFlareCandidate->gsc_bin_index <  curGscBinIndex) {
        curGscBinIndex = pFlareCandidate->gsc_bin_index;
        startGscBinIndex = pFlareCandidate->gsc_bin_index;
      }
      if (pFlareCandidate->gsc_bin_index >  endGscBinIndex) {
        endGscBinIndex = pFlareCandidate->gsc_bin_index;
      }
      
    }
    if (endGscBinIndex <= startGscBinIndex) {
      printf("ERROR: line %5d startGscBinIndex %d is greater than endGscBinIndex %d\n",__LINE__,startGscBinIndex,endGscBinIndex);
      exit(-1);      
    }
    fprintf(consoleHandle,"search_none line %5d flareCandidatesOnly %d interesting_count %2d startGscBinIndex %10d endGscBinIndex %10d\n",
            __LINE__,
            pFileCommon->flareCandidatesOnly,
            pFileCommon->interesting_count,
            startGscBinIndex,
            endGscBinIndex);
#endif
  }



    if (flareCandidatesOnly) {
      pFileCommon->flareCandidatesOnly = 1;
      strcpy(pFileCommon->timestr,timestr2);
      printf("ERROR: flareCandidatesOnly = 1\n");
    }

    MPI_Comm_size(MPI_COMM_WORLD, &pFileCommon->maxproc);
    MPI_Get_processor_name(pFileCommon->processorName, &pFileCommon->namelength);
    if ((charPtr = strstr(pFileCommon->processorName,".rc.fas.harvard.edu")) != NULL) {
      *charPtr = 0;
      pFileCommon->namelength -= 19;
    }

    if (catalogNumber != 0) {
      sprintf(pFileCommon->catalogString,"%d",catalogNumber);
    } else {
      pFileCommon->catalogString[0] = 0;
    }
    if (strlen(pFileCommon->catalogString) > (MAX_CATALOG_STRING-2)) {
      fprintf(consoleHandle,"ERROR: catalogString of size %d is too large\n",strlen(pFileCommon->catalogString));
      exit(-1);
    }
    InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);

    TimeStampPhotometryGlobal(pPhotConnection);

    if (pFileCommon->myrank == 0) {
      fprintf(consoleHandle,"search_none of %s %s, outfile %s databasefile %s platefile %s propermotionfile %s astrometryfile %s bin histogramfile %s magnitudeLimit %d PHOTSTARIMAGE: %d FILECOMMON %d WORKQUEUE %d curGscBinIndex %d minVersionId%d %d minimumGalacticLatitude %f\n",
              __DATE__,__TIME__,outfile,databasefile,platefile,propermotionfile,astrometryfile,binhistogramfile,MAGNITUDE_LIMIT,sizeof(PHOTSTARIMAGE),sizeof(FILECOMMON),sizeof(WORKQUEUE),curGscBinIndex,catalogNumber,pPhotGlobal->minVersionId[catalogNumber],pFileCommon->minimumGalacticLatitude);
      fprintf(outHandle,    "search_none of %s %s, outfile %s databasefile %s platefile %s propermotionfile %s astrometryfile %s bin histogramfile %s magnitudeLimit %d PHOTSTARIMAGE %d FILECOMMON %d WORKQUEUE %d curGscBinIndex %d minVersionId%d %d minimumGalacticLatitude %f\n",
              __DATE__,__TIME__,outfile,databasefile,platefile,propermotionfile,astrometryfile,binhistogramfile,MAGNITUDE_LIMIT,sizeof(PHOTSTARIMAGE),sizeof(FILECOMMON),sizeof(WORKQUEUE),curGscBinIndex,catalogNumber,pPhotGlobal->minVersionId[catalogNumber],pFileCommon->minimumGalacticLatitude);
    }

    time(&startTime);
    oldTime = startTime;
    ptr = localtime(&startTime);
    strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
 
    if (pFileCommon->myrank == 0) {
      fprintf(consoleHandle,"search_none starting at %s\n",timestr);
      fprintf(outHandle,"search_none starting at %s\n",timestr);
      fflush(outHandle);
    }

    if (PopulateGalaxyTree(pGalaxyCommon,galaxy_name,pFileCommon->enableTransientSearch) != 0) {
      exit(-1);
    }


    sprintf(queryString,"update photglobal set keepalive = 'yes';");
    res = ExecuteQuery(pPhotConnection,queryString);
    if (res) {
      exit(-1);
    }




    if (pFileCommon->myrank == 0) {
      /* Find out how many processes there are in the default
         communicator */
#ifdef DEBUG_PARENT
      pFileCommon->maxproc = 31;
#else /* DEBUG_PARENT */
      MPI_Comm_size(MPI_COMM_WORLD, &pFileCommon->maxproc);
#endif /* DEBUG_PARENT */
      fprintf(consoleHandle,"maxproc is %d\n",pFileCommon->maxproc);
      if (pFileCommon->maxproc == 1) {
        pFileCommon->startrank = 0;
        gscIncrement = (endGscBinIndex-startGscBinIndex);
      } else {
        pFileCommon->startrank = 1;
        gscIncrement = (endGscBinIndex-startGscBinIndex)/(pFileCommon->maxproc-1);
      }
      if (binhistogramfile[0] != 0) {
        ReadGscBinHistogram(pFileCommon,consoleHandle,startGscBinIndex,endGscBinIndex,binhistogramfile,binGscIndex,binCount,&binArraySize,&binTotal);
        fprintf(consoleHandle,"binTotal is %lld in %s\n",binTotal,binhistogramfile);
        if (binTotal < pFileCommon->maxproc) {
          fprintf(consoleHandle,"WARNING: dividing bins evenly\n");        
          binhistogramfile[0] = 0;
        }
      }



      /* Now allocate the work queue */
      pFileCommon->workQueueSize = pFileCommon->maxproc;
      pWorkQueue = (PWORKQUEUE)calloc(pFileCommon->workQueueSize,sizeof(WORKQUEUE));
      if (pWorkQueue == NULL) {
        fprintf(consoleHandle,"ERROR: failed to allocate the work queue\n");
        exit(-1);
      }
   
      scriptdir = getenv("DASCH_SCRIPTS");
      if (scriptdir == NULL) {
        printf("ERROR: DASCH_SCRIPTS is not defined\n");
        exit(-1);
      }
      /* Now create all of the file names */
      beginGscIndex = startGscBinIndex;
      for (workIndex = 0; workIndex < pFileCommon->maxproc; workIndex++) {
        pWorkEntry = &pWorkQueue[workIndex];
        pWorkEntry->queueIndex = workIndex;
        sprintf(pWorkEntry->consoleLogname,"%s/search_none_console_%02d.log",scriptdir,workIndex);
#if 0
        printf("Line %d console file for workIndex %d is %s\n",__LINE__,workIndex,pWorkEntry->consoleLogname);
#endif
        sprintf(pWorkEntry->outLogname,"%s_%02d.log",outfile,workIndex);
        sprintf(pWorkEntry->databaseLogname,"%s_%02d.log",databasefile,workIndex);
        sprintf(pWorkEntry->propermotionLogname,"%s_%02d.log",propermotionfile,workIndex);
        sprintf(pWorkEntry->plateLogname,"%s_%02d.log",platefile,workIndex);
        if (idfile[0] != 0) {
          sprintf(pWorkEntry->idLogname,"%s_%02d.log",idfile,workIndex);
        }
        sprintf(pWorkEntry->astrometryLogname,"%s_%02d.log",astrometryfile,workIndex);    
        pWorkEntry->startGscBinIndex = startGscBinIndex;
        pWorkEntry->endGscBinIndex = endGscBinIndex;
        if (workIndex > 0) {
          if (binhistogramfile[0] != 0) {
            /* Change of August 12, 2019 Divide by maxproc-n (last/highest ratio needs to be 0.47 from run of 2019_07_02*/
            /* n = 1    8167076/8072819 = 1.011 (last/highest) */
            /* n = 1.5  4600016/7978949 = 0.576                */
            /* n = 1.6  4036796/7978949 = 0.505                */
            /* n = 1.65 3661316/7978949 = 0.458                */
            /* n = 1.7  3379706/7978949   0.423                */
            /* n = 2.0  1596176/8166689 = 0.195                */
#if 1
            binTargetCount = (100 * binTotal * workIndex)/((100*pFileCommon->maxproc)-165); /* maxproc-1.65 */
#else
            binTargetCount = (binTotal * workIndex)/(pFileCommon->maxproc-1);
#endif
            pWorkEntry->startGscBinIndex = beginGscIndex;
            while (curBinIndex < binArraySize) {
              binCurrentCount += binCount[curBinIndex];
              if (binCurrentCount > binTargetCount) {
                break;
              }
              curBinIndex++;
            }
            if (curBinIndex < binArraySize) {
              if (workIndex == (pFileCommon->maxproc-1)) {
                pWorkEntry->endGscBinIndex = endGscBinIndex;
              } else {
                pWorkEntry->endGscBinIndex = binGscIndex[curBinIndex] -1;
                beginGscIndex =  binGscIndex[curBinIndex];
              }
            } else {
              if (workIndex == (pFileCommon->maxproc-1)) {
                pWorkEntry->endGscBinIndex = endGscBinIndex;
              } else {
                pWorkEntry->endGscBinIndex = endGscBinIndex;
                beginGscIndex =  endGscBinIndex;
              }
            }
            if (pWorkEntry->endGscBinIndex > endGscBinIndex) {
              pWorkEntry->endGscBinIndex = endGscBinIndex;
            }

          } else {

            pWorkEntry->startGscBinIndex = startGscBinIndex+(gscIncrement*(workIndex-1));
            if (workIndex == (pFileCommon->maxproc-1)) {
              pWorkEntry->endGscBinIndex = endGscBinIndex;
            } else {
              pWorkEntry->endGscBinIndex = startGscBinIndex+(gscIncrement*(workIndex));
            }

          }
        }
#if 0
        printf("workIndex %2d startGscBinIndex %10d endGscBinIndex %10d total bins %10lld binTargetCount %13lld\n",workIndex,pWorkEntry->startGscBinIndex,pWorkEntry->endGscBinIndex,pWorkEntry->endGscBinIndex-pWorkEntry->startGscBinIndex,binTargetCount);
#endif
      }
      PlotSymbolKey(timestr2,pFileCommon->enableRematch);
#ifdef DEBUG_PARENT
      printf("ERROR: DEBUG_MASTER exiting in line %d\n",__LINE__);
      exit(-1);
#endif /* DEBUG_PARENT */
      master(consoleHandle,pFileCommon,pWorkQueue,outHandle,databaseHandle,plateHandle,idHandle,astrometryHandle,propermotionfile,qualifier);
    } else {    
      slave(consoleHandle,pFileCommon,qualifier);
    }
    mysql_close(pConnection);
    mysql_close(pPhotConnection);
  
    if (pFileCommon->myrank == 0) {
      FreeFileCommon(pFileCommon,1);
    } else {
      FreeFileCommon(pFileCommon,0);
    }
    if (outHandle != NULL) {
      fclose(outHandle);
    }
    if (databaseHandle != NULL) {
      fclose(databaseHandle);
    }
    if (plateHandle != NULL) {
      fclose(plateHandle);
    }
    if (idHandle != NULL) {
      fclose(idHandle);
    }
    if (flareCandidatesOnly) {
      kd_free(candidate_ptree);
    }

    /* Shut down MPI */
    MPI_Finalize();

    return(0);
  } /* end of main() */
void  GetDiskCounts(FILE *consoleHandle,long long *prxbytes, long long *ptxbytes)
{
  char filename[] = "/proc/diskstats";
  FILE *filehandle = NULL;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  int lineNumber;
  long long  rxbytes;
  long long  txbytes;
  long long  skipval;
  char *charPtr;
  int nvals;
  int verbose = 0;
  *prxbytes = 0;
  *ptxbytes = 0;
  filehandle = fopen(filename,"rt");
  if (filehandle == NULL) {
    fprintf(consoleHandle,"ERROR opening %s\n",filename);
    return;
  }
  lineNumber = 0;
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,filehandle);
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
    lineNumber++;
    if (((charPtr = strstr(inBuffer,"sdc1")) != NULL) ||
        ((charPtr = strstr(inBuffer,"sdb")) != NULL)) { 
      if (verbose) {
        fprintf(consoleHandle,"%s\n",inBuffer);
      }
      charPtr += 4;
      
      nvals = sscanf(charPtr,"%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld",
                     &skipval,&skipval,&rxbytes,&skipval,
                     &skipval,&skipval,&txbytes,&skipval,
                     &skipval,&skipval,&skipval);
      if (nvals != 11) {
        fprintf(consoleHandle,"ERROR: nvals is %d in line %s\n",charPtr);
        
      } else {
        rxbytes = rxbytes * 1024L;
        txbytes = txbytes * 1024L;
        if (verbose) {
          fprintf(consoleHandle,"Disk rx %6lld tx %6lld\n",rxbytes,txbytes);
        }
        *prxbytes += rxbytes;
        *ptxbytes += txbytes;
            

        if (verbose) {
          fprintf(consoleHandle,"\n");
        }
        
      }
    }
  }
  fclose(filehandle);



}
 
void  GetNetworkCounts(FILE *consoleHandle,long long *prxbytes, long long *ptxbytes)
{
  char filename[] = "/proc/net/dev";
  FILE *filehandle = NULL;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  int lineNumber;
  long long  rxbytes;
  long long  txbytes;

  long long  eth_rxbytes = 0;
  long long  eth_txbytes = 0;
  long long  bond_rxbytes = 0;
  long long  bond_txbytes = 0;
  long long  ib_rxbytes = 0;
  long long  ib_txbytes = 0;


  long long  skipval;
  char *charPtr;
  int nvals;
  int verbose = 0;
  *prxbytes = 0;
  *ptxbytes = 0;
  filehandle = fopen(filename,"rt");
  if (filehandle == NULL) {
    fprintf(consoleHandle,"ERROR opening %s\n",filename);
    return;
  }
  lineNumber = 0;
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,filehandle);
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
    lineNumber++;
    if ((charPtr = strstr(inBuffer,":")) != NULL) { 
      if (verbose) {
        fprintf(consoleHandle,"%s\n",inBuffer);
      }
      *charPtr = 0;
      charPtr++;
      
      nvals = sscanf(charPtr,"%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld",
                     &rxbytes,&skipval,&skipval,&skipval,&skipval,&skipval,&skipval,&skipval,
                     &txbytes,&skipval,&skipval,&skipval,&skipval,&skipval,&skipval,&skipval);
      if (nvals != 16) {
        fprintf(consoleHandle,"ERROR: nvals is %d in line %s\n",charPtr);
        
      } else {
        if (verbose) {
          fprintf(consoleHandle,"Link %s rx %6lld tx %6lld\n",inBuffer,rxbytes,txbytes);
        }
        if ((strstr(inBuffer,"eth") != 0)  || (verbose != 0)) {
          eth_rxbytes += rxbytes;
          eth_txbytes += txbytes;         
        } else if ((strstr(inBuffer,"bond") != 0)  || (verbose != 0)) {
          bond_rxbytes += rxbytes;
          bond_txbytes += txbytes;
        } else if ((strstr(inBuffer,"ib") != 0) || (verbose != 0)) {
          ib_rxbytes += rxbytes;
          ib_txbytes += txbytes;
        }

      }
    }
  }
  if ((bond_rxbytes+bond_txbytes) != 0) {    
    *prxbytes = bond_rxbytes + ib_rxbytes;
    *ptxbytes = bond_txbytes + ib_txbytes;
  } else {
    *prxbytes = eth_rxbytes + ib_rxbytes;
    *ptxbytes = eth_txbytes + ib_txbytes;
  }
  fclose(filehandle);



}
void ReadLogFile(FILE *consoleHandle,PWORKQUEUE pWorkEntry) {
  int statResult;
  struct stat filestats;
  int iteration = 10;
  FILE * logFile;
  FILE * altFile;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  char altfilename[MAX_BUFFER];
  int altLogIndex;
  int lineCounter = 0;
  
#if 0
  printf("Line %d console file for workIndex %d is %s\n",__LINE__,pWorkEntry->queueIndex,pWorkEntry->consoleLogname);
#endif
  if (pWorkEntry->queueIndex == 0) {
    /* This is the master acting alone, nothing to read */
    return;
  }

#if 0
  printf("Line %d console file for workIndex %d is %s\n",__LINE__,pWorkEntry->queueIndex,pWorkEntry->consoleLogname);
#endif

  altfilename[0] = 0;
  statResult = -1;
  if (pWorkEntry->errorFlag) {
    fprintf(consoleHandle,"ERROR: errorFlag is set for %s\n",pWorkEntry->consoleLogname);
  } else {
    while (statResult != 0) {
      statResult = stat(pWorkEntry->consoleLogname,&filestats);
      if (statResult !=0) {
        fprintf(consoleHandle,"WARNING: master did not find %s\n",pWorkEntry->consoleLogname);
      } else if (filestats.st_size !=  pWorkEntry->consoleFilesize) {
        fprintf(consoleHandle,"WARNING: master expected size %d got size %d for %s\n",
                pWorkEntry->consoleFilesize,
                filestats.st_size,
                pWorkEntry->consoleLogname);
        statResult = -1;
      } else {
        break;
      }
      sleep(1);
      if (iteration-- <= 0) {
        break;
      }
    
    }
    if (iteration <= 0) {
      fprintf(consoleHandle,"ERROR: Master could not read %s of size %d\n",pWorkEntry->consoleLogname,(int)pWorkEntry->consoleFilesize);
      return;
    }  
    logFile = fopen(pWorkEntry->consoleLogname,"rt");
    if (logFile == NULL) {
      fprintf(consoleHandle,"ERROR: Master could not open %s\n",pWorkEntry->consoleLogname);
      return;
    }
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,logFile);
      if (inBuffer == NULL) {
        break;
      }
      lineCounter++;
      lineLen = strlen(inBuffer);
      if (lineLen > (MAX_BUFFER-10)) {
          printf("ERROR: MAX_BUFFER exceeded in line %d of %s\n",lineCounter,pWorkEntry->consoleLogname);
      }

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
      fprintf(consoleHandle,"%s\n",inBuffer);

    }
    fclose(logFile);
#if 0
    unlink(pWorkEntry->consoleLogname);
#endif
  }
}
void ReadResultFile(FILE *consoleHandle,FILE *resultHandle,int queueIndex,char *resultLogname,int errorFlag,off_t resultFilesize) {
  int statResult;
  struct stat filestats;
  int iteration = 10;
  FILE * logFile;
  FILE * altFile;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  char altfilename[MAX_BUFFER];
  int altLogIndex;
  int lineCounter = 0;

#if 0
  printf("Line %d out file for workIndex %d is %s\n",__LINE__,queueIndex,resultLogname);
#endif
  if (resultLogname[0] == 0) {
    return;
  }


  altfilename[0] = 0;
  statResult = -1;
  if (errorFlag) {
    fprintf(consoleHandle,"ERROR: errorFlag is set for %s\n",resultLogname);
  } else {
    while (statResult != 0) {
      statResult = stat(resultLogname,&filestats);
      if (statResult !=0) {
        fprintf(consoleHandle,"WARNING: master did not find %s\n",resultLogname);
      } else if (filestats.st_size !=  resultFilesize) {
        fprintf(consoleHandle,"WARNING: master expected size %d got size %d for %s\n",
                resultFilesize,
                filestats.st_size,
                resultLogname);
        statResult = -1;
      } else {
        break;
      }
      sleep(1);
      if (iteration-- <= 0) {
        break;
      }
    
    }
    if (iteration <= 0) {
      fprintf(consoleHandle,"ERROR: Master could not read %s of size %d\n",resultLogname,(int)resultFilesize);
      return;
    }  
    logFile = fopen(resultLogname,"rt");
    if (logFile == NULL) {
      fprintf(consoleHandle,"ERROR: Master could not open %s\n",resultLogname);
      return;
    }
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,logFile);
      if (inBuffer == NULL) {
        break;
      }
      lineCounter++;
      lineLen = strlen(inBuffer);
      if (lineLen > (MAX_BUFFER-10)) {
        printf("ERROR: MAX_BUFFER exceeded in line %d of %s\n",lineCounter,resultLogname);
      }


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
      fprintf(resultHandle,"%s\n",inBuffer);

    }
    fclose(logFile);
    unlink(resultLogname);
  }
}
void MergePlateOutlierFile(FILE *consoleHandle,PFILECOMMON pFileCommon,PVECTORSTORE pVectorStore,FILE *plateHandle,int queueIndex,char *plateLogname,int errorFlag,off_t plateFilesize) {
  int statResult;
  struct stat filestats;
  int iteration = 10;
  FILE * logFile;
  FILE * altFile;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  char altfilename[MAX_BUFFER];
  int altLogIndex;
  int lineCounter = 0;
  File plate_handle = NULL;
  char plate_name[MAX_BUFFER];
  TableHead plate_header = NULL;
  PPLATEOUTLIER plate_table = NULL;
  size_t plate_nrecs = 0;
  int plate_index;
  PPLATEOUTLIER pPlateOutlier;
  strcpy(plate_name,plateLogname);
  PPLATECOUNT plateCountTable;

#if 0
  printf("Line %d out file for workIndex %d is %s\n",__LINE__,queueIndex,plateLogname);
#endif
  if (plateLogname[0] == 0) {
    return;
  }


  altfilename[0] = 0;
  statResult = -1;
  if (errorFlag) {
    fprintf(consoleHandle,"ERROR: errorFlag is set for %s\n",plateLogname);
  } else {
    while (statResult != 0) {
      statResult = stat(plateLogname,&filestats);
      if (statResult !=0) {
        fprintf(consoleHandle,"WARNING: master did not find %s\n",plateLogname);
      } else if (filestats.st_size !=  plateFilesize) {
        fprintf(consoleHandle,"WARNING: master expected size %d got size %d for %s\n",
                plateFilesize,
                filestats.st_size,
                plateLogname);
        statResult = -1;
      } else {
        break;
      }
      sleep(1);
      if (iteration-- <= 0) {
        break;
      }
    
    }
    if (iteration <= 0) {
      fprintf(consoleHandle,"ERROR: Master could not read %s of size %d\n",plateLogname,(int)plateFilesize);
      return;
    }
    plate_handle = Open(plate_name,"r");
    if (plate_handle == NULL) {
      fprintf(consoleHandle,"ERROR: No plate outlier file found  %s\n",plate_name);
      errorFlag = 1;
    } else {
      if (pFileCommon->verbose) {
        fprintf(consoleHandle,"Found plate plate outlier file %s\n",plate_name);
      }
    }
    /* Now read in the plate outlier file */
    plate_header = table_header(plate_handle,TABLE_PARSE);
    if (plate_header == NULL) {
      fprintf(consoleHandle,"ERROR: Failed to read header for %s\n",plate_name);
      return;
    }
    plate_table = table_loadva(plate_handle,
                               &plate_header,
                               NULL, /* hbase */
                               NULL, /* rows */
                               NULL,
                               sizeof(PLATEOUTLIER),
                               &plate_nrecs,
                               TblInt,"seriesId",TblOff(PPLATEOUTLIER,seriesId),
                               TblInt,"plateNumber",TblOff(PPLATEOUTLIER,plateNumber),
                               TblInt,"measurementCount",TblOff(PPLATEOUTLIER,measurementCount),
                               TblInt,"matchCount",TblOff(PPLATEOUTLIER,matchCount),
                               TblInt,"lightcurveCount",TblOff(PPLATEOUTLIER,lightcurveCount),
                               TblInt,"outlierCount",TblOff(PPLATEOUTLIER,outlierCount),
                               TblDbl,"matchRatio",TblOff(PPLATEOUTLIER,matchRatio),
                               TblDbl,"outlierRatio",TblOff(PPLATEOUTLIER,outlierRatio),
                               TblBuf,"Plate",TblOff(PPLATEOUTLIER,Plate),MAX_PLATE_NAME,
                               0,"end",0);
    if (plate_table == NULL) {
      fprintf(consoleHandle,"WARNING: Failed to read table for %s\n",plate_name);
      return;
    }
    for (plate_index = 0; plate_index < plate_nrecs; plate_index++) {
      pPlateOutlier = &plate_table[plate_index];
      plateCountTable = pVectorStore->plateCountTable[pPlateOutlier->seriesId];
      if (plateCountTable == NULL) {
        plateCountTable = (PPLATECOUNT)calloc(pFileCommon->maxPlateNumber[pPlateOutlier->seriesId]+1,sizeof(PLATECOUNT));
        if (plateCountTable == NULL) {
          fprintf(consoleHandle,"ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[pPlateOutlier->seriesId]+1);
          exit(-1);
        }
        pVectorStore->plateCountTable[pPlateOutlier->seriesId] = plateCountTable;
      }
      if ((pPlateOutlier->plateNumber <= 0) ||
          (pPlateOutlier->plateNumber > pFileCommon->maxPlateNumber[pPlateOutlier->seriesId])) {
        fprintf(consoleHandle,"ERROR: illegal plate number %d\n",pPlateOutlier->plateNumber);
        exit(-1);
      }
      plateCountTable[pPlateOutlier->plateNumber].measurementCount += pPlateOutlier->measurementCount;
      plateCountTable[pPlateOutlier->plateNumber].matchCount += pPlateOutlier->matchCount;
      plateCountTable[pPlateOutlier->plateNumber].lightcurveCount += pPlateOutlier->lightcurveCount;
      plateCountTable[pPlateOutlier->plateNumber].outlierCount += pPlateOutlier->outlierCount;      
    }

    /* Done with the plate outlier table */
    if (plate_table != NULL) {
      Free(plate_table);
      plate_table = NULL;
    }
    if (plate_header != NULL) {
      table_hdrfree(plate_header);
      plate_header = NULL;
    }
    if (plate_handle != NULL) {
      Close(plate_handle);
      plate_handle = NULL;
    }


    unlink(plateLogname);
  }
}

void MergeDradHistogramFile(FILE *consoleHandle,PFILECOMMON pFileCommon,PVECTORSTORE pVectorStore,FILE *astrometryHandle,int queueIndex,char *astrometryLogname,int errorFlag,off_t astrometryFilesize) {
  int statResult;
  struct stat filestats;
  int iteration = 10;
  FILE * logFile;
  FILE * altFile;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  char altfilename[MAX_BUFFER];
  int altLogIndex;
  int lineCounter = 0;
  FILE * resultHandle = NULL;
  char astrometry_name[MAX_BUFFER];
  int nvals;
  strcpy(astrometry_name,astrometryLogname);

#if 0
  printf("Line %d out file for workIndex %d is %s\n",__LINE__,queueIndex,astrometryLogname);
#endif
  if (astrometryLogname[0] == 0) {
    return;
  }


  altfilename[0] = 0;
  statResult = -1;
  if (errorFlag) {
    fprintf(consoleHandle,"ERROR: errorFlag is set for %s\n",astrometryLogname);
  } else {
    while (statResult != 0) {
      statResult = stat(astrometryLogname,&filestats);
      if (statResult !=0) {
        fprintf(consoleHandle,"WARNING: master did not find %s\n",astrometryLogname);
      } else if (filestats.st_size !=  astrometryFilesize) {
        fprintf(consoleHandle,"WARNING: master expected size %d got size %d for %s\n",
                astrometryFilesize,
                filestats.st_size,
                astrometryLogname);
        statResult = -1;
      } else {
        break;
      }
      sleep(1);
      if (iteration-- <= 0) {
        break;
      }
    
    }
    if (iteration <= 0) {
      fprintf(consoleHandle,"ERROR: Master could not read %s of size %d\n",astrometryLogname,(int)astrometryFilesize);
      return;
    }
    resultHandle = fopen(astrometryLogname,"rt");
    if (resultHandle == NULL) {
      fprintf(consoleHandle,"ERROR: Master could not open %s\n",astrometryLogname);
      return;
    }
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,resultHandle);
      if (inBuffer == NULL) {
        break;
      }
      lineCounter++;
      lineLen = strlen(inBuffer);
      if (lineLen > (MAX_BUFFER-10)) {
        printf("ERROR: MAX_BUFFER exceeded in line %d of %s\n",lineCounter,astrometryLogname);
      }


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
      int magvalue;
      int dradvalue;
      long long draHistogram;
      long long ddecHistogram;
      int tempdraddec;
      if (lineCounter > 2) {
        nvals = sscanf(inBuffer,"%d\t%d\t%lld\t%lld",&magvalue,&dradvalue,&draHistogram,&ddecHistogram);
        if (nvals == 4) {
          if ((magvalue < 0) || (magvalue >= MAX_DRAD_HISTOGRAM_MAG)) {
            fprintf(consoleHandle,"ERROR: range error for magvalue %d\n",magvalue);
            exit(-1);
          }
          tempdraddec = dradvalue + MAX_DRAD_RANGE;
          if ((tempdraddec < 0) || (tempdraddec > (2*(MAX_DRAD_RANGE)+1))) {
            fprintf(consoleHandle,"ERROR: range error for tempdraddec %d\n",tempdraddec);
            exit(-1);
          }
          if (pFileCommon->maxproc == 1) {
            if (pFileCommon->draHistogram[magvalue][tempdraddec] != draHistogram) {
              printf("ERROR: draHistogram mismatch at %d %d %lld\n",magvalue,tempdraddec,draHistogram);
              exit(-1);
            }
            if (pFileCommon->ddecHistogram[magvalue][tempdraddec] != ddecHistogram) {
              printf("ERROR: ddecHistogram mismatch at %d %d %lld\n",magvalue,tempdraddec,ddecHistogram);
              exit(-1);
            }

          } else {
            pFileCommon->draHistogram[magvalue][tempdraddec] += draHistogram;
            pFileCommon->ddecHistogram[magvalue][tempdraddec] += ddecHistogram;
          }
        } else {
          fprintf(consoleHandle,"ERROR: nvals is %d in line %s for %s\n",nvals,inBuffer,astrometryLogname);
          exit(-1);
        }
      }


    }
    fclose(resultHandle);
   


    unlink(astrometryLogname);
  }
}
void UpdateVectorCounterBlock(PVECTORCOUNTERBLOCK pTargetCounterBlock,PVECTORCOUNTERBLOCK pSourceCounterBlock)
{
  if (pSourceCounterBlock->vectorAlloc > pTargetCounterBlock->vectorAlloc) {
    pTargetCounterBlock->vectorAlloc = pSourceCounterBlock->vectorAlloc;
  }
#if 0
  printf("line %d source vectorAlloc %d -> target vectorAlloc %d\n",__LINE__,pSourceCounterBlock->vectorAlloc,pTargetCounterBlock->vectorAlloc);
#endif
  pTargetCounterBlock->outlierCount += pSourceCounterBlock->outlierCount;
  pTargetCounterBlock->lightcurvePointCount += pSourceCounterBlock->lightcurvePointCount;
  pTargetCounterBlock->lightcurveCount += pSourceCounterBlock->lightcurveCount;
  pTargetCounterBlock->qualityRejectCount += pSourceCounterBlock->qualityRejectCount;
  pTargetCounterBlock->pointsRead2 += pSourceCounterBlock->pointsRead2;
  pTargetCounterBlock->pointsRead += pSourceCounterBlock->pointsRead;
  pTargetCounterBlock->pointsRejected += pSourceCounterBlock->pointsRejected;
  pTargetCounterBlock->insufficientPoints += pSourceCounterBlock->insufficientPoints;
  pTargetCounterBlock->highRMSPoints += pSourceCounterBlock->highRMSPoints;
  pTargetCounterBlock->totalStars += pSourceCounterBlock->totalStars;
  pTargetCounterBlock->totalOutliers += pSourceCounterBlock->totalOutliers;
  pTargetCounterBlock->totalPlot2 += pSourceCounterBlock->totalPlot2;
  pTargetCounterBlock->totalNoneCandidates += pSourceCounterBlock->totalNoneCandidates;

}
void UpdateCatalogDist(int* pTarget,int* pSource) 
{
  int index;
  for (index = 0; index < MAX_TC_DIST; index++) {
    pTarget[index] += pSource[index];
  }
}


int master(FILE *consoleHandle,PFILECOMMON pFileCommon,PWORKQUEUE pWorkQueue,FILE *outHandle,FILE *databaseHandle,FILE *plateHandle,FILE *idHandle,FILE *astrometryHandle,char *propermotionfile,char *qualifier)
{
  PWORKQUEUE pWorkEntry;
  PWORKQUEUE pWorkEntry0 = &pWorkQueue[0];
  PWORKQUEUE pWorkEntry2;
  int rank;
  PWORKQUEUE pReceiveWorkEntry;
  int totalSeconds = 0;
  int newTotalSeconds = 0;
  PWORKQUEUE preceive_result;
  char receive_buffer[MAX_MESSAGE_SIZE];
  MPI_Status status;
  int mpitestflag;
  MPI_Request receiveRequest;
  MPI_Request *sendRequestTable = NULL;
  int sendRequestIndex;
  PWORKQUEUE results;
  int mpiresult;
  int numpending = 0;
  int numstarted = 0;
  int numprocessed = 0;
  int workIndex;
  struct stat statbuf;
  int result;
  time_t startTime;
  time_t curTime;
  int sendNewWorkRequest = 0;
  PNODESTATS pstatusmessage;
  double averageTime;
  char rxstring[MAX_BUFFER];
  char txstring[MAX_BUFFER];
  char rxdiskstring[MAX_BUFFER];
  char txdiskstring[MAX_BUFFER];
  int count;
  int cpuIndex;
  int doneFlag = 0;
  int completedFlag = 0;
  int tempCompletedFlag = 0;
  FILE * keepalivefile;
  char *scriptdirectory;
  PNODESTATS pdiemessage;
  VECTORSTORE vectorStore;
  PVECTORSTORE pVectorStore = &vectorStore;
  int seriesIndex;
  int plateIndex;
  PPLATECOUNT plateCountTable;
  int index;
  int dradIndex;

  memset(pVectorStore,0,sizeof(VECTORSTORE));

  if (pFileCommon->verbose) {
    fprintf(consoleHandle,"Have maxproc %d\n",pFileCommon->maxproc);
  }
  sendRequestTable = (MPI_Request *)calloc(sizeof(MPI_Request),pFileCommon->maxproc);
  if (sendRequestTable == NULL) {
    fprintf(consoleHandle,"ERROR: failed to allocate sendRequestTable of size %d\n",pFileCommon->maxproc);
    exit(-1);
  }

  scriptdirectory = getenv("DASCH_SCRIPTS");
  strcpy(pFileCommon->keepalivename,scriptdirectory);
  strcat(pFileCommon->keepalivename,"/search_none_keepalive.txt");
  keepalivefile = fopen(pFileCommon->keepalivename,"wt");
  if (keepalivefile == NULL) {
    fprintf(consoleHandle,"ERROR Failed to open keepalive file %s\n",pFileCommon->keepalivename);
    exit(-1);
  } else {
    fprintf(keepalivefile,"Delete this file to stop pipeline processing\n");
    fclose(keepalivefile);
    keepalivefile = NULL;
  }

  if ((pFileCommon->workQueueSize <= 0) ||
      (pFileCommon->canceledFlag != 0)) {
    /* Nothing to do! */
    return(0);
  }
  pFileCommon->queueIteration++;
  pFileCommon->completeCount = 0;
  time(&startTime);

  /* Seed the slaves by asking them to report their status */

  pFileCommon->statusRequestNumber++;
  for (rank = pFileCommon->startrank; rank < pFileCommon->maxproc; ++rank) {
    pWorkEntry = &pWorkQueue[rank];
    if (pWorkEntry->entryState == STATE_INIT) {
      mpiresult = MPI_Isend(0, 0, MPI_BYTE, rank, STATUSTAG, MPI_COMM_WORLD,&sendRequestTable[rank]);
#if 0
      printf("Line: %d sendRequestTable[%d] 0x%x\n",__LINE__,rank,sendRequestTable[rank]);
#endif
      numpending++;
      if (mpiresult != MPI_SUCCESS) {
        fprintf(consoleHandle,"ERROR: Master sending status request to %d with result %d\n",rank,mpiresult);
        pWorkEntry->entryState = STATE_ERROR;
      } else {
        pWorkEntry->entryState = STATE_STATUS_SENT;
      }
    }
  }


  /* Loop over getting new work requests until there is no more work
     to be done */
  while (doneFlag == 0) {
    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"line %d numpending %d, numstarted %d, numprocessed %d (3)\n",__LINE__,numpending,numstarted,numprocessed);
      for (workIndex = 0; workIndex < pFileCommon->maxproc; ++workIndex) {
        pWorkEntry = &pWorkQueue[workIndex];
        fprintf(consoleHandle,"WorkIndex %d entryState %d %s\n",workIndex,pWorkEntry->entryState,stateString[pWorkEntry->entryState]);
      }
    }
    

    /* Receive results from a slave */
    sleep(1);
    for (sendRequestIndex = 0; sendRequestIndex < pFileCommon->maxproc; sendRequestIndex++) {
      memset(&status,0,sizeof(status));
#if 0
      printf("Line: %d sendRequestTable[%d] 0x%x\n",__LINE__,sendRequestIndex,sendRequestTable[sendRequestIndex]);
#endif
      if (sendRequestTable[sendRequestIndex] != NULL) {
        mpiresult = MPI_Test(&sendRequestTable[sendRequestIndex],&mpitestflag,&status);

        if ((mpitestflag != 0) && (mpiresult != MPI_SUCCESS)) {
          /* Here a send failed */
          fprintf(consoleHandle,"ERROR: MPI send failed with status %d from process %d\n",mpiresult,sendRequestIndex);
          exit(-1);
        }
      }
    }
    mpiresult = MPI_Irecv(&receive_buffer,           /* message buffer */
                          MAX_MESSAGE_SIZE,                 /* one data item */
                          MPI_BYTE,        /* of type double real */
                          MPI_ANY_SOURCE,    /* receive from any sender */
                          MPI_ANY_TAG,       /* any type of message */
                          MPI_COMM_WORLD,    /* default communicator */
                          &receiveRequest);
    if (mpiresult == MPI_SUCCESS) {
      while (1) {
        memset(&status,0,sizeof(status));
        mpiresult = MPI_Test(&receiveRequest,&mpitestflag,&status);
        if (mpiresult != MPI_SUCCESS) {
          break;
        }
        if (mpitestflag) {
          break;
        }
        sleep(1);
      }
    }
    if ((status.MPI_SOURCE < 0) ||
        (status.MPI_SOURCE >= pFileCommon->maxproc)) {
      fprintf(consoleHandle,"ERROR: unknown source. Master received tag %d %s, err %d, count %3d, cancelled %d mpiresult %d from source %d",
              status.MPI_TAG,tagString[status.MPI_TAG],status.MPI_ERROR,count,status._cancelled,mpiresult,status.MPI_SOURCE);
    } else {
      pWorkEntry = &pWorkQueue[status.MPI_SOURCE];
      

      MPI_Get_count(&status,MPI_BYTE,&count);
      if (mpiresult != MPI_SUCCESS) {

        fprintf(consoleHandle,"ERROR: Master receiving tag %d %s, err %d, count %3d, cancelled %d mpiresult %d from source %d",
                status.MPI_TAG,tagString[status.MPI_TAG],status.MPI_ERROR,count,status._cancelled,mpiresult,status.MPI_SOURCE);
        pWorkEntry->entryState = STATE_ERROR; 
      } else {
        if (pFileCommon->verbose) {
          fprintf(consoleHandle,"Master receiving tag %d %s, err %d, count %3d, cancelled %d mpiresult %d from source %d\n",
                  status.MPI_TAG,tagString[status.MPI_TAG],status.MPI_ERROR,count,status._cancelled,mpiresult,status.MPI_SOURCE);
        }

      }
    }

    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"line %d numpending %d, numstarted %d, numprocessed %d (3)\n",__LINE__,numpending,numstarted,numprocessed);
      for (workIndex = 0; workIndex < pFileCommon->maxproc; ++workIndex) {
        pWorkEntry2 = &pWorkQueue[workIndex];
        fprintf(consoleHandle,"WorkIndex %d entryState %d %s\n",workIndex,pWorkEntry2->entryState,stateString[pWorkEntry2->entryState]);
      }
    }

    switch (status.MPI_TAG) {
    case WORKTAG:
      /* Only the master receives this one */
      if (pWorkEntry->entryState == STATE_WORK_SENT) {
        /* Process the work immediately */
        if (sizeof(WORKQUEUE) == count) {        
          pReceiveWorkEntry = (PWORKQUEUE)&receive_buffer;
          if (pFileCommon->verbose) {
            fprintf(consoleHandle," work item %d from itself\n",pReceiveWorkEntry->queueIndex);
          }
          SearchProcess(consoleHandle,pFileCommon,pReceiveWorkEntry,qualifier);
          pWorkEntry->errorFlag = pReceiveWorkEntry->errorFlag;
          pWorkEntry->consoleFilesize = pReceiveWorkEntry->consoleFilesize;
          pWorkEntry->outFilesize = pReceiveWorkEntry->outFilesize;
          pWorkEntry->databaseFilesize = pReceiveWorkEntry->databaseFilesize;
          pWorkEntry->propermotionFilesize = pReceiveWorkEntry->propermotionFilesize;
          pWorkEntry->plateFilesize = pReceiveWorkEntry->plateFilesize;
          pWorkEntry->idFilesize = pReceiveWorkEntry->idFilesize;
          pWorkEntry->astrometryFilesize = pReceiveWorkEntry->astrometryFilesize;
          pWorkEntry->seconds = pReceiveWorkEntry->seconds;
          pWorkEntry->totalStars = pReceiveWorkEntry->totalStars;
          pWorkEntry->recoveredStars = pReceiveWorkEntry->recoveredStars;
          pWorkEntry->badModulusBin = pReceiveWorkEntry->badModulusBin;
          pWorkEntry->binCount = pReceiveWorkEntry->binCount;
          pWorkEntry->groupNumber = pReceiveWorkEntry->groupNumber;
          pWorkEntry->noGscEntryCount = pReceiveWorkEntry->noGscEntryCount;
          pWorkEntry->staleVersionIdCount = pReceiveWorkEntry->staleVersionIdCount;
          pWorkEntry->curGscBinIndex = pReceiveWorkEntry->curGscBinIndex;
          pWorkEntry->totalProcessedMagnitudes = pReceiveWorkEntry->totalProcessedMagnitudes;
          pWorkEntry->totalNoneMagnitudes = pReceiveWorkEntry->totalNoneMagnitudes;
          pWorkEntry->totalBlendNomatchMagnitudes = pReceiveWorkEntry->totalBlendNomatchMagnitudes;
          pWorkEntry->totalMagnitudes = pReceiveWorkEntry->totalMagnitudes;
          pWorkEntry->queryCount = pReceiveWorkEntry->queryCount;
          pWorkEntry->unmatchedStars = pReceiveWorkEntry->unmatchedStars;
          pWorkEntry->goodUnmatchedStars = pReceiveWorkEntry->goodUnmatchedStars;
          memcpy(&pWorkEntry->counterBlock,&pReceiveWorkEntry->counterBlock,sizeof(FILECOUNTERBLOCK));
          memcpy(&pWorkEntry->vectorCounterBlock,&pReceiveWorkEntry->vectorCounterBlock,sizeof(VECTORCOUNTERBLOCK));

#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pReceiveWorkEntry->catalogDistTotal[25],
             pReceiveWorkEntry->catalogDistSimilarMag[25],
             pReceiveWorkEntry->catalogDistDifferentMag[25]);

#endif
          memcpy(&pWorkEntry->catalogDistTotal,&pReceiveWorkEntry->catalogDistTotal,sizeof(pWorkEntry->catalogDistTotal));
          memcpy(&pWorkEntry->catalogDistSimilarMag,&pReceiveWorkEntry->catalogDistSimilarMag,sizeof(pWorkEntry->catalogDistSimilarMag));
          memcpy(&pWorkEntry->catalogDistDifferentMag,&pReceiveWorkEntry->catalogDistDifferentMag,sizeof(pWorkEntry->catalogDistDifferentMag));
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pWorkEntry->catalogDistTotal[25],
             pWorkEntry->catalogDistSimilarMag[25],
             pWorkEntry->catalogDistDifferentMag[25]);

#endif


          pWorkEntry->entryState = STATE_WORK_RECEIVED;
        } else {
          fprintf(consoleHandle,"ERROR: master expected %d bytes, got %d bytes for entry %d\n",sizeof(WORKQUEUE),count,pWorkEntry->queueIndex);
          pWorkEntry->entryState = STATE_ERROR;
        }
      } else {
        fprintf(consoleHandle,"ERROR: WORKTAG received in wrong state %d %s for entry %d\n",pWorkEntry->entryState,stateString[pWorkEntry->entryState],pWorkEntry->queueIndex);
        pWorkEntry->entryState = STATE_ERROR;
      }
      break;
    case RESPONSETAG:
      if (pWorkEntry->entryState == STATE_WORK_SENT) {
        pWorkEntry->entryState = STATE_WORK_RECEIVED;
        preceive_result = (PWORKQUEUE)&receive_buffer;
#if 0
        printf("Line %d console file for workIndex %d is %s\n",__LINE__,workIndex,preceive_result->consoleLogname);
#endif
        if (pFileCommon->verbose) {
          fprintf(consoleHandle,"Response result: value %3d work item %3d entry %d in state %d %s \n",preceive_result->workResult,preceive_result->queueIndex,pWorkEntry->queueIndex,pWorkEntry->entryState,stateString[pWorkEntry->entryState]);
        }
#if 0
        printf("counterBlock.readMagnitudes %lld at line %d\n",preceive_result->counterBlock.readMagnitudes,__LINE__);
#endif
        /* Now copy important information from the receive buffer */
        pWorkEntry->errorFlag = preceive_result->errorFlag;
        pWorkEntry->consoleFilesize = preceive_result->consoleFilesize;
        pWorkEntry->outFilesize = preceive_result->outFilesize;
        pWorkEntry->databaseFilesize = preceive_result->databaseFilesize;
        pWorkEntry->propermotionFilesize = preceive_result->propermotionFilesize;
        pWorkEntry->plateFilesize = preceive_result->plateFilesize;
        pWorkEntry->idFilesize = preceive_result->idFilesize;
        pWorkEntry->astrometryFilesize = preceive_result->astrometryFilesize;
        pWorkEntry->seconds = preceive_result->seconds;
        memcpy(pWorkEntry->histogramTable,preceive_result->histogramTable,sizeof(preceive_result->histogramTable));
        pWorkEntry->totalStars = preceive_result->totalStars;
        pWorkEntry->recoveredStars = preceive_result->recoveredStars;
        pWorkEntry->badModulusBin = preceive_result->badModulusBin;
        pWorkEntry->binCount = preceive_result->binCount;
        pWorkEntry->groupNumber = preceive_result->groupNumber;
        pWorkEntry->noGscEntryCount = preceive_result->noGscEntryCount;
        pWorkEntry->staleVersionIdCount = preceive_result->staleVersionIdCount;
        pWorkEntry->curGscBinIndex = preceive_result->curGscBinIndex;
        pWorkEntry->totalProcessedMagnitudes = preceive_result->totalProcessedMagnitudes;
        pWorkEntry->totalNoneMagnitudes = preceive_result->totalNoneMagnitudes;
        pWorkEntry->totalBlendNomatchMagnitudes = preceive_result->totalBlendNomatchMagnitudes;
        pWorkEntry->totalMagnitudes = preceive_result->totalMagnitudes;
        pWorkEntry->queryCount = preceive_result->queryCount;
        pWorkEntry->unmatchedStars = preceive_result->unmatchedStars;
        pWorkEntry->goodUnmatchedStars = preceive_result->goodUnmatchedStars;
        memcpy(&pWorkEntry->counterBlock,&preceive_result->counterBlock,sizeof(FILECOUNTERBLOCK));
        memcpy(&pWorkEntry->vectorCounterBlock,&preceive_result->vectorCounterBlock,sizeof(VECTORCOUNTERBLOCK));
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             preceive_result->catalogDistTotal[25],
             preceive_result->catalogDistSimilarMag[25],
             preceive_result->catalogDistDifferentMag[25]);

#endif
        memcpy(&pWorkEntry->catalogDistTotal,&preceive_result->catalogDistTotal,sizeof(pWorkEntry->catalogDistTotal));
        memcpy(&pWorkEntry->catalogDistSimilarMag,&preceive_result->catalogDistSimilarMag,sizeof(pWorkEntry->catalogDistSimilarMag));
        memcpy(&pWorkEntry->catalogDistDifferentMag,&preceive_result->catalogDistDifferentMag,sizeof(pWorkEntry->catalogDistDifferentMag));
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pWorkEntry->catalogDistTotal[25],
             pWorkEntry->catalogDistSimilarMag[25],
             pWorkEntry->catalogDistDifferentMag[25]);

#endif
#if 0
        printf("counterBlock.readMagnitudes %lld in pWorkEntry %d at line %d\n",pWorkEntry->counterBlock.readMagnitudes,pWorkEntry->queueIndex,__LINE__);
#endif

        numpending--;
        numprocessed++;
        time(&curTime);
        curTime -= pFileCommon->beginTime;
        fprintf(consoleHandle,"PARENT: Completed subprocess %3d on %3d at %d seconds for bins %9d to %9d\n",numprocessed,rank,curTime,pWorkEntry->startGscBinIndex,pWorkEntry->endGscBinIndex);
        pFileCommon->completeCount++;
        sendNewWorkRequest = 1;
      } else {
        fprintf(consoleHandle,"ERROR: RESPONSETAG received in wrong state %d %s for entry %d\n",pWorkEntry->entryState,stateString[pWorkEntry->entryState],pWorkEntry->queueIndex);
        pWorkEntry->entryState = STATE_ERROR;
      }
      break;
    case STATUSTAG:
      if (pWorkEntry->entryState == STATE_STATUS_SENT) {
        numpending--;
        if (count == 0) {
          pWorkEntry->entryState = STATE_STATUS_RECEIVED;
          fprintf(consoleHandle,"Master receiving STATUSTAG\n");        
        } else if (sizeof(NODESTATS) == count) {
          pWorkEntry->entryState = STATE_STATUS_RECEIVED;
          pstatusmessage = (PNODESTATS)&receive_buffer;
          sprintf(rxstring,"%lld",pstatusmessage->rxrate);
          sprintf(txstring,"%lld",pstatusmessage->txrate);
          sprintf(rxdiskstring,"%lld",pstatusmessage->rxdiskrate);
          sprintf(txdiskstring,"%lld",pstatusmessage->txdiskrate);

          for (cpuIndex = 0; cpuIndex < pFileCommon->remoteCPUTotal; cpuIndex++) {
            if (strcmp(pstatusmessage->nodename,pFileCommon->remoteProcessorName[cpuIndex]) == 0) {
              break;
            }
          }
          if (cpuIndex < MAX_CPUS) {
            if (cpuIndex == pFileCommon->remoteCPUTotal) {
              strcpy(pFileCommon->remoteProcessorName[cpuIndex],pstatusmessage->nodename);
              pFileCommon->remoteCPUTotal++;
            } 
          } else {
            fprintf(consoleHandle,"ERROR: cpuIndex %d greater than MAX_CPUS\n");
            exit(-1);
          }



          fprintf(consoleHandle,"Recstat %d node %3d %s proc: %4d %4d Totalsec %6d %6d Elapsec %6d %6d nrx %6s ntx %6s kB drx %6s dtx %6s kB/sec\n",
                  pFileCommon->statusRequestNumber,
                  status.MPI_SOURCE,
                  pstatusmessage->nodename,
                  pstatusmessage->taskscompleted,
                  pstatusmessage->newTaskscompleted,
                  pstatusmessage->totalSeconds,
                  pstatusmessage->newTotalSeconds,
                  pstatusmessage->elapsedSeconds,
                  pstatusmessage->newElapsedSeconds,
                  rxstring,
                  txstring,
                  rxdiskstring,
                  txdiskstring);
              

        } else {
          pWorkEntry->entryState = STATE_ERROR;
          fprintf(consoleHandle,"ERROR: STATUSTAG with size %d received for entry %d\n",count,pWorkEntry->queueIndex);
        }
      } else {
        fprintf(consoleHandle,"ERROR: STATUS received in wrong state %d %s for entry  %d\n",pWorkEntry->entryState,stateString[pWorkEntry->entryState],pWorkEntry->queueIndex);
        pWorkEntry->entryState = STATE_ERROR;
      }
      break;
    case DIETAG:
      numpending--;
      if (count == 0) {
        printf("Master receiving DIETAG\n");
        pWorkEntry->entryState = STATE_DIE_RECEIVED;
      } else if (sizeof(NODESTATS) == count) {
        pWorkEntry->entryState = STATE_DIE_RECEIVED;
        pdiemessage = (PNODESTATS)&receive_buffer;
        sprintf(rxstring,"%lld",pdiemessage->rxrate);
        sprintf(txstring,"%lld",pdiemessage->txrate);
        sprintf(rxdiskstring,"%lld",pdiemessage->rxdiskrate);
        sprintf(txdiskstring,"%lld",pdiemessage->txdiskrate);


        
        for (cpuIndex = 0; cpuIndex < pFileCommon->remoteCPUTotal; cpuIndex++) {
          if (strcmp(pdiemessage->nodename,pFileCommon->remoteProcessorName[cpuIndex]) == 0) {
            break;
          }
        }
        if (cpuIndex < MAX_CPUS) {
          if (cpuIndex == pFileCommon->remoteCPUTotal) {
            strcpy(pFileCommon->remoteProcessorName[cpuIndex],pdiemessage->nodename);
            pFileCommon->remoteCPUTotal++;
          } 
        } else {
          printf("ERROR: cpuIndex %d greater than MAX_CPUS\n");
          exit(-1);
        }


        printf("Recstat %d node %3d %s proc: %4d %4d Totalsec %6d %6d Elapsec %6d %6d nrx %6s ntx %6s kB drx %6s dtx %6s kB/sec(DIETAG)\n",
               pFileCommon->statusRequestNumber,
               status.MPI_SOURCE,
               pdiemessage->nodename,
               pdiemessage->taskscompleted,
               pdiemessage->newTaskscompleted,
               pdiemessage->totalSeconds,
               pdiemessage->newTotalSeconds,
               pdiemessage->elapsedSeconds,
               pdiemessage->newElapsedSeconds,
               rxstring,
               txstring,
               rxdiskstring,
               txdiskstring);

      } else {
        pWorkEntry->entryState = STATE_ERROR;
        printf("ERROR: DIETAG with size %d received from entry %d\n",count,pWorkEntry->queueIndex);
      }



      /* We are done, having received our own DIETAG */
      break;

    default:
      fprintf(consoleHandle,"ERROR: received unknown tag %d from entry %d\n",status.MPI_TAG,pWorkEntry->queueIndex);
      break;
    }
    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"line %d numpending %d, numstarted %d, numprocessed %d (3)\n",__LINE__,numpending,numstarted,numprocessed);
      for (workIndex = 0; workIndex < pFileCommon->maxproc; ++workIndex) {
        pWorkEntry2 = &pWorkQueue[workIndex];
        fprintf(consoleHandle,"WorkIndex %d entryState %d %s\n",workIndex,pWorkEntry2->entryState,stateString[pWorkEntry2->entryState]);
      }
    }
    /* Send out any responses if necessary */
    tempCompletedFlag = 1;
    doneFlag = 1;
    for (rank = pFileCommon->startrank; rank < pFileCommon->maxproc; ++rank) {
      pWorkEntry = &pWorkQueue[rank];
      switch (pWorkEntry->entryState) {
      case  STATE_INIT:
        fprintf(consoleHandle,"ERROR: entry %d is in state %d %s\n",pWorkEntry->queueIndex,pWorkEntry->entryState,stateString[pWorkEntry->entryState]);
        pWorkEntry->entryState = STATE_ERROR;
        break;
      case  STATE_STATUS_SENT:
        doneFlag = 0;
        tempCompletedFlag = 0;
        break;
      case  STATE_STATUS_RECEIVED:
        /* Here we can send out our work request */
        mpiresult = MPI_Isend(pWorkEntry,             /* message buffer */
                              sizeof(WORKQUEUE),                 /* one data item */
                              MPI_BYTE,           /* data item is an integer */
                              rank,              /* destination process rank */
                              WORKTAG,           /* user chosen message tag */
                              MPI_COMM_WORLD,   /* default communicator */
                              &sendRequestTable[rank]); 
        if (mpiresult != MPI_SUCCESS) {
          fprintf(consoleHandle,"ERROR: Master send to %d with result %d\n",rank,mpiresult);
          pWorkEntry->entryState = STATE_ERROR;
          
        } else {
          doneFlag = 0;
          tempCompletedFlag = 0;
          pWorkEntry->entryState = STATE_WORK_SENT;
          time(&curTime);
          curTime -= pFileCommon->beginTime;
          fprintf(consoleHandle,"PARENT: Created subprocess %3d on %3d at %d seconds for bins %9d to %9d\n",numstarted,rank,curTime,pWorkEntry->startGscBinIndex,pWorkEntry->endGscBinIndex);
          pWorkEntry++;
          numstarted++;
          numpending++;
        }
        break;
      case  STATE_WORK_SENT:
        doneFlag = 0;
        tempCompletedFlag = 0;
        break;
      case  STATE_WORK_RECEIVED:
        doneFlag = 0;
        break;
      case  STATE_DIE_SENT:
        doneFlag = 0;
        break;
      case  STATE_DIE_RECEIVED:
        break;
      case  STATE_ERROR:
        break;

      }
#if 0
      fprintf(consoleHandle,"Line: %d WorkIndex %d entryState %d %s tempCompletedFlag %d doneFlag %d\n",__LINE__,workIndex,pWorkEntry->entryState,stateString[pWorkEntry->entryState],tempCompletedFlag,doneFlag);
#endif
    }
    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"line %d numpending %d, numstarted %d, numprocessed %d (3)\n",__LINE__,numpending,numstarted,numprocessed);
      for (workIndex = 0; workIndex < pFileCommon->maxproc; ++workIndex) {
        pWorkEntry = &pWorkQueue[workIndex];
        fprintf(consoleHandle,"WorkIndex %d entryState %d %s\n",workIndex,pWorkEntry->entryState,stateString[pWorkEntry->entryState]);
      }
    }
    if (pFileCommon->canceledFlag != 0) {
      if (pFileCommon->verbose) {
        fprintf(consoleHandle,"line %d, canceledFlag is set\n",__LINE__);
      }
      tempCompletedFlag = 1;
    }
    if ((tempCompletedFlag == 1) && (completedFlag == 0)) {
      /* Everything is in a terminal state.  Send out our die tags to everyone regardless of state */
      if (pFileCommon->verbose) {
        fprintf(consoleHandle,"line %d, canceledFlag %d completeFlag %d tempCompletedFlag %d\n",__LINE__,pFileCommon->canceledFlag,completedFlag,tempCompletedFlag);
      }

      completedFlag = 1;
      for (rank = pFileCommon->startrank; rank < pFileCommon->maxproc; ++rank) {
        pWorkEntry = &pWorkQueue[rank];
        mpiresult = MPI_Isend(0, 0, MPI_BYTE, rank, DIETAG, MPI_COMM_WORLD,&sendRequestTable[rank]);
        numpending++;
        if (mpiresult != MPI_SUCCESS) {
          fprintf(consoleHandle,"ERROR: Master sending DIETAG to %d with result %d\n",rank,mpiresult);
          pWorkEntry->entryState = STATE_ERROR;
        } else {
          pWorkEntry->entryState = STATE_DIE_SENT;
          doneFlag = 0;
        }
      }
    }
    

    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"line %d numpending %d, numstarted %d, numprocessed %d (3)\n",__LINE__,numpending,numstarted,numprocessed);
      for (workIndex = 0; workIndex < pFileCommon->maxproc; ++workIndex) {
        pWorkEntry = &pWorkQueue[workIndex];
        fprintf(consoleHandle,"WorkIndex %d entryState %d %s\n",workIndex,pWorkEntry->entryState,stateString[pWorkEntry->entryState]);
      }
    }
    result = stat(pFileCommon->keepalivename,&statbuf);

    if (result != 0) {
      fprintf(consoleHandle,"Terminating at line %d because keepalive file %s is missing\n",numstarted+1,pFileCommon->keepalivename);
      pFileCommon->canceledFlag = 1;
    } else {
      if (pFileCommon->verbose) {
        fprintf(consoleHandle,"Stat on %s is %d\n",pFileCommon->keepalivename,result);
      }
    }

  }




  if (pFileCommon->verbose) {
    fprintf(consoleHandle,"Line %d LOOP EXITED: numpending %d, numstarted %d, numprocessed %d (4)\n",__LINE__,numpending,numstarted,numprocessed);
  }
  InitProperMotions(consoleHandle,pFileCommon,propermotionfile,1);
  if (pFileCommon->maxproc > 1) {
    pWorkEntry0->groupNumber = 1;
  }
  for (workIndex = 0; workIndex < pFileCommon->maxproc; ++workIndex) {
    pWorkEntry = &pWorkQueue[workIndex];
    fprintf(consoleHandle,"WorkIndex %d entryState %d %s\n",workIndex,pWorkEntry->entryState,stateString[pWorkEntry->entryState]);
    if (workIndex >= pFileCommon->startrank) {
      /* Open and print out the console buffer */
      ReadLogFile(consoleHandle,pWorkEntry);
      /* Open and append the log file */
      ReadResultFile(consoleHandle,outHandle,pWorkEntry->queueIndex,pWorkEntry->outLogname,pWorkEntry->errorFlag,pWorkEntry->outFilesize);
      /* Open and append the database file */
      ReadResultFile(consoleHandle,databaseHandle,pWorkEntry->queueIndex,pWorkEntry->databaseLogname,pWorkEntry->errorFlag,pWorkEntry->databaseFilesize);
      /* Open and append the id file */
      ReadResultFile(consoleHandle,idHandle,pWorkEntry->queueIndex,pWorkEntry->idLogname,pWorkEntry->errorFlag,pWorkEntry->idFilesize);
      /* Open and append the propermotion file */
      ReadResultFile(consoleHandle,pFileCommon->propermotionHandle,pWorkEntry->queueIndex,pWorkEntry->propermotionLogname,pWorkEntry->errorFlag,pWorkEntry->propermotionFilesize);
      /* Open and merge the plate outlier file */
      MergePlateOutlierFile(consoleHandle,pFileCommon,pVectorStore,plateHandle,pWorkEntry->queueIndex,pWorkEntry->plateLogname,pWorkEntry->errorFlag,pWorkEntry->plateFilesize);
      /* Open and merge the drad histogram file */
      MergeDradHistogramFile(consoleHandle,pFileCommon,pVectorStore,astrometryHandle,pWorkEntry->queueIndex,pWorkEntry->astrometryLogname,pWorkEntry->errorFlag,pWorkEntry->astrometryFilesize);
    }
    if (workIndex > 0) {
      for (index = 0; index < MAX_NONE_HISTOGRAM; index++) {
        pFileCommon->histogramTable[index] += pWorkEntry->histogramTable[index];
       
      }
#if 0
      printf("counterBlock.readMagnitudes %lld in pWorkEntry %d at line %d\n",pWorkEntry->counterBlock.readMagnitudes,pWorkEntry->queueIndex,__LINE__);
#endif
      pWorkEntry0->totalStars += pWorkEntry->totalStars;
      pWorkEntry0->recoveredStars += pWorkEntry->recoveredStars;
      pWorkEntry0->badModulusBin += pWorkEntry->badModulusBin;
      pWorkEntry0->binCount += pWorkEntry->binCount;
      pWorkEntry0->groupNumber += pWorkEntry->groupNumber -1;
      pWorkEntry0->noGscEntryCount += pWorkEntry->noGscEntryCount;
      pWorkEntry0->staleVersionIdCount += pWorkEntry->staleVersionIdCount;
      if (pWorkEntry->curGscBinIndex > pWorkEntry0->curGscBinIndex) {
        pWorkEntry0->curGscBinIndex = pWorkEntry->curGscBinIndex;
      }
      pWorkEntry0->totalProcessedMagnitudes += pWorkEntry->totalProcessedMagnitudes;
      pWorkEntry0->totalNoneMagnitudes += pWorkEntry->totalNoneMagnitudes;
      pWorkEntry0->totalBlendNomatchMagnitudes += pWorkEntry->totalBlendNomatchMagnitudes;
      pWorkEntry0->totalMagnitudes += pWorkEntry->totalMagnitudes;
      pWorkEntry0->queryCount += pWorkEntry->queryCount;
      pWorkEntry0->unmatchedStars += pWorkEntry->unmatchedStars;
      pWorkEntry0->goodUnmatchedStars += pWorkEntry->goodUnmatchedStars;

      UpdateCounterBlock(&pFileCommon->counterBlock,&pWorkEntry->counterBlock);
      UpdateVectorCounterBlock(&pWorkEntry0->vectorCounterBlock,&pWorkEntry->vectorCounterBlock);
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pWorkEntry->catalogDistTotal[25],
             pWorkEntry->catalogDistSimilarMag[25],
             pWorkEntry->catalogDistDifferentMag[25]);

#endif
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pWorkEntry0->catalogDistTotal[25],
             pWorkEntry0->catalogDistSimilarMag[25],
             pWorkEntry0->catalogDistDifferentMag[25]);

#endif
      UpdateCatalogDist(pWorkEntry0->catalogDistTotal,pWorkEntry->catalogDistTotal);
      UpdateCatalogDist(pWorkEntry0->catalogDistSimilarMag,pWorkEntry->catalogDistSimilarMag);
      UpdateCatalogDist(pWorkEntry0->catalogDistDifferentMag,pWorkEntry->catalogDistDifferentMag);
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pWorkEntry0->catalogDistTotal[25],
             pWorkEntry0->catalogDistSimilarMag[25],
             pWorkEntry0->catalogDistDifferentMag[25]);

#endif
#if 0
    printf("counterBlock.readMagnitudes %lld in pWorkEntry %d at line %d\n",pWorkEntry0->counterBlock.readMagnitudes,pWorkEntry0->queueIndex,__LINE__);
#endif

    }

  }
  
  if (pFileCommon->enableTransientSearch > 0) {
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pWorkEntry0->catalogDistTotal[25],
             pWorkEntry0->catalogDistSimilarMag[25],
             pWorkEntry0->catalogDistDifferentMag[25]);

#endif

    DumpCatalogDistanceHistogram(pWorkEntry0->catalogDistTotal,pWorkEntry0->catalogDistSimilarMag,pWorkEntry0->catalogDistDifferentMag);
  }


  pFileCommon->queueIteration++;
  time(&curTime);
  curTime -= startTime;
  if (pFileCommon->completeCount == 0) {
    averageTime = 0.0;
  } else {
    averageTime = (1.0*curTime)/pFileCommon->completeCount;
  }
  fprintf(consoleHandle,"Phase %d Time: %d seconds (%.1f sec/plate) for %d mosaics completed out of %d mosaics %d processors\n",
          pFileCommon->queueIteration,
          curTime,
          averageTime,
          pFileCommon->completeCount,
          pFileCommon->workQueueSize,
          pFileCommon->maxproc);

  for (seriesIndex = 0; seriesIndex <= MAX_SERIES; seriesIndex++) {
    plateCountTable = pVectorStore->plateCountTable[seriesIndex];

    if (plateCountTable == NULL) {
      continue;
    }
    plateCountTable = pVectorStore->plateCountTable[seriesIndex];
    for (plateIndex = 0; plateIndex <= pFileCommon->maxPlateNumber[seriesIndex]; plateIndex++) {
      if ((plateCountTable[plateIndex].measurementCount != 0)  ||
          (plateCountTable[plateIndex].matchCount != 0)  ||
          (plateCountTable[plateIndex].lightcurveCount != 0)  ||
          (plateCountTable[plateIndex].outlierCount != 0)) {
        double outlierRatio = 0.0;
        double matchRatio = 0.0;
        if (plateCountTable[plateIndex].lightcurveCount > 0) {
          outlierRatio = (1.0*plateCountTable[plateIndex].outlierCount)/(1.0*plateCountTable[plateIndex].lightcurveCount);
        }
        if (plateCountTable[plateIndex].measurementCount > 0) {
          matchRatio = (1.0*plateCountTable[plateIndex].matchCount)/(1.0*plateCountTable[plateIndex].measurementCount);
        }
        fprintf(plateHandle,"%s%05d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%f\n",
                GetSeriesString(seriesIndex,1),
                plateIndex,
                seriesIndex,
                plateIndex,
                plateCountTable[plateIndex].measurementCount,
                plateCountTable[plateIndex].matchCount,
                matchRatio,
                plateCountTable[plateIndex].lightcurveCount,
                plateCountTable[plateIndex].outlierCount,
                outlierRatio);
                
      }
    }

  }

  if (astrometryHandle != NULL) {
    int magindex;
    int dradindex;
    int dradvalue;
    for (magindex = 0; magindex < MAX_DRAD_HISTOGRAM_MAG; magindex++) {
      for (dradindex = 0; dradindex < (2*(MAX_DRAD_RANGE)+1); dradindex++) {
        dradvalue = dradindex-MAX_DRAD_RANGE;
        if ((pFileCommon->draHistogram[magindex][dradindex] != 0) ||
            (pFileCommon->ddecHistogram[magindex][dradindex] != 0)) {
          fprintf(astrometryHandle,"%d\t%d\t%lld\t%lld\n",magindex,dradvalue,pFileCommon->draHistogram[magindex][dradindex],pFileCommon->ddecHistogram[magindex][dradindex]);
        }

      }

    }


    fclose(astrometryHandle);
    astrometryHandle = NULL;
  }
  pWorkEntry = &pWorkQueue[0];

  time(&curTime);
  curTime -= pFileCommon->startTime;

#if 0
  printf("totalMagnitudes %d in pWorkEntry %d at line %d\n",pWorkEntry->totalMagnitudes,pWorkEntry->queueIndex,__LINE__);
#endif
  fprintf(consoleHandle,"Execution Time: %d seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d gscBins, %d total magnitudes\n",curTime,pWorkEntry->queryCount,pWorkEntry->totalStars,pWorkEntry->recoveredStars,pWorkEntry->badModulusBin,pWorkEntry->binCount,pWorkEntry->totalMagnitudes);

  fprintf(outHandle,"Execution Time: %d seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d gscBins, %d total magnitudes\n",curTime,pWorkEntry->queryCount,pWorkEntry->totalStars,pWorkEntry->recoveredStars,pWorkEntry->badModulusBin,pWorkEntry->binCount,pWorkEntry->totalMagnitudes);


  fprintf(consoleHandle," queries/sec: %f totalStars/sec: %f recoveredStars+badModulusBin.sec: %f gscBins/sec %f\n",
          (1.0 * pWorkEntry->queryCount)/(1.0*curTime),
          (1.0 * pWorkEntry->totalStars)    /(1.0*curTime),
          (1.0 * (pWorkEntry->recoveredStars + pWorkEntry->badModulusBin) )    /(1.0*curTime),
          (1.0 * pWorkEntry->binCount  )    /(1.0*curTime));
  fprintf(outHandle," queries/sec: %f totalStars/sec: %f recoveredStars+badModulusBin.sec: %f gscBins/sec %f\n",
          (1.0 * pWorkEntry->queryCount)/(1.0*curTime),
          (1.0 * pWorkEntry->totalStars)    /(1.0*curTime),
          (1.0 * (pWorkEntry->recoveredStars+pWorkEntry->badModulusBin) )    /(1.0*curTime),
          (1.0 * pWorkEntry->binCount  )    /(1.0*curTime));
 


  fprintf(consoleHandle,"... %d noGscEntry\n",pWorkEntry->noGscEntryCount);
  fprintf(consoleHandle,"... %d staleVersionId\n",pWorkEntry->staleVersionIdCount);
  fprintf(outHandle,"... %d noGscEntry\n",pWorkEntry->noGscEntryCount);
  fprintf(outHandle,"... %d staleVersionId\n",pWorkEntry->staleVersionIdCount);

  fprintf(consoleHandle,"New curGscBinIndex %d\n",pWorkEntry->curGscBinIndex);
  fprintf(outHandle,"New curGscBinIndex %d\n",pWorkEntry->curGscBinIndex);
#if 0
  printf("counterBlock.readMagnitudes %lld in pWorkEntry %d at line %d\n",pWorkEntry->counterBlock.readMagnitudes,pWorkEntry->queueIndex,__LINE__);
#endif

  fprintf(consoleHandle,"Total Magnitudes %d, Total NONE magnitudes %d Total NOMATCH magnitudes %d Total GSC nonblend magnitudes %d \n",
          pWorkEntry->totalProcessedMagnitudes,pWorkEntry->totalNoneMagnitudes,pWorkEntry->totalBlendNomatchMagnitudes,pWorkEntry->totalProcessedMagnitudes - pWorkEntry->totalNoneMagnitudes - pWorkEntry->totalBlendNomatchMagnitudes);
  fprintf(outHandle,"Total Magnitudes %d, Total NONE magnitudes %d Total NOMATCH magnitudes %d Total GSC nonblend magnitudes %d \n",
          pWorkEntry->totalProcessedMagnitudes,pWorkEntry->totalNoneMagnitudes,pWorkEntry->totalBlendNomatchMagnitudes,pWorkEntry->totalProcessedMagnitudes - pWorkEntry->totalNoneMagnitudes - pWorkEntry->totalBlendNomatchMagnitudes);
  fprintf(consoleHandle,"Total groups %d\n",pWorkEntry->groupNumber-1);
  fprintf(outHandle,"Total groups %d\n",pWorkEntry->groupNumber-1);
  fprintf(consoleHandle,"Max vectorAlloc %d, lightcurveCount %d lightcurvePointCount %d, outlierCount %d, qualityRejectCount %d\n",
          pWorkEntry->vectorCounterBlock.vectorAlloc,
          pWorkEntry->vectorCounterBlock.lightcurveCount,
          pWorkEntry->vectorCounterBlock.lightcurvePointCount,
          pWorkEntry->vectorCounterBlock.outlierCount,
          pWorkEntry->vectorCounterBlock.qualityRejectCount);
  fprintf(consoleHandle,"points read %d %d, pointsRejected %d, insufficientPoints %d, highRMSPoints %d\n",
          pWorkEntry->vectorCounterBlock.pointsRead,
          pWorkEntry->vectorCounterBlock.pointsRead2,
          pWorkEntry->vectorCounterBlock.pointsRejected,
          pWorkEntry->vectorCounterBlock.insufficientPoints,
          pWorkEntry->vectorCounterBlock.highRMSPoints);
  fprintf(consoleHandle,"Unmatch stars %d good unmatched stars %d\n",pWorkEntry->unmatchedStars,pWorkEntry->goodUnmatchedStars);
  fprintf(consoleHandle,"Total Stars %10lld; Total Stars with Outliers %10lld stars plotted %10lld none candidates %10lld\n",
          pWorkEntry->vectorCounterBlock.totalStars,
          pWorkEntry->vectorCounterBlock.totalOutliers,
          pWorkEntry->vectorCounterBlock.totalPlot2,
          pWorkEntry->vectorCounterBlock.totalNoneCandidates);
  OutputRootCounts(consoleHandle,pFileCommon);
  for (dradIndex = 0; dradIndex <= MAX_GROUP_DRAD_HISTOGRAM; dradIndex++) {
    if (pFileCommon->counterBlock.maxGroupHistogram[dradIndex] != 0) {
      fprintf(consoleHandle,"Group Drad: %3d arcsec count: %8d\n",dradIndex,pFileCommon->counterBlock.maxGroupHistogram[dradIndex]);
    }
  }
  for (index = 0; index < MAX_NONE_HISTOGRAM; index++) {
    if (pFileCommon->histogramTable[index] != 0) {
      fprintf(consoleHandle,"Images/Location: %4d  Locations: %4d\n",index,pFileCommon->histogramTable[index]);
      fprintf(outHandle,"Images/Location %4d  Locations: %6d\n",index,pFileCommon->histogramTable[index]);
    }
  }



  FreeVectorStore(pVectorStore);


  return(0);



}

int slave(FILE *consoleHandle,PFILECOMMON pFileCommon,char *qualifier) {
  WORKQUEUE workItem;
  PWORKQUEUE pWorkEntry = &workItem;

  MPI_Status status;
  int mpitestflag;
  MPI_Request receiveRequest;
  PWORKQUEUE results;
  int mpiresult;
  int count;
  struct stat filestats;
  int statResult;
  NODESTATS diemessage;
  time_t startTime;
  time_t newStartTime;
  time_t curTime;
  int totalSeconds = 0;
  int newTotalSeconds = 0;
  long long oldrxbytes;
  long long oldtxbytes;
  long long rxbytes;
  long long txbytes;
  long long denominator;

  long long oldrxdiskbytes;
  long long oldtxdiskbytes;
  long long rxdiskbytes;
  long long txdiskbytes;

  time(&startTime);
  time(&newStartTime);
  GetNetworkCounts(consoleHandle,&oldrxbytes,&oldtxbytes);
  GetDiskCounts(consoleHandle,&oldrxdiskbytes,&oldtxdiskbytes);
#if 0
  fprintf(consoleHandle,"GetNetworkCounts returned rx %lld\n",oldrxbytes);
  fprintf(consoleHandle,"GetNetworkCounts returned tx %lld\n",oldtxbytes);
  sleep(10);
#endif
  memset(&diemessage,0,sizeof(NODESTATS));
  strcpy(diemessage.nodename,pFileCommon->processorName);
  if (pFileCommon->verbose) {
    fprintf(consoleHandle,"Slave %d running on %s\n",pFileCommon->myrank,diemessage.nodename);
  }
    


  while (1) {

    /* Receive a message from the master */

    sleep(1);
    mpiresult = MPI_Irecv(pWorkEntry,           /* message buffer */
                          sizeof(WORKQUEUE),                 /* one data item */
                          MPI_BYTE,        /* of type double real */
                          0,    /* receive from any sender */
                          MPI_ANY_TAG,       /* any type of message */
                          MPI_COMM_WORLD,    /* default communicator */
                          &receiveRequest);
    if (mpiresult == MPI_SUCCESS) {
      while (1) {
        memset(&status,0,sizeof(status));
        mpiresult = MPI_Test(&receiveRequest,&mpitestflag,&status);
        if (mpiresult != MPI_SUCCESS) {
          fprintf(consoleHandle,"ERROR: slave receive failure %d in line %d\n",mpiresult,__LINE__);
          break;
        }
        if (mpitestflag) {
          break;
        }
        sleep(1);
      }
    }
    MPI_Get_count(&status,MPI_BYTE,&count);

    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"Slave receiving tag %d %s, err %d, count %3d, cancelled %d mpiresult %d from source %d\n",
              status.MPI_TAG,tagString[status.MPI_TAG],status.MPI_ERROR,count,status._cancelled,mpiresult,status.MPI_SOURCE);
          
    }

    if (status.MPI_TAG == DIETAG) {
      if (pFileCommon->verbose) {
        fprintf(consoleHandle,"received DIETAG: Slave %d\n",pFileCommon->myrank);
      }
      time(&curTime);
      diemessage.elapsedSeconds = curTime - startTime;
      diemessage.totalSeconds = totalSeconds;
      diemessage.newElapsedSeconds = curTime - newStartTime;
      diemessage.newTotalSeconds = newTotalSeconds;
      GetNetworkCounts(consoleHandle,&rxbytes,&txbytes);
      GetDiskCounts(consoleHandle,&rxdiskbytes,&txdiskbytes);
      if (diemessage.newElapsedSeconds == 0) {
        diemessage.rxrate = 0;
        diemessage.txrate = 0;
        diemessage.rxdiskrate = 0;
        diemessage.txdiskrate = 0;
      } else {
        denominator = NETWORK_MULTIPLIER * diemessage.newElapsedSeconds;
        
        diemessage.rxrate = (rxbytes-oldrxbytes)/denominator;
        diemessage.txrate = (txbytes-oldtxbytes)/denominator; 
        diemessage.rxdiskrate = (rxdiskbytes-oldrxdiskbytes)/denominator;
        diemessage.txdiskrate = (txdiskbytes-oldtxdiskbytes)/denominator; 
      }
#if 0
      fprintf(consoleHandle,"GetNetworkCounts returned rx %lld\n",rxbytes);
      fprintf(consoleHandle,"GetNetworkCounts returned tx %lld\n",txbytes);
      fprintf(consoleHandle,"GetNetworkCounts elapsed %d\n",diemessage.newElapsedSeconds);
      fprintf(consoleHandle,"diemessage.rxrate %lld\n",diemessage.rxrate);
      fprintf(consoleHandle,"diemessage.txrate %lld\n",diemessage.txrate);
      fprintf(consoleHandle,"oldrxbytes %lld\n",oldrxbytes);
      fprintf(consoleHandle,"oldtxbytes %lld\n",oldtxbytes);
      fprintf(consoleHandle,"numerator %lld\n",(rxbytes-oldrxbytes));
      fprintf(consoleHandle,"denominator %lld\n",denominator);
     
#endif
      oldrxbytes = rxbytes;
      oldtxbytes = txbytes;
      oldrxdiskbytes = rxdiskbytes;
      oldtxdiskbytes = txdiskbytes;
      fflush(consoleHandle);
      mpiresult = MPI_Send(&diemessage,sizeof(diemessage), MPI_BYTE, 0, DIETAG, MPI_COMM_WORLD);
      if (mpiresult != MPI_SUCCESS) {
        fprintf(consoleHandle,"ERROR: Slave returning DIETAG from %d with result %d\n",pFileCommon->myrank,mpiresult);
        fflush(consoleHandle);
      }
 
      return(0);
    } else if (status.MPI_TAG == STATUSTAG) {
      if (pFileCommon->verbose) {
        fprintf(consoleHandle,"received STATUSTAG: Slave %d\n",pFileCommon->myrank);
      }
      time(&curTime);
      diemessage.elapsedSeconds = curTime - startTime;
      diemessage.newElapsedSeconds = curTime - newStartTime;
      diemessage.totalSeconds = totalSeconds;
      diemessage.newTotalSeconds = newTotalSeconds;
      GetNetworkCounts(consoleHandle,&rxbytes,&txbytes);
      GetDiskCounts(consoleHandle,&rxdiskbytes,&txdiskbytes);
      if (diemessage.newElapsedSeconds == 0) {
        diemessage.rxrate = 0;
        diemessage.txrate = 0;
        diemessage.rxdiskrate = 0;
        diemessage.txdiskrate = 0;
      } else {
        denominator = NETWORK_MULTIPLIER * diemessage.newElapsedSeconds;
        diemessage.rxrate = (rxbytes-oldrxbytes)/denominator;
        diemessage.txrate = (txbytes-oldtxbytes)/denominator;
        diemessage.rxdiskrate = (rxdiskbytes-oldrxdiskbytes)/denominator;
        diemessage.txdiskrate = (txdiskbytes-oldtxdiskbytes)/denominator; 

      }
#if 0
      fprintf(consoleHandle,"GetNetworkCounts returned rx %lld\n",rxbytes);
      fprintf(consoleHandle,"GetNetworkCounts returned tx %lld\n",txbytes);
      fprintf(consoleHandle,"GetNetworkCounts elapsed %d\n",diemessage.newElapsedSeconds);
      fprintf(consoleHandle,"diemessage.rxrate %lld\n",diemessage.rxrate);
      fprintf(consoleHandle,"diemessage.txrate %lld\n",diemessage.txrate);
#endif
      oldrxbytes = rxbytes;
      oldtxbytes = txbytes;
      oldrxdiskbytes = rxdiskbytes;
      oldtxdiskbytes = txdiskbytes;
      mpiresult = MPI_Send(&diemessage,sizeof(diemessage), MPI_BYTE, 0, STATUSTAG, MPI_COMM_WORLD);
      if (mpiresult != MPI_SUCCESS) {
        fprintf(consoleHandle,"ERROR: Slave returning STATUSTAG from %d with result %d\n",pFileCommon->myrank,mpiresult);
      }
      time(&newStartTime);
      diemessage.newTaskscompleted = 0;
      newTotalSeconds = 0;
      
      

    }  else if (status.MPI_TAG == WORKTAG) {
      diemessage.taskscompleted++;
      diemessage.newTaskscompleted++;
      SearchProcess(consoleHandle,pFileCommon,pWorkEntry,qualifier);
      fflush(consoleHandle);
      fclose(consoleHandle);
      consoleHandle = stdout;
      statResult = stat(pWorkEntry->consoleLogname,&filestats);
      if (statResult != 0) {
        fprintf(consoleHandle,"ERROR: failed to obtain stats on %s\n",pWorkEntry->consoleLogname);
      } else {
        pWorkEntry->consoleFilesize = filestats.st_size;
      }
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pFileCommon->catalogDistTotal[25],
             pFileCommon->catalogDistSimilarMag[25],
             pFileCommon->catalogDistDifferentMag[25]);

      if (pFileCommon->enableTransientSearch > 0) {
        DumpCatalogDistanceHistogram(pFileCommon->catalogDistTotal,pFileCommon->catalogDistSimilarMag,pFileCommon->catalogDistDifferentMag);
      }
#endif
      memcpy(&pWorkEntry->counterBlock,&pFileCommon->counterBlock,sizeof(FILECOUNTERBLOCK));
      memcpy(&pWorkEntry->catalogDistTotal,&pFileCommon->catalogDistTotal,sizeof(pWorkEntry->catalogDistTotal));
      memcpy(&pWorkEntry->catalogDistSimilarMag,&pFileCommon->catalogDistSimilarMag,sizeof(pWorkEntry->catalogDistSimilarMag));
      memcpy(&pWorkEntry->catalogDistDifferentMag,&pFileCommon->catalogDistDifferentMag,sizeof(pWorkEntry->catalogDistDifferentMag));
#if 0 /* Debug only */
      printf("CatalogDistance at line %d for index 25 Total %d Similar %d Different %d\n",
             __LINE__,
             pWorkEntry->catalogDistTotal[25],
             pWorkEntry->catalogDistSimilarMag[25],
             pWorkEntry->catalogDistDifferentMag[25]);

#endif

      mpiresult = MPI_Send(pWorkEntry, sizeof(WORKQUEUE), MPI_BYTE, 0, RESPONSETAG, MPI_COMM_WORLD);
      if (mpiresult != MPI_SUCCESS) {
        fprintf(consoleHandle,"ERROR: sending work item %3d with result %d workResult %8.3f Slave %d \n",pWorkEntry->queueIndex,mpiresult,pWorkEntry->workResult,pFileCommon->myrank);
      }
   

    } else {
      fprintf(consoleHandle,"ERROR: Slave  received unknown tag %d\n",status.MPI_TAG);
      exit(-1);
    }
  }

}


int SearchProcess(FILE *consoleHandle,PFILECOMMON pFileCommon,PWORKQUEUE pWorkEntry,char *qualifier) {
  int curMagnitudes = 0;
  PPHOTSTARIMAGE pMagnitudeTable;
  int magnitudeIndex;
  int keepaliveFlag = 1;
  double lat; /* dec */
  double lon; /* ra  */
  int keepaliveCheckModulus = 0;
  int gotAnswer;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  PMATCHEDCOMMON pMatchedCommon = &pWorkEntry->matchedCommon;
  time_t curTime;
  char timestr[100];
  struct tm *ptr;
  int allocMagnitudes = 0;
  int magnitudeCount;
  int index;
  PPHOTSTARIMAGE pCurStarImage;
  PFILESTARIMAGE pNoneMagnitudeTable = NULL;
  PFILESTARIMAGE tmpNoneMagnitudeTable;
  PFILESTARIMAGE pFileStarImage = NULL;
  PFILESTARIMAGE pNoneImage;
  int noneMagnitudeAlloc = 0;
  PPLATECOUNT plateCountTable;
  int chn_n;
  double p_row;
  double p_coln;
  PPHOTTARGET target_table = NULL;
  PPHOTTARGET pTarget;
  size_t target_nrecs = 0;
  size_t target_alloc = 0;
  size_t target_index;
  long long curPrintMagnitudes = 0;
  int filePrintModulus = 0;
  VECTORSTORE vectorStore;
  PVECTORSTORE pVectorStore = &vectorStore;
  PSKY2KSTAR Sky2KTable = NULL;
  int sky2k_size = 0;
  double maxGSCRadius = 0.0;
  PPHOTSTARIMAGE pGroupImageTable = NULL;
  PPHOTSTARIMAGE pGroupImage;
  int group_image_alloc = 0;
  int group_image_count;
  STARENTRY curStarEntry;
  PSTARENTRY pCurStarEntry = &curStarEntry;
  int refType;
  PARAMETERVECTORSTORE parameterVectorStore;
  PPARAMETERVECTORSTORE pParameterVectorStore = &parameterVectorStore;
  int minGoodStars = MIN_SUMMARY_GOODSTARS;
  int debugMode = 0;
  int excludeSeriesCount = 0; /* Exclude nothing */
  int curMagnitudeIndex;
  int seriesIndex;
  int plateIndex;
  int dradIndex;
  FILE *outHandle = NULL;
  FILE *databaseHandle = NULL;
  FILE *astrometryHandle = NULL;
  FILE *plateHandle = NULL;
  FILE *idHandle = NULL;
  int errorFlag = 0;
  int statResult;
  struct stat filestats;
  int gscBinIndexCurrent = -1;
  PSTARENTRY pStarTable = NULL;
  File catalogHandle = NULL;
  File indexHandle = NULL;
  PGALAXYCOMMON pGalaxyCommon = &pFileCommon->galaxyCommon;
  PFLARECANDIDATE pFlareCandidate = NULL;

  pWorkEntry->totalStars = 0;
  pWorkEntry->recoveredStars = 0;
  pWorkEntry->badModulusBin = 0;
  pWorkEntry->groupNumber = 1;
  pWorkEntry->noGscEntryCount = 0;
  pWorkEntry->staleVersionIdCount = 0;
  pWorkEntry->curGscBinIndex = pWorkEntry->startGscBinIndex;
  pWorkEntry->totalProcessedMagnitudes = 0;
  pWorkEntry->totalNoneMagnitudes = 0;
  pWorkEntry->totalBlendNomatchMagnitudes = 0;
  pWorkEntry->totalMagnitudes = 0;
  pWorkEntry->queryCount = 0;
  pWorkEntry->unmatchedStars = 0;
  pWorkEntry->goodUnmatchedStars = 0;

  fprintf(consoleHandle,"CHILD: SearchProcess starting with rank %d queue index %d start gsc bin %d end gsc bin %d console handle 0x%x console file %s\n",
          pFileCommon->myrank,
          pWorkEntry->queueIndex,
          pWorkEntry->startGscBinIndex,
          pWorkEntry->endGscBinIndex,
          consoleHandle,
          pWorkEntry->consoleLogname);

  memset(pParameterVectorStore,0,sizeof(PARAMETERVECTORSTORE));
  memset(pVectorStore,0,sizeof(VECTORSTORE));
  curMagnitudes = 0;
  pMagnitudeTable = NULL;
  magnitudeIndex = 0;
  curMagnitudes = -1;
  noneMagnitudeAlloc = 100;
  pNoneMagnitudeTable = (PFILESTARIMAGE)calloc(noneMagnitudeAlloc,sizeof(FILESTARIMAGE));
  target_alloc = 100;
  target_table = (PPHOTTARGET)calloc(target_alloc,sizeof(PHOTTARGET));
  if (target_table == NULL) {
    fprintf(consoleHandle,"ERROR: failed to allocate target_table of size %d\n",target_alloc);
    exit(-1);
  }
  if (pFileCommon->enableTransientSearch > 0) {
    pStarTable = (PSTARENTRY) calloc(STAR_LIMIT,sizeof(STARENTRY));
    if (pStarTable == NULL) {
      fprintf(consoleHandle,"ERROR: failed to allocate pStarTable\n");
      exit(-1);
    }


    catalogHandle = Open(catalogname,"r");
    if (catalogHandle == NULL) {
      printf("ERROR Could not open catalog file %s\n",catalogname);
      errorFlag = 1;
    } 
    indexHandle = Open(indexname,"r");
    if (indexHandle == NULL) {
      printf("ERROR Could not open catalog index file %s\n",indexname);
      errorFlag = 1;
    } 
#ifdef  PLOT_LIMITING_MAGNITUDES
		if (OpenGalaxyFiles(pGalaxyCommon,qualifier,NULL) != 0) {
      printf("ERROR Could not open galaxy files\n");
      errorFlag = 1;			
		}
#endif /* PLOT_LIMITING_MAGNITUDES */
  }
  InitProperMotions(consoleHandle,pFileCommon,pWorkEntry->propermotionLogname,0);

  outHandle = fopen(pWorkEntry->outLogname,"a+t");
  if (outHandle == NULL) {
    errorFlag = 1;
    fprintf(consoleHandle,"ERROR: Failed to open the output file %s\n",pWorkEntry->outLogname);
  } else {
    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"Output file %s\n",pWorkEntry->outLogname);
    }
  }


  databaseHandle = fopen(pWorkEntry->databaseLogname,"wt");
  if (databaseHandle == NULL) {
    errorFlag = 1;
    fprintf(consoleHandle,"ERROR: Failed to open the output file %s\n",pWorkEntry->databaseLogname);
  } else {
    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"Database file %s\n",pWorkEntry->databaseLogname);
    }
  }

  astrometryHandle = fopen(pWorkEntry->astrometryLogname,"wt");
  if (astrometryHandle == NULL) {
    errorFlag = 1;
    fprintf(consoleHandle,"ERROR: Failed to open the output file %s\n",pWorkEntry->astrometryLogname);
  } else {
    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"Astrometry file %s\n",pWorkEntry->astrometryLogname);
    }
    fprintf(astrometryHandle,"minmag\tdrad\tdra\tddec\n");
    fprintf(astrometryHandle,"------\t----\t---\t----\n");
  }

  plateHandle = fopen(pWorkEntry->plateLogname,"wt");
  if (plateHandle == NULL) {
    errorFlag = 1;
    fprintf(consoleHandle,"ERROR: Failed to open the output file %s line %d\n",pWorkEntry->plateLogname,__LINE__);
  } else {
    if (pFileCommon->verbose) {
      fprintf(consoleHandle,"Plate file %s\n",pWorkEntry->plateLogname);
    }
    fprintf(plateHandle,"Plate\tseriesId\tplateNumber\tmeasurementCount\tmatchCount\tmatchRatio\tlightcurveCount\toutlierCount\toutlierRatio\n");
    fprintf(plateHandle,"-----\t--------\t-----------\t----------------\t----------\t----------\t---------------\t------------\t------------\n");
  }

  if (pWorkEntry->idLogname[0] != 0) {
    idHandle = fopen(pWorkEntry->idLogname,"wt");
    if (idHandle == NULL) {
      errorFlag = 1;
      fprintf(consoleHandle,"ERROR: Failed to open the output file %s\n",pWorkEntry->idLogname);
    } else {
      if (pFileCommon->verbose) {
        fprintf(consoleHandle,"Id file %s\n",pWorkEntry->idLogname);
      }
    }

  }

  if (errorFlag != 0) {
    return(0);
  }

  if (GetPhotometryGlobal(pFileCommon->pPhotConnection,pPhotGlobal) != 1) {
    fprintf(consoleHandle,"ERROR: failed to get the global photometry table\n");
    exit(-1);
  }
  pMatchedCommon->versionId = pPhotGlobal->currentVersion;
  pMatchedCommon->oldTime = pFileCommon->startTime;

  GetSky2kCatalog(consoleHandle,&Sky2KTable,&sky2k_size,&maxGSCRadius);

  if (pFileCommon->gscBinIndexCount > 0) {
    gscBinIndexCurrent = 0;
  }


  /* Enter the main loop.  We could process each magnitude file individually before going on to the
     next, but then the performance would be extremely slow because we lose the locality of 
     reference into the GSC catalog file.  This more complicated code ensures that we read
     the GSC catalog file in sequence
  */

  while (1) {
    /* First see if we need to read in any data */
    if (keepaliveFlag == 0) {
      break;
    }
    
    if (gscBinIndexCurrent >= 0) {
      /* Select the next gsc bin index from the list */
      if (pFileCommon->flareCandidatesOnly == 0){
        while ((gscBinIndexCurrent < pFileCommon->gscBinIndexCount) &&
               (pFileCommon->gscBinIndexList[gscBinIndexCurrent] < pWorkEntry->curGscBinIndex)) {
          gscBinIndexCurrent++;
        }
        if (gscBinIndexCurrent >= pFileCommon->gscBinIndexCount) {
          pWorkEntry->curGscBinIndex = pGscBin->total_gsc_bins;
        } else {
          pWorkEntry->curGscBinIndex = pFileCommon->gscBinIndexList[gscBinIndexCurrent];
          gscBinIndexCurrent++;
        }
      } else {
        for (pFileCommon->interesting_index = 0; pFileCommon->interesting_index < pFileCommon->interesting_count; pFileCommon->interesting_index) {
          pFlareCandidate = &pFileCommon->candidate_table[pFileCommon->interesting_index];
          if (pFlareCandidate->interesting_peakEvaluation != 0) {
            pFlareCandidate->interesting_peakEvaluation |= (1 << PEAKEVALUATION_PRIMARY_GSC_BIN);
            pWorkEntry->curGscBinIndex =  pFlareCandidate->gsc_bin_index;
            fprintf(consoleHandle,"search_none line %5d flareCandidatesOnly %d interesting_count %2d gsc_bin_index %d interesting_peakEvaluation\n",
                    __LINE__,
                    pFileCommon->flareCandidatesOnly,
                    pFileCommon->interesting_count,
                    pFlareCandidate->gsc_bin_index,
                    pFlareCandidate->interesting_peakEvaluation);
            break;

            
          }      
        }

      }
    }

    if (pWorkEntry->curGscBinIndex >= pGscBin->total_gsc_bins) {
      /* We are done */
      break;
    }

    GetBinCenter(pGscBin,pWorkEntry->curGscBinIndex,&lon,&lat,"search_none");
    if ((lat <  MIN_TC_DECLINATION) ||
        (lat > MAX_TC_DECLINATION) ||
        (lon < MIN_TC_RIGHTASCENSION) ||
        (lon > MAX_TC_RIGHTASCENSION)) {
      pWorkEntry->curGscBinIndex++;
      continue;
    }
    

    if (pFileCommon->minimumGalacticLatitude > -90.0) {
      wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
      if (lat < pFileCommon->minimumGalacticLatitude) {
        pWorkEntry->curGscBinIndex++;
        continue;
      }

    }



    if ((++keepaliveCheckModulus % KEEPALIVE_CHECK_MODULUS) == 0) {
      keepaliveCheckModulus = 0;
      
      
      gotAnswer = GetPhotometryGlobal(pFileCommon->pPhotConnection,pPhotGlobal);
      if ((gotAnswer != 1)  || (pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES)) {
      
        time(&curTime);
        
        ptr = localtime(&curTime);
        strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
        
        fprintf(consoleHandle,"search_none aborting at %s with gotAnswer %d keepalive %d\n",
                timestr,gotAnswer,pPhotGlobal->keepalive);
        fprintf(outHandle,"search_none aborting at %s with gotAnswer %d keepalive %d\n",
                timestr,gotAnswer,pPhotGlobal->keepalive);
        keepaliveFlag = 0;

        break;
      }
    }
    /* Find all of the images in and near this bin */

    LocateNoneImages(pGscBin,pFileCommon,pWorkEntry->curGscBinIndex,&pMagnitudeTable,&allocMagnitudes,&curMagnitudes,pFileCommon->catalogString,1,0);
    if (pFileCommon->enableTransientSearch > 0) {
      pFileCommon->pAllMagnitudeTable = pMagnitudeTable;
      pFileCommon->allMagnitudeCount = curMagnitudes;
    }
#ifdef LOS_DEBUG
    if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
      printf("Searching bin %d, got %d magnitudes<br / > \n",pWorkEntry->curGscBinIndex,curMagnitudes);
      printf("at bin %d\n",LOS_GSC_BIN);
    }
#endif


    curPrintMagnitudes += curMagnitudes;
    if ((++filePrintModulus % FILE_PRINT_MODULUS) == 0) {
      filePrintModulus = 0;
      time(&curTime);
      curTime -= pFileCommon->startTime;
      fprintf(consoleHandle,"Read %12lld magnitudes for curGscBinIndex %10d at %5d seconds\n",curPrintMagnitudes,pWorkEntry->curGscBinIndex,curTime);
      fprintf(outHandle,"Read %12lld magnitudes for curGscBinIndex %10d at %5d seconds\n",curPrintMagnitudes,pWorkEntry->curGscBinIndex,curTime);
      fflush(outHandle);
      curPrintMagnitudes = 0;
    }
#ifdef LOS_DEBUG
    if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
      int losIndex;
      for (losIndex = 0; losIndex < curMagnitudes; losIndex++) {
        pCurStarImage = &pMagnitudeTable[losIndex];
        pFileStarImage = pCurStarImage->pFileStarImage;
        fprintf(consoleHandle,"DEBUG: line %d index %d ra %f dec %f %s%05d_s%d gsc_bin_index %d REF %lld\n",__LINE__,losIndex,pFileStarImage->ra,pFileStarImage->dec,pCurStarImage->series,pFileStarImage->plateNumber,pFileStarImage->solutionNumber,pFileStarImage->gsc_bin_index,pFileStarImage->REFNumber);
        if (pFileStarImage->REFNumber == LOS_REF_NUMBER) {
          fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pFileStarImage->REFNumber);
        }

      }
    }
#endif /* LOS_DEBUG */
    if (keepaliveFlag == 0) {
      break;
    }
    if (curMagnitudes > 0) {
      pWorkEntry->totalMagnitudes += curMagnitudes;
      if (curMagnitudes >= pVectorStore->vectorCounterBlock.vectorAlloc) {
        pVectorStore->vectorCounterBlock.vectorAlloc += curMagnitudes+VECTOR_ALLOC_INCREMENT;
        ReallocVectorStore(consoleHandle,pVectorStore);  
      }
      /* Search for outliers */
      ProcessStars(consoleHandle,pFileCommon,pVectorStore,pMagnitudeTable,curMagnitudes,pWorkEntry->curGscBinIndex,astrometryHandle);

#ifdef LOS_DEBUG
      {
        pCurStarImage = pMagnitudeTable;
        pFileStarImage = pCurStarImage->pFileStarImage;
      
        if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
          fprintf(consoleHandle,"DEBUG: %d curMagnitudeIndex %d, gsc_bin_index %d\n",__LINE__,curMagnitudeIndex,pFileStarImage->gsc_bin_index);
        }
        if (pFileStarImage->REFNumber == LOS_REF_NUMBER) {
          fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pFileStarImage->REFNumber);
        }
      }
#endif /* LOS_DEBUG */


      magnitudeCount = 0;
      /* Because we are selecting images from multiple bins, collate these all into a single file table */
      for (index = 0; index < curMagnitudes; index++) {
        if (magnitudeCount >= noneMagnitudeAlloc) {
          noneMagnitudeAlloc += 1000;
          tmpNoneMagnitudeTable = realloc(pNoneMagnitudeTable,noneMagnitudeAlloc * sizeof(FILESTARIMAGE));
          if (tmpNoneMagnitudeTable == NULL) {
            fprintf(consoleHandle,"ERROR: failed to allocate pNoneMagnitudeTable of size %d\n",noneMagnitudeAlloc);
            exit(-1);
          }
          pNoneMagnitudeTable = tmpNoneMagnitudeTable;
          tmpNoneMagnitudeTable = NULL;
        }
        pCurStarImage = &pMagnitudeTable[index];
        pFileStarImage = pCurStarImage->pFileStarImage;
        if (pFileStarImage->gsc_bin_index == pWorkEntry->curGscBinIndex) {
          plateCountTable = pVectorStore->plateCountTable[pFileStarImage->seriesId];
          if (plateCountTable == NULL) {
            plateCountTable = (PPLATECOUNT)calloc(pFileCommon->maxPlateNumber[pFileStarImage->seriesId]+1,sizeof(PLATECOUNT));
            if (plateCountTable == NULL) {
              fprintf(consoleHandle,"ERROR: failed to allocate plateCountTable of size %d\n",pFileCommon->maxPlateNumber[pFileStarImage->seriesId]+1);
              exit(-1);

            }
            pVectorStore->plateCountTable[pFileStarImage->seriesId] = plateCountTable;
          }
          if ((pFileStarImage->plateNumber <= 0) ||
              (pFileStarImage->plateNumber > pFileCommon->maxPlateNumber[pFileStarImage->seriesId])) {
            fprintf(consoleHandle,"ERROR: illegal plate number %d\n",pFileStarImage->plateNumber);
            exit(-1);
          }
          plateCountTable[pFileStarImage->plateNumber].measurementCount++;
          if (pFileStarImage->REFNumber != 0) {
            plateCountTable[pFileStarImage->plateNumber].matchCount++;

          }
        }


#if 0 /* Change of October 25, 2014 - move all filtering to ProcessNoneImagesX */
#ifdef PERMISSIVE_ID_TABLE
        /* Show everything that is plotted, ignore second quality plates */
        if ((pFileStarImage->REFNumber == 0) &&
            ((pFileStarImage->AFLAGS & (~FILTER_AMASK_PLOT)) != 0)) {
          continue;
        }
#else /* PERMISSIVE_ID_TABLE */

        if (pFileStarImage->REFNumber == 0) {
#if 1
          /* New filtering algorithm of November 8, 2012 */
          if ((pFileStarImage->AFLAGS & FILTER_AMASK_NONE_CANDIDATE) != 0) {
            continue;
          }
#else 
          if (pFileStarImage->AFLAGS != 0) {
            continue;
          }
#endif
          qualityMask = GetPlateQualityMask(pFileCommon,pFileStarImage->seriesId,pFileStarImage->plateNumber,&quality,&pColorterm,&pColorflag,pFileStarImage->REFNumber);
#if 1
          /* New filtering algorithm of November 8, 2012 */
          if ((quality & (QUALITY_MULTIPLE|QUALITY_GRATING|QUALITY_SPECTRA)) != 0) {
            continue;
          }            
#else
          if ((qualityMask & (1 << pFileStarImage->spatial_bin)) != 0) {
            pFileCommon->qualityRejectCount++;
            continue;
          }
#endif
        }
#endif /* PERMISSIVE_ID_TABLE */
#endif /* Change of October 25, 2014 - move all filtering to ProcessNoneImagesX */
        if ((pFileCommon->catalogNumber == CATALOG_KEPLER) && 
            (RA2fpix(pFileStarImage->ra,pFileStarImage->dec,&chn_n,&p_row,&p_coln) == 0)) {
          pFileCommon->counterBlock.chipRejectCount++;
          continue;
        }

        pNoneImage  =&pNoneMagnitudeTable[magnitudeCount];
        memcpy(pNoneImage,pFileStarImage,sizeof(FILESTARIMAGE));
#ifdef LOS_DEBUG
        if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
          fprintf(consoleHandle,"DEBUG2: %d index %d ra %f dec %f %s%05d_s%d gsc_bin_index %d REF %lld\n",__LINE__,magnitudeCount,pFileStarImage->ra,pFileStarImage->dec,pCurStarImage->series,pFileStarImage->plateNumber,pFileStarImage->solutionNumber,pFileStarImage->gsc_bin_index,pFileStarImage->REFNumber);
        }
#endif /* LOS_DEBUG */
        magnitudeCount++;
      }
#ifdef LOS_DEBUG
      if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
        fprintf(consoleHandle,"line %d curMagnitudes %d magnitudeCount %d\n",__LINE__,curMagnitudes,magnitudeCount);
      }
#endif /* LOS_DEBUG */
      target_nrecs = 0;
#if 0
      if ((pWorkEntry->curGscBinIndex >= 6309888) && (pWorkEntry->curGscBinIndex <= 6350848)) {
        printf("\nline %d at curGscBinIndex %d entering\n",__LINE__,pWorkEntry->curGscBinIndex);
      }
#endif
      ProcessNoneImagesX(pGscBin,pFileCommon,&target_table,&target_nrecs,&target_alloc,pNoneMagnitudeTable,magnitudeCount,pFileCommon->verbose,pFileCommon->histogramTable,pFileCommon->minimumImageCount,&pWorkEntry->groupNumber,databaseHandle,0,0,-1,Sky2KTable,sky2k_size,maxGSCRadius,pWorkEntry->curGscBinIndex);
#if 0
      if ((pWorkEntry->curGscBinIndex >= 6309888) && (pWorkEntry->curGscBinIndex <= 6350848)) {
        printf("\nline %d at curGscBinIndex %d exiting\n",__LINE__,pWorkEntry->curGscBinIndex);
      }
#endif

#ifdef LOS_DEBUG
      for (target_index = 0; target_index < target_nrecs; target_index++) {
        pTarget = &target_table[target_index];
        if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
          fprintf(consoleHandle,"line %d target_index %d, REF %s gsc_bin_index %d\n",__LINE__,target_index,pTarget->REF,pTarget->gsc_bin_index);       

        }
      }
#endif /* LOS_DEBUG */ 
      for (target_index = 0; target_index < target_nrecs; target_index++) {
        group_image_count = 0;
        pTarget = &target_table[target_index];
#ifdef REMATCH_BUILD
#ifdef REMATCH_DEBUG
        if (pTarget->bestRematchNumber != 0) {
          printf("At line %d curGscBinIndex %d target_index %d\n",__LINE__,pWorkEntry->curGscBinIndex,target_index);
        }
#endif /* REMATCH_DEBUG */
        if (pTarget->bestRematchNumber != 0) {
          continue;
        }
#endif /* REMATCH_BUILD */
        memset(pCurStarEntry,0,sizeof(STARENTRY));
        strcpy(pCurStarEntry->REF,pTarget->REF);
        GetREFNumber(pTarget->REF,&pCurStarEntry->REFNumber,&refType,1,0);
        pCurStarEntry->ra = pTarget->ra;
        pCurStarEntry->dec = pTarget->dec;
        pCurStarEntry->gsc_bin_index = pTarget->gsc_bin_index;
        pCurStarEntry->DecPM = 0;
        pCurStarEntry->RaPM = 0;
        pCurStarEntry->Stdmag = pTarget->magcal_magdep; /* variable_search.m needs some magnitude to work with */
#ifdef LOS_DEBUG
        if (pCurStarEntry->REFNumber == LOS_REF_NUMBER) {
          fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pCurStarEntry->REFNumber);
        }
#endif /* LOS_DEBUG */
        pTarget->Stdmag = pTarget->magcal_magdep;
        if (pTarget->groupCount >= group_image_alloc) {
          if (pGroupImageTable != NULL) {
            free(pGroupImageTable);
          }
          group_image_alloc = pTarget->groupCount+100;
          pGroupImageTable = (PPHOTSTARIMAGE)calloc(group_image_alloc,sizeof(PHOTSTARIMAGE));
          if (pGroupImageTable == NULL) {
            fprintf(consoleHandle,"ERROR: failed to allocate pGroupImageTable of size %d\n",group_image_alloc);
            exit(-1);
          }
        }
        for (index = 0; index < magnitudeCount; index++) {
          pCurStarImage = &pMagnitudeTable[index];
          pFileStarImage = &pNoneMagnitudeTable[index];
          if (pCurStarEntry->REFNumber == pFileStarImage->REFNumber) {

#ifdef LOS_DEBUG
            if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
              fprintf(consoleHandle,"line %d index %d magnitudeCount %d ra %7.4f dec %7.4f magcal_magdep %6.2f for plate %s%05d_s%d\n",
                      __LINE__,
                      index,
                      magnitudeCount,
                      pFileStarImage->ra,
                      pFileStarImage->dec,
                      pFileStarImage->magcal_magdep,
                      pCurStarImage->series,
                      pFileStarImage->plateNumber,
                      pFileStarImage->solutionNumber);
            }
#endif /* LOS_DEBUG */

            pGroupImage = &pGroupImageTable[group_image_count];
            memcpy(pGroupImage,pCurStarImage,sizeof(PHOTSTARIMAGE));
            pGroupImage->pFileStarImage = pFileStarImage;
            GetREF(pFileStarImage->REFNumber,pGroupImage->REF,1,0);
            group_image_count++;
          }
        }
        if (group_image_count != pTarget->groupCount) {
          fprintf(consoleHandle,"ERROR: group_image_count %d does not agree with groupCount %d ra %f dec %f gsc_bin_index %d curGscBinIndex %d REF %s line %d\n",group_image_count,pTarget->groupCount,pTarget->ra,pTarget->dec,pTarget->gsc_bin_index,pWorkEntry->curGscBinIndex,pTarget->REF,__LINE__);
#if 0
          exit(-1);
#endif
        }

#ifdef LOS_DEBUG
        if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
          int losIndex;
          for (losIndex = 0; losIndex < curMagnitudes; losIndex++) {
            pCurStarImage = &pMagnitudeTable[losIndex];
            pFileStarImage = pCurStarImage->pFileStarImage;
            fprintf(consoleHandle,"DEBUG4: line %d index %d ra %f dec %f %s%05d_s%d gsc_bin_index %d REF %lld\n",__LINE__,losIndex,pFileStarImage->ra,pFileStarImage->dec,pCurStarImage->series,pFileStarImage->plateNumber,pFileStarImage->solutionNumber,pFileStarImage->gsc_bin_index,pFileStarImage->REFNumber);
            if (pFileStarImage->REFNumber == LOS_REF_NUMBER) {
              fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pFileStarImage->REFNumber);
            }
          }
        }
#endif /* LOS_DEBUG */

        pWorkEntry->goodUnmatchedStars += ProcessLightcurveParameters(consoleHandle,pFileCommon,pFileCommon->pPhotConnection,pGroupImageTable,group_image_count,pCurStarEntry,pParameterVectorStore,minGoodStars,pPhotGlobal->currentVersion,idHandle,&pWorkEntry->unmatchedStars,debugMode,&pWorkEntry->staleVersionIdCount,damonSeriesId,excludeSeriesId,excludeSeriesCount,pTarget,pWorkEntry->curGscBinIndex);


      }
#ifdef LOS_DEBUG
      if (pWorkEntry->curGscBinIndex == LOS_GSC_BIN) {
        int losIndex;
        for (losIndex = 0; losIndex < curMagnitudes; losIndex++) {
          pCurStarImage = &pMagnitudeTable[losIndex];
          pFileStarImage = pCurStarImage->pFileStarImage;
          fprintf(consoleHandle,"DEBUG3: %d index %d ra %f dec %f %s%05d_s%d gsc_bin_index %d REF %lld\n",__LINE__,losIndex,pFileStarImage->ra,pFileStarImage->dec,pCurStarImage->series,pFileStarImage->plateNumber,pFileStarImage->solutionNumber,pFileStarImage->gsc_bin_index,pFileStarImage->REFNumber);
        }
      }
#endif /* LOS_DEBUG */
      if (pFileCommon->enableTransientSearch > 0) {
        ProcessMatchedImages(consoleHandle,pFileCommon,pWorkEntry,pParameterVectorStore,pStarTable,outHandle,idHandle,catalogHandle,indexHandle,pWorkEntry->curGscBinIndex);
      }

      pWorkEntry->totalProcessedMagnitudes += magnitudeCount;
      curMagnitudeIndex += magnitudeCount;
    } /* Nonzero curMagnitudes */
      /* On to the next image */
    pWorkEntry->curGscBinIndex++;
    if (pWorkEntry->curGscBinIndex >= pWorkEntry->endGscBinIndex) {
      pWorkEntry->curGscBinIndex = pGscBin->total_gsc_bins+1;
    }
    if (pWorkEntry->curGscBinIndex >= pGscBin->total_gsc_bins) {
      fprintf(consoleHandle,"done curGscBinIndex %d exceeds maximum %d\n",pWorkEntry->curGscBinIndex,pGscBin->total_gsc_bins);
      break;
    }
  }


  time(&curTime);
  curTime -= pFileCommon->startTime;
  
  pWorkEntry->queryCount = GetQueryCount();
  fprintf(consoleHandle,"Rank %2d  CPU %s Time: %d seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d gscBins, %d total magnitudes\n",pFileCommon->myrank,pFileCommon->processorName,curTime,pWorkEntry->queryCount,pWorkEntry->totalStars,pWorkEntry->recoveredStars,pWorkEntry->badModulusBin,pWorkEntry->binCount,pWorkEntry->totalMagnitudes);

  fprintf(outHandle,"Rank %2d  CPU %s Time: %d seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d gscBins, %d total magnitudes\n",pFileCommon->myrank,pFileCommon->processorName,curTime,pWorkEntry->queryCount,pWorkEntry->totalStars,pWorkEntry->recoveredStars,pWorkEntry->badModulusBin,pWorkEntry->binCount,pWorkEntry->totalMagnitudes);



  fprintf(consoleHandle," queries/sec: %f totalStars/sec: %f recoveredStars+badModulusBin.sec: %f gscBins/sec %f\n",
          (1.0 * pWorkEntry->queryCount)/(1.0*curTime),
          (1.0 * pWorkEntry->totalStars)    /(1.0*curTime),
          (1.0 * (pWorkEntry->recoveredStars + pWorkEntry->badModulusBin) )    /(1.0*curTime),
          (1.0 * pWorkEntry->binCount  )    /(1.0*curTime));
  fprintf(outHandle," queries/sec: %f totalStars/sec: %f recoveredStars+badModulusBin.sec: %f gscBins/sec %f\n",
          (1.0 * pWorkEntry->queryCount)/(1.0*curTime),
          (1.0 * pWorkEntry->totalStars)    /(1.0*curTime),
          (1.0 * (pWorkEntry->recoveredStars+pWorkEntry->badModulusBin) )    /(1.0*curTime),
          (1.0 * pWorkEntry->binCount  )    /(1.0*curTime));
 


#if 0
  fprintf(consoleHandle,"... %d noGscEntry\n",pWorkEntry->noGscEntryCount);
  fprintf(consoleHandle,"... %d staleVersionId\n",pWorkEntry->staleVersionIdCount);
  fprintf(outHandle,"... %d noGscEntry\n",pWorkEntry->noGscEntryCount);
  fprintf(outHandle,"... %d staleVersionId\n",pWorkEntry->staleVersionIdCount);


  fprintf(consoleHandle,"New curGscBinIndex %d\n",pWorkEntry->curGscBinIndex);
  fprintf(outHandle,"New curGscBinIndex %d\n",pWorkEntry->curGscBinIndex);

  fprintf(consoleHandle,"Total Magnitudes %d, Total NONE magnitudes %d Total NOMATCH magnitudes %d Total GSC nonblend magnitudes %d \n",
          pWorkEntry->totalProcessedMagnitudes,pWorkEntry->totalNoneMagnitudes,pWorkEntry->totalBlendNomatchMagnitudes,pWorkEntry->totalProcessedMagnitudes - pWorkEntry->totalNoneMagnitudes - pWorkEntry->totalBlendNomatchMagnitudes);
  fprintf(outHandle,"Total Magnitudes %d, Total NONE magnitudes %d Total NOMATCH magnitudes %d Total GSC nonblend magnitudes %d \n",
          pWorkEntry->totalProcessedMagnitudes,pWorkEntry->totalNoneMagnitudes,pWorkEntry->totalBlendNomatchMagnitudes,pWorkEntry->totalProcessedMagnitudes - pWorkEntry->totalNoneMagnitudes - pWorkEntry->totalBlendNomatchMagnitudes);
  fprintf(consoleHandle,"Total groups %d\n",pWorkEntry->groupNumber-1);
  fprintf(outHandle,"Total groups %d\n",pWorkEntry->groupNumber-1);
#endif
  assert(sizeof(pFileCommon->histogramTable) == sizeof(pWorkEntry->histogramTable));
  memcpy(pWorkEntry->histogramTable,pFileCommon->histogramTable,sizeof(pFileCommon->histogramTable));


  if (pMagnitudeTable != NULL) {
    free(pMagnitudeTable);
  }
  if (pNoneMagnitudeTable != NULL) {
    free (pNoneMagnitudeTable);
  }
  if (pFileCommon->enableTransientSearch > 0) {
    free(pStarTable);
    Close(catalogHandle);
    Close(indexHandle);


  }


  for (seriesIndex = 0; seriesIndex <= MAX_SERIES; seriesIndex++) {
    plateCountTable = pVectorStore->plateCountTable[seriesIndex];

    if (plateCountTable == NULL) {
      continue;
    }
    plateCountTable = pVectorStore->plateCountTable[seriesIndex];
    for (plateIndex = 0; plateIndex <= pFileCommon->maxPlateNumber[seriesIndex]; plateIndex++) {
      if ((plateCountTable[plateIndex].measurementCount != 0)  ||
          (plateCountTable[plateIndex].matchCount != 0)  ||
          (plateCountTable[plateIndex].lightcurveCount != 0)  ||
          (plateCountTable[plateIndex].outlierCount != 0)) {
        double outlierRatio = 0.0;
        double matchRatio = 0.0;
        if (plateCountTable[plateIndex].lightcurveCount > 0) {
          outlierRatio = (1.0*plateCountTable[plateIndex].outlierCount)/(1.0*plateCountTable[plateIndex].lightcurveCount);
        }
        if (plateCountTable[plateIndex].measurementCount > 0) {
          matchRatio = (1.0*plateCountTable[plateIndex].matchCount)/(1.0*plateCountTable[plateIndex].measurementCount);
        }
        fprintf(plateHandle,"%s%05d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%f\n",
                GetSeriesString(seriesIndex,1),
                plateIndex,
                seriesIndex,
                plateIndex,
                plateCountTable[plateIndex].measurementCount,
                plateCountTable[plateIndex].matchCount,
                matchRatio,
                plateCountTable[plateIndex].lightcurveCount,
                plateCountTable[plateIndex].outlierCount,
                outlierRatio);
                
      }
    }

  }


  fprintf(consoleHandle,"Max vectorAlloc %d, lightcurveCount %d lightcurvePointCount %d, outlierCount %d, qualityRejectCount %d\n",
          pVectorStore->vectorCounterBlock.vectorAlloc,
          pVectorStore->vectorCounterBlock.lightcurveCount,
          pVectorStore->vectorCounterBlock.lightcurvePointCount,
          pVectorStore->vectorCounterBlock.outlierCount,
          pVectorStore->vectorCounterBlock.qualityRejectCount);
  fprintf(consoleHandle,"points read %d %d, pointsRejected %d, insufficientPoints %d, highRMSPoints %d\n",
          pVectorStore->vectorCounterBlock.pointsRead,
          pVectorStore->vectorCounterBlock.pointsRead2,
          pVectorStore->vectorCounterBlock.pointsRejected,
          pVectorStore->vectorCounterBlock.insufficientPoints,
          pVectorStore->vectorCounterBlock.highRMSPoints);
  fprintf(consoleHandle,"Unmatch stars %d good unmatched stars %d\n",pWorkEntry->unmatchedStars,pWorkEntry->goodUnmatchedStars);
  fprintf(consoleHandle,"Total Stars %10lld; Total Stars with Outliers %10lld stars plotted %10lld none candidates %10lld\n",
          pVectorStore->vectorCounterBlock.totalStars,
          pVectorStore->vectorCounterBlock.totalOutliers,
          pVectorStore->vectorCounterBlock.totalPlot2,
          pVectorStore->vectorCounterBlock.totalNoneCandidates);
#if 0
  fprintf(consoleHandle,"line %d source vectorAlloc %d -> target vectorAlloc %d\n",__LINE__,pVectorStore->vectorCounterBlock.vectorAlloc,pWorkEntry->vectorCounterBlock.vectorAlloc);
#endif
  memcpy(&pWorkEntry->vectorCounterBlock,&pVectorStore->vectorCounterBlock,sizeof(VECTORCOUNTERBLOCK));
#if 0
  fprintf(consoleHandle,"line %d source vectorAlloc %d -> target vectorAlloc %d\n",__LINE__,pVectorStore->vectorCounterBlock.vectorAlloc,pWorkEntry->vectorCounterBlock.vectorAlloc);
#endif


#ifdef OUTPUT_COUNTS
  OutputCounts("AFLAG  ",AflagsTable,NULL,AflagsTableSize,pVectorStore->AFLAGS_COUNTS);
  OutputCounts("BFLAG  ",BflagsTable,NULL,BflagsTableSize,pVectorStore->BFLAGS_COUNTS);
  OutputCounts("QUALITY",NULL,qualityMasks,qualityTableSize,pVectorStore->QUALITY_COUNTS);
#endif /* OUTPUT_COUNTS */
#if 0
  OutputRootCounts(consoleHandle,pFileCommon);
  for (dradIndex = 0; dradIndex <= MAX_GROUP_DRAD_HISTOGRAM; dradIndex++) {
    if (pFileCommon->counterBlock.maxGroupHistogram[dradIndex] != 0) {
      fprintf(consoleHandle,"Group Drad: %3d arcsec count: %8d\n",dradIndex,pFileCommon->counterBlock.maxGroupHistogram[dradIndex]);
    }
  }
#endif
  if (astrometryHandle != NULL) {
    int magindex;
    int dradindex;
    int dradvalue;
    for (magindex = 0; magindex < MAX_DRAD_HISTOGRAM_MAG; magindex++) {
      for (dradindex = 0; dradindex < (2*(MAX_DRAD_RANGE)+1); dradindex++) {
        dradvalue = dradindex-MAX_DRAD_RANGE;
        if ((pFileCommon->draHistogram[magindex][dradindex] != 0) ||
            (pFileCommon->ddecHistogram[magindex][dradindex] != 0)) {
          fprintf(astrometryHandle,"%d\t%d\t%lld\t%lld\n",magindex,dradvalue,pFileCommon->draHistogram[magindex][dradindex],pFileCommon->ddecHistogram[magindex][dradindex]);
        }

      }

    }


    fclose(astrometryHandle);
    astrometryHandle = NULL;
    statResult = stat(pWorkEntry->astrometryLogname,&filestats);
    if (statResult != 0) {
      fprintf(consoleHandle,"ERROR: failed to obtain stats on %s\n",pWorkEntry->astrometryLogname);
    } else {
      pWorkEntry->astrometryFilesize = filestats.st_size;
    }


  }



  FreeVectorStore(pVectorStore);
  FreeParameterVectorStore(pParameterVectorStore);
  if (target_table != NULL) {
    free(target_table);
  }
  if (pGroupImageTable == NULL) {
    free(pGroupImageTable);
  }

  fclose(outHandle);
  statResult = stat(pWorkEntry->outLogname,&filestats);
  if (statResult != 0) {
    fprintf(consoleHandle,"ERROR: failed to obtain stats on %s\n",pWorkEntry->outLogname);
  } else {
    pWorkEntry->outFilesize = filestats.st_size;
  }


  fclose(databaseHandle);
  statResult = stat(pWorkEntry->databaseLogname,&filestats);
  if (statResult != 0) {
    fprintf(consoleHandle,"ERROR: failed to obtain stats on %s\n",pWorkEntry->databaseLogname);
  } else {
    pWorkEntry->databaseFilesize = filestats.st_size;
  }
  fclose(plateHandle);

  statResult = stat(pWorkEntry->plateLogname,&filestats);
  if (statResult != 0) {
    fprintf(consoleHandle,"ERROR: failed to obtain stats on %s\n",pWorkEntry->plateLogname);
  } else {
    pWorkEntry->plateFilesize = filestats.st_size;
  }


  if (idHandle != NULL) {
    fclose(idHandle);
  }
  statResult = stat(pWorkEntry->idLogname,&filestats);
  if (statResult != 0) {
    fprintf(consoleHandle,"ERROR: failed to obtain stats on %s\n",pWorkEntry->idLogname);
  } else {
    pWorkEntry->idFilesize = filestats.st_size;
  }
  pWorkEntry->propermotionFilesize = CloseProperMotions(consoleHandle,pFileCommon,pWorkEntry->propermotionLogname);



  fprintf(consoleHandle,"CHILD: SearchProcess exiting with rank %d queue index %d console handle 0x%x\n",pFileCommon->myrank,pWorkEntry->queueIndex,consoleHandle);
  fflush(consoleHandle);


  return(0);
}

void ProcessMatchedImages(FILE *consoleHandle,PFILECOMMON pFileCommon,PWORKQUEUE pWorkEntry,PPARAMETERVECTORSTORE pParameterVectorStore,PSTARENTRY pStarTable,FILE *outHandle,FILE *idHandle,File catalogHandle,File indexHandle,int curGscBinIndex) {

  int binCount = 0;
  int cur_gsc_bin_index = -1;
  int curStars = 0;
  int debugMode = 0;
  int errorFlag = 0;
  int excludeSeriesCount = 0;
  int goodStars = 0;
  int gscIndex = 0;
  int magnitudeAlloc = 0;
  int magnitudeIndex = 0;
  int magnitudesModulus =  MAGNITUDES_MODULUS;
  int minGoodStars = MIN_SUMMARY_GOODSTARS;
  int noGscEntryCount = 0;
  int oldBinCount = 0;
  int oldGoodStars = 0;
  int oldQueryCount = 0;
  int oldTotalStars = 0;
  int processedCount = 0;
  int readItems = 0;
  int staleVersionIdCount = 0;
  int starIndex = 0 ;
  int starIndex2 = 0;
  int summaryCount = 0;
  int totalStars = 0;
  int usecIndex = 0;
  long long curUsec = 0;
  long long deltaTable[MAX_CHARGE_ENTRY];
  long long oldTotalMagnitudes = 0;
  long long totalMagnitudes = 0;
  long long totalTable[MAX_CHARGE_ENTRY];
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  PFILESTARIMAGE pFileStarImage = NULL;
  PGSCIMAGE pCurGscImage = NULL;
  PPHOTSTARIMAGE pCurSumstarimage;
  PHOTTARGET target;
  PPHOTTARGET pTarget = &target;
  PSTARENTRY pCurStarEntry2 = NULL;
  PSTARENTRY pCurStarEntry = NULL;
  PPHOTSTARIMAGE pMagnitudeTable = NULL;  
  STARINDEX curStarIndex;
  PSTARINDEX pCurStarIndex = &curStarIndex;
  PMATCHEDCOMMON pMatchedCommon = &pWorkEntry->matchedCommon;
  time_t curTime = 0;
  time_t deltaTime = 0;
  int *seriesVersionArray;
  
#ifdef LOS_DEBUG
  int debugPrint = 0;
  if (curGscBinIndex == LOS_GSC_BIN) {
    printf("At line %d curGscBinIndex %d\n",__LINE__,curGscBinIndex);
    debugPrint = 1;
  }

#endif /* LOS_DEBUG */

  /* Find out how many unique REF numbers there are in the table */
  curStars = 0;
  pMagnitudeTable = pFileCommon->pAllMagnitudeTable;
  for (magnitudeIndex = 0; magnitudeIndex < pFileCommon->allMagnitudeCount;magnitudeIndex++) {
    pCurSumstarimage = &pFileCommon->pAllMagnitudeTable[magnitudeIndex];
    pFileStarImage = pCurSumstarimage->pFileStarImage;
#ifdef LOS_DEBUG
    if (pFileStarImage->REFNumber == LOS_REF_NUMBER) {
      fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pFileStarImage->REFNumber);
    }
#endif /* LOS_DEBUG */
#ifdef LOS_DEBUG
    if (debugPrint != 0) {
      printf("At line %d gsc_bin_index %d REFNumber %lld\n",__LINE__,pFileStarImage->gsc_bin_index,pFileStarImage->REFNumber);
    }
#endif /* LOS_DEBUG */
    if (pFileStarImage->REFNumber == 0) {
      continue;
    }

    if (pFileStarImage->gsc_bin_index != curGscBinIndex) {
      continue;
    }
#ifdef LOS_DEBUG
    if ((pFileStarImage->plateNumber == LOS_PLATENUMBER) &&
        (pFileStarImage->seriesId == LOS_SERIESID)) {
      printf("Line %d GetLatestVersionID() At seriesId %d plateNumber %d\n",__LINE__,pFileStarImage->plateNumber,pFileStarImage->seriesId);
    }
#endif /* LOS_DEBUG */

    seriesVersionArray = pFileCommon->seriesVersionTable[pFileStarImage->seriesId];
    if (seriesVersionArray == NULL) {
      seriesVersionArray = (int *)calloc(pFileCommon->maxPlateNumber[pFileStarImage->seriesId]+1,sizeof(int));
      if (seriesVersionArray == NULL) {
        printf("ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[pFileStarImage->seriesId]+1);
        exit(-1);
      }
      pFileCommon->seriesVersionTable[pFileStarImage->seriesId] = seriesVersionArray;
    }
    if (seriesVersionArray[pFileStarImage->plateNumber] < pFileStarImage->versionId) {
      seriesVersionArray[pFileStarImage->plateNumber] = pFileStarImage->versionId;
    }




    for (starIndex = 0; starIndex < curStars; starIndex++) {
      pCurStarEntry = &pStarTable[starIndex];
#ifdef LOS_DEBUG
      if (pCurStarEntry->REFNumber == LOS_REF_NUMBER) {
        fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pCurStarEntry->REFNumber);
      }
#endif /* LOS_DEBUG */
      if (pFileStarImage->REFNumber == pCurStarEntry->REFNumber) {
        break;
      }
    }
#ifdef LOS_DEBUG
    if (debugPrint != 0) {
      printf("At line %d gsc_bin_index %d REFNumber %lld\n",__LINE__,pFileStarImage->gsc_bin_index,pFileStarImage->REFNumber);
    }
#endif /* LOS_DEBUG */
    if (starIndex >= curStars) {
      pCurStarEntry = &pStarTable[curStars];
      memset(pCurStarEntry,0,sizeof(STARENTRY));
      pCurStarEntry->REFNumber = pFileStarImage->REFNumber;
#ifdef LOS_DEBUG
      if (pCurStarEntry->REFNumber == LOS_REF_NUMBER) {
        fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pCurStarEntry->REFNumber);
      }
#endif /* LOS_DEBUG */
      pCurStarEntry->gsc_bin_index = pFileStarImage->gsc_bin_index;
      pCurStarEntry->versionId = pFileStarImage->versionId;
      if ((pFileStarImage->ra_2 < 900) && pFileStarImage->dec_2 <= 90) {
        pCurStarEntry->ra = pFileStarImage->ra_2;
        pCurStarEntry->dec = pFileStarImage->dec_2;
        pCurStarEntry->RaPM = pFileStarImage->RaPM;
        pCurStarEntry->DecPM = pFileStarImage->DecPM;
      }
      pCurStarEntry->updateflag = 1;
      GetREF(pCurStarEntry->REFNumber,pCurStarEntry->REF,0,0);
#if 0
      pCurStarEntry->Stdmag = pFileStarImage->Stdmag;
      pCurStarEntry->color = pFileStarImage->color;
      pCurStarEntry->MAGFlag = pFileStarImage->MAGFlag;
      pCurStarEntry->class = pFileStarImage->class;
      pCurStarEntry->VFlag = pFileStarImage->VFlag;
#endif
       
      curStars++;
      if (curStars > (STAR_LIMIT-1)) {
        fprintf(consoleHandle,"ERROR: curStars %d exceeds STAR_LIMIT %d curMagnitudes %d\n",curStars,STAR_LIMIT,pFileCommon->allMagnitudeCount);
        exit(-1);
      }
    }
  }

  totalStars += curStars;
#if 0
  if (pFileCommon->verbose) {
    time(&curTime);
    curTime -= pFileCommon->startTime;
    fprintf(consoleHandle,"Got %d stars total %d at %d seconds %s\n",curStars,totalStars,curTime,pStarTable->REF);

  }
#endif
  processedCount = 0;
  for (starIndex = 0; starIndex < curStars; starIndex++) {
    if (errorFlag) {
      break;
    }
    pCurStarEntry = &pStarTable[starIndex];      
#ifdef LOS_DEBUG
    if (pCurStarEntry->REFNumber == LOS_REF_NUMBER) {
      fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pCurStarEntry->REFNumber);
    }
#endif /* LOS_DEBUG */
    if (pCurStarEntry->processedflag != 0) {
      /* We processed this one already */
      continue;
    }

    if (pFileCommon->allMagnitudeCount < 0) {
      errorFlag = 1;
      break;
    }

    binCount++;
    for (starIndex2 = starIndex; starIndex2 < curStars; starIndex2++) {
      if (errorFlag) {
        break;
      }
      pCurStarEntry2 = &pStarTable[starIndex2];
#ifdef LOS_DEBUG
      if (pCurStarEntry2->REFNumber == LOS_REF_NUMBER) {
        fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pCurStarEntry2->REFNumber);
      }
#endif /* LOS_DEBUG */
      if (pCurStarEntry2->processedflag != 0) {
        /* We processed this one already */
        continue;
      }
        
      pCurStarEntry2->processedflag = 1;
      processedCount++;




      if (cur_gsc_bin_index != pCurStarEntry2->gsc_bin_index) {
        cur_gsc_bin_index = pCurStarEntry2->gsc_bin_index;
        if (cur_gsc_bin_index >= pGscBin->total_gsc_bins) {
          fprintf(consoleHandle,"ERROR: gsc_bin_index %d exceeds total %d in update_summary2 for REFNUMBER %lld\n",
                  cur_gsc_bin_index,
                  pGscBin->total_gsc_bins,
                  pCurStarEntry2->REFNumber);
          cur_gsc_bin_index = -1;
          memset(pCurStarIndex,0,sizeof(STARINDEX));
        } else {  

          Seek(indexHandle,cur_gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
          readItems = Read(indexHandle,pCurStarIndex,sizeof(STARINDEX),1);
          if (readItems != 1) {
            fprintf(consoleHandle,"ERROR reading star index file errno: %d %s REF %s ra %f dec %f gsc_bin_index %d\n",
                    errno,
                    strerror(errno),
                    pCurStarEntry2->REF,
                    pCurStarEntry2->ra,
                    pCurStarEntry2->dec,
                    pCurStarEntry2->gsc_bin_index);
            cur_gsc_bin_index = -1;
            memset(pCurStarIndex,0,sizeof(STARINDEX));
          } else {
          
            if (pCurStarIndex->numStars > pFileCommon->gscImageAlloc) {
              pFileCommon->gscImageAlloc = pCurStarIndex->numStars+1000;
              pCurGscImage = realloc(pFileCommon->pGscImageTable,pFileCommon->gscImageAlloc*sizeof(GSCIMAGE));
              if (pCurGscImage == NULL) {
                fprintf(consoleHandle,"ERROR: failed to realloc pGscImageTable of size %d\n",pFileCommon->gscImageAlloc*sizeof(GSCIMAGE));
                exit(-1);
              }
              pFileCommon->pGscImageTable = pCurGscImage;
              pCurGscImage = NULL;
            }
            Seek(catalogHandle,pCurStarIndex->offset,SEEK_SET);
            readItems = Read(catalogHandle,pFileCommon->pGscImageTable,sizeof(GSCIMAGE),pCurStarIndex->numStars);
            if (readItems != pCurStarIndex->numStars) {
              fprintf(consoleHandle,"ERROR reading gsc catalog file\n");
              exit(-1);
            }
          }
        }
      }    


      if (pFileCommon->allMagnitudeCount >= minGoodStars) {

        /* Because the ra and dec has been modified for proper motion, get the original ra and dec from the GSC2.3.2 catalog */
        for (gscIndex = 0; gscIndex < pCurStarIndex->numStars; gscIndex++) {
          pCurGscImage = &pFileCommon->pGscImageTable[gscIndex];
          if (pCurGscImage->REFNumber == pCurStarEntry2->REFNumber) {
#ifdef LOS_DEBUG
            if (pCurStarEntry2->REFNumber == LOS_REF_NUMBER) {
              fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pCurStarEntry2->REFNumber);
            }
#endif /* LOS_DEBUG */
            pCurStarEntry2->ra = pCurGscImage->ra;
            pCurStarEntry2->dec = pCurGscImage->dec;
            pCurStarEntry2->MAGFlag = pCurGscImage->MAGFlag;
            pCurStarEntry2->class = pCurGscImage->class;
            pCurStarEntry2->VFlag = pCurGscImage->VFlag;
            pCurStarEntry2->RaPM = pCurGscImage ->RaPM;
            pCurStarEntry2->DecPM = pCurGscImage->DecPM;

            pCurStarEntry2->Stdmag = pCurGscImage->Stdmag;
            pCurStarEntry2->color = pCurGscImage->color;
            pCurStarEntry2->MAGFlag = pCurGscImage->MAGFlag;
            pCurStarEntry2->class = pCurGscImage->class;
            pCurStarEntry2->VFlag = pCurGscImage->VFlag;

            pCurStarEntry->ra = pCurGscImage->ra;
            pCurStarEntry->dec = pCurGscImage->dec;
            pCurStarEntry->MAGFlag = pCurGscImage->MAGFlag;
            pCurStarEntry->class = pCurGscImage->class;
            pCurStarEntry->VFlag = pCurGscImage->VFlag;
            pCurStarEntry->RaPM = pCurGscImage ->RaPM;
            pCurStarEntry->DecPM = pCurGscImage->DecPM;

            pCurStarEntry->Stdmag = pCurGscImage->Stdmag;
            pCurStarEntry->color = pCurGscImage->color;
            pCurStarEntry->MAGFlag = pCurGscImage->MAGFlag;
            pCurStarEntry->class = pCurGscImage->class;
            pCurStarEntry->VFlag = pCurGscImage->VFlag;

            break;
          }

        }
        if (gscIndex == pCurStarIndex->numStars) {
#ifdef LOS_DEBUG
          fprintf(consoleHandle,"WARNING: no gsc entry found in bin %d for star %s ra %f dec %f\n",
                  cur_gsc_bin_index,
                  pCurStarEntry2->REF,
                  pCurStarEntry2->ra,
                  pCurStarEntry2->dec);
#endif /* LOS_DEBUG */
          noGscEntryCount++;
          pCurStarEntry2->Stdmag = 99.0; /* Added Dec 24, 2016 to avoid problems in ProcessTransientCandidate */
          pCurStarEntry2->color = 99.0;
#if 0
          pCurStarEntry2->MAGFlag = 0;
          pCurStarEntry2->class = 0;
          pCurStarEntry2->VFlag = 0;
          pCurStarEntry2->RaPM = 0;
          pCurStarEntry2->DecPM = 0;
#endif
        }

        memset(pTarget,0,sizeof(PHOTTARGET));
        pTarget->Stdmag = pCurStarEntry2->Stdmag;
        pTarget->color = pCurStarEntry2->color;
        pTarget->ra = pCurStarEntry2->ra;
        pTarget->dec = pCurStarEntry2->dec;
        pTarget->MAGFlag = pCurStarEntry2->MAGFlag;
        pTarget->class = pCurStarEntry2->class;
        pTarget->VFlag = pCurStarEntry2->VFlag;
        pTarget->RaPM = pCurStarEntry2->RaPM;
        pTarget->DecPM = pCurStarEntry2->DecPM;
        pTarget->nearestCatalogStar.nearestREFflagGoodMag = -1;
        pTarget->nearestCatalogStar.nearestREFflag = -1;
            
#ifdef LOS_DEBUG
        if (pCurStarEntry2->REFNumber == LOS_REF_NUMBER) {
          fprintf(consoleHandle,"line %d at REFNumber %lld\n",__LINE__,pCurStarEntry2->REFNumber);
        }
#endif /* LOS_DEBUG */

        goodStars += ProcessLightcurveParameters(consoleHandle,pFileCommon,pPhotConnection,pMagnitudeTable,pFileCommon->allMagnitudeCount,pCurStarEntry2,pParameterVectorStore,minGoodStars,pMatchedCommon->versionId,idHandle,&summaryCount,debugMode,&staleVersionIdCount,damonSeriesId,excludeSeriesId,excludeSeriesCount,pTarget,curGscBinIndex);
            
        if (summaryCount >= SUMMARY_LIMIT) {
          summaryCount = 0;

          time(&curTime);
          deltaTime = curTime - pMatchedCommon->oldTime;
          pMatchedCommon->oldTime = curTime;
          curTime -= pFileCommon->startTime;
          fprintf(consoleHandle,"Processing      %d seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",curTime,GetQueryCount(),totalStars,goodStars,magnitudeAlloc,pParameterVectorStore->vectorAlloc,binCount,totalMagnitudes);
          fprintf(outHandle,"Processing      %d seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",curTime,GetQueryCount(),totalStars,goodStars,magnitudeAlloc,pParameterVectorStore->vectorAlloc,binCount,totalMagnitudes);
          if (deltaTime > 0) {
            fprintf(consoleHandle," queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec: %f magnitudes/sec %f\n",
                    (1.0 *(GetQueryCount()-oldQueryCount))/(1.0*deltaTime),
                    (1.0 *(totalStars-oldTotalStars))/(1.0*deltaTime),
                    (1.0 *(goodStars -oldGoodStars ))/(1.0*deltaTime),
                    (1.0 *(binCount - oldBinCount  ))/(1.0*deltaTime),
                    (1.0 *(totalMagnitudes - oldTotalMagnitudes))/(1.0*deltaTime));
            fprintf(outHandle," queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec: %f magnitudes/sec %f\n",
                    (1.0 *(GetQueryCount()-oldQueryCount))/(1.0*deltaTime),
                    (1.0 *(totalStars-oldTotalStars))/(1.0*deltaTime),
                    (1.0 *(goodStars -oldGoodStars ))/(1.0*deltaTime),
                    (1.0 *(binCount - oldBinCount  ))/(1.0*deltaTime),
                    (1.0 *(totalMagnitudes - oldTotalMagnitudes))/(1.0*deltaTime));
            oldQueryCount = GetQueryCount();
            oldTotalStars  = totalStars;
            oldGoodStars   = goodStars;
            oldBinCount    = binCount;
            oldTotalMagnitudes = totalMagnitudes;
          }
          ChargeTime(0,&curUsec,deltaTable,totalTable);
          fprintf(consoleHandle,"Intermediate Time (sec) ");
          fprintf(outHandle,"Intermediate Time (sec) ");
          for (usecIndex = 0; usecIndex < MAX_CHARGE_ENTRY; usecIndex++) {
            double outSeconds;
            outSeconds = (1.0 * deltaTable[usecIndex])/1000000.0;
            fprintf(consoleHandle,"%10.3f ",outSeconds);
            fprintf(outHandle,"%10.3f ",outSeconds);
          }
          fprintf(consoleHandle,"\n");
          fprintf(outHandle,"\n");

          ChargeTime(RESET_DELTA_ENTRY,&curUsec,deltaTable,totalTable);


          fflush(outHandle);

          


        }
      }
              

      
    }
  }

} /* End of ProcessMatchedImages */

