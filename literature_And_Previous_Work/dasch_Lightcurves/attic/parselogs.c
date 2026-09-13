// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* parselogs.c
 *
 *  This program extracts run statistics from a list of log files
 *
 *  Usage:  parselogs -l  <list of logs>  test with run_932.log for i24330
 *          Conversion to openmpi on Aug 8, 2012.  run_673.log
 *          First log containg Recstat is run_678.log
 *          Corrupt (edited) logs run_722a.log
 *          Corrupt (killed) logs run_718.log
 *          Corrupt (unknown) logs run_729.log
 *
 *  TODO: Ambiguity in \update_sextractor if SExtractor has not already been run
 *
 *           parselogs  -v  -l /home/scanner/Pipeline/los.list -o /home/scanner/Pipeline/los.db
 *
 *
           parselogs -l /home/scanner/backup/2014_01_14/production.list -o  /home/scanner/Pipeline/los.db > los9.tmp
           parselogs -r -l /home/scanner/backup/2014_01_14/production.list -o  /home/scanner/Pipeline/los.db > los8.tmp
           parselogs -l /home/scanner/backup/2014_01_14/reprocess.list -o  /home/scanner/Pipeline/los.db > los7.tmp
           parselogs -r -l /home/scanner/backup/2014_01_14/reprocess.list -o  /home/scanner/Pipeline/los.db > los6.tmp

 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm parselogs.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib/mysql  -l mysqlclient -o parselogs
 * 
 *
 * Dec 26, 2013 Edward J. Los - Initial version  
 * Jan 21, 2013 Edward J. Los - Find the mosaic processing time
 * Jan  6, 2020 Edward J. Los - Add '-a' for detailed examination of astrometry.net performance
 *                              the '-o' file will have timeSpent in seconds; plateCount; countCount (normalized cumulative count) , and sequence
 *                              sequence = 1 all Success; 2 all failed; 3 series Success; 4 series failed
 * Mar 30, 2020 Edward J. Los - Add '-p' to parse search_close3 output which looks for shifts in the position of small plates during scanning
 * Apr  6, 2020 Edward J. Los - Add "sectorCount" parsing.
 * Apr  9, 2020 Edward J. Los - Add "clip_count2"
 * Apr 10, 2020 Edward J. Los - Add maxMaskArea, conditionalize on maskArea and maskCenterFlag for individual masks
 * Apr 15, 2020 Edward J. Los - Add <search_close3.list> to -p to record plates actually examined.
 * Apr 18, 2020 Edward J. Los - Expand the SEQUENCE messages to aid in inspection plate selection
 * May  8, 2020 Edward J. Los - Add COMBINE_CENTER_AREA to combine maskCenterFlag in totalMasks sum, and maxCenterMaskArea maxMaskArea
 */   


#include <math.h>
#include <string.h>
#include "table.h"
#include "time.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER     1024
#define MAX_TIME_LENGTH 28
#define MAX_SOLUTION_NUMBER 10
#define COMBINE_CENTER_AREA 1

/* #define DEBUG_SINGLE "i24330" */
#define STEP_UNKNOWN            0  /* Must be zero or conditional leading to "Here we must charge this time to the unknown step" won't work */
#define STEP_MOSAICS            1
#define STEP_ASTROMETRYNET      2  /* find_astrometry.csh */
#define STEP_IMWCS              3
#define STEP_XFER1              4
#define STEP_PREPAREMOSAIC      5
#define STEP_FULLSEXTRACTOR     6  /* run_sextractor.csh - full execution */
#define STEP_SEXTRACTOR1        7  /* run_sextractor.csh - update_sectractor only */
#define STEP_SEARCHCLOSE        8
#define STEP_SCAMPFIRST         9  /* run_scamp.csh */
#define STEP_SEXTRACTOR2       10  /* run_sextractor_second.csh */
#define STEP_SCAMPASTROMETRY   11  /* Note: date brackets solve-field only (find_astrometry01.csh) */
#define STEP_SCAMPSECOND       12  /* run_scamp_01.csh */
#define STEP_SCAMPTHIRD        13  /* run_scamp.csh 2 */
#define STEP_SEXTRACTOR2B      14  /* run_sextractor_second.csh 2 */
#if 0
#define STEP_HEADER            xx  /* No 'seconds' output and no date information -- ignore this one */
#endif
#define STEP_WEDGE             15
#define STEP_MATCH2            16
#define STEP_DEFECT            17
#define STEP_PREPAREOCTAVE     18 
#define STEP_DIVIDEANNULAR     19 
#define STEP_COLORTERM         20
#define STEP_ANNULAR           21
#define STEP_INGESTOCTAVE      22
#define STEP_LOCALCALIBRATION  23
#define STEP_MAGDEPCALIBRATION 24
#define STEP_RECOVERPOINTS     25
#define STEP_FILTERMULTIPLE    26
#define STEP_ASTROMETRYSECOND  27 /*  Note: date brackets solve-field only find_astrometry2.csh */
#if 0
#define STEP_CLEANINTERMEDIATE xx /* No 'seconds' output and no date information -- ignore this one */
#endif
#define STEP_XFER2             28
#define STEP_SEARCHCLOSE3      29
#define MAX_STEP               30  /* Can not be greater than 31 without extending the stepsSeen bitmap */

/* The following steps provide their timings directly from the programs executed and are not included in the parsed "date" intervals */
#define CHARGE_STEP_MASK ((1 << STEP_SEARCHCLOSE) | (1 << STEP_SEXTRACTOR1) | (1 << STEP_SEXTRACTOR2) | (1 << STEP_SEXTRACTOR2B) | (1 << STEP_XFER1) | (1 << STEP_XFER2) | (1 << STEP_MAGDEPCALIBRATION))

#define REPROCESS_MASK ((1 << STEP_XFER1) | (1 << STEP_SEARCHCLOSE) | (1 << STEP_SEXTRACTOR2) | (1 << STEP_WEDGE) | (1 << STEP_MATCH2) | (1 << STEP_DEFECT) | (1 << STEP_PREPAREOCTAVE) | (1 << STEP_DIVIDEANNULAR) | (1 << STEP_COLORTERM) | (1 << STEP_ANNULAR) | (1 << STEP_INGESTOCTAVE) | (1 << STEP_LOCALCALIBRATION) | (1 << STEP_MAGDEPCALIBRATION) | (1 << STEP_RECOVERPOINTS) | (1 << STEP_FILTERMULTIPLE) | (1 << STEP_XFER2))
typedef struct _steptable {
	int sequence;
	int singleFlag;   /* This step is executed as a single thread */
	int dateFlag;     /* This file has a "date" command at the start of its script 1 = start; 2 = start and end */
	int secondsFlag;  /* This step acquires the time from a "seconds" string */
	int fullFlag;     /* This step resides only in full.csh or fulldb.csh (HEADER and ASTROMETRY can be confused with other steps ) */
	int stepmask;
	char *stepname;	char *stepdescription;
} STEPTABLE, *PSTEPTABLE;


STEPTABLE stepTable[] = {
	/*  1 d s f mask                     name                     description */
	{ 0,0,0,0,0,STEP_UNKNOWN           ,"STEP_UNKNOWN"          ,"Unrecorded time"},
	{ 1,1,0,0,0,STEP_MOSAICS           ,"STEP_MOSAICS"          ,"Build mosaic"},
	{ 2,1,0,0,0,STEP_ASTROMETRYNET     ,"STEP_ASTROMETRYNET"    ,"Run astrometry.net"},
	{ 3,1,0,0,0,STEP_IMWCS             ,"STEP_IMWCS"            ,"Run AstrometryWCS"},
	{ 4,0,0,1,1,STEP_XFER1             ,"STEP_XFER1"            ,"Transfer Starbase files to the compute node"},
	{ 5,0,1,0,1,STEP_PREPAREMOSAIC     ,"STEP_PREPAREMOSAIC"    ,"Prepare the mosaic for SCAMP fitting"},
	{ 6,0,1,0,1,STEP_FULLSEXTRACTOR    ,"STEP_FULLSEXTRACTOR"   ,"Run full sextractor on the mosaics"},
	{ 7,0,1,0,1,STEP_SEXTRACTOR1       ,"STEP_SEXTRACTOR1"      ,"Run sextractor (update only) on the mosaics"},
	{ 8,0,0,0,1,STEP_SEARCHCLOSE       ,"STEP_SEARCHCLOSE"      ,"Look for multiple images with small offsets and look for calibration mark"},
	{ 9,0,1,0,1,STEP_SCAMPFIRST        ,"STEP_SCAMPFIRST"       ,"Run SCAMP to obtain a polynomial fit"},
	{10,0,1,0,1,STEP_SEXTRACTOR2       ,"STEP_SEXTRACTOR2"      ,"Run SExtractor on the SCAMP polynomial fit"},
	/*  1 d s f mask                     name                     description */
	{11,0,0,0,0,STEP_SCAMPASTROMETRY   ,"STEP_SCAMPASTROMETRY"  ,"Run astrometry.net on the scamp results"},
	{12,0,1,0,1,STEP_SCAMPSECOND       ,"STEP_SCAMPSECOND"      ,"Run SCAMP a second time to improve the fit for patrol plates (USNO-B data)"},
	{13,0,1,0,1,STEP_SCAMPTHIRD        ,"STEP_SCAMPTHIRD"       ,"Run SCAMP a third time to improve the fit of the patrol plates (UCAC3 data)"},
	{14,0,1,0,1,STEP_SEXTRACTOR2B      ,"STEP_SEXTRACTOR2B"     ,"Run sextractor (update only) on the second SCAMP polynomial fit"},
#if 0
	{xx,0,0,0,0,STEP_HEADER            ,"STEP_HEADER"           ,"Save the header"},
#endif
	{15,0,2,0,1,STEP_WEDGE             ,"STEP_WEDGE"            ,"Find Pickering Wedge objects"},
	{16,0,2,0,1,STEP_MATCH2            ,"STEP_MATCH2"           ,"Match sextractor results with the calibration catalog"},
	{17,0,2,0,1,STEP_DEFECT            ,"STEP_DEFECT"           ,"Filter out plate defects"},
	{18,0,2,0,1,STEP_PREPAREOCTAVE     ,"STEP_PREPAREOCTAVE"    ,"Reformat results for octave routines"},
	{19,0,2,0,1,STEP_DIVIDEANNULAR     ,"STEP_DIVIDEANNULAR"    ,"Run Octave to divide stars into annular bins"},
	{20,0,2,0,1,STEP_COLORTERM         ,"STEP_COLORTERM"        ,"Run Octave to find the colorterm correction for each annular bin"},
	/*  1 d s f mask                     name                     description */
	{21,0,2,0,1,STEP_ANNULAR           ,"STEP_ANNULAR"          ,"Run Octave to find the calibration curve for each bin"},
	{22,0,0,0,1,STEP_INGESTOCTAVE      ,"STEP_INGESTOCTAVE"     ,"Reformat the results into starbase format"},
	{23,0,1,0,1,STEP_LOCALCALIBRATION  ,"STEP_LOCALCALIBRATION" ,"Perform local calibration"},
	{24,0,0,0,1,STEP_MAGDEPCALIBRATION ,"STEP_MAGDEPCALIBRATION","Perform magnitude-dependent calibration"},
	{25,0,2,0,1,STEP_RECOVERPOINTS     ,"STEP_RECOVERPOINTS"    ,"Calculate Magnitudes for all Sextractor Objects"},
	{26,0,2,0,1,STEP_FILTERMULTIPLE    ,"STEP_FILTERMULTIPLE"   ,"Extract unmatched objects for multiple exposure processing"},
	{27,0,0,0,0,STEP_ASTROMETRYSECOND  ,"STEP_ASTROMETRYSECOND" ,"Use astrometry.net to search for multiple exposures"},
#if 0
	{xx,0,0,0,1,STEP_CLEANINTERMEDIATE ,"STEP_CLEANINTERMEDIATE","Remove intermediate files"},
#endif
	{28,0,0,1,1,STEP_XFER2             ,"STEP_XFER2"            ,"Transfer Starbase files to the storage node"},
	{29,0,0,1,1,STEP_SEARCHCLOSE3      ,"STEP_SEARCHCLOSE3"     ,"Look for plates that have shifted during scanning"},

	/*  1 d s f mask                     name                     description */
};
int stepTableLength = sizeof(stepTable)/sizeof(STEPTABLE);

typedef struct _stepstats {
	time_t stepTime;
	int stepCount;
	int startedFlag;
} STEPSTATS,*PSTEPSTATS;

#define DEFAULT_ELLIPTICITY_LIMIT 1.0
#define MAX_SOLVE_FIELD 8
#define ASTROMETRYNETSERIES "ai"
#define MAX_ASTROMETRYNET_SECONDS 1000

typedef struct _astrometrynetstats {
	char series[MAX_SERIES_STRING];
	int plateNumber;
  int mosaicNumber;
  int solveStep; /* solve-field attempt */
  int width;
  int height;
  int num_out;  /* db2xyls output */
  int filterSextractorErrors; /* If true, skip SExtractor errors. */
  int astrometrynetResult;
  double ellipticityLimit; /* Ellipticity limit (default 1.0) */
  int cpu_limit;
} ASTROMETRYNETSTATS,*PASTROMETRYNETSTATS;

typedef struct _astrometrynettotal {
  int solveStep; /* solve-field attempt */
  int num_out;  /* db2xyls output */
  int filterSextractorErrors; /* If true, skip SExtractor errors. */
  double ellipticityLimit; /* Ellipticity limit (default 1.0) */
  int cpu_limit;
  size_t totalPlates;
  size_t totalSeconds;
  double averageSeconds;
} ASTROMETRYNETTOTAL,*PASTROMETRYNETTOTAL;


typedef struct _sizestats {
	int forwardCount;
	size_t aveForwardBytesCopied;
	size_t maxForwardBytesCopied;
	size_t aveForwardBytesPerSecond;
	size_t maxForwardBytesPerSecond;
	int reverseCount;
	size_t aveReverseBytesCopied;
	size_t maxReverseBytesCopied;
	size_t aveReverseBytesPerSecond;
	size_t maxReverseBytesPerSecond;
} SIZESTATS,*PSIZESTATS;


typedef struct _masksector {
  long long sectorFlag;
  int maskSector;
  int tileOverlapCount;
  int convolutionArea; /* Area in pixels of the convolution function */
  double count_med1;
  double count_rms1;
  double count_limit1; /* Original limit for FindMasks() */
  int count_limit2; /* Limit actually used for FindMasks() */
  int maxBinCount0; /* For convolution radius = 0 */
  int maxBinCount1; /* For convolution radius = 1 */
  double maxBinSNR1;
  int totalMasksSector;
  int sectorCount;
  int dimSelectedCount;
  int searchTableSize;
  int curClipCount;
  int maxMaskArea;
  int maxCenterMaskArea;
  double MAG_ISO_med;
  double MAG_ISO_RMS;
} MASKSECTOR,*PMASKSECTOR;

#define MAX_SECTOR_COUNT 20

typedef struct _searchclose3stats {
	char series[MAX_SERIES_STRING];
	int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int seconds;
  double threshold;
  int patternID;
  int maskCenterFlag;
  int multipleMaskSector; /* old search_close.c found this many masks */
  int multipleMaskFlag;   /* old search_close.c found this many normal masks */
  int multipleMaskCenterFlag; /* old search_close.c found this many masks at the center */

  double dual_count_med1;
  double dual_count_rms1;
  double dual_count_limit1;
  int dual_maxBinCount1; /* For convolution radius = 1 */
  int dual_minBinCount1; /* For convolution radius = 1 */
  double dual_maxBinSNR1;
  int dual_totalMasksSector;
  int dual_dimSelectedCount;
  int dual_searchTableSize;
  int dual_count_limit2;
  int dual_maxMaskArea;
  int dual_maxCtrMaskArea;

  double single_count_med1;
  double single_count_rms1;
  double single_count_limit1;  
  int single_maxBinCount1; /* For convolution radius = 1 */
  double single_maxBinSNR1;
  int single_totalMasksSector;
  int single_dimSelectedCount;
  int single_searchTableSize;
  int single_count_limit2;
  int single_maxMaskArea;
  int single_maxCtrMaskArea;
  int sequence;

  int totalMasks;
  int maxMaskSector;
  int minTileOverlapCount;
  int maxTileOverlapCount;
  MASKSECTOR maskSectorTable[MAX_SECTOR_COUNT];
} SEARCHCLOSE3STATS,*PSEARCHCLOSE3STATS;

#define SEARCHCLOSE3_RESULT_UNITIALIZED  0
#define SEARCHCLOSE3_RESULT_INPROGRESS   1
#define SEARCHCLOSE3_RESULT_ABORT        2
#define SEARCHCLOSE3_RESULT_COMPLETE     3
#define SEARCHCLOSE3_RESULT_RECORDED     4

#define SEARCHCLOSE3_MAX_PLATES       8000

typedef struct _platestats {
	time_t totalTime; /* Time reported by the CHILD Exited command */
	time_t beginTime; /* First timestamp */
	time_t endTime; /* Last timestamp */
	time_t chargeTime; /* Time charged to a step sharing a date block */
	time_t totalChargeTime;
  int mosaicNumber; /* Used only for astronometrynetFlag == 1 */
	int solutionNumber;
	int catalogNumber;
	int prevStep;
	int curStep;
	int stepsSeen;             /* Bitmap of steps seen so far */
	int singleFlag;            /* Single thread step has been seen */
	int updateSextractorStep;  /* Current update_sextractor mode, since update_sextractor does not generate a new step */
	int fullScript; /* Set when "copypipeline" is encountered for a full.csh script */
	int xfer1Complete; /* Set to distinguish COPY_XFER1 from COPY_XFER2 */
	char series[MAX_SERIES_STRING];
	int plateNumber;
	STEPSTATS stepStats[MAX_STEP];
  int astrometrynetResult; /* 0 = unknown, 1 = passed, -1 = failed */
  int astrometrynetSolveStep;
  ASTROMETRYNETSTATS astrometrynetstatstable[MAX_SOLVE_FIELD];
  int searchclose3Result; /* 0 = unitialized; 1 = in progress; 2 = abort; 3 complete */
  SEARCHCLOSE3STATS searchclose3statstable;
} PLATESTATS,*PPLATESTATS;
#define SEQUENCE_NONE         0
#define SEQUENCE_OK           1
#define SEQUENCE_RESCAN       2
#define SEQUENCE_INCOMPLETE   3 /* Failed to obtain a getfits image */
#define SEQUENCE_MULTIPLE     4 /* Flagged by search_close.c as a multiple exposure image */
#define SEQUENCE_MASK         5 /* Nonzero maskCount or maskCenterFlag for any zone */
#define SEQUENCE_INSUFFICIENT 6 /* Insufficient points (SEARCH_TABLE_LENGTH) for any zone */

typedef struct _platelistentry {
	char series[MAX_SERIES_STRING];
	int plateNumber;
  int rescanFlag;
  int okFlag;
  int errorFlag;
} PLATELISTENTRY,*PPLATELISTENTRY;




typedef struct _pparsecommon {
	FILE *out_handle;
	time_t startTime;
	time_t totalStepTime;
	time_t totalTime;
	time_t totalChargeTime; 
	time_t totalTimestampTime;
	int verbose;
	int reprocessFlag;
  int astrometrynetFlag; /* Special case for astrometry.net */
  int astrometrynetActive; /* astrometry.net phase detected */
  int searchclose3Flag;    /* Special case for checking shifts in small plates during scanning */
  int searchclose3Active;  /* search_close3 results being checked */
  int searchclose3TotalPlates;  /* total plates processed */
  int searchclose3MaskPlates;   /* plates that have masks */
  int searchclose3SelectedPlates; /* plates selected (dual_totalMasksSector > single_totalMasksSector) */
  int searchclose3MaxSearchTableSize; /* the maximum search table length found on all plates */
  int searchclose3MinSearchTableSize; /* the minimum search table length found in a dual overlap region */
	int maxLineLen;
	PLATESTATS plateStats;
	STEPSTATS totalStepStats[MAX_STEP];
	SIZESTATS totalXferStats[MAX_SOLUTION_NUMBER];
	SIZESTATS totalMatchstarsStats;
	SIZESTATS totalFilterblendedStats;
	SIZESTATS totalRecoverPointsStats;
  ASTROMETRYNETTOTAL allSuccessAstrometrynetTable[MAX_SOLVE_FIELD]; /* All plates */
  ASTROMETRYNETTOTAL allFailedAstrometrynetTable[MAX_SOLVE_FIELD];
  ASTROMETRYNETTOTAL seriesSuccessAstrometrynetTable[MAX_SOLVE_FIELD];  /* plates in ASTROMETRYNETSERIES */
  ASTROMETRYNETTOTAL seriesFailedAstrometrynetTable[MAX_SOLVE_FIELD];
  int allSuccessAstrometrynetSeconds[MAX_ASTROMETRYNET_SECONDS];
  int allFailedAstrometrynetSeconds[MAX_ASTROMETRYNET_SECONDS];
  int seriesSuccessAstrometrynetSeconds[MAX_ASTROMETRYNET_SECONDS];
  int seriesFailedAstrometrynetSeconds[MAX_ASTROMETRYNET_SECONDS];
  int plateListEntryAlloc;
  int plateListEntryCount;
  int plateStatsEntryCount;
  PPLATELISTENTRY pPlateListEntryTable;
  PPLATESTATS pPlateStatsTable;
} PARSECOMMON,*PPARSECOMMON;

char *weekdays[] = {"Sun ","Mon ","Tue ","Wed ","Thu ","Fri ","Sat "};
char numWeekdays = sizeof(weekdays)/sizeof(char *);
char *monthstr[] = {"Jan ","Feb ","Mar ","Apr ","May ","Jun ","Jul ","Aug ","Sep ","Oct ","Nov ","Dec "};
char numMonthstr = sizeof(monthstr)/sizeof(char *);
void	SetNewStep(PPARSECOMMON pParseCommon,int stepNumber,int lineCount,int lineNumber)
{
	PPLATESTATS pPlateStats = &pParseCommon->plateStats;
  PSTEPTABLE pStepTable;
	PSTEPSTATS pStepStats;
	if ((pPlateStats->curStep != STEP_UNKNOWN) &&
			(pPlateStats->curStep != STEP_MAGDEPCALIBRATION)) {
		pPlateStats->prevStep = pPlateStats->curStep;
	}
	pPlateStats->stepsSeen |= (1 << stepNumber);
	pPlateStats->curStep = stepNumber;
	pStepTable = &stepTable[pPlateStats->curStep];
	pStepStats = &pPlateStats->stepStats[pPlateStats->curStep];
	pStepStats->startedFlag++;
	if (pStepTable->singleFlag != 0) {
		pPlateStats->singleFlag = 1;
	}
	if (pStepTable->fullFlag) {
		pPlateStats->fullScript = 1; /* We are running in full.csh or fulldb.csh */
		if (pPlateStats->curStep != STEP_XFER1) {
			pPlateStats->xfer1Complete = 1;
		}
	}
	if (pParseCommon->verbose) {
		printf("SetNewStep step %2d %22s in line %8d parsecommon.c line %d\n",pPlateStats->curStep,pStepTable->stepname,lineCount,lineNumber);
	}
}
void SetStepTime(PPARSECOMMON pParseCommon,time_t stepTime,time_t chargeTime,int overwriteFlag,int lineCount,int lineNumber) 
{
	PSTEPSTATS pStepStats;
	PSTEPTABLE pStepTable;
	PPLATESTATS pPlateStats = &pParseCommon->plateStats;
	
	pPlateStats->totalChargeTime += chargeTime;
	pStepStats = &pPlateStats->stepStats[pPlateStats->curStep];
	pStepTable = &stepTable[pPlateStats->curStep];
	pStepStats->stepCount++;
	if (overwriteFlag) {
		pStepStats->stepTime = stepTime;
	} else {
		pStepStats->stepTime += stepTime;
	}
	if (pPlateStats->curStep == STEP_UNKNOWN) {
		printf("ERROR: line %4d SetStepTime  time %lld seconds chargeTime %5lld overwriteFlag %d step %2d %22s in line %8d parsecommon.c lineCount %d\n",__LINE__,stepTime,chargeTime,overwriteFlag,pPlateStats->curStep,pStepTable->stepname,lineCount,lineNumber);

	} else {
		if (pParseCommon->verbose) {
			printf("SetStepTime time %5lld seconds chargeTime %5lld overwriteFlag %d step %2d %22s in line %8d parsecommon.c line %d\n",stepTime,chargeTime,overwriteFlag,pPlateStats->curStep,pStepTable->stepname,lineCount,lineNumber);
		}
	}
	
	return;
}
void DisplayUpdateSextractorStep(PPARSECOMMON pParseCommon,int lineCount,int lineNumber) {
	PPLATESTATS pPlateStats = &pParseCommon->plateStats;
	PSTEPTABLE pStepTable;
	PSTEPTABLE pUpdateStepTable;
	pStepTable = &stepTable[pPlateStats->curStep];
	pUpdateStepTable = &stepTable[pPlateStats->updateSextractorStep];
	printf("DisplayUpdateSextractorStep for update_sextractor %s curStep %s  in line %8d parsecommon.c line %d\n",pUpdateStepTable->stepname,pStepTable->stepname,lineCount,lineNumber);
	return;
}


/*
 *  Formats time of the form: "Fri Dec 20 13:43:57 EST 2013"
 *  Script started or Script done
 */
int ParseTime(PPARSECOMMON pParseCommon,char *inLine,time_t* timeResult,int lineCount) {
	int day;
	int month;
	char *charPtr = NULL;
	char timebuf[MAX_TIME_LENGTH+1];
	time_t tp;
	struct tm tmstruct;
	struct tm *ptr;
	char timestr[100];
	int nvals;
	*timeResult = 0;
	for (day = 0; day < numWeekdays; day++) {
		charPtr = strstr(inLine,weekdays[day]);
		if (charPtr != NULL) {
			break;
		}
	}
	if (charPtr == NULL) {
		return(-1);
	}
	strncpy(timebuf,charPtr,MAX_TIME_LENGTH);
	for (month = 0; month < numMonthstr; month++) {
		charPtr = strstr(timebuf,monthstr[month]);
		if (charPtr != NULL) {
			break;
		}
	}
	if (charPtr == NULL) {
		printf("ERROR: failed to find month in %s line %d\n",timebuf,lineCount);
		return(-1);
	}

	timebuf[MAX_TIME_LENGTH] = 0;
	memset(&tmstruct,0,sizeof(struct tm));
	if (strstr(timebuf,"EST") != NULL) {
		tmstruct.tm_isdst = 0;
	} else if (strstr(timebuf,"EDT") != NULL) {
		tmstruct.tm_isdst = 1;
	} else {
		if ((strstr(inLine,"Script started") != NULL) ||
				(strstr(inLine,"Script done") != NULL)) {
			return(-1);
		} else {
			printf("ERROR: unrecognized time type in %s line %d\n",inLine,lineCount);
			return(-1);
		}
	}
	if (tmstruct.tm_isdst == 0) {
		nvals = sscanf(&timebuf[8],"%d %d:%d:%d EST %d",&tmstruct.tm_mday,&tmstruct.tm_hour,&tmstruct.tm_min,&tmstruct.tm_sec,&tmstruct.tm_year);
	} else {
		nvals = sscanf(&timebuf[8],"%d %d:%d:%d EDT %d",&tmstruct.tm_mday,&tmstruct.tm_hour,&tmstruct.tm_min,&tmstruct.tm_sec,&tmstruct.tm_year);
	}
	if (nvals != 5) {
		printf("ERROR: failed to decode numbers in %s line %d nvals %d %d %d:%d:%d yr %d\n",inLine,lineCount,nvals,tmstruct.tm_mday,tmstruct.tm_hour,tmstruct.tm_min,tmstruct.tm_sec,tmstruct.tm_year);
		return(-1);
	}
	
	tmstruct.tm_year -= 1900;
	tmstruct.tm_mon = month;
	
	tp = mktime(&tmstruct);
	ptr = localtime(&tp);
	strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
	if (day != tmstruct.tm_wday) {
		printf("ERROR: day disagreement %d %d for %s lineCount %d timestr %s\n",day,tmstruct.tm_wday,inLine,lineCount,timestr);
		return(-1);
	}
	if (pParseCommon->verbose) {
		printf("Found %s  timestr %s lineCount %d\n",inLine,timestr,lineCount);
	}
	*timeResult = tp;
	return(0);
}

void ProcessMaxMemory(PPARSECOMMON pParseCommon,char * inLine,char * log_name) {
	PSIZESTATS pSizeStats;
	char *charPtr;
	int nvals;
	size_t tmpBytesCopied;
	charPtr = strstr(inLine,"properMotionErrLimCount");
	if (charPtr != NULL) {
		pSizeStats = &pParseCommon->totalMatchstarsStats;
	} else  {
		charPtr = strstr(inLine,"Blended stars");
		if (charPtr != NULL) {
			pSizeStats = &pParseCommon->totalFilterblendedStats;
		} else {
			charPtr = strstr(inLine,"Illegal");
			if (charPtr != NULL) {
				pSizeStats = &pParseCommon->totalRecoverPointsStats;
			} else {
				printf("ERROR: unknown maxMemory in line %s of %s\n",inLine,log_name);
				return;
			}
		}
	}
	charPtr = strstr(inLine,"maxMemory");
	if (charPtr == NULL) {
		printf("ERROR: illegal call to ProcessMaxMemory\n");
		exit(-1);
	}
	nvals = sscanf(charPtr,"maxMemory %lld",&tmpBytesCopied);
	if (nvals == 1) {
		pSizeStats->forwardCount++;
		pSizeStats->aveForwardBytesCopied += tmpBytesCopied;
		if (pSizeStats->maxForwardBytesCopied < tmpBytesCopied) {
			pSizeStats->maxForwardBytesCopied = tmpBytesCopied;
		}
	} else {
		printf("ERROR: ProcessMaxMemory can not parse line %s of %s\n",inLine,log_name);
	}

	return;
	
}

int ParseAstrometryNetParameters(PPARSECOMMON pParseCommon,char* copyLine) {
  char inLine[MAX_BUFFER];
  char *charPtr1;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int nvals;
	PPLATESTATS pPlateStats = &pParseCommon->plateStats;
  PASTROMETRYNETSTATS pAstrometrynetStats;
  if (strlen(copyLine) >= MAX_BUFFER) {
    printf("ERROR: ParseAstrometryNetParameters code line %d MAX_BUFFER exceeded with length %d\n",__LINE__,strlen(copyLine));
  }
#if 0
  printf("ParseAstrometryNetParameters code line %d copyLine '%s'\n",__LINE__,copyLine);
#endif
  strcpy(inLine,copyLine);
  charPtr1 = strtok(inLine," ");
  if (charPtr1 == NULL) {
    printf("ERROR line %d for %s\n",__LINE__,copyLine);
    return(-1);
  }
  if (strcmp(charPtr1,"db2xyls") == 0) {
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"-v") == 0) {
       charPtr1 = strtok(NULL," ");
    }

    if (strcmp(charPtr1,"-i") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (ParseFilename(charPtr1,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if ((strcmp(series,pPlateStats->series) != 0) ||
        (plateNumber != pPlateStats->plateNumber)) {
      printf("ERROR line %d plate %s%05d disagrees with  %s\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,copyLine);
      return(-1);
    }
    pPlateStats->mosaicNumber = mosaicNumber;
    if (pPlateStats->astrometrynetSolveStep >= MAX_SOLVE_FIELD) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    pAstrometrynetStats = &pPlateStats->astrometrynetstatstable[pPlateStats->astrometrynetSolveStep];
    pAstrometrynetStats->ellipticityLimit = DEFAULT_ELLIPTICITY_LIMIT;
    strcpy(pAstrometrynetStats->series,pPlateStats->series);
    pAstrometrynetStats->plateNumber = pPlateStats->plateNumber;
    pAstrometrynetStats->mosaicNumber = pPlateStats->mosaicNumber;
    pPlateStats->astrometrynetSolveStep++;
    pAstrometrynetStats->solveStep = pPlateStats->astrometrynetSolveStep; 
 
    while (1) {
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        break;
      }
      if (strcmp(charPtr1,"-w") == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"-h") == 0) {
          /* Here everything is going to fail.  Clean up and exit */
          return(-1);
        }

        nvals = sscanf(charPtr1,"%d",&pAstrometrynetStats->width);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      } else if (strcmp(charPtr1,"-h") == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pAstrometrynetStats->height);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      } else if (strcmp(charPtr1,"-n") == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pAstrometrynetStats->num_out);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      } else if (strcmp(charPtr1,"-e") == 0) {
        pAstrometrynetStats->filterSextractorErrors = 1;
      } else if (strcmp(charPtr1,"-l") == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%lf",&pAstrometrynetStats->ellipticityLimit);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      } 
    }



  } else if (strcmp(charPtr1,"solve-field") == 0) {
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (ParseFilename(charPtr1,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if ((strcmp(series,pPlateStats->series) != 0) ||
        (plateNumber != pPlateStats->plateNumber)) {
      printf("ERROR line %d plate %s%05d disagrees with  %s\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,copyLine);
      return(-1);
    }

    if ((pPlateStats->astrometrynetSolveStep <= -1) ||
        (pPlateStats->astrometrynetSolveStep > (MAX_SOLVE_FIELD))) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    pAstrometrynetStats = &pPlateStats->astrometrynetstatstable[pPlateStats->astrometrynetSolveStep-1];
    if ((strcmp(pAstrometrynetStats->series,series) != 0) ||
        (pAstrometrynetStats->plateNumber != plateNumber) ||
        (pAstrometrynetStats->mosaicNumber != mosaicNumber)) {
      printf("ERROR line %d plate %s%05d disagrees with  %s\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,copyLine);
      return(-1);
    }
    

    while (1) {
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        break;
      }
      if (strcmp(charPtr1,"--cpulimit") == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pAstrometrynetStats->cpu_limit);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      } 
    }



  } else {
    printf("ERROR line %d for %s\n",__LINE__,copyLine);
    return(-1);
  }


  
  return(0);
}
void PrintAstrometryNetStats(PPARSECOMMON pParseCommon,int seriesFlag) {
  ASTROMETRYNETTOTAL successAstrometrynetTable;
  ASTROMETRYNETTOTAL failedAstrometrynetTable;
  PASTROMETRYNETTOTAL pSuccessAstrometrynetTable;
  PASTROMETRYNETTOTAL pFailedAstrometrynetTable;
  PASTROMETRYNETTOTAL pSuccessAstrometrynetTotal;
  PASTROMETRYNETTOTAL pFailedAstrometrynetTotal;
  int *pSuccessAstrometrynetSeconds;
  int *pFailedAstrometrynetSeconds;
  int solveStep;
  int baseSequence; /* sequence = 1 all Success; 2 all failed; 3 series Success; 4 series failed */
  int secondsIndex;
  int cumulativeSuccess = 0;
  int cumulativeFailed = 0;
  int cumulativeSuccessTotal = 0;
  int cumulativeFailedTotal = 0;

  memset(&successAstrometrynetTable,0,sizeof(ASTROMETRYNETTOTAL));
  memset(&failedAstrometrynetTable,0,sizeof(ASTROMETRYNETTOTAL));
  if (seriesFlag == 0) {
    printf("astrometry.net results for all plates\n");
    pSuccessAstrometrynetTable = pParseCommon->allSuccessAstrometrynetTable;
    pFailedAstrometrynetTable = pParseCommon->allFailedAstrometrynetTable;
    pSuccessAstrometrynetSeconds = pParseCommon->allSuccessAstrometrynetSeconds;
    pFailedAstrometrynetSeconds = pParseCommon->allFailedAstrometrynetSeconds;
    baseSequence = 0;
  } else {
    printf("astrometrynet results for '%s' plates\n",ASTROMETRYNETSERIES);
    pSuccessAstrometrynetTable = pParseCommon->seriesSuccessAstrometrynetTable;
    pFailedAstrometrynetTable = pParseCommon->seriesFailedAstrometrynetTable;
    pSuccessAstrometrynetSeconds = pParseCommon->seriesSuccessAstrometrynetSeconds;
    pFailedAstrometrynetSeconds = pParseCommon->seriesFailedAstrometrynetSeconds;
    baseSequence = 2;
    
  }

  for (secondsIndex = 0; secondsIndex < MAX_ASTROMETRYNET_SECONDS; secondsIndex++) {
     cumulativeSuccessTotal += pSuccessAstrometrynetSeconds[secondsIndex];
     cumulativeFailedTotal += pFailedAstrometrynetSeconds[secondsIndex];
  }

  for (secondsIndex = 0; secondsIndex < MAX_ASTROMETRYNET_SECONDS; secondsIndex++) {
    if (pSuccessAstrometrynetSeconds[secondsIndex] != 0) {
      cumulativeSuccess += pSuccessAstrometrynetSeconds[secondsIndex];
      fprintf(pParseCommon->out_handle,"%d\t%d\t%f\t%d\n",secondsIndex,pSuccessAstrometrynetSeconds[secondsIndex],(1.0*cumulativeSuccess)/(1.0*cumulativeSuccessTotal),baseSequence+1);
    }
    if (pFailedAstrometrynetSeconds[secondsIndex] != 0) {
      cumulativeFailed += pFailedAstrometrynetSeconds[secondsIndex];
      fprintf(pParseCommon->out_handle,"%d\t%d\t%f\t%d\n",secondsIndex,pSuccessAstrometrynetSeconds[secondsIndex],(1.0*cumulativeFailed)/(1.0*cumulativeFailedTotal),baseSequence+2);
    }
  }

  for (solveStep = 0; solveStep <= MAX_SOLVE_FIELD; solveStep++) {
    if (solveStep < MAX_SOLVE_FIELD) {
      pSuccessAstrometrynetTotal = &pSuccessAstrometrynetTable[solveStep];
      pFailedAstrometrynetTotal = &pFailedAstrometrynetTable[solveStep];
      successAstrometrynetTable.totalPlates += pSuccessAstrometrynetTotal->totalPlates;
      successAstrometrynetTable.totalSeconds += pSuccessAstrometrynetTotal->totalSeconds;
      failedAstrometrynetTable.totalPlates += pFailedAstrometrynetTotal->totalPlates;
      failedAstrometrynetTable.totalSeconds += pFailedAstrometrynetTotal->totalSeconds;
      
    } else { 
      pSuccessAstrometrynetTotal = &successAstrometrynetTable;
      pFailedAstrometrynetTable = &failedAstrometrynetTable;
    }

    if ((pSuccessAstrometrynetTotal->totalPlates != 0)  ||
        (pFailedAstrometrynetTotal->totalPlates != 0)) {
      if (pSuccessAstrometrynetTotal->totalPlates != 0) {
        pSuccessAstrometrynetTotal->averageSeconds = (1.0 * pSuccessAstrometrynetTotal->totalSeconds)/(1.0 * pSuccessAstrometrynetTotal->totalPlates);
      } else {
        pSuccessAstrometrynetTotal->averageSeconds = 0;
      }
      if (pFailedAstrometrynetTotal->totalPlates != 0) {
        pFailedAstrometrynetTotal->averageSeconds =  (1.0 * pFailedAstrometrynetTotal->totalSeconds)/(1.0 * pFailedAstrometrynetTotal->totalPlates);
      } else {
        pFailedAstrometrynetTotal->averageSeconds = 0;
      }

      printf("solveStep %d success %6d ave sec/plate %7.1f    failed %6d ave sec/plate %7.1f \n",solveStep+1,pSuccessAstrometrynetTotal->totalPlates,pSuccessAstrometrynetTotal->averageSeconds,pFailedAstrometrynetTotal->totalPlates,pFailedAstrometrynetTotal->averageSeconds);

    } 
  }





}
void UpdateAstrometryNetStats(PPARSECOMMON pParseCommon,PASTROMETRYNETTOTAL pAstrometrynetTable,int *pAstrometrynetSeconds)
{
  int solveStep;
	PPLATESTATS pPlateStats = &pParseCommon->plateStats;
  PASTROMETRYNETSTATS pAstrometrynetStats;
  PASTROMETRYNETTOTAL pAstrometrynetTotal;
  int secondsIndex;
  for (solveStep = 0; solveStep < MAX_SOLVE_FIELD; solveStep++) {
    pAstrometrynetStats = &pPlateStats->astrometrynetstatstable[solveStep];
    if (pAstrometrynetStats->solveStep == 0) {
      break;
    }
    if (pAstrometrynetStats->solveStep != (solveStep+1)) {
      printf("ERROR: UpdateAstrometryNetStats code line %d\n",__LINE__);
      exit(-1);
    }
    pAstrometrynetTotal = &pAstrometrynetTable[solveStep];
    if (pAstrometrynetTotal->solveStep == 0) {
      pAstrometrynetTotal->solveStep = pAstrometrynetStats->solveStep;
      pAstrometrynetTotal->num_out = pAstrometrynetStats->num_out;
      pAstrometrynetTotal->filterSextractorErrors = pAstrometrynetStats->filterSextractorErrors;
      pAstrometrynetTotal->ellipticityLimit = pAstrometrynetStats->ellipticityLimit;
      pAstrometrynetTotal->cpu_limit = pAstrometrynetStats->cpu_limit;
    } else {
      if (pAstrometrynetTotal->solveStep != pAstrometrynetStats->solveStep) {
        printf("ERROR: UpdateAstrometryNetStats code line %d\n",__LINE__);
        exit(-1);
      }
      if (pAstrometrynetTotal->num_out != pAstrometrynetStats->num_out) {
        printf("ERROR: UpdateAstrometryNetStats code line %d\n",__LINE__);
        exit(-1);
      }
      if (pAstrometrynetTotal->filterSextractorErrors != pAstrometrynetStats->filterSextractorErrors) {
        printf("ERROR: UpdateAstrometryNetStats code line %d\n",__LINE__);
        exit(-1);
      }
      if (pAstrometrynetTotal->ellipticityLimit != pAstrometrynetStats->ellipticityLimit) {
        printf("ERROR: UpdateAstrometryNetStats code line %d\n",__LINE__);
        exit(-1);
      }
      if (pAstrometrynetTotal->cpu_limit != pAstrometrynetStats->cpu_limit ) {
        printf("ERROR: UpdateAstrometryNetStats code line %d\n",__LINE__);
        exit(-1);
      }
    }
  }
  if(solveStep == 0) {
    printf("ERROR: UpdateAstrometryNetStats code line %d\n",__LINE__);
    exit(-1);
  }
  pAstrometrynetTotal = &pAstrometrynetTable[solveStep-1];
  pAstrometrynetTotal->totalPlates++;
  pAstrometrynetTotal->totalSeconds += pPlateStats->totalTime;
  secondsIndex = pPlateStats->totalTime;
  if (secondsIndex < 0) {
    secondsIndex = 0;
  } else if (secondsIndex >= MAX_ASTROMETRYNET_SECONDS) {
    secondsIndex = MAX_ASTROMETRYNET_SECONDS-1;
  }
  pAstrometrynetSeconds[secondsIndex]++;

}
int ParseSearchClose3Parameters(PPARSECOMMON pParseCommon,int lineCount,char* copyLine) {
  char inLine[MAX_BUFFER];
  char *charPtr1 = NULL;
  char *charPtr2 = NULL;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int nvals;
  int maskSector;
  int maskArea;
  int maskCenterFlag;
	PPLATESTATS pPlateStats = &pParseCommon->plateStats;
  PSEARCHCLOSE3STATS pSearchClose3Stats = &pPlateStats->searchclose3statstable;
  PMASKSECTOR pMaskSector = NULL;
  if (strlen(copyLine) >= MAX_BUFFER) {
    printf("ERROR: ParseSearchClose3Parameters code line %d MAX_BUFFER exceeded with length %d\n",__LINE__,strlen(copyLine));
  }

#if 0
  printf("ParseSearchClose3Parameters code line %d copyLine '%s'\n",__LINE__,copyLine);
#endif
  strcpy(inLine,copyLine);

  charPtr2 = strstr(inLine,"search_close3 sector");
  if (charPtr2 == NULL) {
    charPtr2 = strstr(inLine,"search_close3 Sector");
  }
  if (charPtr2 != NULL) {
    charPtr2 += strlen("search_close3 sector");
    charPtr1 = strtok(charPtr2," ");

    /*  Here we are parsing "line xxxx search_close3 sector  0 overlap  7" */  

    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    nvals = sscanf(charPtr1,"%d",&maskSector);
    if (nvals != 1) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if ((maskSector < 0) ||
        (maskSector >= MAX_SECTOR_COUNT)) {
      printf("ERROR line %d for %s bad MAX_SECTOR_COUNT\n",__LINE__,copyLine);
      return(-1);
    }
    pMaskSector = &pSearchClose3Stats->maskSectorTable[maskSector];
    pMaskSector-> maskSector = maskSector;
    if (maskSector > pSearchClose3Stats->maxMaskSector) {
      pSearchClose3Stats->maxMaskSector = maskSector;
    }

    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"overlap") != 0) {
#if 0 /* do not confuse with "line 2617 search_close3 Sector  0 Min, Max XY: 0.000000 2316.000000 aveXY 1158.000000" */
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
#endif
      return(0);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    nvals = sscanf(charPtr1,"%d",&pMaskSector->tileOverlapCount);
    if (nvals != 1) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }

    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"convolutionArea") == 0) {

      /*  Here we are parsing "line 3023 search_close3 sector  0 overlap  7 convolutionArea 9 count_med1 55.000000 count_rms1 15.143754, count_limit1 282.000000 maxBinCount1 144 maxBinSNR1 5.877010 for b05498_00_01r270ww" */
      if (pPlateStats->searchclose3Result != SEARCHCLOSE3_RESULT_INPROGRESS) {     
        printf("ERROR line %d unrecognized searchclose3Result for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%d",&pMaskSector->convolutionArea);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"count_med1") != 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%lf",&pMaskSector->count_med1);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"count_rms1") != 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%lf",&pMaskSector->count_rms1);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"count_limit1") != 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%lf",&pMaskSector->count_limit1);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"maxBinCount1") != 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%d",&pMaskSector->maxBinCount1);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"maxBinSNR1") != 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%lf",&pMaskSector->maxBinSNR1);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"for") != 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (ParseFilename(charPtr1,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      /* Here we ignore the mosaic number, because it was not set when the child started */
      if ((strcmp(pPlateStats->series,series) != 0)  ||
          (pPlateStats->plateNumber != plateNumber)) {
        printf("ERROR line %d plate %s%05d_%02d is not %s%05d_%02d for %s\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,pPlateStats->mosaicNumber,series,plateNumber,mosaicNumber,copyLine);
        return(-1);
      }
      if (pSearchClose3Stats->series[0] == 0) {
        strcpy(pSearchClose3Stats->series,series);
        pSearchClose3Stats->plateNumber = plateNumber;
        pSearchClose3Stats->mosaicNumber = mosaicNumber;
        pSearchClose3Stats->binning = binning;
        pSearchClose3Stats->rotation = rotation;
        pPlateStats->mosaicNumber = mosaicNumber;
      } else {
        if ((strcmp(pSearchClose3Stats->series,series) != 0)  ||
            (pSearchClose3Stats->plateNumber != plateNumber) ||
            (pSearchClose3Stats->mosaicNumber != mosaicNumber) ||
            (pSearchClose3Stats->binning != binning) ||
            (pSearchClose3Stats->rotation != rotation)) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      } 

    } else if (strcmp(charPtr1,"totalMasksSector") == 0) {

      /*  Here we are parsing "line 3197 search_close3 Sector  0 overlap  7  totalMasksSector      0 totalMasks      0 for b05498_00_01r270ww" */
      /*  alternatively:      "line 3243 search_close3 Sector  0 overlap  7  totalMasksSector      3 totalMasks      3 dimSelectedCount   2495 searchTableSize 2500 curClipCount  15231 MAG_ISO_med -8.008800 MAG_ISO_RMS 0.393396 for br00594_00_01r270ww */
      /*  alternatively:      "line 3243 search_close3 Sector  0 overlap  7  totalMasksSector   1006 totalMasks   1006 sectorCount  10280 dimSelectedCount   1310 searchTableSize 1319 curClipCount   8531 MAG_ISO_med -7.093600 MAG_ISO_RMS 0.618401 for br01064_00_01r90ww" */
      /*  alternatively:      "line 3297 search_close3 Sector  0 overlap  7  totalMasksSector      0 totalMasks      0 sectorCount  55871 dimSelectedCount   1986 searchTableSize 1990 curClipCount  48330 MAG_ISO_med -8.656400 MAG_ISO_RMS 0.581355 count_limit2 18 for br00370_00_01r90ww" */
      /*  alternatively:      "line 1306 for line 3335 search_close3 Sector  0 overlap  7  totalMasksSector      3 totalMasks      3 sectorCount  21457 dimSelectedCount   1984 searchTableSize 1990 curClipCount  18913 MAG_ISO_med -8.833200 MAG_ISO_RMS 0.436182 maxMaskArea   17 maxCenterMaskArea    0 for br00745_00_01r270ww" */


	    if (pPlateStats->searchclose3Result != SEARCHCLOSE3_RESULT_INPROGRESS) {     
        printf("ERROR line %d unrecognized searchclose3Result for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%d",&pMaskSector->totalMasksSector);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"totalMasks") != 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%d",&pSearchClose3Stats->totalMasks);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"sectorCount")  == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pMaskSector->sectorCount);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      }

      if (strcmp(charPtr1,"dimSelectedCount")  == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pMaskSector->dimSelectedCount);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        

        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"searchTableSize") != 0) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pMaskSector->searchTableSize);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }

        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"curClipCount") != 0) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pMaskSector->curClipCount);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }

        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"MAG_ISO_med") != 0) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%lf",&pMaskSector->MAG_ISO_med);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }

        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"MAG_ISO_RMS") != 0) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%lf",&pMaskSector->MAG_ISO_RMS);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }

        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      }
      if (strcmp(charPtr1,"count_limit2") == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pMaskSector->count_limit2);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        

        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }

      }
      if (strcmp(charPtr1,"maxMaskArea") == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pMaskSector->maxMaskArea);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        
        if (strcmp(charPtr1,"maxCenterMaskArea") != 0) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pMaskSector->maxCenterMaskArea);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }

        if (pMaskSector->maxCenterMaskArea > 0) {
           pSearchClose3Stats->maskCenterFlag = 1;
        }

        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      }


      if (strcmp(charPtr1,"for") != 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (ParseFilename(charPtr1,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if ((strcmp(pPlateStats->series,series) != 0)  ||
          (pPlateStats->plateNumber != plateNumber) ||
          (pPlateStats->mosaicNumber != mosaicNumber)) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (pSearchClose3Stats->series[0] == 0) {
        strcpy(pSearchClose3Stats->series,series);
        pSearchClose3Stats->plateNumber = plateNumber;
        pSearchClose3Stats->mosaicNumber = mosaicNumber;
        pSearchClose3Stats->binning = binning;
        pSearchClose3Stats->rotation = rotation;
      } else {
        if ((strcmp(pSearchClose3Stats->series,series) != 0)  ||
            (pSearchClose3Stats->plateNumber != plateNumber) ||
            (pSearchClose3Stats->mosaicNumber != mosaicNumber) ||
            (pSearchClose3Stats->binning != binning) ||
            (pSearchClose3Stats->rotation != rotation)) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      } 
   
    
    } else if (strcmp(charPtr1,"initial") == 0) {
      /*  Here we are parsing "line 3051 search_close3 sector  0 overlap  7 initial count_med0 0.000000 count_rms0 0.000000, count_limit0 0 maxBinCount0 0 maxBinSNR0 0.000000 for br01064_00_01r90ww" */
      /*  alternatively       "line 3051 search_close3 sector  0 overlap  7 initial count_med0 0.000000 count_rms0 0.000000, count_limit0 0 maxBinSNR0 0.000000 for br01064_00_01r90ww" */

      if ((pPlateStats->searchclose3Result != SEARCHCLOSE3_RESULT_UNITIALIZED) && (maskSector == 0)) {     
        printf("ERROR line %d searchclose3 overwrite attempt for plate %s \n",__LINE__,copyLine);
        return(-1);
      } else if ((pPlateStats->searchclose3Result == SEARCHCLOSE3_RESULT_UNITIALIZED) && (maskSector == 0)) {     
        pPlateStats->searchclose3Result = SEARCHCLOSE3_RESULT_INPROGRESS;
      } else if (pPlateStats->searchclose3Result != SEARCHCLOSE3_RESULT_INPROGRESS) {     
        printf("ERROR line %d unrecognized searchclose3Result for %s\n",__LINE__,copyLine);
        return(-1);
      }
      while (1) {
        /* step forward to count_limit0  */
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"count_limit0") == 0) {
          break;
        }
      }
      /* Skip number */
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (strcmp(charPtr1,"maxBinCount0") == 0) {
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        nvals = sscanf(charPtr1,"%d",&pMaskSector->maxBinCount0);
        if (nvals != 1) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
  

      } 
      while (1) {
        /* step forward to for  */
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"for") == 0) {
          break;
        }
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (ParseFilename(charPtr1,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      /* Here we ignore the mosaic number, because it was not set when the child started */
      if ((strcmp(pPlateStats->series,series) != 0)  ||
          (pPlateStats->plateNumber != plateNumber)) {
        printf("ERROR line %d plate %s%05d_%02d is not %s%05d_%02d for %s\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,pPlateStats->mosaicNumber,series,plateNumber,mosaicNumber,copyLine);
        return(-1);
      }
      if (pSearchClose3Stats->series[0] == 0) {
        strcpy(pSearchClose3Stats->series,series);
        pSearchClose3Stats->plateNumber = plateNumber;
        pSearchClose3Stats->mosaicNumber = mosaicNumber;
        pSearchClose3Stats->binning = binning;
        pSearchClose3Stats->rotation = rotation;
        pPlateStats->mosaicNumber = mosaicNumber;
      } else {
        if ((strcmp(pSearchClose3Stats->series,series) != 0)  ||
            (pSearchClose3Stats->plateNumber != plateNumber) ||
            (pSearchClose3Stats->mosaicNumber != mosaicNumber) ||
            (pSearchClose3Stats->binning != binning) ||
            (pSearchClose3Stats->rotation != rotation)) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      }
    } else if (strcmp(charPtr1,"Mask") == 0) {
      /* Here we are processing "line 889 copyLine 'line  646 search_close3 Sector  2 overlap  7 Mask    1 xIndex  -77  -76 yIndex -993 -992 maskArea    1 ratio 1.000 maxBinCount    9 maskCenterFlag 0 maskCenterDistance  995.4 for br00442_00_01r90ww" */
      while (1) {
        /* step forward to maskArea  */
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"maskArea") == 0) {
          break;
        }
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%d",&maskArea);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }

      while (1) {
        /* step forward to maskCenterFlag  */
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"maskCenterFlag") == 0) {
          break;
        }
      }

      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%d",&maskCenterFlag);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (maskCenterFlag == 0) {
        if (maskArea > pMaskSector->maxMaskArea) {
          pMaskSector->maxMaskArea = maskArea;
        }
      } else {
        pSearchClose3Stats->maskCenterFlag = 1;
      }

      while (1) {
        /* step forward to for  */
        charPtr1 = strtok(NULL," ");
        if (charPtr1 == NULL) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
        if (strcmp(charPtr1,"for") == 0) {
          break;
        }
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      if (ParseFilename(charPtr1,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      /* Here we ignore the mosaic number, because it was not set when the child started */
      if ((strcmp(pPlateStats->series,series) != 0)  ||
          (pPlateStats->plateNumber != plateNumber)) {
        printf("ERROR line %d plate %s%05d_%02d is not %s%05d_%02d for %s\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber,pPlateStats->mosaicNumber,series,plateNumber,mosaicNumber,copyLine);
        return(-1);
      }
      if (pSearchClose3Stats->series[0] == 0) {
        strcpy(pSearchClose3Stats->series,series);
        pSearchClose3Stats->plateNumber = plateNumber;
        pSearchClose3Stats->mosaicNumber = mosaicNumber;
        pSearchClose3Stats->binning = binning;
        pSearchClose3Stats->rotation = rotation;
        pPlateStats->mosaicNumber = mosaicNumber;
      } else {
        if ((strcmp(pSearchClose3Stats->series,series) != 0)  ||
            (pSearchClose3Stats->plateNumber != plateNumber) ||
            (pSearchClose3Stats->mosaicNumber != mosaicNumber) ||
            (pSearchClose3Stats->binning != binning) ||
            (pSearchClose3Stats->rotation != rotation)) {
          printf("ERROR line %d for %s\n",__LINE__,copyLine);
          return(-1);
        }
      }
      
 
    }
  }
  /* HERE we are parsing "line 1998 search_close3 QUALITY_MULTIPLE maskCount  1 maskFlag 1 maskCenterFlag 0 for plate br00529_00_01r270ww" */
  charPtr2 = strstr(inLine,"search_close3 QUALITY_MULTIPLE");
  if (charPtr2 != NULL) {
    charPtr2 += strlen("search_close3 QUALITY_MULTIPLE");
    charPtr1 = strtok(charPtr2," ");


    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"maskCount") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    nvals = sscanf(charPtr1,"%d",&pSearchClose3Stats->multipleMaskSector);
    if (nvals != 1) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }

    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"maskFlag") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    nvals = sscanf(charPtr1,"%d",&pSearchClose3Stats->multipleMaskFlag);
    if (nvals != 1) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }

    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"maskCenterFlag") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    nvals = sscanf(charPtr1,"%d",&pSearchClose3Stats->multipleMaskCenterFlag);
    if (nvals != 1) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }

    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"for") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"plate") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (ParseFilename(charPtr1,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    /* Here we ignore the mosaic number, because it was not set when the child started */
    if ((strcmp(pPlateStats->series,series) != 0)  ||
        (pPlateStats->plateNumber != plateNumber)) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (pSearchClose3Stats->series[0] == 0) {
      strcpy(pSearchClose3Stats->series,series);
      pSearchClose3Stats->plateNumber = plateNumber;
      pSearchClose3Stats->mosaicNumber = mosaicNumber;
      pSearchClose3Stats->binning = binning;
      pSearchClose3Stats->rotation = rotation;
      pPlateStats->mosaicNumber = mosaicNumber;
    } else {
      if ((strcmp(pSearchClose3Stats->series,series) != 0)  ||
          (pSearchClose3Stats->plateNumber != plateNumber) ||
          (pSearchClose3Stats->mosaicNumber != mosaicNumber) ||
          (pSearchClose3Stats->binning != binning) ||
          (pSearchClose3Stats->rotation != rotation)) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
    } 
  } 



  /* Here we are parsing "threshold 2.200000 totalMasks 2105 seconds 9444 for c11004_00_01r180ww" */
  /*       alternatively "threshold 2.200000 totalMasks 101 patternID 38 seconds 14 for br00594_00_01r270ww */
  charPtr2 = strstr(inLine,"totalMasks");
  if (charPtr2 != NULL) {
    charPtr2 = strstr(inLine,"threshold");
  }
  if (charPtr2 != NULL) {
    if (pPlateStats->searchclose3Result != SEARCHCLOSE3_RESULT_INPROGRESS) {     
      printf("ERROR line %d unrecognized searchclose3Result for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr2 += strlen("threshold");
    charPtr1 = strtok(charPtr2," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    nvals = sscanf(charPtr1,"%lf",&pSearchClose3Stats->threshold);
    if (nvals != 1) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }

    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"totalMasks") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    nvals = sscanf(charPtr1,"%d",&pSearchClose3Stats->totalMasks);
    if (nvals != 1) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }

    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }

    if (strcmp(charPtr1,"patternID") == 0) {
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      nvals = sscanf(charPtr1,"%d",&pSearchClose3Stats->patternID);
      if (nvals != 1) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
      charPtr1 = strtok(NULL," ");
      if (charPtr1 == NULL) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }


    }
    if (strcmp(charPtr1,"seconds") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    nvals = sscanf(charPtr1,"%d",&pSearchClose3Stats->seconds);
    if (nvals != 1) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (strcmp(charPtr1,"for") != 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    charPtr1 = strtok(NULL," ");
    if (charPtr1 == NULL) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (ParseFilename(charPtr1,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if ((strcmp(pPlateStats->series,series) != 0)  ||
        (pPlateStats->plateNumber != plateNumber) ||
        (pPlateStats->mosaicNumber != mosaicNumber)) {
      printf("ERROR line %d for %s\n",__LINE__,copyLine);
      return(-1);
    }
    if (pSearchClose3Stats->series[0] == 0) {
      strcpy(pSearchClose3Stats->series,series);
      pSearchClose3Stats->plateNumber = plateNumber;
      pSearchClose3Stats->mosaicNumber = mosaicNumber;
      pSearchClose3Stats->binning = binning;
      pSearchClose3Stats->rotation = rotation;
    } else {
      if ((strcmp(pSearchClose3Stats->series,series) != 0)  ||
          (pSearchClose3Stats->plateNumber != plateNumber) ||
          (pSearchClose3Stats->mosaicNumber != mosaicNumber) ||
          (pSearchClose3Stats->binning != binning) ||
          (pSearchClose3Stats->rotation != rotation)) {
        printf("ERROR line %d for %s\n",__LINE__,copyLine);
        return(-1);
      }
    } 
    SetStepTime(pParseCommon,pSearchClose3Stats->seconds,0,1,lineCount,__LINE__);

    /* We are complete */
    pPlateStats->searchclose3Result = SEARCHCLOSE3_RESULT_COMPLETE;
    pPlateStats->chargeTime += pSearchClose3Stats->seconds;

  }




  
  return(0);
}

void PrintSearchClose3Stats(PPARSECOMMON pParseCommon) {
  double percentage;
  if (pParseCommon->searchclose3TotalPlates <= 0) {
    percentage = 0.0;
  } else {
    percentage = 100.0 * 0.00001 * ((100000 * pParseCommon->searchclose3SelectedPlates)/pParseCommon->searchclose3TotalPlates);
  }
  printf("searchclose3 Mask Plates %4d  SelectedPlates %4d  percentage %6.1f %% totalPlates %4d minSearchTableSize %4d maxSearchTableSize %4d\n",
         pParseCommon->searchclose3MaskPlates,
         pParseCommon->searchclose3SelectedPlates,
         percentage,
         pParseCommon->searchclose3TotalPlates,
         pParseCommon->searchclose3MinSearchTableSize,
         pParseCommon->searchclose3MaxSearchTableSize);
         

}

void UpdateSearchClose3Stats(PPARSECOMMON pParseCommon) {
	PPLATESTATS pPlateStats = &pParseCommon->plateStats;
	PPLATESTATS pPlateStats2 = NULL;
  PSEARCHCLOSE3STATS pSearchClose3Stats = &pPlateStats->searchclose3statstable;
  PMASKSECTOR pMaskSector = NULL;
  PPLATELISTENTRY pPlateListEntry;
  int plateListIndex;
  int maskSector;
  int cutoffTileOverlapCount; /* cutoff between single and dual overlap plates */
  int maxSearchTableSize = 0;
  int validMaskPlate = 0;
  int plateStatsIndex = 0;
  int sequence = SEQUENCE_NONE;
#if 0
  if ((strcmp(pPlateStats->series,"b") == 0) &&
      (pPlateStats->plateNumber == 8441)) {
    printf("line %4d at plate %s%05d\n",__LINE__,pPlateStats->series,pPlateStats->plateNumber);
  }
#endif


  if (pPlateStats->searchclose3Result != SEARCHCLOSE3_RESULT_COMPLETE) {
    printf("ERROR line %d UpdateSearchClose3Stats searchclose3Result %d for plate %s%05d_%02d\n",
           __LINE__,
           pPlateStats->searchclose3Result,
           pPlateStats->series,
           pPlateStats->plateNumber,
           pPlateStats->mosaicNumber);
    sequence = SEQUENCE_INCOMPLETE ;
  } else {
    if ((pSearchClose3Stats->multipleMaskSector > 0) ||
        (pSearchClose3Stats->multipleMaskFlag > 0) ||
        (pSearchClose3Stats->multipleMaskCenterFlag > 0)) {
      if (sequence == SEQUENCE_NONE) {
        sequence = SEQUENCE_MULTIPLE;
      }
    }
    /* Find the maximum search table size and be sure that it hasn't changed from an earlier plate */
    for (maskSector = 0; maskSector < pSearchClose3Stats->maxMaskSector; maskSector++) {
      pMaskSector = &pSearchClose3Stats->maskSectorTable[maskSector];
      if (pMaskSector->searchTableSize > maxSearchTableSize) {
        maxSearchTableSize = pMaskSector->searchTableSize;
      }
#ifdef COMBINE_CENTER_AREA
      if (pMaskSector->maxCenterMaskArea > pMaskSector->maxMaskArea) {
        pMaskSector->maxMaskArea = pMaskSector->maxCenterMaskArea;
      }
#endif /* COMBINE_CENTER_AREA */

    }
    if (pParseCommon->searchclose3MaxSearchTableSize == 0) {
      pParseCommon->searchclose3MaxSearchTableSize =  maxSearchTableSize;
      pParseCommon->searchclose3MinSearchTableSize =  maxSearchTableSize;
    } else { 
      if (pParseCommon->searchclose3MaxSearchTableSize != maxSearchTableSize) {
        printf("ERROR: line %4d maxSearchTableSize is %4d, expecting %4d for plate %s%05d_%02d\n",
               __LINE__,
               maxSearchTableSize,
               pParseCommon->searchclose3MaxSearchTableSize,
               pSearchClose3Stats->series,
               pSearchClose3Stats->plateNumber,
               pSearchClose3Stats->mosaicNumber);
      }
      if  (pParseCommon->searchclose3MaxSearchTableSize < maxSearchTableSize) {
        pParseCommon->searchclose3MaxSearchTableSize = maxSearchTableSize;
      }
    
    }


    for (maskSector = 0; maskSector < pSearchClose3Stats->maxMaskSector; maskSector++) {
      pMaskSector = &pSearchClose3Stats->maskSectorTable[maskSector];
      if (pMaskSector->tileOverlapCount == 0) {
        continue;
      }
      if (pSearchClose3Stats->minTileOverlapCount == 0) {
        pSearchClose3Stats->minTileOverlapCount = pMaskSector->tileOverlapCount;
        pSearchClose3Stats->maxTileOverlapCount = pMaskSector->tileOverlapCount;
      } else {
        if (pSearchClose3Stats->minTileOverlapCount > pMaskSector->tileOverlapCount) {
          pSearchClose3Stats->minTileOverlapCount = pMaskSector->tileOverlapCount;
        }
        if (pSearchClose3Stats->maxTileOverlapCount < pMaskSector->tileOverlapCount) {
          pSearchClose3Stats->maxTileOverlapCount = pMaskSector->tileOverlapCount;
        }
      }
    }
    cutoffTileOverlapCount = (pSearchClose3Stats->maxTileOverlapCount/2)+1; /* Decide whether this is a single or dual overlap (See search_close3.c line 2673) */

    for (maskSector = 0; maskSector < pSearchClose3Stats->maxMaskSector; maskSector++) {
      pMaskSector = &pSearchClose3Stats->maskSectorTable[maskSector];
      if (pMaskSector->tileOverlapCount == 0) {
        continue;
      } 
      if (pMaskSector->searchTableSize <  pParseCommon->searchclose3MaxSearchTableSize) {
        if (sequence == SEQUENCE_NONE) {
          sequence = SEQUENCE_INSUFFICIENT;
        }
#if 0
        printf("WARNING: line %4d undersized searchTableSize %4d for plate %s%05d_%02d\n",
               __LINE__,
               pMaskSector->searchTableSize,
               pSearchClose3Stats->series,
               pSearchClose3Stats->plateNumber,
               pSearchClose3Stats->mosaicNumber);
#endif
        if (pMaskSector->searchTableSize <  pParseCommon->searchclose3MinSearchTableSize) {
          pParseCommon->searchclose3MinSearchTableSize = pMaskSector->searchTableSize;
        }
        continue;
      }
      if ((pMaskSector->maxCenterMaskArea > 0) ||
          (pMaskSector->maxMaskArea > 0)) {
        if (sequence == SEQUENCE_NONE) {
          sequence = SEQUENCE_MASK;
        }
      }
      if (pMaskSector->tileOverlapCount > cutoffTileOverlapCount) {
        if (pMaskSector->totalMasksSector > 0) {
          validMaskPlate++; /* This plate has a mask in one of the overlap regions */
        } 
        if (pMaskSector->count_med1 > pSearchClose3Stats->dual_count_med1) {
          pSearchClose3Stats->dual_count_med1 = pMaskSector->count_med1;
        }
        if (pMaskSector->count_rms1 > pSearchClose3Stats->dual_count_rms1) {
          pSearchClose3Stats->dual_count_rms1 = pMaskSector->count_rms1;
        }
        if (pMaskSector->count_limit1 > pSearchClose3Stats->dual_count_limit1) {
          pSearchClose3Stats->dual_count_limit1 = pMaskSector->count_limit1;
        }
        if (pMaskSector->maxBinCount1 > pSearchClose3Stats->dual_maxBinCount1) {
          pSearchClose3Stats->dual_maxBinCount1 = pMaskSector->maxBinCount1;
        }
        if (pMaskSector->maxBinSNR1 > pSearchClose3Stats->dual_maxBinSNR1) {
          pSearchClose3Stats->dual_maxBinSNR1 = pMaskSector->maxBinSNR1;
        }
        if (pMaskSector->totalMasksSector > pSearchClose3Stats->dual_totalMasksSector) {
          pSearchClose3Stats->dual_totalMasksSector = pMaskSector->totalMasksSector;
        }
        if (pMaskSector->maxBinCount1 < pSearchClose3Stats->dual_maxBinCount1) {
          pSearchClose3Stats->dual_minBinCount1 = pMaskSector->maxBinCount1;
        }
        if (pMaskSector->count_limit2 > pSearchClose3Stats->dual_count_limit2) {
          pSearchClose3Stats->dual_count_limit2 = pMaskSector->count_limit2;
        }
        if (pMaskSector->maxMaskArea > pSearchClose3Stats->dual_maxMaskArea) {
          pSearchClose3Stats->dual_maxMaskArea = pMaskSector->maxMaskArea;
        }
        if (pMaskSector->maxCenterMaskArea > pSearchClose3Stats->dual_maxCtrMaskArea) {
          pSearchClose3Stats->dual_maxCtrMaskArea = pMaskSector->maxCenterMaskArea;
        }
      
        /* double overlay */
      } else {
        /* single overlay */
        if (pMaskSector->count_med1 > pSearchClose3Stats->single_count_med1) {
          pSearchClose3Stats->single_count_med1 = pMaskSector->count_med1;
        }
        if (pMaskSector->count_rms1 > pSearchClose3Stats->single_count_rms1) {
          pSearchClose3Stats->single_count_rms1 = pMaskSector->count_rms1;
        }
        if (pMaskSector->count_limit1 > pSearchClose3Stats->single_count_limit1) {
          pSearchClose3Stats->single_count_limit1 = pMaskSector->count_limit1;
        }
        if (pMaskSector->maxBinCount1 > pSearchClose3Stats->single_maxBinCount1) {
          pSearchClose3Stats->single_maxBinCount1 = pMaskSector->maxBinCount1;
        }
        if (pMaskSector->maxBinSNR1 > pSearchClose3Stats->single_maxBinSNR1) {
          pSearchClose3Stats->single_maxBinSNR1 = pMaskSector->maxBinSNR1;
        }
        if (pMaskSector->totalMasksSector > pSearchClose3Stats->single_totalMasksSector) {
          pSearchClose3Stats->single_totalMasksSector = pMaskSector->totalMasksSector;
        }
        if (pMaskSector->count_limit2 > pSearchClose3Stats->single_count_limit2) {
          pSearchClose3Stats->single_count_limit2 = pMaskSector->count_limit2;
        }
        if (pMaskSector->maxMaskArea > pSearchClose3Stats->single_maxMaskArea) {
          pSearchClose3Stats->single_maxMaskArea = pMaskSector->maxMaskArea;
        }
        if (pMaskSector->maxCenterMaskArea > pSearchClose3Stats->single_maxCtrMaskArea) {
          pSearchClose3Stats->single_maxCtrMaskArea = pMaskSector->maxCenterMaskArea;
        }
      }
    }




    printf("line %4d minTileOverlapCount %2d cutoffTileOverlapCount %2d maxTileOverlapCount %2d maxMaskSector %2d series %5s plateNumber %05d totalMasks %5d seconds %3d %3d %3d\n",
           __LINE__,
           pSearchClose3Stats->minTileOverlapCount,
           cutoffTileOverlapCount,
           pSearchClose3Stats->maxTileOverlapCount,
           pSearchClose3Stats->maxMaskSector,
           pPlateStats->series,
           pPlateStats->plateNumber,
           pSearchClose3Stats->totalMasks,
           pSearchClose3Stats->seconds,
           pPlateStats->chargeTime,
           pPlateStats->totalChargeTime);
        
    pParseCommon->searchclose3TotalPlates++;
    if (validMaskPlate > 0) {
      pParseCommon->searchclose3MaskPlates++;
      if (pSearchClose3Stats->dual_maxMaskArea > pSearchClose3Stats->single_maxMaskArea) {
        pParseCommon->searchclose3SelectedPlates++;
      } else {
        printf("WARNING: line %4d plate %s%05d_%02d is not selected\n",
               __LINE__,
               pSearchClose3Stats->series,
               pSearchClose3Stats->plateNumber,
               pSearchClose3Stats->mosaicNumber);
      }
    } else {
      printf("WARNING: line %4d plate %s%05d_%02d is not multiple\n",
             __LINE__,
             pSearchClose3Stats->series,
             pSearchClose3Stats->plateNumber,
             pSearchClose3Stats->mosaicNumber);

    }
    for (plateListIndex = 0; plateListIndex < pParseCommon->plateListEntryCount; plateListIndex++) {
      pPlateListEntry = &pParseCommon->pPlateListEntryTable[plateListIndex];
      if ((pPlateListEntry->plateNumber == pPlateStats->plateNumber) &&
          (strcmp(pPlateListEntry->series,pPlateStats->series) == 0)) {
        break;
      }
    }
    if (plateListIndex < pParseCommon->plateListEntryCount) {
      if (pPlateListEntry->okFlag) {
        sequence = SEQUENCE_OK;
      }
      if (pPlateListEntry->rescanFlag) {
        sequence = SEQUENCE_RESCAN;
      }
    }

  }
  pPlateStats->searchclose3Result = SEARCHCLOSE3_RESULT_RECORDED;         
  pSearchClose3Stats->sequence = sequence;
  for (plateStatsIndex = 0; plateStatsIndex < pParseCommon->plateStatsEntryCount; plateStatsIndex++) {
    pPlateStats2 = &pParseCommon->pPlateStatsTable[plateStatsIndex];
    if ((pPlateStats->plateNumber == pPlateStats2->plateNumber) &&
        (pPlateStats->mosaicNumber == pPlateStats2->mosaicNumber) &&
        (strcmp(pPlateStats->series,pPlateStats2->series) == 0)) {
      break;
    }
  }
  if (plateStatsIndex >= pParseCommon->plateStatsEntryCount) {
    if (plateStatsIndex > (SEARCHCLOSE3_MAX_PLATES-2)) {
      printf("ERROR: line %4d SEARCHCLOSE3_MAX_PLATE exceeded\n",__LINE__);
      exit(-1);
    }
    pPlateStats2 = &pParseCommon->pPlateStatsTable[plateStatsIndex];
    pParseCommon->plateStatsEntryCount++;
  }
  memcpy(pPlateStats2,pPlateStats,sizeof(PLATESTATS));
  

#if 0
  fprintf(pParseCommon->out_handle,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n",
          pPlateStats->series,
          pPlateStats->plateNumber,
          pPlateStats->mosaicNumber,
          pSearchClose3Stats->patternID,
          pSearchClose3Stats->minTileOverlapCount,
          pSearchClose3Stats->maxTileOverlapCount,
          pSearchClose3Stats->maxMaskSector,
          pSearchClose3Stats->totalMasks,    
          pSearchClose3Stats->maskCenterFlag,     /* 9 */
          pSearchClose3Stats->dual_count_med1,
          pSearchClose3Stats->dual_count_rms1,
          pSearchClose3Stats->dual_count_limit1,
          pSearchClose3Stats->dual_maxBinCount1,
          pSearchClose3Stats->dual_minBinCount1,
          pSearchClose3Stats->dual_maxBinSNR1,
          pSearchClose3Stats->dual_totalMasksSector,
          pSearchClose3Stats->dual_count_limit2,
          pSearchClose3Stats->dual_maxMaskArea,   
          pSearchClose3Stats->dual_maxCtrMaskArea,   /* 19 */
          pSearchClose3Stats->single_count_med1,
          pSearchClose3Stats->single_count_rms1,
          pSearchClose3Stats->single_count_limit1,
          pSearchClose3Stats->single_maxBinCount1,
          pSearchClose3Stats->single_maxBinSNR1,
          pSearchClose3Stats->single_totalMasksSector,
          pSearchClose3Stats->single_count_limit2,
          pSearchClose3Stats->single_maxMaskArea, 
          pSearchClose3Stats->single_maxCtrMaskArea,
          sequence);  /* 29 */
#endif

  return;
         
}
void ReadSearchClose3File(PPARSECOMMON pParseCommon,char *searchclose3_name) {
  FILE *searchclose3_handle = NULL;
  int lineCount = 0;
  int lineLen;
  char * inBuffer;
  char copyLine[MAX_BUFFER];
  char inLine[MAX_BUFFER];
  PPLATELISTENTRY pPlateListEntry = NULL;
  int okCount = 0;
  int errorCount = 0;
  int rescanCount = 0;

  if (searchclose3_name[0] == 0) {
    printf("line %4d ERROR: searchclose3_name is not specified\n");
    return;
  }
  searchclose3_handle = fopen(searchclose3_name,"rt");
  if (searchclose3_handle == NULL) {
    printf("ERROR: line %4d Could not open list file %s\n",__LINE__,searchclose3_name);
    return;
  }
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,searchclose3_handle);
    if (inBuffer == NULL) {
      break;
			
    }
    lineLen = strlen(inBuffer);
    if (lineLen > pParseCommon->maxLineLen) {
      pParseCommon->maxLineLen = lineLen;
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
    lineCount++;
  }
  if (searchclose3_handle != NULL) {
    fclose(searchclose3_handle);
  }
  pParseCommon->plateListEntryAlloc = lineCount;
  pParseCommon->pPlateListEntryTable = (PPLATELISTENTRY)calloc(pParseCommon->plateListEntryAlloc,sizeof(PLATELISTENTRY));
  if (pParseCommon->pPlateListEntryTable == NULL) {
    printf("ERROR: line %4d failed to allocate pPlateListEntryTable\n",__LINE__);
    return;
  }

  searchclose3_handle = fopen(searchclose3_name,"rt");
  if (searchclose3_handle == NULL) {
    printf("ERROR: line %4d Could not open list file %s\n",__LINE__,searchclose3_name);
    return;
  }
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,searchclose3_handle);
    if (inBuffer == NULL) {
      break;
			
    }
    lineLen = strlen(inBuffer);
    if (lineLen > pParseCommon->maxLineLen) {
      pParseCommon->maxLineLen = lineLen;
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
    strcpy(copyLine,inBuffer);
    pPlateListEntry = &pParseCommon->pPlateListEntryTable[pParseCommon->plateListEntryCount];
    memset(pPlateListEntry,0,sizeof(PLATELISTENTRY));
    if (ParseFilename2(inBuffer,pPlateListEntry->series,&pPlateListEntry->plateNumber) == 0) {
      printf("line %4d ERROR: ReadSearchClose3 could not parse filename in %s\n",copyLine);
      continue;
    }
    if (strstr(copyLine,"OK") != NULL) {
      okCount++;
      pPlateListEntry->okFlag = 1;
    }
    if (strstr(copyLine,"ERROR") != NULL) {
      errorCount++;
      pPlateListEntry->errorFlag = 1;
    }
    if (strstr(copyLine,"RESCAN") != NULL) {
      rescanCount++;
      pPlateListEntry->rescanFlag = 1;
    }


    pParseCommon->plateListEntryCount++;
  }
  if (searchclose3_handle != NULL) {
    fclose(searchclose3_handle);
  }


  printf("line %4d read %4d of %4d lines okCount %d rescanCount %d errorCount %d in %s\n",__LINE__,pParseCommon->plateListEntryCount,pParseCommon->plateListEntryAlloc,okCount,rescanCount,errorCount,searchclose3_name);
  return;

}
void ProcessPlateStats(PPARSECOMMON pParseCommon, int lineCount, char* log_name) {
  PPLATESTATS pPlateStats = &pParseCommon->plateStats;
  PSTEPSTATS pStepStats;
  PSTEPSTATS pStepStats2;
  PSTEPTABLE pStepTable;
  PSTEPTABLE pStepTable2;
  PSTEPSTATS pTotalStepStats;
  PSTEPSTATS pUnknownStepStats;
  PSTEPSTATS pXFER1Stats = NULL;
  PSTEPSTATS pAstrometryStats = NULL;

  int stepCounter;
  time_t totalStepTime = 0;
  time_t totalTimestampTime = 0;
  int nonChargeStepBitmap = 0;
  int intSeconds;

  if (pParseCommon->searchclose3Active != 0) {
    pParseCommon->searchclose3Active = 0;
    UpdateSearchClose3Stats(pParseCommon);
  }

  if (pPlateStats->series[0] != 0) {
    for (stepCounter = 0; stepCounter < MAX_STEP; stepCounter++) {
      pStepStats = &pPlateStats->stepStats[stepCounter];
      totalStepTime += pStepStats->stepTime;
    }
		
    pStepTable = &stepTable[pPlateStats->curStep];
    pStepStats = &pPlateStats->stepStats[pPlateStats->curStep];
		

    if (pPlateStats->stepsSeen == ((1 << STEP_XFER1) | ( 1 << STEP_SCAMPASTROMETRY))) { 
      pXFER1Stats = &pPlateStats->stepStats[STEP_XFER1];
      pAstrometryStats = &pPlateStats->stepStats[STEP_SCAMPASTROMETRY];
      pXFER1Stats->stepTime = pPlateStats->chargeTime+pPlateStats->totalChargeTime;
      pAstrometryStats->stepTime = pPlateStats->totalTime - pXFER1Stats->stepTime;
      totalStepTime = pPlateStats->totalTime;
    } else if ((pPlateStats->singleFlag) ||
               (pPlateStats->stepsSeen == (1 << STEP_XFER1))) {
      pXFER1Stats = &pPlateStats->stepStats[STEP_XFER1];
      totalStepTime -= pStepStats->stepTime;
      pStepStats->stepTime = 0;
      totalStepTime -= pXFER1Stats->stepTime;
      pXFER1Stats->stepTime = 0;
      SetStepTime(pParseCommon,pPlateStats->totalTime,0,1,lineCount,__LINE__);
      totalStepTime += pPlateStats->totalTime;
    } else {
      if ((pPlateStats->beginTime != 0) &&
          (pPlateStats->endTime != 0)) {
        totalTimestampTime = pPlateStats->endTime - pPlateStats->beginTime;
        if (pPlateStats->curStep == STEP_XFER2) {
          totalTimestampTime += pStepStats->stepTime;
        }
      }
    }
#if 0
    printf("steps seen 0x%x ~CHARGE_STEP_MASK 0x%x conditional %d Total step time %lld total timestamptime %lld\n",
           pPlateStats->stepsSeen,~CHARGE_STEP_MASK,((pPlateStats->stepsSeen | (~CHARGE_STEP_MASK)) == 0),totalStepTime,totalTimestampTime);
#endif

    if (totalTimestampTime < totalStepTime) {
      nonChargeStepBitmap = pPlateStats->stepsSeen & (~CHARGE_STEP_MASK);
      if (nonChargeStepBitmap == 0) {
        /* Here we have had no steps using the entire interval between two "date" markers. -- charge it to STEP_UNKNOWN */
        intSeconds =  pPlateStats->totalTime - totalStepTime;
        pUnknownStepStats = &pPlateStats->stepStats[STEP_UNKNOWN];
        pUnknownStepStats->stepTime += intSeconds;
        pUnknownStepStats->stepCount++;
        totalStepTime += intSeconds;
        totalTimestampTime = totalStepTime;
      } else {
        for (stepCounter = 0; stepCounter < MAX_STEP; stepCounter++) {
          pStepTable2 = &stepTable[stepCounter];
          if (nonChargeStepBitmap == (1 << pStepTable2->stepmask)) {
            /* There is only one other step active than those for which we already have the time, so charge all excess time to this step */
            pStepStats2 = &pPlateStats->stepStats[stepCounter];
            intSeconds = pPlateStats->totalTime - totalStepTime;
            pStepStats2->stepTime += intSeconds;
            pStepStats2->stepCount++;
            totalStepTime += intSeconds;
            totalTimestampTime = totalStepTime;
					 
            break;
          }
        }

      }
    }

    for (stepCounter = 0; stepCounter < MAX_STEP; stepCounter++) {
      pStepStats = &pPlateStats->stepStats[stepCounter];
      pTotalStepStats = &pParseCommon->totalStepStats[stepCounter];
      pStepTable = &stepTable[stepCounter];
      if (pStepStats->stepCount > 0) {
        if (pParseCommon->verbose) {
          printf("Step %2d  %22s count %d time %d\n",stepCounter,pStepTable->stepname,pStepStats->stepCount,pStepStats->stepTime);
        }
        if ((stepCounter == STEP_MOSAICS) && (pStepStats->stepTime > 0)) {
          fprintf(pParseCommon->out_handle,"%d\n",pStepStats->stepTime);
        }

        pTotalStepStats->stepTime += pStepStats->stepTime;
        if (pStepStats->startedFlag > 0) {
          pTotalStepStats->startedFlag++;
        }
        if (pStepStats->stepCount >= 0) {
          pTotalStepStats->stepCount++;
          if ((pParseCommon->reprocessFlag != 0)  &&
              (((1 << stepCounter) & REPROCESS_MASK) == 0)) {
            printf("WARNING: step %s executed while reprocess flag is set in file %s\n",pStepTable->stepname,log_name);
          }
        }
      }
    }
    pParseCommon->totalTime += pPlateStats->totalTime;
    pParseCommon->totalStepTime += totalStepTime;
    pParseCommon->totalChargeTime += pPlateStats->totalChargeTime;
    pParseCommon->totalTimestampTime += totalTimestampTime;
    if ((pParseCommon->astrometrynetActive != 0) &&
        (pPlateStats->astrometrynetResult != 0) &&
        (pPlateStats->astrometrynetSolveStep != 0)) {
      if (pPlateStats->astrometrynetResult > 0) {
        /* Success */
        UpdateAstrometryNetStats(pParseCommon,pParseCommon->allSuccessAstrometrynetTable,pParseCommon->allSuccessAstrometrynetSeconds);
        if (strcmp(pPlateStats->series,ASTROMETRYNETSERIES) == 0) {
          UpdateAstrometryNetStats(pParseCommon,pParseCommon->seriesSuccessAstrometrynetTable,pParseCommon->seriesSuccessAstrometrynetSeconds);              
        }
      } else {
        /* Failed */
        UpdateAstrometryNetStats(pParseCommon,pParseCommon->allFailedAstrometrynetTable,pParseCommon->allFailedAstrometrynetSeconds);
        if (strcmp(pPlateStats->series,ASTROMETRYNETSERIES) == 0) {
          UpdateAstrometryNetStats(pParseCommon,pParseCommon->seriesFailedAstrometrynetTable,pParseCommon->seriesFailedAstrometrynetSeconds);              
        }
      }
    }


    if ((totalTimestampTime != 0) && 
        (totalTimestampTime != totalStepTime)) {
      printf("ERROR: totalStepTime not totalTimestampTime. Finished plate %5s%05d totalTimestampTime %5lld totalChargeTime %5lld totalStepTime %5lld totalTime %5lld line %d\n",pPlateStats->series,pPlateStats->plateNumber,totalTimestampTime,pPlateStats->totalChargeTime,totalStepTime,pPlateStats->totalTime,lineCount);
    } else if (totalStepTime > pPlateStats->totalTime) {
      printf("ERROR: totalStepTime > totalTime Finished plate %5s%05d totalTimestampTime %5lld totalChargeTime %5lld totalStepTime %5lld totalTime %5lld line %d\n",pPlateStats->series,pPlateStats->plateNumber,totalTimestampTime,pPlateStats->totalChargeTime,totalStepTime,pPlateStats->totalTime,lineCount);
    } else {
      if (pParseCommon->verbose) {
        printf("Finished plate %5s%05d totalTimestampTime %5lld totalChargeTime %5lld totalStepTime %5lld totalTime %5lld line %d\n",pPlateStats->series,pPlateStats->plateNumber,totalTimestampTime,pPlateStats->totalChargeTime,totalStepTime,pPlateStats->totalTime,lineCount);
      }
    }
  }

	

  memset(pPlateStats,0,sizeof(PLATESTATS));

}


int ProcessLogBuffer(PPARSECOMMON pParseCommon, char* log_name) {
  char *inBuffer;
  char copyLine[MAX_BUFFER];
  char inLine[MAX_BUFFER];
  char parseLine[MAX_BUFFER];
  char *parsePtr1;
  char *parsePtr2;
  size_t tmpBytesCopied;
  int parseCount;
  int solutionIndex;
  char token[MAX_BUFFER];
  char startToken[MAX_BUFFER];
  int nvals;
  FILE *log_handle = NULL;
  int lineLen;
  int lineCount = 0;
  int recStatCount = 0;
  int parseState = 1; /* 1 = idle; 2 = in progress ; 3 = flush */
  char *charPtr;
  char *charPtr2;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int processCount = 0;
  time_t lastTime = 0;
  time_t prevTime = 0;
  PSTEPSTATS pStepStats = NULL;
  PSTEPTABLE pStepTable;
  PSTEPTABLE pStepTableNext;
  PSTEPTABLE pStepTablePrev;
  PPLATESTATS pPlateStats = &pParseCommon->plateStats;
  double realSeconds;
  int intSeconds;
  int forwardFlag;
  int saveStep;
  int tmpSolutionNumber;
  int catalogNumber;
  PSIZESTATS pSizeStats;

  printf("Processing %s\n",log_name);
  log_handle = fopen(log_name,"rt");
  startToken[0] = 0;
  if (log_handle == NULL) {
    printf("ERROR: line %4d Could not open list file %s\n",__LINE__,log_name);
    return(-1);
  }
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,log_handle);
    if (inBuffer == NULL) {
      break;
			
    }
    lineLen = strlen(inBuffer);
    if (lineLen > pParseCommon->maxLineLen) {
      pParseCommon->maxLineLen = lineLen;
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
    lineCount++;
#if 0
    if (lineCount == 138) {
      printf("At line %d %s\n",lineCount,inLine);
    }
#endif
    if (strstr(inBuffer,"CHILD") != 0) {
      strcpy(copyLine,inBuffer);
      charPtr = inBuffer+strlen("CHILD:");
      while ((*charPtr == ' ') && (*charPtr != 0)) {
        charPtr++;
      }
      charPtr2 = charPtr;
      while ((*charPtr2 != ' ') && (*charPtr2 != 0)) {
        charPtr2++;
      }
      *charPtr2 = 0;
      strcpy(token,charPtr);
      if (parseState == 1) {
        if ((strstr(copyLine,"Exited") == NULL) || (strstr(copyLine," index ") == NULL)) {
          if (strstr(copyLine,"Starting") == NULL) {
            printf("ERROR: 'Starting' expected in line %d %s of %s\n",lineCount,copyLine,log_name);
            fclose(log_handle);
            return(-1);
          }								 
          memset(pPlateStats,0,sizeof(PLATESTATS));
          strcpy(startToken,token);
          if (ParseFilename2(startToken,pPlateStats->series,&pPlateStats->plateNumber) == 0) {
            if ((strstr(copyLine,"FitsLinearity") != 0) ||
                (strstr(copyLine,"FitsMedian") != 0)) {
              parseState = 3; /* Flush rest */
            } else {
              printf("ERROR: Failed to parse startToken %s in line %d %s of %s\n",startToken,lineCount,copyLine,log_name);
              fclose(log_handle);
              return(-1);
            }
          }
          if (pParseCommon->verbose) {
            printf("Have token %s in line %d of %s\n",startToken,lineCount,log_name);
          }
         
#ifdef DEBUG_SINGLE
          if (strstr(startToken,DEBUG_SINGLE) != NULL) {
            printf("%s\n",copyLine);
          }
#endif /* DEBUG_SINGLE */
          parseState = 2;
          processCount++;
          if ((pParseCommon->astrometrynetFlag != 0) && 
              (strstr(copyLine,"find_astrometry.csh") != NULL)) {
            SetNewStep(pParseCommon,STEP_ASTROMETRYNET,lineCount,__LINE__);
            pParseCommon->astrometrynetActive = 1; /* We are in the correct phase */
          }
        }
      } else {
#ifdef DEBUG_SINGLE
        if (strstr(startToken,DEBUG_SINGLE) != NULL) {
          printf("%s\n",copyLine);
        }
#endif /* DEBUG_SINGLE */
        if (strstr(copyLine,"Exited") == NULL) {
          printf("ERROR: 'Exited' expected in line %d %s of %s\n",lineCount,copyLine,log_name);
          fclose(log_handle);
          return(-1);
        }								 
        if (strcmp(token,startToken) != 0) {
          if (strstr(copyLine," index ") == NULL) {
            printf("ERROR: token mismatch %s %s in line %d %s of %s\n",token,startToken,lineCount,copyLine,log_name);
            fclose(log_handle);
            return(-1);
          }
        } else {
          charPtr = strstr(copyLine," at ");
          if (charPtr != NULL) {
            charPtr += strlen(" at ");
            nvals = sscanf(charPtr,"%d",&pPlateStats->totalTime);
            if (nvals != 1) {
              printf("ERROR: failed to decode time in line %d %s of %s\n",lineCount,copyLine,log_name);
              return(-1);
            }
          } else {
            printf("ERROR: failed to find time in line %d %s of %s\n",lineCount,copyLine,log_name);
            return(-1);
          }
						

          if (parseState == 2) {
            parseState = 1;
            startToken[0] = 0;
            ProcessPlateStats(pParseCommon,lineCount,log_name);
          } else if (parseState == 3) {
            parseState = 1;
            ProcessPlateStats(pParseCommon,lineCount,log_name);
            startToken[0] = 0;
          } else {
            printf("ERROR: unrecognized parse state %d in line %d %s of %s\n",parseState,lineCount,copyLine,log_name);
            fclose(log_handle);
            return(-1);
          }
        }

      }
				

    } else {
#ifdef DEBUG_SINGLE
      if (strstr(startToken,DEBUG_SINGLE) != NULL) {
        printf("%s\n",inLine);
      } else {
        if (strstr(inLine,"Recstat") == inLine) {
          recStatCount++;
        }
      }
#else /* DEBUG_SINGLE */
      if (strstr(inLine,"Recstat") == inLine) {
        recStatCount++;
      }
      if (pParseCommon->astrometrynetFlag != 0) {
        if ((pParseCommon->astrometrynetActive != 0) &&
            (strstr(inLine,"Phase") == inLine)) {
          /* We are done with this log */
#if 0
          printf("code line %d exiting astrometrynetActive phase\n",__LINE__);
#endif 
          break;
        } else {
          if ((pParseCommon->astrometrynetActive != 0) && (parseState == 2)) {
            if (strstr(inLine,"db2xyls") == inLine) {
              if (ParseAstrometryNetParameters(pParseCommon,inLine) != 0) {
                parseState = 3;
                continue;
              }
            } else if (strstr(inLine,"solve-field") == inLine) {
              if (ParseAstrometryNetParameters(pParseCommon,inLine) != 0) {
                parseState = 3;
                continue;
              }
            } else if (strstr(inLine,"No result for ") == inLine) {
              pPlateStats->astrometrynetResult = -1;
            } else if (strstr(inLine,"Solved: ") == inLine) {
              pPlateStats->astrometrynetResult = 1;
            } 
          }
          continue;
        }
      }
      if (ParseTime(pParseCommon,inLine,&lastTime,lineCount) == 0) {
        pPlateStats->endTime = lastTime;
        pStepStats = &pPlateStats->stepStats[pPlateStats->curStep];
        pStepTable = &stepTable[pPlateStats->curStep];
        /* Have a valid time */
        if (pPlateStats->beginTime == 0) {
          pPlateStats->beginTime = lastTime;
        } else {
          if (pStepTable->dateFlag != 0) {
            SetStepTime(pParseCommon,lastTime-prevTime-pPlateStats->chargeTime,pPlateStats->chargeTime,0,lineCount,__LINE__);
            pPlateStats->chargeTime = 0;
          } else if ((pPlateStats->curStep == STEP_XFER1) ||
                     ((pPlateStats->curStep > STEP_MOSAICS) &&
                      (pPlateStats->prevStep > STEP_MOSAICS) &&
                      (pPlateStats->curStep < STEP_XFER2))) {
            pStepTablePrev = &stepTable[pPlateStats->prevStep];
            pStepTableNext = &stepTable[pPlateStats->curStep+1];
            if ((pStepTablePrev->dateFlag == 2) && pStepTableNext->dateFlag > 0) {
              SetStepTime(pParseCommon,lastTime-prevTime-pPlateStats->chargeTime,pPlateStats->chargeTime,0,lineCount,__LINE__);
              pPlateStats->chargeTime = 0;
            } else {
              if (pPlateStats->curStep == STEP_SEARCHCLOSE) {
                /* search_close registers only on success, so set the step to the previous step */
                SetNewStep(pParseCommon,pPlateStats->prevStep,lineCount,__LINE__);
                SetStepTime(pParseCommon,lastTime-prevTime-pPlateStats->chargeTime,pPlateStats->chargeTime,0,lineCount,__LINE__);
                pPlateStats->chargeTime = 0;
								
              } else if ((pPlateStats->curStep == STEP_XFER1) ||
                         (pPlateStats->curStep == STEP_SCAMPASTROMETRY) ||
                         (pPlateStats->curStep == STEP_ASTROMETRYSECOND)) {
                SetStepTime(pParseCommon,lastTime-prevTime-pPlateStats->chargeTime,pPlateStats->chargeTime,0,lineCount,__LINE__);
                pPlateStats->chargeTime = 0;
              } else {
                if (pStepTable->singleFlag == 0) {
                  printf("ERROR: valid ParseTime result discarded for %s %d seconds for line %d parselogs.c line %d\n",pStepTable->stepname,lastTime-prevTime,lineCount,__LINE__);
                }
              }
            }
          } else {
            if (pStepTable->singleFlag == 0) {
              if ((pPlateStats->curStep == STEP_UNKNOWN) && 
                  (pPlateStats->stepsSeen <= 1)) {
                /* Here we must charge this time to the unknown step */
                pStepStats->stepTime += lastTime-prevTime-pPlateStats->chargeTime;
                pStepStats->stepCount++;
              } else {
                if (pPlateStats->prevStep != STEP_UNKNOWN) {
                  SetNewStep(pParseCommon,pPlateStats->prevStep,lineCount,__LINE__);
                  SetStepTime(pParseCommon,lastTime-prevTime-pPlateStats->chargeTime,pPlateStats->chargeTime,0,lineCount,__LINE__);
                  pPlateStats->chargeTime = 0;
                } else {
                  printf("ERROR: valid ParseTime result discarded for %s %d seconds for line %d parselogs.c line %d\n",pStepTable->stepname,lastTime-prevTime,lineCount,__LINE__);
                }
              }
            }
          }
        }
        charPtr = strstr(inLine,"Run SCAMP to correct distortions");
        if (charPtr != NULL) {
          if (pPlateStats->curStep == STEP_SCAMPSECOND) {
            SetNewStep(pParseCommon,STEP_SCAMPTHIRD,lineCount,__LINE__);
          } else {
            SetNewStep(pParseCommon,STEP_SCAMPFIRST,lineCount,__LINE__);
          }
        }
        charPtr = strstr(inLine,"Run SCAMP 01 to correct distortions");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_SCAMPSECOND,lineCount,__LINE__);
        }
        prevTime = lastTime;

      } else {
        if (pParseCommon->searchclose3Flag != 0) {
          if ((pParseCommon->searchclose3Active != 0) && (parseState == 2)) {
            if (ParseSearchClose3Parameters(pParseCommon,lineCount,inLine) != 0) {
              parseState = 3;
              continue;
            }
            continue;
          }
        }
        charPtr = strstr(inLine,"FitsMosaic");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_MOSAICS,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"solve-field");
        if ((charPtr != NULL) && (pPlateStats->curStep != STEP_SCAMPASTROMETRY)) {
          charPtr = strstr(inLine,"_full");
          if (charPtr != NULL) {
            SetNewStep(pParseCommon,STEP_SCAMPASTROMETRY,lineCount,__LINE__);
          } else {

            charPtr = strstr(inLine,"_s");
            if (charPtr != NULL) {
              nvals = sscanf(charPtr,"_s%d",&tmpSolutionNumber);
              if (nvals == 1) {
                SetNewStep(pParseCommon,STEP_ASTROMETRYSECOND,lineCount,__LINE__);
              } else {
                SetNewStep(pParseCommon,STEP_ASTROMETRYNET,lineCount,__LINE__);
              }
            } else {
              SetNewStep(pParseCommon,STEP_ASTROMETRYNET,lineCount,__LINE__);
            }
          }
        }
        charPtr = strstr(inLine,"_none.db not found");
        if (charPtr != NULL) {
          charPtr = strstr(inLine,"WARNING");
          if (charPtr != NULL) {
            SetNewStep(pParseCommon,STEP_ASTROMETRYSECOND,lineCount,__LINE__);
          }
        }

        charPtr = strstr(inLine,"AstrometryWCS");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_IMWCS,lineCount,__LINE__);
        }
			
        charPtr = strstr(inLine,"setScampRMS");
        if (charPtr != NULL) {
          charPtr = strstr(inLine,"setfullfail");
          if (charPtr != NULL) {
            SetNewStep(pParseCommon,STEP_SCAMPASTROMETRY,lineCount,__LINE__);
          }

        }
	
        charPtr = strstr(inLine,"Found (");
        if (charPtr == inLine) {
          SetNewStep(pParseCommon,STEP_PREPAREMOSAIC,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"Copy location");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_PREPAREMOSAIC,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"Executing: sex");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_FULLSEXTRACTOR,lineCount,__LINE__);
        }

        /* update_sextractor is in:                          first pass                                             second pass
           run_sextractor.csh            STEP_FULLSEXTRACTOR  (catalog2.tmp) (ww.db)
           run_sextractor.csh            STEP_SEXTRACTOR1      (tnx.db) (ww.fit)                                      (_none.db)
           run_scamp.csh                 STEP_SCAMPFIRST      (.tmp2) (ww.fit)
           run_sextractor_second.csh     STEP_SEXTRACTOR2     (ww_tnx.fit)                                            ww_s1_tnx.hdr
           run_scamp.csh 2               STEP_SCAMPTHIRD      (.tmp2, ww_tnx.fit.provisional and .filtered.head)
           run_sextractor_second.csh 2   STEP_SEXTRACTOR2B    (ww_tnx.fit)                                            ww_s1_tnx.hdr
        */
        charPtr = strstr(inLine,"/update_sextractor ");
        if (charPtr != NULL) {
          charPtr += strlen("/update_sextractor ");
          charPtr = strstr(charPtr,"_s");
          if (charPtr != NULL) {
            nvals = sscanf(charPtr,"_s%d",&tmpSolutionNumber);
            if (nvals == 1) {
              if ((tmpSolutionNumber != 0) && 
                  (pPlateStats->solutionNumber != tmpSolutionNumber)) {
                pPlateStats->solutionNumber = tmpSolutionNumber;
                printf("Setting solution number to %d in line %d parselogs.c line %d\n",pPlateStats->solutionNumber,lineCount,__LINE__);
              }
            }
          }
          if (pPlateStats->solutionNumber == 0) {
            charPtr = strstr(inLine,".tmp2");
            if (charPtr != NULL) {
              charPtr = strstr(inLine,".provisional");
              if (charPtr == NULL) {
                pPlateStats->updateSextractorStep = STEP_SCAMPFIRST;
              } else {
                pPlateStats->updateSextractorStep = STEP_SCAMPTHIRD;
              }
            } else {
              charPtr = strstr(inLine,"catalog2.tmp");
              if (charPtr != NULL) {
                if ((pPlateStats->stepsSeen & (1 << STEP_FULLSEXTRACTOR)) != 0) {
                  pPlateStats->updateSextractorStep = STEP_FULLSEXTRACTOR;
                } else {
                  printf("ERROR: unrecognized update_sextractor  in line %8d parsecommon.c line %d\n",lineCount,__LINE__);
                  pPlateStats->updateSextractorStep = STEP_UNKNOWN;
                }
              } else {
                charPtr = strstr(inLine,"ww_tnx.fit");
                if (charPtr != NULL) {
                  charPtr += 10;
                  if ((charPtr[0] == 0) || (charPtr[0] == ' ')) {
                    if ((pPlateStats->stepsSeen & (1 << STEP_SCAMPTHIRD)) == 0) {
                      pPlateStats->updateSextractorStep = STEP_SEXTRACTOR2;
                    } else {
                      pPlateStats->updateSextractorStep = STEP_SEXTRACTOR2B;
                    }
                  } else {
                    printf("ERROR: unrecognized update_sextractor  in line %8d parsecommon.c line %d\n",lineCount,__LINE__);
                    pPlateStats->updateSextractorStep = STEP_UNKNOWN;
                  }
                } else {
                  charPtr = strstr(inLine,"ww.fit");
                  if (charPtr != NULL) {
                    pPlateStats->updateSextractorStep = STEP_SEXTRACTOR1;
                  } else {
                    printf("ERROR: unrecognized update_sextractor  in line %8d parsecommon.c line %d\n",lineCount,__LINE__);
                    pPlateStats->updateSextractorStep = STEP_UNKNOWN;									
                  }
                }
              }
            }

          } else {
            /* For solutionNumber > 0 */
            charPtr = strstr(inLine,"_none.db");
            if (charPtr != NULL) {
              pPlateStats->updateSextractorStep = STEP_SEXTRACTOR1;
            } else {
              charPtr = strstr(inLine,"_tnx.hdr");
              if (charPtr != NULL) {
                if ((pPlateStats->stepsSeen & (1 << STEP_SCAMPTHIRD)) == 0) {
                  pPlateStats->updateSextractorStep = STEP_SEXTRACTOR2;
                } else {
                  pPlateStats->updateSextractorStep = STEP_SEXTRACTOR2B;
                }
              } else {
                printf("ERROR: unrecognized update_sextractor  in line %8d parsecommon.c line %d\n",lineCount,__LINE__);
                pPlateStats->updateSextractorStep = STEP_UNKNOWN;
              }
            }
          }
          DisplayUpdateSextractorStep(pParseCommon,lineCount,__LINE__);
        }
       
        charPtr = strstr(inLine,"search_close3");
        if (charPtr != NULL) {
          if (pParseCommon->searchclose3Flag != 0) {
            if (pParseCommon->searchclose3Active == 0) {
              SetNewStep(pParseCommon,STEP_SEARCHCLOSE3,lineCount,__LINE__);
              pParseCommon->searchclose3Active = 1;
            }
          }
        } else {
          charPtr = strstr(inLine,"search_close");
          if (charPtr != NULL) {
            SetNewStep(pParseCommon,STEP_SEARCHCLOSE,lineCount,__LINE__);
          }
        }
			 
        charPtr = strstr(inLine,"copypipeline exiting");
        if ((charPtr != NULL) && (pPlateStats->curStep == STEP_UNKNOWN)) {
          SetNewStep(pParseCommon,STEP_XFER1,lineCount,__LINE__);
        }
			 
        charPtr = strstr(inLine,"filter_wedge");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_WEDGE,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"matchstars");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_MATCH2,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"filter_defect");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_DEFECT,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"prepare_octave");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_PREPAREOCTAVE,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"divide_annul9");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_DIVIDEANNULAR,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"colorterm");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_COLORTERM,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"annular9");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_ANNULAR,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"spatial bins found");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_INGESTOCTAVE,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"contour_plot");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_LOCALCALIBRATION,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"recover_points");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_RECOVERPOINTS,lineCount,__LINE__);
        }
				
        charPtr = strstr(inLine,"filter_multiple");
        if (charPtr != NULL) {
          SetNewStep(pParseCommon,STEP_FILTERMULTIPLE,lineCount,__LINE__);
        }
        charPtr = strstr(inLine,"maxMemory");
        if (charPtr != NULL) {
          ProcessMaxMemory(pParseCommon,inLine,log_name);
        }
				
				
        if ((pPlateStats->curStep != STEP_FILTERMULTIPLE) && 
            (pPlateStats->curStep != STEP_ASTROMETRYSECOND)) {
          charPtr = strstr(inLine,"_s");
          if (charPtr != NULL) {
            nvals = sscanf(charPtr,"_s%d",&tmpSolutionNumber);
            if (nvals == 1) {
              if ((tmpSolutionNumber != 0) && 
                  (pPlateStats->solutionNumber != tmpSolutionNumber)) {
                pPlateStats->solutionNumber = tmpSolutionNumber;
                printf("Setting solution number to %d in line %d parselogs.c line %d\n",pPlateStats->solutionNumber,lineCount,__LINE__);
              }
            }
          }
        }
        if (pPlateStats->catalogNumber == 0) {
          for (catalogNumber = 1; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
            charPtr = strstr(inLine,catalogText[catalogNumber]);
            if (charPtr != NULL) {
              pPlateStats->catalogNumber = catalogNumber;
              break;
            }
          }
        }
        pStepStats = &pPlateStats->stepStats[pPlateStats->curStep];
        charPtr = strstr(inLine,"seconds");
        if (charPtr != NULL) {
          if ((pPlateStats->curStep == STEP_SEARCHCLOSE) &
              ((strstr(inLine,"Maximum SNR") != NULL) ||
               (strstr(inLine,"Maximum value") != NULL))) {
            nvals = sscanf(charPtr,"seconds %d\n",&intSeconds);
            if (nvals == 1) {
              saveStep = pPlateStats->prevStep;
              SetStepTime(pParseCommon,intSeconds,0,0,lineCount,__LINE__);
              pPlateStats->chargeTime += intSeconds;
              SetNewStep(pParseCommon,saveStep,lineCount,__LINE__);
            }
          }
          charPtr2 = strstr(inLine,"update_sextractor starsin");
          if (charPtr2 != NULL) {
            nvals = sscanf(charPtr,"seconds %d",&intSeconds);
            if (nvals == 1) {
              DisplayUpdateSextractorStep(pParseCommon,lineCount,__LINE__);
              switch(pPlateStats->updateSextractorStep) {
              case STEP_FULLSEXTRACTOR:
                /* Ignore this time */
                break;
              case STEP_SEXTRACTOR1:
                saveStep = pPlateStats->curStep;
                SetNewStep(pParseCommon,STEP_SEXTRACTOR1,lineCount,__LINE__);
                SetStepTime(pParseCommon,intSeconds,0,0,lineCount,__LINE__);
                pPlateStats->chargeTime += intSeconds;
                SetNewStep(pParseCommon,saveStep,lineCount,__LINE__);							
                break;
              case STEP_SCAMPFIRST:
                /* Ignore this time */
                break;
              case STEP_SEXTRACTOR2:
                saveStep = pPlateStats->curStep;
                SetNewStep(pParseCommon,STEP_SEXTRACTOR2,lineCount,__LINE__);
                SetStepTime(pParseCommon,intSeconds,0,0,lineCount,__LINE__);
                pPlateStats->chargeTime += intSeconds;
                SetNewStep(pParseCommon,saveStep,lineCount,__LINE__);							
                break;
              case STEP_SCAMPTHIRD:
                /* Ignore this time */
                break;
              case STEP_SEXTRACTOR2B:
                saveStep = pPlateStats->curStep;
                SetNewStep(pParseCommon,STEP_SEXTRACTOR2B,lineCount,__LINE__);
                SetStepTime(pParseCommon,intSeconds,0,0,lineCount,__LINE__);
                pPlateStats->chargeTime += intSeconds;
                SetNewStep(pParseCommon,saveStep,lineCount,__LINE__);							
                break;
              default:
                printf("ERROR: unrecognized updateSextractorStep %d line %s lineCount %d\n",pPlateStats->updateSextractorStep,inLine,lineCount);
              }
            }
            pPlateStats->updateSextractorStep = 0;
          }

          nvals = sscanf(inLine,"copypipeline completed in %d",&intSeconds);
          if (nvals == 1) {
            if (pPlateStats->xfer1Complete == 0) {
              SetNewStep(pParseCommon,STEP_XFER1,lineCount,__LINE__);
            } else {
              SetNewStep(pParseCommon,STEP_XFER2,lineCount,__LINE__);
            }
            SetStepTime(pParseCommon,intSeconds,0,0,lineCount,__LINE__);
            pPlateStats->chargeTime += intSeconds;
            /* copypipeline completed in 3 seconds 1 forward 267 checks 74 copies 0 deleted for rh12283_01_01ww s2  441756932 bytes copied */
            /* copypipeline completed in 0 seconds 1 forward 89 checks 0 copies 0 deleted for bm00826_01_01ww s0apass  0 bytes copied */
            nvals = sscanf(inLine,"copypipeline completed in %d seconds %d forward",&intSeconds,&forwardFlag);
            if (nvals == 2) {
              /* Here we have valid bytes copied information */
              strcpy(parseLine,inLine);
              parsePtr1 = parseLine;
              parseCount = 0;
              while (1) {
                parsePtr2 = strtok(parsePtr1," ");
                parsePtr1 = NULL;
                if (parsePtr2 == NULL) {
                  break;
                }
                parseCount++;
                if (parseCount == 16) {
                  nvals = sscanf(parsePtr2,"s%d ",&tmpSolutionNumber);
                  if (nvals == 1) {
                    if (tmpSolutionNumber >= 0) {
                      pPlateStats->solutionNumber = tmpSolutionNumber;
                    }
                  } else {
                    break;
                  }
                }
                if (parseCount == 17) {
                  nvals = sscanf(parsePtr2,"%lld ",&tmpBytesCopied);
                  if (nvals == 1) {
                    if (pPlateStats->solutionNumber >= MAX_SOLUTION_NUMBER) {
                      printf("ERROR: MAX_SOLUTION_NUMBER too small for %d\n",pPlateStats->solutionNumber);
                      exit(-1);
                    }
                    pSizeStats = &pParseCommon->totalXferStats[pPlateStats->solutionNumber];
                    if (forwardFlag == 1) {
                      pSizeStats->forwardCount++;
                      pSizeStats->aveForwardBytesCopied += tmpBytesCopied;
                      if (pSizeStats->maxForwardBytesCopied < tmpBytesCopied) {
                        pSizeStats->maxForwardBytesCopied = tmpBytesCopied;
                      }
                      if (intSeconds > 0) {
                        tmpBytesCopied = tmpBytesCopied/(intSeconds+1);
                        pSizeStats->aveForwardBytesPerSecond += tmpBytesCopied;
                        if (pSizeStats->maxForwardBytesPerSecond < tmpBytesCopied) {
                          pSizeStats->maxForwardBytesPerSecond = tmpBytesCopied;
                        }
                      }
                    } else {
                      pSizeStats->reverseCount++;
                      pSizeStats->aveReverseBytesCopied += tmpBytesCopied;
                      if (pSizeStats->maxReverseBytesCopied < tmpBytesCopied) {
                        pSizeStats->maxReverseBytesCopied = tmpBytesCopied;
                      }
                      if (intSeconds > 0) {
                        tmpBytesCopied = tmpBytesCopied/(intSeconds+1);
                        pSizeStats->aveReverseBytesPerSecond += tmpBytesCopied;
                        if (pSizeStats->maxReverseBytesPerSecond < tmpBytesCopied) {
                          pSizeStats->maxReverseBytesPerSecond = tmpBytesCopied;
                        }
                      }

                    }
											
                  }
                  break;
                }
              }
            }

          }
          /* This format used for only one run! */
          nvals = sscanf(inLine,"copypipeline forward %d completed in %d",&forwardFlag,&intSeconds);
          if (nvals == 2) {
            if (forwardFlag == 1) {
              SetNewStep(pParseCommon,STEP_XFER1,lineCount,__LINE__);
            } else if (forwardFlag == 0) {
              SetNewStep(pParseCommon,STEP_XFER2,lineCount,__LINE__);
            } else {
              printf("ERROR: illegal forwardFlag %d line %s lineCount %d\n",forwardFlag,inLine,lineCount);
              SetNewStep(pParseCommon,STEP_XFER1,lineCount,__LINE__);
            }
            SetStepTime(pParseCommon,intSeconds,0,0,lineCount,__LINE__);
            pPlateStats->chargeTime += intSeconds;
          }

          nvals = sscanf(inLine,"magdepcalibrate seconds %d",&intSeconds);
          if (nvals == 1) {
            saveStep = pPlateStats->curStep;
            SetNewStep(pParseCommon,STEP_MAGDEPCALIBRATION,lineCount,__LINE__);
            SetStepTime(pParseCommon,intSeconds,0,0,lineCount,__LINE__);
            pPlateStats->chargeTime += intSeconds;
            SetNewStep(pParseCommon,saveStep,lineCount,__LINE__);
          }
				 
        }
      }
#endif /* DEBUG_SINGLE */
    }
  }
  pParseCommon->astrometrynetActive = 0;
  if (processCount == 0) {
    printf("ERROR: CHILD not found in %s\n",log_name);
    return(-1);
  }
  printf("Finished processing %s lines %d recStatCount %d\n",log_name,lineCount,recStatCount);
  fclose(log_handle);
  return(0);
}

int main(int argc,char *argv[])
{
  char *argstr;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  char list_name[MAX_INPUT_NAME];
  char out_name[MAX_INPUT_NAME];
  char searchclose3_name[MAX_INPUT_NAME];
  FILE *list_handle = NULL;
  int errorFlag = 0;
  char cmdchar;
  int lineLen;
  time_t curTime;
  PARSECOMMON parsecommon;
  PPARSECOMMON pParseCommon = &parsecommon;
  int stepIndex;
  PSTEPTABLE pStepTable;
  PSTEPSTATS pTotalStepStats;
  PSTEPSTATS pStepStats;
  time_t pTotalStepTime = 0;
  int stepCounter;
  double percentage;
  double percentage2;
  time_t averageTime;
  time_t totalAverageTime = 0;
  int solutionIndex;
  PSIZESTATS pSizeStats;
	PPLATESTATS pPlateStats = NULL;
  PSEARCHCLOSE3STATS pSearchClose3Stats = NULL;
  int plateStatsIndex;

#ifdef SEARCH_TABLE_LENGTH
  printf("ERROR: SEARCH_TABLE_LENGTH is defined at %d\n",SEARCH_TABLE_LENGTH);
#endif /* SEARCH_TABLE_LENGTH */
#ifdef COMBINE_CENTER_AREA
  printf("ERROR: COMBINE_CENTER_AREA is defined\n");
#endif /* COMBINE_CENTER_AREA */

  /* Loop through the arguments */
  list_name[0] = 0;
  out_name[0] = 0;
  memset(pParseCommon,0,sizeof(PARSECOMMON));

  for (stepIndex = 0; stepIndex < stepTableLength; stepIndex++) {
    pStepTable = &stepTable[stepIndex];
    if (pStepTable->sequence != pStepTable->stepmask) {
      printf("ERROR: step table sequence %d does not equal mask %d for '%s'\n",pStepTable->sequence,pStepTable->stepmask,pStepTable->stepname);
      exit(-1);
    }
    if (pStepTable->sequence != stepIndex) {
      printf("ERROR: step table sequence %d does not equal index %d for '%s'\n",pStepTable->sequence,stepIndex,pStepTable->stepname);
      exit(-1);
    }
#if 0
    printf("Step table sequence %2d for '%22s' '%s'\n",pStepTable->sequence,pStepTable->stepname,pStepTable->stepdescription);
#endif
  }

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: stray argument %s\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'v':
        case 'V':
          pParseCommon->verbose = 1;
          break;

        case 'a':
        case 'A':
          pParseCommon->astrometrynetFlag = 1;
          break;

        case 'p': /* search_close3.list file */
        case 'P':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(searchclose3_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            } else {
              pParseCommon->searchclose3Flag = 1;
            }
          }
          break;


        case 'r':
        case 'R':
          pParseCommon->reprocessFlag = 1;
          break;

        case 'l': /* list file name */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(list_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
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
            strncpy(out_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;


        default:
          printf("* illegal command -%c-",cmdchar);
          errorFlag = 1;
          break;
        }
        
      }
    }
  }



  if (list_name[0] == 0) {
    printf("ERROR: List file name not specified\n");
    errorFlag = 1;
  } else {
    list_handle = fopen(list_name,"rt");
    if (list_handle == NULL) {
      errorFlag = 1;
      printf("Could not open list file %s\n",list_name);
    }
  }

  if (out_name[0] == 0) {
    printf("ERROR: Out file name not specified\n");
    errorFlag = 1;
  } else {
    pParseCommon->out_handle = fopen(out_name,"wt");
    if (pParseCommon->out_handle == NULL) {
      errorFlag = 1;
      printf("Could not open out file %s\n",out_name);
    }
  }


  if (errorFlag) {
    printf("Usage: parselogs [-v] -l <list of log files> -o <outfile>\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -a is special processing for first instance of astrometry.net\n");
    printf("       where -r prints summary statistics only for the reprocessing states\n");

    return(-1);
  }

  printf("parselogs of %s %s List Filename %s Output Filename %s CHARGE_STEP_MASK 0x%x\n",
         __DATE__,__TIME__,list_name,out_name,CHARGE_STEP_MASK);
 

  time(&pParseCommon->startTime);

  if (pParseCommon->searchclose3Flag != 0) {
    fprintf(pParseCommon->out_handle,"series\tplateNumber\tmosaicNumber\tpatternID\tminTileOverlapCount\tmaxTileOverlapCount\tmaxMaskSector\ttotalMasks\tmaskCenterFlag\tdual_count_med1\tdual_count_rms1\tdual_count_limit1\tdual_maxBinCount1\tdual_minBinCount1\tdual_maxBinSNR1\tdual_totalMasksSector\tdual_count_limit2\tdual_maxMaskArea\tdual_maxCtrMaskArea\tsingle_count_med1\tsingle_count_rms1\tsingle_count_limit1\tsingle_maxBinCount1\tsingle_maxBinSNR1\tsingle_totalMasksSector\tsingle_count_limit2\tsingle_maxMaskArea\tsingle_maxCtrMaskArea\tsequence\n");
    fprintf(pParseCommon->out_handle,"------\t-----------\t------------\t---------\t-------------------\t-------------------\t-------------\t----------\t--------- ----\t---------------\t---------------\t-----------------\t-----------------\t-----------------\t---------------\t---------------------\t-----------------\t----------------\t-------------------\t-----------------\t-----------------\t-------------------\t-------------------\t-----------------\t-----------------------\t-------------------\t-------------------\t--------------------\t--------\n");
    
    ReadSearchClose3File(pParseCommon,searchclose3_name);
    pParseCommon->pPlateStatsTable = (PPLATESTATS)calloc(SEARCHCLOSE3_MAX_PLATES,sizeof(PLATESTATS));
    if (pParseCommon->pPlateStatsTable == NULL) {
      printf("ERROR: line %4d failed to allocate pPlateStatsTable\n");
      exit(-1);
    }

  }

  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,list_handle);
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
    if (inBuffer[0] != '#') {
      ProcessLogBuffer(pParseCommon,inBuffer);
    }
  }


  if (list_handle != NULL) {
    fclose(list_handle);
  }
  time(&curTime);
  curTime -= pParseCommon->startTime;

  for (stepCounter = 0; stepCounter < MAX_STEP; stepCounter++) {
    if ((pParseCommon->reprocessFlag != 0)  &&
        (((1 << stepCounter) & REPROCESS_MASK) == 0)) {
      /* Skip mosaic building and astrometry steps */
      continue;
    }


    pTotalStepStats = &pParseCommon->totalStepStats[stepCounter];
    if (pTotalStepStats->stepCount == 0) {
      averageTime = 0;
    } else {
      averageTime = pTotalStepStats->stepTime/pTotalStepStats->stepCount;
    }
    totalAverageTime += averageTime;
  }

  for (stepCounter = 0; stepCounter < MAX_STEP; stepCounter++) {
    if ((pParseCommon->reprocessFlag != 0)  &&
        (((1 << stepCounter) & REPROCESS_MASK) == 0)) {
      /* Skip mosaic building and astrometry steps */
      continue;
    }
    pTotalStepStats = &pParseCommon->totalStepStats[stepCounter];
    pStepTable = &stepTable[stepCounter];
    if (pParseCommon->totalStepTime == 0) {
      percentage = 0.0;
    } else {
      percentage = 100.0 * 0.00001 * (100000 * pTotalStepStats->stepTime/pParseCommon->totalStepTime);
    }
    if (pTotalStepStats->stepCount == 0) {
      averageTime = 0;
    } else {
      averageTime = pTotalStepStats->stepTime/pTotalStepStats->stepCount;
    }
    if (totalAverageTime == 0) {
      percentage2 = 0;
    } else {
      percentage2 = 100.0 * 0.00001 * (100000 * averageTime/totalAverageTime);
    }
    printf("Step %2d  %22s started %7d count %7d time %10lld %6.2f ave time %10lld %6.2f %s\n",stepCounter,pStepTable->stepname,pTotalStepStats->startedFlag,pTotalStepStats->stepCount,pTotalStepStats->stepTime,percentage,averageTime,percentage2,pStepTable->stepdescription);
  }
  
  if (pParseCommon->astrometrynetFlag != 0) {
    fprintf(pParseCommon->out_handle,"timeSpent\tplateCount\tcumulativeCount\tsequence\n");
    fprintf(pParseCommon->out_handle,"---------\t----------\t---------------\t--------\n");

    PrintAstrometryNetStats(pParseCommon,0);
    PrintAstrometryNetStats(pParseCommon,1);
  }
 


  for (solutionIndex = 0; solutionIndex < MAX_SOLUTION_NUMBER; solutionIndex++ ) {
    pSizeStats = &pParseCommon->totalXferStats[solutionIndex];
    if ((pSizeStats->forwardCount+pSizeStats->reverseCount) > 0) {
      if (pSizeStats->forwardCount > 0) {
        pSizeStats->aveForwardBytesCopied = pSizeStats->aveForwardBytesCopied/pSizeStats->forwardCount;
        pSizeStats->aveForwardBytesPerSecond = pSizeStats->aveForwardBytesPerSecond/pSizeStats->forwardCount;
        printf("copypipline solution %d forward count %7d forwardBytesCopied %10lld maxBytesCopied %10lld forwardBytes/Second %10lld maxForwardBytes/second %10lld\n",
               solutionIndex,
               pSizeStats->forwardCount,
               pSizeStats->aveForwardBytesCopied,
               pSizeStats->maxForwardBytesCopied,
               pSizeStats->aveForwardBytesPerSecond,
               pSizeStats->maxForwardBytesPerSecond);
      }
      if (pSizeStats->reverseCount > 0) {
        pSizeStats->aveReverseBytesCopied = pSizeStats->aveReverseBytesCopied/pSizeStats->reverseCount;
        pSizeStats->aveReverseBytesPerSecond = pSizeStats->aveReverseBytesPerSecond/pSizeStats->reverseCount;
        printf("copypipline solution %d reverse count %7d reverseBytesCopied %10lld maxBytesCopied %10lld reverseBytes/Second %10lld maxReverseBytes/second %10lld\n",
               solutionIndex,
               pSizeStats->reverseCount,
               pSizeStats->aveReverseBytesCopied,
               pSizeStats->maxReverseBytesCopied,
               pSizeStats->aveReverseBytesPerSecond,
               pSizeStats->maxReverseBytesPerSecond);
      }
    }

  }
  pSizeStats = &pParseCommon->totalMatchstarsStats;
  if (pSizeStats->forwardCount > 0) {
    pSizeStats->aveForwardBytesCopied = pSizeStats->aveForwardBytesCopied/pSizeStats->forwardCount;
    printf("matchstars     count %7d average maxMemory %10lld maximum maxMemory %10lld \n",
           pSizeStats->forwardCount,
           pSizeStats->aveForwardBytesCopied,
           pSizeStats->maxForwardBytesCopied);	
  }
  pSizeStats = &pParseCommon->totalFilterblendedStats;
  if (pSizeStats->forwardCount > 0) {
    pSizeStats->aveForwardBytesCopied = pSizeStats->aveForwardBytesCopied/pSizeStats->forwardCount;
    printf("filterblended  count %7d average maxMemory %10lld maximum maxMemory %10lld \n",
           pSizeStats->forwardCount,
           pSizeStats->aveForwardBytesCopied,
           pSizeStats->maxForwardBytesCopied);	
  }
  pSizeStats = &pParseCommon->totalRecoverPointsStats;
  if (pSizeStats->forwardCount > 0) {
    pSizeStats->aveForwardBytesCopied = pSizeStats->aveForwardBytesCopied/pSizeStats->forwardCount;
    printf("recover_points count %7d average maxMemory %10lld maximum maxMemory %10lld \n",
           pSizeStats->forwardCount,
           pSizeStats->aveForwardBytesCopied,
           pSizeStats->maxForwardBytesCopied);	
  }
  if (pParseCommon->searchclose3Flag != 0) {
    PrintSearchClose3Stats(pParseCommon);
    printf("plateStatsEntryCount is %5d\n",pParseCommon->plateStatsEntryCount);
    for (plateStatsIndex = 0; plateStatsIndex < pParseCommon->plateStatsEntryCount; plateStatsIndex++) {
      pPlateStats = &pParseCommon->pPlateStatsTable[plateStatsIndex];
      pSearchClose3Stats = &pPlateStats->searchclose3statstable;
      fprintf(pParseCommon->out_handle,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n",
              pPlateStats->series,
              pPlateStats->plateNumber,
              pPlateStats->mosaicNumber,
              pSearchClose3Stats->patternID,
              pSearchClose3Stats->minTileOverlapCount,
              pSearchClose3Stats->maxTileOverlapCount,
              pSearchClose3Stats->maxMaskSector,
              pSearchClose3Stats->totalMasks,    
              pSearchClose3Stats->maskCenterFlag,     /* 9 */
              pSearchClose3Stats->dual_count_med1,
              pSearchClose3Stats->dual_count_rms1,
              pSearchClose3Stats->dual_count_limit1,
              pSearchClose3Stats->dual_maxBinCount1,
              pSearchClose3Stats->dual_minBinCount1,
              pSearchClose3Stats->dual_maxBinSNR1,
              pSearchClose3Stats->dual_totalMasksSector,
              pSearchClose3Stats->dual_count_limit2,
              pSearchClose3Stats->dual_maxMaskArea,   
              pSearchClose3Stats->dual_maxCtrMaskArea,   /* 19 */
              pSearchClose3Stats->single_count_med1,
              pSearchClose3Stats->single_count_rms1,
              pSearchClose3Stats->single_count_limit1,
              pSearchClose3Stats->single_maxBinCount1,
              pSearchClose3Stats->single_maxBinSNR1,
              pSearchClose3Stats->single_totalMasksSector,
              pSearchClose3Stats->single_count_limit2,
              pSearchClose3Stats->single_maxMaskArea, 
              pSearchClose3Stats->single_maxCtrMaskArea,
              pSearchClose3Stats->sequence);  /* 29 */
    }

  }


  printf("Execution Time: %lld seconds. maxLineLen %d totalTimestampTime %lld totalChargeTime %lld totalStepTime %lld totalTime %lld totalAverageTime %lld\n",
         curTime,
         pParseCommon->maxLineLen,
         pParseCommon->totalTimestampTime,
         pParseCommon->totalChargeTime,
         pParseCommon->totalStepTime,
         pParseCommon->totalTime,
         totalAverageTime);
  
  if (pParseCommon->pPlateListEntryTable != NULL) {
    free(pParseCommon->pPlateListEntryTable);
  }
  if (pParseCommon->pPlateStatsTable != NULL) {
    free(pParseCommon->pPlateStatsTable);
  }
  if (pParseCommon->out_handle != NULL) {
    fclose(pParseCommon->out_handle);
  }



  return(EXIT_SUCCESS);
}
