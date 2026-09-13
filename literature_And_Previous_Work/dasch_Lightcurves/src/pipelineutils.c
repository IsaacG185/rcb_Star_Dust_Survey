// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* pipelineutils.c
 *
 *  Includes code from statstable.c,
 *                     Starbase Data Tables - An Ascii Database for UNIX
 *                     copyright 1991, 1993, 1995, 1999 John B. Roll jr.
 *  gcc -ggdb -c -O2 -I/usr/include/mysql pipelineutils.c  -o pipelineutils.o
 *
 * For debugging:
 *  gcc -ggdb   -DGENERATE_ARRAY pipelineutils.c -lmysqlclient  -o pipelineutils -I/usr/include/mysql -L/usr/lib${lib64}/mysql  -l mysqlclient -L /dasch/install/lib  -lwcs -lm -ldl -pthread
 *
 *  This subroutine takes the horizontal and vertical width in pixels and a pixel location
 *  and calculates the bin number.   The subroutine is modelled after Sumin Tang's code
 *  in divide_annul9.m
 *
 * Jan 16, 2008 Edward J. Los - Initial Version
 * Apr 24, 2008 Edward J. Los - have CalcMedianAndRMS return the clipped number of stars
 * May 20, 2008 Edward J. Los - Add clipFactor to CalcMedianAndRMS
 * Aug 12, 2008 Edward J. Los - Add zeroFlag to CalcMedianAndRMS
 * Nov 24, 2008 Edward J. Los - Add binFractions global
 * Nov 25, 2008 Edward J. Los - Add MaxDradPixels() with variable tolerance dependent on series
 * Dec  8, 2008 Edward J. Los - Add reject reasons from extract_lightcurves
 * Feb 16, 2009 Edward J. Los - Add CheckBlend and rotateFrame
 * Feb 21, 2009 Edward J. Los - Currect segmententation fault in CalcMedianAndRMS for no valid points
 * Feb 24, 2009 Edward J. Los - Add DecodeAFLAGS
 * Feb 27, 2009 Edward J. Los - Add CalcMemory
 * Mar  9, 2009 Edward J. Los - Add ComputeRange
 * Apr  4, 2009 Edward J. Los - Flag objects without colorterm correction in plots
 * May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 * May 26, 2009 Edward J. Los - Correct GetExposureInfo bug where exposure is not retrieved if timeAccuracy is not available.
 * Jun  3, 2009 Edward J. Los - Tighten distortion tolerances for i plates from 5 to 3 pixels.
 * Jun 23, 2009 Edward J. Los - Add GetMargins
 * Jul 20, 2009 Edward J. Los - Add ChargeTime for performance tracking.
 * Aug 13, 2009 Edward J. Los - Add Kepler support
 * Aug 31, 2009 Edward J. Los - Add solutionNumber
 *                              Move ID conversion routines here
 * Nov 16, 2009 Edward J. Los - Add FILTER_AFLAG_UNCERTAIN_DATE and FILTER_AFLAG_MULTIPLE_BLEND
 * Dec  7, 2009 Edward J. Los - Add FILTER_AFLAG_MULTIPLE_NONE
 * Dec 15, 2009 Edward J. Los - Add SetMosaicFitWCS to keep track of pipeline failures
 * Dec 22, 2009 Edward J. Los - Add SetPlateQuality to keep track of second-quality plates
 * Jan 19, 2010 Edward J. Los - Add elevation to the locations table and GetPlateColor()
 *                            - Implement a new extinction formula for GetExtinctionCoefficient()
 *                                blue   plates: coefficient = 0.4  - (0.2 *elevation/2300m)
 *                                yellow plates: coefficient = 0.25 - (0.11*elevation/2300m)
 *                                red    plates: coefficient = 0.2  - (0.1 *elevation/2300m)
 * Jan 22, 2010 Edward J. Los - Add plate quality support
 * Feb 23, 2010 Edward J. Los - Flag good stars for with sextractor blend or neighbor flags have been set
 * Apr 12, 2010 Edward J. Los - Add GetSeriesInfo for VOTable support
 * May  4, 2010 Edward J. Los - Add GetMysqlError()
 * May 11, 2010 Edward J. Los - Add REF_TYPE fields and GetDASCHCoordinates
 * Jun  4, 2010 Edward J. Los - Move AflagsTable here from showflags
 * Oct 26, 2010 Edward J. Los - Add InaccuratePV flag
 * Nov 29, 2010 Edward J. Los - Implement new SCAMP drad limits.
 * Jan 28, 2011 Edward J. Los - Enhance debug mode to search for smoothing bins in two annular bins.
 * Feb 22, 2011 Edward J. Los - Generate binarray.h and binarray.m
 * Mar 14, 2011 Edward J. Los - Back out binarray.m because of differences with binarray.h calculation
 * Mar 23, 2011 Edward J. Los - Add APASS support
 * Apr  8, 2011 Edward J. Los - Add GetREFType
 * Jun 27, 2011 Edward J. Los - Add SelectBestMosaic
 * Aug 12, 2011 Edward J. Los - Add flatfields table names
 * Aug 19, 2011 Edward J. Los - Add Tycho2 reference numbers
 * Aug 20, 2011 Edward J. Los - Correct incorrect type setting for APASS objects.
 * Aug 31, 2011 Edward J. Los - Return "UNKNOWN" for soft failures to get the reference number.
 * Aug 31, 2011 Edward J. Los - Temporarily remove the REFNumber from the header sanity check (TEMP_REFNUMBER_FIX)
 * Oct 24, 2011 Edward J. Los - Correct median calculation for odd vectors
 * Oct 30, 2011 Edward J. Los - Add magnitude-dependent star correction support.
 * Nov 28, 2011 Edward J. Los - Add REJECT_REASON_SATURATED and SYMBOL_SATURATED
 *                              Change order from "q", "c", "b" to  "q", "s", "c", and "b"
 * Dec 10, 2011 Edward J. Los - Add SYMBOL_NOMAGDEP ("m") and REJECT_REASON_NOMAGDEP.  Order is "q", "m", "s", "c", and "b"
 * Feb 17, 2012 Edward J. Los - Add multiple exposure mask support
 * Feb 24, 2012 Edward J. Los - Add UpdateQuality
 * Mar 13, 2012 Edward J. Los - Add CheckAuthorization
 * Mar 23, 2012 Edward J. Los - Support plot selection by series
 * Jun  4, 2012 Edward J. Los - Add "LARGEA_TEST"
 * Jun 18, 2012 Edward J. Los - Correct SelectBestMosaic for unfitted plates
 * Jun 18, 2012 Edward J. Los - Set illegal catalog declination 99.0 and RA to 999.0
 * Jun 29, 2012 Edward J. Los - Add a flag to GetExposureInfo to get the RA and DEC from the logbook transcriptions.
 * Jul  3, 2012 Edward J. Los - Add high background object support
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Oct  6, 2012 Edward J. Los - Improve GetSubbins performance
 * Oct 19, 2012 Edward J. Los - Always assert nospace for GetREF
 * Nov 23, 2012 Edward J. Los - Correct GetDASCHCoordinates for id table REF format
 * Feb 25, 2013 Edward J. Los   Add support for running astrometry.net on full sized mosaics
 * Jul  3, 2013 Edward J. Los - Correct a ninety degree angle error in converting THETA_J2000 to THETA_IMAGE for CheckBlend()
 * Aug 12, 2013 Edward J. Los - Add the UCAC4 REFNumber type 6.
 * Sep  1, 2013 Edward J. Los - Have CheckAuthorization give priority to the official releases over the calibration fields
 * Jun  3, 2014 Edward J. Los - Support DR3
 * Oct  7, 2014 Edward J. Los - Recogize "KIC" as a Kepler Input Catalog prefix
 * Dec 22, 2014 Edward J. Los   Support DR4
 * Jan 12, 2015 Edward J. Los - In GetExposureInfo, if the exposure is from the card catalog and the date contains the string "00:00:00.0", set the timeAccuracy to 1.0 days.
 * Mar  4, 2015 Edward J. Los - Reorder release fields to progress incrementally from the NGP to the SGP
 *                              Old Release Order:  NGP to SGP  DR1 DR2 DR3 DR4 DR5 DR11 DR12 DR10 DR9 DR8  DR7  DR6
 *                              New Release Order:              DR1 DR2 DR3 DR4 DR5 DR6  DR7  DR8  DR9 DR10 DR11 DR12
 * Mar 10, 2015 Edward J. Los - Correct CheckBlend to test the origin for overlapping similar objects with 45 degrees THETA shift
 * May 15, 2015 Edward J. Los - Add bit numbers to  AflagsTable
 * May 16, 2015 Edward J. Los   Add formatAflags to list individual bits set in AFLAGS
 * Sep 28, 2015 Edward J. Los   Add EXCLUSION table to study transient candidates favoring particular regions of the sky
 * Oct  5, 2015 Edward J. Los   Add local_strlwr
 * Oct 16, 2015 Edward J. Los   Add FITWCS_NOAPASSALLOBJECTS
 * Nov 17, 2015 Edward J. Los   Support DR5
 * Feb 16, 2015 Edward J. Los   Change FILTER_AFLAG_QUALITY description from "Multiple exposure, grating, spectra,
 *                              Pickering Wedge, yellow or red plate (7)" to "Plate-wide issue noted in Quality Flags (7)"
 * Feb 22, 2016 Edward J. Los - In formatAflags, split AFLAGS into two columns: AFLAGS and AFLAGSBits
 * Feb 24, 2016 Edward J. Los - Change "formatAflags" to "FormatFlagsBits" for use with BFLAGS and quality bitmaps and eliminate the "BIT" text
 * Mar  1, 2016 Edward J. Los - Add all of the quality bits to the showflags output
 * Mar  8, 2016 Edward J. Los - Use two separate fields for FormatFlagsBits
 * Sep 30, 2016 Edward J. Los - define REJECTBIT for REJECT_REASON text and add GetRejectReasonText()
 *                              introduce FATAL_REASON and add GetFatalReasonText()
 * Dec 23, 2016 Edward J. Los - Add BlueColorterm
 * Feb 28, 2017 Edward J. Los - correct status return for SetMosaicFitWCS
 * Mar 28, 2017 Edward J. Los   Add peakEvaluation bitmap to indicate reasons for rejection of transient candidate flares
 * Apr  4, 2017 Edward J. Los   Split PEAKEVALUATION_DATERANGE into PEAKEVALUATION_TOOSHORTDATERANGE, PEAKEVALUATION_TOOLONGDATERANGE, and PEAKEVALUATION_PEAKEVALUATION_INSUFFICIENT_POINTS3
 *                              Add PEAKEVALUATION_EXCLUSIONZONE, PEAKEVALUATION_INSUFFICIENT_POINTS2, and PEAKEVALUATION_TOOMANYPOINTS
 * Apr 25, 2017 Edward J. Los   Move WEDGEENTRY from filter_wedge for sharing between filter_wedge.c and update_quality.c
 * May 23, 2017 Edward J. Los   Increase the "ax" series from 4.3 to 7.1 pixel in MaxDradPixels
 * Jun 28, 2017 Edward J. Los   Add maskKernelType to masks
 * Jul 25, 2017 Edward J. Los   Make CalculateBin return 1 spatial bin for 4"x5" or smaller plates
 *                              Add XDmagBins() and YDmagBins() for dynamic bin calculation
 *                              GetMosaicInfo() will now return the number of smoothing bins
 * May 14, 2017 Edward J. Los   Add ReadLinearityFile() from ./FitsWrapper/ScannerCCD.cpp
 * May 21, 2018 Edward J. Los   Add REF_TYPE_GAIA2 for Gaia DR2 (source_id numbers are not guaranteed to be unique across releases)
 * Oct 23, 2018 Edward J. Los   Add REF_TYPE_ATLAS2 for the Atlas 2 catalog
 * Mar  2, 2019 Edward J. Los   Correct DR7 range (-15 to 0) instead of (0 to 15)
 * Jul 12, 2019 Edward J. Los   Add GetPlateCondition and GetPlateEvent for stacks database access
 * Oct  4, 2019 Edward J. Los   Add remaining release level
 * Nov 30, 2019 Edward J. Los   Add GetFlatfieldRecords
 * Dec  3, 2019 Edward J. Los   Add scanDate to GetMosaicInfo()  There are 4 pathological cases where the scanNumber does not reflect the order of scanDate: dny00266, dny00386, j04718, and mc39624.
 * Mar 19, 2019 Edward J. Los   Add patternID to GetMosaicInfo()
 * Apr 22, 2019 Edward J. Los   Add BuildMosaicList() go make a list of the best mosaics
 *                              Add orientation and fittedPlateScale to SERIES_HEADER
 *                              Add checks for incorrect centerSource positions.
 *                              Correct the RA and DEC strings to read "N.A." instead of a bogus "18:36:00.000 +99:00:00.00"
 * Apr  9, 2020 Edward J. Los   Add mosaics.JobId to MOSAIC (MOSAICLIST) for compatabilty with Catalog.h
 * Apr 21, 2020 Edward J. Los   Add MySQL commands to correct inconsistent data
 * May 26, 2020 Edward J. Los   Correct initialization of numExposureRecords in BuildMosaicList()
 * May 30, 2020 Edward J. Los   Correct invocation of GetSeriesInfo from BuildMosaicList() and handling of the "mask" series.
 * May 25, 2020 Edward J. Los   Add LOGBOOK_VERSION for dasch_scat consistency checks
 * Jan 26, 2021 Edward J. Los   Define PEAKEVALUATION_PRIMARY_GSC_BIN and PEAKEVALUATION_SECONDARY_GSC_BIN to
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"
#include "pipelineutils.h"


#define DEBUG_SERIES ""
#define DEBUG_LOWPLATENUMBER     -1
#define DEBUG_HIGHPLATENUMBER    -1
double wcsdist(
               double x1,
               double y1,	/* (RA,Dec) or (Long,Lat) in degrees */
               double x2,
               double y2);	/* (RA,Dec) or (Long,Lat) in degrees */


/* #define LARGEA_TEST 1 */
#ifdef GENERATE_ARRAY
/* #define USE_BINARRAY 1   For debug only */
#else /* GENERATE_ARRAY */
#define USE_BINARRAY 1
#endif /* GENERATE_ARRAY */

#ifdef USE_BINARRAY
#include "binarray.h"
#endif

double fd2jd (char *string);
void ra2str (

             char	*string,	/* Character string (returned) */
             int	lstr,		/* Maximum number of characters in string */
             double	ra,		/* Right ascension in degrees */
             int	ndec);		/* Number of decimal places in seconds */

void
dec2str (

         char	*string,	/* Character string (returned) */
         int	lstr,		/* Maximum number of characters in string */
         double	dec,		/* Declination in degrees */
         int	ndec);		/* Number of decimal places in arcseconds */

double str2dec(		        /* Return Dec in degrees from string */
               const char* in);	/* Character string (dd:mm:ss.sss or dd.dddd) */
double str2ra(		        /* Return RA in degrees from string */
              const char* in);	/* Character string (hh:mm:ss.sss or dd.dddd) */

/* #define WRITE_FILES 1 */

PRANGES pRangeTable = NULL; /* Dump of the range table */
int numRanges = 0;        /* Number of range table entries */

PLOCATIONS pLocationTable = NULL; /* Dump of the location table */
int numLocations = 0;          /* Number of location table entries */

char *releaseFieldText[RELEASE_FIELD_MAX+1] =
  {
    "DR1",
    "DR2",
    "DR3",
    "DR4",
    "DR5",
    "DR6",
    "DR7",
    "DR8",
    "DR9",
    "DR10",
    "DR11",
    "DR12",
    "M44",
    "3C273",
    "BAADE",
    "KEPLER",
    "LMC",
    "Other",
    "ALL"
  };
/* Note: change MAX_QUALIFIER if lengths increase */
char *catalogText[MAX_CATALOG_NUMBER] =
  {"gsc2.3.2",
   "kepler",
   "apass",
   "gaia",
   "atlas",
   "experimental"
  };



/* Note: Any change in the annular bin radii or definitions requires
 *  reruning pipelineutils in debug mode and modification of this table
 */
#ifdef LARGEA_TEST
double binFractions[MAX_SPATIAL_BINS+1] = {
  0.0,
  0.04640,  /* bin 1 */
  0.04800,  /* bin 2 */
  0.04160,  /* bin 3 */
  0.04800,  /* bin 4 */
  0.04400,  /* bin 5 */
  0.04560,  /* bin 6 */
  0.04480,  /* bin 7 */
  0.04800,  /* bin 8 */
  0.63360,  /* bin 9 */
};

#else /* LARGEA_TEST */
double binFractions[MAX_SPATIAL_BINS+1] = {
  0.0,
  0.09025,  /* bin 1 */
  0.09025,  /* bin 2 */
  0.09025,  /* bin 3 */
  0.09025,  /* bin 4 */
  0.09025,  /* bin 5 */
  0.11344,  /* bin 6 */
  0.10973,  /* bin 7 */
  0.14302,  /* bin 8 */
  0.18256,  /* bin 9 */
};
#endif /* LARGEA_TEST */

/* NOTE: the indices into this array must agree with KEPLER_* in pipelineutils.h */
char *keplerSourceText[KEPLER_CQ_MAGFLAG_MAX] =
  {"SCP",
   "2MASS",
   "PHOTO",
   "TYBV",
   "UNCAL",
   "NOCAL"
  };

char *excludeSeriesList[NUM_EXCLUDE_SERIES] = {"ac","ca","am","ax"};

char *pickedCommentString[PICKED_COMMENT_MAX] = {
  "Broken",
  "Masked",
  "Markings",
  "Trailed",
  "Multiple",
  "Discrepancy",
  "Unused2",
  "Out of Focus",
  "Needs Mending",
  "Point Delamination",
  "Yellow Filter",
  "Red Filter",
  "Isochromatic",
  "Panchromatic",
  "Scratched",
  "Scan and Review",
  "Coarse Grating",
  "Fine Grating",
  "Cracked"
};

char *pickedStatusString[PICKED_STATUS_MAX] = {
  "Unknown",
  "Picked",
  "Spectra",
  "Grating XXXXXX", /* Deprecated.  Note: the 'Grating' element can cause false matches with 'Coarse Grating' and 'Fine Grating' */
  "Long Trails",
  "Many Exposures",
  "Missing",
  "Discarded",
  "Too Dark",
  "Too Thick",
  "14x17",
  "10x12",
  "Unused1",
  "Plate JPEG",
  "Jacket JPEG",
  "Scanned",
  "Severe Delamination",
  "Multi Direction Trailing",
  "Rejected",
  "Reseau Grid",
  "On Shelf",
  "Hand Cleaned",
  "Machine Cleaned",
  "Reshelved",
  "Coarse Grating",
  "Fine Grating",
  "Reseau Grid",
  "Temp Location",
  "Perm Location",
  "No Transcription",
  "Flood"
};



#if 0
char *rejectReasonText[REJECT_REASON_MAX] =
  {"Error processing this Plate          ",
   "Object is not on this plate          ",
   "Object in a failed bin               ",
   "Object undetected on this plate      ",
   "Object is a pickering wedge          ",
   "Object has high DRAD                 ",
   "Object failed the defect filter      ",
   "Object in high zout spatial bin      ",
   "Object in glare of brighter object   ",
   "Object is blended                    ",
   "Object is in bin 9                   ",
   "Object has high local calibration rms",
   "Object is close to limiting magnitude",
   "Object is too bright                 ",
   "Object is too close to the horizon   ",
   "Object has high lowess rms           ",
   "Object in a high drad bin            ",
   "Object has case c or d blend         ",
   "Object is a multiple blend           ",
   "Object has uncertain date            ",
   "Unmatched multiple exposure object   ",
   "Plate fails quality criteria         ",
   "Object is saturated                  ",
   "No magnitude-dependent correction",
   "Object has a high background level",
  };
#endif

REJECTBIT rejectReasonBits[REJECT_REASON_MAX] = {
  {REJECT_REASON_ERROR,               "Error processing this Plate             "},
  {REJECT_REASON_OFFPLATE,            "Object is not on this plate             "},
  {REJECT_REASON_BINFAILURE,          "Object in a failed bin                  "},
  {REJECT_REASON_UNDETECTED,          "Object undetected on this plate         "},
  {REJECT_REASON_WEDGE,               "Object is a pickering wedge             "},
  {REJECT_REASON_DRAD,                "Object has high DRAD                    "},
  {REJECT_REASON_DEFECT,              "Object failed the defect filter         "},
  {REJECT_REASON_HIZOUT,              "Object in high zout spatial bin         "},
  {REJECT_REASON_BLEND_NOMATCH,       "Object in glare of brighter object      "},
  {REJECT_REASON_BLEND,               "Object is blended                       "},
  {REJECT_REASON_BIN9,                "Object is in bin 9                      "},
  {REJECT_REASON_LOCAL_RMS,           "Object has high local calibration rms   "},
  {REJECT_REASON_LIMITING_MAG,        "Object is close to limiting magnitude   "},
  {REJECT_REASON_TOO_BRIGHT,          "Object is too bright                    "},
  {REJECT_REASON_LOW_ALTITUDE,        "Object is too close to the horizon      "},
  {REJECT_REASON_ISO_RMS,             "Object has high lowess rms              "},
  {REJECT_REASON_DRADBIN,             "Object in a high drad bin               "},
  {REJECT_REASON_COMPLEX_BLEND,       "Object has case c or d blend            "},
  {REJECT_REASON_MULTIPLE_BLEND,      "Object is a multiple blend              "},
  {REJECT_REASON_UNCERTAIN_DATE,      "Object has uncertain date               "},
  {REJECT_REASON_MULTIPLE_NONE,       "Unmatched multiple exposure object      "},
  {REJECT_REASON_QUALITY,             "Plate fails quality criteria            "},
  {REJECT_REASON_SATURATED,           "Object is saturated                     "},
  {REJECT_REASON_NOMAGDEP,            "No magnitude-dependent correction       "},
  {REJECT_REASON_BACKGROUND,          "Object has a high background level      "},

};


REJECTBIT fatalReasonBits[FATAL_REASON_MAX] = {
  {FATAL_REASON_OFFPLATE,            "Object is not on this plate             "},
  {FATAL_REASON_NOPHOTPLATEENTRY,    "No photplate entry found                "},
  {FATAL_REASON_NOALLOBJECTS,        "No photometry result for this plate     "},
  {FATAL_REASON_ZEROPHOTPLATEVERSION,"Photplate entry has a zero versionId    "},
  {FATAL_REASON_NOMOSAICENTRY,       "No mosaic entry found                   "},
  {FATAL_REASON_NOJULIANDATE,        "No Julian Date found                    "},
  {FATAL_REASON_NOPLATESENTRY,       "No plates table entry found             "},
  {FATAL_REASON_NOSPATIALBINTABLE,   "No spatial bin table found              "},
  {FATAL_REASON_NOSPATIALBINENTRY,   "No spatial bin table entry found        "},
  {FATAL_REASON_ILLEGALSPATIALBIN,   "Illegal spatial bin found               "},
  {FATAL_REASON_NOLOCALBINTABLE,     "No local bin table found                "},
  {FATAL_REASON_NOLOCALBINENTRY,     "No local bin entry found                "},
  {FATAL_REASON_LOCALBINFAILURE,     "Local bin failed                        "},

  {FATAL_REASON_BADERRORFLAG,        "Unexpected set errorFlag                "},
  {FATAL_REASON_BADVALIDPOINT,       "Unexpected clear validPoint             "},
  {FATAL_REASON_NONE,                "No rejection                            "},
};

QUALITYBIT peakEvaluationMasks[PEAKEVALUATION_MAX+1] = {
  { 1 << PEAKEVALUATION_ENTRY,                   "ENTRY"},
  { 1 << PEAKEVALUATION_INSUFFICIENT_POINTS,     "INSUFFICIENT_POINTS"},
  { 1 << PEAKEVALUATION_NOPEAKS,                 "NOPEAKS"},
  { 1 << PEAKEVALUATION_DATERANGE,               "DATERANGE"},
  { 1 << PEAKEVALUATION_MULTIPLE1,               "MULTIPLE1"},
  { 1 << PEAKEVALUATION_SOFTWARE,                "SOFTWARE"},
  { 1 << PEAKEVALUATION_LOW_AMPLITUDE,           "LOW_AMPLITUDE"},
  { 1 << PEAKEVALUATION_NEARBY_DIMMAG,           "NEARBY_DIMMAG"},
  { 1 << PEAKEVALUATION_NEARBY_AVEMAG,           "NEARBY_AVEMAG"},
  { 1 << PEAKEVALUATION_NEARBY_DIMMAG2,          "NEARBY_DIMMAG2"},
  { 1 << PEAKEVALUATION_NEARBY_AVEMAG2,          "NEARBY_AVEMAG2"},
  { 1 << PEAKEVALUATION_NEARTYCHO2,              "NEARTYCHO2"},
  { 1 << PEAKEVALUATION_FINISHED,                "FINISHED"},
  { 1 << PEAKEVALUATION_ACCEPTED,                "ACCEPTED"},
  { 1 << PEAKEVALUATION_TOOSHORTDATERANGE,       "TOOSHORTDATERANGE"},
  { 1 << PEAKEVALUATION_TOOLONGDATERANGE,        "TOOLONGDATERANGE"},
  { 1 << PEAKEVALUATION_INSUFFICIENT_POINTS3,    "INSUFFICIENT_POINTS3"},
  { 1 << PEAKEVALUATION_EXCLUSIONZONE,           "EXCLUSIONZONE"},
  { 1 << PEAKEVALUATION_INSUFFICIENT_POINTS2,    "INSUFFICIENT_POINTS2"},
  { 1 << PEAKEVALUATION_TOOMANYPOINTS,           "TOOMANYPOINTS"},
  { 1 << PEAKEVALUATION_LONG_INSUFFICIENT_POINTS,"LONG_INSUFFICIENT_POINTS"},
  { 1 << PEAKEVALUATION_LONG_ACCEPTED,           "LONG_ACCEPTED"},
  { 1 << PEAKEVALUATION_LONG_TOO_MANY_EVENTS,    "LONG_TOO_MANY_EVENTS"},
  { 1 << PEAKEVALUATION_LONG_REJECTED,           "LONG_REJECTED"},
  { 1 << PEAKEVALUATION_PRIMARY_GSC_BIN,         "PRIMARY_GSC_BIN"},
  { 1 << PEAKEVALUATION_SECONDARY_GSC_BIN,       "SECONDARY_GSC_BIN"},
  { 0,                                           "END"},
};
int peakEvaluationTableSize = sizeof(peakEvaluationMasks)/sizeof(QUALITYBIT);

FITWCSBIT FitWCSMasks[] = {
  {FITWCS_SELECTED,            "Selected"},              /* Selected for fitting */
  {FITWCS_MIKESHAWFIT,         "MikeShawFit"},           /* Previously fitted by Mike Shaw (deprecated) */
  {FITWCS_COMPLETED,           "Completed"},             /* Fitting process successful */
  {FITWCS_COPIED,              "Copied"},                /* Fit results copied to head cluster (deprecated) */
  {FITWCS_MULTFAILEDASTROMETRY,"MultFailedAstrometry"},  /* Fit results failed AstrometryWCS for next exposure */
  {FITWCS_MULTFAILEDSEPARATION,"MultFailedSeparation"},  /* Fit results failed the separation criteria for the next exposure */
  {FITWCS_MULTSUCCEEDED,       "MultSucceeded"},         /* Fit results succeeded for the next exposure */
  {FITWCS_STALESOLUTION,       "StaleSolution"},         /* Better astrometry invalidates this solution */
  {FITWCS_STALEEXPOSURENUMBER, "StaleExposureNumber"},   /* The exposure number of this solution has been changed */
  {FITWCS_WEDGE,               "Wedge"},                 /* This is a Pickering Wedge plate */
  {FITWCS_FILTERBLENDEDTIMEOUT,"FilterblendedTimeout"},  /* The filterblended routine timed out */
  {FITWCS_NOALLOBJECTS,        "NoAllobjects"},          /* The allobjects file for this solution is missing */
  {FITWCS_DUPLICATEREFCOUNT,   "DuplicateREFCount"},     /* This solution shares the same stars with a previous solution */
  {FITWCS_COLORTERMCRASH,      "ColortermCrash"},        /* The colorterm.m script has crashed and needs to be rerun */
  {FITWCS_INACCURATEPV,        "InaccuratePV"},          /* The SCAMP polynomial fit has significant inaccuracy */
  {FITWCS_NOMULTMASK,          "NoMultipleMask"},        /* No Multiple Exposure Mask found */
  {FITWCS_HAVEMULTMASK,        "HaveMultipleMask"},      /* Have a multiple exposure mask */
  {FITWCS_FULLFAILEDASTROMETRY,"FullFailedAstrometry"},  /* Astrometry.net failed on the full size mosaic */
  {FITWCS_FULLFAILEDSCAMPUSNO, "FullFailedScampUSNO"},   /* SCAMP using Astrometry.net catalogs failed on the full size mosaic */
  {FITWCS_FULLFAILEDSCAMPUCAC, "FullFailedScampUCAC"},   /* SCAMP using the UCAC catalog failed on the full size mosaic */
  {FITWCS_FULLFAILEDTOLERANCE, "FullFailedTolerance"},   /* SCAMP did not produce a better solution with the full size mosaic */
  {FITWCS_FULLSUCCEEDED,       "FullSucceeded"},         /* SCAMP succeeded on the full size mosaic */
  {FITWCS_ASTROMETRY2,         "Astrometry2Succeeded"},  /* find_astrometry2 produced a result (for the next solution) */
  {FITWCS_NOAPASSALLOBJECTS,    "NoAPASSAllobjects"},    /* The APASS allobjects file for this solution is missing */
  {0,""}                                                 /* Terminator */
};

QUALITYBIT qualityMasks[] = {
  {QUALITY_MULTIPLE, "multiple"},  /* Multiple exposure */
  {QUALITY_GRATING,  "grating"},   /* Grating exposure */
  {QUALITY_COLOR,    "color"},     /* Color filter used */
  {QUALITY_COLORTERM,"colorterm"}, /* Fails colorterm limits */
  {QUALITY_WEDGE,    "wedge"},     /* Pickering Wedge plate */
  {QUALITY_SPECTRA,  "spectra"},   /* Spectra plate */
  {QUALITY_SATURATED,"saturated"},   /* Saturated images */
  {QUALITY_NOMAGDEP,"nomagdep"}, /* Not magnitude-dependent corrected */
  {QUALITY_PATROL,   "patrol"},   /* nominal scale >= PATROL_PLATE_SCALE */
  {QUALITY_NONPATROL,"nonpatrol"},   /* nominal scale < PATROL_PLATE_SCALE */
  {QUALITY_LIMITING,"limiting"},   /* show limiting magnitudes */
  {QUALITY_UNDETECTED,"undetected"},   /* show points not detected */
  {QUALITY_TRAILED,"trailed"},   /* ELLIPTICITY exceeds MAX_ELLIPTICITY = 0.6 */
  {0,""}                           /* Terminator */
};
int qualityTableSize = sizeof(qualityMasks)/sizeof(QUALITYBIT);


QUALITYBIT webQualityMasks[] = {
  {QUALITY_MULTIPLE, "Multiple Exposure Plates"},  /* Multiple exposure */
  {QUALITY_GRATING,  "Grating Plates"},   /* Grating exposure */
  {QUALITY_COLOR,    "Yellow or Red Plates"},     /* Color filter used */
  {QUALITY_COLORTERM,"Non-Blue Colorterm Spatial Bin"}, /* Fails colorterm limits */
  {QUALITY_WEDGE,    "Pickering Wedge Plates"},     /* Pickering Wedge plate */
  {QUALITY_SPECTRA,  "Spectra Plates"},   /* Spectra plate */
  /* The following fields are not stored in the database but used only by the web site */
  {QUALITY_SATURATED,"Saturated Images"},   /* Saturated images */
  {QUALITY_NOMAGDEP,"Not Magnitude-dependent Corrected"}, /* Not magnitude-dependent corrected */
  {QUALITY_PATROL,   "Wide-Field Patrol Telescopes"},   /* nominal scale >= PATROL_PLATE_SCALE */
  {QUALITY_NONPATROL,"Narrow-Field Telescopes"},   /* nominal scale < PATROL_PLATE_SCALE */
  {QUALITY_LIMITING,"Limiting Magnitudes"},   /* show limiting magnitudes */
  {QUALITY_UNDETECTED,"Undetected"},   /* show points not detected */
  {QUALITY_TRAILED,"Trailed"},   /* ELLIPTICITY exceeds MAX_ELLIPTICITY = 0.6 */
  {0,""}                           /* Terminator */
};
int webQualityTableSize = sizeof(webQualityMasks)/sizeof(QUALITYBIT);


FLAGSENTRY AflagsTable[] = {
  {FILTER_AFLAG_ISO_RMS,0,"FILTER_AFLAG_ISO_RMS","Maximum isophotonic rms exceeded (11)"},
  {FILTER_AFLAG_LOCAL_RMS,0,"FILTER_AFLAG_LOCAL_RMS","Maximum locally smoothed rms exceeded (12)"},
  {FILTER_AFLAG_LIMITING_MAG,0,"FILTER_AFLAG_LIMITING_MAG","Within 0.5 magnitude of the limiting magnitude (13)"},
  {FILTER_AFLAG_BIN9,0,"FILTER_AFLAG_BIN9","In spatial bin 9 (14)"},
  {FILTER_AFLAG_DRADBIN,0,"FILTER_AFLAG_DRADBIN","Bin has high or unknown astrometric error (15)"},
  {GSC_VARIABLE_BIT,0,"GSC_VARIABLE_BIT","Uncertain catalog brightness (19)"},
  {FILTER_AFLAG_CASEB,0,"FILTER_AFLAG_CASEB","Blend of multiple catalog stars (20)"},
  {FILTER_AFLAG_CASEC,0,"FILTER_AFLAG_CASEC","Multiple images for one catalog star (21)"},
  {FILTER_AFLAG_CASED,0,"FILTER_AFLAG_CASED","Complex blend (22)"},
  {FILTER_AFLAG_DRAD,0,"FILTER_AFLAG_DRAD","High astrometric error (23)"},
  {FILTER_AFLAG_WEDGE,0,"FILTER_AFLAG_WEDGE","Pickering Wedge Image (24)"},
  {FILTER_AFLAG_DEFECT,0,"FILTER_AFLAG_DEFECT","Plate Defect (25)"},
  {FILTER_AFLAG_BLEND,0,"FILTER_AFLAG_BLEND","SExtractor blend (26)"},
  {FILTER_AFLAG_BLEND_NOMATCH,0,"FILTER_AFLAG_BLEND_NOMATCH","Rejected blended object (27)"},
  {FILTER_AFLAG_HIZOUT,0,"FILTER_AFLAG_HIZOUT","Excessive smoothing correction (28)"},
  {FILTER_AFLAG_TOO_BRIGHT,0,"FILTER_AFLAG_TOO_BRIGHT","Too bright for calibrated magnitudes (29)"},
  {FILTER_AFLAG_LOW_ALTITUDE,0,"FILTER_AFLAG_LOW_ALTITUDE","Within 23.5 degrees of horizon (30)"},
  {FILTER_AFLAG_MULTIPLE_BLEND,0,"FILTER_AFLAG_MULTIPLE_BLEND","Multiple exposure blend (10)"},
  {FILTER_AFLAG_UNCERTAIN_DATE,0,"FILTER_AFLAG_UNCERTAIN_DATE","Uncertain time for extinction calculation (9)"},
  {FILTER_AFLAG_MULTIPLE_NONE,0,"FILTER_AFLAG_MULTIPLE_NONE","Unmatched object on a multiple exposure plate (8)"},
#if 1
  {FILTER_AFLAG_QUALITY,0,"FILTER_AFLAG_QUALITY","Plate-wide issue noted in Quality Flags (7)"},
#else
  {FILTER_AFLAG_QUALITY,0,"FILTER_AFLAG_QUALITY","Multiple exposure, grating, spectra, Pickering Wedge, yellow or red plate (7)"},
#endif
  {FILTER_AFLAG_BACKGROUND,0,"FILTER_AFLAG_BACKGROUND","High BACKGROUND (6)"}
};
int AflagsTableSize = sizeof(AflagsTable)/sizeof(FLAGSENTRY);

FLAGSENTRY BflagsTable[] = {
  {FILTER_BFLAG_NEIGHBORS,0,"FILTER_BFLAG_NEIGHBORS","Has neighbors"},
  {FILTER_BFLAG_BLEND,0,"FILTER_BFLAG_BLEND","Blended with another"},
  {FILTER_BFLAG_SATURATED,0,"FILTER_BFLAG_SATURATED","Saturated"},
  {FILTER_BFLAG_BOUNDARY,0,"FILTER_BFLAG_BOUNDARY","Too close to boundary"},
  {FILTER_BFLAG_APERTURE_INCOMP,0,"FILTER_BFLAG_APERTURE_INCOMP","Aperture data incomplete"},
  {FILTER_BFLAG_ISOPHOT_INCOMP,0,"FILTER_BFLAG_ISOPHOT_INCOMP","Isophotoal data incomplete"},
  {FILTER_BFLAG_DEBLEND_OVF,0,"FILTER_BFLAG_DEBLEND_OVF","Memory overflow during deblending"},
  {FILTER_BFLAG_EXTRACT_OVF,0,"FILTER_BFLAG_EXTRACT_OVF","Memory overflow during extraction"},
  {FILTER_BFLAG_ADJUST_BLEND,0,"FILTER_BFLAG_ADJUST_BLEND","Magnitude adjusted for blended star"},
  {FILTER_BFLAG_BINDRAD,0,"FILTER_BFLAG_BINDRAD","Low astrometric error in a high astronomic error bin"},
  {FILTER_BFLAG_PSFSATURATED,0,"FILTER_BFLAG_PSFSATURATED","Image PSF is saturated."},
  {FILTER_BFLAG_MAGDEP_MAGCOR,0,"FILTER_BFLAG_MAGDEP_MAGCOR","Magnitude-dependent correction applied."},
#ifndef NEW_COLOR_CALIBRATION
  {GSC_MAGNITUDE_FLAG_BIT,0,  "JMAG substituted for JPGMAG","JMAG substituted for JPGMAG"},
  {GSC_MAGNITUDE_FLAG_BIT+1,0,"BMAG substituted for JPGMAG","BMAG substituted for JPGMAG"},
  {GSC_MAGNITUDE_FLAG_BIT+2,0,"RMAG substituted for FPGMAG","RMAG substituted for FPGMAG"},
  {GSC_MAGNITUDE_FLAG_BIT+3,0,"VMAG substituted for FPGMAG","VMAG substituted for FPGMAG"},
#endif /* NEW_COLOR_CALIBRATION */
  {CAL_FLAG_SHIFT+PIPELINE_FLAG_BIT,0,"PIPELINE_FLAG","Pipeline detected no errors"},
  {CAL_FLAG_SHIFT+LOWESS_CAL_FLAG_BIT,0,"LOWESS_CAL_FLAG","Calibrated with a lowess fit"},
  {CAL_FLAG_SHIFT+LOCAL_CAL_FLAG_BIT,0,"LOCAL_CAL_FLAG","Locally calibrated"},
  {CAL_FLAG_SHIFT+EXTINCTION_FLAG_BIT,0,"EXTINCTION_FLAG","Corrected for extinction"},
  {CAL_FLAG_SHIFT+TOO_BRIGHT_FLAG_BIT,0,"TOO_BRIGHT_FLAG","Too bright"},
  {CAL_FLAG_SHIFT+COLOR_VALID_BIT,0,"COLOR_VALID","Valid color correction"},
  {CAL_FLAG_SHIFT+COLOR_METROPOLIS_BIT,0,"COLOR_METROPOLIS","Metropolis color search suceeded"},
  {FILTER_BFLAG_KEPLER,0,"FILTER_BFLAG_KEPLER","Kepler Input Catalog"},
  {FILTER_BFLAG_LATEMATCH,0,"FILTER_BFLAG_LATEMATCH","Match at end of pipeline"},
  {FILTER_BFLAG_PMERROR,0,"FILTER_BFLAG_PMERROR","High proper motion uncertainty"},
  {FILTER_BFLAG_COLORTERM,0,"FILTER_BFLAG_COLORTERM","Colorterm outside 5 and 95 percentile blue plate range"},
  {FILTER_BFLAG_DRAD_ADJUST,0,"FILTER_BFLAG_DRAD_ADJUST","Position adjusted by bin median errors"},
  {FILTER_BFLAG_BAD_JD,0,"FILTER_BFLAG_BAD_JD","Uncertain Julian Date"},
  {FILTER_BFLAG_PROPER_MOTION,0,"FILTER_BFLAG_PROPER_MOTION","Corrected for proper motion"}
};



int BflagsTableSize = sizeof(BflagsTable)/sizeof(FLAGSENTRY);

/* Exclusion zones come from /dasch/data/scanner/backup/2015_09_24/candidatedistribution.png */
EXCLUSION exclusionTable[] = {
  {241.0,250.0,  0.10,  0.75},
  {241.0,250.5,-11.10,-10.45},
  {231.5,240.5, -9.70, -9.30},
  {230.5,239.5,-20.90,-20.20},
  {214.5,224.0,-20.25,-19.40},
  {214.5,224.5,-31.23,-30.45},
  {185.9,196.0,-28.30,-26.40},
  {184.5,196.0,-38.10,-37.10},
};
int exclusionTableSize = sizeof(exclusionTable)/sizeof(EXCLUSION);

/* The separation below comes from direct measurement of a handle of plates.
 * The meanplussigma_drad will use the following ( for the "I" series)
 *
 * ./find_drad2.csh iplates.list
 * column -a -i find_drad2.db drad | compute 'drad=sqrt(dra^2+ddec^2)' | column drad | statstable
 *  Take the sum of Mean + RMS
 */

WEDGEENTRY wedgeTable[] = {
  {"mc",143.0,2.2},
  {"mc",147.0,2.2},
  {"mc",130.0,2.2},
#if 1
  {"b",266,2.0},
  {"b",260,2.0},
  {"b",270.0,8.9},
#endif
  {"dsy",116.4,3.2},
  {"dsy",115.2,3.2},
  {"dsy",91.2,3.2},
  {"i",287.0,10.8},
  {"i",276.0,10.8},
  {"",0.0,0.0}
};

int wedgeTableSize = sizeof(wedgeTable)/sizeof(WEDGEENTRY);

GSCBIN gscBin01  = {BIN_SIZE_01, DEC_BINS_01, TOTAL_GSC_BINS_01 ,0,NULL};
GSCBIN gscBin02  = {BIN_SIZE_02, DEC_BINS_02, TOTAL_GSC_BINS_02 ,0,NULL};
GSCBIN gscBin04  = {BIN_SIZE_04, DEC_BINS_04, TOTAL_GSC_BINS_04 ,0,NULL};
GSCBIN gscBin08  = {BIN_SIZE_08, DEC_BINS_08, TOTAL_GSC_BINS_08 ,0,NULL};
GSCBIN gscBin16  = {BIN_SIZE_16, DEC_BINS_16, TOTAL_GSC_BINS_16 ,0,NULL};
GSCBIN gscBin32  = {BIN_SIZE_32, DEC_BINS_32, TOTAL_GSC_BINS_32 ,0,NULL};

GSCBIN gscBin64  = {BIN_SIZE_64, DEC_BINS_64, TOTAL_GSC_BINS_64 ,0,NULL};
GSCBIN gscBin128 = {BIN_SIZE_128,DEC_BINS_128,TOTAL_GSC_BINS_128,0,NULL};

#ifdef   DEBUG_BUILDMOSAICLIST
static int debugPrint1 = 0;
static int debugPrint2 = 0;
#endif /* DEBUG_BUILDMOSAICLIST */
#ifdef DEBUG_NOTAVAILABLE
static int debugPrint3 = 0;
#endif /* DEBUG_NOTAVAILABLE */


double sqr(double a) {
  return(a*a);
}

int queryCount = 0;
unsigned int mysqlError = 0;

#if (RELEASE_LEVEL == RELEASE_LEVEL_ALL)
char ReleaseLevel[] = "daschprivate";
#else /* RELEASE_LEVEL */
char ReleaseLevel[] = "daschpublic";
#endif /* RELEASE_LEVEL */

SERIESENTRY seriesTableX[MAX_SERIES];
int seriesTableXInit = 0;


#ifdef USE_BINARRAY
/* The edgeDist parameter is no longer calculated */
int CalculateBin(int width,int height,double xLoc,double yLoc,double *pEdgeDist)
{
  int bin;
  int jx;
  int jy;
  int ix;
  int iy;
  int nx = X_DMAGBINS_NORMAL;
  int ny = Y_DMAGBINS_NORMAL;
  int local_bin;
#ifdef LARGEA_TEST
  static int printFlag = 0;
  if (printFlag == 0) {
    printf("ERROR: LARGEA_TEST is set\n");
    printFlag = 1;
  }
#endif /* LARGEA_TEST */
  *pEdgeDist = 0;
  if ((width*height) < SMALL_PLATE_PIXELS) {
    bin = 1;
    return(bin);
  }

  if ((nx != binarray_nx) ||
      (ny != binarray_ny)) {
    printf("ERROR: binarray.h created with different dimensions ix: %d %d iy: %d %d\n",binarray_nx,nx,binarray_ny,ny);
    exit(-1);
  }

  ix = (1.0*xLoc*nx)/(1.0*width);
  if (ix >= nx) {
    ix = nx-1;
  }
  if (ix < 0) {
    ix = 0;
  }
  iy = (1.0*yLoc*ny)/(1.0*height);
  if (iy >= ny) {
    iy = ny-1;
  }
  if (iy < 0) {
    iy = 0;
  }
  if (width < height) {
    jx = ix;
    jy = iy;
  } else {
    jx = iy;
    jy = ix;
  }

  local_bin = jx + (nx*jy);

  bin = binarray[local_bin];
  return(bin);

}
#else /* USE_BINARRAY */
/* The edgeDist parameter is the distance to the nearest border */
int CalculateBin(int width,int height,double xLoc,double yLoc,double *pEdgeDist)
{
  double dWidth = 1.0 * width;
  double dHeight = 1.0 * height;
  double dXLoc = 1.0 * xLoc;
  double dYLoc = 1.0 * yLoc;
  double xCenter = 0.5 + (dWidth  * 0.5);
  double yCenter = 0.5 + (dHeight * 0.5);
  /* The following variable names follow divide_annul9.m */
  double a0; /* Half shortest dimension */
  double b0; /* Half longest dimension */
  double a5;
  double dr = 0.0;
  double a6;
  double a7;
  double a8;
  double a4;
  double a3;
  double a2;
  double a1;
  double pi = 3.141592654;
  double dedge1 = 0.1;
  double plate_dist = sqrt(sqr(dXLoc-xCenter) + sqr(dYLoc - yCenter));
  int bin = 0;
  double minwidth = 0.5;
  double maxwidth = 0.5 + dWidth;
  double minheight = 0.5;
  double maxheight = 0.5 + dHeight;
  double edgeDist = maxheight;
  double edgeDist2 = maxheight;
  double edgeDist2a;
  double edgeDist2b;
  double edgeDist2c;
  double edgeDist2d;
#ifdef LARGEA_TEST
  static int printFlag = 0;
  if (printFlag == 0) {
    printf("ERROR: LARGEA_TEST is set\n");
    printFlag = 1;
  }
#endif /* LARGEA_TEST */
  *pEdgeDist = 0;

  if (width < height) {
    a0 = dWidth * 0.5;
    b0 = dHeight * 0.5;
  } else {
    b0 = dWidth * 0.5;
    a0 = dHeight * 0.5;
  }
#ifdef LOS_DEBUG_X
  a0 = 4.3268;
  b0 = 5.4316;
#endif /* LOS_DEBUG_X */

  a5 = 0.95 * 2.0 * a0 * sqrt(b0/(a0 * 2.0 * pi));
#ifdef LARGEA_TEST
  a8 = a5 * sqrt(4.0/5.0);
  a7 = a5 * sqrt(3.5/5.0);
  a6 = a5 * sqrt(3.0/5.0);
  a4 = a5 * sqrt(2.0/5.0);
  a3 = a5 * sqrt(1.5/5.0);
  a2 = a5 * sqrt(1.0/5.0);
  a1 = a5 * sqrt(0.5/5.0);
  a5 = a5 * sqrt(2.5/5.0);
#else /* LARGEA_TEST */
  dr = 1.05 * 2.0 * a0 * (((4.0 * b0)/a0) -1)/(30*sqrt((2.0 * pi * b0)/a0));
  a6 = a5 + (1.1*dr);
  a7 = a5 + (2.35 * dr);
  a8 = a5 + (5.0*dr);
  a4 = a5 * sqrt(4.0/5.0);
  a3 = a5 * sqrt(3.0/5.0);
  a2 = a5 * sqrt(2.0/5.0);
  a1 = a5 * sqrt(1.0/5.0);
#endif /* LARGEA_TEST */
#if 0
  printf("a0 %f, b0 %f, a1 %f, a2 %f a3 %f a4 %f a5 %f a6 %f a7 %f a8 %f dr %f\n",a0,b0,a1,a2,a3,a4,a5,a6,a7,a8,dr);
  exit(-1);
#endif

  edgeDist2a = (maxwidth - (dedge1 * a0)) - dXLoc;
  edgeDist2b = dXLoc - (minwidth + (dedge1 * a0));
  edgeDist2c = (maxheight - (dedge1 * a0)) - dYLoc;
  edgeDist2d = dYLoc - (minheight + (dedge1 * a0));

  if ((edgeDist2a >= 0) && (edgeDist2a < edgeDist2)) {
    edgeDist2 = edgeDist2a;
  }
  if ((edgeDist2b >= 0) && (edgeDist2b < edgeDist2)) {
    edgeDist2 = edgeDist2b;
  }
  if ((edgeDist2c >= 0) && (edgeDist2c < edgeDist2)) {
    edgeDist2 = edgeDist2c;
  }
  if ((edgeDist2d >= 0) && (edgeDist2d < edgeDist2)) {
    edgeDist2 = edgeDist2d;
  }


  if (plate_dist < a1) {
    bin = 1;
    edgeDist = a1-plate_dist;
  } else if (plate_dist < a2) {
    bin = 2;
    edgeDist = a2-plate_dist;
    if ((plate_dist-a1) < edgeDist) {
      edgeDist = plate_dist-a1;
    }
  } else if (plate_dist < a3) {
    bin = 3;

    edgeDist = a3-plate_dist;
    if ((plate_dist-a2) < edgeDist) {
      edgeDist = plate_dist-a2;
    }

  } else if (plate_dist < a4) {
    bin = 4;

    edgeDist = a4-plate_dist;
    if ((plate_dist-a3) < edgeDist) {
      edgeDist = plate_dist-a3;
    }

  } else if ((plate_dist < a5) &&
             (dXLoc < (maxwidth - (dedge1 * a0))) &&
             (dXLoc > (minwidth + (dedge1 * a0))) &&
             (dYLoc < (maxheight - (dedge1 * a0))) &&
             (dYLoc > (minheight + (dedge1 * a0)))) {
    bin = 5;

    edgeDist = a5-plate_dist;
    if ((plate_dist-a5) < edgeDist) {
      edgeDist = plate_dist-a4;
    }


  } else if ((plate_dist < a6) &&
             (dXLoc < (maxwidth - (dedge1 * a0))) &&
             (dXLoc > (minwidth + (dedge1 * a0))) &&
             (dYLoc < (maxheight - (dedge1 * a0))) &&
             (dYLoc > (minheight + (dedge1 * a0)))) {
    bin = 6;

    edgeDist = a6-plate_dist;
    if ((plate_dist-a5) < edgeDist) {
      edgeDist = plate_dist-a5;
    }
    if (edgeDist2 < edgeDist) {
      edgeDist = edgeDist2;
    }

  } else if ((plate_dist < a7) &&
             (dXLoc < (maxwidth - (dedge1 * a0))) &&
             (dXLoc > (minwidth + (dedge1 * a0))) &&
             (dYLoc < (maxheight - (dedge1 * a0))) &&
             (dYLoc > (minheight + (dedge1 * a0)))) {
    bin = 7;

    edgeDist = a7-plate_dist;
    if ((plate_dist-a6) < edgeDist) {
      edgeDist = plate_dist-a6;
    }
    if (edgeDist2 < edgeDist) {
      edgeDist = edgeDist2;
    }

  } else if ((plate_dist < a8) &&
             (dXLoc < (maxwidth - (dedge1 * a0))) &&
             (dXLoc > (minwidth + (dedge1 * a0))) &&
             (dYLoc < (maxheight - (dedge1 * a0))) &&
             (dYLoc > (minheight + (dedge1 * a0)))) {
    bin = 8;

    edgeDist = a8-plate_dist;
    if ((plate_dist-a7) < edgeDist) {
      edgeDist = plate_dist-a7;
    }
    if (edgeDist2 < edgeDist) {
      edgeDist = edgeDist2;
    }

  } else {
    bin = 9;
    if ((plate_dist - a8) >= 0) {
      edgeDist = plate_dist - a8;
    }
    edgeDist2a = -((maxwidth - (dedge1 * a0)) - dXLoc);
    edgeDist2b = -(dXLoc - (minwidth + (dedge1 * a0)));
    edgeDist2c = -((maxheight - (dedge1 * a0)) - dYLoc);
    edgeDist2d = -(dYLoc - (minheight + (dedge1 * a0)));

    if ((edgeDist2a >= 0.0) && (edgeDist2a < edgeDist2)) {
      edgeDist2 = edgeDist2a;
    }
    if ((edgeDist2b >= 0.0) && (edgeDist2b < edgeDist2)) {
      edgeDist2 = edgeDist2b;
    }
    if ((edgeDist2c >= 0.0) && (edgeDist2c < edgeDist2)) {
      edgeDist2 = edgeDist2c;
    }
    if ((edgeDist2d >= 0.0) && (edgeDist2d < edgeDist2)) {
      edgeDist2 = edgeDist2d;
    }

    if (edgeDist2 < edgeDist) {
      edgeDist = edgeDist2;
    }

  }

  if ((width*height) < SMALL_PLATE_PIXELS) {
    bin = 1;
  }

  if (edgeDist < 0) {
    printf("Edge distance %f is illegal\n");
  }

  *pEdgeDist = edgeDist;

  return(bin);

}
#endif /* USE_BINARRAY */


// Returns 1 on successful parse, 0 on failure
int ParseFilename2(char *filename,char* series,int *plateNumber)
{
  char *namePtr = filename;
  char *seriesPtr = series;
  char *slashPtr = NULL;
  int seriesLength = 0;
  char charVal;
  int nvals;


  series[0] = 0;
  *plateNumber = 0;

  /* First step to the last directory delimiter */
  while ((charVal = *namePtr) != 0) {
    if (charVal == '/') {
      slashPtr = namePtr;
    }
    namePtr++;
  }
  if (slashPtr == NULL) {
    namePtr = filename;
  } else {
    namePtr = slashPtr+1;
  }

  /* Search for the series */
  while ((charVal = *namePtr) != 0) {
    if ((charVal >= '0') && (charVal <= '9')) {
      break;
    }
    if (seriesLength >= MAX_SERIES_STRING) {
      return(0);
    }
    seriesLength++;
    *seriesPtr = charVal;
    seriesPtr++;
    namePtr++;
  }
  *seriesPtr = 0;
  if (seriesLength == 0) {
    return(0);
  }
  /* Now decode the integers */
  nvals = sscanf(namePtr,"%d_",plateNumber);
  if (nvals == 1) {
    return(1);
  }
  return(0);
}


// Returns 1 on successful parse, 0 on failure
int ParseFilename(char *filename,char* series,int *plateNumber,int *mosaicNumber,int *binning,int *rotation)
{
  char *namePtr = filename;
  char *seriesPtr = series;
  char *slashPtr = NULL;
  int seriesLength = 0;
  char charVal;
  int nvals;


  series[0] = 0;
  *plateNumber = 0;
  *mosaicNumber = 0;
  *rotation = 0;

  /* First step to the last directory delimiter */
  while ((charVal = *namePtr) != 0) {
    if (charVal == '/') {
      slashPtr = namePtr;
    }
    namePtr++;
  }
  if (slashPtr == NULL) {
    namePtr = filename;
  } else {
    namePtr = slashPtr+1;
  }

  /* Search for the series */
  while ((charVal = *namePtr) != 0) {
    if ((charVal >= '0') && (charVal <= '9')) {
      break;
    }
    if (seriesLength >= MAX_SERIES_STRING) {
      return(0);
    }
    seriesLength++;
    *seriesPtr = charVal;
    seriesPtr++;
    namePtr++;
  }
  *seriesPtr = 0;
  if (seriesLength == 0) {
    return(0);
  }
  /* Now decode the integers */
  nvals = sscanf(namePtr,"%d_%d_%dr%dww",plateNumber,mosaicNumber,binning,rotation);
  if (nvals == 4) {
    return(1);
  }
  if (nvals == 3) {
    *rotation = 0;
    return(1);
  }
  return(0);
}


#ifdef GENERATE_ARRAY

int main(int argc,char *argv[])
{
#ifdef WRITE_FILES
  FILE *outFile[10];
  char fileName[100];
#endif /* WRITE_FILES */
  FILE *binHandle;
  FILE *includeHandle;
  int includeCount = 0;
  int width = 17412;
  int height = 22026;
  int x = 22;
  int y = 33;
  int binCounts[11];
  int newBinCounts[11];
  int newBinSum = 0;
  int bin;
  int xindex;
  int yindex;
  int increment = 1;
  int counter = 0;
  int counter2 = 0;
  double edgeDist;
  double minEdgeDist = 1.0 * height;
  double maxEdgeDist = 0.0;
  time_t curTime;
  time_t beginTime;
  int ix;
  int iy;
  int local_bin;
  int nx = X_DMAGBINS;
  int ny = Y_DMAGBINS;
  int local_bin_count[X_DMAGBINS*Y_DMAGBINS];
  int local_bin_average[X_DMAGBINS*Y_DMAGBINS];
  double fraction;
  int intflag;
  double flag;

  memset(local_bin_count,0,sizeof(local_bin_count));
  memset(local_bin_average,0,sizeof(local_bin_average));

#if USE_BINARRAY
  printf("ERROR: USE_BINARRAY is set\n");
#endif /* USE_BINARRAY */

  time(&beginTime);

  memset(binCounts,0,sizeof(binCounts));
#ifdef WRITE_FILES
  for (bin = 1; bin <= 9; bin++) {
    sprintf(fileName,"binpoints_%d.tmp",bin);
    outFile[bin] = fopen(fileName,"wt");
  }
#endif /* WRITE_FILES */
  printf("pipelineutils bin calculation of %s %s with increment %d for width %d and height %d\n",__DATE__,__TIME__,increment,width,height);

#if 0
  bin = CalculateBin(width,height,x,y,&edgeDist);
  printf("Calculate Bin width %d, height %d, x %d, y %d, bin %d dist %f\n",width,height,x,y,bin,edgeDist);
#endif
  for (xindex = 1; xindex <= width; xindex+= increment) {
    for (yindex = 1; yindex <= height; yindex+= increment) {
      bin = CalculateBin(width,height,xindex-0.5,yindex-0.5,&edgeDist);
      if (edgeDist < minEdgeDist) {
        minEdgeDist = edgeDist;
      }
      if (edgeDist > maxEdgeDist) {
        maxEdgeDist = edgeDist;
      }

      binCounts[bin]++;
      counter++;
      ix = (1.0*(xindex-0.5)*nx)/(1.0*width);
      if (ix == nx) {
        ix = nx-1;
      }
      iy = (1.0*(yindex-0.5)*ny)/(1.0*height);
      if (iy == ny) {
        iy = ny-1;
      }
      local_bin = ix + (nx*iy);
      local_bin_count[local_bin]++;
      local_bin_average[local_bin] += bin;


#ifdef WRITE_FILES
      fprintf(outFile[bin],"%d %d\n",xindex,yindex);
#endif /* WRITE_FILES */
#if 0
      printf("Calculate Bin width %d, height %d, x %d, y %d, bin %d local_bin %d edge %f\n",width,height,xindex,yindex,bin,local_bin,edgeDist);
#endif
    }
  }
  binHandle = fopen("pipelineutils.db","wt");
  fprintf(binHandle,"ix\tiy\tlocal_bin\tbin_count\tflag\tfraction\tspatial_bin\n");
  fprintf(binHandle,"--\t--\t---------\t---------\t----\t--------\t-----------\n");
  memset(newBinCounts,0,sizeof(newBinCounts));

  includeHandle = fopen("binarray.h","wt");


  for (local_bin = 0; local_bin < (X_DMAGBINS*Y_DMAGBINS) ; local_bin++) {
    ix = local_bin % nx;
    iy = local_bin /nx;
    if (local_bin_count[local_bin] == 0) {
      flag = -1;
    } else {
      flag = (1.0*local_bin_average[local_bin])/(1.0*local_bin_count[local_bin]);
    }
    intflag = flag;
    fraction = flag-intflag;
    if (fraction > 0.6) {
      fraction = 1.0-fraction;
      intflag++;
      newBinSum++;
      newBinCounts[intflag]++;
    } else {
      newBinCounts[intflag]++;
      newBinSum++;
    }
#if 0
    printf("ix %2d iy %2d localbin %4d, spatial_bin %d newBinSum %4d\n",ix,iy,local_bin,intflag,newBinSum);
#endif
    if (local_bin == 0) {
      fprintf(includeHandle,"/* binarray.h \n");
      fprintf(includeHandle," * This file is generated automatically by pipelineutils when compiled and run with the GENERATE_ARRAY constant \n");
      fprintf(includeHandle," * Created on %s %s with increment %d for width %d and height %d\n",__DATE__,__TIME__,increment,width,height);
      fprintf(includeHandle," * ERROR: LARGEA_TEST is set \n");
      fprintf(includeHandle," */ \n");
      fprintf(includeHandle,"int binarray_nx = %d;\n",X_DMAGBINS);
      fprintf(includeHandle,"int binarray_ny = %d;\n",Y_DMAGBINS);

      fprintf(includeHandle,"int binarray[X_DMAGBINS*Y_DMAGBINS] = { \n    ");





    }

    if (((local_bin+1) % X_DMAGBINS) == 0) {

      if (local_bin == ((X_DMAGBINS*Y_DMAGBINS)-1)) {
        fprintf(includeHandle,"%d\n};\n",intflag);


      } else {
        fprintf(includeHandle,"%d,\n    ",intflag);
      }


    } else {

      fprintf(includeHandle,"%d,",intflag);

    }

    fprintf(binHandle,"%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%d\n",ix,iy,local_bin,local_bin_count[local_bin],flag,fraction,intflag);
  }

  fclose(binHandle);
  fclose(includeHandle);
  newBinSum = 0;
  for (bin = 1; bin <= 9; bin++) {
    newBinSum += newBinCounts[bin];

    printf("bin %d, count %9d, fraction %.5f newbinfraction %.5f ratio %.5f\n",bin,binCounts[bin],(1.0*binCounts[bin])/(1.0 * counter),(1.0*newBinCounts[bin]/2500.0),((1.0*newBinCounts[bin]/2500.0))/((1.0*binCounts[bin])/(1.0 * counter)));
    counter2 += binCounts[bin];
#ifdef WRITE_FILES
    fclose(outFile[bin]);
#endif /* WRITE_FILES */
  }

  time(&curTime);
  curTime = curTime - beginTime;
  printf("Total count %9d %9d seconds: %d min,max edge distance %f %f newBinSum %d\n",counter,counter2,curTime,minEdgeDist,maxEdgeDist,newBinSum);
  return(0);

}
#endif /* GENERATE_ARRAY */

int ExecuteQuery(void *connection,char *queryString)
{
  MYSQL *pConnection = connection;
  int res;
#if 0
  if (strlen(queryString) < 1000) {
    printf("%s: for %s\n",__FUNCTION__,queryString);
  } else {
    printf("%s: NOT PRINTING string of size  %d\n",__FUNCTION__,strlen(queryString));
  }
#endif
  res = mysql_query(pConnection,queryString);
  if (res || mysql_errno(pConnection)) {
    if (mysqlError == 0) {
      mysqlError = mysql_errno(pConnection);
    }
    printf("ERROR res is %d error is %d %s for %s \n",res,mysql_errno(pConnection),mysql_error(pConnection),queryString);
  }
  queryCount++;
  return(res);
}
unsigned int GetMysqlError()
{
  unsigned int retval = mysqlError;
  mysqlError = 0;
  return(retval);
}
void SetQueryCount(int count)
{
  queryCount = count;
}
int GetQueryCount()
{
  return(queryCount);
}
/* Returns 1 on success; 0 on failure */
int GetMosaicInfo(void *connection,
                  char *series,
                  int plateNumber,
                  int mosaicNumber,
                  int solutionNumber,
                  PMOSAIC pMosaic)
{
  MYSQL *pConnection = connection;
#if 1
  int mosaicListResult = 0;
  int mosaicInfoFlag = 1;
  int selectAllFlag = 0;
  int pendingFlag = 1;
  int verbose = 0;
  int numMosaicRecords = 0;
  PMOSAICLIST pMosaicTable = NULL;
  int totalSolution0Count = 0;
  int totalSolutionsCount = 0;



  mosaicListResult = BuildMosaicList(pConnection,
                                     series,     /* if NULL or zero length, consider all series */
                                     plateNumber,  /* if -1, consider all plates */
                                     mosaicNumber,  /* if -1, consider all mosaics */
                                     solutionNumber, /* if -1, consider all solutions */
                                     pendingFlag, /* if 1, then return unscanned plates */
                                     mosaicInfoFlag, /* if 1, then being called from GetMosaicInfo: no exposure or plates table info */
                                     selectAllFlag,  /* if 1, then return all mosaic records: good, deleted, and stale */
                                     verbose,         /* if 1, display warnings */
                                     &totalSolutionsCount,
                                     &totalSolution0Count,
                                     &numMosaicRecords,
                                     &pMosaicTable);
  if (mosaicListResult == 0) {
    return(0);
  }
  memcpy(pMosaic,pMosaicTable,sizeof(MOSAIC));
  if (pMosaicTable != NULL) {
    mosaicNumber = pMosaicTable->mosaicNumber;
    free(pMosaicTable);
    pMosaicTable = NULL;
  }

  if (numMosaicRecords != 1) {
    printf("ERROR: BuildMosaicList returned %d records  in pipelineutils.c line %d\n",numMosaicRecords,__LINE__);
    return(0);
  }
  return(1);

#else /* old code */
#ifdef   DEBUG_BUILDMOSAICLIST
  if (debugPrint1 == 0) {
    printf("ERROR: line %4d DEBUG_BUILDMOSAICLIST upgrade needed for %s\n",__LINE__,__FUNCTION__);
    debugPrint1 = 1;
  }



#endif /* DEBUG_BUILDMOSAICLIST */


  memset(pMosaic,0,sizeof(MOSAIC));
  strcpy(pMosaic->series,series);
  pMosaic->plateNumber = plateNumber;
  pMosaic->mosaicNumber = mosaicNumber;
  pMosaic->solutionNumber = solutionNumber;

  sprintf(queryString,"SELECT scanNumber,exposureNumber,binning+0,mosaicDate,naxis1,naxis2,ctype1,ctype2,crval1,crval2,crpix1,crpix2,cd1_1,cd1_2,cd2_1,cd2_2,FitWCS+0,rotation,WCSSource+0,mosaics.diskLocation,mosaicComment,transform+0,solutionNumber,scanDate,patternID from mosaics INNER JOIN scans using (series,plateNumber,scanNumber) where series = '%s' and plateNumber = %d and mosaicNumber = %d and solutionNumber = %d\n",
          series,plateNumber,mosaicNumber,solutionNumber);
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (gotAnswer) {
          printf("ERROR: two mosaic entries for %s%05d %2d",series,plateNumber,mosaicNumber);
        }
        gotAnswer = 1;
        pMosaic->solutionNumber = solutionNumber;
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pMosaic->scanNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for scanNumber",nvals);
          }
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%d",&pMosaic->exposureNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for exposureNumber",nvals);
          }
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pMosaic->binning);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for binning",nvals);
          }
        }
        if (sqlrow[3]) {
          strcpy(pMosaic->mosaicDate,sqlrow[3]);
        }
        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&pMosaic->naxis1);
        }
        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%d",&pMosaic->naxis2);
        }
        pMosaic->nx = XDmagBins(pMosaic->naxis1,pMosaic->naxis2);
        pMosaic->ny = YDmagBins(pMosaic->naxis1,pMosaic->naxis2);
        if (sqlrow[6]) {
          strcpy(pMosaic->ctype1,sqlrow[6]);
        }
        if (sqlrow[7]) {
          strcpy(pMosaic->ctype2,sqlrow[7]);
        }
        if (sqlrow[8]) {
          nvals = sscanf(sqlrow[8],"%lf",&pMosaic->crval1);
        }
        if (sqlrow[9]) {
          nvals = sscanf(sqlrow[9],"%lf",&pMosaic->crval2);
        }
        if (sqlrow[10]) {
          nvals = sscanf(sqlrow[10],"%lf",&pMosaic->crpix1);
        }
        if (sqlrow[11]) {
          nvals = sscanf(sqlrow[11],"%lf",&pMosaic->crpix2);
        }
        if (sqlrow[12]) {
          nvals = sscanf(sqlrow[12],"%lf",&pMosaic->cd1_1);
        }
        if (sqlrow[13]) {
          nvals = sscanf(sqlrow[13],"%lf",&pMosaic->cd1_2);
        }
        if (sqlrow[14]) {
          nvals = sscanf(sqlrow[14],"%lf",&pMosaic->cd2_1);
        }
        if (sqlrow[15]) {
          nvals = sscanf(sqlrow[15],"%lf",&pMosaic->cd2_2);
        }

        if (sqlrow[16]) {
          nvals = sscanf(sqlrow[16],"%d",&pMosaic->FitWCS);
        }
        if (sqlrow[17]) {
          nvals = sscanf(sqlrow[17],"%d",&pMosaic->rotation);
        }
        if (sqlrow[18]) {
          nvals = sscanf(sqlrow[18],"%d",&pMosaic->WCSSource);
        }
        if (sqlrow[19]) {
          nvals = sscanf(sqlrow[19],"%d",&pMosaic->diskLocation);
        }

        if (sqlrow[20]) {
          strcpy(pMosaic->mosaicComment,sqlrow[20]);
        }
        if (sqlrow[21]) {
          nvals = sscanf(sqlrow[21],"%d",&pMosaic->transform);
          if (nvals != 1) {
            pMosaic->transform = -1;
          }
        } else {
          pMosaic->transform = -1;
        }
        if (sqlrow[23]) {
          strcpy(pMosaic->scanDate,sqlrow[23]);
        } else {
          pMosaic->scanDate[0] = 0;
        }
        if (sqlrow[24]) {
          nvals = sscanf(sqlrow[24],"%d",&pMosaic->patternID);
          if (nvals != 1) {
            pMosaic->patternID = -1;
          }
        } else {
          pMosaic->patternID = -1;
        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
    }
  }
  return(gotAnswer);
#endif /* old code */
}


int GetExposureInfo(void *connection,
                    char *series,
                    int plateNumber,
                    int exposureNumber,
                    PEXPOSURE pExposure,
                    int logbookFlag)
{
  MYSQL *pConnection = (MYSQL *)connection;
  int nvals;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int gotAnswer = 0;
  char queryString[MAX_QUERY_STRING];

  memset(pExposure,0,sizeof(EXPOSURE));
  strcpy(pExposure->series,series);
  pExposure->plateNumber = plateNumber;
  pExposure->exposureNumber = exposureNumber;

  sprintf(queryString,"SELECT dRightAscension,dDeclination,date,centerSource,timeSource,timeAccuracy,exposure,notes from exposures where series = '%s' and plateNumber = %d and exposureNumber = %d\n",
          series,plateNumber,exposureNumber);
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (gotAnswer) {
          printf("ERROR: two exposure entries for %s%05d %2d\n",series,plateNumber,exposureNumber);
        }
        gotAnswer = 1;
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%lf",&pExposure->dRightAscension);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for dRightAscension\n", nvals);
            pExposure->dRightAscension = 999.0;
            pExposure->rightAscension[0] = 0;
          } else {
            ra2str(pExposure->rightAscension,MAX_RIGHTASCENSION_STRING,pExposure->dRightAscension,3);
          }
        } else {
          pExposure->dRightAscension = 999.0;
          pExposure->rightAscension[0] = 0;
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%lf",&pExposure->dDeclination);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for dDeclination\n", nvals);
            pExposure->dDeclination = 99.0;
            pExposure->declination[0] = 0;
          } else {
            dec2str(pExposure->declination,MAX_DECLINATION_STRING,pExposure->dDeclination,2);
          }
        } else {
          pExposure->dDeclination = 99.0;
          pExposure->declination[0] = 0;
        }
        if (sqlrow[2]) {
          strcpy(pExposure->date,sqlrow[2]);
        }
        if (sqlrow[3]) {
          strcpy(pExposure->centerSource,sqlrow[3]);
        }
        if (sqlrow[4]) {
          strcpy(pExposure->timeSource,sqlrow[4]);
        }

        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%lf",&pExposure->timeAccuracy);
          if (nvals != 1) {

            printf("ERROR: nvals is %d for timeAccuracy",nvals);
            pExposure->timeAccuracy = DEFAULT_TIME_ACCURACY;
          }
        } else {
          pExposure->timeAccuracy = DEFAULT_TIME_ACCURACY;
        }
        /* Fix of Jan 12, 2015.  If this is a card catalog entry with a only a date, set the time accuracy to one day */
        if ((pExposure->timeAccuracy == DEFAULT_TIME_ACCURACY) &&
            (strcmp(pExposure->timeSource,"Catalog") == 0) &&
            (strstr(pExposure->date,"T00:00:00.0") != NULL)) {
          pExposure->timeAccuracy = 1.0;
        }

        if (sqlrow[6]) {
          strcpy(pExposure->exposure,sqlrow[6]);
        }
        if (sqlrow[7]) {
          strcpy(pExposure->notes,sqlrow[7]);
        }

      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
    }
  }

  if (logbookFlag != 0) {
    sprintf(queryString,"SELECT dRightAscension,dDeclination from logbook where series = '%s' and plateNumber = %d and exposureNumber = %d\n",
            series,plateNumber,exposureNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%lf",&pExposure->dRightAscension);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for dRightAscension\n", nvals);
              pExposure->dRightAscension = 999.0;
              pExposure->rightAscension[0] = 0;
            } else {
              ra2str(pExposure->rightAscension,MAX_RIGHTASCENSION_STRING,pExposure->dRightAscension,3);
            }
          } else {
            pExposure->dRightAscension = 999.0;
            pExposure->rightAscension[0] = 0;
          }
          if (sqlrow[1]) {
            nvals = sscanf(sqlrow[1],"%lf",&pExposure->dDeclination);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for dDeclination\n", nvals);
              pExposure->dDeclination = 99.0;
              pExposure->declination[0] = 0;
            } else {
              dec2str(pExposure->declination,MAX_DECLINATION_STRING,pExposure->dDeclination,2);
            }
          } else {
            pExposure->dDeclination = 99.0;
            pExposure->declination[0] = 0;
          }

        }
        mysql_free_result(res_ptr);
      }
    }

  }
  return(gotAnswer) ;
}

void GetRangeTable(void *connection)
{
  MYSQL *pConnection = connection;
  PRANGES rangeTablePtr;
  int tmpNumRanges;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char queryString[MAX_QUERY_STRING];
  int curCount = 0;
  int nvals;


  if (numRanges != 0) {
    return;
  }

  sprintf(queryString,"SELECT series,firstPlate,lastPlate,locationId from ranges");
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      tmpNumRanges = (int)mysql_num_rows(res_ptr);
      if (tmpNumRanges == 0) {
        return;
      }
      rangeTablePtr = (PRANGES) malloc(sizeof(RANGES) * tmpNumRanges);
      pRangeTable = rangeTablePtr;
      if (rangeTablePtr) {
        memset(rangeTablePtr,0,sizeof(RANGES)*tmpNumRanges);
      } else {
        return;
      }
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        strcpy(rangeTablePtr->series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&rangeTablePtr->firstPlate);
        if (nvals != 1) {

          printf("ERROR: nvals is %d for firstPlate\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[2],"%d",&rangeTablePtr->lastPlate);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for lastPlate\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[3],"%d",&rangeTablePtr->locationId);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for locationId\n",nvals);
          continue;
        }
        rangeTablePtr++;
        curCount++;

      }
      numRanges = curCount;
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: %s: mysql_store_result failed\n",__FUNCTION__);
    }
  } else {
    printf("ERROR: %s: mysql_query failed\n",__FUNCTION__);
  }

}

int
GetLocationId(
  void *connection,
  char *series,
  int plateNumber
) {
  MYSQL *pConnection = connection;
  PRANGES rangePtr;
  int index, res;
  char queryString[MAX_QUERY_STRING];

  if (numRanges == 0) {
    GetRangeTable(pConnection);
  }

  for (index = 0; index < numRanges; index++) {
    rangePtr = &pRangeTable[index];

    if (
      strcmp(rangePtr->series,series) == 0 &&
      rangePtr->firstPlate <= plateNumber &&
      rangePtr->lastPlate >= plateNumber
    ) {
      return rangePtr->locationId;
    }
  }

  // If not in the `ranges` table, fall back to querying the plate record itself.

  sprintf(
    queryString,
    "SELECT locationId FROM plates WHERE series = '%s' AND plateNumber = %d;",
    series,
    plateNumber
  );

  res = ExecuteQuery(pConnection, queryString);

  if (!res) {
    MYSQL_RES *res_ptr = mysql_store_result(pConnection);

    if (res_ptr) {
      MYSQL_ROW sqlrow;

      while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
        if (sqlrow[0]) {
          int locid;

          if (sscanf(sqlrow[0], "%d", &locid) == 1) {
            return locid;
          }
        }
      }
    }
  }

  // Plate-level location ID is null or something; oh well.
  return -1;
}

void GetLocationTable(MYSQL *pConnection)
{
  PLOCATIONS locationTablePtr;
  int tmpNumLocations;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char queryString[MAX_QUERY_STRING];
  int curCount = 0;
  int nvals;


  if (numLocations != 0) {
    return;
  }
  sprintf(queryString,"SELECT locationId,wLongitude,latitude,elevation,name from locations \n");
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      tmpNumLocations = (int)mysql_num_rows(res_ptr);
      if (tmpNumLocations == 0) {
        return;
      }
      locationTablePtr = (PLOCATIONS) malloc(sizeof(LOCATIONS) * tmpNumLocations);
      pLocationTable = locationTablePtr;
      if (locationTablePtr) {
        memset(locationTablePtr,0,sizeof(LOCATIONS)*tmpNumLocations);
      } else {
        return;
      }
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        nvals = sscanf(sqlrow[0],"%d",&locationTablePtr->locationId);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for locationId\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[1],"%lf",&locationTablePtr->wLongitude);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for wLongitude\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[2],"%lf",&locationTablePtr->latitude);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for latitude\n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[3],"%lf",&locationTablePtr->elevation);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for elevation\n",nvals);
          continue;
        }
        if (sqlrow[4] != NULL) {
          strncpy(locationTablePtr->name,sqlrow[4],MAX_LOCATION_NAME);
        } else {
          printf("ERROR: location %d has no observatory name\n",locationTablePtr->locationId);
        }

        locationTablePtr++;
        curCount++;

      }
      numLocations = curCount;
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: %s: mysql_store_result failed\n",__FUNCTION__);
    }
  } else {
    printf("ERROR: %s: mysql_query failed\n",__FUNCTION__);
  }

}



int GetLongitude(void *connection,
                 int locationId,
                 double *longitude,
                 double *latitude,
                 double *elevation,
                 char *name,
                 int nameSize)
{
  PLOCATIONS pLocation;
  int index;
  MYSQL *pConnection  = (MYSQL *)connection;
  if (numLocations == 0) {
    GetLocationTable(pConnection);
  }

  *longitude = 0.0;
  *latitude = 0.0;
  *elevation = 0.0;
  for (index = 0; index < numLocations; index++) {
    pLocation = &pLocationTable[index];
    if (pLocation->locationId == locationId) {
      *longitude = pLocation->wLongitude;
      *latitude = pLocation->latitude;
      *elevation = pLocation->elevation;
      if (name != NULL) {
        strncpy(name,pLocation->name,nameSize);
      }
      return(1);
    }
  }
  return(0);

}
int dcmp(const void *a,
         const void *b)
{
  double *pa = (double *) a;
  double *pb = (double *) b;

  if ( *pa  > *pb ) return  1;
  if ( *pa == *pb ) return  0;

  return -1;
}


/* Sorting routine borrowed from statstable.c by John Roll */
/* zeroFlag is 0 for normal median and RMS
 *             1 to force the median to zero for a zero-based RMS calculation
 */


int CalcMedianAndRMS(int curObsCount,int minGoodStars,double *vector,double *med,double *rms,int doClip,double clipFactor,int zeroFlag)
{
  int index;
  int tmpObsCount;
  int pass;
  double magcal_local;
  double curMedian;
  double curRMS;
  double magcal_diff;
  double square_val;
  *med = 0;
  *rms = 0;
#if 0
  if (curObsCount == 17) {
    printf("curObsCount %d\n",curObsCount);
  }
#endif
  for (pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      /* Clip the existing data set at three sigma */
      tmpObsCount = 0;
      for (index = 0; index < curObsCount; index++) {
        magcal_local = vector[index];
        magcal_diff = curMedian - magcal_local;
        if (magcal_diff < 0.0) {
          magcal_diff = - magcal_diff;
        }
        if (magcal_diff <= (clipFactor*curRMS)) {
          vector[tmpObsCount] = magcal_local;
          tmpObsCount++;
        }
      }
#if 0
      if (curObsCount != tmpObsCount) {
        printf("curObsCount %d tmpObsCount %d\n",curObsCount,tmpObsCount);
      }
#endif
      curObsCount = tmpObsCount;
    }



    if ((curObsCount >= minGoodStars) && (curObsCount > 0)) {
      double aveMagnitude = 0.0;
      double sqrMagnitude = 0.0;
      for (index = 0; index < curObsCount; index++) {
        magcal_local = vector[index];
        aveMagnitude += magcal_local;
        sqrMagnitude += (magcal_local * magcal_local);
      }

      /* Median and Standard deviation borrowed from John Roll's statstable.c */
      double squ;
      double mean;
      double sum;
      double nval;
      int nred;
      nval = curObsCount;
      nred = nval;
      sum = aveMagnitude;
      squ = sqrMagnitude;
      mean = sum/nval;
      qsort(vector,nred,sizeof(double),dcmp);
      if ( nred == 1 ) {
        curMedian = vector[0];
      } else {
        if ( nred % 2 ) {
          curMedian =  vector[(nred/2)]; /* Bugfix of 10/24/2011 */
        } else {
          curMedian = (vector[(nred/2)-1] + vector[nred/2])/2;
        }
      }
      if (zeroFlag) {
        curMedian = 0.0;
        mean = 0.0;
      }
      if (nval == 1) {
        curRMS = 99.0;
      } else {
        square_val = (squ - 2 * mean * sum + nval * mean * mean);
        if (square_val < 0.0) {
          square_val = 0.0; /* assume round-off error */
        }
        curRMS = sqrt((square_val)/(nval-1));
      }
    } else {
      return(0);
    }
    if (doClip == 0) {
      break;
    }
  }
  *med = curMedian;
  *rms = curRMS;
  return(curObsCount);
}
/* Binary GSC catalog routines */

int GetDecBin(PGSCBIN pGscBin,double dec,char *outputname) {
  int dec_bin;
  dec_bin = (dec+90.0)/pGscBin->bin_size;
  if ((dec_bin < 0) || (dec_bin >= pGscBin->dec_bins)) {
    printf("ERROR: Illegal declination %f in GetDecBin for %s\n",dec,outputname);
    exit(-1);
  }
  return(dec_bin);
}
void InitBinIndex(PGSCBIN pGscBin)
{
  int dec_bin;
  int ra_bin;
  int ra_sum = 0;

  double declination;
  PBININDEX pBinIndex;
  if (pGscBin->pBinMasterIndex != NULL) {
    return;
  }
  pGscBin->pBinMasterIndex = (PBININDEX)calloc(pGscBin->dec_bins,sizeof(BININDEX));
  if (pGscBin->pBinMasterIndex == NULL) {
    printf("ERROR allocating bin index\n");
    exit(-1);
  }
  pGscBin->max_ra_bins = 0;
  for (dec_bin = 0; dec_bin < pGscBin->dec_bins; dec_bin++) {
    pBinIndex = &pGscBin->pBinMasterIndex[dec_bin];
    declination = (dec_bin * pGscBin->bin_size) - 90.0;
    ra_bin = (360.0/pGscBin->bin_size)*cos((declination + (pGscBin->bin_size/2))*DEGREES_TO_RAD);
    pBinIndex->declination = declination;
    pBinIndex->startBin = ra_sum;
    pBinIndex->numBins = ra_bin;
    if (pBinIndex->numBins > pGscBin->max_ra_bins) {
      pGscBin->max_ra_bins = pBinIndex->numBins;
    }
    ra_sum += ra_bin;
#if 0
    printf("dec_bin %3d, declination %11.7f, ra_bin %5d tot_bin %5d\n",
           dec_bin,declination,ra_bin,ra_sum);
#endif
  }
  if (ra_sum != pGscBin->total_gsc_bins) {
    printf("ERROR: Calculated number of GSC bins %d is not expected %d\n",ra_sum,pGscBin->total_gsc_bins);
    exit(-1);
  }

#if 0
  printf("Number of bins for GSC2.3.2 catalog: %d\n",ra_sum);
#endif
}
int GetGSCBin(PGSCBIN pGscBin,double ra,double dec,int *pDecBin,int *pRaBin, char *outputname)
{
  int dec_bin;
  PBININDEX pBinIndex;
  int ra_bin;
  int gsc_bin;
  *pDecBin = 0;
  *pRaBin = 0;

  InitBinIndex(pGscBin);
  dec_bin = GetDecBin(pGscBin,dec,outputname);
  pBinIndex = &pGscBin->pBinMasterIndex[dec_bin];
  while (ra < 0.0) {
    ra = ra+360.0;
  }
  while (ra >= 360.0) {
    ra = ra -360.0;
  }

  ra_bin = ((ra * pBinIndex->numBins) / 360.0);
  gsc_bin = ra_bin + pBinIndex->startBin;
#if 0
  printf("dec_bin %d ra_bin %d gsc_bin %d ",dec_bin,ra_bin,gsc_bin);
#endif
  *pDecBin = dec_bin;
  *pRaBin = ra_bin;
  return(gsc_bin);
}


PBININDEX GetSubbins(PGSCBIN pGscBin,int binNumber,int *pRaBin,int *pDecBin, char* outputname)
{
  PBININDEX pBinIndex;
  int min_bin = 0;
  int test_bin;
  int max_bin = pGscBin->dec_bins;
  InitBinIndex(pGscBin);
  if ((binNumber < 0) || (binNumber >= pGscBin->total_gsc_bins)) {
    printf("ERROR: GetSubBins request with an illegal bin number %d for %s\n",binNumber,outputname);
    exit(-1);
  }
  while (1) {
    test_bin = (min_bin+max_bin)/2;
    if ((test_bin < 0) ||
        (test_bin >= pGscBin->dec_bins)) {
      printf("ERROR: GetSubBins has an illegal bin number %d %d %d %d for %s\n",min_bin,max_bin,binNumber,test_bin,outputname);
      exit(-1);
      return(NULL);
    }
    pBinIndex =  &pGscBin->pBinMasterIndex[test_bin];
    if ((pBinIndex->startBin <= binNumber) &
        (binNumber < (pBinIndex->startBin + pBinIndex->numBins))) {
      *pRaBin = binNumber - pBinIndex->startBin;
      *pDecBin = test_bin;

      return(pBinIndex);
    } else if (pBinIndex->startBin > binNumber) {
      max_bin = test_bin;
    } else {
      min_bin = test_bin;
    }
  }
  printf("ERROR: GetSubBins has an illegal bin number %d for %s\n",binNumber,outputname);
  exit(-1);
  return(NULL);
}


/* Return the ra and dec of the bin center */
int GetBinCenter(PGSCBIN pGscBin,int gsc_bin_index,double *ra,double *dec,char *outputname) {
  double ctrRa;
  double ctrDec;
  int ra_bin;
  int dec_bin;
  PBININDEX pBinIndex;
  *ra = 0;
  *dec = 0;
  pBinIndex = GetSubbins(pGscBin,gsc_bin_index,&ra_bin,&dec_bin,outputname);
  if (pBinIndex == NULL) {
    return(-1);
  }
  ctrDec = pBinIndex->declination + pGscBin->bin_size/2;
  ctrRa = ((360.00 * ((1.0*ra_bin)+0.5))/pBinIndex->numBins);
  *ra = ctrRa;
  *dec = ctrDec;
  return(0);
}
/*
 *   Find least mean square fit of Y = (A*X) + B
 */

int LinearFit(int nPoints,double *X,double *Y,double * pA, double * pB, double * pMean, double * pStd, double * pResidual)
{
  int index;
  double denom;
  double A;
  double B;
  double mean;
  double std;
  double residual;
  double C = 0.0;
  double D = 0.0;
  double E = 0.0;
  double F = 0.0;
  double G = 0.0;
  *pStd = 0;
  *pMean = 0;
  *pA = 0;
  *pB = 0;
  *pResidual = 0;
  for (index = 0 ;index < nPoints; index++) {
    C += X[index];
    D += (X[index]) * (X[index]);
    E += Y[index];
    F += (Y[index]) * (Y[index]);
    G += (X[index]) * (Y[index]);
  }
  if (nPoints >= 2) {
    denom = (D * nPoints) - (C * C);
    if (denom != 0.0) {

      A = ((G * nPoints ) - (E * C))/denom;
      B = ((E * D)  - (G * C))/denom;
      mean = E / nPoints;
      std = sqrt((F/nPoints) - (mean * mean));
      residual = sqrt((F + (A*A*D) + (B * B * nPoints) - 2 *((A * G) + (B * E) - (A * B * C)))/nPoints);
#if 0
      printf("A: %g, B: %g mean %g, std %g residual %g\n",A,B,mean,std,residual);
#endif

      *pStd = std;
      *pMean = mean;
      *pA = A;
      *pB = B;
      *pResidual = residual;

      return(1);
    } else {
      return(0);
    }

  } else {
    return(0);
  }

}
/* MaxDradPixels was originally 5, but set to 2 on Aug 8, 2008
 * On Nov 25, 2008, the tolerance became a function of plate series
 * and was loosened by 1 pixel:
 *
 *  ac: 10 pixels, i: 5 pixels, all others: 3 pixels
 *  NOTE: filterblended uses the "none" series to get the smallest
 *  result.  The smallest result should be the default.
 *
 * On Jun 3, 2009, the i series tolerance was reduced to 3 pixels
 * in recognition that the astrometry errors in outer bins exist because
 * the stars are so highly distorted that sextractor thinks that they are
 * multiple images in a line.
 *
 * On November 29, 2010, these values were adjusted in accordance with SCAMP
 * rms matching values.
 *
 * On May 23, 2017, increase the "ax" series from 4.3 to 7.1 pixels
 */
double MaxDradPixels(char *series,int plateNumber)
{
  double result;
  if ((strcmp(series,"ac") == 0) ||
      (strcmp(series,"am") == 0) ||
      (strcmp(series,"ax") == 0) ||
      (strcmp(series,"ca") == 0)) {
    result = 7.1;
  } else if  ((strcmp(series,"ay") == 0) ||
              (strcmp(series,"x") == 0)) {
    result = 4.3;
  } else {
    result = 3.0;
  }

#if 0
  {
    static int printFlag = 1;
    if (printFlag) {
      printf("MaxDradPixels for %s%05d pixels %f\n",series,plateNumber,result);
      printFlag = 0;
    }
  }
#endif
  return(result);
}
/* Rotate x,y by an angle, which is positive in the counterclockwise direction */

void rotateFrame(double *coords,double angle)
{
  double x1 = coords[0];
  double y1 = coords[1];
  double cosval = cos(angle*DEGREES_TO_RAD);
  double sinval = sin(angle*DEGREES_TO_RAD);
  double x2 = (x1*cosval) - (y1*sinval);
  double y2 = (y1*cosval) + (x1*sinval);
  coords[0] = x2;
  coords[1] = y2;
  return;

}


/* Check blend takes the image dimensions and models them as rectangles.  If one
   of the verticies of rectangle 2 is inside rectangle 1, then we have a blend.

   NOTE: This algorithm will not work for certain pathological cases.  ellipticities
   are assumed to be comparable.

   All dimensions are in pixels

*/
int CheckBlend(double centerX1, /* center of image 1 */
               double centerY1,
               double theta1,   /* Angle of a axis with respect to the x axis for image 1 */
               double aLength1, /* Half length of image 1 */
               double bLength1,
               double centerX2, /* center of image 2 */
               double centerY2,
               double theta2,   /* Angle of b axis with respect to the x axis for image 2 */
               double aLength2, /* Half length of image 2 */
               double bLength2)
{

  /* Image A is the larger of the two images */
  /* Lower case "a" is the Sextractor "A" axis */
  double centerXA; /* center of image A */
  double centerYA;
  double thetaA;   /* Angle of a axis with respect to the x axis for image A */
  double aLengthA; /* Half length of image A */
  double bLengthA;
  double centerXB; /* center of image B */
  double centerYB;
  double thetaB;   /* Angle of b axis with respect to the x axis for image B */
  double aLengthB; /* Half length of image B */
  double bLengthB;
  double upperLeftB[2] = { 0.0, 0.0 };
  double upperRightB[2] = { 0.0, 0.0 };
  double lowerLeftB[2] = { 0.0, 0.0 };
  double lowerRightB[2] = { 0.0, 0.0 };
  double centerB[2];
  int result = 0;
  double xDist;
  double yDist;
  if (aLength1 > aLength2) {
    centerXA = centerX1;
    centerYA = centerY1;
    thetaA = theta1;
    aLengthA = aLength1;
    bLengthA = bLength1;
    centerXB = centerX2;
    centerYB = centerY2;
    thetaB = theta2;
    aLengthB = aLength2;
    bLengthB = bLength2;
  } else {
    centerXA = centerX2;
    centerYA = centerY2;
    thetaA = theta2;
    aLengthA = aLength2;
    bLengthA = bLength2;
    centerXB = centerX1;
    centerYB = centerY1;
    thetaB = theta1;
    aLengthB = aLength1;
    bLengthB = bLength1;

  }

  /* Start with a quick check */
  xDist = centerX1 - centerX2;
  if (xDist < 0.0) {
    xDist = -xDist;
  }
  if (xDist > (aLength1+aLength2)) {
    return(0);
  }
  yDist = centerY1 - centerY2;
  if (yDist < 0.0) {
    yDist = - yDist;
  }
  if (yDist > (aLength1+aLength2)) {
    return(0);
  }
#if 0
  if ((theta1 > 10.0) && (theta2 > 10.0)) {
    printf("theta1 %f, theta2 %f\n",theta1,theta2);
  }
#endif

  /* Get the corner points of the B rectangle
     in the frame of the B ellipse */
  upperLeftB[0]  = - aLengthB;
  upperLeftB[1]  = + bLengthB;

  if (aLengthB != 0) {
    upperRightB[0] = + aLengthB;
    upperRightB[1] = + bLengthB;

    lowerLeftB[0]  = - aLengthB;
    lowerLeftB[1]  = - bLengthB;
    lowerRightB[0] = + aLengthB;
    lowerRightB[1] = - bLengthB;
  }
  /* Transform 1: Covert the B rectangle to
     the frame of the mosaic */

  rotateFrame(upperLeftB,+thetaB);
  if (aLengthB != 0) {
    rotateFrame(upperRightB,+thetaB);
    rotateFrame(lowerLeftB,+thetaB);
    rotateFrame(lowerRightB,+thetaB);
  }
  /* Transform 2: Change the origin from the center of
     the B rectangle to the center of the A rectangle */

  upperLeftB[0] -=  (centerXA - centerXB);
  upperLeftB[1] -=  (centerYA - centerYB);
  if (aLengthB != 0) {
    upperRightB[0] -= (centerXA - centerXB);
    upperRightB[1] -= (centerYA - centerYB);

    lowerLeftB[0] -= (centerXA - centerXB);
    lowerLeftB[1] -= (centerYA - centerYB);
    lowerRightB[0] -= (centerXA - centerXB);
    lowerRightB[1] -= (centerYA - centerYB);
  }
  /* Transform 3: Rotate to the
     frame of the A rectangle */

  rotateFrame(upperLeftB,-thetaA);
  if (aLengthB != 0) {
    rotateFrame(upperRightB,-thetaA);
    rotateFrame(lowerLeftB,-thetaA);
    rotateFrame(lowerRightB,-thetaA);
  }
  /* Now see if any of the B corners lie within the A
     rectangle */


  if ((upperLeftB[0] >= -aLengthA) &&
      (upperLeftB[0] <= aLengthA) &&
      (upperLeftB[1] >= -bLengthA) &&
      (upperLeftB[1] <= bLengthA)) {
    result = 1;
  }
  if ((result == 0) && (aLengthB != 0)) {


    if ((upperRightB[0] >= -aLengthA) &&
        (upperRightB[0] <= aLengthA) &&
        (upperRightB[1] >= -bLengthA) &&
        (upperRightB[1] <= bLengthA)) {
      result = 1;
    } else if ((lowerLeftB[0] >= -aLengthA) &&
               (lowerLeftB[0] <= aLengthA) &&
               (lowerLeftB[1] >= -bLengthA) &&
               (lowerLeftB[1] <= bLengthA)) {
      result = 1;
    } else if ((lowerRightB[0] >= -aLengthA) &&
               (lowerRightB[0] <= aLengthA) &&
               (lowerRightB[1] >= -bLengthA) &&
               (lowerRightB[1] <= bLengthA)) {
      result = 1;
    }  else {
      /* Now see if the center of B lies within the A rectangle - added on Mar 3, 2015*/
      centerB[0] = (upperLeftB[0]+upperRightB[0]+lowerLeftB[0]+lowerRightB[0])/4.0;
      centerB[1] = (upperLeftB[1]+upperRightB[1]+lowerLeftB[1]+lowerRightB[1])/4.0;
      if ((centerB[0] >= -aLengthA) &&
          (centerB[0] <= aLengthA) &&
          (centerB[1] >= -bLengthA) &&
          (centerB[1] <= bLengthA)) {
        result = 1;
      }
    }
  }
  return(result);

}

/* This routine decodes the AFLAGS keywords and provides
 * a rejectReason bitmap and a plotting symbol.
 *
 * Returns:  1 if the point is accepted for plotting
 *           0 if the points is rejected for plotting
 *
 * rejectReason1 is the total summary of reasons for an object
 * rejectResaon2 is the summary of reasons meriting rejection for an object.
 *
 * If any bit AFLAGSMASK is set, then the bit is not used
 * to reject any point for plotting
 *
 *  For extract_lightcurves, AFLAGSMASK should be ((1<<FILTER_AFLAG_DEFECT)|FILTER_AMASK_BLEND)
 */
int DecodeAFLAGS(int AFLAGS,int AFLAGSMASK,int BFLAGS,int BFLAGSMASK,int *pRejectReason1,int *pRejectReason2,double *pSymbol)
{
  int rejectReason1 = 0;
  int rejectReason2 = 0;
  double badSymbol = 0;
  double goodSymbol = 0;
  int status = 0;
  int AFLAGScopy = AFLAGS;
  int BFLAGScopy = BFLAGS;

  AFLAGS &= ~AFLAGSMASK;
  BFLAGS &= ~BFLAGSMASK;

  /* Do not test blends for defect and drad */
  if ((AFLAGScopy & FILTER_AMASK_BLEND) != 0) {
    AFLAGS &= ~((1<<FILTER_AFLAG_DEFECT)|(1<<FILTER_AFLAG_DRAD));
  }

  /* Check for blend nomatch prior to checking for blends. These
     objects plot limiting_mag_local instead of magcal_local */

  if ((AFLAGScopy & (1<<FILTER_AFLAG_BLEND_NOMATCH)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_BLEND_NOMATCH;
    }
    rejectReason1 |= (1<<REJECT_REASON_BLEND_NOMATCH);
    if ((AFLAGS & (1<<FILTER_AFLAG_BLEND_NOMATCH)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_BLEND_NOMATCH;
      }
      rejectReason2 |= (1<<REJECT_REASON_BLEND_NOMATCH);
    }
  }
  if ((AFLAGScopy & (1<<FILTER_AFLAG_WEDGE)) != 0) {
    rejectReason1 |= (1<<REJECT_REASON_WEDGE);
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_WEDGE;
    }
    if ((AFLAGS & (1<<FILTER_AFLAG_WEDGE)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_WEDGE;
      }
      rejectReason2 |= (1<<REJECT_REASON_WEDGE);
    }
  }

  if ((AFLAGScopy & (1<<FILTER_AFLAG_DRAD)) != 0) {
    rejectReason1 |= (1<<REJECT_REASON_DRAD);
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_DRAD;
    }
    if ((AFLAGS & (1<<FILTER_AFLAG_DRAD)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_DRAD;
      }
      rejectReason2 |= (1<<REJECT_REASON_DRAD);
    }
  }
  if ((AFLAGScopy & (1<<FILTER_AFLAG_DEFECT)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_DEFECT;
    }
    rejectReason1 |= (1<<REJECT_REASON_DEFECT);
    if ((AFLAGS & (1<<FILTER_AFLAG_DEFECT)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_DEFECT;
      }
      rejectReason2 |= (1<<REJECT_REASON_DEFECT);
    }
  }

  if ((AFLAGScopy & (1<<FILTER_AFLAG_HIZOUT)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_HIZOUT;
    }
    rejectReason1 |= (1<<REJECT_REASON_HIZOUT);
    if ((AFLAGS & (1<<FILTER_AFLAG_HIZOUT)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_HIZOUT;
      }
      rejectReason2 |= (1<<REJECT_REASON_HIZOUT);
    }
  }

  if ((AFLAGScopy & (1<<FILTER_AFLAG_TOO_BRIGHT)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_TOO_BRIGHT;
    }
    rejectReason1 |= (1<<REJECT_REASON_TOO_BRIGHT);
    if ((AFLAGS & (1<<FILTER_AFLAG_TOO_BRIGHT)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_TOO_BRIGHT;
      }
      rejectReason2 |= (1<<REJECT_REASON_TOO_BRIGHT);
    }
  }

  if ((AFLAGScopy & (1<<FILTER_AFLAG_LOW_ALTITUDE)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_LOW_ALTITUDE;
    }
    rejectReason1 |= (1<<REJECT_REASON_LOW_ALTITUDE);
    if ((AFLAGS & (1<<FILTER_AFLAG_LOW_ALTITUDE)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_LOW_ALTITUDE;
      }
      rejectReason2 |= (1<<REJECT_REASON_LOW_ALTITUDE);
    }
  }


  if ((AFLAGScopy & (1<<FILTER_AFLAG_UNCERTAIN_DATE)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_LOW_ALTITUDE;
    }
    rejectReason1 |= (1<<REJECT_REASON_UNCERTAIN_DATE);
    if ((AFLAGS & (1<<FILTER_AFLAG_UNCERTAIN_DATE)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_LOW_ALTITUDE;
      }
      rejectReason2 |= (1<<REJECT_REASON_UNCERTAIN_DATE);
    }
  }

  if ((AFLAGScopy & (1<<FILTER_AFLAG_MULTIPLE_NONE)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_LOW_ALTITUDE;
    }
    rejectReason1 |= (1<<REJECT_REASON_MULTIPLE_NONE);
    if ((AFLAGS & (1<<FILTER_AFLAG_MULTIPLE_NONE)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_LOW_ALTITUDE;
      }
      rejectReason2 |= (1<<REJECT_REASON_MULTIPLE_NONE);
    }
  }


  if ((AFLAGScopy & ((1<<FILTER_AFLAG_CASEC) | (1<<FILTER_AFLAG_CASED))) != 0) {
    /* Always reject case c and d blended objects  */
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_BLEND;
    }
    rejectReason1 |= (1<<REJECT_REASON_COMPLEX_BLEND);
    if ((AFLAGS & ((1<<FILTER_AFLAG_CASEC) | (1<<FILTER_AFLAG_CASED))) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_BLEND;
      }
      rejectReason2 |= (1<<REJECT_REASON_COMPLEX_BLEND);
    }
  }
  if ((AFLAGScopy & (1<<FILTER_AFLAG_LOCAL_RMS)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_LOCAL_RMS;
    }
    rejectReason1 |= (1<<REJECT_REASON_LOCAL_RMS);
    if ((AFLAGS & (1<<FILTER_AFLAG_LOCAL_RMS)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_LOCAL_RMS;
      }
      rejectReason2 |= (1<<REJECT_REASON_LOCAL_RMS);
    }
  }
  if ((AFLAGScopy & (1<<FILTER_AFLAG_ISO_RMS)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_ISO_RMS;
    }
    rejectReason1 |= (1<<REJECT_REASON_ISO_RMS);
    if ((AFLAGS & (1<<FILTER_AFLAG_ISO_RMS)) != 0) {
      rejectReason2 |= (1<<REJECT_REASON_ISO_RMS);
      if (badSymbol == 0) {
        badSymbol = SYMBOL_ISO_RMS;
      }
    }
  }
  if ((AFLAGScopy & (1<<FILTER_AFLAG_DRADBIN)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_DRAD;
    }
    rejectReason1 |= (1<<REJECT_REASON_DRADBIN);
    if ((AFLAGS & (1<<FILTER_AFLAG_DRADBIN)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_DRAD;
      }
      rejectReason2 |= (1<<REJECT_REASON_DRADBIN);
    }
  }

  if ((AFLAGScopy & FILTER_AMASK_BLEND) != 0) {
    if ((goodSymbol == 0) ||
        (goodSymbol == SYMBOL_DEFECT) ||
        (goodSymbol == SYMBOL_DRAD)) {
      goodSymbol = SYMBOL_BLEND;
    }
    rejectReason1 |= (1<<REJECT_REASON_BLEND);
    if ((AFLAGS & FILTER_AMASK_BLEND) != 0) {
      badSymbol = SYMBOL_BLEND;
      if (badSymbol == 0) {
        rejectReason2 |= (1<<REJECT_REASON_BLEND);
      }
    }
  }

  if ((AFLAGScopy & FILTER_AFLAG_MULTIPLE_BLEND) != 0) {
    if ((goodSymbol == 0) ||
        (goodSymbol == SYMBOL_DEFECT) ||
        (goodSymbol == SYMBOL_BLEND) ||
        (goodSymbol == SYMBOL_DRAD)) {
      goodSymbol = SYMBOL_BLEND;
    }
    rejectReason1 |= (1<<REJECT_REASON_MULTIPLE_BLEND);
    if ((AFLAGS & FILTER_AFLAG_MULTIPLE_BLEND) != 0) {
      badSymbol = SYMBOL_BLEND;
      if (badSymbol == 0) {
        rejectReason2 |= (1<<REJECT_REASON_MULTIPLE_BLEND);
      }
    }
  }


  if ((AFLAGScopy & (1<<FILTER_AFLAG_QUALITY)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_QUALITY;
    }
    rejectReason1 |= (1<<REJECT_REASON_QUALITY);
    if ((AFLAGS & (1<<FILTER_AFLAG_QUALITY)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_QUALITY;
      }
      rejectReason2 |= (1<<REJECT_REASON_QUALITY);
    }
  }


  if ((BFLAGScopy & (1<<FILTER_BFLAG_PSFSATURATED)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_SATURATED;
    }
    rejectReason1 |= (1<<REJECT_REASON_SATURATED);
    if ((BFLAGS & (1<<FILTER_BFLAG_PSFSATURATED)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_SATURATED;
      }
      rejectReason2 |= (1<<REJECT_REASON_SATURATED);
    }
  }



  /* Warning!  We set NOMAGDEP only if the bit is clear */
  if ((BFLAGScopy & (1<<FILTER_BFLAG_MAGDEP_MAGCOR)) == 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_NOMAGDEP;
    }
    rejectReason1 |= (1<<REJECT_REASON_NOMAGDEP);

    if (((BFLAGSMASK  & (1<<FILTER_BFLAG_MAGDEP_MAGCOR)) == 0) &&
        ((BFLAGScopy & (1<<FILTER_BFLAG_MAGDEP_MAGCOR)) == 0)) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_NOMAGDEP;
      }
      rejectReason2 |= (1<<REJECT_REASON_NOMAGDEP);
    }
  }


  if ((AFLAGScopy & (1<<FILTER_AFLAG_BACKGROUND)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_BACKGROUND;
    }
    rejectReason1 |= (1<<REJECT_REASON_BACKGROUND);
    if ((AFLAGS & (1<<FILTER_AFLAG_BACKGROUND)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_BACKGROUND;
      }
      rejectReason2 |= (1<<REJECT_REASON_BACKGROUND);
    }
  }


  /* This one goes last because it it equivalent to a good symbol */
  if ((AFLAGScopy & (1<<FILTER_AFLAG_LIMITING_MAG)) != 0) {
    if (goodSymbol == 0) {
      goodSymbol = SYMBOL_LIMITING_MAG;
    }
    rejectReason1 |= (1<<REJECT_REASON_LIMITING_MAG);
    if ((AFLAGS & (1<<FILTER_AFLAG_LIMITING_MAG)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_LIMITING_MAG;
      }
      rejectReason2 |= (1<<REJECT_REASON_LIMITING_MAG);
    }
  }

  /* At this point, if there is no other condition, set in the symbol for a good
     object */
  if (goodSymbol == 0) {
    goodSymbol = SYMBOL_GOOD;
  }

  if ((AFLAGScopy & (1<<FILTER_AFLAG_BIN9)) != 0) {
    if ((goodSymbol == SYMBOL_GOOD) ||
        (goodSymbol == SYMBOL_LOCAL_RMS)) {
      goodSymbol += SYMBOL_BIN9;
    }
    rejectReason1 |= (1<<REJECT_REASON_BIN9);
    if ((AFLAGS & (1<<FILTER_AFLAG_BIN9)) != 0) {
      if (badSymbol == 0) {
        badSymbol = SYMBOL_GOOD+SYMBOL_BIN9;
      }
      rejectReason2 |= (1<<REJECT_REASON_BIN9);
    }
  }

  /* If there was no color correction for an otherwise good symbol, change
     it to a "c" */
  if ((goodSymbol == SYMBOL_GOOD) && ((BFLAGScopy & (1 << (CAL_FLAG_SHIFT+COLOR_VALID_BIT))) == 0)) {
    goodSymbol = SYMBOL_NOCOLOR;
  }
  /* If sextractor thinks that this star is a blend, but we do not, change the symbol to a "b" */
  if ((goodSymbol == SYMBOL_GOOD) && ((BFLAGScopy & FILTER_BLEND_MASK_SEXTRACTOR) != 0)) {
    goodSymbol = SYMBOL_NEIGHBORS;
  }



  /* Make sure a blend nomatch is not counted as a blend */
  if ((rejectReason1 & (1<<REJECT_REASON_BLEND_NOMATCH)) != 0) {
    rejectReason1 &= ~(1<<REJECT_REASON_BLEND);
  }
  if ((rejectReason2 & (1<<REJECT_REASON_BLEND_NOMATCH)) != 0) {
    rejectReason2 &= ~(1<<REJECT_REASON_BLEND);
  }

  if (badSymbol == 0) {
    status = 1;
    *pSymbol = goodSymbol;
  } else {
    status = 0;
    *pSymbol = badSymbol;
  }
  *pRejectReason1 = rejectReason1;
  *pRejectReason2 = rejectReason2;
  return(status);
}


void CalcMemory(int factor,signed long long newMemory,signed long long *curMemory,signed long long *maxMemory)
{
  if (factor > 0) {
    *curMemory += newMemory;
  } else {
    *curMemory -= newMemory;
  }
  if (*curMemory > *maxMemory) {
    *maxMemory = *curMemory;
  }
  return;
}
void  ComputeRange(int ngood,
                   double *vector,
                   double *min_local,
                   double *min_local2,
                   double *max_local2,
                   double *max_local)

{
  *min_local = 0;
  *max_local = 0;
  *min_local2 = 0;
  *max_local2 = 0;
  if (ngood <= 1) {
    return;
  }
  qsort(vector,ngood,sizeof(double),dcmp);

  *min_local = vector[0];
  *min_local2 = vector[1];
  *max_local2 = vector[ngood-2];
  *max_local = vector[ngood-1];
  return;
}

/* NOTE: this routine uses duplicate code from getlocation.c which should be removed from getlocation.c  */
int GetMargins(void *connection,
               char *series,
               int plateNumber,
               int mosaicNumber,
               int *pLeftMargin,
               int *pRightMargin,
               int *pBottomMargin,
               int *pTopMargin)
{
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  MYSQL *pConnection = connection;
  char queryString[MAX_QUERY_STRING];
  MOSAIC mosaicInfo;
  PMOSAIC pMosaic = &mosaicInfo;
  int totalCount = 0;
  int nvals;
  int patternID;
  int diskLocation;
  int foundEntryScans = 0;
  int foundEntryPattern = 0;
  int leftMargin;
  int rightMargin;
  int topMargin;
  int bottomMargin;
  int leftMargin2;
  int rightMargin2;
  int topMargin2;
  int bottomMargin2;


  *pLeftMargin = 0;
  *pRightMargin = 0;
  *pTopMargin = 0;
  *pBottomMargin = 0;

  if (GetMosaicInfo(pConnection,
                    series,
                    plateNumber,
                    mosaicNumber,
                    0,
                    pMosaic) != 1) {
    printf("ERROR: GetMargins can not find mosaic %s%05d_%02d\n",series,plateNumber,mosaicNumber);
    return(-1);
  }

  sprintf(queryString,"SELECT patternID,diskLocation  FROM scans where series = '%s' and plateNumber = %d and scanNumber = %d;",series,plateNumber,pMosaic->scanNumber);
  if (strlen(queryString) > MAX_QUERY_STRING) {
    printf("ERROR: GetMargins MAX_QUERY_STRING exceeded %lu\n",strlen(queryString));
    return(-1);
  }

  res = mysql_query(pConnection,queryString);
  if (!res) {
#if 0
    printf("insert ID: %lu rows res %d\n",(unsigned long)mysql_affected_rows(pConnection),res);
#endif
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
#if 0
      printf("Retrieved %lu rows\n",(unsigned long)mysql_num_rows(res_ptr));
#endif
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        totalCount++;
        if (sqlrow[0] == NULL) {
          printf("ERROR: GetMargins  No patternID found for %s%05d_%02d\n",series,plateNumber,mosaicNumber);
          continue;
        }
        nvals = sscanf(sqlrow[0],"%d",&patternID);
        if (nvals != 1) {
          printf("ERROR: GetMargins  nvals is %d for patternID\n",nvals);
          continue;
        }
        if (sqlrow[1] == NULL) {
#if 0
          printf("ERROR: GetMargins  No tile diskLocation found for %s%05d_%02d\n",series,plateNumber,mosaicNumber);
#endif
          diskLocation = 99999;
        } else {
          nvals = sscanf(sqlrow[1],"%d",&diskLocation);
          if (nvals != 1) {
            printf("ERROR: GetMargins  nvals is %d for diskLocation\n",nvals);
            continue;
          }
        }
        foundEntryScans++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
    }

  }


  if (foundEntryScans != 1) {
    printf("ERROR: GetMargins  %d scan database entries found for %s%05d_%02d\n",foundEntryScans,series,plateNumber,mosaicNumber);
    return(-1);
  }
  /* Next get the margins from the pattern file */
  sprintf(queryString,"SELECT leftMargin,rightMargin,bottomMargin,topMargin FROM scanPatterns where patternID = %d;",patternID);
  if (strlen(queryString) > MAX_QUERY_STRING) {
    printf("ERROR: GetMargins  MAX_QUERY_STRING exceeded %lu\n",strlen(queryString));
    return(-1);
  }
  if (pMosaic->transform < 0) {
    printf("ERROR: GetMargins  mosaic transform is null for Tiles are on disk %d (%s/%05d_%02d)\n",diskLocation,series,plateNumber,pMosaic->scanNumber);
    return(-1);
  }

  res = mysql_query(pConnection,queryString);
  if (!res) {
#if 0
    printf("insert ID: %lu rows res %d\n",(unsigned long)mysql_affected_rows(pConnection),res);
#endif
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
#if 0
      printf("Retrieved %lu rows\n",(unsigned long)mysql_num_rows(res_ptr));
#endif
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        totalCount++;
        if (sqlrow[0] == NULL) {
          leftMargin = 0;
        } else {
          nvals = sscanf(sqlrow[0],"%d",&leftMargin);
          if (nvals != 1) {
            printf("ERROR: GetMargins  nvals is %d for leftMargin\n",nvals);
            leftMargin = 0;
          }
        }
        if (sqlrow[1] == NULL) {
          rightMargin = 0;
        } else {
          nvals = sscanf(sqlrow[1],"%d",&rightMargin);
          if (nvals != 1) {
            printf("ERROR: GetMargins  nvals is %d for rightMargin\n",nvals);
            rightMargin = 0;
          }
        }
        if (sqlrow[2] == NULL) {
          bottomMargin = 0;
        } else {
          nvals = sscanf(sqlrow[2],"%d",&bottomMargin);
          if (nvals != 1) {
            printf("ERROR: GetMargins  nvals is %d for bottomMargin\n",nvals);
            bottomMargin = 0;
          }
        }
        if (sqlrow[3] == NULL) {
          topMargin = 0;
        } else {
          nvals = sscanf(sqlrow[3],"%d",&topMargin);
          if (nvals != 1) {
            printf("ERROR: GetMargins  nvals is %d for topMargin\n",nvals);
            topMargin = 0;
          }
        }


        foundEntryPattern++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
    }

  }


  if (foundEntryPattern != 1) {
    printf("ERROR: GetMargins  %d pattern database entries found for %s%05d_%02d\n",foundEntryPattern,series,plateNumber,mosaicNumber);
    return(-1);
  }

  /* Perform any necessary transformations */
  if ((pMosaic->transform & TRANSFORM_FLIP) != 0) {
    topMargin2 = bottomMargin;
    bottomMargin2 = topMargin;
    topMargin = topMargin2;
    bottomMargin = bottomMargin2;
  }
  if ((pMosaic->transform & TRANSFORM_MIRROR) != 0) {
    leftMargin2 = rightMargin;
    rightMargin2 = leftMargin;
    leftMargin = leftMargin2;
    rightMargin = rightMargin2;
  }
  switch(pMosaic->rotation) {
  case 0:
    break;
  case 90:
    topMargin2    = rightMargin;
    leftMargin2   = topMargin;
    bottomMargin2 = leftMargin;
    rightMargin2  = bottomMargin;
    topMargin     = topMargin2;
    leftMargin    = leftMargin2;
    bottomMargin  = bottomMargin2;
    rightMargin   = rightMargin2;
    break;
  case 180:
    topMargin2    = bottomMargin;
    leftMargin2   = rightMargin;
    bottomMargin2 = topMargin;
    rightMargin2  = leftMargin;
    topMargin     = topMargin2;
    leftMargin    = leftMargin2;
    bottomMargin  = bottomMargin2;
    rightMargin   = rightMargin2;
    break;
  case 270:
    topMargin2    = leftMargin;
    leftMargin2   = bottomMargin;
    bottomMargin2 = rightMargin;
    rightMargin2  = topMargin;
    topMargin     = topMargin2;
    leftMargin    = leftMargin2;
    bottomMargin  = bottomMargin2;
    rightMargin   = rightMargin2;
    break;
  default:
    printf("ERROR: illegal rotation %d for %s%05d_%02d\n",pMosaic->rotation,series,plateNumber,mosaicNumber);
    return(-1);
    break;
  }



  *pLeftMargin = leftMargin;
  *pRightMargin = rightMargin;
  *pTopMargin = topMargin;
  *pBottomMargin = bottomMargin;

  return(0);

}
void ChargeTime(int entry,long long *pOldUsec,long long *deltaTable,long long *totalTable)
{
  struct timeval tv;
  struct timezone tz;
  int index;
  long long curUsec;
  long long oldUsec = *pOldUsec;

  gettimeofday(&tv,&tz);
  curUsec = (1000000LL * tv.tv_sec) + (tv.tv_usec);

  if (entry == RESET_ALL_ENTRY) {
    for (index = 0; index < MAX_CHARGE_ENTRY; index++) {
      deltaTable[index]  = 0;
      totalTable[index] = 0;
    }
  } else if (entry == RESET_DELTA_ENTRY) {
    for (index = 0; index < MAX_CHARGE_ENTRY; index++) {
      deltaTable[index]  = 0;
    }
  } else {
    if ((entry >= 0) && (entry < MAX_CHARGE_ENTRY)) {
      deltaTable[entry] += curUsec-oldUsec;
      totalTable[entry] += curUsec-oldUsec;


    } else {
      printf("ERROR: Charge time invalid entry %d\n",entry);
      exit(-1);
    }
  }

  *pOldUsec = curUsec;

}
int GetDASCHNumber(double ra,double dec,long long *pREFNumber,int fatal,int refType)
{
  int nvals;
  char tempREF[MAX_REF+1];
  char REFString[MAX_REF+1];
  long long REFNumber = 0;
  if ((ra < 0.0) || (ra >= 360.0)) {
    printf("ERROR: Illegal ra %f in GetDASCHNumber\n",ra);
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
  }
  if ((dec < -90.0) || (dec > 90.0)) {
    printf("ERROR: Illegal dec %f in GetDASCHNumber\n",dec);
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
  }
  switch (refType) {
  case REF_TYPE_DASCH:
    REFString[0] = '3';
    break;
  case REF_TYPE_APASS:
    REFString[0] = '4';
    break;
  default:
    printf("ERROR: Illegal refType %d in GetDASCHNumber\n",refType);
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
    break;
  }

  ra2str(tempREF,MAX_REF+1,ra,1);
  REFString[1] = tempREF[0];
  REFString[2] = tempREF[1];
  REFString[3] = tempREF[3];
  REFString[4] = tempREF[4];
  REFString[5] = tempREF[6];
  REFString[6] = tempREF[7];
  REFString[7] = tempREF[9];
  dec2str(tempREF,MAX_REF+1,dec,0);
  if (dec >= 0) {
    REFString[8] = '1';
  } else {
    REFString[8] = '2';
  }
  REFString[9] = tempREF[1];
  REFString[10] = tempREF[2];
  REFString[11] = tempREF[4];
  REFString[12] = tempREF[5];
  REFString[13] = tempREF[7];
  REFString[14] = tempREF[8];
  REFString[15] = 0;
  nvals = sscanf(REFString,"%lld",&REFNumber);
  if ((nvals != 1) || (REFNumber == 0)) {
    printf("ERROR: failed to decode REFString %s in GetDASCHNumber \n",REFString);
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
  }
  *pREFNumber = REFNumber;
  return(0);
}

int GetDASCHCoordinates(char *REF,double *ra,double *dec, int verbose,int refType) {
  long long REFNumber;
  int tmpRefType;
  char REF2[MAX_REF+1];
  char tempREF[MAX_REF+1];
  double signVal;
  *ra = 0;
  *dec = 0;

  /* First verify that this really is a DASCH reference */
  if ((GetREFNumber(REF,&REFNumber,&tmpRefType,0,verbose) != 0) ||
      (tmpRefType != refType)) {
    return(-1);
  }
  GetREF(REFNumber,REF2,0,1);

  if (REF2[15] == '+') {
    signVal = 1.0;
  } else if (REF2[15] == '-') {
    signVal = -1.0;
  } else {
    return(-1);
  }
  if ((REF2[7] != '0') ||
      (REF2[8] != '0') ||
      (REF2[9] != '0') ||
      (REF2[10] != '0') ||
      (REF2[11] != '0') ||
      (REF2[12] != '0') ||
      (REF2[14] != '0')) {

    tempREF[0] = REF2[7];
    tempREF[1] = REF2[8];
    tempREF[2] = ':';
    tempREF[3] = REF2[9];
    tempREF[4] = REF2[10];
    tempREF[5] = ':';
    tempREF[6] = REF2[11];
    tempREF[7] = REF2[12];
    tempREF[8] = '.';
    tempREF[9] = REF2[14];
    tempREF[10] = 0;
    *ra = str2ra(tempREF);
    if (*ra == 0.0) {
      return(-1);
    }
  }
  if ((REF2[16] != '0') ||
      (REF2[17] != '0') ||
      (REF2[18] != '0') ||
      (REF2[19] != '0') ||
      (REF2[20] != '0') ||
      (REF2[21] != '0')) {
    tempREF[0] = REF2[16];
    tempREF[1] = REF2[17];
    tempREF[2] = ':';
    tempREF[3] = REF2[18];
    tempREF[4] = REF2[19];
    tempREF[5] = ':';
    tempREF[6] = REF2[20];
    tempREF[7] = REF2[21];
    tempREF[8] = 0;


    *dec = str2dec(tempREF);
    if (*dec == 0.0) {
      return(-1);
    }
    *dec = *dec * signVal;
  }

  return(0);
}
int GetREFType(long long REFNumber) {
  char tempREF[MAX_REF+1];
  int refType;
  sprintf(tempREF,"%lld",REFNumber);
  if (tempREF[0] == '0') {
    refType = REF_TYPE_NONE;
  } else if (tempREF[0] == '1') {
    refType = REF_TYPE_GSC;
  } else if (tempREF[0] == '2') {
    refType = REF_TYPE_KEPLER;
  } else if (tempREF[0] == '3') {
    refType = REF_TYPE_DASCH;
  } else if (tempREF[0] == '4') {
    refType = REF_TYPE_APASS;
  } else if (tempREF[0] == '5') {
    refType = REF_TYPE_TYCHO2;
  } else if (tempREF[0] == '6') {
    refType = REF_TYPE_UCAC4;
  } else if (tempREF[0] == '7') {
    refType = REF_TYPE_GAIA1;
  } else if (tempREF[0] == '8') {
    refType = REF_TYPE_GAIA2;
  } else if (tempREF[0] == '9') {
    refType = REF_TYPE_ATLAS2;
  } else {
    printf("ERROR: Unrecognized REFNumber %lld\n",REFNumber);
    exit(-1);
  }
  return(refType);
}

int GetREFNumber(char *REF,long long *pREFNumber,int *pRefType,int fatal,int verbose)
{
  int nvals;
  char tempREF[2*MAX_REF];
  long long REFNumber = 0;
  char *pString;
  int refType = REF_TYPE_NONE;
  *pREFNumber = 0;
  *pRefType = refType;
  if (strlen(REF) >= MAX_REF) {
    if (verbose) {
      printf("ERROR: Reference number %s length %lu exceeds maximum string length %d\n",REF,strlen(REF),MAX_REF);
    }
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
  }
  if ((REF[1] == 'O') && (strcmp(REF,"NONE") == 0)) {
    return(0);
  }
  if (REF[0] == 'N') {
    refType = REF_TYPE_GSC;
    strcpy(&tempREF[1],REF);
    tempREF[1] = '1';
    tempREF[0] = '1';
    pString = &tempREF[2];
  } else if (REF[0] == 'S') {
    refType = REF_TYPE_GSC;
    strcpy(&tempREF[1],REF);
    tempREF[1] = '2';
    tempREF[0] = '1';
    pString = &tempREF[2];
  } else if (REF[0] == 'K') {
    refType = REF_TYPE_KEPLER;
    if ((REF[1] == 'I') && (REF[2] == 'C')) {
      strcpy(tempREF,&REF[2]);
    } else {
      strcpy(tempREF,REF);
    }
      tempREF[0] = '2';
    pString = &tempREF[1];
  } else if (REF[0] == 'D') {
    /* Note: APASS designations below use the identical format */
    refType = REF_TYPE_DASCH;
    if (strlen(REF) == 19) {
      /* This could be a form generated by the id table script.  Replace "DASCH" with "3" and convert to integer */
      if ((pString = strstr(REF,"DASCH")) != NULL) {
        strcpy(&tempREF[1],&REF[5]);
        pString = &tempREF[0];
        *pString = '3';
      } else {
        if (verbose) {
          printf("ERROR: String does not begin with DASCH %s\n",REF);
        }
        if (fatal) {
          exit(-1);
        } else {
          return(-1);
        }

      }

    } else {
      /* Look for the official format */
        strcpy(tempREF,REF);
        if (strlen(REF) != 22) {
        if (verbose) {
          printf("ERROR: Bad DASCH length %lu in %s\n",strlen(REF),REF);
        }
        if (fatal) {
          exit(-1);
        } else {
          return(-1);
        }
      }

      tempREF[7] = 0;
      if ((strcmp(tempREF,"DASCH J") != 0) &&
          (strcmp(tempREF,"DASCH_J") != 0)) {
        if (verbose) {
          printf("ERROR: Unknown DASCH prefix  in %s\n",REF);
        }
        if (fatal) {
          exit(-1);
        } else {
          return(-1);
        }
      }
      if ((REF[13] != '.') ||
          ((REF[15] != '+') &&
           (REF[15] != '-'))) {
        if (verbose) {
          printf("ERROR: decimal or sign missing in %s\n",REF);
        }
        if (fatal) {
          exit(-1);
        } else {
          return(-1);
        }
      }
      pString = tempREF;
      tempREF[0] = '3';
      tempREF[1] = REF[7];
      tempREF[2] = REF[8];
      tempREF[3] = REF[9];
      tempREF[4] = REF[10];
      tempREF[5] = REF[11];
      tempREF[6] = REF[12];
      tempREF[7] = REF[14];
      if (REF[15] == '+') {
        tempREF[8] = '1';
      } else {
        tempREF[8] = '2';
      }
      tempREF[9] = REF[16];
      tempREF[10] = REF[17];
      tempREF[11] = REF[18];
      tempREF[12] = REF[19];
      tempREF[13] = REF[20];
      tempREF[14] = REF[21];
      tempREF[15] = 0;
    }

  } else if ((REF[0] == 'A') &&
             (REF[1] == 'P')) {
    /* Note: DASCH designations above use the identical format */
    refType = REF_TYPE_APASS;
    if (strlen(REF) == 19) {
      /* This could be a form generated by the id table script.  Replace "APASS" with "4" and convert to integer */
      if ((pString = strstr(REF,"APASS")) != NULL) {
        strcpy(&tempREF[1],&REF[5]);
        pString = &tempREF[0];
        *pString = '4';
      } else {
        if (verbose) {
          printf("ERROR: String does not begin with APASS %s\n",REF);
        }
        if (fatal) {
          exit(-1);
        } else {
          return(-1);
        }

      }

    } else {
      /* Look for the official format */
      strcpy(tempREF,REF);
      if (strlen(REF) != 22) {
        if (verbose) {
          printf("ERROR: Bad APASS length %lu in %s\n",strlen(REF),REF);
        }
        if (fatal) {
          exit(-1);
        } else {
          return(-1);
        }
      }

      tempREF[7] = 0;
      if ((strcmp(tempREF,"APASS J") != 0) &&
          (strcmp(tempREF,"APASS_J") != 0)) {
        if (verbose) {
          printf("ERROR: Unknown APASS prefix  in %s\n",REF);
        }
        if (fatal) {
          exit(-1);
        } else {
          return(-1);
        }
      }
      if ((REF[13] != '.') ||
          ((REF[15] != '+') &&
           (REF[15] != '-'))) {
        if (verbose) {
          printf("ERROR: decimal or sign missing in %s\n",REF);
        }
        if (fatal) {
          exit(-1);
        } else {
          return(-1);
        }
      }
      pString = tempREF;
      tempREF[0] = '4';
      tempREF[1] = REF[7];
      tempREF[2] = REF[8];
      tempREF[3] = REF[9];
      tempREF[4] = REF[10];
      tempREF[5] = REF[11];
      tempREF[6] = REF[12];
      tempREF[7] = REF[14];
      if (REF[15] == '+') {
        tempREF[8] = '1';
      } else {
        tempREF[8] = '2';
      }
      tempREF[9] = REF[16];
      tempREF[10] = REF[17];
      tempREF[11] = REF[18];
      tempREF[12] = REF[19];
      tempREF[13] = REF[20];
      tempREF[14] = REF[21];
      tempREF[15] = 0;
    }

  }  else if (REF[0] == 'T') {
    refType = REF_TYPE_TYCHO2;
    strcpy(tempREF,REF);
    tempREF[0] = '5';
    pString = &tempREF[1];
  }  else if (REF[0] == 'U') {
    refType = REF_TYPE_UCAC4;
    strcpy(tempREF,REF);
    tempREF[0] = '6';
    pString = &tempREF[1];
  } else if ((REF[0] == 'G') &&
             (REF[1] == 'A') &&
             (REF[2] == 'I') &&
             (REF[3] == 'A') &&
             (REF[4] == '1') &&
             (REF[5] == '_')) {
    refType = REF_TYPE_GAIA1;
    tempREF[0] = '7';
    strcpy(&tempREF[1],&REF[6]);
    pString = &tempREF[1];
  } else if ((REF[0] == 'G') &&
             (REF[1] == 'A') &&
             (REF[2] == 'I') &&
             (REF[3] == 'A') &&
             (REF[4] == '2') &&
             (REF[5] == '_')) {
    refType = REF_TYPE_GAIA2;
    tempREF[0] = '8';
    strcpy(&tempREF[1],&REF[6]);
    pString = &tempREF[1];
  } else if ((REF[0] == 'A') &&
             (REF[1] == 'T') &&
             (REF[2] == 'L') &&
             (REF[3] == 'A') &&
             (REF[4] == 'S') &&
             (REF[5] == '2') &&
             (REF[6] == '_')) {
    refType = REF_TYPE_ATLAS2;
    tempREF[0] = '9';
    strcpy(&tempREF[1],&REF[7]);
    pString = &tempREF[0];
  } else {
    if (verbose) {
      printf("ERROR: Unknown REF designation %s in GetREFNumber\n",REF);
    }
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
  }
  /* Check integrity of base number */
  nvals = sscanf(pString,"%lld",&REFNumber);
  if ((nvals != 1) || (REFNumber == 0)) {
    if (verbose) {
      printf("ERROR: failed to decode REF %s in GetREFNumber\n",REF);
    }
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
  }
  if (refType == REF_TYPE_GAIA1) {
    if ((REFNumber % GAIA_MODULUS) != 0) {
      printf("ERROR GetREFNumber GAIA REF %s is not a multiple of GAIA_MODULUS %d\n",REF,GAIA_MODULUS);
      if (fatal) {
        exit(-1);
      } else {
        return(-1);
      }
    }
    REFNumber = REFNumber/GAIA_MODULUS;
    sprintf(&tempREF[1],"%lld",REFNumber);
  }
  if (refType == REF_TYPE_GAIA2) {
    if ((REFNumber % GAIA_MODULUS) != 0) {
      printf("ERROR GetREFNumber GAIA REF %s is not a multiple of GAIA_MODULUS %d\n",REF,GAIA_MODULUS);
      if (fatal) {
        exit(-1);
      } else {
        return(-1);
      }
    }
    REFNumber = REFNumber/GAIA_MODULUS;
    sprintf(&tempREF[1],"%lld",REFNumber);
  }
  if ((tempREF[0] == '3') &&
      ((REFNumber < 300000000000000L) ||
       (REFNumber >=  400000000000000L))) {
    if (verbose) {
      printf("ERROR: unexpected text in REF %s in GetREFNumber\n",REF);
    }
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }


  }
  if ((tempREF[0] == '4') &&
      ((REFNumber <   400000000000000L) ||
       (REFNumber >=  500000000000000L))) {
    if (verbose) {
      printf("ERROR: unexpected text in REF %s in GetREFNumber\n",REF);
    }
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
  }
  if ((tempREF[0] == '9') &&
      ((REFNumber <   9000000000L) ||
       (REFNumber >   9992637835L))) {
    if (verbose) {
      printf("ERROR: unexpected text in REF %s in GetREFNumber\n",REF);
    }
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }


  }

#if 1
  nvals = sscanf(tempREF,"%lld",&REFNumber);
  if ((nvals != 1) || (REFNumber == 0)) {
    if (verbose) {
      printf("ERROR: failed to decode REF (2) %s in GetREFNumber\n",REF);
    }
    if (fatal) {
      exit(-1);
    } else {
      return(-1);
    }
  }
#endif
  *pREFNumber = REFNumber;
  *pRefType = refType;
  return(0);

}
int GetREF(long long REFNumber,char *REF,int nospace,int fatal)
{
  int nvals;
  char tempREF[2*MAX_REF];
  long long REFNumber2;
  char tempREF2[2*MAX_REF];
  if (REFNumber == 0) {
    strcpy(REF,"NONE");
    return(0);
  }
  sprintf(tempREF,"%lld",REFNumber);
  if (tempREF[0] == '1') {
    if (tempREF[1] == '1') {

      tempREF[1] = 'N';
      strcpy(REF,&tempREF[1]);
      return(0);
    } else if (tempREF[1] == '2') {
      tempREF[1] = 'S';
      strcpy(REF,&tempREF[1]);
      return(0);
    }
  } else if (tempREF[0] == '2') {
    tempREF[0] = 'K';
    strcpy(REF,&tempREF[0]);
    return(0);
  } else if (tempREF[0] == '3') {
    /* DASCH uses the same format as APASS below */
    if (strlen(tempREF) != 15) {
      printf("ERROR:  GetREF DASCH string length is %lu for %lld\n",strlen(tempREF),REFNumber);
      if (fatal) {
        exit(-1);
      } else {
        strcpy(REF,"UNKNOWN");
        return(-1);
      }
    }
    strcpy(REF,"DASCH J");
#if 1
      REF[5] = '_';
#else
    if (nospace) {
      REF[5] = '_';
    }
#endif
    REF[7]  = tempREF[1];
    REF[8]  = tempREF[2];
    REF[9]  = tempREF[3];
    REF[10] = tempREF[4];
    REF[11] = tempREF[5];
    REF[12] = tempREF[6];
    REF[13] = '.';
    REF[14] = tempREF[7];
    if (tempREF[8] == '1') {
      REF[15] = '+';
    } else if (tempREF[8] == '2') {
      REF[15] = '-';
    } else {
      printf("ERROR: GetREF unable to decode DASCH sign %c in REFNumber %lld\n",tempREF[8],REFNumber);
      if (fatal) {
        exit(-1);
      } else {
        strcpy(REF,"UNKNOWN");
        return(-1);
      }
    }

    REF[16] = tempREF[9];
    REF[17] = tempREF[10];
    REF[18] = tempREF[11];
    REF[19] = tempREF[12];
    REF[20] = tempREF[13];
    REF[21] = tempREF[14];
    REF[22] = 0;
    return(0);
  } else if (tempREF[0] == '4') {
    /* APASS uses the same format as DASCH above */
    if (strlen(tempREF) != 15) {
      printf("ERROR:  GetREF DASCH string length is %lu for %lld\n",strlen(tempREF),REFNumber);
      if (fatal) {
        exit(-1);
      } else {
        strcpy(REF,"UNKNOWN");
        return(-1);
      }
    }
    strcpy(REF,"APASS J");
#if 1
      REF[5] = '_';

#else
    if (nospace) {
      REF[5] = '_';
    }
#endif
    REF[7]  = tempREF[1];
    REF[8]  = tempREF[2];
    REF[9]  = tempREF[3];
    REF[10] = tempREF[4];
    REF[11] = tempREF[5];
    REF[12] = tempREF[6];
    REF[13] = '.';
    REF[14] = tempREF[7];
    if (tempREF[8] == '1') {
      REF[15] = '+';
    } else if (tempREF[8] == '2') {
      REF[15] = '-';
    } else {
      printf("ERROR: GetREF unable to decode APASS sign %c in REFNumber %lld\n",tempREF[8],REFNumber);
      if (fatal) {
        exit(-1);
      } else {
        strcpy(REF,"UNKNOWN");
        return(-1);
      }
    }

    REF[16] = tempREF[9];
    REF[17] = tempREF[10];
    REF[18] = tempREF[11];
    REF[19] = tempREF[12];
    REF[20] = tempREF[13];
    REF[21] = tempREF[14];
    REF[22] = 0;
    return(0);
  } else if (tempREF[0] == '5') {
    tempREF[0] = 'T';
    strcpy(REF,&tempREF[0]);
    return(0);
  } else if (tempREF[0] == '6') {
    tempREF[0] = 'U';
    strcpy(REF,&tempREF[0]);
    return(0);
  } else if (tempREF[0] == '7') {
    strcpy(REF,"GAIA1_");
    /* Here we need to multiply the number by GAIA_MODULUS */
    nvals = sscanf(&tempREF[1],"%lld",&REFNumber2);
    if (nvals != 1) {
      printf("ERROR: GetREF can not decode number in %lld\n",REFNumber);
      if (fatal) {
        exit(-1);
      } else {
        strcpy(REF,"GAIAUNKNOWN");
        return(-1);
      }
    }
    REFNumber2 *= GAIA_MODULUS;
    sprintf(tempREF2,"%lld",REFNumber2);
    strcpy(&REF[6],tempREF2);
    return(0);
  } else if (tempREF[0] == '8') {
    strcpy(REF,"GAIA2_");
    /* Here we need to multiply the number by GAIA_MODULUS */
    nvals = sscanf(&tempREF[1],"%lld",&REFNumber2);
    if (nvals != 1) {
      printf("ERROR: GetREF can not decode number in %lld\n",REFNumber);
      if (fatal) {
        exit(-1);
      } else {
        strcpy(REF,"GAIAUNKNOWN");
        return(-1);
      }
    }
    REFNumber2 *= GAIA_MODULUS;
    sprintf(tempREF2,"%lld",REFNumber2);
    strcpy(&REF[6],tempREF2);
    return(0);
  } else if (tempREF[0] == '9') {
    strcpy(REF,"ATLAS2_");
    strcpy(&REF[7],&tempREF[1]);
    return(0);
  }

  printf("ERROR: GetREF unable to decode REFNumber %lld\n",REFNumber);

  if (fatal) {
    exit(-1);
  } else {
    strcpy(REF,"UNKNOWN");
    return(-1);
  }
}
/* This routine gets the solution date twice: first assuming that
   the solution number is correct and next assuming that the solution number
   is incorrect as it may be after all of the multiple exposures are resolved

   Note: changes to this routine may need to be copied in GetTimeAccuracy();
*/

int GetSolutionJulianDate(void *connection,
                          char *series,
                          int plateNumber,
                          int mosaicNumber,
                          int solutionNumber,
                          int *pExposureNumber,
                          int *pNumExposures,
                          double *pSingleJulianDate,
                          double *pSingleTolerance,
                          double *pAllJulianDate,
                          double *pAllTolerance)
{
  MYSQL *pConnection = (MYSQL *)connection;
  EXPOSURE exposure[MAX_CATALOG_EXPOSURES];
  int numExposures = 0;
  PEXPOSURE pExposure;
  double maxTolerance = 0;
  double minJulianDate;
  double maxJulianDate;
  double singleJulianDate = 0;
  double singleTolerance = 0;
  double allJulianDate;
  double allTolerance;
  MOSAIC mosaicInfo;
  PMOSAIC pMosaic = &mosaicInfo;
  int exposureNumber;


  *pNumExposures = 0;
  *pSingleJulianDate = 0;
  *pSingleTolerance = 0;
  *pAllJulianDate = 0;
  *pAllTolerance = 0;
  *pExposureNumber = -1;

  if (GetMosaicInfo(pConnection,series,plateNumber,mosaicNumber,solutionNumber,pMosaic) != 1) {
    printf("ERROR: GetSolutionJulianDate can not find mosaic %s%05d_%02d s%d\n",series,plateNumber,mosaicNumber,solutionNumber);
    return(-1);
  }
  *pExposureNumber = pMosaic->exposureNumber;

  while (1) {
    pExposure = &exposure[numExposures];
    if (GetExposureInfo(pConnection,series,plateNumber,numExposures,pExposure,0)) {
      numExposures++;
    } else {
      break;
    }
  }
  *pNumExposures = numExposures;
  if (numExposures == 0) {
    return(-1);
  }
  /* See if we know our exposure accuracy to within two hours */
  pExposure = &exposure[0];
  maxTolerance = pExposure->timeAccuracy;
  minJulianDate = fd2jd(pExposure->date);
  maxJulianDate = fd2jd(pExposure->date);
  if (pExposure->exposureNumber == pMosaic->exposureNumber) {
    singleJulianDate = fd2jd(pExposure->date);
    singleTolerance = pExposure->timeAccuracy;
  }
  for (exposureNumber = 1; exposureNumber < numExposures; exposureNumber++) {
    pExposure = &exposure[exposureNumber];
    if (pExposure->exposureNumber == pMosaic->exposureNumber) {
      singleJulianDate = fd2jd(pExposure->date);
      singleTolerance = pExposure->timeAccuracy;
    }

    if (maxTolerance < pExposure->timeAccuracy) {
      maxTolerance = pExposure->timeAccuracy;
    }
    if (minJulianDate > fd2jd(pExposure->date)) {
      minJulianDate = fd2jd(pExposure->date);
    }
    if (maxJulianDate < fd2jd(pExposure->date)) {
      maxJulianDate = fd2jd(pExposure->date);
    }
  }
  if ((maxJulianDate - minJulianDate) > maxTolerance) {
    maxTolerance = maxJulianDate - minJulianDate;
  }
  allJulianDate = (minJulianDate + maxJulianDate)/2;
  allTolerance = maxTolerance;
  if (singleJulianDate == 0) {
    singleJulianDate = allJulianDate;
    singleTolerance = allTolerance;
  }

  *pSingleJulianDate = singleJulianDate;
  *pSingleTolerance = singleTolerance;
  *pAllJulianDate = allJulianDate;
  *pAllTolerance = allTolerance;

  return(0);

}
int SetMosaicFitWCS(void *connection,char* bitName,char *fileroot,int solutionNumber,int setFlag,int readFlag,int *pFitWCS)
{
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int rotation;
  int binning;
  MYSQL *pConnection = (MYSQL *)connection;
  int res = 0;
  char queryString[MAX_QUERY_STRING];
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int FitWCS;
  int tempFitWCS;
  int nvals;
  int gotRow = 0;
  PFITWCSBIT pFitWCSBit = FitWCSMasks;
  if (ParseFilename(fileroot,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
    printf("ERROR SetMosaicFitWCS Failed to parse filename %s as %s%05d_%02d_%02d rotation %d\n",
           fileroot,series,plateNumber,mosaicNumber,binning,rotation);
    return(-1);
  }

  while (pFitWCSBit->FitWCSMask != 0) {
    if (strcmp(bitName,pFitWCSBit->FitWCSDescr) == 0) {
      break;
    }
    pFitWCSBit++;
  }
  if (pFitWCSBit->FitWCSMask == 0) {
    printf("ERROR: SetMosaicFitWCS bit %s is undefined\n",bitName);
    return(-1);
  }

  if (solutionNumber < 0) {
    return(-1);
  }
  sprintf(queryString,"SELECT FitWCS+0 from mosaics where series = '%s' and plateNumber = %d and mosaicNumber = %d and solutionNumber = %d;",series,plateNumber,mosaicNumber,solutionNumber);
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        gotRow = 1;
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&FitWCS);
          if (nvals != 1) {
            FitWCS = 0;
            printf("ERROR: SetMosaicFitWCS nvals is %d for FitWCS\n",nvals);
          }
        } else {
          FitWCS = 0;
        }
        *pFitWCS = FitWCS;
        if (readFlag == 0) {
          tempFitWCS = FitWCS;
          if (setFlag) {
            tempFitWCS |= pFitWCSBit->FitWCSMask;
          } else {
            tempFitWCS &= ~pFitWCSBit->FitWCSMask;
          }

          if (tempFitWCS != FitWCS) {
            /* We need to update this one */
            sprintf(queryString,"UPDATE mosaics SET FitWCS = %d where series = '%s' AND plateNumber = %d AND mosaicNumber = %d and solutionNumber = %d;",tempFitWCS,series,plateNumber,mosaicNumber,solutionNumber);
            res = ExecuteQuery(pConnection,queryString);


          }
        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: SetMosaicFitWCS mysql_store_result failed\n");
      return(-1);
    }
  }
  if (gotRow == 0) {
    printf("ERROR: SetMosaicFitWCS failed to find solution %d for %s\n",solutionNumber,fileroot);
    return(-1);
  }
  return(res);


}

int SetMosaicFitWCS2(
  char* bitName,
  char *fileroot,
  int solutionNumber,
  int setFlag,
  int readFlag,
  int *pFitWCS
) {
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  int result;

  dasch_init_scandb(pConnection);
  result = SetMosaicFitWCS(pConnection, bitName, fileroot, solutionNumber, setFlag, readFlag, pFitWCS);
  mysql_close(pConnection);
  return result;
}

int SetPlateQuality(void *connection,char* bitName,char *fileroot,int setFlag,int readFlag,int *pQuality)
{
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int rotation;
  int binning;
  MYSQL *pConnection = (MYSQL *)connection;
  int res = 0;
  char queryString[MAX_QUERY_STRING];
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int quality;
  int tempQuality;
  int nvals;
  int gotRow = 0;
  PQUALITYBIT pQualityBit = qualityMasks;
  if (ParseFilename(fileroot,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
    printf("ERROR SetPlateQuality Failed to parse filename %s as %s%05d_%02d_%02d rotation %d\n",
           fileroot,series,plateNumber,mosaicNumber,binning,rotation);
    return(-1);
  }

  while (pQualityBit->qualityMask != 0) {
    if (strcmp(bitName,pQualityBit->qualityDescr) == 0) {
      break;
    }
    pQualityBit++;
  }
  if (pQualityBit->qualityMask == 0) {
    printf("ERROR: SetPlateQuality bit %s is undefined\n",bitName);
    return(0);
  }

  sprintf(queryString,"SELECT quality+0 from plates where series = '%s' and plateNumber = %d;",series,plateNumber);
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        gotRow = 1;
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&quality);
          if (nvals != 1) {
            quality = 0;
            printf("ERROR: SetPlateQuality nvals is %d for quality\n",nvals);
          }
        } else {
          quality = 0;
        }
        *pQuality = quality;
        if (readFlag == 0) {
          tempQuality = quality;
          if (setFlag) {
            tempQuality |= pQualityBit->qualityMask;
          } else {
            tempQuality &= ~pQualityBit->qualityMask;
          }

          if (tempQuality != quality) {
            /* We need to update this one */
            sprintf(queryString,"UPDATE plates SET quality = %d where series = '%s' AND plateNumber = %d;",tempQuality,series,plateNumber);
            res = ExecuteQuery(pConnection,queryString);


          }
        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: SetPlateQuality mysql_store_result failed\n");

    }
  }
  if (gotRow == 0) {
    printf("ERROR: SetPlateQuality failed to find plate for %s\n",fileroot);
  }
  return(res);
}


int
SetPlateQuality2(char* bitName, char *fileroot, int setFlag, int readFlag, int *pQuality)
{
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  int result;

  dasch_init_scandb(pConnection);
  result = SetPlateQuality(pConnection, bitName, fileroot, setFlag, readFlag, pQuality);
  mysql_close(pConnection);
  return result;
}

/* Return the plate color, or PLATECOLOR_BLUE if unknown.
 *  Note: the only instances of L3 as of January 2010 were for Hydrogen alpha or 0.670 interference filters, so call them red
 *  Call L4 class plates red also, although no known instances currently exist.
 */
int GetPlateColor(void *connection,char *series,int plateNumber)
{
  MYSQL *pConnection = (MYSQL *)connection;
  char plateClass[MAX_CLASS_STRING+1];
  int result = PLATECOLOR_BLUE;
  int res = 0;
  char queryString[MAX_QUERY_STRING];
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int gotRow = 0;

  plateClass[0] = 0;
  if ((strcmp(series,"dnr") == 0) ||
      (strcmp(series,"dsr") == 0)) {
    return(PLATECOLOR_RED);
  }
  if ((strcmp(series,"dny") == 0) ||
      (strcmp(series,"dsy") == 0)) {
    return(PLATECOLOR_YELLOW);
  }
  sprintf(queryString,"SELECT class from plates where series = '%s' and plateNumber = %d;",series,plateNumber);
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        gotRow = 1;
        if (sqlrow[0]) {
          if (strlen(sqlrow[0]) < MAX_CLASS_STRING) {
            strcpy(plateClass,sqlrow[0]);
          } else {
            printf("ERROR: Plate class length %lu exceeds maximum %d for %s%05d %s\n",strlen(sqlrow[0]),MAX_CLASS_STRING,series,plateNumber,sqlrow[0]);
            exit(-1);
          }
        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: GetPlateColor mysql_store_result failed for %s%05d\n",series,plateNumber);
      exit(-1);
    }
  }
  if (gotRow == 0) {
    printf("ERROR: GetPlateColor failed to find plate %s%05d\n",series,plateNumber);
    exit(-1);
  }


  if ((strcmp(plateClass,"L2") == 0) ||
      (strcmp(plateClass,"L3") == 0) ||
      (strcmp(plateClass,"L4") == 0) ||
      (strcmp(plateClass,"M2") == 0) ||
      (strcmp(plateClass,"C2") == 0) ||
      (strcmp(plateClass,"S2") == 0) ||
      (strcmp(plateClass,"D2") == 0) ||
      (strcmp(plateClass,"d2") == 0) ||
      (strcmp(plateClass,"F2") == 0) ||
      (strcmp(plateClass,"M2S") == 0) ||
      (strcmp(plateClass,"L2S") == 0) ||
      (strcmp(plateClass,"MS2") == 0)) {
    return(PLATECOLOR_RED);
  }

  if ((strcmp(plateClass,"L1") == 0) ||
      (strcmp(plateClass,"M1") == 0) ||
      (strcmp(plateClass,"C1") == 0) ||
      (strcmp(plateClass,"S1") == 0) ||
      (strcmp(plateClass,"D1") == 0) ||
      (strcmp(plateClass,"d1") == 0) ||
      (strcmp(plateClass,"F1") == 0) ||
      (strcmp(plateClass,"FD1") == 0)) {
    return(PLATECOLOR_YELLOW);
  }


  return(result);
}
/* These extinction coefficents were proposed by Sumin Tang in
 *  her memorandum of 12/10/09 9:12 PM
 */
double GetExtinctionCoefficient(int plateColor,double elevation) {
  double extinctionCoefficient;
  switch(plateColor) {
  case PLATECOLOR_BLUE:
    extinctionCoefficient = 0.4 - (0.2*(elevation/2300));
    break;
  case PLATECOLOR_YELLOW:
    extinctionCoefficient = 0.25 - (0.11*(elevation/2300));
    break;
  case PLATECOLOR_RED:
    extinctionCoefficient = 0.2 - (0.1*(elevation/2300));
    break;
  default:
    printf("ERROR: GetExtinction has unknown color %d\n",plateColor);
    exit(-1);
    break;
  }
  if ((elevation < 0.0) ||
      (extinctionCoefficient < 0.0)) {
    printf("ERROR: GetExtinction has illegal elevation %f or extinctionCoefficient %f for plateColor %d\n",
           elevation,extinctionCoefficient,plateColor);
    exit(-1);
  }
  return(extinctionCoefficient);

}
/* If series is null, then use the seriesId to identify the series */
int GetSeriesInfo(void *connection,
                  char *series,
                  int seriesId,
                  PSERIES_HEADER pSeriesHeader)

{
  int nvals;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int gotAnswer = 0;
  char queryString[MAX_QUERY_STRING];
  MYSQL *pConnection  = (MYSQL *)connection;

  memset(pSeriesHeader,0,sizeof(SERIES_HEADER));
  if ((series != NULL) && (strlen(series) > 0)) {
    sprintf(queryString,"SELECT description,aperture,nominalPlateScale,seriesId,sequestered+0,series,orientation+0,fittedPlateScale from series where series = '%s';",series);
  } else {
    sprintf(queryString,"SELECT description,aperture,nominalPlateScale,seriesId,sequestered+0,series,orientation+0,fittedPlateScale from series where seriesId = %d;",seriesId);
  }
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (gotAnswer) {
          printf("ERROR: two series entries for %s",series);
        }
        gotAnswer = 1;
        if (sqlrow[0]) {
          strncpy(pSeriesHeader->description,sqlrow[0],MAX_DESCRIPTION_STRING);
        }

        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%lf",&pSeriesHeader->aperture);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for aperture",nvals);
          }
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%lf",&pSeriesHeader->nominalPlateScale);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for nominalPlateScale",nvals);
          }
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pSeriesHeader->seriesId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for seriesId",nvals);
          }
        }
        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&pSeriesHeader->sequestered);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for sequestered",nvals);
          }
        }
        if (sqlrow[5]) {
          strncpy(pSeriesHeader->series,sqlrow[5],MAX_SERIES_STRING);
        }
        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%d",&pSeriesHeader->orientation);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for orientation",nvals);
          }
        }
        if (sqlrow[7]) {
          nvals = sscanf(sqlrow[7],"%lf",&pSeriesHeader->fittedPlateScale);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for fittedPlateScale",nvals);
          }
        }

      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
    }
  }
  return(gotAnswer) ;
}
int GetPlateInfo(void *connection,
                 char *series,
                 int plateNumber,
                 PPLATE_ENTRY pPlateEntry)
{
  int nvals;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int gotAnswer = 0;
  char queryString[MAX_QUERY_STRING];
  MYSQL *pConnection  = (MYSQL *)connection;

  memset(pPlateEntry,0,sizeof(PLATE_ENTRY));
  strcpy(pPlateEntry->series,series);
  pPlateEntry->plateNumber = plateNumber;

  sprintf(queryString,"SELECT locationId,class,quality,floodJacket+0 from plates where series = '%s' and plateNumber = %d;\n",
          series,plateNumber);
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (gotAnswer) {
          printf("ERROR: two plate entries for %s%05d",series,plateNumber);
        }
        gotAnswer = 1;
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pPlateEntry->locationId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for locationId\n", nvals);
          }
        }
        if (sqlrow[1]) {
          strcpy(pPlateEntry->plateClass,sqlrow[1]);
        }
        if (sqlrow[2]) {
          strncpy(pPlateEntry->quality,sqlrow[2],MAX_QUALITY_STRING);
        } else {
          pPlateEntry->quality[0] = 0;
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pPlateEntry->floodJacket);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for floodJacket\n", nvals);
          }
        }


      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d",__LINE__);
    }
  }
  return(gotAnswer) ;
}
int GetCatalogNumber(char *catalogName)
{
  int catalogNumber;
  for (catalogNumber = 0; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
    if (strstr(catalogName,catalogText[catalogNumber])) {
      return(catalogNumber);
    }
  }
  printf("ERROR: catalog must be one of ");
  for (catalogNumber = 0; catalogNumber < MAX_CATALOG_NUMBER; catalogNumber++) {
    printf("%s ",catalogText[catalogNumber]);
  }
  printf("\n");
  return(-1);

}
/* Given a plate number and a series, select the best mosaic candidate.  This source must match scanread.c and Catalog.cpp */
int SelectBestMosaic(void *pConnection,char *series,int plateNumber,int* pMosaicNumber, int* pRotation)
{
  int mosaicNumber = -1;
  int rotation = -1;
  int mosaicListResult = 0;
  int solutionNumber = 0;
  int mosaicInfoFlag = 1;
  int selectAllFlag = 0;
  int pendingFlag = 1;
  int verbose = 0;
  int numMosaicRecords = 0;
  PMOSAICLIST pMosaicTable = NULL;
  int totalSolution0Count = 0;
  int totalSolutionsCount = 0;

  *pMosaicNumber = -1;
  *pRotation = -1;


#if 1 /* New mosaic selection code */




  mosaicListResult = BuildMosaicList(pConnection,
                                     series,     /* if NULL or zero length, consider all series */
                                     plateNumber,  /* if -1, consider all plates */
                                     mosaicNumber,  /* if -1, consider all mosaics */
                                     solutionNumber, /* if -1, consider all solutions */
                                     pendingFlag, /* if 1, then return unscanned plates */
                                     mosaicInfoFlag, /* if 1, then being called from GetMosaicInfo: no exposure or plates table info */
                                     selectAllFlag,  /* if 1, then return all mosaic records: good, deleted, and stale */
                                     verbose,         /* if 1, display warnings */
                                     &totalSolutionsCount,
                                     &totalSolution0Count,
                                     &numMosaicRecords,
                                     &pMosaicTable);
  if (mosaicListResult == 0) {
    printf("ERROR: BuildMosaicList failed in pipelineutils.c line %d\n",__LINE__);
    return(0);
  }

  if (pMosaicTable != NULL) {
    mosaicNumber = pMosaicTable->mosaicNumber;
    rotation = pMosaicTable->rotation;
    free(pMosaicTable);
    pMosaicTable = NULL;
  }

  if (numMosaicRecords != 1) {
    printf("ERROR: BuildMosaicList returned %d records  in pipelineutils.c line %d\n",numMosaicRecords,__LINE__);
  }




#else /* old code */

#ifdef   DEBUG_BUILDMOSAICLIST
  if (debugPrint2 == 0) {
    printf("ERROR: line %4d DEBUG_BUILDMOSAICLIST upgrade needed for %s\n",__LINE__,__FUNCTION__);
    debugPrint2 = 1;
  }

#endif /* DEBUG_BUILDMOSAICLIST */

  sprintf(queryString,"SELECT mosaicNumber,rotation,FitWCS+0 FROM mosaics  WHERE series = '%s' and plateNumber = %d and solutionNumber = 0 and FitWCS regexp('Selected') and ((mosaicComment IS NULL) or (not(mosaicComment regexp('deleted')))) ORDER BY exposureNumber,scanNumber,mosaicNumber; ",series,plateNumber);
  res = ExecuteQuery(pConnection,queryString);

  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        tmpMosaicNumber = -1;
        tmpRotation = -1;
        tmpFitWCS = -1;
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&tmpMosaicNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for plateNumber",nvals);
            continue;
          }
        } else {
          continue;
        }

        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%d",&tmpRotation);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for rotation",nvals);
            continue;
          }
        } else {
          tmpRotation = 0;
        }

        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&tmpFitWCS);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for FitWCS",nvals);
            continue;
          }
        } else {
          continue;
        }
        if (mosaicNumber < 0) {
          mosaicNumber = tmpMosaicNumber;
          rotation = tmpRotation;
          FitWCS = tmpFitWCS;

        } else {
          if  ((tmpFitWCS & (FITWCS_SELECTED|FITWCS_COMPLETED)) >= (FitWCS & (FITWCS_SELECTED|FITWCS_COMPLETED))) {
            mosaicNumber = tmpMosaicNumber;
            rotation = tmpRotation;
            FitWCS = tmpFitWCS;
          }


        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
    }
  }
#endif /* End of old selection code */
  if (mosaicNumber >= 0) {
    *pMosaicNumber = mosaicNumber;
    *pRotation = rotation;
    return(1);
  }

  return(0);


}
/* Note: The following routines must agree with /dasch/mysql/classes.txt and Catalog.cpp */

int IsMultipleExposure(char *series,char *plateClass)
{
  if (strcmp(series,"b") == 0) {
    if (strcmp(plateClass,"K") == 0) {
      return(1);
    }
  } else {
    if ((strcmp(plateClass,"M") == 0) ||
        (strcmp(plateClass,"M1") == 0) ||
        (strcmp(plateClass,"M2") == 0) ||
        (strcmp(plateClass,"C") == 0) ||
        (strcmp(plateClass,"P") == 0) ||
        (strcmp(plateClass,"T1") == 0) ||
        (strcmp(plateClass,"M2S") == 0) ||
        (strcmp(plateClass,"Pan") == 0) ||
        (strcmp(plateClass,"CS") == 0) ||
        (strcmp(plateClass,"MG") == 0) ||
        (strcmp(plateClass,"MS2") == 0) ||
        (strcmp(plateClass,"Ms") == 0) ||
        (strcmp(plateClass,"PG") == 0) ||
        (strcmp(plateClass,"SM") == 0) ||
        (strcmp(plateClass,"TM") == 0) ||
        (strcmp(plateClass,"PG") == 0) ||
        (strcmp(plateClass,"PG") == 0) ||
        (strcmp(plateClass,"PG") == 0)) {
      return(1);
    }
  }

  return(0);
}
int IsGratingExposure(char *series,char *plateClass)
{
  if ((strcmp(plateClass,"G") == 0) ||
      (strcmp(plateClass,"MG") == 0) ||
      (strcmp(plateClass,"PG") == 0)) {
    return(1);
  }

  return(0);
}
int IsColorFilterPlate(char *series,char *plateClass)
{
  if ((strcmp(series,"dnr") == 0) ||
      (strcmp(series,"dny") == 0) ||
      (strcmp(series,"dsr") == 0) ||
      (strcmp(series,"dsy") == 0)) {
    return(1);
  }
  if ((strcmp(plateClass,"L1") == 0) ||
      (strcmp(plateClass,"L2") == 0) ||
      (strcmp(plateClass,"L3") == 0) ||
      (strcmp(plateClass,"L4") == 0) ||
      (strcmp(plateClass,"M1") == 0) ||
      (strcmp(plateClass,"M2") == 0) ||
      (strcmp(plateClass,"C1") == 0) ||
      (strcmp(plateClass,"C2") == 0) ||
      (strcmp(plateClass,"S1") == 0) ||
      (strcmp(plateClass,"S2") == 0) ||
      (strcmp(plateClass,"D1") == 0) ||
      (strcmp(plateClass,"D2") == 0) ||
      (strcmp(plateClass,"d1") == 0) ||
      (strcmp(plateClass,"d2") == 0) ||
      (strcmp(plateClass,"F1") == 0) ||
      (strcmp(plateClass,"F2") == 0) ||
      (strcmp(plateClass,"FD1") == 0) ||
      (strcmp(plateClass,"M2S") == 0) ||
      (strcmp(plateClass,"L2S") == 0) ||
      (strcmp(plateClass,"MS2") == 0)) {
    return(1);
  }


  return(0);
}
int IsSpectraPlate(char *series,char *plateClass)
{
  if (strcmp(series,"b") == 0) {
    if ((strcmp(plateClass,"E") == 0) ||
        (strcmp(plateClass,"F") == 0) ||
        (strcmp(plateClass,"G") == 0) ||
        (strcmp(plateClass,"I") == 0) ||
        (strcmp(plateClass,"K") == 0)) {
      return(1);
    }
  } else {


    if ((strcmp(plateClass,"D") == 0) ||
        (strcmp(plateClass,"D1") == 0) ||
        (strcmp(plateClass,"D2") == 0) ||
        (strcmp(plateClass,"d") == 0) ||
        (strcmp(plateClass,"d1") == 0) ||
        (strcmp(plateClass,"d2") == 0) ||
        (strcmp(plateClass,"O") == 0) ||
        (strcmp(plateClass,"I") == 0) ||
        (strcmp(plateClass,"OS") == 0) ||
        (strcmp(plateClass,"U") == 0) ||
        (strcmp(plateClass,"FD1") == 0)) {
      return(1);
    }

  }
  return(0);
}


void UpdateQuality(void *connection,char *series,int plateNumber)
{
  /* Check for the following bits
     QUALITY_MULTIPLE    Multiple exposure      - selected by class, logbook entries, scan comments, search_close, and mosaic solutions
     QUALITY_GRATING     Grating exposure       - selected by class and scan comments
     QUALITY_COLOR       Color filter used      - selected from plate class or Damon plate series.
     QUALITY_SPECTRA     Spectra plate          - selected by class and scan comments

     Not checked here -
     QUALITY_WEDGE       Pickering Wedge plate  - set by filter_wedge.c
     QUALITY_COLORTERM   Fails colorterm limits - depends on spatial bin

  */
  MYSQL *pConnection  = (MYSQL *)connection;
  int gotRow = 0;
  int res = 0;
  char queryString[MAX_QUERY_STRING];
  int curQuality;
  int newQuality;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int nvals;
  char plateClass[MAX_CLASS_STRING];
  EXPOSURE tempExposure;
  PEXPOSURE pExposure = &tempExposure;
  int ignoreFields = QUALITY_WEDGE | QUALITY_COLORTERM;
  int mosaicNumber;
  int solutionNumber;
  int FitWCS;
  int oldMosaicNumber;
  int haveMultipleMask = 0;
  int maskCount;
  int maskIndex;
  PMULTMASK pMultipleMask = NULL;
  PMULTMASK pMultipleTable = NULL;

  sprintf(queryString,"SELECT quality+0,class from plates where series = '%s' and plateNumber = %d;",series,plateNumber);
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        gotRow = 1;
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&curQuality);
          if (nvals != 1) {
            curQuality = 0;
            printf("ERROR: UpdateQuality nvals is %d for quality\n",nvals);
          }
        } else {
          curQuality = 0;
        }
        if (sqlrow[1]) {
          if (strlen(sqlrow[1]) < MAX_CLASS_STRING) {
            strcpy(plateClass,sqlrow[1]);
          } else {
            printf("ERROR: UpdateQuality class string %s is too long for %s%05d\n",sqlrow[1],series,plateNumber);
            plateClass[0] = 0;
          }
        } else {
          plateClass[0] = 0;
        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: UpdateQuality mysql_store_result failed\n");
    }
  }
  if (gotRow == 0) {
    printf("ERROR: UpdateQuality failed to find plate for %s%05d\n",series,plateNumber);
    return;
  }
  /* Copy over fields that we will not be checking */
  newQuality = curQuality & ignoreFields;


  if (IsMultipleExposure(series,plateClass) != 0)  {
    newQuality |= QUALITY_MULTIPLE;
  } else if  (GetExposureInfo(pConnection,series,plateNumber,1,pExposure,0) == 1) {
    /* Check to see if there is a multiple exposure in the logbook */
    newQuality |= QUALITY_MULTIPLE;
  }
  if (IsGratingExposure(series,plateClass) != 0) {
    newQuality |= QUALITY_GRATING;
  }
  if (IsColorFilterPlate(series,plateClass) != 0) {
    newQuality |= QUALITY_COLOR;
  }
  if (IsSpectraPlate(series,plateClass) != 0) {
    newQuality |= QUALITY_SPECTRA;
  }


  if ((newQuality & QUALITY_SPECTRA) == 0) {
    sprintf(queryString,"SELECT scanComment from scans where scanComment is not null and scanComment regexp('spectra') and series = '%s' and plateNumber = %d",series,plateNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          newQuality |= QUALITY_SPECTRA;
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: UpdateQuality mysql_store_result failed\n");
      }
    }

  }

  if ((newQuality & QUALITY_GRATING) == 0) {
    sprintf(queryString,"SELECT scanComment from scans where scanComment is not null and scanComment regexp('grating') and series = '%s' and plateNumber = %d",series,plateNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          newQuality |= QUALITY_GRATING;
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: UpdateQuality mysql_store_result failed\n");
      }
    }

  }

  if ((newQuality & QUALITY_MULTIPLE) == 0) {
    sprintf(queryString,"SELECT scanComment from scans where scanComment is not null and (scanComment regexp('multiple') or scanComment regexp('double') or scanComment regexp('triple')) and series = '%s' and plateNumber = %d",series,plateNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          newQuality |= QUALITY_MULTIPLE;
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: UpdateQuality mysql_store_result failed\n");
      }
    }

  }


  if ((newQuality & QUALITY_MULTIPLE) == 0) {
    oldMosaicNumber = -1;
    sprintf(queryString,"SELECT mosaicNumber,solutionNumber,FitWCS+0 from mosaics where series = '%s' and plateNumber = %d and WCSSource = 'imWCS' and FitWCS regexp('Selected') and FitWCS regexp('Completed') and ((mosaicComment IS NULL) OR (mosaicComment NOT regexp('Deleted'))) order by mosaicNumber,solutionNumber;",series,plateNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {


          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&mosaicNumber);
            if (nvals != 1) {
              printf("ERROR: UpdateQuality nvals is %d for mosaicNumber\n",nvals);
              continue;
            }
          } else {
            continue;
          }
          if (sqlrow[1]) {
            nvals = sscanf(sqlrow[1],"%d",&solutionNumber);
            if (nvals != 1) {
              printf("ERROR: UpdateQuality nvals is %d for solutionNumber\n",nvals);
              continue;
            }
          } else {
            continue;
          }
          if (sqlrow[2]) {
            nvals = sscanf(sqlrow[2],"%d",&FitWCS);
            if (nvals != 1) {
              printf("ERROR: UpdateQuality nvals is %d for FitWCS\n",nvals);
              FitWCS = 0;
            }
          } else {
            FitWCS = 0;
          }
#if 0
          printf("%5s%05d %2d %2d\n",series,plateNumber,mosaicNumber,solutionNumber);
#endif
          if (mosaicNumber == oldMosaicNumber) {
#if 0
          printf("STALE: %5s%05d %2d %2d\n",series,plateNumber,mosaicNumber,solutionNumber);
#endif
            continue;
          }
          if ((FitWCS & FITWCS_HAVEMULTMASK) != 0) {
            haveMultipleMask = 1;
          }
          if (((FitWCS & FITWCS_MULTFAILEDSEPARATION) != 0) ||
              ((FitWCS & FITWCS_MULTSUCCEEDED) != 0)) {
            newQuality |= QUALITY_MULTIPLE;
          }
          if ((FitWCS & FITWCS_MULTFAILEDASTROMETRY) != 0) {
            /* Further solutions are stale for this mosaic */
            oldMosaicNumber = mosaicNumber;
          }
        }

        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: UpdateQuality mysql_store_result failed\n");
      }
    }

  }


  if (((newQuality & QUALITY_MULTIPLE) == 0) &&
      (haveMultipleMask == 1)) {
    GetMultipleMask(pConnection,series,plateNumber,-1,&maskCount,&pMultipleTable);
    if (maskCount > 0) {
      for (maskIndex = 0; maskIndex < maskCount; maskIndex++) {
        pMultipleMask = &pMultipleTable[maskIndex];
        if (pMultipleMask->maskCenterFlag != 0) {
#if 0
          printf(".");
#endif
          continue;
        }
        if (pMultipleMask->maskCount > MAX_MASKS) {
#if 0
          printf("x");
#endif
          continue;
        }
        newQuality |= QUALITY_MULTIPLE;
        break;

      }
      if (pMultipleTable != NULL) {
        free(pMultipleTable);
        pMultipleTable= NULL;
      }
    }
  }
  if (curQuality != newQuality) {
    printf("UpdateQuality: curQuality 0x%02x newQuality 0x%02x for plate %5s%05d \n",curQuality,newQuality,series,plateNumber);
    sprintf(queryString,"UPDATE plates set quality = %d where series = '%s' and plateNumber = %d;",newQuality,series,plateNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (res) {
      exit(-1);
    }
  }

#if 0
  if ((curQuality & ~QUALITY_MULTIPLE) != (newQuality & ~QUALITY_MULTIPLE)) {
    /* Currently 42 errors from 1895 errors */

    printf("ERROR: UpdateQuality curQuality 0x%02x newQuality 0x%02x for plate %5s%05d \n",curQuality,newQuality,series,plateNumber);
  }
#endif

}

void CheckQuality(void *connection)
{

  MYSQL *pConnection  = (MYSQL *)connection;
  int res = 0;
  char queryString[MAX_QUERY_STRING];
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int nvals;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int counter = 0;
  int totalCount;


  sprintf(queryString,"SELECT series,plateNumber from plates order by series,plateNumber;");
  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {

      totalCount = (int)mysql_num_rows(res_ptr);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (sqlrow[0]) {
          if (strlen(sqlrow[0]) < MAX_SERIES_STRING) {
            strcpy(series,sqlrow[0]);
          } else {
            printf("ERROR: CheckQuality series string %s is too long for %s\n",sqlrow[0],series);
            strcpy(series,"UNK");
          }
        } else {
          strcpy(series,"UNK");
        }

        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%d",&plateNumber);
          if (nvals != 1) {
            plateNumber = 0;
            printf("ERROR: CheckQuality nvals is %d for quality\n",nvals);
          }
        } else {
          plateNumber = 0;
        }
        counter++;
        UpdateQuality(pConnection,series,plateNumber);
        if ((counter % 10000) == 0) {
          printf("counter %6d of %6d\n",counter,totalCount);
        }

      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: CheckQuality mysql_store_result failed\n");
    }
  }

}
int GetMultipleMask(void *connection,char *series,int plateNumber,int mosaicNumber,int *pMaskCount,PMULTMASK *ppMaskTable)
{
  MYSQL *pConnection  = (MYSQL *)connection;
  int res = 0;
  char queryString[MAX_QUERY_STRING];
  char tempString[MAX_QUERY_STRING];
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int nvals;
  int maskCount;
  int curMaskCount = 0;
  PMULTMASK pMaskTable = NULL;
  PMULTMASK pMultipleMask = NULL;

  *pMaskCount = 0;
  pMaskTable = *ppMaskTable;
  if (pMaskTable != NULL) {
    free(pMaskTable);
    pMaskTable = NULL;
    *ppMaskTable = NULL;
  }

  sprintf(queryString,"SELECT mosaicNumber,maskCount,maskIndex,maskXMin,maskXMax,maskYMin,maskYMax,maskArea,maxBinCount,maskCenterFlag,convolutionRadius,maskCenterDistance,maskAreaRatio,maskKernelType from masks where series = '%s' and plateNumber = %d",series,plateNumber);
  if (mosaicNumber >= 0) {
    sprintf(tempString," and mosaicNumber = %d",mosaicNumber);
    strcat(queryString,tempString);
  }
  strcat(queryString,";");


  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      maskCount =  (int)mysql_num_rows(res_ptr);
      if (maskCount == 0) {
        return(0);
      }
      pMaskTable = (PMULTMASK)calloc(maskCount,sizeof(MULTMASK));
      if (pMaskTable == NULL) {
        printf("ERROR: GetMultipleMask: failed to allocate %d MULTMASK entries\n",maskCount);
        exit(-1);
      }

      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (curMaskCount == maskCount) {
          printf("ERROR: GetMultipleMask: currentMaskCount %d maskCount %d for %s%5d\n",curMaskCount,maskCount,series,plateNumber);
          exit(-1);
        }

        pMultipleMask = &pMaskTable[curMaskCount];

        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&pMultipleMask->mosaicNumber);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for mosaicNumber %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL mosaicNumber %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1],"%d",&pMultipleMask->maskCount);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskCount %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskCount %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%d",&pMultipleMask->maskIndex);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskIndex %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskIndex %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%d",&pMultipleMask->maskXMin);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskXMin %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskXMin %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&pMultipleMask->maskXMax);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskXMax %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskXMax %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%d",&pMultipleMask->maskYMin);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskYMin %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskYMin %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[6]) {
          nvals = sscanf(sqlrow[6],"%d",&pMultipleMask->maskYMax);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskYMax %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskYMax %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[7]) {
          nvals = sscanf(sqlrow[7],"%d",&pMultipleMask->maskArea);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskArea %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskArea %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[8]) {
          nvals = sscanf(sqlrow[8],"%d",&pMultipleMask->maxBinCount);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maxBinCount %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maxBinCount %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[9]) {
          nvals = sscanf(sqlrow[9],"%d",&pMultipleMask->maskCenterFlag);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskCenterFlag %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskCenterFlag %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[10]) {
          nvals = sscanf(sqlrow[10],"%d",&pMultipleMask->convolutionRadius);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for convolutionRadius %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL convolutionRadius %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[11]) {
          nvals = sscanf(sqlrow[11],"%lf",&pMultipleMask->maskCenterDistance);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskCenterDistance %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskCenterDistance %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[12]) {
          nvals = sscanf(sqlrow[12],"%lf",&pMultipleMask->maskAreaRatio);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskAreaRatio %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskAreaRatio %s%5d \n",series,plateNumber);
            continue;
        }
        if (sqlrow[13]) {
           nvals = sscanf(sqlrow[13],"%d",&pMultipleMask->maskKernelType);
          if (nvals != 1) {
            printf("ERROR: GetMultipleMask: nvals is %d for maskKernelType %s%5d \n",nvals,series,plateNumber);
            continue;
          }
        } else {
            printf("ERROR: GetMultipleMask: NULL maskKernelType %s%5d \n",series,plateNumber);
            continue;
        }
        strcpy(pMultipleMask->series,series);
        pMultipleMask->plateNumber = plateNumber;
        curMaskCount++;


      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: GetMultpleMask mysql_store_result failed\n");
    }
  }
  if (curMaskCount > 0) {
    *pMaskCount = curMaskCount;
    *ppMaskTable = pMaskTable;
    return(1);
  }
  if (pMaskTable != NULL) {
    free(pMaskTable);
  }
  return(0);

}
/* NOTE: A change here must also be copied to /home/scanner/linux/apps/FitsWrapper/Coverage.cpp */
int CheckAuthorization(char *regionFlag,double ra, double dec,int *pReleaseField)
{
  int result = 0;
  char regionString[40];
  double lat;
  double lon;
  double ra0;
  double dec0;
  double dist;
  int releaseField = RELEASE_FIELD_OTHER;
  *pReleaseField = releaseField;

  if ((regionFlag != NULL) && (strcmp(regionFlag,"all") == 0)) {
    /* Authorized for everything */
#if 0
    printf("Authorized with %s\n",regionFlag);
#endif
    return(1);
  }
  strcpy(regionString,releaseFieldText[RELEASE_FIELD_OTHER]);


  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR1)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if (lat > 75.0) {
      releaseField = RELEASE_FIELD_DR1;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR1]);
      result = 1;
    }
  }


  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR2)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= 75.0) && (lat > 60.0)) {
      releaseField = RELEASE_FIELD_DR2;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR2]);
      result = 1;
    }
  }

  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR3)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= 60.0) && (lat > 45.0)) {
      releaseField = RELEASE_FIELD_DR3;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR3]);
      result = 1;
    }
  }

  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR4)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= 45.0) && (lat > 30.0)) {
      releaseField = RELEASE_FIELD_DR4;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR4]);
      result = 1;
    }
  }

  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR5)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= 30.0) && (lat > 15.0)) {
      releaseField = RELEASE_FIELD_DR5;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR5]);
      result = 1;
    }
  }

  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR6)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= 15.0) && (lat > 0.0)) {
      releaseField = RELEASE_FIELD_DR6;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR6]);
      result = 1;
    }
  }

  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR7)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= 0) && (lat > -15.0)) {
      releaseField = RELEASE_FIELD_DR7;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR7]);
      result = 1;
    }
  }
  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR8)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= -15.0) && (lat > -30.0)) {
      releaseField = RELEASE_FIELD_DR8;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR8]);
      result = 1;
    }
  }
  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR9)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= -30.0) && (lat > -45.0)) {
      releaseField = RELEASE_FIELD_DR9;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR9]);
      result = 1;
    }
  }
  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR10)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= -45.0) && (lat > -60.0)) {
      releaseField = RELEASE_FIELD_DR10;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR10]);
      result = 1;
    }
  }
  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR11)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= -60.0) && (lat > -75.0)) {
      releaseField = RELEASE_FIELD_DR11;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR11]);
      result = 1;
    }
  }
  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR12)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if (lat <= -75.0) {
      releaseField = RELEASE_FIELD_DR12;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR12]);
      result = 1;
    }
  }
  if (RELEASE_LEVEL >= RELEASE_LEVEL_M44) {

    if (result == 0) {
      ra0 = 130.0924;
      dec0 = 19.67206;
      if (wcsdist(ra0,dec0,ra,dec) < AUTHORIZE_ANGLE) {
        releaseField = RELEASE_FIELD_M44;
        strcpy(regionString,releaseFieldText[RELEASE_FIELD_M44]);
        result = 1; /* M44 */
      }
    }
  }
  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR1)) {
    ra0 = 290.73;
    dec0 = 44.498;
    if (wcsdist(ra0,dec0,ra,dec) < 8.0) {
      releaseField = RELEASE_FIELD_KEPLER;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_KEPLER]);
      result = 1; /* KEPLER */
    }
  }

  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR1)) {
    ra0 = 187.2778;
    dec0 = 2.0523;
    if (wcsdist(ra0,dec0,ra,dec) < AUTHORIZE_ANGLE) {
      releaseField = RELEASE_FIELD_3C273;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_3C273]);
      result = 1; /* 3C273 */
    }
  }


  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR1)) {
    ra0 = 270.8800;
    dec0 = -30.0200;
    dist = wcsdist(ra0,dec0,ra,dec);
    if (dist < AUTHORIZE_ANGLE) {
      releaseField = RELEASE_FIELD_BAADE;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_BAADE]);
      result = 1; /* BAADE */
    }
  }

  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_DR1)) {
    ra0 = 80.89386;
    dec0 = -69.7561;
    if (wcsdist(ra0,dec0,ra,dec) < AUTHORIZE_ANGLE) {
      releaseField = RELEASE_FIELD_LMC;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_LMC]);
      result = 1; /* LMC */
    }
  }


  if ((result == 0) && (RELEASE_LEVEL >= RELEASE_LEVEL_ALL)) {
    /* Now check galactic coordinates */
    lon = ra;
    lat = dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lon,&lat,2000.0);
    if ((lat <= 45.0) && (lat > 30.0)) {
      releaseField = RELEASE_FIELD_DR4;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR4]);
      result = 1;
    } else if ((lat <= 30.0) && (lat > 15.0)) {
      releaseField = RELEASE_FIELD_DR5;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR5]);
      result = 1;
    } else if ((lat <= +15.0) && (lat > 0.0)) {
      releaseField = RELEASE_FIELD_DR6;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR6]);
      result = 1;
    } else if ((lat <= 0.0) && (lat > -15.0)) {
      releaseField = RELEASE_FIELD_DR7;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR7]);
      result = 1;
    } else if ((lat <= -15.0) && (lat > -30.0)) {
      releaseField = RELEASE_FIELD_DR8;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR8]);
      result = 1;
    } else if ((lat <= -30.0) && (lat > -45.0)) {
      releaseField = RELEASE_FIELD_DR9;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR9]);
      result = 1;
    } else if ((lat <= -45.0) && (lat > -60.0)) {
      releaseField = RELEASE_FIELD_DR10;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR10]);
      result = 1;
    } else if ((lat <= -60) && (lat > -75)) {
      releaseField = RELEASE_FIELD_DR11;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR11]);
      result = 1;
    } else if ((lat <= -75.0) && (lat >= -90.0)) {
      releaseField = RELEASE_FIELD_DR12;
      strcpy(regionString,releaseFieldText[RELEASE_FIELD_DR12]);
      result = 1;
    }
    result = 1;
  }


#if 0
  printf("Authorize %d for %s with %s\n",result,regionString,regionFlag);
#endif
  *pReleaseField = releaseField;
  return(result);
}
/* Get a list of contiguous gsc bins around the current gsc bin */
void FindAdjacentBins(PGSCBIN pGscBin,int gsc_bin_index,int *gsc_bin_list,int *pgsc_bin_count) {
  int gsc_bin_count = 0;
  PBININDEX pBinIndex;
  PBININDEX pBinIndex2;
  int raBin;
  int decBin;
  double minRa;
  double maxRa;
  double deltaRa;

  double minRa2;
  double maxRa2;
  double deltaRa2;

  int tmp_gsc_bin_index;

  int raIndex;
  int raIndex2;
  int raIndex3;
  int decIndex;
  pBinIndex = GetSubbins(pGscBin,gsc_bin_index,&raBin,&decBin,(char *)"FindAdjacentBins");
  deltaRa = (360.00 /pBinIndex->numBins);
  minRa = ((360.00 * (1.0*raBin))/pBinIndex->numBins);
  maxRa = minRa + deltaRa;
  for (decIndex = (decBin -1); decIndex <= (decBin+1); decIndex++) {
    if ((decIndex < 0) || (decIndex >= pGscBin->dec_bins)) {
      continue;
    }
    pBinIndex2 =  &pGscBin->pBinMasterIndex[decIndex];
    deltaRa2 = (360.00 /pBinIndex2->numBins);
    if ((decIndex == 0) || (decIndex == (pGscBin->dec_bins-1))) {
      for (raIndex = 0; raIndex < pBinIndex2->numBins; raIndex++) {
        tmp_gsc_bin_index = pBinIndex2->startBin+raIndex;
        if (tmp_gsc_bin_index != gsc_bin_index) {
          gsc_bin_list[gsc_bin_count] = tmp_gsc_bin_index;
          gsc_bin_count++;
          if (gsc_bin_count >= MAX_ADJACENT_BINS) {
            printf("ERROR: MAX_ADJACENT_BINS exceeded\n");
            exit(-1);
          }
        }
      }
    } else {
      raIndex = (((minRa+maxRa)/2.0) * pBinIndex2->numBins)/360.0;
      for (raIndex2 = (raIndex-2); raIndex2 <= (raIndex+2); raIndex2++) {
        if (raIndex2 < 0) {
          raIndex3 = raIndex2 + pBinIndex2->numBins;
        } else if (raIndex2 >=  pBinIndex2->numBins) {
          raIndex3 = raIndex2 - pBinIndex2->numBins;
        } else {
          raIndex3 = raIndex2;
        }
        minRa2 = ((360.00 * (1.0*raIndex3))/pBinIndex2->numBins);
        maxRa2 = minRa2 + deltaRa2;
        if ((minRa2 - minRa) > 180.0) {
          minRa2 -= 360.0;
          maxRa2 -= 360.0;
        } else if ((minRa - minRa2) > 180.0) {
          minRa2 += 360.0;
          maxRa2 += 360.0;

        }
        if ((maxRa2 > (minRa-(deltaRa2/2))) && (minRa2 < (maxRa+(deltaRa2/2)))) {
          tmp_gsc_bin_index = pBinIndex2->startBin+raIndex3;
          if (tmp_gsc_bin_index != gsc_bin_index) {
            gsc_bin_list[gsc_bin_count] = tmp_gsc_bin_index;
            gsc_bin_count++;
            if (gsc_bin_count >= MAX_ADJACENT_BINS) {
              printf("ERROR: MAX_ADJACENT_BINS exceeded\n");
              exit(-1);
            }
          }
        }
      }
    }
  }
  *pgsc_bin_count =  gsc_bin_count;
}
int GetReleaseLevel()
{
  return(RELEASE_LEVEL);

}
void FormatFlagsBits(int AFLAGS, char *result,char *resultbits,int resultSize,int qualityFlag) {
  int length;
  char tmpchar[10];
  int index;
  int firstFlag = 0;
  /* Do not display individual bits for an unitialized quality bitmap */
  if ((qualityFlag != 0) && (AFLAGS == QUALITY_UNINITIALIZED)) {
    sprintf(result,"0");
    sprintf(resultbits,"UNINIT");
    return;
  }

  sprintf(result,"%d",AFLAGS);
  resultbits[0] = 0;
  for (index = 0; index < 31; index++) {
    if ((AFLAGS & 1 << index) != 0) {
      if (firstFlag == 0) {
        sprintf(tmpchar,"%d",index);
        firstFlag = 1;
      } else {
        sprintf(tmpchar,",%d",index);
      }
      length = strlen(resultbits)+strlen(tmpchar)+1;
      if (length > resultSize) {
        printf("ERROR: fatal size %d in formatAFlags for buffer size %d\n",length,resultSize);
        exit(-1);
      } else {
        strcat(resultbits,tmpchar);
      }
    }
  }
  return;

}
char * GetSeriesString(int seriesId,int fatal) {
  static char unknown[] = "Unknown";
#ifdef SKIP_MYSQL
  return("unknown");
#else /* SKIP_MYSQL */
  if (seriesTableXInit == 0) {
    printf("ERROR: Series table not initialized with a call to InitSeriesTable\n");
    if (fatal == 1) {
      exit(-1);
    } else {
      return(unknown);
    }
  }
  if ((seriesId <= 0) || (seriesId >= MAX_SERIES)) {
    printf("ERROR: Illegal series id %d in GetSeriesString\n",seriesId);
    if (fatal == 1) {
      exit(-1);
    } else {
      return(unknown);
    }
  }
  return(seriesTableX[seriesId].series);
#endif /* SKIP_MYSQL */
}
void local_strlwr(char *strPtr)
{
  char *curPtr = strPtr;
  char curChar;
  if (strPtr == NULL) {
    return;
  }
  while ((curChar = *curPtr) != 0) {
    if ((curChar >= 'A') && (curChar <= 'Z')) {
      curChar = curChar - 'A' + 'a';
      *curPtr = curChar;
    }
    curPtr++;

  }

}
const char* GetRejectReasonText(int index)
{
  PREJECTBIT rejectBitPointer;
  static const char badval[] = "Unknown";
  int index2;
  if ((index < 0) || (index >= REJECT_REASON_MAX)) {
    return(badval);
  }
  for (index2 = 0; index2 < REJECT_REASON_MAX; index2++) {
    rejectBitPointer = &rejectReasonBits[index2];
    if (rejectBitPointer->rejectBit == index) {
      return(rejectBitPointer->rejectDescr);
    }
  }
  return(badval);
}
const char* GetFatalReasonText(int index)
{
  PREJECTBIT fatalBitPointer;
  static const char badval[] = "Unknown";
  int index2;
  if ((index < 0) || (index >= FATAL_REASON_MAX)) {
    return(badval);
  }
  for (index2 = 0; index2 < FATAL_REASON_MAX; index2++) {
    fatalBitPointer = &fatalReasonBits[index2];
    if (fatalBitPointer->rejectBit == index) {
      return(fatalBitPointer->rejectDescr);
    }
  }
  return(badval);
}
int BlueColorterm(int catalogNumber,double colorterm)
{
  int result = 0; /* assume plate is not blue */
  switch(catalogNumber) {
  case CATALOG_GSC232:
    if ((colorterm >= MIN_BLUE_COLORTERM) &&
        (colorterm <= MAX_BLUE_COLORTERM)) {
      result = 1; /* a blue plate */
    }
    break;
  case CATALOG_KEPLER:
    if ((colorterm >= MIN_KEPLER_BLUE_COLORTERM) &&
        (colorterm <= MAX_KEPLER_BLUE_COLORTERM)) {
      result = 1; /* a blue plate */
    }
    break;
  case CATALOG_GAIA:
    if ((colorterm >= MIN_GAIA_BLUE_COLORTERM) &&
        (colorterm <= MAX_GAIA_BLUE_COLORTERM)) {
      result = 1; /* a blue plate */
    }
    break;
  case CATALOG_APASS:
    if ((colorterm >= MIN_APASS_BLUE_COLORTERM) &&
        (colorterm <= MAX_APASS_BLUE_COLORTERM)) {
      result = 1; /* a blue plate */

    }
    break;
  case CATALOG_ATLAS:
    if ((colorterm >= MIN_ATLAS_BLUE_COLORTERM) &&
        (colorterm <= MAX_ATLAS_BLUE_COLORTERM)) {
      result = 1; /* a blue plate */

    }
    break;
  case CATALOG_EXPERIMENTAL:
    if ((colorterm >= MIN_EXPERIMENTAL_BLUE_COLORTERM) &&
        (colorterm <= MAX_EXPERIMENTAL_BLUE_COLORTERM)) {
      result = 1; /* a blue plate */

    }
    break;
  default:
    printf("ERROR: pipelineutils line %d unknown catalog %d\n",__LINE__,catalogNumber);
    exit(-1);
    break;
  }

  return result;
}

int XDmagBins(int width,int height)
{
  if ((width*height) >= SMALL_PLATE_PIXELS) {
    return(X_DMAGBINS_NORMAL);
  }
  if (width > height) {
    return(LONG_DMAGBINS_SMALL);
  } else {
    return(SHORT_DMAGBINS_SMALL);
  }

}
int YDmagBins(int width,int height)
{
  if ((width*height) >= SMALL_PLATE_PIXELS) {
    return(Y_DMAGBINS_NORMAL);
  }
  if (width <= height) {
    return(LONG_DMAGBINS_SMALL);
  } else {
    return(SHORT_DMAGBINS_SMALL);
  }

}
void InitLinearityCommon(PLINCOMMON pLinearityCommon)
{
  if (pLinearityCommon->fileOpen == 1) {
    if (pLinearityCommon->binHandle) {
      fclose(pLinearityCommon->binHandle);
      pLinearityCommon->binHandle = NULL;
    }
    if (pLinearityCommon->forwardBuffer) {
      free(pLinearityCommon->forwardBuffer);
      pLinearityCommon->forwardBuffer = NULL;
    }
    if (pLinearityCommon->reverseBuffer) {
      free(pLinearityCommon->reverseBuffer);
      pLinearityCommon->reverseBuffer = NULL;
    }
  }
  memset(pLinearityCommon,0,sizeof(LINCOMMON));
  pLinearityCommon->algorithm = LIN_NO_ALGORITHM;
  pLinearityCommon->delay = -1;
  pLinearityCommon->offset = -1;
  pLinearityCommon->linVersion = 0;
}

/* NOTE: this routine is duplicated in ./FitsWrapper/ScannerCCD.cpp */
#define DELTA_VALUE 400 /* Amount to step back for final interpolation */
int ReadLinearityFile(PLINCOMMON pLinearityCommon,char *binFileName)
{
  LINHDRV2 linHdr;
  LINPOINT linPoint;
  int bytesRead;
  int tmpVal;
  int baseVal;
  int quadrant;
  int index;
  int index2;
  int index3;
  double increment;
  if ((pLinearityCommon->fileOpen != 0) ||
      (pLinearityCommon->binHandle != NULL) ||
      (pLinearityCommon->forwardBuffer != NULL) ||
      (pLinearityCommon->reverseBuffer != NULL)) {
    /* File is already open */
    return(-4);
  }
  InitLinearityCommon(pLinearityCommon);
  pLinearityCommon->binHandle = fopen(binFileName,"rb");
  if (pLinearityCommon->binHandle == NULL) {
    InitLinearityCommon(pLinearityCommon);
    return(-1);

  }
  memset(&linHdr,0,sizeof(LINHDRV2));
  bytesRead = fread(&linHdr,1,sizeof(LINHDRV1),pLinearityCommon->binHandle);
  if (bytesRead != sizeof(LINHDRV1)) {
    InitLinearityCommon(pLinearityCommon);
    return(-2);
  }

  if (linHdr.flag != LIN_FLAG) {
    InitLinearityCommon(pLinearityCommon);
    return(-3);
  }
  if (linHdr.version == LIN_VERSION_V1) {
    /* Invalidate unavailable fields */
    linHdr.delay = -1;
    linHdr.offset = -1;
    linHdr.linVersion = 0;
  } else if (linHdr.version == LIN_VERSION_V2) {
    /* Need to read in the rest of the header */
    bytesRead = fread(&(linHdr.delay),1,sizeof(LINHDRV2)-sizeof(LINHDRV1),pLinearityCommon->binHandle);
    if (bytesRead != (sizeof(LINHDRV2) - sizeof(LINHDRV1))) {
      InitLinearityCommon(pLinearityCommon);
      return(-20);
    }
  } else {
    /* Unrecognized version */
    InitLinearityCommon(pLinearityCommon);
    return(-4);
  }
  pLinearityCommon->algorithm = linHdr.algorithm;
  pLinearityCommon->maxValue = linHdr.size;
  pLinearityCommon->scanTime = linHdr.flatTime;
  pLinearityCommon->maxQuadrant = linHdr.quadrants;
  pLinearityCommon->delay = linHdr.delay;
  pLinearityCommon->linVersion = linHdr.linVersion;
  pLinearityCommon->offset = linHdr.offset;
  pLinearityCommon->forwardBuffer = (int*)malloc(pLinearityCommon->maxValue * pLinearityCommon->maxQuadrant * sizeof(int));
  if (pLinearityCommon->forwardBuffer == NULL) {
    InitLinearityCommon(pLinearityCommon);
    return(-5);
  }
  pLinearityCommon->reverseBuffer = (int*)malloc(pLinearityCommon->maxValue * pLinearityCommon->maxQuadrant * sizeof(int));
  if (pLinearityCommon->reverseBuffer == NULL) {
    InitLinearityCommon(pLinearityCommon);
    return(-6);
  }
  memset(&linPoint,0,sizeof(LINPOINT));
  for (quadrant = 0; quadrant < pLinearityCommon->maxQuadrant; quadrant++) {
    for (index = 0; index < pLinearityCommon->maxValue; index++) {
      bytesRead = fread(&linPoint,1,sizeof(LINPOINT),pLinearityCommon->binHandle);
      if (bytesRead != sizeof(LINPOINT)) {
        InitLinearityCommon(pLinearityCommon);
        return(-7);
      }
      if (linPoint.quadrant != quadrant) {
        InitLinearityCommon(pLinearityCommon);
        return(-8);
      }
      if (linPoint.expValue != index) {
        InitLinearityCommon(pLinearityCommon);
        return(-9);
      }
      pLinearityCommon->forwardBuffer[index + (pLinearityCommon->maxValue * quadrant)] = linPoint.actValue;
      pLinearityCommon->reverseBuffer[index + (pLinearityCommon->maxValue * quadrant)] = -1;
    }
  }

  /* We have the forward array, now construct the reverse array */
  for (quadrant = 0; quadrant < pLinearityCommon->maxQuadrant; quadrant++) {
    for (index = 0; index < pLinearityCommon->maxValue; index++) {
      tmpVal = pLinearityCommon->forwardBuffer[index + (pLinearityCommon->maxValue * quadrant)];
      if (tmpVal < 0) {
        InitLinearityCommon(pLinearityCommon);
        return(-10);
      }
      if (tmpVal < pLinearityCommon->maxValue) {
        pLinearityCommon->reverseBuffer[tmpVal + (pLinearityCommon->maxValue * quadrant)] = index;
      }
    }
    /* Now we need to interpolate for any missing values */
    baseVal = 0;
    for (index = 0; index < pLinearityCommon->maxValue; index++) {
      tmpVal = pLinearityCommon->reverseBuffer[index + (pLinearityCommon->maxValue * quadrant)];
      if (tmpVal >= 0) {
        baseVal = tmpVal;
      } else {
        for (index2 = (index+1); index2 < pLinearityCommon->maxValue; index2++) {
          tmpVal = pLinearityCommon->reverseBuffer[index2 + (pLinearityCommon->maxValue * quadrant)];
          if (tmpVal >= 0) {
            break;
          }
        }
        if (index2 == pLinearityCommon->maxValue) {
          /* No endpoint specified */
          if (index < (2 * DELTA_VALUE)) {
            InitLinearityCommon(pLinearityCommon);
            return(-11);
          }
          tmpVal = pLinearityCommon->reverseBuffer[(index - 1 - DELTA_VALUE) + (pLinearityCommon->maxValue * quadrant)];
          increment = 1.0 * (baseVal - tmpVal)/(DELTA_VALUE);
          tmpVal = (int)((increment *(index2-index+1)) + baseVal);
        }
        increment = 1.0 * (tmpVal - baseVal)/(index2 - index + 1);
        for (index3 = index; index3 < index2; index3++) {
          tmpVal = (int)((increment *(index3-index+1)) + baseVal);
          if (tmpVal >= pLinearityCommon->maxValue) {
            tmpVal = pLinearityCommon->maxValue - 1;
          }
          pLinearityCommon->reverseBuffer[index3 + (pLinearityCommon->maxValue * quadrant)] = tmpVal;
        }
      }
    }
    for (index = 0; index < pLinearityCommon->maxValue; index++) {
      tmpVal = pLinearityCommon->reverseBuffer[index + (pLinearityCommon->maxValue * quadrant)];
#ifdef DEBUG_PLOT_LINEARITY
      fprintf(gmtFile,"%4d %4d\n",index,tmpVal-index);
#endif /* DEBUG_PLOT_LINEARITY */
      if (tmpVal < 0) {
        InitLinearityCommon(pLinearityCommon);
        return(-12);
      }
    }


  }
#ifdef DEBUG_PLOT_LINEARITY
  fclose(gmtFile);
#endif /* DEBUG_PLOT_LINEARITY */

  pLinearityCommon->fileOpen = 1;
  return(0);
}
int GetPlateCondition(void *connection,char *series,int plateNumber,char *editdate,PPLATECONDITION* pPlateConditionTable,FILE *logHandle)
{
  MYSQL *pStacksConnection = (MYSQL *)connection;
  PPLATECONDITION plateConditionTable = *pPlateConditionTable;
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  char queryString[MAX_QUERY_STRING];
  int pickedCommentIndex;
  PPLATECONDITION pPlateCondition;
  int conditionIndex = 0;
  int numConditionRows = 0;


  if (plateConditionTable != NULL) {
    printf("ERROR: GetPlateCondition called with non-null plateConditionTable\n");
    return(-1);
  }
  sprintf(queryString,"SELECT versionId,pickedComment,pickednotes,valid from platecondition where series = '%s' and plateNumber = %d and date = '%s';",series,plateNumber,editdate);
  res = ExecuteQuery(pStacksConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pStacksConnection);
    if (res_ptr) {
      numConditionRows =  mysql_affected_rows(pStacksConnection);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (plateConditionTable == NULL) {
          plateConditionTable = (PPLATECONDITION)calloc(numConditionRows,sizeof(PLATECONDITION));
          if (plateConditionTable == NULL) {
            printf("ERROR: failed to allocate plateConditionTable of size %d\n",numConditionRows);
            fprintf(logHandle,"ERROR: failed to allocate plateConditionTable of size %d\n",numConditionRows);
            return(-1);
          }
          *pPlateConditionTable = plateConditionTable;
        }
        if (conditionIndex >= numConditionRows) {
          printf("ERROR: conditionIndex %d exceeds numConditionRows %d\n",conditionIndex,numConditionRows);
          fprintf(logHandle,"ERROR: conditionIndex %d exceeds numConditionRows %d\n",conditionIndex,numConditionRows);
          return(-1);
        }
        pPlateCondition = &plateConditionTable[conditionIndex];
        strcpy(pPlateCondition->series,series);
        pPlateCondition->plateNumber = plateNumber;
        strcpy(pPlateCondition->date,editdate);
        if (sqlrow[0] != NULL) {
          nvals = sscanf(sqlrow[0],"%d",&pPlateCondition->versionId);
          if (nvals != 1) {
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[1] != NULL) {
          for (pickedCommentIndex = 0; pickedCommentIndex < PICKED_COMMENT_MAX; pickedCommentIndex++) {
            if (strstr(sqlrow[1],pickedCommentString[pickedCommentIndex]) != NULL) {
              pPlateCondition->pickedComment |= 1 << pickedCommentIndex;
            }
          }
        } else {
          pPlateCondition->pickedComment = 0;
        }
        if (sqlrow[2] != NULL) {
          strncpy(pPlateCondition->pickednotes,sqlrow[2],MAX_STACKS_BUFFER);
          pPlateCondition->pickednotes[MAX_STACKS_BUFFER-1] = 0;
        } else {
          pPlateCondition->pickednotes[0] = 0;
        }
        if (sqlrow[3] != NULL) {
          if (strcmp(sqlrow[3],"yes") == 0) {
            pPlateCondition->valid = 1;
          } else if (strcmp(sqlrow[3],"no") == 0) {
            pPlateCondition->valid = 0;
          } else {
            printf("ERROR: pPlateCondition->valid %s is not 'yes' nor 'no'\n",sqlrow[3]);
          }

        } else {
          pPlateCondition->valid = -1;
        }

        conditionIndex++;
      }
      mysql_free_result(res_ptr);
      if (conditionIndex != numConditionRows) {
        printf("ERROR: numConditionRows %d does not match conditionIndex %d\n",numConditionRows,conditionIndex);
        fprintf(logHandle,"ERROR: numConditionRows %d does not match conditionIndex %d\n",numConditionRows,conditionIndex);
        return(-1);
      }
    } else {
      printf("ERROR: mysql_store_result failed at line %d\n",__LINE__);
      fprintf(logHandle,"ERROR: mysql_store_result failed at line %d\n",__LINE__);
      return(-1);
    }

  }
  return(conditionIndex);
}
int GetPlateEvents(void *connection,char *series,int plateNumber,char *editdate,PPLATEEVENTS* pPlateEventsTable,FILE *logHandle)
{
  MYSQL *pStacksConnection = (MYSQL *)connection;
  PPLATEEVENTS plateEventsTable = *pPlateEventsTable;
  int pickedStatusIndex;
  int nvals;
  PPLATEEVENTS pPlateEvents;
  int eventsIndex = 0;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  int numEventsRows = 0;
  char queryString[MAX_QUERY_STRING];


  if (plateEventsTable != NULL) {
    printf("ERROR: GetPlateEvents exited with non-null events table\n");
    return(-1);
  }
  sprintf(queryString,"SELECT versionId,pickedStatus,eventnotes,valid,stackLocation from plateevents where series = '%s' and plateNumber = %d and date = '%s';",series,plateNumber,editdate);
  res = ExecuteQuery(pStacksConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pStacksConnection);
    if (res_ptr) {
      numEventsRows =  mysql_affected_rows(pStacksConnection);
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (plateEventsTable == NULL) {
          plateEventsTable = (PPLATEEVENTS)calloc(numEventsRows,sizeof(PLATEEVENTS));
          if (plateEventsTable == NULL) {
            printf("ERROR: failed to allocate plateEventsTable of size %d\n",numEventsRows);
            fprintf(logHandle,"ERROR: failed to allocate plateEventsTable of size %d\n",numEventsRows);
            return(-1);
          }
          *pPlateEventsTable =  plateEventsTable;
        }
        if (eventsIndex >= numEventsRows) {
          printf("ERROR: eventsIndex %d exceeds numEventsRows %d\n",eventsIndex,numEventsRows);
          fprintf(logHandle,"ERROR: eventsIndex %d exceeds numEventsRows %d\n",eventsIndex,numEventsRows);
          return(-1);
        }
        pPlateEvents = &plateEventsTable[eventsIndex];
        strcpy(pPlateEvents->series,series);
        pPlateEvents->plateNumber = plateNumber;
        strcpy(pPlateEvents->date,editdate);
        if (sqlrow[0] != NULL) {
          nvals = sscanf(sqlrow[0],"%d",&pPlateEvents->versionId);
          if (nvals != 1) {
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[1] != NULL) {
          for (pickedStatusIndex = 0; pickedStatusIndex < PICKED_STATUS_MAX; pickedStatusIndex++) {
            if (strstr(sqlrow[1],pickedStatusString[pickedStatusIndex]) != NULL) {
              pPlateEvents->pickedStatus = pickedStatusIndex;
              break;
            }
          }
        } else {
          pPlateEvents->pickedStatus = 0;
        }
        if (sqlrow[2] != NULL) {
          strncpy(pPlateEvents->eventnotes,sqlrow[2],MAX_STACKS_BUFFER);
          pPlateEvents->eventnotes[MAX_STACKS_BUFFER-1] = 0;
        } else {
          pPlateEvents->eventnotes[0] = 0;
        }
        if (sqlrow[3] != NULL) {
          if (strcmp(sqlrow[3],"yes") == 0) {
            pPlateEvents->valid = 1;
          } else if (strcmp(sqlrow[3],"no") == 0) {
            pPlateEvents->valid = 0;
          } else {
            printf("ERROR: pPlateEvents->valid %s is not 'yes' nor 'no'\n",sqlrow[3]);
          }

        } else {
          pPlateEvents->valid = -1;
        }
        if (sqlrow[4] != NULL) {
          strncpy(pPlateEvents->stackLocation,sqlrow[4],MAX_STACK_LOCATION);
          pPlateEvents->stackLocation[MAX_STACK_LOCATION-1] = 0;
        } else {
          pPlateEvents->stackLocation[0] = 0;
        }

        eventsIndex++;
      }
      mysql_free_result(res_ptr);
      if (eventsIndex != numEventsRows) {
        printf("ERROR: numEventsRows %d does not match eventsIndex %d\n",numEventsRows,eventsIndex);
        fprintf(logHandle,"ERROR: numEventsRows %d does not match eventsIndex %d\n",numEventsRows,eventsIndex);
        return(-1);
      }
    } else {
      printf("ERROR: mysql_store_result failed at line %d\n",__LINE__);
      fprintf(logHandle,"ERROR: mysql_store_result failed at line %d\n",__LINE__);
      return(-1);
    }

  }
  return(eventsIndex);

}
int GetFlatfieldRecords(void *connection, char *startCalibrationDate,char* endCalibrationDate,int calibrationType, int calibrationStatus,int diskLocation, int* pJobId, int* pDelay,char *suffix,PFLATFIELDS *ppFlatfields)
{
   MYSQL *pConnection = connection;
  PFLATFIELDS pFlatfields = *ppFlatfields;
  PFLATFIELDS FlatfieldTable = NULL;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char queryString[MAX_QUERY_STRING];
  char tempString[MAX_QUERY_STRING];
  int recordCount = 0;
  int whereFlag = 0;
  int totalRecords;
  int nvals;

  if (pFlatfields != NULL) {
    printf("ERROR: GetFlatfieldRecords needs to be called with a null ppFlatfields field\n");
    return(-1);
  }



  sprintf(queryString,"SELECT calibrationDate,calibrationType+0,diskLocation,JobId,calibrationStatus+0,delay,alternateDate,suffix from flatfields");
  if ((startCalibrationDate != NULL) && (endCalibrationDate == NULL)) {
    if (whereFlag == 0) {
      sprintf(tempString," where calibrationDate = '%s'",startCalibrationDate);
      whereFlag = 1;
    } else {
      sprintf(tempString," and calibrationDate = '%s'",startCalibrationDate);
    }
    strcat(queryString,tempString);
  } else if ((startCalibrationDate != NULL) && (endCalibrationDate != NULL)) {

    if (whereFlag == 0) {
      sprintf(tempString," where calibrationDate >  '%s' and calibrationDate < '%s'",startCalibrationDate,endCalibrationDate);
      whereFlag = 1;
    } else {
      sprintf(tempString," and calibrationDate >  '%s' and calibrationDate < '%s'",startCalibrationDate,endCalibrationDate);
    }
    strcat(queryString,tempString);

  }

  if (calibrationType >= 0) {
    if (whereFlag == 0) {
      sprintf(tempString," where calibrationType = %d",calibrationType);
      whereFlag = 1;
    } else {
      sprintf(tempString," and calibrationType = %d",calibrationType);
    }
    strcat(queryString,tempString);
  }

  if (calibrationStatus >= 0) {
    if (whereFlag == 0) {
      sprintf(tempString," where calibrationStatus = %d",calibrationStatus);
      whereFlag = 1;
    } else {
      sprintf(tempString," and calibrationStatus = %d",calibrationStatus);
    }
    strcat(queryString,tempString);
  }

  if (diskLocation > 0) {
    if (whereFlag == 0) {
      sprintf(tempString," where diskLocation = %d",diskLocation);
      whereFlag = 1;
    } else {
      sprintf(tempString," and diskLocation = %d",diskLocation);
    }
    strcat(queryString,tempString);
  }

  if (pJobId != NULL) {
    if (whereFlag == 0) {
      sprintf(tempString," where JobId = %d",*pJobId);
      whereFlag = 1;
    } else {
      sprintf(tempString," and JobId = %d",*pJobId);
    }
    strcat(queryString,tempString);
  }

  if (pDelay != NULL) {
    if (whereFlag == 0) {
      sprintf(tempString," where delay = %d",*pDelay);
      whereFlag = 1;
    } else {
      sprintf(tempString," and delay = %d",*pDelay);
    }
    strcat(queryString,tempString);
  }
  if (suffix != NULL) {
    if (whereFlag == 0) {
      sprintf(tempString," where suffix = '%s'",suffix);
      whereFlag = 1;
    } else {
      sprintf(tempString," and suffix = '%s'",suffix);
    }
    strcat(queryString,tempString);
  }
  strcat(queryString," order by calibrationDate,calibrationType,suffix;");
  if (strlen(queryString) > MAX_QUERY_STRING) {
    printf("ERROR: MAX_QUERY_STRING exceeded by %lu in GetFlatfieldRecords\n",strlen(queryString));
    exit(-1);
  }

  res = ExecuteQuery(pConnection,queryString);

  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      totalRecords = (int)mysql_num_rows(res_ptr);
      if (totalRecords > 0) {
        FlatfieldTable = (PFLATFIELDS)calloc(totalRecords,sizeof(FLATFIELDS));
        if (FlatfieldTable == NULL) {
          printf("ERROR: failed to allocate FlatfieldTable of size %d\n",totalRecords);
          exit(-1);
        }
      }
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        pFlatfields = &FlatfieldTable[recordCount];
        memset(pFlatfields,0,sizeof(FLATFIELDS));
        if ((sqlrow[0] != NULL) && (strlen(sqlrow[0]) < MAX_DATE_STRING)) {
          strcpy(pFlatfields->calibrationDate,sqlrow[0]);
        } else {
          continue;
        }
        if (sqlrow[1] != NULL) {
          nvals = sscanf(sqlrow[1],"%d",&pFlatfields->calibrationType);
          if (nvals != 1) {
            printf("ERROR: illegal calibrationType %s in GetFlatfieldRecords\n",sqlrow[1]);
            pFlatfields->calibrationType = 0;
          }
        } else {
          pFlatfields->calibrationType = 0;
        }
        if (sqlrow[2] != NULL) {
          nvals = sscanf(sqlrow[2],"%d",&pFlatfields->diskLocation);
          if (nvals != 1) {
            printf("ERROR: illegal diskLocation %s in GetFlatfieldRecords\n",sqlrow[2]);
            pFlatfields->diskLocation = -1;
          }
        } else {
          pFlatfields->diskLocation = -1;
        }

        if (sqlrow[3] != NULL) {
          nvals = sscanf(sqlrow[3],"%d",&pFlatfields->JobId);
          if (nvals != 1) {
            printf("ERROR: illegal JobId %s in GetFlatfieldRecords\n",sqlrow[3]);
            pFlatfields->JobId = -1;
          }
        } else {
          pFlatfields->JobId = -1;
        }


        if (sqlrow[4] != NULL) {
          nvals = sscanf(sqlrow[4],"%d",&pFlatfields->calibrationStatus);
          if (nvals != 1) {
            printf("ERROR: illegal calibrationStatus %s in GetFlatfieldRecords\n",sqlrow[4]);
            pFlatfields->calibrationStatus = 0;
          }
        } else {
          pFlatfields->calibrationStatus = 0;
        }
        if (sqlrow[5] != NULL) {
          nvals = sscanf(sqlrow[5],"%d",&pFlatfields->delay);
          if (nvals != 1) {
            printf("ERROR: illegal delay %s in GetFlatfieldRecords\n",sqlrow[5]);
            pFlatfields->delay = 0;
          }
        } else {
          pFlatfields->delay = 0;
        }
        if ((sqlrow[6] != NULL) && (strlen(sqlrow[6]) < MAX_DATE_STRING)) {
          strcpy(pFlatfields->alternateDate,sqlrow[6]);
        }
        if ((sqlrow[7] != NULL) && (strlen(sqlrow[7]) < MAX_FLATFIELDS_SUFFIX)) {
          strcpy(pFlatfields->suffix,sqlrow[7]);
        }
        recordCount++;

      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: pipelineutils line %d: mysql_store_result failed",__LINE__);
    }
  } else {
    printf("ERROR: pipelineutils line %d mysql_query_failed",__LINE__);
  }
  *ppFlatfields = FlatfieldTable;
  return(recordCount);

}
/* If there is no valid mosaic associated with this plate, this routine clears all stale mosaic info */
void ResetMosaicList(PMOSAICLIST pMosaicList) {

  pMosaicList->cd1_1 = -1;
  pMosaicList->cd1_2 = -1;
  pMosaicList->cd2_1 = -1;
  pMosaicList->cd2_2 = -1;
  pMosaicList->crpix1 = -1;
  pMosaicList->crpix2 = -1;
  pMosaicList->crval1 = -1;
  pMosaicList->crval2 = -1;
  pMosaicList->diskLocation = -1;
  pMosaicList->mosaicNumber = -1;
  pMosaicList->naxis1 = -1;
  pMosaicList->naxis2 = -1;
  pMosaicList->nx = -1;
  pMosaicList->ny = -1;
  pMosaicList->patternID = -1;
  pMosaicList->rotation = -1;
  pMosaicList->scanNumber = -1;
  pMosaicList->solutionNumber = -1;

  /* Bitmaps */
  pMosaicList->FitWCS  = 0;
  pMosaicList->WCSSource = 0;
  pMosaicList->binning = 0;
  pMosaicList->transform = 0;
  /* Strings */
  pMosaicList->ctype1[0] = 0;
  pMosaicList->ctype2[0] = 0;
  pMosaicList->mosaicComment[0] = 0;
  pMosaicList->scanDate[0] = 0;
  pMosaicList->mosaicDate[0] = 0;

}


/* Returns 1 on success; 0 on failure */
int
BuildMosaicList(
  void *connection,
  char *series,     /* if NULL or zero length, consider all series */
  int plateNumber,  /* if -1, consider all plates */
  int mosaicNumber,  /* if -1, consider all mosaics */
  int solutionNumber, /* if -1, consider all solutions */
  int pendingFlag, /* if 1, then return unscanned plates */
  int mosaicInfoFlag,  /* if 1, being called grom GetMosaicInfo() */
  int selectAllFlag,   /* if 1, then return all mosaic records: good, deleted, and stale */
  int verbose,         /* if 1, display warnings */
  int *pTotalSolutionsCount,
  int *pTotalSolution0Count,
  int *pMosaicListSize,
  PMOSAICLIST *ppMosaicTable
) {
  int nvals;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char queryString[MAX_QUERY_STRING];
  char tempString[MAX_QUERY_STRING];
  MYSQL *pConnection = connection;
  MOSAICLIST mosaicList;
  PMOSAICLIST pMosaicList = &mosaicList;
  int numMosaicRecords = 0;
  int mosaicRecordIndex = 0;
  int mosaicRecordIndex3;
  int oldPlateNumber = -1;
  int oldMosaicNumber = -1;
  int oldMosaicRecordIndex = -1;

  PMOSAICLIST pMosaicTableEntry;
  PMOSAICLIST pMosaicTableEntry2;
  PMOSAICLIST pMosaicTableEntry3;
  PMOSAICLIST pMosaicTable = NULL;
  PMOSAICLIST pMosaicTableEntryDeleted;

  int numExposureRecords = 0;
  int exposureRecordIndex = 0;
  PMOSAICLIST pExposureTable = NULL;
  PMOSAICLIST pExposureTableEntry;

  int numOutputRecords = 0;
  int numOutputAlloc;
  int outputRecordIndex = 0;
  int outputRecordIndex2 = 0;
  int curOutputRecordIndex = 0;

  PMOSAICLIST pOutputTable = NULL;
  PMOSAICLIST pOutputTableEntry;
  PMOSAICLIST pOutputTableEntry2;

  int queryLength;
  char oldseries[MAX_SERIES_STRING]; /* Series */
  SERIES_HEADER seriesHeader;
  PSERIES_HEADER pSeriesHeader = &seriesHeader;

  int validEntry;
  int newValidEntry;
  int maxMosaicNumber;
  int strcmpResult;
  int totalSolutionsCount = 0;
  int totalSolution0Count = 0;

  *pMosaicListSize = 0;
  *pTotalSolutionsCount = 0;
  *pTotalSolution0Count = 0;

  if (*ppMosaicTable != NULL) {
    printf("ERROR: line %2d pipelineutils mosaicTable already allocated\n",__LINE__);
    return(0);
  }

  sprintf(
    queryString,
    "SELECT series,exposureNumber,binning+0,mosaicDate,naxis1,naxis2,ctype1,"
    "ctype2,crval1,crval2,crpix1,crpix2,cd1_1,cd1_2,cd2_1,cd2_2,FitWCS+0,"
    "rotation,WCSSource+0,mosaics.diskLocation,mosaicComment,transform+0,"
    "solutionNumber,scanDate,patternID,scanNumber,plateNumber,mosaicNumber,"
    "solutionNumber,mosaics.JobId "
    "FROM mosaics INNER JOIN scans using (series,plateNumber,scanNumber)"
  );

  if (series != NULL && *series) {
    sprintf(tempString, " WHERE series = '%s'", series);
    strcat(queryString, tempString);
  } else {
    sprintf(tempString, " WHERE series != ''");
    strcat(queryString, tempString);
  }

  if (plateNumber > 0) {
    sprintf(tempString, " AND plateNumber = %d", plateNumber);
    strcat(queryString, tempString);
  }

  if (mosaicNumber >= 0) {
    sprintf(tempString, " AND mosaicNumber = %d", mosaicNumber);
    strcat(queryString, tempString);
  }

  if (solutionNumber >= 0) {
    sprintf(tempString, " AND solutionNumber = %d", solutionNumber);
    strcat(queryString, tempString);
  }

  strcat(queryString, " ORDER BY series,plateNumber,scanDate,mosaicNumber,exposureNumber,solutionNumber;");

  queryLength = strlen(queryString);
  if (queryLength >= MAX_QUERY_STRING) {
    printf("ERROR: line %d pipelineutilsMAX_QUERY_STRING length exceeded (2) from %d\n", __LINE__, queryLength);
    exit(1);
  }

  res = ExecuteQuery(pConnection, queryString);

  if (!res) {
    res_ptr = mysql_store_result(pConnection);

    if (res_ptr) {
      numMosaicRecords = (int) mysql_num_rows(res_ptr);
      pMosaicTable = (PMOSAICLIST) calloc(numMosaicRecords + 1, sizeof(MOSAICLIST));
      if (pMosaicTable == NULL) {
        printf("line %4d ERROR: failed to allocate pMosaicTable\n", __LINE__);
        exit(1);
      }

      mosaicRecordIndex = 0;
      oldseries[0] =  0;
      oldMosaicRecordIndex = -1;
      oldPlateNumber = -1;
      oldMosaicNumber = -1;
      pMosaicTableEntry2 = NULL;

      while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
        memset(pMosaicList, 0, sizeof(MOSAICLIST));
        pMosaicList->solutionNumber = solutionNumber;

        if (sqlrow[0]) {
          strcpy(pMosaicList->series, sqlrow[0]);
        } else {
          pMosaicList->series[0] = 0;
        }

        if (strcmp(pMosaicList->series, oldseries) != 0) {
          strcpy(oldseries,pMosaicList->series);
          if (GetSeriesInfo(pConnection,pMosaicList->series,-1,pSeriesHeader) != 1) {
            if (strcmp(pMosaicList->series,"mask") == 0) {
              memset(pSeriesHeader,0,sizeof(SERIES_HEADER));
              strcpy(pSeriesHeader->series,"mask");
              pSeriesHeader->seriesId = -1;
              pSeriesHeader->orientation = -1;
              pSeriesHeader->sequestered = SEQUESTERED_NO;
              pSeriesHeader->nominalPlateScale = -1;
              pSeriesHeader->fittedPlateScale = -1;
            } else {
              printf("ERROR: line %4d BuildMosaicList failed to get series information for %s\n",__LINE__,series);
              exit(-1);
            }
          }
        }

        if (pSeriesHeader->sequestered != SEQUESTERED_NO) {
          continue;
        }

        pMosaicList->plateScale = pSeriesHeader->nominalPlateScale;

        if (pSeriesHeader->fittedPlateScale != 0) {
          pMosaicList->plateScale = pSeriesHeader->fittedPlateScale;
        }

        pMosaicList->orientation = pSeriesHeader->orientation;
        pMosaicList->seriesId = pSeriesHeader->seriesId;

        if (sqlrow[1]) {
          nvals = sscanf(sqlrow[1], "%d", &pMosaicList->exposureNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for exposureNumber", nvals);
          }
        }

        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2], "%d", &pMosaicList->binning);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for binning", nvals);
          }
        }

        if (sqlrow[3]) {
          strcpy(pMosaicList->mosaicDate, sqlrow[3]);
        }

        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4], "%d", &pMosaicList->naxis1);
        }

        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5], "%d", &pMosaicList->naxis2);
        }

        pMosaicList->nx = XDmagBins(pMosaicList->naxis1, pMosaicList->naxis2);
        pMosaicList->ny = YDmagBins(pMosaicList->naxis1, pMosaicList->naxis2);

        if (sqlrow[6]) {
          strcpy(pMosaicList->ctype1, sqlrow[6]);
        }

        if (sqlrow[7]) {
          strcpy(pMosaicList->ctype2, sqlrow[7]);
        }

        if (sqlrow[8]) {
          nvals = sscanf(sqlrow[8], "%lf", &pMosaicList->crval1);
        }

        if (sqlrow[9]) {
          nvals = sscanf(sqlrow[9], "%lf", &pMosaicList->crval2);
        }

        if (sqlrow[10]) {
          nvals = sscanf(sqlrow[10], "%lf", &pMosaicList->crpix1);
        }

        if (sqlrow[11]) {
          nvals = sscanf(sqlrow[11], "%lf", &pMosaicList->crpix2);
        }

        if (sqlrow[12]) {
          nvals = sscanf(sqlrow[12], "%lf", &pMosaicList->cd1_1);
        }

        if (sqlrow[13]) {
          nvals = sscanf(sqlrow[13], "%lf", &pMosaicList->cd1_2);
        }

        if (sqlrow[14]) {
          nvals = sscanf(sqlrow[14], "%lf", &pMosaicList->cd2_1);
        }

        if (sqlrow[15]) {
          nvals = sscanf(sqlrow[15], "%lf", &pMosaicList->cd2_2);
        }

        if (sqlrow[16]) {
          nvals = sscanf(sqlrow[16], "%d", &pMosaicList->FitWCS);
        }

        if (sqlrow[17]) {
          nvals = sscanf(sqlrow[17], "%d", &pMosaicList->rotation);
        }

        if (sqlrow[18]) {
          nvals = sscanf(sqlrow[18], "%d", &pMosaicList->WCSSource);
        }

        if (sqlrow[19]) {
          nvals = sscanf(sqlrow[19], "%d", &pMosaicList->diskLocation);
        }

        if (sqlrow[20]) {
          strcpy(pMosaicList->mosaicComment,sqlrow[20]);
          if (selectAllFlag == 0) {
            if ((strstr(pMosaicList->mosaicComment, "Deleted") != NULL) ||
                (strstr(pMosaicList->mosaicComment, "deleted") != NULL)) {
              pMosaicList->deleted = 1;
            }
          }
        }

        if (sqlrow[21]) {
          nvals = sscanf(sqlrow[21], "%d", &pMosaicList->transform);
          if (nvals != 1) {
            pMosaicList->transform = -1;
          }
        } else {
          pMosaicList->transform = -1;
        }

        if (sqlrow[23]) {
          strcpy(pMosaicList->scanDate,sqlrow[23]);
        } else {
          pMosaicList->scanDate[0] = 0;
        }

        if (sqlrow[24]) {
          nvals = sscanf(sqlrow[24], "%d", &pMosaicList->patternID);
          if (nvals != 1) {
            pMosaicList->patternID = -1;
          }
        } else {
          pMosaicList->patternID = -1;
        }

        if (sqlrow[25]) {
          nvals = sscanf(sqlrow[25], "%d", &pMosaicList->scanNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for scanNumber",nvals);
          }
        }

        if (sqlrow[26]) {
          nvals = sscanf(sqlrow[26], "%d", &pMosaicList->plateNumber);
          if (nvals != 1) {
            pMosaicList->plateNumber = -1;
          }
        } else {
          pMosaicList->plateNumber = -1;
        }

        if (sqlrow[27]) {
          nvals = sscanf(sqlrow[27], "%d", &pMosaicList->mosaicNumber);
          if (nvals != 1) {
            pMosaicList->mosaicNumber = -1;
          }
        } else {
          pMosaicList->mosaicNumber = -1;
        }

        if (sqlrow[28]) {
          nvals = sscanf(sqlrow[28], "%d", &pMosaicList->solutionNumber);
          if (nvals != 1) {
            pMosaicList->solutionNumber = -1;
          }
        } else {
          pMosaicList->solutionNumber = -1;
        }

        if (sqlrow[29]) {
          nvals = sscanf(sqlrow[29], "%d", &pMosaicList->JobId);
          if (nvals != 1) {
            pMosaicList->JobId = -1;
          }
        } else {
          pMosaicList->JobId = -1;
        }

        /* If we got this far, save the result */
        pMosaicTableEntry = &pMosaicTable[mosaicRecordIndex];
        memcpy(pMosaicTableEntry,pMosaicList,sizeof(MOSAICLIST));

        /* A second fixup checks the the selected mosaicNumber's mosaicComment to be sure that "Deleted" is set for all instances */
        if (oldPlateNumber <  0) {
          oldPlateNumber = pMosaicTableEntry->plateNumber;
          oldMosaicNumber = pMosaicTableEntry->mosaicNumber;
          oldMosaicRecordIndex = mosaicRecordIndex;
          pMosaicTableEntry2 = pMosaicTableEntry;
        } else {
          if ((oldPlateNumber == pMosaicTableEntry->plateNumber) &&
              (oldMosaicNumber == pMosaicTableEntry->mosaicNumber)) {
            if (pMosaicTableEntry->deleted != pMosaicTableEntry2->deleted) {
              if (verbose) {
                if (pMosaicTableEntry->deleted != 0) {
                  pMosaicTableEntryDeleted = pMosaicTableEntry;
                } else {
                  pMosaicTableEntryDeleted = pMosaicTableEntry2;
                }

                printf("line %4d WARNING: inconsistent deletion for %5s%05d_%02d s%02d exposure %02d mosaicComment '%s' '%s'\n",
                       __LINE__,
                       pMosaicTableEntry->series,
                       pMosaicTableEntry->plateNumber,
                       pMosaicTableEntry->mosaicNumber,
                       pMosaicTableEntry->solutionNumber,
                       pMosaicTableEntry->exposureNumber,
                       pMosaicTableEntry->mosaicComment,
                       pMosaicTableEntry2->mosaicComment);
                printf("SELECT series,plateNumber,mosaicNumber,solutionNumber,exposureNumber,mosaicComment from mosaics where series = '%s' and plateNumber = %d and mosaicNumber = %d;\n",
                       pMosaicTableEntry->series,
                       pMosaicTableEntry->plateNumber,
                       pMosaicTableEntry->mosaicNumber);
                printf("UPDATE mosaics set mosaicComment = '%s' where series = '%s' and plateNumber = %d and mosaicNumber = %d;\n",
                       pMosaicTableEntryDeleted->mosaicComment,
                       pMosaicTableEntry->series,
                       pMosaicTableEntry->plateNumber,
                       pMosaicTableEntry->mosaicNumber);
              }

              /* Now fix things up */
              if (pMosaicTableEntry2->deleted == 0) {
                if (strlen(pMosaicTableEntry2->mosaicComment) == 0) {
                  strcpy(pMosaicTableEntry2->mosaicComment,pMosaicTableEntry->mosaicComment);
                  pMosaicTableEntry2->deleted = 1;
                } else {
                  strcat(pMosaicTableEntry2->mosaicComment,", ");
                  strcat(pMosaicTableEntry2->mosaicComment,pMosaicTableEntry->mosaicComment);
                  pMosaicTableEntry2->deleted = 1;
                  pMosaicTableEntry2 = pMosaicTableEntry;
                }

                if (pMosaicTableEntry2->cd1_1 != 0) {
                  totalSolutionsCount--;
                  if (pMosaicTableEntry2->solutionNumber == 0) {
                    totalSolution0Count--;
                  }
                }
              } else {
                if (strlen(pMosaicTableEntry->mosaicComment) == 0) {
                  strcpy(pMosaicTableEntry->mosaicComment,pMosaicTableEntry2->mosaicComment);
                  pMosaicTableEntry->deleted = 1;
                } else {
                  strcat(pMosaicTableEntry->mosaicComment,", ");
                  strcat(pMosaicTableEntry->mosaicComment,pMosaicTableEntry2->mosaicComment);
                  pMosaicTableEntry->deleted = 1;
                }
              }

              for (mosaicRecordIndex3 = oldMosaicRecordIndex; mosaicRecordIndex3 < mosaicRecordIndex; mosaicRecordIndex3++) {
                 pMosaicTableEntry3 = &pMosaicTable[mosaicRecordIndex3];
                 if (pMosaicTableEntry3->deleted == 0) {
                   if (strlen(pMosaicTableEntry3->mosaicComment) == 0) {
                     strcpy(pMosaicTableEntry3->mosaicComment,pMosaicTableEntry2->mosaicComment);
                     pMosaicTableEntry3->deleted = 1;
                   } else {
                     strcat(pMosaicTableEntry3->mosaicComment,", ");
                     strcat(pMosaicTableEntry3->mosaicComment,pMosaicTableEntry2->mosaicComment);
                     pMosaicTableEntry3->deleted = 1;
                   }
                 }
              }
            }
          } else {
            oldPlateNumber = pMosaicTableEntry->plateNumber;
            oldMosaicNumber = pMosaicTableEntry->mosaicNumber;
            oldMosaicRecordIndex = mosaicRecordIndex;
            pMosaicTableEntry2 = pMosaicTableEntry;
          }
        }

        if ((pMosaicTableEntry->deleted == 0) && (pMosaicTableEntry->cd1_1 != 0)) {
          totalSolutionsCount++;
          if (pMosaicTableEntry->solutionNumber == 0) {
            totalSolution0Count++;
          }
        }

        mosaicRecordIndex++;
        if (mosaicRecordIndex > numMosaicRecords) {
          printf("ERROR: pipelineutils line %4d mosaicRecordIndex %d %d exceeded\n",__LINE__,mosaicRecordIndex,numMosaicRecords);
          exit(-1);
        }
      }

      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
    }
  }

  if (mosaicInfoFlag == 0 && selectAllFlag == 0) {
    sprintf(
      queryString,
      "SELECT series,exposureNumber,dRightAscension,dDeclination,date,class,"
      "exposure,plateNumber,centerSource+0 "
      "FROM exposures INNER JOIN plates USING (series,plateNumber)"
    );

    if (series != NULL || series[0] != 0) {
      sprintf(tempString, " WHERE series = '%s'", series);
      strcat(queryString, tempString);
    }

    if (plateNumber > 0) {
      sprintf(tempString," AND plateNumber = %d", plateNumber);
      strcat(queryString, tempString);
    }

    strcat(queryString," ORDER BY series,plateNumber,exposureNumber;");

    res = ExecuteQuery(pConnection, queryString);

    if (!res) {
      res_ptr = mysql_store_result(pConnection);

      if (res_ptr) {
        numExposureRecords = (int) mysql_num_rows(res_ptr);
        pExposureTable = (PMOSAICLIST) calloc(numExposureRecords + 1, sizeof(MOSAICLIST));
        if (pExposureTable == NULL) {
          printf("line %4d ERROR: failed to allocate pExposureTable\n",__LINE__);
          exit(1);
        }

        exposureRecordIndex = 0;
        oldseries[0] = 0;

        while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
          memset(pMosaicList,0,sizeof(MOSAICLIST));
          ResetMosaicList(pMosaicList);

          if (sqlrow[0]) {
            strcpy(pMosaicList->series,sqlrow[0]);
          } else {
            pMosaicList->series[0] = 0;
          }

          if (strcmp(pMosaicList->series,oldseries) != 0) {
            strcpy(oldseries,pMosaicList->series);
            if (GetSeriesInfo(pConnection,pMosaicList->series,-1,pSeriesHeader) != 1) {
              printf("ERROR: line %4d BuildMosaicList failed to get series information for %s\n",__LINE__,series);
              exit(-1);
            }
          }

          if (pSeriesHeader->sequestered != SEQUESTERED_NO) {
            continue;
          }

          pMosaicList->plateScale = pSeriesHeader->nominalPlateScale;

          if (pSeriesHeader->fittedPlateScale != 0) {
            pMosaicList->plateScale = pSeriesHeader->fittedPlateScale;
          }

          pMosaicList->orientation = pSeriesHeader->orientation;
          pMosaicList->seriesId = pSeriesHeader->seriesId;

          if (sqlrow[1]) {
            nvals = sscanf(sqlrow[1],"%d",&pMosaicList->exposureNumber);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for exposureNumber",nvals);
            }
          }

          if (sqlrow[2]) {
            nvals = sscanf(sqlrow[2],"%lf",&pMosaicList->dRightAscension);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for dRightAscension\n", nvals);
              pMosaicList->dRightAscension = 999.0;
              strcpy(pMosaicList->rightAscension,DASCH_NOTAVAILABLE);
            }
          } else {
            pMosaicList->dRightAscension = 999.0;
            strcpy(pMosaicList->rightAscension,DASCH_NOTAVAILABLE);
          }

          if (sqlrow[3]) {
            nvals = sscanf(sqlrow[3],"%lf",&pMosaicList->dDeclination);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for dDeclination\n", nvals);
              pMosaicList->dDeclination = 99.0;
              strcpy(pMosaicList->declination,DASCH_NOTAVAILABLE);
            }
          } else {
            pMosaicList->dDeclination = 99.0;
            strcpy(pMosaicList->declination, DASCH_NOTAVAILABLE);
          }

          if (sqlrow[4]) {
            strcpy(pMosaicList->date,sqlrow[4]);
          } else {
            pMosaicList->date[0] = 0;
          }

          if (sqlrow[5]) {
            strcpy(pMosaicList->plateClass,sqlrow[5]);
          } else {
            pMosaicList->plateClass[0] = 0;
          }

          if (sqlrow[6]) {
            strncpy(pMosaicList->exposure,sqlrow[6], sizeof(pMosaicList->exposure) - 1);
          } else {
            pMosaicList->exposure[0] = 0;;
          }

          if (sqlrow[7]) {
            nvals = sscanf(sqlrow[7],"%d",&pMosaicList->plateNumber);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for exposureNumber",nvals);
            }
          }

          if (sqlrow[8]) {
            nvals = sscanf(sqlrow[8],"%d",&pMosaicList->centerSource);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for exposureNumber",nvals);
            }
          }

          /* If we got this far, save the result */
          pExposureTableEntry = &pExposureTable[exposureRecordIndex];
          memcpy(pExposureTableEntry,pMosaicList,sizeof(MOSAICLIST));

          exposureRecordIndex++;
          if (exposureRecordIndex > numExposureRecords) {
            printf("ERROR: pipelineutils line %4d exposureRecordIndex %d %d exceeded\n",__LINE__,exposureRecordIndex,numExposureRecords);
            exit(-1);
          }
        }
      }

      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
    }
  }  /* End of exposure table processing */

  numOutputRecords = numMosaicRecords + numExposureRecords;

  if (numOutputRecords == 0) {
    if (pMosaicTable != NULL) {
      free(pMosaicTable);
    }

    if (pExposureTable != NULL) {
      free(pExposureTable);
    }

    return 0;
  }

  numMosaicRecords = mosaicRecordIndex;
  numExposureRecords = exposureRecordIndex;

  /* Now select the best mosaic */

  pOutputTable = (PMOSAICLIST) calloc(numOutputRecords + 1, sizeof(MOSAICLIST));
  if (pOutputTable == NULL) {
    printf("line %4d ERROR: failed to allocate pOutputTable\n",__LINE__);
    exit(1);
  }

  numOutputAlloc = numOutputRecords;
  outputRecordIndex = 0;
  curOutputRecordIndex = 0;

  /* Now go through the mosaic list and select the best entry: highest mosaic
   * number with the latest scan date. In addition, find the current highest
   * mosaic number */

  memset(pMosaicList,0,sizeof(MOSAICLIST));
  maxMosaicNumber = -1;
  validEntry = 0;

  for (mosaicRecordIndex = 0; mosaicRecordIndex < numMosaicRecords; mosaicRecordIndex++) {
    pMosaicTableEntry = &pMosaicTable[mosaicRecordIndex];

    if (
      (pMosaicTableEntry->plateNumber != pMosaicList->plateNumber) ||
      (strcmp(pMosaicTableEntry->series,pMosaicList->series) != 0)
    ) {
      for (outputRecordIndex2 = curOutputRecordIndex; outputRecordIndex2 < outputRecordIndex; outputRecordIndex2++) {
        pOutputTableEntry = &pOutputTable[outputRecordIndex2];
        pOutputTableEntry->maxMosaicNumber = maxMosaicNumber;
      }

      validEntry = 0;
      maxMosaicNumber = -1;
      memcpy(pMosaicList,pMosaicTableEntry,sizeof(MOSAICLIST));
      curOutputRecordIndex = outputRecordIndex; /* Save our place */
    }

    if ( maxMosaicNumber < pMosaicTableEntry->mosaicNumber) {
      maxMosaicNumber = pMosaicTableEntry->mosaicNumber;
    }

    if ((selectAllFlag == 0) &&
        ((strstr(pMosaicTableEntry->mosaicComment,"Deleted") != NULL) ||
         (strstr(pMosaicTableEntry->mosaicComment,"deleted") != NULL))) {
      pMosaicTableEntry->deleted = 1;
      newValidEntry = 0;
    } else {
      newValidEntry = 1;
    }

    if (validEntry == 0) {
      memcpy(pMosaicList,pMosaicTableEntry,sizeof(MOSAICLIST));
    }

    if (newValidEntry == 0) {
      continue;
    }

    validEntry = newValidEntry;
    strcmpResult = strcmp(pMosaicList->scanDate,pMosaicTableEntry->scanDate);

    if (strcmpResult < 0 && selectAllFlag == 0) {
      memcpy(pMosaicList,pMosaicTableEntry,sizeof(MOSAICLIST));
      outputRecordIndex = curOutputRecordIndex;       /* New mosaic with scanDate in the future. Restore the old index */
    }

    if (strcmpResult > 0) {
      printf("ERROR: line %d BuildMosaicList scanDate %s %s\n",__LINE__,pMosaicList->scanDate,pMosaicTableEntry->scanDate);
      exit(1);
    }

    if (pMosaicList->scanNumber != pMosaicTableEntry->scanNumber && selectAllFlag == 0) {
      printf("ERROR: line %d BuildMosaicList scanNumber %d %d\n",__LINE__,pMosaicList->scanNumber,pMosaicTableEntry->scanNumber);
      exit(1);
    }

    if (pMosaicList->mosaicNumber != pMosaicTableEntry->mosaicNumber && selectAllFlag == 0) {
      memcpy(pMosaicList,pMosaicTableEntry,sizeof(MOSAICLIST));
      outputRecordIndex = curOutputRecordIndex; /* New mosaic for same scan date, restore old index */
    }

    memcpy(pMosaicList,pMosaicTableEntry,sizeof(MOSAICLIST));
    pOutputTableEntry = &pOutputTable[outputRecordIndex];
    memcpy(pOutputTableEntry,pMosaicList,sizeof(MOSAICLIST));

    if (outputRecordIndex > numOutputRecords) {
      printf("ERROR: pipelineutils line %4d outputRecordIndex %d %d exceeded\n",__LINE__,outputRecordIndex,numOutputRecords);
      exit(1);
    }

    if (outputRecordIndex > numOutputRecords) {
      printf("ERROR: pipelineutils line %4d exposureRecordIndex %d %d exceeded\n",__LINE__,exposureRecordIndex,numExposureRecords);
      exit(1);
    }

    outputRecordIndex++;
  }

  for (outputRecordIndex2 = curOutputRecordIndex; outputRecordIndex2 < outputRecordIndex; outputRecordIndex2++) {
    pOutputTableEntry = &pOutputTable[outputRecordIndex2];
    pOutputTableEntry->maxMosaicNumber = maxMosaicNumber;
  }

  if (validEntry == 0) {
    outputRecordIndex = curOutputRecordIndex;
  }

  numOutputRecords = outputRecordIndex;

  /* The following is a special fixup: the default MySQL exposureNumber is 0 and
   * it should have been -1. As a result, we can have two entries for two
   * solutionNumbers with the same exposureNumber .  In Apr, 2020, there were 75
   * cases where this situation held.  In every case, there was no output for
   * the secondSolution number.  Set the second exposureNumber = 0 to -1 and
   * swap the two entries to keep the sort order correct.
   */

  for (outputRecordIndex = 0; outputRecordIndex < numOutputRecords;outputRecordIndex++) {
    pOutputTableEntry = &pOutputTable[outputRecordIndex];

    if (outputRecordIndex < (numOutputRecords -1)) {
      pOutputTableEntry2 = &pOutputTable[outputRecordIndex+1];
      if ((pOutputTableEntry->plateNumber == pOutputTableEntry2->plateNumber) &&
          (pOutputTableEntry->mosaicNumber == pOutputTableEntry2->mosaicNumber) &&
          (pOutputTableEntry->exposureNumber >= 0) &&
          (pOutputTableEntry->exposureNumber == pOutputTableEntry2->exposureNumber) &&
          (strcmp(pOutputTableEntry->series,pOutputTableEntry->series) == 0)) {

        if (verbose) {
          printf("line %4d WARNING: duplicate entries for %5s%04d_%02d s%02d exposure %02d\n",
                 __LINE__,
                 pOutputTableEntry->series,
                 pOutputTableEntry->plateNumber,
                 pOutputTableEntry->mosaicNumber,
                 pOutputTableEntry->solutionNumber,
                 pOutputTableEntry->exposureNumber);
        }

        pOutputTableEntry2->exposureNumber = -1;
        memcpy(pMosaicList,pOutputTableEntry2,sizeof(MOSAICLIST));
        memcpy(pOutputTableEntry2,pOutputTableEntry,sizeof(MOSAICLIST));
        memcpy(pOutputTableEntry,pMosaicList,sizeof(MOSAICLIST));
      }
    }
  }

  if ((numOutputRecords + numExposureRecords) == 0) {
    if (pMosaicTable != NULL) {
      free(pMosaicTable);
    }

    if (pExposureTable != NULL) {
      free(pExposureTable);
    }

    return 0;
  }

  if (numExposureRecords == 0) {
    *pMosaicListSize = numOutputRecords;
    *pTotalSolutionsCount = totalSolutionsCount;
    *pTotalSolution0Count = totalSolution0Count;
    *ppMosaicTable = pOutputTable;
    return 1;
  }

  /* Here we must merge the exposures table with the output table. */

  /* Free the output table by copying it to the mosaic table */
  for (outputRecordIndex = 0; outputRecordIndex < numOutputRecords;outputRecordIndex++) {
    pOutputTableEntry = &pOutputTable[outputRecordIndex];
    pMosaicTableEntry = &pMosaicTable[outputRecordIndex];
    memcpy(pMosaicTableEntry,pOutputTableEntry,sizeof(MOSAICLIST));
  }

  numMosaicRecords = numOutputRecords;
  outputRecordIndex = 0;
  mosaicRecordIndex = 0;
  pMosaicTableEntry = &pMosaicTable[mosaicRecordIndex];

  for (exposureRecordIndex = 0; exposureRecordIndex < numExposureRecords; exposureRecordIndex++) {
    pExposureTableEntry = &pExposureTable[exposureRecordIndex];

    while (1) {
      if ((numMosaicRecords <= 0) || (pMosaicTableEntry == NULL)) {
        break;
      }

      strcmpResult = strcmp(pMosaicTableEntry->series,pExposureTableEntry->series);

      if ((strcmpResult > 0) || (pMosaicTableEntry->plateNumber > pExposureTableEntry->plateNumber)) {
        break;
      }

      if ((strcmpResult == 0) &&
          (pMosaicTableEntry->plateNumber == pExposureTableEntry->plateNumber) &&
          (pMosaicTableEntry->exposureNumber > pExposureTableEntry->exposureNumber)) {
        break;
      }

      if ((strcmpResult < 0) ||
          (pMosaicTableEntry->plateNumber < pExposureTableEntry->plateNumber) ||
          (pMosaicTableEntry->exposureNumber < pExposureTableEntry->exposureNumber)) {

        if (pMosaicTableEntry->outputFlag == 0) {
          if ((strcmpResult == 0) &&
              (pMosaicTableEntry->plateNumber == pExposureTableEntry->plateNumber)) {
            /* Here, everything but the exposure number matches */
            strcpy(pMosaicTableEntry->date,pExposureTableEntry->date);
            strcpy(pMosaicTableEntry->plateClass,pExposureTableEntry->plateClass);
            strcpy(pMosaicTableEntry->exposure,pExposureTableEntry->exposure);
          }

          pMosaicTableEntry->outputFlag = 1;
          pOutputTableEntry = &pOutputTable[outputRecordIndex];
          memcpy(pOutputTableEntry,pMosaicTableEntry,sizeof(MOSAICLIST));
          outputRecordIndex++;
        }

        mosaicRecordIndex++;

        if (mosaicRecordIndex < numMosaicRecords) {
          pMosaicTableEntry = &pMosaicTable[mosaicRecordIndex];
          continue;
        } else {
          /* This no-effect statement is surely a mistake but fixing it yields a segfault below! */
          /*pMosaicTableEntry == NULL;*/
          break;
        }
      }

      break;
    }

    if (outputRecordIndex > numOutputAlloc) {
      printf("ERROR: pipelineutils line %4d exposureRecordIndex %d %d exceeded\n",__LINE__,exposureRecordIndex,numExposureRecords);
      exit(1);
    }

    if (
      (pMosaicTableEntry->plateNumber == pExposureTableEntry->plateNumber) &&
      (pMosaicTableEntry->exposureNumber == pExposureTableEntry->exposureNumber) &&
      (strcmp(pMosaicTableEntry->series,pExposureTableEntry->series) == 0)
    ) {
      /* We are going to output a mosaic record */
      if (pMosaicTableEntry->outputFlag == 0) {
        pMosaicTableEntry->outputFlag = 1;
        pOutputTableEntry = &pOutputTable[outputRecordIndex];
        memcpy(pOutputTableEntry,pMosaicTableEntry,sizeof(MOSAICLIST));
        outputRecordIndex++;

        pOutputTableEntry->exposureNumber =   pExposureTableEntry->exposureNumber;
        pOutputTableEntry->dRightAscension =   pExposureTableEntry->dRightAscension;
        pOutputTableEntry->dDeclination =   pExposureTableEntry->dDeclination;
        pOutputTableEntry->centerSource =  pExposureTableEntry->centerSource;

        strcpy(pOutputTableEntry->date,pExposureTableEntry->date);
        strcpy(pOutputTableEntry->plateClass,pExposureTableEntry->plateClass);
        strncpy(pOutputTableEntry->exposure,pExposureTableEntry->exposure,sizeof(pOutputTableEntry->exposure) - 1);
        strcpy(pOutputTableEntry->rightAscension,pExposureTableEntry->rightAscension);
        strcpy(pOutputTableEntry->declination,pExposureTableEntry->declination);

        if (verbose) {
          if ((pOutputTableEntry->centerSource == CENTERSOURCE_IMWCS) &&
              (pOutputTableEntry->WCSSource != WCSSOURCE_IMWCS)) {
            printf("line %4d WARNING: stale WCSFIT location for  %5s%04d_%02d s%02d exposure %02d WCSSource %d centerSource %d\n",
                   __LINE__,
                   pOutputTableEntry->series,
                   pOutputTableEntry->plateNumber,
                   pOutputTableEntry->mosaicNumber,
                   pOutputTableEntry->solutionNumber,
                   pOutputTableEntry->exposureNumber,
                   pOutputTableEntry->WCSSource,
                   pOutputTableEntry->centerSource);
            if (pOutputTableEntry->WCSSource == WCSSOURCE_LOGBOOK) {
              printf("SELECT series,plateNumber,exposureNumber,centerSource,version FROM exposures INNER JOIN logbook using (series,plateNumber,exposureNumber) where series = '%s' and plateNumber = %d and exposureNumber = %d;\n",
                     pOutputTableEntry->series,
                     pOutputTableEntry->plateNumber,
                     pOutputTableEntry->exposureNumber);
              printf("UPDATE logbook set version = '%s' where series = '%s' and plateNumber = %d and exposureNumber = %d;\n",
                     LOGBOOK_VERSION,
                     pOutputTableEntry->series,
                     pOutputTableEntry->plateNumber,
                     pOutputTableEntry->exposureNumber);
              printf("UPDATE exposures set centerSource = 'Logbook' where series = '%s' and plateNumber = %d and exposureNumber = %d;\n",
                     pOutputTableEntry->series,
                     pOutputTableEntry->plateNumber,
                     pOutputTableEntry->exposureNumber);
            }
          }
        }
      }
    } else {
      /* Output the exposure record */
      pOutputTableEntry = &pOutputTable[outputRecordIndex];
      memcpy(pOutputTableEntry,pExposureTableEntry,sizeof(MOSAICLIST));
      outputRecordIndex++;
    }
  }

  numOutputRecords = outputRecordIndex;

  if (pMosaicTable != NULL) {
    free(pMosaicTable);
  }

  if (pExposureTable != NULL) {
    free(pExposureTable);
  }

  *pMosaicListSize = numOutputRecords;
  *pTotalSolutionsCount = totalSolutionsCount;
  *pTotalSolution0Count = totalSolution0Count;
  *ppMosaicTable = pOutputTable;
  return 1;
}
