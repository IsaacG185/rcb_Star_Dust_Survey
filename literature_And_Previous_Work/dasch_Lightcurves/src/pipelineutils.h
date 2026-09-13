// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* pipelineutils.h
 *
 * Jan 16, 2008 Edward J. Los - Initial Version
 * Jul  1, 2008 Edward J. Los - add FWHM_WORLD_RADIUS_FACTOR and PIXEL_RADIUS_FACTOR
 *                              to set the GSC and blended object search radius
 *                              add BAD_FLAGS_LIMIT
 * Jul  8, 2008 Edward J. Los - Add FILTER_FLAG_HIZOUT for a local smoothing quality check
 * Jul 28, 2008 Edward J. Los - Add MINIMUM_ZOUT and MAXIMUM_ZOUT
 *                              Add FILTER_FLAG_TOO_BRIGHT
 * Aug  6, 2008 Edward J. Los - Add MAX_DRAD_PIXELS
 * Aug  8, 2008 Edward J. Los - Tie PIXEL_RADIUS_FACTOR to MAX_DRAD_PIXELS (decreasing it from 5 to 2)
 * Aug 11, 2008 Edward J. Los - Add FILTER_FLAG_ADJUST_BLEND flag
 *                              Add zeroFlag to CalcMedianAndRMS
 *                              Add DMAG_REJECT flags
 * Aug 19, 2008 Edward J. Los - Add MAX_LOCAL_RMS
 * Aug 22, 2008 Edward J. Los - Add THRESHOLD_FACTOR
 * Aug 26, 2008 Edward J. Los - Add FILTER_FLAG_BINDRAD
 * Sep 29, 2008 Edward J. Los - Add CCMAP_THRESHOLD
 * Oct 11, 2008 Edward J. Los - Reduce CCMAP_THRESHOLD to THRESHOLD_FACTOR
 * Oct 20, 2008 Edward J. Los - Add duplicate GSC2.3 CLASS flag; Add CLASS_MASK
 * Oct 21, 2008 Edward J. Los - Change MAXIMUM_ZOUT from 1.0 to 0.5
 *                              Add MAXIMUM_ERROUT
 * Nov 25, 2008 Edward J. Los - Replace MAX_DRAD_PIXELS with a series dependent function
 * Dec  8, 2008 Edward J. Los - Add reject reasons from extract_lightcurves
 * Dec 22, 2008 Edward J. Los - Add FILTER_FLAG_LOW_ALTITUDE; Add extinction parameters
 * Dec 27, 2008 Edward J. Los - Increas MAX_REF to 42 to handle new HCO designations.
 * Jan  2, 2009 Edward J. Los - Add WCSSource definitions, increase MAX_QUERY_STRING
 * Jan 22, 2009 Edward J. Los - Add RaPM and DecPM, the proper motions in mas/yr.
 * Jan 28, 2009 Edward J. Los - Move the CAL_FLAGS (now BFLAGS field) here
 *                              Split FLAGS into AFLAGS and BFLAGS
 *                              Add FILTER_BFLAG_BAD_JD and FILTER_BFLAG_PROPER_MOTION
 * Feb  6, 2009 Edward J. Los - Replace FILTER_BFLAG_BLEND with FILTER_AFLAG_BLEND
 *                              Add NEW_BLEND_DEFINITION
 * Feb 10, 2009 Edward J. Los - Add MAX_CATALOG_EXPOSURES
 * Feb 16, 2009 Edward J. Los - Add CheckBlend and rotateFrame
 *                              Add X_DMAGBINS and Y_DMAGBINS
 *                              Add FILTER_BFLAG_DRAD_ADJUST
 *                              Add DRAD_REJECT_FACTOR
 *                              Add MAX_ISO_RMS
 * Feb 24, 2009 Edward J. Los - Add DecodeAFLAGS
 * Feb 27, 2009 Edward J. Los - Add CalcMemory
 * Mar  3, 2009 Edward J. Los - Add MIN_GOODSTARS
 * Mar  9, 2009 Edward J. Los - Add ComputeRange
 * Mar 30, 2009 Edward J. Los - Add COLOR_RATIO
 * Apr  4, 2009 Edward J. Los - Flag objects without colorterm correction in plots
 * May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 * Jun 15, 2009 Edward J. Los - Add binned MAG_ISO vs Stdmag for initial magnitude filtering
 * Jun 29, 2009 Edward J. Los - Back out MAG_ISO_REJECT_FACTOR
 *                              Increase MAX_QUERY_STRING to 3000 to handle InsertNewMagnitude
 * Jun 30, 2009 Edward J. Los - Add REMOVE_NOMATCH to turn off FILTER_AFLAG_BLEND_NOMATCH functionality
 * Jul 20, 2009 Edward J. Los - Add ChargeTime for performance tracking
 * Aug 17, 2009 Edward J. Los - Add Kepler Input Catalog support
 * Aug 31, 2009 Edward J. Los - Add solutionNumber
 *                              Move ID conversion routines here
 * Nov 16, 2009 Edward J. Los - Add UNCERTAIN_DATE and MULTIPLE_BLEND flags and EXT_TIME_TOLERANCE
 * Nov 23, 2009 Edward J. Los - Add UNKNOWN_EXPOSURE
 * Dec  4, 2009 Edward J. Los - Add FILTER_AFLAG_MULTIPLE_NONE
 * Dec 11, 2009 Edward J. Los - Add new color transformation support
 * Dec 15, 2009 Edward J. Los - Add new FitWCS flags
 * Dec 15, 2009 Edward J. Los - Add SetMosaicFitWCS to keep track of pipeline failures
 * Dec 22, 2009 Edward J. Los - Add SetPlateQuality to keep track of second-quality plates
 * Jan 19, 2010 Edward J. Los - Add elevation to the locations table, PLATECOLOR and GetPlateColor();
 * Jan 22, 2010 Edward J. Los - Add plate quality support: FILTER_AFLAG_QUALITY, FILTER_BFLAG_COLORTERM,
 *                              REJECT_REASON_QUALITY, and SYMBOL_QUALITY.  Sumin Tang memorandum
 *                              of 12/11/09 1:07 AM
 * Feb 23, 2009 Edward J. Los - Flag good stars for with sextractor blend or neighbor flags have been set
 * Apr 12, 2010 Edward J. Los - Add GetSeriesInfo and GetPlateInfo, location name for VOTABLE support.
 * Apr 16, 2010 Edward J. Los - Add option to exclude series with known fitting problems
 * May  4, 2010 Edward J. Los - Add GetMysqlError()
 * May 11, 2010 Edward J. Los - Add REF_TYPE fields and GetDSACHCoordinates
 * May 18, 2010 Edward J. Los - Add FILTER_AMASK_PLOT for blends, defects, low altitude, undertain date, and second quality plates.
 * Jun  4, 2010 Edward J. Los - Move flagsentry here from showflags
 * Jul 19, 2010 Edward J. Los - Add DEFAULT_TIME_ACCURACY for exposures with a null timeAccuracy (usually card catalog entries).
 * Aug 24, 2010 Edward J. Los - Add GSC_VARIABLE_BIT to the default list of plotted points
 *                              Add FILTER_AFLAG_MULTIPLE_BLEND to FILTER_AMASK_BLEND and FILTER_AMASK_PLOT
 * Oct 26, 2010 Edward J. Los - Add InaccuratePV flag
 * Dec 21, 2010 Edward J. Los - Remove the grating plates from the ALL_QUALITY_MASK
 * Mar  7, 2011 Edward J. Los - Move catalog definitions here from photometryutils.h
 *                              Add color limits for the catalogs.
 * Mar 23, 2011 Edward J. Los - Add AAVSO APASS support
 * Apr  8, 2011 Edward J. Los - Add GetREFType
 * Jun 27, 2011 Edward J. Los - Add SelectBestMosaic
 * Jul 18, 2011 Edward J. Los - Change THUMBNAIL_SIZEUNITS to THUMBNAIL_PIXELS
 * Aug 10, 2011 Edward J. Los - Add fields from the "flatfields" table
 * Aug 19, 2011 Edward J. Los - Add Tycho2 reference numbers
 * Sep  7, 2011 Edward J. Los - Back out deprecated COLOR_RATIO
 * Sep 17, 2011 Edward J. Los - Expand APASS color range for B and V
 * Sep 28, 2011 Edward J. Los - Add FILTER_BFLAG_PSFSATURATED
 * Oct 30, 2011 Edward J. Los - Add magnitude-dependent star correction support (FILTER_BFLAG_MAGDEP_MAGCOR)
 * Nov 25, 2011 Edward J. Los - Add MAGDEP_V2 changes described in Sumin Tang memorandum of Fri 11/18/11 5:27 PM
 * Nov 28, 2011 Edward J. Los - Add MAGDEP_V3 changes described in Sumin Tang memorandum of Mon 11/28/11 1:47 PM
 *                              Add SYMBOL_SATURATED and REJECT_REASON_SATURATED
 * Dec  5, 2011 Edward J. Los - Increase THUMBNAIL_IMAGESIZE to 300 in conjunction with reducing scaleFactor in web_plot.c.
 * Dec 10, 2011 Edward J. Los - Add SYMBOL_NOMAGDEP and REJECT_REASON_NOMAGDEP and QUALITY_NOMAGDEP
 * Dec 16, 2011 Edward J. Los - Increase MAX_LIMITING_MAG from 0.5 to 0.75 in accordance with Sumin Tang's memorandum of 12/15/11 1:55 PM
 * Dec 19, 2011 Edward J. Los - Add THUMBNAIL_SCALED_PIXELS
 * Feb 17, 2012 Edward J. Los - Add multiple exposure mask support
 * Feb 24, 2012 Edward J. Los - Add UpdateQuality
 * Mar 13, 2012 Edward J. Los - Support password-based authentication
 * Mar 23, 2012 Edward J. Los - Support plot selection by series
 * Apr 25, 2012 Edward J. Los - Add lightcurve estimation symbols
 * Jun 29, 2012 Edward J. Los - Add GetExposureInfo
 * Jul  3, 2012 Edward J. Los - Add high background object support
 * Jul  9, 2012 Edward J. Los - Remove THRESHOLD_FACTOR for the new defect filter of  Fri 6/29/12 9:43 AM
 * Nov  9, 2012 Edward J. Los - Revise the colorterm limits
 * Dec 19, 2012 Edward J. Los - Define FILTER_AMASK_PLOT2 to be ((~ ((1<<FILTER_AFLAG_ISO_RMS) | (1<<FILTER_AFLAG_LOCAL_RMS))) & 0x7FFFFFFF)
 * Dec 26, 2012 Edward J. Los - Support nearby galaxy and variable star finder.
 * Feb 11, 2013 Edward J. Los - Add FILTER_AFLAG_LIMITING_MAG to FILTER_AMASK_PLOT2
 * Feb 15, 2013 Edward J. Los - Reduce MAX_LIMITING_MAG to 0.5 mag (See Josh Grindlay Tue 2/12/13 9:13 PM mail)
 *                              Increase BIN_EXPANSION_RADIUS from 0.010 to 1/64 degree (0.0156 degrees)
 * Feb 25, 2013 Edward J. Los   Add support for running astrometry.net on full sized mosaics
 * Apr 22, 2013 Edward J. Los   Add release levels
 * Jun 18, 2013 Edward J. Los   Add the background flag from FILTER_AMASK_PLOT2
 * Jul 16, 2013 Edward J. Los   Introduce pseudo-quality QUALITY_TRAILED for ELLIPTICITY > MAX_ELLIPTICITY = 0.6
 * Jul 23, 2013 Edward J. Los   Define MIN_SCAMP_ASTROMETRY_SCALE
 * Aug  3, 2013 Edward J. Los   Add FILTER_MASK_RMS_REJECT for the calculation of rms values of stars near the limiting magnitude
 * Aug 12, 2013 Edward J. Los   Add proper motion standard deviations to GSCIMAGE
 *                              Reformat GSCIMAGE, eliminating GSCIMAGEX; Add UCAC4 reference type
 * Apr 15, 2014 Edward J. Los   Introduce GSCAPASSIMAGE with new columns for APASS colors gmag,rmag,imag and errors vmagerr,bmagerr,gmagerror,rmagerr,imagerr
 * Jun  3, 2014 Edward J. Los   Support DR3
 * Sep 30, 2014 Edward J. Los   Add LINEARITY_DELAY_USEC from /apps/FitsWrapper/Catalog.h
 * Oct 25, 2014 Edward J. Los   Add FILTER_AFLAG_WEDGE to FILTER_AMASK_LOWDRAD, FILTER_AMASK_PLOT and FILTER_AMASK_PLOT2
 * Nov 10, 2014 Edward J. Los   Add FILTER_AMASK_NONE_CANDIDATE2 which adds FILTER_AFLAG_LIMITING_MAG to FILTER_AMASK_LOWDRAD
 * Dec 22, 2014 Edward J. Los   Support DR4
 * May 16, 2015 Edward J. Los   Add formatAflags to list individual bits set in AFLAGS
 * Sep 28, 2015 Edward J. Los   Add EXCLUSION to study transient candidates favoring particular regions of the sky
 * Oct  5, 2015 Edward J. Los   Add local_strlwr
 * Oct 16, 2015 Edward J. Los   Add FITWCS_NOAPASSALLOBJECTS
 * Nov 17, 2015 Edward J. Los   Support DR5
 * Feb 24, 2016 Edward J. Los - Change "formatAflags" to "FormatFlagsBits" for use with BFLAGS and quality bitmaps
 * Mar  1, 2016 Edward J. Los - Define a QUALITY_UNITIALIZED field for invalid quality fields.
 * Mar  8, 2016 Edward J. Los - Use two separate fields for FormatFlagsBits
 * Apr 22, 2016 Edward J. Los - Increase MAX_SERIES from 90 to 100
 * Sep 30, 2016 Edward J. Los - define REJECTBIT for REJECT_REASON text; add FATAL_REASON
 * Nov 25, 2016 Edward J. Los - support GAIA catalog
 * Dec 23, 2016 Edward J. Los - Add BlueColorterm()
 * Feb 22, 2017 Edward J. Los   Add INCLUDE_TRANSIENT_LIMITING_MAG_BACKGROUND to include limiting magnitudes and high background in flare searches
 * Feb 28, 2017 Edward J. Los   Add MAX_QUALIFIER
 * Mar 15, 2017 Edward J. Los   Rename INCLUDE_TRANSIENT_LIMITING_MAG_BACKGROUND to INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND
 * Mar 28, 2017 Edward J. Los   Add peakEvaluation bitmap to indicate reasons for rejection of transient candidate flares
 * Apr  4, 2017 Edward J. Los   Split PEAKEVALUATION_DATERANGE into PEAKEVALUATION_TOOSHORTDATERANGE, PEAKEVALUATION_TOOLONGDATERANGE, and PEAKEVALUATION_PEAKEVALUATION_INSUFFICIENT_POINTS3
 *                              Add PEAKEVALUATION_EXCLUSIONZONE, PEAKEVALUATION_INSUFFICIENT_POINTS2, and PEAKEVALUATION_TOOMANYPOINTS
 * Apr 25, 2017 Edward J. Los   Move WEDGEENTRY from filter_wedge for sharing between filter_wedge.c and update_quality.c
 * Jun 27, 2017 Edward J. Los   Add kernelType to MULTMASK
 * Jul 25, 2017 Edward J. Los   Add SMALL_PLATE_PIXELS to handle processing for 5"x7" or smaller plates
 *                              Redefine X_DMAGBINS to X_DMAGBINS_NORMAL and Y_DMAGBINS to Y_DMAGBINS_NORMAL
 *                              Add SHORT_DMAGBINS_SMALL and LONG_DMAG_BINS_SMALL for small plates
 *                              Add XDmagBins() and YDmagBins() for dynamic bin calculation.
 *                              Move X_EXTBINS and Y_EXTBINS from prepare_octave.c here and link them with X_DMAGBINS_NORMAL and Y_DMAGBINS_NORMAL
 *                              GetMosaicInfo() will now return the number of smoothing bins
 * Jan 22, 2018 Edward J. Los   Add RELEASE_EXPERIMENTAL
 * Jan 29, 2018 Edward J. Los   Add FILTER_AFLAG_DRADBIN to FILTER_AMASK_PLOT2
 * Feb 19, 2018 Edward J. Los   Support the full GAIA format by using the high bit of REFNumber
 * May 14, 2018 Edward J. Los   Add ReadLinearityFile() from ./FitsWrapper/ScannerCCD.cpp
 * May 21, 2018 Edward J. Los   Add REF_TYPE_GAIA2 for Gaia DR2 (source_id numbers are not guaranteed to be unique across releases)
 * May 25, 2018 Edward J. Los   Add CATALOG_GAIA to the catalog choices
 * Jun 11, 2018 Edward J. Los   Increase MAX_CATALOG_EXPOSURES from 33 to 52
 * Jul  9, 2018 Edward J. Los   Define MIN_METEOR_ASTROMETRY_SCALE to separate patrol from meteor plates.
 * Jul 23, 2018 Edward J. Los   Increase MAX_SERIES from 100 to 127 (after this the binary file format must change).
 * Oct 23, 2018 Edward J. Los   Add REF_TYPE_ATLAS2 for the Atlas 2 catalog
 * Nov 23, 2018 Edward J. Los   Add the floodJacket field to the plates structure
 * Nov 27, 2018 Edward J. Los   Define ATLAS refcat2 and GAIA dr2 colorterms
 * Dec 10, 2018 Edward J. Los   Define RELEASE_LEVEL_DR7
 * Jul 12, 2019 Edward J. Los   Add GetPlateCondition and GetPlateEvent for stacks database access
 * Oct  4, 2019 Edward J. Los   Add remaining release levels
 * Nov 29, 2019 Edward J. Los   Add GetFlatfieldRecords
 * Dec  3, 2019 Edward J. Los   Add scanDate to GetMosaicInfo()
 *                              Define MAX_SCAN_MOSAIC_MOSAIC_NUMBER for GetMosaicInfo search loops in runpipeline.c
 * Feb  9, 2020 Edward J. Los   Add SCAN_PATTERN... from Catalog.h
 * Mar 19, 2020 Edward J. Los   Add patternID to GetMosaicInfo() and CCD_PIXEL_SIZE (must match ~/apps/Catalog.h)
 * Mar 22, 2020 Edward J. Los   Move RESULTTABLE to pipelineutils and rename it to MOSAICLIST
 *                              Move SEQUESTERED_YES and SEQUESTERED_NO here from photometryutils.h and scanread.h
 *                              Move WCSSOURCE definitions from scanread.h
 *                              Add CENTERSOURCE definitions
 *                              Add orientation and fittedPlateScale to SERIES_HEADER
 *                              Add definitions  DEBUG_BUILDMOSAICLIST and DASCH_NOTAVAILABLE for the "N.A." string
 * Apr  5, 2020 Edward J. Los   Replace MOSAIC with MOSAICLIST, leaving the latter as an alias
 * Apr  9, 2020 Edward J. Los   Add mosaics.JobId to MOSAIC (MOSAICLIST) for compatabilty with Catalog.h
 *                              Rename class to plateClass to avoid cpp errors
 * May 25, 2020 Edward J. Los   Add LOGBOOK_VERSION for dasch_scat consistency checks
 * Jan 26, 2021 Edward J. Los   Define PEAKEVALUATION_PRIMARY_GSC_BIN and PEAKEVALUATION_SECONDARY_GSC_BIN to
 *                               support display of search results for Josh's list of flare candidates
 */

#ifndef _PIPELINEUTILS_H
#define _PIPELINEUTILS_H 1
#define MAX_SERIES_STRING 20
#define MAX_COMMENT_STRING 132
#define MAX_QUERY_STRING 3000
#define MAX_SERIES 127
#define MAX_CLASS_STRING 8
#define MAX_CENTERSOURCE_STRING 10
#define MAX_RIGHTASCENSION_STRING 13
#define MAX_DECLINATION_STRING 16
#define MAX_DATE_STRING 25
#define MAX_FILENAME_LEN 256
#define MAX_CTYPE_STRING 10
#define MAX_EXPOSURE_STRING 8
#define MAX_NOTES_STRING 80
#define MAX_DESCRIPTION_STRING 80
#define MAX_LOCATION_NAME 80
#define MAX_QUALITY_STRING 64
#define MAX_QUALIFIER 15 /* Maximum size of char *catalogText string in pipelineutils.c */
#define PI_VALUE 3.141592654
#define RAD_TO_DEGREES  (180.0/3.141592654)
#define DEGREES_TO_RAD  (3.141592654/180.0)
#define NOMINAL_MM_PER_PIXEL   0.011 /* nominal mm/pixel  Results of May 22, 2006: 0.010988 mm/pixel */
#define FWHM_WORLD_RADIUS_FACTOR 0.5
#define DRAD_REJECT_FACTOR 3.0
#define MAX_PLATE_NUMBER 100000
#define MAX_BRIGHT_ADJUST 0.01 /* Adjust the maximum bright magnitude by this amount */
#define PATROL_PLATE_SCALE (350.*NOMINAL_MM_PER_PIXEL/3600.0) /* patrol plate scale limit in degrees/pixel.  Plates with higher scales are patrol plates */
#define MIN_SCAMP_ASTROMETRY_SCALE   (400 * NOMINAL_MM_PER_PIXEL)  /* Second pass astrometry candidates.  All patrol plates except for rb/rh series.  arcsec/pixel */
#define MIN_PATROL_ASTROMETRY_SCALE   (400 * NOMINAL_MM_PER_PIXEL)  /* Patrol and meteor plates except for rb/rh arcsec/pixel (must agree with ScannerCCD.h)
                                                                       This value for astrometry fitting and excludes rb/rh plates while PATROL_PLATE_SCALE is for
                                                                       the website plotter and included rb/rh plates */
#define CCD_PIXEL_SIZE    4096 /* width and height of the CCD in pixels.  Must match value in ~/apps/Catalog.h */
#define MIN_METEOR_ASTROMETRY_SCALE  (900 * NOMINAL_MM_PER_PIXEL)  /* Meteor plates arcsec/pixel (must agree with ScannerCCD.h) */
#define MAX_REGIONFLAG 10
#define MAX_NEARBY_OBJECTS_STRING 175 /* Last increased Jun 16, 2015 */
#define MATCH_RADIUS_BINS 4 /* Used by matchstars and filterblended to set the maximum match radius - equivalent to 37 pixels for 600 arcsec/mm patrol plates */
#define LINEARITY_DELAY_USEC 35 /* Default delay used by the linearity file processing */
#define SMALL_PLATE_PIXELS  200000000 /* Product of naxis1 and naxis2. (<=  95480580 for the 4"x5" lwla plates and <= 171704836 for the 5"x7" hale plates) */
#define LARGE_PLATE_PIXELS 1000000000 /* Product of naxis1 and naxis2 for 14" x 17" plates  (also in DiskUtilization.h) */
#define MAX_SCAN_MOSAIC_MOSAIC_NUMBER 40 /* Used by runpipeline loop to check all possible mosaics */
/* #define DEBUG_BUILDMOSAICLIST 1 *//* Give warnings of routines which need to be modified. */
/* #define DEBUG_NOTAVAILABLE 1 */   /* If defined, resorts to not using the "N.A." string */
#define DASCH_NOTAVAILABLE "    N.A.    "

/* The reference type is the leading digit of the reference designation */
#define REF_TYPE_NONE   0
#define REF_TYPE_GSC    1
#define REF_TYPE_KEPLER 2
#define REF_TYPE_DASCH  3
#define REF_TYPE_APASS  4
#define REF_TYPE_TYCHO2 5
#define REF_TYPE_UCAC4  6
#define REF_TYPE_GAIA1  7  /* Gaia DR1 */
#define REF_TYPE_GAIA2  8  /* Gaia DR2 */
#define REF_TYPE_ATLAS2 9  /* Atlas refcat 2 */

#define GAIA_MODULUS 128  /* Divide GAIA numbers by this to get at least one extra decimal place */
/* Note: the following bit numbers are 1 based! */
#define GAIA_HEALPIX   0x7ffffff800000000L  /* HEALPix number, bits 36-63  (28 bits) */
#define GAIA_PROCCTR   0x0000000600000000L  /* Data processing center, bits 34-35 (2 bits) */
#define GAIA_SPARE     0x0000000100000000L  /* spare bit bit 33 (1 bit) */
#define GAIA_RUNSEQ    0x00000000ffffff80L  /* running sequence, bits 8-32 (25 bits) */
#define GAIA_COMPONENT 0x000000000000007fL  /* Component number bits 1-7 (7 bits) */

#define GAIA_HEALPIX_BIT   35
#define GAIA_PROCCTR_BIT   33
#define GAIA_SPARE_BIT     32
#define GAIA_RUNSEQ_BIT     7
#define GAIA_COMPONENT_BIT  0

#define GAIA_HEALPIX_MASK   0xfffffffL  /* 28 bits */
#define GAIA_PROCCTR_MASK   0x3L        /*  2 bits */
#define GAIA_SPARE_MASK     0x1L        /*  1 bit  */
#define GAIA_RUNSEQ_MASK    0x1ffffffL  /* 25 bits */
#define GAIA_COMPONENT_MASK 0x7fL       /*  7 bits */

#define MINIMUM_ZOUT -0.5
#define MAXIMUM_ZOUT  0.5
#define MAXIMUM_ERROUT 0.7
#define MAX_LIMITING_MAG 0.5
#define MAX_LOCAL_RMS 1.0
#define MAX_ISO_RMS 1.0

#define MAX_SRC_LENGTH 80 /* The maximum name assigned to a lightcurve, including blends */
#define MAX_BLEND_COUNT 3 /* The maximum number of blends combinded for a lightcurve */

#define MIN_GOODSTARS 10
#define MIN_SUMMARY_GOODSTARS 2
#define REMOVE_NOMATCH 1

#define COLORFLAG_NONE 0
#define COLORFLAG_LINEAR 1
#define COLORFLAG_METROPOLIS 2
#define COLORFLAG_COPY 3
#define COLORFLAG_MAX 4

/* Size of thumbnail for web plotting */
#define THUMBNAIL_IMAGESIZE 300 /* This is a radius */
#define THUMBNAIL_PIXELS    "pixels"
#define THUMBNAIL_ARCSEC    "arcsec"
#define THUMBNAIL_MAX_PIXELS 2400
#define THUMBNAIL_SCALED_PIXELS 300 /* Size of the thumbnail when variable "arcsec" is requested */

/* Defines for the ChargeTime routine */
#define MAX_CHARGE_ENTRY 7
#define RESET_DELTA_ENTRY -1
#define RESET_ALL_ENTRY -2

#define AUTHORIZE_ANGLE 5.0 /* Degrees around point of interest for searches */
#define AUTHORIZE_RADIUS 600.0 /* Largest radius for extracted images */

/* These are the major partitions of the database */
#define CATALOG_GSC232       0
#define CATALOG_KEPLER       1
#define CATALOG_APASS        2
#define CATALOG_GAIA         3
#define CATALOG_ATLAS        4
#define CATALOG_EXPERIMENTAL 5
#define MAX_CATALOG_NUMBER   6

extern char *catalogText[MAX_CATALOG_NUMBER];


#define RELEASE_EXPERIMENTAL 1  /* If 1, then make the experimental calibration available; if 0, do not */
#define RELEASE_EXPERIMENTAL_CAT0 CATALOG_GSC232
#define RELEASE_EXPERIMENTAL_CAT1 CATALOG_APASS
#define RELEASE_LEVEL RELEASE_LEVEL_ALL /* Release level */
#define RELEASE_LEVEL_M44   0 /*   M44 only */
#define RELEASE_LEVEL_DR1   1 /*   DR1 and calibration fields */
#define RELEASE_LEVEL_DR2   2 /*   DR2 */
#define RELEASE_LEVEL_DR3   3 /*   DR3 */
#define RELEASE_LEVEL_DR4   4 /*   DR4 */
#define RELEASE_LEVEL_DR5   5 /*   DR5 */
#define RELEASE_LEVEL_DR6   6 /*   DR6 */
#define RELEASE_LEVEL_DR7   7 /*   DR7 */
#define RELEASE_LEVEL_DR8   8 /*   DR8 */
#define RELEASE_LEVEL_DR9   9 /*   DR9 */
#define RELEASE_LEVEL_DR10 10 /*  DR10 */
#define RELEASE_LEVEL_DR11 11 /*  DR11 */
#define RELEASE_LEVEL_DR12 12 /*  DR12 */
#define RELEASE_LEVEL_ALL 13 /* Everything */
#define RELEASE_FIELD_SEARCH_LIMIT 10 /* Limit lightcurve searches to 10 fields */

#define RELEASE_FIELD_DR1     0
#define RELEASE_FIELD_DR2     1
#define RELEASE_FIELD_DR3     2
#define RELEASE_FIELD_DR4     3
#define RELEASE_FIELD_DR5     4
#define RELEASE_FIELD_DR6     5
#define RELEASE_FIELD_DR7     6
#define RELEASE_FIELD_DR8     7
#define RELEASE_FIELD_DR9     8
#define RELEASE_FIELD_DR10    9
#define RELEASE_FIELD_DR11   10
#define RELEASE_FIELD_DR12   11
#define RELEASE_FIELD_M44    12
#define RELEASE_FIELD_3C273  13
#define RELEASE_FIELD_BAADE  14
#define RELEASE_FIELD_KEPLER 15
#define RELEASE_FIELD_LMC    16
#define RELEASE_FIELD_OTHER  17
#define RELEASE_FIELD_MAX    18
/*
 * The GSC2.3.2 catalog uses several different color systems:
 *    JpgMag  385-540nm (462 nm ave)
 *    FpgMag  610-690nm (650 nm ave)
 *    B        440 nm  (From Henden & Kaitchuck "Astronomical Photometry")
 *    V        550 nm   same
 *    R        700 nm   same
 *
 *   Since FpgMag-JpgMag is a 188 nm span and B-V is a 110 nm span, the COLOR_RATIO is
 *   approximately (FpgMag-JpgMag/(B-V) = 1.7.  Because B is within the JpgMag range,
 *   apply COLOR_RATIO only when the V magnitude bit is set.
 *
 *   A survey of SDSS stars suggested a COLOR_RATIO of 0.9
 *
 *   A series of experiments showed that a COLOR_RATIO is a minimum at 0.6.
 *
 */
/* #define COLOR_RATIO 0.6 This should agree with colorterm.m and annular9.m*/

#define MAX_CATALOG_EXPOSURES 52 /* WARNING - change CFitsExposure::GetExposureIndex if more are needed */

#define UNKNOWN_EXPOSURE -250 /* This must fit in a signed char */

#define CCMAP_THRESHOLD  2.3
/* #define EXTINCTION_COEFF 0.18 */ /* Magnitudes per airmass - Deprecated. Use GetExtinctionCoefficient*/
#define MAX_AIRMASS 2.5 /* Maximum allowable airmass */


#define NEW_BLEND_DEFINITION 1 /* New definition of image blending */

#define MIN_BLEND_MAGNITUDE 0.1 /* Minimum acceptable difference between brightest star and its blend */
#define SEXTRACTOR_NISO 8      /* Number of sextractor ISO levels */

/* This produces 2500 bins for normal plates.  If more bins are needed, check all instances in the code to be sure
 *  sufficient space is allocated.  Must agree with Catalog.h
 */
#define X_DMAGBINS_NORMAL 50                  /* local calibration bins */
#define Y_DMAGBINS_NORMAL X_DMAGBINS_NORMAL   /* local calibration bins */
/* Extinction bins must match local calibration bins because both are included in the localbins photometry tables */
#define X_EXTBINS_NORMAL X_DMAGBINS_NORMAL /* X-axis extinction bins */
#define Y_EXTBINS_NORMAL Y_DMAGBINS_NORMAL /* Y-axis extinction bins */

/*
 * For small plates 20x25 binning produces 500 bins or 20% of normal bins for 25% plate area and 16x20
 * produces 320 bins or 12.8% of normal bins.  These changes make most bins square.
 * Because the small plates are usually for longer focal length scopes, we can expect fewer
 * stars per plate area.  Selection of the 16x20 pattern also coincides with a proposed rewrite of the
 * search_close correlation algorithm to investigate the changes in multiple exposure correlation as
 * a function of position on the plates
 */
#define LONG_DMAGBINS_SMALL  20       /* local calibration bins */
#define SHORT_DMAGBINS_SMALL 16           /* local calibration bins */
#define STARS_PER_DMAGBIN   20  /* Minimum accepted stars in a local calibration bin */

#define TOTAL_DMAGBINS_NORMAL (X_DMAGBINS_NORMAL*Y_DMAGBINS_NORMAL)
#define TOTAL_DMAGBINS_SMALL  (SHORT_DMAGBINS_SMALL*LONG_DMAGBINS_SMALL)

/* Of the 57 closest stars in the RASC "Observers Handbook 2000"
 * Barnards star has the highest proper motion of 10358 mas/year.
 * Since the earliest plate dates to 1870, the maximum movement
 * from 1870 to 2000 is 10358*(2000-1870)/(1000*3600) = 0.374
 * degrees.  After applying the proper motion correction, we will
 * need to expand the search area by this amount. */
#define PIPELINE_MIN_DATE 1870
#define PIPELINE_MAX_DATE 2000

#define PROPER_MOTION_DEGREES 0.374 /* what is this magic number??? */
#define PROPER_MOTION_SIGMA 2.0  /* Number of standard deviations for the proper motion error */
#define PROPER_MOTION_ERROR_LIMIT 1.0  /* Limit in the proper motion error in arcsec for the FILTER_BFLAG_PMERROR bit to be set */
#define INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND 1
/*  INCLUDE_TRANSIENT_LOCAL_RMS_BACKGROUND includes detections close to the limiting magnitude and high background
    limiting magnitude points added Feb 12, 2017
 */

/*
 *  BFLAGS keyword (nonfatal flags)
 *
 *     Bit   Decimal Hexadecimal
 *                             SEXTRACTOR FLAGS
 *      0          1          1         Object has neighbors
 *      1          2          2         Object was blended with another
 *      2          4          4         At least one pixel is saturated
 *      3          8          8         Object is too close to the image boundary
 *      4         16         10         Object aperture data incomplete or corrupted
 *      5         32         20         Object isophotal data incomplete or corrupted
 *      6         64         40         Memory overflow during deblending
 *      7        128         80         Memory overflow during extraction
 *                              Miscellaneous flags
 *      8        256        100 Magnitude adjusted for a blended star (currently unused)
 *      9        512        200 Object has low drad but it is in a high drad bin (currently unused)
 *     10       1024        400 Object PSF considered saturated
 *     11       2048        800 Magnitude-dependent correction has been applied
 *                              GSC2.3 Magnitude substitutions
 *                                      NOTE: The meaning of the next for bits depend on bit 23
 *     12       4096       1000         JMag substituted for JpgMag  NOTE! Changes for these four flags affect xcolorterm.m
 *     13       8192       2000         BMag substituted for JpgMag  NOTE! These colorflags are deprecated.  See below
 *     14      16384       4000         RMag substituted for FpgMag
 *     15      32768       8000         VMag substituted for FpgMag
 *     16      65536      10000         Good star processed by the pipeline
 *     17     131072      20000         Star has been calibrated with the lowess fit
 *     18     262144      40000         Star has been locally calibrated
 *     19     524288      80000         Magnitude has been corrected for extinction
 *     20    1048576     100000         Magnitude calibration not possible because the
 *                                      star is too bright
 *     21    2097152     200000         Bin has a valid color correction
 *     22    4194304     400000         Color correction used the metropolis algorighm
 *     23    8388608     800000         Interpret bits 12-15 as Kepler Input Catalog Source bits
 *                                         (See KEPLER_CQ... below);
 *     27  134217728    8000000         Spatial bin fails QUALITY_COLORTERM below
 *     28  268435456   10000000         Ra and dec have been adjusted by the bin medians
 *     29  536870912   20000000         This plate has an uncertain Julian Date
 *     30 1073741824   40000000         Catalog position corrected for proper motion
 *
 *  AFLAGS keyword (flags for bad points)
 *      6         64         40         Object has high BACKGROUND level
 *      7        128         80         Plate fails quality criterion (see QUALITY_* below)
 *      8        256        100         Unmatched object on a multiple-exposure plate
 *      9        512        200         Object has an uncertain time for extinction calculation
 *     10       1024        400         Object is a multiple exposure blend
 *     11       2048        800         Maximum isophotonic rms exceeded
 *     12       4096       1000         Maximum locally smoothed rms exceeded
 *     13       8192       2000         Star too close to the limiting magnitude
 *     14      16384       4000         Star is in bin 9
 *     15      32768       8000 Object rejected because bin has unknown drad.
 *                             GSC2.3 CLASS FLAG
 *   16-18     65536      10000         Galaxy         (not in GSC 2.3.2)
 *   16-18    131072      20000         Blend          (not in GSC 2.3.2)
 *   16-18    196608      30000         Object is Non-star
 *   16-18    262144      40000         Unclassified   (not in GSC 2.3.2)
 *   16-18    327680      50000         Defect         (not in GSC 2.3.2)
 *   16-18    458752      70000         Duplicate Star (Added for DASCH project)
 *     19     524288      80000         Variable or Uncertain Catalog magnitude
 *                             Pipeline filterblended flags
 *     20    1048576     100000         Case B - blended
 *     21    2097152     200000         Case C - overlapping Sextractor hits for one GSC star
 *     22    4194304     400000         Case D - mixture of cases B and C
 *
 *     23    8388608     800000 Object has a drad three times the bin drad or is in a bad spatial bin or local bin
 *     24   16777216    1000000 Object is a Pickering Wedge object
 *     25   33554432    2000000 Object fails the defect filter
 *     26   67108864    4000000 Copy of the Sextractor blended flag
 *     27  134217728    8000000 Rejected blended object - avoid plotting a limiting magnitude.
 *     28  268435456   10000000 Rejected because of high smoothing correction
 *     29  536870912   20000000 Object is too bright for accurate magnitudes
 *     30 1073741824   40000000 Object is within 23.5 degrees of horizon
 */
/* The following definitions must agree with cnvgsc2.cpp and Gsc23Dasch.cpp */
#define NEW_COLOR_CALIBRATION 1

#ifndef NEW_COLOR_CALIBRATION

#define JMAG_FOR_JPGMAG     1
#define BMAG_FOR_JPGMAG     2
#define RMAG_FOR_FPGMAG     4
#define VMAG_FOR_FPGMAG     8
#else /* NEW_COLOR_CALIBRATION */
/* Definitions added on December 11, 2009, making the above definitions obsolete.
 *  NOTE: bit 3 is not to be used unless COLOR_RATIO corrections are backed out of the
 *  Pipeline
 *
 *  This memorandum was sent on  09/23/09 10:29
 *  The comparison between results for the GSC and KIC catalogs shows that our photometry suffers for the brighter star range where the GSC uses Tycho2 magnitudes.  This memorandum proposes new color transformations to put all of the magnitude sources from the GSC2.3.2 catalog on the same standard.
 *
 *  Initial studies showed that the Tycho2 magnitudes were in three separate populations, but Doug Mink correctly pointed out that these populations are probably separate color classes.  Indeed, a histogram of Tycho2 sources shows three distinct peaks, one at  color = 0.2, the second at color = 0.33 and the third at color = 0.52.  The transformations
 *  reported below should hopefully take into account these color differences.
 *
 *  I have identified 6 separate star populations in the GSC2.3.2 catalog. The "00", "02", "08" and "10" suffixes come from the MAGFlags bitmap that the Pipeline uses.
 *
 *  "gsc00"  - Magnitudes with Valid JpgMag and FpgMag
 *  "gsc02"  - Magnitudes with Bmag and Fpgmag
 *  "gsc08"  - Magnitudes with JpgMag and Vmag
 *  "gsc10"  - Magnitudes with Bmag and Vmag
 *  "tycho2" - Tycho2 substituted magnitudes with Bmag and Vmag
 *  "skymap" - Skymap substituted magnitudes with Bmag and Vmag
 *
 *  With the KIC, we now have a seventh standard star population.  I performed a least-mean-square fit to the formulas
 *               Xmag   = D + B*Tmag + C*Tcol
 *               Xcolor = D + B*Tmag + C*Tcol
 *  where Xmag,Xcolor apply to one of the six populations listed above and Tmag,Tcol is from the KIC.
 *
 *   Xmag   = D + B*Tmag + C*Tcol
 *     "gsc00"  Records  707378, dval: 0.635098, bval: 0.97327,      cval : 0.036179  std: 0.189454
 *     "gsc02"  Records    1196, dval: 0.865432, bval: 0.967795,     cval : -0.322052 std: 0.360847
 *     "gsc08"  Records    2073, dval: 1.00462,  bval: 0.941062,     cval : 0.147403  std: 0.209364
 *     "gsc10"  Records     178, dval: 0.704539, bval: 0.958592,     cval : 0.293712  std: 0.323132
 *     "tycho2" Records   13459, dval: 1.01557,  bval: 0.925296,     cval : 0.626995  std: 0.222655
 *     "skymap" Records    2014, dval: 0.207705, bval: 0.995833,     cval : 0.454896  std: 0.056671
 *
 *   Xcol = D + B*Tmag + C*Tcol
 *     "gsc00"  Records  707378, dval: 0.806126, bval: -0.0202782,   cval : 0.989274  std: 0.172765
 *     "gsc02"  Records    1196, dval: 0.799553, bval: -0.00915018,  cval : 0.612219  std: 0.339703
 *     "gsc08"  Records    2073, dval: 1.09247,  bval: -0.0562771,   cval : 0.812733  std: 0.245912
 *     "gsc10"  Records     178, dval: 0.875439, bval: -0.0444796,   cval : 1.05125   std: 0.27832
 *     "tycho2" Records   13459, dval: 0.635069, bval: -0.0366587,   cval : 1.0445    std: 0.231810
 *     "skymap" Records    2014, dval: 0.20672,  bval: -0.000498154, cval : 1.00922   std: 0.054470
 *
 *   Note that only the "skymap" population matches the KIC with any reasonable accuracy.
 *
 *  Given the above relationships to the KIC, we can solve for the following transformations of all of the star populations to the "gsc00" reference. The first entry is a math check which shows that the gsc00 -> gsc00 transformation reduces to a unity matrix.
 *
 *  gsc00  (COLOR_JPGMAG_FPGMAG)
 *  Gmag = -0.00 + ( 1.00*Tmag) + (-0.00*Tcol)
 *  Gcol = -0.00 + ( 0.00*Tmag) + ( 1.00*Tcol)
 *
 *  gsc02  (COLOR_BMAG_FOR_JPGMAG)
 *  Gmag = -0.71 + ( 1.01*Tmag) + ( 0.59*Tcol)
 *  Gcol = -0.48 + (-0.01*Tmag) + ( 1.61*Tcol)
 *
 *  gsc08  (COLOR_VMAG_FOR_FPGMAG)
 *  Gmag = -0.24 + ( 1.03*Tmag) + (-0.14*Tcol)
 *  Gcol = -0.56 + ( 0.05*Tmag) + ( 1.21*Tcol)
 *
 *  gsc10  (COLOR_VMAG_BMAG_FOR_JPGMAG_FPGMAG)
 *  Gmag =  0.14 + ( 1.00*Tmag) + (-0.25*Tcol)
 *  Gcol = -0.03 + ( 0.02*Tmag) + ( 0.93*Tcol)
 *
 *  tycho  (COLOR_TYCHO)
 *  Gmag = -0.04 + ( 1.03*Tmag) + (-0.58*Tcol)
 *  Gcol =  0.19 + ( 0.02*Tmag) + ( 0.94*Tcol)
 *
 *  skymap (COLOR_SKYMAP)
 *  Gmag =  0.52 + ( 0.98*Tmag) + (-0.40*Tcol)
 *  Gcol =  0.61 + (-0.02*Tmag) + ( 0.99*Tcol)
 */
#define COLOR_JPGMAG_FPGMAG               0
#define COLOR_BMAG_FOR_JPGMAG             3
#define COLOR_VMAG_FOR_FPGMAG             4
#define COLOR_VMAG_BMAG_FOR_JPGMAG_FPGMAG 5
#define COLOR_TYCHO                       6
#define COLOR_SKYMAP                      7
#define COLOR_OBSOLETE_VMAG_FOR_FPGMAG    8
#endif /* NEW_COLOR_CALIBRATION */

#define MAX_MAGNITUDE_FLAG 16
#define MAGNITUDE_MASK    0xf
#define CLASS_MASK        0x7
#define SAME_LOCATION      16
#define SAME_ENTRY         32
#define DUPLICATE_PAIR     64

#define GSC_CLASS_STAR         0
#define GSC_CLASS_GALAXY       1
#define GSC_CLASS_BLEND        2
#define GSC_CLASS_NONSTAR      3
#define GSC_CLASS_UNCLASSIFIED 4
#define GSC_CLASS_DEFECT       5
#define GSC_CLASS_UNDEFINED    6 /* Not part of GSC2.3.2 classification */
#define GSC_CLASS_DUPLICATE    7 /* Not part of GSC2.3.2 classification */

/* MAGFlag has different definitions for the Kepler input catalog */
/* NOTE: do not change these fields without changing keplerSourceText in pipelineutils.c */
#define KEPLER_CQ_SCP         0
#define KEPLER_CQ_2MASS       1
#define KEPLER_CQ_PHOTO       2
#define KEPLER_CQ_TYBV        3
#define KEPLER_CQ_UNCAL       4
#define KEPLER_CQ_NOCAL       5
#define KEPLER_CQ_MAGFLAG_MAX 6
/* This flag is set to distinguish from the gsc2.3.2 format */
#define KEPLER_MAGNITUDE_FLAG 16




#define FILTER_BFLAG_NEIGHBORS       0
#define FILTER_BFLAG_BLEND           1
#define FILTER_BFLAG_SATURATED       2
#define FILTER_BFLAG_BOUNDARY        3
#define FILTER_BFLAG_APERTURE_INCOMP 4
#define FILTER_BFLAG_ISOPHOT_INCOMP  5
#define FILTER_BFLAG_DEBLEND_OVF     6
#define FILTER_BFLAG_EXTRACT_OVF     7


#define FILTER_BMASK_SEXTRACTOR ((1 << FILTER_BFLAG_NEIGHBORS) | (1 << FILTER_BFLAG_BLEND) | (1 << FILTER_BFLAG_SATURATED) | (1 << FILTER_BFLAG_BOUNDARY) | (1 << FILTER_BFLAG_APERTURE_INCOMP) | (1 << FILTER_BFLAG_ISOPHOT_INCOMP) | (1 << FILTER_BFLAG_DEBLEND_OVF) | (1 << FILTER_BFLAG_EXTRACT_OVF))

#define FILTER_BLEND_MASK_SEXTRACTOR ((1 << FILTER_BFLAG_NEIGHBORS) | (1 << FILTER_BFLAG_BLEND))



#define FILTER_BFLAG_ADJUST_BLEND    8
#define FILTER_BFLAG_BINDRAD         9
#define FILTER_BFLAG_PSFSATURATED   10 /* Image PSF is saturated.   ((ISO5/ISO4 > 0.7) && (BACKGROUND+FLUX_MAX > 7500)) */
#define FILTER_BFLAG_MAGDEP_MAGCOR  11 /* Magnitude-dependent correction has been applied */
#define GSC_MAGNITUDE_FLAG_BIT      12
/* Formerly CAL_FLAGS in recover_points.c */
#define CAL_FLAG_SHIFT 16
#define PIPELINE_FLAG_BIT            0 /* Bit 16 Good star processed by the pipeline */
#define LOWESS_CAL_FLAG_BIT          1 /* Bit 17 Star has been calibrated with the lowess fit */
#define LOCAL_CAL_FLAG_BIT           2 /* Bit 18 Star has been locally calibrated */
#define EXTINCTION_FLAG_BIT          3 /* Bit 19 Magnitude has been corrected for extinction */
#define TOO_BRIGHT_FLAG_BIT          4 /* Bit 20 Magnitude calculation not possible because object is too bright */
#define COLOR_VALID_BIT              5 /* Bit 21 Bin has a valid color calculation */
#define COLOR_METROPOLIS_BIT         6 /* Bit 22 If set, color correction used the metropolis algorithm */
#define FILTER_BFLAG_KEPLER         23 /* If set, interpret bits 12-15 as KEPLER_CQ values */
#define FILTER_BFLAG_LATEMATCH      24 /* If set, the image was matched with a catalog object at the end of the pipeline */
#define FILTER_BFLAG_PMERROR        25 /* If set, this image has a 2 sigma proper motion uncertainty exceeding 1 arcsec */
#define FILTER_BFLAG_COLORTERM      27 /* If set, FILTER_AFLAG_QUALITY has been set because the spatial bin colorterm exceeds limits */
#define FILTER_BFLAG_DRAD_ADJUST    28 /* If set, the ra and dec has been adjusted by the bin medians */
#define FILTER_BFLAG_BAD_JD         29 /* If set, this plate has an uncertain Julian Date (replaced by FILTER_AFLAG_UNCERTAIN_DATE) */
#define FILTER_BFLAG_PROPER_MOTION  30 /* If set, the catalog position has been corrected for proper motion */

#define FILTER_BMASK_GSC ((MAGNITUDE_MASK << GSC_MAGNITUDE_FLAG_BIT)


#define PIPELINE_FLAG      1 /* Good star processed by the pipeline */
#define LOWESS_CAL_FLAG    2 /* Star has been calibrated with the lowess fit */
#define LOCAL_CAL_FLAG     4 /* Star has been locally calibrated */
#define EXTINCTION_FLAG    8 /* Magnitude has been corrected for extinction */
#define TOO_BRIGHT_FLAG   16 /* Magnitude calculation not possible because object is too bright */
#define COLOR_VALID       32 /* Bin has a valid color calculation */
#define COLOR_METROPOLIS  64 /* If set, color correction used the metropolis algorithm */

#define AFLAG_BIT_COUNT 32
#define FILTER_AFLAG_BACKGROUND      6
#define FILTER_AFLAG_QUALITY         7
#define FILTER_AFLAG_MULTIPLE_NONE   8
#define FILTER_AFLAG_UNCERTAIN_DATE  9
#define FILTER_AFLAG_MULTIPLE_BLEND 10
#define FILTER_AFLAG_ISO_RMS        11
#define FILTER_AFLAG_LOCAL_RMS      12
#define FILTER_AFLAG_LIMITING_MAG   13
#define FILTER_AFLAG_BIN9           14
#define FILTER_AFLAG_DRADBIN        15
#define GSC_CLASS_BIT               16
#define GSC_VARIABLE_BIT            19
#define FILTER_AFLAG_CASEB          20
#define FILTER_AFLAG_CASEC          21
#define FILTER_AFLAG_CASED          22
#define FILTER_AFLAG_DRAD           23
#define FILTER_AFLAG_WEDGE          24
#define FILTER_AFLAG_DEFECT         25
#define FILTER_AFLAG_BLEND          26
#define FILTER_AFLAG_BLEND_NOMATCH  27
#define FILTER_AFLAG_HIZOUT         28
#define FILTER_AFLAG_TOO_BRIGHT     29
#define FILTER_AFLAG_LOW_ALTITUDE   30

#define FILTER_AFLAG_BAD        1  /* FLAGS >= FILTER_FLAG_BAD is not passed to Octave */
#define FILTER_AMASK_BLEND   ((1<<FILTER_AFLAG_CASEB) | (1<<FILTER_AFLAG_CASEC) | (1<<FILTER_AFLAG_CASED) | (1<<FILTER_AFLAG_BLEND) | (1<<FILTER_AFLAG_BLEND_NOMATCH) | (1<<FILTER_AFLAG_MULTIPLE_BLEND))

#define FILTER_BMASK_BLEND      0
#define FILTER_AMASK_PLOT      ((1<<FILTER_AFLAG_CASEB) | (1<<FILTER_AFLAG_DEFECT) | (1<<FILTER_AFLAG_BLEND) | (1<<FILTER_AFLAG_LOW_ALTITUDE) | (1<<FILTER_AFLAG_WEDGE) | (1<<FILTER_AFLAG_UNCERTAIN_DATE) | (1<<FILTER_AFLAG_QUALITY) | (1<<FILTER_AFLAG_BACKGROUND) | (1<<GSC_VARIABLE_BIT) | (1<<FILTER_AFLAG_MULTIPLE_BLEND))
/* FILTER_AMASK_LIMITING is used in the resort_magfiles plate summary to find the deepest limiting magnitude on a plate */
#define FILTER_AMASK_LIMITING ((1 << FILTER_AFLAG_DRADBIN) | (1 << FILTER_AFLAG_BIN9) | (1 << FILTER_AFLAG_ISO_RMS) | (1 << FILTER_AFLAG_LOCAL_RMS) | (1 << FILTER_AFLAG_HIZOUT))

/* On Dec 19, 2012, these are now the default points NOT plotted by the web browser.
 * Modified on Feb 11, 2013 to include the limiting magnitude flag
 * Modified on Jun 18, 2013 to include the background flag
 * Modified on Oct 25, 2014 to include wedge images (together with FILTER_AMASK_LOWDRAD and FILTER_AMASK_PLOT)
 */
#define FILTER_AMASK_PLOT2     ((~ ((1<<FILTER_AFLAG_ISO_RMS) | (1<<FILTER_AFLAG_LOCAL_RMS) | (1<<FILTER_AFLAG_BACKGROUND) | (1<<FILTER_AFLAG_WEDGE) | (1<<FILTER_AFLAG_LIMITING_MAG) | (1<<FILTER_AFLAG_DRADBIN ))) & 0x7FFFFFFF)

#define FILTER_AMASK_LOWDRAD   ((1<<FILTER_AFLAG_MULTIPLE_BLEND ) | (1<<FILTER_AFLAG_DRADBIN ) | (1<<FILTER_AFLAG_CASEB ) | (1<<FILTER_AFLAG_CASEC ) | (1<<FILTER_AFLAG_CASED ) | (1<<FILTER_AFLAG_WEDGE) | (1<<FILTER_AFLAG_DRAD))
#ifndef PERMISSIVE_ID_TABLE
#define FILTER_AMASK_NONE_CANDIDATE  (FILTER_AMASK_LOWDRAD | (1<<FILTER_AFLAG_DEFECT) | (1<<FILTER_AFLAG_WEDGE))
#endif /* PERMISSIVE_ID_TABLE */
/* New definition for Transient Candidate filtering adopted Nov 6, 2014 */
/* FILTER_AFLAG_LIMITING_MAG Backed out for flare searches Feb 13, 2017 but restored on March 15, 2017 */
#define FILTER_AMASK_NONE_CANDIDATE2  (FILTER_AMASK_LOWDRAD | (1<<FILTER_AFLAG_LIMITING_MAG))

/* The following mask eliminates stars with FILTER_AFLAG_LIMITING_MAG set from the calculation of the RMS value in recover_points.c
 * At a photometry meeting of August 1, 2013, it was decided to reject blends, variables, and non-stars
 */
#define FILTER_MASK_RMS_REJECT  (FILTER_AMASK_BLEND | (CLASS_MASK << GSC_CLASS_BIT) | (1 << GSC_VARIABLE_BIT))

/* The DMAG_REJECT flags are in the rejectFlag column of the <mosaic name>_dmagcor.grid file */
#define DMAG_REJECT_HIZOUT  1  /* local correction bin rejected because ZOUT is out of range */
#define DMAG_REJECT_MEDIAN  2  /* local correction bin rejected because median brightness dimmer than the limiting magnitude */
#define DMAG_REJECT_DRAD    4  /* local correction bin rejected because majority of stars have high drad */


#define PSFSATURATED_ISO   0.7  /* ISO5/ISO4 > this value to be saturated */
#define PSFSATURATED_FLUX 7500  /* BACKGROUND+FLUX_MAX > this value to be saturated */


/* transform bitmask definitions from the mosaics MySQL table */
#define TRANSFORM_FLIP 1
#define TRANSFORM_MIRROR 2

#define WCSSOURCE_NONE 0
#define WCSSOURCE_LOGBOOK 1
#define WCSSOURCE_IMWCS   2

#define CENTERSOURCE_NONE    0
#define CENTERSOURCE_CATALOG 1
#define CENTERSOURCE_LOGBOOK 2
#define CENTERSOURCE_IMWCS   3

#define LOGBOOK_VERSION "2020-05-25" /* used by dasch_scat consistency checks on centerSource */




#define MAX_SPATIAL_BINS 9
extern double binFractions[MAX_SPATIAL_BINS+1];
extern char *keplerSourceText[KEPLER_CQ_MAGFLAG_MAX];

/* The following definitions define the binary catalog format */
#define GSC_EQUINOX 2000

#define MAX_REF 42

/*                                               dec    hex  */
#define REJECT_REASON_ERROR            0 /*        1       1  Error processing this Plate */
#define REJECT_REASON_OFFPLATE         1 /*        2       2  Object not on plate */
#define REJECT_REASON_BINFAILURE       2 /*        4       4  No data for this spatial bin */
#define REJECT_REASON_UNDETECTED       3 /*        8       8  Object undetected on this plate */
/* End of fatal reasons */
#define REJECT_REASON_WEDGE            4 /*       16      10  Object is a pickering wedge */
#define REJECT_REASON_DRAD             5 /*       32      20  Object has high DRAD */
#define REJECT_REASON_DEFECT           6 /*       64      40  Object failed the defect filter */
#define REJECT_REASON_HIZOUT           7 /*      128      80  Object in bad spatial bin */
#define REJECT_REASON_BLEND_NOMATCH    8 /*      256     100  Object in glare of brighter object */
#define REJECT_REASON_BLEND            9 /*      512     200  Object is blended */
#define REJECT_REASON_BIN9            10 /*     1024     400  Object is in bin 9 */
#define REJECT_REASON_LOCAL_RMS       11 /*     2048     800  Object has high locally calibrated RMS */
#define REJECT_REASON_LIMITING_MAG    12 /*     4096    1000  Object within 0.5 mag of the limiting magnitude */
#define REJECT_REASON_TOO_BRIGHT      13 /*     8192    2000  Object is too bright */
#define REJECT_REASON_LOW_ALTITUDE    14 /*    16384    4000  Object is too close to the horizon */
#define REJECT_REASON_ISO_RMS         15 /*    32768    8000  Object has high lowess RMS */
#define REJECT_REASON_DRADBIN         16 /*   131072   20000  Object is in a bin with too high RMS drad */
#define REJECT_REASON_COMPLEX_BLEND   17 /*   262144   40000  Object has a case c or d blend */
#define REJECT_REASON_MULTIPLE_BLEND  18 /*   524288   80000  Object is a multiple blend */
#define REJECT_REASON_UNCERTAIN_DATE  19 /*  1048576  100000  Object has an uncertain time for extinction calculation */
#define REJECT_REASON_MULTIPLE_NONE   20 /*  2097152  200000  Unmatched object on a multiple exposure plate */
#define REJECT_REASON_QUALITY         21 /*  4194304  400000  Plate fails quality criterion - see QUALITY_* below */
#define REJECT_REASON_SATURATED       22 /*  8388608  800000  Object is saturated */
#define REJECT_REASON_NOMAGDEP        23 /* 16777216 1000000  Object has no magnitude-dependent calibration */
#define REJECT_REASON_BACKGROUND      24 /* 33554432 2000000  Object has a high background level */
#define REJECT_REASON_MAX             25 /*        Total reject reasons */

#define FATAL_REASON_OFFPLATE               0
#define FATAL_REASON_NOPHOTPLATEENTRY       1
#define FATAL_REASON_NOALLOBJECTS           2
#define FATAL_REASON_ZEROPHOTPLATEVERSION   3
#define FATAL_REASON_NOMOSAICENTRY          4
#define FATAL_REASON_NOJULIANDATE           5
#define FATAL_REASON_NOPLATESENTRY          6
#define FATAL_REASON_NOSPATIALBINTABLE      7
#define FATAL_REASON_NOSPATIALBINENTRY      8
#define FATAL_REASON_ILLEGALSPATIALBIN      9
#define FATAL_REASON_NOLOCALBINTABLE       10
#define FATAL_REASON_NOLOCALBINENTRY       11
#define FATAL_REASON_LOCALBINFAILURE       12
#define FATAL_REASON_BADERRORFLAG          13
#define FATAL_REASON_BADVALIDPOINT         14
#define FATAL_REASON_NONE                  15
#define FATAL_REASON_MAX                   16

/* These flags are used in the flare transient candidate search in photometryutils.c.  They
   are defined here so that showflags can use them */
#define PEAKEVALUATION_ENTRY                     0
#define PEAKEVALUATION_INSUFFICIENT_POINTS       1
#define PEAKEVALUATION_NOPEAKS                   2
#define PEAKEVALUATION_DATERANGE                 3 /* deprecated */
#define PEAKEVALUATION_MULTIPLE1                 4
#define PEAKEVALUATION_SOFTWARE                  5
#define PEAKEVALUATION_LOW_AMPLITUDE             6
#define PEAKEVALUATION_NEARBY_DIMMAG             7
#define PEAKEVALUATION_NEARBY_AVEMAG             8
#define PEAKEVALUATION_NEARBY_DIMMAG2            9
#define PEAKEVALUATION_NEARBY_AVEMAG2           10
#define PEAKEVALUATION_NEARTYCHO2               11
#define PEAKEVALUATION_FINISHED                 12
#define PEAKEVALUATION_ACCEPTED                 13
#define PEAKEVALUATION_TOOSHORTDATERANGE        14
#define PEAKEVALUATION_TOOLONGDATERANGE         15
#define PEAKEVALUATION_INSUFFICIENT_POINTS3     16
#define PEAKEVALUATION_EXCLUSIONZONE            17
#define PEAKEVALUATION_INSUFFICIENT_POINTS2     18
#define PEAKEVALUATION_TOOMANYPOINTS            19
#define PEAKEVALUATION_LONG_INSUFFICIENT_POINTS 20
#define PEAKEVALUATION_LONG_ACCEPTED            21
#define PEAKEVALUATION_LONG_TOO_MANY_EVENTS     22
#define PEAKEVALUATION_LONG_REJECTED            23
#define PEAKEVALUATION_PRIMARY_GSC_BIN          24
#define PEAKEVALUATION_SECONDARY_GSC_BIN        25
#define PEAKEVALUATION_MAX                      26


#define SYMBOL_BIN_FAILURE    1.0  /* Dot (invisible on the graphs) */
#define SYMBOL_DEFECT         3.0  /* Asterisk */
#define SYMBOL_WEDGE          6.0  /* Open square */
#define SYMBOL_DRAD           7.0  /* Triangle */
#define SYMBOL_HIZOUT        10.5  /* Cusped Square */
#define SYMBOL_BLEND         11.5  /* Diamond */
#define SYMBOL_LOW_ALTITUDE  16.0  /* Solid Square */
#define SYMBOL_LIMITING_MAG  17.0  /* Solid Circle (17.5 for bin 9) */
#define SYMBOL_GOOD          17.0  /* Solid Circle  (17.5 for bin 9) */
#define SYMBOL_LOCAL_RMS     22.0  /* Open Circle  (22.5 for bin 9) */
#define SYMBOL_ISO_RMS       22.0  /* Open Circle  (22.5 for bin 9) */
#define SYMBOL_TOO_BRIGHT    30.0  /* Arrow pointing up */
#define SYMBOL_BLEND_NOMATCH 30.0  /* Arrow Pointing up */
#define SYMBOL_NOT_FOUND     31.0  /* Arrow pointing down (31.5 for bin 9) */
#define SYMBOL_BIN9           0.5  /* Adder for 17, 22, and 31) */
#define SYMBOL_NEIGHBORS     98.5  /* Letter "b" */
#define SYMBOL_NOCOLOR       99.0  /* Letter "c" */
#define SYMBOL_QUALITY      113.5  /* Letter "q" */
#define SYMBOL_SATURATED    115.8  /* Letter "s" */
#define SYMBOL_NOMAGDEP     109.8  /* Letter "m" */
#define SYMBOL_BACKGROUND   108.8  /* Letter "l" */
/* The following symbols are used for Topcat plotting of limiting algorithms.  These may be redefined */
#define SYMBOL_LIMITING1     40.0  /* magnitude of MAG_ISO_med+MAG_ISO_rms */
#define SYMBOL_LIMITING2     41.0  /* magnitude of dimmest star in plate sequence */

#define SCAN_PATTERN_NONE   0
#define SCAN_PATTERN_LEFT   1
#define SCAN_PATTERN_RIGHT  2
#define SCAN_PATTERN_CENTER 3


/*
 * Bins for initial magnitude filtering.  The goal is to find the limiting MAG_ISO and Stdmag as
 * early as possible to reject matches with objects that are too dim.
 *
 *  From a sample of 4176 plates on June 14, 2009, the Stdmag range is 0.03 to 19.0 and the
 *  MAG_ISO range is -22.1 to +0.91.  Use 0.1 magnitude bins from 0 < Stdmag < 19.0 and -22 < MAG_ISO < 0
 *
 *  The index into magnitudeTable is [ ((MAG_ISO - MIN_MAG_ISO)*MAG_ISO_BINS/(MAX_MAG_ISO-MIN_MAG_ISO))
 *                                     + (MAG_ISO_BINS * ((Stdmag-MIN_STDMAG)*STDMAG_BINS/(MAX_STDMAG-MIN_STDMAG)))
 *                                     + ((spatial_bin - 1) * MAGNITUDE_BINS) ]
 *
 *  The index into MAG_ISO_Table is  [((MAG_ISO - MIN_MAG_ISO)*MAG_ISO_BINS/(MAX_MAG_ISO-MIN_MAG_ISO))
 *                                     + ((spatial_bin - 1) * MAG_ISO_BINS) ]
 *
 */
#define MIN_STDMAG   0.0
#define MAX_STDMAG  19.0
#define MIN_MAG_ISO -22.0
#define MAX_MAG_ISO  0.0
#define STDMAG_BINS 190
#define MAG_ISO_BINS 220
#define MAGNITUDE_BINS (STDMAG_BINS*MAG_ISO_BINS)
/* The MAG_ISO_REJECT_FACTOR is the number of match candidates we are willing to accept for each MAG_ISO bin for each Sectractor object */

#define MAG_ISO_REJECT_FACTOR 3000

#define EXT_TIME_TOLERANCE 0.046 /* 1.1 hour time tolerance for extinction calculations */
#define DEFAULT_TIME_ACCURACY 0.0075 /* Leonid Bernikov reported (private communication) an overall accuracy of 0.003 days as a result of his O-C research.  Alison Doane notes an uncertainty in the card catalog entrys as whether geocentric or heliocentric Julian Days were in use for some years.  For an ac plate, the geo/helio uncertainty is 2 * (9 min) * sin (21 degrees) = 6.45 min = 0.0045 days.  Add the two until these issues are resolved (Must agree with Catalog.h) */
/* Note: According to Jason Eastman et al, "Achieving Better Than 1 Minute Accuracy in the Heliocentric and Barycentric Julian Dates", PASp 122:935-946, 2010 August., the maximum HJD accuracy is 8 seconds */

/* Magnitude-dependent star magnitude correction support */
/* #define MAGDEP_V2 1 */ /* Implement changes of memorandum of Fri 11/18/11 5:27 PM by Sumin Tang */
#define MAGDEP_V3 1  /* Implement changes of memorandum of Fri 11/28/11 1:47 PM by Sumin Tang */

#ifdef MAGDEP_V2
#define MAGCAL_MAGDEP_QUALITY_LIMIT 2.5
#else /* MAGDEP_V2 */
#ifdef MAGDEP_V3
#define MAGCAL_MAGDEP_QUALITY_LIMIT 2.0
#else /* MAGDEP_V3 */
#define MAGCAL_MAGDEP_QUALITY_LIMIT 1.0
#endif /* MAGDEP_V3 */
#endif /* MAGDEP_V2 */


typedef struct ranges {
  char series[MAX_SERIES_STRING]; /* Series */
  int firstPlate;
  int lastPlate;
  int locationId;
} RANGES,*PRANGES;

typedef struct locations {
  int locationId;
  double wLongitude;
  double latitude;
  double elevation;
  char name[MAX_LOCATION_NAME+2];
} LOCATIONS, *PLOCATIONS;


typedef struct exposure_struct {
  char series[MAX_SERIES_STRING]; /* Series */
  int plateNumber;
  int exposureNumber;
  char exposure[MAX_EXPOSURE_STRING];
  double dRightAscension;
  double dDeclination;
  char rightAscension[MAX_RIGHTASCENSION_STRING];
  char declination[MAX_DECLINATION_STRING];
  char date[MAX_DATE_STRING];
  char centerSource[MAX_CENTERSOURCE_STRING];
  char timeSource[MAX_CENTERSOURCE_STRING];
  char notes[MAX_NOTES_STRING];
  double timeAccuracy;
} EXPOSURE, *PEXPOSURE;

typedef struct mosaic_struct {
#if 0 /* original version */
  char series[MAX_SERIES_STRING]; /* Series */
  int plateNumber;
  int mosaicNumber;
  int scanNumber;
  int exposureNumber;
  int solutionNumber;
  int binning;
  char scanDate[MAX_DATE_STRING];
  char mosaicDate[MAX_DATE_STRING];
  int nx; /* Number of smoothing bins */
  int ny; /* Number of smoothing bins */
  int naxis1;
  int naxis2;
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double cd1_1;
  double cd1_2;
  double cd2_1;
  double cd2_2;
  int FitWCS;
  int rotation;
  int WCSSource;
  int diskLocation;
  int transform;
  int patternID;
  char ctype1[MAX_CTYPE_STRING];
  char ctype2[MAX_CTYPE_STRING];
  char mosaicComment[MAX_COMMENT_STRING];
#else /* new defintions */
  char series[MAX_SERIES_STRING]; /* Series */
  int seriesId;
  int plateNumber;
  int scanNumber;
  int mosaicNumber;
  int maxMosaicNumber; /* the highest allocated mosaic number */
  int outputFlag; /* If set, the mosaic has already been written */
  int exposureNumber;
  int solutionNumber;
  double plateScale;
  double dRightAscension;
  double dDeclination;
  int centerSource;
  char date[MAX_DATE_STRING];
  char scanDate[MAX_DATE_STRING];
  char mosaicDate[MAX_DATE_STRING];
  char plateClass[MAX_CLASS_STRING];
  char exposure[MAX_EXPOSURE_STRING];
  int WCSSource;
  int FitWCS;
  int nx; /* Number of smoothing bins */
  int ny; /* Number of smoothing bins */
  int naxis1;
  int naxis2;
  int rotation;
  int deleted;
  int diskLocation;
  int JobId; /* mosaics.JobId */
  int transform;
  int patternID;
  int binning;
  int orientation;
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double cd1_1;
  double cd1_2;
  double cd2_1;
  double cd2_2;
  double julianDate; /* Geocentric or heliocentric julian day for sorting */
  double epoch; /* Epoch of the exposure */
  char rightAscension[MAX_RIGHTASCENSION_STRING];
  char declination[MAX_DECLINATION_STRING];
  char ctype1[MAX_CTYPE_STRING];
  char ctype2[MAX_CTYPE_STRING];
  char mosaicComment[MAX_COMMENT_STRING];
#endif
} MOSAIC, *PMOSAIC;

/* GSC bin structures */


#define BIN_SIZE_01 (1.0) /* bin size in degrees (1/64 to avoid floating round off errors) */
#define DEC_BINS_01  180 /* Number of declination bins (2*90/BIN_SIZE) */
#define TOTAL_GSC_BINS_01  41164 /* We should get this number after computing the index */

#define BIN_SIZE_02 (0.5) /* bin size in degrees (1/64 to avoid floating round off errors) */
#define DEC_BINS_02  360 /* Number of declination bins (2*90/BIN_SIZE) */
#define TOTAL_GSC_BINS_02  164828 /* We should get this number after computing the index */

#define BIN_SIZE_04 (0.25) /* bin size in degrees (1/64 to avoid floating round off errors) */
#define DEC_BINS_04  720 /* Number of declination bins (2*90/BIN_SIZE) */
#define TOTAL_GSC_BINS_04  659676 /* We should get this number after computing the index */

#define BIN_SIZE_08 (0.125) /* bin size in degrees (1/64 to avoid floating round off errors) */
#define DEC_BINS_08 1440 /* Number of declination bins (2*90/BIN_SIZE) */
#define TOTAL_GSC_BINS_08  2639480 /* We should get this number after computing the index */

#define BIN_SIZE_16 (0.0625) /* bin size in degrees (1/64 to avoid floating round off errors) */
#define DEC_BINS_16  2880 /* Number of declination bins (2*90/BIN_SIZE) */
#define TOTAL_GSC_BINS_16  10559296 /* We should get this number after computing the index */

#define BIN_SIZE_32 (0.03125) /* bin size in degrees (1/64 to avoid floating round off errors) */
#define DEC_BINS_32  5760 /* Number of declination bins (2*90/BIN_SIZE) */
#define TOTAL_GSC_BINS_32  42240184 /* We should get this number after computing the index */

#define BIN_SIZE_64 (0.015625) /* bin size in degrees (1/64 to avoid floating round off errors) */
#define DEC_BINS_64  11520 /* Number of declination bins (2*90/BIN_SIZE) */
#define TOTAL_GSC_BINS_64  168966386 /* We should get this number after computing the index */

#define BIN_SIZE_128 (0.0078125) /* bin size in degrees (1/128 to avoid floating round off errors) */
#define DEC_BINS_128  23040 /* Number of declination bins (2*90/BIN_SIZE) */
#define TOTAL_GSC_BINS_128 675877030 /* We should get this number after computing the index */


#define MAX_ADJACENT_BINS 10 /* Range is 7-9 average is 8.000108 bins */
/* Expansion of bins in degrees.
 * Examinination of APASS_J075731.1+201735 (the first DASCH symbiotic)
 * suggests a scatter range of +/-18" or 3 ac telescope pixels.  The above
 * value is double this amount, or 32 arcsec
 * at a meeting on Feb 14, 2013, Prof. Grindlay requested that this radius be expanded to the width of a gsc bin
 */
#define BIN_EXPANSION_RADIUS (1.0/64.0)

/* FitWCS bitmask definitions */
#define FITWCS_SELECTED                   1  /* Selected for fitting */
#define FITWCS_MIKESHAWFIT                2  /* Previously fitted by Mike Shaw (deprecated) */
#define FITWCS_COMPLETED                  4  /* Fitting process successful */
#define FITWCS_COPIED                     8  /* Fit results copied to head cluster (deprecated) */
#define FITWCS_MULTFAILEDASTROMETRY      16  /* Fit results failed AstrometryWCS for next exposure */
#define FITWCS_MULTFAILEDSEPARATION      32  /* Fit results failed the separation criteria for the next exposure */
#define FITWCS_MULTSUCCEEDED             64  /* Fit results succeeded for the next exposure */
#define FITWCS_STALESOLUTION            128  /* Better astrometry invalidates this solution */
#define FITWCS_STALEEXPOSURENUMBER      256  /* The exposure number of this solution has been changed */
#define FITWCS_WEDGE                    512  /* This is a Pickering Wedge plate */
#define FITWCS_FILTERBLENDEDTIMEOUT    1024  /* The filterblended routine timed out */
#define FITWCS_NOALLOBJECTS            2048  /* The allobjects file for this solution is missing */
#define FITWCS_DUPLICATEREFCOUNT       4096  /* This solution shares the same stars with a previous solution */
#define FITWCS_COLORTERMCRASH          8192  /* The colorterm.m script has crashed and needs to be rerun */
#define FITWCS_INACCURATEPV           16384  /* The SCAMP polynomial fit has significant inaccuracy */
#define FITWCS_NOMULTMASK             32768  /* No multiple exposure mask found */
#define FITWCS_HAVEMULTMASK           65536  /* Have a multiple exposure mask */
#define FITWCS_FULLFAILEDASTROMETRY  131072  /* Astrometry.net failed on the full size mosaic */
#define FITWCS_FULLFAILEDSCAMPUSNO   262144  /* SCAMP using Astrometry.net catalogs failed on the full size mosaic */
#define FITWCS_FULLFAILEDSCAMPUCAC   524288  /* SCAMP using the UCAC catalog failed on the full size mosaic */
#define FITWCS_FULLFAILEDTOLERANCE  1048576  /* SCAMP did not produce a better solution with the full size mosaic */
#define FITWCS_FULLSUCCEEDED        2097152  /* SCAMP succeeded on the full size mosaic */
#define FITWCS_ASTROMETRY2          4194304  /* find_astrometry2 produced a result (for the next solution) */
#define FITWCS_NOAPASSALLOBJECTS    8388608  /* The APASS allobjects file for this solution is missing */

#define FITWCS_MASK_FULLFAILURE (FITWCS_FULLFAILEDASTROMETRY|FITWCS_FULLFAILEDSCAMPUSNO|FITWCS_FULLFAILEDSCAMPUCAC|FITWCS_FULLFAILEDTOLERANCE|FITWCS_FULLSUCCEEDED)

#define MAX_ELLIPTICITY 0.6                /* Adopted at photometry meeting of July 11, 2013 */
/* Plate Quality mask definitions */
#define QUALITY_MULTIPLE                1  /* bit  0 0x00000001 Multiple exposure */
#define QUALITY_GRATING                 2  /* bit  1 0x00000002 Grating exposure */
#define QUALITY_COLOR                   4  /* bit  2 0x00000004 Color filter used */
#define QUALITY_COLORTERM               8  /* bit  3 0x00000008 Fails colorterm limits (-0.5 to +0.25 for original GSC2.3.2 color calibration) */
                                           /*                        (-0.067 to 0.283 for recalibrated GSC2.3.2) */
                                           /*                        See MIN_BLUE_COLORTERM and MAX_BLUE_COLORTERM below */
#define QUALITY_WEDGE                  16  /* bit  4 0x00000010 Pickering Wedge plate */
#define QUALITY_SPECTRA                32  /* bit  5 0x00000020 Spectra plate */
/* The following bits are for web display only and because they are not saved in the database, can change at will */
/* NOTE: changing the postions of  QUALITY_PATROL and QUALITY_NONPATROL require changes to the mode=2 curquality setting in lightcurve_plot.php redisplay() */
#define QUALITY_SATURATED              64  /* bit  6 0x00000040 Saturated images */
#define QUALITY_NOMAGDEP              128  /* bit  7 0x00000080 No magnitude-dependent correction */
#define QUALITY_PATROL                256  /* bit  8 0x00000100 Patrol plates  (nominal scale >= PATROL_PLATE_SCALE) */
#define QUALITY_NONPATROL             512  /* bit  9 0x00000200 Non-patrol plates (nominal scale < PATROL_PLATE_SCALE) */
#define QUALITY_LIMITING             1024  /* bit 10 0x00000400 Limiting magnitude */
#define QUALITY_UNDETECTED           2048  /* bit 11 0x00000800 Not found on plate */
#define QUALITY_TRAILED              4096  /* bit 12 0x00001000 ELLIPTICITY > MAX_ELLIPTICITY = 0.6 or width = 0.4 * length */

/* New multiple exposure source quality bits */
#define QUALITY_MULTIPLE_CLASS       8192  /* bit 13 0x00002000 Class designation is multiple */
#define QUALITY_MULTIPLE_LOGBOOK    16384  /* bit 14 0x00004000 More than one logbook entry */
#define QUALITY_MULTIPLE_COMMENT    32768  /* bit 15 0x00008000 a comment by stacks personal */
#define QUALITY_MULTIPLE_FITWCS     65536  /* bit 16 0x00010000 multiple astronomy.net solutions */
#define QUALITY_MULTIPLE_MASK      131072  /* bit 17 0x00020000 is close correlation with search_close.c */
#define QUALITY_WEDGE_COMMENT      262144  /* bit 28 0x00040000 is noted as a wedge plate */


#define DATABASE_QUALITY_MASK (QUALITY_MULTIPLE|QUALITY_GRATING|QUALITY_COLOR|QUALITY_COLORTERM|QUALITY_WEDGE|QUALITY_SPECTRA)

#define ALL_QUALITY_MASK (QUALITY_MULTIPLE|QUALITY_GRATING|QUALITY_COLOR|QUALITY_COLORTERM|QUALITY_WEDGE|QUALITY_SPECTRA|QUALITY_NOMAGDEP|QUALITY_SATURATED|QUALITY_PATROL|QUALITY_NONPATROL|QUALITY_LIMITING|QUALITY_UNDETECTED|QUALITY_TRAILED)
#define DEFAULT_QUALITY_MASK (QUALITY_MULTIPLE|QUALITY_GRATING|QUALITY_COLOR|QUALITY_COLORTERM|QUALITY_WEDGE|QUALITY_NOMAGDEP|QUALITY_SATURATED|QUALITY_PATROL|QUALITY_NONPATROL|QUALITY_TRAILED)
#define QUALITY_UNINITIALIZED  (0x40000000 & (~ALL_QUALITY_MASK)) /* This field has not yet been initialized */
#if 0
/* The following colorterm limits were requested by Sumin Tang on 12/11/2009 1:07 AM.
   They are used for calculating the statstics in the id table */
#define ORIGINAL_V354_COLORTERMS 1
#define MIN_BLUE_COLORTERM -0.5
#define MAX_BLUE_COLORTERM  0.25

#else
#if 0
/* The following colorterm limits are based on the recalibrated GSC2.3.2 catalog
 * (New magnitudes and colors for Tycho and Skymap objects) and are the 5%
 * histogram limits for the colorterm distribution as agreed upon at a
 * photometry meeting of January 14, 2010. The colorterm extract from 708 M44 plates is in ~/scanner/backup/2010_01_07/colorterm_M44_new.txt and graphed in
 * http://hea-www.harvard.edu/DASCH/Internal/private/Photometry/KeplerInputCatalog/colorterm_recalibrated_gsc.gif file
*/
#define MIN_BLUE_COLORTERM -0.067 /*  5% limit */
#define MAX_BLUE_COLORTERM  0.283 /* 95% limit */
/* The following limits are from a sample 15517 values of the kepler dataset on January 23 (../2010_01_21/colorterm_kepler.txt) */
/* Use these values also for APASS */
#define MIN_KEPLER_BLUE_COLORTERM -0.067  /*  5% limit */
#define MAX_KEPLER_BLUE_COLORTERM  0.421  /* 95% limit */
#else
/* The following colorterm limits were calculated on November 9, 2012 using only metropolis colorterm values.
 * Extracts from the spatialbin tables are in /dasch/scanner/backup/2012_11_08/colorterm*.db
 *                 Plates        bins          5% limit  95% limit
 *        gsc       34956       283522       -0.097      +0.41
 *        kepler     3695        26191       -0.102      +0.409
 *        apass     32913       255054       -0.146      +0.195
 */
#define MIN_BLUE_COLORTERM                -0.097  /*  5% limit */
#define MAX_BLUE_COLORTERM                 0.410  /* 95% limit */
#define MIN_KEPLER_BLUE_COLORTERM         -0.102  /*  5% limit */
#define MAX_KEPLER_BLUE_COLORTERM          0.409  /* 95% limit */
#define MIN_APASS_BLUE_COLORTERM          -0.4    /* apass_colorterm.png October, 2018 */
#define MAX_APASS_BLUE_COLORTERM           0.4    /* < 1000 points in histogram */
#define MIN_GAIA_BLUE_COLORTERM            0.2    /*  gaiadr2_colorterm.png, October, 2018 */
#define MAX_GAIA_BLUE_COLORTERM            1.1    /*   */
#define MIN_ATLAS_BLUE_COLORTERM            0.00  /* atlas_colorterm.png, October, 2018 */
#define MAX_ATLAS_BLUE_COLORTERM            1.0  /*   */
/* the experimental dataset uses GSC colorterms */
#define MIN_EXPERIMENTAL_BLUE_COLORTERM   -0.097  /*  5% limit */
#define MAX_EXPERIMENTAL_BLUE_COLORTERM    0.410  /* 95% limit */

#endif
#endif

/*
 * From Sumin Tang, memorandum of Thu 2/03/11 4:04 PM
 * According to Allen's book, the B-R color for a O5 star is -0.48, and B-R color for a M5 supergiant (reddest) is 3.98.GSC
 * magnitude error is about 0.2 mag, color error is about 0.3mag.
 *
 * Therefore, I suggest to use a cuttoff of -0.8 & 4.3. For any star bluer than -0.8 or redder than 4.3, they are likely
 * with wrong magnitudes, and we exclude them for calibration (i.e. the calibration curve fitting and local calibration).
 * Also, no color-term corretion for these stars since the colors are wrong. However, we do wanna derive their lightcurves.
 */

#define MIN_GSC_COLOR -0.8
#define MAX_GSC_COLOR  4.3

/* From Sumin Tang, memorandum of Thu 3/03/11 3:53 PM
 * According to Covey et al. 2007 (AJ 134, 2398), g-r color is -0.62 for O5 star, and is 1.72 for M6. Some carbon stars could
 * be redder. Therefore, I suggest to use a cutoff of -1 & 2.5.  For any star bluer than -1 or redder than 2.5, they are likely
 * with wrong magnitudes, and we exclude them for calibration (i.e. the calibration curve fitting and local calibration). Also,
 * no color-term corretion for these stars since the colors are wrong. However, we do wanna derive their lightcurves.
 */
#define MIN_KEPLER_COLOR -1.0
#define MAX_KEPLER_COLOR  2.5
#if 1
/* Switch to B and V:
 * From Sumin Tang, memorandum of Wed 10/12/11 1:21 AM
 * According to Pickles 1998, PASP, 110, 863, the B-V color for a O5 MS
 * star (bluest in table 2) is -0.38, and B-V color for a M2 supergiant
 * (reddest in table 2) is 1.73 (carbon stars could be redder though).
 * Therefore, to allow some uncertainties, I would suggest to use a
 * cutoff of -0.8 & 2.5.
 */

#define MIN_APASS_COLOR -0.8
#define MAX_APASS_COLOR  2.5
#else
/* Since we are using SDSS g and r APASS magnitudes, make them the same as that of KEPLER */
#define MIN_APASS_COLOR -1.0
#define MAX_APASS_COLOR  2.5
 #endif

  /* Gaia DR2 placeholders */
#define MIN_GAIA_COLOR -99.0
#define MAX_GAIA_COLOR  99.0

  /* ATLAS refcat2 placeholders */
#define MIN_ATLAS_COLOR -99.0
#define MAX_ATLAS_COLOR  99.0

/* Return values for GetPlateColor */
#define PLATECOLOR_BLUE    0
#define PLATECOLOR_YELLOW  1
#define PLATECOLOR_RED     2

/* Excluded series list */
#define NUM_EXCLUDE_SERIES 4
extern char *excludeSeriesList[];

/* Fields in the "flatfields" table */
#define CALIBRATION_TYPE_NULL      0
#define CALIBRATION_TYPE_FLAT      1
#define CALIBRATION_TYPE_DARK      2
#define CALIBRATION_TYPE_LINEARITY 3
#define CALIBRATION_MAX_TYPE       4

#define CALIBRATION_STATUS_NULL      0
#define CALIBRATION_STATUS_ACQUIRED  1
#define CALIBRATION_STATUS_PROCESSED 2
#define CALIBRATION_STATUS_FAILED    3
#define CALIBRATION_STATUS_DELETED   4
#define CALIBRATION_MAX_STATUS       5

#define MAX_FLATFIELDS_SUFFIX 10

typedef struct _flatfields {
  char calibrationDate[MAX_DATE_STRING];
  char alternateDate[MAX_DATE_STRING];
  int calibrationType;
  int diskLocation;
  int JobId;             /* Bacula primary backup job id */
  int calibrationStatus;
  int delay;
  char suffix[MAX_FLATFIELDS_SUFFIX];
} FLATFIELDS,*PFLATFIELDS;


/* Polygon jacket picture */
#define FLOODJACKET_YES   1
#define FLOODJACKET_NO    2


/* The bin index allows us to convert a right ascension and
   declination to a GSC232 bin. */
typedef struct _binindex {
  double declination;
  int startBin;
  int numBins;
} BININDEX,*PBININDEX;



typedef struct _gsc_bin {
  double bin_size;  /* bin size in degrees (1/2**N) to avoid floating round-off errors */
  int dec_bins;     /* Number of declination bins (2 * 90 /bin_size) */
  int total_gsc_bins; /* We should get this number after computing the index */
  int max_ra_bins;  /* The maximum number of ra bins for a declination bin */
  PBININDEX pBinMasterIndex; /* Index into the bins */
} GSCBIN,*PGSCBIN;




/* The starindex structure defines the record structure of the index file */

typedef struct _starindex {
  off_t offset;  /* Offset into the output file */
  int binNumber; /* Number of this bin */
  int numStars;  /* Number of stars in this bin */
} STARINDEX,*PSTARINDEX;



typedef struct _gscimage {
  long long REFNumber; /* Encoded catalog reference number */
  double ra;           /* Right Ascension in degrees */
  double dec;          /* Declination in degrees */
  float Stdmag;        /* Blue magnitude */
  float color;         /* color */
  float RaPM;          /* Right ascension proper motion in mas/yr */
  float DecPM;         /* Declination proper motion in mas/yr */
  float RaSigmaPM;     /* Error of the Right Ascension proper motion */
  float DecSigmaPM;    /* Error of the Declination proper motion */
  char class;          /* class */
  char VFlag;          /* variable flag */
  char MAGFlag;        /* magnitude flag */
  char flag;           /* padding - used by check_photometry and formatgsc */
} GSCIMAGE,*PGSCIMAGE;

typedef struct _gscapassimage {
  long long REFNumber; /* Encoded catalog reference number */
  double ra;           /* Right Ascension in degrees */
  double dec;          /* Declination in degrees */
  float Stdmag;        /* Blue magnitude */
  float color;         /* color (bmag - vmag) */
  float RaPM;          /* Right ascension proper motion in mas/yr */
  float DecPM;         /* Declination proper motion in mas/yr */
  float RaSigmaPM;     /* Error of the Right Ascension proper motion */
  float DecSigmaPM;    /* Error of the Declination proper motion */
  float gmag;          /* SDSS g color */
  float rmag;          /* SDSS r color */
  float imag;          /* SDSS i color */
  float vmagerr;       /* magnitude errors */
  float bmagerr;
  float gmagerr;
  float rmagerr;
  float imagerr;
  char class;          /* class */
  char VFlag;          /* variable flag */
  char MAGFlag;        /* magnitude flag */
  char flag;           /* padding - used by check_photometry and formatgsc */
} GSCAPASSIMAGE,*PGSCAPASSIMAGE;

/* The following structure is used by formatgsc */
typedef struct _gscimagec {
  char REFNumber[2*MAX_REF]; /* Encoded catalog reference number as an ASCII number */
  double ra;           /* Right Ascension in degrees */
  double dec;          /* Declination in degrees */
  float Stdmag;        /* Blue magnitude */
  float color;         /* color */
  float RaPM;          /* Right ascension proper motion in mas/yr */
  float DecPM;         /* Declination proper motion in mas/yr */
  float RaSigmaPM;     /* Error of the Right Ascension proper motion */
  float DecSigmaPM;    /* Error of the Declination proper motion */
  char class;          /* class */
  char VFlag;          /* variable flag */
  char MAGFlag;        /* magnitude flag */
  char flag;           /* padding - used by check_photometry and formatgsc */
  /* These colors appear only in the apass catalog */
  float gmag;          /* SDSS g color */
  float rmag;          /* SDSS r color */
  float imag;          /* SDSS i color */
  float vmagerr;       /* magnitude errors */
  float bmagerr;
  float gmagerr;
  float rmagerr;
  float imagerr;
} GSCIMAGEC,*PGSCIMAGEC;

typedef struct _REJECTBIT
{
  int rejectBit;
  const char *rejectDescr;
} REJECTBIT,*PREJECTBIT;


typedef struct _FITWCSBIT
{
  int FitWCSMask;
  const char *FitWCSDescr;
} FITWCSBIT,*PFITWCSBIT;

typedef struct _QUALITYBIT
{
  int qualityMask;
  const char *qualityDescr;
} QUALITYBIT,*PQUALITYBIT;

#define SEQUESTERED_YES 1
#define SEQUESTERED_NO  2



typedef struct series_header {
  char series[MAX_SERIES_STRING+1];
  char description[MAX_DESCRIPTION_STRING+1];
  double aperture;     /* Telescope aperture in meters */
  double nominalPlateScale;  /* Nominal Plate scale in arc seconds per millimeter */
  double fittedPlateScale;  /* Nominal Plate scale in arc seconds per millimeter */
  int seriesId;           /* Series ID for this catalog */
  int orientation;        /* Plate orientation */
  int sequestered;        /* 1 = sequestered, 2 = not */
} SERIES_HEADER, *PSERIES_HEADER;

typedef struct _plate_entry {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int locationId;
  char plateClass[MAX_CLASS_STRING];
  int jacketJpeg;
  int plateJpeg;
  int source;
  char quality[MAX_QUALITY_STRING+2];
  int floodJacket;
} PLATE_ENTRY,*PPLATE_ENTRY;

typedef struct _flagsentry {
  int flagbit;
  int flagmask;
  char *flagname; /* C Language bit definition */
  char *webname;  /* Web-Friendly name */
} FLAGSENTRY,*PFLAGSENTRY;

/* Magnitude-dependent star calibration structures */

typedef struct _magdepcorrection {
  int ixb;        /* x index */
  int iyb;        /* y index */
  int imagb;      /* magnitude index */
  int nstar_magdep;     /* number of stars in the bin */
  int magdep_bin_size;  /* Bin expansion factor */
  int flag;       /* Internal consistency check flag */
  int magdep_bin;  /* Bin number  (ixb + (kxb*iyb) * (kxb*kyb*imagb)) */
  double xcoord_magdep;  /* x pixel of bin center */
  double ycoord_magdep;  /* y pixel of bin center */
  double magdep_bin_edge;   /* upper edge of magnitude bin */
  double magdep_bin_median;      /* center magnitude of the magnitude bin */
  double magdep_bin_magcor;  /* Value of the correction */
  double magcal_magdep_rms;  /* Final estimate of magnitude-dependent rms correction */
  double magdep_bin_quality; /* (magdep_bin_magcor * sqrt(nstar_magdep))/magcal_magdep_rms */
} MAGDEPCORRECTION,*PMAGDEPCORRECTION;

typedef struct _magdeplimits {
  int kxb;        /* Maximum x index */
  int kyb;        /* Maximum y index */
  int kmagb;      /* Maximum magnitude index */
  int nrecs;      /* Total records = kxb * kyb * (kmagb-1); */
  double dx;      /* X bin increment */
  double dy;      /* Y bin increment */
  double *xcoord_magdep; /* X bin centers (length kxb) */
  double *ycoord_magdep; /* Y bin centers (length kyb) */
  double *magdep_bin_edge;  /* Edge magnitudes of the magnitude bins (length kmagb).  The first magnitude is always zero */
} MAGDEPLIMITS,*PMAGDEPLIMITS;

#define MAX_MASKS 50  /* More than this count, assume that the results are bogus */
#define MASK_KERNEL_TYPE_DEFAULT 0 /* Default algorithm */
#define MASK_KERNEL_TYPE_DOUBLE  1 /* Two masks in straight line reflected around the center, no center image */
#define MASK_KERNEL_TYPE_MAX 2

typedef struct multiple_mask {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber; /* Mosaic number */
  int maskKernelType;
  int maskCount;  /* Total masks */
  int maskIndex; /* Index of the mask (one-based) */
  int maskXMin;
  int maskXMax;
  int maskYMin;
  int maskYMax;
  int maskArea;
  int maxBinCount;
  int maskCenterFlag;
  int convolutionRadius;
  double maskCenterDistance; /* Distance of center from the origin */
  double maskAreaRatio;    /* fraction of bounding box filled */
} MULTMASK,*PMULTMASK;

/* High background level table written by search_close and read by update_sextractor and recover_points */
typedef struct _highbackground {
  int NUMBER; /* SExtractor number */
  int AFLAGS; /* FILTER_AFLAG_BACKGROUND */
  double background_ratio; /* Background clip ratio */
} HIGHBACKGROUND,*PHIGHBACKGROUND;


typedef struct _Seriesentry {
  char series[MAX_SERIES_STRING];
  double nominalPlateScale;
  double fittedPlateScale;
  int sequestered;
} SERIESENTRY,*PSERIESENTRY;

/* Exclusion zones come from /dasch/data/scanner/backup/2015_09_24/candidatedistribution.png */
typedef struct _exclusion {
  double ra1;
  double ra2;
  double dec1;
  double dec2;
} EXCLUSION,*PEXCLUSION;


typedef struct _WEDGEENTRY
{
  char series[MAX_SERIES_STRING];
  double separation;   /* Separation between primary and wedge image in arcsec */
  double meanplussigma_drad;  /* Median astrometric error for bin 9 */
} WEDGEENTRY,*PWEDGEENTRY;

extern WEDGEENTRY wedgeTable[];
extern int wedgeTableSize;



extern PRANGES pRangeTable; /* Dump of the range table */
extern int numRanges;        /* Number of range table entries */

extern PLOCATIONS pLocationTable; /* Dump of the location table */
extern int numLocations;          /* Number of location table entries */


/* the following is duplicated in ./FitsWrapper/ScannerCCD.h */
#define LIN_VERSION_V1       1 /* Original file format version */
#define LIN_VERSION_V2       2 /* Expanded header */
#define LIN_NO_ALGORITHM     0 /* No algorithm specified */
#define LIN_ORIG_ALGORITHM   1 /* Original separate fit algorithm */
#define LIN_COMMON_ALGORITHM 2 /* Common fit algorithm. */
#define LIN_MEDIAN_ALGORITHM 3 /* Median fit algorithm. */
#define LIN_FLAG 0x6c696e66  /* 'linf' */
typedef struct _linCommon {
  FILE *binHandle;
  short algorithm;
  int maxValue;
  time_t scanTime;
  int maxQuadrant;
  int fileOpen;             /* True if we have an open file */
  int delay;
  int linVersion;
  int *forwardBuffer;
  int *reverseBuffer;
  int offset;
} LINCOMMON,*PLINCOMMON;

typedef struct _linPoint {
  int quadrant;
  int expValue;
  int actValue;
} LINPOINT, *PLINPOINT;
typedef struct _linHdrV1 {
  int flag;
  short version;
  short algorithm;
  int quadrants;
  int size;
  int flatTime; /* WARNING - not compliant for 2017 */
} LINHDRV1, *PLIN_HDRV1;
typedef struct _linHdrV2 {
  int flag;
  short version;
  short algorithm;
  int quadrants;
  int size;
  int flatTime; /* WARNING - not compliant for 2017 */
  int delay; /* Delay in microseconds */
  int offset; /* Offset in microseconds */
  int linVersion; /* Linearity program version */
  int dummy[3]; /* Extra space */
} LINHDRV2, *PLIN_HDRV2;

#define MAX_STACKS_BUFFER 500
#define MAX_STACK_LOCATION 10
#define PICKED_STATUS_MAX 31
#define PICKED_COMMENT_MAX 19

extern char  *pickedCommentString[PICKED_COMMENT_MAX];
extern char *pickedStatusString[PICKED_STATUS_MAX];

/* The platecondition table is in the stacks database */
typedef struct platecondition {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int versionId;
  char date[MAX_DATE_STRING];
  int pickedComment; /* This is a bitmap of entries in the pickedCommentString */
  int valid; /* -1 = unknown; 0 = "no" 1 = "yes" */
  char pickednotes[MAX_STACKS_BUFFER];
} PLATECONDITION,*PPLATECONDITION;

/* The plateevents table is in the stacks database */

typedef struct plateevents {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int versionId;
  char date[MAX_DATE_STRING];
  int pickedStatus; /* This is the index into the pickedStatusString */
  char stackLocation[MAX_STACK_LOCATION];
  int valid; /* -1 = unknown; 0 = "no" 1 = "yes" */
  char eventnotes[MAX_STACKS_BUFFER];
} PLATEEVENTS,*PPLATEEVENTS;

#if 0
typedef struct mosaiclist {
  char series[MAX_SERIES_STRING]; /* Series */
  int seriesId;
  int plateNumber;
  int scanNumber;
  int mosaicNumber;
  int maxMosaicNumber; /* the highest allocated mosaic number */
  int outputFlag; /* If set, the mosaic has already been written */
  int exposureNumber;
  int solutionNumber;
  double plateScale;
  double dRightAscension;
  double dDeclination;
  int centerSource;
  char date[MAX_DATE_STRING];
  char scanDate[MAX_DATE_STRING];
  char mosaicDate[MAX_DATE_STRING];
  char plateClass[MAX_CLASS_STRING];
  char exposure[MAX_EXPOSURE_STRING];
  int WCSSource;
  int FitWCS;
  int nx; /* Number of smoothing bins */
  int ny; /* Number of smoothing bins */
  int naxis1;
  int naxis2;
  int rotation;
  int deleted;
  int diskLocation;
  int JobId;    /* mosaics.JobId */
  int transform;
  int patternID;
  int binning;
  int orientation;
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double cd1_1;
  double cd1_2;
  double cd2_1;
  double cd2_2;
  double julianDate; /* Geocentric or heliocentric julian day for sorting */
  double epoch; /* Epoch of the exposure */
  char rightAscension[MAX_RIGHTASCENSION_STRING];
  char declination[MAX_DECLINATION_STRING];
  char ctype1[MAX_CTYPE_STRING];
  char ctype2[MAX_CTYPE_STRING];
  char mosaicComment[MAX_COMMENT_STRING];
} MOSAICLIST,*PMOSAICLIST;
#else
#define MOSAICLIST MOSAIC
#define PMOSAICLIST PMOSAIC
#endif

double sqr(double first);
int CalculateBin(int width,int height,double xLoc,double yLoc,double *edgeDist);
int XDmagBins(int width,int height);
int YDmagBins(int width,int height);
int ParseFilename(char *filename,char* series,int *plateNumber,int *mosaicNumber,int *binning,int *rotation);
int ParseFilename2(char *filename,char* series,int *plateNumber);
int dcmp(const void *a,const void * b);
int CalcMedianAndRMS(int curObsCount,int minGoodStars,double *vector,double *med,double *rms,int doClip,double clipFactor,int zeroFlag);

int GetDecBin(PGSCBIN pGscBin,double dec,char *outputname);
int GetGSCBin(PGSCBIN pGscBin,double ra,double dec,int *pDecBin,int *pRaBin,char *outputname);
int GetBinCenter(PGSCBIN pGscBin,int gsc_bin_index,double *ra,double *dec,char *outputname);

PBININDEX GetSubbins(PGSCBIN pGscBin,int binNumber,int *pRaBin,int *pDecBin, char* outputname);
int LinearFit(int nPoints,double *X,double* Y,double * pA, double * pB, double * pMean, double * pStd, double * pResidual);
double MaxDradPixels(char *series,int plateNumber);

void SetQueryCount(int queryCount);
int GetQueryCount();

void rotateFrame(double *coords,double angle);
int CheckBlend(double centerX1, /* center of image 1 */
               double centerY1,
               double theta1,   /* Angle of a axis with respect to the x axis for image 1 */
               double aLength1, /* Half length of image 1 */
               double bLength1,
               double centerX2, /* center of image 2 */
               double centerY2,
               double theta2,   /* Angle of b axis with respect to the x axis for image 2 */
               double aLength2, /* Half length of image 2 */
               double bLength2);

int DecodeAFLAGS(int AFLAGS,int AFLAGSMASK,int BFLAGS,int BFLAGSMASK,int *pRejectReason1,int *pRejectReason2,double *pSymbol);
void CalcMemory(int factor,signed long long newMemory,signed long long *curMemory,signed long long *maxMemory);
void  ComputeRange(int ngood,
                   double *vector,
                   double *min_local,
                   double *min_local2,
                   double *max_local2,
                   double *max_local);
void ChargeTime(int entry,long long *pOldUsec,long long *deltaTable,long long *totalTable);
int GetREFNumber(char *REF,long long *pREFNumber,int *refType,int fatal,int verbose);
int GetREF(long long REFNumber,char *REF,int nospace,int fatal);
int GetSolutionJulianDate(void *pConnection,
                          char *series,
                          int plateNumber,
                          int mosaicNumber,
                          int solutionNumber,
                          int *pExposureNumber,
                          int *pNumExposures,
                          double *pSingleJulianDate,
                          double *pSingleTolerance,
                          double *pAllJulianDate,
                          double *pAllTolerance);

int SetMosaicFitWCS(void *pConnection,char* bitName,char *fileroot,int solutionNumber,int setFlag,int readFlag,int *pFitWCS);
int SetMosaicFitWCS2(char* bitName,char *fileroot,int solutionNumber,int setFlag,int readFlag,int *pFitWCS);
int SetPlateQuality(void *pConnection,char* bitName,char *fileroot,int setFlag,int readFlag,int *pQuality);
int SetPlateQuality2(char* bitName,char *fileroot,int setFlag,int readFlag,int *pQuality);
int GetDASCHNumber(double ra,double dec,long long *pREFNumber,int fatal,int refType);
int GetDASCHCoordinates(char *REF,double *ra,double *dec,int verbose,int refType);
int GetPlateColor(void *pConnection,char *series,int plateNumber);
double GetExtinctionCoefficient(int plateColor,double elevation);
int GetSeriesInfo(void *pConnection,char *series,int seriesId,PSERIES_HEADER pSeriesHeader);
int GetPlateInfo(void *connection,char *series,int plateNumber,PPLATE_ENTRY pPlateEntry);
int GetLongitude(void *connection,int locationId,double *longitude,double *latitude,double *elevation,char *name,int nameSize);
unsigned int GetMysqlError();
int GetCatalogNumber(char *catalogName);
int GetREFType(long long REFNumber);
int SelectBestMosaic(void *pConnection,char *series,int platNumber,int* mosaicNumber, int* rotation);
void UpdateQuality(void *connection,char *series,int plateNumber);
void CheckQuality(void *connection);
int IsMultipleExposure(char *series,char *plateClass);
int IsGratingExposure(char *series,char *plateClass);
int IsColorFilterPlate(char *series,char *plateClass);
int IsSpectraPlate(char *series,char *plateClass);
int GetMultipleMask(void *connection,char *series,int plateNumber,int mosaicNumber,int *pMaskCount,PMULTMASK *ppMaskTable);
int GetExposureInfo(void *connection,char *series,int plateNumber,int exposureNumber,PEXPOSURE pExposure,int logbookFlag);


/* This routine should be in magdeputil.h, but the compiler will do the wrong thing if the user leaves out magdeputil.h */
double GetMagdepBinMagcor(PMAGDEPLIMITS pMagdepLimits,PMAGDEPCORRECTION pMagdepTable,double X_IMAGE,double Y_IMAGE,double Stdmag,int *pmagdep_bin,double *pmagcal_magdep_rms,double *pmagdep_bin_magcor,int rejectFlag);

void FreeMagdepSubarrays(PMAGDEPLIMITS pMagdepLimits);
int GetMagdepBin(PMAGDEPLIMITS pMagdepLimits,double X_IMAGE,double Y_IMAGE,double Stdmag,int *pix,int *piy,int *pimag);
int CheckAuthorization(char *regionFlag,double ra, double dec,int *pReleaseField);
void FindAdjacentBins(PGSCBIN pGscBin,int gsc_bin_index,int *gsc_bin_list,int *pgsc_bin_count);
int GetReleaseLevel();
#define MAX_BITMAP_SIZE 100 /* Minimum size of the FormatFlagsBits result buffer */
void FormatFlagsBits(int AFLAGS, char *result,char *resultbits,int resultSize,int qualityFlag);
int ExecuteQuery(void *connection,char *queryString);
void GetRangeTable(void *connection);
char * GetSeriesString(int seriesId,int fatal);
int GetMosaicInfo(void *connection,
                  char *series,
                  int plateNumber,
                  int mosaicNumber,
                  int solutionNumber,
                  PMOSAIC pMosaic);
void InitBinIndex(PGSCBIN pGscBin);
int GetMargins(void *connection,
               char *series,
               int plateNumber,
               int mosaicNumber,
               int *pLeftMargin,
               int *pRightMargin,
               int *pBottomMargin,
               int *pTopMargin);
int GetLocationId(void *connection,
                  char *series,
                  int plateNumber);

int RA2fpix(double inputRA, double inputDec,int *pChn_n,double *pRow,double *pColn);  /* in RA2fpic.c */
char *urlencode(char *string, char *outstring, int outlength); /* in urlencode.c */

void local_strlwr(char *strPtr);
const char* GetRejectReasonText(int index);
const char* GetFatalReasonText(int index);
int BlueColorterm(int catalogNumber,double colorterm);
void InitLinearityCommon(PLINCOMMON pLinearityCommon);
int ReadLinearityFile(PLINCOMMON pLinearityCommon,char *binFileName);
int GetPlateCondition(void *connection,char *series,int plateNumber,char *editdate,PPLATECONDITION* pPlateConditionTable,FILE *logHandle);
int GetPlateEvents(void *connection,char *series,int plateNumber,char *editdate,PPLATEEVENTS* pPlateEventsTable,FILE *logHandle);
int GetFlatfieldRecords(void *connection, char *startCalibrationDate,char* endCalibrationDate,int calibrationType, int calibrationStatus,int diskLocation, int* pJobId, int* pDelay,char *suffix,PFLATFIELDS *ppFlatfields);
int BuildMosaicList(void *connection,
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
                    PMOSAICLIST *ppMosaicTable);



#endif /* _PIPELINEUTILS_H */
