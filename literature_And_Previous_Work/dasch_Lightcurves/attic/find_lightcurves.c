// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* find_lightcurves.c
 *
 *  Search the master output database for individual stars, and for each star found, output lc_<src_name>.db
 *  with the following columns:
 *   1: sequence
 *   2: heliocentric Julian Date
 *   3: spatial_bin (replaces plate_dist)
 *   4: flag   17   Solid Circle              for a good point
 *             17.5 Small solid circle        for a good point, bin 9
 *             22   Open Circle               for a large error point
 *             22.5 Small open circle         for a large error point, bin 9
 *             31   Arrow pointing down       for not found
 *             31.5 Small arrow pointing down for not found, bin 9
 *              3   Asterisk                  for a plate defect
 *              7   Triangle                  for a high drad object
 *              6   Open Square               for a Pickering Wedge object
 *             16   Solid Square              for a high extinction or uncertain date for calculating extinction
 *             30   Arrow pointing up         for a rejected blend object or a too bright object
 *             11.5 Diamond                   for a blended object (including multiple exposure blends)
 *             10.5 Cusped Square             for a high local correction object.
 *             98.5 Letter "b"                for a good object in which Sextractor blend and/or neighbors flag is set
 *             99.5 Letter "c"                for a good object without color correction
 *            113.5 Letter "h"                for a plate that fails QUALITY_* conditions
 *            115.5 Letter "s"                for saturated stars
 *            109.8 Letter "m"                for stars without magnitude-dependent correction
 *            108.8 Letter "l"                for stars with high BACKGROUND
 *   5: magcal_magdep
 *   6: magcal_local_rms
 *   7: nlocal
 *   8: Plate
 *   9: magcal_iso_rms
 *  10: magcal_local_error
 *  11: magcal_iso
 *  12: limiting_mag_local
 *  13: SFLAGS    - sectractor flags
 *  14: FWHM_IMAGE - full width half maximum in pixels 
 *  15: magcal_local
 *
 *  Also write a file lc_<src_name>_<qualifier>.txt with the following values for use in plotting:
 *   1.  clipped median magcal_magdep
 *   2.  clipped magcal_magdep rms
 *   3.  total number of points (including non-detections)
 *   4.  number of good points, not in bin 9 and magcal_local_rms < 1.0
 *   5.  begining year
 *   6.  ending year
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include -I/usr/include/mysql -L /dasch/install/lib -lm  -L/usr/lib64/mysql -lmysqlclient find_lightcurves.c pipelineutils.a -ltable -lutil -lwcs -o find_lightcurves 
 *los.csh
 * /dasch/Pipeline/find_lightcurves -m  -v  -s -f  -t /dasch/Pipeline/M44_REF.db -d /dasch/Pipeline/lightcurves/tmp
 *                           
 *  For all types of errors: -g 9abcdefghlmpqrstuvw 
 *
 * /dasch/Pipeline/find_lightcurves -i -g dbhquvtc -v   -s -j -f  -t /dasch/Pipeline/Malmquist_REF.db -d /dasch/Pipeline/lightcurves/Malmquist
 * Find specific stars: /dasch/Pipeline/find_lightcurves   -s -j -f  -g dbhquvtc   -d /dasch/Pipeline/lightcurves/tmp  -c S2300002478 N2313102243 N2211330177 N2230030699
 * Find from a table    /dasch/Pipeline/find_lightcurves  -i -g dbhquvtc   -s -j -f    -t /dasch/Pipeline/3c273_REF.db  -d /dasch/Pipeline/lightcurves/tmp -b N120013341 -b N120013339
 * Cone Search:         /dasch/Pipeline/find_lightcurves  -g dbhquvtc  -s -j -f  -r 00:00:10 -p 08:40:24.000 +19:41:00.00  -d /dasch/Pipeline/lightcurves/tmp
 * NONE object:         /dasch/Pipeline/find_lightcurves  -g dbhquvtc  -s -j -f  -r 00:01:00 -p  8:39:50.819 19:33:02.93 -d /dasch/Pipeline/lightcurves/tmp
 * Cone Search:         /dasch/Pipeline/find_lightcurves  -g dbhquvtc  -s -j -f  -r 00:01:00 -c N0303302214 -d /dasch/Pipeline/lightcurves/tmp  (THIS IS NOT SUPPORTED)
 *
 *  ON DELL:ALL:  Output stars 8 from 1612 plates seconds 197 MySQL queries 10405 spatialBinCount 1576 wcsCount 1584 localbinCount 2374
 *          GOOD: Output stars 6 from 1612 plates seconds  83 MySQL queries  8031 spatialBinCount 1576 wcsCount 1584 localbinCount 0
 * 
 * May  4, 2009 Edward J. Los - Adapted from extract_lightcurves.c 
 * May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 * Jul  6, 2009 Edward J. Los - Merge with current version of extract_lightcurves:
 *                              Add allLightcurves, but with no current functionality
 *                              Add GSC catalog data to the bottom of the plot.  Use the gsc_bin_index
 *                              to locate the catalog information
 *                              Add startJD and endJD to the text file to make the axis agree with the start and end ephemeris year.
 *                              Add 0.1% margins to ensure that all points get plotted
 * Jul  7, 2009 Edward J. Los - Add support for reading the plate array from the MySQL database.
 * Jul 19, 2009 Edward J. Los - Add support for individual GSC2.3.2 reference objects.
 * Jul 20, 2009 Edward J. Los - Check plate margins before before flagging a limiting magnitude.
 * Aug 18, 2009 Edward J. Los - Add versionId checking
 * Aug 20, 2009 Edward J. Los - Add Kepler Input Catalog Support
 *                              Make errors in the coverage table non-fatal
 * Oct 27, 2009 Edward J. Los - Print gsc bin of the cone search
 * Nov  3, 2009 Edward J. Los - Remove dependence on the mysql keepalive flag
 * Nov 16, 2009 Edward J. Los - Add FILTER_AFLAG_MULTIPLE_BLEND and FILTER_AFLAG_UNCERTAIN_DATE
 * Nov 24, 2009 Edward J. Los - Implement multiple exposure support
 * Dec 14, 2009 Edward J. Los - Add magnitude file support
 * Jan 22, 2010 Edward J. Los - Add plate quality support
 * Feb 23, 2010 Edward J. Los - Flag good stars for with sextractor blend or neighbor flags have been set
 *                              Add unmatched catalog objects to the cone searches
 * Mar  1, 2010 Edward J. Los - Add SFLAGS, the sextractor flags to the output file
 * Apr  9, 2010 Edward J. Los - Correct target_table reallocation issue with stale pTarget points.
 * Apr 16, 2010 Edward J. Los - Add FWHM_IMAGE
 *                              Add option to exclude four series with known fitting problems (-e qualifier).
 * Apr 20, 2010 Edward J. Los - Write full plate name for extractimage.csh
 *                              Write a region file
 *                              Add the -f flag to include all known information in the db file
 * Jun  3, 2010 Edward J. Los - Correct  "WritePoint already recorded plate" error for plates with errorFlag set
 * Jun 22, 2010 Edward J. Los - Add base_target_nrecs to fix "coverage_table_refcount -1 is negative" error when
 *                              new targets are appended to the target table
 * Aug 24, 2010 Edward J. Los - Add 'v' to the plotting recommondation for "good" lightcurve points.
 * Aug 30, 2010 Edward J. Los - Correct case for radius specified with a single catalog object
 * Sep 17, 2010 Edward J. Los - Add 'x' to avoid writing any files when searching for valid bins.
 * Feb 14, 2011 Edward J. Los - Remover NO_REJECT_FLAG message
 * Mar 15, 2011 Edward J. Los - Correct local_bin_index calculation
 * Apr  1, 2011 Edward J. Los - Make a gsc_bin_index error nonfatal
 * Nov  4, 2011 Edward J. Los - Use magcal_magdep instead of magcal_local
 * Nov 14, 2011 Edward J. Los - Add magnitude-dependent calibration for limiting magnitudes
 * Nov 28, 2011 Edward J. Los - Flag saturated stars with an "s"
 *                              Plot magcal_local as well as magcal_magdep
 * Dec 10, 2011 Edward J. Los - Add SYMBOL_NOMAGDEP and REJECT_REASON_NOMAGDEP
 * May 22, 2012 Edward J. Los - Add option to force the error bar factor to zero.
 * Jul  3, 2012 Edward J. Los - Add option for showing high background objects
 * Jul 18, 2012 Edward J. Los - Make sure that all calls to GetRef do not return spaces in object names
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Aug  1, 2012 Edward J. Los - Add colorterm, errorcolor, and colorflags colorflag to the output table
 * Aug  7, 2012 Edward J. Los - Add plate class to the output table.
 * Sep 10, 2012 Edward J. Los - Reorganize code to optimize magnitude file support.
 *                              Set Stdmag for unmatched stars to 99.0
 * Oct 12, 2012 Edward J. Los - Make the magnitude file readOnly
 * Mar 11, 2013 Edward J. Los - Add RaPM, DecPM, ra_2 and dec_2 to the magnitude file (Version 5);
 * Dec 23, 2014 Edward J. Los - V6 data format: add A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex
 * Aug 17, 2015 Edward J. Los - Remove extra "------" header from the second line of the database file.
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * Oct  2, 2015 Edward J. Los - Add FATAL_REJECT flags and new _fatal_reason.txt file
 * Oct  4, 2015 Edward J. Los - Correct errors for FATAL_REASON_NOSPATIALBINTABLE and FATAL_REASON_NOLOCALBINTABLE
 * Mar 24, 2017 Edward J. Los - Change GetStarEntry to GetStarEntry2 to move from the stars table to the starcatalog table.
 * May 29, 2018 Edward J. Los - Support the GAIA calibration
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 * Jan  4, 2019 Edward J. Los   Add verbose to LocateNoneImages and LoadNoneImages to display the files searched.
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include "magdeputil.h"
#include "searchgsc.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "daschunistd.h"
#define MAX_BUFFER 100
#define MARGIN_WIDTH 0.001
#define MAX_FILENAME 512
/* #define ACCEPT_LOW_ALTITUDE */



struct WorldCoor *
wcskinit (
          int	naxis1,		/* Number of pixels along x-axis */
          int	naxis2,		/* Number of pixels along y-axis */
          char	*ctype1,	/* FITS WCS projection for axis 1 */
          char	*ctype2,	/* FITS WCS projection for axis 2 */
          double crpix1, 
          double crpix2,	/* Reference pixel coordinates */
          double crval1, 
          double crval2,	/* Coordinates at reference pixel in degrees */
          double *cd,		/* Rotation matrix, used if not NULL */
          double cdelt1, 
          double cdelt2,	/* scale in degrees/pixel, ignored if cd is not NULL */
          double crota,		/* Rotation angle in degrees, ignored if cd is not NULL */
          int 	 equinox, /* Equinox of coordinates, 1950 and 2000 supported */
          double epoch);	/* Epoch of coordinates, used for FK4/FK5 conversion
                           * no effect if 0 */
double str2dec(		        /* Return Dec in degrees from string */
               const char* in);	/* Character string (dd:mm:ss.sss or dd.dddd) */
double str2ra(		        /* Return RA in degrees from string */
              const char* in);	/* Character string (hh:mm:ss.sss or dd.dddd) */
 
double wcsdist(	/* Compute angular distance between 2 sky positions */
               double ra1,	/* First longitude/right ascension in degrees */
               double dec1,	/* First latitude/declination in degrees */
               double ra2,	/* Second longitude/right ascension in degrees */
               double dec2);	/* Second latitude/declination in degrees */


double jd2hjd (double	dj,	/* Julian date (geocentric) */
               double	ra,	/* Right ascension (degrees) */
               double	dec,	/* Declination (degrees) */
               int	sys);	/* J2000, B1950, GALACTIC, ECLIPTIC */

void
wcs2pix (
         struct WorldCoor *wcs,	/* World coordinate system structure */
         double	xpos,double ypos,	/* World coordinates in degrees */
         double	*xpix,double *ypix,	/* Image coordinates in pixels */
         int	*offscl);	/* 0 if within bounds, else off scale */

void
dec2str (

         char	*string,	/* Character string (returned) */
         int	lstr,		/* Maximum number of characters in string */
         double	dec,		/* Declination in degrees */
         int	ndec);		/* Number of decimal places in arcseconds */

void ra2str (

             char	*string,	/* Character string (returned) */
             int	lstr,		/* Maximum number of characters in string */
             double	ra,		/* Right ascension in degrees */
             int	ndec);		/* Number of decimal places in seconds */


typedef struct _spatial_bin {
  int spatial_bin;    /* Needed */
  double limiting_mag;  /* Needed */
  double max_bright_mag;
  double limiting_iso;
  double max_bright_iso;
  double upper_limit;  /* Needed */
  double colorterm;
  double errorcolor;
  int colorflag;
} SPATIAL_BIN,*PSPATIAL_BIN;



/* coverage table in the auxscanner database */
typedef struct _coverage
{
  int coverage_bin_index;
  char series[MAX_SERIES_STRING]; /* Series */
  int plateNumber;
  int solutionNumber;
} COVERAGE,*PCOVERAGE;






typedef struct _plate {
  char Plate[MAX_PLATE_NAME];
  char series[MAX_SERIES_STRING];
  int seriesId;
  int plateNumber;
  int mosaicNumber;
  int solutionNumber;
  int errorFlag;
  double geoJulianDate;
  EXPOSURE exposure;
  PMOSAIC pMosaic;
  struct WorldCoor *wcs;
  PSPATIAL_BIN spatial_bin_table;
  int nx;                   /* dmagcor_table width */
  int ny;                   /* dmagcor_table height */
  PLOCALBIN local_bin_table; /* Local bin correction table */
  PMAGDEPLIMITS pMagdepLimits; /* Magnitude-dependent bin correction table */
  PMAGDEPCORRECTION pMagdepTable;
  int starDetected;     /* If non-zero, this object has been detected */
  int validPoint;     /* If non-zero, this point is valid */
  /* Information copied for each point */
  double magcal_magdep; /* Local magnitude, or local limiting magnitude */
  double magcal_local;
  double Date;  /* Heliocentric Julian Date */
  double magcal_local_rms;
  double magcal_iso;
  double magcal_iso_rms;
  double limiting_mag_local;
  double Stdmag;
  double magcal_magdep_rms;
  int spatial_bin;
  int local_bin_index;
  int NUMBER;
  int AFLAGS;
  int AFLAGSCOPY;
  int BFLAGS;
  int magdep_bin;
  /* Information calculated for output */
  double dradRMS2;    /* drad RMS for this bin */
  double FWHM_IMAGE;
  int npoints_local;  /* npout */
  int rejectFlag;     
  double extinction;
  double magcor_local;       /* zout */
  double magcal_local_error; /* errout */


  double clipmed;
  double flag;     /* wip flag symbol (see comment at the top of this file) */
  int reject_reason1;  /* Total  Reject reason mask.  See pipelineutils.h for defintions */
  int reject_reason2;  /* Actual Reject reason mask.  See pipelineutils.h for defintions */
  int fatal_reject_reason;
  PPHOTTARGET pTarget;
  int gsc_bin_index;
  int coverage_table_refcount;
  int coverage_table_select;
  int mosaicWidth;
  int mosaicHeight;
  int leftMargin;
  int rightMargin;
  int bottomMargin;
  int topMargin;
  int versionId;
  /* The following fields are used by the fullFlag option */
  double BACKGROUND;
  double Blendedmag;
  double dec;                 /* Declination in degrees */
  double ELLIPTICITY;
  double FLUX_ISO;
  double FLUX_MAX;
  double FWHM_WORLD;
  double ISOAREA_WORLD;
  double KRON_RADIUS;
  double MAG_APER;
  double MAG_AUTO;
  double MAG_ISO;
  double plate_dist;
  double ra;                  /* Right Ascension in degrees */
  double THETA_J2000;
  double X_IMAGE;
  double Y_IMAGE;
  int exposureNumber;
  int ISO0;
  int ISO1;
  int ISO2;
  int ISO3;
  int ISO4;
  int ISO5;
  int ISO6;
  int ISO7;
  long long REFNumber;        /* Translated reference number */

  double colorterm;
  double errorcolor;
  int colorflag;
  char plateClass[MAX_CLASS_STRING];
  double ra_2;
  double dec_2;
  double RaPM;
  double DecPM;
  int A2FLAGS;                                     /* Sextractor and pipeline flags */
  int B2FLAGS;                                     /* Sextractor and pipeline flags */
  float timeAccuracy;                              /* Estimated time accuracy in days */
  unsigned short maskIndex;                        /* Multiple exposure mask number (using search_close) */


} PLATE,*PPLATE;


extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
extern GSCBIN gscBin01;
PGSCBIN pGscBin01 = &gscBin01;
static int ndec = 3;		/* Number of decimal places in non-angles */
int binSearchCount = 0;
int cacheHits = 0;
int cacheMisses = 0;
static int printFlag = 0;

void CleanupPlate(PPLATE pPlate)
{
  pPlate->coverage_table_select = 0;
  if (pPlate->local_bin_table != NULL) {
    free(pPlate->local_bin_table);
    pPlate->local_bin_table = NULL;
  }
  
  if (pPlate->pMagdepLimits != NULL) {
    FreeMagdepSubarrays(pPlate->pMagdepLimits);
    free(pPlate->pMagdepLimits);
    pPlate->pMagdepLimits = NULL;    
  }

  if (pPlate->pMagdepTable != NULL) {
    free(pPlate->pMagdepTable);
    pPlate->pMagdepTable = NULL;    
  }

  if (pPlate->pMosaic != NULL) {
    free(pPlate->pMosaic);
    pPlate->pMosaic = NULL;
  }
  if (pPlate->wcs != NULL) {
    wcsfree(pPlate->wcs);
    pPlate->wcs = NULL;
  }
  if (pPlate->spatial_bin_table != NULL) {
    free(pPlate->spatial_bin_table);
    pPlate->spatial_bin_table = NULL;
  }
  if (pPlate->local_bin_table != NULL) {
    free(pPlate->local_bin_table);
    pPlate->local_bin_table = NULL;
  }

}




/* Sort routine to put the plates in Julian Day order */
int JulianDayCompare(const void *first, const void *second) 
{
  double numberFirst = ((PPLATE)first)->geoJulianDate;
  double numberSecond = ((PPLATE)second)->geoJulianDate;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}

int InitPlateSpatialBinTable(MYSQL* pPhotConnection,PPLATE pPlate,int verbose,int *spatialBinCount,int dumpMySQL,char *catalogString,int solutionNumberFix)
{
  PPHOT_SPATIAL_BIN spatial_bin_table = NULL;
  PPHOT_SPATIAL_BIN pSpatialBinEntry;
  int spatial_bin_nrecs = 0;
  int spatial_bin_index;
  int spatial_bin;
  PSPATIAL_BIN pSpatial_bin;
  if (pPlate->spatial_bin_table != NULL) {
    return(0);
  }
  spatial_bin_nrecs = ReadPhotSpatialBin(pPhotConnection,&spatial_bin_table,pPlate->series,pPlate->plateNumber,pPlate->solutionNumber,0,catalogString,solutionNumberFix);
  if (spatial_bin_nrecs == 0) {
#if 0
    printf("WARNING: Failed to read spatial bin table for %s%05d s%d\n",pPlate->series,pPlate->plateNumber,pPlate->solutionNumber);
#endif
    if (dumpMySQL) {
      printf("SELECT * FROM spatialbin%s WHERE seriesId = %d and plateNumber = %d and exposureNumber = %d;\n",catalogString,pPlate->seriesId,pPlate->plateNumber,pPlate->solutionNumber);
    } 
    if (pPlate->fatal_reject_reason < 0) {
      pPlate->fatal_reject_reason = (FATAL_REASON_NOSPATIALBINTABLE);
    }    
    pPlate->errorFlag = 1;
    
    return(1);
  }
#if 0
  if (verbose) {
    printf("Read %d records from the spatial bin table for plate %s%05d s%d\n",spatial_bin_nrecs,pPlate->series,pPlate->plateNumber,pPlate->solutionNumber);
  }
#endif
  pPlate->spatial_bin_table = (PSPATIAL_BIN)calloc(MAX_SPATIAL_BINS+1,sizeof(SPATIAL_BIN));
  if (pPlate->spatial_bin_table == NULL) {
    printf("ERROR: Failed to allocate a spatial_bin_table\n");
    exit(-1);
  }
  for (spatial_bin_index = 0; spatial_bin_index < spatial_bin_nrecs; spatial_bin_index++) {
    pSpatialBinEntry = &spatial_bin_table[spatial_bin_index];
    spatial_bin = pSpatialBinEntry->spatial_bin;
    if ((spatial_bin > 0) && (spatial_bin <= MAX_SPATIAL_BINS)) {
      pSpatial_bin = &pPlate->spatial_bin_table[spatial_bin];
      pSpatial_bin->limiting_mag = pSpatialBinEntry->limiting_mag;
      pSpatial_bin->max_bright_mag = pSpatialBinEntry->max_bright_mag;
      pSpatial_bin->limiting_iso = pSpatialBinEntry->limiting_iso;
      pSpatial_bin->max_bright_iso = pSpatialBinEntry->max_bright_iso;
      pSpatial_bin->spatial_bin = spatial_bin;
      pSpatial_bin->upper_limit = pSpatialBinEntry->upper_limit;
      pSpatial_bin->colorterm = pSpatialBinEntry->colorterm;
      pSpatial_bin->errorcolor = pSpatialBinEntry->errorcolor;
      pSpatial_bin->colorflag = pSpatialBinEntry->colorflag;
    } else {
      printf("ERROR: Illegal spatial bin %d\n",spatial_bin);
      exit(-1);
    }
  }
  if (spatial_bin_table != NULL) {
    free(spatial_bin_table);
  }
  (*spatialBinCount)++;
  return(0);

}


void InitPlateEntry(MYSQL *pConnection,MYSQL* pPhotConnection,PPLATE pPlate,int verbose,int ignoreBin,PPHOTGLOBAL pPhotGlobal,int dumpMySQL,char *catalogString, int solutionNumberFix) 
{
  int plateIndex;
  int index;
  char solutionString[MAX_PLATE_NAME];
  int local_bin_index;
  int ix;
  int iy;
  char *underscorePtr;
  double cd[4];
  PPHOT_SPATIAL_BIN spatial_bin_table = NULL;
  PPHOT_SPATIAL_BIN pSpatialBinEntry;
  PHOTPLATES basePhotPlates;
  PPHOTPLATES pPhotPlates = &basePhotPlates;
  PMOSAIC pMosaic;

  int spatial_bin_nrecs = 0;
  int spatial_bin_index;
  int spatial_bin;
  PSPATIAL_BIN pSpatial_bin;
  
  int exposureNumber;
  int numExposures;
  double singleJulianDate;
  double singleTolerance;
  double allJulianDate;
  double allTolerance;
  PLATE_ENTRY plateEntry;
  PPLATE_ENTRY pPlateEntry = &plateEntry;



#if 0
  if ((strcmp(pPlate->series,"rh") == 0) &&
      (pPlate->plateNumber == 1146)) {
    printf("At plate %s%05d\n",pPlate->series,pPlate->plateNumber);
  }

#endif
  sprintf(pPlate->Plate,"%s%05d",pPlate->series,pPlate->plateNumber);


  if (GetPhotPlate(pPhotConnection,pPlate->series,pPlate->plateNumber,pPhotPlates,catalogString,0)) {
    printf("WARNING: Failed to get photplate entry for %s%05d\n",pPlate->series,pPlate->plateNumber);
    if (pPlate->fatal_reject_reason < 0) {
      pPlate->fatal_reject_reason = ( FATAL_REASON_NOPHOTPLATEENTRY);
    }
    pPlate->errorFlag = 1;
    return;
  }
  if (pPhotPlates->versionId == 0) {
#if 0
    printf("WARNING: Failed to get the photplate entry for %s%05d\n",pPlate->series,pPlate->plateNumber);
#endif
    if (dumpMySQL) {
      printf("SELECT * from photplates where seriesId = %d and plateNumber = %d;\n",pPlate->seriesId,pPlate->plateNumber);
      printf("DELETE from coverage where series = '%s' and plateNumber = %d;\n",pPlate->series,pPlate->plateNumber);
    }
    if (pPlate->fatal_reject_reason < 0) {
      pPlate->fatal_reject_reason = ( FATAL_REASON_ZEROPHOTPLATEVERSION);
    }
    pPlate->errorFlag = 1;
    return;

  }
  pPlate->mosaicNumber = pPhotPlates->mosaicNumber;
  pPlate->mosaicWidth = pPhotPlates->mosaicWidth;
  pPlate->mosaicHeight = pPhotPlates->mosaicHeight;
  pPlate->leftMargin = pPhotPlates->leftMargin;
  pPlate->rightMargin = pPhotPlates->rightMargin;
  pPlate->topMargin = pPhotPlates->topMargin;
  pPlate->bottomMargin = pPhotPlates->bottomMargin;
  pPlate->versionId = pPhotPlates->versionId;


  pMosaic = (PMOSAIC)calloc(sizeof(MOSAIC),1);
  if (pMosaic == NULL) {
    printf("ERROR: Failed to allocate mosaic entry\n");
    exit(-1);
  }
  pPlate->pMosaic = pMosaic;


  if (GetMosaicInfo(pConnection,pPlate->series,pPlate->plateNumber,pPlate->mosaicNumber,pPlate->solutionNumber,pMosaic) != 1) {
    printf("ERROR: Failed to get mosaic entry for %s%05d_%02d s%d\n",pPlate->series,pPlate->plateNumber,pPlate->mosaicNumber,pPlate->solutionNumber);
    if (pPlate->fatal_reject_reason < 0) {
      pPlate->fatal_reject_reason = ( FATAL_REASON_NOMOSAICENTRY);
    }
    pPlate->errorFlag = 1;
    return;
  }
  pPlate->nx = pMosaic->nx;
  pPlate->ny = pMosaic->ny;
  
  if (pMosaic->rotation == 0) {
    sprintf(pPlate->Plate,"%s%05d_%02d_01ww",pPlate->series,pPlate->plateNumber,pPlate->mosaicNumber);
  } else {
    sprintf(pPlate->Plate,"%s%05d_%02d_01r%dww",pPlate->series,pPlate->plateNumber,pPlate->mosaicNumber,pMosaic->rotation);
  }
  if (pPlate->solutionNumber != 0) {
    sprintf(solutionString,"_s%d",pPlate->solutionNumber);
    strcat(pPlate->Plate,solutionString);
  }



  pPlate->seriesId = GetSeriesId(pPlate->series,1);



  if ((verbose > 0) && (((plateIndex+1) % 100) == 0)) {
    printf("Reading plate calibration files for plate %5d %s\n",plateIndex+1,pPlate->Plate);
  }

  if (GetSolutionJulianDate(pConnection,pPlate->series,pPlate->plateNumber,pPlate->mosaicNumber,pPlate->solutionNumber,&exposureNumber,&numExposures,&singleJulianDate,&singleTolerance,&allJulianDate,&allTolerance)) {
    printf("ERROR: Failed to get exposure record for %s\n",pPlate->Plate);
    if (pPlate->fatal_reject_reason < 0) {
      pPlate->fatal_reject_reason = ( FATAL_REASON_NOJULIANDATE);
    }
    pPlate->errorFlag = 1;
    return;
  }
  pPlate->geoJulianDate = singleJulianDate;



  if (spatial_bin_table != NULL) {
    free(spatial_bin_table);
  }

  if (GetPlateInfo(pConnection,pPlate->series,pPlate->plateNumber,pPlateEntry) == 1) {
    strcpy(pPlate->plateClass,pPlateEntry->plateClass);
  } else {
    printf("ERROR: Failed to get plate entry for %s%05d_%02d s%d\n",pPlate->series,pPlate->plateNumber,pPlate->mosaicNumber,pPlate->solutionNumber);
    if (pPlate->fatal_reject_reason < 0) {
      pPlate->fatal_reject_reason = ( FATAL_REASON_NOPLATESENTRY);
    }
    pPlate->errorFlag = 1;

  }




  return;
}

int InitPlateWCS(MYSQL *pConnection,PPLATE pPlate,int verbose,int *wcsCount) 
{
  PMOSAIC pMosaic;
  double cd[4];
  char solutionString[MAX_PLATE_NAME];
#if 0
  if ((strcmp(pPlate->series,"rh") == 0) &&
      (pPlate->plateNumber == 1146)) {
    printf("At plate %s%05d\n",pPlate->series,pPlate->plateNumber);
  }

#endif
  pMosaic = pPlate->pMosaic;
  if (pPlate->pMosaic == NULL) {
    printf("ERROR: InitPlateWCS called without allocating a mosaic\n");
    exit(-1);
  }
  if (pPlate->wcs != NULL) {
    return(0);
  }


  /* Initialize the wcs structure for this plate */
  if (strstr(pMosaic->ctype1,"DEC")) {
    char tmpPtr[MAX_CTYPE_STRING];
    double dtmp;
    int itmp;
    strcpy(tmpPtr,pMosaic->ctype2);
    strcpy(pMosaic->ctype2,pMosaic->ctype1);
    strcpy(pMosaic->ctype1,tmpPtr);
    dtmp = pMosaic->crval1;
    pMosaic->crval1 = pMosaic->crval2;
    pMosaic->crval2 = dtmp;
            


    cd[0] = pMosaic->cd2_1;
    cd[1] = pMosaic->cd2_2;
    cd[2] = pMosaic->cd1_1;
    cd[3] = pMosaic->cd1_2;

  } else {

    cd[0] = pMosaic->cd1_1;
    cd[1] = pMosaic->cd1_2;
    cd[2] = pMosaic->cd2_1;
    cd[3] = pMosaic->cd2_2;

  }



  pPlate->wcs = wcskinit(pMosaic->naxis1,
                         pMosaic->naxis2,
                         pMosaic->ctype1,
                         pMosaic->ctype2,
                         pMosaic->crpix1,
                         pMosaic->crpix2,
                         pMosaic->crval1,
                         pMosaic->crval2,
                         cd,
                         0,  /* cdelt1 */
                         0,  /* cdelt2 */
                         0,  /* crota */
                         2000, /* equinox */
                         0);   /* epoch */


  (*wcsCount)++;
  return(0);
}



/* Writepoint returns +1 for a good entry, 0 for no entry, and -1 for a fatal error */

int WritePoint(PPLATE plateArray,int plateListLength,PPHOTSTARIMAGE pMaster,int *debugFlag,int *staleEntryCount, int fullFlag) 
{
  PPLATE pPlate;
  int plateIndex;
  int iy;
  int ix;
  PFILESTARIMAGE pFileStarImage = pMaster->pFileStarImage;
#if 0
  if ((strcmp(pFileStarImage->REF,"K10561482") == 0) &&
      (strcmp(pMaster->series,"i") == 0) &&
      (pFileStarImage->plateNumber == 1761)) {
    printf("At plate %s for %s NUMBER %d\n",pMaster->Plate,pMaster->REF,pFileStarImage->NUMBER);
  }
#endif


  if (*debugFlag) {
    printf("Saving point for %s ",pMaster->REF);
    *debugFlag = 0;
  }

  for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
    pPlate = &plateArray[plateIndex];
    if (pPlate->coverage_table_select == 0) {
      continue;
    }
    if (strcmp(pMaster->Plate,pPlate->Plate) == 0) {
      /* We have a match */

      if (pFileStarImage->versionId < pPlate->versionId) {
        (*staleEntryCount)++;
        return(0);
      }

      if (pPlate->starDetected) {
        printf("ERROR: WritePoint already recorded plate %s for star %s NUMBER: %d %d\n",pMaster->Plate,pMaster->REF,pFileStarImage->NUMBER,pPlate->NUMBER);
        return(0);
      }
      pPlate->magcal_magdep = pFileStarImage->magcal_magdep;
      pPlate->magcal_local = pFileStarImage->magcal_local;
      pPlate->Date = pFileStarImage->Date;
      pPlate->spatial_bin = pFileStarImage->spatial_bin;
      pPlate->local_bin_index = pFileStarImage->local_bin_index;
      pPlate->magcal_local_rms = pFileStarImage->magcal_local_rms;
      pPlate->Stdmag = pMaster->Stdmag;
      pPlate->magcal_iso_rms = pFileStarImage->magcal_iso_rms;
      pPlate->magcal_iso = pFileStarImage->magcal_iso;
      pPlate->limiting_mag_local = pFileStarImage->limiting_mag_local;
      pPlate->npoints_local = pFileStarImage->npoints_local;
      pPlate->AFLAGS = pFileStarImage->AFLAGS;
      pPlate->AFLAGSCOPY = pMaster->AFLAGSCOPY;
      pPlate->BFLAGS = pFileStarImage->BFLAGS;
      pPlate->gsc_bin_index = pFileStarImage->gsc_bin_index;
      pPlate->starDetected = 1;
      pPlate->dradRMS2 = pFileStarImage->dradRMS2;
      pPlate->FWHM_IMAGE = pFileStarImage->FWHM_IMAGE;
      pPlate->rejectFlag = pFileStarImage->rejectFlag;
      pPlate->extinction = pFileStarImage->extinction;
      pPlate->magcor_local = pFileStarImage->magcor_local;
      pPlate->magcal_local_error = pFileStarImage->magcal_local_error;
      
      pPlate->NUMBER = pFileStarImage->NUMBER;
      if (fullFlag) {
        pPlate->BACKGROUND = pFileStarImage->BACKGROUND;
        pPlate->Blendedmag = pFileStarImage->Blendedmag;
        pPlate->dec = pFileStarImage->dec;
        pPlate->ELLIPTICITY = pFileStarImage->ELLIPTICITY;
        pPlate->FLUX_ISO = pFileStarImage->FLUX_ISO;
        pPlate->FLUX_MAX = pFileStarImage->FLUX_MAX;
        pPlate->FWHM_WORLD = pFileStarImage->FWHM_WORLD;
        pPlate->ISOAREA_WORLD = pFileStarImage->ISOAREA_WORLD;
        pPlate->KRON_RADIUS = pFileStarImage->KRON_RADIUS;
        pPlate->MAG_APER = pFileStarImage->MAG_APER;
        pPlate->MAG_AUTO = pFileStarImage->MAG_AUTO;
        pPlate->MAG_ISO = pFileStarImage->MAG_ISO;
        pPlate->plate_dist = pFileStarImage->plate_dist;
        pPlate->ra = pFileStarImage->ra;
        pPlate->THETA_J2000 = pFileStarImage->THETA_J2000;
        pPlate->X_IMAGE = pFileStarImage->X_IMAGE;
        pPlate->Y_IMAGE = pFileStarImage->Y_IMAGE;
        pPlate->exposureNumber = pFileStarImage->exposureNumber;
        pPlate->ISO0 = pFileStarImage->ISO0;
        pPlate->ISO1 = pFileStarImage->ISO1;
        pPlate->ISO2 = pFileStarImage->ISO2;
        pPlate->ISO3 = pFileStarImage->ISO3;
        pPlate->ISO4 = pFileStarImage->ISO4;
        pPlate->ISO5 = pFileStarImage->ISO5;
        pPlate->ISO6 = pFileStarImage->ISO6;
        pPlate->ISO7 = pFileStarImage->ISO7;
        pPlate->magdep_bin = pFileStarImage->magdep_bin;
        pPlate->magcal_magdep_rms = pFileStarImage->magcal_magdep_rms;
        pPlate->REFNumber = pFileStarImage->REFNumber;
        pPlate->colorterm = pMaster->colorterm;
        pPlate->errorcolor = pMaster->errorcolor;
        pPlate->colorflag = pMaster->colorflag;

        pPlate->ra_2 = pFileStarImage->ra_2;
        pPlate->dec_2 = pFileStarImage->dec_2;
        pPlate->RaPM = pFileStarImage->RaPM;
        pPlate->DecPM = pFileStarImage->DecPM;

        pPlate->A2FLAGS = pFileStarImage->A2FLAGS;
        pPlate->B2FLAGS = pFileStarImage->B2FLAGS;
        pPlate->timeAccuracy = pFileStarImage->timeAccuracy;
        pPlate->maskIndex = pFileStarImage->maskIndex;

      }

      return(1);
    }
  }
#if 1
  return(0);
#else
  printf("ERROR: WritePoint Failed to find plate %s in plate list\n",pMaster->Plate);
  return(-1);
#endif
}
void MergeBlendArray(char *blendStrings,int blendCount,PPLATE *blendArray,int *savedBlendCount, PPHOTTARGET pTarget, int plateListLength) {
  int blendIndex;
  PPLATE plateArray1;
  PPLATE plateArray2;
  PPLATE pPlate1;
  PPLATE pPlate2;
  double clipmed;
  int plateIndex;

  memset(pTarget,0,sizeof(PHOTTARGET));
  for (blendIndex = 0; blendIndex < blendCount; blendIndex++) {
    plateArray2 = blendArray[blendIndex];
    strcat(pTarget->REF,plateArray2->pTarget->REF);
    strcat(pTarget->REF,"_");
    strcat(pTarget->src_name,plateArray2->pTarget->src_name);
    strcat(pTarget->src_name,"_");



    if (blendIndex == 0) {
      plateArray1 = plateArray2;
      continue;
    }
    clipmed = plateArray2->clipmed;
    if (clipmed == 0.0) {
      for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
        pPlate2 = &plateArray2[plateIndex];
        if (pPlate2->starDetected != 0) {
          clipmed = pPlate2->Stdmag;
          break;
        }

      }
    }
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate1 = &plateArray1[plateIndex];
      pPlate2 = &plateArray2[plateIndex];
      pPlate1->AFLAGS = pPlate1->AFLAGSCOPY;
      pPlate2->AFLAGS = pPlate2->AFLAGSCOPY;
      pPlate1->validPoint = 0;
#if 0
      if (strcmp(pPlate1->Plate,"rb09369_00_01r180ww") == 0) {
        printf("At plate %s\n",pPlate1->Plate);
      }
#endif
      if ((pPlate1->errorFlag == 0) &&
          (pPlate1->starDetected != 0)) {
        double tmpFlux1;
        double tmpFlux2;
        double totalFlux;
        /* Star of interest has been detected.  See if the companion has also been detected,
           but be sure it is not an unmatched blended object */
        if ((pPlate2->errorFlag == 0) &&
            (pPlate2->starDetected != 0) &&
            ((pPlate2->AFLAGS & (1<<FILTER_AFLAG_BLEND_NOMATCH)) == 0) &&
            (pPlate2->magcal_magdep < 90.0) &&
            (pPlate2->magcal_iso < 90.0)) {          
          if (((pPlate1->AFLAGS & (1<<FILTER_AFLAG_BLEND)) == 0) &&
              ((pPlate2->AFLAGS & (1<<FILTER_AFLAG_BLEND)) == 0) &&
              ((pPlate1->AFLAGS & (1<<FILTER_AFLAG_CASED)) != 0) &&
              ((pPlate2->AFLAGS & (1<<FILTER_AFLAG_CASED)) != 0)) {
            /* Here both stars are part of a complex blend, yet neither one has been
               detected as blended by Sextractor.  Assume that they are not really
               blended but marked as such because of elliptical rather than circular trails */
            /* Clear the blend flag */
            pPlate1->AFLAGS &= ~(1<<FILTER_AFLAG_CASED);
            pPlate2->AFLAGS &= ~(1<<FILTER_AFLAG_CASED);
            if ((pPlate1->AFLAGS < FILTER_AFLAG_BAD) && (pPlate2->AFLAGS < FILTER_AFLAG_BAD)) {
              printf("WARNING: %s and %s for plate %s are not really blended\n",
                     pTarget->REF,plateArray2->pTarget->REF,pPlate1->Plate);
            }
          }


          /* Other star has been detected, sum up the magnitudes */
          tmpFlux1 = exp10(-pPlate1->magcal_magdep/2.5);
          tmpFlux2 = exp10(-pPlate2->magcal_magdep/2.5);
          totalFlux = tmpFlux1+tmpFlux2;
          if (totalFlux == 0) {
            pPlate1->magcal_magdep = 99.0;
            pPlate1->magcal_local_rms = 99.0;
          } else {
            pPlate1->magcal_magdep = -2.5*log10(totalFlux);
          }
          tmpFlux1 = exp10(-pPlate1->magcal_local/2.5);
          tmpFlux2 = exp10(-pPlate2->magcal_local/2.5);
          totalFlux = tmpFlux1+tmpFlux2;
          if (totalFlux == 0) {
            pPlate1->magcal_local = 99.0;
          } else {
            pPlate1->magcal_local = -2.5*log10(totalFlux);
          }
          tmpFlux1 = exp10(-pPlate1->magcal_iso/2.5);
          tmpFlux2 = exp10(-pPlate2->magcal_iso/2.5);
          totalFlux = tmpFlux1+tmpFlux2;
          if (totalFlux == 0) {
            pPlate1->magcal_iso = 99.0;
            pPlate1->magcal_iso_rms = 99.0;
          } else {
            pPlate1->magcal_iso = -2.5*log10(totalFlux);
          }
        } else {
          /* Other star has not been detected */
          if (((pPlate2->AFLAGS & (1<<FILTER_AFLAG_BLEND_NOMATCH)) == 0) ||
              (pPlate2->starDetected == 0)) {
            /* Other star is not part of primary star */
            if (clipmed != 0.0) {
              /* Use the DASCH median magnitude for the undetected star */
              tmpFlux1 = exp10(-pPlate1->magcal_magdep/2.5);
              tmpFlux2 = exp10(-clipmed/2.5);
              totalFlux = tmpFlux1+tmpFlux2;
              if (totalFlux == 0) {
                pPlate1->magcal_magdep = 99.0;
                pPlate1->magcal_local_rms = 99.0;
              } else {
                pPlate1->magcal_magdep = -2.5*log10(totalFlux);
              }
              tmpFlux1 = exp10(-pPlate1->magcal_local/2.5);
              tmpFlux2 = exp10(-clipmed/2.5);
              totalFlux = tmpFlux1+tmpFlux2;
              if (totalFlux == 0) {
                pPlate1->magcal_local = 99.0;
              } else {
                pPlate1->magcal_local = -2.5*log10(totalFlux);
              }
              tmpFlux1 = exp10(-pPlate1->magcal_iso/2.5);
              tmpFlux2 = exp10(-clipmed/2.5);
              totalFlux = tmpFlux1+tmpFlux2;
              if (totalFlux == 0) {
                pPlate1->magcal_iso = 99.0;
                pPlate1->magcal_iso_rms = 99.0;
              } else {
                pPlate1->magcal_iso = -2.5*log10(totalFlux);
              }

            }

          }

        }
      }
    }
  }
  pTarget->REF[strlen(pTarget->REF)-1] = 0;
  pTarget->src_name[strlen(pTarget->src_name)-1] = 0;
}

/* If the current object matches the blend list, then save the array */
void SaveLightcurve(PPLATE plateArray,int plateListLength,PPHOTTARGET pTarget,char *blendStrings,int blendCount,PPLATE *blendArray,int *savedBlendCount,double clipmed)
{
  int blendIndex;
  int saveFlag = 0;
  int curBlendCount = *savedBlendCount;
  if ((curBlendCount < 0) || (curBlendCount >= MAX_BLEND_COUNT)) {
    printf("ERROR: curBlendCount %d too high for MAX_BLEND_COUNT %d\n",curBlendCount,MAX_BLEND_COUNT);
  }
  if (blendCount == 0) {
    return;
  }
  for (blendIndex = 0; blendIndex < blendCount; blendIndex++) {
    if (strcmp(pTarget->REF,&blendStrings[blendIndex*MAX_BUFFER]) == 0) {
      saveFlag = 1;
      break;
    }

  }
  if (saveFlag == 0) {
    return;
  }
  blendArray[blendIndex] = (PPLATE)calloc(plateListLength,sizeof(PLATE));
  if (blendArray[blendIndex] == NULL) {
    printf("ERROR: Failed to allocate blend array of size %d\n",plateListLength*sizeof(PLATE));
    exit(-1);
  }
  plateArray->clipmed = clipmed;
  memcpy(blendArray[blendIndex],plateArray,plateListLength*sizeof(PLATE));
  blendArray[blendIndex]->pTarget = pTarget;
  curBlendCount++;
  *savedBlendCount = curBlendCount;
  printf("Saved blend object %s %s index %d\n",pTarget->src_name,pTarget->REF,blendIndex);
  return;
}

int WriteLightcurve(MYSQL *pConnection,MYSQL* pPhotConnection,PPLATE plateArray,int plateListLength,int match_count,PPHOTTARGET pTarget,char *lightcurveDirectory,double *vector1,double *vector2,int verbose,int goodFlag,int printReject,double *pClipMed, int AFLAGSMASK,int BFLAGSMASK,int ignoreBin,int allLightcurves,int dumpMySQL,int *spatialBinCount,int *wcsCount,int *localbinCount,char *catalogString,int solutionNumberFix,int fullFlag,int forceErrorBarFactor)
{
  PPLATE pPlate;
  int plateIndex;
  int curveWritten = 0;
  int sequence = 0;
  int validPoints = 0;
  int goodPoints = 0;
  int errorPoints = 0;
  int limitingPoints1 = 0;
  int limitingPoints2 = 0;
  int badPoints = 0;
  int bin9Points = 0;
  int bothPoints = 0;
  double X_IMAGE;
  double Y_IMAGE;
  int offscl;
  double edgeDist;
  int ix;
  int iy;
  int eix;
  int eiy;
  int irec;
  double rawmed = 0;
  double rawrms = 0;
  double clipmed = 0;
  double cliprms = 0;
  double rawerrrms = 0;
  double rawerrmed = 0;
  char lightcurveName[MAX_BUFFER];
  FILE *lightcurveHandle;
  double startYear = 0.0;
  double endYear;
  double startJD = 0.0;
  double endJD;
  double rangeJD;
  int wedgeCount = 0;
  int highDradCount = 0;
  int defectCount = 0;
  int hiZoutCount = 0;
  int tooBrightCount = 0;
  int lowAltitudeCount = 0;
  int blendCount = 0;
  int localBinFailure = 0;
  int reject_count_all1[REJECT_REASON_MAX];
  int reject_count_only1[REJECT_REASON_MAX];
  int reject_count_all2[REJECT_REASON_MAX];
  int reject_count_only2[REJECT_REASON_MAX];
  int fatal_reject_count_plate_s0[FATAL_REASON_MAX];
  int fatal_reject_count_plate_sn[FATAL_REASON_MAX];
  int totalPointsWritten = 0;
  int total_s0PointsWritten = 0;
  int total_snPointsWritten = 0;
  int fatal_s0PointCount = 0;
  int fatal_snPointCount = 0;
  int total_s0Count = 0;
  int total_snCount = 0;
  double fatal_s0Total;
  double fatal_snTotal;

  int rejectIndex;
  double error_bar_factor = 1.0;
  PMOSAIC pMosaic;
  char rstr[32], dstr[32];
  PLOCALBIN pLocalBin;
  int local_bin_index;
  int marginFlag;
  char MAGFlagString[MAX_BUFFER];
  double minMag = 99.0; 
  double maxMag = 0;

  *pClipMed = 0;
  if (printReject) {
    for (rejectIndex = 0; rejectIndex < REJECT_REASON_MAX; rejectIndex++) {
      reject_count_all1[rejectIndex] = 0;
      reject_count_only1[rejectIndex] = 0;
      reject_count_all2[rejectIndex] = 0;
      reject_count_only2[rejectIndex] = 0;
      
    }

    for (rejectIndex = 0; rejectIndex < FATAL_REASON_MAX; rejectIndex++) {
      fatal_reject_count_plate_s0[rejectIndex] = 0;
      fatal_reject_count_plate_sn[rejectIndex] = 0;
      totalPointsWritten = 0;
      total_s0PointsWritten = 0;
      total_snPointsWritten = 0;
      fatal_s0PointCount = 0;
      fatal_snPointCount = 0;
      total_s0Count = 0;
      total_snCount = 0;
    }
  }

  /* Go through the array and set in the correct flags value.  Also
     calculate and set in the limiting magnitude if the star was not
     detected */
  for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
    pPlate = &plateArray[plateIndex];
    if (fullFlag) {
      if (InitPlateSpatialBinTable(pPhotConnection,pPlate,verbose,spatialBinCount,dumpMySQL,catalogString,solutionNumberFix) == 0) {

        pPlate->colorterm = pPlate->spatial_bin_table[pPlate->spatial_bin].colorterm;
        pPlate->errorcolor = pPlate->spatial_bin_table[pPlate->spatial_bin].errorcolor;
        pPlate->colorflag = pPlate->spatial_bin_table[pPlate->spatial_bin].colorflag;
      }
    }


    pMosaic = NULL;
    marginFlag = 0;
    if (pPlate->coverage_table_select == 0) {
      continue;
    }
    pPlate->reject_reason1 = 0;
    pPlate->reject_reason2 = 0;


    localBinFailure = 0;

#if 0
    if ((strcmp(pTarget->REF,"N2313102243") == 0) &&
        (strcmp(pPlate->Plate,"i52702") == 0)) {
      printf("At plate %s for %s \n",pPlate->Plate,pTarget->REF);
    }
#endif


#if 0
    if (strcmp(pTarget->REF,"N2313102243") == 0) {
      printf("Writing lightcurve for %s, plate %s\n",pTarget->REF,pPlate->Plate);
    }
#endif
#if 0
    if (strcmp(pPlate->Plate,"dnb03443") == 0) {
      printf("Writing lightcurve for %s, plate %s\n",pTarget->REF,pPlate->Plate);
    }
#endif
    if (pPlate->errorFlag != 0) {
      pPlate->reject_reason1 |= (1 << REJECT_REASON_ERROR);
      pPlate->reject_reason2 |= (1 << REJECT_REASON_ERROR);
    } else {
      if (pPlate->starDetected == 0) {
        if (pPlate->wcs == NULL) {
          if (InitPlateWCS(pConnection,pPlate,verbose,wcsCount) != 0) {
            pPlate->reject_reason1 |= (1 << REJECT_REASON_ERROR);
            pPlate->reject_reason2 |= (1 << REJECT_REASON_ERROR);
           

          }
        }
        if (pPlate->errorFlag == 0) {
          wcs2pix(pPlate->wcs,pTarget->ra,pTarget->dec,&X_IMAGE,&Y_IMAGE,&offscl);
          if (!offscl) {
            /* Check plate margins */
            if ((X_IMAGE < pPlate->leftMargin) ||
                (Y_IMAGE < pPlate->bottomMargin) ||
                ((pPlate->mosaicWidth-X_IMAGE) < pPlate->rightMargin) ||
                ((pPlate->mosaicHeight-Y_IMAGE) < pPlate->topMargin)) {
#if 0
              printf("In Margin for plate %s X_IMAGE %f Y_IMAGE %f width %d height %d margins left %d right %d bottom %d top %d\n",
                     pPlate->Plate,
                     X_IMAGE,Y_IMAGE,
                     pPlate->mosaicWidth,
                     pPlate->mosaicHeight,
                     pPlate->leftMargin,
                     pPlate->rightMargin,
                     pPlate->bottomMargin,
                     pPlate->topMargin);
#endif
              marginFlag = 1;
            }
          }
          if (offscl || marginFlag) {
            pPlate->reject_reason1 |= (1 << REJECT_REASON_OFFPLATE);
            pPlate->reject_reason2 |= (1 << REJECT_REASON_OFFPLATE);
            if (pPlate->fatal_reject_reason < 0) {
              pPlate->fatal_reject_reason = ( FATAL_REASON_OFFPLATE);
            }    
          } else {
            pMosaic = pPlate->pMosaic;
            pPlate->spatial_bin = CalculateBin(pMosaic->naxis1,pMosaic->naxis2,X_IMAGE,Y_IMAGE,&edgeDist);
            ix = (X_IMAGE * pMosaic->nx) /(1.0 * pMosaic->naxis1);
            iy = (Y_IMAGE * pMosaic->ny) /(1.0 * pMosaic->naxis2);
            pPlate->local_bin_index = ix + (pMosaic->nx * iy);

            if ((pPlate->spatial_bin < 1) || (pPlate->spatial_bin > MAX_SPATIAL_BINS)) {
              pPlate->reject_reason1 |= (1 << REJECT_REASON_OFFPLATE);          
              pPlate->reject_reason2 |= (1 << REJECT_REASON_OFFPLATE);          
              if (pPlate->fatal_reject_reason < 0) {
                pPlate->fatal_reject_reason = (FATAL_REASON_ILLEGALSPATIALBIN);
              }    
            } else {
              if (pPlate->spatial_bin == 9) {
                pPlate->reject_reason1 |= (1 << REJECT_REASON_BIN9);
                pPlate->reject_reason2 |= (1 << REJECT_REASON_BIN9);
              }
              if (ignoreBin) {
                pPlate->reject_reason1 |= (1 << REJECT_REASON_UNDETECTED);
                pPlate->reject_reason2 |= (1 << REJECT_REASON_UNDETECTED);
              } else { /* ignoreBin */
                if (InitPlateSpatialBinTable(pPhotConnection,pPlate,verbose,spatialBinCount,dumpMySQL,catalogString,solutionNumberFix) != 0) {
                  /* No spatial bin table to be found */
                  pPlate->reject_reason1 |= (1 << REJECT_REASON_BINFAILURE);
                  pPlate->reject_reason2 |= (1 << REJECT_REASON_BINFAILURE);

                } else {
                  if (pPlate->spatial_bin != pPlate->spatial_bin_table[pPlate->spatial_bin].spatial_bin) {
                    if (pPlate->fatal_reject_reason < 0) {
                      pPlate->fatal_reject_reason = ( FATAL_REASON_NOSPATIALBINENTRY);
                    }    
                    
                    pPlate->reject_reason1 |= (1 << REJECT_REASON_BINFAILURE);
                    pPlate->reject_reason2 |= (1 << REJECT_REASON_BINFAILURE);
                  } else {
                    pPlate->reject_reason1 |= (1 << REJECT_REASON_UNDETECTED);
                    pPlate->reject_reason2 |= (1 << REJECT_REASON_UNDETECTED);
                  }
                }
#if 0 /* Can not do this if we never detected a star */
                /* Now see if there is a local correction for this magnitude */
                if ((pPlate->AFLAGS & (1<<FILTER_AFLAG_HIZOUT)) != 0) {
                  pPlate->reject_reason1 |= (1 << REJECT_REASON_HIZOUT);
                  pPlate->reject_reason2 |= (1 << REJECT_REASON_HIZOUT);
                }
#endif        
            
            
              } /* ignoreBin */
            }
          }
        }
      } /* Star detected */ 
    } /* Error flag clear */
    
    if (pPlate->errorFlag == 0) {
      if (pPlate->starDetected) {
       
        pPlate->validPoint = DecodeAFLAGS(pPlate->AFLAGS,AFLAGSMASK,pPlate->BFLAGS,BFLAGSMASK,&pPlate->reject_reason1,&pPlate->reject_reason2,&pPlate->flag);
        if (pPlate->flag == SYMBOL_WEDGE) {
          wedgeCount++;
        } else if (pPlate->flag == SYMBOL_DRAD) {
          highDradCount++;
        } else if (pPlate->flag == SYMBOL_HIZOUT) {
          hiZoutCount++;
        } else if (pPlate->flag == SYMBOL_TOO_BRIGHT) {
          if ((pPlate->reject_reason2 & (1<<REJECT_REASON_TOO_BRIGHT)) != 0) {
            tooBrightCount++;
          }         
        } else if (pPlate->flag == SYMBOL_LOW_ALTITUDE) {
          lowAltitudeCount++;
        } else if (pPlate->flag == SYMBOL_DEFECT) {
          defectCount++;
        }
        if ((pPlate->reject_reason2 & (1<<REJECT_REASON_BLEND)) != 0) {
          blendCount++;
        }
        if ((pPlate->reject_reason2 & (1<<REJECT_REASON_BIN9)) != 0) {
          bin9Points++;
        }
        if ((pPlate->reject_reason2 & (1<<REJECT_REASON_LOCAL_RMS)) != 0) {
          badPoints++;
        }
        if ((pPlate->flag == SYMBOL_GOOD) || (pPlate->flag == SYMBOL_NOCOLOR)  || (pPlate->flag == SYMBOL_NEIGHBORS) || (pPlate->flag == SYMBOL_SATURATED)  || (pPlate->flag == SYMBOL_NOMAGDEP)) {
          /* This currently includes stars near the limiting magnitude */
          vector1[goodPoints] = pPlate->magcal_magdep;
          vector2[goodPoints] = pPlate->magcal_local_rms;
#if 0
          printf("DEBUG: point %4d AFLAGS %08x plate %25s target %s\n",goodPoints,pPlate->AFLAGS,pPlate->Plate,pTarget->REF);
#endif
          goodPoints++;
          if (pPlate->validPoint == 0) {
            limitingPoints2++;
          }
        }
        if (pPlate->validPoint != 0) {
          validPoints++;
        }
          
      

      } else {
        if (ignoreBin) {
          continue;
        }
        /* Star has not been detected. Figure out where this
           star is on the plate */
        if (pPlate->wcs == NULL) {
          printf("ERROR: no wcs in line %d\n",__LINE__);
          exit(-1);
        }
        wcs2pix(pPlate->wcs,pTarget->ra,pTarget->dec,&X_IMAGE,&Y_IMAGE,&offscl);
        
        if (offscl || marginFlag) {
          if (pPlate->fatal_reject_reason < 0) {
            printf("ERROR line %d for errorPoints %d\n",__LINE__,errorPoints);
          }
          errorPoints++;
          continue;
        }
#if 0
        if ((strcmp(pTarget->REF,"DASCH_J192827.0+431358") == 0) && (strcmp(pPlate->Plate,"ac02838_00_01ww") == 0)) {
          printf("At plate %s for %s \n",pPlate->Plate,pTarget->REF);
        }
#endif
        pPlate->X_IMAGE = X_IMAGE;
        pPlate->Y_IMAGE = Y_IMAGE;
        if ((pMosaic == NULL) || 
            (pPlate->spatial_bin <= 0) ||
            (pPlate->spatial_bin_table == NULL) ||
            (pPlate->local_bin_index < 0)) {
          printf("ERROR: no mosaic or spatial bin in line %d\n",__LINE__);
          exit(-1);
        }
        if ((pPlate->spatial_bin < 1) || (pPlate->spatial_bin > MAX_SPATIAL_BINS)) {
          if (pPlate->fatal_reject_reason < 0) {
            printf("ERROR line %d for errorPoints %d\n",__LINE__,errorPoints);
          }
          errorPoints++;
          continue;
        }
        if (pPlate->spatial_bin != pPlate->spatial_bin_table[pPlate->spatial_bin].spatial_bin) {
          if (pPlate->fatal_reject_reason < 0) {
            printf("ERROR line %d for errorPoints %d\n",__LINE__,errorPoints);
          }
          errorPoints++;
          continue;
        }
        limitingPoints1++;
        if (goodFlag == 1) {
          continue;
        }
        /* Start with the iso upper limit */
        pPlate->magcal_iso = 0;
        pPlate->magcal_iso_rms = 0;
        pPlate->limiting_mag_local = pPlate->spatial_bin_table[pPlate->spatial_bin].limiting_mag;
        pPlate->magcal_local_rms = 0;
        pPlate->npoints_local = 0;
        pPlate->magcal_magdep = pPlate->spatial_bin_table[pPlate->spatial_bin].upper_limit;
        pPlate->magcal_local = pPlate->magcal_magdep;
        pPlate->magcal_iso = pPlate->magcal_magdep;
        
        if (pPlate->local_bin_table == NULL) {
          /* Need to allocate a local bin table now */
          pPlate->local_bin_table = (PLOCALBIN)calloc(pPlate->nx*pPlate->ny,sizeof(LOCALBIN));
          if (pPlate->local_bin_table == NULL) {
            printf("ERROR: Failed to allocate the local bin table \n");
            exit(-1);
          }
          for (local_bin_index = 0; local_bin_index < (pPlate->nx*pPlate->ny); local_bin_index++) {
            pLocalBin = &pPlate->local_bin_table[local_bin_index];
            pLocalBin->local_bin_index = -2;
          }
        }
        pLocalBin = &pPlate->local_bin_table[pPlate->local_bin_index];
        if (pLocalBin->local_bin_index == pPlate->local_bin_index) {
          pPlate->dradRMS2 = pLocalBin->dradRMS2;
          pPlate->npoints_local = pLocalBin->npoints_local;
          pPlate->rejectFlag = pLocalBin->rejectFlag;
          pPlate->extinction = pLocalBin->extinction;
          pPlate->magcor_local = pLocalBin->magcor_local;
          pPlate->magcal_local_error = pLocalBin->magcal_local_error;
          
        } else if (pLocalBin->local_bin_index == -2) {

          /* Here we need to read a bin of in the local bin table */
          if (GetPhotLocalBin(pPhotConnection,pLocalBin,1,pPlate->series,pPlate->plateNumber,pPlate->solutionNumber,pPlate->local_bin_index,catalogString,solutionNumberFix,0)) {
            printf("WARNING: Failed to read local bin table for %s%05d  seriesId = %d plateNumber = %d solutionNumber = %d\n",pPlate->series,pPlate->plateNumber,pPlate->seriesId,pPlate->plateNumber,pPlate->solutionNumber);
            pPlate->magcal_local_error = 99.0;
            pLocalBin->local_bin_index = -1;
            if (dumpMySQL) {
              printf("SELECT * from localbin where seriesId = %d and plateNumber = %d and local_bin_index = %d and exposureNumber = %d;\n",pPlate->seriesId,pPlate->plateNumber,pPlate->local_bin_index,pPlate->solutionNumber);

            } else {
              if (pPlate->fatal_reject_reason < 0) {
                pPlate->fatal_reject_reason = (FATAL_REASON_NOLOCALBINENTRY);
              }    
              pPlate->errorFlag = 1;

            }
          } else {
            (*localbinCount)++;
            pPlate->dradRMS2 = pLocalBin->dradRMS2;
            pPlate->npoints_local = pLocalBin->npoints_local;
            pPlate->rejectFlag = pLocalBin->rejectFlag;
            pPlate->extinction = pLocalBin->extinction;
            pPlate->magcor_local = pLocalBin->magcor_local;
            pPlate->magcal_local_error = pLocalBin->magcal_local_error;

          }
        } else {
          if (pPlate->fatal_reject_reason < 0) {
            pPlate->fatal_reject_reason = (FATAL_REASON_NOLOCALBINENTRY);
          }    
          pPlate->errorFlag = 1;
          pPlate->magcal_local_error = 99.0;
        }



        /* Now see if there is a local correction for this magnitude */
        if ((pPlate->magcal_local_error >= 90.0) ||
            (pPlate->rejectFlag != 0)) {
          if (pPlate->fatal_reject_reason < 0) {
            pPlate->fatal_reject_reason = (FATAL_REASON_LOCALBINFAILURE);
          }    
          localBinFailure = 1;
        } else {
          pPlate->limiting_mag_local -= pPlate->magcor_local;
          pPlate->magcal_magdep -= pPlate->magcor_local;
          pPlate->magcal_local -= pPlate->magcor_local;
        }
        /* Apply the magnitude-dependent correction for this plate */
        if (pPlate->pMagdepLimits == NULL) {
          /* We need to retrieve the magnitude-dependent correction table */
          pPlate->pMagdepLimits = (PMAGDEPLIMITS)calloc(1,sizeof(MAGDEPLIMITS));
          if (pPlate->pMagdepLimits == NULL) {
            printf("Error: failed to allocate pMagdepLimits\n");
            exit(-1);
          }
          PhotLoadMagdepCorrections(pPhotConnection,pPlate->seriesId,pPlate->plateNumber,pPlate->solutionNumber,catalogString,pPlate->pMagdepLimits,&pPlate->pMagdepTable,verbose);
          
        }
        if (pPlate->pMagdepLimits->nrecs > 0) {
          double temp_magdep_rms;
          double temp_magdep_magcor;
          int temp_magdep_bin;
          pPlate->limiting_mag_local -= GetMagdepBinMagcor(pPlate->pMagdepLimits,pPlate->pMagdepTable,pPlate->X_IMAGE,pPlate->Y_IMAGE,pPlate->limiting_mag_local,&temp_magdep_bin,&temp_magdep_rms,&temp_magdep_magcor,1);
          pPlate->magcal_magdep -= GetMagdepBinMagcor(pPlate->pMagdepLimits,pPlate->pMagdepTable,pPlate->X_IMAGE,pPlate->Y_IMAGE,pPlate->magcal_magdep,&temp_magdep_bin,&temp_magdep_rms,&temp_magdep_magcor,1);
          pPlate->magcal_local -= GetMagdepBinMagcor(pPlate->pMagdepLimits,pPlate->pMagdepTable,pPlate->X_IMAGE,pPlate->Y_IMAGE,pPlate->magcal_local,&temp_magdep_bin,&temp_magdep_rms,&temp_magdep_magcor,1);
        }


        if (pPlate->extinction > 0) {
          pPlate->magcal_magdep -= pPlate->extinction;
          pPlate->magcal_iso -= pPlate->extinction;
          pPlate->magcal_local -= pPlate->extinction;
          pPlate->limiting_mag_local -= pPlate->extinction;
        }
         

        if (localBinFailure) {
          pPlate->flag = 1.0;
          pPlate->reject_reason1 |= (1 << REJECT_REASON_BINFAILURE);
          pPlate->reject_reason2 |= (1 << REJECT_REASON_BINFAILURE);
          continue;     /* Here do not allow the points to clutter the plot */
        } else if (pPlate->spatial_bin == 9) {
          pPlate->flag = 31.5;
        } else {
          pPlate->flag = 31.0;
        }

        pPlate->magcal_iso_rms = 0.0;

        pPlate->Date = jd2hjd(pPlate->geoJulianDate,pTarget->ra,pTarget->dec,WCS_J2000);
        pPlate->validPoint = 1;
        validPoints++;
      }
    } else {
      if (pPlate->fatal_reject_reason < 0) {
        printf("ERROR line %d for errorPoints %d\n",__LINE__,errorPoints);
      }
      errorPoints++;
    }

  }
  if ((allLightcurves == 0) && (validPoints > 1)) {

    if (CalcMedianAndRMS(goodPoints,2,vector1,&rawmed,&rawrms,0,3.0,0) == 0) {
      rawmed = 0.0;
      rawrms = 99.0;
    }
    if (CalcMedianAndRMS(goodPoints,2,vector1,&clipmed,&cliprms,1,3.0,0) == 0) {
      clipmed = 0.0;
      cliprms = 99.0;
    }
    *pClipMed = clipmed;
    /* Find the zero-based clipped rms of our error bars */
    if (CalcMedianAndRMS(goodPoints,2,vector2,&rawerrmed,&rawerrrms,0,3.0,1) == 0) {
      rawerrmed = 0.0;
      rawerrrms = 99.0;
    }
    if (forceErrorBarFactor) {
      error_bar_factor = 1.0;
    } else {
      if ((rawrms == 0.0) || 
          (rawrms == 99.0) ||
          (rawerrrms == 0.0) || 
          (rawerrrms == 99.0)) {
        error_bar_factor = 1.0;
      } else {
        error_bar_factor = rawrms/rawerrrms;
      }
      if (error_bar_factor > 1.0) {
        error_bar_factor = 1.0;
      }
    }

    /* We have at least two good points, write the curve */
    strcpy(lightcurveName,lightcurveDirectory);
    strcat(lightcurveName,"/lc_");
    strcat(lightcurveName,pTarget->src_name);
    if (goodFlag) {
      strcat(lightcurveName,"_good");
    } else {
      strcat(lightcurveName,"_all");
    }
    strcat(lightcurveName,".db");
    lightcurveHandle = fopen(lightcurveName,"wt");
    if (lightcurveHandle == NULL) {
      printf("ERROR: Failed to open %s\n",lightcurveName);
      return(-1);
    }
    fprintf(lightcurveHandle,"sequence\thJD\tspatial_bin\tflag\tmagcal_magdep\tmagcal_local_rms\tnlocal\tPlate\tmagcal_iso_rms\tmagcal_local_error\tmagcal_iso\tlimiting_mag_local\tAFLAGS\treject_reason1\treject_reason2\tdradRMS2\tSFLAGS\tFWHM_IMAGE\tmagcal_local");
    if (fullFlag) {
      fprintf(lightcurveHandle,"\tBACKGROUND\tBlendedmag\tdec\tELLIPTICITY\textinction\tFLUX_ISO\tFLUX_MAX\tFWHM_WORLD\tISOAREA_WORLD\tKRON_RADIUS\tMAG_APER\tMAG_AUTO\tmagcor_local\tMAG_ISO\tplate_dist\tra\tStdmag\tTHETA_J2000\tX_IMAGE\tY_IMAGE\tBFLAGS\texposureNumber\tgsc_bin_index\tISO0\tISO1\tISO2\tISO3\tISO4\tISO5\tISO6\tISO7\tlocal_bin_index\tNUMBER\tplateNumber\trejectFlag\tseriesId\tsolutionNumber\tversionId\tmagdep_bin\tmagcal_magdep_rms\tcolorterm\terrorcolor\tcolorflag\tplateClass\tra_2\tdec_2\tRaPM\tDecPM\tA2FLAGS\tB2FLAGS\ttimeAccuracy\tmaskIndex\tREFNumber");

    }
    fprintf(lightcurveHandle,"\n");

    fprintf(lightcurveHandle,"--------\t---\t-----------\t----\t------------\t----------------\t------\t-----\t--------------\t-----------------\t----------\t------------------\t-----\t--------------\t--------------\t--------\t------\t----------\t------------");
    if (fullFlag) {
      fprintf(lightcurveHandle,"\t----------\t----------\t---\t-----------\t----------\t--------\t--------\t----------\t-------------\t-----------\t--------\t--------\t------------\t-------\t----------\t--\t------\t-----------\t-------\t-------\t------\t--------------\t-------------\t----\t----\t----\t----\t----\t----\t----\t----\t---------------\t------\t-----------\t----------\t--------\t--------------\t---------\t----------\t-----------------\t---------\t----------\t---------\t---------\t----\t-----\t----\t-----\t-----\t-------\t-------\t-----------\t----------\t---------");
    }
    fprintf(lightcurveHandle,"\n");


    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (pPlate->coverage_table_select == 0) {
        continue;
      }
      if ((pPlate->errorFlag == 0) &&
          (pPlate->validPoint != 0)) {
        endJD = pPlate->Date;
        if (startJD == 0.0) {
          startJD = endJD;
        }
        sequence++;
        if (pPlate->magcal_magdep < minMag) {
          minMag = pPlate->magcal_magdep;
        }
        if (pPlate->magcal_magdep > maxMag) {
          maxMag = pPlate->magcal_magdep;
        }
        if (pPlate->magcal_iso < minMag) {
          minMag = pPlate->magcal_iso;
        }
        if (pPlate->magcal_iso > maxMag) {
          maxMag = pPlate->magcal_iso;
        }
        if (pPlate->magcal_local < minMag) {
          minMag = pPlate->magcal_local;
        }
        if (pPlate->magcal_local > maxMag) {
          maxMag = pPlate->magcal_local;
        }
        totalPointsWritten++;
        if (pPlate->solutionNumber == 0) {
          total_s0PointsWritten++;
        } else {
          total_snPointsWritten++;
        }

        fprintf(lightcurveHandle,"%d\t%f\t%d\t%.1f\t%.2f\t%.2f\t%d\t%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%d\t%8.2f\t%.2f",
                sequence,
                pPlate->Date,
                pPlate->spatial_bin,
                pPlate->flag,
                pPlate->magcal_magdep,
                pPlate->magcal_local_rms*error_bar_factor,
                pPlate->npoints_local,
                pPlate->Plate,
                pPlate->magcal_iso_rms*error_bar_factor,
                pPlate->magcal_local_error*error_bar_factor,
                pPlate->magcal_iso,
                pPlate->limiting_mag_local,
                pPlate->AFLAGSCOPY,
                pPlate->reject_reason1,
                pPlate->reject_reason2,
                pPlate->dradRMS2,
                (pPlate->BFLAGS & FILTER_BMASK_SEXTRACTOR),
                pPlate->FWHM_IMAGE,
                pPlate->magcal_local);
        if (fullFlag) {
          fprintf(lightcurveHandle,"\t%12.7f\t%.2f\t%11.7f\t%8.3f\t%.3f\t%12.7f\t%12.7f\t%12.7f\t%12.7f\t%5.2f\t%8.4f\t%8.4f\t%f\t%8.4f\t%.4f\t%11.7f\t%.2f\t%6.2f\t%10.3f\t%10.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%f\t%f\t%d\t%s\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%d\t%lld",
                  pPlate->BACKGROUND,
                  pPlate->Blendedmag,
                  pPlate->dec,
                  pPlate->ELLIPTICITY,
                  pPlate->extinction,
                  pPlate->FLUX_ISO,
                  pPlate->FLUX_MAX,
                  pPlate->FWHM_WORLD,
                  pPlate->ISOAREA_WORLD,
                  pPlate->KRON_RADIUS,
                  pPlate->MAG_APER,
                  pPlate->MAG_AUTO,
                  pPlate->magcor_local,
                  pPlate->MAG_ISO,
                  pPlate->plate_dist,
                  pPlate->ra,
                  pPlate->Stdmag,
                  pPlate->THETA_J2000,
                  pPlate->X_IMAGE,
                  pPlate->Y_IMAGE,
                  pPlate->BFLAGS,
                  pPlate->exposureNumber,
                  pPlate->gsc_bin_index,
                  pPlate->ISO0,
                  pPlate->ISO1,
                  pPlate->ISO2,
                  pPlate->ISO3,
                  pPlate->ISO4,
                  pPlate->ISO5,
                  pPlate->ISO6,
                  pPlate->ISO7,
                  pPlate->local_bin_index,
                  pPlate->NUMBER,
                  pPlate->plateNumber,
                  pPlate->rejectFlag,
                  pPlate->seriesId,
                  pPlate->solutionNumber,
                  pPlate->versionId,
                  pPlate->magdep_bin,
                  pPlate->magcal_magdep_rms,
                  pPlate->colorterm,
                  pPlate->errorcolor,
                  pPlate->colorflag,
                  pPlate->plateClass,
                  pPlate->ra_2,
                  pPlate->dec_2,
                  pPlate->RaPM,
                  pPlate->DecPM,
                  pPlate->A2FLAGS,
                  pPlate->B2FLAGS,
                  pPlate->timeAccuracy,
                  pPlate->maskIndex,
                  pPlate->REFNumber
	  
                  );


        }
        fprintf(lightcurveHandle,"\n");

      }

    }
    fclose(lightcurveHandle);
    /* Now write out the summary text */
    rangeJD = (endJD-startJD);
    endJD += MARGIN_WIDTH * rangeJD;
    startJD -= MARGIN_WIDTH *rangeJD;


    startYear = jd2ep(startJD);
    endYear = jd2ep(endJD);
    strcpy(lightcurveName,lightcurveDirectory);
    strcat(lightcurveName,"/lc_");
    strcat(lightcurveName,pTarget->src_name);
    if (goodFlag) {
      strcat(lightcurveName,"_good");
    } else {
      strcat(lightcurveName,"_all");
    }
    strcat(lightcurveName,".txt");
    lightcurveHandle = fopen(lightcurveName,"wt");
    if (lightcurveHandle == NULL) {
      printf("ERROR: Failed to open %s\n",lightcurveName);
      return(-1);
    }
    fprintf(lightcurveHandle,"#raw median mag, raw median rms, clipped median mag, clipped median rms, total points, good points, start year, end year,error_bar_factor, start JD,end JD,REF,ra,dec,Stdmag,color,VFlag,MAGFlag,RaPM,DecPM,minMag,maxMag\n");

    fprintf(lightcurveHandle,"%.2f %.2f %.2f %.2f %d %d %f %f %f %.0f %.0f ",
            rawmed,rawrms,clipmed,cliprms,validPoints,goodPoints,startYear,endYear,error_bar_factor,startJD,endJD);

    if (pTarget->gsc_bin_index != 0) {
      ra2str (rstr, 16, pTarget->ra, ndec);
      dec2str (dstr, 16, pTarget->dec, ndec-1);    
      if (pTarget->updateflag == UPDATEFLAG_NO_UPDATE) {

        if ((pTarget->MAGFlag & KEPLER_MAGNITUDE_FLAG) == 0) {
          sprintf(MAGFlagString,"%d",pTarget->MAGFlag);
        } else {
          if ((pTarget->MAGFlag & MAGNITUDE_MASK) < KEPLER_CQ_MAGFLAG_MAX) {
            strcpy(MAGFlagString,keplerSourceText[pTarget->MAGFlag & MAGNITUDE_MASK]);
          } else {
            sprintf(MAGFlagString,"%d??",pTarget->MAGFlag);
          }
        }


        fprintf(lightcurveHandle,"%s %s %s %.2f %.2f %d %s %.1f %.1f ",
                pTarget->REF,
                rstr,
                dstr,
                pTarget->Stdmag,
                pTarget->color,
                pTarget->VFlag,
                MAGFlagString,
                pTarget->RaPM,
                pTarget->DecPM);
      } else {
        fprintf(lightcurveHandle,"%s %s %s UNK UNK UNK UNK UNK UNK ",
                pTarget->REF,
                rstr,
                dstr);

      }
    } else {
      fprintf(lightcurveHandle,"UNK UNK UNK UNK UNK UNK UNK UNK UNK ");

    }
    fprintf(lightcurveHandle,"%.2f %.2f\n",minMag,maxMag);

    fclose(lightcurveHandle);
    if (allLightcurves == 0) {
      printf("Lightcurve written goodFlag %d good: %4d err: %4d bad: %4d limiting %4d %4d bin9 %4d both %4d wedge %4d high drad %4d defect %4d hiZout %4d tooBright %d too low %d blend %4d errbarfactor %f for %s %s modulus %03d gsc_bin_index %d\n",goodFlag,goodPoints,errorPoints,badPoints,limitingPoints1,limitingPoints2,bin9Points,bothPoints,wedgeCount,highDradCount,defectCount,hiZoutCount,tooBrightCount,lowAltitudeCount,blendCount,error_bar_factor,pTarget->src_name,pTarget->REF,pTarget->gsc_bin_index % MAGNITUDES_MODULUS, pTarget->gsc_bin_index);
    }
    curveWritten = 1;



  } else {
    if (verbose & (allLightcurves == 0)) {
      printf("No valid points found for %s\n",pTarget->src_name);
    }
  }
             
  if (printReject) {
    char rejectFileName[MAX_BUFFER];
    FILE *rejectFileHandle = NULL;
    char undetectedName[MAX_BUFFER];
    FILE *undetectedHandle = NULL;
    char fatalCountFileName[MAX_BUFFER];
    FILE *fatalCountFileHandle = NULL;
    char fatalPlateFileName[MAX_BUFFER];
    FILE *fatalPlateFileHandle = NULL;
    double fatal_s0Percent;
    double fatal_snPercent;
    strcpy(rejectFileName,lightcurveDirectory);
    strcat(rejectFileName,"/lc_");
    strcat(rejectFileName,pTarget->src_name);
    strcat(rejectFileName,"_reject_reason.txt");
    rejectFileHandle = fopen(rejectFileName,"wt");
    if (rejectFileHandle == NULL) {
      printf("ERROR: Failed to open %s\n",rejectFileName);
      return(-1);
    }
    strcpy(undetectedName,lightcurveDirectory);
    strcat(undetectedName,"/lc_");
    strcat(undetectedName,pTarget->src_name);
    strcat(undetectedName,"_undetected.txt");
    undetectedHandle = fopen(undetectedName,"wt");
    if (undetectedHandle == NULL) {
      printf("ERROR: Failed to open %s\n",undetectedName);
      return(-1);
    }
    strcpy(fatalCountFileName,lightcurveDirectory);
    strcat(fatalCountFileName,"/lc_");
    strcat(fatalCountFileName,pTarget->src_name);
    strcat(fatalCountFileName,"_fatal_reason.txt");
    fatalCountFileHandle = fopen(fatalCountFileName,"wt");
    if (fatalCountFileHandle == NULL) {
      printf("ERROR: Failed to open %s\n",fatalCountFileName);
      return(-1);
    }

    strcpy(fatalPlateFileName,lightcurveDirectory);
    strcat(fatalPlateFileName,"/lc_");
    strcat(fatalPlateFileName,pTarget->src_name);
    strcat(fatalPlateFileName,"_fatal_plate.txt");
    fatalPlateFileHandle = fopen(fatalPlateFileName,"wt");
    if (fatalPlateFileHandle == NULL) {
      printf("ERROR: Failed to open %s\n",fatalPlateFileName);
      return(-1);
    }



    fprintf(undetectedHandle,"Plate\treject_reason\n");
    fprintf(undetectedHandle,"-----\t-------------\n");

    /* Tally up the individual counts */
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (pPlate->coverage_table_select == 0) {
        continue;
      }
      if (pPlate->reject_reason2 & ((1<<REJECT_REASON_ERROR)|
                                    (1<<REJECT_REASON_OFFPLATE)|
                                    (1<<REJECT_REASON_BINFAILURE)|
                                    (1<<REJECT_REASON_UNDETECTED))) {
        fprintf(undetectedHandle,"%s\t%d\n",pPlate->Plate,pPlate->reject_reason2);
      }
      for (rejectIndex = 0; rejectIndex < REJECT_REASON_MAX; rejectIndex++) {
        if ((pPlate->reject_reason1 & (1 << rejectIndex)) != 0) {
          reject_count_all1[rejectIndex]++;
          if (pPlate->reject_reason1 == (1 << rejectIndex)) {
            reject_count_only1[rejectIndex]++;
          }
        }
        if ((pPlate->reject_reason2 & (1 << rejectIndex)) != 0) {
          reject_count_all2[rejectIndex]++;
          if (pPlate->reject_reason2 == (1 << rejectIndex)) {
            reject_count_only2[rejectIndex]++;
          }
        }
      }
      if (pPlate->solutionNumber == 0) {
        total_s0Count++;
      } else {
        total_snCount++;
      }

      if ((pPlate->fatal_reject_reason < 0) &&
          (pPlate->errorFlag != 0)) {
        pPlate->fatal_reject_reason = FATAL_REASON_BADERRORFLAG;
      }
      if ((pPlate->fatal_reject_reason < 0) &&
          (pPlate->validPoint == 0)) {
        pPlate->fatal_reject_reason = FATAL_REASON_BADVALIDPOINT ;
      }

      if (pPlate->fatal_reject_reason < 0) {
        if (pPlate->solutionNumber == 0) {
          fatal_reject_count_plate_s0[FATAL_REASON_NONE]++;          
        } else {
          fatal_reject_count_plate_sn[FATAL_REASON_NONE]++;
        }
      } else {
        if (pPlate->solutionNumber == 0) {
          fatal_reject_count_plate_s0[pPlate->fatal_reject_reason]++;
          fatal_s0PointCount++;
        } else {
          fatal_reject_count_plate_sn[pPlate->fatal_reject_reason]++;
          fatal_snPointCount++;
        }
      }
    }
    fprintf(rejectFileHandle,"Reject reasons of %d plates for %s %s\n",
            plateListLength,
            pTarget->REF,
            pTarget->src_name);
    fprintf(rejectFileHandle,"Column 1: Reject reason index (See REJECT_REASON in pipelineutils.h)\n");
    fprintf(rejectFileHandle,"Column 2: Total of all plates affected by the condition\n");
    fprintf(rejectFileHandle,"Column 3: Total of all plates exclusively affected by the condition\n");
    fprintf(rejectFileHandle,"Column 4: Number of plates plotted in 'good' lightcurve affected by the condition\n");
    fprintf(rejectFileHandle,"Column 5: Number of plates plotted in 'good' lightcurve exclusively affected by the condition\n");
    fprintf(rejectFileHandle,"ignoreBin %d, AFLAGSMASK %d BFLAGSMASK %d (Use showflags <aflagsmask> to show included error bits\n",ignoreBin,AFLAGSMASK,BFLAGSMASK);


    for (rejectIndex = 0; rejectIndex < REJECT_REASON_MAX; rejectIndex++) {
      fprintf(rejectFileHandle,"%4d %4d %4d %4d %4d %s\n",
              rejectIndex,
              reject_count_all1[rejectIndex],
              reject_count_only1[rejectIndex],
              reject_count_all2[rejectIndex],
              reject_count_only2[rejectIndex],
              GetRejectReasonText(rejectIndex));
    }
    fatal_s0Total = 1.0*(total_s0Count-fatal_reject_count_plate_s0[FATAL_REASON_OFFPLATE]);
    fatal_snTotal = 1.0*(total_snCount-fatal_reject_count_plate_sn[FATAL_REASON_OFFPLATE]);

    for (rejectIndex = 0; rejectIndex < FATAL_REASON_MAX; rejectIndex++) {
      if ((fatal_s0Total == 0) || (rejectIndex == FATAL_REASON_OFFPLATE)) {
        fatal_s0Percent = 0;
      } else {
        fatal_s0Percent = (100.0*fatal_reject_count_plate_s0[rejectIndex])/(fatal_s0Total);
      }
      if ((fatal_snTotal == 0) || (rejectIndex == FATAL_REASON_OFFPLATE)) {
        fatal_snPercent = 0;
      } else {
        fatal_snPercent = (100.0*fatal_reject_count_plate_sn[rejectIndex])/(fatal_snTotal);
      }
      fprintf(fatalCountFileHandle,"Fatal reject index %2d count s0: %6d %5.1f sn: %6d %5.1f %s\n",
              rejectIndex,
              fatal_reject_count_plate_s0[rejectIndex],
              fatal_s0Percent,
              fatal_reject_count_plate_sn[rejectIndex],
              fatal_snPercent,
              GetFatalReasonText(rejectIndex));
      printf("Fatal reject index %2d count s0: %6d %5.1f sn: %6d %5.1f %s\n",
              rejectIndex,
              fatal_reject_count_plate_s0[rejectIndex],
              fatal_s0Percent,
              fatal_reject_count_plate_sn[rejectIndex],
              fatal_snPercent,
              GetFatalReasonText(rejectIndex));
    }
    fprintf(fatalCountFileHandle,"totalExposureCount s0: %d sn: %d fatalPointCount s0: %d sn: %d, totalPointsWritten s0: %d sn: %d total %d\n",
            total_s0Count,
            total_snCount,
            fatal_s0PointCount,
            fatal_snPointCount,
            total_s0PointsWritten,
            total_snPointsWritten,
            totalPointsWritten);
    printf("totalExposureCount s0: %d sn: %d fatalPointCount s0: %d sn: %d, totalPointsWritten s0: %d sn: %d total %d\n",
            total_s0Count,
            total_snCount,
            fatal_s0PointCount,
            fatal_snPointCount,
            total_s0PointsWritten,
            total_snPointsWritten,
            totalPointsWritten);

    for (rejectIndex = 0; rejectIndex < FATAL_REASON_MAX; rejectIndex++) {
      for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
        pPlate = &plateArray[plateIndex];
        if (pPlate->fatal_reject_reason != rejectIndex) {
          continue;
        }
        if ((pPlate->fatal_reject_reason >= 0) &&
            (pPlate->fatal_reject_reason != FATAL_REASON_NONE) &&
            (pPlate->fatal_reject_reason != FATAL_REASON_OFFPLATE)) {
          

          fprintf(fatalPlateFileHandle,"%s Reason %d %s",pPlate->Plate,pPlate->fatal_reject_reason,GetFatalReasonText(pPlate->fatal_reject_reason));
          if (pPlate->pMosaic == NULL) {
            fprintf(fatalPlateFileHandle," NO PHOTOMETRY PLATE ENTRY");
          } else if ((pPlate->pMosaic->FitWCS & FITWCS_NOALLOBJECTS) != 0) {
            fprintf(fatalPlateFileHandle," NO GSC PHOTOMETRY RESULTS");
          }
          if (pPlate->fatal_reject_reason == FATAL_REASON_LOCALBINFAILURE) {
            fprintf(fatalPlateFileHandle," Local bin (0-49) RA: %2d Dec: %2d",pPlate->local_bin_index % pPlate->nx,pPlate->local_bin_index /pPlate->ny);

          }
          fprintf(fatalPlateFileHandle,"\n");
        }
      }
    }

    fclose(rejectFileHandle);
    fclose(fatalCountFileHandle);
    fclose(undetectedHandle);
    fclose(fatalPlateFileHandle);
  }



  return(curveWritten);
}
int TargetCompare(const void *first, const void *second) 
{
  int gsc_bin_indexFirst = ((PPHOTTARGET)first)->gsc_bin_index;
  int gsc_bin_indexSecond = ((PPHOTTARGET)second)->gsc_bin_index;
  if (gsc_bin_indexFirst > gsc_bin_indexSecond) {
    return(1);
  } else if (gsc_bin_indexFirst < gsc_bin_indexSecond) {
    return(-1);
  } else {
    return(0);
  }
}



int GscBinCompare(const void *first, const void *second) 
{
  int coverage_bin_indexFirst = *((int *)first);
  int coverage_bin_indexSecond = *((int *)second);
  if (coverage_bin_indexFirst > coverage_bin_indexSecond) {
    return(1);
  } else if (coverage_bin_indexFirst < coverage_bin_indexSecond) {
    return(-1);
  } else {
    return(0);
  }
}
int ConeSearch(int **pgsc_bin_table,int *pgsc_bin_alloc,int *pgsc_bin_count,double ra0,double dec0,double rad0,int verbose)
{
#if 0
intMYSQL *pPhotConnection,PPHOTTARGET *pTarget_table,size_t *pTarget_nrecs,int dumpMySQL,char *catalogString
#endif
  int totalBins = 0;
  int *gsc_bin_table = NULL;
  int gsc_bin_index;
  int alloc_bins;
  double rad1;
  int minDecBin;
  int maxDecBin;
  double dec;
  double ra;
  double decHigh;
  double decLow;
  double raHigh;
  double raLow;
  int raBin;
  int raBin0;
  int decBin;
  int wrapFlag;
  int binIndex;
  PBININDEX pBinIndex;
#if 0

  PPHOTTARGET target_table = NULL;
  PPHOTTARGET tmp_target_table = NULL;
  PPHOTTARGET pTarget = NULL;
  int target_alloc = 0;
  size_t target_nrecs = 0;
  int curStars;
  int starIndex;
  PSTARENTRY pStarTable = NULL;
  PSTARENTRY pStarEntry = NULL;
  double radactual;
  char tmp_src_name[MAX_SRC_LENGTH];
  int decCBin;
  int raCBin;
  long long REFNumber;
#endif

  *pgsc_bin_table = NULL;
  *pgsc_bin_alloc = 0;
  *pgsc_bin_count = 0;
  

  rad1 = rad0 + (BIN_EXPANSION_RADIUS * 3600) + (2*pGscBin->bin_size*3600);

  alloc_bins = 2.0*PI_VALUE*sqr(2.0+(rad1/(3600.0*pGscBin->bin_size)));
  gsc_bin_table = (int *)calloc(alloc_bins,sizeof(int));
  


  if (gsc_bin_table == NULL) {
    printf("ERROR: Failed to allocate gsc_bin_table for %d bins\n",alloc_bins);
    exit(-1);
  }
  *pgsc_bin_table = gsc_bin_table;
  *pgsc_bin_alloc = alloc_bins;

  gsc_bin_index = GetGSCBin(pGscBin,ra0,dec0,&minDecBin,&raBin,"ConeSearch");
  printf("Cone Search centered around GSC bin %d\n",gsc_bin_index);

  /* Find the lower search limit */

  dec = dec0-(rad1/3600.);
  if (dec < -90.0) {
    dec = -90.0;
  }
  gsc_bin_index = GetGSCBin(pGscBin,ra0,dec,&minDecBin,&raBin,"ConeSearch");

  /* Now find the uppper limit */
  dec = dec0+(rad1/3600.);
  if (dec >= +90.0) {
    dec = +90.0-(0.1*pGscBin->bin_size);
  }
  gsc_bin_index = GetGSCBin(pGscBin,ra0,dec,&maxDecBin,&raBin,"ConeSearch");
  

  /* Search from South to North */
  for (decBin = minDecBin; decBin <= maxDecBin; decBin++) {
    pBinIndex = &pGscBin->pBinMasterIndex[decBin];
    decLow = (decBin * pGscBin->bin_size) - 90.0;
    decHigh = ((decBin+1) * pGscBin->bin_size) - 90.0;
    raBin0 = ((ra0 * pBinIndex->numBins) / 360.0);
    /* Search to the West */
    wrapFlag = 0;
    raBin = raBin0;
    while(1) {
      if (raBin >= pBinIndex->numBins) {
        raBin = 0;
      }
      raLow = (raBin * 360.0)/pBinIndex->numBins;
      raHigh = ((raBin+1) * 360.0)/pBinIndex->numBins;
      if (((wcsdist(ra0,dec0,raLow ,decLow )*3600) < rad1) ||
          ((wcsdist(ra0,dec0,raHigh,decLow )*3600) < rad1) ||
          ((wcsdist(ra0,dec0,raLow ,decHigh)*3600) < rad1) ||
          ((wcsdist(ra0,dec0,raHigh,decHigh)*3600) < rad1)) {
        /* Add this one to our table */
        if (totalBins < alloc_bins) {
          gsc_bin_table[totalBins] = pBinIndex->startBin+raBin;
          totalBins++;
        } else {
          printf("ERROR totalBins %d exceeds allocate %d\n",totalBins,alloc_bins);
          exit(-1);
        }


      } else {
        break;
      }
      raBin++;
      if (raBin == raBin0) {
        /* We wrapped! */
        wrapFlag = 1;
        break;
      }
    }
     
    /* Search to the East */
    if (wrapFlag == 0) {
      raBin = raBin0 -1;
      while(1) {
        if (raBin < 0) {
          raBin = pBinIndex->numBins -1;
        }
        raLow = (raBin * 360.0)/pBinIndex->numBins;
        raHigh = ((raBin+1) * 360.0)/pBinIndex->numBins;
        if (((wcsdist(ra0,dec0,raLow ,decLow )*3600) < rad1) ||
            ((wcsdist(ra0,dec0,raHigh,decLow )*3600) < rad1) ||
            ((wcsdist(ra0,dec0,raLow ,decHigh)*3600) < rad1) ||
            ((wcsdist(ra0,dec0,raHigh,decHigh)*3600) < rad1)) {
          /* Add this one to our table */
          if (totalBins < alloc_bins) {
            gsc_bin_table[totalBins] = pBinIndex->startBin+raBin;
            totalBins++;
          } else {
            printf("ERROR totalBins %d exceeds allocate %d\n",totalBins,alloc_bins);
            exit(-1);
          }


        } else {
          break;
        }
        raBin--;
      }
    }
  }
  
  if (verbose) {
    printf("Searching through %d gsc bins\n",totalBins);
  }


  /* For performance, sort in increasing bin order */
  qsort((void*)gsc_bin_table,totalBins,sizeof(int),GscBinCompare);
  for (binIndex = 1; binIndex < totalBins; binIndex++) {
    if (gsc_bin_table[binIndex-1] == gsc_bin_table[binIndex]) {
      printf("ERROR: duplicate index in gsc_bin_table\n");
      exit(-1);
    }
  }
  *pgsc_bin_count = totalBins;
  return(0);
}

int GetCoverage(MYSQL * pAuxConnection,PCOVERAGE *pcoverage_table,int *pcoverage_table_size,int coverage_bin_index,int excludeSeriesCount)
{
  char queryString[MAX_QUERY_STRING];
  char tmpString[MAX_QUERY_STRING];
  int nvals;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  PCOVERAGE coverage_table = *pcoverage_table;
  int coverage_table_size = *pcoverage_table_size;
  int result_size;
  int foundCount = 0;
  PCOVERAGE pCoverage;
  int excludeSeriesIndex;

  *pcoverage_table_size = 0;
  sprintf(queryString,"SELECT series,plateNumber,solutionNumber from coverage where coverage_bin_index = %d",coverage_bin_index);
  for (excludeSeriesIndex = 0; excludeSeriesIndex < excludeSeriesCount; excludeSeriesIndex++) {
    sprintf(tmpString," and series != '%s'",excludeSeriesList[excludeSeriesIndex]);
    strcat(queryString,tmpString);
  }
  strcat(queryString," order by series,plateNumber,solutionNumber;");
  res = ExecuteQuery(pAuxConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pAuxConnection);
    if (res_ptr) {
      result_size = (int)mysql_num_rows(res_ptr);
      if (result_size > 0) {
        if (coverage_table != NULL) {
          free(coverage_table);
        }
        coverage_table = (PCOVERAGE)calloc(result_size,sizeof(COVERAGE));
        *pcoverage_table = coverage_table;
        if (coverage_table == NULL) {
          printf("ERROR: allocation failure for coverage_table of size %d\n",result_size);
          exit(-1);
        }

      }
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        pCoverage = &coverage_table[foundCount];
        if (sqlrow[0] == NULL) {
          printf("ERROR: series is NULL in coverage table\n");
          continue;
        }
        strcpy(pCoverage->series,sqlrow[0]);
        nvals = sscanf(sqlrow[1],"%d",&pCoverage->plateNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for plateNumber in coverage table \n",nvals);
          continue;
        }
        nvals = sscanf(sqlrow[2],"%d",&pCoverage->solutionNumber);
        if (nvals != 1) {
          printf("ERROR: nvals is %d for solutionNumber in coverage table \n",nvals);
          continue;
        }
        pCoverage->coverage_bin_index = coverage_bin_index;
        foundCount++;
#if 0
        if ((pCoverage->plateNumber == 1146) &&
            (strcmp("rh",pCoverage->series) == 0)) {
          printf("At plate %s%05d coverage_bin_index %d\n",pCoverage->series,pCoverage->plateNumber,coverage_bin_index);
        }
#endif
        if (foundCount > result_size) {
          printf("ERROR: foundCount %d exceeds allocCount %d\n",foundCount,result_size);
          exit(-1);
        }
      }
      mysql_free_result(res_ptr);
    }
  } else {
    printf("Select error %d: %s res %d\n",mysql_errno(pAuxConnection),mysql_error(pAuxConnection),res);
    return(0);
  } 
  if (foundCount != result_size) {
    printf("ERROR: foundCount %d not equal to allocCount %d\n",foundCount,result_size);
    return(0);
  }
  *pcoverage_table_size = result_size;
  return(1);
}

void ReInitializePlate(PPLATE pPlate) {
  pPlate->starDetected = 0;
  pPlate->validPoint = 0;
  pPlate->reject_reason1 = 0;
  pPlate->reject_reason2 = 0;
  pPlate->spatial_bin = 0;
  pPlate->local_bin_index = -1;
  pPlate->fatal_reject_reason= -1;
  pPlate->BACKGROUND = 0;
  pPlate->Blendedmag = 0;
  pPlate->dec = 0;
  pPlate->ELLIPTICITY = 0;
  pPlate->FLUX_ISO = 0;
  pPlate->FLUX_MAX = 0;
  pPlate->FWHM_WORLD = 0;
  pPlate->ISOAREA_WORLD = 0;
  pPlate->KRON_RADIUS = 0;
  pPlate->MAG_APER = 0;
  pPlate->MAG_AUTO = 0;
  pPlate->MAG_ISO = 0;
  pPlate->plate_dist = 0;
  pPlate->ra = 0;
  pPlate->THETA_J2000 = 0;
  pPlate->X_IMAGE = 0;
  pPlate->Y_IMAGE = 0;
  pPlate->exposureNumber = 0;
  pPlate->ISO0 = 0;
  pPlate->ISO1 = 0;
  pPlate->ISO2 = 0;
  pPlate->ISO3 = 0;
  pPlate->ISO4 = 0;
  pPlate->ISO5 = 0;
  pPlate->ISO6 = 0;
  pPlate->ISO7 = 0;
  pPlate->magdep_bin = 0;
  pPlate->magcal_magdep_rms = 0;
  pPlate->REFNumber = 0;
  pPlate->colorterm = 0;
  pPlate->errorcolor = 0;
  pPlate->colorflag = 0;

  pPlate->magcal_magdep = 0;
  pPlate->magcal_local = 0;
  pPlate->Date = 0;
  pPlate->magcal_local_rms = 0;
  pPlate->Stdmag = 0;
  pPlate->magcal_iso_rms = 0;
  pPlate->magcal_iso = 0;
  pPlate->limiting_mag_local = 0;
  pPlate->npoints_local = 0;
  pPlate->AFLAGS = 0;
  pPlate->AFLAGSCOPY = 0;
  pPlate->BFLAGS = 0;
  pPlate->gsc_bin_index = 0;
  pPlate->dradRMS2 = 0;
  pPlate->FWHM_IMAGE = 0;
  pPlate->rejectFlag = 0;
  pPlate->extinction = 0;
  pPlate->magcor_local = 0;
  pPlate->magcal_local_error = 0;
  pPlate->NUMBER = 0;
  pPlate->flag = 0;
  pPlate->ra_2 = 999.0;
  pPlate->dec_2 = 99.0;
  pPlate->RaPM = NULL_PROPER_MOTION;
  pPlate->DecPM = NULL_PROPER_MOTION;
  pPlate->A2FLAGS = 0;
  pPlate->B2FLAGS = 0;
  pPlate->timeAccuracy = -1;
  pPlate->maskIndex = 0;
  /* The following fields appear to be unused here */
#if 0
  EXPOSURE exposure;
  int nx;              
  int ny;              
  double clipmed;
  PPHOTTARGET pTarget; 
#endif

}

int main(int argc,char *argv[])
{
  int nvals;
  char *argstr2;
  char cmdchar2;
  char *argstr;
  char cmdchar;
  char lightcurveDirectory[MAX_BUFFER];
  char rmsName[MAX_BUFFER];
  FILE *rmsHandle = NULL;
  char regionFileName[MAX_BUFFER];
  FILE *regionFileHandle = NULL;

  char *dotPos;
  char *slashPos;
  char *lastSlashPos;
 
  time_t startTime;
  time_t curTime;
  int verbose = 0;
  int forceErrorBarFactor = 0;
  int skipWrite = 0;
  int fullFlag = 0;
  int dumpMySQL= 0;
  int printReject = 0;
  int debugFlag = 0;
  int outputPlotCount = 0;
  int spatialBinCount = 0;
  int localbinCount = 0;
  int staleEntryCount = 0;
  int wcsCount = 0;
  int outputResult;
  int goodFlag = 0;
  
  int lineLen;
  int nlines = 0;
  char *inBuffer;
  PPLATE plateArray = NULL;
  PPLATE tmpPlateArray;
  double *vector1 = NULL;
  double *vector2 = NULL;
  int plateListLength = 0;
  int duplicatePlateEntries = 0;
  int plateListAlloc = 0;
  int compareFlag;
  int match_count = 0;
  PPLATE pPlate;
  PPLATE pPrevPlate;
  int plateIndex;
  double clipmed;
  PFILESTARIMAGE pNoneMagnitudeTable = NULL;
  PFILESTARIMAGE pNoneImage;
  int noneMagnitudeAlloc = 0;


  int ignoreBin = 0;

  File target_handle = NULL;
  char target_name[MAX_BUFFER];
  char catalog_name[MAX_BUFFER];
  FILE *catalog_handle = NULL;
  TableHead target_header = NULL;
  PPHOTTARGET target_table = NULL; 
  PPHOTTARGET tmp_target_table = NULL;
  size_t target_nrecs = 0;
  size_t target_alloc = 0;
  size_t base_target_nrecs = 0;
  int target_index;
  int target_index2;
  int target_index3;
  PPHOTTARGET pTarget;
  PPHOTTARGET pTarget2;
  PPHOTTARGET pTarget3;
  PHOTTARGET blendTarget;

  PPHOTTARGET pLastTarget = NULL;


  int master_nrecs = 0;
  int master_index;
  PPHOTSTARIMAGE pMaster = NULL;

  int errorFlag = 0;

  int blendCount = 0;
  char blendStrings[MAX_BLEND_COUNT*MAX_BUFFER];
  PPLATE blendArray[MAX_BLEND_COUNT];
  int savedBlendCount = 0;
  int blendIndex;
  int writePointResult;
  int decBin;
  int raBin;
  int decCBin;
  int raCBin;
  
  char *mysqlhost;
  char *username;
  char *password;
  char *auxscanner;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL aux_connection;
  MYSQL *pAuxConnection = &aux_connection;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;

  int curStars;
  STARENTRY starTable;
  PSTARENTRY pStarTable = &starTable;
  char tmp_src_name[MAX_SRC_LENGTH];
  double radactual;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  int gotAnswer;
  int curMagnitudes;
  int ignoreVersionId = 1;
  int staleVersionIdCount = 0;


  PPHOTSTARIMAGE pMagnitudeTable = NULL;
  PPHOTSTARIMAGE pCurSumstarimage;
  PFILESTARIMAGE pFileStarImage = NULL;
  int magnitudeAlloc = 0;
  int magnitudeIndex;
  int numMagnitudes;
  int vectorAlloc = 0;

  int allLightcurves = 0;
  int old_coverage_bin_index;
  int old_gsc_bin_index;
  int totalCoverageBins = 0;
  int coverageBinIndex;
  PCOVERAGE *coverage_table;
  PCOVERAGE pCoverage;
  int *coverage_table_size;
  int coverageIndex;
  PLOCALBIN pLocalBin;
  int local_bin_index;
  int coverageTableCount = 0;
  int maxref = 0;
  int curref;
  int coneSearchFlag = 0;
  long long *REFNumberArray = NULL;
  long long tmpREFNumber;
  char tmpREF[MAX_REF];
  double rad0 = -1.0; /* Search radius in arcsec */
  double ra0 = -99.0;	/* Initial center RA in degrees */
  double dec0 = -99.0;	/* Initial center Dec in degrees */
  int cFlag = 0;
  int pFlag = 0;
  int tFlag = 0;
  char indexname[MAX_BUFFER];
  File indexHandle;
  char* catalogdir;
  char* slashPtr;
  char catalogname[MAX_BUFFER];
  File catalogHandle;
  GSC_CACHE_ENTRY pGscCache[GSC_CACHE_ENTRIES];
  int newGscBinIndex;
  PGSCIMAGE pCurGscImage;
  char * dotPtr;
  int statResult;
  struct stat statbuf;
  int AFLAGSMASK = -1;
  int BFLAGSMASK = -1;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int coverageErrorCount = 0;
  char magnitudeFileName[MAX_FILENAME];
  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon = &fileCommon;
  int solutionNumberFix = 0;
  int minimumImageCount = 1;
  int histogramTable[MAX_NONE_HISTOGRAM];
  int groupNumber = 1;
  int excludeSeriesCount = 0;
  int excludeSeriesIndex;
  int excludeSeriesFlag;
  int imageExcludeCount = 0;
  int refType;
  int readOnly = 1;
  int *gsc_bin_table = NULL;
  int *coverage_bin_for_gsc_bin_table = NULL;
  int *coverage_bin_table = NULL;
  int gsc_bin_alloc = 0;
  int gsc_bin_count = 0;
  int gsc_bin_iterator;
  int gsc_bin_index;
  int coverage_bin_index;
  int coverage_bin_iterator;
  
  double raCenter;
  double decCenter;
  double ra;
  double dec;


  memset(pGscCache,0,sizeof(pGscCache));

  time(&startTime);
  catalogString[0] = 0;
  target_name[0] = 0;
  catalog_name[0] = 0;
  lightcurveDirectory[0] = 0;
  SetQueryCount(0);
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




        case 'b': /* blended object name */
        case 'B':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            if (blendCount < MAX_BLEND_COUNT) {
              strncpy(&blendStrings[blendCount*MAX_BUFFER],*++argv,MAX_BUFFER-2);
              blendCount++;
              if (strlen(*argv) >= MAX_BUFFER-2) {
                printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
                errorFlag = 1;
              }
            } else {
              printf("ERROR: More than %d blend stars\n",blendCount);
              errorFlag = 1;                          
            }
          }
          break;

        case 't': /* target catalog file name */
        case 'T':
          tFlag = 1;
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(target_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'n': /* new catalog file name */
        case 'N':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(catalog_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
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
            strncpy(lightcurveDirectory,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'c':  /* GSC2.3.2 reference object list */
        case 'C':
          cFlag = 1;
          maxref = argc-1;
          for (curref = 1; curref <= maxref; curref++) {
            if (argv[curref][0] == '-') {
              curref--;
              break;
            }
          }
          if (curref == (maxref+1)) {
            curref--;
          }
          if (curref > 0) {
            REFNumberArray = (long long *)calloc(curref,sizeof(long long));
            if (REFNumberArray == NULL) {
              printf("ERROR: Failed to allocate REFNumberArray of size %d\n",curref);
              exit(-1);
            }
            maxref = curref;
            for (curref = 0; curref < maxref; curref++) {
              argc--;
              if ((GetREFNumber(*++argv,&REFNumberArray[curref],&refType,1,1) != 0) ||
                  (REFNumberArray[curref] == 0)) {
                printf("ERROR: %s is an unrecognized or illegal catalog reference value\n",*argv);
                errorFlag = 1;
              }
            }
            
          } else {
            printf("ERROR: No star id values follow '-c'\n");
            argc--;
            ++argv;
            errorFlag = 1;
          }


          break;

        case 'g': /* write good points only */
        case 'G': 
          goodFlag =  1;
          AFLAGSMASK = 0;
          BFLAGSMASK = 0;
          if ((argc > 1) && ((*(argv+1))[0] != '-')) {
            argc--;
            ++argv;
            argstr2 = *argv;
            while ((cmdchar2 = *(argstr2++)) != 0) {
              switch(cmdchar2) {
              case 'd':
                AFLAGSMASK |= (1 << FILTER_AFLAG_DEFECT);
                break;
              case 'b':
                AFLAGSMASK |= (1 << FILTER_AFLAG_BLEND) | (1 << FILTER_AFLAG_CASEB);
                break;
              case 'h':
                AFLAGSMASK |= (1 << FILTER_AFLAG_LOW_ALTITUDE);
                break;
              case 'u':
                AFLAGSMASK |= (1 << FILTER_AFLAG_UNCERTAIN_DATE);
                break;
              case 'q':
                AFLAGSMASK |= (1 << FILTER_AFLAG_QUALITY);
                break;
              case 'f':
                AFLAGSMASK |= (1 << FILTER_AFLAG_MULTIPLE_BLEND);
                break;
              case 'l':
                AFLAGSMASK |= (1 << FILTER_AFLAG_HIZOUT);
                break;
              case 'w':
                AFLAGSMASK |= (1 << FILTER_AFLAG_WEDGE);
                break;
              case 'a':
                AFLAGSMASK |= (1 << FILTER_AFLAG_DRAD);
                break;
              case 'z':
                AFLAGSMASK |= (1 << FILTER_AFLAG_DRADBIN);
                break;
              case 'm':
                AFLAGSMASK |= (1 << FILTER_AFLAG_TOO_BRIGHT);
                break;
              case 'v':
                AFLAGSMASK |= (1 << GSC_VARIABLE_BIT);
                break;
              case 's':
                AFLAGSMASK |= (CLASS_MASK << GSC_CLASS_BIT);
                break;
              case '9':
                AFLAGSMASK |= (1 << FILTER_AFLAG_BIN9);
                break;
              case 'g':
                AFLAGSMASK |= (1 << FILTER_AFLAG_LIMITING_MAG);
                break;
              case 'r':
                AFLAGSMASK |= (1 << FILTER_AFLAG_ISO_RMS);
                break;
              case 'e':
                AFLAGSMASK |= (1 << FILTER_AFLAG_LOCAL_RMS);
                break;
              case 'c':
                AFLAGSMASK |= (1 << FILTER_AFLAG_BACKGROUND);
                break;
              case 'y':
                AFLAGSMASK |= (1 << FILTER_AFLAG_CASEC) | (1 << FILTER_AFLAG_CASED);
                break;
              case 't':
                BFLAGSMASK |= (1 << FILTER_BFLAG_PSFSATURATED) ;
                break;
              case 'p':
                BFLAGSMASK |= (1 << FILTER_BFLAG_MAGDEP_MAGCOR) ;
                break;
              default:
                printf("ERROR: Unknown good plot accept flag %c\n",cmdchar2);
                errorFlag = 1;
              }

            }

          }



          break;

        case 'e': /* Exclude bad astrometry series */
        case 'E':
          excludeSeriesCount = NUM_EXCLUDE_SERIES;
          break;

        case 'v': /* verbose */
          verbose += 1;
          break;

        case 'j': /* force error bar factor to 1.0 */
          forceErrorBarFactor = 1;
          break;

        case 'x': /* Do not write files */
          skipWrite = 1;
          break;

        case 'f': /* full dump of photometry data */
          fullFlag = 1;
          break;

        case 'm': /* MySQL dump */
          dumpMySQL= 1;
          break;

        case 'i': /* Ignore limiting magnitudes */
          ignoreBin = 1;
          break;


        case 's': /* show rejection statistics */
          printReject  = 1;
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




        case 'r': /* Use a search radius in arcsec */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strchr(*argv,':')) {
              rad0 = 3600.0 * str2dec(*argv);
            } else {
              rad0 = atof(*argv);
            }
            if (rad0 <= 0.0) {
              printf("ERROR: Illegal radius %s specified\n",*argv);
              errorFlag = 1;
            }
          }

          break;

        case 'V': /* even more verbose */
          verbose += 2;
          break;

        case 'p': /* RA and DEC */
          pFlag = 1;
          if (argc < 3) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            argc--;
            ++argv;
            ra0 = str2ra (*argv);
            argc--;
            ++argv;
            dec0 = str2dec (*argv);
            if ((ra0 < 0.0) || 
                (ra0 >= 360.0) ||
                (dec0 < -90.0) ||
                (dec0 > +90.0)) {
              printf("ERROR: Illegal ra %d or dec %d\n",ra0,dec0);
            }
          }

          break;


        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */
  if ((target_name[0] == 0) && 
      (maxref == 0) &&
      (ra0 < 0.0) &&
      (dec0 < -90.0)) {
    printf("ERROR: No target catalog file was specified\n");
    errorFlag = 1;
  }
  if ((target_name[0] != 0) && (maxref > 0)) {
    printf("ERROR: Only one of -c or -t may be specified\n");
    errorFlag = 1;
  }

  if (lightcurveDirectory[0] == 0) {
    printf("ERROR: No lightcurve directory was specified\n");
    errorFlag = 1;
  } else {
    statResult = stat(lightcurveDirectory,&statbuf);
    if (statResult != 0) {
      printf("ERROR: Lightcurve directory %s does not exist\n",lightcurveDirectory);
      errorFlag = 1;
    }
  }
  /* Open and read in a new catalog */
  if (catalog_name[0] != 0) {
    catalog_handle = fopen(catalog_name,"wt");
    if (catalog_handle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to find the catalog catalog file %s\n",catalog_name);
    } else {
      if (verbose) {
        printf("Found catalog catalog file %s\n",catalog_name);
      }
      fprintf(catalog_handle,"src_name\tra\tdec\tREF\n");
      fprintf(catalog_handle,"--------\t--\t---\t---\n");
    }
  }


  /* Open and read in the target catalog file */
  if (target_name[0] != 0) {
    target_handle = Open(target_name,"r");
    if (target_handle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to find the target catalog file %s\n",target_name);
    } else {
      if (verbose) {
        printf("Found target catalog file %s\n",target_name);
      }
      /* Now read in the target catalog */

      target_header = table_header(target_handle,TABLE_PARSE);
      if (target_header == NULL) {
        printf("ERROR: Failed to read header for %s\n",target_name);
        return(-1);
      }
      target_table = table_loadva(target_handle,
                                  &target_header,
                                  NULL, /* hbase */
                                  NULL, /* rows */
                                  NULL,
                                  sizeof(PHOTTARGET),
                                  &target_nrecs,
                                  TblBuf,"REF"    ,TblOff(PPHOTTARGET,REF),MAX_REF,
                                  TblBuf,"src_name"    ,TblOff(PPHOTTARGET,src_name),MAX_SRC_LENGTH,
                                  TblDbl,"ra",TblOff(PPHOTTARGET,ra),
                                  TblDbl,"dec",TblOff(PPHOTTARGET,dec),
                                  0,"end",0);
      
      if (target_table == NULL) {
        printf("ERROR: Failed to read table for %s\n",target_name);
        errorFlag = 1;
      } else {
        target_alloc = target_nrecs;
        for (target_index = 0; target_index < target_nrecs; target_index++) {
          pTarget = &target_table[target_index];
#if 0
          if (target_index == 0) {
            printf("Index (1) %d REF %s ra %f dec %f\n",target_index,pTarget->REF,pTarget->ra,pTarget->dec);
          }
#endif
          pTarget->haveLocation = 1;
          pTarget->updateflag = 0;
          pTarget->processed = 0;
          pTarget->daschFlag = 0;
          pTarget->gsc_bin_index = 0;
          pTarget->coverage_bin_index = 0;
          pTarget->coverage_table_index = 0;
          pTarget->groupCount = 0;
          pTarget->Stdmag = 0;
          pTarget->color = 0;
          pTarget->class = 0;
          pTarget->VFlag = 0;
          pTarget->MAGFlag = 0;
          pTarget->RaPM = 0;
          pTarget->DecPM = 0;
          pTarget->magcal_magdep = 0;
        }
        if (verbose) {
          printf("read %d records for %s\n",target_nrecs,target_name);
        }
        if (rad0 > 0) {
          if ((dec0 >= -90) || (ra0 >= 0)) {
            printf("ERROR: can not specify -p and -t simultaneously\n");
            errorFlag = 1;
          } else {
            ra0 = target_table[0].ra;
            dec0 = target_table[0].dec;
          }

        }
      }

    }
  }


  if ((rad0 > 0) && 
      ((target_nrecs > 1) ||
       (maxref > 1))) {   
    printf("ERROR: Only one target star may be used as the center of a radius search \n");
    errorFlag = 1;
  }

 

  if ((cFlag+pFlag+tFlag) != 1) {
    printf("ERROR: must specifiy only one of  -c, -p or -t\n");
    errorFlag = 1;
  }
  

 

  if (errorFlag) {
    printf("Usage: find_lightcurves   {-t <target catalog file> | -c <gsc2.3.2 reference ID ... \\\n");
    printf("                                  | -p <J2000 RA hh:mm:ss.sss or dd.dddd> <J2000 DEC dd:mm:ss.sss or dd.dddd>}\n");
    printf("                           -d <lightcurve output directory> \n");
    printf("                          [-s] print rejection statistics\n");
    printf("                          [-b  <blend reference ID>] (Can specify up to ten times)\n");
    printf("                          [-i]  do not plot limiting magnitudes \n");
    printf("                          [-j]  force error_bar_factor to 1.0 \n");
    printf("                          [-v]  verbose\n");
    printf("                          [-f]  include all photometry data in the star database file\n");
    printf("                          [-V  even more verbose]\n");
    printf("                          [-x  skip writing files]\n");
    printf("                          [-m] print MySQL data retrieval commands for debugging\n");
    printf("                          [-n <table>] generate a target catalog file\n");
    printf("                          [-e] excludes series\n");
    printf("                          [-r] search radius in arcsec\n");
    printf("                          [-q <gaia|kepler|apass|atlas|experimental] catalog and filename qualifier ");
    for (excludeSeriesIndex = 0; excludeSeriesIndex < NUM_EXCLUDE_SERIES; excludeSeriesIndex++) {
      printf(" %s",excludeSeriesList[excludeSeriesIndex]);
    }
    
    printf("\n");
    printf("                          [-g  [<includeflags>] write good stars only]\n");
    printf("                               where includeflags adds one or more of the following rejected star types to the plot:\n");
    printf("                               d = defects\n");
    printf("                               b = blends exept multiple exposure blends (Note -b implies this b also)\n");
    printf("                               h = within 23.5 degrees of horizon\n");
    printf("                               l = high smoothing correction required\n");
    printf("                               w = Pickering wedge objects\n");
    printf("                               a = high drad object, spatial bin, or local bin\n");
    printf("                               m = object is too bright\n");
    printf("                               v = flagged as variable in the catalog\n");
    printf("                               s = flagged as non-star in the catalog\n");
    printf("                               9 = star in bin 9\n");
    printf("                               g = too close to the limiting magnitude\n");
    printf("                               r = isophotonic rms is too high\n");
    printf("                               e = local smoothing rms is too high\n");
    printf("                               u = extinction may be in error because the date is uncertain\n");
    printf("                               f = object is a multiple exposure blend\n");
    printf("                               q = plate fails Sumin Tangs quality criteria of 12/11/09 1:07 AM\n");
    printf("                               t = star is saturated\n");
    printf("                               p = star has no magnitude-dependent calibration\n");
    printf("                               c = star has high background\n");
    printf("  Note: until July 21, 2009 a good lightcurve used -g db\n");


    return(-1);
  }

  printf("find_lightcurves of %s %s, ignoreBin %d goodFlag %d AFLAGSMASK %d 0x%08x BFLAGSMASK %d %08x ra %.2f dec %.2f arcsec %.2f excludedSeries %d\n",
         __DATE__,__TIME__,ignoreBin,goodFlag,AFLAGSMASK,AFLAGSMASK,BFLAGSMASK,BFLAGSMASK,ra0,dec0,rad0,excludeSeriesCount);



  /* Connect to the database */


  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    printf("ERROR: DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    printf("ERROR: DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    printf("ERROR: DASCH_PASSWORD is not defined\n");
    return(-1);
  }
                   
  catalogdir = getenv("DASCH_CATALOG");
  if (catalogdir == NULL) {
    fprintf(stderr,"DASCH_CATALOG is not defined\n");
    printf("0\n");
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



  mysql_init(pConnection);


  if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
    }
    return(-1);
  }



  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    printf("DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }

  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    printf("DASCH_PHOT_USERNAME is not defined\n");
    return(-1);
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    printf("DASCH_PHOT_PASSWORD is not defined\n");
    return(-1);
  }
  mysql_init(pPhotConnection);


  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
    }
    return(-1);
  }

  /* Connect to the auxiliary scanner database to find out what plates cover the area we are interested in */

  auxscanner = getenv("DASCH_AUXSCANNER");
  if (password == NULL) {
    printf("ERROR: DASCH_AUXSCANNER is not defined\n");
    return(-1);
  }

  mysql_init(pAuxConnection);


  if (!mysql_real_connect(pAuxConnection,mysqlhost,username,password,auxscanner,0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pAuxConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pAuxConnection),mysql_error(pAuxConnection));
    }
    return(-1);
  }




  InitSeriesTable(pConnection,pPhotConnection);
  gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
  if (gotAnswer != 1)  {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }
  
  if (pPhotGlobal->magnitudeFile == PHOT_MAGNITUDEFILE_YES) {
    InitFileCommon(stdout,pFileCommon,pConnection,pPhotConnection,catalogNumber,1);

#if 0
    printf("ERROR: rad0 trace is enabled\n");
    if (rad0 > 0) {
      pFileCommon->ra0 = ra0;
      pFileCommon->dec0 = dec0;
      pFileCommon->rad0 = rad0;
    }
#endif


    InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);


  } else {
    printf("ERROR: the magnitudes table is no longer supported\n");
    exit(-1);
  }
  if (pPhotGlobal->solutionNumber == PHOT_SOLUTIONNUMBER_YES) {
    solutionNumberFix = 1;
  }



  if (target_name[0] == 0) {
    if (rad0 > 0) {
      if (maxref == 1) {
        /* Need to convert our star into a searchable position */
        refType = GetREFType(REFNumberArray[0]);
        switch(refType) {
        case REF_TYPE_DASCH:
        case REF_TYPE_APASS:
          GetREF(REFNumberArray[0],tmpREF,1,1);
          if (GetDASCHCoordinates(tmpREF,&ra0,&dec0,1,refType) == 0) {
            /* This is a DASCH number */
            break;
          } else {
            printf("ERROR: find_lightcurves line %d failed to decode coordinates from %s\n",__LINE__,tmpREF);
            exit(-1);
          }
          break;

        default:

          curStars = GetStarEntry2(pPhotConnection,pStarTable,tmpREF,1,1); 
          if (curStars == 1) {
            ra0 = pStarTable->ra;
            dec0 = pStarTable->dec;
          } else {
            printf("ERROR: The object %s is not in the photometry database\n",tmpREF);
            exit(-1);
          }
          break;
        }
      }

      /* Perform a cone search */
      if (ConeSearch(&gsc_bin_table,&gsc_bin_alloc,&gsc_bin_count,ra0,dec0,rad0,verbose)) {
        exit(-1);
      }
      coneSearchFlag = 1;
    } else {
      target_nrecs = maxref;
      target_alloc = target_nrecs;
      target_table = (PPHOTTARGET)calloc(target_nrecs,sizeof(PHOTTARGET));
      if (target_table == NULL) {
        printf("STDERR: failed to allocate target_table\n");
        exit(-1);
      }
      for (target_index = 0; target_index < target_nrecs; target_index++) {
        pTarget = &target_table[target_index];
        GetREF(REFNumberArray[target_index],pTarget->REF,1,1);
        strcpy(pTarget->src_name,pTarget->REF);
      }
    }

  }
  if (coneSearchFlag == 0) {
    for (target_index = 0; target_index < target_nrecs; target_index++) {
      pTarget = &target_table[target_index];
#if 0
      if (target_index == 0) {
        printf("Index (2) %d REF %s ra %f dec %f\n",target_index,pTarget->REF,pTarget->ra,pTarget->dec);
      }
#endif
      pTarget->processed = 0;
      pTarget->class = 0;
      pTarget->VFlag = 0;
      pTarget->MAGFlag = 0;
      pTarget->RaPM = 0;
      pTarget->DecPM = 0;
      GetREFNumber(pTarget->REF,&tmpREFNumber,&refType,1,1);
      switch(refType) {
      case REF_TYPE_DASCH:
      case REF_TYPE_APASS:
        if (GetDASCHCoordinates(pTarget->REF,&pTarget->ra,&pTarget->dec,1,refType) == 0) {
          /* This is a DASCH or APASS number */
          pTarget->haveLocation = 1;        
          pTarget->Stdmag = 99.0;
          pTarget->color = 99.0;
          pTarget->gsc_bin_index = GetGSCBin(pGscBin,pTarget->ra,pTarget->dec,&decBin,&raBin,pTarget->REF);
          pTarget->coverage_bin_index = GetGSCBin(pGscBin01,pTarget->ra,pTarget->dec,&decCBin,&raCBin,pTarget->REF);
          pTarget->daschFlag = 1;
          pTarget->updateflag = UPDATEFLAG_NO_UPDATE;
        } else {
          printf("ERROR: find_lightcurves line %d failed to decode coordinates from %s\n",__LINE__,tmpREF);
          exit(-1);
        }
        break;
      default:   

        pTarget->Stdmag = 0;
        pTarget->color = 0;
        pTarget->updateflag = UPDATEFLAG_NEED_UPDATE;
        if (pTarget->haveLocation) {
          pTarget->gsc_bin_index = GetGSCBin(pGscBin,pTarget->ra,pTarget->dec,&decBin,&raBin,pTarget->REF);
          pTarget->coverage_bin_index = GetGSCBin(pGscBin01,pTarget->ra,pTarget->dec,&decCBin,&raCBin,pTarget->REF);
        }
        curStars = GetStarEntry2(pPhotConnection,pStarTable,pTarget->REF,1,1);
        if (curStars == 0) {
          if (dumpMySQL) {
            long long REFNumber;
            GetREFNumber(pTarget->REF,&REFNumber,&refType,1,1);
            printf("SELECT * from stars%s where REFNumber = %lld;\n",catalogString,REFNumber);
          }
          printf("WARNING: No star entry found in the database for %s\n",pTarget->REF);
        } else {
          if ((target_name[0] != 0) && 
              (pTarget->gsc_bin_index != pStarTable->gsc_bin_index)) {
            printf("WARNING: gsc bin indices %d and %d do not agree for star %s\n",
                   pTarget->gsc_bin_index,
                   pStarTable->gsc_bin_index,
                   pTarget->REF);
            pTarget->gsc_bin_index =  pStarTable->gsc_bin_index;
          }
        }

        pTarget->Stdmag = pStarTable->Stdmag;
        pTarget->color  = pStarTable->color;
        pTarget->class = pStarTable->class;
        pTarget->VFlag = pStarTable->VFlag;
        pTarget->MAGFlag = pStarTable->MAGFlag;
        pTarget->RaPM = pStarTable->RaPM;
        pTarget->DecPM = pStarTable->DecPM;
        pTarget->updateflag = pStarTable->updateflag;
        if (target_name[0] == 0) {
          pTarget->ra = pStarTable->ra;
          pTarget->dec = pStarTable->dec;
          pTarget->haveLocation = 1;
          pTarget->gsc_bin_index = GetGSCBin(pGscBin,pTarget->ra,pTarget->dec,&decBin,&raBin,pTarget->REF);
          pTarget->coverage_bin_index = GetGSCBin(pGscBin01,pTarget->ra,pTarget->dec,&decCBin,&raCBin,pTarget->REF);
        }
        if (dumpMySQL) {
          GetMagnitudeFileName(catalogNumber,pTarget->gsc_bin_index,magnitudeFileName,MAX_FILENAME);
          printf("Need magnitude file %s\n",magnitudeFileName);
        }
        if (verbose) {
          printf("Found %s at ra %f dec %f Stdmag %f color %f  modulus %03d gsc_bin_index %d\n",pTarget->REF,pTarget->ra,pTarget->dec,pTarget->Stdmag,pTarget->color,pTarget->gsc_bin_index % MAGNITUDES_MODULUS ,pTarget->gsc_bin_index);
        }
        break;
      }
    } /* Loop on all targets */
  
  }
  if (coneSearchFlag == 0) {
    /* Now check all of the targets and update GSC catalog information if necessary */
    for (target_index = 0; target_index < target_nrecs; target_index++) {
      pTarget = &target_table[target_index];
      if ((pTarget->processed != 0) || 
          (pTarget->updateflag == UPDATEFLAG_NO_UPDATE)) {
        continue;
      }
      if (pTarget->daschFlag != 0) {
        /* DASCH objects will not be in the catalog */
        continue;
      }
      pCurGscImage = SearchStar(pGscBin,pGscCache,indexHandle,catalogHandle,pTarget->gsc_bin_index,&newGscBinIndex,pTarget->REF);
      if (pCurGscImage == NULL) {
        printf("WARNING: no gsc entry found in bin %d for star %s\n",pTarget->gsc_bin_index,pTarget->REF);
      } else {
        pTarget->gsc_bin_index = newGscBinIndex;
        pTarget->coverage_bin_index = GetGSCBin(pGscBin01,pCurGscImage->ra,pCurGscImage->dec,&decCBin,&raCBin,pTarget->REF);
        pTarget->ra = pCurGscImage->ra;
        pTarget->dec = pCurGscImage->dec;
        pTarget->Stdmag = pCurGscImage->Stdmag;
        pTarget->color = pCurGscImage->color;
        pTarget->MAGFlag = pCurGscImage->MAGFlag;
        pTarget->class = pCurGscImage->class;
        pTarget->VFlag = pCurGscImage->VFlag;
        pTarget->RaPM = pCurGscImage ->RaPM;
        pTarget->DecPM = pCurGscImage->DecPM;
        pTarget->updateflag = UPDATEFLAG_NO_UPDATE;
      }

    }


    /* Now sort the target table into increasing gsc bin slots */
    qsort((void*)target_table,target_nrecs,sizeof(PHOTTARGET),TargetCompare);

    /* Count up the gsc bins that we need to examine */
    old_gsc_bin_index = -1;
    gsc_bin_table = (int *)calloc(target_nrecs,sizeof(int));
    if (gsc_bin_table == NULL) {
      printf("ERROR: failed to allocate the gsc bin table \n");
      exit(-1);
    }
    for (target_index = 0; target_index < target_nrecs; target_index++) {
      pTarget = &target_table[target_index];
      if (pTarget->processed) {
        continue;
      }
      if (pTarget->gsc_bin_index != old_gsc_bin_index) {
        gsc_bin_table[gsc_bin_count] = pTarget->gsc_bin_index;
        gsc_bin_count++;
      }
   
      old_gsc_bin_index = pTarget->gsc_bin_index;
    }
  } 
  base_target_nrecs = target_nrecs;
  /*   Find the coverage bins from the gsc bin list */
  coverage_bin_table = (int *)calloc(gsc_bin_count,sizeof(int));
  if (coverage_bin_table == NULL) {
    printf("ERROR: failed to allocate the coverage bin table \n");
    exit(-1);
  }
  coverage_bin_for_gsc_bin_table = (int *)calloc(gsc_bin_count,sizeof(int));
  if (coverage_bin_for_gsc_bin_table == NULL) {
    printf("ERROR: failed to allocate the coverage bin table \n");
    exit(-1);
  }
  old_coverage_bin_index = -1;
  for (gsc_bin_iterator = 0; gsc_bin_iterator < gsc_bin_count; gsc_bin_iterator++) {
    if(GetBinCenter(pGscBin,gsc_bin_table[gsc_bin_iterator],&raCenter,&decCenter,"find_lightcurves")) {
      printf("ERROR: GetBinCenter failed for bin %d\n",gsc_bin_table[gsc_bin_iterator]);
      exit(-1);
    }
    coverage_bin_index = GetGSCBin(pGscBin01,raCenter,decCenter,&decCBin,&raCBin,"find_lightcurves");
#if 0
    printf("gsc_bin_index %d coverage_bin_index %d\n",gsc_bin_table[gsc_bin_iterator],coverage_bin_index);
#endif
    if (coverage_bin_index != old_coverage_bin_index) {
      coverage_bin_table[totalCoverageBins] = coverage_bin_index;
      coverage_bin_for_gsc_bin_table[gsc_bin_iterator] = coverage_bin_index;
      totalCoverageBins++;
    }
    old_coverage_bin_index = coverage_bin_index;
  }
  qsort((void*)coverage_bin_table,totalCoverageBins,sizeof(int),GscBinCompare);
  if (totalCoverageBins > 0) {
    coverage_bin_iterator = 0;
    for (coverage_bin_index = 1; coverage_bin_index < totalCoverageBins; coverage_bin_index++) {
      if (coverage_bin_table[coverage_bin_iterator] != coverage_bin_table[coverage_bin_index]) {
        coverage_bin_iterator++;
        coverage_bin_table[coverage_bin_iterator] = coverage_bin_table[coverage_bin_index];
      }
    }
    totalCoverageBins = coverage_bin_iterator+1;
  }

  printf("Total coverage bins to process: %d first bin: %d\n",totalCoverageBins,old_coverage_bin_index);
  coverage_table = (PCOVERAGE *)calloc(totalCoverageBins,sizeof(PCOVERAGE));
  coverage_table_size = (int *)calloc(totalCoverageBins,sizeof(int));
  if ((coverage_table == NULL)  || (coverage_table_size == NULL)) {
    printf("ERROR: Failed to allocate coverage_table and coverage_table_size\n");
    exit(-1);
  }

  /* Get the plates we are interested in and initialize their data */
  for (coverage_bin_iterator = 0; coverage_bin_iterator < totalCoverageBins; coverage_bin_iterator++) {    
    if (GetCoverage(pAuxConnection,&coverage_table[coverage_bin_iterator],&coverage_table_size[coverage_bin_iterator],coverage_bin_table[coverage_bin_iterator],excludeSeriesCount)) {
      for (coverageIndex = 0; coverageIndex < coverage_table_size[coverage_bin_iterator]; coverageIndex++) {
        pCoverage = &((coverage_table[coverage_bin_iterator])[coverageIndex]);
#if 0
        printf("pCoverage entry %d: %s%05d s%d %8d\n",coverageIndex,pCoverage->series,pCoverage->plateNumber,pPlate->solutionNumber,pCoverage->coverage_bin_index);
#endif
#if 0
        if ((pCoverage->plateNumber == 1146) &&
            (strcmp("rh",pCoverage->series) == 0)) {
          printf("At plate %s%05d coverage_bin_index %d\n",pCoverage->series,pCoverage->plateNumber,pCoverage->coverage_bin_index);
        }
#endif
        for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
          pPlate = &plateArray[plateIndex];
          if ((pPlate->plateNumber == pCoverage->plateNumber) &&
              (pPlate->solutionNumber == pCoverage->solutionNumber) &&
              (strcmp(pPlate->series,pCoverage->series) == 0)) {
              
            break;
              
          }
        }
        if (plateIndex == plateListLength) {
          if ((plateListLength+1) >= plateListAlloc) {
            /* We need a bigger table */
            plateListAlloc += 100;


            tmpPlateArray = realloc(plateArray,plateListAlloc * sizeof(PLATE));
            if (tmpPlateArray == NULL) {
              printf("ERROR: Failed to reallocate plateArray of size %d\n",plateListAlloc * sizeof(PLATE));
              exit(-1);
            }
            plateArray = tmpPlateArray;
          }
          pPlate = &plateArray[plateListLength];
          memset(pPlate,0,sizeof(PLATE));
          pPlate->local_bin_index = -1;
          pPlate->fatal_reject_reason= -1;
          strcpy(pPlate->series,pCoverage->series);
          pPlate->seriesId = GetSeriesId(pPlate->series,1);
          pPlate->plateNumber = pCoverage->plateNumber;
          pPlate->solutionNumber = pCoverage->solutionNumber;
          plateListLength++;


        } else {
          printf("ERROR: duplicate plate entry\n");
        }
        pPlate->coverage_table_refcount++;
        if (plateListLength > 1) {
          pPrevPlate = &plateArray[plateListLength - 2];
          if ((pPlate->plateNumber == pPrevPlate->plateNumber) &&
              (pPlate->seriesId = pPrevPlate->seriesId)) {
            duplicatePlateEntries++;
          }
        }
      }
       

    }
  }
  printf("Reading data for %d exposures from %d plates\n",plateListLength,plateListLength-duplicatePlateEntries);
  
#if 0
  for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
    pPlate = &plateArray[plateIndex];
    printf("plateIndex %d %6s%05d_s%d\n",index,pPlate->series,pPlate->plateNumber,pPlate->solutionNumber);
  }
#endif

  for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
    pPlate = &plateArray[plateIndex];
#if 0
    if ((pCoverage->plateNumber == 1146) &&
        (strcmp("rh",pCoverage->series) == 0)) {
      printf("At plate %s%05d coverage_bin_index %d\n",pCoverage->series,pCoverage->plateNumber,pCoverage->coverage_bin_index);
    }
#endif
    InitPlateEntry(pConnection,pPhotConnection,pPlate,verbose,ignoreBin,pPhotGlobal,dumpMySQL,catalogString,solutionNumberFix);
  }
  /* Sort the plate array in the order of increasing Julian Day */
  qsort((void*)plateArray,plateListLength,sizeof(PLATE),JulianDayCompare);
      
  time(&curTime);
  curTime -= startTime;
  printf("Reading magnitudes at %d seconds\n",curTime);

  /* Now read in all of the magnitude measurements */
  for (coverage_bin_iterator = 0; coverage_bin_iterator < totalCoverageBins; coverage_bin_iterator++) {
    coverage_bin_index = coverage_bin_table[coverage_bin_iterator];
    coverageTableCount = 0;

    /* Go through the plates, select new ones and clean up any that will no longer be referenced */
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (pPlate->coverage_table_select) {
        pPlate->coverage_table_select = 0;
        pPlate->coverage_table_refcount--;
        if (pPlate->coverage_table_refcount < 0) {
          printf("ERROR: coverage_table_refcount %d is negative for plate %s %s%05d\n",pPlate->coverage_table_refcount,pPlate->Plate,pPlate->series,pPlate->plateNumber);
#if 0
          exit(-1);
#else
          pPlate->coverage_table_refcount = 5000000; /* Just keep the plate indefinitely until we understand this error */
#endif
        } else if (pPlate->coverage_table_refcount == 0) {
          /* We are no longer using this plate */
          if (savedBlendCount != 0) {
            printf("ERROR: attempting to clean up with a savedBlendCount of %d\n",savedBlendCount);
            exit(-1);
          }
          CleanupPlate(pPlate);
        }
      }
    }
    /* Now see if which plates are in the new list */
    for (coverageIndex = 0; coverageIndex < coverage_table_size[coverage_bin_iterator]; coverageIndex++) {
      pCoverage = &((coverage_table[coverage_bin_iterator])[coverageIndex]);
      for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
        pPlate = &plateArray[plateIndex];
#if 0
        if ((pCoverage->plateNumber == 1146) &&
            (strcmp("rh",pCoverage->series) == 0)) {
          printf("At plate %s%05d coverage_bin_index %d\n",pCoverage->series,pCoverage->plateNumber,pCoverage->coverage_bin_index);
        }
#endif
        if ((pPlate->plateNumber == pCoverage->plateNumber) &&
            (pPlate->solutionNumber == pCoverage->solutionNumber) &&
            (strcmp(pPlate->series,pCoverage->series) == 0)) {
          pPlate->coverage_table_select = 1;
          coverageTableCount++;
          break;
            
        }
      }
    }
    if (coverageTableCount != coverage_table_size[coverage_bin_iterator]) {
      printf("ERROR: coverage table counts %d %d do not agree for index %d\n",coverageTableCount,coverage_table_size[coverage_bin_iterator],coverage_bin_iterator);
#if 0
      exit(-1);
#endif
    }
    if (verbose) {
      printf("NOTE: coverage bin index %d for Number of plates %d\n",coverage_bin_index,coverageTableCount);
    }
     
    for (gsc_bin_iterator = 0; gsc_bin_iterator < gsc_bin_count; gsc_bin_iterator++) {
      gsc_bin_index = gsc_bin_table[gsc_bin_iterator];
      if(GetBinCenter(pGscBin,gsc_bin_table[gsc_bin_iterator],&raCenter,&decCenter,"find_lightcurves")) {
        printf("ERROR: GetBinCenter failed for bin %d\n",gsc_bin_table[gsc_bin_iterator]);
        exit(-1);
      }
      if (coverage_bin_index != GetGSCBin(pGscBin01,raCenter,decCenter,&decCBin,&raCBin,"find_lightcurves")) {
        continue;
      }
      /* Get the images in and near this bin */
      LocateNoneImages(pGscBin,pFileCommon,gsc_bin_index,&pMagnitudeTable,&magnitudeAlloc,&numMagnitudes,catalogString,1,0);
      if (noneMagnitudeAlloc < numMagnitudes) {
        if (pNoneMagnitudeTable != NULL) {
          free(pNoneMagnitudeTable);
        }
        noneMagnitudeAlloc = numMagnitudes+1000;
        pNoneMagnitudeTable = (PFILESTARIMAGE)calloc(noneMagnitudeAlloc,sizeof(FILESTARIMAGE));
        if (pNoneMagnitudeTable == NULL) {
          printf("ERROR: failed to allocate pNoneMagnitudeTable of size %d\n",noneMagnitudeAlloc);
          exit(-1);
        }
      }

      for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes;magnitudeIndex++) {
        pCurSumstarimage = &pMagnitudeTable[magnitudeIndex];
        pNoneImage = &pNoneMagnitudeTable[magnitudeIndex];
        pFileStarImage = pCurSumstarimage->pFileStarImage;
        memcpy(pNoneImage,pFileStarImage,sizeof(FILESTARIMAGE));

        if (pFileStarImage->gsc_bin_index == gsc_bin_index) {

          if (pFileStarImage->REFNumber != 0) {
            /* See if we need to expand our target array */
            if ((pLastTarget != NULL) && (strcmp(pCurSumstarimage->REF,pLastTarget->REF) == 0)) {
              pCurSumstarimage->Stdmag = pLastTarget->Stdmag;    
            } else {
              /* Search our target array for this entry */
              for (target_index = 0; target_index < target_nrecs; target_index++) {
                pTarget = &target_table[target_index]; 
                if (strcmp(pCurSumstarimage->REF,pTarget->REF) == 0) {
                  pCurSumstarimage->Stdmag = pTarget->Stdmag;
                  pLastTarget = pTarget;
                  break;
                }
              }
              if ((target_index == target_nrecs) && (coneSearchFlag != 0)) {
                /* Since we are doing a cone search, add this entry to the target table */
                if (target_nrecs >= target_alloc) {
                  target_alloc += 100;
                  tmp_target_table = realloc(target_table,target_alloc*sizeof(PHOTTARGET));
                  if (tmp_target_table == NULL) {
                    printf("ERROR: failed to realocate target_table of size %d\n",target_alloc*sizeof(PHOTTARGET));
                  }
                  target_table = tmp_target_table;
                  tmp_target_table = NULL;
                }
                pTarget = &target_table[target_nrecs];
                memset(pTarget,0,sizeof(PHOTTARGET));
                strcpy(pTarget->REF,pCurSumstarimage->REF);
                GetREFNumber(pTarget->REF,&tmpREFNumber,&refType,1,1);
                switch(refType) {
                case REF_TYPE_DASCH:
                case REF_TYPE_APASS:
                  if (GetDASCHCoordinates(pTarget->REF,&pTarget->ra,&pTarget->dec,1,refType) == 0) {
                    /* This is a DASCH number */
                    pTarget->haveLocation = 1;        
                    pTarget->Stdmag = 99.0;
                    pTarget->color = 99.0;
                    pTarget->gsc_bin_index = GetGSCBin(pGscBin,pTarget->ra,pTarget->dec,&decBin,&raBin,pTarget->REF);
                    pTarget->coverage_bin_index = GetGSCBin(pGscBin01,pTarget->ra,pTarget->dec,&decCBin,&raCBin,pTarget->REF);
                    pTarget->daschFlag = 1;
                    pTarget->updateflag = UPDATEFLAG_NO_UPDATE;
                    target_nrecs++;
                    break;
                  } else {
                    printf("ERROR: find_lightcurves line %d failed to decode coordinates from %s\n",__LINE__,tmpREF);
                    exit(-1);
                  }
                  break;

                default:

              
                  curStars = GetStarEntry2(pPhotConnection,pStarTable,pTarget->REF,1,1);
                  if (curStars <= 0) {
                    printf("ERROR: Can not find star table entry for %s\n",pTarget->REF);
                    pTarget->processed = 1;
                  } else {
                    radactual = wcsdist(ra0,dec0,pStarTable->ra,pStarTable->dec) * 3600;
                    if ((rad0 < 0) || (radactual > rad0)) {
                      /* Not interested, but use it as a placeholder */
                      pTarget->processed = 1;
                    } else {
                      pTarget->Stdmag = pStarTable->Stdmag;
                      pTarget->color  = pStarTable->color;
                      pTarget->class = pStarTable->class;
                      pTarget->VFlag = pStarTable->VFlag;
                      pTarget->MAGFlag = pStarTable->MAGFlag;
                      pTarget->RaPM = pStarTable->RaPM;
                      pTarget->DecPM = pStarTable->DecPM;
                      pTarget->updateflag = pStarTable->updateflag;
                      pTarget->ra = pStarTable->ra;
                      pTarget->dec = pStarTable->dec;
                      pTarget->haveLocation = 1;
                      pTarget->gsc_bin_index = gsc_bin_index;
                      pTarget->coverage_bin_index = coverage_bin_index;
                      strcpy(pTarget->REF,pStarTable->REF);
                      strcpy(pTarget->src_name,pStarTable->REF);
                      sprintf(tmp_src_name,"_%06.0f",radactual);
                      strcat(pTarget->src_name,tmp_src_name);
                    }
                  }
          
                  target_nrecs++;
                }
              }
            }
          }
        }
        pCurSumstarimage->AFLAGSCOPY = pFileStarImage->AFLAGS;
        
        if (fullFlag) {

          /* Note: color, dra, ddec, and THRESHOLD do not exist in _filestarimage */
          pCurSumstarimage->colorterm = 99.0;
          pCurSumstarimage->errorcolor = 99.0;
          pCurSumstarimage->colorflag = -1;
        }

      }
    
      if (skipWrite != 0) {
        numMagnitudes = 0;
      }
      /* Now assign DASCH numbers to these objects */

      ProcessNoneImagesX(pGscBin,pFileCommon,&target_table,&target_nrecs,&target_alloc,pNoneMagnitudeTable,numMagnitudes,verbose,histogramTable,minimumImageCount,&groupNumber,NULL,ra0,dec0,rad0,NULL,0,0,gsc_bin_index);
      for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes;magnitudeIndex++) {
        pCurSumstarimage = &pMagnitudeTable[magnitudeIndex];
        pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
        if (strcmp(pCurSumstarimage->REF,"NONE") == 0) {
          GetREF(pFileStarImage->REFNumber,pCurSumstarimage->REF,1,1);
        }

      }
      /* Check to see if we have a better position for any new unmatched objects */
      if ((coneSearchFlag == 0) && (target_nrecs > base_target_nrecs)) {
        for (target_index2 = base_target_nrecs; target_index2 < target_nrecs; target_index2++) {
          pTarget2 = &target_table[target_index2];
          if (strstr(pTarget2->REF,"DASCH") == NULL) {
            continue;
          }
          for (target_index3 = 0; target_index3 < base_target_nrecs; target_index3++) {
            pTarget3 = &target_table[target_index3];
            if (strcmp(pTarget2->REF,pTarget3->REF) == 0) {
              /* We have a match.  Copy over the ra and dec */
              pTarget3->ra = pTarget2->ra;
              pTarget3->dec = pTarget2->dec;
            }
          }
        }
        target_nrecs = base_target_nrecs; /* Done with additional table entries */
      }
      master_nrecs += numMagnitudes;
      for (target_index2 = 0; target_index2 < target_nrecs; target_index2++) {
        pTarget2 = &target_table[target_index2];
        if (pTarget2->processed) {
          continue;
        }
#if 0

        if (strcmp(pTarget2->REF,"DASCH_J192901.4+431030") == 0) {
          printf("At %s in line %d\n",pTarget2->REF,__LINE__);
          printFlag = 1;
        } else {
          printFlag = 0;
        }
#endif

        if (pTarget2->gsc_bin_index == gsc_bin_index) {
          /* This star is in our magnitudes table.  Process it. */
          if (verbose) {
            printf("Processing %s %s\n",pTarget2->REF,pTarget2->src_name);
          }
          match_count = 0;
          for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
            pMaster = &pMagnitudeTable[magnitudeIndex];
            pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
            excludeSeriesFlag = 0;
            for (excludeSeriesIndex = 0; excludeSeriesIndex < excludeSeriesCount; excludeSeriesIndex++) {

              if (strcmp(pMaster->series,excludeSeriesList[excludeSeriesIndex]) == 0) {
                excludeSeriesFlag = 1;
                imageExcludeCount++;
                break;
              }
            }
            if (excludeSeriesFlag) {
              continue;
            }
	  
#if 0
            if (strstr(pMaster->REF,"DASCH") != 0) {
              printf("pMaster->REF is %s at line %d\n",pMaster->REF,__LINE__);
            }
            if (printFlag) {
              printf("pMaster->REF is %s at line %d\n",pMaster->REF,__LINE__);
            }
#endif

            if (strcmp(pMaster->REF,pTarget2->REF) == 0) {
              for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
                pPlate = &plateArray[plateIndex];
#if 0
                if ((pCoverage->plateNumber == 1146) &&
                    (strcmp("rh",pCoverage->series) == 0)) {
                  printf("At plate %s%05d coverage_bin_index %d\n",pCoverage->series,pCoverage->plateNumber,pCoverage->coverage_bin_index);
                }
#endif
                if (pPlate->coverage_table_select == 0) {
                  continue;
                }
                if ((pPlate->plateNumber == pFileStarImage->plateNumber) &&
                    (pPlate->solutionNumber == pFileStarImage->solutionNumber) &&
                    (strcmp(pPlate->series,pMaster->series) == 0)) {
                  break;
                }
              }
              if (plateIndex == plateListLength) {
                if (coverageErrorCount == 0) {
                  printf("ERROR: can not find plate %s%05d s%d from coverage table\n",
                         pMaster->series,pFileStarImage->plateNumber,pPlate->solutionNumber);
                }
                coverageErrorCount++;
                continue;
              }

              if (strlen(pPlate->Plate) > 0) {
                strcpy(pMaster->Plate,pPlate->Plate);
              }
              pMaster->Stdmag = pTarget2->Stdmag;
              pMaster->AFLAGSCOPY = pFileStarImage->AFLAGS;
              if (pPlate->local_bin_table == NULL) {
                /* Need to allocate a local bin table now */
                pPlate->local_bin_table = (PLOCALBIN)calloc(X_DMAGBINS_NORMAL*Y_DMAGBINS_NORMAL,sizeof(LOCALBIN));
                if (pPlate->local_bin_table == NULL) {
                  printf("ERROR: Failed to allocate the local bin table \n");
                  exit(-1);
                }
                for (local_bin_index = 0; local_bin_index < (X_DMAGBINS_NORMAL*Y_DMAGBINS_NORMAL); local_bin_index++) {
                  pLocalBin = &pPlate->local_bin_table[local_bin_index];
                  pLocalBin->local_bin_index = -2;
                }
              }
            
              if (pFileStarImage->npoints_local < 0) {
                pLocalBin = &pPlate->local_bin_table[pFileStarImage->local_bin_index];
                if (pLocalBin->local_bin_index == pFileStarImage->local_bin_index) {
                  pFileStarImage->dradRMS2 = pLocalBin->dradRMS2;
                  pFileStarImage->npoints_local = pLocalBin->npoints_local;
                  pFileStarImage->rejectFlag = pLocalBin->rejectFlag;
                  pFileStarImage->extinction = pLocalBin->extinction;
                  pFileStarImage->magcor_local = pLocalBin->magcor_local;
                  pFileStarImage->magcal_local_error = pLocalBin->magcal_local_error;
                } else if (pLocalBin->local_bin_index == -2) {

                  /* Here we need to read a bin of in the local bin table */
                  if (GetPhotLocalBin(pPhotConnection,pLocalBin,1,pPlate->series,pPlate->plateNumber,pPlate->solutionNumber,pFileStarImage->local_bin_index,catalogString,solutionNumberFix,0)) {
                    printf("WARNING: Failed to read local bin table for %s%05d  seriesId = %d plateNumber = %d solutionNumber = %d\n",pPlate->series,pPlate->plateNumber,pPlate->seriesId,pPlate->plateNumber,pPlate->solutionNumber);
                    pLocalBin->local_bin_index = -1;
                    if (dumpMySQL) {
                      printf("SELECT * from localbin where seriesId = %d and plateNumber = %d and exposureNumber = %d and local_bin_index = %d;\n",pPlate->seriesId,pPlate->plateNumber,pPlate->solutionNumber,pFileStarImage->local_bin_index);

                    } else {
                      if (pPlate->fatal_reject_reason < 0) {
                        pPlate->fatal_reject_reason = (FATAL_REASON_NOLOCALBINENTRY);
                      }    
                      pPlate->errorFlag = 1;

                    }
                  } else {
                    localbinCount++;
                    pFileStarImage->dradRMS2 = pLocalBin->dradRMS2;
                    pFileStarImage->npoints_local = pLocalBin->npoints_local;
                    pFileStarImage->rejectFlag = pLocalBin->rejectFlag;
                    pFileStarImage->extinction = pLocalBin->extinction;
                    pFileStarImage->magcor_local = pLocalBin->magcor_local;
                    pFileStarImage->magcal_local_error = pLocalBin->magcal_local_error;

                  }
                } else {
                  if (pPlate->fatal_reject_reason < 0) {
                    pPlate->fatal_reject_reason = (FATAL_REASON_NOLOCALBINENTRY);
                  }    
                  pPlate->errorFlag = 1;
                }
              } else {
                /* Save this data in the local bin table for the next read */
                pLocalBin = &pPlate->local_bin_table[pFileStarImage->local_bin_index];
                pLocalBin->local_bin_index = pFileStarImage->local_bin_index;
                pLocalBin->dradRMS2 = pFileStarImage->dradRMS2;
                pLocalBin->npoints_local = pFileStarImage->npoints_local;
                pLocalBin->rejectFlag = pFileStarImage->rejectFlag;
                pLocalBin->extinction = pFileStarImage->extinction;
                pLocalBin->magcor_local = pFileStarImage->magcor_local;
                pLocalBin->magcal_local_error = pFileStarImage->magcal_local_error;
              }
#if 0

              if (strcmp(pCurSumstarimage->REF,"DASCH_J192901.4+431030") == 0) {
                printf("At %s in line %d\n",pTarget2->REF,__LINE__);
              }
#endif

              if ((writePointResult = WritePoint(plateArray,plateListLength,pMaster,&debugFlag,&staleEntryCount,fullFlag)) < 0) {
                return(-1);
              }
              match_count+= writePointResult;
            }
          }
          if (vectorAlloc < plateListAlloc) {
            if (vector1 != NULL) {
              free(vector1);
            }
            if (vector2 != NULL) {
              free(vector2);
            }
            vector1 = (double*)calloc(plateListAlloc,sizeof(double));
            vector2 = (double*)calloc(plateListAlloc,sizeof(double));
            if ((vector1 == NULL) || 
                (vector2 == NULL)) {
              printf("ERROR: failed to allocate vectors of length %d\n",plateListLength);
              exit(-1);
            }
            vectorAlloc = plateListAlloc;
          }
          if (strstr(pTarget2->src_name,"PLACEHOLDER") == NULL) {
            if (skipWrite) {
              outputResult = 0;
            } else { 
#if 0
              if (strcmp(pTarget2->REF,"DASCH_J192901.4+431030") == 0) {
                printf("At %s in line %d\n",pTarget2->REF,__LINE__);
              }
#endif
              outputResult = WriteLightcurve(pConnection,pPhotConnection,plateArray,plateListLength,match_count,pTarget2,lightcurveDirectory,vector1,vector2,verbose,goodFlag,printReject,&clipmed,AFLAGSMASK,BFLAGSMASK,ignoreBin,allLightcurves,dumpMySQL,&spatialBinCount,&wcsCount,&localbinCount,catalogString,solutionNumberFix,fullFlag,forceErrorBarFactor);
              SaveLightcurve(plateArray,plateListLength,pTarget2,blendStrings,blendCount,blendArray,&savedBlendCount,clipmed);
            }
          }
          pTarget2->processed = 1;
          /* Re-initialize the plateArray for the next star */
          for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
            pPlate = &plateArray[plateIndex];
            ReInitializePlate(pPlate);
          }
          if (verbose > 1) {
            debugFlag = 1;
          }
          if (outputResult < 0) {
            return(-1);
          } else {
            outputPlotCount += outputResult;
            match_count = 0;
          }

          /* If we need to produce a blend graph, do it now */
          if (blendCount > 0) {
            if (blendCount == savedBlendCount) {
              match_count = 0;
              MergeBlendArray(blendStrings,blendCount,blendArray,&savedBlendCount,&blendTarget,plateListLength);
              if (skipWrite) {
                outputResult = 0;
              } else {
                outputResult = WriteLightcurve(pConnection,pPhotConnection,blendArray[0],plateListLength,match_count,&blendTarget,lightcurveDirectory,vector1,vector2,verbose,goodFlag,printReject,&clipmed,AFLAGSMASK,BFLAGSMASK,ignoreBin,allLightcurves,dumpMySQL,&spatialBinCount,&wcsCount,&localbinCount,catalogString,solutionNumberFix,fullFlag,forceErrorBarFactor);
              }
              if (outputResult < 0) {
                return(-1);
              } else {
                outputPlotCount += outputResult;
                match_count = 0;
              }
              savedBlendCount = 0;
            } 
          }

          for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
            pPlate = &plateArray[plateIndex];
            ReInitializePlate(pPlate);
          }
      
        }
      }
    } /* gsc bin iterator */
  } /* coverage bin iterator */
    /* Write out a region file */
  strcpy(regionFileName,lightcurveDirectory);
  strcat(regionFileName,"/stars.reg");
  regionFileHandle = fopen(regionFileName,"wt");
  if (regionFileHandle == NULL) {
    printf("ERROR: Failed to open %s\n",regionFileName);
    exit(-1);
  }

  for (target_index = 0; target_index < target_nrecs; target_index++) {
    char typeString[40];
    char magnitudeString[40];
    pTarget = &target_table[target_index];
#if 0
    if ((pTarget->ra == 0.0) && (pTarget->dec == 0)) {
      printf("ERROR: zero ra and dec for haveLocation %d\n",pTarget->haveLocation);
    }
#endif
    if (pTarget->haveLocation == 0) {
      continue;
    }

    if (strstr(pTarget->src_name,"PLACEHOLDER") == NULL) {
      if (strstr(pTarget->src_name,"DASCH")) {
        strcpy(typeString,"x");
        magnitudeString[0] = 0;
      } else {
        strcpy(typeString,"circle");
        sprintf(magnitudeString," text = {%.1f}",pTarget->Stdmag);
      }
      fprintf(regionFileHandle,"j2000;point(%f,%f) # point = %s%s\n",pTarget->ra,pTarget->dec,typeString,magnitudeString);
    }
  }
  fclose(regionFileHandle);
  if (catalog_handle != NULL) {
    for (target_index = 0; target_index < target_nrecs; target_index++) {
      pTarget = &target_table[target_index];
      if (pTarget->haveLocation != 0) {
        fprintf(catalog_handle,"%s\t%f\t%f\t%s\n",pTarget->src_name,pTarget->ra,pTarget->dec,pTarget->REF);
      }
    }
  }
  if (printReject) {
    char goodBinsFileName[MAX_BUFFER];
    FILE *goodBinsFileHandle = NULL;
    char *seriesTable = (char *)calloc(MAX_SERIES_STRING*(MAX_SERIES+1),1);
    int numSeries = 0;
    int seriesIndex;
    int *spatialBinTable = NULL;
    char series[MAX_SERIES_STRING];
    int plateNumber;
    int mosaicNumber;
    int binning;
    int rotation;
    int spatial_bin;
    /* Print an overall summary of successful bins as a function of plate series */
    if (seriesTable == NULL) {
      printf("ERROR: seriesTable allocation failure\n");
      exit(-1);
    }
    /* First figure out how many series we have */
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (ParseFilename2(pPlate->Plate,series,&plateNumber) != 0) {
        for (seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
          if (strcmp(series,&seriesTable[seriesIndex*MAX_SERIES_STRING]) == 0) {
            break;
          }
        }
        if (seriesIndex >= numSeries) {
          strcpy(&seriesTable[numSeries*MAX_SERIES_STRING],series);
          numSeries++;
          if (numSeries >= MAX_SERIES) {
            printf("ERROR: MAX_SERIES is too small\n");
            exit(-1);
          }
        }
      }
    }
    /* O.K. To allocate the bin table */
    strcpy(&seriesTable[numSeries*MAX_SERIES_STRING],"TOTAL");
    spatialBinTable = calloc((numSeries+1)*(MAX_SPATIAL_BINS+1),sizeof(int));
    if (spatialBinTable == NULL) {
      printf("ERROR: failed to allocate spatialBinTable\n");
      exit(-1);
    }
    /* Now count up our successful bins */
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (ParseFilename2(pPlate->Plate,series,&plateNumber) != 0) {
        for (seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
          if (strcmp(series,&seriesTable[seriesIndex*MAX_SERIES_STRING]) == 0) {
            break;
          }
        }
        if (seriesIndex < numSeries) {
          spatialBinTable[seriesIndex*(MAX_SPATIAL_BINS+1)]++; /* Entry total is the total count of plates in the series */
          spatialBinTable[numSeries*(MAX_SPATIAL_BINS+1)]++; /* Entry total is the total count of plates in the series */
          if (pPlate->spatial_bin_table != NULL) {
            for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
              if (spatial_bin == pPlate->spatial_bin_table[spatial_bin].spatial_bin) {
                spatialBinTable[(seriesIndex*(MAX_SPATIAL_BINS+1))+spatial_bin]++;
                spatialBinTable[(numSeries*(MAX_SPATIAL_BINS+1))+spatial_bin]++;
              }
            }
            
          }
        }
      }
    }
    

    strcpy(goodBinsFileName,lightcurveDirectory);
    strcat(goodBinsFileName,"/binsummary.txt");
    goodBinsFileHandle = fopen(goodBinsFileName,"wt");
    if (goodBinsFileHandle == NULL) {
      printf("ERROR: Failed to open %s\n",goodBinsFileName);
      exit(-1);
    }


    fprintf(goodBinsFileHandle,"  LEGEND: \n");
    fprintf(goodBinsFileHandle,"    Solid Circle (17)                for a good point \n");
    fprintf(goodBinsFileHandle,"    Small solid circle (17.5)        for a good point, bin 9 \n");
    fprintf(goodBinsFileHandle,"    Open Circle (22)                 for a large error point \n");
    fprintf(goodBinsFileHandle,"    Small open circle (22.5)         for a large error point, bin 9 \n");
    fprintf(goodBinsFileHandle,"    Arrow pointing down (31)         for not found \n");
    fprintf(goodBinsFileHandle,"    Small arrow pointing down (31.5) for not found, bin 9 \n");
    fprintf(goodBinsFileHandle,"    Asterisk (3)                     for a plate defect \n");
    fprintf(goodBinsFileHandle,"    Triangle (7)                     for a high drad object \n");
    fprintf(goodBinsFileHandle,"    Open Square (6)                  for a Pickering Wedge object \n");
    fprintf(goodBinsFileHandle,"    Solid Square (16)                for a high extinction \n");
    fprintf(goodBinsFileHandle,"    Arrow pointing up (30)           for a too bright object \n");
    fprintf(goodBinsFileHandle,"    Diamond (11.5)                   for a blended object \n");
    fprintf(goodBinsFileHandle,"    Cusped Square (10.5)             for a high local correction object \n");
    fprintf(goodBinsFileHandle,"    Letter 'b' (98.5)                for a good object with sextractor blend or neighbors flags \n");
    fprintf(goodBinsFileHandle,"    Letter 'c' (99.5)                for a good object without color correction \n");
    fprintf(goodBinsFileHandle,"    Letter 's' (115.5)               for a good object which is saturated \n");
    fprintf(goodBinsFileHandle,"    Letter 'q' (113.5)               for a good object in a plate with:\n");
    fprintf(goodBinsFileHandle,"                                         multiple exposures\n");
    fprintf(goodBinsFileHandle,"                                         gratings\n");
    fprintf(goodBinsFileHandle,"                                         red or yellow filters\n");
    fprintf(goodBinsFileHandle,"                                         Pickering wedge\n");
    fprintf(goodBinsFileHandle,"                                         Spectra\n");
    fprintf(goodBinsFileHandle,"                                         Fails colorterm metropolis algorithm\n");
    fprintf(goodBinsFileHandle,"                                         Fails colorterm limits %f to %f for GSC colors\n",MIN_BLUE_COLORTERM,MAX_BLUE_COLORTERM);
    fprintf(goodBinsFileHandle,"                                         Fails colorterm limits %f to %f for KIC\n",MIN_KEPLER_BLUE_COLORTERM,MAX_KEPLER_BLUE_COLORTERM);
    fprintf(goodBinsFileHandle,"                                         Fails colorterm limits %f to %f for APASS\n",MIN_APASS_BLUE_COLORTERM,MAX_APASS_BLUE_COLORTERM);
    fprintf(goodBinsFileHandle,"                                         Fails colorterm limits %f to %f for GAIA\n",MIN_GAIA_BLUE_COLORTERM,MAX_GAIA_BLUE_COLORTERM);
    fprintf(goodBinsFileHandle,"\n");
    fprintf(goodBinsFileHandle,"  Series Plates  bin1  bin2  bin3  bin4  bin5  bin6  bin7  bin8  bin9\n");
    for (seriesIndex = 0; seriesIndex <= numSeries; seriesIndex++) {
      fprintf(goodBinsFileHandle,"%8s ",&seriesTable[seriesIndex*MAX_SERIES_STRING]);
      for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
        fprintf(goodBinsFileHandle," %5d",spatialBinTable[(seriesIndex*(MAX_SPATIAL_BINS+1))+spatial_bin]);
      }
      fprintf(goodBinsFileHandle,"\n");
    }
    fclose(goodBinsFileHandle);
    
    free(spatialBinTable);
    free(seriesTable);

  }
  if ((blendCount > 0) && (savedBlendCount != 0)) {
    printf("ERROR: blendCount %d does not equal found blend objects %d\n",blendCount,savedBlendCount);
  }


  if (verbose) {
    printf("read %d records from the magnitude table\n",master_nrecs);
  }


  /* All done.  Clean up */
 

  if (target_table != NULL) {
    free(target_table);
  }
  if (target_header != NULL) {
    table_hdrfree(target_header);
  }
  if (target_handle != NULL) {
    Close(target_handle);
  }


  if(plateArray != NULL) {
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (pPlate->coverage_table_select) {
        pPlate->coverage_table_refcount--;
        if (pPlate->coverage_table_refcount < 0) {
          printf("ERROR: coverage_table_refcount is negative %d for plate %s %s%05d\n",pPlate->coverage_table_refcount,pPlate->Plate,pPlate->series,pPlate->plateNumber);
        } else if (pPlate->coverage_table_refcount == 0) {
          /* We are no longer using this plate */
          pPlate->coverage_table_select = 0;
          CleanupPlate(pPlate);
        } else {
          printf("ERROR: coverage_table_refcount is %d for plate %s %s%05d\n",pPlate->coverage_table_refcount,pPlate->Plate,pPlate->series,pPlate->plateNumber);
        }
      }
      /* Do it anyways to be safe */
      CleanupPlate(pPlate);

    }
  
    free(plateArray);
  }
  if (vector1 != NULL) {
    free(vector1);
  }
  if (vector2 != NULL) {
    free(vector2);
  }
  if (rmsHandle != NULL) {
    fclose(rmsHandle);
  }

  if (pMagnitudeTable != NULL) {
    free(pMagnitudeTable);
  }
  if (pNoneMagnitudeTable != NULL) {
    free(pNoneMagnitudeTable);
  }


  for (blendIndex = 0; blendIndex < savedBlendCount; blendIndex++) {
    free(blendArray[blendIndex]);
    blendArray[blendIndex] = NULL;
 
  }

  for (coverageBinIndex = 0; coverageBinIndex < totalCoverageBins; coverageBinIndex++) {
    if (coverage_table[coverageBinIndex] != NULL) {
      free(coverage_table[coverageBinIndex]);
      coverage_table[coverageBinIndex] = NULL;
    }

  }
  if (coverage_table != NULL) {
    free(coverage_table);
  }
  if (coverage_table_size != NULL) {
    free(coverage_table_size);
  }
  if (REFNumberArray != NULL) {
    free(REFNumberArray);
  }
  if (gsc_bin_table != NULL) {
    free(gsc_bin_table);
  }
  if (coverage_bin_table != NULL) {
    free(coverage_bin_table);
  }
  if (coverage_bin_for_gsc_bin_table != NULL) {
    free(coverage_bin_for_gsc_bin_table);
  }
  if (catalog_handle != NULL) {
    fclose(catalog_handle);
  }


  time(&curTime);
  curTime -= startTime;

  FreeFileCommon(pFileCommon,1);

  mysql_close(pPhotConnection);
  mysql_close(pConnection);
  mysql_close(pAuxConnection);


  fprintf(stdout,"Output stars %d from %d plates seconds %d MySQL queries %d spatialBinCount %d wcsCount %d localbinCount %d,staleEntryCount %d imagesExcluded %d\n",
          outputPlotCount,plateListLength,curTime,GetQueryCount(),spatialBinCount,wcsCount,localbinCount,staleEntryCount,imageExcludeCount);
    
  if (coverageErrorCount > 0) {
    fprintf(stdout,"ERROR: The coverage table in the MySQL database is out of date, %d point(s) were not plotted\n",coverageErrorCount);
  }

  return(0);
}

