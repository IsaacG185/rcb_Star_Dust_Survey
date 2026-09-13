// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#ifndef _PHOTOMETRY_H
#define _PHOTOMETRY_H 1

#include <mysql.h>

#include "galaxyutils.h"
#include <semaphore.h>

// photdb.c:
//
// This function will connect to a photometry SQL database, aborting the process
// if anything goes wrong. The choice of database is set by the file
// `$DASCH_PHOT_ROOT/dbname.txt`, and the database connection parameters come
// from $DASCH_MYSQLHOST, $DASCH_USERNAME, and $DASCH_PASSWORD. The function
// *sets* the environment variables $DASCH_PHOT_MAGNITUDES@KEY@ to
// subdirectories of $DASCH_PHOT_ROOT, e.g. $DASCH_PHOT_MAGNITUDES2 =
// $DASCH_PHOT_ROOT/apass.
//
// If the first character of the `dbname.txt` file is an asterisk, that
// indicates that the photometry "mag" files use the new directory layout
// scheme. The asterisk is ignored for the database name.
extern void dasch_init_photdb(MYSQL *conn);
extern int dasch_photdb_magfiles_use_new_layout(void);

/* These are the major partitions of the database */
#define MAX_NONE_HISTOGRAM 3000
#define SEARCH_NONE_DRAD_FACTOR 2.5
#define SEARCH_NONE_MIN_DRAD 0.5 /* Pixels */

/* The following minumum arcsecond parameter combines distinct images on high resolution plates
 * into a single blend for comparison with images on patrol plates to avoid splitting transient
 * candidates into multiple groups
 */
#define SEARCH_NONE_MIN_ARCSEC  0.0
#define SEARCH_NONE_MAX_ARCSEC 18.0  /* A test on March 13 to reduce this did not improve KV UMa */
#define PMSTATS_MAX_RMS_ERROR  0.0007  /* error in degrees 2.52 arcsec */

#define MAX_VERSION_NAME 24
#define MAX_TIMESTAMP_STRING MAX_VERSION_NAME /* These two must be equal */
#define MAX_OVERLAP_TABLE 20
#define PHOT_KEEPALIVE_YES 1
#define PHOT_KEEPALIVE_NO  2

#define MAX_GROUP_DRAD_HISTOGRAM 100 /* group histogram size in arcsec */
#define BIN_HISTOGRAM_SIZE 1800

/*
 * The following permissive flag affects the id_ tables for matched and unmatched,
 * and the search_none tables and summary tables.
 *
 * It includes all plotted objects including:
 *  AFLAGS: blends, defects,
 *  BFLAGS: no magnitude-dependent correction, sextractor blends,
 *  all quality flags: multiple,grating, color, colorterm, wedge, spectra plates.
 *
 * Does not affect the plate_outlier_*.db files
 */

#define PERMISSIVE_ID_TABLE 1

#define PHOT_MAGNITUDEFILE_YES 1  /* A magnitude file is in use */
#define PHOT_MAGNITUDEFILE_NO  2

#define PHOT_SOLUTIONNUMBER_YES 1 /* The exposureNumber column has been renamed to solutionNumber */
#define PHOT_SOLUTIONNUMBER_NO  2

#define COMBINE_ALGORITHM_BOTH 1      /* Display both images */
#define COMBINE_ALGORITHM_BEST 2      /* Display the best image */
#define COMBINE_ALGORITHM_AVERAGE 3   /* Display the average of the images */
#define COMBINE_ALGORITHM_DEFAULT COMBINE_ALGORITHM_AVERAGE

#define PASSBIT_GSC232 1
#define MAGNITUDES_MODULUS 200
#define MAX_PLATE_NAME 32

#define MAX_GRID_ENTRY 4000 /* The maximum size of the calibration table for a spatial bin */

/* Note database file semaphores may be found in /dev/shm/sem.dasch_phot_magnitudes_%d where %d is the catalog number */
#define SEM_NAME_STRING "/dasch_phot_magnitudes_%d"
#define SEM_TIMEOUT 600 /* 10 minutes */

#define KEPLERFIELD_YES 1 /* Object is on a Kepler CCD chip */
#define KEPLERFIELD_NO  2

#define NULL_PROPER_MOTION 999999.0 /* Invalid proper motion mas/yr */
#define MAX_PROPER_MOTION  990000.0 /* Maximum proper motion mas/yr */
#define ENABLE_REMATCH 1  /* set enableRematch to 1 in affected modules */

typedef struct _photglobal {
  double epoch;
  double equinox;
  int maxSpatialBins;
  int keepalive;
  int magnitudeFile;
  int solutionNumber;
  int currentVersion;
  int minVersionId[MAX_CATALOG_NUMBER];
  char timeStamp[MAX_TIMESTAMP_STRING+1];
  char versionName[MAX_VERSION_NAME+1];
  char versionDate[MAX_TIMESTAMP_STRING+1];
} PHOTGLOBAL,*PPHOTGLOBAL;

typedef struct _photplates {
  char series[MAX_SERIES_STRING]; /* Series */
  int seriesId; /* Series ID */
  int plateNumber;
  int mosaicNumber;
  int versionId;
  double THRESHOLD;
  double scale;
  int mosaicWidth;
  int mosaicHeight;
  int leftMargin;
  int rightMargin;
  int bottomMargin;
  int topMargin;
  int quality;
  int unmatched; /* Number of good unmatched images */
  int exposures;  /* Exposures fitted for this plate */
} PHOTPLATES, *PPHOTPLATES;

typedef struct _phot_spatial_bin {
  char series[MAX_SERIES_STRING]; /* Series */
  int seriesId;                   /* Series ID */
  int plateNumber;
  int solutionNumber;
  int versionId;
  /* The following come from the /ingest/<plate>.out.spatial_bins.db file */
  int spatial_bin;
  double rms;
  int n1;
  int n2;
  double limiting_mag;
  double limiting_iso;
  double max_bright_mag;
  double max_bright_iso;
  double upper_limit;
  double colorterm;
  double errorcolor;
  int colorflag;

} PHOT_SPATIAL_BIN,*PPHOT_SPATIAL_BIN;

/* The update flag in the starentry structure lets update_summary know which stars need
   to have their DASCH magnitudes recomputed */
#define UPDATEFLAG_IGNORE      0
#define UPDATEFLAG_NEED_UPDATE 1
#define UPDATEFLAG_NO_UPDATE   2

typedef struct _starentry {
  long long REFNumber;
  char REF[MAX_REF];          /* GSC2.3.2 reference number */
  int versionId;
  double Stdmag;
  double color;
  double dec;
  double ra;
  int MAGFlag;
  int class;
  int VFlag;
  double RaPM;
  double DecPM;
  int gsc_bin_index;
  int updateflag;
  /* The following entries are not in the database */
  int processedflag;         /* If set, this star has been processed */
} STARENTRY,*PSTARENTRY;

typedef struct _localbin
{
  char series[MAX_SERIES_STRING]; /* Series */
  int seriesId;
  int plateNumber;
  int solutionNumber;
  int local_bin_index;
  int versionId;
  double altitude;
  double extinction;
  double magcor_local;       /* zout */
  double magcal_local_error; /* errout */
  int npoints_local;
  int rejectFlag;
  int drad_bin_count;
  int drad_bin_size;
  int drad_reject_count;
  double draMedian;
  double draRMS;
  double ddecMedian;
  double ddecRMS;
  int drad_bin_count2;
  int drad_bin_size2;
  int drad_reject_count2;
  double draMedian2;
  double draRMS2;
  double ddecMedian2;
  double ddecRMS2;
  double dradRMS2;
} LOCALBIN, *PLOCALBIN;

/* Magnitude file header versions */
#define MAG_HEADER_VERSION_ONE   1 /* Version 1 had header bugs */
#define MAG_HEADER_VERSION_TWO   2
#define MAG_HEADER_VERSION_THREE 3 /* Version 3 adds rejectFlag (Error here: it should have been the data version below */
#define MAG_HEADER_VERSION_FOUR  4 /* Version 4 adds magdep_bin, magcal_magdep, and magcal_magdep_rms */
#define MAG_HEADER_VERSION_FIVE  5 /* Version 5 adds dra and ddec  */
#define MAG_HEADER_VERSION       6 /* Version 6 adds A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex   */
#define MAG_HEADER_DATA   0x68656164 /* 'head' */

/* Magnitude file star data versions */
#define MAG_FORMAT_VERSION_ONE   1
#define MAG_FORMAT_VERSION_TWO   1
#define MAG_FORMAT_VERSION_THREE 1
#define MAG_FORMAT_VERSION_FOUR  4  /* Version 4 adds magdep_bin, magcal_magdep, and magcal_magdep_rms */
#define MAG_FORMAT_VERSION_FIVE  5  /* Version 5 adds dra and ddec  */
#define MAG_FORMAT_VERSION       6  /* Version 5 adds A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex  */
#define MAG_FORMAT_DATA   0x64617461 /* 'data' */
#define MAG_SUBDIR_MODULUS 1000

#define MAG_FILE_MODULUS 1024

#define CACHE_MAGNITUDE_ENTRIES 15
#define POSITION_REALLOC_INCREMENT 1000
#define STATS_REALLOC_INCREMENT 1000

/* WARNING: be sure to swap bytes in photometryutils.c for any new fields */
typedef struct _filestarheader {
  long long versionTag;                            /* Version tag */
  off_t curFileSize;
  off_t fileOffset[MAG_FILE_MODULUS];   /* Offset to the start of the bin */
  int starCount[MAG_FILE_MODULUS];      /* Number of stars in the bin */
  int gsc_base_index;                   /* Base gsc bin index */
  int headersize;                       /* Header size */
  int starimagesize;                    /* Star Image size */
  int unused;                           /* Pad to 8 byte multiple */
} FILESTARHEADER,*PFILESTARHEADER;

/* Note: Needs a rejectFlag, quality, and qualityMask */
typedef struct _filestarimage_v2 {
  long long versionTag;                            /* Version tag */
  long long REFNumber;                             /* Star reference designation - first digit is catalog */
                                                   /*    For GSC2.3.2, replace a "N" with 11 and a "S" with 12 */
                                                   /*    For Kepler, replace "K" with 2 */
                                                   /*    NONE - 0 */
  double X_IMAGE;                                  /* X pixel location */
  double Y_IMAGE;                                  /* Y pixel location */
  double MAG_ISO;                                  /* Sextractor isophotal estimate */
  double ra;                                       /* Object Right Ascension in degrees */
  double dec;                                      /* Object Declination in degrees */
  double magcal_iso;                               /* Lowess magnitude estimate */
  double magcal_iso_rms;                           /* Lowess magnitude estimate error */
  double magcal_local;                             /* Local magnitude estimate */
  double magcal_local_rms;                         /* Local magnitude estimate error */
  double Date;                                     /* Julian date of exposure */
  double FLUX_ISO;                                 /* Peak Sextractor flux */
  double MAG_APER;                                 /* Sextractor aperture magnitude */
  double MAG_AUTO;                                 /* Sextractor auto magnitude */
  double KRON_RADIUS;                              /* Sextractor Kron radius */
  double BACKGROUND;                               /* Sextractor background estimate */
  double FLUX_MAX;                                 /* Sextractor maximum flux in ADU */
  double THETA_J2000;                              /* Sextractor major axis angle */
  double ELLIPTICITY;                              /* Sextractor ellipticity */
  double ISOAREA_WORLD;                            /* Sextractor area in square degrees */
  double FWHM_IMAGE;                               /* Width of the image in pixels */
  double FWHM_WORLD;                               /* Width of the image in degrees */
  double plate_dist;                               /* Plate distance from center in degrees */
  double Blendedmag;                               /* Blended magnitude */
  double limiting_mag_local;                       /* Locally corrected limiting magnitude */
  double magcal_local_error;                       /* + Local calibration correction error */
  double dradRMS2;                                 /* + clipped astrometric error for the local bin */
  double magcor_local;                             /* + Local correction magnitude */
  double extinction;                               /* + Extinction correction */
  int gsc_bin_index;                               /* GSC bin index */
  int plateNumber;                                 /* Plate number */
  int NUMBER;                                      /* Sextractor object number */
  int versionId;                                   /* Pipeline version Id */
  int AFLAGS;                                      /* Sextractor and pipeline flags */
  int BFLAGS;                                      /* Sextractor and pipeline flags */
  int ISO0;                                        /* Isophotal level 0 area in sq pixels */
  int ISO1;                                        /* Isophotal level 1 area in sq pixels */
  int ISO2;                                        /* Isophotal level 2 area in sq pixels */
  int ISO3;                                        /* Isophotal level 3 area in sq pixels */
  int ISO4;                                        /* Isophotal level 4 area in sq pixels */
  int ISO5;                                        /* Isophotal level 5 area in sq pixels */
  int ISO6;                                        /* Isophotal level 6 area in sq pixels */
  int ISO7;                                        /* Isophotal level 7 area in sq pixels */
  int npoints_local;                               /* + Points used for local calibration */
  unsigned short passBits;                         /* Pass identifier; currently only bit 1: GSC2.3.2 pass1 is recognized */
  unsigned short local_bin_index;                  /* Local bin index: ix+(nx*iy) where ix is the horizontal zero-based index */
  unsigned char seriesId;                          /* Plate series */
  char exposureNumber;                             /* + Exposure number (zero based), negative numbers for unmatched solutions. */
  unsigned char solutionNumber;                    /* Solution number */
  unsigned char spatial_bin;                       /* Spatial bin */
  unsigned char pad[4];                            /* Padding to make the sizes agree on 32 and 64 bit systems */
} FILESTARIMAGEV2,*PFILESTARIMAGEV2;

/* Version 3 adds a rejectFlag and increases padding to 8  */
/* sizeof FILESTARHEADER 12320 sizeof FILESTARIMAGE 320 */
/* WARNING: be sure to swap bytes in photometryutils.c for any new fields */
typedef struct _filestarimage_v3 {
  long long versionTag;                            /* Version tag */
  long long REFNumber;                             /* Star reference designation - first digit is catalog */
                                                   /*    For GSC2.3.2, replace a "N" with 11 and a "S" with 12 */
                                                   /*    For Kepler, replace "K" with 2 */
                                                   /*    NONE - 0 */
  double X_IMAGE;                                  /* X pixel location */
  double Y_IMAGE;                                  /* Y pixel location */
  double MAG_ISO;                                  /* Sextractor isophotal estimate */
  double ra;                                       /* Object Right Ascension in degrees */
  double dec;                                      /* Object Declination in degrees */
  double magcal_iso;                               /* Lowess magnitude estimate */
  double magcal_iso_rms;                           /* Lowess magnitude estimate error */
  double magcal_local;                             /* Local magnitude estimate */
  double magcal_local_rms;                         /* Local magnitude estimate error */
  double Date;                                     /* Julian date of exposure */
  double FLUX_ISO;                                 /* Peak Sextractor flux */
  double MAG_APER;                                 /* Sextractor aperture magnitude */
  double MAG_AUTO;                                 /* Sextractor auto magnitude */
  double KRON_RADIUS;                              /* Sextractor Kron radius */
  double BACKGROUND;                               /* Sextractor background estimate */
  double FLUX_MAX;                                 /* Sextractor maximum flux in ADU */
  double THETA_J2000;                              /* Sextractor major axis angle */
  double ELLIPTICITY;                              /* Sextractor ellipticity */
  double ISOAREA_WORLD;                            /* Sextractor area in square degrees */
  double FWHM_IMAGE;                               /* Width of the image in pixels */
  double FWHM_WORLD;                               /* Width of the image in degrees */
  double plate_dist;                               /* Plate distance from center in degrees */
  double Blendedmag;                               /* Blended magnitude */
  double limiting_mag_local;                       /* Locally corrected limiting magnitude */
  double magcal_local_error;                       /* + Local calibration correction error */
  double dradRMS2;                                 /* + clipped astrometric error for the local bin */
  double magcor_local;                             /* + Local correction magnitude */
  double extinction;                               /* + Extinction correction */
  int gsc_bin_index;                               /* GSC bin index */
  int plateNumber;                                 /* Plate number */
  int NUMBER;                                      /* Sextractor object number */
  int versionId;                                   /* Pipeline version Id */
  int AFLAGS;                                      /* Sextractor and pipeline flags */
  int BFLAGS;                                      /* Sextractor and pipeline flags */
  int ISO0;                                        /* Isophotal level 0 area in sq pixels */
  int ISO1;                                        /* Isophotal level 1 area in sq pixels */
  int ISO2;                                        /* Isophotal level 2 area in sq pixels */
  int ISO3;                                        /* Isophotal level 3 area in sq pixels */
  int ISO4;                                        /* Isophotal level 4 area in sq pixels */
  int ISO5;                                        /* Isophotal level 5 area in sq pixels */
  int ISO6;                                        /* Isophotal level 6 area in sq pixels */
  int ISO7;                                        /* Isophotal level 7 area in sq pixels */
  int npoints_local;                               /* + Points used for local calibration */
  int rejectFlag;                                  /* + localbin reject flag */
  unsigned short passBits;                         /* Pass identifier; currently only bit 1: GSC2.3.2 pass1 is recognized */
  unsigned short local_bin_index;                  /* Local bin index: ix+(nx*iy) where ix is the horizontal zero-based index */
  unsigned char seriesId;                          /* Plate series */
  char exposureNumber;                             /* + Exposure number (zero based), negative numbers for unmatched solutions. */
  unsigned char solutionNumber;                    /* Solution number */
  unsigned char spatial_bin;                       /* Spatial bin */
  unsigned char pad[8];                            /* Padding to make the sizes agree on 32 and 64 bit systems */
} FILESTARIMAGEV3,*PFILESTARIMAGEV3;

/* Version 4 adds magnitude-dependent variables  */
/* sizeof FILESTARHEADER 12320 sizeof FILESTARIMAGE 336 */
/* WARNING: be sure to swap bytes in photometryutils.c for any new fields */
typedef struct _filestarimage_v4 {
  long long versionTag;                            /* Version tag */
  long long REFNumber;                             /* Star reference designation - first digit is catalog */
                                                   /*    For GSC2.3.2, replace a "N" with 11 and a "S" with 12 */
                                                   /*    For Kepler, replace "K" with 2 */
                                                   /*    NONE - 0 */
  double X_IMAGE;                                  /* X pixel location */
  double Y_IMAGE;                                  /* Y pixel location */
  double MAG_ISO;                                  /* Sextractor isophotal estimate */
  double ra;                                       /* Object Right Ascension in degrees */
  double dec;                                      /* Object Declination in degrees */
  double magcal_iso;                               /* Lowess magnitude estimate */
  double magcal_iso_rms;                           /* Lowess magnitude estimate error */
  double magcal_local;                             /* Local magnitude estimate */
  double magcal_local_rms;                         /* Local magnitude estimate error */
  double Date;                                     /* Julian date of exposure */
  double FLUX_ISO;                                 /* Peak Sextractor flux */
  double MAG_APER;                                 /* Sextractor aperture magnitude */
  double MAG_AUTO;                                 /* Sextractor auto magnitude */
  double KRON_RADIUS;                              /* Sextractor Kron radius */
  double BACKGROUND;                               /* Sextractor background estimate */
  double FLUX_MAX;                                 /* Sextractor maximum flux in ADU */
  double THETA_J2000;                              /* Sextractor major axis angle */
  double ELLIPTICITY;                              /* Sextractor ellipticity */
  double ISOAREA_WORLD;                            /* Sextractor area in square degrees */
  double FWHM_IMAGE;                               /* Width of the image in pixels */
  double FWHM_WORLD;                               /* Width of the image in degrees */
  double plate_dist;                               /* Plate distance from center in degrees */
  double Blendedmag;                               /* Blended magnitude */
  double limiting_mag_local;                       /* Locally corrected limiting magnitude */
  double magcal_local_error;                       /* + Local calibration correction error */
  double dradRMS2;                                 /* + clipped astrometric error for the local bin */
  double magcor_local;                             /* + Local correction magnitude */
  double extinction;                               /* + Extinction correction */
  double magcal_magdep;                            /* magnitude-corrected local magnitude */
  double magcal_magdep_rms;                        /* RMS of the magnitude-corrected local magnitude */
  int gsc_bin_index;                               /* GSC bin index */
  int plateNumber;                                 /* Plate number */
  int NUMBER;                                      /* Sextractor object number */
  int versionId;                                   /* Pipeline version Id */
  int AFLAGS;                                      /* Sextractor and pipeline flags */
  int BFLAGS;                                      /* Sextractor and pipeline flags */
  int ISO0;                                        /* Isophotal level 0 area in sq pixels */
  int ISO1;                                        /* Isophotal level 1 area in sq pixels */
  int ISO2;                                        /* Isophotal level 2 area in sq pixels */
  int ISO3;                                        /* Isophotal level 3 area in sq pixels */
  int ISO4;                                        /* Isophotal level 4 area in sq pixels */
  int ISO5;                                        /* Isophotal level 5 area in sq pixels */
  int ISO6;                                        /* Isophotal level 6 area in sq pixels */
  int ISO7;                                        /* Isophotal level 7 area in sq pixels */
  int npoints_local;                               /* + Points used for local calibration */
  int rejectFlag;                                  /* + localbin reject flag */
  int magdep_bin;                                  /* Magnituded-dependent calibration bin (-1) if undefined */
  unsigned short passBits;                         /* Pass identifier; currently only bit 1: GSC2.3.2 pass1 is recognized */
  unsigned short local_bin_index;                  /* Local bin index: ix+(nx*iy) where ix is the horizontal zero-based index */
  unsigned char seriesId;                          /* Plate series */
  char exposureNumber;                             /* + Exposure number (zero based), negative numbers for unmatched solutions. */
  unsigned char solutionNumber;                    /* Solution number */
  unsigned char spatial_bin;                       /* Spatial bin */
  unsigned char pad[4];                            /* Padding to make the sizes agree on 32 and 64 bit systems */
} FILESTARIMAGEV4,*PFILESTARIMAGEV4;

/* Version 5 adds dra and ddec  */
/* sizeof FILESTARHEADER 12320 sizeof FILESTARIMAGE 320 */
/* WARNING: be sure to swap bytes in photometryutils.c for any new fields */
typedef struct _filestarimage_v5 {
  long long versionTag;                            /* Version tag */
  long long REFNumber;                             /* Star reference designation - first digit is catalog */
                                                   /*    For GSC2.3.2, replace a "N" with 11 and a "S" with 12 */
                                                   /*    For Kepler, replace "K" with 2 */
                                                   /*    NONE - 0 */
  double X_IMAGE;                                  /* X pixel location */
  double Y_IMAGE;                                  /* Y pixel location */
  double MAG_ISO;                                  /* Sextractor isophotal estimate */
  double ra;                                       /* Object Right Ascension in degrees */
  double dec;                                      /* Object Declination in degrees */
  double Date;                                     /* Julian date of exposure */
  double FLUX_ISO;                                 /* Peak Sextractor flux */
  double MAG_APER;                                 /* Sextractor aperture magnitude */
  double MAG_AUTO;                                 /* Sextractor auto magnitude */
  double KRON_RADIUS;                              /* Sextractor Kron radius */
  double BACKGROUND;                               /* Sextractor background estimate */
  double FLUX_MAX;                                 /* Sextractor maximum flux in ADU */
  double THETA_J2000;                              /* Sextractor major axis angle */
  double ELLIPTICITY;                              /* Sextractor ellipticity */
  double ISOAREA_WORLD;                            /* Sextractor area in square degrees */
  double FWHM_IMAGE;                               /* Width of the image in pixels */
  double FWHM_WORLD;                               /* Width of the image in degrees */
  double plate_dist;                               /* Plate distance from center in degrees */
  double Blendedmag;                               /* Blended magnitude */
  double dradRMS2;                                 /* + clipped astrometric error for the local bin */
  double ra_2;                                     /* Catalog star Right Ascension in degrees corrected for proper motion */
  double dec_2;                                    /* Catalog star Declination in degrees corrected for proper motion */
  float magcal_iso;                               /* Lowess magnitude estimate */
  float magcal_iso_rms;                           /* Lowess magnitude estimate error */
  float magcal_local;                             /* Local magnitude estimate */
  float magcal_local_rms;                         /* Local magnitude estimate error */
  float limiting_mag_local;                       /* Locally corrected limiting magnitude */
  float magcal_local_error;                       /* + Local calibration correction error */
  float magcor_local;                             /* + Local correction magnitude */
  float extinction;                               /* + Extinction correction */
  float magcal_magdep;                            /* magnitude-corrected local magnitude */
  float magcal_magdep_rms;                        /* RMS of the magnitude-corrected local magnitude */
  float RaPM;                                      /* Right Ascension proper motion in mas/yr */
  float DecPM;                                     /* Declination proper motion in mas/yr */
  int gsc_bin_index;                               /* GSC bin index */
  int plateNumber;                                 /* Plate number */
  int NUMBER;                                      /* Sextractor object number */
  int versionId;                                   /* Pipeline version Id */
  int AFLAGS;                                      /* Sextractor and pipeline flags */
  int BFLAGS;                                      /* Sextractor and pipeline flags */
  int ISO0;                                        /* Isophotal level 0 area in sq pixels */
  int ISO1;                                        /* Isophotal level 1 area in sq pixels */
  int ISO2;                                        /* Isophotal level 2 area in sq pixels */
  int ISO3;                                        /* Isophotal level 3 area in sq pixels */
  int ISO4;                                        /* Isophotal level 4 area in sq pixels */
  int ISO5;                                        /* Isophotal level 5 area in sq pixels */
  int ISO6;                                        /* Isophotal level 6 area in sq pixels */
  int ISO7;                                        /* Isophotal level 7 area in sq pixels */
  int npoints_local;                               /* + Points used for local calibration */
  int rejectFlag;                                  /* + localbin reject flag */
  int magdep_bin;                                  /* Magnituded-dependent calibration bin (-1) if undefined */
  unsigned short passBits;                         /* Pass identifier; currently only bit 1: GSC2.3.2 pass1 is recognized */
  unsigned short local_bin_index;                  /* Local bin index: ix+(nx*iy) where ix is the horizontal zero-based index */
  unsigned char seriesId;                          /* Plate series */
  char exposureNumber;                             /* + Exposure number (zero based), negative numbers for unmatched solutions. */
  unsigned char solutionNumber;                    /* Solution number */
  unsigned char spatial_bin;                       /* Spatial bin */
  unsigned char pad[4];                            /* Padding to make the sizes agree on 32 and 64 bit systems */
} FILESTARIMAGEV5,*PFILESTARIMAGEV5;

/* Version 6 adds A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex  */
/* sizeof FILESTARHEADER 12320 sizeof FILESTARIMAGE 336 */
/* WARNING: be sure to swap bytes in photometryutils.c for any new fields */
typedef struct _filestarimage_v6 {
  long long versionTag;                            /* Version tag */
  long long REFNumber;                             /* Star reference designation - first digit is catalog */
                                                   /*    For GSC2.3.2, replace a "N" with 11 and a "S" with 12 */
                                                   /*    For Kepler, replace "K" with 2 */
                                                   /*    NONE - 0 */
  double X_IMAGE;                                  /* X pixel location */
  double Y_IMAGE;                                  /* Y pixel location */
  double MAG_ISO;                                  /* Sextractor isophotal estimate */
  double ra;                                       /* Object Right Ascension in degrees */
  double dec;                                      /* Object Declination in degrees */
  double Date;                                     /* Julian date of exposure */
  double FLUX_ISO;                                 /* Peak Sextractor flux */
  double MAG_APER;                                 /* Sextractor aperture magnitude */
  double MAG_AUTO;                                 /* Sextractor auto magnitude */
  double KRON_RADIUS;                              /* Sextractor Kron radius */
  double BACKGROUND;                               /* Sextractor background estimate */
  double FLUX_MAX;                                 /* Sextractor maximum flux in ADU */
  double THETA_J2000;                              /* Sextractor major axis angle */
  double ELLIPTICITY;                              /* Sextractor ellipticity */
  double ISOAREA_WORLD;                            /* Sextractor area in square degrees */
  double FWHM_IMAGE;                               /* Width of the image in pixels */
  double FWHM_WORLD;                               /* Width of the image in degrees */
  double plate_dist;                               /* Plate distance from center in degrees */
  double Blendedmag;                               /* Blended magnitude */
  double dradRMS2;                                 /* + clipped astrometric error for the local bin */
  double ra_2;                                     /* Catalog star Right Ascension in degrees corrected for proper motion */
  double dec_2;                                    /* Catalog star Declination in degrees corrected for proper motion */
  float magcal_iso;                               /* Lowess magnitude estimate */
  float magcal_iso_rms;                           /* Lowess magnitude estimate error */
  float magcal_local;                             /* Local magnitude estimate */
  float magcal_local_rms;                         /* Local magnitude estimate error */
  float limiting_mag_local;                       /* Locally corrected limiting magnitude */
  float magcal_local_error;                       /* + Local calibration correction error */
  float magcor_local;                             /* + Local correction magnitude */
  float extinction;                               /* + Extinction correction */
  float magcal_magdep;                            /* magnitude-corrected local magnitude */
  float magcal_magdep_rms;                        /* RMS of the magnitude-corrected local magnitude */
  float RaPM;                                      /* Right Ascension proper motion in mas/yr */
  float DecPM;                                     /* Declination proper motion in mas/yr */
  float timeAccuracy;                              /* Estimated time accuracy in days */
  int gsc_bin_index;                               /* GSC bin index */
  int plateNumber;                                 /* Plate number */
  int NUMBER;                                      /* Sextractor object number */
  int versionId;                                   /* Pipeline version Id */
  int AFLAGS;                                      /* Sextractor and pipeline flags */
  int A2FLAGS;                                     /* Sextractor and pipeline flags */
  int BFLAGS;                                      /* Sextractor and pipeline flags */
  int B2FLAGS;                                     /* Sextractor and pipeline flags */
  int ISO0;                                        /* Isophotal level 0 area in sq pixels */
  int ISO1;                                        /* Isophotal level 1 area in sq pixels */
  int ISO2;                                        /* Isophotal level 2 area in sq pixels */
  int ISO3;                                        /* Isophotal level 3 area in sq pixels */
  int ISO4;                                        /* Isophotal level 4 area in sq pixels */
  int ISO5;                                        /* Isophotal level 5 area in sq pixels */
  int ISO6;                                        /* Isophotal level 6 area in sq pixels */
  int ISO7;                                        /* Isophotal level 7 area in sq pixels */
  int npoints_local;                               /* + Points used for local calibration */
  int rejectFlag;                                  /* + localbin reject flag */
  int magdep_bin;                                  /* Magnituded-dependent calibration bin (-1) if undefined */
  unsigned short passBits;                         /* Pass identifier; currently only bit 1: GSC2.3.2 pass1 is recognized */
  unsigned short local_bin_index;                  /* Local bin index: ix+(nx*iy) where ix is the horizontal zero-based index */
  unsigned short maskIndex;                        /* Multiple exposure mask number (using search_close) */
  unsigned char seriesId;                          /* Plate series */
  char exposureNumber;                             /* + Exposure number (zero based), negative numbers for unmatched solutions. */
  unsigned char solutionNumber;                    /* Solution number */
  unsigned char spatial_bin;                       /* Spatial bin */
  unsigned char catalogNumber;                     /* Only for the experimental catalog merge: catalog number */
  unsigned char pad[5];                            /* Padding to make the sizes agree on 32 and 64 bit systems */
} FILESTARIMAGE,*PFILESTARIMAGE;

/* Extension used by the web-based plotter */
typedef struct _filestarimageext {
  FILESTARIMAGE filestarimage;
  PHOTPLATES photplates;
  int mosaicNumber;
  int quality;
  int plateVersionId;
} FILESTARIMAGEEXT,*PFILESTARIMAGEEXT;

/* Parallel array with non-filebased additions */
typedef struct _starimage {
  PFILESTARIMAGE pFileStarImage;
  char series[MAX_SERIES_STRING]; /* Series */
  char REF[MAX_REF];              /* GSC2.3.2 reference number */
  char Plate[MAX_PLATE_NAME];     /* Plate string with solution number */
  /* The following are used by find_lightcurves */
  double Stdmag;              /* GSC2.3.2 magnitude */
  int AFLAGSCOPY;
  double colorterm;
  double errorcolor;
  int colorflag;
  int quality;                /* Used in transient searches */
  int badcolorflag;
} PHOTSTARIMAGE,*PPHOTSTARIMAGE;

/* Table used in the search for unmatched objects */
typedef struct _nonestats {
  long long REFNumber;
  double curdrad;
  double scale;  /* Degrees per pixel */
  double drad;   /* Degrees */
  double drad2; /* This drad does not go below SEARCH_NONE_MIN_DRAD pixels */
  double aLength; /* Ellipse long axis semi-length  (degrees) */
	double bLength; /* Ellipse short axis semi-length (degrees) */
	double imageArea; /* area, including drad in (square degrees) */
	double epoch;
	double blendAngle;
  double MAG_ISO; /* Tiebreaker: choose the brighter of two candidates with equal drad */
  double colorterm;
  int selected;
	int flink;  /* Index of the next group member */
  int groupNumber;
	int magnitudeIndex;
  int colorflag;
	/* The following section for late matches */
	double magcal_magdep;
	int blendCandidate;
	int blendGroupNumber;
	int blendFlink;
	int blendRefcnt;
} NONESTATS,*PNONESTATS;

/* Table used for sorting NONE objects in increasing effective area */
typedef struct _AreaSort {
	int magnitudeIndex;
	PNONESTATS pNoneStats;
} AREASORT,*PAREASORT;

/* Table used to calculate proper motions */
typedef struct _pmstats {
  long long REFNumber;
	double ra;
	double dec;
	double catra;
	double catdec;
	double ra_2;
	double dec_2;
	double RaPM;
	double DecPM;
	double epoch;   /* Year */
	double ELLIPTICITY; /* Used to eliminate trailed images */
	int count;      /* Entry in this group */
	int link;       /* Link to next entry for this REFNumber */
	int AFLAGS;
	int BFLAGS;     /* Used to eliminate saturated points */
} PMSTATS,*PPMSTATS;

/* Table used in the search for unmatched objects */
typedef struct _groupstats {
	long long REFNumber;
  char nearbyObjects[MAX_NEARBY_OBJECTS_STRING];
  int groupCount;
  int groupIndex;
	int groupNumber;
	int flink; /* Index of the first/next group member */
  int GSCBlendFlag;
  int magcal_magdep_count; /* Number of stars for magcal_magdep estimate */
  int gsc_bin_index; /* GSC bin index of this group */
  double ra; /* Position of the lowest drad star in this group */
  double dec;
	double factor; /* cosine of declination */
	double drad2;
  double magcal_magdep;  /* Estimate of average magcal_magdep */
	/* The following fields are for comparisons with nearby catalog objects that may have high proper motions.
		 we take the average of the earliest and latest observation to calculate the catalog object position rather
		 than the J2000 position */
	double minEpoch; /* earliest date of the group */
	double maxEpoch; /* latest date of this group */
	double maxDec; /* Maximum declination in this group */
	double minDec; /* Minimum declination in this group */
	double maxRA; /* Maximum Right Ascension in this group */
	double minRA; /* Minimum Right Ascension in this group */
  NEARESTCATALOGSTAR nearestCatalogStar;
	int blendGroupNumber;
	int blendFlink;
	int blendGroupMemberCount;
} GROUPSTATS,*PGROUPSTATS;

/* Table used in the search for unmatched objects */
typedef struct _nearbycatalog {
  long long REFNumber;
  /* Values for good objects */
  double ra;
  double dec;
  double dradRMS2;
  double magcal_magdep;
  int magcal_magdep_count;
  /* Values for objects objects above the limiting magnitude */
  double raGoodMag;
  double decGoodMag;
  double dradRMS2GoodMag;
  double magcal_magdepGoodMag;
	double RaPM;
	double DecPM;
  int magcal_magdepGoodMag_count;
	int gsc_bin_index;
	int firstProperMotionIndex; /* Chain of proper motion values */
	int lastProperMotionIndex;
	int properMotionCount;
} NEARBYCATALOG, *PNEARBYCATALOG;

#define MAX_DRAD_HISTOGRAM_MAG 20
#define MAX_DRAD_RANGE  3600 /* 0.1 arcsec steps */
#define WEB_TIMEOUT 31536000   /* Timeout for web activity (Disable by setting this to 1 year) */
#define DEFAULT_GAL_LATITUDE 40.0 /* For search_none processing */
#define MAX_CATALOG_STRING 40
#define MAX_CPUS 10
#define NETWORK_MULTIPLIER 1000L
#define DASCH_MAX_PROCESSOR_NAME 257 /* Must be greater than MPI_MAX_PROCESSOR_NAME (enforced by an assert()) */
#define DASCH_MAX_FILENAME  512
#define TIME_ACCURACY_ALLOCATION 1000
#define MAX_TC_DIST 200 /* Size of catalogDistance histogram in arcsec */

typedef struct _filecounterblock {
  off_t magnitudeBufferAlloc;
  long long returnedMagnitudes; /* Number of magnitudes returned */
  long long readMagnitudes;     /* Number of magnitudes read */
  long long currentMagnitudes;  /* Number of valid versionId checked magnitudes */
	double maxPMrarms;  /* Maximum rms in proper motion correction */
	double maxPMdecrms; /* Maximum rms in proper motion correction */
  /* Cache counters */
  int cache_access_count;
  int cache_hit_count;
  int cache_miss_count;
  int cache_magnitude_read_count;
  int cache_magnitude_reject_count;
  /* Counters */
  int staleVersionCount;
  int versionIdQueries;
  int versionIdHits;
  int binRejectCount;
  int nomatchRejectCount;
  int flagsRejectCount;
  int refRejectCount;
  int sequesteredCount;
  int binZeroErrors;
	int properMotionErrorCount;
  int noHeaderCount; /* File has no header */
  int versionOneFile; /* Version one had a bug */
  int corruptStarCount; /* star count table is incorrect */
  int newDataCount; /* New data appended to the file */
  int noValidImages; /* No images found in the file */
  int qualityRejectCount; /* Number of images rejected because of plate quality */
  int chipRejectCount; /* Number of images not on the Kepler chips */
  int subroutineCalls; /* Times GetFileSummaryMagnitudes has been called */
  int excludeSeriesCount; /* Number of images rejected because of bad astrometry series */
  int filteredVersionCount; /* Number of images rejected by request */
  int mixedCatalogCount; /* Number of images rejected because of mixed catalogs */
  int staleREFCount; /* Number of matched catalog stars no longer in the catalog (typically a new APASS version) */
  int rejectedTCGroup; /* Number of Transient Candidate blend groups rejected because of no suitable root */
  int limitingmag1; /* new limiting magnitude definition of Sumin's memo of Thu 12/15/11 1:52 PM */
  int limitingmag2; /* new limiting magnitude definition of Sumin's memo of Thu 12/15/11 1:52 PM */
	double maxGroupDrad; /* Maximum group extent found */
	int AFLAGS_ZERO_COUNT;       /* count of root Transient Candidate objects with zero AFLAGS */
	int AFLAGS_NONZERO_COUNT;    /* count of root TC objects with nonzero AFLAGS */
	int AFLAGS_MULT_COUNT[AFLAG_BIT_COUNT];   /* count of root TC AFLAGS per-bit settings*/
	int AFLAGS_UNIQUE_COUNT[AFLAG_BIT_COUNT]; /* unique count of root TC AFLAGS per-bit settings */
	int maxGroupHistogram[MAX_GROUP_DRAD_HISTOGRAM+1];
} FILECOUNTERBLOCK,*PFILECOUNTERBLOCK;

typedef struct _timeAccuracyEntry {
  double singleJulianDate;
  double singleTolerance;
  double allTolerance;
  double returnedTolerance;
  int flink;  /* flink is an index into the seriesTimeAccuracyTable[0], but is never zero */
  char valid;   /* Entry is valid if nonzero */
} TIMEACCURACYENTRY, *PTIMEACCURACYENTRY;
#define MAX_TIMESTR_SIZE 26

#define MAX_FLARECOMMENT_SIZE 1000
#define MAX_FLARECANDIDATE_SIZE  700
#define MAX_FLARESUBDIR_NAME 10
#define MAX_FLAREBUFFER 2048
#define MAX_FLAREFILE_DISTANCE 14

typedef struct _flarecandidate {
  long long REFNumber;
  double ra;
  double dec;
  int refType;
  int unrecognized;
  int interestingFlag; /* if 1, on Josh's list */
  int previousFlag; /* if 1, on the flarecandidate_name list */
  int input_index; /* Pointer to the matching object in the input table */
  int flarecandidate_index;
  int interesting_peakEvaluation;  /* PEAKEVALUATION bitmap defined in photoutils.h and photutils.c that shows reasons for transient rejections */
  int prevPeakCount;
  int newPeakCount;
  int gsc_bin_index;
  char typeFlag;
  char allNameOutputFlag; /* Set if this flarecandidate was written to the <all_name> table */
  char runTimestamp[MAX_FLARESUBDIR_NAME+1]; /* This is an old runDate which may or may not be equal to the current runDate */
  char name[MAX_FLARECANDIDATE_SIZE+1];
  char comment[MAX_FLARECOMMENT_SIZE+1];
} FLARECANDIDATE,*PFLARECANDIDATE;

typedef struct _filecommon {
	time_t startTime;
  MYSQL* pConnection;
  MYSQL* pPhotConnection;
  sem_t *fileSemaphore;
  FILE *debugHandle;
	FILE *propermotionHandle;
  FILECOUNTERBLOCK counterBlock;
  GALAXYCOMMON galaxyCommon;
  /* File read buffer, output buffer */
  PFILESTARIMAGE magnitudeBuffer;
  int curMagnitudes;            /* Number of valid magnitudes in the magnitudeBuffer */

  /* Access to catalog data */
  int catalogFD;                /* Catalog data file descriptor number */
  int indexFD;                  /* Catalog index file descriptor number */
  int gsc_bin_index;            /* Catalog index in cache */
  int catalogTableSize;         /* Number of catalog records */
  int catalogTableAlloc;        /* Number of allocated catalog records */
  PGSCIMAGE  catalogTable;      /* cached catalog data, sorted by REFNumber */

  /* Series version information */
  int seriesTimeAccuracyCount;
  int seriesTimeAccuracyAlloc;
  PTIMEACCURACYENTRY seriesTimeAccuracyTable[MAX_SERIES+1];
  int *seriesVersionTable[MAX_SERIES+1];
  int *seriesQualityValueTable[MAX_SERIES+1];
  int *seriesQualityMaskTable[MAX_SERIES+1];
  double *seriesColortermTable[MAX_SERIES+1];
  int *seriesColorflagTable[MAX_SERIES+1];
  float **seriesLimitingMagnitudeArray;
  int maxPlateNumber[MAX_SERIES+1];
  int catalogNumber;
  int minVersionId;
  int currentMagnitudesBin; /* Current gsc bin to which to add currentMagnitudes */
  /* The Magnitude cache */
  int cache_current_sequence;
  int cache_gsc_bin_index[CACHE_MAGNITUDE_ENTRIES];  /* The smallest gsc bin index in this entry */
  int cache_bin_magnitudes[CACHE_MAGNITUDE_ENTRIES]; /* Number of valid magnitudes in this entry */
  int cache_access_sequence[CACHE_MAGNITUDE_ENTRIES]; /* Access sequence for this entry (oldest gets deleted) */
  PPHOTSTARIMAGE pCacheMagnitudeTable[CACHE_MAGNITUDE_ENTRIES];            /* List of valid magnitudes for this entry */
  PFILESTARIMAGE pCacheFileMagnitudeTable[CACHE_MAGNITUDE_ENTRIES];        /* Associated file entries for the magnitudes table */

  double ra0; /* Right ascension (for debugging) */
  double dec0; /* Declination (for debugging) */
  double rad0; /* Search radius in arcsec (for debugging); */

  /* Vectors for calculating the median of catalog stars for search_none */
  PNONESTATS none_stats_table;
	PPMSTATS pm_stats_table; /* Proper motion statistics */
	PAREASORT area_sort_table; /* Area sort pointers */
  PGROUPSTATS group_stats_table;
  PNEARBYCATALOG nearby_catalog_table; /* Used only by PlotTransientCandidate() */
  PPHOTSTARIMAGE pAllMagnitudeTable;   /* Used only by PlotTransientCandidate() */

  int allMagnitudeCount;
  int catalogDistTotal[MAX_TC_DIST];  /* Histogram of the distances of DASCH objects from nearby catalog objects */
  int catalogDistSimilarMag[MAX_TC_DIST]; /* Histogram of the distances of DASCH objects from nearby catalog objects within MAX_TC_NEARBY_MAG magnitudes */
  int catalogDistDifferentMag[MAX_TC_DIST]; /* Histogram of the distances of DASCH objects from nearby catalog objects differing by more than MAX_TC_NEARBY_MAG magnitudes */
  int nearbyCatalogCount;
  int statsAllocCount;
	/* Vectors for finding miscellaneous proper motion statistics */
  /* Jul 13, 2015 - used also to find the new TC REFNumber */
	double *vector1;
	double *vector2;
	double *vector3;
	int vectorAllocCount;

	/* For the following histograms the second index is in 0.1 arcsec units.  [5][MAX_DRAD_RANGE] is 0.0-0.1 arcsec for magnitudes 5-6.0 */
	long long  draHistogram[MAX_DRAD_HISTOGRAM_MAG][2*(MAX_DRAD_RANGE)+1];
	long long  ddecHistogram[MAX_DRAD_HISTOGRAM_MAG][2*(MAX_DRAD_RANGE)+1];

  int numGroups;
	int blendGroupCount;
  /* Open MPI fields */
  int myrank;             /* Process rank (0 = master) */
  char processorName[DASCH_MAX_PROCESSOR_NAME];  /* CPU name */
  int namelength;                     /* Length of CPU name */
  /* Fields used by MPI search_none */
  double minimumGalacticLatitude;
  int verbose;
  int provideWebSummary; /* used by update_summary2.c for monthly web statistics */
  int histogramTable[MAX_NONE_HISTOGRAM];
  int minimumImageCount;
  int maxproc;         /* Number of processes active */
  int workQueueSize;   /* Size of the work queue */
  int startrank;       /* Number of our process */
  time_t beginTime;  /* Overall execution time */
  int canceledFlag;      /* If nonzero, the user is requesting an abort */
  int completeCount;     /* Number of successfully competed entries in an iteration */
  int queueIteration;  /* Work queue iteration */
  int remoteCPUTotal;
  int statusRequestNumber; /* Number of status request cycles */
  int extraColumnFlag;        /* When enabled, add extra columns not in the MySQL starsummary table */
  int enableTransientSearch;  /* When enabled, change observation selection parameters to assist transient search */
  int enableRematch;          /* When enabled, repeat star matching to optimize location of transients */

  /* The following are for Josh's flare candidate list */
  int flareCandidatesOnly;    /* When enabled, search only the regions near flare candidates */
  int interesting_count;
  int interesting_index;
  PFLARECANDIDATE candidate_table;
  /* End of flare candidate list handling */

  int gscBinIndexCount; /* Number of gsc bins to search */
  int *gscBinIndexList; /* List of gsc bins to search */
  char remoteProcessorName[MAX_CPUS][DASCH_MAX_PROCESSOR_NAME];  /* CPU name */
  char keepalivename[DASCH_MAX_FILENAME];
  char timestr[MAX_TIMESTR_SIZE]; /* Time stamp for transient candidate search */
  /* The following is used by ProcessMatchedImages() */
  int gscImageAlloc;
  PGSCIMAGE pGscImageTable;
} FILECOMMON,*PFILECOMMON;

typedef struct _phottarget {
  char REF[MAX_REF*MAX_BLEND_COUNT];   /* GSC2.3.2 reference number */
  char src_name[MAX_SRC_LENGTH*MAX_BLEND_COUNT];   /* GSC2.3.2 reference number */
  char nearestREF[MAX_REF]; /* Nearest matched REFNumber */
  char nearbyObjects[MAX_NEARBY_OBJECTS_STRING];
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
  int haveLocation;   /* Nonzero if the above ra and dec are valid */
  int updateflag;
  int processed;
  int daschFlag;      /* If non zero, REF_TYPE_DASCH */
  int gsc_bin_index;
  int coverage_bin_index;
  int coverage_table_index;
  int groupCount; /* Written by ProcessNoneImagesX and web_query */
  double Stdmag;
  double color;
  char class;  /* This is gscclass in the MySQL database */
  char VFlag;
  char MAGFlag;
  double RaPM;
  double DecPM;
  double magcal_magdep; /* Written by ProcessNoneImagesX and web_query */
#if 0
  double nearestREFarcsec; /* Nearest match REFNumber distance */
  double nearestREFmag; /* Nearest match REFNumber magnitude */
  int nearestREFflag;   /* Set if this object is has questionable astrometry */
#endif
  NEARESTCATALOGSTAR nearestCatalogStar;
} PHOTTARGET,*PPHOTTARGET;

typedef struct _sky2kstar {
  double ra;           /* Right ascension */
  double dec;          /* Declination */
  double magnitude;    /* Magnitude */
  double radius;       /* Radius of incomplete GSC2.3.2 ares in degrees */
} SKY2KSTAR,*PSKY2KSTAR;

/* The following is used for the series buttons on the lightcurve plotter */
#define START_SERIES_BIT 0 /* Start with bit one */
#define MAX_SERIES_BIT  30 /* Maximum series in the bitmask */
typedef struct serieslist {
  int numPlates;
  int bitPosition;      /* Bit position of this plate */
  int bitMask;         /* Bit mask */
} SERIESLIST,*PSERIESLIST;

/* The following structure is used in ProcessLightcurveParameters() */
typedef struct _sumstarimage_old {
  char series[MAX_SERIES_STRING]; /* Series */
  int seriesId;
  int plateNumber;
  int solutionNumber;
#if 0
  int magdep_bin;
  int local_bin_index;
#endif
  char REF[MAX_REF];          /* GSC2.3.2 reference number */
  long long REFNumber;        /* Translated reference number */
#if 0
  int NUMBER;                 /* Sextractor reference number */
#endif
  int passBits;               /* Pass identifier */
#if 0
  double X_IMAGE;
  double Y_IMAGE;
#endif
  int AFLAGS;                 /* Flags keyword. See header of pipelineutils.h */
  int BFLAGS;                 /* Flags keyword. See header of pipelineutils.h */
#if 0
  int npoints_local;
  double MAG_ISO;
#endif
  double ra;                  /* Right Ascension in degrees */
  double dec;                 /* Declination in degrees */
  double magcal_iso;          /* lowess magnitude estimate */
  double magcal_iso_rms;      /* local  error */
  double magcal_local_error;
  double magcal_local_rms;    /* Overall error */
  int spatial_bin;            /* Lowess magnitude spatial bin  number */
  double Date;                /* Heliocentric Julian Date  */
  double limiting_mag_local;  /* Limiting magnitude */
  double magcal_magdep;       /* Magnitude-dependent corrected magnitude */
#if 0
  double magcal_local;        /* Local magnitude calibration */
  double magcal_magdep_rms;   /* Error of the magnitude-dependent correction */
  double extinction;          /* extinction correction */
  double Stdmag;              /* GSC2.3.2 magnitude */
  double color;               /* GSC2.3.2 color */
  double dra;
  double ddec;
  double FLUX_ISO;
  double MAG_APER;
  double MAG_AUTO;
  double KRON_RADIUS;
  double BACKGROUND;
  double THRESHOLD;
  double FLUX_MAX;
  double THETA_J2000;
  double ELLIPTICITY;
  double ISOAREA_WORLD;
  double FWHM_IMAGE;
  double FWHM_WORLD;
  int ISO0;
  int ISO1;
  int ISO2;
  int ISO3;
  int ISO4;
  int ISO5;
  int ISO6;
  int ISO7;
  double plate_dist;
  double Blendedmag;
#endif
  int gsc_bin_index;
  int versionId;
  double dradRMS2;
} SUMSTARIMAGE_OLD,*PSUMSTARIMAGE_OLD;
/* The following structure is used in ProcessLightcurveParameters() */

#define PARAMETERVECTOR_ALLOC_INCREMENT 500
#define PARAMETERVECTOR_MAGCAL_MAGDEP          0 /* Vector for locally corrected magnitudes (magcal_magdep)   old vector1  */
#define PARAMETERVECTOR_MAGCAL_MAGDEP2         1 /* Use this one for clipping and sorting */
#define PARAMETERVECTOR_MAGCAL_MAGDEP3         2 /* For smoothing, with max and min replaced by median */
#define PARAMETERVECTOR_MAGCAL_MAGDEP4         3 /* This one is sorted */
#define PARAMETERVECTOR_JULIAN_DAY            4 /* Vector for julian days                                   old vector2  */
#define PARAMETERVECTOR_MAGCAL_ISO            5 /* Vector for lowess magnitudes (magcal_iso)                old vector3  */
#define PARAMETERVECTOR_MAGCAL_LOCAL_RMS      6 /* Vector for magcal_local_rms                              old vector4  */
#define PARAMETERVECTOR_DRAD                  7 /* Vector for drad estimate                                 old vector5  */
#define PARAMETERVECTOR_BDRAD                 8 /* Vector for dradB estimate                                old vector5b */
#define PARAMETERVECTOR_NONDAMON              9 /* Vector for non-Damon                                     old vector7  */
#define PARAMETERVECTOR_DAMON                10 /* Vector for Damon                                         old vector8  */
#define PARAMETERVECTOR_RA                   11 /* Vector for Right Ascension */
#define PARAMETERVECTOR_DEC                  12 /* Vector for Declination */
#define PARAMETERVECTOR_ERRREAL              13 /* ave of magcal_local_error and magcal_iso_rms added in quadrature */
#define PARAMETERVECTOR_PMAG                 14 /* Estimate of upper error bar magnitudes */
#define PARAMETERVECTOR_MMAG                 15 /* Estimate of lower error bar magnitudes */
#define PARAMETERVECTOR_LIMITING_MAG         16 /* Limiting magnitude */
#define PARAMETERVECTOR_DRADRMS2             17 /* rms astrometric error of the local smoothing bin */
#define PARAMETERVECTOR_XVALS                18 /* Integer run for smoothing */
#define PARAMETERVECTOR_DATEREAL             19 /* Modified hJd */
#define PARAMETERVECTOR_MAGS                 20 /* Lowess smoothed, span of 0.4 */
#define PARAMETERVECTOR_MAGSB                21 /* Lowess smoothed, span of 0.8 */
#define PARAMETERVECTOR_MAGSC                22 /* Sgolay smoothed, span of 10 */
#define PARAMETERVECTOR_MAGSD                23 /* loess  smoothed, span of 15 */
#define PARAMETERVECTOR_IDTREND              24 /* Used to calculate lightcurverms<n> */
#define PARAMETERVECTOR_SCATTER              25 /* Used to calculate lightcurverms<n> */
#define PARAMETERVECTOR_MAGCAL_MAGDEP5       26 /* Used for finding the median of the temp_ngood population */
#define MAX_PARAMETERVECTOR                  27

#define PARAMETERVECTOR_MALMQUIST        0  /* Vector for the Malmquist_factor                         old vector6  */
#define PARAMETERVECTOR_BMALMQUIST       1  /* Vector for the Malmquist_factor                         old vector6b */
#define PARAMETERVECTOR_SELECTION        2  /* Selection vector */
#define PARAMETERVECTOR_SLOPE_ALL        3  /* Used to calculate slope_all */
#define PARAMETERVECTOR_SLOPE_60         4  /* Used to calculate slope_all_60 */
#define PARAMETERVECTOR_PEAKVALUES       5  /* Used by ProcessTransientCandidates to search for peaks in the lightcurve */
#define PARAMETERVECTOR_NPOINTS          6  /* Used by ProcessTransientCandidates to calculate peakExtra */
#define MAX_MALMQUIST_PARAMETER          7

#define PARAMETERVECTOR_IBURST3          0
#define PARAMETERVECTOR_IDIP3            1
#define PARAMETERVECTOR_IBURST4          2
#define PARAMETERVECTOR_IDIP4            3
#define PARAMETERVECTOR_IBURST5          4
#define PARAMETERVECTOR_IDIP5            5
#define MAX_IPARAMETERVECTOR             6
#define NEW_LIMITING_MAG 1

/* Selection flags in the MALMQUIST structure */
#define SELECT_TEMP_GOOD   1 /* Meets AFLAGS criteria */
#define SELECT_TEMP_BLEND  2 /* Meets blend criteria */
#define SELECT_GOOD        4 /* Final good image criteria */
#define SELECT_BLEND       8 /* Final blend criteria */
#define SELECT_TRANSIENT  16 /* Checked for a Transient Candidate */

/* The Malmquist structure is used see whether mag_local is a function of limiting_mag_local
 * The structure was expanded for the additional parameters requested by Sumin Tang on 7/03/10 10:53
 */
#define MALMQUIST_STARS 20
#define DAMON_STARS 20
#define MIN_SLOPESTARS 20

typedef struct _Malmquist {
  double limiting_mag_local;
  double magcal_magdep;
  double Date;  /* Heliocentric Julian Day */
  double drad;
  double dradRMS2;
  int imageIndex;  /* Index of this star in the magnitudes table */
  int secondaryIndex; /* Secondary index referencing another Malmquist structure */
  int selectionFlags; /* Selection flags listed above */
  int AFLAGS;
} MALMQUIST,*PMALMQUIST;

/* the following structure is used by ProcessTransientCandidates and parseidtable.c*/

#define PLOT_LIMITING_MAGNITUDES 1
#define MIN_TC_POINTS 3 /* Minimum to be included in search_none */
#define LOW_TC_POINTS 4  /* Minimum to avoid exclusion in parseidtable.c */
#define MAX_TC_POINTS 50 /* George Miller's criteria for npoints */
#define MIN_TC_DAYS 1.0  /* Jul 13, 2015 now for all of the points, not just the points after the peak. */
#define MAX_TC_DAYS 100.0
#define MIN_TC_MAGNITUDE 0.25  /* Selection of possible flare peaks */
#define LARGE_TC_MAGNITUDE 0.5  /* Selection of largest flares */
#define LONG_TC_PRE_DAYS  (5.0*365.0) /* Averaging interval for the peakLongOutburst  baseline */
#define LONG_TC_FLARE_DAYS 365.0      /* Average interval to search for peakLongOutburst */
#define LONG_TC_SKIP_DAYS (2.0*365.0) /* Days to skip before searching for a second flare */
#define LONG_TC_MIN_DAYS (LONG_TC_FLARE_DAYS+LONG_TC_SKIP_DAYS-0.5)
#define LONG_TC_POINTS    10    /* Minimum points for peakLongOutburst average */
#define LONG_TC_MAGNITUDE  1.0  /* Minimum magnitude increase for peakLongOutburst average */
#define LONG_TC_MAX_EVENTS 100   /* Disable this until we understand the function better */
#define MIN_TC_DRAD 9 /* arcsec */
#define MIN_TC_LIMITING_YEARS 22   /* Minimum number of years with at least one magnitude observation or limiting magnitude observation (run of 3/17/17)  */
                                   /* (17 points suggested by run of 3/13/17 */
#define MIN_TC_LIMITING_POINTS 50  /* Minimum number of points, either magnitude observation or limiting magnitude observation (run of 3/17/17) */
#define MAX_PEAK_RATIO 0.4
#define MAX_TC_BRIGHT_ARCSEC 360 /* Limit to reject all 10th mag or brighter stars from the Tycho2 catalog */
#define MAX_TC_BRIGHT_MAG     10
#define MAX_TC_NEARBY_ARCSEC  18 /* Limit to reject all catalog stars within 2 magnitudes of the transient */
#define MAX_TC_NEARBY_MAG      2 /* Limit for MAX_TC_NEARBY_ARCSEC and also for Stdmag of catalog stars */
#define MIN_TC_CLIP_NGOOD 100 /* Points needed to establish a minimum magnitude for the star.  Median points is now 31.  30% of the lightcurves have more than 100 points */
#define MAX_TC_EXCESS_POINTS 10 /* Points outside the MAX_TC_DAYS detection window */
#define MAX_TC_MULTIPLE_POINTS 1 /* In effect when INCLUDE_TRANSIENT_ONE_MULTIPLE is defined */
#if 1
#define MIN_TC_DECLINATION -90.0
#define MAX_TC_DECLINATION  90.0
#define MIN_TC_RIGHTASCENSION 0.0
#define MAX_TC_RIGHTASCENSION 360.0
#else
#define MIN_TC_DECLINATION -90.0
#define MAX_TC_DECLINATION  15.0
#define MIN_TC_RIGHTASCENSION   0.0
#define MAX_TC_RIGHTASCENSION 285.0
#endif
#define MIN_TC_DATERANGE 0.0001141 /* Minimum allowable difference in date between two TC candidates for the 'sameday' flag  1 hour/((24 hours/day) * 365.25 days per year) = 0.0001141 years */
#define MF_TC_SERIESID 27  /* MF series */
#define MF_TC_PLATE_HEIGHT  21953 /* Pattern id 10, 11, and 12 */
#define MF_TC_EXCLUSION_PIXELS 1691 /*  3737 = 1.9 degrees / (1.83 arcsec/pixel); 1691 for 0.86 degree mean width */
#define MAX_TC_PEAKDRADRMS3 5.0  /* Reduced from 10 to 5 on April 16, 2019 */ /* See ~/backup/2016_12_20/peakDradRMS3.png */
#define PLOT_TC_MAG_MARGIN 0.5 /* Margins of the upper two panels in magnitudes */
#define PLOT_TC_MAG_RANGE  3.0 /* Range of the upper two panels measured from the minimum magnitude */

typedef struct _transientCandidateCommon {
   /* Input values */
  int npoints;
  int full_ngood;  /* ngood */
  int minGoodStars; /* Minimum good stars requested by user */
  double clip_med; /* clip_median_local */
  double clip_rms; /* clip_rms_local */
  int clip_ngood;
  /* Return Values */
  double peakDays; /* Number days after the peak */
  double peakYear; /* Year of the peak */
  double peakMag;  /* Brightness of the peak */
  double peakSlope; /* Best fit slope from peak to end (mag/day) */
  double peakSlopeRMS; /* Residual RMS of the best fit slope from peak to end (mag/day) */
  double peakRMS; /* Residual RMS of the best linear fit (mag) */
  double peakDradRMS; /* Positional RMS of selected points in the peak */
  double peakDradRMS2; /* dradRMS2 points added in quadrature for the selected points in the peak */
  double peakDradRMS3; /* Positional RMS of selected points in the peak referred to peakRA and peakDec */
  double peakRA;  /* 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window */
  double peakDec; /* 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window */
  double peakNearbyDistance; /* distance to nearby stars that are within 2 magnitudes of the brightest or average magnitude of the transient candidate */
  int peakCount;   /* Number of good lightcurve points inside the selected transient window */
  int peakCountNF;  /* Number of good lightcurve points inside the selected transient window from narrow field telescopes  */
  int peakCountWF;  /* Number of good lightcurve points inside the selected transient window from wide field telescopes */
  int peakOutside;  /* Number of good lightcurve points outside the selected transient window (ngood - peakCount) */
  int peakOutsideNF;  /* Number of good lightcurve points outside the selected transient window from narrow field telescopes */
  int peakOutsideWF;  /* Number of good lightcurve points outside the selected transient window from wide field telescopes */
  int peakEventCount; /* number (order) of peaks that trigger an event */
  int peakEventCount2; /* number of peaks with at least one point in a given MAX_TC_DAYS period  */
  int peakEventCount3; /* total number of peaks */
  int peakBadColorCount; /* total number with bad or undertain color */
  int peakNumber;  /* The brightness order of this peak */
  int peakUpperCount; /* Number of points above (peakMaxMag+peakMinMag)/2 */
  int nearbyREFflag;
  int peakExtra;   /* Number of points outside the flare candidate */
  int peakDefectCount; /* the number of lightcurve defects in the selected transient window */
  int peakMultipleCount; /* the number of multiple exposure plates in the selected transient window */
  int peakLimitingYears;  /* The number of years of coverage at or below peakMag (including the transient year)*/
  int peakLimitingPoints; /* The number of plates of coverage at or below peakMag (including the transient points) */
  int peakEvaluation; /* PEAKEVALUATION bitmap defined in photoutils.h and photutils.c that shows reasons for transient rejections */
  int peakLongOutburst; /* The number of times in which the average mag for 10 points during a year increases by 1 mag for the next year */
  int gsc_bin_index; /* current bin of the transient */
  double peakMaxMag; /* Maximum (brightest) B magnitude */
  double peakNoDefectMag; /* Maximum (brightest) B magnitude excluding defect points */
  double peakMinMag; /* Minimum (dimmest) B magnitude */
} TRANCOMMON,*PTRANCOMMON;

typedef struct _parametervectorstorage {
  int vectorAlloc;   /* Size of the vectors */
  double     *vector[MAX_PARAMETERVECTOR];
  PMALMQUIST mvector[MAX_MALMQUIST_PARAMETER];
  int        *ivector[MAX_IPARAMETERVECTOR];
} PARAMETERVECTORSTORE,*PPARAMETERVECTORSTORE;

/* This table is used to replace all other catalog identifiers with APASS catalog identifiers for the combined GSC/APASS catalog */
typedef struct _substitution {
  long long badREFNumber;
  long long goodREFNumber;
  int count;
  int valid;
} SUBSTITUTION,*PSUBSTITUTION;

#define MAX_SERIES_CHAR 26
#define SERIES_ALLOC_INCREMENT 150

typedef struct _seriesTree {
  int level;    /* Level of the tree, designates character in the series string */
  int seriesId; /* Index to the previous character */
  int treeIndex; /* Index to this tree */
  int charIndex[MAX_SERIES_CHAR]; /* If 0, no such series */
                                  /* If > 0, index to the next tree */
                                  /* If < 0 the negative of the seriesId */
  char charVal;
} SERIESTREE,*PSERIESTREE;

#define NUM_DAMON_SERIES 6

extern int damonSeriesId[NUM_DAMON_SERIES]; /* Series id for the damon series */
extern int excludeSeriesId[NUM_EXCLUDE_SERIES]; /* Series to exclude */

extern char* colortable[COLORFLAG_MAX];

char *GetPhotFileBase(char *catalogString);
int GetPhotPlate(MYSQL *pPhotConnection,char *series,int plateNumber,PPHOTPLATES pPhotPlates,char *catalogString,int readOnly);
int GetPhotometryGlobal(MYSQL *pPhotConnection,PPHOTGLOBAL pPhotGlobal);
int UpdatePhotometryGlobal(MYSQL *pPhotConnection,PPHOTGLOBAL pPhotGlobal,int new_version,char *versionname);
int TimeStampPhotometryGlobal(MYSQL *pPhotConnection);
int GetPlateScale(MYSQL *pConnection,char * series,int plateNumber,int mosaicNumber,int solutionNumber,int *pMosaicWidth,int *pMosaicHeight,double *pScale);
int GetStarEntry(MYSQL *pPhotConnection,PSTARENTRY *ppStarEntry,int maxstars,int updateFlag,char *REF,int gsc_bin_index,char *catalogString,int verbose,int nospace);
int GetStarEntry2(MYSQL *pPhotConnection,PSTARENTRY pCurStarEntry,char *REF, int verbose,int nospace);
int GetStarEntry3(MYSQL *pPhotConnection,char *catalogname,double *pra,double *pdeclination);
int  ReadPhotSpatialBin(MYSQL *pPhotConnection,PPHOT_SPATIAL_BIN *ppSpatialBinTable,char *series,int plateNumber,int solutionNumber,int spatial_bin,char *catalogString,int solutionNumberFix);
int GetSeriesId(char *series,int fatal);
void InitSeriesTable(MYSQL *pConnection,MYSQL *pPhotConnection);
int  GetPhotLocalBin(MYSQL *pPhotConnection,LOCALBIN *pLocalBinTable,int localBinTableSize, char *series,int plateNumber,int solutionNumber,int local_bin_index,char *catalogString,int solutionNumberFix,int debugMode);
double GetFittedPlateScale(int seriesId,int plateNumber);
int GetSequesteredFlag(int seriesId);
PSERIESENTRY GetSeriesEntry(int seriesId);
void InitFileCommon(FILE *consoleHandle,PFILECOMMON pFileCommon,MYSQL * pConnection,MYSQL * pPhotConnection,int catalogNumber,int verbose);
int InitGscBinList(PFILECOMMON pFileCommon,char *gscbinfile);
void FreeFileCommon(PFILECOMMON pFileCommon,int verbose);
int GetFileSummaryMagnitudes(PGSCBIN pGscBin,PFILECOMMON pFileCommon,int gsc_bin_index_low,int gsc_bin_index_high,int zeroAFLAGS,int REFonly,int magnitudeLimit,int dumpAllFlag,FILE *headerHandle,char *catalogString,int verbose,FILE * repairHandle,int readOnly,int nosequestered);
int GetLatestVersionId(PFILECOMMON pFileCommon,unsigned char seriesId,int plateNumber,int versionId,long long REFNumber,int fatal);

void ProcessNoneImagesX(
  PGSCBIN pGscBin,
  PFILECOMMON pFileCommon,
  PPHOTTARGET *pTarget_table,
  size_t *pTarget_nrecs,
  size_t *pTarget_alloc,
  PFILESTARIMAGE pNoneMagnitudeTable,
  int magnitudeCount,
  int verbose,
  int *histogramTable,
  int *pGroupNumber,
  int curGscBinIndex
);

void GetSky2kCatalog(FILE *consoleHandle,PSKY2KSTAR *pSky2KTable,int *pSky2k_size,double *pMaxGSCRadius);
void InitMaxPlateNumber(MYSQL *pConnection,int *maxPlateNumber);
void WriteStarbaseRecord(PFILESTARIMAGE pFileStarImage,PPHOTPLATES pPhotPlates,FILE *dbHandle,char *title, int* pWriteHeader,int catalogNumber);
void ChangeGroup(char *filename);
void LocateNoneImages(PGSCBIN pGscBin,PFILECOMMON pFileCommon,int gsc_bin_index,PPHOTSTARIMAGE* ppMagnitudeTable,int *pAllocMagnitudes,int *pCurMagnitudes,char *catalogString,int readOnly,int verbose);
int ProcessLightcurveParameters(FILE *consoleHandle,PFILECOMMON pFileCommon,MYSQL *pPhotConnection,PPHOTSTARIMAGE pMagnitudeTable,int numMagnitudes,PSTARENTRY pCurStarEntry,PPARAMETERVECTORSTORE pParametervectorstore,int minGoodStars,int versionId,FILE *summary_handle,int *summaryCount,int debugMode,int *staleVersionIdCount,int *damonSeriesId,int *excludeSeriesId,int excludeSeriesCount,PPHOTTARGET pTarget,int gsc_bin_index);
void FreeParameterVectorStore(PPARAMETERVECTORSTORE pVectorStore);
void InitCatalogAccess(PFILECOMMON pFileCommon);

void PlotTransientCandidates(PFILECOMMON pFileCommon,PTRANCOMMON pTranCommon,PPARAMETERVECTORSTORE pVectorStore,PPHOTSTARIMAGE pMagnitudeTable,PGALAXYCOMMON pGalaxyCommon,PSTARENTRY pCurStarEntry,int refType,char *nearbyObjects,int AFLAGSMASK1,int QUALITYMASK,PMALMQUIST pMalmquist1,PMALMQUIST pMalmquist3);
void AllocFileCommonVectors(PFILECOMMON pFileCommon,int vectorSize);
void GetFullQuality(PFILESTARIMAGE pFileStarImage,int oldquality,int *pQuality);
void PlotSymbolKey();
void CombineMultipleTables(PPHOTGLOBAL pPhotGlobal,
                           PPHOTSTARIMAGE *ppMagnitudeTable,
                           PFILESTARIMAGE *ppFileStarImageTable,
                           PPHOTSTARIMAGE pMagnitudeTable1,
                           PPHOTSTARIMAGE pMagnitudeTable2,
                           PSUBSTITUTION* pSubstitutionTable,
                           int numMagnitudes1,
                           int numMagnitudes2,
                           int catalogNumber1,
                           int catalogNumber2,
                           int* pAllocMagnitudes,
                           int* pNumMagnitudes,
                           int combineAlgorithm,
                           int gsc_bin_index);

#endif /* _PHOTOMETRY_H */
