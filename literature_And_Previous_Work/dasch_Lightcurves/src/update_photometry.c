// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* Populate the photometry database.
 *
 * Inputs:
 *
 * - $DASCH_BINOUTPUT/${series}${id}_${scannum}_01${rot}ww_a$N.grid
 * - $DASCH_CATALOGALL/${series}${id}_${scannum}_01${rot}ww_extinction.db
 * - $DASCH_INGEST/${series}${id}_${scannum}_01${rot}ww.out.spatial_bins.db
 * - $DASCH_INGEST/${series}${id}_${scannum}_01${rot}ww_magdepcalibrate.db
 * - $DASCH_INGEST/${series}${id}_${scannum}_01${rot}ww_dmagcor.db
 * - $DASCH_INGEST/${series}${id}_${scannum}_01${rot}ww_allobjects.db
 * - plate list file in `-l` argument; entries are of the form `a15529_00_01r180ww`
 *
 * Outputs:
 *
 * - One of:
 *   - $DASCH_PHOT_ROOT/{refcat}/{nnn}/{nnnnnnnnnn}.dat  (new structure)
 *   - $DASCH_PHOT_ROOT/{refcat}/magNNN/magNNN/magNNN/magNNNNNNNNN.dat  (old structure)
 * - Logfile associated with `-o` argument
 *
 * Database updates:
 *
 * - New rows in `photometry.calibration{key}`
 * - New rows in `photometry.localbin{key}`
 * - New rows in `photometry.magdepcalibrate{key}`
 * - Update `photometry.photglobal` timestamp
 * - New rows in `photmetry.photinsert` for each plate
 * - New rows in `photmetry.photplates` for each plate
 *
 * Environment variables:
 *
 * - DASCH_BINOUTPUT
 * - DASCH_CATALOGALL
 * - DASCH_GID
 * - DASCH_INGEST
 * - DASCH_PHOT_ROOT
 * - DASCH_SCRATCH
 */

#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>
#include <errno.h>
#include <time.h>

#include <table.h>

#include <mysql.h>
#include <errmsg.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include "magdeputil.h"
#include "daschunistd.h"


/* Make GCC happy about potential overflows */
#define BUFPRINTF(BUF, FMT, ARGS...) \
  do { \
    int rv = snprintf((BUF), sizeof(BUF), (FMT), ARGS); \
    if (rv < 0 || rv >= sizeof(BUF)) \
      abort(); \
  } while (0)

#define MAX_FILENAME 256
#define MAX_BUFFER 256
#define MAX_LIST_STRING 25
#define MAX_QUERY_COUNT 200
#define MAGNITUDES_PLATES 10 /* Dump the magnitudes every 10 plates */
#define MAGNITUDES_LIMIT 5000000 /* Dump the magnitudes whever this limit is reached */
#define GRID_STATE_INIT 0
#define GRID_STATE_BEGIN_INTERPOLATE 1
#define GRID_STATE_NORMAL 2
#define GRID_STATE_END_INTERPOLATE 3


extern GSCBIN gscBin64;


static int invokeCmd(char *cmdStr,int *timeSpent);

static char *globalQuery1String;
static size_t globalQuery1StringSize = 0;
static size_t globalQuery1StringInitial = 0;
static size_t globalQuery1StringAlloc = 0;
static int curQuery1Count = 0;

static char *globalQuery2String;
static size_t globalQuery2StringSize = 0;
static size_t globalQuery2StringInitial = 0;
static size_t globalQuery2StringAlloc = 0;
static int curQuery2Count = 0;

static PGSCBIN pGscBin = &gscBin64;

static int magfiles_use_new_layout = 0;

typedef struct _plateEntry {
  char series[MAX_SERIES_STRING];
  char listName[MAX_LIST_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int fileExists;
  int mosaicWidth;  /* Mosaic width in pixels */
  int mosaicHeight; /* Mosaic height in pixels */
  double scale;     /* Plate scale in arcsec/pixel */
  int exposures;    /* Number of exposures */
  int quality;      /* Plate quality */
} PLATEENTRY,*PPLATEENTRY;

typedef struct _qualityEntry {
  char series[MAX_SERIES_STRING];
  int seriesId;
  int plateNumber;
  int quality;
} QUALITYENTRY,*PQUALITYENTRY;

typedef struct _dmagcor {
  int nx;             /* Total bins in width */
  int ny;             /* Total bins in height */
  int ix;             /* width bin */
  int iy;             /* height bin */
  double zout;        /* Magnitude correction */
  double errout;      /* RMS of magnitude correction */
  int npout;          /* Number of points used for correction */
  int rejectFlag;     /* Reason, if any, for rejecting this bin */
  double bright_magcor_local; /* magcor_local for brightest stars */
  double bright_magcor_error;      /* magcal_local_error for brightest stars */
  int bright_npoints_local;   /* npoints_local for brightest stars */
  /* First pass astrometry check */
  int drad_bin_count;   /* Stars in bin */
  int drad_bin_size;    /* Number of local smoothing bins used */
  int drad_reject_count; /* Stars rejected for high drad */
  double draMedian;     /* Median right ascension error in arcsec*/
  double draRMS;        /* RMS of draMedian */
  double ddecMedian;    /* Median declination error in arcsec */
  double ddecRMS;       /* RMS of ddecRMS */
  /* Second pass astrometry check */
  int drad_bin_count2;   /* Stars in bin */
  int drad_bin_size2;    /* Number of local smoothing bins used */
  int drad_reject_count2; /* Stars rejected for high drad */
  double draMedian2;     /* Median right ascension error in arcsec*/
  double draRMS2;        /* RMS of draMedian */
  double ddecMedian2;    /* Median declination error in arcsec */
  double ddecRMS2;       /* RMS of ddecRMS */
  double dradRMS2;       /* Drad RMS of the bin, after applying corrections */
} DMAGCOR,*PDMAGCOR;


typedef struct _extinction {
  int enx;             /* Total bins in width */
  int eny;             /* Total bins in height */
  int eix;             /* width bin */
  int eiy;             /* height bin */
  double altitude;    /* altitude in degrees */
  double extinction;  /* extinction in magnitudes */
} EXTINCTION,*PEXTINCTION;


typedef struct _magfilecache {
  sem_t *fileSemaphore;
  int enabled;
  int curSequence;
  int cacheHits;
  int cacheMisses;
  int fileOpens;
  int directoryCreates;
  int cur_base_index;
  PFILESTARIMAGE cacheStarImages;
  int cacheStarCount;
  int cacheStarAlloc;
  int catalogNumber;
} MAGFILECACHE, *PMAGFILECACHE;


typedef struct _updateCommon {
  char stars_name[MAX_BUFFER];
  int verbose;
  long long curUsec;
  long long deltaTable[MAX_CHARGE_ENTRY];
  long long totalTable[MAX_CHARGE_ENTRY];
  char catalogString[MAX_BUFFER];
  MYSQL *pConnection;
  MYSQL *pPhotConnection;
  time_t curTime;
  time_t keepaliveTime;
  time_t maxKeepalive;
  time_t maxInsertTime;
  time_t plateStartTime;

  char *photfilebase;
  time_t startTime;
  FILE *outHandle;
  int plateIndex;
  int curPlateCount;
  int magnitudeCount;
  int starCount;
  int overlapCount;
  int badGscBinIndexCount;
  int rejectedBlendNomatch;
  PMAGFILECACHE pMagFileCache;
} UPDATECOMMON, *PUPDATECOMMON;


/* Write remainder of multiple mysql query string. This routine helps
   consolidate multiple mysql requests into one for an improvement in speed */
static void
CompleteQuery(
  MYSQL *pPhotConnection,
  char *globalQueryString,
  int globalQueryStringSize,
  int *pCurQueryCount
) {
  int curQueryCount = *pCurQueryCount;
  int res;

  if (curQueryCount > 0) {
    globalQueryString[globalQueryStringSize-1] = 0;
    strcat(globalQueryString,";\n");
    res = ExecuteQuery(pPhotConnection,globalQueryString);
    if (res) {
      exit(1);
    }
    curQueryCount = 0;
  }
  *pCurQueryCount = curQueryCount;
}


/* Append a new query string to an in-progress query.  This routine helps
   consolidate multiple mysql requests into one for an improvement in speed */
static void
AddToQuery(
  MYSQL *pPhotConnection,
  char *queryString,
  char **pGlobalQueryString,
  size_t *pGlobalQueryStringSize,
  size_t *pGlobalQueryStringAlloc,
  int *pCurQueryCount,
  int globalQueryStringInitial,
  int debugSequence,
  int maxQueryCount
) {
  size_t globalQueryStringSize = *pGlobalQueryStringSize;
  size_t globalQueryStringAlloc = *pGlobalQueryStringAlloc;
  int curQueryCount = *pCurQueryCount;
  char *globalQueryString = *pGlobalQueryString;
  char *tmpGlobalQueryString;
  size_t querySize;

  querySize = strlen(queryString);
  if (querySize > MAX_QUERY_STRING) {
    printf("ERROR: querySize (%zu) exceeded MAX_QUERY_STRING (%d)\n", querySize, MAX_QUERY_STRING);
    exit(1);
  }

  if (globalQueryStringSize + querySize > globalQueryStringAlloc - 5) {
    // Need to increase the size of our query string
    tmpGlobalQueryString = realloc(globalQueryString, globalQueryStringSize + 2 * querySize);
    if (tmpGlobalQueryString == NULL) {
      printf(
        "ERROR: Failed to reallocate tmpGlobalQueryString of size %zu\n",
        globalQueryStringSize + 2 * querySize
      );
      exit(1);
    }

    globalQueryString = tmpGlobalQueryString;
    globalQueryStringAlloc = globalQueryStringSize + 2 * querySize;
    tmpGlobalQueryString = NULL;
  }

  strcat(globalQueryString, queryString);
  globalQueryStringSize += querySize;

  if (++curQueryCount > maxQueryCount) {
    CompleteQuery(
      pPhotConnection,
      globalQueryString,
      globalQueryStringSize,
      &curQueryCount
    );

    globalQueryString[globalQueryStringInitial] = '\0';
    globalQueryStringSize = globalQueryStringInitial;
    curQueryCount = 0;
  }

  *pGlobalQueryStringSize = globalQueryStringSize;
  *pGlobalQueryStringAlloc = globalQueryStringAlloc;
  *pCurQueryCount = curQueryCount;
  *pGlobalQueryString = globalQueryString;
}


static void
FlushCache(PMAGFILECACHE pMagFileCache, char *photfilebase)
{
  size_t written;
  char filename[MAX_FILENAME];
  int fd, i;

  if (pMagFileCache->cacheStarCount == 0) {
    return;
  }

  pMagFileCache->fileOpens++;

  if (magfiles_use_new_layout) {
    // cur_base_index is a multiple of 1024, so to hash effectively:
    int key = (pMagFileCache->cur_base_index >> 10) % 400;
    sprintf(filename, "%s/%03d/%010d.dat", photfilebase, key, pMagFileCache->cur_base_index);

    fd = open(filename, O_CREAT|O_RDWR|O_APPEND, 0664);

    if (fd < 0) {
      if (errno == ENOENT) {
        // Probably need to do this:
        char dirname[MAX_FILENAME];

        sprintf(dirname, "%s/%03d", photfilebase, key);

        if (mkdir(dirname, 0775 | S_ISGID)) {
          printf("ERROR: failed to create directory %s\n", dirname);
          exit(1);
        }

        pMagFileCache->directoryCreates++;

        fd = open(filename, O_CREAT|O_RDWR|O_APPEND, 0664);
        if (fd < 0) {
          printf("ERROR: Exiting (1) with errno %d %s %s\n", errno, strerror(errno), filename);
          exit(1);
        }
      } else {
        printf("ERROR: Exiting (2) with errno %d %s %s\n", errno, strerror(errno), filename);
        exit(1);
      }
    }
  } else {
    int dir1 = pMagFileCache->cur_base_index / (MAG_SUBDIR_MODULUS * MAG_SUBDIR_MODULUS);
    int dir2 = pMagFileCache->cur_base_index / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);

    sprintf(filename, "%s/mag%03d/mag%03d/mag%09d.dat", photfilebase, dir1, dir2, pMagFileCache->cur_base_index);

    fd = open(filename, O_CREAT|O_RDWR|O_APPEND, S_IRWXU|S_IRWXG|S_IROTH);

    if (fd < 0) {
      if (errno == ENOENT) {
        char chgrpname[MAX_FILENAME];
        char cmdstr[MAX_FILENAME];
        int result;
        int timeSpent;

        /* We need to create the directory */
        sprintf(cmdstr, "mkdir -p %s/mag%03d/mag%03d", photfilebase, dir1, dir2);
        result = invokeCmd(cmdstr, &timeSpent);

        if (result != 0) {
          printf("ERROR: failed to create directory %s\n", cmdstr);
          exit(1);
        }

        sprintf(chgrpname, "%s/mag%03d", photfilebase, dir1);
        ChangeGroup(chgrpname);
        sprintf(chgrpname, "%s/mag%03d/mag%03d", photfilebase, dir1, dir2);
        ChangeGroup(chgrpname);

        pMagFileCache->directoryCreates++;
        fd = open(filename, O_CREAT|O_RDWR|O_APPEND, S_IRWXU|S_IRWXG|S_IROTH);
        if (fd < 0) {
          printf("ERROR: Exiting (1) with result %d errno %d %s %s\n", fd, errno, strerror(errno), filename);
          exit(1);
        }
      } else {
        printf("ERROR: Exiting (2) with result %d errno %d %s %s\n", fd, errno, strerror(errno), filename);
        exit(1);
      }
    }
  }

  // Write the data item-by-item so that they will atomically append to the file
  // without stomping on other potential writes. This is only guaranteed for
  // writes of size less than PIPE_BUF.

  for (i = 0; i < pMagFileCache->cacheStarCount; i++) {
    written = write(fd, &pMagFileCache->cacheStarImages[i], sizeof(FILESTARIMAGE));

    if (written != sizeof(FILESTARIMAGE)) {
      printf("Wrote only %zu bytes to for gsc_base_index %d\n", written, pMagFileCache->cur_base_index);
      break;
    }
  }

  close(fd);

  if (!magfiles_use_new_layout) {
    ChangeGroup(filename);
  }

  pMagFileCache->cacheStarCount = 0;
}


static int
invokeCmd(char *cmdStr,int *timeSpent) {
  int result;
  time_t beginTime;
  time_t curTime;
  time(&beginTime);
  *timeSpent = 0;
  result = system(cmdStr);
  time(&curTime);
  beginTime = curTime - beginTime;
  *timeSpent = (int)beginTime;
  if (WIFSIGNALED(result) &&
      (WTERMSIG(result) == SIGINT || WTERMSIG(result) == SIGQUIT)) {
    printf("Terminated with signal %d at %d seconds\n",WTERMSIG(result),(int)curTime);
    exit(1);
  }

  return result;
}


#define STAR_CACHE_INCREMENT 1024

static int
StoreMagnitudeFileEntry(
  PMAGFILECACHE pMagFileCache,
  PFILESTARIMAGE pFileStarImage,
  char *photfilebase
) {
  int gsc_base_index = pFileStarImage->gsc_bin_index - (pFileStarImage->gsc_bin_index % MAG_FILE_MODULUS);
  PFILESTARIMAGE pCacheStarImage;

  pFileStarImage->versionTag = MAG_FORMAT_VERSION;
  pFileStarImage->versionTag = pFileStarImage->versionTag << 32;
  pFileStarImage->versionTag = pFileStarImage->versionTag | MAG_FORMAT_DATA;

  if (pMagFileCache->enabled == 0) {
    pMagFileCache->curSequence = 0;
    pMagFileCache->enabled = 1;
    pMagFileCache->cacheHits = 0;
    pMagFileCache->cacheMisses = 0;
    pMagFileCache->fileOpens = 0;
    pMagFileCache->cur_base_index = gsc_base_index;
    pMagFileCache->cacheStarImages = (PFILESTARIMAGE) calloc(STAR_CACHE_INCREMENT, sizeof(FILESTARIMAGE));
    if (pMagFileCache->cacheStarImages == NULL) {
      printf("ERROR: failed on first allocation of cacheStarImages\n");
      exit(1);
    }
    pMagFileCache->cacheStarCount = 0;
    pMagFileCache->cacheStarAlloc = STAR_CACHE_INCREMENT;
  }

  if (pMagFileCache->cur_base_index != gsc_base_index) {
    pMagFileCache->cacheMisses++;
    FlushCache(pMagFileCache, photfilebase);
  } else {
    pMagFileCache->cacheHits++;
  }

  pMagFileCache->cur_base_index = gsc_base_index;
  pMagFileCache->curSequence++;

  if (pMagFileCache->cacheStarCount >= pMagFileCache->cacheStarAlloc) {
    // Need a larger cache
    pMagFileCache->cacheStarAlloc += STAR_CACHE_INCREMENT;

    pMagFileCache->cacheStarImages = realloc(
      pMagFileCache->cacheStarImages,
      pMagFileCache->cacheStarAlloc * sizeof(FILESTARIMAGE)
    );

    if (pMagFileCache->cacheStarImages == NULL) {
      printf(
        "ERROR: failed to reallocate cacheStarImages of size %d %zu\n",
        pMagFileCache->cacheStarAlloc,
        pMagFileCache->cacheStarAlloc * sizeof(FILESTARIMAGE)
      );
      exit(1);
    }
  }

  pCacheStarImage = &pMagFileCache->cacheStarImages[pMagFileCache->cacheStarCount];
  memcpy(pCacheStarImage, pFileStarImage, sizeof(FILESTARIMAGE));
  pMagFileCache->cacheStarCount++;
  return 1;
}


typedef struct _allobjectsrecord {
  char series[MAX_SERIES_STRING]; /* Series */
  int seriesId;
  int plateNumber;
  int exposureNumber;
  int solutionNumber;
  int local_bin_index;
  int magdep_bin;             /* Magnituded-dependent calibration bin (-1) if undefined */
  char REF[MAX_REF];          /* GSC2.3.2 reference number */
  long long REFNumber;        /* Translated reference number */
  int NUMBER;                 /* Sextractor reference number */
  int passBits;               /* Pass identifier */
  double X_IMAGE;
  double Y_IMAGE;
  int AFLAGS;                 /* Flags keyword. See header of pipelineutils.h */
  int BFLAGS;                 /* Flags keyword. See header of pipelineutils.h */
  int npoints_local;
  int spatial_bin;            /* Lowess magnitude spatial bin  number */
  double MAG_ISO;
  double ra;                  /* Right Ascension in degrees */
  double dec;                 /* Declination in degrees */
  double magcal_iso;          /* lowess magnitude estimate */
  double magcal_iso_rms;      /* local  error */
  double magcal_local;        /* Local magnitude calibration */
  double magcal_local_error;
  double magcal_local_rms;    /* Overall error */
  double magcal_magdep;  /* magnitude-corrected local magnitude */
  double magcal_magdep_rms; /* RMS of the magnitude-corrected local magnitude */
  double Date;                /* Heliocentric Julian Date  */
  double limiting_mag_local;  /* Limiting magnitude */
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
  int gsc_bin_index;
  int versionId;
  /* The following fields come from other tables and are used in find_lightcurves.c */
  double colorterm;
  double errorcolor;
  int colorflag;
  int AFLAGSCOPY;
  int rejectFlag;
  double dradRMS2;
  double magcor_local;       /* zout */
  char Plate[MAX_PLATE_NAME];
} ALLOBJECTSRECORD,*PALLOBJECTSRECORD;


/* Sort routine based on the GSC bin index */
static int
AllobjectsCompare(const void *first, const void *second)
{
  int numberFirst = ((PALLOBJECTSRECORD)first)->gsc_bin_index;
  int numberSecond = ((PALLOBJECTSRECORD)second)->gsc_bin_index;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  } else {
    return 0;
  }
}


static void
WriteFinalStatistics(PUPDATECOMMON pUpdateCommon)
{
  int usecIndex;
  time(&pUpdateCommon->curTime);
  pUpdateCommon->curTime -= pUpdateCommon->startTime;

  printf(
    "Summary: badGscBinIndex=%d, BlendNomatch=%d maxKeepalive=%ld maxInsertTime=%ld\n",
    pUpdateCommon->badGscBinIndexCount,
    pUpdateCommon->rejectedBlendNomatch,
    pUpdateCommon->maxKeepalive,
    pUpdateCommon->maxInsertTime
  );

  ChargeTime(0,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
  printf("Time Breakdown (sec) ");
  fprintf(pUpdateCommon->outHandle,"Time Breakdown (sec) ");
  for (usecIndex = 0; usecIndex < MAX_CHARGE_ENTRY; usecIndex++) {
    double outSeconds;
    outSeconds = (1.0 * pUpdateCommon->totalTable[usecIndex])/1000000.0;
    printf("%10.3f ",outSeconds);
    fprintf(pUpdateCommon->outHandle,"%10.3f ",outSeconds);
  }
  printf("\n");
  fprintf(pUpdateCommon->outHandle,"\n");

  FlushCache(pUpdateCommon->pMagFileCache, pUpdateCommon->photfilebase);
  printf(
    "curSequence=%d cacheHits=%d cacheMisses=%d fileOpens=%d newDirs=%d cacheAlloc=%d\n",
    pUpdateCommon->pMagFileCache->curSequence,
    pUpdateCommon->pMagFileCache->cacheHits,
    pUpdateCommon->pMagFileCache->cacheMisses,
    pUpdateCommon->pMagFileCache->fileOpens,
    pUpdateCommon->pMagFileCache->directoryCreates,
    pUpdateCommon->pMagFileCache->cacheStarAlloc
  );
}


int
main(int argc, char *argv[])
{
  UPDATECOMMON updateCommon;
  PUPDATECOMMON pUpdateCommon = &updateCommon;
  char *argstr;
  char listname[MAX_FILENAME];
  char versionname[MAX_VERSION_NAME+1];
  FILE * listFile;
  int errorFlag = 0;
  int numPlates = 0;
  int maxString = 0;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  int lineLen;
  char* listStrings = NULL;
  int curline = 0;
  time_t curPlateTime;
  int remoteDirectory = 1;
  int useDiskLocation = 0;
  int diskLocation = -1;
  int noNewVersion = 0;
  char cmdchar;
  PPLATEENTRY pPlateEntry;
  PPLATEENTRY pPlateTable = NULL;
  char outfile[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  char *scriptDir = NULL;
  char *ingestDirectory = NULL;
  char *binoutputDirectory = NULL;
  char *catalogallDirectory = NULL;

  char ingestRemoteDirectory[MAX_BUFFER];
  char binoutputRemoteDirectory[MAX_BUFFER];
  char catalogallRemoteDirectory[MAX_BUFFER];

  char *scratchDirectory;

  int magnitudeFlag = 0;
  int magdepFlag = 0;
  int platesFlag = 0;
  int localBinFlag = 0;
  int lowessGridFlag = 0;
  int spatialBinFlag = 0;
  int magnitudeLimitCount = 0;
  int gridWorkDone;

  int irec;
  int local_bin_index;

  int decBin;
  int raBin;
  int tmp_gsc_bin_index;
  int nx;                   /* dmagcor_table width */
  int ny;                   /* dmagcor_table height */
  int ix;
  int iy;

  char timestr[100];
  struct tm *ptr;
  int gotAnswer;

  MYSQL my_connection;
  MYSQL my_phot_connection;
  int spatial_bin;

  File spatial_bins_handle = NULL;
  char spatial_bins_name[MAX_BUFFER];
  char solutionString[MAX_BUFFER];
  char spatial_bin_string[MAX_BUFFER];
  PPHOT_SPATIAL_BIN tmp_spatial_bin_table = NULL;
  TableHead spatial_bins_header = NULL;
  PPHOT_SPATIAL_BIN pSpatial_bin;
  PPHOT_SPATIAL_BIN pSpatial_bin2;
  size_t spatial_bin_nrecs;
  int index;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  PHOTPLATES basePhotPlates;
  PPHOTPLATES pPhotPlates = &basePhotPlates;

  PHOT_SPATIAL_BIN spatial_bin_table[MAX_SPATIAL_BINS+1];

  File dmagcor_handle = NULL;
  char dmagcor_name[MAX_BUFFER];
  TableHead dmagcor_header = NULL;
  size_t dmagcor_nrecs;
  PDMAGCOR dmagcor_table[MAX_CATALOG_EXPOSURES];
  PDMAGCOR pDmagcor;

  File extinction_handle = NULL;
  char extinction_name[MAX_BUFFER];
  char grid_name[MAX_BUFFER];
  FILE *grid_file = NULL;
  TableHead extinction_header = NULL;
  PEXTINCTION extinction_table[MAX_CATALOG_EXPOSURES];
  double timeAccuracy[MAX_CATALOG_EXPOSURES];
  size_t extinction_nrecs = 0;
  PEXTINCTION pExtinction;

  int nvals;
  int gridCounter;
  int gridState;
  float maggrid[MAX_GRID_ENTRY];
  float isogrid[MAX_GRID_ENTRY];
  float griderr[MAX_GRID_ENTRY];
  int flaggrid[MAX_GRID_ENTRY];
  float lastmaggrid;
  float lastisogrid;
  float lastgriderr;
  int lastflaggrid;

  File allobjects_handle = NULL;
  char allobjects_name[MAX_BUFFER];
  TableHead allobjects_header = NULL;
  PALLOBJECTSRECORD allobjects_table = NULL;
  size_t allobjects_nrecs = 0;
  int allobjects_index;
  PALLOBJECTSRECORD pAllobjects = NULL;

  double THRESHOLD;

  int new_version = 0; /* If true, create a new version */
  int res;
  char queryString[MAX_QUERY_STRING];

  int plateMagnitudeCount = 0;
  int plateQueryCount = 0;

  MOSAIC mosaic;
  PMOSAIC pMosaic = &mosaic;
  int fullSucceededFlag = 0;
  int lastFitWCS = 0;

  FILESTARIMAGE fileStarImage;
  PFILESTARIMAGE pFileStarImage = &fileStarImage;
  MAGFILECACHE magFileCache;
  int columnIndex;
  int gotExposureNumber = 0;
  int gotSolutionNumber = 0;
  int solutionNumber;
  int maxSolutionNumber;
  int printMaxSolutionError = 1;

  /* The following four variables are used for ChargeTime()
   * Index 0: non database work
   *       1: update localbin
   *       2: update magnitudes
   *       3: update stars
   *       4: update photplates
   *       5: update calibration
   *       6: update magnitude dependent calibration
   */

  char solutionNumberString[20];
  char dasch_phot_magnitudes[40];
  int refType;
  int binPrintFlag = 1;

  double singleJulianDate;
  double singleTolerance;
  double allJulianDate;
  double allTolerance;
  int numExposures = 0;
  int exposureNumber;
  PMAGDEPCORRECTION pMagdepTable = NULL;
  PMAGDEPCORRECTION pMagdep;
  MAGDEPLIMITS magdepLimits;
  PMAGDEPLIMITS pMagdepLimits = &magdepLimits;
  char magdep_name[MAX_BUFFER];
  int magdep_load_status;
  int magdep_index;
  int sanityErrorCount = 0;
  int timeAccuracyError = 0;
  int seriesId = -1;

  memset(pUpdateCommon, 0, sizeof(UPDATECOMMON));
  pUpdateCommon->pConnection = &my_connection;
  pUpdateCommon->pPhotConnection = &my_phot_connection;
  pUpdateCommon->pMagFileCache = &magFileCache;

  for (solutionNumber = 0; solutionNumber < MAX_CATALOG_EXPOSURES; solutionNumber++) {
    extinction_table[solutionNumber] = NULL;
    dmagcor_table[solutionNumber] = NULL;
    timeAccuracy[solutionNumber] = -1;
  }

  memset(pMagdepLimits, 0, sizeof(MAGDEPLIMITS));

  memset(pUpdateCommon->pMagFileCache, 0, sizeof(MAGFILECACHE));

  pUpdateCommon->catalogString[0] = 0;
  ChargeTime(RESET_ALL_ENTRY, &pUpdateCommon->curUsec, pUpdateCommon->deltaTable, pUpdateCommon->totalTable);

  SetQueryCount(0);
  globalQuery1String = (char *) malloc(MAX_QUERY_STRING);
  if (globalQuery1String == NULL) {
    printf("Error allocating globalQuery1String\n");
    exit(1);
  }
  globalQuery1StringAlloc = MAX_QUERY_STRING;

  globalQuery2String = (char *) malloc(MAX_QUERY_STRING);
  if (globalQuery2String == NULL) {
    printf("Error allocating globalQuery2String\n");
    exit(1);
  }
  globalQuery2StringAlloc = MAX_QUERY_STRING;

  listname[0] = 0;
  outfile[0] = 0;
  qualifier[0] = 0;
  versionname[0] = 0;

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
        case 'l': /* list file name */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for list file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(listname,*++argv,MAX_BUFFER-2);
            if (strlen(listname) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(outfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            pUpdateCommon->pMagFileCache->catalogNumber = GetCatalogNumber(*++argv);

            if (pUpdateCommon->pMagFileCache->catalogNumber < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              sprintf(pUpdateCommon->catalogString,"%d",pUpdateCommon->pMagFileCache->catalogNumber);
              sprintf(qualifier,"_%s",catalogText[pUpdateCommon->pMagFileCache->catalogNumber]);
            }
          }
          break;

        case 'n': /* new version */
        case 'N':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(versionname,*++argv,MAX_VERSION_NAME-2);
            if (strlen(versionname) >= MAX_VERSION_NAME-3) {
              printf("ERROR: MAX_VERSION NAME %d exceeded %zu\n",MAX_VERSION_NAME,strlen(versionname));
              errorFlag = 1;
            } else {
              new_version = 1;
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          pUpdateCommon->verbose = 1;
          break;

        case 'u': /* No new version */
        case 'U':
          noNewVersion = 1;
          break;

        case 'a': /* All Flags */
        case 'A':
          magnitudeFlag = 1;
          platesFlag = 1;
          localBinFlag = 1;
          spatialBinFlag = 1;
          lowessGridFlag = 1;
          magdepFlag = 1;
          break;

        case 'f': /* magnitude-dependent calibration */
        case 'F':
          magdepFlag = 1;
          break;

        case 'm': /* magnitude flag */
        case 'M':
          magnitudeFlag = 1;
          break;

        case 'p': /* plates flag */
        case 'P':
          platesFlag = 1;
          break;

        case 'c': /* local bin flag */
        case 'C':
          localBinFlag = 1;
          break;

        case 'd': /* lowess calibration flag */
          lowessGridFlag = 1;
          break;

        case 'b': /* spatial bin flag */
        case 'B':
          spatialBinFlag = 1;
          break;

        case 'g': /* Synchronize quality columns */
        case 'G':
          fprintf(stderr, "fatal: -g \"synchronize quality column\" no longer implemented\n");
          return 1;

        case 'r': /* Use local directory instead of remote */
        case 'R':
          remoteDirectory = 0;
          break;

        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&diskLocation);
            if (nvals != 1) {
              printf("ERROR: Unable to decode disk location %s\n",*argv);
              errorFlag = 1;
            } else {
              useDiskLocation = 1;
            }
          }

          break;

        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  if (useDiskLocation == 1 && remoteDirectory == 0) {
    printf("ERROR: using both -D and -r is invalid\n");
    errorFlag = 1;
  }

  scratchDirectory = getenv("DASCH_SCRATCH");
  if (scratchDirectory == NULL) {
    printf("ERROR: DASCH_SCRATCH is not defined\n");
    return 1;
  }

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pUpdateCommon->pConnection);
  dasch_init_photdb(pUpdateCommon->pPhotConnection);
  magfiles_use_new_layout = dasch_photdb_magfiles_use_new_layout();

  InitSeriesTable(pUpdateCommon->pConnection, pUpdateCommon->pPhotConnection);

  if (GetPhotometryGlobal(pUpdateCommon->pPhotConnection, pPhotGlobal) != 1) {
    printf("ERROR: failed to get the global photometry table\n");
    exit(1);
  }

  if (pPhotGlobal->magnitudeFile != PHOT_MAGNITUDEFILE_YES) {
    fprintf(stderr, "fatal: magnitudeFile must be true\n");
    exit(1);
  }

  if (pPhotGlobal->solutionNumber == PHOT_SOLUTIONNUMBER_YES) {
    strcpy(solutionNumberString, "solutionNumber");
  } else {
    strcpy(solutionNumberString, "exposureNumber");
  }

  /* Validate arguments */

  pUpdateCommon->outHandle = fopen(outfile,"a+t");
  if (pUpdateCommon->outHandle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to open the output file %s\n", outfile);
  } else if (pUpdateCommon->verbose) {
    printf("Output file %s\n",outfile);
  }

  if (new_version == 0 || listname[0] != '\0') {
    if (outfile[0] == '\0') {
      printf("ERROR: No output filename was specified\n");
      errorFlag = 1;
    }

    /* Attempt to open the list */

    if (listname[0] == 0) {
      printf("ERROR: No plate list specified\n");
      errorFlag = 1;
    }

    if (listname[0] != 0) {
      listFile = fopen(listname,"rt");
      if (listFile == NULL) {
        printf("ERROR: Could not open file %s\n",listname);
        errorFlag = 1;
      } else {
        /* Count the number of records in the file */
        while (1) {
          inBuffer = fgets(inLine,MAX_BUFFER,listFile);
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

          if (lineLen > 5) {
            numPlates++;
            if (maxString < lineLen) {
              maxString = lineLen;
            }
          }
        }

        fclose(listFile);

        if (numPlates == 0) {
          printf("ERROR: no lines in the list file\n");
          errorFlag = 1;
        }
      }
    }
  }

  if (errorFlag) {
    printf("Usage: update_photometry options\n");
    printf("  options: -v verbose\n");
    printf("           -l <Plate list file>\n");
    printf("           -o <output file>\n");
    printf("           -n create a new version \n");
    printf("           -a write all tables \n");
    printf("           -m write magnitudes table\n");
    printf("           -p write plates table\n");
    printf("           -c write local bin table\n");
    printf("           -f write magnitude-dependent calibration table\n");
    printf("           -b write spatial bin table\n");
    printf("           -d write lowess calibration tables\n");
    printf("           -q input catalog and filename qualifer\n");
    printf("           -u disable automatic version (does not override -n)\n");
    printf("           -g synchronize the quality column with that of the plates database\n");
    printf("           -r read files locally rather than from the system with the scripts directory\n");
    printf("           Starbase files are on DASCH_BINOUTPUT, DASCH_CATALOGBIN, and  DASCH_INGEST unless: \n");
    printf("           -D nnn use starbase files on /dasch/raidnnn/Pipeline\n");
    printf("           -r use starbase files on DASCH_SCRIPTS\n");
    return 1;
  }

  BUFPRINTF(dasch_phot_magnitudes, "DASCH_PHOT_MAGNITUDES%s", pUpdateCommon->catalogString);

  pUpdateCommon->photfilebase = getenv(dasch_phot_magnitudes);
  if (pUpdateCommon->photfilebase == NULL) {
    printf("%s is not defined\n",dasch_phot_magnitudes);
    return 1;
  }

  if (useDiskLocation != 0) {
    sprintf(ingestRemoteDirectory,"/dasch/raid%03d/Pipeline/ingest",diskLocation);
    ingestDirectory = ingestRemoteDirectory;

    sprintf(binoutputRemoteDirectory,"/dasch/raid%03d/Pipeline/bin9",diskLocation);
    binoutputDirectory = binoutputRemoteDirectory;

    sprintf(catalogallRemoteDirectory,"/dasch/raid%03d/Pipeline/catalogall",diskLocation);
    catalogallDirectory = catalogallRemoteDirectory;
  } else if (remoteDirectory == 0) {
    catalogallDirectory = getenv("DASCH_CATALOGALL");
    if (catalogallDirectory == NULL) {
      printf("ERROR: DASCH_CATALOGALL is not defined\n");
      return 1;
    }

    ingestDirectory = getenv("DASCH_INGEST");
    if (ingestDirectory == NULL) {
      printf("ERROR: DASCH_INGEST is not defined\n");
      return 1;
    }

    binoutputDirectory = getenv("DASCH_BINOUTPUT");
    if (binoutputDirectory == NULL) {
      printf("ERROR: DASCH_BINOUTPUT is not defined\n");
      return 1;
    }
  } else {
    scriptDir = getenv("DASCH_SCRIPTS");
    if (scriptDir == NULL) {
      printf("ERROR: DASCH_SCRIPTS is not defined\n");
      return 0;
    }

    sprintf(ingestRemoteDirectory,"%s/ingest",scriptDir);
    ingestDirectory = ingestRemoteDirectory;

    sprintf(binoutputRemoteDirectory,"%s/bin9",scriptDir);
    binoutputDirectory = binoutputRemoteDirectory;

    sprintf(catalogallRemoteDirectory,"%s/catalogall",scriptDir);
    catalogallDirectory = catalogallRemoteDirectory;
  }

  printf(
    "update_photometry of %s %s, list file %s entries %d, max string %d catalog %s filebase %s\n",
    __DATE__,
    __TIME__,
    listname,
    numPlates,
    maxString,
    catalogText[pUpdateCommon->pMagFileCache->catalogNumber],
    pUpdateCommon->photfilebase
  );
  printf(
    "Reading Starbase files from %s %s and %s\n",
    binoutputDirectory,
    catalogallDirectory,
    ingestDirectory
  );

  /* The following relation is used for MAG_FORMAT_DATA and  MAG_FORMAT_VERSION tags */
  if (sizeof(int) * 2 > sizeof(long long)) {
    printf("ERROR: int size %zu not half of long long size %zu\n",sizeof(int),sizeof(long long));
    exit(1);
  }

  if (UpdatePhotometryGlobal(pUpdateCommon->pPhotConnection, pPhotGlobal, new_version, versionname) != 1) {
    return 1;
  }

  printf(
    "update_photometry is now at version %d %s timestamp %s\n",
    pPhotGlobal->currentVersion,
    pPhotGlobal->versionName,
    pPhotGlobal->timeStamp
  );

  TimeStampPhotometryGlobal(pUpdateCommon->pPhotConnection);

  if (new_version == 1 && listname[0] == '\0') {
    fclose(pUpdateCommon->outHandle);
    mysql_close(pUpdateCommon->pPhotConnection);
    mysql_close(pUpdateCommon->pConnection);
    exit(0);
  }

  /* Now allocate space for all of the records */
  if (listname[0] != 0) {
    maxString += 5;

    listStrings = calloc((numPlates + 1) * maxString, sizeof(char));
    if (listStrings == NULL) {
      printf("Failed to allocate listStrings\n");
      return 1;
    }

    strcpy(&listStrings[numPlates * maxString], "UNKNOWN");

    pPlateTable = (PPLATEENTRY) calloc(numPlates, sizeof(PLATEENTRY));
    if (pPlateTable == NULL) {
      printf("Failed to allocate pPlateTable\n");
      return 1;
    }

    /* Attempt to open the list */
    listFile = fopen(listname,"rt");
    if (listFile == NULL) {
      printf("ERROR: Could not open file %s\n", listname);
      return 1;
    } else {
      /* Read the file records */
      while (1) {
        inBuffer = fgets(inLine,MAX_BUFFER,listFile);
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

        if (lineLen > 5) {
          strcpy(&listStrings[curline*maxString],inBuffer);
          curline++;
        }
      }

      fclose(listFile);
    }
  }

  time(&pUpdateCommon->startTime);
  time(&pUpdateCommon->keepaliveTime);

  ptr = localtime(&pUpdateCommon->startTime);
  strftime(timestr, 25, "%Y-%m-%dT%H-%M-%S", ptr);
  printf("update_photometry starting at %s\n", timestr);

  sprintf(queryString,"update photglobal set keepalive = 'yes';");
  res = ExecuteQuery(pUpdateCommon->pPhotConnection, queryString);

  // For each plate ...

  for (pUpdateCommon->plateIndex = 0; pUpdateCommon->plateIndex < numPlates; pUpdateCommon->plateIndex++) {
    time(&pUpdateCommon->plateStartTime);
    plateMagnitudeCount = 0;
    plateQueryCount = GetQueryCount();
    gotAnswer = GetPhotometryGlobal(pUpdateCommon->pPhotConnection, pPhotGlobal);

    if (gotAnswer != 1 || pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES) {
      ptr = localtime(&pUpdateCommon->plateStartTime);
      strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
      printf("update_photometry aborting at %s with plates %d  gotAnswer %d keepalive %d\n",
             timestr,pUpdateCommon->plateIndex,gotAnswer,pPhotGlobal->keepalive);
      fprintf(pUpdateCommon->outHandle,"update_photometry aborting at %s with plates %d  gotAnswer %d keepalive %d\n",
              timestr,pUpdateCommon->plateIndex,gotAnswer,pPhotGlobal->keepalive);
      break;
    }

    // Load up this plate's basic information

    pPlateEntry = &pPlateTable[pUpdateCommon->plateIndex];
    strcpy(pPlateEntry->listName,&listStrings[pUpdateCommon->plateIndex*maxString]);
    if (ParseFilename(pPlateEntry->listName,pPlateEntry->series,&pPlateEntry->plateNumber,&pPlateEntry->mosaicNumber,&pPlateEntry->binning,&pPlateEntry->rotation) == 0) {
      printf("ERROR: Failed to parse the name  %s\n",pPlateEntry->listName);
      continue;
    }

    if (strcmp(pPlateEntry->series,"na") == 0) {
      printf("ERROR: excluding 'na' series plates %s\n",pPlateEntry->listName);
      continue;
    }

    if (GetPlateScale(pUpdateCommon->pConnection,pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->mosaicNumber,0,&pPlateEntry->mosaicWidth,&pPlateEntry->mosaicHeight,&pPlateEntry->scale)) {
      printf("ERROR: Failed to get plate scale for %s\n",pPlateEntry->listName);
      continue;
    }

    nx = XDmagBins(pPlateEntry->mosaicWidth,pPlateEntry->mosaicHeight);
    ny = YDmagBins(pPlateEntry->mosaicWidth,pPlateEntry->mosaicHeight);

    /* Read the plate quality here */
    if ((SetPlateQuality(pUpdateCommon->pConnection,"multiple",pPlateEntry->listName,-1,1,&pPlateEntry->quality)) != 0) {
      exit(1);
    }

    if (GetPhotPlate(pUpdateCommon->pPhotConnection,pPlateEntry->series,pPlateEntry->plateNumber,pPhotPlates,pUpdateCommon->catalogString,0)) {
      exit(1);
    }

    pPhotPlates->quality = pPlateEntry->quality;
    pPhotPlates->mosaicNumber = pPlateEntry->mosaicNumber;

    // Do we need to increment the photometry version?

    if (
      (pPhotPlates->versionId != 0) &&
      (pPhotPlates->versionId == pPhotGlobal->currentVersion) &&
      (new_version == 0)
    ) {
      // We have a duplicate plate version with new data. Update the version.

      if (noNewVersion) {
        fprintf(stderr, "FATAL: need to update phot version from %d but we are forbidden\n", pPhotGlobal->currentVersion);
        exit(1);
      }

      new_version = 1;
      if (GetPhotometryGlobal(pUpdateCommon->pPhotConnection, pPhotGlobal) != 1) {
        printf("ERROR reading photglobal\n");
        exit(1);
      }

      strcpy(versionname, pPhotGlobal->timeStamp);

      if (UpdatePhotometryGlobal(pUpdateCommon->pPhotConnection, pPhotGlobal, new_version, versionname) != 1) {
        return 1;
      }

      printf("New version %d %s because of a duplicate plate\n", pPhotGlobal->currentVersion, versionname);
    }

    // Log the insert with a failure flag; we'll clear it if/when the insert fully succeeds.

    seriesId = GetSeriesId(pPlateEntry->series,1);
    sprintf(queryString,"INSERT IGNORE INTO photinsert (seriesId,plateNumber,catalogNumber) VALUES (%d,%d,%d);\n",
            seriesId,
            pPlateEntry->plateNumber,
            pUpdateCommon->pMagFileCache->catalogNumber);
    res = ExecuteQuery(pUpdateCommon->pPhotConnection,queryString);
    if (res != 0) {
      exit(1);
    }

    sprintf(queryString,"UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
            pPlateEntry->mosaicNumber,
            pPhotGlobal->currentVersion,
            1,
            seriesId,
            pPlateEntry->plateNumber,
            pUpdateCommon->pMagFileCache->catalogNumber);
    res = ExecuteQuery(pUpdateCommon->pPhotConnection,queryString);
    if (res != 0) {
      exit(1);
    }

    // Free buffers from previous iteration, and other cleanup

    timeAccuracyError = 0;
    for (solutionNumber = 0; solutionNumber < MAX_CATALOG_EXPOSURES; solutionNumber++) {
      if (extinction_table[solutionNumber] != NULL) {
        Free(extinction_table[solutionNumber]);
        extinction_table[solutionNumber] = NULL;
      }
      if (dmagcor_table[solutionNumber] != NULL) {
        Free(dmagcor_table[solutionNumber]);
        dmagcor_table[solutionNumber] = NULL;
      }
      timeAccuracy[solutionNumber] = -1;
    }

    if (pMagdepTable != NULL) {
      free(pMagdepTable);
      FreeMagdepSubarrays(pMagdepLimits);
      pMagdepTable = NULL;
    }

    maxSolutionNumber = 0;
    fullSucceededFlag = 0;

    // For each solution ...

    for (solutionNumber = 0; solutionNumber < MAX_CATALOG_EXPOSURES; solutionNumber++) {
      if (solutionNumber == 0) {
        solutionString[0] = 0;
      } else {
        sprintf(solutionString,"_s%d",solutionNumber);
      }

      if (tmp_spatial_bin_table != NULL) {
        Free(tmp_spatial_bin_table);
        tmp_spatial_bin_table = NULL;
      }

      if (spatial_bins_header != NULL) {
        table_hdrfree(spatial_bins_header);
        spatial_bins_header = NULL;
      }

      if (spatial_bins_handle != NULL) {
        Close(spatial_bins_handle);
        spatial_bins_handle = NULL;
      }

      if (extinction_header != NULL) {
        table_hdrfree(extinction_header);
        extinction_header = NULL;
      }

      if (extinction_handle != NULL) {
        Close(extinction_handle);
        extinction_handle = NULL;
      }

      if (dmagcor_header != NULL) {
        table_hdrfree(dmagcor_header);
        dmagcor_header = NULL;
      }

      if (dmagcor_handle != NULL) {
        Close(dmagcor_handle);
        dmagcor_handle = NULL;
      }

      // Load more metadata

      errorFlag = 0;
      numExposures = 0;

      if (GetMosaicInfo(pUpdateCommon->pConnection,pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->mosaicNumber,solutionNumber,pMosaic) != 1) {
        errorFlag = 1;
      } else {
        if (GetSolutionJulianDate(pUpdateCommon->pConnection,pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->mosaicNumber,solutionNumber,&exposureNumber,&numExposures,&singleJulianDate,&singleTolerance,&allJulianDate,&allTolerance)) {
          timeAccuracy[solutionNumber] = -1;
        } else {
          timeAccuracy[solutionNumber] = singleTolerance;
        }

        if (strstr(pMosaic->mosaicComment,"Rejected for Markings") != NULL) {
          /* Change of Feb 18, 2020 */
          printf("ERROR: Rejected for Markings plate %s\n",pPlateEntry->listName);
          errorFlag = 1;
        }

        if ((fullSucceededFlag != 0) && ((pMosaic->FitWCS & FITWCS_FULLSUCCEEDED) == 0)) {
          /* Here we have a mixture of new and old solutions. Do not go further */
          errorFlag = 1;
        }

        fullSucceededFlag = pMosaic->FitWCS & FITWCS_FULLSUCCEEDED;

        if ((solutionNumber > 0) && ((lastFitWCS & (FITWCS_MULTFAILEDASTROMETRY|FITWCS_MULTFAILEDSEPARATION)) != 0)) {
          errorFlag = 1;
        }

        lastFitWCS = pMosaic->FitWCS;
      }

      // Load "spatial bin" calibration info

      if (!errorFlag && spatialBinFlag) {
        strcpy(spatial_bins_name,ingestDirectory);
        strcat(spatial_bins_name,"/");
        strcat(spatial_bins_name,pPlateEntry->listName);
        strcat(spatial_bins_name,solutionString);
        strcat(spatial_bins_name,qualifier);
        strcat(spatial_bins_name,".out.spatial_bins.db");
        spatial_bins_handle = Open(spatial_bins_name,"r");
        if (spatial_bins_handle == NULL) {
          if (solutionNumber == 0) {
            printf("ERROR: Failed to find spatial bins file %s\n",spatial_bins_name);
          }
          errorFlag = 1;

        } else {
          if (pUpdateCommon->verbose) {
            printf("Found spatial bins file %s\n",spatial_bins_name);
          }
        }

        /* Now read in the spatial bins file */
        if (!errorFlag) {
          spatial_bins_header = table_header(spatial_bins_handle,TABLE_PARSE);
          if (spatial_bins_header == NULL) {
            printf("ERROR: Failed to read header for %s\n",spatial_bins_name);
            errorFlag = 1;

          }
        }

        if (!errorFlag) {
          tmp_spatial_bin_table = table_loadva(spatial_bins_handle,
                                               &spatial_bins_header,
                                               NULL, /* hbase */
                                               NULL, /* rows */
                                               NULL,
                                               sizeof(PHOT_SPATIAL_BIN),
                                               &spatial_bin_nrecs,
                                               TblInt,"spatial_bin",TblOff(PPHOT_SPATIAL_BIN,spatial_bin),
                                               TblDbl,"limiting_mag",TblOff(PPHOT_SPATIAL_BIN,limiting_mag),
                                               TblDbl,"rms",TblOff(PPHOT_SPATIAL_BIN,rms),
                                               TblInt,"n1",TblOff(PPHOT_SPATIAL_BIN,n1),
                                               TblInt,"n2",TblOff(PPHOT_SPATIAL_BIN,n2),
                                               TblDbl,"max_bright_mag",TblOff(PPHOT_SPATIAL_BIN,max_bright_mag),
                                               TblDbl,"max_bright_iso",TblOff(PPHOT_SPATIAL_BIN,max_bright_iso),
                                               TblDbl,"limiting_iso",TblOff(PPHOT_SPATIAL_BIN,limiting_iso),
                                               TblDbl,"upper_limit",TblOff(PPHOT_SPATIAL_BIN,upper_limit),
                                               TblDbl,"colorterm",TblOff(PPHOT_SPATIAL_BIN,colorterm),
                                               TblDbl,"errorcolor",TblOff(PPHOT_SPATIAL_BIN,errorcolor),
                                               TblInt,"colorflag",TblOff(PPHOT_SPATIAL_BIN,colorflag),
                                               0,"end",0);
          if (tmp_spatial_bin_table == NULL) {
            printf("ERROR: Failed to read table for %s\n",spatial_bins_name);
            errorFlag = 1;

          } else {
            if (pUpdateCommon->verbose) {
              printf("read %zu records for %s\n",spatial_bin_nrecs,spatial_bins_name);
            }
          }
        }

        if (!errorFlag) {
          memset(spatial_bin_table,0,sizeof(spatial_bin_table));
          for (index = 0; index <= MAX_SPATIAL_BINS; index++) {
            spatial_bin_table[index].colorflag = 0;
          }

          /* Now copy the spatial bin entries into their proper slots */
          for (index = 0; index < spatial_bin_nrecs; index++) {
            spatial_bin = tmp_spatial_bin_table[index].spatial_bin;
            pSpatial_bin2 = &tmp_spatial_bin_table[index];
            if (isinf(pSpatial_bin2->errorcolor)) {
              pSpatial_bin2->errorcolor = 99.0;
            }

            if ((spatial_bin > 0) && (spatial_bin <= MAX_SPATIAL_BINS)) {
              pSpatial_bin = &spatial_bin_table[spatial_bin];

              memcpy(pSpatial_bin,&tmp_spatial_bin_table[index],sizeof(PHOT_SPATIAL_BIN));
              if ((pSpatial_bin->colorflag < COLORFLAG_NONE) || (pSpatial_bin->colorflag >= COLORFLAG_MAX)) {
                printf("ERROR: Illegal color flag %d in spatial_bin %d\n",pSpatial_bin->colorflag,spatial_bin);
                exit(1);
              }
            } else {
              printf("ERROR: Illegal spatial bin %d\n",spatial_bin);
              exit(1);
            }
          }
        }

        if (!errorFlag && spatialBinFlag) {
          for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
            pSpatial_bin = &spatial_bin_table[spatial_bin];

            if (pSpatial_bin->spatial_bin == spatial_bin) {
              sprintf(queryString,"REPLACE into spatialbin%s (seriesId,plateNumber,%s,spatial_bin,versionId,rms,n1,n2,limiting_mag,limiting_iso,max_bright_mag,max_bright_iso,upper_limit,colorterm,errorcolor,colorflag) values (%d,%d,%d,%d,%d,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,'%s');\n",
                      pUpdateCommon->catalogString,
                      solutionNumberString,
                      pPhotPlates->seriesId,
                      pPhotPlates->plateNumber,
                      solutionNumber,
                      spatial_bin,
                      pPhotGlobal->currentVersion,
                      pSpatial_bin->rms,
                      pSpatial_bin->n1,
                      pSpatial_bin->n2,
                      pSpatial_bin->limiting_mag,
                      pSpatial_bin->limiting_iso,
                      pSpatial_bin->max_bright_mag,
                      pSpatial_bin->max_bright_iso,
                      pSpatial_bin->upper_limit,
                      pSpatial_bin->colorterm,
                      pSpatial_bin->errorcolor,
                      colortable[pSpatial_bin->colorflag]);

              if (strlen(queryString) > MAX_QUERY_STRING) {
                printf("ERROR: MAX_QUERY_STRING (2) exceeded %zu\n",strlen(queryString));
                return 1;
              }

              res = ExecuteQuery(pUpdateCommon->pPhotConnection,queryString);
              if (res) {
                exit(1);
              }
            }
          }
        }
      }

      // Load local, magnitude-dependent calibration. This can be partial or missing.

      if (magdepFlag && !errorFlag) {
        strcpy(magdep_name,ingestDirectory);
        strcat(magdep_name,"/");
        strcat(magdep_name,pPlateEntry->listName);
        strcat(magdep_name,solutionString);
        strcat(magdep_name,qualifier);
        strcat(magdep_name,"_magdepcalibrate.db");

        if((magdep_load_status = LoadMagdepCorrections(magdep_name,pMagdepLimits,&pMagdepTable,pUpdateCommon->verbose)) == 0) {
          sprintf(globalQuery2String,"REPLACE into magdepcalibrate%s (seriesId,plateNumber,magdep_bin,solutionNumber,versionId,ixb,iyb,imagb,nstar_magdep,magdep_bin_size,xcoord_magdep,ycoord_magdep,magdep_bin_edge,magdep_bin_median,magdep_bin_magcor,magcal_magdep_rms,magdep_bin_quality) VALUES",pUpdateCommon->catalogString);
          curQuery2Count = 0;
          globalQuery2StringSize = strlen(globalQuery2String);
          globalQuery2StringInitial = globalQuery2StringSize;

          for (magdep_index = 0; magdep_index < pMagdepLimits->nrecs; magdep_index++) {
            pMagdep = &pMagdepTable[magdep_index];
            sprintf(queryString," (%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f),",
                    pPhotPlates->seriesId,
                    pPhotPlates->plateNumber,
                    pMagdep->magdep_bin,
                    solutionNumber,
                    pPhotGlobal->currentVersion,
                    pMagdep->ixb,
                    pMagdep->iyb,
                    pMagdep->imagb,
                    pMagdep->nstar_magdep,
                    pMagdep->magdep_bin_size,
                    pMagdep->xcoord_magdep,
                    pMagdep->ycoord_magdep,
                    pMagdep->magdep_bin_edge,
                    pMagdep->magdep_bin_median,
                    pMagdep->magdep_bin_magcor,
                    pMagdep->magcal_magdep_rms,
                    pMagdep->magdep_bin_quality);
            AddToQuery(pUpdateCommon->pPhotConnection,
                       queryString,
                       &globalQuery2String,
                       &globalQuery2StringSize,
                       &globalQuery2StringAlloc,
                       &curQuery2Count,
                       globalQuery2StringInitial,
                       1,
                       MAX_QUERY_COUNT);
          }

          ChargeTime(0,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
          CompleteQuery(pUpdateCommon->pPhotConnection,
                        globalQuery2String,
                        globalQuery2StringSize,
                        &curQuery2Count);
          ChargeTime(6,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
          free(pMagdepTable);
          FreeMagdepSubarrays(pMagdepLimits);
          pMagdepTable = NULL;
        }
      }

      // "Lowess" spatial-binned photometry calibration

      if (lowessGridFlag && !errorFlag) {
        gridWorkDone = 0;
        sprintf(globalQuery2String,"REPLACE into calibration%s (seriesId,plateNumber,solutionNumber,spatial_bin,versionId,maggrid,isogrid,griderr,flaggrid) VALUES",pUpdateCommon->catalogString);
        curQuery2Count = 0;
        globalQuery2StringSize = strlen(globalQuery2String);
        globalQuery2StringInitial = globalQuery2StringSize;

        for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
          strcpy(grid_name,binoutputDirectory);
          strcat(grid_name,"/");
          strcat(grid_name,pPlateEntry->listName);
          strcat(grid_name,solutionString);
          strcat(grid_name,qualifier);
          sprintf(spatial_bin_string,"_a%d.grid",spatial_bin);
          strcat(grid_name,spatial_bin_string);
          grid_file = fopen(grid_name,"rt");

          if (grid_file != NULL) {
            gridCounter = 0;
            gridState = GRID_STATE_INIT;

            while (1) {
              inBuffer = fgets(inLine,MAX_BUFFER,grid_file);
              if (inBuffer == NULL) {
                break;
              }

              nvals = sscanf(inLine,"%f %f %f %d",&maggrid[gridCounter],&isogrid[gridCounter],&griderr[gridCounter],&flaggrid[gridCounter]);

              if (nvals == 4) {
                if (gridState == GRID_STATE_INIT) {
                  if (flaggrid[gridCounter] == 0) {
                    gridState = GRID_STATE_BEGIN_INTERPOLATE;
                    gridCounter++;
                  } else {
                    printf("ERROR: Format 1 error for lowess calibration file %s\n",grid_name);
                    break;
                  }
                } else if (gridState == GRID_STATE_BEGIN_INTERPOLATE) {
                  if (flaggrid[gridCounter] == 0) {
                    lastmaggrid = maggrid[gridCounter];
                    lastisogrid = isogrid[gridCounter];
                    lastgriderr = griderr[gridCounter];
                    lastflaggrid = flaggrid[gridCounter];
                  } else if (flaggrid[gridCounter] == 1) {
                    maggrid[gridCounter+1] = maggrid[gridCounter];
                    isogrid[gridCounter+1] = isogrid[gridCounter];
                    griderr[gridCounter+1] = griderr[gridCounter];
                    flaggrid[gridCounter+1] = flaggrid[gridCounter];
                    maggrid[gridCounter] = lastmaggrid;
                    isogrid[gridCounter] = lastisogrid;
                    griderr[gridCounter] = lastgriderr;
                    flaggrid[gridCounter] = lastflaggrid;
                    gridState = GRID_STATE_NORMAL;
                    gridCounter += 2;
                  } else {
                    printf("ERROR: Format 2 error for lowess calibration file %s\n",grid_name);
                    break;
                  }
                } else if (gridState == GRID_STATE_NORMAL) {
                  if (flaggrid[gridCounter] == 1) {
                    gridCounter++;
                  } else if (flaggrid[gridCounter] == -1) {
                    gridCounter++;
                    gridState = GRID_STATE_END_INTERPOLATE;
                  } else {
                    printf("ERROR: Format 3 error for lowess calibration file %s\n",grid_name);
                    break;
                  }
                } else if (gridState == GRID_STATE_END_INTERPOLATE) {
                  if (flaggrid[gridCounter] != -1) {
                    printf("ERROR: Format 4 error for lowess calibration file %s\n",grid_name);
                    break;
                  }
                }
              }

              if (gridCounter+5 >= MAX_GRID_ENTRY) {
                printf("ERROR: lowess calibration file %s has count %d greater than max count %d\n",grid_name,gridCounter,MAX_GRID_ENTRY);
                break;
              }
            }

            gridCounter++;

            if (gridState == GRID_STATE_END_INTERPOLATE) {
              gridWorkDone = 1;

              for (index = 0; index < gridCounter; index++) {
                sprintf(queryString," (%d,%d,%d,%d,%d,%f,%f,%f,%d),",
                        pPhotPlates->seriesId,
                        pPhotPlates->plateNumber,
                        solutionNumber,
                        spatial_bin,
                        pPhotGlobal->currentVersion,
                        maggrid[index],
                        isogrid[index],
                        griderr[index],
                        flaggrid[index]);
                AddToQuery(pUpdateCommon->pPhotConnection,
                           queryString,
                           &globalQuery2String,
                           &globalQuery2StringSize,
                           &globalQuery2StringAlloc,
                           &curQuery2Count,
                           globalQuery2StringInitial,
                           1,
                           MAX_QUERY_COUNT);
              }
            }

            fclose(grid_file);
          }
        }

        if (gridWorkDone == 0) {
          errorFlag = 1;

          if (solutionNumber == 0) {
            printf("ERROR: Failed to find the lowess calibration file %s for solution %d\n",pPlateEntry->listName,solutionNumber);
          }
        } else {
          ChargeTime(0,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
          CompleteQuery(pUpdateCommon->pPhotConnection,
                        globalQuery2String,
                        globalQuery2StringSize,
                        &curQuery2Count);
          ChargeTime(5,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
        }
      }

      // Extinction and "local bin" data

      if (!errorFlag) {
        strcpy(extinction_name,catalogallDirectory);
        strcat(extinction_name,"/");
        strcat(extinction_name,pPlateEntry->listName);
        strcat(extinction_name,solutionString);
        strcat(extinction_name,qualifier);
        strcat(extinction_name,"_extinction.db");

        extinction_handle = Open(extinction_name,"r");
        if (extinction_handle == NULL) {
          printf("ERROR: Failed to find the extinction file %s\n",extinction_name);
          errorFlag = 1;
        } else {
          if (pUpdateCommon->verbose > 0) {
            printf("Found extinction file %s\n",extinction_name);
          }
        }
      }

      if (!errorFlag) {
        /* Now read in the extinction file */
        extinction_header = table_header(extinction_handle,TABLE_PARSE);
        if (extinction_header == NULL) {
          printf("ERROR: Failed to read header for %s\n",extinction_name);
          errorFlag = 1;
        }
      }

      if (!errorFlag) {
        extinction_table[solutionNumber] = table_loadva(extinction_handle,
                                                        &extinction_header,
                                                        NULL, /* hbase */
                                                        NULL, /* rows */
                                                        NULL,
                                                        sizeof(EXTINCTION),
                                                        &extinction_nrecs,
                                                        TblInt,"enx",TblOff(PEXTINCTION,enx),
                                                        TblInt,"eny",TblOff(PEXTINCTION,eny),
                                                        TblInt,"eix",TblOff(PEXTINCTION,eix),
                                                        TblInt,"eiy",TblOff(PEXTINCTION,eiy),
                                                        TblDbl,"altitude",TblOff(PEXTINCTION,altitude),
                                                        TblDbl,"extinction",TblOff(PEXTINCTION,extinction),
                                                        0,"end",0);

        if (extinction_table[solutionNumber] == NULL) {
          printf("ERROR: Failed to read table for %s\n",extinction_name);
          errorFlag = 1;
        }

        if (pUpdateCommon->verbose > 0) {
          printf("read %zu records for %s\n",extinction_nrecs,extinction_name);
        }
      }

      /* Find out which local smoothing bins are valid for this file */
      if (!errorFlag) {
        strcpy(dmagcor_name,ingestDirectory);
        strcat(dmagcor_name,"/");
        strcat(dmagcor_name,pPlateEntry->listName);
        strcat(dmagcor_name,solutionString);
        strcat(dmagcor_name,qualifier);
        strcat(dmagcor_name,"_dmagcor.db");
        dmagcor_handle = Open(dmagcor_name,"r");

        if (dmagcor_handle == NULL) {
          printf("ERROR: Failed to find the local calibration file %s\n",dmagcor_name);
          errorFlag = 1;
        } else {
          if (pUpdateCommon->verbose) {
            printf("Found local calibration file %s\n",dmagcor_name);
          }
        }
      }

      /* Now read in the local calibration file */

      if (!errorFlag) {
        dmagcor_header = table_header(dmagcor_handle,TABLE_PARSE);
        if (dmagcor_header == NULL) {
          printf("ERROR: Failed to read header for %s\n",dmagcor_name);
          errorFlag = 1;
        }
      }

      if (!errorFlag) {
        dmagcor_table[solutionNumber] = table_loadva(dmagcor_handle,
                                                      &dmagcor_header,
                                                      NULL, /* hbase */
                                                      NULL, /* rows */
                                                      NULL,
                                                      sizeof(DMAGCOR),
                                                      &dmagcor_nrecs,
                                                      TblInt,"nx",TblOff(PDMAGCOR,nx),
                                                      TblInt,"ny",TblOff(PDMAGCOR,ny),
                                                      TblInt,"ix",TblOff(PDMAGCOR,ix),
                                                      TblInt,"iy",TblOff(PDMAGCOR,iy),
                                                      TblInt,"rejectFlag",TblOff(PDMAGCOR,rejectFlag),
                                                      TblInt,"npout",TblOff(PDMAGCOR,npout),
                                                      TblDbl,"zout",TblOff(PDMAGCOR,zout),
                                                      TblDbl,"errout",TblOff(PDMAGCOR,errout),
                                                      TblDbl,"bright_magcor_local",TblOff(PDMAGCOR,bright_magcor_local),
                                                      TblDbl,"bright_magcor_error",TblOff(PDMAGCOR,bright_magcor_error),
                                                      TblInt,"bright_npoints_local",TblOff(PDMAGCOR,bright_npoints_local),
                                                      TblInt,"drad_bin_count",TblOff(PDMAGCOR,drad_bin_count),
                                                      TblInt,"drad_bin_size",TblOff(PDMAGCOR,drad_bin_size),
                                                      TblInt,"drad_reject_count",TblOff(PDMAGCOR,drad_reject_count),
                                                      TblDbl,"draMedian",TblOff(PDMAGCOR,draMedian),
                                                      TblDbl,"draRMS",TblOff(PDMAGCOR,draRMS),
                                                      TblDbl,"ddecMedian",TblOff(PDMAGCOR,ddecMedian),
                                                      TblDbl,"ddecRMS",TblOff(PDMAGCOR,ddecRMS),
                                                      TblInt,"drad_bin_count2",TblOff(PDMAGCOR,drad_bin_count2),
                                                      TblInt,"drad_bin_size2",TblOff(PDMAGCOR,drad_bin_size2),
                                                      TblInt,"drad_reject_count2",TblOff(PDMAGCOR,drad_reject_count2),
                                                      TblDbl,"draMedian2",TblOff(PDMAGCOR,draMedian2),
                                                      TblDbl,"draRMS2",TblOff(PDMAGCOR,draRMS2),
                                                      TblDbl,"ddecMedian2",TblOff(PDMAGCOR,ddecMedian2),
                                                      TblDbl,"ddecRMS2",TblOff(PDMAGCOR,ddecRMS2),
                                                      TblDbl,"dradRMS2",TblOff(PDMAGCOR,dradRMS2),
                                                      0,"end",0);

        if (dmagcor_table[solutionNumber] == NULL) {
          printf("ERROR: Failed to read table for %s\n",dmagcor_name);
          errorFlag = 1;
        }
      }

      if (!errorFlag && localBinFlag) {
        sprintf(globalQuery1String,"REPLACE into localbin%s (seriesId,plateNumber,%s,local_bin_index,versionId,altitude,extinction,magcor_local,magcal_local_error,npoints_local,rejectFlag,drad_bin_count,drad_bin_size,drad_reject_count,draMedian,draRMS,ddecMedian,ddecRMS,drad_bin_count2,drad_bin_size2,drad_reject_count2,draMedian2,draRMS2,ddecMedian2,ddecRMS2,dradRMS2,bright_magcor_local,bright_magcor_error,bright_npoints_local) VALUES",pUpdateCommon->catalogString,solutionNumberString);
        curQuery1Count = 0;
        globalQuery1StringSize = strlen(globalQuery1String);
        globalQuery1StringInitial = globalQuery1StringSize;

        /* We are now ready to write the local bins table */
        if (!errorFlag) {
          nx = XDmagBins(pPlateEntry->mosaicWidth,pPlateEntry->mosaicHeight);
          ny = YDmagBins(pPlateEntry->mosaicWidth,pPlateEntry->mosaicHeight);
          for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
              irec = iy + (ny*ix);
              local_bin_index = ix + (nx*iy);
              pExtinction = &extinction_table[solutionNumber][irec];

              if ((pExtinction->enx != nx) ||
                  (pExtinction->eny != ny) ||
                  (pExtinction->eix != ix) ||
                  (pExtinction->eiy != iy)) {
                printf("ERROR Mismatch of nx %d %d, ny %d %d, ix %d %d, or iy %d %d for %s\n",
                        pExtinction->enx,nx,
                        pExtinction->eny,ny,
                        pExtinction->eix,ix,
                        pExtinction->eiy,iy,
                        extinction_name);
                return 1;
              }

              pDmagcor = &dmagcor_table[solutionNumber][irec];

              if ((pDmagcor->nx != nx) ||
                  (pDmagcor->ny != ny) ||
                  (pDmagcor->ix != ix) ||
                  (pDmagcor->iy != iy)) {
                printf("ERROR Mismatch of nx %d %d, ny %d %d, ix %d %d, or iy %d %d for %s\n",
                        pDmagcor->nx,nx,
                        pDmagcor->ny,ny,
                        pDmagcor->ix,ix,
                        pDmagcor->iy,iy,
                        dmagcor_name);
                return 1;
              }

              sprintf(queryString," (%d,%d,%d,%d,%d,%f,%f,%f,%f,%d,%d,%d,%d,%d,%f,%f,%f,%f,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%d),",
                      pPhotPlates->seriesId,
                      pPhotPlates->plateNumber,
                      solutionNumber,
                      local_bin_index,
                      pPhotGlobal->currentVersion,
                      pExtinction->altitude,
                      pExtinction->extinction,
                      pDmagcor->zout,
                      pDmagcor->errout,
                      pDmagcor->npout,
                      pDmagcor->rejectFlag,
                      pDmagcor->drad_bin_count,
                      pDmagcor->drad_bin_size,
                      pDmagcor->drad_reject_count,
                      pDmagcor->draMedian,
                      pDmagcor->draRMS,
                      pDmagcor->ddecMedian,
                      pDmagcor->ddecRMS,
                      pDmagcor->drad_bin_count2,
                      pDmagcor->drad_bin_size2,
                      pDmagcor->drad_reject_count2,
                      pDmagcor->draMedian2,
                      pDmagcor->draRMS2,
                      pDmagcor->ddecMedian2,
                      pDmagcor->ddecRMS2,
                      pDmagcor->dradRMS2,
                      pDmagcor->bright_magcor_local,
                      pDmagcor->bright_magcor_error,
                      pDmagcor->bright_npoints_local);
              AddToQuery(pUpdateCommon->pPhotConnection,
                          queryString,
                          &globalQuery1String,
                          &globalQuery1StringSize,
                          &globalQuery1StringAlloc,
                          &curQuery1Count,
                          globalQuery1StringInitial,
                          1,
                          MAX_QUERY_COUNT);
            }
          }

          ChargeTime(0,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
          CompleteQuery(pUpdateCommon->pPhotConnection,
                        globalQuery1String,
                        globalQuery1StringSize,
                        &curQuery1Count);
          ChargeTime(1,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
        }
      }

      if (errorFlag) {
        if (solutionNumber) {
          errorFlag = 0;
          maxSolutionNumber = solutionNumber - 1;
        }

        break;
      }
    }

    // We've loaded all of the calibration data for every solution. Now actually
    // load up the `allobjects.db` source table.

    if (magnitudeFlag || platesFlag) {
      if (!errorFlag) {
        if (allobjects_header != NULL) {
          table_hdrfree(allobjects_header);
          allobjects_header = NULL;
        }

        if (allobjects_handle != NULL) {
          Close(allobjects_handle);
          allobjects_handle = NULL;
        }

        if (allobjects_table != NULL) {
          Free(allobjects_table);
        }

        /* Now build and open the allobjects file */
        strcpy(allobjects_name,ingestDirectory);
        strcat(allobjects_name,"/");
        strcat(allobjects_name,pPlateEntry->listName);
        strcat(allobjects_name,qualifier);
        strcat(allobjects_name,"_allobjects.db");

        allobjects_handle = Open(allobjects_name,"r");
        if (allobjects_handle == NULL) {
          printf("ERROR: Could not open file %s\n",allobjects_name);
          errorFlag = 1;
        }
      }

      if (!errorFlag) {
        /* Now read in the allobjects file */
        allobjects_header = table_header(allobjects_handle,TABLE_PARSE);
        if (allobjects_header == NULL) {
          printf("ERROR: Failed to read header for %s\n",allobjects_name);
          errorFlag = 1;
        }
      }

      if (!errorFlag) {
        /* At this point, check the file to see if the exposureNumber and solutionNumber fields are present */
        gotExposureNumber = 0;
        gotSolutionNumber = 0;

        for (columnIndex = 1; columnIndex <= allobjects_header->header->ncol; columnIndex++) {
          if (strcmp(allobjects_header->header->column[columnIndex],"exposureNumber") == 0) {
            gotExposureNumber = 1;
          }
          if (strcmp(allobjects_header->header->column[columnIndex],"solutionNumber") == 0) {
            gotSolutionNumber = 1;
          }
        }

        if (gotExposureNumber && gotSolutionNumber) {
          allobjects_table = table_loadva(allobjects_handle,
                                          &allobjects_header,
                                          NULL, /* hbase */
                                          NULL, /* rows */
                                          NULL,
                                          sizeof(ALLOBJECTSRECORD),
                                          &allobjects_nrecs,
                                          TblBuf,"REF"    ,TblOff(PALLOBJECTSRECORD,REF),MAX_REF,
                                          TblInt,"NUMBER",TblOff(PALLOBJECTSRECORD,NUMBER),
                                          TblDbl,"X_IMAGE"     ,TblOff(PALLOBJECTSRECORD,X_IMAGE),
                                          TblDbl,"Y_IMAGE"     ,TblOff(PALLOBJECTSRECORD,Y_IMAGE),
                                          TblInt,"AFLAGS",TblOff(PALLOBJECTSRECORD,AFLAGS),
                                          TblInt,"BFLAGS",TblOff(PALLOBJECTSRECORD,BFLAGS),
                                          TblInt,"npoints_local",TblOff(PALLOBJECTSRECORD,npoints_local),
                                          TblDbl,"MAG_ISO"     ,TblOff(PALLOBJECTSRECORD,MAG_ISO),
                                          TblDbl,"ra"     ,TblOff(PALLOBJECTSRECORD,ra),
                                          TblDbl,"dec"         ,TblOff(PALLOBJECTSRECORD,dec),
                                          TblDbl,"magcal_iso"    ,TblOff(PALLOBJECTSRECORD,magcal_iso),
                                          TblDbl,"magcal_iso_rms"    ,TblOff(PALLOBJECTSRECORD,magcal_iso_rms),
                                          TblDbl,"magcal_local"  ,TblOff(PALLOBJECTSRECORD,magcal_local),
                                          TblDbl,"magcal_local_error"    ,TblOff(PALLOBJECTSRECORD,magcal_local_error),
                                          TblDbl,"magcal_local_rms"    ,TblOff(PALLOBJECTSRECORD,magcal_local_rms),
                                          TblInt,"spatial_bin",TblOff(PALLOBJECTSRECORD,spatial_bin),
                                          TblInt,"local_bin",TblOff(PALLOBJECTSRECORD,local_bin_index),
                                          TblDbl,"Date"    ,TblOff(PALLOBJECTSRECORD,Date),
                                          TblDbl,"limiting_mag_local"    ,TblOff(PALLOBJECTSRECORD,limiting_mag_local),
                                          TblDbl,"extinction"    ,TblOff(PALLOBJECTSRECORD,extinction),
                                          TblDbl,"Stdmag"    ,TblOff(PALLOBJECTSRECORD,Stdmag),
                                          TblDbl,"color"    ,TblOff(PALLOBJECTSRECORD,color),
                                          TblDbl,"dra"     ,TblOff(PALLOBJECTSRECORD,dra),
                                          TblDbl,"ddec"    ,TblOff(PALLOBJECTSRECORD,ddec),
                                          TblDbl,"FLUX_ISO"     ,TblOff(PALLOBJECTSRECORD,FLUX_ISO),
                                          TblDbl,"MAG_APER"     ,TblOff(PALLOBJECTSRECORD,MAG_APER),
                                          TblDbl,"MAG_AUTO"     ,TblOff(PALLOBJECTSRECORD,MAG_AUTO),
                                          TblDbl,"KRON_RADIUS"     ,TblOff(PALLOBJECTSRECORD,KRON_RADIUS),
                                          TblDbl,"BACKGROUND"     ,TblOff(PALLOBJECTSRECORD,BACKGROUND),
                                          TblDbl,"THRESHOLD"     ,TblOff(PALLOBJECTSRECORD,THRESHOLD),
                                          TblDbl,"FLUX_MAX"     ,TblOff(PALLOBJECTSRECORD,FLUX_MAX),
                                          TblDbl,"THETA_J2000"     ,TblOff(PALLOBJECTSRECORD,THETA_J2000),
                                          TblDbl,"ELLIPTICITY"     ,TblOff(PALLOBJECTSRECORD,ELLIPTICITY),
                                          TblDbl,"ISOAREA_WORLD"     ,TblOff(PALLOBJECTSRECORD,ISOAREA_WORLD),
                                          TblDbl,"FWHM_IMAGE"     ,TblOff(PALLOBJECTSRECORD,FWHM_IMAGE),
                                          TblDbl,"FWHM_WORLD"     ,TblOff(PALLOBJECTSRECORD,FWHM_WORLD),
                                          TblInt,"ISO0",TblOff(PALLOBJECTSRECORD,ISO0),
                                          TblInt,"ISO1",TblOff(PALLOBJECTSRECORD,ISO1),
                                          TblInt,"ISO2",TblOff(PALLOBJECTSRECORD,ISO2),
                                          TblInt,"ISO3",TblOff(PALLOBJECTSRECORD,ISO3),
                                          TblInt,"ISO4",TblOff(PALLOBJECTSRECORD,ISO4),
                                          TblInt,"ISO5",TblOff(PALLOBJECTSRECORD,ISO5),
                                          TblInt,"ISO6",TblOff(PALLOBJECTSRECORD,ISO6),
                                          TblInt,"ISO7",TblOff(PALLOBJECTSRECORD,ISO7),
                                          TblDbl,"plate_dist"     ,TblOff(PALLOBJECTSRECORD,plate_dist),
                                          TblDbl,"Blendedmag"     ,TblOff(PALLOBJECTSRECORD,Blendedmag),
                                          TblInt,"gsc_bin_index",TblOff(PALLOBJECTSRECORD,gsc_bin_index),
                                          TblInt,"exposureNumber",TblOff(PALLOBJECTSRECORD,exposureNumber),
                                          TblInt,"solutionNumber",TblOff(PALLOBJECTSRECORD,solutionNumber),
                                          TblInt,"magdep_bin",TblOff(PALLOBJECTSRECORD,magdep_bin),
                                          TblDbl,"magcal_magdep"     ,TblOff(PALLOBJECTSRECORD,magcal_magdep),
                                          TblDbl,"magcal_magdep_rms"     ,TblOff(PALLOBJECTSRECORD,magcal_magdep_rms),
                                          0,"end",0);
        } else {
          allobjects_table = table_loadva(allobjects_handle,
                                          &allobjects_header,
                                          NULL, /* hbase */
                                          NULL, /* rows */
                                          NULL,
                                          sizeof(ALLOBJECTSRECORD),
                                          &allobjects_nrecs,
                                          TblBuf,"REF"    ,TblOff(PALLOBJECTSRECORD,REF),MAX_REF,
                                          TblInt,"NUMBER",TblOff(PALLOBJECTSRECORD,NUMBER),
                                          TblDbl,"X_IMAGE"     ,TblOff(PALLOBJECTSRECORD,X_IMAGE),
                                          TblDbl,"Y_IMAGE"     ,TblOff(PALLOBJECTSRECORD,Y_IMAGE),
                                          TblInt,"AFLAGS",TblOff(PALLOBJECTSRECORD,AFLAGS),
                                          TblInt,"BFLAGS",TblOff(PALLOBJECTSRECORD,BFLAGS),
                                          TblInt,"npoints_local",TblOff(PALLOBJECTSRECORD,npoints_local),
                                          TblDbl,"MAG_ISO"     ,TblOff(PALLOBJECTSRECORD,MAG_ISO),
                                          TblDbl,"ra"     ,TblOff(PALLOBJECTSRECORD,ra),
                                          TblDbl,"dec"         ,TblOff(PALLOBJECTSRECORD,dec),
                                          TblDbl,"magcal_iso"    ,TblOff(PALLOBJECTSRECORD,magcal_iso),
                                          TblDbl,"magcal_iso_rms"    ,TblOff(PALLOBJECTSRECORD,magcal_iso_rms),
                                          TblDbl,"magcal_local"  ,TblOff(PALLOBJECTSRECORD,magcal_local),
                                          TblDbl,"magcal_local_error"    ,TblOff(PALLOBJECTSRECORD,magcal_local_error),
                                          TblDbl,"magcal_local_rms"    ,TblOff(PALLOBJECTSRECORD,magcal_local_rms),
                                          TblInt,"spatial_bin",TblOff(PALLOBJECTSRECORD,spatial_bin),
                                          TblInt,"local_bin",TblOff(PALLOBJECTSRECORD,local_bin_index),
                                          TblDbl,"Date"    ,TblOff(PALLOBJECTSRECORD,Date),
                                          TblDbl,"limiting_mag_local"    ,TblOff(PALLOBJECTSRECORD,limiting_mag_local),
                                          TblDbl,"extinction"    ,TblOff(PALLOBJECTSRECORD,extinction),
                                          TblDbl,"Stdmag"    ,TblOff(PALLOBJECTSRECORD,Stdmag),
                                          TblDbl,"color"    ,TblOff(PALLOBJECTSRECORD,color),
                                          TblDbl,"dra"     ,TblOff(PALLOBJECTSRECORD,dra),
                                          TblDbl,"ddec"    ,TblOff(PALLOBJECTSRECORD,ddec),
                                          TblDbl,"FLUX_ISO"     ,TblOff(PALLOBJECTSRECORD,FLUX_ISO),
                                          TblDbl,"MAG_APER"     ,TblOff(PALLOBJECTSRECORD,MAG_APER),
                                          TblDbl,"MAG_AUTO"     ,TblOff(PALLOBJECTSRECORD,MAG_AUTO),
                                          TblDbl,"KRON_RADIUS"     ,TblOff(PALLOBJECTSRECORD,KRON_RADIUS),
                                          TblDbl,"BACKGROUND"     ,TblOff(PALLOBJECTSRECORD,BACKGROUND),
                                          TblDbl,"THRESHOLD"     ,TblOff(PALLOBJECTSRECORD,THRESHOLD),
                                          TblDbl,"FLUX_MAX"     ,TblOff(PALLOBJECTSRECORD,FLUX_MAX),
                                          TblDbl,"THETA_J2000"     ,TblOff(PALLOBJECTSRECORD,THETA_J2000),
                                          TblDbl,"ELLIPTICITY"     ,TblOff(PALLOBJECTSRECORD,ELLIPTICITY),
                                          TblDbl,"ISOAREA_WORLD"     ,TblOff(PALLOBJECTSRECORD,ISOAREA_WORLD),
                                          TblDbl,"FWHM_IMAGE"     ,TblOff(PALLOBJECTSRECORD,FWHM_IMAGE),
                                          TblDbl,"FWHM_WORLD"     ,TblOff(PALLOBJECTSRECORD,FWHM_WORLD),
                                          TblInt,"ISO0",TblOff(PALLOBJECTSRECORD,ISO0),
                                          TblInt,"ISO1",TblOff(PALLOBJECTSRECORD,ISO1),
                                          TblInt,"ISO2",TblOff(PALLOBJECTSRECORD,ISO2),
                                          TblInt,"ISO3",TblOff(PALLOBJECTSRECORD,ISO3),
                                          TblInt,"ISO4",TblOff(PALLOBJECTSRECORD,ISO4),
                                          TblInt,"ISO5",TblOff(PALLOBJECTSRECORD,ISO5),
                                          TblInt,"ISO6",TblOff(PALLOBJECTSRECORD,ISO6),
                                          TblInt,"ISO7",TblOff(PALLOBJECTSRECORD,ISO7),
                                          TblDbl,"plate_dist"     ,TblOff(PALLOBJECTSRECORD,plate_dist),
                                          TblDbl,"Blendedmag"     ,TblOff(PALLOBJECTSRECORD,Blendedmag),
                                          TblInt,"gsc_bin_index",TblOff(PALLOBJECTSRECORD,gsc_bin_index),
                                          TblInt,"magdep_bin",TblOff(PALLOBJECTSRECORD,magdep_bin),
                                          TblDbl,"magcal_magdep"     ,TblOff(PALLOBJECTSRECORD,magcal_magdep),
                                          TblDbl,"magcal_magdep_rms"     ,TblOff(PALLOBJECTSRECORD,magcal_magdep_rms),
                                          0,"end",0);
        }

        if (allobjects_table == NULL) {
          printf("ERROR: Failed to  readr %s\n",allobjects_name);
          errorFlag = 1;
        }

        if (pUpdateCommon->verbose) {
          printf("read %zu records for %s\n",allobjects_nrecs,allobjects_name);
        }

        if (!gotExposureNumber || !gotSolutionNumber) {
          for (allobjects_index = 0; allobjects_index < allobjects_nrecs; allobjects_index++) {
            pAllobjects = &allobjects_table[allobjects_index];
            pAllobjects->exposureNumber = 0;
            pAllobjects->solutionNumber = 0;
          }
        }
      }

      if (!errorFlag) {
        THRESHOLD = -1.0;

        // Sort in order of increasing gsc_bin_index
        qsort((void*)allobjects_table, allobjects_nrecs, sizeof(ALLOBJECTSRECORD), AllobjectsCompare);

        pPhotPlates->unmatched = 0;
        binPrintFlag = 1;
        sanityErrorCount = 0;

        // For each source ...

        for (allobjects_index = 0; allobjects_index < allobjects_nrecs; allobjects_index++) {
          pAllobjects = &allobjects_table[allobjects_index];

          // Load and homogenize metadata for this measurement

          if (GetREFNumber(pAllobjects->REF,&pAllobjects->REFNumber,&refType,0,1) != 0) {
            printf("ERROR: Bad REF  in catalog %s for NUMBER %d in %s\n",catalogText[pUpdateCommon->pMagFileCache->catalogNumber],pAllobjects->NUMBER,allobjects_name);
            errorFlag = 1;
            break;
          }

          if ((pUpdateCommon->pMagFileCache->catalogNumber == CATALOG_KEPLER) &&
              ((refType ==  REF_TYPE_GSC) || (refType ==  REF_TYPE_APASS))) {
            printf("ERROR: Found refType %d in catalog %s for NUMBER %d in %s\n",refType,catalogText[pUpdateCommon->pMagFileCache->catalogNumber],pAllobjects->NUMBER,allobjects_name);
            exit(1);
          } else if ((pUpdateCommon->pMagFileCache->catalogNumber == CATALOG_APASS) && (refType ==  REF_TYPE_KEPLER)) {
            printf("ERROR: Found refType %d in catalog %s for NUMBER %d in %s\n",refType,catalogText[pUpdateCommon->pMagFileCache->catalogNumber],pAllobjects->NUMBER,allobjects_name);
            exit(1);
          } else if ((pUpdateCommon->pMagFileCache->catalogNumber == CATALOG_EXPERIMENTAL) &&
                     ((refType ==  REF_TYPE_GSC) || (refType ==  REF_TYPE_KEPLER))) {
            printf("ERROR: Found refType %d in catalog %s for NUMBER %d in %s\n",refType,catalogText[pUpdateCommon->pMagFileCache->catalogNumber],pAllobjects->NUMBER,allobjects_name);
            exit(1);
          } else if ((pUpdateCommon->pMagFileCache->catalogNumber == CATALOG_GAIA) && (refType !=  REF_TYPE_GAIA2) && (refType != REF_TYPE_NONE)) {
            printf("ERROR: Found refType %d in catalog %s for NUMBER %d in %s\n",refType,catalogText[pUpdateCommon->pMagFileCache->catalogNumber],pAllobjects->NUMBER,allobjects_name);
            exit(1);
          } else if ((pUpdateCommon->pMagFileCache->catalogNumber == CATALOG_ATLAS) && (refType !=  REF_TYPE_ATLAS2) && (refType != REF_TYPE_NONE)) {
            printf("ERROR: Found refType %d in catalog %s for NUMBER %d in %s\n",refType,catalogText[pUpdateCommon->pMagFileCache->catalogNumber],pAllobjects->NUMBER,allobjects_name);
            exit(1);
          }

          /* Set the QUALITY fit if needed */
          if ((pAllobjects->AFLAGS & (1<<FILTER_AFLAG_QUALITY)) == 0) {
            if (pPlateEntry->quality != 0) {
              pAllobjects->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
            } else {
              switch(pUpdateCommon->pMagFileCache->catalogNumber) {
              case CATALOG_GSC232:
                if ((spatial_bin_table[pAllobjects->spatial_bin].colorflag != COLORFLAG_METROPOLIS) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm > MAX_BLUE_COLORTERM) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm < MIN_BLUE_COLORTERM)) {
                  pAllobjects->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
                }
                break;

              case CATALOG_KEPLER:
                if ((spatial_bin_table[pAllobjects->spatial_bin].colorflag != COLORFLAG_METROPOLIS) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm > MAX_KEPLER_BLUE_COLORTERM) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm < MIN_KEPLER_BLUE_COLORTERM)) {
                  pAllobjects->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);

                }
                break;

              case CATALOG_APASS:
                if ((spatial_bin_table[pAllobjects->spatial_bin].colorflag != COLORFLAG_METROPOLIS) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm > MAX_APASS_BLUE_COLORTERM) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm < MIN_APASS_BLUE_COLORTERM)) {
                  pAllobjects->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
                }
                break;

              case CATALOG_GAIA:
                if ((spatial_bin_table[pAllobjects->spatial_bin].colorflag != COLORFLAG_METROPOLIS) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm > MAX_GAIA_BLUE_COLORTERM) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm < MIN_GAIA_BLUE_COLORTERM)) {
                  pAllobjects->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
                }
                break;

              case CATALOG_ATLAS:
                if ((spatial_bin_table[pAllobjects->spatial_bin].colorflag != COLORFLAG_METROPOLIS) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm > MAX_ATLAS_BLUE_COLORTERM) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm < MIN_ATLAS_BLUE_COLORTERM)) {
                  pAllobjects->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
                }
                break;

              case CATALOG_EXPERIMENTAL:
                if ((spatial_bin_table[pAllobjects->spatial_bin].colorflag != COLORFLAG_METROPOLIS) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm > MAX_EXPERIMENTAL_BLUE_COLORTERM) ||
                    (spatial_bin_table[pAllobjects->spatial_bin].colorterm < MIN_EXPERIMENTAL_BLUE_COLORTERM)) {
                  pAllobjects->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
                }
                break;

              default:
                printf("ERROR: unknown catalog %d\n",pUpdateCommon->pMagFileCache->catalogNumber);
                exit(1);
                break;
              }
            }
          }

          if ((pAllobjects->REFNumber == 0)  && (pAllobjects->AFLAGS == 0)) {
            /* We are looking for high numbers of unmatched objects */
            pPhotPlates->unmatched++;
          }

          if ((pAllobjects->AFLAGS & (1<<FILTER_AFLAG_BLEND_NOMATCH)) != 0) {
            pUpdateCommon->rejectedBlendNomatch++;
            continue;
          }

          if (THRESHOLD < 0) {
            THRESHOLD = pAllobjects->THRESHOLD;
            pPlateEntry->exposures = 1;
            if ((magnitudeFlag == 0) && (platesFlag != 0)) {
              /* The threshold is all we need if we are not writing the magnitudes or overlap table  */
              break;
            }
          } else {
            if (THRESHOLD != pAllobjects->THRESHOLD) {
              printf("ERROR: THRESHOLD %f and %f disagree\n",THRESHOLD,pAllobjects->THRESHOLD);
            }
          }

          nx = XDmagBins(pPlateEntry->mosaicWidth,pPlateEntry->mosaicHeight);
          ny = YDmagBins(pPlateEntry->mosaicWidth,pPlateEntry->mosaicHeight);
          ix = (pAllobjects->X_IMAGE * nx) /(1.0 * pPlateEntry->mosaicWidth);
          iy = (pAllobjects->Y_IMAGE * ny) /(1.0 * pPlateEntry->mosaicHeight);

          if (ix < 0) {
            ix = 0;
          }

          if (ix >= nx) {
            ix = nx-1;
          }

          if (iy < 0) {
            iy = 0;
          }

          if (iy >= ny) {
            iy = ny-1;
          }

          /* These records are reversed from the local bin */
          irec = iy + (ny*ix);
          local_bin_index = ix + (nx*iy);

          if (pAllobjects->local_bin_index != local_bin_index) {
            printf("ERROR: local_bin_index %d %d mismatch for %d in %s\n",pAllobjects->local_bin_index,local_bin_index,pAllobjects->NUMBER,allobjects_name);
            errorFlag = 1;
            break;
          }

          // Construct the magnitude record in our buffer

          if (magnitudeFlag) {
            solutionNumber = pAllobjects->solutionNumber;

            if (solutionNumber > maxSolutionNumber) {
              if (printMaxSolutionError) {
                printMaxSolutionError = 0;
                printf("ERROR: solutionNumber %d exceeds maxSolutionNumber %d for %s\n",solutionNumber,maxSolutionNumber,allobjects_name);
              }
              solutionNumber = 0;
            }

            pDmagcor = &dmagcor_table[solutionNumber][irec];
            if ((pDmagcor->iy + (pDmagcor->ny * pDmagcor->ix)) != irec) {
              printf("ERROR: local bin table out of order for %s\n",dmagcor_name);
              exit(1);
            }

            pExtinction = &extinction_table[solutionNumber][irec];
            if ((pExtinction->eiy + (pExtinction->eny * pExtinction->eix)) != irec) {
              printf("ERROR: extinction table out of order for %s\n",extinction_name);
              exit(1);
            }

            memset(pFileStarImage, 0, sizeof(FILESTARIMAGE));

            pFileStarImage->gsc_bin_index =  pAllobjects->gsc_bin_index;
            pFileStarImage->seriesId =  pPhotPlates->seriesId;
            pFileStarImage->plateNumber =  pPhotPlates->plateNumber;
            pFileStarImage->NUMBER =  pAllobjects->NUMBER;
            pFileStarImage->passBits = PASSBIT_GSC232;
            pFileStarImage->versionId =  pPhotGlobal->currentVersion;
            pFileStarImage->exposureNumber =  pAllobjects->exposureNumber;
            pFileStarImage->solutionNumber =  pAllobjects->solutionNumber;
            pFileStarImage->spatial_bin =  pAllobjects->spatial_bin;
            pFileStarImage->local_bin_index =  local_bin_index;
            pFileStarImage->REFNumber =  pAllobjects->REFNumber;
            pFileStarImage->X_IMAGE =  pAllobjects->X_IMAGE;
            pFileStarImage->Y_IMAGE =  pAllobjects->Y_IMAGE;
            pFileStarImage->AFLAGS =  pAllobjects->AFLAGS;
            pFileStarImage->BFLAGS =  pAllobjects->BFLAGS;
            pFileStarImage->MAG_ISO =  pAllobjects->MAG_ISO;
            pFileStarImage->ra =  pAllobjects->ra;
            pFileStarImage->dec =  pAllobjects->dec;
            pFileStarImage->magcal_iso =  pAllobjects->magcal_iso;
            pFileStarImage->magcal_iso_rms =  pAllobjects->magcal_iso_rms;
            pFileStarImage->magcal_local =  pAllobjects->magcal_local;
            pFileStarImage->magcal_local_rms =  pAllobjects->magcal_local_rms;
            pFileStarImage->Date =  pAllobjects->Date;
            pFileStarImage->FLUX_ISO =  pAllobjects->FLUX_ISO;
            pFileStarImage->MAG_APER =  pAllobjects->MAG_APER;
            pFileStarImage->MAG_AUTO =  pAllobjects->MAG_AUTO;
            pFileStarImage->KRON_RADIUS =  pAllobjects->KRON_RADIUS;
            pFileStarImage->BACKGROUND =  pAllobjects->BACKGROUND;
            pFileStarImage->FLUX_MAX =  pAllobjects->FLUX_MAX;
            pFileStarImage->THETA_J2000 =  pAllobjects->THETA_J2000;
            pFileStarImage->ELLIPTICITY =  pAllobjects->ELLIPTICITY;
            pFileStarImage->ISOAREA_WORLD =  pAllobjects->ISOAREA_WORLD;
            pFileStarImage->FWHM_IMAGE =  pAllobjects->FWHM_IMAGE;
            pFileStarImage->FWHM_WORLD =  pAllobjects->FWHM_WORLD;
            pFileStarImage->ISO0 =  pAllobjects->ISO0;
            pFileStarImage->ISO1 =  pAllobjects->ISO1;
            pFileStarImage->ISO2 =  pAllobjects->ISO2;
            pFileStarImage->ISO3 =  pAllobjects->ISO3;
            pFileStarImage->ISO4 =  pAllobjects->ISO4;
            pFileStarImage->ISO5 =  pAllobjects->ISO5;
            pFileStarImage->ISO6 =  pAllobjects->ISO6;
            pFileStarImage->ISO7 =  pAllobjects->ISO7;
            pFileStarImage->magdep_bin = pAllobjects->magdep_bin;
            pFileStarImage->magcal_magdep = pAllobjects->magcal_magdep;
            pFileStarImage->magcal_magdep_rms = pAllobjects->magcal_magdep_rms;
            pFileStarImage->plate_dist =  pAllobjects->plate_dist;
            pFileStarImage->Blendedmag =  pAllobjects->Blendedmag;
            pFileStarImage->limiting_mag_local =  pAllobjects->limiting_mag_local;
            /* New V6 fields */
            pFileStarImage->A2FLAGS = 0;
            pFileStarImage->B2FLAGS = 0;
            pFileStarImage->maskIndex = 0;

            if ((pAllobjects->solutionNumber >= 0)  && (pAllobjects->solutionNumber < MAX_CATALOG_EXPOSURES)) {
              pFileStarImage->timeAccuracy = timeAccuracy[pAllobjects->solutionNumber];
            } else {
              pFileStarImage->timeAccuracy = -1;
            }

            if (pFileStarImage->timeAccuracy == 0) {
              if (timeAccuracyError == 0) {
                printf("ERROR: zero timeAccuracy for %s\n",pPlateEntry->listName);
                timeAccuracyError = 1;
              }
              pFileStarImage->timeAccuracy = 1;
            }

            /* Catalog position */
            pFileStarImage->RaPM = NULL_PROPER_MOTION;
            pFileStarImage->DecPM = NULL_PROPER_MOTION;
            pFileStarImage->ra_2 = 999.0;
            pFileStarImage->dec_2 = 99.0;

            if (pFileStarImage->AFLAGS == 0 && (pAllobjects->magcal_local_error > 90.0 || pDmagcor->errout > 90.0)) {
              sanityErrorCount++;

              if (sanityErrorCount == 1) {
                printf("ERROR: AFLAGS %d and magcal_local_error %f %f ix %d iy %d local_bin %d are inconsistent for %d in %s\n",
                        pFileStarImage->AFLAGS,
                        pDmagcor->errout,
                        pAllobjects->magcal_local_error,
                        pDmagcor->ix,
                        pDmagcor->iy,
                        pFileStarImage->local_bin_index,
                        pFileStarImage->NUMBER,
                        pPlateEntry->listName);
              }
              continue;
            }

            if (fabs(pAllobjects->magcal_local_error - pDmagcor->errout) > 0.01) {
              sanityErrorCount++;

              if (sanityErrorCount == 1) {
                printf("ERROR: AFLAGS %d and magcal_local_err or %f %f ix %d iy %d local_bin %d are inconsistent for %d in %s\n",
                        pFileStarImage->AFLAGS,
                        pDmagcor->errout,
                        pAllobjects->magcal_local_error,
                        pDmagcor->ix,
                        pDmagcor->iy,
                        pFileStarImage->local_bin_index,
                        pFileStarImage->NUMBER,
                        pPlateEntry->listName);
              }
              continue;
            }

            pFileStarImage->magcal_local_error = pDmagcor->errout;
            pFileStarImage->dradRMS2 = pDmagcor->dradRMS2;
            pFileStarImage->magcor_local = pDmagcor->zout;
            pFileStarImage->extinction = pExtinction->extinction;
            pFileStarImage->npoints_local = pDmagcor->npout;
            pFileStarImage->rejectFlag = pDmagcor->rejectFlag;

            if (pFileStarImage->gsc_bin_index == 0) {
              tmp_gsc_bin_index = GetGSCBin(pGscBin,pFileStarImage->ra,pFileStarImage->dec,&decBin,&raBin,"binCheck");

              if (tmp_gsc_bin_index != pFileStarImage->gsc_bin_index) {
                if (binPrintFlag != 0) {
                  printf("ERROR: gsc_bin_index %d is not %d for NUMBER %d in  %s\n",tmp_gsc_bin_index,pFileStarImage->gsc_bin_index,pAllobjects->NUMBER,allobjects_name);
                  binPrintFlag = 0;
                }
                continue;
              }
            }

            ChargeTime(0,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
            StoreMagnitudeFileEntry(pUpdateCommon->pMagFileCache,pFileStarImage,pUpdateCommon->photfilebase); // result ignored
            ChargeTime(2,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
            magnitudeLimitCount++;
            pUpdateCommon->magnitudeCount++;
            plateMagnitudeCount++;
          }
        }

        if (sanityErrorCount > 0) {
          printf("ERROR: sanity check inconsistent: %d in %s\n",sanityErrorCount,pPlateEntry->listName);
          errorFlag = 1;
        }

        // Emit any remaining buffered magnitudes to disk (StoreMagnitudeFileEntry() might also flush)
        FlushCache(pUpdateCommon->pMagFileCache, pUpdateCommon->photfilebase);
      }
    }

    // Emit information about this plate, if requested

    if (platesFlag) {
      if (!errorFlag) {
        if (GetMargins(pUpdateCommon->pConnection,
                       pPhotPlates->series,
                       pPhotPlates->plateNumber,
                       pPhotPlates->mosaicNumber,
                       &pPhotPlates->leftMargin,
                       &pPhotPlates->rightMargin,
                       &pPhotPlates->bottomMargin,
                       &pPhotPlates->topMargin) != 0) {
          pPhotPlates->leftMargin = -1;
          pPhotPlates->rightMargin = -1;
          pPhotPlates->bottomMargin = -1;
          pPhotPlates->topMargin = 1;
        }

        sprintf(queryString,"INSERT IGNORE into photplates (seriesId,plateNumber) values (%d,%d);\n",
                pPhotPlates->seriesId,pPhotPlates->plateNumber);
        res = ExecuteQuery(pUpdateCommon->pPhotConnection,queryString);
        if (res) {
          return 1;
        }

        sprintf(queryString,"UPDATE photplates SET mosaicNumber = %d,versionId%s = %d,THRESHOLD = %f,scale = %f ,mosaicWidth = %d,mosaicHeight = %d,exposures = %d,leftMargin = %d,rightMargin = %d,topMargin = %d,bottomMargin = %d,quality = %d,unmatched = %d WHERE seriesid = %d and plateNumber = %d;\n",
                pPhotPlates->mosaicNumber,
                pUpdateCommon->catalogString,
                pPhotGlobal->currentVersion,
                THRESHOLD,
                pPlateEntry->scale,
                pPlateEntry->mosaicWidth,
                pPlateEntry->mosaicHeight,
                pPlateEntry->exposures,
                pPhotPlates->leftMargin,
                pPhotPlates->rightMargin,
                pPhotPlates->topMargin,
                pPhotPlates->bottomMargin,
                pPhotPlates->quality,
                pPhotPlates->unmatched,
                pPhotPlates->seriesId,
                pPhotPlates->plateNumber);
        ChargeTime(0,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);
        res = ExecuteQuery(pUpdateCommon->pPhotConnection,queryString);
        ChargeTime(4,&pUpdateCommon->curUsec,pUpdateCommon->deltaTable,pUpdateCommon->totalTable);

        if (res == 0) {
          seriesId = GetSeriesId(pPlateEntry->series,1);
          sprintf(queryString,"UPDATE photinsert set mosaicNumber = %d,versionId = %d,photFailFlag = %d WHERE seriesId = %d and plateNumber = %d and catalogNumber = %d;\n",
                  pPlateEntry->mosaicNumber,
                  pPhotGlobal->currentVersion,
                  0,
                  seriesId,
                  pPlateEntry->plateNumber,
                  pUpdateCommon->pMagFileCache->catalogNumber);
          res = ExecuteQuery(pUpdateCommon->pPhotConnection,queryString);
          if (res != 0) {
            exit(1);
          }
        }
      }
    }

    // We're done processing this plate; emit a status update

    time(&curPlateTime);
    pUpdateCommon->curTime = curPlateTime - pUpdateCommon->startTime;
    curPlateTime -= pUpdateCommon->plateStartTime;
    plateQueryCount = GetQueryCount() - plateQueryCount;

    if (!errorFlag) {
      double magRate = 0.0;
      double queryRate = 0.0;
      double plateRate = 0.0;

      pUpdateCommon->curPlateCount++;

      if (pUpdateCommon->curTime > 0) {
        magRate = 1.0 * pUpdateCommon->magnitudeCount / pUpdateCommon->curTime;
        queryRate = 1.0 * GetQueryCount() / pUpdateCommon->curTime;
        plateRate = 1.0 * pUpdateCommon->curPlateCount / pUpdateCommon->curTime;
      }

      printf(
        "plate %d %19s mag=%d query=%d sec=%ld | "
        "tot: mag=%d query=%d sec=%ld | "
        "rate: mag=%.0f query=%.1f plate=%.3f\n",
        pUpdateCommon->plateIndex + 1,
        pPlateEntry->listName,
        plateMagnitudeCount,
        plateQueryCount,
        curPlateTime,
        pUpdateCommon->magnitudeCount,
        GetQueryCount(),
        pUpdateCommon->curTime,
        magRate,
        queryRate,
        plateRate
      );
      ChargeTime(0, &pUpdateCommon->curUsec, pUpdateCommon->deltaTable, pUpdateCommon->totalTable);
      ChargeTime(RESET_DELTA_ENTRY, &pUpdateCommon->curUsec, pUpdateCommon->deltaTable, pUpdateCommon->totalTable);
      fflush(pUpdateCommon->outHandle);
    }
  }

  // Done with all plates. Cleanup.

  WriteFinalStatistics(pUpdateCommon);

  mysql_close(pUpdateCommon->pPhotConnection);
  mysql_close(pUpdateCommon->pConnection);

  if (listStrings != NULL) {
    free(listStrings);
  }

  if (pPlateTable != NULL) {
    free(pPlateTable);
  }

  if (pUpdateCommon->outHandle != NULL) {
    fclose(pUpdateCommon->outHandle);
  }

  if (tmp_spatial_bin_table != NULL) {
    Free(tmp_spatial_bin_table);
    tmp_spatial_bin_table = NULL;
  }

  if (spatial_bins_header != NULL) {
    table_hdrfree(spatial_bins_header);
    spatial_bins_header = NULL;
  }

  if (spatial_bins_handle != NULL) {
    Close(spatial_bins_handle);
    spatial_bins_handle = NULL;
  }

  for (solutionNumber = 0; solutionNumber < MAX_CATALOG_EXPOSURES; solutionNumber++) {
    if (extinction_table[solutionNumber] != NULL) {
      Free(extinction_table[solutionNumber]);
      extinction_table[solutionNumber] = NULL;
    }

    if (dmagcor_table[solutionNumber] != NULL) {
      Free(dmagcor_table[solutionNumber]);
      dmagcor_table[solutionNumber] = NULL;
    }
  }

  if (extinction_header != NULL) {
    table_hdrfree(extinction_header);
    extinction_header = NULL;
  }

  if (extinction_handle != NULL) {
    Close(extinction_handle);
    extinction_handle = NULL;
  }

  if (dmagcor_header != NULL) {
    table_hdrfree(dmagcor_header);
    dmagcor_header = NULL;
  }

  if (dmagcor_handle != NULL) {
    Close(dmagcor_handle);
    dmagcor_handle = NULL;
  }

  if (allobjects_header != NULL) {
    table_hdrfree(allobjects_header);
    allobjects_header = NULL;
  }

  if (allobjects_handle != NULL) {
    Close(allobjects_handle);
    allobjects_handle = NULL;
  }

  if (allobjects_table != NULL) {
    Free(allobjects_table);
  }

  if (globalQuery1String != NULL) {
    free(globalQuery1String);
  }

  if (globalQuery2String != NULL) {
    free(globalQuery2String);
  }

  return EXIT_SUCCESS;
}


/*
 * Dec 18, 2008 Edward J. Los - Original Version
 * Jan 16, 2008 Edward J. Los - Flush output file to monitor job progress
 *                              Implement the keepalive mechanism to kill the job cleanly
 * Jan 21, 2008 Edward J. Los - Add additional statistics, remove scale calculations.
 * Jan 23, 2008 Edward J. Los - Add skip flag for stars, magnitudes, and overlaps.
 * Jan 28, 2009 Edward J. Los - Split FLAGS into AFLAGS and BFLAGS, remove CAL_FLAG
 * Feb  6, 2009 Edward J. Los - Revise flag handling
 * Feb 11, 2009 Edward J. Los - Limit the drad error to the series matching limit
 * Feb 14, 2009 Edward J. Los - Presort the overlap table
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 * Mar  2, 2009 Edward J. Los - Remove overlap table support
 *                              Remove drad table and dradAverage since all values are now in dmagcor.
 *                              Read gsc_bin_index prior to precession,
 * May  4, 2009 Edward J. Los - Implement automatic version update
 * May  5, 2009 Edward J. Los - Add mosaicNubmer to photplates
 * May 12, 2009 Edward J. Los - Add a check for bad gsc_bin_index
 * Jun  3, 2009 Edward J. Los - Replace series with seriesId
 * Jun 23, 2009 Edward J. Los - Partition the magnitude table.
 *                              Add plate margins and plate quality.
 * Jun 30, 2009 Edward J. Los - Check for an illegal GSC bin index
 * Jul 20, 2009 Edward J. Los - Add detailed statistics.
 *                              Do not store FILTER_AFLAG_BLEND_NOMATCH entries in the database
 * Jul 27, 2009 Edward J. Los - Add a keepalive database call prior to long MySQL transactions.
 *                            - Do not bump the version if only the platesFlag is set
 * Aug 12, 2009 Edward J. Los - Correct an "Inf" in the spatial bin error color by setting it to +99
 * Aug 20, 2009 Edward J. Los - Add Kepler Input Catalog support
 * Sep 19, 2009 Edward J. Los - Use files instead of MySQL magnitude tables
 * Sep 30, 2009 Edward J. Los - Read the entire allobjects table at once for performance.
 * Oct 15, 2009 Edward J. Los - Add DASCH_PHOT_MAGNITUDES for the photometry files.
 * Nov 16, 2009 Edward J. Los - add noNewVersion feature
 * Nov 20, 2009 Edward J. Los - Add multiple exposure support
 * Dec  2, 2009 Edward J. Los - Store the magnitudeFileFlag in the photglobal table
 * Dec 29, 2009 Edward J. Los - Sync the quality column with the same column in the plates table of the scanner database
 * Jan 27, 2010 Edward J. Los - Add rejectFlag to the file-based database
 * Feb 10, 2010 Edward J. Los - Fix solutionNumber column.
 *                              Count the number of good unmatched stars
 * Feb 12, 2010 Edward J. Los - Allow update_photometry to run from the scanweb account
 *                              Add temporary quality support
 * Feb 23, 2010 Edward J. Los - Improve error reporting
 * Apr  2, 2010 Edward J. Los - Support kepler magnitude files
 * May  5, 2010 Edward J. Los - Add error recovery from a lost server
 * Nov  5, 2010 Edward J. Los - Change group ownership of files to the scanner group
 * Dec 18, 2010 Edward J. Los - Invoke ChangeGroup to change the group name
 * Feb  2, 2011 Edward J. Los - Write the lowess calibration tables
 * Feb 28, 2011 Edward J. Los - add bright local calibration parameters
 * Mar 14, 2011 Edward J. Los - Correct local bin indices.
 * Mar 23, 2011 Edward J. Los - Add apass catalog support
 * Apr  8, 2011 Edward J. Los - Add a sanity check to avoid mixing catalogs.
 * Nov  8, 2011 Edward J. Los - Add magnitude-dependent photometry correction
 * Jul 23, 2012 Edward J. Los - Allow the pipeline files to be on the machine with the DASCH_SCRIPTS directory.
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Aug 16, 2012 Edward J. Los - Add sanity checks for possible corruption
 * Sep 10, 2012 Edward J. Los - Reorganize code to optimize magnitude file support.
 * Nov  9, 2012 Edward J. Los - Correct colorterm limit code.
 * Nov 23, 2012 Edward J. Los - Implement the star hash table for performance optimisation
 * Mar  4, 2013 Edward J. Los - Add dra and ddec to the magnitude file (Version 5);
 * Mar  5, 2013 Edward J. Los - Check the mosaics table for stale solutions
 * Mar 11, 2013 Edward J. Los - Add RaPM, DecPM, ra_2 and dec_2 to the magnitude file (Version 5);
 * Mar 25, 2013 Edward J. Los - Use a different semaphores for each catalog database.
 * May 31, 2013 Edward J. Los - Correct a bug in handling multiple exposures
 * Jun 29, 2013 Edward J. Los - Correct initialization error for the fullSucceededFlag
 * Jan 20, 2015 Edward J. Los - V6 data format: add A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex
 * Mar  3, 2015 Edward J. Los - Add a trap for zero timeAccuracy
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * Mar 21, 2015 Edward J. Los - Deprecate the stars table.  Replacing with starcatalog;
 * May 29, 2018 Edward J. Los - Support gaia
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 * Jan 29, 2019 Edward J. Los - Use inputs from a different disk for Odyssey queue processing
 * Mar  4, 2019 Edward J. Los - Introduce '-D' to specify the raid array for the input starbase files
 * May 29, 2019 Edward J. Los - Handle GetREFNumber error with better descriptive material
 * Feb 18, 2020 Edward J. Los - Do not insert rejected mosaics with markings
 * Jun 13, 2020 Edward J. Los - Always store the new version in the logs when '-n' is specified
 *                              Add photinsert table: set photFailFlag when before processing a plate and clear photFailFlag if processing succeeds
 * Nov  9, 2020 Edward J. Los - Exclude 'na' series plates from insertion
 */
