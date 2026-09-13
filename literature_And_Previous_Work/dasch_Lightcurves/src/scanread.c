// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* Adapted from libwcs/scanread.c */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <fcntl.h>

#include "mysql.h"

#include "libwcs/wcs.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcscat.h"

#include "scandb.h"
#include "pipelineutils.h"
#include "scanread.h"

extern int getdateform(void);
int CheckAuthorization(char *regionFlag,double ra, double dec,int *pReleaseField);
int dcmp(const void *a,const void * b);

#ifndef O_BINARY
#define O_BINARY 0
#endif

#define STATUS_STRING_SIZE 200
#define MAX_MOSAICID_STRING 50
#define WEB_SIZE_LIMIT 80000 /* High enough so we never see a problem (3000 plates) */

/* Guarantee that the object is on the plate by using the shorter length.
   Take 90% of this length to ensure that the object is not in bin 9 */
#define PESSIMISTIC_HALFWIDTH10X10  96.52 /* Assume a  7.6" x  7.6"  plate */
#define PESSIMISTIC_HALFWIDTH17X17 168.90 /* Assume a 13.3" x 13.3"  plate */
/* Normal search - large area.  These sizes must match ProcLogbook.h */
#define OPTIMISTIC_HALFWIDTH10X10 127.0 /* Assume a 10"x10" plate */
#define OPTIMISTIC_HALFWIDTH17X17 215.9 /* Assume a 17"x17" plate */

static char *searchPlateSeries = NULL;
static char *searchPlateNumber = NULL;
static char *searchPlateClass = NULL;
static char searchType[10] = {"wcsfit"};
static char *webtable = NULL;
static char *webextract = NULL;
static char *tarballextract = NULL;
static char *webfiletype = NULL;
static char *webarchivetype = NULL;
static char *regionflag = NULL;
static FILE *outfile = NULL;
static double epoch1 = 0.0;
static double epoch2 = 0.0;
static int mosaicFlag = 0;
static int pessimisticFlag = 0;
static int singleColumnFlag = 0;
static int degout = 0;

static int query_tool_mode = 0;

void scanclose();

extern struct WorldCoor *GetFITSWCS();
struct StarCat *scanopen (char *bincat);

static   MYSQL *sqlConnection = NULL;  /* MySQL database connection */
static   PSERIESTABLE sqlSeriesTable = NULL;   /* MySQL database series table */

/* from catutil.c */

void
DS_setlimdeg (int degoutx)
{
  degout = degoutx;
}


void
printTableHeader(
  int showDistance,
  int showHeliocentric,
  int publicRelease
) {
  if (!publicRelease) {
    printf(
      "<p><i>Sorry, cutout images are not available for this source due to "
      "<a target=\"_blank\" href=\"https://dasch.cfa.harvard.edu/data-access/#restrictions\">temporary "
      "data access restrictions</a>.</i></p>\n"
    );
  }

  printf("<table class=\"bordered\" summary=\"A table of search results\">\n");
  printf("<caption>Search Results</caption>\n");
  printf(
    "<tr>"
    "<th>Series</th>"
    "<th>Plate Number (Select to view metadata)</th>"
    "<th>Class</th>"
    "<th>Exposure Number</th>"
    "<th>Right Ascension</th>"
    "<th>Declination</th>"
    "<th>Exposure (Minutes)</th>"
  );

  if (getdateform() == EP_JD) {
    if (showHeliocentric) {
      printf("<th>Heliocentric Julian Date</th>");
    } else {
      printf("<th>Geocentric Julian Date</th>");
    }
  } else {
    printf("<th>Geocentric Date</th>");
  }

  if (showDistance) {
    printf("<th>cm from center</th><th>cm from edge</th>");
  }

  printf("<th>Status");
  if (publicRelease) {
    printf(" (click for cutout)");
  }
  printf("</th>");
  printf("</tr>\n");
}


int
StoreResult(
  PMOSAICLIST pResult,
  PSERIESTABLE seriesTablePtr,
  char *series,
  double	cra,		/* Search center J2000 right ascension in degrees */
  double	cdec,		/* Search center J2000 declination in degrees */
  double  dra,        /* Right ascension half width in degrees */
  double  ddec,       /* Declination half width in degrees */
  double  drad,           /* Extraction radius in degrees */
  double	*tnum,		/* Array of catalog numbers (returned) */
  double	*tra,		/* Array of right ascensions (returned) */
  double	*tdec,		/* Array of declinations (returned) */
  double	**tmag,		/* 2-D Array of magnitudes (returned) */
  char	**tclass,	 /* Array of spectral classes (returned) */
  PSORTTABLE pSortTable, /* Table for sorting by Julian Day */
  int *pSortTableSize,      /* size of the sort table */
  int nstar,               /* Index for the above arrays */
  int nstarmax,              /* Allocated size of the arrays */
  int publicRelease,       /* True if data can be released */
  int *selectWCSFit,       /* Number of images enabled for extraction */
  double *selectWCSSize,   /* Size of extracted regions in bytes */
  int *selectWCSStrlen   /* Size of the plate identification string */
) {
  char statusString[STATUS_STRING_SIZE];
  char distanceString[64];
  int naxis1;
  int naxis2;
  char ctype1[20];
  char ctype2[20];
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double cd1_1;
  double cd1_2;
  double cd2_1;
  double cd2_2;
  double cd[4];
  struct WorldCoor *wcs;
  double xpix;
  double ypix;
  int offscl;
  double halfwidth;
  int showDistance = 0;
  double centerDistance = 100.0;
  double edgeDistance = 100.0;
  int showHeliocentric = 0;
  PSORTTABLE pSortEntry;
  int sortTableSize = *pSortTableSize;
  int have_pos_filter = (cra > -95.0) && (cdec > -95.0);

  if (pResult->plateNumber == 0) {
    return(0);
  }

  if (nstar > nstarmax) {  /* Bugfix of Nov 20, 2019 */
    return(0);
  }

  if (have_pos_filter && getdateform() == EP_JD) {
    showHeliocentric = 1;
  }

  distanceString[0] = 0;
  statusString[0] = 0;

  /* Verify that this object is really on the plate */
  if (have_pos_filter) {
    if (pessimisticFlag == 1) {
      halfwidth = PESSIMISTIC_HALFWIDTH10X10; /* halfwidth mm for a 10" x 10" plate */
      if (strcmp(seriesTablePtr->series,"a") == 0) {
        halfwidth = PESSIMISTIC_HALFWIDTH17X17; /* halfwidth mm for a 17" x 17" plate */
      }
    } else {
      halfwidth = OPTIMISTIC_HALFWIDTH10X10; /* halfwidth mm for a 10" x 10" plate */
      if (strcmp(seriesTablePtr->series,"a") == 0) {
        halfwidth = OPTIMISTIC_HALFWIDTH17X17; /* halfwidth mm for a 17" x 17" plate */
      }
    }

    if (pResult->naxis1 == 0 || pResult->WCSSource != WCSSOURCE_IMWCS) {
      naxis1 = 2.0*halfwidth*PIXELS_PER_MM;
      naxis2 = 2.0*halfwidth*PIXELS_PER_MM;
      strcpy(ctype1,"RA---TAN");
      strcpy(ctype2,"DEC--TAN");
      crval1 = pResult->dRightAscension;
      crval2 = pResult->dDeclination;
      crpix1 = ((1.0 * naxis1)/2.0) + 0.5;
      crpix2 = ((1.0 * naxis2)/2.0) + 0.5;
      cd2_2 = (2.0 * dra)/(1.0 * naxis1);
      cd1_2 = 0.0;
      cd2_1 = 0.0;
      cd1_1 = -cd2_2;
    } else {
      naxis1 = pResult->naxis1;
      naxis2 = pResult->naxis2;
      strcpy(ctype1,pResult->ctype1);
      strcpy(ctype2,pResult->ctype2);
      crval1 = pResult->crval1;
      crval2 = pResult->crval2;
      crpix1 = pResult->crpix1;
      crpix2 = pResult->crpix2;
      cd1_1 = pResult->cd1_1;
      cd1_2 = pResult->cd1_2;
      cd2_1 = pResult->cd2_1;
      cd2_2 = pResult->cd2_2;
    }

    /* If the first axis is DEC, swap things around to ensure conformance
       with the standard */
    if (strstr(ctype1,"DEC")) {
      char tmpctype[20];
      double dtmp;
      strcpy(tmpctype,ctype2);
      strcpy(ctype2,ctype1);
      strcpy(ctype1,tmpctype);
      dtmp = crval1;
      crval1 = crval2;
      crval2 = dtmp;
      dtmp = cd2_1;
      cd2_1 = cd1_1;
      cd1_1 = dtmp;
      dtmp = cd2_2;
      cd2_2 = cd1_2;
      cd1_2 = dtmp;
    }

    cd[0] = cd1_1;
    cd[1] = cd1_2;
    cd[2] = cd2_1;
    cd[3] = cd2_2;

    wcs = wcskinit(
      naxis1,
      naxis2,
      ctype1,
      ctype2,
      crpix1,
      crpix2,
      crval1,
      crval2,
      cd,
      0.0,  /* cdelt1 */
      0.0,  /* cdelt2 */
      0.0,  /* crota */
      2000, /* equinox */
      0.0   /* epoch */
    );

    wcs2pix(wcs, cra, cdec, &xpix, &ypix, &offscl);

    if (webtable != NULL || query_tool_mode || pessimisticFlag == 2) {
      double centerx = xpix - (((1.0 * naxis1)/2.0) + 0.5);
      double centery = ypix - (((1.0 * naxis2)/2.0) + 0.5);

      showDistance = 1;
      centerDistance = sqrt((centerx * centerx) + (centery * centery));
      centerDistance = centerDistance /(10.0 * PIXELS_PER_MM);
      edgeDistance = xpix;

      if (ypix < edgeDistance) {
        edgeDistance = ypix;
      }

      if (((1.0 * naxis1) - xpix) < edgeDistance) {
        edgeDistance = (1.0 * naxis1) - xpix;
      }

      if (((1.0 * naxis2) - ypix) < edgeDistance) {
        edgeDistance = (1.0 * naxis2) - ypix;
      }

      edgeDistance = edgeDistance / (10.0 * PIXELS_PER_MM);
      sprintf(distanceString, "<td>%4.1f</td><td>%4.1f</td>", centerDistance, edgeDistance);
    }

    wcsfree(wcs);

    if (pessimisticFlag == 2) {
      if ((centerDistance/(centerDistance+fabs(edgeDistance))) > 0.3) {
        return(0);
      }
    } else {
      if (offscl != 0) {
        return(0);
      }
    }
  }

  pResult->julianDate = fd2jd(pResult->date);
  if (showHeliocentric) {
    pResult->julianDate = jd2hjd(pResult->julianDate,cra,cdec,WCS_J2000);
  }

  pResult->epoch = jd2ep(pResult->julianDate);

  if (query_tool_mode) {
    char *wcssource;
    float exptime;

    if (nstar == 0) {
      // Print the header
      printf(
        "series\tplatenum\tscannum\tmosnum\texpnum\tsolnum\tclass\t"
        "ra\tdec\texptime\tjd\tepoch\twcssource\t"
        "scale\tscandate\tmosdate\trotation\tbinflags"
      );

      if (have_pos_filter) {
        printf("\tcenterdist\tedgedist");
      }

      printf("\n");
    }

    if (sscanf(pResult->exposure, "%f", &exptime) != 1) {
      exptime = -1.0;
    }

    switch (pResult->WCSSource) {
    case WCSSOURCE_NONE:
    default:
      wcssource = "none";
      break;

    case WCSSOURCE_LOGBOOK:
      wcssource = "logbook";
      break;

    case WCSSOURCE_IMWCS:
      wcssource = "imwcs";
      break;
    }

    printf(
      "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%.8f\t%.8f\t%.1f\t%12.5f\t%.5f\t%s\t%.3f\t%s\t%s\t%d\t%d",
      series,
      pResult->plateNumber,
      pResult->scanNumber,
      pResult->mosaicNumber,
      pResult->exposureNumber,
      pResult->solutionNumber,
      pResult->plateClass,
      pResult->dRightAscension,
      pResult->dDeclination,
      exptime,
      pResult->julianDate,
      pResult->epoch,
      wcssource,
      pResult->plateScale,
      pResult->scanDate,
      pResult->mosaicDate,
      pResult->rotation,
      pResult->binning
    );

    if (have_pos_filter) {
      printf(
        "\t%.1f\t%.1f",
        centerDistance,
        edgeDistance
      );
    }

    printf("\n");
  } else if (webtable == NULL) {
    if (mosaicFlag) {
      int rotationCode = pResult->rotation/90;
      if (pResult->mosaicNumber < 0) {
        /* No mosaic exists */
        rotationCode = 5; /* Flag exposure */
        tnum[nstar] = (1.0 * pResult->plateNumber) + ((1.0 * MAX_PLATE_NUMBER) * seriesTablePtr->seriesId) + ((1.0 * pResult->exposureNumber)/MAX_EXPOSURE_NUMBER)  + ((1.0 * rotationCode)/(MAX_EXPOSURE_NUMBER*MAX_ROTATION));
      } else if (pResult->WCSSource != WCSSOURCE_IMWCS) {
        rotationCode = 4; /* Flag for no WCS fit */
        tnum[nstar] = (1.0 * pResult->plateNumber) + ((1.0 * MAX_PLATE_NUMBER) * seriesTablePtr->seriesId) + ((1.0 * pResult->mosaicNumber)/MAX_EXPOSURE_NUMBER)  + ((1.0 * rotationCode)/(MAX_EXPOSURE_NUMBER*MAX_ROTATION));
      } else {
        /* Store the mosaic number and the rotation */
        tnum[nstar] = (1.0 * pResult->plateNumber) + ((1.0 * MAX_PLATE_NUMBER) * seriesTablePtr->seriesId) + ((1.0 * pResult->mosaicNumber)/MAX_EXPOSURE_NUMBER) + ((1.0 * rotationCode)/(MAX_EXPOSURE_NUMBER*MAX_ROTATION));
      }
    } else {
      /* Store the exposure number */
      tnum[nstar] = (1.0 * pResult->plateNumber) + ((1.0 * MAX_PLATE_NUMBER) * seriesTablePtr->seriesId) + ((1.0 * pResult->exposureNumber)/MAX_EXPOSURE_NUMBER);
    }

    if ((pResult->dRightAscension < 999.0) && (pResult->dDeclination < 99.0)) {
      tra[nstar] = str2ra(pResult->rightAscension);
      tdec[nstar] = str2dec(pResult->declination);
    } else {
      tra[nstar] = pResult->dRightAscension;
      tdec[nstar] = pResult->dDeclination;
    }

    tmag[0][nstar] = fd2ep(pResult->date);
    if (strlen(pResult->plateClass) > 0) {
      tclass[nstar] = (char *)calloc(strlen(pResult->plateClass)+1,1);
      if (tclass[nstar]) {
        strcpy(tclass[nstar],pResult->plateClass);
      }
    } else {
      tclass[nstar] = NULL;
    }
  } else {
    /* Display results as we find them so the web client does not have to wait. */

    /* As best I can tell, WCSSource is more about scan status than WCS status, as follows: */
    switch (pResult->WCSSource) {
    default:
    case WCSSOURCE_NONE:
      strcpy(statusString, "not scanned");
      break;

    case WCSSOURCE_LOGBOOK:
      strcpy(statusString, "scanned, WCS from logbook (Astrometry.Net failed)");
      break;

    case WCSSOURCE_IMWCS:
      if ((webextract == NULL) ||
          (webfiletype == NULL) ||
          (webarchivetype == NULL) ||
          (tarballextract == NULL) ||
          (publicRelease == 0) ||
          (drad == 0) ||
          !have_pos_filter) {
        strcpy(statusString, "scanned, WCS from Astrometry.Net");
      } else {
        int result;
        float plateScale;
        double totalPixels;
        char tempBuff[MAX_MOSAICID_STRING];
        int tempLength;

        result = snprintf(
          statusString,
          sizeof(statusString),
          "<a href=\"%s?filetype=%s&series=%s&plateNumber=%d&mosaicNumber=%d&ra=%f&dec=%f&system=J2000&radius=%f\" "
          "target=\"search\">scanned, WCS from Astrometry.Net</a>",
          webextract,
          webfiletype,
          seriesTablePtr->series,
          pResult->plateNumber,
          pResult->mosaicNumber,
          cra,
          cdec,
          drad
        );

        if (result < 0) {
          printf("buffer overflow in statusString\n");
          exit(1);
        }

        (*selectWCSFit)++;

        plateScale = sqrt((pResult->cd2_1 * pResult->cd2_1) + (pResult->cd2_2 * pResult->cd2_2));
        if (plateScale > 0.0) {
          totalPixels = (2.0*drad/plateScale)*(2.0*drad/plateScale);
          (*selectWCSSize) += totalPixels;
        }

        /* Prepare the tarball command list, which will be colon separated for each plate and
           space separated for the series, plate number, and mosaic number */
        sprintf(
          tempBuff,
          "%s %d %d %9.4f;",
          seriesTablePtr->series,
          pResult->plateNumber,
          pResult->mosaicNumber,
          pResult->epoch
        );

        tempLength = strlen(tempBuff);

        (*selectWCSStrlen) += tempLength;

        tclass[nstar] = (char *) calloc(tempLength + 1, 1);
        if (tclass[nstar]) {
          strcpy(tclass[nstar],tempBuff);
          pSortEntry = &pSortTable[sortTableSize];
          pSortEntry->julianDate = pResult->julianDate;
          pSortEntry->listString = tclass[nstar];
          sortTableSize++;
          *pSortTableSize = sortTableSize;
        }
      }

      break;
    }

    if (nstar == 0) {
      printTableHeader(showDistance, showHeliocentric, publicRelease);
    }

    printf(
      "<tr><td>%s</td><td><a href=\"%s?series=%s&plateNumber=%d",
      series,
      webtable,
      series,
      pResult->plateNumber
    );

    /* Add extraction information if requested */
    if (
      (webextract != NULL) &&
      (webfiletype != NULL) &&
      (webarchivetype != NULL) &&
      (tarballextract != NULL) &&
      (drad > 0.0) &&
      have_pos_filter
    ) {
      printf(
        "&filetype=%s&ra=%f&dec=%f&system=J2000&radius=%f",
        webfiletype,
        cra,
        cdec,
        drad
      );
    }

    if (getdateform() == EP_JD) {
      sprintf(pResult->date, "%12.5f", pResult->julianDate);
    }

    printf(
      "\" target=\"search\">%d</a></td>"
      "<td>%s</td>"
      "<td>%d</td>"
      "<td>%s</td>"
      "<td>%s</td>"
      "<td>%s</td>"
      "<td>%s</td>"
      "%s"
      "<td>%s</td>"
      "</tr>\n",
      pResult->plateNumber,
      pResult->plateClass,
      pResult->exposureNumber,
      pResult->rightAscension,
      pResult->declination,
      pResult->exposure,
      pResult->date,
      distanceString,
      statusString
    );
  }

  /* Mark the table as already printed */
  memset(pResult, 0, sizeof(MOSAICLIST));
  return 1;
}


/* SCANREAD -- Read binary catalog sources + names in specified region */

int
scanread(
  char *bincat, /* Name of reference star catalog file */
  double cra,  /* Search center J2000 right ascension in degrees */
  double cdec,  /* Search center J2000 declination in degrees */
  double drad,           /* Search half-width and half-height in degrees */
  int sysout,  /* Search coordinate system */
  double eqout,  /* Search coordinate equinox */
  double epout,  /* Proper motion epoch (0.0 for no proper motion) */
  int nstarmax, /* Maximum number of sources to be returned */
  struct StarCat **starcat, /* Star catalog data structure */
  double *tnum,  /* Array of catalog numbers (returned) */
  double *tra,  /* Array of right ascensions (returned) */
  double *tdec,  /* Array of declinations (returned) */
  double **tmag,  /* 2-D Array of magnitudes (returned) */
  char **tclass,  /* Array of spectral classes (returned) */
  int  nlog
) {
  double ra1 = 0, ra2 = 0;  /* Limiting output right ascensions of region */
  double dec1 = 0, dec2 = 0;  /* Limiting output declinations of region */
  struct StarCat *sc;  /* Star catalog data structure */
  int nstar;
  int nplates;
  int istar = 0;
  int verbose;
  MYSQL *pConnection;
  PSERIESTABLE seriesTablePtr;
  char series[MAX_SERIES_STRING];
  int nvals;
  int oldPlateNumber = -1;
  int plateNumberCopy2 = -1;
  int oldPlateNumber2 = -1; /* For the plate count */
  int poleflag = 0;
  int tempResult;
  double halfwidth;  /* Half width in mm */
  double dra = 0;        /* Right ascension half width in degrees */
  double ddec = 0;       /* Declination half width in degrees */
  int naxis1;
  int naxis2;
  char ctype1[20];
  char ctype2[20];
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double cd1_1;
  double cd1_2;
  double cd2_1;
  double cd2_2;
  double xmin;
  double ymin;
  double xmax;
  double ymax;
  double sra;
  double sdec;
  int wrap1;
  int wrap2;
  double cd[4];
  MOSAICLIST bestResultTable;
  PMOSAICLIST pBestResult = &bestResultTable;
  int nResult;
  int selectWCSFit = 0;
  double selectWCSSize = 0.0;
  int selectWCSStrlen = 0;
  struct WorldCoor *wcs;
  char archiveString[40];
  int pendingSearch = 0;
  int allSearch = 0;
  int totalSearch = 0; /* Search for wcsfit with the highest good mosaic number */
  int wcsfitSearch = 0; /* Search for all wcsfit solutions */
  int scannedSearch = 0; /* Search for all scanned plates with at least one matched exposureNumber */
  PMOSAICLIST pStoredResultTable;
  PSORTTABLE pSortTable = NULL;
  int sortTableSize = 0;
  int publicRelease = 0;
  int releaseField;
  int pendingFlag;
  int numMosaicRecords = 0;
  int mosaicRecordIndex = 0;
  int mosaicRecordIndex2 = 0;
  PMOSAICLIST pMosaicTable = NULL;
  PMOSAICLIST pMosaicTableEntry;
  PMOSAICLIST pMosaicTableEntry2;
  int pendingPlateNumber = -1;
  int pendingScannedFlag = 0;
  int plateNumber = -1;
  int mosaicNumber = -1;
  int solutionNumber = -1;
  int mosaicInfoFlag = 0;
  int selectAllFlag = 0;
  char *dateStr1 = NULL;
  char *dateStr2 = NULL;
  int totalSolution0Count = 0;
  int totalSolutionsCount = 0;

  memset(pBestResult, 0, sizeof(MOSAICLIST));

  if (webtable) {
    outfile = stdout;
  } else {
    outfile = stderr;
  }

  sc = *starcat;

  if (nlog > 0)
    verbose = 1;
  else
    verbose = 0;

  /* Open catalog */
  if (sc == NULL)
    sc = scanopen (bincat);

  *starcat = sc;
  if (sc == NULL)
    return 0;

  nstar = 0;
  nplates = 0;

  /* If our center is B1950, convert it to J2000 to match the database */
  wcscon(sysout, WCS_J2000, eqout, 2000.0, &cra, &cdec, 2000.0);

  publicRelease = CheckAuthorization(regionflag, cra, cdec, &releaseField);

  pStoredResultTable = (PMOSAICLIST) calloc(MAX_EXPOSURE_NUMBER, sizeof(MOSAICLIST));
  if (pStoredResultTable == NULL) {
    printf("ERROR: Failed to allocate pStoredResultTable\n");
    exit(1);
  }

  pSortTable = (PSORTTABLE) calloc(nstarmax, sizeof(SORTTABLE));
  if (pSortTable == NULL) {
    printf("ERROR: Failed to allocate pSortTable\n");
    exit(1);
  }

  /* The outermost loop goes through the series because each series has a different plate scale */

  pConnection = sqlConnection;
  if (sqlConnection == NULL) {
    printf("ERROR: null sqlConnection in scanread\n");
    exit(1);
  }

  seriesTablePtr = sqlSeriesTable;
  if (seriesTablePtr == NULL) {
    printf("ERROR: null sqlConnection in scanread\n");
    exit(1);
  }

  while (seriesTablePtr->seriesId > 0) {
    strcpy(series, seriesTablePtr->series);

    /* Skip plate scale records. If we are searching on a given series, skip all other records */
    if (
      (seriesTablePtr->plateScale > 0) &&
      (seriesTablePtr->sequestered == SEQUESTERED_NO) &&
      ((searchPlateSeries == NULL) || (strcmp(searchPlateSeries, series) == 0))
    ) {
      pendingFlag = 0;
      pendingSearch = 0;
      allSearch = 0;

      if (strstr(searchType, "pending")) {
        pendingSearch = 1;
        pendingFlag = 1;
      } else if (strstr(searchType, "all")) {
        allSearch = 1;
        pendingFlag = 1;
      } else if (strstr(searchType, "total")) {
        totalSearch = 1;
      } else if (strstr(searchType, "wcsfit")) {
        wcsfitSearch = 1;
      }

      if (strstr(searchType,"scanned")) {
        scannedSearch = 1;
      }

      if (searchPlateNumber != NULL) {
        nvals = sscanf(searchPlateNumber, "%d", &plateNumber);
        if (nvals != 1) {
          printf("ERROR: bad plate number %s in scanread\n", searchPlateNumber);
          exit(1);
        }
      }

      BuildMosaicList(
        pConnection,
        series,     /* if NULL or zero length, consider all series */
        plateNumber,  /* if -1, consider all plates */
        mosaicNumber,  /* if -1, consider all mosaics */
        solutionNumber, /* if -1, consider all solutions */
        pendingFlag, /* if 1, then return unscanned plates */
        mosaicInfoFlag, /* if 1, then being called from GetMosaicInfo: no exposure or plates table info */
        selectAllFlag,   /* if 1, then return all mosaic records: good, deleted, and stale */
        verbose,         /* if 1, display warnings */
        &totalSolutionsCount,
        &totalSolution0Count,
        &numMosaicRecords,
        &pMosaicTable
      );

      /* Now go through the table and exclude anything not requested by the remaining search parameters */

      if (epoch1 > 0.0 && epoch2 > 0.0) {
        dateStr1 = ep2fd(epoch1);
        dateStr2 = ep2fd(epoch2);
      }

      /* See if we are searching on a plate center */

      if (cra > -95.0 && cdec > -95.0) {
        /* Now that we know the scale of the plate, we can compute the search
         * limits. Assume a square plate with North at 0, 90, 180, or 270
         * degrees (not true for at least one J series plate), and then refine
         * the search for each candidate plate once we have the actual WCS. */

        if (pessimisticFlag == 1) {
          halfwidth = PESSIMISTIC_HALFWIDTH10X10; /* halfwidth mm for a 10" x 10" plate */
          if (strcmp(seriesTablePtr->series, "a") == 0) {
            halfwidth = PESSIMISTIC_HALFWIDTH17X17; /* halfwidth mm for a 17" x 17" plate */
          }
        } else {
          halfwidth = OPTIMISTIC_HALFWIDTH10X10; /* halfwidth mm for a 10" x 10" plate */
          if (strcmp(seriesTablePtr->series, "a") == 0) {
            halfwidth = OPTIMISTIC_HALFWIDTH17X17; /* halfwidth mm for a 17" x 17" plate */
          }
        }

        dra = halfwidth * seriesTablePtr->plateScale / 3600.0; /* halfwidth in degrees for RA */

        ddec = dra; /* Assume square plate */

        dec1 = cdec - ddec;
        dec2 = cdec + ddec;

        /* Search zones which include the poles cover 360 degrees in RA */
        poleflag = 0;

        if (cdec - ddec < -90.0) {
          dec1 = -90.0;
          poleflag = 1;
        }

        if (cdec + ddec > 90.0) {
          dec2 = +90.0;
          poleflag = 1;
        }

        /* Now refine the search by checking the corners of a square plate */
        naxis1 = 2.0 * halfwidth * PIXELS_PER_MM;
        naxis2 = 2.0 * halfwidth * PIXELS_PER_MM;
        strcpy(ctype1, "RA---TAN");
        strcpy(ctype2, "DEC--TAN");
        crval1 = cra;
        crval2 = cdec;
        crpix1 = ((1.0 * naxis1)/2.0) + 0.5;
        crpix2 = ((1.0 * naxis2)/2.0) + 0.5;
        cd2_2 = (2.0 * dra)/(1.0 * naxis1);
        cd1_2 = 0.0;
        cd2_1 = 0.0;
        cd1_1 = -cd2_2;

        cd[0] = cd1_1;
        cd[1] = cd1_2;
        cd[2] = cd2_1;
        cd[3] = cd2_2;

        wcs = wcskinit(
          naxis1,
          naxis2,
          ctype1,
          ctype2,
          crpix1,
          crpix2,
          crval1,
          crval2,
          cd,
          0.0,  /* cdelt1 */
          0.0,  /* cdelt2 */
          0.0,  /* crota */
          2000, /* equinox */
          0.0   /* epoch */
        );

        xmin = 0.5;
        ymin = 0.5;
        xmax = 0.5 + (1.0* naxis1);
        ymax = 0.5 + (1.0* naxis2);

        /* Lower left corner */
        pix2wcs(wcs,xmin,ymin,&sra,&sdec);
        ra1 = sra;
        if (sra < cra) {
          /* Flag if we crossed the equinox */
          wrap1 = 1;
        }
        if (sdec < dec1) {
          dec1 = sdec;
        }

        /* Lower right corner */
        pix2wcs(wcs,xmax,ymin,&sra,&sdec);
        ra2 = sra;
        if (sra > cra) {
          /* Flag if we crossed the equinox */
          wrap2 = 1;
        }
        if (sdec < dec1) {
          dec1 = sdec;
        }

        /* Upper left corner */
        pix2wcs(wcs,xmin,ymax,&sra,&sdec);
        if (sra < cra) {
          /* Here we wrapped */
          if (wrap1) {
            /* Both wrapped. Keep greater */
            if (sra > ra1) {
              ra1 = sra;
            }
          } else {
            /* Upper wrapped but lower did not.  Use upper */
            ra1 = sra;
            wrap1 = 1;
          }
        } else {
          /* Here there is no wrap */
          if (wrap1) {
            /* Lower wrapped but upper did not.  Keep lower */
          } else {
            /* Both didn't wrap. Keep greater */
            if (sra > ra1) {
              ra1 = sra;
            }
          }
        }

        if (sdec > dec2) {
          dec2 = sdec;
        }

        /* Upper right corner */
        pix2wcs(wcs,xmax,ymax,&sra,&sdec);
        if (sra > cra) {
          /* Here we wrapped */
          if (wrap2) {
            /* Both wrapped. Keep lesser */
            if (sra < ra2) {
              ra2 = sra;
            }
          } else {
            /* Upper wrapped but lower did not.  Use upper */
            ra2 = sra;
            wrap2 = 1;
          }
        } else {
          /* Here there is no wrap */
          if (wrap2) {
            /* Lower wrapped but upper did not.  Keep lower */
          } else {
            /* Both didn't wrap. Keep lesser */
            if (sra < ra2) {
              ra2 = sra;
            }
          }
        }
        if (sdec > dec2) {
          dec2 = sdec;
        }

        wcsfree(wcs);
      }

      for (mosaicRecordIndex = 0; mosaicRecordIndex < numMosaicRecords; mosaicRecordIndex++) {
        pMosaicTableEntry = &pMosaicTable[mosaicRecordIndex];

        /* A pending search must look forward to see if any entry for this plate has been scanned */
        if (pendingSearch == 1) {
          if (pMosaicTableEntry->plateNumber == pendingPlateNumber && pendingScannedFlag != 0) {
            continue;
          }

          if (pMosaicTableEntry->mosaicNumber >= 0) {
            pendingPlateNumber = pMosaicTableEntry->plateNumber;
            pendingScannedFlag = 1;
            continue;
          }

          pendingPlateNumber = -1;
          pendingScannedFlag = 0;

          for (mosaicRecordIndex2 = mosaicRecordIndex + 1; mosaicRecordIndex2 < numMosaicRecords; mosaicRecordIndex2++) {
            pMosaicTableEntry2 = &pMosaicTable[mosaicRecordIndex2];

            if (pMosaicTableEntry2->plateNumber != pMosaicTableEntry->plateNumber) {
              break;
            }

            if (pMosaicTableEntry2->mosaicNumber >= 0) {
              pendingPlateNumber = pMosaicTableEntry->plateNumber;
              pendingScannedFlag = 1;
              break;
            }
          }

          if (pMosaicTableEntry->plateNumber == pendingPlateNumber && pendingScannedFlag != 0) {
            continue;
          }

          if (pMosaicTableEntry->centerSource == CENTERSOURCE_IMWCS) {
            printf("line %4d WARNING: stale WCSFIT location for  %5s%04d_%02d s%02d exposure %02d WCSSource %d centerSource %d\n",
                    __LINE__,
                    pMosaicTableEntry->series,
                    pMosaicTableEntry->plateNumber,
                    pMosaicTableEntry->mosaicNumber,
                    pMosaicTableEntry->solutionNumber,
                    pMosaicTableEntry->exposureNumber,
                    pMosaicTableEntry->WCSSource,
                    pMosaicTableEntry->centerSource);
            printf("SELECT series,plateNumber,exposureNumber,centerSource,version FROM exposures INNER JOIN logbook using (series,plateNumber,exposureNumber) where series = '%s' and plateNumber = %d and exposureNumber = %d;\n",
                      pMosaicTableEntry->series,
                      pMosaicTableEntry->plateNumber,
                      pMosaicTableEntry->exposureNumber);
            printf("UPDATE logbook set version = '%s' where series = '%s' and plateNumber = %d and exposureNumber = %d;\n",
                      LOGBOOK_VERSION,
                      pMosaicTableEntry->series,
                      pMosaicTableEntry->plateNumber,
                      pMosaicTableEntry->exposureNumber);
            printf("UPDATE exposures set centerSource = 'Logbook' where series = '%s' and plateNumber = %d and exposureNumber = %d;\n",
                      pMosaicTableEntry->series,
                      pMosaicTableEntry->plateNumber,
                      pMosaicTableEntry->exposureNumber);
          }
        }

        if (strstr(pMosaicTableEntry->mosaicComment,"Deleted") != NULL) {
          continue;
        }

        /* Remove rejected dates */
        if (dateStr1 != NULL && dateStr1 != NULL) {
          if (strcmp(pMosaicTableEntry->date, dateStr1) < 0) {
            continue;
          }

          if (strcmp(pMosaicTableEntry->date, dateStr2) > 0) {
            continue;
          }
        }

        /* Remove unwanted classes */
        if (searchPlateClass != NULL) {
          if (strcmp(pMosaicTableEntry->plateClass, searchPlateClass) != 0) {
            continue;
          }
        }

        /* Include only fitted mosaics */
        if (wcsfitSearch == 1 || totalSearch == 1) {
          if (pMosaicTableEntry->WCSSource != WCSSOURCE_IMWCS) {
            continue;
          }
        } else if (scannedSearch == 1) {
          if (pMosaicTableEntry->mosaicNumber < 0) {
            continue;
          }
        }

        /* Include only unscanned plates? */
        if (pendingSearch == 1) {
          if (pMosaicTableEntry->mosaicNumber >= 0) {
            continue;
          }
        }

        /* For wcsfit and scanned, include only solutions matched to a logbook exposure */
        if (
          (wcsfitSearch == 1 || scannedSearch == 1 || allSearch == 1) &&
          pMosaicTableEntry->exposureNumber < 0
        ) {
          continue;
        }

        if (cra > -95.0 && cdec > -95.0) {
          if (pendingFlag == 0) {
            if ((pMosaicTableEntry->dDeclination < dec1) ||
                (pMosaicTableEntry->dDeclination > dec2)) {
              continue;
            }
          }

          if (pendingFlag == 0 && poleflag == 0) {
            if (wrap1 != 0 || wrap2 != 0) {
              if ((pMosaicTableEntry->dRightAscension < ra2) && (pMosaicTableEntry->dRightAscension > ra1)) {
                continue;
              }
            } else {
              if ((pMosaicTableEntry->dRightAscension < ra2) || (pMosaicTableEntry->dRightAscension > ra1)) {
                continue;
              }
            }
          }
        }

        /* If the output is in another system, convert the string value to the
         * other system for display and output, but leave the decimal value in
         * J2000 */
        {
          double dra2 = pMosaicTableEntry->dRightAscension;
          double ddec2 = pMosaicTableEntry->dDeclination;

          if ((dra2 < 999.0) && (ddec2 < 99.0)) {
            if ((sysout != WCS_J2000) || (eqout != 2000.0)) {
              wcscon(WCS_J2000,sysout,2000.0,eqout,&dra2,&ddec2,epout);
            }

            if (degout) {
              sprintf(pMosaicTableEntry->rightAscension,"%.3f",dra2);
              sprintf(pMosaicTableEntry->declination,"%.2f",ddec2);
            } else {
              ra2str(pMosaicTableEntry->rightAscension,MAX_RIGHTASCENSION_STRING,dra2,3);
              dec2str(pMosaicTableEntry->declination,MAX_DECLINATION_STRING,ddec2,2);
            }
          }
        }

        /* output only one plate */
        if (totalSearch == 1 || singleColumnFlag != 0) {
          if (pMosaicTableEntry->plateNumber == oldPlateNumber) {
            continue;
          } else {
            oldPlateNumber = pMosaicTableEntry->plateNumber;
          }
        }

        plateNumberCopy2 = pMosaicTableEntry->plateNumber;
        tempResult = StoreResult(
          pMosaicTableEntry,
          seriesTablePtr,
          series,
          cra,
          cdec,
          dra,
          ddec,
          drad,
          tnum,
          tra,
          tdec,
          tmag,
          tclass,
          pSortTable,
          &sortTableSize,
          nstar,
          nstarmax,
          publicRelease,
          &selectWCSFit,
          &selectWCSSize,
          &selectWCSStrlen
        );

        if (tempResult) {
          nstar += tempResult;

          if (plateNumberCopy2 != oldPlateNumber2) {
            nplates++;
            oldPlateNumber2 = plateNumberCopy2;
          }
        }
      } /* End of mosaic selection loop */

      if (epoch1 > 0.0 && epoch2 > 0.0) {
        free(dateStr1);
        dateStr1 = NULL;
        free(dateStr2);
        dateStr2 = NULL;
      }

      if (pMosaicTable != NULL) {
        free(pMosaicTable);
        pMosaicTable = NULL;
      }
    }

    seriesTablePtr++;
  }

  /* Summarize search */
  if (query_tool_mode) {
    // No summary
  } else if (webtable == NULL) {
    if (nlog > 0) {
      fprintf(outfile, "SCANREAD: Catalog %s : %d / %d found\n", bincat, istar, sc->nstars);
      if (nstar > nstarmax)
        fprintf(outfile,"SCANREAD: %d stars found; only %d returned\n", nstar, nstarmax);
    }
  } else {
    nResult = nstar;

    if (nResult > nstarmax) {
      nResult = nstarmax;
    }

    if (nstar > 0) {
      printf("</table>\n");
      printf("<p>%d exposure(s) and %d plate(s) listed out of %d candidates.</p>", nResult, nplates, nstar);

      if (publicRelease != 0) {
        if (selectWCSFit > 0) {
          char *listString = calloc(selectWCSStrlen+10,1);
          int tempIndex;
          int curTempIndex = 0;

          if (listString) {
            qsort(pSortTable,sortTableSize,sizeof(SORTTABLE),dcmp);

            if (strstr(webarchivetype,"ZIP")) {
              strcpy(archiveString,"zip file");
            } else {
              strcpy(archiveString,"tarball");
            }

            printf(
              "<p>Create %s for a total of %d extracted images containing a total of %.3f megapixels:</p>\n",
              archiveString,
              selectWCSFit,
              selectWCSSize / 1000000.0
            );

            while (curTempIndex < sortTableSize) {
              printf("<form action=\"%s\" method=POST target=\"tarball\">\n", tarballextract);

              for (tempIndex = curTempIndex; tempIndex < sortTableSize; tempIndex++) {
                PSORTTABLE pSortEntry = &pSortTable[tempIndex];
                if (pSortEntry->listString != NULL) {
                  strcat(listString,pSortEntry->listString);
                }

                if (strlen(listString) > WEB_SIZE_LIMIT) {
                  tempIndex++;
                  break;
                }
              }

              if (strlen(listString) > selectWCSStrlen) {
                printf("Software error on string length %d\n",selectWCSStrlen);
              }

              printf("<input type=\"hidden\" name=\"filetype\" value = \"%s\" />\n", webfiletype);
              printf("<input type=\"hidden\" name=\"archivetype\" value = \"%s\"/>\n", webarchivetype);
              printf("<input type=\"hidden\" name=\"ra\" value = \"%f\"/>\n", cra);
              printf("<input type=\"hidden\" name=\"dec\" value = \"%f\"/>\n", cdec);
              printf("<input type=\"hidden\" name=\"system\" value = \"J2000\"/>\n");
              printf("<input type=\"hidden\" name=\"imagesize\" value = \"%.0f\"/>\n", 3600 * drad);
              printf("<input type=\"hidden\" name=\"sizeunits\" value = \"arcsec\"/>\n");
              printf("<input type=\"hidden\" name=\"selection\" value = \"%s\"/>\n", listString);
              printf("<input type=\"submit\" name=\"submit\" value=\"%s-%d-%d\">\n", archiveString, curTempIndex + 1, tempIndex);
              printf("</form>");
              curTempIndex = tempIndex;
              listString[0] = 0;
            }

            free(listString);
          }
        }
      }
    } else {
      printf("<p><strong>ERROR: No database entries found</strong></p>\n");
    }
  }

  if (pStoredResultTable != NULL) {
    free(pStoredResultTable);
  }
  if (pSortTable != NULL) {
    free(pSortTable);
  }

  if (pMosaicTable != NULL) {
    free(pMosaicTable);
  }

  return (nstar);
} /* End of scanread() */


/* SCANRNUM -- Read binary catalog stars with specified numbers */

int
scanrnum (bincat, nnum, sysout, eqout, epout, match,
          tnum,tra,tdec,tpra,tpdec,tmag,tpeak,tobj,nlog)

     char	*bincat;	/* Name of reference star catalog file */
     int	nnum;		/* Number of stars to look for */
     int	sysout;		/* Search coordinate system */
     double	eqout;		/* Search coordinate equinox */
     double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
     int	match;		/* If 1, match number exactly, else number is sequence*/
     double	*tnum;		/* Array of star numbers to look for */
     double	*tra;		/* Array of right ascensions (returned) */
     double	*tdec;		/* Array of declinations (returned) */
     double  *tpra;		/* Array of right ascension proper motions (returned) */
     double  *tpdec;		/* Array of declination proper motions (returned) */
     double	**tmag;		/* 2-D Array of magnitudes (returned) */
     int	*tpeak;		/* Array of peak counts (returned) */
     char	**tobj;		/* Array of object names (returned) */
     int	nlog;
{
  if (webtable) {
    outfile = stdout;
  } else {
    outfile = stderr;
  }

  return(0);
}


/* SCANBIN -- Fill FITS WCS image with stars from binary catalog */

int
scanbin (bincat, wcs, header, image, mag1, mag2, sortmag, magscale, nlog)

     char	*bincat;	/* Name of reference star catalog file */
     struct WorldCoor *wcs;	/* World coordinate system for image */
     char	*header;	/* FITS header for output image */
     char	*image;		/* Output FITS image */
     double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
     int	sortmag;	/* Magnitude by which to sort (1 to nmag) */
     double	magscale;	/* Scaling factor for magnitude to pixel flux
                       * (number of catalog objects per bin if 0) */
     int	nlog;
{
  if (webtable) {
    outfile = stdout;
  } else {
    outfile = stderr;
  }

  return (0);
}


/* ReadSeries -- Read in the series table */
static void *
ReadSeries(MYSQL * pConnection)
{
  PSERIESTABLE seriesTablePtr;
  PSERIESTABLE baseSeriesPtr = NULL;
  int numSeries;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char queryString[MAX_QUERY_STRING];
  int nvals;
  double plateScale;
  sprintf(queryString,"SELECT seriesId,series,nominalPlateScale,fittedPlateScale,orientation+0,sequestered+0 from series;");

  res = ExecuteQuery(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
      numSeries = (int)mysql_num_rows(res_ptr);
      if (numSeries == 0) {
        return(NULL);
      }
      seriesTablePtr = (PSERIESTABLE) calloc(numSeries+1,sizeof(SERIESTABLE));
      baseSeriesPtr = seriesTablePtr;
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        nvals = sscanf(sqlrow[0],"%d",&seriesTablePtr->seriesId);
        if (nvals != 1) {
          fprintf(outfile,"ERROR: nvals is %d for firstPlate\n", nvals);
          continue;
        }
        strcpy(seriesTablePtr->series,sqlrow[1]);
        if (sqlrow[2]) {
          nvals = sscanf(sqlrow[2],"%lf",&plateScale);
          if (nvals != 1) {
            fprintf(outfile,"ERROR: nvals is %d for nominalPlateScale\n", nvals);
          } else {
            seriesTablePtr->plateScale = plateScale;
          }

        }
        if (sqlrow[3]) {
          nvals = sscanf(sqlrow[3],"%lf",&plateScale);
          if (nvals != 1) {
            fprintf(outfile,"ERROR: nvals is %d for fittedPlateScale\n", nvals);
          } else {
            seriesTablePtr->plateScale = plateScale;
          }

        }
        if (sqlrow[4]) {
          nvals = sscanf(sqlrow[4],"%d",&seriesTablePtr->orientation);
          if (nvals != 1) {
            fprintf(outfile,"ERROR: nvals is %d for orientation\n", nvals);
          }
        }
        if (sqlrow[5]) {
          nvals = sscanf(sqlrow[5],"%d",&seriesTablePtr->sequestered);
          if (nvals != 1) {
            fprintf(outfile,"ERROR: nvals is %d for sequestered\n", nvals);
          }
        }
        seriesTablePtr++;
      }
      mysql_free_result(res_ptr);
    } else {
      fprintf(outfile,"%s: mysql_store_result failed\n",__FUNCTION__);
    }
  } else {
    fprintf(outfile,"%s: mysql_query failed\n",__FUNCTION__);
  }

  return(baseSeriesPtr);
}



/* SCANOPEN -- Open scanner catalog, returning number of entries */

struct StarCat *
scanopen (char *bincat /* Binary catalog file name */)
{
  struct StarCat *sc;
  MYSQL * pConnection;

  if (webtable) {
    outfile = stdout;
  } else {
    outfile = stderr;
  }

  sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));

  if (sc == NULL) {
    return(NULL);
  }
  if (sqlConnection != NULL) {
    printf("ERROR: Attempt to connect to the MySQL database twice in scanopen\n");
    exit(1);
  }

  pConnection = (MYSQL *)calloc(1,sizeof(MYSQL));
  if (pConnection == NULL) {
    scanclose(sc);
    return(NULL);
  }

  dasch_init_scandb(pConnection);

  sqlConnection = (void *)pConnection;

  sqlSeriesTable = ReadSeries(pConnection);
  if (sqlSeriesTable == NULL) {
    scanclose(sc);
    return(NULL);
  }

  sc->nstars = SCANCATSTARS;
  sc->inform = 'H';
  sc->coorsys = WCS_J2000;
  sc->epoch = 2000.0;
  sc->equinox = 2000.0;
  sc->nmag = 1;
  sc->stnum = 1;
  sc->nepoch = 1;
  strcpy(sc->isname,"Harvard Scanner MySQL Database");
  strcpy(sc->keymag[sc->nmag-1], "epoch");
  if (mosaicFlag) {
    sc->nnfld = 24; /* mosaic name: ssssnnnnn_mm_01rrrrww  where m is the mosaic and rrrr is the rotation */
  } else {
    sc->nnfld = 13; /* four characters for series, 5 for plate number, 1 space and 3 for exposure */
  }
  sc->stnum = 1;

  /* Set number of decimal places in star numbers. This will be the exposure number */
  sc->nndec = 3;



  sc->refcat = SCANCAT;
  return (sc);
}


void
scanclose (struct StarCat *sc /* Star catalog descriptor */)
{
  MYSQL * pConnection;

  if (sc == NULL) {
    return;
  }
  pConnection = sqlConnection;
  if (pConnection) {

    mysql_close(pConnection);
    free((void *)pConnection);
    sqlConnection = NULL;
  }

  if (sqlSeriesTable != NULL) {
    free(sqlSeriesTable);
    sqlSeriesTable = NULL;
  }

  if (sc->catline != NULL) {
    free ((void *)sc->catline);
  }
  if (sc != NULL) {
    free ((void *)sc);
  }
  return;
}


/* SCANSTAR -- Get binary catalog entry for one star;
   return 0 if successful */

int
scanstar (sc, st, istar)

     struct StarCat *sc;	/* Star catalog descriptor */
     struct Star *st;	/* Current star entry */
     int istar;	/* Star sequence number in binary catalog */
{
  if (webtable) {
    outfile = stdout;
  } else {
    outfile = stderr;
  }

  return (0);
}


void
setMosaicFlag(int flag)
{
  mosaicFlag = flag;
}


int
getMosaicFlag()
{
  return mosaicFlag;
}


void
setPessimisticFlag(int flag)
{
  pessimisticFlag = flag;
}


void
setSingleColumnFlag(int flag)
{
  singleColumnFlag = flag;
}


int
getSingleColumnFlag()
{
  return(singleColumnFlag);
}


void
setPlateSeries(char *strPtr)
{
  int length;

  if (searchPlateSeries) {
    free(searchPlateSeries);
    searchPlateSeries = NULL;
  }
  length = strlen(strPtr);
  if (length > 0) {
    searchPlateSeries = (char *)calloc(length+1,1);
  }
  if (searchPlateSeries) {
    strcpy(searchPlateSeries,strPtr);
    local_strlwr(searchPlateSeries); /* Fix of 2019-08-05 */
  }
}


void
setPlateNumber(char *strPtr)
{
  int length;

  if (searchPlateNumber) {
    free(searchPlateNumber);
    searchPlateNumber = NULL;
  }
  length = strlen(strPtr);
  if (length > 0) {
    searchPlateNumber = (char *)calloc(length+1,1);
  }
  if (searchPlateNumber) {
    strcpy(searchPlateNumber,strPtr);
  }
}


void
setPlateClass(char *strPtr)
{
  int length;

  if (searchPlateClass) {
    free(searchPlateClass);
    searchPlateClass = NULL;
  }
  length = strlen(strPtr);
  if (length > 0) {
    searchPlateClass = (char *)calloc(length+1,1);
  }
  if (searchPlateClass) {
    strcpy(searchPlateClass,strPtr);
  }
}


void
setSearchType(char *strPtr)
{
  searchType[0] = 0;
  if (strstr(strPtr,"scanned")) {
    strcpy(searchType,"scanned");
  } else if (strstr(strPtr,"pending")) {
    strcpy(searchType,"pending");
  } else if (strstr(strPtr,"wcsfit")) {
    strcpy(searchType,"wcsfit");
  } else if (strstr(strPtr,"all")) {
    strcpy(searchType,"all");
  } else if (strstr(strPtr,"total")) {
    strcpy(searchType,"total");
  } else {
    printf("ERROR scanread line %d searchType %s is undefined\n",__LINE__,strPtr);
  }
}


void
setWebTable(char *strPtr)
{
  int length;

  if (webtable) {
    free(webtable);
    webtable = NULL;
  }
  length = strlen(strPtr);
  if (length > 0) {
    webtable = (char *)calloc(length+1,1);
  }
  if (webtable) {
    strcpy(webtable,strPtr);
  }
}


void
setWebExtract(char *strPtr)
{
  int length;

  if (webextract) {
    free(webextract);
    webextract = NULL;
  }
  if (strstr(strPtr,"NULL") == NULL) {
    length = strlen(strPtr);
    if (length > 0) {
      webextract = (char *)calloc(length+1,1);
    }
    if (webextract) {
      strcpy(webextract,strPtr);
    }
  }
}


void
setTarballExtract(char *strPtr)
{
  int length;

  if (tarballextract) {
    free(tarballextract);
    tarballextract = NULL;
  }
  if (strstr(strPtr,"NULL") == NULL) {
    length = strlen(strPtr);
    if (length > 0) {
      tarballextract = (char *)calloc(length+1,1);
    }
    if (tarballextract) {
      strcpy(tarballextract,strPtr);
    }
  }
}


void
setWebFileType(char *strPtr)
{
  int length;

  if (webfiletype) {
    free(webfiletype);
    webfiletype = NULL;
  }
  length = strlen(strPtr);
  if (length > 0) {
    webfiletype = (char *)calloc(length+1,1);
  }
  if (webfiletype) {
    strcpy(webfiletype,strPtr);
  }
}


void
setWebArchiveType(char *strPtr)
{
  int length;

  if (webarchivetype) {
    free(webarchivetype);
    webarchivetype = NULL;
  }
  length = strlen(strPtr);
  if (length > 0) {
    webarchivetype = (char *)calloc(length+1,1);
  }
  if (webarchivetype) {
    strcpy(webarchivetype,strPtr);
  }
}


void
setRegionFlag(char *strPtr)
{
  int length;

  if (regionflag) {
    free(regionflag);
    regionflag = NULL;
  }
  length = strlen(strPtr);
  if (length > 0) {
    regionflag = (char *)calloc(length+1,1);
  }
  if (regionflag) {
    strcpy(regionflag,strPtr);
  }
}


void
ScanCatNum (struct StarCat *sc,int nnfld,int nndec,double dnum,char * numstr)
{
  PSERIESTABLE seriesTablePtr;
  int seriesNumber = (dnum + 0.5)/MAX_PLATE_NUMBER;
  int plateNumber = (dnum+0.5) - (1.0 * MAX_PLATE_NUMBER * seriesNumber);
  int rotation = dnum + 0.5;
  int exposureNumber; /* Can also be the mosaic number */

  rotation = ((dnum-rotation) * (MAX_EXPOSURE_NUMBER*MAX_ROTATION))+ 0.005;
  exposureNumber = rotation/MAX_ROTATION;
  rotation = (rotation - (exposureNumber*MAX_ROTATION))* 90;
  seriesTablePtr = sqlSeriesTable;

  if (seriesTablePtr != NULL) {
    while (seriesTablePtr->seriesId > 0) {
      if (seriesTablePtr->seriesId == seriesNumber) {
        if (mosaicFlag) {
          if (rotation == 450) {
            sprintf(numstr,"%4s%05d   %3d      ",seriesTablePtr->series,plateNumber,exposureNumber);
          } else if (rotation == 360) {
            sprintf(numstr,"%4s%05d_%02d_01      ",seriesTablePtr->series,plateNumber,exposureNumber);
          } else {
            if (rotation == 0) {
              sprintf(numstr,"%4s%05d_%02d_01ww    ",seriesTablePtr->series,plateNumber,exposureNumber);
            } else {
              if (rotation == 90) {
                sprintf(numstr,"%4s%05d_%02d_01r%dww ",seriesTablePtr->series,plateNumber,exposureNumber,rotation);
              } else {
                sprintf(numstr,"%4s%05d_%02d_01r%dww",seriesTablePtr->series,plateNumber,exposureNumber,rotation);
              }
            }
          }
        } else {
          sprintf(numstr,"%4s%05d %3d",seriesTablePtr->series,plateNumber,exposureNumber);
        }
        return;
      }
      seriesTablePtr++;
    }
  }

  sprintf(numstr,"%3d %5d %3d",seriesNumber,plateNumber,exposureNumber);
}


void
setSearchParameters(double date1,double date2)
{
  epoch1 = date1;
  epoch2 = date2;
}

void
dasch_scanread_set_query_tool_mode(int mode)
{
  query_tool_mode = mode;
}

/*
 * 2008-03-13: Divide tarballs into separate units to avoid exceeding the 8000 byte web request limit
 * 2008-03-26: Correct 90 degree rotation file names (r90ww instead of r090ww).
 * 2008-05-05: Give preference to mosaics that have been selected. i.e. FitWCS regexp("Selected")
 *             Correct mosaic ordering to select latest scan then latest mosaic.
 * 2008-05-12: Do not show deleted mosaics
 * 2008-05-14: Correct file name for unfitted mosaics.
 * 2008-06-06  Port to wcstools 3.7.3
 * 2008-08-04: If the webextract or tarballextract strings are NULL, do not save the string.
 * 2008-08-05: Restore LEFT JOIN to MySQL query removed on 2008-05-12
 * 2009-06-02: Use pop-up windows to avoid repeating the search
 * 2009-06-11: Reduce the halfwidth
 * 2009-11-03: Correct wraps near the equinox: use OR instead of AND on the right ascension
 *             Make the width of the search variable.
 * 2010-10-18: Convert tarball request from GET to POST
 * 2010-09-06: Fix Solaris build warnings
 * 2011-04-04: Add "pending" to search unscanned plates
 *             Clear pBestResult before query to avoid false positives
 * 2011-09-13  Add a new pessimistic search, so the center is in  bin 1
 * 2011-11-29  Sequester non-Harvard series.
 * 2011-12-13  Correct "pending" searches for the multiple exposure case
 * 2012-03-27  Sort the tarball string by Julian Date and prepend the year
 * 2012-06-18  Set illegal catalog declination 99.0 and RA to 999.0
 * 2012-12-17  Support decimal output of RA and Declination
 * 2013-03-30  Support data release
 * 2014-12-23  Remove SOLARIS_BUILD support
 * 2015-05-06  Make plate lists public in non-released areas (Josh's memorandum of 5/05/15 11:52)
 * 2018-07-23  Fix infinite loop when WEB_SIZE_LIMIT is exceeded
 * 2019-08-05  Make the series lowercase in setPlateSeries
 * 2019-10-04  Make "all" and "pending" search types accept deleted mosaics.  Reqires adding "deleted" flag to RESULTTABLE,
 * 2019-10-15  Define a "total" search to build total.list with the highest good mosaic number
 * 2019-11-20  Prevent StoreResult from storing more than requested number of stars
 * 2020-04-22  Rewrite scanread() to use a common pipelineutils interface and BuildMosaicList()
 *             Move RESULTTABLE to pipelineutils and rename it to MOSAICLIST
 *             Force default search type to be "wcsfit"
 *             Add setSingleColumnFlag()and getSingleColumn() flag to avoid using an awk script for the website "Show Plates" function
 *             Pass through 999 error flag for RA and 99 error flag for DEC
*/
