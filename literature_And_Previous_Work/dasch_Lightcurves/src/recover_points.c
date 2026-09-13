// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// Applying photometric calibrations.
//
// Inputs:
//
// - spatial bin file, `${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}.out.spatial_bins.db`
// - "maggrid" files, `${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_a${bin}_grid.tmp`
// - defect file, `${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_defect.db`
// - background file, `${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_background.db`
// - "drad" file, `${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_drad.db`
// - local cal / "dmagcor" backup, `${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_dmagcor_copy.db`
// - mag-dep-cal file, `${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_magdepcalibrate.db` (optional)
// - extinction file, `${DASCH_CATALOGALL}/${series}${id}_${scannum}_01${rot}ww${qual}_extinction.db`
// - SExtractor source list, `${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_tnx.db`
// - match file, `${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_tnx_u.db`
// - blend file, `${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_tnx_b.db`
// - "local" match file, `${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}.out.local.db`
// - scanner database info:
//   - => JD of solution,
//   - exposure number
//   - quality flags
//   - observatory location
//   - plate color
//
// Outputs:
//
// - calibrated table, `${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_allobjects.db`
// - local cal table, `${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_dmagcor.db`
//
// See the Git history for definitions of "RMS algorithms" A through D. We now hardcode "algorithm E".

#include <math.h>
#include <time.h>

// Starbase
#include <table.h>

// Mysql
#include <mysql.h>

// wcstools
#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "magdeputil.h"
#include "daschunistd.h"


#define MAX_BUFFER 512
#define RMSVALUE_INITIAL_ALLOC 1000 /* Initial allocation of the RMS table */
#define RMSCALIB_INCREMENT 20  /* Number of stars for which to calculate the RMS value */


typedef struct _spatial_bin {
  int spatial_bin;
  double limiting_mag;
  double max_bright_mag;
  double limiting_iso;
  double max_bright_iso;
  double colorterm;
  double errorcolor;
  int colorflag;
} SPATIAL_BIN, *PSPATIAL_BIN;

typedef struct _defect {
  int NUMBER;
  int AFLAGS;
} DEFECT, *PDEFECT;

typedef struct _rmsValue {
  double magcal_local;
  double diffval;
} RMSVALUE, *PRMSVALUE;

typedef struct _rmsCalib {
  double magcal_local;
  double rms;
} RMSCALIB, *PRMSCALIB;

typedef struct _dmagcor {
  int nx;             /* Total bins in width */
  int ny;             /* Total bins in height */
  int ix;             /* width bin */
  int iy;             /* height bin */
  int irec;           /* file index */
  double zout;        /* Magnitude correction */
  double errout;      /* RMS of magnitude correction */
  int npout;          /* Number of points used for correction */
  int rejectFlag;     /* Reason, if any, for rejecting this bin */
  double bright_magcor_local; /* magcor_local for brightest stars (obsolete) */
  double bright_magcor_error; /* magcal_local_error for brightest stars (obsolete) */
  int bright_npoints_local;   /* npoints_local for brightest stars (obsolete) */

  /* Drad calculations done by filterblended */
  int drad_bin_count;      /* Count in this bin (after clipping) */
  int drad_bin_size;       /* Size of this bin (expanded for insufficient points) */
  int drad_reject_count;   /* Stars rejected because of high drad*/
  double draMedian;
  double draRMS;
  double ddecMedian;
  double ddecRMS;
  /* Drad calculations done using good stars only */
  int drad_bin_count2;      /* Count in this bin (after clipping) */
  int drad_bin_size2;       /* Size of this bin (expanded for insufficient points) */
  int drad_reject_count2;   /* Stars rejected because of high drad*/
  double draMedian2;
  double draRMS2;
  double ddecMedian2;
  double ddecRMS2;
  double dradRMS2;

  struct _starimage *flink; /* Flink to next member of this spatial bin */
} DMAGCOR, *PDMAGCOR;

typedef struct _extinction {
  int enx;             /* Total bins in width */
  int eny;             /* Total bins in height */
  int eix;             /* width bin */
  int eiy;             /* height bin */
  double altitude;     /* altitude in degrees */
  double extinction;   /* extinction in magnitudes */
  double Date;         /* Exposure heliocentric Julian date */
  double timeAccuracy; /* Exposure time accuracy */
} EXTINCTION, *PEXTINCTION;

typedef struct _maggrid {
  double maggrid;        /* Stdmag */
  double isogrid;        /* MAG_ISO */
  double griderr;       /* Stdmag rms */
  int flaggrid;          /* Grid flag  0 = too bright
                          *            1 = good magnitude
                          *           -1 = too dim */
} MAGGRID, *PMAGGRID;

typedef struct _starimage {
  int NUMBER;         /* Sextractor reference number NOTE: Both _starimage and _blend must have NUMBER as the first entry */
  int AFLAGS;          /* flags word */
  int BFLAGS;          /* flags word */
  int selected;       /* This entry has been selected for output (begins with cal_local flag) */
  int outputFlag;     /* This entry is accepted for output */
  int npoints_local;  /* Number of points used for local calibration */
  int spatial_bin;    /* Lowess magnitude spatial bin  number */
  int local_bin;      /* Local calibration bin number */
  int magdep_bin;     /* Magnituded-dependent calibration bin (-1) if undefined */
  double X_IMAGE;        /* Sextractor X location in pixels */
  double Y_IMAGE;        /* Sextractor Y location in pixels */
  double MAG_ISO;     /* Sextractor isophotonic magnitude */
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
  double magcal_iso;  /* Lowess magnitude calibration */
  double magcal_iso_rms; /* Lowess magnitude calibration error */
  double magcal_local;   /* Local magnitude calibration */
  double magcal_local_error; /* Local magnitude calibration error */
  double magcal_magdep;  /* magnitude-corrected local magnitude */
  double magcal_magdep_rms; /* RMS of the magnitude-corrected local magnitude */
  double magcal_local_rms;   /* Total magnitude calibration error */
  double extinction;    /* Extinction correction */
  /* The following items are for GSC2.3.2 identified objects */
  double Stdmag;       /* GSC2.3.2 magnitude */
  double color;        /* GSC2.3.2 color */
  double dra;          /* arcsec error in right ascension */
  double ddec;         /* arcsec error in declination */
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  /* Additional fields added on April 17, 2008 */
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
  /* Blended magnitude for blended images */
  double Blendedmag;   /* Blended magnitude */
  struct _starimage *flink; /* Flink to next member of this spatial bin */
  int gsc_bin_index; /* gsc bin index */
  int irec;          /* dmag file index */
  /* Fields added Apr 16, 2013 */
  double heliocentricJD;
  double limiting_mag_local;
  double dradRMS2;
} STARIMAGE, *PSTARIMAGE;

typedef struct _blend {
  int NUMBER;          /* Sextractor reference number NOTE: Both _starimage and _blend must have NUMBER as the first entry */
  double Stdmag;       /* GSC2.3.2 magnitude */
  double Blendedmag;   /* Blended magnitude */
} BLEND, *PBLEND;

typedef struct _drad {
  int nx;             /* Total bins in width */
  int ny;             /* Total bins in height */
  int ix;             /* width bin */
  int iy;             /* height bin */
  /* Drad calculations done by filterblended */
  int drad_bin_count;      /* Count in this bin (after clipping) */
  int drad_bin_size;       /* Size of this bin (expanded for insufficient points) */
  int drad_reject_count;   /* Stars rejected because of high drad*/
  double draMedian;
  double draRMS;
  double ddecMedian;
  double ddecRMS;
} DRAD, *PDRAD;

extern GSCBIN gscBin64;

static PGSCBIN pGscBin = &gscBin64;


/* Sort routine based on the Sextractor reference number */
static int
ImageCompare(const void *first, const void *second)
{
  int numberFirst = ((PSTARIMAGE)first)->NUMBER;
  int numberSecond = ((PSTARIMAGE)second)->NUMBER;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  } else {
    return 0;
  }
}


/* Sort routine based on the GSC2.3.2 magnitude */
static int
RmsValueCompare(const void *first, const void *second)
{
  double valueFirst = ((PRMSVALUE)first)->magcal_local;
  double valueSecond = ((PRMSVALUE)second)->magcal_local;

  if (valueFirst > valueSecond) {
    return 1;
  } else if (valueFirst < valueSecond) {
    return -1;
  } else {
    return 0;
  }
}


// Table lookup for the rms calculation
static double
CalcRMS_E(
  int *binCalibCount,
  PRMSCALIB *pRmsCalibTable,
  int spatial_bin,
  double magcal_local
) {
  PRMSCALIB pHighRmsCalib;
  PRMSCALIB pLowRmsCalib = NULL;
  int index;

  if (binCalibCount[spatial_bin] <= 0) {
    return 99.0;
  }

  pLowRmsCalib = &(pRmsCalibTable[spatial_bin][0]);
  if (magcal_local < pLowRmsCalib->magcal_local) {
    return pLowRmsCalib->rms;
  }

  for (index = 1; index < binCalibCount[spatial_bin]; index++) {
    pHighRmsCalib = &(pRmsCalibTable[spatial_bin][index]);

    if (magcal_local > pLowRmsCalib->magcal_local && magcal_local <= pHighRmsCalib->magcal_local) {
      return pHighRmsCalib->rms;
    }

    pLowRmsCalib = pHighRmsCalib;
  }

  /* Off the top - no RMS */
  return 99.0;
}


// Given the MAG_ISO, return the magcal_iso and magcal_iso_rms
static int
CalculateLowessMagnitude(
  double MAG_ISO,
  double *magcal_iso,
  double *magcal_iso_rms,
  PMAGGRID maggrid_table,
  size_t maggrid_nrecs
) {
  int index;
  PMAGGRID pMaggrid1; /* Grid entry just before point */
  PMAGGRID pMaggrid2; /* Grid entry just beyond point */
  PMAGGRID pMaggrid3; /* Grid entry two steps before point because of annular9.m round-off error */
  int result = 0;

  *magcal_iso = 0;
  *magcal_iso_rms = 0;

  for (index = 1; index < maggrid_nrecs; index++) {
    pMaggrid2 = &maggrid_table[index];

    if (MAG_ISO < pMaggrid2->isogrid) {
      pMaggrid1 = &maggrid_table[index-1];

      if (index >= 2) {
        pMaggrid3 = &maggrid_table[index-2];
      } else {
        pMaggrid3 = pMaggrid2;
      }

      if (pMaggrid1->flaggrid == 1 || pMaggrid2->flaggrid == 1 || pMaggrid3->flaggrid == 1) {
        if (pMaggrid2->isogrid == pMaggrid1->isogrid) {
          *magcal_iso = (pMaggrid1->maggrid + pMaggrid2->maggrid) / 2;
        } else {
          *magcal_iso = pMaggrid1->maggrid + ((pMaggrid2->maggrid-pMaggrid1->maggrid) * (MAG_ISO - pMaggrid1->isogrid) / (pMaggrid2->isogrid-pMaggrid1->isogrid));
        }

        if (pMaggrid1->flaggrid == 1 && pMaggrid2->flaggrid == 1) {
          if (pMaggrid2->isogrid == pMaggrid1->isogrid) {
            *magcal_iso_rms = (pMaggrid1->griderr + pMaggrid2->griderr) / 2;
          } else {
            *magcal_iso_rms = pMaggrid1->griderr + ((pMaggrid2->griderr-pMaggrid1->griderr) * (MAG_ISO - pMaggrid1->isogrid) / (pMaggrid2->isogrid-pMaggrid1->isogrid));
          }
        } else if (pMaggrid1->flaggrid == 1) {
          *magcal_iso_rms = pMaggrid1->griderr;
        } else if (pMaggrid2->flaggrid == 1) {
          *magcal_iso_rms = pMaggrid2->griderr;
        } else {
          *magcal_iso_rms = pMaggrid3->griderr;
        }

        result = 1;
      }

      break;
    }
  }

  return result;
}


// The annular9.m procedure uses cubic interpolation which can occasionally blow
// up if there are a sparse number of points. Look for negative rms values and
// force them to be the no less than the last peak value
static void
CheckMaggridTable(
  PMAGGRID maggrid_table,
  size_t maggrid_nrecs,
  char *fileroot,
  int spatial_bin
) {
  size_t index;
  size_t index2;
  size_t negativeIndex = 0; /* First negative griderr encountered */
  size_t startIndex = 0; /* First griderr less than 90.0 encountered */
  size_t peakIndex = 0;
  PMAGGRID pMaggrid1;
  PMAGGRID pMaggrid2;
  PMAGGRID pMaggrid3;
  int errorFlag = 0;
  int startFlag = 0;
  double peakGriderr = 0;
  double dmagold = 0;
  double dmag;

  for (index = 0; index < maggrid_nrecs; index++) {
    pMaggrid1 = &maggrid_table[index];

    if (pMaggrid1->flaggrid == 1 && startFlag == 0) {
      startFlag = 1;
      startIndex = index;
    }

    if (index > 0) {
      /* Record the last peak of griderr */
      pMaggrid2 = &maggrid_table[index-1];

      if (pMaggrid1->flaggrid == 1 && pMaggrid2->flaggrid == 1) {
        dmag = pMaggrid1->griderr - pMaggrid2->griderr;

        if (dmag < 0 && dmagold >= 0) {
          peakGriderr = pMaggrid2->griderr;
          peakIndex = index-1;
        }

        dmagold = dmag;
      }
    }

    if (pMaggrid1->griderr < 0) {
      if (errorFlag == 0) {
        negativeIndex = index;
        errorFlag = 1;

        if (negativeIndex != startIndex) {
          pMaggrid3 = &maggrid_table[peakIndex];
          printf(
            "WARNING: CheckMaggridTable negative magcal_iso_rms %f peakGriderr %f at maggrid %f for spatial bin %d of %s\n",
            pMaggrid1->griderr,
            peakGriderr,
            pMaggrid3->maggrid,
            spatial_bin,
            fileroot
          );
        }
      }
    }
  }

  if (errorFlag == 1 && negativeIndex == startIndex) {
    /* Here the curve goes negative at the start of the trace */
    negativeIndex = maggrid_nrecs - 1; /* First negative griderr encountered */
    startIndex = maggrid_nrecs - 1; /* First griderr less than 90.0 encountered */
    startFlag = 0;
    peakGriderr = 0;
    peakIndex = maggrid_nrecs - 1;
    dmagold = 0;
    errorFlag = 0;

    for (index = 0; index < maggrid_nrecs; index++) {
      index2 = maggrid_nrecs-1-index;
      pMaggrid1 = &maggrid_table[index2];

      if (pMaggrid1->flaggrid == 1 && startFlag == 0) {
        startFlag = 1;
        startIndex = index2;
      }

      if (index2 < maggrid_nrecs - 1) {
        /* Record the first peak of griderr */
        pMaggrid2 = &maggrid_table[index2 + 1];

        if (pMaggrid1->flaggrid == 1 && pMaggrid2->flaggrid == 1) {
          dmag = pMaggrid1->griderr - pMaggrid2->griderr;

          if (dmag < 0 && dmagold >= 0) {
            peakGriderr = pMaggrid2->griderr;
            peakIndex = index2 - 1;
          }

          dmagold = dmag;
        }
      }

      if (pMaggrid1->griderr < 0) {
        if (errorFlag == 0) {
          negativeIndex = index2;
          errorFlag = 1;
          pMaggrid3 = &maggrid_table[peakIndex];
          printf(
            "WARNING: CheckMaggridTable negative magcal_iso_rms (2) %f peakGriderr %f at maggrid %f for spatial bin %d of %s\n",
            pMaggrid1->griderr,
            peakGriderr,
            pMaggrid3->maggrid,
            spatial_bin,
            fileroot
          );
        }
      }
    }

    if (errorFlag == 1) {
      /* Here clean up on the bright magnitude side */
      if (peakIndex < startIndex && startFlag == 1 && negativeIndex < peakIndex && peakGriderr > 0) {
        /* It is o.k. to fix things up by using a flat griderr set to the last peak */
        for (index = 0; index < maggrid_nrecs; index++) {
          index2 = maggrid_nrecs - 1 - index;
          pMaggrid1 = &maggrid_table[index2];

          if (pMaggrid1->flaggrid == 1 && index2 <= peakIndex) {
            pMaggrid1->griderr = peakGriderr;
          }
        }
      } else {
        /* Problem here, just set all negative values to 99.0 so as not to contaminate the database */
        pMaggrid1 = &maggrid_table[negativeIndex];

        printf(
          "ERROR: CheckMaggridTable negative magcal_iso_rms (2) %f not recoverable for "
          "spatial bin %d of %s startFlag %d peakGridError %f startIndex %zu peakIndex %zu negativeIndex %zu\n",
          pMaggrid1->griderr,
          spatial_bin,fileroot,
          startFlag,
          peakGriderr,
          startIndex,
          peakIndex,
          negativeIndex
        );

        for (index = 0; index < maggrid_nrecs; index++) {
          pMaggrid1 = &maggrid_table[index];

          if (pMaggrid1->flaggrid == 1 && pMaggrid1->griderr <= 0) {
            pMaggrid1->griderr = 99.0;
          }
        }
      }
    }
  } else {
    /* Here clean up on the dim magnitude side */
    if (errorFlag == 1) {
      if (peakIndex > startIndex && startFlag == 1 && negativeIndex > peakIndex && peakGriderr > 0) {
        /* It is o.k. to fix things up by using a flat griderr set to the last peak */
        for (index = 0; index < maggrid_nrecs; index++) {
          pMaggrid1 = &maggrid_table[index];

          if (pMaggrid1->flaggrid == 1 && index >= peakIndex) {
            pMaggrid1->griderr = peakGriderr;
          }
        }
      } else {
        /* Problem here, just set all negative values to 99.0 so as not to contaminate the database */
        pMaggrid1 = &maggrid_table[negativeIndex];

        printf(
          "ERROR: CheckMaggridTable negative magcal_iso_rms %f not recoverable for "
          "spatial bin %d of %s startFlag %d peakGridError %f startIndex %zu peakIndex %zu negativeIndex %zu\n",
          pMaggrid1->griderr,
          spatial_bin,fileroot,
          startFlag,
          peakGriderr,
          startIndex,
          peakIndex,
          negativeIndex
        );

        for (index = 0; index < maggrid_nrecs; index++) {
          pMaggrid1 = &maggrid_table[index];

          if (pMaggrid1->flaggrid == 1 && pMaggrid1->griderr <= 0) {
            pMaggrid1->griderr = 99.0;
          }
        }
      }
    }
  }
}


static void
ComputeDrad2(
  int nx,
  int ny,
  PDMAGCOR dmagcor_table,
  double maxDradArcsec,
  int *dradCount1,
  signed long long *pCurMemory,
  signed long long *pMaxMemory
) {
  DMAGCOR *pDmagcor;
  DMAGCOR *pDmagcor2 = NULL;
  STARIMAGE *pCurEntry;
  double *draVector = NULL;
  double *tmpdraVector;
  double *ddecVector = NULL;
  double *tmpddecVector;
  double *dradVector = NULL;
  double dradMedian2;
  int vectorAlloc = 0;
  int vectorCount;
  int ix;
  int iy;
  int irec;
  int draCount;
  int ddecCount;
  int dx;
  int dy;
  int iix;
  int iiy;

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      vectorCount = 0;
      irec = iy + (ix * ny);
      pDmagcor = &dmagcor_table[irec];
      pCurEntry = pDmagcor->flink;

      while (pCurEntry != NULL) {
        if (vectorCount >= vectorAlloc) {
          vectorAlloc += 100;
          tmpdraVector = realloc(draVector, vectorAlloc * sizeof(double));
          CalcMemory(1,100 * sizeof(double), pCurMemory, pMaxMemory);

          if (tmpdraVector == NULL) {
            printf("ERROR allocating tmpdravector of size %d\n",vectorAlloc);
            exit(1);
          }

          draVector = tmpdraVector;
          tmpdraVector = NULL;

          tmpddecVector = realloc(ddecVector, vectorAlloc * sizeof(double));
          CalcMemory(1,100 * sizeof(double), pCurMemory, pMaxMemory);

          if (tmpddecVector == NULL) {
            printf("ERROR allocating tmpddecvector of size %d\n",vectorAlloc);
            exit(1);
          }

          ddecVector = tmpddecVector;
          tmpddecVector = NULL;
        }

        draVector[vectorCount] = pCurEntry->dra;
        ddecVector[vectorCount] = pCurEntry->ddec;
        vectorCount++;
        pCurEntry = pCurEntry->flink;
      }

      draCount = CalcMedianAndRMS(vectorCount, STARS_PER_DMAGBIN, draVector, &pDmagcor->draMedian2, &pDmagcor->draRMS2, 0, 3.0, 0);
      ddecCount = CalcMedianAndRMS(vectorCount, STARS_PER_DMAGBIN, ddecVector, &pDmagcor->ddecMedian2, &pDmagcor->ddecRMS2, 0, 3.0, 0);

      if (draCount == 0 || ddecCount == 0) {
        /* We did not have sufficient stars here to calculate drad. Expand to all of the adjacent bins */
        vectorCount = 0;

        for (dx = -1; dx <= 1; dx++) {
          for (dy = -1; dy <= 1; dy++) {
            iix = dx + ix;
            if (iix < 0 || iix >= nx) {
              continue;
            }

            iiy = dy + iy;
            if (iiy < 0 || iiy >= ny) {
              continue;
            }

            pDmagcor->drad_bin_size2++;
            irec = iiy + (iix * ny);
            pDmagcor2 = &dmagcor_table[irec];
            pCurEntry = pDmagcor2->flink;

            while (pCurEntry != NULL) {
              if (vectorCount >= vectorAlloc) {
                vectorAlloc += 100;
                tmpdraVector = realloc(draVector, vectorAlloc * sizeof(double));
                CalcMemory(1, 100 * sizeof(double), pCurMemory, pMaxMemory);

                if (tmpdraVector == NULL) {
                  printf("ERROR allocating tmpdravector of size %d\n",vectorAlloc);
                  exit(1);
                }

                draVector = tmpdraVector;
                tmpdraVector = NULL;

                tmpddecVector = realloc(ddecVector, vectorAlloc * sizeof(double));
                CalcMemory(1,100 * sizeof(double), pCurMemory, pMaxMemory);

                if (tmpddecVector == NULL) {
                  printf("ERROR allocating tmpddecvector of size %d\n",vectorAlloc);
                  exit(1);
                }

                ddecVector = tmpddecVector;
                tmpddecVector = NULL;
              }

              draVector[vectorCount] = pCurEntry->dra;
              ddecVector[vectorCount] = pCurEntry->ddec;
              vectorCount++;
              pCurEntry = pCurEntry->flink;
            }
          }
        }

        draCount = CalcMedianAndRMS(vectorCount, 0, draVector, &pDmagcor->draMedian2, &pDmagcor->draRMS2, 0, 3.0, 0);
        ddecCount = CalcMedianAndRMS(vectorCount, 0, ddecVector, &pDmagcor->ddecMedian2, &pDmagcor->ddecRMS2, 0, 3.0, 0);

        if (draCount != 0 && ddecCount != 0) {
          /* Now find the bin median drad */
          dradVector = draVector;
          vectorCount = 0;

          for (dx = -1; dx <= 1; dx++) {
            for (dy = -1; dy <= 1; dy++) {
              iix = dx + ix;
              if (iix < 0 || iix >= nx) {
                continue;
              }

              iiy = dy + iy;
              if (iiy < 0 || iiy >= ny) {
                continue;
              }

              irec = iiy + (iix * ny);
              pDmagcor2 = &dmagcor_table[irec];
              pCurEntry = pDmagcor2->flink;

              while (pCurEntry != NULL) {
                if (vectorCount >= vectorAlloc) {
                  printf("ERROR: vectorCount %d exceeds vectorAlloc %d (2) in ComputeDrad2\n", vectorCount, vectorAlloc);
                  exit(1);
                }

                dradVector[vectorCount] = sqrt(sqr(pCurEntry->dra - pDmagcor->draMedian2) + sqr(pCurEntry->ddec - pDmagcor->ddecMedian2));
                vectorCount++;
                pCurEntry = pCurEntry->flink;
              }
            }
          }

          if (CalcMedianAndRMS(vectorCount, 0, dradVector, &dradMedian2, &pDmagcor->dradRMS2, 0, 3.0, 1) == 0) {
            pDmagcor->dradRMS2 = 99.0;
          }
        } else {
          pDmagcor->dradRMS2 = 99.0;
        }
      } else {
        pDmagcor->drad_bin_size2 = 1;

        /* Now find the bin median drad */
        dradVector = draVector;
        vectorCount = 0;
        pCurEntry = pDmagcor->flink;

        while (pCurEntry != NULL) {
          if (vectorCount >= vectorAlloc) {
            printf("ERROR: vectorCount %d exceeds vectorAlloc %d (1) in ComputeDrad2\n", vectorCount, vectorAlloc);
            exit(1);
          }

          dradVector[vectorCount] = sqrt(sqr(pCurEntry->dra - pDmagcor->draMedian2) + sqr(pCurEntry->ddec - pDmagcor->ddecMedian2));
          vectorCount++;
          pCurEntry = pCurEntry->flink;
        }

        if (CalcMedianAndRMS(vectorCount, 0, dradVector, &dradMedian2, &pDmagcor->dradRMS2, 0, 3.0, 1) == 0) {
          pDmagcor->dradRMS2 = 99.0;
        }
      }

      if (ddecCount == 0) {
        pDmagcor->ddecRMS2 = 99.0;
      }

      if (draCount == 0) {
        pDmagcor->draRMS2 = 99.0;
      }

      if (ddecCount < draCount) {
        pDmagcor->drad_bin_count2 = ddecCount;
      } else {
        pDmagcor->drad_bin_count2 = draCount;
      }
    }
  }

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      vectorCount = 0;
      irec = iy + (ix * ny);
      pDmagcor = &dmagcor_table[irec];

      /* If the local bin RMS is too high, flag everything in the bin as bad */
      if (pDmagcor->dradRMS2 > 90.0 || pDmagcor->dradRMS2 > maxDradArcsec) {
        pDmagcor->rejectFlag |= DMAG_REJECT_DRAD;
        pCurEntry = pDmagcor2->flink;

        while (pCurEntry != NULL) {
          /* This is a bad local smoothing bin.  Flag everything here as high drad */
          if (pCurEntry->AFLAGS < FILTER_AFLAG_BAD) {
            pDmagcor->drad_reject_count2++;
            *dradCount1 = *dradCount1 + 1;
          }

          pCurEntry->AFLAGS |= (1 << FILTER_AFLAG_DRADBIN);
          pCurEntry = pCurEntry->flink;
        }
      }

      /* Clean up the flinks in case will be reused later */
      pDmagcor->flink = NULL;
    }
  }

  free(draVector);
  CalcMemory(-1, vectorAlloc * sizeof(double), pCurMemory, pMaxMemory);
  free(ddecVector);
  CalcMemory(-1, vectorAlloc * sizeof(double), pCurMemory, pMaxMemory);
}


int
main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  char *ingestDirectory;
  char *matchDirectory;
  char *binDirectory;
  int numBins = 0;
  char outfile[MAX_BUFFER];
  char fileroot[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  FILE *outHandle = NULL;
  int printExtinctionError = 1;
  int printStdmagMatchError = 1;

  char magdep_name[MAX_BUFFER];
  int magdep_load_status;

  File spatial_bins_handle = NULL;
  char spatial_bins_name[MAX_BUFFER];
  PSPATIAL_BIN tmp_spatial_bin_table = NULL;
  TableHead spatial_bins_header = NULL;
  PSPATIAL_BIN spatial_bin_table = NULL;
  PSPATIAL_BIN pSpatial_bin;

  File maggrid_handle[MAX_SPATIAL_BINS+1];
  char tmp_name[MAX_BUFFER];
  char maggrid_name[MAX_BUFFER];
  TableHead maggrid_header[MAX_SPATIAL_BINS+1];
  size_t maggrid_nrecs[MAX_SPATIAL_BINS+1];
  PMAGGRID maggrid_table[MAX_SPATIAL_BINS+1];
  PMAGGRID pMaggrid;
  PMAGGRID pMaggrid2;

  File defect_handle = NULL;
  char defect_name[MAX_BUFFER];
  TableHead defect_header = NULL;
  PDEFECT defect_table = NULL;
  PDEFECT pDefect;
  PDEFECT pDefectPrev;
  size_t defect_nrecs;
  int defect_index;

  File background_handle = NULL;
  char background_name[MAX_BUFFER];
  TableHead background_header = NULL;
  PHIGHBACKGROUND background_table = NULL;
  PHIGHBACKGROUND pBackground;
  PHIGHBACKGROUND pBackgroundPrev;
  size_t background_nrecs;
  int background_index;

  File dmagcor_handle1 = NULL;
  FILE *dmagcor_handle2 = NULL;
  char dmagcor_name1[MAX_BUFFER];
  char dmagcor_name2[MAX_BUFFER];
  TableHead dmagcor_header = NULL;
  PDMAGCOR dmagcor_table = NULL;
  PDMAGCOR pDmagcor;

  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  PSTARIMAGE sextractor_table = NULL;
  size_t sextractor_nrecs = 0;
  int sextractor_index;
  PSTARIMAGE pSextractor;

  File match_handle = NULL;
  char match_name[MAX_BUFFER];
  TableHead match_header = NULL;
  PSTARIMAGE match_table = NULL;
  size_t match_nrecs = 0;
  int match_index;
  PSTARIMAGE pMatch;

  File blend_handle = NULL;
  char blend_name[MAX_BUFFER];
  TableHead blend_header = NULL;
  PBLEND blend_table = NULL;
  size_t blend_nrecs = 0;
  int blend_index;
  PBLEND pBlend;

  File drad_handle = NULL;
  char drad_name[MAX_BUFFER];
  TableHead drad_header = NULL;
  PDRAD drad_table = NULL;
  size_t drad_nrecs = 0;
  PDRAD pDrad;

  int maxSextractorNUMBER = 0;

  File local_handle = NULL;
  char local_name[MAX_BUFFER];
  TableHead local_header = NULL;
  PSTARIMAGE local_table = NULL;
  size_t local_nrecs = 0;
  int local_index;
  PSTARIMAGE pLocal;

  File extinction_handle = NULL;
  char extinction_name[MAX_BUFFER];
  TableHead extinction_header = NULL;
  PEXTINCTION extinction_table = NULL;
  size_t extinction_nrecs = 0;
  PEXTINCTION pExtinction;

  int x_dmagbins;
  int y_dmagbins;
  int nx;
  int ny;
  int ix;
  int iy;
  int local_bin;
  int enx;
  int eny;
  int eix;
  int eiy;
  int spatial_bin;
  int effective_spatial_bin;
  int irec = 0;
  size_t spatial_bin_nrecs = 0;
  size_t dmagcor_nrecs = 0;
  int maxLength = 0;
  int maxLengthIndex = 0;
  int curLength;
  int index;
  int index2;
  int verbose = 0;
  double geocentricJD = 0.0;
  double heliocentricJD = 0.0;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int illegalBinCount = 0;
  int missingBinCount = 0;
  int tooDimCount = 0;
  int tooBrightCount = 0;
  int tooBrightFlag = 0;
  int noLocalPoints = 0;
  int noRefCount = 0;
  int outputCount = 0;
  int singleBinCount = 0;
  int highZoutCount = 0;
  int dmagHiZoutCount = 0; /* Rejected because zout out of limits */
  int dmagMedianCount = 0; /* Rejected because median is too dim */
  int dmagBothCount = 0;   /* Rejected for both reasons */
  int dmagDradCount = 0; /* Rejected because of high Drad */
  int maxVectorCount = 0;
  int blendTotalCount = 0;
  int blendCorrectedCount = 0;
  int lowAltitudeCount = 0;
  int belowLimitingMagCount = 0;
  int goodLimitingMagRangeCount = 0;
  int secondLimitingMagRangeCount = 0;
  int allLimitingMagRangeCount = 0;
  int dradCount1 = 0;
  int dradCount2 = 0;
  int dradCount3 = 0;
  int dradCount4 = 0;

  double limiting_iso;
  double max_bright_iso;
  double magcal_iso;  /* Lowess magnitude calibration */
  double magcal_iso_rms; /* Lowess magnitude calibration error */
  double magcal_local; /* Local magnitude calibration */

  double min_X_IMAGE;
  double max_X_IMAGE;
  double min_Y_IMAGE;
  double max_Y_IMAGE;
  double limiting_mag_local;

  time_t startTime;
  time_t curTime;
  double edgeDist;

  double scale = -1.0; /* Plate scale in arcsec per pixel */
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int solutionNumber = -1;
  int binning;
  int rotation;
  double stdmagdiff;
  signed long long curMemory = 0;
  signed long long maxMemory = 0;
  int tmpMAGFlag;

  int gotDateColumn = 0;
  int gotTimeColumn = 0;
  int columnIndex;
  double extinctionFileTimeAccuracy;
  double extinctionFileDate;
  double singleJulianDate;
  double singleTolerance;
  double allJulianDate;
  double allTolerance;
  double dateDifference;
  double accuracyDifference;
  int numExposures = 0;
  int exposureNumber;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  int locationID;
  double latitude;
  double longitude;
  double elevation;
  int plateColor;
  double extinctionCoefficient;
  double maxAltitude = 0.0;
  double minExtinction = 0.0;
  double z1;
  double extinctionDifference;
  double dradlimit; /* Drad limit in arcsec */
  int printTimeAccuracy = 1;
  int catalogNumber =  CATALOG_GSC232;
  int tmp_spatial_bin;

  PMAGDEPCORRECTION pMagdepTable = NULL;
  MAGDEPLIMITS magdepLimits;
  PMAGDEPLIMITS pMagdepLimits = &magdepLimits;
  int temp_magdep_bin;
  double temp_magdep_rms;
  double temp_magdep_magcor;
  char *charPtr;
  char algorithm[] = "E";
  int rmsValueAlloc[MAX_SPATIAL_BINS+1];
  int rmsValueCount[MAX_SPATIAL_BINS+1];
  int binCalibCount[MAX_SPATIAL_BINS+1];
  PRMSVALUE pRmsValueTable[MAX_SPATIAL_BINS+1];
  PRMSCALIB pRmsCalibTable[MAX_SPATIAL_BINS+1];
  PRMSCALIB pRmsCalib;
  PRMSVALUE pRmsValue;
  int rmsCalibIndex;
  int quality = 0;  /* QUALITY_* plate quality flags */
  int filterMaskRMSReject = FILTER_MASK_RMS_REJECT;

  time(&startTime);

  local_name[0] = 0;
  match_name[0] = 0;
  blend_name[0] = 0;
  sextractor_name[0] = 0;
  extinction_name[0] = 0;
  outfile[0] = 0;
  fileroot[0] = 0;
  qualifier[0] = 0;

  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    maggrid_handle[spatial_bin] = NULL;
    maggrid_header[spatial_bin] = NULL;
    maggrid_table[spatial_bin] = NULL;
    maggrid_nrecs[spatial_bin] = 0;
  }

  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    printf("ERROR: DASCH_INGEST is not defined\n");
    return 1;
  }

  matchDirectory = getenv("DASCH_MATCH");
  if (ingestDirectory == NULL) {
    printf("ERROR: DASCH_MATCH is not defined\n");
    return 1;
  }

  binDirectory = getenv("DASCH_BINOUTPUT");
  if (binDirectory == NULL) {
    printf("ERROR: DASCH_BINOUTPUT is not defined\n");
    return 1;
  }

  // Handle arguments

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: argument %s does not have a qualifier\n", argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'q': /* file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(qualifier,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
              errorFlag = 1;
            }
            if (strstr(qualifier,"kepler") != NULL) {
              catalogNumber = CATALOG_KEPLER;
            }
            if (strstr(qualifier,"apass") != NULL) {
              catalogNumber = CATALOG_APASS;
            }
            if (strstr(qualifier,"gaia") != NULL) {
              catalogNumber = CATALOG_GAIA;
            }
            if (strstr(qualifier,"atlas") != NULL) {
              catalogNumber = CATALOG_ATLAS;
            }
            if (strstr(qualifier,"experimental") != NULL) {
              catalogNumber = CATALOG_EXPERIMENTAL;
            }

          }
          break;

        case 'l': /* local file name */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(local_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'e': /* extinction output file name */
        case 'E':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(extinction_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'm': /* match file name */
        case 'M':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(match_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'n': /* no match file name */
        case 'N':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            // This input is no longer used.
            argv++;
          }
          break;

        case 'd': /* blend file name */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(blend_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'r': /* Plate scale */
        case 'R':
          argc--;
          nvals = sscanf(*++argv,"%lf",&scale);
          if (nvals != 1) {
            printf("ERROR: Can not decode the plate scale\n");
            errorFlag = 1;
          }
          break;

        case 's': /* sextractor file name */
        case 'S':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sextractor_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'p': /* file root */
        case 'P':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(fileroot,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'b': /* Number of bins */
        case 'B':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&numBins);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the number of bins %s\n",*argv);
              errorFlag = 1;
            } else {
              if ((numBins != 1) && (numBins != 9)) {
                printf("ERROR: Unrecognized DASCH_NUMBINS value %d\n",numBins);
                errorFlag = 1;
              }
              if (numBins == 1) {
                printf("recover_points single bin case\n");
              }
            }
          }
          break;

        case 'c': /* Catalog RMS. This is in the form "catrms=<rms>" */
        case 'C':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            // This argument is no longer used.
            argv++;
          }
          break;

        case 'a': /* solution number */
        case 'A':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&solutionNumber);
            if (nvals != 1) {
              printf("ERROR: Unable to decode solution Number %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic width %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* mosaic height */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  // Verify that we have everything

  if (numBins == 0) {
    printf("ERROR: No number of bins specified\n");
    errorFlag = 1;
  }

  if (solutionNumber < 0) {
    printf("ERROR: No solution number specified\n");
    errorFlag = 1;
  }

  if (mosaicWidth == 0) {
    printf("ERROR: No mosaic width specified\n");
    errorFlag = 1;
  }

  if (mosaicHeight == 0) {
    printf("ERROR: No mosaic height specified\n");
    errorFlag = 1;
  }

  x_dmagbins = XDmagBins(mosaicWidth, mosaicHeight);
  y_dmagbins = YDmagBins(mosaicWidth, mosaicHeight);

  if (outfile[0] == 0) {
    printf("ERROR: No output filename was specified\n");
    errorFlag = 1;
  }

  if (local_name[0] == 0) {
    printf("ERROR: No local filename was specified\n");
    errorFlag = 1;
  }

  if (match_name[0] == 0) {
    printf("ERROR: No match filename was specified\n");
    errorFlag = 1;
  }

  if (blend_name[0] == 0) {
    printf("ERROR: No blend filename was specified\n");
    errorFlag = 1;
  }

  if (extinction_name[0] == 0) {
    printf("ERROR: No extinction filename was specified\n");
    errorFlag = 1;
  }

  if (sextractor_name[0] == 0) {
    printf("ERROR: No sextractor filename was specified\n");
    errorFlag = 1;
  }

  if (scale <= 0.0) {
    printf("ERROR: No plate scale\n");
    errorFlag = 1;
  }

  if (fileroot[0] == 0) {
    printf("ERROR: No file root was specified\n");
    errorFlag = 1;
  }

  // Ready to get going! Open the spatial bin file.

  strcpy(spatial_bins_name, ingestDirectory);
  strcat(spatial_bins_name, "/");
  strcat(spatial_bins_name, fileroot);
  strcat(spatial_bins_name, qualifier);
  strcat(spatial_bins_name, ".out.spatial_bins.db");

  spatial_bins_handle = Open(spatial_bins_name, "r");

  if (spatial_bins_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find spatial bins file %s\n",spatial_bins_name);
  } else {
    if (verbose) {
      printf("Found spatial bins file %s\n",spatial_bins_name);
    }
  }

  // Defect file

  strcpy(defect_name, matchDirectory);
  strcat(defect_name, "/");
  strcat(defect_name, fileroot);

  if (strstr(qualifier,"kepler")) {
    strcat(defect_name,qualifier);
  }

  if (strstr(qualifier,"apass")) {
    strcat(defect_name,qualifier);
  }

  if (strstr(qualifier,"gaia")) {
    strcat(defect_name,qualifier);
  }

  if (strstr(qualifier,"atlas")) {
    strcat(defect_name,qualifier);
  }

  if (strstr(qualifier,"experimental")) {
    strcat(defect_name,qualifier);
  }

  strcat(defect_name,"_defect.db");

  defect_handle = Open(defect_name,"r");

  if (defect_handle == NULL) {
    printf("No defect file found  %s\n",defect_name);
  } else {
    if (verbose) {
      printf("Found plate defect file %s\n",defect_name);
    }
  }

  // Background file

  strcpy(background_name, matchDirectory);
  strcat(background_name, "/");
  strcat(background_name, fileroot);

  /* background files apply to every exposure in a multiple exposure set */
  charPtr = strstr(background_name,"_s");
  if (charPtr != NULL) {
    *charPtr = 0;
  }

  strcat(background_name,"_background.db");
  background_handle = Open(background_name,"r");

  if (background_handle == NULL) {
    printf("ERROR: No background file found  %s\n",background_name);
  } else {
    if (verbose) {
      printf("Found plate background file %s\n",background_name);
    }
  }

  // drad file

  strcpy(drad_name, matchDirectory);
  strcat(drad_name, "/");
  strcat(drad_name, fileroot);

  if (strstr(qualifier,"kepler")) {
    strcat(drad_name,qualifier);
  }

  if (strstr(qualifier,"apass")) {
    strcat(drad_name,qualifier);
  }

  if (strstr(qualifier,"gaia")) {
    strcat(drad_name,qualifier);
  }

  if (strstr(qualifier,"atlas")) {
    strcat(drad_name,qualifier);
  }

  if (strstr(qualifier,"experimental")) {
    strcat(drad_name,qualifier);
  }

  strcat(drad_name,"_drad.db");

  drad_handle = Open(drad_name,"r");

  if (drad_handle == NULL) {
    printf("ERROR: No drad file found  %s\n",drad_name);
    errorFlag = 1;
  } else {
    if (verbose) {
      printf("Found plate drad file %s\n",drad_name);
    }
  }

  // Local calibration files

  strcpy(dmagcor_name1, ingestDirectory);
  strcat(dmagcor_name1, "/");
  strcat(dmagcor_name1, fileroot);
  strcat(dmagcor_name1, qualifier);
  strcat(dmagcor_name1, "_dmagcor_copy.db");

  strcpy(dmagcor_name2, ingestDirectory);
  strcat(dmagcor_name2, "/");
  strcat(dmagcor_name2, fileroot);
  strcat(dmagcor_name2, qualifier);
  strcat(dmagcor_name2, "_dmagcor.db");

  dmagcor_handle1 = Open(dmagcor_name1,"r");

  if (dmagcor_handle1 == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the local calibration file %s\n",dmagcor_name1);
  } else {
    if (verbose) {
      printf("Found local calibration file %s\n",dmagcor_name1);
    }
  }

  // Extinction file

  extinction_handle = Open(extinction_name,"r");

  if (extinction_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the local calibration file %s\n",extinction_name);
  } else {
    if (verbose) {
      printf("Found local calibration file %s\n",extinction_name);
    }
  }

  // SExtractor source list file

  sextractor_handle = Open(sextractor_name,"r");

  if (sextractor_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the sextractor file %s\n",sextractor_name);
  } else {
    if (verbose) {
      printf("Found sextractor file %s\n",sextractor_name);
    }
  }

  // Match file

  match_handle = Open(match_name,"r");

  if (match_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the match file %s\n",match_name);
  } else {
    if (verbose) {
      printf("Found match file %s\n",match_name);
    }
  }

  // Blend file

  blend_handle = Open(blend_name,"r");

  if (blend_handle == NULL) {
    printf("WARNING: Failed to find the blend file %s\n",blend_name);
    blend_nrecs = 0;
  } else {
    if (verbose) {
      printf("Found blend file %s\n",blend_name);
    }
  }

  // "local" matched sources file

  local_handle = Open(local_name,"r");
  if (local_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the local file %s\n",local_name);
  } else {
    if (verbose) {
      printf("Found local file %s\n",local_name);
    }
  }

  // Parse output name for plate ID info that we need for the DB lookup

  if (ParseFilename(outfile, series, &plateNumber, &mosaicNumber, &binning, &rotation) == 0) {
    printf("ERROR: Failed to parse %s\n",outfile);
    errorFlag = 1;
  }

  if (errorFlag) {
    printf("Usage: recover_points -p <file root> \n");
    printf("                      -w <mosaic width in pixels> \n");
    printf("                      -h <mosaic height in pixels> \n");
    printf("                      -l <local calibration file> \n");
    printf("                      -m <match file> \n");
    printf("                      -n <no-match file> \n");
    printf("                      -d <blend file name>\n");
    printf("                      -b <Number of bins>\n");
    printf("                      -e <extinction file> \n");
    printf("                      -s <sextractor file> \n");
    printf("                      -o <output file>\n");
    printf("                      -c catrms=<rms> where rms is the catalog error in magnitudes\n");
    printf("                      -r <plate scale arcsec/pixel\n");
    printf("                      -a <solution number>\n");
    printf("                      -v verbose\n");
    return 1;
  }

  // So far so good ...

  printf("recover_points of %s %s RMS algorithm '%s'\n", __DATE__, __TIME__, algorithm);

  dradlimit = MaxDradPixels("none", plateNumber) * scale;

  dasch_init_scandb(pConnection);

  // Read the spatial bins data

  spatial_bins_header = table_header(spatial_bins_handle,TABLE_PARSE);
  if (spatial_bins_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",spatial_bins_name);
    return 1;
  }

  tmp_spatial_bin_table = table_loadva(
    spatial_bins_handle,
    &spatial_bins_header,
    NULL, /* hbase */
    NULL, /* rows */
    NULL,
    sizeof(SPATIAL_BIN),
    &spatial_bin_nrecs,
    TblInt, "spatial_bin", TblOff(PSPATIAL_BIN, spatial_bin),
    TblDbl, "limiting_mag", TblOff(PSPATIAL_BIN, limiting_mag),
    TblDbl, "max_bright_mag", TblOff(PSPATIAL_BIN, max_bright_mag),
    TblDbl, "limiting_iso", TblOff(PSPATIAL_BIN, limiting_iso),
    TblDbl, "max_bright_iso", TblOff(PSPATIAL_BIN, max_bright_iso),
    TblDbl, "colorterm", TblOff(PSPATIAL_BIN, colorterm),
    TblDbl, "errorcolor", TblOff(PSPATIAL_BIN, errorcolor),
    TblInt, "colorflag", TblOff(PSPATIAL_BIN, colorflag),
    0,"end",0
  );

  if (tmp_spatial_bin_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",spatial_bins_name);
    return 1;
  }

  CalcMemory(1, spatial_bin_nrecs * sizeof(SPATIAL_BIN), &curMemory, &maxMemory);

  if (verbose) {
    printf("read %zu records for %s\n", spatial_bin_nrecs, spatial_bins_name);
  }

  // Copy the spatial bin entries into their proper slots

  spatial_bin_table = (PSPATIAL_BIN) calloc(MAX_SPATIAL_BINS + 1, sizeof(SPATIAL_BIN));
  CalcMemory(1, (MAX_SPATIAL_BINS + 1) * sizeof(SPATIAL_BIN), &curMemory, &maxMemory);
  if (spatial_bin_table == NULL) {
    printf("ERROR: Failed to allocate the spatial bin table \n");
    return 1;
  }

  for (index = 0; index < spatial_bin_nrecs; index++) {
    spatial_bin = tmp_spatial_bin_table[index].spatial_bin;

    if (spatial_bin > 0 && spatial_bin <= MAX_SPATIAL_BINS) {
      memcpy(&spatial_bin_table[spatial_bin], &tmp_spatial_bin_table[index], sizeof(SPATIAL_BIN));

      // Open the associated "grid" file and read it in
      strcpy(maggrid_name, ingestDirectory);
      strcat(maggrid_name, "/");
      strcat(maggrid_name, fileroot);
      strcat(maggrid_name, qualifier);
      sprintf(tmp_name, "_a%d_grid.tmp", spatial_bin);
      strcat(maggrid_name, tmp_name);

      maggrid_handle[spatial_bin] = Open(maggrid_name,"r");
      if (maggrid_handle[spatial_bin] == NULL) {
        printf("ERROR: Failed to find magnitude grid  file %s\n",maggrid_name);
        return 1;
      } else {
        if (verbose) {
          printf("Found magnitude grid file %s\n",maggrid_name);
        }
      }

      maggrid_header[spatial_bin] = table_header(maggrid_handle[spatial_bin], TABLE_PARSE);
      if (maggrid_header[spatial_bin] == NULL) {
        printf("ERROR: Failed to read header for %s\n",maggrid_name);
        return 1;
      }

      maggrid_table[spatial_bin] = table_loadva(
        maggrid_handle[spatial_bin],
        &maggrid_header[spatial_bin],
        NULL, /* hbase */
        NULL, /* rows */
        NULL,
        sizeof(MAGGRID),
        &maggrid_nrecs[spatial_bin],
        TblDbl, "maggrid", TblOff(PMAGGRID, maggrid),
        TblDbl, "isogrid", TblOff(PMAGGRID, isogrid),
        TblDbl, "griderr", TblOff(PMAGGRID, griderr),
        TblInt, "flaggrid", TblOff(PMAGGRID, flaggrid),
        0, "end", 0
      );

      if (maggrid_table[spatial_bin] == NULL) {
        printf("ERROR: Failed to read table for %s\n",maggrid_name);
        return 1;
      }

      CheckMaggridTable(maggrid_table[spatial_bin], maggrid_nrecs[spatial_bin], fileroot, spatial_bin);
      CalcMemory(1, maggrid_nrecs[spatial_bin] * sizeof(MAGGRID), &curMemory, &maxMemory);

      if (verbose) {
        printf("read %zu records for %s\n",maggrid_nrecs[spatial_bin],maggrid_name);
      }

      // Validate the table, looking for out-of order fields

      for (index2 = 1; index2 < maggrid_nrecs[spatial_bin]; index2++) {
        pMaggrid2 = &maggrid_table[spatial_bin][index2];
        pMaggrid = &maggrid_table[spatial_bin][index2 - 1];

        if (pMaggrid2->maggrid < pMaggrid->maggrid) {
          printf("ERROR: maggrid %f %f out of order for entry %d of %s\n", pMaggrid->maggrid, pMaggrid2->maggrid, index2, maggrid_name);
        }

        if (pMaggrid2->isogrid < pMaggrid->isogrid) {
          if ((pMaggrid->isogrid - pMaggrid2->isogrid) > 0.03) {
            printf("ERROR: isogrid %f %f out of order for entry %d of %s\n", pMaggrid->isogrid, pMaggrid2->isogrid, index2, maggrid_name);
          } else {
            printf("WARNING: isogrid %f %f out of order for entry %d of %s\n", pMaggrid->isogrid, pMaggrid2->isogrid, index2, maggrid_name);
          }
        }
      }
    } else {
      printf("ERROR: Illegal spatial bin %d\n", spatial_bin);
    }
  }

  if (tmp_spatial_bin_table != NULL) {
    Free(tmp_spatial_bin_table);
  }

  CalcMemory(-1, spatial_bin_nrecs * sizeof(SPATIAL_BIN), &curMemory, &maxMemory);

  // Read the drad file

  drad_header = table_header(drad_handle,TABLE_PARSE);
  if (drad_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",drad_name);
    return 1;
  }

  drad_table = table_loadva(
    drad_handle,
    &drad_header,
    NULL, /* hbase */
    NULL, /* rows */
    NULL,
    sizeof(DRAD),
    &drad_nrecs,
    TblInt, "nx", TblOff(PDRAD, nx),
    TblInt, "ny", TblOff(PDRAD, ny),
    TblInt, "ix", TblOff(PDRAD, ix),
    TblInt, "iy", TblOff(PDRAD, iy),
    TblInt, "drad_bin_count", TblOff(PDRAD, drad_bin_count),
    TblInt, "drad_bin_size", TblOff(PDRAD, drad_bin_size),
    TblInt, "drad_reject_count", TblOff(PDRAD, drad_reject_count),
    TblDbl, "draMedian", TblOff(PDRAD, draMedian),
    TblDbl, "draRMS", TblOff(PDRAD, draRMS),
    TblDbl, "ddecMedian", TblOff(PDRAD, ddecMedian),
    TblDbl, "ddecRMS", TblOff(PDRAD, ddecRMS),
    0, "end", 0
  );

  if (drad_table == NULL) {
    printf("ERROR: Failed to read table for %s\n", drad_name);
    return 1;
  }

  CalcMemory(1, drad_nrecs * sizeof(DRAD), &curMemory, &maxMemory);

  if (verbose) {
    printf("read %zu records for %s\n", drad_nrecs, drad_name);
  }

  // validate the table
  if (drad_table[0].nx != x_dmagbins || drad_table[0].ny != y_dmagbins || drad_nrecs != (x_dmagbins*y_dmagbins)) {
    printf(
      "ERROR: drad table has wrong size nx: %d %d, ny %d %d, total %zu %d\n",
      drad_table[0].nx,
      x_dmagbins,
      drad_table[0].ny,
      y_dmagbins,
      drad_nrecs,
      x_dmagbins * y_dmagbins
    );
    exit(1);
  }

  nx = x_dmagbins;
  ny = y_dmagbins;
  irec = 0;

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      pDrad = &drad_table[irec];

      if (pDrad->nx != nx && pDrad->ny != ny && pDrad->ix != ix && pDrad->iy != iy) {
        printf(
          "ERROR Mismatch of nx %d %d, ny %d %d, ix %d %d, or iy %d %d\n",
          pDrad->nx,
          nx,
          pDrad->ny,
          ny,
          pDrad->ix,
          ix,
          pDrad->iy,
          iy
        );
        return 1;
      }

      irec++;
    }
  }

  // Read in the local calibration file

  dmagcor_header = table_header(dmagcor_handle1, TABLE_PARSE);
  if (dmagcor_header == NULL) {
    printf("ERROR: Failed to read header for %s\n", dmagcor_name1);
    return 1;
  }

  dmagcor_table = table_loadva(
    dmagcor_handle1,
    &dmagcor_header,
    NULL, /* hbase */
    NULL, /* rows */
    NULL,
    sizeof(DMAGCOR),
    &dmagcor_nrecs,
    TblInt, "nx", TblOff(PDMAGCOR, nx),
    TblInt, "ny", TblOff(PDMAGCOR, ny),
    TblInt, "ix", TblOff(PDMAGCOR, ix),
    TblInt, "iy", TblOff(PDMAGCOR, iy),
    TblInt, "npout", TblOff(PDMAGCOR, npout),
    TblDbl, "zout", TblOff(PDMAGCOR, zout),
    TblDbl, "errout", TblOff(PDMAGCOR, errout),
    TblDbl, "bright_magcor_local", TblOff(PDMAGCOR, bright_magcor_local),
    TblDbl, "bright_magcor_error", TblOff(PDMAGCOR, bright_magcor_error),
    TblInt, "bright_npoints_local", TblOff(PDMAGCOR, bright_npoints_local),
    0,"end",0
  );

  if (dmagcor_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",dmagcor_name1);
    return 1;
  }

  CalcMemory(1, dmagcor_nrecs * sizeof(DMAGCOR), &curMemory, &maxMemory);

  if (verbose) {
    printf("read %zu records for %s\n", dmagcor_nrecs, dmagcor_name1);
  }

  // Validate
  nx = dmagcor_table[0].nx;
  ny = dmagcor_table[0].ny;
  if (dmagcor_nrecs != nx * ny) {
    printf("ERROR: Number of records %zu does not match product of %d and %d\n", dmagcor_nrecs, nx, ny);
    return 1;
  }

  if (nx != x_dmagbins || ny != y_dmagbins) {
    printf("ERROR: Mismatch in local bins between drad %d %d and dmagcor %d %d tables\n", x_dmagbins, y_dmagbins, nx, ny);
    exit(1);
  }

  irec = 0;

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      // See if we need to flag all stars in this bin because of excessively high or low zout
      pDmagcor = &dmagcor_table[irec];
      pDmagcor->rejectFlag = 0;
      pDmagcor->irec = irec;
      pDmagcor->flink = NULL;

      if (pDmagcor->zout < MINIMUM_ZOUT || pDmagcor->zout > MAXIMUM_ZOUT || pDmagcor->errout > MAXIMUM_ERROUT) {
        pDmagcor->rejectFlag |= DMAG_REJECT_HIZOUT;
        dmagHiZoutCount++;
      }

      // Combine the drad and dmagcor tables
      pDrad = &drad_table[irec];

      if (pDrad->ix != pDmagcor->ix || pDrad->iy != pDmagcor->iy) {
        printf(
          "ERROR: pDrad ix %d iy %d does not match pDmagcor ix %d iy %d\n",
          pDrad->ix,
          pDrad->iy,
          pDmagcor->ix,
          pDmagcor->iy
        );
        exit(1);
      }

      pDmagcor->drad_bin_count = pDrad->drad_bin_count;
      pDmagcor->drad_bin_size = pDrad->drad_bin_size;
      pDmagcor->drad_reject_count = pDrad->drad_reject_count;
      pDmagcor->draMedian = pDrad->draMedian;
      pDmagcor->draRMS = pDrad->draRMS;
      pDmagcor->ddecMedian = pDrad->ddecMedian;
      pDmagcor->ddecRMS = pDrad->ddecRMS;
      pDmagcor->drad_bin_count2 = 0;
      pDmagcor->drad_bin_size2 = 0;
      pDmagcor->drad_reject_count2 = 0;
      pDmagcor->draMedian2 = 0;
      pDmagcor->draRMS2 = 0;
      pDmagcor->ddecMedian2 = 0;
      pDmagcor->ddecRMS2 = 0;

      if (pDmagcor->nx != nx && pDmagcor->ny != ny && pDmagcor->ix != ix && pDmagcor->iy != iy) {
        printf(
          "ERROR Mismatch of nx %d %d, ny %d %d, ix %d %d, or iy %d %d\n",
          pDmagcor->nx,
          nx,
          pDmagcor->ny,
          ny,
          pDmagcor->ix,
          ix,
          pDmagcor->iy,
          iy
        );
        return 1;
      }

      irec++;
    }
  }

  // Done with drad

  if (drad_table != NULL) {
    Free(drad_table);
    drad_table = NULL;
  }

  CalcMemory(-1, drad_nrecs * sizeof(DRAD), &curMemory, &maxMemory);

  if (drad_header != NULL) {
    table_hdrfree(drad_header);
    drad_header = NULL;
  }

  if (drad_handle != NULL) {
    Close(drad_handle);
    drad_handle = NULL;
  }

  // Read in the extinction file

  extinction_header = table_header(extinction_handle, TABLE_PARSE);
  if (extinction_header == NULL) {
    printf("ERROR: Failed to read header for %s\n", extinction_name);
    return 1;
  }

  for (columnIndex = 1; columnIndex <= extinction_header->header->ncol; columnIndex++) {
    if (strcmp(extinction_header->header->column[columnIndex], "Date") == 0) {
      gotDateColumn = 1;
    }

    if (strcmp(extinction_header->header->column[columnIndex], "timeAccuracy") == 0) {
      gotTimeColumn = 1;
    }
  }

  if (gotDateColumn && gotTimeColumn) {
    extinction_table = table_loadva(
      extinction_handle,
      &extinction_header,
      NULL, /* hbase */
      NULL, /* rows */
      NULL,
      sizeof(EXTINCTION),
      &extinction_nrecs,
      TblInt, "enx", TblOff(PEXTINCTION, enx),
      TblInt, "eny", TblOff(PEXTINCTION, eny),
      TblInt, "eix", TblOff(PEXTINCTION, eix),
      TblInt, "eiy", TblOff(PEXTINCTION, eiy),
      TblDbl, "altitude", TblOff(PEXTINCTION, altitude),
      TblDbl, "extinction", TblOff(PEXTINCTION, extinction),
      TblDbl, "Date", TblOff(PEXTINCTION, Date),
      TblDbl, "timeAccuracy", TblOff(PEXTINCTION, timeAccuracy),
      0, "end", 0
    );
  } else {
    extinction_table = table_loadva(
      extinction_handle,
      &extinction_header,
      NULL, /* hbase */
      NULL, /* rows */
      NULL,
      sizeof(EXTINCTION),
      &extinction_nrecs,
      TblInt, "enx", TblOff(PEXTINCTION, enx),
      TblInt, "eny", TblOff(PEXTINCTION, eny),
      TblInt, "eix", TblOff(PEXTINCTION, eix),
      TblInt, "eiy", TblOff(PEXTINCTION, eiy),
      TblDbl, "altitude", TblOff(PEXTINCTION, altitude),
      TblDbl, "extinction", TblOff(PEXTINCTION, extinction),
      0, "end", 0
    );
  }

  if (extinction_table == NULL) {
    printf("ERROR: Failed to read table for %s\n", extinction_name);
    return 1;
  }

  CalcMemory(1, extinction_nrecs * sizeof(EXTINCTION), &curMemory, &maxMemory);
  if (verbose) {
    printf("read %zu records for %s\n", extinction_nrecs, extinction_name);
  }

  // Validate this table. If all of the extinction values are zero, then we know
  // that the plate has uncertain extinction and we must mark every image as
  // "low altitiude".

  enx = extinction_table[0].enx;
  eny = extinction_table[0].eny;

  if (extinction_nrecs !=  enx * eny) {
    printf("ERROR: Number of records %zu does not match product of %d and %d\n", extinction_nrecs, enx, eny);
    return 1;
  }

  irec = 0;

  for (eix = 0; eix < enx; eix++) {
    for (eiy = 0; eiy < eny; eiy++) {
      if (extinction_table[irec].altitude > maxAltitude) {
        maxAltitude = extinction_table[irec].altitude;
        minExtinction = extinction_table[irec].extinction;
      }

      if (
        extinction_table[irec].enx != enx &&
        extinction_table[irec].eny != eny &&
        extinction_table[irec].eix != eix &&
        extinction_table[irec].eiy != eiy
      ) {
        printf(
          "ERROR Mismatch of nx %d %d, ny %d %d, ix %d %d, or iy %d %d\n",
          extinction_table[irec].enx,
          enx,
          extinction_table[irec].eny,
          eny,
          extinction_table[irec].eix,
          eix,
          extinction_table[irec].eiy,
          eiy
        );
        return 1;
      }

      irec++;
    }
  }

  if (!(gotDateColumn && gotTimeColumn)) {
    for (irec = 0; irec < extinction_nrecs; irec++) {
      extinction_table[irec].Date = 0.0;
      extinction_table[irec].timeAccuracy = 99.0;
    }

    extinctionFileTimeAccuracy = 99.0;
    extinctionFileDate = 0.0;
  } else {
    extinctionFileTimeAccuracy = extinction_table[0].timeAccuracy;
    extinctionFileDate = extinction_table[0].Date;

    for (irec = 0; irec < extinction_nrecs; irec++) {
      if (
        extinction_table[irec].Date != extinctionFileDate ||
        extinction_table[irec].timeAccuracy != extinctionFileTimeAccuracy
      ) {
        if (printTimeAccuracy) {
          printf(
            "ERROR Mismatch of Julian Date %f %f or timeAccuracy %f %f in record %d of %s\n",
            extinction_table[irec].Date,
            extinctionFileDate,
            extinction_table[irec].timeAccuracy,
            extinctionFileTimeAccuracy,
            irec,
            extinction_name
          );
          printTimeAccuracy = 0;
        }

        return 1;
      }
    }
  }

  // Go to the database and get the exposure records for this mosaic, and check
  // them against the extinction file

  if (
    GetSolutionJulianDate(
      pConnection,
      series,
      plateNumber,
      mosaicNumber,
      solutionNumber,
      &exposureNumber,
      &numExposures,
      &singleJulianDate,
      &singleTolerance,
      &allJulianDate,
      &allTolerance
    )
  ) {
    printf("ERROR: failed to get exposure information for %s\n", fileroot);
    return 1;
  }

  // Update the plate quality bits if necessary
  UpdateQuality(pConnection, series, plateNumber);
  if (SetPlateQuality(pConnection, "multiple", fileroot, 0, 1, &quality) != 0) {
    printf("ERROR: failed to get plate quality information for %s\n", fileroot);
    return 1;
  }

  locationID = GetLocationId(pConnection, series, plateNumber);

  if (locationID <= 0) {
    printf("ERROR: failed to get a location ID for %s\n", fileroot);
    return 1;
  } else if (GetLongitude(pConnection, locationID, &longitude, &latitude, &elevation, NULL, 0) == 0) {
    printf("ERROR: failed to get a longitude for %s\n",fileroot);
    return 1;
  }

  plateColor = GetPlateColor(pConnection, series, plateNumber);
  extinctionCoefficient = GetExtinctionCoefficient(plateColor, elevation);

  // Done with the MySQL connection
  mysql_close(pConnection);

  if (!(gotDateColumn && gotTimeColumn)) {
    // Older version of prepare_octave.  We are o.k. if there is only one exposure record
    if (numExposures != 1) {
      singleJulianDate = allJulianDate;
      singleTolerance = allTolerance;
    }
  } else {
    dateDifference = singleJulianDate - extinctionFileDate;
    accuracyDifference = singleTolerance - extinctionFileTimeAccuracy;

    if (minExtinction > 0.0 && maxAltitude > 0.0) {
      z1 = sin(maxAltitude * DEGREES_TO_RAD);

      if (1 / z1 <= MAX_AIRMASS) {
        // Don't do this check if all of the stars are going to be flagged low altitiude
        extinctionDifference = minExtinction - extinctionCoefficient / z1;

        if (extinctionDifference < 0.0) {
          extinctionDifference = -extinctionDifference;
        }

        if (extinctionDifference > 0.001) {
          printf(
            "ERROR: stale extinction algorithm. extinction difference %f is too large in %s\n",
            extinctionDifference,
            extinction_name
          );
          exit(1);
        }
      }
    }

    // NOTE: the following uses a 1/10 second accuracy to match ProcLogbook.cpp, PrepareExposure routine
    if (
      dateDifference < -0.000001 ||
      dateDifference >  0.000001 ||
      accuracyDifference < -0.000001 ||
      accuracyDifference >  0.000001
    ) {
      printf(
        "ERROR: stale extinction file. date difference %f or time difference %f is too large in %s\n",
        dateDifference,
        accuracyDifference,
        extinction_name
      );
      singleJulianDate = allJulianDate;
      singleTolerance = allTolerance;
    }
  }

  geocentricJD = singleJulianDate;

  // Read in the sextractor file

  sextractor_header = table_header(sextractor_handle, TABLE_PARSE);
  if (sextractor_header == NULL) {
    printf("ERROR: Failed to read header for %s\n", sextractor_name);
    return 1;
  }

  sextractor_table = table_loadva(
    sextractor_handle,
    &sextractor_header,
    NULL, /* hbase */
    NULL, /* rows */
    NULL,
    sizeof(STARIMAGE),
    &sextractor_nrecs,
    TblInt, "NUMBER", TblOff(PSTARIMAGE, NUMBER),
    TblDbl, "X_IMAGE", TblOff(PSTARIMAGE, X_IMAGE),
    TblDbl, "Y_IMAGE", TblOff(PSTARIMAGE, Y_IMAGE),
    TblInt, "AFLAGS", TblOff(PSTARIMAGE, AFLAGS),
    TblInt, "BFLAGS", TblOff(PSTARIMAGE, BFLAGS),
    TblDbl, "MAG_ISO", TblOff(PSTARIMAGE, MAG_ISO),
    TblDbl, "ra", TblOff(PSTARIMAGE, ra),
    TblDbl, "dec", TblOff(PSTARIMAGE, dec),
    TblDbl, "FLUX_ISO", TblOff(PSTARIMAGE, FLUX_ISO),
    TblDbl, "MAG_APER", TblOff(PSTARIMAGE, MAG_APER),
    TblDbl, "MAG_AUTO", TblOff(PSTARIMAGE, MAG_AUTO),
    TblDbl, "KRON_RADIUS", TblOff(PSTARIMAGE, KRON_RADIUS),
    TblDbl, "BACKGROUND", TblOff(PSTARIMAGE, BACKGROUND),
    TblDbl, "THRESHOLD", TblOff(PSTARIMAGE, THRESHOLD),
    TblDbl, "FLUX_MAX", TblOff(PSTARIMAGE, FLUX_MAX),
    TblDbl, "THETA_J2000", TblOff(PSTARIMAGE, THETA_J2000),
    TblDbl, "ELLIPTICITY", TblOff(PSTARIMAGE, ELLIPTICITY),
    TblDbl, "ISOAREA_WORLD", TblOff(PSTARIMAGE, ISOAREA_WORLD),
    TblDbl, "FWHM_IMAGE", TblOff(PSTARIMAGE, FWHM_IMAGE),
    TblDbl, "FWHM_WORLD", TblOff(PSTARIMAGE, FWHM_WORLD),
    TblInt, "ISO0", TblOff(PSTARIMAGE, ISO0),
    TblInt, "ISO1", TblOff(PSTARIMAGE, ISO1),
    TblInt, "ISO2", TblOff(PSTARIMAGE, ISO2),
    TblInt, "ISO3", TblOff(PSTARIMAGE, ISO3),
    TblInt, "ISO4", TblOff(PSTARIMAGE, ISO4),
    TblInt, "ISO5", TblOff(PSTARIMAGE, ISO5),
    TblInt, "ISO6", TblOff(PSTARIMAGE, ISO6),
    TblInt, "ISO7", TblOff(PSTARIMAGE, ISO7),
    TblDbl, "plate_dist", TblOff(PSTARIMAGE, plate_dist),
    0, "end", 0
  );

  if (sextractor_table == NULL) {
    printf("ERROR: Failed to read table for %s\n", sextractor_name);
    return 1;
  }

  CalcMemory(1, sextractor_nrecs * sizeof(STARIMAGE), &curMemory, &maxMemory);

  if (verbose) {
    printf("read %zu records for %s\n", sextractor_nrecs, sextractor_name);
  }

  // Sort the table according to sextractor NUMBER */
  qsort((void *) sextractor_table, sextractor_nrecs, sizeof(STARIMAGE), ImageCompare);

  for (index = 0; index < sextractor_nrecs; index++) {
    sextractor_table[index].REF[0] = 0;
    sextractor_table[index].selected = 0;
    sextractor_table[index].npoints_local = 0;
    sextractor_table[index].spatial_bin = 0;
    ix = (sextractor_table[index].X_IMAGE * nx) / (1.0 * mosaicWidth);
    iy = (sextractor_table[index].Y_IMAGE * ny) / (1.0 * mosaicHeight);
    sextractor_table[index].local_bin = ix + nx * iy;
    sextractor_table[index].magcal_iso = 0.0;
    sextractor_table[index].magcal_iso_rms = 0.0;
    sextractor_table[index].magcal_local = 0.0;
    sextractor_table[index].magcal_local_error = 0.0;
    sextractor_table[index].magcal_magdep = 0.0;
    sextractor_table[index].magcal_magdep_rms = 99.0;
    sextractor_table[index].magdep_bin = -1;
    sextractor_table[index].magcal_local_rms = 0.0;
    sextractor_table[index].extinction = 0.0;
    sextractor_table[index].Stdmag = 0.0;
    sextractor_table[index].color = 0.0;
    sextractor_table[index].Blendedmag = 0.0;
    sextractor_table[index].outputFlag = 0;
    sextractor_table[index].heliocentricJD = 0;
    sextractor_table[index].limiting_mag_local = 0;
    sextractor_table[index].dradRMS2 = 0;
  }

  // Read in the "local" file

  local_header = table_header(local_handle, TABLE_PARSE);
  if (local_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",local_name);
    return 1;
  }

  local_table = table_loadva(
    local_handle,
    &local_header,
    NULL, /* hbase */
    NULL, /* rows */
    NULL,
    sizeof(STARIMAGE),
    &local_nrecs,
    TblInt, "NUMBER", TblOff(PSTARIMAGE, NUMBER),
    TblDbl, "X_IMAGE", TblOff(PSTARIMAGE, X_IMAGE),
    TblDbl, "Y_IMAGE", TblOff(PSTARIMAGE, Y_IMAGE),
    TblInt, "BFLAGS", TblOff(PSTARIMAGE, BFLAGS),
    TblInt, "cal_local", TblOff(PSTARIMAGE, selected),
    TblInt, "npoints_local", TblOff(PSTARIMAGE, npoints_local),
    TblInt, "spatial_bin", TblOff(PSTARIMAGE, spatial_bin),
    TblInt, "local_bin", TblOff(PSTARIMAGE, local_bin),
    TblDbl, "refmag", TblOff(PSTARIMAGE, Stdmag),
    TblDbl, "color", TblOff(PSTARIMAGE, color),
    TblDbl, "MAG_ISO", TblOff(PSTARIMAGE, MAG_ISO),
    TblDbl, "ra", TblOff(PSTARIMAGE, ra),
    TblDbl, "dec", TblOff(PSTARIMAGE, dec),
    TblDbl, "magcal_iso", TblOff(PSTARIMAGE, magcal_iso),
    TblDbl, "magcal_iso_rms", TblOff(PSTARIMAGE, magcal_iso_rms),
    TblDbl, "magcal_local", TblOff(PSTARIMAGE, magcal_local),
    TblDbl, "magcal_local_error", TblOff(PSTARIMAGE, magcal_local_error),
    TblDbl, "extinction", TblOff(PSTARIMAGE, extinction),
    TblBuf, "REF", TblOff(PSTARIMAGE, REF), MAX_REF,
    0, "end", 0
  );

  if (local_table == NULL) {
    printf("ERRORD: Failed to read table for %s\n", local_name);
    return 1;
  }

  CalcMemory(1, local_nrecs * sizeof(STARIMAGE), &curMemory, &maxMemory);
  if (verbose) {
    printf("read %zu records for %s\n", local_nrecs, local_name);
  }

  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
    pSextractor->AFLAGS = 0;
    pSextractor->Blendedmag = 0;
    pSextractor->magdep_bin = -1;
    pSextractor->magcal_magdep = 0;
    pSextractor->magcal_magdep_rms = 99.0;

    // "V3.6.5 - discard unreasonable colors"

    switch (catalogNumber) {
    case CATALOG_GSC232:
    case CATALOG_EXPERIMENTAL:
      if (pSextractor->color < MIN_GSC_COLOR || pSextractor->color > MAX_GSC_COLOR) {
        printf("ERROR: Stale GSC color value %f\n", pSextractor->color);
        pSextractor->color = 99.0;
        exit(1);
      }
      break;

    case CATALOG_KEPLER:
      if (pSextractor->color < MIN_KEPLER_COLOR || pSextractor->color > MAX_KEPLER_COLOR) {
        printf("ERROR: Stale KIC color value %f\n", pSextractor->color);
        pSextractor->color = 99.0;
        exit(1);
      }
      break;

    case CATALOG_APASS:
      if (pSextractor->color < MIN_APASS_COLOR || pSextractor->color > MAX_APASS_COLOR) {
        printf("ERROR: Stale APASS color value %f\n", pSextractor->color);
        pSextractor->color = 99.0;
        exit(1);
      }
      break;

    case CATALOG_GAIA:
      if (pSextractor->color < MIN_GAIA_COLOR || pSextractor->color > MAX_GAIA_COLOR) {
        printf("ERROR: Stale GAIA color value %f\n", pSextractor->color);
        pSextractor->color = 99.0;
        exit(1);
      }
      break;

    case CATALOG_ATLAS:
      if (pSextractor->color < MIN_ATLAS_COLOR || pSextractor->color > MAX_ATLAS_COLOR) {
        printf("ERROR: Stale ATLAS color value %f\n", pSextractor->color);
        pSextractor->color = 99.0;
        exit(1);
      }
      break;

    default:
      printf("ERROR: Unknown catalog in prepare_octave\n");
      exit(1);
      break;
    }
  }

  // Read in the magnitude-dependent correction table, if present, and apply the
  // magnitude-dependent correction to magcal_local.

  strcpy(magdep_name, ingestDirectory);
  strcat(magdep_name, "/");
  strcat(magdep_name, fileroot);
  strcat(magdep_name, qualifier);
  strcat(magdep_name, "_magdepcalibrate.db");

  if((magdep_load_status = LoadMagdepCorrections(magdep_name, pMagdepLimits, &pMagdepTable, verbose)) != 0) {
    // No magnitude-dependent calibration table exists. Set magcal_magdep to
    // magcal_local and magcal_magdep_rms to 99.
    for (index = 0; index < local_nrecs; index++) {
      pSextractor = &local_table[index];
      pSextractor->magcal_magdep = pSextractor->magcal_local;
      pSextractor->magcal_magdep_rms = 99.0;
    }
  } else {
    // Apply the magnitude-dependent correction .
    for (index = 0; index < local_nrecs; index++) {
      pSextractor = &local_table[index];

      pSextractor->magcal_magdep = pSextractor->magcal_local - GetMagdepBinMagcor(
        pMagdepLimits,
        pMagdepTable,
        pSextractor->X_IMAGE,
        pSextractor->Y_IMAGE,
        pSextractor->magcal_local,
        &pSextractor->magdep_bin,
        &pSextractor->magcal_magdep_rms,
        &temp_magdep_magcor,
        1
      );

      if (pSextractor->magcal_magdep_rms < 90) {
        pSextractor->BFLAGS |= (1 << FILTER_BFLAG_MAGDEP_MAGCOR);
      }
    }
  }

  // Initialize the local structures

  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    pRmsValueTable[spatial_bin] = (PRMSVALUE) calloc(RMSVALUE_INITIAL_ALLOC, sizeof(RMSVALUE));
    CalcMemory(1, RMSVALUE_INITIAL_ALLOC * sizeof(RMSVALUE), &curMemory, &maxMemory);
    if (pRmsValueTable[spatial_bin] == NULL) {
      printf("ERROR: Failed to allocate pRmsValueTable for spatial bin %d\n", spatial_bin);
      exit(1);
    }

    pRmsCalibTable[spatial_bin] = NULL;
    rmsValueAlloc[spatial_bin] = RMSVALUE_INITIAL_ALLOC;
    rmsValueCount[spatial_bin] = 0;
    binCalibCount[spatial_bin] = 0;
  }

  // Populate the difference value per spatial bin

  for (index = 0; index < local_nrecs; index++) {
    pSextractor = &local_table[index];
    spatial_bin = pSextractor->spatial_bin;
    if (spatial_bin < 1 || spatial_bin > MAX_SPATIAL_BINS) {
      printf("ERROR: Illegal spatial bin %d for NUMBER %d\n", pSextractor->spatial_bin, pSextractor->NUMBER);
      exit(1);
    }

    if (numBins == 1) {
      /* Use the bin 1 calibration for all stars in single bin mode */
      effective_spatial_bin = 1;
    } else {
      effective_spatial_bin = spatial_bin;
    }

    if (rmsValueCount[effective_spatial_bin] >= rmsValueAlloc[effective_spatial_bin]) {
      /* Need to get more spatial bin entries */
      rmsValueAlloc[effective_spatial_bin] += RMSVALUE_INITIAL_ALLOC;
      pRmsValue = realloc(pRmsValueTable[effective_spatial_bin], rmsValueAlloc[effective_spatial_bin] * sizeof(RMSVALUE));
      CalcMemory(1, RMSVALUE_INITIAL_ALLOC * sizeof(RMSVALUE), &curMemory, &maxMemory);
      if (pRmsValue == NULL) {
        printf("ERROR: failed to reallocate pRmsValueTable for spatial bin %d size %d\n", effective_spatial_bin, rmsValueAlloc[effective_spatial_bin]);
        exit(1);
      }

      pRmsValueTable[effective_spatial_bin] = pRmsValue;
    }

    pRmsValue = &(pRmsValueTable[effective_spatial_bin][rmsValueCount[effective_spatial_bin]]);
    pRmsValue->magcal_local = pSextractor->magcal_local;
    pRmsValue->diffval = pSextractor->magcal_local - pSextractor->Stdmag;
    rmsValueCount[effective_spatial_bin]++;
  }

  // Estimate the RMS as a function of magnitude for each spatial bin

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    double *vector;

    if (rmsValueCount[spatial_bin] <= 0) {
      continue;
    }

    qsort((void *)pRmsValueTable[spatial_bin], rmsValueCount[spatial_bin], sizeof(RMSVALUE), RmsValueCompare);

    binCalibCount[spatial_bin] = rmsValueCount[spatial_bin] / RMSCALIB_INCREMENT;
    if (binCalibCount[spatial_bin] == 0) {
      binCalibCount[spatial_bin]++;
    }

    pRmsCalib = (PRMSCALIB)calloc(binCalibCount[spatial_bin], sizeof(RMSCALIB));
    CalcMemory(1, binCalibCount[spatial_bin] * sizeof(RMSCALIB), &curMemory, &maxMemory);
    if (pRmsCalib == NULL) {
      printf("ERROR: failed to allocate pRmsCAlib for spatial bin %d, size %d\n", spatial_bin, binCalibCount[spatial_bin]);
      exit(1);
    }

    pRmsCalibTable[spatial_bin] = pRmsCalib;
    vector = (double *)calloc(rmsValueCount[spatial_bin], sizeof(double));
    CalcMemory(1,rmsValueCount[spatial_bin] * sizeof(double), &curMemory, &maxMemory);
    if (vector == NULL) {
      printf("ERROR: failed to allocate vector for spatial bin %d, size %d\n", spatial_bin, rmsValueCount[spatial_bin]);
      exit(1);
    }

    for (rmsCalibIndex = 0; rmsCalibIndex < binCalibCount[spatial_bin]; rmsCalibIndex++) {
      int minValueIndex = rmsCalibIndex * RMSCALIB_INCREMENT;
      int maxValueIndex = minValueIndex + RMSCALIB_INCREMENT;
      int valueIndex;
      double median;
      double rms;
      pRmsValue = NULL;

      if (rmsCalibIndex == binCalibCount[spatial_bin] - 1) {
        maxValueIndex = rmsValueCount[spatial_bin];
      }

      for (valueIndex = minValueIndex; valueIndex < maxValueIndex; valueIndex++) {
        pRmsValue = &(pRmsValueTable[spatial_bin][valueIndex]);
        vector[valueIndex - minValueIndex] = pRmsValue->diffval;
      }

      if (CalcMedianAndRMS(maxValueIndex - minValueIndex, 0, vector, &median, &rms, 1, 3.0, 1) > 0) {
        pRmsCalib = &(pRmsCalibTable[spatial_bin][rmsCalibIndex]);
        pRmsCalib->rms = rms;
        pRmsCalib->magcal_local = pRmsValue->magcal_local;
      } else {
        printf(
          "ERROR: Failed to calculate a rms curve for spatial bin %d number of stars %d for %s\n",
          spatial_bin,
          maxValueIndex - minValueIndex,
          fileroot
        );
        exit(1);
      }
    }

    free(vector);
    CalcMemory(-1, rmsValueCount[spatial_bin] * sizeof(double), &curMemory, &maxMemory);
    vector = NULL;
  }

  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (pRmsValueTable[spatial_bin] != NULL) {
      free(pRmsValueTable[spatial_bin]);
      CalcMemory(-1, rmsValueAlloc[spatial_bin] * sizeof(RMSVALUE), &curMemory, &maxMemory);
      rmsValueAlloc[spatial_bin] = 0;
      pRmsValueTable[spatial_bin] = NULL;
    }
  }

  // Sort the table according to sextractor NUMBER
  qsort((void *) local_table, local_nrecs, sizeof(STARIMAGE), ImageCompare);

  maxLength  = 0;
  min_X_IMAGE = 1.0 * mosaicWidth;
  max_X_IMAGE = 0.0;
  min_Y_IMAGE = 1.0 * mosaicHeight;
  max_Y_IMAGE = 0.0;

  for (index = 0; index < local_nrecs; index++) {
    if (local_table[index].X_IMAGE > max_X_IMAGE) {
      max_X_IMAGE = local_table[index].X_IMAGE;
    }

    if (local_table[index].X_IMAGE < min_X_IMAGE) {
      min_X_IMAGE = local_table[index].X_IMAGE;
    }

    if (local_table[index].Y_IMAGE > max_Y_IMAGE) {
      max_Y_IMAGE = local_table[index].Y_IMAGE;
    }

    if (local_table[index].Y_IMAGE < min_Y_IMAGE) {
      min_Y_IMAGE = local_table[index].Y_IMAGE;
    }

    curLength = strlen(local_table[index].REF);
    if (curLength > maxLength) {
      maxLength = curLength;
      maxLengthIndex = index;
    }

    if (local_table[index].selected != PIPELINE_FLAG) {
      printf("ERROR: cal_local is %d instead of %d for %s\n", local_table[index].selected, PIPELINE_FLAG, local_table[index].REF);
      return 1;
    }

    // Correct the local calibration RMS with the catalog RMS in quadrature
    // using "RMS algorithm E"

    if (local_table[index].magcal_local_error > 90.0) {
      // No local calibration here
      local_table[index].selected = PIPELINE_FLAG | LOWESS_CAL_FLAG;
      local_table[index].magcal_local_rms = 99.0;
      local_table[index].magcal_local_rms = local_table[index].magcal_iso_rms;
      local_table[index].magcal_local = local_table[index].magcal_iso;
    } else {
      local_table[index].selected = PIPELINE_FLAG | LOWESS_CAL_FLAG | LOCAL_CAL_FLAG;

      if (numBins == 1) {
        // Use the bin 1 calibration for all stars in single bin mode
        effective_spatial_bin = 1;
      } else {
        effective_spatial_bin = local_table[index].spatial_bin;
      }

      local_table[index].magcal_local_rms = CalcRMS_E(
        binCalibCount,
        pRmsCalibTable,
        effective_spatial_bin,
        local_table[index].magcal_local
      );
    }

    // If this object is in a bad spatial bin, flag it now because we need to
    // evaluate each local bin. We will flag everything else in the
    // sextractor_table later.

    // Always check the local bin and spatial bin values

    ix = (local_table[index].X_IMAGE * nx) / (1.0 * mosaicWidth);
    iy = (local_table[index].Y_IMAGE * ny) / (1.0 * mosaicHeight);

    local_bin = ix + nx * iy;
    if (local_bin != local_table[index].local_bin) {
      printf(
        "ERROR: local bin is %d instead of expected %d for %s\n",
        local_bin,
        local_table[index].local_bin,
        local_table[index].REF
      );
      return 1;
    }

    spatial_bin = CalculateBin(mosaicWidth, mosaicHeight, local_table[index].X_IMAGE, local_table[index].Y_IMAGE, &edgeDist);
    if (spatial_bin != local_table[index].spatial_bin) {
      printf(
        "ERROR: spatial bin is %d instead of expected %d dist %f for %s\n",
        spatial_bin,
        local_table[index].spatial_bin,
        edgeDist,
        local_table[index].REF
      );
      exit(1);
    }
  }

  if (verbose) {
    printf(
      "Bins %d Mosaic X range %.1f %.1f for  width %d  Y range %.1f %.1f for height %d\n",
      numBins,
      min_X_IMAGE,
      max_X_IMAGE,
      mosaicWidth,
      min_Y_IMAGE,
      max_Y_IMAGE,
      mosaicHeight
    );
  }

  if (maxLength >= MAX_REF) {
    printf(
      "ERROR: maxLength %d is greater than or equal to MAX_REF %d for %s\n",
      maxLength,
      MAX_REF,
      local_table[maxLengthIndex].REF
    );
    return 1;
  } else if (verbose) {
    printf("Maximum REF length is %d for %s\n", maxLength, local_table[maxLengthIndex].REF);
  }

  // Step through the sextractor table and replace the sextractor entries with
  // any locally calibrated entries.

  local_index = 0;
  pLocal = &local_table[local_index];

  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];

    if (pSextractor->NUMBER > maxSextractorNUMBER) {
      maxSextractorNUMBER = pSextractor->NUMBER;
    }

    while (pLocal->NUMBER < pSextractor->NUMBER) {
      local_index++;

      if (local_index >= local_nrecs) {
        break;
      }

      pLocal = &local_table[local_index];
    }

    if (pSextractor->NUMBER == pLocal->NUMBER) {
      // First copy over items not in pLocal
      pLocal->FLUX_ISO = pSextractor->FLUX_ISO;
      pLocal->MAG_APER = pSextractor->MAG_APER;
      pLocal->MAG_AUTO = pSextractor->MAG_AUTO;
      pLocal->KRON_RADIUS = pSextractor->KRON_RADIUS;
      pLocal->BACKGROUND = pSextractor->BACKGROUND;
      pLocal->THRESHOLD = pSextractor->THRESHOLD;
      pLocal->FLUX_MAX = pSextractor->FLUX_MAX;
      pLocal->THETA_J2000 = pSextractor->THETA_J2000;
      pLocal->ELLIPTICITY = pSextractor->ELLIPTICITY;
      pLocal->ISOAREA_WORLD = pSextractor->ISOAREA_WORLD;
      pLocal->FWHM_IMAGE = pSextractor->FWHM_IMAGE;
      pLocal->FWHM_WORLD = pSextractor->FWHM_WORLD;
      pLocal->ISO0 = pSextractor->ISO0;
      pLocal->ISO1 = pSextractor->ISO1;
      pLocal->ISO2 = pSextractor->ISO2;
      pLocal->ISO3 = pSextractor->ISO3;
      pLocal->ISO4 = pSextractor->ISO4;
      pLocal->ISO5 = pSextractor->ISO5;
      pLocal->ISO6 = pSextractor->ISO6;
      pLocal->ISO7 = pSextractor->ISO7;
      pLocal->plate_dist = pSextractor->plate_dist;

      // Always use the sextractor ra and dec
      pLocal->ra = pSextractor->ra;
      pLocal->dec = pSextractor->dec;

      memcpy(pSextractor, pLocal, sizeof(STARIMAGE));
    }

    // Calculate the spatial bin

    tmp_spatial_bin = CalculateBin(
      mosaicWidth,
      mosaicHeight,
      pSextractor->X_IMAGE,
      pSextractor->Y_IMAGE,
      &edgeDist
    );

    if (pSextractor->spatial_bin == 0) {
      pSextractor->spatial_bin = tmp_spatial_bin;
    } else if (tmp_spatial_bin != pSextractor->spatial_bin) {
      printf(
        "ERROR: spatial bin is %d %d for NUMBER %d for %s\n",
        tmp_spatial_bin,
        pSextractor->spatial_bin,
        pSextractor->NUMBER,
        fileroot
      );
      exit(1);
    }
  }

  // Done with the local table

  if (local_table != NULL) {
    Free(local_table);
    local_table = NULL;
  }

  CalcMemory(-1, local_nrecs * sizeof(STARIMAGE), &curMemory, &maxMemory);

  if (local_header != NULL) {
    table_hdrfree(local_header);
    local_header = NULL;
  }

  if (local_handle != NULL) {
    Close(local_handle);
    local_handle = NULL;
  }

  pLocal = NULL;

  // Read in the match file

  match_header = table_header(match_handle, TABLE_PARSE);
  if (match_header == NULL) {
    printf("ERROR: Failed to read header for %s\n", match_name);
    return 1;
  }

  match_table = table_loadva(
    match_handle,
    &match_header,
    NULL, /* hbase */
    NULL, /* rows */
    NULL,
    sizeof(STARIMAGE),
    &match_nrecs,
    TblInt, "NUMBER", TblOff(PSTARIMAGE, NUMBER),
    TblDbl, "X_IMAGE", TblOff(PSTARIMAGE, X_IMAGE),
    TblDbl, "Y_IMAGE", TblOff(PSTARIMAGE, Y_IMAGE),
    TblInt, "AFLAGS", TblOff(PSTARIMAGE, AFLAGS),
    TblInt, "BFLAGS", TblOff(PSTARIMAGE, BFLAGS),
    TblDbl, "MAG_ISO", TblOff(PSTARIMAGE, MAG_ISO),
    TblDbl, "ra_2", TblOff(PSTARIMAGE, ra),
    TblDbl, "dec_2", TblOff(PSTARIMAGE, dec),
    TblDbl, "Stdmag", TblOff(PSTARIMAGE, Stdmag),
    TblDbl, "color", TblOff(PSTARIMAGE, color),
    TblDbl, "dra", TblOff(PSTARIMAGE, dra),
    TblDbl, "ddec", TblOff(PSTARIMAGE, ddec),
    TblBuf, "REF", TblOff(PSTARIMAGE, REF), MAX_REF,
    TblInt, "gsc_bin_index", TblOff(PSTARIMAGE, gsc_bin_index),
    0, "end", 0
  );

  if (match_table == NULL) {
    printf("ERROR: Failed to read table for %s\n", match_name);
    return 1;
  }

  CalcMemory(1, match_nrecs * sizeof(STARIMAGE), &curMemory, &maxMemory);
  if (verbose) {
    printf("read %zu records for %s\n", match_nrecs, match_name);
  }

  // Sort the table according to sextractor NUMBER
  qsort((void *) match_table, match_nrecs, sizeof(STARIMAGE), ImageCompare);

  for (index = 0; index < match_nrecs; index++) {
    curLength = strlen(match_table[index].REF);

    if (curLength > maxLength) {
      maxLength = curLength;
      maxLengthIndex = index;
    }
  }

  if (maxLength >= MAX_REF) {
    printf(
      "ERROR: maxLength %d is greater than or equal to MAX_REF %d for %s\n",
      maxLength,
      MAX_REF,
      match_table[maxLengthIndex].REF
    );
    return 1;
  } else if (verbose) {
    printf("Maximum REF length is %d for %s\n", maxLength, match_table[maxLengthIndex].REF);
  }

  for (index = 0; index < match_nrecs; index++) {
    match_table[index].selected = 0;
    match_table[index].npoints_local = 0;
    match_table[index].spatial_bin = 0;
    match_table[index].local_bin = 0;
    match_table[index].magcal_iso = 0.0;
    match_table[index].magcal_iso_rms = 0.0;
    match_table[index].magcal_local = 0.0;
    match_table[index].magcal_local_error = 0.0;
    match_table[index].magcal_magdep = 0.0;
    match_table[index].magcal_magdep_rms = 99.0;
    match_table[index].magdep_bin = -1;
    match_table[index].magcal_local_rms = 0.0;
    match_table[index].extinction = 0.0;

    // "V3.6.5 - discard unreasonable colors"

    switch (catalogNumber) {
    case CATALOG_GSC232:
    case CATALOG_EXPERIMENTAL:
      if (match_table[index].color < MIN_GSC_COLOR || match_table[index].color > MAX_GSC_COLOR) {
        match_table[index].color = 99.0;
      }
      break;

    case CATALOG_KEPLER:
      if (match_table[index].color < MIN_KEPLER_COLOR || match_table[index].color > MAX_KEPLER_COLOR) {
        match_table[index].color = 99.0;
      }
      break;

    case CATALOG_APASS:
      if (match_table[index].color < MIN_APASS_COLOR || match_table[index].color > MAX_APASS_COLOR) {
        match_table[index].color = 99.0;
      }
      break;

    case CATALOG_GAIA:
      if (match_table[index].color < MIN_GAIA_COLOR || match_table[index].color > MAX_GAIA_COLOR) {
        match_table[index].color = 99.0;
      }
      break;

    case CATALOG_ATLAS:
      if (match_table[index].color < MIN_ATLAS_COLOR || match_table[index].color > MAX_ATLAS_COLOR) {
        match_table[index].color = 99.0;
      }
      break;

    default:
      printf("ERROR: Unknown catalog in prepare_octave\n");
      exit(1);
      break;
    }
  }

  // Step through the sextractor table and replace any missing GSC2.3.2 ID's from the match table

  match_index = 0;
  pMatch = &match_table[match_index];

  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];
    pSextractor->gsc_bin_index = -1;

    while (pMatch->NUMBER < pSextractor->NUMBER) {
      match_index++;

      if (match_index >= match_nrecs) {
        break;
      }

      pMatch = &match_table[match_index];
    }

    if (pSextractor->NUMBER == pMatch->NUMBER) {
      if (strlen(pSextractor->REF) != 0 && strcmp(pSextractor->REF, pMatch->REF) != 0) {
        printf(
          "ERROR: NUMBER %d has a stale ref %s which should be %s for %s\n",
          pSextractor->NUMBER,
          pSextractor->REF,
          pMatch->REF,
          fileroot
        );
        exit(1);
      }

      pSextractor->dra = pMatch->dra;
      pSextractor->ddec = pMatch->ddec;
      pSextractor->gsc_bin_index = pMatch->gsc_bin_index; // Move the object to the "J2000" bin of the gsc coordinates

      if (pSextractor->Stdmag != 0) {
        stdmagdiff = pSextractor->Stdmag - pSextractor->extinction - pMatch->Stdmag;

        if (numBins == 1) {
          /* Use the bin 1 calibration for all stars in single bin mode */
          effective_spatial_bin = 1;
        } else {
          effective_spatial_bin = pSextractor->spatial_bin;
        }

        if (spatial_bin_table[effective_spatial_bin].colorflag > 0 && pSextractor->color < 99.0) {
          stdmagdiff -= pMatch->color * spatial_bin_table[effective_spatial_bin].colorterm;
        }

        if ((stdmagdiff < -0.02 || stdmagdiff > +0.02) && printStdmagMatchError) {
          printf(
            "ERROR: Stdmag %f %f NUMBER %d mismatch for %s \n",
            pSextractor->Stdmag,
            pMatch->Stdmag,
            pSextractor->NUMBER,
            fileroot
          );
          printStdmagMatchError = 0;
        }
      }

      // It is very important to replace all Stdmag values from the locally
      // corrected magnitude table with the values from the match table, because
      // the former have been adjusted for extinction and color.

      pSextractor->Stdmag = pMatch->Stdmag;
      pSextractor->color  = pMatch->color;

      if (pSextractor->selected == 0) {
        strcpy(pSextractor->REF, pMatch->REF);
        pSextractor->AFLAGS = pMatch->AFLAGS;
        pSextractor->BFLAGS = pMatch->BFLAGS;
      } else if ((pMatch->BFLAGS & (1 << FILTER_BFLAG_KEPLER)) != 0) {
          // We have cleared the kepler source bits and must now put them back
          tmpMAGFlag = pMatch->BFLAGS >> GSC_MAGNITUDE_FLAG_BIT;
          tmpMAGFlag &= MAGNITUDE_MASK;
          pSextractor->BFLAGS |= (tmpMAGFlag << GSC_MAGNITUDE_FLAG_BIT);
      }
    }

    if (pSextractor->gsc_bin_index < 0) {
      // Find the gsc catalog bin for all unmatched stars
      int decBin;
      int raBin;
      pSextractor->gsc_bin_index = GetGSCBin(pGscBin, pSextractor->ra, pSextractor->dec, &decBin, &raBin, fileroot);
    }
  }

  // Done with the match table.

  if (match_table != NULL) {
    Free(match_table);
    match_table = NULL;
  }

  CalcMemory(-1, match_nrecs * sizeof(STARIMAGE), &curMemory, &maxMemory);

  if (match_header != NULL) {
    table_hdrfree(match_header);
    match_header = NULL;
  }

  if (match_handle != NULL) {
    Close(match_handle);
    match_handle = NULL;
  }

  pMatch = NULL;

  // Repeat the calculation of the local calibration bin drad for using all of
  // the stars that are still considered good stars

  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];

    ix = (pSextractor->X_IMAGE * nx) / (1.0 * mosaicWidth);
    iy = (pSextractor->Y_IMAGE * ny) / (1.0 * mosaicHeight);
    if (ix < 0) {
      ix = 0;
    }

    if (ix >= nx) {
      ix = nx - 1;
    }

    if (iy < 0) {
      iy = 0;
    }

    if (iy >= ny) {
      iy = ny - 1;
    }

    // These records are reversed from the local bin
    pSextractor->irec = iy + ny * ix;

    if (
      pSextractor->AFLAGS < FILTER_AFLAG_BAD &&
      (pSextractor->BFLAGS & FILTER_BMASK_BLEND) == 0 &&
      (pSextractor->BFLAGS & (1 << FILTER_BFLAG_PMERROR)) == 0 &&
      pSextractor->selected == (PIPELINE_FLAG | LOWESS_CAL_FLAG | LOCAL_CAL_FLAG)
    ) {
      // Add this one to the local calibration table
      pDmagcor = &dmagcor_table[pSextractor->irec];
      pSextractor->flink = pDmagcor->flink;
      pDmagcor->flink = pSextractor;
    }
  }

  // Tighten drad tolerance by rejecting bad spatial bins
  ComputeDrad2(nx, ny, dmagcor_table, dradlimit, &dradCount1, &curMemory, &maxMemory);

  // At this point, we need to find out which local calibration bins fail the
  // criterion: `median GSC2.3.2 magnitude > median (limiting_lowess magnitude)
  // + 0.5`

  {
    double *vector1;
    double *vector2;
    double medStdmag;
    double rmsStdmag;
    double medLimMag;
    double rmsLimMag;

    vector1 = (double *) calloc(local_nrecs, sizeof(double));
    vector2 = (double *) calloc(local_nrecs, sizeof(double));
    CalcMemory(1, local_nrecs * 2 * sizeof(double), &curMemory, &maxMemory);

    if (vector1 == NULL || vector2 == NULL) {
      printf("ERROR: failed to allocate sort vectors of size %zu for %s\n", local_nrecs, fileroot);
      return 1;
    }

    // For each star of interest, populate the queue beginning with the DMAGCOR table

    for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
      pSextractor = &sextractor_table[sextractor_index];

      if (pSextractor->REF[0] != 0) {
        double newdrad;

        pDmagcor = &dmagcor_table[pSextractor->irec];

        // Correct for draMedian2 and ddecMedian2
        if (pDmagcor->drad_bin_count2 > 0) {
          double factor;

          pSextractor->dec -= pDmagcor->ddecMedian2 / 3600.;
          pSextractor->ddec -= pDmagcor->ddecMedian2;

          factor = cos(DEGREES_TO_RAD * pSextractor->dec);
          if (factor != 0.0) {
            pSextractor->ra -= (pDmagcor->draMedian2 / 3600.) / factor;
            pSextractor->dra -= pDmagcor->draMedian2;
          }

          pSextractor->BFLAGS |= (1 << FILTER_BFLAG_DRAD_ADJUST);
          newdrad = sqrt(sqr(pSextractor->dra) + sqr(pSextractor->ddec));

          if (newdrad > pDmagcor->dradRMS2 * DRAD_REJECT_FACTOR || newdrad > dradlimit) {
            if (pSextractor->AFLAGS < FILTER_AFLAG_BAD) {
              pDmagcor->drad_reject_count2++;
              dradCount4++;
            }

            pSextractor->AFLAGS |= (1 << FILTER_AFLAG_DRAD);
          }
        }

        if (numBins == 1) {
          /* Use the bin 1 calibration for all stars in single bin mode */
          effective_spatial_bin = 1;
        } else {
          if (pSextractor->spatial_bin == 0) {
            pSextractor->spatial_bin = CalculateBin(
              mosaicWidth,
              mosaicHeight,
              pSextractor->X_IMAGE,
              pSextractor->Y_IMAGE,
              &edgeDist
            );
          }

          effective_spatial_bin = pSextractor->spatial_bin;
        }

        pSpatial_bin = &spatial_bin_table[effective_spatial_bin];

        if (pSextractor->selected == (PIPELINE_FLAG | LOWESS_CAL_FLAG | LOCAL_CAL_FLAG)) {
          pSextractor->flink = pDmagcor->flink;
          pDmagcor->flink = pSextractor;
        }
      }
    }

    for (local_bin = 0; local_bin < nx * ny; local_bin++) {
      int vectorCount = 0;

      pDmagcor = &dmagcor_table[local_bin];
      pSextractor = pDmagcor->flink;

      while (pSextractor != NULL) {
        spatial_bin = CalculateBin(
          mosaicWidth,
          mosaicHeight,
          pSextractor->X_IMAGE,
          pSextractor->Y_IMAGE,
          &edgeDist
        );

        vector1[vectorCount] = pSextractor->Stdmag;
        pSpatial_bin = &spatial_bin_table[spatial_bin];
        vector2[vectorCount] = pSpatial_bin->limiting_mag;
        vectorCount++;

        if (vectorCount >= local_nrecs) {
          printf("ERROR: zout vectors exceed %zu for %s\n", local_nrecs, fileroot);
          return 1;
        }

        pSextractor = pSextractor->flink;
      }

      if (vectorCount > 0) {
        if (vectorCount > maxVectorCount) {
          maxVectorCount = vectorCount;
        }

        CalcMedianAndRMS(vectorCount, 0, vector1, &medStdmag, &rmsStdmag, 0, 0, 0);
        CalcMedianAndRMS(vectorCount, 0, vector2, &medLimMag, &rmsLimMag, 0, 0, 0);

        if (medStdmag > medLimMag + 0.5) {
          // This is a bad bin.
          pDmagcor->rejectFlag |= DMAG_REJECT_MEDIAN;
          dmagMedianCount++;

          if ((pDmagcor->rejectFlag & DMAG_REJECT_HIZOUT) != 0) {
            dmagBothCount++;
          }
        }
      }
    }

    free(vector1);
    free(vector2);
    CalcMemory(-1, local_nrecs * 2 * sizeof(double), &curMemory, &maxMemory);
  }

  fprintf(
    stdout,
    "Number of local bins rejected for high zout %d, dim Stdmag %d, both %d drad %d, maxVectorCount %d for %s\n",
    dmagHiZoutCount,
    dmagMedianCount,
    dmagBothCount,
    dmagDradCount,
    maxVectorCount,
    fileroot
  );

  // Read in the blend file, if present

  if (blend_handle != NULL) {
    blend_header = table_header(blend_handle, TABLE_PARSE);
    if (blend_header == NULL) {
      printf("ERROR: Failed to read header for %s\n", blend_name);
      return 1;
    }

    blend_table = table_loadva(
      blend_handle,
      &blend_header,
      NULL, /* hbase */
      NULL, /* rows */
      NULL,
      sizeof(BLEND),
      &blend_nrecs,
      TblInt, "NUMBER", TblOff(PBLEND, NUMBER),
      TblDbl, "Stdmag", TblOff(PBLEND, Stdmag),
      TblDbl, "Blendedmag", TblOff(PBLEND, Blendedmag),
      0, "end", 0
    );

    if (blend_table == NULL) {
      /* assume this is an empty table */
      printf("WARNING: Failed to read table for %s\n", blend_name);
      blend_nrecs = 0;
    }
  }

  CalcMemory(1, blend_nrecs * sizeof(BLEND), &curMemory, &maxMemory);

  if (verbose) {
    printf("read %zu records for %s\n", blend_nrecs, blend_name);
  }

  if (blend_nrecs > 0) {
    // Sort the table according to sextractor NUMBER
    qsort((void *) blend_table, blend_nrecs, sizeof(BLEND), ImageCompare);

    // Go through the image table and adjust the magnitudes of any blended
    // stars.
    blend_index = 0;
    pBlend = &blend_table[blend_index];

    for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
      pSextractor = &sextractor_table[sextractor_index];

      while (pBlend->NUMBER < pSextractor->NUMBER) {
        blend_index++;

        if (blend_index >= blend_nrecs) {
          break;
        }

        pBlend = &blend_table[blend_index];
      }

      if (pSextractor->NUMBER == pBlend->NUMBER) {
        if (pSextractor->Stdmag - pBlend->Stdmag > 0.01 || pSextractor->Stdmag - pBlend->Stdmag < -0.01) {
          // Skip stars that have been purged because of an attempt to write out
          // the REFnumber twice.

          if (pSextractor->Stdmag != 0.0 || pSextractor->selected != 0 || pSextractor->REF[0] != 0) {
            printf(
              "ERROR: Stdmag %f in input file does not agree with blend Stdmag %f for %s\n",
              pSextractor->Stdmag,
              pBlend->Stdmag,
              fileroot
            );
            exit(1);
          }
        } else {
          pSextractor->Blendedmag = pBlend->Blendedmag;
        }
      }
    }
  }

  // Done with the blend table.

  if (blend_table != NULL) {
    Free(blend_table);
    blend_table = NULL;
  }

  CalcMemory(-1, blend_nrecs * sizeof(BLEND), &curMemory, &maxMemory);

  if (blend_header != NULL) {
    table_hdrfree(blend_header);
    blend_header = NULL;
  }

  if (blend_handle != NULL) {
    Close(blend_handle);
    blend_handle = NULL;
  }

  pBlend = NULL;

  // Read the defect file, if available

  if (defect_handle != NULL) {
    defect_header = table_header(defect_handle, TABLE_PARSE);
    if (defect_header == NULL) {
      printf("ERROR: Failed to read header for %s\n", defect_name);
      return 1;
    }

    defect_table = table_loadva(
      defect_handle,
      &defect_header,
      NULL, /* hbase */
      NULL, /* rows */
      NULL,
      sizeof(DEFECT),
      &defect_nrecs,
      TblInt, "NUMBER", TblOff(PDEFECT, NUMBER),
      TblInt, "AFLAGS", TblOff(PDEFECT, AFLAGS),
      0, "end", 0
    );

    if (defect_table == NULL) {
      printf("ERROR: Failed to read table for %s\n", defect_name);
      return 1;
    }

    CalcMemory(1, defect_nrecs * sizeof(DEFECT), &curMemory, &maxMemory);

    if (verbose) {
      printf("read %zu records for %s\n", defect_nrecs, defect_name);
    }

    // Verify that the defect table is sorted numerically by NUMBER
    for (defect_index = 1; defect_index < defect_nrecs; defect_index++) {
      pDefect = &defect_table[defect_index];
      pDefectPrev = &defect_table[defect_index - 1];

      if (pDefect->NUMBER <= pDefectPrev->NUMBER) {
        printf(
          "ERROR: defect table is not sorted by NUMBER at %d %d %d for %s\n",
          defect_index,
          pDefectPrev->NUMBER,
          pDefect->NUMBER,
          fileroot
        );
        return 1;
      }
    }

    // Step through the sextractor table and set the FILTER_AFLAG_DEFECT bit for
    // any plate defects

    defect_index = 0;
    pDefect = &defect_table[defect_index];

    for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
      pSextractor = &sextractor_table[sextractor_index];

      while (pDefect->NUMBER < pSextractor->NUMBER) {
        defect_index++;

        if (defect_index >= defect_nrecs) {
          break;
        }

        pDefect = &defect_table[defect_index];
      }

      if (pDefect->NUMBER == pSextractor->NUMBER) {
        if ((pDefect->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) != 0) {
          pSextractor->AFLAGS |= (1 << FILTER_AFLAG_DEFECT);
        } else if ((pSextractor->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) != 0) {
          pSextractor->AFLAGS &= ~(1 << FILTER_AFLAG_DEFECT) ;
          printf(
            "ERROR: stale defect flag found for NUMBER %d for %s\n",
            pSextractor->NUMBER,
            fileroot
          );
        }
      }
    }
  }

  // Done with the defect table.

  if (defect_table != NULL) {
    free(defect_table);
    defect_table = NULL;
  }

  CalcMemory(-1, defect_nrecs * sizeof(DEFECT), &curMemory, &maxMemory);

  if (defect_header != NULL) {
    table_hdrfree(defect_header);
    defect_header = NULL;
  }

  if (defect_handle != NULL) {
    Close(defect_handle);
    defect_handle = NULL;
  }

  pDefect = NULL;
  pDefectPrev = NULL;

  // Read in the background file, if available

  if (background_handle != NULL) {
    background_header = table_header(background_handle, TABLE_PARSE);
    if (background_header == NULL) {
      printf("ERROR: Background header read error for %s\n", background_name);
      return 1;
    }

    background_table = table_loadva(
      background_handle,
      &background_header,
      NULL, /* hbase */
      NULL, /* rows */
      NULL,
      sizeof(HIGHBACKGROUND),
      &background_nrecs,
      TblInt, "NUMBER", TblOff(PHIGHBACKGROUND, NUMBER),
      TblInt, "AFLAGS", TblOff(PHIGHBACKGROUND, AFLAGS),
      0, "end", 0
    );

    if (background_table == NULL) {
      // This is o.k. if the plate has a good background
      printf("WARNING: Failed to read table for %s\n", background_name);
    } else {
      CalcMemory(1, background_nrecs * sizeof(HIGHBACKGROUND), &curMemory, &maxMemory);

      if (verbose) {
        printf("read %zu records for %s\n", background_nrecs, background_name);
      }

      // Verify that the background table is sorted numerically by NUMBER
      for (background_index = 1; background_index < background_nrecs; background_index++) {
        pBackground = &background_table[background_index];
        pBackgroundPrev = &background_table[background_index - 1];

        if (pBackground->NUMBER <= pBackgroundPrev->NUMBER) {
          printf(
            "ERROR: background table is not sorted by NUMBER at %d %d %d for %s\n",
            background_index,
            pBackgroundPrev->NUMBER,
            pBackground->NUMBER,
            fileroot
          );
          return 1;
        }
      }

      // Step through the sextractor table and set the FILTER_AFLAG_BACKGROUND
      // bit for any objects with high background

      background_index = 0;
      pBackground = &background_table[background_index];

      for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
        pSextractor = &sextractor_table[sextractor_index];

        while (pBackground->NUMBER < pSextractor->NUMBER) {
          background_index++;

          if (background_index >= background_nrecs) {
            break;
          }

          pBackground = &background_table[background_index];
        }

        if (pBackground->NUMBER == pSextractor->NUMBER) {
          if ((pBackground->AFLAGS & (1 << FILTER_AFLAG_BACKGROUND)) != 0) {
            pSextractor->AFLAGS |= (1 << FILTER_AFLAG_BACKGROUND);
          } else if ((pSextractor->AFLAGS & (1 << FILTER_AFLAG_BACKGROUND)) != 0) {
            pSextractor->AFLAGS &= ~(1 << FILTER_AFLAG_BACKGROUND);
            printf(
              "ERROR: stale background flag found for NUMBER %d for %s\n",
              pSextractor->NUMBER,
              fileroot
            );
          }
        }
      }
    }
  }

  // Done with the background table.

  if (background_table != NULL) {
    free(background_table);
    background_table = NULL;
  }

  CalcMemory(-1, background_nrecs * sizeof(HIGHBACKGROUND), &curMemory, &maxMemory);

  if (background_header != NULL) {
    table_hdrfree(background_header);
    background_header = NULL;
  }

  if (background_handle != NULL) {
    Close(background_handle);
    background_handle = NULL;
  }

  pBackground = NULL;
  pBackgroundPrev = NULL;

  // We are finally ready to process all of the sextractor entries that have not
  // been assigned magnitudes

  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];

    tooBrightFlag = 0;

    // Get the magnitude correction record for this bin. These records are
    // reversed from the local bin.

    pDmagcor = &dmagcor_table[pSextractor->irec];
    pSextractor->magcal_local_error = pDmagcor->errout;
    pSextractor->npoints_local = pDmagcor->npout;

    if (pDmagcor->zout < MINIMUM_ZOUT || pDmagcor->zout > MAXIMUM_ZOUT || pDmagcor->errout > MAXIMUM_ERROUT) {
      pDmagcor->rejectFlag |= DMAG_REJECT_HIZOUT;
    }

    if ((pDmagcor->rejectFlag & (DMAG_REJECT_HIZOUT | DMAG_REJECT_MEDIAN)) != 0) {
      pSextractor->AFLAGS |= (1 << FILTER_AFLAG_HIZOUT);
      highZoutCount++;
    } else if ((pDmagcor->rejectFlag & DMAG_REJECT_DRAD) != 0) {
      if (pSextractor->AFLAGS < FILTER_AFLAG_BAD) {
        dradCount2++;
      }

      pSextractor->AFLAGS |= (1 << FILTER_AFLAG_DRADBIN);
    }

    // Get the extinction record for this bin

    eix = (pSextractor->X_IMAGE * enx) / (1.0 * mosaicWidth);
    eiy = (pSextractor->Y_IMAGE * eny) / (1.0 * mosaicHeight);

    if (eix < 0) {
      eix = 0;
    }

    if (eix >= enx) {
      eix = enx - 1;
    }

    if (eiy < 0) {
      eiy = 0;
    }

    if (eiy >= eny) {
      eiy = eny - 1;
    }

    // These records are reversed from the local bin
    irec = eiy + eny * eix;
    pExtinction = &extinction_table[irec];

    // First apply lowess calibration

    spatial_bin = CalculateBin(
      mosaicWidth,
      mosaicHeight,
      pSextractor->X_IMAGE,
      pSextractor->Y_IMAGE,
      &edgeDist
    );

    if (spatial_bin < 1 || spatial_bin > MAX_SPATIAL_BINS) {
      illegalBinCount++;
      continue;
    }

    if (numBins == 1) {
      /* Use the bin 1 calibration for all stars in single bin mode */
      effective_spatial_bin = 1;
    } else {
      effective_spatial_bin = spatial_bin;
    }

    pSpatial_bin = &spatial_bin_table[effective_spatial_bin];
    if (pSpatial_bin->spatial_bin != effective_spatial_bin) {
      missingBinCount++; // No lowess fit exists for this bin
      continue;
    }

    if (pSextractor->selected == 0) {
      limiting_iso = pSpatial_bin->limiting_iso;
      max_bright_iso = pSpatial_bin->max_bright_iso;

      if (pSextractor->MAG_ISO > limiting_iso) {
        tooDimCount++; /* This star is too dim for processing */
        continue;
      }

      if (pSextractor->MAG_ISO < max_bright_iso + MAX_BRIGHT_ADJUST) {
        tooBrightCount++; /* This star is too bright for processing */
        tooBrightFlag = 1;

        pSextractor->magcal_iso = pSpatial_bin->max_bright_mag;
        pSextractor->magcal_iso_rms = 99.0;
        pSextractor->magcal_local = pSpatial_bin->max_bright_mag;
        pSextractor->magcal_local_rms = 99.0;
        pSextractor->magcal_magdep = pSextractor->magcal_local;
        pSextractor->magcal_magdep_rms = 99.0;
        pSextractor->selected |= TOO_BRIGHT_FLAG;
        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_TOO_BRIGHT);
      }

      if (tooBrightFlag == 0) {
        if (
          CalculateLowessMagnitude(
            pSextractor->MAG_ISO,
            &magcal_iso,
            &magcal_iso_rms,
            maggrid_table[effective_spatial_bin],
            maggrid_nrecs[effective_spatial_bin]
          )
        ) {
          if (magcal_iso_rms < 0.0) {
            printf(
              "ERROR: line %d negative magcal_iso_rms %f for NUMBER %d for %s\n",
              __LINE__,
              magcal_iso_rms,
              pSextractor->NUMBER,
              fileroot
            );
            magcal_iso_rms = 99.0;
          }

          pSextractor->selected |= LOWESS_CAL_FLAG;
          pSextractor->magcal_iso = magcal_iso;
          pSextractor->magcal_iso_rms = magcal_iso_rms;
          pSextractor->magcal_local_rms = magcal_iso_rms;
        } else {
          printf(
            "ERRORD: CalculateLowessMagnitude failed for MAG_ISO %f spatial_bin %d star %s plate %s\n",
            pSextractor->MAG_ISO,
            effective_spatial_bin,
            pSextractor->REF,
            fileroot
          );
        }

        if (pDmagcor->errout > 90.0 || pSextractor->selected == 0) {
          // No local calibration here
          pSextractor->magcal_local_rms = 99.0;
          pSextractor->magcal_local_rms = pSextractor->magcal_iso_rms;
          pSextractor->magcal_local_error = pDmagcor->errout;
          pSextractor->npoints_local = pDmagcor->npout;
          magcal_local = pSextractor->magcal_iso;

          if (pDmagcor->errout > 90.0) {
            noLocalPoints++;
          }
        } else {
          magcal_local = pSextractor->magcal_iso - pDmagcor->zout;
          pSextractor->magcal_local_error = pDmagcor->errout;
          pSextractor->npoints_local = pDmagcor->npout;
          pSextractor->selected |= LOCAL_CAL_FLAG;
          pSextractor->magcal_local_rms = CalcRMS_E(binCalibCount, pRmsCalibTable, effective_spatial_bin, magcal_local);
        }

        pSextractor->magcal_local = magcal_local;

        // Apply the magnitude-dependent local calibration
        if ((pSextractor->BFLAGS & (1 << FILTER_BFLAG_MAGDEP_MAGCOR)) != 0) {
          // We shouln't be here for already selected objects
          printf(
            "ERROR: double magnitude-dependent correction for NUMBER %d of %s\n",
            pSextractor->NUMBER,
            fileroot
          );
          exit(1);
        }

        if (magdep_load_status == 0)  {
          pSextractor->magcal_magdep = pSextractor->magcal_local - GetMagdepBinMagcor(
            pMagdepLimits, pMagdepTable,
            pSextractor->X_IMAGE,
            pSextractor->Y_IMAGE,
            pSextractor->magcal_local,
            &pSextractor->magdep_bin,
            &pSextractor->magcal_magdep_rms,
            &temp_magdep_magcor,
            1
          );

          if (pSextractor->magcal_magdep_rms < 90) {
            pSextractor->BFLAGS |= (1 << FILTER_BFLAG_MAGDEP_MAGCOR);
          }
        } else {
          pSextractor->magcal_magdep = pSextractor->magcal_local;
          pSextractor->magcal_magdep_rms = 99.0;
        }
      }
    }

    // In single bin mode, do not write out any uncalibrated objects beyond the
    // limits of the input data set

    if (
      numBins == 1 &&
      pSextractor->selected != 0 &&
      (pSextractor->selected & PIPELINE_FLAG) == 0 &&
      (
        pSextractor->X_IMAGE < min_X_IMAGE ||
        pSextractor->X_IMAGE > max_X_IMAGE ||
        pSextractor->Y_IMAGE < min_Y_IMAGE ||
        pSextractor->Y_IMAGE > max_Y_IMAGE
      )
    ) {
      pSextractor->selected = 0;
      singleBinCount++;
    }

    // Write out this entry if it is of interest
    if (pSextractor->selected != 0) {
      if (pSextractor->REF[0] == 0) {
        noRefCount++;
        strcpy(pSextractor->REF, "NONE");
        pSextractor->Stdmag = 99.0;
        pSextractor->color = 0;
        pSextractor->dra = 0.0;
        pSextractor->ddec = 0.0;
      }

      heliocentricJD = jd2hjd(geocentricJD, pSextractor->ra, pSextractor->dec, WCS_J2000);
      limiting_mag_local = pSpatial_bin->limiting_mag;

      if (pDmagcor->errout < 90.0) {
        // We have an extinction correction in this bin. Apply it also to the
        // limiting magnitude
        limiting_mag_local -= pDmagcor->zout;
      }

      if (magdep_load_status == 0)  {
        // The magnitude-dependent correction for the limiting magnitude is limiting_mag_iso
        limiting_mag_local -= GetMagdepBinMagcor(
          pMagdepLimits,
          pMagdepTable,
          pSextractor->X_IMAGE,
          pSextractor->Y_IMAGE,
          limiting_mag_local,
          &temp_magdep_bin,
          &temp_magdep_rms,
          &temp_magdep_magcor,
          1
        );
      }

      // Now apply the color correction to both the output magnitude and to the
      // limiting magnitude. Subtract the extinction correction.

      if (pSpatial_bin->colorflag > 0 && pSextractor->color < 90.0) {
        pSextractor->selected |= COLOR_VALID;

        if (pSpatial_bin->colorflag > 1) {
          pSextractor->selected |= COLOR_METROPOLIS;
        }

        limiting_mag_local -= pSextractor->color * pSpatial_bin->colorterm;
        pSextractor->magcal_local -= pSextractor->color * pSpatial_bin->colorterm;
        pSextractor->magcal_magdep -= pSextractor->color * pSpatial_bin->colorterm;
      }

      // Now apply the extinction correction to both the output magnitude and to
      // the limiting magnitude. Subtract the extinction to make everything
      // brighter.

      if (pExtinction->extinction > 0) {
        limiting_mag_local -= pExtinction->extinction;
        pSextractor->magcal_local -= pExtinction->extinction;
        pSextractor->magcal_magdep -= pExtinction->extinction;
        pSextractor->selected |= EXTINCTION_FLAG;

        if (pSextractor->extinction > 0.0 && printExtinctionError > 0) {
          double extinctionDiff = pSextractor->extinction - pExtinction->extinction;

          if (extinctionDiff < 0) {
            extinctionDiff = -extinctionDiff;
          }

          if (extinctionDiff > 0.001) {
            printf(
              "ERROR: grid extinction %f does not agree with star extinction %f for NUMBER %d of %s\n",
              pExtinction->extinction,
              pSextractor->extinction,
              pSextractor->NUMBER,
              fileroot
            );
            printExtinctionError = 0;

            // We probably have a mismatch between the old extinction algorithm
            // and the algorithm of Jan 19, 2009.  Delete the output file and
            // exit now.

            if (outHandle != NULL) {
              fclose(outHandle);
              outHandle = NULL;
              remove(outfile);
            }

            exit(1);
          }
        }
      }

      if (singleTolerance > EXT_TIME_TOLERANCE) {
        pSextractor->BFLAGS |= (1 << FILTER_BFLAG_BAD_JD);
        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_UNCERTAIN_DATE);
        lowAltitudeCount++;
      } else if (pExtinction->extinction == 0.0 || pExtinction->extinction > extinctionCoefficient * MAX_AIRMASS) {
        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_LOW_ALTITUDE);
        lowAltitudeCount++;
      }

      if (pDmagcor->rejectFlag & DMAG_REJECT_DRAD) {
        // This is a bad spatial bin. Flag everything here as high drad.
        if (pSextractor->AFLAGS < FILTER_AFLAG_BAD) {
          dradCount3++;
        }

        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_DRADBIN);
      }

      // Set additional flags that used to be considered in computemag.c and
      // extract_lightcurves.c.

      if (spatial_bin == 9) {
        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_BIN9);
      }

      if (pSextractor->magcal_iso_rms > MAX_ISO_RMS) {
        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_ISO_RMS);
      }

      if (pSextractor->magcal_local_rms > MAX_LOCAL_RMS) {
        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_LOCAL_RMS);
      }

      if (pSextractor->magcal_local_rms < 0) {
        printf(
          "ERRORD: line %d negative magcal_local_rms %f for NUMBER %d for %s\n",
          __LINE__,
          pSextractor->magcal_local_rms,
          pSextractor->NUMBER,
          fileroot
        );
        pSextractor->magcal_local_rms = 99.0;
      }

      if (limiting_mag_local - pSextractor->magcal_magdep < MAX_LIMITING_MAG) {
        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_LIMITING_MAG);
      }

      if (quality != 0) {
        pSextractor->AFLAGS |= (1 << FILTER_AFLAG_QUALITY);
      } else {
        // If the colorterm was not computed or the colorterm exceeds limits,
        // this is a spatial bin quality issue.

        switch(catalogNumber) {
        case CATALOG_GSC232:
          if ((pSpatial_bin->colorflag != COLORFLAG_METROPOLIS) ||
              (pSpatial_bin->colorterm < MIN_BLUE_COLORTERM) ||
              (pSpatial_bin->colorterm > MAX_BLUE_COLORTERM)) {
            pSextractor->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
            pSextractor->BFLAGS |= (1<<FILTER_BFLAG_COLORTERM);
          }
          break;

        case CATALOG_KEPLER:
          if ((pSpatial_bin->colorflag != COLORFLAG_METROPOLIS) ||
              (pSpatial_bin->colorterm < MIN_KEPLER_BLUE_COLORTERM) ||
              (pSpatial_bin->colorterm > MAX_KEPLER_BLUE_COLORTERM)) {
            pSextractor->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
            pSextractor->BFLAGS |= (1<<FILTER_BFLAG_COLORTERM);
          }
          break;

        case CATALOG_APASS:
          if ((pSpatial_bin->colorflag != COLORFLAG_METROPOLIS) ||
              (pSpatial_bin->colorterm < MIN_APASS_BLUE_COLORTERM) ||
              (pSpatial_bin->colorterm > MAX_APASS_BLUE_COLORTERM)) {
            pSextractor->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
            pSextractor->BFLAGS |= (1<<FILTER_BFLAG_COLORTERM);
          }
          break;

        case CATALOG_GAIA:
          if ((pSpatial_bin->colorflag != COLORFLAG_METROPOLIS) ||
              (pSpatial_bin->colorterm < MIN_GAIA_BLUE_COLORTERM) ||
              (pSpatial_bin->colorterm > MAX_GAIA_BLUE_COLORTERM)) {
            pSextractor->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
            pSextractor->BFLAGS |= (1<<FILTER_BFLAG_COLORTERM);
          }
          break;

        case CATALOG_ATLAS:
          if ((pSpatial_bin->colorflag != COLORFLAG_METROPOLIS) ||
              (pSpatial_bin->colorterm < MIN_ATLAS_BLUE_COLORTERM) ||
              (pSpatial_bin->colorterm > MAX_ATLAS_BLUE_COLORTERM)) {
            pSextractor->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
            pSextractor->BFLAGS |= (1<<FILTER_BFLAG_COLORTERM);
          }
          break;

        case CATALOG_EXPERIMENTAL:
          if ((pSpatial_bin->colorflag != COLORFLAG_METROPOLIS) ||
              (pSpatial_bin->colorterm < MIN_EXPERIMENTAL_BLUE_COLORTERM) ||
              (pSpatial_bin->colorterm > MAX_EXPERIMENTAL_BLUE_COLORTERM)) {
            pSextractor->AFLAGS |= (1<<FILTER_AFLAG_QUALITY);
            pSextractor->BFLAGS |= (1<<FILTER_BFLAG_COLORTERM);
          }
          break;

        default:
          printf("ERROR: unknown catalog %d\n", catalogNumber);
          exit(1);
          break;
        }
      }

      if ((pSextractor->AFLAGS & (1 << FILTER_AFLAG_BLEND_NOMATCH)) == 0) {
        if (limiting_mag_local - pSextractor->magcal_magdep < 0) {
          belowLimitingMagCount++;
          continue;
        }
      } else {
        // Object is buried in a blend. Allow NOMTCH_MAGNITUDE_ERROR for
        // GSC2.3.2 magnitude errors.
        continue;
      }

      outputCount++;

      if ((pSextractor->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) != 0) {
        allLimitingMagRangeCount++;
      }

      // Sanity check: Be sure that every object has a bright bin defined.
      if ((pSextractor->magdep_bin < 0 && pSextractor->magcal_magdep_rms < 99.0) || pSextractor->magcal_magdep == 0) {
        printf(
          "ERROR: Star %d has magdep_bin %d and magcal_magdep %lf magcal_magdep_rms %lf for %s\n",
          pSextractor->NUMBER,
          pSextractor->magdep_bin,
          pSextractor->magcal_magdep,
          pSextractor->magcal_magdep_rms,
          fileroot
        );
        exit(1);
      }

      if (outHandle == NULL) {
        outHandle = fopen(outfile, "wt");

        if (outHandle == NULL) {
          printf("ERROR: Failed to open the output file %s\n", outfile);
          exit(1);
        } else if (verbose) {
          printf("Output file %s\n", outfile);
        }

        fprintf(
          outHandle,
          "REF\t"
          "NUMBER\t"
          "X_IMAGE\t"
          "Y_IMAGE\t"
          "AFLAGS\t"
          "BFLAGS\t"
          "npoints_local\t"
          "MAG_ISO\t"
          "ra\t"
          "dec\t"
          "magcal_iso\t"
          "magcal_iso_rms\t"
          "magcal_local\t"
          "magcal_local_error\t"
          "magcal_local_rms\t"
          "spatial_bin\t"
          "Date\t"
          "limiting_mag_local\t"
          "extinction\t"
          "Stdmag\t"
          "color\t"
          "dra\t"
          "ddec\t"
          "FLUX_ISO\t"
          "MAG_APER\t"
          "MAG_AUTO\t"
          "KRON_RADIUS\t"
          "BACKGROUND\t"
          "THRESHOLD\t"
          "FLUX_MAX\t"
          "THETA_J2000\t"
          "ELLIPTICITY\t"
          "ISOAREA_WORLD\t"
          "FWHM_IMAGE\t"
          "FWHM_WORLD\t"
          "ISO0\t"
          "ISO1\t"
          "ISO2\t"
          "ISO3\t"
          "ISO4\t"
          "ISO5\t"
          "ISO6\t"
          "ISO7\t"
          "plate_dist\t"
          "Blendedmag\t"
          "gsc_bin_index\t"
          "local_bin\t"
          "dradRMS2\t"
          "magdep_bin\t"
          "magcal_magdep\t"
          "magcal_magdep_rms\n"
        );
        fprintf(
          outHandle,
          "---\t"
          "------\t"
          "-------\t"
          "-------\t"
          "------\t"
          "------\t"
          "-------------\t"
          "-------\t"
          "--\t"
          "---\t"
          "----------\t"
          "--------------\t"
          "------------\t"
          "------------------\t"
          "----------------\t"
          "-----------\t"
          "----\t"
          "------------------\t"
          "----------\t"
          "------\t"
          "-----\t"
          "---\t"
          "----\t"
          "--------\t"
          "--------\t"
          "--------\t"
          "-----------\t"
          "----------\t"
          "---------\t"
          "--------\t"
          "-----------\t"
          "-----------\t"
          "-------------\t"
          "----------\t"
          "----------\t"
          "----\t"
          "----\t"
          "----\t"
          "----\t"
          "----\t"
          "----\t"
          "----\t"
          "----\t"
          "----------\t"
          "----------\t"
          "-------------\t"
          "---------\t"
          "--------\t"
          "----------\t"
          "-------------\t"
          "-----------------\n"
        );
      }

      iy = pSextractor->irec % ny;
      ix = pSextractor->irec / ny;

      if (pSextractor->local_bin != ix + nx * iy) {
        printf("ERROR: local_bin mismatch %d %d for %d in %s\n", pSextractor->local_bin, ix + (nx * iy), pSextractor->NUMBER, fileroot);
        exit(1);
      }

      {
        int tmpix;
        int tmpiy;
        PDMAGCOR tmppDmagcor;
        int tmpirec;

        tmpix = (pSextractor->X_IMAGE * x_dmagbins) / (1.0 * mosaicWidth);
        tmpiy = (pSextractor->Y_IMAGE * y_dmagbins) / (1.0 * mosaicHeight);
        tmpirec = tmpiy + y_dmagbins * tmpix;
        tmppDmagcor = &dmagcor_table[tmpirec];

        if (tmppDmagcor->nx != x_dmagbins || tmppDmagcor->ny != y_dmagbins || tmppDmagcor->ix != tmpix || tmppDmagcor->iy != tmpiy) {
          printf(
            "ERROR inconsistent nx %d %d, ny %d %d, ix %d %d, or iy %d %d for %s\n",
            tmppDmagcor->nx,x_dmagbins,
            tmppDmagcor->ny,y_dmagbins,
            tmppDmagcor->ix,tmpix,
            tmppDmagcor->iy,tmpiy,
            dmagcor_name1
          );
          exit (-1);
        }

        if (fabs(tmppDmagcor->errout - pSextractor->magcal_local_error) > 0.001) {
          printf(
            "ERROR: errout %f and magcal_local_error %f are inconsistent for %d in %s\n",
            tmppDmagcor->errout,
            pSextractor->magcal_local_error,
            pSextractor->NUMBER,
            fileroot
          );
          exit(1);
        }

        if (pSextractor->AFLAGS == 0 && pSextractor->magcal_local_error > 90.0) {
          printf(
            "ERROR: AFLAGS %d and magcal_local_error %f are inconsistent for %d in %s\n",
            pSextractor->AFLAGS,
            pSextractor->magcal_local_error,
            pSextractor->NUMBER,
            fileroot
          );
          exit(1);
        }
      }

      pSextractor->spatial_bin = spatial_bin;
      pSextractor->heliocentricJD = heliocentricJD;
      pSextractor->limiting_mag_local = limiting_mag_local;
      pSextractor->extinction = pExtinction->extinction;
      pSextractor->dradRMS2 = pDmagcor->dradRMS2;
      pSextractor->outputFlag = 1;
    }
  }

  // At this point, we want to estimate the rms values for those points which
  // have FILTER_AFLAG_LIMITING_MAG set.  We will rebuild the rms calibration
  // tables using points with FILTER_AFLAG_LIMITING_MAG and a valid catalog
  // magnitude.
  //
  // Initialize the local structures.

  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (pRmsCalibTable[spatial_bin] != NULL) {
      free(pRmsCalibTable[spatial_bin]);
      CalcMemory(1, rmsValueCount[spatial_bin] * sizeof(double), &curMemory, &maxMemory);
      pRmsCalibTable[spatial_bin] = NULL;
    }

    rmsValueCount[spatial_bin] = 0;
    binCalibCount[spatial_bin] = 0;
  }

  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];

    if (
      pSextractor->outputFlag == 0 ||
      pSextractor->Stdmag > 90.0 ||
      pSextractor->magcal_local > 90.0 ||
      (pSextractor->AFLAGS & filterMaskRMSReject) != 0 ||
      (pSextractor->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) == 0
    ) {
      continue;
    }

    goodLimitingMagRangeCount++;
    spatial_bin = pSextractor->spatial_bin;

    if (spatial_bin < 1 || spatial_bin > MAX_SPATIAL_BINS) {
      printf("ERROR: Illegal spatial bin %d for NUMBER %d\n", pSextractor->spatial_bin, pSextractor->NUMBER);
      exit(1);
    }

    if (numBins == 1) {
      /* Use the bin 1 calibration for all stars in single bin mode */
      effective_spatial_bin = 1;
    } else {
      effective_spatial_bin = spatial_bin;
    }

    if (rmsValueCount[effective_spatial_bin] >= rmsValueAlloc[effective_spatial_bin]) {
      // Need to get more spatial bin entries
      rmsValueAlloc[effective_spatial_bin] += RMSVALUE_INITIAL_ALLOC;
      pRmsValue = realloc(pRmsValueTable[effective_spatial_bin], rmsValueAlloc[effective_spatial_bin] * sizeof(RMSVALUE));
      CalcMemory(1, RMSVALUE_INITIAL_ALLOC * sizeof(RMSVALUE), &curMemory, &maxMemory);
      if (pRmsValue == NULL) {
        printf(
          "ERROR: failed to reallocate pRmsValueTable for spatial bin %d size %d\n",
          effective_spatial_bin,
          rmsValueAlloc[effective_spatial_bin]
        );
        exit(1);
      }

      pRmsValueTable[effective_spatial_bin] = pRmsValue;
    }

    pRmsValue = &(pRmsValueTable[effective_spatial_bin][rmsValueCount[effective_spatial_bin]]);
    pRmsValue->magcal_local = pSextractor->magcal_local;
    pRmsValue->diffval = pSextractor->magcal_local - pSextractor->Stdmag;
    rmsValueCount[effective_spatial_bin]++;
  }

  // Estimate the RMS as a function of magnitude for each spatial bin.

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    double *vector;

    if (rmsValueCount[spatial_bin] <= 0) {
      continue;
    }

    qsort((void *) pRmsValueTable[spatial_bin], rmsValueCount[spatial_bin], sizeof(RMSVALUE), RmsValueCompare);

    binCalibCount[spatial_bin] = rmsValueCount[spatial_bin] / RMSCALIB_INCREMENT;
    if (binCalibCount[spatial_bin] == 0) {
      binCalibCount[spatial_bin]++;
    }

    pRmsCalib = (PRMSCALIB) calloc(binCalibCount[spatial_bin], sizeof(RMSCALIB));
    CalcMemory(1, binCalibCount[spatial_bin] * sizeof(RMSCALIB), &curMemory, &maxMemory);
    if (pRmsCalib == NULL) {
      printf("ERROR: failed to allocate pRmsCAlib for spatial bin %d, size %d\n", spatial_bin, binCalibCount[spatial_bin]);
      exit(1);
    }

    pRmsCalibTable[spatial_bin] = pRmsCalib;
    vector = (double *)calloc(rmsValueCount[spatial_bin], sizeof(double));
    CalcMemory(1, rmsValueCount[spatial_bin] * sizeof(double), &curMemory, &maxMemory);
    if (vector == NULL) {
      printf("ERROR: failed to allocate vector for spatial bin %d, size %d\n", spatial_bin, rmsValueCount[spatial_bin]);
      exit(1);
    }

    for (rmsCalibIndex = 0; rmsCalibIndex < binCalibCount[spatial_bin]; rmsCalibIndex++) {
      int minValueIndex = rmsCalibIndex * RMSCALIB_INCREMENT;
      int maxValueIndex = minValueIndex + RMSCALIB_INCREMENT;
      int valueIndex;
      double median;
      double rms;

      pRmsValue = NULL;

      if (rmsCalibIndex == binCalibCount[spatial_bin] - 1) {
        maxValueIndex = rmsValueCount[spatial_bin];
      }

      for (valueIndex = minValueIndex; valueIndex < maxValueIndex; valueIndex++) {
        pRmsValue = &(pRmsValueTable[spatial_bin][valueIndex]);
        vector[valueIndex - minValueIndex] = pRmsValue->diffval;
      }

      if (CalcMedianAndRMS(maxValueIndex - minValueIndex, 0, vector, &median, &rms, 1, 3.0, 1) > 0) {
        pRmsCalib = &(pRmsCalibTable[spatial_bin][rmsCalibIndex]);
        pRmsCalib->rms = rms;
        pRmsCalib->magcal_local = pRmsValue->magcal_local;
      } else {
        printf(
          "ERROR: Failed to calculate a rms curve for spatial bin %d number of stars %d for %s\n",
          spatial_bin,
          maxValueIndex - minValueIndex,
          fileroot
        );
        exit(1);
      }
    }

    free(vector);
    CalcMemory(-1, rmsValueCount[spatial_bin] * sizeof(double), &curMemory, &maxMemory);
    vector = NULL;
  }

  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (pRmsValueTable[spatial_bin] != NULL) {
      free(pRmsValueTable[spatial_bin]);
      CalcMemory(-1, rmsValueAlloc[spatial_bin] * sizeof(RMSVALUE), &curMemory, &maxMemory);
      pRmsValueTable[spatial_bin] = NULL;
    }
  }

  // Now go back to all of our selected stars and insert the new RMS value. We
  // do not use filterMaskRMSReject here because the image can be assigned a
  // valid RMS even if it does not qualify for RMS calibration.

  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];

    if (
      pSextractor->outputFlag == 0 ||
      pSextractor->Stdmag > 90.0 ||
      pSextractor->magcal_local > 90.0 ||
      pSextractor->magcal_local_rms < 90.0 ||
      (pSextractor->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) == 0
    ) {
      continue;
    }

    secondLimitingMagRangeCount++;

    spatial_bin = pSextractor->spatial_bin;
    if (spatial_bin < 1 || spatial_bin > MAX_SPATIAL_BINS) {
      printf("ERROR: Illegal spatial bin %d for NUMBER %d\n", pSextractor->spatial_bin, pSextractor->NUMBER);
      exit(1);
    }

    if (numBins == 1) {
      /* Use the bin 1 calibration for all stars in single bin mode */
      effective_spatial_bin = 1;
    } else {
      effective_spatial_bin = spatial_bin;
    }

    pSextractor->magcal_local_rms = CalcRMS_E(
      binCalibCount,
      pRmsCalibTable,
      effective_spatial_bin,
      pSextractor->magcal_local
    );

    if (pSextractor->magcal_local_rms > MAX_LOCAL_RMS) {
      pSextractor->AFLAGS |= (1 << FILTER_AFLAG_LOCAL_RMS);
    } else {
      pSextractor->AFLAGS &= ~(1 << FILTER_AFLAG_LOCAL_RMS);
    }
  }

  // Finally ready to emit!

  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];

    if (pSextractor->outputFlag == 1) {
      fprintf(
        outHandle,
        "%s\t"
        "%d\t"
        "%f\t"
        "%f\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%d\t"
        "%.6f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%f\t"
        "%f\t"
        "%d\t"
        "%d\t"
        "%f\t"
        "%d\t"
        "%f\t"
        "%f\n",
        pSextractor->REF,
        pSextractor->NUMBER,
        pSextractor->X_IMAGE,
        pSextractor->Y_IMAGE,
        pSextractor->AFLAGS,
        pSextractor->BFLAGS | ((pSextractor->selected) << CAL_FLAG_SHIFT),
        pSextractor->npoints_local,
        pSextractor->MAG_ISO,
        pSextractor->ra,
        pSextractor->dec,
        pSextractor->magcal_iso,
        pSextractor->magcal_iso_rms,
        pSextractor->magcal_local,
        pSextractor->magcal_local_error,
        pSextractor->magcal_local_rms,
        pSextractor->spatial_bin,
        pSextractor->heliocentricJD,
        pSextractor->limiting_mag_local,
        pSextractor->extinction,
        pSextractor->Stdmag,
        pSextractor->color,
        pSextractor->dra,
        pSextractor->ddec,
        pSextractor->FLUX_ISO,
        pSextractor->MAG_APER,
        pSextractor->MAG_AUTO,
        pSextractor->KRON_RADIUS,
        pSextractor->BACKGROUND,
        pSextractor->THRESHOLD,
        pSextractor->FLUX_MAX,
        pSextractor->THETA_J2000,
        pSextractor->ELLIPTICITY,
        pSextractor->ISOAREA_WORLD,
        pSextractor->FWHM_IMAGE,
        pSextractor->FWHM_WORLD,
        pSextractor->ISO0,
        pSextractor->ISO1,
        pSextractor->ISO2,
        pSextractor->ISO3,
        pSextractor->ISO4,
        pSextractor->ISO5,
        pSextractor->ISO6,
        pSextractor->ISO7,
        pSextractor->plate_dist,
        pSextractor->Blendedmag,
        pSextractor->gsc_bin_index,
        pSextractor->local_bin,
        pSextractor->dradRMS2,
        pSextractor->magdep_bin,
        pSextractor->magcal_magdep,
        pSextractor->magcal_magdep_rms
      );
    }
  }

  // Write out the rejectFlag and dradAverage with for the dmagcor.grid file

  if (dmagcor_handle1 != NULL) {
    Close(dmagcor_handle1);
    dmagcor_handle1 = NULL;
  }

  if (dmagcor_header != NULL) {
    table_hdrfree(dmagcor_header);
    dmagcor_header = NULL;
  }

  pDmagcor = NULL;

  dmagcor_handle2 = fopen(dmagcor_name2, "wt");

  if (dmagcor_handle2 == NULL) {
    printf("ERROR: Failed to open %s for writing\n", dmagcor_name2);
    exit(1);
  }

  if (verbose) {
    printf("Writing updated %s\n", dmagcor_name2);
  }

  fprintf(
    dmagcor_handle2,
    "nx\t"
    "ny\t"
    "ix\t"
    "iy\t"
    "zout\t"
    "errout\t"
    "npout\t"
    "rejectFlag\t"
    "drad_bin_count\t"
    "drad_bin_size\t"
    "drad_reject_count\t"
    "draMedian\t"
    "draRMS\t"
    "ddecMedian\t"
    "ddecRMS\t"
    "drad_bin_count2\t"
    "drad_bin_size2\t"
    "drad_reject_count2\t"
    "draMedian2\t"
    "draRMS2\t"
    "ddecMedian2\t"
    "ddecRMS2\t"
    "dradRMS2\t"
    "bright_magcor_local\t"
    "bright_magcor_error\t"
    "bright_npoints_local\n"
  );
  fprintf(
    dmagcor_handle2,
    "--\t"
    "--\t"
    "--\t"
    "--\t"
    "----\t"
    "------\t"
    "-----\t"
    "----------\t"
    "--------------\t"
    "-------------\t"
    "-----------------\t"
    "---------\t"
    "------\t"
    "----------\t"
    "-------\t"
    "---------------\t"
    "--------------\t"
    "------------------\t"
    "----------\t"
    "-------\t"
    "-----------\t"
    "--------\t"
    "--------\t"
    "------------------\t"
    "-------------\t"
    "----------------\n"
  );

  irec = 0;

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      pDmagcor = &dmagcor_table[irec];

      fprintf(
        dmagcor_handle2,
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%.6f\t"
        "%.6f\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%f\t"
        "%.6f\t"
        "%.6f\t"
        "%d\n",
        pDmagcor->nx,
        pDmagcor->ny,
        pDmagcor->ix,
        pDmagcor->iy,
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
        pDmagcor->bright_npoints_local
      );

      irec++;
    }
  }

  if (dmagcor_handle2 != NULL) {
    fclose(dmagcor_handle2);
  }

  if (dmagcor_table != NULL) {
    Free(dmagcor_table);
    dmagcor_table = NULL;
  }

  if (dmagcor_header != NULL) {
    table_hdrfree(dmagcor_header);
    dmagcor_header = NULL;
  }

  pDmagcor = NULL;

  // All done! Clean up.

  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (pRmsCalibTable[spatial_bin] != NULL) {
      free(pRmsCalibTable[spatial_bin]);
      pRmsCalibTable[spatial_bin] = NULL;
    }
  }

  if (spatial_bin_table != NULL) {
    free(spatial_bin_table);
  }

  if (spatial_bins_header != NULL) {
    table_hdrfree(spatial_bins_header);
  }

  if (spatial_bins_handle != NULL) {
    Close(spatial_bins_handle);
  }

  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (maggrid_handle[spatial_bin] != NULL) {
      Close(maggrid_handle[spatial_bin]);
    }

    if (maggrid_header[spatial_bin] != NULL) {
      table_hdrfree(maggrid_header[spatial_bin]);
    }

    if (maggrid_table[spatial_bin] != NULL) {
      Free(maggrid_table[spatial_bin]);
    }
  }

  if (extinction_table != NULL) {
    Free(extinction_table);
  }

  if (extinction_header != NULL) {
    table_hdrfree(extinction_header);
  }

  if (extinction_handle != NULL) {
    Close(extinction_handle);
  }

  if (sextractor_table != NULL) {
    Free(sextractor_table);
  }

  if (sextractor_header != NULL) {
    table_hdrfree(sextractor_header);
  }

  if (sextractor_handle != NULL) {
    Close(sextractor_handle);
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }

  if (pMagdepTable != NULL) {
    free(pMagdepTable);
  }

  FreeMagdepSubarrays(pMagdepLimits);

  time(&curTime);
  curTime -= startTime;

  fprintf(
    stdout,
    "Illegal %d missing %d too dim %d too bright %d highZout %d too low %d no local %d"
    " noRef %d output %d %f single bin %d numBins %d blendCorr %d blendUncorr %d "
    "drad %d %d %d %d belowLimit %d seconds %ld maxMemory %lld for %s%s\n",
    illegalBinCount,
    missingBinCount,
    tooDimCount,
    tooBrightCount,
    highZoutCount,
    lowAltitudeCount,
    noLocalPoints,
    noRefCount,
    outputCount,
    (1.0 * noRefCount) / outputCount,
    singleBinCount,
    numBins,
    blendCorrectedCount,
    blendTotalCount - blendCorrectedCount,
    dradCount1,
    dradCount2,
    dradCount3,
    dradCount4,
    belowLimitingMagCount,
    curTime,
    maxMemory,
    fileroot,
    qualifier
  );

  fprintf(
    stdout,
    "allLimitingMagRangeCount %d goodLimitingMagRangeCount %d secondLimitingMagRangeCount %d for %s%s\n",
    allLimitingMagRangeCount,
    goodLimitingMagRangeCount,
    secondLimitingMagRangeCount,
    fileroot,
    qualifier
  );

  return 0;
}

/* Jan 22, 2008 Edward J. Los - Initial version
 * Feb 15, 2008 Edward J. Los - add verbose flag
 *                              Allow zero slope in the calibration curve
 *                              Handle round-off errors in annular9.m by allowing
 *                              data slightly beyond the limiting magnitude of the calibration curve
 *                              Relax the check on the monotonicity of the calibration curve.
 * Feb 18, 2008 Edward J. Los - Correct pass 2 by adding a qualifier
 *                            - Add Heliocentric Julian Day support
 * Feb 25, 2008 Edward J. Los - Add single-bin support: Take only those points within the range
 *                              of properly calibrated points.
 *                              If the standard deviation for local calibration is above 90, disable
 *                              local calibration for this star
 * Mar  5, 2008 Edward J. Los - Add the limiting_mag to the output table
 *                              Add the number of bins as a qualifier
 *                              Add extinction
 * Mar 10, 2008 Edward J. Los - Correct single bin case error introduced on Mar 5.  Only stars in the
 *                              original bin 1 were being output.
 * Mar 28, 2008 Edward J. Los - Use ra_2 and dec_2 instead of ra_1 and dec_1
 * Mar 30, 2008 Edward J. Los - For GSC2.3.2 stars, add Stdmag, dra and ddec to help look for
 *                              unusual stars and endure that they are not double exposures.
 * Apr 22, 2008 Edward J. Los - Add the following fields: FLUX_ISO MAG_APER MAG_AUTO KRON_RADIUS BACKGROUND THRESHOLD
 *                              FLUX_MAX THETA_J2000 ELLIPTICITY ISOAREA_WORLD FWHM_IMAGE FWHM_WORLD ISO0 ISO1 ISO2
 *                              ISO3 ISO4 ISO5 ISO6 ISO7 plate_dist
 * Apr 23, 2008 Edward J. Los - Fix extinction and dmgacor error for objects at the edge of the plate.
 *                              Recover FLAGS from the match or local calibration tables
 * Apr 28, 2008 Edward J. Los - Output the plate RA and DEC, not the GSC2.3.2 RA and DEC for all stars.
 *                              Ouptut dra and ddec for all matched stars
 * Apr 30, 2008 Edward J. Los - Fix Stdmag for stars which made it through the pipeline.
 * Jun 10, 2008 Edward J. Los - Use the plate defect file to update the Sextractor FLAGS value.
 * Jun 13, 2008 Edward J. Los - Add new algorithm (RMS_ALGORITHM_B) for calculating magcal_local_rms
 * Jun 18, 2008 Edward J. Los - Add new algorithm (RMS_ALGORITHM_C) for calculating magcal_local_rms
 * Jul  1, 2008 Edward J. Los - Add new algorithm (RMS_ALGORITHM_D) for calculating magcal_local_rms
 * Jul  7, 2008 Edward J. Los - Read the match_*_n.db file to get entries that are flagged FILTER_FLAG_BLEND_NOMATCH
 *                              add these to the table but with phoney Sextractor NUMBER values
 * Jul  8, 2008 Edward J. Los - Add FILTER_FLAG_HIZOUT for a dmag smoothing sanity check.
 * Jul 28, 2008 Edward J. Los - Add FILTER_FLAG_TOO_BRIGHT for objects that are to bright to be calibrated.
 * Aug 11, 2008 Edward J. Los - Add magnitude adjustment for blended stars and use FILTER_FLAG_ADJUST_BLEND to flag these stars
 *                              Add new algorithm (RMS_ALGORITHM_E) for calculating magcal_local_rms
 *                              Rewrite the <mosaic>_dmagcor.grid file so the plotting routines can see the rejectFlag
 *                              Read in the _drad.db file and if the bin drad is too high, flag every point in the bin
 *                              as a high drad point.
 *                              Reject all stars in a local calibration bin if the number of high drad stars exceeds
 *                               the number of low drad stars (all test stars must be 0.5 mag brighter than the limiting mag)
 * Aug 18, 2008 Edward J. Los - Correct the blended magnitude adjustment to skip the adjustment if the resulting flux is negative
 * Aug 22, 2008 Edward J. Los - Correct the high drad local calibration bin rejection algorithm by including matched
 *                              objects that did not pass through the pipeline use the FLUX_MAX > THRESHOLD_FACTOR * THRESHOLD
 *                              Back out correction of case b blended magnitudes, but save the Blendedmag field
 *                              Compute the mean local bin drad value, dradAverage, and save it in the local bin table, and add it to the output file
 * Oct  1, 2008 Edward J. Los - Read the saturationFlux from the drad file.  Do not consider objects with high flux for drad processing.
 * Oct 20, 2008 Edward J. Los - Reject local bins if the errout is too high
 * Nov 24, 2008 Edward J. Los - Back out saturation flux
 * Nov 25, 2008 Edward J. Los - Replace MaxDradPixels(series,plateNumber) with MaxDradPixels() function
 * Dec 19, 2008 Edward J. Los - If CalculateLowessMagnitude fails, do not attempt to apply a local
 *                              correction.
 *                              Add MAX_BRIGHT_ADJUST to the max_bright_iso to avoid CalculateLowessMagnitude failures
 *                              near the beginning of the bright range of the star.
 * Dec 22, 2008 Edward J. Los - Flag objects with high extinction or missing extinction
 * Jan 12, 2009 Edward J. Los - Add colorterm support
 * Jan 27, 2009 Edward J. Los - split FLAGS into AFLAGS and BFLAGS.  The CAL_FLAG now becomes part of BFLAGS
 *                              Correct a bug in setting the defect bit.
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 * Feb 18, 2009 Edward J. Los - Change _dmagcor.grid to _dmagcor.db
 *                              Back out saturation flux
 *                              Implement new drad algorithm (NEW_BLEND_DEFINITION); remove dradAverage
 *                              add local_bin_index and gsc_bin_index
 *                              Set AFLAGS bits for conditions that used to be filtered in computmag.c and extract_lightcurves.c
 * Feb 27, 2009 Edward J. Los - Add CalcMemory
 * Mar  2, 2009 Edward J. Los - Add dradRMS2 to allobjects.db
 * Mar 17, 2009 Edward J. Los - Correct gsc_bin_index for nomtch stars
 * Mar 30, 2009 Edward J. Los - Add COLOR_RATIO to the colorterm calibration
 * May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 * May 19, 2009 Edward J. Los - Avoid reallocating the large sextractor table.
 *                              Do not output images that are below the limiting magnitude
 * Jun 30, 2009 Edward J. Los - Remove FILTER_AFLAG_BLEND_NOMATCH functionality
 * Aug 20, 2009 Edward J. Los - Add Kepler Input Catalog support
 * Sep  3, 2009 Edward J. Los - Correct drad file name for the kepler case
 * Oct 21, 2009 Edward J. Los - Create the output file only if we have data to write
 * Nov 16, 2009 Edward J. Los - Add multiple exposure support, particularly with respect to getting the JulianDate
 *                              If the JulianDate is not accurate, flag it here rather than in prepare_octave
 * Dec  3, 2009 Edward J. Los - Correct garbage Blendedmag entries
 * Jan 19, 2010 Edward J. Los - Implement a new extinction formula using GetExtinctionCoefficient
 * Jan 22, 2010 Edward J. Los - Add plate quality support
 * Apr 17, 2010 Edward J. Los - Tighten the tolerance for rejecting bad spatial bins by using a 3 pixel limit for all
 *                              plate series.
 * Apr 22, 2010 Edward J. Los - Turn off the stale extinction error temporarily
 * Feb 28, 2011 Edward J. Los - add bright entries to the _dmagcor.db file
 * Mar  7, 2011 Edward J. Los - Do not use colorterm processing for colors outside of reasonable bounds
 * Mar 14, 2011 Edward J. Los - make local bin checking mandatory
 *                              avoid using irec for extinction and dmag file indices
 *                              Correct local_bin_index error in the output file
 * Apr 20, 2011 Edward J. Los - Add the qualifier to the final report
 * Aug 29, 2011 Edward J. Los - Skip false ERROR: Stdmag 0.000000 in input file does not agree with blend Stdmag when it occurs
 *                              because filterblended attempted to write out a REF twice.
 * Sep  7, 2011 Edward J. Los - Back out deprecated COLOR_RATIO
 * Sep 17, 2011 Edward J. Los - Adjust color limits for B-V apass catalog
 * Oct 30, 2011 Edward J. Los - Add magnitude-dependent star correction support (FILTER_BFLAG_MAGDEP_MAGCOR)
 *                              Introduce magcal_magdep: magnitude-corrected local magnitude.
 *                              New output values:  magdep_bin
 *                                                  magcal_magdep
 *                                                  magcal_magdep_rms
 *                              In applying the correction, use magcal_iso which should represent not Stdmag, but Stdmag with
 *                                 extinction and color corrections applied prior to running annular9.m
 * Nov 20, 2011 Edward J. Los - Correct error in which extinction and color correction was not applied to magcal_magdep
 * Feb 27, 2012 Edward J. Los - Add UpdateQuality to look for multiple exposures
 *                              Ignore errors reading the blend table, assuming that there were no blends to be read.
 * Apr  6, 2012 Edward J. Los - Add diagnostic printout for individual reference objects
 * Jul  3, 2012 Edward J. Los - Add high background object support
 * Jul  9, 2012 Edward J. Los - Remove THRESHOLD_FACTOR for the new defect filter of  Fri 6/29/12 9:43 AM
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Aug 13, 2012 Edward J. Los - Avoid corrupted errout (magcal_local_error) and ensure that FILTER_AFLAG_HIZOUT is set.
 * Aug 16, 2012 Edward J. Los - Distinguish between input and output dmagcor.db files
 * Nov  9, 2012 Edward J. Los - Correct colorterm limit code.
 * Jan  7, 2013 Edward J. Los - Correct background file name for multiple exposures.
 * Mar  8, 2013 Edward J. Los - Correct dra
 * Apr 16, 2013 Edward J. Los - Count the stars in the MAX_LIMITING_MAG range
 *                              Eliminate references to REMOVE_NOMATCH
 * Jun 24, 2013 Edward J. Los - Make an unreadable, but not a missing background file benign, because the background
 *                              file may contain no points.
 * Jul 30, 2013 Edward J. Los - Correct RMS_ALGORITHM_E to use magcal_local instead of magcal_magdep
 *                              Reuse RMS_ALGORITHM_E to estimate the rms of all stars with FILTER_AFLAG_LIMITING set
 * Aug  3, 2013 Edward J. Los - use FILTER_MASK_RMS_REJECT to eliminate blends, variables, and non-stars from the rms
 *                              estimate for stars with FILTER_AFLAG_LIMITING set.
 * Sep  3, 2013 Edward J. Los - Do not use high proper motion error objects to compute drad
 * Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
 * Jun  6, 2014 Edward J. Los - Do not perform the stale extinction check if all of the stars will be marked below the horizon.
 * Sep 15, 2014 Edward J. Los - Flag negative values of magcal_local_rms
 * Dec  8, 2014 Edward J. Los - Add CheckMaggridTable to set negative griderr values to the previous peak error on the dim side of the curve
 * Dec 10, 2014 Edward J. Los - Modify CheckMaggridTable to handle negative griderr values on the bright side of the curve
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * Nov 28, 2016 Edward J. Los - Correct bug in CalcRMS_E for a magcal_local equal to a pRmsCalibTable entry
 * Jul 26, 2017 Edward J. Los - Handle 4"x5" plates.   Correct Dmagcor consistency check for nx != ny
 * May 28, 2018 Edward J. Los - Support GAIA
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 */
