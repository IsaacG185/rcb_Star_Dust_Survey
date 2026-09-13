// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// dasch-querycat -- query a DASCH binary catalog
//
// Derived from matchstars.c.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Starlink -- needed for File, Seek, etc.
#include <table.h>

// libwcs
#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

// local
#include "pipelineutils.h"

#define MAX_BUFFER 512

extern GSCBIN gscBin64;


typedef struct _gscimagebin {
  GSCIMAGE gscimage;
  int gsc_bin_index;
} GSCIMAGEBIN;

// This structure is silly for this particular program, but it's what we've
// inherited from matchstars, and there's value in maintaining comparability
// with that code.
typedef struct _inputdata {
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
} INPUTDATA;

typedef struct _cacherow {
  int numEntries;  /* Number of entries in this row */
  int curIndex;    /* Current index into this row */
  INPUTDATA *inputdata; /*  row contents */
} CACHEROW;


static GSCBIN *pGscBin = &gscBin64;
static double search_radius_deg = 20. / 3600.;  // default: 20 arcsec


/*
 * Return 0 if the gsc RA is within range,
 *        1 if the gsc RA is higher
 *       -1 if the gsc RA is lower
 */
static inline int CheckRA(
  GSCIMAGEBIN *pGscImage,
  double minImageRa,
  double maxImageRa
) {
  double gscRa = pGscImage->gscimage.ra;

  minImageRa -= pGscImage->gscimage.RaSigmaPM;
  maxImageRa += pGscImage->gscimage.RaSigmaPM;

  if (minImageRa < 0.0) {
    /* Here the low end of the range croses the first point of Aries */
    if (gscRa >= minImageRa + 360.0 || gscRa <= maxImageRa) {
      return 0;
    }

    if (gscRa - maxImageRa < minImageRa + 360.0 - gscRa) {
      return 1;
    }

    return -1;
  } else if (maxImageRa >= 360.0) {
    /* Here the upper end of the range crosses the first point of Aries */
    if (gscRa >= minImageRa || gscRa <= (maxImageRa - 360.0)) {
      return 0;
    }

    if (minImageRa - gscRa < gscRa - (maxImageRa - 360.0)) {
      return -1;
    }

    return 1;
  }

  if (gscRa < minImageRa) {
    return -1;
  } else if (gscRa > maxImageRa) {
    return 1;
  }

  return 0;
}


/* Sort routine based on the right ascension */
static int
GscImageCompare(const void *first, const void *second)
{
  double numberFirst = ((GSCIMAGEBIN *) first)->gscimage.ra;
  double numberSecond = ((GSCIMAGEBIN *) second)->gscimage.ra;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  }

  return 0;
}


// Returns the number of stars in this bin
static int
ReadCatalogDecBin(
  File catalog_handle,
  File index_handle,
  STARINDEX **pStarIndex,
  int *pStarIndexAlloc,
  int *pStarIndexEntries,
  GSCIMAGEBIN **pGscImage,
  GSCIMAGE **pGscStorage,
  int *pGscImageAlloc,
  int *pGscImageEntries,
  int curCatalogBin,
  double *minRa,
  double *maxRa,
  double limitingMagnitude,
  double epoch_delta,
  int n_ra_chunks
) {
  STARINDEX *pCurStarIndex = *pStarIndex;
  STARINDEX *tempStarIndex;
  int raIndex;
  int firstBin;
  int lastBin;
  int numBins;
  int totalBins = 0;
  int decBin;
  int raBin;
  double declination = (curCatalogBin * pGscBin->bin_size) - 90.0;
  int readItems;
  int totalStars = 0;
  int tempStars;
  int numStars;
  int catIndex;
  int starIndex;
  int starIndex2;
  GSCIMAGEBIN *pCurGscImage = *pGscImage;
  GSCIMAGE *pCurGscStorage = *pGscStorage;
  GSCIMAGE *tempGscStorage = NULL;
  GSCIMAGEBIN *tempGscImage;
  GSCIMAGEBIN *tempGscImage2;
  double factor;
  double raPrecession = 0.0;
  double tmpRa;
  int precessionFlag = 0;

  factor = cos(DEGREES_TO_RAD * declination);

  if (factor != 0.0) {
    raPrecession = PROPER_MOTION_DEGREES / factor;
  } else {
    raPrecession = PROPER_MOTION_DEGREES;
  }

  *pStarIndexEntries = 0;
  *pGscImageEntries = 0;

  for (raIndex = 0; raIndex < n_ra_chunks; raIndex++) {
    // Find out which star index bins we need to read
    tmpRa = minRa[raIndex] - raPrecession - pGscBin->bin_size;

    if (tmpRa < 0.0) {
      tmpRa = 0.0;
    }

    firstBin = GetGSCBin(pGscBin, tmpRa, declination, &decBin, &raBin, "dasch-querycat");
    GetSubbins(pGscBin, firstBin, &raBin, &decBin, "dasch-querycat");
    tmpRa = maxRa[raIndex] + raPrecession + pGscBin->bin_size;

    if (tmpRa >= 360.0) {
      tmpRa = 359.9999;
    }

    lastBin = GetGSCBin(pGscBin, tmpRa, declination, &decBin, &raBin, "dasch-querycat");
    numBins = lastBin - firstBin + 1;

    if (numBins <= 0) {
      fprintf(stderr, "error: numBins is %d in ReadCatalogDecBin\n", numBins);
      exit(1);
    }

    if ((numBins + totalBins) >= *pStarIndexAlloc) {
      *pStarIndexAlloc = numBins + totalBins + 100;
      tempStarIndex = realloc(*pStarIndex, (*pStarIndexAlloc) * sizeof(STARINDEX));

      if (tempStarIndex == NULL) {
        fprintf(stderr, "error: allocating star index of size %d\n", *pStarIndexAlloc);
        exit(1);
      }

      *pStarIndex = tempStarIndex;
      pCurStarIndex = &tempStarIndex[totalBins];
    }

    totalBins += numBins;
    Seek(index_handle, firstBin * sizeof(STARINDEX), SEEK_SET);
    readItems = Read(index_handle, pCurStarIndex, sizeof(STARINDEX), numBins);

    if (readItems != numBins) {
      fprintf(stderr, "error: reading star index file\n");
      exit(1);
    }

    /* Figure out how many bins we need to read */
    numStars = 0;

    for (catIndex = 0; catIndex < numBins; catIndex++) {
      tempStarIndex = &pCurStarIndex[catIndex];
      numStars += tempStarIndex->numStars;
    }

    if (numStars == 0) {
      continue;
    }

    if ((totalStars + numStars) >= *pGscImageAlloc) {
      *pGscImageAlloc = totalStars + numStars + 100;
      tempGscImage = realloc(*pGscImage, (*pGscImageAlloc) * sizeof(GSCIMAGEBIN));

      if (tempGscImage == NULL) {
        fprintf(stderr, "error: allocating gsc index of size %d\n", *pGscImageAlloc);
        exit(1);
      }

      *pGscImage = tempGscImage;
      pCurGscImage = &tempGscImage[totalStars];

      tempGscStorage = realloc(*pGscStorage, (*pGscImageAlloc) * sizeof(GSCIMAGE));

      if (tempGscStorage == NULL) {
        fprintf(stderr, "error: allocating pCurGscStorage of size %d\n", *pStarIndexAlloc);
        exit(1);
      }

      *pGscStorage = tempGscStorage;
      pCurGscStorage = tempGscStorage;
    }

    /* Go to the first bin in the index */

    Seek(catalog_handle, pCurStarIndex->offset, SEEK_SET);
    totalStars += numStars;
    /* Now read the required number of stars */
    readItems = Read(catalog_handle, pCurGscStorage, sizeof(GSCIMAGE), numStars);

    if (readItems != numStars) {
      fprintf(stderr, "error: reading GSC catalog file\n");
      exit(1);
    }

    /* Now copy everything into the larger GSCIMAGEBIN area */
    for (starIndex = 0; starIndex < numStars; starIndex++) {
      tempGscImage = &pCurGscImage[starIndex];
      tempGscStorage = &pCurGscStorage[starIndex];
      memcpy(&tempGscImage->gscimage, tempGscStorage, sizeof(GSCIMAGE));
      tempGscImage->gsc_bin_index = 0;
    }

    pCurStarIndex = &(*pStarIndex)[totalBins];
    pCurGscImage = &(*pGscImage)[totalStars];
  }

  if (limitingMagnitude != 0.0) {
    /* If we have a limiting magnitude, get rid of all dimmer
       stars */
    tempStars = totalStars;
    totalStars = 0;
    pCurGscImage = *pGscImage;

    for (starIndex = 0; starIndex < tempStars; starIndex++) {
      tempGscImage2 = &pCurGscImage[starIndex];

      if (tempGscImage2->gscimage.Stdmag <= limitingMagnitude) {
        if (starIndex != totalStars) {
          tempGscImage = &pCurGscImage[totalStars];
          memcpy(tempGscImage, tempGscImage2, sizeof(GSCIMAGEBIN));
        }

        totalStars++;
      }
    }
  }

  /* Get rid of all Stdmag == 0 stars.  These appear to be
      catalog errors */
  tempStars = totalStars;
  totalStars = 0;
  pCurGscImage = *pGscImage;

  for (starIndex = 0; starIndex < tempStars; starIndex++) {
    tempGscImage2 = &pCurGscImage[starIndex];

    if (tempGscImage2->gscimage.Stdmag != 0) {
      if (starIndex != totalStars) {
        tempGscImage = &pCurGscImage[totalStars];
        memcpy(tempGscImage, tempGscImage2, sizeof(GSCIMAGEBIN));
      }
      totalStars++;
    }
  }

  /* Now save the gsc bin indices in the gsc bin table prior to changing
     the ra and dec because of precession */
  pCurGscImage = *pGscImage;

  for (starIndex = 0; starIndex < totalStars; starIndex++) {
    tempGscImage = &pCurGscImage[starIndex];
    tempGscImage->gsc_bin_index = GetGSCBin(
      pGscBin,
      tempGscImage->gscimage.ra,
      tempGscImage->gscimage.dec,
      &decBin,
      &raBin,
      "dasch-querycat"
    );
  }

  // Apply proper motions, and overload {Ra,Dec}SigmaPM to express RA/Dec
  // uncertainty.

  pCurGscImage = *pGscImage;

  for (starIndex = 0; starIndex < totalStars; starIndex++) {
    tempGscImage2 = &pCurGscImage[starIndex];

    if (tempGscImage2->gscimage.RaPM != 0 || tempGscImage2->gscimage.DecPM != 0) {
      precessionFlag = 1;
      tempGscImage2->gscimage.dec += (tempGscImage2->gscimage.DecPM * epoch_delta) / (3600.0 * 1000.0);

      if ((tempGscImage2->gscimage.dec < 90.0) && (tempGscImage2->gscimage.dec > -90.0)) {
        tempGscImage2->gscimage.ra += ((tempGscImage2->gscimage.RaPM * epoch_delta) / (3600.0 * 1000.0)) / cos(tempGscImage2->gscimage.dec * DEGREES_TO_RAD);
      }
    }

    if (tempGscImage2->gscimage.RaSigmaPM != 0.0 || tempGscImage2->gscimage.DecSigmaPM != 0.0) {
      double factor2 = cos(DEGREES_TO_RAD * tempGscImage2->gscimage.dec);

      if (factor2 == 0.0) {
        factor2 = 1.0;
      }

      tempGscImage2->gscimage.RaSigmaPM *= fabs((PROPER_MOTION_SIGMA * epoch_delta) / (3600.0 * 1000.0 * factor2)); /* Convert to degrees in R.A. */
      tempGscImage2->gscimage.DecSigmaPM *= fabs((PROPER_MOTION_SIGMA * epoch_delta) / (3600.0 * 1000.0)); /* Convert to degrees */
    }
  }

  if (precessionFlag) {
    /* We need to resort this bin so that the RA's are back in order */
    qsort((void*)pCurGscImage, totalStars, sizeof(GSCIMAGEBIN), GscImageCompare);
  }

  if (minRa[1] <= maxRa[1]) {
    /* There may be duplicates here that we need to remove */
    if (precessionFlag == 0) {
      qsort((void*)pCurGscImage, totalStars, sizeof(GSCIMAGEBIN), GscImageCompare);
    }

    starIndex = 1;

    while (starIndex < totalStars) {
      tempGscImage = &pCurGscImage[starIndex-1];
      tempGscImage2 = &pCurGscImage[starIndex];

      if (tempGscImage->gscimage.REFNumber == tempGscImage2->gscimage.REFNumber) {
        /* Duplicate entry - eliminate it, but keep sorting order */

        for (starIndex2 = starIndex; starIndex2 < (totalStars-1); starIndex2++) {
          tempGscImage = &pCurGscImage[starIndex2];
          tempGscImage2 = &pCurGscImage[starIndex2+1];
          memcpy(tempGscImage, tempGscImage2, sizeof(GSCIMAGEBIN));
        }

        totalStars--;
      } else {
        starIndex++;
      }
    }
  }

  *pStarIndexEntries = totalBins;
  *pGscImageEntries = totalStars;
  return totalStars;
}


static void
output_match(
  GSCIMAGEBIN *pGscImage,
  INPUTDATA *pInput,
  double limitingMagnitude,
  double epoch_delta
) {
  double ddec;
  double dra;
  double factor;
  double RaSigmaPM;
  double DecSigmaPM;
  char ref_text[MAX_REF];
  double deltara;

  /* If the user requested a limiting magnitude,
     do not output match if we are above the limit */
  if (limitingMagnitude != 0.0 && pGscImage->gscimage.Stdmag > limitingMagnitude) {
    return;
  }

  // Calculate position difference between this source and the query center, in
  // arcsec.

  factor = cos(DEGREES_TO_RAD * ((pInput->dec + pGscImage->gscimage.dec) / 2.0));
  deltara = pInput->ra - pGscImage->gscimage.ra;

  if (deltara < -180.0) {
    deltara = 360. + deltara;
  } else if (deltara > 180.0) {
    deltara = deltara - 360.0;
  }

  dra = 3600 * factor * deltara;
  ddec = 3600 * (pInput->dec - pGscImage->gscimage.dec);

  // During processing the gscimage values were overloaded to represent ra/dec uncertainties.
  // Undo that, to get us back to the 2-sigma proper motion uncerts in mas/yr.

  RaSigmaPM = fabs(
    (3600.0 * 1000.0 * factor / (PROPER_MOTION_SIGMA * epoch_delta)) * pGscImage->gscimage.RaSigmaPM
  );
  DecSigmaPM = fabs(
    (3600.0 * 1000.0 / (PROPER_MOTION_SIGMA * epoch_delta)) * pGscImage->gscimage.DecSigmaPM
  );

  GetREF(pGscImage->gscimage.REFNumber, ref_text, 0, 1);

  printf(
    "%s\t%lld\t%d\t%lf\t%lf\t%f\t%f\t%.3lf\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\n",
    ref_text,
    pGscImage->gscimage.REFNumber,
    pGscImage->gsc_bin_index,
    pGscImage->gscimage.ra,
    pGscImage->gscimage.dec,
    dra,
    ddec,
    epoch_delta + GSC_EQUINOX,
    pGscImage->gscimage.RaPM,
    pGscImage->gscimage.DecPM,
    RaSigmaPM,
    DecSigmaPM,
    pGscImage->gscimage.Stdmag,
    pGscImage->gscimage.color,
    pGscImage->gscimage.VFlag,
    pGscImage->gscimage.MAGFlag,
    pGscImage->gscimage.class
  );
}


static void
SearchCatalogDecBin(
  GSCIMAGEBIN *pGscImageTable,
  int gscImageEntries,
  CACHEROW **inputCache,
  int numCacheLines,
  int decBin,
  double limitingMagnitude,
  double epoch_delta
) {
  int cacheIndex;
  int imageIndex;
  CACHEROW *pCacheRow;
  INPUTDATA *pInput;
  double minImageRa;
  double maxImageRa;
  double minImageDec;
  double maxImageDec;
  double minCatalogDec = (decBin * pGscBin->bin_size) - 90.0 - PROPER_MOTION_DEGREES;
  double maxCatalogDec = ((decBin + 1) * pGscBin->bin_size) - 90.0 + PROPER_MOTION_DEGREES;
  double minCatalogRa;
  double maxCatalogRa;
  GSCIMAGEBIN *pGscImage;

  if (gscImageEntries == 0) {
    /* Nothing to do here */
    return;
  }

  if (PROPER_MOTION_DEGREES != 0.0) {
    int i;
    minCatalogDec = 90.0;
    maxCatalogDec = -90.0;

    for (i = 0; i < gscImageEntries; i++) {
      pGscImage = &pGscImageTable[i];

      if (pGscImage->gscimage.dec > maxCatalogDec) {
        maxCatalogDec = pGscImage->gscimage.dec + pGscImage->gscimage.DecSigmaPM;
      }

      if (pGscImage->gscimage.dec < minCatalogDec) {
        minCatalogDec = pGscImage->gscimage.dec - pGscImage->gscimage.DecSigmaPM;
      }
    }
  }

  pGscImage = &pGscImageTable[0];
  minCatalogRa = pGscImage->gscimage.ra - pGscImage->gscimage.RaSigmaPM;
  pGscImage = &pGscImageTable[gscImageEntries - 1];
  maxCatalogRa = pGscImage->gscimage.ra + pGscImage->gscimage.DecSigmaPM;
  pGscImage = NULL; /* For safety */

  for (cacheIndex = 0; cacheIndex < numCacheLines; cacheIndex++) {
    pCacheRow = inputCache[cacheIndex];

    if (pCacheRow != NULL && pCacheRow->numEntries > 0) {
      pCacheRow->curIndex = 0;

      for (imageIndex = 0; imageIndex < pCacheRow->numEntries; imageIndex++) {
        int numSearched = 1;
        int curIndex = pCacheRow->curIndex;
        int curIncrementFlag;
        int originalCase;
        int insideFlag = 0; /* nonzero if we are inside the match area */
        int wrapFlag; /* Nonzero if we just wrapped around the first point of Ares */

        pInput = &pCacheRow->inputdata[imageIndex];

        maxImageDec = pInput->dec + search_radius_deg;
        minImageDec = pInput->dec - search_radius_deg;

        if (minImageDec > maxCatalogDec || maxImageDec < minCatalogDec) {
          /* Image not in this catalog slice */
          continue;
        }

        if (pInput->dec != 90.0) {
          double searchRadiusRa = search_radius_deg / cos(pInput->dec * DEGREES_TO_RAD);
          maxImageRa = pInput->ra + searchRadiusRa;
          minImageRa = pInput->ra - searchRadiusRa;

          if (maxImageRa < 360.0 && minImageRa > 0.0) {
            if (maxImageRa < minCatalogRa || minImageRa > maxCatalogRa) {
              continue;
            }
          }
        } else {
          minImageRa = 0.0;
          maxImageRa = 359.99999;
        }

        // The GSC slice is sorted in increasing RA.

        pGscImage = &pGscImageTable[curIndex];

        switch (CheckRA(pGscImage, minImageRa, maxImageRa)) {
        case -1:
          /* Search forward */
          originalCase = -1;
          curIncrementFlag = +1;
          break;

        case 0:
          // Search forward until we are above the region of interest. Since we
          // have to search again in the opposite direction, do not set the
          // insideFlag.
          originalCase = 0;
          curIncrementFlag = +1;

          if (
            (pGscImage->gscimage.dec <= (maxImageDec + pGscImage->gscimage.DecSigmaPM)) &&
            (pGscImage->gscimage.dec >= (minImageDec - pGscImage->gscimage.DecSigmaPM))
          ) {
            output_match(
              pGscImage,
              pInput,
              limitingMagnitude,
              epoch_delta
            );
          }
          break;

        case 1:
          // Search backward
          originalCase = 1;
          curIncrementFlag = -1;
          break;
        }

        while (numSearched < gscImageEntries) {
          numSearched++;
          curIndex += curIncrementFlag;
          wrapFlag = 0;

          if (curIndex >= gscImageEntries) {
            curIndex = 0;
            wrapFlag = 1;
          } else if (curIndex < 0) {
            curIndex = gscImageEntries -1;
            wrapFlag = 1;
          }

          pGscImage = &pGscImageTable[curIndex];

          switch (CheckRA(pGscImage, minImageRa, maxImageRa)) {
          case -1:
            if (insideFlag == 1 || (wrapFlag == 0 && originalCase == 1)) {
              /* We found the region searching backward and are now
                 below the region.  Terminate the search */
              numSearched = gscImageEntries;
            } else if (originalCase == 0) {
              /* We were searching forward, and are
                 now outside the region.  Go back to the beginning and
                 reverse direction */
              curIndex = pCacheRow->curIndex;
              curIncrementFlag = -1;
              insideFlag = 1; /* Do this only once */
            } else {
              /* Advance the saved index until we find something */
              pCacheRow->curIndex = curIndex;
            }
            break;

          case 0:
            /* Keep going in the same direction until we are outside the region */
            if (originalCase != 0) {
              insideFlag = 1;
            }

            if (
              (pGscImage->gscimage.dec <= (maxImageDec + pGscImage->gscimage.DecSigmaPM)) &&
              (pGscImage->gscimage.dec >= (minImageDec - pGscImage->gscimage.DecSigmaPM))
            ) {
              output_match(
                pGscImage,
                pInput,
                limitingMagnitude,
                epoch_delta
              );
            }
            break;

          case 1:
            if (insideFlag == 1 || (wrapFlag == 0 && originalCase == -1)) {
              // We searched through the region of interest and are done
              numSearched = gscImageEntries;
            } else if (originalCase == 0) {
              // We were searching forward, and are now above the region. Go
              // back to the beginning and reverse direction
              curIndex = pCacheRow->curIndex;
              curIncrementFlag = -1;
              insideFlag = 1; /* Do this only once */
            }
            break;
          }
        } /* Loop over all catalog entries */
      } /* Loop over all images in the cache */
    } /* Got a good cache line */
  } /* loop over all cache lines */
}


int
main(int argc, char *argv[])
{
  int errorFlag = 0;
  double limitingMagnitude = 0.0;
  File index_handle;
  File catalog_handle;
  STARINDEX *pStarIndex = NULL;
  int starIndexAlloc = 0;
  int starIndexEntries;
  GSCIMAGEBIN *pGscImage = NULL;
  GSCIMAGE *pGscStorage = NULL;
  int gscImageAlloc = 0;
  int gscImageEntries;
  char *dotPtr;

  // This default value is J2000.0 minus half a day (!); wcstools' jd2ep() of
  // this yields 2000.000000. I would have thought that 2451545.0 would have
  // done so, but no -- but, it does for wcstools' jd2ep*j*(). I don't fully
  // understand what's going on here but there's no way that 12 hours of proper
  // motion is affecting anything in the DASCH analysis.
  double geocentricJD = 2451544.5;
  double userepoch;
  double epoch_delta;

  double ra_deg = -99.0;
  double dec_deg = -99.0;
  char catalog_path[MAX_BUFFER];
  char index_path[MAX_BUFFER];
  CACHEROW input;
  CACHEROW *input_list[1];

  catalog_path[0] = '\0';

  for (argv++; --argc > 0; argv++) {
    char *argstr = *argv;

    if (strcmp(argstr, "--help") == 0) {
      errorFlag = 1;
      break;
    } else if (argstr[0] != '-') {
      errorFlag = 1;
      fprintf(stderr, "error: unhandled argument \"%s\"\n", argstr);
    } else {
      char cmdchar;

      while ((cmdchar = *++argstr) != '\0') {
        switch(cmdchar) {

        case '?':
        case 'H':
        case 'h':
          errorFlag = 1;
          break;

        // -j <geocentric JD>
        case 'j':
          argc--;
          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else if (sscanf(*++argv, "%lf", &geocentricJD) != 1) {
            fprintf(stderr, "error: Unable to decode geocentric Julian Day %s\n", *argv);
            errorFlag = 1;
          }
          break;

        // -c <catalog path>
        case 'c':
          argc--;
          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            strncpy(catalog_path, *++argv, MAX_BUFFER - 2);
            if (strlen(*argv) >= MAX_BUFFER - 2) {
              fprintf(stderr, "error: buffer overflow in catalog_path\n");
            }
          }
          break;

        // -e <RA(deg)> <dec(deg)>
        case 'e':
          argc--;
          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else if (sscanf(*++argv, "%lf %lf", &ra_deg, &dec_deg) != 2) {
            fprintf(stderr, "error: Unable to decode RA/dec in decimal degrees: \"%s\"\n", *argv);
            errorFlag = 1;
          }
          break;

        // -m <limiting mag>
        case 'm':
          argc--;
          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else if (sscanf(*++argv, "%lf", &limitingMagnitude) != 1) {
            fprintf(stderr, "error: unable to decode limiting magnitude argument \"%s\"\n", *argv);
            errorFlag = 1;
          }
          break;

        // -r <search radius (arcsec)>
        case 'r':
          argc--;

          if (argc < 1) {
            fprintf(stderr, "error: insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            double arcsec;

            if (sscanf(*++argv, "%lf", &arcsec) != 1) {
              fprintf(stderr, "error: unable to decode search radius argument \"%s\"\n", *argv);
              errorFlag = 1;
            } else if (arcsec <= 0.) {
              fprintf(stderr, "error: illegal search radius value \"%s\"\n", *argv);
              errorFlag = 1;
            } else {
              search_radius_deg = arcsec / 3600;
            }
          }
          break;

        default:
          fprintf(stderr, "error: illegal command -%c-", cmdchar);
          errorFlag = 1;
          break;
        }
      }
    }
  }

  if (catalog_path[0] == '\0') {
    fprintf(stderr, "error: the `-c <catpath>` argument is required\n");
    errorFlag = 1;
  }

  if (dec_deg < -90.0) {
    fprintf(stderr, "error: the `-e <RA(deg) dec(deg)>` argument is required\n");
    errorFlag = 1;
  }

  if (errorFlag) {
    fprintf(
      stderr,
      "usage: dasch-querycat [options...]\n\n"
      "  -c <catpath>             Path to the input catalog `.dat` file (required)\n"
      "  -e <RA(deg) dec(deg)>    Equatorial coordinates of search position as decimal degrees (required)\n"
      "  -r <radius(arcsec)>      Search radius (default: 20 arcsec)\n"
      "  -j <geocentricJD>        Geocentric Julian Date (timescale\?\?) epoch of position for results (default: J2000.0)\n"
      "  -m <limiting-mag>        Limiting magnitude filter to apply to catalog (default: no filtering)\n"
      "  -h, -H, -?, --help       Show this help output\n"
    );
    return 1;
  }

  // Further setup

  strcpy(index_path, catalog_path);
  dotPtr = strrchr(index_path, '.');
  if (dotPtr) {
    *dotPtr = '\0';
  }
  strcat(index_path, ".idx");

  catalog_handle = Open(catalog_path, "r");
  if (catalog_handle == NULL) {
    fprintf(stderr, "error: Could not open catalog file %s\n", catalog_path);
    return 1;
  }

  index_handle = Open(index_path, "r");
  if (index_handle == NULL) {
    fprintf(stderr, "error: Could not open catalog index file %s\n", index_path);
    return 1;
  }

  userepoch = jd2ep(geocentricJD);
  epoch_delta = userepoch - GSC_EQUINOX;

  input_list[0] = &input;

  input.numEntries = 1;
  input.curIndex = 0;
  input.inputdata = malloc(sizeof(INPUTDATA));
  input.inputdata->ra = ra_deg;
  input.inputdata->dec = dec_deg;

  // Do it!

  printf(
    "ref_text\t"
    "ref_number\t"
    "gscBinIndex\t"
    "raDeg\t"
    "decDeg\t"
    "draAsec\t"
    "ddecAsec\t"
    "posEpoch\t"
    "pmRaMasyr\t"
    "pmDecMasyr\t"
    "uPMRaMasyr\t"
    "uPMDecMasyr\t"
    "stdmag\t"
    "color\t"
    "vFlag\t"
    "magFlag\t"
    "class\n"
  );
  printf(
    "--------\t"
    "----------\t"
    "-----------\t"
    "-----\t"
    "------\t"
    "-------\t"
    "--------\t"
    "--------\t"
    "---------\t"
    "----------\t"
    "----------\t"
    "-----------\t"
    "------\t"
    "-----\t"
    "-----\t"
    "-------\t"
    "-----\n"
  );

  {
    double minRa[2];
    double maxRa[2];
    int n_ra_chunks = 1;
    double minDec = fmax(dec_deg - search_radius_deg, -90.);
    double maxDec = fmin(dec_deg + search_radius_deg, 89.999);  // getDecBin will error out if we go to exactly 90
    int bin0 = GetDecBin(pGscBin, minDec, "dasch-querycat");
    int bin1 = GetDecBin(pGscBin, maxDec, "dasch-querycat");
    double cos_dec = fmin(cos(minDec * DEGREES_TO_RAD), cos(maxDec * DEGREES_TO_RAD));
    int ibin;

    if (cos_dec <= 0) {
      minRa[0] = 0.;
      maxRa[0] = 360.;
    } else {
      double search_radius_ra = search_radius_deg / cos_dec;

      minRa[0] = ra_deg - search_radius_ra;
      maxRa[0] = ra_deg + search_radius_ra;

      if (minRa[0] < 0 && maxRa[0] > 360) {
        // We cover all RA's, which might happen with a reasonable radius if
        // we're right at the poles. This is OK.
        minRa[0] = 0.;
        maxRa[0] = 360.;
      } else if (minRa[0] < 0) {
        // We need to break our search into two RA chunks.
        // {min, max}[0] = (0, naive-max)
        // {min, max}[1] = (wrapped-naive-min, 360)
        minRa[1] = minRa[0] + 360;
        minRa[0] = 0.;
        maxRa[1] = 360.;
        n_ra_chunks = 2;
      } else if (maxRa[0] > 360) {
        // Analogous to the previous case.
        maxRa[1] = maxRa[0] - 360;
        maxRa[0] = 360;
        minRa[1] = 0.;
        n_ra_chunks = 2;
      }
    }

    for (ibin = bin0; ibin <= bin1; ibin++) {
      if (
        ReadCatalogDecBin(
          catalog_handle,
          index_handle,
          &pStarIndex,
          &starIndexAlloc,
          &starIndexEntries,
          &pGscImage,
          &pGscStorage,
          &gscImageAlloc,
          &gscImageEntries,
          ibin,
          minRa,
          maxRa,
          limitingMagnitude,
          epoch_delta,
          n_ra_chunks
        )
      ) {
        SearchCatalogDecBin(
          pGscImage,
          gscImageEntries,
          (CACHEROW **) &input_list,
          1, // numCacheLines
          ibin,
          limitingMagnitude,
          epoch_delta
        );
      }
    }
  }

  // All done. Clean up.

  if (pStarIndex != NULL) {
    free(pStarIndex);
  }

  if (pGscImage != NULL) {
    free(pGscImage);
  }

  if (pGscStorage != NULL) {
    free(pGscStorage);
  }

  Close(catalog_handle);
  Close(index_handle);

  return 0;
}
