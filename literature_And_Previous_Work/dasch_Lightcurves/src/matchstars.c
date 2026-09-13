// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// "Rewrite of John Roll's search utility using a binary GSC2.3.2 catalog."

#include <math.h>
#include <time.h>
#include <errno.h>

#include "table.h"

#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "libwcs/fitshead.h"

#include "pipelineutils.h"


#define MAX_BUFFER 512


/* We need to add the gsc_bin_index to the GSCIMAGE structure */
typedef struct _gscimagebin {
  GSCIMAGE gscimage;  /* The normal GSCIMAGE structure */
  int gsc_bin_index;  /* The gsc_bin_index for this structure */
} GSCIMAGEBIN;

/* reformatted Sextractor results */
typedef struct _inputdata {
  double MAG_ISO;     /* Sextractor isophotonic magnitude */
  double X_IMAGE;        /* Sextractor X location in pixels */
  double Y_IMAGE;        /* Sextractor Y location in pixels */
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
  double FWHM_IMAGE;  /* Image full width half maximum in pixels */
  double FWHM_WORLD;  /* Image full width half maximum in degrees */
  double THRESHOLD;   /* Sextractor threshold */
  double FLUX_MAX;    /* Sextractor maximum flux */
  int NUMBER;         /* Sextractor reference number */
  int AFLAGS;          /* Sextractor flags */
  int BFLAGS;          /* Sextractor flags */
  double aLength;      /* Long axis in pixels */
  double bLength;      /* Short axis in pixels */
  double THETA_J2000;  /* Image orientation in degrees */
  int spatial_bin;     /* Spatial bin */
} INPUTDATA;

typedef struct _cacherow {
  int lineNumber;  /* Line Number */
  int numEntries;  /* Number of entries in this row */
  int allocEntries; /* Entries allocated for this row */
  int curIndex;    /* Current index into this row */
  int decBin;      /* Declination bin */
  double minDec;   /* minimum declination */
  double maxDec;   /* maximum declination */
  int overlapIndex;  /* Index for equinox overlap */
  INPUTDATA *inputdata; /*  row contents */
} CACHEROW;


extern GSCBIN gscBin64;

static double properMotionDegrees = PROPER_MOTION_DEGREES;
static double maxpmerror = 0;
static long long maxpmerrorREFNumber = 0;
static PGSCBIN pGscBin = &gscBin64;
static int zeroMagnitudeCount = 0;

static void WriteMagnitudeTable(int *magnitudeTable, int *MAG_ISO_Table, FILE *estimateHandle, char *estimatename);


/*
 * Return 0 if the gsc RA is within range,
 *        1 if the gsc RA is higher
 *       -1 if the gsc RA is lower
 */
static inline int
CheckRA(
  GSCIMAGEBIN *pGscImage,
  double minImageRa,
  double maxImageRa,
  char* outputname
) {
  double gscRa = pGscImage->gscimage.ra;

  minImageRa -= pGscImage->gscimage.RaSigmaPM;
  maxImageRa += pGscImage->gscimage.RaSigmaPM;

  if (minImageRa < 0.0) {
    /* Here the low end of the range croses the first point of Aries */
    if (gscRa >= minImageRa + 360.0 || gscRa <= maxImageRa) {
      return 0;
    }

    if ((gscRa - maxImageRa) < (minImageRa + 360.0 - gscRa)) {
      return 1;
    }

    return -1;
  }

  if (maxImageRa >= 360.0) {
    /* Here the upper end of the range crosses the first point of Aries */
    if (gscRa >= minImageRa || gscRa <= (maxImageRa - 360.0)) {
      return 0;
    }

    if ((minImageRa - gscRa) < (gscRa - (maxImageRa - 360.0))) {
      return -1;
    }

    return 1;
  }

  if (gscRa < minImageRa) {
    return -1;
  } else if (gscRa > maxImageRa) {
    return 1;
  } else {
    return 0;
  }
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
  } else {
    return 0;
  }
}


/* Sort routine based on the right ascension */
static int
ImageCompare(const void *first, const void *second)
{
  double numberFirst = ((INPUTDATA *) first)->ra;
  double numberSecond = ((INPUTDATA *) second)->ra;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  } else {
    return 0;
  }
}


static int
ReadCatalogDecBin(
  File catalogHandle,
  File indexHandle,
  PSTARINDEX *pStarIndex,
  int *pStarIndexAlloc,
  int *pStarIndexEntries,
  GSCIMAGEBIN **pGscImage,
  PGSCIMAGE *pGscStorage,
  int *pGscImageAlloc,
  int *pGscImageEntries,
  int curCatalogBin,
  double *minRa,
  double *maxRa,
  double limitingMagnitude,
  char *outputname,
  double plateepoch,
  int *properMotionErrLimCount,
  int *duplicateCount,
  int catalogNumber,
  signed long long *pCurMemory,
  signed long long *pMaxMemory
) {
  PSTARINDEX pCurStarIndex = *pStarIndex;
  PSTARINDEX tempStarIndex;

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
  PGSCIMAGE pCurGscStorage = *pGscStorage;
  PGSCIMAGE tempGscStorage = NULL;
  GSCIMAGEBIN *tempGscImage;
  GSCIMAGEBIN *tempGscImage2;
  double factor;
  double raMotion;
  double tmpRa;
  int motionFlag = 0;

  // Compute RA motion avoiding explosions around the poles. The cutoff value
  // of 0.01 corresponds to decs of around +/- 89.5 degrees.
  factor = cos(DEGREES_TO_RAD * declination);
  factor = fmax(factor, 0.01);
  raMotion = properMotionDegrees / factor;

  *pStarIndexEntries = 0;
  *pGscImageEntries = 0;

  for (raIndex = 0; raIndex < 2; raIndex++) {
    if (minRa[raIndex] >= maxRa[raIndex]) {
      /* No work here */
      continue;
    }

    /* Find out which star index bins we need to read */
    tmpRa = minRa[raIndex] - raMotion - pGscBin->bin_size;
    if (tmpRa < 0.0) {
      tmpRa = 0.0;
    }

    firstBin = GetGSCBin(pGscBin, tmpRa, declination, &decBin, &raBin, outputname);
    GetSubbins(pGscBin, firstBin, &raBin, &decBin, outputname);

    tmpRa = maxRa[raIndex] + raMotion + pGscBin->bin_size;
    if (tmpRa >= 360.0) {
      tmpRa = 359.9999;
    }

    lastBin = GetGSCBin(pGscBin, tmpRa, declination, &decBin, &raBin, outputname);
    numBins = lastBin - firstBin + 1;
    if (numBins <= 0) {
      printf(
        "ERROR: numBins is %d in ReadCatalogDecBin for %s (%d, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf)\n",
        numBins,
        outputname,
        raIndex,
        minRa[raIndex],
        maxRa[raIndex],
        raMotion,
        pGscBin->bin_size,
        declination
      );
      exit(1);
    }

    if ((numBins + totalBins) >= *pStarIndexAlloc) {
      CalcMemory(1, (numBins + totalBins + 100 - *pStarIndexAlloc) * sizeof(STARINDEX) ,pCurMemory, pMaxMemory);

      *pStarIndexAlloc = numBins + totalBins + 100;
      tempStarIndex = realloc(*pStarIndex, (*pStarIndexAlloc) * sizeof(STARINDEX));
      if (tempStarIndex == NULL) {
        printf("ERROR allocating star index of size %d for %s\n",*pStarIndexAlloc,outputname);
        exit(1);
      }

      *pStarIndex = tempStarIndex;
      pCurStarIndex = &tempStarIndex[totalBins];
    }

    totalBins += numBins;
    Seek(indexHandle,firstBin * sizeof(STARINDEX),SEEK_SET);
    readItems = Read(indexHandle,pCurStarIndex,sizeof(STARINDEX),numBins);
    if (readItems != numBins) {
      printf("ERROR reading star index file for %s\n",outputname);
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
      CalcMemory(1,(totalStars + numStars + 100 - *pGscImageAlloc)*(sizeof(STARINDEX)+sizeof(int)),pCurMemory,pMaxMemory);
      *pGscImageAlloc = totalStars + numStars + 100;
      tempGscImage = realloc(*pGscImage,(*pGscImageAlloc) * sizeof(GSCIMAGEBIN));
      if (tempGscImage == NULL) {
        printf("ERROR allocating gsc index of size %d for %s\n",*pGscImageAlloc,outputname);
        exit(1);
      }

      *pGscImage = tempGscImage;
      pCurGscImage = &tempGscImage[totalStars];

      tempGscStorage = realloc(*pGscStorage,(*pGscImageAlloc)*sizeof(GSCIMAGE));
      if (tempGscStorage == NULL) {
        printf("ERROR allocating pCurGscStorage of size %d for %s\n",*pStarIndexAlloc,outputname);
        exit(1);
      }

      *pGscStorage = tempGscStorage;
      pCurGscStorage = tempGscStorage;
    }

    /* Go to the first bin in the index */

    Seek(catalogHandle,pCurStarIndex->offset,SEEK_SET);
    totalStars += numStars;
    /* Now read the required number of stars */
    readItems = Read(catalogHandle,pCurGscStorage,sizeof(GSCIMAGE),numStars);
    if (readItems != numStars) {
      printf("ERROR reading GSC catalog file for %s\n",outputname);
      exit(1);
    }

    /* Now copy everything into the larger GSCIMAGEBIN area */
    for (starIndex = 0; starIndex < numStars; starIndex++) {
      tempGscImage = &pCurGscImage[starIndex];
      tempGscStorage = &pCurGscStorage[starIndex];
      memcpy(&tempGscImage->gscimage,tempGscStorage,sizeof(GSCIMAGE));
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
          memcpy(tempGscImage,tempGscImage2,sizeof(GSCIMAGEBIN));
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
        memcpy(tempGscImage,tempGscImage2,sizeof(GSCIMAGEBIN));
      }
      totalStars++;
    } else {
      zeroMagnitudeCount++;
    }
  }

  /* Now save the gsc bin indices in the gsc bin table prior to changing
     the ra and dec because of proper motion */
  pCurGscImage = *pGscImage;

  for (starIndex = 0; starIndex < totalStars; starIndex++) {
    tempGscImage = &pCurGscImage[starIndex];

    if (tempGscImage->gscimage.dec > 90. || tempGscImage->gscimage.dec < -90.) {
      fprintf(stderr, "about to crash: TNX.db source has illegal decl.; suggests corrupt input file\n");
      fflush(stderr);
    }

    tempGscImage->gsc_bin_index = GetGSCBin(
      pGscBin,
      tempGscImage->gscimage.ra,
      tempGscImage->gscimage.dec,
      &decBin,
      &raBin,
      outputname
    );
  }

  // Move any stars that need it

  pCurGscImage = *pGscImage;

  for (starIndex = 0; starIndex < totalStars; starIndex++) {
    tempGscImage2 = &pCurGscImage[starIndex];

    if (tempGscImage2->gscimage.RaPM != 0 || tempGscImage2->gscimage.DecPM != 0) {
      motionFlag = 1;
      tempGscImage2->gscimage.dec += (tempGscImage2->gscimage.DecPM * plateepoch) / (3600.0 * 1000.0);

      if (tempGscImage2->gscimage.dec < 90.0 && tempGscImage2->gscimage.dec > -90.0) {
        tempGscImage2->gscimage.ra +=
          ((tempGscImage2->gscimage.RaPM * plateepoch) / (3600.0 * 1000.0))
          / cos(tempGscImage2->gscimage.dec * DEGREES_TO_RAD);
      }
    }

    if (tempGscImage2->gscimage.RaSigmaPM != 0.0 || tempGscImage2->gscimage.DecSigmaPM != 0.0) {
      double factor2 = cos(DEGREES_TO_RAD * tempGscImage2->gscimage.dec);
      factor2 = fmax(factor2, 0.0001);

      // Overload the meaning of these fields; convert to degrees
      tempGscImage2->gscimage.RaSigmaPM *= fabs((PROPER_MOTION_SIGMA * plateepoch) / (3600.0 * 1000.0 * factor2));
      tempGscImage2->gscimage.DecSigmaPM *= fabs((PROPER_MOTION_SIGMA * plateepoch) / (3600.0 * 1000.0));
    }
  }

  if (motionFlag) {
    /* We need to resort this bin so that the RA's are back in order */
    qsort((void*) pCurGscImage, totalStars, sizeof(GSCIMAGEBIN), GscImageCompare);
  }

  if (minRa[1] < maxRa[1]) {
    /* There may be duplicates here that we need to remove */
    if (!motionFlag) {
      qsort((void*) pCurGscImage, totalStars, sizeof(GSCIMAGEBIN), GscImageCompare);
    }

    starIndex = 1;
    while (starIndex < totalStars) {
      tempGscImage = &pCurGscImage[starIndex - 1];
      tempGscImage2 = &pCurGscImage[starIndex];

      if (tempGscImage->gscimage.REFNumber == tempGscImage2->gscimage.REFNumber) {
        /* Duplicate entry - eliminate it, but keep sorting order */
        (*duplicateCount)++;
        for (starIndex2 = starIndex; starIndex2 < (totalStars-1); starIndex2++) {
          tempGscImage = &pCurGscImage[starIndex2];
          tempGscImage2 = &pCurGscImage[starIndex2+1];
          memcpy(tempGscImage,tempGscImage2,sizeof(GSCIMAGEBIN));
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
OutputImage(
  FILE *outputHandle,
  GSCIMAGEBIN *pGscImage,
  INPUTDATA *pInput,
  double limitingMagnitude,
  double plateepoch,
  int *properMotionErrLimCount,
  int *pOutputCount,
  int *magnitudeTable,
  int catalogNumber
) {
  double ddec;
  double dra;
  double drad;
  double factor;
  double pm;
  double RaSigmaPM;
  double DecSigmaPM;
  double pmerror;
  int BFLAGS; /* Make a copy because otherwise, the next match may have the wrong BFLAGS */
  char REF[MAX_REF];
  double deltara;

  if (limitingMagnitude != 0.0 && pGscImage->gscimage.Stdmag > limitingMagnitude) {
    return;
  }

  (*pOutputCount)++;

  factor = cos(DEGREES_TO_RAD * ((pInput->dec + pGscImage->gscimage.dec) / 2.0));
  deltara = pInput->ra - pGscImage->gscimage.ra;

  if (deltara < -180.0) {
    deltara = 360. + deltara;
  } else if (deltara > 180.0) {
    deltara = deltara - 360.0;
  }

  dra = 3600 * factor * deltara;
  ddec = 3600 * (pInput->dec - pGscImage->gscimage.dec);
  drad = sqrt((dra * dra) + (ddec * ddec));

  BFLAGS = pInput->BFLAGS;
  pm = sqrt(sqr(pGscImage->gscimage.RaPM) + sqr(pGscImage->gscimage.DecPM));
  if (pm > 0.0) {
    BFLAGS |= (1 << FILTER_BFLAG_PROPER_MOTION);
  }

  // before we write things out, restore RaSigmaPM and DecSigmaPM to their original values (mas/yr)
  RaSigmaPM = fabs(((3600.0 * 1000.0 * factor) / (PROPER_MOTION_SIGMA * plateepoch)) * pGscImage->gscimage.RaSigmaPM);
  DecSigmaPM = fabs(((3600.0 * 1000.0) / (PROPER_MOTION_SIGMA * plateepoch)) * pGscImage->gscimage.DecSigmaPM);

  pmerror = fabs(PROPER_MOTION_SIGMA * plateepoch * (sqrt(sqr(RaSigmaPM) * sqr(DecSigmaPM))) / 1000.0);
  if (pmerror > maxpmerror) {
    maxpmerror = pmerror;
    maxpmerrorREFNumber = pGscImage->gscimage.REFNumber;
  }

  if (pmerror > PROPER_MOTION_ERROR_LIMIT) {
    BFLAGS |= (1 << FILTER_BFLAG_PMERROR);
    (*properMotionErrLimCount)++;
  }

  GetREF(pGscImage->gscimage.REFNumber, REF, 0, 1);

  fprintf(
    outputHandle,
    "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n",
    pInput->NUMBER,
    pInput->MAG_ISO,
    pInput->X_IMAGE,
    pInput->Y_IMAGE,
    pInput->ra,
    pInput->dec,
    pInput->FWHM_IMAGE,
    pInput->FWHM_WORLD,
    pInput->AFLAGS,
    BFLAGS,
    REF,
    pGscImage->gscimage.ra,
    pGscImage->gscimage.dec,
    RaSigmaPM,
    DecSigmaPM,
    pGscImage->gscimage.Stdmag,
    pGscImage->gscimage.color,
    pGscImage->gscimage.VFlag,
    pGscImage->gscimage.MAGFlag,
    pGscImage->gscimage.class,
    dra,
    ddec,
    drad,
    pInput->THRESHOLD,
    pInput->FLUX_MAX,
    pInput->aLength,
    pInput->bLength,
    pInput->THETA_J2000,
    pGscImage->gsc_bin_index,
    pInput->spatial_bin
  );

  if (magnitudeTable != NULL && pGscImage->gscimage.Stdmag < 90.0 && pInput->spatial_bin > 0) {
    int iStdmag = (pGscImage->gscimage.Stdmag - MIN_STDMAG) * STDMAG_BINS / (MAX_STDMAG - MIN_STDMAG);
    int iMAG_ISO = (pInput->MAG_ISO - MIN_MAG_ISO) * MAG_ISO_BINS / (MAX_MAG_ISO - MIN_MAG_ISO);

    if (iStdmag < 0) {
      iStdmag = 0;
    } else if (iStdmag >= STDMAG_BINS) {
      iStdmag = STDMAG_BINS - 1;
    }

    if (iMAG_ISO < 0) {
      iMAG_ISO = 0;
    } else if (iMAG_ISO >= MAG_ISO_BINS) {
      iMAG_ISO = MAG_ISO_BINS - 1;
    }

    magnitudeTable[iMAG_ISO + (MAG_ISO_BINS * iStdmag) + ((pInput->spatial_bin - 1) * MAGNITUDE_BINS)]++;
  }
}


static void
SearchCatalogDecBin(
  FILE *outputHandle,
  GSCIMAGEBIN *pGscImageTable,
  int gscImageEntries,
  CACHEROW **inputCache,
  int numCacheLines,
  double scale,
  int decBin,
  double limitingMagnitude,
  double plateepoch,
  int *pOutputCount,
  int *properMotionErrLimCount,
  char *outputname,
  int *magnitudeTable,
  int catalogNumber
) {
  int cacheIndex;
  int imageIndex;
  CACHEROW *pCacheRow;
  INPUTDATA *pInput;
  double minImageRa;
  double maxImageRa;
  double minImageDec;
  double maxImageDec;
  double searchRadiusDec;
  double searchRadiusRa;
  double minCatalogDec = (decBin * pGscBin->bin_size) - 90.0 - properMotionDegrees;
  double maxCatalogDec = ((decBin + 1) * pGscBin->bin_size) - 90.0 + properMotionDegrees;
  double minCatalogRa;
  double maxCatalogRa;
  GSCIMAGEBIN *pGscImage;
  char series[MAX_SERIES_STRING];
  int plateNumber;

  if (gscImageEntries == 0) {
    return;
  }

  if (properMotionDegrees != 0.0) {
    int curIndex;

    minCatalogDec = 90.0;
    maxCatalogDec = -90.0;

    for (curIndex = 0; curIndex < gscImageEntries; curIndex++) {
      pGscImage = &pGscImageTable[curIndex];

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

  if (ParseFilename2(outputname, series, &plateNumber) == 0) {
    printf("ERROR: SearchCatalogDecBin failed to parse filename %s\n", outputname);
    exit(1);
  }

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

        /* A change here requires a change in filterblended.c */
        /* THETA_J2000 is measured in degrees clockwise from North. East is 90 degrees, and West is 270 degrees */
        double cosTheta = cos(pInput->THETA_J2000 * DEGREES_TO_RAD);
        double sinTheta = sin(pInput->THETA_J2000 * DEGREES_TO_RAD);

        if (cosTheta < 0.0) {
          cosTheta = -cosTheta;
        }

        if (sinTheta < 0.0) {
          sinTheta = -sinTheta;
        }

        double raLength = (pInput->aLength * sinTheta) + (pInput->bLength * cosTheta); /* ra extent in pixels */
        double decLength = (pInput->aLength * cosTheta) + (pInput->bLength * sinTheta); /* dec extent in pixels */

        searchRadiusDec = (decLength + MaxDradPixels(series, plateNumber)) * (scale / 3600.0); /* Search radius in degrees */
        searchRadiusRa  = (raLength + MaxDradPixels(series, plateNumber)) * (scale / 3600.0); /* Search radius in degrees */

        maxImageDec = pInput->dec + searchRadiusDec;
        minImageDec = pInput->dec - searchRadiusDec;

        if (minImageDec > maxCatalogDec || maxImageDec < minCatalogDec) {
          continue; /* Image not in this catalog slice */
        }

        if (pInput->dec != 90.0) {
          searchRadiusRa = searchRadiusRa / cos(pInput->dec * DEGREES_TO_RAD);
          maxImageRa = pInput->ra + searchRadiusRa;
          minImageRa = pInput->ra - searchRadiusRa;

          if (pCacheRow->overlapIndex == 0 && maxImageRa < 360.0 && minImageRa > 0.0) {
            if (maxImageRa < minCatalogRa || minImageRa > maxCatalogRa) {
              continue;
            }
          }
        } else {
          minImageRa = 0.0;
          maxImageRa = 359.99999;
        }

        /* The gsc slice is sorted in increasing RA.  */

        pGscImage = &pGscImageTable[curIndex];

        switch (CheckRA(pGscImage, minImageRa, maxImageRa, outputname)) {
        case -1:
          /* Search forward */
          originalCase = -1;
          curIncrementFlag = +1;
          break;

        case 0:
          /* Search forward until we are above the region of interest */
          /* Since we have to search again in the opposite direction, do not set the insideFlag */
          originalCase = 0;
          curIncrementFlag = +1;
          if (
            pGscImage->gscimage.dec <= (maxImageDec + pGscImage->gscimage.DecSigmaPM) &&
            pGscImage->gscimage.dec >= (minImageDec - pGscImage->gscimage.DecSigmaPM)
          ) {
            OutputImage(
              outputHandle,
              pGscImage,
              pInput,
              limitingMagnitude,
              plateepoch,
              properMotionErrLimCount,
              pOutputCount,
              magnitudeTable,
              catalogNumber
            );
          }
          break;

        case 1:
          /* Search Backward */
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
            curIndex = gscImageEntries - 1;
            wrapFlag = 1;
          }

          pGscImage = &pGscImageTable[curIndex];

          switch (CheckRA(pGscImage, minImageRa, maxImageRa, outputname)) {
          case -1:
            if (insideFlag == 1 || (wrapFlag == 0 && originalCase == 1)) {
              // We found the region searching backward and are now below the
              // region. Terminate the search
              numSearched = gscImageEntries;
            } else if (originalCase == 0) {
              // We were searching forward, and are now outside the region. Go
              // back to the beginning and reverse direction.
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
              pGscImage->gscimage.dec <= (maxImageDec + pGscImage->gscimage.DecSigmaPM) &&
              pGscImage->gscimage.dec >= (minImageDec - pGscImage->gscimage.DecSigmaPM)
            ) {
              OutputImage(
                outputHandle,
                pGscImage,
                pInput,
                limitingMagnitude,
                plateepoch,
                properMotionErrLimCount,
                pOutputCount,
                magnitudeTable,
                catalogNumber
              );
            }
            break;

          case 1:
            if (insideFlag == 1 || (wrapFlag == 0 && originalCase == -1)) {
              // We searched through the region of interest and are done
              numSearched = gscImageEntries;
            } else {
              if (originalCase == 0) {
                // We were searching forward, and are now above the region. Go
                // back to the beginning and reverse direction
                curIndex = pCacheRow->curIndex;
                curIncrementFlag = -1;
                insideFlag = 1; /* Do this only once */
              }
            }
            break;
          }
        }
      }
    }
  }
}


static int
ReadInputData(
  INPUTDATA *pInput,
  File input_handle,
  TableHead input_header,
  TblDescriptor input_descriptor,
  TableRow* input_row,
  int mosaicWidth,
  int mosaicHeight,
  int *MAG_ISO_Table
) {
  double edgeDist;
  int iMAG_ISO;

  memset(pInput, 0, sizeof(INPUTDATA));

  *input_row = table_rowget(input_handle, input_header, *input_row, NULL, NULL, 0);
  if (*input_row == NULL) {
    return 0;
  }

  if (!table_loadrow(input_handle, input_header, *input_row, input_descriptor, (char *) pInput)) {
    printf("ERROR: Read Input table_loadrow failed\n");
    return 0;
  }

  pInput->spatial_bin = CalculateBin(mosaicWidth, mosaicHeight, pInput->X_IMAGE, pInput->Y_IMAGE, &edgeDist);

  iMAG_ISO = (pInput->MAG_ISO - MIN_MAG_ISO) * MAG_ISO_BINS / (MAX_MAG_ISO - MIN_MAG_ISO);
  if (iMAG_ISO < 0) {
    iMAG_ISO = 0;
  } else if (iMAG_ISO >= MAG_ISO_BINS) {
    iMAG_ISO = MAG_ISO_BINS - 1;
  }

  if (MAG_ISO_Table != NULL) {
    MAG_ISO_Table[iMAG_ISO + MAG_ISO_BINS * (pInput->spatial_bin - 1)]++;
  }

  return 1;
}


static void
GetRanges(
  CACHEROW **inputCache,
  int numCacheLines,
  double *pMinDec,
  double *pMaxDec,
  double *pMinRa1,
  double *pMaxRa1,
  double *pMinRa2,
  double *pMaxRa2,
  char *outputname
) {
  int cacheIndex;
  CACHEROW *pCacheRow;
  INPUTDATA *pInput;
  double curMinDec;
  double curMaxDec;
  double curMinRa;
  double curMaxRa;

  *pMinDec = 90.0;
  *pMaxDec = -90.0;
  *pMinRa1 = 360.0;
  *pMaxRa1 = 0.0;
  *pMinRa2 = 360.0;
  *pMaxRa2 = 0.0;

  int firstOverlap = 0;
  double minDiff1;
  double minDiff2;

  if (inputCache == NULL) {
    return;
  }

  for (cacheIndex = 0; cacheIndex < numCacheLines; cacheIndex++) {
    pCacheRow = inputCache[cacheIndex];

    if (pCacheRow != NULL && pCacheRow->numEntries > 0) {
      curMinDec = (pCacheRow->decBin * pGscBin->bin_size) - 90.0;
      curMaxDec = curMinDec + pGscBin->bin_size;

      if (curMinDec < *pMinDec) {
        *pMinDec = curMinDec;
      }

      if (curMaxDec > *pMaxDec) {
        *pMaxDec = curMaxDec;
      }

      if (pCacheRow->overlapIndex == 0) {
        /* Normal case, no overlap */
        pInput = pCacheRow->inputdata;
        curMinRa = pInput->ra;
        pInput = &pCacheRow->inputdata[pCacheRow->numEntries - 1];
        curMaxRa = pInput->ra;

        if (firstOverlap == 0) {
          /* Normal case with one range active */
          if (curMinRa < *pMinRa1) {
            *pMinRa1 = curMinRa;
          }
          if (curMaxRa > *pMaxRa1) {
            *pMaxRa1 = curMaxRa;
          }
        } else {
          /* Here we have two ranges - see which range this one is closest to */
          minDiff1 = curMinRa - *pMinRa1;

          if (minDiff1 < 0.0) {
            minDiff1 = -minDiff1;
          }

          minDiff2 = curMinRa - *pMinRa2;

          if (minDiff2 < 0.0) {
            minDiff2 = -minDiff2;
          }

          if (minDiff1 < minDiff2) {
            /* This is the low range */
            if (curMinRa < *pMinRa1) {
              *pMinRa1 = curMinRa;
            }

            if (curMaxRa > *pMaxRa1) {
              *pMaxRa1 = curMaxRa;
            }
          } else {
            /* This is the high range */
            if (curMinRa < *pMinRa2) {
              *pMinRa2 = curMinRa;
            }

            if (curMaxRa > *pMaxRa2) {
              *pMaxRa2 = curMaxRa;
            }
          }
        }
      } else {
        /* Overlap case - two ranges */
        pInput = pCacheRow->inputdata;
        curMinRa = pInput->ra;
        pInput = &pCacheRow->inputdata[pCacheRow->overlapIndex - 1];
        curMaxRa = pInput->ra;

        if (firstOverlap == 0) {
          minDiff1 = curMinRa - *pMinRa1;

          if (minDiff1 < 0.0) {
            minDiff1 = -minDiff1;
          }

          if (minDiff1 > 180.0) {
            /* Our existing range is the high range - move it to the second range */
            *pMinRa2 = *pMinRa1;
            *pMaxRa2 = *pMaxRa1;
            *pMinRa1 = 360.0;
            *pMaxRa1 = 0.0;
          }
        }

        firstOverlap = 1;

        if (curMinRa < *pMinRa1) {
          *pMinRa1 = curMinRa;
        }

        if (curMaxRa > *pMaxRa1) {
          *pMaxRa1 = curMaxRa;
        }

        pInput = &pCacheRow->inputdata[pCacheRow->overlapIndex];
        curMinRa = pInput->ra;
        pInput = &pCacheRow->inputdata[pCacheRow->numEntries - 1];
        curMaxRa = pInput->ra;

        if (curMinRa < *pMinRa2) {
          *pMinRa2 = curMinRa;
        }

        if (curMaxRa > *pMaxRa2) {
          *pMaxRa2 = curMaxRa;
        }
      }
    }
  }
}


static int
ReadNextLine(
  INPUTDATA *pOldInput,
  File input_handle,
  TableHead input_header,
  TblDescriptor input_descriptor,
  TableRow *input_row,
  CACHEROW **inputCache,
  int numCacheLines,
  size_t *input_nrecs,
  int verbose,
  char *outputname,
  signed long long *pCurMemory,
  signed long long *pMaxMemory,
  int mosaicWidth,
  int mosaicHeight,
  int *MAG_ISO_Table
) {
  int cacheIndex;
  CACHEROW *pCacheRow;
  CACHEROW *pCurCacheRow;
  int curCacheIndex = 0;
  int curDecBin;
  int oldDecBin = -1;
  INPUTDATA *pInput;
  INPUTDATA *pPrevInput;
  int index;
  int readResult = 0;

  /* Start by finding a stale cache descriptor */
  pCurCacheRow = NULL;

  for (cacheIndex = 0; cacheIndex < numCacheLines; cacheIndex++) {
    pCacheRow = inputCache[cacheIndex];

    if (pCacheRow == NULL) {
      curCacheIndex = cacheIndex;
      pCurCacheRow = NULL;
      break;
    }

    if (pCurCacheRow == NULL) {
      pCurCacheRow = pCacheRow;
      curCacheIndex = cacheIndex;
    } else if (pCurCacheRow->lineNumber > pCacheRow->lineNumber) {
      pCurCacheRow = pCacheRow;
      curCacheIndex = cacheIndex;
    }
  }

  if (pCurCacheRow == NULL) {
    pCurCacheRow = (CACHEROW *) calloc(1, sizeof(CACHEROW));

    if (pCurCacheRow == NULL) {
      printf("ERROR allocating CACHEROW\n");
      exit(1);
    }

    CalcMemory(1, sizeof(CACHEROW), pCurMemory, pMaxMemory);
    inputCache[curCacheIndex] = pCurCacheRow;
    pCurCacheRow->allocEntries = 100;
    pCurCacheRow->inputdata = (INPUTDATA *) calloc(pCurCacheRow->allocEntries, sizeof(INPUTDATA));
    CalcMemory(1, pCurCacheRow->allocEntries * sizeof(INPUTDATA), pCurMemory, pMaxMemory);
  }

  pCurCacheRow->lineNumber = (*input_nrecs) + 1;
  pCurCacheRow->curIndex = 0;
  pCurCacheRow->numEntries = 0;
  pCurCacheRow->decBin = 0;
  pCurCacheRow->minDec = 0.0;
  pCurCacheRow->maxDec = 0.0;
  pCurCacheRow->overlapIndex = 0;
  pInput = pCurCacheRow->inputdata;

  if (pOldInput->NUMBER != 0) {
    /* We have one star left over from the last invocation.*/
    memcpy(pInput, pOldInput, sizeof(INPUTDATA));
    memset(pOldInput, 0, sizeof(INPUTDATA));
    oldDecBin = GetDecBin(pGscBin, pInput->dec, outputname);

    pCurCacheRow->decBin = oldDecBin;
    pCurCacheRow->minDec = (oldDecBin * pGscBin->bin_size) - 90.0;
    pCurCacheRow->maxDec = ((oldDecBin + 1) * pGscBin->bin_size) - 90.0;
    pCurCacheRow->numEntries++;
    pInput++;
  }

  while (1) {
    if (pCurCacheRow->numEntries >= pCurCacheRow->allocEntries) {
      CalcMemory(1, 100 * sizeof(INPUTDATA), pCurMemory, pMaxMemory);
      pCurCacheRow->allocEntries += 100;

      pInput = realloc(pCurCacheRow->inputdata, pCurCacheRow->allocEntries * sizeof(INPUTDATA));
      if (pInput == NULL) {
        printf("ERROR reallocating INPUTDATA\n");
        exit(1);
      }

      pCurCacheRow->inputdata = pInput;
      pInput = &pCurCacheRow->inputdata[pCurCacheRow->numEntries];
    }

    if (!(
      readResult = ReadInputData(
        pInput,
        input_handle,
        input_header,
        input_descriptor,
        input_row,
        mosaicWidth,
        mosaicHeight,
        MAG_ISO_Table
      )
    )) {
      if (verbose) {
        printf("At end of data in record %zu\n", *input_nrecs);
      }
      break;
    }

    *input_nrecs += 1;
    curDecBin = GetDecBin(pGscBin, pInput->dec, outputname);

    if (curDecBin < oldDecBin) {
      printf("ERROR: Input sextractor file is not sorted by declination for %s\n", outputname);
    } else if (oldDecBin != -1 && curDecBin > oldDecBin) {
      // Next bin. Temporarily save the results for the next invocation.
      memcpy(pOldInput, pInput, sizeof(INPUTDATA));
      break;
    }

    oldDecBin = curDecBin;
    pCurCacheRow->decBin = curDecBin;
    pCurCacheRow->minDec = (curDecBin * pGscBin->bin_size) - 90.0;
    pCurCacheRow->maxDec = ((curDecBin + 1) * pGscBin->bin_size) - 90.0;
    pCurCacheRow->numEntries++;
    pInput++;
  }

  /* We have a full row */
  if (pCurCacheRow->numEntries <= 0) {
    return readResult;
  }

  /* Now sort the row in increasing right ascension */
  qsort((void *) pCurCacheRow->inputdata, pCurCacheRow->numEntries, sizeof(INPUTDATA), ImageCompare);

  /* Now look for an equinox overlap */
  pPrevInput = pCurCacheRow->inputdata;
  for (index = 1; index < pCurCacheRow->numEntries; index++) {
    pInput = &pCurCacheRow->inputdata[index];

    if ((pInput->ra - pPrevInput->ra) > 180.0) {
      pCurCacheRow->overlapIndex = index;
      break;
    }

    pPrevInput = pInput;
  }

  return 1;
}


int
main(int argc, char *argv[])
{
  char *argstr;
  char cmdchar;
  int nvals;
  int row_result;
  char outputname[MAX_BUFFER];
  char estimatename[MAX_BUFFER];
  FILE * outputHandle;
  FILE * estimateHandle = NULL;

  /* NOTE: Most filtering according to MAG_ISO was effectively removed on Jan 29, 2009 by increasing
     the MAG_ISO_REJECT_FACTOR to a high value of 3000 */
  /* The MAG_ISO_Table is a histogram of all input points as a function of MAG_ISO and spatial_bin */
  int *MAG_ISO_Table = NULL;
  /* The magnitudeTable is a histogram of all output matches as a function of (MAG_ISO,Stdmag) and spatial_bin */
  int *magnitudeTable = NULL;
  int errorFlag = 0;
  int verbose = 0;
  double limitingMagnitude = 0.0;
  double matchRadius = MATCH_RADIUS_BINS * (pGscBin->bin_size);
  time_t startTime;
  time_t curTime;
  int numCacheLines; /* Number of input file lines to keep at any one time */
  CACHEROW **inputCache; /* The cache itself */
  CACHEROW *pCacheRow;
  int numRows = 0;
  int memoryAllocated = 0;
  int index;
  int maxStarsPerRow = 0;
  double minDec;
  double maxDec = 0;
  double minRa[2];
  double maxRa[2];
  int curCatalogDecBin = -1;
  int oldCatalogDecBin = -1;
  int tmpCatalogDecBin;
  double curCatalogDec = 0;
  double maxCatalogDec;
  double minCatalogDec;
  double scale = 0.0; /* Plate scale in arcsec per pixel */

  char indexname[MAX_BUFFER];
  File indexHandle;
  char* catalogdir;
  char* slashPtr;
  char qualifier[MAX_BUFFER];
  char catalogname[MAX_BUFFER];
  int catalogNumber = 0;
  File catalogHandle;
  PSTARINDEX pStarIndex = NULL;
  int starIndexAlloc = 0;
  int starIndexEntries;
  GSCIMAGEBIN *pGscImage = NULL;
  PGSCIMAGE pGscStorage = NULL;
  int gscImageAlloc = 0;
  int gscImageEntries;
  char * dotPtr;
  int matchedStars = 0;
  double geocentricJD = 0.0;
  double plateepoch;

  File input_handle = NULL;
  char input_name[MAX_BUFFER];
  TableHead input_header = NULL;
  size_t input_nrecs = 0;
  TblDescriptor input_descriptor = NULL;
  TableRow input_row = NULL;
  INPUTDATA input_block;
  INPUTDATA *pInput = &input_block;
  signed long long curMemory = 0;
  signed long long maxMemory = 0;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int properMotionErrLimCount = 0;
  int duplicateCount = 0;

  memset(pInput, 0, sizeof(INPUTDATA));

  outputname[0] = 0;
  input_name[0] = 0;
  estimatename[0] = 0;
  qualifier[0] = 0;

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: unqualified argument %s\n", argstr);
    } else {
      while ((cmdchar = *++argstr) != '\0') {
        switch(cmdchar) {
        case 'v':
        case 'V':
          verbose = 1;
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

        case 's': // plate scale
        case 'S':
          argc--;
          nvals = sscanf(*++argv,"%lf",&scale);
          if (nvals != 1) {
            printf("ERROR: Can not decode the plate scale\n");
            errorFlag = 1;
          }
          break;

        case 'j': /* geocentric Julian Day */
        case 'J':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&geocentricJD);
            if (nvals != 1) {
              printf("ERROR: Unable to decode geocentric Julian Day %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'q': // file name qualifier (= reference catalog selection)
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(qualifier,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
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
            strncpy(outputname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'i': /* input file name */
        case 'I':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(input_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'e': /* MAG_ISO estimate file name */
        case 'E':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(estimatename,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'm': /* limiting magnitude */
        case 'M':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&limitingMagnitude);
            if (nvals != 1) {
              printf("ERROR: Unable to decode geocentric Julian Day %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        default:
          printf("ERROR: illegal command -%c-", cmdchar);
          errorFlag = 1;
          break;
        }
      }
    }
  }

  /* Verify that we have everything */
  if (outputname[0] == 0) {
    printf("ERROR: No output file specified\n");
    errorFlag = 1;
  }

  if (input_name[0] == 0) {
    printf("ERROR: No input file specified\n");
    errorFlag = 1;
  }

  if (geocentricJD == 0.0) {
    printf("ERROR: No geocentric Julian Date specified\n");
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

  if (scale <= 0.0) {
    printf("ERROR: no scale in matchstars for %s\n", outputname);
    errorFlag = 1;
  }

  plateepoch = jd2ep(geocentricJD) - GSC_EQUINOX;

  /* Now verify that we have all of the files */
  input_handle = Open(input_name, "rt");

  if (input_handle == NULL) {
    printf("ERROR: Could not open file %s\n", input_name);
    errorFlag = 1;
  } else {
    input_header = table_header(input_handle,TABLE_PARSE);
    if (input_header == NULL) {
      printf("ERROR: Failed to read header for %s\n",input_name);
      errorFlag = 1;
    } else {
      input_descriptor = table_create_descrip(
        &input_nrecs,
        TblDbl, "MAG_ISO", TblOff(INPUTDATA *, MAG_ISO),
        TblDbl, "X_IMAGE", TblOff(INPUTDATA *, X_IMAGE),
        TblDbl, "Y_IMAGE", TblOff(INPUTDATA *, Y_IMAGE),
        TblDbl, "ra", TblOff(INPUTDATA *, ra),
        TblDbl, "dec", TblOff(INPUTDATA *, dec),
        TblDbl, "FWHM_IMAGE", TblOff(INPUTDATA *, FWHM_IMAGE),
        TblDbl, "FWHM_WORLD", TblOff(INPUTDATA *, FWHM_WORLD),
        TblDbl, "THRESHOLD", TblOff(INPUTDATA *, THRESHOLD),
        TblDbl, "FLUX_MAX", TblOff(INPUTDATA *, FLUX_MAX),
        TblInt, "NUMBER", TblOff(INPUTDATA *, NUMBER),
        TblInt, "AFLAGS", TblOff(INPUTDATA *, AFLAGS),
        TblInt, "BFLAGS", TblOff(INPUTDATA *, BFLAGS),
        TblDbl, "aLength", TblOff(INPUTDATA *, aLength),
        TblDbl, "bLength", TblOff(INPUTDATA *, bLength),
        TblDbl, "THETA_J2000", TblOff(INPUTDATA *, THETA_J2000),
        0, "end", 0
      );

      if (input_descriptor == NULL) {
        printf("ERROR: Failed to allocate descriptor for %s\n", input_name);
        return -1;
      }

      table_loadmap(input_header, input_descriptor);
    }
  }

  catalogdir = getenv("DASCH_CATALOG");
  if (catalogdir == NULL) {
    printf("DASCH_CATALOG is not defined\n");
    return 1;
  }

  if (qualifier[0] != '\0') {
    catalogNumber = GetCatalogNumber(qualifier);
  }

  strcpy(catalogname, catalogdir);

  if (strstr(qualifier,"kepler")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"kepler.dat");
  } else if (strstr(qualifier,"apass")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"apass.dat");
  } else if (strstr(qualifier,"atlas")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"atlas.dat");
  } else if (strstr(qualifier,"gaia")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"gaiadr2.dat");
  } else if (strstr(qualifier,"experimental")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"experimental.dat");
  }

  strcpy(indexname, catalogname);
  dotPtr = strrchr(indexname, '.');
  if (dotPtr != NULL) {
    *dotPtr = 0;
  }
  strcat(indexname, ".idx");

  catalogHandle = Open(catalogname, "r");
  if (catalogHandle == NULL) {
    printf("ERROR Could not open catalog file %s\n", catalogname);
    errorFlag = 1;
  }

  indexHandle = Open(indexname, "r");
  if (indexHandle == NULL) {
    printf("ERROR Could not open catalog index file %s\n", indexname);
    errorFlag = 1;
  }

  outputHandle = fopen(outputname, "wt");
  if (outputHandle == NULL) {
    printf("ERROR Failed to open the output file %s\n", outputname);
    errorFlag = 1;
  }

  // Default buffer size is 8192 -- make it much larger
  setvbuf(outputHandle, NULL, _IOFBF, 1048576);

  if (estimatename[0] != 0) {
    estimateHandle = fopen(estimatename, "wt");
    if (outputHandle == NULL) {
      printf("ERROR Failed to open the estimate file %s\n", estimatename);
      errorFlag = 1;
    }

    magnitudeTable = (int *) calloc(MAGNITUDE_BINS * MAX_SPATIAL_BINS, sizeof(int));
    CalcMemory(1, sizeof(int) * MAGNITUDE_BINS * MAX_SPATIAL_BINS, &curMemory, &maxMemory);

    if (magnitudeTable == NULL) {
      printf("ERROR: Failed to allocate the magnitudesTable\n");
      errorFlag = 1;
    }

    MAG_ISO_Table = (int *) calloc(MAG_ISO_BINS * MAX_SPATIAL_BINS, sizeof(int));
    CalcMemory(1, sizeof(int) * MAG_ISO_BINS * MAX_SPATIAL_BINS, &curMemory, &maxMemory);

    if (MAG_ISO_Table == NULL) {
      printf("ERROR: Failed to allocate the MAG_ISO_Table\n");
      errorFlag = 1;
    }
  }

  if (errorFlag) {
    printf("Usage: matchstars -i <input name> \n");
    printf("                  -o <output name> \n");
    printf("                  -r <match radius dd:mm:ss.sss> \n");
    printf("                  -m limiting magnitude (default: all magnitudes)\n");
    printf("                  -s scale (arcsec/pixel) (default: match on match radius alone) \n");
    printf("                  -w is the mosaic width in pixels\n");
    printf("                  -h is the mosaic height in pixels\n");
    printf("                  -e <limiting magnitude estimate file\n");
    printf("                  -q <catalog>\n");
    printf("                  -j Julian date of the plate");
    printf("                  -v verbose\n");
    return 1;
  }

  printf(
    "matchstars of %s %s \n"
    " Input Filename %s\n"
    " Output Filename %s\n"
    " Estimate Filename %s\n"
    " Catalog file %s\n"
    " Limiting Magnitude %5.2f, Plate scale %f arcsec/pixel\n",
    __DATE__,
    __TIME__,
    input_name,
    outputname,
    estimatename,
    catalogname,
    limitingMagnitude,
    scale
  );

  numCacheLines = ((2 * (matchRadius + properMotionDegrees)) / pGscBin->bin_size) + 3;
  printf("Lines of input file to cache in memory: %d\n", numCacheLines);
  inputCache = (CACHEROW **) calloc(numCacheLines, sizeof(CACHEROW *));
  CalcMemory(1, numCacheLines * sizeof(CACHEROW *), &curMemory, &maxMemory);

  time(&startTime);

  fprintf(
    outputHandle,
    "NUMBER\t"
    "MAG_ISO\t"
    "X_IMAGE\t"
    "Y_IMAGE\t"
    "ra_1\t"
    "dec_1\t"
    "FWHM_IMAGE\t"
    "FWHM_WORLD\t"
    "AFLAGS\t"
    "BFLAGS\t"
    "REF\t"
    "ra_2\t"
    "dec_2\t"
    "RaSigmaPM\t"
    "DecSigmaPM\t"
    "Stdmag\t"
    "color\t"
    "VFlag\t"
    "MAGFlag\t"
    "class\t"
    "dra\t"
    "ddec\t"
    "drad\t"
    "THRESHOLD\t"
    "FLUX_MAX\t"
    "aLength\t"
    "bLength\t"
    "THETA_J2000\t"
    "gsc_bin_index\t"
    "spatial_bin\n"
  );
  fprintf(
    outputHandle,
    "------\t"
    "-------\t"
    "-------\t"
    "-------\t"
    "----\t"
    "-----\t"
    "----------\t"
    "----------\t"
    "------\t"
    "------\t"
    "---\t"
    "----\t"
    "-----\t"
    "---------\t"
    "----------\t"
    "------\t"
    "-----\t"
    "-----\t"
    "-------\t"
    "-----\t"
    "---\t"
    "----\t"
    "----\t"
    "---------\t"
    "--------\t"
    "-------\t"
    "-------\t"
    "-----------\t"
    "-------------\t"
    "-----------\n"
  );

  // Start reading

  while (
    (row_result = ReadNextLine(
      pInput,
      input_handle,
      input_header,
      input_descriptor,
      &input_row,
      inputCache,
      numCacheLines,
      &input_nrecs,
      verbose,
      outputname,
      &curMemory,
      &maxMemory,
      mosaicWidth,
      mosaicHeight,
      MAG_ISO_Table)
    ) != 0
  ) {
    numRows++;
    if (verbose && (numRows % 10) == 0) {
      printf("At record %zu, row %d\n", input_nrecs, numRows);
    }

    /* Figure out what part of the sky we have */
    GetRanges(inputCache, numCacheLines, &minDec, &maxDec, &minRa[0], &maxRa[0], &minRa[1], &maxRa[1], outputname);
    curCatalogDec = minDec - matchRadius - properMotionDegrees;
    if (curCatalogDec < -90.0) {
      curCatalogDec = -90.0;
    }

    if (numRows < numCacheLines) {
      /* Be sure we fill up our cache */
      continue;
    }

    if (curCatalogDecBin == -1) {
      curCatalogDecBin = GetDecBin(pGscBin, curCatalogDec, outputname);
    } else {
      curCatalogDecBin = oldCatalogDecBin + 1;
      tmpCatalogDecBin = GetDecBin(pGscBin, curCatalogDec, outputname);
      if (tmpCatalogDecBin > curCatalogDecBin) {
        /* Here there was nothing in the image for a number of bins */
        curCatalogDecBin = tmpCatalogDecBin;
      }
    }

    while (1) {
      maxCatalogDec = ((curCatalogDecBin + 1) * pGscBin->bin_size) - 90.0;
      if (maxCatalogDec > 90.0) {
        maxCatalogDec = 90.0;
      }

      if ((maxCatalogDec + properMotionDegrees + matchRadius) > maxDec) {
        /* Need to read in another bin from the sextractor file */
        break;
      }

      if (
        ReadCatalogDecBin(
          catalogHandle,
          indexHandle,
          &pStarIndex,
          &starIndexAlloc,
          &starIndexEntries,
          &pGscImage,
          &pGscStorage,
          &gscImageAlloc,
          &gscImageEntries,
          curCatalogDecBin,
          minRa,
          maxRa,
          limitingMagnitude,
          outputname,
          plateepoch,
          &properMotionErrLimCount,
          &duplicateCount,
          catalogNumber,
          &curMemory,
          &maxMemory
        )
      ) {
        SearchCatalogDecBin(
          outputHandle,
          pGscImage,
          gscImageEntries,
          inputCache,
          numCacheLines,
          scale,
          curCatalogDecBin,
          limitingMagnitude,
          plateepoch,
          &matchedStars,
          &properMotionErrLimCount,
          input_name,
          magnitudeTable,
          catalogNumber
        );
        oldCatalogDecBin = curCatalogDecBin;
      }

      /* Try to read in another row from the catalog */
      curCatalogDecBin++;
    }
  }

  /* At this point, we are out of new sextractor rows, so just continue with the catalog */

  if (curCatalogDecBin == -1) {
    curCatalogDecBin = GetDecBin(pGscBin, curCatalogDec, outputname);
  } else {
    curCatalogDecBin = oldCatalogDecBin;
  }

  while (1) {
    /* Go to the next catalog bin */
    curCatalogDecBin++;
    minCatalogDec = ((curCatalogDecBin) * pGscBin->bin_size) - 90.0;
    if (minCatalogDec >= 90.0 || (minCatalogDec - properMotionDegrees) >= (maxDec + matchRadius)) {
      /* We are done */
      break;
    }

    if (
      ReadCatalogDecBin(
        catalogHandle,
        indexHandle,
        &pStarIndex,
        &starIndexAlloc,
        &starIndexEntries,
        &pGscImage,
        &pGscStorage,
        &gscImageAlloc,
        &gscImageEntries,
        curCatalogDecBin,
        minRa,
        maxRa,
        limitingMagnitude,
        outputname,
        plateepoch,
        &properMotionErrLimCount,
        &duplicateCount,
        catalogNumber,
        &curMemory,
        &maxMemory
      )
    ) {
      SearchCatalogDecBin(
        outputHandle,
        pGscImage,
        gscImageEntries,
        inputCache,
        numCacheLines,
        scale,
        curCatalogDecBin,
        limitingMagnitude,
        plateepoch,
        &matchedStars,
        &properMotionErrLimCount,
        input_name,
        magnitudeTable,
        catalogNumber
      );
    }
  }

  time(&curTime);
  curTime -= startTime;

  if (inputCache != NULL) {
    memoryAllocated = numCacheLines * sizeof(CACHEROW *);

    for (index = 0; index < numCacheLines; index++) {
      pCacheRow = inputCache[index];

      if (pCacheRow != NULL) {
        memoryAllocated += sizeof(CACHEROW);

        if (pCacheRow->inputdata != NULL) {
          if (pCacheRow->allocEntries > maxStarsPerRow) {
            maxStarsPerRow = pCacheRow->allocEntries;
          }

          memoryAllocated += pCacheRow->allocEntries * sizeof(INPUTDATA);
          free(pCacheRow->inputdata);
        }

        free(pCacheRow);
      }
    }

    free(inputCache);
  }

  if (pStarIndex != NULL) {
    memoryAllocated += starIndexAlloc * sizeof(STARINDEX);
    free(pStarIndex);
  }

  if (pGscImage != NULL) {
    memoryAllocated += gscImageAlloc * sizeof(GSCIMAGEBIN);
    free(pGscImage);
  }

  if (pGscStorage != NULL) {
    memoryAllocated += gscImageAlloc * sizeof(GSCIMAGE);
    free(pGscStorage);
  }

  if (magnitudeTable != NULL) {
    WriteMagnitudeTable(magnitudeTable, MAG_ISO_Table, estimateHandle, estimatename);

    memoryAllocated += sizeof(int) * MAGNITUDE_BINS * MAX_SPATIAL_BINS;
    free(magnitudeTable);

    memoryAllocated += sizeof(int) * MAG_ISO_BINS * MAX_SPATIAL_BINS;
    free(MAG_ISO_Table);

    if (estimateHandle != NULL) {
      fclose(estimateHandle);
    }
  }

  printf(
    "Execution time %ld seconds. Memory allocated %d, input stars %zu, input rows %d "
    "max stars per row %d matchedStars %d zero magnitude stars %d "
    "properMotionErrLimCount %d duplicateCount %d maxMemory %lld for %s\n",
    curTime,
    memoryAllocated,
    input_nrecs,
    numRows,
    maxStarsPerRow,
    matchedStars,
    zeroMagnitudeCount,
    properMotionErrLimCount,
    duplicateCount,
    maxMemory,
    outputname
  );

  if (input_header != NULL) {
    table_hdrfree(input_header);
  }

  if (input_handle != NULL) {
    Close(input_handle);
  }

  if (input_descriptor != NULL) {
    Free(input_descriptor);
  }

  if (input_row != NULL) {
    table_rowfree(input_row);
  }

  Close(catalogHandle);
  Close(indexHandle);

  if (outputHandle != NULL) {
    fclose(outputHandle);
  }

  return 0;
}


static void
WriteMagnitudeTable(int *magnitudeTable, int *MAG_ISO_Table, FILE *estimateHandle, char *estimatename)
{
  int iStdmag;
  int iMAG_ISO;
  int binCount;
  double MAG_ISO;
  double Stdmag;
  int spatial_bin = 1;
  char *suffixPtr;
  char limit1Name[MAX_BUFFER];
  FILE *limit1Handle;
  char limit2Name[MAX_BUFFER];
  FILE *limit2Handle;
  char limit3Name[MAX_BUFFER];
  FILE *limit3Handle;
  int MAG_ISO_Limit;
  int MAG_ISO_Current;
  int MAG_ISO_Saved;
  int MAG_ISO_Accept = 0;
  int MAG_ISO_Reject = 0;
  double MAX_Stdmag;
  double StdmagTable1[MAG_ISO_BINS * MAX_SPATIAL_BINS];
  double StdmagTable2[MAG_ISO_BINS * MAX_SPATIAL_BINS];

  strcpy(limit1Name, estimatename);
  suffixPtr = strstr(limit1Name, "_lim.db");
  if (suffixPtr == NULL) {
    printf("ERROR: estimatename %s has unknown suffix\n", estimatename);
    exit(1);
  }

  *suffixPtr = 0;
  strcat(limit1Name, "_lim1.db");
  limit1Handle = fopen(limit1Name, "wt");
  if (limit1Handle == NULL) {
    printf("ERROR: failed to open %s\n", limit1Name);
    exit(1);
  }

  strcpy(limit2Name, estimatename);
  suffixPtr = strstr(limit2Name, "_lim.db");
  if (suffixPtr == NULL) {
    printf("ERROR: estimatename %s has unknown suffix\n", estimatename);
    exit(1);
  }

  *suffixPtr = 0;
  strcat(limit2Name, "_lim2.db");
  limit2Handle = fopen(limit2Name, "wt");
  if (limit2Handle == NULL) {
    printf("ERROR: failed to open %s\n", limit2Name);
    exit(1);
  }

  strcpy(limit3Name, estimatename);
  suffixPtr = strstr(limit3Name, "_lim.db");
  if (suffixPtr == NULL) {
    printf("ERROR: estimatename %s has unknown suffix\n", estimatename);
    exit(1);
  }

  *suffixPtr = 0;
  strcat(limit3Name, "_lim3.db");
  limit3Handle = fopen(limit3Name, "wt");
  if (limit3Handle == NULL) {
    printf("ERROR: failed to open %s\n", limit3Name);
    exit(1);
  }

  // For each MAG_ISO, find the Stdmag at which we reject all dimmer refcat stars

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    for (iMAG_ISO = 0; iMAG_ISO < MAG_ISO_BINS; iMAG_ISO++) {
      MAG_ISO_Limit = MAG_ISO_Table[iMAG_ISO + (MAG_ISO_BINS * (spatial_bin - 1))] * MAG_ISO_REJECT_FACTOR;
      MAG_ISO_Current = 0;
      MAG_ISO_Saved = 0;
      MAX_Stdmag = MAX_STDMAG;

      for (iStdmag = 0; iStdmag < STDMAG_BINS; iStdmag++) {
        binCount = magnitudeTable[iMAG_ISO + (MAG_ISO_BINS * iStdmag) + ((spatial_bin - 1) * MAGNITUDE_BINS)];
        MAG_ISO_Current += binCount;

        if (MAG_ISO_Current > MAG_ISO_Limit) {
          MAG_ISO_Reject += binCount;
          magnitudeTable[iMAG_ISO + (MAG_ISO_BINS * iStdmag) + ((spatial_bin - 1) * MAGNITUDE_BINS)] = 0;
        } else {
          Stdmag = MIN_STDMAG + ((1.0 * (iStdmag + 1)) * ((MAX_STDMAG - MIN_STDMAG) / STDMAG_BINS));
          MAX_Stdmag = Stdmag;
          MAG_ISO_Accept += binCount;
          MAG_ISO_Saved += binCount;
        }
      }

      StdmagTable1[iMAG_ISO + (MAG_ISO_BINS * (spatial_bin - 1))] = MAX_Stdmag;

      /* Now repeat the procedure, but save the median Stdmag of the remaining stars */
      StdmagTable2[iMAG_ISO + (MAG_ISO_BINS * (spatial_bin - 1))] = MAX_STDMAG;
      MAG_ISO_Current = 0;

      for (iStdmag = 0; iStdmag < STDMAG_BINS; iStdmag++) {
        binCount = magnitudeTable[iMAG_ISO + (MAG_ISO_BINS * iStdmag) + ((spatial_bin - 1) * MAGNITUDE_BINS)];
        MAG_ISO_Current += binCount;

        if (MAG_ISO_Current > MAG_ISO_Saved / 2) {
          Stdmag = MIN_STDMAG + ((1.0 * (iStdmag + 1)) * ((MAX_STDMAG - MIN_STDMAG) / STDMAG_BINS));
          StdmagTable2[iMAG_ISO + (MAG_ISO_BINS * (spatial_bin - 1))] = Stdmag;
          break;
        }
      }
    }
  }

  spatial_bin = 1;

  for (iMAG_ISO = 0; iMAG_ISO < MAG_ISO_BINS; iMAG_ISO++) {
    MAG_ISO = MIN_MAG_ISO + ((1.0 * iMAG_ISO) * ((MAX_MAG_ISO - MIN_MAG_ISO) / MAG_ISO_BINS));

    for (iStdmag = 0; iStdmag < STDMAG_BINS; iStdmag++) {
      Stdmag = MIN_STDMAG + ((1.0 * iStdmag) * ((MAX_STDMAG - MIN_STDMAG) / STDMAG_BINS));
      binCount = magnitudeTable[iMAG_ISO + (MAG_ISO_BINS * iStdmag) + ((spatial_bin - 1) * MAGNITUDE_BINS)];
      fprintf(limit1Handle, "%f %f %d\n", Stdmag, MAG_ISO, binCount);
    }

    fprintf(limit2Handle, "%f %f\n", StdmagTable1[iMAG_ISO], MAG_ISO);
    fprintf(limit3Handle, "%f %f\n", StdmagTable2[iMAG_ISO], MAG_ISO);
  }

  fprintf(estimateHandle, "spatial_bin\tiMAG_ISO\tStdmag\tMAG_ISO\tmedianStdmag\n");
  fprintf(estimateHandle, "-----------\t--------\t------\t-------\t------------\n");

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    for (iMAG_ISO = 0; iMAG_ISO < MAG_ISO_BINS; iMAG_ISO++) {
      MAG_ISO = MIN_MAG_ISO + ((1.0 * iMAG_ISO) * ((MAX_MAG_ISO - MIN_MAG_ISO) / MAG_ISO_BINS));
      fprintf(
        estimateHandle,
        "%d\t%d\t%f\t%f\t%f\n",
        spatial_bin,
        iMAG_ISO,
        StdmagTable1[iMAG_ISO],
        MAG_ISO,
        StdmagTable2[iMAG_ISO]
      );
    }
  }

  fclose(limit1Handle);
  fclose(limit2Handle);
  fclose(limit3Handle);

  printf(
    "MAG_ISO_Accept %d MAG_ISO_Reject %d MAG_ISO_REJECT_FACTOR is %d for %s\n",
    MAG_ISO_Accept,
    MAG_ISO_Reject,
    MAG_ISO_REJECT_FACTOR,
    estimatename
  );
}

/*
 *  Mar  7, 2008 Edward J. Los - Initial version
 *  Mar 14, 2008 Edward J. Los - Do not test for magnitude if it has not been specified
 *                             - Take the smaller of FWHM_WORLD and matchRadius
 *  Mar 15, 2008 Edward J. Los - Use zero scale for match on matchRadius alone.
 *  Jun 10, 2008 Edward J. Los - Get rid of all Stdmag == 0.0 GSC2.3.2 stars.  These appear to be a catalog error.
 *  Jul  1, 2008 Edward J. Los - Sextractor believes are blended
 *                               Change the search radius from FWHM_WORLD + 1 pixel
 *                                                      to (FWHM_WORLD/2) + 5 pixels
 *                               Add THRESHOLD and FLUX_MAX for drad filtering
 *  Aug 14, 2008 Edward J. Los - Correct wrap test of curIndex in SearchCatalogDecBin
 *  Sep 29, 2008 Edward J. Los - Add the plate text to the Aries error message
 *                               Bugfixes for the first point of Aries.
 *  Nov 25, 2008 Edward J. Los - Replace PIXEL_RADIUS_FACTOR with MaxDradPixels() function
 *  Jan 22, 2009 Edward J. Los - Add RaPM and DecPM, the proper motions in mas/yr.  Adjust the
 *                               GSC catalog position by these values.
 *  Jan 27, 2009 Edward J. Los - Split FLAGS into AFLAGS and BFLAGS
 *                               Expand search area to account for proper motion.
 *  Feb  4, 2009 Edward J. Los - Fix bug in variable search radius near top of plate
 *  Feb 10, 2009 Edward J. Los - Use new blend detection algorithm
 *  Feb 16, 2009 Edward J. Los - Add gsc_bin_index prior to calculating precession
 *  Feb 27, 2009 Edward J. Los - Add CalcMemory
 *  Mar  5, 2009 Edward J. Los - Correct bin sizing error "numBins is -12370 in ReadCatalogDecBin"
 *  May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 *  Jun  1, 2009 Edward J. Los - Correct gsc_bin_index error.
 *  Jun  5, 2009 Edward J. Los - Add the spatial_bin for preliminary analysis of the limiting magnitude.
 *  Jun 15, 2009 Edward J. Los - Add an estimate of the limiting magnitudes for this plate.
 *                               (This code was subsequently disabled in dasch_match.csh)
 *  Aug 19, 2009 Edward J. Los - Add kepler catalog support
 *  Mar 28, 2011 Edward J. Los - Correct bug in assignment of gsc_bin_index.
 *                               Correct improper setting of FILTER_BFLAG_PROPER_MOTION from
 *                               a previous match
 *                               Rewrite the search algorithm to correct false duplicate values
 *                               Add apass catalog support
 *  Apr 12, 2011 Edward J. Los   Fix maxPM calculation
 *  Jul 30, 2012 Edward J. Los   Add experimental catalog support
 *  Aug 27, 2013 Edward J. Los   Remove matchRadius, replacing it with three GSC bins
 *                               Do not filter GSC proper motions, all proper motions are now ucac4.
 *                               Add proper motion error support
 *  Nov 19, 2013 Edward J. Los   Correct meridian wrap
 *  May 29, 2018 Edward J. Los   Support gaia
 *  Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 */
