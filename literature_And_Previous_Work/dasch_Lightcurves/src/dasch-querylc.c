// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// dasch-querylc -- extract lightcurve data
//
// Originally derived from web_plot.c, with *lots* of changes.
//
// This needs the following environment variables:
//
// - DASCH_MYSQLHOST
// - DASCH_USERNAME
// - DASCH_PASSWORD
// - DASCH_CATALOG (../gsc232bin.dat)
// - DASCH_PHOT_ROOT

#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

// Starbase
#include <table.h>

// libwcs
#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>
#include <libwcs/wcscat.h>

// local
#include "scandb.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include "photometryutils.h"
#include "searchgsc.h"
#include "daschunistd.h"


#define MAX_BUFFER 500


extern char *catalogText[MAX_CATALOG_NUMBER];
extern GSCBIN gscBin64;

static PGSCBIN pGscBin = &gscBin64;


static int
PlateCompare(const void *first, const void *second)
{
  PFILESTARIMAGEEXT pPlateRecFirst = (PFILESTARIMAGEEXT) first;
  PFILESTARIMAGEEXT pPlateRecSecond = (PFILESTARIMAGEEXT) second;

  if (pPlateRecFirst->filestarimage.seriesId > pPlateRecSecond->filestarimage.seriesId) {
    return 1;
  } else if (pPlateRecFirst->filestarimage.seriesId < pPlateRecSecond->filestarimage.seriesId) {
    return -1;
  }

  if (pPlateRecFirst->filestarimage.plateNumber > pPlateRecSecond->filestarimage.plateNumber) {
    return 1;
  } else if (pPlateRecFirst->filestarimage.plateNumber < pPlateRecSecond->filestarimage.plateNumber) {
    return -1;
  }

  if (pPlateRecFirst->mosaicNumber > pPlateRecSecond->mosaicNumber) {
    return 1;
  } else if (pPlateRecFirst->mosaicNumber < pPlateRecSecond->mosaicNumber) {
    return -1;
  }

  if (pPlateRecFirst->filestarimage.solutionNumber > pPlateRecSecond->filestarimage.solutionNumber) {
    return 1;
  } else if (pPlateRecFirst->filestarimage.solutionNumber < pPlateRecSecond->filestarimage.solutionNumber) {
    return -1;
  }

  return 0;
}


static int
compare_lim_and_image(PLATELIMITINGREC *plc, FILESTARIMAGEEXT *fsie)
{
  if (plc->seriesId < fsie->filestarimage.seriesId)
    return -2;
  else if (plc->seriesId > fsie->filestarimage.seriesId)
    return 2;

  if (plc->plateNumber < fsie->filestarimage.plateNumber)
    return -2;
  else if (plc->plateNumber > fsie->filestarimage.plateNumber)
    return 2;

  if (plc->mosaicNumber < fsie->mosaicNumber)
    return -1;
  else if (plc->mosaicNumber > fsie->mosaicNumber)
    return 1;

  if (plc->solutionNumber < fsie->filestarimage.solutionNumber)
    return -1;
  else if (plc->solutionNumber > fsie->filestarimage.solutionNumber)
    return 1;

  return 0;
}


static int
JulianDateCompare(const void *first, const void *second)
{
  double dateFirst = ((PFILESTARIMAGEEXT) first)->filestarimage.Date;
  double dateSecond = ((PFILESTARIMAGEEXT) second)->filestarimage.Date;

  if (dateFirst > dateSecond) {
    return 1;
  } else if (dateFirst < dateSecond) {
    return -1;
  }

  return 0;
}


static void
AddLimitingMagnitudeEntry(
  PPLATELIMITINGREC pPlateLimitingRec,
  PFILESTARIMAGEEXT db_table,
  int *pNumMagnitudes,
  double centerRa,
  double centerDec,
  int ref_type
) {
  int numMagnitudes = *pNumMagnitudes;
  PFILESTARIMAGEEXT pFileStarImageExt = &db_table[numMagnitudes];
  PFILESTARIMAGE pFileStarImage = &pFileStarImageExt->filestarimage;
  PPHOTPLATES pPhotPlates = &pFileStarImageExt->photplates;

  memset(pFileStarImageExt, 0, sizeof(FILESTARIMAGEEXT));

  pFileStarImageExt->mosaicNumber = pPlateLimitingRec->mosaicNumber;
  pFileStarImageExt->quality = QUALITY_UNDETECTED;
  pFileStarImageExt->plateVersionId = pPlateLimitingRec->versionId;
  pFileStarImage->Date = pPlateLimitingRec->geoJulianDate;
  pFileStarImage->seriesId = pPlateLimitingRec->seriesId;
  pFileStarImage->plateNumber = pPlateLimitingRec->plateNumber;
  pFileStarImage->solutionNumber = pPlateLimitingRec->solutionNumber;
  pFileStarImage->limiting_mag_local = pPlateLimitingRec->limiting_mag_local;
  pFileStarImage->versionId = pPlateLimitingRec->versionId;
  pFileStarImage->timeAccuracy = -1;
  pPhotPlates->mosaicNumber = pPlateLimitingRec->mosaicNumber;
  pPhotPlates->quality = QUALITY_UNDETECTED;
  pPhotPlates->versionId = pPlateLimitingRec->versionId;

  if (ref_type == REF_TYPE_NONE) {
    pFileStarImage->ra = centerRa;
    pFileStarImage->dec = centerDec;
  }

  numMagnitudes++;
  *pNumMagnitudes = numMagnitudes;
}


int
main(int argc, char *argv[])
{
  char ref_text[MAX_REF];
  long long ref_number;
  int ref_type;
  long long ref_number_start;
  int errorFlag = 0;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int gsc_bin_index = -1;
  PFILESTARIMAGEEXT db_table = NULL;
  PFILESTARIMAGEEXT pFileStarImageExt;

  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;

  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon0 = &fileCommon;

  int numMagnitudes;
  int newNumMagnitudes;
  int magnitudeIndex;
  int curMagnitudeIndex;
  PFILESTARIMAGE pNoneMagnitudeTable = NULL;
  PPHOTPLATES pPhotPlates = NULL;
  GALAXYCOMMON galaxycommon;
  PGALAXYCOMMON pGalaxyCommon = &galaxycommon;
  double centerRa;
  double centerDec;
  PPHOTSTARIMAGE pMagnitudeTable = NULL;
  int allocMagnitudes = 0;
  int quality = 0;
  int plateIndex;
  int writeHeader = 1;

  memset(pGalaxyCommon, 0, sizeof(GALAXYCOMMON));

  SetQueryCount(0);
  qualifier[0] = 0;
  catalogString[0] = 0;

  // Parse arguments

  for (argv++; --argc > 0; argv++) {
    char *argstr = *argv;

    if (argstr[0] != '-') {
      errorFlag = 1;
      fprintf(stderr, "error: unqualified argument %s argc: %d\n", argstr, argc);
    } else {
      char cmdchar;

      while ((cmdchar = *++argstr) != 0) {
        switch (cmdchar) {

        case 'g': /* gsc_bin_index */
          argc--;

          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            if (sscanf(*++argv, "%d", &gsc_bin_index) != 1) {
              fprintf(stderr, "error: Unable to decode the gsc_bin_index %s\n", *argv);
              errorFlag = 1;
            }
          }
          break;

        case 'q': /* Catalog and file name qualifier */
          argc--;
          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              fprintf(stderr, "error: Illegal catalog name %s\n", *argv);
              errorFlag = 1;
            } else {
              if (catalogNumber > 0) {
                sprintf(catalogString,"%d",catalogNumber);
                strcpy(qualifier,*argv);
              }
            }
          }
          break;

        case 'r': /* Object designation */
          argc--;
          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_REF-1)) {
              fprintf(stderr, "error: designation length %zu for %s is too long\n", strlen(*argv), *argv);
              errorFlag = 1;
            } else {
              strcpy(ref_text,*argv);
              GetREFNumber(ref_text, &ref_number, &ref_type, 1, 1);
            }
          }
          break;

        default:
          fprintf(stderr, "error:  unknown command -%c\n", cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  if (gsc_bin_index < 0) {
    fprintf(stderr, "error: gsc_bin_index is not specified\n");
    errorFlag = 1;
  }

  if (errorFlag == 1) {
    printf("usage: dasch-querylc -r <REF> -q <catalog> -g <gsc_bin_index>\n");
    return 1;
  }

  if (GetBinCenter(pGscBin, gsc_bin_index, &centerRa, &centerDec, "web_plot") != 0) {
    fprintf(stderr, "error: auth failure??\n");
    return 1;
  }

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pConnection);
  dasch_init_photdb(pPhotConnection);

  // More init

  InitSeriesTable(pConnection,pPhotConnection);
  if (GetPhotometryGlobal(pPhotConnection, pPhotGlobal) != 1)  {
    fprintf(stderr, "error: failed to get the global photometry table\n");
    return 1;
  }

  if (pPhotGlobal->magnitudeFile != PHOT_MAGNITUDEFILE_YES) {
    fprintf(stderr, "error: database uses obsolute magnitude table\n");
    return 1;
  }

  InitFileCommon(stderr, pFileCommon0, pConnection, pPhotConnection, catalogNumber, 0);
  InitMaxPlateNumber(pConnection, pFileCommon0->maxPlateNumber);

  if (OpenGalaxyFiles(pGalaxyCommon, qualifier, NULL, GetPhotFileBase(catalogString)) == 0) {
    LoadGalaxyTable(pGalaxyCommon, centerRa, centerDec);
  }

  // "Find all of the images in this bin"

  if (ref_number == 0 && ref_type == REF_TYPE_NONE) {
    numMagnitudes = 0;
  } else {
    LocateNoneImages(
      pGscBin,
      pFileCommon0,
      gsc_bin_index,
      &pMagnitudeTable,
      &allocMagnitudes,
      &numMagnitudes,
      catalogString,
      1,
      0
    );
  }

  if (numMagnitudes > 0) {
    int noneMagnitudeAlloc = numMagnitudes + 1000;

    if (pNoneMagnitudeTable != NULL) {
      free(pNoneMagnitudeTable);
    }

    pNoneMagnitudeTable = (PFILESTARIMAGE) calloc(noneMagnitudeAlloc, sizeof(FILESTARIMAGE));
    if (pNoneMagnitudeTable == NULL) {
      fprintf(stderr, "error: failed to allocate pNoneMagnitudeTable of size %d\n", noneMagnitudeAlloc);
      return 1;
    }

    for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
      PHOTSTARIMAGE *pCurStarImage = &pMagnitudeTable[magnitudeIndex];
      FILESTARIMAGE *pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
      memcpy(pFileStarImage, pCurStarImage->pFileStarImage, sizeof(FILESTARIMAGE));
    }
  } else {
    if (ref_number != 0 || ref_type != REF_TYPE_NONE) {
      fprintf(stderr, "error: no magnitudes found at gsc bin %d\n", gsc_bin_index);
      return 1;
    }
  }

  // "Now search through the results and extract everything that applies to our
  // reference number". I don't understand this code well enough to know if
  // we can get rid of this or what :-(
  {
    int histogramTable[MAX_NONE_HISTOGRAM];
    int groupNumber = 1;

    ProcessNoneImagesX(
      pGscBin,
      pFileCommon0,
      NULL,
      NULL,
      NULL,
      pNoneMagnitudeTable,
      numMagnitudes,
      0, // verbose
      histogramTable,
      &groupNumber,
      gsc_bin_index
    );
  }

  if (ref_type == REF_TYPE_DASCH) {
    double ra;
    double dec;
    double radminimum = 90.0;

    // This is a DASCH reference number.  We need to perform astrometric checks
    // to get everything in the vicinity that overlaps

    if (GetDASCHCoordinates(ref_text, &ra, &dec, 1, REF_TYPE_DASCH) != 0) {
      return 1;
    }

    // Now assign DASCH numbers to everything in the magnitude table

    ref_number_start = 0;

    // Now find the group which is closest to our desired location

    for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
      char tmp_ref_text[MAX_REF];
      long long tmp_ref_number;
      int tmp_ref_type;
      FILESTARIMAGE *pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];

      GetREF(pFileStarImage->REFNumber, tmp_ref_text, 1, 1);
      GetREFNumber(tmp_ref_text, &tmp_ref_number, &tmp_ref_type, 1, 1);

      if (pFileStarImage->REFNumber == ref_number) {
        // We found the one we were looking for, so just exit
        ref_number_start = ref_number;
        break;
      }

      if (tmp_ref_type == REF_TYPE_DASCH) {
        double radactual = wcsdist(ra, dec, pFileStarImage->ra, pFileStarImage->dec);

        if (radactual < radminimum) {
          radminimum = radactual;
          ref_number_start = pFileStarImage->REFNumber;
        }
      }
    }

    if (ref_number_start == 0) {
      fprintf(stderr, "error: No objects found for %s\n", ref_text);
      return 1;
    }
  } else {
    // Non DASCH object. Just use its reference number
    ref_number_start = ref_number;
  }

  // Now extract all objects with our desired reference number

  curMagnitudeIndex = 0;
  db_table = calloc(numMagnitudes + pGalaxyCommon->plateCount, sizeof(FILESTARIMAGEEXT));

  if (db_table == NULL) {
    fprintf(stderr, "error: failed to allocate db_table of size %d\n", numMagnitudes);
    return 1;
  }

  for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
    FILESTARIMAGE *pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];

    if (pFileStarImage->REFNumber == ref_number_start) {
      if (curMagnitudeIndex != magnitudeIndex) {
        FILESTARIMAGE *pCurFileStarImage = &pNoneMagnitudeTable[curMagnitudeIndex];
        memcpy(pCurFileStarImage, pFileStarImage, sizeof(FILESTARIMAGE));
      }

      curMagnitudeIndex++;
    }
  }

  numMagnitudes = curMagnitudeIndex;

  for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
    FILESTARIMAGE *pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
    pFileStarImageExt = &db_table[magnitudeIndex];
    memcpy(&pFileStarImageExt->filestarimage, pFileStarImage, sizeof(FILESTARIMAGE));
    pPhotPlates = &pFileStarImageExt->photplates;

    if (
      GetPhotPlate(
        pPhotConnection,
        GetSeriesString(pFileStarImage->seriesId, 1),
        pFileStarImage->plateNumber,
        pPhotPlates,
        catalogString,
        1
      ) != 0
    ) {
      fprintf(
        stderr,
        "warning: failed to get photplates record for %s%05d\n",
        GetSeriesString(pFileStarImage->seriesId, 0),
        pFileStarImage->plateNumber
      );
      pPhotPlates->mosaicNumber = 99;
      pPhotPlates->quality = QUALITY_UNINITIALIZED;
      pPhotPlates->versionId = 0;
    }

    if (pPhotPlates->versionId == 0) {
      pPhotPlates->mosaicNumber = 99;
      pPhotPlates->quality = QUALITY_UNINITIALIZED;
      pPhotPlates->versionId = 0;
    }

    if (pPhotPlates->versionId != pFileStarImage->versionId) {
      fprintf(
        stderr,
        "warning: versionId mismatch for plate %s%05d: vid %d/%d, gsc_bin_index %d\n",
        GetSeriesString(pFileStarImage->seriesId, 0),
        pFileStarImage->plateNumber,
        pPhotPlates->versionId,
        pFileStarImage->versionId,
        pFileStarImage->gsc_bin_index
      );
      pPhotPlates->mosaicNumber = 99;
      pPhotPlates->quality = QUALITY_UNINITIALIZED;
      pPhotPlates->versionId = 0;
    }

    GetFullQuality(&pFileStarImageExt->filestarimage, pPhotPlates->quality, &quality);

    pPhotPlates->quality = quality;
    pFileStarImageExt->mosaicNumber = pPhotPlates->mosaicNumber;
    pFileStarImageExt->quality = pPhotPlates->quality;
    pFileStarImageExt->plateVersionId = pPhotPlates->versionId;
  }

  // Sort the results by series and plateNumber so that we can add limiting
  // magnitude records. Limiting magnitude records are already sorted by
  // seriesId and plateNumber

  qsort(db_table, numMagnitudes, sizeof(FILESTARIMAGEEXT), PlateCompare);

  int magnitudeIndex2 = 0;
  pFileStarImageExt = &db_table[magnitudeIndex2];
  plateIndex = 0;
  newNumMagnitudes = numMagnitudes;

  while (plateIndex < pGalaxyCommon->plateCount) {
    PLATELIMITINGREC *pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
    int code = compare_lim_and_image(pPlateLimitingRec, pFileStarImageExt);

    if (code >= -1 && code <= 1 && pPlateLimitingRec->versionId != pFileStarImageExt->filestarimage.versionId) {
      fprintf(
        stderr,
        "warning: versionId mismatch for lim/image %s%05d: vid %d/%d, mosnum %d/%d, gsc_bin_index %d\n",
        GetSeriesString(pFileStarImageExt->filestarimage.seriesId, 0),
        pFileStarImageExt->filestarimage.plateNumber,
        pPlateLimitingRec->versionId,
        pFileStarImageExt->filestarimage.versionId,
        pPlateLimitingRec->mosaicNumber,
        pFileStarImageExt->mosaicNumber,
        pFileStarImageExt->filestarimage.gsc_bin_index
      );
    }

    if (code == 0) {
      // Have a match, so we already know the limiting magnitude
      plateIndex++;

      if (plateIndex >= pGalaxyCommon->plateCount) {
        break;
      }

      pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
      magnitudeIndex2++;

      if (magnitudeIndex2 >= numMagnitudes) {
        break;
      }

      pFileStarImageExt = &db_table[magnitudeIndex2];
    } else if (code < 0) {
      AddLimitingMagnitudeEntry(pPlateLimitingRec, db_table, &newNumMagnitudes, centerRa, centerDec, ref_type);
      plateIndex++;

      if (plateIndex >= pGalaxyCommon->plateCount) {
        break;
      }

      pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
    } else {
      magnitudeIndex2++;

      if (magnitudeIndex2 >= numMagnitudes) {
        break;
      }

      pFileStarImageExt = &db_table[magnitudeIndex2];
    }
  }

  while (plateIndex < pGalaxyCommon->plateCount) {
    PLATELIMITINGREC *pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
    AddLimitingMagnitudeEntry(pPlateLimitingRec, db_table, &newNumMagnitudes, centerRa, centerDec, ref_type);
    plateIndex++;
  }

  numMagnitudes = newNumMagnitudes;

  // Now sort by JD and report

  qsort(db_table, numMagnitudes, sizeof(FILESTARIMAGEEXT), JulianDateCompare);

  for (magnitudeIndex2 = 0; magnitudeIndex2 < numMagnitudes; magnitudeIndex2++) {
    pFileStarImageExt = &db_table[magnitudeIndex2];
    FILESTARIMAGE *pFileStarImage = &pFileStarImageExt->filestarimage;
    pPhotPlates = &pFileStarImageExt->photplates;

    if (pFileStarImageExt->mosaicNumber != 99) {
      WriteStarbaseRecord(
        pFileStarImage,
        pPhotPlates,
        stdout,
        NULL, // no title
        &writeHeader,
        catalogNumber
      );
    }
  }

  // Clean up

  FreeFileCommon(pFileCommon0, 0);
  free(pMagnitudeTable);
  free(pNoneMagnitudeTable);
  mysql_close(pPhotConnection);
  mysql_close(pConnection);
  free(db_table);
  return 0;
}
