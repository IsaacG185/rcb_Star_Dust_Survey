// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// dasch-queryunmatched -- extract photometry data for unmatched sources
//
// Derived from dasch-querylc.
//
// This needs the following environment variables:
//
// - DASCH_MYSQLHOST
// - DASCH_USERNAME
// - DASCH_PASSWORD
// - DASCH_CATALOG (../gsc232bin.dat)
// - DASCH_PHOT_ROOT
//
// Test case: if we search "V* IM Nor", there is DASCH_J153924.8-521912 at a
// separation of 17 arcsec, with 1 measurement from MF32611 with a reported
// APASS mag of 13.41. A query that should recover this is:
//
// `dasch-queryunmatched -q apass -j 234.8603 -52.3217 -r 30`

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Starbase
#include <table.h>

// libwcs
#include <libwcs/wcs.h>

// local
#include "scandb.h"
#include "pipelineutils.h"
#include "photometryutils.h"


extern GSCBIN gscBin64;

static PGSCBIN pGscBin = &gscBin64;


static int
compare_jd(const void *first, const void *second)
{
  double dateFirst = ((FILESTARIMAGEEXT *) first)->filestarimage.Date;
  double dateSecond = ((FILESTARIMAGEEXT *) second)->filestarimage.Date;

  if (dateFirst > dateSecond) {
    return 1;
  } else if (dateFirst < dateSecond) {
    return -1;
  }

  return 0;
}


static double
sep_arcsec(const double ra0, const double dec0, const double ra1, const double dec1)
{
  // This is super naive, but should be OK for our purposes. Surely wcslib
  // provides a better function but I don't know what it is. All arguments in
  // degrees; return value in arcsec.

  double factor = cos(dec0 * M_PI / 180);
  double ddec = dec1 - dec0;
  double dra = factor * (ra1 - ra0);
  return hypot(dra, ddec) * 3600;
}


int
main(int argc, char *argv[])
{
  int i, j;
  int error_flag = 0;
  int catalog_number = 0;
  char catalog_tag[12] = "";

  int gsc_bin_index = -1, gsc_ra_bin, gsc_dec_bin;
  double ra0 = 999., dec0 = 99.;
  double radius = -1.;

  MYSQL my_connection, my_phot_connection;

  PHOTGLOBAL phot_global;
  FILECOMMON file_common;

  int num_phot, alloc_phot;
  PHOTSTARIMAGE *phot_table = NULL;

  FILESTARIMAGEEXT *table = NULL;
  int num_nones;

  int write_header = 1;

  // Parse arguments

  for (argv++; --argc > 0; argv++) {
    char *argstr = *argv;

    if (argstr[0] != '-') {
      error_flag = 1;
      fprintf(stderr, "error: unqualified argument %s argc: %d\n", argstr, argc);
    } else {
      char cmdchar;

      while ((cmdchar = *++argstr) != 0) {
        switch (cmdchar) {

        // search center
        case 'j':
          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            error_flag = 1;
          } else {
            ra0 = str2ra(argv[1]);
            dec0 = str2dec(argv[2]);
            argv += 2;
            argc -= 2;
          }
          break;

        // search radius, arcsec
        case 'r':
          argc--;

          if (argc < 1) {
            fprintf(stderr, "error: insufficient arguments for -%c\n", cmdchar);
            error_flag = 1;
          } else if (sscanf(*++argv, "%lf", &radius) != 1) {
            fprintf(stderr, "error: unable to decode search radius argument \"%s\"\n", *argv);
            error_flag = 1;
          } else if (radius <= 0.) {
            fprintf(stderr, "error: illegal search radius value \"%s\"\n", *argv);
            error_flag = 1;
          }
          break;

        case 'q': /* Catalog and file name qualifier */
          argc--;
          if (argc < 1) {
            fprintf(stderr, "error: Insufficient arguments for -%c\n", cmdchar);
            error_flag = 1;
          } else {
            catalog_number = GetCatalogNumber(*++argv);
            if (catalog_number < 0) {
              fprintf(stderr, "error: Illegal catalog name %s\n", *argv);
              error_flag = 1;
            } else {
              if (catalog_number > 0) {
                sprintf(catalog_tag,"%d",catalog_number);
              }
            }
          }
          break;

        default:
          fprintf(stderr, "error:  unknown command -%c\n", cmdchar);
          error_flag = 1;
        }
      }
    }
  }

  if (ra0 == 999.) {
    fprintf(stderr, "error: -j argument is not specified\n");
    error_flag = 1;
  }

  if (radius <= 0.) {
    fprintf(stderr, "error: -r argument invalid or unspecified\n");
    error_flag = 1;
  }

  if (error_flag == 1) {
    printf("usage: dasch-queryunmatched -q <catalog> -j <RA> <dec> -r <radius(arcsec)>\n");
    return 1;
  }

  gsc_bin_index = GetGSCBin(pGscBin, ra0, dec0, &gsc_dec_bin, &gsc_ra_bin, "dasch-queryunmatched");

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(&my_connection);
  dasch_init_photdb(&my_phot_connection);

  // More init

  InitSeriesTable(&my_connection, &my_phot_connection);

  if (GetPhotometryGlobal(&my_phot_connection, &phot_global) != 1)  {
    fprintf(stderr, "error: failed to get the global photometry table\n");
    return 1;
  }

  if (phot_global.magnitudeFile != PHOT_MAGNITUDEFILE_YES) {
    fprintf(stderr, "error: database uses obsolute magnitude table\n");
    return 1;
  }

  InitFileCommon(stderr, &file_common, &my_connection, &my_phot_connection, catalog_number, 0);
  InitMaxPlateNumber(&my_connection, file_common.maxPlateNumber);

  // Load up full info

  LocateNoneImages(
    pGscBin,
    &file_common,
    gsc_bin_index,
    &phot_table,
    &alloc_phot,
    &num_phot,
    catalog_tag,
    1, // readonly
    0 // verbose
  );

  num_nones = 0;

  for (i = 0; i < num_phot; i++) {
    if (phot_table[i].pFileStarImage->REFNumber == 0) {
      num_nones++;
    }
  }

  table = (FILESTARIMAGEEXT *) calloc(num_nones, sizeof(FILESTARIMAGEEXT));
  if (table == NULL) {
    fprintf(stderr, "error: failed to allocate table of size %d\n", num_nones);
    return 1;
  }

  for (i = 0, j = 0; i < num_phot; i++) {
    PHOTSTARIMAGE *src = &phot_table[i];
    FILESTARIMAGEEXT *dest = &table[j];

    if (
      src->pFileStarImage->REFNumber == 0 &&
      sep_arcsec(ra0, dec0, src->pFileStarImage->ra, src->pFileStarImage->dec) <= radius
    ) {
      memcpy(&dest->filestarimage, src->pFileStarImage, sizeof(FILESTARIMAGE));
      j++;
    }
  }

  num_nones = j;

  // Augment

  for (i = 0; i < num_nones; i++) {
    int quality;
    FILESTARIMAGEEXT *rec = &table[i];

    if (
      GetPhotPlate(
        &my_phot_connection,
        GetSeriesString(rec->filestarimage.seriesId, 1),
        rec->filestarimage.plateNumber,
        &rec->photplates,
        catalog_tag,
        1
      ) != 0
    ) {
      fprintf(
        stderr,
        "warning: failed to get photplates record for %s%05d\n",
        GetSeriesString(rec->filestarimage.seriesId, 0),
        rec->filestarimage.plateNumber
      );
      rec->photplates.mosaicNumber = 99;
      rec->photplates.quality = QUALITY_UNINITIALIZED;
      rec->photplates.versionId = 0;
    }

    if (rec->photplates.versionId == 0) {
      rec->photplates.mosaicNumber = 99;
      rec->photplates.quality = QUALITY_UNINITIALIZED;
      rec->photplates.versionId = 0;
    }

    if (rec->photplates.versionId != rec->filestarimage.versionId) {
      fprintf(
        stderr,
        "warning: versionId mismatch for plate %s%05d: vid %d/%d, gsc_bin_index %d\n",
        GetSeriesString(rec->filestarimage.seriesId, 0),
        rec->filestarimage.plateNumber,
        rec->photplates.versionId,
        rec->filestarimage.versionId,
        rec->filestarimage.gsc_bin_index
      );
      rec->photplates.mosaicNumber = 99;
      rec->photplates.quality = QUALITY_UNINITIALIZED;
      rec->photplates.versionId = 0;
    }

    GetFullQuality(&rec->filestarimage, rec->photplates.quality, &quality);

    rec->photplates.quality = quality;
    rec->mosaicNumber = rec->photplates.mosaicNumber;
    rec->quality = rec->photplates.quality;
    rec->plateVersionId = rec->photplates.versionId;
  }

  // Now sort by JD and report

  qsort(table, num_nones, sizeof(FILESTARIMAGEEXT), compare_jd);

  for (i = 0; i < num_nones; i++) {
    FILESTARIMAGEEXT *rec = &table[i];

    if (rec->mosaicNumber != 99) {
      WriteStarbaseRecord(
        &rec->filestarimage,
        &rec->photplates,
        stdout,
        NULL, // no title
        &write_header,
        catalog_number
      );
    }
  }

  // Clean up

  FreeFileCommon(&file_common, 0);
  free(phot_table);
  free(table);
  mysql_close(&my_phot_connection);
  mysql_close(&my_connection);
  return 0;
}
