// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* This file takes a TNX FITS file and restores it to the TAN FITS using
 * the WCS values saved in the MySQL database
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mysql.h"

#include "fitsio.h"
#include "longnam.h"

#include "scandb.h"
#include "pipelineutils.h"

#define MAXPV 100
#define NAXISPV 2

char *scampKeywords[] = {
  "FGROUPNO",
  "ASTIRMS1",
  "ASTIRMS2",
  "ASTRRMS1",
  "ASTRRMS2",
  "ASTINST",
  "FLXSCALE",
  "MAGZEROP",
  "PHOTIRMS",
  "PHOTINST",
  "PHOTLINK",
  "CUNIT1",
  "CUNIT2"
};

const int numScampKeywords = sizeof(scampKeywords) / sizeof(char *);

int main(int argc, char *argv[])
{
  int res;
  MYSQL my_connection;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  double cd1_1;
  double cd2_2;
  double cd2_1;
  double cd1_2;
  int mosaicNumber;
  int rotation;
  int expRotation;
  int binning;
  int nvals;
  char filename[MAX_FILENAME_LEN];
  char queryString[MAX_QUERY_STRING];
  int index;
  int naxis1;
  int naxis2;
  char oldctype1[FLEN_KEYWORD];
  char ctype1[FLEN_KEYWORD];
  char ctype2[FLEN_KEYWORD];
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  int totalCount = 0;
  fitsfile *fptr = NULL;
  int status = 0;
  int iomode = READWRITE;
  char err_text[FLEN_ERRMSG];
  int foundEntry = 0;
  int watIndex;
  int pvIndex;
  int axisIndex;
  char keywordstr[FLEN_KEYWORD];
  int tnxFlag = 0; /* Set to 1 if a tnx file */
  int nfound = 0;
  long naxes[2];

  if (argc < 2) {
    fprintf(stderr, "Usage: restoretan <filename>\n");
    return -1;
  }

  strncpy(filename, argv[1], MAX_FILENAME_LEN);
  filename[MAX_FILENAME_LEN-1] = 0;

  if (strstr(filename, "_tnx")) {
    fprintf(stderr, "Operation not permitted with files that have '_tnx' in the name\n");
    return -1;
  }

  if (ParseFilename(filename, series, &plateNumber, &mosaicNumber, &binning, &rotation) == 0) {
    fprintf(
      stderr,
      "Failed to parse filename %s as %s%05d_%02d_%02d rotation %d\n",
      filename,
      series,
      plateNumber,
      mosaicNumber,
      binning,
      rotation
    );
    return -1;
  }

  if (binning != 1) {
    fprintf(stderr, "File binning must be 1\n");
    return -1;
  }

  dasch_init_scandb(&my_connection);

  sprintf(
    queryString,
    "SELECT cd1_1, cd1_2, cd2_1, cd2_2, naxis1, naxis2, ctype1, ctype2, crval1, crval2, crpix1, crpix2, rotation FROM mosaics where series = '%s' and plateNumber = %d and mosaicNumber = %d and solutionNumber = 0;",
    series,
    plateNumber,
    mosaicNumber
  );
  if (strlen(queryString) > MAX_QUERY_STRING) {
    fprintf(stderr, "ERROR: MAX_QUERY_STRING exceeded %zu\n", strlen(queryString));
    return -1;
  }

  res = mysql_query(&my_connection, queryString);
  if (!res) {
    res_ptr = mysql_store_result(&my_connection);

    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        totalCount++;
        if (sqlrow[0] == NULL) {
          fprintf(stderr, "No WCS found for %s\n", filename);
          continue;
        }

        nvals = sscanf(sqlrow[0], "%lf", &cd1_1);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for cd1_1\n", nvals);
          continue;
        }

        nvals = sscanf(sqlrow[1], "%lf", &cd1_2);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for cd1_2\n", nvals);
          continue;
        }

        nvals = sscanf(sqlrow[2], "%lf", &cd2_1);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for cd2_1\n", nvals);
          continue;
        }

        nvals = sscanf(sqlrow[3], "%lf", &cd2_2);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for cd2_2\n", nvals);
          continue;
        }

        nvals = sscanf(sqlrow[4], "%d", &naxis1);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for naxis1\n", nvals);
          continue;
        }

        nvals = sscanf(sqlrow[5], "%d", &naxis2);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for naxis2\n", nvals);
          continue;
        }

        strncpy(ctype1, sqlrow[6], FLEN_KEYWORD - 1);
        strncpy(ctype2, sqlrow[7], FLEN_KEYWORD - 1);

        nvals = sscanf(sqlrow[8], "%lf", &crval1);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for crval1\n", nvals);
          continue;
        }

        nvals = sscanf(sqlrow[9], "%lf", &crval2);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for crval2\n", nvals);
          continue;
        }

        nvals = sscanf(sqlrow[10], "%lf", &crpix1);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for crpix1\n", nvals);
          continue;
        }

        nvals = sscanf(sqlrow[11], "%lf", &crpix2);
        if (nvals != 1) {
          fprintf(stderr, "ERROR: nvals is %d for crpix2\n", nvals);
          continue;
        }

        if (sqlrow[12] == NULL) {
            fprintf(stderr, "ERROR: rotation is null\n");
            continue;
        } else {
          nvals = sscanf(sqlrow[12], "%d", &expRotation);
          if (nvals != 1) {
            fprintf(stderr, "ERROR: nvals is %d for rotation\n", nvals);
            continue;
          }
        }

        foundEntry++;
      }

      mysql_free_result(res_ptr);
    }
  } else {
    fprintf(stderr, "Select error %d: %s res %d\n", mysql_errno(&my_connection), mysql_error(&my_connection), res);
  }

  mysql_close(&my_connection);

  if (foundEntry != 1) {
    fprintf(stderr, "ERROR: %d database entries found\n", foundEntry);
    return -1;
  }

  if (rotation != expRotation) {
    fprintf(stderr, "ERROR: filename rotation %d does not agree with database rotation %d for %s\n", rotation, expRotation, filename);
    return -1;
  }

  /* Open the file and find out if it is TAN or TNX */
  fits_open_file(&fptr, filename, iomode, &status);

  if (status != 0) {
    fprintf(stderr, "Cannot read FITS file %s\n", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
  if (status != 0) {
    fprintf(stderr, "Cannot read NAXIS keywords in %s\n", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  if (nfound != 2) {
    fprintf(stderr, "Found only %d axes in %s\n", nfound, filename);
    return -1;
  }

  fits_read_key(fptr, TSTRING, "CTYPE1", oldctype1, NULL, &status);
  if (status != 0) {
    fprintf(stderr, "Cannot read CTYPE1 keyword in %s\n", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  if (strstr(oldctype1, "TNX") != NULL) {
    tnxFlag = 1;
  }

  /* Restore the WCS keywords */
  fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, "X-axis coordinate type", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update ctype1\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, "Y-axis coordinate type", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update ctype2\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1, "X-axis coordinate value", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update crval1\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2, "Y-axis coordinate value", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update crval2\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "X-axis reference pixel", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update crpix1\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TDOUBLE, "CRPIX2", &crpix2, "Y-axis reference pixel", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update crpix2\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TDOUBLE, "CD1_1", &cd1_1, "Change in RA---TAN along X axis", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update cd1_1\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TDOUBLE, "CD1_2", &cd1_2, "Change in RA---TAN along Y axis", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update cd1_2\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TDOUBLE, "CD2_1", &cd2_1, "Change in DEC--TAN along X axis", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update cd2_1\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  fits_update_key(fptr, TDOUBLE, "CD2_2", &cd2_2, "Change in DEC--TAN anong Y axis", &status);
  if (status != 0) {
    fprintf(stderr, "Cannot update cd2_2\n in %s", filename);
    fits_get_errstatus(status, err_text);
    fprintf(stderr, "CFITSIO ERROR %d: %s\n", status, err_text);
    return -1;
  }

  if (tnxFlag) {
    /* Now delete residual TNX keys */
    fits_delete_key(fptr, "ORIGIN", &status);
    if (status != 0) {
      fits_get_errstatus(status, err_text);
      fprintf(stderr, "Cannot delete origin %s %d %s\n", filename, status, err_text);
      status = 0;
    }

    fits_delete_key(fptr, "DATE", &status);
    if (status != 0) {
      fits_get_errstatus(status, err_text);
      fprintf(stderr, "Cannot delete date %s %d %s\n", filename, status, err_text);
      status = 0;
    }

    fits_delete_key(fptr, "IRAF-TLM", &status);
    if (status != 0) {
      fits_get_errstatus(status, err_text);
      fprintf(stderr, "Cannot delete iraf-tlm %s %d %s\n", filename, status, err_text);
      status = 0;
    }

    fits_delete_key(fptr, "WCSDIM", &status);
    if (status != 0) {
      fits_get_errstatus(status, err_text);
      fprintf(stderr, "Cannot delete wcsdim %s %d %s\n", filename, status, err_text);
      status = 0;
    }

    fits_delete_key(fptr, "LTM1_1", &status);
    if (status != 0) {
      fits_get_errstatus(status, err_text);
      fprintf(stderr, "Cannot delete ltm1_1 %s %d %s\n", filename, status, err_text);
      status = 0;
    }

    fits_delete_key(fptr, "LTM2_2", &status);
    if (status != 0) {
      fits_get_errstatus(status, err_text);
      fprintf(stderr, "Cannot delete ltm2_2 %s %d %s\n", filename, status, err_text);
      status = 0;
    }

    for (watIndex = 0; watIndex < 3; watIndex++) {
      for (index = 0; index < 30; index++) {
        status = 0;
        sprintf(keywordstr, "WAT%d_%03d", watIndex, index);
        fits_delete_key(fptr, keywordstr, &status);
        if (status != 0 && status != KEY_NO_EXIST) {
          fits_get_errstatus(status, err_text);
          fprintf(stderr, "Error deleting key %s %d %s\n", keywordstr, status, err_text);
        }
      }
    }
  } else {
    /* This is a TAN file. Delete PV keywords if they are present */
    for (axisIndex = 0; axisIndex < NAXISPV; axisIndex++) {
      for (pvIndex = 0; pvIndex < MAXPV; pvIndex++) {
        sprintf(keywordstr, "PV%d_%d ", axisIndex+1, pvIndex);
        status = 0;
        fits_delete_key(fptr, keywordstr, &status);
        if (status != 0 && status != KEY_NO_EXIST) {
          fits_get_errstatus(status, err_text);
          fprintf(stderr, "Error deleting key %s %d %s\n", keywordstr, status, err_text);
        }
      }
    }

    /* Now delete other keywords generated by SCAMP */
    for (index = 0; index < numScampKeywords; index++) {
      status = 0;
      strcpy(keywordstr, scampKeywords[index]);
      fits_delete_key(fptr, keywordstr, &status);
      if (status != 0 && status != KEY_NO_EXIST) {
        fits_get_errstatus(status, err_text);
        fprintf(stderr, "Error deleting key %s %d %s\n", keywordstr, status, err_text);
      }
    }
  }

  status = 0;
  fits_close_file(fptr, &status);
  return EXIT_SUCCESS;
}
