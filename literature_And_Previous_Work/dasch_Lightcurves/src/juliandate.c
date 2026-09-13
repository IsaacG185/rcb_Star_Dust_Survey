// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// "WARNING: This routine uses only exposure zero and should be used only by
// matchstars where the time accuracy is good enough for proper motion
// corrections."

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"
#include "pipelineutils.h"


#define MAX_COMMENT_STRING 132
#define MAX_CENTERSOURCE_STRING 10
#define MAX_RIGHTASCENSION_STRING 13
#define MAX_DECLINATION_STRING 16
#define MAX_DATE_STRING 25
#define MAX_FILENAME_LEN 256
#define MAX_NOTES_STRING 80
#define MAX_EXPOSURE_STRING 8


int
main(int argc, char *argv[])
{
  MYSQL my_connection;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  char filename[MAX_FILENAME_LEN];
  double julianDate;
  EXPOSURE exposureTable;
  PEXPOSURE pExposure = &exposureTable;

  if (argc < 2) {
    fprintf(stderr, "Usage: juliandate <filename>\n");
    printf("0\n");
    return 1;
  }

  strncpy(filename, argv[1], MAX_FILENAME_LEN);
  filename[MAX_FILENAME_LEN - 1] = 0;

  if (strstr(filename, "_tnx")) {
    fprintf(stderr, "Operation not permitted with files that have '_tnx' in the name\n");
    printf("0\n");
    return 1;
  }

  if (ParseFilename2(filename, series, &plateNumber) == 0) {
    fprintf(stderr, "Failed to parse filename %s as %s%05d\n", filename, series, plateNumber);
    printf("0\n");
    return 1;
  }

  dasch_init_scandb(&my_connection);

  if (GetExposureInfo(&my_connection, series, plateNumber, 0, pExposure, 0)) {
    julianDate = fd2jd(pExposure->date);
    printf("%12.4f %.6f %s%05d\n", julianDate, pExposure->timeAccuracy, series, plateNumber);
  } else {
    fprintf(stderr, "Failed to get exposure record for %s\n",filename);
    printf("0\n");
  }

  mysql_close(&my_connection);
  return 0;
}

/* 2008-02-19 Edward J. Los - add the time accuracy
 * 2008-06-21 Edward J. Los - add plate number
 */
