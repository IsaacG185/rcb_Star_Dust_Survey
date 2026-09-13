// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// Do some basic mosaic geometry calculations.
//
// Usage: getlocation [-vdgas] [-r <search pixels>] [-e <solution_number>] [-t <arcsec>] <filename>
//
// Input specifiers:
//
// -e          Solution number
// <filename>  Mosaic filename
//
// Modes:
//
// (default)     Print: centerRA(deg) centerDec(deg) maxRad(arcmin) scale(arcsec/px) searchRadius(sexagesimal)
// -a            Print: naxis1 naxis2 platescale(arsec) marginLeft marginRight marginBot marginTop
// -b            Print: (small|normal|large) series plateNum mosaicNum filename
// -s            Print: series scale(asec/px) binA1Radius(deg)
// -g            Print a `gsc23dasch` query command that could be run
// -n            Exit 0 if this plate is a candidate for a second pass at astrometry (= large scale patrol), 1 otherwise
// -t <r(asec)>  Print size in pixels: 2 * <r> / (3600 * scale)
//
// Other flags:
//
// -v       Verbose: report input filename
// -d       Output angles in decimal degrees
// -l       Output galactic coordinates in decimal degrees
// -p       Don't check whether filename rotation and database rotation agree
// -r <px>  A search radius in pixels for the search-radius output

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <mysql.h>

// cfitsio
#include <fitsio.h>
#include <longnam.h>

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
#define A1_FACTOR 0.23666  /* From annular9.m */
#define DEFAULT_SEARCH_PIXELS 40


static int ndec = 3;		/* Number of decimal places in non-angles */


int
main(int argc, char *argv[])
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
  double scale;
  int mosaicNumber;
  int rotation;
  int expRotation;
  int binning;
  char filename[MAX_FILENAME_LEN];
  char queryString[MAX_QUERY_STRING];
  int naxis1;
  int naxis2;
  char *ctype1;
  char *ctype2;
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double cra;
  double cdec;
  double lontemp;
  double lattemp;
  double era;
  double edec;
  double dist;
  double pixdist;
  double xmin;
  double ymin;
  double xctr;
  double yctr;
  struct WorldCoor *wcs;

  int totalCount = 0;
  int foundEntry = 0;
  int foundEntryScans = 0;
  int foundEntryPattern = 0;
  double cd[4];
  int degout = 0;
  int errorFlag = 0;
  char *argstr;
  char cmdchar;
  int verboseFlag = 0;
  int gsc23Flag = 0;
  int dradFlag = 0;
  int patrolFlag = 0;
  int ignoreRotation = 0;
  char rstr[32], dstr[32];
  char radstr[32];
  int axisFlag = 0;
  int plateSizeFlag = 0;
  double dRightAscension = 999;
  double dDeclination = 99;
  double nominalPlateScale = -1;
  double fittedPlateScale = -1;
  int scanNumber;
  int transform = 0;
  int WCSSource;
  int pixelRadius = DEFAULT_SEARCH_PIXELS;
  int patternID;
  int leftMargin;
  int rightMargin;
  int topMargin;
  int bottomMargin;
  int leftMargin2;
  int rightMargin2;
  int topMargin2;
  int bottomMargin2;
  int nullMessage = 0;
  int solutionNumber = 0;
  int diskLocation;
  double extractArcsec = 0; /* Extract radius for getfits */
  double extractPixels; /* Pixel width for getfits */
  double galacticout = 0;

  // Handle arguments

  filename[0] = 0;

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    if (argstr[0] != '-') {
      /* This must be the list of plates */
      if (strlen(filename) > 0) {
        errorFlag = 1;
        fprintf(stderr, "ERROR: list %s is being overwritten by %s\n", filename, argstr);
      } else {
        strncpy(filename, argstr, MAX_FILENAME_LEN - 1);
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'a': /* show only the axes and the plate scale*/
          axisFlag = 1;
          break;

        case 'b': /* show the size of the plate */
          plateSizeFlag = 1;
          break;

        case 'g': /* GSC2.3 flag */
          gsc23Flag = 1;
          break;

        case 'l': /* Generate galactic coordinates */
          galacticout = 1;
          degout = 1;
          break;

        case 'p': /* ignore rotation */
          ignoreRotation = 1;
          break;

        case 's':  /* Series information */
          dradFlag = 1;
          degout = 1;
          break;

        case 'v': /* verbose flag */
          verboseFlag = 1;
          break;

        case 'n': /* patrol flag */
          patrolFlag = 1;
          break;

        case 'd': /* degrees out */
          degout++;
          break;

        case 'r': /* search radius in pixels */
          argc--;
          if (argc < 1) {
            fprintf(stderr, "ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else if (sscanf(*++argv, "%d", &pixelRadius) != 1) {
            fprintf(stderr, "ERROR: Unable to decode the search radius %s\n", *argv);
            errorFlag = 1;
          }
          break;

        case 't': /* getfits radius in arcsec */
          argc--;
          if (argc < 1) {
            fprintf(stderr, "ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else if (sscanf(*++argv, "%lf", &extractArcsec) != 1) {
            fprintf(stderr, "ERROR: Unable to decode the getfits extraction radius %s\n", *argv);
            errorFlag = 1;
          }
          break;

        case 'e': /* solution number */
          argc--;
          if (argc < 1) {
            fprintf(stderr, "ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else if (sscanf(*++argv, "%d", &solutionNumber) != 1) {
            fprintf(stderr, "ERROR: Unable to decode the solutionNumber %s\n", *argv);
            errorFlag = 1;
          }
          break;

        default:
          fprintf(stderr, "* illegal option -%c-", cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  if (filename[0] == 0) {
    fprintf(stderr, "ERROR: No filename was specified\n");
    errorFlag = 1;
  }

  /* We do not need the plate center if the "-a" or the "-s" or "-b" qualifier is specified */
  if ((dradFlag || axisFlag) && solutionNumber) {
    solutionNumber = 0;
  }

  if (errorFlag) {
    printf("Usage: getlocation [-v] [-d] [-g] [-a] [-s] [-r <search pixels>] [-e <solution_number>] [-t <arcsec>] <filename>\n");
    printf(" Returns: center ra (degrees)\n");
    printf("          center dec (degrees)\n");
    printf("          max radius (arcmin)\n");
    printf("          scale (arcsec/pixel)\n");
    printf("          search radius (<search pixels>*scale) in sexagesimal\n");
    printf("          -a display the axis length in pixels\n");
    printf("             (width in pixels, height in pixels, arcsec/pixel, scan pattern id\n");
    printf("             left, right, bottom, top margins in pixels\n");
    printf("          -b display 'small', 'normal', or 'large'; series, plateNumber, mosaicNumber and filename\n");
    printf("             where small = 4in x 5in, normal = 8in x 10in, and large = 14in x 17in\n");
    printf("          -d outputs in degrees\n");
    printf("          -l outputs galactic coordinates (lon lat) in degrees\n");
    printf("          -g generate gsc2.3 query\n");
    printf("          -r number of pixels in search radius (default is 40)\n");
    printf("          -s display the series, arcsec/pixel, bin a1 radius in degrees\n");
    printf("          -v is the verbosity flag\n");
    printf("          -p ignore rotation\n");
    printf("          -e <solution number>\n");
    printf("          -t <arcsec> displays pixel width only for getfits\n");
    return 1;
  }

  dasch_init_scandb(&my_connection);

  // Parse filename

  if (ParseFilename(filename, series, &plateNumber, &mosaicNumber, &binning, &rotation) == 0) {
    if (ParseFilename2(filename, series, &plateNumber) == 0) {
      fprintf(
        stderr,
        "ERROR: getlocation Failed to parse filename %s as %s%05d_%02d_%02d rotation %d\n",
        filename,
        series,
        plateNumber,
        mosaicNumber,
        binning,
        rotation
      );
      return 1;
    }

    /* Here we have the series and plateNumber but not the mosaic number.  We set the binning to 1 */
    if (SelectBestMosaic(&my_connection, series, plateNumber, &mosaicNumber, &rotation) == 0) {
      fprintf(
        stderr,
        "ERROR: getlocation Failed to find a candidate mosaic for plate filename %s, series %s plateNumber %d\n",
        filename,
        series,
        plateNumber
      );
      return 1;
    }

    binning = 1;
    fprintf(
      stderr,
      "WARNING: Selecting mosaic for %s as %s%05d_%02d_%02d rotation %d\n",
      filename,
      series,
      plateNumber,
      mosaicNumber,
      binning,
      rotation
    );
  }

  if (binning != 1) {
    fprintf(stderr, "ERROR: getlocation File binning must be 1\n");
    return 1;
  }

  // Query

  if (solutionNumber > 0) {
    sprintf(
      queryString,
      "SELECT cd1_1,cd1_2,cd2_1,cd2_2,naxis1,naxis2,ctype1,ctype2,crval1,crval2,crpix1,crpix2,WCSSource+0,rotation,scanNumber,transform+0 "
      "FROM mosaics where series = '%s' and plateNumber = %d and mosaicNumber = %d and solutionNumber = %d;",
      series,
      plateNumber,
      mosaicNumber,
      solutionNumber
    );
  } else {
    sprintf(
      queryString,
      "SELECT cd1_1,cd1_2,cd2_1,cd2_2,naxis1,naxis2,ctype1,ctype2,crval1,crval2,crpix1,crpix2,WCSSource+0,rotation,scanNumber,transform+0,dRightAscension,dDeclination,nominalPlateScale,fittedPlateScale "
      "FROM mosaics INNER JOIN exposures USING (series,plateNumber,exposureNumber) INNER JOIN series using (series) "
      "where series = '%s' and plateNumber = %d and mosaicNumber = %d and solutionNumber = %d;",
      series,
      plateNumber,
      mosaicNumber,
      solutionNumber
    );
  }

  if (strlen(queryString) > MAX_QUERY_STRING) {
    fprintf(stderr, "ERROR: getlocation MAX_QUERY_STRING exceeded %zu\n", strlen(queryString));
    return 1;
  }

  res = mysql_query(&my_connection, queryString);

  if (!res) {
    res_ptr = mysql_store_result(&my_connection);

    if (res_ptr) {
      while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
        totalCount++;
        if (sqlrow[0] == NULL) {
          fprintf(stderr, "ERROR: No WCS found for %s\n",filename);
          continue;
        }

        if (sscanf(sqlrow[0], "%lf", &cd1_1) != 1) {
          fprintf(stderr, "ERROR: failed to parse cd1_1\n");
          continue;
        }

        if (sscanf(sqlrow[1], "%lf", &cd1_2) != 1) {
          fprintf(stderr, "ERROR: failed to parse cd1_2\n");
          continue;
        }

        if (sscanf(sqlrow[2], "%lf", &cd2_1) != 1) {
          fprintf(stderr, "ERROR: failed to parse cd2_1\n");
          continue;
        }

        if (sscanf(sqlrow[3], "%lf", &cd2_2) != 1) {
          fprintf(stderr, "ERROR: failed to parse cd2_2\n");
          continue;
        }

        if (sscanf(sqlrow[4], "%d", &naxis1) != 1) {
          fprintf(stderr, "ERROR: failed to parse naxis1\n");
          continue;
        }

        if (sscanf(sqlrow[5], "%d", &naxis2) != 1) {
          fprintf(stderr, "ERROR: failed to parse naxis2\n");
          continue;
        }

        ctype1 = sqlrow[6];
        ctype2 = sqlrow[7];

        if (sscanf(sqlrow[8], "%lf", &crval1) != 1) {
          fprintf(stderr, "ERROR: failed to parse crval1\n");
          continue;
        }

        if (sscanf(sqlrow[9], "%lf", &crval2) != 1) {
          fprintf(stderr, "ERROR: failed to parse crval2\n");
          continue;
        }

        if (sscanf(sqlrow[10], "%lf", &crpix1) != 1) {
          fprintf(stderr, "ERROR: failed to parse crpix1\n");
          continue;
        }

        if (sscanf(sqlrow[11], "%lf", &crpix2) != 1) {
          fprintf(stderr, "ERROR: failed to parse crpix2\n");
          continue;
        }

        if (sscanf(sqlrow[12], "%d", &WCSSource) != 1) {
          fprintf(stderr, "ERROR: failed to parse WCSSource\n");
          continue;
        }

        if (sqlrow[13] != NULL) {
          if (sscanf(sqlrow[13], "%d", &expRotation) != 1) {
            fprintf(stderr, "ERROR: failed to parse rotation\n");
            continue;
          }
        } else if (WCSSource < WCSSOURCE_IMWCS) {
          expRotation = 0;
        } else {
          fprintf(stderr, "ERROR: rotation field is null\n");
          continue;
        }

        if (sscanf(sqlrow[14], "%d", &scanNumber) != 1) {
          fprintf(stderr, "ERROR: failed to parse scanNumber\n");
          continue;
        }

        if (sqlrow[15] != NULL) {
          if (sscanf(sqlrow[15], "%d", &transform) != 1) {
            fprintf(stderr, "ERROR: failed to parse transform\n");
            continue;
          }
        } else {
          transform = 0;
          nullMessage = 1;
        }

        if (solutionNumber == 0 && plateSizeFlag == 0){
          if (sqlrow[16] == NULL) {
            dRightAscension = 999;
          } else if (sscanf(sqlrow[16], "%lf", &dRightAscension) != 1) {
            dRightAscension = 999;
          }

          if (sqlrow[17] == NULL) {
            dDeclination = 99;
          } else if (sscanf(sqlrow[17], "%lf", &dDeclination) != 1) {
            dDeclination = 99;
          }

          if (sqlrow[18]) {
            if (sscanf(sqlrow[18], "%lf", &nominalPlateScale) != 1) {
              nominalPlateScale = -1;
            }
          } else {
            nominalPlateScale = -1;
          }

          if (sqlrow[19]) {
            if (sscanf(sqlrow[19], "%lf", &fittedPlateScale) != 1) {
              fittedPlateScale = -1;
            }
          } else {
            fittedPlateScale = -1;
          }
        } else {
          dRightAscension = 0;
          dDeclination = 0;
        }

        foundEntry++;
      }

      mysql_free_result(res_ptr);
    }
  } else {
    fprintf(stderr, 
      "ERROR: getlocation Select error %d: %s res %d\n",
      mysql_errno(&my_connection),
      mysql_error(&my_connection),
      res
    );
  }

  // If `-a` mode, we need to find the scan pattern ID and then the margins for the mosaic.

  if (axisFlag) {
    sprintf(
      queryString,
      "SELECT patternID,diskLocation FROM scans where series = '%s' and plateNumber = %d and scanNumber = %d;",
      series,
      plateNumber,
      scanNumber
    );

    if (strlen(queryString) > MAX_QUERY_STRING) {
      fprintf(stderr, "ERROR: getlocation MAX_QUERY_STRING exceeded %zu\n",strlen(queryString));
      return 1;
    }

    res = mysql_query(&my_connection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(&my_connection);

      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
          totalCount++;
          if (sqlrow[0] == NULL) {
            fprintf(stderr, "ERROR: No patternID found for %s\n",filename);
            continue;
          }

          if (sscanf(sqlrow[0], "%d", &patternID) != 1) {
            fprintf(stderr, "ERROR: failed to parse patternID\n");
            continue;
          }

          if (sqlrow[1] == NULL) {
            diskLocation = 99999;
          } else if (sscanf(sqlrow[1], "%d", &diskLocation) != 1) {
            fprintf(stderr, "ERROR: failed to parse diskLocation\n");
            continue;
          }

          foundEntryScans++;
        }

        mysql_free_result(res_ptr);
      }
    }

    if (foundEntryScans != 1) {
      if (foundEntryScans != 0 || solutionNumber == 0) {
        fprintf(stderr, "ERROR: getlocation %d scan database entries found for %s\n", foundEntryScans, filename);
      }

      return 1;
    }

    // Next get the margins from the pattern table

    sprintf(
      queryString,
      "SELECT leftMargin,rightMargin,bottomMargin,topMargin FROM scanPatterns where patternID = %d;",
      patternID
    );

    if (strlen(queryString) > MAX_QUERY_STRING) {
      fprintf(stderr, "ERROR: getlocation MAX_QUERY_STRING exceeded %zu\n", strlen(queryString));
      return 1;
    }

    if (nullMessage) {
      fprintf(stderr, 
        "ERROR: mosaic transform is null for %s. Tiles are on disk %d (%s/%05d_%02d)\n",
        filename,
        diskLocation,
        series,
        plateNumber,
        scanNumber
      );
    }

    res = mysql_query(&my_connection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(&my_connection);

      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
          totalCount++;

          if (sqlrow[0] == NULL) {
            leftMargin = 0;
          } else if (sscanf(sqlrow[0], "%d", &leftMargin) != 1) {
            fprintf(stderr, "ERROR: failed to parse leftMargin\n");
            leftMargin = 0;
          }

          if (sqlrow[1] == NULL) {
            rightMargin = 0;
          } else if (sscanf(sqlrow[1], "%d", &rightMargin) != 1) {
            fprintf(stderr, "ERROR: failed to parse rightMargin\n");
            rightMargin = 0;
          }

          if (sqlrow[2] == NULL) {
            bottomMargin = 0;
          } else if (sscanf(sqlrow[2], "%d", &bottomMargin) != 1) {
            fprintf(stderr, "ERROR: failed to parse bottomMargin\n");
            bottomMargin = 0;
          }

          if (sqlrow[3] == NULL) {
            topMargin = 0;
          } else if (sscanf(sqlrow[3], "%d", &topMargin) != 1) {
            fprintf(stderr, "ERROR: failed to parse topMargin\n");
            topMargin = 0;
          }

          foundEntryPattern++;
        }

        mysql_free_result(res_ptr);
      }
    }

    if (foundEntryPattern != 1) {
      fprintf(stderr, "ERROR: getlocation %d pattern database entries found for %s\n",foundEntryPattern,filename);
      return 1;
    }
  }

  mysql_close(&my_connection);

  if (foundEntry != 1) {
    if (foundEntry != 0 || solutionNumber == 0) {
      fprintf(stderr, "ERROR: getlocation %d database entries found for %s\n", foundEntry, filename);
    }

    return 1;
  }

  if (rotation != expRotation && !ignoreRotation && !plateSizeFlag) {
    fprintf(stderr, 
      "ERROR: getlocation filename rotation %d does not agree with database rotation %d for %s\n",
      rotation,
      expRotation,
      filename
    );
    return 1;
  }

  // Compute WCS stuff

  if (strstr(ctype1, "DEC")) {
    char *tmpPtr;
    double dtmp;

    tmpPtr = ctype2;
    ctype2 = ctype1;
    ctype1 = tmpPtr;

    dtmp = crval1;
    crval1 = crval2;
    crval2 = dtmp;

    cd[0] = cd2_1;
    cd[1] = cd2_2;
    cd[2] = cd1_1;
    cd[3] = cd1_2;
  } else {
    cd[0] = cd1_1;
    cd[1] = cd1_2;
    cd[2] = cd2_1;
    cd[3] = cd2_2;
  }

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
    0,  /* cdelt1 */
    0,  /* cdelt2 */
    0,  /* crota */
    2000, /* equinox */
    0   /* epoch */
  );

  xctr = 0.5 + (0.5 * naxis1);
  yctr = 0.5 + (0.5 * naxis2);
  pix2wcs(wcs, xctr, yctr, &cra, &cdec);

  xmin = 0.5;
  ymin = 0.5;
  pix2wcs(wcs, xmin, ymin, &era, &edec);

  dist = wcsdist(cra, cdec, era, edec);
  pixdist = sqrt(((xctr - xmin) * (xctr - xmin)) + ((yctr - ymin) * (yctr - ymin)));

  if (dist > 0.0 && pixdist > 0.0) {
    scale = dist / pixdist;
  } else if (fittedPlateScale > 0) {
    scale = fittedPlateScale * NOMINAL_MM_PER_PIXEL / 3600.;
  } else {
    scale = nominalPlateScale * NOMINAL_MM_PER_PIXEL / 3600.;
  }

  if (scale < 0) {
    scale = -1;
  }

  // `-n` mode?

  if (patrolFlag) {
    if (scale < 0) {
      fprintf(stderr, 
        "ERROR: getlocation no scale is available for for %s%05d_%02d\n",
        series,
        plateNumber,
        mosaicNumber
      );
      exit(1);
    }

    if (3600.0 * scale > MIN_SCAMP_ASTROMETRY_SCALE) {
      exit(0);
    } else {
      exit(1);
    }
  }

  // Start working on what we're going to report

  if (WCSSource < WCSSOURCE_IMWCS) {
    /* Current logbook position is more reliable */
    cra = dRightAscension;
    cdec = dDeclination;
  }

  if (degout == 0) {
    ra2str(rstr, 16, cra, ndec);
    dec2str(dstr, 16, cdec, ndec - 1);
  } else {
    if (galacticout) {
      lontemp = cra;
      lattemp = cdec;
      wcscon(WCS_J2000, WCS_GALACTIC, 2000.0, 2000.0, &lontemp, &lattemp, 2000.0);
      cra = lontemp;
      cdec = lattemp;
    }

    num2str(rstr, cra, 0, ndec);
    num2str(dstr, cdec, 0, ndec);
  }

  if (scale > 0) {
    dec2str(radstr, 16, scale * (1.0 * pixelRadius), ndec);
  } else {
    strcpy(radstr, DASCH_NOTAVAILABLE);
  }

  // Which mode?

  if (axisFlag) {
    /* Perform any necessary transformations */
    if ((transform & TRANSFORM_FLIP) != 0) {
      topMargin2 = bottomMargin;
      bottomMargin2 = topMargin;
      topMargin = topMargin2;
      bottomMargin = bottomMargin2;
    }

    if ((transform & TRANSFORM_MIRROR) != 0) {
      leftMargin2 = rightMargin;
      rightMargin2 = leftMargin;
      leftMargin = leftMargin2;
      rightMargin = rightMargin2;
    }

    switch(expRotation) {
    case 0:
      break;

    case 90:
      topMargin2    = rightMargin;
      leftMargin2   = topMargin;
      bottomMargin2 = leftMargin;
      rightMargin2  = bottomMargin;
      topMargin     = topMargin2;
      leftMargin    = leftMargin2;
      bottomMargin  = bottomMargin2;
      rightMargin   = rightMargin2;
      break;

    case 180:
      topMargin2    = bottomMargin;
      leftMargin2   = rightMargin;
      bottomMargin2 = topMargin;
      rightMargin2  = leftMargin;
      topMargin     = topMargin2;
      leftMargin    = leftMargin2;
      bottomMargin  = bottomMargin2;
      rightMargin   = rightMargin2;
      break;

    case 270:
      topMargin2    = leftMargin;
      leftMargin2   = bottomMargin;
      bottomMargin2 = rightMargin;
      rightMargin2  = topMargin;
      topMargin     = topMargin2;
      leftMargin    = leftMargin2;
      bottomMargin  = bottomMargin2;
      rightMargin   = rightMargin2;
      break;

    default:
      fprintf(stderr, "ERROR: illegal rotation %d for %s\n",expRotation,filename);
      return 1;
      break;
    }

    if (scale > 0) {
      printf("%d %d %8.3f %d %d %d %d %d\n",naxis1,naxis2,3600*scale,patternID,leftMargin,rightMargin,bottomMargin,topMargin);
    } else {
      printf("%d %d %s %d %d %d %d %d\n",naxis1,naxis2,DASCH_NOTAVAILABLE,patternID,leftMargin,rightMargin,bottomMargin,topMargin);
    }
  } else if (plateSizeFlag) {
    long long mosaicSize = naxis1;

    mosaicSize = mosaicSize * naxis2;

    if (mosaicSize < SMALL_PLATE_PIXELS) {
      printf("small %s %d %d %s\n", series, plateNumber, mosaicNumber, filename);
    } else if (mosaicSize > LARGE_PLATE_PIXELS) {
      printf("large %s %d %d %s\n", series, plateNumber, mosaicNumber, filename);
    } else {
      printf("normal %s %d %d %s\n", series, plateNumber, mosaicNumber, filename);
    }
  } else if (extractArcsec != 0) {
    if (scale > 0) {
      extractPixels = 2 * extractArcsec / (3600.0 * scale);
    } else {
      fprintf(stderr, "ERROR: illegal scale %f for %s\n", scale, filename);
      extractPixels = 0;
    }

    printf("%.0f\n", extractPixels);
  } else if (dradFlag) {
    if (scale > 0) {
      printf("%s %8.3f %8.3f\n", series, 3600 * scale, dist * A1_FACTOR);
    } else {
      printf("%s %s %s\n", series, DASCH_NOTAVAILABLE, DASCH_NOTAVAILABLE);
    }
  } else if (gsc23Flag) {
    if (scale > 0) {
      printf("gsc23dasch -ra %s -dec %s -r2 %8.3f\n", rstr, dstr, 60 * dist);
    } else {
      printf("gsc23dasch -ra %s -dec %s -r2 %s\n", rstr, dstr, DASCH_NOTAVAILABLE);
    }
  } else {
    if (scale > 0) {
      printf("%s %s %8.3f %12.7f %s", rstr, dstr, 60 * dist, 3600 * scale, radstr);
    } else {
      printf("%s %s %8.3f %s %s", rstr, dstr, 60 * dist, DASCH_NOTAVAILABLE, radstr);
    }

    if (verboseFlag) {
      printf(" for %s", filename);
    }

    printf("\n");
  }

  return EXIT_SUCCESS;
}

/* Dec  5, 2007 Edward J. Los - Initial version
 * Jan 16, 2008 Edward J. Los - Add axis display
 * May 13, 2008 Edward J. Los - Fix behavior for unsolved
 *                              mosaics, using the exposures table
 *                              for the best plate centers.
 * Sep 22, 2008 Edward J. Los - Allow specification of the search
 *                              radius in pixels
 * Dec  1, 2008 Edward J. Los - Current segmentation fault for no disk
 *                              location
 * Oct  5, 2009 Edward J. Los - Add multiple exposure support
 * Apr 20, 2010 Edward J. Los - Add -t to convert arcsec radius to pixel width for getfits
 * Jul 22, 2010 Edward J. Los - Ignore the solution number if "-a" or "-s" is specified.
 * Jun 27, 2011 Edward J. Los - Add SelectBestMosaic to work with only a series and plateNumber
 * Apr 15, 2012 Edward J. Los - Add option to ignore rotation
 * Nov 24, 2012 Edward J. Los - Add galactic coordinates (approx 122 plates/sec)
 * Jul 23, 2013 Edward J. Los - Return -1 status if the scale is > 400 arcsec/mm (MIN_SCAMP_ASTROMETRY_SCALE) and the parameter is "-n"
 * Aug  3, 2017 Edward J. Los - Add -b to return the plate size, series, plateNumber, mosaicNumber, and input file name
 * Sep 11, 2017 Edward J. Los - Ignore miscellaneous errors when using the -b qualifier
 * May 11, 2017 Edward J. Los - Improve the robustness of getlocation: return "N.A." if the RA and DEC are not available
 *                              Get the scale from the series table if there is no valid wcs scale
 */
