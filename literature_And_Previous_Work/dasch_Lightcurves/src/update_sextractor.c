// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* Update the positions in the sextractor file with the image header */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Starbase
#include <table.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>
#include <libwcs/fitswcs.h>

// Astrometry.Net
#include <astrometry/kdtree.h>

#include "scandb.h"
#include "pipelineutils.h"

/* Sigh, wcstools doesn't export a prototype for this: */
extern struct WorldCoor *GetFITSWCS(
  char    *filename, // FITS or IRAF file name
  char    *header,   // Image FITS header
  int     verbose,   // Extra printing if =1
  double  *cra,      // Center right ascension in degrees (returned)
  double  *cdec,     // Center declination in degrees (returned)
  double  *dra,      // Right ascension half-width in degrees (returned)
  double  *ddec,     // Declination half-width in degrees (returned)
  double  *secpix,   // Arcseconds per pixel (returned)
  int     *wp,       // Image width in pixels (returned)
  int     *hp,       // Image height in pixels (returned)
  int     *sysout,   // Coordinate system to return (0=image, returned)
  double  *eqout     // Equinox to return (0=image, returned)
);


#define MAX_BUFFER 512
#define SHORT_ALLOCATION_INCREMENT 10000
#define ONE_YEAR (1.)
#define CLIP_RMS_FACTOR 1.0 /* This value should match filter_wedge.c */


static int verbose = 0;		/* verbose/debugging flag */
static double eqim = 0.0;
static int sysim = 0;


/* All we want from the filter_multiple "none" table is the sextractor number */
typedef struct _none {
  int NUMBER;
} NONE, *PNONE;

typedef struct _starimage {
  int NUMBER;
  double FLUX_ISO;
  double MAG_ISO;
  double MAGERR_ISO;
  double MAG_APER;
  double MAGERR_APER;
  double MAG_AUTO;
  double MAGERR_AUTO;
  double KRON_RADIUS;
  double BACKGROUND;
  double THRESHOLD;
  double FLUX_MAX;
  int ISOAREA_IMAGE;
  double ISOAREA_WORLD;
  double X_IMAGE;
  double Y_IMAGE;
  double ra;
  double dec;
  double THETA_IMAGE;
  double THETA_J2000;
  double ELLIPTICITY;
  double ERRTHETA_J2000;
  double FWHM_IMAGE;
  double FWHM_WORLD;
  int ISO[SEXTRACTOR_NISO];
  int AFLAGS;
  int BFLAGS;
  double plate_dra;
  double plate_ddec;
  double plate_dist;
  /* New sextractor parameters for Scamp */
  double ERRA_IMAGE;
  double ERRB_IMAGE;
  double ERRTHETA_IMAGE;
  double FLUXERR_ISO;
} SEXTRACTOR, *PSEXTRACTOR;

typedef struct _shortimage {
  int NUMBER;
  double THRESHOLD;
  double FLUX_MAX;
  double X_IMAGE;
  double Y_IMAGE;
  double THETA_IMAGE;
  double ELLIPTICITY;
  int ISO[SEXTRACTOR_NISO];
  int AFLAGS;
  int BFLAGS;
  /* Derived values: half image sizes in pixels at half brightness */
  double aLength;
  double bLength;
} SHORT, *PSHORT;


static int
ReadSextractor(
  PSEXTRACTOR pSextractor,
  File sextractor_handle,
  TableHead sextractor_header,
  TblDescriptor sextractor_descriptor,
  TableRow* sextractor_row
) {
  *sextractor_row = table_rowget(
    sextractor_handle,
    sextractor_header,
    *sextractor_row,
    NULL,
    NULL,
    0
  );

  if (*sextractor_row == NULL) {
    return 0;
  }

  if (!table_loadrow(
    sextractor_handle,
    sextractor_header,
    *sextractor_row,
    sextractor_descriptor,
    (char *) pSextractor
  )) {
    printf("ERROR: Read Sextractor table_loadrow failed\n");
    return 0;
  }

  /* Keep only the original sextractor flags */
  pSextractor->AFLAGS = 0;
  pSextractor->BFLAGS &= FILTER_BMASK_SEXTRACTOR;

  /* Now set the PSF saturated flag */
  if (
    pSextractor->ISO[4] > 0 &&
    ((1.0*pSextractor->ISO[5])/(1.0*pSextractor->ISO[4])) > PSFSATURATED_ISO &&
    (pSextractor->BACKGROUND + pSextractor->FLUX_MAX) > PSFSATURATED_FLUX
  ) {
    pSextractor->BFLAGS |= (1 << FILTER_BFLAG_PSFSATURATED);
  }

  return 1;
}


static int
ReadShort(
  PSHORT pShort,
  File short_handle,
  TableHead short_header,
  TblDescriptor short_descriptor,
  TableRow* short_row
) {
  *short_row = table_rowget(short_handle, short_header, *short_row, NULL, NULL, 0);
  if (*short_row == NULL) {
    return 0;
  }

  if (!table_loadrow(short_handle, short_header, *short_row, short_descriptor, (char *) pShort)) {
    printf("ERROR: Read Short table_loadrow failed\n");
    return 0;
  }

  /* Keep only the original sextractor flags */
  pShort->AFLAGS = 0;
  pShort->BFLAGS &= FILTER_BMASK_SEXTRACTOR;
  return 1;
}


/* Sort routine based on Y_IMAGE */
static int
ImageCompare(const void *first, const void *second)
{
  double numberFirst = ((PSHORT)first)->Y_IMAGE;
  double numberSecond = ((PSHORT)second)->Y_IMAGE;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  } else {
    return 0;
  }
}


static int
NoneCompare(const void *first, const void *second)
{
  int numberFirst = ((PNONE)first)->NUMBER;
  int numberSecond = ((PNONE)second)->NUMBER;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  } else {
    return 0;
  }
}


static int
GetEllipseDimensions(
  double FLUX_MAX,  /* Maxiumum flux of the object */
  double THRESHOLD, /* Object threshold */
  double ELLIPTICITY, /* Object ellipticity */
  int *ISO,      /* Object's ISO array */
  double *aLength,  /* Long axis semi-length (returned) */
  double *bLength  /* Short axis semi-length (returned) */
) {
  int isoArea;

  // New algorithm of July 18, 2013 because a reading of the SExtractor code
  // shows that the ISO levels are linearly spaced in flux. ISO[3.5] will be
  // half way in flux between THRESHOLD and FLUX_MAX.
  isoArea = (ISO[3] + ISO[4]) / 2;

  // Some tables have ellipticities that even reach 1.0, which blows up
  // `maxLength` in the blend search, causing every point to be tested against
  // every other point. That's a bit ridiculous.
  if (ELLIPTICITY > 0.9)
    ELLIPTICITY = 0.9;

  *aLength = sqrt((1.0 * isoArea) / (PI_VALUE * (1.0 - ELLIPTICITY)));
  *bLength = *aLength * (1.0 - ELLIPTICITY);
  return 0;
}


static void
SaveMAG_ISO(
  char *magisofile,
  char *fileroot,
  double MAG_ISO_med,
  double MAG_ISO_rms
) {
  FILE *magisoHandle = NULL;

  magisoHandle = fopen(magisofile, "wt");
  if (magisoHandle == NULL) {
    printf("ERROR: Failed to open %s\n", magisofile);
    exit(1);
  }

  fprintf(magisoHandle, "MAG_ISO_med\tMAG_ISO_rms\tPlate\n");
  fprintf(magisoHandle, "-----------\t-----------\t-----\n");
  fprintf(magisoHandle, "%.2f\t%.2f\t%s\n", MAG_ISO_med, MAG_ISO_rms, fileroot);
  fclose(magisoHandle);
}


int
main(int argc, char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  char *matchDirectory;
  char outfile[MAX_BUFFER];
  char mosaicfile[MAX_BUFFER];
  FILE *outHandle = NULL;
  int outputCount = 0;
  char magisofile[MAX_BUFFER];
  char fileroot[MAX_BUFFER];
  char tmpstring[MAX_BUFFER];
  char charVal;
  char *slashPtr = NULL;
  char *dotPtr;
  char *namePtr;

  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  size_t sextractor_nrecs = 0;
  SEXTRACTOR sextractor_record;
  PSEXTRACTOR pSextractor = &sextractor_record;
  TblDescriptor sextractor_descriptor = NULL;
  TableRow sextractor_row = NULL;

  File background_handle = NULL;
  char background_name[MAX_BUFFER];
  TableHead background_header = NULL;
  PHIGHBACKGROUND background_table = NULL;
  PHIGHBACKGROUND pBackground;
  size_t background_nrecs;
  int background_index;

  File none_handle = NULL;
  char none_name[MAX_BUFFER];
  TableHead none_header = NULL;
  PNONE none_table = NULL;
  size_t none_nrecs = 0;
  size_t none_index;
  PNONE pNone1;
  PNONE pNone2;

  File short_handle = NULL;
  char short_name[MAX_BUFFER];
  TableHead short_header = NULL;
  PSHORT short_table = NULL;
  PSHORT tmp_short_table = NULL;
  int short_allocation = 0;
  size_t short_nrecs = 0;
  int short_total_nrecs = 0;
  int short_blend_nrecs = 0;
  int short_index;
  SHORT short_record;
  PSHORT pShort = &short_record;
  TblDescriptor short_descriptor = NULL;
  TableRow short_row = NULL;
  int *AFLAGS_ARRAY = NULL;

  struct WorldCoor *wcs = NULL;
  time_t startTime;
  time_t curTime;
  char * header = NULL;
  double cra;
  double cdec;
  double dra;
  double ddec;
  double secpix;
  int wp;
  int hp;
  int skipFlag = 0;
  int skipPVUpdateFlag = 0;
  int skipMargin = 0;
  int leftMargin = -1;
  int rightMargin = -1;
  int topMargin = -1;
  int bottomMargin = -1;
  double leftMarginD;
  double rightMarginD;
  double bottomMarginD;
  double topMarginD;
  int marginCount = 0;
  int noneCount = 0;
  int updateBlendFlag = 0;
  double maxLength = 0; /* Maximum image half-length in pixels */
  int maxNUMBER = 0;
  int FitWCS;
  MYSQL my_connection;
  int solutionNumber = -1;
  EXPOSURE exposureTable;
  PEXPOSURE pExposure = &exposureTable;
  double epoch;
  char series[MAX_SERIES_STRING];
  int plateNumber;

  double *vector = NULL;
  int vector_nrecs = 0;
  int curClipCount;
  double MAG_ISO_med;
  double MAG_ISO_rms;

  time(&startTime);

  sextractor_name[0] = 0;
  short_name[0] = 0;
  outfile[0] = 0;
  background_name[0] = 0;
  mosaicfile[0] = 0;
  none_name[0] = 0;
  magisofile[0] = 0;
  fileroot[0] = 0;

  matchDirectory = getenv("DASCH_MATCH");
  if (matchDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_MATCH is not defined\n");
    return 1;
  }

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      fprintf(stderr,"ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'C': /* background file name */
        case 'c':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(background_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'm': /* mosaic fits file name */
        case 'M':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(mosaicfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'f': /* median and RMS MAG_ISO file */
        case 'F':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(magisofile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
              errorFlag = 1;
            }
            strcpy(tmpstring,magisofile);
            namePtr = tmpstring;
            /* First step to the last directory delimiter */
            while ((charVal = *namePtr) != 0) {
              if (charVal == '/') {
                slashPtr = namePtr;
              }
              namePtr++;
            }
            if (slashPtr == NULL) {
              namePtr = tmpstring;
            } else {
              namePtr = slashPtr+1;
            }
            dotPtr = strstr(namePtr,".");
            if (dotPtr != NULL) {
              *dotPtr = 0;
            }
            dotPtr = strstr(namePtr,"_magiso");
            if (dotPtr != NULL) {
              *dotPtr = 0;
            }
            strcpy(fileroot,namePtr);
            if (fileroot[0] == 0) {
              printf("ERROR: Failed to obtain fileroot from %s\n",magisofile);
              errorFlag = 1;
            }
          }
          break;

        case 'e': /* Solution Number */
        case 'E':
          argc--;
          nvals = sscanf(*++argv,"%d",&solutionNumber);
          if (nvals != 1) {
            printf("ERROR: Can not decode solutionNumber\n");
            errorFlag = 1;
          }
          break;

        case 'i': /* sextractor file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sextractor_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'n': /* filter_multiple "none"  file name */
        case 'N':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(none_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'l': /* left margin */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&leftMargin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode left margin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'r': /* right margin */
        case 'R':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&rightMargin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode right margin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'b': /* bottom margin */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&bottomMargin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode bottom margin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 't': /* top margin */
        case 'T':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&topMargin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode top margin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 's': /* skip flag */
        case 'S':
          skipFlag = 1;
          break;

        case 'p': /* skip flag */
        case 'P':
          skipPVUpdateFlag = 1;
          break;

        case 'u': /* update blend flag */
        case 'U':
          updateBlendFlag = 1;
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        default:
          fprintf(stderr,"ERROR: unknown command -%c\n",cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  /* Verify that we have everything */

  if (solutionNumber < 0) {
    fprintf(stderr,"ERROR: No solution number\n");
    errorFlag = 1;
  }

  if (outfile[0] == 0) {
    fprintf(stderr,"ERROR: No output filename was specified\n");
    errorFlag = 1;
  }

  if (mosaicfile[0] == 0) {
    fprintf(stderr,"ERROR: No mosaic filename was specified\n");
    errorFlag = 1;
  }

  if (sextractor_name[0] == 0) {
    fprintf(stderr,"ERROR: No sextractor filename was specified\n");
    errorFlag = 1;
  }

  strcpy(short_name, sextractor_name);

  if (
    leftMargin  == -1 &&
    rightMargin == -1 &&
    bottomMargin == -1 &&
    topMargin == -1
  ) {
    skipMargin = 1;
  } else if (
    leftMargin == -1 ||
    rightMargin == -1 ||
    bottomMargin == -1 ||
    topMargin == -1
  ) {
    errorFlag = 1;
    fprintf(
      stderr,
      "ERROR: One of leftMargin %d, rightMargin %d, bottomMargin %d, or topMargin %d not set\n",
      leftMargin,
      rightMargin,
      bottomMargin,
      topMargin
    );
  }

  if (updateBlendFlag) {
    short_handle = Open(short_name,"r");

    if (short_handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to find the short file %s\n",short_name);
    } else if (verbose) {
      fprintf(stderr,"Found short file %s\n",short_name);
    }
  }

  /* Open the sextractor file */
  sextractor_handle = Open(sextractor_name,"r");
  if (sextractor_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the sextractor file %s\n",sextractor_name);
  } else if (verbose) {
    fprintf(stderr,"Found sextractor file %s\n",sextractor_name);
  }

  if (none_name[0] != 0) {
    none_handle = Open(none_name,"r");

    if (none_handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to find the none file %s\n",none_name);
    } else {
      if (verbose) {
        fprintf(stderr,"Found none file %s\n",none_name);
      }
    }
  }

  if ((header = GetFITShead(mosaicfile, verbose)) == NULL) {
    fprintf(stderr,"ERROR: Failed to read mosaic file header %s\n",mosaicfile);
    errorFlag = 1;
  } else {
    dasch_init_scandb(&my_connection);

    /* Read the fits file header */
    wcs = GetFITSWCS(
      mosaicfile,
      header,
      verbose,
      &cra,
      &cdec,
      &dra,
      &ddec,
      &secpix,
      &wp,
      &hp,
      &sysim,
      &eqim
    );

    if (nowcs(wcs)) {
      wcsfree(wcs);
      wcs = NULL;
      fprintf(stderr,"ERROR: No valid WCS in the mosaic file %s\n",mosaicfile);
      errorFlag = 1;
    } else if (strstr(wcs->ctype[0],"RA") == NULL) {
      fprintf(stderr,"ERROR: ctype1 is not RA in the mosaic file %s\n",mosaicfile);
      errorFlag = 1;
    }

    if (skipPVUpdateFlag == 0) {
      if (wcs->pvfail) {
        printf("ERRORD: Significant inaccuracy likely to occur in projection for %s\n",mosaicfile);
        SetMosaicFitWCS(&my_connection,(char *)"InaccuratePV",mosaicfile,solutionNumber,1,0,&FitWCS);
      } else if (wcs->inv_x != NULL || wcs->lngcor != NULL) {
        SetMosaicFitWCS(&my_connection,(char *)"InaccuratePV",mosaicfile,solutionNumber,0,0,&FitWCS);
      }
    }

    /* Now check the epoch against the latest exposure date */
    if (fileroot[0] == 0) {
      strcpy(tmpstring, mosaicfile);
    } else {
      strcpy(tmpstring, fileroot);
    }

    if (ParseFilename2(tmpstring, series, &plateNumber) > 0) {
      if (GetExposureInfo(&my_connection, series, plateNumber, 0, pExposure, 0) > 0) {
        epoch = fd2ep(pExposure->date);

        if ((epoch - wcs->epoch) > ONE_YEAR || (wcs->epoch - epoch) > ONE_YEAR) {
          printf("ERROR: mosaic epoch %f not equal to logbook epoch %f for %s\n", wcs->epoch, epoch, mosaicfile);
        }
      } else {
        printf("ERROR: GetExposureInfo failed for %s\n", mosaicfile);
      }
    } else {
      printf("ERROR: ParseFilename2 failed for %s\n", mosaicfile);
    }

    mysql_close(&my_connection);
  }

  if (skipMargin == 0) {
    /* Set our margin limits */
    leftMarginD = 1.0 * leftMargin;
    rightMarginD = 1.0 * (wp - rightMargin);
    bottomMarginD = 1.0 * bottomMargin;
    topMarginD = 1.0 * (hp - topMargin);
  }

  outHandle = fopen(outfile,"wt");
  if (outHandle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
  } else if (verbose) {
    fprintf(stderr,"Output file %s\n",outfile);
  }

  if (background_name[0] != 0) {
    background_handle = Open(background_name,"r");
    if (background_handle == NULL) {
      fprintf(stderr,"No background file found  %s\n",background_name);
    } else if (verbose) {
      fprintf(stderr,"Found plate background file %s\n",background_name);
    }
  }

  if (errorFlag) {
    fprintf(stderr,"Usage: update_sextractor -i <sextractor file> \n");
    fprintf(stderr,"                    -o <output file> \n");
    fprintf(stderr,"                    -e <solution number> \n");
    fprintf(stderr,"                    -m <fits mosaic file> \n");
    fprintf(stderr,"                    -n <filter_multiple none file\n");
    fprintf(stderr,"                    -f <output of median and rms MAG_ISO for use by estimate_lightcurves\n");
    fprintf(stderr,"                    -s skip ra and dec conversion \n");
    fprintf(stderr,"                    -l <left margin in pixels> \n");
    fprintf(stderr,"                    -r <right margin in pixels> \n");
    fprintf(stderr,"                    -b <bottom margin in pixels> \n");
    fprintf(stderr,"                    -t <top margin in pixels> \n");
    fprintf(stderr,"                    -u Update the blend flag\n");
    fprintf(stderr,"                    -p Skip update of InaccuratePV flag\n");
    fprintf(stderr,"                    -c <background file>\n");
    fprintf(stderr,"                    -v  verbose\n");
    return 1;
  }

  printf(
    "update_sextractor of %s %s skipFlag %d updateBlendFlag %d skipMargin %d Margins %d %d %d %d BMASK %x solution %d\n",
    __DATE__,
    __TIME__,
    skipFlag,
    updateBlendFlag,
    skipMargin,
    leftMargin,
    rightMargin,
    bottomMargin,
    topMargin,
    FILTER_BMASK_SEXTRACTOR,
    solutionNumber
  );

  fprintf(
    outHandle,
    "NUMBER\t"
    "FLUX_ISO\t"
    "MAG_ISO\t"
    "MAGERR_ISO\t"
    "MAG_APER\t"
    "MAGERR_APER\t"
    "MAG_AUTO\t"
    "MAGERR_AUTO\t"
    "KRON_RADIUS\t"
    "BACKGROUND\t"
    "THRESHOLD\t"
    "FLUX_MAX\t"
    "ISOAREA_IMAGE\t"
    "ISOAREA_WORLD\t"
    "X_IMAGE\t"
    "Y_IMAGE\t"
    "ra\t"
    "dec\t"
    "THETA_J2000\t"
    "ELLIPTICITY\t"
    "ERRTHETA_J2000\t"
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
    "AFLAGS\t"
    "BFLAGS\t"
    "plate_dra\t"
    "plate_ddec\t"
    "plate_dist\t"
    "aLength\t"
    "bLength\t"
    "THETA_IMAGE\t"
    "ERRA_IMAGE\t"
    "ERRB_IMAGE\t"
    "ERRTHETA_IMAGE\t"
    "FLUXERR_ISO\n"
  );
  fprintf(
    outHandle,
    "------\t"
    "--------\t"
    "-------\t"
    "----------\t"
    "--------\t"
    "-----------\t"
    "--------\t"
    "-----------\t"
    "-----------\t"
    "----------\t"
    "---------\t"
    "--------\t"
    "-------------\t"
    "-------------\t"
    "-------\t"
    "-------\t"
    "--\t"
    "---\t"
    "-----------\t"
    "-----------\t"
    "--------------\t"
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
    "------\t"
    "------\t"
    "---------\t"
    "----------\t"
    "----------\t"
    "-------\t"
    "-------\t"
    "-----------\t"
    "----------\t"
    "----------\t"
    "--------------\t"
    "-----------\n"
  );

  /* If the filter_multiple "none" file exists, read it in now */
  if (none_name[0] != 0) {
    none_header = table_header(none_handle,TABLE_PARSE);
    if (none_header == NULL) {
      fprintf(stderr,"ERROR: Failed to read header for %s\n",none_name);
      return 1;
    }

    none_table = table_loadva(
      none_handle,
      &none_header,
      NULL, /* hbase */
      NULL, /* rows */
      NULL,
      sizeof(NONE),
      &none_nrecs,
      TblInt, "NUMBER", TblOff(PNONE, NUMBER),
      0, "end", 0
    );

    if (none_table == NULL) {
      fprintf(stderr,"ERROR: Failed to read table for %s\n",none_name);
      return 1;
    }

    if (verbose) {
      fprintf(stderr,"read %zu records for %s\n",none_nrecs,none_name);
    }

    if (none_nrecs < 2) {
      fprintf(stderr,"ERROR: only %zu records read from %s\n",none_nrecs,none_name);
      return 1;
    }

    /* Sort the NONE table in order of increasing Sextractor number */
    qsort((void *) none_table, none_nrecs, sizeof(NONE), NoneCompare);
  }

  none_index = 0;
  pNone1 = &none_table[none_index];
  pNone2 = &none_table[none_index+1];

  if (updateBlendFlag) {
    printf("Update blend ...\n");

    short_header = table_header(short_handle,TABLE_PARSE);
    if (short_header == NULL) {
      fprintf(stderr,"ERROR: Failed to read header for %s\n",short_name);
      return 1;
    }

    short_descriptor = table_create_descrip(
      &short_nrecs,
      TblInt,"NUMBER", TblOff(PSHORT, NUMBER),
      TblDbl,"THRESHOLD", TblOff(PSHORT, THRESHOLD),
      TblDbl,"FLUX_MAX", TblOff(PSHORT, FLUX_MAX),
      TblDbl,"X_IMAGE", TblOff(PSHORT, X_IMAGE),
      TblDbl,"Y_IMAGE", TblOff(PSHORT, Y_IMAGE),
      TblDbl,"THETA_IMAGE", TblOff(PSHORT, THETA_IMAGE),
      TblDbl,"ELLIPTICITY", TblOff(PSHORT, ELLIPTICITY),
      TblInt,"ISO0", TblOff(PSHORT, ISO[0]),
      TblInt,"ISO1", TblOff(PSHORT, ISO[1]),
      TblInt,"ISO2", TblOff(PSHORT, ISO[2]),
      TblInt,"ISO3", TblOff(PSHORT, ISO[3]),
      TblInt,"ISO4", TblOff(PSHORT, ISO[4]),
      TblInt,"ISO5", TblOff(PSHORT, ISO[5]),
      TblInt,"ISO6", TblOff(PSHORT, ISO[6]),
      TblInt,"ISO7", TblOff(PSHORT, ISO[7]),
      TblInt,"BFLAGS", TblOff(PSHORT, BFLAGS),
      0, "end", 0
    );

    if (short_descriptor == NULL) {
      fprintf(stderr, "ERROR: Failed to allocate descriptor for %s\n", short_name);
      return 1;
    }

    table_loadmap(short_header, short_descriptor);
    short_nrecs = 0;

    while (1) {
      memset(pShort, 0, sizeof(SHORT));
      if (ReadShort(pShort, short_handle, short_header, short_descriptor, &short_row) == 0) {
        break;
      }

      if (pShort->NUMBER > maxNUMBER) {
        maxNUMBER = pShort->NUMBER;
      }

      short_total_nrecs++;

      /* Keep only the objects flagged as blended */
      if (pShort->BFLAGS & FILTER_BFLAG_BLEND) {
        if (short_allocation < short_nrecs + 1) {
          short_allocation += SHORT_ALLOCATION_INCREMENT;
          tmp_short_table = realloc(short_table, short_allocation * sizeof(SHORT));

          if (tmp_short_table == NULL) {
            fprintf(stderr,"ERROR: Failed to reallocate tmp_short_table of size %d\n",short_allocation);
            exit(1);
          }

          short_table = tmp_short_table;
          tmp_short_table = NULL;
        }

        memcpy(&short_table[short_nrecs], pShort, sizeof(SHORT));
        short_nrecs++;
      }
    }

    if (verbose) {
      fprintf(stderr, "read %zu records for %s\n", short_nrecs, short_name);
    }

    if (short_total_nrecs > 0) {
      time_t t0, elapsed;
      kdtree_t *kdt;
      kdtree_qres_t *kdi = NULL;
      double *kd_data, kd_lim_lo[2], kd_lim_hi[2];
      double maxDist2;
      const int n_leaf = 32; // a recommended value
      const int n_dim = 2;
      unsigned long n_tests = 0;

      if (magisofile[0] != 0) {
        vector = (double *) calloc(short_total_nrecs, sizeof(double));
        if (vector == NULL) {
          printf("ERROR: failed to allocate vector of size %d\n", short_total_nrecs);
          exit(1);
        }
      }

      /* Sort this table in increasing Y_IMAGE */
      qsort((void *) short_table, short_nrecs, sizeof(SHORT), ImageCompare);

      printf("- pass 1 ...\n");
      kd_data = calloc(n_dim * short_nrecs, sizeof(double));

      for (short_index = 0; short_index < short_nrecs; short_index++) {
        pShort = &short_table[short_index];

        if (GetEllipseDimensions(
          pShort->FLUX_MAX,
          pShort->THRESHOLD,
          pShort->ELLIPTICITY,
          pShort->ISO,
          &pShort->aLength,
          &pShort->bLength
        )) {
          fprintf(stderr, "ERROR: Illegal ISO index for NUMBER %d of %s\n", pShort->NUMBER, short_name);
          exit(1);
        }

        if (pShort->aLength > maxLength) {
          maxLength = pShort->aLength;
        }

        if (pShort->bLength > maxLength) {
          maxLength = pShort->bLength;
        }

        kd_data[2 * short_index] = pShort->X_IMAGE;
        kd_data[2 * short_index + 1] = pShort->Y_IMAGE;

        if (short_index == 0) {
          kd_lim_lo[0] = pShort->X_IMAGE;
          kd_lim_lo[1] = pShort->Y_IMAGE;
          kd_lim_hi[0] = pShort->X_IMAGE;
          kd_lim_hi[1] = pShort->Y_IMAGE;
        } else {
          kd_lim_lo[0] = Min(kd_lim_lo[0], pShort->X_IMAGE);
          kd_lim_lo[1] = Min(kd_lim_lo[1], pShort->Y_IMAGE);
          kd_lim_hi[0] = Max(kd_lim_lo[0], pShort->X_IMAGE);
          kd_lim_hi[1] = Max(kd_lim_lo[1], pShort->Y_IMAGE);
        }
      }

      printf("- maxLength: %.1f\n", maxLength);
      printf("- building tree (%zu recs) ...\n", short_nrecs);
      kdt = kdtree_new(short_nrecs, n_dim, n_leaf);
      kdtree_set_limits(kdt, kd_lim_lo, kd_lim_hi);
      kdt = kdtree_build(kdt, kd_data, short_nrecs, n_dim, n_leaf, KDTT_DOUBLE, KD_BUILD_BBOX | KD_BUILD_SPLIT);
      maxDist2 = 4 * maxLength * maxLength;

      if (verbose) {
        fprintf(stderr, "Max Length %f, maxNUMBER %d %s\n", maxLength, maxNUMBER, short_name);
      }

      AFLAGS_ARRAY = (int *) calloc(maxNUMBER + 1, sizeof(int));
      if (AFLAGS_ARRAY == NULL) {
        fprintf(stderr,"ERROR: Can not allocate AFLAGS_ARRAY of size %d for %s\n",maxNUMBER+1,short_name);
        exit(1);
      }

      printf("- pass 2 ...\n");
      t0 = time(NULL);

      for (short_index = 0; short_index < short_nrecs; short_index++) {
        u32 j;

        // We could potentially search using a size customized to this
        // particular source, rather than determining the global `maxLength`
        // parameter. If we do that, note that we'll have to skip the `idx2 <=
        // short_index` check below.
        pShort = &short_table[short_index];
        kdi = kdtree_rangesearch_options_reuse(
          kdt,
          kdi,
          &(pShort->X_IMAGE),
          maxDist2,
          KD_OPTIONS_NO_RESIZE_RESULTS
        );

        for (j = 0; j < kdi->nres; j++) {
          int idx2 = kdi->inds[j];
          SHORT *pShort2;

          if (idx2 <= short_index)
            continue; // already did this pair

          pShort2 = &short_table[idx2];
          n_tests++;

          // Note that aLength and bLength here are set by GetEllipseDimensions,
          // so out handling of ELLIPTICITY ~ 1 will already be in effect.
          if (CheckBlend(
            pShort->X_IMAGE,
            pShort->Y_IMAGE,
            pShort->THETA_IMAGE,
            pShort->aLength,
            pShort->bLength,
            pShort2->X_IMAGE,
            pShort2->Y_IMAGE,
            pShort2->THETA_IMAGE,
            pShort2->aLength,
            pShort2->bLength
          )) {
            /* This is a blend! */
            if ((pShort->AFLAGS & (1 << FILTER_AFLAG_BLEND)) == 0) {
              pShort->AFLAGS |= (1 << FILTER_AFLAG_BLEND);
              AFLAGS_ARRAY[pShort->NUMBER] = pShort->AFLAGS;
              short_blend_nrecs++;
            }

            if ((pShort2->AFLAGS & (1 << FILTER_AFLAG_BLEND)) == 0) {
              pShort2->AFLAGS |= (1 << FILTER_AFLAG_BLEND);
              AFLAGS_ARRAY[pShort2->NUMBER] = pShort2->AFLAGS;
              short_blend_nrecs++;
            }
          }
        }
      }

      if (verbose) {
        fprintf(stderr, "revised blend records %d %s\n", short_blend_nrecs, short_name);
      }

      elapsed = time(NULL) - t0;
      printf(
        "- ... done: %lld s, %lu tests, %.0lf recs/s, %.0lf tests/s\n",
        (long long) elapsed,
        n_tests,
        (double) short_nrecs / elapsed,
        (double) n_tests / elapsed
      );
      kdtree_free_query(kdi);
      kdtree_free(kdt);
    } /* short_nrecs test */
  } /* update blend calculation */

  /* Free up memory */

  if (short_header != NULL) {
    table_hdrfree(short_header);
    short_header = NULL;
  }

  if (short_handle != NULL) {
    Close(short_handle);
    short_handle = NULL;
  }

  if (short_descriptor != NULL) {
    Free(short_descriptor);
    short_descriptor = NULL;
  }

  if (short_row != NULL) {
    table_rowfree(short_row);
    short_row = NULL;
  }

  if (short_table != NULL) {
    free(short_table);
    short_table = NULL;
  }

  /* Now read in the background file */
  if (background_handle != NULL) {
    background_header = table_header(background_handle, TABLE_PARSE);

    if (background_header == NULL) {
      fprintf(stderr, "ERROR: Failed to read header for %s\n", background_name);
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
      0,"end",0
    );

    if (background_table == NULL) {
      fprintf(stderr, "ERROR: Failed to read table for %s\n", background_name);
      return 1;
    }

    if (verbose) {
      fprintf(stderr,"read %zu records for %s\n",background_nrecs,background_name);
    }

    /* Now set the appropriate bits in the AFLAGS array */

    for (background_index = 0; background_index < background_nrecs; background_index++) {
      pBackground = &background_table[background_index];

      if (pBackground->NUMBER <= maxNUMBER) {
        AFLAGS_ARRAY[pBackground->NUMBER] |= pBackground->AFLAGS;
      } else {
        printf("ERROR: background NUMBER %d exceeds %d\n",pBackground->NUMBER,maxNUMBER);
        exit(1);
      }
    }

    /* We are done with the background table.  Clean up memory */

    if (background_table != NULL) {
      free(background_table);
      background_table = NULL;
    }

    if (background_header != NULL) {
      table_hdrfree(background_header);
      background_header = NULL;
    }

    if (background_handle != NULL) {
      Close(background_handle);
      background_handle = NULL;
    }

    pBackground = NULL;
  }

  /* Now read in the sextractor file */

  sextractor_header = table_header(sextractor_handle, TABLE_PARSE);
  if (sextractor_header == NULL) {
    fprintf(stderr, "ERROR: Failed to read header for %s\n", sextractor_name);
    return 1;
  }

  sextractor_descriptor = table_create_descrip(
    &sextractor_nrecs,
    TblInt, "NUMBER", TblOff(PSEXTRACTOR, NUMBER),
    TblDbl, "FLUX_ISO", TblOff(PSEXTRACTOR, FLUX_ISO),
    TblDbl, "MAG_ISO", TblOff(PSEXTRACTOR, MAG_ISO),
    TblDbl, "MAGERR_ISO", TblOff(PSEXTRACTOR, MAGERR_ISO),
    TblDbl, "MAG_APER", TblOff(PSEXTRACTOR, MAG_APER),
    TblDbl, "MAGERR_APER", TblOff(PSEXTRACTOR, MAGERR_APER),
    TblDbl, "MAG_AUTO", TblOff(PSEXTRACTOR, MAG_AUTO),
    TblDbl, "MAGERR_AUTO", TblOff(PSEXTRACTOR, MAGERR_AUTO),
    TblDbl, "KRON_RADIUS", TblOff(PSEXTRACTOR, KRON_RADIUS),
    TblDbl, "BACKGROUND", TblOff(PSEXTRACTOR, BACKGROUND),
    TblDbl, "THRESHOLD", TblOff(PSEXTRACTOR, THRESHOLD),
    TblDbl, "FLUX_MAX", TblOff(PSEXTRACTOR, FLUX_MAX),
    TblInt, "ISOAREA_IMAGE", TblOff(PSEXTRACTOR, ISOAREA_IMAGE),
    TblDbl, "ISOAREA_WORLD", TblOff(PSEXTRACTOR, ISOAREA_WORLD),
    TblDbl, "X_IMAGE", TblOff(PSEXTRACTOR, X_IMAGE),
    TblDbl, "Y_IMAGE", TblOff(PSEXTRACTOR, Y_IMAGE),
    TblDbl, "ra", TblOff(PSEXTRACTOR, ra),
    TblDbl, "dec", TblOff(PSEXTRACTOR, dec),
    TblDbl, "THETA_IMAGE", TblOff(PSEXTRACTOR, THETA_IMAGE),
    TblDbl, "THETA_J2000", TblOff(PSEXTRACTOR, THETA_J2000),
    TblDbl, "ELLIPTICITY", TblOff(PSEXTRACTOR, ELLIPTICITY),
    TblDbl, "ERRTHETA_J2000", TblOff(PSEXTRACTOR, ERRTHETA_J2000),
    TblDbl, "FWHM_IMAGE", TblOff(PSEXTRACTOR, FWHM_IMAGE),
    TblDbl, "FWHM_WORLD", TblOff(PSEXTRACTOR, FWHM_WORLD),
    TblInt, "ISO0", TblOff(PSEXTRACTOR, ISO[0]),
    TblInt, "ISO1", TblOff(PSEXTRACTOR, ISO[1]),
    TblInt, "ISO2", TblOff(PSEXTRACTOR, ISO[2]),
    TblInt, "ISO3", TblOff(PSEXTRACTOR, ISO[3]),
    TblInt, "ISO4", TblOff(PSEXTRACTOR, ISO[4]),
    TblInt, "ISO5", TblOff(PSEXTRACTOR, ISO[5]),
    TblInt, "ISO6", TblOff(PSEXTRACTOR, ISO[6]),
    TblInt, "ISO7", TblOff(PSEXTRACTOR, ISO[7]),
    TblInt, "BFLAGS", TblOff(PSEXTRACTOR, BFLAGS),
    TblDbl, "plate_dra", TblOff(PSEXTRACTOR, plate_dra),
    TblDbl, "plate_ddec", TblOff(PSEXTRACTOR, plate_ddec),
    TblDbl, "plate_dist", TblOff(PSEXTRACTOR, plate_dist),
    TblDbl, "ERRA_IMAGE", TblOff(PSEXTRACTOR, ERRA_IMAGE),
    TblDbl, "ERRB_IMAGE", TblOff(PSEXTRACTOR, ERRB_IMAGE),
    TblDbl, "ERRTHETA_IMAGE", TblOff(PSEXTRACTOR, ERRTHETA_IMAGE),
    TblDbl, "FLUXERR_ISO", TblOff(PSEXTRACTOR, FLUXERR_ISO),
    0,"end", 0
  );

  if (sextractor_descriptor == NULL) {
    fprintf(stderr,"ERROR: Failed to allocate descriptor for %s\n",sextractor_name);
    return 1;
  }

  table_loadmap(sextractor_header, sextractor_descriptor);

  while (1) {
    double aLength, bLength;

    memset(pSextractor, 0, sizeof(SEXTRACTOR));

    if(
      ReadSextractor(
        pSextractor,
        sextractor_handle,
        sextractor_header,
        sextractor_descriptor,
        &sextractor_row
      ) == 0
    ) {
      break;
    }

    sextractor_nrecs++;

    if (skipMargin == 0) {
      if (
        pSextractor->X_IMAGE < leftMarginD ||
        pSextractor->X_IMAGE > rightMarginD ||
        pSextractor->Y_IMAGE < bottomMarginD ||
        pSextractor->Y_IMAGE > topMarginD
      ) {
        marginCount++;
        continue;
      }
    }

    /* Reject this entry if NUMBER is not in the NONE table */

    if (none_table != NULL) {
      while (pSextractor->NUMBER < pNone1->NUMBER) {
        if (none_index <= 0) {
          break;
        }

        none_index--;
        pNone1 = &none_table[none_index];
        pNone2 = &none_table[none_index + 1];
      }

      while (pSextractor->NUMBER >= pNone2->NUMBER) {
        if (none_index + 1 >= none_nrecs) {
          break;
        }

        none_index++;
        pNone1 = &none_table[none_index];
        pNone2 = &none_table[none_index + 1];
      }

      if (pSextractor->NUMBER != pNone1->NUMBER) {
        noneCount++;
        continue;
      }
    }

    if (skipFlag == 0) {
      pix2wcs(
        wcs,
        pSextractor->X_IMAGE,
        pSextractor->Y_IMAGE,
        &pSextractor->ra,
        &pSextractor->dec
      );

      // Update derived values. Equations here mirror the ones in the pipeline
      // scripts that compute the original values.
      //
      // When this tool updates a SExtractor derived from a WW mosaic with WCS
      // from a TNX mosaic, the following columns should be updated but aren't:
      // ISOAREA_WORLD, THETA_J2000, ERRTHETA_J2000, FWHM_WORLD. ("Should" in
      // the sense that if you ran SExtractor on the TNX file to start, you
      // would get different results.) Out of those columns, THETA_J2000 and
      // FWHM_WORLD are used in the pipeline analysis.

      pSextractor->plate_ddec = pSextractor->dec - cdec;
      pSextractor->plate_dra = (pSextractor->ra - cra) * cos(((pSextractor->dec + cdec) / 2) / 57.29577951);
      pSextractor->plate_dist = sqrt(
        pSextractor->plate_dra * pSextractor->plate_dra + pSextractor->plate_ddec * pSextractor->plate_ddec
      );
    }

    if (
      GetEllipseDimensions(
        pSextractor->FLUX_MAX,
        pSextractor->THRESHOLD,
        pSextractor->ELLIPTICITY,
        pSextractor->ISO,
        &aLength,
        &bLength
      )
    ) {
      fprintf(stderr, "ERROR: Illegal ISO index for NUMBER %d of %s\n", pShort->NUMBER, short_name);
      exit(1);
    }

    if (pSextractor->NUMBER <= maxNUMBER) {
      pSextractor->AFLAGS = AFLAGS_ARRAY[pSextractor->NUMBER];
    }

    fprintf(
      outHandle,
      "%d\t"
      "%12.7g\t"
      "%.4f\t"
      "%.4f\t"
      "%.4f\t"
      "%.4f\t"
      "%.4f\t"
      "%.4f\t"
      "%.2f\t"
      "%12.7g\t"
      "%12.7g\t"
      "%12.7g\t"
      "%d\t"
      "%12.7g\t"
      "%.3f\t"
      "%.3f\t"
      "%.7f\t"
      "%.7f\t"
      "%.2f\t"
      "%.3f\t"
      "%.1f\t"
      "%.2f\t"
      "%12.7g\t"
      "%d\t"
      "%d\t"
      "%d\t"
      "%d\t"
      "%d\t"
      "%d\t"
      "%d\t"
      "%d\t"
      "%d\t"
      "%d\t"
      "%.6f\t"
      "%.4f\t"
      "%.4f\t"
      "%.4f\t"
      "%.4f\t"
      "%f\t"
      "%.4f\t"
      "%.4f\t"
      "%.1f\t"
      "%.7f\n",
      pSextractor->NUMBER,
      pSextractor->FLUX_ISO,
      pSextractor->MAG_ISO,
      pSextractor->MAGERR_ISO,
      pSextractor->MAG_APER,
      pSextractor->MAGERR_APER,
      pSextractor->MAG_AUTO,
      pSextractor->MAGERR_AUTO,
      pSextractor->KRON_RADIUS,
      pSextractor->BACKGROUND,
      pSextractor->THRESHOLD,
      pSextractor->FLUX_MAX,
      pSextractor->ISOAREA_IMAGE,
      pSextractor->ISOAREA_WORLD,
      pSextractor->X_IMAGE,
      pSextractor->Y_IMAGE,
      pSextractor->ra,
      pSextractor->dec,
      pSextractor->THETA_J2000,
      pSextractor->ELLIPTICITY,
      pSextractor->ERRTHETA_J2000,
      pSextractor->FWHM_IMAGE,
      pSextractor->FWHM_WORLD,
      pSextractor->ISO[0],
      pSextractor->ISO[1],
      pSextractor->ISO[2],
      pSextractor->ISO[3],
      pSextractor->ISO[4],
      pSextractor->ISO[5],
      pSextractor->ISO[6],
      pSextractor->ISO[7],
      pSextractor->AFLAGS,
      pSextractor->BFLAGS,
      pSextractor->plate_dra,
      pSextractor->plate_ddec,
      pSextractor->plate_dist,
      aLength,
      bLength,
      pSextractor->THETA_IMAGE,
      pSextractor->ERRA_IMAGE,
      pSextractor->ERRB_IMAGE,
      pSextractor->ERRTHETA_IMAGE,
      pSextractor->FLUXERR_ISO
    );

    outputCount++;

    if (magisofile[0] != 0) {
      if (vector_nrecs >= short_total_nrecs) {
        printf("ERROR: vector_nrecs %d exceeds short_total_nrecs %d\n",vector_nrecs,short_total_nrecs);
        exit(1);
      }

      vector[vector_nrecs] = pSextractor->MAG_ISO;
      vector_nrecs++;
    }
  }

  if (magisofile[0] != 0) {
    curClipCount = vector_nrecs;
    curClipCount = CalcMedianAndRMS(
      curClipCount,
      10,
      vector,
      &MAG_ISO_med,
      &MAG_ISO_rms,
      1,
      CLIP_RMS_FACTOR,
      0
    );

    if (verbose) {
      printf("curClipCount %d MAG_ISO median %f, MAG_ISO rms %f\n", curClipCount, MAG_ISO_med, MAG_ISO_rms);
    }

    SaveMAG_ISO(magisofile, fileroot, MAG_ISO_med, MAG_ISO_rms);
  }

  if (verbose) {
    fprintf(stderr, "read %zu records for %s\n", sextractor_nrecs, sextractor_name);
  }

  /* We are done with the sextractor table.  Clear memory */

  if (sextractor_header != NULL) {
    table_hdrfree(sextractor_header);
    sextractor_header = NULL;
  }

  if (sextractor_handle != NULL) {
    Close(sextractor_handle);
    sextractor_handle = NULL;
  }

  if (sextractor_descriptor != NULL) {
    Free(sextractor_descriptor);
    sextractor_descriptor = NULL;
  }

  if (sextractor_row != NULL) {
    table_rowfree(sextractor_row);
    sextractor_row = NULL;
  }

  /* Done with the none table */

  if (none_table != NULL) {
    Free(none_table);
    none_table = NULL;
  }

  if (none_header != NULL) {
    table_hdrfree(none_header);
    none_header = NULL;
  }

  if (none_handle != NULL) {
    Close(none_handle);
    none_handle = NULL;
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }

  if (background_handle != NULL) {
    fclose(background_handle);
  }

  if (wcs != NULL) {
    wcsfree (wcs);
    wcs = NULL;
  }

  if (AFLAGS_ARRAY != NULL) {
    free(AFLAGS_ARRAY);
  }

  if (vector != NULL) {
    free(vector);
    vector = NULL;
  }

  time(&curTime);
  curTime -= startTime;

  fprintf(
    stdout,
    "update_sextractor starsin %zu starsout %d margin stars %d blend in %zu blend out %d noneRejects %d seconds %ld for %s\n",
    sextractor_nrecs,
    outputCount,
    marginCount,
    short_nrecs,
    short_blend_nrecs,
    noneCount,
    curTime,
    outfile
  );

  return 0;
}

/*
 * Sep  1, 2008 Edward J. Los - Initial version
 * Nov  4, 2008 Edward J. Los - Add filtering of mosaic edge patterns
 * Jan 27, 2009 Edward J. Los - Add THETA_IMAGE Convert the sextractor FLAGS parameter to BFLAGS
 * Feb  6, 2009 Edward J. Los - Add -u to recalculate the blend flag
 * Sep 25, 2009 Edward J. Los - Write THETA_IMAGE to the output file
 * Oct  5, 2009 Edward J. Los - Update for multiple exposures
 * Sep 20, 2010 Edward J. Los - Add new sextractractor parameters for scamp support
 * Sep 29, 2010 Edward J. Los - Add a check to ensure that ctype1 is RA
 * Oct 26, 2010 Edward J. Los - Add InaccuratePV flag
 * Nov 24, 2010 Edward J. Los - Check mosaic DATE-OBS against final logbook date.
 * Sep 28, 2011 Edward J. Los - Add FILTER_BFLAG_PSFSATURATED to BFLAGS
 * Apr 15, 2012 Edward J. Los - add option to skip the update of the InaccuratePV flag
 * Apr 24, 2012 Edward J. Los - add option to output the median and rms clipped MAG_ISO
 * May 14, 2012 Edward J. Los - get the fileroot from magisofile if available
 * Jul  3, 2012 Edward J. Los - Add high background object support
 * Jul 18, 2012 Edward J. Los - Change GetEllipseDimensions to use the average of ISO3 and ISO4 because a reading of the SExtractor sources
 *                              shows that the ISO array has equal flux intervals, not magnitude intervals.
 * Jun 30, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
 */
