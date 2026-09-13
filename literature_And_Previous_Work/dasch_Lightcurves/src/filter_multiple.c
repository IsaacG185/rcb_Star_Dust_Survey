// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// This routine takes the source extractor file and any number of _allobjects
// files and outputs a file without matched objects for astrometry.net
// processing.
//
// Usage:
//   filter_multiple
//     -e <solutionNumber, >= 1>
//     -q <catalog qualifier, ignored>
//     -s <input sextractor file>
//     -a <combined allobjects output file>
//     -o <output file for unmatched ("none") objects>
//     <previous allobjects files ...>

#include <math.h>
#include <time.h>

#include <table.h>

#include <mysql.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "daschunistd.h"

#define MAX_BUFFER 512

#define STAR_SELECTED  1  /* Image is present in the sextractor file */
#define NONE_STAR      2  /* Image found to be unmatched */
#define REF_STAR       4  /* Image matched with at least one catalog star */
#define MULTIPLE_STAR  8  /* Image matched with multiple catalog stars */
#define MULTIPLE_NONE 16  /* Multiple unmatched solutions for this object */

#define SEXTRACTOR_ALLOCATION 10000

typedef struct _allobjects {
  char REF[MAX_REF];
  double X_IMAGE;
  double Y_IMAGE;
  double MAG_ISO;
  double ra;
  double dec;
  double magcal_iso;
  double magcal_iso_rms;
  double magcal_local;
  double magcal_local_error;
  double magcal_local_rms;
  double Date;
  double limiting_mag_local;
  double extinction;
  double Stdmag;
  double color;
  double dra;
  double ddec;
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
  double plate_dist;
  double Blendedmag;
  double dradRMS2;
  double magcal_magdep;  /* magnitude-corrected local magnitude */
  double magcal_magdep_rms; /* RMS of the magnitude-corrected local magnitude */
  int NUMBER;
  int AFLAGS;
  int BFLAGS;
  int npoints_local;
  int spatial_bin;
  int ISO0;
  int ISO1;
  int ISO2;
  int ISO3;
  int ISO4;
  int ISO5;
  int ISO6;
  int ISO7;
  int gsc_bin_index;
  int local_bin;
  int solutionNumber;
  int exposureNumber;
  int nextAllobject;
  int magdep_bin;     /* Magnituded-dependent calibration bin (-1) if undefined */
} ALLOBJECTS, *PALLOBJECTS;

typedef struct _starimage {
  struct _allobjects allObjects;
  int NUMBER;
  int BFLAGS;
  double MAG_ISO;
  double X_IMAGE;
  double Y_IMAGE;
  double ELLIPTICITY;
  int spatial_bin;
  int selectFlag;
  int selectCount;
} SEXTRACTOR, *PSEXTRACTOR;


static int
ReadSextractor(
  PSEXTRACTOR pSextractor,
  File sextractor_handle,
  TableHead sextractor_header,
  TblDescriptor sextractor_descriptor,
  TableRow* sextractor_row
) {
  *sextractor_row = table_rowget(sextractor_handle,sextractor_header,*sextractor_row,NULL,NULL,0);
  if (*sextractor_row == NULL) {
    return 0;
  }

  if (!table_loadrow(sextractor_handle,sextractor_header,*sextractor_row,sextractor_descriptor,(char *)pSextractor)) {
    printf("ERROR: filter_multiple  Read Sextractor table_loadrow failed\n");
    return 0;
  }

  memset(&pSextractor->allObjects,0,sizeof(ALLOBJECTS));
  pSextractor->selectFlag = 0;
  pSextractor->selectCount = 0;
  pSextractor->spatial_bin = 0;
  return 1;
}


static int
ReadAllobjects(
  PALLOBJECTS pAllobjects,
  File allobjects_handle,
  TableHead allobjects_header,
  TblDescriptor allobjects_descriptor,
  TableRow* allobjects_row
) {
  *allobjects_row = table_rowget(allobjects_handle, allobjects_header, *allobjects_row, NULL, NULL, 0);
  if (*allobjects_row == NULL) {
    return 0;
  }

  if (!table_loadrow(allobjects_handle, allobjects_header, *allobjects_row, allobjects_descriptor, (char *) pAllobjects)) {
    printf("ERROR: filter_multiple  Read Allobjects table_loadrow failed\n");
    return 0;
  }

  pAllobjects->solutionNumber = 0;
  pAllobjects->exposureNumber = 0;
  pAllobjects->nextAllobject = 0;
  return 1;
}


int
main(int argc, char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  char outfile[MAX_BUFFER];
  char debugname[MAX_BUFFER];
  char combinedfile[MAX_BUFFER];
  int numFiles = 0;
  int fileIndex;
  char fileroot[MAX_BUFFER];

  FILE *outHandle = NULL;
  FILE *debugHandle = NULL;
  FILE *combinedHandle = NULL;
  int verbose = 0;
  char *charPtr;
  int pixelBinning = 0;

  int spatial_bin;
  int acceptBin[MAX_SPATIAL_BINS+1];
  int binCount[MAX_SPATIAL_BINS+1];

  time_t startTime;
  time_t curTime;

  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  PSEXTRACTOR sextractor_table = NULL;
  PSEXTRACTOR tmp_sextractor_table;
  size_t sextractor_nrecs = 0;   /* Number of sextractor images actually read */
  int sextractor_totrecs = 0; /* Additional images for multiple matches */
  int sextractor_index;
  int sextractor_alloc = 0;
  int sextractor_new_alloc;
  int sextractorMaxNUMBER = 0;
  SEXTRACTOR sextractor_record;
  PSEXTRACTOR pSextractor = &sextractor_record;
  PSEXTRACTOR pSextractor2;
  PSEXTRACTOR pSextractor3;
  PSEXTRACTOR pSextractorSel;
  TblDescriptor sextractor_descriptor = NULL;
  TableRow sextractor_row = NULL;

  File allobjects_handle[MAX_CATALOG_EXPOSURES];
  char allobjects_name[MAX_CATALOG_EXPOSURES][MAX_BUFFER];
  TableHead allobjects_header = NULL;
  size_t allobjects_nrecs = 0;
  ALLOBJECTS allobjects_record;
  PALLOBJECTS pAllobjects = &allobjects_record;
  PALLOBJECTS pSavedAllobjects = NULL;
  TblDescriptor allobjects_descriptor = NULL;
  TableRow allobjects_row = NULL;
  PALLOBJECTS allobjects_table = NULL;
  int duplicateFlag = 0;

  int matchedObjectCount = 0;
  int outputCount = 0;
  int duplicateREFCount = 0;
  int multipleNoneCount = 0;
  int initialNoneCount = 0;
  int replacedNoneCount = 0;
  int initialREFCount = 0;
  int multipleREFCount = 0;
  int galaxyCount = 0;
  int combinedCount = 0;
  int multipleBlendGalaxyCount = 0;
  int multipleBlendCount = 0;

  int exposureNumber;
  double singleJulianDate;
  double singleTolerance;
  double allJulianDate;
  double allTolerance;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int rotation;
  int binning;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  int numExposures;

  int galaxyFlag = 0;
  int nextAllobject;
  PALLOBJECTS pAllobjectsSel;
  double tmpFlux;
  double blendedFlux = 0.0;
  double blendedMag;
  int curNoneMultiple = 0;
  int maxNoneMultiple = 0;
  int maxNoneNUMBER = 0;
  int solutionNumber = 0;
  int readFlag = 0;
  int FitWCS;

  time(&startTime);

  sextractor_name[0] = 0;
  debugname[0] = 0;
  outfile[0] = 0;
  combinedfile[0] = 0;

  for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
    allobjects_handle[fileIndex] = NULL;
  }

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    acceptBin[spatial_bin] = 1;
  }

  // Process arguments

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    if (argstr[0] != '-') {
      /* This must be a filename */
      strncpy(allobjects_name[numFiles],argstr,MAX_BUFFER-1);

      if (numFiles >= MAX_CATALOG_EXPOSURES) {
        fprintf(stdout,"ERROR: filter_multiple  Too many files, maximum is %d\n",MAX_CATALOG_EXPOSURES);
        return 1;
      }

      numFiles++;
    } else {
      while ((cmdchar = *++argstr) != '\0') {
        switch(cmdchar) {

        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: filter_multiple  Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: filter_multiple  MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'a': /* combined output file name */
        case 'A':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: filter_multiple  Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(combinedfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: filter_multiple  MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'e': /* solution number */
        case 'E':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&solutionNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the solutionNumber %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'q': /* file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          }

          // We don't actually do anything with this argument anymore, but preserve
          // it for compatibility
          argv++;
          break;

        case 's': /* sextractor file name */
        case 'S':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: filter_multiple  Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sextractor_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: filter_multiple  MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'b': /* pixel binning */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: filter_multiple  Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&pixelBinning);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: filter_multiple  Unable to decode pixel binning %s\n",*argv);
              errorFlag = 1;
            } else {
              if (pixelBinning < 1) {
                fprintf(stderr,"ERROR: filter_multiple  Invalid pixel binning value %d\n",pixelBinning);
                errorFlag = 1;
              }
            }
          }
          break;

        case 'd': /* debug file name */
        case 'D':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: filter_multiple  Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(debugname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: filter_multiple  MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        default:
          fprintf(stderr,"ERROR: filter_multiple   unknown command -%c\n",cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  // Verify that we have everything

  if (outfile[0] == 0) {
    fprintf(stderr,"ERROR: filter_multiple  No output filename was specified\n");
    errorFlag = 1;
  }

  if (combinedfile[0] == 0) {
    fprintf(stderr,"ERROR: filter_multiple  No combined allobjects output filename was specified\n");
    errorFlag = 1;
  }

  if (sextractor_name[0] == 0) {
    fprintf(stderr,"ERROR: filter_multiple  No sextractor filename was specified\n");
    errorFlag = 1;
  }

  if (numFiles == 0) {
    fprintf(stderr,"ERROR: filter_multiple  No allobjects files specified\n");
    errorFlag = 1;
  }

  if (pixelBinning == 0) {
    fprintf(stderr,"ERROR: filter_multiple  No pixel binning specified\n");
    errorFlag = 1;
  }

  charPtr = strrchr(sextractor_name, '/');
  if (charPtr == NULL) {
    strcpy(fileroot, sextractor_name);
  } else {
    charPtr++;
    strcpy(fileroot, charPtr);
  }

  if ((charPtr = strstr(fileroot, "_tnx")) != NULL) {
    *charPtr = 0;
  } else if ((charPtr = strstr(fileroot, ".")) != NULL) {
    *charPtr = 0;
  }

  if (debugname[0] != 0) {
    debugHandle = fopen(debugname, "wt");
    if (debugHandle == NULL) {
      printf("ERROR: failed to open %s\n",debugname);
    }
  }

  sextractor_handle = Open(sextractor_name, "r");
  if (sextractor_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr, "ERROR: filter_multiple  Failed to find the sextractor file %s\n", sextractor_name);
  } else if (verbose) {
    fprintf(stderr, "Found sextractor file %s\n", sextractor_name);
  }

  if (ParseFilename(sextractor_name, series, &plateNumber, &mosaicNumber, &binning, &rotation) == 0) {
    fprintf(
      stderr,
      "ERROR Failed to parse filename %s as %s%05d_%02d_%02d rotation %d\n",
      sextractor_name,
      series,
      plateNumber,
      mosaicNumber,
      binning,
      rotation
    );
    errorFlag = 1;
  }

  for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
    allobjects_handle[fileIndex] = Open(allobjects_name[fileIndex], "r");

    if (allobjects_handle[fileIndex] == NULL) {
      errorFlag = 1;
      fprintf(stderr, "ERRORD: filter_multiple  Failed to find the allobjects file %s\n", allobjects_name[fileIndex]);
    } else if (verbose) {
      fprintf(stderr, "Found allobjects file %s\n", allobjects_name[fileIndex]);
    }
  }

  if (errorFlag) {
    fprintf(stderr, "Usage: filter_multiple <allobjects file 1> [<allobject file 2>...] options\n");
    fprintf(stderr, "                    -e <solutionNumber>\n");
    fprintf(stderr, "                    -s <sextractor file> \n");
    fprintf(stderr, "                    -o <output file> \n");
    fprintf(stderr, "                    -a <combined allobjects output file\n");
    fprintf(stderr, "                    -d  <debug file>\n");
    fprintf(stderr, "                    -q <catalog qualifier>  (no longer used)\n");
    fprintf(stderr, "                    -v  verbose\n");
    return 1;
  }

  // Ready to really get going here

  printf(
    "filter_multiple of %s %s writing %s %s\n",
    __DATE__,
    __TIME__,
    outfile,
    combinedfile
  );

  dasch_init_scandb(pConnection);

  // Read in the sextractor file

  sextractor_header = table_header(sextractor_handle, TABLE_PARSE);
  if (sextractor_header == NULL) {
    fprintf(stderr,"ERROR: filter_multiple  Failed to read header for %s\n", sextractor_name);
    return 1;
  }

  sextractor_descriptor = table_create_descrip(
    &sextractor_nrecs,
    TblInt, "NUMBER", TblOff(PSEXTRACTOR, NUMBER),
    TblDbl, "X_IMAGE", TblOff(PSEXTRACTOR, X_IMAGE),
    TblDbl, "Y_IMAGE", TblOff(PSEXTRACTOR, Y_IMAGE),
    TblInt, "BFLAGS", TblOff(PSEXTRACTOR, BFLAGS),
    TblDbl, "MAG_ISO", TblOff(PSEXTRACTOR, MAG_ISO),
    TblDbl, "ELLIPTICITY", TblOff(PSEXTRACTOR, ELLIPTICITY),
    0, "end", 0
  );

  if (sextractor_descriptor == NULL) {
    fprintf(stderr,"ERROR: filter_multiple  Failed to allocate descriptor for %s\n",sextractor_name);
    return 1;
  }

  table_loadmap(sextractor_header, sextractor_descriptor);

  while (1) {
    memset(pSextractor, 0, sizeof(SEXTRACTOR));

    if(ReadSextractor(pSextractor, sextractor_handle, sextractor_header, sextractor_descriptor, &sextractor_row) == 0) {
      break;
    }

    sextractor_nrecs++;

    // Insert this record in the correct position of the sextractor table

    if (pSextractor->NUMBER > sextractorMaxNUMBER) {
      sextractorMaxNUMBER = pSextractor->NUMBER;
    }

    if (pSextractor->NUMBER > sextractor_alloc - 1) {
      sextractor_new_alloc = (pSextractor->NUMBER + SEXTRACTOR_ALLOCATION + 1) / SEXTRACTOR_ALLOCATION;
      sextractor_new_alloc = sextractor_new_alloc * SEXTRACTOR_ALLOCATION;
      tmp_sextractor_table = realloc(sextractor_table, sextractor_new_alloc * sizeof(SEXTRACTOR));
      if (tmp_sextractor_table == NULL) {
        fprintf(stderr,"ERROR: filter_multiple  Failed to reallocate sextractor table of size %d\n",sextractor_new_alloc);
        exit(-1);
      }

      sextractor_table = tmp_sextractor_table;
      tmp_sextractor_table = NULL;

      for (sextractor_index = sextractor_alloc; sextractor_index <  sextractor_new_alloc; sextractor_index++) {
        pSextractor2 = &sextractor_table[sextractor_index];
        pSextractor2->selectFlag = 0;
        pSextractor2->selectCount = 0;
        memset(&pSextractor2->allObjects,0,sizeof(ALLOBJECTS));
      }

      sextractor_alloc = sextractor_new_alloc;
    }

    pSextractor->selectFlag = STAR_SELECTED;
    pSextractor2 = &sextractor_table[pSextractor->NUMBER];
    memcpy(pSextractor2, pSextractor, sizeof(SEXTRACTOR));
  }

  if (verbose) {
    fprintf(
      stderr,
      "read %zu records, maximum NUMBER %d for %s\n",
      sextractor_nrecs,
      sextractorMaxNUMBER,
      sextractor_name
    );
  }

  sextractor_totrecs = sextractorMaxNUMBER + 1; /* Start allocating for multiple REF matches for a single image */
  if (sextractor_nrecs <= 0) {
    fprintf(stderr, "ERROR: filter_multiple read no records from %s\n", sextractor_name);
    exit(-1);
  }

  // Now open and process each allobjects file in turn

  for (fileIndex = 0; fileIndex < numFiles; fileIndex++) {
    /* Get the exposure number */
    if (
      GetSolutionJulianDate(
        pConnection,
        series,
        plateNumber,
        mosaicNumber,
        fileIndex,
        &exposureNumber,
        &numExposures,
        &singleJulianDate,
        &singleTolerance,
        &allJulianDate,
        &allTolerance
      )
    ) {
      fprintf(
        stderr,
        "ERROR: filter_multiple  failed to get exposure information for %s s%d\n",
        fileroot,
        fileIndex
      );
      return 1;
    }

    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      binCount[spatial_bin] = 0;
    }

    allobjects_header = table_header(allobjects_handle[fileIndex], TABLE_PARSE);
    if (allobjects_header == NULL) {
      fprintf(stderr,"ERROR: filter_multiple  Failed to read header for %s\n",allobjects_name[fileIndex]);
      return 1;
    }

    allobjects_descriptor = table_create_descrip(
      &allobjects_nrecs,
      TblBuf,"REF"    ,TblOff(PALLOBJECTS,REF),MAX_REF,
      TblInt,"NUMBER",TblOff(PALLOBJECTS,NUMBER),
      TblInt,"AFLAGS",TblOff(PALLOBJECTS,AFLAGS),
      TblInt,"BFLAGS",TblOff(PALLOBJECTS,BFLAGS),
      TblInt,"npoints_local",TblOff(PALLOBJECTS,npoints_local),
      TblInt,"spatial_bin",TblOff(PALLOBJECTS,spatial_bin),
      TblInt,"ISO0",TblOff(PALLOBJECTS,ISO0),
      TblInt,"ISO1",TblOff(PALLOBJECTS,ISO1),
      TblInt,"ISO2",TblOff(PALLOBJECTS,ISO2),
      TblInt,"ISO3",TblOff(PALLOBJECTS,ISO3),
      TblInt,"ISO4",TblOff(PALLOBJECTS,ISO4),
      TblInt,"ISO5",TblOff(PALLOBJECTS,ISO5),
      TblInt,"ISO6",TblOff(PALLOBJECTS,ISO6),
      TblInt,"ISO7",TblOff(PALLOBJECTS,ISO7),
      TblInt,"gsc_bin_index",TblOff(PALLOBJECTS,gsc_bin_index),
      TblInt,"local_bin",TblOff(PALLOBJECTS,local_bin),
      TblInt,"magdep_bin",TblOff(PALLOBJECTS,magdep_bin),
      TblDbl,"X_IMAGE",TblOff(PALLOBJECTS,X_IMAGE),
      TblDbl,"Y_IMAGE",TblOff(PALLOBJECTS,Y_IMAGE),
      TblDbl,"MAG_ISO",TblOff(PALLOBJECTS,MAG_ISO),
      TblDbl,"ra",TblOff(PALLOBJECTS,ra),
      TblDbl,"dec",TblOff(PALLOBJECTS,dec),
      TblDbl,"magcal_iso",TblOff(PALLOBJECTS,magcal_iso),
      TblDbl,"magcal_iso_rms",TblOff(PALLOBJECTS,magcal_iso_rms),
      TblDbl,"magcal_local",TblOff(PALLOBJECTS,magcal_local),
      TblDbl,"magcal_local_error",TblOff(PALLOBJECTS,magcal_local_error),
      TblDbl,"magcal_local_rms",TblOff(PALLOBJECTS,magcal_local_rms),
      TblDbl,"Date",TblOff(PALLOBJECTS,Date),
      TblDbl,"limiting_mag_local",TblOff(PALLOBJECTS,limiting_mag_local),
      TblDbl,"extinction",TblOff(PALLOBJECTS,extinction),
      TblDbl,"Stdmag",TblOff(PALLOBJECTS,Stdmag),
      TblDbl,"color",TblOff(PALLOBJECTS,color),
      TblDbl,"dra",TblOff(PALLOBJECTS,dra),
      TblDbl,"ddec",TblOff(PALLOBJECTS,ddec),
      TblDbl,"FLUX_ISO",TblOff(PALLOBJECTS,FLUX_ISO),
      TblDbl,"MAG_APER",TblOff(PALLOBJECTS,MAG_APER),
      TblDbl,"MAG_AUTO",TblOff(PALLOBJECTS,MAG_AUTO),
      TblDbl,"KRON_RADIUS",TblOff(PALLOBJECTS,KRON_RADIUS),
      TblDbl,"BACKGROUND",TblOff(PALLOBJECTS,BACKGROUND),
      TblDbl,"THRESHOLD",TblOff(PALLOBJECTS,THRESHOLD),
      TblDbl,"FLUX_MAX",TblOff(PALLOBJECTS,FLUX_MAX),
      TblDbl,"THETA_J2000",TblOff(PALLOBJECTS,THETA_J2000),
      TblDbl,"ELLIPTICITY",TblOff(PALLOBJECTS,ELLIPTICITY),
      TblDbl,"ISOAREA_WORLD",TblOff(PALLOBJECTS,ISOAREA_WORLD),
      TblDbl,"FWHM_IMAGE",TblOff(PALLOBJECTS,FWHM_IMAGE),
      TblDbl,"FWHM_WORLD",TblOff(PALLOBJECTS,FWHM_WORLD),
      TblDbl,"plate_dist",TblOff(PALLOBJECTS,plate_dist),
      TblDbl,"Blendedmag",TblOff(PALLOBJECTS,Blendedmag),
      TblDbl,"dradRMS2",TblOff(PALLOBJECTS,dradRMS2),
      TblDbl,"magcal_magdep",TblOff(PALLOBJECTS,magcal_magdep),
      TblDbl,"magcal_magdep_rms",TblOff(PALLOBJECTS,magcal_magdep_rms),
      0,"end",0
    );

    if (allobjects_descriptor == NULL) {
      fprintf(stderr,"ERROR: filter_multiple  Failed to allocate descriptor for %s\n",allobjects_name[fileIndex]);
      return 1;
    }

    table_loadmap(allobjects_header,allobjects_descriptor);

    memset(pAllobjects,0,sizeof(ALLOBJECTS));
    allobjects_nrecs = 0;

    while (1) {
      if(
        ReadAllobjects(
          pAllobjects,
          allobjects_handle[fileIndex],
          allobjects_header,
          allobjects_descriptor,
          &allobjects_row
        ) == 0
      ) {
        break;
      }

      allobjects_nrecs++;

      if ((pAllobjects->spatial_bin > 0)  && (pAllobjects->spatial_bin <= MAX_SPATIAL_BINS)) {
        binCount[pAllobjects->spatial_bin]++;
      } else {
        fprintf(stderr,"ERROR: filter_multiple  Illegal spatial bin %d for NUMBER %d in %s\n",pAllobjects->spatial_bin,pAllobjects->NUMBER,allobjects_name[fileIndex]);
      }

      pAllobjects->solutionNumber = fileIndex; /* Save the solution of this file */
      pAllobjects->exposureNumber = exposureNumber;

      if (strcmp(pAllobjects->REF,"NONE") == 0) {
        /* Unmatched object, retain the selected flag and mark the fact that we processed it */
        if (pAllobjects->NUMBER <= sextractorMaxNUMBER) {
          pSextractor2 = &sextractor_table[pAllobjects->NUMBER];
          if (((pSextractor2->selectFlag & REF_STAR) != 0) ||
              ((pSextractor2->selectFlag & NONE_STAR) != 0)) {
            /* star either has been matched or already has an unmatched star.  Save this record because we need to write out all matches */
            if ((pSextractor2->selectFlag & NONE_STAR) != 0) {
              multipleNoneCount++;
              if (sextractor_totrecs >= (sextractor_alloc-1)) {
                /* Here we need to expand our array */
                sextractor_alloc  += SEXTRACTOR_ALLOCATION;
                tmp_sextractor_table = realloc(sextractor_table,sextractor_alloc * sizeof(SEXTRACTOR));
                if (tmp_sextractor_table == NULL) {
                  fprintf(stderr,"ERROR: filter_multiple  Failed to reallocate sextractor table of size %d\n",sextractor_alloc);
                  exit(-1);
                }
                sextractor_table = tmp_sextractor_table;
                tmp_sextractor_table = NULL;
                pSextractor2 = &sextractor_table[pAllobjects->NUMBER];

              }
              pSextractor2->selectFlag |= MULTIPLE_NONE;
              pSavedAllobjects = &pSextractor2->allObjects;
              pSavedAllobjects->nextAllobject = sextractor_totrecs;
              pSextractor3 =  &sextractor_table[sextractor_totrecs];
              memset(pSextractor3,0,sizeof(SEXTRACTOR));
              pSextractor3->NUMBER = sextractor_totrecs;
              pSextractor3->selectFlag |= NONE_STAR | MULTIPLE_NONE;
              pSavedAllobjects = &pSextractor3->allObjects;
              memcpy(pSavedAllobjects,pAllobjects,sizeof(ALLOBJECTS));
              sextractor_totrecs++;
            }
            continue;
          }
          /* Save the first unmatched image */
          pSavedAllobjects = &pSextractor2->allObjects;
          pSextractor2->selectFlag |= NONE_STAR;
          pSextractor2->spatial_bin = pAllobjects->spatial_bin;
          memcpy(pSavedAllobjects,pAllobjects,sizeof(ALLOBJECTS));
          initialNoneCount++;
          continue;

        }
        continue;
      }
      matchedObjectCount++;
      if (pAllobjects->NUMBER > sextractorMaxNUMBER) {
        /* Fake object for FILTER_AFLAG_BLEND_NOMATCH objects */
        continue;
      }


      /* This object has been matched to a GSC object. */
      pSextractor2 = &sextractor_table[pAllobjects->NUMBER];

      if (pSextractor2->selectFlag == 0) {
        /* Not in sextractor file */
        fprintf(stderr,"ERROR, allobject file %s does not match sextractor file %s\n",allobjects_name[fileIndex],sextractor_name);
        exit(-1);
      }

      if ((pSextractor2->selectFlag & NONE_STAR) !=  0) {
        replacedNoneCount++;
        /* Image has been matched, substitute its entry for the "NONE" entry and deselect it */
        pSavedAllobjects = &pSextractor2->allObjects;
        pSextractor2->selectFlag |= REF_STAR;
        pSextractor2->selectFlag &= ~(NONE_STAR|MULTIPLE_NONE);
        pSextractor2->spatial_bin = pAllobjects->spatial_bin;
        memcpy(pSavedAllobjects,pAllobjects,sizeof(ALLOBJECTS));
      } else if ((pSextractor2->selectFlag & REF_STAR) ==  0) {
        /* Image not yet matched, just store it */
        initialREFCount++;
        pSavedAllobjects = &pSextractor2->allObjects;
        pSextractor2->selectFlag |= REF_STAR;
        pSextractor2->spatial_bin = pAllobjects->spatial_bin;
        memcpy(pSavedAllobjects,pAllobjects,sizeof(ALLOBJECTS));

      } else if ((pSextractor2->selectFlag & REF_STAR) !=  0) {
        /* Reject this image if we have a duplicate REF.  This means that the two solution are really overlapping */
        duplicateFlag = 0;
        pSextractor3 = pSextractor2;
        while (1) {
          pSavedAllobjects = &pSextractor3->allObjects;
          if (strcmp(pSavedAllobjects->REF,pAllobjects->REF) == 0) {
            duplicateFlag = 1;
            break;
          }
          if (pSavedAllobjects->nextAllobject == 0) {
            break;
          }
          pSextractor3 = &sextractor_table[pSavedAllobjects->nextAllobject];
        }
        if (duplicateFlag) {
          duplicateREFCount++;
          continue;
        }

        multipleREFCount++;
        /* Here we have multiple exposure blend, save the results and link to it */
        if (sextractor_totrecs >= (sextractor_alloc-1)) {
          /* Here we need to expand our array */
          sextractor_alloc  += SEXTRACTOR_ALLOCATION;
          tmp_sextractor_table = realloc(sextractor_table,sextractor_alloc * sizeof(SEXTRACTOR));
          if (tmp_sextractor_table == NULL) {
            fprintf(stderr,"ERROR: filter_multiple  Failed to reallocate sextractor table of size %d\n",sextractor_alloc);
            exit(-1);
          }
          sextractor_table = tmp_sextractor_table;
          tmp_sextractor_table = NULL;
          pSextractor2 = &sextractor_table[pAllobjects->NUMBER];

        }
        pSextractor2->selectFlag |= MULTIPLE_STAR;
        pSavedAllobjects = &pSextractor2->allObjects;
        pSavedAllobjects->nextAllobject = sextractor_totrecs;
        pSextractor3 =  &sextractor_table[sextractor_totrecs];
        memset(pSextractor3,0,sizeof(SEXTRACTOR));
        pSextractor3->NUMBER = sextractor_totrecs;
        pSextractor3->selectFlag |= REF_STAR | MULTIPLE_STAR;
        pSavedAllobjects = &pSextractor3->allObjects;
        memcpy(pSavedAllobjects,pAllobjects,sizeof(ALLOBJECTS));
        sextractor_totrecs++;

      } else {
        /* Just ignore subsequent "NONE" images */
        fprintf(stderr,"ERROR: filter_multiple  NONE star in unexpected position for NUMBER %d file %s\n",pAllobjects->NUMBER,allobjects_name[fileIndex]);
        exit(-1);
      }


    }
    /* Reject any spatial bins not in every file */
    if (verbose) {
      fprintf(stderr,"Bin, binCount, acceptBin: ");
    }
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      if(binCount[spatial_bin] == 0) {
        acceptBin[spatial_bin] = 0;
      }
      if (verbose) {
        fprintf(stderr,"%d %6d %d |",spatial_bin,binCount[spatial_bin],acceptBin[spatial_bin]);
      }
    }
    if (verbose) {
      fprintf(stderr,"\n");
    }

    if (verbose) {
      fprintf(stderr,"Read %zu records matchedObjectCount %d duplicateREFCount %d multipleNoneCount %d initialNoneCount %d initialREFCount %d multipleREFCount %d galaxyCount %d multipleBlendGalaxyCount %d multipleBlendCount %d from %s\n",allobjects_nrecs,matchedObjectCount,duplicateREFCount,multipleNoneCount,initialNoneCount,initialREFCount,multipleREFCount,galaxyCount,multipleBlendGalaxyCount,multipleBlendCount,allobjects_name[fileIndex]);
    }
    /* Free up memory */
    if (allobjects_header != NULL) {
      table_hdrfree(allobjects_header);
      allobjects_header = NULL;
    }
    if (allobjects_handle[fileIndex] != NULL) {
      Close(allobjects_handle[fileIndex]);
      allobjects_handle[fileIndex] = NULL;
    }
    if (allobjects_descriptor != NULL) {
      Free(allobjects_descriptor);
      allobjects_descriptor = NULL;
    }
    if (allobjects_row != NULL) {
      table_rowfree(allobjects_row);
      allobjects_row = NULL;
    }
    if (allobjects_table != NULL) {
      free(allobjects_table);
      allobjects_table = NULL;
    }


  }

  /* Now write the output file */
  for (sextractor_index = 0; sextractor_index < (sextractorMaxNUMBER+1); sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];
    if (((pSextractor->selectFlag == (STAR_SELECTED|NONE_STAR)) ||
         (pSextractor->selectFlag == (STAR_SELECTED|NONE_STAR|MULTIPLE_NONE))) &&
        (acceptBin[pSextractor->spatial_bin] != 0)) {
      if (debugHandle != NULL) {
        fprintf(debugHandle,"ACCEPT selectFlag %d 0x%x\n",pSextractor->selectFlag,pSextractor->selectFlag);
      }

      if (outHandle == NULL) {
        outHandle = fopen(outfile,"wt");
        if (outHandle == NULL) {
          fprintf(stderr,"ERROR: filter_multiple  Failed to open the output file %s\n",outfile);
          exit(-1);
        } else {
          fprintf(outHandle,"X_IMAGE\tY_IMAGE\tMAG_ISO\tFLAGS\tELLIPTICITY\tNUMBER\n");
          fprintf(outHandle,"-------\t-------\t-------\t-----\t-----------\t------\n");
          if (verbose) {
            fprintf(stderr,"Output file %s\n",outfile);
          }
        }
      }


      fprintf(outHandle,"%f\t%f\t%f\t%d\t%f\t%d\n",
              pSextractor->X_IMAGE/pixelBinning,
              pSextractor->Y_IMAGE/pixelBinning,
              pSextractor->MAG_ISO,
              pSextractor->BFLAGS & FILTER_BMASK_SEXTRACTOR,
              pSextractor->ELLIPTICITY,
              pSextractor->NUMBER);
      outputCount++;
    } else {
      if ((debugHandle != NULL) && (acceptBin[pSextractor->spatial_bin] != 0)) {
        fprintf(debugHandle,"REJECT selectFlag %d 0x%x\n",pSextractor->selectFlag,pSextractor->selectFlag);
      }

    }
    if ((pSextractor->selectFlag & (NONE_STAR|REF_STAR)) != 0) {

      /* We need write this one to the combined file */
      if ((pSextractor->selectFlag & MULTIPLE_STAR) != 0) {
        galaxyFlag = 0;
        blendedFlux = 0.0;
        /* We need to select the best image and possibly set the multiple exposure blend flag */
        /* Start by searching for the brightest object */
        pSextractorSel = pSextractor;
        pAllobjectsSel = &pSextractorSel->allObjects;

        if (pAllobjectsSel->nextAllobject == 0) {
          fprintf(stderr,"ERROR: filter_multiple  filter_multiple software inconsistancy at line %d\n",__LINE__);
          exit(-1);
        }
        nextAllobject = pAllobjectsSel->NUMBER;
        pSavedAllobjects = pAllobjectsSel;
        while (nextAllobject != 0) {
          pSextractor3 = &sextractor_table[nextAllobject];
          pSavedAllobjects = &pSextractor3->allObjects;
          if ((((pSavedAllobjects->AFLAGS) >> GSC_CLASS_BIT) & CLASS_MASK) != 0) {
            galaxyCount++;
            galaxyFlag = 1;
          }
          if (pSavedAllobjects->Stdmag < 99.0) {
            tmpFlux = exp10(-pSavedAllobjects->Stdmag/2.5);
            blendedFlux += tmpFlux;

          }
          if (pSextractorSel != pSextractor3) {
            if (pAllobjectsSel->Stdmag > pSavedAllobjects->Stdmag) {
              pSextractorSel = pSextractor3;
              pAllobjectsSel = pSavedAllobjects;
            }
          }
          nextAllobject = pSavedAllobjects->nextAllobject;
        }
        /* Copy over our brightest */
        if (pSextractor != pSextractorSel) {
          memcpy(pSextractor,pSextractorSel,sizeof(SEXTRACTOR));
          pAllobjectsSel = &pSextractorSel->allObjects;
        }
        if (galaxyFlag) {
          /* If there is a galaxy in the mix, flag it unconditionally */
          pAllobjectsSel->AFLAGS &= ~(1<<FILTER_AFLAG_DEFECT); /* all blends are presumed to be defects */
          pAllobjectsSel->AFLAGS |=  (1<<FILTER_AFLAG_MULTIPLE_BLEND);
          multipleBlendGalaxyCount++;
        } else {
          blendedMag = -2.5*log10(blendedFlux);
          if ((pAllobjectsSel->Stdmag > 99.0) ||
              ((pAllobjectsSel->Stdmag - blendedMag) > MIN_BLEND_MAGNITUDE)) {
            multipleBlendCount++;
            pAllobjectsSel->AFLAGS &= ~(1<<FILTER_AFLAG_DEFECT); /* all blends are presumed to be defects */
            pAllobjectsSel->AFLAGS |= (1 << FILTER_AFLAG_MULTIPLE_BLEND);
          }

        }

      }
      if (combinedHandle == NULL) {
        combinedHandle = fopen(combinedfile,"wt");
        if (combinedHandle == NULL) {
          fprintf(stderr,"ERROR: filter_multiple  Failed to open the combined output file %s\n",combinedfile);
          exit(-1);
        } else {
          /* Note: this code must agree with recover_points.c */
          fprintf(combinedHandle,"REF\tNUMBER\tX_IMAGE\tY_IMAGE\tAFLAGS\tBFLAGS\tnpoints_local\tMAG_ISO\tra\tdec\tmagcal_iso\tmagcal_iso_rms\tmagcal_local\tmagcal_local_error\tmagcal_local_rms\tspatial_bin\tDate\tlimiting_mag_local\textinction\tStdmag\tcolor\tdra\tddec\tFLUX_ISO\tMAG_APER\tMAG_AUTO\tKRON_RADIUS\tBACKGROUND\tTHRESHOLD\tFLUX_MAX\tTHETA_J2000\tELLIPTICITY\tISOAREA_WORLD\tFWHM_IMAGE\tFWHM_WORLD\tISO0\tISO1\tISO2\tISO3\tISO4\tISO5\tISO6\tISO7\tplate_dist\tBlendedmag\tgsc_bin_index\tlocal_bin\tdradRMS2\texposureNumber\tsolutionNumber\tmagdep_bin\tmagcal_magdep\tmagcal_magdep_rms\n");
          fprintf(combinedHandle,"---\t------\t-------\t-------\t------\t------\t-------------\t-------\t--\t---\t----------\t--------------\t------------\t------------------\t----------------\t-----------\t----\t------------------\t----------\t------\t-----\t---\t----\t--------\t--------\t--------\t-----------\t----------\t---------\t--------\t-----------\t-----------\t-------------\t----------\t----------\t----\t----\t----\t----\t----\t----\t----\t----\t----------\t----------\t-------------\t---------------\t--------\t--------------\t--------------\t----------\t-------------\t-----------------\n");

          if (verbose) {
            fprintf(stderr,"Combined output file %s\n",combinedfile);
          }
        }
      }
      pSavedAllobjects = &pSextractor->allObjects;

      if (((pSextractor->selectFlag & NONE_STAR) != 0) &&
          ((pSextractor->selectFlag & MULTIPLE_NONE) != 0)) {
        /* Because we do not know which of the "NONE" positions are correct, we need to flag them and write them all out */
        pSavedAllobjects->AFLAGS |= (1<<FILTER_AFLAG_MULTIPLE_NONE);

      }


      fprintf(combinedHandle,"%s\t%d\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%.6f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%f\t%d\t%d\t%d\t%f\t%f\n",
              pSavedAllobjects->REF,
              pSavedAllobjects->NUMBER,
              pSavedAllobjects->X_IMAGE,
              pSavedAllobjects->Y_IMAGE,
              pSavedAllobjects->AFLAGS,
              pSavedAllobjects->BFLAGS,
              pSavedAllobjects->npoints_local,
              pSavedAllobjects->MAG_ISO,
              pSavedAllobjects->ra,
              pSavedAllobjects->dec,
              pSavedAllobjects->magcal_iso,
              pSavedAllobjects->magcal_iso_rms,
              pSavedAllobjects->magcal_local,
              pSavedAllobjects->magcal_local_error,
              pSavedAllobjects->magcal_local_rms,
              pSavedAllobjects->spatial_bin,
              pSavedAllobjects->Date,
              pSavedAllobjects->limiting_mag_local,
              pSavedAllobjects->extinction,
              pSavedAllobjects->Stdmag,
              pSavedAllobjects->color,
              pSavedAllobjects->dra,
              pSavedAllobjects->ddec,
              pSavedAllobjects->FLUX_ISO,
              pSavedAllobjects->MAG_APER,
              pSavedAllobjects->MAG_AUTO,
              pSavedAllobjects->KRON_RADIUS,
              pSavedAllobjects->BACKGROUND,
              pSavedAllobjects->THRESHOLD,
              pSavedAllobjects->FLUX_MAX,
              pSavedAllobjects->THETA_J2000,
              pSavedAllobjects->ELLIPTICITY,
              pSavedAllobjects->ISOAREA_WORLD,
              pSavedAllobjects->FWHM_IMAGE,
              pSavedAllobjects->FWHM_WORLD,
              pSavedAllobjects->ISO0,
              pSavedAllobjects->ISO1,
              pSavedAllobjects->ISO2,
              pSavedAllobjects->ISO3,
              pSavedAllobjects->ISO4,
              pSavedAllobjects->ISO5,
              pSavedAllobjects->ISO6,
              pSavedAllobjects->ISO7,
              pSavedAllobjects->plate_dist,
              pSavedAllobjects->Blendedmag,
              pSavedAllobjects->gsc_bin_index,
              pSavedAllobjects->local_bin,
              pSavedAllobjects->dradRMS2,
              pSavedAllobjects->exposureNumber,
              pSavedAllobjects->solutionNumber,
              pSavedAllobjects->magdep_bin,
              pSavedAllobjects->magcal_magdep,
              pSavedAllobjects->magcal_magdep_rms);


      combinedCount++;

      if (((pSextractor->selectFlag & NONE_STAR) != 0) &&
          ((pSextractor->selectFlag & MULTIPLE_NONE) != 0)) {
        /* Because we do not know which of the "NONE" positions are correct, we need to flag them and write them all out */
        pSextractorSel = pSextractor;
        pAllobjectsSel = &pSextractorSel->allObjects;
        if (pAllobjectsSel->nextAllobject == 0) {
          fprintf(stderr,"ERROR: filter_multiple  filter_multiple software inconsistancy at line %d\n",__LINE__);
          exit(-1);
        }
        nextAllobject = pAllobjectsSel->nextAllobject;
        curNoneMultiple = 1;
        while (nextAllobject != 0) {
          pSextractor3 = &sextractor_table[nextAllobject];
          pSavedAllobjects = &pSextractor3->allObjects;
          pSavedAllobjects->AFLAGS |= (1<<FILTER_AFLAG_MULTIPLE_NONE);
          curNoneMultiple++;
          fprintf(combinedHandle,"%s\t%d\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%.6f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%f\t%d\t%d\t%d\t%f\t%f\n",
                  pSavedAllobjects->REF,
                  pSavedAllobjects->NUMBER,
                  pSavedAllobjects->X_IMAGE,
                  pSavedAllobjects->Y_IMAGE,
                  pSavedAllobjects->AFLAGS,
                  pSavedAllobjects->BFLAGS,
                  pSavedAllobjects->npoints_local,
                  pSavedAllobjects->MAG_ISO,
                  pSavedAllobjects->ra,
                  pSavedAllobjects->dec,
                  pSavedAllobjects->magcal_iso,
                  pSavedAllobjects->magcal_iso_rms,
                  pSavedAllobjects->magcal_local,
                  pSavedAllobjects->magcal_local_error,
                  pSavedAllobjects->magcal_local_rms,
                  pSavedAllobjects->spatial_bin,
                  pSavedAllobjects->Date,
                  pSavedAllobjects->limiting_mag_local,
                  pSavedAllobjects->extinction,
                  pSavedAllobjects->Stdmag,
                  pSavedAllobjects->color,
                  pSavedAllobjects->dra,
                  pSavedAllobjects->ddec,
                  pSavedAllobjects->FLUX_ISO,
                  pSavedAllobjects->MAG_APER,
                  pSavedAllobjects->MAG_AUTO,
                  pSavedAllobjects->KRON_RADIUS,
                  pSavedAllobjects->BACKGROUND,
                  pSavedAllobjects->THRESHOLD,
                  pSavedAllobjects->FLUX_MAX,
                  pSavedAllobjects->THETA_J2000,
                  pSavedAllobjects->ELLIPTICITY,
                  pSavedAllobjects->ISOAREA_WORLD,
                  pSavedAllobjects->FWHM_IMAGE,
                  pSavedAllobjects->FWHM_WORLD,
                  pSavedAllobjects->ISO0,
                  pSavedAllobjects->ISO1,
                  pSavedAllobjects->ISO2,
                  pSavedAllobjects->ISO3,
                  pSavedAllobjects->ISO4,
                  pSavedAllobjects->ISO5,
                  pSavedAllobjects->ISO6,
                  pSavedAllobjects->ISO7,
                  pSavedAllobjects->plate_dist,
                  pSavedAllobjects->Blendedmag,
                  pSavedAllobjects->gsc_bin_index,
                  pSavedAllobjects->local_bin,
                  pSavedAllobjects->dradRMS2,
                  pSavedAllobjects->exposureNumber,
                  pSavedAllobjects->solutionNumber,
                  pSavedAllobjects->magdep_bin,
                  pSavedAllobjects->magcal_magdep,
                  pSavedAllobjects->magcal_magdep_rms);
          nextAllobject = pSavedAllobjects->nextAllobject;
        }
        if (curNoneMultiple > maxNoneMultiple) {
          maxNoneMultiple = curNoneMultiple;
          maxNoneNUMBER = pSavedAllobjects->NUMBER;
        }
      }
    }
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


  if (outHandle != NULL) {
    fclose(outHandle);
  }
  time(&curTime);
  curTime -= startTime;
  if (duplicateREFCount > 0) {
    printf("ERRORD: filter_multiple  duplicateREFCount is %d for %s s%d\n",duplicateREFCount,fileroot,solutionNumber);
    SetMosaicFitWCS(pConnection,"DuplicateREFCount",fileroot,solutionNumber-1,1,readFlag,&FitWCS);
  } else {
    SetMosaicFitWCS(pConnection,"DuplicateREFCount",fileroot,solutionNumber-1,0,readFlag,&FitWCS);
  }


  mysql_close(pConnection);

  fprintf(stdout,"filter_multiple starsin %zu matchedObjectCount %d duplicateREFCount %d multipleNoneCount %d initialNoneCount %d initialREFCount %d multipleREFCount %d galaxyCount %d multipleBlendGalaxyCount %d multipleBlendCount %d none output %d combined output %d maxNoneCount %d maxNoneNUMBER %d seconds %ld for %s\n",
          sextractor_nrecs,
          matchedObjectCount,
          duplicateREFCount,
          multipleNoneCount,
          initialNoneCount,
          initialREFCount,
          multipleREFCount,
          galaxyCount,
          multipleBlendGalaxyCount,
          multipleBlendCount,
          outputCount,
          combinedCount,
          maxNoneMultiple,
          maxNoneNUMBER,
          curTime,
          fileroot);

  if (debugHandle != NULL) {
    fclose(debugHandle);
  }

  return 0;
}


 /*
 * Feb 10, 2009 Edward J. Los - Initial version
 * Oct  5, 2009 Edward J. Los - Add the NUMBER column to match this file with Sextractor files
 * Oct 20, 2009 Edward J. Los - Count the number of REF matches
 * Nov  3, 2009 Edward J. Los - Write a new allobjects file with the combined result
 * Dec  7, 2009 Edward J. Los - List completely unmatched objects with all sets of coordinates and
 *                              add FILTER_AFLAG_MULTIPLE_NONE to note the uncertainty in time and position.
 * Nov  4, 2011 Edward J. Los - Add magnitude-dependent parameters (magdep_bin, magcal_magdep, and magcal_magdep_rms).
 * Feb 27, 2012 Edward J. Los - Support pass2 qualifiers; do not write an output file for pass 2 qualifiers.
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support.
 * Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * May 28, 2018 Edward J. Los - Add gaia support
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 *
 * TODO: Need to search for missing emulsion
 * TODO: Need to recognize non-adjacent spatial bins
 * NOTE: There is no need to recognize Pickering Wedge objects because they will be picked up as a double exposure.
 */
