// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 * "This routine is a version of update_summary which does not write the
 * starsummary MySQL tables, but writes directly to a Starbase file". In turn,
 * update_summary will "compute the DASCH median magnitude for all stars that
 * need updating".
 *
 * Inputs:
 *
 * - $DASCH_PHOT_MAGNITUDES@KEY@/magNNN/magNNN/magNNN/magNNNNNNNNN.db
 *   - This tool iterates over every single potential magfile name, so it is
 *     very slow even if the magnitude database is small.
 *
 * Outputs:
 *
 * - Logfile associated with `-o` argument
 * - "ID database" associated with `-u` argument
 *   - traditionally `$DASCH_INGEST/id_${catname}_${datecode}.db`
 *
 * Database updates:
 *
 * - TBC?
 *
 * Environment variables:
 *
 * - DASCH_CATALOG
 * - DASCH_PHOT_MAGNITUDES@KEY@ with @KEY@ a digit or ""
 * - DASCH_PHOT_MYSQLHOST
 * - DASCH_PHOT_PASSWORD
 * - DASCH_PHOT_USERNAME
 * - DASCH_MYSQLHOST
 * - DASCH_PASSWORD
 * - DASCH_USERNAME
 */

#include <math.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include <table.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "photometryutils.h"


#define MAX_BUFFER 256
#define STAR_LIMIT     50000
#define SUMMARY_LIMIT  50000
#define MAX_ALLOCATION 0x10000000
#define AVE_STARS 1000 /* Expected average number of objects per bin */
#define PRINT_INTERVAL 21600 /* seconds; = 6 hr */


extern GSCBIN gscBin64;


static PGSCBIN pGscBin = &gscBin64;


int
main(int argc, char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  time_t deltaTime;
  time_t oldTime;
  time_t printTime = 0;
  int verbose = 0;
  int provideWebSummary = 0;
  int debugMode = 0;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  char gscbinfile[MAX_BUFFER];
  int nvals;
  FILE *outHandle = NULL;
  char timestr[100];
  struct tm *ptr;
  int curMagnitudes;
  long long totalMagnitudes = 0;
  int staleVersionIdCount = 0;
  char *dotPtr;
  int processedCount = 0;
  int minGoodStars = MIN_SUMMARY_GOODSTARS;

  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;

  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  PSTARENTRY pStarTable;
  PSTARENTRY pCurStarEntry = NULL;
  PSTARENTRY pCurStarEntry2 = NULL;

  PPHOTSTARIMAGE pMagnitudeTable = NULL;
  PPHOTSTARIMAGE pCurSumstarimage;
  PFILESTARIMAGE pFileStarImage = NULL;

  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon = &fileCommon;
  PGALAXYCOMMON pGalaxyCommon = &pFileCommon->galaxyCommon;

  int seriesId;
  int plateNumber;
  int *seriesVersionArray;
  int magnitudeIndex;
  int magnitudeAlloc = 0;
  int goodStars = 0;
  int binCount = 0;
  int noGscEntryCount = 0;

  char debug_name[MAX_BUFFER];
  FILE *debug_handle = NULL;

  int curStars;
  int totalStars = 0;
  int starIndex;
  int starIndex2;
  int summaryCount = 0;

  int oldQueryCount = 0;
  int oldTotalStars = 0;
  int oldGoodStars = 0;
  int oldBinCount = 0;
  long long oldTotalMagnitudes = 0;
  PARAMETERVECTORSTORE vectorStore;
  PPARAMETERVECTORSTORE pVectorStore = &vectorStore;
  int enableRematch = 0;
  int curGscBinIndex = 0;
  int maxGscBinIncrement;
  char indexname[MAX_BUFFER];
  File indexHandle;
  char* catalogdir;
  char* slashPtr;
  char galaxy_name[MAX_BUFFER];
  char catalogname[MAX_BUFFER];
  File catalogHandle;
  int cur_gsc_bin_index = -1;
  STARINDEX curStarIndex;
  PSTARINDEX pCurStarIndex = &curStarIndex;
  PGSCIMAGE pGscImageTable = NULL;
  int gscImageAlloc = 0;
  PGSCIMAGE pCurGscImage;
  int readItems;
  int gscIndex;
  int excludeSeriesCount = 0;
  int magnitudesModulus = MAGNITUDES_MODULUS;
  int dumpAllFlag = 0;
  FILE *headerHandle = NULL;
  double lat;
  double lon;

  /* The following four values are used for ChargeTime
   * Index 0: non database work
   *       1: GetStarEntry
   *       2: GetSummaryMagnitudes
   *       3: LOAD starsummary
   *       4: LOAD stars
   */
  int usecIndex;
  long long curUsec;
  long long deltaTable[MAX_CHARGE_ENTRY];
  long long totalTable[MAX_CHARGE_ENTRY];

  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];

  double minimumGalacticLatitude = -90.0;
  char idfile[MAX_BUFFER];
  FILE *idHandle = NULL;
  PHOTTARGET target;
  PPHOTTARGET pTarget = &target;
  PPHOTTARGET target_table = NULL;
  int gscBinIndexCurrent = -1;

  memset(pVectorStore, 0, sizeof(PARAMETERVECTORSTORE));

  catalogString[0] = 0;
  ChargeTime(RESET_ALL_ENTRY, &curUsec, deltaTable, totalTable);
  SetQueryCount(0);

  outfile[0] = 0;
  qualifier[0] = 0;
  idfile[0] = 0;
  gscbinfile[0] = 0;

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d\n",argstr,argc);
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
            if (strlen(outfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'g': /* limited bin list to search */
        case 'G':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(gscbinfile,*++argv,MAX_BUFFER-2);
            if (strlen(gscbinfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'n': /* Minimum number of good stars */
        case 'N':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&minGoodStars);
            if (nvals != 1) {
              printf("ERROR: Unable to decode minGoodStars %s\n",*argv);
              errorFlag = 1;
            } else {
              if (minGoodStars < 1) {
                printf("ERROR: illegal value for minGoodStars: %d  Using default %d\n",minGoodStars,MIN_SUMMARY_GOODSTARS);
                minGoodStars = MIN_SUMMARY_GOODSTARS;
              }
            }
          }
          break;

        case 'l': /* Minimum galactic latitude */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&minimumGalacticLatitude);
            if (nvals != 1) {
              printf("ERROR: Unable to decode minimum galactic latitude %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'u': /* id file name */
        case 'U':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(idfile,*++argv,MAX_BUFFER-2);
            if (strlen(idfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              sprintf(catalogString,"%d",catalogNumber);
              strcpy(qualifier,*argv);
            }
          }
          break;

        case 'e': /* Exclude bad astrometry series */
        case 'E':
          excludeSeriesCount = NUM_EXCLUDE_SERIES;
          break;

        case 'w': /* Brief run to get summary totals only */
        case 'W':
          provideWebSummary = 1;
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'd': /* debug mode */
        case 'D':
          verbose = 1;
          debugMode = 1;
          break;

        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  /* Validate arguments */

  if (outfile[0] == 0) {
    printf("ERROR: No output filename was specified\n");
    errorFlag = 1;
  }

  if (provideWebSummary == 0) {
    outHandle = fopen(outfile, "a+t");
    if (outHandle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to open the output file %s\n", outfile);
    } else if (verbose) {
      printf("Output file %s\n",outfile);
    }

    if (idfile[0] != 0) {
      idHandle = fopen(idfile, "wt");
      if (idHandle == NULL) {
        errorFlag = 1;
        printf("ERROR: Failed to open the output file %s\n", idfile);
      } else {
        if (verbose) {
          printf("Id file %s\n",idfile);
        }

        /* This header implies that pFileCommon->extraColumnFlag == 1 */
        fprintf(idHandle,"REF\tdrad\tdradB\tmax_drad\tmax_dradB\tyrbegin\tyrend\tnpoints\tmin_local\tmax_local\trange_local\tmin_local2\tmax_local2\trange_local2\tmedian_local\trms_local\tngood\tngoodB\tclip_median_local\tclip_rms_local\tclip_ngood\tmedian_iso\trange_iso\terror_bar_factor\tMalmquist_factor\tMalmquist_factorB\tDamon_factor\tSextractor_Blend\tnblend\tnNonDamonBlue\tmedNonDamonBlue\tnonDamonBluerms\tnDamonBlue\tmedDamonBlue\tdamonBlueRms\tmagvslimitingcorr\tmagvsracorr\tmagvsdeccorr\trmsdradrms2\trarms\tdecrms\tnburst\tnburst2\tnburst3\tnburst4\tndip\tndip2\tndip3\tndip4\tndev3\tndev2\tadjacentburstdip\tadjacentburstdip2\tadjacentburstdip3\tlightcurverms1\tlightcurverms2\tlightcurverms3\tlightcurverms4\tlightcurverms5\tslope_all\tslope_all_err\tnslope_all\tslope_60\tslope_60_err\tnslope_60\trms_factor\tdmagcatalog\tStdmag\tcolor\tra\tdeclination\tMAGFlag\tgscclass\tVFlag\tRaPM\tDecPM\tkeplerField\tversionId\tnearbyObjects\tlat\tnearbyREFflag\treleaseField\tpeakDays\tpeakYear\tpeakMag\tpeakSlope\tpeakRMS\tpeakCount\tpeakUpperCount\tpeakNumber\tpeakMaxMag\tpeakMinMag\tpeakEventCount\tpeakDradRMS\tpeakEventCount2\tpeakEventCount3\tpeakExtra\tpeakCountNF\tpeakCountWF\tpeakOutside\tpeakOutsideNF\tpeakOutsideWF\tpeakDefectCount\tpeakRA\tpeakDec\tpeakNearbyDistance\tpeakMultipleCount\tgsc_bin_index\tpeakDradRMS2\tpeakDradRMS3\tpeakBadColorCount\tpeakLimitingYears\tpeakLimitingPoints\tpeakEvaluation\tpeakNoDefectMag\tpeakLongOutburst\n");
        fprintf(idHandle,"---\t----\t-----\t--------\t---------\t-------\t-----\t-------\t---------\t---------\t-----------\t----------\t----------\t------------\t------------\t---------\t-----\t------\t-----------------\t--------------\t----------\t----------\t---------\t----------------\t----------------\t-----------------\t------------\t----------------\t------\t-------------\t---------------\t---------------\t----------\t------------\t------------\t-----------------\t-----------\t------------\t-----------\t-----\t------\t------\t-------\t-------\t-------\t----\t-----\t-----\t-----\t-----\t-----\t----------------\t-----------------\t-----------------\t--------------\t--------------\t--------------\t--------------\t--------------\t---------\t-------------\t----------\t--------\t------------\t---------\t----------\t-----------\t------\t-----\t--\t-----------\t-------\t--------\t-----\t----\t-----\t-----------\t---------\t-------------\t---\t-------------\t------------\t--------\t--------\t-------\t---------\t-------\t---------\t--------------\t----------\t----------\t----------\t--------------\t-----------\t---------------\t---------------\t----------\t-----------\t-----------\t-----------\t-------------\t-------------\t--------------\t------\t-------\t------------------\t-----------------\t-------------\t------------\t------------\t-----------------\t-----------------\t------------------\t--------------\t---------------\t----------------\n");
      }
    } else {
      printf("ERROR: No id filename was specified\n");
      errorFlag = 1;
    }
  }

  catalogdir = getenv("DASCH_CATALOG");
  if (catalogdir == NULL) {
    printf("DASCH_CATALOG is not defined\n");
    return -1;
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
  }

  if (strstr(qualifier,"apass")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"apass.dat");
  }

  if (strstr(qualifier,"gaia")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"gaiadr2.dat");
  }

  if (strstr(qualifier,"atlas")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"atlas.dat");
  }

  if (strstr(qualifier,"experimental")) {
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
  strcat(indexname,".idx");

  catalogHandle = Open(catalogname,"r");
  if (catalogHandle == NULL) {
    printf("ERROR Could not open catalog file %s\n", catalogname);
    errorFlag = 1;
  }

  indexHandle = Open(indexname,"r");
  if (indexHandle == NULL) {
    printf("ERROR Could not open catalog index file %s\n", indexname);
    errorFlag = 1;
  }

  if (errorFlag) {
    printf("Usage: update_summary2 options\n");
    printf("  options: -v verbose\n");
    printf("           -d debugging mode\n");
    printf("           -o <output file>\n");
    printf("           -n <minimum stars per lightcurve>\n");
    printf("           -q input catalog qualifer\n");
    printf("           -e excludes ac,ca,ax, and am series\n");
    printf("           -l <minimum galactic latitude>, default %f\n", DEFAULT_GAL_LATITUDE);
    printf("           -u <id summary file>\n");
    printf("           -g gsc bin list for limited execution.\n");
    printf("           -w provide only monthly website statistics\n");
    return -1;
  }

  maxGscBinIncrement = magnitudesModulus * (MAX_ALLOCATION / (AVE_STARS * (sizeof(PHOTSTARIMAGE) + sizeof(FILESTARIMAGE))));

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pConnection);
  dasch_init_photdb(pPhotConnection);

  InitSeriesTable(pConnection, pPhotConnection);

  if (GetPhotometryGlobal(pPhotConnection, pPhotGlobal) != 1) {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }

  if (pPhotGlobal->magnitudeFile == PHOT_MAGNITUDEFILE_YES) {
    printf("File Header size %zu File Image size %zu\n", sizeof(FILESTARHEADER), sizeof(FILESTARIMAGE));
    if (outHandle != NULL) {
      fprintf(outHandle, "File Header size %zu File Image size %zu\n", sizeof(FILESTARHEADER), sizeof(FILESTARIMAGE));
    }
  } else {
    printf("ERROR: magnitudes stored in the MySQL database are no longer supported\n");
    exit(-1);
  }

  InitFileCommon(stdout, pFileCommon, pConnection, pPhotConnection, catalogNumber, 1);
  InitMaxPlateNumber(pConnection, pFileCommon->maxPlateNumber);
  pFileCommon->provideWebSummary = provideWebSummary;
  pFileCommon->extraColumnFlag = 1;

  if (enableRematch) {
    pFileCommon->enableRematch = 1;
  }

  if (InitGscBinList(pFileCommon,gscbinfile) != 0) {
    exit(-1);
  }

  if (pFileCommon->gscBinIndexCount > 0) {
    gscBinIndexCurrent = 0;
  }

  mysql_close(pConnection);

  printf(
    "update_summary2 of %s %s, outfile %s minimum stars %d starLimit %d summaryLimit %d sizeof(PHOTSTARIMAGE) %zu maxGscBinIncrement %d excludeSeriesCount %d magnitudeFile %d minVersionId%d %d minGalacticLatitude %f\n",
    __DATE__,
    __TIME__,
    outfile,
    minGoodStars,
    STAR_LIMIT,
    SUMMARY_LIMIT,
    sizeof(PHOTSTARIMAGE),
    maxGscBinIncrement,
    excludeSeriesCount,
    1,
    catalogNumber,
    pPhotGlobal->minVersionId[catalogNumber],
    minimumGalacticLatitude
  );

  if (outHandle != NULL) {
    fprintf(
      outHandle,
      "update_summary2 of %s %s, outfile %s minimum stars %d starLimit %d summaryLimit %d sizeof(PHOTSTARIMAGE) %zu maxGscBinIncrement %d excludeSeriesCount %d magnitudeFile %d minVersionId%d %d minGalacticLatitude %f\n",
      __DATE__,
      __TIME__,
      outfile,
      minGoodStars,
      STAR_LIMIT,
      SUMMARY_LIMIT,
      sizeof(PHOTSTARIMAGE),
      maxGscBinIncrement,
      excludeSeriesCount,
      1,
      catalogNumber,
      pPhotGlobal->minVersionId[catalogNumber],
      minimumGalacticLatitude
    );
  }

  time(&startTime);
  oldTime = startTime;
  printTime = startTime;
  ptr = localtime(&startTime);
  strftime(timestr, 25, "%Y-%m-%dT%H-%M-%S", ptr);

  printf("update_summary2 starting at %s\n", timestr);
  if (outHandle != NULL) {
    fprintf(outHandle,"update_summary2 starting at %s\n", timestr);
    fflush(outHandle);
  }

  pStarTable = (PSTARENTRY) calloc(STAR_LIMIT, sizeof(STARENTRY));
  if (pStarTable == NULL) {
    fprintf(stderr,"ERROR: failed to allocate pStarTable\n");
    exit(-1);
  }

  strcpy(galaxy_name, catalogdir);
  slashPtr = strrchr(galaxy_name, '/');
  if (slashPtr != NULL) {
    slashPtr++;
  } else {
    slashPtr = galaxy_name;
  }
  *slashPtr = 0;
  strcat(galaxy_name, "galaxy.dat");

  if (provideWebSummary == 0) {
    if (PopulateGalaxyTree(pGalaxyCommon, galaxy_name, 0) != 0) {
      exit(-1);
    }
  }

  /* Enter the main loop */

  while (1) {
    time(&curTime);
    curTime -= printTime;

    if (curTime > PRINT_INTERVAL) {
      time(&curTime);
      printTime = curTime;
      curTime -= startTime;

      if (pCurStarEntry2 != NULL) {
        if (outHandle != NULL) {
          fprintf(
            outHandle,
            "At gsc_bin_index %d and REF %s and ra %5.1f declination %5.1f\n",
            pCurStarEntry2->gsc_bin_index,
            pCurStarEntry2->REF,
            pCurStarEntry2->ra,
            pCurStarEntry2->dec
          );
        }
      }

      if (outHandle != NULL) {
        fprintf(
          outHandle,
          "Elapsed Time: %ld seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins, %lld magnitudes\n",
          curTime,
          GetQueryCount(),
          totalStars,
          goodStars,
          magnitudeAlloc,
          pVectorStore->vectorAlloc,
          binCount,
          totalMagnitudes
        );

        fprintf(
          outHandle,
          " queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec %f  magnitude/sec:%f\n",
          (1.0 * GetQueryCount()) / (1.0 * curTime),
          (1.0 * totalStars) / (1.0 * curTime),
          (1.0 * goodStars ) / (1.0 * curTime),
          (1.0 * binCount  ) / (1.0 * curTime),
          (1.0 * totalMagnitudes) / (1.0 * curTime)
        );

        fprintf(outHandle, "Intermediate Time (sec) ");
        for (usecIndex = 0; usecIndex < MAX_CHARGE_ENTRY; usecIndex++) {
          double outSeconds;
          outSeconds = (1.0 * totalTable[usecIndex]) / 1000000.0;
          fprintf(outHandle, "%10.3f ", outSeconds);
        }
        fprintf(outHandle, "\n");

        fflush(outHandle);
      }
    }

    if (gscBinIndexCurrent >= 0) {
      /* Select the next gsc bin index from the list */
      while (gscBinIndexCurrent < pFileCommon->gscBinIndexCount && pFileCommon->gscBinIndexList[gscBinIndexCurrent] < curGscBinIndex) {
        gscBinIndexCurrent++;
      }

      if (gscBinIndexCurrent >= pFileCommon->gscBinIndexCount) {
        curGscBinIndex = pGscBin->total_gsc_bins;
      } else {
        curGscBinIndex = pFileCommon->gscBinIndexList[gscBinIndexCurrent];
        gscBinIndexCurrent++;
      }
    }

    if (curGscBinIndex >= pGscBin->total_gsc_bins) {
      /* We are done */
      break;
    }

    if (debugMode) {
      if (debug_handle != NULL) {
        fclose(debug_handle);
        debug_handle = NULL;
      }

      debug_handle = fopen(debug_name, "wt");
      if (debug_handle == NULL) {
        printf("ERROR: Failed to open summary file %s\n", debug_name);
        exit(-1);
      }
    }

    processedCount = 0;

    if (minimumGalacticLatitude > -90.0) {
      GetBinCenter(pGscBin, curGscBinIndex, &lon, &lat, "search_none");
      wcscon(WCS_J2000, WCS_GALACTIC, 2000.0, 2000.0, &lon, &lat, 2000.0);
      if (lat < minimumGalacticLatitude) {
        curGscBinIndex++;
        continue;
      }
    }

    ChargeTime(0, &curUsec, deltaTable, totalTable);

    curMagnitudes = GetFileSummaryMagnitudes(
      pGscBin,
      pFileCommon,
      curGscBinIndex,
      curGscBinIndex,
      0,
      1,
      0,
      dumpAllFlag,
      headerHandle,
      catalogString,
      0,
      NULL,
      0,
      1
    );

    /* For now, just copy the results over */

    if (magnitudeAlloc < curMagnitudes) {
      free(pMagnitudeTable);
      magnitudeAlloc = curMagnitudes;
      pMagnitudeTable = (PPHOTSTARIMAGE) calloc(magnitudeAlloc, sizeof(PHOTSTARIMAGE));
      if (pMagnitudeTable == NULL) {
        printf("ERROR: failed to allocate %d PHOTSTARIMAGEs in GetSummaryMagnitudes\n", magnitudeAlloc);
        exit(-1);
      }
    }

    for (magnitudeIndex = 0; magnitudeIndex < curMagnitudes; magnitudeIndex++) {
      pCurSumstarimage = &pMagnitudeTable[magnitudeIndex];
      pFileStarImage = &pFileCommon->magnitudeBuffer[magnitudeIndex];

      memset(pCurSumstarimage, 0, sizeof(PHOTSTARIMAGE));
      pCurSumstarimage->pFileStarImage = pFileStarImage;
      GetREF(pFileStarImage->REFNumber, pCurSumstarimage->REF, 0, 1);
      strcpy(pCurSumstarimage->series, GetSeriesString(pFileStarImage->seriesId, 1));

      if (GetREF(pFileStarImage->REFNumber, pCurSumstarimage->REF, 0, 0) != 0) {
        printf("ERROR: corrupt REFNumber at magnitude index %d magnitudeIndex\n", magnitudeIndex);
        exit(-1);
      }
    }

    totalMagnitudes += curMagnitudes;
    ChargeTime(2, &curUsec, deltaTable, totalTable);

    for (magnitudeIndex = 0; magnitudeIndex < curMagnitudes; magnitudeIndex++) {
      pCurSumstarimage = &pMagnitudeTable[magnitudeIndex];
      pFileStarImage = &pFileCommon->magnitudeBuffer[magnitudeIndex];

      seriesId = pFileStarImage->seriesId;
      if (seriesId < 0 || seriesId > MAX_SERIES) {
        printf(
          "ERROR: illegal seriesId for seriesId %d, plateNumber %d, NUMBER %lld\n",
          seriesId,
          pFileStarImage->plateNumber,
          pFileStarImage->REFNumber
        );
        exit(-1);
      }

      plateNumber = pFileStarImage->plateNumber;
      if (plateNumber < 0 || plateNumber > pFileCommon->maxPlateNumber[seriesId]) {
        printf(
          "ERROR: illegal plateNumber for seriesId %d, plateNumber %d, NUMBER %lld\n",
          seriesId,
          pFileStarImage->plateNumber,
          pFileStarImage->REFNumber
        );
        exit(-1);
      }

      seriesVersionArray = pFileCommon->seriesVersionTable[seriesId];
      if (seriesVersionArray == NULL) {
        seriesVersionArray = (int *) calloc(pFileCommon->maxPlateNumber[seriesId] + 1, sizeof(int));
        if (seriesVersionArray == NULL) {
          printf("ERROR: failed to allocate seriesVersion Array of size %d\n", pFileCommon->maxPlateNumber[seriesId] + 1);
          exit(-1);
        }

        pFileCommon->seriesVersionTable[seriesId] = seriesVersionArray;
      }

      if (seriesVersionArray[plateNumber] < pFileStarImage->versionId) {
        seriesVersionArray[plateNumber] = pFileStarImage->versionId;
      }
    }

    /* Find out how many unique REF numbers there are in the table */

    curStars = 0;

    for (magnitudeIndex = 0; magnitudeIndex < curMagnitudes; magnitudeIndex++) {
      pCurSumstarimage = &pMagnitudeTable[magnitudeIndex];
      pFileStarImage = &pFileCommon->magnitudeBuffer[magnitudeIndex];
      if (pFileStarImage->REFNumber == 0) {
        continue;
      }

      for (starIndex = 0; starIndex < curStars; starIndex++) {
        pCurStarEntry = &pStarTable[starIndex];
        if (pFileStarImage->REFNumber == pCurStarEntry->REFNumber) {
          break;
        }
      }

      if (starIndex >= curStars) {
        pCurStarEntry = &pStarTable[curStars];
        memset(pCurStarEntry, 0, sizeof(STARENTRY));
        pCurStarEntry->REFNumber = pFileStarImage->REFNumber;
        pCurStarEntry->gsc_bin_index = pFileStarImage->gsc_bin_index;
        pCurStarEntry->versionId = pFileStarImage->versionId;

        if (pFileStarImage->ra_2 < 900 && pFileStarImage->dec_2 <= 90) {
          pCurStarEntry->ra = pFileStarImage->ra_2;
          pCurStarEntry->dec = pFileStarImage->dec_2;
          pCurStarEntry->RaPM = pFileStarImage->RaPM;
          pCurStarEntry->DecPM = pFileStarImage->DecPM;
        }

        pCurStarEntry->updateflag = 1;
        GetREF(pCurStarEntry->REFNumber, pCurStarEntry->REF, 0, 0);

        curStars++;
        if (curStars > STAR_LIMIT - 1) {
          printf("ERROR: curStars %d exceeds STAR_LIMIT %d curMagnitudes %d\n", curStars, STAR_LIMIT, curMagnitudes);
          exit(-1);
        }
      }
    }

    totalStars += curStars;

    if (debugMode) {
      if (debug_handle != NULL) {
        fclose(debug_handle);
        debug_handle = NULL;
      }
      debug_handle = fopen(debug_name, "wt");
      if (debug_handle == NULL) {
        fprintf(stderr,"ERROR: Failed to open summary file  %s\n", debug_name);
        exit(-1);
      }
    }

    processedCount = 0;

    for (starIndex = 0; starIndex < curStars; starIndex++) {
      if (errorFlag) {
        break;
      }

      pCurStarEntry = &pStarTable[starIndex];
      if (pCurStarEntry->processedflag != 0) {
        /* We processed this one already */
        continue;
      }

      if (debugMode) {
        time(&curTime);
        curTime -= startTime;
        fprintf(debug_handle, "starIndex %d finished getting magnitudes for %d at %ld\n", starIndex, pCurStarEntry->gsc_bin_index, curTime);
      }

      if (curMagnitudes < 0) {
        errorFlag = 1;
        break;
      }

      binCount++;

      for (starIndex2 = starIndex; starIndex2 < curStars; starIndex2++) {
        if (errorFlag) {
          break;
        }

        pCurStarEntry2 = &pStarTable[starIndex2];
        if (pCurStarEntry2->processedflag != 0) {
          /* We processed this one already */
          continue;
        }

        pCurStarEntry2->processedflag = 1;
        processedCount++;

        if (cur_gsc_bin_index != pCurStarEntry2->gsc_bin_index) {
          cur_gsc_bin_index = pCurStarEntry2->gsc_bin_index;
          if (cur_gsc_bin_index >= pGscBin->total_gsc_bins) {
            printf(
              "ERROR: gsc_bin_index %d exceeds total %d in update_summary2 for REFNUMBER %lld\n",
              cur_gsc_bin_index,
              pGscBin->total_gsc_bins,
              pCurStarEntry2->REFNumber
            );
            cur_gsc_bin_index = -1;
            memset(pCurStarIndex, 0, sizeof(STARINDEX));
          } else {
            Seek(indexHandle,cur_gsc_bin_index * sizeof(STARINDEX), SEEK_SET);
            readItems = Read(indexHandle, pCurStarIndex, sizeof(STARINDEX), 1);

            if (readItems != 1) {
              printf(
                "ERROR reading star index file errno: %d %s REF %s ra %f dec %f gsc_bin_index %d\n",
                errno,
                strerror(errno),
                pCurStarEntry2->REF,
                pCurStarEntry2->ra,
                pCurStarEntry2->dec,
                pCurStarEntry2->gsc_bin_index
              );
              cur_gsc_bin_index = -1;
              memset(pCurStarIndex, 0, sizeof(STARINDEX));
            } else {
              if (pCurStarIndex->numStars > gscImageAlloc) {
                gscImageAlloc = pCurStarIndex->numStars + 1000;
                pCurGscImage = realloc(pGscImageTable, gscImageAlloc * sizeof(GSCIMAGE));
                if (pCurGscImage == NULL) {
                  printf("ERROR: failed to realloc pGscImageTable of size %zu\n", gscImageAlloc * sizeof(GSCIMAGE));
                  exit(-1);
                }

                pGscImageTable = pCurGscImage;
                pCurGscImage = NULL;
              }

              Seek(catalogHandle, pCurStarIndex->offset, SEEK_SET);
              readItems = Read(catalogHandle, pGscImageTable, sizeof(GSCIMAGE), pCurStarIndex->numStars);
              if (readItems != pCurStarIndex->numStars) {
                printf("ERROR reading gsc catalog file\n");
                exit(-1);
              }
            }
          }
        }

        if (curMagnitudes >= minGoodStars) {
          /* Because the ra and dec has been modified for proper motion, get the original ra and dec from the GSC2.3.2 catalog */
          for (gscIndex = 0; gscIndex < pCurStarIndex->numStars; gscIndex++) {
            pCurGscImage = &pGscImageTable[gscIndex];

            if (pCurGscImage->REFNumber == pCurStarEntry2->REFNumber) {
              pCurStarEntry2->ra = pCurGscImage->ra;
              pCurStarEntry2->dec = pCurGscImage->dec;
              pCurStarEntry2->MAGFlag = pCurGscImage->MAGFlag;
              pCurStarEntry2->class = pCurGscImage->class;
              pCurStarEntry2->VFlag = pCurGscImage->VFlag;
              pCurStarEntry2->RaPM = pCurGscImage ->RaPM;
              pCurStarEntry2->DecPM = pCurGscImage->DecPM;

              pCurStarEntry2->Stdmag = pCurGscImage->Stdmag;
              pCurStarEntry2->color = pCurGscImage->color;
              pCurStarEntry2->MAGFlag = pCurGscImage->MAGFlag;
              pCurStarEntry2->class = pCurGscImage->class;
              pCurStarEntry2->VFlag = pCurGscImage->VFlag;

              pCurStarEntry->ra = pCurGscImage->ra;
              pCurStarEntry->dec = pCurGscImage->dec;
              pCurStarEntry->MAGFlag = pCurGscImage->MAGFlag;
              pCurStarEntry->class = pCurGscImage->class;
              pCurStarEntry->VFlag = pCurGscImage->VFlag;
              pCurStarEntry->RaPM = pCurGscImage ->RaPM;
              pCurStarEntry->DecPM = pCurGscImage->DecPM;

              pCurStarEntry->Stdmag = pCurGscImage->Stdmag;
              pCurStarEntry->color = pCurGscImage->color;
              pCurStarEntry->MAGFlag = pCurGscImage->MAGFlag;
              pCurStarEntry->class = pCurGscImage->class;
              pCurStarEntry->VFlag = pCurGscImage->VFlag;

              break;
            }
          }

          if (gscIndex == pCurStarIndex->numStars) {
            noGscEntryCount++;
          }

          memset(pTarget, 0, sizeof(PHOTTARGET));
          pTarget->Stdmag = pCurStarEntry2->Stdmag;
          pTarget->color = pCurStarEntry2->color;
          pTarget->ra = pCurStarEntry2->ra;
          pTarget->dec = pCurStarEntry2->dec;
          pTarget->MAGFlag = pCurStarEntry2->MAGFlag;
          pTarget->class = pCurStarEntry2->class;
          pTarget->VFlag = pCurStarEntry2->VFlag;
          pTarget->RaPM = pCurStarEntry2->RaPM;
          pTarget->DecPM = pCurStarEntry2->DecPM;
          pTarget->nearestCatalogStar.nearestREFflagGoodMag = -1;
          pTarget->nearestCatalogStar.nearestREFflag = -1;

          goodStars += ProcessLightcurveParameters(
            stdout,
            pFileCommon,
            pPhotConnection,
            pMagnitudeTable,
            curMagnitudes,
            pCurStarEntry2,
            pVectorStore,
            minGoodStars,
            pPhotGlobal->currentVersion,
            idHandle,
            &summaryCount,
            debugMode,
            &staleVersionIdCount,
            damonSeriesId,
            excludeSeriesId,
            excludeSeriesCount,
            pTarget,
            curGscBinIndex
          );

          if (summaryCount >= SUMMARY_LIMIT) {
            summaryCount = 0;

            time(&curTime);
            deltaTime = curTime - oldTime;
            oldTime = curTime;
            curTime -= startTime;

            printf(
              "Processing      %ld seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",
              curTime,
              GetQueryCount(),
              totalStars,
              goodStars,
              magnitudeAlloc,
              pVectorStore->vectorAlloc,
              binCount,
              totalMagnitudes
            );

            if (outHandle != NULL) {
              fprintf(
                outHandle,
                "Processing      %ld seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",
                curTime,
                GetQueryCount(),
                totalStars,
                goodStars,
                magnitudeAlloc,
                pVectorStore->vectorAlloc,
                binCount,
                totalMagnitudes
              );
            }

            if (deltaTime > 0) {
              printf(
                " queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec: %f magnitudes/sec %f\n",
                (1.0 * (GetQueryCount() - oldQueryCount)) / (1.0 * deltaTime),
                (1.0 * (totalStars - oldTotalStars)) / (1.0 * deltaTime),
                (1.0 * (goodStars - oldGoodStars)) / (1.0 * deltaTime),
                (1.0 * (binCount - oldBinCount)) / (1.0 * deltaTime),
                (1.0 * (totalMagnitudes - oldTotalMagnitudes)) / (1.0 * deltaTime)
              );

              if (outHandle != NULL) {
                fprintf(
                  outHandle,
                  " queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec: %f magnitudes/sec %f\n",
                  (1.0 * (GetQueryCount() - oldQueryCount)) / (1.0 * deltaTime),
                  (1.0 * (totalStars - oldTotalStars)) / (1.0 * deltaTime),
                  (1.0 * (goodStars - oldGoodStars)) / (1.0 * deltaTime),
                  (1.0 * (binCount - oldBinCount)) / (1.0 * deltaTime),
                  (1.0 * (totalMagnitudes - oldTotalMagnitudes)) / (1.0 * deltaTime)
                );
              }

              oldQueryCount = GetQueryCount();
              oldTotalStars = totalStars;
              oldGoodStars = goodStars;
              oldBinCount = binCount;
              oldTotalMagnitudes = totalMagnitudes;
            }

            ChargeTime(0, &curUsec, deltaTable, totalTable);

            printf("Intermediate Time (sec) ");
            if (outHandle != NULL) {
              fprintf(outHandle,"Intermediate Time (sec) ");
            }

            for (usecIndex = 0; usecIndex < MAX_CHARGE_ENTRY; usecIndex++) {
              double outSeconds;
              outSeconds = (1.0 * deltaTable[usecIndex]) / 1000000.0;
              printf("%10.3f ", outSeconds);
              if (outHandle != NULL) {
                fprintf(outHandle, "%10.3f ", outSeconds);
              }
            }

            printf("\n");
            if (outHandle != NULL) {
              fprintf(outHandle,"\n");
            }

            ChargeTime(RESET_DELTA_ENTRY, &curUsec, deltaTable, totalTable);

            if (outHandle != NULL) {
              fflush(outHandle);
            }
          }
        }

        if (debugMode) {
          time(&curTime);
          curTime -= startTime;
          fprintf(
            debug_handle,
            "Finished processing %s summaryCount %d starIndex %d starIndex2 %d, processedCount %d at %ld\n",
            pCurStarEntry2->REF,
            summaryCount,
            starIndex,
            starIndex2,
            processedCount,
            curTime
          );
        }
      }
    }

    curGscBinIndex++;
  }

  time(&curTime);
  curTime -= startTime;

  printf(
    "Execution Time: %ld seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",
    curTime,
    GetQueryCount(),
    totalStars,
    goodStars,
    magnitudeAlloc,
    pVectorStore->vectorAlloc,
    binCount,
    totalMagnitudes
  );

  if (outHandle != NULL) {
    fprintf(
      outHandle,
      "Execution Time: %ld seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",
      curTime,
      GetQueryCount(),
      totalStars,
      goodStars,
      magnitudeAlloc,
      pVectorStore->vectorAlloc,
      binCount,
      totalMagnitudes
    );
  }

  printf(
    " queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec %f magnitudes/sec %f\n",
    (1.0 * GetQueryCount()) / (1.0 * curTime),
    (1.0 * totalStars) / (1.0 * curTime),
    (1.0 * goodStars ) / (1.0 * curTime),
    (1.0 * binCount  ) / (1.0 * curTime),
    (1.0 * totalMagnitudes) / (1.0 * curTime)
  );

  if (outHandle != NULL) {
    fprintf(
      outHandle,
      " queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec %f magnitudes/sec %f\n",
      (1.0 * GetQueryCount()) / (1.0 * curTime),
      (1.0 * totalStars) / (1.0 * curTime),
      (1.0 * goodStars ) / (1.0 * curTime),
      (1.0 * binCount  ) / (1.0 * curTime),
      (1.0 * totalMagnitudes) / (1.0 * curTime)
    );
  }

  ChargeTime(0, &curUsec, deltaTable, totalTable);
  printf("Intermediate Time (sec) ");
  if (outHandle != NULL) {
    fprintf(outHandle, "Intermediate Time (sec) ");
  }

  for (usecIndex = 0; usecIndex < MAX_CHARGE_ENTRY; usecIndex++) {
    double outSeconds;
    outSeconds = (1.0 * totalTable[usecIndex]) / 1000000.0;
    printf("%10.3f ", outSeconds);
    if (outHandle != NULL) {
      fprintf(outHandle, "%10.3f ", outSeconds);
    }
  }

  printf("\n");
  if (outHandle != NULL) {
    fprintf(outHandle,"\n");
  }

  printf("... %d noGscEntry\n", noGscEntryCount);
  printf("... %d staleVersionId\n", staleVersionIdCount);
  if (outHandle != NULL) {
    fprintf(outHandle, "... %d noGscEntry\n", noGscEntryCount);
    fprintf(outHandle, "... %d staleVersionId\n", staleVersionIdCount);
  }

  mysql_close(pPhotConnection);
  if (pStarTable != NULL) {
    free(pStarTable);
  }

  if (pMagnitudeTable != NULL) {
    free(pMagnitudeTable);
  }

  FreeFileCommon(pFileCommon,1);
  FreeParameterVectorStore(pVectorStore);

  if (pGscImageTable != NULL) {
    free(pGscImageTable);
  }

  if (debug_handle != NULL) {
    fclose(debug_handle);
    debug_handle = NULL;
  }

  if (target_table != NULL) {
    free(target_table);
  }

  Close(catalogHandle);
  Close(indexHandle);

  if (outHandle != NULL) {
    fclose(outHandle);
  }

  if (idHandle != NULL) {
    fclose(idHandle);
  }

  return 0;
}


/* Mar 31, 2015 Edward J. Los - adapted from update_summary
 * Apr 13, 2015 Edward J. Los - add galactic latitude as an extra table entry.
 * May  6, 2015 Edward J. Los -  Add nearbyREFflag if a TC has a catalog star within FWHM/2 and the catalog object is at least 0.5 mag dimmer than the TC (modified Josh request of Tue 5/05/15 9:15 AM)
 *                               Also add releaseField, peakDays, peakYear, peakMag, peakSlope, peakSlopeRMS, and peakCount
 * Jun 15, 2015 Edward J. Los - Change peakSlopeRMS to peakRMS, add peakUpperCount, peakNumber, peakMaxMag, and peakMinMag
 * Jun 24, 2015 Edward J. Los - Add peakEventCount
 * Jun 29, 2015 Edward J. Los - Add peakDradRMS
 * Jul 13, 2015 Edward J. Los - Add peakEventCount2
 * Jul 17, 2015 Edward J. Los - Add peakEventCount3
 * Aug  1, 2015 Edward J. Los - Add -t parameter to enable transient flare searches
 * Aug 18, 2015 Edward J. Los - Add peakExcess = npoints outside of flare
 *                              Add "catalogDistance" histogram of distances of DASCH objects from catalog objects
 *                              Correct flare search for matched transient candidates
 * Sep 15, 2015 Edward J. Los - Rename "peakExcess" to "peakExtra"
 *                              Add peakCountNF, peakCountWF, peakOutside, peakOutsideNF, and peakOutsideWF
 *                              Add peakRA and peakDec with 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window
 *                              Add peakDefectCount for the number of lightcurve defects in the selected transient window.
 * Sep 24, 2015 Edward J. Los - Add peakNearbyDistance for the distance to nearby stars that are within 2 magnitudes of the brightest or average magnitude of the transient candidate
 * Oct  4, 2015 Edward J. Los   Add peakMultipleCount for the number of multiple exposure plates in the selected transient window
 * Oct 26, 2015 Edward J. Los - Add gsc_bin_index for the current gsc bin of the transient
 *                              Add peakDradRMS2 for the dradRMS2 values added in quadrature
 *                              Add -g qualifier for testing and repeat execution
 * Oct 29, 2015 Edward J. Los - Add peakDradRMS3 for the drad rms of points in the selected 90 day window relative to peakRA and peakDec
 *                              Do not generate the timestamp from the current time - have the user pick the string.
 * Dec 23, 2015 Edward J. Los - Correct Seek bug in reading the star catalog
 * Mar 21, 2017 Edward J. Los - Remover the transient search.  search_none now handles both matched and unmatched objects
 * Jul  1, 2017 Edward J. Los - Introduce -w to provide only the web summary counts
 * Mar  5, 2018 Edward J. Los - Add peakBadColorCount, peakLimitingYears, peakLimitingPoints, peakEvaluation, and peakNoDefectMag to the id table
 * Apr  2, 2018 Edward J. Los - Add peakLongOutburst to the id table
 * May 29, 2018 Edward J. Los - support gaia
 * Aug 13, 2018 Edward J. Los   Define enableRematch to optimize location of transients
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 */
