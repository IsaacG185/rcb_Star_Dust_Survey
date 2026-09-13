// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* Step through the entire photometry magnitudes database to resort the files
 * and get an accurate count of the magnitudes
 *
 * Modified in-place:
 *
 * - $DASCH_PHOT_ROOT/{refcat}/magNNN/magNNN/magNNN/magNNNNNNNNN.db
 *
 * Outputs:
 *
 * - Logfile associated with `-o` argument
 * - Bin histogram associated with `-g` argument
 * - Plate list file in `-p` argument
 *
 * Database updates:
 *
 * - ???
 *
 * Environment variables:
 *
 * - DASCH_CATALOG
 * - DASCH_PHOT_ROOT
 * - DASCH_MYSQLHOST
 * - DASCH_PASSWORD
 * - DASCH_USERNAME
 * - SKY2K_PATH
 */

#include <math.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>

#include "mysql.h"

#include "table.h"

#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "libwcs/wcscat.h"

#include "scandb.h"
#include "photometryutils.h"
#include "pipelineutils.h"
#include "searchgsc.h"

#define MAX_BUFFER 256
#define VECTOR_ALLOC_INCREMENT 500
#if 0
#define MAGNITUDE_LIMIT 2000
#else
#define MAGNITUDE_LIMIT 50000
#endif
#define RESET_DELTA_ENTRY -1
#define RESET_ALL_ENTRY -2
/* #define LOS_DEBUG 1 */
/* #define LOS_GSC_BIN 112757560 */
#define MIN_YEAR 1880
#define MAX_YEAR 2000
#define TOTAL_YEARS (MAX_YEAR-MIN_YEAR)
#define SEQUENCE_PATROL     1
#define SEQUENCE_NONPATROL  2
#define SEQUENCE_TOTAL      3

#define MAX_MAG 20
#define LIMITING_FACTOR 100L

/* The following for N03102131317 at 311.22600006 7.72723982 raPM = 4.6 decPM = -9.7
   median measured 311.229      7.72543
   nearby object   311.234      7.71949

   It shows the malmquist effect
   unknown object                    311.228      7.7252
   RA -.001 degrees in 70 years raPM = -51  mas/yr  (unmatched  slightly positive when matched)
   DEC +0.0004 degrees in 70 years  = +20  mas/yr  (unmatched, slightly negative when matched)
   unknown object
*/
#if 0
#define START_GSC_BIN_INDEX 138240000
#define END_GSC_BIN_INDEX   143360000
#endif

#define FILE_PRINT_MODULUS 5000
#define MIN_OUTLIER_SEARCH 20 /* Minimum number of stars for an outlier search */

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
/*
 *  Search the spatial bin for good matched stars with good points.  If one is found with more than 10 good points, then look for outliers
 *
 */
#define VECTOR_MAGCAL_MAGDEP1 0
#define VECTOR_MAGCAL_MAGDEP2 1  /* Necessary because the median calculation sorts the vectors */
#define MAX_VECTOR          2


#define VECTOR_PROCESSED 0
#define VECTOR_INDEX     1
#define MAX_IVECTOR      2

typedef struct _vectorstorage {
  int vectorAlloc;   /* Size of the vectors */
  int outlierCount;
  int lightcurvePointCount;
  int lightcurveCount;
  int qualityRejectCount;
  int pointsRead2;
  int pointsRead;
  int pointsRejected;
  int insufficientPoints;
  int highRMSPoints;

  double     *vector[MAX_VECTOR];
  int        *ivector[MAX_IVECTOR];
  int *seriesLightcurveCountTable[MAX_SERIES+1];
  int *seriesOutlierCountTable[MAX_SERIES+1];


} VECTORSTORE,*PVECTORSTORE;

long long  binHistogram[BIN_HISTOGRAM_SIZE];

void ReallocVectorStore(PVECTORSTORE pVectorStore)
{
  int index;
  double *tmpVector;
  int *tmpivector;
  for (index = 0; index < MAX_VECTOR; index++) {
    tmpVector = realloc(pVectorStore->vector[index],pVectorStore->vectorAlloc * sizeof(double));
    if (tmpVector == NULL) {
      fprintf(stderr,"ERROR: Failed to reallocate vector %d of size %zu\n",index,pVectorStore->vectorAlloc*sizeof(double));
      exit(-1);
    }
    pVectorStore->vector[index] = tmpVector;

  }

  for (index = 0; index < MAX_IVECTOR; index++) {
    tmpivector = realloc(pVectorStore->ivector[index],pVectorStore->vectorAlloc * sizeof(int));
    if (tmpivector == NULL) {
      fprintf(stderr,"ERROR: Failed to reallocate ivector %d of size %zu\n",index,pVectorStore->vectorAlloc*sizeof(double));
      exit(-1);
    }
    pVectorStore->ivector[index] = tmpivector;

  }

}

void FreeVectorStore(PVECTORSTORE pVectorStore)
{
  int index;
  double *vector;
  int * ivector;

  for (index = 0; index < MAX_VECTOR; index++) {
    vector = pVectorStore->vector[index];
    if (vector != NULL) {
      free(vector);
    }
  }
  for (index = 0; index < MAX_IVECTOR; index++) {
    ivector = pVectorStore->ivector[index];
    if (ivector != NULL) {
      free(ivector);
    }
  }
  for (index = 0; index <= MAX_SERIES; index++) {
    ivector = pVectorStore->seriesLightcurveCountTable[index];
    if (ivector != NULL) {
      free(ivector);
    }
    ivector = pVectorStore->seriesOutlierCountTable[index];
    if (ivector != NULL) {
      free(ivector);
    }
  }

  memset(pVectorStore,0,sizeof(VECTORSTORE));
}

int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int verbose = 0;
  int publicReleaseFlag = 0;
  int releaseField;
  int authorizedBinCount = 0;
  int totalBinCount = 0;
  double centerRA;
  double centerDec;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char platefile[MAX_BUFFER];
  char binhistogramfile[MAX_BUFFER];
  char magcountvsyear[MAX_BUFFER];
  char magcountvsmag[MAX_BUFFER];
  char limvsyear[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int nvals;
  FILE *outHandle = NULL;
  FILE *magcountvsyearHandle = NULL;
  FILE *magcountvsmagHandle = NULL;
  FILE *limvsyearHandle = NULL;
  FILE *plateHandle = NULL;
  FILE *binhistogramHandle = NULL;
  char timestr[100];
  struct tm *ptr;
  int gotAnswer;
  int currentIndex;
  int curMagnitudes = 0;
  int totalProcessedMagnitudes = 0;
  int totalNoneMagnitudes = 0;
  int totalBlendNomatchMagnitudes = 0;
  int staleVersionIdCount = 0;
  char *dotPtr;
  int index;

  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;

  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;

  PPHOTSTARIMAGE pMagnitudeTable;
  PFILESTARIMAGE pNoneMagnitudeTable = NULL;

  int noneMagnitudeAlloc = 0;
  int recoveredStars = 0;
  int badModulusBin = 0;
  int binCount = 0;
  int noGscEntryCount = 0;
  int totalStars = 0;

  int res;
  char queryString[MAX_QUERY_STRING];
  int magnitudeIndex;
#ifdef START_GSC_BIN_INDEX
  int curGscBinIndex = START_GSC_BIN_INDEX;
#else /* START_GSC_BIN_INDEX */
  int curGscBinIndex = 0;
#endif /* START_GSC_BIN_INDEX */
  char indexname[MAX_BUFFER];
  char* catalogdir;
  char* slashPtr;
  char catalogname[MAX_BUFFER];

  SetQueryCount(0);

  int keepaliveFlag = 1;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  int filePrintModulus = 0;
  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon = &fileCommon;
  PFILESTARIMAGE pFileStarImage = NULL;
  long long curPrintMagnitudes = 0;
  long long totalMagnitudes = 0;
  int histogramTable[MAX_NONE_HISTOGRAM];
  int groupNumber = 1;
  PSKY2KSTAR Sky2KTable = NULL;
  int sky2k_size = 0;
  double maxGSCRadius = 0.0;
  VECTORSTORE vectorStore;
  PVECTORSTORE pVectorStore = &vectorStore;
  int goodPlateCount = 0;
  int *seriesVersionArray;
  PTIMEACCURACYENTRY seriesTimeAccuracyArray;
  PTIMEACCURACYENTRY seriesTimeAccuracyEntry;
  double timeAccuracy;
  int seriesId;
  int plateNumber;
  char* series;
  int magFileModulus = MAG_FILE_MODULUS;
  int dateyear;
  long long patrolYearHistogram[TOTAL_YEARS+1];
  long long nonPatrolYearHistogram[TOTAL_YEARS+1];
  long long patrolMagHistogram[MAX_MAG+1];
  long long nonPatrolMagHistogram[MAX_MAG+1];
  long long patrolMagCount = 0;
  long long nonPatrolMagCount = 0;
  long long patrolImageCount = 0;
  long long nonPatrolImageCount = 0;
  long long limitingSum;
  long long patrolLimitingHistogram[TOTAL_YEARS+1];
  long long nonPatrolLimitingHistogram[TOTAL_YEARS+1];
  long long patrolLimitingCount[TOTAL_YEARS+1];
  long long nonPatrolLimitingCount[TOTAL_YEARS+1];
  int limitingMagnitude;
  double aveLimitingMagnitude;
  PSERIESENTRY pSeriesEntry;
  long long totalSequesteredMagnitudes = 0;
  long long matchedMagnitudes = 0;
  long long unMatchedMagnitudes = 0;
  long long zerodatecount = 0;
#ifdef START_GSC_BIN_INDEX
  printf("ERROR: START_GSC_BIN_INDEX is set\n");
#endif /* START_GSC_BIN_INDEX */
#ifdef END_GSC_BIN_INDEX
  printf("ERROR: END_GSC_BIN_INDEX is set\n");
#endif /* END_GSC_BIN_INDEX */
  long long histogramIndex;


  memset(pVectorStore,0,sizeof(VECTORSTORE));
  memset(patrolYearHistogram,0,sizeof(patrolYearHistogram));
  memset(nonPatrolYearHistogram,0,sizeof(patrolYearHistogram));
  memset(patrolMagHistogram,0,sizeof(patrolMagHistogram));
  memset(nonPatrolMagHistogram,0,sizeof(patrolMagHistogram));


  memset(patrolLimitingHistogram,0,sizeof(patrolLimitingHistogram));
  memset(nonPatrolLimitingHistogram,0,sizeof(patrolLimitingHistogram));

  memset(patrolLimitingCount,0,sizeof(patrolLimitingCount));
  memset(nonPatrolLimitingCount,0,sizeof(patrolLimitingCount));


  GetSky2kCatalog(stdout,&Sky2KTable,&sky2k_size,&maxGSCRadius);

#if 0 /* Testing of DASCH reference designation support */
  {
    double dec;
    double ra;
    long long REFNumber;
    int result;
    int refType;
    char REF[MAX_REF];
    dec = 14.120293;
    ra = 127.660377;
    result = GetDASCHNumber(ra,dec,&REFNumber,1,REF_TYPE_DASCH);
    printf("result is %d for ra %f dec %f, REFNumber %lld\n",result,ra,dec,REFNumber);
    GetREF(REFNumber,REF,0,1);
    printf("REF is %s for REFNumber %lld\n",REF,REFNumber);
    result = GetREFNumber(REF,&REFNumber,&refType,1);
    printf("result %d for %s is REFNumber %lld\n",result,REF,REFNumber);
    exit(-1);

  }


#endif

#if 0 /* Testing of FindAdjacentBins. Takes 1054 seconds with -O2 */
  {
    int gsc_bin;
    int bin_count;
    int bin_index;
    int max_bin_count = 0;
    int min_bin_count;
    int ave_bin_count;
    int bin_list[MAX_ADJACENT_BINS];
    time(&startTime);
    for (gsc_bin = 0; gsc_bin < pGscBin->total_gsc_bins; gsc_bin++) {
      FindAdjacentBins(pGscBin,gsc_bin,bin_list,&bin_count);
      if (bin_count > max_bin_count) {
        max_bin_count = bin_count;
      }
      if (gsc_bin == 0) {
        min_bin_count = bin_count;
      } else if (bin_count < min_bin_count) {
        min_bin_count = bin_count;
      }
      ave_bin_count += bin_count;
      if ((gsc_bin % 10000000) == 0) {
        time(&curTime);
        curTime -= startTime;
        printf("%10d of %10d bins min_bin_count %d ave_bin_count %f max_bin_count is %d in %d seconds\n",gsc_bin,pGscBin->total_gsc_bins,min_bin_count,(1.0*ave_bin_count)/(gsc_bin+1),max_bin_count,curTime);
        for (bin_index = 0; bin_index < bin_count; bin_index++) {
          printf("%d ",bin_list[bin_index]);
        }
        printf("\n");
      }
    }
    time(&curTime);
    curTime -= startTime;
    printf("min_bin_count %d ave_bin_count %f max_bin_count is %d in %d seconds\n",min_bin_count,(1.0*ave_bin_count)/gsc_bin,max_bin_count,curTime);
  }
  exit(-1);
#endif


  catalogString[0] = 0;
  memset(histogramTable,0,sizeof(histogramTable));


  outfile[0] = 0;
  binhistogramfile[0] = 0;
  platefile[0] = 0;
  qualifier[0] = 0;
  magcountvsyear[0] = 0;
  magcountvsmag[0] = 0;
  limvsyear[0] = 0;

  noneMagnitudeAlloc = 100;
  pNoneMagnitudeTable = (PFILESTARIMAGE)calloc(noneMagnitudeAlloc,sizeof(FILESTARIMAGE));

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d\n",argstr,argc);
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
            if (strlen(outfile) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'p': /* plate file name */
        case 'P':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(platefile,*++argv,MAX_BUFFER-2);
            if (strlen(platefile) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'g': /* declination bin histogram name */
        case 'G':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(binhistogramfile,*++argv,MAX_BUFFER-2);
            if (strlen(binhistogramfile) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          InitBinIndex(pGscBin);
          break;


        case 'y': /* magnitude count vs year file*/
        case 'Y':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(magcountvsyear,*++argv,MAX_BUFFER-2);
            if (strlen(magcountvsyear) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'm': /* magnitude count vs magnitude */
        case 'M':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(magcountvsmag,*++argv,MAX_BUFFER-2);
            if (strlen(magcountvsmag) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'l': /* limiting magnitude vs year */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(limvsyear,*++argv,MAX_BUFFER-2);
            if (strlen(limvsyear) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'i':  /* Ignore the versionId */
        case 'I':
          break;

        case 's': /* Starting GSC bin index */
        case 'S':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&curGscBinIndex);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode initial gsc bin index %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'r': /* public release */
        case 'R':
          publicReleaseFlag = 1;
          magFileModulus = 1;
          break;

        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              fprintf(stderr,"ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              sprintf(catalogString,"%d",catalogNumber);
              strcpy(qualifier,*argv);
            }
          }
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
    fprintf(stderr,"ERROR: No output filename was specified\n");
    errorFlag = 1;
  }

  outHandle = fopen(outfile,"a+t");
  if (outHandle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
  } else {
    if (verbose) {
      fprintf(stderr,"Output file %s\n",outfile);
    }
  }

  if (platefile[0] != 0) {
    plateHandle = fopen(platefile,"wt");
    if (plateHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the plate file %s\n",platefile);
    } else {

      if (verbose) {
        fprintf(stderr,"Plate file %s\n",platefile);
      }
    }
  }


  if (magcountvsyear[0] != 0) {
    magcountvsyearHandle = fopen(magcountvsyear,"wt");
    if (magcountvsyearHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the magnitude count vs year file %s\n",magcountvsyear);
    } else {
      if (verbose) {
        fprintf(stderr,"Magnitude count vs year file %s\n",magcountvsyear);
      }
    }
  }

  if (magcountvsmag[0] != 0) {
    magcountvsmagHandle = fopen(magcountvsmag,"wt");
    if (magcountvsmagHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the magnitude count vs magnitude file %s\n",magcountvsmag);
    } else {
      if (verbose) {
        fprintf(stderr,"magnitude count vs magnitude file %s\n",magcountvsmag);
      }
    }
  }

  if (limvsyear[0] != 0) {
    limvsyearHandle = fopen(limvsyear,"wt");
    if (limvsyearHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the limiting magnitude vs year file %s\n",limvsyear);
    } else {
      if (verbose) {
        fprintf(stderr,"limiting magnitude vs year file %s\n",limvsyear);
      }
    }
  }



  catalogdir = getenv("DASCH_CATALOG");
  if (catalogdir == NULL) {
    fprintf(stderr,"DASCH_CATALOG is not defined\n");
    printf("0\n");
    return(-1);
  }
  strcpy(catalogname,catalogdir);
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



  strcpy(indexname,catalogname);
  dotPtr = strstr(indexname,".");
  if (dotPtr != NULL) {
    *dotPtr = 0;
  }
  strcat(indexname,".idx");



  if (errorFlag) {
    printf("Usage: resort_magfiles options\n");
    printf("  options: -v verbose\n");
    printf("           -i ignore version id\n");
    printf("           -r do publicly released plates only\n");
    printf("           -q <catalog>\n");
    printf("           -c <minimum image count for magnitude file>\n");
    printf("           -s <starting gsc bin index>\n");
    printf("           -o <output file>\n");
    printf("           -d <database file>\n");
    printf("           -p <list of plates found>\n");
    printf("           -y <table of magnitude measurements vs year>\n");
    printf("           -m <table of magnitude measurements vs magnitude>\n");
    printf("           -l <table of limiting magnitudes vs year\n");
    printf("           -g <bin histogram file>\n");

    return(-1);
  }

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pConnection);
  dasch_init_photdb(pPhotConnection);

  InitSeriesTable(pConnection,pPhotConnection);

  if (GetPhotometryGlobal(pPhotConnection,pPhotGlobal) != 1) {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }
  if (pPhotGlobal->magnitudeFile == PHOT_MAGNITUDEFILE_NO ) {
    printf("ERROR: magnitudes stored in a MySQL table is no longer supported\n");
    exit(-1);
  }

  InitFileCommon(stdout,pFileCommon,pConnection,pPhotConnection,catalogNumber,1);
  InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);
  InitCatalogAccess(pFileCommon);

  TimeStampPhotometryGlobal(pPhotConnection);
#if 0
  histogramIndex = BIN_HISTOGRAM_SIZE;
  histogramIndex *= pGscBin->total_gsc_bins;
  printf("resort_magfiles BIN_HISTOGRAM_SIZE*pGscBin->total_gsc_bins %lld\n",histogramIndex);
#endif
  printf("resort_magfiles of %s %s, outfile %s platefile %s magnitudeLimit %d PHOTSTARIMAGE: %zu curGscBinIndex %d minVersionId%d %d publicRelease %d release level %d\n",
         __DATE__,__TIME__,outfile,platefile,MAGNITUDE_LIMIT,sizeof(PHOTSTARIMAGE),curGscBinIndex,catalogNumber,pPhotGlobal->minVersionId[catalogNumber],publicReleaseFlag,RELEASE_LEVEL);
  fprintf(outHandle,"resort_magfiles of %s %s, outfile %s platefile %s magnitudeLimit %d PHOTSTARIMAGE %zu curGscBinIndex %d minVersionId%d %d publicRelease %d release level %d\n",
          __DATE__,__TIME__,outfile,platefile,MAGNITUDE_LIMIT,sizeof(PHOTSTARIMAGE),curGscBinIndex,catalogNumber,pPhotGlobal->minVersionId[catalogNumber],publicReleaseFlag,RELEASE_LEVEL);


  time(&startTime);
  ptr = localtime(&startTime);
  strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);

  printf("resort_magfiles starting at %s\n",timestr);
  fprintf(outHandle,"resort_magfiles starting at %s\n",timestr);
  fflush(outHandle);

  if (plateHandle != NULL) {
    pFileCommon->seriesLimitingMagnitudeArray = (float **)calloc(MAX_SERIES+1,sizeof(float *));
    if (pFileCommon->seriesLimitingMagnitudeArray == NULL) {
      printf("ERROR: failed to allocate seriesLimitingMagnitudeArray of size %zu\n",(MAX_SERIES+1)*sizeof(float *));
    }
    fprintf(plateHandle,"Plate\tversionId\tcatalogNumber\tpatrolPlate\tlimiting_mag_local\ttimeAccuracy\n");
    fprintf(plateHandle,"-----\t---------\t-------------\t-----------\t------------------\t------------\n");
  }
  memset(binHistogram,0,sizeof(binHistogram));

  sprintf(queryString,"update photglobal set keepalive = 'yes';");
  res = ExecuteQuery(pPhotConnection,queryString);
  if (res) {
    exit(-1);
  }
  /* Enter the main loop */
  curMagnitudes = 0;
  pMagnitudeTable = NULL;
  magnitudeIndex = 0;
  curMagnitudes = -1;


  /* Enter the main loop.  We could process each magnitude file individually before going on to the
     next, but then the performance would be extremely slow because we lose the locality of
     reference into the GSC catalog file.  This more complicated code ensures that we read
     the GSC catalog file in sequence
  */

  while (1) {
    /* First see if we need to read in any data */
    if (keepaliveFlag == 0) {
      break;
    }
    if (curGscBinIndex >= pGscBin->total_gsc_bins) {
      /* We are done */
      break;
    }
    if ((++filePrintModulus % FILE_PRINT_MODULUS) == 0) {
      filePrintModulus = 0;


      gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
      if ((gotAnswer != 1)  || (pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES)) {

        time(&curTime);

        ptr = localtime(&curTime);
        strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);

        printf("resort_magfiles aborting at %s with gotAnswer %d keepalive %d\n",
               timestr,gotAnswer,pPhotGlobal->keepalive);
        fprintf(outHandle,"resort_magfiles aborting at %s with gotAnswer %d keepalive %d\n",
                timestr,gotAnswer,pPhotGlobal->keepalive);
        keepaliveFlag = 0;

        break;
      }
    }
    if ((filePrintModulus % FILE_PRINT_MODULUS) == 0) {
      time(&curTime);
      curTime -= startTime;
      printf("Read %10lld magnitudes for curGscBinIndex %10d at %5ld seconds\n",curPrintMagnitudes,curGscBinIndex,curTime);
      fprintf(outHandle,"Read %10lld magnitudes for curGscBinIndex %10d at %5ld seconds\n",curPrintMagnitudes,curGscBinIndex,curTime);
      fflush(outHandle);
      curPrintMagnitudes = 0;
    }

    totalBinCount++;
    if (publicReleaseFlag) {
      if (GetBinCenter(pGscBin,curGscBinIndex,&centerRA,&centerDec,"resort_magfiles") != 0) {
        printf("ERROR: GetBinCenter failure\n");
        exit(-1);
      }
      if (CheckAuthorization("none",centerRA,centerDec,&releaseField) == 0) {
        curGscBinIndex += magFileModulus;
#ifdef END_GSC_BIN_INDEX
        if (curGscBinIndex > END_GSC_BIN_INDEX) {
          curGscBinIndex = pGscBin->total_gsc_bins+1;
        }
#endif /* END_GSC_BIN_INDEX */
        if (curGscBinIndex >= pGscBin->total_gsc_bins) {
          printf("done curGscBinIndex %d exceeds maximum %d\n",curGscBinIndex,pGscBin->total_gsc_bins);
          break;
        }

        continue;
      }
      authorizedBinCount++;
    } else {
      authorizedBinCount++;
    }

    curMagnitudes = GetFileSummaryMagnitudes(pGscBin,pFileCommon,curGscBinIndex,curGscBinIndex+magFileModulus-1,0,0,0,0,NULL,catalogString,0,NULL,0,0);
    curPrintMagnitudes += curMagnitudes;

    if ((curMagnitudes != 0) && (binhistogramfile[0] != 0)) {
      histogramIndex = curGscBinIndex;
      histogramIndex *= BIN_HISTOGRAM_SIZE;
      histogramIndex = histogramIndex/pGscBin->total_gsc_bins;
      if (histogramIndex < 0) {
        histogramIndex = 0;
      }
      if (histogramIndex >= pGscBin->total_gsc_bins) {
        histogramIndex = pGscBin->total_gsc_bins;
      }
      binHistogram[histogramIndex] += curMagnitudes;
    }
#if 0
    if (curMagnitudes != 0) {
      printf("Have %d magnitudes at line %d\n",curMagnitudes,__LINE__);
    }
#endif

    if ((platefile[0] != 0) ||
        (magcountvsyear[0] != 0) ||
        (magcountvsmag[0] != 0) ||
        (limvsyear[0] != 0)) {
      for (currentIndex = 0; currentIndex < curMagnitudes; currentIndex++) {
#if 0
        if ((currentIndex % 10000) == 0) {
          printf("at %d in line %d\n",currentIndex,__LINE__);
        }
#endif
        pFileStarImage = &pFileCommon->magnitudeBuffer[currentIndex];
        pSeriesEntry = GetSeriesEntry(pFileStarImage->seriesId);
        if ((pSeriesEntry == NULL) ||
            (pSeriesEntry->sequestered == SEQUESTERED_YES)) {
          /* Skip sequestered series */
          totalSequesteredMagnitudes++;
          continue;
        }
        if (pFileStarImage->REFNumber == 0) {
          unMatchedMagnitudes++;
        } else {
          matchedMagnitudes++;
        }
        if (pFileStarImage->Date == 0) {
          zerodatecount++;
          dateyear = TOTAL_YEARS;
        } else {
          dateyear = jd2ep(pFileStarImage->Date) - MIN_YEAR;
          if ((dateyear < 0) || (dateyear >= MAX_YEAR)) {
          dateyear = TOTAL_YEARS;
          }
        }
        magnitudeIndex = pFileStarImage->magcal_magdep;
        if ((magnitudeIndex < 1) || (magnitudeIndex > MAX_MAG)) {
          magnitudeIndex = MAX_MAG;
        }
        limitingMagnitude = pFileStarImage->limiting_mag_local;
        if ((limitingMagnitude < 1) || (limitingMagnitude > MAX_MAG)) {
          limitingMagnitude = MAX_MAG;
        }


        if (GetFittedPlateScale(pFileStarImage->seriesId,pFileStarImage->plateNumber) >= PATROL_PLATE_SCALE) {
          /* Patrol plate counters */
          patrolYearHistogram[dateyear]++;
          patrolImageCount++;
          patrolMagHistogram[magnitudeIndex]++;
          if (magnitudeIndex != MAX_MAG) {
            patrolMagCount++;
          }
          if (limitingMagnitude != MAX_MAG) {
            limitingSum = (pFileStarImage->limiting_mag_local*LIMITING_FACTOR); /* A direct addition would peg at 34359738368 */
            patrolLimitingHistogram[dateyear] += limitingSum;
            patrolLimitingCount[dateyear]++;
#if 0
            /* NOTE: cut down the size with: grep "XXXX" limvsyear_apass_DR4.db | awk '((NR%10000)==0){print $0}' > los.db */
            if (dateyear == (1940-MIN_YEAR)) {
              if ((debugPrintCounter++ % 10000) == 0) {
                fprintf(limvsyearHandle,"XXXX%d\t%d\t%d\t%f\t%lld\t%lld\t%lld\t%lld\n",pFileStarImage->seriesId,pFileStarImage->plateNumber,pFileStarImage->NUMBER,pFileStarImage->limiting_mag_local,pFileStarImage->REFNumber,patrolLimitingHistogram[dateyear],patrolLimitingCount[dateyear],patrolLimitingHistogram[dateyear]/patrolLimitingCount[dateyear]);
              }
            }
#endif
          }


        } else {
          /* Non-patrol plate counters */
          nonPatrolYearHistogram[dateyear]++;
          nonPatrolImageCount++;
          nonPatrolMagHistogram[magnitudeIndex]++;
          if (magnitudeIndex != MAX_MAG) {
            nonPatrolMagCount++;
          }
          if (limitingMagnitude != MAX_MAG) {
            limitingSum = (pFileStarImage->limiting_mag_local*LIMITING_FACTOR);
            nonPatrolLimitingHistogram[dateyear] += limitingSum;
            nonPatrolLimitingCount[dateyear]++;
          }
        }

        if (RELEASE_EXPERIMENTAL == 0 || catalogNumber != CATALOG_EXPERIMENTAL) {
          GetLatestVersionId(pFileCommon, pFileStarImage->seriesId, pFileStarImage->plateNumber, 0, pFileStarImage->REFNumber, 0);
        }
      }
    }


    totalMagnitudes += curMagnitudes;


#ifdef LOS_DEBUG
    {
      int losIndex;
      for (losIndex = 0; losIndex < curMagnitudes; losIndex++) {
        pCurStarImage = &pMagnitudeTable[losIndex];
        pFileStarImage = pCurStarImage->pFileStarImage;
        printf("DEBUG: index %d ra %f dec %f %s%05d gsc_bin_index %d\n",losIndex,pFileStarImage->ra,pFileStarImage->dec,pCurStarImage->series,pFileStarImage->plateNumber,pFileStarImage->gsc_bin_index);

      }
    }


#endif /* LOS_DEBUG */
    if (keepaliveFlag == 0) {
      break;
    }

    /* On to the next image */
    curGscBinIndex += magFileModulus;
#ifdef END_GSC_BIN_INDEX
    if (curGscBinIndex > END_GSC_BIN_INDEX) {
      curGscBinIndex = pGscBin->total_gsc_bins+1;
    }
#endif /* END_GSC_BIN_INDEX */
    if (curGscBinIndex >= pGscBin->total_gsc_bins) {
      printf("done curGscBinIndex %d exceeds maximum %d\n",curGscBinIndex,pGscBin->total_gsc_bins);
      break;
    }
  }


  time(&curTime);
  curTime -= startTime;
  if (platefile[0] != 0) {
    for (seriesId = 0; seriesId < MAX_SERIES; seriesId++) {
      seriesVersionArray = pFileCommon->seriesVersionTable[seriesId];
      if (seriesVersionArray == NULL) {
        continue;
      }
      seriesTimeAccuracyArray = pFileCommon->seriesTimeAccuracyTable[seriesId];
      series = GetSeriesString(seriesId,1);
      for (plateNumber = 0; plateNumber <= pFileCommon->maxPlateNumber[seriesId]; plateNumber++) {
        int patrolPlate;
        float limiting_mag_local;
        timeAccuracy = -1;
        if (seriesVersionArray[plateNumber] != 0) {
          float **seriesLimitingMagnitudeArray = pFileCommon->seriesLimitingMagnitudeArray;
          float *seriesLimitingMagnitudeTable;
          if (GetFittedPlateScale(seriesId,plateNumber) >= PATROL_PLATE_SCALE) {
            patrolPlate = 1;
          } else {
            patrolPlate = 0;
          }
          if (seriesLimitingMagnitudeArray != NULL) {
            seriesLimitingMagnitudeTable = seriesLimitingMagnitudeArray[seriesId];
            if (seriesLimitingMagnitudeTable != NULL) {
              limiting_mag_local = seriesLimitingMagnitudeTable[plateNumber];
            } else {
              limiting_mag_local = 99.0;
            }
          } else {
            limiting_mag_local = 99.0;
          }
          if (seriesTimeAccuracyArray != NULL) {
            seriesTimeAccuracyEntry = &seriesTimeAccuracyArray[plateNumber];
            if (seriesTimeAccuracyEntry->valid != 0) {
              timeAccuracy = seriesTimeAccuracyEntry->returnedTolerance;
            }
          }

          fprintf(plateHandle,"%s%05d\t%d\t%d\t%d\t%0.2f\t%f\n",series,plateNumber,seriesVersionArray[plateNumber],catalogNumber,patrolPlate,limiting_mag_local,timeAccuracy);
          goodPlateCount++;

        } else {
          if (seriesTimeAccuracyArray != NULL) {
            seriesTimeAccuracyEntry = &seriesTimeAccuracyArray[plateNumber];
            if (seriesTimeAccuracyEntry->valid != 0) {
              timeAccuracy = seriesTimeAccuracyEntry->returnedTolerance;
              patrolPlate = -1;
              limiting_mag_local = 99.0;
              fprintf(plateHandle,"%s%05d\t%d\t%d\t%d\t%0.2f\t%f\n",series,plateNumber,seriesVersionArray[plateNumber],catalogNumber,patrolPlate,limiting_mag_local,timeAccuracy);

            }

          }
        }
      }
    }
    printf("Found %d good plates\n",goodPlateCount);
    fprintf(outHandle,"Found %d good plates\n",goodPlateCount);
  }



  if (binhistogramfile[0] != 0) {
    binhistogramHandle = fopen(binhistogramfile,"wt");
    if (binhistogramHandle == NULL) {
      fprintf(stderr,"ERROR: Failed to open the binhistogram file %s\n",binhistogramfile);
    } else {
      fprintf(stderr,"Binhistogram file %s\n",binhistogramfile);
    }
  }

  if (binhistogramHandle != NULL) {
    fprintf(binhistogramHandle,"gsc_bin_index\tcount2\n");
    fprintf(binhistogramHandle,"-------------\t------\n");
    for (histogramIndex = 0; histogramIndex < BIN_HISTOGRAM_SIZE; histogramIndex++) {
      if (binHistogram[histogramIndex] != 0) {
        fprintf(binhistogramHandle,"%lld\t%lld\n",(histogramIndex * (pGscBin->total_gsc_bins/BIN_HISTOGRAM_SIZE)),binHistogram[histogramIndex]);
      }
    }

  }



  printf("Execution Time: %ld seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d gscBins, %lld total magnitudes %lld currentMagnitudes\n",curTime,GetQueryCount(),totalStars,recoveredStars,badModulusBin,binCount,totalMagnitudes,pFileCommon->counterBlock.currentMagnitudes);

  fprintf(outHandle,"Execution Time: %ld seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d gscBins, %lld total magnitudes %lld currentMagnitudes\n",curTime,GetQueryCount(),totalStars,recoveredStars,badModulusBin,binCount,totalMagnitudes,pFileCommon->counterBlock.currentMagnitudes);
  printf(" queries/sec: %f totalStars/sec: %f recoveredStars+badModulusBin.sec: %f gscBins/sec %f\n",
         (1.0 * GetQueryCount())/(1.0*curTime),
         (1.0 * totalStars)    /(1.0*curTime),
         (1.0 * (recoveredStars + badModulusBin) )    /(1.0*curTime),
         (1.0 * binCount  )    /(1.0*curTime));
  fprintf(outHandle," queries/sec: %f totalStars/sec: %f recoveredStars+badModulusBin.sec: %f gscBins/sec %f\n",
          (1.0 * GetQueryCount())/(1.0*curTime),
          (1.0 * totalStars)    /(1.0*curTime),
          (1.0 * (recoveredStars+badModulusBin) )    /(1.0*curTime),
          (1.0 * binCount  )    /(1.0*curTime));



  printf("... %d noGscEntry\n",noGscEntryCount);
  printf("... %d staleVersionId\n",staleVersionIdCount);
  fprintf(outHandle,"... %d noGscEntry\n",noGscEntryCount);
  fprintf(outHandle,"... %d staleVersionId\n",staleVersionIdCount);

  printf("New curGscBinIndex %d\n",curGscBinIndex);
  fprintf(outHandle,"New curGscBinIndex %d\n",curGscBinIndex);

  printf("Total Magnitudes %d, Total NONE magnitudes %d Total NOMATCH magnitudes %d Total GSC nonblend magnitudes %d \n",
         totalProcessedMagnitudes,totalNoneMagnitudes,totalBlendNomatchMagnitudes,totalProcessedMagnitudes - totalNoneMagnitudes - totalBlendNomatchMagnitudes);
  fprintf(outHandle,"Total Magnitudes %d, Total NONE magnitudes %d Total NOMATCH magnitudes %d Total GSC nonblend magnitudes %d \n",
          totalProcessedMagnitudes,totalNoneMagnitudes,totalBlendNomatchMagnitudes,totalProcessedMagnitudes - totalNoneMagnitudes - totalBlendNomatchMagnitudes);
  printf("Total groups %d\n",groupNumber-1);
  fprintf(outHandle,"Total groups %d\n",groupNumber-1);


  for (index = 0; index < MAX_NONE_HISTOGRAM; index++) {
    if (histogramTable[index] != 0) {
      printf("Images/Location: %4d  Locations: %4d\n",index,histogramTable[index]);
      fprintf(outHandle,"Images/Location %4d  Locations: %6d\n",index,histogramTable[index]);
    }
  }

  mysql_close(pConnection);
  mysql_close(pPhotConnection);
  if (pMagnitudeTable != NULL) {
    free(pMagnitudeTable);
  }
  if (pNoneMagnitudeTable != NULL) {
    free (pNoneMagnitudeTable);
  }

  if ((platefile[0] != 0) ||
      (magcountvsyear[0] != 0) ||
      (magcountvsmag[0] != 0) ||
      (limvsyear[0] != 0)) {
    printf("Matched magnitudes %lld  Unmatched Magnitudes %lld\n",matchedMagnitudes,unMatchedMagnitudes);
    printf("Valid Patrol images %lld Valid non-Patrol images %lld\n",patrolImageCount,nonPatrolImageCount);
    printf("Valid Patrol magnitudes %lld Valid non-Patrol magnitudes %lld\n",patrolMagCount,nonPatrolMagCount);
    if (zerodatecount > 0) {
      printf("ERROR: zerodatecount %lld\n",zerodatecount);
    }
    if (magcountvsyearHandle != NULL) {
      fprintf(magcountvsyearHandle,"sequence\tyear\tcatalogNumber\tcount\n");
      fprintf(magcountvsyearHandle,"--------\t----\t-------------\t-----\n");
      for (dateyear = 0; dateyear <= TOTAL_YEARS; dateyear++) {
        if (patrolYearHistogram[dateyear] != 0) {
          fprintf(magcountvsyearHandle,"%d\t%d\t%d\t%lld\n",SEQUENCE_PATROL,dateyear+MIN_YEAR,catalogNumber,patrolYearHistogram[dateyear]);
        }
        if (nonPatrolYearHistogram[dateyear] != 0) {
          fprintf(magcountvsyearHandle,"%d\t%d\t%d\t%lld\n",SEQUENCE_NONPATROL,dateyear+MIN_YEAR,catalogNumber,nonPatrolYearHistogram[dateyear]);
        }
        if ((patrolYearHistogram[dateyear] + nonPatrolYearHistogram[dateyear]) != 0) {
          fprintf(magcountvsyearHandle,"%d\t%d\t%d\t%lld\n",SEQUENCE_TOTAL,dateyear+MIN_YEAR,catalogNumber,nonPatrolYearHistogram[dateyear]+patrolYearHistogram[dateyear]);
        }

      }

    }
    if (magcountvsmagHandle != NULL) {
      fprintf(magcountvsmagHandle,"sequence\tmagcal_magdep\tcatalogNumber\tcount\n");
      fprintf(magcountvsmagHandle,"--------\t-------------\t-------\t-----\n");
      for (magnitudeIndex = 0; magnitudeIndex <= MAX_MAG; magnitudeIndex++) {
        if (patrolMagHistogram[magnitudeIndex] != 0) {
          fprintf(magcountvsmagHandle,"%d\t%d\t%d\t%lld\n",SEQUENCE_PATROL,magnitudeIndex,catalogNumber,patrolMagHistogram[magnitudeIndex]);
        }
        if (nonPatrolMagHistogram[magnitudeIndex] != 0) {
          fprintf(magcountvsmagHandle,"%d\t%d\t%d\t%lld\n",SEQUENCE_NONPATROL,magnitudeIndex,catalogNumber,nonPatrolMagHistogram[magnitudeIndex]);
        }
        if ((patrolMagHistogram[magnitudeIndex] + nonPatrolMagHistogram[magnitudeIndex]) != 0) {
          fprintf(magcountvsmagHandle,"%d\t%d\t%d\t%lld\n",SEQUENCE_TOTAL,magnitudeIndex,catalogNumber,nonPatrolMagHistogram[magnitudeIndex]+patrolMagHistogram[magnitudeIndex]);
        }

      }

    }
    if (limvsyearHandle != NULL) {
      fprintf(limvsyearHandle,"sequence\tyear\tcatalogNumber\tlimiting_mag_local\n");
      fprintf(limvsyearHandle,"--------\t----\t-------------\t------------------\n");
      for (dateyear = 0; dateyear <= TOTAL_YEARS; dateyear++) {
        if (patrolLimitingCount[dateyear] != 0) {

          aveLimitingMagnitude = (1.0*patrolLimitingHistogram[dateyear])/(1.0*patrolLimitingCount[dateyear]*LIMITING_FACTOR);
#if 0
          aveLongLimitingMagnitude = patrolLimitingHistogram[dateyear]/patrolLimitingCount[dateyear];
          fprintf(limvsyearHandle,"ZZZZ (1) %.2f %lld %lld %lld\n",aveLimitingMagnitude,patrolLimitingHistogram[dateyear],patrolLimitingCount[dateyear],aveLongLimitingMagnitude);
#endif
          fprintf(limvsyearHandle,"%d\t%d\t%d\t%.2f\n",SEQUENCE_PATROL,dateyear+MIN_YEAR,catalogNumber,aveLimitingMagnitude);
        }
        if (nonPatrolLimitingCount[dateyear] != 0) {
          aveLimitingMagnitude = (1.0*nonPatrolLimitingHistogram[dateyear])/(1.0*nonPatrolLimitingCount[dateyear]*LIMITING_FACTOR);
#if 0
          aveLongLimitingMagnitude = nonPatrolLimitingHistogram[dateyear]/nonPatrolLimitingCount[dateyear];
          fprintf(limvsyearHandle,"ZZZZ (2) %.2f %lld %lld %lld\n",aveLimitingMagnitude,nonPatrolLimitingHistogram[dateyear],nonPatrolLimitingCount[dateyear],aveLongLimitingMagnitude);
#endif
          fprintf(limvsyearHandle,"%d\t%d\t%d\t%.2f\n",SEQUENCE_NONPATROL,dateyear+MIN_YEAR,catalogNumber,aveLimitingMagnitude);
        }
        if ((patrolLimitingCount[dateyear] + nonPatrolLimitingCount[dateyear]) != 0) {
          aveLimitingMagnitude = (1.0*(patrolLimitingHistogram[dateyear]+nonPatrolLimitingHistogram[dateyear]))/(1.0*(patrolLimitingCount[dateyear] + nonPatrolLimitingCount[dateyear])*LIMITING_FACTOR);
#if 0
          aveLongLimitingMagnitude = (patrolLimitingHistogram[dateyear]+nonPatrolLimitingHistogram[dateyear])/(patrolLimitingCount[dateyear] + nonPatrolLimitingCount[dateyear]);
          fprintf(limvsyearHandle,"ZZZZ (3) %.2f %lld %lld %lld\n",aveLimitingMagnitude,patrolLimitingHistogram[dateyear]+nonPatrolLimitingHistogram[dateyear],patrolLimitingCount[dateyear] + nonPatrolLimitingCount[dateyear],aveLongLimitingMagnitude);
#endif
          fprintf(limvsyearHandle,"%d\t%d\t%d\t%.2f\n",SEQUENCE_TOTAL,dateyear+MIN_YEAR,catalogNumber,aveLimitingMagnitude);
        }

      }

    }
  }

  printf("Max vectorAlloc %d, lightcurveCount %d lightcurvePointCount %d, outlierCount %d, qualityRejectCount %d\n",
         pVectorStore->vectorAlloc,
         pVectorStore->lightcurveCount,
         pVectorStore->lightcurvePointCount,
         pVectorStore->outlierCount,
         pVectorStore->qualityRejectCount);
  printf("points read %d %d, pointsRejected %d, insufficientPoints %d, highRMSPoints %d\n",
         pVectorStore->pointsRead,
         pVectorStore->pointsRead2,
         pVectorStore->pointsRejected,
         pVectorStore->insufficientPoints,
         pVectorStore->highRMSPoints);
  if (publicReleaseFlag) {
    printf("total bins %d of %d, authorized bins %d\n",totalBinCount,pGscBin->total_gsc_bins,authorizedBinCount);
  }
  if (pFileCommon->seriesTimeAccuracyAlloc != 0) {
    printf("seriesTimeAccuracyAlloc %d, seriesTimeAccuracyCount %d\n",pFileCommon->seriesTimeAccuracyAlloc,pFileCommon->seriesTimeAccuracyCount);
  }

  FreeFileCommon(pFileCommon,1);
  FreeVectorStore(pVectorStore);

  fclose(outHandle);
  if (plateHandle != NULL) {
    fclose(plateHandle);
  }
  if (binhistogramHandle != NULL) {
    fclose(binhistogramHandle);
  }
  if (magcountvsyearHandle != NULL) {
    fclose(magcountvsyearHandle);
  }
  if (magcountvsmagHandle != NULL) {
    fclose(magcountvsmagHandle);
  }
  if (limvsyearHandle != NULL) {
    fclose(limvsyearHandle);
  }
  return(0);
}

/*
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   resort_magfiles.c pipelineutils.a -ltable -lutil  -lwcs -o resort_magfiles -L/usr/lib64/mysql -L/usr/lib/mysql -lmysqlclient
 *
 *  Nov  5, 2012 Edward J. Los - Adapted from search_none.c
 *  Jun 29, 2013 Edward J. Los - Add option to print a list of plates found
 *  Sep 17, 2013 Edward J. Los - Add statistical measurments for publicly released data.
 *  Apr 19, 2014 Edward J. Los - Obtain gsc bin star histogram for use in search_none.c
 *  Feb  2, 2015 Edward J. Los - Support update of the timeAccuracy value.
 *  May 29, 2018 Edward J. Los - Support gaia
 *  Sep 22, 2018 Edward J. Los - look for Date == 0 in the binary files
 *  Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 *
 *  gsc
 *   resort_magfiles  -o /dasch/Pipeline/ingest/resort_magfiles.log
 *
 *  Kepler field starts at GSC bin about 141505545
 *   resort_magfiles -q kepler  -o /dasch/Pipeline/ingest/resort_magfiles_kepler.log
 *
 *  apass
 *   resort_magfiles -q apass  -o /dasch/Pipeline/ingest/resort_magfiles_apass.log
 *
 *  experimental
 *   resort_magfiles -q experimental  -o /dasch/Pipeline/ingest/resort_magfiles_experimental.log
 *
 *  Statistics:
 resort_magfiles -r  -o /dasch/Pipeline/ingest/resort_magfiles_gsc.log -p /dasch/Pipeline/ingest/plates_gsc.txt -y /dasch/Pipeline/ingest/countvsyear_gsc.db -m /dasch/Pipeline/ingest/countvsmag_gsc.db -l /dasch/Pipeline/ingest/limvsyear_gsc.db

 resort_magfiles -r -q apass -o /dasch/Pipeline/ingest/resort_magfiles_apass.log -p /dasch/Pipeline/ingest/plates_apass.txt -y /dasch/Pipeline/ingest/countvsyear_apass.db -m /dasch/Pipeline/ingest/countvsmag_apass.db -l /dasch/Pipeline/ingest/limvsyear_apass.db
 *
 */
