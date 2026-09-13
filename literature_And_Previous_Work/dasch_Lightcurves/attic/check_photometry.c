// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* check_photometry.c
 *
 * This routine steps through the entire photometry magnitudes table while performing miscellaneous actions.  At the momemnt, we are
 * looking for gsc bin errors.
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   check_photometry.c pipelineutils.a -ltable -lutil  -lwcs -o check_photometry -L/usr/lib64/mysql -lmysqlclient
 *
 * May 19, 2009 Edward J. Los - Original Version
 * Jun 19, 2009 Edward J. Los - Add automatic repair
 * Jun 26, 2009 Edward J. Los - Update for multiple tables
 * Jul 11, 2009 Edward J. Los - Count magnitudes, and FILTER_AFLAG_BLEND_NOMATCH objects
 * Aug 24, 2009 Edward J. Los - Add kepler input catalog support
 * Dec  2, 2009 Edward J. Los - Add magnitude file support.
 * Mar 13, 2011 Edward J. Los - Add apass support
 * Nov  4, 2011 Edward J. Los - Add magnitude-dependent fields (magdep_bin, magcal_magdep, and magcal_magdep_rms)
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Sep 10, 2012 Edward J. Los - Reorganize code to optimize magnitude file support.
 * Mar  3, 2015 Edward J. Los - Look for zero timeAccuracy values
 * May 29, 2015 Edward J. Los - support the GAIA catalog
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 * Jan 26, 2020 Edward J. Los - Add '-p' to create a histogram of ISO<n> with FILTER_AFLAG_LIMITING_MAG set
 *
 *   check_photometry -r -v -d -i -o /dasch/Pipeline/ingest/check_photometry.log
 *   echo "check_photometry -v -r -i -o /dasch/Pipeline/ingest/check_photometry.log" | at now
 *
 *   search for timeAccuracy issues: check_photometry -q apass -s -i -o /dasch/Pipeline/ingest/apass_check_photometry.log 
 *   make an ISO0 histogram          check_photometry -q apass -p -i -o /dasch/Pipeline/ingest/apass_check_photometry.log 
 *         NOTE: ISO0 results in /backup/2020_01_23/atlas_check_photometry2020_01_27a.log and atlas_ISO0_size_and_grain_noise.png
 *    
 *   echo "check_photometry  -r  -v -d -i -o /dasch/Pipeline/ingest/check_photometry.log" | at now
 *   echo "check_photometry  -r -i -o /dasch/Pipeline/ingest/check_photometry.log" | at now
 *   echo "check_photometry  -q kepler -r -i -o /dasch/Pipeline/ingest/kepler_check_photometry.log" | at now
 *        340 seconds on dell with dell.list
 * WARNING: no gsc entry in bin 73331039 should be 73331039 for star S2002233116 ra 193.610156 dec -7.592091
 * WARNING: no gsc entry in bin 82491895 should be 82491895 for star S2031000177 ra 202.293584 dec -1.357030
 * REF	ra	dec	Stdmag	color	class	VFlag	MAGFlag	RaPM	DecPM
 * ---	--	---	------	-----	-----	-----	-------	----	-----
 * S2002233116	 193.61211169	  -7.61044510	 8.78	 0.31	0	0	10	  18.3	 -31.1
 * S2031000177	 202.31214768	  -1.36442529	 7.86	 1.32	0	0	10	 -42.6	 -58.9
 *
 *  S2002233116 found in bin 73285365
 *  S2031000177 found in bin 82468864 
 *
 *
 */


#include <math.h>
#include <errno.h>
#include "table.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include <sys/time.h>
#include "pipelineutils.h"
#include "mysql.h"
#include "photometryutils.h"
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
#define FILE_PRINT_MODULUS 5000
#define ISO_MAX_VALUE 15000
#define FILTER_AMASK_REJECT   ((1<<FILTER_AFLAG_CASEB) | (1<<FILTER_AFLAG_CASEC) | (1<<FILTER_AFLAG_CASED) | (1<<FILTER_AFLAG_BLEND) | (1<<FILTER_AFLAG_BLEND_NOMATCH) | (1<<FILTER_AFLAG_MULTIPLE_BLEND) | (1<<FILTER_AFLAG_DEFECT) | (1<<FILTER_AFLAG_LOW_ALTITUDE) | (1<<FILTER_AFLAG_WEDGE) | (1<<FILTER_AFLAG_BACKGROUND)  | (1<<FILTER_AFLAG_MULTIPLE_BLEND) | (1 << FILTER_AFLAG_DRADBIN) | (1 << FILTER_AFLAG_BIN9) | (1 << FILTER_AFLAG_ISO_RMS) | (1 << FILTER_AFLAG_LOCAL_RMS) | (1 << FILTER_AFLAG_HIZOUT))



extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

int binSearchCount = 0;
int cacheHits = 0;
int cacheMisses = 0;


/*
 * Given a GSC bin number and a position within the bin, this routine returns the
 * distance in degrees to the nearest edge.
 */
double GetGSCBorder(PGSCBIN pGscBin,int gsc_bin_index,double ra,double dec,char* outputname) 
{
  PBININDEX pBinIndex;
  double decval;
  double raval;
  int raBin;
  int decBin;
  double borderDist = 180.0;
  double tempDist;

  pBinIndex = GetSubbins(pGscBin,gsc_bin_index,&raBin,&decBin,outputname);
  raval = (360.0 * raBin)/(pBinIndex->numBins) + (180.0/(pBinIndex->numBins));
  decval = (decBin * pGscBin->bin_size) - 90.0 +  (pGscBin->bin_size/2.0);

  /* Bottom Edge */
  tempDist = fabs(dec - decval);
  if (tempDist < borderDist) {
    borderDist = tempDist;
  }

  /* Top Edge */
  decval += pGscBin->bin_size;
  tempDist = fabs(dec - decval);
  if (tempDist < borderDist) {
    borderDist = tempDist;
  }

  /* Left Edge */
  decval -= (pGscBin->bin_size/2.0);
  tempDist = fabs(ra - raval) * cos(decval*DEGREES_TO_RAD);
  if (tempDist > 180.0) {
    tempDist = 360.0-tempDist;
  }
  if (tempDist < borderDist) {
    borderDist = tempDist;
  }

  /* Right Edge */
  if (raBin == (pBinIndex->numBins-1)) {
    raBin = 0;
  } else {
    raBin++;
  }
  raval = (360.0 * raBin)/(pBinIndex->numBins) + (180.0/(pBinIndex->numBins));
  tempDist = fabs(ra - raval) * cos(decval*DEGREES_TO_RAD);
  if (tempDist > 180.0) {
    tempDist = 360.0-tempDist;
  }
  if (tempDist < borderDist) {
    borderDist = tempDist;
  }

  return(borderDist);


}


int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  time_t deltaTime;
  time_t oldTime;
  int suppressErrorMessages = 0;
  int limitingISOHistogram  = 0;
  int verbose = 0;
  int debugMode = 0;
  int ignoreVersionId = 0;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int nvals;
  FILE *outHandle;
  char timestr[100];
  struct tm *ptr;
  int gotAnswer;
  int curMagnitudes;
  int totalMagnitudes = 0;
  int totalProcessedMagnitudes = 0;
  int totalNoneMagnitudes = 0;
  int totalBlendNomatchMagnitudes = 0;
  int minGscBinIndex;
  int maxGscBinIndex;
  int staleVersionIdCount = 0;
  char *dotPtr;
  int processedCount = 0;
  int magnitudeFileFlag = 0;

  char *mysqlhost;
  char *username;
  char *password;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;

  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  STARENTRY curStarEntry;
  PSTARENTRY pCurStarEntry = &curStarEntry;

  PPHOTSTARIMAGE pMagnitudeTable;
  PPHOTSTARIMAGE pCurMagnitudeTable;
  PPHOTSTARIMAGE pCurStarImage;
  PFILESTARIMAGE pFileStarImage = NULL;


  int magnitudeAlloc;
  int totalMagnitudeAlloc = 0;
  int recoveredStars = 0;
  int badModulusBin = 0;
  int binCount = 0;
  int noGscEntryCount = 0;
  int repairMagnitudes = 0;

  FILE *stars_handle = NULL;
  FILE *summary_handle = NULL;
  char stars_name[MAX_BUFFER];
  char summary_name[MAX_BUFFER];
  char debug_name[MAX_BUFFER];
  FILE *debug_handle = NULL;
  char *ingestDirectory;
 
  int curStars;
  int totalStars = 0;
  int starIndex;
  int starIndex2;
  int summaryCount = 0;

  int res;
  char queryString[MAX_QUERY_STRING];
  int querySize;
  int oldQueryCount = 0;
  int oldTotalStars = 0;
  int oldGoodStars = 0;
  int oldBinCount = 0;
  int fileMagnitudeIndex;
  int magnitudeIndex;
  int oldGscBinIndex;
  int newGscBinIndex;
  int cacheIndex;
  PGSC_CACHE_ENTRY curCacheEntry;
  int curGscBinIndex = 0;

  char indexname[MAX_BUFFER];
  File indexHandle;
  char* catalogdir;
  char* slashPtr;
  char catalogname[MAX_BUFFER];
  File catalogHandle;
  GSC_CACHE_ENTRY pGscCache[GSC_CACHE_ENTRIES];

  PGSCIMAGE pCurGscImage;
  int readItems;
  int gscIndex;

  int usecIndex;
  SetQueryCount(0);
  double borderDist;
  double maxBorderDist = 0;
  int keepaliveFlag = 1;
  int workToDo;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon = &fileCommon;
  int filePrintModulus = 0;
  int curPrintMagnitudes = 0;
  int lastFileGscBinIndex = -MAG_FILE_MODULUS;
  int readOnly = 1;
  int iso_index;
  double iso_ratio;
  long long* iso0_badhistogram = NULL;
  long long iso0_badcount = 0;
  long long iso0_badcumulativeCount = 0;
  int iso0_badindex;
  long long* iso0_goodhistogram = NULL;
  long long iso0_goodcount = 0;
  long long iso0_goodcumulativeCount = 0;
  int iso0_goodindex;



  catalogString[0] = 0;

  memset(pGscCache,0,sizeof(pGscCache));


  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    printf("ERROR: DASCH_INGEST is not defined\n");
    return(-1);
  }

  outfile[0] = 0;
  qualifier[0] = 0;

 
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
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(outfile) >= MAX_BUFFER-3) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'i':  /* Ignore the versionId */
        case 'I':
          ignoreVersionId = 1;
          break;

        case 'r':  /* Repair the database */
        case 'R':
          repairMagnitudes = 1;
          readOnly = 0;
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

        case 's': /* suppress gsc catalog error messages */
        case 'S':
          suppressErrorMessages = 1;
          break;

        case 'p': /* create a limitingISO histogram */
        case 'P':
          limitingISOHistogram = 1;
          suppressErrorMessages = 1;
          readOnly = 1;
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


  outHandle = fopen(outfile,"a+t");
  if (outHandle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to open the output file %s\n",outfile);
  } else {
    if (verbose) {
      printf("Output file %s\n",outfile);
    }
  }


  catalogdir = getenv("DASCH_CATALOG");
  if (catalogdir == NULL) {
    printf("DASCH_CATALOG is not defined\n");
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



  catalogHandle = Open(catalogname,"r");
  if (catalogHandle == NULL) {
    printf("ERROR Could not open catalog file %s\n",catalogname);
    errorFlag = 1;
  } 
  indexHandle = Open(indexname,"r");
  if (indexHandle == NULL) {
    printf("ERROR Could not open catalog index file %s\n",indexname);
    errorFlag = 1;
  } 


  
  if (errorFlag) {
    printf("Usage: check_photometry options\n");
    printf("  options: -v verbose\n");    
    printf("           -r repair the database\n");
    printf("           -i ignore version id\n");
    printf("           -d debugging mode\n");
    printf("           -o <output file>\n");
    printf("           -q input catalog qualifer\n");
    printf("           -s Suppress GSC entry errors\n");
    printf("           -p Create an ISO<n> histogram for FILTER_AFLAG_LIMITING_MAG set\n");

    return(-1);
  }

  printf("check_photometry of %s %s, outfile %s magnitudeLimit %d PHOTSTARIMAGE: %d\n",
         __DATE__,__TIME__,outfile,MAGNITUDE_LIMIT,sizeof(PHOTSTARIMAGE));
  fprintf(outHandle,"check_photometry of %s %s, outfile %s magnitudeLimit %d PHOTSTARIMAGE %d\n",
          __DATE__,__TIME__,outfile,MAGNITUDE_LIMIT,sizeof(PHOTSTARIMAGE));


  if (limitingISOHistogram) {
    iso0_badhistogram = (long long *) calloc(ISO_MAX_VALUE,sizeof(long long));
    if (iso0_badhistogram == NULL) {
      printf("ERROR: failed to allocate iso0_badhistogram\n");
      exit(-1);
    }
    iso0_goodhistogram = (long long *) calloc(ISO_MAX_VALUE,sizeof(long long));
    if (iso0_goodhistogram == NULL) {
      printf("ERROR: failed to allocate iso0_goodhistogram\n");
      exit(-1);
    }



  }
#if 0
  readItems = 1;
#define MAX_ALLOCATED_BIN 70
  {
    int maxAlloc = 0;
    STARINDEX curStarIndex;
    PSTARINDEX pCurStarIndex = &curStarIndex;
    GSCIMAGE pGscImageTable[MAX_ALLOCATED_BIN];
    int tmpGscBin;
    int tmpDecBin;
    int tmpRaBin;
    while (readItems > 0) {
      readItems = Read(indexHandle,pCurStarIndex,sizeof(STARINDEX),1);
      if (readItems != 1) {
        break;
      }
      if(pCurStarIndex->numStars > maxAlloc) {
        maxAlloc = pCurStarIndex->numStars;
      }
      
      if ((pCurStarIndex->numStars) && (pCurStarIndex->numStars <= MAX_ALLOCATED_BIN)) {

        readItems = Read(catalogHandle,pGscImageTable,sizeof(GSCIMAGE),pCurStarIndex->numStars);

        if (readItems != pCurStarIndex->numStars) {
          break;
        }
        
        for (gscIndex = 0; gscIndex < pCurStarIndex->numStars; gscIndex++) {
          pCurGscImage = &pGscImageTable[gscIndex];
          if (strcmp(pCurGscImage->REF,"N2312232541") == 0) {

            tmpGscBin = GetGSCBin(pGscBin,pCurGscImage->ra,pCurGscImage->dec,&tmpDecBin,&tmpRaBin,(char *)"SearchStar2");

            printf("%s found in bin %d calculated bin is %d\n",pCurGscImage->REF,pCurStarIndex->binNumber,tmpGscBin);      
          }
        }
      }
    }
    printf("MaxAlloc is %d\n",maxAlloc);
  }
#endif



  /* Connect to the database */

  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    printf("ERROR: DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    printf("ERROR: DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    printf("ERROR: DASCH_PASSWORD is not defined\n");
    return(-1);
  }
                   
  mysql_init(pConnection);


  if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
    }
    return(-1);
  }


  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    printf("DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }

  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    printf("DASCH_PHOT_USERNAME is not defined\n");
    return(-1);
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    printf("DASCH_PHOT_PASSWORD is not defined\n");
    return(-1);
  }
  mysql_init(pPhotConnection);


  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
    }
    return(-1);
  }

  InitSeriesTable(NULL,pPhotConnection);

  if (GetPhotometryGlobal(pPhotConnection,pPhotGlobal) != 1) {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }
  if (pPhotGlobal->magnitudeFile == PHOT_MAGNITUDEFILE_YES) {
    magnitudeFileFlag = 1;
  } else {
    printf("ERROR: magnitudeFileFlag = 0 is no longer supported\n");
    exit(-1);
  }


  InitFileCommon(stdout,pFileCommon,pConnection,pPhotConnection,catalogNumber,1);
  InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);


  TimeStampPhotometryGlobal(pPhotConnection);



  time(&startTime);
  oldTime = startTime;
  ptr = localtime(&startTime);
  strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
 
  printf("check_photometry starting at %s\n",timestr);
  fprintf(outHandle,"check_photometry starting at %s\n",timestr);
  fflush(outHandle);

  strcpy(stars_name,ingestDirectory);
  strcat(stars_name,"/");
  strcat(stars_name,"check_photometry_stars.txt");

  strcpy(summary_name,ingestDirectory);
  strcat(summary_name,"/");
  strcat(summary_name,"check_photometry.txt");

  strcpy(debug_name,ingestDirectory);
  strcat(debug_name,"/");
  strcat(debug_name,"check_photometry.debug");
 
  summary_handle = fopen(summary_name,"wt");
  if (summary_handle == NULL) {
    printf("ERROR: Failed to open summary file  %s\n",summary_name);
    exit(-1);
  }


  sprintf(queryString,"update photglobal set keepalive = 'yes';");
  res = ExecuteQuery(pPhotConnection,queryString);
  if (res) {
    exit(-1);
  }
  /* Enter the main loop */
  minGscBinIndex = -1;
  maxGscBinIndex = -1;
  pMagnitudeTable  = NULL;
  curMagnitudes = -1;
  magnitudeIndex = 0;
  magnitudeAlloc = 0;
  curGscBinIndex = 0;

  /* Enter the main loop.  We could process each magnitude file individually before going on to the
     next, but then the performance would be extremely slow because we lose the locality of 
     reference into the GSC catalog file.  This more complicated code ensures that we read
     the GSC catalog file in sequence
  */

  while (1) {
    /* First see if we need to read in any data */
    if ((curMagnitudes < 0) || 
        (curGscBinIndex > maxGscBinIndex &&
         (curMagnitudes != 0))) {

      if (((++filePrintModulus % FILE_PRINT_MODULUS) == 0) || (magnitudeFileFlag == 0)) {
        filePrintModulus = 0;
	
        gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
        if ((gotAnswer != 1)  || (pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES)) {
      
          time(&curTime);

          ptr = localtime(&curTime);
          strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
 
          printf("check_photometry aborting at %s with gotAnswer %d keepalive %d\n",
                 timestr,gotAnswer,pPhotGlobal->keepalive);
          fprintf(outHandle,"check_photometry aborting at %s with gotAnswer %d keepalive %d\n",
                  timestr,gotAnswer,pPhotGlobal->keepalive);
          keepaliveFlag = 0;

          break;
        }
      }
      minGscBinIndex = curGscBinIndex;
      maxGscBinIndex = curGscBinIndex+MAG_FILE_MODULUS-1;
      curMagnitudes = GetFileSummaryMagnitudes(pGscBin,pFileCommon,curGscBinIndex,curGscBinIndex+MAG_FILE_MODULUS-1,0,0,0,0,NULL,catalogString,0,NULL,readOnly,0);
      /* For now, just copy the results over */
      if (magnitudeAlloc < curMagnitudes) {
        free(pMagnitudeTable);
        magnitudeAlloc = curMagnitudes;
        pMagnitudeTable = (PPHOTSTARIMAGE)calloc(magnitudeAlloc,sizeof(PHOTSTARIMAGE));
        if (pMagnitudeTable == NULL) {
          printf("ERROR: failed to allocate %d PHOTSTARIMAGEs in GetSummaryMagnitudes\n",magnitudeAlloc);
          exit(-1);
        }
      }
      pCurMagnitudeTable = pMagnitudeTable;
      for (fileMagnitudeIndex = 0; fileMagnitudeIndex < curMagnitudes;fileMagnitudeIndex++) {
        pCurStarImage = &pCurMagnitudeTable[fileMagnitudeIndex];
        pFileStarImage = &pFileCommon->magnitudeBuffer[fileMagnitudeIndex];
        if (limitingISOHistogram) {
          /* Because MAX_LIMITING_MAG was changed on Feb 15, 2013, this code enforces the new setting */
          if ((pFileStarImage->limiting_mag_local < 90.) &&
              (pFileStarImage->magcal_magdep < 90.) &&
              ((pFileStarImage->AFLAGS & FILTER_AMASK_REJECT) == 0)) {
            
            if ((pFileStarImage->limiting_mag_local - pFileStarImage->magcal_magdep) < MAX_LIMITING_MAG) {
              pFileStarImage->AFLAGS |= (1<<FILTER_AFLAG_LIMITING_MAG);
            } else {
              pFileStarImage->AFLAGS &= (~(1<<FILTER_AFLAG_LIMITING_MAG));
            }
         
            if ((pFileStarImage->AFLAGS & (1 << FILTER_AFLAG_LIMITING_MAG)) != 0) {
              iso0_badindex = pFileStarImage->ISO0;
              if (iso0_badindex < 0) {
                iso0_badindex = 0;
              } else if (iso0_badindex >= ISO_MAX_VALUE) {
                iso0_badindex = ISO_MAX_VALUE-1;
              }
              iso0_badhistogram[iso0_badindex]++;
              iso0_badcount++;


            } else {
              iso0_goodindex = pFileStarImage->ISO0;
              if (iso0_goodindex < 0) {
                iso0_goodindex = 0;
              } else if (iso0_goodindex >= ISO_MAX_VALUE) {
                iso0_goodindex = ISO_MAX_VALUE-1;
              }
              iso0_goodhistogram[iso0_goodindex]++;
              iso0_goodcount++;


            }
          }
          continue;
        }


        if (pFileStarImage->timeAccuracy == 0) {
          printf("ERROR: timeAccuracy is zero for %s%05d s%d exposure %d\n",GetSeriesString(pFileStarImage->seriesId,0),pFileStarImage->plateNumber,pFileStarImage->solutionNumber,pFileStarImage->exposureNumber);
        }
        memset(pCurStarImage,0,sizeof(PHOTSTARIMAGE));
        pCurStarImage->pFileStarImage = pFileStarImage;
        strcpy(pCurStarImage->series,GetSeriesString(pFileStarImage->seriesId,1));
        GetREF(pFileStarImage->REFNumber,pCurStarImage->REF,0,1);

      }

        
   	  
      
    
#if 0
      if (curMagnitudes == 0) {
        printf("curMagnitudes is %d \n",curMagnitudes);
      }
#endif
      curPrintMagnitudes += curMagnitudes;
      if ((minGscBinIndex-lastFileGscBinIndex) != MAG_FILE_MODULUS) {
        printf("ERROR: files not searched increment %d: %d-%d\n",
               MAG_FILE_MODULUS,
               lastFileGscBinIndex,
               minGscBinIndex);
      }
      lastFileGscBinIndex = minGscBinIndex;
      if (verbose) {
	  
        if (((filePrintModulus % FILE_PRINT_MODULUS) == 0) || (magnitudeFileFlag == 0)) {
          time(&curTime);
          curTime -= startTime;
          printf("Read %10d magnitudes for curGscBinIndex %10d at %5d seconds\n",curPrintMagnitudes,minGscBinIndex,curTime);
          curPrintMagnitudes = 0;
        }
      }
      if (curMagnitudes > 0) {
        magnitudeIndex = 0;
        totalMagnitudes += curMagnitudes;
        pCurStarImage = pCurMagnitudeTable;
        if (limitingISOHistogram == 0) {
          pFileStarImage = pCurStarImage->pFileStarImage;
          minGscBinIndex = pFileStarImage->gsc_bin_index;

          pCurStarImage = &pCurMagnitudeTable[curMagnitudes-1];
          pFileStarImage = pCurStarImage->pFileStarImage;
        }
      }
    }
    
    workToDo = 0;
    /* Now select our next object to read */
    if (curMagnitudes > 0) {
      if (workToDo == 0) {
        curGscBinIndex = minGscBinIndex;
      } else {
        if (minGscBinIndex < curGscBinIndex) {
          curGscBinIndex = minGscBinIndex;
        }
      }
      workToDo = 1;
    }

    if ((curMagnitudes == 0) || ( limitingISOHistogram != 0)) {
      curGscBinIndex += MAG_FILE_MODULUS;
      curMagnitudes = -1;
      if (curGscBinIndex > pGscBin->total_gsc_bins) {
        printf("curGscBinIndex %d exceeds maximum %d\n",curGscBinIndex,pGscBin->total_gsc_bins);
        break;
      }
	
      continue;
    }
    
    pCurMagnitudeTable = pMagnitudeTable;
    pCurStarImage = &pCurMagnitudeTable[magnitudeIndex];
    pFileStarImage = pCurStarImage->pFileStarImage;
 
    totalProcessedMagnitudes++;
    if (pFileStarImage->REFNumber != 0) {
      if ((pFileStarImage->AFLAGS & (1 << FILTER_AFLAG_BLEND_NOMATCH)) != 0) {
        totalBlendNomatchMagnitudes++;
      }
      /* Search for our star */
      pCurGscImage = SearchStar(pGscBin,pGscCache,indexHandle,catalogHandle,pFileStarImage->gsc_bin_index,&newGscBinIndex,pCurStarImage->REF);

      if (pCurGscImage == NULL) {
        curStars = GetStarEntry(pPhotConnection,&pCurStarEntry,1,UPDATEFLAG_IGNORE,pCurStarImage->REF,0,catalogString,1,0);
        if (suppressErrorMessages == 0) {
          if (curStars != 1) {
            printf("WARNING: no gsc entry found in bin %d for star %s plate %s%05d AFLAGS %d 0x%08x\n",
                   pFileStarImage->gsc_bin_index,
                   pCurStarImage->REF,
                   pCurStarImage->series,
                   pFileStarImage->plateNumber,
                   pFileStarImage->AFLAGS,
                   pFileStarImage->AFLAGS);
          } else {
            printf("WARNING: no gsc entry in bin %d should be %d for star %s ra %f dec %f plate %s%05d  AFLAGS %d 0x%08x\n",
                   pFileStarImage->gsc_bin_index,
                   pCurStarEntry->gsc_bin_index,
                   pCurStarImage->REF,
                   pCurStarEntry->ra,
                   pCurStarEntry->dec,
                   pCurStarImage->series,
                   pFileStarImage->plateNumber,
                   pFileStarImage->AFLAGS,
                   pFileStarImage->AFLAGS);
          }
        }
        noGscEntryCount++;
      } else {
        if (newGscBinIndex != pFileStarImage->gsc_bin_index) {
          recoveredStars++;
          

#if 0
          borderDist = GetGSCBorder(pGscBin,newGscBinIndex,pCurGscImage->ra,pCurGscImage->dec,pCurStarImage->REF);
          if (borderDist > maxBorderDist) {
            maxBorderDist = borderDist;
          }

          printf("WARNING: The gsc entry in bin %d should be %d for star %s ra %f dec %f plate %s%05d distance %f degrees  AFLAGS %d 0x%08x \n",
                 pFileStarImage->gsc_bin_index,
                 newGscBinIndex,
                 pCurStarImage->REF,
                 pCurGscImage->ra,
                 pCurGscImage->dec,
                 pCurStarImage->series,
                 pFileStarImage->plateNumber,
                 borderDist,
                 pFileStarImage->AFLAGS,
                 pFileStarImage->AFLAGS);
#endif
          if (repairMagnitudes) {
            if (magnitudeFileFlag) {
              printf("ERROR: magnitude repair is not implemented for the magnitude file format\n");
              exit(-1);
            }
          }

        }
        if (pCurGscImage->flag == 0 ) {
          pCurGscImage->flag = 1;
          totalStars++;

        }

      }
    } else {
      totalNoneMagnitudes++;
    }

#if 0
    /* On to the next star */
    if (curGscBinIndex == 129971741) {
      printf("At gsc_bin_index %d\n",curGscBinIndex);
    }
#endif

    magnitudeIndex++;
    if (magnitudeIndex < curMagnitudes) {
      pCurStarImage = &pCurMagnitudeTable[magnitudeIndex];
      pFileStarImage = pCurStarImage->pFileStarImage;
      minGscBinIndex = pFileStarImage->gsc_bin_index;
    } else {
      curMagnitudes = -1;
      curGscBinIndex = maxGscBinIndex+1;
      

    }

  }
  totalMagnitudeAlloc += magnitudeAlloc;
  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d magnitudeAlloc %d gscBins, %d total magnitudes\n",curTime,GetQueryCount(),totalStars,recoveredStars,badModulusBin,totalMagnitudeAlloc,binCount,totalMagnitudes);
  if (limitingISOHistogram) {
    printf("ISO0\tcatalogNumber\tcountBad\tcumulativeCountBad\tcountGood\tcumulativeCountGood\tratio\n");
    printf("----\t-------------\t--------\t------------------\t---------\t-------------------\t-----\n");
    for (iso_index = 0; iso_index < ISO_MAX_VALUE; iso_index++) {
      if ((iso0_badhistogram[iso_index] != 0) ||
          (iso0_goodhistogram[iso_index] != 0)) {
        iso0_badcumulativeCount += iso0_badhistogram[iso_index];
        iso0_goodcumulativeCount += iso0_goodhistogram[iso_index];
        iso_ratio = (1.0*iso0_goodcumulativeCount)/(1.0*(iso0_badcumulativeCount+iso0_goodcumulativeCount));
        printf("%d\t%d\t%lld\t%f\t%lld\t%f\t%f\n",iso_index,catalogNumber,iso0_badhistogram[iso_index],(1.0*iso0_badcumulativeCount)/(1.0*iso0_badcount),iso0_goodhistogram[iso_index],(1.0*iso0_goodcumulativeCount)/(1.0*iso0_goodcount),iso_ratio);
      }

    }
    printf("ISO0\tdone\n");

  }


  fprintf(outHandle,"Execution Time: %d seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d magnitudeAlloc %d gscBins, %d total magnitudes\n",curTime,GetQueryCount(),totalStars,recoveredStars,badModulusBin,totalMagnitudeAlloc,binCount,totalMagnitudes);
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


  printf("binSearchCount %d, cacheHits %d, cacheMisses %d, maxBorderDist %f\n",binSearchCount,cacheHits,cacheMisses,maxBorderDist);
  fprintf(outHandle,"binSearchCount %d, cacheHits %d, cacheMisses %d, maxBorderDist %f\n",binSearchCount,cacheHits,cacheMisses,maxBorderDist);

  printf("Total Magnitudes %d, Total NONE magnitudes %d Total NOMATCH magnitudes %d Total GSC nonblend magnitudes %d \n",
         totalProcessedMagnitudes,totalNoneMagnitudes,totalBlendNomatchMagnitudes,totalProcessedMagnitudes - totalNoneMagnitudes - totalBlendNomatchMagnitudes);
  fprintf(outHandle,"Total Magnitudes %d, Total NONE magnitudes %d Total NOMATCH magnitudes %d Total GSC nonblend magnitudes %d \n",
          totalProcessedMagnitudes,totalNoneMagnitudes,totalBlendNomatchMagnitudes,totalProcessedMagnitudes - totalNoneMagnitudes - totalBlendNomatchMagnitudes);

  mysql_close(pConnection);
  mysql_close(pPhotConnection);
  if (pMagnitudeTable != NULL) {
    free(pMagnitudeTable);
  }
  for (cacheIndex = 0; cacheIndex < GSC_CACHE_ENTRIES; cacheIndex++) {
    curCacheEntry = &pGscCache[cacheIndex];
    printf("Freeing cache entry %d with %d entries\n",cacheIndex,curCacheEntry->gscImageAlloc);
    if (curCacheEntry->pGscImageTable != NULL) {
      free(curCacheEntry->pGscImageTable);
    }

  }
  if (debug_handle != NULL) {
    fclose(debug_handle);
    debug_handle = NULL;
  }
  if (iso0_badhistogram != NULL) {
    free(iso0_badhistogram);
  }
  FreeFileCommon(pFileCommon,1);
  Close(catalogHandle);
  Close(indexHandle);
  fclose(outHandle);
  return(0);
}

