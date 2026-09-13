// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* combine_photometry.c
 *
 * Combine the APASS and GSC2.3.2 catalogs to produce a new experimental catalog
 * 
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   combine_photometry.c pipelineutils.a -ltable -lutil  -lwcs -o combine_photometry -L/usr/lib64/mysql -lmysqlclient
 *
 * Jan 19, 2018 Edward J. Los - Initial version
 * Mar 19, 2018 Edward J. Los - Remove support for combined binary photometry files.  Allow only combining the limiting magnitude files.
 * Jun  1, 2018 Edward J. Los - Correct sanity check for two different mosaic numbers
 *
     combine_photometry -v -i -p -l -q apass gsc2.3.2 experimental -o  /dasch/Pipeline/ingest/combine_photometry.log

 *
 */


#include <math.h>
#include <errno.h>
#include "table.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <fcntl.h>
#include "pipelineutils.h"
#include "mysql.h"
#include "photometryutils.h"
#include "searchgsc.h"
#include "daschunistd.h"

#define MAX_BUFFER 256
#define VECTOR_ALLOC_INCREMENT 5003
#define RESET_DELTA_ENTRY -1
#define RESET_ALL_ENTRY -2
#define FILE_PRINT_MODULUS 5000
#define MAX_FILENAME 256
#define LIMITING_BUFFER_INCREMENT 10000


extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

extern GSCBIN gscBin01;
PGSCBIN pGscBin01 = &gscBin01;

int binSearchCount = 0;
int cacheHits = 0;
int cacheMisses = 0;
int invokeCmd(char *cmdStr,int *timeSpent) {
  int result;
  time_t beginTime;
  time_t curTime;
  time(&beginTime);
  *timeSpent = 0;
  result = system(cmdStr);
  time(&curTime);
  beginTime = curTime - beginTime;
  *timeSpent = (int)beginTime;
  if (WIFSIGNALED(result) &&
      (WTERMSIG(result) == SIGINT || WTERMSIG(result) == SIGQUIT)) {
    printf("Terminated with signal %d at %d seconds\n",WTERMSIG(result),(int)curTime);
    exit(-1);
  }
#if 0
  if (WIFEXITED(result)) {
    printf("Exited with status %d at %d sec delta %d sec\n",WEXITSTATUS(result),(int)curTime,(int)beginTime);
  } else {
    printf("Unknown result %d at %d sec delta %d \n",result,(int)curTime,(int)beginTime);
  }
  printf("result: %d\n",result);
#endif

  return(result);
}

/* Now sort by gsc_bin_number, seriesId,plateNumber,NUMBER and catalogNumber */
int OldCatalogImageCompare(const void *first, const void *second) 
{
  PFILESTARIMAGE pFirst = (PFILESTARIMAGE)first;
  PFILESTARIMAGE pSecond = (PFILESTARIMAGE)second;
  if ((pFirst->catalogNumber == CATALOG_GAIA) || 
      (pFirst->catalogNumber == CATALOG_GAIA)) {
    printf("ERROR in line %d of combine_photometry - CATALOG_GAIA is not supported\n",__LINE__);
    exit(-1);
  }
  if ((pFirst->catalogNumber == CATALOG_ATLAS) || 
      (pFirst->catalogNumber == CATALOG_ATLAS)) {
    printf("ERROR in line %d of combine_photometry - CATALOG_ATLAS is not supported\n",__LINE__);
    exit(-1);
  }

  if (pFirst->gsc_bin_index > pSecond->gsc_bin_index) {
    return(1);
  } else if (pFirst->gsc_bin_index < pSecond->gsc_bin_index) {
    return(-1);
  } else {

    if (pFirst->seriesId > pSecond->seriesId) {
      return(1);
    } else if (pFirst->seriesId < pSecond->seriesId) {
      return(-1);
    } else {

      if (pFirst->plateNumber > pSecond->plateNumber) {
        return(1);
      } else if (pFirst->plateNumber < pSecond->plateNumber) {
        return(-1);
      } else {

        if (pFirst->NUMBER > pSecond->NUMBER) {
          return(1);
        } else if (pFirst->NUMBER < pSecond->NUMBER) {
          return(-1);
        } else {
#if 0
          if (((pFirst->REFNumber == 120132012267L) && (pSecond->REFNumber == 405144112674112L)) ||
              ((pSecond->REFNumber == 120132012267L) && (pFirst->REFNumber == 405144112674112L))) {
            printf("DUMP line %4d gsc_bin_index %d seriesId %2d %3s plateNumber %d NUMBER %d catalogNumber %d REFNumber %lld\n",
                   __LINE__,
                   pFirst->gsc_bin_index,
                   pFirst->seriesId,
                   GetSeriesString(pFirst->seriesId,0),
                   pFirst->plateNumber,
                   pFirst->NUMBER,
                   pFirst->catalogNumber,
                   pFirst->REFNumber);
          }
#endif
          /* Favor the APASS catalog */
          if ((pFirst->catalogNumber == CATALOG_APASS) && pSecond->catalogNumber != CATALOG_APASS) {
            return(1);
          } else if ((pFirst->catalogNumber != CATALOG_APASS) && pSecond->catalogNumber == CATALOG_APASS) {
            return(-1);
          } else {
            if (pFirst->catalogNumber > pSecond->catalogNumber) {
              return(1);
            } else if (pFirst->catalogNumber < pSecond->catalogNumber) {
              return(-1);
            } else {
              return(0);
            }
          }

        }

      }

    }

  }
    

}
/* Sort by geoJulianDate,plate number and deepest limiting magnitude */
int LimitingMagCompare(const void *first, const void *second) 
{
  PPLATELIMITINGREC pFirst = (PPLATELIMITINGREC)first;
  PPLATELIMITINGREC pSecond = (PPLATELIMITINGREC)second;
  if (pFirst->geoJulianDate > pSecond->geoJulianDate) {
    return(1);
  } else if (pFirst->geoJulianDate < pSecond->geoJulianDate) {
    return(-1);
  } else {
    if (pFirst->seriesId > pSecond->seriesId) {
      return(1);
    } else if (pFirst->seriesId < pSecond->seriesId) {
      return(-1);
    } else {
      if (pFirst->plateNumber > pSecond->plateNumber) {
        return(1);
      } else if (pFirst->plateNumber < pSecond->plateNumber) {
        return(-1);
      } else {
        if (pFirst->mosaicNumber > pSecond->mosaicNumber) {
          return(1);
        } else if (pFirst->mosaicNumber < pSecond->mosaicNumber) {
          return(-1);
        } else {
          /* Return the deepest limiting magnitude found */
          if (pFirst->limiting_mag_local > pSecond->limiting_mag_local) {
            return(1);
          } else if (pFirst->limiting_mag_local < pSecond->limiting_mag_local) {
            return(-1);
          } else {
            return(0);
          }
        }
      }
    }
  }    
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
  int verbose = 0;
  int combineBinaryFiles = 0;
  int combineLimitingFiles = 0;
  int debugMode = 0;
  int ignoreVersionId = 0;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char qualifier[3][MAX_BUFFER];
  int nvals;
  FILE *outHandle;
  char timestr[100];
  struct tm *ptr;
  int gotAnswer;
  int curMagnitudes0;
  int curMagnitudes1;
  int totalMagnitudes = 0;
  int totalProcessedMagnitudes = 0;
  int totalNoneMagnitudes = 0;
  int totalBlendNomatchMagnitudes = 0;
  int minGscBinIndex = 0;
  int maxGscBinIndex = 0;
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
  PFILESTARIMAGE pFileStarImage2 = NULL;


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
  int curGscBinIndex = 0;

  char indexname[3][MAX_BUFFER];
  File indexHandle[3];
  char* catalogdir;
  char* slashPtr;
  char catalogname[3][MAX_BUFFER];
  File catalogHandle[3];
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
  int catalogNumber[3];
  char catalogString[3][MAX_BUFFER];
  FILECOMMON fileCommon0;
  FILECOMMON fileCommon1;
  FILECOMMON fileCommon2;
  PFILECOMMON pFileCommon0 = &fileCommon0;
  PFILECOMMON pFileCommon1 = &fileCommon1;
  PFILECOMMON pFileCommon2 = &fileCommon2;
  GALAXYCOMMON galaxyCommon0;
  GALAXYCOMMON galaxyCommon1;
  GALAXYCOMMON galaxyCommon2;
  PGALAXYCOMMON pGalaxyCommon0 = &galaxyCommon0;
  PGALAXYCOMMON pGalaxyCommon1 = &galaxyCommon1;
  PGALAXYCOMMON pGalaxyCommon2 = &galaxyCommon2;
  int haveGalaxyTable0;
  int haveGalaxyTable1;
  double limitingra;
  double limitingdec;
  PPLATELIMITINGREC pPlateLimitingRec0;
  PPLATELIMITINGREC pPlateLimitingRec1;
  PPLATELIMITINGREC pPlateLimitingRec2;
  PPLATELIMITINGREC pPlateLimitingRec;
  PPLATELIMITINGREC pPlateLimitingBuffer = NULL;
  int plateIndex;
  int plateIndexRead;
  int plateIndexWrite;
  int plateIndexSaved;
  PPLATELIMITINGREC pPlateLimitingRecRead;
  PPLATELIMITINGREC pPlateLimitingRecSaved;
  PPLATELIMITINGREC pPlateLimitingRecWrite;
  double savedJulianDate;
  int savedSeriesId;
  int savedPlateNumber;
  int savedMosaicNumber;
  int savedlimiting_mag_local;

  int totalPlateCount;
  int allocPlateCount = 0;

  int filePrintModulus = 0;
  int curPrintMagnitudes = 0;
  int lastFileGscBinIndex = -MAG_FILE_MODULUS;
  int catalogIndex;

  PSUBSTITUTION substitutionTable = NULL;
  PSUBSTITUTION pSubstitution;
  int substitutionCount = 0;
  int substitutionIndex;
  PFILESTARIMAGE cacheStarImages = NULL;
  int cacheStarCount = 0;
  int cacheStarAlloc = 0;
  int cacheStarIndex = 0;
  int REFNumberReplacementCount = 0;
  char *photfilebase;
  char dasch_phot_magnitudes[40]; 
  sem_t *fileSemaphore = NULL;
  int directoryCreates = 0;

  int dir1;
  int dir2;
  size_t writeBytes;
  char filename[MAX_FILENAME];
  char chgrpname[MAX_FILENAME];
  char cmdstr[MAX_FILENAME];
  int timeSpent;
  int result;
  int fileNumber;
  int statResult;
  struct stat filestats;
  int debugPrint = 1;
  int debugCount = 0;
  STARINDEX starIndexRecord;
  PSTARINDEX pStarIndex = &starIndexRecord;
  off_t platerec_bytes = 0;
  int totalGalaxyWrite;
  off_t byteCount;
  int objectIndex;
  int writeItems;
  PGALAXYREC pGalaxyRec;


  for (catalogIndex = 0; catalogIndex < 3; catalogIndex++) {

    catalogNumber[catalogIndex] = 0;
    catalogString[catalogIndex][0] = 0;
    catalogHandle[catalogIndex] = NULL;
    indexHandle[catalogIndex] = NULL;
    qualifier[catalogIndex][0] = 0;
  }

  memset(pGscCache,0,sizeof(pGscCache));


  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    printf("ERROR: DASCH_INGEST is not defined\n");
    return(-1);
  }

  outfile[0] = 0;
 
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

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'd': /* debug mode */
        case 'D':
          verbose = 1;
          debugMode = 1;
          break;

        case 'p': /* do the binary magnitude files */
        case 'P':
          combineBinaryFiles = 1;
          printf("ERROR: combining photometry binary files is no longer supported\n");
          errorFlag = 1;
          break;

        case 'l': /* do the limiting magnitude */
        case 'L':
          combineLimitingFiles = 1;
          break;

 

        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          if (argc < 4) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else { 
            argc--;
            catalogNumber[0] = GetCatalogNumber(*++argv);
            if (catalogNumber[0] < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              if (catalogNumber[0] == 0) {
                catalogString[0][0] = 0;
              } else {
                sprintf(catalogString[0],"%d",catalogNumber[0]);
                strcpy(qualifier[0],*argv);
              }
            }
            argc--;
            catalogNumber[1] = GetCatalogNumber(*++argv);
            if (catalogNumber[1] < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              if (catalogNumber[1] == 0) {
                catalogString[1][0] = 0;
              } else {
                sprintf(catalogString[1],"%d",catalogNumber[1]);
                strcpy(qualifier[1],*argv);
              }
            }
            argc--;
            catalogNumber[2] = GetCatalogNumber(*++argv);
            if (catalogNumber[2] != CATALOG_EXPERIMENTAL) {
              printf("ERROR: output catalog %s is not %s\n",
                     catalogText[catalogNumber[2]],
                     catalogText[CATALOG_EXPERIMENTAL]);
              exit(-1);
            }
            if (catalogNumber[2] < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              if (catalogNumber[2] == 0) {
                catalogString[2][0] = 0;
              } else {
                sprintf(catalogString[2],"%d",catalogNumber[2]);
                strcpy(qualifier[2],*argv);
              }
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

  if ((combineBinaryFiles == 0) && (combineLimitingFiles == 0)) {
    printf("Error: neither -p nor -l options are specified\n");
    errorFlag = 1;
  }

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
  
  
  for (catalogIndex = 0; catalogIndex < 3; catalogIndex++) {

    strcpy(catalogname[catalogIndex],catalogdir);
    if (strstr(qualifier[catalogIndex],"kepler")) {
      slashPtr = strrchr(catalogname[catalogIndex],'/');
      if (slashPtr != NULL) {
        slashPtr++;
      } else {
        slashPtr = catalogname[catalogIndex];
      }
      *slashPtr = 0;
      strcat(catalogname[catalogIndex],"keplerlimiting.dat");
    }

    if (strstr(qualifier[catalogIndex],"apass")) {
      slashPtr = strrchr(catalogname[catalogIndex],'/');
      if (slashPtr != NULL) {
        slashPtr++;
      } else {
        slashPtr = catalogname[catalogIndex];
      }
      *slashPtr = 0;
      strcat(catalogname[catalogIndex],"apasslimiting.dat");
    }

    if (strstr(qualifier[catalogIndex],"experimental")) {
      slashPtr = strrchr(catalogname[catalogIndex],'/');
      if (slashPtr != NULL) {
        slashPtr++;
      } else {
        slashPtr = catalogname[catalogIndex];
      }
      *slashPtr = 0;
      strcat(catalogname[catalogIndex],"experimentallimiting.dat");
    }


    strcpy(indexname[catalogIndex],catalogname[catalogIndex]);
    dotPtr = strstr(indexname[catalogIndex],".");
    if (dotPtr != NULL) {
      *dotPtr = 0;
    }
    strcat(indexname[catalogIndex],".idx");

#if 0
    catalogHandle[catalogIndex] = Open(catalogname[catalogIndex],"r");
    if (catalogHandle[catalogIndex] == NULL) {
      printf("ERROR Could not open catalog file %s\n",catalogname[catalogIndex]);
      errorFlag = 1;
    } 
    indexHandle[catalogIndex] = Open(indexname[catalogIndex],"r");
    if (indexHandle[catalogIndex] == NULL) {
      printf("ERROR Could not open catalog index file %s\n",indexname[catalogIndex]);
      errorFlag = 1;
    } 
#endif

  }
  sprintf(dasch_phot_magnitudes,"DASCH_PHOT_MAGNITUDES%s",catalogString[2]);
  photfilebase = getenv(dasch_phot_magnitudes);
  
  if (errorFlag) {
    printf("Usage: combine_photometry options\n");
    printf("  options: -v verbose\n");    
    printf("           -i ignore version id\n");
    printf("           -o <output file>\n");
    printf("           -q <input1> <input2> <output> catalog qualifers\n");
    printf("           -p combine the binary photometry files\n");
    printf("           -l combine the limiting magnitude files\n");
    return(-1);
  }

  printf("combine_photometry of %s %s, combining %s and %s to produce %s outfile %s magnitudeLimit %d PHOTSTARIMAGE: %d\n",
         __DATE__,__TIME__,qualifier[0],qualifier[1],qualifier[2],outfile,sizeof(PHOTSTARIMAGE));
  fprintf(outHandle,"combine_photometry of %s %s, combining %s and %s to produce %s outfile %s magnitudeLimit %d PHOTSTARIMAGE: %d\n",
          __DATE__,__TIME__,qualifier[0],qualifier[1],qualifier[2],outfile,sizeof(PHOTSTARIMAGE));




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

  InitSeriesTable(pConnection,pPhotConnection);

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


  InitFileCommon(stdout,pFileCommon0,pConnection,pPhotConnection,catalogNumber[0],1);
  InitMaxPlateNumber(pConnection,pFileCommon0->maxPlateNumber);
  InitFileCommon(stdout,pFileCommon1,pConnection,pPhotConnection,catalogNumber[1],1);
  InitMaxPlateNumber(pConnection,pFileCommon1->maxPlateNumber);
  InitFileCommon(stdout,pFileCommon2,pConnection,pPhotConnection,catalogNumber[2],1);
  InitMaxPlateNumber(pConnection,pFileCommon2->maxPlateNumber);





  time(&startTime);
  oldTime = startTime;
  ptr = localtime(&startTime);
  strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
 
  printf("combine_photometry starting at %s\n",timestr);
  fprintf(outHandle,"combine_photometry starting at %s\n",timestr);
  fflush(outHandle);

  strcpy(stars_name,ingestDirectory);
  strcat(stars_name,"/");
  strcat(stars_name,"combine_photometry_stars.txt");

  strcpy(summary_name,ingestDirectory);
  strcat(summary_name,"/");
  strcat(summary_name,"combine_photometry.txt");

  strcpy(debug_name,ingestDirectory);
  strcat(debug_name,"/");
  strcat(debug_name,"combine_photometry.debug");
 
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
#if 0   /* /dasch/data/junk/filetest3/mag006/mag326/mag006326272.dat */
  /* RA,Dec 79.3151 -67.6802   APASS_J052142.5-674046 has 2184 points   APASS_J051258.3-674024 has 3308 APASS_J051441.1-674112 has 3761 points, but GSC has a separate id */
  /* APASS_J051441.1-674112 (405144112674112)  is confused with S0132012267   (120132012267) */
  minGscBinIndex = 6326272;
  maxGscBinIndex = 6326272+1024;
#endif
#if 1
  minGscBinIndex = 0;
  maxGscBinIndex = pGscBin->total_gsc_bins;
#endif
  pMagnitudeTable  = NULL;
  magnitudeIndex = 0;
  magnitudeAlloc = 0;
  curGscBinIndex = 0;

  if (combineBinaryFiles) {
    printf("Combining the binary photometry magnitude files\n");
  
    for (curGscBinIndex = minGscBinIndex; curGscBinIndex < maxGscBinIndex; curGscBinIndex += MAG_FILE_MODULUS) {
      if (((++filePrintModulus % FILE_PRINT_MODULUS) == 0) || (magnitudeFileFlag == 0)) {
        filePrintModulus = 0;
	
        gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
        if ((gotAnswer != 1)  || (pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES)) {
      
          time(&curTime);

          ptr = localtime(&curTime);
          strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
 
          printf("combine_photometry aborting at %s with gotAnswer %d keepalive %d\n",
                 timestr,gotAnswer,pPhotGlobal->keepalive);
          fprintf(outHandle,"combine_photometry aborting at %s with gotAnswer %d keepalive %d\n",
                  timestr,gotAnswer,pPhotGlobal->keepalive);
          keepaliveFlag = 0;

          break;
        }
      }
      curMagnitudes0 = GetFileSummaryMagnitudes(pGscBin,pFileCommon0,curGscBinIndex,curGscBinIndex+MAG_FILE_MODULUS-1,0,0,0,0,NULL,catalogString[0],0,NULL,1,1);
      curMagnitudes1 = GetFileSummaryMagnitudes(pGscBin,pFileCommon1,curGscBinIndex,curGscBinIndex+MAG_FILE_MODULUS-1,0,0,0,0,NULL,catalogString[1],0,NULL,1,1);

      if (cacheStarAlloc < (curMagnitudes0+curMagnitudes1)) {
        if (cacheStarImages != NULL) {
          free(cacheStarImages);
          cacheStarImages = NULL;
        }
        if (substitutionTable != NULL) {
          free(substitutionTable);
          substitutionTable = NULL;
        }

        cacheStarAlloc = (curMagnitudes0+curMagnitudes1)+1000;
        cacheStarImages = (PFILESTARIMAGE)calloc(cacheStarAlloc,sizeof(FILESTARIMAGE));
        if (cacheStarImages == NULL) {
          printf("ERROR: failed to allocate cacheStarImages for %d stars\n",cacheStarAlloc);
        }

        substitutionTable = (PSUBSTITUTION)calloc(cacheStarAlloc,sizeof(SUBSTITUTION));
        if (substitutionTable == NULL) {
          printf("ERROR: failed to allocate substitutionTable for %d stars\n",cacheStarAlloc);
        }
      }
      cacheStarCount = 0;
      substitutionCount = 0;
      for (fileMagnitudeIndex = 0; fileMagnitudeIndex < curMagnitudes0;fileMagnitudeIndex++) {
        pFileStarImage = &pFileCommon0->magnitudeBuffer[fileMagnitudeIndex];
        pFileStarImage2 = &cacheStarImages[cacheStarCount];
        memcpy(pFileStarImage2,pFileStarImage,sizeof(FILESTARIMAGE));
        pFileStarImage2->catalogNumber = catalogNumber[0];
#if 0
        pFileStarImage2->versionId = pPhotGlobal->currentVersion;
#endif
        cacheStarCount++;
        if (cacheStarCount >= cacheStarAlloc) {
          printf("ERROR: line %d cacheStarCount %d less than cacheStarAlloc %d\n",__LINE__,cacheStarCount,cacheStarAlloc);
          exit(-1);
        }
      }
      for (fileMagnitudeIndex = 0; fileMagnitudeIndex < curMagnitudes1;fileMagnitudeIndex++) {
        pFileStarImage = &pFileCommon1->magnitudeBuffer[fileMagnitudeIndex];
        pFileStarImage2 = &cacheStarImages[cacheStarCount];
        memcpy(pFileStarImage2,pFileStarImage,sizeof(FILESTARIMAGE));
        pFileStarImage2->catalogNumber = catalogNumber[1];
        pFileStarImage2->versionId = pPhotGlobal->currentVersion;
        cacheStarCount++;
        if (cacheStarCount >= cacheStarAlloc) {
          printf("ERROR: line %d cacheStarCount %d less than cacheStarAlloc %d\n",__LINE__,cacheStarCount,cacheStarAlloc);
          exit(-1);
        }
      }
      /* Now sort by gsc_bin_number, seriesId,plateNumber,NUMBER, and catalogNumber */
      qsort((void*)cacheStarImages,cacheStarCount,sizeof(FILESTARIMAGE),OldCatalogImageCompare);
      /* Now we want to replace any non-APASS REFNumbers with APASS refnumbers */
      for (cacheStarIndex = 1; cacheStarIndex < cacheStarCount; cacheStarIndex++) {
        pFileStarImage =  &cacheStarImages[cacheStarIndex-1];

        pFileStarImage2 =  &cacheStarImages[cacheStarIndex];
#if 0
        if (pFileStarImage2->REFNumber == 405144112674112L) {
          printf("At REFNumber %lld\n",pFileStarImage2->REFNumber);
          debugPrint = 1;
        }
#endif

        if ((pFileStarImage->gsc_bin_index == pFileStarImage2->gsc_bin_index) &&
            (pFileStarImage->seriesId == pFileStarImage2->seriesId) &&
            (pFileStarImage->plateNumber == pFileStarImage2->plateNumber) &&
            (pFileStarImage->NUMBER == pFileStarImage2->NUMBER) &&
            (pFileStarImage->catalogNumber != CATALOG_APASS) &&
            (pFileStarImage2->catalogNumber == CATALOG_APASS)) {
          if ((pFileStarImage->REFNumber != pFileStarImage2->REFNumber) &&
              (GetREFType(pFileStarImage->REFNumber) != GetREFType(pFileStarImage2->REFNumber)) &&
              (pFileStarImage->REFNumber != REF_TYPE_NONE) &&
              (pFileStarImage2->REFNumber != REF_TYPE_NONE)) {
            /* Here we want to replace pFileStarImage->catalogNumber with pFileStarImage2->catalogNumber */
#if 0
            if (substitutionCount == 79) {
              printf("line %4d At %d count %3d\n",__LINE__,substitutionCount,debugCount);
              if (debugCount == 312) {
                printf("At debugCount %d\n",debugCount);
              }
              debugCount++;
            }
#endif
            for (substitutionIndex = 0; substitutionIndex < substitutionCount; substitutionIndex++) {
              pSubstitution = &substitutionTable[substitutionIndex];
              if ((pSubstitution->badREFNumber == pFileStarImage->REFNumber) &&
                  (pSubstitution->goodREFNumber == pFileStarImage2->REFNumber)) {
                break; /* Already in the table */
              }
            }
            if (substitutionIndex >= substitutionCount) {
              if ((substitutionCount+1) >= cacheStarAlloc) {
                printf("ERROR: substitutionTable overrun\n");
                exit(-1);
              }
              pSubstitution = &substitutionTable[substitutionCount];
              pSubstitution->badREFNumber = pFileStarImage->REFNumber;
              pSubstitution->goodREFNumber = pFileStarImage2->REFNumber;
              substitutionCount++;
            }
          }

#if 0
          pFileStarImage2->REFNumber = pFileStarImage->REFNumber;
          REFNumberReplacementCount++;
#endif
        }
#if 0
        if (debugPrint != 0) {
          int coverage_bin_index;
          int DecBin;
          int RABin;
          coverage_bin_index = GetGSCBin(pGscBin01,pFileStarImage->ra,pFileStarImage->dec,&DecBin,&RABin,"combine_photometry");

          printf("DUMP line %4d gsc_bin_index %d seriesId %2d %3s plateNumber %d NUMBER %d catalogNumber %d, coverage_bin %d change %5d REFNumber %lld\n",
                 __LINE__,
                 pFileStarImage->gsc_bin_index,
                 pFileStarImage->seriesId,
                 GetSeriesString(pFileStarImage->seriesId,0),
                 pFileStarImage->plateNumber,
                 pFileStarImage->NUMBER,
                 pFileStarImage->catalogNumber,
                 coverage_bin_index,
                 REFNumberReplacementCount,
                 pFileStarImage->REFNumber);
        }      
#endif

      }
#if 0
      for (substitutionIndex = 0; substitutionIndex < substitutionCount; substitutionIndex++) {
        pSubstitution = &substitutionTable[substitutionIndex];
        printf("substitutionIndex %5d badREFNumber %16lld goodREFNumber %16lld\n",
               substitutionIndex,
               pSubstitution->badREFNumber,
               pSubstitution->goodREFNumber);
      }

#endif
      /* Now perform the substitutions */
      for (cacheStarIndex = 0; cacheStarIndex < cacheStarCount; cacheStarIndex++) {
        pFileStarImage =  &cacheStarImages[cacheStarIndex];
        for (substitutionIndex = 0; substitutionIndex < substitutionCount; substitutionIndex++) {
          pSubstitution = &substitutionTable[substitutionIndex];
          if (pSubstitution->badREFNumber == pFileStarImage->REFNumber)  {
            pFileStarImage->REFNumber = pSubstitution->goodREFNumber;
            REFNumberReplacementCount++;
            break;
          }
        }
      }
      

      /* Write out the file */
      dir1 = curGscBinIndex / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
      dir2 = curGscBinIndex / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
      sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,curGscBinIndex);
      statResult = stat(filename,&filestats);
      if (statResult == 0) {
        remove(filename);
      }


      FileSemaphoreWait(&fileSemaphore,catalogNumber[2]);
      fileNumber = open(filename,O_CREAT|O_RDWR|O_APPEND,S_IRWXU|S_IRWXG|S_IROTH);
      if (fileNumber < 0) {
        if (errno == ENOENT) {
          /* We need to create the directory */
          sprintf(cmdstr,"mkdir -p %s/mag%03d/mag%03d",photfilebase,dir1,dir2);
          result = invokeCmd(cmdstr,&timeSpent);
          if (result != 0) {
            printf("ERROR: failed to create directory %s\n",cmdstr);
            FileSemaphorePost(fileSemaphore);

            exit(-1);
          }
          sprintf(chgrpname,"%s/mag%03d",photfilebase,dir1);
          ChangeGroup(fileSemaphore,chgrpname);
          sprintf(chgrpname,"%s/mag%03d/mag%03d",photfilebase,dir1,dir2);
          ChangeGroup(fileSemaphore,chgrpname);

          directoryCreates++;
          fileNumber = open(filename,O_CREAT|O_RDWR|O_APPEND,S_IRWXU|S_IRWXG|S_IROTH);
          if (fileNumber < 0) {
            printf("ERROR: Exiting (1) with result %d errno %d %s %s\n",fileNumber,errno,strerror(errno),filename);
            FileSemaphorePost(fileSemaphore);

            exit(-1);
          }

        } else {

          printf("ERROR: Exiting (2) with result %d errno %d %s %s\n",fileNumber,errno,strerror(errno),filename);
          FileSemaphorePost(fileSemaphore);
          exit(-1);
        }
      }
   
      writeBytes = write(fileNumber,cacheStarImages,sizeof(FILESTARIMAGE)*cacheStarCount);
      if (writeBytes != (sizeof(FILESTARIMAGE)*cacheStarCount)) {
        printf("Wrote only %d bytes to for gsc_base_index %d\n",writeBytes,curGscBinIndex);
      }
      close(fileNumber);
      ChangeGroup(fileSemaphore,filename);

      FileSemaphorePost(fileSemaphore);
    }
  }
  if (combineLimitingFiles) {
    printf("Combining the limiting magnitude files\n");
    memset(pGalaxyCommon0,0,sizeof(GALAXYCOMMON));
    memset(pGalaxyCommon1,0,sizeof(GALAXYCOMMON));
    memset(pGalaxyCommon2,0,sizeof(GALAXYCOMMON));
    if (OpenGalaxyFiles(pGalaxyCommon0,qualifier[0],NULL) != 0) {
      printf("ERROR failed to open galaxy common 0\n");
      exit(-1);
    }
    if (OpenGalaxyFiles(pGalaxyCommon1,qualifier[1],NULL) != 0) {
      printf("ERROR failed to open galaxy common 1\n");
      exit(-1);
    }
    if ((pGalaxyCommon0->pLimitingGscBin->total_gsc_bins != pGalaxyCommon1->pLimitingGscBin->total_gsc_bins)) {
      printf("ERROR catalog size %d and %d disagree\n",pGalaxyCommon0->pLimitingGscBin->bin_size,pGalaxyCommon1->pLimitingGscBin->bin_size);
      exit(-1);
    }
    catalogIndex = 2;
    catalogHandle[catalogIndex] = Open(catalogname[catalogIndex],"w");
    if (catalogHandle[catalogIndex] == NULL) {
      printf("ERROR Could not open catalog file %s\n",catalogname[catalogIndex]);
      errorFlag = 1;
    } 
    indexHandle[catalogIndex] = Open(indexname[catalogIndex],"w");
    if (indexHandle[catalogIndex] == NULL) {
      printf("ERROR Could not open catalog index file %s\n",indexname[catalogIndex]);
      errorFlag = 1;
    } 


    /* Now enter the main loop */
    minGscBinIndex = 0;
    maxGscBinIndex = pGalaxyCommon0->pLimitingGscBin->total_gsc_bins;
    
    for (curGscBinIndex = minGscBinIndex; curGscBinIndex < maxGscBinIndex; curGscBinIndex ++) {
      if (GetBinCenter(pGalaxyCommon0->pLimitingGscBin,curGscBinIndex,&limitingra,&limitingdec,"combine_photometry")) {
        printf("ERROR: GetBinCenter failed for index %d\n",curGscBinIndex);
        exit(-1);       
      }
      haveGalaxyTable0 = LoadGalaxyTable(pGalaxyCommon0,limitingra,limitingdec);
      haveGalaxyTable1 = LoadGalaxyTable(pGalaxyCommon1,limitingra,limitingdec);
      if ((haveGalaxyTable0 == 0)  || (haveGalaxyTable1 == 0)) {
        printf("Failed to load galaxy table %d %d for index %d\n",haveGalaxyTable0,haveGalaxyTable1,curGscBinIndex);
        exit(-1);
      }
      totalPlateCount = pGalaxyCommon0->plateCount+pGalaxyCommon1->plateCount;
      if (totalPlateCount >= allocPlateCount) {
        if (pPlateLimitingBuffer != NULL) {
          free(pPlateLimitingBuffer);
        }
        allocPlateCount = totalPlateCount + LIMITING_BUFFER_INCREMENT;
        pPlateLimitingBuffer = (PPLATELIMITINGREC)calloc(allocPlateCount,sizeof(PLATELIMITINGREC));
        if (pPlateLimitingBuffer == NULL) {
          printf("ERROR: failed to allocate pPlateLimitingBuffer of size %d\n",allocPlateCount);
          exit(-1);
        }
      }
      totalPlateCount = 0;
      for (plateIndex = 0; plateIndex < pGalaxyCommon0->plateCount; plateIndex++) {
        pPlateLimitingRec0 =  &pGalaxyCommon0->plateLimitingBuffer[plateIndex];
        pPlateLimitingRec2 = &pPlateLimitingBuffer[totalPlateCount];
        memcpy(pPlateLimitingRec2,pPlateLimitingRec0,sizeof(PLATELIMITINGREC));
        totalPlateCount++;
        if (totalPlateCount >= allocPlateCount) {
          printf("ERROR: line %d insufficient allocation\n",__LINE__);
          exit(-1);
        }
      }
      for (plateIndex = 0; plateIndex < pGalaxyCommon1->plateCount; plateIndex++) {
        pPlateLimitingRec1 =  &pGalaxyCommon1->plateLimitingBuffer[plateIndex];
        pPlateLimitingRec2 = &pPlateLimitingBuffer[totalPlateCount];
        memcpy(pPlateLimitingRec2,pPlateLimitingRec1,sizeof(PLATELIMITINGREC));
        totalPlateCount++;
        if (totalPlateCount >= allocPlateCount) {
          printf("ERROR: line %d insufficient allocation\n",__LINE__);
          exit(-1);
        }
      }
      if (totalPlateCount > 1) {
#if 0
        if (curGscBinIndex == 56) {
          printf("At bin %d\n",curGscBinIndex);
          debugPrint = 1;
        }
        if (debugPrint) {
          printf("line %4d curGscBinIndex %d totalPlateCount %6d %6d %6d\n",__LINE__,curGscBinIndex,pGalaxyCommon0->plateCount,pGalaxyCommon1->plateCount,totalPlateCount);
        }
#endif
        qsort((void *)pPlateLimitingBuffer,totalPlateCount,sizeof(PLATELIMITINGREC),LimitingMagCompare);
        /* Now eliminate duplicates, always choosing the deepest limiting magnitude */
        plateIndexRead = 0;
        plateIndexWrite = 0;
        plateIndexSaved = 0;
        pPlateLimitingRecRead = &pPlateLimitingBuffer[plateIndexRead];
        pPlateLimitingRecSaved = pPlateLimitingRecRead;
        pPlateLimitingRecWrite = &pPlateLimitingBuffer[plateIndexWrite];
        savedJulianDate = pPlateLimitingRecRead->geoJulianDate;
        savedSeriesId = pPlateLimitingRecRead->seriesId;
        savedPlateNumber = pPlateLimitingRecRead->plateNumber;
        savedMosaicNumber = pPlateLimitingRecRead->mosaicNumber;
        savedlimiting_mag_local = pPlateLimitingRecRead->limiting_mag_local;
#if 0
        if (debugPrint) {
          printf("line %4d geoJulianDate %f seriesId %2d plateNumber %5d mosaicNumber %5d limiting_mag_local %4.1f\n",
                 __LINE__,
                 pPlateLimitingRecRead->geoJulianDate,
                 pPlateLimitingRecRead->seriesId,
                 pPlateLimitingRecRead->plateNumber,
                 pPlateLimitingRecRead->mosaicNumber,
                 pPlateLimitingRecRead->limiting_mag_local);
        }
#endif
        while (plateIndexRead < (totalPlateCount-1)) {
          plateIndexRead++;
          pPlateLimitingRecRead = &pPlateLimitingBuffer[plateIndexRead];
#if 0
          if (debugPrint) {
            printf("line %4d geoJulianDate %f seriesId %2d plateNumber %5d mosaicNumber %5d limiting_mag_local %4.1f\n",
                   __LINE__,
                   pPlateLimitingRecRead->geoJulianDate,
                   pPlateLimitingRecRead->seriesId,
                   pPlateLimitingRecRead->plateNumber,
                   pPlateLimitingRecRead->mosaicNumber,
                   pPlateLimitingRecRead->limiting_mag_local);
          }
#endif
          if ((pPlateLimitingRecRead->geoJulianDate == savedJulianDate) &&
              (pPlateLimitingRecRead->seriesId == savedSeriesId) &&
              (pPlateLimitingRecRead->plateNumber == savedPlateNumber) &&
              (pPlateLimitingRecRead->mosaicNumber == savedMosaicNumber)) {
#if 0
            if (debugPrint) {
              printf("line %4d MATCH!\n",__LINE__);
            }
#endif
            if (pPlateLimitingRecRead->limiting_mag_local > savedlimiting_mag_local) {
              plateIndexSaved =  plateIndexRead;
#if 0
              if (debugPrint) {
                printf("line %4d Deeper plateIndexSaved %d\n",__LINE__,plateIndexSaved);
              }
#endif
              pPlateLimitingRecSaved = pPlateLimitingRecRead;
              savedJulianDate = pPlateLimitingRecRead->geoJulianDate;
              savedSeriesId = pPlateLimitingRecRead->seriesId;
              savedPlateNumber = pPlateLimitingRecRead->plateNumber;
              savedMosaicNumber = pPlateLimitingRecRead->mosaicNumber;
              savedlimiting_mag_local = pPlateLimitingRecRead->limiting_mag_local;
            }
          } else {
            /* new plate, save old */
            if (plateIndexSaved != plateIndexWrite) {
              pPlateLimitingRecWrite = &pPlateLimitingBuffer[plateIndexWrite];
              memcpy(pPlateLimitingRecWrite,pPlateLimitingRecSaved,sizeof(PLATELIMITINGREC));
            }
#if 0
            if ((plateIndexWrite >= 123) && (curGscBinIndex == 1389)) {
              printf("At plateIndexWrite %d, bin %d\n",plateIndexWrite,curGscBinIndex);
            }
            if (debugPrint) {
              printf("line %4d geoJulianDate %f seriesId %2d plateNumber %5d mosaicNumber %5d limiting_mag_local %4.1f plateIndexWrite %5d\n",
                     __LINE__,
                     pPlateLimitingRecWrite->geoJulianDate,
                     pPlateLimitingRecWrite->seriesId,
                     pPlateLimitingRecWrite->plateNumber,
                     pPlateLimitingRecWrite->mosaicNumber,
                     pPlateLimitingRecWrite->limiting_mag_local,
                     plateIndexWrite);
            }
#endif
            plateIndexWrite++;
            plateIndexSaved = plateIndexRead;
            pPlateLimitingRecSaved = pPlateLimitingRecRead;
            savedJulianDate = pPlateLimitingRecRead->geoJulianDate;
            savedSeriesId = pPlateLimitingRecRead->seriesId;
            savedPlateNumber = pPlateLimitingRecRead->plateNumber;
            savedMosaicNumber = pPlateLimitingRecRead->mosaicNumber;
            savedlimiting_mag_local = pPlateLimitingRecRead->limiting_mag_local;
          }
        }
        /* save the last */
        if (plateIndexSaved != plateIndexWrite) {
          pPlateLimitingRecWrite = &pPlateLimitingBuffer[plateIndexWrite];
          memcpy(pPlateLimitingRecWrite,pPlateLimitingRecSaved,sizeof(PLATELIMITINGREC));
        }
#if 0
        if (debugPrint) {
          printf("line %4d geoJulianDate %f seriesId %2d plateNumber %5d mosaicNumber %5d limiting_mag_local %4.1f plateIndexWrite %5d\n",
                 __LINE__,
                 pPlateLimitingRecWrite->geoJulianDate,
                 pPlateLimitingRecWrite->seriesId,
                 pPlateLimitingRecWrite->plateNumber,
                 pPlateLimitingRecWrite->mosaicNumber,
                 pPlateLimitingRecWrite->limiting_mag_local,
                 plateIndexWrite);
        }
#endif
        plateIndexWrite++;
        totalPlateCount = plateIndexWrite;
        /* Sanity check */
        for (plateIndex = 1; plateIndex < totalPlateCount; plateIndex++) {
          pPlateLimitingRecRead = &pPlateLimitingBuffer[plateIndex];
          pPlateLimitingRecSaved = &pPlateLimitingBuffer[plateIndex-1];
#if 0
          if (debugPrint) {
            printf("line %4d geoJulianDate %f seriesId %2d plateNumber %5d mosaicNumber %5d limiting_mag_local %4.1f plateIndexRead %5d\n",
                   __LINE__,
                   pPlateLimitingRecRead->geoJulianDate,
                   pPlateLimitingRecRead->seriesId,
                   pPlateLimitingRecRead->plateNumber,
                   pPlateLimitingRecRead->mosaicNumber,
                   pPlateLimitingRecRead->limiting_mag_local,
                   plateIndex);
          }
#endif
          if ((pPlateLimitingRecRead->geoJulianDate ==pPlateLimitingRecSaved->geoJulianDate ) &&
              (pPlateLimitingRecRead->seriesId == pPlateLimitingRecSaved->seriesId ) &&
              (pPlateLimitingRecRead->plateNumber == pPlateLimitingRecSaved->plateNumber) &&
              (pPlateLimitingRecRead->mosaicNumber == pPlateLimitingRecSaved->mosaicNumber)) {
            printf("ERROR: line %d unexpected duplicate plate date %f seriesId %d plateNumber %d \n",
                   __LINE__,
                   pPlateLimitingRecRead->geoJulianDate,
                   pPlateLimitingRecRead->seriesId,
                   pPlateLimitingRecRead->plateNumber);
            exit(-1);
          }
        }

      }
      /* Now write out the results */
      /* First write out all of the records for variables and galaxies */
      totalGalaxyWrite = 0;
      byteCount = 0;
      for (objectIndex = 0; objectIndex < pGalaxyCommon0->galaxyCount; objectIndex++) {
        pGalaxyRec = &pGalaxyCommon0->galaxyInputBuffer[objectIndex];
        writeItems = Write(catalogHandle[catalogIndex],pGalaxyRec,sizeof(GALAXYREC),1);
				if (writeItems != 1) {
					printf("ERROR writing a galaxy/variable record\n");
					exit(-1);
				}
        totalGalaxyWrite++;
        byteCount += sizeof(GALAXYREC);
      }
      writeItems = Write(catalogHandle[catalogIndex],pPlateLimitingBuffer,sizeof(PLATELIMITINGREC),totalPlateCount);
      if (writeItems != totalPlateCount) {
        printf("ERROR writing %d of %d items to the limiting magnitude table\n",writeItems,totalPlateCount);
        exit(-1);
      }
      pStarIndex->offset = platerec_bytes;
      pStarIndex->binNumber = curGscBinIndex;
      platerec_bytes += (writeItems*sizeof(PLATELIMITINGREC)) + (totalGalaxyWrite *sizeof(GALAXYREC));
      pStarIndex->numStars = (writeItems*sizeof(PLATELIMITINGREC)) + (totalGalaxyWrite *sizeof(GALAXYREC));

   		writeItems = Write(indexHandle[catalogIndex],pStarIndex,sizeof(STARINDEX),1);
      if (writeItems != 1) {
        fprintf(stderr,"ERROR writing the index file\n");
        exit(-1);
      }
   

      
    }        
  }
 
  /* STOP HERE */

  totalMagnitudeAlloc += magnitudeAlloc;
  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds, %d queries, %d totalStars, %d recoveredStars, %d badModulusBin, %d magnitudeAlloc %d gscBins, %d total magnitudes\n",curTime,GetQueryCount(),totalStars,recoveredStars,badModulusBin,totalMagnitudeAlloc,binCount,totalMagnitudes);

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

  printf("substitutionCount %d REFNumberReplacementCount %d directoryCreates %d\n",substitutionCount,REFNumberReplacementCount,directoryCreates);
  fprintf(outHandle,"substitutionCount %d REFNumberReplacementCount %d directoryCreates %d\n",substitutionCount,REFNumberReplacementCount,directoryCreates);


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
  if (debug_handle != NULL) {
    fclose(debug_handle);
    debug_handle = NULL;
  }
  FreeFileCommon(pFileCommon0,1);
  FreeFileCommon(pFileCommon1,1);
  FreeFileCommon(pFileCommon2,1);
  for (catalogIndex = 0; catalogIndex < 3; catalogIndex++) {
    if (catalogHandle[catalogIndex] != NULL) {
      Close(catalogHandle[catalogIndex]);
    }
    if (indexHandle[catalogIndex] != NULL) {
      Close(indexHandle[catalogIndex]);
    }
  }


  fclose(outHandle);
  return(0);
}

