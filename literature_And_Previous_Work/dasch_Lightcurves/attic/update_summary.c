// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* update_summary.c
 *
 * This routine supercedes computemag.c to compute the DASCH median magnitude for all stars that need updating.
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   update_summary.c pipelineutils.a -ltable -lutil  -lwcs -o update_summary -L/usr/lib64/mysql -lmysqlclient
 *
 * Mar  2, 2009 Edward J. Los - Original Version
 * May  1, 2009 Edward J. Los - Improve error handling
 * May  3, 2009 Edward J. Los - Do not abort on gsc catalog read errors
 * May  8, 2009 Edward J. Los - Implement a Malmquist test
 * May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 * Jun  8, 2009 Edward J. Los - Convert to new database format
 * Jul  1, 2009 Edward J. Los - Add MAGflag, class, RaPM, DecPM to STARENTRY
 * Jul 15, 2009 Edward J. Los - Back out fatal error if no GSC entry can be found 
 * Jul 19, 2009 Edward J. Los - Correct ChargeTime indices
 * Jul 28, 2009 Edward J. Los - Enhance performance by reading only the necessary data for each magnitude
 *                              and by processing multiple stars with each MySQL access.
 * Aug 18, 2009 Edward J. Los - Add Damon_factor
 *                              Always test for version ID
 * Oct 19, 2009 Edward J. Los - Covert the magnitudes tables to files
 * Dec  2, 2009 Edward J. Los - Store the magnitudeFileFlag in the photglobal table
 * Dec 29, 2009 Edward J. Los - Accept only points from plates that have no quality bits set
 *                              and that have a colorterm within range
 * Feb 26, 2010 Edward J. Los - Add Sextractor_Blend, the number of good points with the Sextractor blend or neighbors bits set
 * Apr 16, 2010 Edward J. Los - Add option to exclude four series with known fitting problems (-e qualifier).
 * Apr 30, 2010 Edward J. Los - Correct the drad calculation to include proper motion
 * May 18, 2010 Edward J. Los - Do not clip the dradmed and dradrms calculation.
 *                              Add max_drad for good points
 *                              Add ngoodB, dradB, max_dradB, and Malmquist_factorB for points ignoring FILTER_AMASK_PLOT in AFLAGS
 * Aug  4, 2010 Edward J. Los - Add new parameters from Sumin Tang request of  7/03/10 10:53.  Stars with high dra and ddec and Sextractor blends are now rejected.
 *                              Ignore the GSC_VARIABLE_BIT when creating star summary statistics
 * Dec 31, 2010 Edward J. Los - display minVersionId
 * Mar 23, 2011 Edward J. Los - Add apass catalog support
 * May 10, 2011 Edward J. Los - Add slope parameters from Sumin Tang's memorandum of Thu 4/28/11 4:37 PM
 * Jun  3, 2011 Edward J. Los - Preserve the versionId in the stars table.
 * Jun 13, 2011 Edward J. Los - Move the spatial bin plate quality check before the blended image selection.
 *                              Preserve nDamonBlue and nNonDamonBlue counts during clipping
 * Nov  9, 2011 Edward J. Los - Add magnitude-dependent parameters (magdep_bin, magcal_magdep, and magcal_magdep_rms).
 * Nov 14, 2011 Edward J. Los - make minVersionId dependent on the catalog
 * Dec  9, 2011 Edward J. Los - Reject observations that do not have magnitude-dependent correction
 * Dec 16, 2011 Edward J. Los - Test new limiting magnitude criteria of Sumin's memo of Thu 12/15/11 1:52 PM
 * Jul 18, 2012 Edward J. Los - Make sure that all calls to GetRef do not return spaces in object names
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Aug 13, 2012 Edward J. Los - Change "Execution" to "Intermediate" for interim report.
 * Sep 28, 2012 Edward J. Los - Move ProcessStar and associated routes to photometryutils and rename it ProcessLightcurveParameters
 *                            - Calculate lightcurve parameters for unmatched objects.
 * Dec  8, 2012 Edward J. Los - Make the star selection more permissive
 * Mar 13, 2013 Edward J. Los - Make the routine MP safe by adding the qualifier to the output files
 * Mar 21, 2017 Edward J. Los - Deprecate this program in favor of update_summary2.
 * May 28, 2018 Edward J. Los - support gaia
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 *
 *   update_summary   -n 2  -i -o /dasch/Pipeline/ingest/file_update_summary.log
 *   echo "update_summary   -n 2  -i -o /dasch/Pipeline/ingest/file_update_summary.log " | at now
 *
 *    update_summary  -q kepler  -n 2  -i -o /dasch/Pipeline/ingest/file_kepler_update_summary.log
 *    echo "update_summary  -q kepler  -n 2  -i -o /dasch/Pipeline/ingest/file_kepler_update_summary.log" | at now
 *
 *    update_summary  -q apass  -n 2  -i -o /dasch/Pipeline/ingest/file_apass_update_summary.log
 *    echo "update_summary  -q apass  -n 2  -i -o /dasch/Pipeline/ingest/file_apass_update_summary.log " | at now
 *
 *    update_summary  -q experimental  -n 2  -i -o /dasch/Pipeline/ingest/file_experimental_update_summary.log
 *    echo "update_summary  -q experimental  -n 2  -i -o /dasch/Pipeline/ingest/file_experimental_update_summary.log" | at now
 *
 *   UPDATE stars set updateflag = "yes";  
 *    Query OK, 94125314 rows affected (2 min 43.86 sec)
 *
 *   UPDATE stars1 set updateflag = "yes";
 *    Query OK, 3119478 rows affected (8 min 23.35 sec)
 *
 *   UPDATE stars2 set updateflag = "yes";
 *    Query OK, 88294065 rows affected (1 hour 42 min 37.20 sec)
 *
 *   UPDATE stars3 set updateflag = "yes";
 *    Query OK, 2729761 rows affected (6 min 0.48 sec)
 *
 *   Sumin's parameter test
     UPDATE stars1 set updateflag = "yes" where REFNumber = 210311340;
     UPDATE stars1 set updateflag = "yes" where REFNumber = 210311993;
     UPDATE stars1 set updateflag = "yes" where REFNumber = 210312235;  
 *   SELECT * FROM starsummary1 where REFNumber = 210561482 or  REFNumber = 211007282 or  REFNumber = 211103080;
 *
 *   On DELL: (dell.list) -n 2
 *             Execution Time: 1225 seconds, 310572 queries, 339387 totalStars, 11699 goodstars, 40 magnitudeAlloc 500 vectorAlloc, 310429 gscBins
 *
 */


#include <math.h>
#include <errno.h>
#include "table.h"
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include "pipelineutils.h"
#include "photometryutils.h"
/* #define DUMP_VALUES 1 */


#define EPS 2.2204e-16
#define MAX_BUFFER 256
#if 0
#define STAR_LIMIT     10000
#define SUMMARY_LIMIT  100
#else
#define STAR_LIMIT     50000
#define SUMMARY_LIMIT  50000
#endif
/* #define SKIP_REMOVAL 1 */
#define GSC_BIN_SPACING (250*magnitudesModulus) /* Read GSC bins in this size block */
#define MAX_ALLOCATION (0x10000000)
#define AVE_STARS 1000 /* Expected average number of objects per bin */
#define MAX_FILENAME 512

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;













#ifdef DUMP_VALUES
void DumpVector(int curCount,double *vector) {
  int index;
  for (index = 0; index < curCount; index++) {
    printf("Index %d, value %f\n",index,vector[index]);
  }

}
void DumpMVector(int curCount,PMALMQUIST vector) {
  int index;
  for (index = 0; index < curCount; index++) {
    PMALMQUIST pVector = &vector[index];
    printf("IndexM %d, value %f\n",index,pVector->magcal_magdep);
  }

}
#endif /* DUMP_VALUES */



int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  time_t deltaTime;
  time_t oldTime;
  time_t printTime = 0;
  int verbose = 0;
  int debugMode = 0;
  int magnitudeFileFlag = 0;
  int ignoreVersionId = 1;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int nvals;
  FILE *outHandle;
  char timestr[100];
  struct tm *ptr;
  int gotAnswer;
  int curMagnitudes;
  long long totalMagnitudes = 0;
  int staleVersionIdCount = 0;
  char *dotPtr;
  int processedCount = 0;

  int minGoodStars = MIN_SUMMARY_GOODSTARS;

  char *mysqlhost;
  char *username;
  char *password;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;
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
  long long oldTotalMagnitudes = 0;
  PARAMETERVECTORSTORE vectorStore;
  PPARAMETERVECTORSTORE pVectorStore = &vectorStore;
#if 0
  double *vector1 = NULL; /* Vector for locally corrected magnitudes (magcal_magdep) */
  double *vector2 = NULL; /* Vector for julian days */
  double *vector3 = NULL; /* Vector for lowess magnitudes (magcal_iso) */
  double *pVectorStore->vector[VECTOR_MAGCAL_LOCAL_RMS] = NULL; /* Vector for magcal_local_rms */
  double *vector5 = NULL; /* Vector for drad estimate */
  double *vector5B = NULL; /* Vector for dradB estimate */
  PMALMQUIST vector6 = NULL; /* Vector for the Malmquist_factor */
  PMALMQUIST vector6B = NULL; /* Vector for the Malmquist_factor */
  double *vector7 = NULL; /* Vector for non-Damon */
  double *vector8 = NULL; /* Vector for Damon */
#endif
  int minGscBinIndex;
  int maxGscBinIndex;
  int highGscBinIndex;
  int curGscBinIndex;
  int curGscModulus;
  int maxGscBinIncrement;
#if 0
  int vectorAlloc = 0;
#endif
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
  
  int excludeSeriesIndex;
  int magnitudesModulus =  MAGNITUDES_MODULUS;
  int dumpAllFlag = 0;
  FILE *headerHandle = NULL;
#if 1
  printf("ERROR: update_summary is deprecated because the stars table is no longer updated.  Use update_summary2 instead\n");
  exit(-1);
#endif


  memset(pVectorStore,0,sizeof(PARAMETERVECTORSTORE));

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

#ifdef NEW_LIMITING_MAG
#if 0
  printf("ERROR: NEW_LIMITING_MAG is set\n");
#endif
#endif /* NEW_LIMITING_MAG */
  catalogString[0] = 0;
  ChargeTime(RESET_ALL_ENTRY,&curUsec,deltaTable,totalTable);
  SetQueryCount(0);

  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_INGEST is not defined\n");
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
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(outfile) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        


        case 'n': /* Minimum number of good stars */
        case 'N':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&minGoodStars);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode minGoodStars %s\n",*argv);
              errorFlag = 1;
            } else {
              if (minGoodStars < 1) {
                fprintf(stderr,"ERROR: illegal value for minGoodStars: %d  Using default %d\n",minGoodStars,MIN_SUMMARY_GOODSTARS);
                minGoodStars = MIN_SUMMARY_GOODSTARS;
              }
            }
          }
          break;
          
        case 'i':  /* Ignore the versionId */
        case 'I':
          ignoreVersionId = 1;
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



        case 'e': /* Exclude bad astrometry series */
        case 'E':
          excludeSeriesCount = NUM_EXCLUDE_SERIES;
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
    printf("Usage: update_summary options\n");
    printf("  options: -v verbose\n");
    printf("           -d debugging mode\n");
    printf("           -o <output file>\n");
    printf("           -n <minimum stars per lightcurve>\n");
    printf("           -q input catalog qualifer\n");
    printf("           -e excludes ac,ca,ax, and am series\n");


    return(-1);
  }

  maxGscBinIncrement  = magnitudesModulus * (MAX_ALLOCATION / (AVE_STARS * (sizeof(PHOTSTARIMAGE)+sizeof(FILESTARIMAGE))));

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
    fprintf(stderr,"DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }

  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    fprintf(stderr,"DASCH_PHOT_USERNAME is not defined\n");
    return(-1);
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    fprintf(stderr,"DASCH_PHOT_PASSWORD is not defined\n");
    return(-1);
  }
  mysql_init(pPhotConnection);


  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
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
  }

  


  if (magnitudeFileFlag) {
    printf("File Header size %d File Image size %d\n",sizeof(FILESTARHEADER),sizeof(FILESTARIMAGE));
    fprintf(outHandle,"File Header size %d File Image size %d\n",sizeof(FILESTARHEADER),sizeof(FILESTARIMAGE));
  } else {
    printf("ERROR: magnitudes stored in the MySQL database are no longer supported\n");
    exit(-1);

  }


  InitFileCommon(stdout,pFileCommon,pConnection,pPhotConnection,catalogNumber,1);
  InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);

  mysql_close(pConnection);
                    
  printf("update_summary of %s %s, outfile %s minimum stars %d starLimit %d summaryLimit %d sizeof(PHOTSTARIMAGE) %d maxGscBinIncrement %d excludeSeriesCount %d magnitudeFile %d minVersionId%d %d\n",
         __DATE__,__TIME__,outfile,minGoodStars,STAR_LIMIT,SUMMARY_LIMIT,sizeof(PHOTSTARIMAGE),maxGscBinIncrement,excludeSeriesCount,magnitudeFileFlag,catalogNumber,pPhotGlobal->minVersionId[catalogNumber]);
  fprintf(outHandle,"update_summary of %s %s, outfile %s minimum stars %d starLimit %d summaryLimit %d sizeof(PHOTSTARIMAGE) %d maxGscBinIncrement %d excludeSeriesCount %d magnitudeFile %d minVersionId%d %d\n",
          __DATE__,__TIME__,outfile,minGoodStars,STAR_LIMIT,SUMMARY_LIMIT,sizeof(PHOTSTARIMAGE),maxGscBinIncrement,excludeSeriesCount,magnitudeFileFlag,catalogNumber,pPhotGlobal->minVersionId[catalogNumber]);

  TimeStampPhotometryGlobal(pPhotConnection);


  time(&startTime);
  oldTime = startTime;
  printTime = startTime;
  ptr = localtime(&startTime);
  strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
 
  printf("update_summary starting at %s\n",timestr);
  fprintf(outHandle,"update_summary starting at %s\n",timestr);
  fflush(outHandle);
  pStarTable = (PSTARENTRY) calloc(STAR_LIMIT,sizeof(STARENTRY));
  if (pStarTable == NULL) {
    fprintf(stderr,"ERROR: failed to allocate pStarTable\n");
    exit(-1);
  }

  strcpy(stars_name,ingestDirectory);
  strcat(stars_name,"/");
  strcat(summary_name,qualifier);
  strcat(stars_name,"update_summary_stars.txt");

  strcpy(summary_name,ingestDirectory);
  strcat(summary_name,"/");
  strcat(summary_name,qualifier);
  strcat(summary_name,"update_summary.txt");

  strcpy(debug_name,ingestDirectory);
  strcat(debug_name,"/");
  strcat(summary_name,qualifier);
  strcat(debug_name,"update_summary.debug");
 
  summary_handle = fopen(summary_name,"wt");
  if (summary_handle == NULL) {
    fprintf(stderr,"ERROR: Failed to open summary file  %s\n",summary_name);
    exit(-1);
  }


  strcpy(galaxy_name,catalogdir);
  slashPtr = strrchr(galaxy_name,'/');
  if (slashPtr != NULL) {
    slashPtr++;
  } else {
    slashPtr = galaxy_name;
  }
  *slashPtr = 0;
  strcat(galaxy_name,"galaxy.dat");


  if (PopulateGalaxyTree(pGalaxyCommon,galaxy_name,0) != 0) {
    exit(-1);
  }


  sprintf(queryString,"update photglobal set keepalive = 'yes';");
  res = ExecuteQuery(pPhotConnection,queryString);
  if (res) {
    exit(-1);
  }
  /* Enter the main loop */
  while (1) {
    gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
    time(&curTime);
    if ((gotAnswer != 1)  || (pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES)) {
      

      ptr = localtime(&curTime);
      strftime(timestr,25,"%Y-%m-%dT%H-%M-%S", ptr);
 
      printf("update_summary aborting at %s with gotAnswer %d keepalive %d\n",
             timestr,gotAnswer,pPhotGlobal->keepalive);
      fprintf(outHandle,"update_summary aborting at %s with gotAnswer %d keepalive %d\n",
              timestr,gotAnswer,pPhotGlobal->keepalive);
      break;
    }
    curTime -= printTime;
    if (curTime > 21600) {
      time(&curTime);
      printTime = curTime;
      curTime -= startTime;
      if (pCurStarEntry2 != NULL) {
        fprintf(outHandle,"At gsc_bin_index %d and REF %s and ra %5.1f declination %5.1f\n",pCurStarEntry2->gsc_bin_index,pCurStarEntry2->REF,pCurStarEntry2->ra,pCurStarEntry2->dec);
      }
      fprintf(outHandle,"Elapsed Time: %d seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins, %lld magnitudes\n",curTime,GetQueryCount(),totalStars,goodStars,magnitudeAlloc,pVectorStore->vectorAlloc,binCount,totalMagnitudes);
      fprintf(outHandle," queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec %f  magnitude/sec:%f\n",
              (1.0 * GetQueryCount())/(1.0*curTime),
              (1.0 * totalStars)    /(1.0*curTime),
              (1.0 * goodStars )    /(1.0*curTime),
              (1.0 * binCount  )    /(1.0*curTime),
              (1.0 * totalMagnitudes)/(1.0*curTime));
 
      fprintf(outHandle,"Intermediate Time (sec) ");
      for (usecIndex = 0; usecIndex < MAX_CHARGE_ENTRY; usecIndex++) {
        double outSeconds;
        outSeconds = (1.0 * totalTable[usecIndex])/1000000.0;
        fprintf(outHandle,"%10.3f ",outSeconds);
      }
      fprintf(outHandle,"\n");
      fflush(outHandle);
    }

    ChargeTime(0,&curUsec,deltaTable,totalTable);
    curStars = GetStarEntry(pPhotConnection,&pStarTable,STAR_LIMIT,UPDATEFLAG_NEED_UPDATE,NULL,0,catalogString,1,0);
    
    ChargeTime(1,&curUsec,deltaTable,totalTable);
    if (curStars <= 0) {
      printf("GetStarEntry curStars is %d\n",curStars);
      break;
    }

    totalStars += curStars;
    if (verbose) {
      time(&curTime);
      curTime -= startTime;
      fprintf(stderr,"Got %d stars total %d at %d seconds %s\n",curStars,totalStars,curTime,pStarTable->REF);
      if (debugMode) {
        fprintf(outHandle,"Got %d stars total %d at %d seconds %s\n",curStars,totalStars,curTime,pStarTable->REF);
        fflush(outHandle);
      }

    }

    if (debugMode) {
      if (debug_handle != NULL) {
        fclose(debug_handle);
        debug_handle = NULL;
      }
      debug_handle = fopen(debug_name,"wt");
      if (debug_handle == NULL) {
        fprintf(stderr,"ERROR: Failed to open summary file  %s\n",debug_name);
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
      if (pCurStarEntry->gsc_bin_index <= 0) {
        pCurStarEntry->processedflag = 1;
        continue;
      } 

      if (debugMode) {
        time(&curTime);
        curTime -= startTime;
        fprintf(debug_handle,"starIndex %d processedCount %d getting magnitudes for %d at %d\n",starIndex,processedCount,pCurStarEntry->gsc_bin_index,curTime);
      } 
      /* Here we look for contiguous runs of the gsc_bin_index */

      minGscBinIndex = pCurStarEntry->gsc_bin_index;
      curGscBinIndex = minGscBinIndex;
      curGscModulus = minGscBinIndex % magnitudesModulus;
      if (magnitudeFileFlag == 0) {

       

        for (starIndex2 = starIndex; starIndex2 < curStars; starIndex2++) {
          if (errorFlag) {
            break;
          }
          pCurStarEntry2 = &pStarTable[starIndex2];
          if (pCurStarEntry2->processedflag != 0) {
            /* We processed this one already */
            continue;
          }
          if ((pCurStarEntry2->gsc_bin_index % magnitudesModulus) != 0) {
            continue;
          }
          if ((pCurStarEntry2->gsc_bin_index - GSC_BIN_SPACING) > curGscBinIndex) {
            break;
          }
          curGscBinIndex = pCurStarEntry2->gsc_bin_index;
        }
        maxGscBinIndex = curGscBinIndex;


        if ((maxGscBinIndex - minGscBinIndex) > maxGscBinIncrement) {
          maxGscBinIndex = minGscBinIndex + maxGscBinIncrement;
        }


        curGscModulus = minGscBinIndex % magnitudesModulus;
        if (maxGscBinIndex == minGscBinIndex) {
          highGscBinIndex = 0;
        } else {
          highGscBinIndex = maxGscBinIndex;
        }

      } else {
        maxGscBinIndex = ((minGscBinIndex + MAG_FILE_MODULUS) / MAG_FILE_MODULUS);
        maxGscBinIndex = (maxGscBinIndex * MAG_FILE_MODULUS) -1;
        highGscBinIndex = maxGscBinIndex;
      }

      ChargeTime(0,&curUsec,deltaTable,totalTable);
      if (magnitudeFileFlag == 0) {
        exit(-1);
      } else {
        curMagnitudes = GetFileSummaryMagnitudes(pGscBin,pFileCommon,minGscBinIndex,highGscBinIndex,0,1,0,dumpAllFlag,headerHandle,catalogString,0,NULL,0,1);
        /* For now, just copy the results over */
        if (magnitudeAlloc < curMagnitudes) {
          free(pMagnitudeTable);
          magnitudeAlloc = curMagnitudes;
          pMagnitudeTable = (PPHOTSTARIMAGE)calloc(magnitudeAlloc,sizeof(PHOTSTARIMAGE));
          if (pMagnitudeTable == NULL) {
            fprintf(stderr,"ERROR: failed to allocate %d PHOTSTARIMAGEs in GetSummaryMagnitudes\n",magnitudeAlloc);
            exit(-1);
          }
        }
        for (magnitudeIndex = 0; magnitudeIndex < curMagnitudes;magnitudeIndex++) {
          pCurSumstarimage = &pMagnitudeTable[magnitudeIndex];
          pFileStarImage = &pFileCommon->magnitudeBuffer[magnitudeIndex];
          memset(pCurSumstarimage,0,sizeof(PHOTSTARIMAGE));
          pCurSumstarimage->pFileStarImage = pFileStarImage;
          GetREF(pFileStarImage->REFNumber,pCurSumstarimage->REF,0,1);
          strcpy(pCurSumstarimage->series,GetSeriesString(pFileStarImage->seriesId,1));
          if (GetREF(pFileStarImage->REFNumber,pCurSumstarimage->REF,0,0) != 0) {
            fprintf(stderr,"ERROR: corrupt REFNumber at magnitude index %d magnitudeIndex,bin %d to %d\n",
                    magnitudeIndex,minGscBinIndex,highGscBinIndex);
            exit(-1);
          }
        }

      }
      totalMagnitudes += curMagnitudes;
      ChargeTime(2,&curUsec,deltaTable,totalTable);

      for (magnitudeIndex = 0; magnitudeIndex < curMagnitudes;magnitudeIndex++) {
        pCurSumstarimage = &pMagnitudeTable[magnitudeIndex];
        pFileStarImage = &pFileCommon->magnitudeBuffer[magnitudeIndex];
        seriesId = pFileStarImage->seriesId;
        if ((seriesId < 0) || (seriesId > MAX_SERIES)) {
          printf("ERROR: illegal seriesId for seriesId %d, plateNumber %d, NUMBER %lld\n",seriesId,pFileStarImage->plateNumber,pFileStarImage->REFNumber);
          exit(-1);
        }
        plateNumber = pFileStarImage->plateNumber;
        if ((plateNumber < 0) || (plateNumber > pFileCommon->maxPlateNumber[seriesId])) {
          printf("ERROR: illegal plateNumber for seriesId %d, plateNumber %d, NUMBER %lld\n",seriesId,pFileStarImage->plateNumber,pFileStarImage->REFNumber);
          exit(-1);
        }
        seriesVersionArray = pFileCommon->seriesVersionTable[seriesId];
        if (seriesVersionArray == NULL) {
          seriesVersionArray = (int *)calloc(pFileCommon->maxPlateNumber[seriesId]+1,sizeof(int));
          if (seriesVersionArray == NULL) {
            printf("ERROR: failed to allocate seriesVersion Array of size %d\n",pFileCommon->maxPlateNumber[seriesId]+1);
            exit(-1);
          }
          pFileCommon->seriesVersionTable[seriesId] = seriesVersionArray;
        }
        if (seriesVersionArray[plateNumber] < pFileStarImage->versionId) {
          seriesVersionArray[plateNumber] = pFileStarImage->versionId;
        }
      }


      if (debugMode) {
        time(&curTime);
        curTime -= startTime;
        fprintf(debug_handle,"starIndex %d finished getting magnitudes for %d at %d\n",starIndex,pCurStarEntry->gsc_bin_index,curTime);
      }
      if (curMagnitudes < 0) {
        errorFlag = 1;
        break;
      }

      binCount++;
      curGscBinIndex = minGscBinIndex;
      for (starIndex2 = starIndex; starIndex2 < curStars; starIndex2++) {
        if (errorFlag) {
          break;
        }
        pCurStarEntry2 = &pStarTable[starIndex2];
        if (pCurStarEntry2->processedflag != 0) {
          /* We processed this one already */
          continue;
        }
        
        if ((magnitudeFileFlag == 0) && ((pCurStarEntry2->gsc_bin_index % magnitudesModulus) != curGscModulus)) {
          /* Not in our partition of the database */
          continue;
        }
        if (pCurStarEntry2->gsc_bin_index > maxGscBinIndex) {
          /* Need to read another block of magnitudes */
          break;
        } else {

          pCurStarEntry2->processedflag = 1;
          processedCount++;




          if (cur_gsc_bin_index != pCurStarEntry2->gsc_bin_index) {
            cur_gsc_bin_index = pCurStarEntry2->gsc_bin_index;
            if (cur_gsc_bin_index >= pGscBin->total_gsc_bins) {
              printf("ERROR: gsc_bin_index %d exceeds total %d in update_summary for REFNUMBER %lld\n",
                     cur_gsc_bin_index,
                     pGscBin->total_gsc_bins,
                     pCurStarEntry2->REFNumber);
              cur_gsc_bin_index = -1;
              memset(pCurStarIndex,0,sizeof(STARINDEX));
            } else {  

              Seek(indexHandle,cur_gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
              readItems = Read(indexHandle,pCurStarIndex,sizeof(STARINDEX),1);
              if (readItems != 1) {
                fprintf(stderr,"ERROR reading star index file errno: %d %s REF %s ra %f dec %f gsc_bin_index %d\n",
                        errno,
                        strerror(errno),
                        pCurStarEntry2->REF,
                        pCurStarEntry2->ra,
                        pCurStarEntry2->dec,
                        pCurStarEntry2->gsc_bin_index);
                cur_gsc_bin_index = -1;
                memset(pCurStarIndex,0,sizeof(STARINDEX));
              } else {
          
                if (pCurStarIndex->numStars > gscImageAlloc) {
                  gscImageAlloc = pCurStarIndex->numStars+1000;
                  pCurGscImage = realloc(pGscImageTable,gscImageAlloc*sizeof(GSCIMAGE));
                  if (pCurGscImage == NULL) {
                    fprintf(stderr,"ERROR: failed to realloc pGscImageTable of size %d\n",gscImageAlloc*sizeof(GSCIMAGE));
                    exit(-1);
                  }
                  pGscImageTable = pCurGscImage;
                  pCurGscImage = NULL;
                }
                Seek(catalogHandle,pCurStarIndex->offset,SEEK_SET);
                readItems = Read(catalogHandle,pGscImageTable,sizeof(GSCIMAGE),pCurStarIndex->numStars);
                if (readItems != pCurStarIndex->numStars) {
                  fprintf(stderr,"ERROR reading gsc catalog file\n");
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
                break;
              }

            }
            if (gscIndex == pCurStarIndex->numStars) {
#if 0
              fprintf(stderr,"WARNING: no gsc entry found in bin %d for star %s ra %f dec %f\n",
                      cur_gsc_bin_index,
                      pCurStarEntry2->REF,
                      pCurStarEntry2->ra,
                      pCurStarEntry2->dec);
              exit(-1);
#endif
              noGscEntryCount++;
#if 0
              pCurStarEntry2->ra = -1.0;
              pCurStarEntry2->dec = -91.0;
              pCurStarEntry2->MAGFlag = 0;
              pCurStarEntry2->class = 0;
              pCurStarEntry2->VFlag = 0;
              pCurStarEntry2->RaPM = 0;
              pCurStarEntry2->DecPM = 0;
#endif
            }



            goodStars += ProcessLightcurveParameters(stdout,pFileCommon,pPhotConnection,pMagnitudeTable,curMagnitudes,pCurStarEntry2,pVectorStore,minGoodStars,pPhotGlobal->currentVersion,summary_handle,&summaryCount,debugMode,&staleVersionIdCount,damonSeriesId,excludeSeriesId,excludeSeriesCount,NULL,cur_gsc_bin_index);
            if (summaryCount >= SUMMARY_LIMIT) {
              fclose(summary_handle);
              if (debugMode) {
                fprintf(outHandle,"LOAD DATA LOCAL INFILE '%s' REPLACE INTO TABLE starsummary%s;",summary_name,catalogString);
                fflush(outHandle);
              }
              sprintf(queryString,"LOAD DATA LOCAL INFILE '%s' REPLACE INTO TABLE starsummary%s;",summary_name,catalogString);
              ChargeTime(0,&curUsec,deltaTable,totalTable);
              res = ExecuteQuery(pPhotConnection,queryString);
              ChargeTime(3,&curUsec,deltaTable,totalTable);
              if (debugMode) {
                fprintf(outHandle,"LOAD DATA LOCAL into starsummary%s completed\n",catalogString);
                fflush(outHandle);
              }
              if (res) {
                errorFlag = 1;
                break;
              } 
#if SKIP_REMOVAL
              fprintf(stderr,"ERROR: removal of %s skipped\n",summary_name);
#else /* SKIP_REMOVAL */
              remove(summary_name);
#endif /* SKIP_REMOVAL */
              summary_handle = fopen(summary_name,"wt");
              if (summary_handle == NULL) {
                fprintf(stderr,"ERROR: Failed to open summary file  %s\n",summary_name);
                exit(-1);
              }
              
              summaryCount = 0;

              time(&curTime);
              deltaTime = curTime - oldTime;
              oldTime = curTime;
              curTime -= startTime;
              printf("Processing      %d seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",curTime,GetQueryCount(),totalStars,goodStars,magnitudeAlloc,pVectorStore->vectorAlloc,binCount,totalMagnitudes);
              fprintf(outHandle,"Processing      %d seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",curTime,GetQueryCount(),totalStars,goodStars,magnitudeAlloc,pVectorStore->vectorAlloc,binCount,totalMagnitudes);
              if (deltaTime > 0) {
                printf(" queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec: %f magnitudes/sec %f\n",
                       (1.0 *(GetQueryCount()-oldQueryCount))/(1.0*deltaTime),
                       (1.0 *(totalStars-oldTotalStars))/(1.0*deltaTime),
                       (1.0 *(goodStars -oldGoodStars ))/(1.0*deltaTime),
                       (1.0 *(binCount - oldBinCount  ))/(1.0*deltaTime),
                       (1.0 *(totalMagnitudes - oldTotalMagnitudes))/(1.0*deltaTime));
                fprintf(outHandle," queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec: %f magnitudes/sec %f\n",
                        (1.0 *(GetQueryCount()-oldQueryCount))/(1.0*deltaTime),
                        (1.0 *(totalStars-oldTotalStars))/(1.0*deltaTime),
                        (1.0 *(goodStars -oldGoodStars ))/(1.0*deltaTime),
                        (1.0 *(binCount - oldBinCount  ))/(1.0*deltaTime),
                        (1.0 *(totalMagnitudes - oldTotalMagnitudes))/(1.0*deltaTime));
                oldQueryCount = GetQueryCount();
                oldTotalStars  = totalStars;
                oldGoodStars   = goodStars;
                oldBinCount    = binCount;
                oldTotalMagnitudes = totalMagnitudes;
              }
              ChargeTime(0,&curUsec,deltaTable,totalTable);
              printf("Intermediate Time (sec) ");
              fprintf(outHandle,"Intermediate Time (sec) ");
              for (usecIndex = 0; usecIndex < MAX_CHARGE_ENTRY; usecIndex++) {
                double outSeconds;
                outSeconds = (1.0 * deltaTable[usecIndex])/1000000.0;
                printf("%10.3f ",outSeconds);
                fprintf(outHandle,"%10.3f ",outSeconds);
              }
              printf("\n");
              fprintf(outHandle,"\n");

              ChargeTime(RESET_DELTA_ENTRY,&curUsec,deltaTable,totalTable);


              fflush(outHandle);

            }


          }
          if (debugMode) {
            time(&curTime);
            curTime -= startTime;
            fprintf(debug_handle,"Finished processing %s summaryCount %d starIndex %d starIndex2 %d, processedCount %d at %d\n",
                    pCurStarEntry2->REF,summaryCount,starIndex,starIndex2,processedCount,curTime);
          }
              

        }
      }
    } 
    if (errorFlag) {
      break;
    }

  





    stars_handle = fopen(stars_name,"wt");
    if (stars_handle == NULL) {
      fprintf(stderr,"ERROR: Failed to open stars file  %s\n",stars_name);
      exit(-1);
    }
    for (starIndex = 0; starIndex < curStars; starIndex++) {
      pCurStarEntry = &pStarTable[starIndex];
      fprintf(stars_handle,"%lld\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%d\t%s\n",
              pCurStarEntry->REFNumber,
              pCurStarEntry->versionId,
              pCurStarEntry->Stdmag,
              pCurStarEntry->color,
              pCurStarEntry->ra,
              pCurStarEntry->dec,
              pCurStarEntry->MAGFlag,
              pCurStarEntry->class,
              pCurStarEntry->VFlag,
              pCurStarEntry->RaPM,
              pCurStarEntry->DecPM,
              pCurStarEntry->gsc_bin_index,
              "no");
    }
    fclose(stars_handle);
    if (debugMode) {
      fprintf(outHandle,"LOAD DATA LOCAL INFILE '%s' REPLACE INTO TABLE stars%s;\n",stars_name,catalogString);
      fflush(outHandle);
    }
    sprintf(queryString,"LOAD DATA LOCAL INFILE '%s' REPLACE INTO TABLE stars%s;",stars_name,catalogString);
    ChargeTime(0,&curUsec,deltaTable,totalTable);
    res = ExecuteQuery(pPhotConnection,queryString);
    ChargeTime(4,&curUsec,deltaTable,totalTable);

    if (debugMode) {
      fprintf(outHandle,"LOAD DATA LOCAL into stars%s completed\n",catalogString);
      fflush(outHandle);
    }
    if (res) {
      errorFlag = 1;
      break;
    }
#if SKIP_REMOVAL
    fprintf(stderr,"ERROR: removal of %s skipped\n",stars_name);
#else /* SKIP_REMOVAL */
    remove(stars_name);
#endif /* SKIP_REMOVAL */
  }

  if ((errorFlag == 0) && (goodStars > 0)) {
    fclose(summary_handle);
    if (debugMode) {
      fprintf(outHandle,"LOAD DATA LOCAL INFILE '%s' REPLACE INTO TABLE starsummary%s;",summary_name,catalogString);
      fflush(outHandle);
    }
    sprintf(queryString,"LOAD DATA LOCAL INFILE '%s' REPLACE INTO TABLE starsummary%s;",summary_name,catalogString);
    ChargeTime(0,&curUsec,deltaTable,totalTable);
    res = ExecuteQuery(pPhotConnection,queryString);
    ChargeTime(3,&curUsec,deltaTable,totalTable);

    if (debugMode) {
      fprintf(outHandle,"LOAD DATA LOCAL into starsummary%s completed\n",catalogString);
      fflush(outHandle);
    }
    if (!res) {
#if SKIP_REMOVAL
      fprintf(stderr,"ERROR: removal of %s skipped\n",summary_name);
#else /* SKIP_REMOVAL */
      remove(summary_name);
#endif /* SKIP_REMOVAL */

    }
  }

  


  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins\n",curTime,GetQueryCount(),totalStars,goodStars,magnitudeAlloc,pVectorStore->vectorAlloc,binCount);

  fprintf(outHandle,"Execution Time: %d seconds, %d queries, %d totalStars, %d goodstars, %d magnitudeAlloc %d vectorAlloc, %d gscBins %lld totalMagnitudes\n",curTime,GetQueryCount(),totalStars,goodStars,magnitudeAlloc,pVectorStore->vectorAlloc,binCount,totalMagnitudes);
  printf(" queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec %f magnitudes/sec %f\n",
         (1.0 * GetQueryCount())/(1.0*curTime),
         (1.0 * totalStars)    /(1.0*curTime),
         (1.0 * goodStars )    /(1.0*curTime),
         (1.0 * binCount  )    /(1.0*curTime),
         (1.0 * totalMagnitudes)/(1.0*curTime));
  fprintf(outHandle," queries/sec: %f totalStars/sec: %f goodstars.sec: %f gscBins/sec %f magnitudes/sec %f\n",
          (1.0 * GetQueryCount())/(1.0*curTime),
          (1.0 * totalStars)    /(1.0*curTime),
          (1.0 * goodStars )    /(1.0*curTime),
          (1.0 * binCount  )    /(1.0*curTime),
          (1.0 * totalMagnitudes)/(1.0*curTime));
 
  ChargeTime(0,&curUsec,deltaTable,totalTable);
  printf("Intermediate Time (sec) ");
  fprintf(outHandle,"Intermediate Time (sec) ");
  for (usecIndex = 0; usecIndex < MAX_CHARGE_ENTRY; usecIndex++) {
    double outSeconds;
    outSeconds = (1.0 * totalTable[usecIndex])/1000000.0;
    printf("%10.3f ",outSeconds);
    fprintf(outHandle,"%10.3f ",outSeconds);
  }
  printf("\n");
  fprintf(outHandle,"\n");


  printf("... %d noGscEntry\n",noGscEntryCount);
  printf("... %d staleVersionId\n",staleVersionIdCount);
  fprintf(outHandle,"... %d noGscEntry\n",noGscEntryCount);
  fprintf(outHandle,"... %d staleVersionId\n",staleVersionIdCount);



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


  Close(catalogHandle);
  Close(indexHandle);
  fclose(outHandle);
  return(0);
}

