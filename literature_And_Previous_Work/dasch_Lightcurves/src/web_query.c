// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 * Usage: web_query designation [designation2 ...]
 *          -n <number of points>
 *          -d <temporary subdirectory name>
 *          -b <full temporary directory name>
 *          -e <executable binary image directory>
 *          -l <gzip binary image directory>
 *          -v <verbose>
 *          -q <catalog and file name qualifier>
 *          -r <search radius in arcsec>
 *          -f <authentication region flag>
 */

#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>

#include <table.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/fitshead.h>
#include <libwcs/wcs.h>
#include <libwcs/wcscat.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include "photometryutils.h"
#include "daschunistd.h"

#define MAX_BUFFER 1024
#define MAX_PATH   512
#define URL_SIZE 256
#define GSC_ALLOC_INCREMENT 1000

extern char *catalogText[MAX_CATALOG_NUMBER];
extern GSCBIN gscBin64;

extern char *urlencode(char *string, char *outstring, int outlength);
extern int parse_coo(char *coi, double *a, double *d);

typedef struct _search_center {
  char *queryName;
  char *queryNameurlencode;
  long long REFNumber;
  double ra;
  double dec;
  int haveCoordinates; /* Have valid coordinates */
  int catalogError;    /* Designation not valid for this catalog */
  int notFound;        /* Object not found */
  int unrecognized;    /* Unrecognized designation */
  int authorized;      /* Authorized access if nonzero */
} SEARCHCENTER,*PSEARCHCENTER;

/* Stars for which we have lightcurves */
typedef struct _objtable
{
  long long REFNumber;
  double ra; /* Right Ascension */
  double dec; /* Declination */
  double raGood; /* Right Ascension of good points */
  double decGood; /* Declination of good points */
  double magcal_magdep; /* Magnitude estimate */
  double magcal_magdepGood; /* Magnitude estimate of good points */
  double searchDistance; /* Search distance in arcsec */
  int searchIndex; /* Index of closest object */
  int numPoints; /* Number of points found */
  int numPointsGood; /* Number of good points found */
  int numPointsPlot; /* Number of default points to plot*/
  int catalogFlag; /* If nonzero, this is a matched object */
  int gsc_bin_index; /* GSC2.3.2 spatial bin */
  char REF[MAX_REF];          /* GSC2.3.2 reference number */
  NEARESTCATALOGSTAR nearestCatalogStar;
} OBJENTRY, *POBJENTRY;

static PGSCBIN pGscBin64 = &gscBin64;


static int
GscBinCompare(const void *first, const void *second)
{
  int coverage_bin_indexFirst = *((int *)first);
  int coverage_bin_indexSecond = *((int *)second);

  if (coverage_bin_indexFirst > coverage_bin_indexSecond) {
    return 1;
  } else if (coverage_bin_indexFirst < coverage_bin_indexSecond) {
    return -1;
  } else {
    return 0;
  }
}

static int
ObjTableCompare(const void *first, const void *second)
{
  double searchDistanceFirst = ((POBJENTRY)first)->searchDistance;
  double searchDistanceSecond = ((POBJENTRY)second)->searchDistance;

  if (searchDistanceFirst > searchDistanceSecond) {
    return 1;
  } else if (searchDistanceFirst < searchDistanceSecond) {
    return -1;
  } else {
    return 0;
  }
}

/* Given a center and a radius, get a list of prospective GSC2.3.2 bins to search.  rad0 is in arcsec */
static int
GetBinList(
  int **pgsc_bin_list,
  int *pgsc_bin_alloc,
  int *pgsc_bin_count,
  double ra0,
  double dec0,
  double rad0
) {
  int *gsc_bin_list = *pgsc_bin_list;
  int gsc_bin_alloc = *pgsc_bin_alloc;
  int gsc_bin_count = *pgsc_bin_count;
  int *tmp_gsc_bin_list = NULL;
  double rad1;

  int minDecBin;
  int maxDecBin;
  int raBin;
  int raBin0;
  int decBin;

  double dec;
  double decHigh;
  double decLow;
  double raHigh;
  double raLow;
  PBININDEX pBinIndex;
  int wrapFlag;
  int binIndex;
  int binIndex2;

  if (gsc_bin_alloc == 0) {
    gsc_bin_alloc += GSC_ALLOC_INCREMENT;
    gsc_bin_list = (int *)calloc(gsc_bin_alloc,sizeof(int));
    if (gsc_bin_list == NULL) {
      printf("ERROR: Failed to allocate gsc_bin_list for %d items\n",gsc_bin_alloc);
      exit(-1);
    }
  }

  rad1 = 3600.0*(rad0 + (2*pGscBin64->bin_size));
  GetGSCBin(pGscBin64,ra0,dec0,&minDecBin,&raBin,"ConeSearch");

  /* Find the lower search limit */
  dec = dec0-(rad1/3600.);
  if (dec < -90.0) {
    dec = -90.0;
  }
  GetGSCBin(pGscBin64,ra0,dec,&minDecBin,&raBin,"ConeSearch");

  /* Now find the uppper limit */
  dec = dec0+(rad1/3600.);
  if (dec >= +90.0) {
    dec = +90.0-(0.1*pGscBin64->bin_size);
  }
  GetGSCBin(pGscBin64,ra0,dec,&maxDecBin,&raBin,"ConeSearch");

  /* Search from South to North */
  for (decBin = minDecBin; decBin <= maxDecBin; decBin++) {
    pBinIndex = &pGscBin64->pBinMasterIndex[decBin];
    decLow = (decBin * pGscBin64->bin_size) - 90.0;
    decHigh = ((decBin+1) * pGscBin64->bin_size) - 90.0;
    raBin0 = ((ra0 * pBinIndex->numBins) / 360.0);
    /* Search to the West */
    wrapFlag = 0;
    raBin = raBin0;
    while(1) {
      if (raBin >= pBinIndex->numBins) {
        raBin = 0;
      }
      raLow = (raBin * 360.0)/pBinIndex->numBins;
      raHigh = ((raBin+1) * 360.0)/pBinIndex->numBins;
      if (((wcsdist(ra0,dec0,raLow ,decLow )*3600) < rad1) ||
          ((wcsdist(ra0,dec0,raHigh,decLow )*3600) < rad1) ||
          ((wcsdist(ra0,dec0,raLow ,decHigh)*3600) < rad1) ||
          ((wcsdist(ra0,dec0,raHigh,decHigh)*3600) < rad1)) {
        /* Add this one to our table */
        if (gsc_bin_count+1 >= gsc_bin_alloc) {
          gsc_bin_alloc += GSC_ALLOC_INCREMENT;
          tmp_gsc_bin_list = (int *)realloc(gsc_bin_list,gsc_bin_alloc * sizeof(int));
          if (tmp_gsc_bin_list == NULL) {
            printf("ERROR: failed to reallocate gsc_bin_list of size %d\n",gsc_bin_alloc);
            exit(-1);
          }
          gsc_bin_list = tmp_gsc_bin_list;
          tmp_gsc_bin_list = NULL;

        }
        gsc_bin_list[gsc_bin_count] = pBinIndex->startBin+raBin;
        gsc_bin_count++;
      } else {
        break;
      }

      raBin++;
      if (raBin == raBin0) {
        /* We wrapped! */
        wrapFlag = 1;
        break;
      }
    }

    /* Search to the East */
    if (wrapFlag == 0) {
      raBin = raBin0 -1;
      while(1) {
        if (raBin < 0) {
          raBin = pBinIndex->numBins -1;
        }
        raLow = (raBin * 360.0)/pBinIndex->numBins;
        raHigh = ((raBin+1) * 360.0)/pBinIndex->numBins;
        if (((wcsdist(ra0,dec0,raLow ,decLow )*3600) < rad1) ||
            ((wcsdist(ra0,dec0,raHigh,decLow )*3600) < rad1) ||
            ((wcsdist(ra0,dec0,raLow ,decHigh)*3600) < rad1) ||
            ((wcsdist(ra0,dec0,raHigh,decHigh)*3600) < rad1)) {
          /* Add this one to our table */
          if (gsc_bin_count+1 >= gsc_bin_alloc) {
            gsc_bin_alloc += GSC_ALLOC_INCREMENT;
            tmp_gsc_bin_list = (int *)realloc(gsc_bin_list,gsc_bin_alloc * sizeof(int));
            if (tmp_gsc_bin_list == NULL) {
              printf("ERROR: failed to reallocate gsc_bin_list of size %d\n",gsc_bin_alloc);
              exit(-1);
            }
            gsc_bin_list = tmp_gsc_bin_list;
            tmp_gsc_bin_list = NULL;

          }
          gsc_bin_list[gsc_bin_count] = pBinIndex->startBin+raBin;
          gsc_bin_count++;
        } else {
          break;
        }
        raBin--;
      }
    }
  }

  /* For performance, sort in increasing bin order */
  qsort((void*)gsc_bin_list,gsc_bin_count,sizeof(int),GscBinCompare);
  for (binIndex = 1; binIndex < gsc_bin_count; binIndex++) {
    if (gsc_bin_list[binIndex-1] == gsc_bin_list[binIndex]) {
      /* This is a duplicate, shift everything down */
      for (binIndex2 = binIndex; binIndex2 < (gsc_bin_count-1); binIndex2++) {
        gsc_bin_list[binIndex2] = gsc_bin_list[binIndex2+1];
      }
      gsc_bin_count--;

    }
  }

  gsc_bin_list = *pgsc_bin_list = gsc_bin_list;
  gsc_bin_alloc = *pgsc_bin_alloc = gsc_bin_alloc;
  gsc_bin_count = *pgsc_bin_count = gsc_bin_count;;
  return 0;
}


int
main(int argc, char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;

  int verbose = 0;
  int verbose2 = 0;

  int lineLen;

  PPHOTTARGET target_table = NULL;
  size_t target_nrecs = 0;
  size_t target_alloc = 0;
  int target_index;
  PPHOTTARGET pTarget;

  int decBin;
  int raBin;
  int decBin2;
  int raBin2;

  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;

  int curStars;
  STARENTRY starTable;
  PSTARENTRY pStarTable = &starTable;

  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  int gotAnswer;
  char *tmpdir = NULL;
  char *fulltmpdir = NULL;

  PFILESTARIMAGE pFileStarImage = NULL;
  int magnitudeIndex;
  PPHOTSTARIMAGE pMagnitudeTable = NULL;
  PFILESTARIMAGE pFileStarImageTable = NULL;
  PPHOTSTARIMAGE pMagnitudeTable1 = NULL;
  PPHOTSTARIMAGE pMagnitudeTable2 = NULL;
  int allocMagnitudes = 0;
  int allocMagnitudes1 = 0;
  int allocMagnitudes2 = 0;

  PPHOTSTARIMAGE pCurStarImage;
  PFILESTARIMAGE pNoneMagnitudeTable = NULL;
  int noneMagnitudeAlloc = 0;
  int numMagnitudes;
  int numMagnitudes1;
  int numMagnitudes2;

  double rad0 = 0; /* Search radius in arcsec */
  int statResult;
  struct stat statbuf;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  char catalogString1[MAX_BUFFER];
  char catalogString2[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  char source[MAX_BUFFER];
  FILECOMMON fileCommon0;
  PFILECOMMON pFileCommon0 = &fileCommon0;
  FILECOMMON fileCommon1;
  PFILECOMMON pFileCommon1 = &fileCommon1;
  FILECOMMON fileCommon2;
  PFILECOMMON pFileCommon2 = &fileCommon2;
  int histogramTable[MAX_NONE_HISTOGRAM];
  int groupNumber = 1;
  int numPoints = 0;

  PSEARCHCENTER search_center_table = NULL;
  PSEARCHCENTER tmp_search_center_table;
  PSEARCHCENTER pSearchCenter;
  int searchCenterSize = 0;
  int searchCenterAlloc = 100;

  POBJENTRY object_table = NULL;
  POBJENTRY tmp_object_table;
  POBJENTRY pObjectEntry;
  int objTableAlloc = 100;
  int objTableSize = 0;
  int object_index;

  int index;
  char* pDesignation;
  int refType;
  int searchStarsTable;
  char rstr[32], dstr[32];
  int lobj;
  char url[URL_SIZE];
  int lbuff;
  char *buff;
  int i;
  int *gsc_bin_list = NULL;
  int gsc_bin_alloc = 0;
  int gsc_bin_count = 0;
  int bin_list_index;
  int gsc_bin_index = 0;
  int printFlag;
  int printedObjects = 0;
  char REFurlencode[4*MAX_REF];          /* GSC2.3.2 reference number */
  char cmdStr[MAX_BUFFER];
  char table_file_name[MAX_PATH];
  char vo_file_name[MAX_PATH];
  char vo_file_name_gz[MAX_PATH];
  FILE* table_file_handle = NULL;
  char object_file_name[MAX_PATH];
  char object_xml_name[MAX_PATH];
  FILE* object_file_handle = NULL;
  int result;
  char binaries[MAX_BUFFER];
  char lsbin[MAX_BUFFER];
  char* regionFlag = NULL;
  double tmpra;
  double tmpdec;
  GALAXYCOMMON galaxycommon;
  PGALAXYCOMMON pGalaxyCommon = &galaxycommon;
  int haveGalaxyTable = 1;
  char nearbyObjects[MAX_NEARBY_OBJECTS_STRING];
  char nearbyObjects2[MAX_NEARBY_OBJECTS_STRING + 8];
  int releaseField;
  PHOTPLATES photPlates;
  PPHOTPLATES pPhotPlates = &photPlates;
  int writeHeader = 1;
  char *pTabChar;
  int combineAlgorithm = COMBINE_ALGORITHM_DEFAULT;
  char rematchString[10];
  int enableRematch = 0;
  char starbase_title[MAX_BUFFER];

  PSUBSTITUTION substitutionTable = NULL;

  memset(pGalaxyCommon,0,sizeof(GALAXYCOMMON));
  memset(pPhotPlates,0,sizeof(PHOTPLATES));
  pPhotPlates->mosaicNumber = 99;
  pPhotPlates->quality = QUALITY_UNINITIALIZED;

  table_file_name[0] = 0;
  object_file_name[0] = 0;

  search_center_table = (PSEARCHCENTER)calloc(searchCenterAlloc,sizeof(SEARCHCENTER));
  if (search_center_table == NULL) {
    printf("ERROR: failed to allocate search_center_table of size %d\n",searchCenterAlloc);
    exit(-1);
  }

  object_table = (POBJENTRY)calloc(objTableAlloc,sizeof(OBJENTRY));
  if (object_table == NULL) {
    printf("ERROR: failed to allocate object_table of size %d\n",objTableAlloc);
    exit(-1);
  }

  catalogString[0] = 0;
  catalogString1[0] = 0;
  catalogString2[0] = 0;
  SetQueryCount(0);
  qualifier[0] = 0;
  source[0] = 0;
  binaries[0] = 0;
  lsbin[0] = 0;

  /* Loop through the arguments. If there's a usage error, report it but plow on
   * ahead; since we're a web interface script, we want to try as hard as
   * possible to succeed, but we also don't want to leak information about this
   * program. */

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be a designation */
      pDesignation = strtok(argstr,"\n");
      while (pDesignation != NULL) {
        if ((RELEASE_LEVEL <= RELEASE_LEVEL_DR1) &&
            (searchCenterSize >= RELEASE_FIELD_SEARCH_LIMIT)) {
          break;
        }

        if ((searchCenterSize+1) >= searchCenterAlloc) {
          searchCenterAlloc += 100;
          tmp_search_center_table = (PSEARCHCENTER)realloc(search_center_table,searchCenterAlloc*sizeof(SEARCHCENTER));
          if (tmp_search_center_table == NULL) {
            printf("ERROR: failed to realloc search_center_table of size %d\n",searchCenterAlloc);
            exit(-1);
          }
          search_center_table = tmp_search_center_table;
          tmp_search_center_table = NULL;
        }

        pSearchCenter = &search_center_table[searchCenterSize];
        memset(pSearchCenter,0,sizeof(SEARCHCENTER));
        pSearchCenter->authorized = 1;
        pSearchCenter->queryName = (char *)calloc(strlen(pDesignation)+2,sizeof(char));
        pSearchCenter->queryNameurlencode = (char *)calloc(4*(strlen(pDesignation)+2),sizeof(char));
        strcpy(pSearchCenter->queryName,pDesignation);
        lineLen = strlen(pSearchCenter->queryName);
        /* Trim off the carriage return */
        if (pSearchCenter->queryName[lineLen-1] == 10) {
          pSearchCenter->queryName[lineLen-1] = 0;
        }
        /* Trim off the line feed */
        lineLen = strlen(pSearchCenter->queryName);
        if (pSearchCenter->queryName[lineLen-1] == 13) {
          pSearchCenter->queryName[lineLen-1] = 0;
        }
        /* Now convert tabs to spaces */
        while (1) {
          pTabChar = strstr(pSearchCenter->queryName,"\t");
          if (pTabChar == NULL) {
            break;
          }
          *pTabChar = ' ';
        }

        urlencode(pSearchCenter->queryName,pSearchCenter->queryNameurlencode,4*(strlen(pDesignation)));

        searchCenterSize++;
        pDesignation = strtok(NULL,"\n");
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'n': /* number of points */
        case 'N':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            nvals = sscanf(*++argv,"%d",&numPoints);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the numPoints %s\n",*argv);
            }
          }
          break;


        case 'd': /* temporary subdirectory name */
        case 'D':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            tmpdir = *++argv;
          }
          break;

        case 'b': /* full temporary directory name */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            fulltmpdir = *++argv;
          }
          break;

        case 'e': /* Executable binary image directory */
        case 'E':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: binaries directory length %zu for %s is too long\n",strlen(*argv),*argv);
            } else {
              strcpy(binaries,*argv);
            }
          }
          break;

        case 'l': /* gzip binary image directory */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: lsbin directory length %zu for %s is too long\n",strlen(*argv),*argv);
            } else {
              strcpy(lsbin,*argv);
            }
          }
          break;

        case 'v': /* verbose */
          verbose += 1;
          verbose2 = -1; /* For LocateNoneImages only */
          break;

        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              fprintf(stderr,"ERROR: Illegal catalog name %s\n",*argv);
            } else {
              strcpy(source,*argv);
              if (catalogNumber > 0) {
                sprintf(catalogString,"%d",catalogNumber);
                strcpy(qualifier,*argv);
              }
            }
          }
          break;

        case 'r': /* Use a search radius in arcsec */
        case 'R': /* Use a search radius in arcsec */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            ++argv;
            if (strchr(*argv,':')) {
              rad0 = 3600*str2dec(*argv);
            } else {
              rad0 = atof(*argv);
            }
            if (rad0 == 0.0) {
              printf("ERROR: Illegal radius %s specified\n",*argv);
            }
          }
          break;

        case 'f': /* Region Flag */
        case 'F':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            regionFlag = *++argv;
          }
          break;

        case 'O':
          enableRematch = 1;
          break;

        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
        }
      }
    }
  }

  if (fulltmpdir != NULL) {
    sprintf(table_file_name,"%s/table_%s.db",fulltmpdir,tmpdir);
    sprintf(vo_file_name,"%s/table_%s.xml",fulltmpdir,tmpdir);
    strcpy(vo_file_name_gz,vo_file_name);
    strcat(vo_file_name_gz,".gz");
    sprintf(object_file_name,"%s/object_%s.db",fulltmpdir,tmpdir);
    sprintf(object_xml_name,"%s/object_%s.xml.gz",fulltmpdir,tmpdir);
  }

  /* Now display our banner */
  if (enableRematch == 0) {
    printf("<H3>DASCH (%s) Catalog Query Results (%.0f\")</H3>",catalogText[catalogNumber],rad0);
    rematchString[0] = 0;
  } else {
    printf("<H3>DASCH (%s) Transient Optimization Results (%.0f\")</H3>",catalogText[catalogNumber],rad0);
    strcpy(rematchString,"checked");
  }
  printf("<pre>");

  if (strlen(binaries) == 0) {
    printf("ERROR: binaries directory is not specified\n");
  }

  if (strlen(lsbin) == 0) {
    printf("ERROR: lsbin directory is not specified\n");
  }

  if (regionFlag == NULL) {
    printf("ERROR: regionFlag is not specified\n");
    return -1;
  }

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pConnection);
  dasch_init_photdb(pPhotConnection);

  if (OpenGalaxyFiles(pGalaxyCommon, qualifier, NULL, GetPhotFileBase(catalogString)) != 0) {
    haveGalaxyTable = 0;
  }

  InitSeriesTable(pConnection,pPhotConnection);
  gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
  if (gotAnswer != 1)  {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }

  if (pPhotGlobal->magnitudeFile != PHOT_MAGNITUDEFILE_YES) {
    printf("ERROR: database uses obsolute magnitude table\n");
    exit(-1);
  }

  if ((RELEASE_EXPERIMENTAL != 0) && (catalogNumber == CATALOG_EXPERIMENTAL)) {
    InitFileCommon(stdout,pFileCommon0,pConnection,pPhotConnection,catalogNumber,0);
    InitMaxPlateNumber(pConnection,pFileCommon0->maxPlateNumber);
    InitFileCommon(stdout,pFileCommon1,pConnection,pPhotConnection,RELEASE_EXPERIMENTAL_CAT0,0);
    InitMaxPlateNumber(pConnection,pFileCommon1->maxPlateNumber);
    if (RELEASE_EXPERIMENTAL_CAT0 != 0) {
      sprintf(catalogString1,"%d",RELEASE_EXPERIMENTAL_CAT0);
    }
    InitFileCommon(stdout,pFileCommon2,pConnection,pPhotConnection,RELEASE_EXPERIMENTAL_CAT1,0);
    InitMaxPlateNumber(pConnection,pFileCommon2->maxPlateNumber);
    if (RELEASE_EXPERIMENTAL_CAT1 != 0) {
      sprintf(catalogString2,"%d",RELEASE_EXPERIMENTAL_CAT1);
    }
  } else {
    InitFileCommon(stdout,pFileCommon0,pConnection,pPhotConnection,catalogNumber,0);
    InitMaxPlateNumber(pConnection,pFileCommon0->maxPlateNumber);
  }

  if (enableRematch) {
    pFileCommon0->enableRematch = 1;
    pFileCommon1->enableRematch = 1;
    pFileCommon2->enableRematch = 1;
  }

  /* Now get the center RA and DEC for all of these strings */
  for (index = 0; index < searchCenterSize; index++) {
    pSearchCenter = &search_center_table[index];
    searchStarsTable = 0;

    if(parse_coo(pSearchCenter->queryName,&pSearchCenter->ra,&pSearchCenter->dec)== 0) {
      pSearchCenter->haveCoordinates = 1;
    } else if (GetREFNumber(pSearchCenter->queryName,&pSearchCenter->REFNumber,&refType,0,0) == 0) {
      switch(refType) {
      case REF_TYPE_NONE:
        pSearchCenter->unrecognized = 1;
        break;
      case REF_TYPE_GSC:
      case REF_TYPE_TYCHO2:
        if ((catalogNumber != CATALOG_GSC232) &&
            (catalogNumber != CATALOG_APASS)){
          pSearchCenter->catalogError = 1;
        } else {
          searchStarsTable = 1;
        }
        break;
      case REF_TYPE_KEPLER:
        if (catalogNumber != CATALOG_KEPLER) {
          pSearchCenter->catalogError = 1;
        } else {
          searchStarsTable = 1;
        }
        break;
      case REF_TYPE_ATLAS2:
        if (catalogNumber != CATALOG_ATLAS) {
          pSearchCenter->catalogError = 1;
        } else {
          searchStarsTable = 1;
        }
        break;
      case REF_TYPE_GAIA2:
        if (catalogNumber != CATALOG_GAIA) {
          pSearchCenter->catalogError = 1;
        } else {
          searchStarsTable = 1;
        }
        break;
      case REF_TYPE_DASCH:
        if (GetDASCHCoordinates(pSearchCenter->queryName,
                                &pSearchCenter->ra,
                                &pSearchCenter->dec,0,REF_TYPE_DASCH) == 0) {
          pSearchCenter->haveCoordinates = 1;
        } else {
          pSearchCenter->unrecognized = 1;
        }

        break;
      case REF_TYPE_APASS:
        if (GetDASCHCoordinates(pSearchCenter->queryName,
                                &pSearchCenter->ra,
                                &pSearchCenter->dec,0,REF_TYPE_APASS) == 0) {
          pSearchCenter->haveCoordinates = 1;
        } else {
          pSearchCenter->unrecognized = 1;
        }

        break;
      default:
        break;
      }

    } else {

      /* Query the "identifier" database */
      if (GetStarEntry3(pPhotConnection,pSearchCenter->queryName,&pSearchCenter->ra,&pSearchCenter->dec) == 1) {
        pSearchCenter->haveCoordinates = 1;
      } else {
        char *posline;
        char *posra;
        char *posdec;
        /* Ask Simbad about this one, using Doug Mink's code wcstools (simpos.c) */
        /* Replace underscores and spaces with plusses */
        lobj = strlen (pSearchCenter->queryName);
        for (i = 0; i < lobj; i++) {
          if (pSearchCenter->queryName[i] == '_') {
            pSearchCenter->queryName[i] = '+';
          }
          if (pSearchCenter->queryName[i] == ' ') {
            pSearchCenter->queryName[i] = '+';
          }
        }

        strcpy (url,"http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oI/SNV?"); /* Change made Feb 4, 2019 */
        strncat (url,pSearchCenter->queryNameurlencode,URL_SIZE-strlen(url)-2);
        url[URL_SIZE-1] = 0;
        buff = webbuff (url, verbose, &lbuff);
        if (buff == NULL) {
          pSearchCenter->unrecognized = 1;

        } else if ((posline = strsrch (buff, "%J ")) != NULL) {
          posra = posline + 3;
          while (*posra == ' ')
            posra++;
          pSearchCenter->ra = atof (posra);
          posdec = strchr (posra, ' ');
          while (*posdec == ' ')
            posdec++;
          posline = strchr (posdec, ' ');
          pSearchCenter->dec = atof (posdec);
          pSearchCenter->haveCoordinates = 1;
        } else {
          pSearchCenter->unrecognized = 1;

        }

        free(buff);

        for (i = 0; i < lobj; i++) {
          if (pSearchCenter->queryName[i] == '+') {
            pSearchCenter->queryName[i] = ' ';
          }
        }
      }
    }

    if (searchStarsTable) {
      if ((curStars = GetStarEntry2(pPhotConnection,pStarTable,pSearchCenter->queryName,0,1)) > 0) {
        pSearchCenter->ra = pStarTable->ra;
        pSearchCenter->dec = pStarTable->dec;
        pSearchCenter->haveCoordinates = 1;
      } else {
        pSearchCenter->notFound = 1;
      }
    }

    if (pSearchCenter->haveCoordinates) {

      if (CheckAuthorization(regionFlag,pSearchCenter->ra,pSearchCenter->dec,&releaseField)) {
        GetBinList(&gsc_bin_list,&gsc_bin_alloc,&gsc_bin_count,pSearchCenter->ra,pSearchCenter->dec,(rad0/3600)+BIN_EXPANSION_RADIUS);
      } else {
        pSearchCenter->authorized = 0;
      }
    }
  }

  /* Now search for objects */

  bin_list_index = 0;
  while (bin_list_index < gsc_bin_count) {
    gsc_bin_index = gsc_bin_list[bin_list_index];
    bin_list_index++;
    /* Find all of the images in this bin */
    if ((RELEASE_EXPERIMENTAL != 0) && (catalogNumber == CATALOG_EXPERIMENTAL)) {
      LocateNoneImages(pGscBin64,pFileCommon1,gsc_bin_index,&pMagnitudeTable1,&allocMagnitudes1,&numMagnitudes1,catalogString1,1,verbose2);
      LocateNoneImages(pGscBin64,pFileCommon2,gsc_bin_index,&pMagnitudeTable2,&allocMagnitudes2,&numMagnitudes2,catalogString2,1,verbose2);
      CombineMultipleTables(pPhotGlobal,&pMagnitudeTable,&pFileStarImageTable,pMagnitudeTable1,pMagnitudeTable2,&substitutionTable,numMagnitudes1,numMagnitudes2,RELEASE_EXPERIMENTAL_CAT0,RELEASE_EXPERIMENTAL_CAT1,&allocMagnitudes,&numMagnitudes,combineAlgorithm,gsc_bin_index);

    } else {
      LocateNoneImages(pGscBin64,pFileCommon0,gsc_bin_index,&pMagnitudeTable,&allocMagnitudes,&numMagnitudes,catalogString,1,verbose2);
    }

    if (numMagnitudes > 0) {
      if (numMagnitudes >= noneMagnitudeAlloc) {
        noneMagnitudeAlloc = numMagnitudes + 1000;
        if (pNoneMagnitudeTable != NULL) {
          free(pNoneMagnitudeTable);
        }
        pNoneMagnitudeTable = (PFILESTARIMAGE)calloc(noneMagnitudeAlloc,sizeof(FILESTARIMAGE));
        if (pNoneMagnitudeTable == NULL) {
          printf("ERROR: failed to allocate pNoneMagnitudeTable of size %d\n",noneMagnitudeAlloc);
          exit(-1);
        }
      }
      for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
        pCurStarImage = &pMagnitudeTable[magnitudeIndex];
        pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
        memcpy(pFileStarImage,pCurStarImage->pFileStarImage,sizeof(FILESTARIMAGE));
      }

      target_nrecs = 0;

      ProcessNoneImagesX(
        pGscBin64,
        pFileCommon0,
        &target_table,
        &target_nrecs,
        &target_alloc,
        pNoneMagnitudeTable,
        numMagnitudes,
        verbose,
        histogramTable,
        &groupNumber,
        gsc_bin_index
      );

      for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
        pCurStarImage = &pMagnitudeTable[magnitudeIndex];
        pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];

        if ((pFileStarImage->gsc_bin_index == gsc_bin_index) &&
            (pFileStarImage->magcal_magdep < 90.0) &&
            (pFileStarImage->REFNumber != 0)) {
          if ((object_file_handle == NULL) && (object_file_name[0] != 0)) {
            object_file_handle = fopen(object_file_name,"wt");
          }
          if (object_file_handle != NULL) {
            pPhotPlates->versionId = pFileStarImage->versionId;
            if (enableRematch == 0) {
              strcpy(starbase_title,"Catalog Query Data");
            } else {
              strcpy(starbase_title,"Catalog Query Data Optimized for Transients");
            }
            WriteStarbaseRecord(pFileStarImage,pPhotPlates,object_file_handle,starbase_title,&writeHeader,catalogNumber);
          }

          for (object_index = 0; object_index < objTableSize; object_index++) {
            pObjectEntry = &object_table[object_index];
            if (pFileStarImage->REFNumber == pObjectEntry->REFNumber) {
              break;
            }
          }

          if (object_index >= objTableSize) {
            if ((objTableSize+1) >= objTableAlloc) {
              /* Time to reallocate the table */
              objTableAlloc += 100;
              tmp_object_table = (POBJENTRY)realloc(object_table,objTableAlloc*sizeof(OBJENTRY));
              if (tmp_object_table == NULL) {
                printf("ERROR: failed to realloc object_table of size %d\n",objTableAlloc);
                exit(-1);
              }
              object_table = tmp_object_table;
              tmp_object_table = NULL;

            }
            pObjectEntry = &object_table[objTableSize];
            memset(pObjectEntry,0,sizeof(OBJENTRY));
            pObjectEntry->dec = 99.0;
            pObjectEntry->ra = 999.0;
            pObjectEntry->gsc_bin_index = -1;
            InitNearestCatalogStar(&pObjectEntry->nearestCatalogStar);
            pObjectEntry->REFNumber = pFileStarImage->REFNumber;
            GetREF(pObjectEntry->REFNumber,pObjectEntry->REF,1,1);

            if (pFileStarImage->DecPM < MAX_PROPER_MOTION) {
              double epoch;
              double factor;
              epoch = jd2ep(pFileStarImage->Date) - GSC_EQUINOX;
              pObjectEntry->dec =  pFileStarImage->dec_2 - ((pFileStarImage->DecPM * epoch)/(3600.0*1000.0));
              if ((pObjectEntry->dec < 90.0) && (pObjectEntry->dec > -90.0)) {
                factor = cos(pObjectEntry->dec * DEGREES_TO_RAD);
                pObjectEntry->ra = pFileStarImage->ra_2 - (((pFileStarImage->RaPM *epoch)/(3600.0*1000.0))/factor);
              } else {
                pObjectEntry->dec = pFileStarImage->ra_2;
                factor = 1.0;
              }
              pObjectEntry->raGood = pObjectEntry->ra;
              pObjectEntry->decGood = pObjectEntry->dec;
            }

            pObjectEntry->magcal_magdep = pFileStarImage->magcal_magdep;
            pObjectEntry->numPoints = 1;
            pObjectEntry->gsc_bin_index = pFileStarImage->gsc_bin_index;
            pObjectEntry->catalogFlag = 1;
            if ((pFileStarImage->AFLAGS & FILTER_AMASK_LOWDRAD) == 0) {
              pObjectEntry->magcal_magdepGood = pFileStarImage->magcal_magdep;
              pObjectEntry->numPointsGood = 1;

            }
            if ((pFileStarImage->AFLAGS & (~FILTER_AMASK_PLOT2)) == 0) {
              pObjectEntry->numPointsPlot = 1;
            }
            objTableSize++;
          } else { /* Got a new object */
            /* Got an existing object, just increment the count */
            pObjectEntry->numPoints++;
            pObjectEntry->magcal_magdep += pFileStarImage->magcal_magdep;

            if ((pFileStarImage->AFLAGS & (~FILTER_AMASK_PLOT2)) == 0) {
              pObjectEntry->numPointsPlot++;
            }

            if ((pFileStarImage->AFLAGS & FILTER_AMASK_LOWDRAD) == 0) {
              pObjectEntry->numPointsGood++;
              pObjectEntry->magcal_magdepGood += pFileStarImage->magcal_magdep;
            }
          }
        } /* Got a legal REFNumber */
      }

      if (target_nrecs > 0) {
        /* Have some new objects, transfer them to our object list */
        for (target_index = 0; target_index < target_nrecs; target_index++) {
          long long REFNumber;
          int RefType;

          pTarget = &target_table[target_index];
          if (pTarget->gsc_bin_index != gsc_bin_index) {
            continue;
          }

          if (GetREFNumber(pTarget->REF,&REFNumber,&RefType,0,0) == 0) {
            for (object_index = 0; object_index < objTableSize; object_index++) {
              pObjectEntry = &object_table[object_index];
              if (REFNumber == pObjectEntry->REFNumber) {
                break;
              }
            }

            if (object_index >= objTableSize) {
              if ((objTableSize+1) >= objTableAlloc) {
                /* Time to reallocate the table */
                objTableAlloc += 100;
                tmp_object_table = (POBJENTRY)realloc(object_table,objTableAlloc*sizeof(OBJENTRY));
                if (tmp_object_table == NULL) {
                  printf("ERROR: failed to realloc object_table of size %d\n",objTableAlloc);
                  exit(-1);
                }
                object_table = tmp_object_table;
                tmp_object_table = NULL;
              }

              pObjectEntry = &object_table[objTableSize];
              memset(pObjectEntry,0,sizeof(OBJENTRY));
              pObjectEntry->gsc_bin_index = pTarget->gsc_bin_index;
              InitNearestCatalogStar(&pObjectEntry->nearestCatalogStar);
              pObjectEntry->REFNumber = REFNumber;
              GetREF(pObjectEntry->REFNumber,pObjectEntry->REF,1,1);

              if (GetDASCHCoordinates(pObjectEntry->REF,&pObjectEntry->ra,&pObjectEntry->dec,1,REF_TYPE_DASCH) != 0) {
                printf("ERROR: Unable to get DASCH coordinates for %s\n",pObjectEntry->REF);
              }

              pObjectEntry->magcal_magdep = pTarget->magcal_magdep;
              memcpy(&pObjectEntry->nearestCatalogStar,&pTarget->nearestCatalogStar,sizeof(NEARESTCATALOGSTAR));
              for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
                pCurStarImage = &pMagnitudeTable[magnitudeIndex];
                pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
                if (pObjectEntry->REFNumber != pFileStarImage->REFNumber) {
                  continue;
                }
                pObjectEntry->numPoints++;
                if ((pFileStarImage->AFLAGS & FILTER_AMASK_LOWDRAD) == 0) {
                  pObjectEntry->numPointsGood++;
                }
                if ((pFileStarImage->AFLAGS & (~FILTER_AMASK_PLOT2)) == 0) {
                  pObjectEntry->numPointsPlot++;
                }
              }
              objTableSize++;
            } /* Got a new object */
          } /* Got a legal REFNumber */
        } /* End of object table loop */
      } /* Got valid none objects */
    } /* Got a bin with magnitudes */
  } /* End of bin index loop */

  /* Now for each object, find the source that it is closest to */

  for (object_index = 0; object_index < objTableSize; object_index++) {
    pObjectEntry = &object_table[object_index];
    pObjectEntry->searchIndex = -1;
    pObjectEntry->searchDistance = 0;

    if ((pObjectEntry->catalogFlag != 0) && (pObjectEntry->numPoints > 0)) {
      if (pObjectEntry->numPointsGood > 0) {
        pObjectEntry->magcal_magdep = pObjectEntry->magcal_magdepGood/pObjectEntry->numPointsGood;
      } else {
        pObjectEntry->magcal_magdep = pObjectEntry->magcal_magdep/pObjectEntry->numPoints;
      }

      refType = GetREFType(pObjectEntry->REFNumber);

      switch(refType) {
      case REF_TYPE_APASS:
      case REF_TYPE_DASCH:
        if (GetDASCHCoordinates(pObjectEntry->REF,
                                &tmpra,
                                &tmpdec,0,refType) == 0) {
          pObjectEntry->ra = tmpra;
          pObjectEntry->dec = tmpdec;
          if (pObjectEntry->gsc_bin_index < 0) {
            pObjectEntry->gsc_bin_index = GetGSCBin(pGscBin64,tmpra,tmpdec,&decBin,&raBin,"ConeSearch");
          }
        }
        break;
      default:
        break;
      }
    }

    for (index = 0; index < searchCenterSize; index++) {
      pSearchCenter = &search_center_table[index];
      if ((pSearchCenter->haveCoordinates) && (pSearchCenter->authorized)) {
        double rad1;
        rad1 = wcsdist(pObjectEntry->ra,pObjectEntry->dec,pSearchCenter->ra,pSearchCenter->dec) * 3600;
        if ((pObjectEntry->searchIndex < 0) ||
            (rad1 < pObjectEntry->searchDistance)) {
          pObjectEntry->searchIndex = index;
          pObjectEntry->searchDistance = rad1;

        }
      }
    }
  }

  /* At this point, sort all of the targets with respect to search distance */
  qsort((void*)object_table,objTableSize,sizeof(OBJENTRY),ObjTableCompare);

  for (index = 0; index < searchCenterSize; index++) {
    pSearchCenter = &search_center_table[index];

    if ((pSearchCenter->haveCoordinates) && (pSearchCenter->authorized)) {
      printFlag = 1;
      printedObjects = 0;
      ra2str(rstr,16,pSearchCenter->ra,1);
      dec2str(dstr,16,pSearchCenter->dec,0);
      if (haveGalaxyTable) {
        haveGalaxyTable = LoadGalaxyTable(pGalaxyCommon,pSearchCenter->ra,pSearchCenter->dec);
      }
      printf("<B>'%s'  ra: %s  dec: %s",pSearchCenter->queryName,rstr,dstr);
      if (haveGalaxyTable) {
        printf(" approx. plates: %d\n",pGalaxyCommon->plateCount);
      }
      printf("</B>\n");

      /* Now list points to the objects */;
      for (object_index = 0; object_index < objTableSize; object_index++) {
        pObjectEntry = &object_table[object_index];

        if ((pObjectEntry->searchIndex == index) &&
            (pObjectEntry->searchDistance <= rad0) &&
            (pObjectEntry->numPoints >= numPoints)) {
          printedObjects++;
          if (printFlag) {
            printf(
              "arcsec  Nobs Nplot   mag  id "
              "(<a target=\"_blank\" href=\"https://dasch.cfa.harvard.edu/data-guide/#nearbyObjects_ext\">nearbyObjects</a>"
              " [(mag)] [:object type] @ arcsec distance)\n"
            );
            printFlag = 0;
          }

          /* Escape the plus sign in DASCH references */
          urlencode(pObjectEntry->REF,REFurlencode,4*MAX_REF);
          FindNearestGalaxy(pGalaxyCommon,pObjectEntry->ra,pObjectEntry->dec,&pObjectEntry->nearestCatalogStar,nearbyObjects,NULL);

          /* At this point, open a tabular form of the search table */

          if ((table_file_handle == NULL) && (table_file_name[0] != 0)) {
            table_file_handle = fopen(table_file_name,"wt");
            if (table_file_handle) {
              fprintf(table_file_handle,"src_name\tcra\tcdec\tREF\tdrad\tnpoints\tnplot\tmagcal_magdep\tra\tdec\tnearbyObjects\n");
              fprintf(table_file_handle,"--------\t---\t----\t---\t----\t-------\t-----\t-------------\t--\t---\t-------------\n");
            }
          }

          if (table_file_handle != NULL) {
            fprintf(table_file_handle,"%s\t%.4f\t%.4f\t%s\t%.0f\t%d\t%d\t%.2f\t%.4f\t%.4f\t%s\n",
                    pSearchCenter->queryName,
                    pSearchCenter->ra,
                    pSearchCenter->dec,
                    pObjectEntry->REF,
                    pObjectEntry->searchDistance,
                    pObjectEntry->numPoints,
                    pObjectEntry->numPointsPlot,
                    pObjectEntry->magcal_magdep,
                    pObjectEntry->ra,
                    pObjectEntry->dec,
                    nearbyObjects);
          }

          if (strlen(nearbyObjects) > 0) {
            sprintf(nearbyObjects2," (%s)",nearbyObjects);
          } else {
            nearbyObjects2[0] = 0;
          }

          printf("%5.0f  %5d %5d %5.2f  <a href=\"lightcurve_plot_size.php?REF=%s&center=%s&distance=%.0f&gsc_bin_index=%d&tmpdir=%s&source=%s&rematch=%s&AFLAGS=%d&quality=%d&imagesize=%d&sizeunits=%s\" target=\"data_plot\" >%s</a>%s\n",
                 pObjectEntry->searchDistance,
                 pObjectEntry->numPoints,
                 pObjectEntry->numPointsPlot,
                 pObjectEntry->magcal_magdep,
                 REFurlencode,
                 pSearchCenter->queryNameurlencode,
                 pObjectEntry->searchDistance,
                 pObjectEntry->gsc_bin_index,
                 tmpdir,
                 source,
                 rematchString,
                 FILTER_AMASK_PLOT2,
                 DEFAULT_QUALITY_MASK,
                 THUMBNAIL_IMAGESIZE,
                 THUMBNAIL_ARCSEC,
                 pObjectEntry->REF,
                 nearbyObjects2);
        }
      }

      if (printedObjects == 0) {
        if ((haveGalaxyTable != 0) &&
            (pGalaxyCommon->plateCount > 0) &&
            (pSearchCenter->haveCoordinates != 0)) {
          gsc_bin_index = GetGSCBin(pGscBin64,pSearchCenter->ra,pSearchCenter->dec,&decBin2,&raBin2,"web_query");

          printf("<a href=\"lightcurve_plot_size.php?REF=%s&center=%s&distance=%.0f&gsc_bin_index=%d&tmpdir=%s&source=%s&rematch=%s&AFLAGS=%d&quality=%d&imagesize=%d&sizeunits=%s\" target=\"data_plot\" >No Objects Found. Show Nearby Limiting Magnitudes</a>\n",
                 "NONE",
                 pSearchCenter->queryNameurlencode,
                 0.0,
                 gsc_bin_index,
                 tmpdir,
                 source,
                 rematchString,
                 0x7FFFFFFF,
                 0x7FFFFFFF,
                 THUMBNAIL_IMAGESIZE,
                 THUMBNAIL_ARCSEC);
        } else {
          printf("                   No Plates Scanned\n");
        }
      }
    } else if (pSearchCenter->authorized == 0) {
      printf("Sorry, these data are not available due to <a target=\"_blank\" href=\"https://dasch.cfa.harvard.edu/data-access/#restrictions\">temporary data access restrictions</a><br/>\n");
    }  else if (pSearchCenter->catalogError) {
      printf("<B>'%s' Incorrect source catalog </B>\n",pSearchCenter->queryName);
    } else if (pSearchCenter->notFound) {
      printf("<B>'%s' NOT FOUND </B>\n",pSearchCenter->queryName);
    } else if (pSearchCenter->unrecognized) {
      printf("<B>'%s' Unrecognized Identifier </B>\n",pSearchCenter->queryName);
    } else {
      printf("<B>'%s' SOFTWARE ERROR </B>\n",pSearchCenter->queryName);
    }
  }

  if (table_file_handle) {
    fclose(table_file_handle);

    printf("\n\n<a href=\"tmp/%s/table_%s.db\">Display this table as a text file</a>\n",tmpdir,tmpdir);

    /* Now write a VOTable to vo__file_name */
    sprintf(cmdStr,"%s/votable -i %s -o %s\n",binaries,table_file_name,vo_file_name);
    result = system(cmdStr);
    if (result != 0) {
      printf("ERROR: result %d(d) %x(x) %s creating the votable\n",result,result,strerror(errno));
    } else {
      statResult = stat(vo_file_name_gz,&statbuf);
      if (statResult == 0) {
        unlink(vo_file_name_gz);
      }

      sprintf(cmdStr,"%s/gzip %s\n",lsbin,vo_file_name);
      result = system(cmdStr);
      printf("<a href=\"tmp/%s/table_%s.xml.gz\">Display this table as a VOTable</a>\n",tmpdir,tmpdir);
    }
  }

  if (object_file_handle) {
    fclose(object_file_handle);
    if (printedObjects != 0) {
      printf("<a href=\"lightcurve_data.php?dbfilename=/tmp/%s/object_%s.db&vofilename=/tmp/%s/object_%s.xml.gz&tmpdir=%s&REF=object&source=%s&rematch=%s&GetData=Download+Data+in+tab-delimited+or+VOtable+format\" target=\"results_frame\">Download all points in table form.</a>\n",tmpdir,tmpdir,tmpdir,tmpdir,tmpdir,source,rematchString);
    }
  }

  printf("</pre>\n");

  for (index = 0; index < searchCenterSize; index++) {
    pSearchCenter = &search_center_table[index];
    free(pSearchCenter->queryName);
    free(pSearchCenter->queryNameurlencode);
  }

  free(search_center_table);
  if (target_table != NULL) {
    free(target_table);
  }
  if (pMagnitudeTable != NULL) {
    free(pMagnitudeTable);
  }
  if (pFileStarImageTable != NULL) {
    free(pFileStarImageTable);
  }
  if (pMagnitudeTable1 != NULL) {
    free(pMagnitudeTable1);
  }
  if (pMagnitudeTable2 != NULL) {
    free(pMagnitudeTable2);
  }

  if (pNoneMagnitudeTable != NULL) {
    free(pNoneMagnitudeTable);
  }
  if (substitutionTable != NULL) {
    free(substitutionTable);
    substitutionTable = NULL;
  }

  FreeGalaxyCommon(pGalaxyCommon,0);
  if ((RELEASE_EXPERIMENTAL != 0) && (catalogNumber == CATALOG_EXPERIMENTAL)) {
    FreeFileCommon(pFileCommon0,0);
    FreeFileCommon(pFileCommon1,0);
    FreeFileCommon(pFileCommon2,0);
  } else {
    FreeFileCommon(pFileCommon0,0);
  }

  mysql_close(pPhotConnection);
  mysql_close(pConnection);
  return 0;
}

/* May  5, 2010 Edward J. Los - Adapted from find_lightcurves.c and asas_cat_input provided by
 *                              Grzegorz Pojmanski <gp@astrouw.edu.pl>  ASAS source author on Apr 26, 2010
 * Jun 15, 2010 Edward J. Los - Correct cone search degree-arcsec problem and correct the check for the high gsc_bin_index.
 * Jun 16, 2010 Edward J. Los - Improve the positional estimate for all objects.
 * Sep  6, 2010 Edward J. Los - Correct for Solaris build.  Allow the photometry database to be in the scanner database
 * Jul 18, 2011 Edward J. Los - Provide the search list in tabular and VOtable format
 * Jul 18, 2011 Edward J. Los - Change THUMBNAIL_SIZEUNITS to THUMBNAIL_PIXELS
 *                              Make THUMBNAIL_ARCSEC the default
 * Nov 29, 2011 Edward J. Los - Do not display sequestered series
 * Feb 21, 2012 Edward J. Los - Correct errors when no objects are found.
 *                              Display the number of good points
 * Mar 13, 2012 Edward J. Los - Support region-specific authentication
 * Jul 18, 2012 Edward J. Los - Make sure that all calls to GetRef do not return spaces in object names
 * Sep 10, 2012 Edward J. Los - Reorganize code to optimize magnitude file support.
 * Oct 30, 2012 Edward J. Los - For APASS and DASCH objects, use the coordinates embedded in REF.
 * Nov  9, 2012 Edward J. Los - Allow unmatched objects to have only one point regardless of what the user requests.
 * Nov 23, 2012 Edward J. Los - Support the new APASS catalog with GSC and Tycho-2 placeholders.
 * Dec 19, 2012 Edward J. Los - Change FILTER_AMASK_PLOT to FILTER_AMASK_PLOT2
 * Feb 11, 2013 Edward J. Los - Add the magnitude of the nearest catalog star.
 * Mar 30, 2013 Edward J. Los - Support new data release mechanism
 * Apr  5, 2013 Edward J. Los - Correct search to use the J2000 position for catalog objects.
 * Jul 17, 2013 Edward J. Los - Accept new objects with only one point for completeness
 * Oct  4, 2013 Edward J. Los - Add code to dump magnitude measurements (SUPPORT_LIGHTCURVE_TABLE)
 * Jan 15, 2014 Edward J. Los - Allow tabs to be coordinate delimiters.
 * Mar 24, 2014 Edward J. Los - Covert tabs to spaces for coordinates - the previous fix is incorrect because the tabs foul up the table_*.db starbase file
 * Sep 26, 2014 Edward J. Los - use the urlencode version of the object name for Simbad queries.
 * Sep 30, 2014 Edward J. Los - do not over-write the gsc_bin_index for DASCH objects to avoid round-off errors.
 * Feb 16, 2014 Edward J. Los - Add timeAccuracy to the short form lightcurve table
 * Feb 26, 2014 Edward J. Los - Remove SUPPORT_SHORT_LIGHTCURVE_TABLE
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * Feb 29, 2016 Edward J. Los - Add "source" to the lightcurve_data invocation.
 * Mar 24, 2017 Edward J. Los - Change GetStarEntry to GetStarEntry2 to move from the stars table to the starcatalog table.
 * Dec 20, 2017 Edward J. Los - Call GetStarEntry3 to query the new identifer table for coordinates before querying Simbad
 * Jan 23, 2018 Edward J. Los - Support the merged experimental table: add catalogNumber to WriteStarbaseRecord.
 * Mar 19, 2018 Edward J. Los - Support direct read of the "experimental" binary database.
 * May  8, 2018 Edward J. Los - correct gsc_bin_index when only limiting magnitudes are available.
 * Aug 13, 2018 Edward J. Los   Define enableRematch to optimize location of transients
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 * Jan  4, 2019 Edward J. Los   Add verbose to LocateNoneImages and LoadNoneImages to display the files searched.
 * Feb  4, 2019 Edward J. Los   Switch simpos queries to the CFA address.
 */
