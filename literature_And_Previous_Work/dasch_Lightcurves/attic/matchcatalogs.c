// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* matchcatalogs.c
 *
 *  Compare the GSC2.3.2 catalog with the Kepler Input Catalog
 *  Optionally create a starbase table the joined results
 *  Optionally create a merge catalog with Kepler entries substituted by GSC2.3.2 entries
 *
 *  APASS          matchcatalogs -c a -t 2.0 -k -m -w -f /dasch/Pipeline/catalogs/ucac3.dat -s /dasch/Pipeline/catalogs/apass_temp.dat
 *  YB6            matchcatalogs -c y -t 1.0 -f /dasch/Pipeline/catalogs/yb6.dat -s /dasch/Pipeline/catalogs/kepler.dat
 *  UCAC3          matchcatalogs -c u -t 1.0 -f /dasch/Pipeline/catalogs/ucac3.dat -s /dasch/Pipeline/catalogs/kepler.dat
 *  GSC            matchcatalogs -c g -t 0.2 -f /dasch/Pipeline/catalogs/gsc232bin.dat -s /dasch/Pipeline/catalogs/kepler.dat
 *  KEPLER & APASS matchcatalogs -c a -t 2.0 -k -m -w -f /dasch/Pipeline/catalogs/kepler.dat -s /dasch/Pipeline/catalogs/apass.dat
 *  APASS & KEPLER matchcatalogs -c a -t 2.0 -k -m -w -f /dasch/Pipeline/catalogs/kepler.dat -s /dasch/Pipeline/catalogs/apass.dat
 *  EXPERIMENTAL   matchcatalogs -c e -t 2.0 -k -m -w -f /dasch/Pipeline/catalogs/gsc232bin.dat -s /dasch/Pipeline/catalogs/apass.dat
 *                 The experimental apass catalog contains GSC colorterms.
 *  UCAC3&UCAC4    matchcatalogs -c a -t 0.2 -k -m  -f /dasch/Pipeline/catalogs/ucac3.dat -s /dasch/Pipeline/catalogs/ucac4.dat
 *  GSC p.motion   matchcatalogs -c g -t 1.0 -p -m -w -f /dasch/Pipeline/catalogs/ucac4.dat -s /dasch/Pipeline/catalogs/gsc232bin.dat
 *  KEPLER p.mot   matchcatalogs -c k -t 1.0 -p -m -w -f /dasch/Pipeline/catalogs/ucac4.dat -s /dasch/Pipeline/catalogs/kepler.dat
 *  APASS pos+mot  matchcatalogs -c a -t 2.0 -l -k -m -w -f /dasch/Pipeline/catalogs/ucac4.dat -s /dasch/Pipeline/catalogs/apass_temp2.dat
 * 
 *  gcc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm  matchcatalogs.c pipelineutils.a -ltable -lutil -lwcs  -L/usr/lib/mysql -lmysqlclient -o matchcatalogs -lpthread -ldl
 * 
 *  matchcatalogs
 *
 *  Aug 14, 2009 Edward J. Los - Initial version
 *  Aug 28, 2009 Edward J. Los - produce a merge catalog
 *  Nov 13, 2010 Edward J. Los - Optionally replace the GSC catalog with the UCAC3 catalog.  Replace the Kepler proper motions with the UCAC3 proper motions
 *  Jul 28, 2011 Edward J. Los - Parameterize this route and remove compile flags.
 *  Jan 28, 2013 Edward J. Los - Support ucac4, add proper motion qualifier plus position qualifier
 *  Aug  7, 2013 Edward J. Los - Move the ucac4 matching function to matchucac4 to support APASS
 *  Aug 26, 2013 Edward J. Los - Use unique names for output files.
 *  Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 *                               
 */   

#include <math.h>
#include <time.h>
#include <errno.h>
#include "table.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"


#if 0
#define DEBUG_START_INDEX 110623744
#define DEBUG_END_INDEX   110624768
#endif

#define CATALOG_TYPE_NONE         0
#define CATALOG_TYPE_APASS        1
#define CATALOG_TYPE_YB6          2
#define CATALOG_TYPE_UCAC3        3
#define CATALOG_TYPE_GSC          4
#define CATALOG_TYPE_EXPERIMENTAL 5
#define CATALOG_TYPE_KEPLER       6
#define CATALOG_TYPE_GAIA         7
#define CATALOG_TYPE_ATLAS        8
#define CATALOG_TYPE_MAX          9
#define MAX_BUFFER 512
#define INDEX_READ_SIZE 1000000
char *catalogTypeString[CATALOG_TYPE_MAX] = {"none","apass","yb6","ucac3","gsc232","experimental","kepler","gaia","atlas"};

/* The observed Dec RMS is 0.000069 degrees in Dec and 0.000062 degrees in ra for 2985 matches */
extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;


int main(int argc,char *argv[])
{
  char *argstr; 
  int catalogType = CATALOG_TYPE_NONE;
  char typechar;
  double tolerance = 0;
  char firstdatname[MAX_BUFFER];
  char firstidxname[MAX_BUFFER];
  char seconddatname[MAX_BUFFER];
  char secondidxname[MAX_BUFFER];
  char *charPtr;

  char mergedatname[MAX_BUFFER] = "/dasch/Pipeline/catalogs/merge.dat";
  char mergeidxname[MAX_BUFFER] = "/dasch/Pipeline/catalogs/merge.idx";
  char outname[MAX_BUFFER] = "/dasch/Pipeline/catalogs/matchgsckepler.db";
  FILE * outHandle;
  File gscCatalogHandle;
  File gscIndexHandle;
  File keplerCatalogHandle;
  File keplerIndexHandle;
  File mergeCatalogHandle;
  File mergeIndexHandle;
  int errorFlag = 0;
#ifdef DEBUG_START_INDEX
#define IGNORE_OFFSET_CHECK
  int base_gsc_bin_index = DEBUG_START_INDEX;
#else /* DEBUG_START_INDEX */
	int base_gsc_bin_index = 0;
#endif /* DEBUG_START_INDEX */
  int gsc_bin_index = 0;
  int populatedBinCount = 0;
  int bin_index;
  PSTARINDEX pKeplerIndexTable = NULL;
  PSTARINDEX pKeplerIndex;
  PSTARINDEX pGscIndexTable = NULL;
  PSTARINDEX pGscIndex;

  PGSCIMAGE pGscImageTable = NULL;
  PGSCIMAGE pGscImage;
  PGSCIMAGE pBestGscImage;
  double bestDrad;
  double diffPM;
  PGSCIMAGE pKeplerImageTable = NULL;
  PGSCIMAGE pKeplerImage;
  int imageAllocationCount = 0;
  int matchCount = 0;
  int keplerCount = 0;

  int keplerReadItems;
  int keplerReadTotal = 0;
  int gscReadItems;
  int gscReadTotal = 0;
  int mergeWriteItems;
  int keplerReadCount;
  int gscReadCount;
  int mergeWriteCount;

  off_t mergeSeekAddress = 0;
  int keplerIndex;
  int gscIndex;
  double deltaDec;
  double deltaRA;
  double drad;
  double factor;
  time_t startTime;
  time_t curTime;
  char cmdchar;
  int nvals;
  int writeMergeCatalog = 0;
  int writeMatchStarbase = 0;
  int includeKeplerColor = 0;
  int includeProperMotions = 0;
  int includePositions = 0;
	char keplerREF[MAX_REF];
	char gscREF[MAX_REF];
	int debugPrintFlag = 0;

#if (defined(DEBUG_START_INDEX) || defined(DEBUG_END_INDEX))
	printf("ERROR: DEBUG_START_INDEX or DEBUG_END_INDEX is defined\n");
#endif
  firstdatname[0] = 0;
  seconddatname[0] = 0;

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: stray argument %s\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {



        case 'c': /* catalog type */
        case 'C':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            typechar = *(*++argv);
            switch (typechar) {
            case 'a':
              catalogType = CATALOG_TYPE_APASS;
              break;
            case 'y':
              catalogType = CATALOG_TYPE_YB6;
              break;
            case 'u':
              catalogType = CATALOG_TYPE_UCAC3;
              break;
            case 'g':
              catalogType = CATALOG_TYPE_GSC;
              break;
            case 'e':
              catalogType = CATALOG_TYPE_EXPERIMENTAL;
              break;
            case 'k':
              catalogType = CATALOG_TYPE_KEPLER;
              break;
            case 't':
              catalogType = CATALOG_TYPE_GAIA;
              break;
            case 'l':
              catalogType = CATALOG_TYPE_ATLAS;
              break;
            default:
              fprintf(stderr,"ERROR: unrecognized catalog type %c\n",typechar);
              break;
            }
          }
          break;


        case 't': /* Tolerance in arcsec */
        case 'T':
          argc--;
          nvals = sscanf(*++argv,"%lf",&tolerance);
          if (nvals != 1) {
            fprintf(stderr,"ERROR: Can not decode the tolerance\n");
            errorFlag = 1;
          } else {
            tolerance = tolerance/3600.;
          }
          break;


        case 'f': /* first catalog name */
        case 'F':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(firstdatname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 's': /* second catalog name */
        case 'S':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(seconddatname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'w':
        case 'W':
          writeMergeCatalog = 1;
          break;

        case 'm':
        case 'M':
          writeMatchStarbase = 1;
          break;

        case 'k':
        case 'K':
          includeKeplerColor = 1;
          break;

        case 'p':
        case 'P':
          includeProperMotions = 1;
          break;

        case 'l':
        case 'L':
          includePositions = 1;
          break;


        default:
          printf("* illegal command -%c-",cmdchar);
          errorFlag = 1;
          break;
        }
        
      }
    }
  }
  if (catalogType == CATALOG_TYPE_NONE) {
    printf("ERROR: missing catalog type\n");
    errorFlag = 1;
  }
  if (tolerance == 0) {
    printf("ERROR: tolerance is not specified\n");
    errorFlag = 1;
  }
  if (firstdatname[0] == 0) {
    printf("ERROR: first catalog name is not specified\n");
    errorFlag  = 1;
  } else {
    strcpy(firstidxname,firstdatname);
    charPtr = strstr(firstidxname,".dat");
    if (charPtr == NULL) {
      printf("ERROR: first catalog name does not end in '.dat'\n");
      errorFlag = 1;
    } else {
      *charPtr = 0;
      strcat(firstidxname,".idx");
    }
  }
  if (seconddatname[0] == 0) {
    printf("ERROR: second catalog name is not specified\n");
    errorFlag  = 1;
  } else {
    strcpy(secondidxname,seconddatname);
    charPtr = strstr(secondidxname,".dat");
    if (charPtr == NULL) {
      printf("ERROR: second catalog name does not end in '.dat'\n");
      errorFlag = 1;
    } else {
      *charPtr = 0;
      strcat(secondidxname,".idx");
    }
  }
    


  if (errorFlag) {
    printf("Usage: matchcatalogs \n");
    printf("       where -c <catalogtype> where catalogtype = 'a' for APASS\n");
    printf("                                                  'y; for YB6\n");
    printf("                                                  'u' for UCAC3\n");
    printf("                                                  'g' for GSC2.3.2\n");
    printf("                                                  'k' for KEPLER\n");
    printf("                                                  't' for Tycho-Gaia Astrometric Solution\n");
    printf("                                                  'l' for ATLAS refcat2\n");
    printf("                                                  'e' for experimental (APASS w/GSC color)\n");
    printf("             -k for INCLUDE_KEPLER_COLOR \n");
    printf("             -p for include proper motions \n");
    printf("             -l for include positions from ucac4 unless class is nonzero \n");
    printf("             -m for WRITE_MATCH_STARBASE (/dasch/Pipeline/catalogs/matchgsckepler.db)\n");
    printf("             -w for WRITE_MERGE_CATALOG (/dasch/Pipeline/catalogs/merge.dat)  \n");
    printf("             -t <tolerance> where tolerance is the match tolerance in arcsec\n");    
    printf("             -f for first catalog name  (G), eg /dasch/Pipeline/catalogs/ucac3.dat (GSC catalog)\n");
    printf("             -s for second catalog name (K), eg /dasch/Pipeline/catalogs/apass_temp.dat (Kepler catalog)\n");
    printf("             Note: the second (kepler) catalog determines the output count and merged index\n");

    return(-1);
  }


  time(&startTime);

  pKeplerIndexTable = (PSTARINDEX)calloc(INDEX_READ_SIZE,sizeof(STARINDEX));
  pGscIndexTable = (PSTARINDEX)calloc(INDEX_READ_SIZE,sizeof(STARINDEX));


  keplerCatalogHandle = Open(seconddatname,"r");
  if (keplerCatalogHandle == NULL) {
    printf("ERROR Could not open catalog file %s\n",seconddatname);
    errorFlag = 1;
  } 
  keplerIndexHandle = Open(secondidxname,"r");
  if (keplerIndexHandle == NULL) {
    printf("ERROR Could not open catalog index file %s\n",secondidxname);
    errorFlag = 1;
  } 

  gscCatalogHandle = Open(firstdatname,"r");
  if (gscCatalogHandle == NULL) {
    printf("ERROR Could not open catalog file %s\n",firstdatname);
    errorFlag = 1;
  } 
  gscIndexHandle = Open(firstidxname,"r");
  if (gscIndexHandle == NULL) {
    printf("ERROR Could not open catalog index file %s\n",firstidxname);
    errorFlag = 1;
  } 
	strcpy(mergedatname,"/dasch/Pipeline/catalogs/merge");
	strcat(mergedatname,catalogTypeString[catalogType]);
	strcpy(mergeidxname,mergedatname);

	strcat(mergedatname,".dat");
	strcat(mergeidxname,".idx");
 
	strcpy(outname,"/dasch/Pipeline/catalogs/match");
	strcat(outname,catalogTypeString[catalogType]);
	strcat(outname,"ucac4.db");

  if (writeMergeCatalog) {
    mergeCatalogHandle = Open(mergedatname,"w");
    if (mergeCatalogHandle == NULL) {
      printf("Could not open file %s\n",mergedatname);
      exit(-1);
    } 

    mergeIndexHandle = Open(mergeidxname,"w");
    if (mergeIndexHandle == NULL) {
      printf("Could not open file %s\n",mergeidxname);
      exit(-1);
    } 

  }


  if (writeMatchStarbase) {
    outHandle = fopen(outname,"wt");
    if (outHandle == NULL) {
      printf("ERROR count not open output file %s\n",outname);
    }
    fprintf(outHandle,"KREF\tKra\tKdec\tKStdmag\tKcolor\tKRaPM\tKDecPM\tKclass\tKVFlag\tKMAGFlag\tGREF\tGra\tGdec\tGStdmag\tGcolor\tGRaPM\tGDecPM\tGclass\tGVFlag\tGMAGFlag\tdegDrad\tdiffPM\n");
    fprintf(outHandle,"----\t---\t----\t-------\t------\t-----\t------\t------\t------\t--------\t----\t---\t----\t-------\t------\t-----\t------\t------\t------\t--------\t-------\t------\n");

  }

  if (errorFlag) {
    exit(-1);
  }

  printf("matchcatalogs of %s %s, GSCIMAGE size %d\n matching %s and %s \n output %s total bins %d tolerance %f degrees\n",
         __DATE__,__TIME__,sizeof(GSCIMAGE),firstdatname,seconddatname,outname,pGscBin->total_gsc_bins,tolerance);

  while (1) {
    /* Read in the star index */
#ifdef DEBUG_END_INDEX
    if (base_gsc_bin_index >= DEBUG_END_INDEX) {
      break;
    }
#endif
    time(&curTime);
    curTime -= startTime;

    printf("Reading index %9d populatedBinCount %d matches %d at %d seconds\n",base_gsc_bin_index,populatedBinCount,matchCount,curTime);
    Seek(keplerIndexHandle,base_gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
    keplerReadItems = Read(keplerIndexHandle,pKeplerIndexTable,sizeof(STARINDEX),INDEX_READ_SIZE);
    if (keplerReadItems  < 1) {
      break;
    }
    Seek(gscIndexHandle,base_gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
    gscReadItems = Read(gscIndexHandle,pGscIndexTable,sizeof(STARINDEX),INDEX_READ_SIZE);
    if (gscReadItems  != keplerReadItems) {
      printf("ERROR: gscReadItems %d does not equal keplerReadItems %d at gsc_bin_index %d\n",
             gscReadItems,keplerReadItems,base_gsc_bin_index);
    }
    if (writeMergeCatalog) {

      /* We are using only the objects in the Kepler catalog so the Kepler index becomes the
         merge index */
      mergeWriteItems = Write(mergeIndexHandle,pKeplerIndexTable,sizeof(STARINDEX),keplerReadItems);
      if (mergeWriteItems != keplerReadItems) {
        fprintf(stderr,"ERROR writing the index file\n");
        exit(-1);
      }

    }

    for (bin_index = 0; bin_index < keplerReadItems; bin_index++) {
      pKeplerIndex = &pKeplerIndexTable[bin_index];
      pGscIndex = &pGscIndexTable[bin_index];
      gsc_bin_index = bin_index+base_gsc_bin_index;
      if (pKeplerIndex->numStars > 0) {
        keplerCount += pKeplerIndex->numStars;
        if ((pKeplerIndex->numStars > imageAllocationCount) ||
            (pGscIndex->numStars > imageAllocationCount)) {
          free(pGscImageTable);
          free(pKeplerImageTable);
          imageAllocationCount = pKeplerIndex->numStars;
          if (pGscIndex->numStars > imageAllocationCount) {
            imageAllocationCount = pGscIndex->numStars;
          }
          imageAllocationCount += 100;
          pKeplerImageTable = (PGSCIMAGE)calloc(imageAllocationCount,sizeof(GSCIMAGE));
          pGscImageTable    = (PGSCIMAGE)calloc(imageAllocationCount,sizeof(GSCIMAGE));
        }
        Seek(keplerCatalogHandle,pKeplerIndex->offset,SEEK_SET);
        keplerReadCount = Read(keplerCatalogHandle,pKeplerImageTable,sizeof(GSCIMAGE),pKeplerIndex->numStars);
        keplerReadTotal += keplerReadCount;
        if (keplerReadCount != pKeplerIndex->numStars) {
          printf("ERROR keplerReadCount %d %d does not agree for bin %d\n",
                 keplerReadCount,pKeplerIndex->numStars,gsc_bin_index);
        }
        if  (pGscIndex->numStars > 0) {
          populatedBinCount++;
#if 0
          printf("Kepler stars %d, gsc stars %d for bin %9d\n",pKeplerIndex->numStars,pGscIndex->numStars,gsc_bin_index);
#endif
          Seek(gscCatalogHandle,pGscIndex->offset,SEEK_SET);
          gscReadCount = Read(gscCatalogHandle,pGscImageTable,sizeof(GSCIMAGE),pGscIndex->numStars);
          gscReadTotal += gscReadCount;
          if (gscReadCount != pGscIndex->numStars) {
            printf("ERROR gscReadCount %d %d does not agree for bin %d\n",
                   gscReadCount,pGscIndex->numStars,gsc_bin_index);
          }
          for (keplerIndex = 0; keplerIndex < keplerReadCount; keplerIndex++) {
            pKeplerImage = &pKeplerImageTable[keplerIndex];
#if 0
						GetREF(pKeplerImage->REFNumber,keplerREF,1,1);
            if (strcmp(keplerREF,"N12332221614") == 0) {  /* Occurs at base index 139000000 */
              printf("At %s\n",keplerREF);
							debugPrintFlag = 1;
            } else {
							debugPrintFlag = 0;
						}
#endif
            bestDrad = 9999.0;
            pBestGscImage = &pGscImageTable[0];
            for (gscIndex = 0; gscIndex < gscReadCount; gscIndex++) {
              pGscImage = &pGscImageTable[gscIndex];
              factor = cos(DEGREES_TO_RAD* pKeplerImage->dec);
              deltaDec = pKeplerImage->dec - pGscImage->dec;
              if (deltaDec < 0) {
                deltaDec = - deltaDec;
              }
              deltaRA  = factor *(pKeplerImage->ra  - pGscImage->ra);
              if (deltaRA < 0) {
                deltaRA = -deltaRA;
              }
              if ((deltaDec < tolerance) &&
                  (deltaRA < tolerance)) {
                drad = sqrt(sqr(deltaDec)+sqr(deltaRA));
#if 0
								if (debugPrintFlag) {
									GetREF(pGscImage->REFNumber,gscREF,1,1);
									if (strcmp(gscREF,"U4125686288") == 0) {  /* Occurs at base index 139000000 */
										printf("At %s drad %f\n",gscREF,drad);
									}
								}
#endif						


                if (drad < tolerance) {
                  if (drad < bestDrad) {
                    pBestGscImage = pGscImage;
                    bestDrad = drad;
                  }
                }
              }
            }
            if (bestDrad < 9999.0) {
              matchCount++;
              pGscImage = pBestGscImage;
              if (writeMatchStarbase) {
								GetREF(pKeplerImage->REFNumber,keplerREF,0,1);
								GetREF(pGscImage->REFNumber,gscREF,0,1);

                diffPM = sqrt(sqr(pKeplerImage->RaPM-pGscImage->RaPM)+sqr(pKeplerImage->DecPM-pGscImage->DecPM));
                fprintf(outHandle,"%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\n",
                        keplerREF,
                        pKeplerImage->ra,
                        pKeplerImage->dec,
                        pKeplerImage->Stdmag,
                        pKeplerImage->color,
                        pKeplerImage->RaPM,
                        pKeplerImage->DecPM,
                        pKeplerImage->class,
                        pKeplerImage->VFlag,
                        pKeplerImage->MAGFlag,
                        gscREF,
                        pGscImage->ra,
                        pGscImage->dec,
                        pGscImage->Stdmag,
                        pGscImage->color,
                        pGscImage->RaPM,
                        pGscImage->DecPM,
                        pGscImage->class,
                        pGscImage->VFlag,
                        pGscImage->MAGFlag,
                        bestDrad,
                        diffPM);
              }
            }
            if (writeMergeCatalog) {
              switch (catalogType) {
              case CATALOG_TYPE_UCAC3:
              case CATALOG_TYPE_YB6:
              case CATALOG_TYPE_APASS:
                /* Discard Kepler proper motions */
                pKeplerImage->RaPM = 0;
                pKeplerImage->DecPM = 0;
								pKeplerImage->RaSigmaPM = 0;
								pKeplerImage->DecSigmaPM = 0;
                break;
              case CATALOG_TYPE_GSC:
              case CATALOG_TYPE_EXPERIMENTAL:
              case CATALOG_TYPE_KEPLER:
              case CATALOG_TYPE_GAIA:
              default:
                break;
              }
              if (includeProperMotions) {
                /* Discard old proper motions */
                pKeplerImage->RaPM = 0;
                pKeplerImage->DecPM = 0;
								pKeplerImage->RaSigmaPM = 0;
								pKeplerImage->DecSigmaPM = 0;
              }
              if (bestDrad < 9999.0) {
                pGscImage = pBestGscImage;
                switch (catalogType) {
                case CATALOG_TYPE_UCAC3:
                case CATALOG_TYPE_YB6:
                case CATALOG_TYPE_APASS:
                case CATALOG_TYPE_GAIA:
                  if ((includePositions != 0) &&
                      (pGscImage->class == 0)) {
                    pKeplerImage->ra = pGscImage->ra;
                    pKeplerImage->dec = pGscImage->dec;
                    
                  }

                  /* Replace only the proper motions */
                  pKeplerImage->RaPM = pGscImage->RaPM;
                  pKeplerImage->DecPM = pGscImage->DecPM;
									pKeplerImage->RaSigmaPM = pGscImage->RaSigmaPM;
									pKeplerImage->DecSigmaPM = pGscImage->DecSigmaPM;
                  break;
                case CATALOG_TYPE_EXPERIMENTAL:
                  /* Replace only the color */
                  pKeplerImage->color = pGscImage->color;
                  break;
                case CATALOG_TYPE_GSC:
                case CATALOG_TYPE_KEPLER:
                default:

                  if (includeKeplerColor) {
                    /* Replace only GSC ra, dec, Stdmag, proper motion, but not
                       color, REF, and MAGFlag (The latter two interpret the color */

                    pKeplerImage->RaPM = pGscImage->RaPM;
                    pKeplerImage->DecPM = pGscImage->DecPM;
										pKeplerImage->RaSigmaPM = pGscImage->RaSigmaPM;
										pKeplerImage->DecSigmaPM = pGscImage->DecSigmaPM;
                    if (catalogType != CATALOG_TYPE_APASS) {
                      pKeplerImage->ra = pGscImage->ra;
                      pKeplerImage->dec = pGscImage->dec;
                      pKeplerImage->Stdmag = pGscImage->Stdmag;
                      pKeplerImage->class = pGscImage->class;
                      pKeplerImage->VFlag = pGscImage->VFlag;
                    }

                  } else { /* includeKeplerColor */
                    if (includeProperMotions) {
                      /* Replace only the proper motions */
                      pKeplerImage->RaPM = pGscImage->RaPM;
                      pKeplerImage->DecPM = pGscImage->DecPM;
											pKeplerImage->RaSigmaPM = pGscImage->RaSigmaPM;
											pKeplerImage->DecSigmaPM = pGscImage->DecSigmaPM;
                    } else {
                      /* Replace the Kepler entry with the GSC entry */
                      pKeplerImage->REFNumber,pGscImage->REFNumber;
                      pKeplerImage->ra = pGscImage->ra;
                      pKeplerImage->dec = pGscImage->dec;
                      pKeplerImage->Stdmag = pGscImage->Stdmag;
                      pKeplerImage->color = pGscImage->color;
                      pKeplerImage->RaPM = pGscImage->RaPM;
                      pKeplerImage->DecPM = pGscImage->DecPM;
											pKeplerImage->RaSigmaPM = pGscImage->RaSigmaPM;
											pKeplerImage->DecSigmaPM = pGscImage->DecSigmaPM;
                      pKeplerImage->class = pGscImage->class;
                      pKeplerImage->VFlag = pGscImage->VFlag;
                      pKeplerImage->MAGFlag = pGscImage->MAGFlag;
                    }
                  }
                }
              }
            }
          }
        }
        if (writeMergeCatalog) {
          /* Make sure we are writing to the correct offset */
#ifndef IGNORE_OFFSET_CHECK
          if (pKeplerIndex->offset != mergeSeekAddress) {
            fprintf(stderr,"ERROR expect offset %d and got %d\n",pKeplerIndex->offset,mergeSeekAddress);
            exit(-1);
          }
#endif /* IGNORE_OFFSET_CHECK */
          mergeSeekAddress += (sizeof(GSCIMAGE) * pKeplerIndex->numStars);
          /* Write out the modified Kepler catalog */
          mergeWriteCount = Write(mergeCatalogHandle,pKeplerImageTable,sizeof(GSCIMAGE),pKeplerIndex->numStars);
          if (mergeWriteCount != pKeplerIndex->numStars) {
            fprintf(stderr,"ERROR writing the catalog file\n");
            exit(-1);
          }
        }
    
      }

    }


    base_gsc_bin_index += INDEX_READ_SIZE;
    


  }
  time(&curTime);
  curTime -= startTime;

  printf("Done  populatedBinCount %d matches %d for %d kepler stars at %d seconds\n",populatedBinCount,matchCount,keplerCount,curTime);
  printf("gscReadTotal (first) %d  keplerReadTotal (second) %d\n",gscReadTotal,keplerReadTotal);

  Close(gscCatalogHandle);
  Close(gscIndexHandle);
  Close(keplerCatalogHandle);
  Close(keplerIndexHandle);
  if (writeMatchStarbase) {
    fclose(outHandle);
  }
  if (writeMergeCatalog) {
    Close(mergeCatalogHandle);
    Close(mergeIndexHandle);
  }
  free(pKeplerIndexTable);
  free(pGscIndexTable);
  free(pKeplerImageTable);
  free(pGscImageTable);
  return(0);
}
