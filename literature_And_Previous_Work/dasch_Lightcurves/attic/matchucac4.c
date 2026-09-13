// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* matchcatalogs.c
 *
 * Adds proper motion to the GSC, Kepler, or APASS catalogs
 *
 *  APASS pos+mot  matchucac4 -c a -t 3.63 -l -k -m -w -f /dasch/Pipeline/catalogs/ucac4.dat -s /dasch/Pipeline/catalogs/apass_temp4.dat
 *                                    maxGhostImages 53 maxGhostSequence 103032 maxGhostBin 6060 maxGhostREFNumber 64000201050 in 130 seconds
 *                                    ucac4 Alloc 170000 Read 113780093 Flush 114001545
 *                                    proper_motion_sequence 180356
 *  
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include   matchucac4.c pipelineutils.a -L /dasch/install/lib -lm  -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql -lmysqlclient -ldl -pthread -o matchucac4
 * 
 *                  
 *  Support ability to match uncertain APASS dates (2009-2014) with UCAC4 proper motions.  
 *  The maximum UCAC4 proper motion in DEC is 10.33"/yr, or 145" in 14 years or 2.57 1/64 degree declination bins distant.  Potential matching stars can be 3 declination bins distant
 *  The APASS pixel size is 2.5", we had been using a 3.0" match radius.
 *  The maximum possible proper motion error is 45.0 mas/yr or 0.63" in 14 years.  Handle this by increasing the match radius to 3.63"
 *  Model high proper motion stars with a ghost image every 1.0 arcsec.  Increase the match radius to 3.66" to handle any hypothenuse effects.  Use multiple images to model proper motions > 0.5" in 14 years.
 *
 *  matchucac4
 *
 *  Aug  7, 2013 Edward J. Los - Adopted from matchcatalogs.c 
 *  Aug 31, 2013 Edward J. Los - Always include a zero proper motion ghost image to handle non APASS catalog objects
 *  Jul  1, 2014 Edward J. Los - Support all APASS colors and color errors.  Write a starbase file for everything with color.
 */   

#include <math.h>
#include <time.h>
#include <errno.h>
#include "table.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "kdtree.h"
#include <sys/stat.h>
#if 0
#define DEBUG_START_INDEX 6050
#define DEBUG_END_INDEX   6070
#endif
#if 0
#define DEBUG_START_INDEX    0
#define DEBUG_END_INDEX    200
#endif
/* #define DEBUG_REFNUMBER  12110000012418L */
/* #define DEBUG_WRITTEN 1 */
/* #define DEBUG_MATCH   1 */

#define DEC_BIN_OFFSET 3 /* Maximum extent of bins to search */

#define APASS_START_YEAR 2009.0
#define APASS_END_YEAR 2014.0
#define EXTENT_BIN_SIZE 1.0  /* Number of arcsec for an extent */
#define EXTENSTION_ALLOC_INCREMENT 10000 /* Extension allocation limit */

#define CATALOG_TYPE_NONE         0
#define CATALOG_TYPE_APASS        1
#define CATALOG_TYPE_YB6          2
#define CATALOG_TYPE_UCAC3        3
#define CATALOG_TYPE_GSC          4
#define CATALOG_TYPE_EXPERIMENTAL 5
#define CATALOG_TYPE_KEPLER       6
#define MAX_BUFFER 512
#define INDEX_READ_SIZE 1000000

/* The observed Dec RMS is 0.000069 degrees in Dec and 0.000062 degrees in ra for 2985 matches */
extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
#if 0
typedef struct _pgscext {
	GSCIMAGE gscImage;  /* Copy of the catalog index */
	double ra;           /* precessed right ascension */
	double dec;          /* precessed declination */
	int proper_motion_sequence;  /* Unique identifier for each proper motion trace */
	int proper_motion_count; /* Total number of proper motion images */
	int proper_motion_index; /* Index of this proper motion image */
	int gsc_bin_index;       /* GSC bin index */
	int dec_bin_index;       /* Original declination bin index */
	int matchFlag;           /* Nonzero if already matched, otherwise contains a count of stars matched */
	int matchSequence;       /* Unique match sequence identifying a single UCAC4 stars matching multiple APASS stars */
	int rejected;            /* Nonzero for a rejected duplicate */
} GSCEXT,*PGSCEXT;
#endif

typedef struct _pgscapassext {
	GSCAPASSIMAGE gscApassImage;  /* Copy of the catalog index */
	double ra;           /* precessed right ascension */
	double dec;          /* precessed declination */
	int proper_motion_sequence;  /* Unique identifier for each proper motion trace */
	int proper_motion_count; /* Total number of proper motion images */
	int proper_motion_index; /* Index of this proper motion image */
	int gsc_bin_index;       /* GSC bin index */
	int dec_bin_index;       /* Original declination bin index */
	int matchFlag;           /* Nonzero if already matched, otherwise contains a count of stars matched */
	int matchSequence;       /* Unique match sequence identifying a single UCAC4 stars matching multiple APASS stars */
	int rejected;            /* Nonzero for a rejected duplicate */
} GSCAPASSEXT,*PGSCAPASSEXT;


typedef struct _matchcommon {
	long long maxGhostREFNumber;
	long long mergeCatalogOffset;
	double tolerance; /* matching tolerance in degrees */
	double maxdupRA;
	double maxdupDec;
	int maxdupCount;
	int ucac4_loaded_bin;
	int cur_dec_bin;
	int cur_write_bin;
	int last_written_gsc_bin;
	int proper_motion_sequence; /* Number of high proper motion stars (> 1" between 2009 and 2014) */
	int matchSequence;         /* Unique match sequence identifying a single UCAC4 stars matching multiple APASS stars */
	int maxGhostImages;
	int maxGhostSequence;
	int maxGhostBin;
  char ucac4datname[MAX_BUFFER];
  char ucac4idxname[MAX_BUFFER];
  File ucac4CatalogHandle;
  File ucac4IndexHandle;
  File apassCatalogHandle;
  File apassIndexHandle;
	File mergeCatalogHandle;
	File mergeIndexHandle;
  int ucac4recordsize;
  PSTARINDEX pIndexTable;
  PGSCIMAGE pImageTable;
	int imageTableAlloc;
  PGSCAPASSIMAGE pImageTableApass;
	int imageTableApassAlloc;
	PGSCAPASSEXT pUcac4ExtTable;
  int ucac4CurCount;
 	int ucac4AllocCount;
	int ucac4ReadCount;
	int ucac4FlushCount;
	PGSCAPASSEXT pApassExtTable;
	int apassCurCount;
	int apassAllocCount;
	int apassReadCount;
	int apassFlushCount;
	int duplicateMatchCount;
	int matchCount;
	int maxMatchCount;
	int writeMatchStarbase;
	int writeMergeCatalog;
	int mergeWriteCount;
	int includeKeplerColor;
	int includeProperMotions;
	int includePositions;
  int apassColorFlag; /* 0 for GSCIMAGE; 1 for GSCAPASSIMAGE */

	FILE * outHandle;
  FILE * resultHandle;
	int catalogType;
	void *ptree; /* kdtree for everything in the buffer */

} MATCHCOMMON,*PMATCHCOMMON;

int ImageCompare(const void *first, const void *second) 
{
  int numberFirst = ((PGSCAPASSEXT)first)->gsc_bin_index;
  int numberSecond = ((PGSCAPASSEXT)second)->gsc_bin_index;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}

#ifdef DEBUG_WRITTEN
void DumpImageTable(PMATCHCOMMON pMatchCommon,char *reportString,int numStars)
{
	int index;
	PGSCIMAGE pApassImage;
	PGSCAPASSIMAGE pApassImageColor;

	for (index = 0; index < numStars;index++) {
    if (pMatchCommon->apassColorFlag == 0) {
      pApassImage = &pMatchCommon->pImageTable[index];
#ifdef DEBUG_REFNUMBER
      if ((pApassImage->REFNumber == DEBUG_REFNUMBER) && (strstr(reportString,"f Written") != NULL)) {
        printf("At REFNumber (20) index %d %lld \n",index,pApassImage->REFNumber);
      }
#endif


      printf("%s: %lld\n",reportString,pApassImage->REFNumber);
    } else {
      pApassImageColor = &pMatchCommon->pImageTableApass[index];
#ifdef DEBUG_REFNUMBER
      if ((pApassImageColor->REFNumber == DEBUG_REFNUMBER) && (strstr(reportString,"f Written") != NULL)) {
        printf("At REFNumber (21) index %d %lld \n",index,pApassImageColor->REFNumber);
      }
#endif
      printf("%s: %lld\n",reportString,pApassImageColor->REFNumber);

    }
	}
}
#endif /* DEBUG_WRITTEN */
void CheckDuplicates(PMATCHCOMMON pMatchCommon,PGSCIMAGE pImageTable,int numStars)
{

	/* This is a quick and dirty check for duplicates within a single bin */
	PGSCIMAGE pGscImage1;
	PGSCIMAGE pGscImage2;
	int index1;
	int index2;
	double factor;
	double deltaDec;
	double deltaRA;
	double drad;
	if (numStars < 2) {
		return;
	}
	for (index1 = 0; index1 < numStars; index1++) {
		pGscImage1 = &pImageTable[index1];
		factor = cos(DEGREES_TO_RAD * pGscImage1->dec);
		for (index2 = 1 ; index2 < numStars; index2++) {
			pGscImage2 = &pImageTable[index2];
			deltaDec = pGscImage1->dec - pGscImage2->dec;
			deltaRA = factor *(pGscImage1->ra - pGscImage2->ra);
			drad = sqrt(sqr(deltaDec)+sqr(deltaRA));
			if ((drad*3600.0) < pMatchCommon->tolerance) {
				printf("ERROR: drad is %d arcsec for %lld and %lld\n",drad,pGscImage1->REFNumber,pGscImage2->REFNumber);
				exit(-1);
			}
		}
	}
	if (pMatchCommon->maxdupCount < numStars) {
		pMatchCommon->maxdupCount = numStars;
		pMatchCommon->maxdupRA = pGscImage1->ra;
		pMatchCommon->maxdupDec = pGscImage1->dec;
	}


}



void WriteMatchedStars(PMATCHCOMMON pMatchCommon)
{
	STARINDEX apassIndex;
	PSTARINDEX pApassIndex = &apassIndex;
	PGSCAPASSEXT pApassExt;
	PGSCAPASSEXT pApassExt2;
	PGSCAPASSEXT pBestApassExt;
	PBININDEX pBinIndex;
	PGSCIMAGE pGscImage;
	int index1;
	int index2;
	int mergeWriteCount;
	int mergeWriteItems;
	int max_gsc_bin_index;
	int imageTableSize = 0;
	int debugPrintFlag = 0;
#ifdef DEBUG_REFNUMBER
	for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
		pApassExt = &pMatchCommon->pApassExtTable[index1];
		if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
			printf("At REFNumber (10) index %d %lld sequence %d\n",index1,pApassExt->gscApassImage.REFNumber,pApassExt->matchSequence);
		}
	}
#endif


	if ((pMatchCommon->cur_write_bin < 0) ||
			(pMatchCommon->cur_write_bin >= pGscBin->dec_bins)) {
		return; /* Nothing to load here */
	}
	pBinIndex = &pGscBin->pBinMasterIndex[pMatchCommon->cur_write_bin];
	max_gsc_bin_index = pBinIndex->startBin + pBinIndex->numBins-1;
	memset(pApassIndex,0,sizeof(STARINDEX));
	if ((pMatchCommon->writeMergeCatalog != 0) &&
			(pMatchCommon->apassCurCount > 0)) {
		/* Sort all of our results by gsc_bin_index */
		qsort((void*)pMatchCommon->pApassExtTable,pMatchCommon->apassCurCount,sizeof(GSCAPASSEXT),ImageCompare);
#ifdef DEBUG_REFNUMBER
		for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
			pApassExt = &pMatchCommon->pApassExtTable[index1];
			if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
				printf("At REFNumber (11) index %d %lld apassCurCount %d sequence %d\n",index1,pApassExt->gscApassImage.REFNumber,pMatchCommon->apassCurCount,pApassExt->matchSequence);
				debugPrintFlag = 1;
			}
		}
#endif
		pApassIndex->offset = pMatchCommon->mergeCatalogOffset;
		pApassIndex->binNumber = pMatchCommon->last_written_gsc_bin+1;
		pApassIndex->numStars = 0;
		
		for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
			pApassExt = &pMatchCommon->pApassExtTable[index1];

#if 0
			if (debugPrintFlag) {
				printf("At REFNumber (12) index %d %lld apassCurCount %d\n",index1,pApassExt->gscApassImage.REFNumber,pMatchCommon->apassCurCount);

			}
#endif
#ifdef DEBUG_REFNUMBER
			if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
				printf("At REFNumber (8) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
			}
#endif

			if (pApassExt->rejected) {
				continue;
			}
			if (pApassExt->matchFlag > 1) {
#ifdef DEBUG_REFNUMBER
				if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
					printf("At REFNumber (3) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
				}
#endif
				/* Here we have duplicate matches.  Decide on the best match and reject all others */
				pBestApassExt = pApassExt;
				for (index2 = index1+1; index2 < pMatchCommon->apassCurCount; index2++) {
					pApassExt2 = &pMatchCommon->pApassExtTable[index2];

					if ((pApassExt2->rejected == 0) && 
							(pApassExt2->matchSequence == pBestApassExt->matchSequence)) {
#if 0
						if (pBestApassExt->matchSequence == 1501) {
							printf("At sequence %d\n",pBestApassExt->matchSequence);
						}
#endif				

						/* Select the one with a valid magnitude */
						if ((pBestApassExt->gscApassImage.Stdmag < 90.0) &&
								(pApassExt2->gscApassImage.Stdmag > 90.0)) {
							pApassExt2->rejected = 1;
#if (defined(DEBUG_WRITTEN) || defined(DEBUG_MATCH))
							printf("d Rejected: %lld sequence %d\n",pApassExt2->gscApassImage.REFNumber,pApassExt2->matchSequence);
#endif /* DEBUG_WRITTEN */
							
						} else if ((pBestApassExt->gscApassImage.Stdmag > 90.0) &&
											 (pApassExt2->gscApassImage.Stdmag < 90.0)) {
							pBestApassExt->rejected = 1;
#if (defined(DEBUG_WRITTEN) || defined(DEBUG_MATCH))
							printf("d Rejected: %lld sequence %d\n",pBestApassExt->gscApassImage.REFNumber,pBestApassExt->matchSequence);
#endif /* DEBUG_WRITTEN */
							pBestApassExt = pApassExt2;
						} else if (GetREFType(pBestApassExt->gscApassImage.REFNumber) == REF_TYPE_APASS) {
							/* Select the APASS reference type over TYCHO and GSC reference types */
#if (defined(DEBUG_WRITTEN) || defined(DEBUG_MATCH))
							printf("d Rejected: %lld sequence %d\n",pApassExt2->gscApassImage.REFNumber,pApassExt2->matchSequence);
#endif /* DEBUG_WRITTEN */
							pApassExt2->rejected = 1;
						} else if (GetREFType(pApassExt2->gscApassImage.REFNumber) == REF_TYPE_APASS) {
#if (defined(DEBUG_WRITTEN) || defined(DEBUG_MATCH))
							printf("d Rejected: %lld sequence %d\n",pBestApassExt->gscApassImage.REFNumber,pBestApassExt->matchSequence);
#endif /* DEBUG_WRITTEN */
							pBestApassExt->rejected = 1;
							pBestApassExt = pApassExt2;
						} else {
							/* By default, just reject the duplicate */
#if (defined(DEBUG_WRITTEN) || defined(DEBUG_MATCH))
							printf("d Rejected: %lld sequence %d\n",pApassExt2->gscApassImage.REFNumber,pApassExt2->matchSequence);
#endif /* DEBUG_WRITTEN */
							pApassExt2->rejected = 1;
						}
					}
				}
			}
			if (pApassExt->rejected) {
				continue;
			}

			if (pApassExt->gsc_bin_index <= pMatchCommon->last_written_gsc_bin) {
				printf("ERROR: bin index %d already written (1), last %d\n",pApassExt->gsc_bin_index,pMatchCommon->last_written_gsc_bin);
				exit(-1);
			}
			if ((pMatchCommon->last_written_gsc_bin+1) != pApassExt->gsc_bin_index) {
				/* Write null index entries if necessary */
				while (pApassIndex->binNumber < pApassExt->gsc_bin_index) {
					if (pApassIndex->binNumber > max_gsc_bin_index) {
						break;
					}
					mergeWriteItems = Write(pMatchCommon->mergeIndexHandle,pApassIndex,sizeof(STARINDEX),1);
					if (mergeWriteItems != 1) {
						printf("ERROR writing the merge index file\n");
						exit(-1);
					}
					if (pApassIndex->numStars > 0) {
#ifdef DEBUG_WRITTEN
						DumpImageTable(pMatchCommon,"f Written",pApassIndex->numStars);
#endif /* DEBUG_WRITTEN */

#ifdef DEBUG_REFNUMBER
            {
              int index3;
              for (index3 = 0; index3 < pApassIndex->numStars; index3++) {
                PGSCIMAGE pGscImage3 = &pMatchCommon->pImageTable[index3];


                if (pGscImage3->REFNumber == DEBUG_REFNUMBER) {
                  printf("At REFNumber (23) index %d %lld\n",index1,pGscImage3->REFNumber);
                }
              }
            }
#endif /* DEBUG_REFNUMBER */

						
						mergeWriteCount = Write(pMatchCommon->mergeCatalogHandle,pMatchCommon->pImageTable,sizeof(GSCIMAGE),pApassIndex->numStars);
						if (mergeWriteCount != pApassIndex->numStars) {
							fprintf(stderr,"ERROR writing the catalog file\n");
							exit(-1);
						}
						pMatchCommon->mergeWriteCount += mergeWriteCount;
						imageTableSize = 0;
					}

					pMatchCommon->last_written_gsc_bin++;
#if 0
					if (pMatchCommon->last_written_gsc_bin == 91163789) {
						printf("At last_written_gsc_bin %d\n",pMatchCommon->last_written_gsc_bin);
					}
#endif
					pMatchCommon->mergeCatalogOffset += (pApassIndex->numStars * sizeof(GSCIMAGE));
					pApassIndex->offset = pMatchCommon->mergeCatalogOffset;
					pApassIndex->binNumber = pMatchCommon->last_written_gsc_bin+1;
					pApassIndex->numStars = 0;
					
				}
			}
			if (pApassIndex->binNumber > max_gsc_bin_index) {
				break;
			}
			/* If we got here, save the result in our image buffer */
			if (imageTableSize+1 > pMatchCommon->imageTableAlloc) {
				printf("ERROR: image table size %d exceeds allocation %d\n",imageTableSize);
				exit(-1);
			}
			pGscImage = &pMatchCommon->pImageTable[imageTableSize];
#ifdef DEBUG_REFNUMBER
			if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
				printf("At REFNumber (24) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
			}
#endif
			/* memcpy(pGscImage,&pApassExt->gscApassImage,sizeof(GSCIMAGE)); Bug! */
      pGscImage->REFNumber = pApassExt->gscApassImage.REFNumber;
      pGscImage->ra = pApassExt->gscApassImage.ra;
      pGscImage->dec = pApassExt->gscApassImage.dec;
      pGscImage->Stdmag = pApassExt->gscApassImage.Stdmag;
      pGscImage->color = pApassExt->gscApassImage.color;
      pGscImage->RaPM = pApassExt->gscApassImage.RaPM;
      pGscImage->DecPM = pApassExt->gscApassImage.DecPM;
      pGscImage->RaSigmaPM = pApassExt->gscApassImage.RaSigmaPM;
      pGscImage->DecSigmaPM = pApassExt->gscApassImage.DecSigmaPM;
      pGscImage->class = pApassExt->gscApassImage.class;
      pGscImage->VFlag = pApassExt->gscApassImage.VFlag;
      pGscImage->MAGFlag = pApassExt->gscApassImage.MAGFlag;
      if (pMatchCommon->resultHandle != NULL) {
        char REF[MAX_REF];
        long long REFNumber;
        int refType;
        GetREF(pGscImage->REFNumber,REF,1,1);
        GetREFNumber(REF,&REFNumber,&refType,1,0);
        if (refType == REF_TYPE_APASS) {
          fprintf(pMatchCommon->resultHandle,"%s\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                  REF,
                  pGscImage->ra,
                  pGscImage->dec,
                  pGscImage->Stdmag,
                  pGscImage->color,
                  pGscImage->RaPM,
                  pGscImage->DecPM,
                  pGscImage->RaSigmaPM,
                  pGscImage->DecSigmaPM,
                  pGscImage->class,
                  pGscImage->VFlag,
                  pGscImage->MAGFlag,
                  pGscImage->flag,
                  pApassExt->gscApassImage.gmag,
                  pApassExt->gscApassImage.rmag,
                  pApassExt->gscApassImage.imag,
                  pApassExt->gscApassImage.vmagerr,
                  pApassExt->gscApassImage.bmagerr,
                  pApassExt->gscApassImage.gmagerr,
                  pApassExt->gscApassImage.rmagerr,
                  pApassExt->gscApassImage.imagerr);
        }



      }

#ifdef DEBUG_WRITTEN
			printf("e Copied: %lld\n",pApassExt->gscApassImage.REFNumber);
#endif /* DEBUG_WRITTEN */
			pApassIndex->numStars++;
			imageTableSize++;			
		}
	}
	if ((pMatchCommon->writeMergeCatalog != 0) &&
			(pApassIndex->numStars > 0)) {
		/* Do a final sanity check for duplicates */
		CheckDuplicates(pMatchCommon,pMatchCommon->pImageTable,pApassIndex->numStars);
		/* We need to write out the last buffer */
		mergeWriteItems = Write(pMatchCommon->mergeIndexHandle,pApassIndex,sizeof(STARINDEX),1);
		if (mergeWriteItems != 1) {
			printf("ERROR writing the merge index file\n");
			exit(-1);
		}

#ifdef DEBUG_REFNUMBER
    {
      int index3;
      for (index3 = 0; index3 < pApassIndex->numStars; index3++) {
        PGSCIMAGE pGscImage3 = &pMatchCommon->pImageTable[index3];


        if (pGscImage3->REFNumber == DEBUG_REFNUMBER) {
          printf("At REFNumber (22) index %d %lld\n",index1,pGscImage3->REFNumber);
        }
      }
    }
#endif /* DEBUG_REFNUMBER */


		mergeWriteCount = Write(pMatchCommon->mergeCatalogHandle,pMatchCommon->pImageTable,sizeof(GSCIMAGE),pApassIndex->numStars);
		if (mergeWriteCount != pApassIndex->numStars) {
			fprintf(stderr,"ERROR writing the catalog file\n");
			exit(-1);
		}
		pMatchCommon->mergeWriteCount += mergeWriteCount;
		pMatchCommon->last_written_gsc_bin++;
#if 0
		if (pMatchCommon->last_written_gsc_bin == 91163789) {
			printf("At last_written_gsc_bin %d\n",pMatchCommon->last_written_gsc_bin);
		}
#endif
		pMatchCommon->mergeCatalogOffset += (pApassIndex->numStars * sizeof(GSCIMAGE));
	}
 

	/* Now flush the old stars */
	index1 = 0;
	while (index1 < pMatchCommon->apassCurCount) {
		pApassExt = &pMatchCommon->pApassExtTable[index1];
#if 0
		if (debugPrintFlag) {
			printf("At REFNumber (13) index %d %lld apassCurCount %d\n",index1,pApassExt->gscApassImage.REFNumber,pMatchCommon->apassCurCount);	
		}
#endif
		if (pApassExt->gsc_bin_index <= max_gsc_bin_index) {
#ifdef DEBUG_REFNUMBER
			if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
				printf("At REFNumber (5) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
			}
#endif
			
			if (index1 < (pMatchCommon->apassCurCount-1)) {
				pApassExt2 =  &pMatchCommon->pApassExtTable[pMatchCommon->apassCurCount-1];
				memcpy(pApassExt,pApassExt2,sizeof(GSCAPASSEXT));
			}
			pMatchCommon->apassFlushCount++;
			pMatchCommon->apassCurCount--;
		} else {
#ifdef DEBUG_REFNUMBER
			if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
				printf("At REFNumber (14) index %d %lld sequence %d\n",index1,pApassExt->gscApassImage.REFNumber,pApassExt->matchSequence);
			}
#endif
			if (pApassExt->gsc_bin_index <= pMatchCommon->last_written_gsc_bin) {
				printf("ERROR: bin index %d already written (2), last %d\n",pApassExt->gsc_bin_index,pMatchCommon->last_written_gsc_bin);
				exit(-1);
			}
			index1++;
		}
	}
}
void FlushUcac4Stars(PMATCHCOMMON pMatchCommon)
{
	PGSCAPASSEXT pUcac4Ext;
	PGSCAPASSEXT pUcac4Ext2;
	PBININDEX pBinIndex;
	double minDeclination;
	int index;
	int index1;
#ifdef DEBUG_REFNUMBER
	for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
		PGSCAPASSEXT pApassExt;
		pApassExt = &pMatchCommon->pApassExtTable[index1];
		if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
			printf("At REFNumber (17) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
		}
	}
#endif
	if ((pMatchCommon->cur_write_bin < 0) ||
			(pMatchCommon->cur_write_bin >= pGscBin->dec_bins)) {
		return; /* Nothing to load here */
	}
	pBinIndex = &pGscBin->pBinMasterIndex[pMatchCommon->cur_write_bin];
	minDeclination = pBinIndex->declination+pGscBin->bin_size;
	index = 0;
	while (index < pMatchCommon->ucac4CurCount) {
		pUcac4Ext = &pMatchCommon->pUcac4ExtTable[index];
		if (pUcac4Ext->dec < minDeclination) {
			if (index < (pMatchCommon->ucac4CurCount-1)) {
				pUcac4Ext2 =  &pMatchCommon->pUcac4ExtTable[pMatchCommon->ucac4CurCount-1];
				memcpy(pUcac4Ext,pUcac4Ext2,sizeof(GSCAPASSEXT));
			}
			pMatchCommon->ucac4FlushCount++;
			pMatchCommon->ucac4CurCount--;
		} else {
			index++;
		}
	}

#ifdef DEBUG_REFNUMBER
	for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
		PGSCAPASSEXT pApassExt;
		pApassExt = &pMatchCommon->pApassExtTable[index1];
		if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
			printf("At REFNumber (19) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
		}
	}
#endif
}
void LoadAPASSStars(PMATCHCOMMON pMatchCommon)
{
	int raBin;
	int decBin;
	PBININDEX pBinIndex;
	int readItems;
	int starCount;
	int binIndex;
	int index1;
	PSTARINDEX pStarIndex;
	PGSCAPASSEXT tmpApassExtTable;
	int readCount;
	int starIndex;
	PGSCIMAGE pGscImage = NULL;
  PGSCAPASSIMAGE pGscApassImage = NULL;
	PGSCAPASSEXT pApassExt;
#ifdef DEBUG_END_INDEX
  if (pMatchCommon->cur_dec_bin >= DEBUG_END_INDEX) {
    return;
  }
#endif
#ifdef DEBUG_REFNUMBER
	for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
		pApassExt = &pMatchCommon->pApassExtTable[index1];
		if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
			printf("At REFNumber (15) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
		}
	}
#endif
	if ((pMatchCommon->cur_dec_bin < 0) ||
			(pMatchCommon->cur_dec_bin >= pGscBin->dec_bins)) {
		return; /* Nothing to load here */
	}
	pBinIndex = &pGscBin->pBinMasterIndex[pMatchCommon->cur_dec_bin];

	if (pBinIndex->numBins > pGscBin->max_ra_bins) {
		printf("ERROR numBins exceeds max_ra_bins\n");
		exit(-1);
	}
	Seek(pMatchCommon->apassIndexHandle,pBinIndex->startBin * sizeof(STARINDEX),SEEK_SET);
	readItems = Read(pMatchCommon->apassIndexHandle,pMatchCommon->pIndexTable,sizeof(STARINDEX),pBinIndex->numBins);
	if (readItems != pBinIndex->numBins) {
		printf("ERROR reading the apass index file\n");
		exit(-1);
	}
	/* Find out how many stars we need to read */
	starCount = 0;
	for (binIndex = 0; binIndex < pBinIndex->numBins; binIndex++) {
		pStarIndex = &pMatchCommon->pIndexTable[binIndex];
		starCount += pStarIndex->numStars;
	}
	if (starCount == 0) {
		/* Nothing to read, we are done */
		return;
	}

  if (pMatchCommon->apassColorFlag == 0) {
    if (starCount > pMatchCommon->imageTableAlloc) {
      if (pMatchCommon->pImageTable != 0) {
        free(pMatchCommon->pImageTable);
      }
      pMatchCommon->pImageTable = (PGSCIMAGE)calloc(starCount+EXTENSTION_ALLOC_INCREMENT,sizeof(GSCIMAGE));
      if (pMatchCommon->pImageTable == NULL) {
        printf("ERROR allocating image buffer of size %d\n",starCount);
      }
      pMatchCommon->imageTableAlloc = starCount+EXTENSTION_ALLOC_INCREMENT;
    }
  } else {
    if (starCount > pMatchCommon->imageTableApassAlloc) {
      if (pMatchCommon->pImageTableApass != 0) {
        free(pMatchCommon->pImageTableApass);
      }
      pMatchCommon->pImageTableApass = (PGSCAPASSIMAGE)calloc(starCount+EXTENSTION_ALLOC_INCREMENT,sizeof(GSCAPASSIMAGE));
      if (pMatchCommon->pImageTableApass == NULL) {
        printf("ERROR allocating image buffer of size %d\n",starCount);
      }
      pMatchCommon->imageTableApassAlloc = starCount+EXTENSTION_ALLOC_INCREMENT;
    }
  }


	if ((pMatchCommon->apassCurCount+starCount+1) >= pMatchCommon->apassAllocCount) {
		/* Need to reallocate the APASS extenstion list */
		pMatchCommon->apassAllocCount += starCount+EXTENSTION_ALLOC_INCREMENT;
		tmpApassExtTable = (PGSCAPASSEXT)realloc(pMatchCommon->pApassExtTable,pMatchCommon->apassAllocCount*sizeof(GSCAPASSEXT));
		if (tmpApassExtTable == NULL) {
			printf("ERROR: failed to reallocate pMatchCommon->pApassExtTable of size %d\n",pMatchCommon->apassAllocCount);
			exit(-1);
		}
		pMatchCommon->pApassExtTable = tmpApassExtTable;
		tmpApassExtTable = NULL;
	}

	Seek(pMatchCommon->apassCatalogHandle,pMatchCommon->pIndexTable->offset,SEEK_SET);
  if (pMatchCommon->apassColorFlag == 0) {
    readCount = Read(pMatchCommon->apassCatalogHandle,pMatchCommon->pImageTable,sizeof(GSCIMAGE),starCount);
  } else {
    readCount = Read(pMatchCommon->apassCatalogHandle,pMatchCommon->pImageTableApass,sizeof(GSCAPASSIMAGE),starCount);
  }

#ifdef DEBUG_WRITTEN
	DumpImageTable(pMatchCommon,"a Read",starCount);
#endif /* DEBUG_WRITTEN */

	if (readCount != starCount) {
		printf("ERROR: readCount %d %d does not agree\n",
					 readCount,starCount);
		exit(-1);		
	}
	pMatchCommon->apassReadCount += readCount;
	for (starIndex = 0; starIndex < readCount; starIndex++) {
    if (pMatchCommon->apassColorFlag == 0) {
      pGscImage = &pMatchCommon->pImageTable[starIndex];
      pGscImage->flag = 0;
#ifdef DEBUG_REFNUMBER
      if (pGscImage->REFNumber == DEBUG_REFNUMBER) {
        printf("At REFNumber (2) index %d %lld\n",pMatchCommon->apassCurCount,pGscImage->REFNumber);
      }
#endif
    } else {
      pGscApassImage = &pMatchCommon->pImageTableApass[starIndex];
      pGscApassImage->flag = 0;
#ifdef DEBUG_REFNUMBER
      if (pGscApassImage->REFNumber == DEBUG_REFNUMBER) {
        printf("At REFNumber (2) index %d %lld\n",pMatchCommon->apassCurCount,pGscApassImage->REFNumber);
      }
#endif
    }
		pApassExt = &pMatchCommon->pApassExtTable[pMatchCommon->apassCurCount];
		memset(pApassExt,0,sizeof(GSCAPASSEXT));
    if (pMatchCommon->apassColorFlag == 0) {
      pApassExt->gscApassImage.REFNumber = pGscImage->REFNumber;
      pApassExt->gscApassImage.ra = pGscImage->ra;
      pApassExt->gscApassImage.dec = pGscImage->dec;
      pApassExt->gscApassImage.Stdmag = pGscImage->Stdmag;
      pApassExt->gscApassImage.color = pGscImage->color;
      pApassExt->gscApassImage.RaPM = pGscImage->RaPM;
      pApassExt->gscApassImage.DecPM = pGscImage->DecPM;
      pApassExt->gscApassImage.RaSigmaPM = pGscImage->RaSigmaPM;
      pApassExt->gscApassImage.DecSigmaPM = pGscImage->DecSigmaPM;
      pApassExt->gscApassImage.class = pGscImage->class;
      pApassExt->gscApassImage.VFlag = pGscImage->VFlag;
      pApassExt->gscApassImage.MAGFlag = pGscImage->MAGFlag;
      pApassExt->gscApassImage.flag = pGscImage->flag;
      pApassExt->gscApassImage.gmag = 99.0;
      pApassExt->gscApassImage.rmag = 99.0;
      pApassExt->gscApassImage.imag = 99.0;
      pApassExt->gscApassImage.vmagerr = 9.999;
      pApassExt->gscApassImage.bmagerr = 9.999;
      pApassExt->gscApassImage.gmagerr = 9.999;
      pApassExt->gscApassImage.rmagerr = 9.999;
      pApassExt->gscApassImage.imagerr = 9.999;
      pApassExt->ra = pGscImage->ra;
      pApassExt->dec = pGscImage->dec;
    } else {
      memcpy(&pApassExt->gscApassImage,pGscApassImage,sizeof(GSCAPASSIMAGE));
      pApassExt->ra = pGscApassImage->ra;
      pApassExt->dec = pGscApassImage->dec;
    
    }
    pApassExt->dec_bin_index = pMatchCommon->cur_dec_bin;
    pApassExt->gsc_bin_index = GetGSCBin(pGscBin,pApassExt->ra,pApassExt->dec,&decBin,&raBin,"LoadAPASSStars");

    pMatchCommon->apassCurCount++;
		
	}

#ifdef DEBUG_REFNUMBER
	for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
		pApassExt = &pMatchCommon->pApassExtTable[index1];
		if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
			printf("At REFNumber (16) index %d %lld sequence %d\n",index1,pApassExt->gscApassImage.REFNumber,pApassExt->matchSequence);
		}
	}
#endif

	return;
}
void MatchStars(PMATCHCOMMON pMatchCommon)
{
	PGSCAPASSEXT pUcac4Ext;
	PGSCAPASSEXT pUcac4Ext2;
	PGSCAPASSEXT pBestUcac4Ext;
	PGSCAPASSEXT pApassExt;
	PGSCAPASSEXT pApassExt2;
	PGSCAPASSIMAGE pApassImage;
	PGSCAPASSIMAGE pUcac4Image;
	double x;
	double y;
	double z;
	double cosx;
	int index1;
	int index2;
	double pt[3];  
	double pos[3];
	struct kdres *presults;
	int resultSize;
	double curDistance;
	double factor;
	double bestDrad;
	double deltaDec;
	double deltaRA;
	double drad;
	char apassREF[MAX_REF];
	char ucac4REF[MAX_REF];
	int curMatchCount;
	long long tempREFNumber;
	int raBin;
	int decBin;
	if (pMatchCommon->ucac4CurCount == 0) {
		/* No Ucac4 stars, nothing to do */
		return;
	}
	/* First create a kdtree with all active ucac4 objects */
	kd_free(pMatchCommon->ptree);
	pMatchCommon->ptree = kd_create(3);

	for (index1 = 0; index1 < pMatchCommon->ucac4CurCount; index1++) {
		pUcac4Ext = &pMatchCommon->pUcac4ExtTable[index1];
		z = sin(DEGREES_TO_RAD*pUcac4Ext->dec);
		cosx = cos(DEGREES_TO_RAD*pUcac4Ext->dec);
		x = cosx*cos(DEGREES_TO_RAD*pUcac4Ext->ra);
		y = cosx*sin(DEGREES_TO_RAD*pUcac4Ext->ra);
		if (kd_insert3(pMatchCommon->ptree,x,y,z,pUcac4Ext) != 0) {
			printf("ERROR: fatal return from kd_insert3\n");
			exit(-1);
		}
	}
	/* Now go through all of the apass entries and attempt to match each one with a ucac4 object */
	for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
		pApassExt = &pMatchCommon->pApassExtTable[index1];
		pApassImage = &pApassExt->gscApassImage;
		GetREF(pApassImage->REFNumber,apassREF,0,1);
		if ((pApassExt->dec_bin_index == pMatchCommon->cur_dec_bin) &&
				(pApassExt->matchFlag == 0)) {
			pApassImage->RaPM = 0;
			pApassImage->DecPM = 0;
			pApassImage->RaSigmaPM = 0;
			pApassImage->DecSigmaPM = 0;

			z = sin(DEGREES_TO_RAD*pApassImage->dec);
			cosx = cos(DEGREES_TO_RAD*pApassImage->dec);
			x = cosx*cos(DEGREES_TO_RAD*pApassImage->ra);
			y = cosx*sin(DEGREES_TO_RAD*pApassImage->ra);
			pt[0] = x;
			pt[1] = y;
			pt[2] = z;
			presults = kd_nearest_range(pMatchCommon->ptree, pt, pMatchCommon->tolerance*DEGREES_TO_RAD );
			resultSize = kd_res_size(presults);
			bestDrad = 9999.0;
			while( !kd_res_end( presults ) ) {
				/* get the data and position of the current result item */
				pUcac4Ext = (PGSCAPASSEXT)kd_res_item( presults, pos );
				factor = cos(DEGREES_TO_RAD* pApassImage->dec);
				deltaDec = pApassImage->dec - pUcac4Ext->dec;
				if (deltaDec < 0) {
					deltaDec = - deltaDec;
				}
				deltaRA  = factor *(pApassImage->ra  - pUcac4Ext->ra);
				if (deltaRA < 0) {
					deltaRA = -deltaRA;
				}
				if ((deltaDec < pMatchCommon->tolerance) &&
						(deltaRA < pMatchCommon->tolerance)) {
					drad = sqrt(sqr(deltaDec)+sqr(deltaRA));
					if (drad < pMatchCommon->tolerance) {
						if (drad < bestDrad) {
							pBestUcac4Ext = pUcac4Ext;
							bestDrad = drad;
						}
					}
				}
				/* go to the next entry */
				kd_res_next( presults );
			}
			
			if (bestDrad < 9999.0) {
				pMatchCommon->matchCount++;
				pUcac4Ext = pBestUcac4Ext;

				pUcac4Image = &pBestUcac4Ext->gscApassImage;
				if (pMatchCommon->writeMatchStarbase) {
					GetREF(pUcac4Image->REFNumber,ucac4REF,0,1);



					fprintf(pMatchCommon->outHandle,"%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
									apassREF,
									pApassImage->ra,
									pApassImage->dec,
									pApassImage->Stdmag,
									pApassImage->color,
									pApassImage->RaPM,
									pApassImage->DecPM,
									pApassImage->class,
									pApassImage->VFlag,
									pApassImage->MAGFlag,
									ucac4REF,
									pUcac4Image->ra,
									pUcac4Image->dec,
									pUcac4Image->Stdmag,
									pUcac4Image->color,
									pUcac4Image->RaPM,
									pUcac4Image->DecPM,
									pUcac4Image->class,
									pUcac4Image->VFlag,
									pUcac4Image->MAGFlag,
									bestDrad,
									pUcac4Ext->matchFlag,
                  pApassImage->gmag,
                  pApassImage->rmag,
                  pApassImage->imag,
                  pApassImage->vmagerr,
                  pApassImage->bmagerr,
                  pApassImage->gmagerr,
                  pApassImage->rmagerr,
                  pApassImage->imagerr);
				}
#ifdef DEBUG_REFNUMBER
				if (pApassImage->REFNumber == DEBUG_REFNUMBER) {
					printf("At REFNumber (1) index %d %lld\n",index1,pApassImage->REFNumber);
				}
#endif
				if (pUcac4Image->class == 0) {
					/* Accept the Ucac4 position as good */
					pApassImage->ra = pUcac4Image->ra;
					pApassImage->dec = pUcac4Image->dec;
					/* Calcuate a new current bin based on the new position */
					pApassExt->gsc_bin_index = GetGSCBin(pGscBin,pApassImage->ra,pApassImage->dec,&decBin,&raBin,"MatchStars");
					pApassExt->dec_bin_index = GetDecBin(pGscBin,pApassImage->dec,"matchucac4");
					/* If this is an APASS object, calculate a new position */
					if (GetREFType(pApassImage->REFNumber) == REF_TYPE_APASS) {
#ifdef DEBUG_REFNUMBER
						if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
							printf("At REFNumber (6) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
						}
#endif
#if (defined(DEBUG_WRITTEN) || defined(DEBUG_MATCH))
						tempREFNumber = pApassImage->REFNumber;
#endif /* DEBUG_WRITTEN */
						GetDASCHNumber(pApassImage->ra,pApassImage->dec,&pApassImage->REFNumber,1,REF_TYPE_APASS);
#ifdef DEBUG_REFNUMBER
						if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
							printf("At REFNumber (7) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
						}
#endif
#if (defined(DEBUG_WRITTEN) || defined(DEBUG_MATCH))
						if (tempREFNumber != pApassImage->REFNumber) {
							printf("b Oldname: %lld\n",tempREFNumber);
							printf("c Newname: %lld\n",pApassImage->REFNumber);
						}
#endif /* DEBUG_WRITTEN */
					}
				}
				/* Replace the proper motions */
				pApassImage->RaPM = pUcac4Image->RaPM;
				pApassImage->DecPM = pUcac4Image->DecPM;
				pApassImage->RaSigmaPM = pUcac4Image->RaSigmaPM;
				pApassImage->DecSigmaPM = pUcac4Image->DecSigmaPM;

				if (pUcac4Ext->matchFlag == 0) {

					++pMatchCommon->matchSequence;
					pUcac4Ext->matchFlag = 1;
					pUcac4Ext->matchSequence = pMatchCommon->matchSequence;
					pApassExt->matchFlag = 1;
					pApassExt->matchSequence = pMatchCommon->matchSequence;
#ifdef DEBUG_MATCH
					printf("First match %lld and %lld sequence %d\n",pApassImage->REFNumber,pUcac4Ext->gscImage.REFNumber,pMatchCommon->matchSequence);
#endif /* DEBUG_MATCH */
					if (pUcac4Ext->proper_motion_count > 1) {
						/* Having matched one of the ghost images, make sure that the matchFlag is set for the remaining ghost images */
						for (index2 = 0; index2 < pMatchCommon->ucac4CurCount; index2++) {
							pUcac4Ext2 = &pMatchCommon->pUcac4ExtTable[index2];
							if (pUcac4Ext->proper_motion_sequence == pUcac4Ext2->proper_motion_sequence) {
								pUcac4Ext2->matchFlag = 1;
								pUcac4Ext2->matchSequence = pMatchCommon->matchSequence;
							}
						}	
					}
				} else {
					pMatchCommon->duplicateMatchCount++;
					/* Search for the previous match and update the match count */
					curMatchCount = pUcac4Ext->matchFlag + 1;
					if (curMatchCount > pMatchCommon->maxMatchCount) {
						pMatchCommon->maxMatchCount = curMatchCount;
					}
					pApassExt->matchSequence = pUcac4Ext->matchSequence;
					for (index2 = 0; index2 < pMatchCommon->apassCurCount; index2++) {
						pApassExt2 = &pMatchCommon->pApassExtTable[index2];
						if (pApassExt2->matchSequence == pUcac4Ext->matchSequence) {
#ifdef DEBUG_REFNUMBER
							if (pApassExt2->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
								printf("At REFNumber (4) index %d %lld\n",index2,pApassExt2->gscApassImage.REFNumber);
							}
#endif
#ifdef DEBUG_MATCH
							printf("Subsequent match %lld and %lld sequence %d\n",pApassExt2->gscApassImage.REFNumber,pUcac4Ext->gscImage.REFNumber,pMatchCommon->matchSequence);
#endif /* DEBUG_MATCH */
							pApassExt2->matchFlag = curMatchCount;
							pApassExt2->matchSequence = pUcac4Ext->matchSequence;
						}
					}

				}
			}

			kd_res_free(presults);
 

		}

	}
#ifdef DEBUG_REFNUMBER
	for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
		pApassExt = &pMatchCommon->pApassExtTable[index1];
		if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
			printf("At REFNumber (9) index %d %lld sequence %d\n",index1,pApassExt->gscApassImage.REFNumber,pApassExt->matchSequence);
		}
	}
#endif
		
}



void LoadUcac4Stars(PMATCHCOMMON pMatchCommon)
{
	PBININDEX pBinIndex;
	PSTARINDEX pStarIndex;
	int readItems;
	int starCount;
	int binIndex;
	int readCount;
	int starIndex;
	PGSCIMAGE pGscImage;
	double factor;
	double extent;
	int numGhostImages;
	PGSCAPASSEXT tmpUcac4ExtTable;
	double plateepoch;
	int ghostIndex;
	PGSCAPASSEXT pUcac4Ext;
	int result;
	off_t offset;
	int index1;
#ifdef DEBUG_REFNUMBER
	for (index1 = 0; index1 < pMatchCommon->apassCurCount; index1++) {
		PGSCAPASSEXT pApassExt;
		pApassExt = &pMatchCommon->pApassExtTable[index1];
		if (pApassExt->gscApassImage.REFNumber == DEBUG_REFNUMBER) {
			printf("At REFNumber (18) index %d %lld\n",index1,pApassExt->gscApassImage.REFNumber);
		}
	}
#endif
	if ((pMatchCommon->ucac4_loaded_bin < 0) ||
			(pMatchCommon->ucac4_loaded_bin >= pGscBin->dec_bins)) {
		return; /* Nothing to load here */
	}
	pBinIndex = &pGscBin->pBinMasterIndex[pMatchCommon->ucac4_loaded_bin];

	if (pBinIndex->numBins > pGscBin->max_ra_bins) {
		printf("ERROR numBins exceeds max_ra_bins\n");
		exit(-1);
	}
	Seek(pMatchCommon->ucac4IndexHandle,pBinIndex->startBin * sizeof(STARINDEX),SEEK_SET);
	readItems = Read(pMatchCommon->ucac4IndexHandle,pMatchCommon->pIndexTable,sizeof(STARINDEX),pBinIndex->numBins);
	if (readItems != pBinIndex->numBins) {
		printf("ERROR reading the ucac4 index file\n");
		exit(-1);
	}
	/* Find out how many stars we need to read */
	starCount = 0;
	for (binIndex = 0; binIndex < pBinIndex->numBins; binIndex++) {
		pStarIndex = &pMatchCommon->pIndexTable[binIndex];
		starCount += pStarIndex->numStars;
	}
#if 0
	printf("%d stars for bin %d\n",starCount,pMatchCommon->ucac4_loaded_bin);
#endif
	if (starCount == 0) {
		/* Nothing to read, we are done */
		return;
	}
	if (starCount > pMatchCommon->imageTableAlloc) {
		if (pMatchCommon->pImageTable != 0) {
			free(pMatchCommon->pImageTable);
		}
		pMatchCommon->pImageTable = (PGSCIMAGE)calloc(starCount+EXTENSTION_ALLOC_INCREMENT,sizeof(GSCIMAGE));
		if (pMatchCommon->pImageTable == NULL) {
			printf("ERROR allocating image buffer of size %d\n",starCount);
		}
		pMatchCommon->imageTableAlloc = starCount+EXTENSTION_ALLOC_INCREMENT;
	}
	result = Seek(pMatchCommon->ucac4CatalogHandle,pMatchCommon->pIndexTable->offset,SEEK_SET);
	offset = Tell(pMatchCommon->ucac4CatalogHandle);
	readCount = Read(pMatchCommon->ucac4CatalogHandle,pMatchCommon->pImageTable,sizeof(GSCIMAGE),starCount);
	if (readCount != starCount) {
		printf("ERROR: readCount %d %d does not agree\n",
					 readCount,starCount);
		exit(-1);		
	}
	pMatchCommon->ucac4ReadCount += starCount;
	for (starIndex = 0; starIndex < readCount; starIndex++) {
		pGscImage = &pMatchCommon->pImageTable[starIndex];
		factor = cos(DEGREES_TO_RAD*pGscImage->dec);
		extent = (APASS_END_YEAR-APASS_START_YEAR)*sqrt(sqr(pGscImage->RaPM)+sqr(pGscImage->DecPM))/1000.0;
		if (extent < EXTENT_BIN_SIZE/2.0) {
			numGhostImages = 1; /* Use only the starting image */
		} else {
			numGhostImages = 1+((extent+2.0)/(EXTENT_BIN_SIZE));
		}
		if (numGhostImages > 1) {
			pMatchCommon->proper_motion_sequence++;
			if (numGhostImages > pMatchCommon->maxGhostImages) {
				pMatchCommon->maxGhostREFNumber = pGscImage->REFNumber;
				pMatchCommon->maxGhostImages = numGhostImages;
				pMatchCommon->maxGhostSequence = pMatchCommon->proper_motion_sequence;
				pMatchCommon->maxGhostBin = pMatchCommon->ucac4_loaded_bin;
			}
#if 0
			if (pMatchCommon->proper_motion_sequence == 67915) {
				printf("At proper_motion_sequence %d\n",pMatchCommon->proper_motion_sequence);
			}
#endif
		}

		for (ghostIndex = 0; ghostIndex < numGhostImages; ghostIndex++) {
			if ((pMatchCommon->ucac4CurCount+1) >= pMatchCommon->ucac4AllocCount) {
				/* Need to reallocate the UCAC4 extenstion list */
				pMatchCommon->ucac4AllocCount += EXTENSTION_ALLOC_INCREMENT;
				tmpUcac4ExtTable = (PGSCAPASSEXT)realloc(pMatchCommon->pUcac4ExtTable,pMatchCommon->ucac4AllocCount*sizeof(GSCAPASSEXT));
				if (tmpUcac4ExtTable == NULL) {
					printf("ERROR: failed to reallocate pMatchCommon->pUcac4ExtTable of size %d\n",pMatchCommon->ucac4CurCount);
					exit(-1);
				}
				pMatchCommon->pUcac4ExtTable = tmpUcac4ExtTable;
				tmpUcac4ExtTable = NULL;
			}
			pUcac4Ext = &pMatchCommon->pUcac4ExtTable[pMatchCommon->ucac4CurCount];
			memset(pUcac4Ext,0,sizeof(GSCAPASSEXT));
      pUcac4Ext->gscApassImage.REFNumber = pGscImage->REFNumber;
      pUcac4Ext->gscApassImage.ra = pGscImage->ra;
      pUcac4Ext->gscApassImage.dec = pGscImage->dec;
      pUcac4Ext->gscApassImage.Stdmag = pGscImage->Stdmag;
      pUcac4Ext->gscApassImage.color = pGscImage->color;
      pUcac4Ext->gscApassImage.RaPM = pGscImage->RaPM;
      pUcac4Ext->gscApassImage.DecPM = pGscImage->DecPM;
      pUcac4Ext->gscApassImage.RaSigmaPM = pGscImage->RaSigmaPM;
      pUcac4Ext->gscApassImage.DecSigmaPM = pGscImage->DecSigmaPM;
      pUcac4Ext->gscApassImage.class = pGscImage->class;
      pUcac4Ext->gscApassImage.VFlag = pGscImage->VFlag;
      pUcac4Ext->gscApassImage.MAGFlag = pGscImage->MAGFlag;
      pUcac4Ext->gscApassImage.flag = pGscImage->flag;
      pUcac4Ext->gscApassImage.gmag = 99.0;
      pUcac4Ext->gscApassImage.rmag = 99.0;
      pUcac4Ext->gscApassImage.imag = 99.0;
      pUcac4Ext->gscApassImage.vmagerr = 9.999;
      pUcac4Ext->gscApassImage.bmagerr = 9.999;
      pUcac4Ext->gscApassImage.gmagerr = 9.999;
      pUcac4Ext->gscApassImage.rmagerr = 9.999;
      pUcac4Ext->gscApassImage.imagerr = 9.999;
			if (numGhostImages > 1) {
				pUcac4Ext->proper_motion_sequence = pMatchCommon->proper_motion_sequence;
				pUcac4Ext->proper_motion_count = numGhostImages;
				pUcac4Ext->proper_motion_index = ghostIndex;
			}
			if (numGhostImages == 1) {
				plateepoch = (APASS_START_YEAR - GSC_EQUINOX);
			} else {
				if (ghostIndex == 0) {
					plateepoch = 0.0;
				} else {
					plateepoch = (APASS_START_YEAR - GSC_EQUINOX) + ((1.0*ghostIndex)*(APASS_END_YEAR-APASS_START_YEAR)/(1.0*numGhostImages-1));
				}
			}
			pUcac4Ext->ra = pGscImage->ra + (((pGscImage->RaPM *plateepoch)/(3600.0*1000.0))/factor);
			while (pUcac4Ext->ra > 360.0) {
				pUcac4Ext->ra -= 360.0;
			}
			while (pUcac4Ext->ra < 0.0) {
				pUcac4Ext->ra += 360.0;
			}

			pUcac4Ext->dec = pGscImage->dec +  ((pGscImage->DecPM * plateepoch)/(3600.0*1000.0));
			if (pUcac4Ext->dec > 90.0) {
				pUcac4Ext->dec = 180.0-pUcac4Ext->dec;
				pUcac4Ext->ra += 180.0;
				if (pUcac4Ext->ra >= 360.0) {
					pUcac4Ext->ra -= 360.0;
				}
			}
			if (pUcac4Ext->dec < -90.0) {
				pUcac4Ext->dec = -180.0-pUcac4Ext->dec;
				pUcac4Ext->ra += 180.0;
				if (pUcac4Ext->ra >= 360.0) {
					pUcac4Ext->ra -= 360.0;
				}
			}
				

			pMatchCommon->ucac4CurCount++;
#if 0
			if (numGhostImages > 1) {
				printf("%d\t%d\t%d\t%f\t%f\t%f\t%f\n",
							 pUcac4Ext->proper_motion_sequence,
							 pUcac4Ext->proper_motion_count,
							 pUcac4Ext->proper_motion_index,
							 extent,
							 plateepoch,
							 pUcac4Ext->ra,
							 pUcac4Ext->dec);
			}
#endif
		}

	}
	return;
}
int FindRecordSize(char *filename,File indexHandle)
{
  int statResult;
  struct stat filestats;
  off_t filesize;
  STARINDEX indexentry;
  PSTARINDEX pIndexEntry = &indexentry;
  int readItems;
  int iteration = 0;
  int gsc_bin_index = pGscBin->total_gsc_bins-1;
  int recordSize = -1;
  

  /* Get the file size */
  statResult = stat(filename,&filestats);
  if (statResult != 0) {
    printf("ERROR: FindRecordSize failed to stat %s\n",filename);
    exit(-1);
  }
  filesize = filestats.st_size;
  /* Now get the last index in the file */
  while (gsc_bin_index > 0) {

    Seek(indexHandle,gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
    readItems = Read(indexHandle,pIndexEntry,sizeof(STARINDEX),1);
    if (readItems != 1) {
      printf("ERROR: FindRecordSize failed to stat %s\n",filename);
      exit(-1);
    }
    if (pIndexEntry->numStars > 0) {
      recordSize = (filesize-pIndexEntry->offset)/pIndexEntry->numStars;
      return(recordSize);
    }
    iteration++;
    gsc_bin_index--;
  }
  printf("ERROR: FindRecordSize found no stars in %s\n",filename);
  exit(-1);
  return(recordSize);
  
}



int main(int argc,char *argv[])
{
	char apassdatname[MAX_BUFFER];
	char *argstr; 
	int errorFlag = 0;
	char cmdchar;
	int nvals;
	char *charPtr;
	char apassidxname[MAX_BUFFER];
	time_t startTime;
	time_t curTime;
	int mergeWriteItems;
	char typechar;
	char mergedatname[] = "/dasch/Pipeline/catalogs/mergeapass.dat";
	char mergeidxname[] = "/dasch/Pipeline/catalogs/mergeapass.idx";
	char outname[] = "/dasch/Pipeline/catalogs/matchapassucac4.db";
  char resultname[] =  "/dasch/Pipeline/catalogs/matchapass_final.db";
	int printFlag = 1;

	MATCHCOMMON matchCommon;
	PMATCHCOMMON pMatchCommon = &matchCommon;
	STARINDEX starIndex;
	PSTARINDEX pApassIndex = &starIndex;
  int recordSize;

	memset(pApassIndex,0,sizeof(STARINDEX));
	memset(pMatchCommon,0,sizeof(MATCHCOMMON));
	pMatchCommon->catalogType = CATALOG_TYPE_NONE;
	pMatchCommon->last_written_gsc_bin = -1;

#if (defined(DEBUG_START_INDEX) || defined(DEBUG_END_INDEX))
	printf("ERROR: DEBUG_START_INDEX %d or DEBUG_END_INDEX %d is defined\n",DEBUG_START_INDEX,DEBUG_END_INDEX);
#endif
	pMatchCommon->ucac4datname[0] = 0;
	apassdatname[0] = 0;

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
							pMatchCommon->catalogType = CATALOG_TYPE_APASS;
							break;
						case 'y':
							pMatchCommon->catalogType = CATALOG_TYPE_YB6;
							break;
						case 'u':
							pMatchCommon->catalogType = CATALOG_TYPE_UCAC3;
							break;
						case 'g':
							pMatchCommon->catalogType = CATALOG_TYPE_GSC;
							break;
						case 'e':
							pMatchCommon->catalogType = CATALOG_TYPE_EXPERIMENTAL;
							break;
						case 'k':
							pMatchCommon->catalogType = CATALOG_TYPE_KEPLER;
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
					nvals = sscanf(*++argv,"%lf",&pMatchCommon->tolerance);
					if (nvals != 1) {
						fprintf(stderr,"ERROR: Can not decode the tolerance\n");
						errorFlag = 1;
					} else {
						pMatchCommon->tolerance = pMatchCommon->tolerance/3600.;
					}
					break;


				case 'f': /* ucac4 catalog name */
				case 'F':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(pMatchCommon->ucac4datname,*++argv,MAX_BUFFER-2);
						if (strlen(*argv) >= MAX_BUFFER-2) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
						}
					}
					break;

				case 's': /* apass catalog name */
				case 'S':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(apassdatname,*++argv,MAX_BUFFER-2);
						if (strlen(*argv) >= MAX_BUFFER-2) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
						}
					}
					break;

				case 'w':
				case 'W':
					pMatchCommon->writeMergeCatalog = 1;
					break;

				case 'm':
				case 'M':
					pMatchCommon->writeMatchStarbase = 1;
					break;

				case 'k':
				case 'K':
					pMatchCommon->includeKeplerColor = 1;
					break;

				case 'p':
				case 'P':
					pMatchCommon->includeProperMotions = 1;
					break;

				case 'l':
				case 'L':
					pMatchCommon->includePositions = 1;
					break;


				default:
					printf("* illegal command -%c-",cmdchar);
					errorFlag = 1;
					break;
				}
        
			}
		}
	}
	if (pMatchCommon->catalogType == CATALOG_TYPE_NONE) {
		printf("ERROR: missing catalog type\n");
		errorFlag = 1;
	}
	if (pMatchCommon->tolerance == 0) {
		printf("ERROR: tolerance is not specified\n");
		errorFlag = 1;
	}
	if (pMatchCommon->ucac4datname[0] == 0) {
		printf("ERROR: ucac4 catalog name is not specified\n");
		errorFlag  = 1;
	} else {
		strcpy(pMatchCommon->ucac4idxname,pMatchCommon->ucac4datname);
		charPtr = strstr(pMatchCommon->ucac4idxname,".dat");
		if (charPtr == NULL) {
			printf("ERROR: ucac4 catalog name does not end in '.dat'\n");
			errorFlag = 1;
		} else {
			*charPtr = 0;
			strcat(pMatchCommon->ucac4idxname,".idx");
		}
	}
	if (apassdatname[0] == 0) {
		printf("ERROR: apass catalog name is not specified\n");
		errorFlag  = 1;
	} else {
		strcpy(apassidxname,apassdatname);
		charPtr = strstr(apassidxname,".dat");
		if (charPtr == NULL) {
			printf("ERROR: apass catalog name does not end in '.dat'\n");
			errorFlag = 1;
		} else {
			*charPtr = 0;
			strcat(apassidxname,".idx");
		}
	}
    


	if (errorFlag) {
		printf("Usage: matchucac4 \n");
		printf("       where -c <catalogtype> where catalogtype = 'a' for APASS\n");
		printf("                                                  'y; for YB6\n");
		printf("                                                  'u' for UCAC3\n");
		printf("                                                  'g' for GSC2.3.2\n");
		printf("                                                  'k' for KEPLER\n");
		printf("                                                  'e' for experimental (APASS w/GSC color)\n");
		printf("             -k for INCLUDE_KEPLER_COLOR \n");
		printf("             -p for include proper motions \n");
		printf("             -l for include positions from ucac4 unless class is nonzero \n");
		printf("             -m for WRITE_MATCH_STARBASE (/dasch/Pipeline/catalogs/matchgsckepler.db)\n");
		printf("             -w for WRITE_MERGE_CATALOG (/dasch/Pipeline/catalogs/merge.dat)  \n");
		printf("             -t <tolerance> where tolerance is the match tolerance in arcsec\n");    
		printf("             -f for ucac4 catalog name  (G), eg /dasch/Pipeline/catalogs/ucac3.dat (GSC catalog)\n");
		printf("             -s for apass catalog name (K), eg /dasch/Pipeline/catalogs/apass_temp.dat (Kepler catalog)\n");
		printf("             Note: the apass (kepler) catalog determines the output count and merged index\n");

		return(-1);
	}


	time(&startTime);
  printf("matchucac4 of %s %s sizeof(STARINDEX) %d sizeof(GSCIMAGE) %d sizeof(GSCAPASSEXT) %d\n",__DATE__,__TIME__,sizeof(STARINDEX),sizeof(GSCIMAGE),sizeof(GSCAPASSEXT));


	pMatchCommon->apassCatalogHandle = Open(apassdatname,"r");
	if (pMatchCommon->apassCatalogHandle == NULL) {
		printf("ERROR Could not open catalog file %s\n",apassdatname);
		errorFlag = 1;
	} 
	pMatchCommon->apassIndexHandle = Open(apassidxname,"r");
	if (pMatchCommon->apassIndexHandle == NULL) {
		printf("ERROR Could not open catalog index file %s\n",apassidxname);
		errorFlag = 1;
	} 

  recordSize = FindRecordSize(apassdatname,pMatchCommon->apassIndexHandle);
  if (recordSize == sizeof(GSCIMAGE)) {
    pMatchCommon->apassColorFlag = 0;
  } else if (recordSize == sizeof(GSCAPASSIMAGE)) {
    pMatchCommon->apassColorFlag = 1;
  } else {
    printf("ERROR: recordsize %d for %s is not GSCIMAGE %d or GSCAPASSIMAGE %d\n",recordSize,apassdatname,sizeof(GSCIMAGE),sizeof(GSCAPASSIMAGE));
    exit(-1);
  }


	pMatchCommon->ucac4CatalogHandle = Open(pMatchCommon->ucac4datname,"r");
	if (pMatchCommon->ucac4CatalogHandle == NULL) {
		printf("ERROR Could not open catalog file %s\n",pMatchCommon->ucac4datname);
		errorFlag = 1;
	} 
	pMatchCommon->ucac4IndexHandle = Open(pMatchCommon->ucac4idxname,"r");
	if (pMatchCommon->ucac4IndexHandle == NULL) {
		printf("ERROR Could not open catalog index file %s\n",pMatchCommon->ucac4idxname);
		errorFlag = 1;
	} 

  recordSize = FindRecordSize(pMatchCommon->ucac4datname,pMatchCommon->ucac4IndexHandle);
  if (recordSize != sizeof(GSCIMAGE)) {
    printf("ERROR: wrong record size %d %d in file %s\n",recordSize,sizeof(GSCIMAGE),pMatchCommon->ucac4datname);
    exit(-1);
  }


	if (pMatchCommon->writeMergeCatalog) {
		pMatchCommon->mergeCatalogHandle = Open(mergedatname,"w");
		if (pMatchCommon->mergeCatalogHandle == NULL) {
			printf("Could not open file %s\n",mergedatname);
			exit(-1);
		} 

		pMatchCommon->mergeIndexHandle = Open(mergeidxname,"w");
		if (pMatchCommon->mergeIndexHandle == NULL) {
			printf("Could not open file %s\n",mergeidxname);
			exit(-1);
		} 

	}


	if (pMatchCommon->writeMatchStarbase) {
		pMatchCommon->outHandle = fopen(outname,"wt");
		if (pMatchCommon->outHandle == NULL) {
			printf("ERROR count not open output file %s\n",outname);
      exit(-1);
		}
		fprintf(pMatchCommon->outHandle,"KREF\tKra\tKdec\tKStdmag\tKcolor\tKRaPM\tKDecPM\tKclass\tKVFlag\tKMAGFlag\tGREF\tGra\tGdec\tGStdmag\tGcolor\tGRaPM\tGDecPM\tGclass\tGVFlag\tGMAGFlag\tdegDrad\tflag\tgmag\trmag\timag\tvmagerr\tbmagerr\tgmagerr\trmagerr\timagerr\n");
		fprintf(pMatchCommon->outHandle,"----\t---\t----\t-------\t------\t-----\t------\t------\t------\t--------\t----\t---\t----\t-------\t------\t-----\t------\t------\t------\t--------\t-------\t----\t----\t----\t----\t-------\t-------\t-------\t-------\t-------\n");

	}

  pMatchCommon->resultHandle = fopen(resultname,"wt");
  if (pMatchCommon->resultHandle == NULL) {
    printf("ERROR count not open result file %s\n",resultname);
    exit(-1);
  }

  fprintf(pMatchCommon->resultHandle,"REF\tra\tdec\tStdmag\tcolor\tRaPM\tDecPM\tRaSigmaPM\tDecSigmaPM\tclass\tVFlag\tMAGFlag\tflag\tgmag\trmag\timag\tvmagerr\tbmagerr\tgmagerr\trmagerr\timagerr\n");
  fprintf(pMatchCommon->resultHandle,"---\t--\t---\t------\t-----\t----\t-----\t---------\t----------\t-----\t-----\t-------\t----\t----\t----\t----\t-------\t-------\t-------\t-------\t-------\n");


	if (errorFlag) {
		exit(-1);
	}

	printf("matchucac4 of %s %s, GSCIMAGE size %d\n matching %s and %s \n output %s total bins %d tolerance %f degrees\n",
				 __DATE__,__TIME__,sizeof(GSCIMAGE),pMatchCommon->ucac4datname,apassdatname,outname,pGscBin->total_gsc_bins,pMatchCommon->tolerance);
	InitBinIndex(pGscBin);

#ifdef DEBUG_START_INDEX
	pMatchCommon->ucac4_loaded_bin = DEBUG_START_INDEX;
#endif

	pMatchCommon->cur_dec_bin = pMatchCommon->ucac4_loaded_bin-DEC_BIN_OFFSET;
	pMatchCommon->cur_write_bin = pMatchCommon->ucac4_loaded_bin -(1+ (2*DEC_BIN_OFFSET));
	pMatchCommon->pIndexTable = (PSTARINDEX)calloc(pGscBin->max_ra_bins,sizeof(STARINDEX));
	if (pMatchCommon->pIndexTable == NULL) {
		printf("ERROR: failed to allocate the star index table of size %d\n",pGscBin->max_ra_bins);
		exit(-1);
	}
	while (1) {
		/* Load in the next declination range of UCAC4 stars */
		time(&curTime);
		curTime -= startTime;
		if ((curTime % 60) == 0) {
			if (printFlag == 1) {
				printf(" At %5.1f at %3d seconds\n",(pGscBin->bin_size*pMatchCommon->cur_dec_bin)-90.0,curTime);
				fflush(stdout);
				printFlag = 0;
			}
		} else {
			printFlag = 1;
		}

		LoadUcac4Stars(pMatchCommon);
		LoadAPASSStars(pMatchCommon);
		MatchStars(pMatchCommon);
		WriteMatchedStars(pMatchCommon);
		FlushUcac4Stars(pMatchCommon);
		pMatchCommon->ucac4_loaded_bin++;
		pMatchCommon->cur_dec_bin++;
		pMatchCommon->cur_write_bin++;
#ifdef DEBUG_END_INDEX
		if (pMatchCommon->cur_write_bin > DEBUG_END_INDEX) {
			break;
		}
#endif
		if (pMatchCommon->cur_write_bin > pGscBin->dec_bins) {
			break;
		}
	}
#ifdef DEBUG_END_INDEX
	time(&curTime);
	curTime -= startTime;
	printf("Writing remainder of index at %d seconds\n",curTime);
#endif /* DEBUG_END_INDEX */
	while ((pMatchCommon->last_written_gsc_bin+1) <  pGscBin->total_gsc_bins) {
		pApassIndex->offset = pMatchCommon->mergeCatalogOffset;
		pApassIndex->binNumber = pMatchCommon->last_written_gsc_bin+1;
		pApassIndex->numStars = 0;
		mergeWriteItems = Write(pMatchCommon->mergeIndexHandle,pApassIndex,sizeof(STARINDEX),1);
		if (mergeWriteItems != 1) {
			printf("ERROR writing the merge index file\n");
			exit(-1);
		}
		pMatchCommon->last_written_gsc_bin++;
	}
	time(&curTime);
	curTime -= startTime;

	printf("Finished %d matches with %d duplicates for %d written stars at %d seconds\n",pMatchCommon->matchCount,pMatchCommon->duplicateMatchCount,pMatchCommon->mergeWriteCount,curTime);
	printf("maxdupCount %d maxdupRA %f maxdupDec %f\n",pMatchCommon->maxdupCount,pMatchCommon->maxdupRA,pMatchCommon->maxdupDec);
	printf("maxGhostImages %d maxGhostSequence %d maxGhostBin %d maxGhostREFNumber %lld\nucac4 Alloc %d Read %d Flush %d\napass Alloc %d Read %d Flush %d\nproper_motion_sequence %d matchSequence %d maxMatchCount %d\n",
				 pMatchCommon->maxGhostImages,
				 pMatchCommon->maxGhostSequence,
				 pMatchCommon->maxGhostBin,
				 pMatchCommon->maxGhostREFNumber,
				 pMatchCommon->ucac4AllocCount,
				 pMatchCommon->ucac4ReadCount,
				 pMatchCommon->ucac4FlushCount,
				 pMatchCommon->apassAllocCount,
				 pMatchCommon->apassReadCount,
				 pMatchCommon->apassFlushCount,
				 pMatchCommon->proper_motion_sequence,
				 pMatchCommon->matchSequence,
				 pMatchCommon->maxMatchCount);

	kd_free(pMatchCommon->ptree);
	Close(pMatchCommon->ucac4CatalogHandle);
	Close(pMatchCommon->ucac4IndexHandle);
	Close(pMatchCommon->apassCatalogHandle);
	Close(pMatchCommon->apassIndexHandle);
	if (pMatchCommon->writeMatchStarbase) {
		fclose(pMatchCommon->outHandle);
	}
  if (pMatchCommon->resultHandle != NULL) {
    fclose(pMatchCommon->resultHandle);
  }

	if (pMatchCommon->writeMergeCatalog) {
		Close(pMatchCommon->mergeCatalogHandle);
		Close(pMatchCommon->mergeIndexHandle);
	}
	if (pMatchCommon->pImageTable != NULL) {
		free(pMatchCommon->pImageTable);
	}
	if (pMatchCommon->pImageTableApass != NULL) {
		free(pMatchCommon->pImageTableApass);
	}
	if (pMatchCommon->pIndexTable != NULL) {
		free(pMatchCommon->pIndexTable);
	}
	if (pMatchCommon->pApassExtTable != NULL) {
		free(pMatchCommon->pApassExtTable);
	}
	if (pMatchCommon->pUcac4ExtTable != NULL) {
		free(pMatchCommon->pUcac4ExtTable);
	}
	

	return(0);
}
