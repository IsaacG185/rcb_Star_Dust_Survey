// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* searchgsc.c
 *
 *   Search the gsc catalog for a given object in a given gsc_bin_index or its neighbors
 *  
 * cc -ggdb -c -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include searchgsc.c -o searchgsc.o
ar csr pipelineutils.a pipelineutils.o photometryutils.o searchgsc.o
 *   
 * Jun 19, 2009 Edward J. Los - Initial Version
 * Aug 12, 2013 Edward J. Los   Reformat GSCIMAGE, eliminating GSCIMAGEX
 */


#include <math.h>
#include <errno.h>
#include "table.h"
#include "pipelineutils.h"
#include "searchgsc.h"
#include "mysql.h"

extern int binSearchCount;
extern int cacheHits;
extern int cacheMisses;

PGSCIMAGE SearchStarSingleBin(PGSCBIN pGscBin,PGSC_CACHE_ENTRY pGscCache,File indexHandle,File catalogHandle,int gsc_bin_index,int *trueGscBinIndex,char *REF) 
{
  int cacheIndex;
  int  static curSequence = 1;
  PGSC_CACHE_ENTRY curCacheEntry;
  PGSC_CACHE_ENTRY nextCacheEntry;
  PSTARINDEX pCurStarIndex;
  int readItems;
  PGSCIMAGE pCurGscImage;
  int gscIndex;
	int refType;
	long long REFNumber;

	if (GetREFNumber(REF,&REFNumber,&refType,0,0) != 0) {
		printf("ERROR: SearchStarSingleBin invalid REF %s at gsc_bin_index %d \n",REF,gsc_bin_index);
		exit(-1);
	}



  *trueGscBinIndex = -1;

  binSearchCount++;
  for (cacheIndex = 0; cacheIndex < GSC_CACHE_ENTRIES; cacheIndex++) {
    curCacheEntry = &pGscCache[cacheIndex];
    if (curCacheEntry->gsc_bin_index == gsc_bin_index) {
      cacheHits++;
      break;
    }

  }
  if (cacheIndex == GSC_CACHE_ENTRIES) {
    cacheMisses++;
    nextCacheEntry = &pGscCache[0];
    for (cacheIndex = 1; cacheIndex < GSC_CACHE_ENTRIES; cacheIndex++) {
      curCacheEntry = &pGscCache[cacheIndex];
      if (curCacheEntry->sequence < nextCacheEntry->sequence) {
        nextCacheEntry = curCacheEntry;
      }
    }
    curCacheEntry = nextCacheEntry;
  }
  pCurStarIndex = &curCacheEntry->gscStarIndex;
  if (curCacheEntry->sequence != curSequence) {
    curCacheEntry->sequence = ++curSequence;
    if (curSequence > 0x7fffffff) {
      fprintf(stderr,"ERROR: curSequence is wrapping in SearchStarSingleBin\n");
      exit(-1);
    }
  }
  if (curCacheEntry->gsc_bin_index != gsc_bin_index) {
    /* We need to read in these stars */
    if (gsc_bin_index >= pGscBin->total_gsc_bins) {
      printf("ERROR: gsc_bin_index %d exceeds total %d in SearchStarSingleBin\n",gsc_bin_index,pGscBin->total_gsc_bins);
      curCacheEntry->gsc_bin_index = -1;
      curCacheEntry->sequence = 0;
      memset(pCurStarIndex,0,sizeof(STARINDEX));
    } else {
      curCacheEntry->gsc_bin_index = gsc_bin_index;
      Seek(indexHandle,curCacheEntry->gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
      readItems = Read(indexHandle,pCurStarIndex,sizeof(STARINDEX),1);
      if (readItems != 1) {
        fprintf(stderr,"ERROR reading star index file errno: %d %s REF %s gsc_bin_index %d\n",
                errno,
                strerror(errno),
                REF,
                gsc_bin_index);
        curCacheEntry->gsc_bin_index = -1;
        curCacheEntry->sequence = 0;
        memset(pCurStarIndex,0,sizeof(STARINDEX));
      } else {
          
        if (pCurStarIndex->numStars > curCacheEntry->gscImageAlloc) {
          curCacheEntry->gscImageAlloc = pCurStarIndex->numStars+1000;
          pCurGscImage = realloc(curCacheEntry->pGscImageTable,curCacheEntry->gscImageAlloc*sizeof(GSCIMAGE));
          if (pCurGscImage == NULL) {
            fprintf(stderr,"ERROR: failed to realloc curCacheEntry->pGscImageTable of size %d\n",curCacheEntry->gscImageAlloc*sizeof(GSCIMAGE));
            exit(-1);
          }
          curCacheEntry->pGscImageTable = pCurGscImage;
          pCurGscImage = NULL;
        }
        Seek(catalogHandle,pCurStarIndex->offset,SEEK_SET);
        readItems = Read(catalogHandle,curCacheEntry->pGscImageTable,sizeof(GSCIMAGE),pCurStarIndex->numStars);
        if (readItems != pCurStarIndex->numStars) {
          fprintf(stderr,"ERROR reading gsc catalog file\n");
          exit(-1);
        }
      }
  

    }

  }
  

  for (gscIndex = 0; gscIndex < pCurStarIndex->numStars; gscIndex++) {
    pCurGscImage = &curCacheEntry->pGscImageTable[gscIndex];
    if (pCurGscImage->REFNumber == REFNumber) {
      *trueGscBinIndex = gsc_bin_index;
      
      return(pCurGscImage);
    }
  }
  return(NULL);

}


PGSCIMAGE SearchStar(PGSCBIN pGscBin,PGSC_CACHE_ENTRY pGscCache,File indexHandle,File catalogHandle,int gsc_bin_index,int *trueGscBinIndex,char *REF) 
{
  PBININDEX pBinIndex;
  PBININDEX pOrigBinIndex;
  int origRaBin;
  int origDecBin;
  double raval;
  double decval;
  int raBin;
  int decBin;
  int next_bin_index;
  int center_bin_index;
  PGSCIMAGE pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,gsc_bin_index,trueGscBinIndex,REF);

#if 0
  if (strcmp(REF,"N2312232541") == 0) {
    printf("At %s\n",REF);      
  }
  if (strcmp(REF,"S2031000177") == 0) {
    printf("At %s\n",REF);      
  }
#endif


  if (pGscImage != NULL) {
    return(pGscImage);
  }
  /* We need to expand our search to adjacent RA bins */
  pBinIndex = GetSubbins(pGscBin,gsc_bin_index,&raBin,&decBin,(char *)"SearchStar");
  if (pBinIndex == NULL) {
    fprintf(stderr,"ERROR: failed to get GSC bin index for bin %d\n",gsc_bin_index);
    return(NULL);
  }
  pOrigBinIndex = pBinIndex;
  origRaBin = raBin;
  origDecBin = decBin;
  /* Search the next highest bin */
  if (raBin == (pBinIndex->numBins-1)) {
    next_bin_index = pBinIndex->startBin;
  } else {
    next_bin_index = gsc_bin_index+1;
  }
  pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

  if (pGscImage != NULL) {
    return(pGscImage);
  }

  if (raBin == 0) {
    next_bin_index = pBinIndex->startBin + pBinIndex->numBins-1;
  } else {
    next_bin_index = gsc_bin_index-1;
  }
  pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

  if (pGscImage != NULL) {
    return(pGscImage);
  }
  
  /* Need to search in the next highest declination bin */
  
  raval = (360.0 * origRaBin)/(pOrigBinIndex->numBins) + (180.0/(pOrigBinIndex->numBins));
  decval = (origDecBin * pGscBin->bin_size) - 90.0 +  (pGscBin->bin_size/2.0);
  decval += pGscBin->bin_size;
  if (decval < 90.0) {
    center_bin_index = GetGSCBin(pGscBin,raval,decval,&decBin,&raBin,(char *)"SearchStar2");
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,center_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }
    pBinIndex = GetSubbins(pGscBin,center_bin_index,&raBin,&decBin,(char *)"SearchStar");
    if (pBinIndex == NULL) {
      fprintf(stderr,"ERROR: failed to get GSC bin index for bin %d\n",gsc_bin_index);
      return(NULL);
    }
    /* Search the next highest bin */
    if (raBin == (pBinIndex->numBins-1)) {
      next_bin_index = pBinIndex->startBin;
    } else {
      next_bin_index = center_bin_index+1;
    }
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }


    if ((raBin+2) >= pBinIndex->numBins) {
      next_bin_index = pBinIndex->startBin+raBin+2-pBinIndex->numBins;
    } else {
      next_bin_index = center_bin_index+2;
    }
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }

    if (raBin == 0) {
      next_bin_index = pBinIndex->startBin + pBinIndex->numBins-1;
    } else {
      next_bin_index = center_bin_index-1;
    }
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }


    if (raBin-2 < 0) {
      next_bin_index = pBinIndex->startBin + pBinIndex->numBins+raBin-2;
    } else {
      next_bin_index = center_bin_index-2;
    }
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }

  }


  /* Need to search in the next lowest declination bin */
  
  raval = (360.0 * origRaBin)/(pOrigBinIndex->numBins) + (180.0/(pOrigBinIndex->numBins));
  decval = (origDecBin * pGscBin->bin_size) - 90.0 +  (pGscBin->bin_size/2.0);
  decval -= pGscBin->bin_size;
  if (decval > -90.0) {
    center_bin_index = GetGSCBin(pGscBin,raval,decval,&decBin,&raBin,(char *)"SearchStar2");
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,center_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }
    pBinIndex = GetSubbins(pGscBin,center_bin_index,&raBin,&decBin,(char *)"SearchStar");
    if (pBinIndex == NULL) {
      fprintf(stderr,"ERROR: failed to get GSC bin index for bin %d\n",gsc_bin_index);
      return(NULL);
    }
    /* Search the next highest bin */
    if (raBin == (pBinIndex->numBins-1)) {
      next_bin_index = pBinIndex->startBin;
    } else {
      next_bin_index = center_bin_index+1;
    }
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }

    if ((raBin+2) >= pBinIndex->numBins) {
      next_bin_index = pBinIndex->startBin+raBin+2-pBinIndex->numBins;
    } else {
      next_bin_index = center_bin_index+2;
    }
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }

    if (raBin == 0) {
      next_bin_index = pBinIndex->startBin + pBinIndex->numBins-1;
    } else {
      next_bin_index = center_bin_index-1;
    }
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }

    if (raBin-2 < 0) {
      next_bin_index = pBinIndex->startBin + pBinIndex->numBins+raBin-2;
    } else {
      next_bin_index = center_bin_index-2;
    }
    pGscImage = SearchStarSingleBin(pGscBin,pGscCache,indexHandle,catalogHandle,next_bin_index,trueGscBinIndex,REF);

    if (pGscImage != NULL) {
      return(pGscImage);
    }

  }


  return(NULL);

}
