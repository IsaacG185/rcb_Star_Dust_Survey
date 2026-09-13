// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* searchgsc.h
 *
 * Jun 19, 2009 Edward J. Los - Initial Version
 */

#ifndef _SEARCHGSC_H
#define _SEARCHGSC_H 1
#define GSC_CACHE_ENTRIES 12
typedef struct _gsc_cache_entry {
  int gscImageAlloc;
  int gsc_bin_index;
  int sequence;
  PGSCIMAGE pGscImageTable;
  STARINDEX gscStarIndex;
} GSC_CACHE_ENTRY,*PGSC_CACHE_ENTRY;

PGSCIMAGE SearchStarSingleBin(PGSCBIN pGscBin,PGSC_CACHE_ENTRY pGscCache,File indexHandle,File catalogHandle,int gsc_bin_index,int *trueGscBinIndex,char *REF) ;
PGSCIMAGE SearchStar(PGSCBIN pGscBin,PGSC_CACHE_ENTRY pGscCache,File indexHandle,File catalogHandle,int gsc_bin_index,int *trueGscBinIndex,char *REF);
#endif /* _SEARCHGSC_H */
