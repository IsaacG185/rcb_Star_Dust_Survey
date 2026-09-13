// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* galaxyutils.h
 *
 *  Add  support for external variable star and galaxy catalogs and for plate limiting magnitudes
 *
 * Dec 28, 2012 Edward J. Los - Initial Version
 * Feb 23, 2013 Edward J. Los - Add PNEARESTCATALOGSTAR to support printing of additional nearby objects.
 * May  6, 2013 Edward J. Los - Split the limiting magnitudes file into a "nearby" and "plates" file.
 * May 26, 2015 Edward J. Los - Add pGalaxyResult to support the SDSS catalog
 * Jun 29, 2015 Edward J. Los - Use mmap() and munmap() to map the catalogs and avoid exhausting system resources
 * Jul  6, 2015 Edward J. Los - Use gsc bin indexed SDSS catalog
 * Aug  3, 2015 Edward J. Los - Read in the Tycho2 catalog to check transient candidates against bright stars 
 *                              Remove USE_MMAP_FUNCTION because it didn't work.  Enable USE_INDEXED_SDSS permanently
 * Feb 24, 2017 Edward J. Los - Add limitingDataHandle and limitingIndexHandle to avoid conflict
 *                              Add sanity checks to prevent double initilization of galaxy common
 * Apr  4, 2018 Edward J. Los - Add SearchLimiting
 */

#ifndef _GALAXYUTILS_H
#define _GALAXYUTILS_H 1

#include "pipelineutils.h" // STARINDEX, etc.

/* Support in formatvsx and update_limiting to find nearby galaxies and stars */

  /* Of 939,546 PGC Galaxies with radii 
   *     10 > 960 "
   *     20 > 600 "
   *    100 > 240 "
   *   1000 >  90 "
   * Use a search radius of 600"
   */
#define GALAXY_SEARCH_RADIUS (600.0/3600.0)


#define MAX_GALAXY_NAME 32
#define MAX_GALAXY_TYPE 31
#define NEAREST_MAG_DIFFERENCE 1.0

#define GALAXY_DEFAULT_RADIUS 6.0 /* Default galaxy radius in arcsec */
#define GALAXY_FLAG 0x786c6167 /* "galx" */
#define GALAXY_VERSION  1
#define GALAXY_INITIALIZATION_FLAG 0x696E6967 /* gini */
typedef struct _galaxyrec {
  int galaxyflag;
  int galaxyversion;
  double ra;
  double dec;
  double catalogmag;
  double radius; /* Radius in arcsec */
  char catalogname[MAX_GALAXY_NAME];
  char galaxytype[MAX_GALAXY_TYPE];
  char variableFlag;
} GALAXYREC,*PGALAXYREC;

/* Support for plate limiting magnitudes */

#define PLATE_FLAG 0x74616c70 /* "plat" */
#define PLATE_VERSION  1
typedef struct _platelimitingrec {
  int galaxyflag;
  int galaxyversion;
  double limiting_mag_local;
  double geoJulianDate;
  int seriesId;
  int plateNumber;
  int mosaicNumber;
  int solutionNumber;
  int versionId;
	int unused;
} PLATELIMITINGREC,*PPLATELIMITINGREC;


/* Common storage area for galaxy retrieval routines */
typedef struct _galaxycommon
{
  void * nearbyDataHandle; /* (File) */
  void * nearbyIndexHandle;  /* (File) */
  void * platesDataHandle; /* (File) */
  void * platesIndexHandle;  /* (File) */
  void * limitingDataHandle; /* (File) */
  void * limitingIndexHandle;  /* (File) */
  off_t galaxy_count;
  off_t sdss_count;
  off_t tycho2_count;
  unsigned long initialized;
  int inputBufferAlloc;
  int nearbyBufferAlloc;
  int galaxyCount;
  int plateCount;
	int nearbyEnabled;
  int enableTransientSearch;
  /* The following accesses an indexed SDSS catalog */
  int sdssCatalogFD;
  int sdssIndexFD;
  int sdssCatalogAlloc;
  PGALAXYREC sdssBuffer;
  int sdss_bin_index;  /* currently cached gsc_bin_index */
  STARINDEX cursdssStarIndex;
  int sdssGscBinReadCount;
  int sdssGscCacheHitCount;
  /* The following accesses an indexec Tycho2 catalog */
  int tycho2CatalogFD;
  int tycho2IndexFD;
  int tycho2CatalogAlloc;
  PGSCIMAGE tycho2Buffer;
  int tycho2_bin_index;  /* currently cached gsc_bin_index */
  STARINDEX curtycho2StarIndex;
  int tycho2GscBinReadCount;
  int tycho2GscCacheHitCount;
  /* End of Tycho2 access */

  char *inputBuffer;
  char *nearbyBuffer;
  void *ptree;
	PGSCBIN pLimitingGscBin;
	PGSCBIN pNearbyGscBin;
  PGALAXYREC galaxybuffer;
  PGALAXYREC galaxyInputBuffer;
  PPLATELIMITINGREC plateLimitingBuffer;
} GALAXYCOMMON,*PGALAXYCOMMON;


/* Supporting return data for FindNearestGalaxy */
typedef struct _galaxyresult
{
  int variableFoundFlag;
  int galaxyFoundFlag;
  double minVariableDistance; /* Distance in arcsec */
  double minGalaxyDistance; /* Distance in arcsec */
  GALAXYREC galaxyRecVariable;
  GALAXYREC galaxyRecGalaxy;
 } GALAXYRESULT,*PGALAXYRESULT;

/* Input table for FindNearestGalaxy */
typedef struct _nearestCatalogStar {
  char nearestREF[MAX_REF]; /* Nearest matched REFNumber */
  double nearestREFarcsec; /* Nearest match REFNumber distance */
  double nearestREFmag; /* Nearest match REFNumber magnitude */
  int nearestREFflag;   /*  -1: No position; 0: all stars; +1: good stars */
  char nearestREFGoodMag[MAX_REF]; /* Nearest matched REFNumber */
  double nearestREFarcsecGoodMag; /* Nearest match REFNumber distance for good magnitudes */
  double nearestREFmagGoodMag; /* Nearest match REFNumber magnitude for good magnitudes */
  int nearestREFflagGoodMag;   /* -1: No position; 0: all stars; +1: good stars */
} NEARESTCATALOGSTAR,*PNEARESTCATALOGSTAR;

int OpenGalaxyFiles(PGALAXYCOMMON pGalaxyCommon, char *qualifier, char *prefix, char *photfilebase);
int LoadGalaxyTable(PGALAXYCOMMON pGalaxyCommon,double ra,double dec);
void FindNearestGalaxy(PGALAXYCOMMON pGalaxyCommon,double ra,double dec,PNEARESTCATALOGSTAR pNearestCatalogStar,char *nearbyObjects,PGALAXYRESULT pGalaxyResult);
void FreeGalaxyCommon(PGALAXYCOMMON pGalaxyCommon,int verbose);
int PopulateGalaxyTree(PGALAXYCOMMON pGalaxyCommon,char *galaxy_name,int enableTransientSearch);
void InitNearestCatalogStar(PNEARESTCATALOGSTAR pNearestCatalogStar);
void *GalaxyOpen(char *name,char *mode);
void DumpLimitingTable(PGALAXYCOMMON pGalaxyCommon);
void InitGalaxyCommon(PGALAXYCOMMON pGalaxyCommon);
int FindNearestTycho2Star(PGALAXYCOMMON pGalaxyCommon,double ra,double dec,double tycho2mag,double tycho2arcsec);
void SearchLimiting(PGALAXYCOMMON  pGalaxyCommon,int printLimit,char* coveragefile);

#endif /* _GALAXYUTILS_H */
