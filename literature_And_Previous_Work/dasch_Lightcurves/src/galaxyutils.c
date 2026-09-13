// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#include <table.h>
#include <math.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <time.h>

#include "mysql.h"

#include "kdtree.h" /* Astrometry.Net */

#include "pipelineutils.h"
#include "galaxyutils.h"
#include "daschunistd.h"

#define MAX_BUFFER 500
#define GSC_BIN_INCREMENT 1000

extern GSCBIN gscBin01;
extern GSCBIN gscBin02;
extern GSCBIN gscBin04;
extern GSCBIN gscBin08;
extern GSCBIN gscBin16;
extern GSCBIN gscBin32;
extern GSCBIN gscBin64;


void
InitGalaxyCommon(PGALAXYCOMMON pGalaxyCommon)
{
  if (pGalaxyCommon->initialized == GALAXY_INITIALIZATION_FLAG) {
    printf("ERROR: InitGalaxyCommon already initialized\n");
    exit(-1);
  } else if (pGalaxyCommon->initialized != 0) {
    printf("ERROR: InitGalaxyCommon: corrupt or not pre-zeroed\n");
    exit(-1);
  }

  memset(pGalaxyCommon, 0, sizeof(GALAXYCOMMON));

  pGalaxyCommon->sdssCatalogFD = -1;
  pGalaxyCommon->sdssIndexFD = -1;
  pGalaxyCommon->sdss_bin_index = -1;
  pGalaxyCommon->tycho2CatalogFD = -1;
  pGalaxyCommon->tycho2IndexFD = -1;
  pGalaxyCommon->tycho2_bin_index = -1;
  pGalaxyCommon->initialized = GALAXY_INITIALIZATION_FLAG;
}


int
LoadGalaxyTable(PGALAXYCOMMON pGalaxyCommon,double ra,double dec)
{
  int gsc_bin_index;
  int decBin;
  int raBin;
  int readItems;
  STARINDEX curStarIndex;
  PSTARINDEX pCurStarIndex = &curStarIndex;
  PGALAXYREC pGalaxyRec = NULL;
  PPLATELIMITINGREC pPlateRec;
  int byteCount;

	if (pGalaxyCommon->pLimitingGscBin == NULL) {
		printf("ERROR: LoadGalaxyTable called without setting pLimitingGscBin\n");
		exit(-1);
	}

  gsc_bin_index = GetGSCBin(pGalaxyCommon->pLimitingGscBin,ra,dec,&decBin,&raBin,"LoadGalaxyTable");
  if (pGalaxyCommon->limitingIndexHandle == NULL || pGalaxyCommon->limitingDataHandle == NULL) {
    return 0;
  }

  Seek((File)pGalaxyCommon->limitingIndexHandle, gsc_bin_index * sizeof(STARINDEX), SEEK_SET);
  readItems = Read((File)pGalaxyCommon->limitingIndexHandle,pCurStarIndex,sizeof(STARINDEX),1);
  if (readItems != 1) {
    printf("Error: LoadGalaxyTable can not read bin %d\n",gsc_bin_index);
    exit(-1);
  }

  /* Here "numStars" really means "numBytes" */
  if (pCurStarIndex->numStars > pGalaxyCommon->inputBufferAlloc) {
    if (pGalaxyCommon->inputBuffer != NULL) {
      free(pGalaxyCommon->inputBuffer);
    }

    pGalaxyCommon->inputBuffer = (char *)calloc(pCurStarIndex->numStars,sizeof(char));

    if (pGalaxyCommon->inputBuffer == NULL) {
      printf("ERROR: failed to allocate the galaxy input buffer of size %d\n",pCurStarIndex->numStars);
      exit(-1);
    }

    pGalaxyCommon->inputBufferAlloc = pCurStarIndex->numStars;
  }

  Seek((File)pGalaxyCommon->limitingDataHandle,pCurStarIndex->offset,SEEK_SET);
  readItems = Read((File)pGalaxyCommon->limitingDataHandle,pGalaxyCommon->inputBuffer,sizeof(char),pCurStarIndex->numStars);

  if (readItems != pCurStarIndex->numStars) {
    printf("ERROR reading gsc catalog file\n");
    exit(-1);
  }

  byteCount = 0;
  pGalaxyCommon->galaxyCount = 0;
  pGalaxyCommon->galaxyInputBuffer = (PGALAXYREC)pGalaxyCommon->inputBuffer;

  while(byteCount < readItems) {
    pGalaxyRec = &pGalaxyCommon->galaxyInputBuffer[pGalaxyCommon->galaxyCount];

    if (pGalaxyRec->galaxyflag == GALAXY_FLAG &&
        pGalaxyRec->galaxyversion == GALAXY_VERSION &&
        (readItems-byteCount) >= sizeof(GALAXYREC)
    ) {
			pGalaxyCommon->galaxyCount++;
			byteCount += sizeof(GALAXYREC);
    } else {
      break;
    }
  }

  pGalaxyCommon->plateCount = 0;
  pGalaxyCommon->plateLimitingBuffer = (PPLATELIMITINGREC)pGalaxyRec;

  while(byteCount < readItems) {
    pPlateRec = &pGalaxyCommon->plateLimitingBuffer[pGalaxyCommon->plateCount];

    if (pPlateRec->galaxyflag == PLATE_FLAG &&
        pPlateRec->galaxyversion == PLATE_VERSION &&
        (readItems-byteCount) >= sizeof(PLATELIMITINGREC)
    ) {
			pGalaxyCommon->plateCount++;
			byteCount += sizeof(PLATELIMITINGREC);
    } else {
      break;
    }
  }

  if (byteCount != readItems) {
    printf(
      "ERROR parsing the galaxy catalog got %d actual %d PLATELIMITINGREC %zu GALAXYREC %zu\n",
      byteCount,
      readItems,
      sizeof(PLATELIMITINGREC),
      sizeof(GALAXYREC)
    );
    return 0;
  }

	if (pGalaxyCommon->nearbyEnabled) {
		/* There is an update list of galaxies available, discard the old list and read the new one in */
		if (pGalaxyCommon->pNearbyGscBin == NULL) {
			printf("ERROR: LoadGalaxyTable called without setting pNearbyGscBin\n");
			exit(-1);
		}

		gsc_bin_index = GetGSCBin(pGalaxyCommon->pNearbyGscBin,ra,dec,&decBin,&raBin,"LoadGalaxyTable");
		if ((pGalaxyCommon->nearbyIndexHandle == NULL) ||
				(pGalaxyCommon->nearbyDataHandle == NULL)) {
			return 0;
		}

		Seek((File)pGalaxyCommon->nearbyIndexHandle,gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
		readItems = Read((File)pGalaxyCommon->nearbyIndexHandle,pCurStarIndex,sizeof(STARINDEX),1);
		if (readItems != 1) {
			printf("Error: LoadGalaxyTable can not read bin %d\n",gsc_bin_index);
			exit(-1);
		}

		/* Here "numStars" really means "numBytes" */
		if (pCurStarIndex->numStars > pGalaxyCommon->nearbyBufferAlloc) {
			if (pGalaxyCommon->nearbyBuffer != NULL) {
				free(pGalaxyCommon->nearbyBuffer);
			}

			pGalaxyCommon->nearbyBuffer = (char *)calloc(pCurStarIndex->numStars,sizeof(char));

			if (pGalaxyCommon->nearbyBuffer == NULL) {
				printf("ERROR: failed to allocate the galaxy input buffer of size %d\n",pCurStarIndex->numStars);
				exit(-1);
			}
			pGalaxyCommon->nearbyBufferAlloc = pCurStarIndex->numStars;
		}

		Seek((File)pGalaxyCommon->nearbyDataHandle,pCurStarIndex->offset,SEEK_SET);
		readItems = Read((File)pGalaxyCommon->nearbyDataHandle,pGalaxyCommon->nearbyBuffer,sizeof(char),pCurStarIndex->numStars);

		if (readItems != pCurStarIndex->numStars) {
			printf("ERROR reading gsc catalog file\n");
			exit(-1);
		}

		byteCount = 0;
		pGalaxyCommon->galaxyCount = 0;
		pGalaxyCommon->galaxyInputBuffer = (PGALAXYREC)pGalaxyCommon->nearbyBuffer;

		while(byteCount < readItems) {
			pGalaxyRec = &pGalaxyCommon->galaxyInputBuffer[pGalaxyCommon->galaxyCount];

			if (pGalaxyRec->galaxyflag == GALAXY_FLAG &&
					pGalaxyRec->galaxyversion == GALAXY_VERSION &&
					(readItems-byteCount) >= sizeof(GALAXYREC)
      ) {
				pGalaxyCommon->galaxyCount++;
				byteCount += sizeof(GALAXYREC);
			} else {
				break;
			}
		}

		if (byteCount != readItems) {
			printf(
        "ERROR parsing the nearby galaxy catalog got %d actual %d PLATELIMITINGREC %zu GALAXYREC %zu\n",
        byteCount,
        readItems,
        sizeof(PLATELIMITINGREC),
        sizeof(GALAXYREC)
      );
			return 0;
		}

	}

  return 1;
}


/* This routine assumes that the nearbyObjects string has already been cleared or use to save the nearestREF */
void
FindNearestGalaxy(
  PGALAXYCOMMON pGalaxyCommon,
  double ra,
  double dec,
  PNEARESTCATALOGSTAR pNearestCatalogStar,
  char *nearbyObjects,
  PGALAXYRESULT pGalaxyResult
) {
  int objectIndex;
  PGALAXYREC pGalaxyRec;
  PGALAXYREC pGalaxyRecVariable = NULL;
  PGALAXYREC pGalaxyRecGalaxy = NULL;
  double minVariableDistance = BIN_EXPANSION_RADIUS * BIN_EXPANSION_RADIUS;
  double minGalaxyDistance = GALAXY_SEARCH_RADIUS * GALAXY_SEARCH_RADIUS;
  double factor;
  double curDistance;
  double radius;
  char tempBuffer[MAX_NEARBY_OBJECTS_STRING];
  double x;
  double y;
  double z;
  double cosx;
  double pt[3];
  double pos[3];
  struct kdres *presults;
  char questionable[2];
  int duplicateFlag = 0;
  int raBin;
  int decBin;
  int readBytes;
  int readItems;
  int readIndex;
  int cur_bin_index;

  if (pGalaxyResult != NULL) {
    memset(pGalaxyResult, 0, sizeof(GALAXYRESULT));
  }

  nearbyObjects[0] = 0;
  factor = cos(DEGREES_TO_RAD*dec);

  if (pGalaxyCommon->platesIndexHandle != NULL && pGalaxyCommon->platesDataHandle != NULL) {
    for (objectIndex = 0; objectIndex < pGalaxyCommon->galaxyCount; objectIndex++) {
      pGalaxyRec = &pGalaxyCommon->galaxyInputBuffer[objectIndex];
      curDistance = sqr(pGalaxyRec->dec-dec) + sqr(factor*(pGalaxyRec->ra-ra));

      if (pGalaxyRec->variableFlag) {
        if (curDistance < minVariableDistance) {
          minVariableDistance = curDistance;
          pGalaxyRecVariable = pGalaxyRec;
        }
      } else {
        if ((radius = pGalaxyRec->radius) == 0.0) {
          radius = GALAXY_DEFAULT_RADIUS;
        }

        if ((curDistance < sqr(GALAXY_SEARCH_RADIUS)) &&
            (curDistance < minGalaxyDistance) &&
            (curDistance < sqr(2.*radius/3600.0))) {
          minGalaxyDistance = curDistance;
          pGalaxyRecGalaxy = pGalaxyRec;
        }
      }
    }
  } else if (pGalaxyCommon->ptree != NULL) {
    z = sin(DEGREES_TO_RAD * dec);
    cosx = cos(DEGREES_TO_RAD * dec);
    x = cosx * cos(DEGREES_TO_RAD * ra);
    y = cosx * sin(DEGREES_TO_RAD * ra);
    pt[0] = x;
    pt[1] = y;
    pt[2] = z;

    presults = kd_nearest_range(pGalaxyCommon->ptree, pt, GALAXY_SEARCH_RADIUS * DEGREES_TO_RAD);

    while (!kd_res_end(presults)) {
      /* get the data and position of the current result item */
      pGalaxyRec = (PGALAXYREC) kd_res_item(presults, pos);
      curDistance = sqr(pGalaxyRec->dec - dec) + sqr(factor * (pGalaxyRec->ra - ra));

      if (pGalaxyRec->variableFlag) {
        if (curDistance < minVariableDistance) {
          minVariableDistance = curDistance;
          pGalaxyRecVariable = pGalaxyRec;
        }
      } else {
        if ((radius = pGalaxyRec->radius) == 0.0) {
          radius = GALAXY_DEFAULT_RADIUS;
        }

        if ((curDistance < sqr(GALAXY_SEARCH_RADIUS)) &&
            (curDistance < minGalaxyDistance) &&
            (curDistance < sqr(2.*radius/3600.0))) {
          minGalaxyDistance = curDistance;
          pGalaxyRecGalaxy = pGalaxyRec;
        }
      }

      /* go to the next entry */
      kd_res_next(presults);
    }

    kd_res_free(presults);

    if (pGalaxyCommon->sdssCatalogFD >= 0 && pGalaxyCommon->sdssIndexFD >= 0) {
      cur_bin_index = GetGSCBin(&gscBin64, ra, dec, &decBin, &raBin, "FindNearestGalaxy");

      if (cur_bin_index != pGalaxyCommon->sdss_bin_index) {
        pGalaxyCommon->sdss_bin_index = cur_bin_index;
        pGalaxyCommon->sdssGscBinReadCount++;

        lseek(pGalaxyCommon->sdssIndexFD, pGalaxyCommon->sdss_bin_index * sizeof(STARINDEX), SEEK_SET);
        readBytes = read(pGalaxyCommon->sdssIndexFD, &pGalaxyCommon->cursdssStarIndex, sizeof(STARINDEX));
        readItems = readBytes / sizeof(STARINDEX);

        if (readItems != 1) {
          printf(
            "ERROR FindNearestGalaxy reading star index file errno: %d %s gsc_bin_index %d\n",
            errno,
            strerror(errno),
            pGalaxyCommon->sdss_bin_index
          );
          exit(-1);
        }

        if (pGalaxyCommon->cursdssStarIndex.numStars > 0) {
          if (pGalaxyCommon->cursdssStarIndex.binNumber != pGalaxyCommon->sdss_bin_index) {
            printf(
              "ERROR: FindNearestGalaxy gsc_bin_index %d %d do not agree with catalog index\n",
              pGalaxyCommon->sdss_bin_index,
              pGalaxyCommon->cursdssStarIndex.binNumber
            );
            exit(-1);
          }

          if (pGalaxyCommon->cursdssStarIndex.numStars > pGalaxyCommon->sdssCatalogAlloc) {
            if (pGalaxyCommon->sdssBuffer != NULL) {
              free(pGalaxyCommon->sdssBuffer);
            }

            pGalaxyCommon->sdssCatalogAlloc = pGalaxyCommon->cursdssStarIndex.numStars;
            pGalaxyCommon->sdssBuffer = (PGALAXYREC) calloc(pGalaxyCommon->sdssCatalogAlloc, sizeof(GALAXYREC));

            if (pGalaxyCommon->sdssBuffer == NULL) {
              printf(
                "ERROR: FindNearestGalaxy failed to allocate sdssBuffer of size %d\n",
                pGalaxyCommon->sdssCatalogAlloc
              );
              exit(-1);
            }
          }

          lseek(pGalaxyCommon->sdssCatalogFD,pGalaxyCommon->cursdssStarIndex.offset,SEEK_SET);

          readBytes = read(pGalaxyCommon->sdssCatalogFD,pGalaxyCommon->sdssBuffer,sizeof(GALAXYREC) * pGalaxyCommon->cursdssStarIndex.numStars);
          readItems = readBytes/sizeof(GALAXYREC);

          if (readItems != pGalaxyCommon->cursdssStarIndex.numStars) {
            fprintf(
              stderr,
              "ERROR reading gsc catalog file at index %d\n",
              pGalaxyCommon->sdss_bin_index
            );
            exit(-1);
          }
        }
      } else {
        pGalaxyCommon->sdssGscCacheHitCount++;
      }

      if (pGalaxyCommon->cursdssStarIndex.numStars > 0) {
        for (readIndex = 0; readIndex < pGalaxyCommon->cursdssStarIndex.numStars; readIndex++) {
          pGalaxyRec = &pGalaxyCommon->sdssBuffer[readIndex];
          curDistance = sqr(pGalaxyRec->dec - dec) + sqr(factor * (pGalaxyRec->ra - ra));

          if (curDistance < sqr(BIN_EXPANSION_RADIUS) && curDistance < minGalaxyDistance) {
            minGalaxyDistance = curDistance;
            pGalaxyRecGalaxy = pGalaxyRec;
          }
        }
      }
    }
  } else if (pGalaxyCommon->limitingIndexHandle != NULL && pGalaxyCommon->limitingDataHandle != NULL) {
    for (objectIndex = 0; objectIndex < pGalaxyCommon->galaxyCount; objectIndex++) {
      pGalaxyRec = &pGalaxyCommon->galaxyInputBuffer[objectIndex];
      curDistance = sqr(pGalaxyRec->dec - dec) + sqr(factor * (pGalaxyRec->ra - ra));

      if (pGalaxyRec->variableFlag) {
        if (curDistance < minVariableDistance) {
          minVariableDistance = curDistance;
          pGalaxyRecVariable = pGalaxyRec;
        }
      } else {
        if ((radius = pGalaxyRec->radius) == 0.0) {
          radius = GALAXY_DEFAULT_RADIUS;
        }

        if ((curDistance < sqr(GALAXY_SEARCH_RADIUS)) &&
            (curDistance < minGalaxyDistance) &&
            (curDistance < sqr(2.*radius/3600.0))) {
          minGalaxyDistance = curDistance;
          pGalaxyRecGalaxy = pGalaxyRec;
        }
      }
    }
  }

  if (pGalaxyRecVariable != NULL) {
    pGalaxyRec = pGalaxyRecVariable;
    sprintf(
      tempBuffer,
      "%s:%s @ %.1f",
      pGalaxyRec->catalogname,
      pGalaxyRec->galaxytype,
      3600. * sqrt(minVariableDistance)
    );

    if (strlen(tempBuffer) + strlen(nearbyObjects) + 5 < MAX_NEARBY_OBJECTS_STRING) {
      if (strlen(nearbyObjects) > 0) {
        strcat(nearbyObjects, "; ");
      }

      strcat(nearbyObjects, tempBuffer);
    } else {
      printf(
        "ERROR: MAX_NEARBY_OBJECTS_STRING %d needs to be %zu\n",
        MAX_NEARBY_OBJECTS_STRING,
        strlen(tempBuffer) + strlen(nearbyObjects) + 5
      );
    }

    if (pGalaxyResult != NULL) {
      pGalaxyResult->variableFoundFlag = 1;
      pGalaxyResult->minVariableDistance = 3600.*sqrt(minVariableDistance);
      memcpy(&pGalaxyResult->galaxyRecVariable,pGalaxyRec,sizeof(GALAXYREC));
    }
  }

  if (pGalaxyRecGalaxy != NULL) {
    pGalaxyRec = pGalaxyRecGalaxy;
    sprintf(
      tempBuffer,
      "%s:%s @ %.1f",
      pGalaxyRec->catalogname,
      pGalaxyRec->galaxytype,
      3600. * sqrt(minGalaxyDistance)
    );

    if (strlen(tempBuffer) + strlen(nearbyObjects) + 5 < MAX_NEARBY_OBJECTS_STRING) {
      if (strlen(nearbyObjects) > 0) {
        strcat(nearbyObjects,"; ");
      }

      strcat(nearbyObjects,tempBuffer);
    } else {
      printf(
        "ERROR: MAX_NEARBY_OBJECTS_STRING %d needs to be %zu\n",
        MAX_NEARBY_OBJECTS_STRING,
        strlen(tempBuffer) + strlen(nearbyObjects) + 5
      );
    }

    if (pGalaxyResult != NULL) {
      pGalaxyResult->galaxyFoundFlag = 1;
      pGalaxyResult->minGalaxyDistance = 3600.*sqrt(minGalaxyDistance);
      memcpy(&pGalaxyResult->galaxyRecGalaxy,pGalaxyRec,sizeof(GALAXYREC));
    }
  }

  if (
    pNearestCatalogStar->nearestREFflagGoodMag >= 0 &&
    strlen(pNearestCatalogStar->nearestREFGoodMag) > 0 &&
    pNearestCatalogStar->nearestREFarcsecGoodMag < 9990.0
  ) {
    questionable[0] = '*';
    questionable[1] = 0;
    sprintf(
      tempBuffer,
      "%s (%.1f) @ %.1f%s",
      pNearestCatalogStar->nearestREFGoodMag,
      pNearestCatalogStar->nearestREFmagGoodMag,
      pNearestCatalogStar->nearestREFarcsecGoodMag,
      questionable
    );

    if (strlen(tempBuffer) + strlen(nearbyObjects) + 5 < MAX_NEARBY_OBJECTS_STRING) {
      if (strlen(nearbyObjects) > 0) {
        strcat(nearbyObjects,"; ");
      }

      strcat(nearbyObjects,tempBuffer);
      duplicateFlag = 1;
    } else {
      printf(
        "ERROR: MAX_NEARBY_OBJECTS_STRING %d needs to be %zu\n",
        MAX_NEARBY_OBJECTS_STRING,
        strlen(tempBuffer) + strlen(nearbyObjects) + 5
      );
    }
  }

  if (
    pNearestCatalogStar->nearestREFflag >= 0 &&
    strlen(pNearestCatalogStar->nearestREF) > 0 &&
    pNearestCatalogStar->nearestREFarcsec < 9990.0
  ) {
    if (duplicateFlag == 0 || strcmp(pNearestCatalogStar->nearestREF,pNearestCatalogStar->nearestREFGoodMag) != 0) {
      questionable[0] = 0;
      sprintf(
        tempBuffer,
        "%s (%.1f) @ %.1f%s",
        pNearestCatalogStar->nearestREF,
        pNearestCatalogStar->nearestREFmag,
        pNearestCatalogStar->nearestREFarcsec,
        questionable
      );

      if (strlen(tempBuffer) + strlen(nearbyObjects) + 5 < MAX_NEARBY_OBJECTS_STRING) {
        if (strlen(nearbyObjects) > 0) {
          strcat(nearbyObjects,"; ");
        }

        strcat(nearbyObjects,tempBuffer);
      } else {
        printf(
          "ERROR: MAX_NEARBY_OBJECTS_STRING %d needs to be %zu\n",
          MAX_NEARBY_OBJECTS_STRING,
          strlen(tempBuffer) + strlen(nearbyObjects) + 5
        );
      }
    }
  }
}


void
FreeGalaxyCommon(PGALAXYCOMMON pGalaxyCommon, int verbose)
{
  if (verbose > 0 && pGalaxyCommon->sdssCatalogFD >= 0) {
    printf(
      "FreeGalaxyCommon sdssCatalogAlloc %d sdssGscBinReadCount %d sdssGscCacheHitCount %d\n",
      pGalaxyCommon->sdssCatalogAlloc,
      pGalaxyCommon->sdssGscBinReadCount,
      pGalaxyCommon->sdssGscCacheHitCount
    );
  }

  if (pGalaxyCommon->sdssCatalogFD >= 0) {
    close(pGalaxyCommon->sdssCatalogFD);
    pGalaxyCommon->sdssCatalogFD = -1;
  }

  if (pGalaxyCommon->sdssIndexFD >= 0) {
    close(pGalaxyCommon->sdssIndexFD);
    pGalaxyCommon->sdssIndexFD = -1;
  }

  if (pGalaxyCommon->sdssBuffer != NULL) {
    free(pGalaxyCommon->sdssBuffer);
    pGalaxyCommon->sdssBuffer = NULL;
    pGalaxyCommon->sdss_bin_index = 1;
  }

  if (verbose > 0 && pGalaxyCommon->tycho2CatalogFD >= 0) {
    printf(
      "FreeGalaxyCommon tycho2CatalogAlloc %d tycho2GscBinReadCount %d tycho2GscCacheHitCount %d\n",
      pGalaxyCommon->tycho2CatalogAlloc,
      pGalaxyCommon->tycho2GscBinReadCount,
      pGalaxyCommon->tycho2GscCacheHitCount
    );
  }

  if (pGalaxyCommon->tycho2CatalogFD >= 0) {
    close(pGalaxyCommon->tycho2CatalogFD);
    pGalaxyCommon->tycho2CatalogFD = -1;
  }

  if (pGalaxyCommon->tycho2IndexFD >= 0) {
    close(pGalaxyCommon->tycho2IndexFD);
    pGalaxyCommon->tycho2IndexFD = -1;
  }

  if (pGalaxyCommon->tycho2Buffer != NULL) {
    free(pGalaxyCommon->tycho2Buffer);
    pGalaxyCommon->tycho2Buffer = NULL;
    pGalaxyCommon->tycho2_bin_index = 1;
  }

  if (pGalaxyCommon->galaxybuffer != NULL) {
    free(pGalaxyCommon->galaxybuffer);
  }

  if (pGalaxyCommon->inputBuffer != NULL) {
    free(pGalaxyCommon->inputBuffer);
  }

  if (pGalaxyCommon->nearbyBuffer != NULL) {
    free(pGalaxyCommon->nearbyBuffer);
  }

	if (pGalaxyCommon->nearbyDataHandle != NULL) {
		Close(pGalaxyCommon->nearbyDataHandle);
		pGalaxyCommon->nearbyDataHandle = NULL;
	}

	if (pGalaxyCommon->platesDataHandle != NULL) {
		Close(pGalaxyCommon->platesDataHandle);
		pGalaxyCommon->platesDataHandle = NULL;
	}

  if (pGalaxyCommon->limitingDataHandle != NULL) {
		Close(pGalaxyCommon->limitingDataHandle);
		pGalaxyCommon->limitingDataHandle = NULL;
	}

	if (pGalaxyCommon->nearbyIndexHandle != NULL) {
		Close(pGalaxyCommon->nearbyIndexHandle);
		pGalaxyCommon->nearbyIndexHandle = NULL;
	}

	if (pGalaxyCommon->platesIndexHandle != NULL) {
		Close(pGalaxyCommon->platesIndexHandle);
		pGalaxyCommon->platesIndexHandle = NULL;
	}

	if (pGalaxyCommon->limitingIndexHandle != NULL) {
		Close(pGalaxyCommon->limitingIndexHandle);
		pGalaxyCommon->limitingIndexHandle = NULL;
	}

  memset(pGalaxyCommon, 0, sizeof(GALAXYCOMMON));
}


int
PopulateGalaxyTree(
  PGALAXYCOMMON pGalaxyCommon,
  char *galaxy_name,
  int enableTransientSearch
) {
  int galaxyStatResult;
  struct stat galaxy_statbuf;
  File galaxy_handle = NULL;
  off_t galaxy_readItems = 0;
  off_t galaxy_index = 0;
  char *slashPtr;
  char sdss_name[MAX_BUFFER];
  char sdss_index_name[MAX_BUFFER];
  int sdssStatResult;
  struct stat sdss_statbuf;
  char tycho2_name[MAX_BUFFER];
  char tycho2_index_name[MAX_BUFFER];
  int tycho2StatResult;
  struct stat tycho2_statbuf;
  off_t total_count;
  char *dotPtr;
  PGALAXYREC pGalaxyRec;
  double x;
  double y;
  double z;
  double cosx;

  sdss_name[0] = 0;
  tycho2_name[0] = 0;

  if (pGalaxyCommon->ptree != NULL) {
    printf("ERROR: PopulateGalaxyTree ptree already populated\n");
    exit(-1);
  }

  if (pGalaxyCommon->initialized == 0) {
    InitGalaxyCommon(pGalaxyCommon);
  }

  pGalaxyCommon->enableTransientSearch = enableTransientSearch;
  pGalaxyCommon->ptree = kd_create(3);

  galaxyStatResult = stat(galaxy_name, &galaxy_statbuf);
  if (galaxyStatResult != 0) {
    printf("Could not find binary file %s\n",galaxy_name);
    return -1;
  }

  if (enableTransientSearch > 0) {
    if (strlen(galaxy_name) >= MAX_BUFFER + 1) {
      printf(
        "ERROR PopulateGalaxyTree galaxy_name size %zu exceeds %d\n",
        strlen(galaxy_name),
        MAX_BUFFER
      );
      exit(-1);
    }

    strcpy(sdss_name, galaxy_name);
    slashPtr = strrchr(sdss_name, '/');

    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = sdss_name;
    }

    *slashPtr = 0;
    strcat(sdss_name, "sdss.dat");

    if (strlen(sdss_name) >= MAX_BUFFER + 1) {
      printf(
        "ERROR PopulateGalaxyTree sdss_name size %zu exceeds %d\n",
        strlen(sdss_name),
        MAX_BUFFER
      );
      exit(-1);
    }

    sdssStatResult = stat(sdss_name, &sdss_statbuf);
    if (sdssStatResult != 0) {
      printf("Could not find binary file %s\n", sdss_name);
      return -1;
    }

    strcpy(tycho2_name, galaxy_name);
    slashPtr = strrchr(tycho2_name, '/');

    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = tycho2_name;
    }

    *slashPtr = 0;
    strcat(tycho2_name, "tycho2.dat");

    if (strlen(tycho2_name) >= MAX_BUFFER + 1) {
      printf(
        "ERROR PopulateGalaxyTree tycho2_name size %zu exceeds %d\n",
        strlen(tycho2_name),
        MAX_BUFFER
      );
      exit(-1);
    }

    tycho2StatResult = stat(tycho2_name, &tycho2_statbuf);
    if (tycho2StatResult != 0) {
      printf("Could not find binary file %s\n", tycho2_name);
      return -1;
    }
  } else {
    memset(&sdss_statbuf, 0, sizeof(sdss_statbuf));
    memset(&tycho2_statbuf, 0, sizeof(tycho2_statbuf));
  }

  pGalaxyCommon->galaxy_count = galaxy_statbuf.st_size / sizeof(GALAXYREC);
  pGalaxyCommon->sdss_count = sdss_statbuf.st_size / sizeof(GALAXYREC);
  pGalaxyCommon->tycho2_count = tycho2_statbuf.st_size / sizeof(GSCIMAGE);

  /* Here we just open the files for reading */
  if (enableTransientSearch > 0) {
    pGalaxyCommon->sdssCatalogFD = open(sdss_name, O_RDONLY, S_IRWXU|S_IRGRP);
    if (pGalaxyCommon->sdssCatalogFD < 0) {
      printf("ERROR Could not open catalog binary file %s\n", sdss_name);
      pGalaxyCommon->sdssCatalogFD = -1;
      exit(-1);
    }

    if (strlen(sdss_name) > MAX_BUFFER - 2) {
      printf("ERROR: PopulateGalaxyTree sdss_name %s exceeds MAX_BUFFER\n", sdss_name);
      exit(-1);
    }

    strcpy(sdss_index_name, sdss_name);
    dotPtr = strrchr(sdss_index_name, '.');

    if (dotPtr == NULL) {
      printf("ERROR: PopulateGalaxyTree sdss_name %s has no suffix\n", sdss_name);
      exit(-1);
    }

    *dotPtr = 0;
    strcat(sdss_index_name, ".idx");

    pGalaxyCommon->sdssIndexFD = open(sdss_index_name, O_RDONLY, S_IRWXU|S_IRGRP);
    if (pGalaxyCommon->sdssIndexFD < 0) {
      printf("ERROR Could not open catalog index file %s\n", sdss_index_name);
      pGalaxyCommon->sdssIndexFD = -1;
      exit(-1);
    }

    pGalaxyCommon->tycho2CatalogFD = open(tycho2_name, O_RDONLY, S_IRWXU|S_IRGRP);
    if (pGalaxyCommon->tycho2CatalogFD < 0) {
      printf("ERROR Could not open catalog binary file %s\n", tycho2_name);
      pGalaxyCommon->tycho2CatalogFD = -1;
      exit(-1);
    }

    if (strlen(tycho2_name) > MAX_BUFFER - 2) {
      printf("ERROR: PopulateGalaxyTree tycho2_name %s exceeds MAX_BUFFER\n", tycho2_name);
      exit(-1);
    }

    strcpy(tycho2_index_name, tycho2_name);
    dotPtr = strrchr(tycho2_index_name, '.');

    if (dotPtr == NULL) {
      printf("ERROR: PopulateGalaxyTree tycho2_name %s has no suffix\n", tycho2_name);
      exit(-1);
    }

    *dotPtr = 0;
    strcat(tycho2_index_name, ".idx");

    pGalaxyCommon->tycho2IndexFD = open(tycho2_index_name, O_RDONLY, S_IRWXU|S_IRGRP);
    if (pGalaxyCommon->tycho2IndexFD < 0) {
      printf("ERROR Could not open catalog index file %s\n", tycho2_index_name);
      pGalaxyCommon->tycho2IndexFD = -1;
      exit(-1);
    }
  }

  total_count = pGalaxyCommon->galaxy_count;

  pGalaxyCommon->galaxybuffer = (PGALAXYREC) calloc(galaxy_statbuf.st_size, 1);
  if (pGalaxyCommon->galaxybuffer == NULL) {
    printf("ERROR: failed to allocate galaxy buffer of size %ld\n", galaxy_statbuf.st_size);
    return -1;
  }

  galaxy_handle = Open(galaxy_name, "r");
  if (galaxy_handle == NULL) {
    printf("Could not open binary file %s\n", galaxy_name);
    return -1;
  }

  galaxy_readItems = Read(
    galaxy_handle,
    pGalaxyCommon->galaxybuffer,
    sizeof(GALAXYREC),
    pGalaxyCommon->galaxy_count
  );

  if (galaxy_readItems != pGalaxyCommon->galaxy_count) {
    printf(
      "ERROR: Read only %ld of %ld items from %s\n",
      galaxy_readItems,
      pGalaxyCommon->galaxy_count,
      galaxy_name
    );
    return -1;
  }

  for (galaxy_index = 0; galaxy_index < total_count; galaxy_index++) {
    pGalaxyRec = &pGalaxyCommon->galaxybuffer[galaxy_index];
    z = sin(DEGREES_TO_RAD * pGalaxyRec->dec);
    cosx = cos(DEGREES_TO_RAD * pGalaxyRec->dec);
    x = cosx * cos(DEGREES_TO_RAD * pGalaxyRec->ra);
    y = cosx * sin(DEGREES_TO_RAD * pGalaxyRec->ra);

    if (kd_insert3(pGalaxyCommon->ptree, x, y, z, pGalaxyRec) != 0) {
      printf("ERROR: fatal return from kd_insert3\n");
      exit(-1);
    }
  }

  return 0;
}


void
InitNearestCatalogStar(PNEARESTCATALOGSTAR pNearestCatalogStar)
{
  memset(pNearestCatalogStar, 0, sizeof(NEARESTCATALOGSTAR));
  pNearestCatalogStar->nearestREF[0] = 0;
  pNearestCatalogStar->nearestREFarcsec = 9999.0;
  pNearestCatalogStar->nearestREFflag = -1;
  pNearestCatalogStar->nearestREFGoodMag[0] = 0;
  pNearestCatalogStar->nearestREFarcsecGoodMag = 9999.0;
  pNearestCatalogStar->nearestREFflagGoodMag = -1;
}


int
OpenGalaxyFiles(PGALAXYCOMMON pGalaxyCommon, char *qualifier, char *prefix, char *photfilebase)
{
  char *catalogdir;
	int catalogNumber;
	char catalogname[MAX_BUFFER];
  char indexname[MAX_BUFFER];
  char *slashPtr;
  char *dotPtr;
	struct stat filestats;
  int result;

  if (pGalaxyCommon->limitingDataHandle != NULL || pGalaxyCommon->limitingIndexHandle != NULL) {
		printf("ERROR: OpenGalaxyFiles already invoked\n");
		exit(-1);
	}

  if (pGalaxyCommon->initialized == 0) {
    InitGalaxyCommon(pGalaxyCommon);
  }

  if (prefix == NULL) {
    prefix = "";
  }

  /* Load up the limiting database associated with the phot files */

  result = snprintf(catalogname, sizeof(catalogname), "%s/%slimiting.dat", photfilebase, prefix);
  if (result < 0) {
    printf("ERROR: buffer overrun for catalogname\n");
    exit(1);
  }

	if (*prefix) {
		printf("Opening limiting file %s\n", catalogname);
	}

	pGalaxyCommon->limitingDataHandle = Open(catalogname, "r");
	if (pGalaxyCommon->limitingDataHandle == NULL) {
		printf("ERROR: OpenGalaxyFiles could not open %s\n", catalogname);
		return -1;
	}

  result = snprintf(indexname, sizeof(indexname), "%s/%slimiting.idx", photfilebase, prefix);
  if (result < 0) {
    printf("ERROR: buffer overrun for indexname\n");
    exit(1);
  }

	pGalaxyCommon->limitingIndexHandle = Open(indexname, "r");
	if (pGalaxyCommon->limitingIndexHandle == NULL) {
		printf("ERROR: OpenGalaxyFiles could not open %s\n", indexname);
		return -1;
	}

	/* Figure out what binning we have here */
	result = stat(indexname, &filestats);

	if (result != 0) {
		printf("ERROR: failed to find size of file %s\n", indexname);
		exit(-1);
	}

	if (filestats.st_size == sizeof(STARINDEX) * gscBin01.total_gsc_bins) {
		pGalaxyCommon->pLimitingGscBin = &gscBin01;
		pGalaxyCommon->pNearbyGscBin = &gscBin01;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin02.total_gsc_bins) {
		pGalaxyCommon->pLimitingGscBin = &gscBin02;
		pGalaxyCommon->pNearbyGscBin = &gscBin02;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin04.total_gsc_bins) {
		pGalaxyCommon->pLimitingGscBin = &gscBin04;
		pGalaxyCommon->pNearbyGscBin = &gscBin04;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin08.total_gsc_bins) {
		pGalaxyCommon->pLimitingGscBin = &gscBin08;
		pGalaxyCommon->pNearbyGscBin = &gscBin08;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin16.total_gsc_bins) {
		pGalaxyCommon->pLimitingGscBin = &gscBin16;
		pGalaxyCommon->pNearbyGscBin = &gscBin16;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin32.total_gsc_bins) {
		pGalaxyCommon->pLimitingGscBin = &gscBin32;
		pGalaxyCommon->pNearbyGscBin = &gscBin32;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin64.total_gsc_bins) {
		pGalaxyCommon->pLimitingGscBin = &gscBin64;
		pGalaxyCommon->pNearbyGscBin = &gscBin64;
	} else {
		printf(
      "ERROR: No expected multiple of index size %zu in file size %ld\n",
      sizeof(STARINDEX),
      filestats.st_size
    );
		exit(-1);
	}

	/* Now look for the "nearby" files, using $DASCH_CATALOG */

	if (strlen(qualifier) > 0) {
		catalogNumber = GetCatalogNumber(qualifier);
		if (catalogNumber < 0) {
			printf("ERROR: illegal qualifier in OpenGalaxyFiles: %s\n", qualifier);
			exit(-1);
		}
	} else {
		catalogNumber = 0;
	}

	catalogdir = getenv("DASCH_CATALOG");
	if (catalogdir == NULL) {
		printf("ERROR: DASCH_CATALOG is not defined in OpenGalaxyFiles\n");
		return -1;
	}

	strcpy(catalogname, catalogdir);
	slashPtr = strrchr(catalogname, '/');

	if (slashPtr != NULL) {
		slashPtr++;
	} else {
		slashPtr = catalogname;
	}

	*slashPtr = 0;

	if (prefix != NULL) {
		strcat(catalogname, prefix);
	}

	strcat(catalogname, "nearbylimiting.dat");
	strcpy(indexname, catalogname);
	dotPtr = strrchr(indexname, '.');

	if (dotPtr != NULL) {
		*dotPtr = 0;
	}

	strcat(indexname, ".idx");

	/* For backward compatability, we must assume that the plates files also has nearby galaxy data */

	pGalaxyCommon->nearbyDataHandle = Open(catalogname, "r");
	if (pGalaxyCommon->nearbyDataHandle == NULL) {
		return 0;
	}

	pGalaxyCommon->nearbyIndexHandle = Open(indexname, "r");
	if (pGalaxyCommon->nearbyIndexHandle == NULL) {
		Close(pGalaxyCommon->nearbyDataHandle);
		pGalaxyCommon->nearbyDataHandle = NULL;
		return 0;
	}

	result = stat(indexname, &filestats);
	if (result != 0) {
		printf("ERROR: failed to find size of file %s\n", indexname);
		exit(-1);
	}

	if (filestats.st_size == sizeof(STARINDEX) * gscBin01.total_gsc_bins) {
		pGalaxyCommon->pNearbyGscBin = &gscBin01;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin02.total_gsc_bins) {
		pGalaxyCommon->pNearbyGscBin = &gscBin02;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin04.total_gsc_bins) {
		pGalaxyCommon->pNearbyGscBin = &gscBin04;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin08.total_gsc_bins) {
		pGalaxyCommon->pNearbyGscBin = &gscBin08;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin16.total_gsc_bins) {
		pGalaxyCommon->pNearbyGscBin = &gscBin16;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin32.total_gsc_bins) {
		pGalaxyCommon->pNearbyGscBin = &gscBin32;
	} else if (filestats.st_size == sizeof(STARINDEX) * gscBin64.total_gsc_bins) {
		pGalaxyCommon->pNearbyGscBin = &gscBin64;
	} else {
		printf(
      "ERROR: No expected multiple of index size %zu in file size %ld\n",
      sizeof(STARINDEX),
      filestats.st_size
    );
		exit(-1);
	}

	pGalaxyCommon->nearbyEnabled = 1;
 	if (prefix != NULL) {
		printf("Opening nearby file %s\n", catalogname);
	}

	return 0;
}


void
DumpLimitingTable(PGALAXYCOMMON pGalaxyCommon)
{
	int plateIndex;
	PPLATELIMITINGREC pPlateLimitingRec;

	if (pGalaxyCommon->plateCount <= 0) {
		printf("No limiting magnitudes are available\n");
		return;
	}

	for (plateIndex = 0; plateIndex < pGalaxyCommon->plateCount; plateIndex++) {
		 pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
		 printf(
      "Entry %d Plate %5s%05d_%02d s%d limiting_mag_local %.2f\n",
      plateIndex,
      GetSeriesString(pPlateLimitingRec->seriesId, 0),
      pPlateLimitingRec->plateNumber,
      pPlateLimitingRec->mosaicNumber,
      pPlateLimitingRec->solutionNumber,
      pPlateLimitingRec->limiting_mag_local
    );
	}
}


/* Returns 1 if a Tycho2 star is within tycho2srcsec and brighter than tycho2mag */
int
FindNearestTycho2Star(
  PGALAXYCOMMON pGalaxyCommon,
  double ra,
  double dec,
  double tycho2mag,
  double tycho2arcsec
) {
  int cur_bin_index;
  int decBin;
  int raBin;
  int readItems;
  int readBytes;
  PGSCIMAGE pGscImage;
  double curDistance;
  int readIndex;
  double factor;

  if (pGalaxyCommon->tycho2CatalogFD >= 0 && pGalaxyCommon->tycho2IndexFD >= 0) {
    cur_bin_index = GetGSCBin(&gscBin08, ra, dec, &decBin, &raBin, "FindNearestGalaxy");

    if (cur_bin_index != pGalaxyCommon->tycho2_bin_index) {
      pGalaxyCommon->tycho2_bin_index = cur_bin_index;
      pGalaxyCommon->tycho2GscBinReadCount++;

      lseek(pGalaxyCommon->tycho2IndexFD, pGalaxyCommon->tycho2_bin_index * sizeof(STARINDEX), SEEK_SET);
      readBytes = read(pGalaxyCommon->tycho2IndexFD, &pGalaxyCommon->curtycho2StarIndex, sizeof(STARINDEX));
      readItems = readBytes / sizeof(STARINDEX);
      if (readItems != 1) {
        printf(
          "ERROR FindNearestGalaxy reading star index file errno: %d %s gsc_bin_index %d\n",
          errno,
          strerror(errno),
          pGalaxyCommon->tycho2_bin_index
        );
        exit(-1);
      }

      if (pGalaxyCommon->curtycho2StarIndex.numStars > 0) {
        if (pGalaxyCommon->curtycho2StarIndex.binNumber != pGalaxyCommon->tycho2_bin_index) {
          printf(
            "ERROR: FindNearestGalaxy gsc_bin_index %d %d do not agree with catalog index\n",
            pGalaxyCommon->tycho2_bin_index,
            pGalaxyCommon->curtycho2StarIndex.binNumber
          );
          exit(-1);
        }

        if (pGalaxyCommon->curtycho2StarIndex.numStars > pGalaxyCommon->tycho2CatalogAlloc) {
          if (pGalaxyCommon->tycho2Buffer != NULL) {
            free(pGalaxyCommon->tycho2Buffer);
          }

          pGalaxyCommon->tycho2CatalogAlloc = pGalaxyCommon->curtycho2StarIndex.numStars;
          pGalaxyCommon->tycho2Buffer = (PGSCIMAGE) calloc(pGalaxyCommon->tycho2CatalogAlloc, sizeof(GSCIMAGE));
          if (pGalaxyCommon->tycho2Buffer == NULL) {
            printf(
              "ERROR: FindNearestGalaxy failed to allocate tycho2Buffer of size %d\n",
              pGalaxyCommon->tycho2CatalogAlloc
            );
            exit(-1);
          }
        }

        lseek(pGalaxyCommon->tycho2CatalogFD, pGalaxyCommon->curtycho2StarIndex.offset, SEEK_SET);
        readBytes = read(
          pGalaxyCommon->tycho2CatalogFD,
          pGalaxyCommon->tycho2Buffer,
          sizeof(GSCIMAGE) * pGalaxyCommon->curtycho2StarIndex.numStars
        );
        readItems = readBytes / sizeof(GSCIMAGE);

        if (readItems != pGalaxyCommon->curtycho2StarIndex.numStars) {
          fprintf(
            stderr,
            "ERROR reading gsc catalog file at index %d\n",
            pGalaxyCommon->tycho2_bin_index
          );
          exit(-1);
        }
      }
    } else {
      pGalaxyCommon->tycho2GscCacheHitCount++;
    }

    if (pGalaxyCommon->curtycho2StarIndex.numStars > 0) {
      factor = cos(DEGREES_TO_RAD * dec);

      for (readIndex = 0; readIndex < pGalaxyCommon->curtycho2StarIndex.numStars; readIndex++) {
        pGscImage = &pGalaxyCommon->tycho2Buffer[readIndex];
        if (pGscImage->Stdmag > tycho2mag || pGscImage->Stdmag <= 0) {
          continue;
        }

        curDistance = 3600. * sqrt(sqr(pGscImage->dec - dec) + sqr(factor * (pGscImage->ra - ra)));
        if (curDistance < tycho2arcsec) {
          return 1;
        }
      }
    }
  }

  return 0;
}


/* Sort based on numStars */
int
StarIndexCompare(const void *first, const void *second)
{
  int numberFirst = ((PSTARINDEX) first)->numStars;
  int numberSecond = ((PSTARINDEX) second)->numStars;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  } else {
    return 0;
  }
}


void
SearchLimiting(
  PGALAXYCOMMON pGalaxyCommon,
  int printLimit,
  char* coveragefile
) {
  int cur_read_bin;
  PSTARINDEX resultTable;
  int result_index;
  int result_count = 0;
  PSTARINDEX readBuffer;
  int read_index;
  int readItems;
  PSTARINDEX pResult;
  PSTARINDEX pSmallestResult = NULL;
  PSTARINDEX pLimiting;
  double ra;
  double dec;
  FILE *coveragehandle = NULL;
  int coverageFlag = 0;
  int binIndex;
  PBININDEX pBinIndex;
  int decBinIndex;
  int raBinIndex;
  PGSCBIN pGscBin = &gscBin01;
  double decval;
  double raval;

  if (coveragefile != NULL) {
    if (coveragefile[0] != 0) {
      coveragehandle = fopen(coveragefile, "wt");

      if (coveragehandle == NULL) {
        printf("ERROR: Failed to open the coverage file %s\n", coveragefile);
        exit(-1);
      } else {
        printf("coverage file %s\n", coveragefile);
        coverageFlag = 1;
      }
    }
  }

  resultTable = (PSTARINDEX) calloc(printLimit, sizeof(STARINDEX));
  if (resultTable == NULL) {
    printf("ERROR: failed to allocate resultTable in SearchLimiting\n");
    exit(-1);
  }

  if (coverageFlag == 0) {
    readBuffer = (PSTARINDEX) calloc(GSC_BIN_INCREMENT, sizeof(STARINDEX));
    if (readBuffer == NULL) {
      printf("ERROR: failed to allocate readBuffer in SearchLimiting\n");
      exit(-1);
    }

    for (cur_read_bin = 0; cur_read_bin < pGalaxyCommon->pLimitingGscBin->total_gsc_bins; cur_read_bin += GSC_BIN_INCREMENT) {
      Seek((File) pGalaxyCommon->limitingIndexHandle, cur_read_bin * sizeof(STARINDEX), SEEK_SET);
      readItems = Read(
        (File) pGalaxyCommon->limitingIndexHandle,
        readBuffer,
        sizeof(STARINDEX),
        GSC_BIN_INCREMENT
      );
      if (readItems <= 0) {
        printf("Error: SearchLimiting can not read bin %d\n", cur_read_bin);
        exit(-1);
      }

      for (read_index = 0; read_index < readItems; read_index++) {
        pLimiting = &readBuffer[read_index];

        if (result_count < printLimit) {
          pResult = &resultTable[result_count];
          memcpy(pResult, pLimiting, sizeof(STARINDEX));
          result_count++;

          if (pSmallestResult == NULL) {
            pSmallestResult = pResult;
          } else if (pSmallestResult->numStars > pResult->numStars) {
            pSmallestResult = pResult;
          }
        } else if (pLimiting->numStars > pSmallestResult->numStars) {
          memcpy(pSmallestResult, pLimiting, sizeof(STARINDEX));

          for (result_index = 0; result_index < printLimit; result_index++) {
            pResult = &resultTable[result_index];

            if (pResult->numStars < pSmallestResult->numStars) {
              pSmallestResult = pResult;
            }
          }
        }
      }
    }

    qsort((void *)resultTable, printLimit, sizeof(STARINDEX), StarIndexCompare);

    for (result_index = 0; result_index < printLimit; result_index++) {
      pResult = &resultTable[result_index];
      GetBinCenter(pGalaxyCommon->pLimitingGscBin, pResult->binNumber, &ra, &dec, "SearchLimiting");
      printf("binNumber %8d numStars %8d ra,dec %9.4f %8.4f\n", pResult->binNumber, pResult->numStars, ra, dec);
    }
  } else {
    for (binIndex = 0; binIndex < pGscBin->total_gsc_bins; binIndex++) {
      pBinIndex = GetSubbins(pGscBin, binIndex, &raBinIndex, &decBinIndex, "galaxyutils");
      raval = (360.0 * raBinIndex) / pBinIndex->numBins + 180.0 / pBinIndex->numBins;
      decval = decBinIndex * pGscBin->bin_size - 90.0 + pGscBin->bin_size / 2.0;

      if (raval < 0 || raval > 360.0) {
        printf("raval error\n");
      }

      if (decval < -90.0 || decval > 90.0) {
        printf("decval error\n");
      }

      LoadGalaxyTable(pGalaxyCommon, raval, decval);
      fprintf(
        coveragehandle,
        "%f %f %d\n",
        360.0 - raval,
        decval,
        pGalaxyCommon->plateCount
      );
    }
  }

  if (coveragehandle != NULL) {
    fclose(coveragehandle);
  }
}

/* Dec 28, 2012 Edward J. Los - Initial Version
 * Jan  4, 2013 Edward J. Los - kdtree search support
 * Feb 11, 2013 Edward J. Los - Add nearestREFmag support.
 * Feb 23, 2013 Edward J. Los - Add PNEARESTCATALOGSTAR to support printing of additional nearby objects.
 * May 26, 2015 Edward J. Los - Add pGalaxyResult to support the SDSS catalog
 * Jun 29, 2015 Edward J. Los - Use mmap() and munmap() to map the catalogs and avoid exhausting system resources
 * Jul  6, 2015 Edward J. Los - Convert sdss catalog access to an indexed file
 *                              Correct minGalaxyDistance to get the closest galaxy
 * Aug  3, 2015 Edward J. Los - Read in the Tycho2 catalog to check transient candidates against bright stars
 *                              Remove USE_MMAP_FUNCTION because it didn't work.  Enable USE_INDEXED_SDSS permanently
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflict
 * Feb 24, 2017 Edward J. Los - Add limitingDataHandle and limitingIndexHandle to avoid conflict
 *                              Add sanity checks to prevent double initilization of galaxy common
 * Apr  4, 2018 Edward J. Los - Add SearchLimiting
 */
