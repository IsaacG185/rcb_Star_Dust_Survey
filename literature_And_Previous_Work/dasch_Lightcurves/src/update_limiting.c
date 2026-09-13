// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* This program reads the photometry database to find the limiting magnitude for
 * each plate and each square degree of the sky. The results are limiting
 * magnitude maps of the sky as a function of limiting magnitude and a series of
 * binary files which can be read to find the number of plates scanned for a
 * given square degree of sky and what the limiting magnitude of those plates
 * happens to be.
 *
 * Note that this tool is constructed to do partial updates -- it only processes
 * plates where `photplates.versionId@KEY@ >= photglobal.minVersionId@KEY@`.
 * AFAICT, nothing in the pipeline updates the global setting, though! Also,
 * this of course means that plates must be logged in photplates in order to be
 * included in processing.
 *
 * Inputs:
 *
 * - $DASCH_PHOT_MAGNITUDES@KEY@/magNNN/magNNN/magNNN/magNNNNNNNNN.db
 * - $(dirname $DASCH_CATALOG)/galaxy.dat (if processing in "nearby" mode)
 *
 * Modified in-place:
 *
 * - Main output specified with `-o` argument, traditionally
 *   `{refcat}limiting.dat`
 * - Index file: `-o` argument with `.dat` replaced with  `.idx`
 *
 * Outputs:
 *
 * - Logfile associated with `-l` argument
 *
 * Database updates:
 *
 * - Update `photometry.photglobal` timestamp
 *
 * Environment variables:
 *
 * - DASCH_CATALOG
 * - DASCH_PHOT_ROOT
 * - Database vars:
 *   - DASCH_MYSQLHOST
 *   - DASCH_PASSWORD
 *   - DASCH_USERNAME
 */

#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>

/* Starlink */
#include "table.h"

/* MySQL / MariaDB */
#include "mysql.h"
#include "errmsg.h"

/* libwcs */
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

/* Astrometry.net */
#include "kdtree.h"

/* DASCH pipeline local */
#include "scandb.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include "photometryutils.h"
#include "magdeputil.h"

#define MAXIMUM_BUFFER_MEMORY 2000000000 /* Need twice this to make things work ! */
#define MAX_BUFFER 256
#define GSCBINREC_ALLOC_COUNT 500
#define PLATE_TABLE_INCREMENT 4000
#define AC_SCALE_FACTOR 0.98

extern GSCBIN gscBin01;
extern GSCBIN gscBin02;
extern GSCBIN gscBin04;
extern GSCBIN gscBin08;
extern GSCBIN gscBin16;
extern GSCBIN gscBin32;
extern GSCBIN gscBin64;

static PGSCBIN pGscBin = &gscBin02;

extern int binarray[X_DMAGBINS_NORMAL*Y_DMAGBINS_NORMAL];
extern char *releaseFieldText[RELEASE_FIELD_MAX+1];


typedef struct _gsclimitingbinrec {
	double limiting_mag_local;
	int gsc_bin_index;
	int local_bin_index;
} GSCLIMITINGBINREC, *PGSCLIMITINGBINREC;

typedef struct _plateEntry {
	double geoJulianDate;
	int errorFlag;
	int versionId;
	int seriesId;
	int plateNumber;
	int mosaicNumber;
	int rotation;
	int solutionNumber;
	int leftMargin;
	int rightMargin;
	int topMargin;
	int bottomMargin;
	int maxBin;
	int minBin;
	int nx; /* X-axis local bins */
	int ny; /* Y-axis local bins */
	int totalGoodBins;
	int totalPartialBins;
	int gscLimitingBinCount;
	int gscLimitingBinAlloc;
	PGSCLIMITINGBINREC gscLimitingBinBuffer;
	MOSAIC mosaic;
	char series[MAX_SERIES_STRING];
	char Plate[MAX_PLATE_NAME];
} PLATEENTRY, *PPLATEENTRY;

typedef struct _limitingmag {
	double limiting_mag_local;
	int gsc_bin_index;
	int local_bin_index;
	int valid;
} LIMITINGMAG, *PLIMITINGMAG;

/* Can be as low as 32 bytes, but padded to equal PPLATELIIMITINGREC */
typedef struct _file_limiting_record {
	int galaxyflag; /* Padding */
	int galaxyversion; /* Padding*/
	double unused2; /* Padding */
	double geoJulianDate;
	double limiting_mag_local;
	int plateNumber;
	int versionId;
	unsigned char seriesId;
	unsigned char mosaicNumber;
	unsigned char solutionNumber;
	unsigned char unused;
	int gsc_bin_index;
} FILELIMITINGREC, *PFILELIMITINGREC;

typedef struct _updateCommon {
	long long currentMemory;
	long long maxMemory;
	int verbose;
	int catalogNumber;
	char catalogString[MAX_BUFFER];
	MYSQL *pConnection;
	MYSQL *pPhotConnection;
	time_t curTime;
	time_t startTime;
	time_t plateStartTime;
	int numPlates;

	char *photfilebase;
	FILE *logHandle;
	File outputHandle;
	File indexHandle;
	File temporaryHandle;

	int magnitudeFileFlag;
	int plateTableAlloc;
	int nearbyFlag;
	int platesFlag;
	PPLATEENTRY pPlateTable;
	size_t fileWrittenCount;
	size_t fileHandledCount;

	int *file_written_table; /* Number of records written to the file for this bin */
	int *file_handled_table; /* Number of records re-written for this bin */
	PFILELIMITINGREC pFileBuffer1;
	PFILELIMITINGREC pFileBuffer2;
	size_t file_alloc_count;
	off_t platerec_bytes; /* Size written to the limiting magnitude tables */
	int maxGalaxyBinCount;
	int aveGalaxyBinCount;
	int countRelease[RELEASE_FIELD_MAX];
	int maxPlateBinCount;
	int avePlateBinCount;
	int zeroBinFlushCount;
	int staleBinFlushCount;
	int staleBinMaximum;
	int lastBinFlushed;
	double searchRadius; /* Search radius in degrees */
	char *checkedTable;        /* 1 = queued; 2 = checked */
	char *releaseMappingTable; /* Maps bin numbers to release field numbers */
	int *binWorkQueue;       /* Queue of bins to be checked */
	int workQueueSize;       /* Current size of the bin_work_queue */
	PHOTGLOBAL photGlobal;
} UPDATECOMMON, *PUPDATECOMMON;


static void
BuildPlatesTable(PUPDATECOMMON pUpdateCommon)
{
	int nvals;
	MYSQL_RES *res_ptr;
	MYSQL_ROW sqlrow;
	int res;
	char queryString[MAX_QUERY_STRING];
	int result_size;
	int plateIndex = 0;
	PPLATEENTRY pPlate;

	sprintf(
		queryString,
		"SELECT seriesId,plateNumber,versionId%s,leftMargin,rightMargin,topMargin,bottomMargin FROM photplates WHERE versionId%s IS NOT NULL and versionId%s >= %d;",
		pUpdateCommon->catalogString,
		pUpdateCommon->catalogString,
		pUpdateCommon->catalogString,
		pUpdateCommon->photGlobal.minVersionId[pUpdateCommon->catalogNumber]
	);
	res = ExecuteQuery(pUpdateCommon->pPhotConnection, queryString);

	if (!res) {
		res_ptr = mysql_store_result(pUpdateCommon->pPhotConnection);

		if (res_ptr) {
			result_size = (int) mysql_num_rows(res_ptr);
			pUpdateCommon->plateTableAlloc = result_size + PLATE_TABLE_INCREMENT;
			pUpdateCommon->pPlateTable = (PPLATEENTRY) calloc(pUpdateCommon->plateTableAlloc,sizeof(PLATEENTRY));
			pUpdateCommon->currentMemory += pUpdateCommon->plateTableAlloc*sizeof(PLATEENTRY);

			if (pUpdateCommon->currentMemory > pUpdateCommon->maxMemory) {
				pUpdateCommon->maxMemory = pUpdateCommon->currentMemory;
			}

			while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
				pPlate = &pUpdateCommon->pPlateTable[plateIndex];
				memset(pPlate, 0, sizeof(PLATEENTRY));

				if (sqlrow[0]) {
					nvals = sscanf(sqlrow[0],"%d",&pPlate->seriesId);
					if (nvals != 1) {
						continue;
					}
				} else {
					continue;
				}

				if (sqlrow[1]) {
					nvals = sscanf(sqlrow[1],"%d",&pPlate->plateNumber);
					if (nvals != 1) {
						continue;
					}
				} else {
					continue;
				}

				if (sqlrow[2]) {
					nvals = sscanf(sqlrow[2],"%d",&pPlate->versionId);
					if (nvals != 1) {
						continue;
					}
				} else {
					continue;
				}

				if (sqlrow[3]) {
					nvals = sscanf(sqlrow[3],"%d",&pPlate->leftMargin);
					if (nvals != 1) {
						continue;
					}
				} else {
					continue;
				}

				if (sqlrow[4]) {
					nvals = sscanf(sqlrow[4],"%d",&pPlate->rightMargin);
					if (nvals != 1) {
						continue;
					}
				} else {
					continue;
				}

				if (sqlrow[5]) {
					nvals = sscanf(sqlrow[5],"%d",&pPlate->topMargin);
					if (nvals != 1) {
						continue;
					}
				} else {
					continue;
				}

				if (sqlrow[6]) {
					nvals = sscanf(sqlrow[6],"%d",&pPlate->bottomMargin);
					if (nvals != 1) {
						continue;
					}
				} else {
					continue;
				}

				strcpy(pPlate->series,GetSeriesString(pPlate->seriesId,0));
				plateIndex++;
				if (plateIndex > result_size) {
					printf("ERROR: plateIndex %d greater than result_size %d at line %d\n", plateIndex, result_size, __LINE__);
				}
			}

			mysql_free_result(res_ptr);
		} else {
			printf("ERROR: BuildPlatesTable failed\n");
		}
	}

	printf(
		"GetPlatesTable read %d plates with versionId%s >= %d\n",
		plateIndex,
		pUpdateCommon->catalogString,
		pUpdateCommon->photGlobal.minVersionId[pUpdateCommon->catalogNumber]
	);
	pUpdateCommon->numPlates = plateIndex;
}


static void
StoreLimitingMagnitude(
	PUPDATECOMMON pUpdateCommon,
	PPLATEENTRY pPlate,
	int gsc_bin_index,
	int gsc_bin_count,
	double limiting_mag_local_average,
	int local_bin_index,
	int local_bin_index2
) {
	PGSCLIMITINGBINREC pGscLimitingBin;
	int index;
	double limiting_mag_local;

	if (gsc_bin_index < 0 || gsc_bin_count <= 0) {
		return;
	}

	limiting_mag_local = limiting_mag_local_average / (1.0 * gsc_bin_count);

	if (pPlate->gscLimitingBinCount > 0 && pPlate->gscLimitingBinCount >= pPlate->gscLimitingBinAlloc) {
		printf("ERROR: line %d gscLimitingBinCount %d >= gscLimitingBinAlloc %d\n", __LINE__, pPlate->gscLimitingBinCount, pPlate->gscLimitingBinAlloc);
		exit(-1);
	}

	for (index = 0; index < pPlate->gscLimitingBinCount; index++) {
		pGscLimitingBin = &pPlate->gscLimitingBinBuffer[index];

		if (pGscLimitingBin->gsc_bin_index == gsc_bin_index) {
			pGscLimitingBin->limiting_mag_local = limiting_mag_local;
			break;
		}
	}
}


static void
SaveLimitingMagnitudes(
	PUPDATECOMMON pUpdateCommon,
	PPLATEENTRY pPlate,
	PLIMITINGMAG limiting_mag_table,
	int *missingBinCount,
	int *errorPrintFlag
) {
	PGSCLIMITINGBINREC pGscLimitingBin;
	PLIMITINGMAG pLimitingMag;
	FILELIMITINGREC fileRec;
	PFILELIMITINGREC pFileRec = &fileRec;
	off_t writeItems = 0;
	int index;

	if (pPlate->gscLimitingBinCount > 0 && pPlate->gscLimitingBinCount >= pPlate->gscLimitingBinAlloc) {
		printf("ERROR: line %d gscLimitingBinCount %d >= gscLimitingBinAlloc %d\n", __LINE__, pPlate->gscLimitingBinCount,pPlate->gscLimitingBinAlloc);
		exit(-1);
	}

	for (index = 0; index < pPlate->gscLimitingBinCount; index++) {
		pGscLimitingBin = &pPlate->gscLimitingBinBuffer[index];
		if (pGscLimitingBin->limiting_mag_local == 0) {
			printf("ERROR: limiting_mag_local (3) is %f\n", pGscLimitingBin->limiting_mag_local);
			exit(-1);
		}

		if (pGscLimitingBin->limiting_mag_local < 0.0) {
			/* We missed this bin in the averaging process */
			if (pGscLimitingBin->local_bin_index >= TOTAL_DMAGBINS_NORMAL) {
				printf("ERROR: line %d local_bin_index %d exceeds %d\n", __LINE__, pGscLimitingBin->local_bin_index, TOTAL_DMAGBINS_NORMAL);
				exit(-1);
			}

			pLimitingMag = &limiting_mag_table[pGscLimitingBin->local_bin_index];
			if (pLimitingMag->valid != 0 && pLimitingMag->limiting_mag_local < 90.0) {
				if (pLimitingMag->local_bin_index != pGscLimitingBin->local_bin_index) {
					printf("ERROR:local bin indexing problem in line SaveLimitingMagnitudes\n");
					exit(-1);
				}

				(*missingBinCount)++;
				pGscLimitingBin->limiting_mag_local = pLimitingMag->limiting_mag_local;

				if (pLimitingMag->limiting_mag_local == 0) {
					printf("ERROR: limiting_mag_local (2) is %f\n", pLimitingMag->limiting_mag_local);
					exit(-1);
				}
			} else {
				continue;
			}
		}

		if (pGscLimitingBin->gsc_bin_index <= pUpdateCommon->lastBinFlushed) {
			pUpdateCommon->staleBinFlushCount++;

			if (pGscLimitingBin->gsc_bin_index == 0) {
				pUpdateCommon->zeroBinFlushCount++;
			} else {
				if (pUpdateCommon->lastBinFlushed-pGscLimitingBin->gsc_bin_index > pUpdateCommon->staleBinMaximum) {
					pUpdateCommon->staleBinMaximum = (pUpdateCommon->lastBinFlushed-pGscLimitingBin->gsc_bin_index);
				}
			}

			return;
		}

		memset(pFileRec, 0, sizeof(FILELIMITINGREC));
		pFileRec->limiting_mag_local = pGscLimitingBin->limiting_mag_local;

		if (pFileRec->limiting_mag_local <= 0) {
			if (*errorPrintFlag < 100) {
				printf("ERROR: limiting_mag_local (4) is %f for gsc_bin_index %d and local_bin_index %d seriesId %d plateNumber %d\n",
					pGscLimitingBin->limiting_mag_local,
					pGscLimitingBin->gsc_bin_index,
					pGscLimitingBin->local_bin_index,
					pPlate->seriesId,
					pPlate->plateNumber);
				(*errorPrintFlag)++;
			}
			pFileRec->limiting_mag_local = 99.0;
		}

		pFileRec->galaxyflag = PLATE_FLAG;  /* Padding -- not really needed here */
		pFileRec->galaxyversion = PLATE_VERSION; /* Padding -- not really needed here */
		pFileRec->gsc_bin_index = pGscLimitingBin->gsc_bin_index;
		pFileRec->geoJulianDate = pPlate->geoJulianDate;
		pFileRec->seriesId = pPlate->seriesId;
		pFileRec->plateNumber = pPlate->plateNumber;
		pFileRec->mosaicNumber = pPlate->mosaicNumber;
		pFileRec->solutionNumber = pPlate->solutionNumber;
		pFileRec->versionId = pPlate->versionId;
		pFileRec->unused = 0;
		pFileRec->unused2 = 0;

		pUpdateCommon->file_written_table[pGscLimitingBin->gsc_bin_index]++;
		writeItems = Write(pUpdateCommon->temporaryHandle, pFileRec, sizeof(FILELIMITINGREC), 1);
		if (writeItems != 1) {
			printf("ERROR writing %ld of %d items to the limiting magnitude table\n", writeItems, 1);
			exit(-1);
		}

		pUpdateCommon->fileWrittenCount++;
	}
}


static int
PlaterecCompare(const void *first, const void *second)
{
	PPLATELIMITINGREC pPlateRecFirst = (PPLATELIMITINGREC) first;
	PPLATELIMITINGREC pPlateRecSecond = (PPLATELIMITINGREC) second;

	if (pPlateRecFirst->seriesId > pPlateRecSecond->seriesId) {
		return 1;
	} else if (pPlateRecFirst->seriesId < pPlateRecSecond->seriesId) {
		return -1;
	}

	if (pPlateRecFirst->plateNumber > pPlateRecSecond->plateNumber) {
		return 1;
    } else if (pPlateRecFirst->plateNumber < pPlateRecSecond->plateNumber) {
		return -1;
    }

	if (pPlateRecFirst->mosaicNumber > pPlateRecSecond->mosaicNumber) {
		return 1;
	} else if (pPlateRecFirst->mosaicNumber < pPlateRecSecond->mosaicNumber) {
		return -1;
	}

	if (pPlateRecFirst->solutionNumber > pPlateRecSecond->solutionNumber) {
		return 1;
	} else if (pPlateRecFirst->solutionNumber < pPlateRecSecond->solutionNumber) {
		return -1;
	}

	return 0;
}


static int
BinrecCompare(const void *first, const void *second)
{
	PGSCLIMITINGBINREC pBinRecFirst = (PGSCLIMITINGBINREC) first;
	PGSCLIMITINGBINREC pBinRecSecond = (PGSCLIMITINGBINREC) second;

	if (pBinRecFirst->local_bin_index >= TOTAL_DMAGBINS_NORMAL) {
		printf("ERROR: line %d local_bin_index %d exceeds %d\n", __LINE__, pBinRecFirst->local_bin_index, TOTAL_DMAGBINS_NORMAL);
		exit(-1);
	}
	if (pBinRecSecond->local_bin_index >= TOTAL_DMAGBINS_NORMAL) {
		printf("ERROR: line %d local_bin_index %d exceeds %d\n", __LINE__, pBinRecSecond->local_bin_index, TOTAL_DMAGBINS_NORMAL);
		exit(-1);
	}

	if (pBinRecFirst->gsc_bin_index > pBinRecSecond->gsc_bin_index) {
		return 1;
	} else if (pBinRecFirst->gsc_bin_index < pBinRecSecond->gsc_bin_index) {
		return -1;
	}

	return 0;
}


static int
GscBinCompare(const void *first, const void *second)
{
	double gsc_bin_indexFirst = ((PFILELIMITINGREC) first)->gsc_bin_index;
	double gsc_bin_indexSecond = ((PFILELIMITINGREC) second)->gsc_bin_index;

	if (gsc_bin_indexFirst > gsc_bin_indexSecond) {
		return 1;
	} else if (gsc_bin_indexFirst < gsc_bin_indexSecond) {
		return -1;
	}

	return 0;
}


static void
FlushBinCache(
	PUPDATECOMMON pUpdateCommon,
	PGALAXYCOMMON pGalaxyCommon,
	int minGscBin,
	int maxGscBin,
	int entryCount
) {
	int gsc_bin_index;
	STARINDEX starIndex;
	PSTARINDEX pStarIndex = &starIndex;
	int numObjects;
	off_t writeItems = 0;
	double ra;
	double dec;
	double x;
	double y;
	double z;
	double cosx;
	double pt[3];
	double pos[3];
	struct kdres *presults;
	PGALAXYREC pGalaxyRec;
	int totalWrite = 0;
	int resultSize;
	int curMinBin;
	int readItems;
	int totalReadItems;
	int readIndex;
	int curHandledCount;
	int curHandledIndex;
	int totalWriteItems = 0;
	PFILELIMITINGREC pFileRec1;
	PFILELIMITINGREC pFileRec2;
	PPLATELIMITINGREC pPlateRec1;
	PPLATELIMITINGREC pPlateRec2;

	time(&pUpdateCommon->curTime);
	pUpdateCommon->curTime -= pUpdateCommon->startTime;

	printf(
		"FlushBinCache for bins %8d to %8d, entries %8d at %ld seconds\n",
		minGscBin,
		maxGscBin,
		entryCount,
		pUpdateCommon->curTime
	);

	/* Write out the bin cache */
	curMinBin = minGscBin;
	if (curMinBin != pUpdateCommon->lastBinFlushed + 1) {
		printf("ERROR: bin flush sequence error\n");
	}

	/* Rewind our temporary file */
	Seek(pUpdateCommon->temporaryHandle, 0, SEEK_SET);

	/* Now read in the requested records */
	curHandledCount = 0;
	totalReadItems = 0;

	while (totalReadItems < pUpdateCommon->fileWrittenCount) {
		readItems = Read(
			pUpdateCommon->temporaryHandle,
			pUpdateCommon->pFileBuffer1,
			sizeof(FILELIMITINGREC),
			pUpdateCommon->file_alloc_count
		);

		if (readItems <= 0) {
			printf("ERROR reading the temporary file\n");
			break;
		}

		totalReadItems += readItems;

		for (readIndex = 0; readIndex < readItems; readIndex++) {
			pFileRec1 = &pUpdateCommon->pFileBuffer1[readIndex];

			if (pFileRec1->gsc_bin_index >= minGscBin && pFileRec1->gsc_bin_index <= maxGscBin) {
				pFileRec2 = &pUpdateCommon->pFileBuffer2[curHandledCount];

				pFileRec2->galaxyflag = PLATE_FLAG;
				pFileRec2->galaxyversion = PLATE_VERSION;
				pFileRec2->limiting_mag_local = pFileRec1->limiting_mag_local;
				pFileRec2->geoJulianDate = pFileRec1->geoJulianDate;
				pFileRec2->seriesId = pFileRec1->seriesId;
				pFileRec2->plateNumber = pFileRec1->plateNumber;
				pFileRec2->mosaicNumber = pFileRec1->mosaicNumber;
				pFileRec2->solutionNumber = pFileRec1->solutionNumber;
				pFileRec2->versionId = pFileRec1->versionId;
				pFileRec2->gsc_bin_index = pFileRec1->gsc_bin_index;
				pFileRec2->unused = 0;

				curHandledCount++;
				if (curHandledCount > entryCount) {
					printf("ERROR: curHandledCount %d exceeds entryCount %d\n", curHandledCount, entryCount);
					exit(-1);
				}
			}
		}
	}

	if (curHandledCount != entryCount) {
		printf("ERROR: curHandledCount %d does not equal entryCount %d\n", curHandledCount, entryCount);
		exit(-1);
	}

	if (curHandledCount > 0) {
		/* Sort these records in increasing GSC bin */
		qsort(pUpdateCommon->pFileBuffer2, curHandledCount, sizeof(FILELIMITINGREC), GscBinCompare);

		/* Now reformat these records into plate limiting magnitude records */
		for (curHandledIndex = 0; curHandledIndex < curHandledCount; curHandledIndex++) {
			pFileRec2 = &pUpdateCommon->pFileBuffer2[curHandledIndex];
			pPlateRec1 = (PPLATELIMITINGREC) &pUpdateCommon->pFileBuffer1[curHandledIndex];

			pPlateRec1->galaxyflag = PLATE_FLAG;
			pPlateRec1->galaxyversion = PLATE_VERSION;
			pPlateRec1->limiting_mag_local = pFileRec2->limiting_mag_local;
			pPlateRec1->geoJulianDate = pFileRec2->geoJulianDate;
			pPlateRec1->seriesId = pFileRec2->seriesId;
			pPlateRec1->plateNumber = pFileRec2->plateNumber;
			pPlateRec1->mosaicNumber = pFileRec2->mosaicNumber;
			pPlateRec1->solutionNumber = pFileRec2->solutionNumber;
			pPlateRec1->versionId = pFileRec2->versionId;
			pPlateRec1->unused = 0;
		}
	}

	/* Now write out all of the records */
	curHandledIndex = 0;
	pFileRec1 = pUpdateCommon->pFileBuffer2; /* Start of read */
	pFileRec2 = pUpdateCommon->pFileBuffer2; /* End of read */
	pPlateRec1 = (PPLATELIMITINGREC)pUpdateCommon->pFileBuffer1; /* Start of write */
	pPlateRec2 = (PPLATELIMITINGREC)pUpdateCommon->pFileBuffer1; /* End of write */

	for (gsc_bin_index = minGscBin; gsc_bin_index <= maxGscBin; gsc_bin_index++) {
		/* First we write out all records for variables and galaxies in the vicinity of the bin */
		totalWrite = 0;

		if (pUpdateCommon->nearbyFlag != 0) {
			if (GetBinCenter(pGscBin, gsc_bin_index, &ra, &dec, "FlushBinCache")) {
				printf("ERROR: GetBinCenter failed for index %d\n", gsc_bin_index);
				exit(-1);
			}

			z = sin(DEGREES_TO_RAD * dec);
			cosx = cos(DEGREES_TO_RAD * dec);
			x = cosx * cos(DEGREES_TO_RAD * ra);
			y = cosx * sin(DEGREES_TO_RAD * ra);
			pt[0] = x;
			pt[1] = y;
			pt[2] = z;

			presults = kd_nearest_range(pGalaxyCommon->ptree, pt, pUpdateCommon->searchRadius * DEGREES_TO_RAD);
			resultSize = kd_res_size(presults);

			if (resultSize > pUpdateCommon->maxGalaxyBinCount) {
				pUpdateCommon->maxGalaxyBinCount = resultSize;
			}

			pUpdateCommon->aveGalaxyBinCount += resultSize;

			while (!kd_res_end(presults)) {
				/* get the data and position of the current result item */
				pGalaxyRec = (PGALAXYREC) kd_res_item(presults, pos);

				/* Write out this record */
				writeItems = Write(pUpdateCommon->outputHandle, pGalaxyRec, sizeof(GALAXYREC), 1);
				if (writeItems != 1) {
					printf("ERROR writing a galaxy/variable record\n");
					exit(-1);
				}

				totalWrite++;
				/* go to the next entry */
				kd_res_next(presults);
			}

			kd_res_free(presults);
			if (totalWrite != resultSize) {
				printf("ERROR galaxy/variable write count %d not equal to search result %d\n", totalWrite, resultSize);
				exit(-1);
			}
		}

		writeItems = 0;

		if (pUpdateCommon->platesFlag != 0) {
			numObjects = 0;

			/* Now figure out how many objects in our buffer apply to this bin index */
			while (pFileRec2->gsc_bin_index == gsc_bin_index) {
				numObjects++;
				pFileRec2++;
				pPlateRec2++;
			}

			/* Now write out the limiting magnitude cache */
			if (numObjects > 0) {
				if (numObjects > pUpdateCommon->maxPlateBinCount) {
					pUpdateCommon->maxPlateBinCount = numObjects;
				}

				pUpdateCommon->avePlateBinCount += numObjects;

				qsort((void *)pPlateRec1, numObjects, sizeof(PLATELIMITINGREC), PlaterecCompare);

				writeItems = Write(
					pUpdateCommon->outputHandle,
					pPlateRec1,
					sizeof(PLATELIMITINGREC),
					numObjects
				);

				if (writeItems != numObjects) {
					printf("ERROR writing %ld of %d items to the limiting magnitude table\n", writeItems, numObjects);
					exit(-1);
				}

				totalWriteItems += writeItems;

				if (totalWriteItems > entryCount) {
					printf("ERROR: totalWriteItems %d is greater than entryCount %d\n", totalWriteItems, entryCount);
					exit(-1);
				}

				pFileRec1 = pFileRec2;
				pPlateRec1 = pPlateRec2;
			}
		}

		pStarIndex->offset = pUpdateCommon->platerec_bytes;
		pStarIndex->binNumber = gsc_bin_index;
		pUpdateCommon->platerec_bytes += (writeItems * sizeof(PLATELIMITINGREC)) + (totalWrite * sizeof(GALAXYREC));
		pStarIndex->numStars = (writeItems * sizeof(PLATELIMITINGREC)) + (totalWrite * sizeof(GALAXYREC));

		writeItems = Write(pUpdateCommon->indexHandle, pStarIndex, sizeof(STARINDEX), 1);
		if (writeItems != 1) {
			fprintf(stderr,"ERROR writing the index file\n");
			exit(-1);
		}

		pUpdateCommon->lastBinFlushed = gsc_bin_index;
	} /* gsc index loop */

	if (totalWriteItems != entryCount) {
		printf("ERROR: totalWriteItems %d does not equal entryCount %d\n", totalWriteItems, entryCount);
		exit(-1);
	}

	pUpdateCommon->fileHandledCount += totalWriteItems;
}


static void
CheckBin(
	struct WorldCoor *wcs,
	int gsc_bin_index,
	int *goodBin,
	int *partialBin,
	double *xpix,
	double *ypix,
	PPLATEENTRY pPlate,
	double leftMarginD,
	double rightMarginD,
	double topMarginD,
	double bottomMarginD
) {
	PBININDEX pBinIndex;
	int decBin;
	int raBin;
	double ra;
	double dec;
	double deltadec;
	double deltara;
	double xpix2;
	double ypix2;
	int offscl;
	int goodCount = 0;

	pBinIndex = GetSubbins(pGscBin, gsc_bin_index, &raBin, &decBin, pPlate->Plate);
	*goodBin = 0;
	*partialBin = 0;

	GetBinCenter(pGscBin, gsc_bin_index, &ra, &dec, pPlate->Plate);
	wcs2pix(wcs, ra, dec, xpix, ypix, &offscl);

	if (offscl == 0 && *xpix > leftMarginD && *xpix < rightMarginD && *ypix > bottomMarginD && *ypix < topMarginD) {
		goodCount++;
	}

	deltadec = pGscBin->bin_size / 2;
	deltara =  0.5 * (360.00 / pBinIndex->numBins);

	/* Now do the bin corners */

	wcs2pix(wcs, ra + deltara, dec - deltadec, &xpix2, &ypix2, &offscl);
	if (offscl == 0 && xpix2 > leftMarginD && xpix2 < rightMarginD && ypix2 > bottomMarginD && ypix2 < topMarginD) {
		goodCount++;
	}

	wcs2pix(wcs, ra - deltara, dec - deltadec, &xpix2, &ypix2, &offscl);
	if (offscl == 0 && xpix2 > leftMarginD && xpix2 < rightMarginD && ypix2 > bottomMarginD && ypix2 < topMarginD) {
		goodCount++;
	}

	wcs2pix(wcs, ra + deltara, dec + deltadec, &xpix2, &ypix2, &offscl);
	if (offscl == 0 && xpix2 > leftMarginD && xpix2 < rightMarginD && ypix2 > bottomMarginD && ypix2 < topMarginD) {
		goodCount++;
	}

	wcs2pix(wcs, ra - deltara, dec + deltadec, &xpix2, &ypix2, &offscl);
	if (offscl == 0 && xpix2 > leftMarginD && xpix2 < rightMarginD && ypix2 > bottomMarginD && ypix2 < topMarginD) {
		goodCount++;
	}

	if (goodCount == 5) {
		*goodBin = 1;
	} else if (goodCount > 0) {
		*partialBin = 1;
	}
}


static void
GetPlateInformation(
	PPHOTGLOBAL pPhotGlobal,
	PUPDATECOMMON pUpdateCommon,
	PPLATEENTRY pPlate,
	int *releaseFieldPlates
) {
	int exposureNumber;
	int numExposures;
	double singleJulianDate;
	double singleTolerance;
	double allJulianDate;
	double allTolerance;
	PMOSAIC pMosaic = &pPlate->mosaic;
	double ra;
	double dec;
	struct WorldCoor *wcs;
	double xpix;
	double ypix;
	double leftMarginD;
	double rightMarginD;
	double bottomMarginD;
	double topMarginD;
	double cd[4];
	double xctr;
	double yctr;
	int decBin;
	int raBin;
	int gsc_bin_index;
	int goodBin;
	int partialBin;
	PGSCLIMITINGBINREC pGscLimitingBin;
	PGSCLIMITINGBINREC pGscLimitingBin2;
	int ix;
	int iy;
	int index;
	int gsc_bin_list[MAX_ADJACENT_BINS];
	int gsc_bin_count;
	int new_gsc_bin;
	int minReleaseField = RELEASE_FIELD_MAX;
	int fieldRelease[RELEASE_FIELD_MAX];
	int releaseIndex;

	memset(fieldRelease, 0, sizeof(fieldRelease));

	if (pPlate->solutionNumber == 0) {
		sprintf(pPlate->Plate, "%s%05d_%02d", pPlate->series, pPlate->plateNumber, pPlate->mosaicNumber);
	} else {
		sprintf(pPlate->Plate, "%s%05d_%02d_s%d", pPlate->series, pPlate->plateNumber, pPlate->mosaicNumber, pPlate->solutionNumber);
	}

	if (GetSolutionJulianDate(
			pUpdateCommon->pConnection,
			pPlate->series,
			pPlate->plateNumber,
			pPlate->mosaicNumber,
			pPlate->solutionNumber,
			&exposureNumber,
			&numExposures,
			&singleJulianDate,
			&singleTolerance,
			&allJulianDate,
			&allTolerance
	)) {
		printf("ERROR: Failed to get exposure record for %s\n", pPlate->Plate);
		pPlate->errorFlag = 1;
		return;
	}

	pPlate->geoJulianDate = singleJulianDate;

	if (strstr(pMosaic->ctype1,"DEC")) {
		char tmpPtr[MAX_CTYPE_STRING];
		double dtmp;

		strcpy(tmpPtr, pMosaic->ctype2);
		strcpy(pMosaic->ctype2, pMosaic->ctype1);
		strcpy(pMosaic->ctype1, tmpPtr);
		dtmp = pMosaic->crval1;
		pMosaic->crval1 = pMosaic->crval2;
		pMosaic->crval2 = dtmp;
		cd[0] = pMosaic->cd2_1;
		cd[1] = pMosaic->cd2_2;
		cd[2] = pMosaic->cd1_1;
		cd[3] = pMosaic->cd1_2;
	} else {
		cd[0] = pMosaic->cd1_1;
		cd[1] = pMosaic->cd1_2;
		cd[2] = pMosaic->cd2_1;
		cd[3] = pMosaic->cd2_2;
	}

	wcs = wcskinit(
		pMosaic->naxis1,
		pMosaic->naxis2,
		pMosaic->ctype1,
		pMosaic->ctype2,
		pMosaic->crpix1,
		pMosaic->crpix2,
		pMosaic->crval1,
		pMosaic->crval2,
		cd,
		0,  /* cdelt1 */
		0,  /* cdelt2 */
		0,  /* crota */
		2000, /* equinox */
		0   /* epoch */
	);

	if (wcs == NULL) {
		printf("ERROR: failed to get the WCS for plate %s\n", pPlate->Plate);
		pPlate->errorFlag = 1;
		return;
	}

	xctr = 0.5 + (0.5 * pMosaic->naxis1);
	yctr = 0.5 + (0.5 * pMosaic->naxis2);
	leftMarginD = 1.0 * pPlate->leftMargin;
	rightMarginD = 1.0 * (pMosaic->naxis1 - pPlate->rightMargin);
	bottomMarginD = 1.0 * pPlate->bottomMargin;
	topMarginD = 1.0 * (pMosaic->naxis2 - pPlate->topMargin);

	/* Start with the center bin and put it on the queue */
	pix2wcs(wcs, xctr, yctr, &ra, &dec);

	/* Now begin our search for all of the bins covered by the plate */
	memset(pUpdateCommon->checkedTable, 0, pGscBin->total_gsc_bins * sizeof(char));
	gsc_bin_index = GetGSCBin(pGscBin, ra, dec, &decBin, &raBin, pPlate->Plate);
	pUpdateCommon->workQueueSize = 1;
	pUpdateCommon->checkedTable[gsc_bin_index] = 1;
	pUpdateCommon->binWorkQueue[0] = gsc_bin_index;

	if (pPlate->gscLimitingBinBuffer != NULL || pPlate->gscLimitingBinCount != 0 || pPlate->gscLimitingBinAlloc != 0) {
		printf("ERROR: update_limiting gscLimitingBin error\n");
		exit(-1);
	}

	pPlate->gscLimitingBinCount = 0;
	pPlate->gscLimitingBinAlloc = 0;
	pPlate->gscLimitingBinBuffer = NULL;
	pPlate->minBin = gsc_bin_index;
	pPlate->maxBin = gsc_bin_index;

	while (pUpdateCommon->workQueueSize != 0) {
		pUpdateCommon->workQueueSize--;
		gsc_bin_index = pUpdateCommon->binWorkQueue[pUpdateCommon->workQueueSize];

		/* Flag this bin as being checked */

		if (pUpdateCommon->checkedTable[gsc_bin_index] != 1) {
			printf("ERROR: checkTable consistency problem with gsc bin %d\n", gsc_bin_index);
			exit(-1);
		}

		pUpdateCommon->checkedTable[gsc_bin_index] = 2;
		CheckBin(wcs, gsc_bin_index, &goodBin, &partialBin, &xpix, &ypix, pPlate, leftMarginD, rightMarginD, topMarginD, bottomMarginD);

		if (goodBin) {
			pPlate->totalGoodBins++;
			fieldRelease[(int) pUpdateCommon->releaseMappingTable[gsc_bin_index]] = 1;

			if (minReleaseField > pUpdateCommon->releaseMappingTable[gsc_bin_index]) {
				minReleaseField = pUpdateCommon->releaseMappingTable[gsc_bin_index];
			}

			if (pPlate->minBin > gsc_bin_index) {
				pPlate->minBin = gsc_bin_index;
			}

			if (pPlate->maxBin < gsc_bin_index) {
				pPlate->maxBin = gsc_bin_index;
			}

			/* Now store this bin */
			if (pPlate->gscLimitingBinCount + 2 >= pPlate->gscLimitingBinAlloc) {
				PGSCLIMITINGBINREC temp_gscLimitingBinBuffer;

				pPlate->gscLimitingBinAlloc += (pPlate->gscLimitingBinCount + 2 + GSCBINREC_ALLOC_COUNT);
				temp_gscLimitingBinBuffer = (PGSCLIMITINGBINREC) realloc(pPlate->gscLimitingBinBuffer, pPlate->gscLimitingBinAlloc * sizeof(GSCLIMITINGBINREC));
				pUpdateCommon->currentMemory += GSCBINREC_ALLOC_COUNT * sizeof(GSCLIMITINGBINREC);

				if (pUpdateCommon->currentMemory > pUpdateCommon->maxMemory) {
					pUpdateCommon->maxMemory = pUpdateCommon->currentMemory;
				}

				if (temp_gscLimitingBinBuffer == NULL) {
					printf("ERROR: failed to reallocate gsc limiting bin table of size %d\n", pPlate->gscLimitingBinAlloc);
				}

				pPlate->gscLimitingBinBuffer = temp_gscLimitingBinBuffer;
				temp_gscLimitingBinBuffer = NULL;
			}

			pGscLimitingBin = &pPlate->gscLimitingBinBuffer[pPlate->gscLimitingBinCount];
			memset(pGscLimitingBin, 0, sizeof(GSCLIMITINGBINREC));
			pGscLimitingBin->limiting_mag_local = -1.0;
			pGscLimitingBin->gsc_bin_index = gsc_bin_index;

			pPlate->nx = XDmagBins(pMosaic->naxis1, pMosaic->naxis2);
			pPlate->ny = YDmagBins(pMosaic->naxis1, pMosaic->naxis2);

			ix = (xpix * pPlate->nx) / (1.0 * pMosaic->naxis1);
			iy = (ypix * pPlate->ny) / (1.0 * pMosaic->naxis2);

			if (ix < 0) {
				ix = 0;
			}

			if (ix >= pPlate->nx) {
				ix = pPlate->nx - 1;
			}

			if (iy < 0) {
				iy = 0;
			}

			if (iy >= pPlate->ny) {
				iy = pPlate->ny - 1;
			}

			pGscLimitingBin->local_bin_index = ix + pPlate->nx * iy;
			pPlate->gscLimitingBinCount++;
		}

		if (partialBin != 0) {
			pPlate->totalPartialBins++;
		}

		if (goodBin != 0 || partialBin != 0) {
			/* We have a hit, so put all of the surrounding bins on the work queue if they haven't yet been checked */
			FindAdjacentBins(pGscBin, gsc_bin_index, gsc_bin_list, &gsc_bin_count);

			for (index = 0; index < gsc_bin_count; index++) {
				new_gsc_bin = gsc_bin_list[index];

				if (pUpdateCommon->checkedTable[new_gsc_bin] == 0) {
					/* Not yet checked -- put it on the list */
					pUpdateCommon->checkedTable[new_gsc_bin] = 1;
					pUpdateCommon->binWorkQueue[pUpdateCommon->workQueueSize] = new_gsc_bin;
					pUpdateCommon->workQueueSize++;

					if (pUpdateCommon->workQueueSize >= pGscBin->total_gsc_bins) {
						printf("ERROR: too many bins on the work queue\n");
						exit(-1);
					}
				}
			}
		}
	}

	if (pPlate->solutionNumber == 0) {
		if (minReleaseField < RELEASE_FIELD_MAX) {
			releaseFieldPlates[RELEASE_FIELD_MAX]++;
			releaseFieldPlates[minReleaseField]++;
		}

		for (releaseIndex = 0; releaseIndex < RELEASE_FIELD_MAX; releaseIndex++) {
			if (fieldRelease[releaseIndex] != 0) {
				pUpdateCommon->countRelease[releaseIndex]++;
			}
		}
	}

	/* Now sort these legal bins in increasing gsc_bin_index */

	if (pPlate->gscLimitingBinCount > 1) {
		if (pPlate->gscLimitingBinCount >= pPlate->gscLimitingBinAlloc) {
			printf("ERROR: line %d gscLimitingBinCount %d >= gscLimitingBinAlloc %d\n", __LINE__, pPlate->gscLimitingBinCount, pPlate->gscLimitingBinAlloc);
			exit(-1);
		}

		for (index = 0; index < (pPlate->gscLimitingBinCount); index++) {
			pGscLimitingBin = &pPlate->gscLimitingBinBuffer[index];

			if (pGscLimitingBin->local_bin_index >= TOTAL_DMAGBINS_NORMAL) {
				printf("ERROR: line %d local_bin_index %d exceeds %d\n", __LINE__, pGscLimitingBin->local_bin_index, TOTAL_DMAGBINS_NORMAL);
				exit(-1);
			}
		}

		qsort((void *)pPlate->gscLimitingBinBuffer, pPlate->gscLimitingBinCount, sizeof(GSCLIMITINGBINREC), BinrecCompare);

		if (pPlate->gscLimitingBinCount >= pPlate->gscLimitingBinAlloc) {
			printf("ERROR: line %d gscLimitingBinCount %d >= gscLimitingBinAlloc %d\n", __LINE__, pPlate->gscLimitingBinCount, pPlate->gscLimitingBinAlloc);
			exit(-1);
		}

		for (index = 0; index < pPlate->gscLimitingBinCount - 1; index++) {
			pGscLimitingBin = &pPlate->gscLimitingBinBuffer[index];

			if (pGscLimitingBin->local_bin_index >= TOTAL_DMAGBINS_NORMAL) {
				printf("ERROR: line %d local_bin_index %d exceeds %d\n", __LINE__, pGscLimitingBin->local_bin_index, TOTAL_DMAGBINS_NORMAL);
				exit(-1);
			}

			pGscLimitingBin2 = &pPlate->gscLimitingBinBuffer[index+1];

			if (pGscLimitingBin2->local_bin_index >= TOTAL_DMAGBINS_NORMAL) {
				printf("ERROR: line %d local_bin_index %d exceeds %d\n", __LINE__, pGscLimitingBin2->local_bin_index, TOTAL_DMAGBINS_NORMAL);
				exit(-1);
			}

			if (pGscLimitingBin->gsc_bin_index == pGscLimitingBin2->gsc_bin_index) {
				printf("ERROR: duplicate gsc bin index %d\n", pGscLimitingBin->gsc_bin_index);
				exit(-1);
			}
		}
	}

	if (pUpdateCommon->verbose) {
		printf(
			"plate %5s%05d totalGoodBins %5d totalPartialBins %5d ratio %f\n",
			pPlate->series,
			pPlate->plateNumber,
			pPlate->totalGoodBins,
			pPlate->totalPartialBins,
			(1.0 * pPlate->totalGoodBins) / (1.0 * (pPlate->totalGoodBins + pPlate->totalPartialBins))
		);
	}

	wcsfree(wcs);
}


static int
LimitingMagCompare(const void *first, const void *second)
{
	double gsc_bin_indexFirst = ((PLIMITINGMAG) first)->gsc_bin_index;
	double gsc_bin_indexSecond = ((PLIMITINGMAG) second)->gsc_bin_index;

	if (gsc_bin_indexFirst > gsc_bin_indexSecond) {
		return 1;
	} else if (gsc_bin_indexFirst < gsc_bin_indexSecond) {
		return -1;
	}

	return 0;
}


static void
GetReleaseRegions(PUPDATECOMMON pUpdateCommon)
{
	int gsc_bin_index;
	double ra;
	double dec;
	int releaseField;

	pUpdateCommon->binWorkQueue = (int *) calloc(pGscBin->total_gsc_bins, sizeof(int));
	if (pUpdateCommon->binWorkQueue == NULL) {
		printf("ERROR: failed to allocate binWorkQueue of size %d\n", pGscBin->total_gsc_bins);
		exit(-1);
	}

	pUpdateCommon->checkedTable = (char *) calloc(pGscBin->total_gsc_bins, sizeof(char));
	if (pUpdateCommon->checkedTable == NULL) {
		printf("ERROR: failed to allocate checkedTable of size %d\n", pGscBin->total_gsc_bins);
		exit(-1);
	}

	pUpdateCommon->releaseMappingTable = (char *) calloc(pGscBin->total_gsc_bins, sizeof(char));
	if (pUpdateCommon->releaseMappingTable == NULL) {
		printf("ERROR: failed to allocate releaseMappingTable of size %d\n", pGscBin->total_gsc_bins);
		exit(-1);
	}

	pUpdateCommon->currentMemory += pGscBin->total_gsc_bins + (sizeof(int *) + sizeof(char *) + sizeof(char*));
	if (pUpdateCommon->currentMemory > pUpdateCommon->maxMemory) {
		pUpdateCommon->maxMemory = pUpdateCommon->currentMemory;
	}

	/* Now get the release regions for every GSC bin. */
	for (gsc_bin_index = 0; gsc_bin_index < pGscBin->total_gsc_bins; gsc_bin_index++) {
		if (GetBinCenter(pGscBin, gsc_bin_index, &ra, &dec, "update_limiting")) {
			printf("ERROR: GetBinCenter failed for index %d\n", gsc_bin_index);
			exit(-1);
		}

		CheckAuthorization(NULL, ra, dec, &releaseField);
		if (releaseField >= 0 && releaseField < RELEASE_FIELD_MAX) {
			pUpdateCommon->releaseMappingTable[gsc_bin_index] = releaseField;
		}
	}
}


int main(int argc, char *argv[])
{
	UPDATECOMMON updateCommon;
	PUPDATECOMMON pUpdateCommon = &updateCommon;
	GALAXYCOMMON galaxycommon;
	PGALAXYCOMMON pGalaxyCommon = &galaxycommon;
	PPHOTGLOBAL pPhotGlobal = &pUpdateCommon->photGlobal;

	char *argstr;
	int errorFlag = 0;
	int maxString = 0;
	char cmdchar;
	char logfilename[MAX_BUFFER];
	char qualifier[MAX_BUFFER];
	char outputname[MAX_BUFFER];
	char indexname[MAX_BUFFER];
	char temporaryname[MAX_BUFFER];
	char* catalogdir;
	char galaxy_name[MAX_BUFFER];

	int plateIndex;
	char timestr[100];
	struct tm *ptr;
	int gotAnswer;

	MYSQL my_connection;
	MYSQL my_phot_connection;
	PPHOT_SPATIAL_BIN pSpatial_bin;
	int spatial_bin_nrecs;
	int good_bin_nrecs;
	PPLATEENTRY pPlate;
	PPLATEENTRY pPlate2;
	char queryString[MAX_QUERY_STRING];
	int plateCount = 0;
	int solutionCount = 0;
	int solutionNumber;
	int solutionNumberFix = 0;
	char solutionNumberString[20];
	char dasch_phot_magnitudes[800];
	MOSAIC mosaicInfo;
	PMOSAIC pMosaic = &mosaicInfo;
	PHOT_SPATIAL_BIN spatial_bin_table[MAX_SPATIAL_BINS+1];
	PPHOT_SPATIAL_BIN spatial_bin_table2 = NULL;
	PPHOT_SPATIAL_BIN pSpatialBinEntry;
	int curPlateCount;
	int spatial_bin_index;
	int spatial_bin;
	LOCALBIN local_bin_table[TOTAL_DMAGBINS_NORMAL];
	PLOCALBIN pLocalBin;
	LIMITINGMAG limiting_mag_table[TOTAL_DMAGBINS_NORMAL];
	LIMITINGMAG limiting_mag_table_unsorted[TOTAL_DMAGBINS_NORMAL];
	PLIMITINGMAG pLimitingMag;
	double cd[4];
	struct WorldCoor *wcs;
	int local_bin_index;
	int lb_index;
	int ix;
	int iy;
	double xctr;
	double yctr;
	double cra;
	double cdec;
	int raBin;
	int decBin;
	MAGDEPLIMITS magdepLimits;
	PMAGDEPLIMITS pMagdepLimits = &magdepLimits;
	PMAGDEPCORRECTION pMagdepTable = NULL;
	int goodPlateCount = 0;
	int gsc_bin_index = 0;
	char *dotPtr;
	char *slashPtr;
	time_t curTime;
	int nvals;
	int binningFactor = 2;
	int missingBinCount = 0;
	int releaseFieldPlates[RELEASE_FIELD_MAX+1];
	int releaseField;
	int minGscBin = 0;
	size_t writtenEntryCount = 0;
	int maxGscBin = 0;
	int currentEntryCount = 0;
	int errorPrintFlag = 0;

	memset(releaseFieldPlates, 0, sizeof(releaseFieldPlates));
	memset(pUpdateCommon, 0, sizeof(UPDATECOMMON));
	memset(pGalaxyCommon, 0, sizeof(GALAXYCOMMON));

	pUpdateCommon->pConnection = &my_connection;
	pUpdateCommon->pPhotConnection = &my_phot_connection;
	pUpdateCommon->catalogString[0] = 0;
	pUpdateCommon->lastBinFlushed = -1;
	pUpdateCommon->platesFlag = 1;
	pUpdateCommon->nearbyFlag = 1;

	/* Take the galaxy radius cutoff and add to it half the bin diagonal and add
	 * an extra half bin for objects near the edge of the bin */
	pUpdateCommon->searchRadius = GALAXY_SEARCH_RADIUS + (0.5 * pGscBin->bin_size * sqrt(2.0)) + (0.5 * pGscBin->bin_size);

	logfilename[0] = 0;
	qualifier[0] = 0;
	outputname[0] = 0;
	indexname[0] = 0;
	temporaryname[0] = 0;

	/* Loop through the arguments */

	for (argv++; --argc > 0; argv++) {
		argstr = *argv;

		if (argstr[0] != '-') {
			errorFlag = 1;
			printf("ERROR: unqualified argument %s argc: %d\n", argstr, argc);
		} else {
			while ((cmdchar = *++argstr) != 0) {
				switch(cmdchar) {
				case 'l': /* log file name */
				case 'L':
					argc--;
					if (argc < 1) {
						printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(logfilename,*++argv,MAX_BUFFER-2);
						if (strlen(logfilename) >= MAX_BUFFER-3) {
							printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
						}
					}
					break;

				case 'b': /* binning factor */
				case 'B':
					argc--;
					if (argc < 1) {
						printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						nvals = sscanf(*++argv,"%d",&binningFactor);
						if (nvals != 1) {
							printf("ERROR: Unable to decode the binningFactor %s\n",*argv);
							errorFlag = 1;
						} else {
							switch(binningFactor) {
							case 1:
								pGscBin = &gscBin01;
								break;
							case 2:
								pGscBin = &gscBin02;
								break;
							case 4:
								pGscBin = &gscBin04;
								break;
							case 8:
								pGscBin = &gscBin08;
								break;
							case 16:
								pGscBin = &gscBin16;
								break;
							case 32:
								pGscBin = &gscBin32;
								break;
							case 64:
								pGscBin = &gscBin64;
								break;
							default:
								printf("ERROR: binning factor must be one of 1,2,4,8,16,32, or 64\n");
								errorFlag = 1;
							}
						}
					}
					break;

				case 'o': /* binary output name */
				case 'O':
					argc--;
					if (argc < 1) {
						printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(outputname,*++argv,MAX_BUFFER-2);
						if (strlen(outputname) >= MAX_BUFFER-3) {
							printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
						}
					}
					break;

				case 'q': /* Catalog and file name qualifier */
				case 'Q':
					argc--;
					if (argc < 1) {
						printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						pUpdateCommon->catalogNumber = GetCatalogNumber(*++argv);

						if (pUpdateCommon->catalogNumber < 0) {
							printf("ERROR: Illegal catalog name %s\n",*argv);
							errorFlag = 1;
						} else {
							sprintf(pUpdateCommon->catalogString,"%d",pUpdateCommon->catalogNumber);
							sprintf(qualifier,"_%s",catalogText[pUpdateCommon->catalogNumber]);
						}
					}
					break;

				case 'v': /* verbose */
				case 'V':
					pUpdateCommon->verbose = 1;
					break;

				case 'n': /* nearby galaxies only */
				case 'N':
					pUpdateCommon->platesFlag = 0;
					break;

				case 'p': /* plates only */
				case 'P':
					pUpdateCommon->nearbyFlag = 0;
					break;

				default:
					printf("ERROR: * illegal command -%c-\n",cmdchar);
					errorFlag = 1;
				}
			}
		}
	}

	catalogdir = getenv("DASCH_CATALOG");
	if (catalogdir == NULL) {
		fprintf(stderr,"DASCH_CATALOG is not defined\n");
		return -1;
	}

	strcpy(galaxy_name, catalogdir);
	slashPtr = strrchr(galaxy_name,'/');
	if (slashPtr != NULL) {
		slashPtr++;
	} else {
		slashPtr = galaxy_name;
	}
	*slashPtr = 0;
	strcat(galaxy_name, "galaxy.dat");

	// Connect to databases. These will abort the process if any unsolvable
	// problems occur.
	dasch_init_scandb(pUpdateCommon->pConnection);
	dasch_init_photdb(pUpdateCommon->pPhotConnection);

	InitSeriesTable(pUpdateCommon->pConnection, pUpdateCommon->pPhotConnection);

	if (GetPhotometryGlobal(pUpdateCommon->pPhotConnection, pPhotGlobal) != 1) {
		printf("ERROR: failed to get the global photometry table\n");
		exit(-1);
	}

	if (pPhotGlobal->magnitudeFile == PHOT_MAGNITUDEFILE_YES) {
		pUpdateCommon->magnitudeFileFlag = 1;
	}

	if (pPhotGlobal->solutionNumber == PHOT_SOLUTIONNUMBER_YES) {
		solutionNumberFix = 1;
		strcpy(solutionNumberString,"solutionNumber");
	} else {
		strcpy(solutionNumberString,"exposureNumber");
	}

	/* Validate arguments */

	if (logfilename[0] == 0) {
		printf("ERROR: No log filename was specified\n");
		errorFlag = 1;
	} else {
		pUpdateCommon->logHandle = fopen(logfilename,"a+t");
		if (pUpdateCommon->logHandle == NULL) {
			errorFlag = 1;
			printf("ERROR: Failed to open the output file %s\n",logfilename);
		} else if (pUpdateCommon->verbose) {
			printf("Output file %s\n",logfilename);
		}
	}

	if (outputname[0] == 0) {
		printf("ERROR: no output name specified\n");
		errorFlag = 1;
	} else {
		dotPtr = strrchr(outputname, '.');
		if (dotPtr != NULL) {
			*dotPtr = 0;
		}
		strcat(outputname, ".dat");

		strcpy(indexname, outputname);
		dotPtr = strrchr(indexname, '.');
		if (dotPtr != NULL) {
			*dotPtr = 0;
		}
		strcat(indexname, ".idx");

		strcpy(temporaryname, outputname);
		dotPtr = strrchr(temporaryname, '.');
		if (dotPtr != NULL) {
			*dotPtr = 0;
		}
		strcat(temporaryname, ".tmp");
	}

	if (errorFlag) {
		printf("Usage: update_limiting options\n");
		printf("  options: -v verbose\n");
		printf("           -l <log file>\n");
		printf("           -o <output file>\n");
		printf("           -b <binning factor> (default 2)\n");
		printf("           -n nearby galaxies only\n");
		printf("           -p plate limiting magnitudes only\n");
		printf("           -q input catalog and filename qualifer\n");
		return -1;
	}

	pUpdateCommon->outputHandle = Open(outputname,"w");
	if (pUpdateCommon->outputHandle == NULL) {
		printf("Could not open file %s\n",outputname);
		exit(-1);
	}

	pUpdateCommon->indexHandle = Open(indexname,"w");
	if (pUpdateCommon->indexHandle == NULL) {
		printf("Could not open file %s\n",indexname);
		exit(-1);
	}

	pUpdateCommon->temporaryHandle = Open(temporaryname,"w");
	if (pUpdateCommon->temporaryHandle == NULL) {
		printf("Could not open file %s\n",temporaryname);
		exit(-1);
	}

	sprintf(dasch_phot_magnitudes, "DASCH_PHOT_MAGNITUDES%s", pUpdateCommon->catalogString);

	pUpdateCommon->photfilebase = getenv(dasch_phot_magnitudes);
	if (pUpdateCommon->photfilebase == NULL) {
		printf("%s is not defined\n", dasch_phot_magnitudes);
		return -1;
	}

	assert(sizeof(PLATELIMITINGREC) == sizeof(FILELIMITINGREC));

	printf(
		"update_limiting of %s %s, entries %d, max string %d catalog %s filebase %s nearbyFlag %d\n",
		__DATE__,
		__TIME__,
		pUpdateCommon->numPlates,
		maxString,
		catalogText[pUpdateCommon->catalogNumber],
		pUpdateCommon->photfilebase,
		pUpdateCommon->nearbyFlag
	);

	fprintf(
		pUpdateCommon->logHandle,
		"update_limiting of %s %s, entries %d, max string %d catalog %s filebase %s nearbyFlag %d\n",
		__DATE__,
		__TIME__,
		pUpdateCommon->numPlates,
		maxString,
		catalogText[pUpdateCommon->catalogNumber],
		pUpdateCommon->photfilebase,
		pUpdateCommon->nearbyFlag
	);

	printf(
		"sizeof STARINDEX %zu sizeof GALAXYREC %zu PLATELIMITINGREC %zu FILELIMITINGREC %zu\n",
		sizeof(STARINDEX),
		sizeof(GALAXYREC),
		sizeof(PLATELIMITINGREC),
		sizeof(FILELIMITINGREC)
	);

	time(&pUpdateCommon->startTime);

	ptr = localtime(&pUpdateCommon->startTime);
	strftime(timestr, 25, "%Y-%m-%dT%H-%M-%S", ptr);

	printf("update_limiting starting at %s\n", timestr);
	fprintf(pUpdateCommon->logHandle, "update_limiting starting at %s\n", timestr);
	fflush(pUpdateCommon->logHandle);

	sprintf(queryString,"update photglobal set keepalive = 'yes';");
	ExecuteQuery(pUpdateCommon->pPhotConnection,queryString);

	if (pUpdateCommon->nearbyFlag != 0) {
		if (PopulateGalaxyTree(pGalaxyCommon, galaxy_name, 0) != 0) {
			exit(-1);
		}
	}

	if (pUpdateCommon->platesFlag != 0) {
		GetReleaseRegions(pUpdateCommon);
		time(&pUpdateCommon->curTime);
		pUpdateCommon->curTime -= pUpdateCommon->startTime;
		printf("Time to get release regions for %d bins is %ld seconds\n", pGscBin->total_gsc_bins, pUpdateCommon->curTime);

		BuildPlatesTable(pUpdateCommon);

		/* First loop through the plate table.  Get the mosaic info so we can sort the table by bin index */
		curPlateCount = pUpdateCommon->numPlates;

		for (plateIndex = 0; plateIndex < pUpdateCommon->numPlates; plateIndex++) {
			time(&pUpdateCommon->plateStartTime);

			pPlate = &pUpdateCommon->pPlateTable[plateIndex];

			if (((plateIndex + 1) % 10000) == 0) {
				time(&curTime);
				curTime -= pUpdateCommon->startTime;
				printf("Line %d At plate %6d %s%05d in %ld seconds\n", __LINE__, plateIndex + 1, pPlate->series, pPlate->plateNumber, curTime);
				gotAnswer = GetPhotometryGlobal(pUpdateCommon->pPhotConnection, pPhotGlobal);

				if (gotAnswer != 1 || pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES) {
					ptr = localtime(&pUpdateCommon->plateStartTime);
					strftime(timestr, 25, "%Y-%m-%dT%H-%M-%S", ptr);

					printf(
						"update_limiting aborting at %s with plates %d  gotAnswer %d keepalive %d\n",
						timestr,
						plateIndex,
						gotAnswer,
						pPhotGlobal->keepalive
					);
					fprintf(
						pUpdateCommon->logHandle,
						"update_limiting aborting at %s with plates %d  gotAnswer %d keepalive %d\n",
						timestr,
						plateIndex,
						gotAnswer,
						pPhotGlobal->keepalive
					);
					break;
				}
			}

			pPlate = &pUpdateCommon->pPlateTable[plateIndex];

			if (SelectBestMosaic(
				pUpdateCommon->pConnection,
				pPlate->series,
				pPlate->plateNumber,
				&pPlate->mosaicNumber,
				&pPlate->rotation
			) == 0) {
				printf("Failed to find a candidate mosaic for plate %s%05d\n", pPlate->series, pPlate->plateNumber);
				pPlate->errorFlag = 1;
				continue;
			}

			for (solutionNumber = 0; solutionNumber < MAX_CATALOG_EXPOSURES; solutionNumber++) {
				if (GetMosaicInfo(
					pUpdateCommon->pConnection,
					pPlate->series,
					pPlate->plateNumber,
					pPlate->mosaicNumber,
					solutionNumber,
					pMosaic
				) != 1) {
					break;
				}

				if (strcmp(pPlate->series,"ac") == 0) {
					pMosaic->cd1_1 = AC_SCALE_FACTOR * pMosaic->cd1_1;
					pMosaic->cd2_2 = AC_SCALE_FACTOR * pMosaic->cd2_2;
					pMosaic->cd1_2 = AC_SCALE_FACTOR * pMosaic->cd1_2;
					pMosaic->cd2_1 = AC_SCALE_FACTOR * pMosaic->cd2_1;
				}

				if (solutionNumber == 0) {
					memcpy(&pPlate->mosaic, pMosaic, sizeof(MOSAIC));
				} else {
					/* Need to allocate a new plate structure for this solution */
					if ((curPlateCount+1) >= pUpdateCommon->plateTableAlloc) {
						pUpdateCommon->plateTableAlloc += PLATE_TABLE_INCREMENT;
						PPLATEENTRY tmpPlateTable = (PPLATEENTRY) realloc(pUpdateCommon->pPlateTable, pUpdateCommon->plateTableAlloc * sizeof(PLATEENTRY));
						pUpdateCommon->currentMemory += PLATE_TABLE_INCREMENT * sizeof(PLATEENTRY);

						if (pUpdateCommon->currentMemory > pUpdateCommon->maxMemory) {
							pUpdateCommon->maxMemory = pUpdateCommon->currentMemory;
						}

						if (tmpPlateTable == NULL) {
							printf("ERROR: failed to reallocate plate table of size %d\n", pUpdateCommon->plateTableAlloc);
						}

						pUpdateCommon->pPlateTable = tmpPlateTable;
						tmpPlateTable = NULL;
					}

					pPlate2 = &pUpdateCommon->pPlateTable[plateIndex];
					pPlate = &pUpdateCommon->pPlateTable[curPlateCount];
					memcpy(pPlate, pPlate2, sizeof(PLATEENTRY));
					pPlate->gscLimitingBinCount = 0;
					pPlate->gscLimitingBinAlloc = 0;
					pPlate->totalGoodBins = 0;
					pPlate->totalPartialBins = 0;
					pPlate->gscLimitingBinBuffer = NULL;
					pPlate->solutionNumber = solutionNumber;
					memcpy(&pPlate->mosaic, pMosaic, sizeof(MOSAIC));
					curPlateCount++;
				}

				solutionCount++;
			}

			plateCount++;
		}

		printf(
			"Found %d additional solutions PLATE_TABLE_INCREMENT is %d\n",
			curPlateCount - pUpdateCommon->numPlates,
			PLATE_TABLE_INCREMENT
		);
		printf(
			"Max memory %lld, current memory %lld\n",
			pUpdateCommon->maxMemory,
			pUpdateCommon->currentMemory
		);
		pUpdateCommon->numPlates = curPlateCount;

		/* Now we can process all of the plates in order */

		pUpdateCommon->file_written_table = (int *) calloc(pGscBin->total_gsc_bins, sizeof(int));
		pUpdateCommon->file_handled_table = (int *) calloc(pGscBin->total_gsc_bins, sizeof(int));
		pUpdateCommon->currentMemory += 2 * sizeof(int);

		if (pUpdateCommon->currentMemory > pUpdateCommon->maxMemory) {
			pUpdateCommon->maxMemory = pUpdateCommon->currentMemory;
		}

		printf(
			"Max memory %lld, current memory %lld\n",
			pUpdateCommon->maxMemory,
			pUpdateCommon->currentMemory
		);

		for (plateIndex = 0; plateIndex < pUpdateCommon->numPlates; plateIndex++) {
			time(&pUpdateCommon->plateStartTime);
			pPlate = &pUpdateCommon->pPlateTable[plateIndex];

			if ((plateIndex + 1) % 1000 == 0) {
				time(&curTime);
				curTime -= pUpdateCommon->startTime;
				printf("Line %d At plate %6d %s%05d in %ld seconds\n", __LINE__, plateIndex + 1, pPlate->series, pPlate->plateNumber, curTime);

				gotAnswer = GetPhotometryGlobal(pUpdateCommon->pPhotConnection, pPhotGlobal);

				if (gotAnswer != 1 || pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES) {
					ptr = localtime(&pUpdateCommon->plateStartTime);
					strftime(timestr, 25, "%Y-%m-%dT%H-%M-%S", ptr);

					printf(
						"update_limiting aborting at %s with plates %d  gotAnswer %d keepalive %d\n",
						timestr,
						plateIndex,
						gotAnswer,
						pPhotGlobal->keepalive
					);
					fprintf(
						pUpdateCommon->logHandle,
						"update_limiting aborting at %s with plates %d  gotAnswer %d keepalive %d\n",
						timestr,
						plateIndex,
						gotAnswer,
						pPhotGlobal->keepalive
					);
					break;
				}
			}

			if (pPlate->errorFlag) {
				continue;
			}

			GetPlateInformation(pPhotGlobal, pUpdateCommon, pPlate, releaseFieldPlates);

			if (pPlate->errorFlag) {
				continue;
			}

			/* Read in the spatial bin table */

			spatial_bin_nrecs = ReadPhotSpatialBin(
				pUpdateCommon->pPhotConnection,
				&spatial_bin_table2,
				pPlate->series,
				pPlate->plateNumber,
				pPlate->solutionNumber,
				0,
				pUpdateCommon->catalogString,
				solutionNumberFix
			);
			if (spatial_bin_nrecs == 0) {
				pPlate->errorFlag = 1;
				continue;
			}

			good_bin_nrecs = 0;

			for (spatial_bin_index = 0; spatial_bin_index < spatial_bin_nrecs; spatial_bin_index++) {
				pSpatialBinEntry = &spatial_bin_table2[spatial_bin_index];
				spatial_bin = pSpatialBinEntry->spatial_bin;

				if (spatial_bin > 0 && spatial_bin <= MAX_SPATIAL_BINS) {
					pSpatial_bin = &spatial_bin_table[spatial_bin];
					if (pSpatialBinEntry->versionId >= pPlate->versionId) {
						memcpy(pSpatial_bin, pSpatialBinEntry, sizeof(PHOT_SPATIAL_BIN));
						pSpatial_bin->spatial_bin = spatial_bin;
						good_bin_nrecs++;
					}
				} else {
					printf("ERROR: Illegal spatial bin %d\n", spatial_bin);
					exit(-1);
				}
			}

			if (spatial_bin_table2 != NULL) {
				free(spatial_bin_table2);
				spatial_bin_table2 = NULL;
			}

			if (good_bin_nrecs == 0) {
				pPlate->errorFlag = 1;
				continue;
			}

			/* Now read in the local bin table */

			if (GetPhotLocalBin(
				pUpdateCommon->pPhotConnection,
				local_bin_table,
				pPlate->nx * pPlate->ny,
				pPlate->series,
				pPlate->plateNumber,
				pPlate->solutionNumber,
				-1,
				pUpdateCommon->catalogString,
				solutionNumberFix,
				0 /* debug mode */
			)) {
				printf(
					"WARNING: Failed to read local bin table for %s%05d  seriesId = %d plateNumber = %d solutionNumber = %d\n",
					pPlate->series,
					pPlate->plateNumber,
					pPlate->seriesId,
					pPlate->plateNumber,
					pPlate->solutionNumber
				);
				pPlate->errorFlag = 1;
				continue;
			}

			pMosaic = &pPlate->mosaic;

			if (strstr(pMosaic->ctype1,"DEC")) {
				char tmpPtr[MAX_CTYPE_STRING];
				double dtmp;

				strcpy(tmpPtr,pMosaic->ctype2);
				strcpy(pMosaic->ctype2,pMosaic->ctype1);
				strcpy(pMosaic->ctype1,tmpPtr);
				dtmp = pMosaic->crval1;
				pMosaic->crval1 = pMosaic->crval2;
				pMosaic->crval2 = dtmp;
				cd[0] = pMosaic->cd2_1;
				cd[1] = pMosaic->cd2_2;
				cd[2] = pMosaic->cd1_1;
				cd[3] = pMosaic->cd1_2;
			} else {
				cd[0] = pMosaic->cd1_1;
				cd[1] = pMosaic->cd1_2;
				cd[2] = pMosaic->cd2_1;
				cd[3] = pMosaic->cd2_2;
			}

			wcs = wcskinit(
				pMosaic->naxis1,
				pMosaic->naxis2,
				pMosaic->ctype1,
				pMosaic->ctype2,
				pMosaic->crpix1,
				pMosaic->crpix2,
				pMosaic->crval1,
				pMosaic->crval2,
				cd,
				0,  /* cdelt1 */
				0,  /* cdelt2 */
				0,  /* crota */
				2000, /* equinox */
				0   /* epoch */
			);

			if (wcs == NULL) {
				printf("ERROR: failed to get the WCS for plate %s\n", pPlate->Plate);
				pPlate->errorFlag = 1;
				continue;
			}

			memset(pMagdepLimits, 0, sizeof(MAGDEPLIMITS));

			if (PhotLoadMagdepCorrections(
				pUpdateCommon->pPhotConnection,
				pPlate->seriesId,
				pPlate->plateNumber,
				pPlate->solutionNumber,
				pUpdateCommon->catalogString,
				pMagdepLimits,
				&pMagdepTable,
				pUpdateCommon->verbose
			) != 0) {
				/* Nothing - missing magdepcal is allowed */
			}

			memset(limiting_mag_table, 0, sizeof(limiting_mag_table));

			/* Now find the local bin index for each of these bins */

			for (local_bin_index = 0; local_bin_index < pPlate->nx * pPlate->ny; local_bin_index++) {
				if (local_bin_index >= TOTAL_DMAGBINS_NORMAL) {
					printf("ERROR: line %d local_bin_index %d exceeds %d\n", __LINE__, local_bin_index, TOTAL_DMAGBINS_NORMAL);
					exit(-1);
				}

				pLocalBin = &local_bin_table[local_bin_index];
				if (pLocalBin->local_bin_index != local_bin_index) {
					printf("ERROR in local bin table indices\n");
					exit(-1);
				}

				if (local_bin_index >= TOTAL_DMAGBINS_NORMAL) {
					printf("ERROR: line %d local_bin_index %d exceeds %d\n", __LINE__, local_bin_index, TOTAL_DMAGBINS_NORMAL);
					exit(-1);
				}

				pLimitingMag = &limiting_mag_table[local_bin_index];
				pLimitingMag->limiting_mag_local = 99.0;
				pLimitingMag->local_bin_index = local_bin_index;
				pLimitingMag->gsc_bin_index = -1;
				if (pLocalBin->rejectFlag != 0) {
					continue;
				}

				if (pMosaic->naxis1 * pMosaic->naxis2 < SMALL_PLATE_PIXELS) {
					spatial_bin = 1;
				} else {
					spatial_bin = binarray[local_bin_index];
				}

				pSpatial_bin = &spatial_bin_table[spatial_bin];
				if (pSpatial_bin->spatial_bin != spatial_bin) {
					continue;
				}

				ix = local_bin_index % pPlate->nx;
				iy = local_bin_index / pPlate->nx;
				xctr = (1.0 * pMosaic->naxis1) * (1.0 * ix + 0.5) / (1.0 * pPlate->nx);
				yctr = (1.0 * pMosaic->naxis2) * (1.0 * iy + 0.5) / (1.0 * pPlate->ny);
				pix2wcs(wcs, xctr, yctr, &cra, &cdec);
				pLimitingMag->gsc_bin_index = GetGSCBin(pGscBin, cra, cdec, &decBin, &raBin, pPlate->Plate);
				pLimitingMag->limiting_mag_local = pSpatial_bin->limiting_mag;

				if (pLocalBin->magcal_local_error < 90.0) {
					pLimitingMag->limiting_mag_local -= pLocalBin->magcor_local;
				}

				if (pMagdepLimits->nrecs > 0) {
					double temp_magdep_rms;
					double temp_magdep_magcor;
					int temp_magdep_bin;

					pLimitingMag->limiting_mag_local -= GetMagdepBinMagcor(
						pMagdepLimits,
						pMagdepTable,
						xctr,
						yctr,
						pLimitingMag->limiting_mag_local,
						&temp_magdep_bin,
						&temp_magdep_rms,
						&temp_magdep_magcor,
						1
					);
				}

				if (pLocalBin->extinction > 0) {
					pLimitingMag->limiting_mag_local -= pLocalBin->extinction;
				}

				pLimitingMag->valid = 1;
			}

			goodPlateCount++;
			wcsfree(wcs);

			if (pMagdepLimits->nrecs > 0) {
				FreeMagdepSubarrays(pMagdepLimits);
			}

			if (pMagdepTable != NULL) {
				free(pMagdepTable);
				pMagdepTable = NULL;
			}

			/* Now sort the results in increasing GSC bin number */

			int last_gsc_bin_index = -1;
			int last_local_bin_index = -1;
			int last_local_bin_index2 = -1;
			int gsc_bin_count = 0;
			double limiting_mag_local_average = 0.0;

			memcpy(limiting_mag_table_unsorted, limiting_mag_table, sizeof(limiting_mag_table_unsorted));
			qsort((void *) limiting_mag_table, pPlate->nx * pPlate->ny, sizeof(LIMITINGMAG), LimitingMagCompare);

			for (lb_index = 0; lb_index < pPlate->nx * pPlate->ny; lb_index++) {
				if (lb_index >= TOTAL_DMAGBINS_NORMAL) {
					printf("ERROR: line %d local_bin_index %d exceeds %d\n",__LINE__,lb_index,TOTAL_DMAGBINS_NORMAL);
					exit(-1);
				}

				pLimitingMag = &limiting_mag_table[lb_index];
				local_bin_index = pLimitingMag->local_bin_index;

				if (pLimitingMag->valid) {
					if (pLimitingMag->gsc_bin_index != last_gsc_bin_index) {
						StoreLimitingMagnitude(
							pUpdateCommon,
							pPlate,
							last_gsc_bin_index,
							gsc_bin_count,
							limiting_mag_local_average,
							last_local_bin_index,
							last_local_bin_index2
						);
						gsc_bin_count = 0;
						limiting_mag_local_average = 0.0;
						last_local_bin_index2 = -1;
					}

					last_gsc_bin_index = pLimitingMag->gsc_bin_index;
					if (last_local_bin_index2 < 0) {
						last_local_bin_index2 = local_bin_index;
					}

					last_local_bin_index = local_bin_index;
					gsc_bin_count++;
					limiting_mag_local_average += pLimitingMag->limiting_mag_local;
				}
			}

			/* Store the last local bin */
			StoreLimitingMagnitude(
				pUpdateCommon,
				pPlate,
				last_gsc_bin_index,
				gsc_bin_count,
				limiting_mag_local_average,
				last_local_bin_index,
				last_local_bin_index2
			);

			/* Now store everything for the plate */
			SaveLimitingMagnitudes(
				pUpdateCommon,
				pPlate,
				limiting_mag_table_unsorted,
				&missingBinCount,
				&errorPrintFlag
			);

			/* Can flush the plate memory structures */
			free(pPlate->gscLimitingBinBuffer);
			pPlate->gscLimitingBinBuffer = NULL;
			pPlate->gscLimitingBinCount = 0;
			pPlate->gscLimitingBinAlloc = 0;
			pUpdateCommon->currentMemory -= pPlate->gscLimitingBinAlloc * sizeof(GSCLIMITINGBINREC);
		} /* Bottom of the plate loop */

		printf("Max memory %lld, current memory %lld\n", pUpdateCommon->maxMemory, pUpdateCommon->currentMemory);

		/* Now we can read everything back in */

		pUpdateCommon->file_alloc_count = MAXIMUM_BUFFER_MEMORY / sizeof(FILELIMITINGREC);
		if (pUpdateCommon->file_alloc_count > pUpdateCommon->fileWrittenCount) {
			pUpdateCommon->file_alloc_count = pUpdateCommon->fileWrittenCount;
		}

		time(&pUpdateCommon->curTime);
		pUpdateCommon->curTime -= pUpdateCommon->startTime;
		printf(
			"Allocating memory for %zu records total for %zu records allocated, %zu bytes at %ld seconds\n",
			pUpdateCommon->fileWrittenCount,
			pUpdateCommon->file_alloc_count,
			2 * (pUpdateCommon->file_alloc_count * sizeof(FILELIMITINGREC)),
			pUpdateCommon->curTime
		);

		pUpdateCommon->pFileBuffer1 = (PFILELIMITINGREC) calloc(pUpdateCommon->file_alloc_count, sizeof(FILELIMITINGREC));
		if (pUpdateCommon->pFileBuffer1 == NULL) {
			printf("ERROR: failed to allocate pFileBuffer1\n");
			exit(-1);
		}

		pUpdateCommon->pFileBuffer2 = (PFILELIMITINGREC) calloc(pUpdateCommon->file_alloc_count, sizeof(FILELIMITINGREC));
		if (pUpdateCommon->pFileBuffer2 == NULL) {
			printf("ERROR: failed to allocate pFileBuffer2\n");
			exit(-1);
		}

		pUpdateCommon->currentMemory += 2 * (pUpdateCommon->file_alloc_count * sizeof(FILELIMITINGREC));
		if (pUpdateCommon->currentMemory > pUpdateCommon->maxMemory) {
			pUpdateCommon->maxMemory = pUpdateCommon->currentMemory;
		}
		printf("Max memory %lld, current memory %lld\n", pUpdateCommon->maxMemory, pUpdateCommon->currentMemory);

		Close(pUpdateCommon->temporaryHandle);

		pUpdateCommon->temporaryHandle = Open(temporaryname,"r");
		if (pUpdateCommon->temporaryHandle == NULL) {
			printf("Could not open file %s\n",temporaryname);
			exit(-1);
		}

		while (writtenEntryCount < pUpdateCommon->fileWrittenCount && maxGscBin < (pGscBin->total_gsc_bins -1)) {
			maxGscBin = minGscBin;
			currentEntryCount = 0;

			for (gsc_bin_index = minGscBin; gsc_bin_index < pGscBin->total_gsc_bins; gsc_bin_index++) {
				if (currentEntryCount + pUpdateCommon->file_written_table[gsc_bin_index] > pUpdateCommon->file_alloc_count) {
					break;
				}

				currentEntryCount += pUpdateCommon->file_written_table[gsc_bin_index];
			}

			if (gsc_bin_index == minGscBin) {
				printf("ERROR: could not write %d records for gsc_bin_index %d\n", pUpdateCommon->file_written_table[gsc_bin_index], gsc_bin_index);
				exit(-1);
			}

			maxGscBin = gsc_bin_index-1;
			FlushBinCache(pUpdateCommon, pGalaxyCommon, minGscBin, maxGscBin, currentEntryCount);
			minGscBin = maxGscBin+1;
			writtenEntryCount += currentEntryCount;
		}

		if (writtenEntryCount != pUpdateCommon->fileWrittenCount || writtenEntryCount != pUpdateCommon->fileHandledCount || maxGscBin != (pGscBin->total_gsc_bins-1)) {
			printf(
				"ERROR: inconsistent results writtenEntryCount %zu, fileWrittenCount %zu, fileHandledCount %ld, maxGscBin %d, total_gsc_bins %d\n",
				writtenEntryCount,
				pUpdateCommon->fileWrittenCount,
				pUpdateCommon->fileHandledCount,
				maxGscBin,
				pGscBin->total_gsc_bins
			);
		}
 	} else {
		/* Here for nearby galaxies only */
		FlushBinCache(pUpdateCommon, pGalaxyCommon, 0, pGscBin->total_gsc_bins - 1, 0);
	}

	mysql_close(pUpdateCommon->pPhotConnection);
	mysql_close(pUpdateCommon->pConnection);

	if (pUpdateCommon->logHandle != NULL) {
		fclose(pUpdateCommon->logHandle);
	}

	if (pUpdateCommon->platesFlag != 0) {
		if (pUpdateCommon->file_written_table != NULL) {
			free(pUpdateCommon->file_written_table);
		}
		if (pUpdateCommon->file_handled_table != NULL) {
			free(pUpdateCommon->file_handled_table);
		}
	}

	if (pUpdateCommon->checkedTable != NULL) {
		free(pUpdateCommon->checkedTable);
	}
	if (pUpdateCommon->releaseMappingTable != NULL) {
		free(pUpdateCommon->releaseMappingTable);
	}
	if (pUpdateCommon->binWorkQueue != NULL) {
		free(pUpdateCommon->binWorkQueue);
	}

	for (plateIndex = 0; plateIndex < pUpdateCommon->numPlates; plateIndex++) {
		pPlate = &pUpdateCommon->pPlateTable[plateIndex];

		if (pPlate->gscLimitingBinBuffer != NULL) {
			free(pPlate->gscLimitingBinBuffer);
			pPlate->gscLimitingBinBuffer = NULL;
			pPlate->gscLimitingBinCount = 0;
			pPlate->gscLimitingBinAlloc = 0;
		}
	}

	if (pUpdateCommon->pPlateTable != NULL) {
		free(pUpdateCommon->pPlateTable);
	}

	FreeGalaxyCommon(pGalaxyCommon, 1);
	Close(pUpdateCommon->outputHandle);
	Close(pUpdateCommon->indexHandle);
	Close(pUpdateCommon->temporaryHandle);

	/* If there's an error, ignore it. */
	unlink(temporaryname);

	time(&pUpdateCommon->curTime);
	pUpdateCommon->curTime -= pUpdateCommon->startTime;

	printf(
		"Maximum galaxies/variables in a bin %d Average %f; Maximum plates in a bin %d Average %f; Total bins %d\n",
		pUpdateCommon->maxGalaxyBinCount,
		(1.0 * pUpdateCommon->aveGalaxyBinCount) / (1.0 * pGscBin->total_gsc_bins),
		pUpdateCommon->maxPlateBinCount,
		(1.0 * pUpdateCommon->avePlateBinCount) / (1.0 * pGscBin->total_gsc_bins),
		pGscBin->total_gsc_bins
	);
	printf(
		"update_limiting got %d plates and %d solutions %d queries %d good solutions in Execution Time %ld seconds\n",
		plateCount,
		solutionCount,
		GetQueryCount(),
		goodPlateCount,
		pUpdateCommon->curTime
	);
	printf("Max memory %lld, current memory %lld\n", pUpdateCommon->maxMemory, pUpdateCommon->currentMemory);
	printf(
		"zeroBinFlushCount %d staleBinFlushCount %d staleBinMaximum %d missingBinCount %d\n",
		pUpdateCommon->zeroBinFlushCount,
		pUpdateCommon->staleBinFlushCount,
		pUpdateCommon->staleBinMaximum,
		missingBinCount
	);

	for (releaseField = 0; releaseField < RELEASE_FIELD_MAX; releaseField++) {
		printf("%6d plates for field (multiple counting) %s\n", pUpdateCommon->countRelease[releaseField], releaseFieldText[releaseField]);
	}
	printf("\n");

	for (releaseField = 0; releaseField <= RELEASE_FIELD_MAX; releaseField++) {
		printf("%6d plates for field (uniquely assigned) %s\n", releaseFieldPlates[releaseField], releaseFieldText[releaseField]);
	}

	return EXIT_SUCCESS;
}


/*
 * Dec 21, 2012 Edward J. Los - Original Version, adapted from update_photometry.c
 * May  6, 2013 Edward J. Los - Split the limiting magnitudes file into a "nearby" and "plates" file.
 * May 16, 2013 Edward J. Los - Make counts of plate centers as a function of relase field
 * May 19, 2013 Edward J. Los - Correct bin search algorithm
 *                            - Shrink AC scales by 2% to allow for polynomial fitting
 * Jun 18, 2013 Edward J. Los - Use slice algorithm to reduce memory
 *                            - Account for the cropping zones around the plate
 * Jun 21, 2013 Edward J. Los - Rewrite to use files instead of memory.
 * Oct  4, 2013 Edward J. Los - Count only solution 0 to get the plate totals.
 * Mar 22, 2016 Edward J. Los - Add SINGLE_SERIES_DEBUG
 * Jul  7, 2017 Edward J. Los - Correct false error messages for ir02344.
 * Aug  3, 2017 Edward J. Los - Add -d for debugging mode to reduce error printout and -e to limit plate selection
 * Oct 16, 2018 Edward J. Los - Correct large file size overflow errors
 * Jul  5, 2018 Edward J. Los - Add debug code for PhotLoadMagdepCorrections
 */
