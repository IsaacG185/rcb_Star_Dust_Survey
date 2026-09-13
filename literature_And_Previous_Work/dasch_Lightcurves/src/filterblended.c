// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#include <math.h>
#include <errno.h>
#include <time.h>

#include <table.h> // Starbase

#include <libwcs/fitsfile.h>  // jd2ep

#include "pipelineutils.h"
#include "daschunistd.h"

#define CHECK_TOTAL_COUNT 2000000 /* Number of stars in a blend where checking of time begins */
#define MAX_TOTAL_COUNT  99000000 /* Maximum number of stars in a blend */
#define MAX_INPUTNAME 256
#define MAX_BUFFER 1000
#define BLEND_MAGNITUDE_DIFFERENCE 5.0 /* The magnitude difference between the brightest star in a blend and
                                          a dimmer blend star which changes the overall brightness by 0.01 mag */

/* We will set the Sextractor NEW_BLEND_FLAG and change Stdmag if the blended
 *  magnitude difference is greater then MIN_BLEND_MAGNITUDE (THIS CODE BACKED OUT!)
 */

#define DUPLICATE_REF 1  /* Duplicate GSC2.2 ref found */
#define DUPLICATE_NUMBER 2 /* Duplicate Sextracter reference found */

#define CASE_3     0
#define CASE_A     1
#define CASE_B     2
#define CASE_C     3
#define CASE_DA    4
#define CASE_DB    5
#define CASE_DC    6
#define CASE_TOTAL 7
#define CASE_MAX   8

extern GSCBIN gscBin64;

static int rejectedCount[CASE_MAX] = {0,0,0,0,0,0,0,0};
static double rejectedMagnitude[CASE_MAX] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static int acceptedCount[CASE_MAX] = {0,0,0,0,0,0,0,0};
static double acceptedMagnitude[CASE_MAX] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static int blendedAccept = 0;
static int blendedReject = 0;
static PGSCBIN pGscBin = &gscBin64;

static char *caseString[CASE_MAX] = {
  "Step 3 ",
  "Case a ",
  "Case b ",
  "Case c ",
  "Case da",
  "Case db",
  "Case dc",
  "Total  "
};


/* Full contents of the match record */
typedef struct _matchfull {
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  double MAG_ISO;
  double X_IMAGE;
  double Y_IMAGE;
  double ra_1;
  double dec_1;
  double FWHM_IMAGE;
  double FWHM_WORLD;
  double ra_2;
  double dec_2;
  double Stdmag;
  double color;
  double dra;
  double ddec;
  double drad;
  int NUMBER;
  int AFLAGS;
  int BFLAGS;
  int gsc_bin_index;
} MATCHFULL, *PMATCHFULL;

/* Short record for initial read */
typedef struct _matchshort {
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  double MAG_ISO;
  double X_IMAGE;
  double Y_IMAGE;
  double Stdmag;
  double FWHM_WORLD;
  double drad;
  double THRESHOLD;
  double FLUX_MAX;
  double dra;
  double ddec;
  double aLength;      /* Long axis in pixels */
  double bLength;      /* Short axis in pixels */
  double THETA_J2000;  /* Image orientation in degrees */
  double dec_1;
  float RaSigmaPM;     /* mas/yr */
  float DecSigmaPM;    /* mas/yr */
  int NUMBER;
  int VFlag;
  int MAGFlag;
  int AFLAGS;
  int BFLAGS;
  int class;
  int spatial_bin;
  off_t seekPosition;
} MATCHSHORT, *PMATCHSHORT;

typedef struct _check_repeat {
  long long REFNumber;                /* GSC2ID */
  int NUMBER;
  int refCount1;     /* Times REF actually written out */
  int numberCount1;  /* Times NUMBER actually written out */
} CHECK_REPEAT, *PCHECK_REPEAT;

typedef struct _match_entry {
  int linenumber;
  struct _match_entry *REFflink;  /* forward Link for duplicate GSC2 REF's */
  struct _match_entry *REFblink;  /* backward Link for duplicate GSC2 REF's */
  struct _match_entry *NUMBERflink; /* forward Link for duplicate Sextractor numbers */
  struct _match_entry *NUMBERblink; /* backward Link for duplicate Sextractor numbers */
  struct _match_entry *chainFlink;  /* Link to the next chain entry */
  long long REFNumber;                /* GSC2ID */
  double Stdmag;                /* blue magnitude */
  double MAG_ISO;               /* isophotal magnitude */
  double X_IMAGE;
  double Y_IMAGE;
  double Blendedmag;            /* blended magnitude (Stdmag if not blended) */
  double FWHM_WORLD;            /* full width half maximum in arcsec */
  int irec;
  double dra;
  double ddec;
  double aLength;      /* Long axis in pixels */
  double bLength;      /* Short axis in pixels */
  double THETA_J2000;  /* Image orientation in degrees */
  float RaSigmaPM;     /* mas/yr */
  float DecSigmaPM;    /* mas/yr */
  double drad;                  /* separation in arcsec */
  double FLUX_MAX;
  int NUMBER;                   /* Sextractor reference */
  int AFLAGS;                   /* Sextractor flags */
  int BFLAGS;                   /* Sextractor flags */
  int VFlag;                    /* GSC2 variable flag */
  int MAGFlag;                  /* GSC2 magnitude substition */
  int class;                    /* GSC2 class */
  int Afilter_flag;             /* See FILTER_AFLAG_* in pipelineutils.h  */
  int duplicateFlag;            /* See DUPLICATE_* above */
  off_t seekPosition;            /* location of this record in the input file */
  int spatial_bin;
  int chainID;
  int queueIndex;               /* Index into a queue for creating chains */
  int thisIndex;                /* Index into the table for this entry */
  struct _match_entry *flink;   /* Flink to next member of this local calibration bin */

  char selected1;               /* Used during chain searches.  If true, the
                                 * object is already part of the chain */
  char selected2;               /* Used for case d stars to partition the chain into groups
                                   for recursive processing. */
  char selected3;               /* If non zero, the star has been selected for output */
  char selected4;               /* If non zero, the star is in a blend and at least 5 mags dimmer than the brightest star
                                   in the blend */
  char rejected1;               /* If non zero, the star has been marked for a purge */
} MATCH_ENTRY, *PMATCH_ENTRY;

/* drad statistics table */
typedef struct _drad {
  int drad_bin_count;      /* Count in this bin (after clipping) */
  int drad_bin_size;       /* Size of this bin (expanded for insufficient points) */
  int drad_reject_count;   /* Stars rejected because of high drad*/
  double draMedian;
  double draRMS;
  double ddecMedian;
  double ddecRMS;
  int vectorAlloc;
  int vectorCount;
  double *draVector;
  double *ddecVector;
} DRAD, *PDRAD;

/* Index used for sorting */
typedef struct _ref_index {
  long long REFNumber;                /* GSC2ID */
  int index;
} REF_INDEX, *PREF_INDEX;

typedef struct _number_index {
  int NUMBER;                   /* Sextractor reference */
  int index;
} NUMBER_INDEX, *PNUMBER_INDEX;

/* Index entry used for sorting */
typedef struct _chain_entry {
  struct _match_entry *chainFlink;  /* forward Link to chain members */
  int chainID;                  /* ID of this chain */
  int inBand;                   /* If one, this chain is entirely within or below the current band */
  int duplicateFlag;            /* See DUPLICATE_* above */
  int totalCount;               /* Total members in the chain */
  int REFCount;                 /* Duplicate GSC2ID members */
  int NUMBERCount;              /* Duplicate Sextractor members */
} CHAIN_ENTRY, *PCHAIN_ENTRY;

/* This is a table of limiting magnitude estimates used for pre-filtering */
typedef struct _estimate {
  int spatial_bin;
  int iMAG_ISO;
  double Stdmag;
  double MAG_ISO;
} ESTIMATE, *PESTIMATE;


static void
DumpChain(PCHAIN_ENTRY pCurChain,char * rootname)
{
   PMATCH_ENTRY pCurEntry;

  printf("Type %d REF count: %d, NUMBER count %d, total count %d for chain %6d\n",
         pCurChain->duplicateFlag,
         pCurChain->REFCount,
         pCurChain->NUMBERCount,
         pCurChain->totalCount,
         pCurChain->chainID);
  pCurEntry = pCurChain->chainFlink;

  while (pCurEntry != NULL) {
    printf("     Chain %6d member %13lld Stdmag %6.2f %7d MAG_ISO %6.2f Blended %6.2f flag %3d X %6.0f Y %6.0f\n",
           pCurEntry->chainID,
           pCurEntry->REFNumber,
           pCurEntry->Stdmag,
           pCurEntry->NUMBER,
           pCurEntry->MAG_ISO,
           pCurEntry->Blendedmag,
           pCurEntry->Afilter_flag,
           pCurEntry->X_IMAGE,
           pCurEntry->Y_IMAGE);

    pCurEntry = pCurEntry->chainFlink;
  }
}


static int
ReadMatchShort(
  PMATCHSHORT pMatchShort,
  File match_handle,
  TableHead match_header,
  TblDescriptor match_descriptor,
  TableRow* match_row,
  char *rootname
) {
  pMatchShort->seekPosition = Tell(match_handle);
  *match_row = table_rowget(match_handle,match_header,*match_row,NULL,NULL,0);
  if (*match_row == NULL) {
    return 0;
  }

  if (!table_loadrow(match_handle,match_header,*match_row,match_descriptor,(char *)pMatchShort)) {
    printf("ERROR: Read Match table_loadrow failed for %s\n",rootname);
    return 0;
  }

  return 1;
}


static int
ReadMatchFull(
  PMATCHFULL pMatchFull,
  File match_handle,
  TableHead match_header,
  TblDescriptor match_descriptor,
  TableRow* match_row,
  char *rootname
) {
  *match_row = table_rowget(match_handle,match_header,*match_row,NULL,NULL,0);
  if (*match_row == NULL) {
    return 0;
  }

  if (!table_loadrow(match_handle,match_header,*match_row,match_descriptor,(char *)pMatchFull)) {
    printf("ERROR: Read Match table_loadrow failed for %s\n",rootname);
    return 0;
  }

  return 1;
}


static int
NumberCompare(const void *first, const void *second)
{
  int numberFirst = ((PNUMBER_INDEX)first)->NUMBER;
  int numberSecond = ((PNUMBER_INDEX)second)->NUMBER;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  } else {
    return 0;
  }
}


static int
RefCompare(const void *first, const void *second)
{
  long long numberFirst = ((PREF_INDEX)first)->REFNumber;
  long long numberSecond = ((PREF_INDEX)second)->REFNumber;

  if (numberFirst > numberSecond) {
    return 1;
  } else if (numberFirst < numberSecond) {
    return -1;
  } else {
    return 0;
  }
}


static int
InBandCandidate(
  int nx,
  int ny,
  int iyCur,
  PMATCH_ENTRY pCurEntry,
  PCHAIN_ENTRY pChainList
) {
  PCHAIN_ENTRY pCurChain = NULL;

  if (pCurEntry->duplicateFlag != 0) {
    pCurChain = &pChainList[pCurEntry->chainID-1];

    if (pCurChain->inBand == 0) {
      return 0;
    }
  } else if ((pCurEntry->irec % ny) > iyCur) {
    return 0;
  }

  return 1;
}


/* NOTE: case a is the simple case handled within the main routine. */
/*
 * case b:
 * Multiple GSC2.3 stars matching one Sextractor object: this is the classic
 * "blended" object.  Accept only the brightest star, only if the blended
 * magnitude does not change the brightest magnitude by more than
 * MIN_BLEND_MAGNITUDE
 */
static void SelectCaseB(
  int nx,
  int ny,
  PCHAIN_ENTRY pCurChain,
  FILE *blendHandle,
  int iyCur,
  PCHAIN_ENTRY pChainList,
  char *rootname
) {
  PMATCH_ENTRY pCurEntry = pCurChain->chainFlink;
  PMATCH_ENTRY pSelectedEntry = pCurEntry;
  pCurEntry = pCurEntry->chainFlink;
  double tmpFlux;
  double blendedFlux = 0.0;

  /* Find the brightest object */

  while (pCurEntry != NULL) {
    if ((pCurEntry->Stdmag < 99.0) &&
        (pCurEntry->Stdmag < pSelectedEntry->Stdmag)) {
      pSelectedEntry = pCurEntry;
    } else if ((pCurEntry->Stdmag == pSelectedEntry->Stdmag) &&
               (pCurEntry->REFNumber < pSelectedEntry->REFNumber)) {
      pSelectedEntry = pCurEntry;

    }
    pCurEntry = pCurEntry->chainFlink;
  }
  pSelectedEntry->selected3 = 1;
  acceptedCount[CASE_B]++;
  acceptedMagnitude[CASE_B] += pSelectedEntry->Stdmag;

  /* Get the blended magnitude and compute statistics */

  pCurEntry = pCurChain->chainFlink;
  while (pCurEntry != NULL) {
    if (pCurEntry->Stdmag < 90.0) {
      tmpFlux = exp10(-pCurEntry->Stdmag/2.5);
      blendedFlux += tmpFlux;
    }

    if (pCurEntry->selected3 == 0) {
      rejectedCount[CASE_B]++;
      rejectedMagnitude[CASE_B] += pCurEntry->Stdmag;
      if ((pCurEntry->Stdmag-pSelectedEntry->Stdmag) > BLEND_MAGNITUDE_DIFFERENCE) {
        /* Too dim. Toss this one prematurely */
        pCurEntry->selected4 = 1;
      }
    }
    pCurEntry = pCurEntry->chainFlink;
  }

  if (blendedFlux == 0) {
    pSelectedEntry->Blendedmag = 99.0;
  } else {
    pSelectedEntry->Blendedmag = -2.5*log10(blendedFlux);
  }

  if ((pSelectedEntry->Stdmag > 90.0) ||
      ((pSelectedEntry->Stdmag - pSelectedEntry->Blendedmag) > MIN_BLEND_MAGNITUDE)) {
    pSelectedEntry->Afilter_flag |= (1 << FILTER_AFLAG_CASEB);
    blendedReject++;

    if (InBandCandidate(nx,ny,iyCur,pSelectedEntry,pChainList)) {
      fprintf(blendHandle,"%d\t%.2f\t%.2f\n",pSelectedEntry->NUMBER,pSelectedEntry->Stdmag,pSelectedEntry->Blendedmag);
    } else {
      blendedAccept++;
    }
  }
}


/* case c:
 * One GSC2.3 star matching multiple Sextractor objects: Assign this star
 * to the brightest of the Sextractor objects.
 */
static void
SelectCaseC(PCHAIN_ENTRY pCurChain, char *rootname)
{
  PMATCH_ENTRY pCurEntry = pCurChain->chainFlink;
  PMATCH_ENTRY pSelectedEntry = pCurEntry;
  pCurEntry = pCurEntry->chainFlink;

  /* Find the brightest object */
  while (pCurEntry != NULL) {
    if (pCurEntry->MAG_ISO < pSelectedEntry->MAG_ISO) {
      pSelectedEntry = pCurEntry;
    } else if ((pCurEntry->MAG_ISO == pSelectedEntry->MAG_ISO) &&
               (pCurEntry->REFNumber < pSelectedEntry->REFNumber)) {
      pSelectedEntry = pCurEntry;
    }

    pCurEntry = pCurEntry->chainFlink;
  }

  pSelectedEntry->selected3 = 1;
  pSelectedEntry->Afilter_flag |= (1 << FILTER_AFLAG_CASEC);
  acceptedCount[CASE_C]++;
  acceptedMagnitude[CASE_C] += pSelectedEntry->Stdmag;

  /* Compute statistics */
  pCurEntry = pCurChain->chainFlink;
  while (pCurEntry != NULL) {
    if (pCurEntry->selected3 == 0) {
      rejectedCount[CASE_C]++;
      rejectedMagnitude[CASE_C] += pCurEntry->Stdmag;
    }

    pCurEntry = pCurEntry->chainFlink;
  }
}


/*
 * case d:
 * Multiple GSC2.3 stars matching multiple Sextractor objects.  This is a
 * problematic case which for now will be dealt with arbitrarily:  Start with the
 * brightest Sextractor object.  If it has multiple GSC2.3 objects within its
 * (factor*FWHM_WORLD) then treat these objects as case b (case db below).
 * Otherwise, treat it as case c (case dc below).  Repeat this process for all the
 * remaining stars in the group.  Any leftover individual pairs are treated as
 * case a (case da below).  Another option for this case might be to reject the
 * entries altogether.
 */
static void
SelectCaseD(
  int nx,
  int ny,
  PCHAIN_ENTRY pCurChain,
  FILE *blendHandle,
  int iyCur,
  PCHAIN_ENTRY pChainList,
  char *rootname
) {
  PMATCH_ENTRY pCurEntry = pCurChain->chainFlink;
  PMATCH_ENTRY pSelectedEntry = NULL;
  PMATCH_ENTRY pTmpChain;
  int tmpChainCount = 0;
  int tmpCount;
  double tmpMagnitude;
  double blendedFlux = 0.0;
  double tmpFlux;

  /* Find the brightest Sextractor object */
  while (pCurEntry != NULL) {
    if (pCurEntry->selected2 == 0) {
      if (pSelectedEntry == NULL) {
        pSelectedEntry = pCurEntry;
      } else {
        if (pCurEntry->MAG_ISO < pSelectedEntry->MAG_ISO) {
          pSelectedEntry = pCurEntry;
        } else if ((pCurEntry->MAG_ISO == pSelectedEntry->MAG_ISO) &&
                   (pCurEntry->REFNumber < pSelectedEntry->REFNumber)) {
          pSelectedEntry = pCurEntry;
        }
      }
    }

    pCurEntry = pCurEntry->chainFlink;
  }

  if (pSelectedEntry == NULL) {
    /* No more stars */
    return;
  }

  /* Bugfix of April 29, 2008 - for all objects with the
     brightest MAG_ISO, select the one with the brightest Stdmag */
  pCurEntry = pCurChain->chainFlink;
  while (pCurEntry != NULL) {
    if (pCurEntry->selected2 == 0) {
      if ((pCurEntry->MAG_ISO == pSelectedEntry->MAG_ISO) &&
          (pCurEntry->Stdmag < 90.0) &&
          (pCurEntry->Stdmag < pSelectedEntry->Stdmag)) {
        pSelectedEntry = pCurEntry;
      } else if ((pCurEntry->MAG_ISO == pSelectedEntry->MAG_ISO) &&
                 (pCurEntry->Stdmag == pSelectedEntry->Stdmag) &&
                 (pCurEntry->REFNumber < pSelectedEntry->REFNumber)) {
        pSelectedEntry = pCurEntry;

      }
    }
    pCurEntry = pCurEntry->chainFlink;
  }

  /* Select this entry */
  pSelectedEntry->selected3 = 1;
  if ((pSelectedEntry->NUMBERflink != NULL) ||
      (pSelectedEntry->NUMBERblink != NULL)) {
    /* This is case d-b or d-a. Reject other instances of this
       GSC2.3 object */
    /* Go to the start of the REF chain */
    tmpMagnitude = 0.0;
    tmpCount = 0;
    pTmpChain = pSelectedEntry;
    while (pTmpChain->REFblink) {
      pTmpChain = pTmpChain->REFblink;
    }

    pCurEntry = pTmpChain;

    while (pCurEntry) {
      if ((pCurEntry->selected2 == 0) &&
          (pCurEntry->selected3 == 0)) {
        pCurEntry->selected2 = 1;
        tmpCount++;
        tmpMagnitude += pCurEntry->Stdmag;
      }
      pCurEntry = pCurEntry->REFflink;
    }

    /* Go to the start of the NUMBER chain */
    pTmpChain = pSelectedEntry;

    while (pTmpChain->NUMBERblink) {
      pTmpChain = pTmpChain->NUMBERblink;
    }

    /* Now find the brightest GSC2 object in the chain */

    pCurEntry = pTmpChain;

    while (pCurEntry) {
      if (pCurEntry->selected2 == 0) {
        pCurEntry->selected2 = 1;
        tmpChainCount++;

        if (pCurEntry->Stdmag < 90.0) {
          tmpFlux = exp10(-pCurEntry->Stdmag/2.5);
          blendedFlux += tmpFlux;
        }

        if (pCurEntry->selected3 == 0) {
          tmpCount++;
          tmpMagnitude += pCurEntry->Stdmag;
          if ((pCurEntry->Stdmag-pSelectedEntry->Stdmag) > BLEND_MAGNITUDE_DIFFERENCE) {
            /* Too dim. Toss this one prematurely */
            pCurEntry->selected4 = 1;
          }
        }
      }

      pCurEntry = pCurEntry->NUMBERflink;
    }

    if (blendedFlux == 0) {
      pSelectedEntry->Blendedmag = 99.0;
    } else {
      pSelectedEntry->Blendedmag = -2.5*log10(blendedFlux);
    }

    if ((pSelectedEntry->Stdmag > 90.0) ||
        ((pSelectedEntry->Stdmag - pSelectedEntry->Blendedmag) > MIN_BLEND_MAGNITUDE)) {
      pSelectedEntry->Afilter_flag |= (1 << FILTER_AFLAG_CASED);
      blendedReject++;

      if (InBandCandidate(nx,ny,iyCur,pSelectedEntry,pChainList)) {
        fprintf(blendHandle,"%d\t%.2f\t%.2f\n",pSelectedEntry->NUMBER,pSelectedEntry->Stdmag,pSelectedEntry->Blendedmag);
      } else {
        blendedAccept++;
      }
    }

    if (tmpChainCount == 1) {
      /* This is case da */
      acceptedCount[CASE_DA]++;
      acceptedMagnitude[CASE_DA] += pSelectedEntry->Stdmag;
      rejectedCount[CASE_DA] += tmpCount;
      rejectedMagnitude[CASE_DA] += tmpMagnitude;
    } else {
      /* This is case db */
      acceptedCount[CASE_DB]++;
      acceptedMagnitude[CASE_DB] += pSelectedEntry->Stdmag;
      rejectedCount[CASE_DB] += tmpCount;
      rejectedMagnitude[CASE_DB] += tmpMagnitude;

    }
  } else {
    tmpCount = 0;
    tmpMagnitude = 0.0;
    tmpChainCount = 0;

    pSelectedEntry->Afilter_flag |= (1 << FILTER_AFLAG_CASED);
    /* This is case d-c or d-a: go to the start of the NUMBER
       chain and reject all other instances */
    pTmpChain = pSelectedEntry;

    while (pTmpChain->NUMBERblink) {
      pTmpChain = pTmpChain->NUMBERblink;
    }

    pCurEntry = pTmpChain;

    while (pCurEntry) {
      if ((pCurEntry->selected2 == 0) &&
          (pCurEntry->selected3 == 0)) {
        pCurEntry->selected2 = 1;
        tmpCount++;
        tmpMagnitude += pCurEntry->Stdmag;
      }
      pCurEntry = pCurEntry->NUMBERflink;
    }

    /* go to the start of the chain */
    pTmpChain = pSelectedEntry;

    while (pTmpChain->REFblink) {
      pTmpChain = pTmpChain->REFblink;
    }

    /* Now find the brightest GSC2 object in the chain */

    pCurEntry = pTmpChain;

    while (pCurEntry) {
      if (pCurEntry->selected2 == 0) {
        pCurEntry->selected2 = 1;
        tmpChainCount++;
        if (pCurEntry->selected3 == 0) {
          tmpCount++;
          tmpMagnitude += pCurEntry->Stdmag;
        }
      }

      pCurEntry = pCurEntry->REFflink;
    }

    if (tmpChainCount == 1) {
      /* This is case da */
      acceptedCount[CASE_DA]++;
      acceptedMagnitude[CASE_DA] += pSelectedEntry->Stdmag;
      rejectedCount[CASE_DA] += tmpCount;
      rejectedMagnitude[CASE_DA] += tmpMagnitude;
    } else {
      /* This is case dc */
      acceptedCount[CASE_DC]++;
      acceptedMagnitude[CASE_DC] += pSelectedEntry->Stdmag;
      rejectedCount[CASE_DC] += tmpCount;
      rejectedMagnitude[CASE_DC] += tmpMagnitude;
    }
  }

  /* Do it again for the remaining stars in the group */
  SelectCaseD(nx,ny,pCurChain,blendHandle,iyCur,pChainList,rootname);
}


/* Place the entry on a fifo queue */
static void
AddToQueue(PMATCH_ENTRY pCurEntry, PMATCH_ENTRY pEntryList, int *pQueueSize)
{
  int queueSize = *pQueueSize;
  PMATCH_ENTRY pQueueEntry;

  if (pCurEntry->selected1) {
    return;
  }

  pQueueEntry = &pEntryList[queueSize];
  pQueueEntry->queueIndex = pCurEntry->thisIndex;
  queueSize++;
  *pQueueSize = queueSize;
}


static void
CountRefs(
  PCHAIN_ENTRY pCurChain,
  PMATCH_ENTRY pEntryList,
  int *pQueueSize,
  int *pMaxTotalCount,
  int solutionNumber,
  char *rootname,
  char * qualifier,
  int maxTimeout,
  int skipPVUpdateFlag,
  time_t startTime
) {
  int queueSize = *pQueueSize;
  PMATCH_ENTRY pQueueEntry;
  int index;
  PMATCH_ENTRY pCurEntry;
  int readFlag = 0;
  int FitWCS;
  static int checkCounter = 0;
  time_t curTime;

  queueSize--;

  if (queueSize < 0) {
    printf("ERROR: queueSize is %d in CountRefs\n",queueSize);
    exit(1);
  }

  *pQueueSize = queueSize;
  pQueueEntry = &pEntryList[queueSize];
  index = pQueueEntry->queueIndex;
  pCurEntry = &pEntryList[index];

  if (pCurEntry->selected1) {
    return;
  }

  pCurEntry->selected1 = 1;
  pCurEntry->chainID = pCurChain->chainID;
  pCurEntry->chainFlink = pCurChain->chainFlink;
  pCurChain->chainFlink = pCurEntry;
  pCurChain->duplicateFlag |= pCurEntry->duplicateFlag;
  pCurChain->totalCount += 1;

  if ((maxTimeout > 0) &&
      (*pMaxTotalCount > CHECK_TOTAL_COUNT)) {
    if (++checkCounter >= CHECK_TOTAL_COUNT) {
      checkCounter = 0;
      time(&curTime);
      if ((curTime - startTime) > maxTimeout) {
        if ((qualifier[0] == 0) && (skipPVUpdateFlag == 0)) {
          SetMosaicFitWCS2("FilterblendedTimeout",rootname,solutionNumber,1,readFlag,&FitWCS);
        }
        printf("ERRORD: filterblended maxTimeout of %ld/%d seconds maxTotalCount %d reached for %s\n",
                curTime-startTime,
                maxTimeout,
                *pMaxTotalCount,
                rootname);
        exit(1);
      }
    }
  }

  if (pCurChain->totalCount > *pMaxTotalCount) {
    *pMaxTotalCount = pCurChain->totalCount;

    if (*pMaxTotalCount > MAX_TOTAL_COUNT) {
      if ((qualifier[0] == 0)  && (skipPVUpdateFlag == 0)) {
        SetMosaicFitWCS2("FilterblendedTimeout",rootname,solutionNumber,1,readFlag,&FitWCS);
      }
      printf("ERROR: maxTotalCount %d exceeds %d for %s\n",*pMaxTotalCount,MAX_TOTAL_COUNT,rootname);
      exit(1);
    }
  }

  if ((pCurEntry->REFflink) || (pCurEntry->REFblink)) {
    pCurChain->REFCount += 1;
    if (pCurEntry->REFflink) {
      AddToQueue(pCurEntry->REFflink,pEntryList,pQueueSize);
    }

    if (pCurEntry->REFblink) {
      AddToQueue(pCurEntry->REFblink,pEntryList,pQueueSize);
    }
  }

  if ((pCurEntry->NUMBERflink) || (pCurEntry->NUMBERblink)) {
    pCurChain->NUMBERCount += 1;

    if (pCurEntry->NUMBERflink) {
      AddToQueue(pCurEntry->NUMBERflink,pEntryList,pQueueSize);
    }

    if (pCurEntry->NUMBERblink) {
      AddToQueue(pCurEntry->NUMBERblink,pEntryList,pQueueSize);
    }
  }
}


static void
AddSelectedEntry(
  int nx,
  int ny,
  int iyCur,
  PDRAD dradTable,
  PMATCH_ENTRY pSelectedEntry,
  double THRESHOLD,
  signed long long *pCurMemory,
  signed long long *pMaxMemory
) {
  PDRAD pDrad;
  double *tmpdraVector;
  double *tmpddecVector;
  int iy = pSelectedEntry->irec % ny;

  if (((iyCur == 0) && (iy <= 1)) || (iy == (iyCur+1))) {
    /* Skip observations that have issues or that have a high proper motion error */

    if  ((pSelectedEntry->AFLAGS < FILTER_AFLAG_BAD) &&
         ((pSelectedEntry->BFLAGS & (1 << FILTER_BFLAG_PMERROR)) == 0)) {

      pDrad = &dradTable[pSelectedEntry->irec];
      if (pDrad->vectorCount >= pDrad->vectorAlloc) {
        pDrad->vectorAlloc += 100;
        tmpdraVector = realloc(pDrad->draVector,pDrad->vectorAlloc*sizeof(double));
        CalcMemory(1,100 * sizeof(double),pCurMemory,pMaxMemory);

        if (tmpdraVector == NULL) {
          printf("ERROR allocating tmpdravector of size %d\n",pDrad->vectorAlloc);
          exit(1);
        }

        pDrad->draVector = tmpdraVector;
        tmpdraVector = NULL;
        tmpddecVector = realloc(pDrad->ddecVector,pDrad->vectorAlloc*sizeof(double));
        CalcMemory(1,100 * sizeof(double),pCurMemory,pMaxMemory);

        if (tmpddecVector == NULL) {
          printf("ERROR allocating tmpddecvector of size %d\n",pDrad->vectorAlloc);
          exit(1);
        }

        pDrad->ddecVector = tmpddecVector;
        tmpddecVector = NULL;
      }


      pDrad->draVector[pDrad->vectorCount] = pSelectedEntry->dra;
      pDrad->ddecVector[pDrad->vectorCount] = pSelectedEntry->ddec;
      pDrad->vectorCount++;
    }
  }
}


/* Go through the chains and select the one with the brightest Stdmag and the
 * brightest MAG_ISO.  Put this one in the dradTable
 */
static void
AddBrightest(
  int nx,
  int ny,
  int iyCur,
  PDRAD dradTable,
  PMATCH_ENTRY pEntryList,
  int goodlines,
  PCHAIN_ENTRY pChainList,
  int maxChainID,
  double THRESHOLD,
  signed long long *pCurMemory,
  signed long long *pMaxMemory,
  char *rootname
) {
  PCHAIN_ENTRY pCurChain;
  PMATCH_ENTRY pCurEntry;
  PMATCH_ENTRY pSelectedEntry;
  int chainIndex;

  for (chainIndex = 0; chainIndex < maxChainID; chainIndex++) {
    pCurChain = &pChainList[chainIndex];
    pCurEntry = pCurChain->chainFlink;
    pSelectedEntry = pCurEntry;

    while (pCurEntry != NULL) {
      if (pCurEntry->Stdmag < pSelectedEntry->Stdmag) {
        pSelectedEntry = pCurEntry;
      }

      pCurEntry = pCurEntry->chainFlink;
    }

    /* Now for the brightest GSC object, select the brightest image */
    pCurEntry = pCurChain->chainFlink;

    while (pCurEntry != NULL) {
      if ((pCurEntry->MAG_ISO <= pSelectedEntry->MAG_ISO) &&
          (pCurEntry->Stdmag == pSelectedEntry->Stdmag)) {
        pSelectedEntry = pCurEntry;
      }

      pCurEntry = pCurEntry->chainFlink;
    }

    if ((pSelectedEntry->Stdmag < 90.0) && (pSelectedEntry->MAG_ISO < 90.0)) {
      AddSelectedEntry(nx,ny,iyCur,dradTable,pSelectedEntry,THRESHOLD,pCurMemory,pMaxMemory);
    }
  }
}


/*
 * Compute Drad
 */
static void
ComputeDrad(
  int nx,
  int ny,
  int iyCur,
  PDRAD dradTable,
  double THRESHOLD,
  signed long long *pCurMemory,
  signed long long *pMaxMemory
) {
  PDRAD pDrad;
  PDRAD pDrad2;
  double *draVector = NULL;
  double *tmpdraVector;
  double *ddecVector = NULL;
  double *tmpddecVector;
  int vectorAlloc = 0;
  int vectorCount;
  int ix;
  int iy;
  int irec;
  int draCount;
  int ddecCount;
  int dx;
  int dy;
  int iix;
  int iiy;
  int index;
  iy = iyCur;
  int vectorAllocOld;

  for (ix = 0; ix < nx; ix++) {
    vectorCount = 0;
    irec = iy + (ix * ny);
    pDrad = &dradTable[irec];

    /* We must copy the vectors since clipping will change their contents */
    if ((pDrad->vectorCount) >= vectorAlloc) {
      vectorAllocOld = vectorAlloc;
      vectorAlloc = pDrad->vectorCount+100;
      tmpdraVector = realloc(draVector,vectorAlloc*sizeof(double));
      CalcMemory(1,(vectorAlloc-vectorAllocOld) * sizeof(double),pCurMemory,pMaxMemory);

      if (tmpdraVector == NULL) {
        printf("ERROR allocating tmpdravector of size %d\n",vectorAlloc);
        exit(1);
      }

      draVector = tmpdraVector;
      tmpdraVector = NULL;
      tmpddecVector = realloc(ddecVector,vectorAlloc*sizeof(double));
      CalcMemory(1,(vectorAlloc-vectorAllocOld) * sizeof(double),pCurMemory,pMaxMemory);

      if (tmpddecVector == NULL) {
        printf("ERROR allocating tmpddecvector of size %d\n",vectorAlloc);
        exit(1);
      }

      ddecVector = tmpddecVector;
      tmpddecVector = NULL;
    }

    for (index = 0; index < pDrad->vectorCount; index++) {
      draVector[index] = pDrad->draVector[index];
      ddecVector[index] = pDrad->ddecVector[index];
      vectorCount++;
    }

    draCount  = CalcMedianAndRMS(pDrad->vectorCount,STARS_PER_DMAGBIN,draVector,&pDrad->draMedian,&pDrad->draRMS,1,3.0,0);
    ddecCount = CalcMedianAndRMS(pDrad->vectorCount,STARS_PER_DMAGBIN,ddecVector,&pDrad->ddecMedian,&pDrad->ddecRMS,1,3.0,0);

    if ((draCount == 0) || (ddecCount == 0)) {
      /* We did not have sufficient stars here to calculate drad.  Expand to all of the adjacent bins */
      vectorCount = 0;

      for (dx = -1; dx <= 1; dx++) {
        for (dy = -1; dy <= 1; dy++) {
          iix = dx+ix;
          if ((iix < 0) || (iix >= nx)) {
            continue;
          }

          iiy = dy+iy;
          if ((iiy < 0) || (iiy >= ny)) {
            continue;
          }

          pDrad->drad_bin_size++;
          irec = iiy + (iix * ny);
          pDrad2 = &dradTable[irec];

          if (pDrad2->vectorCount > 0) {
            if ((vectorCount+pDrad2->vectorCount) >= vectorAlloc) {
              vectorAllocOld = vectorAlloc;
              vectorAlloc = vectorCount+pDrad2->vectorCount+100;
              tmpdraVector = realloc(draVector,vectorAlloc*sizeof(double));
              CalcMemory(1,(vectorAlloc-vectorAllocOld) * sizeof(double),pCurMemory,pMaxMemory);
              if (tmpdraVector == NULL) {
                printf("ERROR allocating tmpdravector of size %d\n",vectorAlloc);
                exit(1);
              }

              draVector = tmpdraVector;
              tmpdraVector = NULL;
              tmpddecVector = realloc(ddecVector,vectorAlloc*sizeof(double));
              CalcMemory(1,(vectorAlloc-vectorAllocOld) * sizeof(double),pCurMemory,pMaxMemory);
              if (tmpddecVector == NULL) {
                printf("ERROR allocating tmpddecvector of size %d\n",vectorAlloc);
                exit(1);
              }

              ddecVector = tmpddecVector;
              tmpddecVector = NULL;
            }

            for (index = 0; index < pDrad2->vectorCount; index++) {
              draVector[vectorCount] = pDrad2->draVector[index];
              ddecVector[vectorCount] = pDrad2->ddecVector[index];
              vectorCount++;
            }
          }
        }
      }

      draCount  = CalcMedianAndRMS(vectorCount,STARS_PER_DMAGBIN,draVector,&pDrad->draMedian,&pDrad->draRMS,1,3.0,0);
      ddecCount = CalcMedianAndRMS(vectorCount,STARS_PER_DMAGBIN,ddecVector,&pDrad->ddecMedian,&pDrad->ddecRMS,1,3.0,0);
    } else {
      pDrad->drad_bin_size = 1;
    }

    if (ddecCount == 0) {
      pDrad->ddecRMS = 99.0;
    }

    if (draCount == 0) {
      pDrad->draRMS = 99.0;
    }

    if (ddecCount < draCount) {
      pDrad->drad_bin_count = ddecCount;
    } else {
      pDrad->drad_bin_count = draCount;
    }
  }

  if (draVector != NULL) {
    free(draVector);
    CalcMemory(-1,vectorAlloc * sizeof(double),pCurMemory,pMaxMemory);
  }

  if (ddecVector != NULL) {
    free(ddecVector);
    CalcMemory(-1,vectorAlloc * sizeof(double),pCurMemory,pMaxMemory);
  }

  /* Clean up vectors that we no longer need */
  iy = iyCur - 1;

  if (iy >= 0) {
    for (ix = 0; ix < nx; ix++) {
      irec = iy + (ix * ny);
      pDrad = &dradTable[irec];

      if (pDrad->vectorAlloc > 0) {
        free(pDrad->draVector);
        free(pDrad->ddecVector);
        pDrad->draVector = NULL;
        pDrad->ddecVector = NULL;
        CalcMemory(-1,2*pDrad->vectorAlloc * sizeof(double),pCurMemory,pMaxMemory);
        pDrad->vectorCount = 0;
        pDrad->vectorAlloc = 0;
      }
    }
  }
}


static void
PurgeHighDrad(
  int nx,
  int ny,
  int iyCur,
  PDRAD dradTable,
  PMATCH_ENTRY pEntryList,
  int *pGoodlines,
  int *pHighDradCount,
  double maxDrad,
  double scale,
  char *rootname,
  FILE *regionHandle
) {
  int goodlines = *pGoodlines;
  int highDradCount = 0;
  PMATCH_ENTRY pCurEntry = pEntryList;
  PMATCH_ENTRY pLastEntry;
  PDRAD pDrad;
  int curIndex = 0;
  double dra;
  double ddec;
  double draderror;
  double aLengthMax;
  double bLengthMax;
  double blendAngle;

  *pHighDradCount = 0;

  while (curIndex < goodlines) {
    pDrad = &dradTable[pCurEntry->irec];

    if ((pCurEntry->irec % ny) > iyCur) {
      curIndex++;
      pCurEntry++;
      continue; /* We have not computed these drads yet */
    }

    if (pDrad->drad_bin_count >= STARS_PER_DMAGBIN) {
      /* We have enough stars to adjust the dra and ddec.  Otherwise, check against the maximum drad error */
      dra = pCurEntry->dra - pDrad->draMedian;
      ddec = pCurEntry->ddec - pDrad->ddecMedian;
      draderror = sqrt(sqr(pDrad->draRMS) + sqr(pDrad->ddecRMS)) * DRAD_REJECT_FACTOR;
      if (draderror > maxDrad) {
        draderror = maxDrad;
      }
    } else {
      dra = pCurEntry->dra;
      ddec = pCurEntry->ddec;
      draderror = maxDrad;
    }

    aLengthMax =  (pCurEntry->aLength*scale);
    if (aLengthMax < draderror) {
      aLengthMax = draderror;
    }

    bLengthMax =  (pCurEntry->bLength*scale);
    if (bLengthMax < draderror) {
      bLengthMax = draderror;
    }

    /* Convert J2000 to IMAGE */
    blendAngle = pCurEntry->THETA_J2000 + 90.;
    if (blendAngle > 90) {
      blendAngle -= 180;
    }

    if (!CheckBlend(-dra, /* x axis is positive */
                    ddec,
                    blendAngle,
                    aLengthMax,
                    bLengthMax,
                    0.0,
                    0.0,
                    0.0,
                    pCurEntry->RaSigmaPM,
                    pCurEntry->DecSigmaPM)) {
      /* This image fails, toss it and replace it with the last entry on the list */
      highDradCount++;
      goodlines--;

      if (curIndex != goodlines) {
        pLastEntry = &pEntryList[goodlines];
        memcpy(pCurEntry,pLastEntry,sizeof(MATCH_ENTRY));
      }
    } else {
      /* Accept this image */
      if (regionHandle != NULL) {
        /* Here we are not as accuracte the mosaic is not well aligned to the North */
        fprintf(regionHandle,"box %f %f %f %f %f # color=green\n",
                pCurEntry->X_IMAGE,
                pCurEntry->Y_IMAGE,
                2*pCurEntry->aLength,
                2*pCurEntry->bLength,
                pCurEntry->THETA_J2000);
        fprintf(regionHandle,"box %f %f %f %f %f # color=red\n",
                pCurEntry->X_IMAGE,
                pCurEntry->Y_IMAGE,
                2*aLengthMax/scale,
                2*bLengthMax/scale,
                pCurEntry->THETA_J2000);
      } /* regionHandle test */

      curIndex++;
      pCurEntry++;

    }
  }

  *pGoodlines = goodlines;
  *pHighDradCount = highDradCount;
}


static void
ReadBand(
  int nx,
  int ny,
  int iyCur,
  int mosaicWidth,
  int mosaicHeight,
  PMATCH_ENTRY *ppEntryList,
  int *pmatch_alloc,
  size_t *pmatch_nrecs,
  int *pgoodlines,
  int *pbandlines,
  int *phighdradcount,
  File match_handle,
  TableHead match_header,
  TblDescriptor match_descriptor_short,
  TableRow *match_row,
  char *rootname,
  PMATCH_ENTRY pFWHM_WORLD_Entry,
  signed long long *pCurMemory,
  signed long long *pMaxMemory,
  double maxDradPixels,
  double plateepoch,
  double scale,
  double *pTHRESHOLD,
  double verbose,
  long long *pMaxREF,
  double *pOldY_IMAGE,
  double *poldiy,
  double *Stdmag_Table,
  int *estimateFilterCount
) {
  PMATCH_ENTRY pEntryList = *ppEntryList;
  int match_alloc = *pmatch_alloc;
  int match_nrecs = *pmatch_nrecs;
  MATCHSHORT match_record_short;
  PMATCHSHORT pMatchShort = &match_record_short;
  int highdradflag;
  int highdradcount = *phighdradcount;
  int goodlines = *pgoodlines;
  int ix;
  int iy = 0;
  int bandlines = 0;
  PMATCH_ENTRY pCurEntry;
  int index;
  int refType;
  double Stdmag;
  int iMAG_ISO;
  double factor;
  double blendAngle;

  /* If we have existing records, see how many of them are in the current band */
  for (index = 0; index < goodlines; index++) {
    pCurEntry = &pEntryList[index];

    iy = (pCurEntry->Y_IMAGE * ny) /(1.0 * mosaicHeight);
    if (iy < 0) {
      iy = 0;
    }

    if (iy >= ny) {
      iy = ny-1;
    }

    if (iy <= iyCur) {
      bandlines++;
    }
  }

  iy = 0;

  while (1) {
    if (iy > (iyCur+2)) {
      /* Exit if we have reached the top of the second highest band */
      break;
    }

    if (goodlines >= match_alloc) {
      match_alloc += 1000;
      ReAlloc(pEntryList,match_alloc * sizeof(MATCH_ENTRY));
      CalcMemory(1,1000*sizeof(MATCH_ENTRY),pCurMemory,pMaxMemory);
      if (pEntryList == NULL) {
        printf("ERROR: Failed to allocate pEntryList size %d %zu %s\n",match_alloc,match_alloc * sizeof(MATCH_ENTRY),rootname);
        exit(1);
      }
    }

    pCurEntry = &pEntryList[goodlines];

    memset(pMatchShort,0,sizeof(MATCHSHORT));
    if(ReadMatchShort(pMatchShort,match_handle,match_header,match_descriptor_short,match_row,rootname) == 0) {
      break;
    }
    match_nrecs++;

    if (pMatchShort->Y_IMAGE < *pOldY_IMAGE) {
      printf("ERROR: Y_IMAGE decreases from %f to %f for record %d\n",pMatchShort->Y_IMAGE,*pOldY_IMAGE,match_nrecs);
      exit(1);
    }

    *pOldY_IMAGE = pMatchShort->Y_IMAGE;

    if ((Stdmag_Table != NULL) &&
        (pMatchShort->spatial_bin != 0) &&
        (pMatchShort->Stdmag < 90.0)) {
      /* Here are are filtering out objects from over-populated MAG_ISO bins */
      iMAG_ISO = (pMatchShort->MAG_ISO - MIN_MAG_ISO)*MAG_ISO_BINS/(MAX_MAG_ISO - MIN_MAG_ISO);
      if (iMAG_ISO < 0) {
        iMAG_ISO = 0;
      } else if (iMAG_ISO >= MAG_ISO_BINS) {
        iMAG_ISO = MAG_ISO_BINS -1;
      }
      Stdmag = Stdmag_Table[iMAG_ISO + (MAG_ISO_BINS * (pMatchShort->spatial_bin-1))];
      if (pMatchShort->Stdmag > Stdmag) {
        (*estimateFilterCount)++;
        continue;
      }
    }

    // Check to see if dra and ddec lie outside the rectangle defining the
    // object.  For this check, expand the object dimensions by MaxDradPixels
    //
    // Note: CheckBlend is written for THETA_IMAGE where the first parameter
    // defines horizontal, the second vertical, and the third the angle measured
    // conterclockwise from horizontal. Since THETA_J2000 is measured from North
    // and the angle is clockwise from North, reverse horizontal and vertical
    // and change the direction of the angle.

    // Convert the proper motion errors to arcseconds and expand the GSC object
    // by these errors (proper motion has already been applied to the matched
    // image
    factor = cos(DEGREES_TO_RAD *  pMatchShort->dec_1);
    if (factor == 0.0) {
      factor = 1.0;
    }

    pMatchShort->RaSigmaPM *= fabs((PROPER_MOTION_SIGMA*plateepoch)/(1000.0*factor)); /* Convert to arcseconds in R.A. */
    pMatchShort->DecSigmaPM *= fabs((PROPER_MOTION_SIGMA*plateepoch)/(1000.0)); /* Convert to arcseconds */

    /* Convert J2000 to IMAGE */
    blendAngle = pMatchShort->THETA_J2000 + 90.;
    if (blendAngle > 90) {
      blendAngle -= 180;
    }

    highdradflag = !CheckBlend(pMatchShort->dra, /* Center of image  with respect to the GSC object */
                               pMatchShort->ddec, /* Center of image  with respect to the GSC object */
                               blendAngle, /* Angle of the axis with respect the X axis */
                               (pMatchShort->aLength + maxDradPixels)*scale,
                               (pMatchShort->bLength + maxDradPixels)*scale,
                               0.0, /* Location of the GSC object */
                               0.0,
                               0.0, /* Angle of the GSC object */
                               pMatchShort->RaSigmaPM, /* The GSC object size is the uncertainty in proper motion */
                               pMatchShort->DecSigmaPM);

    if (highdradflag) {
      highdradcount++;
      rejectedCount[CASE_3]++;
      rejectedMagnitude[CASE_3] += pMatchShort->Stdmag;
      continue;
    } else {
      acceptedCount[CASE_3]++;
      acceptedMagnitude[CASE_3] += pMatchShort->Stdmag;
    }

    memset(pCurEntry,0,sizeof(MATCH_ENTRY));

    /* Decode the GSC2 ID */
    GetREFNumber(pMatchShort->REF,&pCurEntry->REFNumber,&refType,1,1);

    if (*pTHRESHOLD == 0.0) {
      *pTHRESHOLD = pMatchShort->THRESHOLD;
    } else {
      if (*pTHRESHOLD != pMatchShort->THRESHOLD) {
        printf("ERROR: Threshold %f for NUMBER %d is not equal to %f for %s\n",
               pMatchShort->THRESHOLD,pMatchShort->NUMBER,*pTHRESHOLD,rootname);
        exit(1);
      }
    }

    pCurEntry->NUMBER = pMatchShort->NUMBER;
    pCurEntry->AFLAGS = pMatchShort->AFLAGS;
    pCurEntry->BFLAGS = pMatchShort->BFLAGS;
    pCurEntry->class = pMatchShort->class;
    pCurEntry->VFlag = pMatchShort->VFlag;
    pCurEntry->MAGFlag = pMatchShort->MAGFlag;
    pCurEntry->Stdmag = pMatchShort->Stdmag;
    pCurEntry->MAG_ISO = pMatchShort->MAG_ISO;
    pCurEntry->X_IMAGE = pMatchShort->X_IMAGE;
    pCurEntry->Y_IMAGE = pMatchShort->Y_IMAGE;
    pCurEntry->FWHM_WORLD = pMatchShort->FWHM_WORLD * 3600.;
    pCurEntry->dra = pMatchShort->dra;
    pCurEntry->ddec = pMatchShort->ddec;
    pCurEntry->drad = pMatchShort->drad;
    pCurEntry->FLUX_MAX = pMatchShort->FLUX_MAX;
    pCurEntry->seekPosition = pMatchShort->seekPosition;
    pCurEntry->linenumber = match_nrecs;
    pCurEntry->aLength = pMatchShort->aLength;
    pCurEntry->bLength = pMatchShort->bLength;
    pCurEntry->THETA_J2000 = pMatchShort->THETA_J2000;
    pCurEntry->RaSigmaPM = pMatchShort->RaSigmaPM;
    pCurEntry->DecSigmaPM = pMatchShort->DecSigmaPM;

    pCurEntry->flink = NULL;
    ix = (pCurEntry->X_IMAGE * nx) /(1.0 * mosaicWidth);
    iy = (pCurEntry->Y_IMAGE * ny) /(1.0 * mosaicHeight);

    if (ix < 0) {
      ix = 0;
    }

    if (ix >= nx) {
      ix = nx-1;
    }

    if (iy < 0) {
      iy = 0;
    }

    if (iy >= ny) {
      iy = ny-1;
    }

    if (iy <= iyCur) {
      bandlines++;
    }

    if (iy < *poldiy) {
      printf("ERROR: iy decreases from %d to %lf\n", iy, *poldiy);
      exit(1);
    }

    *poldiy = iy;

    /* These records are reversed from the local bin */
    pCurEntry->irec = iy + (ny*ix);

    if (verbose) {
      if ((pCurEntry->AFLAGS == 0) &&
          (pCurEntry->FWHM_WORLD > pFWHM_WORLD_Entry->FWHM_WORLD)) {
        memcpy(pFWHM_WORLD_Entry,pCurEntry,sizeof(MATCH_ENTRY));
      }
    }

    goodlines++;
    pCurEntry++;
  }

  *ppEntryList =  pEntryList;
  *pmatch_alloc = match_alloc;
  *pmatch_nrecs =  match_nrecs;
  *phighdradcount =  highdradcount;
  *pgoodlines = goodlines;
  *pbandlines = bandlines;
}


static void
CreateChains(
  int nx,
  int ny,
  int iyCur,
  int mosaicWidth,
  int mosaicHeight,
  PMATCH_ENTRY pEntryList,
  PCHAIN_ENTRY pChainList,
  PREF_INDEX pRefIndexList,
  PNUMBER_INDEX pNumberIndexList,
  PDRAD dradTable,
  int goodlines,
  double THRESHOLD,
  signed long long *pCurMemory,
  signed long long *pMaxMemory,
  double verbose,
  int *pMaxChainID,
  int *pDuplicateREF,
  int *pDuplicateNUMBER,
  int *pInBandChains,
  int *pAcceptedCaseACount,
  double *pAcceptedCaseAMagnitude,
  int pass, /* Pass 1 adds to the drad vectors, pass 2 does not */
  char *rootname,
  int *pMaxTotalCount,
  int solutionNumber,
  char *qualifier,
  int maxTimeout,
  int skipPVUpdateFlag,
  time_t startTime
) {
  PMATCH_ENTRY pSearchEntry;
  PREF_INDEX pRefIndex;
  PNUMBER_INDEX pNumberIndex;
  PMATCH_ENTRY pCurEntry;
  int queueSize;
  int maxChainID = 0;
  int duplicateREF = 0;
  int duplicateNUMBER = 0;
  int inBandChains = 0;
  int index;
  PCHAIN_ENTRY pCurChain;
  int acceptedCaseACount = 0;
  int acceptedCaseAMagnitude = 0.0;

  /* Now sort the tables so we can look for duplicate ID's */
  for (index = 0; index < goodlines; index++) {
    pSearchEntry = &pEntryList[index];
    pSearchEntry->thisIndex = index;
    pSearchEntry->queueIndex = 0;
    pRefIndex = &pRefIndexList[index];
    pNumberIndex = &pNumberIndexList[index];
    pRefIndex->index = index;
    pRefIndex->REFNumber = pSearchEntry->REFNumber;
    pNumberIndex->index = index;
    pNumberIndex->NUMBER = pSearchEntry->NUMBER;
  }

  qsort(pRefIndexList,goodlines,sizeof(REF_INDEX),RefCompare);
  qsort(pNumberIndexList,goodlines,sizeof(NUMBER_INDEX),NumberCompare);

  /* Now search for duplicate REF ID's */
  for (index = 0; index < (goodlines -1); index++) {
    pRefIndex = &pRefIndexList[index];
    pSearchEntry = &pEntryList[pRefIndex->index];
    pRefIndex = &pRefIndexList[index+1];
    pCurEntry = &pEntryList[pRefIndex->index];

    if (pSearchEntry->REFNumber == pCurEntry->REFNumber) {
      pSearchEntry->REFflink = pCurEntry;
      pCurEntry->REFblink = pSearchEntry;
      duplicateREF++;
      pCurEntry->duplicateFlag |= DUPLICATE_REF;
      pSearchEntry->duplicateFlag |= DUPLICATE_REF;
    }
  }

  /* Now search for duplicate NUMBER ID's */
  for (index = 0; index < (goodlines -1); index++) {
    pNumberIndex = &pNumberIndexList[index];
    pSearchEntry = &pEntryList[pNumberIndex->index];
    pNumberIndex = &pNumberIndexList[index+1];
    pCurEntry = &pEntryList[pNumberIndex->index];

    if (pSearchEntry->NUMBER == pCurEntry->NUMBER) {
      pSearchEntry->NUMBERflink = pCurEntry;
      pCurEntry->NUMBERblink = pSearchEntry;
      duplicateNUMBER++;
      pCurEntry->duplicateFlag |= DUPLICATE_NUMBER;
      pSearchEntry->duplicateFlag |= DUPLICATE_NUMBER;
    }
  }

  /* Find all of the common chains  */

  for (index = 0; index < goodlines; index++) {
    pCurEntry = &pEntryList[index];

    if (pCurEntry->duplicateFlag == 0) {
      // Not in a chain, this is Case A: one GSC2.3 star and one Sextractor
      // object.  Since the pair is uniquely defined, select it.
      if ((pCurEntry->irec % ny) <= iyCur) {
        pCurEntry->selected3 = 1;
        acceptedCaseACount++;
        acceptedCaseAMagnitude += pCurEntry->Stdmag;
      }

      // We will use all accepted values to estimate drad.  Since this entry
      // will be written out immediately, accept it
      if (pass == 1) {
        AddSelectedEntry(nx,ny,iyCur,dradTable,pCurEntry,THRESHOLD,pCurMemory,pMaxMemory);
      }
      continue;
    }

    if (pCurEntry->selected1 != 0) {
      continue;
    }

    pCurChain = &pChainList[maxChainID];
    maxChainID++;
    pCurChain->chainID = maxChainID;

    // Initialize the queue. This queue is necessary because the original
    // implementation using recursion ran off the end of the stack for Baade
    // window plates
    queueSize = 0;
    AddToQueue(pCurEntry,pEntryList,&queueSize);

    while(queueSize > 0) {
      CountRefs(pCurChain,pEntryList,&queueSize,pMaxTotalCount,solutionNumber,rootname,qualifier,maxTimeout,skipPVUpdateFlag,startTime);
    }
  }

  /* At this point, scan through the chains and flag the ones that are in or below the current band */
  for (index = 0; index < maxChainID; index++) {
    pCurChain = &pChainList[index];
    pCurChain->inBand = 1; /* Assume we are in band */
    pCurEntry = pCurChain->chainFlink;

    while (pCurEntry != NULL) {
      if ((pCurEntry->irec % ny) > iyCur) {
        pCurChain->inBand = 0;
        break;
      }

      pCurEntry = pCurEntry->chainFlink;
    }

    if (pCurChain->inBand != 0) {
      inBandChains++;
    }
  }

  *pMaxChainID = maxChainID;
  *pDuplicateREF = duplicateREF;
  *pDuplicateNUMBER = duplicateNUMBER;
  *pInBandChains = inBandChains;

  *pAcceptedCaseACount = acceptedCaseACount;
  *pAcceptedCaseAMagnitude = acceptedCaseAMagnitude;
}


int
main(int argc, char *argv[])
{
  char *argstr;
  char outputname[MAX_INPUTNAME];
  FILE * outputHandle = NULL;
  char REF[MAX_REF];
  char nomatchname[MAX_INPUTNAME];

  File estimate_handle = NULL;
  char estimate_name[MAX_INPUTNAME];
  char qualifier[MAX_BUFFER];
  TableHead estimate_header = NULL;
  PESTIMATE estimate_table = NULL;
  size_t estimate_nrecs = 0;
  int estimate_index;
  PESTIMATE pEstimate;
  double *Stdmag_Table = NULL;
  int Stdmag_index;
  int estimateFilterCount = 0;

  int solutionNumber = 0;
  char blendname[MAX_INPUTNAME];
  FILE * blendHandle = NULL;

  char dradname[MAX_INPUTNAME];
  char regionname[MAX_INPUTNAME];
  char rootname[MAX_INPUTNAME];
  char *rootPtr;
  FILE * dradHandle;
  FILE * regionHandle = NULL;
  int errorFlag = 0;
  int goodlines = 0;
  int allocGoodLines = 0;
  int sumGoodLines = 0;
  int bandlines = 0;
  int maxChainID = 0;
  int duplicateREF = 0;
  int tempDuplicateREF;
  int duplicateNUMBER = 0;
  int tempDuplicateNUMBER;
  int inBandChains = 0;

  int acceptedCaseACount;
  double acceptedCaseAMagnitude;

  int nvals;
  time_t startTime;
  time_t curTime;
  int index;
  long long maxREF = 0;
  PCHAIN_ENTRY pChainList = NULL;
  PMATCH_ENTRY pEntryList = NULL;
  PCHAIN_ENTRY pCurChain;
  MATCH_ENTRY FWHM_WORLD_Entry;
  PMATCH_ENTRY pFWHM_WORLD_Entry = &FWHM_WORLD_Entry;
  PMATCH_ENTRY pCurEntry;
  PMATCH_ENTRY pSearchEntry;
  PMATCH_ENTRY pLastEntry;

  PREF_INDEX pRefIndexList = NULL;
  PNUMBER_INDEX pNumberIndexList = NULL;

  PDRAD dradTable = NULL;
  PDRAD pDrad = NULL;

  char cmdchar;
  int highdradcount = 0;
  int highdradcount2 = 0;
  int tempdradcount;
  int verbose = 0;
  int dumpFlag = 0;
  int skipPVUpdateFlag = 0;
  int regionFlag = 0;
  int totalSingleCount1 = 0;
  int totalREFCount1 = 0;
  int totalNUMBERCount1 = 0;
  int totalBOTHCount1 =0;

  int totalREFCount2 = 0;
  int totalNUMBERCount2 = 0;
  int totalBOTHCount2 =0;
  int maxTotalCount = 0;

  int totalBlendNomatch = 0;
  double scale = 0.0; /* Plate scale in arcsec per pixel */
  char *dbPtr;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int nx = 0;
  int ny = 0;
  int ix;
  int iy;
  int iyCur; /* Current local smoothing bin */
  int irec;
  int totalDradReject = 0;

  int totalHighFluxAccept = 0;
  double THRESHOLD = 0;

  File match_handle_short = NULL;
  File match_handle_full = NULL;
  char match_name[MAX_INPUTNAME];
  TableHead match_header_short = NULL;
  TableHead match_header_full = NULL;
  size_t match_nrecs = 0;
  int match_alloc = 0;
  MATCHFULL match_record_full;
  PMATCHFULL pMatchFull = &match_record_full;
  TblDescriptor match_descriptor_short = NULL;
  TblDescriptor match_descriptor_full = NULL;
  TableRow match_row = NULL;
  double maxDradPixels;
  double maxDradPixels2; /* The smallest blend error after applying the bin drad correction */

  char series[MAX_SERIES_STRING];
  int plateNumber;
  signed long long curMemory = 0;
  signed long long maxMemory = 0;
  double matchRadius =  MATCH_RADIUS_BINS * (pGscBin->bin_size);
  double  searchPixels;
  double bandPixels;
  double oldY_IMAGE = -99.;
  double oldiy = -1;
  int maxChainLength = 0;
  int maxTimeout = 0;
  int tmpMAGFlag;
  int readFlag = 0;
  int FitWCS;

  PCHECK_REPEAT pCheckRefTable = NULL;
  PCHECK_REPEAT pCheckRefEntry = NULL;
  PCHECK_REPEAT pTmpCheckRefTable = NULL;
  int checkRefAllocation = 0;
  int checkRefIndex;
  int checkRefAcceptCount = 0;

  int totalCheckRefReject = 0;
  int checkRefCount = 0;
  double geocentricJD = 0.0;
  double plateepoch;

  /* Loop through the arguments */
  match_name[0] = 0;
  outputname[0] = 0;
  nomatchname[0] = 0;
  blendname[0] = 0;
  dradname[0] = 0;
  rootname[0] = 0;
  regionname[0] = 0;
  estimate_name[0] = 0;
  qualifier[0] = 0;

  memset(pFWHM_WORLD_Entry,0,sizeof(MATCH_ENTRY));

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      if (strlen(match_name) > 0) {
        errorFlag = 1;
        printf("ERROR: list %s is being overwritten by %s\n",match_name,argstr);
      } else {
        strncpy(match_name,argstr,MAX_INPUTNAME-1);
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'v':
        case 'V':
          verbose = 1;
          break;

        case 'j': /* geocentric Julian Day */
        case 'J':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&geocentricJD);
            if (nvals != 1) {
              printf("ERROR: Unable to decode geocentric Julian Day %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'o': /* sextractor region file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(regionname,*++argv,MAX_INPUTNAME-2);
            if (strlen(*argv) >= MAX_INPUTNAME-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            } else {
              regionFlag = 1;
            }
          }
          break;

        case 'p': /* skip flag */
        case 'P':
          skipPVUpdateFlag = 1;
          break;

        case 'd':
        case 'D':
          dumpFlag = 1;
          break;

        case 'e': /* MAG_ISO estimate file name */
        case 'E':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(estimate_name,*++argv,MAX_INPUTNAME-2);
            if (strlen(*argv) >= MAX_INPUTNAME-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'f': /* solution number */
        case 'F':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&solutionNumber);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the solutionNumber %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic width %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* mosaic height */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'q': /* file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(qualifier,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 's':
        case 'S':
          argc--;
          nvals = sscanf(*++argv,"%lf",&scale);
          if (nvals != 1) {
            printf("ERROR: Can not decode the plate scale\n");
            errorFlag = 1;
          }
          break;

        case 't': /* Timeout */
        case 'T':
          argc--;
          nvals = sscanf(*++argv,"%d",&maxTimeout);
          if (nvals != 1) {
            printf("ERROR: Can not decode maximum timeout\n");
            errorFlag = 1;
          }
          break;

        default:
          printf("ERROR: illegal command -%c-",cmdchar);
          errorFlag = 1;
          break;
        }
      }
    }
  }

  if (mosaicWidth == 0) {
    printf("ERROR: No mosaic width specified\n");
    errorFlag = 1;
  }

  if (mosaicHeight == 0) {
    printf("ERROR: No mosaic height specified\n");
    errorFlag = 1;
  }

  nx = XDmagBins(mosaicWidth,mosaicHeight);
  ny = YDmagBins(mosaicWidth,mosaicHeight);

  if (scale == 0.0) {
    printf("ERROR: No plate scale specified\n");
    errorFlag = 1;
  }

  if (geocentricJD == 0.0) {
    printf("ERROR: No geocentric Julian Date specified\n");
    errorFlag = 1;
  }

  /* Attempt to open the list */
  if (match_name[0] == '\0') {
    printf("ERROR: Input file name not specified\n");
    errorFlag = 1;
  }

  strcpy(outputname, match_name);
  dbPtr = strrchr(outputname, '.');

  if (dbPtr == NULL) {
    errorFlag = 1;
    printf("ERROR: expecting match_name %s to contain a '.'\n", match_name);
  } else {
    *dbPtr = '\0';
    strcat(outputname, "_u.db");
  }

  strcpy(blendname, match_name);
  dbPtr = strrchr(blendname, '.');

  if (dbPtr == NULL) {
    errorFlag = 1;
    printf("ERROR: expecting blendname %s to contain a '.'\n", blendname);
  } else {
    *dbPtr = '\0';
    strcat(blendname, "_b.db");
  }

  blendHandle = fopen(blendname, "wt");

  if (blendHandle == NULL) {
    printf("ERROR: Failed to open output file %s\n",blendname);
    errorFlag = 1;
  } else {
    fprintf(blendHandle,"NUMBER\tStdmag\tBlendedmag\n");
    fprintf(blendHandle,"------\t-------\t---------\n");
  }

  rootPtr = strstr(match_name, "match_");

  if (rootPtr == NULL) {
    errorFlag = 1;
    printf("ERROR: expecting match_name %s to begin with 'match_'\n",match_name);
  } else {
    rootPtr += 6;
    strcpy(rootname, rootPtr);
    dbPtr = strrchr(rootname, '.');

    if (dbPtr == NULL) {
      errorFlag = 1;
      printf("ERROR: expecting rootname %s to contain a '.'\n", rootname);
    } else {
      *dbPtr = '\0';
      dbPtr = strstr(rootname, "_tnx");

      if (dbPtr != NULL) {
        *dbPtr = '\0';
      }
    }
  }

  strcpy(dradname, match_name);
  dbPtr = strstr(dradname, "match_");

  if (dbPtr == NULL) {
    errorFlag = 1;
    printf("ERROR: expecting match_name %s to begin with 'match_'\n",match_name);
  } else {
    *dbPtr = '\0';
    strcat(dradname, rootname);

    if (strlen(qualifier) > 0) {
      strcat(dradname, "_");
      strcat(dradname, qualifier);
    }

    strcat(dradname, "_drad.db");
  }

  if (regionFlag) {
    regionHandle = fopen(regionname,"wt");
    if (regionHandle == NULL) {
      printf("ERROR: Failed to open region file %s\n",regionname);
      return 1;
    }
  }

  if (estimate_name[0] != 0) {
    /* Open the no estimate file */
    estimate_handle = Open(estimate_name,"r");

    if (estimate_handle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to find the no-match file %s\n",estimate_name);
    } else {
      if (verbose) {
        printf("Found estimate file %s\n",estimate_name);
      }
    }
  }

  if (ParseFilename2(dradname,series,&plateNumber) == 0) {
    printf("ERROR: Failed to parse %s\n",dradname);
    errorFlag = 1;
  }

  maxDradPixels = MaxDradPixels(series,plateNumber);
  maxDradPixels2 = MaxDradPixels((char *)"none",1);

  if (errorFlag) {
    printf("Usage: filterblended <match_name> -w <width> -h <height> [-v] [-d] [-s <scale (arcsec/pixel)>]\n");
    printf("       where -v is the verbose flag\n");
    printf("             -d dumps the chains\n");
    printf("             -w is the mosaic width in pixels\n");
    printf("             -h is the mosaic height in pixels\n");
    printf("             -e <limiting magnitude estimate file\n");
    printf("             -f <solutionNumber>\n");
    printf("             -r <match radius>\n");
    printf("             -t <timeout in seconds\n");
    printf("             -o <region file>\n");
    printf("             -j Julian date of the plate\n");
    return 1;
  }

  searchPixels = (3600.0 * (matchRadius + (3*pGscBin->bin_size)))/scale;
  plateepoch = jd2ep(geocentricJD)- GSC_EQUINOX;
  bandPixels = mosaicHeight/ny;

  printf("filterblended of %s %s \n Input Filename %s\n Output Filename %s\n nomatch Filename %s\ndrad Filename %s\nPlate scale %12.7f width %d height %d\nSearch Pixels %f Band Pixels %f\n",
         __DATE__,__TIME__,match_name,outputname,nomatchname,dradname,scale,mosaicWidth,mosaicHeight,searchPixels,bandPixels);

  if (searchPixels > bandPixels) {
    printf("WARNING: searchPixels %f is greater than bandPixels %f for %s\n",searchPixels,bandPixels,rootname);
  }

  if ((qualifier[0] == 0)  && (skipPVUpdateFlag == 0)) {
    SetMosaicFitWCS2("FilterblendedTimeout",rootname,solutionNumber,0,readFlag,&FitWCS);
  }

  time(&startTime);

  /* Read in the estimate table */

  if (estimate_handle != NULL) {
    estimate_header = table_header(estimate_handle,TABLE_PARSE);
    if (estimate_header == NULL) {
      printf("ERROR: Failed to read header for %s\n",estimate_name);
      return 1;
    }

    estimate_table = table_loadva(estimate_handle,
                                  &estimate_header,
                                  NULL, /* hbase */
                                  NULL, /* rows */
                                  NULL,
                                  sizeof(ESTIMATE),
                                  &estimate_nrecs,
                                  TblInt,"spatial_bin" ,TblOff(PESTIMATE,spatial_bin),
                                  TblInt,"iMAG_ISO",TblOff(PESTIMATE,iMAG_ISO),
                                  TblDbl,"Stdmag" ,TblOff(PESTIMATE,Stdmag),
                                  TblDbl,"MAG_ISO",TblOff(PESTIMATE,MAG_ISO),
                                  0,"end",0);
    if (estimate_table == NULL) {
      printf("ERROR: Failed to read table for %s\n",estimate_name);
      return 1;
    }

    if (verbose) {
      printf("read %zu records for %s\n",estimate_nrecs,estimate_name);
    }

    if (estimate_nrecs != (MAG_ISO_BINS*MAX_SPATIAL_BINS)) {
      printf("ERROR: read only %zu of expected %d records from %s\n",estimate_nrecs,MAG_ISO_BINS*MAX_SPATIAL_BINS,estimate_name);
      return 1;
    }

    Stdmag_Table = (double *)calloc(MAG_ISO_BINS*MAX_SPATIAL_BINS,sizeof(double));
    CalcMemory(+1,MAG_ISO_BINS*MAX_SPATIAL_BINS*sizeof(double),&curMemory,&maxMemory);

    if (Stdmag_Table == NULL) {
      printf("ERROR: failed to allocate Stdmag_table for %s\n",estimate_name);
      return 1;
    }

    for (estimate_index = 0; estimate_index < estimate_nrecs; estimate_index++) {
      pEstimate = &estimate_table[estimate_index];
      Stdmag_index = pEstimate->iMAG_ISO + (MAG_ISO_BINS *(pEstimate->spatial_bin-1));

      if ((Stdmag_index < 0) || (Stdmag_index >= (MAG_ISO_BINS*MAX_SPATIAL_BINS))) {
        printf("ERROR: Illegal iMAG_ISO %d or spatial_bin %d in record %d of %s\n",
                pEstimate->iMAG_ISO,
                pEstimate->spatial_bin,
                estimate_index,
                estimate_name);
        return 1;
      }

      if (Stdmag_Table[Stdmag_index] != 0) {
        printf("ERROR: duplicate iMAG_ISO %d or spatial_bin %d in record %d of %s\n",
                pEstimate->iMAG_ISO,
                pEstimate->spatial_bin,
                estimate_index,
                estimate_name);
        return 1;
      }

      Stdmag_Table[Stdmag_index] = pEstimate->Stdmag;
    }

    /* We are done with the estimate table.  Clear memory */
    if (estimate_table != NULL) {
      Free(estimate_table);
      estimate_table = NULL;
    }

    if (estimate_header != NULL) {
      table_hdrfree(estimate_header);
      estimate_header = NULL;
    }

    if (estimate_handle != NULL) {
      Close(estimate_handle);
      estimate_handle = NULL;
    }

    pEstimate = NULL;
  }

  /* Attempt to open the list */
  match_handle_short = Open(match_name,"rt");
  if (match_handle_short == NULL) {
    printf("ERROR: Could not open file %s\n",match_name);
    return 1;
  }

  match_handle_full = Open(match_name,"rt");
  if (match_handle_full == NULL) {
    printf("ERROR: Could not open file %s\n",match_name);
    return 1;
  }

  /* Now read in the match file one record at a time */

  match_header_short = table_header(match_handle_short,TABLE_PARSE);
  if (match_header_short == NULL) {
    printf("ERROR: Failed to read header for %s\n",match_name);
    return 1;
  }

  match_descriptor_short = table_create_descrip(&match_nrecs,
                                                TblBuf,"REF"    ,TblOff(PMATCHSHORT,REF),MAX_REF,
                                                TblDbl,"MAG_ISO"  ,TblOff(PMATCHSHORT,MAG_ISO),
                                                TblDbl,"X_IMAGE"  ,TblOff(PMATCHSHORT,X_IMAGE),
                                                TblDbl,"Y_IMAGE"  ,TblOff(PMATCHSHORT,Y_IMAGE),
                                                TblDbl,"Stdmag"  ,TblOff(PMATCHSHORT,Stdmag),
                                                TblDbl,"FWHM_WORLD"  ,TblOff(PMATCHSHORT,FWHM_WORLD),
                                                TblDbl,"drad"  ,TblOff(PMATCHSHORT,drad),
                                                TblDbl,"THRESHOLD"  ,TblOff(PMATCHSHORT,THRESHOLD),
                                                TblDbl,"FLUX_MAX"  ,TblOff(PMATCHSHORT,FLUX_MAX),
                                                TblDbl,"dra"   ,TblOff(PMATCHSHORT,dra),
                                                TblDbl,"ddec"  ,TblOff(PMATCHSHORT,ddec),
                                                TblDbl,"aLength"  ,TblOff(PMATCHSHORT,aLength),
                                                TblDbl,"bLength"  ,TblOff(PMATCHSHORT,bLength),
                                                TblDbl,"THETA_J2000"  ,TblOff(PMATCHSHORT,THETA_J2000),
                                                TblDbl,"dec_1"  ,TblOff(PMATCHSHORT,dec_1),
                                                TblFlt,"RaSigmaPM"  ,TblOff(PMATCHSHORT,RaSigmaPM),
                                                TblFlt,"DecSigmaPM"  ,TblOff(PMATCHSHORT,DecSigmaPM),
                                                TblInt,"NUMBER",TblOff(PMATCHSHORT,NUMBER),
                                                TblInt,"VFlag",TblOff(PMATCHSHORT,VFlag),
                                                TblInt,"MAGFlag",TblOff(PMATCHSHORT,MAGFlag),
                                                TblInt,"AFLAGS",TblOff(PMATCHSHORT,AFLAGS),
                                                TblInt,"BFLAGS",TblOff(PMATCHSHORT,BFLAGS),
                                                TblInt,"class",TblOff(PMATCHSHORT,class),
                                                TblInt,"spatial_bin",TblOff(PMATCHSHORT,spatial_bin),
                                                0,"end",0);
  if (match_descriptor_short == NULL) {
    printf("ERROR: Failed to allocate descriptor for %s\n",match_name);
    return 1;
  }

  table_loadmap(match_header_short,match_descriptor_short);

  match_header_full = table_header(match_handle_full,TABLE_PARSE);
  if (match_header_full == NULL) {
    printf("ERROR: Failed to read header for %s\n",match_name);
    return 1;
  }

  /* Now change our descriptor allocation to the long form */
  match_descriptor_full = table_create_descrip(&match_nrecs,
                                               TblBuf,"REF"    ,TblOff(PMATCHFULL,REF),MAX_REF,
                                               TblDbl,"MAG_ISO"  ,TblOff(PMATCHFULL,MAG_ISO),
                                               TblDbl,"X_IMAGE"  ,TblOff(PMATCHFULL,X_IMAGE),
                                               TblDbl,"Y_IMAGE"  ,TblOff(PMATCHFULL,Y_IMAGE),
                                               TblDbl,"ra_1"  ,TblOff(PMATCHFULL,ra_1),
                                               TblDbl,"dec_1"  ,TblOff(PMATCHFULL,dec_1),
                                               TblDbl,"FWHM_IMAGE"  ,TblOff(PMATCHFULL,FWHM_IMAGE),
                                               TblDbl,"FWHM_WORLD"  ,TblOff(PMATCHFULL,FWHM_WORLD),
                                               TblDbl,"ra_2"  ,TblOff(PMATCHFULL,ra_2),
                                               TblDbl,"dec_2"  ,TblOff(PMATCHFULL,dec_2),
                                               TblDbl,"Stdmag"  ,TblOff(PMATCHFULL,Stdmag),
                                               TblDbl,"color"  ,TblOff(PMATCHFULL,color),
                                               TblDbl,"dra"  ,TblOff(PMATCHFULL,dra),
                                               TblDbl,"ddec"  ,TblOff(PMATCHFULL,ddec),
                                               TblDbl,"drad"  ,TblOff(PMATCHFULL,drad),
                                               TblInt,"NUMBER",TblOff(PMATCHFULL,NUMBER),
                                               TblInt,"AFLAGS",TblOff(PMATCHFULL,AFLAGS),
                                               TblInt,"AFLAGS",TblOff(PMATCHFULL,AFLAGS),
                                               TblInt,"gsc_bin_index",TblOff(PMATCHFULL,gsc_bin_index),
                                               0,"end",0);
  if (match_descriptor_full == NULL) {
    printf("ERROR: Failed to allocate full descriptor for %s\n",match_name);
    return 1;
  }

  table_loadmap(match_header_full,match_descriptor_full);

  /* Allocate statistics table for drad bins */
  dradTable = (PDRAD)calloc(TOTAL_DMAGBINS_NORMAL,sizeof(DRAD));
  if (dradTable == NULL) {
    printf("ERROR: failed to allocate dradTable\n");
    exit(1);
  }

  CalcMemory(1,TOTAL_DMAGBINS_NORMAL*sizeof(DRAD),&curMemory,&maxMemory);

  outputHandle = fopen(outputname,"wt");
  if (outputHandle == NULL) {
    printf("ERROR: Failed to open output file %s\n",outputname);
    return 1;
  }

  fprintf(outputHandle,"NUMBER\tMAG_ISO\tX_IMAGE\tY_IMAGE\tra_1\tdec_1\tFWHM_IMAGE\tFWHM_WORLD\tAFLAGS\tBFLAGS\tREF\tra_2\tdec_2\tStdmag\tcolor\tdra\tddec\tdrad\tgsc_bin_index\n");
  fprintf(outputHandle,"------\t-------\t-------\t-------\t----\t-----\t----------\t----------\t------\t------\t---\t----\t-----\t------\t-----\t---\t----\t----\t-------------\n");

  /* The main loop is on a band of local smoothing bins */
  for (iyCur = 0; iyCur < ny; iyCur++) {
    if (maxTimeout > 0) {
      time(&curTime);
      if ((curTime - startTime) > maxTimeout) {
        if ((qualifier[0] == 0) && (skipPVUpdateFlag == 0)) {
          SetMosaicFitWCS2("FilterblendedTimeout",rootname,solutionNumber,1,readFlag,&FitWCS);
        }
        printf("ERROR: filterblended maxTimeout of %ld/%d seconds maxTotalCount %d reached for %s\n",
                curTime-startTime,
                maxTimeout,
                maxTotalCount,
                rootname);
        exit(1);
      }

    }

    sumGoodLines -= goodlines;
    if (verbose) {
      printf("Processing band %d of %d goodlines %d bandlines %d sumGoodLines %d\n",iyCur,ny,goodlines,bandlines,sumGoodLines);
    }

    ReadBand(nx,ny,iyCur,mosaicWidth,mosaicHeight,&pEntryList,&match_alloc,&match_nrecs,&goodlines,&bandlines,&highdradcount,match_handle_short,match_header_short,match_descriptor_short,&match_row,rootname,pFWHM_WORLD_Entry,&curMemory,&maxMemory,maxDradPixels,plateepoch,scale,&THRESHOLD,verbose,&maxREF,&oldY_IMAGE,&oldiy,Stdmag_Table,&estimateFilterCount);
    sumGoodLines += goodlines;

    if (bandlines <= 0) {
      /* Nothing to process in this band, continue with the next band */
      continue;
    }

    /* Now allocate space for all of the records */
    if (goodlines > allocGoodLines) {
      if (allocGoodLines > 0) {
        CalcMemory(-1,allocGoodLines*(sizeof(CHAIN_ENTRY)+sizeof(REF_INDEX)+sizeof(NUMBER_INDEX)),&curMemory,&maxMemory);

        if (pChainList != NULL) {
          free(pChainList);
          pChainList = NULL;
        }

        if (pRefIndexList != NULL) {
          free(pRefIndexList);
          pRefIndexList = NULL;
        }

        if (pNumberIndexList != NULL) {
          free(pNumberIndexList);
          pNumberIndexList = NULL;
        }
      }

      CalcMemory(1,goodlines*(sizeof(CHAIN_ENTRY)+sizeof(REF_INDEX)+sizeof(NUMBER_INDEX)),&curMemory,&maxMemory);
      pChainList = (PCHAIN_ENTRY)calloc(goodlines,sizeof(CHAIN_ENTRY));
      pRefIndexList = (PREF_INDEX)calloc(goodlines,sizeof(REF_INDEX));
      pNumberIndexList = (PNUMBER_INDEX)calloc(goodlines,sizeof(NUMBER_INDEX));
      allocGoodLines = goodlines;
    }

    if (pEntryList == NULL) {
      printf("ERROR: Failed to allocate pEntryList for %s\n",rootname);
      return 1;
    }

    if (pChainList == NULL) {
      printf("ERROR: Failed to allocate pChainList for %s\n",rootname);
      return 1;
    }

    if (pRefIndexList == NULL) {
      printf("ERROR: Failed to allocate pRefIndexList for %s\n",rootname);
      return 1;
    }

    memset(pChainList,0,goodlines*sizeof(CHAIN_ENTRY));
    CreateChains(nx,ny,iyCur,
                 mosaicWidth,
                 mosaicHeight,
                 pEntryList,
                 pChainList,
                 pRefIndexList,
                 pNumberIndexList,
                 dradTable,
                 goodlines,
                 THRESHOLD,
                 &curMemory,
                 &maxMemory,
                 verbose,
                 &maxChainID,
                 &tempDuplicateREF,
                 &tempDuplicateNUMBER,
                 &inBandChains,
                 &acceptedCaseACount,
                 &acceptedCaseAMagnitude,
                 1,
                 rootname,
                 &maxTotalCount,
                 solutionNumber,
                 qualifier,
                 maxTimeout,
                 skipPVUpdateFlag,
                 startTime);

    if (inBandChains > 0) {
      /* At this point, add the brightest members of the each chain list to the dradTable */
      AddBrightest(nx,ny,iyCur,dradTable,pEntryList,goodlines,pChainList,maxChainID,THRESHOLD,&curMemory,&maxMemory,rootname);

      /* Now compute the DRAD values */
      ComputeDrad(nx,ny,iyCur,dradTable,THRESHOLD,&curMemory,&maxMemory);

      /* Now use the drad table to purge all high drad matches */
      tempdradcount = 0;
      PurgeHighDrad(nx,ny,iyCur,dradTable,pEntryList,&goodlines,&tempdradcount,maxDradPixels2*scale,scale,rootname,regionHandle);
      highdradcount2 += tempdradcount;

      if (tempdradcount > 0) {
        /* Here we have purged some entries, clear the chains and recalculate */
        for (index = 0; index < goodlines; index++) {
          pCurEntry = &pEntryList[index];
          pCurEntry->flink = NULL;
          pCurEntry->REFflink = NULL;
          pCurEntry->REFblink = NULL;
          pCurEntry->NUMBERflink = NULL;
          pCurEntry->NUMBERblink = NULL;
          pCurEntry->chainFlink = NULL;
          pCurEntry->selected1 = 0;
          pCurEntry->selected3 = 0;
          pCurEntry->selected4 = 0;
          pCurEntry->rejected1 = 0;
          pCurEntry->chainID = 0;
          pCurEntry->duplicateFlag = 0;
        }

        memset(pChainList,0,goodlines*sizeof(CHAIN_ENTRY));
        CreateChains(nx,ny,iyCur,
                     mosaicWidth,
                     mosaicHeight,
                     pEntryList,
                     pChainList,
                     pRefIndexList,
                     pNumberIndexList,
                     dradTable,
                     goodlines,
                     THRESHOLD,
                     &curMemory,
                     &maxMemory,
                     verbose,
                     &maxChainID,
                     &tempDuplicateREF,
                     &tempDuplicateNUMBER,
                     &inBandChains,
                     &acceptedCaseACount,
                     &acceptedCaseAMagnitude,
                     2,
                     rootname,
                     &maxTotalCount,
                     solutionNumber,
                     qualifier,
                     maxTimeout,
                     skipPVUpdateFlag,
                     startTime);
      }

      duplicateREF += tempDuplicateREF;
      duplicateNUMBER += tempDuplicateNUMBER;
      acceptedCount[CASE_A] += acceptedCaseACount;
      acceptedMagnitude[CASE_A] += acceptedCaseAMagnitude;

      /* Now count everything up and display it */
      pCurEntry = pEntryList;

      for (index = 0; index < goodlines; index++) {
        if (InBandCandidate(nx,ny,iyCur,pCurEntry,pChainList)) {
          switch(pCurEntry->duplicateFlag) {
          case 0:
            totalSingleCount1++;
            break;
          case DUPLICATE_REF:
            totalREFCount1++;
            break;
          case DUPLICATE_NUMBER:
            totalNUMBERCount1++;
            break;
          case (DUPLICATE_REF|DUPLICATE_NUMBER):
            totalBOTHCount1++;
            break;
          default:
            printf("ERROR: Illegal value %d in duplicateFlag for %s\n",pCurEntry->duplicateFlag,rootname);
            return 1;
          }
        }

        pCurEntry++;
      }

      pCurChain = pChainList;

      for (index = 0; index < maxChainID; index++) {
        if (pCurChain->inBand != 0) {
          if (pCurChain->totalCount > maxChainLength) {
            maxChainLength = pCurChain->totalCount;
          }

          switch(pCurChain->duplicateFlag) {
          case DUPLICATE_REF:
            totalREFCount2++;
            SelectCaseC(pCurChain,rootname);
            break;

          case DUPLICATE_NUMBER:
            SelectCaseB(nx,ny,pCurChain,blendHandle,iyCur,pChainList,rootname);
            totalNUMBERCount2++;
            break;

          case (DUPLICATE_REF|DUPLICATE_NUMBER):
            SelectCaseD(nx,ny,pCurChain,blendHandle,iyCur,pChainList,rootname);
            totalBOTHCount2++;
            break;

          default:
            printf("ERROR: Illegal value %d in duplicateFlag for %s\n",pCurChain->duplicateFlag,rootname);
            return 1;
          }

          if (dumpFlag) {
            DumpChain(pCurChain,rootname);
          }
        }

        pCurChain++;
      }

      /* Now populate the table */
      for (index = 0; index < goodlines; index++) {
        pCurEntry = &pEntryList[index];

        if (((pCurEntry->irec % ny) > iyCur) ||
            (pCurEntry->selected3 == 0)) {
          continue;
        }

        /* Update the FLAGS value
         * While we are here, make a copy of the Sextractor blend flag so that these stars will
         * be rejected as blended objects
         */

        pCurEntry->AFLAGS = pCurEntry->AFLAGS |
          (pCurEntry->VFlag << GSC_VARIABLE_BIT) |
          (pCurEntry->class << GSC_CLASS_BIT) |
          pCurEntry->Afilter_flag;

        tmpMAGFlag = pCurEntry->MAGFlag;

        if (tmpMAGFlag & KEPLER_MAGNITUDE_FLAG) {
          pCurEntry->BFLAGS |= (1 << FILTER_BFLAG_KEPLER);
          tmpMAGFlag &= ~KEPLER_MAGNITUDE_FLAG;
          pCurEntry->BFLAGS = pCurEntry->BFLAGS | (tmpMAGFlag << GSC_MAGNITUDE_FLAG_BIT);
        } else {
          pCurEntry->BFLAGS = pCurEntry->BFLAGS | (pCurEntry->MAGFlag << GSC_MAGNITUDE_FLAG_BIT);
        }
      }

      /* Now flag any stars with a drad greater than DRAD_REJECT_FACTOR times the median */
      for (index = 0; index < goodlines; index++) {
        pCurEntry = &pEntryList[index];

        if (((pCurEntry->irec % ny) > iyCur) || (pCurEntry->selected3 == 0)) {
          continue;
        }

        {
          double dra;
          double ddec;
          double draderror;
          double drad;
          pDrad = &dradTable[pCurEntry->irec];

          if (pDrad->drad_bin_count >= STARS_PER_DMAGBIN) {
            /* We have good drad values for this bin */
            dra = pCurEntry->dra - pDrad->draMedian;
            ddec = pCurEntry->ddec - pDrad->ddecMedian; /* = changed to - on Feb 24, 2010 */
            drad = sqrt(sqr(dra) + sqr(ddec));
            draderror = sqrt(sqr(pDrad->draRMS) + sqr(pDrad->ddecRMS)) * DRAD_REJECT_FACTOR;

            /* Bugfix of Apr 17, 2010: change maxDradPixels to maxDradPixels2 */
            if (draderror > (maxDradPixels2 * scale)) {
              draderror = (maxDradPixels2 * scale);
            }

            if (drad > draderror) {
              pCurEntry->Afilter_flag |= (1 << FILTER_AFLAG_DRAD);
              if  (((pCurEntry->Afilter_flag & FILTER_AMASK_BLEND) == 0) &&
                   ((pCurEntry->BFLAGS & FILTER_BMASK_BLEND) == 0) &&
                   ((pCurEntry->BFLAGS & (1 << FILTER_BFLAG_PMERROR)) == 0)) {
                pDrad->drad_reject_count++;
                totalDradReject++;
              }
            }
          } else {
            /* We have no idea what the drad is for this bin, so reject everything */
            pCurEntry->Afilter_flag |= (1 << FILTER_AFLAG_DRADBIN);
            if  (((pCurEntry->Afilter_flag & FILTER_AMASK_BLEND) == 0) &&
                 ((pCurEntry->BFLAGS & FILTER_BMASK_BLEND) == 0) &&
                 ((pCurEntry->BFLAGS & (1 << FILTER_BFLAG_PMERROR)) == 0)) {
              pDrad->drad_reject_count++;
              totalDradReject++;
            }
          }
        }
      }

      // Modification of Feb 17, 2010.  For any selected entry in band iycur-1,
      // write it out anyways, but be sure to reject all subsequent instances of
      // REFNumber

      for (index = 0; index < goodlines; index++) {
        pCurEntry = &pEntryList[index];

        if ((pCurEntry->selected3 != 0) &&
            (pCurEntry->rejected1 == 0) &&
            ((pCurEntry->irec % ny) < iyCur)) {

          /* Traverse the REF queue and mark all entries rejected */
          pSearchEntry = pCurEntry->REFflink;
          while (pSearchEntry) {
            pSearchEntry->rejected1 = 1;
            pSearchEntry = pSearchEntry->REFflink;
          }

          pSearchEntry = pCurEntry->REFblink;
          while (pSearchEntry) {
            pSearchEntry->rejected1 = 1;
            pSearchEntry = pSearchEntry->REFblink;
          }

          pSearchEntry = pCurEntry->NUMBERflink;
          while (pSearchEntry) {
            pSearchEntry->rejected1 = 1;
            pSearchEntry = pSearchEntry->NUMBERflink;
          }

          pSearchEntry = pCurEntry->NUMBERblink;
          while (pSearchEntry) {
            pSearchEntry->rejected1 = 1;
            pSearchEntry = pSearchEntry->NUMBERblink;
          }

          /* Now make sure we have not written this one out yet */
          for (checkRefIndex = 0; checkRefIndex < checkRefCount; checkRefIndex++) {
            pCheckRefEntry = &pCheckRefTable[checkRefIndex];

            if (pCurEntry->REFNumber == pCheckRefEntry->REFNumber) {
              printf("WARNING: Attempt to write out REF %lld twice for %s\n",pCurEntry->REFNumber,rootname);
              pCurEntry->rejected1 = 1;
              break;
            }

            if (pCurEntry->NUMBER == pCheckRefEntry->NUMBER) {
              printf("ERROR: Attempt to write out NUMBER %d twice for %s\n",pCurEntry->NUMBER,rootname);
              exit(1);
            }
          }

          if (pCurEntry->rejected1 != 0) {
            continue;
          }

          /* Add this one to the table */
          if ((checkRefCount+1) >= checkRefAllocation) {
            checkRefAllocation += 100;
            pTmpCheckRefTable = realloc(pCheckRefTable,checkRefAllocation * sizeof(CHECK_REPEAT));

            if (pTmpCheckRefTable == NULL) {
              printf("ERROR allocating pCheckRefTable of size %d\n",checkRefAllocation);
              exit(1);
            }

            CalcMemory(1,100 * sizeof(CHECK_REPEAT),&curMemory,&maxMemory);
            pCheckRefTable = pTmpCheckRefTable;
            pTmpCheckRefTable = NULL;
          }

          pCheckRefEntry = &pCheckRefTable[checkRefCount];
          memset(pCheckRefEntry,0,sizeof(CHECK_REPEAT));
          pCheckRefEntry->REFNumber = pCurEntry->REFNumber;
          pCheckRefEntry->NUMBER = pCurEntry->NUMBER;
          checkRefCount++;
        }
      }

      /* Write out the results */

      for (index = 0; index < goodlines; index++) {
        pCurEntry = &pEntryList[index];

        if ((pCurEntry->selected3 == 0) ||
            (pCurEntry->rejected1 != 0))
            {
              continue;
            }

        if (InBandCandidate(nx,ny,iyCur,pCurEntry,pChainList) == 0) {
          if ((pCurEntry->irec % ny) >= iyCur) {
            continue;
          } else {
            checkRefAcceptCount++;
          }
        }

        /* Now make sure we have not written this one out yet */
        for (checkRefIndex = 0; checkRefIndex < checkRefCount; checkRefIndex++) {
          pCheckRefEntry = &pCheckRefTable[checkRefIndex];

          if (pCurEntry->REFNumber == pCheckRefEntry->REFNumber) {
            if (pCheckRefEntry->refCount1 > 0) {
              printf("WARNING: Attempt (2) to write out REF %lld twice for %s\n",pCurEntry->REFNumber,rootname);
              pCurEntry->rejected1 = 1;
              break;
            }

            pCheckRefEntry->refCount1++;
            break;
          }

          if (pCurEntry->NUMBER == pCheckRefEntry->NUMBER) {
            if (pCheckRefEntry->numberCount1 > 0) {
              printf("ERROR: Attempt (2) to write out NUMBER %d twice for %s\n",pCurEntry->NUMBER,rootname);
              exit(1);
            }

            pCheckRefEntry->numberCount1++;
            break;
          }
        }

        if (pCurEntry->rejected1 != 0) {
          continue;
        }

        Seek(match_handle_full,pCurEntry->seekPosition,SEEK_SET);
        memset(pMatchFull,0,sizeof(MATCHFULL));

        if(ReadMatchFull(pMatchFull,match_handle_full,match_header_full,match_descriptor_full,&match_row,rootname) == 0) {
          printf("ERROR reading at position %ld %d %s %s\n",pCurEntry->seekPosition,errno,strerror(errno),match_name);
          return 1;
        }

        GetREF(pCurEntry->REFNumber,REF,1,1);

        if ((pMatchFull->NUMBER != pCurEntry->NUMBER) || (strcmp(pMatchFull->REF,REF) != 0)) {
          printf("ERROR: full %d %s and short %d %s NUMBER or REF disagree\n",pMatchFull->NUMBER,pMatchFull->REF,pCurEntry->NUMBER,REF);
          return 1;
        }

        /* Create the expanded FLAGS variable */
        pMatchFull->AFLAGS = pCurEntry->AFLAGS |
          (pCurEntry->class << GSC_CLASS_BIT) |
          (pCurEntry->VFlag << GSC_VARIABLE_BIT) |
          pCurEntry->Afilter_flag;

        tmpMAGFlag = pCurEntry->MAGFlag;

        if (tmpMAGFlag & KEPLER_MAGNITUDE_FLAG) {
          pMatchFull->BFLAGS |= (1 << FILTER_BFLAG_KEPLER);
          tmpMAGFlag &= ~KEPLER_MAGNITUDE_FLAG;
          pMatchFull->BFLAGS = pCurEntry->BFLAGS | (tmpMAGFlag << GSC_MAGNITUDE_FLAG_BIT);
        } else {
          pMatchFull->BFLAGS = pCurEntry->BFLAGS | (pCurEntry->MAGFlag << GSC_MAGNITUDE_FLAG_BIT);
        }

        if ((pMatchFull->AFLAGS & (1 << FILTER_AFLAG_BLEND_NOMATCH)) == 0) {
          fprintf(outputHandle,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",
                  pMatchFull->NUMBER,
                  pMatchFull->MAG_ISO,
                  pMatchFull->X_IMAGE,
                  pMatchFull->Y_IMAGE,
                  pMatchFull->ra_1,
                  pMatchFull->dec_1,
                  pMatchFull->FWHM_IMAGE,
                  pMatchFull->FWHM_WORLD,
                  pMatchFull->AFLAGS,
                  pMatchFull->BFLAGS,
                  pMatchFull->REF,
                  pMatchFull->ra_2,
                  pMatchFull->dec_2,
                  pMatchFull->Stdmag,
                  pMatchFull->color,
                  pMatchFull->dra,
                  pMatchFull->ddec,
                  pMatchFull->drad,
                  pMatchFull->gsc_bin_index);
        }
      }
    }

    /* At this point, discard all of the points we just wrote */
    index = 0;
    pCurEntry = pEntryList;

    while (index < goodlines) {
      if ((pCurEntry->rejected1 != 0) ||
          (InBandCandidate(nx,ny,iyCur,pCurEntry,pChainList)))
      {
        /* Toss this image */
        if (pCurEntry->rejected1 != 0) {
          totalCheckRefReject++;
        }

        goodlines--;
        if (index != goodlines) {
          pLastEntry = &pEntryList[goodlines];
          memcpy(pCurEntry,pLastEntry,sizeof(MATCH_ENTRY));
        }
      } else {
        /* Keep this image for further consideration */
        index++;
        pCurEntry++;
      }
    }

    /* Reset the entries */
    for (index = 0; index < goodlines; index++) {
      pCurEntry = &pEntryList[index];
      pCurEntry->flink = NULL;
      pCurEntry->REFflink = NULL;
      pCurEntry->REFblink = NULL;
      pCurEntry->NUMBERflink = NULL;
      pCurEntry->NUMBERblink = NULL;
      pCurEntry->chainFlink = NULL;
      pCurEntry->selected1 = 0;
      pCurEntry->selected3 = 0;
      pCurEntry->selected4 = 0;
      pCurEntry->rejected1 = 0;
      pCurEntry->chainID = 0;
      pCurEntry->duplicateFlag = 0;
      pCurEntry->AFLAGS &= ~(1 << FILTER_AFLAG_BLEND_NOMATCH);
      pCurEntry->AFLAGS &= ~(pCurEntry->Afilter_flag);
      pCurEntry->Afilter_flag = 0;
    }
  }

  CalcMemory(-1,allocGoodLines*(sizeof(CHAIN_ENTRY)+sizeof(REF_INDEX)+sizeof(NUMBER_INDEX)),&curMemory,&maxMemory);

  if (pChainList != NULL) {
    free(pChainList);
    pChainList = NULL;
  }

  if (pRefIndexList != NULL) {
    free(pRefIndexList);
    pRefIndexList = NULL;
  }

  if (pNumberIndexList != NULL) {
    free(pNumberIndexList);
    pNumberIndexList = NULL;
  }

  if (sumGoodLines == 0) {
    printf("ERROR: No lines read in from %s\n",match_name);
    return 1;
  }

  if (verbose) {
    printf("Max full width half maximum is %f for %lld\n",
           pFWHM_WORLD_Entry->FWHM_WORLD,pFWHM_WORLD_Entry->REFNumber);
    printf("match_alloc %d, match_nrecs %zu, good lines %d\n",match_alloc,match_nrecs,sumGoodLines);
  }

  if (verbose) {
    printf("Total pairs: %d, single: %d duplicate REF: %d duplicate NUMBER: %d, both: %d\n",
           totalSingleCount1+totalREFCount1+totalNUMBERCount1+totalBOTHCount1,
           totalSingleCount1,totalREFCount1,totalNUMBERCount1,totalBOTHCount1);
    printf("Total chains %d\n",maxChainID);
  }

  if (verbose) {
    printf("Total chains %d duplicate REF chains: %d duplicate NUMBER chains: %d, both: %d\n",
           totalREFCount2+totalNUMBERCount2+totalBOTHCount2,
           totalREFCount2,totalNUMBERCount2,totalBOTHCount2);
  }

  /* Write out the drad file */
  dradHandle = fopen(dradname,"wt");
  if (dradHandle == NULL) {
    printf("ERROR: Failed to open drad file %s\n",dradname);
    return 1;
  }

  fprintf(dradHandle,"nx\tny\tix\tiy\tdrad_bin_count\tdrad_bin_size\tdrad_reject_count\tdraMedian\tdraRMS\tddecMedian\tddecRMS\n");
  fprintf(dradHandle,"--\t--\t--\t--\t--------------\t-------------\t-----------------\t---------\t------\t----------\t-------\n");

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      irec = iy + (ix * ny);
      pDrad = &dradTable[irec];
      fprintf(dradHandle,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",
              nx,
              ny,
              ix,
              iy,
              pDrad->drad_bin_count,
              pDrad->drad_bin_size,
              pDrad->drad_reject_count,
              pDrad->draMedian,
              pDrad->draRMS,
              pDrad->ddecMedian,
              pDrad->ddecRMS);
    }
  }

  fclose(dradHandle);
  fclose(outputHandle);
  fclose(blendHandle);

  if (regionHandle != NULL) {
    fclose(regionHandle);
  }

  if (match_header_short != NULL) {
    table_hdrfree(match_header_short);
  }

  if (match_header_full != NULL) {
    table_hdrfree(match_header_full);
  }

  if (match_handle_short != NULL) {
    Close(match_handle_short);
  }

  if (match_handle_full != NULL) {
    Close(match_handle_full);
  }

  if (match_descriptor_short != NULL) {
    Free(match_descriptor_short);
  }

  if (match_descriptor_full != NULL) {
    Free(match_descriptor_full);
  }

  if (match_row != NULL) {
    table_rowfree(match_row);
  }

  CalcMemory(-1,match_alloc*sizeof(MATCH_ENTRY),&curMemory,&maxMemory);

  if (pEntryList != NULL) {
    Free(pEntryList);
  }

  for (iy = 0; iy < ny; iy++) {
    for (ix = 0; ix < nx; ix++) {
      irec = iy + (ix * ny);
      pDrad = &dradTable[irec];
      if (pDrad->vectorAlloc > 0) {
        free(pDrad->draVector);
        free(pDrad->ddecVector);
        pDrad->draVector = NULL;
        pDrad->ddecVector = NULL;
        CalcMemory(-1,2*pDrad->vectorAlloc * sizeof(double),&curMemory,&maxMemory);
        pDrad->vectorCount = 0;
        pDrad->vectorAlloc = 0;
      }
    }
  }

  CalcMemory(-1,TOTAL_DMAGBINS_NORMAL*sizeof(DRAD),&curMemory,&maxMemory);
  if (dradTable != NULL) {
    free(dradTable);
  }

  if (Stdmag_Table != NULL) {
    CalcMemory(-1,MAG_ISO_BINS*MAX_SPATIAL_BINS*sizeof(double),&curMemory,&maxMemory);
    free(Stdmag_Table);
  }

  if (pCheckRefTable != NULL) {
    free(pCheckRefTable);
    CalcMemory(-1,checkRefAllocation*sizeof(CHECK_REPEAT),&curMemory,&maxMemory);
    pCheckRefTable = NULL;
  }

  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %ld seconds for %d input stars for %s\n",curTime,sumGoodLines,rootname);

  if (verbose) {
    printf("Stars rejected because drad is greater than FWHM %d %d for %s\n",highdradcount,highdradcount2,rootname);
    printf("Maximum REF: %lld\n",maxREF);
    printf("Duplicate GSC2ID %d, duplicate Sextractor id %d for %s\n",duplicateREF,duplicateNUMBER,rootname);
  }

  printf(
    "Blended stars accepted as normal stars: %d; rejected: %d\n"
    "  DRAD rejections: %d\n"
    "  Unmatched blends: %d\n"
    "  High Flux Accept: %d\n"
    "  maxChainLength: %d\n"
    "  estimateFilterCount: %d\n"
    "  maxTotalCount: %d\n"
    "  checkRefCount: %d\n"
    "  totalCheckRefReject: %d\n"
    "  maxMemory: %lld; curMemory: %lld\n",
    blendedAccept,
    blendedReject,
    totalDradReject,
    totalBlendNomatch,
    totalHighFluxAccept,
    maxChainLength,
    estimateFilterCount,
    maxTotalCount,
    checkRefCount,
    totalCheckRefReject,
    maxMemory,
    curMemory
  );

  for (index = 0; index < CASE_MAX; index++) {
    if (index != CASE_TOTAL) {
      rejectedCount[CASE_TOTAL] += rejectedCount[index];
      rejectedMagnitude[CASE_TOTAL] += rejectedMagnitude[index];
    }

    if ((index >= CASE_A) && (index < CASE_TOTAL)) {
      acceptedCount[CASE_TOTAL] += acceptedCount[index];
      acceptedMagnitude[CASE_TOTAL] += acceptedMagnitude[index];
    }

    if (rejectedCount[index] != 0) {
      rejectedMagnitude[index] = rejectedMagnitude[index]/rejectedCount[index];
    }

    if (acceptedCount[index] != 0) {
      acceptedMagnitude[index] = acceptedMagnitude[index]/acceptedCount[index];
    }

    printf("%s Rejected %5d Ave Mag %6.2f; Accepted %5d Ave Mag %6.2f\n",
           caseString[index],
           rejectedCount[index],
           rejectedMagnitude[index],
           acceptedCount[index],
           acceptedMagnitude[index]);
  }

  return EXIT_SUCCESS;
}


/*
 *  Nov 23, 2007 Edward J. Los - Initial version
 *  Dec 12, 2007 Edward J. Los -         (a) flag all case b, c, and d stars for eventual rejection
 *                                       (b) do not used blended magnitudes
 *                                       (c) Add astrometric uncertainty to the search radius
 *                                       (d) flag GSC2.3 non-stellar stars and variable stars for eventual
 *                                           rejection and remove the "class" and "VFlag" columns
 *  Dec 20, 2007 Edward J. Los -          Add magnitudeFlag (MAGFlag) to FLAGS
 *                                        Accept blended stars as normal stars if the magnitude difference is too small
 *                                        Validate the column headers.
 *  Feb 20, 2008 Edward J. Los -          Rewrite search algorithm for performance
 *  Mar 27, 2008 Edward J. Los -          Convert to starbase format
 *  Apr 29, 2008 Edward J. Los -          Correct selection of stars for case DB
 *                                        Make case DB conform to case B: If the combined magnitude does not change
 *                                        significantly, do not call the star blended.
 *  May  2, 2008 Edward J. Los - V3.4.2   Perform additional filtering on objects with high drad.
 *  May 28, 2008 Edward J. Los -          Move FLAGS definition to pipelineutils.h
 *  May 29, 2008 Edward J. Los - V3.4.4   Change MIN_BLEND_MAGNITUDE to 0.05  on
 *                                        recommendation of Sumin Tang.
 *  Jul  1, 2008 Edward J. Los - V3.4.7   Set FILTER_FLAG_BLEND_COPY to reject objects that
 *                                        Sextractor believes are blended
 *                                        Change the search radius from FWHM_WORLD + 1 pixel
 *                                                               to (FWHM_WORLD/2) + 5 pixels
 *                                        Set FILTER_FLAG_BLEND_NOMATCH to include objects that
 *                                        would otherwise be plotted as objects with bad limiting
 *                                        magnitudes.  Write these objects to a different file to avoid
 *                                        issues with duplicate Sextractor NUMBER fields.
 *                                        Perform drad analysis only for no bad FLAGS and for FLUX_MAX > 3 * THRESHOLD
 *                                                    Perform 3 sigma clipping to get the median drad per bin.
 *                                        Do not do drad filtering if no median drad can be calculated for the bin.
 *  Aug  8, 2008 Edward J. Los - V3.4.9   Set a drad ceiling at MAX_DRAD_PIXELS
 *                                        Do not set the FILTER_FLAG_DRAD for blended stars
 *                                        Write the match_<plate>_tnx_b.db file for blended stars so we can correct
 *                                        the blend magnitude after
 *  Sep 22, 2008 Edward J. Los - V3.4.12  Change MIN_BLEND_MAGNITUDE to 0.1 on recommendation of scanner meeting
 *  Sep 26, 2008 Edward J. Los - V3.4.12  Check and flag blended stars for high drad, but still do not count
 *                                        them in the drad statistics file.
 *  Oct  1, 2008 Edward J. Los - V3.4.13  Improve error reporting
 *                                        Read in the saturationFlux from the *_ccmap file, write it out to the
 *                                        drad summary file and use this value to exclude setting drad flags for
 *                                        saturated images.
 *  Oct 10, 2008 Edward J. Los -          Make the ccmap file optional.
 *  Nov 24, 2008 Edward J. Los -          Back out saturation flux checks
 *  Nov 25, 2008 Edward J. Los -          Replace MAX_DRAD_PIXELS and PIXEL_RADIUS_FACTOR with MaxDradPixels() function
 *  Jan 27, 2009 Edward J. Los -          Split flags into AFLAGS and BFLAGS
 *  Feb  6, 2009 Edward J. Los -          Use new blended star algorithm to replace the Sextractor blended flag.
 *                                        Ignore the Sextractor blend flag in favor of a more restrictive blend definition.
 *  Feb 10, 2009 Edward J. Los -          Use the new blend detection algorithm no longer based on FWHM_WORLD
 *  Feb 16, 2009 Edward J. Los -          Add gsc_bin_index to output
 *  Feb 16, 2009 Edward J. Los -          Use size_t for the number of records in a table to avoid crashes on 64 bit systems
 *                                        when the table size exceeds 2GB
 *                                        Make drad filtering local smoothing bin dependent instead of annular bin
 *                                        dependent, but do not make the final adjustment of positions
 *                                        because of high rms values.
 *                                        Back out saturationFlux (ccmap) code
 * Feb 27, 2009 Edward J. Los -           Add CalcMemory
 * Mar  3, 2009 Edward J. Los -           Operate in Y_IMAGE raster mode for performance
 * Mar 13, 2009 Edward J. Los -           Correct AddSelectedEntry comparison error: images matching the same GSC star will all have
 *                                        the same Stdmag
 *                                        If magnitudes are equal in Case B, C, and D, select the lowest REF
 *                                        Rewrite chain creation to avoid a recursion failure running off the stack
 *                                        Add timeout support.
 * May  1, 2009 Edward J. Los -           Tighten the drad error radius after applying the bin drad correction.
 *                                        Correct bug in checking the drad error radius
 * May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 * May 25, 2009 Edward J. Los - Make CountRefs more efficient
 * Jun 15, 2009 Edward J. Los - Add an estimate of the limiting magnitude for this plate
 *                              (NOTE: this code is now disabled with a change to dasch_match.csh)
 * Jun 30, 2009 Edward J. Los - Add REMOVE_NOMATCH to turn off FILTER_AFLAG_BLEND_NOMATCH functionality
 *                              To improve performance, remove FILTER_AFLAG_BLEND_NOMATCH objects early if the brightest star
 *                              is 5 mags above the potential matching star
 * Jul 21, 2009 Edward J. Los - Demote the searchPixels > bandPixels error to a warning.  We may miss some high proper motion
 *                              stars.
 * Aug 31, 2009 Edward J. Los - Covert GSC and Kepler ID conversion to photometry utilities
 * Sep 28, 2009 Edward J. Los - Add maxTotalCount and other infinite loop debugging information.
 * Dec 15, 2009 Edward J. Los - Flag plates that time out in the database.
 * Feb 17, 2010 Edward J. Los - Use rejected1 to write out images in band iyCur-1 even if they are part of a out-of-band chain
 *                              Correct seekPosition data type
 * Feb 24, 2010 Edward J. Los - Double the maximum objects in a blend
 *                              Correct calculation of ddec when setting the FILTER_AFLAG_DRAD flag
 * Feb 25, 2010 Edward J. Los - Start checking execution times when blends get large
 * Mar 16, 2010 Edward J. Los - Change duplicate REF error to a warning when change of Feb 17 above is active.
 * Apr 17, 2010 Edward J. Los - Impose a tighter draderror criterion for ac, am, ax, and ca plates.
 * May 27, 2010 Edward J. Los - Add option to write a region file
 *                              Correct a bug: Failure to copy THETA_J2000
 *                              For high drad bins, check to tighter tolerances.  The bin will most likely be rejected anyways.
 *                              In PurgeHighDrad do not add draderror to aLengthMax and bLengthMax.  Take the minimum of the two.
 * Apr 12, 2011 Edward J. Los - Remove unused keplerFlag
 * Apr 15, 2012 Edward J. Los - add option to skip the update of the InaccuratePV flag
 * Jul  9, 2012 Edward J. Los - Remove THRESHOLD_FACTOR for the new defect filter of  Fri 6/29/12 9:43 AM
 * Jul  3, 2013 Edward J. Los - Correct CheckBlend angle error in PurgeHighDrad   (Assume that the CheckBlend in ReadBand is o.k.)
 * Aug 27, 2013 Edward J. Los   Remove matchRadius, replacing it with three GSC bins
 * Sep  2, 2013 Edward J. Los - Make the CheckBlend in ReadBand conform to the routine definition.
 *                            - Expand CheckBlend size by the error in proper motion, do not use high proper motion error drad objects to compute blend angles.
 * May  4, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 */