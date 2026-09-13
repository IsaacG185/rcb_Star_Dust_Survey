// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// setScampRMS
//
// This tool has evolved into a generic database-checking utility.
//
// Commmon arguments:
//
//   -p <platename> -- specifies plate of interest
//   -e <solnum>    -- specifies solution number of interest
//
// Modes:
//
//   -n <PolyRefStars> -r <PolyRARMS> -d <PolyDecRMS> -- "original"; save SCAMP info to database
//   -c                                               -- "check new algorithm success"
//   -i                                               -- "Check success of new procedure with previous solution"
//   -s                                               -- "Check success of new procedure with stale solution"
//   -f                                               -- "check and set new algorithm failure"
//   -a                                               -- "declare new algorithm astrometry.wcs success"
//   -g                                               -- "check readiness to run the SCAMP USNO procedure"
//   -q                                               -- "declare new algorithm SCAMP USNO success"
//   -b                                               -- "declare new algorithm SCAMP multiple exposure astrometry.wcs success"
//   -k                                               -- "check readiness to run SCAMP UCAC procedure"
//   -h <hdr> -n <prefstars> -r <pRARMS> -d <pDecRMS> -- "declare new algorithm SCAMP UCAC success providing that the new header is acceptable"
//   -l                                               -- "declare new algorithm SCAMP UCAC success for multiple exposures"
//   -m <header>                                      -- "check header output for scamp for inaccuracy"
//   -y                                               -- "reset full processing status completely"
//
// Virtually all of these modes are still used in the current pipeline.

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>
#include <libwcs/fitswcs.h>

// sigh, relying on un-headered libwcs functions ...
extern struct WorldCoor *GetFITSWCS (
  char *filename,
  char *header,
  int verbose,
  double *cra,
  double *cdec,
  double *dra,
  double *ddec,
  double *secpix,
  int *wp,
  int *hp,
  int *sysout,
  double *eqout
);

#include "scandb.h"
#include "pipelineutils.h"

static int verbose = 0;		/* verbose/debugging flag */
static double eqim = 0.0;
static int sysim = 0;

#define MAX_BUFFER 512

/* Parameter validation */
#define ARG_POLYREFSTARS              1  /* -n */
#define ARG_POLYRARMS                 2  /* -r */
#define ARG_POLYDECRMS                4  /* -d */
#define ARG_FULLSUCCESS               8  /* -c */
#define ARG_FULLFAILURE              16  /* -f */
#define ARG_FULLASTROMETRY           32  /* -a */
#define ARG_FULLSCAMPUSNO            64  /* -b */
#define ARG_FULLSCAMPUCAC           128  /* -h */
#define ARG_FULLRESETFAIL           256  /* -y */
#define ARG_FULLCHECKASTROM         512  /* -g */
#define ARG_FULLCHECKUSNO          1024  /* -k */
#define ARG_CHECKHEADER            2048  /* -m */
#define ARG_FULLCHECKASTROMSECOND  4096  /* -q */
#define ARG_FULLSCAMPUCACSECOND    8192  /* -l */
#define ARG_CHECKPREV             16384  /* -i */
#define ARG_CHECKSTALE            32768  /* -s */

#define FUNCTION_ORIGINAL               (ARG_POLYREFSTARS|ARG_POLYRARMS|ARG_POLYDECRMS)
#define FUNCTION_CHECKNEW               ARG_FULLSUCCESS
#define FUNCTION_CHECKPREV              ARG_CHECKPREV
#define FUNCTION_CHECKSTALE              ARG_CHECKSTALE
#define FUNCTION_SETFULLFAIL            ARG_FULLFAILURE
#define FUNCTION_SETFULLASTROMETRY      ARG_FULLASTROMETRY
#define FUNCTION_SETFULLSCAMPUSNO       ARG_FULLSCAMPUSNO
#define FUNCTION_SETFULLSCAMPUCAC       (ARG_FULLSCAMPUCAC|ARG_POLYREFSTARS|ARG_POLYRARMS|ARG_POLYDECRMS)
#define FUNCTION_RESETFULLFAIL          ARG_FULLRESETFAIL
#define FUNCTION_CHECKFULLASTROM        ARG_FULLCHECKASTROM
#define FUNCTION_CHECKFULLUSNO          ARG_FULLCHECKUSNO
#define FUNCTION_CHECKHEADER            ARG_CHECKHEADER
#define FUNCTION_CHECKFULLASTROMSECOND  ARG_FULLCHECKASTROMSECOND
#define FUNCTION_SETFULLSCAMPUCACSECOND ARG_FULLSCAMPUCACSECOND

#define DEFAULT_POLYREFSTARS -999999

static int functionlist[] = {
  FUNCTION_ORIGINAL,
  FUNCTION_CHECKNEW,
  FUNCTION_CHECKPREV,
  FUNCTION_CHECKSTALE,
  FUNCTION_SETFULLFAIL,
  FUNCTION_SETFULLASTROMETRY,
  FUNCTION_SETFULLSCAMPUSNO,
  FUNCTION_SETFULLSCAMPUCAC,
  FUNCTION_SETFULLSCAMPUCACSECOND,
  FUNCTION_RESETFULLFAIL,
  FUNCTION_CHECKFULLASTROM,
  FUNCTION_CHECKFULLASTROMSECOND,
  FUNCTION_CHECKFULLUSNO,
  FUNCTION_CHECKHEADER
};

static char *functionNames[] = {
  "original",
  "-c checknew",
  "-i checkprev",
  "-i checkstale",
  "-f setfullfail",
  "-a setfullastrometry",
  "-b setfullscampusno",
  "-h setfullscampucac",
  "-l setfullscampucacsecond",
  "-y resetfullfail",
  "-g checkfullastrom",
  "-q checkfullastromsecond",
  "-k checkfullusno",
  "-m checkheader",
  "UNKNOWN"
};

static int functionlistsize = sizeof(functionlist) / sizeof(int);

typedef struct _checkCommon {
  char platename[MAX_BUFFER];
  char headername[MAX_BUFFER];
  double polyRARMS;
  double polyDecRMS;
  double oldPolyRARMS;
  double oldPolyDecRMS;
  int solutionNumber;
  int polyRefStars;
  int oldPolyRefStars;
  int prevFitWCS;
  int oldFitWCS;
  int newFitWCS;
} CHECKCOMMON, *PCHECKCOMMON;


static int
CheckScampHeader(PCHECKCOMMON pCheckCommon, int declareError)
{
  struct WorldCoor *wcs = NULL;
  char *header = NULL;
  double cra;
  double cdec;
  double dra;
  double ddec;
  double secpix;
  int wp;
  int hp;

  if ((header = GetFITShead (pCheckCommon->headername, verbose)) == NULL) {
    printf(
      "ERROR: setScampRMS line %d Unable to read header  PolyRefStars %d->%d PolyRARMS %.2f->%.2f PolyDecRMS %.2f->%.2f for %s s%d\n",
      __LINE__,
      pCheckCommon->oldPolyRefStars,
      pCheckCommon->polyRefStars,
      pCheckCommon->oldPolyRARMS,
      pCheckCommon->polyRARMS,
      pCheckCommon->oldPolyDecRMS,
      pCheckCommon->polyDecRMS,
      pCheckCommon->platename,
      pCheckCommon->solutionNumber
    );
    return -1;
  }

  wcs = GetFITSWCS (pCheckCommon->headername, header, verbose, &cra, &cdec, &dra, &ddec,
                    &secpix, &wp, &hp, &sysim, &eqim);
  if (nowcs (wcs)) {
    wcsfree (wcs);
    wcs = NULL;
    printf("ERROR: setScampRMS line %d Unable to get wcs PolyRefStars %d->%d PolyRARMS %.2f->%.2f PolyDecRMS %.2f->%.2f for %s s%d\n",
           __LINE__,
           pCheckCommon->oldPolyRefStars,
           pCheckCommon->polyRefStars,
           pCheckCommon->oldPolyRARMS,
           pCheckCommon->polyRARMS,
           pCheckCommon->oldPolyDecRMS,
           pCheckCommon->polyDecRMS,
           pCheckCommon->platename,
           pCheckCommon->solutionNumber);
    if (header != NULL) {
      free(header);
    }
    return -1;
  } else {
    if (wcs->pvfail) {
      printf("WARNING: setScampRMS line %d Significant inaccuracy PolyRefStars %d->%d PolyRARMS %.2f->%.2f PolyDecRMS %.2f->%.2f for %s s%d\n",
             __LINE__,
             pCheckCommon->oldPolyRefStars,
             pCheckCommon->polyRefStars,
             pCheckCommon->oldPolyRARMS,
             pCheckCommon->polyRARMS,
             pCheckCommon->oldPolyDecRMS,
             pCheckCommon->polyDecRMS,
             pCheckCommon->platename,
             pCheckCommon->solutionNumber);

      wcsfree (wcs);
      wcs = NULL;
      if (declareError != 0) {
        pCheckCommon->newFitWCS &= ~FITWCS_FULLFAILEDSCAMPUCAC;
        pCheckCommon->newFitWCS |= FITWCS_FULLFAILEDTOLERANCE;
      }
      if (header != NULL) {
        free(header);
      }

      return -1;
    }
  }
  wcsfree (wcs);
  wcs = NULL;
  if (header != NULL) {
    free(header);
  }
  return 0;

}


static int
CheckFullScampUCACResult(PCHECKCOMMON pCheckCommon)
{
  /* First make sure that we are improving matters */
  /* Removed: ((pCheckCommon->oldFitWCS & FITWCS_INACCURATEPV) == 0) */
  if (pCheckCommon->oldPolyRefStars != DEFAULT_POLYREFSTARS) {
    if ((pCheckCommon->oldPolyRefStars > pCheckCommon->polyRefStars) ||
        ((sqr(pCheckCommon->oldPolyRARMS)+sqr(pCheckCommon->oldPolyDecRMS)) <
         (sqr(pCheckCommon->polyRARMS)+sqr(pCheckCommon->polyDecRMS)))) {
      printf("ERRORD: setScampRMS line %d No astrometry improvement  PolyRefStars %d->%d PolyRARMS %.2f->%.2f PolyDecRMS %.2f->%.2f for %s s%d\n",
             __LINE__,
             pCheckCommon->oldPolyRefStars,
             pCheckCommon->polyRefStars,
             pCheckCommon->oldPolyRARMS,
             pCheckCommon->polyRARMS,
             pCheckCommon->oldPolyDecRMS,
             pCheckCommon->polyDecRMS,
             pCheckCommon->platename,
             pCheckCommon->solutionNumber);
      pCheckCommon->newFitWCS &= ~FITWCS_FULLFAILEDSCAMPUCAC;
      pCheckCommon->newFitWCS |= FITWCS_FULLFAILEDTOLERANCE;
      return -1;
    }
  }

  if (CheckScampHeader(pCheckCommon,1) != 0) {
    return -1;
  }

  /* All checks pass */
  pCheckCommon->newFitWCS &= ~FITWCS_FULLFAILEDSCAMPUCAC;
  pCheckCommon->newFitWCS |= FITWCS_FULLSUCCEEDED;
  printf("setScampRMS line %d succeeded PolyRefStars %d->%d PolyRARMS %.2f->%.2f PolyDecRMS %.2f->%.2f for %s s%d\n",
         __LINE__,
         pCheckCommon->oldPolyRefStars,
         pCheckCommon->polyRefStars,
         pCheckCommon->oldPolyRARMS,
         pCheckCommon->polyRARMS,
         pCheckCommon->oldPolyDecRMS,
         pCheckCommon->polyDecRMS,
         pCheckCommon->platename,
         pCheckCommon->solutionNumber);
  return 0;
}


int
main(int argc, char *argv[])
{
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  int nvals;
  CHECKCOMMON checkCommon;
  PCHECKCOMMON pCheckCommon = &checkCommon;
  int function = 0;
  int functionIndex;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  int result;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int rotation;
  int binning;
  int res = 0;
  char queryString[MAX_QUERY_STRING];
  int gotAnswer = 0;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;

  memset(pCheckCommon, 0, sizeof(CHECKCOMMON));

  pCheckCommon->solutionNumber = -1;
  pCheckCommon->polyRefStars = DEFAULT_POLYREFSTARS;
  pCheckCommon->oldPolyRefStars = DEFAULT_POLYREFSTARS;

  // Handle arguments

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    if (argstr[0] != '-') {
      printf("ERROR: setScampRMS line %d  Stray argument %s\n",__LINE__,argstr);
      errorFlag = 1;
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'p': /* plate name */
        case 'P':
          argc--;
          if (argc < 1) {
            printf("ERROR: setScampRMS line %d  Insufficient arguments for -%c\n",__LINE__,cmdchar);
            errorFlag = 1;
          } else {
            strncpy(pCheckCommon->platename,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: setScampRMS line %d  MAX_BUFFER too small for %s\n",__LINE__,*argv);
            }
          }
          break;

        case 'e': /* solution number */
        case 'E':
          argc--;
          if (argc < 1) {
            printf("ERROR: setScampRMS line %d  Insufficient arguments for -%c\n",__LINE__,cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&pCheckCommon->solutionNumber);
            if (nvals != 1) {
              printf("ERROR: setScampRMS line %d  Unable to decode the solutionNumber %s\n",__LINE__,*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'n': /* PolyRefStars */
        case 'N':
          argc--;
          if (argc < 1) {
            printf("ERROR: setScampRMS line %d  Insufficient arguments for -%c\n",__LINE__,cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&pCheckCommon->polyRefStars);
            if (nvals != 1) {
              printf("ERROR: setScampRMS line %d  Unable to decode PolyRefStars %s\n",__LINE__,*argv);
              errorFlag = 1;
            } else {
              function |= ARG_POLYREFSTARS;
            }
          }
          break;

        case 'd': /* PolyDecRMS */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: setScampRMS line %d  Insufficient arguments for -%c\n",__LINE__,cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&pCheckCommon->polyDecRMS);
            if (nvals != 1) {
              printf("ERROR: setScampRMS line %d  Unable to decode PolyDecRMS %s\n",__LINE__,*argv);
              errorFlag = 1;
            } else {
              function |= ARG_POLYDECRMS;
            }
          }
          break;

        case 'r': /* PolyRARMS */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: setScampRMS line %d  Insufficient arguments for -%c\n",__LINE__,cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&pCheckCommon->polyRARMS);
            if (nvals != 1) {
              printf("ERROR: setScampRMS line %d  Unable to decode PolyRARMS %s\n",__LINE__,*argv);
              errorFlag = 1;
            } else {
              function |= ARG_POLYRARMS;
            }
          }
          break;

        case 'c': /* Check success of new procedure */
        case 'C':
          function |= ARG_FULLSUCCESS;
          break;

        case 'i': /* Check success of new procedure with previous solution */
        case 'I':
          function |= ARG_CHECKPREV;
          break;

        case 's': /* Check success of new procedure with previous solution */
        case 'S':
          function |= ARG_CHECKSTALE;
          break;

        case 'f': /* Check failure of new procedure */
        case 'F':
          function |= ARG_FULLFAILURE;
          break;

        case 'a': /* Declare success of astrometry.net on full mosaic */
        case 'A':
          function |= ARG_FULLASTROMETRY;
          break;

        case 'g': /* Check success of astrometry.net on full mosaic */
        case 'G':
          function |= ARG_FULLCHECKASTROM;
          break;

        case 'q': /* Check success of astrometry.net on full mosaic for multiple exposures */
        case 'Q':
          function |= ARG_FULLCHECKASTROMSECOND;
          break;

        case 'b': /* Declare success of SCAMP USNO-B fitting on full mosaic */
        case 'B':
          function |= ARG_FULLSCAMPUSNO;
          break;

        case 'k': /* Check success of SCAMP USNO-B fitting on full mosaic */
        case 'K':
          function |= ARG_FULLCHECKUSNO;
          break;

        case 'y': /* Reset previous failure bits */
        case 'Y':
          function |= ARG_FULLRESETFAIL;
          break;

        case 'h': /* Declare success of SCAMP UCAC4 fitting on the full mosaic */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: setScampRMS line %d  Insufficient arguments for -%c\n",__LINE__,cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) < (MAX_BUFFER-5)) {
              strcpy(pCheckCommon->headername,*argv);
              function |= ARG_FULLSCAMPUCAC;
            } else {
              printf("ERROR: setScampRMS line %d  header file name length %zu is too long\n",__LINE__,strlen(*argv));
              errorFlag = 1;
            }
          }
          break;

        case 'l': /* Declare success of SCAMP UCAC4 fitting for multiple exposures */
        case 'L':
          function |= ARG_FULLSCAMPUCACSECOND;
          break;

        case 'm': /* Check the mosaic headef for InaccuratePV */
        case 'M':
          argc--;
          if (argc < 1) {
            printf("ERROR: setScampRMS line %d  Insufficient arguments for -%c\n",__LINE__,cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) < (MAX_BUFFER-5)) {
              strcpy(pCheckCommon->headername,*argv);
              function |= ARG_CHECKHEADER;
            } else {
              printf("ERROR: setScampRMS line %d  header file name length %zu is too long\n",__LINE__,strlen(*argv));
              errorFlag = 1;
            }
          }
          break;

        default:
          printf("ERROR: setScampRMS line %d   unknown command -%c\n",__LINE__,cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  // Verify that we have everything

  if (pCheckCommon->platename[0] == 0) {
    printf("ERROR: setScampRMS line %d  No platename was specified\n",__LINE__);
    errorFlag = 1;
  }

  if (pCheckCommon->solutionNumber < 0) {
    printf("ERROR: setScampRMS line %d  No solutionNumber specified\n",__LINE__);
    errorFlag = 1;
  }

  for (functionIndex = 0; functionIndex < functionlistsize; functionIndex++) {
    if (function == functionlist[functionIndex]) {
      break;
    }
  }

  if (functionIndex == functionlistsize) {
    printf("ERROR: setScampRMS line %d  illegal/unrecognized combination of function flags\n",__LINE__);
    errorFlag = 1;
  }

  if (isnan(pCheckCommon->polyRARMS)) {
    printf("ERROR: setScampRMS line %d  polyRARMS is nan for %s\n",__LINE__,pCheckCommon->platename);
    pCheckCommon->polyRARMS = 0;
  }

  if (isnan(pCheckCommon->polyDecRMS)) {
    printf("ERROR: setScampRMS line %d  pCheckCommon->polyDecRMS is nan for %s\n",__LINE__,pCheckCommon->platename);
    pCheckCommon->polyDecRMS = 0;
  }

  if (errorFlag) {
    printf("Usage for original SCAMP algorithm: \n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -n <PolyRefStars> -r <PolyRARMS> -d <PolyDecRMS>\n");
    printf("Usage to check new algorithm success:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -c\n");
    printf("Usage to check and set new algorithm failure:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -f\n");
    printf("Usage to declare new algorithm astrometry.wcs success:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -a\n");
    printf("Usage to check readiness to run the SCAMP USNO procedure:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -g\n");
    printf("Usage to declare new algorithm SCAMP USNO success:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -q\n");
    printf("Usage to declare new algorithm SCAMP multiple exposure astrometry.wcs success:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -b\n");
    printf("Usage to check readiness to run SCAMP UCAC procedure:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -k\n");
    printf("Usage to declare new algorithm SCAMP UCAC success providing that the new header is acceptable:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -h <New mosaic header> -n <PolyRefStars> -r <PolyRARMS> -d <PolyDecRMS>\n");
    printf("Usage to declare new algorithm SCAMP UCAC success for multiple exposures\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -l \n");
    printf("Usage to check header output for scamp for inaccuracy:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -m <New mosaic header>\n");
    printf("Usage to reset full processing status completely:\n");
    printf("  setScampRMS -e <solutionNumber> -p <plate name>   -y\n");
    return -1;
  }

  printf(
    "setScampRMS of %s %s plate %s solution %d function %s stars %d RA RMS %f pixels; Dec RMS %f pixels\n",
    __DATE__,
    __TIME__,
    pCheckCommon->platename,
    pCheckCommon->solutionNumber,
    functionNames[functionIndex],
    pCheckCommon->polyRefStars,
    pCheckCommon->polyRARMS,
    pCheckCommon->polyDecRMS
  );

  dasch_init_scandb(pConnection);

  if (ParseFilename(pCheckCommon->platename,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
    printf("ERROR: setScampRMS line %d Failed to parse filename %s as %s%05d_%02d_%02d rotation %d\n",__LINE__,
           pCheckCommon->platename,series,plateNumber,mosaicNumber,binning,rotation);
    return 1;
  }

  // Do whatever we're supposed to do

  switch(function) {
  case FUNCTION_ORIGINAL:
    break;

  case FUNCTION_CHECKSTALE:
    /* Query the database to get the current status of the previous solution */
    if (pCheckCommon->solutionNumber > 0) {
      gotAnswer = 0;
      sprintf(queryString,"SELECT FitWCS+0 from mosaics where series = '%s' AND plateNumber = %d AND mosaicNumber = %d and solutionNumber = %d;",series,plateNumber,mosaicNumber,pCheckCommon->solutionNumber-1);
      res = ExecuteQuery(pConnection,queryString);
      if (!res) {

        res_ptr = mysql_store_result(pConnection);
        if (res_ptr) {
          while ((sqlrow = mysql_fetch_row(res_ptr))) {
            if (gotAnswer) {
              printf("ERROR: setScampRMS line %d  two mosaic entries for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
            }
            gotAnswer = 1;
            if (sqlrow[0]) {
              nvals = sscanf(sqlrow[0],"%d",&pCheckCommon->prevFitWCS);
              if (nvals != 1) {
                printf("ERROR: setScampRMS line %d failed to decode FitWCS %s for %s s%d\n",__LINE__,sqlrow[0],pCheckCommon->platename,pCheckCommon->solutionNumber);
                return 1;
              }
            } else {
              printf("ERROR: setScampRMS line %d null FitWCS for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
              return 1;
            }
          }
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: setScampRMS line %d  mysql_store_result failed in pipelineutils.c for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
        return 1;
      }
      if (gotAnswer == 0)  {
        printf("ERROR: setScampRMS line %d  no mosaics found for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
        return 1;
      }
    }
    // fall through!

  default:
    /* Query the database to get the current status of FitWCS */
    gotAnswer = 0;
    sprintf(queryString,"SELECT FitWCS+0,PolyRefStars,PolyRARMS,PolyDecRMS from mosaics where series = '%s' AND plateNumber = %d AND mosaicNumber = %d and solutionNumber = %d;",series,plateNumber,mosaicNumber,pCheckCommon->solutionNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {

      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          if (gotAnswer) {
            printf("ERROR: setScampRMS line %d  two mosaic entries for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
          }
          gotAnswer = 1;
          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&pCheckCommon->oldFitWCS);
            if (nvals != 1) {
              printf("ERROR: setScampRMS line %d failed to decode FitWCS %s for %s s%d\n",__LINE__,sqlrow[0],pCheckCommon->platename,pCheckCommon->solutionNumber);
              return 1;
            }
            pCheckCommon->newFitWCS = pCheckCommon->oldFitWCS;
          } else {
            printf("ERROR: setScampRMS line %d null FitWCS for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
            return 1;
          }
          if (sqlrow[1]) {
            nvals = sscanf(sqlrow[1],"%d",&pCheckCommon->oldPolyRefStars);
            if (nvals != 1) {
              printf("ERROR: setScampRMS line %d failed to decode PolyRefStars %s for %s s%d\n",__LINE__,sqlrow[1],pCheckCommon->platename,pCheckCommon->solutionNumber);
              return 1;
            }
          } else {
            pCheckCommon->oldPolyRefStars = DEFAULT_POLYREFSTARS;
          }

          if (sqlrow[2]) {
            nvals = sscanf(sqlrow[2],"%lf",&pCheckCommon->oldPolyRARMS);
            if (nvals != 1) {
              printf("ERROR: setScampRMS line %d failed to decode PolyRARMS %s for %s s%d\n",__LINE__,sqlrow[2],pCheckCommon->platename,pCheckCommon->solutionNumber);
              return 1;
            }
          } else {
            pCheckCommon->oldPolyRARMS = 0;
          }
          if (sqlrow[3]) {
            nvals = sscanf(sqlrow[3],"%lf",&pCheckCommon->oldPolyDecRMS);
            if (nvals != 1) {
              printf("ERROR: setScampRMS line %d failed to decode PolyDecRMS %s for %s s%d\n",__LINE__,sqlrow[3],pCheckCommon->platename,pCheckCommon->solutionNumber);
              return 1;
            }
          } else {
            pCheckCommon->oldPolyDecRMS = 0;
          }
        }
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: setScampRMS line %d  mysql_store_result failed in pipelineutils.c for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      return 1;
    }
    if (gotAnswer == 0)  {
      printf("ERROR: setScampRMS line %d  no mosaics found for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      return 1;
    }
  }

  if (pCheckCommon->polyRefStars == DEFAULT_POLYREFSTARS) {
    pCheckCommon->polyRefStars = pCheckCommon->oldPolyRefStars;
    pCheckCommon->polyRARMS = pCheckCommon->oldPolyRARMS;
    pCheckCommon->polyDecRMS = pCheckCommon->oldPolyDecRMS;
  }

  result = -1;

  switch(function) {
  case FUNCTION_CHECKNEW:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERRORD: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    } else {
      result = 0;
    }
    break;

  case FUNCTION_CHECKPREV:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) == 0) {
      printf("WARNING: setScampRMS line %d found wcs has not been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    } else {
      result = 0;
    }
    break;

  case FUNCTION_CHECKSTALE:
    if (pCheckCommon->solutionNumber == 0) {
      result = 0;
    } else {
      if (((pCheckCommon->prevFitWCS & FITWCS_FULLSUCCEEDED) != 0) &&
          ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) == 0)) {
      printf("ERROR: stale setScampRMS line %d found wcs has not been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
      } else {
        result = 0;
      }
    }
    break;

  case FUNCTION_SETFULLFAIL:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_MASK_FULLFAILURE) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been attempted with the full mosaic algorithm FitWCS %d for %s s%d\n",__LINE__,pCheckCommon->oldFitWCS,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    pCheckCommon->newFitWCS |= FITWCS_FULLFAILEDASTROMETRY;
    result = 0;
    break;

  case FUNCTION_SETFULLASTROMETRY:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDASTROMETRY) == 0) {
      printf("ERROR: setScampRMS line %d found wcs astrometry has not been attempted with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    pCheckCommon->newFitWCS &= ~FITWCS_FULLFAILEDASTROMETRY;
    pCheckCommon->newFitWCS |= FITWCS_FULLFAILEDSCAMPUSNO;
    result = 0;
    break;

  case FUNCTION_CHECKFULLASTROM:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDASTROMETRY) != 0) {
      printf("ERRORD: setScampRMS line %d found wcs astrometry failed with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDSCAMPUSNO) == 0) {
      printf("ERRORD: setScampRMS line %d found wcs astrometry not attempted with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    result = 0;
    break;

  case FUNCTION_CHECKFULLASTROMSECOND:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDASTROMETRY) != 0) {
      printf("ERROR: setScampRMS line %d found wcs astrometry failed with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDSCAMPUSNO) == 0) {
      printf("ERROR: setScampRMS line %d found wcs astrometry not attempted with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    pCheckCommon->newFitWCS &= ~FITWCS_FULLFAILEDSCAMPUSNO;
    pCheckCommon->newFitWCS |= FITWCS_FULLFAILEDSCAMPUCAC;
    result = 0;
    break;

  case FUNCTION_SETFULLSCAMPUSNO:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDSCAMPUSNO) == 0) {
      printf("ERROR: setScampRMS line %d found wcs SCAMP USNO-B fit has not been attempted with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    pCheckCommon->newFitWCS &= ~FITWCS_FULLFAILEDSCAMPUSNO;
    pCheckCommon->newFitWCS |= FITWCS_FULLFAILEDSCAMPUCAC;
    result = 0;
    break;

  case FUNCTION_CHECKFULLUSNO:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDSCAMPUSNO) != 0) {
      printf("ERROR: setScampRMS line %d found wcs SCAMP USNO-B fit failed for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDSCAMPUCAC) == 0) {
      printf("ERROR: setScampRMS line %d found wcs SCAMP USNO-B fit not attempted for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    result = 0;
    break;

  case FUNCTION_CHECKHEADER:
    result = CheckScampHeader(pCheckCommon,0);
    break;

  case FUNCTION_SETFULLSCAMPUCAC:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDSCAMPUCAC) == 0) {
      printf("ERROR: setScampRMS line %d found wcs SCAMP UCAC fit has not been attempted with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    result = CheckFullScampUCACResult(pCheckCommon);
    break;

  case FUNCTION_SETFULLSCAMPUCACSECOND:
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLSUCCEEDED) != 0) {
      printf("ERROR: setScampRMS line %d found wcs has already been solved with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    if ((pCheckCommon->oldFitWCS & FITWCS_FULLFAILEDSCAMPUCAC) == 0) {
      printf("ERROR: setScampRMS line %d found wcs SCAMP UCAC fit has not been attempted with the full mosaic algorithm for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
      break;
    }
    pCheckCommon->newFitWCS &= ~FITWCS_FULLFAILEDSCAMPUCAC;
    pCheckCommon->newFitWCS |= FITWCS_FULLSUCCEEDED;
    result = 0;
    break;

  case FUNCTION_RESETFULLFAIL:
    if ((pCheckCommon->oldFitWCS && FITWCS_FULLSUCCEEDED) != 0) {
      /* Set fit parameters to zero rather than null so we know that they have been reset */
      pCheckCommon->polyRefStars = 0;
      pCheckCommon->polyRARMS = 0;
      pCheckCommon->polyDecRMS = 0;
    }

    pCheckCommon->newFitWCS &= ~(FITWCS_MASK_FULLFAILURE);
    result = 0;
    break;

 case FUNCTION_ORIGINAL:
    result = 0;
    break;

 default:
    printf("ERROR: setScampRMS line %d unrecognized function %d for %s s%d\n",__LINE__,function,pCheckCommon->platename,pCheckCommon->solutionNumber);
    break;
  }

  if ((result == 0) && (pCheckCommon->newFitWCS != pCheckCommon->oldFitWCS)) {
    sprintf(queryString,"UPDATE mosaics SET FitWCS = %d where series = '%s' AND plateNumber = %d AND mosaicNumber = %d and solutionNumber = %d;",pCheckCommon->newFitWCS,series,plateNumber,mosaicNumber,pCheckCommon->solutionNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      result = 0;
    } else {
      result = -1;
      printf("ERROR: setScampRMS line %d WCSFit update failed for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
    }
  }

  if ((result == 0) &&
      ((pCheckCommon->polyRefStars != pCheckCommon->oldPolyRefStars) ||
       (pCheckCommon->polyRARMS != pCheckCommon->oldPolyRARMS) ||
       (pCheckCommon->polyDecRMS != pCheckCommon->oldPolyDecRMS))) {
    sprintf(queryString,"UPDATE mosaics SET PolyRefStars = %d, PolyRARMS = %.2f, PolyDecRMS = %.2f where series = '%s' AND plateNumber = %d AND mosaicNumber = %d and solutionNumber = %d;",pCheckCommon->polyRefStars,pCheckCommon->polyRARMS,pCheckCommon->polyDecRMS,series,plateNumber,mosaicNumber,pCheckCommon->solutionNumber);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      result = 0;
    } else {
      printf("ERROR: setScampRMS line %d update failed for %s s%d\n",__LINE__,pCheckCommon->platename,pCheckCommon->solutionNumber);
    }
  }

  mysql_close(pConnection);
  return result;
}

/* Dec 16, 2009 Edward J. Los - Initial version
 * Feb 25, 2013 Edward J. Los - Add support for full sized mosaic fitting with astrometry.net and a two-stage SCAMP procedure,
 *                              by supporting the  FullFailedAstrometry, FullFailedScampUSNO, FullFailedScampUCAC, FullFailedTolerance,
 *                              and FullSucceeded FitWCS keywords
 * Jun 27, 2013 Edward J. Los - Add a check for stale source lists, i.e. FullSucceeded for previous but not for current.
 * Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
 */
