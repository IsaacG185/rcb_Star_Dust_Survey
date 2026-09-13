// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* Get unfinished mosaics and write appropriate commands 
 *
 * Compile with rpm installation
 *
 * gcc -ggdb -O0  -I/usr/include/mysql -I /dasch/install/include  -L /dasch/install/lib getscale.c pipelineutils.a  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -lcfitsio -lwcs -lm  -o getscale
 *
 *  Other inputs: DASCH_USERNAME - MySQL user name
 *                DASCH_PASSWORD - MySQL password
 *
 *  Jan 17, 2008 Edward J. Los - split am series into two.
 *  May 23, 2008 Edward J. Los - remove am split because there appears to be no difference.
 *  May  3, 2010 Edward J. Los - Add multiple "ac" series lenses
 *  May  7, 2010 Edward J. Los - Add additional sub-series
 *  Aug 18, 2020 Edward J. Los - Use capital letters for sub-series to avoid confusion with plate numbers
 *  Nov 13, 2010 Edward J. Los - Count the "InaccuratePV" instances
 *  Nov 29, 2010 Edward J. Los - Compute Poly statistics
 *  Mar 28, 2012 Edward J. Los - Move from the mysql directory
 *  Jun  4, 2012 Edward J. Los - Add DUPLICATEREFCOUNT to INACCURATEPV and rename to InaccurateFit
 *  Oct 29, 2013 Edward J. Los - target.txt to starbase
 *  Jun 20, 2018 Edward J. Los - Add PolyRefStars to plate-specific printout
 *                             - Add SINGLE_SERIES qualifier
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "pipelineutils.h"
#define USE_WCSKINIT 1
#define RMS_FLAG_SCALEA 1 /* Scale exceeds the 3 sigma limit */
#define RMS_FLAG_SCALEB 2
#define RMS_FLAG_ANGLE  4 /* Angle exceeds the 3 sigma limit */
#define RMS_LIMIT 3.0

#define RMS_FLAG_BADFIT 8 /* bad astrometric fit */
#define RMS_FLAG_BADPOLY 16 /* bad polynomial fit rms */
#define POLY_PIXEL_LIMIT 3.0 /* polynomial pixel limit (1.414 * 3.0) */
#define NUM_MOSAICS_LIMIT 20 /* need at least 20 plates in the series to find rms */
#define LARGE_NAXIS 30000 /* Size limit for 14" x 17" plates (13 inches) */
#define MAX_BUFFER 80
/* #define SINGLE_SERIES "ax" */

/* FitWCS bitmask definitions */


struct WorldCoor *
wcskinit (
          int	naxis1,		/* Number of pixels along x-axis */
          int	naxis2,		/* Number of pixels along y-axis */
          char	*ctype1,	/* FITS WCS projection for axis 1 */
          char	*ctype2,	/* FITS WCS projection for axis 2 */
          double crpix1, 
          double crpix2,	/* Reference pixel coordinates */
          double crval1, 
          double crval2,	/* Coordinates at reference pixel in degrees */
          double *cd,		/* Rotation matrix, used if not NULL */
          double cdelt1, 
          double cdelt2,	/* scale in degrees/pixel, ignored if cd is not NULL */
          double crota,		/* Rotation angle in degrees, ignored if cd is not NULL */
          int 	 equinox, /* Equinox of coordinates, 1950 and 2000 supported */
          double epoch);	/* Epoch of coordinates, used for FK4/FK5 conversion
                                 * no effect if 0 */





typedef struct mosaic_entry {
  char series[MAX_SERIES_STRING]; /* Series */
  char origseries[MAX_SERIES_STRING]; /* Series */
  int plateNumber;
  double cd1_1;
  double cd1_2;
  double cd2_1;
  double cd2_2;
  int scanNumber;
  int solutionNumber;
  int naxis1;
  int naxis2;
  char ctype1[MAX_CTYPE_STRING];
  char ctype2[MAX_CTYPE_STRING];
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double scaleA;
  double scaleB;
  double dotangle;
  double PolyRARMS;
  double PolyDecRMS;
  double PolyRMS;
  int rmsFlag;
  int FitWCS;
  int PolyRefStars;
  int inaccuratefit;
  int rotation;
  int mosaicNumber;
} MOSAIC_ENTRY,*PMOSAIC_ENTRY;



typedef struct series_entry {
  char series[MAX_SERIES_STRING]; /* Series */
  int numMosaics; /* Number of mosaics found */
  int numWCSFit;  /* Number of WCS fits compared */
  int numInaccuratefit; /* Number of inaccurate PV entries */
  int numBadFlags; /* Number with bad flags */
  int numPoly; /* Number of Poly plates without inaccurate PV */
  double minScale; /* arcsec per pixel for x16 bin */
  double maxScale; /* arcsec per pixel for x16 bin */
  double aveScale;
  double squaredScale;
  double rmsScale;
  double minScaleRatio; /* horiz vs vert scale */
  double maxScaleRatio; /* horiz vs vert scale */
  double aveAngle;   /* Average absolutes angle, degrees */
  double maxAngle;   /* Maximum angle, degrees */
  double aveRAError; /* Average RA error of plate center */
  double aveDECError; /* Average DEC error of plate center */
  double aveScale2;
  double squaredScale2;
  double rmsScale2;
  double avedotangle;
  double squaredotangle;
  double rmsdotangle;
  double rmsPolyRefStars;
  double rmsPolyRMS; /* The rms of sqrt(sqr(PolyRARMS**2) + sqr(PolyDecRMS)); */
} SERIES_ENTRY, *PSERIES_ENTRY;


SERIES_ENTRY seriesList[MAX_SERIES];

struct WorldCoor *
GetFITSWCS (char *filename, char *header, int verbose,double * cra, double *cdec,double* dra, double *ddec, double* secpix,int* wp,int* hp,
	    int *sysout, double* eqout);

static int verbose = 0;		/* verbose/debugging flag */
static int dss = 0;		/* Flag to drop extra stuff for DSS */
static int dssc = 0;		/* Flag to drop extra stuff for DSS */
static int degout = 0;		/* Flag to print center in degrees */
static double eqout = 0.0;
static double eqim = 0.0;
static int sysout0 = 0;
static int sysim = 0;
static int printepoch = 0;
static int printrange = 0;	/* Flag to print range rather than center */
static int version = 0;		/* If 1, print only program name and version */
static int ndec = 3;		/* Number of decimal places in non-angles */

int SeriesCompare(const void *first, const void *second) 
{
  PSERIES_ENTRY pFirst  = (PSERIES_ENTRY)first;
  PSERIES_ENTRY pSecond = (PSERIES_ENTRY)second;
  return(strcmp(pFirst->series,pSecond->series));
}

int PlateCompare(const void *first, const void *second) 
{
  int result;
  PMOSAIC_ENTRY pFirst  = (PMOSAIC_ENTRY)first;
  PMOSAIC_ENTRY pSecond = (PMOSAIC_ENTRY)second;
  result = strcmp(pFirst->series,pSecond->series);
  if (result != 0) {
    return(result);
  }
  if (pFirst->plateNumber > pSecond->plateNumber) {
    return(1);
  } else if (pFirst->plateNumber < pSecond->plateNumber) {
    return(-1);
  } else {
    return(0);
  }
}


/* 
 * Query the database for an exposure table entry.  Returns 1 if
 * an entry was found
 */



void GetScanComment(MYSQL *pConnection,char *series,int plateNumber,int scanNumber, char* buffer,int bufferLength)
{
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  buffer[0] = 0;
  char queryString[MAX_QUERY_STRING];
  sprintf(queryString,"SELECT scanComment from scans where series = '%s' and plateNumber = %d and scanNumber = %d\n",series,plateNumber,scanNumber);
  res = mysql_query(pConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pConnection);
    if (res_ptr) {
#if 0
      printf("%s: Retrieved %lu rows\n",__FUNCTION__,(unsigned long)mysql_num_rows(res_ptr));
#endif
      sqlrow = mysql_fetch_row(res_ptr);
      if (sqlrow && sqlrow[0]) {
        strncpy(buffer,sqlrow[0],bufferLength);
      } else {
#if 0
        printf("%s: mysql_fetch_row failed %x\n",__FUNCTION__,sqlrow);
#endif
      }
      mysql_free_result(res_ptr);
    } else {
      printf("%s: mysql_store_result failed\n",__FUNCTION__);
    }
  } else {
    printf("%s: mysql_query failed\n",__FUNCTION__);
  }
}
#ifndef USE_WCSKINIT

char * CreateFitsHeader(int naxis1,
                        int naxis2,
                        char *ctype1,
                        char *ctype2,
                        double crval1,
                        double crval2,
                        double crpix1,
                        double crpix2,
                        double cd1_1,
                        double cd1_2,
                        double cd2_1,
                        double cd2_2)
{
  char *fitsheader;
  int nblock;
  int nbfits;
  char endline[81];
  int i;

  nblock = 2;
  nbfits = (nblock + 5) * 2880 + 4;
  fitsheader = (char *) calloc (nbfits,1);
  /* Set up last line of FITS header */
  (void)strncpy (endline,"END", 3);
  for (i = 3; i < 80; i++)
    endline[i] = ' ';
  endline[80] = 0;
  hlength(fitsheader,nbfits);
  strncpy(fitsheader,endline,80);
  hputl (fitsheader, "SIMPLE", 1);
  hputi4 (fitsheader,"BITPIX",16);
  hputi4 (fitsheader,"NAXIS",2);
  hputi4 (fitsheader,"NAXIS1",naxis1);
  hputi4 (fitsheader,"NAXIS2",naxis2);
  hputi4 (fitsheader,"EPOCH", 2000);
  hputi4 (fitsheader,"EQUINOX", 2000);
  hputs  (fitsheader,"CTYPE1",ctype1);
  hputs  (fitsheader,"CTYPE2",ctype2);
  hputnr8(fitsheader,"CRVAL1", 9, crval1);
  hputnr8(fitsheader,"CRVAL2", 9, crval2);
  hputnr8(fitsheader,"CRPIX1", 4, crpix1);
  hputnr8(fitsheader,"CRPIX2", 4, crpix2);
  hputnr8(fitsheader,"CD1_1", 12, cd1_1);
  hputnr8(fitsheader,"CD1_2", 12, cd1_2);
  hputnr8(fitsheader,"CD2_1", 12, cd2_1);
  hputnr8(fitsheader,"CD2_2", 12, cd2_2);
  return(fitsheader);

}

#endif /* USE_WCSKINIT */
int main(int argc,char *argv[])
{
  int res;
  MYSQL my_connection;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  char series[MAX_SERIES_STRING+1];
  char origseries[MAX_SERIES_STRING+1];
  char fullname[MAX_BUFFER];
  int plateNumber;
  double cd1_1;
  double cd2_2;
  double cd2_1;
  double cd1_2;
  double scaleA;
  double scaleB;
  /* double scale; */
  double scaleRatio;
  double angle;
  double angle2;
  double cosdot;
  int scanNumber;
  int nvals;
  char filename[256];
  int numSeries = 0;
  int index;
  PSERIES_ENTRY pSeriesEntry;
  int naxis1;
  int naxis2;
  char *ctype1;
  char *ctype2;
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  int solutionNumber;
  double x;
  double dx;
  double y;
  double dy;
  double cra;
  double cdec;
  double edec;
  double era;
  EXPOSURE exposureTable;
  PEXPOSURE pExposure = &exposureTable;
  FILE *gmtHandle1;
  FILE *gmtHandle2;
  char gmtFilename1[] = "/dasch/scanner/linux/gmt/target.txt";
  char gmtFilename2[] = "/dasch/scanner/linux/gmt/aplates.txt";
  struct WorldCoor *wcs;
  char *header;		/* FITS image header */
  int verbose = 0;
  double dra;
  double ddec;
  double secpix;
  int wp;
  int hp;
  int sysout;
  double eqout;
  char *username;
  char *password;
  int rmsCheckCount = 0;
  int rmsErrorCount = 0;
  double totalAveRAError = 0; /* Average RA error of plate center */
  double totalAveDECError = 0; /* Average DEC error of plate center */
  int totalWCSFit = 0;
  PMOSAIC_ENTRY plateList = NULL;
  PMOSAIC_ENTRY pPlateEntry;
  int plateIndex = 0;
  int totalPlates = 0;
  int seriesIndex;
  int bestMosaicNumber;
  int bestRotation;
  char queryString[MAX_QUERY_STRING];

#ifdef USE_WCSKINIT
  double cd[4];
  double xcenter;
  double ycenter;
#endif /* USE_WCSKINIT */

#ifdef SINGLE_SERIES
  printf("ERROR: SINGLE_SERIES is defined for series %s\n",SINGLE_SERIES);
#endif /* SINGLE_SERIES */

  memset(seriesList,0,sizeof(seriesList));
  gmtHandle1 = fopen(gmtFilename1,"wt");

  if (gmtHandle1 == NULL) {
    printf("Failed to open %s\n",gmtFilename1);
    exit(-1);
  }
	fprintf(gmtHandle1,"dra\tddec\n");
	fprintf(gmtHandle1,"---\t----\n");
  gmtHandle2 = fopen(gmtFilename2,"wt");

  if (gmtHandle2 == NULL) {
    printf("Failed to open %s\n",gmtFilename2);
    exit(-1);
  }

  char scanComment[MAX_COMMENT_STRING];

  for (index = 0; index < MAX_SERIES; index++) {
    pSeriesEntry = &seriesList[index];
    pSeriesEntry->minScale = 5000.0;
    pSeriesEntry->maxScale = 0.0;
    pSeriesEntry->minScaleRatio = 5000.0;
    pSeriesEntry->maxScaleRatio = 0.0;
  }

  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    fprintf(stderr,"DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    fprintf(stderr,"DASCH_PASSWORD is not defined\n");
    return(-1);
  }


  mysql_init(&my_connection);


  if (mysql_real_connect(&my_connection,"localhost.localdomain",username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    printf("Connection success\n");
#ifdef SINGLE_SERIES
    sprintf(queryString,"SELECT series,plateNumber,cd1_1,cd1_2,cd2_1,cd2_2,scanNumber,naxis1,naxis2,ctype1,ctype2,crval1,crval2,crpix1,crpix2,solutionNumber,FitWCS+0,PolyRefStars,PolyRARMS,PolyDecRMS,mosaicNumber,rotation FROM mosaics where series = '%s' and WCSSource = 'imWCS'",SINGLE_SERIES);
    if (strlen(queryString) >= MAX_QUERY_STRING) {
      printf("ERROR: MAX_QUERY_STRING exceeded in line %d\n",__LINE__);
      exit(-1);
    }
    res = mysql_query(&my_connection,queryString);
#else /* SINGLE_SERIES */
    res = mysql_query(&my_connection,"SELECT series,plateNumber,cd1_1,cd1_2,cd2_1,cd2_2,scanNumber,naxis1,naxis2,ctype1,ctype2,crval1,crval2,crpix1,crpix2,solutionNumber,FitWCS+0,PolyRefStars,PolyRARMS,PolyDecRMS,mosaicNumber,rotation FROM mosaics where series != 'mask' and WCSSource = 'imWCS'");
#endif /* SINGLE_SERIES */
    if (!res) {
      printf("insert ID: %lu rows res %d\n",(unsigned long)mysql_affected_rows(&my_connection),res);
      res_ptr = mysql_store_result(&my_connection);
      if (res_ptr) {
        MYSQL_FIELD *field_ptr;
        totalPlates = mysql_num_rows(res_ptr);
        printf("Retrieved %d rows\n",totalPlates);
        plateList = (PMOSAIC_ENTRY)calloc(totalPlates,sizeof(MOSAIC_ENTRY));
        if (plateList == NULL) {
          printf("ERROR: failed to allocate plate list of size %d\n",totalPlates);
          exit(-1);
        }

#if 0
        printf("Column details: \n");
#endif
        while((field_ptr = mysql_fetch_field(res_ptr)) != NULL) {
#if 0
          printf("\t Name: %s\t Type: %d\n",field_ptr->name,field_ptr->type);
#endif
        }
#if 0
        printf("Fetched data:\n");
#endif
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          if (plateIndex >= totalPlates) {
            printf("ERROR: plateIndex %d exceeds totalPlates %d\n",plateIndex,totalPlates);
            exit(-1);
          }
	 
          pPlateEntry = &plateList[plateIndex];
          memset(pPlateEntry,0,sizeof(MOSAIC_ENTRY));
          strncpy(series,sqlrow[0],MAX_SERIES_STRING);
          strncpy(origseries,sqlrow[0],MAX_SERIES_STRING);
          strncpy(pPlateEntry->origseries,sqlrow[0],MAX_SERIES_STRING);
          nvals = sscanf(sqlrow[1],"%d",&plateNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for plateNumber\n");
            continue;
          }
          pPlateEntry->plateNumber = plateNumber;

          if (sqlrow[20] != NULL) {
            nvals = sscanf(sqlrow[20],"%d",&pPlateEntry->mosaicNumber);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for mosaicNumber\n");
              pPlateEntry->mosaicNumber = 0;
              exit(-1);
            }
          } else {
            pPlateEntry->mosaicNumber = 0;
          }
          if (SelectBestMosaic(&my_connection,series,plateNumber,&bestMosaicNumber,&bestRotation) == 1) {
            if (pPlateEntry->mosaicNumber != bestMosaicNumber) {
              continue;
            }

          } else {
            continue;
          }

          nvals = sscanf(sqlrow[6],"%d",&scanNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for scanNumber\n");
            continue;
          }
          pPlateEntry->scanNumber = scanNumber;
          GetScanComment(&my_connection,origseries,plateNumber,scanNumber,scanComment,MAX_COMMENT_STRING);


          if (strcmp(series,"mc") == 0) {
            if  ((plateNumber < 4166) &&
                 (plateNumber >= 3500)) {
              strcat(series,"B"); 
            } else if (plateNumber < 3499) {
              strcat(series,"A");
            } else {
              strcat(series,"C"); /* Vacuum "subsequent to 4171 */
            }
          }
          if (strcmp(series,"ac") == 0) {
            if ((plateNumber >=   1641) && (plateNumber <=   9169)) {
              strcat(series,"B"); /* Cooke #4665 */
            } else if (((plateNumber >=  9170) && (plateNumber <=  9207)) ||
                       ((plateNumber >= 12361) && (plateNumber <= 12365))) {
              strcat(series,"C"); /* Cooke #832 */
            } else if (((plateNumber >=  9208) && (plateNumber <= 12360)) ||
                       ((plateNumber >= 12366) && (plateNumber <= 13050)) ||
                       ((plateNumber >= 13051) && (plateNumber <= 13066))) {
              strcat(series,"D"); /* Cooke #17486 */
            } else if (plateNumber >= 13067) {
              strcat(series,"E"); /* Cooke #V30401 (or #17486 ??) */
            } else {
              strcat(series,"A"); /* Cooke no number */
            }

          }
          if (strcmp(series,"am") == 0) {
            if (plateNumber <= 16722) {
              strcat(series,"A"); /* 1 inch lens */
            } else {
              strcat(series,"B"); /* 1.5 inch lens */
            }
          }
          if (strcmp(series,"mb") == 0) {
            if (plateNumber <= 327) {
              strcat(series,"A"); /* 4 in Cooke lens  191"/mm */
            } else if ((plateNumber >= 328) && (plateNumber <= 1776)) {
              strcat(series,"B"); /* 6 in             354"/mm */
            } else {
              strcat(series,"C"); /* 3 in Ross Lundin 390"/mm */
            }
          }
          strcpy(pPlateEntry->series,series);
          /* NOTE: aco series may have used different lenses */
          /* The "x" series had a reversible crown element for visual or blue focus */

          nvals = sscanf(sqlrow[2],"%lf",&cd1_1);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for cd1_1\n");
            continue;
          }
          pPlateEntry->cd1_1 = cd1_1;
          nvals = sscanf(sqlrow[3],"%lf",&cd1_2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for cd1_2\n");
            continue;
          }
          pPlateEntry->cd1_2 = cd1_2;

          nvals = sscanf(sqlrow[4],"%lf",&cd2_1);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for cd2_1\n");
            continue;
          }
          pPlateEntry->cd2_1 = cd2_1;
          nvals = sscanf(sqlrow[5],"%lf",&cd2_2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for cd2_2\n");
            continue;
          }
          pPlateEntry->cd2_2 = cd2_2;

          nvals = sscanf(sqlrow[7],"%d",&naxis1);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for naxis1\n");
            continue;
          }
          pPlateEntry->naxis1 = naxis1;
          nvals = sscanf(sqlrow[8],"%d",&naxis2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for naxis2\n");
            continue;
          }
          pPlateEntry->naxis2 = naxis2;
          
          if (strcmp(pPlateEntry->series,"a") == 0) {
            if ((pPlateEntry->naxis1 > LARGE_NAXIS) && 
                (pPlateEntry->naxis2 > LARGE_NAXIS)) {
              strcat(series,"A"); 
            } else {
              strcat(series,"B");
            }
            strcpy(pPlateEntry->series,series);
          } else {
            if ((pPlateEntry->naxis1 > LARGE_NAXIS) && 
                (pPlateEntry->naxis2 > LARGE_NAXIS)) {
              printf("ERROR: plate %s%05d is 11x17\n",pPlateEntry->series,pPlateEntry->plateNumber);
            }
          }

        
          ctype1 = sqlrow[9];
          strcpy(pPlateEntry->ctype1,ctype1);
          ctype2 = sqlrow[10];
          strcpy(pPlateEntry->ctype2,ctype2);


          nvals = sscanf(sqlrow[11],"%lf",&crval1);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for crval1\n");
            continue;
          }
          pPlateEntry->crval1 = crval1;
          nvals = sscanf(sqlrow[12],"%lf",&crval2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for crval2\n");
            continue;
          }
          pPlateEntry->crval2 = crval2;

          nvals = sscanf(sqlrow[13],"%lf",&crpix1);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for crpix1\n");
            continue;
          }
          pPlateEntry->crpix1 = crpix1;

          nvals = sscanf(sqlrow[14],"%lf",&crpix2);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for crpix2\n");
            continue;
          }
          pPlateEntry->crpix2 = crpix2;

          nvals = sscanf(sqlrow[15],"%d",&solutionNumber);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for solutionNumber\n");
            continue;
          }
          pPlateEntry->solutionNumber = solutionNumber;

          if (sqlrow[16] != NULL) {
            nvals = sscanf(sqlrow[16],"%d",&pPlateEntry->FitWCS);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for FitWCS\n");
              pPlateEntry->FitWCS = 0;
              exit(-1);
            }
          } else {
            pPlateEntry->FitWCS = 0;
          }



          if (sqlrow[17] != NULL) {
            nvals = sscanf(sqlrow[17],"%d",&pPlateEntry->PolyRefStars);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for PolyRefStars\n");
              pPlateEntry->PolyRefStars = 0;
              exit(-1);
            }
          } else {
            pPlateEntry->PolyRefStars = 0;
          }


          if (sqlrow[18] != NULL) {
            nvals = sscanf(sqlrow[18],"%lf",&pPlateEntry->PolyRARMS);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for PolyRARMS\n");
              pPlateEntry->PolyRARMS = 0;
              exit(-1);
            }
          } else {
            pPlateEntry->PolyRARMS = 0;
          }

          if (sqlrow[19] != NULL) {
            nvals = sscanf(sqlrow[19],"%lf",&pPlateEntry->PolyDecRMS);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for PolyDecRMS\n");
              pPlateEntry->PolyDecRMS = 0;
              exit(-1);
            }
          } else {
            pPlateEntry->PolyDecRMS = 0;
          }

          pPlateEntry->PolyRMS = sqrt((pPlateEntry->PolyRARMS * pPlateEntry->PolyRARMS)+(pPlateEntry->PolyDecRMS * pPlateEntry->PolyDecRMS));
          if ((pPlateEntry->PolyRARMS > POLY_PIXEL_LIMIT) ||
              (pPlateEntry->PolyDecRMS > POLY_PIXEL_LIMIT)) {
            pPlateEntry->rmsFlag |= RMS_FLAG_BADPOLY;
          }
          if (sqlrow[21] != NULL) {
            nvals = sscanf(sqlrow[21],"%d",&pPlateEntry->rotation);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for rotation\n");
              pPlateEntry->rotation = 0;
              exit(-1);
            }
          } else {
            pPlateEntry->rotation = 0;
          }

	  
          plateIndex++;
	

          
  
          for (index = 0; index < numSeries; index++) {
            pSeriesEntry = &seriesList[index];
            if (strcmp(series,pSeriesEntry->series) == 0) {
              break;
            }
          }
          if ((index == numSeries) || (numSeries == 0)) {
            /* New series */
            pSeriesEntry = &seriesList[numSeries];
            strcpy(pSeriesEntry->series,series);
            if (++numSeries >= MAX_SERIES) {
              printf("ERROR: MAX_SERIES is too small\n");
              exit(-1);
            }
          }
          pSeriesEntry->numMosaics++;


          if ((pPlateEntry->FitWCS & (FITWCS_INACCURATEPV|FITWCS_DUPLICATEREFCOUNT)) != 0) {
            pPlateEntry->inaccuratefit = 1;
            pPlateEntry->rmsFlag |= RMS_FLAG_BADFIT;
            pSeriesEntry->numInaccuratefit++;
          } else {
            if ((pPlateEntry->PolyRefStars != 0) &&
                (pPlateEntry->PolyRARMS != 0) &&
                (pPlateEntry->PolyDecRMS != 0)) {
              pSeriesEntry->numPoly++;
              pSeriesEntry->rmsPolyRefStars += 1.0 * pPlateEntry->PolyRefStars * pPlateEntry->PolyRefStars;
              pSeriesEntry->rmsPolyRMS  += (pPlateEntry->PolyRARMS * pPlateEntry->PolyRARMS);
              pSeriesEntry->rmsPolyRMS  += (pPlateEntry->PolyDecRMS * pPlateEntry->PolyDecRMS);
	      
            }


          }




#if 0
          if ((pPlateEntry->plateNumber == 7701) &&
              (strcmp(pPlateEntry->series,"x") == 0)) {
            printf("At %s%05d\n",pPlateEntry->series,pPlateEntry->plateNumber);
          }
#endif

          scaleA = sqrt((cd1_1*cd1_1) + (cd2_1 * cd2_1)) * 3600.0 *16.0;
          if (scaleA < 0.0) {
            scaleA = - scaleA;
          }
          pPlateEntry->scaleA = scaleA;
          pSeriesEntry->aveScale2 += scaleA;
          pSeriesEntry->squaredScale2 += (scaleA*scaleA);
          if (scaleA < pSeriesEntry->minScale) {
            pSeriesEntry->minScale = scaleA;
          }
          if (scaleA > pSeriesEntry->maxScale) {
            pSeriesEntry->maxScale = scaleA;
          }
          scaleB = sqrt((cd1_2 * cd1_2) + (cd2_2 * cd2_2)) * 3600.0 *16.0;
          if (scaleB < 0.0) {
            scaleB = - scaleB;
          }
          pPlateEntry->scaleB = scaleB;
          pSeriesEntry->aveScale2 += scaleB;
          pSeriesEntry->squaredScale2 += (scaleB*scaleB);


          if (scaleB < pSeriesEntry->minScale) {
            pSeriesEntry->minScale = scaleB;
          }
          if (scaleB > pSeriesEntry->maxScale) {
            pSeriesEntry->maxScale = scaleB;
          }
          pSeriesEntry->aveScale += (scaleA+scaleB)/2;
          pSeriesEntry->squaredScale += ((scaleA+scaleB)/2)*((scaleA+scaleB)/2);
          scaleRatio = sqrt((cd1_1 * cd1_1) + (cd2_1 * cd2_1))/sqrt((cd1_2 * cd1_2) + (cd2_2 * cd2_2));
          if (scaleRatio < 0.0) {
            scaleRatio = - scaleRatio;
          }
          if (scaleRatio < pSeriesEntry->minScaleRatio) {
            pSeriesEntry->minScaleRatio = scaleRatio;
          }
          if (scaleRatio > pSeriesEntry->maxScaleRatio) {
            pSeriesEntry->maxScaleRatio = scaleRatio;
          }

          angle = cd1_2/cd1_1;
          if (angle < 0.0) {
            angle = -angle;
          }
          angle2 = cd2_1/cd2_2;
          if (angle2 < 0.0) {
            angle2 = - angle;
          }
          pSeriesEntry->aveAngle = angle + angle2;
          if (angle2 > angle) {
            angle = angle2;
          }
          if (angle > pSeriesEntry->maxAngle) {
            pSeriesEntry->maxAngle = angle;
          }

          angle = atan(angle) * 180.0 / 3.14159;
          cosdot = 3600.0*16.0*3600.0*16.0*((cd1_1*cd2_1)+(cd1_2*cd2_2))/(scaleA*scaleB);
          pPlateEntry->dotangle = acos(cosdot) * 180.0 /3.14159;
          pSeriesEntry->avedotangle += pPlateEntry->dotangle;
          pSeriesEntry->squaredotangle += pPlateEntry->dotangle * pPlateEntry->dotangle;

#if 0
          printf("%s%05d scale %5.3f angle %5.2f cd1_1: %10.7f cd1_2: %13.10f cd2_1: %13.10f cd2_2: %10.7f %s\n",
                 series,plateNumber,scaleRatio,angle,cd1_1,cd1_2,cd2_1,cd2_2,scanComment);
#endif
#ifdef USE_WCSKINIT
          if (strstr(ctype1,"DEC")) {
            char *tmpPtr;
            double dtmp;
            int itmp;
            tmpPtr = ctype2;
            ctype2 = ctype1;
            ctype1 = tmpPtr;
            dtmp = crval1;
            crval1 = crval2;
            crval2 = dtmp;
            dtmp = cd2_1;
            cd2_1 = cd1_1;
            cd1_1 = dtmp;
            dtmp = cd2_2;
            cd2_2 = cd1_2;
            cd1_2 = dtmp;
          }
          cd[0] = cd1_1;
          cd[1] = cd1_2;
          cd[2] = cd2_1;
          cd[3] = cd2_2;


          wcs = wcskinit(naxis1,
                         naxis2,
                         ctype1,
                         ctype2,
                         crpix1,
                         crpix2,
                         crval1,
                         crval2,
                         cd,
                         0,  /* cdelt1 */
                         0,  /* cdelt2 */
                         0,  /* crota */
                         2000, /* equinox */
                         0);   /* epoch */

          xcenter = 0.5 + (0.5 * naxis1);
          ycenter = 0.5 + (0.5 * naxis2);
          pix2wcs(wcs,xcenter,ycenter,&cra,&cdec);

 
#else /* USE_WCSKINIT */


          header = CreateFitsHeader(naxis1,
                                    naxis2,
                                    ctype1,
                                    ctype2,
                                    crval1,
                                    crval2,
                                    crpix1,
                                    crpix2,
                                    cd1_1,
                                    cd1_2,
                                    cd2_1,
                                    cd2_2);
          filename[0] = 0;
          wcs = GetFITSWCS (filename, header, verbose, &cra, &cdec, &dra, &ddec,
                            &secpix, &wp, &hp, &sysim, &eqim);

#endif /* WCSKINIT */
          if (nowcs(wcs)) {
            wcsfree(wcs);
            return(-1);
          }
          wcsfree(wcs);
#ifndef USE_WCSKINIT
          free(header);
#endif /* USE_WCSKINIT */
          if (strcmp(series,"aA") == 0) {
            fprintf(gmtHandle2,"%.5f %.5f %s%05d\n",cra,cdec,origseries,plateNumber);
          }
          if (GetExposureInfo(&my_connection,origseries,plateNumber,0,pExposure,1)) {
            

            edec = str2dec(pExposure->declination);
            era = str2ra(pExposure->rightAscension);
            
            if ((edec == 0) || 
                (era == 0) ||
                (strlen(pExposure->declination) == 0) ||
                (strlen(pExposure->rightAscension )== 0)) {
              printf("Bad exposure info for %4s%05d\n",series,plateNumber);
            } else {
              /*  if (pExposure->timeAccuracy < 0.0010) { */
              double errorDEC;
              double errorRA;
							double radifference;
							radifference = (cra-era);
							if (radifference > 180.0) {
								radifference = 360.0 - radifference;
							}
							if (radifference < -180.0) {
								radifference = -360.0 - radifference;
							}
#if 1
              fprintf(gmtHandle1,"%f\t%f\n",3600.*radifference,3600.*(cdec-edec));
#else
              fprintf(gmtHandle1,"%10.5f %10.5f\n",radifference,cdec-edec);
#endif
              errorDEC = cdec-edec;
              if (errorDEC < 0) {
                errorDEC = -errorDEC;
              }
              errorRA = cra-era;
              if (errorRA < 0) {
                errorRA = - errorRA;
              }
              pSeriesEntry->aveRAError += errorRA;
              pSeriesEntry->aveDECError += errorDEC;
              pSeriesEntry->numWCSFit++;
              totalAveRAError += errorRA;
              totalAveDECError += errorDEC;
              totalWCSFit++;
              /* } */
            }

          }

        }

        mysql_free_result(res_ptr); 
      }
      printf("NOTE: The scale below is in arcseconds/pixel for the bin x 16 plate\n");
      printf("Multiply by 1/16 and divide by 0.011 mm/pixel to get arcsec/mm for an unbinned plate\n");

      /* At this point, sort the series list for readability */
      qsort((void *)seriesList,numSeries,sizeof(SERIES_ENTRY),SeriesCompare);

      /* Also sort the plate list */
      totalPlates = plateIndex;
      qsort((void *)plateList,totalPlates,sizeof(MOSAIC_ENTRY),PlateCompare);

      for (index = 0; index < numSeries; index++) {
        pSeriesEntry = &seriesList[index];
        scaleRatio = pSeriesEntry->maxScale/pSeriesEntry->minScale;
        pSeriesEntry->aveScale =  ( pSeriesEntry->aveScale / pSeriesEntry->numMosaics);
        pSeriesEntry->rmsScale = (pSeriesEntry->squaredScale/pSeriesEntry->numMosaics) - (pSeriesEntry->aveScale * pSeriesEntry->aveScale);
        if (pSeriesEntry->rmsScale > 0) {
          pSeriesEntry->rmsScale = sqrt(pSeriesEntry->rmsScale);
        } else {
          pSeriesEntry->rmsScale = 0;
        }
        pSeriesEntry->aveScale2 =  ( pSeriesEntry->aveScale2 / (2*pSeriesEntry->numMosaics));
        pSeriesEntry->rmsScale2 = (pSeriesEntry->squaredScale2/(2*pSeriesEntry->numMosaics)) - (pSeriesEntry->aveScale2 * pSeriesEntry->aveScale2);
        if (pSeriesEntry->rmsScale2 > 0) {
          pSeriesEntry->rmsScale2 = sqrt(pSeriesEntry->rmsScale2);
        } else {
          pSeriesEntry->rmsScale2 = 0;
        }

        pSeriesEntry->avedotangle =  ( pSeriesEntry->avedotangle / pSeriesEntry->numMosaics);
        pSeriesEntry->rmsdotangle = (pSeriesEntry->squaredotangle/pSeriesEntry->numMosaics) - (pSeriesEntry->avedotangle * pSeriesEntry->avedotangle);
        if (pSeriesEntry->rmsdotangle > 0) {
          pSeriesEntry->rmsdotangle = sqrt(pSeriesEntry->rmsdotangle);
        } else {
          pSeriesEntry->rmsdotangle = 0;
        }


        angle = atan(0.5* (pSeriesEntry->aveAngle/pSeriesEntry->numMosaics)) * 180.0/3.14159;
        angle2 = atan(pSeriesEntry->maxAngle) * 180.0/3.14159;

        printf("%2d:  Series %4s, entries %4d scale %6.2f %6.2f %6.2f %6.4f  ratio %7.5f %7.5f angle %5.2f %5.2f scale2 %6.2f rms %8.4f dot %6.2f rms %8.4f\n",
               index,
               pSeriesEntry->series,
               pSeriesEntry->numMosaics,
               pSeriesEntry->minScale,
               pSeriesEntry->aveScale,
               pSeriesEntry->maxScale,
               scaleRatio,
               pSeriesEntry->minScaleRatio,
               pSeriesEntry->maxScaleRatio,
               angle,angle2,
               pSeriesEntry->aveScale2,
               pSeriesEntry->rmsScale2,
               pSeriesEntry->avedotangle,
               pSeriesEntry->rmsdotangle);
      }




      for (index = 0; index < numSeries; index++) {
        pSeriesEntry = &seriesList[index];

        printf("Series %-3s, entries %4d unbinned arcsec/mm %8.2f rms %8.2f\n",
               pSeriesEntry->series,
               pSeriesEntry->numMosaics,
               pSeriesEntry->aveScale/(16.0 * 0.011),
               pSeriesEntry->rmsScale/(16.0 * 0.011));
      }


      for (index = 0; index < numSeries; index++) {
        pSeriesEntry = &seriesList[index];
        if (pSeriesEntry->numWCSFit == 0) {
          continue;
        }
        pSeriesEntry->aveRAError = pSeriesEntry->aveRAError/pSeriesEntry->numWCSFit;
        pSeriesEntry->aveDECError = pSeriesEntry->aveDECError/pSeriesEntry->numWCSFit;

        printf("%2d:  Series %4s, entries %4d ave RA Error %9.4f ave DEC error %9.4f\n",
               index,
               pSeriesEntry->series,
               pSeriesEntry->numWCSFit,
               pSeriesEntry->aveRAError,
               pSeriesEntry->aveDECError);
      }

    } else {
      printf("Insert error %d: %s res %d\n",mysql_errno(&my_connection),mysql_error(&my_connection),res);
    }



    mysql_close(&my_connection);
  } else {
    printf("Connection failed\n");
    if (mysql_errno(&my_connection)) {
      printf("Connection error %d: %s\n",mysql_errno(&my_connection),mysql_error(&my_connection));
    }
  }
  seriesIndex = 0;
  pSeriesEntry = &seriesList[0];
  for (plateIndex = 0; plateIndex < totalPlates; plateIndex++) {
    double minrmsScale;
    double maxrmsScale;
    double minrmsAngle;
    double maxrmsAngle;
    double stdScaleA;
    double stdScaleB;
    double stdDotAngle;
    pPlateEntry = &plateList[plateIndex];
#if 0
    if ((pPlateEntry->plateNumber == 7701) &&
        (strcmp(pPlateEntry->series,"x") == 0)) {
      printf("At %s%05d\n",pPlateEntry->series,pPlateEntry->plateNumber);
    }
#endif


    while (strcmp(pPlateEntry->series,pSeriesEntry->series) != 0) {
      pSeriesEntry++;
      seriesIndex++;
      if (seriesIndex >= numSeries) {
        printf("ERROR: seriesIndex %d exceeds numSeries %d\n",seriesIndex,numSeries);
        exit(-1);
      }
    }
    if (pSeriesEntry->rmsScale2 > 0) {
      stdScaleA = fabs((pPlateEntry->scaleA - pSeriesEntry->aveScale2)/pSeriesEntry->rmsScale2);
      stdScaleB = fabs((pPlateEntry->scaleB - pSeriesEntry->aveScale2)/pSeriesEntry->rmsScale2);
    } else {
      stdScaleA = 0.0;
      stdScaleB = 0.0;
    }

    if (pSeriesEntry->rmsdotangle > 0) {
      stdDotAngle = fabs((pPlateEntry->dotangle - pSeriesEntry->avedotangle)/pSeriesEntry->rmsdotangle);
    } else {
      stdDotAngle = 0.0;
    }


    minrmsScale = pSeriesEntry->aveScale2 - (RMS_LIMIT * pSeriesEntry->rmsScale2);
    maxrmsScale = pSeriesEntry->aveScale2 + (RMS_LIMIT * pSeriesEntry->rmsScale2);
    minrmsAngle = pSeriesEntry->avedotangle - (RMS_LIMIT * pSeriesEntry->rmsdotangle);
    maxrmsAngle = pSeriesEntry->avedotangle + (RMS_LIMIT * pSeriesEntry->rmsdotangle);
    if ((pPlateEntry->scaleA < minrmsScale) ||
        (pPlateEntry->scaleA > maxrmsScale)) {
      pPlateEntry->rmsFlag |= RMS_FLAG_SCALEA;
    }
    if ((pPlateEntry->scaleB < minrmsScale) ||
        (pPlateEntry->scaleB > maxrmsScale)) {
      pPlateEntry->rmsFlag |= RMS_FLAG_SCALEB;
    }
    if ((pPlateEntry->dotangle < minrmsAngle) ||
        (pPlateEntry->dotangle > maxrmsAngle)) {
      pPlateEntry->rmsFlag |= RMS_FLAG_ANGLE;
    }
    if (pPlateEntry->rmsFlag != 0) {
      pSeriesEntry->numBadFlags ++;
    }
    if ((pPlateEntry->solutionNumber == 0) && 
        (pSeriesEntry->numMosaics >= NUM_MOSAICS_LIMIT)) {

      if (pPlateEntry->rotation == 0) {
        sprintf(fullname,"%s%05d_%02d_01ww",
                pPlateEntry->origseries,
                pPlateEntry->plateNumber,
                pPlateEntry->mosaicNumber);
      } else {
        sprintf(fullname,"%s%05d_%02d_01r%dww",
                pPlateEntry->origseries,
                pPlateEntry->plateNumber,
                pPlateEntry->mosaicNumber,
                pPlateEntry->rotation);

      }

      printf("%s%05d %d scale %5.3f %5.3f dot product %5.3f Ratio %4.2f %4.2f %4.2f flag %2d fit %d PolyRefStars %5d PolyRMS %4.1f %s \n",
             pPlateEntry->series,
             pPlateEntry->plateNumber,
             pPlateEntry->solutionNumber,
             pPlateEntry->scaleA,
             pPlateEntry->scaleB,
             pPlateEntry->dotangle,
             stdScaleA,
             stdScaleB,
             stdDotAngle,
             pPlateEntry->rmsFlag,
             pPlateEntry->inaccuratefit,
             pPlateEntry->PolyRefStars,
             pPlateEntry->PolyRMS,
             fullname);
      if (pPlateEntry->rmsFlag != 0) {
        rmsErrorCount++;
      }
      rmsCheckCount++;
    }

  }
 
  for (index = 0; index < numSeries; index++) {
    pSeriesEntry = &seriesList[index];
    if (pSeriesEntry->numPoly > 0) {
      pSeriesEntry->rmsPolyRefStars = sqrt(pSeriesEntry->rmsPolyRefStars/pSeriesEntry->numPoly);
      pSeriesEntry->rmsPolyRMS = sqrt(pSeriesEntry->rmsPolyRMS/pSeriesEntry->numPoly);
    } else {
      pSeriesEntry->rmsPolyRefStars = 0;
      pSeriesEntry->rmsPolyRMS = 0;	  
    }

    if (pSeriesEntry->numMosaics > 0) {
      printf("Series %4s, entries %5d, Inaccuratefit %5d, Ratio %5.1f percent,bad flags %5d, Ratio %5.1f percent, Scale %8.2f unbinned arcsec/mm numPoly %5d rmsPolyRefStars %6.1f rmsPolyRMS %4.2f\n",
             pSeriesEntry->series,
             pSeriesEntry->numMosaics,
             pSeriesEntry->numInaccuratefit,
             100.0 * pSeriesEntry->numInaccuratefit/(1.0*pSeriesEntry->numMosaics),
             pSeriesEntry->numBadFlags,
             100.0 * pSeriesEntry->numBadFlags/(1.0*pSeriesEntry->numMosaics),
             pSeriesEntry->aveScale/(16.0 * 0.011),
             pSeriesEntry->numPoly,
             pSeriesEntry->rmsPolyRefStars,
             pSeriesEntry->rmsPolyRMS);
		 
    }
  }

  
  printf("Total plates %5d, Wcs fit plates %5d, Checked Plates %5d Error Plates %5d totalAveRAError (deg) %9.4f totalAveDECError %9.4f (deg)\n",totalPlates,totalWCSFit,rmsCheckCount,rmsErrorCount,totalAveRAError/totalWCSFit,totalAveDECError/totalWCSFit);
  if (plateList != NULL) {
    free(plateList);
  }

  fclose(gmtHandle1);
  fclose(gmtHandle2);

  return(EXIT_SUCCESS);
}
