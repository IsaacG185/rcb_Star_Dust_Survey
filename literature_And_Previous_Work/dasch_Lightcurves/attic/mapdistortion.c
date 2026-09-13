// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* mapdistortion.c
 *
 *   Experimental non-production program to find image distortion. 
 *   All of the input parameters are currently hard-coded in.
 *   Preliminary results for the "ax" series show that the distortion is is insignificant, but the TNX fitting in imwsky2k and imwty2 "-n 8" is too unconstrained
 *
 *  gcc -ggdb -O0  -I/usr/include/mysql -I /dasch/install/include  -L /dasch/install/lib mapdistortion.c pipelineutils.a  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -lwcs -lm -o mapdistortion
 *
 * Calculations are borrowed from Doug Mink's imwcs.
 *
 *  Other inputs: DASCH_USERNAME - MySQL user name
 *                DASCH_PASSWORD - MySQL password
 *
 * Jun 15, 2018 Edward J. Los - Initial version
 * Jul  2, 2018 Edward J. Los - Add DEMO_PLOT to show difference between fitting with different TAN reference pixels
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mysql.h"
#include "fitsio.h"
#include "longnam.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include <sys/stat.h>
#define MAX_COMMENT_STRING 132
#define MAX_CENTERSOURCE_STRING 10
#define MAX_RIGHTASCENSION_STRING 13
#define MAX_DECLINATION_STRING 16
#define MAX_DATE_STRING 25
#define MAX_RMS 60 /* pixels */
#define MAX_FILENAME_LEN 256
#define MAX_BUFFER 500
#if 1
#define MAX_OFFSET 0.001
#else
#define MAX_OFFSET 0.05 /* Five percent of the plate width */
#endif
#define A1_FACTOR 0.23666  /* From annular9.m */
#define DEFAULT_SEARCH_PIXELS 40
/* #define PLATE_LIMIT 30  */
#define MAX_PIXDIST  15000 /* sqrt(4**2 + 5**2)*24.5 mm/inch/0.011 mm/pixel = 14261 pixels */
#define SUMMARY_PLOT 1
#define MAX_VECTOR_COUNT 6000
/* #define DEMO_PLOT 1 */

char *
GetFITShead (char *filename,	/* FITS or IRAF file filename */
             int verbose);	/* Print error messages if nonzero */



void ra2str (

             char	*string,	/* Character string (returned) */
             int	lstr,		/* Maximum number of characters in string */
             double	ra,		/* Right ascension in degrees */
             int	ndec);		/* Number of decimal places in seconds */

void
dec2str (

         char	*string,	/* Character string (returned) */
         int	lstr,		/* Maximum number of characters in string */
         double	dec,		/* Declination in degrees */
         int	ndec);		/* Number of decimal places in arcseconds */


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



struct WorldCoor *
GetFITSWCS (char *filename, char *header, int verbose,double * cra, double *cdec,double* dra, double *ddec, double* secpix,int* wp,int* hp,
	    int *sysout, double* eqout);
void
pix2wcs (

         struct WorldCoor *wcs,		/* World coordinate system structure */
         double	xpix, double ypix,	/* x and y image coordinates in pixels */
         double	*xpos,double *ypos);	/* RA and Dec in degrees (returned) */
double wcsdist(
               double x1,
               double y1,	/* (RA,Dec) or (Long,Lat) in degrees */
               double x2,
               double y2);	/* (RA,Dec) or (Long,Lat) in degrees */


void
wcs2pix (
struct WorldCoor *wcs,	/* World coordinate system structure */
double	xpos,double ypos,	/* World coordinates in degrees */
double	*xpix,double *ypix,	/* Image coordinates in pixels */
int	*offscl);	/* 0 if within bounds, else off scale */


void
num2str (

char	*string,	/* Character string (returned) */
double	num,		/* Number */
int	field,		/* Number of characters in output field (0=any) */
int	ndec);		/* Number of decimal places in degree string */




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

typedef struct _plate {
  char Plate[MAX_PLATE_NAME];
  char platename[MAX_PLATE_NAME];
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int flag;
  int fit; 
  int PolyRefStars;
  double PolyRMS;
  int seriesId;
  int rejectFlag[TOTAL_DMAGBINS_NORMAL];
  int drad_bin_size2[TOTAL_DMAGBINS_NORMAL];
#if 0
  int solutionNumber;
  int errorFlag;
  double geoJulianDate;
  EXPOSURE exposure;
  PMOSAIC pMosaic;
  struct WorldCoor *wcs;
  int nx;                   /* dmagcor_table width */
  int ny;                   /* dmagcor_table height */
  PLOCALBIN local_bin_table;
  int starDetected;     /* If non-zero, this object has been detected */
  int validPoint;     /* If non-zero, this point is valid */
  /* Information copied for each point */
  double magcal_magdep; /* Local magnitude, or local limiting magnitude */
  double Date;  /* Heliocentric Julian Date */
  double magcal_local_rms;
  double magcal_iso;
  double magcal_iso_rms;
  double limiting_mag_local;
  double Stdmag;
  int spatial_bin;
  int local_bin_index;
  int NUMBER;
  int AFLAGS;
  int AFLAGSCOPY;
  int BFLAGS;
  /* Information calculated for output */
  double dradRMS2;    /* drad RMS for this bin */
  double FWHM_IMAGE;
  int npoints_local;  /* npout */
  int rejectFlag;     
  double extinction;
  double magcor_local;       /* zout */
  double magcal_local_error; /* errout */


  double clipmed;
  int reject_reason1;  /* Total  Reject reason mask.  See pipelineutils.h for defintions */
  int reject_reason2;  /* Actual Reject reason mask.  See pipelineutils.h for defintions */
  PPHOTTARGET pTarget;
  int gsc_bin_index;
  int coverage_table_refcount;
  int coverage_table_select;
  int leftMargin;
  int rightMargin;
  int bottomMargin;
  int topMargin;
  int versionId;
  /* The following fields are used by the fullFlag option */
  double BACKGROUND;
  double Blendedmag;
  double dec;                 /* Declination in degrees */
  double ELLIPTICITY;
  double FLUX_ISO;
  double FLUX_MAX;
  double FWHM_WORLD;
  double ISOAREA_WORLD;
  double KRON_RADIUS;
  double MAG_APER;
  double MAG_AUTO;
  double MAG_ISO;
  double plate_dist;
  double ra;                  /* Right Ascension in degrees */
  double THETA_J2000;
  double X_IMAGE;
  double Y_IMAGE;
  double dra;
  double ddec;
  int exposureNumber;
  int ISO0;
  int ISO1;
  int ISO2;
  int ISO3;
  int ISO4;
  int ISO5;
  int ISO6;
  int ISO7;
  long long REFNumber;        /* Translated reference number */
 
#endif
} PLATE,*PPLATE;

typedef struct _MAPCOMMON 
{
	FILE *debugFile;
	int idsw;
	int igmprd;
  int irotat;
  int iresol;
} MAPCOMMON, *PMAPCOMMON;

typedef struct _pixdist {
  int count;
  double value;
  double xOffset;
  double yOffset;
  double median;
  double rms;
  int vectorCount;
  double *vector;
} PIXDIST,*PPIXDIST;


/* Print an error message */
void plterr(int errnum,double fxx, double fyy)
{
	printf("ERROR: plterr number %d, fxx = %g, fyy = %g\n",errnum,fxx,fyy);
	return;
}

/* 
 * resolv - converts from cartesian to polar
 *
 *  x - length along x axis
 *  y - length along y axis
 *  z - length along z axis
 *  dec - angle (degrees) above the x-y plane
 *  rha - angle (degrees) in the x-y plane from the x axis
 */
void
resolv(PMAPCOMMON pCommon,
	   double x,
	   double y,
	   double z,
	   double * dec,
	   double * rha)
{
	double temp = (x*x)+(y*y);
	if (temp == 0)
	{
		*dec = 90.;
		if (z < 0)
			*dec = -90.;
	}
	else
	{
		*dec = atan(z/sqrt(temp))*RAD_TO_DEGREES;
	}
	if (x != 0)
	{
		*rha = atan(y/x)*RAD_TO_DEGREES;
		if ((y < 0) && (x < 0))
			*rha = *rha - 180.;
		if ((y >= 0) && (x < 0))
			*rha = *rha + 180.;
	}
	else
	{
		*rha = 90.;
		if (y < 0.0)
			*rha = - 90.;
	}
	if ((pCommon->idsw != 0) && (pCommon->iresol <= 10))
	{
		pCommon->iresol++;
		fprintf(pCommon->debugFile," in reolv: iresol %d, x %g y %g z %g dec %g rha %g\n",
			pCommon->iresol,x,y,z,*dec,*rha);
	}
}

/* create a rotation matrix 
	ind = axis index  1 = about x axis
                    2 = about y axis
					          3 = about z axis
    rth = rotation angle (degrees)
          (use the right hand rule for positive direction 
    gmat = rotation matrix.  Nine elements to convert x,y,z to xp, yp, zp
*/

void
rotate(PMAPCOMMON pCommon,
       int ind,		
       double rth,
       double *gmat)
{
	double pi = 3.141592654;
	double ttmp = (rth*pi)/180.;
	double ci = cos(ttmp);
	double rt;
	double si;
	int i;
	for (i = 0; i < 9; i++) {
		gmat[i] = 0.0;
	}
	for (i = 1; i <= 3; i++) {
		rt = ci;
		if (ind == i)
			rt = 1.;
		gmat[(4*i) - 4] = rt;
	}
	si = sin(ttmp);
	if (ind == 1) {
		gmat[7] = si;
		gmat[5] = - si;
	} else if (ind == 2)	{
		gmat[2] = si;
		gmat[6] = - si;
	} else if (ind == 3) {
		gmat[3] = si;
		gmat[1] = -si;
	} else {
		plterr(50,0.0,0.0);
	}
	if (pCommon->idsw != 0) {
		pCommon->irotat += 1;
		if (pCommon->irotat <= 10) {
			fprintf(pCommon->debugFile," in rotate: irotat %d, ind %d, rth %g \n",
              pCommon->irotat,ind,rth);
			for (i = 0; i < 9; i+=3) {
				fprintf(pCommon->debugFile," gmat %10.5f %10.5f %10.5f\n",
                gmat[i+0],gmat[i+1],gmat[i+2]);
			}
		}
	}
}


/* gmprd, a canned matrix product routine
 * 
 * a = name of the first input matrix
 * b = name of the second input matrix
 * c = name of the output matrix
 * n = number of rows in a
 * m = number of columns in a and rows  in b
 * l = number of columns in b
 *
 * Therefore:
 *  
 * [R] = [A][B] where [R] is l x n
 *                    [A] is n x m
 *                    [B] is m x l
 */
void
gmprd(PMAPCOMMON pCommon,
      double *a,
      double *b,
      double *r,
      int n,
      int m,
      int l)
{
	int ir = 0;
	int ik = -m;
	int k;
	int ji;
	int ib;
	int i;
	int j;


	for (k = 1; k <= l; k++) {
		ik = ik + m;
		for (j = 1; j <= n; j++) {
			ir = ir + 1;
			ji = j - n;
			ib = ik;
			r[ir-1] = 0.;
			for (i = 1; i <= m; i++) {
				ji = ji + n;
				ib = ib + 1;
				r[ir-1] = r[ir-1] + (a[ji-1]*b[ib-1]);
			}


		}
	}
	if (pCommon->idsw != 0) {
		pCommon->igmprd += 1;
		if (pCommon->igmprd <= 10) {
			fprintf(pCommon->debugFile," in gmprod igmprd %d, n %d, m %d, l %d \n",
              pCommon->igmprd,n,m,l);
			for (k = 1; k <= 9; k++) {
				fprintf(pCommon->debugFile,"a[%d] %10.5f b[%d] %10.5f r[%d] %10.5f\n",
                k,a[k-1],k,b[k-1],k,r[k-1]);
			}
		}
	}
}
/* Reverse sort by PolyRefStars */
int PlateCompare(const void *first, const void *second) 
{
  PPLATE pPlateRecFirst = (PPLATE)first;
  PPLATE pPlateRecSecond = (PPLATE)second;
  
  if (pPlateRecFirst->PolyRefStars > pPlateRecSecond->PolyRefStars) {
    return(-1);
  } else if (pPlateRecFirst->PolyRefStars < pPlateRecSecond->PolyRefStars) {
    return(1);
  } else {
  }
}



int main(int argc,char *argv[])
{
  int res;
  MYSQL my_connection;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  double cd1_1;
  double cd2_2;
  double cd2_1;
  double cd1_2;
  double scale = 0;
  double scaleRatio;
  double angle;
  double angle2;
  int expRotation;
  int nvals;
  char filename[MAX_FILENAME_LEN];
  char outname[MAX_FILENAME_LEN];
  char scalename[MAX_FILENAME_LEN];
  FILE *outhandle = NULL;
  FILE *scalehandle = NULL;
  char queryString[MAX_QUERY_STRING];
  int numSeries = 0;
  int index;
  int naxis1;
  int naxis2;
  char *ctype1;
  char *ctype2;
  double crval1;
  double crval2;
  double crpix1;
  double crpix2;
  double x;
  double dx;
  double y;
  double dy;
  double cra;
  double cdec;
  double lontemp;
  double lattemp;
  double era;
  double edec;
  double dist;
  double pixdist;
  double sra[4];
  double sdec[4];
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double xctr;
  double yctr;
  struct WorldCoor *wcs;
  struct WorldCoor *tan_wcs;
  char *header;		/* FITS image header */
  char *image;	      /* Image pixels */
  
  double dra;
  double ddec;
  double secpix;
  int wp;
  int hp;
  int sysout;
  double eqout;
  int mosaicCount = 0;
  int noWCSCount = 0;
  int totalCount = 0;
  int zeroCrossingCount = 0;
  char *mysqlhost;
  char *username;
  char *password;
  char nullPassword = 0;
  fitsfile *fptr = NULL;
  int status = 0;
  int iomode = READWRITE;
  char err_text[FLEN_ERRMSG];
  int foundEntry = 0;
  int foundEntryScans = 0;
  int foundEntryPattern = 0;
  int watIndex;
  char keywordstr[FLEN_KEYWORD];
  double cd[4];
  int errorFlag = 0;
  char *argstr;
  char cmdchar;
  int verboseFlag = 0;
  int gsc23Flag = 0;
  int dradFlag = 0;
	int patrolFlag = 0;
  int ignoreRotation = 0;
  char rstr[32], dstr[32];
  char radstr[32];
  int axisFlag = 0;
  int plateSizeFlag = 0;
  double dRightAscension;
  double dDeclination;
  int scanNumber;
  int transform = 0;
  int WCSSource;
  int pixelRadius = DEFAULT_SEARCH_PIXELS;
  int patternID;
  int leftMargin;
  int rightMargin;
  int topMargin;
  int bottomMargin;
  int leftMargin2;
  int rightMargin2;
  int topMargin2;
  int bottomMargin2;
  int nullMessage = 0;
  int solutionNumber = 0;
  int diskLocation;
  double extractArcsec = 0; /* Extract radius for getfits */
  double extractPixels; /* Pixel width for getfits */
  double galacticout = 0;
  MAPCOMMON map_common;
  PMAPCOMMON pCommon = &map_common;
  int nx = X_DMAGBINS_NORMAL;
  int ny = Y_DMAGBINS_NORMAL;
  int ix;
  int iy;
  int local_bin;
  double vertang;
  double pixDist;
  double radius;
  double drad;
	double vecA[9];
	double vecB[9];
  double vecC[9];
  double vecD[9];
	double matA[9];
  double matB[9];
  double matC[9];
  double cosx;
  double rainv;
  double decinv;
  double denominator;
  double cra2;
  double cdec2;
  double decinv2;
  double rainv2;
  double decinv3;
  double rainv3;
  double verticalang;
  double xPixelsPerBin;
  double yPixelsPerBin;
#define CENTER_POINT TOTAL_DMAGBINS_NORMAL /* The center point if the mosaic */
#define TEST_POINT (TOTAL_DMAGBINS_NORMAL+1)   /* One pixel above the center point along NAXIS2 */
  double xloc[TOTAL_DMAGBINS_NORMAL+2];
  double yloc[TOTAL_DMAGBINS_NORMAL+2];
  double actualra[TOTAL_DMAGBINS_NORMAL+2];
  double actualdec[TOTAL_DMAGBINS_NORMAL+2];
  double mappedra[TOTAL_DMAGBINS_NORMAL+2];
  double mappeddec[TOTAL_DMAGBINS_NORMAL+2];
  double val1[TOTAL_DMAGBINS_NORMAL+2][9];
  double val2[TOTAL_DMAGBINS_NORMAL+2][9];
  double val3[TOTAL_DMAGBINS_NORMAL+2][9];
  int lineLen;
  int lineCounter = 0;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
  char *charPtr;
  char *charPtr2;
  int maxString = 0;
  PPLATE plateTable = NULL;
  PPLATE pPlate = NULL;
  int maxPlates = 0;
  int plateCount = 0;
  int platesProcessed = 0;
  int plateIndex;
  char *headersDirectory;
  struct stat statbuf;
  int statResult;
  time_t curTime;
  time_t startTime;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  int local_bin_index;
  int numRows;
  int curRow = 0;
  int goodLocalBins = 0;
  int pixIndex;
  PPIXDIST pixDistTable = NULL;
  PPIXDIST pixDistSummaryTable = NULL;
  PPIXDIST pPixDist;
  PPIXDIST pPixDistSummary;
  double xOffset;
  double yOffset;
  double fullOffset;
  int maxPixDistCount = 0;
  int maxPixDistSummaryCount = 0;
  int pixDistVectorAllocCount = 0;
  int pixDistVectorError = 0;
  int pixDistCountError = 0;
  int vectorIndex;
  int bestCount;

#ifdef DEMO_PLOT
  int demoCounter = 0;
  int demoIteration = 0; /* Set to 100 for changes with scale with respect to reference pixel */
  double xplatectr;
  double yplatectr;
  double ra1;
  double dec1;
  double ra2;
  double dec2;
  double ra3;
  double dec3;
  double ra4;
  double dec4;
  double scale1;
  double scale2;
#endif /* DEMO_PLOT */
#ifdef PLATE_LIMIT
  printf("ERROR: plateCount limited to %d\n",PLATE_LIMIT);
#endif /* PLATE_LIMIT */

#ifdef DEMO_PLOT
  printf("ERROR: DEMO_PLOT is enabled\n");
#endif /* DEMO_PLOT */
#ifndef MAX_OFFSET
  printf("ERROR: MAX_OFFSET is not defined\n");
#endif /* MAX_OFFSET */

  outname[0] = 0;
  scalename[0] = 0;
  time(&startTime);

  memset(pCommon,0,sizeof(MAPCOMMON));
  pCommon->idsw = 0;
  pCommon->debugFile = stdout;

  headersDirectory = getenv("DASCH_HEADERS");
  if (headersDirectory == NULL) {
    printf("ERROR: DASCH_HEADERS is not defined\n");
    errorFlag = 1;
  }

  pixDistTable = (PPIXDIST)calloc(MAX_PIXDIST,sizeof(PIXDIST));
  pixDistSummaryTable = (PPIXDIST)calloc(MAX_PIXDIST,sizeof(PIXDIST));
  if ((pixDistTable == NULL) || 
      (pixDistSummaryTable == NULL)) {
    printf("ERROR: failed to allocate pixDist tables %lld %lld\n",pixDistTable,pixDistSummaryTable);
    exit(-1);
  }

#if 0
  /* Loop through the arguments */
  filename[0] = 0;
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      if (strlen(filename) > 0) {
        errorFlag = 1;
        printf("ERROR: list %s is being overwritten by %s\n",filename,argstr);
      } else {
        strncpy(filename,argstr,MAX_FILENAME_LEN-1);
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'a': /* show only the axes and the plate scale*/
        case 'A':
          axisFlag = 1;
          break;

        case 'b': /* show the size of the plate */
        case 'B':
          plateSizeFlag = 1;
          break;
       
        case 'g': /* GSC2.3 flag */
        case 'G':
          gsc23Flag = 1;
          dradFlag = 0;
          degout = 0;
          break;
        case 'l': /* Generate galactic coordinates */
        case 'L':
          galacticout = 1;
          degout++;
          break;
        case 'p': /* ignore rotation */
        case 'P':
          ignoreRotation = 1;
          break;
        case 's':  /* Series information */
        case 'S':
          gsc23Flag = 0;
          dradFlag = 1;
          degout = 1;
          break;
        case 'v': /* verbose flag */
        case 'V':
          verboseFlag = 1;
          break;
        case 'n': /* patrol flag */
        case 'N':
          patrolFlag = 1;
          break;
        case 'd': /* degrees out */
        case 'D':
          degout++;
          break;

        case 'r': /* search radius in pixels */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&pixelRadius);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the search radius %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 't': /* getfits radius in arcsec */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&extractArcsec);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the getfits extraction radius %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'e': /* solution number */
        case 'E':
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


        default:
          printf("* illegal command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  if (filename[0] == 0) {
    printf("ERROR: No filename was specified\n");
    errorFlag = 1;
  }
  /* We do not need the plate center if the "-a" or the "-s" or "-b" qualifier is specified */
  if (((axisFlag != 0) || (dradFlag != 0) || (axisFlag)) && 
      (solutionNumber != 0)) {
    solutionNumber = 0;
  }
  
  if (outname[0] == 0) {
    printf("ERROR: No output file specified\n");
    errorFlag = 1;
  }
			 

  if (errorFlag) {
    printf("Usage: mapdistortion [-v] [-d] [-g] [-a] [-s] [-r <search pixels>] [-e <solution_number>] [-t <arcsec>] <filename>\n");
    printf(" Returns: center ra (degrees)\n");
    printf("          center dec (degrees)\n");
    printf("          max radius (arcmin)\n");
    printf("          scale (arcsec/pixel)\n");
    printf("          search radius (<search pixels>*scale) in sexagesimal\n");
    printf("          -a display the axis length in pixels\n");
    printf("             (width in pixels, height in pixels, arcsec/pixel, scan pattern id\n"); 
    printf("             left, right, bottom, top margins in pixels\n");
    printf("          -b display 'small', 'normal', or 'large'; series, plateNumber, mosaicNumber and filename\n");
    printf("             where small = 4in x 5in, normal = 8in x 10in, and large = 14in x 17in\n");
    printf("          -d outputs in degrees\n");
    printf("          -l outputs galactic coordinates (lon lat) in degrees\n");
    printf("          -g generate gsc2.3 query\n");
    printf("          -r number of pixels in search radius (default is 40)\n");
    printf("          -s display the series, arcsec/pixel, bin a1 radius in degrees\n");
    printf("          -v is the verbosity flag\n");
    printf("          -p ignore rotation\n");
    printf("          -e <solution number>\n");
    printf("          -t <arcsec> displays pixel width only for getfits\n");
    return(-1);
  }
#endif

  mysql_init(&my_connection);

  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    printf("ERROR: DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    printf("ERROR: DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    printf("ERROR: DASCH_PASSWORD is not defined\n");
    return(-1);
  }
                   

  if (mysql_real_connect(&my_connection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
#if 0
    printf("Connection success\n");
#endif
  } else {
    printf("ERROR: mapdistortion Connection failed\n");
    if (mysql_errno(&my_connection)) {
      printf("ERROR: mapdistortion Connection error %d: %s\n",mysql_errno(&my_connection),mysql_error(&my_connection));
    }
    return(-1);
  }


  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    fprintf(stderr,"DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }

  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    fprintf(stderr,"DASCH_PHOT_USERNAME is not defined\n");
    return(-1);
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    fprintf(stderr,"DASCH_PHOT_PASSWORD is not defined\n");
    return(-1);
  }
  mysql_init(pPhotConnection);
  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
    }
    return(-1);
  }



  strcpy(outname,"/home/scanner/Pipeline/losmap.db");
  strcpy(scalename,"/home/scanner/backup/2018_06_14/getscale.log");
#if 1
  strcpy(filename,"/dasch/data/Pipeline/headers/ax04044_01_01ww_tnx.hdr"); /* This is the very best fit so far */
#endif
#if 0
  strcpy(filename,"/dasch/data/Pipeline/headers/ax04832_00_01r180ww_tnx.hdr");
#endif
#if 0
  strcpy(filename,"/dasch/data/Pipeline/headers/ax04832_00_01r180ww_TAN.hdr");
#endif


  outhandle = fopen(outname,"wt");
  if (outhandle == NULL) {
    printf("ERROR: failed to open output file %s\n",outname);
    exit(-1);
  }
#ifdef SUMMARY_PLOT
  fprintf(outhandle,"xIndex\tpixDist\tdrad\tscale\tPlate\n");
  fprintf(outhandle,"------\t-------\t----\t-----\t-----\n");
#else /* SUMMARY_PLOT */
  fprintf(outhandle,"xIndex\tX_IMAGE\tY_IMAGE\tra\tdec\tpixDist\tradius\tdrad\tscale\tPlate\n");
  fprintf(outhandle,"------\t-------\t-------\t--\t---\t-------\t------\t----\t-----\t-----\n");
#endif /* SUMMARY_PLOT */

  scalehandle = fopen(scalename,"rt");
  if (scalehandle == NULL) {
    printf("ERROR: failed to open input file %s line %d\n",scalename,__LINE__);
    exit(-1);
  }
  

  /* Count the number of records in the file */
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,scalehandle);
    if (inBuffer == NULL) {
      break;
    }
    lineCounter++;
    lineLen = strlen(inBuffer);
    /* Trim off the carriage return */
    if (inBuffer[lineLen-1] == 10) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    /* Trim off the line feed */
    if (inBuffer[lineLen-1] == 13) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    if (lineLen > 5) {
      if (strstr(inBuffer,"dot product") != NULL) {
        maxPlates++;
      }
      charPtr = strstr(inBuffer,"unbinned arcsec/mm");
      if (charPtr != NULL) {
        if (strstr(inBuffer,"unbinned arcsec/mm numPoly") == NULL) {
          charPtr += strlen("unbinned arcsec/mm");
          nvals = sscanf(charPtr,"%lf",&scale);
          if (nvals == 1) {
            scale = scale * NOMINAL_MM_PER_PIXEL; /* scale in arcsec per unbinned pixel */
          } else {
            scale = 0;
            printf("ERROR: could not decode scale in line %d of %s\n",lineCounter,scalename);
            exit(-1);
          }
        }
      }
      if (maxString < lineLen) {
        maxString = lineLen;
      }
    }
  }
  fclose(scalehandle);
  printf("Found %d plates and a scale %f arcsec/pixel maxString %d\n",maxPlates,scale,maxString);
  plateTable = (PPLATE)calloc(maxPlates,sizeof(PLATE));
  if (plateTable == NULL) {
    printf("ERROR: failed to allocate plateTable\n");
    exit(-1);
  }

  /* Now reopen the file and parse the mosaic specific information */


  scalehandle = fopen(scalename,"rt");
  if (scalehandle == NULL) {
    printf("ERROR: failed to open input file %s line %d\n",scalename,__LINE__);
    exit(-1);
  }
  
  lineCounter = 0;
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,scalehandle);
    if (inBuffer == NULL) {
      break;
    }
    lineCounter++;
    lineLen = strlen(inBuffer);
    /* Trim off the carriage return */
    if (inBuffer[lineLen-1] == 10) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    /* Trim off the line feed */
    if (inBuffer[lineLen-1] == 13) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    /* Trim trailing spaces */
    if (inBuffer[lineLen-1] == ' ') {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    if (lineLen > 5) {
      if (strstr(inBuffer,"dot product") != NULL) {
        pPlate = &plateTable[plateCount];
        memset(pPlate,0,sizeof(PLATE));

        strcpy(copyLine,inBuffer);
        charPtr = strchr(inBuffer,' ');
        if (charPtr == NULL) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }
        *charPtr = 0;
        if (strlen(inLine) >= MAX_PLATE_NAME) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }
        strcpy(pPlate->Plate,inBuffer);
        charPtr++;
        charPtr2 = strstr(charPtr,"flag");
        if (charPtr2 == NULL) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }
        charPtr2 += strlen("flag");
        nvals = sscanf(charPtr2,"%d",&pPlate->flag);
        if (nvals != 1) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }

        charPtr = charPtr2;
        charPtr2 = strstr(charPtr,"fit");
        if (charPtr2 == NULL) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }
        charPtr2 += strlen("fit");
        nvals = sscanf(charPtr2,"%d",&pPlate->fit);
        if (nvals != 1) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }

        charPtr = charPtr2;
        charPtr2 = strstr(charPtr,"PolyRefStars");
        if (charPtr2 == NULL) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }
        charPtr2 += strlen("PolyRefStars");
        nvals = sscanf(charPtr2,"%d",&pPlate->PolyRefStars);
        if (nvals != 1) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }

        charPtr = charPtr2;
        charPtr2 = strstr(charPtr,"PolyRMS");
        if (charPtr2 == NULL) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }
        charPtr2 += strlen("PolyRMS");
        nvals = sscanf(charPtr2,"%lf",&pPlate->PolyRMS);
        if (nvals != 1) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }

        charPtr = charPtr2;
        charPtr2 = strstr(charPtr,pPlate->Plate);
        if (charPtr2 == NULL) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }
        if (strlen(charPtr2) >= MAX_PLATE_NAME) {
          printf("ERROR: format error in software line %d at line %d of %s\n",__LINE__,lineCounter,scalename);
          continue;
        }
        strcpy(pPlate->platename,charPtr2);
        if (plateCount >= maxPlates) {
          printf("ERROR: plateCount %d exceeds maxPlates %d\n",plateCount,maxPlates);
          exit(-1);
        }

        plateCount++;
      }
    }
  }
  fclose(scalehandle);
  printf("Found %d of %d plates and a scale %f arcsec/pixel maxString %d\n",plateCount,maxPlates,scale,maxString);
  qsort((void*)plateTable,plateCount,sizeof(PLATE),PlateCompare);

  memset(pixDistSummaryTable,0,MAX_PIXDIST*sizeof(PIXDIST));
  for (plateIndex = 0; plateIndex < plateCount; plateIndex++) {
    pPlate = &plateTable[plateIndex];
    strcpy(filename,headersDirectory);
    strcat(filename,"/");
    strcat(filename,pPlate->platename);
    strcat(filename,"_tnx.hdr");
    statResult = stat(filename,&statbuf);
    if (statResult < 0) {
      printf("ERROR: failed to find file %s\n",filename);
      continue;
    }

    if (ParseFilename(filename,pPlate->series,&pPlate->plateNumber,&pPlate->mosaicNumber,&pPlate->binning,&pPlate->rotation) == 0) {
      printf("ERROR: failed to parse %s\n",filename);
      continue;
    }
    if (pPlate->binning != 1) {
      printf("ERROR: mapdistortion File binning must be 1\n");
      return(-1);
    }
    for (local_bin_index = 0; local_bin_index < TOTAL_DMAGBINS_NORMAL; local_bin_index++) {
      pPlate->rejectFlag[local_bin_index] = -1;
    }
    sprintf(queryString,"SELECT seriesId,local_bin_index,rejectFlag,drad_bin_size2 FROM localbin INNER JOIN photseries USING(seriesId) where series = '%s' and plateNumber = %d and solutionNumber = 0;",
            pPlate->series,
            pPlate->plateNumber);
    if (strlen(queryString) > (MAX_QUERY_STRING-2)) {
      printf("ERROR: queryString exceeded  %d\n",strlen(queryString));
      exit(-1);
    }
    res = ExecuteQuery(pPhotConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pPhotConnection);
      if (res_ptr) {
        numRows = mysql_affected_rows(pPhotConnection);
        if (numRows != TOTAL_DMAGBINS_NORMAL) {
          printf("ERROR: Plate %s%05d expects %d localbin rows and has %d localbin rows\n",pPlate->series,pPlate->plateNumber,TOTAL_DMAGBINS_NORMAL,numRows);
        } else {
    
          while ((sqlrow = mysql_fetch_row(res_ptr))) {
            if (sqlrow[0]) {
              nvals = sscanf(sqlrow[0],"%d",&pPlate->seriesId);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for seriesId in GetPhotLocalBin\n");
                continue;
              }
            } else {
              printf("ERROR: NULL seriesId in GetPhotLocalBin\n");
              continue;
            }
            if (sqlrow[1]) {
              nvals = sscanf(sqlrow[1],"%d",&local_bin_index);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for local_bin_index\n");
                continue;
              }
            } else {
              printf("ERROR: NULL local_bin_index in GetPhotLocalBin\n");
              continue;
            }
            if ((local_bin_index < 0) || (local_bin_index >= TOTAL_DMAGBINS_NORMAL)) {
              printf("ERROR illegal local_bin_index %d\n",local_bin_index);
              continue;
            }
            if (sqlrow[2]) {
              nvals = sscanf(sqlrow[2],"%d",&pPlate->rejectFlag[local_bin_index]);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for rejectFlag\n");
                continue;
              }
            }
            if (sqlrow[3]) {
              nvals = sscanf(sqlrow[3],"%d",&pPlate->drad_bin_size2[local_bin_index]);
              if (nvals != 1) {
                printf("ERROR: nvals is %d for drad_bin_size2\n");
                continue;
              }
            }
  
            curRow++;
          }
        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
      }
    }




    if ((header = GetFITShead (filename, verbose)) == NULL) {
      printf("Failed to read the header for %s\n",filename);
      return(-1);
    }
    wcs = GetFITSWCS (filename, header, verbose, &cra, &cdec, &dra, &ddec,
                      &secpix, &wp, &hp, &sysim, &eqim);
    if (nowcs (wcs)) {
      wcsfree (wcs);
      wcs = NULL;
      fprintf(stderr,"ERROR: no wcs for %s\n",filename);
      continue;
    }
    if (wcs->pvfail) {
      printf("ERROR: Significant inaccuracy likely to occur in projection for %s\n",filename);
      wcsfree(wcs);
      continue;
    }
    naxis1 = wcs->nxpix;
    naxis2 = wcs->nypix;
    xctr = 0.5 + ((1.0*wcs->nxpix)/2.0);
    yctr = 0.5 + ((1.0*wcs->nypix)/2.0);
    pix2wcs(wcs,xctr,yctr,&cra,&cdec);
    if (wcs->nxpix < wcs->nypix) {
      denominator = wcs->nxpix;
    } else {
      denominator = wcs->nypix;
    }

    xOffset = (wcs->crpix[0]-xctr)/denominator;
    yOffset = (wcs->crpix[1]-yctr)/denominator;
    fullOffset = sqrt(sqr(xOffset)+sqr(yOffset));
#ifdef MAX_OFFSET
    if (fullOffset < MAX_OFFSET) {
      wcsfree(wcs);
      continue;
    }
#endif /* MAX_OFFSET */

#ifdef DEMO_PLOT
    xplatectr = 0.5 + ((1.0*wcs->nxpix)/2.0);
    yplatectr = 0.5 + ((1.0*wcs->nypix)/2.0);
    while (1) {
      if (demoIteration > 0) {
        yctr = 0.5 + ((1.0*wcs->nypix)/2.0);
        yctr += (1.0*demoIteration*naxis1)/200.;
      } else {
        if (demoCounter == 0) {
          yctr = 0.5 + ((1.0*wcs->nypix)/2.0);
          yctr += (1.0*naxis1)/20.;
        }
      }
      wcsfree(wcs); 
      cd[0] = -scale/3600.;
      cd[1] = 0.0;
      cd[2] = 0.0;
      cd[3] = scale/3600.;
      wcs =  wcskinit(naxis1,
                      naxis2,
                      "RA---TAN",
                      "DEC--TAN",
                      xctr,
                      yctr,
                      cra,
                      cdec,
                      cd,
                      0,  /* cdelt1 */
                      0,  /* cdelt2 */
                      0,  /* crota */
                      2000, /* equinox */
                      0);   /* epoch */
      if (demoIteration > 0) { 
        pix2wcs(wcs,xctr,yctr,&ra1,&dec1);
        pix2wcs(wcs,xctr,yctr+1.,&ra2,&dec2);
        scale1 = wcsdist(ra1,dec1,ra2,dec2)* 3600.;
        pix2wcs(wcs,xplatectr,yplatectr,&ra3,&dec3);
        pix2wcs(wcs,xplatectr,yplatectr+1,&ra4,&dec4);
        scale2 = wcsdist(ra3,dec3,ra4,dec4)* 3600.;
        printf("iteration %3d yctr %f height %f offset scale %f  center scale %f ratio %f\n",demoIteration,yctr,1.0*wcs->nypix,scale1,scale2,scale1/scale2);
        demoIteration--;
      }
      if (demoIteration < 0) {
        break;
      }
    }
  
#endif /* DEMO_PLOT */

    printf("plateIndex %4d of %4d  xctr %f yctr %f cra %f cdec %f for %s\n",plateIndex,plateCount,xctr,yctr,cra,cdec,filename);
    memset(vecA,0,sizeof(vecA));
    memset(vecB,0,sizeof(vecB));
    memset(vecC,0,sizeof(vecC));
    memset(vecD,0,sizeof(vecD));
    memset(matA,0,sizeof(matA));
    memset(matB,0,sizeof(matB));
    memset(matC,0,sizeof(matC));

    cosx = cos(DEGREES_TO_RAD*cdec);
    vecA[0] = cosx*cos(DEGREES_TO_RAD*cra);
    vecA[1] = cosx*sin(DEGREES_TO_RAD*cra);
    vecA[2] = sin(DEGREES_TO_RAD*cdec);

    resolv(pCommon,vecA[0],vecA[1],vecA[2],&decinv,&rainv);
    if (verboseFlag) {
      printf("plateIndex %4d rainv %f decinv %f\n",plateIndex,rainv,decinv);
    }
 
    /* Rotate to the first point of Aries */
    rotate(pCommon,3,cra,matA);
    gmprd(pCommon,matA,vecA,vecB,3,3,1);
    resolv(pCommon,vecB[0],vecB[1],vecB[2],&decinv,&rainv);
    if (verboseFlag) {
      printf("plateIndex %4d rainv %f decinv %f\n",plateIndex,rainv,decinv);
    }
    /* Now rotate to the celestial equator */
    rotate(pCommon,2,-cdec,matB);
    gmprd(pCommon,matB,vecB,vecC,3,3,1);
    resolv(pCommon,vecC[0],vecC[1],vecC[2],&decinv,&rainv);
    if (verboseFlag) {
      printf("plateIndex %4d rainv %f decinv %f\n",plateIndex,rainv,decinv);
    }
    pix2wcs(wcs,xctr,yctr+1,&cra2,&cdec2);
    if (verboseFlag) {
      printf("plateIndex %4d xctr %f yctr %f cra2 %f cdec2 %f for %s\n",plateIndex,xctr,yctr,cra2,cdec2,filename);
    }
    cosx = cos(DEGREES_TO_RAD*cdec2);
    vecA[0] = cosx*cos(DEGREES_TO_RAD*cra2);
    vecA[1] = cosx*sin(DEGREES_TO_RAD*cra2);
    vecA[2] = sin(DEGREES_TO_RAD*cdec2);
    gmprd(pCommon,matA,vecA,vecB,3,3,1);
    gmprd(pCommon,matB,vecB,vecC,3,3,1);
    resolv(pCommon,vecC[0],vecC[1],vecC[2],&decinv2,&rainv2);
    if (verboseFlag) {
      printf("plateIndex %4d rainv2 %f decinv2 %f\n",plateIndex,rainv2,decinv2);
    }
    if (decinv2 == 0) {
      if (rainv2 > 0) {
        verticalang = 90.0;
      } else {
        verticalang = -90.0;
      }
    } else {
      verticalang = atan(rainv2/decinv2)*RAD_TO_DEGREES;
    }
    rotate(pCommon,1,-verticalang,matC);
    gmprd(pCommon,matC,vecC,vecD,3,3,1);
    resolv(pCommon,vecD[0],vecD[1],vecD[2],&decinv3,&rainv3);
    if (verboseFlag) {
      printf("plateIndex %4d rainv3 %f decinv3 %f verticalang %f\n",plateIndex,rainv3,decinv3,verticalang);
    }

    /* First find our angles */
    xloc[CENTER_POINT] = 0.5 + ((1.0*wcs->nxpix)/2.0);
    yloc[CENTER_POINT] = 0.5 + ((1.0*wcs->nypix)/2.0);
    pix2wcs(wcs,xloc[CENTER_POINT],yloc[CENTER_POINT],&actualra[CENTER_POINT],&actualdec[CENTER_POINT]);
    if (verboseFlag) {
      printf("plateIndex %4d xctr[CENTER_POINT] %f yloc[CENTER_POINT] %f cra %f cdec %f for %s\n",plateIndex,xloc[CENTER_POINT],yloc[CENTER_POINT],actualra[CENTER_POINT],actualdec[CENTER_POINT],filename);
    }
    cosx = cos(DEGREES_TO_RAD*actualdec[CENTER_POINT]);
    val1[CENTER_POINT][0] = cosx*cos(DEGREES_TO_RAD*actualra[CENTER_POINT]);
    val1[CENTER_POINT][1] = cosx*sin(DEGREES_TO_RAD*actualra[CENTER_POINT]);
    val1[CENTER_POINT][2] = sin(DEGREES_TO_RAD*actualdec[CENTER_POINT]);


    resolv(pCommon,val1[CENTER_POINT][0],val1[CENTER_POINT][1],val1[CENTER_POINT][2],&decinv,&rainv);
    if (verboseFlag) {
      printf("plateIndex %4d rainv %f decinv %f\n",plateIndex,rainv,decinv);
    }
 
    /* Rotate to the first point of Aries */
    rotate(pCommon,3,actualra[CENTER_POINT],matA);
    gmprd(pCommon,matA,val1[CENTER_POINT],val2[CENTER_POINT],3,3,1);
    resolv(pCommon,val2[CENTER_POINT][0],val2[CENTER_POINT][1],val2[CENTER_POINT][2],&decinv,&rainv);
    if (verboseFlag) {
      printf("plateIndex %4d rainv %f decinv %f\n",plateIndex,rainv,decinv);
    }
    /* Now rotate to the celestial equator */
    rotate(pCommon,2,-actualdec[CENTER_POINT],matB);
    gmprd(pCommon,matB,val2[CENTER_POINT],val1[CENTER_POINT],3,3,1);
    resolv(pCommon,val1[CENTER_POINT][0],val1[CENTER_POINT][1],val1[CENTER_POINT][2],&decinv,&rainv);
    if (verboseFlag) {
      printf("plateIndex %4d rainv[CENTER_POINT] %f decinv[CENTER_POINT] %f\n",plateIndex,rainv,decinv);
    }

    /* Now do the same for our test point */
    xloc[TEST_POINT] = xloc[CENTER_POINT];
    yloc[TEST_POINT] = yloc[CENTER_POINT]+1;
    pix2wcs(wcs,xloc[TEST_POINT],yloc[TEST_POINT],&actualra[TEST_POINT],&actualdec[TEST_POINT]);
    if (verboseFlag) {
      printf("plateIndex %4d xloc[TEST_POINT] %f yloc[TEST_POINT] %f cra2 %f cdec2 %f for %s\n",plateIndex,xloc[TEST_POINT],yloc[TEST_POINT],actualra[TEST_POINT],actualdec[TEST_POINT],filename);
    }
    cosx = cos(DEGREES_TO_RAD*actualdec[TEST_POINT]);
    val1[TEST_POINT][0] = cosx*cos(DEGREES_TO_RAD*actualra[TEST_POINT]);
    val1[TEST_POINT][1] = cosx*sin(DEGREES_TO_RAD*actualra[TEST_POINT]);
    val1[TEST_POINT][2] = sin(DEGREES_TO_RAD*actualdec[TEST_POINT]);
    gmprd(pCommon,matA,val1[TEST_POINT],val2[TEST_POINT],3,3,1);
    gmprd(pCommon,matB,val2[TEST_POINT],val1[TEST_POINT],3,3,1);
    resolv(pCommon,val1[TEST_POINT][0],val1[TEST_POINT][1],val1[TEST_POINT][2],&decinv2,&rainv2);
    if (verboseFlag) {
      printf("plateIndex %4d rainv2 %f decinv2 %f\n",plateIndex,rainv2,decinv2);
    }
    if (decinv2 == 0) {
      if (rainv2 > 0) {
        vertang = 90.0;
      } else {
        vertang = -90.0;
      }
    } else {
      vertang = atan(rainv2/decinv2)*RAD_TO_DEGREES;
    }
    rotate(pCommon,1,-vertang,matC);
    gmprd(pCommon,matC,val1[TEST_POINT],val2[TEST_POINT],3,3,1);
    resolv(pCommon,val2[TEST_POINT][0],val2[TEST_POINT][1],val2[TEST_POINT][2],&mappeddec[TEST_POINT],&mappedra[TEST_POINT]);
    if (verboseFlag) {
      printf("plateIndex %4d mappedra[TEST_POINT] %f mappeddec[TEST_POINT] %f vertang %f\n",plateIndex,mappedra[TEST_POINT],mappeddec[TEST_POINT],vertang);
    } 
    /* Now check the center point again */
  
    if (verboseFlag) {
      printf("plateIndex %4d val1[CENTER_POINT][0] %f val1[CENTER_POINT][1] %f val1[CENTER_POINT][2] %f line %d\n",plateIndex,val1[CENTER_POINT][0],val1[CENTER_POINT][1],val1[CENTER_POINT][2],__LINE__);
      printf("plateIndex %4d rainv[CENTER_POINT] %f decinv[CENTER_POINT] %f\n",plateIndex,rainv,decinv);
    }
    gmprd(pCommon,matC,val1[CENTER_POINT],val2[CENTER_POINT],3,3,1);
    resolv(pCommon,val1[CENTER_POINT][0],val1[CENTER_POINT][1],val1[CENTER_POINT][2],&decinv,&rainv);
    resolv(pCommon,val2[CENTER_POINT][0],val2[CENTER_POINT][1],val2[CENTER_POINT][2],&mappeddec[CENTER_POINT],&mappedra[CENTER_POINT]);
    if (verboseFlag) {
      printf("plateIndex %4d mappedra[CENTER_POINT] %f mappeddec[CENTER_POINT] %f vertang %f\n",plateIndex,mappedra[CENTER_POINT],mappeddec[CENTER_POINT],vertang);
    }
    radius = wcsdist(mappedra[CENTER_POINT],mappeddec[CENTER_POINT],mappedra[TEST_POINT],mappeddec[TEST_POINT]);
    if (verboseFlag) {
      printf("plateIndex %4d Plate scale is %f arcsec/pixel %f degrees/pixel\n",plateIndex,3600*scale,scale); 
    }
    cd[0] = -scale;
    cd[1] = 0.0;
    cd[2] = 0.0;
    cd[3] = scale;
    tan_wcs =  wcskinit(wcs->nxpix,
                        wcs->nypix,
                        "RA---TAN",
                        "DEC--TAN",
                        xloc[CENTER_POINT],
                        yloc[CENTER_POINT],
                        0.0,
                        0.0,
                        cd,
                        0,  /* cdelt1 */
                        0,  /* cdelt2 */
                        0,  /* crota */
                        2000, /* equinox */
                        0);   /* epoch */


    /* Now handle all of the smoothing bins */
    for (pixIndex = 0; pixIndex < MAX_PIXDIST; pixIndex++) {
      pPixDist = &pixDistTable[pixIndex];
      pPixDist->count = 0;
      pPixDist->value = 0;
      pPixDist->xOffset = 0;
      pPixDist->yOffset = 0;
      pPixDist->vectorCount = 0;
    }
    
    for (local_bin = 0; local_bin < TOTAL_DMAGBINS_NORMAL; local_bin++) {
      /* Find the bin centers */
      ix = (local_bin % nx);
      iy = (local_bin / ny);
      xPixelsPerBin = (1.0*wcs->nxpix)/(1.0*nx);
      yPixelsPerBin = (1.0*wcs->nypix)/(1.0*nx);

      xloc[local_bin] =  xloc[CENTER_POINT]+(xPixelsPerBin*(0.5+(1.0*(ix-(nx/2)))));
      yloc[local_bin] =  yloc[CENTER_POINT]+(yPixelsPerBin*(0.5+(1.0*(iy-(ny/2)))));
      pix2wcs(wcs,xloc[local_bin],yloc[local_bin],&actualra[local_bin],&actualdec[local_bin]);
      cosx = cos(DEGREES_TO_RAD*actualdec[local_bin]);
      val1[local_bin][0] = cosx*cos(DEGREES_TO_RAD*actualra[local_bin]);
      val1[local_bin][1] = cosx*sin(DEGREES_TO_RAD*actualra[local_bin]);
      val1[local_bin][2] = sin(DEGREES_TO_RAD*actualdec[local_bin]);
      gmprd(pCommon,matA,val1[local_bin],val2[local_bin],3,3,1);
      gmprd(pCommon,matB,val2[local_bin],val1[local_bin],3,3,1);
      gmprd(pCommon,matC,val1[local_bin],val2[local_bin],3,3,1);
      resolv(pCommon,val2[local_bin][0],val2[local_bin][1],val2[local_bin][2],&mappeddec[local_bin],&mappedra[local_bin]);
      pixDist = sqrt(sqr(xloc[local_bin]-xloc[CENTER_POINT])+sqr(yloc[local_bin]-yloc[CENTER_POINT]));
      radius = wcsdist(mappedra[CENTER_POINT],mappeddec[CENTER_POINT],mappedra[local_bin],mappeddec[local_bin]);
      drad = RAD_TO_DEGREES*(cos((90.0-radius)*DEGREES_TO_RAD)/sin((90.0-radius)*DEGREES_TO_RAD));
#ifdef DEMO_PLOT
      pPlate->rejectFlag[local_bin] = 0;
      pPlate->drad_bin_size2[local_bin] = 1;
#endif /* DEMO_PLOT */

      if ((pPlate->rejectFlag[local_bin] == 0) && (pPlate->drad_bin_size2[local_bin] == 1)) {

        goodLocalBins++;
#ifdef SUMMARY_PLOT
        pixIndex = pixDist;
        if ((pixIndex < 0) || (pixIndex >= MAX_PIXDIST)) {
          printf("ERROR: illegal pixIndex %d\n",pixIndex);
          exit(-1);
        }
        pPixDist = &pixDistTable[pixIndex];
        pPixDist->count++;
        if (pPixDist->count > maxPixDistCount) {
          maxPixDistCount = pPixDist->count;
        }
        pPixDist->value += ((drad/(scale/3600.))-pixDist);
        if (pPixDist->vector == NULL) {
          pPixDist->vector = (double *)calloc(MAX_VECTOR_COUNT,sizeof(double));
          pixDistVectorAllocCount++;
        }
        if (pPixDist->vector == NULL) {
          if (pixDistVectorError == 0) {
            printf("ERROR: failed to allocate pixDist vector at allocation %d\n",pixDistVectorAllocCount);
            pixDistVectorError++;
          }
        } else {
          if (pPixDist->vectorCount > (MAX_VECTOR_COUNT-1)) {
            if (pixDistCountError == 0) {
              printf("ERROR: MAX_VECTOR_COUNT exceeded\n");
              pixDistCountError++;
            } 
          } else {
            pPixDist->vector[pPixDist->vectorCount] = ((drad/(scale/3600.))-pixDist);
            pPixDist->vectorCount++;
          }
            
        }

        xOffset =  xloc[local_bin]-xloc[CENTER_POINT];
        yOffset =  yloc[local_bin]-yloc[CENTER_POINT];
        pPixDist->xOffset += xOffset;
        pPixDist->yOffset += yOffset;
#if 0
        printf("local_bin %d pixIndex %d, count %d value %f %f, xOffset %f yOffset %f\n",local_bin,pixIndex,pPixDist->count,((drad/(scale/3600.))-pixDist),pPixDist->value,xOffset,yOffset);
#endif


#else /* SUMMARY_PLOT */
        fprintf(outhandle,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",plateIndex,xloc[local_bin],yloc[local_bin],mappedra[local_bin],mappeddec[local_bin],pixDist,radius,drad,scale,pPlate->Plate);
#endif /* SUMMARY_PLOT */

      }
      /* plot   ((drad/0.002079061754146236)-pixDist) as a function of radius */
      /* plot   ((drad/(scale/3600.))-pixDist) as a function of radius */

    }
    wcsfree(wcs);
    wcsfree(tan_wcs);
    platesProcessed++;
#ifdef SUMMARY_PLOT
    for (pixIndex = 0; pixIndex < MAX_PIXDIST; pixIndex++) {
      pPixDist = &pixDistTable[pixIndex];
      pPixDistSummary = &pixDistSummaryTable[pixIndex];

      if (pPixDist->count == 0) {
        continue;
      } else if ((pPixDist->count % 4 == 0) &&
                 (fabs(pPixDist->xOffset) < 1.0) && 
                 (fabs(pPixDist->yOffset) < 1.0)) {        
#if 0
        printf("pixIndex %d, count %d value %f, xOffset %f yOffset %f\n",pixIndex,pPixDist->count,pPixDist->value,pPixDist->xOffset,pPixDist->yOffset);
#endif

        fprintf(outhandle,"%d\t%d\t%f\t%f\t%s\n",plateIndex,pixIndex, pPixDist->value/(1.0*pPixDist->count),scale,pPlate->Plate);
        
        pPixDistSummary->count += pPixDist->count;
        pPixDistSummary->value += pPixDist->value;

        if (pPixDistSummary->count > maxPixDistSummaryCount) {
          maxPixDistSummaryCount = pPixDistSummary->count;
        }
        if (pPixDistSummary->vector == NULL) {
          pPixDistSummary->vector = (double *)calloc(MAX_VECTOR_COUNT,sizeof(double));
          pixDistVectorAllocCount++;
        }
        if (pPixDistSummary->vector == NULL) {
          if (pixDistVectorError == 0) {
            printf("ERROR: failed to allocate pixDist vector at allocation %d\n",pixDistVectorAllocCount);
            pixDistVectorError++;
          }
        } else {
          if ((pPixDistSummary->vectorCount+pPixDist->vectorCount) > (MAX_VECTOR_COUNT-1)) {
            if (pixDistCountError == 0) {
              printf("ERROR: MAX_VECTOR_COUNT exceeded\n");
              pixDistCountError++;
            } 
          } else {
            for (vectorIndex = 0; vectorIndex < pPixDist->vectorCount; vectorIndex++) { 
              pPixDistSummary->vector[pPixDistSummary->vectorCount] = pPixDist->vector[vectorIndex];
              pPixDistSummary->vectorCount++;
            }
            pPixDist->vectorCount = 0;
          }
            
        }

      } else if (pPixDist->count >= 16) {
        printf("ERROR: illegal pixDistCount %d\n",pPixDist->count);
        exit(-1);
      }
    }
#endif /* SUMMARY_PLOT */


#ifdef PLATE_LIMIT
    if (platesProcessed > PLATE_LIMIT) {
      break;
    }
#endif /* POLYREFSTARS_LIMIT */
#ifdef DEMO_PLOT
    demoCounter++;
    if (demoCounter > 0) {
      break;
    }
#endif /* DEMO_PLOT */

  }

#ifdef SUMMARY_PLOT
  for (pixIndex = 0; pixIndex < MAX_PIXDIST; pixIndex++) {
    pPixDistSummary = &pixDistSummaryTable[pixIndex];
    if (pPixDistSummary->count == 0) {
      continue;
    } else {        
      bestCount = CalcMedianAndRMS(pPixDistSummary->vectorCount,10,pPixDistSummary->vector,&pPixDistSummary->median,&pPixDistSummary->rms,0,3.0,0);
      if (bestCount == 0) {
        pPixDistSummary->median = -10000.0;
        pPixDistSummary->rms = -99.0;
      } else {
        if (pPixDistSummary->rms < MAX_RMS) {
          fprintf(outhandle,"%d\t%d\t%f\t%f\t%s\n",-1,pixIndex, pPixDistSummary->value/(1.0*pPixDistSummary->count),scale,"summary");
          fprintf(outhandle,"%d\t%d\t%f\t%f\t%s\n",-2,pixIndex, pPixDistSummary->median,scale,"median");
          fprintf(outhandle,"%d\t%d\t%f\t%f\t%s\n",-3,pixIndex, pPixDistSummary->rms,scale,"rms");
        }
      }
    }
  }
#endif /* SUMMARY_PLOT */


  mysql_close(&my_connection);
  mysql_close(pPhotConnection);

  if (plateTable != NULL) {
    free(plateTable);
  }
  for (pixIndex = 0; pixIndex < MAX_PIXDIST; pixIndex++) {
    pPixDist = &pixDistTable[pixIndex];
    if (pPixDist->vector != NULL ) {
      free(pPixDist->vector);
      pPixDist->vector = NULL;
    }
    pPixDistSummary = &pixDistSummaryTable[pixIndex];
    if (pPixDistSummary->vector != NULL ) {
      free(pPixDistSummary->vector);
      pPixDistSummary->vector = NULL;
    }
  }


  time(&curTime);
  curTime = curTime - startTime;
  printf("Finished goodLocalBins %d maxPixDistCount %d maxPixDistSummaryCount %d pixDistVectorAllocCount %d MAX_RMS %d in %d seconds\n",goodLocalBins,maxPixDistCount,maxPixDistSummaryCount,pixDistVectorAllocCount,MAX_RMS,curTime);
  if ((pixDistVectorError != 0) || (pixDistCountError != 0)) {
    printf("ERROR pixDistVectorError %d pixDistCountError %d\n",pixDistVectorError,pixDistCountError);
  }
  return(EXIT_SUCCESS);

} 
