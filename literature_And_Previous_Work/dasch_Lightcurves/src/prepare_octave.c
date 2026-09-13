// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* prepare_octave.c
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I /dasch/install/include  -L /dasch/install/lib -lm  -L/usr/lib/mysql -L/usr/lib64/mysql  -l mysqlclient prepare_octave.c pipelineutils.a -ltable -lutil -lwcs -o prepare_octave
 *
 *   Fields from match_${table}
 *    REF  ra_2   dec_2  Stdmag  color  MAG_ISO  X_IMAGE  Y_IMAGE ddec dra NUMBER FLAGS 
 *
 *
 *   Output Fields: X_WIDTH and Y_WIDTH are the width and height in pixels of the mosaic.
 *   1: REF 2: ra_2   3: dec_2  4: Stdmag+extinction  5: color  6: MAG_ISO  7: NUMBER 8: FLAGS
 *   9: X_IMAGE 10: Y_IMAGE 11: ddec 12: dra 13: extinction 14: spatial_bin 15: local_bin 16: X_WIDTH  17: Y_WIDTH 18: PLATE_SCALE * 
 *
 * Mar  5, 2008 Edward J. Los - Initial version
 * Mar 10, 2008 Edward J. Los - Increase time tolerance to 1.1 hours.
 * Mar 28, 2008 Edward J. Los - Remove MAG_APER, ISOAREA, and plate_dist for performance.  Switch to ra_2 and dec_2
 * Mar 31, 2008 Edward J. Los - Correct column bug: drad should be dra
 *                              Use native octave mode for performance.
 * Apr 16, 2008 Edward J. Los - Correct extinction error when writing native mode files
 * Jul  3, 2008 Edward J. Los - Add FILTER_FLAG_BAD
 * Aug 13, 2008 Edward J. Los - Correct bug in updating Stdmag for extinction.
 * Dec 22, 2008 Edward J. Los - Eliminate High extinction objects (> 2.5 airmass)
 * Jan 27, 2009 Edward J. Los - split FLAGS into AFLAGS and BFLAGS
 * Feb  4, 2009 Edward J. Los - Add FILTER_BFLAG_BAD_JD for plates with uncertain Julian Days 
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 * Aug 20, 2009 Edward J. Los - Add Kepler Input Catalog support
 * Nov  4, 2009 Edward J. Los - Add multiple exposure support. Add the JulianDate and timeAccuracy to the extinction file to be sure
 *                              the solution does not change as more exposures are solved.  
 * Dec 15, 2009 Edward J. Los - Clear the StaleExposureNumber bit because we are recalculating the extinction
 * Jan 19, 2010 Edward J. Los - Implement a new extinction formula using GetExtinctionCoefficient
 * Mar  7, 2011 Edward J. Los - Reject stars with colors outside reasonable bounds.
 * Mar 14, 2011 Edward J. Los - Calculate the spatial bin index and local bin index here
 * Mar 23, 2011 Edward J. Los - Add apass catalog support
 * Jul 30, 2011 Edward J. Los - Add experimental catalog support
 * May  7, 2014 Edward J. Los - If all of the of the local bins have more than MAX_AIRMASS, then do not set the FILTER_AFLAG_LOW_ALTITUDE until the recover_points stage.
 * Jun  6, 2014 Edward J. Los - Change "setting extinction to zero" to "no points have an airmass greater than than MAX_AIRMASS"
 * Jul 25, 2017 Edward J. Los - Make the number of extinction bins variable and identical to the number of local bins as required by the localbins tables
 *                              in the photometry database.
 * May 28, 2018 Edward J. Los - Support the gaia calibration.
 * Oct 27, 2018 Edward J. Los - Support the atlas calibration.
 */


#include <math.h>
#include <time.h>

#include "table.h"

#include "mysql.h"

#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

#include "scandb.h"
#include "pipelineutils.h"

#define NATIVE_OCTAVE 1
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


void
pix2wcs (

         struct WorldCoor *wcs,		/* World coordinate system structure */
         double	xpix, double ypix,	/* x and y image coordinates in pixels */
         double	*xpos,double *ypos);	/* RA and Dec in degrees (returned) */

double fd2jd (char *string);
double wcsdist(
               double x1,
               double y1,	/* (RA,Dec) or (Long,Lat) in degrees */
               double x2,
               double y2);	/* (RA,Dec) or (Long,Lat) in degrees */

double jd2lst (double dj);
void setlongitude (double longitude0);

/* #define MAX_CATALOG_EXPOSURES 33 */
#define MAX_BUFFER 512
#define MAX_COMMENT_STRING 132
#define MAX_CENTERSOURCE_STRING 10
#define MAX_RIGHTASCENSION_STRING 13
#define MAX_DECLINATION_STRING 16
#define MAX_DATE_STRING 25



typedef struct _starimage {
  /* Fields read in from the match table */
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  double ra;            /* Catalog Right Ascension in degrees */
  double dec;           /* Catalog Declination in degrees */
  double Stdmag;        /* Catalog magnitude */
  double color;         /* Catalog color */
  double MAG_ISO;       /* Sextractor isophotonic magnitude */
#if 0
  double MAG_APER;      /* Sextractor aperture magnitude */
  int ISOAREA_IMAGE; /* Sextractor isophotonic image area */
  double plate_dist;    /* Distance in degrees from plate center */
#endif
  double ddec;          /* declination deviation in arcsec */
  double dra;           /* ra deviation in arcsec */
  int NUMBER;           /* Sextractor reference number */
  int AFLAGS;            /* flags word */
  int BFLAGS;            /* flags word */
  double X_IMAGE;       /* Sextractor X location in pixels */
  double Y_IMAGE;       /* Sextractor Y location in pixels */
  /* Additional output fields */
  double X_WIDTH;        /* Width of the mosaic in pixels */
  double Y_WIDTH;        /* Height of the mosaic in pixels */
  double PLATE_SCALE;    /* Plate scale in arcsec */
  double extinction;     /* Extinction in magnitudes for this star */
#ifdef NATIVE_OCTAVE
  int selected;          /* If nonzero, selected for output */
#endif /* NATIVE_OCTAVE */
  int spatial_bin;       /* Annular bin index */
  int local_bin;         /* Smoothing bin index */
} STARIMAGE,*PSTARIMAGE;

typedef struct _catalog {
  /* Fields read in from the match table */
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  double median_local;        /* DASCH median magnitude */
} CATALOG,*PCATALOG;




int MatchCompare(const void *first, const void *second) 
{
  char *REFFirst = ((PSTARIMAGE)first)->REF;
  char *REFSecond = ((PSTARIMAGE)second)->REF;
  return(strcmp(REFFirst,REFSecond));
}

int CatalogCompare(const void *first, const void *second) 
{
  char *REFFirst = ((PCATALOG)first)->REF;
  char *REFSecond = ((PCATALOG)second)->REF;
  return(strcmp(REFFirst,REFSecond));
}




int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char infile[MAX_BUFFER];
  char cmdchar;
  int errorFlag = 0;
  char outfile[MAX_BUFFER];
  char fileroot[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  char extinctionfile[MAX_BUFFER];
  FILE *outHandle = NULL;
  FILE *extinctionHandle = NULL;

  File match_handle = NULL;
  char match_name[MAX_BUFFER];
  TableHead match_header = NULL;
  PSTARIMAGE match_table = NULL;
  size_t match_nrecs = 0;
  int match_index;
  PSTARIMAGE pMatch;

  File catalog_handle = NULL;
  char catalog_name[MAX_BUFFER];
  TableHead catalog_header = NULL;
  PCATALOG catalog_table = NULL;
  size_t catalog_nrecs = 0;
  int catalog_index;
  PCATALOG pCatalog;

  int verbose = 0;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int rotation;
  int binning;
  MYSQL my_connection;
  MOSAIC mosaicInfo;
  PMOSAIC pMosaic = &mosaicInfo;
  int exposureNumber;
  int solutionNumber = -1;
  double singleJulianDate;
  double singleTolerance;
  double allJulianDate;
  double allTolerance;
  int numExposures = 0;
#if 0
  PEXPOSURE pExposure;
  EXPOSURE exposure[MAX_CATALOG_EXPOSURES];
  double maxTolerance = 0;
  double minJulianDate;
  double maxJulianDate;

#endif

  double cd[4];
  struct WorldCoor *wcs;
  double xctr;
  double yctr;
  double cra;
  double cdec;
  double era;
  double edec;
  double dist;
  double pixdist;
  double scale;
  double xmin;
  double ymin;;
  int outCount = 0;
  int badFlagCount = 0;
  int keplerRejectCount = 0;
  int noColorCount = 0;
  int colorRejectCount = 0;
  int noMedianCount = 0;
  int lowAltitudeCount = 0;
  int cmpResult;
  int doExtinction;
  double dLST; /* LST in degrees */
  double minExtinction;
  double maxExtinction;
  double minAltitude;
  double maxAltitude;
  double extinction = 0;
  time_t startTime;
  time_t curTime;
  int maxREF = 0;
  int indexREF;
  char tmpREF[MAX_REF];
  int tmpMAGFlag;
  int readFlag = 0;
  int FitWCS;
  double hacorr = 0.0; /* normally zero, but set to force the center ha to zero when the timeAccuracy is too uncertain */

  double extinctionArray[X_EXTBINS_NORMAL][Y_EXTBINS_NORMAL]; /* Extinction array in magnitudes */
  double altitudeArray[X_EXTBINS_NORMAL][Y_EXTBINS_NORMAL];   /* Altitude array in degrees */
  int ix;
  int iy;
  int locationID;
  double latitude;
  double longitude;
  double elevation;
  int plateColor;
  double extinctionCoefficient;

  int catalogNumber;
  int nx = 0;
  int ny = 0;
  double edgeDist;
  int goodExtinctionCount = 0;
  int badExtinctionCount = 0;
  int allBelowHorizon = 0;
  int enx = 0;
  int eny = 0;

  time(&startTime);

  match_name[0] = 0;
  catalog_name[0] = 0;
  outfile[0] = 0;
  extinctionfile[0] = 0;
  fileroot[0] = 0;
  qualifier[0] = 0;
  catalogNumber =  CATALOG_GSC232;

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      printf("ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'e': /* extinction output file name */
        case 'E':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(extinctionfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
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


        case 'm': /* match file name */
        case 'M':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(match_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'p': /* file root */
        case 'P':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(fileroot,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'c': /* Pass 2 catalog file */
        case 'C':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(catalog_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 's': /* Solution Number */
        case 'S':
          argc--;
          nvals = sscanf(*++argv,"%d",&solutionNumber);
          if (nvals != 1) {
            printf("ERROR: Can not decode solutionNumber\n");
            errorFlag = 1;
          }
          break;


        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }

  if (outfile[0] == 0) {
    printf("ERROR: No output filename was specified\n");
    errorFlag = 1;
  }
  if (extinctionfile[0] == 0) {
    printf("ERROR: No extinction filename was specified\n");
    errorFlag = 1;
  }
  if (match_name[0] == 0) {
    printf("ERROR: No match filename was specified\n");
    errorFlag = 1;
  }
  if (fileroot[0] == 0) {
    printf("ERROR: No file root was specified\n");
    errorFlag = 1;
  }

  if ((qualifier[0] != 0) || (catalog_name[0] != 0)) {
    if (qualifier[0] == 0) {
      printf("ERROR: No qualifier was specified with an input catalog\n");
      errorFlag = 1;
    }
    if (catalog_name[0] == 0) {
      if (strstr(qualifier,"kepler") == NULL) {
        if (strstr(qualifier,"apass") == NULL) {
          if (strstr(qualifier,"gaia") == NULL) {
            if (strstr(qualifier,"atlas") == NULL) {
              if (strstr(qualifier,"experimental") == NULL) {
                printf("ERROR: A qualifier was specified without a catalog\n");
                errorFlag = 1;
              } else {
                catalogNumber = CATALOG_EXPERIMENTAL;
              }
            } else {
              catalogNumber = CATALOG_ATLAS;
            }
          } else {
            catalogNumber = CATALOG_GAIA;
          }
        } else {
          catalogNumber = CATALOG_APASS;
        }
      } else {
        catalogNumber = CATALOG_KEPLER;
      }
    } else {
      if (strstr(catalog_name,"kepler") != NULL) {
        catalogNumber = CATALOG_KEPLER;
      } else if (strstr(catalog_name,"apass") != NULL) {
        catalogNumber = CATALOG_APASS;
      } else if (strstr(catalog_name,"gaia") != NULL) {
        catalogNumber = CATALOG_GAIA;
      } else if (strstr(catalog_name,"atlas") != NULL) {
        catalogNumber = CATALOG_ATLAS;
      } else if (strstr(catalog_name,"experimental") != NULL) {
        catalogNumber = CATALOG_EXPERIMENTAL;
      }
    }

  }


  if (ParseFilename(fileroot,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
    printf("ERROR Failed to parse filename %s%s as %s%05d_%02d_%02d rotation %d\n",
           fileroot,qualifier,series,plateNumber,mosaicNumber,binning,rotation);
    errorFlag = 1;
  }


  /* Open the match file */
  match_handle = Open(match_name,"r");
  if (match_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the match file %s\n",match_name);
  } else {
    if (verbose) {
      printf("Found match file %s\n",match_name);
    }
  }

  if (strlen(catalog_name) > 0) {
    /* Open the catalog file */
    catalog_handle = Open(catalog_name,"r");
    if (catalog_handle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to find the catalog file %s\n",catalog_name);
    } else {
      if (verbose) {
        printf("Found catalog file %s\n",catalog_name);
      }
    }

  }

  outHandle = fopen(outfile,"wt");
  if (outHandle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to open the output file %s\n",outfile);
  } else {
    if (verbose) {
      printf("Output file %s\n",outfile);
    }
  }

  extinctionHandle = fopen(extinctionfile,"wt");
  if (extinctionHandle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to open extinction file %s\n",extinctionfile);
  }



  if (errorFlag) {
    printf("Usage: prepare_octave -p <file root> \n");
    printf("                      -m <match file> \n");
    printf("                      -o  <output file>\n");
    printf("                      -e  <extinction output file>\n");
    printf("                      -c  <pass 2 calibration file>\n");
    printf("                      -q  <pass 2 identifier>\n");
    printf("                      -s  <solution number>\n");
    printf("                      -v  verbose\n");

    return(-1);
  }

  doExtinction = 1; /* Assume we are caluclation extinction */
  dasch_init_scandb(&my_connection);

  /* Query the database for all the parameters we will need */
  if (GetSolutionJulianDate(&my_connection,series,plateNumber,mosaicNumber,solutionNumber,&exposureNumber,&numExposures,&singleJulianDate,&singleTolerance,&allJulianDate,&allTolerance)) {
    printf("ERROR: failed to get exposure information for %s%s\n",fileroot,qualifier);
    return(-1);
  }
    
#if 0 /* As of Nov 16, 2009, we will calculate the extinction anyways */
  if (singleTolerance > EXT_TIME_TOLERANCE) {
    printf("Extinction time tolerance %f days exceeds limit %f days for %s%s\n",singleTolerance,EXT_TIME_TOLERANCE,fileroot,qualifier);
    doExtinction = 0;
  }
#endif

  if (GetMosaicInfo(&my_connection,series,plateNumber,mosaicNumber,solutionNumber,pMosaic) != 1) {
    printf("ERROR: failed to get mosaic information for %s%s\n",fileroot,qualifier);
    return(-1);
  }

  if (pMosaic->rotation != rotation) {
    printf("ERROR: filename rotation %d does not agree with database rotation %d for %s%s\n",rotation,pMosaic->rotation,fileroot,qualifier);
    return(-1);
  }

  locationID = GetLocationId(&my_connection,series,plateNumber);

  if (locationID <= 0) {
    printf("ERROR: failed to get a location ID for %s%s\n",fileroot,qualifier);
    return(-1);
  } else {
    if (GetLongitude(&my_connection,locationID,&longitude,&latitude,&elevation,NULL,0) == 0) {
      printf("ERROR: failed to get a longitude for %s%s\n",fileroot,qualifier);
      return(-1);
    }
  }
  plateColor = GetPlateColor(&my_connection,series,plateNumber);
  extinctionCoefficient = GetExtinctionCoefficient(plateColor,elevation);

  /* Clear the StaleSolution bit */
  if (qualifier[0] == 0) {
    SetMosaicFitWCS(&my_connection,(char *)"StaleExposureNumber",fileroot,solutionNumber,0,readFlag,&FitWCS);
  }
  /* Done with the MySQL connection */
  mysql_close(&my_connection);





  /* Figure out the plate scale */
  if (strstr(pMosaic->ctype1,"DEC")) {
    char tmpPtr[MAX_CTYPE_STRING];
    double dtmp;
    int itmp;
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



  wcs = wcskinit(pMosaic->naxis1,
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
                 0);   /* epoch */

  nx = XDmagBins(pMosaic->naxis1,pMosaic->naxis2);
  ny = YDmagBins(pMosaic->naxis1,pMosaic->naxis2);
  enx = nx;
  eny = ny;

  xctr = 0.5 + (0.5 *pMosaic->naxis1);
  yctr = 0.5 + (0.5 *pMosaic->naxis2);
  pix2wcs(wcs,xctr,yctr,&cra,&cdec);
  xmin = 0.5;
  ymin = 0.5;
  pix2wcs(wcs,xmin,ymin,&era,&edec);

  dist = wcsdist(cra,cdec,era,edec);

  pixdist = sqrt(((xctr-xmin)*(xctr-xmin)) + ((yctr-ymin)*(yctr-ymin)));
  if (pixdist > 0.0) {
    scale = 3600.0*(dist/pixdist);
  } else {
    printf("ERROR, scale is unavailable for %s%s\n",fileroot,qualifier);
    return(-1);
  }



  /* Initialize the extinction array */
  for (ix = 0; ix < enx;ix++) {
    for (iy = 0; iy < eny;iy++) {
      extinctionArray[ix][iy] = 0.0;
      altitudeArray[ix][iy] = 0.0;
    }
  }
  minExtinction = 100000.00;
  maxExtinction = 0.00;
  minAltitude = 90.0;
  maxAltitude = -90.0;

  /* Calculate the extinction for all grid centers */
  if (doExtinction) {
    setlongitude(longitude);
    dLST = jd2lst(singleJulianDate) * 360.0/86400.0;
    /* Bugfix of Jan 15, 2015.  If the time is uncertain, set the scope at the meridian */
    if (singleTolerance > EXT_TIME_TOLERANCE) {
      hacorr = dLST-cra;
      printf("WARNING prepare_octave forcing hour angle to zero for filename %s%s_s%d \n",fileroot,qualifier,solutionNumber);
      
    }
    
    for (ix = 0; ix < enx;ix++) {
      for (iy = 0; iy < eny;iy++) {
        double xPixel;
        double yPixel;
        double dHA;
        double dRA;
        double rdecl;
        double z0;
        double x0;
        double z1;
        double altitude;
#if 0
        if ((ix == (enx/2)) && (iy == (eny/2))) {
          printf("At center\n");
        }
#endif
        xPixel = (0.5 + (1.0*ix)) * pMosaic->naxis1 / enx;
        yPixel = (0.5 + (1.0*iy)) * pMosaic->naxis2 / eny;
        pix2wcs(wcs,xPixel,yPixel,&dRA,&rdecl);

        dHA = dLST - dRA - hacorr;
        while (dHA < -180.0) {
          dHA = dHA + 360.;
        }
        while (dHA >= 180.) {
          dHA = dHA - 360.;
        }

#if 0
        printf("dHA %f for ix %d iy %d\n",dHA,ix,iy);
#endif


        /* cartesian coordinates in the frame of the celestial sphere */
        z0 = sin(rdecl *  DEGREES_TO_RAD);
        x0 = cos(rdecl * DEGREES_TO_RAD) * cos(dHA * DEGREES_TO_RAD);
        /* cartesian coordinates in the frame of the observer */
        z1 = (x0 * sin((90.0-latitude) * DEGREES_TO_RAD)) + (z0 * cos((90.0-latitude) * DEGREES_TO_RAD));
        altitude = RAD_TO_DEGREES * asin(z1);
        if (z1 > 0.0) {
          extinction = extinctionCoefficient/z1;
          if (1/z1 > MAX_AIRMASS) {
            badExtinctionCount++;
          } else {
            goodExtinctionCount++;
          }
          extinctionArray[ix][iy] = extinction; 
          if (extinction < minExtinction) {
            minExtinction = extinction;
          }
          if (extinction > maxExtinction) {
            maxExtinction = extinction;
          }
        } else {
          badExtinctionCount++;
          extinction = 0.0;
          extinctionArray[ix][iy] = extinction; 
        }
        altitudeArray[ix][iy] = altitude;
        if (altitude < minAltitude) {
          minAltitude = altitude;
        }
        if (altitude > maxAltitude) {
          maxAltitude = altitude;
        }

      }
    }
    printf("Extinction range is %f to %f (%f)  Altitude range is %f to %f (%f) for %s%s\n",
           minExtinction,maxExtinction,maxExtinction-minExtinction,
           minAltitude,maxAltitude,maxAltitude-minAltitude,
           fileroot,qualifier);
    if ((goodExtinctionCount+badExtinctionCount) != (enx*eny)) {
      printf("ERROR: prepare_octave extinction sanity check\n");
      exit(-1);
    }
    if (goodExtinctionCount == 0) {
      printf("WARNING: no points have an airmass greater than than MAX_AIRMASS %.2f %s%s\n",MAX_AIRMASS,fileroot,qualifier);
      allBelowHorizon = 1;
#if 0
      for (ix = 0; ix < enx;ix++) {
        for (iy = 0; iy < eny;iy++) {
          extinctionArray[ix][iy] = 0.0; 
        
        }
      }
#endif
    }

  }

  /* Now write out the extinction information to the extinction file 
   * Make it starbase compatible 
   */
  fprintf(extinctionHandle,"enx\teny\teix\teiy\taltitude\textinction\tDate\ttimeAccuracy\n");
  fprintf(extinctionHandle,"---\t---\t---\t---\t--------\t----------\t----\t------------\n");


  for (ix = 0; ix < enx;ix++) {
    for (iy = 0; iy < eny;iy++) {
      fprintf(extinctionHandle,"%d\t%d\t%d\t%d\t%.3f\t%.3f\t%f\t%f\n",
              enx,eny,ix,iy,altitudeArray[ix][iy],extinctionArray[ix][iy],singleJulianDate,singleTolerance);
    }
  }



  if (verbose) {
    printf("Read mosaic information for %s%s scale: %f arcsec/pixel\n",fileroot,qualifier,scale);
  }


  /* Now read in the match file */

  match_header = table_header(match_handle,TABLE_PARSE);
  if (match_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",match_name);
    return(-1);
  }

  match_table = table_loadva(match_handle,
                             &match_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(STARIMAGE),
                             &match_nrecs,
                             TblBuf,"REF"    ,TblOff(PSTARIMAGE,REF),MAX_REF,
                             TblDbl,"ra_2"   ,TblOff(PSTARIMAGE,ra),
                             TblDbl,"dec_2"  ,TblOff(PSTARIMAGE,dec),
                             TblDbl,"Stdmag",TblOff(PSTARIMAGE,Stdmag),
                             TblDbl,"color",TblOff(PSTARIMAGE,color),
                             TblDbl,"MAG_ISO",TblOff(PSTARIMAGE,MAG_ISO),
#if 0
                             TblDbl,"MAG_APER",TblOff(PSTARIMAGE,MAG_APER),
                             TblInt,"ISOAREA_IMAGE",TblOff(PSTARIMAGE,ISOAREA_IMAGE),
                             TblDbl,"plate_dist",TblOff(PSTARIMAGE,plate_dist),
#endif
                             TblDbl,"ddec",TblOff(PSTARIMAGE,ddec),
                             TblDbl,"dra",TblOff(PSTARIMAGE,dra),
                             TblInt,"NUMBER" ,TblOff(PSTARIMAGE,NUMBER),
                             TblInt,"AFLAGS"  ,TblOff(PSTARIMAGE,AFLAGS),
                             TblInt,"BFLAGS"  ,TblOff(PSTARIMAGE,BFLAGS),
                             TblDbl,"X_IMAGE",TblOff(PSTARIMAGE,X_IMAGE),
                             TblDbl,"Y_IMAGE",TblOff(PSTARIMAGE,Y_IMAGE),
                             0,"end",0);
  if (match_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",match_name);
    return(-1);
  }
  if (verbose) {
    printf("read %d records for %s\n",match_nrecs,match_name);
  }



  /* Now read in the catalog file */
  if (catalog_handle != NULL) {
    catalog_header = table_header(catalog_handle,TABLE_PARSE);
    if (catalog_header == NULL) {
      printf("ERROR: Failed to read header for %s\n",catalog_name);
      return(-1);
    }

    catalog_table = table_loadva(catalog_handle,
                                 &catalog_header,
                                 NULL, /* hbase */
                                 NULL, /* rows */
                                 NULL,
                                 sizeof(CATALOG),
                                 &catalog_nrecs,
                                 TblBuf,"REF"    ,TblOff(PCATALOG,REF),MAX_REF,
                                 TblDbl,"median_local"   ,TblOff(PCATALOG,median_local),
                                 0,"end",0);
    if (catalog_table == NULL) {
      printf("ERROR: Failed to read table for %s\n",catalog_name);
      return(-1);
    }
    if (verbose) {
      printf("read %d records for %s\n",catalog_nrecs,catalog_name);
    }

    /* Now sort both tables by REF */
    if (verbose) {
      time(&curTime);
      curTime -= startTime;
      printf("Begin catalog sort at %d seconds\n",curTime);
    }
    qsort((void*)catalog_table,catalog_nrecs,sizeof(CATALOG),CatalogCompare);
    if (verbose) {
      time(&curTime);
      curTime -= startTime;
      printf("Begin match table sort at %d seconds\n",curTime);
    }
    qsort((void*)match_table,match_nrecs,sizeof(STARIMAGE),MatchCompare);
    if (verbose) {
      time(&curTime);
      curTime -= startTime;
      printf("End match table sort at %d seconds\n",curTime);
    }
    catalog_index = 0;
    pCatalog = &catalog_table[catalog_index];

  }



  for (match_index = 0; match_index < match_nrecs; match_index++) {
    pMatch = &match_table[match_index];
#if 0
    if (pMatch->NUMBER == 100) {
      printf("At NUMBER %d\n",pMatch->NUMBER);
    }
#endif
    pMatch->selected = 0;
    if (catalog_handle != NULL) {
      while ((cmpResult = strcmp(pMatch->REF,pCatalog->REF)) > 0) {
        catalog_index++;
        if (catalog_index >= catalog_nrecs) {
          /* Can not find a match */
          break;
        }
        pCatalog = &catalog_table[catalog_index];
      }
      if (cmpResult == 0) {
        pMatch->Stdmag = pCatalog->median_local;
      } else {
        noMedianCount++;
        continue;
      }
    } 

    /* Find the extinction for this file */
    ix = (pMatch->X_IMAGE * enx)/pMosaic->naxis1;
    iy = (pMatch->Y_IMAGE * eny)/pMosaic->naxis2;
    if (ix < 0) {
      ix = 0;
    }
    if (ix >= enx) {
      ix = enx-1;
    }
    if (iy < 0) {
      iy = 0;
    }
    if (iy >= eny) {
      iy = eny-1;
    }
    /* We add extinction to the GSC2.3.2 magnitudes to make them dimmer */
    extinction = extinctionArray[ix][iy];

#if 0
    if (doExtinction == 0) {
      pMatch->BFLAGS |= (1 << FILTER_BFLAG_BAD_JD);
    }
#endif
#if 0
    printf("AFLAGS %d NUMBER %d\n",pMatch->AFLAGS,pMatch->NUMBER);
#endif
    if ((doExtinction != 0) && ((extinction == 0.0) || (extinction > (extinctionCoefficient * MAX_AIRMASS)))) {
      if (allBelowHorizon == 0) {
        pMatch->AFLAGS |= (1 << FILTER_AFLAG_LOW_ALTITUDE);
      }
      lowAltitudeCount++;
    }

    if (pMatch->AFLAGS >= FILTER_AFLAG_BAD) {
      
      badFlagCount++;
      continue;
    }
  
    /* For the Kepler Input Catalog, toss all but SCP and Tycho-2 values. 
       Set the MAGFlag to zero to avoid confusion with color corrections */
    if ((pMatch->BFLAGS & (1<<FILTER_BFLAG_KEPLER)) != 0) {
      tmpMAGFlag = pMatch->BFLAGS >> GSC_MAGNITUDE_FLAG_BIT;
      tmpMAGFlag &= MAGNITUDE_MASK;
      if ((tmpMAGFlag != KEPLER_CQ_SCP) &&
          (tmpMAGFlag != KEPLER_CQ_TYBV)) {
        keplerRejectCount++;
        continue;
      } else {
        pMatch->BFLAGS &= ~(MAGNITUDE_MASK <<  GSC_MAGNITUDE_FLAG_BIT);
      }
	
    }

    /* Reject all stars without a color or outside the color limits */
    if (pMatch->color > 90.0) {
      noColorCount++;
      continue;
    } else {
      
      switch (catalogNumber) {
      case CATALOG_GSC232:
      case CATALOG_EXPERIMENTAL:
        if ((pMatch->color < MIN_GSC_COLOR) ||
            (pMatch->color > MAX_GSC_COLOR)) {
          colorRejectCount++;
          continue;
        }
        break;
      case CATALOG_KEPLER:
        if ((pMatch->color < MIN_KEPLER_COLOR) ||
            (pMatch->color > MAX_KEPLER_COLOR)) {
          colorRejectCount++;
          continue;
        }
        break;
      case CATALOG_APASS:
        if ((pMatch->color < MIN_APASS_COLOR) ||
            (pMatch->color > MAX_APASS_COLOR)) {
          colorRejectCount++;
          continue;
        }
        break;
      case CATALOG_GAIA:
        if ((pMatch->color < MIN_GAIA_COLOR) ||
            (pMatch->color > MAX_GAIA_COLOR)) {
          colorRejectCount++;
          continue;
        }
        break;
      case CATALOG_ATLAS:
        if ((pMatch->color < MIN_ATLAS_COLOR) ||
            (pMatch->color > MAX_ATLAS_COLOR)) {
          colorRejectCount++;
          continue;
        }
        break;
      default:
        printf("ERROR: Unknown catalog in prepare_octave\n");
        exit(-1);
        break;
      }
    }
    /* Calculate the spatial bin index and the local bin index here */
#if 0
    if (pMatch->NUMBER == 162463) {
      printf("At NUMBER %d\n",pMatch->NUMBER);
    }
#endif
    pMatch->spatial_bin = CalculateBin(pMosaic->naxis1,pMosaic->naxis2,pMatch->X_IMAGE,pMatch->Y_IMAGE,&edgeDist);
    ix = (pMatch->X_IMAGE * nx) /(1.0 * pMosaic->naxis1);
    iy = (pMatch->Y_IMAGE * ny) /(1.0 * pMosaic->naxis2);
    pMatch->local_bin = ix + (nx*iy);

    outCount++;
#ifdef NATIVE_OCTAVE
    pMatch->selected = 1;
    pMatch->extinction = extinction;
    if (strlen(pMatch->REF) > maxREF) {
      maxREF = strlen(pMatch->REF);
    }

#else /* NATIVE_OCTAVE */
    fprintf(outHandle,"%s\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%f\n",
            pMatch->REF,
            pMatch->ra,
            pMatch->dec,
            pMatch->Stdmag+extinction,
            pMatch->color,
            pMatch->MAG_ISO,
            pMatch->NUMBER,
            pMatch->BFLAGS,
            pMatch->X_IMAGE,
            pMatch->Y_IMAGE,
            pMatch->ddec,
            pMatch->dra,
            extinction,
            pMatch->spatial_bin,
            pMatch->local_bin,
            pMosaic->naxis1,
            pMosaic->naxis2,
            scale);
#endif /* NATIVE_OCTAVE */
  }

#ifdef NATIVE_OCTAVE
  fprintf(outHandle,"# Created by Octave 2.9.8, Mon Mar 31 19:59:22 2008 EDT <scanner@localhost.localdomain>\n");
  fprintf(outHandle,"# name: inputREF\n");
  fprintf(outHandle,"# type: sq_string\n");
  fprintf(outHandle,"# elements: %d\n",outCount);
  if (maxREF > MAX_REF) {
    printf("ERROR: maxREF %d exceeds MAX_REF %d\n",maxREF,MAX_REF);
    return(-1);
  }
  for (match_index = 0; match_index < match_nrecs; match_index++) {
    pMatch = &match_table[match_index];
#if 0
    if (pMatch->NUMBER == 100) {
      printf("At NUMBER %d\n",pMatch->NUMBER);
    }
#endif
    if (pMatch->selected) {
      fprintf(outHandle,"# length: %d\n",maxREF);
      strcpy(tmpREF,pMatch->REF);
      for (indexREF = strlen(pMatch->REF); indexREF < MAX_REF; indexREF++) {
        tmpREF[indexREF] = ' ';
      }
      tmpREF[maxREF] = 0;
      fprintf(outHandle,"%s\n",tmpREF);
    }
  }
  fprintf(outHandle,"# name: inputx\n");
  fprintf(outHandle,"# type: matrix\n");
  fprintf(outHandle,"# rows: %d\n",outCount);
  fprintf(outHandle,"# columns: 18\n");
  for (match_index = 0; match_index < match_nrecs; match_index++) {
    pMatch = &match_table[match_index];
#if 1
    if (pMatch->NUMBER == 100) {
      printf("At NUMBER %d\n",pMatch->NUMBER);
    }
#endif
    if (pMatch->selected) {
      fprintf(outHandle," 0 %f %f %f %f %f %d %d %f %f %f %f %f %d %d %d %d %f\n",
              pMatch->ra,
              pMatch->dec,
              pMatch->Stdmag+pMatch->extinction,
              pMatch->color,
              pMatch->MAG_ISO,
              pMatch->NUMBER,
              pMatch->BFLAGS,
              pMatch->X_IMAGE,
              pMatch->Y_IMAGE,
              pMatch->ddec,
              pMatch->dra,
              pMatch->extinction,
              pMatch->spatial_bin,
              pMatch->local_bin,
              pMosaic->naxis1,
              pMosaic->naxis2,
              scale);
    }
  }


#endif /* NATIVE_OCTAVE */


  /* All done.  Clean up */


  if (match_header != NULL) {
    table_hdrfree(match_header);
  }
  if (match_handle != NULL) {
    Close(match_handle);
  }

  if (catalog_header != NULL) {
    table_hdrfree(catalog_header);
  }
  if (catalog_handle != NULL) {
    Close(catalog_handle);
  }


  if (outHandle != NULL) {
    fclose(outHandle);
  }
  if (extinctionHandle != NULL) {
    fclose(extinctionHandle);
  }
  if (pRangeTable != NULL) {
    free(pRangeTable);
  }
  if (pLocationTable != NULL) {
    free(pLocationTable);
  }


  time(&curTime);
  curTime -= startTime;
  printf("prepare_octave wrote %d records.  Bad flags %d No median %d too low %d noColor %d colorReject %d in %d seconds for %s%s\n",
         outCount,
         badFlagCount,
         noMedianCount,
         lowAltitudeCount,
         noColorCount,
         colorRejectCount,
         curTime,
         fileroot,
         qualifier);
    


  return(0);
}
