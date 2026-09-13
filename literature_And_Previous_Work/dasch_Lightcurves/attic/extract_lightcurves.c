// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* extract_lightcurves.c
 *
 *  Search the master output database for individual stars, and for each star found, output lc_<src_name>_<qualifier>.db
 *  with the following columns:
 *   1: sequence
 *   2: heliocentric Julian Date
 *   3: spatial_bin (replaces plate_dist)
 *   4: flag   17   Solid Circle              for a good point
 *             17.5 Small solid circle        for a good point, bin 9
 *             22   Open Circle               for a large error point
 *             22.5 Small open circle         for a large error point, bin 9
 *             31   Arrow pointing down       for not found
 *             31.5 Small arrow pointing down for not found, bin 9
 *              3   Asterisk                  for a plate defect
 *              7   Triangle                  for a high drad object
 *              6   Open Square               for a Pickering Wedge object
 *             16   Solid Square              for a high extinction
 *             30   Arrow pointing up         for a rejected blend object or a too bright object
 *             11.5 Diamond                   for a blended object
 *             10.5 Cusped Square             for a high local correction object.
 *             99.5 Letter "c"                for a good object without color correction
 *   5: magcal_local
 *   6: magcal_local_rms
 *   7: nlocal
 *   8: Plate
 *   9: magcal_iso_rms
 *  10: magcal_local_error
 *  11: magcal_iso
 *  12: limiting_mag_local
 *
 *  Also write a file lc_<src_name>_<qualifier>.txt with the following values for use in plotting:
 *   1.  clipped median magcal_local
 *   2.  clipped magcal_local rms
 *   3.  total number of points (including non-detections)
 *   4.  number of good points, not in bin 9 and magcal_local_rms < 1.0
 *   5.  begining year
 *   6.  ending year
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include -I/usr/include/mysql -L /dasch/install/lib -lm  -L/usr/lib64/mysql -lmysqlclient extract_lightcurves.c pipelineutils.a -ltable -lutil -lwcs -o extract_lightcurves 
 *
 * /dasch/Pipeline/extract_lightcurves -i -l /dasch/Pipeline/3c273_total.list -t /dasch/Pipeline/3c273_REF.db -d /dasch/Pipeline/lightcurves/3c273
 *
 *
 * /dasch/Pipeline/extract_lightcurves -g -v  -l /dasch/Pipeline/3c273_total.list -t /dasch/Pipeline/3c273_REF.db  -d /dasch/Pipeline/lightcurves/3c273 -b N120013341 -b N120013339
 *
 *  To extract everything from M44 
 *  echo "extract_lightcurves -i -g -v   -f -l /dasch/Pipeline/M44_total.list  -a 3 -d /dasch/junk/M44_all" | at now
 *  echo "extract_lightcurves -i -g -v   -f -l /dasch/Pipeline/M44_total.list  -a 100 -d /dasch/raid011/Backup/M44_all" | at now
 * 
 * Mar 11, 2008 Edward J. Los - Initial version
 * Apr 24, 2008 Edward J. Los - Revice CalcMedianAndRMS return test
 * May 29, 2008 Edward J. Los - Remove high drad and Pickering wedge objects.
 * Jun 10, 2008 Edward J. Los - Remove plate defects
 * Jun 13, 2008 Edward J. Los - add magcal_local_error and magcal_iso_rms
 * Jun 18, 2008 Edward J. Los - restore magcal_iso
 * Jun 20, 2008 Edward J. Los - correct length of src_name on table read
 * Jul  1, 2008 Edward J. Los - add limiting_mag_local.
 *                              Plot blended stars as diamonds.
 *                              Plot unmatched blended stars as up arrows at the limiting magnitude
 * Jul  8, 2008 Edward J. Los - Plot stars with high zout as a square with concave sides.
 * Jul 15, 2008 Edward J. Los - Add "-s" and "-c" qualifiers to test rms implementations
 * Jul 28, 2008 Edward J. Los - Do not plot undetected points where local calibration has failed.
 *                              Plot too bright stars as upward pointing arrows.
 *                              Use smaller diamonds for blended stars, and high zout stars
 * Aug 15, 2008 Edward J. Los - Restore blended objects to favor 
 *                              Do not consider an image within 0.5 mag of the limiting magnitude a good point
 * Aug 18, 2008 Edward J. Los - Reject blended stars that fail any other tests, except possibly the defect test
 *                              Add IGNORE_LIMITING for home debugging.
 * Aug 19, 2008 Edward J. Los - Add FLAGS and reject_region to the output file
 *                              Correct plot for objects close to the limiting magnitude
 * Aug 22, 2008 Edward J. Los - Enable ADD_DEFECT_TO_BLEND to ignore the defect flag on blended objects.
 *                              Add error_bar_factor, the factor by which magcal_local_rms and magcal_iso_rms has been
 *                              reduced so that the clipped median RMS is equal to the clipped zero-based RMS of 
 *                              magcal_local_rms for good stars.
 * Aug 25, 2008 Edward J. Los - Correct omission to plot TOO_BRIGHT objects correctly.
 *                              Change test for rms1 statsticstics
 * Sep 13, 2008 Edward J. Los - Add capability to merge multiple stars in a single graph to handle blended images
 * Sep 15, 2008 Edward J. Los - Reject case C and D blended stars
 *                              If a merge companion can not be found, use the DASCH median magnitude.  Else use
 *                              Stdmag
 *                              Correct calculation of dasch median rms by excluding stars close to the limiting magnitude.
 * Oct  7, 2008 Edward J. Los - Make a copy of the FLAGS integer and if it is a blended object, clear the FILTER_FLAG_DRAD
 *                              and FILTER_FLAG_DEFECT bits in the original FLAGS entry.  However, always display
 *                              the original FLAGS integer
 *                              Correct magcal_iso and limiting_mag_local for objects not found.
 * Oct 10, 2008 Edward J. Los - Add OFFPLATE and BINFAILURE reject reasons, create summaries of rejected and undetected objects.
 * Oct 14, 2008 Edward J. Los - Add -f to ignore the defect flag
 * Oct 23, 2008 Edward J. Los - Plot upper limit instead of limiting magnitude for undetected stars
 * Nov 13, 2008 Edward J. Los - Correct IGNORE_LIMITING bug in finding the reject reason
 * Nov 23, 2008 Edward J. Los - Note whether undetected stars lie in a failed spatial bin
 * Dec  8, 2008 Edward J. Los - Move reject reason to pipelineutils.
 * Dec 15, 2008 Edward J. Los - Correct erroneous REJECT_REASON_HIZOUT reject reason flag
 * Dec 22, 2008 Edward J. Los - Add reject reason for low altitude (high extinction objects)
 * Jan 28, 2009 Edward J. Los - Split FLAGS into AFLAGS and BFLAGS.
 * Feb  6, 2009 Edward J. Los - Replace FILTER_BFLAG_BLEND with FILTER_AFLAG_BLEND
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 * Feb 23, 2009 Edward J. Los - Remove rms statistics code
 * Feb 24, 2009 Edward J. Los - Use common DecodeAFLAGS routine
 *                              Add dradRMS2 to the output files.
 * Feb 25, 2009 Edward J. Los - make duplicate GSC2.3.2 ID's for a plate a soft error.
 * Apr  4, 2009 Edward J. Los - Flag objects without colorterm correction in plots
 * May  4, 2009 Edward J. Los - Make IGNORE_LIMITING an input parameter
 * Jun 12, 2009 Edward J. Los - Add the -a qualifier to extract all lightcurves.
 * Jun 16, 2009 Edward J. Los - Add GSC catalog data to the bottom of the plot.  Use the gsc_bin_index
 *                              to locate the catalog information
 * Jun 22, 2009 Edward J. Los - Add startJD and endJD to the text file to make the axis agree with the start and end ephemeris year.
 * Jun 29, 2009 Edward J. Los - Add 0.1% margins to ensure that all points get plotted
 * Aug 20, 2009 Edward J. Los - Add Kepler Input Catalog Support
 * Mar 13, 2011 Edward J. Los - Add apass catalog support
 * Jul 30, 2012 Edward J. Los - Add experimental catalog support
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * May 29, 2018 Edward J. Los - Add gaia support
 * Oct 28, 2018 Edward J. Los - Add atlas refcat2 support
 *
 * WARNING: This file is obsolete and will not handle multiple exposures correctly
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "searchgsc.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "daschunistd.h"
#define MAX_PLATE_NAME 32
#define MAX_BUFFER 100
#define MARGIN_WIDTH 0.001
/* #define ACCEPT_LOW_ALTITUDE */


extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;


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

double jd2hjd (double	dj,	/* Julian date (geocentric) */
               double	ra,	/* Right ascension (degrees) */
               double	dec,	/* Declination (degrees) */
               int	sys);	/* J2000, B1950, GALACTIC, ECLIPTIC */

void
wcs2pix (
         struct WorldCoor *wcs,	/* World coordinate system structure */
         double	xpos,double ypos,	/* World coordinates in degrees */
         double	*xpix,double *ypix,	/* Image coordinates in pixels */
         int	*offscl);	/* 0 if within bounds, else off scale */

void
dec2str (

         char	*string,	/* Character string (returned) */
         int	lstr,		/* Maximum number of characters in string */
         double	dec,		/* Declination in degrees */
         int	ndec);		/* Number of decimal places in arcseconds */

void ra2str (

             char	*string,	/* Character string (returned) */
             int	lstr,		/* Maximum number of characters in string */
             double	ra,		/* Right ascension in degrees */
             int	ndec);		/* Number of decimal places in seconds */

typedef struct _spatial_bin {
  int spatial_bin;
  double limiting_mag;
  double max_bright_mag;
  double limiting_iso;
  double max_bright_iso;
  double upper_limit;
} SPATIAL_BIN,*PSPATIAL_BIN;

typedef struct _dmagcor {
  int nx;             /* Total bins in width */
  int ny;             /* Total bins in height */
  int ix;             /* width bin */
  int iy;             /* height bin */
  double zout;        /* Magnitude correction */
  double errout;      /* RMS of magnitude correction */
  int npout;          /* Number of points used for correction */
  int rejectFlag;     /* Reasons for rejecting this bin */ 
  double dradRMS2;    /* drad RMS for this bin */
} DMAGCOR,*PDMAGCOR;

typedef struct _extinction {
  int enx;             /* Total bins in width */
  int eny;             /* Total bins in height */
  int eix;             /* width bin */
  int eiy;             /* height bin */
  double altitude;    /* altitude in degrees */
  double extinction;  /* extinction in magnitudes */


} EXTINCTION,*PEXTINCTION;

typedef struct _target {
  char REF[MAX_REF*MAX_BLEND_COUNT];   /* GSC2.3.2 reference number */
  char src_name[MAX_SRC_LENGTH*MAX_BLEND_COUNT];   /* GSC2.3.2 reference number */
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
  /* The following data is extracted from the GSC catalog */
  int gsc_bin_index;
  double Stdmag;
  double color;
  char class;
  char VFlag;
  char MAGFlag;
  double RaPM;
  double DecPM;

} TARGET,*PTARGET;

typedef struct _master {
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  double magcal_local;
  double Date;
  double magcal_iso;
  double magcal_iso_rms;
  double magcal_local_error;
  double magcal_local_rms;
  double limiting_mag_local;
  double Stdmag;
  double ra;
  double dec;
  int npoints_local;
  int spatial_bin;
  int AFLAGS;
  int AFLAGSCOPY;
  int BFLAGS;
  int local_bin_index;
  int gsc_bin_index;
  char Plate[MAX_PLATE_NAME];
} MASTER,*PMASTER;


int binSearchCount = 0;
int cacheHits = 0;
int cacheMisses = 0;

typedef struct _plate {
  char Plate[MAX_PLATE_NAME];
  char PlateStub[MAX_PLATE_NAME];
  int errorFlag;
  double geoJulianDate;
  EXPOSURE exposure;
  MOSAIC mosaic;
  struct WorldCoor *wcs;
  PSPATIAL_BIN spatial_bin_table;
  int nx;                   /* dmagcor_table width */
  int ny;                   /* dmagcor_table height */
  PDMAGCOR dmagcor_table;
  int enx;                   /* extinction width */
  int eny;                   /* extinction height */
  PEXTINCTION extinction_table;
  int starDetected;     /* If non-zero, this object has been detected */
  int validPoint;     /* If non-zero, this point is valid */
  /* Information copied for each point */
  double magcal_local; /* Local magnitude, or local limiting magnitude */
  double Date;  /* Heliocentric Julian Date */
  double magcal_local_rms;
  double magcal_local_error;
  double magcal_iso;
  double magcal_iso_rms;
  double limiting_mag_local;
  double Stdmag;
  int spatial_bin;
  int npoints_local;
  int AFLAGS;
  int AFLAGSCOPY;
  int BFLAGS;
  /* Information calculated for output */
  double dradRMS2;    /* drad RMS for this bin */
  double clipmed;
  double flag;     /* wip flag symbol (see comment at the top of this file) */
  int reject_reason1;  /* Total  Reject reason mask.  See pipelineutils.h for defintions */
  int reject_reason2;  /* Actual Reject reason mask.  See pipelineutils.h for defintions */
  int gsc_bin_index;
  PTARGET pTarget;
} PLATE,*PPLATE;

static int ndec = 3;		/* Number of decimal places in non-angles */
int ReadMaster(PMASTER pMaster,
               File master_handle,
               TableHead master_header,
               TblDescriptor master_descriptor,
               TableRow* master_row)
{
  *master_row = table_rowget(master_handle,master_header,*master_row,NULL,NULL,0);
  if (*master_row == NULL) {
    return(0);
  }
  if (!table_loadrow(master_handle,master_header,*master_row,master_descriptor,(char *)pMaster)) {
    printf("ERROR: Read Master table_loadrow failed\n");
    return(0);
  }
  pMaster->AFLAGSCOPY = pMaster->AFLAGS;
  return(1);
}

/* Sort routine to put the plates in Julian Day order */
int JulianDayCompare(const void *first, const void *second) 
{
  double numberFirst = ((PPLATE)first)->geoJulianDate;
  double numberSecond = ((PPLATE)second)->geoJulianDate;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}


int InitPlateList(PPLATE plateArray,int plateListLength,char *ingestDirectory,char* catalogallDirectory,char *qualifier,int verbose,int ignoreLimiting) 
{
  MYSQL my_connection;
  char *mysqlhost;
  char *username;
  char *password;
  PPLATE pPlate;
  int plateIndex;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int index;
  int irec;
  int ix;
  int iy;
  int eix;
  int eiy;
  char *underscorePtr;
  int mosaicNumber;
  int binning;
  int rotation;
  double cd[4];

  File spatial_bins_handle = NULL;
  char spatial_bins_name[MAX_BUFFER];
  PSPATIAL_BIN tmp_spatial_bin_table = NULL;
  TableHead spatial_bins_header = NULL;
  PSPATIAL_BIN pSpatial_bin;
  int spatial_bin;
  size_t spatial_bin_nrecs;

  File dmagcor_handle = NULL;
  char dmagcor_name[MAX_BUFFER];
  TableHead dmagcor_header = NULL;
  PDMAGCOR pDmagcor;
  size_t dmagcor_nrecs;

  File extinction_handle = NULL;
  char extinction_name[MAX_BUFFER];
  TableHead extinction_header = NULL;
  size_t extinction_nrecs = 0;
  int extinction_index;
  PEXTINCTION pExtinction;

  /* Connect to the database */
  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    fprintf(stderr,"DASCH_MYSQLHOST is not defined\n");
    return(-1);
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


  if (!mysql_real_connect(&my_connection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(&my_connection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(&my_connection),mysql_error(&my_connection));
    }
    return(-1);
  }

  for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {


    pPlate = &plateArray[plateIndex];
    pPlate->errorFlag = 0;
    pPlate->starDetected = 0;
    pPlate->validPoint = 0;
    pPlate->reject_reason1 = 0;
    pPlate->reject_reason2 = 0;
    if ((verbose > 0) && (((plateIndex+1) % 100) == 0)) {
      fprintf(stderr,"Reading plate calibration files for plate %5d %s\n",plateIndex+1,pPlate->PlateStub);
    }
    strcpy(pPlate->PlateStub,pPlate->Plate);
    underscorePtr = strstr(pPlate->PlateStub,"_");
    if (underscorePtr != NULL) {
      *underscorePtr = 0;
    }


    if (ParseFilename(pPlate->Plate,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
      fprintf(stderr,"ERROR: Failed to parse plate name %s\n",pPlate->Plate);
      pPlate->errorFlag = 1;
      continue;
    }
#if 0
    if ((strcmp(series,"dsb") == 0) &&
        (plateNumber == 2822)) {
      printf("At %s%05d\n",series,plateNumber);
    }
#endif
    if(GetExposureInfo(&my_connection,
                       series,
                       plateNumber,
                       0,
                       &pPlate->exposure,
                       0)) {
      pPlate->geoJulianDate = fd2jd(pPlate->exposure.date);
    } else {
      fprintf(stderr,"Failed to get exposure record for %s\n",pPlate->Plate);
      pPlate->errorFlag = 1;
      continue;
    }

    if (GetMosaicInfo(&my_connection,series,plateNumber,mosaicNumber,0,&pPlate->mosaic) != 1) {
      fprintf(stderr,"ERROR: failed to get mosaic information for %s\n",pPlate->Plate);
      pPlate->errorFlag = 1;
      continue;
    }
    /* Initialize the wcs structure for this plate */
    if (strstr(pPlate->mosaic.ctype1,"DEC")) {
      char tmpPtr[MAX_CTYPE_STRING];
      double dtmp;
      int itmp;
      strcpy(tmpPtr,pPlate->mosaic.ctype2);
      strcpy(pPlate->mosaic.ctype2,pPlate->mosaic.ctype1);
      strcpy(pPlate->mosaic.ctype1,tmpPtr);
      dtmp = pPlate->mosaic.crval1;
      pPlate->mosaic.crval1 = pPlate->mosaic.crval2;
      pPlate->mosaic.crval2 = dtmp;
            


      cd[0] = pPlate->mosaic.cd2_1;
      cd[1] = pPlate->mosaic.cd2_2;
      cd[2] = pPlate->mosaic.cd1_1;
      cd[3] = pPlate->mosaic.cd1_2;

    } else {

      cd[0] = pPlate->mosaic.cd1_1;
      cd[1] = pPlate->mosaic.cd1_2;
      cd[2] = pPlate->mosaic.cd2_1;
      cd[3] = pPlate->mosaic.cd2_2;

    }



    pPlate->wcs = wcskinit(pPlate->mosaic.naxis1,
                           pPlate->mosaic.naxis2,
                           pPlate->mosaic.ctype1,
                           pPlate->mosaic.ctype2,
                           pPlate->mosaic.crpix1,
                           pPlate->mosaic.crpix2,
                           pPlate->mosaic.crval1,
                           pPlate->mosaic.crval2,
                           cd,
                           0,  /* cdelt1 */
                           0,  /* cdelt2 */
                           0,  /* crota */
                           2000, /* equinox */
                           0);   /* epoch */


    if (ignoreLimiting == 0) {
      /* Now read in the spatial bin file */

      strcpy(spatial_bins_name,ingestDirectory);
      strcat(spatial_bins_name,"/");
      strcat(spatial_bins_name,pPlate->Plate);
      strcat(spatial_bins_name,qualifier);
      strcat(spatial_bins_name,".out.spatial_bins.db");
      spatial_bins_handle = Open(spatial_bins_name,"r");
      if (spatial_bins_handle == NULL) {
        fprintf(stderr,"ERROR: Failed to find spatial bins file %s\n",spatial_bins_name);
        pPlate->errorFlag = 1;
        continue;
      } else {
        if (verbose > 1) {
          fprintf(stderr,"Found spatial bins file %s\n",spatial_bins_name);
        }
      }

      /* Now read in the spatial bins file */
      spatial_bins_header = table_header(spatial_bins_handle,TABLE_PARSE);
      if (spatial_bins_header == NULL) {
        fprintf(stderr,"ERROR: Failed to read header for %s\n",spatial_bins_name);
        pPlate->errorFlag = 1;
        continue;
      }

      tmp_spatial_bin_table = table_loadva(spatial_bins_handle,
                                           &spatial_bins_header,
                                           NULL, /* hbase */
                                           NULL, /* rows */
                                           NULL,
                                           sizeof(SPATIAL_BIN),
                                           &spatial_bin_nrecs,
                                           TblInt,"spatial_bin",TblOff(PSPATIAL_BIN,spatial_bin),
                                           TblDbl,"limiting_mag",TblOff(PSPATIAL_BIN,limiting_mag),
                                           TblDbl,"max_bright_mag",TblOff(PSPATIAL_BIN,max_bright_mag),
                                           TblDbl,"limiting_iso",TblOff(PSPATIAL_BIN,limiting_iso),
                                           TblDbl,"max_bright_iso",TblOff(PSPATIAL_BIN,max_bright_iso),
                                           TblDbl,"upper_limit",TblOff(PSPATIAL_BIN,upper_limit),
                                           0,"end",0);
      if (tmp_spatial_bin_table == NULL) {
        fprintf(stderr,"ERROR: Failed to read table for %s\n",spatial_bins_name);
        pPlate->errorFlag = 1;
        continue;
      }
      if (verbose > 1) {
        fprintf(stderr,"read %d records for %s\n",spatial_bin_nrecs,spatial_bins_name);
      }
      /* Now copy the spatial bin entries into their proper slots */
      pPlate->spatial_bin_table = (PSPATIAL_BIN)calloc(MAX_SPATIAL_BINS+1,sizeof(SPATIAL_BIN));
      if (pPlate->spatial_bin_table == NULL) {
        fprintf(stderr,"ERROR: Failed to allocate the spatial bin table \n");
        return(-1);
      }
      for (index = 0; index < spatial_bin_nrecs; index++) {
        spatial_bin = tmp_spatial_bin_table[index].spatial_bin;
        if ((spatial_bin > 0) && (spatial_bin <= MAX_SPATIAL_BINS)) {
          memcpy(&pPlate->spatial_bin_table[spatial_bin],&tmp_spatial_bin_table[index],sizeof(SPATIAL_BIN));
        } else {
          fprintf(stderr,"ERROR: Illegal spatial bin %d\n",spatial_bin);
        }
      }
         
      if (tmp_spatial_bin_table != NULL) {
        Free(tmp_spatial_bin_table);
      }
  
      if (spatial_bins_header != NULL) {
        table_hdrfree(spatial_bins_header);
      }
      if (spatial_bins_handle != NULL) {
        Close(spatial_bins_handle);
      }

      /* Open the local calibration file */
      strcpy(dmagcor_name,ingestDirectory);
      strcat(dmagcor_name,"/");
      strcat(dmagcor_name,pPlate->Plate);
      strcat(dmagcor_name,qualifier);
      strcat(dmagcor_name,"_dmagcor.db");
      dmagcor_handle = Open(dmagcor_name,"r");
      if (dmagcor_handle == NULL) {
        fprintf(stderr,"ERROR: Failed to find the local calibration file %s\n",dmagcor_name);
        pPlate->errorFlag = 1;
        continue;
      } else {
        if (verbose > 1) {
          fprintf(stderr,"Found local calibration file %s\n",dmagcor_name);
        }
      }
      /* Now read in the local calibration file */

      dmagcor_header = table_header(dmagcor_handle,TABLE_PARSE);
      if (dmagcor_header == NULL) {
        fprintf(stderr,"ERROR: Failed to read header for %s\n",dmagcor_name);
        pPlate->errorFlag = 1;
        continue;
      }

      pPlate->dmagcor_table = table_loadva(dmagcor_handle,
                                           &dmagcor_header,
                                           NULL, /* hbase */
                                           NULL, /* rows */
                                           NULL,
                                           sizeof(DMAGCOR),
                                           &dmagcor_nrecs,
                                           TblInt,"nx",TblOff(PDMAGCOR,nx),
                                           TblInt,"ny",TblOff(PDMAGCOR,ny),
                                           TblInt,"ix",TblOff(PDMAGCOR,ix),
                                           TblInt,"iy",TblOff(PDMAGCOR,iy),
                                           TblInt,"npout",TblOff(PDMAGCOR,npout),
                                           TblInt,"rejectFlag",TblOff(PDMAGCOR,rejectFlag),
                                           TblDbl,"zout",TblOff(PDMAGCOR,zout),
                                           TblDbl,"errout",TblOff(PDMAGCOR,errout),
                                           TblDbl,"dradRMS2",TblOff(PDMAGCOR,dradRMS2),
                                           0,"end",0);
      if (pPlate->dmagcor_table == NULL) {
        fprintf(stderr,"ERROR: Failed to read table for %s\n",dmagcor_name);
        pPlate->errorFlag = 1;
        continue;
      }
      if (verbose > 1) {
        fprintf(stderr,"read %d records for %s\n",dmagcor_nrecs,dmagcor_name);
      }
      /* Now validate this table */
      pPlate->nx = pPlate->dmagcor_table[0].nx;
      pPlate->ny = pPlate->dmagcor_table[0].ny;
      if (dmagcor_nrecs != (pPlate->nx*pPlate->ny)) {
        fprintf(stderr,"ERROR: Number of records %d does not match product of %d and %d\n",dmagcor_nrecs,pPlate->nx,pPlate->ny);
        return(-1);
      }
      irec = 0;
      for (ix = 0; ix < pPlate->nx; ix++) {
        for (iy = 0; iy < pPlate->ny; iy++) {
          if ((pPlate->dmagcor_table[irec].nx != pPlate->nx) &&
              (pPlate->dmagcor_table[irec].ny != pPlate->ny) &&
              (pPlate->dmagcor_table[irec].ix != ix) &&
              (pPlate->dmagcor_table[irec].iy != iy)) {
            fprintf(stderr,"ERROR Mismatch of nx %d %d, ny %d %d, ix %d %d, or iy %d %d\n",
                    pPlate->dmagcor_table[irec].nx,pPlate->nx,
                    pPlate->dmagcor_table[irec].ny,pPlate->ny,
                    pPlate->dmagcor_table[irec].ix,ix,
                    pPlate->dmagcor_table[irec].iy,iy);
            return(-1);
          }
          irec++;
        }
      }

      if (dmagcor_header != NULL) {
        table_hdrfree(dmagcor_header);
      }
      if (dmagcor_handle != NULL) {
        Close(dmagcor_handle);
      }


      /* Open the extinction file */
      strcpy(extinction_name,catalogallDirectory);
      strcat(extinction_name,"/");
      strcat(extinction_name,pPlate->Plate);
      strcat(extinction_name,qualifier);
      strcat(extinction_name,"_extinction.db");


      extinction_handle = Open(extinction_name,"r");
      if (extinction_handle == NULL) {
        fprintf(stderr,"ERROR: Failed to find the extinction file %s\n",extinction_name);
        pPlate->errorFlag = 1;
        continue;
      } else {
        if (verbose > 1) {
          fprintf(stderr,"Found extinction file %s\n",extinction_name);
        }
      }

      /* Now read in the extinction file */

      extinction_header = table_header(extinction_handle,TABLE_PARSE);
      if (extinction_header == NULL) {
        fprintf(stderr,"ERROR: Failed to read header for %s\n",extinction_name);
        pPlate->errorFlag = 1;
        continue;
      }

      pPlate->extinction_table = table_loadva(extinction_handle,
                                              &extinction_header,
                                              NULL, /* hbase */
                                              NULL, /* rows */
                                              NULL,
                                              sizeof(EXTINCTION),
                                              &extinction_nrecs,
                                              TblInt,"enx",TblOff(PEXTINCTION,enx),
                                              TblInt,"eny",TblOff(PEXTINCTION,eny),
                                              TblInt,"eix",TblOff(PEXTINCTION,eix),
                                              TblInt,"eiy",TblOff(PEXTINCTION,eiy),
                                              TblDbl,"altitude",TblOff(PEXTINCTION,altitude),
                                              TblDbl,"extinction",TblOff(PEXTINCTION,extinction),
                                              0,"end",0);
      if (pPlate->extinction_table == NULL) {
        fprintf(stderr,"ERROR: Failed to read table for %s\n",extinction_name);
        pPlate->errorFlag = 1;
        continue;
      }
      if (verbose > 1) {
        fprintf(stderr,"read %d records for %s\n",extinction_nrecs,extinction_name);
      }
      /* Now validate this table */
      pPlate->enx = pPlate->extinction_table[0].enx;
      pPlate->eny = pPlate->extinction_table[0].eny;
      if (extinction_nrecs != (pPlate->enx*pPlate->eny)) {
        fprintf(stderr,"ERROR: Number of records %d does not match product of %d and %d\n",extinction_nrecs,pPlate->enx,pPlate->eny);
        return(-1);
      }
      irec = 0;
      for (eix = 0; eix < pPlate->enx; eix++) {
        for (eiy = 0; eiy < pPlate->eny; eiy++) {
          if ((pPlate->extinction_table[irec].enx != pPlate->enx) &&
              (pPlate->extinction_table[irec].eny != pPlate->eny) &&
              (pPlate->extinction_table[irec].eix != eix) &&
              (pPlate->extinction_table[irec].eiy != eiy)) {
            fprintf(stderr,"ERROR Mismatch of nx %d %d, ny %d %d, ix %d %d, or iy %d %d\n",
                    pPlate->extinction_table[irec].enx,pPlate->enx,
                    pPlate->extinction_table[irec].eny,pPlate->eny,
                    pPlate->extinction_table[irec].eix,eix,
                    pPlate->extinction_table[irec].eiy,eiy);
            return(-1);
          }
          irec++;
        }
      }


      if (extinction_header != NULL) {
        table_hdrfree(extinction_header);
      }
      if (extinction_handle != NULL) {
        Close(extinction_handle);
      }
    } /* ignoreLimiting */


  }


  mysql_close(&my_connection);

  /* Now sort the plate list in the order of increasing Julian Day */

  qsort((void*)plateArray,plateListLength,sizeof(PLATE),JulianDayCompare);
  return(0);
}
/* Writepoint returns +1 for a good entry, 0 for no entry, and -1 for a fatal error */

int WritePoint(PPLATE plateArray,int plateListLength,PMASTER pMaster,int *debugFlag) 
{
  PPLATE pPlate;
  PDMAGCOR pDmagcor;
  int plateIndex;
#if 0
  if ((strcmp(pMaster->REF,"N120013341") == 0) &&
      (strcmp(pMaster->Plate,"dsb02822") == 0)) {
    printf("At plate %s for %s\n",pMaster->Plate,pMaster->REF);
  }
#endif

  if (*debugFlag) {
    fprintf(stderr,"Saving point for %s ",pMaster->REF);
    *debugFlag = 0;
  }

  for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
    pPlate = &plateArray[plateIndex];
    if (strcmp(pMaster->Plate,pPlate->PlateStub) == 0) {
      /* We have a match */
      if (pPlate->starDetected) {
        fprintf(stderr,"ERROR: WritePoint already recorded plate %s for star %s\n",pMaster->Plate,pMaster->REF);
        return(0);
      }
      pPlate->magcal_local = pMaster->magcal_local;
      pPlate->Date = pMaster->Date;
      pPlate->spatial_bin = pMaster->spatial_bin;
      pPlate->magcal_local_rms = pMaster->magcal_local_rms;
      pPlate->magcal_local_error = pMaster->magcal_local_error;
      pPlate->Stdmag = pMaster->Stdmag;
      pPlate->magcal_iso_rms = pMaster->magcal_iso_rms;
      pPlate->magcal_iso = pMaster->magcal_iso;
      pPlate->limiting_mag_local = pMaster->limiting_mag_local;
      pPlate->npoints_local = pMaster->npoints_local;
      pPlate->AFLAGS = pMaster->AFLAGS;
      pPlate->AFLAGSCOPY = pMaster->AFLAGSCOPY;
      pPlate->BFLAGS = pMaster->BFLAGS;
      pPlate->gsc_bin_index = pMaster->gsc_bin_index;
      pPlate->starDetected = 1;
      if (pPlate->dmagcor_table != NULL) {
        pDmagcor = &pPlate->dmagcor_table[pMaster->local_bin_index];
        pPlate->dradRMS2 = pDmagcor->dradRMS2;
      } else {
        pPlate->dradRMS2 = 99.0;
      }
      return(1);
    }
  }
#if 1
  return(0);
#else
  fprintf(stderr,"ERROR: WritePoint Failed to find plate %s in plate list\n",pMaster->Plate);
  return(-1);
#endif
}
void MergeBlendArray(char *blendStrings,int blendCount,PPLATE *blendArray,int *savedBlendCount, PTARGET pTarget, int plateListLength) {
  int blendIndex;
  PPLATE plateArray1;
  PPLATE plateArray2;
  PPLATE pPlate1;
  PPLATE pPlate2;
  double clipmed;
  int plateIndex;

  memset(pTarget,0,sizeof(TARGET));
  for (blendIndex = 0; blendIndex < blendCount; blendIndex++) {
    plateArray2 = blendArray[blendIndex];
    strcat(pTarget->REF,plateArray2->pTarget->REF);
    strcat(pTarget->REF,"_");
    strcat(pTarget->src_name,plateArray2->pTarget->src_name);
    strcat(pTarget->src_name,"_");



    if (blendIndex == 0) {
      plateArray1 = plateArray2;
      continue;
    }
    clipmed = plateArray2->clipmed;
    if (clipmed == 0.0) {
      for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
        pPlate2 = &plateArray2[plateIndex];
        if (pPlate2->starDetected != 0) {
          clipmed = pPlate2->Stdmag;
          break;
        }

      }
    }
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate1 = &plateArray1[plateIndex];
      pPlate2 = &plateArray2[plateIndex];
      pPlate1->AFLAGS = pPlate1->AFLAGSCOPY;
      pPlate2->AFLAGS = pPlate2->AFLAGSCOPY;
      pPlate1->validPoint = 0;
#if 0
      if (strcmp(pPlate1->Plate,"rb09369_00_01r180ww") == 0) {
        printf("At plate %s\n",pPlate1->Plate);
      }
#endif
#if 0
      printf("Index %d, errorFlag %d starDetected %d\n",plateIndex,pPlate1->errorFlag,pPlate1->starDetected);
#endif
      if ((pPlate1->errorFlag == 0) &&
          (pPlate1->starDetected != 0)) {
        double tmpFlux1;
        double tmpFlux2;
        double totalFlux;
        /* Star of interest has been detected.  See if the companion has also been detected,
           but be sure it is not an unmatched blended object */
        if ((pPlate2->errorFlag == 0) &&
            (pPlate2->starDetected != 0) &&
            ((pPlate2->AFLAGS & (1<<FILTER_AFLAG_BLEND_NOMATCH)) == 0) &&
            (pPlate2->magcal_local < 90.0) &&
            (pPlate2->magcal_iso < 90.0)) {          
          if (((pPlate1->AFLAGS & (1<<FILTER_AFLAG_BLEND)) == 0) &&
              ((pPlate2->AFLAGS & (1<<FILTER_AFLAG_BLEND)) == 0) &&
              ((pPlate1->AFLAGS & (1<<FILTER_AFLAG_CASED)) != 0) &&
              ((pPlate2->AFLAGS & (1<<FILTER_AFLAG_CASED)) != 0)) {
            /* Here both stars are part of a complex blend, yet neither one has been
               detected as blended by Sextractor.  Assume that they are not really
               blended but marked as such because of elliptical rather than circular trails */
            /* Clear the blend flag */
            pPlate1->AFLAGS &= ~(1<<FILTER_AFLAG_CASED);
            pPlate2->AFLAGS &= ~(1<<FILTER_AFLAG_CASED);
            if ((pPlate1->AFLAGS < FILTER_AFLAG_BAD) && (pPlate2->AFLAGS < FILTER_AFLAG_BAD)) {
              fprintf(stderr,"WARNING: %s and %s for plate %s are not really blended\n",
                      pTarget->REF,plateArray2->pTarget->REF,pPlate1->Plate);
            }
          }


          /* Other star has been detected, sum up the magnitudes */
          tmpFlux1 = exp10(-pPlate1->magcal_local/2.5);
          tmpFlux2 = exp10(-pPlate2->magcal_local/2.5);
          totalFlux = tmpFlux1+tmpFlux2;
          if (totalFlux == 0) {
            pPlate1->magcal_local = 99.0;
            pPlate1->magcal_local_rms = 99.0;
          } else {
            pPlate1->magcal_local = -2.5*log10(totalFlux);
          }
          tmpFlux1 = exp10(-pPlate1->magcal_iso/2.5);
          tmpFlux2 = exp10(-pPlate2->magcal_iso/2.5);
          totalFlux = tmpFlux1+tmpFlux2;
          if (totalFlux == 0) {
            pPlate1->magcal_iso = 99.0;
            pPlate1->magcal_iso_rms = 99.0;
          } else {
            pPlate1->magcal_iso = -2.5*log10(totalFlux);
          }
        } else {
          /* Other star has not been detected */
          if (((pPlate2->AFLAGS & (1<<FILTER_AFLAG_BLEND_NOMATCH)) == 0) ||
              (pPlate2->starDetected == 0)) {
            /* Other star is not part of primary star */
            if (clipmed != 0.0) {
              /* Use the DASCH median magnitude for the undetected star */
              tmpFlux1 = exp10(-pPlate1->magcal_local/2.5);
              tmpFlux2 = exp10(-clipmed/2.5);
              totalFlux = tmpFlux1+tmpFlux2;
              if (totalFlux == 0) {
                pPlate1->magcal_local = 99.0;
                pPlate1->magcal_local_rms = 99.0;
              } else {
                pPlate1->magcal_local = -2.5*log10(totalFlux);
              }
              tmpFlux1 = exp10(-pPlate1->magcal_iso/2.5);
              tmpFlux2 = exp10(-clipmed/2.5);
              totalFlux = tmpFlux1+tmpFlux2;
              if (totalFlux == 0) {
                pPlate1->magcal_iso = 99.0;
                pPlate1->magcal_iso_rms = 99.0;
              } else {
                pPlate1->magcal_iso = -2.5*log10(totalFlux);
              }

            }

          }

        }
      }
    }
  }
  pTarget->REF[strlen(pTarget->REF)-1] = 0;
  pTarget->src_name[strlen(pTarget->src_name)-1] = 0;
}

/* If the current object matches the blend list, then save the array */
void SaveLightcurve(PPLATE plateArray,int plateListLength,PTARGET pTarget,char *blendStrings,int blendCount,PPLATE *blendArray,int *savedBlendCount,double clipmed)
{
  int blendIndex;
  int saveFlag = 0;
  int curBlendCount = *savedBlendCount;
  if ((curBlendCount < 0) || (curBlendCount >= MAX_BLEND_COUNT)) {
    fprintf(stderr,"ERROR: curBlendCount %d too high for MAX_BLEND_COUNT %d\n",curBlendCount,MAX_BLEND_COUNT);
  }
  if (blendCount == 0) {
    return;
  }
  for (blendIndex = 0; blendIndex < blendCount; blendIndex++) {
    if (strcmp(pTarget->REF,&blendStrings[blendIndex*MAX_BUFFER]) == 0) {
      saveFlag = 1;
      break;
    }

  }
  if (saveFlag == 0) {
    return;
  }
  blendArray[blendIndex] = (PPLATE)calloc(plateListLength,sizeof(PLATE));
  if (blendArray[blendIndex] == NULL) {
    fprintf(stderr,"ERROR: Failed to allocate blend array of size %d\n",plateListLength*sizeof(PLATE));
    exit(-1);
  }
  plateArray->clipmed = clipmed;
  memcpy(blendArray[blendIndex],plateArray,plateListLength*sizeof(PLATE));
  blendArray[blendIndex]->pTarget = pTarget;
  curBlendCount++;
  *savedBlendCount = curBlendCount;
  printf("Saved blend object %s %s index %d\n",pTarget->src_name,pTarget->REF,blendIndex);
  return;
}
int WriteLightcurve(PGSC_CACHE_ENTRY pGscCache,File indexHandle, File catalogHandle,PPLATE plateArray,int plateListLength,int match_count,PTARGET pTarget,char *lightcurveDirectory,double *vector1,double *vector2,int verbose,int goodFlag,int printReject,double *pClipMed, int ignoreDefect,int ignoreLimiting,int allLightcurves)
{
  PPLATE pPlate;
  int plateIndex;
  int curveWritten = 0;
  int sequence = 0;
  int validPoints = 0;
  int goodPoints = 0;
  int errorPoints = 0;
  int limitingPoints1 = 0;
  int limitingPoints2 = 0;
  int badPoints = 0;
  int bin9Points = 0;
  int bothPoints = 0;
  double X_IMAGE;
  double Y_IMAGE;
  int offscl;
  double edgeDist;
  int ix;
  int iy;
  int eix;
  int eiy;
  int irec;
  double rawmed = 0;
  double rawrms = 0;
  double clipmed = 0;
  double cliprms = 0;
  double rawerrrms = 0;
  double rawerrmed = 0;
  PDMAGCOR pDmagcor;
  PEXTINCTION pExtinction;
  char lightcurveName[MAX_BUFFER];
  char MAGFlagString[MAX_BUFFER];
  FILE *lightcurveHandle;
  double startYear = 0.0;
  double endYear;
  double startJD = 0.0;
  double endJD;
  double rangeJD;
  int wedgeCount = 0;
  int highDradCount = 0;
  int defectCount = 0;
  int hiZoutCount = 0;
  int tooBrightCount = 0;
  int lowAltitudeCount = 0;
  int blendCount = 0;
  int localBinFailure = 0;
  int reject_count_all1[REJECT_REASON_MAX];
  int reject_count_only1[REJECT_REASON_MAX];
  int reject_count_all2[REJECT_REASON_MAX];
  int reject_count_only2[REJECT_REASON_MAX];
  int rejectIndex;
  double error_bar_factor = 1.0;
  int AFLAGSMASK;
  int BFLAGSMASK = -1;
  PGSCIMAGE pGscImage = NULL;
  int newGscBinIndex;
  int tmpDecBin;
  int tmpRaBin;
  int plateIndex2;
  char rstr[32], dstr[32];

  *pClipMed = 0;
  if (printReject) {
    for (rejectIndex = 0; rejectIndex < REJECT_REASON_MAX; rejectIndex++) {
      reject_count_all1[rejectIndex] = 0;
      reject_count_only1[rejectIndex] = 0;
      reject_count_all2[rejectIndex] = 0;
      reject_count_only2[rejectIndex] = 0;
    }
  }

  /* Go through the array and set in the correct flags value.  Also
     calculate and set in the limiting magnitude if the star was not
     detected */
  for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
    pPlate = &plateArray[plateIndex];
    pPlate->reject_reason1 = 0;
    pPlate->reject_reason2 = 0;

#if 0
  if ((strcmp(pTarget->REF,"N120013341") == 0) &&
      (strcmp(pPlate->Plate,"dsb02822_00_01ww") == 0)) {
    printf("At plate %s for %s\n",pPlate->Plate,pTarget->REF);
  }
#endif

    localBinFailure = 0;


#if 0
    if (strcmp(pTarget->REF,"N120013341") == 0) {
      printf("Writing lightcurve for %s, plate %s\n",pTarget->REF,pPlate->Plate);
    }
#endif
#if 0
    if (strcmp(pPlate->PlateStub,"dnb03443") == 0) {
      printf("Writing lightcurve for %s, plate %s\n",pTarget->REF,pPlate->Plate);
    }
#endif
    if (pPlate->errorFlag != 0) {
      pPlate->reject_reason1 |= (1 << REJECT_REASON_ERROR);
      pPlate->reject_reason2 |= (1 << REJECT_REASON_ERROR);
    } else {
      if (pPlate->starDetected == 0) {
        wcs2pix(pPlate->wcs,pTarget->ra,pTarget->dec,&X_IMAGE,&Y_IMAGE,&offscl);
        if (offscl) {
          pPlate->reject_reason1 |= (1 << REJECT_REASON_OFFPLATE);
          pPlate->reject_reason2 |= (1 << REJECT_REASON_OFFPLATE);
        } else {
          pPlate->spatial_bin = CalculateBin(pPlate->mosaic.naxis1,pPlate->mosaic.naxis2,X_IMAGE,Y_IMAGE,&edgeDist);
          if ((pPlate->spatial_bin < 1) || (pPlate->spatial_bin > MAX_SPATIAL_BINS)) {
            pPlate->reject_reason1 |= (1 << REJECT_REASON_OFFPLATE);          
            pPlate->reject_reason2 |= (1 << REJECT_REASON_OFFPLATE);          
          } else {
            if (pPlate->spatial_bin == 9) {
              pPlate->reject_reason1 |= (1 << REJECT_REASON_BIN9);
              pPlate->reject_reason2 |= (1 << REJECT_REASON_BIN9);
            }
            if (ignoreLimiting) {
              pPlate->reject_reason1 |= (1 << REJECT_REASON_UNDETECTED);
              pPlate->reject_reason2 |= (1 << REJECT_REASON_UNDETECTED);
            } else { /* ignoreLimiting */
              if (pPlate->spatial_bin != pPlate->spatial_bin_table[pPlate->spatial_bin].spatial_bin) {
                pPlate->reject_reason1 |= (1 << REJECT_REASON_BINFAILURE);
                pPlate->reject_reason2 |= (1 << REJECT_REASON_BINFAILURE);
              } else {
                pPlate->reject_reason1 |= (1 << REJECT_REASON_UNDETECTED);
                pPlate->reject_reason2 |= (1 << REJECT_REASON_UNDETECTED);
              }

              /* Now see if there is a local correction for this magnitude */
              if ((pPlate->mosaic.naxis1 != 0) &&
                  (pPlate->mosaic.naxis2 != 0)) {
                ix = (X_IMAGE * pPlate->nx) /(1.0 * pPlate->mosaic.naxis1);
                iy = (Y_IMAGE * pPlate->ny) /(1.0 * pPlate->mosaic.naxis2);
                /* These records are reversed from the local bin */
                irec = iy + (pPlate->ny*ix);
                if (pPlate->dmagcor_table != NULL) { 
                  pDmagcor = &pPlate->dmagcor_table[irec];
                  pPlate->dradRMS2 = pDmagcor->dradRMS2;
                  if ((pDmagcor->errout >= 90.0) ||
                      (pDmagcor->rejectFlag != 0)) {
                    pPlate->reject_reason1 |= (1 << REJECT_REASON_HIZOUT);
                    pPlate->reject_reason2 |= (1 << REJECT_REASON_HIZOUT);
                  }
                }
              }
            } /* ignoreLimiting */
          }
        }
      
      } /* Star detected */ 
    } /* Error flag clear */

    if (pPlate->errorFlag == 0) {
      if (pPlate->starDetected) {

        if (goodFlag == 0) {
          AFLAGSMASK = 0xffffffff;
        } else {
          AFLAGSMASK = (1<<FILTER_AFLAG_CASEB) | (1<<FILTER_AFLAG_BLEND);
          if (ignoreDefect) {
            AFLAGSMASK |= (1<<FILTER_AFLAG_DEFECT);
          } 
        }
#ifdef ACCEPT_LOW_ALTITUDE
        AFLAGSMASK |= (1<<FILTER_AFLAG_LOW_ALTITUDE);
#endif /* ACCEPT_LOW_ALTITUDE */

        pPlate->validPoint = DecodeAFLAGS(pPlate->AFLAGS,AFLAGSMASK,pPlate->BFLAGS,BFLAGSMASK,&pPlate->reject_reason1,&pPlate->reject_reason2,&pPlate->flag);
        if (pPlate->flag == SYMBOL_WEDGE) {
          wedgeCount++;
        } else if (pPlate->flag == SYMBOL_DRAD) {
          highDradCount++;
        } else if (pPlate->flag == SYMBOL_HIZOUT) {
          hiZoutCount++;
        } else if (pPlate->flag == SYMBOL_TOO_BRIGHT) {
          if ((pPlate->reject_reason2 & (1<<REJECT_REASON_TOO_BRIGHT)) != 0) {
            tooBrightCount++;
          }         
        } else if (pPlate->flag == SYMBOL_LOW_ALTITUDE) {
          lowAltitudeCount++;
        } else if (pPlate->flag == SYMBOL_DEFECT) {
          defectCount++;
        }
        if ((pPlate->reject_reason2 & (1<<REJECT_REASON_BLEND)) != 0) {
          blendCount++;
        }
        if ((pPlate->reject_reason2 & (1<<REJECT_REASON_BIN9)) != 0) {
          bin9Points++;
        }
        if ((pPlate->reject_reason2 & (1<<REJECT_REASON_LOCAL_RMS)) != 0) {
          badPoints++;
        }
        if ((pPlate->flag == SYMBOL_GOOD) || (pPlate->flag == SYMBOL_NOCOLOR) || (pPlate->flag == SYMBOL_NEIGHBORS) || (pPlate->flag == SYMBOL_SATURATED) || (pPlate->flag == SYMBOL_NOMAGDEP)) {
          /* This currently includes stars near the limiting magnitude */
          vector1[goodPoints] = pPlate->magcal_local;
          vector2[goodPoints] = pPlate->magcal_local_rms;
#if 0
          printf("DEBUG: point %4d AFLAGS %08x plate %25s target %s\n",goodPoints,pPlate->AFLAGS,pPlate->Plate,pTarget->REF);
#endif
          goodPoints++;
          if (pPlate->validPoint == 0) {
            limitingPoints2++;
          }
        }
        if (pPlate->validPoint != 0) {
          validPoints++;
        }
          
      

      } else {
        if (ignoreLimiting) {
          continue;
        }
        /* Star has not been detected. Figure out where this
           star is on the plate */
        wcs2pix(pPlate->wcs,pTarget->ra,pTarget->dec,&X_IMAGE,&Y_IMAGE,&offscl);
        if (offscl) {
          errorPoints++;
          continue;
        }
        pPlate->spatial_bin = CalculateBin(pPlate->mosaic.naxis1,pPlate->mosaic.naxis2,X_IMAGE,Y_IMAGE,&edgeDist);
        if ((pPlate->spatial_bin < 1) || (pPlate->spatial_bin > MAX_SPATIAL_BINS)) {
          errorPoints++;
          continue;
        }
        if (pPlate->spatial_bin != pPlate->spatial_bin_table[pPlate->spatial_bin].spatial_bin) {
          errorPoints++;
          continue;
        }
        limitingPoints1++;
        if (goodFlag == 1) {
          continue;
        }
        /* Start with the iso upper limit */
        pPlate->magcal_iso = 0;
        pPlate->magcal_iso_rms = 0;
        pPlate->limiting_mag_local = pPlate->spatial_bin_table[pPlate->spatial_bin].limiting_mag;
        pPlate->magcal_local_rms = 0;
        pPlate->npoints_local = 0;
        pPlate->magcal_local = pPlate->spatial_bin_table[pPlate->spatial_bin].upper_limit;
        pPlate->magcal_iso = pPlate->magcal_local;

        /* Now see if there is a local correction for this magnitude */
        ix = (X_IMAGE * pPlate->nx) /(1.0 * pPlate->mosaic.naxis1);
        iy = (Y_IMAGE * pPlate->ny) /(1.0 * pPlate->mosaic.naxis2);
        /* These records are reversed from the local bin */
        irec = iy + (pPlate->ny*ix);
        pDmagcor = &pPlate->dmagcor_table[irec];
        if ((pDmagcor->errout < 90.0) &&
            (pDmagcor->rejectFlag == 0)) {
          pPlate->npoints_local = pDmagcor->npout;
          pPlate->magcal_local -= pDmagcor->zout;
          pPlate->limiting_mag_local -= pDmagcor->zout;

        } else {
          localBinFailure = 1;
        }


        /* Finally, see if there is an extinction correction */
        eix = (X_IMAGE * pPlate->enx) /(1.0 * pPlate->mosaic.naxis1);
        eiy = (Y_IMAGE * pPlate->eny) /(1.0 * pPlate->mosaic.naxis2);
        /* These records are reversed from the local bin */
        irec = eiy + (pPlate->eny*eix);
        pExtinction = &pPlate->extinction_table[irec];
        if (pExtinction->extinction > 0) {
          pPlate->magcal_local -= pExtinction->extinction;
          pPlate->magcal_iso -= pExtinction->extinction;
          pPlate->limiting_mag_local -= pExtinction->extinction;
        }
         

        if (localBinFailure) {
          pPlate->flag = 1.0;
          continue;     /* Here do not allow the points to clutter the plot */
        } else if (pPlate->spatial_bin == 9) {
          pPlate->flag = 31.5;
        } else {
          pPlate->flag = 31.0;
        }

        pPlate->magcal_iso_rms = 0.0;

        pPlate->Date = jd2hjd(pPlate->geoJulianDate,pTarget->ra,pTarget->dec,WCS_J2000);
        pPlate->validPoint = 1;
        validPoints++;
      }
    } else {
      errorPoints++;
    }

  }
  if (((allLightcurves == 0) && (validPoints > 1)) ||
      ((allLightcurves > 0) && (validPoints >= allLightcurves))) {

    if (CalcMedianAndRMS(goodPoints,2,vector1,&rawmed,&rawrms,0,3.0,0) == 0) {
      rawmed = 0.0;
      rawrms = 99.0;
    }
    if (CalcMedianAndRMS(goodPoints,2,vector1,&clipmed,&cliprms,1,3.0,0) == 0) {
      clipmed = 0.0;
      cliprms = 99.0;
    }
    *pClipMed = clipmed;
    /* Find the zero-based clipped rms of our error bars */
    if (CalcMedianAndRMS(goodPoints,2,vector2,&rawerrmed,&rawerrrms,0,3.0,1) == 0) {
      rawerrmed = 0.0;
      rawerrrms = 99.0;
    }
    if ((rawrms == 0.0) || 
        (rawrms == 99.0) ||
        (rawerrrms == 0.0) || 
        (rawerrrms == 99.0)) {
      error_bar_factor = 1.0;
    } else {
      error_bar_factor = rawrms/rawerrrms;
    }
    if (error_bar_factor > 1.0) {
      error_bar_factor = 1.0;
    }


    /* We have at least two good points, write the curve */
    strcpy(lightcurveName,lightcurveDirectory);
    strcat(lightcurveName,"/lc_");
    strcat(lightcurveName,pTarget->src_name);
    if (goodFlag) {
      strcat(lightcurveName,"_good");
    } else {
      strcat(lightcurveName,"_all");
    }
    strcat(lightcurveName,".db");
    lightcurveHandle = fopen(lightcurveName,"wt");
    if (lightcurveHandle == NULL) {
      fprintf(stderr,"ERROR: Failed to open %s\n",lightcurveName);
      return(-1);
    }
    fprintf(lightcurveHandle,"sequence\thJD\tspatial_bin\tflag\tmagcal_local\tmagcal_local_rms\tnlocal\tPlate\tmagcal_iso_rms\tmagcal_local_error\tmagcal_iso\tlimiting_mag_local\tAFLAGS\treject_reason1\treject_reason2\tdradRMS2");
    fprintf(lightcurveHandle,"\n");

    fprintf(lightcurveHandle,"--------\t---\t-----------\t----\t------------\t----------------\t------\t-----\t--------------\t-----------------\t----------\t------------------\t-----\t--------------\t--------------\t--------");
    fprintf(lightcurveHandle,"\n");


    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if ((pPlate->errorFlag == 0) &&
          (pPlate->validPoint != 0)) {
        endJD = pPlate->Date;
        if (startJD == 0.0) {
          startJD = endJD;
        }
        sequence++;
        fprintf(lightcurveHandle,"%d\t%f\t%d\t%.1f\t%.2f\t%.2f\t%d\t%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f",
                sequence,
                pPlate->Date,
                pPlate->spatial_bin,
                pPlate->flag,
                pPlate->magcal_local,
                pPlate->magcal_local_rms*error_bar_factor,
                pPlate->npoints_local,
                pPlate->PlateStub,
                pPlate->magcal_iso_rms*error_bar_factor,
                pPlate->magcal_local_error*error_bar_factor,
                pPlate->magcal_iso,
                pPlate->limiting_mag_local,
                pPlate->AFLAGSCOPY,
                pPlate->reject_reason1,
                pPlate->reject_reason2,
                pPlate->dradRMS2);
        fprintf(lightcurveHandle,"\n");

      }

    }
    fclose(lightcurveHandle);
    /* Now write out the summary text */
    rangeJD = (endJD-startJD);
    endJD += MARGIN_WIDTH * rangeJD;
    startJD -= MARGIN_WIDTH *rangeJD;

    startYear = jd2ep(startJD);
    endYear = jd2ep(endJD);
    strcpy(lightcurveName,lightcurveDirectory);
    strcat(lightcurveName,"/lc_");
    strcat(lightcurveName,pTarget->src_name);
    if (goodFlag) {
      strcat(lightcurveName,"_good");
    } else {
      strcat(lightcurveName,"_all");
    }
    if (pTarget->gsc_bin_index == 0) {
      /* We need to get get the catalog information */
      pGscImage = NULL;
      for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
        pPlate = &plateArray[plateIndex];
        if (pPlate->gsc_bin_index != 0) {
          pGscImage = SearchStar(pGscBin,pGscCache,indexHandle,catalogHandle,pPlate->gsc_bin_index,&newGscBinIndex,pTarget->REF);
          if (pGscImage != NULL) {
            break;
          } else {
            for (plateIndex2 = plateIndex+1; plateIndex2 < plateListLength; plateIndex2++) {
              if (plateArray[plateIndex2].gsc_bin_index == pPlate->gsc_bin_index) {
                plateArray[plateIndex2].gsc_bin_index = 0;
              }
            }
        
          }
        }
      }
      if (pGscImage == NULL) {
        int tmpGscBinIndex;
        /* Try to get the gsc_bin_index from the star ra and dec */
        tmpGscBinIndex = GetGSCBin(pGscBin,pTarget->ra,pTarget->dec,&tmpDecBin,&tmpRaBin,(char *)"extract_lightcurves");
        pGscImage = SearchStar(pGscBin,pGscCache,indexHandle,catalogHandle,tmpGscBinIndex,&newGscBinIndex,pTarget->REF);
        


      }
      if (pGscImage != NULL) {
        pTarget->gsc_bin_index = newGscBinIndex;
        pTarget->Stdmag = pGscImage->Stdmag;
        pTarget->color = pGscImage->color;
        pTarget->VFlag = pGscImage->VFlag;
        pTarget->MAGFlag = pGscImage->MAGFlag;
        pTarget->RaPM = pGscImage->RaPM;
        pTarget->DecPM = pGscImage->DecPM;
        
      }
      
    }
 

    strcat(lightcurveName,".txt");
    lightcurveHandle = fopen(lightcurveName,"wt");

    fprintf(lightcurveHandle,"#raw median mag, raw median rms, clipped median mag, clipped median rms, total points, good points, start year, end year,error_bar_factor, start JD,end JD,REF,ra,dec,Stdmag,color,VFlag,MAGFlag,RaPM,DecPM\n");

    fprintf(lightcurveHandle,"%.2f %.2f %.2f %.2f %d %d %f %f %f %.0f %.0f ",
            rawmed,rawrms,clipmed,cliprms,validPoints,goodPoints,startYear,endYear,error_bar_factor,startJD,endJD);
    if (pTarget->gsc_bin_index != 0) {
      ra2str (rstr, 16, pTarget->ra, ndec);
      dec2str (dstr, 16, pTarget->dec, ndec-1);    

      if ((pTarget->MAGFlag & KEPLER_MAGNITUDE_FLAG) == 0) {
        sprintf(MAGFlagString,"%d",pTarget->MAGFlag);
      } else {
        if ((pTarget->MAGFlag & MAGNITUDE_MASK) < KEPLER_CQ_MAGFLAG_MAX) {
          strcpy(MAGFlagString,keplerSourceText[pTarget->MAGFlag & MAGNITUDE_MASK]);
        } else {
          sprintf(MAGFlagString,"%d??",pTarget->MAGFlag);
        }
      }

      fprintf(lightcurveHandle,"%s %s %s %.2f %.2f %d %s %.1f %.1f\n",
              pTarget->REF,
              rstr,
              dstr,
              pTarget->Stdmag,
              pTarget->color,
              pTarget->VFlag,
              MAGFlagString,
              pTarget->RaPM,
              pTarget->DecPM);
    } else {
      fprintf(lightcurveHandle,"UNK UNK UNK UNK UNK UNK UNK UNK UNK\n");

    }

    if (lightcurveHandle == NULL) {
      fprintf(stderr,"ERROR: Failed to open %s\n",lightcurveName);
      return(-1);
    }
    fclose(lightcurveHandle);
    if (allLightcurves == 0) {
      fprintf(stderr,"Lightcurve written goodFlag %d good: %4d err: %4d bad: %4d limiting %4d %4d bin9 %4d both %4d wedge %4d high drad %4d defect %4d hiZout %4d tooBright %d too low %d blend %4d errorbarfactor %f for %s %s\n",goodFlag,goodPoints,errorPoints,badPoints,limitingPoints1,limitingPoints2,bin9Points,bothPoints,wedgeCount,highDradCount,defectCount,hiZoutCount,tooBrightCount,lowAltitudeCount,blendCount,error_bar_factor,pTarget->src_name,pTarget->REF);
    }
    curveWritten = 1;



  } else {
    if (verbose & (allLightcurves == 0)) {
     
      fprintf(stderr,"No valid points found for %s\n",pTarget->src_name);
    }
  }
             
  if (printReject) {
    char rejectFileName[MAX_BUFFER];
    FILE *rejectFileHandle = NULL;
    char undetectedName[MAX_BUFFER];
    FILE *undetectedHandle = NULL;
    strcpy(rejectFileName,lightcurveDirectory);
    strcat(rejectFileName,"/lc_");
    strcat(rejectFileName,pTarget->src_name);
    strcat(rejectFileName,"_reject_reason.txt");
    rejectFileHandle = fopen(rejectFileName,"wt");
    if (rejectFileHandle == NULL) {
      fprintf(stderr,"ERROR: Failed to open %s\n",rejectFileName);
      return(-1);
    }
    strcpy(undetectedName,lightcurveDirectory);
    strcat(undetectedName,"/lc_");
    strcat(undetectedName,pTarget->src_name);
    strcat(undetectedName,"_undetected.txt");
    undetectedHandle = fopen(undetectedName,"wt");
    if (undetectedHandle == NULL) {
      fprintf(stderr,"ERROR: Failed to open %s\n",undetectedName);
      return(-1);
    }
    fprintf(undetectedHandle,"Plate\treject_reason\n");
    fprintf(undetectedHandle,"-----\t-------------\n");

    /* Tally up the individual counts */
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (pPlate->reject_reason2 & ((1<<REJECT_REASON_ERROR)|
                                    (1<<REJECT_REASON_OFFPLATE)|
                                    (1<<REJECT_REASON_BINFAILURE)|
                                    (1<<REJECT_REASON_UNDETECTED))) {
        fprintf(undetectedHandle,"%s\t%d\n",pPlate->PlateStub,pPlate->reject_reason2);
      }
      for (rejectIndex = 0; rejectIndex < REJECT_REASON_MAX; rejectIndex++) {
        if ((pPlate->reject_reason1 & (1 << rejectIndex)) != 0) {
          reject_count_all1[rejectIndex]++;
          if (pPlate->reject_reason1 == (1 << rejectIndex)) {
            reject_count_only1[rejectIndex]++;
          }
        }
        if ((pPlate->reject_reason2 & (1 << rejectIndex)) != 0) {
          reject_count_all2[rejectIndex]++;
          if (pPlate->reject_reason2 == (1 << rejectIndex)) {
            reject_count_only2[rejectIndex]++;
          }
        }
      }
    }
    fprintf(rejectFileHandle,"Reject reasons of %d plates for %s %s\n",
            plateListLength,
            pTarget->REF,
            pTarget->src_name);
    for (rejectIndex = 0; rejectIndex < REJECT_REASON_MAX; rejectIndex++) {
      fprintf(rejectFileHandle,"%4d %4d %4d %4d %4d %s\n",
              rejectIndex,
              reject_count_all1[rejectIndex],
              reject_count_only1[rejectIndex],
              reject_count_all2[rejectIndex],
              reject_count_only2[rejectIndex],
              GetRejectReasonText(rejectIndex));
    }
    fclose(rejectFileHandle);
    fclose(undetectedHandle);
  }



  return(curveWritten);
}






int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  char list_name[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  char lightcurveDirectory[MAX_BUFFER];
  char rmsName[MAX_BUFFER];
  FILE *rmsHandle = NULL;

  char *dotPos;
  char *slashPos;
  char *lastSlashPos;
 
  FILE *list_handle;
  time_t startTime;
  time_t curTime;
  int verbose = 0;
  int allLightcurves = 0;
  int printReject = 0;
  int debugFlag = 0;
  int outputPlotCount = 0;
  int outputResult;
  int goodFlag = 0;
  char *ingestDirectory;
  char *catalogallDirectory;
  
  char inLine[MAX_BUFFER];
  int lineLen;
  int nlines = 0;
  char *inBuffer;
  PPLATE plateArray = NULL;
  double *vector1 = NULL;
  double *vector2 = NULL;
  int plateListLength = 0;
  int compareFlag;
  int match_count = 0;
  PPLATE pPlate;
  int plateIndex;
  double clipmed;


  int ignoreLimiting = 0;

  File target_handle = NULL;
  char target_name[MAX_BUFFER];
  TableHead target_header = NULL;
  PTARGET target_table = NULL;
  size_t target_nrecs = 0;
  int target_index;
  PTARGET pTarget;
  TARGET oldTarget;
  PTARGET pOldTarget = &oldTarget;
  TARGET blendTarget;

  File master_handle = NULL;
  char master_name[MAX_BUFFER];
  TableHead master_header = NULL;
  int master_nrecs = 0;
  int master_index;
  MASTER old_master_record;
  PMASTER pOldMaster = &old_master_record;
  MASTER master_record;
  PMASTER pMaster = &master_record;
  TblDescriptor master_descriptor = NULL;
  TableRow master_row = NULL;

  int errorFlag = 0;

  int blendCount = 0;
  char blendStrings[MAX_BLEND_COUNT*MAX_BUFFER];
  PPLATE blendArray[MAX_BLEND_COUNT];
  int savedBlendCount = 0;
  int blendIndex;
  int ignoreDefect = 0;
  int writePointResult;
  GSC_CACHE_ENTRY pGscCache[GSC_CACHE_ENTRIES];
  char indexname[MAX_BUFFER];
  File indexHandle;
  char* catalogdir;
  char* slashPtr;
  char catalogname[MAX_BUFFER];
  File catalogHandle;
  char *dotPtr;


#ifdef ACCEPT_LOW_ALTITUDE
  printf("ERROR: FILTER_AFLAG_LOW_ALTITUDE stars are being plotted\n");
#endif /* ACCEPT_LOW_ALTITUDE */


  printf("ERROR: This file is obsolete and will not handle multiple exposures correctly\n");


  memset(pGscCache,0,sizeof(pGscCache));
  time(&startTime);
  list_name[0] = 0;
  target_name[0] = 0;
  qualifier[0] = 0;
  master_name[0] = 0;
  lightcurveDirectory[0] = 0;

  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_INGEST is not defined\n");
    return(-1);
  }
  catalogallDirectory = getenv("DASCH_CATALOGALL");
  if (catalogallDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_CATALOGALL is not defined\n");
    return(-1);
  }



  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      fprintf(stderr,"ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {


        case 'l': /* list file name */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(list_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'q': /* file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(qualifier,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'b': /* blended object name */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            if (blendCount < MAX_BLEND_COUNT) {
              strncpy(&blendStrings[blendCount*MAX_BUFFER],*++argv,MAX_BUFFER-2);
              blendCount++;
              if (strlen(*argv) >= MAX_BUFFER-2) {
                fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
                errorFlag = 1;
              }
            } else {
              fprintf(stderr,"ERROR: More than %d blend stars\n",blendCount);
              errorFlag = 1;                          
            }
          }
          break;

        

        case 't': /* target catalog file name */
        case 'T':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(target_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'd': /* lightcurve output directory */
        case 'D':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(lightcurveDirectory,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'a': /* Handle all lightcurves */
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&allLightcurves);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode mosaic width %s\n",*argv);
              errorFlag = 1;
            }
            if (allLightcurves <= 2) {
              fprintf(stderr,"ERROR: value %d of -a is too small\n",allLightcurves);
              errorFlag = 1;
            }
          }
          break;


        case 'g': /* write good points only */
        case 'G': 
          goodFlag =  1;
          break;

        case 'v': /* verbose */
          verbose += 1;
          break;

        case 'i': /* Ignore limiting magnitudes */
          ignoreLimiting = 1;;
          break;


        case 'r': /* show rejection statistics */
          printReject  = 1;
          break;

        case 'V': /* even more verbose */
          verbose += 2;
          break;

        case 'f': /* Ignore defect filter */
        case 'F': 
          ignoreDefect = 1;
          break;


        default:
          fprintf(stderr,"ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */
  if (list_name[0] == 0) {
    fprintf(stderr,"ERROR: No list file was specified\n");
    errorFlag = 1;
  }
 

  if ((target_name[0] == 0) && (allLightcurves == 0)) {
    fprintf(stderr,"ERROR: No target catalog file was specified\n");
    errorFlag = 1;
  }

  if ((target_name[0] != 0) && (allLightcurves > 0)) {
    fprintf(stderr,"ERROR: Specify only one of -t or -a \n");
    errorFlag = 1;
  }
  

  if (lightcurveDirectory[0] == 0) {
    fprintf(stderr,"ERROR: No lightcurve directory was specified\n");
    errorFlag = 1;
  }

  
  if (allLightcurves == 0) {
    /* Open the target catalog file */
    target_handle = Open(target_name,"r");
    if (target_handle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to find the target catalog file %s\n",target_name);
    } else {
      if (verbose) {
        fprintf(stderr,"Found target catalog file %s\n",target_name);
      }
    }
  }

  /* Open the master photometry file */

  strcpy(master_name,ingestDirectory);
  strcat(master_name,"/master_all_");
  lastSlashPos = list_name;
  while ((slashPos = strstr(lastSlashPos,"/")) != NULL) {
    lastSlashPos = ++slashPos;
  }
  strcat(master_name,lastSlashPos);
  dotPos = strstr(master_name,".");
  if (dotPos != NULL) {
    *dotPos = 0;
  }
  if (qualifier[0] != 0) {
    strcat(master_name,qualifier);
  }
  strcat(master_name,".db");
  master_handle = Open(master_name,"r");
  if (master_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the master catalog file %s\n",master_name);
  } else {
    if (verbose) {
      fprintf(stderr,"Found master catalog file %s\n",master_name);
    }
  }

  /* Open the list file */
  list_handle = fopen(list_name,"rt");
  if (list_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the list file %s\n",list_name);
  } else {
    if (verbose) {
      fprintf(stderr,"Found list file %s\n",list_name);
    }
  }
  if (ignoreLimiting) {
    fprintf(stderr,"ERROR: IGNORE_LIMITING is set in extract_lightcurves.c\n");
  }
  

  catalogdir = getenv("DASCH_CATALOG");
  if (catalogdir == NULL) {
    fprintf(stderr,"DASCH_CATALOG is not defined\n");
    printf("0\n");
    return(-1);
  }
  strcpy(catalogname,catalogdir);
  if (strstr(qualifier,"kepler")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"kepler.dat");
  }

  if (strstr(qualifier,"apass")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"apass.dat");
  }

  if (strstr(qualifier,"atlas")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"atlas.dat");
  }

  if (strstr(qualifier,"gaia")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"gaiadr2.dat");
  }

  if (strstr(qualifier,"experimental")) {
    slashPtr = strrchr(catalogname,'/');
    if (slashPtr != NULL) {
      slashPtr++;
    } else {
      slashPtr = catalogname;
    }
    *slashPtr = 0;
    strcat(catalogname,"experimental.dat");
  }



  strcpy(indexname,catalogname);
  dotPtr = strstr(indexname,".");
  if (dotPtr != NULL) {
    *dotPtr = 0;
  }
  strcat(indexname,".idx");



  catalogHandle = Open(catalogname,"r");
  if (catalogHandle == NULL) {
    printf("ERROR Could not open catalog file %s\n",catalogname);
    errorFlag = 1;
  } 
  indexHandle = Open(indexname,"r");
  if (indexHandle == NULL) {
    printf("ERROR Could not open catalog index file %s\n",indexname);
    errorFlag = 1;
  } 



  if (errorFlag) {
    fprintf(stderr,"Usage: extract_lightcurves -l <plate list> \n");
    fprintf(stderr,"                           -t <target catalog file> \n");
    fprintf(stderr,"                           -a [npoints] extract all stars with at least this number of points \n");
    fprintf(stderr,"                           -d <lightcurve output directory> \n");
    fprintf(stderr,"                          [-g  write good stars only]\n");
    fprintf(stderr,"                          [-q <qualifier>] \n");
    fprintf(stderr,"                          [-r] print rejection statistics\n");
    fprintf(stderr,"                          [-b  <blend reference ID>] (Can specify up to ten times)\n");
    fprintf(stderr,"                          [-i]  do not plot limiting magnitudes \n");
    fprintf(stderr,"                          [-f]  ignore the defect bit \n");
    fprintf(stderr,"                          [-v]  verbose\n");
    fprintf(stderr,"                          [-V  even more verbose]\n");
    return(-1);
  }
  if (ignoreDefect) {
    fprintf(stderr,"ERROR: defect bits are being ignored\n");
  }


  /* Find out how many records are in the list file */
  if (list_handle != NULL) {
    while(1) {
      inBuffer = fgets(inLine,MAX_BUFFER,list_handle);
      if (inBuffer == NULL) {
        break;
      }
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
        if (lineLen > (MAX_PLATE_NAME-2)) {
          fprintf(stderr,"ERROR: MAX_PLATE_NAME %d exceeded with %d\n",MAX_PLATE_NAME-2,lineLen);
        }
        nlines++;
      }
    }
    fclose(list_handle);
    plateArray = (PPLATE)calloc(nlines,sizeof(PLATE));
    vector1 = (double*)calloc(nlines,sizeof(double));
    vector2 = (double*)calloc(nlines,sizeof(double));

    /* Now read in and store the plates */
    if (plateArray == NULL) {
      fprintf(stderr,"ERROR: allocation failed on plateArray\n");
      return(-1);
    }
    if (vector1 == NULL) {
      fprintf(stderr,"ERROR: allocation failed on vector1\n");
      return(-1);
    }
    if (vector2 == NULL) {
      fprintf(stderr,"ERROR: allocation failed on vector2\n");
      return(-1);
    }
    list_handle = fopen(list_name,"rt");
    if (list_handle == NULL) {
      fprintf(stderr,"ERROR: Failed to find the list file %s\n",list_name);
      return(-1);
    }


    while(1) {
      inBuffer = fgets(inLine,MAX_BUFFER,list_handle);
      if (inBuffer == NULL) {
        break;
      }
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
        if ((lineLen+2) > MAX_PLATE_NAME) {
          fprintf(stderr,"ERROR: lineLen %d exceeds maxLineLen %d\n",lineLen,MAX_PLATE_NAME);
          return(-1);
        }
        if (plateListLength >= nlines) {
          fprintf(stderr,"ERROR: plateListLength %d exceeds nlines %d\n",plateListLength,nlines);
          return(-1);
        }
        strcpy(plateArray[plateListLength].Plate,inBuffer);
        plateListLength++;
      }
    }
    fclose(list_handle);
    list_handle = NULL;
    if (verbose) {
      fprintf(stderr,"Read a total of %d lines from %s\n",plateListLength,list_name);
    }
  }
  /* Now access the database for these plates */
  if (InitPlateList(plateArray,plateListLength,ingestDirectory,catalogallDirectory,qualifier,verbose,ignoreLimiting) != 0) {
    return(-1);
  }
  if (printReject) {
    char goodBinsFileName[MAX_BUFFER];
    FILE *goodBinsFileHandle = NULL;
    char *seriesTable = (char *)calloc(MAX_SERIES_STRING*(MAX_SERIES+1),1);
    int numSeries = 0;
    int seriesIndex;
    int *spatialBinTable = NULL;
    char series[MAX_SERIES_STRING];
    int plateNumber;
    int mosaicNumber;
    int binning;
    int rotation;
    int spatial_bin;
    /* Print an overall summary of successful bins as a function of plate series */
    if (seriesTable == NULL) {
      fprintf(stderr,"ERROR: seriesTable allocation failure\n");
      exit(-1);
    }
    /* First figure out how many series we have */
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (ParseFilename(pPlate->Plate,series,&plateNumber,&mosaicNumber,&binning,&rotation) != 0) {
        for (seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
          if (strcmp(series,&seriesTable[seriesIndex*MAX_SERIES_STRING]) == 0) {
            break;
          }
        }
        if (seriesIndex >= numSeries) {
          strcpy(&seriesTable[numSeries*MAX_SERIES_STRING],series);
          numSeries++;
          if (numSeries >= MAX_SERIES) {
            fprintf(stderr,"ERROR: MAX_SERIES is too small\n");
            exit(-1);
          }
        }
      }
    }
    /* O.K. To allocate the bin table */
    strcpy(&seriesTable[numSeries*MAX_SERIES_STRING],"TOTAL");
    spatialBinTable = calloc((numSeries+1)*(MAX_SPATIAL_BINS+1),sizeof(int));
    if (spatialBinTable == NULL) {
      fprintf(stderr,"ERROR: failed to allocate spatialBinTable\n");
      exit(-1);
    }
    /* Now count up our successful bins */
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (ParseFilename(pPlate->Plate,series,&plateNumber,&mosaicNumber,&binning,&rotation) != 0) {
        for (seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
          if (strcmp(series,&seriesTable[seriesIndex*MAX_SERIES_STRING]) == 0) {
            break;
          }
        }
        if (seriesIndex < numSeries) {
          spatialBinTable[seriesIndex*(MAX_SPATIAL_BINS+1)]++; /* Entry total is the total count of plates in the series */
          spatialBinTable[numSeries*(MAX_SPATIAL_BINS+1)]++; /* Entry total is the total count of plates in the series */
          if (pPlate->spatial_bin_table != NULL) {
            for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
              if (spatial_bin == pPlate->spatial_bin_table[spatial_bin].spatial_bin) {
                spatialBinTable[(seriesIndex*(MAX_SPATIAL_BINS+1))+spatial_bin]++;
                spatialBinTable[(numSeries*(MAX_SPATIAL_BINS+1))+spatial_bin]++;
              }
            }
            
          }
        }
      }
    }
    

    strcpy(goodBinsFileName,lightcurveDirectory);
    strcat(goodBinsFileName,"/binsummary.txt");
    goodBinsFileHandle = fopen(goodBinsFileName,"wt");
    if (goodBinsFileHandle == NULL) {
      fprintf(stderr,"ERROR: Failed to open %s\n",goodBinsFileName);
      exit(-1);
    }

    fprintf(goodBinsFileHandle,"  Series Plates  bin1  bin2  bin3  bin4  bin5  bin6  bin7  bin8  bin9\n");
    for (seriesIndex = 0; seriesIndex <= numSeries; seriesIndex++) {
      fprintf(goodBinsFileHandle,"%8s ",&seriesTable[seriesIndex*MAX_SERIES_STRING]);
      for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
	fprintf(goodBinsFileHandle," %5d",spatialBinTable[(seriesIndex*(MAX_SPATIAL_BINS+1))+spatial_bin]);
      }
      fprintf(goodBinsFileHandle,"\n");
    }
    fclose(goodBinsFileHandle);
    
    free(spatialBinTable);
    free(seriesTable);

  }
  /* Now read in the target catalog */

  if (allLightcurves == 0) {
    target_header = table_header(target_handle,TABLE_PARSE);
    if (target_header == NULL) {
      fprintf(stderr,"ERROR: Failed to read header for %s\n",target_name);
      return(-1);
    }
    target_table = table_loadva(target_handle,
                                &target_header,
                                NULL, /* hbase */
                                NULL, /* rows */
                                NULL,
                                sizeof(TARGET),
                                &target_nrecs,
                                TblBuf,"REF"    ,TblOff(PTARGET,REF),MAX_REF,
                                TblBuf,"src_name"    ,TblOff(PTARGET,src_name),MAX_SRC_LENGTH,
                                TblDbl,"ra",TblOff(PTARGET,ra),
                                TblDbl,"dec",TblOff(PTARGET,dec),
                                0,"end",0);
    if (target_table == NULL) {
      fprintf(stderr,"ERROR: Failed to read table for %s\n",target_name);
      return(-1);
    }
    if (verbose) {
      fprintf(stderr,"read %d records for %s\n",target_nrecs,target_name);
    }
    for (target_index = 0; target_index < target_nrecs; target_index++) {
      pTarget = &target_table[target_index];
      pTarget->gsc_bin_index = 0;
      pTarget->color = 0;
      pTarget->class = 0;
      pTarget->VFlag = 0;
      pTarget->MAGFlag = 0;
      pTarget->RaPM = 0;
      pTarget->DecPM = 0;
     
    }


  } else {
    target_table = (PTARGET)calloc(1,sizeof(TARGET));
    if (target_table == NULL) {
      fprintf(stderr,"ERROR: failed to allocate target_table\n");
    }
  }
  /* Read the header of the master table */
  master_header = table_header(master_handle,TABLE_PARSE);
  if (master_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",master_name);
    return(-1);
  }


  master_descriptor = table_create_descrip(&master_nrecs,
                                           TblDbl,"magcal_local",TblOff(PMASTER,magcal_local),
                                           TblDbl,"Stdmag",TblOff(PMASTER,Stdmag),
                                           TblDbl,"Date",TblOff(PMASTER,Date),
                                           TblDbl,"ra",TblOff(PMASTER,ra),
                                           TblDbl,"dec",TblOff(PMASTER,dec),
                                           TblDbl,"limiting_mag_local",TblOff(PMASTER,limiting_mag_local),
                                           TblDbl,"magcal_local_rms",TblOff(PMASTER,magcal_local_rms),
                                           TblDbl,"magcal_iso",TblOff(PMASTER,magcal_iso),
                                           TblDbl,"magcal_iso_rms",TblOff(PMASTER,magcal_iso_rms),
                                           TblDbl,"magcal_local_error",TblOff(PMASTER,magcal_local_error),
                                           TblInt,"spatial_bin",TblOff(PMASTER,spatial_bin),
                                           TblInt,"npoints_local",TblOff(PMASTER,npoints_local),
                                           TblInt,"local_bin_index",TblOff(PMASTER,local_bin_index),
                                           TblInt,"AFLAGS",TblOff(PMASTER,AFLAGS),
                                           TblInt,"BFLAGS",TblOff(PMASTER,BFLAGS),
                                           TblBuf,"REF"    ,TblOff(PMASTER,REF),MAX_REF,
                                           TblBuf,"Plate"    ,TblOff(PMASTER,Plate),MAX_PLATE_NAME,
                                           TblInt,"gsc_bin_index",TblOff(PMASTER,gsc_bin_index),
                                           0,"end",0);

  if (master_descriptor == NULL) {
    fprintf(stderr,"ERROR: Failed to allocate descriptor for %s\n",master_name);
    return(-1);
  }
  table_loadmap(master_header,master_descriptor);


  if (verbose) {
    fprintf(stderr,"Reading %s\n",master_name);
  }
  target_index = 0;
  pTarget = &target_table[0];
  memset(pOldMaster,0,sizeof(MASTER));
  memset(pMaster,0,sizeof(MASTER));
  if (verbose > 1) {
    debugFlag = 1;
  }
  if(ReadMaster(pMaster,master_handle,master_header,master_descriptor,&master_row) == 0) {
    fprintf(stderr,"ERROR: Failed to read first record for %s\n",master_name);
    return(-1);
  }
  master_nrecs++;


  while(1) {
    int nextMaster = 0;
    compareFlag = strcmp(pTarget->REF,pMaster->REF);
    if (compareFlag == 0) {
      if ((allLightcurves > 0) && (strcmp(pMaster->REF,"NONE") == 0)) {
        match_count = 0;
      } else {
        if ((writePointResult = WritePoint(plateArray,plateListLength,pMaster,&debugFlag)) < 0) {
          return(-1);
        }
        match_count+= writePointResult;
      }
      nextMaster = 1;
    } else if (compareFlag < 0) {
      outputResult = WriteLightcurve(pGscCache,indexHandle,catalogHandle,plateArray,plateListLength,match_count,pTarget,lightcurveDirectory,vector1,vector2,verbose,goodFlag,printReject,&clipmed,ignoreDefect,ignoreLimiting,allLightcurves);
      SaveLightcurve(plateArray,plateListLength,pTarget,blendStrings,blendCount,blendArray,&savedBlendCount,clipmed);
      /* Re-initialize the plateArray for the next star */
      for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
        pPlate = &plateArray[plateIndex];
        pPlate->starDetected = 0;
        pPlate->validPoint = 0;
        pPlate->reject_reason1 = 0;
        pPlate->reject_reason2 = 0;
      }
      if (verbose > 1) {
        debugFlag = 1;
      }
      if (outputResult < 0) {
        return(-1);
      } else {
        outputPlotCount += outputResult;
        match_count = 0;
      }
      memcpy(pOldTarget,pTarget,sizeof(TARGET));
      if (allLightcurves > 0) {
        strcpy(pTarget->REF,pMaster->REF);
        strcpy(pTarget->src_name,pMaster->REF);
        pTarget->ra = pMaster->ra;
        pTarget->dec = pMaster->dec;
      } else {
        target_index++;
        if (target_index >= target_nrecs) {
          break;
        }
        pTarget = &target_table[target_index];
        if (strcmp(pOldTarget->REF,pTarget->REF) > 0) {
          fprintf(stderr,"ERROR: target file is out of order at %s %s\n",pOldTarget->REF,pTarget->REF);
          return(-1);
        }
      }
    } else {
      nextMaster = 1;
    }
    if (nextMaster != 0) {
      memcpy(pOldMaster,pMaster,sizeof(MASTER));
      if (ReadMaster(pMaster,master_handle,master_header,master_descriptor,&master_row) == 0) {
        break;
      }
      master_nrecs++;
      if (strcmp(pOldMaster->REF,pMaster->REF) > 0) {
        fprintf(stderr,"ERROR: master file is out of order at %s %s\n",pOldMaster->REF,pMaster->REF);
        return(-1);
      }
    }
  }
  if (strcmp(pOldTarget->REF,pTarget->REF) != 0) {
    outputResult = WriteLightcurve(pGscCache,indexHandle,catalogHandle,plateArray,plateListLength,match_count,pTarget,lightcurveDirectory,vector1,vector2,verbose,goodFlag,printReject,&clipmed,ignoreDefect,ignoreLimiting,allLightcurves);
    SaveLightcurve(plateArray,plateListLength,pTarget,blendStrings,blendCount,blendArray,&savedBlendCount,clipmed);
    /* Re-initialize the plateArray for the next star */
    for (plateIndex = 0; plateIndex < plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      pPlate->starDetected = 0;
      pPlate->validPoint = 0;
      pPlate->reject_reason1 = 0;
      pPlate->reject_reason2 = 0;
    }
    if (outputResult < 0) {
      return(-1);
    } else {
      outputPlotCount += outputResult;
      match_count = 0;
    }
  }
  if (verbose) {
    fprintf(stderr,"read %d records for %s\n",master_nrecs,master_name);
  }
  /* If we need to produce a blend graph, do it now */
  if (blendCount > 0) {
    if (blendCount == savedBlendCount) {
      match_count = 0;
      MergeBlendArray(blendStrings,blendCount,blendArray,&savedBlendCount,&blendTarget,plateListLength);
      outputResult = WriteLightcurve(pGscCache,indexHandle,catalogHandle,blendArray[0],plateListLength,match_count,&blendTarget,lightcurveDirectory,vector1,vector2,verbose,goodFlag,printReject,&clipmed,ignoreDefect,ignoreLimiting,allLightcurves);
    
      if (outputResult < 0) {
        return(-1);
      } else {
        outputPlotCount += outputResult;
        match_count = 0;
      }
    } else {
      fprintf(stderr,"ERROR: blendCount %d does not equal found blend objects %d\n",blendCount,savedBlendCount);
    }
  }




  /* All done.  Clean up */
 
  if (master_header != NULL) {
    table_hdrfree(master_header);
  }
  if (master_handle != NULL) {
    Close(master_handle);
  }
  if (master_descriptor != NULL) {
    Free(master_descriptor);
  }
  if (master_row != NULL) {
    table_rowfree(master_row);
  }

  if (target_table != NULL) {
    free(target_table);
  }
  if (target_header != NULL) {
    table_hdrfree(target_header);
  }
  if (target_handle != NULL) {
    Close(target_handle);
  }

  if (list_handle != NULL) {
    fclose(list_handle);
  }

  if(plateArray != NULL) {
    for (plateIndex = 0; plateIndex <  plateListLength; plateIndex++) {
      pPlate = &plateArray[plateIndex];
      if (pPlate->spatial_bin_table != NULL) {
        free(pPlate->spatial_bin_table);
      }
      if (pPlate->dmagcor_table != NULL) {
        Free(pPlate->dmagcor_table);
      }
      if (pPlate->extinction_table != NULL) {
        Free(pPlate->extinction_table);
      }

    }
    free(plateArray);
  }
  if (vector1 != NULL) {
    free(vector1);
  }
  if (vector2 != NULL) {
    free(vector2);
  }
  if (rmsHandle != NULL) {
    fclose(rmsHandle);
  }
  for (blendIndex = 0; blendIndex < savedBlendCount; blendIndex++) {
    free(blendArray[blendIndex]);
    blendArray[blendIndex] = NULL;
 
  }

  time(&curTime);
  curTime -= startTime;
  fprintf(stdout,"Output stars %d seconds %d\n",
          outputPlotCount,curTime);
  Close(catalogHandle);
  Close(indexHandle);
    


  return(0);
}

