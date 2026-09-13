// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* web_data.c - support display of multiple object data files
 *
 * -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm  -L/usr/lib/mysql -L/usr/lib64/mysql  -l mysqlclient  web_data.c pipelineutils.a -ltable -lutil -lgd -lwcs -o web_data
 *
 *
 * Oct  7, 2013 Edward J. Los - Initial version
 * Jan 20, 2015 Edward J. Los - V6 data format: add A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex
 * Feb 16, 2014 Edward J. Los - Add timeAccuracy to the short form lightcurve table
 * May 16, 2015 Edward J. Los - List the bits in AFLAGS for the short output formats
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * Feb 22, 2016 Edward J. Los - In formatAflags, split AFLAGS into two columns: AFLAGS and AFLAGSBits
 * Feb 24, 2016 Edward J. Los - Change "formatAflags" to "FormatFlagsBits" and add BFLAGS and quality bitmaps
 * Feb 29, 2016 Edward J. Los - Add a source qualifier to for accurate quality bitmaps.
 * Mar  8, 2016 Edward J. Los - Use two separate fields for FormatFlagsBits
 * Jan 23, 2018 Edward J. Los   Support the merged experimental table: add catalogNumber to WriteStarbaseRecord.
 * Aug 13, 2018 Edward J. Los   Define enableRematch to optimize location of transients (Redefine -O) to support this function
 * Dec 10, 2018 Edward J. Los   Add NUMBER to the short plot for positive identification of the image
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <math.h>

#include <table.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include "photometryutils.h"
#include "searchgsc.h"
#include "daschunistd.h"

extern FLAGSENTRY AflagsTable[];
extern int AflagsTableSize;
extern QUALITYBIT webQualityMasks[];
extern int webQualityTableSize;
double fmax(double x, double y);
double fmin(double x, double y);

/* Make GCC happy about potential overflows */
#define BUFPRINTF(BUF, FMT, ARGS...) \
  do { \
    int rv = snprintf((BUF), sizeof(BUF), (FMT), ARGS); \
    if (rv < 0 || rv >= sizeof(BUF)) \
      abort(); \
  } while (0)

#define MAX_BUFFER 500
#define MAX_FILENAME 512
#define GSC_ALLOC_INCREMENT 1000
#define LIST_STRING_ALLOC 2000
#define LIST_SUBSTR_SIZE 28
#define MAX_MOSAICID_STRING 50
#define ERROR_BAR_CLIP 0.4
/* #define dasch_proto 1 Experimental version enabled */
/* #define ACCEPT_LOW_ALTITUDE */
/* #define DEBUG_LIMITING 1 */

#define PLOT_AVE	0x01
#define PLOT_FLIP	0x02
#define PLOT_LINE	0x04
#define PLOT_EXACT	0x08
#define PLOT_UPPER_L	0x10
#define MAX_LABEL 64
#define DOTSZ 2
#define UPPER_LIMIT 19.0;
#define USE_BUTTONS 1

#define COLOR_BLUE 1  /* too bright */
#define COLOR_RED  2  /* second quality */
#define COLOR_GRAY 5  /* undected */

/* #define LOS_DEBUG */
#define LOS_NUMBER 125406
#define LOS_PLATENUMBER  5289
#define LOS_SERIESID 22
#define LOS_GSC_BIN 113776398
#define LOS_REF_NUMBER  309313291201748L
#define LOS_REF  "DASCH_J093132.9+201748"

typedef struct _object_table {
    long long REFNumber;
    double cra;   /* Search center */
    double cdec;  /* Search center */
    double drad;  /* Distance in arcsec from the search center */
    double magcal_magdep; /* Average magnitude of the lightcurve */
    double ra;    /* object center */
    double dec;   /* object center */
    int npoints;  /* Total points in the lightcurve */
    int nplot;    /* Plotted points in the lightcurve */
    int refType;
    int pointsfound; /* Number of points found */
    int written;
    char src_name[MAX_SRC_LENGTH];
    char nearbyObjects[MAX_NEARBY_OBJECTS_STRING];
    char REF[MAX_REF];
} OBJECTTABLE,*POBJECTTABLE;

/* Extension used by the web-based plotter */
typedef struct _starimageext {
  FILESTARIMAGE filestarimage;
  PHOTPLATES photplates;
  int mosaicNumber;
  int quality;
  int plateVersionId;
    int object_index;
    char REFNumberC[2*MAX_REF]; /* Encoded catalog reference number as an ASCII number */
} STARIMAGEEXT,*PSTARIMAGEEXT;

extern char *catalogText[MAX_CATALOG_NUMBER];
int gif=1;

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
double str2dec(		        /* Return Dec in degrees from string */
    const char* in);	/* Character string (dd:mm:ss.sss or dd.dddd) */
double str2ra(		        /* Return RA in degrees from string */
    const char* in);	/* Character string (hh:mm:ss.sss or dd.dddd) */

double wcsdist(	/* Compute angular distance between 2 sky positions */
               double ra1,	/* First longitude/right ascension in degrees */
               double dec1,	/* First latitude/declination in degrees */
               double ra2,	/* Second longitude/right ascension in degrees */
               double dec2);	/* Second latitude/declination in degrees */


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
char *webbuff(	/* Read URL into buffer across the web */
          char *url,	/* URL to read */
          int diag,	/* 1 to print diagnostic messages */
          int *lbuff);	/* Length of buffer (returned) */


typedef struct sortTable {
  double julianDate; /* This item must be first for sorting */
  char listString[LIST_SUBSTR_SIZE+1]; /* Pointer to the list string in tclass */
} SORTTABLE,*PSORTTABLE;




extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
int binSearchCount = 0;
int cacheHits = 0;
int cacheMisses = 0;
int votablerun = 0;

SERIESLIST seriesList[MAX_SERIES+1];
int seriesCount = 0;
int allSeriesMask = 0;
int OutputTableCompare(const void *first, const void *second)
{
    PSTARIMAGEEXT pImageFirst = (PSTARIMAGEEXT)first;
    PSTARIMAGEEXT pImageSecond = (PSTARIMAGEEXT)second;
    int result;
    if (pImageFirst->object_index > pImageSecond->object_index) {
        return(1);
    } else if (pImageFirst->object_index < pImageSecond->object_index) {
        return(-1);
    } else {
        result = strcmp(pImageFirst->REFNumberC,pImageSecond->REFNumberC);
        if (result != 0) {
            return(result);
        } else {
            if (pImageFirst->filestarimage.Date > pImageSecond->filestarimage.Date) {
                return(1);
            } else if (pImageFirst->filestarimage.Date < pImageSecond->filestarimage.Date) {
                return(-1);
            } else {
                return(0);
            }
        }
    }
}
void WriteOutputLinks(char *description,char *lsbin,char *binaries,char* short_name,char *long_name,char *b_name,char *c_name,char *d_name, char * f_name,char* g_name,char *h_name)
{
  int statResult;
  struct stat statbuf;
    char *slashPtr;
    char *charPtr;
    char prezip_name[MAX_FILENAME];
  char postzip_name[MAX_FILENAME];
  char zip_name[MAX_FILENAME];
  char txt_name[MAX_FILENAME];
    char votable_name[MAX_FILENAME];
    char votable_zip_name[MAX_FILENAME];
  char cmdStr[MAX_BUFFER];
    char *dotPtr;
    /*
   *   Short table: List in the browser window
     */
    printf("<br /><b>%s</b> ",description);
    statResult = stat(short_name,&statbuf);
    if (statResult != 0) {
        printf("ERROR: short file %s does not exist\n",short_name);
        exit(-1);
    }
    slashPtr = strrchr(short_name,'/');
    if (slashPtr == NULL) {
        slashPtr = short_name;
    } else {
        slashPtr++;
    }
    charPtr = strstr(short_name,"/tmp");
    if (charPtr == NULL) {
        charPtr = short_name;
    } else {
        charPtr++;
    }
    printf("<br />A: <a href=\"%s\">%s</a> ",charPtr,slashPtr);
    /*
   *   Short Table: Download db as a zip file
     */
    strcpy(zip_name,short_name);
    dotPtr = strstr(zip_name,".db");
    if (dotPtr != NULL) {
        *dotPtr = 0;
    }
  strcpy(txt_name,zip_name);
    strcat(zip_name,".db.gz");
  strcat(txt_name,".txt");
    strcpy(prezip_name,short_name);
    dotPtr = strstr(prezip_name,".db");
    if (dotPtr != NULL) {
        *dotPtr = 0;
    }
    strcat(prezip_name,"_zip.db");
  strcpy(postzip_name,prezip_name);
  strcat(postzip_name,".gz");


    statResult = stat(zip_name,&statbuf);
    if (statResult == 0) {
        unlink(zip_name);
    }
    statResult = stat(postzip_name,&statbuf);
    if (statResult == 0) {
        unlink(postzip_name);
    }
    sprintf(cmdStr,"%s/cp %s %s\n",lsbin,short_name,prezip_name);
    system(cmdStr);
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,short_name);
    system(cmdStr);

    slashPtr = strrchr(txt_name,'/');
    if (slashPtr == NULL) {
        slashPtr = txt_name;
    } else {
        slashPtr++;
    }
    charPtr = strstr(txt_name,"/tmp");
    if (charPtr == NULL) {
        charPtr = txt_name;
    } else {
        charPtr++;
    }
    printf("<br />B: <a href=\"%s\">%s</a> ",charPtr,slashPtr);
    strcpy(b_name,slashPtr); /* short db */

    slashPtr = strrchr(zip_name,'/');
    sprintf(cmdStr,"%s/mv %s %s\n",lsbin,prezip_name,short_name);
    system(cmdStr);
    if (slashPtr == NULL) {
        slashPtr = zip_name;
    } else {
        slashPtr++;
    }
    charPtr = strstr(zip_name,"/tmp");
    if (charPtr == NULL) {
        charPtr = zip_name;
    } else {
        charPtr++;
    }
    printf("<br />C: <a href=\"%s\">%s</a> ",charPtr,slashPtr);
    strcpy(c_name,slashPtr); /* short db */
    /*
     * Short table download as a zipped VOTable file
     */
    strcpy(votable_name,short_name);
    dotPtr = strstr(votable_name,".db");
    if (dotPtr != NULL) {
        *dotPtr = 0;
    }
    strcat(votable_name,".xml");
    strcpy(votable_zip_name,votable_name);
    strcat(votable_zip_name,".gz");

    statResult = stat(votable_name,&statbuf);
    if (statResult == 0) {
        unlink(votable_name);
    }
    statResult = stat(votable_zip_name,&statbuf);
    if (statResult == 0) {
        unlink(votable_zip_name);
    }
    sprintf(cmdStr,"%s/votable -i %s -o %s\n",binaries,short_name,votable_name);
    system(cmdStr);
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,votable_name);
    system(cmdStr);
    slashPtr = strrchr(votable_zip_name,'/');
    if (slashPtr == NULL) {
        slashPtr = votable_zip_name;
    } else {
        slashPtr++;
    }
    charPtr = strstr(votable_zip_name,"/tmp");
    if (charPtr == NULL) {
        charPtr = votable_zip_name;
    } else {
        charPtr++;
    }
    printf("<br />D: <a href=\"%s\">%s</a> ",charPtr,slashPtr);
     strcpy(d_name,slashPtr); /* short xml */

    /*
   *   Full table: show in browser window
     */
    statResult = stat(long_name,&statbuf);
    if (statResult != 0) {
        printf("ERROR: short file %s does not exist\n",long_name);
        exit(-1);
    }
    slashPtr = strrchr(long_name,'/');
    if (slashPtr == NULL) {
        slashPtr = long_name;
    } else {
        slashPtr++;
    }
    charPtr = strstr(long_name,"/tmp");
    if (charPtr == NULL) {
        charPtr = long_name;
    } else {
        charPtr++;
    }
    printf("<br />E: <a href=\"%s\">%s</a> ",charPtr,slashPtr);
    /*
   *   Full table: download db as a zip file
     */
    strcpy(zip_name,long_name);
  strcat(zip_name,".gz");
    strcpy(prezip_name,long_name);
    dotPtr = strstr(prezip_name,".db");
    if (dotPtr != NULL) {
        *dotPtr = 0;
    }
  strcpy(txt_name,prezip_name);
    strcat(prezip_name,"_zip.db");
  strcat(txt_name,".txt");
    statResult = stat(zip_name,&statbuf);
    if (statResult == 0) {
        unlink(zip_name);
    }
  strcpy(postzip_name,prezip_name);
  strcat(postzip_name,".gz");
    statResult = stat(postzip_name,&statbuf);
    if (statResult == 0) {
        unlink(postzip_name);
    }

    sprintf(cmdStr,"%s/cp %s %s\n",lsbin,long_name,prezip_name);
    system(cmdStr);
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,long_name);
    system(cmdStr);
    sprintf(cmdStr,"%s/mv %s %s\n",lsbin,prezip_name,long_name);
    system(cmdStr);

    slashPtr = strrchr(txt_name,'/');
    if (slashPtr == NULL) {
        slashPtr = txt_name;
    } else {
        slashPtr++;
    }
    charPtr = strstr(txt_name,"/tmp");
    if (charPtr == NULL) {
        charPtr = txt_name;
    } else {
        charPtr++;
    }
    printf("<br />F: <a href=\"%s\">%s</a> ",charPtr,slashPtr);
     strcpy(f_name,slashPtr); /* long db */

    slashPtr = strrchr(zip_name,'/');
    if (slashPtr == NULL) {
        slashPtr = zip_name;
    } else {
        slashPtr++;
    }
    charPtr = strstr(zip_name,"/tmp");
    if (charPtr == NULL) {
        charPtr = zip_name;
    } else {
        charPtr++;
    }
    printf("<br />G: <a href=\"%s\">%s</a> ",charPtr,slashPtr);
     strcpy(g_name,slashPtr); /* long db */

    /*
     * full table download as a zipped VOTable file
     */
    strcpy(votable_name,long_name);
    dotPtr = strstr(votable_name,".db");
    if (dotPtr != NULL) {
        *dotPtr = 0;
    }
    strcat(votable_name,".xml");
    strcpy(votable_zip_name,votable_name);
    strcat(votable_zip_name,".gz");

    statResult = stat(votable_name,&statbuf);
    if (statResult == 0) {
        unlink(votable_name);
    }
    statResult = stat(votable_zip_name,&statbuf);
    if (statResult == 0) {
        unlink(votable_zip_name);
    }
    sprintf(cmdStr,"%s/votable -i %s -o %s\n",binaries,long_name,votable_name);
    system(cmdStr);
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,votable_name);
    system(cmdStr);
    slashPtr = strrchr(votable_zip_name,'/');
    if (slashPtr == NULL) {
        slashPtr = votable_zip_name;
    } else {
        slashPtr++;
    }
    charPtr = strstr(votable_zip_name,"/tmp");
    if (charPtr == NULL) {
        charPtr = votable_zip_name;
    } else {
        charPtr++;
    }
    printf("<br />H: <a href=\"%s\">%s</a> ",charPtr,slashPtr);
    strcpy(h_name,slashPtr); /* long xml */
    printf("<br />\n");
}

/* Wrapper for GetPhotPlate() */
void GetPhotPlateExt(MYSQL *pPhotConnection,PSTARIMAGEEXT pFileStarImageExt,char *catalogString)
{
  PPHOTPLATES pPhotPlates;
  PFILESTARIMAGE pFileStarImage = &pFileStarImageExt->filestarimage;
  int quality;
  pPhotPlates = &pFileStarImageExt->photplates;
  memset(pPhotPlates,0,sizeof(PHOTPLATES));
  pPhotPlates->mosaicNumber = 99;
  pPhotPlates->quality = QUALITY_UNINITIALIZED;
  if (GetPhotPlate(pPhotConnection,
                   GetSeriesString(pFileStarImage->seriesId,1),
                   pFileStarImage->plateNumber,
                   pPhotPlates,
                   catalogString,1) != 0) {
    printf("ERROR: Failed to get photplates record for %s%05d\n",GetSeriesString(pFileStarImage->seriesId,0),pFileStarImage->plateNumber);
    pPhotPlates->mosaicNumber=99;
    pPhotPlates->quality = QUALITY_UNINITIALIZED;
    pPhotPlates->versionId = 0;
  }
  if (pPhotPlates->versionId == 0) {
#if 0
    printf("ERROR: Null photplates record for %s%05d\n",GetSeriesString(pFileStarImage->seriesId,0),pFileStarImage->plateNumber);
    printf("SELECT * from photplates where seriesId = %d and plateNumber = %d;\n",pFileStarImage->seriesId,pFileStarImage->plateNumber);
#endif
    pPhotPlates->mosaicNumber=99;
    pPhotPlates->quality = QUALITY_UNINITIALIZED;
    pPhotPlates->versionId = 0;

  }
#ifdef DEBUG_BADSTALEID
  pFileStarImage->versionId = pPhotPlates->versionId;
#endif /* DEBUG_BADSTALEID */
  if (pPhotPlates->versionId != pFileStarImage->versionId) {
#if 0
    printf("ERROR: versionId mismatch %d %d for %s%05d\n",pPhotPlates->versionId,pFileStarImage->versionId,GetSeriesString(pFileStarImage->seriesId,0),pFileStarImage->plateNumber);
#else
    printf("*");
#endif
#if 0
    printf("SELECT * from photplates where seriesId = %d and plateNumber = %d;\n",pFileStarImage->seriesId,pFileStarImage->plateNumber);
#endif

    pPhotPlates->mosaicNumber=99;
    pPhotPlates->quality = QUALITY_UNINITIALIZED;
    pPhotPlates->versionId = 0;
  }
  GetFullQuality(&pFileStarImageExt->filestarimage,pPhotPlates->quality,&quality);
  pPhotPlates->quality = quality;
  pFileStarImageExt->mosaicNumber = pPhotPlates->mosaicNumber;
  pFileStarImageExt->quality = pPhotPlates->quality;
  pFileStarImageExt->plateVersionId = pPhotPlates->versionId;


}


int main(int argc,char *argv[])
{

  int nvals;
  char *argstr;
  char cmdchar;
  char REF[MAX_REF];          /* GSC2.3.2 reference number */
  int errorFlag = 0;
  char *tmpdir = NULL;
  char source[MAX_BUFFER];
  char binaries[MAX_BUFFER];
  char lsbin[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
    char dbfilename[MAX_BUFFER];
  char textfilename[MAX_BUFFER]; /* -t qualifier */
  char origtmpdir[MAX_BUFFER];
  char *charPtr;
    char object_name[MAX_FILENAME];
  char db_name[MAX_FILENAME];
  char zipdb_name[MAX_FILENAME];
  char zipdb_name_gz[MAX_FILENAME];
  char gz_name[MAX_FILENAME];
  char vo_name[MAX_FILENAME];
  char short_db_name[MAX_FILENAME];
  char short_txt_name[MAX_FILENAME];
  char short_zipdb_name[MAX_FILENAME];
  char short_zipdb_name_gz[MAX_FILENAME];
  char short_gz_name[MAX_FILENAME];
  char short_vo_name[MAX_FILENAME];
  char short_vo_name_gz[MAX_FILENAME];
    char output_db[MAX_FILENAME];
  char output_txt_db[MAX_FILENAME];
    char short_tar_db_name[MAX_FILENAME];
    char short_tar_txt_name[MAX_FILENAME];
    char tar_db_name[MAX_FILENAME];
    char tar_txt_name[MAX_FILENAME];
    char short_tar_xml_name[MAX_FILENAME];
    char tar_xml_name[MAX_FILENAME];
    FILE *outputhandle = NULL;
    FILE *outputtxthandle = NULL;
    FILE *shortoutputhandle = NULL;
    FILE *shorttxtoutputhandle = NULL;
    FILE *objectoutputhandle = NULL;
    FILE *objecttxtoutputhandle = NULL;
    FILE *shortobjectoutputhandle = NULL;
    FILE *shorttxtobjectoutputhandle = NULL;
    char objectrootname[MAX_FILENAME];
    char objectdescription[MAX_FILENAME];
    char objectfilename[MAX_FILENAME];
    char objecttxtfilename[MAX_FILENAME];
    char shortobjectfilename[MAX_FILENAME];
    char shorttxtobjectfilename[MAX_FILENAME];
    char b_name[MAX_FILENAME]; /* short txt */
    char c_name[MAX_FILENAME]; /* short db */
    char d_name[MAX_FILENAME]; /* short xml */
    char f_name[MAX_FILENAME]; /* full txt */
    char g_name[MAX_FILENAME]; /* full db */
     char h_name[MAX_FILENAME]; /* full xml */
    char* full_b_name = NULL;
    char* full_c_name = NULL;
    char* full_d_name = NULL;
    char* full_f_name = NULL;
    char* full_g_name = NULL;
    char* full_h_name = NULL;
  char starbase_title[MAX_BUFFER];
  char starbase_txt_title[MAX_BUFFER];

  File db_handle = NULL;
  TableHead db_header = NULL;
  PSTARIMAGEEXT db_table = NULL;
  PSTARIMAGEEXT pStarImageExt;
    PSTARIMAGEEXT pOutputImageExt;
    PSTARIMAGEEXT pNextOutputImageExt;
  PSTARIMAGEEXT pPrevOutputImageExt;
    PSTARIMAGEEXT output_table = NULL;
    size_t output_nrecs = 0;
    size_t output_index;
    size_t output_index2;

  size_t db_nrecs = 0;

  File object_handle = NULL;
  TableHead object_header = NULL;
  POBJECTTABLE object_table = NULL;
    POBJECTTABLE pObject = NULL;
  size_t object_nrecs = 0;
    size_t object_index = 0;

  int *db_indices = NULL;
  int *short_indices = NULL;

  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon = &fileCommon;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  int gotAnswer;
  int verbose = 0;
  int m44release = 0;
  char *dotloc;
  int statResult;
  struct stat statbuf;
  int numMagnitudes;
  int magnitudeIndex;
    PPHOTPLATES pPhotPlates;
  PPHOTPLATES pPrevPhotPlates;
  int dbCallocFlag = 0;
  char *listString = NULL;
  PSORTTABLE pSortTable = NULL;
  char* slashPtr;
  GALAXYCOMMON galaxycommon;
  PGALAXYCOMMON pGalaxyCommon = &galaxycommon;

    int writeObjectHeader = 1;
    int writeTxtObjectHeader = 1;
    int writeOutputHeader = 1;
    int writeTxtOutputHeader = 1;
  char AFLAGSBuffer[MAX_BITMAP_SIZE];
  char BFLAGSBuffer[MAX_BITMAP_SIZE];
  char qualityBuffer[MAX_BITMAP_SIZE];
  char AFLAGSbitsBuffer[MAX_BITMAP_SIZE];
  char BFLAGSbitsBuffer[MAX_BITMAP_SIZE];
  char qualitybitsBuffer[MAX_BITMAP_SIZE];
  char catalogString[MAX_BUFFER];
  int catalogNumber = 0;
  int tmpCatalogNumber;
  int enableRematch = 0;

  memset(seriesList,0,sizeof(seriesList));


#ifdef dasch_proto
  printf("dasch_proto is active\n");
#endif /* dasch_proto */

  memset(pGalaxyCommon,0,sizeof(GALAXYCOMMON));
  SetQueryCount(0);
  qualifier[0] = 0;
  source[0] = 0;
  dbfilename[0] = 0;
  binaries[0] = 0;
  lsbin[0] = 0;

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d\n",argstr,argc);

    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {


        case 'd': /* Temporary directory */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            tmpdir = *++argv;
          }
          break;

        case 'e': /* Executable binary image directory */
        case 'E':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: binaries directory length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(binaries,*argv);
            }
          }
          break;

        case 'l': /* gzip binary image directory */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: lsbin directory length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(lsbin,*argv);
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
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              strcpy(source,*argv);
              if (catalogNumber > 0) {
                sprintf(catalogString,"%d",catalogNumber);
                strcpy(qualifier,*argv);
              }
            }
          }
          break;

        case 'v': /* verbose */
          verbose += 1;
          break;

        case 'o': /* M44 Release */
          m44release += 1;
          break;

        case 'O':
          enableRematch = 1;
          break;




        case 'r': /* Object designation */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_REF-1)) {
              printf("ERROR: designation length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(REF,*argv);

            }
          }
          break;





        case 't': /* Text file name */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: text file name length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(dbfilename,*argv);
            }
          }
          break;







        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }

      }

    }
  }
  if (strlen(dbfilename) == 0) {
    printf("ERROR: dbfilename is not specified\n");
    errorFlag = 1;
  }
  if (tmpdir == NULL) {
    printf("ERROR: plot directory is not specified\n");
    errorFlag = 1;
  }

  if (strlen(binaries) == 0) {
    printf("ERROR: binaries directory is not specified\n");
    errorFlag = 1;
  }

  if (strlen(lsbin) == 0) {
    printf("ERROR: lsbin directory is not specified\n");
    errorFlag = 1;
  }

  if (strlen(dbfilename)+strlen(tmpdir)+5 > MAX_BUFFER) {
    printf("ERROR: full name of the text file is too long \n");
    errorFlag = 1;
  }



  /* We are done with the requested plot size.  If the user checks on a box, then we need to rescale the Y-axis.  There is no need to rescale the X-axis, but make sure we
     display the full range */











  if (errorFlag == 1) {
    printf("Usage: web_data -r <REF> -q <catalog>  -d <directory> -t <text name>  [-v]\n");
    printf("       where -v is the verbose flag\n");
    printf("             -r is the object reference\n");
    printf("             -d is the plot and text directory\n");
    printf("             -t is the text file name\n");
    printf("             -p is the plot file name\n");
    printf("             -e is the directory containing the executable binaries\n");
    printf("             -l is the directory containing the gzip binary\n");
    printf("             -O optimize location of transients\n");

    exit(-1);
  }

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pConnection);
  dasch_init_photdb(pPhotConnection);

  InitSeriesTable(pConnection,pPhotConnection);
  gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
  if (gotAnswer != 1)  {
    printf("ERROR: failed to get the global photometry table\n");
    exit(-1);
  }

  if (pPhotGlobal->magnitudeFile != PHOT_MAGNITUDEFILE_YES) {
    printf("ERROR: database uses obsolute magnitude table\n");
    exit(-1);
  }
  InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);
  if (enableRematch) {
    pFileCommon->enableRematch = 1;
  }


  /* Now see if we have already written out a result file for this object */
    slashPtr = strrchr(dbfilename,'/');
    if (slashPtr == NULL) {
        strcpy(textfilename,dbfilename);
    } else{
        slashPtr++;
        strcpy(textfilename,slashPtr);
    }
  strcpy(db_name,tmpdir);
  strcat(db_name,"/");
  strcat(db_name,textfilename);

  charPtr = strstr(tmpdir,"tmp/");
  if (charPtr == NULL) {
    printf("ERROR: invalid format for tmpdir %s\n",tmpdir);
    exit(-1);
  }
  charPtr += 4;
  strcpy(origtmpdir,charPtr);

  strcpy(object_name,tmpdir);
  strcat(object_name,"/table_");
  strcat(object_name,origtmpdir);
    strcat(object_name,".db");

  strcpy(output_db,tmpdir);
  strcat(output_db,"/output_");
  strcat(output_db,origtmpdir);
  strcpy(output_txt_db,output_db);
    strcat(output_db,".db");
  strcat(output_txt_db,".txt");

  strcpy(short_db_name,tmpdir);
  strcat(short_db_name,"/short_output_");
  strcat(short_db_name,origtmpdir);
  strcpy(short_txt_name,short_db_name);
    strcat(short_db_name,".db");
  strcat(short_txt_name,".txt");

  strcpy(short_tar_db_name,tmpdir);
  strcat(short_tar_db_name,"/short_all_db_");
  strcat(short_tar_db_name,origtmpdir);
    strcat(short_tar_db_name,".tar");


  strcpy(short_tar_txt_name,tmpdir);
  strcat(short_tar_txt_name,"/short_all_txt_");
  strcat(short_tar_txt_name,origtmpdir);
    strcat(short_tar_txt_name,".tar");

  strcpy(short_tar_xml_name,tmpdir);
  strcat(short_tar_xml_name,"/short_all_xml_");
  strcat(short_tar_xml_name,origtmpdir);
    strcat(short_tar_xml_name,".tar");

  strcpy(tar_db_name,tmpdir);
  strcat(tar_db_name,"/all_db_");
  strcat(tar_db_name,origtmpdir);
    strcat(tar_db_name,".tar");

  strcpy(tar_txt_name,tmpdir);
  strcat(tar_txt_name,"/all_txt_");
  strcat(tar_txt_name,origtmpdir);
    strcat(tar_txt_name,".tar");

  strcpy(tar_xml_name,tmpdir);
  strcat(tar_xml_name,"/all_xml_");
  strcat(tar_xml_name,origtmpdir);
    strcat(tar_xml_name,".tar");



  dotloc = strrchr(db_name,'.');
  if (dotloc) {
    *dotloc = 0;
    strcpy(vo_name,db_name);
    strcpy(zipdb_name,db_name);
    strcat(db_name,".db");
    strcat(vo_name,".xml");
    strcat(zipdb_name,"_zip.db");
    strcpy(zipdb_name_gz,zipdb_name);
    strcat(zipdb_name_gz,".gz");
    strcpy(gz_name,db_name);
    strcat(gz_name,".gz");
  } else {
    printf("ERROR: %s has no suffix\n",textfilename);
    exit(-1);
  }

  dotloc = strrchr(short_db_name,'.');
  if (dotloc) {
    *dotloc = 0;
    strcpy(short_vo_name,short_db_name);
    strcpy(short_zipdb_name,short_db_name);
    strcat(short_db_name,".db");
    strcat(short_vo_name,".xml");
        strcpy(short_vo_name_gz,short_vo_name);
        strcat(short_vo_name_gz,".gz");
    strcat(short_zipdb_name,"_zip.db");
        strcpy(short_zipdb_name_gz,short_zipdb_name);
        strcat(short_zipdb_name_gz,".gz");
    strcpy(short_gz_name,short_db_name);
    strcat(short_gz_name,".gz");
  } else {
    printf("ERROR: %s has no suffix\n",textfilename);
    exit(-1);
  }
  /* Clear out files from a previous invocation */
  statResult = stat(short_db_name,&statbuf);
  if (statResult == 0) {
    unlink(short_db_name);
  }
  statResult = stat(short_txt_name,&statbuf);
  if (statResult == 0) {
    unlink(short_txt_name);
  }
  statResult = stat(short_zipdb_name,&statbuf);
  if (statResult == 0) {
    unlink(short_zipdb_name);
  }
  statResult = stat(short_zipdb_name_gz,&statbuf);
  if (statResult == 0) {
    unlink(short_zipdb_name_gz);
  }
  statResult = stat(short_gz_name,&statbuf);
  if (statResult == 0) {
    unlink(short_gz_name);
  }
  statResult = stat(short_vo_name,&statbuf);
  if (statResult == 0) {
    unlink(short_vo_name);
  }
  statResult = stat(short_vo_name_gz,&statbuf);
  if (statResult == 0) {
    unlink(short_vo_name_gz);
  }


  statResult = stat(db_name,&statbuf);
  if (statResult != 0) {
        printf("ERROR: failed to find file %s\n",db_name);
        exit(-1);
  }
    /* Here we have already prepared the data.  Read it in */
    db_handle = Open(db_name,"r");
    if (db_handle == NULL) {
        errorFlag = 1;
        printf("ERROR: Failed to find the db file %s\n",db_name);
    } else {
        if (verbose) {
            printf("Found db file %s\n",db_name);
        }
    }


    db_table = table_loadva(db_handle,
                                                    &db_header,
                                                    NULL, /* hbase */
                                                    NULL, /* rows */
                                                    NULL,
                                                    sizeof(STARIMAGEEXT),
                                                    &db_nrecs,
                                                    TblDbl,"X_IMAGE",TblOff(PSTARIMAGEEXT,filestarimage.X_IMAGE),
                                                    TblDbl,"Y_IMAGE",TblOff(PSTARIMAGEEXT,filestarimage.Y_IMAGE),
                                                    TblDbl,"MAG_ISO",TblOff(PSTARIMAGEEXT,filestarimage.MAG_ISO),
                                                    TblDbl,"ra",TblOff(PSTARIMAGEEXT,filestarimage.ra),
                                                    TblDbl,"dec",TblOff(PSTARIMAGEEXT,filestarimage.dec),
                                                    TblDbl,"Date",TblOff(PSTARIMAGEEXT,filestarimage.Date),
                                                    TblDbl,"FLUX_ISO",TblOff(PSTARIMAGEEXT,filestarimage.FLUX_ISO),
                                                    TblDbl,"MAG_APER",TblOff(PSTARIMAGEEXT,filestarimage.MAG_APER),
                                                    TblDbl,"MAG_AUTO",TblOff(PSTARIMAGEEXT,filestarimage.MAG_AUTO),
                                                    TblDbl,"KRON_RADIUS",TblOff(PSTARIMAGEEXT,filestarimage.KRON_RADIUS),
                                                    TblDbl,"BACKGROUND",TblOff(PSTARIMAGEEXT,filestarimage.BACKGROUND),
                                                    TblDbl,"FLUX_MAX",TblOff(PSTARIMAGEEXT,filestarimage.FLUX_MAX),
                                                    TblDbl,"THETA_J2000",TblOff(PSTARIMAGEEXT,filestarimage.THETA_J2000),
                                                    TblDbl,"ELLIPTICITY",TblOff(PSTARIMAGEEXT,filestarimage.ELLIPTICITY),
                                                    TblDbl,"ISOAREA_WORLD",TblOff(PSTARIMAGEEXT,filestarimage.ISOAREA_WORLD),
                                                    TblDbl,"FWHM_IMAGE",TblOff(PSTARIMAGEEXT,filestarimage.FWHM_IMAGE),
                                                    TblDbl,"FWHM_WORLD",TblOff(PSTARIMAGEEXT,filestarimage.FWHM_WORLD),
                                                    TblDbl,"plate_dist",TblOff(PSTARIMAGEEXT,filestarimage.plate_dist),
                                                    TblDbl,"Blendedmag",TblOff(PSTARIMAGEEXT,filestarimage.Blendedmag),
                                                    TblDbl,"dradRMS2",TblOff(PSTARIMAGEEXT,filestarimage.dradRMS2),
                                                    TblDbl,"ra_2",TblOff(PSTARIMAGEEXT,filestarimage.ra_2),
                                                    TblDbl,"dec_2",TblOff(PSTARIMAGEEXT,filestarimage.dec_2),

                                                    TblFlt,"magcal_iso",TblOff(PSTARIMAGEEXT,filestarimage.magcal_iso),
                                                    TblFlt,"magcal_iso_rms",TblOff(PSTARIMAGEEXT,filestarimage.magcal_iso_rms),
                                                    TblFlt,"magcal_local",TblOff(PSTARIMAGEEXT,filestarimage.magcal_local),
                                                    TblFlt,"magcal_local_rms",TblOff(PSTARIMAGEEXT,filestarimage.magcal_local_rms),
                                                    TblFlt,"limiting_mag_local",TblOff(PSTARIMAGEEXT,filestarimage.limiting_mag_local),
                                                    TblFlt,"magcal_local_error",TblOff(PSTARIMAGEEXT,filestarimage.magcal_local_error),
                                                    TblFlt,"magcor_local",TblOff(PSTARIMAGEEXT,filestarimage.magcor_local),
                                                    TblFlt,"extinction",TblOff(PSTARIMAGEEXT,filestarimage.extinction),
                                                    TblFlt,"magcal_magdep",TblOff(PSTARIMAGEEXT,filestarimage.magcal_magdep),
                                                    TblFlt,"magcal_magdep_rms",TblOff(PSTARIMAGEEXT,filestarimage.magcal_magdep_rms),
                                                    TblFlt,"RaPM",TblOff(PSTARIMAGEEXT,filestarimage.RaPM),
                                                    TblFlt,"DecPM",TblOff(PSTARIMAGEEXT,filestarimage.DecPM),
                          TblInt,"A2FLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.A2FLAGS),
                          TblInt,"B2FLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.B2FLAGS),
                          TblFlt,"timeAccuracy",TblOff(PFILESTARIMAGEEXT,filestarimage.timeAccuracy),
                          TblInt,"maskIndex",TblOff(PFILESTARIMAGEEXT,filestarimage.maskIndex),

                                                    TblInt,"gsc_bin_index",TblOff(PSTARIMAGEEXT,filestarimage.gsc_bin_index),
                                                    TblInt,"plateNumber",TblOff(PSTARIMAGEEXT,filestarimage.plateNumber),
                                                    TblInt,"NUMBER",TblOff(PSTARIMAGEEXT,filestarimage.NUMBER),
                                                    TblInt,"versionId",TblOff(PSTARIMAGEEXT,filestarimage.versionId),
                                                    TblInt,"AFLAGS",TblOff(PSTARIMAGEEXT,filestarimage.AFLAGS),
                                                    TblInt,"BFLAGS",TblOff(PSTARIMAGEEXT,filestarimage.BFLAGS),
                                                    TblInt,"ISO0",TblOff(PSTARIMAGEEXT,filestarimage.ISO0),
                                                    TblInt,"ISO1",TblOff(PSTARIMAGEEXT,filestarimage.ISO1),
                                                    TblInt,"ISO2",TblOff(PSTARIMAGEEXT,filestarimage.ISO2),
                                                    TblInt,"ISO3",TblOff(PSTARIMAGEEXT,filestarimage.ISO3),
                                                    TblInt,"ISO4",TblOff(PSTARIMAGEEXT,filestarimage.ISO4),
                                                    TblInt,"ISO5",TblOff(PSTARIMAGEEXT,filestarimage.ISO5),
                                                    TblInt,"ISO6",TblOff(PSTARIMAGEEXT,filestarimage.ISO6),
                                                    TblInt,"ISO7",TblOff(PSTARIMAGEEXT,filestarimage.ISO7),
                                                    TblInt,"npoints_local",TblOff(PSTARIMAGEEXT,filestarimage.npoints_local),
                                                    TblInt,"rejectFlag",TblOff(PSTARIMAGEEXT,filestarimage.rejectFlag),
                                                    TblSht,"local_bin_index",TblOff(PSTARIMAGEEXT,filestarimage.local_bin_index),
                                                    TblByt,"seriesId",TblOff(PSTARIMAGEEXT,filestarimage.seriesId),
                                                    TblByt,"exposureNumber",TblOff(PSTARIMAGEEXT,filestarimage.exposureNumber),
                                                    TblByt,"solutionNumber",TblOff(PSTARIMAGEEXT,filestarimage.solutionNumber),
                                                    TblByt,"spatial_bin",TblOff(PSTARIMAGEEXT,filestarimage.spatial_bin),
                                                    TblInt,"mosaicNumber",TblOff(PSTARIMAGEEXT,mosaicNumber),
                                                    TblInt,"quality",TblOff(PSTARIMAGEEXT,quality),
                                                    TblInt,"plateVersionId",TblOff(PSTARIMAGEEXT,plateVersionId),
                                                    TblInt,"magdep_bin",TblOff(PSTARIMAGEEXT,filestarimage.magdep_bin),
                                                    TblBuf,"REFNumber",TblOff(PSTARIMAGEEXT,REFNumberC),2*MAX_REF,
                                                    0,"end",0);
    if (db_table == NULL) {
        printf("ERROR: Failed to read table for %s\n",db_name);
        return(-1);
    }
    numMagnitudes = db_nrecs;

    if (db_handle != NULL) {
        Close(db_handle);
    }
    if (db_header != NULL) {
        table_hdrfree(db_header);
    }


    /* Read in the object list */
    object_handle = Open(object_name,"r");
    if (object_handle == NULL) {
        errorFlag = 1;
        printf("ERROR: Failed to find the db file %s\n",object_name);
    } else {
        if (verbose) {
            printf("Found db file %s\n",object_name);
        }
    }


    object_table = table_loadva(object_handle,
                                                            &object_header,
                                                            NULL, /* hbase */
                                                            NULL, /* rows */
                                                            NULL,
                                                            sizeof(OBJECTTABLE),
                                                            &object_nrecs,
                                                            TblDbl,"ra",TblOff(POBJECTTABLE,ra),
                                                            TblDbl,"dec",TblOff(POBJECTTABLE,dec),
                                                            TblDbl,"cra",TblOff(POBJECTTABLE,cra),
                                                            TblDbl,"cdec",TblOff(POBJECTTABLE,cdec),
                                                            TblDbl,"drad",TblOff(POBJECTTABLE,drad),
                                                            TblDbl,"magcal_magdep",TblOff(POBJECTTABLE,magcal_magdep),
                                                            TblInt,"npoints",TblOff(POBJECTTABLE,npoints),
                                                            TblInt,"nplot",TblOff(POBJECTTABLE,nplot),
                                                            TblBuf,"REF",TblOff(POBJECTTABLE,REF),MAX_REF,
                                                            TblBuf,"src_name",TblOff(POBJECTTABLE,src_name),MAX_SRC_LENGTH,
                                                            TblBuf,"nearbyObjects",TblOff(POBJECTTABLE,nearbyObjects),MAX_NEARBY_OBJECTS_STRING,
                                                            0,"end",0);
    if (object_table == NULL) {
        printf("ERROR: Failed to read table for %s\n",object_name);
        return(-1);
    }
    for (object_index = 0; object_index < object_nrecs; object_index++) {
        pObject = &object_table[object_index];
        while ((charPtr = strchr(pObject->src_name,' ')) != NULL) {
            *charPtr = '_';
        }
        pObject->pointsfound = 0;
        pObject->written = 0;
        GetREFNumber(pObject->REF,&pObject->REFNumber,&pObject->refType,1,0);
    }
    if (object_nrecs > 1) {
        full_b_name = (char *)calloc(object_nrecs*MAX_FILENAME,sizeof(char));
        if (full_b_name == NULL) {
            printf("ERROR: failed to allocate full_b_name\n");
            exit(-1);
        }
        sprintf(full_b_name,"cd %s;tar -cf %s ",tmpdir,short_tar_txt_name);
        full_c_name = (char *)calloc(object_nrecs*MAX_FILENAME,sizeof(char));
        if (full_c_name == NULL) {
            printf("ERROR: failed to allocate full_c_name\n");
            exit(-1);
        }
        sprintf(full_c_name,"cd %s;tar -cf %s ",tmpdir,short_tar_db_name);
        full_d_name = (char *)calloc(object_nrecs*MAX_FILENAME,sizeof(char));
        if (full_d_name == NULL) {
            printf("ERROR: failed to allocate full_d_name\n");
            exit(-1);
        }
        sprintf(full_d_name,"cd %s;tar -cf %s ",tmpdir,short_tar_xml_name);

        full_f_name = (char *)calloc(object_nrecs*MAX_FILENAME,sizeof(char));
        if (full_f_name == NULL) {
            printf("ERROR: failed to allocate full_f_name\n");
            exit(-1);
        }
        sprintf(full_f_name,"cd %s;tar -cf %s ",tmpdir,tar_txt_name);

        full_g_name = (char *)calloc(object_nrecs*MAX_FILENAME,sizeof(char));
        if (full_g_name == NULL) {
            printf("ERROR: failed to allocate full_g_name\n");
            exit(-1);
        }
        sprintf(full_g_name,"cd %s;tar -cf %s ",tmpdir,tar_db_name);



        full_h_name = (char *)calloc(object_nrecs*MAX_FILENAME,sizeof(char));
        if (full_h_name == NULL) {
            printf("ERROR: failed to allocate full_h_name\n");
            exit(-1);
        }
        sprintf(full_h_name,"cd %s;tar -cf %s ",tmpdir,tar_xml_name);
    }

    if (object_handle != NULL) {
        Close(object_handle);
    }
    if (object_header != NULL) {
        table_hdrfree(object_header);
    }
    output_table = (PSTARIMAGEEXT)calloc(numMagnitudes,sizeof(STARIMAGEEXT));

    /* Now purge the database file of any objects that are not in the object list */
  for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
    pStarImageExt = &db_table[magnitudeIndex];
    pPhotPlates = &pStarImageExt->photplates;
    memset(pPhotPlates,0,sizeof(PHOTPLATES));
    pPhotPlates->mosaicNumber = 99;
    pPhotPlates->quality = QUALITY_UNINITIALIZED;

        /* Covert the REFNumber to a REF */
        nvals = sscanf(pStarImageExt->REFNumberC,"%lld",&pStarImageExt->filestarimage.REFNumber);
        if (nvals != 1) {
            printf("ERROR: failed to decode REFNumber %s\n",pStarImageExt->REFNumberC);
            exit(-1);
        }

        pOutputImageExt = &output_table[output_nrecs];
        for (object_index = 0; object_index < object_nrecs; object_index++) {
            pObject = &object_table[object_index];
            if (pObject->REFNumber == pStarImageExt->filestarimage.REFNumber) {
                memcpy(pOutputImageExt,pStarImageExt,sizeof(STARIMAGEEXT));
                pObject->pointsfound++;
                pOutputImageExt->object_index = object_index;
                output_nrecs++;
            }
        }
    }
    if (output_nrecs == 0) {
        printf("ERROR: no output records\n");
        exit(-1);
    }
    /* Sort everything by REF number first and date second */
  qsort((void*)output_table,output_nrecs,sizeof(STARIMAGEEXT),OutputTableCompare);
  /* At this point, we need to get the plate quality bits */
    for (output_index = 0; output_index < output_nrecs; output_index++) {
        pOutputImageExt = &output_table[output_index];
    pPhotPlates = &pOutputImageExt->photplates;
    /* First search for a previous invocation of GetPhotPlate */
    for (output_index2 = 0; output_index2 < output_index; output_index2++) {
      pPrevOutputImageExt = &output_table[output_index2];
      pPrevPhotPlates = &pPrevOutputImageExt->photplates;
      if ((pPrevOutputImageExt->filestarimage.plateNumber == pOutputImageExt->filestarimage.plateNumber) &&
          (pPrevOutputImageExt->filestarimage.seriesId == pOutputImageExt->filestarimage.seriesId)) {
        /* We already got this information, just copy it over */
        memcpy(pPhotPlates,pPrevPhotPlates,sizeof(PHOTPLATES));
        pOutputImageExt->mosaicNumber = pPhotPlates->mosaicNumber;
        pOutputImageExt->quality = pPhotPlates->quality;
        pOutputImageExt->plateVersionId = pPhotPlates->versionId;
        break;
      }
    }
    if (output_index2 >= output_index) {
      GetPhotPlateExt(pPhotConnection,pOutputImageExt,catalogString);
    }
  }


    outputhandle = fopen(output_db,"wt");
    if (outputhandle == NULL) {
        printf("ERROR: failed to open output file %s\n",output_db);
    }
    outputtxthandle = fopen(output_txt_db,"wt");
    if (outputtxthandle == NULL) {
        printf("ERROR: failed to open output file %s\n",output_txt_db);
    }
    shortoutputhandle = fopen(short_db_name,"wt");
    if (shortoutputhandle == NULL) {
        printf("ERROR: failed to open output file %s\n",short_db_name);
    }
    shorttxtoutputhandle = fopen(short_txt_name,"wt");
    if (shorttxtoutputhandle == NULL) {
        printf("ERROR: failed to open output file %s\n",short_txt_name);
    }
  sprintf(starbase_title," Photometry data using the %s calibration catalog",source);
  if (enableRematch != 0) {
    strcat(starbase_title,", optimized for transients.");
  }
#if 0
  printf("line %d web_data enableRematch %d starbase_title %s\n",__LINE__,enableRematch,starbase_title);
#endif

  /* Hack here!  Adding an ASCII 0x7 = "BEL" = "\a" turns the file into a binary file in the eyes of firefox and chrome */
  strcpy(starbase_txt_title,starbase_title);
  strcat(starbase_txt_title,"\a");
#if (RELEASE_EXPERIMENTAL == 1)
  fprintf(shortoutputhandle,"title %s\n",starbase_title);
    fprintf(shortoutputhandle,"REF\tDate\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\tcatalogNumber\n");
    fprintf(shortoutputhandle,"---\t----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\t-------------\n");

  fprintf(shorttxtoutputhandle,"title %s\n",starbase_txt_title);
    fprintf(shorttxtoutputhandle,"REF\tDate\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\tcatalogNumber\n");
    fprintf(shorttxtoutputhandle,"---\t----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\t-------------\n");

#else /* (RELEASE_EXPERIMENTAL == 1) */


  fprintf(shortoutputhandle,"title %s\n",starbase_title);
    fprintf(shortoutputhandle,"REF\tDate\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\n");
    fprintf(shortoutputhandle,"---\t----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\n");

  fprintf(shorttxtoutputhandle,"title %s\n",starbase_txt_title);
    fprintf(shorttxtoutputhandle,"REF\tDate\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\n");
    fprintf(shorttxtoutputhandle,"---\t----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\n");
#endif /* (RELEASE_EXPERIMENTAL == 1) */

    writeOutputHeader = 1;
    writeTxtOutputHeader = 1;


    for (output_index = 0; output_index < output_nrecs; output_index++) {
        double short_year;
        char REF[MAX_REF];
        char short_plate[MAX_FILENAME];
        pOutputImageExt = &output_table[output_index];
        pPhotPlates->versionId = pOutputImageExt->plateVersionId;
        WriteStarbaseRecord(&pOutputImageExt->filestarimage,pPhotPlates,outputhandle,starbase_title,&writeOutputHeader,catalogNumber);
        WriteStarbaseRecord(&pOutputImageExt->filestarimage,pPhotPlates,outputtxthandle,starbase_txt_title,&writeTxtOutputHeader,catalogNumber);
        GetREF(pOutputImageExt->filestarimage.REFNumber,REF,1,1);
        if (pOutputImageExt->filestarimage.solutionNumber == 0) {
            sprintf(short_plate,"%s%05d",GetSeriesString(pOutputImageExt->filestarimage.seriesId,0),pOutputImageExt->filestarimage.plateNumber);
        } else {
            sprintf(short_plate,"%s%05d_s%d",GetSeriesString(pOutputImageExt->filestarimage.seriesId,0),pOutputImageExt->filestarimage.plateNumber,pOutputImageExt->filestarimage.solutionNumber);
        }
        short_year = jd2ep(pOutputImageExt->filestarimage.Date);

    FormatFlagsBits(pOutputImageExt->filestarimage.AFLAGS,AFLAGSBuffer,AFLAGSbitsBuffer,sizeof(AFLAGSBuffer),0);
    FormatFlagsBits(pOutputImageExt->quality,qualityBuffer,qualitybitsBuffer,sizeof(qualityBuffer),1);
    FormatFlagsBits(pOutputImageExt->filestarimage.BFLAGS,BFLAGSBuffer,BFLAGSbitsBuffer,sizeof(BFLAGSBuffer),0);

    tmpCatalogNumber = catalogNumber;
    if ((RELEASE_EXPERIMENTAL != 0) && (catalogNumber == CATALOG_EXPERIMENTAL)) {
      tmpCatalogNumber = pOutputImageExt->filestarimage.catalogNumber;
    } else {
      tmpCatalogNumber = catalogNumber;
    }

#if (RELEASE_EXPERIMENTAL == 1)
        fprintf(shortoutputhandle,"%s\t%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\t%d\n",
                        REF,
                        pOutputImageExt->filestarimage.Date,
                        short_year,
                        pOutputImageExt->filestarimage.magcal_magdep,
                        pOutputImageExt->filestarimage.magcal_local_rms,
                        pOutputImageExt->filestarimage.limiting_mag_local,
                        pOutputImageExt->filestarimage.ra,
                        pOutputImageExt->filestarimage.dec,
                        pOutputImageExt->filestarimage.THETA_J2000,
                        pOutputImageExt->filestarimage.ELLIPTICITY,
                        short_plate,
                        pOutputImageExt->filestarimage.versionId,
                        AFLAGSBuffer,
            BFLAGSBuffer,
                        pOutputImageExt->filestarimage.NUMBER,
            qualityBuffer,
                        pOutputImageExt->filestarimage.timeAccuracy,
            AFLAGSbitsBuffer,
            BFLAGSbitsBuffer,
            qualitybitsBuffer,
            tmpCatalogNumber);
        fprintf(shorttxtoutputhandle,"%s\t%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\t%d\n",
                        REF,
                        pOutputImageExt->filestarimage.Date,
                        short_year,
                        pOutputImageExt->filestarimage.magcal_magdep,
                        pOutputImageExt->filestarimage.magcal_local_rms,
                        pOutputImageExt->filestarimage.limiting_mag_local,
                        pOutputImageExt->filestarimage.ra,
                        pOutputImageExt->filestarimage.dec,
                        pOutputImageExt->filestarimage.THETA_J2000,
                        pOutputImageExt->filestarimage.ELLIPTICITY,
                        short_plate,
                        pOutputImageExt->filestarimage.versionId,
                        AFLAGSBuffer,
            BFLAGSBuffer,
                        pOutputImageExt->filestarimage.NUMBER,
            qualityBuffer,
                        pOutputImageExt->filestarimage.timeAccuracy,
            AFLAGSbitsBuffer,
            BFLAGSbitsBuffer,
            qualitybitsBuffer,
            tmpCatalogNumber);

#else /* (RELEASE_EXPERIMENTAL == 1) */

        fprintf(shortoutputhandle,"%s\t%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\n",
                        REF,
                        pOutputImageExt->filestarimage.Date,
                        short_year,
                        pOutputImageExt->filestarimage.magcal_magdep,
                        pOutputImageExt->filestarimage.magcal_local_rms,
                        pOutputImageExt->filestarimage.limiting_mag_local,
                        pOutputImageExt->filestarimage.ra,
                        pOutputImageExt->filestarimage.dec,
                        pOutputImageExt->filestarimage.THETA_J2000,
                        pOutputImageExt->filestarimage.ELLIPTICITY,
                        short_plate,
                        pOutputImageExt->filestarimage.versionId,
                        AFLAGSBuffer,
            BFLAGSBuffer,
                        pOutputImageExt->filestarimage.NUMBER,
            qualityBuffer,
                        pOutputImageExt->filestarimage.timeAccuracy,
            AFLAGSbitsBuffer,
            BFLAGSbitsBuffer,
            qualitybitsBuffer);
        fprintf(shorttxtoutputhandle,"%s\t%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\n",
                        REF,
                        pOutputImageExt->filestarimage.Date,
                        short_year,
                        pOutputImageExt->filestarimage.magcal_magdep,
                        pOutputImageExt->filestarimage.magcal_local_rms,
                        pOutputImageExt->filestarimage.limiting_mag_local,
                        pOutputImageExt->filestarimage.ra,
                        pOutputImageExt->filestarimage.dec,
                        pOutputImageExt->filestarimage.THETA_J2000,
                        pOutputImageExt->filestarimage.ELLIPTICITY,
                        short_plate,
                        pOutputImageExt->filestarimage.versionId,
                        AFLAGSBuffer,
            BFLAGSBuffer,
                        pOutputImageExt->filestarimage.NUMBER,
            qualityBuffer,
                        pOutputImageExt->filestarimage.timeAccuracy,
            AFLAGSbitsBuffer,
            BFLAGSbitsBuffer,
            qualitybitsBuffer);
#endif /* (RELEASE_EXPERIMENTAL == 1) */
        if ((objectoutputhandle == NULL) && (shortobjectoutputhandle == NULL)) {
            object_index = pOutputImageExt->object_index;
            pObject = &object_table[object_index];

            if (pObject->written != 0) {
                printf("ERROR: double write for object %zu %s %s\n",object_index,pObject->src_name,pObject->REF);
                exit(-1);
            }
            pObject->written = 1;
            sprintf(objectrootname,"%s_%s_%04.0f",pObject->src_name,pObject->REF,pObject->drad);
            strcpy(objectfilename,tmpdir);
            strcat(objectfilename,"/");
            strcat(objectfilename,objectrootname);
      strcpy(objecttxtfilename,objectfilename);
            strcat(objectfilename,".db");
      strcat(objecttxtfilename,".txt");
            strcpy(shortobjectfilename,tmpdir);
            strcat(shortobjectfilename,"/short_");
            strcat(shortobjectfilename,objectrootname);
      strcpy(shorttxtobjectfilename,shortobjectfilename);
            strcat(shortobjectfilename,".db");
      strcat(shorttxtobjectfilename,".txt");
            objectoutputhandle = fopen(objectfilename,"wt");
            if (objectoutputhandle == NULL) {
                printf("ERROR: failed to open %s\n",objectfilename);
                exit(-1);
            }
            objecttxtoutputhandle = fopen(objecttxtfilename,"wt");
            if (objecttxtoutputhandle == NULL) {
                printf("ERROR: failed to open %s\n",objecttxtfilename);
                exit(-1);
            }
            shortobjectoutputhandle = fopen(shortobjectfilename,"wt");
            if (shortobjectoutputhandle == NULL) {
                printf("ERROR: failed to open %s\n",shortobjectfilename);
                exit(-1);
            }
            shorttxtobjectoutputhandle = fopen(shorttxtobjectfilename,"wt");
            if (shorttxtobjectoutputhandle == NULL) {
                printf("ERROR: failed to open %s\n",shorttxtobjectfilename);
                exit(-1);
            }
      writeObjectHeader = 1;
      writeTxtObjectHeader = 1;
      sprintf(starbase_title," Photometry data for %s from search centered on %s using the %s calibration catalog",pObject->REF,pObject->src_name,catalogText[catalogNumber]);
      if (enableRematch != 0) {
        strcat(starbase_title,", optimized for transients.");
      }
#if 0
      printf("line %d web_data enableRematch %d starbase_title %s\n",__LINE__,enableRematch,starbase_title);
#endif
      /* Hack here!  Adding an ASCII 0x7 = "BEL" = "\a" turns the file into a binary file in the eyes of firefox and chrome */


      strcpy(starbase_txt_title,starbase_title);
      strcat(starbase_txt_title,"\a");

#if (RELEASE_EXPERIMENTAL == 1)
      fprintf(shortobjectoutputhandle,"title %s\n",starbase_title);
            fprintf(shortobjectoutputhandle,"REF\tDate\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\tcatalogNumber\n");
            fprintf(shortobjectoutputhandle,"---\t----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\t-------------\n");

      fprintf(shorttxtobjectoutputhandle,"title %s\n",starbase_txt_title);
            fprintf(shorttxtobjectoutputhandle,"REF\tDate\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\tcatalogNumber\n");
            fprintf(shorttxtobjectoutputhandle,"---\t----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\t-------------\n");
#else /* (RELEASE_EXPERIMENTAL == 1) */
      fprintf(shortobjectoutputhandle,"title %s\n",starbase_title);
            fprintf(shortobjectoutputhandle,"REF\tDate\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\n");
            fprintf(shortobjectoutputhandle,"---\t----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\n");

      fprintf(shorttxtobjectoutputhandle,"title %s\n",starbase_txt_title);
            fprintf(shorttxtobjectoutputhandle,"REF\tDate\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\n");
            fprintf(shorttxtobjectoutputhandle,"---\t----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\n");
#endif /* (RELEASE_EXPERIMENTAL == 1) */

        }
        WriteStarbaseRecord(&pOutputImageExt->filestarimage,pPhotPlates,objectoutputhandle,starbase_title,&writeObjectHeader,catalogNumber);
        WriteStarbaseRecord(&pOutputImageExt->filestarimage,pPhotPlates,objecttxtoutputhandle,starbase_txt_title,&writeTxtObjectHeader,catalogNumber);
    FormatFlagsBits(pOutputImageExt->filestarimage.AFLAGS,AFLAGSBuffer,AFLAGSbitsBuffer,sizeof(AFLAGSBuffer),0);
    FormatFlagsBits(pOutputImageExt->quality,qualityBuffer,qualitybitsBuffer,sizeof(qualityBuffer),1);
    FormatFlagsBits(pOutputImageExt->filestarimage.BFLAGS,BFLAGSBuffer,BFLAGSbitsBuffer,sizeof(BFLAGSBuffer),0);

    tmpCatalogNumber = catalogNumber;
    if ((RELEASE_EXPERIMENTAL != 0) && (catalogNumber == CATALOG_EXPERIMENTAL)) {
      tmpCatalogNumber = pOutputImageExt->filestarimage.catalogNumber;
    } else {
      tmpCatalogNumber = catalogNumber;
    }
#if (RELEASE_EXPERIMENTAL == 1)

        fprintf(shortobjectoutputhandle,"%s\t%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\t%d\n",
                        REF,
                        pOutputImageExt->filestarimage.Date,
                        short_year,
                        pOutputImageExt->filestarimage.magcal_magdep,
                        pOutputImageExt->filestarimage.magcal_local_rms,
                        pOutputImageExt->filestarimage.limiting_mag_local,
                        pOutputImageExt->filestarimage.ra,
                        pOutputImageExt->filestarimage.dec,
                        pOutputImageExt->filestarimage.THETA_J2000,
                        pOutputImageExt->filestarimage.ELLIPTICITY,
                        short_plate,
                        pOutputImageExt->filestarimage.versionId,
                        AFLAGSBuffer,
            BFLAGSBuffer,
                        pOutputImageExt->filestarimage.NUMBER,
            qualityBuffer,
                        pOutputImageExt->filestarimage.timeAccuracy,
            AFLAGSbitsBuffer,
            BFLAGSbitsBuffer,
            qualitybitsBuffer,
            tmpCatalogNumber);
        fprintf(shorttxtobjectoutputhandle,"%s\t%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\t%d\n",
                        REF,
                        pOutputImageExt->filestarimage.Date,
                        short_year,
                        pOutputImageExt->filestarimage.magcal_magdep,
                        pOutputImageExt->filestarimage.magcal_local_rms,
                        pOutputImageExt->filestarimage.limiting_mag_local,
                        pOutputImageExt->filestarimage.ra,
                        pOutputImageExt->filestarimage.dec,
                        pOutputImageExt->filestarimage.THETA_J2000,
                        pOutputImageExt->filestarimage.ELLIPTICITY,
                        short_plate,
                        pOutputImageExt->filestarimage.versionId,
                        AFLAGSBuffer,
            BFLAGSBuffer,
                        pOutputImageExt->filestarimage.NUMBER,
            qualityBuffer,
                        pOutputImageExt->filestarimage.timeAccuracy,
            AFLAGSbitsBuffer,
            BFLAGSbitsBuffer,
            qualitybitsBuffer,
            tmpCatalogNumber);
#else /* (RELEASE_EXPERIMENTAL == 1) */

        fprintf(shortobjectoutputhandle,"%s\t%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\n",
                        REF,
                        pOutputImageExt->filestarimage.Date,
                        short_year,
                        pOutputImageExt->filestarimage.magcal_magdep,
                        pOutputImageExt->filestarimage.magcal_local_rms,
                        pOutputImageExt->filestarimage.limiting_mag_local,
                        pOutputImageExt->filestarimage.ra,
                        pOutputImageExt->filestarimage.dec,
                        pOutputImageExt->filestarimage.THETA_J2000,
                        pOutputImageExt->filestarimage.ELLIPTICITY,
                        short_plate,
                        pOutputImageExt->filestarimage.versionId,
                        AFLAGSBuffer,
            BFLAGSBuffer,
                        pOutputImageExt->filestarimage.NUMBER,
            qualityBuffer,
                        pOutputImageExt->filestarimage.timeAccuracy,
            AFLAGSbitsBuffer,
            BFLAGSbitsBuffer,
            qualitybitsBuffer);
        fprintf(shorttxtobjectoutputhandle,"%s\t%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\n",
                        REF,
                        pOutputImageExt->filestarimage.Date,
                        short_year,
                        pOutputImageExt->filestarimage.magcal_magdep,
                        pOutputImageExt->filestarimage.magcal_local_rms,
                        pOutputImageExt->filestarimage.limiting_mag_local,
                        pOutputImageExt->filestarimage.ra,
                        pOutputImageExt->filestarimage.dec,
                        pOutputImageExt->filestarimage.THETA_J2000,
                        pOutputImageExt->filestarimage.ELLIPTICITY,
                        short_plate,
                        pOutputImageExt->filestarimage.versionId,
                        AFLAGSBuffer,
            BFLAGSBuffer,
                        pOutputImageExt->filestarimage.NUMBER,
            qualityBuffer,
                        pOutputImageExt->filestarimage.timeAccuracy,
            AFLAGSbitsBuffer,
            BFLAGSbitsBuffer,
            qualitybitsBuffer);
#endif /* (RELEASE_EXPERIMENTAL == 1) */

        pNextOutputImageExt = &output_table[output_index+1];
        if ((output_index+1 == output_nrecs) ||
                (pOutputImageExt->object_index != pNextOutputImageExt->object_index)) {
            if (objectoutputhandle != NULL) {
                fclose(objectoutputhandle);
                objectoutputhandle = NULL;
            }
            if (objecttxtoutputhandle != NULL) {
                fclose(objecttxtoutputhandle);
                objectoutputhandle = NULL;
            }
            if (shortobjectoutputhandle  != NULL) {
                fclose(shortobjectoutputhandle);
                shortobjectoutputhandle = NULL;
            }
            if (shorttxtobjectoutputhandle  != NULL) {
                fclose(shorttxtobjectoutputhandle);
                shorttxtobjectoutputhandle = NULL;
            }
            writeObjectHeader = 1;
            BUFPRINTF(objectdescription,"%s with %d points at %.0f arcsec from search center",objectrootname,pObject->npoints,pObject->drad);
            WriteOutputLinks(objectdescription,lsbin,binaries,shortobjectfilename,objectfilename,b_name,c_name,d_name,f_name,g_name,h_name);
            if (object_nrecs > 1) {
                strcat(full_b_name," ");
                strcat(full_b_name,b_name);

                strcat(full_c_name," ");
                strcat(full_c_name,c_name);

                strcat(full_d_name," ");
                strcat(full_d_name,d_name);

                strcat(full_f_name," ");
                strcat(full_f_name,f_name);

                strcat(full_g_name," ");
                strcat(full_g_name,g_name);

                strcat(full_h_name," ");
                strcat(full_h_name,h_name);


            }


        }

    }
    fclose(outputhandle);
    outputhandle = NULL;
    fclose(shortoutputhandle);
    shortoutputhandle = NULL;
    fclose(shorttxtoutputhandle);
    shorttxtoutputhandle = NULL;
    if (object_nrecs > 1) {
        WriteOutputLinks("Table of all the data:",lsbin,binaries,short_db_name,output_db,b_name,c_name,d_name,f_name,g_name,h_name);
        printf("<br /><b>Collections of the above Files</b> ");

        strcat(full_b_name,"\n");
        system(full_b_name);
        slashPtr = strrchr(short_tar_txt_name,'/');
        if (slashPtr == NULL) {
            slashPtr = short_tar_txt_name;
        } else {
            slashPtr++;
        }
        charPtr = strstr(short_tar_txt_name,"/tmp");
        if (charPtr == NULL) {
            charPtr = short_tar_txt_name;
        } else {
            charPtr++;
        }
        printf("<br />B: <a href=\"%s\">%s</a> ",charPtr,slashPtr);

        strcat(full_c_name,"\n");
        system(full_c_name);
        slashPtr = strrchr(short_tar_db_name,'/');
        if (slashPtr == NULL) {
            slashPtr = short_tar_db_name;
        } else {
            slashPtr++;
        }
        charPtr = strstr(short_tar_db_name,"/tmp");
        if (charPtr == NULL) {
            charPtr = short_tar_db_name;
        } else {
            charPtr++;
        }
        printf("<br />C: <a href=\"%s\">%s</a> ",charPtr,slashPtr);

        strcat(full_d_name,"\n");
        system(full_d_name);
        slashPtr = strrchr(short_tar_xml_name,'/');
        if (slashPtr == NULL) {
            slashPtr = short_tar_xml_name;
        } else {
            slashPtr++;
        }
        charPtr = strstr(short_tar_xml_name,"/tmp");
        if (charPtr == NULL) {
            charPtr = short_tar_xml_name;
        } else {
            charPtr++;
        }
        printf("<br />D: <a href=\"%s\">%s</a> ",charPtr,slashPtr);

        strcat(full_f_name,"\n");
        system(full_f_name);
        slashPtr = strrchr(tar_txt_name,'/');
        if (slashPtr == NULL) {
            slashPtr = tar_txt_name;
        } else {
            slashPtr++;
        }
        charPtr = strstr(tar_txt_name,"/tmp");
        if (charPtr == NULL) {
            charPtr = tar_txt_name;
        } else {
            charPtr++;
        }
        printf("<br />F: <a href=\"%s\">%s</a> ",charPtr,slashPtr);

        strcat(full_g_name,"\n");
        system(full_g_name);
        slashPtr = strrchr(tar_db_name,'/');
        if (slashPtr == NULL) {
            slashPtr = tar_db_name;
        } else {
            slashPtr++;
        }
        charPtr = strstr(tar_db_name,"/tmp");
        if (charPtr == NULL) {
            charPtr = tar_db_name;
        } else {
            charPtr++;
        }
        printf("<br />G: <a href=\"%s\">%s</a> ",charPtr,slashPtr);

        strcat(full_h_name,"\n");
        system(full_h_name);
        slashPtr = strrchr(tar_xml_name,'/');
        if (slashPtr == NULL) {
            slashPtr = tar_xml_name;
        } else {
            slashPtr++;
        }
        charPtr = strstr(tar_xml_name,"/tmp");
        if (charPtr == NULL) {
            charPtr = tar_xml_name;
        } else {
            charPtr++;
        }
        printf("<br />H: <a href=\"%s\">%s</a> ",charPtr,slashPtr);


    }

  if (db_table != NULL) {
    if (dbCallocFlag) {
      free(db_table);
      db_table = NULL;
    } else {
      Free(db_table);
      db_table = NULL;
    }
  }
    if (output_table != NULL) {
        free(output_table);
    }
  if (db_indices != NULL) {
    free(db_indices);
  }
  if (short_indices != NULL) {
    free(short_indices);
  }

  if (listString != NULL) {
    free(listString);
  }
  if (pSortTable != NULL) {
    free(pSortTable);
  }

  return(0);
}

