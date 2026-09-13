// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* find_distortion.c
 *
 * This program has several modes of operation:
 * MODE_ORIGINAL: Given a list of files, find the difference between the TNX and TAN fit so we can estimate the maximum distortion.
 * MODE_CREATE_DISTORTION: For all spatial bins with valid lowess fit solutions, output a grid database of the distortion per bin.
 * MODE_READ_DISTORTION:   Read the grid database and generate per-series mean, rms, and count database files.
 * MODE_FIND_DISTORTION:   Combination of the above two steps with no distortion file printed.
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio -lmysqlclient  find_distortion.c pipelineutils.a -ltable -lutil  -lwcs -o find_distortion
 *
 * Sep 16, 2008 Edward J. Los - Initial version
 * Oct 15, 2008 Edward J. Los - Add MODE_CREATE_DISTORTION
 * Oct 24, 2008 Edward J. Los - Add MODE_READ_DISTORTION, MODE_FIND_DISTORTION
 *              TODO: Need mc and ac sub-series
 *              TODO: Need full plate rotation and transformation support
 * Dec 15, 2008 Edward J. Los - Read in the dmagcor corrections and eliminate any
 *              bins with errors
 * Dec 16, 2008 Edward J. Los - Zero starbase structures after deleting them.
 *                              Delete spatial bins starbase structures
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 *
 * 
 * find_distortion -v -m 1 -l /dasch/Pipeline/good_scans.debug -o /dasch/Pipeline/los1.txt -p /dasch/Pipeline/los.txt -s /dasch/Pipeline/los2.txt
 *
 * find_distortion -v -m 1 -l /dasch/Pipeline/total.list -o /dasch/backup/2008_10_16/distortion1.txt -p /dasch/backup/2008_10_16/distortion.txt -s /dasch/backup/2008_10_16/distortion2.txt
 *
 * find_distortion -v -m 3 -g /dasch/Pipeline/match/distortion.db -s /dasch/Pipeline/match/distortion_summary.txt
 *
 * find_distortion -v -m 4 -l /dasch/Pipeline/total.list -g /dasch/Pipeline/match/total.db -s /dasch/Pipeline/match/total_summary.txt
 *
 */


#include <math.h>
#include <sys/stat.h>
#include "table.h"
#include "mysql.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "time.h"
#include "pipelineutils.h"
#define MAX_FILENAME 256
#define MAX_BUFFER 256
#define MAX_LIST_STRING 25
#define SAMPLE_POINTS 49
#define ALLOC_INCREMENT 100
#if 1
#define MINIMUM_COUNT 20
#else
#define MINIMUM_COUNT 2
#endif

/* Values for programMode */
#define MODE_ORIGINAL          1  /* Original mode: output file is max distortion, plot file is distortion vs plate distance */
#define MODE_CREATE_DISTORTION 2  /* Read a plate list and create distortion.db as the output file */
#define MODE_READ_DISTORTION   3  /* Read distortion.db and create distortion_<series>_median.db  (median distortion map) 
                                                                   distortion_<series>_rms.db     (rms distortion map)
                                                                   distortion_<series>_count.db   (count of plates map) */
#define MODE_FIND_DISTORTION   4 /* Both MODE_CREATE_DISTORTION and MODE_FIND_DISTORTION without writing distortion.db */
#define MODE_MAX 4



char *
GetFITShead (char *filename,	/* FITS or IRAF file filename */
             int verbose);	/* Print error messages if nonzero */
struct WorldCoor *
GetFITSWCS (char *filename, char *header, int verbose,double * cra, double *cdec,double* dra, double *ddec, double* secpix,int* wp,int* hp,
	    int *sysout, double* eqout);

void
wcs2pix (
struct WorldCoor *wcs,	/* World coordinate system structure */
double	xpos,double ypos,	/* World coordinates in degrees */
double	*xpix,double *ypix,	/* Image coordinates in pixels */
int	*offscl);	/* 0 if within bounds, else off scale */


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

double wcsdist(
               double x1,
               double y1,	/* (RA,Dec) or (Long,Lat) in degrees */
               double x2,
               double y2);	/* (RA,Dec) or (Long,Lat) in degrees */


typedef struct _spatial_bin {
  int spatial_bin;
  double limiting_mag;
  double max_bright_mag;
  double limiting_iso;
  double max_bright_iso;
} SPATIAL_BIN,*PSPATIAL_BIN;


typedef struct _plateEntry {
  char series[MAX_SERIES_STRING];
  char listName[MAX_LIST_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int fileExists;
} PLATEENTRY,*PPLATEENTRY;

typedef struct _gridentry {
  char series[MAX_SERIES_STRING+1];
  int plateNumber;
  int spatial_bin;
  int dnx;
  int dny;
  int dix;
  int diy;
  int width;
  int height;
  double xval;
  double yval;
  double dxval;
  double dyval;
  double dPlateNumber;
} GRIDENTRY,*PGRIDENTRY;

typedef struct _seriestable {
  char gridname[MAX_FILENAME];
  int seriesCount;
  int plateCount;
  int seriesPlateCount;
  int allocCount[(SAMPLE_POINTS+1)*(SAMPLE_POINTS+1)];
  int curCount[(SAMPLE_POINTS+1)*(SAMPLE_POINTS+1)];
  double *dxvector[(SAMPLE_POINTS+1)*(SAMPLE_POINTS+1)];
  double *dyvector[(SAMPLE_POINTS+1)*(SAMPLE_POINTS+1)];
  char seriesStrings[MAX_SERIES_STRING*(MAX_SERIES+1)];
  int seriesPlates[MAX_SERIES+1];
  double maxmedian[MAX_SERIES+1];
  double maxrms[MAX_SERIES+1];
} SERIESTABLE, *PSERIESTABLE;

typedef struct _dmagcor {
  int nx;             /* Total bins in width */
  int ny;             /* Total bins in height */
  int ix;             /* width bin */
  int iy;             /* height bin */
  double zout;        /* Magnitude correction */
  double errout;      /* RMS of magnitude correction */
  int npout;          /* Number of points used for correction */
  int rejectFlag;     /* Reason, if any, for rejecting this bin */
  double dradAverage; /* drad average */
} DMAGCOR,*PDMAGCOR;


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

/* wf_gsclose -- procedure to free the surface descriptor */

static void
wf_gsclose (struct IRAFsurface *sf /* the surface descriptor */)

{
    if (sf != NULL) {
	if (sf->xbasis != NULL)
	    free (sf->xbasis);
	if (sf->ybasis != NULL)
	    free (sf->ybasis);
	if (sf->coeff != NULL)
	    free (sf->coeff);
	free (sf);
	}
    return;
}

int ReadGridEntry(PGRIDENTRY pGridEntry,
                   File grid_handle,
                   TableHead grid_header,
                   TblDescriptor grid_descriptor,
                   TableRow* grid_row,
                   char *rootname)
{
  *grid_row = table_rowget(grid_handle,grid_header,*grid_row,NULL,NULL,0);
  if (*grid_row == NULL) {
    return(0);
  }
  if (!table_loadrow(grid_handle,grid_header,*grid_row,grid_descriptor,(char *)pGridEntry)) {
    printf("ERROR: ReadGridEntry table_loadrow failed for %s\n",rootname);
    return(0);
  }
  /* This necessary because leading zeroes force conversion to octal! */
  pGridEntry->plateNumber = pGridEntry->dPlateNumber+0.0001;
  return(1);
}
void ProcessSeries(PSERIESTABLE pSeriesTable, PGRIDENTRY pOldGridEntry,int verbose)
{
  int dix;
  int diy;
  int irec;
  int result;
  char *dotPtr;
  char filename[2*MAX_FILENAME];
  FILE *medianHandle = NULL;
  FILE *countHandle = NULL;
  FILE *rmsHandle = NULL;
  double xmedian;
  double ymedian;
  double xrms;
  double yrms;
  /* Make up our three file names and open them */
  strcpy(filename,pSeriesTable->gridname);
  dotPtr = strstr(filename,".");
  if (dotPtr) {
    *dotPtr = 0;
  }
  strcat(filename,"_");
  strcat(filename,pOldGridEntry->series);
  strcat(filename,"_median.db");
  medianHandle = fopen(filename,"wt");
  if (medianHandle == NULL) {
    fprintf(stderr,"ERROR: failed to open file %s\n",filename);
    exit(-1);
  }
  strcpy(filename,pSeriesTable->gridname);
  dotPtr = strstr(filename,".");
  if (dotPtr) {
    *dotPtr = 0;
  }
  strcat(filename,"_");
  strcat(filename,pOldGridEntry->series);
  strcat(filename,"_count.db");
  countHandle = fopen(filename,"wt");
  if (countHandle == NULL) {
    fprintf(stderr,"ERROR: failed to open file %s\n",filename);
    exit(-1);
  }
  strcpy(filename,pSeriesTable->gridname);
  dotPtr = strstr(filename,".");
  if (dotPtr) {
    *dotPtr = 0;
  }
  strcat(filename,"_");
  strcat(filename,pOldGridEntry->series);
  strcat(filename,"_rms.db");
  rmsHandle = fopen(filename,"wt");
  if (rmsHandle == NULL) {
    fprintf(stderr,"ERROR: failed to open file %s\n",filename);
    exit(-1);
  }
  fprintf(medianHandle,"series\tdnx\tdny\tdix\tdiy\tdxmedian\tdymedian\n");
  fprintf(medianHandle,"------\t---\t---\t---\t---\t--------\t--------\n");
  fprintf(rmsHandle,"series\tdnx\tdny\tdix\tdiy\tdxrms\tdyrms\n");
  fprintf(rmsHandle,"------\t---\t---\t---\t---\t-----\t-----\n");
  fprintf(countHandle,"series\tdnx\tdny\tdix\tdiy\tcount\n");
  fprintf(countHandle,"------\t---\t---\t---\t---\t-----\n");


  for (dix = 0; dix <= SAMPLE_POINTS; dix++) {
    for (diy = 0; diy <= SAMPLE_POINTS; diy++) {
      irec = dix + (diy * (SAMPLE_POINTS+1));
      fprintf(countHandle,"%s\t%d\t%d\t%d\t%d\t%d\n",
              pOldGridEntry->series,
              SAMPLE_POINTS+1,
              SAMPLE_POINTS+1,
              dix,
              diy,
              pSeriesTable->curCount[irec]);
      if (pSeriesTable->curCount[irec] >= MINIMUM_COUNT) {
        double tmpmedian;
        result = CalcMedianAndRMS(pSeriesTable->curCount[irec],0,pSeriesTable->dxvector[irec],&xmedian,&xrms,1,3.0,0);
        result = CalcMedianAndRMS(pSeriesTable->curCount[irec],0,pSeriesTable->dyvector[irec],&ymedian,&yrms,1,3.0,0);
        tmpmedian = xmedian;
        if (tmpmedian < 0.0) {
          tmpmedian = -tmpmedian;
        }
        if (tmpmedian > pSeriesTable->maxmedian[pSeriesTable->seriesCount]) {
          pSeriesTable->maxmedian[pSeriesTable->seriesCount] = tmpmedian;
        }
        tmpmedian = ymedian;
        if (tmpmedian < 0.0) {
          tmpmedian = -tmpmedian;
        }
        if (tmpmedian > pSeriesTable->maxmedian[pSeriesTable->seriesCount]) {
          pSeriesTable->maxmedian[pSeriesTable->seriesCount] = tmpmedian;
        }
        if (xrms >  pSeriesTable->maxrms[pSeriesTable->seriesCount]) {
          pSeriesTable->maxrms[pSeriesTable->seriesCount] = xrms;
        }
        if (yrms >  pSeriesTable->maxrms[pSeriesTable->seriesCount]) {
          pSeriesTable->maxrms[pSeriesTable->seriesCount] = yrms;
        }
       

      } else {
        xmedian = 0.0;
        ymedian = 0.0;
        xrms = 99.0;
        yrms = 99.0;
      }
      fprintf(medianHandle,"%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n",
              pOldGridEntry->series,
              SAMPLE_POINTS+1,
              SAMPLE_POINTS+1,
              dix,
              diy,
              xmedian,
              ymedian);
      fprintf(rmsHandle,"%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n",
              pOldGridEntry->series,
              SAMPLE_POINTS+1,
              SAMPLE_POINTS+1,
              dix,
              diy,
              xrms,
              yrms);
    }
  }


  /* Now clear the counts for the next series */

  for (dix = 0; dix <= SAMPLE_POINTS; dix++) {
    for (diy = 0; diy <= SAMPLE_POINTS; diy++) {
      irec = dix + (diy * (SAMPLE_POINTS+1));
      pSeriesTable->curCount[irec] = 0;
    }
  }

  fclose(rmsHandle);
  fclose(medianHandle);
  fclose(countHandle);
}

void ProcessGridEntry(PSERIESTABLE pSeriesTable, PGRIDENTRY pOldGridEntry,PGRIDENTRY pGridEntry,int nrecs,int verbose)
{
  int series_result;
  int dix;
  int diy;
  int irec;
  series_result = strcmp(pOldGridEntry->series,pGridEntry->series);
  if (series_result < 0) {
    fprintf(stderr,"New series %s at line %d\n",pGridEntry->series,nrecs);
    if (strlen(pOldGridEntry->series) > 0) {
      if (strlen(pOldGridEntry->series) > (MAX_SERIES_STRING-1)) {
        fprintf(stderr,"ERROR: series length %d exceeds MAX_SERIES_STRING %d\n",
                strlen(pOldGridEntry->series),
                MAX_SERIES_STRING);
        exit(-1);
      }
      if (pSeriesTable->seriesCount > MAX_SERIES) {
        fprintf(stderr,"ERROR: MAX_SERIES exceeded\n");
        exit(-1);
      }
      strcpy(&pSeriesTable->seriesStrings[MAX_SERIES_STRING * pSeriesTable->seriesCount],pOldGridEntry->series);
      pSeriesTable->seriesPlates[pSeriesTable->seriesCount] = pSeriesTable->seriesPlateCount;
      pSeriesTable->seriesPlateCount = 0;
      ProcessSeries(pSeriesTable,pOldGridEntry,verbose);
      pSeriesTable->seriesCount++;
    }
    
    series_result = 0;
    pOldGridEntry->plateNumber = 0;

  } else if (series_result > 0) {
    fprintf(stderr,"ERROR: Series %s and %s are not in alphabetical order at record %d\n",
            pOldGridEntry->series,pGridEntry->series,nrecs);
    exit(-1);
  }
  if (series_result == 0) {
    if (pOldGridEntry->plateNumber != pGridEntry->plateNumber) {
      pSeriesTable->plateCount++;
      pSeriesTable->seriesPlateCount++;
      if (verbose) {
        fprintf(stderr,"New plate: %s%05d at record %d\n",pGridEntry->series,pGridEntry->plateNumber,nrecs);
      }
    }
  }
  if (strcmp(pGridEntry->series,"zzzzz") == 0) {
    /* We are done if we encounter a phoney series */
    return;
  }
  /* Here we validate the entry and store the results in the plate array */
  if ((pGridEntry->dnx != (SAMPLE_POINTS+1)) ||
      (pGridEntry->dny != (SAMPLE_POINTS+1)) ||
      (pGridEntry->dix < 0) ||
      (pGridEntry->dix > (SAMPLE_POINTS)) ||
      (pGridEntry->diy < 0) ||
      (pGridEntry->diy > (SAMPLE_POINTS))) {
    fprintf(stderr,"Entry at line %d has invalid dnx %d, dny %d,dix %d, or diy %d\n",nrecs,
            pGridEntry->dnx,
            pGridEntry->dny,
            pGridEntry->dix,
            pGridEntry->diy);
    exit(-1);
  }
  /* Now we transform the matrix to the original plate position.  For the momement, just do a simple 90 degree rotation */
  if (pGridEntry->width < pGridEntry->height) {
    dix = pGridEntry->dix;
    diy = pGridEntry->diy;
  } else {
    dix = pGridEntry->dix;
    diy = pGridEntry->diy;
  }
  irec = dix + (diy * (SAMPLE_POINTS+1));
  if (pSeriesTable->allocCount[irec] == 0) {
    pSeriesTable->dxvector[irec]  = (double *)calloc(ALLOC_INCREMENT,sizeof(double));
    pSeriesTable->dyvector[irec]  = (double *)calloc(ALLOC_INCREMENT,sizeof(double));
    if ((pSeriesTable->dxvector[irec] == NULL) || (pSeriesTable->dyvector[irec] == NULL)) {
      fprintf(stderr,"ERROR: failed to allocate dxvector or dyvector\n");
      exit(-1);
    }
    pSeriesTable->allocCount[irec] = ALLOC_INCREMENT;
  }


  if (pSeriesTable->curCount[irec] >= pSeriesTable->allocCount[irec]) {
    double * tmpvector;
    /* Need to realloc vector space */
    pSeriesTable->allocCount[irec] += ALLOC_INCREMENT;
    tmpvector = realloc(pSeriesTable->dxvector[irec],pSeriesTable->allocCount[irec]*sizeof(double));
    if (tmpvector == NULL) {
      fprintf(stderr,"ERROR failed to reallocate dxvector at size %d\n",pSeriesTable->allocCount[irec]);
      exit(-1);
    }
    pSeriesTable->dxvector[irec] = tmpvector;
    tmpvector = realloc(pSeriesTable->dyvector[irec],pSeriesTable->allocCount[irec]*sizeof(double));
    if (tmpvector == NULL) {
      fprintf(stderr,"ERROR failed to reallocate dyvector at size %d\n",pSeriesTable->allocCount[irec]);
      exit(-1);
    }
    pSeriesTable->dyvector[irec] = tmpvector;
  }
  pSeriesTable->dxvector[irec][pSeriesTable->curCount[irec]] = pGridEntry->dxval;
  pSeriesTable->dyvector[irec][pSeriesTable->curCount[irec]] = pGridEntry->dyval;
  pSeriesTable->curCount[irec]++;
  return;

}


int main(int argc,char *argv[])
{
  char *argstr;
  char listname[MAX_FILENAME];
  char mosaicname[MAX_FILENAME];
  char inputname[MAX_FILENAME];
  char *header = NULL;
  struct stat statbuf;
  FILE * listFile;
  FILE *outHandle = NULL;
  FILE *plotHandle = NULL;
  FILE *summaryHandle = NULL;
  int errorFlag = 0;
  int errorFlag2 = 0;
  int numPlates = 0;
  int numSeries = 0;
  int maxString = 0;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  int nvals;
  int lineLen;
  char* listStrings = NULL;
  int curline = 0;
  time_t startTime;
  time_t curTime;
  int verbose = 0;
  int verbose2 = 0;
  char cmdchar;
  int plateIndex;
  PPLATEENTRY pPlateEntry;
  PPLATEENTRY pPlateTable = NULL;
  char outfile[MAX_BUFFER];
  char plotfile[MAX_BUFFER];
  char summaryfile[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  struct WorldCoor *tanwcs = NULL;
  struct WorldCoor *tnxwcs = NULL;
  char *ingestDirectory;

  double xval;
  double yval;
  double xval2;
  double yval2;
  int offscl;
  int dix;
  int diy;
  int irec;
  double tnxra;
  double tnxdec;
  double tanra;
  double tandec;
  double xmin;
  double ymin;
  double xctr;
  double yctr;
  double dra;
  double ddec;
  int wp;
  int hp;
  double cra;
  double cdec;
  double era;
  double edec;
  double dist;
  double pixdist;
  double scale;
  double secpix;
  double cd[4];
  int programMode = MODE_ORIGINAL;
  int nx;                   /* dmagcor_table width */
  int ny;                   /* dmagcor_table height */
  int ix;
  int iy;
  int maxBin = MAX_SPATIAL_BINS;


  char *mysqlhost;
  char *username;
  char *password;
  MYSQL my_connection;
  MOSAIC mosaicInfo;
  PMOSAIC pMosaic = &mosaicInfo;
  int spatial_bin;
  double edgeDist;
  double maxDist;
  double distPixels;
  double maxDistPixels[MAX_SPATIAL_BINS+1];
  int sampleCount[MAX_SPATIAL_BINS+1];

  File spatial_bins_handle = NULL;
  char spatial_bins_name[MAX_BUFFER];
  PSPATIAL_BIN tmp_spatial_bin_table = NULL;
  TableHead spatial_bins_header = NULL;
  PSPATIAL_BIN pSpatial_bin;
  size_t spatial_bin_nrecs;
  int index;

  SPATIAL_BIN spatial_bin_table[MAX_SPATIAL_BINS+1];




  File grid_handle = NULL;
  char gridname[MAX_FILENAME];
  TableHead grid_header = NULL;
  int grid_nrecs = 0;
  int grid_index;
  GRIDENTRY grid_record;
  PGRIDENTRY pGridEntry = &grid_record;
  GRIDENTRY old_grid_record;
  PGRIDENTRY pOldGridEntry = &old_grid_record;
  TblDescriptor grid_descriptor = NULL;
  TableRow grid_row = NULL;

  File dmagcor_handle = NULL;
  char dmagcor_name[MAX_BUFFER];
  TableHead dmagcor_header = NULL;
  size_t dmagcor_nrecs;
  PDMAGCOR dmagcor_table = NULL;
  PDMAGCOR pDmagcor;


  PSERIESTABLE pSeriesTable = NULL;

  listname[0] = 0;
  outfile[0] = 0;
  plotfile[0] = 0;
  summaryfile[0] = 0;
  qualifier[0] = 0;
  gridname[0] = 0;
  memset(pOldGridEntry,0,sizeof(GRIDENTRY));
 
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      errorFlag = 1;
      printf("ERROR: unqualified argument %s\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'l': /* list file name */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for list file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(listname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'g': /* grid file name */
        case 'G':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(gridname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'm': /* program mode */
        case 'M':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&programMode);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode program mode %s\n",*argv);
              errorFlag = 1;
            } else {
              if ((programMode < MODE_ORIGINAL) || (programMode > MODE_MAX)) {
                errorFlag = 1;
                fprintf(stderr,"ERROR: Program mode must be an integer from %d to %d\n",MODE_ORIGINAL,MODE_MAX);
              }
            }
          }
          break;

        case 'b': /* maximum spatial bin */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&maxBin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode maxbin %s\n",*argv);
              errorFlag = 1;
            } else {
              if ((maxBin < 1) || (maxBin > MAX_SPATIAL_BINS)) {
                errorFlag = 1;
                fprintf(stderr,"ERROR: maximum bin must be an integer from 1 to %d\n",MAX_SPATIAL_BINS);
              }
            }
          }
          break;

        case 'i': /* input file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(inputname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'p': /* plot file name */
        case 'P':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(plotfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 's': /* summary file name */
        case 'S':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(summaryfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Validate arguments */

  if ((programMode == MODE_ORIGINAL) ||
      (programMode == MODE_CREATE_DISTORTION)) {
    if (outfile[0] == 0) {
      fprintf(stderr,"ERROR: No output filename was specified\n");
      errorFlag = 1;
    }


    outHandle = fopen(outfile,"wt");
    if (outHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output file %s\n",outfile);
      }
    }

  }

  if (plotfile[0] != 0) {
    plotHandle = fopen(plotfile,"wt");
    if (plotHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the plot file %s\n",plotfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Plot file %s\n",plotfile);
      }
    }
  }

  if ((programMode == MODE_READ_DISTORTION) ||
      (programMode == MODE_FIND_DISTORTION)) {
    if (summaryfile[0] != 0) {
      summaryHandle = fopen(summaryfile,"wt");
      if (summaryHandle == NULL) {
        errorFlag = 1;
        fprintf(stderr,"ERROR: Failed to open the summary file %s\n",summaryfile);
      } else {
        if (verbose) {
          fprintf(stderr,"Summary file %s\n",summaryfile);
        }
      }
    }
  }
  /* Attempt to open the list */

  if ((programMode == MODE_ORIGINAL) ||
      (programMode == MODE_CREATE_DISTORTION) ||
      (programMode == MODE_FIND_DISTORTION)) {
    if (listname[0] == 0) {
      printf("ERROR: No plate list specified\n");
      errorFlag = 1;
    }

    if (listname[0] != 0) {
      listFile = fopen(listname,"rt");
      if (listFile == NULL) {
        printf("Could not open file %s\n",listname);
        errorFlag = 1;
      } else {


        /* Count the number of records in the file */
        while (1) {
          inBuffer = fgets(inLine,MAX_BUFFER,listFile);
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
            numPlates++;
            if (maxString < lineLen) {
              maxString = lineLen;
            }
          }
  

        }
        fclose(listFile);
        if (numPlates == 0) {
          printf("Error: no lines in the list file\n");
          errorFlag = 1;
        }
      }
    } 
  }
  
  if (errorFlag) {
    printf("Usage: find_distortion -m <program mode> options\n");
    printf("  where program mode = 1: output file is maximum distortion, plot file is distortion vs plate distance\n");
    printf("        program mode = 2: output file is distortion grid\n");
    printf("        program mode = 3: input file is distortion grid\n");
    printf("        program mode = 4: combine steps 2 and 3\n");
    printf("  options: -v verbose\n");
    printf("           -l <Plate list file>\n");
    printf("           -i <input file>\n");
    printf("           -o <output file\n");
    printf("           -p <plot file>\n");
    printf("           -s <summary file>\n");
    printf("           -b <maximum bin>\n");
    return(-1);
  }

  printf("find_distortion of %s %s, mode %d maxBin %d minCount %d list file %s entries %d, max string %d\n",
         __DATE__,__TIME__,programMode,maxBin,MINIMUM_COUNT,listname,numPlates,maxString);
  fprintf(summaryHandle,"find_distortion of %s %s, mode %d maxBin %d minCount %d list file %s entries %d, max string %d\n",
          __DATE__,__TIME__,programMode,maxBin,MINIMUM_COUNT,listname,numPlates,maxString);
  if (programMode == MODE_CREATE_DISTORTION) {
    fprintf(outHandle,"series\tplateNumber\tspatial_bin\tdnx\tdny\tdix\tdiy\twidth\theight\txval\tyval\tdxval\tdyval\n");
    fprintf(outHandle,"------\t-----------\t-----------\t---\t---\t---\t---\t-----\t------\t----\t----\t-----\t-----\n");
  }
  if (programMode == MODE_ORIGINAL) {
    fprintf(outHandle,"series\tplateNumber\tmaxDist1\tmaxDist2\tmaxDist3\tmaxDist4\tmaxDist5\tmaxDist6\tmaxDist7\tmaxDist8\tmaxDist9\n");
    fprintf(outHandle,"------\t-----------\t--------\t--------\t--------\t--------\t--------\t--------\t--------\t--------\t--------\n");
  }
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



  /* Now allocate space for all of the records */
  if (listname[0] != 0) {
    maxString += 5;
    listStrings = calloc((numPlates+1)*maxString,sizeof(char));
    strcpy(&listStrings[numPlates*maxString],"UNKNOWN");
    if (listStrings == NULL) {
      printf("Failed to allocate listStrings\n");
    }
    pPlateTable = (PPLATEENTRY)calloc(numPlates,sizeof(PLATEENTRY));
    if (pPlateTable == NULL) {
      printf("Failed to allocate pPlateTable");
    }


    /* Attempt to open the list */
    listFile = fopen(listname,"rt");
    if (listFile == NULL) {
      printf("Could not open file %s\n",listname);
      return(-1);
    } else {

      /* Read the file records */
      while (1) {
        inBuffer = fgets(inLine,MAX_BUFFER,listFile);
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
          strcpy(&listStrings[curline*maxString],inBuffer);
          curline++;
        }
  
      }
      fclose(listFile);
    }
  }
  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_INGEST is not defined\n");
    return(-1);
  }
  time(&startTime);
  if ((programMode == MODE_READ_DISTORTION) ||
      (programMode == MODE_FIND_DISTORTION)) {
    pSeriesTable = calloc(1,sizeof(SERIESTABLE));
    if (pSeriesTable == NULL) {
      fprintf(stderr,"ERROR: Failed to allocate the series table \n");
      exit(-1);
    }

    strcpy(pSeriesTable->gridname,gridname);
    if (gridname[0] == 0) {
      fprintf(stderr,"ERROR: Input grid file name not specified\n");
      errorFlag = 1;
    }
  }


  for (plateIndex = 0; plateIndex < numPlates; plateIndex++) {

    pPlateEntry = &pPlateTable[plateIndex];
    strcpy(pPlateEntry->listName,&listStrings[plateIndex*maxString]);
    if (ParseFilename(pPlateEntry->listName,pPlateEntry->series,&pPlateEntry->plateNumber,&pPlateEntry->mosaicNumber,&pPlateEntry->binning,&pPlateEntry->rotation) == 0) {
      fprintf(stderr,"ERROR: Failed to parse the name  %s\n",pPlateEntry->listName);
      continue;
    }
    /* Note: we are using only solution 0, which presumably is the most accurate */
    if (GetMosaicInfo(&my_connection,pPlateEntry->series,pPlateEntry->plateNumber,pPlateEntry->mosaicNumber,0,pMosaic) != 1) {
      fprintf(stderr,"ERROR: failed to get mosaic information for %s\n",pPlateEntry->listName);
      continue;
    }
    if (pMosaic->rotation != pPlateEntry->rotation) {
      fprintf(stderr,"ERROR: filename rotation %d does not agree with database rotation %d for %s\n",pPlateEntry->rotation,pMosaic->rotation,pPlateEntry->listName);
      continue;
    }
    if (pMosaic->rotation == 0) {
      sprintf(mosaicname,"/dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d_01ww_tnx.fit",
              pMosaic->diskLocation,
              pPlateEntry->series,
              pPlateEntry->plateNumber,
              pPlateEntry->mosaicNumber,
              pPlateEntry->series,
              pPlateEntry->plateNumber,
              pPlateEntry->mosaicNumber);
    } else {
      sprintf(mosaicname,"/dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d_01r%dww_tnx.fit",
              pMosaic->diskLocation,
              pPlateEntry->series,
              pPlateEntry->plateNumber,
              pPlateEntry->mosaicNumber,
              pPlateEntry->series,
              pPlateEntry->plateNumber,
              pPlateEntry->mosaicNumber,
              pPlateEntry->rotation);

    }
    errorFlag = 0;
    errorFlag2 = 0;
    if ((header = GetFITShead (mosaicname, verbose2)) == NULL) {
      fprintf(stderr,"ERROR: Failed to read mosaic file header %s\n",mosaicname);
      errorFlag = 1;
      continue;
    } else {
     
      tnxwcs = GetFITSWCS (mosaicname, header, verbose2, &cra, &cdec, &dra, &ddec,
                           &secpix, &wp, &hp, &sysim, &eqim);
      if (nowcs (tnxwcs)) {
        fprintf(stderr,"ERROR: No valid WCS in the mosaic file %s\n",mosaicname);
        errorFlag = 1;
        continue;
      }
    }
    if ((programMode == MODE_CREATE_DISTORTION) ||
        (programMode == MODE_FIND_DISTORTION)) {
      /* Find out which spatial bins are valid for this file */
      if ((errorFlag == 0) && (summaryfile[0] != 0)) {
        for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
          spatial_bin_table[spatial_bin].spatial_bin = 0;
        }


        strcpy(spatial_bins_name,ingestDirectory);
        strcat(spatial_bins_name,"/");
        strcat(spatial_bins_name,pPlateEntry->listName);
        strcat(spatial_bins_name,qualifier);
        strcat(spatial_bins_name,".out.spatial_bins.db");
        spatial_bins_handle = Open(spatial_bins_name,"r");
        if (spatial_bins_handle == NULL) {
          fprintf(stderr,"ERROR: Failed to find spatial bins file %s\n",spatial_bins_name);
          errorFlag2 = 1;
         
        } else {
          if (verbose) {
            fprintf(stderr,"Found spatial bins file %s\n",spatial_bins_name);
          }
        }
      

        /* Now read in the spatial bins file */
        if (errorFlag2 == 0) {
          spatial_bins_header = table_header(spatial_bins_handle,TABLE_PARSE);
          if (spatial_bins_header == NULL) {
            fprintf(stderr,"ERROR: Failed to read header for %s\n",spatial_bins_name);
            errorFlag2 = 1;
          
          }
        }
        if (errorFlag2 == 0) {
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
                                               0,"end",0);
          if (tmp_spatial_bin_table == NULL) {
            fprintf(stderr,"ERROR: Failed to read table for %s\n",spatial_bins_name);
            errorFlag2 = 1;
          
          } else {
            if (verbose) {
              fprintf(stderr,"read %d records for %s\n",spatial_bin_nrecs,spatial_bins_name);
            }
          }
        }
        if (errorFlag2 == 0) {
          /* Now copy the spatial bin entries into their proper slots */
          for (index = 0; index < spatial_bin_nrecs; index++) {
            spatial_bin = tmp_spatial_bin_table[index].spatial_bin;
            if ((spatial_bin > 0) && (spatial_bin <= MAX_SPATIAL_BINS)) {
              memcpy(&spatial_bin_table[spatial_bin],&tmp_spatial_bin_table[index],sizeof(SPATIAL_BIN));
            } else {
              fprintf(stderr,"ERROR: Illegal spatial bin %d\n",spatial_bin);
            }
          }
        }
        if (tmp_spatial_bin_table != NULL) {
          Free(tmp_spatial_bin_table);
          tmp_spatial_bin_table = NULL;
        }
  
        if (spatial_bins_header != NULL) {
          table_hdrfree(spatial_bins_header);
          spatial_bins_header = NULL;
        }
        if (spatial_bins_handle != NULL) {
          Close(spatial_bins_handle);
          spatial_bins_handle = NULL;
        }
     



      }

    }
    if ((programMode == MODE_CREATE_DISTORTION) ||
        (programMode == MODE_FIND_DISTORTION)) {
      /* Find out which local smoothing bins are valid for this file */
      if ((errorFlag == 0) && (summaryfile[0] != 0)) {
        if (dmagcor_table != NULL) {
          Free(dmagcor_table);
          dmagcor_table = NULL;
        }
        if (dmagcor_header != NULL) {
          table_hdrfree(dmagcor_header);
          dmagcor_header = NULL;
        }

        if (dmagcor_handle != NULL) {
          Close(dmagcor_handle);
          dmagcor_handle = NULL;
        }

        strcpy(dmagcor_name,ingestDirectory);
        strcat(dmagcor_name,"/");
        strcat(dmagcor_name,pPlateEntry->listName);
        strcat(dmagcor_name,qualifier);
        strcat(dmagcor_name,"_dmagcor.grid");
        dmagcor_handle = Open(dmagcor_name,"r");
        if (dmagcor_handle == NULL) {
          fprintf(stderr,"ERROR: Failed to find the local calibration file %s\n",dmagcor_name);
          errorFlag2 = 1;
        } else {
          if (verbose) {
            fprintf(stderr,"Found local calibration file %s\n",dmagcor_name);
          }
        }
    
        /* Now read in the local calibration file */

        if (errorFlag2 == 0) {
          dmagcor_header = table_header(dmagcor_handle,TABLE_PARSE);
          if (dmagcor_header == NULL) {
            fprintf(stderr,"ERROR: Failed to read header for %s\n",dmagcor_name);
      
            errorFlag2 = 1;
          }
        }
        if (errorFlag2 == 0) {
          dmagcor_table = table_loadva(dmagcor_handle,
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
                                       TblInt,"rejectFlag",TblOff(PDMAGCOR,rejectFlag),
                                       TblInt,"npout",TblOff(PDMAGCOR,npout),
                                       TblDbl,"zout",TblOff(PDMAGCOR,zout),
                                       TblDbl,"errout",TblOff(PDMAGCOR,errout),
                                       TblDbl,"dradAverage",TblOff(PDMAGCOR,dradAverage),

                                       0,"end",0);
          if (dmagcor_table == NULL) {
            fprintf(stderr,"ERROR: Failed to read table for %s\n",dmagcor_name);
            errorFlag2 = 1;
          }
 
        }
      }
    }
    if (errorFlag == 0) {



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



      tanwcs = wcskinit(pMosaic->naxis1,
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


      if (nowcs (tanwcs)) {
        fprintf(stderr,"ERROR: No valid WCS in the database for %s\n",pPlateEntry->listName);
        errorFlag = 1;
      }
      if (tanwcs->lngcor != NULL) {
        wf_gsclose(tanwcs->lngcor);
        tanwcs->lngcor = NULL;
      }
      if (tanwcs->latcor != NULL) {
        wf_gsclose(tanwcs->latcor);
        tanwcs->latcor = NULL;
      }

    }
    if (errorFlag == 0) {
      xctr = 0.5 + (0.5 *pMosaic->naxis1);
      yctr = 0.5 + (0.5 *pMosaic->naxis2);
      pix2wcs(tanwcs,xctr,yctr,&cra,&cdec);
      xmin = 0.5;
      ymin = 0.5;
      pix2wcs(tanwcs,xmin,ymin,&era,&edec);

      dist = wcsdist(cra,cdec,era,edec);

      pixdist = sqrt(((xctr-xmin)*(xctr-xmin)) + ((yctr-ymin)*(yctr-ymin)));
      if (pixdist > 0.0) {
        scale = dist/pixdist;
      } else {
        fprintf(stderr,"ERROR: No plate scale for %s\n",pPlateEntry->listName);
        errorFlag = 1;
      }
    }

    if (errorFlag == 0) {
      for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
        maxDistPixels[spatial_bin] = 0.0;
        sampleCount[spatial_bin] = 0;
      }

      for (dix = 0; dix <= SAMPLE_POINTS; dix++) {
        xval = (pMosaic->naxis1/SAMPLE_POINTS)*dix;
        for (diy = 0; diy <= SAMPLE_POINTS; diy++) {
          yval = (pMosaic->naxis2/SAMPLE_POINTS)*diy;
          if ((programMode == MODE_CREATE_DISTORTION) ||
              (programMode == MODE_FIND_DISTORTION)) {
            spatial_bin = CalculateBin(pMosaic->naxis1,pMosaic->naxis2,xval,yval,&edgeDist);



            if (spatial_bin == 0) {
              fprintf(stderr,"ERROR: spatial bin is %d for %f %f\n",spatial_bin,xval,yval);
              exit(-1);
            }
            if (errorFlag2 != 0) {
              /* No spatial bin or dmagcor tables available */
              continue;
            }
            if (spatial_bin != spatial_bin_table[spatial_bin].spatial_bin) {
              /* Not in a legal spatial bin */
              continue;
            }
            if (spatial_bin > maxBin) {
              /* Excluded bin */
              continue;
            }
            nx = dmagcor_table[0].nx;
            ny = dmagcor_table[0].ny;
            ix = (xval * nx) /(1.0 * pMosaic->naxis1);
            iy = (yval * ny) /(1.0 * pMosaic->naxis2);
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
            /* These records are reversed from the local bin */
            irec = iy + (ny*ix);
            pDmagcor = &dmagcor_table[irec];
            if ((pDmagcor->nx != nx) ||
                (pDmagcor->ny != ny) ||
                (pDmagcor->ix != ix) ||
                (pDmagcor->iy != iy)) {
              fprintf(stderr,"ERROR: invalid dmagcor table ix %d iy %d nx %d ny %d for %s\n",
                      ix,iy,nx,ny,dmagcor_name);
            }
            if ((pDmagcor->errout >= 90.0) ||
                (pDmagcor->rejectFlag != 0)) {
              /* This is a bad local smoothing bin */
              continue;
            }   

            sampleCount[spatial_bin]++;
            /* Find the TNX RA AND DEC */
            pix2wcs(tnxwcs,xval,yval,&tnxra,&tnxdec);
            /* Now find the TAN reverse values */
            wcs2pix(tanwcs,tnxra,tnxdec,&xval2,&yval2,&offscl);
            if (offscl) {
              /* Can not get inverse coordinate */
              continue;
            }
            xval2 -= xval;
            yval2 -= yval;
            
            if (programMode == MODE_CREATE_DISTORTION) {
              fprintf(outHandle,"%s\t%05d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n",
                      pPlateEntry->series,
                      pPlateEntry->plateNumber,
                      spatial_bin,
                      SAMPLE_POINTS+1,
                      SAMPLE_POINTS+1,
                      dix,
                      diy,
                      pMosaic->naxis1,
                      pMosaic->naxis2,
                      xval,
                      yval,
                      xval2,
                      yval2);
            }
            if (programMode == MODE_FIND_DISTORTION) {
              memset(pGridEntry,0,sizeof(GRIDENTRY));

              strcpy(pGridEntry->series,pPlateEntry->series);
              
              pGridEntry->plateNumber = pPlateEntry->plateNumber;
              pGridEntry->spatial_bin = spatial_bin;
              pGridEntry->dnx         = SAMPLE_POINTS+1;
              pGridEntry->dny         = SAMPLE_POINTS+1;
              pGridEntry->dix         = dix;
              pGridEntry->diy         = diy;
              pGridEntry->width       = pMosaic->naxis1;
              pGridEntry->height      = pMosaic->naxis2;
              pGridEntry->xval        = xval;
              pGridEntry->yval        = yval;
              pGridEntry->dxval       = xval2;
              pGridEntry->dyval       = yval2;

              ProcessGridEntry(pSeriesTable,pOldGridEntry,pGridEntry,grid_nrecs,verbose);
              memcpy(pOldGridEntry,pGridEntry,sizeof(GRIDENTRY));
            }

          }
          if (programMode == MODE_ORIGINAL) {
            pix2wcs(tanwcs,xval,yval,&tanra,&tandec);
            pix2wcs(tnxwcs,xval,yval,&tnxra,&tnxdec);
            dist = wcsdist(tanra,tandec,tnxra,tnxdec);
            distPixels = dist/scale;
            pixdist = sqrt(((xctr-xval)*(xctr-xval)) + ((yctr-yval)*(yctr-yval)));

            spatial_bin = CalculateBin(pMosaic->naxis1,pMosaic->naxis2,xval,yval,&edgeDist);
            if (spatial_bin == 0) {
              fprintf(stderr,"ERROR: spatial bin is %d for %f %f\n",spatial_bin,xval,yval);
            } else {
              sampleCount[spatial_bin]++;
              if (distPixels > maxDistPixels[spatial_bin]) {
                maxDistPixels[spatial_bin] = distPixels;
              }
            }


            if (plotHandle != 0) {
              fprintf(plotHandle,"%.0f %f\n",pixdist,distPixels);


            }

          }
        
        }
      }

      if (programMode == MODE_ORIGINAL) {
        fprintf(outHandle,"%s\t%05d",pPlateEntry->series,pPlateEntry->plateNumber);
        for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
          fprintf(outHandle,"\t%.2f",maxDistPixels[spatial_bin]);
        }
        fprintf(outHandle,"\n");
      }
      pPlateEntry->fileExists = 1;
      printf("Spatial Bin Counts");
      for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
        printf(" %4d",sampleCount[spatial_bin]);
      }
     
      printf(" for Plate: %s\n",pPlateEntry->listName);

    }
    if (tnxwcs != NULL) {
      wcsfree(tnxwcs);
      tnxwcs = NULL;
    }
    if (header != NULL) {
      free(header);
    }
    if (tanwcs != NULL) {
      wcsfree(tanwcs);
      tanwcs = NULL;
    }
  }
  

  if (programMode == MODE_READ_DISTORTION) {
    grid_handle = Open(gridname,"rt");
    if (grid_handle == NULL) {
      fprintf(stderr,"ERROR: Could not open file %s\n",gridname);
      return(-1);
    } 
  

    /* Now read in the grid file one record at a time */

    grid_header = table_header(grid_handle,TABLE_PARSE);
    if (grid_header == NULL) {
      fprintf(stderr,"ERROR: Failed to read header for %s\n",gridname);
      return(-1);
    }

    grid_descriptor = table_create_descrip(&grid_nrecs,
                                           TblBuf,"series"    ,TblOff(PGRIDENTRY,series),MAX_SERIES_STRING,
                                           TblDbl,"plateNumber",TblOff(PGRIDENTRY,dPlateNumber),
                                           TblInt,"spatial_bin",TblOff(PGRIDENTRY,spatial_bin),
                                           TblInt,"dnx",TblOff(PGRIDENTRY,dnx),
                                           TblInt,"dny",TblOff(PGRIDENTRY,dny),
                                           TblInt,"dix",TblOff(PGRIDENTRY,dix),
                                           TblInt,"diy",TblOff(PGRIDENTRY,diy),
                                           TblInt,"width",TblOff(PGRIDENTRY,width),
                                           TblInt,"height",TblOff(PGRIDENTRY,height),
                                           TblDbl,"xval"  ,TblOff(PGRIDENTRY,xval),
                                           TblDbl,"yval"  ,TblOff(PGRIDENTRY,yval),
                                           TblDbl,"dxval"  ,TblOff(PGRIDENTRY,dxval),
                                           TblDbl,"dyval"  ,TblOff(PGRIDENTRY,dyval),
                                           0,"end",0);
    if (grid_descriptor == NULL) {
      fprintf(stderr,"ERROR: Failed to allocate descriptor for %s\n",gridname);
      return(-1);
    }
    table_loadmap(grid_header,grid_descriptor);






    while (1) {
      memset(pGridEntry,0,sizeof(GRIDENTRY));
      if(ReadGridEntry(pGridEntry,grid_handle,grid_header,grid_descriptor,&grid_row,gridname) == 0) {
        break;
      }
      grid_nrecs++;
      ProcessGridEntry(pSeriesTable,pOldGridEntry,pGridEntry,grid_nrecs,verbose);
      memcpy(pOldGridEntry,pGridEntry,sizeof(GRIDENTRY));
    }
    
  }
  
  if ((programMode == MODE_READ_DISTORTION) ||
      (programMode == MODE_FIND_DISTORTION)) {
    /* Handle the last record with a phoney series */
    strcpy(pGridEntry->series,"zzzzz");
    pGridEntry->plateNumber = 1;
    ProcessGridEntry(pSeriesTable,pOldGridEntry,pGridEntry,grid_nrecs,verbose);
    numPlates = pSeriesTable->plateCount;
    numSeries = pSeriesTable->seriesCount;
  }
  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds for %d mosaics %d series\n",curTime,numPlates,numSeries);
  /* Done with the MySQL connection */
  mysql_close(&my_connection);


  if ((programMode == MODE_READ_DISTORTION) ||
      (programMode == MODE_FIND_DISTORTION)) {
    for (index = 0; index < pSeriesTable->seriesCount; index++) {
      fprintf(stderr,"%3d series %5s  plates %5d maxmedian %6.2f maxrms %7.2f\n",index,
              &pSeriesTable->seriesStrings[MAX_SERIES_STRING * index],
              pSeriesTable->seriesPlates[index],
              pSeriesTable->maxmedian[index],
              pSeriesTable->maxrms[index]);
      if (summaryHandle != NULL) {
        fprintf(summaryHandle,"%3d series %5s  plates %5d maxmedian %6.2f maxrms %7.2f\n",index,
                &pSeriesTable->seriesStrings[MAX_SERIES_STRING * index],
                pSeriesTable->seriesPlates[index],
                pSeriesTable->maxmedian[index],
                pSeriesTable->maxrms[index]);

      }
    }
  }
  if (listStrings != NULL) {
    free(listStrings);
  }
  if (pPlateTable != NULL) {
    free(pPlateTable);
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }
  if (plotHandle != NULL) {
    fclose(plotHandle);
  }
  if (summaryHandle != NULL) {
    fclose(summaryHandle);
  }

  if (grid_header != NULL) {
    table_hdrfree(grid_header);
  }
  if (grid_handle != NULL) {
    Close(grid_handle);
  }
  if (grid_descriptor != NULL) {
    Free(grid_descriptor);
  }
  if (grid_row != NULL) {
    table_rowfree(grid_row);
  }
  if (dmagcor_table != NULL) {
    Free(dmagcor_table);
    dmagcor_table = NULL;
  }
  if (dmagcor_header != NULL) {
    table_hdrfree(dmagcor_header);
    dmagcor_header = NULL;
  }

  if (dmagcor_handle != NULL) {
    Close(dmagcor_handle);
    dmagcor_handle = NULL;
  }
  
  for (dix = 0; dix <= SAMPLE_POINTS; dix++) {
    for (diy = 0; diy <= SAMPLE_POINTS; diy++) {
      irec = dix + (diy * (SAMPLE_POINTS+1));
      if (pSeriesTable->dxvector[irec] != NULL) {
        free(pSeriesTable->dxvector[irec]);
      }
      if (pSeriesTable->dyvector[irec] != NULL) {
        free(pSeriesTable->dyvector[irec]);
      }
    }
  }
  if (pSeriesTable != NULL) {
    free(pSeriesTable);
    pSeriesTable = NULL;
  }

  return(EXIT_SUCCESS);
}

