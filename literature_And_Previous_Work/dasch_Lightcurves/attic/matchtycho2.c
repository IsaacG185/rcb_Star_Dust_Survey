// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* matchtycho2.c
 *
 *  Match entries in the sextractor table with entries in the tycho-2 catalog
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio -lmysqlclient  matchtycho2.c pipelineutils.a -ltable -lutil  -lwcs -o matchtycho2
 *
 *  imty2  -w -d -n 5000000 /dasch/data/ExposureData/Mosaics/dnb/03393_00/dnb03393_00_01ww.fit 
 *  mv dnb03393_00_01ww.fit.tycho2 /dasch/Pipeline/match
 *  matchtycho2 -v -w 17493 -h 21953 -c /dasch/Pipeline/match/dnb03393_00_01ww.fit.tycho2 -i /dasch/Pipeline/match/dnb03393_00_01ww.db -o /dasch/Pipeline/match/match_dnb03393_00_01ww_u.db -p 20 -l 101 -r 333 -b 130 -t 130 -s 5.902 -n 5000
 *  rm  /dasch/Pipeline/match/dnb03393_00_01ww.fit.tycho2
 *
 *
 *  imty2  -w -d -n 5000000 /dasch/data/ExposureData/Mosaics/ac/42227_00/ac42227_00_01ww.fit 
 *  mv ac42227_00_01ww.fit.tycho2 /dasch/Pipeline/match
 *  matchtycho2 -v -w 17493 -h 21953 -c /dasch/Pipeline/match/ac42227_00_01ww.fit.tycho2 -i /dasch/Pipeline/match/ac42227_00_01ww.db -o /dasch/Pipeline/match/match_ac42227_00_01ww_u.db -p 50 -l 101 -r 333 -b 130 -t 130 -s 6.302 -n 5000
 *  rm  /dasch/Pipeline/match/ac42227_00_01ww.fit.tycho2
 *
 *
 * Nov 17, 2008  Edward J. Los - Initial version
 * Dec  1, 2008  Edward J. Los - Correct segmentation fault in the search algorithm
 * Jan 23, 2008  Edward J. Los - split FLAGS into AFLAGS and BFLAGS
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 */


#include <math.h>
#include "table.h"
#include "time.h"
#include "pipelineutils.h"
#define MAX_FILENAME 256
#define MAX_BUFFER 256
#define MAX_LIST_STRING 25
#define DUPLICATE_FACTOR 2 /* ratio of sextractor stars to be matched with tycho2 stars */
/* #define VERIFY_TABLE 1 */


typedef struct _starimage {
  int reference_count; /* Number of times this image is referenced by a catalog object */
  int NUMBER;
  int spatial_bin;
  double FLUX_ISO;
  double MAG_ISO;
  double MAGERR_ISO;
  double MAG_APER;
  double MAGERR_APER;
  double MAG_AUTO;
  double MAGERR_AUTO;
  double KRON_RADIUS;
  double BACKGROUND;
  double THRESHOLD;
  double FLUX_MAX;
  int ISOAREA_IMAGE;
  double ISOAREA_WORLD;
  double X_IMAGE;
  double Y_IMAGE;
  double ra;
  double dec;
  double THETA_J2000;
  double ELLIPTICITY;
  double ERRTHETA_J2000;
  double FWHM_IMAGE;
  double FWHM_WORLD;
  int ISO0;
  int ISO1;
  int ISO2;
  int ISO3;
  int ISO4;
  int ISO5;
  int ISO6;
  int ISO7;
  int AFLAGS;
  int BFLAGS;
  double plate_dra;
  double plate_ddec;
  double plate_dist;

} SEXTRACTOR,*PSEXTRACTOR;

typedef struct _catalog {
  double ra;
  double dec;
  double magb;
  double x;
  double y;
  int spatial_bin;
  PSEXTRACTOR pSextractor;
} CATALOG,*PCATALOG;




int ReadSextractor(PSEXTRACTOR pSextractor,
               File sextractor_handle,
               TableHead sextractor_header,
               TblDescriptor sextractor_descriptor,
               TableRow* sextractor_row)
{
  *sextractor_row = table_rowget(sextractor_handle,sextractor_header,*sextractor_row,NULL,NULL,0);
  if (*sextractor_row == NULL) {
    return(0);
  }
  if (!table_loadrow(sextractor_handle,sextractor_header,*sextractor_row,sextractor_descriptor,(char *)pSextractor)) {
    printf("ERROR: Read Sextractor table_loadrow failed\n");
    return(0);
  }
  return(1);
}

/* Sort in ascending MAG_ISO */
int SextractorImageCompare(const void *first, const void *second) 
{
  double mag_isoFirst = ((PSEXTRACTOR)first)->MAG_ISO;
  double mag_isoSecond = ((PSEXTRACTOR)second)->MAG_ISO;
  if (mag_isoFirst > mag_isoSecond) {
    return(1);
  } else if (mag_isoFirst < mag_isoSecond) {
    return(-1);
  } else {
    return(0);
  }

}
int CatalogImageCompare(const void *first, const void *second) 
{
  double magbFirst = ((PCATALOG)first)->magb;
  double magbSecond = ((PCATALOG)second)->magb;
  if (magbFirst > magbSecond) {
    return(1);
  } else if (magbFirst < magbSecond) {
    return(-1);
  } else {
    return(0);
  }

}


/* Sort in ascending Y_IMAGE */
int SextractorYCompare(const void *first, const void *second) 
{
  double y_imageFirst = ((PSEXTRACTOR)first)->Y_IMAGE;
  double y_imageSecond = ((PSEXTRACTOR)second)->Y_IMAGE;
  if (y_imageFirst > y_imageSecond) {
    return(1);
  } else if (y_imageFirst < y_imageSecond) {
    return(-1);
  } else {
    return(0);
  }

}
int CatalogYCompare(const void *first, const void *second) 
{
  double yFirst = ((PCATALOG)first)->y;
  double ySecond = ((PCATALOG)second)->y;
  if (yFirst > ySecond) {
    return(1);
  } else if (yFirst < ySecond) {
    return(-1);
  } else {
    return(0);
  }

}


int main(int argc,char *argv[])
{
  char *argstr;
  FILE *outHandle = NULL;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int verbose = 0;
  char cmdchar;
  char outfile[MAX_BUFFER];
  char rootname[MAX_BUFFER];
  char *dotPtr;
  char *matchDirectory;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int pixelRadius = 0;
  int maxStars = 0;
  int nvals;
  int skipMargin = 0;
  int leftMargin = -1;
  int rightMargin = -1;
  int topMargin = -1;
  int bottomMargin = -1;
  double leftMarginD;
  double rightMarginD;
  double bottomMarginD;
  double topMarginD;
  int bin_output_count[MAX_SPATIAL_BINS+1];
  double *binVector[MAX_SPATIAL_BINS+1];
  int catalog_bin_count[MAX_SPATIAL_BINS+1];
  int new_catalog_bin_count[MAX_SPATIAL_BINS+1];
  int catalog_bin_index[MAX_SPATIAL_BINS+1];
  int sextractor_bin_count[MAX_SPATIAL_BINS+1];
  int new_sextractor_bin_count[MAX_SPATIAL_BINS+1];
  int sextractor_bin_index[MAX_SPATIAL_BINS+1]; /* Highest index value for a given bin */
  int double_match_bin[MAX_SPATIAL_BINS+1];
  int double_reference_bin[MAX_SPATIAL_BINS+1];
  int spatial_bin;
  int maxCatalogBinCount = 0;
  double edgeDist;
  int marginCount = 0;
  int doubleMatchCount = 0;
  int outCount = 0;
  int doubleReferenceCount = 0;
  int trialCount = 0;
  double factor;
  double dra;
  double ddec;
  double drad;
  double binMedian;
  double binRms;

  File catalog_handle = NULL;
  char catalog_name[MAX_BUFFER];
  TableHead catalog_header = NULL;
  size_t catalog_nrecs;
  int new_catalog_nrecs = 0;
  int catalog_index;
  int catalog_index2;
  PCATALOG catalog_table = NULL;
  PCATALOG pCatalog;
  PCATALOG pCatalog2;
  CATALOG catalogEntry;
  PCATALOG pCatalog3 = &catalogEntry;
  
  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  PSEXTRACTOR sextractor_table = NULL;
  size_t sextractor_nrecs = 0;
  int new_sextractor_nrecs = 0;
  int sextractor_index;
  int sextractor_index2;
  PSEXTRACTOR pSextractor;
  PSEXTRACTOR pSextractor2;
  SEXTRACTOR sextractorEntry;
  PSEXTRACTOR pSextractor3 = &sextractorEntry;
  TblDescriptor sextractor_descriptor = NULL;
  TableRow sextractor_row = NULL;

  int foundReplacement = 1;
  int new_bin_count;
  int max_bin_count;

  int sextractor_min_index;
  double xDist;
  double yDist;
  double totalDist;
  double arcsecPerPixel = 0.0;

  sextractor_name[0] = 0;
  outfile[0] = 0;
  catalog_name[0] = 0;

  time(&startTime);

  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    bin_output_count[spatial_bin] = 0;
    binVector[spatial_bin] = NULL;
    catalog_bin_count[spatial_bin] = 0;
    new_catalog_bin_count[spatial_bin] = 0;
    catalog_bin_index[spatial_bin] = -1;
    sextractor_bin_count[spatial_bin] = 0;
    new_sextractor_bin_count[spatial_bin] = 0;
    sextractor_bin_index[spatial_bin] = -1;
    double_match_bin[spatial_bin] = 0;
    double_reference_bin[spatial_bin] = 0;
  }

  matchDirectory = getenv("DASCH_MATCH");
  if (matchDirectory == NULL) {
    fprintf(stderr,"ERROR: DASCH_MATCH is not defined\n");
    return(-1);
  }


  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      fprintf(stderr,"ERROR Unrecognized argument %s\n",argstr);
      errorFlag = 1;
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
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

        case 'c': /* catalog file name */
        case 'C':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(catalog_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'i': /* sextractor file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sextractor_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 's': /* arcsec/pixel */
        case 'S':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&arcsecPerPixel);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Can not decode number of arcsec per pixel\n");
              errorFlag = 1;
            }
          }
          break;



        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicWidth);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode mosaic width %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* mosaic height */
        case 'H':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicHeight);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'p': /* search radius in pixels */
        case 'P':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&pixelRadius);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the search radius %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'n': /* maximum number of stars to match */
        case 'N':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&maxStars);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the search radius %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;
        case 'l': /* left margin */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&leftMargin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode left margin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'r': /* right margin */
        case 'R':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&rightMargin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode right margin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'b': /* bottom margin */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&bottomMargin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode bottom margin %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 't': /* top margin */
        case 'T':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&topMargin);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode top margin %s\n",*argv);
              errorFlag = 1;
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

  if (outfile[0] == 0) {
    fprintf(stderr,"ERROR: No output filename was specified\n");
    errorFlag = 1;
  }

  if (catalog_name[0] == 0) {
    fprintf(stderr,"ERROR: No catalog filename was specified\n");
    errorFlag = 1;
  }

  if (sextractor_name[0] == 0) {
    fprintf(stderr,"ERROR: No sextractor filename was specified\n");
    errorFlag = 1;
  }
  if (arcsecPerPixel <= 0.0) {
    fprintf(stderr,"ERROR: No pixel scale was specified\n");
    errorFlag = 1;
  }

  if (mosaicWidth <= 0) {
    fprintf(stderr,"ERROR: No mosaic width was specified\n");
    errorFlag = 1;
  }
  if (mosaicHeight <= 0) {
    fprintf(stderr,"ERROR: No mosaic height was specified\n");
    errorFlag = 1;
  }
  if (pixelRadius <= 0) {
    fprintf(stderr,"ERROR: No pixel radius was specified\n");
    errorFlag = 1;
  }
  if (maxStars <= 0) {
    fprintf(stderr,"ERROR: No maximum number of stars was specified\n");
    errorFlag = 1;
  }
    
  if ((leftMargin  == -1) &&
      (rightMargin == -1) &&
      (bottomMargin == -1) &&
      (topMargin == -1)) {
    skipMargin = 1;
  } else {
    if ((leftMargin  == -1) ||
        (rightMargin == -1) ||
        (bottomMargin == -1) ||
        (topMargin == -1)) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: One of leftMargin %d, rightMargin %d, bottomMargin %d, or topMargin %d not set\n",
              leftMargin,rightMargin,bottomMargin,topMargin);
    }
  }

  if (skipMargin == 0) {
    /* Set our margin limits */
    leftMarginD = 1.0*leftMargin;
    rightMarginD = 1.0*(mosaicWidth-rightMargin);
    bottomMarginD = 1.0*bottomMargin;
    topMarginD = 1.0 *(mosaicHeight-topMargin);

  }



  /* Open the sextractor file */
  sextractor_handle = Open(sextractor_name,"r");
  if (sextractor_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the sextractor file %s\n",sextractor_name);
  } else {
    if (verbose) {
      fprintf(stderr,"Found sextractor file %s\n",sextractor_name);
    }
  }

  /* Open the catalog file */
  catalog_handle = Open(catalog_name,"r");
  if (catalog_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the catalog file %s\n",catalog_name);
  } else {
    if (verbose) {
      fprintf(stderr,"Found catalog file %s\n",catalog_name);
    }
  }


  outHandle = fopen(outfile,"wt");
  if (outHandle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
  } else {
    fprintf(outHandle,"NUMBER\tMAG_ISO\tX_IMAGE\tY_IMAGE\tra_1\tdec_1\tFWHM_IMAGE\tFWHM_WORLD\tAFLAGS\tBFLAGS\tREF\tra_2\tdec_2\tStdmag\tcolor\tVFlag\tMAGFlag\tclass\tdra\tddec\tdrad\tTHRESHOLD\tFLUX_MAX\n");
    fprintf(outHandle,"------\t-------\t-------\t-------\t----\t-----\t----------\t----------\t------\t------\t---\t----\t-----\t------\t-----\t-----\t-------\t-----\t---\t----\t----\t---------\t--------\n");
    if (verbose) {
      fprintf(stderr,"Output file %s\n",outfile);
    }
  }

  dotPtr = strrchr(sextractor_name,'/');
  if (dotPtr == NULL) {
    strcpy(rootname,sextractor_name);
  } else {
    dotPtr++;
    strcpy(rootname,dotPtr);
  }
  dotPtr = strstr(rootname,".");
  if (dotPtr != NULL) {
    *dotPtr = 0;
  }
  
  if (errorFlag) {
    fprintf(stderr,"Usage: matchtycho2 [arguments]\n");
    fprintf(stderr,"       where [-v] verbose\n");
    fprintf(stderr,"             -c <catalog file from imty2>\n");
    fprintf(stderr,"             -o <match output file>\n");
    fprintf(stderr,"             -i <sextractor file> \n");
    fprintf(stderr,"             -w <mosaic width in pixels>\n");
    fprintf(stderr,"             -h <mosaic height in pixels\n");
    fprintf(stderr,"             -l <left margin in pixels> \n");
    fprintf(stderr,"             -r <right margin in pixels> \n");
    fprintf(stderr,"             -b <bottom margin in pixels> \n");
    fprintf(stderr,"             -t <top margin in pixels> \n");
    fprintf(stderr,"             -p <match radius in pixels\n"); 
    fprintf(stderr,"             -n <maximum number of stars to match\n"); 
    fprintf(stderr,"             -s <pixel scale in arcsec/pixel>\n");
    return(-1);
  }
  printf("matchtycho2 width %d height %d radius %d maxStars %d skipMargin %d margins %d %d %d %d for %s\n",mosaicWidth,mosaicHeight,pixelRadius,maxStars,skipMargin,leftMargin,rightMargin,bottomMargin,topMargin,rootname);

    
  /* Now read in the catalog file */

  catalog_header = table_header(catalog_handle,TABLE_PARSE);
  if (catalog_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",catalog_name);
    exit(-1);
  }

  catalog_table = table_loadva(catalog_handle,
                               &catalog_header,
                               NULL, /* hbase */
                               NULL, /* rows */
                               NULL,
                               sizeof(CATALOG),
                               &catalog_nrecs,
                               TblDbl,"ra",TblOff(PCATALOG,ra),
                               TblDbl,"dec",TblOff(PCATALOG,dec),
                               TblDbl,"magb",TblOff(PCATALOG,magb),
                               TblDbl,"x",TblOff(PCATALOG,x),
                               TblDbl,"y",TblOff(PCATALOG,y),
                               0,"end",0);

  if (catalog_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",catalog_name);
    exit(-1);
  }


  if (verbose) {
    fprintf(stderr,"read %d records for %s\n",catalog_nrecs,catalog_name);
  }


  /* Now read in the sextractor file */

  sextractor_header = table_header(sextractor_handle,TABLE_PARSE);
  if (sextractor_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",sextractor_name);
    exit(-1);
  }

  sextractor_table = table_loadva(sextractor_handle,
                                  &sextractor_header,
                                  NULL, /* hbase */
                                  NULL, /* rows */
                                  NULL,
                                  sizeof(SEXTRACTOR),
                                  &sextractor_nrecs,
                                  TblInt,"NUMBER",TblOff(PSEXTRACTOR,NUMBER),
#ifndef MINIMUM_DEBUG
                                  TblDbl,"FLUX_ISO",TblOff(PSEXTRACTOR,FLUX_ISO),
                                  TblDbl,"MAG_ISO",TblOff(PSEXTRACTOR,MAG_ISO),
                                  TblDbl,"MAGERR_ISO",TblOff(PSEXTRACTOR,MAGERR_ISO),
                                  TblDbl,"MAG_APER",TblOff(PSEXTRACTOR,MAG_APER),
                                  TblDbl,"MAGERR_APER",TblOff(PSEXTRACTOR,MAGERR_APER),
                                  TblDbl,"MAG_AUTO",TblOff(PSEXTRACTOR,MAG_AUTO),
                                  TblDbl,"MAGERR_AUTO",TblOff(PSEXTRACTOR,MAGERR_AUTO),
                                  TblDbl,"KRON_RADIUS",TblOff(PSEXTRACTOR,KRON_RADIUS),
                                  TblDbl,"BACKGROUND",TblOff(PSEXTRACTOR,BACKGROUND),
                                  TblDbl,"THRESHOLD",TblOff(PSEXTRACTOR,THRESHOLD),
                                  TblDbl,"FLUX_MAX",TblOff(PSEXTRACTOR,FLUX_MAX),
                                  TblInt,"ISOAREA_IMAGE",TblOff(PSEXTRACTOR,ISOAREA_IMAGE),
                                  TblDbl,"ISOAREA_WORLD",TblOff(PSEXTRACTOR,ISOAREA_WORLD),
#endif /* MINIMUM_DEBUG */
                                  TblDbl,"X_IMAGE",TblOff(PSEXTRACTOR,X_IMAGE),
                                  TblDbl,"Y_IMAGE",TblOff(PSEXTRACTOR,Y_IMAGE),
                                  TblDbl,"ra",TblOff(PSEXTRACTOR,ra),
                                  TblDbl,"dec",TblOff(PSEXTRACTOR,dec),
#ifndef MINIMUM_DEBUG
                                  TblDbl,"THETA_J2000",TblOff(PSEXTRACTOR,THETA_J2000),
                                  TblDbl,"ELLIPTICITY",TblOff(PSEXTRACTOR,ELLIPTICITY),
                                  TblDbl,"ERRTHETA_J2000",TblOff(PSEXTRACTOR,ERRTHETA_J2000),
                                  TblDbl,"FWHM_IMAGE",TblOff(PSEXTRACTOR,FWHM_IMAGE),
                                  TblDbl,"FWHM_WORLD",TblOff(PSEXTRACTOR,FWHM_WORLD),
                                  TblInt,"ISO0",TblOff(PSEXTRACTOR,ISO0),
                                  TblInt,"ISO1",TblOff(PSEXTRACTOR,ISO1),
                                  TblInt,"ISO2",TblOff(PSEXTRACTOR,ISO2),
                                  TblInt,"ISO3",TblOff(PSEXTRACTOR,ISO3),
                                  TblInt,"ISO4",TblOff(PSEXTRACTOR,ISO4),
                                  TblInt,"ISO5",TblOff(PSEXTRACTOR,ISO5),
                                  TblInt,"ISO6",TblOff(PSEXTRACTOR,ISO6),
                                  TblInt,"ISO7",TblOff(PSEXTRACTOR,ISO7),
                                  TblInt,"AFLAGS",TblOff(PSEXTRACTOR,AFLAGS),
                                  TblInt,"BFLAGS",TblOff(PSEXTRACTOR,BFLAGS),
                                  TblDbl,"plate_dra",TblOff(PSEXTRACTOR,plate_dra),
                                  TblDbl,"plate_ddec",TblOff(PSEXTRACTOR,plate_ddec),
                                  TblDbl,"plate_dist",TblOff(PSEXTRACTOR,plate_dist),
#endif /* MINIMUM_DEBUG */
                                  0,"end",0);

  if (sextractor_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",sextractor_name);
    exit(-1);
  }


  if (verbose) {
    fprintf(stderr,"read %d records for %s\n",sextractor_nrecs,sextractor_name);
  }

  /* Now  divide the catalog file according to spatial bins.  Include stars that may be
     in the margin.  This division is loose so we can pare down the sextractor file (the x and y 
     catalog locations do not agree with the final plate locations */
  for (catalog_index = 0; catalog_index < catalog_nrecs; catalog_index++) {
    pCatalog = &catalog_table[catalog_index];
    pCatalog->pSextractor = NULL;
    pCatalog->spatial_bin = CalculateBin(mosaicWidth,mosaicHeight,pCatalog->x,pCatalog->y,&edgeDist);
    if (pCatalog->magb == 0) {
      pCatalog->spatial_bin = 0;
    }
    catalog_bin_count[pCatalog->spatial_bin]++;
  }

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (catalog_bin_count[spatial_bin] > maxCatalogBinCount) {
      maxCatalogBinCount = catalog_bin_count[spatial_bin];
    }
  }
  maxCatalogBinCount = maxCatalogBinCount * DUPLICATE_FACTOR;
  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    binVector[spatial_bin] = (double *)calloc(maxCatalogBinCount,sizeof(double));
  }




  /* Do the same for the sextractor file, but consider the margins */
  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];
    pSextractor->reference_count = 0;
    pSextractor->spatial_bin = CalculateBin(mosaicWidth,mosaicHeight,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,&edgeDist);
    sextractor_bin_count[pSextractor->spatial_bin]++;
    if (skipMargin == 0) {
      if ((pSextractor->X_IMAGE < leftMarginD) ||
          (pSextractor->X_IMAGE > rightMarginD) ||
          (pSextractor->Y_IMAGE < bottomMarginD) ||
          (pSextractor->Y_IMAGE > topMarginD)) {
        marginCount++;
        pSextractor->spatial_bin = 0;
      }
    }
  }
 
  /* At this point, sort the sextractor table according to spatial bin. 
     We hand code this because the search is so specialized */

  sextractor_index2 = 0;
  pSextractor2 = &sextractor_table[sextractor_index2];
  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    sextractor_index = sextractor_nrecs -1;
    pSextractor = &sextractor_table[sextractor_index];
#ifdef VERIFY_TABLE
    if (verbose) {
      printf("spatial_bin %d, sextractor_index2 %5d\n",spatial_bin,sextractor_index2);
    }
#endif /* VERIFY_TABLE */
    foundReplacement = 1;
    while (foundReplacement) {
#if 0
      if ((spatial_bin == 9) && 
          (sextractor_index == (sextractor_nrecs-1)) &&
          (sextractor_index2 > (sextractor_nrecs-3))) {
        printf("At spatial_bin %d, sextractor_index %d, sextractor_index2 %d\n",spatial_bin,sextractor_index,sextractor_index2);
      }
#endif
      foundReplacement = 0;
      if (pSextractor2->spatial_bin == spatial_bin) {
        sextractor_bin_index[spatial_bin] = sextractor_index2;
        if (sextractor_index2 < (sextractor_nrecs-1)) {
          foundReplacement = 1;
          sextractor_index2++;
          pSextractor2 = &sextractor_table[sextractor_index2];
        } else {
          foundReplacement = 0;
        }
      } else {
        /* Search for a suitable replacement */
        while (sextractor_index > sextractor_index2) {
          if (pSextractor->spatial_bin == spatial_bin) {
            foundReplacement = 1;
            sextractor_bin_index[spatial_bin] = sextractor_index2;
            memcpy(pSextractor3,pSextractor2,sizeof(SEXTRACTOR));
            memcpy(pSextractor2,pSextractor,sizeof(SEXTRACTOR));
            memcpy(pSextractor,pSextractor3,sizeof(SEXTRACTOR));
            break;
          } else {
            sextractor_index--;
            pSextractor = &sextractor_table[sextractor_index];            
          }
        }
        
      }
    }
  }

#ifdef VERIFY_TABLE
  spatial_bin = 0;
  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];
    if (pSextractor->spatial_bin < spatial_bin) {
      fprintf(stderr,"ERROR: table error for bin %d at index %d\n",spatial_bin,pSextractor->spatial_bin,sextractor_index);
      exit(-1);
    }
    if (pSextractor->spatial_bin > spatial_bin) {
      fprintf(stderr,"verify spatial bin %d %d at index %5d %5d\n",spatial_bin,pSextractor->spatial_bin,sextractor_index-1,sextractor_index);
    }
    spatial_bin = pSextractor->spatial_bin;
  }

#endif /* VERIFY_TABLE */

  /* Now we pare down the sextractor catalog size until it matches DUPLICATE_FACTOR * min(the tycho2 size for each bin,maxStars for each bin)
     Here we toss anything in spatial bin 0 */
  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    new_bin_count = DUPLICATE_FACTOR * catalog_bin_count[spatial_bin];
    max_bin_count = DUPLICATE_FACTOR * maxStars * binFractions[spatial_bin];
    if (new_bin_count > max_bin_count) {
      new_bin_count = max_bin_count;
    }
    if (sextractor_bin_count[spatial_bin] < new_bin_count) {
#if 0
      fprintf(stderr,"ERROR spatial bin %d has %d stars, needs %d stars for %s\n",spatial_bin,sextractor_bin_count[spatial_bin],new_bin_count,rootname);
#endif
      new_bin_count = sextractor_bin_count[spatial_bin];
    } else {
      /* Sort everything in the designated spatial bin in descending brightness */
      qsort((void*)&sextractor_table[sextractor_bin_index[spatial_bin-1]+1],sextractor_bin_count[spatial_bin],sizeof(SEXTRACTOR),SextractorImageCompare);
    }
    new_sextractor_bin_count[spatial_bin] = new_bin_count;
    /* Now take the first new_bin_count entries and copy them up */
    if ((sextractor_bin_index[spatial_bin-1]+1) > new_sextractor_nrecs) {
      pSextractor = &sextractor_table[new_sextractor_nrecs];
      pSextractor2 = &sextractor_table[sextractor_bin_index[spatial_bin-1]+1];
      for (sextractor_index = 0; sextractor_index < new_bin_count; sextractor_index++) { 
        memcpy(pSextractor,pSextractor2,sizeof(SEXTRACTOR));
        pSextractor++;
        pSextractor2++;
        new_sextractor_nrecs++;
      }
    } else {
      new_sextractor_nrecs += new_bin_count;
    }

#ifdef VERIFY_TABLE
    if (verbose) {
      printf("spatial_bin %d, new_bin_count %5d, new_sextractor_nrecs %5d \n",spatial_bin,new_bin_count,new_sextractor_nrecs);
    }
#endif /* VERIFY_TABLE */

  }
  sextractor_nrecs = new_sextractor_nrecs;

  /* At this point, sort the catalog table according to spatial bin. 
     We hand code this because the search is so specialized */

  catalog_index2 = 0;
  pCatalog2 = &catalog_table[catalog_index2];
  for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    catalog_index = catalog_nrecs -1;
    pCatalog = &catalog_table[catalog_index];
#ifdef VERIFY_TABLE
    if (verbose) {
      printf("spatial_bin %d, catalog_index2 %5d\n",spatial_bin,catalog_index2);
    }
#endif /* VERIFY_TABLE */
    foundReplacement = 1;
    while (foundReplacement) {
      foundReplacement = 0;
      if (pCatalog2->spatial_bin == spatial_bin) {
        catalog_bin_index[spatial_bin] = catalog_index2;
        if (catalog_index2 < (catalog_nrecs-1)) {
          foundReplacement = 1;
          catalog_index2++;
          pCatalog2 = &catalog_table[catalog_index2];
        } else {
          foundReplacement = 0;
        }
      } else {
        /* Search for a suitable replacement */
        while (catalog_index > catalog_index2) {
          if (pCatalog->spatial_bin == spatial_bin) {
            foundReplacement = 1;
            catalog_bin_index[spatial_bin] = catalog_index2;
            memcpy(pCatalog3,pCatalog2,sizeof(CATALOG));
            memcpy(pCatalog2,pCatalog,sizeof(CATALOG));
            memcpy(pCatalog,pCatalog3,sizeof(CATALOG));
            break;
          } else {
            catalog_index--;
            pCatalog = &catalog_table[catalog_index];            
          }
        }
        
      }
    }
  }

#ifdef VERIFY_TABLE
  spatial_bin = 0;
  for (catalog_index = 0; catalog_index < catalog_nrecs; catalog_index++) {
    pCatalog = &catalog_table[catalog_index];
    if (pCatalog->spatial_bin < spatial_bin) {
      fprintf(stderr,"ERROR: table error for bin %d at index %d\n",spatial_bin,pCatalog->spatial_bin,catalog_index);
      exit(-1);
    }
    if (pCatalog->spatial_bin > spatial_bin) {
      fprintf(stderr,"verify spatial bin %d %d at index %5d %5d\n",spatial_bin,pCatalog->spatial_bin,catalog_index-1,catalog_index);
    }
    spatial_bin = pCatalog->spatial_bin;
  }

#endif /* VERIFY_TABLE */

  /* Now we pare down the catalog catalog size until it matches the sextractor stars for each bin dividided by the duplicate factor
     Here we toss anything in spatial bin 0 */
  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    new_bin_count = (1.0 * new_sextractor_bin_count[spatial_bin])/DUPLICATE_FACTOR;
    if (catalog_bin_count[spatial_bin] < new_bin_count) {
      new_bin_count = catalog_bin_count[spatial_bin];
    } else {
      /* Sort everything in the designated spatial bin in descending brightness */
      qsort((void*)&catalog_table[catalog_bin_index[spatial_bin-1]+1],catalog_bin_count[spatial_bin],sizeof(CATALOG),CatalogImageCompare);
    }
    new_catalog_bin_count[spatial_bin] = new_bin_count;
    /* Now take the first new_bin_count entries and copy them up */
    if ((catalog_bin_index[spatial_bin-1]+1) > new_catalog_nrecs) {
      pCatalog = &catalog_table[new_catalog_nrecs];
      pCatalog2 = &catalog_table[catalog_bin_index[spatial_bin-1]+1];
      for (catalog_index = 0; catalog_index < new_bin_count; catalog_index++) { 
        memcpy(pCatalog,pCatalog2,sizeof(CATALOG));
        pCatalog++;
        pCatalog2++;
        new_catalog_nrecs++;
      }
    } else {
      new_catalog_nrecs += new_bin_count;
    }

#ifdef VERIFY_TABLE
    if (verbose) {
      printf("spatial_bin %d, new_bin_count %5d, new_catalog_nrecs %5d \n",spatial_bin,new_bin_count,new_catalog_nrecs);
    }
#endif /* VERIFY_TABLE */

  }
  catalog_nrecs = new_catalog_nrecs;

  


  /* We are now ready to begin our match search.  Sort both tables by the height */
  qsort((void*)sextractor_table,sextractor_nrecs,sizeof(SEXTRACTOR),SextractorYCompare);
  qsort((void*)catalog_table,catalog_nrecs,sizeof(CATALOG),CatalogYCompare);

  sextractor_min_index = 0;
  sextractor_index = 0;
  time(&curTime);
  curTime -= startTime;
  if (verbose) {
    fprintf(stderr,"Begin search for matches at %d seconds\n",curTime);
  }
  for (catalog_index = 0; catalog_index < catalog_nrecs; catalog_index++) {
    pCatalog = &catalog_table[catalog_index];
#if 0
    if ((verbose != 0)  && 
        ((catalog_index % 1000) == 0)) {
      time(&curTime);
      curTime -= startTime;
      printf("Searching catalog_index %d, sextractor_min_index %d, sextractor_index %d, trialCount %d at %d seconds\n",catalog_index,sextractor_min_index,sextractor_index,trialCount,curTime);
    }
#endif 
    sextractor_index = sextractor_min_index;
    while (sextractor_index < sextractor_nrecs) {
      pSextractor = &sextractor_table[sextractor_index];
      yDist = pSextractor->Y_IMAGE - pCatalog->y;
      if (yDist < -pixelRadius) {
        /* Never going to look at this record again */
        sextractor_min_index++;
      }
      if (yDist > pixelRadius) {
        /* Stop searching for matches for this object */
        break;
      }
      if (yDist < 0.0) {
        yDist = -yDist;
      }
      xDist = pSextractor->X_IMAGE - pCatalog->x;
      if (xDist < 0.0) {
        xDist = - xDist;
      }
      if ((xDist < pixelRadius) &&
          (yDist < pixelRadius)) {
        trialCount++;
        /* Potential match. Calculate total distance */
        totalDist = sqrt(sqr(xDist)+sqr(yDist));
        if (totalDist < pixelRadius) {
          /* We have a match. */
          if (pCatalog->pSextractor != NULL ) {
            /* Already matched something else, choose the
               brighter object */
            if (pSextractor->MAG_ISO < pCatalog->pSextractor->MAG_ISO) {
              pCatalog->pSextractor->reference_count--;
              pCatalog->pSextractor = pSextractor;
              pCatalog->pSextractor->reference_count++;
              doubleMatchCount++;
              double_match_bin[pCatalog->spatial_bin]++;
            }
          } else {
            pCatalog->pSextractor = pSextractor;
            pSextractor->reference_count++;
          }
        }


      }
      sextractor_index++;

    }
  }
  time(&curTime);
  curTime -= startTime;
  if (verbose) {
    fprintf(stderr,"Begin search for duplicates at %d seconds\n",curTime);
  }

  /* At this point, search the sextractor table for double references and keep the one with the brightest catalog entry */
  for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
    pSextractor = &sextractor_table[sextractor_index];
    if (pSextractor->reference_count > 1) {
      doubleReferenceCount++;
      double_reference_bin[pSextractor->spatial_bin]++;
      pCatalog2 = NULL; /* Our desired catalog object */
      for (catalog_index = 0; catalog_index < catalog_nrecs; catalog_index++) {
        pCatalog = &catalog_table[catalog_index];
        if (pCatalog->pSextractor == pSextractor) {
          if (pCatalog2 == NULL) {
            pCatalog2 = pCatalog;
          } else {
            if (pCatalog->magb < pCatalog2->magb) {
              /* New one is brighter, get rid of the old one */
              pCatalog2->pSextractor->reference_count--;
              pCatalog2->pSextractor = NULL;
              pCatalog2 = pCatalog;
            } else {
              /* Old one is brighter.  Get rid of the new one */
              pCatalog->pSextractor->reference_count--;
              pCatalog->pSextractor = NULL;

            }

            if (pSextractor->reference_count <= 1) {
              /* Found them all, exit now */
              break;
            }
          }

        }
      
      
      }
      if ((pCatalog2 == NULL) || (pSextractor->reference_count != 1)) {
        fprintf(stderr,"ERROR: No duplicate found or reference count bad for index %d, pointer %08x, reference count %d for %s\n",sextractor_index,
                pCatalog2,pSextractor->reference_count,rootname);
        exit(-1);
      
      }

    }

  }

  if (verbose) {
    fprintf(stderr,"Writing Results\n");
  }

  /* We can now write everything out */
  for (catalog_index = 0; catalog_index < catalog_nrecs; catalog_index++) {
    pCatalog = &catalog_table[catalog_index];
    pSextractor = pCatalog->pSextractor;
    if (pSextractor != NULL) {
      outCount++;

      factor = cos(DEGREES_TO_RAD*((pSextractor->dec + pCatalog->dec)/2.0));
      dra = 3600*factor*(pSextractor->ra - pCatalog->ra);
      ddec = 3600*(pSextractor->dec - pCatalog->dec);
      drad = sqrt((dra*dra) + (ddec*ddec));
      spatial_bin = pSextractor->spatial_bin;
      binVector[spatial_bin][bin_output_count[spatial_bin]] = drad;
#if 0
      printf("binVector %d %5.2f\n",spatial_bin,drad);
#endif
      bin_output_count[spatial_bin]++;
      if (bin_output_count[spatial_bin] >= maxCatalogBinCount) {
        fprintf(stderr,"ERROR: bin_output_count %d for bin %d is greater than %d\n",bin_output_count[spatial_bin],spatial_bin,maxCatalogBinCount);
        exit(-1);
      }


      fprintf(outHandle,"%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%s\t%.6f\t%.6f\t%.3f\t%.2f\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
              pSextractor->NUMBER,
              pSextractor->MAG_ISO,
              pSextractor->X_IMAGE,
              pSextractor->Y_IMAGE,
              pSextractor->ra,
              pSextractor->dec,
              pSextractor->FWHM_IMAGE,
              pSextractor->FWHM_WORLD,
              pSextractor->AFLAGS,
              pSextractor->BFLAGS,
              "TYCHO", /* pCatalog->REF, */
              pCatalog->ra,
              pCatalog->dec,
              pCatalog->magb,
              0.0, /* pCatalog->color, */
              0,   /* pCatalog->VFlag, */
              0,   /* pCatalog->MAGFlag, */
              0,   /* pCatalog->class, */
              dra,
              ddec,
              drad,
              pSextractor->THRESHOLD,
              pSextractor->FLUX_MAX);

    }
  }


#if 1
  {
    FILE *tmpHandle;
    tmpHandle = fopen("los.reg","wt");
    for (catalog_index = 0; catalog_index < catalog_nrecs; catalog_index++) {
      pCatalog = &catalog_table[catalog_index];
      pSextractor = pCatalog->pSextractor;
      if (pSextractor != NULL) {
        outCount++;
        fprintf(tmpHandle,"j2000;box(%f,%f,0.02,0.02,0.0);\n",pSextractor->ra,pSextractor->dec);
        fprintf(tmpHandle,"j2000;box(%f,%f,0.02,0.02,45.0);\n",pCatalog->ra,pCatalog->dec);
      }
    }
    fclose(tmpHandle);
  }

#endif


  time(&curTime);
  curTime -= startTime;
  if (verbose) {
    for (spatial_bin = 0; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      binMedian = 0.0;
      binRms = 0.0;
      if (spatial_bin > 0) {
        CalcMedianAndRMS(bin_output_count[spatial_bin],0,binVector[spatial_bin],&binMedian,&binRms,1,3.0,0);
      }
      printf("matchtycho2: bin %d, catalog %6d new %4d, image %6d, new %6d dmatch %6d dref %6d output %6d median %5.1f drad rms %5.1f for %s\n",spatial_bin,
             catalog_bin_count[spatial_bin],
             new_catalog_bin_count[spatial_bin],
             sextractor_bin_count[spatial_bin],
             new_sextractor_bin_count[spatial_bin],
             double_match_bin[spatial_bin],
             double_reference_bin[spatial_bin],
             bin_output_count[spatial_bin],
             binMedian/arcsecPerPixel,
             binRms/arcsecPerPixel,
             rootname);
    }
  


  }

  printf("matchtycho2: %d seconds margin %d double match %d double reference %d trialCount %d output %d cat %d tycho2 %d pixelRadius %d for %s\n",
         curTime,
         marginCount,
         doubleMatchCount,
         doubleReferenceCount,
         trialCount,
         outCount,
         sextractor_nrecs,
         catalog_nrecs,
         pixelRadius,
         rootname);

  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (binVector[spatial_bin] != NULL) {
      free(binVector[spatial_bin]);
    }
  }


  if (catalog_handle != NULL) {
    Close(catalog_handle);
    catalog_handle = NULL;
  }
  if (catalog_table != NULL) {
    Free(catalog_table);
  }
  if (catalog_header != NULL) {
    table_hdrfree(catalog_header);
  }

  if (sextractor_handle != NULL) {
    Close(sextractor_handle);
    sextractor_handle = NULL;
  }
  if (sextractor_table != NULL) {
    Free(sextractor_table);
  }
  if (sextractor_header != NULL) {
    table_hdrfree(sextractor_header);
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }


  return(EXIT_SUCCESS);
}
