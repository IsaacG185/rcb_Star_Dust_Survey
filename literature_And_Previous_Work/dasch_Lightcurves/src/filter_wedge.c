// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* filter_wedge.c
 *
 *  Search for Pickering Wedge images in a limited number of series
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm -L/usr/lib64/mysql -L/usr/lib/mysql -lmysqlclient filter_wedge.c pipelineutils.a -ltable -lutil -lwcs -o filter_wedge 
 *
 *  filter_wedge -p -v -s 1.786 -w 17412 -h 22026 -r i31013_00_01r180ww -i /dasch/Pipeline/match/i31013_00_01r180ww_tnx.tmp -o /dasch/Pipeline/match/i31013_00_01r180ww_wedge.db
 *
 *  filter_wedge -p -v  -s 1.787 -w 17412 -h 22026 -r i31090_00_01r180ww -i /dasch/Pipeline/match/i31090_00_01r180ww_tnx.tmp -o /dasch/Pipeline/match/i31090_00_01r180ww_wedge.db
 *
 * filter_wedge -v  -p  -s 1.023 -w 17412 -h 22026 -r mc00380_03_01ww -i /dasch/Pipeline/match/mc00380_03_01ww_tnx.tmp -o /dasch/Pipeline/match/mc00380_03_01ww_wedge.db
 *
 *
 * May  6, 2008 Edward J. Los - Inital version
 * May 29, 2008 Edward J. Los - add B plate result
 * May 31, 2008 Edward J. Los - Switch from median to mean + RMS
 * Jan 27, 2008 Edward J. Los - split FLAGS into AFLAGS and BFLAGS
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 * Dec 28, 2009 Edward J. Los - Set or clear the Wedge bit in the FitWCS column of the mosaics table
 *                              Set or clear the wedge bit in the quality column of the plates table
 * Feb  4, 2010 Edward J. Los - look for dsy ghost images
 * Mar  7, 2010 Edward J. Los - add second test for "i" series wedge images using i23336 as a prototype with 276" separation.
 * Feb  6, 2012 Edward J. Los - Add two suspected mc series at 130" and 147"
 * Jul 23, 2012 Edward J. Los - Add high background filter support
 * Aug  4, 2012 Edward J. Los - Correct Y_IMAGE sort position.
 * Jun  2, 2014 Edward J. Los - Add new dsy distance found by the search_close algorithm.
 * Apr  8, 2015 Edward J. Los - Revise the B plate parameters
 * Feb  4, 2017 Edward J. Los - Make failure to read an empty background table a warning.
 * Apr 25, 2017 Edward J. Los   Move WEDGEENTRY from filter_wedge for sharing between filter_wedge.c and update_quality.c
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"

#define MAX_BUFFER 512
#define SEARCH_TABLE_LENGTH 300
#define BORDER_EDGE 0.025
#define NUM_GMT_FILES 7
#define MATCH_ARRAY_ALLOC 1000000
#define SQUARE_BIN_OBJECTS 10
#define CLIP_RMS_FACTOR 1.0
#define MAG_ISO_CRIT_MARGIN 1.0
#define EXPECTED_PEAKS 1
#define EXCESS_LIMIT 10.0






typedef struct _starimage {
  int NUMBER;         /* Sextractor reference number */
  int AFLAGS;          /* flags word */
  int MATCH_NUMBER;   /* If nonzero, matching bright star */
  double MAG_ISO;
  double FWHM_WORLD;
  double X_IMAGE;
  double Y_IMAGE;
  
} STARIMAGE,*PSTARIMAGE;

typedef struct _searchentry {
  int sextractorIndex;          /* Index into the sextractor table */
  double MAG_ISO;     /* Brightness of this object */
  double X_IMAGE;
  double Y_IMAGE;
} SEARCHENTRY,*PSEARCHENTRY;

SEARCHENTRY searchTable[SEARCH_TABLE_LENGTH];

typedef struct _matchentry {
  double xValue;
  double yValue;
  double magsum;  /* dim + bright iso magnitude */
  double magdiff; /* dim - bright iso magnitude */
  double drad;    /* Distance from centroid */
  int brightIndex;
  int dimIndex;
  int selected;
  int binNumber;
} MATCHENTRY,*PMATCHENTRY;

typedef struct _magsum_index {
  double magsum;
  int index;
} MAGSUM_INDEX,*PMAGSUM_INDEX;


/* Given an angle and a radius, compute the limits */
void GetLimits(double binAngle,double radius,double *pMinX,double *pMaxX,double *pMinY,double *pMaxY)
{
  double xval = radius * cos(binAngle*DEGREES_TO_RAD);
  double yval = radius * sin(binAngle*DEGREES_TO_RAD);
  if (xval < *pMinX) {
    *pMinX = xval;
  }
  if (xval > *pMaxX) {
    *pMaxX = xval;
  }
  if (yval < *pMinY) {
    *pMinY = yval;
  }
  if (yval > *pMaxY) {
    *pMaxY = yval;
  }
#if 0
  printf("binAngle %f radius %f, xval %f yval %f\n",binAngle,radius,xval,yval);
#endif
  return;
}

/* NOTE: The following is a reverse sort */
int magsumCompare(const void *first, const void *second) 
{
  double magsumFirst = ((PMAGSUM_INDEX)first)->magsum;
  double magsumSecond = ((PMAGSUM_INDEX)second)->magsum;
  if (magsumFirst < magsumSecond) {
    return(1);
  } else if (magsumFirst > magsumSecond) {
    return(-1);
  } else {
    return(0);
  }

}

int sextractorCompare(const void *first, const void *second) 
{
  int NUMBERFirst = ((PSTARIMAGE)first)->NUMBER;
  int NUMBERSecond = ((PSTARIMAGE)second)->NUMBER;
  if (NUMBERFirst < NUMBERSecond) {
    return(-1);
  } else if (NUMBERFirst > NUMBERSecond) {
    return(1);
  } else {
    return(0);
  }

}

int sextractorYCompare(const void *first, const void *second) 
{
  double Y_IMAGEFirst = ((PSTARIMAGE)first)->Y_IMAGE;
  double Y_IMAGESecond = ((PSTARIMAGE)second)->Y_IMAGE;
  if (Y_IMAGEFirst < Y_IMAGESecond) {
    return(-1);
  } else if (Y_IMAGEFirst > Y_IMAGESecond) {
    return(1);
  } else {
    return(0);
  }

}

int backgroundCompare(const void *first, const void *second) 
{
  int NUMBERFirst = ((PHIGHBACKGROUND)first)->NUMBER;
  int NUMBERSecond = ((PHIGHBACKGROUND)second)->NUMBER;
  if (NUMBERFirst < NUMBERSecond) {
    return(-1);
  } else if (NUMBERFirst > NUMBERSecond) {
    return(1);
  } else {
    return(0);
  }

}


int main(int argc,char *argv[])
{
  char *argstr;
  char cmdchar;
  int nvals;
  int errorFlag = 0;
  char fileroot[MAX_BUFFER];
  char outfile[MAX_BUFFER];
  FILE *outHandle = NULL;
  FILE *gmtHandle[NUM_GMT_FILES];
  char gmtname[MAX_BUFFER];
  char tmpname[MAX_BUFFER];
  int gmtIndex;
  char *slashPtr;
  char *curPtr;
  int index;
  int searchIndex;
  int insertIndex;
  int seriesIndex = -1;
  PWEDGEENTRY pWedgeEntry = NULL;
  PSEARCHENTRY pSearchEntry1 = NULL;
  PSEARCHENTRY pSearchEntry2 = NULL;
  int searchTableSize = 0;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int candidateCount = 0;


  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  PSTARIMAGE sextractor_table = NULL;
  size_t sextractor_nrecs = 0;
  int sextractorIndex;
  int sextractor_index;
  int baseIndex;
  PSTARIMAGE pSextractor1 = NULL;
  PSTARIMAGE pSextractor2 = NULL;
  PSTARIMAGE pSextractorBright = NULL;
  PSTARIMAGE pSextractorDim = NULL;
 

  double minX;
  double maxX;
  double minY;
  double maxY;
  double binAngleSize;
  int angleBins;
  int *angleBinTable = NULL;
  int angleBinNumber;
  int averageBin;
  double sigmaSum = 0;
  double sigmaSqr = 0;
  double binMean;
  double binStd;
  double binValue;
  double binAngle;

  int verbose = 0;
  int doPlots = 0;
  time_t startTime;
  time_t curTime;

  char series[MAX_SERIES_STRING];
  int plateNumber;
  double excess;
  double maxExcess = -1000000.0;
  int maxIndex1 = 0;
  int maxIndex2 = 0;
  int index1;
  int index2;

  double minDec;
  double maxDec;
  double minRa;
  double maxRa;
  double drad;
  double xValue;
  double yValue;
  double angle;
  PMATCHENTRY pMatchArray = NULL;
  PMATCHENTRY pSelectedMatchArray = NULL;
  PMATCHENTRY pNewMatchArray;
  int matchArrayAlloc = 0;
  int matchArraySize = 0;
  int selectedCount = 0;

  PMATCHENTRY pMatchEntry;
  PMATCHENTRY pSelectedMatchEntry;

  int totSqrBins; /* Total number of square bins */
  int numSqrBins; /* Number of sqare bins per axis */
  int xSqrBin;
  int ySqrBin;
  int iSqrBin;
  int matchIndex;
  int binIndex;
  int maxSqrBinValue;
  int maxSqrBinIndex;
  int aveBinValue;
  int startIndex;
  int endIndex;
  double brightestMagIso;
  int brightestMagIsoIndex;
  PMAGSUM_INDEX pMagsumIndexList = NULL;
  PMAGSUM_INDEX pMagsumIndex;
  int plotCount = 0;
  double arcsecPerPixel;
  double degreesPerPixel;
 
  double * vector = NULL;
  int curClipCount;
  double MAG_ISO_med;
  double MAG_ISO_rms;
  double MAG_ISO_crit;
  int numPeaks = 0;
  int abovePeak;
  int initAbovePeak;
  int finalAbovePeak;
  int solutionNumber = 0;
  int readFlag = 0;
  int FitWCS;
  int quality;
  int wedgeFlag = -1;

  File background_handle = NULL;
  char background_name[MAX_BUFFER];
  TableHead background_header = NULL;
  PHIGHBACKGROUND background_table = NULL;
  PHIGHBACKGROUND pBackground;
  PHIGHBACKGROUND pBackgroundPrev;
  size_t background_nrecs;
  int background_index;

  time(&startTime);



  background_name[0] = 0;

  sextractor_name[0] = 0;
  fileroot[0] = 0;
  arcsecPerPixel = 0.0;
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

        case 'C': /* background file name */
        case 'c':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(background_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'r': /* file root */
        case 'R':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(fileroot,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
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

        case 'e': /* solution number */
        case 'E':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&solutionNumber);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the solutionNumber %s\n",*argv);
              errorFlag = 1;
            }
          }
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

 
        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'p': /* do plotting */
        case 'P':
          doPlots = 1;
          break;





        default:
          fprintf(stderr,"ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }

  if (mosaicWidth == 0) {
    fprintf(stderr,"ERROR: No mosaic width specified\n");
    errorFlag = 1;
  }
  if (mosaicHeight == 0) {
    fprintf(stderr,"ERROR: No mosaic height specified\n");
    errorFlag = 1;
  }
  if (arcsecPerPixel == 0.0) {
    fprintf(stderr,"ERROR: No arcsec per pixel plate scale was specified\n");
    errorFlag = 1;
  }
  degreesPerPixel = arcsecPerPixel/3600.0;

  if (outfile[0] == 0) {
    fprintf(stderr,"ERROR: No output filename was specified\n");
    errorFlag = 1;
  }
  if (sextractor_name[0] == 0) {
    fprintf(stderr,"ERROR: No sextractor filename was specified\n");
    errorFlag = 1;
  }
  if (fileroot[0] == 0) {
    fprintf(stderr,"ERROR: No file root was specified\n");
    errorFlag = 1;
  }
  if (ParseFilename2(sextractor_name,series,&plateNumber) == 0) {
    fprintf(stderr,"ERROR: Failed to parse filename %s\n",sextractor_name);
    errorFlag = 1;
  }

  for (index = 0; index < wedgeTableSize; index++) {
    pWedgeEntry = &wedgeTable[index];
    if (strlen(pWedgeEntry->series) == 0) {
      break;
    }
    if (strcmp(pWedgeEntry->series,series) == 0) {
      break;
    }
  }
  if (strlen(pWedgeEntry->series) == 0) {
    fprintf(stderr,"ERROR: Series %s is not supported by this program\n",series);
    errorFlag = 1;
  }

  

  sextractor_handle = Open(sextractor_name,"r");
  if (sextractor_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the sextractor file %s\n",sextractor_name);
  } else {
    if (verbose) {
      fprintf(stderr,"Found sextractor file %s\n",sextractor_name);
    }
  }

   
  if (doPlots) {
    for (gmtIndex = 0; gmtIndex < NUM_GMT_FILES; gmtIndex++) {
      
      strcpy(gmtname,outfile);
      curPtr = gmtname;
      while ((slashPtr = strstr(curPtr,"/")) != NULL) {
        curPtr = slashPtr+1;
      }
      if (curPtr != NULL) {
        *curPtr = 0;
      }
      strcat(gmtname,fileroot);
      sprintf(tmpname,"_gmt%d.txt",gmtIndex);
      strcat(gmtname,tmpname);
      gmtHandle[gmtIndex] = fopen(gmtname,"wt");
      if (gmtHandle[gmtIndex] == NULL) {
        errorFlag = 1;
        fprintf(stderr,"ERROR: Failed to open the gmt file %s\n",gmtname);
      } else {
        if (verbose) {
          fprintf(stderr,"GMT file %s\n",gmtname);
        }
      }
    }
  }

  if (background_name[0] != 0) {
    background_handle = Open(background_name,"r");
    if (background_handle == NULL) {
      fprintf(stderr,"No background file found  %s\n",background_name);
    } else {
      if (verbose) {
        fprintf(stderr,"Found plate background file %s\n",background_name);
      }
    }
  }


  if (errorFlag) {
    fprintf(stderr,"Usage: filter_wedge -r <file root> \n");
    fprintf(stderr,"                      -i <sextractor file> \n");
    fprintf(stderr,"                      -o  <output file>\n");
    fprintf(stderr,"                      -p  write GMT plotting files\n");
    fprintf(stderr,"                      -s <arcsec per pixel>\n");
    fprintf(stderr,"                      -e <solutionNumber>\n");
    fprintf(stderr,"                      -c <background file>\n");
    fprintf(stderr,"                      -v  verbose\n");

    return(-1);
  }
  /* Now read in the sextractor file */

  sextractor_header = table_header(sextractor_handle,TABLE_PARSE);
  if (sextractor_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",sextractor_name);
    return(-1);
  }

  sextractor_table = table_loadva(sextractor_handle,
                                  &sextractor_header,
                                  NULL, /* hbase */
                                  NULL, /* rows */
                                  NULL,
                                  sizeof(STARIMAGE),
                                  &sextractor_nrecs,
                                  TblInt,"NUMBER",TblOff(PSTARIMAGE,NUMBER),
                                  TblInt,"AFLAGS",TblOff(PSTARIMAGE,AFLAGS),
                                  TblDbl,"MAG_ISO",TblOff(PSTARIMAGE,MAG_ISO),
                                  TblDbl,"FWHM_WORLD",TblOff(PSTARIMAGE,FWHM_WORLD),
                                  TblDbl,"X_IMAGE",TblOff(PSTARIMAGE,X_IMAGE),
                                  TblDbl,"Y_IMAGE",TblOff(PSTARIMAGE,Y_IMAGE),
                                  0,"end",0);
  if (sextractor_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",sextractor_name);
    return(-1);
  }
  if (verbose) {
    fprintf(stderr,"read %d records for %s\n",sextractor_nrecs,sextractor_name);
  }
  /* Sort this table in increasing NUMBER */
  qsort(sextractor_table,sextractor_nrecs,sizeof(STARIMAGE),sextractorCompare);
  /* Now read in the background file */
  if (background_handle != NULL) {
    background_header = table_header(background_handle,TABLE_PARSE);
    if (background_header == NULL) {
      fprintf(stderr,"ERROR: Failed to read header for %s\n",background_name);
      return(-1);
    }

    background_table = table_loadva(background_handle,
                                    &background_header,
                                    NULL, /* hbase */
                                    NULL, /* rows */
                                    NULL,
                                    sizeof(HIGHBACKGROUND),
                                    &background_nrecs,
                                    TblInt,"NUMBER" ,TblOff(PHIGHBACKGROUND,NUMBER),
                                    TblInt,"AFLAGS"  ,TblOff(PHIGHBACKGROUND,AFLAGS),
                                    0,"end",0);
    if (background_table == NULL) {
      /* This is o.k. if there are no entries in the file */
      fprintf(stderr,"WARNING: Failed to read table for %s\n",background_name);
    } else {
      if (verbose) {
        fprintf(stderr,"read %d records for %s\n",background_nrecs,background_name);
      }
      /* Sort this table in increasing NUMBER */
      qsort(background_table,background_nrecs,sizeof(HIGHBACKGROUND),backgroundCompare);
      /* Now we step through the sextractor table and set the FILTER_AFLAG_BACKGROUND bit for any objects with high background */
      background_index = 0;
      pBackground = &background_table[background_index];
      for (sextractor_index = 0; sextractor_index < sextractor_nrecs; sextractor_index++) {
        pSextractor1 = &sextractor_table[sextractor_index];
        /* Be sure the background flag is clear */
        pSextractor1->AFLAGS &= ~(1<<FILTER_AFLAG_BACKGROUND) ;
        while (pBackground->NUMBER < pSextractor1->NUMBER) {
          background_index++;
          if (background_index >= background_nrecs) {
            break;
          }
          pBackground = &background_table[background_index];
        }
 
        if (pBackground->NUMBER == pSextractor1->NUMBER) {
          if ((pBackground->AFLAGS & (1<<FILTER_AFLAG_BACKGROUND)) != 0) {
            /* Set our background flag */
            pSextractor1->AFLAGS |= (1<<FILTER_AFLAG_BACKGROUND) ;
          }
        }

      }

  
  



      /* We are done with the background table.  Clean up memory */

      if (background_table != NULL) {
        free(background_table);
        background_table = NULL;
      }
    }
    if (background_header != NULL) {
      table_hdrfree(background_header);
      background_header = NULL;
    }
    if (background_handle != NULL) {
      Close(background_handle);
      background_handle = NULL;
    }
    pBackground = NULL;
    pBackgroundPrev = NULL;
  }


  /* Sort this table in increasing Y_IMAGE */
  qsort(sextractor_table,sextractor_nrecs,sizeof(STARIMAGE),sextractorYCompare);


  matchArrayAlloc = MATCH_ARRAY_ALLOC;
  pMatchArray = (PMATCHENTRY)calloc(matchArrayAlloc,sizeof(MATCHENTRY));
  if (pMatchArray == NULL) {
    fprintf(stderr,"ERROR Failed to allocate pMatchArray\n");
    exit(-1);
  }
  seriesIndex = -1;
  fprintf(stderr,"filter_wedge of %s %s \n Input Filename %s\n Output Filename %s\n Series %s separation %f drad %f mosaic width %d mosaic height %d scale %f arcsec/pixel\n",
          __DATE__,__TIME__,sextractor_name,outfile,
          pWedgeEntry->series,
          pWedgeEntry->separation,
          pWedgeEntry->meanplussigma_drad,
          mosaicWidth,
          mosaicHeight,
          arcsecPerPixel);


  while (1) {
    seriesIndex++;
    pWedgeEntry  = &wedgeTable[seriesIndex];
    if (strlen(pWedgeEntry->series) == 0) {
      break;
    }
    if (strcmp(pWedgeEntry->series,series) != 0) {
      continue;
    }
    candidateCount = 0;
    matchArraySize = 0;
    searchTableSize = 0;
    selectedCount = 0;
    maxIndex1 = 0;
    maxIndex2 = 0;
    maxExcess = -1000000.0;
    numPeaks = 0;
    sigmaSum = 0;
    sigmaSqr = 0;
    
    memset(searchTable,0,sizeof(searchTable));
    for (searchIndex = 0; searchIndex < SEARCH_TABLE_LENGTH; searchIndex++) {
      pSearchEntry1 = &searchTable[searchIndex];
      pSearchEntry1->sextractorIndex = -1;
    }



    /* Convert arcsec to degrees */
    pWedgeEntry->separation = pWedgeEntry->separation/3600.;
    pWedgeEntry->meanplussigma_drad = pWedgeEntry->meanplussigma_drad/3600.;
  
    binAngleSize = atan(DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad/pWedgeEntry->separation)*RAD_TO_DEGREES;
    angleBins = 360.0/binAngleSize;
    angleBins = angleBins * 2; /* Double the angle bins (this number will always be even) to that
                                  we chan check across boundaries */
    fprintf(stderr,"Angle bin size %f degrees, Number of angle bins %d for %s\n",binAngleSize,angleBins/2,fileroot);
    if (angleBinTable != NULL) {
      free(angleBinTable);
      angleBinTable = NULL;
    }
    angleBinTable = (int *)calloc(angleBins,sizeof(int));

    if (angleBins < 2) {
      fprintf(stderr,"ERROR: angleBins %d is less than 2\n",angleBins);
      exit(-1);
    }

    minX = mosaicWidth * BORDER_EDGE;
    minY = mosaicHeight * BORDER_EDGE;
    maxX = mosaicWidth - minX;
    maxY = mosaicHeight - minY;
    fprintf(stderr,"Min, Max X: %f %f, Min Max Y %f %f\n",minX,maxX,minY,maxY);


 
    /* We assume that most of the objects are dim.  Calculate the clipped median to
       get the limiting dim MAG_ISO */
  
    if (vector != NULL) {
      free(vector);
      vector = NULL;
    }
    vector = (double *)calloc(sextractor_nrecs,sizeof(double));
    for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
      pSextractorBright = &sextractor_table[sextractorIndex];
      vector[sextractorIndex] = pSextractorBright->MAG_ISO;
    }

    curClipCount = sextractor_nrecs;
    curClipCount = CalcMedianAndRMS(curClipCount,10,vector,&MAG_ISO_med,&MAG_ISO_rms,1,CLIP_RMS_FACTOR,0);
    if (verbose) {
      printf("curClipCount %d MAG_ISO median %f, MAG_ISO rms %f\n",curClipCount,MAG_ISO_med,MAG_ISO_rms);
    }


    /* At this point, we will produce a sorted list of the brightest objects on the plate
       with a reasonable FWHM */

    for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
      pSextractorBright = &sextractor_table[sextractorIndex];
#if 0
      fprintf(stderr,"Object %6d MAG_ISO %10f FWHM_WORLD %10f\n",
              sextractorIndex,
              pSextractorBright->MAG_ISO,
              pSextractorBright->FWHM_WORLD);
#endif
      if ((pSextractorBright->FWHM_WORLD / 2.0) > pWedgeEntry->separation) {
#if 0
        /* This object is so large it will cover any Pickering Wedge object */
        fprintf(stderr,"LARGE FWHM_WORLD: %6d MAG_ISO %10f FWHM_WORLD %10f deg %10f arcsec > %10f  arcsec\n",
                sextractorIndex,
                pSextractorBright->MAG_ISO,
                pSextractorBright->FWHM_WORLD,
                pSextractorBright->FWHM_WORLD* 3600.0/2.0,
                pWedgeEntry->separation*3600.0
                );
#endif
        continue;
      }
      if ((pSextractorBright->X_IMAGE < minX) ||
          (pSextractorBright->X_IMAGE > maxX) ||
          (pSextractorBright->Y_IMAGE < minY) ||
          (pSextractorBright->Y_IMAGE > maxY)) {
#if 0
        fprintf(stderr,"IMAGE ON EDGE %10f %10f\n",
                pSextractorBright->X_IMAGE,
                pSextractorBright->Y_IMAGE);
#endif
        continue;
      }

      /* Insert this object in the search table */
      if ((searchTableSize < SEARCH_TABLE_LENGTH) ||
          (pSextractorBright->MAG_ISO < searchTable[SEARCH_TABLE_LENGTH-1].MAG_ISO)) {
        for (searchIndex = 0; searchIndex < SEARCH_TABLE_LENGTH; searchIndex++) {
          pSearchEntry1 = &searchTable[searchIndex];
          if (pSearchEntry1->sextractorIndex < 0) {
            pSearchEntry1->sextractorIndex = sextractorIndex;
            pSearchEntry1->MAG_ISO = pSextractorBright->MAG_ISO;
            pSearchEntry1->X_IMAGE     = pSextractorBright->X_IMAGE;
            pSearchEntry1->Y_IMAGE     = pSextractorBright->Y_IMAGE;
            searchTableSize++;
            break;
          } else {
            if (pSextractorBright->MAG_ISO < pSearchEntry1->MAG_ISO) {
              /* Our entry goes here - move everything else down */
              for (insertIndex = (SEARCH_TABLE_LENGTH-1); insertIndex > searchIndex; insertIndex--) {
                memcpy(&searchTable[insertIndex],&searchTable[insertIndex-1],sizeof(SEARCHENTRY));
              }
              pSearchEntry1->sextractorIndex = sextractorIndex;
              pSearchEntry1->MAG_ISO = pSextractorBright->MAG_ISO;
              pSearchEntry1->X_IMAGE     = pSextractorBright->X_IMAGE;
              pSearchEntry1->Y_IMAGE     = pSextractorBright->Y_IMAGE;
              if (searchTableSize < SEARCH_TABLE_LENGTH) {
                searchTableSize++;
              }
              break;
            }
          }
        } /* Insertion loop */
#if 0
        for (insertIndex = 0; insertIndex < SEARCH_TABLE_LENGTH; insertIndex++) {
          pSearchEntry2 = &searchTable[insertIndex];
          if (pSearchEntry2->sextractorIndex >= 0) {
            pSextractor2 = &sextractor_table[pSearchEntry2->sextractorIndex];
            fprintf(stderr,"Entry %3d  Index %7d MAG_ISO %10f NUMBER %6d X_IMAGE %10f Y_IMAGE %10f\n",
                    insertIndex,
                    pSearchEntry2->sextractorIndex,
                    pSearchEntry2->MAG_ISO,
                    pSextractor2->NUMBER,
                    pSextractor2->X_IMAGE,
                    pSextractor2->Y_IMAGE);

          }
        }
#endif
      } /* Star bright enough for insertion */
#if 0
      if (sextractorIndex > (2 * SEARCH_TABLE_LENGTH)) {
        fprintf(stderr,"ERROR: early abort\n");
        exit(-1);
      }
#endif
    } /* Sextractor table index */
#if 0
    for (insertIndex = 0; insertIndex < SEARCH_TABLE_LENGTH; insertIndex++) {
      pSearchEntry2 = &searchTable[insertIndex];
      if (pSearchEntry2->sextractorIndex >= 0) {
        pSextractor2 = &sextractor_table[pSearchEntry2->sextractorIndex];
        fprintf(stderr,"Entry %3d  Index %7d MAG_ISO %10f NUMBER %6d X_IMAGE %10f Y_IMAGE %10f\n",
                insertIndex,
                pSearchEntry2->sextractorIndex,
                pSearchEntry2->MAG_ISO,
                pSextractor2->NUMBER,
                pSextractor2->X_IMAGE,
                pSextractor2->Y_IMAGE);

      }
    }
#endif
    /* Now go through the sextractor list for each of these objects and save likely candidates */
    for (searchIndex = 0; searchIndex < searchTableSize; searchIndex++) {
      pSearchEntry1 = &searchTable[searchIndex];
      minDec = (pSearchEntry1->Y_IMAGE  * degreesPerPixel) - pWedgeEntry->separation - (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad);
      maxDec = (pSearchEntry1->Y_IMAGE  * degreesPerPixel) + pWedgeEntry->separation + (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad);
      minRa =  (pSearchEntry1->X_IMAGE  * degreesPerPixel) - pWedgeEntry->separation - (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad);
      maxRa =  (pSearchEntry1->X_IMAGE  * degreesPerPixel) + pWedgeEntry->separation + (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad);
      for (sextractorIndex = 0; sextractorIndex < sextractor_nrecs; sextractorIndex++) {
        pSextractorDim = &sextractor_table[sextractorIndex];

        if ((sextractorIndex > 0) &&
            (pSextractorDim->Y_IMAGE < pSextractor2->Y_IMAGE)) {
          fprintf(stderr,"ERROR: sextractor table is not sorted by Y_IMAGE %f %f %s\n",
                  pSextractorDim->Y_IMAGE,
                  pSextractor2->Y_IMAGE,
                  sextractor_name);
          exit(-1);
        }


        pSextractor2 = pSextractorDim;


        if ((pSextractorDim->Y_IMAGE * degreesPerPixel) > maxDec) {
          break;
        }
        if (((pSextractorDim->Y_IMAGE * degreesPerPixel) < minDec) ||
            ((pSextractorDim->X_IMAGE * degreesPerPixel) < minRa) ||
            ((pSextractorDim->X_IMAGE * degreesPerPixel) > maxRa)) {
          continue;
        }
        /* Object is close, perform a more refined search */
        drad = sqrt(sqr((pSextractorDim->X_IMAGE * degreesPerPixel) - (pSearchEntry1->X_IMAGE * degreesPerPixel))
                    +(sqr((pSextractorDim->Y_IMAGE * degreesPerPixel) - (pSearchEntry1->Y_IMAGE * degreesPerPixel))));
        if ((drad > pWedgeEntry->separation + (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad)) ||
            (drad < pWedgeEntry->separation - (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad))) {
          continue;

        }
        /* Here we have a potential matching object */
        yValue =  (pSextractorDim->Y_IMAGE * degreesPerPixel) - (pSearchEntry1->Y_IMAGE * degreesPerPixel);
        xValue =  (pSextractorDim->X_IMAGE * degreesPerPixel) - (pSearchEntry1->X_IMAGE * degreesPerPixel);

        if (xValue == 0.0) {
          if (yValue > 0) {
            angle = 90.0;
          } else {
            angle = -90.0;
          }
        } else {
          angle = atan(yValue/xValue)*RAD_TO_DEGREES;
          if ((xValue < 0) && (yValue < 0)) {
            angle +=  180.0; 
          } else if ((xValue < 0) && (yValue > 0)) {
            angle += 180.0; 
          } else if ((xValue > 0) && (yValue < 0)) {
            angle += 360.0;
          } 
        }
#if 0
        fprintf(stderr,"Candidate at separation %15f angle %15f %15f %15f\n",3600*drad,3600*xValue,3600*yValue,angle);
#endif
        if (doPlots) {
          fprintf(gmtHandle[0],"%f %f\n",3600*xValue,3600*yValue);
        }
        angleBinNumber = angleBins * angle/360.0;
        if (angleBinNumber < 0) {
          angleBinNumber = 0;
        }
        if (angleBinNumber >= angleBins) {
          angleBinNumber = angleBins-1;
        }
        angleBinTable[angleBinNumber]++;
        candidateCount++;
      }
    }
    averageBin = (2*candidateCount)/angleBins;

    for (index = 0; index < angleBins; index+=2) {
      binValue = (angleBinTable[index] + angleBinTable[index+1]);
      binAngle = ((index+1) * 360.0)/angleBins;
      if (doPlots) {
        fprintf(gmtHandle[1],"%f %f\n",binAngle,binValue);
      }
      sigmaSum += 1.0 * binValue;
      sigmaSqr += 1.0 * binValue * binValue;
    }
    binMean = 2*sigmaSum/angleBins;
    binStd = sqrt((2*sigmaSqr/angleBins) - (binMean * binMean));
    fprintf(stderr,"averageBin %d, binMean %f binStd %f\n",averageBin,binMean,binStd);
    /* Now count the number of peaks */
    abovePeak = 0;
    for (index1 = 0; index1 < angleBins; index1+=2) {
      index2 = index1+1;
      binValue = (angleBinTable[index1] + angleBinTable[index2]);      
      excess = binValue-averageBin - (3*binStd);
      if (excess > 0) {
        abovePeak = 1;
      } else {
        if (abovePeak) {
          numPeaks++;
        }
        abovePeak = 0;
      }
      if (index1 == 0) {
        initAbovePeak = abovePeak;
      } else if (index1 == (angleBins -2)) {
        finalAbovePeak = abovePeak;
      }
    }
    /* check for wrap */
    if (finalAbovePeak == 1) {
      if (initAbovePeak == 0) {
        numPeaks;
      }
    }



    /* We are now going to perform two staggered searches to be sure that the image is not
       on a bin boundary */

    for (index1 = 0; index1 < angleBins; index1+=2) {
      index2 = index1+1;
      binValue = (angleBinTable[index1] + angleBinTable[index2]);      
      excess = binValue-averageBin - (3*binStd);
      if (excess > maxExcess) {
        maxExcess = excess;
        maxIndex1 = index1;
        maxIndex2 = index2;
      }
    }
    for (index1 = 1; index1 < angleBins; index1+=2) {
      index2 = index1+1;
      if (index2 >= angleBins) {
        index2 = index2-angleBins;
      }
      binValue = (angleBinTable[index1] + angleBinTable[index2]);      
      double excess = binValue-averageBin - (3*binStd);
      if (excess > maxExcess) {
        maxExcess = excess;
        maxIndex1 = index1;
        maxIndex2 = index2;
      }
    }
    

  
  

    /* Print our highest value */
    fprintf(stderr,"Angle Bin %3d, Angle %6.1f, Count %4d  Ave %5.1f std %5.1f Excess %5.1f Peaks %d for %s\n",
            (maxIndex1+1),
            ((maxIndex1+1)*360.0)/angleBins,
            (angleBinTable[maxIndex1]+angleBinTable[maxIndex2]),
            binMean,
            binStd,
            maxExcess,
            numPeaks,
            fileroot);
  

    if ((numPeaks == EXPECTED_PEAKS) && (maxExcess > EXCESS_LIMIT)) {
    

      /* The next step is to go through the entire file looking for candidates.  We
         define a search area that covers at least the peak bin and its two neighbors */
      minX = pWedgeEntry->separation + (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad);
      maxX = -minX;
      minY = minX;
      maxY = maxX;
      binAngle = ((maxIndex1-1)) * 360.0/angleBins;
      if (binAngle < 0.0) {
        binAngle += 360.0;
      }

      GetLimits(binAngle,pWedgeEntry->separation + (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad),&minX,&maxX,&minY,&maxY);
      GetLimits(binAngle,pWedgeEntry->separation - (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad),&minX,&maxX,&minY,&maxY);

      binAngle = ((maxIndex1+2)) * 360.0/angleBins;
      if (binAngle >= 360.0) {
        binAngle -= 360.0;
      }

      GetLimits(binAngle,pWedgeEntry->separation + (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad),&minX,&maxX,&minY,&maxY);
      GetLimits(binAngle,pWedgeEntry->separation - (DRAD_REJECT_FACTOR * pWedgeEntry->meanplussigma_drad),&minX,&maxX,&minY,&maxY);

      if (verbose) {
        printf("Search range is ra: %f %f dec %f %f \n",3600*minX,3600*maxX,3600*minY,3600*maxY);
      }
      /* At this point, we now repeat the search for all objects, looking for shadow images in the
         restricted search range.  Take advantage of the declination sorted table to reduce the
         search time */
  
      baseIndex = 0;

      for (searchIndex = 0; searchIndex < sextractor_nrecs; searchIndex++) {
        pSextractorBright = &sextractor_table[searchIndex];
#if 0
        if ((pSextractorBright->FWHM_WORLD / 2.0) > pWedgeEntry->separation) {
          continue;
        }
        if (pSextractorBright->MAG_ISO > searchTable[SEARCH_TABLE_LENGTH-1].MAG_ISO) {
          continue;
        }
#endif


        startIndex = matchArraySize;  /* Start of list for this object */

        minDec = (pSextractorBright->Y_IMAGE * degreesPerPixel) + minY;
        maxDec = (pSextractorBright->Y_IMAGE * degreesPerPixel) + maxY;


        sextractorIndex = baseIndex;
        while (sextractorIndex < sextractor_nrecs) {
          pSextractorDim = &sextractor_table[sextractorIndex];


          sextractorIndex++;

          if ((pSextractorDim->Y_IMAGE * degreesPerPixel) < minDec) {
            /* Declination is too small, step our base index forward
               for the next search */
            baseIndex = sextractorIndex;
            continue;
          }
          if ((pSextractorDim->Y_IMAGE * degreesPerPixel) > maxDec) {
            /* Exit the loop if our declination is too large */
            break;
          }
          if (pSextractorDim->MAG_ISO <= pSextractorBright->MAG_ISO) {
            /* Object is brighter than ghost image, skip it */
            continue;
          }
          yValue =  (pSextractorDim->Y_IMAGE * degreesPerPixel) - (pSextractorBright->Y_IMAGE * degreesPerPixel);
          xValue =  (pSextractorDim->X_IMAGE * degreesPerPixel) - (pSextractorBright->X_IMAGE * degreesPerPixel);




          if ((xValue < minX) ||
              (xValue > maxX) ||
              (yValue < minY) ||
              (yValue > maxY)) {
            /* Outside the box, skip it */
            continue;
          }
#if 0
          fprintf(stderr,"Candidate at x %15f y %15f \n",3600*xValue,3600*yValue);
#endif

          if (matchArraySize >= matchArrayAlloc) {
            matchArrayAlloc += MATCH_ARRAY_ALLOC;
            pNewMatchArray = realloc(pMatchArray,(matchArrayAlloc * sizeof(MATCHENTRY)));
            if (pNewMatchArray == NULL) {
              fprintf(stderr,"ERROR pMatchArray realloc failed\n");
            }
            pMatchArray = pNewMatchArray;
        
          }
          pMatchEntry = &pMatchArray[matchArraySize];
          memset(pMatchEntry,0,sizeof(MATCHENTRY));
          pMatchEntry->brightIndex = searchIndex;
          pMatchEntry->dimIndex = sextractorIndex-1;
          pMatchEntry->xValue = xValue;
          pMatchEntry->yValue = yValue;
#if 0
          pMatchEntry->magsum = pSextractorDim->MAG_ISO - pSextractorBright->MAG_ISO;
#else
          pMatchEntry->magsum = - pSextractorDim->MAG_ISO - pSextractorBright->MAG_ISO;
          pMatchEntry->magdiff = pSextractorDim->MAG_ISO - pSextractorBright->MAG_ISO;
#endif
          matchArraySize++;
#if 0
          if  ((sextractor_table[pMatchEntry->dimIndex].MAG_ISO -
                sextractor_table[pMatchEntry->brightIndex].MAG_ISO) < 0.0) {
            printf("Bright %8d MAG_ISO: %15f Dim %8d MAG_ISO %15f Difference %15f\n",
                   pMatchEntry->brightIndex,
                   sextractor_table[pMatchEntry->brightIndex].MAG_ISO,
                   pMatchEntry->dimIndex,
                   sextractor_table[pMatchEntry->dimIndex].MAG_ISO,
                   sextractor_table[pMatchEntry->dimIndex].MAG_ISO -
                   sextractor_table[pMatchEntry->brightIndex].MAG_ISO);
          }
#endif
        }
        endIndex = matchArraySize;
        if (endIndex > startIndex) {
          /* Find the brightest companion for this star */
          pMatchEntry = &pMatchArray[startIndex];
          pSextractorDim = &sextractor_table[pMatchEntry->dimIndex];
          brightestMagIsoIndex = startIndex;
          brightestMagIso = pSextractorDim->MAG_ISO;
          for (matchIndex = startIndex+1; matchIndex < endIndex; matchIndex++) {
            pMatchEntry = &pMatchArray[matchIndex];
            pSextractorDim = &sextractor_table[pMatchEntry->dimIndex];
            if (pSextractorDim->MAG_ISO < brightestMagIso) {
              brightestMagIso = pSextractorDim->MAG_ISO;
              brightestMagIsoIndex = matchIndex;
            }
          }
          /* Select this star for further consideration */
          pMatchEntry = &pMatchArray[brightestMagIsoIndex];
          pMatchEntry->selected = 1;
          selectedCount++;
#if 0
          if (doPlots) {
            fprintf(gmtHandle[2],"%f %f\n",3600*pMatchEntry->xValue,3600*pMatchEntry->yValue);
            fprintf(gmtHandle[3],"%f %f\n",sextractor_table[pMatchEntry->brightIndex].MAG_ISO,pMatchEntry->magdiff);
          }
#endif      
        }



      }
      if (verbose) {
        printf("Found %d objects in the region of interest\n",matchArraySize);
      }

      /* Now sort on the largest magnitude differences */
      if (pMagsumIndexList != NULL) {
        free(pMagsumIndexList);
        pMagsumIndexList = NULL;
      }


      pMagsumIndexList = (PMAGSUM_INDEX)calloc(matchArraySize,sizeof(MAGSUM_INDEX));
      for (index = 0; index < matchArraySize; index++) {
        pMatchEntry  = &pMatchArray[index];
        pMagsumIndex = &pMagsumIndexList[index];
        pMagsumIndex->index = index;
        pMagsumIndex->magsum = pMatchEntry->magsum;
      }



      qsort(pMagsumIndexList,matchArraySize,sizeof(MAGSUM_INDEX),magsumCompare);
      if (pSelectedMatchArray != NULL) {
        free(pSelectedMatchArray);
        pSelectedMatchArray = NULL;
      }

      pSelectedMatchArray = (PMATCHENTRY)calloc(selectedCount,sizeof(MATCHENTRY));
      {
        double xAve;
        double yAve;
        double xAve2;
        double yAve2;
        int selectCount2;
        int selectCount3;  
        double dradAve = 0;
        double dradAve2 = 0;
        double dradSqr = 0;
        double minDrad;
        int minDradIndex = -1;
        int selectCount;
        int selectIndex;
        double aveMagDiff = 0;
        double aveMAG_ISO = 0;
        selectCount = 0;
        for (matchIndex = 0; matchIndex < matchArraySize; matchIndex++) {
          pMagsumIndex = &pMagsumIndexList[matchIndex];
          pMatchEntry = &pMatchArray[pMagsumIndex->index];
          if (pMatchEntry->selected == 0) {
            continue;
          }
          pSextractorBright = &sextractor_table[pMatchEntry->brightIndex];
      
          pSelectedMatchEntry = &pSelectedMatchArray[selectCount];
          memcpy(pSelectedMatchEntry,pMatchEntry,sizeof(MATCHENTRY));
          selectCount++;


          if (selectCount < 10) {
            continue;
          } 

      

          xAve = 0;
          yAve = 0;
          /* Now we compute the mean position of selected objects */
          for (index = 0; index < selectCount; index++) {
            pSelectedMatchEntry = &pSelectedMatchArray[index];
            xAve += pSelectedMatchEntry->xValue;
            yAve += pSelectedMatchEntry->yValue;
          }
          xAve = xAve/selectCount;
          yAve = yAve/selectCount;
          /* Now compute the radial difference */
          dradAve = 0;
          for (index = 0; index < selectCount; index++) {
            pSelectedMatchEntry = &pSelectedMatchArray[index];
            pSelectedMatchEntry->drad = sqrt(sqr(pSelectedMatchEntry->xValue - xAve) + sqr(pSelectedMatchEntry->yValue - yAve));
            dradAve += pSelectedMatchEntry->drad;
            dradSqr += sqr(pSelectedMatchEntry->drad);
          }
          dradAve = dradAve/selectCount;
          dradSqr = sqrt(dradSqr/selectCount);
          if ((minDradIndex < 0) || (dradAve < minDrad)) {
            minDrad = dradAve;
            minDradIndex = selectCount;  
          }
#if 0
          printf("Count %5d, xAve %7.1f yAve %7.1f, drad %7.1f, drad RMS_0 %7.1f\n",
                 selectCount,3600*xAve,3600*yAve,3600*dradAve,3600*dradSqr);
#endif
          if ((selectCount >= 1000) && 
              (minDradIndex != selectCount)) {
            break;
          }

        }
        /* Now cast out the outliers and recompute the centroid */
        xAve2 = 0;
        yAve2 = 0;
        selectCount2 = 0;
        for (selectIndex = 0; selectIndex <= minDradIndex; selectIndex++) {
          pSelectedMatchEntry = &pSelectedMatchArray[selectIndex];
          if (pSelectedMatchEntry->drad > dradAve) {
            continue;
          }
          xAve2 += pSelectedMatchEntry->xValue;
          yAve2 += pSelectedMatchEntry->yValue;
          aveMagDiff += pSelectedMatchEntry->magdiff;
          aveMAG_ISO += sextractor_table[pSelectedMatchEntry->brightIndex].MAG_ISO;
          selectCount2++;

        }


        xAve = xAve2/selectCount2;
        yAve = yAve2/selectCount2;
        aveMagDiff = aveMagDiff/selectCount2;
        aveMAG_ISO = aveMAG_ISO/selectCount2;
        /* We find the intersection of these two curves */
        /*  magdiff = MAG_ISO_med - MAG_ISO_crit (-1 slope assumed) */
        /*  magdiff = K + MAG_ISO_crit           (+1 slope assumed for a nonlinear relationship) */
        /*  K is found by the ave point calculated above: */
        /*  aveMagDiff = K + aveMAG_ISO */


        MAG_ISO_crit = (MAG_ISO_med - aveMagDiff + aveMAG_ISO)/2.0;

        /* Two straight lines showing the critical parameters */

        if (doPlots) {
          fprintf(gmtHandle[6],"%f %f\n",MAG_ISO_med,0.0);
          fprintf(gmtHandle[6],"%f %f\n",MAG_ISO_crit,(MAG_ISO_med - MAG_ISO_crit));
          fprintf(gmtHandle[6],"%f %f\n",aveMAG_ISO,aveMagDiff);
        }
#if 0
        printf("New centroid x: %7.1f  y: %7.1f dist %7.1f from %d of %d objects drad=max( %7.1f %7.1f ) MAG_ISO median %6.2f rms %6.2f aveMagDiff %6.2f aveMAG_ISO %6.2f critical MAG_ISO %6.2f\n",
               3600*xAve,
               3600*yAve,
               3600*sqrt(sqr(xAve)+sqr(yAve)),
               selectCount2,
               selectCount,
               3600*minDrad,
               3600*pWedgeEntry->meanplussigma_drad,
               MAG_ISO_med,
               MAG_ISO_rms,
               aveMagDiff,
               aveMAG_ISO,
               MAG_ISO_crit);
#endif
        if (minDrad < pWedgeEntry->meanplussigma_drad) {
          minDrad = pWedgeEntry->meanplussigma_drad;
        }
        /* At this point, we search for everything within three times our new drad */
        minDrad = minDrad * DRAD_REJECT_FACTOR;

        selectCount3 = 0;
        baseIndex = 0;

        for (searchIndex = 0; searchIndex < sextractor_nrecs; searchIndex++) {
          pSextractorBright = &sextractor_table[searchIndex];
          pSextractorBright->MATCH_NUMBER = 0;
        }
        for (searchIndex = 0; searchIndex < sextractor_nrecs; searchIndex++) {
          pSextractorBright = &sextractor_table[searchIndex];
    

          minDec = (pSextractorBright->Y_IMAGE  * degreesPerPixel) + yAve - minDrad;
          maxDec = (pSextractorBright->Y_IMAGE  * degreesPerPixel) + yAve + minDrad;



          sextractorIndex = baseIndex;
          while (sextractorIndex < sextractor_nrecs) {
            pSextractorDim = &sextractor_table[sextractorIndex];


            sextractorIndex++;

            if ((pSextractorDim->Y_IMAGE * degreesPerPixel) < minDec) {
              /* Declination is too small, step our base index forward
                 for the next search */
              baseIndex = sextractorIndex;
              continue;
            }
            if ((pSextractorDim->Y_IMAGE * degreesPerPixel) > maxDec) {
              /* Exit the loop if our declination is too large */
              break;
            }
            if (pSextractorDim->MAG_ISO <= pSextractorBright->MAG_ISO) {
              /* Object is brighter than ghost image, skip it */
              continue;
            }
            minRa =  (pSextractorBright->X_IMAGE  * degreesPerPixel) + xAve - minDrad;
            maxRa =  (pSextractorBright->X_IMAGE  * degreesPerPixel) + xAve + minDrad;
            if (((pSextractorDim->X_IMAGE * degreesPerPixel) < minRa) ||
                ((pSextractorDim->X_IMAGE * degreesPerPixel) > maxRa)) {
              continue;
            }

            /* Object is close, perform a more refined search */
            yValue =  (pSextractorDim->Y_IMAGE * degreesPerPixel) - (pSextractorBright->Y_IMAGE * degreesPerPixel);
            xValue =  (pSextractorDim->X_IMAGE * degreesPerPixel) - (pSextractorBright->X_IMAGE * degreesPerPixel);
            drad = sqrt(sqr(xValue-xAve) + sqr(yValue-yAve));
            if ((drad < minDrad) && 
                (pSextractorBright->MAG_ISO < (MAG_ISO_crit+MAG_ISO_CRIT_MARGIN)) &&
                (pSextractorDim->MAG_ISO < (MAG_ISO_med - MAG_ISO_rms))) {
              /* We have a match */
              selectCount3++;
              dradAve2 += drad;
              pSextractorDim->MATCH_NUMBER = pSextractorBright->NUMBER;
#if 0
              printf("Index %5d Bright MAG_ISO %6.2f Dim MAG_ISO %6.2f magsum %6.2f  X %7.1f  Y %7.1f \n",
                     selectCount3,
                     pSextractorBright->MAG_ISO,
                     pSextractorDim->MAG_ISO,
                     - pSextractorDim->MAG_ISO - pSextractorBright->MAG_ISO,
                     3600*xValue,
                     3600*yValue);
#endif
              if (doPlots) {
                fprintf(gmtHandle[4],"%f %f\n",3600*xValue,3600*yValue);
                fprintf(gmtHandle[5],"%f %f\n",pSextractorBright->MAG_ISO,pSextractorDim->MAG_ISO - pSextractorBright->MAG_ISO);
              }
            } else {
              if (doPlots) {
                fprintf(gmtHandle[2],"%f %f\n",3600*xValue,3600*yValue);
                fprintf(gmtHandle[3],"%f %f\n",pSextractorBright->MAG_ISO,pSextractorDim->MAG_ISO - pSextractorBright->MAG_ISO);
              }

            }
                  



          }
        }
        plotCount = selectCount3;
        if (plotCount > 0) {
          dradAve2 = dradAve2/plotCount;
        }

        printf("New centroid x: %7.1f  y: %7.1f dist %7.1f from %d of %d objects drad=max( %7.1f %7.1f ) final drad %7.1f MAG_ISO median %6.2f rms %6.2f aveMagDiff %6.2f aveMAG_ISO %6.2f critical MAG_ISO %6.2f final drad %7.1f\n",
               3600*xAve,
               3600*yAve,
               3600*sqrt(sqr(xAve)+sqr(yAve)),
               selectCount2,
               selectCount,
               3600*minDrad,
               3600*pWedgeEntry->meanplussigma_drad,
               3600*dradAve2,
               MAG_ISO_med,
               MAG_ISO_rms,
               aveMagDiff,
               aveMAG_ISO,
               MAG_ISO_crit);

      }


      /* Now write the output file - use same order as the input file for 
         a table join operation */
      outHandle = fopen(outfile,"wt");
      if (outHandle == NULL) {
        fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
        exit(-1);
      } else {
        if (verbose) {
          fprintf(stderr,"Output file %s\n",outfile);
        }
      }
      fprintf(outHandle,"NUMBER\tMATCH_NUMBER\tAFLAGS\n");
      fprintf(outHandle,"------\t------------\t------\n");

      for (searchIndex = 0; searchIndex < sextractor_nrecs; searchIndex++) {
        pSextractor1 = &sextractor_table[searchIndex];

        if (pSextractor1->MATCH_NUMBER == 0) {
          /* Normal star, be sure that the wedge flag bit is clear */
          pSextractor1->AFLAGS &= (~(1<<FILTER_AFLAG_WEDGE));
        } else {

          /* A wedge star, set the wedge bit */
          pSextractor1->AFLAGS |= (1<<FILTER_AFLAG_WEDGE);

        }
        fprintf(outHandle,"%d\t%d\t%d\n",pSextractor1->NUMBER,pSextractor1->MATCH_NUMBER,pSextractor1->AFLAGS);
      }
      wedgeFlag = 1;
      break;
    } else {
      wedgeFlag = 0;
    }
  }

  if (wedgeFlag == 1) {
    SetMosaicFitWCS2("Wedge",fileroot,solutionNumber,1,readFlag,&FitWCS);
    SetPlateQuality2("wedge",fileroot,1,readFlag,&quality);
    printf("Pickering Wedge found for %s\n",fileroot);
  } else if (wedgeFlag == 0) {
    SetMosaicFitWCS2("Wedge",fileroot,solutionNumber,0,readFlag,&FitWCS);
    SetPlateQuality2("wedge",fileroot,0,readFlag,&quality);
    printf("No Pickering Wedge detected for %s\n",fileroot);
  }

  /* All done.  Clean up */
  if (angleBinTable != NULL) {
    free(angleBinTable);
  }
 
  if (sextractor_table != NULL) {
    Free(sextractor_table);
  }
  if (sextractor_header != NULL) {
    table_hdrfree(sextractor_header);
  }
  if (sextractor_handle != NULL) {
    Close(sextractor_handle);
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }
  if (doPlots) {
    for (gmtIndex = 0; gmtIndex < NUM_GMT_FILES; gmtIndex++) {
      if (gmtHandle[gmtIndex] != NULL) {
        fclose(gmtHandle[gmtIndex]);
      }
    }
  }
  if (pMatchArray != NULL) {
    free(pMatchArray);
  }
  if (pSelectedMatchArray != NULL) {
    free(pSelectedMatchArray);
  }
  if (pMagsumIndexList != NULL) {
    free(pMagsumIndexList);
  }
  if (vector != NULL) {
    free(vector);
  }

  time(&curTime);
  curTime -= startTime;
  fprintf(stdout,"Candidates: %d MatchArray: %d selected: %d plotted %d seconds %d for %s\n",
          candidateCount,
          matchArraySize,
          selectedCount,
          plotCount,
          curTime,
          fileroot);
    


  return(0);
}

