// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* matchgsc.c
 *
 *  Match an extract from the SDSS catalog with the GSC2.3.2 catalog
 * 
 *   gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm  -L/usr/lib64/mysql -lmysqlclient matchgsc.c pipelineutils.a -ltable -lutil -lwcs -o matchgsc 
 * 
 *  matchgsc -v -i /dasch/backup/2009_03_26/sdss.db -o /dasch/Pipeline/sdss.tmp
 *
 *  Mar 31, 2009 Edward J. Los  Initial Version
 *  Apr  4, 2009 Edward J. Los  Add flags statistics
 *  May 17, 2009 Edward J. Los - Allow multiple GSC bin index sizes
 *
 */   

#include <math.h>
#include <errno.h>
#include "table.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#define MAX_SDSS_ID_LENGTH 25
#define MAX_BUFFER 100
#define CLIP_DRAD 1.0 /* Clipping radius in arcsec */
#define CLIP_MAGDIFF 1.5 /* Clipping magnitude difference */
#define STAR_TYPE 6
#define GALAXY_TYPE 3
#define MAX_FLAGS 16

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

typedef struct _sdssdata {
  char objid[MAX_SDSS_ID_LENGTH];  /* SDSS ID */
  int  type;          /* Type 3 = galaxy and 6 = star */
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
  double gmag;        /* Green magnitude */
  double rmag;        /* Red magnitude */
} SDSSDATA,*PSDSSDATA;
int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char outfile[MAX_BUFFER];
  char cmdchar;
  int errorFlag = 0;
  int verbose = 0;
  FILE *outHandle = NULL;
  int outCount = 0;
  double magdiff;
  int sdssStarCount = 0;
  int sdssStarOutCount = 0;

  char indexname[MAX_BUFFER];
  File indexHandle;
  char* catalogname;
  File catalogHandle;

  int gsc_bin_index;
  int decBin;
  int raBin;
  int cur_gsc_bin_index = -1;
  STARINDEX curStarIndex;
  PSTARINDEX pCurStarIndex = &curStarIndex;


  PGSCIMAGE pGscImageTable = NULL;
  int gscImageAlloc = 0;
  PGSCIMAGE pCurGscImage = NULL;
  PGSCIMAGE pBestGscImage;
  int readItems;
  int gscIndex;
  double bestDrad;
  double curDrad;
  double factor;
  char *dotPtr;
	char REF[MAX_REF];

  File sdss_handle = NULL;
  char sdss_name[MAX_BUFFER];
  TableHead sdss_header = NULL;
  PSDSSDATA sdss_table = NULL;
  size_t sdss_nrecs = 0;
  int sdss_index;
  PSDSSDATA pSDSSEntry;
  int index;
  int flagsCount[MAX_FLAGS];

  for (index = 0; index < MAX_FLAGS; index++) {
    flagsCount[index] = 0;
  }

  outfile[0] = 0;
  sdss_name[0] = 0;

  catalogname = getenv("DASCH_CATALOG");
  if (catalogname == NULL) {
    fprintf(stderr,"DASCH_CATALOG is not defined\n");
    printf("0\n");
    return(-1);
  }

  strcpy(indexname,catalogname);
  dotPtr = strstr(indexname,".");
  if (dotPtr != NULL) {
    *dotPtr = 0;
  }
  strcat(indexname,".idx");



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
        case 'i': /* input file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sdss_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;





        default:
          fprintf(stderr,"ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }

  if (sdss_name[0] == 0) {
    fprintf(stderr,"ERROR: No match filename was specified\n");
    errorFlag = 1;
  }
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
  fprintf(outHandle,"REF\tdrad\tmagdiff\tsdsscolor\tgsccolor\tMAGFlag\tgmag\n");
  fprintf(outHandle,"---\t----\t-------\t---------\t--------\t-------\t----\n");

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



  /* Open the match file */
  sdss_handle = Open(sdss_name,"r");
  if (sdss_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the match file %s\n",sdss_name);
  } else {
    if (verbose) {
      fprintf(stderr,"Found match file %s\n",sdss_name);
    }
  }


  if (errorFlag) {
    fprintf(stderr,"Usage: matchgsc -i <sdss file> -o <output file> \n");

    return(-1);
  }

  printf("matchgsc of %s %s \n",
         __DATE__,__TIME__);




  sdss_header = table_header(sdss_handle,TABLE_PARSE);
  if (sdss_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",sdss_name);
    return(-1);
  }

  sdss_table = table_loadva(sdss_handle,
                            &sdss_header,
                            NULL, /* hbase */
                            NULL, /* rows */
                            NULL,
                            sizeof(SDSSDATA),
                            &sdss_nrecs,
                            TblInt,"type" ,TblOff(PSDSSDATA,type),
                            TblDbl,"ra"   ,TblOff(PSDSSDATA,ra),
                            TblDbl,"dec"  ,TblOff(PSDSSDATA,dec),
                            TblDbl,"g" ,TblOff(PSDSSDATA,gmag),
                            TblDbl,"r"  ,TblOff(PSDSSDATA,rmag),
                            TblBuf,"objid"    ,TblOff(PSDSSDATA,objid),MAX_SDSS_ID_LENGTH,
                            0,"end",0);
  if (sdss_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",sdss_name);
    return(-1);
  }
  if (verbose) {
    fprintf(stderr,"read %d records for %s\n",sdss_nrecs,sdss_name);
  }

  for (sdss_index = 0; sdss_index < sdss_nrecs; sdss_index++) {
    pSDSSEntry = &sdss_table[sdss_index];
    if (pSDSSEntry->type == STAR_TYPE) {
      sdssStarCount++;
    }
    gsc_bin_index = GetGSCBin(pGscBin,pSDSSEntry->ra,pSDSSEntry->dec,&decBin,&raBin,"matchgsc");
    if (gsc_bin_index != cur_gsc_bin_index) {
      cur_gsc_bin_index = gsc_bin_index;
      Seek(indexHandle,cur_gsc_bin_index * sizeof(STARINDEX),SEEK_SET);
      readItems = Read(indexHandle,pCurStarIndex,sizeof(STARINDEX),1);
      if (readItems != 1) {
        fprintf(stderr,"ERROR reading star index file\n");
        exit(-1);
      }
      if (pCurStarIndex->numStars > gscImageAlloc) {
        gscImageAlloc = pCurStarIndex->numStars+1000;
        pCurGscImage = realloc(pGscImageTable,gscImageAlloc*sizeof(GSCIMAGE));
        if (pCurGscImage == NULL) {
          fprintf(stderr,"ERROR: failed to realloc pGscImageTable of size %d\n",gscImageAlloc*sizeof(GSCIMAGE));
          exit(-1);
        }
        pGscImageTable = pCurGscImage;
        pCurGscImage = NULL;
      }
      Seek(catalogHandle,pCurStarIndex->offset,SEEK_SET);
      readItems = Read(catalogHandle,pGscImageTable,sizeof(GSCIMAGE),pCurStarIndex->numStars);
      if (readItems != pCurStarIndex->numStars) {
        fprintf(stderr,"ERROR reading gsc catalog file\n");
        exit(-1);
      }
      pBestGscImage = NULL;
      double bestDrad = 360.0;
      factor = cos(DEGREES_TO_RAD*(pSDSSEntry->dec));
      for (gscIndex = 0; gscIndex < pCurStarIndex->numStars; gscIndex++) {
        pCurGscImage = &pGscImageTable[gscIndex];
        curDrad = sqrt(sqr(factor*(pSDSSEntry->ra-pCurGscImage->ra)) + sqr(pSDSSEntry->dec-pCurGscImage->dec));
        if (curDrad < bestDrad) {
          bestDrad = curDrad;
          pBestGscImage = pCurGscImage;
        }
      }
      if (pBestGscImage != NULL) {
        bestDrad = bestDrad * 3600.0;
        magdiff = pBestGscImage->Stdmag-pSDSSEntry->gmag;
        if (magdiff < 0) {
          magdiff = - magdiff;
        }
        if ((bestDrad < CLIP_DRAD) && 
            (magdiff < CLIP_MAGDIFF) &&
            (pBestGscImage->color < 90.0)) {
          outCount++;
#if 0
          printf("REF %20s bestDrad %.5f, magDiff %.2f\n",
                 pBestGscImage->REF,
                 bestDrad,
                 magdiff);
#endif
          if (pSDSSEntry->type == STAR_TYPE) {
            sdssStarOutCount++;
          }
					GetREF(pBestGscImage->REFNumber,REF,0,1);
          fprintf(outHandle,"%s\t%f\t%f\t%f\t%f\t%d\t%f\n",
                  REF,
                  bestDrad,
                  magdiff,
                  pSDSSEntry->gmag-pSDSSEntry->rmag,
                  pBestGscImage->color,
                  pBestGscImage->MAGFlag,
                  pSDSSEntry->gmag);
          if ((pBestGscImage->MAGFlag >= 0)  &&
              (pBestGscImage->MAGFlag < MAX_FLAGS)) {
            flagsCount[pBestGscImage->MAGFlag]++;
          }
  
        }
      }
    }
  }



  /* We are done with the match table.  Clear memory */
  if (sdss_table != NULL) {
    Free(sdss_table);
    sdss_table = NULL;
  }
  if (sdss_header != NULL) {
    table_hdrfree(sdss_header);
    sdss_header = NULL;
  }
  if (sdss_handle != NULL) {
    Close(sdss_handle);
    sdss_handle = NULL;
  }
  pSDSSEntry = NULL;
  Close(catalogHandle);
  Close(indexHandle);
  if (pGscImageTable != NULL) {
    free(pGscImageTable);
  }
  printf("SDSS objects read %d, stars %d  GSC objects output %d stars output %d\n",
         sdss_nrecs,
         sdssStarCount,
         outCount,
         sdssStarOutCount);
  for (index = 0; index < MAX_FLAGS; index++) {
    if (flagsCount[index] != 0) {
      printf("Flag %2d, count %10d\n",index,flagsCount[index]);
    }
  }


  fclose(outHandle);
}
