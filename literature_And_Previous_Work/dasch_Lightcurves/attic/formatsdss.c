// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatsdss.c
 *
 * Read in an extract of the SDSS catalog in accordance with Josh's memorandum of Mon 5/18/15 10:29 AM and put it into GALAXYREC format
    http://tdc-www.harvard.edu/computing/sdss/
    http://tdc-www.harvard.edu/computing/sdss/photo.html
 *
 *  Condition 1: select   mag range g = 17 - 20 and color range g - r = 0.5 - 1.0 objects: 
 *  date;row ' ((g - extinction_g) < 20)  && ((g - extinction_g) > 17) && (((g - extinction_g)-(r - extinction_r)) > 0.5) && (((g - extinction_g)-(r - extinction_r)) < 1.0)' < /data/astrocat2/SDSS-dr7/PhotoPrim.db | column  objID ra dec u g r i z extinction_u extinction_g extinction_r exti
nction_i extinction_z petroMag_g_r type status flags flags2 > sdssdr7.db;date
 *  Condition 2: Use only type = 6 (star) records
 *  Condition 3: Put in a requirement that the stars examined be at least 1arcmin from any bright galaxy with total mag (if SDSS lists this) with g <14. by using the existing galaxy.dat database.
 *
 *
 *    formatsdss -o /dasch/Pipeline/catalogs/sdss.db -b  /dasch/Pipeline/catalogs/sdsstemp.dat
 *
 *    formatapass of May 26 2015 11:05:13 Output Filename /dasch/Pipeline/catalogs/sdss.db Binary Filename /dasch/Pipeline/catalogs/sdss.dat
 *    Size of BININDEX is 16.  Size of STARINDEX is 16. Size of SDSSIMAGE is 136 Size of GALAXYREC is 104 
 *    Read 16224192 records from /dasch/Pipeline/catalogs/sdssdr7/sdssdr7.db at 112 seconds
 *    Allocated 16224192 records 112 seconds
 *    Read the galaxy common at 113 seconds
 *    Prepared 12971119 records in 202 seconds
 *    Writing output table at 202 seconds
 *    Wrote 12971119 items of size 104 to /dasch/Pipeline/catalogs/sdss.dat
 *    nonStarCount 3202547 nearbyGalaxyCount 50526
 *    Execution Time: 213 seconds; total entries 16224192 stars written: 12971119
 *
 *    
 *    [elos@dasch8 sdssdr7]$ formatsdss -o /dasch/Pipeline/catalogs/sdss.db -b  /dasch/Pipeline/catalogs/sdsstemp.dat
 *    formatapass of Jun 22 2015 09:40:47 Output Filename /dasch/Pipeline/catalogs/sdss.db Binary Filename /dasch/Pipeline/catalogs/sdss.dat
 *    Size of BININDEX is 16.  Size of STARINDEX is 16. Size of SDSSIMAGE is 136 Size of GALAXYREC is 104 
 *    Read 157137849 records from /dasch/Pipeline/catalogs/sdssdr7/sdssdr7.db at 960 seconds
 *    Allocated 157137849 records 960 seconds
 *    Read the galaxy common at 963 seconds
 *    Prepared 156238055 records in 1980 seconds
 *    Writing output table at 1980 seconds
 *    Wrote 156238055 items of size 104 to /dasch/Pipeline/catalogs/sdss.dat
 *    nonStarCount 0 nearbyGalaxyCount 899794
 *    Execution Time: 1992 seconds; total entries 157137849 stars written: 156238055
 *    
 *    New indexed format:
 *    [scanner@localhost ~/Pipeline]$ formatsdss -o /dasch/Pipeline/catalogs/sdss.db -b  /dasch/Pipeline/catalogs/sdsstemp.dat
 *    formatapass of Jul  6 2015 15:38:16 Output Filename /dasch/Pipeline/catalogs/sdss.db Binary Filename /dasch/Pipeline/catalogs/sdsstemp.dat Sorted Binary Filename /dasch/Pipeline/catalogs/sdss.dat Index filename /dasch/Pipeline/catalogs/sdss.idx
 *    Size of BININDEX is 16.  Size of STARINDEX is 16. Size of SDSSIMAGE is 136 Size of GALAXYREC is 104 
 *    [scanner@localhost ~/Pipeline]$ 
 *    [scanner@localhost ~/Pipeline]$ formatsdss -o /dasch/Pipeline/catalogs/sdss.db -b  /dasch/Pipeline/catalogs/sdsstemp.dat
 *    formatapass of Jul  6 2015 15:38:16 Output Filename /dasch/Pipeline/catalogs/sdss.db Binary Filename /dasch/Pipeline/catalogs/sdsstemp.dat Sorted Binary Filename /dasch/Pipeline/catalogs/sdss.dat Index filename /dasch/Pipeline/catalogs/sdss.idx
 *    Size of BININDEX is 16.  Size of STARINDEX is 16. Size of SDSSIMAGE is 136 Size of GALAXYREC is 104 
 *    Read 16224192 records from /dasch/Pipeline/catalogs/sdssdr7/sdssdr7.db at 111 seconds
 *    Allocated 16224192 records 112 seconds
 *    Read the galaxy common at 113 seconds
 *    Prepared 12971119 records in 210 seconds
 *    Writing output table at 210 seconds
 *    Wrote 12971119 items of size 104 to /dasch/Pipeline/catalogs/sdsstemp.dat
 *    Allocating sort index at 223 seconds
 *    Preparing sort index at 223 seconds
 *    Sorting 116740071 items at 229 seconds
 *    Writing 116740071 items at 245 seconds
 *    nonStarCount 3202547 nearbyGalaxyCount 50526
 *    Execution Time: 428 seconds; total entries 16224192 stars written: 12971119 index entries written 116740071
 *    
 *
 *    [elos@dasch8 ~/Pipeline]$ formatsdss -o /dasch/Pipeline/catalogs/sdss.db -b /dasch/Pipeline/catalogs/sdsstemp.dat
 *    formatapass of Jul  6 2015 16:58:07 Output Filename /dasch/Pipeline/catalogs/sdss.db Binary Filename /dasch/Pipeline/catalogs/sdsstemp.dat Sorted Binary Filename /dasch/Pipeline/catalogs/sdss.dat Index filename /dasch/Pipeline/catalogs/sdss.idx
 *    Size of BININDEX is 16.  Size of STARINDEX is 16. Size of SDSSIMAGE is 136 Size of GALAXYREC is 104 
 *    Read 157137849 records from /dasch/Pipeline/catalogs/sdssdr7/sdssdr7.db at 953 seconds
 *    Allocated 157137849 records 953 seconds
 *    Read the galaxy common at 954 seconds
 *    Prepared 156238055 records in 1985 seconds
 *    Writing output table at 1985 seconds
 *    Wrote 156238055 items of size 104 to /dasch/Pipeline/catalogs/sdsstemp.dat
 *    Allocating sort index at 1995 seconds
 *    Preparing sort index at 1995 seconds
 *    Sorting 1406142495 items at 2059 seconds
 *    Writing 1406142495 items at 2270 seconds
 *    nonStarCount 0 nearbyGalaxyCount 899794
 *    Execution Time: 3095 seconds; total entries 157137849 stars written: 156238055 index entries written 1406142495
 *    
 *
 * May 25, 2015 Edward J. Los - Initial version * Jun 15, 2015 Edward J. Los - Add g and (g-r) to galaxyType so it reads SDSSgXX.XcYY.Y 
 * Jul  6, 2015 Edward J. Los - Because of the expansion of the SDSS catalog, need to store it on disk in the same format as the other catalogs
 *                              We will create sdss.dat sorted by gsc bin index and sdss.idx.  The original binary will be sdsstemp.dat.  Unlike
 *                              the other catalogs, each entry will have a copy in an adjacent bin for the purposes of sorting.
 */ 
  

#define _GNU_SOURCE
#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include <stddef.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <assert.h>
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "kdtree.h"

#define MAX_INPUT_NAME 512
#define MAX_BUFFER 512
#define INPUT_FILES 1
#define MAX_FILTER_SIZE 100
#define GALAXY_DEFAULT_RADIUS 6.0 /* Default galaxy radius in arcsec */

#define SDSS_TYPE_UNKNOWN   0
#define SDSS_TYPE_COSMICRAY 1
#define SDSS_TYPE_DEFECT    2
#define SDSS_TYPE_GALAXY    3
#define SDSS_TYPE_GHOST     4
#define SDSS_TYPE_KNOWNOBJ  5
#define SDSS_TYPE_STAR      6
#define SDSS_TYPE_TRAIL     7
#define SDSS_TYPE_SKY       8

/* #define DEBUG_TYPE 1 */
/* #define DEBUG_SDSSIMAGE 1 */

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
extern char *releaseFieldText[RELEASE_FIELD_MAX+1];
/* #define LOS_DEBUG 1 */  

#define MAX_OBJID_SIZE 19
typedef struct _sdssimage {
  char objID[MAX_OBJID_SIZE+2];
  double ra;
  double dec;
  double u;
  double g;
  double r;
  double i;
  double z;
  double extinction_u;
  double extinction_g;
  double extinction_r;
  double extinction_i;
  double extinction_z;
  int type;
  int status;
  int flags;
  int flags2;
} SDSSIMAGE,*PSDSSIMAGE;


typedef struct _sortindex {
  int gsc_bin_index;
  int record_number;
} SORTINDEX,*PSORTINDEX;

int IndexCompare(const void *first, const void *second)  {
  PSORTINDEX indexFirst = ((PSORTINDEX)first);
  PSORTINDEX indexSecond = ((PSORTINDEX)second);
  if (indexFirst->gsc_bin_index > indexSecond->gsc_bin_index) {
    return(1);
  } else if (indexFirst->gsc_bin_index < indexSecond->gsc_bin_index) {
    return(-1);
  } else {
    return(0);
  }
}

int main(int argc,char *argv[])
{
  char *argstr;
  char output_name[MAX_INPUT_NAME];
  char sdss_name[MAX_INPUT_NAME];
  char sdss_binary_name[MAX_INPUT_NAME];
  char sdss_index_name[MAX_INPUT_NAME];
  FILE* output_handle = NULL;
  File galaxy_handle = NULL;
  File galaxy_binary_handle = NULL;
  File galaxy_index_handle = NULL;
  int totalCount = 0;
  int finalCount = 0;
  int nonStarCount = 0;
  int curRecord = 0;
  int nearbyGalaxyCount = 0;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  char cmdchar;
  int verbose = 0;
  

  SDSSIMAGE sdssimage;
  PSDSSIMAGE pSdssImage = &sdssimage;

  
  PGALAXYREC galaxy_table = NULL;
  GALAXYREC galaxyRec;
  PGALAXYREC pGalaxyRec;
  int writeItems;
  int readItems;
  PSORTINDEX sort_index_table = NULL;
  PSORTINDEX pSortIndex;
  int sortIndexSize = 0;
  int curSortIndex;

  File sdss_output_handle = NULL;
  File sdss_input_handle = NULL;
  char sdss_input_name[] = "/dasch/Pipeline/catalogs/sdssdr7/sdssdr7.db";


  TableHead sdss_header = NULL;
  size_t sdss_nrecs;
  size_t sdss_index;
  PSDSSIMAGE sdss_table;
  PSDSSIMAGE pSdss;
  GALAXYCOMMON galaxyCommon;
  PGALAXYCOMMON pGalaxyCommon = &galaxyCommon;
  GALAXYRESULT galaxyResult;
  PGALAXYRESULT pGalaxyResult = &galaxyResult;
	NEARESTCATALOGSTAR nullCatalogStar;
	PNEARESTCATALOGSTAR pNullCatalogStar = &nullCatalogStar;
	char nearbyObjects[MAX_NEARBY_OBJECTS_STRING];
  char galaxyType[MAX_BUFFER];

  char* catalogdir;
  char* slashPtr;
  char galaxy_name[MAX_BUFFER];

  int gsc_bin_index;
  int cur_bin_index;
  int raBin;
  int decBin;
  int bin_count;
  int bin_index;
  int bin_list[MAX_ADJACENT_BINS+1];
	STARINDEX sdssIndex;
	PSTARINDEX pSdssIndex = &sdssIndex;
  off_t curOffset;


 	InitNearestCatalogStar(pNullCatalogStar);
  assert(MAX_BUFFER > sizeof(SDSSIMAGE));
  assert(MAX_BUFFER > sizeof(SDSSIMAGE));
  assert(MAX_GALAXY_NAME > sizeof(MAX_OBJID_SIZE+2));
#ifdef DEBUG_TYPE
  printf("ERROR: DEBUG_TYPE is set\n");
#endif

  memset(pGalaxyCommon,0,sizeof(GALAXYCOMMON));
  /* Loop through the arguments */
  output_name[0] = 0;
  sdss_name[0] = 0;

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: stray argument %s\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'v':
        case 'V':
          verbose = 1;
          break;



        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(output_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 'b': /* binary file name */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sdss_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;



        default:
          printf("* illegal command -%c-",cmdchar);
          errorFlag = 1;
          break;
        }
        
      }
    }
  }

  /* Attempt to open the list */
  if (output_name[0] == 0) {
    printf("ERROR: Output file name not specified\n");
    errorFlag = 1;
  } else {
    output_handle = fopen(output_name,"wt");
    if (output_handle == NULL) {
      errorFlag = 1;
      printf("Could not open output file %s\n",output_name);
    }
  }




  sdss_input_handle = Open(sdss_input_name,"r");
  if (sdss_input_handle == NULL) {
    printf("ERROR: Failed to find Sdss file %s\n",sdss_input_name);
    errorFlag = 1;
  } else {
    if (verbose) {
      printf("Found the Sdss file file %s\n",sdss_input_name);
    }
  }

  strcpy(sdss_binary_name,sdss_name);
  slashPtr = strrchr(sdss_binary_name,'/');
  if (slashPtr == NULL) {
    sdss_binary_name[0] = 0;
  } else {
    slashPtr++;
    *slashPtr = 0;
  }
  strcpy(sdss_index_name,sdss_binary_name);
  strcat(sdss_binary_name,"sdss.dat");
  strcat(sdss_index_name,"sdss.idx");


  if (errorFlag) {
    printf("Usage: formatapass -o {output name} \n");
		printf("       where -b <filename> is the binary file name \n");
    printf("       where -v is the verbose flag\n");
    return(-1);
  }

  printf("formatapass of %s %s Output Filename %s Binary Filename %s Sorted Binary Filename %s Index filename %s\n",
         __DATE__,__TIME__,output_name,sdss_name,sdss_binary_name,sdss_index_name);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of SDSSIMAGE is %d Size of GALAXYREC is %d \n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(SDSSIMAGE),sizeof(GALAXYREC));
 

  time(&startTime);
  
#ifdef DEBUG_TYPE
  fprintf(output_handle,"catalogname\tra\tdec\tsequence\tcatalogmag\tradius\n");
#else /* DEBUG_TYPE */
  fprintf(output_handle,"catalogname\tra\tdec\ttype\tcatalogmag\tradius\n");
#endif /* DEBUG_TYPE */
  fprintf(output_handle,"-----------\t--\t---\t----\t----------\t------\n");

		

  /* Now read in the Sdss catalog */
  sdss_header = table_header(sdss_input_handle,TABLE_PARSE);
  if (sdss_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",sdss_input_name);
    
    exit(-1);
  }
  sdss_table = table_loadva(sdss_input_handle,
                            &sdss_header,
                            NULL, /* hbase */
                            NULL, /* rows */
                            NULL,
                            sizeof(SDSSIMAGE),
                            &sdss_nrecs,
                            TblBuf,"objID"       ,TblOff(PSDSSIMAGE,objID),MAX_OBJID_SIZE,
                            TblDbl,"ra"          ,TblOff(PSDSSIMAGE,ra),
                            TblDbl,"dec"         ,TblOff(PSDSSIMAGE,dec),
                            TblDbl,"u"           ,TblOff(PSDSSIMAGE,u),
                            TblDbl,"g"           ,TblOff(PSDSSIMAGE,g),
                            TblDbl,"r"           ,TblOff(PSDSSIMAGE,r),
                            TblDbl,"i"           ,TblOff(PSDSSIMAGE,i),
                            TblDbl,"z"           ,TblOff(PSDSSIMAGE,z),
                            TblDbl,"extinction_u",TblOff(PSDSSIMAGE,extinction_u),
                            TblDbl,"extinction_g",TblOff(PSDSSIMAGE,extinction_g),
                            TblDbl,"extinction_r",TblOff(PSDSSIMAGE,extinction_r),
                            TblDbl,"extinction_i",TblOff(PSDSSIMAGE,extinction_i),
                            TblDbl,"extinction_z",TblOff(PSDSSIMAGE,extinction_z),
                            TblInt,"type"        ,TblOff(PSDSSIMAGE,type),
                            TblInt,"status"      ,TblOff(PSDSSIMAGE,status),
                            TblInt,"flags"       ,TblOff(PSDSSIMAGE,flags),
                            TblInt,"flags2"      ,TblOff(PSDSSIMAGE,flags2),
                            0,"end",0);
  if (sdss_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",sdss_input_name);
    exit(-1);
  }
	time(&curTime);
	curTime -= startTime;
  printf("Read %lld records from %s at %d seconds\n",sdss_nrecs,sdss_input_name,curTime);
  totalCount = sdss_nrecs;
  if (sdss_name[0] != 0) {
    galaxy_handle = Open(sdss_name,"w");
    if (galaxy_handle == NULL) {
      printf("Could not open binary file %s\n",sdss_name);
      exit(-1);
    } else {
      galaxy_table = (PGALAXYREC)calloc(totalCount,sizeof(GALAXYREC));
      if (galaxy_table == NULL) {
        printf("ERROR: failed to allocate binary table of size %d\n",totalCount);
        exit(-1);
      }
    }
    galaxy_binary_handle = Open(sdss_binary_name,"w");
    if (galaxy_binary_handle == NULL) {
      printf("Could not open sorted binary file %s\n",sdss_binary_name);
      exit(-1);
    }    
    galaxy_index_handle = Open(sdss_index_name,"w");
    if (galaxy_index_handle == NULL) {
      printf("Could not open index file %s\n",sdss_index_name);
      exit(-1);
    }
  }
	time(&curTime);
	curTime -= startTime;
  printf("Allocated %d records %d seconds\n",totalCount,curTime);


  catalogdir = getenv("DASCH_CATALOG");
  if (catalogdir == NULL) {
    fprintf(stderr,"DASCH_CATALOG is not defined\n");
    printf("0\n");
    return(-1);
  }

  strcpy(galaxy_name,catalogdir);
  slashPtr = strrchr(galaxy_name,'/');
  if (slashPtr != NULL) {
    slashPtr++;
  } else {
    slashPtr = galaxy_name;
  }
  *slashPtr = 0;
  strcat(galaxy_name,"galaxy.dat");

  if (PopulateGalaxyTree(pGalaxyCommon,galaxy_name,0) != 0) {
    exit(-1);
  }
	time(&curTime);
	curTime -= startTime;
  printf("Read the galaxy common at %d seconds\n",curTime);

  for (sdss_index = 0; sdss_index < sdss_nrecs; sdss_index++) {
    pSdss = &sdss_table[sdss_index];
    /* Accept only stars */
    if (pSdss->type != SDSS_TYPE_STAR) {
      nonStarCount++;
      continue;
    }
    /* Reject any galaxies or NGC objects that are too close */
    FindNearestGalaxy(pGalaxyCommon,pSdss->ra,pSdss->dec,pNullCatalogStar,nearbyObjects,pGalaxyResult);
    if (pGalaxyResult->galaxyFoundFlag != 0) {
      nearbyGalaxyCount++;
      continue;
    }
    sprintf(galaxyType,"SDSSg%.1fc%.1f",pSdss->g-pSdss->extinction_g,pSdss->g-pSdss->extinction_g-(pSdss->r-pSdss->extinction_r));
    if (strlen(galaxyType) > (MAX_GALAXY_TYPE-2)) {
      printf("ERROR: MAX_GALAXY_TYPE exceeded for %s\n",galaxyType);
      exit(-1);
    }

    fprintf(output_handle,"%s\t%.5f\t%.5f\t%s\t%.2f\t%f\n",pSdss->objID,pSdss->ra,pSdss->dec,galaxyType,pSdss->g-pSdss->extinction_g,GALAXY_DEFAULT_RADIUS);
    if (galaxy_table != NULL) {
      pGalaxyRec = &galaxy_table[finalCount];
      pGalaxyRec->galaxyflag       = GALAXY_FLAG;
      pGalaxyRec->galaxyversion    = GALAXY_VERSION;
      pGalaxyRec->ra         = pSdss->ra;
      pGalaxyRec->dec        = pSdss->dec;
      pGalaxyRec->catalogmag = pSdss->g-pSdss->extinction_g;
      pGalaxyRec->radius     = GALAXY_DEFAULT_RADIUS;
      strcpy(pGalaxyRec->catalogname,pSdss->objID);
      strcpy(pGalaxyRec->galaxytype,galaxyType);
      
    }
    finalCount++;

  }
  if (sdss_table != NULL) {
    Free(sdss_table);
    sdss_table = NULL;
  }

  if (sdss_header != NULL) {
    table_hdrfree(sdss_header);
    sdss_header = NULL;
  }

  if (sdss_input_handle != NULL) {
    Close(sdss_input_handle);
    sdss_input_handle = NULL;
  }
 
  time(&curTime);
  curTime -= startTime;
  printf("Prepared %d records in %d seconds\n",
         finalCount,
         curTime);
  

	time(&curTime);
	curTime -= startTime;
  printf("Writing output table at %d seconds\n",curTime);



  if (galaxy_table != NULL) {
    writeItems = Write(galaxy_handle,galaxy_table,sizeof(GALAXYREC),finalCount);
    if (writeItems != finalCount) {
      printf("ERROR writing %s items %d\n",sdss_name,writeItems);
      exit(-1);
    } else {
      printf("Wrote %d items of size %d to %s\n",writeItems,sizeof(GALAXYREC),sdss_name);
    }
    Close(galaxy_handle);
    galaxy_handle = Open(sdss_name,"r");
    if (galaxy_handle == NULL) {
      printf("Could not open binary file %s for reading\n",sdss_name);
      exit(-1);
    }
  
    time(&curTime);
    curTime -= startTime;
    printf("Allocating sort index at %d seconds\n",curTime);


    sort_index_table = (PSORTINDEX)calloc(finalCount*MAX_ADJACENT_BINS,sizeof(SORTINDEX));
    if (sort_index_table == NULL) {
      printf("ERROR: failed to allocate the sort index table of size %d\n",totalCount*MAX_ADJACENT_BINS);
      exit(-1);
    }
    time(&curTime);
    curTime -= startTime;
    printf("Preparing sort index at %d seconds\n",curTime);
    for (curRecord = 0; curRecord < finalCount; curRecord++) {
      pGalaxyRec = &galaxy_table[curRecord];
      gsc_bin_index = GetGSCBin(pGscBin,pGalaxyRec->ra,pGalaxyRec->dec,&decBin,&raBin,"formatsdss");
      FindAdjacentBins(pGscBin,gsc_bin_index,bin_list,&bin_count);
      bin_list[bin_count] = gsc_bin_index;
      bin_count++;
      for (bin_index = 0; bin_index < bin_count; bin_index++) {
        pSortIndex = &sort_index_table[sortIndexSize];
        pSortIndex->gsc_bin_index = bin_list[bin_index];
        pSortIndex->record_number = curRecord;
        sortIndexSize++;
        if (sortIndexSize >= finalCount*MAX_ADJACENT_BINS) {
          printf("ERROR: sortIndexSize %d exceeds %d\n",sortIndexSize,finalCount*MAX_ADJACENT_BINS);
          exit(-1);
        }
      }
    }
    time(&curTime);
    curTime -= startTime;
    printf("Sorting %d items at %d seconds\n",sortIndexSize,curTime);
    qsort((void *)sort_index_table,sortIndexSize,sizeof(SORTINDEX),IndexCompare);
    time(&curTime);
    curTime -= startTime;
    printf("Writing %d items at %d seconds\n",sortIndexSize,curTime);
  
    curOffset = 0;
    curSortIndex = 0;
    pSortIndex = &sort_index_table[curSortIndex];
    gsc_bin_index = pSortIndex->gsc_bin_index;
    pGalaxyRec = &galaxyRec;
    for (cur_bin_index = 0; cur_bin_index < pGscBin->total_gsc_bins; cur_bin_index++) {
      pSdssIndex->offset = curOffset;
      pSdssIndex->binNumber = cur_bin_index;
      pSdssIndex->numStars = 0;
      while (1) {
        if ((cur_bin_index < gsc_bin_index) ||
            (cur_bin_index >= pGscBin->total_gsc_bins)) {
          break;
        }
        if (cur_bin_index == gsc_bin_index) {
          Seek(galaxy_handle,pSortIndex->record_number*sizeof(GALAXYREC),SEEK_SET);
          readItems = Read(galaxy_handle,pGalaxyRec,sizeof(GALAXYREC),1);
          if (readItems != 1) {
            printf("ERROR: reading the binary data file\n");
            exit(-1);
          }
          writeItems = Write(galaxy_binary_handle,pGalaxyRec,sizeof(GALAXYREC),1);
          if (writeItems != 1) {
            printf("ERROR: writing binary file\n");
            exit(-1);
          }
          curOffset += sizeof(GALAXYREC);
          pSdssIndex->numStars++;
          curSortIndex++;
          if (curSortIndex < sortIndexSize) {
            pSortIndex = &sort_index_table[curSortIndex];
            gsc_bin_index =  pSortIndex->gsc_bin_index;
          } else {
            gsc_bin_index =  pGscBin->total_gsc_bins;
            curRecord = finalCount;
            break;
          }
        }
      }
      writeItems = Write(galaxy_index_handle,pSdssIndex,sizeof(STARINDEX),1);
      if (writeItems != 1) {
        printf("ERROR writing index file\n");
        exit(-1);
      }
    }
  }
  free(galaxy_table);
  free(sort_index_table);
	fclose(output_handle);
  if (galaxy_handle != NULL) {
    Close(galaxy_handle);
    Close(galaxy_binary_handle);
    Close(galaxy_index_handle);
  }
	time(&curTime);
	curTime -= startTime;



	time(&curTime);
	curTime -= startTime;
  printf("nonStarCount %d nearbyGalaxyCount %d\n",nonStarCount,nearbyGalaxyCount);
	printf("Execution Time: %d seconds; total entries %d stars written: %d index entries written %d\n",curTime,totalCount,finalCount,sortIndexSize);

	return(EXIT_SUCCESS);
}
