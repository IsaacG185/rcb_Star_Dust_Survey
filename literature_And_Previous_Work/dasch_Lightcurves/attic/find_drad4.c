// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* find_drad4.c
 *
 * Given a list of files, go to the ingest directory and find all the <filename>_dmagcor.grid entries
 * and compile statistics from them.
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio -lmysqlclient  find_drad4.c pipelineutils.a -ltable -lutil  -lwcs -o find_drad4
 *
 * Aug 26, 2008 Edward J. Los - Initial version
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 * Mar 26, 2010 Edward J. Los - Correct the dmagcor name and remove dradAverage
 */


#include <math.h>
#include "table.h"
#include "time.h"
#include "pipelineutils.h"
#include <sys/stat.h>
#define MAX_FILENAME 256
#define MAX_BUFFER 256
#define MAX_LIST_STRING 25

typedef struct _plateEntry {
  char series[MAX_SERIES_STRING];
  char listName[MAX_LIST_STRING];
  int plateNumber;
  int fileExists;
  int fileObsolete;
  int goodBins;
  int goodLocalBins;
  int goodZoutBins;
  int goodMedianBins;
  int goodDradBins;
} PLATEENTRY,*PPLATEENTRY;

typedef struct _dmagcor {
  int nx;             /* Total bins in width */
  int ny;             /* Total bins in height */
  int ix;             /* width bin */
  int iy;             /* height bin */
  double zout;        /* Magnitude correction */
  double errout;      /* RMS of magnitude correction */
  int npout;          /* Number of points used for correction */
  int rejectFlag;     /* Reason, if any, for rejecting this bin */
}  DMAGCOR,*PDMAGCOR;

int main(int argc,char *argv[])
{
  char *argstr;
  char listname[MAX_FILENAME];
  struct stat statbuf;
  FILE * listFile;
  FILE *outHandle = NULL;
  int errorFlag = 0;
  int numPlates = 0;
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
  char cmdchar;
  int plateIndex;
  PPLATEENTRY pPlateEntry;
  PPLATEENTRY pPlateTable = NULL;
  int badBin = 0;
  char outfile[MAX_BUFFER];
  listname[0] = 0;
  outfile[0] = 0;

  File dmagcor_handle = NULL;
  char dmagcor_name[MAX_BUFFER];
  TableHead dmagcor_header = NULL;
  size_t dmagcor_nrecs;
  PDMAGCOR dmagcor_table = NULL;
  PDMAGCOR pDmagcor;
  char *ingestDirectory;
  
  int nx;
  int ny;
  int irec;
  int ix;
  int iy;
  int obsoleteCount = 0;
  int existsCount = 0;


  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      if (strlen(listname) > 0) {
        errorFlag = 1;
        printf("ERROR: list %s is being overwritten by %s\n",listname,argstr);
      } else {
        strncpy(listname,argstr,MAX_FILENAME-1);
      }
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

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Validate arguments */
  if (listname[0] == 0) {
    printf("ERROR: No plate list specified\n");
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

  /* Attempt to open the list */


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

  
  if (errorFlag) {
    printf("Usage: find_drad4 [-v] -o <outname>  <listname>\n");
    return(-1);
  }
  printf("find_drad4 of %s %s, list file %s entries %d, max string %d\n",
         __DATE__,__TIME__,listname,numPlates,maxString);

#ifdef CHILD_DEBUG
  printf("ERROR: CHILD_DEBUG is set\n");
#endif /* CHILD_DEBUG */
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
  for (plateIndex = 0; plateIndex < numPlates; plateIndex++) {
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

    pPlateEntry = &pPlateTable[plateIndex];
    strcpy(pPlateEntry->listName,&listStrings[plateIndex*maxString]);
    if (ParseFilename2(pPlateEntry->listName,pPlateEntry->series,&pPlateEntry->plateNumber) == 0) {
      fprintf(stderr,"ERROR: Failed to parse the name  %s\n",pPlateEntry->listName);
      continue;
    }
    /* Open the local calibration file */
    strcpy(dmagcor_name,ingestDirectory);
    strcat(dmagcor_name,"/");
    strcat(dmagcor_name,&listStrings[plateIndex*maxString]);
    strcat(dmagcor_name,"_dmagcor.db");
   
    dmagcor_handle = Open(dmagcor_name,"r");
    if (dmagcor_handle == NULL) {
      fprintf(stderr,"ERROR: Failed to find the local calibration file %s\n",dmagcor_name);
      continue;
    } else {
      if (verbose) {
        fprintf(stderr,"Found local calibration file %s\n",dmagcor_name);
      }
    }
    
    /* Now read in the local calibration file */

    dmagcor_header = table_header(dmagcor_handle,TABLE_PARSE);
    if (dmagcor_header == NULL) {
      fprintf(stderr,"ERROR: Failed to read header for %s\n",dmagcor_name);
      
      continue;
    }

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

                                 0,"end",0);
    if (dmagcor_table == NULL) {
      fprintf(stderr,"ERROR: Failed to read table for %s\n",dmagcor_name);
      pPlateEntry->fileObsolete = 1;
      obsoleteCount++;
      continue;
    }
    pPlateEntry->fileExists = 1;
    existsCount++;
    if (verbose) {
      fprintf(stderr,"read %d records for %s\n",dmagcor_nrecs,dmagcor_name);
    }

    /* Now validate this table */
    nx = dmagcor_table[0].nx;
    ny = dmagcor_table[0].ny;
    if (dmagcor_nrecs != (nx*ny)) {
      fprintf(stderr,"ERROR: Number of records %d does not match product of %d and %d\n",dmagcor_nrecs,nx,ny);
      return(-1);
    }
    irec = 0;
    for (ix = 0; ix < nx; ix++) {
      for (iy = 0; iy < ny; iy++) {
        pDmagcor = &dmagcor_table[irec];

        if ((pDmagcor->nx != nx) &&
            (pDmagcor->ny != ny) &&
            (pDmagcor->ix != ix) &&
            (pDmagcor->iy != iy)) {
          fprintf(stderr,"ERROR Mismatch of nx %d %d, ny %d %d, ix %d %d, or iy %d %d\n",
                  pDmagcor->nx,nx,
                  pDmagcor->ny,ny,
                  pDmagcor->ix,ix,
                  pDmagcor->iy,iy);
          return(-1);
        }
        badBin = 0;
        if (pDmagcor->errout < 90.0) {
          pPlateEntry->goodLocalBins++;
        } else {
          badBin = 1;
        }
#if 0
        if (pDmagcor->rejectFlag != 0) {
          printf("rejectFlag: %d\n",pDmagcor->rejectFlag);
        }
        if ((ix == 1) && (iy == 46)) {
          printf("At ix  %d iy %d rejectFlag %d\n",ix,iy,pDmagcor->rejectFlag);
        }
#endif
        if ((pDmagcor->rejectFlag & DMAG_REJECT_HIZOUT) == 0) {
          pPlateEntry->goodZoutBins++;
        } else {
          badBin = 1;
        }
        if ((pDmagcor->rejectFlag & DMAG_REJECT_MEDIAN) == 0) {
          pPlateEntry->goodMedianBins++;
        } else {
          badBin = 1;
        }
        if ((pDmagcor->rejectFlag & DMAG_REJECT_DRAD) == 0) {
          pPlateEntry->goodDradBins++;
        } else {
          badBin = 1;
        }
        if (badBin == 0) {
          pPlateEntry->goodBins++;
        } 
        

        irec++;
      }
    }


  }
  fprintf(outHandle,"series\tplateNumber\tfileExists\tfileObsolete\tgoodBins\tgoodLocalBins\tgoodZoutBins\tgoodMedianBins\tgoodDradBins\tlistName\n");
  fprintf(outHandle,"------\t-----------\t----------\t------------\t--------\t-------------\t------------\t--------------\t------------\t--------\n");

  /* Now print everything out */
  for (plateIndex = 0; plateIndex < numPlates; plateIndex++) {
    pPlateEntry = &pPlateTable[plateIndex];
#if 0
    if (pPlateEntry->fileObsolete != 0) {
      continue;
    }
#endif
    fprintf(outHandle,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
            pPlateEntry->series,
            pPlateEntry->plateNumber,
            pPlateEntry->fileExists,
            pPlateEntry->fileObsolete,
            pPlateEntry->goodBins,
            pPlateEntry->goodLocalBins,
            pPlateEntry->goodZoutBins,
            pPlateEntry->goodMedianBins,
            pPlateEntry->goodDradBins,
            pPlateEntry->listName);

  }
  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds for %d mosaics\n",curTime,numPlates);
  printf("Existing Mosaics %d  Obsolete Mosaics %d\n",existsCount,obsoleteCount);
  if (listStrings != NULL) {
    free(listStrings);
  }
  if (pPlateTable != NULL) {
    free(pPlateTable);
  }
  if (dmagcor_handle != NULL) {
    Close(dmagcor_handle);
    dmagcor_handle = NULL;
  }
  if (dmagcor_table != NULL) {
    Free(dmagcor_table);
  }
  if (dmagcor_header != NULL) {
    table_hdrfree(dmagcor_header);
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }


  return(EXIT_SUCCESS);
}
