// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatcatalog.c
 *
 *  Creates starcatalog.txt for insertion into the starcatalog table of the MySQL database.  This table converts catalog names that do not have recognizable positions into RA and Dec positions.
 *
 *  gcc -g -O0 -D_FILE_OFFSET_BITS=64 formatcatalog.c  -o formatcatalog -I/usr/include/mysql -I /dasch/install/include -L /dasch/install/lib -L/usr/lib${lib64}/mysql pipelineutils.a -lwcs -lmysqlclient -lm 
 * 
 *  formatcatalog -v -o /dasch/Pipeline/catalogs/starcatalog.txt
 *
 *  see /home/scanner/backup/(*)/formatcatalog*.bat and /home/scanner/backup/(*)/formatcatalog*.log
 * 
 * Mar 14, 2017 Edward J. Los - Initial version
 * Mar 17, 2017 Edward J. Los - Add subtable to remove duplicates
 * Mar 20, 2017 Edward J. Los - Add gsc2.3.2 duplicates: N0200023379 and N020211312460
 * Mar 24, 2017 Edward J. Los - Add gsc_bin_index
 * May 29, 2018 Edward J. Los - Support GAIA
 * Oct 23, 2018 Edward J. Los   Add REF_TYPE_ATLAS2 for the Atlas 2 catalog
 */ 
  
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_INPUT_RECORDS 1000

typedef struct _catalogentry {
  int refType;
  char *catalogName;
  size_t totalRecords;
  size_t acceptedRecords;
  size_t rejectedRecords;
} CATALOGENTRY,*PCATALOGENTRY;

CATALOGENTRY catalogTable[] = {
  {REF_TYPE_GSC,"gsc232bin.dat",0,0}, /* GSC2.3.2 must come before APASS for subTable to handle GSC duplicates in APASS properly */
  {REF_TYPE_TYCHO2,"apass.dat",0,0},
  {REF_TYPE_KEPLER,"kepler.dat",0,0},
  {REF_TYPE_GAIA2,"gaiadr2.dat",0,0},
  {REF_TYPE_ATLAS2,"atlas.dat",0,0},
};
int catalogTableSize = sizeof(catalogTable)/sizeof(CATALOGENTRY);

typedef struct _subtable {
  long long starid;
  int matches; /* Number of matches found */
  int refType; /* Active refType */
} SUBTABLE,*PSUBTABLE;

SUBTABLE subTable[] = {
{5358902821L,0,REF_TYPE_TYCHO2},
{5358902841L,0,REF_TYPE_TYCHO2},
{5358902851L,0,REF_TYPE_TYCHO2},
{5358903531L,0,REF_TYPE_TYCHO2},
{5358903931L,0,REF_TYPE_TYCHO2},
{5358906421L,0,REF_TYPE_TYCHO2},
{5358909461L,0,REF_TYPE_TYCHO2},
{5358909701L,0,REF_TYPE_TYCHO2},
{5358916891L,0,REF_TYPE_TYCHO2},
{5358917731L,0,REF_TYPE_TYCHO2},
{5358917841L,0,REF_TYPE_TYCHO2},
{5358919001L,0,REF_TYPE_TYCHO2},
{5358921131L,0,REF_TYPE_TYCHO2},
{5570701241L,0,REF_TYPE_TYCHO2},
{5570701471L,0,REF_TYPE_TYCHO2},
{5570701731L,0,REF_TYPE_TYCHO2},
{5570702261L,0,REF_TYPE_TYCHO2},
{5570703291L,0,REF_TYPE_TYCHO2},
{5790001841L,0,REF_TYPE_TYCHO2},
{5790003151L,0,REF_TYPE_TYCHO2},
{5790003691L,0,REF_TYPE_TYCHO2},
{5790004131L,0,REF_TYPE_TYCHO2},
{5790006621L,0,REF_TYPE_TYCHO2},
{5790007401L,0,REF_TYPE_TYCHO2},
{5790008941L,0,REF_TYPE_TYCHO2},
{5790009381L,0,REF_TYPE_TYCHO2},
{5790012381L,0,REF_TYPE_TYCHO2},
{5790013121L,0,REF_TYPE_TYCHO2},
{5790015431L,0,REF_TYPE_TYCHO2},
{5790016221L,0,REF_TYPE_TYCHO2},
{110200023379L,0,REF_TYPE_GSC},
{11020211312460L,0,REF_TYPE_GSC},
};

int numSubTable = sizeof(subTable)/sizeof(SUBTABLE);



extern GSCBIN gscBin64;  
PGSCBIN pGscBin = &gscBin64; 

int main(int argc,char *argv[])
{
  char *argstr;
  char output_name[MAX_INPUT_NAME];
  char binary_name[MAX_INPUT_NAME];
  
  FILE* output_handle = NULL;
  FILE* binary_handle = NULL;
  char* slashPtr;
  int outputCount = 0;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  char *charPtr;
  char cmdchar;
  int verbose = 0;
  int catalogIndex;
  PCATALOGENTRY pCatalogEntry;
  PGSCIMAGE pGscTable = NULL;
  PGSCIMAGE pGscEntry;
  size_t inputRecords;
  size_t inputIndex;
  size_t totalRecords = 0;
  size_t acceptedRecords = 0;
  size_t rejectedRecords = 0;
  char REF[MAX_REF];
  long long REFNumber;
  int refType;
  int gsc_bin_index;
  int minDecBin;
  int raBin;

  pGscTable = calloc( MAX_INPUT_RECORDS,sizeof(GSCIMAGE));
  if (pGscTable == NULL) {
    printf("ERROR: failed to allocate pGscTable of size %d\n",MAX_INPUT_RECORDS);
    exit(-1);
  }


  /* Loop through the arguments */
  output_name[0] = 0;


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

  if (errorFlag) {
    printf("Usage: formatcatalog -o {output name} [=v]\n");
    printf("       where -v is the verbose flag\n");
    return(-1);
  }

  printf("formatcatalog of %s %s Output Filename %s\n",
         __DATE__,__TIME__,output_name);
  printf("Size of GSCIMAGE is %d\n",sizeof(GSCIMAGE));
 

  time(&startTime);
#if 0
  fprintf(output_handle,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tgsc_bin_index\n");
  fprintf(output_handle,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t-------------\n");
#endif
  for (catalogIndex = 0; catalogIndex < catalogTableSize; catalogIndex++) {
    pCatalogEntry = &catalogTable[catalogIndex];
    strcpy(binary_name,output_name);
    slashPtr = strrchr(binary_name,'/');
    if (slashPtr == NULL) {
      binary_name[0] = 0;
    } else {
      slashPtr++;
      *slashPtr = 0;
    }
    strcat(binary_name,pCatalogEntry->catalogName);

    binary_handle = fopen(binary_name,"rb");
    if (binary_handle == NULL) {
      printf("Could not open sorted binary file %s\n",binary_name);
      exit(-1);
    } else {
      printf("Opened binary file %s\n",binary_name);
    }

    while (1) {
      inputRecords = fread(pGscTable,sizeof(GSCIMAGE),MAX_INPUT_RECORDS,binary_handle);
      if (inputRecords <= 0) {
        break;
      }
      pCatalogEntry->totalRecords += inputRecords;
      for (inputIndex = 0; inputIndex < inputRecords; inputIndex++) {
        pGscEntry = &pGscTable[inputIndex];
#if 1
        if (pGscEntry->REFNumber == 59502000411L) {
          printf("At REF %lld\n",pGscEntry->REFNumber);
        }
#endif
        GetREF(pGscEntry->REFNumber,(char *)&REF,1,1);
        GetREFNumber(REF,&REFNumber,&refType,1,1);
        /* Look for duplicates in the catalog */
        int subTableIndex;
        for (subTableIndex = 0; subTableIndex < numSubTable; subTableIndex++) {
          if  (subTable[subTableIndex].starid == pGscEntry->REFNumber) {
            subTable[subTableIndex].matches++;
            REFNumber = pGscEntry->REFNumber;
            if (subTable[subTableIndex].refType == REF_TYPE_TYCHO2) {
              if (subTable[subTableIndex].matches == 2) {
                pGscEntry->REFNumber += 3; /* Bump the reference number to avoid a dup */
                break;
              } else if (subTable[subTableIndex].matches > 2) {
                printf("ERROR: more than 1 match for %lld\n", pGscEntry->REFNumber);
                exit(-1);
              }
            } else if (subTable[subTableIndex].refType == REF_TYPE_GSC) {
              if (subTable[subTableIndex].matches > 1) {
                /* Here we use the formate 114ssmm where mm is the number of matches and ss is the sub table index */
                if ((subTable[subTableIndex].matches >= 100) ||
                    (subTableIndex >= 100)) {
                  printf("ERROR matches %d or index %d greater than 100\n",subTable[subTableIndex].matches,subTableIndex);
                  exit(-1);
                }
                pGscEntry->REFNumber = 1140000 + (100*subTableIndex) + subTable[subTableIndex].matches;
              }
            }
            printf("WARNING: subTable matches %2d for %lld -> %lld\n",subTable[subTableIndex].matches,REFNumber,pGscEntry->REFNumber);
          }
        }

        totalRecords++;
        if (refType == pCatalogEntry->refType) {
          pCatalogEntry->acceptedRecords++;
          acceptedRecords++;
          gsc_bin_index = GetGSCBin(pGscBin,pGscEntry->ra,pGscEntry->dec,&minDecBin,&raBin,"formatcatalog");
          fprintf(output_handle,"%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%d\n",
                  pGscEntry->REFNumber,
                  pGscEntry->ra,
                  pGscEntry->dec,
                  pGscEntry->Stdmag,
                  pGscEntry->color,
                  pGscEntry->class,
                  pGscEntry->VFlag,
                  pGscEntry->MAGFlag,
                  pGscEntry->RaPM,
                  pGscEntry->DecPM,
                  gsc_bin_index);
        } else {
          pCatalogEntry->rejectedRecords++;
          rejectedRecords++;
        }
        if ((totalRecords % 1000000) == 0) {
          time(&curTime);
          curTime -= startTime;
          printf("At record %9lld in %lld seconds\n",totalRecords,curTime);
        }
      }
    }
    time(&curTime);
    curTime -= startTime;
    printf("File %s input records %lld accepted records %lld rejected records %lld at %lld seconds\n",pCatalogEntry->catalogName,pCatalogEntry->totalRecords,pCatalogEntry->acceptedRecords,pCatalogEntry->rejectedRecords,curTime);
    fclose(binary_handle);
  }
  time(&curTime);
  curTime -= startTime;
  fclose(output_handle);
  printf("Total records %lld accepted records %lld rejected records %lld Execution Time: %d seconds\n",totalRecords,acceptedRecords,rejectedRecords,curTime);

  return(0);
}
