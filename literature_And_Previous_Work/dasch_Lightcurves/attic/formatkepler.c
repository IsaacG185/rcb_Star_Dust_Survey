// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatkepler.c
 *
 * Read a starbase-formatted kepler input catalog and put it into the gsc232bin.db format
 *
 *  1.  Use this program to reformat the Kepler Input Catalog
 *      formatkepler -v -i /Kepler/temp.db -o /dasch/Pipeline/catalogs/kepler.db
 *  2.  Sort the results into 0.01 degree bins:
 *      
 *      index -mb -n kepler.db dec
 *
 *  4.  Run formatgsc to put this file in binary format and to index it
 *      
 *      formatgsc  /dasch/Pipeline/catalogs/kepler.db
 * 
 *      For statistics only:
 *      echo "formatgsc -s /dasch/Pipeline/catalogs/kepler.db" | at now
 *       
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm formatkepler.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib64/mysql -lmysqlclient -o formatkepler
 * 
 *
 *   Execution Time: 815 seconds for 13161645 records; header 616 noCQCount 0 maxcqLength 5 for ----- maxdeltacolor 0.000000 for K4046
 *  
 *  Source:          SCP   2MASS   PHOTO    TYBV   UNCAL   NOCAL
 *  Count:       2106821  417229 9092682   16414 1490571   37312
 *  No GMAG:      115304  834458       0       0 1495869   74624
 *  No RMAG:           0       0       0       0       0       0
 *  Variable:          0       8     434       0      54       0
 *  Galaxy:            0       0       0       0       0   37312
 *  No PMRA:      178744       0       0       0       0   37312
 *  No PMDEC:     178744       0       0       0       0   37312
 *  Mag   0            0       0       0       0       0       0
 *  Mag   1            0       0       0       0       0       0
 *  Mag   2            0       0       0       1       0       0
 *  Mag   3            0       0       0       1       0       0
 *  Mag   4            0       0       0       6       0       0
 *  Mag   5            0       0       0      29       0       0
 *  Mag   6            7       0       0      51       0       0
 *  Mag   7           93       0       0     137       0       0
 *  Mag   8          327       0       3     292       2       0
 *  Mag   9         1090       0     173     739     155       0
 *  Mag  10         3229       0     694    1843      81       0
 *  Mag  11         8404       0    2437    5117      84       0
 *  Mag  12        21223       0   12199    6753     639       0
 *  Mag  13        49905       0   40516    1377    1717       0
 *  Mag  14       108463       0   98529      57    3554       0
 *  Mag  15       218130       0  210084       1    6902       0
 *  Mag  16       412936       0  421464       0   14204       0
 *  Mag  17       661019       0  831364      10   26880       0
 *  Mag  18       497115       0 1736509       0   45698       0
 *  Mag  19        11553       0 5738710       0  330761       0
 *  
 * Aug 17, 2009 Edward J. Los - Initial version
 */   


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512

#define MAGBINS  19   /* Number of magnitude bins */

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  

typedef struct _keplerimage {
  int kepler_id;        /* Kepler ID */
  double ra;           /* Right Ascension in degrees */
  double dec;          /* Declination in degrees */
  double gmag;         /* green magnitude */
  double rmag;         /* reg magnitude */
  double color;        /* color */
  double pmra;        /* ra proper motion */
  double pmdec;           /* dec proper motion*/
  int cqtype;          /* cq variable type */
  char galaxy;         /* galaxy flag*/
  char variable;       /* variable flag */
  char cq[MAX_REF];   /* quality string */
  
} KEPLERIMAGE,*PKEPLERIMAGE;


int ReadInput(PKEPLERIMAGE pInput,
              File input_handle,
              TableHead input_header,
              TblDescriptor input_descriptor,
              TableRow* input_row)
{
  *input_row = table_rowget(input_handle,input_header,*input_row,NULL,NULL,0);
  if (*input_row == NULL) {
    return(0);
  }
  if (!table_loadrow(input_handle,input_header,*input_row,input_descriptor,(char *)pInput)) {
    printf("ERROR: Read Input table_loadrow failed\n");
    return(0);
  }
  return(1);
}



int main(int argc,char *argv[])
{
  char *argstr;
  char output_name[MAX_INPUT_NAME];
  File output_handle = NULL;
  int outputCount = 0;
  int zeromagCount = 0;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int linecounter = 0;
  int index;
  char *charPtr;
  char cmdchar;
  int verbose = 0;
  char *dotPtr = NULL;
  int compareResult;
  int binNumber;
  int oldBinNumber = -1;
  int dec_bin;
  int dec_bin2;
  int oldDecBin = 0;
  int savedDecBin = -2;
  int ra_bin;
  int *decBinCounts = NULL;
  int maxDecBin = 0;
  double maxRaError = 0;
  double maxDecError = 0;
  int maxWriteCount = 0;
  int totalBinCount = 0;

  File input_handle = NULL;
  char input_name[MAX_INPUT_NAME];
  TableHead input_header = NULL;
  int input_nrecs = 0;
  int input_index;
  KEPLERIMAGE input_record;
  PKEPLERIMAGE pInput = &input_record;
  KEPLERIMAGE old_input_record;
  PKEPLERIMAGE pOldInput = &old_input_record;
  GSCIMAGE gsc_record;
  PGSCIMAGE pGscImage = &gsc_record;

  TblDescriptor input_descriptor = NULL;
  TableRow input_row = NULL;
  PKEPLERIMAGE pSortInput = NULL;
  int sortArrayAlloc = 0;
  int sortArraySize = 0;
 
  double tmpRaPM;
  double tmpDecPM;
  double tmpPM;
  double maxRaPM = 0;
  double maxDecPM = 0;
  double maxPM = 0;
  int countPM = 0;
  int tempmag;

  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int keplerCQCount[KEPLER_CQ_MAGFLAG_MAX];
  int keplerNoGMAG[KEPLER_CQ_MAGFLAG_MAX];
  int keplerNoRMAG[KEPLER_CQ_MAGFLAG_MAX];
  int keplerVariable[KEPLER_CQ_MAGFLAG_MAX];
  int keplerGalaxy[KEPLER_CQ_MAGFLAG_MAX];
  int keplerNoPMRA[KEPLER_CQ_MAGFLAG_MAX];
  int keplerNoPMDEC[KEPLER_CQ_MAGFLAG_MAX];
  int magcount[KEPLER_CQ_MAGFLAG_MAX][MAGBINS+1];
  int keplerNoCQCount = 0;
  int keplerHeaderCount = 0;
  int maxcqLength = 0;
  char maxcq[MAX_REF];
  double deltacolor;
  double maxdeltacolor = 0.0;
  int maxkeplercolorid;
	char REF[MAX_REF];
	int refType;

  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    keplerCQCount[index] = 0;
    keplerNoGMAG[index]= 0;
    keplerNoRMAG[index]= 0;
    keplerVariable[index]= 0;
    keplerGalaxy[index]= 0;
    keplerNoPMRA[index]= 0;
    keplerNoPMDEC[index]= 0;
  }


  decBinCounts = (int *)calloc(pGscBin->dec_bins,sizeof(int));

  /* Loop through the arguments */
  input_name[0] = 0;
  output_name[0] = 0;


  memset(magcount,0,sizeof(magcount));
  

  maxcq[0] = 0;
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
        case 'i': /* input file name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(input_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
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
  if (input_name[0] == 0) {
    printf("ERROR: Input file name not specified\n");
    errorFlag = 1;
  } else {

    input_handle = Open(input_name,"r");
    if (input_handle == NULL) {
      printf("Could not open input file %s\n",input_name);
      errorFlag = 1;
    } 
  }
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
    printf("Usage: formatkepler <input_name> [-v][-s]\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s provides statistics only and skips writing output files\n");

    return(-1);
  }

  printf("formatkepler of %s %s \n Input Filename %s\n Output Filename %s\n",
         __DATE__,__TIME__,input_name,output_name);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of KEPLERIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(KEPLERIMAGE));
 

  time(&startTime);
  
  
  fprintf(output_handle,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\n");
  fprintf(output_handle,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\n");

  /* Now read in the input file */

  input_header = table_header(input_handle,TABLE_PARSE);
  if (input_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",input_name);
    return(-1);
  }

  input_descriptor = table_create_descrip(&input_nrecs,
                                          TblInt,"KEPLER_ID",TblOff(PKEPLERIMAGE,kepler_id),
                                          TblDbl,"RA"       ,TblOff(PKEPLERIMAGE,ra),
                                          TblDbl,"DEC"      ,TblOff(PKEPLERIMAGE,dec),
                                          TblDbl,"GMAG"     ,TblOff(PKEPLERIMAGE,gmag),
                                          TblDbl,"RMAG"     ,TblOff(PKEPLERIMAGE,rmag),
                                          TblDbl,"GRCOLOR"    ,TblOff(PKEPLERIMAGE,color),
                                          TblDbl,"PMRA"     ,TblOff(PKEPLERIMAGE,pmra),
                                          TblDbl,"PMDEC"    ,TblOff(PKEPLERIMAGE,pmdec),
                                          TblByt,"GALAXY"   ,TblOff(PKEPLERIMAGE,galaxy),
                                          TblByt,"VARIABLE" ,TblOff(PKEPLERIMAGE,variable),
                                          TblBuf,"CQ"       ,TblOff(PKEPLERIMAGE,cq),MAX_REF,
                                          0,"end",0);
  if (input_descriptor == NULL) {
    fprintf(stderr,"ERROR: Failed to allocate descriptor for %s\n",input_name);
    return(-1);
  }
  table_loadmap(input_header,input_descriptor);
  sortArrayAlloc = 1000;
  pSortInput = (PKEPLERIMAGE)calloc(sortArrayAlloc,sizeof(KEPLERIMAGE));

  input_nrecs = 0;
  while (1) {
    memset(pInput,0,sizeof(KEPLERIMAGE));
    memset(pGscImage,0,sizeof(GSCIMAGE));
    if(ReadInput(pInput,input_handle,input_header,input_descriptor,&input_row) == 0) {
      break;
    }
    input_nrecs++;
    if (verbose) {
      if (((input_nrecs+1) % 500000) == 0) {
        time(&curTime);
        curTime -= startTime;
        printf("At record %7d, time %5d sec\n",input_nrecs,curTime);
      }
    }

    if (strlen(pInput->cq) > maxcqLength) {
      maxcqLength = strlen(pInput->cq);
      if (maxcqLength < (MAX_REF-1)) {
        strcpy(maxcq,pInput->cq);
      }
    }
    
    /* Now decode the PQ field to see what we have */
    for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
      if (strstr(pInput->cq,keplerSourceText[index])) {
        break;
      }
    }
    pInput->cqtype = index;
    if (pInput->cqtype < KEPLER_CQ_MAGFLAG_MAX) {
      keplerCQCount[index]++;

      if (pInput->gmag < 0) {
        keplerNoGMAG[index]++;
      } else {
        tempmag = pInput->gmag-0.001;
        if (tempmag < 0) {
          tempmag = 0;
        }
        if (tempmag >= MAGBINS) {
          tempmag = MAGBINS;
        }
        magcount[index][tempmag]++;

      }
      if (pInput->rmag < 0) {
        keplerNoGMAG[index]++;
      }
      if (pInput->variable > 0) {
        keplerVariable[index]++;
      }
      if (pInput->galaxy > 0) {
        keplerGalaxy[index]++;
      }

      if (pInput->pmra < -90.) {
        keplerNoPMRA[index]++;
      }
      if (pInput->pmdec < -90.) {
        keplerNoPMDEC[index]++;
      }
      /* Now fill in our GSC record */
      sprintf(REF,"K%d",pInput->kepler_id); /* A preceeding 'K' flags kepler ids */
			if (GetREFNumber(REF,&pGscImage->REFNumber,&refType,0,0) != 0) {
				printf("ERROR: formatkepler invalid REF %s\n",REF);
				exit(-1);
			}
			pGscImage->ra = pInput->ra * 15.0; /* Convert hours to degrees */
      pGscImage->dec = pInput->dec;
      if (pInput->gmag < -90) {
        pGscImage->Stdmag = 99.0;
      } else {
        pGscImage->Stdmag = pInput->gmag;
      }
      if ((pInput->gmag < -90) || (pInput->rmag < -90)) {
        pGscImage->color = 99.0;
      } else {
        pGscImage->color = pInput->gmag - pInput->rmag;
        deltacolor = pGscImage->color - pInput->color;
        if (deltacolor < 0) {
          deltacolor = -deltacolor;
        }
        if (deltacolor > maxdeltacolor) {
          maxdeltacolor = deltacolor;
          maxkeplercolorid = pInput->kepler_id;
        }
      }
      if (pInput->pmra < -90) {
        pGscImage->RaPM = 0;
      } else {
        pGscImage->RaPM = pInput->pmra * 1000.0;
      }
      if (pInput->pmdec < -90) {
        pGscImage->DecPM = 0;
      } else {
        pGscImage->DecPM = pInput->pmdec * 1000.0;
      }
      if (pInput->galaxy) {
        pGscImage->class = 3;
      }
      if (pInput->variable) {
        pGscImage->VFlag = 1;
      }
      pGscImage->MAGFlag = index | KEPLER_MAGNITUDE_FLAG;
      if (((pInput->gmag > -90.0) &&
           (pInput->gmag <= 19.0)) ||
          ((pInput->rmag > -90.0) &&
           (pInput->rmag <= 19.0))) {

      
        outputCount++;
        fprintf(output_handle,"%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\n",
                pGscImage->REFNumber,
                pGscImage->ra,
                pGscImage->dec,
                pGscImage->Stdmag,
                pGscImage->color,
                pGscImage->class,
                pGscImage->VFlag,
                pGscImage->MAGFlag,
                pGscImage->RaPM,
                pGscImage->DecPM);
      }

    } else {
      if (strstr(pInput->cq,"CQ"))  {
        keplerHeaderCount++;
      } else if (strstr(pInput->cq,"---")) {
        keplerHeaderCount++;
      } else {
        printf("ERROR: Unrecognized CQ flag %s\n",pInput->cq);
        keplerNoCQCount++;
      }
    }
      
  }


  if (input_header != NULL) {
    table_hdrfree(input_header);
  }
  if (input_handle != NULL) {
    Close(input_handle);
  }
  if (input_descriptor != NULL) {
    Free(input_descriptor);
  }
  if (input_row != NULL) {
    table_rowfree(input_row);
  }

  Close(output_handle);
 
  time(&curTime);
  curTime -= startTime;

  printf("Execution Time: %d seconds for %d records; output %d header %d noCQCount %d maxcqLength %d for %s maxdeltacolor %f for K%d\n",curTime,input_nrecs,outputCount,keplerHeaderCount,keplerNoCQCount,maxcqLength,maxcq,maxdeltacolor,maxkeplercolorid);

  printf("Source:     ");
  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    printf("%8s",keplerSourceText[index]);
  }
  printf("\n");

  printf("Count:      ");
  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    printf("%8d",keplerCQCount[index]);
  }
  printf("\n");


  printf("No GMAG:    ");
  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    printf("%8d",keplerNoGMAG[index]);
  }
  printf("\n");

  printf("No RMAG:    ");
  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    printf("%8d",keplerNoRMAG[index]);
  }
  printf("\n");

  printf("Variable:   ");
  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    printf("%8d",keplerVariable[index]);
  }
  printf("\n");

  printf("Galaxy:     ");
  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    printf("%8d",keplerGalaxy[index]);
  }
  printf("\n");

  printf("No PMRA:    ");
  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    printf("%8d",keplerNoPMRA[index]);
  }
  printf("\n");

  printf("No PMDEC:   ");
  for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
    printf("%8d",keplerNoPMDEC[index]);
  }
  printf("\n");
  for (tempmag = 0; tempmag <= MAGBINS; tempmag++) {
    printf("Mag %3d     ",tempmag);
    for (index = 0; index < KEPLER_CQ_MAGFLAG_MAX; index++) {
      printf("%8d",magcount[index][tempmag]);
    }
    printf("\n");
  }    



  if (decBinCounts != NULL) {
    free(decBinCounts);
  }
  if (pSortInput != NULL) {
    free(pSortInput);
  }

  return(EXIT_SUCCESS);
}
