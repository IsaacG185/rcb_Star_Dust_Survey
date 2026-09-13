// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* dumpxyls.c
 *
 * convert an xyls file to a starbase db file
 * 
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm  dumpxyls.c pipelineutils.a -ltable -lutil  -lcfitsio -o dumpxyls 
 * 
 *    dumpxyls  -v -i /dasch/junk/test/mc39048_01_16r270ww.xyls -o /dasch/junk/test/los.db
 *
 *  Apr 21, 2008 Edward J. Los - Initial version
 *  May 12, 2008 Edward J. Los - Correct verbose flag
 */   

#include <math.h>
#include <errno.h>
#include <time.h>
#include "table.h"
#include "fitsio.h"
#include "longnam.h"
#include "pipelineutils.h"



#define MAX_BUFFER 512
/* #define LOS_DEBUG 1 */


/* reformatted Sextractor results */
typedef struct _inputdata {
  double MAG_ISO;     /* Sextractor isophotonic magnitude */
  double X_IMAGE;        /* Sextractor X location in pixels */
  double Y_IMAGE;        /* Sextractor Y location in pixels */
} INPUTDATA,*PINPUTDATA;


/* REVERSE Sort routine based on MAG_ISO */
int ImageCompare(const void *first, const void *second) 
{
  double numberFirst = ((PINPUTDATA)first)->MAG_ISO;
  double numberSecond = ((PINPUTDATA)second)->MAG_ISO;
  /* REVERSE COMPARISON */
  if (numberFirst > numberSecond) {
    return(1);
    /* REVERSE_COMPARISON */
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }
}
void gsc_report_error(int line,
                      fitsfile *fptr,
                      int status)
{
  if (fptr == NULL) {
    fprintf(stderr,"Error in line %d. NO FILE\n",
            line);
  } else {
    fprintf(stderr,"Error in line %d, file %s\n",
            line,fptr->Fptr->filename);
  }
  fits_report_error(stderr,status);
  exit(-1);
  return;
}


int main(int argc,char *argv[])
{
  char *argstr;
  char cmdchar;
  int nvals;
  char output_name[MAX_BUFFER];
  FILE * outputHandle;
  int errorFlag = 0;
  int verbose = 0;
  time_t startTime;
  time_t curTime;
  int binning = 1;
  int num_out = 0;
  int width = 0;
  int height = 0;
  fitsfile *fptr = NULL;
  int status = 0;  
  int iomode = READWRITE;
  char err_text[FLEN_ERRMSG];
  char keywordstr[FLEN_KEYWORD];
  int bitpix   =  BYTE_IMG; /* number of bits per data pixel */
  long naxis    =   0; /* 0 - dimensional image                            */ 
  double *outTable = NULL;
  int firstRow;
  int hdutype;
  int inull = -1;
  int anynull = 0;
  
  File input_handle = NULL;
  char input_name[MAX_BUFFER];
  TableHead input_header = NULL;
  PINPUTDATA input_table = NULL;
  long input_nrecs = 0;
  long input_index;
  PINPUTDATA pInput = NULL;
  PINPUTDATA pInput2 = NULL;

  int bfields = 3;
  char bextname[] = "SOURCES";             /* extension name */
  

  char *btype[] = { "X","Y","FLUX"};
  char *bform[] = { "D"     ,"D"       ,"D"          };
  char *bunit[] = { "pixel"    ,"pixel","mag"  };

  input_name[0] = 0;
  output_name[0] = 0;

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      fprintf(stderr,"ERROR: unqualified argument %s\n",argstr);
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
            strncpy(output_name,*++argv,MAX_BUFFER-2);
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
            strncpy(input_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;
        


        default:
          fprintf(stderr,"ERROR: illegal command -%c-",cmdchar);
          errorFlag = 1;
          break;
        }
        
      }
    }
  }

  /* Verify that we have everything */
  if (output_name[0] == 0) {
    fprintf(stderr,"ERROR: No output file specified\n");
    errorFlag = 1;
  }
  if (input_name[0] == 0) {
    fprintf(stderr,"ERROR: No input file specified\n");
    errorFlag = 1;
  }
  /* Now verify that we have all of the files */
  input_handle = Open(input_name,"rt");
  if (input_handle == NULL) {
    fprintf(stderr,"ERROR: Could not open file %s\n",input_name);
    errorFlag = 1;
  }


  /* Open the input file */
  fits_open_file(&fptr, input_name, READONLY, &status);
  if (status != 0) {
    printf("Error Opening %s\n",input_name);
    fits_report_error(stderr, status);
    errorFlag = 1;
  }



  outputHandle = fopen(output_name,"wt");
  if (outputHandle == NULL) {
    fprintf(stderr,"ERROR Failed to open the output file %s\n",output_name);
    errorFlag = 1;
  }
    




  if (errorFlag) {
    fprintf(stderr,"Usage: dumpxyls -i <input name> \n");
    fprintf(stderr,"                -o <output name> \n");
    fprintf(stderr,"                -v verbose \n");
    return(-1);
  }

  if (verbose) {
    printf("dumpxyls of %s %s \n Input Filename %s\n Output Filename %s\n\n",
         __DATE__,__TIME__,input_name,output_name);


    printf("Size of INPUTDATA %d\n",sizeof(INPUTDATA));
  }

  time(&startTime);

  // Move to the second Binary table extension
  fits_movabs_hdu (fptr, 2, &hdutype, &status);
  if ( status != 0 ) {
    fits_report_error(stderr, status);
    return(-1);
  }

  fits_get_num_rows (fptr,&input_nrecs, &status);
  if ( status != 0 ) {
    fits_report_error(stderr, status);
    return(-1);
  }
  outTable = (double *)calloc(input_nrecs,sizeof(double));
  if (outTable == NULL) {
    fprintf(stderr,"ERROR allocating outTable\n");
    exit(-1);
  }
  input_table = (PINPUTDATA)calloc(input_nrecs,sizeof(INPUTDATA));
  if (input_table == NULL) {
    fprintf(stderr,"ERROR allocating input_table\n");
    exit(-1);
  }

  firstRow = 1;
  fits_read_col(fptr, TDOUBLE, 1, firstRow, 1L, input_nrecs, &inull, outTable,
                &anynull, &status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr,status);
    return(-1);
  }

  if (anynull != 0) {
    fprintf(stderr,"ERROR: null found in column 1\n");
    return(-1);
  }
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
    pInput->X_IMAGE = outTable[input_index];
  }

  fits_read_col(fptr, TDOUBLE, 2, firstRow, 1L, input_nrecs, &inull, outTable,
                &anynull, &status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr, status);
    return(-1);
  }

  if (anynull != 0) {
    fprintf(stderr,"ERROR: null found in column 2\n");
    return(-1);
  }
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
    pInput->Y_IMAGE = outTable[input_index];
  }

  fits_read_col(fptr, TDOUBLE, 3, firstRow, 1L, input_nrecs, &inull, outTable,
                &anynull, &status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr,status);
    return(-1);
  }

  if (anynull != 0) {
    fprintf(stderr,"ERROR: null found in column 1\n");
    return(-1);
  }
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
    pInput->MAG_ISO = outTable[input_index];
  }

  fprintf(outputHandle,"X_IMAGE\tY_IMAGE\tMAG_ISO\n");
  fprintf(outputHandle,"-------\t-------\t-------\n");



  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
    fprintf(outputHandle,"%f\t%f\t%f\n",
            pInput->X_IMAGE,
            pInput->Y_IMAGE,
            pInput->MAG_ISO);

  }




  time(&curTime);
  curTime -= startTime;


  if (verbose) {
    printf("Execution time %d seconds. Stars in %d\n",curTime,input_nrecs);
  }

  status = 0;
  fits_close_file(fptr,&status);
  if (outTable != NULL) {
    free(outTable);
  }
  if (input_table != NULL) {
    free(input_table);
  }
  if (outputHandle != NULL) {
    fclose(outputHandle);
  }


}
