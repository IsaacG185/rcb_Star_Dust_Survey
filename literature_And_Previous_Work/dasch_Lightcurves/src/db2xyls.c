// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* db2xyls.c
 *
 * Convert a starbase file to an xyls file for use by astrometry.net 
 * 
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm  db2xyls.c pipelineutils.a -ltable -lutil  -lcfitsio -o db2xyls 
 * 
 *  db2xyls -v -w 1376 -h 1088  -i /dasch/junk/test/los.db -o /dasch/junk/test/los.xyls
 *
 *  Apr 21, 2008 Edward J. Los - Initial version
 *  May 10, 2008 Edward J. Los - Retain original column labels
 *  May 12, 2008 Edward J. Los - Filter stars near plate edges
 *  May 13, 2008 Edward J. Los - Add experimental fields and experimental output file
 *                               Filter out sextractor errors
 *                               Add KRON_RADIUS filtering (which has a random effect)
 *                               Add ELLIPTICITY filtering
 *  Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                               exceeds 2GB
 *  Apr 21, 2009 Edward J. Los - Correct error messages
 *  Mar  5, 2012 Edward J. Los - Add mask support
 *                               Redefine the binning factor to relate the mask to the input file
 *  Mar 22, 2013 Edward J. Los - Introduce reverse binning to undo the effects of filter_multiple
 */   

#include <math.h>
#include <errno.h>
#include <time.h>
#include "table.h"
#include "fitsio.h"
#include "longnam.h"
#include "pipelineutils.h"

/* #define ALL_FIELDS 1  */

#define MAX_BUFFER 512
#define REGION_FILE 1 /* Output the debug file as a region file */
/* #define LOS_DEBUG 1 */
/* #define REVERSE_KRON 1 */
/* #define USE_KRON_RADIUS 1 */
/* Mask table */
typedef struct _maskdata {
  int NUMBER; /* Sextractor NUMBER */
  int maskIndex; /* Mask index */
  int X_IMAGE; /* X_IMAGE as an integer */
  int Y_IMAGE; /* Y_IMAGE as an integer */
} MASKDATA,*PMASKDATA;

/* reformatted Sextractor results */
typedef struct _inputdata {
  int FLAGS;            /* Sextractor flags field */
  int NUMBER;           /* Sextractor object number */
  double MAG_ISO;        /* Sextractor isophotonic magnitude */
  double X_IMAGE;        /* Sextractor X location in pixels */
  double Y_IMAGE;        /* Sextractor Y location in pixels */
  double ELLIPTICITY;
#ifdef USE_KRON_RADIUS
  double KRON_RADIUS;
#endif /* USE_KRON_RADIUS */
#ifdef ALL_FIELDS
  double FLUX_ISO;
  double BACKGROUND;
  double FLUX_MAX;
  double THETA_J2000;
  double FWHM_IMAGE;
  int ISOAREA_IMAGE;
#endif /* ALL_FIELDS */
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


/* Forward sort on NUMBER */
int ImageNUMBERCompare(const void *first, const void *second) 
{
  int numberFirst = ((PINPUTDATA)first)->NUMBER;
  int numberSecond = ((PINPUTDATA)second)->NUMBER;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }
}
int MaskNUMBERCompare(const void *first, const void *second) 
{
  int numberFirst = ((PMASKDATA)first)->NUMBER;
  int numberSecond = ((PMASKDATA)second)->NUMBER;
  if (numberFirst > numberSecond) {
    return(1);
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
  char *charPtr;
  int nvals;
  char output_name[MAX_BUFFER];
  char debug_name[MAX_BUFFER];
  FILE * debugHandle = NULL;
  int errorFlag = 0;
  int verbose = 0;
  time_t startTime;
  time_t curTime;
  int binning = 1;

  double kronRadiusLimit = 0.0;
  double ellipticityLimit = 1.0;
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
  double fraction = 0.0;
  int firstrow;
  int filterSextractorErrors = 0;
  int minimumMaskIndex = -1;
  int maskRejected = 0;

  File input_handle = NULL;
  char input_name[MAX_BUFFER];
  TableHead input_header = NULL;
  PINPUTDATA input_table = NULL;
  size_t input_nrecs = 0;
  int input_index;
  int input_index2;
  PINPUTDATA pInput = NULL;
  PINPUTDATA pInput2 = NULL;
  int cur_out = 0;

  int bfields = 3;
  char bextname[] = "SOURCES";             /* extension name */
  int columnIndex;
  int gotFLAGSColumn = 0;
  int gotBFLAGSColumn  = 0;
  int gotNUMBERColumn = 0;

  char *btype[] = { "X_IMAGE","Y_IMAGE","MAG_ISO"};
  char *bform[] = { "D"     ,"D"       ,"D"          };
  char *bunit[] = { "pixel"    ,"pixel","mag"  };

  File mask_handle = NULL;
  char mask_name[MAX_BUFFER];
  TableHead mask_header = NULL;
  PMASKDATA mask_table = NULL;
  size_t mask_nrecs = 0;
  int mask_index;
  PMASKDATA pMask = NULL;


  input_name[0] = 0;
  output_name[0] = 0;
  debug_name[0] = 0;
  mask_name[0] = 0;

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

        case 'e':
        case 'E':
          filterSextractorErrors = 1;
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

        case 'd': /* debug file name */
        case 'D':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(debug_name,*++argv,MAX_BUFFER-2);
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
        case 'r':  
        case 'R':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&binning);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the number of bins %s\n",*argv);
              errorFlag = 1;
            } else {
              if (binning <= 0) {
                fprintf(stderr,"ERROR: binning %d must be a positive integer\n",binning);
                errorFlag = 1;
              }
            }
          }
          break;

        case 'm':  /* Minimum mask index */  
        case 'M':
          argc--;
          if (argc < 2) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&minimumMaskIndex);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the number of bins %s\n",*argv);
              errorFlag = 1;
            } else { 
              charPtr = argv[1];
              if (*charPtr == '-') {
                printf("ERROR: Second argument to -m must be a filename but is %s\n",charPtr);
                errorFlag = 1;
              } else if (strlen(charPtr) > (MAX_BUFFER-2)) {
                printf("ERROR: Second argument to -m is too long for a filename\n");
                errorFlag = 1;
              } else {
                argc--;
                strcpy(mask_name,*++argv);
              }
            }

          }
              
          
          break;

#ifdef USE_KRON_RADIUS
        case 'k': /* Kron Radius limit */
        case 'K':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&kronRadiusLimit);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the Kron Radius Limit %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;
#endif /* USE_KRON_RADIUS */

        case 'l': /* Ellipticity limit */
        case 'L':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&ellipticityLimit);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the Ellipticity Limit %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        
        case 'f': /* Border fraction */
        case 'F':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&fraction);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the border fraction %s\n",*argv);
              errorFlag = 1;
            } else {
              if ((fraction < 0.0) || (fraction >= 0.5)) {
                fprintf(stderr,"ERROR: border fraction %f must be between 0 and 0.5 \n",fraction);
                errorFlag = 1;
              }
            }
          }
          break;

        
        case 'n': /* Number to output */
        case 'N':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&num_out);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the number to output %s\n",*argv);
              errorFlag = 1;
            } else {
              if (num_out <= 0) {
                fprintf(stderr,"ERROR: num_out %d must be a positive integer\n",num_out);
                errorFlag = 1;
              }
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
            nvals = sscanf(*++argv,"%d",&width);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the mosaic width %s\n",*argv);
              errorFlag = 1;
            } else {
              if (width <= 0) {
                fprintf(stderr,"ERROR: width %d must be a positive integer\n",width);
                errorFlag = 1;
              }
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
            nvals = sscanf(*++argv,"%d",&height);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode the mosaic height %s\n",*argv);
              errorFlag = 1;
            } else {
              if (height <= 0) {
                fprintf(stderr,"ERROR: height %d must be a positive integer\n",height);
                errorFlag = 1;
              }
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
  if (width == 0) {
    fprintf(stderr,"ERROR: No width specified\n");
    errorFlag = 1;
  }
  if (height == 0) {
    fprintf(stderr,"ERROR: No height specified\n");
    errorFlag = 1;
  }
  /* Now verify that we have all of the files */
  input_handle = Open(input_name,"rt");
  if (input_handle == NULL) {
    fprintf(stderr,"ERROR: Could not open file %s\n",input_name);
    errorFlag = 1;
  }
  remove(output_name); // Remove the old file if it exists
  fits_create_file(&fptr,output_name,&status);
  if (status != 0) {
    printf("Failed to create %s\n",output_name);
    fits_report_error(stderr,status);
    errorFlag = 1;
  }

  if (strlen(debug_name) > 0) {
    debugHandle = fopen(debug_name,"wt");
    if (debugHandle == NULL) {
      printf("Failed to open the debug file %s\n",debug_name);
      errorFlag = 1;

    }
  }



  if (errorFlag) {
    fprintf(stderr,"Usage: db2xyls -i <input name> \n");
    fprintf(stderr,"                  -o <output name> \n");
    fprintf(stderr,"                  -n <number to write> \n");
    
    fprintf(stderr,"                 [-r <reverse binning factor>] \n");
    fprintf(stderr,"                 [-f <fraction>}\n");
#ifdef USE_KRON_RADIUS
    fprintf(stderr,"                 [-k <Kron radius limit>}\n");
#endif /* USE_KRON_RADIUS */
    fprintf(stderr,"                 [-l <ellipticity limit>}\n");
    fprintf(stderr,"                  -m <mask number> <mask filename> enable mask support\n");
    fprintf(stderr,"                  -v verbose\n");
    fprintf(stderr,"                  -e Remove sextractor errors (FLAGS > 7)\n");
    fprintf(stderr,"                  -d <debug file for plotting>\n");

    
    fprintf(stderr,"NOTE: the fraction determines a square central region\n");

    return(-1);
  }

  if (verbose) {
    printf("db2xyls of %s %s \n Input Filename %s\n Output Filename %s\n width %d height %d number to write %d border %f Kron %f ellipticity %f\n",
           __DATE__,__TIME__,input_name,output_name,width,height,num_out,fraction,kronRadiusLimit,ellipticityLimit);
  }
#ifdef REVERSE_KRON

  printf("ERROR: KRON RADIUS CHECK IS REVERSED\n");
#endif /* REVERSE_KRON */

  if (verbose) {
    printf("Size of INPUTDATA %d\n",sizeof(INPUTDATA));
  }

  time(&startTime);

  /* Read in the input file */

  input_header = table_header(input_handle,TABLE_PARSE);
  if (input_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",input_name);
    return(-1);
  }

  /* See if FLAGS or BFLAGS is present */
  for (columnIndex = 1; columnIndex <= input_header->header->ncol; columnIndex++) {
    if (strcmp(input_header->header->column[columnIndex],"FLAGS") == 0) {
      gotFLAGSColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"BFLAGS") == 0) {
      gotBFLAGSColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"NUMBER") == 0) {
      gotNUMBERColumn = 1;
    }
  }
  if ((gotNUMBERColumn == 0) && 
      ((minimumMaskIndex >= 0) || 
       (mask_name[0] != 0))) {
    printf("ERROR: an input file with a NUMBER column is needed when processing masks %s\n",input_name);
    return(-1);
  }

  if (gotFLAGSColumn) {

    if (gotNUMBERColumn) {
      input_table = table_loadva(input_handle,
                                 &input_header,
                                 NULL, /* hbase */
                                 NULL, /* rows */
                                 NULL,
                                 sizeof(INPUTDATA),
                                 &input_nrecs,
                                 TblInt,"FLAGS",TblOff(PINPUTDATA,FLAGS),
                                 TblInt,"NUMBER",TblOff(PINPUTDATA,NUMBER),
                                 TblDbl,"X_IMAGE",TblOff(PINPUTDATA,X_IMAGE),
                                 TblDbl,"Y_IMAGE",TblOff(PINPUTDATA,Y_IMAGE),
                                 TblDbl,"MAG_ISO",TblOff(PINPUTDATA,MAG_ISO),
#ifdef USE_KRON_RADIUS
                                 TblDbl,"KRON_RADIUS",TblOff(PINPUTDATA,KRON_RADIUS),
#endif /* USE_KRON_RADIUS */
                                 TblDbl,"ELLIPTICITY",TblOff(PINPUTDATA,ELLIPTICITY),
#ifdef ALL_FIELDS
                                 TblDbl,"FLUX_ISO",TblOff(PINPUTDATA,FLUX_ISO),
                                 TblDbl,"BACKGROUND",TblOff(PINPUTDATA,BACKGROUND),
                                 TblDbl,"FLUX_MAX",TblOff(PINPUTDATA,FLUX_MAX),
                                 TblDbl,"THETA_J2000",TblOff(PINPUTDATA,THETA_J2000),
                                 TblDbl,"FWHM_IMAGE",TblOff(PINPUTDATA,FWHM_IMAGE),
                                 TblInt,"ISOAREA_IMAGE",TblOff(PINPUTDATA,ISOAREA_IMAGE),

#endif /* ALL_FIELDS */
                                 0,"end",0);
    } else {

      input_table = table_loadva(input_handle,
                                 &input_header,
                                 NULL, /* hbase */
                                 NULL, /* rows */
                                 NULL,
                                 sizeof(INPUTDATA),
                                 &input_nrecs,
                                 TblInt,"FLAGS",TblOff(PINPUTDATA,FLAGS),
                                 TblDbl,"X_IMAGE",TblOff(PINPUTDATA,X_IMAGE),
                                 TblDbl,"Y_IMAGE",TblOff(PINPUTDATA,Y_IMAGE),
                                 TblDbl,"MAG_ISO",TblOff(PINPUTDATA,MAG_ISO),
#ifdef USE_KRON_RADIUS
                                 TblDbl,"KRON_RADIUS",TblOff(PINPUTDATA,KRON_RADIUS),
#endif /* USE_KRON_RADIUS */
                                 TblDbl,"ELLIPTICITY",TblOff(PINPUTDATA,ELLIPTICITY),
#ifdef ALL_FIELDS
                                 TblDbl,"FLUX_ISO",TblOff(PINPUTDATA,FLUX_ISO),
                                 TblDbl,"BACKGROUND",TblOff(PINPUTDATA,BACKGROUND),
                                 TblDbl,"FLUX_MAX",TblOff(PINPUTDATA,FLUX_MAX),
                                 TblDbl,"THETA_J2000",TblOff(PINPUTDATA,THETA_J2000),
                                 TblDbl,"FWHM_IMAGE",TblOff(PINPUTDATA,FWHM_IMAGE),
                                 TblInt,"ISOAREA_IMAGE",TblOff(PINPUTDATA,ISOAREA_IMAGE),

#endif /* ALL_FIELDS */
                                 0,"end",0);

    }

  } else if (gotBFLAGSColumn) {


    input_table = table_loadva(input_handle,
                               &input_header,
                               NULL, /* hbase */
                               NULL, /* rows */
                               NULL,
                               sizeof(INPUTDATA),
                               &input_nrecs,
                               TblDbl,"X_IMAGE",TblOff(PINPUTDATA,X_IMAGE),
                               TblDbl,"Y_IMAGE",TblOff(PINPUTDATA,Y_IMAGE),
                               TblDbl,"MAG_ISO",TblOff(PINPUTDATA,MAG_ISO),
                               TblInt,"BFLAGS",TblOff(PINPUTDATA,FLAGS),
                               TblInt,"NUMBER",TblOff(PINPUTDATA,NUMBER),
#ifdef USE_KRON_RADIUS
                               TblDbl,"KRON_RADIUS",TblOff(PINPUTDATA,KRON_RADIUS),
#endif /* USE_KRON_RADIUS */
                               TblDbl,"ELLIPTICITY",TblOff(PINPUTDATA,ELLIPTICITY),
#ifdef ALL_FIELDS
                               TblDbl,"FLUX_ISO",TblOff(PINPUTDATA,FLUX_ISO),
                               TblDbl,"BACKGROUND",TblOff(PINPUTDATA,BACKGROUND),
                               TblDbl,"FLUX_MAX",TblOff(PINPUTDATA,FLUX_MAX),
                               TblDbl,"THETA_J2000",TblOff(PINPUTDATA,THETA_J2000),
                               TblDbl,"FWHM_IMAGE",TblOff(PINPUTDATA,FWHM_IMAGE),
                               TblInt,"ISOAREA_IMAGE",TblOff(PINPUTDATA,ISOAREA_IMAGE),

#endif /* ALL_FIELDS */
                               0,"end",0);

    

  } else {
    printf("ERROR: No FLAGS nor BFLAGS in %s\n",input_name);
    return(-1);
  }

  if (input_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",input_name);
    return(-1);
  }
  if (verbose) {
    fprintf(stderr,"read %d records for %s\n",input_nrecs,input_name);
  }
  if (gotBFLAGSColumn) {
    /* If we read in full Sextractor file, keep only the original Sextractor flags */
    for (input_index = 0; input_index < input_nrecs; input_index++) {
      pInput = &input_table[input_index];
      pInput->FLAGS = pInput->FLAGS & FILTER_BMASK_SEXTRACTOR;
#if 0
      printf("INPUT: %d\n",pInput->NUMBER);
#endif
    }
  }
  if (binning != 1) {
    for (input_index = 0; input_index < input_nrecs; input_index++) {
      pInput = &input_table[input_index];
      pInput->X_IMAGE  = pInput->X_IMAGE * binning;
      pInput->Y_IMAGE  = pInput->Y_IMAGE * binning;
    }
  }

  if ((minimumMaskIndex >= 0) && (mask_name[0] != 0)) {
    /* We have a mask table.  Read it in and use it to filter the Sextractor table */
    mask_handle = Open(mask_name,"rt");
    if (mask_handle == NULL) {
      fprintf(stderr,"ERROR: Could not open file %s\n",mask_name);
      errorFlag = 1;
    }
    mask_header = table_header(mask_handle,TABLE_PARSE);
    if (mask_header == NULL) {
      fprintf(stderr,"ERROR: Failed to read header for %s\n",mask_name);
      return(-1);
    }

    mask_table = table_loadva(mask_handle,
                              &mask_header,
                              NULL, /* hbase */
                              NULL, /* rows */
                              NULL,
                              sizeof(MASKDATA),
                              &mask_nrecs,
                              TblInt,"X_IMAGE",TblOff(PMASKDATA,X_IMAGE),
                              TblInt,"Y_IMAGE",TblOff(PMASKDATA,Y_IMAGE),
                              TblInt,"maskIndex",TblOff(PMASKDATA,maskIndex),
                              TblInt,"NUMBER",TblOff(PMASKDATA,NUMBER),
                              0,"end",0);

    if (mask_table == NULL) {
      fprintf(stderr,"ERROR: Failed to read table for %s\n",mask_name);
      return(-1);
    }
    if (verbose) {
      fprintf(stderr,"read %d records for %s\n",mask_nrecs,mask_name);
    }
    /* Sort both the Sextractor table and the mask table in ascending Sextractor NUMBER */
    qsort((void*)input_table,input_nrecs,sizeof(INPUTDATA),ImageNUMBERCompare);
    qsort((void*)mask_table,mask_nrecs,sizeof(MASKDATA),MaskNUMBERCompare);
    
    mask_index = 0;
    pMask = mask_table;
    input_index2 = 0;
    pInput2 = input_table;
    for (input_index = 0; input_index < input_nrecs; input_index++) {
      pInput = &input_table[input_index];
      if ((pInput->NUMBER < pMask->NUMBER) ||
          (mask_index == mask_nrecs)) {
        if (input_index2 != input_index) {
          pInput2 = &input_table[input_index2];
          memcpy(pInput2,pInput,sizeof(INPUTDATA));
        }
        input_index2++;
        continue;
      }
      if (pInput->NUMBER == pMask->NUMBER) {
        if (pMask->maskIndex < minimumMaskIndex) {
          if (input_index2 != input_index) {
            pInput2 = &input_table[input_index2];
            memcpy(pInput2,pInput,sizeof(INPUTDATA));
          }
          input_index2++;
        } else {
#if 0
          printf("MASKED: %d %d %d %d\n",pMask->NUMBER,pMask->maskIndex,pMask->X_IMAGE,pMask->Y_IMAGE);
#endif
        }
    
      } 

      while ((pInput->NUMBER >= pMask->NUMBER) && (mask_index < mask_nrecs)) {
        pMask++;
        mask_index++;
      }
    }
    maskRejected = input_nrecs - input_index2;
    input_nrecs = input_index2;

  }

#if 0
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
    printf("OUTPUT: %d\n",pInput->NUMBER);
  }
#endif
  /* Reverse Sort the table according to MAG_ISO */
  qsort((void*)input_table,input_nrecs,sizeof(INPUTDATA),ImageCompare);
#ifdef LOS_DEBUG
  for (input_index = 1; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
    pInput2 = &input_table[input_index-1];
    fprintf(stderr,"Index %6d  MAG_ISO %f X_INDEX %f Y_INDEX %f\n",
            input_index,
            pInput->MAG_ISO,
            pInput->X_IMAGE,
            pInput->Y_IMAGE);
    if (pInput->MAG_ISO < pInput2->MAG_ISO) {
      fprintf(stderr,"ERROR: MAG_ISO not properly sorted %d %f %d %f\n",
              input_index,
              pInput->MAG_ISO,
              input_index-1,
              pInput2->MAG_ISO);
      exit(-1);
    }
  }
#endif /* LOS_DEBUG */
  if (fraction > 0.0) {
    double subwidthX = ((1.0*width)/2.0) - (fraction*width);
    double subwidthY = ((1.0*height)/2.0) - (fraction*height);
    double minsubwidth;
    double minXindex;
    double maxXindex;
    double minYindex;
    double maxYindex;

    if (subwidthX < subwidthY) {
      minsubwidth = subwidthX;
    } else {
      minsubwidth = subwidthY;
    }
    minXindex = ((1.0*width)/2) - minsubwidth;
    maxXindex = ((1.0*width)/2) + minsubwidth;
    minYindex = ((1.0*height)/2) - minsubwidth;
    maxYindex = ((1.0*height)/2) + minsubwidth;

    for (input_index = 0; input_index < input_nrecs; input_index++) {
      pInput = &input_table[input_index];
      if ((pInput->X_IMAGE > minXindex) &&
          (pInput->X_IMAGE < maxXindex) &&
          (pInput->Y_IMAGE > minYindex) &&
          (pInput->Y_IMAGE < maxYindex)) {
        if ((filterSextractorErrors) && 
            (pInput->FLAGS > 7)) {
          continue;
        }
#ifdef USE_KRON_RADIUS
#ifdef REVERSE_KRON
        if ((kronRadiusLimit > 0.0) &&
            ((pInput->KRON_RADIUS == 0) ||
             (pInput->KRON_RADIUS < kronRadiusLimit))) {
          continue;
        }

#else /* REVERSE_KRON */
        if ((kronRadiusLimit > 0.0) &&
            (pInput->KRON_RADIUS > kronRadiusLimit)) {
          continue;
        }
#endif /* REVERSE_KRON */
#endif /* USE_KRON_RADIUS */

        if (pInput->ELLIPTICITY > ellipticityLimit) {
          continue;
        }

        if (input_index != cur_out) {
          pInput2 = &input_table[cur_out];
          memcpy(pInput2,pInput,sizeof(INPUTDATA));
        }
        cur_out++;        
      }
    }
    
  } else {
    cur_out = input_nrecs;
  }

  if ((num_out == 0) ||
      (num_out > cur_out)) {
    num_out = cur_out;
  }

  outTable = (double *)calloc(num_out,sizeof(double));
  if (outTable == NULL) {
    fprintf(stderr,"ERROR allocating outTable\n");
    exit(-1);
  }
  /* Write a debug file suitable for plotting */
  if (debugHandle != NULL) {
#ifndef REGION_FILE
    fprintf(debugHandle,"MAG_ISO\tELLIPTICITY\tX_IMAGE\tY_IMAGE\n");
    fprintf(debugHandle,"-------\t-----------\t-------\t-------\n");
#endif /* REGION_FILE */    

    for (input_index = 0; input_index < num_out; input_index++) {
      pInput = &input_table[input_index];
#ifdef REGION_FILE
      fprintf(debugHandle,"CIRCLE(%.0f,%.0f,10) # text = {%.2f %.2f}\n",pInput->X_IMAGE,pInput->Y_IMAGE,pInput->MAG_ISO,pInput->ELLIPTICITY);

#else /* REGION_FILE */
      fprintf(debugHandle,"%f\t%f\t%f\t%f\n",pInput->MAG_ISO,1.0*(pInput->ELLIPTICITY),pInput->X_IMAGE,pInput->Y_IMAGE);
#endif /* REGION_FILE */
    }
  }
    
  /* Now create the first HDU */
  fits_create_img(fptr,bitpix,naxis,NULL,&status);
  if (status != 0) {
    gsc_report_error(__LINE__,fptr,status);
  }


  fits_create_tbl(fptr,BINARY_TBL,num_out,bfields,btype,bform,bunit,bextname,&status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr, status);
    return(-1);
  }

  fits_write_key(fptr,TINT,"IMAGEW",&width,"Image width in pixels",&status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr, status);
    return(-1);
  }
  fits_write_key(fptr,TINT,"IMAGEH",&height,"Image width in pixels",&status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr, status);
    return(-1);
  }
  for (input_index = 0; input_index < num_out; input_index++) {
    pInput = &input_table[input_index];
    outTable[input_index] = pInput->X_IMAGE;
  }
  firstrow = 1;
  fits_write_col(fptr, TDOUBLE, 1, firstrow, 1L, num_out,  outTable,  &status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr, status);
    return(-1);
  }


  for (input_index = 0; input_index < num_out; input_index++) {
    pInput = &input_table[input_index];
    outTable[input_index] = pInput->Y_IMAGE;
  }
  fits_write_col(fptr, TDOUBLE, 2, firstrow, 1L, num_out,outTable,  &status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr, status);
    return(-1);
  }


  for (input_index = 0; input_index < num_out; input_index++) {
    pInput = &input_table[input_index];
    outTable[input_index] = pInput->MAG_ISO;
  }
  fits_write_col(fptr, TDOUBLE, 3, firstrow, 1L, num_out,outTable,  &status);
  if ( status != 0 ) {
    gsc_report_error(__LINE__,fptr, status);
    return(-1);
  }


  time(&curTime);
  curTime -= startTime;



  
  printf("Execution time %d seconds. Total %d, center %d  written %d masked %d eFlag %d border %f Kron %f ellip %f  %s\n",
         curTime,
         input_nrecs,
         cur_out,
         num_out,
         maskRejected,
         filterSextractorErrors,
         fraction,
         kronRadiusLimit,
         ellipticityLimit,
         output_name);

  status = 0;
  fits_close_file(fptr,&status);
  if (outTable != NULL) {
    free(outTable);
  }

  if (input_table != NULL) {
    Free(input_table);
  }
  if (input_header != NULL) {
    table_hdrfree(input_header);
  }
  if (input_handle != NULL) {
    Close(input_handle);
  }
  if (debugHandle != NULL) {
    fclose(debugHandle);
  }

}
