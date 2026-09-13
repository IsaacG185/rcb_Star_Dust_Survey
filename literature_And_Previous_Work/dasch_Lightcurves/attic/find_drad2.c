// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* find_drad2.c
 *
 *  Given a match_*.db files, plot a map of distortions as a series of fits filed
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm -lcfitsio find_drad2.c pipelineutils.a -ltable -lutil  -o find_drad2
 *
 * /dasch/Pipeline/find_drad2 -w 17402 -h 21053 -i /dasch/junk/find_drad2.db
 *
 * Feb 8, 2008 Edward J. Los - Initial version
 * Feb 16, 2009 Edward J. Los - Use size_t for the number of records in a table to avoid crashes on 64 bit systems when the table size
 *                              exceeds 2GB
 *
 */


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include "fitsio.h"
#include "longnam.h"
#define MAX_BUFFER 100
#define BIN_FACTOR 64

typedef struct _image {
  double X_IMAGE;
  double Y_IMAGE;
  double dra;
  double ddec;
} STARIMAGE,*PSTARIMAGE;

int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char infile[MAX_BUFFER];
  char cmdchar;
  int errorFlag = 0;
  char *ingestDirectory;
  char *binDirectory;
  char outFile[MAX_BUFFER];
  char fileroot[MAX_BUFFER];
  FILE *outHandle = NULL;
  char *charPtr;
  char *basePtr;
  char charVal;
  char err_text[FLEN_ERRMSG];

  File input_handle = NULL;
  char input_name[MAX_BUFFER];
  TableHead input_header = NULL;
  PSTARIMAGE input_table = NULL;
  size_t input_nrecs = 0;
  int input_index;
  PSTARIMAGE pInput;
  int width;
  int height;
  int xBins;
  int yBins;
  int totalBins;
  int usedBins = 0;
  int index;
  int iX;
  int iY;
  double *draValue = NULL;
  double *draStd = NULL;
  double *ddecValue = NULL;
  double *ddecStd = NULL;
  int *binCount = NULL;

  int totalBinCount;
  double totalDraValue;
  double totalDraStd;
  double totalDdecValue;
  double totalDdecStd;

  int maxBinCount = 0;
  double  minDraValue = 0;
  double  maxDraValue = 0;
  double  maxDraStd = 0;
  double  minDdecValue = 0;
  double  maxDdecValue = 0;
  double  maxDdecStd = 0;
  fitsfile *fptr = NULL;
  int status = 0;
  long naxes[2];

  time_t startTime;
  time_t curTime;

  time(&startTime);

  input_name[0] = 0;
  outFile[0] = 0;


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

        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&width);
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
            nvals = sscanf(*++argv,"%d",&height);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;





        default:
          fprintf(stderr,"ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */

  if (width == 0) {
    fprintf(stderr,"ERROR: No mosaic width specified\n");
    errorFlag = 1;
  }
  if (height == 0) {
    fprintf(stderr,"ERROR: No mosaic height specified\n");
    errorFlag = 1;
  }

  if (input_name[0] == 0) {
    fprintf(stderr,"ERROR: No input filename was specified\n");
    errorFlag = 1;
  }




  /* Open the input file */
  input_handle = Open(input_name,"r");
  if (input_handle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to find the input file %s\n",input_name);
  } else {
    fprintf(stderr,"Found input file %s\n",input_name);
  }


  /* Now name the output files */
  strcpy(fileroot,input_name);
  basePtr = fileroot;
  charPtr = strstr(basePtr,"/");
  while (charPtr != NULL) {
    charPtr++;
    basePtr = charPtr;
    charPtr = strstr(basePtr,"/");
  }
  charPtr = strstr(basePtr,".");
  if (charPtr != NULL) {
    *charPtr = 0;
  }


  if (errorFlag) {
    fprintf(stderr,"Usage: find_drad2 -w <mosaic width in pixels> \n");
    fprintf(stderr,"                  -h <mosaic height in pixels> \n");
    fprintf(stderr,"                  -i <input file> \n");

    return(-1);
  }

  /* Now read in the input file */

  input_header = table_header(input_handle,TABLE_PARSE);
  if (input_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",input_name);
    return(-1);
  }

  input_table = table_loadva(input_handle,
                             &input_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(STARIMAGE),
                             &input_nrecs,
                             TblDbl,"X_IMAGE",TblOff(PSTARIMAGE,X_IMAGE),
                             TblDbl,"Y_IMAGE",TblOff(PSTARIMAGE,Y_IMAGE),
                             TblDbl,"dra",TblOff(PSTARIMAGE,dra),
                             TblDbl,"ddec",TblOff(PSTARIMAGE,ddec),
                             0,"end",0);
  if (input_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",input_name);
    return(-1);
  }
  fprintf(stderr,"read %d records for %s\n",input_nrecs,input_name);

  xBins = (width + BIN_FACTOR-1)/BIN_FACTOR;
  yBins = (height + BIN_FACTOR-1)/BIN_FACTOR;
  totalBins = xBins*yBins;
  
  /* Allocate the necessary arrays */
  draValue = calloc(totalBins,sizeof(double));
  draStd = calloc(totalBins,sizeof(double));
  ddecValue = calloc(totalBins,sizeof(double));
  ddecStd = calloc(totalBins,sizeof(double));
  binCount = calloc(totalBins,sizeof(int));

  if ((draValue == NULL) ||
      (draStd == NULL) ||
      (ddecValue == NULL) ||
      (ddecStd == NULL) ||
      (binCount == NULL)) {
    fprintf(stderr,"ERROR: Failed to allocate arrays \n");
    return(-1);
  }
  /* Now go through the input data and populate the arrays */
  for (index = 0; index < input_nrecs; index++) {
    pInput = &input_table[index];
    iX = pInput->X_IMAGE/BIN_FACTOR;
    if (iX < 1) {
      continue;
    }
    if (iX >= (xBins-1)) {
      continue;
    }
    iY = pInput->Y_IMAGE/BIN_FACTOR;
    if (iY < 1) {
      continue;
    }
    if (iY >= (yBins-1)) {
      continue;
    }
    draValue[iX+(iY*xBins)] += pInput->dra;
    draStd[iX+(iY*xBins)] += (pInput->dra)*(pInput->dra);
    ddecValue[iX+(iY*xBins)] += pInput->ddec;
    ddecStd[iX+(iY*xBins)] += (pInput->ddec)*(pInput->ddec);

    binCount[iX+(iY*xBins)]++;


    totalDraValue += pInput->dra;
    totalDraStd += (pInput->dra)*(pInput->dra);
    totalDdecValue += pInput->ddec;
    totalDdecStd += (pInput->ddec)*(pInput->ddec);
    totalBinCount++;

    binCount[iX+(iY*xBins)]++;
    
  }
  /* Now compute averages, maxima and minima */

  for (iX = 1; iX < (xBins-1); iX++) {
    for (iY = 1; iY < (yBins-1); iY++) {
      index = iX + (iY*xBins);
      if (binCount[index] > 0) {
        usedBins++;
        draValue[index] = draValue[index]/binCount[index];
        draStd[index] = draStd[index]/binCount[index] - (draValue[index] * draValue[index]);
        if (draStd[index] < 0.0) {
          draStd[index] = 0.0;
        }
        draStd[index] = sqrt(draStd[index]);
        ddecValue[index] = ddecValue[index]/binCount[index];
        ddecStd[index] = ddecStd[index]/binCount[index] - (ddecValue[index] * ddecValue[index]);
        if (ddecStd[index] < 0.0) {
          ddecStd[index] = 0.0;
        }
        ddecStd[index] = sqrt(ddecStd[index]);
        if (binCount[index] > maxBinCount) {
          maxBinCount = binCount[index];
        }

        if (draValue[index] > maxDraValue) {
          maxDraValue = draValue[index];
        }
        if (draValue[index] < minDraValue) {
          minDraValue = draValue[index];
        }
        if (draStd[index] > maxDraStd) {
          maxDraStd = draStd[index];
        }

        if (ddecValue[index] > maxDdecValue) {
          maxDdecValue = ddecValue[index];
        }
        if (ddecValue[index] < minDdecValue) {
          minDdecValue = ddecValue[index];
        }
        if (ddecStd[index] > maxDdecStd) {
          maxDdecStd = ddecStd[index];
        }

      }
    }
  }
  

  totalDraValue = totalDraValue/totalBinCount;
  totalDraStd = sqrt(totalDraStd/totalBinCount - (totalDraValue * totalDraValue));
  totalDdecValue = totalDdecValue/totalBinCount;
  totalDdecStd = sqrt(totalDdecStd/totalBinCount - (totalDdecValue * totalDdecValue));

  /* Start with the bin count */
  strcpy(outFile,fileroot);
  strcat(outFile,"_binCount.fit");
  remove(outFile);
  fits_create_file(&fptr,outFile,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  naxes[0] = xBins;
  naxes[1] = yBins;
  fits_create_img(fptr,LONG_IMG,2,naxes, &status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create image %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }

  for (iY = 0; iY < yBins; iY++) {
    fits_write_img(fptr,TUINT,1+(xBins * (yBins-1-iY)),xBins,&binCount[iY*xBins],&status);
    if (status != 0) {
      fits_get_errstatus(status,err_text);
      printf("Failed to write file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
      return(-1);
    }
   

  }
  fits_close_file(fptr,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to close file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  printf("Write: %s\n",outFile);

  /* Next output the draValue */
  strcpy(outFile,fileroot);
  strcat(outFile,"_draValue.fit");
  remove(outFile);
  fits_create_file(&fptr,outFile,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  naxes[0] = xBins;
  naxes[1] = yBins;
  fits_create_img(fptr,DOUBLE_IMG,2,naxes, &status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create image %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }

  for (iY = 0; iY < yBins; iY++) {
    fits_write_img(fptr,TDOUBLE,1+(xBins * (yBins-1-iY)),xBins,&draValue[iY*xBins],&status);
    if (status != 0) {
      fits_get_errstatus(status,err_text);
      printf("Failed to write file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
      return(-1);
    }
   

  }
  fits_close_file(fptr,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to close file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  printf("Write: %s\n",outFile);



  /* Next output the draStd */
  strcpy(outFile,fileroot);
  strcat(outFile,"_draStd.fit");
  remove(outFile);
  fits_create_file(&fptr,outFile,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  naxes[0] = xBins;
  naxes[1] = yBins;
  fits_create_img(fptr,DOUBLE_IMG,2,naxes, &status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create image %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }

  for (iY = 0; iY < yBins; iY++) {
    fits_write_img(fptr,TDOUBLE,1+(xBins * (yBins-1-iY)),xBins,&draStd[iY*xBins],&status);
    if (status != 0) {
      fits_get_errstatus(status,err_text);
      printf("Failed to write file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
      return(-1);
    }
   

  }
  fits_close_file(fptr,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to close file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  printf("Write: %s\n",outFile);



  /* Next output the ddecValue */
  strcpy(outFile,fileroot);
  strcat(outFile,"_ddecValue.fit");
  remove(outFile);
  fits_create_file(&fptr,outFile,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  naxes[0] = xBins;
  naxes[1] = yBins;
  fits_create_img(fptr,DOUBLE_IMG,2,naxes, &status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create image %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }

  for (iY = 0; iY < yBins; iY++) {
    fits_write_img(fptr,TDOUBLE,1+(xBins * (yBins-1-iY)),xBins,&ddecValue[iY*xBins],&status);
    if (status != 0) {
      fits_get_errstatus(status,err_text);
      printf("Failed to write file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
      return(-1);
    }
   

  }
  fits_close_file(fptr,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to close file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  printf("Write: %s\n",outFile);



  /* Next output the ddecStd */
  strcpy(outFile,fileroot);
  strcat(outFile,"_ddecStd.fit");
  remove(outFile);
  fits_create_file(&fptr,outFile,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  naxes[0] = xBins;
  naxes[1] = yBins;
  fits_create_img(fptr,DOUBLE_IMG,2,naxes, &status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to create image %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }

  for (iY = 0; iY < yBins; iY++) {
    fits_write_img(fptr,TDOUBLE,1+(xBins * (yBins-1-iY)),xBins,&ddecStd[iY*xBins],&status);
    if (status != 0) {
      fits_get_errstatus(status,err_text);
      printf("Failed to write file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
      return(-1);
    }
   

  }
  fits_close_file(fptr,&status);
  if (status != 0) {
    fits_get_errstatus(status,err_text);
    printf("Failed to close file %s %d %s at line %d\n",outFile,status,err_text,__LINE__);
    return(-1);
  }
  printf("Write: %s\n",outFile);



  /* All done.  Clean up */

  if (draValue != NULL) {
    free(draValue);
  }
  if (draStd != NULL) {
    free(draStd);
  }
  if (ddecValue != NULL) {
    free(ddecValue);
  }
  if (ddecStd != NULL) {
    free(ddecStd);
  }

  if (binCount != NULL) {
    free(binCount);
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

  if (outHandle != NULL) {
    fclose(outHandle);
  }



  time(&curTime);
  curTime -= startTime;
  fprintf(stdout,"maxBinCount %d\n",maxBinCount);
  fprintf(stdout,"minDraValue %f maxDraValue %f maxDraStd %f\n",minDraValue,maxDraValue,maxDraStd); 
  fprintf(stdout,"minDdecValue %f maxDdecValue %f maxDdecStd %f\n",minDdecValue,maxDdecValue,maxDdecStd); 


  fprintf(stdout,"totalBinCount %d\n",totalBinCount);
  fprintf(stdout,"totalDraValue %f totalDraStd %f\n",totalDraValue,totalDraStd); 
  fprintf(stdout,"totalDdecValue %f totalDdecStd %f\n",totalDdecValue,totalDdecStd); 

  fprintf(stdout,"xBins %d yBins %d, totalBins %d usedBins %d,count/usedBin %d\n",xBins,yBins,totalBins,usedBins,totalBinCount/usedBins);

  fprintf(stdout,"maxseconds %d for %s\n",
          curTime,
          input_name);
    


  return(0);
}

