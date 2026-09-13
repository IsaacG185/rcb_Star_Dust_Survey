// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* comparefits.c
 *
 *  Compare the bodies of two FITS files
 *
 *  Usage comparefits <filename1> <filename2>
 *
 *  gcc -g -O0  comparefits.c  -o comparefits -I /dasch/install/include -L /dasch/install/lib -lcfitsio -lm
 *
 *  2009-04-27 Edward J. Los: Correct for files > 2GB
 *  2011-09-17 Edward J. Los: Add length checks
 *  2015-09-15 Edward J. Los: Add daschunistd.h to test for conflicts with unistd.h
 *  2019-06-11 Edward J. Los: Correct the error printout for a corrupted file
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include "daschunistd.h" /* Test for conflicts */
#include "fitsio.h"
#include "longnam.h"

#define BUFFER_SIZE 2880
int findend(FILE *handle,char *filename)
{
  int counter = 1;
  char buffer[BUFFER_SIZE];
  size_t size;
  while (1) {
    size = fread(buffer,1,BUFFER_SIZE,handle);
    if (size != BUFFER_SIZE) {
      printf("ERROR: comparefits Illegal buffer size while looking for END in %s\n",filename);
      exit(-1);
    }
    buffer[BUFFER_SIZE-1] = 0;
    if (strstr(buffer,"END             ") != 0) {
#if 0
      printf("END found in buffer %d\n",counter);
#endif
      break;
    }
    counter++;
  }

  return(0);
}

int main(int argc,char *argv[])
{
  FILE *handle1;
  FILE *handle2;
  size_t size1;
  size_t size2;
  char buffer1[BUFFER_SIZE];
  char buffer2[BUFFER_SIZE];
  char *buf1Ptr;
  char *buf2Ptr;
  int bufferNum = 0;
  int index;
  off_t filesize1;
  off_t filesize2;
  off_t totalBytes;
  off_t headerSize;
  time_t startTime;
  time_t curTime;
  fitsfile *fptr = NULL;
  int iomode = READONLY;
  char err_text[FLEN_ERRMSG];
  int status = 0;
  long naxes[2];
  int nfound = 0;
  int fpackFlag = 0;
  struct stat filestats;
  int statResult;
  int hdutyp;
  char datasum[FLEN_VALUE];
  char comm[FLEN_COMMENT];
  unsigned long dsum1;
  unsigned long dsum2;
  unsigned long pcount1;
  unsigned long pcount2;
  unsigned long heapbytes1;
  unsigned long heapbytes2;
  double tdouble;
  
  time(&startTime);
  if (argc < 3) {
    printf("ERROR: comparefits Usage: comparefits <filename1> <filename2>\n");
    exit(-1);
  }
 
  if ((strstr(argv[1],".fit.fz") != NULL) &&
      (strstr(argv[2],".fit.fz") != NULL)) {
    fpackFlag = 1;
  }


  statResult = stat(argv[1],&filestats);
  if (statResult != 0) {
    printf("ERROR: comparefits can not stat %s\n",argv[1]);
    exit(-1);
  }
  


  /* Find out how big the first file should be */
  fits_open_file(&fptr,argv[1],iomode,&status);
  if (status != 0) {
    printf("ERROR: comparefits can not read FITS file %s\n",argv[1]);
    fits_get_errstatus(status,err_text);
    printf("CFITSIO ERROR %d: %s\n",status,err_text);
    exit(-1);
  } 
  if (fpackFlag) {
    fits_movrel_hdu(fptr,1,&hdutyp,&status);
    if (status != 0) {
      printf("ERROR: comparefits can not find a second header in  %s\n",argv[1]);
      fits_get_errstatus(status,err_text);
      printf("CFITSIO ERROR %d: %s\n",status,err_text);
      exit(-1);
    }
    if (hdutyp != IMAGE_HDU) {
      printf("ERROR: comparefits found the wrong header type %d in %s\n",hdutyp,argv[1]);
      exit(-1);
    }
    /* Look for the datasum key */
    fits_read_key_str(fptr,"DATASUM",datasum,comm,&status);
    if (status != 0) {
      printf("ERROR: comparefits can not find the datasum keyword  %s\n",argv[1]);
      fits_get_errstatus(status,err_text);
      printf("CFITSIO ERROR %d: %s\n",status,err_text);
      exit(-1);
    }
    tdouble = atof(datasum); /* read as a double as a workaround */
    dsum1 = (unsigned long) tdouble;

    /* Look for the pcount key */
    fits_read_key_str(fptr,"PCOUNT",datasum,comm,&status);
    if (status != 0) {
      printf("ERROR: comparefits can not find the pcount keyword  %s\n",argv[1]);
      fits_get_errstatus(status,err_text);
      printf("CFITSIO ERROR %d: %s\n",status,err_text);
      exit(-1);
    }
    tdouble = atof(datasum); /* read as a double as a workaround */
    pcount1 = (unsigned long) tdouble;
    

  }
    


  fits_read_keys_lng(fptr,"NAXIS",1,2,naxes,&nfound,&status);
  if (status != 0) {
    printf("ERROR: comparefits Can not read NAXIS keywords in  %s\n",argv[1]);
    fits_get_errstatus(status,err_text);
    printf("CFITSIO ERROR %d: %s\n",status,err_text);
    exit(-1);
  }
  if (nfound != 2) {
    printf("ERROR: comparefits Found only %d axes in  %s\n",nfound,argv[1]);
    exit(-1);
  }
  filesize1 = naxes[0];
  filesize1 = filesize1*naxes[1]*2;
#if 0
  printf("File %s has %d bytes of data\n",argv[1],filesize1);
#endif
  status = 0;
  fits_close_file(fptr,&status);
  if (fpackFlag) {
    filesize1 = (((filesize1/2)+BUFFER_SIZE-1)/BUFFER_SIZE) * BUFFER_SIZE;
    heapbytes1 = ((pcount1+BUFFER_SIZE-1)/BUFFER_SIZE) * BUFFER_SIZE;
    filesize1 += heapbytes1;
  }


  headerSize = filestats.st_size - filesize1;
  if (headerSize < 2880) {
    printf("ERROR: comparefits header of %s has only %lldd bytes\n",argv[1],headerSize);
    exit(-1);
  }


  status = 0;

  statResult = stat(argv[1],&filestats);
  if (statResult != 0) {
    printf("ERROR: comparefits can not stat %s\n",argv[1]);
    exit(-1);
  }
  /* Find out how big the second file should be */
  fits_open_file(&fptr,argv[2],iomode,&status);
  if (status != 0) {
    printf("ERROR: comparefits Can not read FITS file %s\n",argv[2]);
    fits_get_errstatus(status,err_text);
    printf("CFITSIO ERROR %d: %s\n",status,err_text);
    exit(-1);
  } 
  if (fpackFlag) {
    fits_movrel_hdu(fptr,1,&hdutyp,&status);
    if (status != 0) {
      printf("ERROR: comparefits can not find a second header in  %s\n",argv[2]);
      fits_get_errstatus(status,err_text);
      printf("CFITSIO ERROR %d: %s\n",status,err_text);
      exit(-1);
    }
    if (hdutyp != IMAGE_HDU) {
      printf("ERROR: comparefits found the wrong header type %d in %s\n",hdutyp,argv[2]);
      exit(-1);
    }
    /* Look for the datasum key */
    fits_read_key_str(fptr,"DATASUM",datasum,comm,&status);
    if (status != 0) {
      printf("ERROR: comparefits can not find the datasum keyword  %s\n",argv[1]);
      fits_get_errstatus(status,err_text);
      printf("CFITSIO ERROR %d: %s\n",status,err_text);
      exit(-1);
    }
    tdouble = atof(datasum); /* read as a double as a workaround */
    dsum2 = (unsigned long) tdouble;


    /* Look for the pcount key */
    fits_read_key_str(fptr,"PCOUNT",datasum,comm,&status);
    if (status != 0) {
      printf("ERROR: comparefits can not find the pcount keyword  %s\n",argv[2]);
      fits_get_errstatus(status,err_text);
      printf("CFITSIO ERROR %d: %s\n",status,err_text);
      exit(-1);
    }
    tdouble = atof(datasum); /* read as a double as a workaround */
    pcount2 = (unsigned long) tdouble;
    

    if (dsum2 != dsum1) {
      printf("ERROR: comparefits data checksums %u and %u do not agree for %s and %s\n",dsum1,dsum2,argv[1],argv[2]);
      exit(-1);
    }

    if (pcount2 != pcount1) {
      printf("ERROR: comparefits heap sizes (pcount) %u and %u do not agree for %s and %s\n",pcount1,pcount2,argv[1],argv[2]);
      exit(-1);
    }


  }





  fits_read_keys_lng(fptr,"NAXIS",1,2,naxes,&nfound,&status);
  if (status != 0) {
    printf("ERROR: comparefits Can not read NAXIS keywords in  %s\n",argv[2]);
    fits_get_errstatus(status,err_text);
    printf("CFITSIO ERROR %d: %s\n",status,err_text);
    exit(-1);
  }
  if (nfound != 2) {
    printf("ERROR: comparefits Found only %d axes in  %s\n",nfound,argv[2]);
    exit(-1);
  }
  filesize2 = naxes[0];
  filesize2 = filesize2*naxes[1]*2;

  if (fpackFlag) {
    filesize2 = (((filesize2/2)+BUFFER_SIZE-1)/BUFFER_SIZE) * BUFFER_SIZE;
    heapbytes2 = ((pcount2+BUFFER_SIZE-1)/BUFFER_SIZE) * BUFFER_SIZE;
    filesize2 += heapbytes2;
  }
#if 0
  printf("File %s has %d bytes of data\n",argv[2],filesize2);
#endif
  status = 0;
  fits_close_file(fptr,&status);
  status = 0;

  headerSize = filestats.st_size - filesize2;
  if (headerSize < 2880) {
    printf("ERROR: comparefits header of %s has only %lld bytes\n",argv[1],headerSize);
    exit(-1);
  }

  if (filesize1 != filesize2) {
    printf("ERROR: comparefits file sizes %lld and %lld do not agree for %s and %s\n",filesize1,filesize2,argv[1],argv[2]);
    exit(-1);
  }


  handle1 = fopen(argv[1],"rb");
  if (handle1 == NULL) {
    printf("ERROR: comparefits Failed to open %s\n",argv[1]);
    exit(-1);
  }
  handle2 = fopen(argv[2],"rb");
  if (handle2 == NULL) {
    printf("ERROR: comparefits Failed to open %s\n",argv[2]);
    exit(-1);
  }
  findend(handle1,argv[1]);
  findend(handle2,argv[2]);
  if (fpackFlag) {
    findend(handle1,argv[1]);
    findend(handle2,argv[2]);

  }

  while (1) {
    bufferNum++;
    if ((bufferNum % 10000) == 0) {
      totalBytes = BUFFER_SIZE;
      totalBytes = bufferNum * totalBytes;
#if 0
      printf("At buffer %7d bytes %10ld\n",bufferNum,totalBytes);
#endif
    }
    size1 = fread(buffer1,1,BUFFER_SIZE,handle1);
    size2 = fread(buffer2,1,BUFFER_SIZE,handle2);

    if (size1 != size2) {
      printf("comparefits size1 %d does not agree with size2 %d, buffer %d file %s\n",size1,size2,bufferNum,argv[2]);
      break;
    }
    if (size1 == 0) {
      break;
    }
    if (size1 != BUFFER_SIZE) {
      printf("comparefits Buffer size %d not %d at buffer %d file %s\n",size1,BUFFER_SIZE,bufferNum,argv[2]);
    }
    buf1Ptr = buffer1;
    buf2Ptr = buffer2;
    for (index = 0; index < size1; index++) {
      if (*buf1Ptr++ != *buf2Ptr++) {
        buf1Ptr--;
        buf2Ptr--;
        printf("ERROR: comparefits value %d not %d at index %d buffer %d for %s\n",*buf1Ptr,*buf2Ptr,index,bufferNum,argv[2]);
	exit(-1);
      }
    }

  }
  totalBytes = BUFFER_SIZE;
  totalBytes = bufferNum * totalBytes;
  if (totalBytes < filesize1) {
    printf("ERROR: comparefits found only %d bytes while expecting %lld bytes in %s\n",totalBytes,filesize1,argv[2]);
    exit(-1);
  }

 
  printf("Successfully compared %7d buffers,  for %s %10lld bytes\n",bufferNum,argv[2],totalBytes);
  fclose(handle1);
  fclose(handle2);
  time(&curTime);
  curTime -= startTime;
#if 0
  printf("Execution Time: %d seconds for %s\n",curTime,argv[2]);
#endif
  exit(0);
}
