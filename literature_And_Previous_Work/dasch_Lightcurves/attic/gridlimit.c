// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* gridlimit.c
 *
 *  gcc -ggdb   -I /dasch/install/include  -L /dasch/install/lib gridlimit.c -lm  -o gridlimit
 * 
 *  Read in a list of annular9.m calibration grid files and determine the limits of these
 *  files.
 */   


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_INPUTNAME 256
#define MAX_BUFFER 512
#define MAX_FILE 4000  /* 3251 */
#define MAX_MAGNITUDE 24
#define BINS_PER_MAGNITUDE 100
#define MAX_BINS 2400 /* This is the maximum expected file size */

/* #define LOS_DEBUG 1 */  
int ReadGridFile(char *filename,double *inputrefmag,double * inputisomag,int *pGridCount)
{
  char inLine[MAX_BUFFER];
  int lineLen;
  FILE * inputHandle;
  char *inBuffer;
  int nlines = 0;
  int gridcount = 0;
  double refmag;
  double isomag;
  double rms;
  int flag;
  int nvals;
  *pGridCount = 0;

  inputHandle = fopen(filename,"rt");
  if (inputHandle == NULL) {
    printf("Could not open file %s\n",filename);
    return(-1);
  }
    
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,inputHandle);
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
      nvals = sscanf(inBuffer,"%lf %lf %lf %d",&refmag,&isomag,&rms,&flag);
      if (nvals == 4) {
        if (flag == 1) {
          
          if (gridcount < MAX_BINS) {
            inputrefmag[gridcount] = refmag;
            inputisomag[gridcount] = isomag;
            gridcount++;
          } else {
            printf("MAX_BIN exceeded\n");
            return(-1);
          }
        }
        

      } else {
        printf("Format error in line %d of %s %s\n",nlines,filename,inBuffer);
      }
      nlines++;
    }
  } 
  fclose(inputHandle);
#if 0
  printf("Number of lines in file: %d, number of points %d in file %s\n",nlines,gridcount,filename);
#endif
  *pGridCount = gridcount;
  return(0);
}


int ReadInputFile(char *filenames,int *pFileCount,char *listname)
{
  char inLine[MAX_BUFFER];
  int lineLen;
  FILE * inputHandle;
  char *inBuffer;
  int nlines = 0;
  int fileCount = 0;
  *pFileCount = 0;

  inputHandle = fopen(listname,"rt");
  if (inputHandle == NULL) {
    printf("Could not open file %s\n",listname);
    return(-1);
  }
    
  while (1) {
    inBuffer = fgets(inLine,MAX_BUFFER,inputHandle);
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
      if (lineLen >= MAX_INPUTNAME) {
        printf("Line length %d exceeds maximum input name %d\n",lineLen,MAX_INPUTNAME);
      }
      if (fileCount < MAX_FILE) {
        strcpy(&filenames[fileCount*MAX_INPUTNAME],inBuffer);
        fileCount++;
      }
      nlines++;
    }
  } 
  fclose(inputHandle);
  printf("Number of files in the list: %d, number used %d\n",nlines,fileCount);

  *pFileCount = fileCount;
  return(0);
}

/* Sorting routine borrowed from statstable.c by John Roll */
int dcmp(const void *a, const void *b)
{
  double *pa = (double *) a;
  double *pb = (double *) b;

  if ( *pa  > *pb ) return  1;
  if ( *pa == *pb ) return  0;
  if ( *pa  < *pb ) return -1;
}

int main(int argc,char *argv[]) {
  char filenames[MAX_INPUTNAME*MAX_FILE];
  int filecount;
  int fileindex;
  int gridcount;
  int maxgridcount = 0;
  int totgridcount = 0;
  time_t startTime;
  time_t curTime;
  double inputrefmag[MAX_BINS];
  double inputisomag[MAX_BINS];
  double miniso = 500;
  double maxiso = -500;
  double minref = 500.0;
  double maxref = -500.0;
  int index;
  double *refarray;
  int *refsize;
  int bin;
  int binerror1 = 0;
  int binerror2 = 0;
  int binerror3 = 0;
  char *minname;
  FILE *outHandle;
  int maxbinsize = 0;

  time(&startTime);

  refarray = (double *)calloc(MAX_FILE*MAX_BINS,sizeof(double));
  refsize = (int *)calloc(MAX_BINS,sizeof(int));
  if ((refarray == NULL) ||
      (refsize == NULL)) {
    printf("Failed to allocate arrays\n");
    return(-1);
  }
  if (argc < 3) {
    printf("Usage: gridlimit <file list> <output file>\n");
    return(-1);
  }
  if (ReadInputFile(filenames,&filecount,argv[1]) != 0) {
    return(-1);
  }
  outHandle = fopen(argv[2],"wt");
  if (outHandle == NULL) {
    printf("Failed to open %s\n",argv[2]);
    return(-1);
  }

  for (fileindex = 0; fileindex < filecount; fileindex++) {
    if (ReadGridFile(&filenames[fileindex*MAX_INPUTNAME],inputrefmag,inputisomag,&gridcount) != 0) {
      return(-1);
    }
    totgridcount += maxgridcount;
    if (gridcount > maxgridcount) {
      maxgridcount = gridcount;
    }
    for (index = 0; index < gridcount; index++) {
      if (inputisomag[index] < miniso) {
        miniso = inputisomag[index];
      }
      if (inputisomag[index] > maxiso) {
        maxiso = inputisomag[index];
      }
      if (inputrefmag[index] < minref) {
        minref = inputrefmag[index];
        minname = &filenames[fileindex*MAX_INPUTNAME];
      }
      if (inputrefmag[index] > maxref) {
        maxref = inputrefmag[index];
      }
      bin = (inputrefmag[index] * BINS_PER_MAGNITUDE) + 0.005;
      if ((bin >= 0) && (bin < MAX_BINS)) {
        if (refsize[bin] >= MAX_FILE) {
          binerror1++;
        } else {
          if (((MAX_FILE * bin) + refsize[bin]) >= (MAX_FILE*MAX_BINS)) {
            binerror2++;
          } else {
            refarray[(MAX_FILE * bin) + refsize[bin]] = inputisomag[index];
            refsize[bin]++;
          }
        }
      } else {
        binerror3++;
      }


    } /* End of point loop within each file */

  } /* End of file loop */
  for (bin = 0; bin < MAX_BINS; bin++) {
    if (refsize[bin] > 0) {
      
      int binsize = refsize[bin];
      if (binsize > maxbinsize) {
        maxbinsize = binsize;
      }
      qsort((void*)&refarray[MAX_FILE*bin],binsize,sizeof(double),dcmp);
      fprintf(outHandle,"%f %f %f %f %f %f %d\n",
              (1.0*bin)/(1.0 * BINS_PER_MAGNITUDE),
              refarray[(MAX_FILE*bin)],
              refarray[(MAX_FILE*bin)+(binsize/10)],
              refarray[(MAX_FILE*bin)+(binsize/2)],
              refarray[(MAX_FILE*bin)+((9*binsize)/10)],
              refarray[(MAX_FILE*bin)+binsize-1],
              binsize);
    }

  }

  printf("Max bin size %d\n",maxbinsize);
  printf("Min, max refmag: %f %f, min,max isomag: %f %f\n",minref,maxref,miniso,maxiso);
  printf("min refmag name %s\n",minname);
  printf("Total points %d, maximum per file %d\n",totgridcount,maxgridcount);
  printf("Error points %d %d %d\n",binerror1,binerror2,binerror3);

  free(refarray);
  free(refsize);
  fclose(outHandle);
  time(&curTime);
  curTime -= startTime;
  printf("Execution time: %d seconds\n",curTime);
 

}

