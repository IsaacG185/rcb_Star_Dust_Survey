// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* 
 * compare_colorterm.c - compares old and new colorterm results
 *
 * gcc -ggdb -O0 compare_colorterm.c -o compare_colorterm
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define MAX_PATH 512
#define MAX_BUFFER 512
#define SECOND_COMPARE 1

char *platelist[] =
  {"bm00008_00_01r180ww",
"bm01926_00_01ww",
"bm03206_00_01ww",
"dnr00153_00_01ww",
"dnr00160_00_01ww",
"i14853_00_01ww",
"i31090_00_01r180ww",
"i53176_00_01r180ww",
"mb02645_00_01r180ww",
"mc00342_01_01ww",
"mc12668_01_01r270ww",
"mc21436_04_01r270ww",
"mc21439_03_01r270ww",
"mc39048_01_01r270ww",
"mc39048_02_01r270ww",
"rb14928_00_01ww",
"rh00134_01_01ww",
"rh00796_00_01ww",
"rh03761_00_01ww",
"rh07499_00_01ww",
"rh15577_00_01ww"};

int plateCount = sizeof(platelist)/sizeof(char *);

int main(int argc,char *argv[])
{
  int index;
  int bin;
#ifdef SECOND_COMPARE
  char gscname[] = "/dasch/scanner/linux/gmt/colorterm2.txt";
#else /* SECOND_COMPARE */
  char gscname[] = "/dasch/scanner/linux/gmt/colorterm.txt";
#endif /* SECOND_COMPARE */
  char filename[MAX_PATH];
  FILE *gscHandle;
  FILE *fileHandle;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  int nvals;
  double origColorterm;
  double origErrcolor;
  int origFlag;
  double newColorterm;
  double newErrcolor;
  int newFlag;
  gscHandle = fopen(gscname,"wt");
  if (gscHandle == NULL) {
    printf("ERROR: opening %s\n",gscname);
    exit(-1);
  }

  printf("INFO: Plate count %d\n",plateCount);
#ifdef SECOND_COMPARE
  printf("SECOND_COMPARE is set\n");
#endif /* SECOND_COMPARE */
  for (index = 0; index < plateCount; index++) {
    printf("INFO: Plate %d is %s\n",index,platelist[index]);
    for (bin = 1; bin <= 9; bin++) {
      origColorterm = 99.0;
      origErrcolor = 99.0;
      origFlag = 99;
#ifdef SECOND_COMPARE
      sprintf(filename,"%s%s%s%d%s","/dasch/junk/colorterm/octave/",platelist[index],"_a",bin,".colorterm.txt");
#else /* SECOND_COMPARE */
      sprintf(filename,"%s%s%s%d%s","/dasch/junk/colorterm/colorterm/",platelist[index],"_a",bin,".out.colorterm.txt");
#endif /* SECOND_COMPARE */
      printf("INFO: opening %s\n",filename);
      fileHandle = fopen(filename,"rt");
      if (fileHandle == NULL) {
        printf("ERROR: Could not open file %s\n",filename);
        
      } else {
        inBuffer = fgets(inLine,MAX_BUFFER,fileHandle);
        fclose(fileHandle);
        if (inBuffer == NULL) {
          printf("ERROR: No input buffer for file %s\n",filename);
        } else {
          nvals = sscanf(inBuffer,"%lf %lf %d",&origColorterm,&origErrcolor,&origFlag);
          if (nvals != 3) {
            printf("ERROR: nvals is %d for file %s\n",filename);
            origColorterm = 99.0;
            origErrcolor = 99.0;
            origFlag = 99;
          } 
        }       
      }
      newColorterm = 99.0;
      newErrcolor = 99.0;
      newFlag = 99;
#ifdef SECOND_COMPARE
      sprintf(filename,"%s%s%s%d%s","/dasch/junk/colorterm/octave2/",platelist[index],"_a",bin,".colorterm.txt");
#else /* SECOND_COMPARE */
      sprintf(filename,"%s%s%s%d%s","/dasch/junk/colorterm/octave/",platelist[index],"_a",bin,".colorterm.txt");
#endif /* SECOND_COMPARE */
      printf("INFO: opening %s\n",filename);
      fileHandle = fopen(filename,"rt");
      if (fileHandle == NULL) {
        printf("ERROR: Could not open file %s\n",filename);
        
      } else {
        inBuffer = fgets(inLine,MAX_BUFFER,fileHandle);
        fclose(fileHandle);
        if (inBuffer == NULL) {
          printf("ERROR: No input buffer for file %s\n",filename);
        } else {
          nvals = sscanf(inBuffer,"%lf %lf %d",&newColorterm,&newErrcolor,&newFlag);
          if (nvals != 3) {
            printf("ERROR: nvals is %d for file %s\n",filename);
            newColorterm = 99.0;
            newErrcolor = 99.0;
            newFlag = 99;
          } 
        }       
      }
      printf("RESULT: colorterm %6.2f %6.2f %6.2f  errcolor %6.2f %6.2f %6.2f  Flag %2d %2d %2d for %s_a%d\n",
             origColorterm,
             newColorterm,
             origColorterm-newColorterm,
             origErrcolor,
             newErrcolor,
             origErrcolor-newErrcolor,
             origFlag,
             newFlag,
             origFlag-newFlag,
             platelist[index],
             bin);
      if ((origColorterm < 90.0) && (newColorterm < 90.0)) {
        fprintf(gscHandle,"%f %f\n",origColorterm,newColorterm);
      }

    } /* Bin index */
  

  } /* Plate Index */
}
