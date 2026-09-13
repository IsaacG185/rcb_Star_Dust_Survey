// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* checksize.c
 *
 * This routine looks for a file's existence and makes sure that the file is not
 * just a starbase or octave header.  If the file is good, return a status of 0
 * Otherwise return a status of -1;
 *
 * cc -g -O0 -o checksize checksize.c
 *
 *   checksize -v -f checksize.c -l 10000
 *   echo $status
 *   result is 0: input error or file o.k.
 *             1: no such file
 *             2: file has zero length
 *             3: can not open file
 *             4: file has too few lines
 *
 *  Apr 13, 2011 Edward J. Los - Initial version
 *
 */
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#define MAX_PATH 512
#define MAX_BUFFER 512
#define LESS_VERBOSE 1
int main(int argc,char *argv[])
{
  char *argstr;
  int errorFlag = 0;
  int minimumLines = 3;
  int nvals;
  int verbose = 0;
  char cmdchar;
  char filename[MAX_PATH];
  int statResult;
  struct stat statbuf;
  int result = 0;
  FILE *fileHandle;
  int numlines = 0;
  char *inBuffer;
  char inLine[MAX_BUFFER];

  /* Loop through all the arguments */
  filename[0] = 0;
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'l': /* Minimum number of lines in the file */
        case 'L':
          if (argc > 1) {
            argc--;
            nvals = sscanf(*++argv,"%d",&minimumLines);
            if (nvals != 1) {
              printf("ERROR: Can not decode the minimum number of lines\n");
              errorFlag = 1;
            }
          } else {
            printf("ERROR: No line number follows -l\n");
            errorFlag = 1;
          }

        
          break;
        case 'f': /* Select mosaics from one series */
        case 'F':
          if (argc > 1) {
            strncpy(filename,*++argv,MAX_PATH-1);
            argc--;
          } else {
            printf("ERROR: No series identifier follows -s\n");
            errorFlag = 1;
          }
          break;
        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;
         
        default:
          printf("ERROR: illegal command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  if (filename[0] == 0) {
    printf("ERROR: No filename specified\n");
    errorFlag = 1;
  }
  if (errorFlag) {
    printf("Usage: checksize -f <filename> [-l <minimum lines>][ -v]\n");
    printf("                  where -f is the filename to check\n");
    printf("                  where -l is the minimum number of lines\n");
    printf("                  where -f is verbose mode\n");
    exit(0);
  }
#ifndef LESS_VERBOSE
  if (verbose) {
    printf("checkstatus of %s %s min lines %d for %s n",
           __DATE__,__TIME__,minimumLines,filename);
  }
#endif
  statResult = stat(filename,&statbuf);
  if (statResult == 0) {
    if (statbuf.st_size > 0) {
      fileHandle = fopen(filename,"rt");
      if (fileHandle != NULL) {
        while (1) {
          inBuffer = fgets(inLine,MAX_BUFFER,fileHandle);
          if (inBuffer == NULL) {
            break;
          }

          numlines++;
          if (numlines >= minimumLines) {
            break;
          }
        }
        if (numlines >= minimumLines) {
          result = 0;
        } else {
          result = 4; /* Not enough lines */
        }
      } else {
        result = 3; /* Can not open file */
      }
    } else {
      result = 2; /* Zero size file */
    }
  } else {
    result = 1; /* File not found */
  }

  if (verbose) {

#ifdef LESS_VERBOSE
    printf("result %d %d for %s\n",result,numlines,filename);
#else /* LESS_VERBOSE */
    printf("numlines %d result: %d\n",numlines,result);
#endif /* LESS_VERBOSE */

  }
  exit(result);
}
