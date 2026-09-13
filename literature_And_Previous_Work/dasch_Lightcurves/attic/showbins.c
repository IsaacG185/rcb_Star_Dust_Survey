// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* showbins.c
 *
 *  Given a pixel value and a mosaic size, this program
 *  shows the spatial bin and the local calibration bin. 
 *
 * gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm -lmysqlclient showbins.c pipelineutils.a -ltable -lutil -lwcs -o showbins 
 * 
 *  Use getfits to translate from ra and dec to pixel value
 *  Use getlocation -a to find width and height
 *
 *  For bm01749
 *  showbins -w 17412 -h 22026 -x 15549 -y 14141
 *
 *  For rh15338
 *  showbins -w 17384 -h 22026 -x 12827 -y 9939
 *
 *  Oct 20, 2008 Edward J. Los - Initial Version
 *  Aug 31, 2009 Edward J. Los - irec is now local_bin_index
 */   

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pipelineutils.h"

#define X_BINS 50
#define Y_BINS 50

int main(int argc,char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  int spatial_bin;
  int verbose = 0;

  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int xPixel = 0;
  int yPixel = 0;
  double edgeDist;
  int nx;             /* Total bins in width */
  int ny;             /* Total bins in height */
  int ix;             /* width bin */
  int iy;             /* height bin */
  int local_bin_index;

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
        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicWidth);
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
            nvals = sscanf(*++argv,"%d",&mosaicHeight);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'x': /* x pixel */
        case 'X':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&xPixel);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode x pixel value %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'y': /* Y pixel  */
        case 'Y':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&yPixel);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode y pixel value %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;


        default:
          fprintf(stderr,"ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */
  if (mosaicWidth == 0) {
    fprintf(stderr,"ERROR: No mosaic width specified\n");
    errorFlag = 1;
  }
  if (mosaicHeight == 0) {
    fprintf(stderr,"ERROR: No mosaic height specified\n");
    errorFlag = 1;
  }
  if (xPixel == 0) {
    fprintf(stderr,"ERROR: No x pixel value specified\n");
    errorFlag = 1;
  }
  if (yPixel == 0) {
    fprintf(stderr,"ERROR: No y pixel value specified\n");
    errorFlag = 1;
  }


  if (errorFlag) {
    fprintf(stderr,"Usage: showbins\n");
    fprintf(stderr,"                    -w <mosaic width in pixels> \n");
    fprintf(stderr,"                    -h <mosaic height in pixels> \n");
    fprintf(stderr,"                    -x <x pixel value> \n");
    fprintf(stderr,"                    -y <y pixel value> \n");
    fprintf(stderr,"                    -v : verbose\n");

    return(-1);
  }

  printf("showbins of %s %s\n",
         __DATE__,__TIME__);
  printf("Mosaic width %d height %d\n",mosaicWidth,mosaicHeight);
  printf("X Pixel value %d Y pixel value %d\n",xPixel,yPixel);
  spatial_bin = CalculateBin(mosaicWidth,mosaicHeight,xPixel,yPixel,&edgeDist);
  printf("spatial bin %d edge distance %.2f\n",spatial_bin,edgeDist);

  nx = X_BINS;
  ny = Y_BINS;

  ix = (xPixel * nx) /(1.0 * mosaicWidth);
  iy = (yPixel * ny) /(1.0 * mosaicHeight);
  if (ix < 0) {
    ix = 0;
  }
  if (ix >= nx) {
    ix = nx-1;
  }
  if (iy < 0) {
    iy = 0;
  }
  if (iy >= ny) {
    iy = ny-1;
  }
  /* These records are reversed from the local bin */
  local_bin_index = iy + (ny*ix);
  printf("X bins %d ix %d Y bins %d iy %d local_bin_index %d\n",nx,ix,ny,iy,local_bin_index);

  return(0);
}
