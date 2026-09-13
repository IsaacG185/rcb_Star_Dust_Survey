// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* setPlateQuality.c
 * 
 *  This routine sets and clears named bits in the quality column of the plates database
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I/usr/include/mysql  -I /dasch/install/include  -L /dasch/install/lib -lm  -L/usr/lib64/mysql -L/usr/lib/mysql -lmysqlclient setPlateQuality.c pipelineutils.a -ltable -lutil -lwcs -o setPlateQuality 
 *
 *  
 *   To set   the bit: setPlateQuality -p i31013_00_01r180ww -b wedge   -s
 *   To clear the bit: setPlateQuality -p i31013_00_01r180ww -b wedge   -c
 *   To read  the bit: setPlateQuality -p i31013_00_01r180ww -b wedge   -r
 *
 * Dec 16, 2009 Edward J. Los - Initial version
 */


#include <math.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#define MAX_BUFFER 512

extern QUALITYBIT qualityMasks[];

int main(int argc,char *argv[])
{
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  int nvals;
  int setFlag = -1;
  int readFlag = 0;
  char bitname[MAX_BUFFER];
  char platename[MAX_BUFFER];
  PQUALITYBIT pQualityBit = qualityMasks;
  int quality;
  platename[0] = 0;
  bitname[0] = 0;
 
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      fprintf(stderr,"ERROR: Stray argument %s\n",argstr);
      errorFlag = 1;

    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'b': /* bitname */
        case 'B':
          argc--;
	if (argc < 1) {
	  fprintf(stderr,"ERROR: setPlateQuality  Insufficient arguments for -%c\n",cmdchar);
	  errorFlag = 1;
	} else {
	  strncpy(bitname,*++argv,MAX_BUFFER-2);
	  if (strlen(*argv) >= MAX_BUFFER-2) {
	    fprintf(stderr,"ERROR: setPlateQuality  MAX_BUFFER too small for %s\n",*argv);
	  }
	}
	break;

        case 'p': /* plate name */
        case 'P':
          argc--;
	if (argc < 1) {
	  fprintf(stderr,"ERROR: setPlateQuality  Insufficient arguments for -%c\n",cmdchar);
	  errorFlag = 1;
	} else {
	  strncpy(platename,*++argv,MAX_BUFFER-2);
	  if (strlen(*argv) >= MAX_BUFFER-2) {
	    fprintf(stderr,"ERROR: setPlateQuality  MAX_BUFFER too small for %s\n",*argv);
	  }
	}
	break;



        case 'r': /* read the bit */
        case 'R':
	  if (setFlag >= 0) {
	    fprintf(stderr,"ERROR: setFlag is already %d\n",setFlag);
	    errorFlag = 1;
	  }
	readFlag = 1;
	break;

        case 's': /* set the bit */
        case 'S':
	  if (setFlag >= 0) {
	    fprintf(stderr,"ERROR: setFlag is already %d\n",setFlag);
	    errorFlag = 1;
	  }
	if (readFlag > 0) {
	  fprintf(stderr,"ERROR: readFlag is already %d\n",readFlag);
	  errorFlag = 1;
	}

	setFlag = 1;
	break;

        case 'c': /* clear the bit */
        case 'C':
	  if (setFlag >= 0) {
	    fprintf(stderr,"ERROR: setFlag is already %d\n",setFlag);
	    errorFlag = 1;
	  }
	if (readFlag > 0) {
	  fprintf(stderr,"ERROR: readFlag is already %d\n",readFlag);
	  errorFlag = 1;
	}
	setFlag = 0;
	break;


        default:
          fprintf(stderr,"ERROR: setPlateQuality   unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */

  if (bitname[0] == 0) {
    fprintf(stderr,"ERROR: setPlateQuality  No bit name  specified\n");
    errorFlag = 1;
  }
  if (platename[0] == 0) {
    fprintf(stderr,"ERROR: setPlateQuality  No platename was specified\n");
    errorFlag = 1;
  }
  if ((setFlag < 0) && (readFlag == 0)) {
    fprintf(stderr,"ERROR: setPlateQuality  Set, clear, or read not specified\n");
    errorFlag = 1;
  }
  if (errorFlag) {
    fprintf(stderr,"Usage: setPlateQuality -b <bitname> -p <platename> -s|-c|-r \n");
    fprintf(stderr,"                    -s sets the bit \n");
    fprintf(stderr,"                    -c clears the bit \n");
    fprintf(stderr,"                    -r read the bit\n");
    while (pQualityBit->qualityMask != 0) {
      fprintf(stderr,"Bitname: %4d 0x%04x %-20s\n",pQualityBit->qualityMask,pQualityBit->qualityMask,pQualityBit->qualityDescr);
      pQualityBit++;
    }
   
    return(-1);
  }

  printf("setPlateQuality of %s %s bit %s plate %s setFlag %d",
         __DATE__,__TIME__,bitname,platename,setFlag);

  SetPlateQuality2(bitname,platename,setFlag,readFlag,&quality);
  if (readFlag) {
    printf(" quality %d 0x%x ",quality,quality);
    while (pQualityBit->qualityMask != 0) {
      if ((pQualityBit->qualityMask & quality) != 0) {
	printf(" %s",pQualityBit->qualityDescr);
      }
      pQualityBit++;
    }
    
  }
  printf("\n");
  

  return(0);
}

