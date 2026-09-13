// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* add_galactic.c - add galactic coordinate columns
 *
 *  Accepts a three column input of the form REF ra dec and produces REF lng lat
 *  where ra is the Right Ascension in degrees
 *        dec is the declination in degrees
 *        lon is the galactic longitude in degrees
 *        lat is the galactic latitude in degrees
 *

 * Without "nearbyObjects"

   sorttable REF < los2.db > los6.db
   column -i los2.db REF ra declination dec | compute 'dec=declination' > los3.db
   sorttable REF  < los3.db > los4.db
   add_galactic -i los4.db -o los5.db
   jointable -j REF los6.db los5.db > los7.db 


* With "nearby Objects "

   For the "id_unmatched" table:

   sorttable REF < los2.db > los6.db
   column -i los2.db REF ra declination dec nearbyObjects | compute 'dec=declination' > los3.db
   sorttable REF  < los3.db > los4.db
   add_galactic -n -i los4.db -o los5.db
   jointable -j REF los6.db los5.db | compute 'nearbyObjects = notes'  > los11.db 
   cp los11.db ~/Pipeline/ingest


 *
 * 
 *  
 * Oct 26, 2012 Edward J. Los - initial version
 * Jan 10, 2012 Edward J. Los - Add nearbyObjects format support.
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 *  
 */


#include <math.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include "searchgsc.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "daschunistd.h"
typedef struct inputrecord {
  double ra;          /* Right Ascension */     
  double dec;         /* Declination */
  double lon;         /* Galactic longitude */
  double lat;         /* Galactic latitude */
  char REF[MAX_REF];  /* catalog reference number */
  char nearbyObjects[MAX_NEARBY_OBJECTS_STRING+1];
} INPUTRECORD,*PINPUTRECORD;

#define MAX_BUFFER 512
#define MAX_FILENAME 512




int main(int argc,char *argv[])
{
  char input_name[MAX_BUFFER]; /* -t qualifier */
  char output_name[MAX_BUFFER]; /* -t qualifier */
  File input_handle = NULL;
  TableHead input_header = NULL;
  PINPUTRECORD input_table = NULL;
  PINPUTRECORD pStar;
  size_t input_nrecs = 0;
  size_t input_index = 0;
  double lontemp = 0;
  double lattemp = 0;
  FILE *outHandle = NULL;
  char *argstr;
  int errorFlag = 0;
  int verbose = 0;
  int nearbyObjectsFlag = 0;
  char cmdchar;
  char *charPtr;
  input_name[0] = 0;
  output_name[0] = 0;
 

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d\n",argstr,argc);
      
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {



        case 'i': /* input */
        case 'I': 
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: text file name length %d for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(input_name,*argv);
            }
          }
          break;

        case 'o': /* output */
        case 'O': 
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: text file name length %d for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(output_name,*argv);
            }
          }
          break;
        case 'V':
        case 'v':
          verbose = 1;
          break;

        case 'N':
        case 'n':
          nearbyObjectsFlag = 1;
          break;


        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }

        
      }

    }
  }

  if (strlen(input_name) == 0) {
    printf("ERROR: input_name is not specified\n");
    errorFlag = 1;
  }
  if (strlen(output_name) == 0) {
    printf("ERROR: output_name is not specified\n");
    errorFlag = 1;
  }

  outHandle = fopen(output_name,"wt");
  if (outHandle == NULL) {
    errorFlag = 1;
    fprintf(stderr,"ERROR: Failed to open the output file %s\n",output_name);
  } else {
    if (nearbyObjectsFlag) {
      fprintf(outHandle,"REF\tlon\tlat\tnotes\n");
      fprintf(outHandle,"---\t---\t---\t-----\n");
    } else {
      fprintf(outHandle,"REF\tlon\tlat\n");
      fprintf(outHandle,"---\t---\t---\n");
    }
  }
  

  if (errorFlag == 1) {
    printf("Usage: add_galactic -i <input file> -o <output file>\n");
    printf("       where the input file is a starbase table with  columns REF,ra,dec  (equatorial J2000)\n");
    printf("             the output file is a starbase table with columns REF,lon,lat (galactic   J2000)\n");
    printf("             -n make sure the nearbyObjects field can be read by octave by converting spaces\n");

    exit(-1);
  }
  
  /* Here we have already prepared the data.  Read it in */
  input_handle = Open(input_name,"r");
  if (input_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the db file %s\n",input_name);
  } else {
    if (verbose) {
      printf("Found db file %s\n",input_name);
    }
  }
  
  input_header = table_header(input_handle,TABLE_PARSE);
  if (input_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",input_name);
    return(-1);
  }

  if (nearbyObjectsFlag) {
    input_table = table_loadva(input_handle,
                               &input_header,
                               NULL, /* hbase */
                               NULL, /* rows */
                               NULL,
                               sizeof(INPUTRECORD),
                               &input_nrecs,
                               TblDbl,"ra",TblOff(PINPUTRECORD,ra),
                               TblDbl,"dec",TblOff(PINPUTRECORD,dec),
                               TblBuf,"REF",TblOff(PINPUTRECORD,REF),MAX_REF,
                               TblBuf,"nearbyObjects",TblOff(PINPUTRECORD,nearbyObjects),MAX_NEARBY_OBJECTS_STRING,
                               0,"end",0);
  } else {
    input_table = table_loadva(input_handle,
                               &input_header,
                               NULL, /* hbase */
                               NULL, /* rows */
                               NULL,
                               sizeof(INPUTRECORD),
                               &input_nrecs,
                               TblDbl,"ra",TblOff(PINPUTRECORD,ra),
                               TblDbl,"dec",TblOff(PINPUTRECORD,dec),
                               TblBuf,"REF",TblOff(PINPUTRECORD,REF),MAX_REF,
                               0,"end",0);

  }
  if (input_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",input_name);
    return(-1);
  }
  
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pStar = &input_table[input_index];
    lontemp = pStar->ra;
    lattemp = pStar->dec;
    wcscon(WCS_J2000,WCS_GALACTIC,2000.0,2000.0,&lontemp,&lattemp,2000.0);
    pStar->lon = lontemp;
    pStar->lat = lattemp;
    if (verbose) {
      printf("REF %s ra %f dec %f lon %f lat %f\n",pStar->REF,pStar->ra,pStar->dec,pStar->lon,pStar->lat);
    }
    if (nearbyObjectsFlag) {
      if (strlen(pStar->nearbyObjects) == 0) {
        /* Octave does not tolerate a null string */
        strcpy(pStar->nearbyObjects,"X");
      } else {
        while ((charPtr = strstr(pStar->nearbyObjects," ")) != NULL) {
          /* Octave uses spaces as well as tabs for separators */
          *charPtr = '_';
        }
      }
      fprintf(outHandle,"%s\t%f\t%f\t%s\n",pStar->REF,pStar->lon,pStar->lat,pStar->nearbyObjects);
    } else {
      fprintf(outHandle,"%s\t%f\t%f\n",pStar->REF,pStar->lon,pStar->lat);
    }
  }
  fclose(outHandle);


  if (input_handle != NULL) {
    Close(input_handle);
  }
  if (input_header != NULL) {
    table_hdrfree(input_header);
  }

  if (input_table != NULL) {
    free(input_table);
  }


  return(0);
}

