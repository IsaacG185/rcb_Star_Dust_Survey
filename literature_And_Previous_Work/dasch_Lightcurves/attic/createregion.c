// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* createregion.c
 *
 *  Generate region files for ds9 from starbase inputs
 * 
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm  createregion.c pipelineutils.a -ltable -lutil -lwcs  -L/usr/lib64/mysql -lmysqlclient -o createregion 
 * 
 *
 *  createregion -i /home/scanner/Pipeline/match/rb00273_00_01r180ww_tnx.db -o /home/scanner/junk/rb00273.reg -r 187.277908333 -d 2.052325000 -a 0.59 -l
 *
 *  Jul  2, 2013 Edward J. Los - Initial version
 *  Feb  1, 2014 Edward J. Los - Add -m qualifier to compute aLength and bLength from ISO0
 */   

#include <math.h>
#include <errno.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"


/* #define SKIP_PROPER_MOTION 1 */
#define FILTER_PROPER_MOTION 1 /* Use only Tycho-2 and Skymap proper motions */
double str2dec(const char* in);


#define MAX_BUFFER 512
	
extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;


/* reformatted Sextractor results */
typedef struct _inputdata {
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
	double X_IMAGE;
	double Y_IMAGE;
	double THETA_IMAGE;
	double aLength;
	double bLength;
	double ELLIPTICITY;
	int NUMBER;
	int ISO0;
	char REF[MAX_REF];
} INPUTDATA,*PINPUTDATA;


int main(int argc,char *argv[])
{
  char *argstr;
  char cmdchar;
  int nvals;
  int row_result;
  char outputname[MAX_BUFFER];
  char estimatename[MAX_BUFFER];
  char radiustext[MAX_BUFFER];
  char extinctionfile[MAX_BUFFER];
  FILE * outputHandle = NULL;
  FILE * estimateHandle = NULL;
	PINPUTDATA pInputData;
	double centerRa = -1;
	double centerDec = -1;
	double factor;
	int outputCount = 0;
  /* NOTE: Most filtering according to MAG_ISO was effectively removed on Jan 29, 2009 by increasing
     the MAG_ISO_REJECT_FACTOR to a high value of 3000 */
  /* The MAG_ISO_Table is a histogram of all input points as a function of MAG_ISO and spatial_bin */
  int *MAG_ISO_Table = NULL;
  /* The magnitudeTable is a histogram of all output matches as a function of (MAG_ISO,Stdmag) and spatial_bin */
  int *magnitudeTable = NULL;
  int errorFlag = 0;
  int verbose = 0;
	int aLengthFlag = 0;
	int ISO0Flag = 0;
  double limitingMagnitude = 0.0;
  double matchRadius = 0.0;
  time_t startTime;
  time_t curTime;
  int numCacheLines; /* Number of input file lines to keep at any one time */
  int numRows = 0;
  int memoryAllocated = 0;
  int index;
  int maxStarsPerRow = 0;
  double minDec;
  double maxDec;
  double minRa[2];
  double maxRa[2];
  int curCatalogDecBin = -1;
  int oldCatalogDecBin = -1;
  int tmpCatalogDecBin;
  double curCatalogDec;
  double maxCatalogDec;
  double minCatalogDec;
  double scale = 0.0; /* Plate scale in arcsec per pixel */

  char indexname[MAX_BUFFER];
  char* catalogdir;
  char* slashPtr;
	char qualifier[MAX_BUFFER];
  char catalogname[MAX_BUFFER];
  int catalogNumber = 0;
  File catalogHandle;
  PSTARINDEX pStarIndex = NULL;
  int starIndexAlloc = 0;
  int starIndexEntries;
  int gscImageAlloc = 0;
  int gscImageEntries;
  char * dotPtr;
  int matchedStars = 0;
  double geocentricJD = 0.0;
  double plateepoch;
	double circleRadius = -1;
	double regionRadius = -1;
  File input_handle = NULL;
  char input_name[MAX_BUFFER];
  TableHead input_header = NULL;
  size_t input_nrecs = 0;
  size_t input_index;

	
  PINPUTDATA input_table = NULL;
  INPUTDATA input_block;
  PINPUTDATA pInput = &input_block;
  signed long long curMemory = 0;
  signed long long maxMemory = 0;
  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int properMotionSkipCount = 0;
  int duplicateCount = 0;
  MYSQL my_connection;
  char *mysqlhost;
  char *username;
  char *password;
  char series[MAX_SERIES_STRING];
  int plateNumber;
	int mosaicNumber;
	int solutionNumber = 0;
	int result;
	int binning;
	int rotation;
	double aLength;
	double bLength;

  memset(pInput,0,sizeof(INPUTDATA));
  /* Loop through the arguments */
  outputname[0] = 0;
  input_name[0] = 0;
  radiustext[0] = 0;
  estimatename[0] = 0;
  qualifier[0] = 0;
#ifdef SKIP_PROPER_MOTION
  printf("WARNING: SKIP_PROPER_MOTION is set in createregion degrees %f\n",properMotionDegrees);
#endif /* SKIP_PROPER_MOTION */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: unqualified argument %s\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'v':
        case 'V':
          verbose = 1;
          break;

        case 'l':
        case 'L':
          aLengthFlag = 1;
          break;

        case 'm':
        case 'M':
          ISO0Flag = 1;
          break;


        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outputname,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'i': /* input file name */
        case 'I':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(input_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'c': /* circle radius */
        case 'C':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&circleRadius);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the circle radius %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'a': /* region radius */
        case 'A':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&regionRadius);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the region radius %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'd': /* center declination */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&centerDec);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the center declination %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'r': /* center Right Ascension */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&centerRa);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the center Right Ascension %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;



        default:
          printf("ERROR: illegal command -%c-",cmdchar);
          errorFlag = 1;
          break;
        }
        
      }
    }
  }

  /* Verify that we have everything */
  if (outputname[0] == 0) {
    printf("ERROR: No output file specified\n");
    errorFlag = 1;
  }
  if (input_name[0] == 0) {
    printf("ERROR: No input file specified\n");
    errorFlag = 1;
  }
	if ((aLengthFlag == 0) && 
			(circleRadius < 0)) {
		printf("ERROR: No circle radius specified\n");
		errorFlag = 1;
	}
	
  if (regionRadius < 0) {
    printf("ERROR: No region radius specified\n");
    errorFlag = 1;
  }
  if (centerRa < 0) {
    printf("ERROR: No center right ascension specified\n");
    errorFlag = 1;
  }
  if (centerDec < 0) {
    printf("ERROR: No center declination specified\n");
    errorFlag = 1;
  }

  /* Now verify that we have all of the files */
  input_handle = Open(input_name,"rt");
  if (input_handle == NULL) {
    printf("ERROR: Could not open file %s\n",input_name);
    errorFlag = 1;
  } else {
    input_header = table_header(input_handle,TABLE_PARSE);
    if (input_header == NULL) {
      printf("ERROR: Failed to read header for %s\n",input_name);
      errorFlag = 1;
    } 
  }


  outputHandle = fopen(outputname,"wt");
  if (outputHandle == NULL) {
    printf("ERROR Failed to open the output file %s\n",outputname);
    errorFlag = 1;
  }


  if (errorFlag) {
    printf("Usage: createregion -i <input name> \n");
    printf("                  -o <output name> \n");
    printf("                  -r <right ascension in degrees> \n");
    printf("                  -d <declination degrees> \n");
    printf("                  -a <radius around ra and dec in degrees> \n");
    printf("                  -c <circle radius in degrees> \n");
		printf("                  -l use aLength and bLength\n");
		printf("                  -m Compute aLength and bLength from ISO0\n");
    printf("                  -v verbose\n");

    return(-1);
  }

  mysql_init(&my_connection);

  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    printf("ERROR: DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    printf("ERROR: DASCH_USERNAME is not defined\n");
    return(-1);
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    printf("ERROR: DASCH_PASSWORD is not defined\n");
    return(-1);
  }
                   

  if (mysql_real_connect(&my_connection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
#if 0
    printf("Connection success\n");
#endif
  } else {
    printf("ERROR: Connection failed\n");
    if (mysql_errno(&my_connection)) {
      printf("ERROR: Connection error %d: %s\n",mysql_errno(&my_connection),mysql_error(&my_connection));
    }
    return(-1);
  }

  if (ParseFilename(input_name,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
		printf("ERROR: Failed to parse filename %s as %s%05d_%02d_%02d rotation %d\n",
						input_name,series,plateNumber,mosaicNumber,binning,rotation);
		return(-1);
	}

	result = GetPlateScale(&my_connection,series,plateNumber,mosaicNumber,solutionNumber,&mosaicWidth,&mosaicHeight,&scale);
	if (result != 0) {
		printf("ERROR: can not get plate scale for %s\n",input_name);
		exit(-1);;
	}
	printf("Plate scale for %s%05d_%02d_%02d rotation %d is %f arcsec/pixel\n",
				 series,plateNumber,mosaicNumber,binning,rotation,scale);


	input_table = table_loadva(input_handle,
														 &input_header,
														 NULL, /* hbase */
														 NULL, /* rows */
														 NULL,
														 sizeof(INPUTDATA),
														 &input_nrecs,
#if 0
														 TblBuf,"REF"    ,TblOff(PINPUTDATA,REF),MAX_REF,
#endif
														 TblInt,"NUMBER",TblOff(PINPUTDATA,NUMBER),
														 TblInt,"ISO0",TblOff(PINPUTDATA,ISO0),
														 TblDbl,"ra"  ,TblOff(PINPUTDATA,ra),
														 TblDbl,"dec"  ,TblOff(PINPUTDATA,dec),
														 TblDbl,"X_IMAGE"  ,TblOff(PINPUTDATA,X_IMAGE),
														 TblDbl,"Y_IMAGE"  ,TblOff(PINPUTDATA,Y_IMAGE),
														 TblDbl,"THETA_IMAGE"  ,TblOff(PINPUTDATA,THETA_IMAGE),
														 TblDbl,"aLength"  ,TblOff(PINPUTDATA,aLength),
														 TblDbl,"bLength"  ,TblOff(PINPUTDATA,bLength),
														 TblDbl,"ELLIPTICITY",TblOff(PINPUTDATA,ELLIPTICITY),
														 0,"end",0);

  if (input_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",input_name);
    return(-1);
  }
	printf("read %d records for %s\n",input_nrecs,input_name);

	factor = cos(centerDec * DEGREES_TO_RAD);
	for (input_index = 0; input_index < input_nrecs; input_index++) {
		pInputData = &input_table[input_index];
#if 0
		if (pInputData->NUMBER == 561886)	 {
			printf("At NUMBER %d\n",pInputData->NUMBER);
		}
#endif
		if ((pInputData->ra < (centerRa+(regionRadius/factor))) &&
				(pInputData->ra > (centerRa-(regionRadius/factor))) &&
				(pInputData->dec < (centerDec + regionRadius)) &&
				(pInputData->dec > (centerDec - regionRadius))) {
			if (ISO0Flag != 0) {
				aLength = sqrt((1.0*pInputData->ISO0)/(PI_VALUE*(1.0 - pInputData->ELLIPTICITY)));
				bLength = aLength *(1.0 - pInputData->ELLIPTICITY);
			} else {
				aLength = pInputData->aLength;
				bLength = pInputData->bLength;
			}


			if (aLengthFlag != 0) {
#if 0
				fprintf(outputHandle,"J2000;BOX(%f,%f,%f,%f,%f) # color = red\n",
								pInputData->ra,
								pInputData->dec,
								2*aLength*scale/3600.,
								0.0,
								pInputData->THETA_IMAGE);

#else
				fprintf(outputHandle,"J2000;BOX(%f,%f,%f,%f,%f) # color = red\n",
								pInputData->ra,
								pInputData->dec,
								2*aLength*scale/3600.,
								2*bLength*scale/3600.,
								pInputData->THETA_IMAGE);
#endif
			} else {
				fprintf(outputHandle,"J2000;CIRCLE(%f,%f,%f) # color = red\n",pInputData->ra,pInputData->dec,circleRadius);
			}
			outputCount++;

		}

	}
	printf("wrote %d records for %s\n",outputCount,outputname);
  if (input_header != NULL) {
    table_hdrfree(input_header);
  }
  if (input_handle != NULL) {
    Close(input_handle);
  }
  if (outputHandle != NULL) {
    fclose(outputHandle);
  }
  mysql_close(&my_connection);


}
