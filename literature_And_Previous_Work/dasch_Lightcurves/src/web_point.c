// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* extract data about a single point, including the image if it is available
 *
 * Invokes the `showflags` subprogram.
 */

#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

#include "table.h"

#include "mysql.h"

#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "libwcs/fitswcs.h"

#include "scandb.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include "searchgsc.h"

#define MAX_BUFFER 512
#define MAX_FILENAME 512

extern QUALITYBIT webQualityMasks[];
extern int webQualityTableSize;

extern char *catalogText[MAX_CATALOG_NUMBER];

/* Provided by libwcs but not declared */
extern struct WorldCoor *GetFITSWCS(
  char *filename,
  char *header,
  int verbose,
  double *cra,
  double *cdec,
  double *dra,
  double *ddec,
  double *secpix,
  int *wp,
  int *hp,
  int *sysout,
  double *eqout
);

static int sysim = 0;
static double eqim = 0.0;
static SERIESLIST seriesList[MAX_SERIES+1];
static int seriesCount = 0;
static int allSeriesMask = 0;


int
StarDateCompare(const void *first, const void *second)
{
  double dateFirst = ((PFILESTARIMAGEEXT)first)->filestarimage.Date;
  double dateSecond = ((PFILESTARIMAGEEXT)second)->filestarimage.Date;

  if (dateFirst > dateSecond) {
    return 1;
  } else if (dateFirst < dateSecond) {
    return -1;
  } else {
    return 0;
  }
}


int
CheckImageLocation(
  char *filename,
  double averageRa,
  double averageDec,
  char *pSpatial_bin,
  int *pOffscl,
  double *pEdgeDistance
) {
	char *header = NULL;
  double cra;
  double cdec;
  double dra;
  double ddec;
  double secpix;
  int wp;
  int hp;
	int verbose = 0;
	double xpix;
	double ypix;
	int offscl = 0;
	double edgeDistance = 0; /* Distance from the plate edge in cm. */
	struct WorldCoor *wcs = NULL;		/* World coordinate system structure */
	int spatial_bin;
	double pNullDist;

	*pOffscl = offscl;
	*pEdgeDistance = 0;

	if ((header = GetFITShead(filename, verbose)) == NULL) {
		return -1;
	}

	wcs = GetFITSWCS(
    filename,
    header,
    verbose,
    &cra,
    &cdec,
    &dra,
    &ddec,
		&secpix,
    &wp,
    &hp,
    &sysim,
    &eqim
  );

	if (nowcs(wcs)) {
		wcsfree(wcs);
		wcs = NULL;
		return -1;
	}

	wcs2pix(wcs, averageRa, averageDec, &xpix, &ypix, &offscl);

	*pOffscl = offscl;

	if (offscl == 0) {
		spatial_bin = CalculateBin(wp, hp, xpix, ypix, &pNullDist);
		*pSpatial_bin = spatial_bin;
		edgeDistance = xpix;

		if (edgeDistance > ypix) {
			edgeDistance = ypix;
		}

		if (edgeDistance > ((1.0*wp) - xpix)) {
			edgeDistance = ((1.0*wp) - xpix);
		}

		if (edgeDistance > ((1.0*hp) - ypix)) {
			edgeDistance = ((1.0*wp) - xpix);
		}

		edgeDistance = edgeDistance * NOMINAL_MM_PER_PIXEL / 10.0; /* convert pixels to cm */
		*pEdgeDistance = edgeDistance;
	}

	if (header != NULL) {
		free(header);
	}

  if (wcs != NULL) {
    wcsfree (wcs);
    wcs = NULL;
  }

	return 0;
}


int
main(int argc, char *argv[])
{
  int nvals;
  char *argstr;
  char cmdchar;
  char REF[MAX_REF]; /* GSC2.3.2 reference number */
  long long REFNumber;
  int RefType;
  int errorFlag = 0;
  char *tmpdir = NULL;
  char tmpdir2[MAX_BUFFER];
  int catalogNumber = 0;
  char source[MAX_BUFFER];
  char binaries[MAX_BUFFER];
  char catalogString[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int reqAFLAGSMASK;
  int reqqualitymask;
  int reqseriesmask;
  double reqstartdate;
  double reqenddate;
  int gsc_bin_index = -1;
  int index;
  PQUALITYBIT pQualityEntry;
  int imagesize = -1; /* Radius of the extracted image in pixels */
  char *sizeunits = NULL;
  char defaultsizeunits[] = THUMBNAIL_PIXELS;
  double imagearcmin;
  int plotWidth = -1;
  int plotHeight = -1;
  int automode = -1;
  char textfilename[MAX_BUFFER]; /* -t qualifier */
  char db_name[MAX_FILENAME];
  char filename[MAX_FILENAME];
  char filename2[2048];
  char rotsuffix[MAX_BUFFER];
  File db_handle = NULL;
  char tempname[MAX_FILENAME];
  FILE *tempHandle = NULL;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  TableHead db_header = NULL;
  PFILESTARIMAGEEXT db_table = NULL;
  size_t db_nrecs = 0;
	size_t db_alloc = 0;
  MOSAIC mosaicInfo;
  PMOSAIC pMosaic = &mosaicInfo;
  double scaleFactor = 1.0;
  int verbose = 0;
  int m44release = 0;
  int ndata = 0;
  int idata;
  int bestidata;
  int *crossindex = NULL;
  double *x = NULL;
  double *y = NULL;
  char *dotloc;
  int statResult;
  struct stat statbuf;
  int magnitudeIndex;
  PFILESTARIMAGEEXT pFileStarImage = NULL;
  char store[MAX_FILENAME] = "";
  char storefilename[MAX_FILENAME] = "";
  char *charPtr;
  int ix = -1; /* x Position on graph of the point */
  int iy = -1; /* y Position on graph of the point */
  int nextPointIndex = -1; /* Next point */
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  PHOTPLATES basePhotPlates;
  PPHOTPLATES pPhotPlates = &basePhotPlates;
  char cmdStr[MAX_BUFFER];
  int result;
  double dateyear;
  /* Grzegorz Pomanski's plotting parameters */
  FILE *fp;
  int sx;
  int x0;
  int x1;
  double fx;
  double fx0;
  int sy;
  int y0;
  int y1;
  double fy;
  double fy0;
  double upper_limit;
	int foldflag;
	double foldperiod;
	double foldcenter;
	double folddate;
	int foldfactor;
  int webquality;
  int webseriesbit;
  char* regionFlag = NULL;
  PSERIESLIST pSeriesList;
  int seriesIndex;
  double averageRa = 999.;
  double averageDec = 99.;
  double reqdimmag;
  double reqbrightmag;
	int releaseField;
	double sizePixels;
	double edgeDistance = 0.0;
	int undetectedFlag = 0;
	int offscl = 0;
  int oldEnableRematch = 0;

  memset(seriesList, 0, sizeof(seriesList));

  SetQueryCount(0);
  qualifier[0] = 0;
  source[0] = 0;
  catalogString[0] = 0;
  textfilename[0] = 0;
  storefilename[0] = 0;
  binaries[0] = 0;

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
        case 'd': /* Temporary directory */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            tmpdir = *++argv;
          }
          break;

        case 'f': /* Region Flag */
        case 'F':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            regionFlag = *++argv;
          }
          break;

        case 'g': /* gsc_bin_index */
        case 'G':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&gsc_bin_index);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the gsc_bin_index %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'i': /* Image width in pixels or arcsec */
        case 'I':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            nvals = sscanf(*++argv,"%d",&imagesize);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the imagesize %s\n",*argv);
            }
          }
          break;

        case 'j': /* Units of image width: pixels or arcsec */
        case 'J':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            sizeunits = *++argv;
          }
          break;

        case 'w': /* plot width */
        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&plotWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the plotWidth %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'a': /* automatic mode */
        case 'A':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&automode);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the automode %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* plotHeight */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&plotHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the plotHeight %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'k': /* next point index */
        case 'K':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&nextPointIndex);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the next point index %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'v': /* verbose */
          verbose += 1;
          break;

        case 'm': /* M44 Release */
          m44release += 1;
          break;

        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              strcpy(source,*argv);
              if (catalogNumber > 0) {
                sprintf(catalogString,"%d",catalogNumber);
                strcpy(qualifier,*argv);
              }
            }
          }
          break;

        case 'r': /* Object designation */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n", cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_REF-1)) {
              printf("ERROR: designation length %zu for %s is too long\n", strlen(*argv), *argv);
              errorFlag = 1;
            } else {
              strcpy(REF,*argv);
              GetREFNumber(REF,&REFNumber,&RefType,1,1);

            }
          }
          break;

        case 't': /* Text file name */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: text file name length %zu for %s is too long\n", strlen(*argv), *argv);
              errorFlag = 1;
            } else {
              strcpy(textfilename,*argv);
            }
          }
          break;

        case 'b': /* Executable binary image directory */
        case 'B':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: binaries directory length %zu for %s is too long\n", strlen(*argv), *argv);
              errorFlag = 1;
            } else {
              strcpy(binaries,*argv);
            }
          }
          break;

        case 'n': /* parameter store filename */
        case 'N':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: parameter store filename length %zu for %s is too long\n", strlen(*argv), *argv);
              errorFlag = 1;
            } else {
              strcpy(storefilename,*argv);
            }
          }
          break;

        case 'c': /* Point coordinates (x,y) */
        case 'C':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if ((charPtr = strtok(*argv,",")) != NULL) {
              nvals = sscanf(charPtr, "%d", &ix);

              if (nvals != 1) {
                printf("ERROR: Unable to decode ix in %s\n", *argv);
                errorFlag = 1;
              } else if ((charPtr = strtok(NULL,",")) != NULL) {
                nvals = sscanf(charPtr,"%d",&iy);
                if (nvals != 1) {
                  printf("ERROR: Unable to decode iy in %s\n", *argv);
                  errorFlag = 1;
                }
              } else {
                printf("ERROR: Failed to parse second coordinate %s\n", *argv);
                errorFlag = 1;
              }
            } else {
              printf("ERROR: Failed to parse first coordinate %s\n", *argv);
              errorFlag = 1;
            }
          }
          break;

        case 'O':
          /*enableRematch = 1;*/
          break;

        default:
          printf("ERROR:  unknown command -%c\n", cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  /* Check validity of arguments */

  if (gsc_bin_index < 0) {
    printf("ERROR: gsc_bin_index is not specified\n");
    errorFlag = 1;
  }

  if (imagesize <= 0) {
    printf("ERROR: imagesize is not specified using %d\n", THUMBNAIL_IMAGESIZE);
    imagesize = THUMBNAIL_IMAGESIZE;
  }

  if ((sizeunits == NULL) ||
      ((strcmp(sizeunits,THUMBNAIL_PIXELS) != 0) &&
       (strcmp(sizeunits,THUMBNAIL_ARCSEC) != 0))) {
    sizeunits = defaultsizeunits;
    printf("ERROR: sizeunits is not specified, using %s\n",sizeunits);
  }

  if (plotWidth < 0) {
    printf("ERROR: plotWidth is not specified\n");
    errorFlag = 1;
  }

  if (plotHeight < 0) {
    printf("ERROR: plotHeight is not specified\n");
    errorFlag = 1;
  }

  if (automode < 0) {
    printf("ERROR: automode is not specified\n");
    errorFlag = 1;
  }

  if (strlen(textfilename) == 0) {
    printf("ERROR: textfilename is not specified\n");
    errorFlag = 1;
  }

  if (strlen(binaries) == 0) {
    printf("ERROR: binaries directory is not specified\n");
    errorFlag = 1;
  }

  if (strlen(storefilename) == 0) {
    printf("ERROR: parameter store filename is not specified\n");
    errorFlag = 1;
  }

  if (tmpdir == NULL) {
    printf("ERROR: plot directory is not specified\n");
    errorFlag = 1;
  }

  if (regionFlag == NULL) {
    printf("ERROR: regionFlag is not specified\n");
    errorFlag = 1;
  }

  if (nextPointIndex < 0) {
    if ((ix < 0) ||
        (ix >= plotWidth) ||
        (iy < 0) ||
        (iy >= plotHeight)) {
      printf("ERROR: invalid coordinates %d and %d\n",ix,iy);
      errorFlag = 1;
    }
  }

  if (strlen(textfilename)+strlen(tmpdir)+5 > MAX_BUFFER) {
    printf("ERROR: full name of the text file is too long \n");
    errorFlag = 1;
  }

  if (strlen(storefilename)+strlen(tmpdir)+5 > MAX_BUFFER) {
    printf("ERROR: full name of the parameter store file is too long \n");
    errorFlag = 1;
  }

  strcpy(tempname,tmpdir);
  strcat(tempname,"/");
  strcat(tempname,REF);
  strcat(tempname,".tmp");

  if (errorFlag == 1) {
    printf("Usage: web_point -r <REF> -q <source> -g <gsc_bin_index> -d <directory> -t <text name> -p <plot name> -s <imagesize> [-v]\n");
    printf("       where -v is the verbose flag\n");
    printf("             -r is the object reference\n");
    printf("             -w is the plot width in pixels\n");
    printf("             -h is the plot height in pixels\n");
    printf("             -d is the plot and text directory\n");
    printf("             -t is the text file name\n");
    printf("             -i is the width of the extracted image in units specified by -j below\n");
    printf("             -j units of width of the extracted images - pixels or arcsec");
    printf("             -c is the x,y coordinate of the requested point\n");
    printf("             -k is the index of the next point\n");
    printf("             -n is the parameter store filename\n");
    printf("             -m limited M44 release - no mosaics\n");
    printf("             -a automatic display mode\n");
    printf("             -q <catalog qualifier>\n");
    printf("             -b <executable binary image directory>\n");
    printf("             -f <region flag>\n");
    exit(1);
  }

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pConnection);
  dasch_init_photdb(pPhotConnection);

  InitSeriesTable(pConnection, pPhotConnection);

  /* Read in the plotting parameters */

  strcpy(store,tmpdir);
  strcat(store,"/");
  strcat(store,storefilename);

  if((fp = fopen(store, "r")) != NULL){
    fscanf(fp, "%d %d %d %lf %lf", &sx, &x0, &x1, &fx, &fx0);
    fscanf(fp, "%d %d %d %lf %lf", &sy, &y0, &y1, &fy, &fy0);
    fscanf(fp, "%lf %d %lf %lf %d",&upper_limit,&foldflag,&foldperiod,&foldcenter,&oldEnableRematch);
    fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf",&reqAFLAGSMASK,&reqqualitymask,&reqseriesmask,&reqstartdate,&reqenddate,&reqdimmag,&reqbrightmag,&averageRa,&averageDec);
    fclose(fp);
  } else {
    printf("Cannot open %s\n", store);
    exit(1);
  }

  charPtr = strrchr(tmpdir, '/');
  if (charPtr == NULL) {
    strcpy(tmpdir2, tmpdir);
  } else {
    charPtr++;
    strcpy(tmpdir2, charPtr);
  }

  /* Now read in the starbase file*/
  strcpy(db_name, tmpdir);
  strcat(db_name, "/");
  strcat(db_name, textfilename);
  dotloc = strrchr(db_name, '.');

  if (dotloc) {
    *dotloc = 0;
    strcat(db_name,".db");
  } else {
    printf("ERROR: %s has no suffix\n",textfilename);
    exit(1);
  }

  /* Here we have already prepared the data.  Read it in */
  db_handle = Open(db_name,"r");
  if (db_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the db file %s\n",db_name);
  } else {
    if (verbose) {
      printf("Found db file %s\n",db_name);
    }
  }

  db_header = table_header(db_handle,TABLE_PARSE);
  if (db_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",db_name);
    return 1;
  }

  db_table = table_loadva(db_handle,
                          &db_header,
                          NULL, /* hbase */
                          NULL, /* rows */
                          NULL,
                          sizeof(FILESTARIMAGEEXT),
                          &db_nrecs,
                          TblDbl,"X_IMAGE",TblOff(PFILESTARIMAGEEXT,filestarimage.X_IMAGE),
                          TblDbl,"Y_IMAGE",TblOff(PFILESTARIMAGEEXT,filestarimage.Y_IMAGE),
                          TblDbl,"MAG_ISO",TblOff(PFILESTARIMAGEEXT,filestarimage.MAG_ISO),
                          TblDbl,"ra",TblOff(PFILESTARIMAGEEXT,filestarimage.ra),
                          TblDbl,"dec",TblOff(PFILESTARIMAGEEXT,filestarimage.dec),
                          TblDbl,"Date",TblOff(PFILESTARIMAGEEXT,filestarimage.Date),
                          TblDbl,"FLUX_ISO",TblOff(PFILESTARIMAGEEXT,filestarimage.FLUX_ISO),
                          TblDbl,"MAG_APER",TblOff(PFILESTARIMAGEEXT,filestarimage.MAG_APER),
                          TblDbl,"MAG_AUTO",TblOff(PFILESTARIMAGEEXT,filestarimage.MAG_AUTO),
                          TblDbl,"KRON_RADIUS",TblOff(PFILESTARIMAGEEXT,filestarimage.KRON_RADIUS),
                          TblDbl,"BACKGROUND",TblOff(PFILESTARIMAGEEXT,filestarimage.BACKGROUND),
                          TblDbl,"FLUX_MAX",TblOff(PFILESTARIMAGEEXT,filestarimage.FLUX_MAX),
                          TblDbl,"THETA_J2000",TblOff(PFILESTARIMAGEEXT,filestarimage.THETA_J2000),
                          TblDbl,"ELLIPTICITY",TblOff(PFILESTARIMAGEEXT,filestarimage.ELLIPTICITY),
                          TblDbl,"ISOAREA_WORLD",TblOff(PFILESTARIMAGEEXT,filestarimage.ISOAREA_WORLD),
                          TblDbl,"FWHM_IMAGE",TblOff(PFILESTARIMAGEEXT,filestarimage.FWHM_IMAGE),
                          TblDbl,"FWHM_WORLD",TblOff(PFILESTARIMAGEEXT,filestarimage.FWHM_WORLD),
                          TblDbl,"plate_dist",TblOff(PFILESTARIMAGEEXT,filestarimage.plate_dist),
                          TblDbl,"Blendedmag",TblOff(PFILESTARIMAGEEXT,filestarimage.Blendedmag),
                          TblDbl,"dradRMS2",TblOff(PFILESTARIMAGEEXT,filestarimage.dradRMS2),
                          TblFlt,"magcal_iso",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_iso),
                          TblFlt,"magcal_iso_rms",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_iso_rms),
                          TblFlt,"magcal_local",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_local),
                          TblFlt,"magcal_local_rms",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_local_rms),
                          TblFlt,"limiting_mag_local",TblOff(PFILESTARIMAGEEXT,filestarimage.limiting_mag_local),
                          TblFlt,"magcal_local_error",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_local_error),
                          TblFlt,"magcor_local",TblOff(PFILESTARIMAGEEXT,filestarimage.magcor_local),
                          TblFlt,"extinction",TblOff(PFILESTARIMAGEEXT,filestarimage.extinction),
                          TblFlt,"magcal_magdep",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_magdep),
                          TblFlt,"magcal_magdep_rms",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_magdep_rms),
                          TblFlt,"RaPM",TblOff(PFILESTARIMAGEEXT,filestarimage.RaPM),
                          TblFlt,"DecPM",TblOff(PFILESTARIMAGEEXT,filestarimage.DecPM),
                          TblInt,"A2FLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.A2FLAGS),
                          TblInt,"B2FLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.B2FLAGS),
                          TblFlt,"timeAccuracy",TblOff(PFILESTARIMAGEEXT,filestarimage.timeAccuracy),
                          TblInt,"maskIndex",TblOff(PFILESTARIMAGEEXT,filestarimage.maskIndex),
                          TblInt,"gsc_bin_index",TblOff(PFILESTARIMAGEEXT,filestarimage.gsc_bin_index),
                          TblInt,"plateNumber",TblOff(PFILESTARIMAGEEXT,filestarimage.plateNumber),
                          TblInt,"NUMBER",TblOff(PFILESTARIMAGEEXT,filestarimage.NUMBER),
                          TblInt,"versionId",TblOff(PFILESTARIMAGEEXT,filestarimage.versionId),
                          TblInt,"AFLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.AFLAGS),
                          TblInt,"BFLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.BFLAGS),
                          TblInt,"ISO0",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO0),
                          TblInt,"ISO1",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO1),
                          TblInt,"ISO2",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO2),
                          TblInt,"ISO3",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO3),
                          TblInt,"ISO4",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO4),
                          TblInt,"ISO5",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO5),
                          TblInt,"ISO6",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO6),
                          TblInt,"ISO7",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO7),
                          TblInt,"npoints_local",TblOff(PFILESTARIMAGEEXT,filestarimage.npoints_local),
                          TblInt,"rejectFlag",TblOff(PFILESTARIMAGEEXT,filestarimage.rejectFlag),
                          TblSht,"local_bin_index",TblOff(PFILESTARIMAGEEXT,filestarimage.local_bin_index),
                          TblByt,"seriesId",TblOff(PFILESTARIMAGEEXT,filestarimage.seriesId),
                          TblByt,"exposureNumber",TblOff(PFILESTARIMAGEEXT,filestarimage.exposureNumber),
                          TblByt,"solutionNumber",TblOff(PFILESTARIMAGEEXT,filestarimage.solutionNumber),
                          TblByt,"spatial_bin",TblOff(PFILESTARIMAGEEXT,filestarimage.spatial_bin),
                          TblInt,"mosaicNumber",TblOff(PFILESTARIMAGEEXT,mosaicNumber),
                          TblInt,"quality",TblOff(PFILESTARIMAGEEXT,quality),
                          TblInt,"plateVersionId",TblOff(PFILESTARIMAGEEXT,plateVersionId),
                          TblInt,"catalogNumber",TblOff(PFILESTARIMAGEEXT,filestarimage.catalogNumber),
                          0,"end",0);

  if (db_table == NULL) {
    fprintf(stderr,"ERROR: Failed to read table for %s\n",db_name);
    return 1;
  }

  if (db_handle != NULL) {
    Close(db_handle);
  }

  if (db_header != NULL) {
    table_hdrfree(db_header);
  }

  /* Now decide how many series we have */

  for (magnitudeIndex = 0; magnitudeIndex < db_nrecs; magnitudeIndex++) {
    pFileStarImage = &db_table[magnitudeIndex];
    if ((pFileStarImage->filestarimage.seriesId > 0) &&
        (pFileStarImage->filestarimage.seriesId <= MAX_SERIES)) {
      pSeriesList = &seriesList[pFileStarImage->filestarimage.seriesId];
      pSeriesList->numPlates++;
    }
  }

  seriesCount = START_SERIES_BIT;

  for (seriesIndex = 0; seriesIndex <= MAX_SERIES; seriesIndex++) {
    pSeriesList = &seriesList[seriesIndex];

    if ((pSeriesList->numPlates > 0) &&
        (pSeriesList->bitMask == 0) &&
        (seriesCount <= MAX_SERIES_BIT)) {
      pSeriesList->bitPosition = seriesCount;
      pSeriesList->bitMask = 1 << seriesCount;
      allSeriesMask |= pSeriesList->bitMask;
      seriesCount++;
    } else {
      pSeriesList->numPlates = 0; /* Can not include this series! */
    }
  }

  /* Now create the plotting array */

  if (db_nrecs > 0) {
		if (foldflag) {
			db_alloc = 2*db_nrecs;
		} else {
			db_alloc = db_nrecs;
		}

    crossindex = (int *)calloc(db_alloc,sizeof(int));
    x = (double *)calloc(db_alloc,sizeof(double));
    y = (double *)calloc(db_alloc,sizeof(double));

    if ((crossindex == NULL) || (x == NULL) || (y == NULL)) {
      printf(
        "ERROR allocation failure: crossindex %p, x %p, or y %p\n",
        crossindex,
        x,
        y
      );
      exit(1);
    }

    /* Sort the table in increasing Julian Date for the movie option */
    qsort((void*)db_table,db_nrecs,sizeof(FILESTARIMAGEEXT),StarDateCompare);

    for (magnitudeIndex = 0; magnitudeIndex < db_nrecs; magnitudeIndex++) {
      pFileStarImage = &db_table[magnitudeIndex];
      webquality = pFileStarImage->quality;

      /* Handle pseudo-quality bits */
      if ((pFileStarImage->filestarimage.BFLAGS & (1 << FILTER_BFLAG_PSFSATURATED)) != 0) {
        webquality |= QUALITY_SATURATED;
      }

      if (((pFileStarImage->filestarimage.BFLAGS & (1 << FILTER_BFLAG_MAGDEP_MAGCOR)) == 0) &&
          ((pFileStarImage->quality & QUALITY_UNDETECTED) == 0)) {
        webquality |= QUALITY_NOMAGDEP;
      }

      if (GetFittedPlateScale(pFileStarImage->filestarimage.seriesId,pFileStarImage->filestarimage.plateNumber) >= PATROL_PLATE_SCALE) {
        webquality |= QUALITY_PATROL;
      } else {
        webquality |= QUALITY_NONPATROL;
      }

			if (pFileStarImage->filestarimage.ELLIPTICITY > MAX_ELLIPTICITY) {
        webquality |= QUALITY_TRAILED;
			}

      webseriesbit = 0;

      if ((pFileStarImage->filestarimage.seriesId > 0) &&
          (pFileStarImage->filestarimage.seriesId <= MAX_SERIES)) {
        pSeriesList = &seriesList[pFileStarImage->filestarimage.seriesId];
        webseriesbit = pSeriesList->bitMask;
      }

      pFileStarImage->quality = webquality;
      dateyear = jd2ep(pFileStarImage->filestarimage.Date);

      if (((pFileStarImage->filestarimage.AFLAGS & ~reqAFLAGSMASK) == 0) &&
          ((webquality & ~reqqualitymask) == 0)&&
          ((webseriesbit & ~reqseriesmask) == 0) &&
          (dateyear >= (reqstartdate-0.001)) &&
          (dateyear <= (reqenddate+0.001))) {
        x[ndata] = dateyear;
        y[ndata] = pFileStarImage->filestarimage.magcal_magdep;

        if ((pFileStarImage->quality & QUALITY_UNDETECTED) != 0) {
          y[ndata] = pFileStarImage->filestarimage.limiting_mag_local;
        }

        crossindex[ndata] = magnitudeIndex;
        ndata++;
      }
    }
  }

	if (foldflag) {
		for (idata = 0; idata < ndata; idata++) {
			pFileStarImage = &db_table[crossindex[idata]];
			foldfactor = (pFileStarImage->filestarimage.Date-foldcenter)/foldperiod;
			folddate = (pFileStarImage->filestarimage.Date-foldcenter) - (foldperiod*foldfactor);
			x[idata] = folddate/foldperiod;

			while (x[idata] < 0.0) {
				x[idata] += 1.0;
			}

			while (x[idata] > 1.0) {
				x[idata] -= 1.0;
			}

			x[idata+ndata] = x[idata] - 1.0;
			y[idata+ndata] = y[idata];
			crossindex[idata+ndata] = crossindex[idata];
		}

		ndata = ndata * 2;
	}

  {
    /* Here is Grzegorz Pojmanski's code from asas_raw_query.c to decode the click location */
    double xpixel; /* Renamed from x */
    double ypixel; /* Renamed from y */
    double d;
    double dmin;

    if (nextPointIndex < 0) {
      if(ix<x0 || ix>x1 || iy<y0 || iy>y1){
        printf("Plot coordinates file: %s<br>\n", store);
        printf(
          "You have clicked outside the valid area: "
          "%d %d\nsx=%d sy=%d\nx0,x1=%d %d fx=%lf fx0=%lf\ny0,y1=%d %d fy=%lf fy0=%lf<br>\n",
          ix, iy, sx, sy, x0, x1, fx, fx0, y0, y1, fy, fy0
        );
        exit(1);
      }

      xpixel = (ix-x0)/fx+fx0;
      ypixel = (iy-y0)/fy+fy0;
      dmin = 1.e32;
      bestidata = -1;

      for (idata = 0; idata < ndata; idata++) {
        if (y[idata] >= 29.99)
          y[idata]=upper_limit;

        xpixel = (x[idata]-fx0)*fx+x0;
        ypixel = ((y[idata]-fy0)*fy+y0);
        d = sqrt((xpixel-ix)*(xpixel-ix)+(ypixel-iy)*(ypixel-iy));

        if (d < dmin) {
          dmin = d;
          bestidata = idata;
        }
      }
    } else {
      bestidata = nextPointIndex;

      if (bestidata >= ndata) {
        printf("ERROR: selected point %d is greater than %d<br>\n", bestidata, ndata);
        exit(1);
      }
    }
  }

  if (bestidata >= 0) {
    magnitudeIndex = crossindex[bestidata];
    pFileStarImage = &db_table[magnitudeIndex];

    if ((pFileStarImage->quality & QUALITY_UNDETECTED) == 0) {
      averageRa = pFileStarImage->filestarimage.ra;
      averageDec = pFileStarImage->filestarimage.dec;
    }

    if (RefType == REF_TYPE_NONE) {
      averageRa = pFileStarImage->filestarimage.ra;
      averageDec = pFileStarImage->filestarimage.dec;
    }

    if (CheckAuthorization(regionFlag, averageRa, averageDec, &releaseField)) {
			double plateScale = GetFittedPlateScale(
        pFileStarImage->filestarimage.seriesId,
        pFileStarImage->filestarimage.plateNumber
      ); /* Degrees per pixel */

      if (strcmp(sizeunits,THUMBNAIL_ARCSEC) == 0) {
				if ((GetReleaseLevel() < RELEASE_LEVEL_ALL) && (imagesize > AUTHORIZE_RADIUS)) {
					imagesize = AUTHORIZE_RADIUS;
				}

        if (plateScale > 0) {
          /* New algorithm - fix the size at 300 pixels */
          scaleFactor = (3600. * THUMBNAIL_SCALED_PIXELS * plateScale)/(2.0 * imagesize * scaleFactor);
          sizePixels = ((2.0 * imagesize * scaleFactor) / 3600.0)/plateScale;
          if (sizePixels > THUMBNAIL_MAX_PIXELS) {
            printf("ERROR: Restricting size (1) to %d pixels <br />",THUMBNAIL_MAX_PIXELS);
            imagesize = (imagesize*THUMBNAIL_MAX_PIXELS)/sizePixels;
          }
        }

        imagearcmin = (2.0 * imagesize)/60.0;
      } else {
				if ((GetReleaseLevel() < RELEASE_LEVEL_ALL) && ((imagesize*plateScale*3600.) > AUTHORIZE_RADIUS)) {
					imagesize = AUTHORIZE_RADIUS/(plateScale*3600);
				}

        if (2*imagesize > THUMBNAIL_MAX_PIXELS) {
          printf("ERROR: Restricting size (2) to %d pixels <br />",THUMBNAIL_MAX_PIXELS);
          imagesize = THUMBNAIL_MAX_PIXELS/2;
        }

        imagearcmin = 2.0 * 60.0 * GetFittedPlateScale(
          pFileStarImage->filestarimage.seriesId,
          pFileStarImage->filestarimage.plateNumber
        ) * imagesize;
      }

      if (
        GetPhotPlate(
          pPhotConnection,
          GetSeriesString(pFileStarImage->filestarimage.seriesId, 1),
          pFileStarImage->filestarimage.plateNumber,
          pPhotPlates,
          "",
          1
        ) != 0
      ) {
        printf(
          "ERROR: Failed to get photplates record for %s%05d s%d\n",
          GetSeriesString(pFileStarImage->filestarimage.seriesId, 0),
          pFileStarImage->filestarimage.plateNumber,
          pFileStarImage->filestarimage.solutionNumber
        );
        exit(1);
      }

      if (
        GetMosaicInfo(
          pConnection,
          GetSeriesString(pFileStarImage->filestarimage.seriesId, 1),
          pFileStarImage->filestarimage.plateNumber,
          pPhotPlates->mosaicNumber,
          pFileStarImage->filestarimage.solutionNumber,
          pMosaic
        ) != 1
      ) {
        printf(
          "ERROR: GetMosaicInfo can not find mosaic %s%05d_%02d\n",
          GetSeriesString(pFileStarImage->filestarimage.seriesId, 0),
          pFileStarImage->filestarimage.plateNumber,
          pPhotPlates->mosaicNumber
        );
        exit(1);
      }

      // Try the new TNX archive bucket location before falling back to the old
      // raidNNN system.

      if (
        snprintf(
          filename,
          sizeof(filename),
          "/n/boslfs02/LABS/dasch_project/buckets/tnx.bucket/data/%02d/%02d/%s%05d_%02d_full.fits",
          pMosaic->plateNumber % 100,
          (pMosaic->plateNumber / 100) % 100,
          pMosaic->series,
          pMosaic->plateNumber,
          pMosaic->mosaicNumber
        ) >= sizeof(filename)
      ) {
        printf(
          "ERROR: buffer overflow for \"filename\"\n"
        );
        exit(1);
      }

      statResult = stat(filename, &statbuf);

      if (statResult < 0) {
        // Fall back.
        strcpy(filename, tmpdir);
        charPtr = strstr(filename, "tmp");

        if (charPtr == NULL) {
          printf("ERROR Invalid format for %s\n", tmpdir);
          exit(1);
        }

        *charPtr = '\0';

        if (pMosaic->rotation == 0) {
          rotsuffix[0] = '\0';
        } else {
          sprintf(rotsuffix, "r%d", pMosaic->rotation);
        }

        sprintf(
          filename2,
          "raid%03d/ExposureData/Mosaics/%s/%05d_%02d/%s%05d_%02d_01%sww_tnx.fit",
          pMosaic->diskLocation,
          pMosaic->series,
          pMosaic->plateNumber,
          pMosaic->mosaicNumber,
          pMosaic->series,
          pMosaic->plateNumber,
          pMosaic->mosaicNumber,
          rotsuffix
        );

        strcat(filename, filename2);
        statResult = stat(filename, &statbuf);
      }

      if (statResult == 0) {
        if (RefType == REF_TYPE_NONE) {
          printf("<b>Limiting Magnitude Data</b><br>\n");
        } else {
          printf("<b>Photometry data for %s (size %.0f\' x %.0f\')</b><br>\n",REF,imagearcmin,imagearcmin);
        }
      } else {
        if (m44release) {
          printf("<b>Preliminary photometry data for %s</b><br>\n",REF);
        } else {
          if (RefType == REF_TYPE_NONE) {
            printf("<b>Limiting Magnitude Data</b><br>\n");
          } else {
            printf("<b>Photometry data for %s</b><br>\n",REF);
          }
        }
      }

			if (pFileStarImage->filestarimage.spatial_bin == 0) {
				undetectedFlag = 1;
				result = CheckImageLocation(
          filename,
          averageRa,
          averageDec,
          (char *) &pFileStarImage->filestarimage.spatial_bin,
          &offscl,
          &edgeDistance
        );

				if (result < 0) {
					printf("ERROR: can not obtain image location for %s\n",filename);
					statResult = -1;
				}
			}

      dateyear = jd2ep(pFileStarImage->filestarimage.Date);

			if (foldflag != 0) {
				if (undetectedFlag == 0) {
					printf("%d Date %f Phase %f Mag %.2f\n", bestidata, dateyear, x[bestidata], y[bestidata]);
				} else {
					printf("%d Date %f Phase %f Limiting Mag %.2f\n", bestidata, dateyear, x[bestidata], y[bestidata]);
				}
			} else {
				if (undetectedFlag == 0) {
					printf("%d Date %f Mag %.2f\n", bestidata, x[bestidata], y[bestidata]);
				} else {
					printf("%d Date %f Limiting Mag %.2f\n", bestidata, x[bestidata], y[bestidata]);
				}
			}

      if ((RELEASE_EXPERIMENTAL != 0) && (catalogNumber == CATALOG_EXPERIMENTAL)) {
        printf("Calibration: %s\n", catalogText[pFileStarImage->filestarimage.catalogNumber]);
      }

      printf(
        "<a href='showplate.php?series=%s&plateNumber=%d&mosaicNumber=%d&filetype=FITS&ra=%f&dec=%f&system=J2000&radius=%f' "
        "target='results_frame'>Plate %s%05d</a>",
        GetSeriesString(pFileStarImage->filestarimage.seriesId, 1),
        pFileStarImage->filestarimage.plateNumber,
        pPhotPlates->mosaicNumber,
        averageRa,
        averageDec,
        imagearcmin / (2.0 * 60.0),
        GetSeriesString(pFileStarImage->filestarimage.seriesId, 1),
        pFileStarImage->filestarimage.plateNumber
      );

			if (offscl == 0) {
				printf(
          " Solution %d spatial_bin %d\n",
          pFileStarImage->filestarimage.solutionNumber,
          pFileStarImage->filestarimage.spatial_bin
        );

				if (undetectedFlag != 0) {
					printf("The image is %.1f cm from the plate edge\n\n", edgeDistance);
				} else {
					printf("Limiting magnitude %.2f\n\n", pFileStarImage->filestarimage.limiting_mag_local);
				}
			} else {
				printf(" Solution %d.\n\n Image is off edge of plate.\n\n", pFileStarImage->filestarimage.solutionNumber);
			}

      if (ndata > 1) {
        if (bestidata > 0) {
          printf("<button onClick='shownext(%d,%d); return true;'>Prev</button>", bestidata - 1, ndata);
        }

        if (bestidata < (ndata-1)) {
          if (automode == 0) {
            printf("<button id='autobutton' onClick='showauto(%d,%d); return true;'>Auto</button>", bestidata, ndata);
          } else {
            printf("<button id='autobutton' onClick='showauto(%d,%d); return true;'>Pause</button>", bestidata, ndata);
          }

          printf("<button onClick='shownext(%d,%d); return true;'>Next</button>", bestidata + 1, ndata);
        }

        printf("<br>");
      }

      if (statResult == 0) {
				if (offscl == 0) {
					printf(
            "<a href=\"extractimage.php?filetype=FITS&series=%s&plateNumber=%d&mosaicNumber=%d&ra=%f&dec=%f&system=J2000&radius=%f&scaleFactor=%f\">"
            "Get this image in FITS format.</a>\n",
            GetSeriesString(pFileStarImage->filestarimage.seriesId, 1),
            pFileStarImage->filestarimage.plateNumber,
            pPhotPlates->mosaicNumber,
            averageRa,
            averageDec,
            imagearcmin / (2.0 * 60.0),
            scaleFactor
          );
					printf(
            "<img src=\"extractimage.php?filetype=JPEG&series=%s&plateNumber=%d&mosaicNumber=%d&ra=%f&dec=%f&system=J2000&radius=%f&scaleFactor=%f\"><br />",
            GetSeriesString(pFileStarImage->filestarimage.seriesId, 1),
            pFileStarImage->filestarimage.plateNumber,
            pPhotPlates->mosaicNumber,
            averageRa,
            averageDec,
            imagearcmin / (2.0 * 60.0),
            scaleFactor
          );
					printf("<form action=\"lightcurve_point.php\" method=\"GET\">");
					printf("<input type=\"submit\" name=\"Replot\" value=\"Replot\" /> with radius ");
					printf("<input type=\"text\" name=\"imagesize\" id=\"imagesizewindow\" maxlength=5 size=5 value=\"%d\" />", imagesize);
					printf("<input type=\"radio\" name=\"sizeunits\" id=\"pixelsbutton\" value=\"pixels\"");

					if (strcmp(sizeunits, THUMBNAIL_PIXELS) == 0) {
						printf("checked=\"checked\"");
					}

					printf("/> Pixels ");
					printf("<input type=\"radio\" name=\"sizeunits\" id=\"pixelsbutton\" value=\"arcsec\"");

					if (strcmp(sizeunits, THUMBNAIL_ARCSEC) == 0) {
						printf("checked=\"checked\"");
					}

					printf("/> Arcsec ");

					if (GetReleaseLevel() < RELEASE_LEVEL_ALL) {
						printf("<br />The image radius must be less than %.0f arcsec\n", AUTHORIZE_RADIUS);
					}

					printf("<input type=\"hidden\" name=\"REF\" value=\"%s\" />", REF);
					printf("<input type=\"hidden\" name=\"source\" value=\"%s\" />", catalogText[catalogNumber]);
					printf("<input type=\"hidden\" name=\"tmpdir\" value=\"%s\" />", tmpdir2);
					printf("<input type=\"hidden\" name=\"gsc_bin_index\" value=\"%d\" />", gsc_bin_index);
					printf("<input type=\"hidden\" name=\"width\" value=\"%d\" />", plotWidth);
					printf("<input type=\"hidden\" name=\"height\" value=\"%d\" />", plotHeight);
					printf("<input type=\"hidden\" name=\"auto\" value=\"0\" />");
					printf("<input type=\"hidden\" name=\"coords\" value=\"");

					if ((ix == -1 && iy == -1) && nextPointIndex >= 0) {
						printf("%d\" />", nextPointIndex);
					} else {
						printf("X %d,%d\" />", ix, iy);
					}

					printf("</form><br>");
				}
      } else {
        printf("<b>WARNING: The file %s is not available</b><br>\n", filename);
      }

      if ((pFileStarImage->quality & QUALITY_UNDETECTED) != 0) {
        printf(
          "<b>No brightness estimate is available for this star on this plate.</b><br>"
          "The limiting magnitude is an average over one square degree and does not reflect "
          "local conditions: plate edges, missing emulsion, or sections of the plate where "
          "reliable magnitude measurements are not possible.<br>\n"
        );
      } else {
        sprintf(
          cmdStr,
          "%s/showflags -w %d >%s",
          binaries,
          pFileStarImage->filestarimage.AFLAGS,
          tempname
        );
        result = system(cmdStr);

        if (result != 0) {
          printf("ERROR: result %d decoding AFLAGS\n", result);
        }

        tempHandle = fopen(tempname, "rt");
        if (tempHandle == NULL) {
          printf("ERROR: Failed to open the temporary file %s\n", tempname);
          exit(1);
        } else if (verbose) {
          printf("Temporary file %s\n",tempname);
        }

        while (1) {
          inBuffer = fgets(inLine, MAX_BUFFER, tempHandle);
          if (inBuffer == NULL) {
            break;
          }

          printf("%s", inBuffer);
        }

        fclose(tempHandle);

        if (pFileStarImage->quality != 0) {
          printf("Plate Quality Flags %d<br />", pFileStarImage->quality);

          for (index = 0; index < webQualityTableSize; index++) {
            pQualityEntry = &webQualityMasks[index];

            if ((pQualityEntry->qualityMask & pFileStarImage->quality) != 0) {
              printf("Bit: %2d  %s<br />", index, pQualityEntry->qualityDescr);
            }
          }
        }

        sprintf(
          cmdStr,
          "%s/showflags -w -b %d >%s",
          binaries,
          pFileStarImage->filestarimage.BFLAGS,
          tempname
        );
        result = system(cmdStr);

        if (result != 0) {
          printf("ERROR: result %d decoding BFLAGS\n", result);
        }

        tempHandle = fopen(tempname, "rt");

        if (tempHandle == NULL) {
          printf("ERROR: Failed to open the temporary file %s\n", tempname);
          exit(1);
        } else if (verbose) {
          printf("Temporary file %s\n", tempname);
        }

        while (1) {
          inBuffer = fgets(inLine, MAX_BUFFER, tempHandle);

          if (inBuffer == NULL) {
            break;
          }

          printf("%s", inBuffer);
        }

        fclose(tempHandle);
      }

      printf(
        "<br>For more information about this plate, click on the above <q>Plate %s%05d</q> link.",
        GetSeriesString(pFileStarImage->filestarimage.seriesId, 1),
        pFileStarImage->filestarimage.plateNumber
      );
    } else {
      /* Authorization failure */
			printf("Sorry, these data are not available due to <a target=\"_blank\" href=\"https://dasch.cfa.harvard.edu/data-access/#restrictions\">temporary data access restrictions</a><br>\n");
		}
  } else {
    printf("ERROR: no object found near selected point\n");
    exit(1);
  }

  mysql_close(pPhotConnection);
  mysql_close(pConnection);

  if (crossindex != NULL) {
    free(crossindex);
  }

  if (x != NULL) {
    free(x);
  }

  if (y != NULL) {
    free(y);
  }

  if (db_table != NULL) {
    free(db_table);
  }

  return 0;
}

/* May 26, 2010 Edward J. Los - Adapted from asas_raw_query.c provided by
 *                              Grzegorz Pojmanski <gp@astrouw.edu.pl>  ASAS source author on Apr 26, 2010
 * Jul  2, 2010 Edward J. Los - Implement the "movie" feature.
 * Sep  6, 2010 Edward J. Los - Update for Solaris build
 * Jul 18, 2011 Edward J. Los - Change THUMBNAIL_SIZEUNITS to THUMBNAIL_PIXELS
 *                              Scale all constant-arcsec plots the same size.
 * Nov 11, 2011 Edward J. Los - Support magnitude-dependent calibration
 * Dec  5, 2011 Edward J. Los - Reduce scaleFactor by 50% in conjunction with doubling THUMBNAIL_IMAGESIZE
 * Dec  9, 2011 Edward J. Los - Filter on magnitude-dependent correction available
 * Mar 13, 2012 Edward J. Los - Support region-specific authentication
 * Mar 23, 2012 Edward J. Los - Support plot selection by series
 * Mar 11, 2013 Edward J. Los - Add RaPM, DecPM, ra_2 and dec_2 to the magnitude file (Version 5);
 * Mar 30, 2013 Edward J. Los - Support new data release mechanism
 * Jun 24, 2013 Edward J. Los - For points without magnitude measurements, check to see if the image is offscale,
 *                              otherwise find the spatial bin and the distance from the plate edge.
 * Jan 28, 2014 Edward J. Los - Add the limiting magnitude of all detections.
 * Jan 20, 2015 Edward J. Los - V6 data format: add A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex
 * Mar 30, 2018 Edward J. Los - Suport the "experimental" catalog combination of apass and gsc2.3.2 results.
 * Aug 13, 2018 Edward J. Los - Define enableRematch to optimize location of transients
 * Dec 13, 2019 Edward J. Los - Add a ds9 warning.
 * Mar 10, 2020 Edward J. Los - Update the ds9 warning
 */
