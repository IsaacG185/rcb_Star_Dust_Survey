// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* web_authorize.c
 *
 *  gcc -O0 -g -o web_authorize -I/usr/include/mysql web_authorize.c pipelineutils.a   -L /dasch/install/lib -L/usr/lib${lib64}/mysql  -lmysqlclient -lm -lwcs
cp web_authorize /dasch/install/bin
 *
 * Mar 12, 2010 Edward J. Los - Initial version
 *
 *  Usage: (not printed for security!)
 *  web_authorize -a <current authorization level, e.g. "guest" or "all">
 *
 *              PASSWORD AUTHORIZATION AND CHANGE REQUEST (login page)
 *                [-u <researcher name>]
 *                [-p <resarecher password>]
 *                [-n <new password>]
 *                [-t <retype password>]
 *
 *              LOCATION-BASED ACCESS REQUEST (all other pages) (returns "Allow" or "Deny");
 *                [-g <gsc_bin_index>]
 *                [-r <rightAscension> -d <declination>]  in decimal degrees
 *                [-m <mosaic reference>]
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"

#define MAX_USER_LENGTH 60
#define MAX_PASSWORD_LENGTH 50
#define MAX_BUFFER 80
#define MAX_EXPIRATION_DATE 25
#define MAX_SOLUTIONS 10
#include "mysql.h"
#include "pipelineutils.h"

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
struct WorldCoor *
wcskinit (
          int	naxis1,		/* Number of pixels along x-axis */
          int	naxis2,		/* Number of pixels along y-axis */
          char	*ctype1,	/* FITS WCS projection for axis 1 */
          char	*ctype2,	/* FITS WCS projection for axis 2 */
          double crpix1,
          double crpix2,	/* Reference pixel coordinates */
          double crval1,
          double crval2,	/* Coordinates at reference pixel in degrees */
          double *cd,		/* Rotation matrix, used if not NULL */
          double cdelt1,
          double cdelt2,	/* scale in degrees/pixel, ignored if cd is not NULL */
          double crota,		/* Rotation angle in degrees, ignored if cd is not NULL */
          int 	 equinox, /* Equinox of coordinates, 1950 and 2000 supported */
          double epoch);	/* Epoch of coordinates, used for FK4/FK5 conversion
                                 * no effect if 0 */

void
pix2wcs (

         struct WorldCoor *wcs,		/* World coordinate system structure */
         double	xpix, double ypix,	/* x and y image coordinates in pixels */
         double	*xpos,double *ypos);	/* RA and Dec in degrees (returned) */


int main(int argc,char *argv[])
{
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  char regionflag[MAX_BUFFER];
  char mosaicString[MAX_BUFFER];
  char researcher[MAX_USER_LENGTH];
  char userpassword[MAX_PASSWORD_LENGTH];
  char newpassword[MAX_PASSWORD_LENGTH];
  char retypepassword[MAX_PASSWORD_LENGTH];
  char **nextarg;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  char queryString[MAX_QUERY_STRING];
  int nvals;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int userId = -1;
  int tmpUserId;
  char tmpRegionFlag[MAX_REGIONFLAG];
  char tmpExpirationDate[MAX_EXPIRATION_DATE];
  char curTimeString[MAX_EXPIRATION_DATE];
  time_t tp;
  struct tm *ptr;
  int gsc_bin_index = -1;
  double rightAscension = -99.0;
  double declination = -99.0;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int rotation;
  MOSAIC mosaicStruct;
  PMOSAIC pMosaic = &mosaicStruct;
  int solutionNumber;
  struct WorldCoor *wcs;
  double xctr;
  double yctr;
  double cra;
  double cdec;

  double centerRa;
  double centerDec;
  int solutionCount;

  int recordNumber = 0;
  double cd[4];
	int releaseField;

  time(&tp);
  ptr = localtime(&tp);
  strftime(curTimeString,MAX_EXPIRATION_DATE, "%Y-%m-%dT%H-%M-%S", ptr);

  mosaicString[0] = 0;
  regionflag[0] = 0;
  researcher[0] = 0;
  userpassword[0] = 0;
  newpassword[0] = 0;
  retypepassword[0] = 0;
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


        case 'a': /* Current authorization level*/
        case 'A':
          if (argc > 1) {
            nextarg = argv;
            nextarg++;
            if (*nextarg[0] != '-') {
              argc--;
              argv++;
              if (strlen(*argv) >= MAX_BUFFER) {
                printf("ERROR: authorization level exceeds maximum buffer length\n");
                errorFlag = 1;
              } else {
                strcpy(regionflag,*argv);
              }
            }
          }
          break;

        case 'u': /* researcher name*/
        case 'U':
          if (argc > 1) {
            nextarg = argv;
            nextarg++;
            if (*nextarg[0] != '-') {
              argc--;
              argv++;
              if (strlen(*argv) >= MAX_USER_LENGTH) {
                printf("ERROR: user name exceeds maximum buffer length\n");
                errorFlag = 1;
              } else {
                strcpy(researcher,*argv);
              }
            }
          }
          break;

        case 'p': /* Password*/
        case 'P':
          if (argc > 1) {
            nextarg = argv;
            nextarg++;
            if (*nextarg[0] != '-') {
              argc--;
              argv++;
              if (strlen(*argv) >= MAX_PASSWORD_LENGTH) {
                printf("ERROR: password exceeds maximum buffer length\n");
                errorFlag = 1;
              } else {
                strcpy(userpassword,*argv);
              }
            }
          }
          break;

        case 'n': /* New password*/
        case 'N':
          if (argc > 1) {
            nextarg = argv;
            nextarg++;
            if (*nextarg[0] != '-') {
              argc--;
              argv++;
              if (strlen(*argv) >= MAX_PASSWORD_LENGTH) {
                printf("ERROR: new password exceeds maximum buffer length\n");
                errorFlag = 1;
              } else {
                strcpy(newpassword,*argv);
              }
            }
          }
          break;

        case 't': /* retype password*/
        case 'T':
          if (argc > 1) {
            nextarg = argv;
            nextarg++;
            if (*nextarg[0] != '-') {
              argc--;
              argv++;
              if (strlen(*argv) >= MAX_PASSWORD_LENGTH) {
                printf("ERROR: retyped exceeds maximum buffer length\n");
                errorFlag = 1;
              } else {
                strcpy(retypepassword,*argv);
              }
            }
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



        case 'r': /* rightAscension */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&rightAscension);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the rightAscension %s\n",*argv);
              errorFlag = 1;
            } else if (strstr(*argv,":") || strstr(*argv,",") || strstr(*argv," ")) {
              printf("ERROR: rightAscension %s is not in decimal degrees\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 'd': /* declination */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&declination);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the declination %s\n",*argv);
              errorFlag = 1;
            } else if (strstr(*argv,":") || strstr(*argv,",") || strstr(*argv," ")) {
              printf("ERROR: declination %s is not in decimal degrees\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 'm': /* mosaic string */
        case 'M':
          if (argc > 1) {
            nextarg = argv;
            nextarg++;
            if (*nextarg[0] != '-') {
              argc--;
              argv++;
              if (strlen(*argv) >= MAX_BUFFER) {
                printf("ERROR: mosaic string exceeds maximum buffer length\n");
                errorFlag = 1;
              } else {
                strcpy(mosaicString,*argv);

                if (ParseFilename2(mosaicString,series,&plateNumber) == 0) {
                  printf("ERROR: mosaic string %s can not be parsed\n",mosaicString);
                  errorFlag = 1;
                }

              }
            }
          }
          break;




        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }

      }

    }
  }
  if (regionflag[0] == 0) {
    strcpy(regionflag,"guest");
  }
  if ((newpassword[0] != 0)  && (strcmp(newpassword,retypepassword) != 0)) {
    printf("ERROR: Retyped password does not agree with new password\n");
    errorFlag = 1;
  }
  if (errorFlag== 1) {
    return(-1);
  }

  if (userpassword[0] != 0) {
    dasch_init_scandb(pConnection);

    sprintf(queryString,"SELECT userId,researcher,regionflag,expirationDate FROM users order by userId;",researcher);
    res = ExecuteQuery(pConnection,queryString);

    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          recordNumber++;
          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&tmpUserId);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for userId",nvals);
              exit(-1);
            }
          } else {
            printf("ERROR: missing userId\n");
            exit(-1);
          }

          if (sqlrow[1]) {
            if (strcmp(sqlrow[1],researcher) != 0) {
              continue;
            }
          } else {
            continue;
          }
          userId = tmpUserId;
          if (sqlrow[2]) {
            if (strlen(sqlrow[2]) < MAX_REGIONFLAG) {
              strcpy(tmpRegionFlag,sqlrow[2]);
            } else {
              printf("ERROR: regionflag is too large\n");
              exit(-1);
            }
          } else {
            printf("ERROR: missing regionflag\n");
            exit(-1);
          }
          if (sqlrow[3]) {
            if (strlen(sqlrow[3]) < MAX_EXPIRATION_DATE) {
              strcpy(tmpExpirationDate,sqlrow[3]);
              if (strcmp(curTimeString,tmpExpirationDate) > 0) {
                printf("guest\n");
                printf("Your account has expired\n");
                exit(0);
              }
            } else {
              printf("ERROR: expiration date is too large\n");
              exit(-1);
            }
          } else {
            tmpExpirationDate[0] = 0;
          }



        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
        exit(-1);
      }
    }
    if (userId < 0) {
      printf("guest\n");
      printf("No such account\n");
      exit(0);
    }
    tmpUserId = -1;
    sprintf(queryString,"SELECT userId FROM passwords where userId = %d and password = PASSWORD('%s');",userId,userpassword);
    res = ExecuteQuery(pConnection,queryString);
    if (!res) {
      res_ptr = mysql_store_result(pConnection);
      if (res_ptr) {
        while ((sqlrow = mysql_fetch_row(res_ptr))) {
          if (sqlrow[0]) {
            nvals = sscanf(sqlrow[0],"%d",&tmpUserId);
            if (nvals != 1) {
              printf("ERROR: nvals is %d for userId",nvals);
              exit(-1);
            }
          } else {
            printf("ERROR: missing userId\n");
            exit(-1);
          }

        }
        mysql_free_result(res_ptr);
      } else {
        printf("ERROR: mysql_store_result failed in pipelineutils.c line %d\n",__LINE__);
        exit(-1);
      }
    }

    if (tmpUserId != userId) {
      printf("guest\n");
      printf("Invalid Password\n");
      exit(0);
    }
    sprintf(queryString,"UPDATE passwords set accessDate = '%s' where userId = %d and password = PASSWORD('%s');",curTimeString,userId,userpassword);
    res = ExecuteQuery(pConnection,queryString);
    if (res) {
      exit(-1);
    }
    if (newpassword[0] != 0) {
      sprintf(queryString,"UPDATE passwords set password = PASSWORD('%s'),changeDate = '%s' where userId = %d and password = PASSWORD('%s');",newpassword,curTimeString,userId,userpassword);
      res = ExecuteQuery(pConnection,queryString);
      if (res) {
        exit(-1);
      }
      printf("%s\n",tmpRegionFlag);
      printf("Password changed\n");

    } else {
      printf("%s\n",tmpRegionFlag);
      printf("Authorization Succeeded\n");
    }


    mysql_close(pConnection);
    exit(0);



  }


  if (gsc_bin_index > 0) {
    /* Return status based on the current authorization flag and bin location */
    if ((GetBinCenter(pGscBin,gsc_bin_index,&centerRa,&centerDec,"web_plot") != 0) ||
        (CheckAuthorization(regionflag,centerRa,centerDec,&releaseField) == 0)) {
      /* Authorization failure */
      printf("%s\n",regionflag);
      printf("Deny\n");
      exit(0);
    } else {
      printf("%s\n",regionflag);
      printf("Allow\n");
      exit(0);
    }
  } else if ((rightAscension >= 0.0) || (declination >= -90.0)) {
    if (CheckAuthorization(regionflag,rightAscension,declination,&releaseField)) {
      printf("%s\n",regionflag);
      printf("Allow\n");
      exit(0);
    } else {
      printf("%s\n",regionflag);
      printf("Deny\n");
      exit(0);
    }

  } else if (mosaicString[0] != 0) {
    dasch_init_scandb(pConnection);

    if (SelectBestMosaic(pConnection,series,plateNumber,&mosaicNumber,&rotation) == 0) {
      /* Here there are no mosaics.  Say "Deny" to be safe, but there is nothing to deny */
      printf("%s\n",regionflag);
      printf("Deny\n");
      mysql_close(pConnection);
      exit(0);
    } else {
      for (solutionCount = 0; solutionCount < MAX_SOLUTIONS; solutionCount++) {
        if (GetMosaicInfo(pConnection,series,plateNumber,mosaicNumber,solutionCount,pMosaic) == 0) {
          break;
        }
        if (strstr(pMosaic->ctype1,"DEC")) {
          char tmpPtr[MAX_CTYPE_STRING];
          double dtmp;
          int itmp;
          strcpy(tmpPtr,pMosaic->ctype2);
          strcpy(pMosaic->ctype2,pMosaic->ctype1);
          strcpy(pMosaic->ctype1,tmpPtr);
          dtmp = pMosaic->crval1;
          pMosaic->crval1 = pMosaic->crval2;
          pMosaic->crval2 = dtmp;



          cd[0] = pMosaic->cd2_1;
          cd[1] = pMosaic->cd2_2;
          cd[2] = pMosaic->cd1_1;
          cd[3] = pMosaic->cd1_2;

        } else {

          cd[0] = pMosaic->cd1_1;
          cd[1] = pMosaic->cd1_2;
          cd[2] = pMosaic->cd2_1;
          cd[3] = pMosaic->cd2_2;

        }



        wcs = wcskinit(pMosaic->naxis1,
                       pMosaic->naxis2,
                       pMosaic->ctype1,
                       pMosaic->ctype2,
                       pMosaic->crpix1,
                       pMosaic->crpix2,
                       pMosaic->crval1,
                       pMosaic->crval2,
                       cd,
                       0,  /* cdelt1 */
                       0,  /* cdelt2 */
                       0,  /* crota */
                       2000, /* equinox */
                       0);   /* epoch */
        if (wcs) {
          xctr = 0.5 + (0.5 *pMosaic->naxis1);
          yctr = 0.5 + (0.5 *pMosaic->naxis2);
          pix2wcs(wcs,xctr,yctr,&cra,&cdec);
          wcsfree(wcs);
          if (CheckAuthorization(regionflag,cra,cdec,&releaseField) == 0) {
            printf("%s\n",regionflag);
            printf("Deny\n");
            mysql_close(pConnection);
            exit(0);
          }


        } else {
          printf("ERROR creating wcs structure for %s%05d\n",series,mosaicNumber);
          mysql_close(pConnection);
          exit(-1);

        }

      }
      printf("%s\n",regionflag);
      printf("Allow\n");
      mysql_close(pConnection);
      exit(0);

    }
    mysql_close(pConnection);

  }



  printf("guest\n");
  printf("You are currently not logged into the website\n");
  return(0);
}
