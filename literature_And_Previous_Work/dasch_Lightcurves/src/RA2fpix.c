/* RA2fpix.c
 *
 *  RA2fpix by David Koch  D.Koch@nasa.gov, version of Wed, 17 Jun 2009
 *  Feb  5, 2009 Edward J. Los - modified for pipeline use
 *  Jul  4, 2011 Edward J. Los - correct infinite loops waiting for console input
 *
 * RA2fpix.c 
 * Purpose is to have a single program that converts RA both to "2" pixels as well as 
 from "f" pixels 
 * version 2.1 added FGS and correct values for chips 7 21 and 29 
 *
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm -lmysqlclient testkepler.c pipelineutils.a -ltable -lutil -lwcs -o testkepler 
 *
 *  ./testkepler -v -i /dasch/junk/id_all_kepler.db  -o /dasch/junk/kepler.tmp
 *
 * 
   column -b -i kepler.tmp ra dec > los.tmp
   psxy  los.tmp -R279/302/36/53 -Y4i -JX6i/5.5i -P -B5:"RA":/5:"DEC":WSne -Sc0.03  > los.ps
   gv los.ps &
 *
 * Jul  1, 2000 Edward J. Los - Initial version
 */


#include <math.h>
#if 0
#include "table.h"
#include "pipelineutils.h"
#endif
#include <stdio.h>
#include <string.h>
#include <time.h>
/* #include <ps.h> /* plotting postscript */
#include <stdlib.h> /* needed for exit */


#define MAX_BUFFER 100
#define nzt 35 /* number of angle pairs in opt_xf */

void sep(double alp1,double del1,double alp2,double del2,double *delta,int units); /* added 2/22/05 */
void RA2(double tra, double tdec, double *lngr, double *latm, int *chip_n, int *chn_m, double *p_row, double *p_coln);
void RAF(int chn_n, double p_row,double p_coln, double *tra, double *tdec);
int RA2fpix(double inputRA, double inputDec,int *pChn_n,double *pRow,double *pColn);
	
double DCM11,DCM12,DCM13,DCM21,DCM22,DCM23,DCM31,DCM32,DCM33; /* direction cosine matrix */	

int chn_num[10][10] ={0,0,56,55,36,35,16,15, 0,0,		/* array of channel numbers */
                      0,87,53,54,33,34,13,14,85,0,		/* added FGS 87 & 85  3/4/05 */
                      75,74,60,59,40,39,17,20, 1,4,
                      76,73,57,58,37,38,18,19, 2,3,
                      79,78,63,62,44,43,21,24, 5,8,		/* note layout of matrix mimics  */
                      80,77,64,61,41,42,22,23, 6,7,
                      83,82,67,66,46,45,26,25, 9,12,		/* actual physical layout of chn numbers */
                      84,81,68,65,47,48,27,28,10,11,		/* on FPA */
                      0,88,70,69,50,49,30,29,86,0,		/* added FGS 88 & 86  3/4/05 */
                      0, 0,71,72,51,52,31,32, 0,0	};

double plate_scale=3.9753235, plate_scale_fgs ; /* arc sec per pixel 2/24/05 */
								
double DCM11c[46], DCM12c[46], DCM13c[46],
  DCM21c[46], DCM22c[46], DCM23c[46], 
  DCM31c[46], DCM32c[46], DCM33c[46];	
	
double chip_offset[46][3]={ /* chip offsets
                               first term is rotation in degrees added to last rotation in chip_trans[][]
                               second term is half gap width from module center in pixels, a row offset
                               last term is a column offset in pixels
                            */
  /* values from RA2pixCoef.out */
  -0.00540,  40.38896,  2.43919,  /* chip   1 */
  -0.02789,  37.47266,  -2.29810,  /* chip   2 */
  -0.00677,  37.74270,  1.69278,  /* chip   3 */
  0.00203,  37.08065,  -0.61181,  /* chip   4 */
  0.04963,  37.83317,  0.67468,  /* chip   5 */
  -0.05890,  39.88030,  0.60112,  /* chip   6 */
  -0.01208,  37.02522,  1.34406,  /* chip   7 */
  -0.01600,  39.81689,  -0.89545,  /* chip   8 */
  0.03726,  39.79497,  1.39803,  /* chip   9 */
  -0.01928,  35.37856,  -0.14501,  /* chip  10 */
  0.00814,  38.00031,  1.31946,  /* chip  11 */
  -0.01825,  37.89786,  -0.87843,  /* chip  12 */
  -0.05113,  38.70468,  1.21300,  /* chip  13 */
  0.00754,  37.09967,  -0.64806,  /* chip  14 */
  0.02746,  39.87765,  0.37930,  /* chip  15 */
  -0.03317,  36.03914,  -0.22458,  /* chip  16 */
  -0.01514,  39.13289,  0.16309,  /* chip  17 */
  -0.00731,  38.04260,  -0.55602,  /* chip  18 */
  -0.00407,  37.05139,  0.25807,  /* chip  19 */
  -0.01971,  37.75357,  0.36864,  /* chip  20 */
  0.04503,  37.55200,  0.06207,  /* chip  21 */
  0.03818,  37.82782,  0.05873,  /* chip  22 */
  -0.04545,  39.63746,  0.77143,  /* chip  23 */
  -0.00438,  37.94519,  -0.93808,  /* chip  24 */
  0.00852,  38.28699,  1.31147,  /* chip  25 */
  -0.00908,  38.07401,  -1.12547,  /* chip  26 */
  0.00533,  37.37255,  2.31620,  /* chip  27 */
  0.00337,  37.09931,  -1.05059,  /* chip  28 */
  -0.00260,  37.88098,  1.20230,  /* chip  29 */
  0.00327,  37.33139,  -0.42545,  /* chip  30 */
  0.02050,  38.75998,  -0.56909,  /* chip  31 */
  -0.02936,  38.24879,  -0.01055,  /* chip  32 */
  0.01456,  39.91342,  -0.23958,  /* chip  33 */
  -0.00597,  38.31758,  -0.49553,  /* chip  34 */
  -0.00872,  37.90025,  0.38379,  /* chip  35 */
  0.01572,  39.84887,  -1.32037,  /* chip  36 */
  0.02500,  36.97141,  0.17216,  /* chip  37 */
  0.01718,  37.20423,  -0.01043,  /* chip  38 */
  0.02684,  38.94127,  0.60935,  /* chip  39 */
  -0.00249,  37.41979,  -1.03399,  /* chip  40 */
  0.03646,  39.72337,  0.67074,  /* chip  41 */
  0.00517,  38.37851,  -1.39885,  /* chip  42 */
  0.11444,  -5.92846,  -0.61846,  /* chip  43 */
  0.10827,  -3.75692,  2.89846,  /* chip  44 */
  -0.12767,  -2.08077,  2.17615,  /* chip  45 */
  0.11151,  -6.13615,  -2.54154  /* chip  46 */
};

double /* changed name so didn't have to comment out 46 lines */ foo, chip_offset_orig[46][3]={ /* chip offsets
                                                                                                   first term is rotation in degrees added to last rotation in chip_trans[][]
                                                                                                   second term is half gap width from module center in pixels, a row offset
                                                                                                   last term is a column offset in pixels
                                                                                                */
  0.00000, 38.884, 0.0, /* chip  1 */
  0.00000, 38.884, 0.0, /* chip  2 */
  0.00000, 38.884, 0.0, /* chip  3 */
  0.00000, 38.884, 0.0, /* chip  4 */
  0.00000, 38.884, 0.0, /* chip  5 */
  0.00000, 38.884, 0.0, /* chip  6 */
  0.00000, 38.884, 0.0, /* chip  7 */
  0.00000, 38.884, 0.0, /* chip  8 */
  0.00000, 38.884, 0.0, /* chip  9 */
  0.00000, 38.884, 0.0, /* chip 10 */
  0.00000, 38.884, 0.0, /* chip 11 */
  0.00000, 38.884, 0.0, /* chip 12 */
  0.00000, 38.884, 0.0, /* chip 13 */
  0.00000, 38.884, 0.0, /* chip 14 */
  0.00000, 38.884, 0.0, /* chip 15 */
  0.00000, 38.884, 0.0, /* chip 16 */
  0.00000, 38.884, 0.0, /* chip 17 */
  0.00000, 38.884, 0.0, /* chip 18 */
  0.00000, 38.884, 0.0, /* chip 19 */
  0.00000, 38.884, 0.0, /* chip 20 */
  0.00000, 38.884, 0.0, /* chip 21 */
  0.00000, 38.884, 0.0, /* chip 22 */
  0.00000, 38.884, 0.0, /* chip 23 */
  0.00000, 38.884, 0.0, /* chip 24 */
  0.00000, 38.884, 0.0, /* chip 25 */
  0.00000, 38.884, 0.0, /* chip 26 */
  0.00000, 38.884, 0.0, /* chip 27 */
  0.00000, 38.884, 0.0, /* chip 28 */
  0.00000, 38.884, 0.0, /* chip 29 */
  0.00000, 38.884, 0.0, /* chip 30 */
  0.00000, 38.884, 0.0, /* chip 31 */
  0.00000, 38.884, 0.0, /* chip 32 */
  0.00000, 38.884, 0.0, /* chip 33 */
  0.00000, 38.884, 0.0, /* chip 34 */
  0.00000, 38.884, 0.0, /* chip 35 */
  0.00000, 38.884, 0.0, /* chip 36 */
  0.00000, 38.884, 0.0, /* chip 37 */
  0.00000, 38.884, 0.0, /* chip 38 */
  0.00000, 38.884, 0.0, /* chip 39 */
  0.00000, 38.884, 0.0, /* chip 40 */
  0.00000, 38.884, 0.0, /* chip 41 */
  0.00000, 38.884, 0.0, /* chip 42 */
  0.00000,  0.000, 0.0, /* fgs 1 added 3/4/05 */
  0.00000,  0.000, 0.0, /* fgs 2 added 3/4/05 */
  0.00000,  0.000, 0.0, /* fgs 3 added 3/4/05 */
  0.00000,  0.000, 0.0  /* fgs 4 added 3/4/05 */
};




/******************  convert from row and column to chip lng and lat *************/

void RAF(int chn_n, double p_row,double p_coln, double *trax, double *tdecx)
{
  double pi=3.141592654, tpi=2.*pi, pi2=pi/2., dtr=pi/180., rtd=180./pi;
  double tra, tdec;	
  int chip_n, quad;
  double lngr, latp, latm, lngm, lpm, mpm, npm, lp, mp, np;
  double cosa, cosb, cosg;
	
  if(chn_n<85) { /* do sci chn */
    chip_n=(chn_n+1)/2-1;  /* chip number index, -1 since DCMc array index is 0-41 */
    quad=chn_n % 4; /* need to find quadrant in module */
    lngr=plate_scale*(1024.0 + chip_offset[chip_n][1] - p_row);
    latp= 1100-p_coln;
    if(quad==0 || quad==2) latp=-latp;
    latm=(latp + chip_offset[chip_n][2])*plate_scale; /* in arc sec */
  } else { /* do FGS chips */
    chip_n=chn_n-43; /* FGS chip numbers */
    lngr=p_row*plate_scale_fgs+chip_offset[chip_n][1]; 
    /* chn_n-42 results in indices 43-46 for the FGS */
    latm=p_coln*plate_scale_fgs+chip_offset[chip_n][2]; /* added 3/4/05 */
  }
  latm=latm/3600./rtd; /*latm in radians */
  lngm=lngr/cos(latm); /* correct for cos effect going from spherical to rectangular coor 
                          one deg of long at one deg of lat is smaller than one deg at zero lat
                          by cos(lat) , amounts to 1/2 arc sec=1/8 pix */
  lngm=lngm/3600./rtd; /* lngm in radians */
  /* get direction cosines */
  lpm=cos(lngm)*cos(latm);
  mpm=sin(lngm)*cos(latm);
  npm=sin(latm);
			
  /* transform from chip coor to FPA coor 
     Do inverse of above transform (swap matrix row/coln indices) */
  lp=DCM11c[chip_n]*lpm+DCM21c[chip_n]*mpm+DCM31c[chip_n]*npm;
  mp=DCM12c[chip_n]*lpm+DCM22c[chip_n]*mpm+DCM32c[chip_n]*npm;
  np=DCM13c[chip_n]*lpm+DCM23c[chip_n]*mpm+DCM33c[chip_n]*npm;
	
  /* Transform from FPA to RA and Dec Again use inverse of transform matrix */
  cosa=DCM11*lp+DCM21*mp+DCM31*np;
  cosb=DCM12*lp+DCM22*mp+DCM32*np;
  cosg=DCM13*lp+DCM23*mp+DCM33*np;
		
  /* convert dir cosines to equatorial coor system */
  tdec=asin(cosg)*rtd; /* transformed Dec in deg */
  tra=atan2(cosb,cosa)*rtd; /* transformed RA in deg */
  if(tra<0.) tra=tra+360.; 
		
  *trax=tra;
  *tdecx=tdec;
  return;
} /* end RAF */	
			

/************************ convert tra and tdec to direction cosines *************/

void RA2(double tra, double tdec, double *lngrx, double *latmx, int *chip_m, int *chn_m, double *p_row, double *p_coln)
{
  double pi=3.141592654, tpi=2.*pi, pi2=pi/2., dtr=pi/180., rtd=180./pi;
  double cosa,cosb,cosg,lp,mp,np, lpm,mpm,npm, lng,lat,latp, lngm, latm, lngr;
  int chn_i, chn_j, chn_n, chip_n, bad=-99;
	
  tra=tra*dtr; /* tra in radians */
  tdec=tdec*dtr; /* tdec in radians */
  /**		printf("tra=%f tdec=%f\n", tra,tdec); /**/
  cosa=cos(tra)*cos(tdec);
  cosb=sin(tra)*cos(tdec);
  cosg=sin(tdec);

  /* now do coor transformation get direction cosines in FPA coor*/
  lp=DCM11*cosa+DCM12*cosb+DCM13*cosg;
  mp=DCM21*cosa+DCM22*cosb+DCM23*cosg;
  np=DCM31*cosa+DCM32*cosb+DCM33*cosg;
  /*		printf("lp=%f mp=%f np=%f\n", lp,mp,np); */
		
  /* convert dir cosines to longitude and lat in FPA coor system */
  lat=asin(np)*rtd; /* transformed lat +Z' in deg */
  lng=atan2(mp,lp)*rtd; /* transformed long +Y' in deg */
  /*	fprintf(stderr,"FPA: lng=%f lat=%f  ",lng,lat); /* */
		
  /* find which chn this falls onto */
  chn_i=floor(lat/1.430)+5;
  chn_j=floor(lng/1.430)+5;
  /*	fprintf(stderr,"lat i=%d lng j=%d\n",chn_i,chn_j);*/
  if(chn_i<0 || chn_i>9 || chn_j<0 || chn_j>9) {
    chn_n=bad; /* coor well out of FOV */	
  } else {
    chn_n=chn_num[chn_i][chn_j]; /* channel number */
  }
  if (chn_n>0 && chn_n<85)  /* added if to deal with fgs case Fixed lower lim 8/25/05 */
    chip_n=(chn_n+1)/2-1;  /* chip number index, -1 since DCMc array index is 0-41 */
  else if (chn_n>84 && chn_n<89)  /* FGS chip Fixed lower lim 8/25/05 */
    chip_n=chn_n-43; /* added for fgs */
  else {/* chn_n==0 */
    chip_n=bad;
    *p_row=bad;
    *p_coln=bad;
  }
  /**	fprintf(stderr,"chn#=%d chip#=%d  Must be >=0 and <=41 \n",chn_n,chip_n); /**/
  /* can now transform to module coordinates Use direction cosine in FPA coor*/
  /* now do transformation to module chip coor*/

  if(chip_n>-1 && chip_n<46) { /* coor falls on a chip Fixed lower lim 8/25/05 */
    lpm=DCM11c[chip_n]*lp+DCM12c[chip_n]*mp+DCM13c[chip_n]*np;
    mpm=DCM21c[chip_n]*lp+DCM22c[chip_n]*mp+DCM23c[chip_n]*np;
    npm=DCM31c[chip_n]*lp+DCM32c[chip_n]*mp+DCM33c[chip_n]*np;
    /* define chip coor as: rotation about the center of the module(field flattener lens) &
       angular row and column from this center
       then column 1100 is angle zero and decreases up and down with increasing angle
       towards readout amp on each corner
       and row 1024 starts after a gap of 39 pixels decreasing with increasing angle */
			
    latm=asin(npm);/* transformed lat +Z' to chip coor in radians */
    lngm=atan2(mpm,lpm)*rtd*3600.; /* transformed long +Y' to chip coor in arc sec */
    lngr=lngm*cos(latm); /* correct for cos effect going from spherical to rectangular coor 
                            one deg of long at one deg of lat is smaller than one deg at zero lat
                            by cos(lat) , amounts to 1/2 arc sec=1/8 pix */
    latm=latm*rtd*3600.; /* latm in arc sec */
    /*	fprintf(stderr,"Chip: lng=%f lat=%f \n",lngm/3600.,latm/3600.);*/
    /* now convert to row and column */
    if(chn_n>0 && chn_n<85) { /* do for sci CCDs added 3/4/05 fixed lower lim 8/25/05 */
      *p_row=1024.0-lngr/plate_scale + chip_offset[chip_n][1];
      latp=latm/plate_scale - chip_offset[chip_n][2]; /* +/-latitude in pixels on chip */
      if(latp>=0.0)
        *p_coln=1100.0-latp;
      else
        *p_coln=1100.0+latp;
				
      /* get correct chn_n, since initial chn_n was a close guess */
      if(latp<0.0) /* then falls on top half of chip */
        chn_n=(int)((chn_n+1)/2) *2;
      /* if +latp side of chip chn_n is even*/
      else
        chn_n=(int)((chn_n-1)/2) *2+1;
      /*   if -latp side of chip chn_n is odd */		
    } else { /* do for fgs added 3/4/05 */
      *p_row=lngr/plate_scale_fgs+chip_offset[chip_n][1]; 
      /* chn_n-42 results in indices 43-46 for the FGS */
      *p_coln=latm/plate_scale_fgs+chip_offset[chip_n][2]; /* added 3/4/05 */
    }
  }
  /* finished conversion to chn_n row and column */
  *chn_m=chn_n;
  *chip_m=chip_n;
  *lngrx=lngr;
  *latmx=latm;
  return;
} /* end RA2 */


/********************************   Calculates angular separation *************************/
	
void sep(double alp1,double del1,double alp2,double del2, double *delta,int units)
{
	
  double pi=3.141592654, tpi=2.*pi, dtr=pi/180., rtd=180./pi;
  double cossep,conv;
	
  if(units==0) {/* convert deg to radians */
    alp1=alp1*dtr;
    del1=del1*dtr; 
    alp2=alp2*dtr;
    del2=del2*dtr;
    conv=180./pi;
  } else 
    conv=1.;
		
  cossep=sin(del1)*sin(del2)+cos(del1)*cos(del2)*cos(alp1-alp2);
	
  if(cossep>=1. && cossep<1.000001) 
    *delta=0.;
  else if(cossep<=-1. && cossep>-1.000001) 
    *delta=pi*conv;
  else if(cossep<1. && cossep>-1.)
    *delta=acos(cossep)*conv;
  else {
    *delta=9999999999.;
    printf("Error in coor sep calculation Halt.\n");
    exit(0 );
  }
  /*printf("del=%f\n",delta);*/
  return;
}

		
		
		


            

int RA2fpix(double inputRA, double inputDec,int *pChn_n,double *pRow,double *pColn)

{
  double pi=3.141592654, tpi=2.*pi, pi2=pi/2., dtr=pi/180., rtd=180./pi;
  FILE *fp1; /* fp1 has the output chn row and column for option s */
  FILE *fp2; /* fp2 has the input RA and dec in deg to convert for option f & s */
  char file2[80], fileN[100];  /* file name and junk data to read over */

  /* Use array opt_xf to find transformation angles so corners come out to the deired value
     Invoke this array using option z 
     when using these values, program will convert codeV angles given here to azim and elev, 
     which changes tan of second angle value by cos of first 
     These data are from Ball-Angela Sparks from DFM.PHT.019 3/11/05 , FGS updated 6/24/05 
     when adding to this array you need to change value of nzt in #define */
  double opt_xf[ nzt ][2]={	
    0.0,       0.0, 		/* zoom 1 */
    -1.173932 , 1.2149095,		/* data are for zoom 1-6 */ 
    1.173932,  1.2149095,	 	/* except for z6 coor are center, upper left, upper right, */
    -1.173932, -1.2149095,		/* lower left, lower right */
    1.173932, -1.2149095,		/* changed values for module 13 3/4/05 */
    0.0,       2.8594335,	/* zoom 2 */	/* added 22 Feb 05 */
    -1.2176889, 4.0333823,		
    1.2176889, 4.0333823,		
    -1.215182 , 1.6854851,		
    1.215182 , 1.6854851,		
    0.0,       5.7189145, /* zoom 3 */
    -1.2235547, 6.89289,
    1.2235547, 6.89289,
    -1.21851  , 4.54494,
    1.21851  , 4.54494,
    -2.859409,  2.859409, /* zoom 4 */ /* needed to change sign on X for this zoom */
    -1.684259,  4.073138,
    -4.0375935, 4.0785615,
    -1.683729,  1.642762,
    -4.032057,  1.6431795,
    2.8593305, 5.7188497, /* zoom 5 */
    1.638719 , 6.8900115,
    4.0857918, 6.8985335,
    1.637894 , 4.541664,
    4.07491  , 4.5452094,
    -2.8593305, 5.7188497, /* zoom 5a Create by reflecting 5 about x axis */
    -1.638719 , 6.8900115,
    -4.0857918, 6.8985335,
    -1.637894 , 4.541664,
    -4.07491  , 4.5452094,
    -4.77637,   4.66802,   /* loc 2  FGS zoom 6  added 3/4/05  refined positions 6/24/05 */
    -4.91982  , 4.80933  , /* loc 6  changed sign on x for FGS just like for z 4 */
    -4.91885  , 4.52759,   /* loc 7  FGS data from DFM.PHT.019a 6/23/05 */
    -4.63383  , 4.8084   , /* loc 9 */
    -4.63303  , 4.5268     /* loc 10  changed y from 4.52759 to 4.5268 since appeared to 
                              be off by 1.5 pix 6/24/05 */
  };
  int nz; /* note in later loop test for nz equal to number of pairs in above array */
  int zi,ci,cc,cx; /*added 2/22/05*/
	
  /***************  TRANFORMATION VALUES *******/

  /*double rot_3=290.6835251;*/ double rot_3=290.66666667; /* 3-rot corresponds to FOV at RA 19h 22m 40s changed 9/8/04 */
  /*double rot_2=-44.48147667; */ double rot_2=-44.5; 	/* 2-rot corresponds to FOV Dec at 44d 30m 00s Note minus sign
    since +dec corresponds to -rot changed 9/8/04 */
								
								
								
  double rot_1; /* depends on season */
  int season_int, convert_it;
  double first_roll=110.; /* Spacecraft Roll angle for summer. Replaces use of season 9/10/04 */
  /* the season determines the last rotation angle */
  /* summer =0, fall =1, winter =2, spring =3, 
     changed 9/8/04 */
  double clocking_angle=13.; /* add to spacecraft roll and n*90 deg to get FPAA roll=rot_1 9/10/04 */

  int i, check;
  char ang; /* units for entering angle d or t 
               OR if z then uses angles in opt_xf for input and rot_1 _2 _3 set to zero*/
	
  double srac,crac,sdec,cdec,srotc,crotc;
	
  double chip_trans[46][3]={ 
    /* starting values are obtained from DFM.PHT.019a for the field=1 of each module
       then they are adjusted so that the corners of the modules match the 
       field positions 5,6,8,9
       DFM gives angles as azim,azim (u,v) need azim, elev (u,w)
       Convert second azim v to w using tan w=tan v * cos u Done in z option code below
						
       These are the 3-2-1 transform values for the 46 chips= (channel+1)/2 
       transform is from center of FPA to center of the module for each chip (DCA) 
       transformed coor will have valid long values negative 0 to -1.1321 deg
       lat will be +/- 1.2161 deg 
       lng corresponds to rows, lat  to columns*/
						
    /* z'-rot   Y'-rot    X'-rot 
       Y-angle  X'-angle   module_rotation 
       where tan X'= tan X * cos Y  since CodeV uses (azim,azim) coor and need (azim,elev) coor   */
    5.71881,   2.84513, 180.14080,    /* mod  2, chip 1  chn 1 & 2    zoom 5 */
    5.71881,   2.84513,   0.14080,    /* mod  2, chip 2  chn 3 & 4    */ 
    5.71893,   0.00000, 180.00000,    /* mod  3, chip 3  chn 5 & 6    zoom 3 */
    5.71893,   0.00000,   0.00000,    /* mod  3, chip 4  chn 7 & 8    */
    5.71881,  -2.84513, 179.85920,    /* mod  4, chip 5  chn 9 & 10   zoom 5 */
    5.71881,  -2.84513,  -0.14080,    /* mod  4, chip 6  chn 11 & 12  */
    /* get next by doing 90 degree rotation of zoom coor */
    /* sin Y''=sin Y * cos X' , sin X''=sin X' / cos Y''  provides 90 coor rotation
     */
    2.85934,   5.71174,  90.14400,    /* mod  6, chip 7  chn 13 & 14  zoom 5+90deg rot */
    2.85934,   5.71174, 270.14400,    /* mod  6, chip 8  chn 15 & 16  */
    2.85942,   2.85586, 180.07160,    /* mod  7, chip 9  chn 17 & 18  zoom 4 */
    2.85942,   2.85586,   0.07160,    /* mod  7, chip 10 chn 19 & 20  */
    2.85945,   0.00000, 180.00000,    /* mod  8, chip 11 chn 21 & 22  zoom 2 */
    2.85945,   0.00000,   0.00000,    /* mod  8, chip 12 chn 23 & 24  */
    2.85942,  -2.85587, 269.92900,    /* mod  9, chip 13 chn 25 & 26  zoom 4 */
    2.85942,  -2.85587,  89.92900,    /* mod  9, chip 14 chn 27 & 28  */
    2.85934,  -5.71174, 269.85600,    /* mod 10, chip 15 chn 29 & 30  zoom 5+90deg rot */
    2.85934,  -5.71174,  89.85600,    /* mod 10, chip 16 chn 31 & 32  */
    0.00000,   5.71893,  90.00000,    /* mod 11, chip 17 chn 33 & 34  zoom 3 */
    0.00000,   5.71893, 270.00000,    /* mod 11, chip 18 chn 35 & 36  */
    0.00000,   2.85945,  90.00000,    /* mod 12, chip 19 chn 37 & 38  zoom 2 */
    0.00000,   2.85945, 270.00000,    /* mod 12, chip 20 chn 39 & 40  */
    0.00000,   0.00000,  90.00000,    /* mod 13, chip 21 chn 41 & 42  zoom 1 */
    0.00000,   0.00000, 270.00000,    /* mod 13, chip 22 chn 43 & 44  */
    0.00000,  -2.85945, 270.00000,    /* mod 14, chip 23 chn 45 & 46  zoom 2 */
    0.00000,  -2.85945,  90.00000,    /* mod 14, chip 24 chn 47 & 48  */
    0.00000,  -5.71893, 270.00000,    /* mod 15, chip 25 chn 49 & 50  zoom 3 */
    0.00000,  -5.71893,  90.00000,    /* mod 15, chip 26 chn 51 & 52  */
    -2.85934,   5.71174,  89.85600,    /* mod 16, chip 27 chn 53 & 54  zoom 5+90 deg rot */
    -2.85934,   5.71174, 269.85600,    /* mod 16, chip 28 chn 55 & 56  */
    -2.85942,   2.85587,  89.92900,    /* mod 17, chip 29 chn 57 & 58  zoom 4 */
    -2.85942,   2.85587, 269.92900,    /* mod 17, chip 30 chn 59 & 60  */
    -2.85945,   0.00000,   0.00000,    /* mod 18, chip 31 chn 61 & 62  zoom 2 */
    -2.85945,   0.00000, 180.00000,    /* mod 18, chip 32 chn 63 & 64  */
    -2.85942,  -2.85586,   0.07160,    /* mod 19, chip 33 chn 65 & 66  zoom 4 */
    -2.85942,  -2.85586, 180.07160,    /* mod 19, chip 34 chn 67 & 68  */
    -2.85934,  -5.71174, 270.14400,    /* mod 20, chip 35 chn 69 & 70  zoom 5+90deg rot */
    -2.85934,  -5.71174,  90.14400,    /* mod 20, chip 36 chn 71 & 72  */
    -5.71881,   2.84513,  -0.14080,    /* mod 22, chip 37 chn 73 & 74  zoom 5 */
    -5.71881,   2.84513, 179.85920,    /* mod 22, chip 38 chn 75 & 76  */
    -5.71893,   0.00000,   0.00000,    /* mod 23, chip 39 chn 77 & 78  zoom 3 */
    -5.71893,   0.00000, 180.00000,    /* mod 23, chip 40 chn 79 & 80  */
    -5.71881,  -2.84513,   0.14080,    /* mod 24, chip 41 chn 81 & 82  zoom 5 */
    -5.71881,  -2.84513, 180.14080,    /* mod 24, chip 42 chn 83 & 84  */
    /* for FGS reference position is 6, an outside corner of the chip , 
       not the center of the module */
    /* FGS 1 and 4 are 90 deg rot of 2 and 3 with appropriate sign changes 
       also added 0.0064 to y in 2&3 and to x in 1&4 to shift by 12 pixels */	
    4.92625,    4.79173,    90.21,	    /* mod 1,  fgs 1, chn 85 added 3/4/05 season 0 */
    4.80933,   -4.90900,   179.79,		/* mod 5,  fgs 2, chn 86 added 3/4/05 season 3 */
    -4.80933,    4.90900,    -0.21,		/* mod 21, fgs 3, chn 87 added 3/4/05 season 1 */
    -4.92625,   -4.79173,   270.21		/* mod 25, fgs 4, chn 88 added 3/4/05 season 2 */	
  };

	
  int chn_n, chip_n, sid;
  double p_row, p_coln;
  float fchn_n,f_row,f_coln;
  double lngr, latm ;
  double alp1, del1, alp2, del2, delta; /* added 2/22/05 */	
  double tra,tdec, tra_sav, tdec_sav;/* target Ra and Dec in degrees and then radians */
  float hr, fmin,sec,dd,mm,ss;
  int ihr, imin, idd, imm;
  int c1,c2;
  int dread, chnn, rown, colnn;

  time_t systime;

#if 1
  *pChn_n = 0;
  *pRow = 0;
  *pColn = 0;
#endif

  systime=time(NULL);
#if 0
  printf("RA2fpix version 2.1\n%s  ",ctime(&systime));
#endif
  plate_scale_fgs=plate_scale*13.000/27.0 ; /* arc sec per pixel added 3/4/05 */
  foo=chip_offset[0][0]; /* added so array is used once and not get warning message */
  /*
    for(i=0;i<5;i++) { // test for print out of chn layout 
    for(j=0;j<10;j++){
    printf("%i %i %i\n",i-5,j-5,chn_num[i][j]); 
    }}
  */
#if 0	
  printf("\nSelect units for RA and DEC. Type d for decimal t for time of f for file: \n");
  printf("Or c to calculate RA and Dec of corners: ");
  scanf("%c",&ang);  /* x  z and m are also acceptable, 
                        /*z option uses the opt_xf angles w/o transforming to center of FOV 
                        z option converts CodeV zoom angles to pixels. THen adjusted offsets until corners aligned */
  /* x is also acceptable for doing corner calculations 
     and then transform back to see if get 0,0 for corners */
  /* m outputs the corner information in mongo format for Latham	*/
  /* f reads an input file in deg RA and Dec added 2/1/08 */
  /* s reads an input file in sexiamgiismal RA and Dec hh mm ss dd mm ss add 11/19/08 */
#else
  ang = 'd';
#endif

	
  check=0;
  while(check==0) {
    if(ang!='d' && ang!='t' && ang!='z' && ang!='c'&& ang!='x' && ang!='m' && ang!='f'&& ang!='s') {
      /* added z option 2/22/05 c option 7/20/05 m option 11/05  f on 2/1/08 */
      printf("\nPlease type  d  t  or  c  and a carriage return: ");
#if 1
      printf("ERROR: ang is %c\n",ang);
      return(0);
#else
      scanf("%c",&ang);
#endif
    } else check=1;
  }
	
  if(ang=='t'){
    printf("\nWarning: when entering negative declination between 0 and -1 deg\n");
    printf("You need to attach the - sign a nonzero value, either the min or sec\n\n");
  }

  if(ang=='f' || ang=='s') {
    convert_it=0;
    printf(" Enter file name:  " );
#if 1
    return(0);
#else
    scanf("%s", file2);
#endif
    printf(" File name %s\n",file2);
    if((fp2=fopen(file2,"r"))==NULL){ 
      printf("Input File doesn't exist!!!\n");
      exit(1);
    }
    if (ang=='s') {
      if((fp1=fopen("SIM-CRC.dat","w"))==NULL){
        printf("Output File doesn't exist!!!\n");
        exit(1);
      }
    }
    if(ang=='f'){
      for(i=0;i<8;i++){ /* read in 8 lines of header stuff */
        fgets(fileN,100,fp2); /* reads to end of line or 100 char */
        printf("-- %s",fileN);
      }
    }

  }


  if(ang=='d'|| ang=='t') {
#if 0
    convert_it=-1;
    while (convert_it<0 || convert_it>1){
      printf("Enter 0 if want to convert TO pixels or 1 to convert FROM pixels : ");
      scanf("%d",&convert_it);
    }
#else 
    convert_it = 0;
#endif
  }
#if 0
  season_int=-1;
  while (season_int<0 || season_int>3){
    printf("\nPick season (FOV rotation), summer =0, fall =1 winter =2, spring =3 : ");
    scanf("%d",&season_int);
  }
#else
  season_int = 0;
#endif
	
  /* Calculate the Direction Cosine Matrix to transform from RA and Dec to FPA coordinates */
		
  if(ang=='z') { /* added z option 2/22/05 */
    rot_1=90.0*season_int;
    rot_2=0.0; /* don't move to FOV center coor */
    rot_3=0.0; /* don't move to FOV center coor */
		
  } else {
    rot_1=first_roll+clocking_angle+season_int*90.; /* changed 9/10/04 */
    if(rot_1 > 360.) rot_1=rot_1-360.;
  }	
#if 0	
  printf("Using FOV center at RA=%f Dec=%f and Rot angle=%f\n",rot_3,-rot_2,rot_1);
  printf("Spacecraft roll angle is %f since FPA clocking adds %f deg to s/c roll\n",
         rot_1-clocking_angle, clocking_angle);
#endif	
  rot_1=rot_1+180.0; /* Need to account for 180 deg rotation of field due to imaging of mirror */
  srac=sin(rot_3*dtr); /* sin phi 3 rotation */
  crac=cos(rot_3*dtr); /* cos phi */
  sdec=sin(rot_2*dtr); /* sin theta 2 rotation 
                          Note 2 rotation is negative of dec in right hand sense */
  cdec=cos(rot_2*dtr); /* cos theta */
  srotc=sin(rot_1*dtr); /* sin psi 1 rotation */
  crotc=cos(rot_1*dtr); /* cos psi */
	
	
  /* DCM for spacecraft & focal plane 3-2-1 rotation, Wertz p764 */
  DCM11=cdec*crac;
  DCM12=cdec*srac;
  DCM13=-sdec;
  DCM21=-crotc*srac+srotc*sdec*crac;
  DCM22=crotc*crac+srotc*sdec*srac;
  DCM23=srotc*cdec;
  DCM31=srotc*srac+crotc*sdec*crac;
  DCM32=-srotc*crac+crotc*sdec*srac;
  DCM33=crotc*cdec;
	
  /* CALCULATE DCM for each chip relative to center of FOV includes FGS 43-46 */
  for (i=0;i<46;i++){ /* step through each chip */
    /*fprintf(stderr,"%d %f %f %f\n",i,chip_trans[i][0],chip_trans[i][1],chip_trans[i][2]); */
    srac=sin(chip_trans[i][0]*dtr); /* sin phi 3 rotation */
    crac=cos(chip_trans[i][0]*dtr); /* cos phi */
    sdec=sin(chip_trans[i][1]*dtr); /* sin theta 2 rotation */
    cdec=cos(chip_trans[i][1]*dtr); /* cos theta */
    srotc=sin((chip_trans[i][2]+chip_offset[i][0])*dtr); /* sin psi 1 rotation includes rotation offset */
    crotc=cos((chip_trans[i][2]+chip_offset[i][0])*dtr); /* cos psi */
		
    /* DCM for a 3-2-1 rotation, Wertz p762 */
    DCM11c[i]=cdec*crac;
    DCM12c[i]=cdec*srac;
    DCM13c[i]=-sdec;
    DCM21c[i]=-crotc*srac+srotc*sdec*crac;
    DCM22c[i]=crotc*crac+srotc*sdec*srac;
    DCM23c[i]=srotc*cdec;
    DCM31c[i]=srotc*srac+crotc*sdec*crac;
    DCM32c[i]=-srotc*crac+crotc*sdec*srac;
    DCM33c[i]=crotc*cdec;
  }
	
  /*********** finished setting up of matrices ***********/	
	
  if(ang=='c' || ang=='x' || ang=='m'){ /*********** find RA and dec of corners ********/
    for (chn_n=1;chn_n<85;chn_n++){ /* loop thru all Science chn */
      p_coln=0.;/* define outside corner of each chn = chip */
      for(c1=0;c1<2;c1++ ){ /* expanded to do row 1024 as well as row 0 */
        p_row=1024.0 *c1;
        /* also did one time run with p_row=512 and p_coln=550 to get centers */
        RAF(chn_n, p_row, p_coln, &tra, &tdec); /* get RA and dec FROM pixels */
				
        if(ang=='m') { /* generate mongo output */
          if(chn_n % 4 ==1){
            printf("relocate %f %f",tra,tdec);
            tra_sav=tra;
            tdec_sav=tdec;
          } else 
            printf("draw %f %f",tra,tdec);
          if(chn_n % 4 == 0) printf("\ndraw %f %f \n",tra_sav,tdec_sav);
				
        } else {/* convert to sexagismal */
          printf("chn=%2d row=%6.1f coln=%6.1f %f %f ",chn_n, p_row, p_coln, tra, tdec);
          ihr=(int) tra/15.;/* convert to decimal hours */
          fmin=60.*(tra/15.-ihr);
          imin= (int) fmin;
          sec=60.*(fmin-imin);
					
          idd= (int) tdec;
          mm=60.*(tdec-idd);
          imm= (int) mm;
          ss=60.*(mm-imm);
          printf(" RA %2dh %2dm %6.3fs  Dec %3dd %2d' %6.3f\" ",ihr,imin,sec,idd,imm,ss);
				
          if(ang=='x') { /* convert back to pixels */
            RA2(tra, tdec, &lngr, &latm, &chip_n, &chn_n, &p_row,&p_coln);
            printf("chn=%2d row=%8.3f coln=%8.3f", chn_n,p_row,p_coln);
          }
        }
        printf("\n");
      }
    }
    /* now do FGS channels */
    for (chn_n=85;chn_n<89;chn_n++){ 	
      for (c1=0;c1<2;c1++){
        for (c2=0;c2<2;c2++){
          p_row=c1*528;
          p_coln=c2*548;
          RAF(chn_n, p_row, p_coln, &tra, &tdec); /* get RA and dec FROM pixels */
          if(ang=='m') { /* generate mongo output */
            if(c1==0 && c2==0){
              printf("relocate %f %f\n",tra,tdec);
              tra_sav=tra;
              tdec_sav=tdec;
            } else 
              printf("draw %f %f\n",tra,tdec);
            if(c1==1 && c2==1) printf("draw %f %f \n\n",tra_sav,tdec_sav);
						
          } else {
            printf("chn=%2d row=%6.1f coln=%6.1f %f %f",chn_n, p_row, p_coln, tra, tdec);
            /* convert to sexagismal */
            tra=tra; 	
            ihr=(int) tra/15.;/* convert to decimal hours */
            fmin=60.*(tra/15.-ihr);
            imin= (int) fmin;
            sec=60.*(fmin-imin);
					
            idd= (int) tdec;
            mm=60.*(tdec-idd);
            imm= (int) mm;
            ss=60.*(mm-imm);
            printf(" RA %2dh %2dm %6.3fs  Dec %3dd %2d' %6.3f\" ",ihr,imin,sec,idd,imm,ss);
				
            if(ang=='x') { /* convert back to pixels */
              RA2(tra, tdec, &lngr, &latm, &chip_n, &chn_n, &p_row,&p_coln);
              printf("chn=%2d row=%8.3f coln=%8.3f", chn_n,p_row,p_coln);
            }
            printf("\n");
          }
        }
      }
    }	
    /* finished c and x option */
	
  } else if(ang=='z') { /************ z option convert zoom angles to pix */

    /* treat opt_xf[][0] as if RA and opt_xt[][1] as if dec */
		
    /* convert angles from CodeV from azim,azim to azim,dec */
    for(zi=0;zi<nzt;zi++){
      opt_xf[zi][1]=rtd*atan( tan(opt_xf[zi][1]*dtr) * cos(opt_xf[zi][0]*dtr) ); /* fixed 3/6/05 */
    }
		
    /* calculate separation between corners added 2/22/05 */
    printf("Calculated width and height of all sides of each module as check\n");
    for (zi=0;zi<7;zi++){ /* loop thru zoom 2-5a */
      for(ci=0;ci<4;ci++) { /*loop thru corners 1 to 4 */
        if(ci==0){ cc=5*zi+1; cx=cc+2; /* corners 8-9 */
        } else if(ci==1) { cc=5*zi+3; cx=cc+1; /* corners 9-6 */
        } else if(ci==2) { cc=5*zi+4; cx=cc-2; /* corners 6-5 */
        } else if(ci==3) { cc=5*zi+2; cx=cc-1; /* corners 5-8 */
        }
				
        alp1=opt_xf[cc][0];
        del1=opt_xf[cc][1];

        alp2=opt_xf[cx][0];
        del2=opt_xf[cx][1];
        sep( alp1, del1, alp2, del2, &delta, 0);
        printf("alp1 %f del1 %f alp2 %f del2 %f \n",alp1,del1,alp2,del2);
        if(zi<5) printf("For Z %d separation=%f¡\n", zi+1,delta);
        else if (zi==5) printf("For Z 5a separation=%f¡\n",delta);  
        /* fixed zoom 5a, since zoom 6 is now FGS 6/24/05 */
        else if (zi==6) printf("For Z 6 separation=%f¡\n",delta);  
        /* fixed zoom 5a, since zoom 6 is now FGS 6/24/05 */
      }
    } 
		
    /* check FGS sep center to four corners  added 7/18/05 */
    alp1=opt_xf[30][0]; /* FGS center */
    del1=opt_xf[30][1];
    for (cx=31;cx<35;cx++){ /* loop thru 4 corners */
      alp2=opt_xf[cx][0]; /* index 1=corner 6, 2 is 7, 3 is 9 and 4 is 10 */ 
      del2=opt_xf[cx][1];
      sep( alp1, del1, alp2, del2, &delta, 0);
      printf("FGS center to corner %d is %f¡\n",cx-30,delta);
    }
    printf("\n\n");
		
    /* now do transforms of zoom angles to pixels */
    for(nz=0;nz<nzt;nz++){
      tra=opt_xf[nz][0];
      tdec=opt_xf[nz][1];
      printf("RA=%f  Dec=%f\n",tra,tdec);
      RA2(tra, tdec, &lngr, &latm, &chip_n, &chn_n, &p_row,&p_coln);
      /**/ printf("chip long=%9.5f lat=%9.5f chip=%2d ",lngr/3600.,latm/3600.,chip_n+1); /**/
      printf("chn=%2d row=%8.3f coln=%8.3f\n", chn_n,p_row,p_coln);
    }
    /* finished z option */

  } else if(convert_it==0){ /************ Do conversion from RA & Dec to Pixels */
    tra=1.; /* so can enter while loop */
    while(tra>=-1000.) {
      if(ang=='d') {
#if 0
        printf("\nEnter RA and dec in degrees: ");
        scanf("%lf %lf",&tra, &tdec);
#else
        tra = inputRA;
        tdec = inputDec;
#endif

        check=0;

        while(check==0){
          if(tra<-180. || tra>360. || tdec<-90. || tdec>90.) {
            printf("\nAngles out of range. Reenter both: ");
#if 1
            printf("ERROR: (1) tra %f tdec %f\n",tra,tdec);
            return(0);
#else
            scanf("%lf %lf",&tra, &tdec);
#endif
          } else check=1;
        }
      } else if(ang=='t') {/* ang in time  */ 
        check=0;
        printf("\nEnter RA and Dec in hh mm ss and deg min sec: ");
#if 1
        return(0);
#else
        scanf("%f %f %f   %f %f %f", &hr, &fmin, &sec, &dd, &mm, &ss);
#endif
        /*		fprintf(stderr,"\n Have input ");
                        fprintf(stderr," %3.0f %3.0f %4.1f   %3.0f %3.0f %4.1f \n",hr, fmin, sec, dd, mm, ss); /**/

        while(check==0){
          tra=15.*(fabs(hr)+(fabs(fmin)+fabs(sec)/60.)/60.);
          if(hr<0. || fmin<0. || sec<0.) tra = -tra;
          tdec=fabs(dd)+(fabs(mm)+fabs(ss)/60.)/60.;
          if(dd<0. || mm<0. || ss<0.) tdec=-tdec;
			
          if(tra<-180. || tra>360. || tdec<-90. || tdec>90.) {
            printf("\nRA or Dec out of range. Reenter both: ");
#if 1
            printf("ERROR: (2) tra %f tdec %f\n",tra,tdec);
            return(0);
#else
            scanf("%f %f %f   %f %f %f", &hr, &fmin, &sec, &dd, &mm, &ss);
#endif
          } else  check=1;
        }			
      } else if (ang=='f'){ /* reading from file produced by corners added 2/1/08 */
				
	
        fgets(fileN,5,fp2);  printf("%s ",fileN); /* reads in 5 char */
        dread=fscanf(fp2, "%d %s %f  ", &chnn, fileN, &rown); 
        dread=fscanf(fp2, "%s %f  ", fileN, &colnn); 
        dread=fscanf(fp2, "%lf %lf  ", &tra, &tdec); 
        fgets(fileN,50,fp2);  /* reads to end of line  or 50 char */
        if(dread==EOF) exit(0);
        printf(" %d %f %f \n", chnn, tra, tdec);
        if (chnn==89) exit(0);
				
			
      } else {/* read from file sexigismal values */
        dread=fscanf(fp2,"%d %f %f %f   %f %f %f", &sid, &hr, &fmin, &sec, &dd, &mm, &ss);
        fprintf(stderr,"%d %3.0f %3.0f %4.1f   %3.0f %3.0f %4.1f \n", sid,hr, fmin, sec, dd, mm, ss); /**/
        if(dread==EOF) {
          fprintf(stderr,"Finished\n");
          exit (0);
        }
        check=0;
        while(check==0){
          tra=15.*(fabs(hr)+(fabs(fmin)+fabs(sec)/60.)/60.);
          if(hr<0. || fmin<0. || sec<0.) tra = -tra;
          tdec=fabs(dd)+(fabs(mm)+fabs(ss)/60.)/60.;
          if(dd<0. || mm<0. || ss<0.) tdec=-tdec;
			
          if(tra<-180. || tra>360. || tdec<-90. || tdec>90.) {
            fprintf(stderr,"%d %f %f %f   %f %f %f \n",sid, hr, fmin, sec, dd, mm, ss); 
            printf("\nRA or Dec out of range. star ID %d %f %f",sid,tra,tdec);
            exit(0);
          } else  check=1;
        }	
			
			
      }
			
      /*	fprintf(stderr," tra=%f tdec=%f\n", tra,tdec );/* */
			
      /* have input tra, tdec now do conversion to chn_n row and column */

      RA2(tra, tdec, &lngr, &latm, &chip_n, &chn_n, &p_row,&p_coln);

      /*printf("\nEnter target RA and Dec in deg (neg. RA ends loop): ");
        scanf("%lf %lf",&tra,&tdec);
        printf("%f %f\n", tra, tdec);*/

      if(ang=='s'){/* print to file */
        fprintf(fp1, "%d %f %f %2d %8.3f %8.3f", sid,tra/15., tdec, chn_n,p_row,p_coln);
        /* if(sid>20) exit(0); /**/
        if(chn_n==-99 || chn_n==0 ) {
          fprintf(fp1, "   *\n");
        } else {
          if(p_row<0. || p_coln<0. || p_row>1024. || p_coln>1100.) 
            fprintf(fp1,"   **\n");
          else if( chn_n>84 &&(p_row>528 || p_coln>548)) /*   added 3/4/05 */
            fprintf(fp1,"   ***\n");  
          /* fgs has fewer rows and columns than a sci chip  added 3/4/05 */
          else fprintf(fp1,"\n");
				
        }
      }else{
        if(chn_n==-99){
#if 0
          printf("Coordinate well out of the FOV ***\n");
#endif
#if 1
          return(0);
#endif
        } else if(chn_n==0) {
#if 0
          printf("Coordinate in Fine Guidance Sensor corner, but off FGS chip ***\n");
#else
          return(0);
#endif
        } else {
#if 0
          /**/ printf("chip long=%9.5f lat=%9.5f chip=%2d ",lngr/3600.,latm/3600.,chip_n+1); 	
          printf("chn=%2d row=%8.3f coln=%8.3f", chn_n,p_row,p_coln);
#endif				
          if(p_row<-10. || p_coln<-10. || p_row>1034. || p_coln>1100.) {
#if 0
            printf("   ***\n");
#else
            return(0);
#endif
          } else if( chn_n>84 &&(p_row>528 || p_coln>548)) {/*   added 3/4/05 */
#if 0
            printf("   ***\n");
#else
            return(0);
#endif
          /* fgs has fewer rows and columns than a sci chip  added 3/4/05 */
          } else {
#if 0
            printf("\n");
#else
#if 1
            *pChn_n = chn_n;
            *pRow = p_row;
            *pColn = p_coln;
#endif

            return(1);
#endif
          }
        }
      }
    }
    /* finsihed convert_it 0 case RA Dec to pix */
#if 1
    printf("RA2fpix ERROR: should not be here!\n");
    exit(-1);
#endif	
  } else { /*************** convert_it==1 do conversion from pixels to RA and Dec */
    printf("To quit, enter a negative chn number\n");
		
    printf("Enter chn # (1-84), row # (0-1024) column # (0-1100): \n");
    printf("For FGS chn # (85-88), row # (0-528) column # (0-548): \n");
		
    printf("You may use decimal values for row from 0.0 to <1024.0 and column from 0.0 to <1100.0 : ");
#if 1
    return(0);
#else
    scanf( "%f %f %f", &fchn_n,  &f_row, &f_coln);
#endif
    chn_n=(int)fchn_n;
    p_row=(double)f_row;
    p_coln=(double)f_coln;
		
    while(chn_n>=1){ /* enter loop */

      RAF(chn_n, p_row, p_coln, &tra, &tdec); /* get RA and dec FROM pixels */
			
      if(ang=='t'){ /* convert to sexagismal */
        tra=tra; /* convert to decimal hours */
        ihr=(int) tra/15.;
        fmin=60.*(tra/15.-ihr);
        imin= (int) fmin;
        sec=60.*(fmin-imin);
				
        idd= (int) tdec;
        mm=60.*(tdec-idd);
        imm= (int) mm;
        ss=60.*(mm-imm);
        printf(" RA %2dh %2dm %6.3fs  Dec %3dd %2d' %6.3f\" \n",ihr,imin,sec,idd,imm,ss);
				
      } else {
        printf("RA=%f Dec=%f  ",tra,tdec);
      }
      f_row=-1; /* to get into next while so can read next value
                   Doesn't check to see if have correct FGS row col */
      while(f_row<0. || f_row>1065. || f_coln<0. || f_coln>1100.) {
        printf("Enter chn # (1-84), row # (0-1024) column # (0-1100): \n");
#if 1
        return(0);
#else
        scanf( "%f %f %f", &fchn_n,  &f_row, &f_coln);
#endif
      }
      chn_n=(int)fchn_n;  /* if entered negative chn_n while loop will end */
      p_row=(double)f_row;
      p_coln=(double)f_coln;
    }
  }
  printf("Finished\n");
  return 0;
}

