// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 * hammer.c - Test of Hammer-Aitoff plotting routines
 * hammer transformation borrowed from Generic Mapping Tools V4.0
 * Copyright (c) 1991-2004 by P. Wessel and W. H. F. Smith
 *
 * 
 * gcc -ggdb -O2  -I/usr/include/plplot -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include plothammer.c pipelineutils.a -L /dasch/install/lib -lplplotd  -lm  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread    -ltable -lutil  -lwcs -lgd  -o plothammer
 *
 * 
 * NOTE: to use plothammer on Odyssey, switch the loaded plplot version:
 * module unload hpc/plplot-5.9.4
 * module load module load  hpc/plplot-5.9.9
 *  
 * plothammer -c -l 500  -t 'Scanned Plates' -i /home/scanner/scanner/linux/gmt/coverage2.txt -o /home/scanner/scanner/linux/gmt/coverage.png
 * plothammer -c -l 500 -g -t 'Scanned Plates' -i /home/scanner/scanner/linux/gmt/coverage2galactic.txt -o /home/scanner/scanner/linux/gmt/coveragegalactic.png
 *
 * Public release: Note: edit LABELBLOCK to add new releases
 * 
 * plothammer -c -l 0 -r    -p -t 'DR3 Release' -i /home/scanner/scanner/linux/gmt/coverage2.txt -o /home/scanner/scanner/linux/gmt/coveragepublic.png
 * plothammer -c -l 0 -r -g -p -t 'DR3 Release' -i /home/scanner/scanner/linux/gmt/coverage2galactic.txt -o /home/scanner/scanner/linux/gmt/coveragegalacticpublic.png
 * plothammer -c -l 0 -g -p -t 'M44 Release' -i /home/scanner/scanner/linux/gmt/coverage2galactic.txt -o /home/scanner/web/dasch/coveragegalacticm44.png
 *
 * Numbered release map  
 *
 * plothammer -c -l 0 -r  -p -m -a -t "DASCH Release Sequence " -i /home/scanner/scanner/linux/gmt/coverage2galactic.txt -o /home/scanner/scanner/linux/gmt/releasepublic.png
 * plothammer -c -l 0 -r -g -p -m -a -t "DASCH Release Sequence " -i /home/scanner/scanner/linux/gmt/coverage2galactic.txt -o /home/scanner/scanner/linux/gmt/releasegalacticpublic.png
 *
 * Apass plots
   plothammer -l 10 -c -s 'stars/deg**2' -u 10000 -t 'APASS DR6' -i /home/scanner/scanner/linux/gmt/apass_temp.txt -o /home/scanner/Pipeline/losdr6.png
   plothammer -l 10 -c -s 'stars/deg**2' -u 10000 -t 'APASS DR7' -i /home/scanner/scanner/linux/gmt/apass_temp4.txt -o /home/scanner/Pipeline/losdr7.png
 *
 * Limiting magnitude tables (On odyssey /n/dasch8/scanner/web/limiting)
 * 
 *  equatorial
 *
   plothammer -l 10 -t "Kepler Input Catalog Calibration " -q kepler -o /home/scanner/web/dasch/newlimiting/limitingkepler.png
   plothammer -l 10 -t "APASS Catalog Calibration " -q apass -o /home/scanner/web/dasch/newlimiting/limitingapass.png
   plothammer -l 10 -t "GSC2.3.2 Catalog Calibration " -q gsc2.3.2 -o /home/scanner/web/dasch/newlimiting/limitinggsc.png
 *
 *  galactic
 *
   plothammer -g -l 10 -t "Kepler Input Catalog Calibration " -q kepler -o /home/scanner/web/dasch/newlimiting/gallimitingkepler.png
   plothammer -g -l 10 -t "APASS Catalog Calibration " -q apass -o /home/scanner/web/dasch/newlimiting/gallimitingapass.png
   plothammer -g -l 10 -t "GSC2.3.2 Catalog Calibration " -q gsc2.3.2 -o /home/scanner/web/dasch/newlimiting/gallimitinggsc.png
 
 *  pngtopnm < gallimitingapass_ALL_15.png | pnmscalefixed -quiet -ysize=100  | cjpeg -rgb -quality 100 > gallimitingapass_ALL_15_thumb.jpg
 *
 *  pngtopnm < coveragepublic.png | cjpeg -rgb -quality 100 > coveragepublic.jpg
 *  pngtopam < coveragepublic.png | pnmquant 256 | pamtogif > coveragepublic.gif
 *  ppm, pgm, bmp, Targa, RLE
 *
 * pngtopnm limitingkepler_ALL_12.png | pnmtops > limitingkepler_ALL_12.ps
 * ~/Pipeline/pscombine.csh limiting*.ps
 * 
 *  grid only

    plothammer -d -t "Galactic Coordinate Grid" -o /home/scanner/junk/galacticgrid.png
    pngtopnm galacticgrid.png | pnmtops  > galacticgrid.ps
    ps2pdf galacticgrid.ps

 *   rematchFlag 

      plothammer  -f -t 'REMATCH' -s 'count' -i /home/scanner/Pipeline/los1.db -o /home/scanner/Pipeline/los1.png

 *
 * Apr 12, 2013 Edward J. Los - Initial version
 * Apr 16, 2013 Edward J. Los - Add upper, lower limits and title
 * Apr 17, 2013 Edward J. Los - Obtain and plot limiting magnitudes
 * May 10, 2013 Edward J. Los - Support APASS plots
 * Jan 21, 2013 Edward J. Los - Correct for Fedora 20 interface change
 * Mar 11, 2013 Edward J. Los - Introduce simple release map colored by release status.
 * Jun  4, 2014 Edward J. Los - Support DR3; Count plates in each release field
 * Dec 22, 2014 Edward J. Los - Support DR4
 * Feb 13, 2015 Edward J. Los - Reorder releases into 15 degree bands from the NGP to SGP
 * Jun  5, 2015 Edward J. Los - add -d qualifier to plot the galactic grid only
 * Jul 31, 2015 Edward J. los - add -e qualifier to plot empty bins
 * Nov 17, 2015 Edward J. Los - Support DR5
 * Nov 27, 2017 Edward J. Los - Support DR6
 * Sep 22, 2018 Edward J. Los - Add -f for Optimized TC (rematch) plots
 * Aug 30, 2019 Edward J. Los - Add FORCE_RELEASE_FIELD for all-sky plots
 * Nov 12, 2019 Edward J. Los - Stop at DR7 for the release*public.png plots
 * Nov 16, 2019 Edward J. Los - Restore the DR8-DR12 labels in release*public.png plots, but keep the yellow color
 */
#define _GNU_SOURCE 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "plplot.h"
#include "plplotP.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

#define NUMDEC 181
#define NUMRA  361
#if 0
#define NSHADES 5
#define NCOLORBAR 6
#else
#define NSHADES 100
#define NCOLORBAR 11
#endif
#define NPOINTDEF 100
#define MAX_BUFFER 512
#define STEP 1
#define MINDEC -90.0
#define MAXDEC +90.0
#define MINRA    0.0
#define MAXRA  360.0
#define DEGREES_TO_RAD  (3.141592654/180.0)
#define RAD_TO_DEGREES  (180.0/3.141592654)
#define MAGNITUDE_MAX 21
/* #define LIMIT_REMATCH 1 */ /* Clean up for Josh 2018-11-07 */
/* #define FORCE_RELEASE_FIELD  RELEASE_FIELD_DR6 */

#if 0
#define LOW_CUTOFF 130.0
#define HIGH_CUTOFF 230.0
#else
#define LOW_CUTOFF 145.0
#define HIGH_CUTOFF 210.0
#endif


extern GSCBIN gscBin01;
PGSCBIN pGscBin01 = &gscBin01;
extern char *releaseFieldText[RELEASE_FIELD_MAX+1];

#define M_SQRT2 1.41421356237309504880
#define GMT_CONV_LIMIT	1.0e-8	/* Fairly tight convergence limit or "close to zero" limit */
#if 1
#define EQ_RAD (1.0/M_SQRT2)
#else
#define EQ_RAD 6371008.7714
#endif
#define M_PI 3.14159265358979323846
#define D2RHAMMER (M_PI / 180.0)
#define R2DHAMMER (180.0 / M_PI)
#define d_asin(x) (fabs (x) >= 1.0 ? copysign (M_PI_2, (x)) : asin (x))


typedef struct _labelblock {
  int releaseField;
	double lon;
	double lat;
	char *labelText;
} LABELBLOCK,*PLABELBLOCK;

typedef struct _imageTotal {
  int group;
  int count;
} IMAGETOTAL,*PIMAGETOTAL;


/* Region labels */
LABELBLOCK labels[] = {	
#if 1
  /* Release order revision of Februrary, 2015 */
 	{RELEASE_FIELD_DR1   ,0.0+255.0, 82.5,"DR1"},
	{RELEASE_FIELD_DR2   ,0.0+255.0, 67.5,"DR2"},
	{RELEASE_FIELD_DR3   ,0.0+255.0, 52.5,"DR3"},
	{RELEASE_FIELD_DR4   ,0.0+255.0, 37.5,"DR4"},
	{RELEASE_FIELD_DR5   ,0.0+255.0, 22.5,"DR5"},
	{RELEASE_FIELD_DR6   ,0.0+255.0,  7.5,"DR6"},
	{RELEASE_FIELD_DR7   ,0.0+255.0, -7.5,"DR7"},
	{RELEASE_FIELD_DR8   ,0.0+255.0,-22.5,"DR8"},
	{RELEASE_FIELD_DR9   ,0.0+255.0,-37.5,"DR9"},
	{RELEASE_FIELD_DR10  ,0.0+255.0,-52.5,"DR10"},
	{RELEASE_FIELD_DR11  ,0.0+255.0,-67.5,"DR11"},
	{RELEASE_FIELD_DR12  ,0.0+255.0,-82.5,"DR12"},

#else
  /* Original release sequence */
 	{RELEASE_FIELD_DR1   ,0.0+255.0, 82.5,"DR1"},
	{RELEASE_FIELD_DR2   ,0.0+255.0, 67.5,"DR2"},
	{RELEASE_FIELD_DR3   ,0.0+255.0, 52.5,"DR3"},
	{RELEASE_FIELD_DR4   ,0.0+255.0, 37.5,"DR4"},
	{RELEASE_FIELD_DR5   ,0.0+255.0, 22.5,"DR5"},
	{RELEASE_FIELD_DR6   ,0.0+255.0,-82.5,"DR6"},
	{RELEASE_FIELD_DR7   ,0.0+255.0,-67.5,"DR7"},
	{RELEASE_FIELD_DR8   ,0.0+255.0,-52.5,"DR8"},
	{RELEASE_FIELD_DR9   ,0.0+255.0,-37.5,"DR9"},
	{RELEASE_FIELD_DR10  ,0.0+255.0,-22.5,"DR10"},
	{RELEASE_FIELD_DR11  ,0.0+255.0, -7.5,"DR11"},
	{RELEASE_FIELD_DR12  ,0.0+255.0,  7.5,"DR12"},
#endif

	{RELEASE_FIELD_M44   ,205.91,32.48,"M44"},
	{RELEASE_FIELD_3C273 ,289.96,64.36+4.0,"3C273"},
	{RELEASE_FIELD_BAADE ,1.03,-3.91,"Baade's Window"},
  {RELEASE_FIELD_KEPLER,76.34,13.45,"Kepler Field"},
	{RELEASE_FIELD_LMC   ,280.47,-32.89+8.0,"LMC"}
};
int labelSize = sizeof(labels)/sizeof(LABELBLOCK);

/* Data Release dividers */
#if 0
double latDivider[] = {-15.0};
#else
double latDivider[] = {75.0,60.0,45.0,30.0,15.0,0.0,-15.0,-30.0,-45.0,-60.0,-75.0};
#endif
int latDividerSize = sizeof(latDivider)/sizeof(double);
void    plot_galactic_grid(char *outfile,char* title);
void plot_rematch_grid(double *raTable,double *decTable,double *countTable,double *sequenceTable,int curCount,double *gscVector,int count,char *outfile,int galacticFlag,double upperLimit,double lowerLimit,char* title,char* ralabel,char* declabel,char* colormap_scale,int releaseField,int catalogNumber);


void GMT_hammer(double lon, double lat, double *x, double *y)
{
	/* Convert lon/lat to Hammer-Aitoff Equal-Area x/y */
	double slat, clat, slon, clon, D;
	
	if (fabs (fabs (lat) - 90.0) < GMT_CONV_LIMIT) {	/* Save time */
		*x = 0.0;
		*y = M_SQRT2 * copysign(EQ_RAD,lat) ;
		return;
	}
	
	while (lon < -180.0) lon += 360.0;
	while (lon > 180.0) lon -= 360.0;
	lat *= D2RHAMMER;
	lon *= (0.5 * D2RHAMMER);
	sincos (lat, &slat, &clat);
	sincos (lon, &slon, &clon);
	
	D = EQ_RAD * sqrt (2.0 / (1.0 + clat * clon));
	*x = 2.0 * D * clat * slon;
	*y = D * slat;
}

void GMT_ihammer (double *lon, double *lat, double x, double y)
{
	/* Convert Hammer_Aitoff Equal-Area x/y to lon/lat */
	double rho, c, angle;
	
	x *= 0.5;
	rho = hypot (x, y);
	
	if (fabs (rho) < GMT_CONV_LIMIT) {
		*lat = 0.0;
		*lon = 0.0;
	}
	else {
		c = 2.0 * d_asin (0.5 * rho * (1.0/EQ_RAD));
		*lat = d_asin (y * sin (c) / rho) * R2DHAMMER;
		if (fabs (c - M_PI_2) < GMT_CONV_LIMIT)
			angle = (fabs (x) < GMT_CONV_LIMIT) ? 0.0 : copysign (180.0, x);
		else
			angle = 2.0 * R2DHAMMER * atan (x * tan (c) / rho);
		*lon =  angle;
	}
}


void hammer_projection(double ra,double dec,double *pXval,double *pYval)
{
	double xval;
	double yval;
	double racosval;
	double rasinval;
	double deccosval;
	double decsinval;
	double denominator;
	double zval;
	double rainv;
	double decinv;
#if 1
	GMT_hammer (ra,dec,&xval,&yval);
	GMT_ihammer(&rainv,&decinv,xval,yval);
	
#else
	racosval = cos(DEGREES_TO_RAD*(ra)/2.0);
	rasinval = sin(DEGREES_TO_RAD*(ra)/2.0);
	deccosval = cos(DEGREES_TO_RAD*(dec));
	decsinval = sin(DEGREES_TO_RAD*(dec));
	denominator = sqrt(1+(racosval*deccosval));
	xval = 2.*sqrt(2.)*deccosval*rasinval/denominator;
	yval = 2.*sqrt(2.)*decsinval/denominator;
	zval = sqrt(1- sqr(xval/4.0) - sqr(yval/2.0));
	rainv = 2*atan((zval*xval)/(2.0*((2.0*sqr(zval)) -1))) * RAD_TO_DEGREES;
	decinv = asin(zval*yval) * RAD_TO_DEGREES;
#endif
	*pXval = xval;
	*pYval = yval;
#if 0
	printf("%f\t%f\t%f\t%f\t%f\t%f\n",ra,dec,xval,yval,rainv,decinv);
#endif
	return;
}

void grid_pltr2( PLFLT x, PLFLT y, PLFLT *tx, PLFLT *ty, PLPointer pltr_data )
{
	if (x > 180.0) {
		*tx = x - 360.0;
	} else {
		*tx = x;
	}
	*ty = y-90.0;
#if 0
	printf("%f\t%f\t%f\t%f\n",x,y,*tx,*ty);
#endif
	return;
}
void grid_pltr3( PLFLT x, PLFLT y, PLFLT *tx, PLFLT *ty, PLPointer pltr_data )
{
  *tx = x;
	*ty = y;
#if 0
	printf("%f\t%f\t%f\t%f\n",x,y,*tx,*ty);
#endif
	return;
}
void grid_pltr( PLFLT x, PLFLT y, PLFLT *tx, PLFLT *ty, PLPointer pltr_data )
{
	*tx = x;
	*ty = y-90.0;
#if 0
	printf("%f\t%f\t%f\t%f\n",x,y,*tx,*ty);
#endif
	return;
}


void temp_pltr( PLFLT x, PLFLT y, PLFLT *tx, PLFLT *ty, PLPointer pltr_data )
{
#if 0 /* Unity */
	*tx = x;
	*ty = y-90.0;
#if 0
	printf("temp_pltr x %f y %f tx %f  ty %f\n",x,y,*tx,*ty);
#endif
	return;
#endif /* Unity */
#if 1 /* From GMT 4.0 */
	{
		double xval;
		double yval;
		GMT_hammer (x-180,y-90.0,&xval,&yval);		
		*tx = (xval+2)*90.0;
		*ty = yval * 90.0;
	}
	return;
#endif

#if 0 /* From Wikipedia */
	double xval;
	double yval;
	double racosval;
	double rasinval;
	double deccosval;
	double decsinval;
	double denominator;
	racosval = cos(DEGREES_TO_RAD*(x-90.0)/2.0);
	rasinval = sin(DEGREES_TO_RAD*(x-90.0)/2.0);
	deccosval = cos(DEGREES_TO_RAD*(y-90.0));
	decsinval = sin(DEGREES_TO_RAD*(y-90.0));
#if 0
	printf("x %f y %f racosval %f rasinval %f deccosval %f decsinval %f\n",x,y,racosval,rasinval,deccosval,decsinval);
#endif
	denominator = sqrt(1+(racosval*deccosval));
	xval = 2.*sqrt(2.)*deccosval*rasinval/denominator;
	yval = 2.*sqrt(2.)*decsinval/denominator;
		
	*tx = (xval)*90.0;
	*ty = yval*22.5;
#if 0
	printf("temp_pltr x %f y %f xval %f yval %f tx %f  ty %f\n",x,y,xval,yval,*tx,*ty);
#endif
#endif /* Wikipedia */

}


void plot_grid(double *raTable,double *decTable,double *countTable,double *gscVector,int count,char *outfile,int galacticFlag,double upperLimit,double lowerLimit,char* title,char* colormap_scale,int releaseField,int catalogNumber)
{
	int index;
	double minCount = 1000000;
	double minCountx = 1000000;
	double maxCount = 0;
	double maxCountx = 0;
	int argc;
	double xmin;
	double xmax;
	double ymin;
	double ymax;

	double xminmm;
	double xmaxmm;
	double yminmm;
	double ymaxmm;
	double mscale;
	double dscale;

	char str[512];
	PLFLT	**scale;
	PLFLT clevel[NSHADES];
	PLFLT shedge[NCOLORBAR];
	int xindex;
	int yindex;

	double plateCount;
	double ra;
	double dec;
	PLFLT cpoint[3];
	PLFLT r[3];
	PLFLT g[3];
	PLFLT b[3];
	PLFLT x[NPOINTDEF+1];
	PLFLT y[NPOINTDEF+1];
	int idec;
	int ira;
	int gsc_bin_index;
	int decBin;
	int raBin;
	int just = 1.0;

	plAlloc2dGrid(&scale, NUMRA, NUMDEC);

#if 0
	printf("ra\tdec\tcount\n");
	printf("--\t---\t-----\n");
#endif
	for (xindex = 0; xindex < NUMRA; xindex++) {
		for (yindex = 0; yindex < NUMDEC; yindex++) {
			ra = 1.0*xindex;
			dec = -90.0 + (1.0*yindex);
			if (dec >= 90.0) {
				dec = 89.99999;
			}
			if (galacticFlag) {
				ra = ra - 180.0;
				if (ra < 0.0) {
					ra += 360.0;
				}
			}

			gsc_bin_index = GetGSCBin(pGscBin01,ra,dec,&decBin,&raBin,"plot_hammer1");
			plateCount = gscVector[gsc_bin_index];


#if 0
			printf("%f\t%f\t%f\n",ra,dec,plateCount);
#endif
			if ((plateCount > 0 ) && (minCount > plateCount)) {
				minCount = plateCount;
			}
	
			if (maxCount < plateCount) {
				maxCount = plateCount;
			}
			if (lowerLimit > -1)  {
				if (plateCount < lowerLimit) {
					plateCount = 0;
				}
			}
			if (upperLimit > -1) {
				if (plateCount >= upperLimit) {
					plateCount = upperLimit-1;
				}
			}
			scale[xindex][yindex] = plateCount;
		}
	}
 
	if (lowerLimit > -1) {
		minCountx = lowerLimit;
	} else {
		minCountx = minCount;
	}
	if (minCountx <= 0) {
		minCountx = 1;
	}
	if (upperLimit > -1) {
		maxCountx = upperLimit;
	} else {
		maxCountx = maxCount;
	}
	
	printf("minCount %f minCountx %f, maxCount %f, maxCountx %f\n",minCount,minCountx,maxCount,maxCountx);

	if ((minCountx == 1) && (maxCountx == 1)) {
		minCountx = 0;
	}

	dscale = (1.0*(maxCountx-minCountx))/(1.0*(NSHADES-1));
	if (dscale <= 0) {
		printf("ERROR: dscale (1) is %f for %s\n",dscale,title);
		return;
	}

	for (index = 0; index < NSHADES; index++) {
		clevel[index] = minCountx + (index * dscale);
#if 0
		printf("index %d clevel %f\n",index,clevel[index]);
#endif
	}
	dscale = (1.0*(maxCountx-minCountx))/(1.0*(NCOLORBAR-1));
	if (dscale <= 0) {
		printf("ERROR: dscale (2) is %f for %s\n",dscale,title);
		return;
	}

	for (index = 0; index < NCOLORBAR; index++) {
		shedge[index] = minCountx + (index * dscale);
#if 0
		printf("index %d shedge %f\n",index,shedge[index]);
#endif
	}


#if 1
	plsdev("png");
  plsetopt("-o",outfile);
#endif
#if 0
	plsdev("psc");
  plsetopt("-o","loshammer.psc");
#endif
#if 0
	plsdev("xwin");
#endif
#if 0
	plsori(1);  /* Set to portrait */
#endif
  plscolbg(255,255,255);	/* Force the background colour to white */  
  plscol0(1, 0,0,0);		/* Force the foreground colour to black */
  plscol0(15,255,0,0);		/* Move red to 15 */
	plinit();
  if (releaseField == RELEASE_FIELD_KEPLER) {
		if (galacticFlag) {
			xmin = 90.0;
			xmax = 65.0;
			ymin =   0.0;
			ymax =  25.0;
			
		} else {
			xmin = 305.0;
			xmax = 275.0;
			ymin =  35.0;
			ymax =  53.0; /* Originally 55.0 */


		}

	} else if (releaseField == RELEASE_FIELD_M44) {
		if (galacticFlag) {
			xmin = 215.0;
			xmax = 195.0;
			ymin =  40.0;
			ymax =  25.0;
			
		} else {
			xmin = 140.0;
			xmax = 120.0;
			ymin =  10.0;
			ymax =  30.0;


		}
	} else if (releaseField == RELEASE_FIELD_3C273) {
		if (galacticFlag) {
			/* Note xmin and xmax chosen so the plot looks reasonably circular */
			xmin = 290.0+22.5;
			xmax = 290.0-22.5;
			ymin =  55.0;
			ymax =  70.0;
			just = 0;
			
		} else {
			xmin = 195.0;
			xmax = 180.0;
			ymin =  -5.0;
			ymax =  10.0;


		}
	}  else if (releaseField == RELEASE_FIELD_BAADE) {
		if (galacticFlag) {
			xmin =  15.0;
			xmax = -15.0;
			ymin = -15.0;
			ymax =   5.0;
			
		} else {
			xmin = 280.0;
			xmax = 260.0;
			ymin = -40.0;
			ymax = -20.0;


		}
	}  else if (releaseField == RELEASE_FIELD_LMC) {
		if (galacticFlag) {
			xmin = 290.0;
			xmax = 270.0;
			ymin = -25.0;
			ymax = -40.0;
			
		} else {
			xmin = 100.0;
			xmax =  60.0;
			ymin = -80.0;
			ymax = -60.0;


		}
	}  else if (releaseField == RELEASE_FIELD_DR1) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin =  70.0;
			ymax =  90.0;
			just = 0;
			
		} else {
			xmin = 215;
			xmax = 175.0;
			ymin =   5.0;
			ymax =  50.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR2) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin =  55.0;
			ymax =  80.0;
			just = 0;
			
		} else {
			xmin = 230;
			xmax = 155;
			ymin = -5.0;
			ymax =  60.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR3) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin =  40.0;
			ymax =  65.0;
			just = 0;
			
		} else {
			xmin = 250;
			xmax = 135;
			ymin = -20.0;
			ymax =  75.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR4) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin =  25.0;
			ymax =  50.0;
			just = 0;
			
		} else {
			xmin = 280.0;
			xmax = 110.0;
			ymin = -40.0;
			ymax =  90.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR5) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin =  10.0;
			ymax =  35.0;
			just = 0;
			
		} else {
			xmin =  360.0;
			xmax =  0.0;
			ymin = -50.0;
			ymax =  90.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR6) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin =  -5.0;
			ymax =  20.0;
			just = 0;
			
		} else {
			xmin =  360.0;
			xmax =  0.0;
			ymin = -70.0;
			ymax =  85.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR7) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin = -20.0;
			ymax =   5.0;
			just = 0;
			
		} else {
			xmin =  360.0;
			xmax =  0.0;
			ymin = -90.0;
			ymax =  70.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR8) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin = -35.0;
			ymax = -10.0;
			just = 0;
			
		} else {
			xmin =  360.0;
			xmax =  0.0;
			ymin = -90.0;
			ymax =  90.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR9) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin = -50.0;
			ymax = -25.0;
			just = 0;
			
		} else {
			xmin =  360.0;
			xmax =  0.0;
			ymin = -90.0;
			ymax =  90.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR10) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin = -65.0;
			ymax = -40.0;
			just = 0;
			
		} else {
			xmin =  360.0;
			xmax =  0.0;
			ymin = -90.0;
			ymax =  90.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR11) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin = -80.0;
			ymax = -55.0;
			just = 0;
			
		} else {
			xmin =  360.0;
			xmax =  0.0;
			ymin = -90.0;
			ymax =  90.0;


		}
	} else if (releaseField == RELEASE_FIELD_DR12) {
	 
		if (galacticFlag) {
			xmin = 360.0;
			xmax =   0.0;
			ymin = -90.0;
			ymax = -70.0;
			just = 0;
			
		} else {
			xmin =  360.0;
			xmax =  0.0;
			ymin = -90.0;
			ymax = -70.0;


		}
	} else {
    if (catalogNumber == CATALOG_KEPLER) {
      if (galacticFlag) {
        xmin = 92.0;
        xmax = 62.0;
        ymin = -4.0;
        ymax =  25.0;
      } else {
        xmin = 310.0;
        xmax = 274.0;
        ymin =  28.0;
        ymax =  56.0;
      }

    } else {
      xmin = 360.0;
      xmax = 0.0;
      ymin = -90.0;
      ymax =  90.0;
    }
	}
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax,just, 2); 
	plwind(xmin,xmax,ymin,ymax);
	plgspa(&xminmm,&xmaxmm,&yminmm,&ymaxmm);
	printf("xminmm %f xmaxmm %f yminmm %f ymaxmm %f\n",xminmm,xmaxmm,yminmm,ymaxmm);
	/*	plwind(xmin,xmax,ymin,ymax); */
	plcol0(15);
#if 0
	plpoin(count,raTable,decTable,1);
#endif
#if 0 /* unsaturated table */
	cpoint[0] = 0.0; r[0] = 1.0; g[0] = 0.5; b[0] = 0.5;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 0.5; g[2] = 0.5; b[2] = 1.0;
#endif
#if 1 /* saturated table */
	cpoint[0] = 0.0; r[0] = 1.0; g[0] = 0.0; b[0] = 0.0;
  cpoint[1] = 0.5; r[1] = 0.0; g[1] = 1.0; b[1] = 0.0;
  cpoint[2] = 1.0; r[2] = 0.0; g[2] = 0.0; b[2] = 1.0;
#endif
#if 0 /* Inverted table */
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;

#endif
  plscmap1l(1, 3, cpoint, r, g, b, NULL);


	printf("xmin %f xmax %f ymin %f ymax %f\n",xmin,xmax,ymin,ymax);
	if ((galacticFlag != 0) && (releaseField == RELEASE_FIELD_BAADE)) {
		plshades((const PLFLT * const *)scale,NUMRA,NUMDEC,NULL,xmin,xmax,ymin,ymax,clevel,NSHADES,1,0,0,plfill,0,grid_pltr2,NULL);
	} else {
		plshades((const PLFLT * const *)scale,NUMRA,NUMDEC,NULL,xmin,xmax,ymin,ymax,clevel,NSHADES,1,0,0,plfill,0,grid_pltr,NULL);
	}

  plcol0(1);
	if (galacticFlag) {
		pllab("Galactic Longitude","Latitude",title);
	} else {
		pllab("J2000 RA","J2000 Dec",title);
	}
	{
		PLFLT colorbar_width;
		PLFLT colorbar_height;
		PLINT cont_color = 0;
		PLINT cont_width = 0;
		// Smaller text
		plschr( 0.0, 0.6);
		// Small ticks on the vertical axis
		plsmaj( 0.0, 0.5 );
		plsmin( 0.0, 0.5 );
#ifndef plwidth  /* Fedora 18 version */
		plcolorbar( &colorbar_width, &colorbar_height,
								PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL |PL_COLORBAR_LABEL_TOP, 0, 
								0.001,0.0,   /* x,y */
								0.02, 0.875,  /* x_length,y_length */
								0, 1, 1,        /* bg_color,bb_color, bb_style */
								0.0, 0.0,       /* low_cap_color, high_cap_color */
								cont_color, cont_width, 0.0, 0, "bcvtm",colormap_scale,
								NCOLORBAR, shedge );
#else /* plwidth Fedora 20 version */
		{
			PLINT      n_labels     = 1;
			PLINT      label_opts[] = {
        PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL | PL_COLORBAR_LABEL_TOP,
			};
			const char *labels[] = {
        colormap_scale,
			};
			const char *axis_opts[] = {
        "bcvtm",
			};
			PLINT      num_values[1];
      PLFLT      *values[1];
			PLFLT      axis_ticks[1] = {
        0.0,
			};
			PLINT      axis_subticks[1] = {
        0,
			};


			values [0] = shedge;
			num_values[0] = NCOLORBAR;

			plcolorbar( &colorbar_width, &colorbar_height,
									PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL |PL_COLORBAR_LABEL_TOP, 0, 
									0.001,0.0,   /* x,y */
									0.02, 0.875,  /* x_length,y_length */
									0, 1, 1,        /* bg_color,bb_color, bb_style */
									0.0, 0.0,       /* low_cap_color, high_cap_color */
									cont_color, cont_width, n_labels, label_opts, 
									labels,                            /* Old arg 19 New arg 18 label */
									1,                                 /*            New arg 19 n_axes */
									axis_opts,                         /* Old arg 18 New arg 20 axis_opts*/
									axis_ticks,                        /*            New arg 21 ticks */
									axis_subticks,                     /*            New arg 22 sub_ticks */
								  num_values,                        /* Old arg 20 New arg 23 n_values */
									(const PLFLT * const *)values );   /* Old arg 21 New arg 24 values*/
		}

#endif /* plwidth */
		// Reset text and tick sizes
		plschr( 0.0, 1.0 );
		plsmaj( 0.0, 1.0 );
		plsmin( 0.0, 1.0 );
		

	}



	plend();
	plFree2dGrid(scale, NUMRA, NUMDEC);

	return;


}
void plot_hammer(double *raTable,double *decTable,double *countTable,double *gscVector,int count,char *outfile,int galacticFlag,double upperLimit,double lowerLimit,char* title,char *colormap_scale,int releaseLabelFlag,int allReleaseLabelFlag,int releaseMapFlag,int missingBinFlag)
{
	int index;
	double minCount = 1000000;
	double minCountx = 1000000;
	double maxCount = 0;
	double maxCountx = 0;
	int argc;
	double xmin;
	double xmax;
	double ymin;
	double ymax;

	double xminmm;
	double xmaxmm;
	double yminmm;
	double ymaxmm;
	double mscale;
	double dscale;

	char str[512];
	PLFLT	**scale;
	PLFLT clevel[NSHADES];
	PLFLT shedge[NCOLORBAR];
	int xindex;
	int yindex;

	double plateCount;
	double ra;
	double dec;
	PLFLT cpoint[3];
	PLFLT r[3];
	PLFLT g[3];
	PLFLT b[3];
	PLFLT x[NPOINTDEF+1];
	PLFLT y[NPOINTDEF+1];
	int idec;
	int ira;
	int gsc_bin_index;
	int decBin;
	int raBin;
	int labelIndex;
	PLABELBLOCK pLabelBlock;
	double lon;
	double lat;

	plAlloc2dGrid(&scale, NUMRA, NUMDEC);

#if 0
	printf("ra\tdec\tcount\n");
	printf("--\t---\t-----\n");
#endif
	for (xindex = 0; xindex < NUMRA; xindex++) {
		for (yindex = 0; yindex < NUMDEC; yindex++) {
			ra = 1.0*xindex;
			dec = -90.0 + (1.0*yindex);
			if (dec >= 90.0) {
				dec = 89.99999;
			}
			gsc_bin_index = GetGSCBin(pGscBin01,ra,dec,&decBin,&raBin,"plot_hammer1");
			plateCount = gscVector[gsc_bin_index];


#if 0
			printf("%f\t%f\t%f\n",ra,dec,plateCount);
#endif
			if ((plateCount > 0 ) && (minCount > plateCount)) {
				minCount = plateCount;
			}
      if (releaseMapFlag) {
        maxCount = 3.0;
      } else {
        if (maxCount < plateCount) {
          maxCount = plateCount;
        }
      }
			if (lowerLimit > -1)  {
				if (plateCount < lowerLimit) {
					plateCount = 0;
				}
			}
			if (upperLimit > -1) {
				if (plateCount >= upperLimit) {
					plateCount = upperLimit-1;
				}
			}
			scale[xindex][yindex] = plateCount;
		}
	}
	if (lowerLimit > -1) {
		minCountx = lowerLimit;
	} else {
		minCountx = minCount;
	}
	if (minCountx <= 0) {
		minCountx = 1;
	}
	if (upperLimit > -1) {
		maxCountx = upperLimit;
	} else {
		maxCountx = maxCount;
	}
	
	printf("minCount %f minCountx %f, maxCount %f, maxCountx %f\n",minCount,minCountx,maxCount,maxCountx);

	if ((minCountx == 1) && (maxCountx == 1)) {
		minCountx = 0;
	}

	dscale = (1.0*(maxCountx-minCountx))/(1.0*(NSHADES-1));
	if (dscale <= 0) {
		printf("ERROR: dscale (3) is %f for %s\n",dscale,title);
		return;
	}

	if (dscale <= 0) {
		printf("ERROR: dscale is %f for %s\n",dscale,title);
		return;
	}

	for (index = 0; index < NSHADES; index++) {
		clevel[index] = minCountx + (index * dscale);
#if 0
		printf("index %d clevel %f\n",index,clevel[index]);
#endif
	}
	dscale = (1.0*(maxCountx-minCountx))/(1.0*(NCOLORBAR-1));
	if (dscale <= 0) {
		printf("ERROR: dscale (4) is %f for %s\n",dscale,title);
		return;
	}
	for (index = 0; index < NCOLORBAR; index++) {
		shedge[index] = minCountx + (index * dscale);
#if 0
		printf("index %d shedge %f\n",index,shedge[index]);
#endif
	}


#if 1
	plsdev("png");
  plsetopt("-o",outfile);
#endif
#if 0
	plsdev("psc");
  plsetopt("-o","loshammer.psc");
#endif
#if 0
	plsdev("xwin");
#endif
#if 0
	plsori(1);  /* Set to portrait */
#endif
  plscolbg(255,255,255);	/* Force the background colour to white */  
  plscol0(1, 0,0,0);		/* Force the foreground colour to black */
  plscol0(15,255,0,0);		/* Move red to 15 */
	plinit();


  xmin = 360.0;
  xmax = 0.0;
  ymin = -90.0;
  ymax =  90.0;
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -2); 
	plwind(xmin,xmax,ymin,ymax);
	plgspa(&xminmm,&xmaxmm,&yminmm,&ymaxmm);
	printf("xminmm %f xmaxmm %f yminmm %f ymaxmm %f\n",xminmm,xmaxmm,yminmm,ymaxmm);
	/*	plwind(xmin,xmax,ymin,ymax); */
	plcol0(15);
#if 0
	plpoin(count,raTable,decTable,1);
#endif
#if 0 /* unsaturated table */
	cpoint[0] = 0.0; r[0] = 1.0; g[0] = 0.5; b[0] = 0.5;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 0.5; g[2] = 0.5; b[2] = 1.0;
#endif
#if 1 /* saturated table */
	cpoint[0] = 0.0; r[0] = 1.0; g[0] = 0.0; b[0] = 0.0;
  cpoint[1] = 0.5; r[1] = 0.0; g[1] = 1.0; b[1] = 0.0;
  cpoint[2] = 1.0; r[2] = 0.0; g[2] = 0.0; b[2] = 1.0;
#endif
#if 0 /* Inverted table */
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;

#endif
  plscmap1l(1, 3, cpoint, r, g, b, NULL);


	printf("xmin %f xmax %f ymin %f ymax %f\n",xmin,xmax,ymin,ymax);
  plshades((const PLFLT * const *)scale,NUMRA,NUMDEC,NULL,xmin,xmax,ymin,ymax,clevel,NSHADES,1,0,0,plfill,0,temp_pltr,NULL);
  plcol0(1);
	if (galacticFlag) {
		pllab("Galactic Longitude","Latitude",title);
	} else {
		pllab("J2000 RA","J2000 Dec",title);
	}

	if (releaseLabelFlag == 1) {
		pllsty(2);
		for (labelIndex = 0; labelIndex < latDividerSize; labelIndex++) {
			for (ira = 0; ira <= NPOINTDEF; ira ++) {
				lat = latDivider[labelIndex];
				lon = (360. * ira)/(1.0*NPOINTDEF);
				if (galacticFlag == 0) {
					wcscon(WCS_GALACTIC,WCS_J2000,2000.0,2000.0,&lon,&lat,2000.0);
					lon += 180.0;
					while (lon > 360.0) {
						lon -= 360;
					}
				}
				GMT_hammer(lon-180.0,lat,&x[ira],&y[ira]);	
				x[ira] = (x[ira]+2.0)*90.0;
				y[ira] = -y[ira]*90.0;
#if 0
				printf("labelIndex %d ira %d, lon %f lat %f\n",labelIndex,ira,x[ira],y[ira]);
#endif
			}
			if (galacticFlag == 0) {
				/* We can not cross the 360 degree line! */
				int baseira = 0;
				ira = 1;
				while (ira <=  NPOINTDEF ) {
					if (((x[ira] < LOW_CUTOFF) && (x[ira-1] >= HIGH_CUTOFF)) || 
							((x[ira] >=  HIGH_CUTOFF) && (x[ira-1] <  LOW_CUTOFF))) {
						plline(ira - baseira,&x[baseira],&y[baseira]);
						baseira = ira+1;
#if 0
						printf("labelIndex %d baseira %d ira %d\n",labelIndex,baseira,ira);
#endif
            plline(ira-baseira,&x[baseira],&y[baseira]);
					}
					ira++;
					
				}
#if 0
				printf("labelIndex %d baseira %d ira %d\n",labelIndex,baseira,ira);
#endif
				plline(ira-baseira-1,&x[baseira],&y[baseira]);


			} else {
				plline(NPOINTDEF+1,x,y);
			}
		}
		pllsty(1);

	}

	/* Now plot the axes */

	for (idec = -60; idec <= +60; idec += 30) {
		for (ira = 0; ira <= NPOINTDEF; ira++) {
			ra = (360. * ira)/(1.0*NPOINTDEF);
			dec = 1.0*idec;
			GMT_hammer(ra-180.0,dec,&x[ira],&y[ira]);
			x[ira] = (x[ira]+2.0)*90.0;
			y[ira] = y[ira]*90.0;
		}
		plline(NPOINTDEF+1,x,y);
	}
	for (ira = 0; ira <= 360; ira += 30) {
		for (idec = 0; idec <= NPOINTDEF; idec++) {
			ra = 1.0*ira;
			dec = ((180. * idec)/(1.0*NPOINTDEF)) - 90.0;
			GMT_hammer(ra-180.0,dec,&x[idec],&y[idec]);
			x[idec] = (x[idec]+2.0)*90.0;
			y[idec] = y[idec]*90.0;
		}
		plline(NPOINTDEF+1,x,y);
	}

	plschr( 0.0, 0.6);

	if (releaseLabelFlag == 1) {
#if 0
		for (index = 0; index < 360; index+= 10) {
			sprintf(str,"%d",index);
			GMT_hammer(1.0*index,82.5,&x[0],&y[0]);
			plptex((x[0]+2.0)*90.0,y[0]*90.0,-1.0,0.0,1.0,str);

		}
#endif
		for (labelIndex = 0; labelIndex < labelSize; labelIndex++) {
			pLabelBlock = &labels[labelIndex];
      if (allReleaseLabelFlag ==0) {
        if ((pLabelBlock->releaseField >= RELEASE_LEVEL) &&  
            (pLabelBlock->releaseField < RELEASE_FIELD_M44)) {
          continue;
        }
      } else {
#if 0   /* Removed    Nov 16, 2018 */
        /* Changed of Nov 12, 2019 */
        if ((pLabelBlock->releaseField >= RELEASE_FIELD_DR8) &&  
            (pLabelBlock->releaseField <= RELEASE_FIELD_DR12)) {
          continue;
        }
#endif
      }
			lon = pLabelBlock->lon;
			lat = pLabelBlock->lat;
			if (galacticFlag == 0) {
				wcscon(WCS_GALACTIC,WCS_J2000,2000.0,2000.0,&lon,&lat,2000.0);
				lon += 180.0;
				if (lon > 360.0) {
					lon -= 360;
				}
			}
			GMT_hammer(lon,lat,&x[0],&y[0]);
			plptex((x[0]+2.0)*90.0,y[0]*90.0,-1.0,0.0,0.5,pLabelBlock->labelText);
#if 0
      printf("line %d label %s\n",__LINE__,pLabelBlock->labelText);
#endif
		}
	}

	for (ira = 0; ira <= 360; ira += 90) {
		double raLabel;
		ra = 1.0*ira;
		if (galacticFlag) {
			raLabel = ra - 180.0;
			if (raLabel < 0.0) {
				raLabel += 360.0;
			}
			sprintf(str,"%.0f",raLabel);
		 
		} else {
			sprintf(str,"%d#uh",ira/15);
		}

		GMT_hammer(ra-180.0,2.0,&x[0],&y[0]);
		x[0] = (x[0]+2.0)*90.0;
		y[0] = y[0]*90.0;
		if (ira == 0) {
			plptex(x[0],y[0],-1.0,0.0,1.0,str);
		} else {
			plptex(x[0],y[0],-1.0,0.0,0.0,str);
		}
	}


	for (idec = -60; idec <= 60; idec += 30) {
		if (idec == 0) {
			continue;
		}
		sprintf(str,"%d#uo",idec);
		GMT_hammer(-180.0,1.0*idec,&x[0],&y[0]);
		x[0] = (x[0]+2.0)*90.0;
		y[0] = y[0]*90.0;
		plptex(x[0]-5.0,y[0],-1.0,0.0,0.0,str);
		GMT_hammer(+180.0,1.0*idec,&x[0],&y[0]);
		x[0] = (x[0]+2.0)*90.0;
		y[0] = y[0]*90.0;
		plptex(x[0]+10.0,y[0],-1.0,0.0,1.0,str);
	}


	{
		PLFLT colorbar_width;
		PLFLT colorbar_height;
		PLINT cont_color = 0;
		PLINT cont_width = 0;
		// Smaller text
		plschr( 0.0, 0.6);
		// Small ticks on the vertical axis
		plsmaj( 0.0, 0.5 );
		plsmin( 0.0, 0.5 );
    if ((releaseMapFlag == 0) && (missingBinFlag == 0)) {
#ifndef plwidth
      plcolorbar( &colorbar_width, &colorbar_height,
                  PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL |PL_COLORBAR_LABEL_TOP, 0, 
                  0.001,0.0,   /* x,y */
                  0.02, 0.875,  /* x_length,y_length */
                  0, 1, 1,        /* bg_color,bb_color, bb_style */
                  0.0, 0.0,       /* low_cap_color, high_cap_color */
                  cont_color, cont_width, 0.0, 0, "bcvtm",colormap_scale,
                  NCOLORBAR, shedge );
#else /* plwidth Fedora 20 version */
      {
        PLINT      n_labels     = 1;
        PLINT      label_opts[] = {
          PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL | PL_COLORBAR_LABEL_TOP,
        };
        const char *labels[] = {
          colormap_scale,
        };
        const char *axis_opts[] = {
          "bcvtm",
        };
        PLINT      num_values[1];
        PLFLT      *values[1];
        PLFLT      axis_ticks[1] = {
          0.0,
        };
        PLINT      axis_subticks[1] = {
          0,
        };

        values [0] = shedge;
        num_values[0] = NCOLORBAR;
        plcolorbar( &colorbar_width, &colorbar_height,
                    PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL |PL_COLORBAR_LABEL_TOP, 0, 
                    0.001,0.0,   /* x,y */
                    0.02, 0.875,  /* x_length,y_length */
                    0, 1, 1,        /* bg_color,bb_color, bb_style */
                    0.0, 0.0,       /* low_cap_color, high_cap_color */
                    cont_color, cont_width, n_labels, label_opts, 
                    labels,                           /* Old arg 19 New arg 18 label */
                    1,                                /*            New arg 19 n_axes */
                    axis_opts,                        /* Old arg 18 New arg 20 axis_opts*/
                    axis_ticks,                       /*            New arg 21 ticks */
                    axis_subticks,                    /*            New arg 22 sub_ticks */
                    num_values,                       /* Old arg 20 New arg 23 n_values */
                    (const PLFLT * const *)values );  /* Old arg 21 New arg 24 values*/
      }

#endif /* plwidth */
    }
		// Reset text and tick sizes
		plschr( 0.0, 1.0 );
		plsmaj( 0.0, 1.0 );
		plsmin( 0.0, 1.0 );
		

	}



	plend();
	plFree2dGrid(scale, NUMRA, NUMDEC);

	return;


}


int main(int argc,char *argv[])
{
	double *raTable = NULL;
	double *decTable = NULL;
	double *countTable = NULL;
  double *sequenceTable = NULL;
	double *gscVector = NULL;
	double ra;
	double dec;
	double count;
	int xindex;
	int yindex;
	int tindex = 0;
	double xval;
	double yval;
  char *argstr;
	char cmdchar;
  int verbose = 0;
  int galacticGrid = 0;
	int galacticFlag = 0;
	int publicFlag = 0;
	int	releaseLabelFlag = 0;
  int allReleaseLabelFlag = 0;
  FILE *input_handle = NULL;
  char input_name[MAX_BUFFER];
  char outfile[MAX_BUFFER]; 
	char magnitudeFile[MAX_BUFFER];
	char tempBuffer[MAX_BUFFER];
	char *pDot;
	int errorFlag = 0;
	int complementFlag = 0;
  int missingBinFlag = 0;
  int rematchFlag = 0;
  char inLine[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
	int lineCount = 0;
	int curCount = 0;
	int authorizedCount = 0;
	int nvals;
	int index;
	int gsc_bin_index;
	int decBin;
	int raBin;
	int unitializedCount = 0;
	int authorized;
	char title[MAX_BUFFER];
  char ralabel[MAX_BUFFER];
  char declabel[MAX_BUFFER];
	char colormap_scale[MAX_BUFFER];
	double upperLimit = -1.0;
	double lowerLimit = -1.0;
	int catalogNumber = -1;
	char source[MAX_BUFFER];
	char catalogString[MAX_BUFFER];
	char qualifier[MAX_BUFFER];
  char* catalogdir;
  GALAXYCOMMON galaxycommon;
  PGALAXYCOMMON pGalaxyCommon = &galaxycommon;
	int releaseField;
	double *pLimitingVector;
	double *limitingVector[RELEASE_FIELD_MAX+1][MAGNITUDE_MAX];
	int plateCountLimiting[RELEASE_FIELD_MAX+1][MAGNITUDE_MAX];
	int plateCountMaximum[RELEASE_FIELD_MAX+1];
	int plateCount;
	int releaseFieldIndex;
	int magnitudeIndex;
	PPLATELIMITINGREC pPlateRec;
	int limitingMagLocal;
	int releaseFieldCount[RELEASE_FIELD_MAX+1];
	int magnitudeCount[MAGNITUDE_MAX];
	int brightMagnitudeLimit;
	double lon;
	double lat;
  int releaseMapFlag = 0;
  int *seriesReleaseTable[MAX_SERIES+1];
  int *seriesReleaseArray;
  int plateNumber;
  int seriesId;
  int releaseFieldTotals[RELEASE_FIELD_MAX+1];
#ifdef FORCE_RELEASE_FIELD
  printf("ERROR: FORCE_RELEASE_FIELD defined as %d\n",FORCE_RELEASE_FIELD);
#endif

  memset(releaseFieldTotals,0,sizeof(releaseFieldTotals));
  memset(seriesReleaseTable,0,sizeof(seriesReleaseTable));
	memset(limitingVector,0,sizeof(limitingVector));
	memset(plateCountLimiting,0,sizeof(plateCountLimiting));
	memset(plateCountMaximum,0,sizeof(plateCountMaximum));
	memset(pGalaxyCommon,0,sizeof(GALAXYCOMMON));
	memset(releaseFieldCount,0,sizeof(releaseFieldCount));
	memset(magnitudeCount,0,sizeof(magnitudeCount));

	title[0] = 0;
	colormap_scale[0] = 0;
	input_name[0] = 0;
	outfile[0] = 0;
	qualifier[0] = 0;
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
				case 'i': /* input file */
				case 'I':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for input file -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(input_name,*++argv,MAX_BUFFER-2);
						if (strlen(input_name) >= MAX_BUFFER-3) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
							exit(-1);
						}
					}
					break;

				case 'o': /* output file name */
				case 'O':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(outfile,*++argv,MAX_BUFFER-2);
						if (strlen(outfile) >= MAX_BUFFER-3) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
							exit(-1);
						}
					}
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
          if (catalogNumber == CATALOG_GAIA) {
            printf("ERROR in line %d of plothammer - CATALOG_GAIA is not supported\n",__LINE__);
            exit(-1);
          }
#if 0
          if (catalogNumber == CATALOG_ATLAS) {
            printf("ERROR in line %d of plothammer - CATALOG_ATLAS is not supported\n",__LINE__);
            exit(-1);
          }
#endif
          break;

				case 't': /* title */
				case 'T':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for title -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(title,*++argv,MAX_BUFFER-2);
						if (strlen(title) >= MAX_BUFFER-3) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
							exit(-1);
						}
					}
					break;

				case 's': /* colormap_scale */
				case 'S':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for title -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(colormap_scale,*++argv,MAX_BUFFER-2);
						if (strlen(title) >= MAX_BUFFER-3) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
							exit(-1);
						}
					}
					break;

        case 'u': /* Upper plotted limit */
        case 'U':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&upperLimit);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the upperLimit %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'l': /* Lower plotted limit */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&lowerLimit);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the lowerLimit %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

				

				case 'v': /* verbose */
				case 'V':
					verbose = 1;
					break;

				case 'd': /* galactic grid */
				case 'D':
					galacticGrid = 1;
					break;

				case 'm': /* release map plot */
				case 'M':
					releaseMapFlag = 1;
					break;

				case 'g': /* galacticFlag */
				case 'G':
					galacticFlag = 1;
					break;

				case 'p': /* publicFlag */
				case 'P':
					publicFlag = 1;
					break;

				case 'r': /* label release regions */
				case 'R':
					releaseLabelFlag = 1;
					break;

				case 'a': /* label release regions */
				case 'A':
					allReleaseLabelFlag = 1;
					break;

				case 'c': /* RA is 360 - column value */
				case 'C':
					complementFlag = 1;
					break;

				case 'e': /* show bins without data */
				case 'E':
					missingBinFlag = 1;
          upperLimit = 100.0;
          lowerLimit = 0.5;
					break;

        case 'f':
        case 'F':
          rematchFlag = 1;
#ifdef LIMIT_REMATCH
          printf("ERROR: REMATCH_FLAG is set to simply the plot\n");
#endif /* LIMIT_REMATCH */
          break;


				default:
					printf("ERROR: * illegal command -%c-",cmdchar);
					errorFlag = 1;

				}
        
      }

    }
  }

	if ((input_name[0] == 0) && (catalogNumber < 0) && (galacticGrid == 0)) {
		printf("ERROR: no input file, limiting magnitude catalog, or galactic grid specified\n");
		errorFlag = 1;
	}
	if (outfile[0] == 0) {
		printf("ERROR: no output file specified\n");
		errorFlag = 1;
	}
	if (title[0] == 0) {
		printf("WARNING: no title specified\n");
	}
	if (colormap_scale[0] == 0) {
		strcpy(colormap_scale,"plates");
	}

	if (errorFlag) {
		printf("Usage: plothammer [-q <catalog>|-i<input file>] options\n");
    printf("  options: -v verbose\n");
    printf("           -i <input file>\n");
		printf("           -q <calibration catalog> to plot limiting magnitudes\n");
    printf("           -o <output file>\n");
    printf("           -g galactic coordinates\n");
    printf("           -p public data only\n");
		printf("           -r label the release regions\n");
		printf("           -a label all release regions\n");
		printf("           -t title\n");
    printf("           -m release map\n");
		printf("           -s colormap scale\n");
		printf("           -u upper limit of colorbar\n");
		printf("           -l lower limit of colorbar\n");
		printf("           -c input RA column has RA = 360-(column value)\n");		
    printf("           -d print the galactic grid only\n");
    printf("           -f show optimized TC match diagnostic plot\n");
    printf("           -e show empty bins only\n");
    return(-1);
	}		
	printf("plothammer of %s %s complementFlag %d\n",__DATE__,__TIME__,complementFlag);

	if (input_name[0] != 0) {
		input_handle = fopen(input_name,"rt");
		if (input_handle == NULL) {
			printf("ERROR: could not open input file %s\n",input_name);
			return(-1);
		}
		while(1) {
	
			/* Read in the input file and see how many lines it has */
			inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
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
      if (rematchFlag) {
        if (strchr(inBuffer,'\t') == NULL) {
          continue;
        }
#if 1
        printf("line %d lineCount %d %s\n",__LINE__,lineCount,inBuffer);
#endif
      }

			lineCount++;
		}

		fclose(input_handle);
		input_handle = NULL;
		if (lineCount < pGscBin01->total_gsc_bins) {
			lineCount = pGscBin01->total_gsc_bins;
		}
		raTable = (double *)calloc(lineCount,sizeof(double));
		decTable = (double *)calloc(lineCount,sizeof(double));
		countTable = (double *)calloc(lineCount,sizeof(double));
		sequenceTable = (double *)calloc(lineCount,sizeof(double));
		gscVector = (double *)calloc(lineCount,sizeof(double));

		/* Now allocate our arrays */
		if ((raTable == NULL) ||
				(decTable == NULL) ||
				(gscVector == NULL) ||
        (sequenceTable == NULL) ||
				(countTable == NULL)) {
			printf("ERROR: failed to allocate tables of size %d\n",lineCount);
			exit(-1);
		}
		for (index = 0; index < pGscBin01->total_gsc_bins; index++) {
			gscVector[index] = -1.0;
		}
		/* Now read in the data */

		input_handle = fopen(input_name,"rt");
		if (input_handle == NULL) {
			printf("ERROR: could not open input file %s\n",input_name);
			return(-1);
		}
		while(1) {
	
			/* Read in the input file and see how many lines it has */
			inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
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
      if (rematchFlag) {
        if (strchr(inBuffer,'\t') == NULL) {
          continue;
        }
        nvals = sscanf(inBuffer,"%lf\t%lf\t%lf\t%lf",&raTable[curCount],&decTable[curCount],&countTable[curCount],&sequenceTable[curCount]);
        if (nvals == 4) {
          curCount++;
        }  
      } else {

        nvals = sscanf(inBuffer,"%lf %lf %lf",&raTable[curCount],&decTable[curCount],&countTable[curCount]);
        if (nvals == 3) {
          if (complementFlag) {
            raTable[curCount] = 360.0 - raTable[curCount];
          }
          if (missingBinFlag) {
            if (countTable[curCount] > 0) {
              countTable[curCount] = 0;
            } else {
              countTable[curCount] = 1;
            }
          }

          gsc_bin_index = GetGSCBin(pGscBin01,raTable[curCount],decTable[curCount],&decBin,&raBin,"plothammer");
          if (publicFlag != 0) {
            if (galacticFlag) {
              /* ra, dec is really lon, lat -- convert to ra and dec. Variable names are reversed! */
              lon = raTable[curCount];
              lat = decTable[curCount];
              lon -= -180;
              if (lon < 0.0) {
                lon += 360.0;
              }
              wcscon(WCS_GALACTIC,WCS_J2000,2000.0,2000.0,&lon,&lat,2000.0);
              authorized = CheckAuthorization("none",lon,lat,&releaseField);


            } else {
              authorized = CheckAuthorization("none",raTable[curCount],decTable[curCount],&releaseField);
            } 
          } else {

            authorized = 1;

          }
          if (releaseMapFlag == 0) {
            if (authorized != 0) {
              authorizedCount++;
              gscVector[gsc_bin_index] = countTable[curCount];
            } else {
              gscVector[gsc_bin_index] = 0;
            }
          } else {
            switch (releaseField) {
            case RELEASE_FIELD_DR1:
            case RELEASE_FIELD_DR2:
            case RELEASE_FIELD_DR3:
            case RELEASE_FIELD_DR4:
            case RELEASE_FIELD_DR5:
            case RELEASE_FIELD_DR6:
            case RELEASE_FIELD_M44:
            case RELEASE_FIELD_3C273:
            case RELEASE_FIELD_BAADE:
            case RELEASE_FIELD_KEPLER:
            case RELEASE_FIELD_LMC:
              gscVector[gsc_bin_index] = 2.0;
              break;
            case RELEASE_FIELD_DR8:
            case RELEASE_FIELD_DR7:  /* Change of Nov 12, 2019 */
            case RELEASE_FIELD_DR9:  /* Change of Nov 12, 2019 */
            case RELEASE_FIELD_DR10: /* Change of Nov 12, 2019 */
            case RELEASE_FIELD_DR11: /* Change of Nov 12, 2019 */
            case RELEASE_FIELD_DR12: /* Change of Nov 12, 2019 */
            default: /* Change of Nov 12, 2019 */
              gscVector[gsc_bin_index] = 1.5;
              break;
#if 0 /* Change of Nov 12, 2019 */
            default:
              gscVector[gsc_bin_index] = 0.0;
              break;
#endif
            }
          }
        }
        curCount++;
      }
			
		}
		printf("Bins read %d, authorized bins %d for release level %d\n",curCount,authorizedCount,RELEASE_LEVEL);

		fclose(input_handle);
		input_handle = NULL;

		if ((rematchFlag == 0) && (curCount != pGscBin01->total_gsc_bins)) {
			printf("WARNING: curCount %d is not equal to total_gsc_bins %d\n",curCount,pGscBin01->total_gsc_bins);
		}
		for (index = 0; index < pGscBin01->total_gsc_bins; index++) {
			if (gscVector[index] < 0) {
				if (index == 0) {
					if (gscVector[index+1] > 0.0) {
						gscVector[index] = gscVector[index+1];
					}
				} else if (index == (pGscBin01->total_gsc_bins-1)) {
					if (gscVector[index-1] > 0.0) {
						gscVector[index] = gscVector[index-1];
					}
				} else if ((gscVector[index-1] > 0) && (gscVector[index+1] > 0)) {
					gscVector[index] = (gscVector[index-1]+gscVector[index+1])/2.0;
				} else if (gscVector[index-1] > 0) {
					gscVector[index] = gscVector[index-1];
				} else if (gscVector[index+1] > 0) {
					gscVector[index] = gscVector[index+1];
				} else {
					unitializedCount++;
				}
			}
		}
		if ((rematchFlag == 0) && (unitializedCount > 0)) {
			printf("ERROR: %d unitialized bins are present\n",unitializedCount);
		}
    if (rematchFlag == 0) {
      plot_hammer(raTable,decTable,countTable,gscVector,curCount,outfile,galacticFlag,upperLimit,lowerLimit,title,colormap_scale,releaseLabelFlag,allReleaseLabelFlag,releaseMapFlag,missingBinFlag);
    } else {
      plot_rematch_grid(raTable,decTable,countTable,sequenceTable,curCount,gscVector,lineCount,outfile,galacticFlag,upperLimit,lowerLimit,tempBuffer,ralabel,declabel,colormap_scale,releaseFieldIndex,catalogNumber);
    }
	} else if (galacticGrid != 0) {
    plot_galactic_grid(outfile,title);
    exit(-1);
  } else  {
		/* Limiting magnitude study */
		/* Open the galaxy/limiting magnitude database and load in the file */
		OpenGalaxyFiles(pGalaxyCommon,source,NULL);
		lineCount = pGscBin01->total_gsc_bins;
		raTable = (double *)calloc(lineCount,sizeof(double));
		decTable = (double *)calloc(lineCount,sizeof(double));
		countTable = (double *)calloc(lineCount,sizeof(double));
		gscVector = (double *)calloc(lineCount,sizeof(double));
		if ((raTable == NULL) ||
				(decTable == NULL) ||
				(gscVector == NULL) ||
				(countTable == NULL)) {
			printf("ERROR: failed to allocate tables of size %d\n",lineCount);
			exit(-1);
		}
		for (releaseFieldIndex = 0; releaseFieldIndex <= RELEASE_FIELD_MAX; releaseFieldIndex++) {
			plateCountMaximum[releaseFieldIndex] = 0;
			for (magnitudeIndex = 0; magnitudeIndex < MAGNITUDE_MAX; magnitudeIndex++) {
				plateCountLimiting[releaseFieldIndex][magnitudeIndex] = 0;
				limitingVector[releaseFieldIndex][magnitudeIndex] = (double *)calloc(lineCount,sizeof(double));
				if (limitingVector[releaseFieldIndex][magnitudeIndex] == NULL) {
					printf("ERROR: failed to allocate limiting vector[%d][%d]\n",releaseFieldIndex,magnitudeIndex);
					exit(-1);
				}
				
			}
		}
		
		for (index = 0; index < pGscBin01->total_gsc_bins; index++) {
			gscVector[index] = -1.0;
		}
		for (curCount = 0; curCount < lineCount; curCount++) {
			if (((curCount+1) % 4000) == 0) {
				printf("At curCount %5d of %d\n",curCount+1,lineCount);
			}
			/* Now consider each gsc bin in turn */
			GetBinCenter(pGscBin01,curCount,&ra,&dec,"plothammer");
			if (galacticFlag) {
				/* Here the gscVector is actually in galactic coordinates but the origin is shifted */
				lon = ra;
				lat = dec;
				lon -= -180;
				if (lon < 0.0) {
					lon += 360.0;
				}
				ra = lon;
				wcscon(WCS_GALACTIC,WCS_J2000,2000.0,2000.0,&ra,&dec,2000.0);
				raTable[curCount] = lon;
				decTable[curCount] = lat;
			} else {
				raTable[curCount] = ra;
				decTable[curCount] = dec;
			}
			CheckAuthorization(NULL,ra,dec,&releaseFieldIndex);
#ifdef FORCE_RELEASE_FIELD
      if (releaseFieldIndex < FORCE_RELEASE_FIELD) {
        releaseFieldIndex = FORCE_RELEASE_FIELD;
      }
#endif /* FORCE_RELEASE_FIELD */
			/* Now read the limiting magnitudes for this bin */
			LoadGalaxyTable(pGalaxyCommon,ra,dec);
			for (plateCount = 0; plateCount < pGalaxyCommon->plateCount; plateCount++) {
				pPlateRec	= &pGalaxyCommon->plateLimitingBuffer[plateCount];
        /* For each plate, we are looking for the minimum release */
        if ((pPlateRec->seriesId > 0) && (pPlateRec->seriesId <= MAX_SERIES)) {
          seriesReleaseArray = seriesReleaseTable[pPlateRec->seriesId];
          if (seriesReleaseArray == NULL) {
            seriesReleaseArray = (int *)calloc(MAX_PLATE_NUMBER,sizeof(int*));
            if (seriesReleaseArray == NULL) {
              printf("ERROR allocating seriesReleaseArray\n");
              exit(-1);
            }
            seriesReleaseTable[pPlateRec->seriesId] = seriesReleaseArray;
            for (plateNumber = 0; plateNumber < MAX_PLATE_NUMBER; plateNumber++) {
              seriesReleaseArray[plateNumber] = RELEASE_FIELD_MAX;
            }
          }
          if ((pPlateRec->plateNumber >= 0) && (pPlateRec->plateNumber < MAX_PLATE_NUMBER)) {

            if (seriesReleaseArray[pPlateRec->plateNumber] > releaseFieldIndex) {
              seriesReleaseArray[pPlateRec->plateNumber] = releaseFieldIndex;
            }
          }
          
        }

				limitingMagLocal = pPlateRec->limiting_mag_local + 0.0001;
				if ((limitingMagLocal >= 0) && (limitingMagLocal < MAGNITUDE_MAX)) {

					for (magnitudeIndex = limitingMagLocal; magnitudeIndex >= 0; magnitudeIndex--) {
						plateCountLimiting[releaseFieldIndex][magnitudeIndex]++;

								
						releaseFieldCount[releaseFieldIndex]++;
						magnitudeCount[magnitudeIndex]++;
						pLimitingVector = limitingVector[releaseFieldIndex][magnitudeIndex];
						pLimitingVector[curCount]++;
						if (plateCountMaximum[releaseFieldIndex] < pLimitingVector[curCount]) {
							plateCountMaximum[releaseFieldIndex] = pLimitingVector[curCount];
						}


					}
				}
			}
			releaseFieldIndex = RELEASE_FIELD_MAX;
			for (plateCount = 0; plateCount < pGalaxyCommon->plateCount; plateCount++) {
				pPlateRec	= &pGalaxyCommon->plateLimitingBuffer[plateCount];
				limitingMagLocal = pPlateRec->limiting_mag_local + 0.0001;
				if ((limitingMagLocal >= 0) && (limitingMagLocal < MAGNITUDE_MAX)) {

					for (magnitudeIndex = limitingMagLocal; magnitudeIndex >= 0; magnitudeIndex--) {
						plateCountLimiting[releaseFieldIndex][magnitudeIndex]++;

								
						releaseFieldCount[releaseFieldIndex]++;
						magnitudeCount[magnitudeIndex]++;
						pLimitingVector = limitingVector[releaseFieldIndex][magnitudeIndex];
						pLimitingVector[curCount]++;
						if (plateCountMaximum[releaseFieldIndex] < pLimitingVector[curCount]) {
							plateCountMaximum[releaseFieldIndex] = pLimitingVector[curCount];
						}


					}
				}
			}
		}
		for (brightMagnitudeLimit = 0; brightMagnitudeLimit < (MAGNITUDE_MAX-1); brightMagnitudeLimit++) {
			if (magnitudeCount[brightMagnitudeLimit] != magnitudeCount[brightMagnitudeLimit+1]) {
				break;
			}
		}

		printf("                           ");
		for (magnitudeIndex = brightMagnitudeLimit; magnitudeIndex < MAGNITUDE_MAX; magnitudeIndex++) {
			printf(" %8d",magnitudeIndex);
		}
		printf("\n");
		for (releaseFieldIndex = 0; releaseFieldIndex <= RELEASE_FIELD_MAX; releaseFieldIndex++) {
			if (releaseFieldCount[releaseFieldIndex] == 0) {
				continue;
			}
			printf("Field %8s Max %8d",releaseFieldText[releaseFieldIndex],plateCountMaximum[releaseFieldIndex]);
			for (magnitudeIndex = brightMagnitudeLimit; magnitudeIndex < MAGNITUDE_MAX; magnitudeIndex++) {
				printf(" %8d",plateCountLimiting[releaseFieldIndex][magnitudeIndex]);
			}
			printf("\n");
				
		}
		/* Now generate the plots */
		for (releaseFieldIndex = 0; releaseFieldIndex <= RELEASE_FIELD_MAX; releaseFieldIndex++) {
			if (releaseFieldCount[releaseFieldIndex] == 0) {
				continue;
			}
			for (magnitudeIndex = brightMagnitudeLimit; magnitudeIndex < MAGNITUDE_MAX; magnitudeIndex++) {
				if (plateCountLimiting[releaseFieldIndex][magnitudeIndex] == 0) {
					continue;
				}
				for (curCount = 0; curCount < lineCount; curCount++) {
					pLimitingVector = limitingVector[releaseFieldIndex][magnitudeIndex];
					gscVector[curCount] = pLimitingVector[curCount];
				}
				strcpy(tempBuffer,outfile);
				pDot = strstr(tempBuffer,".");
				if (pDot != NULL) {
					*pDot = 0;
				}
				sprintf(magnitudeFile,"%s_%s_%02d.png",tempBuffer,releaseFieldText[releaseFieldIndex],magnitudeIndex);
				sprintf(tempBuffer,"%s Field: %s Limiting Magnitude %2d",title,releaseFieldText[releaseFieldIndex],magnitudeIndex);
				printf("Plotting %s with title %s\n",magnitudeFile,tempBuffer);
#if 0
				if ((releaseFieldIndex == RELEASE_FIELD_MAX) && (magnitudeIndex == 10)) {
#if 1
					plot_hammer(raTable,decTable,countTable,gscVector,lineCount,magnitudeFile,galacticFlag,1.0*plateCountMaximum[releaseFieldIndex],lowerLimit,tempBuffer,releaseLabelFlag,allReleaseLabelFlag,releaseMapFlag,missingBinFlag);
#else
					plot_grid(raTable,decTable,countTable,gscVector,lineCount,magnitudeFile,galacticFlag,1.0*plateCountMaximum[releaseFieldIndex],lowerLimit,tempBuffer,colormap_scale,releaseFieldIndex,catalogNumber);
#endif					
				}
#else
        if ((releaseFieldIndex == RELEASE_FIELD_MAX) && (catalogNumber == CATALOG_KEPLER)) {
          plot_grid(raTable,decTable,countTable,gscVector,lineCount,magnitudeFile,galacticFlag,1.0*plateCountMaximum[releaseFieldIndex],lowerLimit,tempBuffer,colormap_scale,releaseFieldIndex,catalogNumber);
        } else {
          if (
#ifdef FORCE_RELEASE_FIELD
              (releaseFieldIndex == FORCE_RELEASE_FIELD) ||
#endif /* FORCE_RELEASE_FIELD */
              (releaseFieldIndex == RELEASE_FIELD_OTHER) || 
                (releaseFieldIndex == RELEASE_FIELD_MAX)) {
            plot_hammer(raTable,decTable,countTable,gscVector,lineCount,magnitudeFile,galacticFlag,1.0*plateCountMaximum[releaseFieldIndex],lowerLimit,tempBuffer,colormap_scale,releaseLabelFlag,allReleaseLabelFlag,releaseMapFlag,missingBinFlag);
          } else {
            plot_grid(raTable,decTable,countTable,gscVector,lineCount,magnitudeFile,galacticFlag,1.0*plateCountMaximum[releaseFieldIndex],lowerLimit,tempBuffer,colormap_scale,releaseFieldIndex,catalogNumber);
          }
        }
#endif
			}
				
		}
    /* Now print out the plate table */
    for (seriesId = 0; seriesId <= MAX_SERIES; seriesId++) {
      seriesReleaseArray = seriesReleaseTable[seriesId];
      if (seriesReleaseArray == NULL) {
        continue;
      }
      for (plateNumber = 0; plateNumber < MAX_PLATE_NUMBER; plateNumber++) {
        releaseFieldIndex = seriesReleaseArray[plateNumber];
        releaseFieldTotals[releaseFieldIndex]++;
      }
    }
    for (releaseFieldIndex = 0; releaseFieldIndex < RELEASE_FIELD_MAX; releaseFieldIndex++) {
      printf("%6d plates for field (uniquely assigned) %s\n",releaseFieldTotals[releaseFieldIndex],releaseFieldText[releaseFieldIndex]);

    }

	} /* End of limiting magnitude study */

	if (raTable != NULL) {
		free(raTable);
	}
	if (decTable != NULL) {
		free(decTable);
	}
	if (countTable != NULL) {
		free(countTable);
	}
  if (sequenceTable != NULL) {
    free(sequenceTable);
  }
	if (gscVector != NULL) {
		free(gscVector);
	}



	



	for (releaseFieldIndex = 0; releaseFieldIndex <= RELEASE_FIELD_MAX; releaseFieldIndex++) {
		for (magnitudeIndex = 0; magnitudeIndex < MAGNITUDE_MAX; magnitudeIndex++) {
			if (limitingVector[releaseFieldIndex][magnitudeIndex] != NULL) {
				free(limitingVector[releaseFieldIndex][magnitudeIndex]);
			}
				
		}
	}


	if (input_handle != NULL) {
		fclose(input_handle);
	}



	return(0);
}

void plot_galactic_grid(char *outfile,char* title)
{

  int galacticFlag = 1;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double xminmm;
	double xmaxmm;
	double yminmm;
	double ymaxmm;
	int idec;
	int ira;
	double ra;
	double dec;
	PLFLT x[NPOINTDEF+1];
	PLFLT y[NPOINTDEF+1];
	char str[512];
#if 0
	int index;
	double minCount = 1000000;
	double minCountx = 1000000;
	double maxCount = 0;
	double maxCountx = 0;
	int argc;

	double mscale;
	double dscale;
  

	int xindex;
	int yindex;

	double plateCount;
	PLFLT r[3];
	PLFLT g[3];
	PLFLT b[3];
	int gsc_bin_index;
	int decBin;
	int raBin;
	int labelIndex;
	PLABELBLOCK pLabelBlock;
	double lon;
	double lat;
#endif


	plsdev("png");
  plsetopt("-o",outfile);
  plsetopt("-geometry","1440x1080");
  plscolbg(255,255,255);	/* Force the background colour to white */  
  plscol0(1, 0,0,0);		/* Force the foreground colour to black */
  plscol0(15,255,0,0);		/* Move red to 15 */
	plinit();


  xmin = 360.0;
  xmax = 0.0;
  ymin = -90.0;
  ymax =  90.0;
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -2); 
	plwind(xmin,xmax,ymin,ymax);
	plgspa(&xminmm,&xmaxmm,&yminmm,&ymaxmm);
	printf("xminmm %f xmaxmm %f yminmm %f ymaxmm %f\n",xminmm,xmaxmm,yminmm,ymaxmm);
	/*	plwind(xmin,xmax,ymin,ymax); */
	plcol0(15);


	printf("xmin %f xmax %f ymin %f ymax %f\n",xmin,xmax,ymin,ymax);
  plcol0(1);
	if (galacticFlag) {
		pllab("Galactic Longitude","Latitude",title);
	} else {
		pllab("J2000 RA","J2000 Dec",title);
	}


	/* Now plot the axes */

	for (idec = -80; idec <= +80; idec += 10) {
		for (ira = 0; ira <= NPOINTDEF; ira++) {
			ra = (360. * ira)/(1.0*NPOINTDEF);
			dec = 1.0*idec;
			GMT_hammer(ra-180.0,dec,&x[ira],&y[ira]);
			x[ira] = (x[ira]+2.0)*90.0;
			y[ira] = y[ira]*90.0;
		}
		plline(NPOINTDEF+1,x,y);
	}
	for (ira = 0; ira <= 360; ira += 15) {
		for (idec = 0; idec <= NPOINTDEF; idec++) {
			ra = 1.0*ira;
			dec = ((180. * idec)/(1.0*NPOINTDEF)) - 90.0;
			GMT_hammer(ra-180.0,dec,&x[idec],&y[idec]);
			x[idec] = (x[idec]+2.0)*90.0;
			y[idec] = y[idec]*90.0;
		}
		plline(NPOINTDEF+1,x,y);
	}

	plschr( 0.0, 0.6);


	for (ira = 0; ira <= 360; ira += 30) {
		double raLabel;
		ra = 1.0*ira;
		if (galacticFlag) {
			raLabel = ra - 180.0;
			if (raLabel < 0.0) {
				raLabel += 360.0;
			}
			sprintf(str,"%.0f",raLabel);
		 
		} else {
			sprintf(str,"%d#uh",ira/15);
		}

		GMT_hammer(ra-180.0,2.0,&x[0],&y[0]);
		x[0] = (x[0]+2.0)*90.0;
		y[0] = y[0]*90.0;
		if (ira == 0) {
			plptex(x[0],y[0],-1.0,0.0,1.0,str);
		} else {
			plptex(x[0],y[0],-1.0,0.0,0.0,str);
		}
	}


	for (idec = -60; idec <= 60; idec += 30) {
		if (idec == 0) {
			continue;
		}
		sprintf(str,"%d#uo",idec);
		GMT_hammer(-180.0,1.0*idec,&x[0],&y[0]);
		x[0] = (x[0]+2.0)*90.0;
		y[0] = y[0]*90.0;
		plptex(x[0]-5.0,y[0],-1.0,0.0,0.0,str);
		GMT_hammer(+180.0,1.0*idec,&x[0],&y[0]);
		x[0] = (x[0]+2.0)*90.0;
		y[0] = y[0]*90.0;
		plptex(x[0]+10.0,y[0],-1.0,0.0,1.0,str);
	}




	plend();

	return;


}
/* Borrowed from x12c.c for use in filling in the boxes */
void
plfbox( PLFLT x0, PLFLT y0 )
{
    PLFLT x[4], y[4];

    x[0] = x0 - 0.5;
    y[0] = y0 - 0.5;
    x[1] = x0 + 0.5;
    y[1] = y0 - 0.5;
    x[2] = x0 + 0.5;
    y[2] = y0 + 0.5;
    x[3] = x0 - 0.5;
    y[3] = y0 + 0.5;
    plfill( 4, x, y );
    plcol0( 1 );
    pllsty( 1 );
    plline( 4, x, y );
}
int IntCompare(const void *first, const void *second) {
  int numberFirst = *((int*)first);
  int numberSecond = *((int*)second);

  
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}


/* This plot takes output generated from photometryutils.c  ProcessNoneImagesX() when enableRematch = 1 and  REMATCH_DEBUG is defined
 * 
 *  The input is tab-separated debugging statements which form a starbase table.
 *  column 4 (sequence) defines the meanings of the other columns.
 *  column 1 (Q0) is the REMATCH_RESOLUTION for sequence == 0; otherwise it is the rematchXIndex
 *  column 2 (Q1) is the REMATCH_BINS for sequence == 0; otherwise it is the rematchYIndex
 *  column 3 (count) for sequence 0 is the gsc_bin_index defined by REMATCH_DEBUG
 *                   for sequence 1 is a star or stars of interest identifier and the image count is the group number
 *                   for sequence 2 is image count in the bin
 *                   for sequence 3 is the group index.
 *           
 */
#define SEQUENCE_PARAMETERS 0
#define SEQUENCE_INTEREST   1
#define SEQUENCE_IMAGECOUNT 2
#define SEQUENCE_GROUPCOUNT 3


void plot_rematch_grid(double *raTable,double *decTable,double *countTable,double *sequenceTable,int totalCount,double *gscVector,int count,char *outfile,int galacticFlag,double upperLimit,double lowerLimit,char* title,char* ralabel,char* declabel, char* colormap_scale,int releaseField,int catalogNumber)
{
	int index;
	double minCount = 1000000;
	double minCountx = 1000000;
	double maxCount = 0;
	double maxCountx = 0;
	int argc;
	double xmin;
	double xmax;
	double ymin;
	double ymax;

	double xminmm;
	double xmaxmm;
	double yminmm;
	double ymaxmm;
	double mscale;
	double dscale;

	char str[512];
	PLFLT	**scale;
	PLFLT clevel[NSHADES];
  PLFLT* contlevel = NULL;
  int* rematchGroupTable = NULL;
  int maxlevel;
  int indexlevel;
  int nlevel = 0;
	PLFLT shedge[NCOLORBAR];
	int xindex;
	int yindex;

	double plateCount;
	double ra;
	double dec;
	PLFLT cpoint[3];
	PLFLT r[3];
	PLFLT g[3];
	PLFLT b[3];
	PLFLT x[NPOINTDEF+1];
	PLFLT y[NPOINTDEF+1];
	int idec;
	int ira;
	int decBin;
	int raBin;
	int just = 1.0;
  PLINT nx = 0;
  PLINT ny = 0;
  PLINT ix = 0;
  PLINT iy = 0;
  int curCount;
  int intCount;
  double rematch_resolution = 0;
  int rematch_bins = 0;
  int gsc_bin_index = 0;
  int sequenceImageCount = 0;
  int sequenceGroupCount = 0;
  int sequenceInterestCount = 0;
  int groupIndex;
  int numGroups = 0;
  PIMAGETOTAL imageTotalTable;
  PIMAGETOTAL pImageTotal;
  int imageTotalIndex;
  int imageTotalCount;
  int maxImageIndex = -1;
  int maxImageCount = 0;
  

  static PLFLT pos[]   = { 0.0, 0.25, 0.5, 0.75, 1.0 };
#if 0 /* Original table */
  static PLFLT red[]   = { 0.0, 0.25, 0.5, 1.0, 1.0 };
  static PLFLT green[] = { 1.0, 0.5, 0.5, 0.5, 1.0 };
  static PLFLT blue[]  = { 1.0, 1.0, 0.5, 0.25, 0.0 };
#endif
#if 1 /* Unsaturated 0.0->0.5  0.25 ->0.625  0.5 -> .75 */
  static PLFLT red[]   = { 0.5, 00.625, 0.75, 1.0, 1.0 };
  static PLFLT green[] = { 1.0, 0.75, 0.75, 0.75, 1.0 };
  static PLFLT blue[]  = { 1.0, 1.0, 0.75, 00.625, 0.5 };
#endif

  /* Now parse the data */
  for (curCount = 0; curCount < totalCount; curCount++) {

    if (sequenceTable[curCount] == SEQUENCE_PARAMETERS) {
      rematch_resolution = raTable[curCount];
      rematch_bins = decTable[curCount];
      gsc_bin_index = countTable[curCount];
      nx = rematch_bins;
      ny = rematch_bins;
      printf("line %d rematch_resolution %.2f, rematch_bins %d, gsc_bin_index %d\n",__LINE__,rematch_resolution,rematch_bins,gsc_bin_index);
#ifdef LIMIT_REMATCH

#else /* LIMIT_REMATCH */
      sprintf(title,"rematch_resolution %f, rematch_bins %d, gsc_bin_index %d",rematch_resolution,rematch_bins,gsc_bin_index);
#endif /* LIMIT_REMATCH */
      break;
    }
  }




	plAlloc2dGrid(&scale, nx+1, ny+1);

  for (curCount = 0; curCount < totalCount; curCount++) {
    if (sequenceTable[curCount] == SEQUENCE_GROUPCOUNT) {
      sequenceGroupCount++;
    }
    if (sequenceTable[curCount] == SEQUENCE_INTEREST) {
      sequenceInterestCount++;
    }

    if (sequenceTable[curCount] != SEQUENCE_IMAGECOUNT) {
        continue;
    }
    sequenceImageCount++;
    intCount = raTable[curCount]+0.1;
    ix = intCount;
    if ((ix < 0) || (ix > nx)) {
      printf("ERROR: line %d  curCount %d raTable %f\n",__LINE__,curCount,raTable[curCount]);
      continue;
    }
    intCount = decTable[curCount]+0.1;

    if ((ix < 0) || (ix > nx)) {
      printf("ERROR: line %d  curCount %d decTable %f\n",__LINE__,curCount,decTable[curCount]);
      continue;
    }

    iy = intCount;
#if 1
    plateCount = countTable[curCount];
#else
    plateCount = log(10*(countTable[curCount]+0.1));
#endif
    if ((plateCount > 0 ) && (minCount > plateCount)) {
      minCount = plateCount;
    }
    if (maxCount < plateCount) {
      maxCount = plateCount;
    }
    if (lowerLimit > -1)  {
      if (plateCount < lowerLimit) {
        plateCount = 0;
      }
    }
    if (upperLimit > -1) {
      if (plateCount >= upperLimit) {
        plateCount = upperLimit-1;
      }
    }
    scale[ix][iy] = plateCount;
#if 1
    printf("line %d index %d ix %d iy %d plateCount %f\n",__LINE__,curCount,ix,iy,plateCount);
#endif
  }



	if (lowerLimit > -1) {
		minCountx = lowerLimit;
	} else {
		minCountx = minCount;
	}
	if (minCountx <= 0) {
		minCountx = 1;
	}
	if (upperLimit > -1) {
		maxCountx = upperLimit;
	} else {
		maxCountx = maxCount;
	}
	
	printf("line %d minCount %f minCountx %f, maxCount %f, maxCountx %f\n",__LINE__,minCount,minCountx,maxCount,maxCountx);

	if ((minCountx == 1) && (maxCountx == 1)) {
		minCountx = 0;
	}

#ifdef DO_SHADES
	dscale = (1.0*(maxCountx-minCountx))/(1.0*(NSHADES-1));
	if (dscale <= 0) {
		printf("ERROR: dscale (1) is %f for %s\n",dscale,title);
		return;
	}

	for (index = 0; index < NSHADES; index++) {
		clevel[index] = minCountx + (index * dscale);
#if 0
		printf("index %d clevel %f\n",index,clevel[index]);
#endif
	}
	dscale = (1.0*(maxCountx-minCountx))/(1.0*(NCOLORBAR-1));
	if (dscale <= 0) {
		printf("ERROR: dscale (2) is %f for %s\n",dscale,title);
		return;
	}
	for (index = 0; index < NCOLORBAR; index++) {
		shedge[index] = minCountx + (index * dscale);
#if 0
		printf("index %d shedge %f\n",index,shedge[index]);
#endif
	}
#endif /* DO_SHADES */

  maxlevel = maxCountx+5;
  contlevel = (PLFLT*)calloc(maxlevel,sizeof(PLFLT));
  if (contlevel == NULL) {
    printf("ERROR: line %d allocation failure\n",__LINE__);
    exit(-1);
  }
  /* Select isobars above and below each point */
  for (curCount = 0; curCount < totalCount; curCount++) {
    if (sequenceTable[curCount] != SEQUENCE_IMAGECOUNT) {
        continue;
    }
#if 0
    indexlevel = countTable[curCount];
    if ((indexlevel >= 0) && (indexlevel < maxlevel)) {
      if  (contlevel[indexlevel] == 0) {
        contlevel[indexlevel] = indexlevel;
      }
    }
    indexlevel = countTable[curCount]+1;
    if ((indexlevel >= 0) && (indexlevel < maxlevel)) {
      if  (contlevel[indexlevel] == 0) {
        contlevel[indexlevel] = indexlevel;
      }
    }
#endif
    indexlevel = countTable[curCount]-1;
    if ((indexlevel >= 0) && (indexlevel < maxlevel)) {
      if  (contlevel[indexlevel] == 0) {
        contlevel[indexlevel] = indexlevel;
      }
    }
  }
  /* Now compress the isobar table */
  for (indexlevel = 0; indexlevel < maxlevel; indexlevel++) {
    if (contlevel[indexlevel] != 0) {
      contlevel[nlevel] = contlevel[indexlevel];
#if 1
      printf("line %d indexlevel %d nlevel %d contlevel %f\n",__LINE__,indexlevel,nlevel,contlevel[nlevel]);
#endif
      nlevel++;
    }
  }

  /* No do the group coloring */
  rematchGroupTable = (int*)calloc(sequenceGroupCount,sizeof(int));
  if (rematchGroupTable == NULL) {
    printf("ERROR allocating the rematchGroupTable\n");
    exit(-1);
  }
  for (curCount = 0; curCount < totalCount; curCount++) {
    if (sequenceTable[curCount] != SEQUENCE_GROUPCOUNT) {
      continue;
    }
    for (groupIndex = 0; groupIndex < numGroups; groupIndex++) {
      if (rematchGroupTable[groupIndex] == countTable[curCount]) {
        break;
      }
    }
    if (groupIndex >= numGroups) {
      rematchGroupTable[numGroups] = countTable[curCount];
      numGroups++;
    }
  }
  qsort((void*)rematchGroupTable,numGroups,sizeof(int),IntCompare);
  printf("Found %d groups in %d bins\n",numGroups,sequenceGroupCount);
#if 1
  for (groupIndex = 0; groupIndex < numGroups; groupIndex++) {
    printf("line %d groupIndex %3d group %3d\n",__LINE__,groupIndex,rematchGroupTable[groupIndex]);
  }
#endif

#if 1
	plsdev("png");
  plsetopt("-o",outfile);
  plsetopt("-geometry","2880x2160");
#endif
#if 0
	plsdev("psc");
  plsetopt("-o","loshammer.psc");
#endif
#if 0
	plsdev("xwin");
#endif
#if 0
	plsori(1);  /* Set to portrait */
#endif
  plscolbg(255,255,255);	/* Force the background colour to white */  
  plscol0(1, 0,0,0);		/* Force the foreground colour to black */
  plscol0(15,255,0,0);		/* Move red to 15 */
	plinit();
  xmin = 0.0;
  xmax = nx;
  ymin = 0.0;
  ymax = ny;


  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax,just, 2); 
	plwind(xmin,xmax,ymin,ymax);
	plgspa(&xminmm,&xmaxmm,&yminmm,&ymaxmm);
	printf("line %d xminmm %f xmaxmm %f yminmm %f ymaxmm %f\n",__LINE__,xminmm,xmaxmm,yminmm,ymaxmm);
	/*	plwind(xmin,xmax,ymin,ymax); */




#if 0
#if 1 /* unsaturated table */
	cpoint[0] = 0.0; r[0] = 1.0; g[0] = 0.5; b[0] = 0.5;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 0.5; g[2] = 0.5; b[2] = 1.0;
#endif
#if 0 /* saturated table */
	cpoint[0] = 0.0; r[0] = 1.0; g[0] = 0.0; b[0] = 0.0;
  cpoint[1] = 0.5; r[1] = 0.0; g[1] = 1.0; b[1] = 0.0;
  cpoint[2] = 1.0; r[2] = 0.0; g[2] = 0.0; b[2] = 1.0;
#endif
#if 0 /* Inverted table */
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;

#endif
  plscmap1l(1, 3, cpoint, r, g, b, NULL);
#endif
  plscmap1l( 1, 5, pos, red, green, blue, NULL );
  for (curCount = 0; curCount < totalCount; curCount++) {
    if (sequenceTable[curCount] != SEQUENCE_GROUPCOUNT) {
      continue;
    }
    for (groupIndex = 0; groupIndex < numGroups; groupIndex++) {
      if (rematchGroupTable[groupIndex] == countTable[curCount]) {
        break;
      }
    }
    if (groupIndex >= numGroups) {
      printf("ERROR line %d\n",__LINE__);
      exit(-1);
    }
    plcol1((1.0*groupIndex)/(1.0*numGroups));
    plfbox(raTable[curCount],decTable[curCount]);

  }



	printf("line %d xmin %f xmax %f ymin %f ymax %f\n",__LINE__,xmin,xmax,ymin,ymax);
#ifdef DO_SHADES
  plshades((const PLFLT * const *)scale,nx,ny,NULL,xmin,xmax,ymin,ymax,clevel,NSHADES,1,0,0,plfill,0,grid_pltr3,NULL);
#endif /* DO_SHADES */
  plcol0(7); /* use grey for the contours */
  plcont((const PLFLT * const *)scale,nx+1,ny+1,1,nx,1,ny,contlevel,nlevel,grid_pltr3,NULL);


	plcol0(1);
#if 1
  imageTotalTable = (PIMAGETOTAL)calloc(sequenceInterestCount,sizeof(IMAGETOTAL));

  for (curCount = 0; curCount < totalCount; curCount++) {
    if (sequenceTable[curCount] != SEQUENCE_INTEREST) {
      continue;
    }
#if 1
    printf("line %d index %d ra %f dec %f group %f\n",__LINE__,curCount,raTable[curCount],decTable[curCount],countTable[curCount]);
#endif
    for (imageTotalIndex = 0; imageTotalIndex < imageTotalCount; imageTotalIndex++) {
      pImageTotal = &imageTotalTable[imageTotalIndex];
      if (pImageTotal->group == countTable[curCount]) {
        break;
      }
    }
    if (imageTotalIndex >= imageTotalCount) {
      if (imageTotalCount >= sequenceInterestCount) {
        printf("ERROR: line %d too many images\n",__LINE__);
        exit(-1);
      }
      pImageTotal = &imageTotalTable[imageTotalCount];
      pImageTotal->group = countTable[curCount];
      imageTotalIndex = imageTotalCount;
      imageTotalCount++;
    }
    pImageTotal->count++;
    if (pImageTotal->count > maxImageCount) {
      maxImageCount = pImageTotal->count;
      maxImageIndex = imageTotalIndex;
    }
  }
  pImageTotal = &imageTotalTable[maxImageIndex];
  printf("Maximum group number %d count %d total groups %d total images %d\n",pImageTotal->group,pImageTotal->count,imageTotalCount,sequenceInterestCount);
  for (curCount = 0; curCount < totalCount; curCount++) {
    if (sequenceTable[curCount] != SEQUENCE_INTEREST) {
      continue;
    }
    ra = raTable[curCount]-0.5;
    dec = decTable[curCount]-0.5;
     if (countTable[curCount] ==  pImageTotal->group) {
       plcol0(9); /* Blue for in group */
     } else {
       plcol0(15); /* Red for in group */
     }
       plpoin(1,&ra,&dec,12);

  }
#endif
  plcol0(1);
  plschr( 0.0, 0.2 );
  for (curCount = 0; curCount < totalCount; curCount++) {
    if (sequenceTable[curCount] != SEQUENCE_IMAGECOUNT) {
        continue;
    }
#ifndef LIMIT_REMATCH
    sprintf(str,"%.0f",countTable[curCount]);
    plptex(raTable[curCount],decTable[curCount],1.0,0,0.5,str);
#endif /* LIMIT_REMATCH*/
  }
  plschr( 0.0, 1.0 );

	if (galacticFlag) {
		pllab("Galactic Longitude","Latitude",title);
	} else {
    sprintf(ralabel,"RA (arcsec/%.2f)",rematch_resolution);
    sprintf(declabel,"Dec (arcsec/%.2f)",rematch_resolution);
		pllab(ralabel,declabel,title);
	}
#ifdef USE_COLORBAR
	{
		PLFLT colorbar_width;
		PLFLT colorbar_height;
		PLINT cont_color = 0;
		PLINT cont_width = 0;
		// Smaller text
		plschr( 0.0, 0.6);
		// Small ticks on the vertical axis
		plsmaj( 0.0, 0.5 );
		plsmin( 0.0, 0.5 );
#ifndef plwidth  /* Fedora 18 version */
		plcolorbar( &colorbar_width, &colorbar_height,
								PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL |PL_COLORBAR_LABEL_TOP, 0, 
								0.001,0.0,   /* x,y */
								0.02, 0.875,  /* x_length,y_length */
								0, 1, 1,        /* bg_color,bb_color, bb_style */
								0.0, 0.0,       /* low_cap_color, high_cap_color */
								cont_color, cont_width, 0.0, 0, "bcvtm",colormap_scale,
								NCOLORBAR, shedge );
#else /* plwidth Fedora 20 version */
		{
			PLINT      n_labels     = 1;
			PLINT      label_opts[] = {
        PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL | PL_COLORBAR_LABEL_TOP,
			};
			const char *labels[] = {
        colormap_scale,
			};
			const char *axis_opts[] = {
        "bcvtm",
			};
			PLINT      num_values[1];
      PLFLT      *values[1];
			PLFLT      axis_ticks[1] = {
        0.0,
			};
			PLINT      axis_subticks[1] = {
        0,
			};


			values [0] = shedge;
			num_values[0] = NCOLORBAR;

			plcolorbar( &colorbar_width, &colorbar_height,
									PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL |PL_COLORBAR_LABEL_TOP, 0, 
									0.001,0.0,   /* x,y */
									0.02, 0.875,  /* x_length,y_length */
									0, 1, 1,        /* bg_color,bb_color, bb_style */
									0.0, 0.0,       /* low_cap_color, high_cap_color */
									cont_color, cont_width, n_labels, label_opts, 
									labels,                            /* Old arg 19 New arg 18 label */
									1,                                 /*            New arg 19 n_axes */
									axis_opts,                         /* Old arg 18 New arg 20 axis_opts*/
									axis_ticks,                        /*            New arg 21 ticks */
									axis_subticks,                     /*            New arg 22 sub_ticks */
								  num_values,                        /* Old arg 20 New arg 23 n_values */
									(const PLFLT * const *)values );   /* Old arg 21 New arg 24 values*/
		}

#endif /* plwidth */
		// Reset text and tick sizes
		plschr( 0.0, 1.0 );
		plsmaj( 0.0, 1.0 );
		plsmin( 0.0, 1.0 );
		

	}
#endif /* USE_COLORBAR */


	plend();
	plFree2dGrid(scale, nx, ny);
  if (contlevel != NULL) {
    free(contlevel);
  }
  if (rematchGroupTable != NULL) {
    free(rematchGroupTable);
  }
  if (imageTotalTable != NULL) {
    free(imageTotalTable);
  }


	return;


}
