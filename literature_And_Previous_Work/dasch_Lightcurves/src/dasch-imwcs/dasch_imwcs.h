// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#ifndef __DASCH_IMWCS_H__
#define __DASCH_IMWCS_H__

// DASCH customizations to upstream wcstools
#define CROP_BORDER
#define LIMIT_SHIFT
#define POLE_FIX

#include <stdio.h> // FILE needed by libwcs headers

#include <libwcs/wcs.h>

// imsetwcs.c:

#define SetWCSFITS dasch_SetWCSFITS
#define settolerance dasch_settolerance
#define setnxydec dasch_setnxydec
#define setirafout dasch_setirafout
#define setreflim dasch_setreflim
#define setfitwcs dasch_setfitwcs
#define setfitplate dasch_setfitplate
#define setminstars dasch_setminstars
#define setnofit dasch_setnofit
#define setfrac dasch_setfrac
#define setmaxcat dasch_setmaxcat
#define setiterate dasch_setiterate
#define setnfiterate dasch_setnfiterate
#define setiteratet dasch_setiteratet
#define setrecenter dasch_setrecenter
#define setsortmag dasch_setsortmag
#define setmagfit dasch_setmagfit
#define settabkrw dasch_settabkrw
#define setcropborder dasch_setcropborder
#define setsquareborder dasch_setsquareborder

extern int dasch_SetWCSFITS(
    char *filename, // image file name
    char *header, // FITS header
    char *image,  // image pixels
    char *refcatname,  // name of reference catalog
    int verbose
);

extern void dasch_settolerance(double tol);
extern void dasch_setnxydec(int ndex);
extern void dasch_setirafout(void);
extern void dasch_setreflim(double lim1, double lim2);
extern void dasch_setfitwcs(int wfit);
extern void dasch_setfitplate(int nc);
extern void dasch_setminstars(int minstars);
extern void dasch_setnofit(void);
extern void dasch_setfrac(double frac0);
extern void dasch_setmaxcat(int ncat);
extern void dasch_setiterate(int iter);
extern void dasch_setnfiterate(int iter);
extern void dasch_setiteratet(int iter);
extern void dasch_setrecenter(int recenter);
extern void dasch_setsortmag(int imag);
extern void dasch_setmagfit(void);
extern void dasch_settabkrw(char *keyword0);
extern void dasch_setcropborder(int borderpix);
extern void dasch_setsquareborder(int squareborder);

// matchstar.c:

#define StarMatch dasch_StarMatch
#define getnfit dasch_getnfit
#define setnfit dasch_setnfit
#define setminbin dasch_setminbin
#define setmaxshift dasch_setmaxshift

extern int dasch_StarMatch(
    int ns, // number of image stars
    double *sx, // image star pixel coordinates
    double *sy,
    int refcat, // reference catalog code
    int ng, // number of reference stars (refstars)
    double *gnum, // refstar catalog numbers
    double *gra, // refstar sky coordinates in degrees
    double *gdec,
    int *goff, // refstar offscale flags
    double *gx, // refstar pixel coordinates
    double *gy,
    double tol, // tolerance in pixels
    struct WorldCoor *wcs, // fit returned here
    int debug
);
extern int dasch_getnfit(void);
extern void dasch_setnfit(int);
extern void dasch_setminbin(int);
#ifdef LIMIT_SHIFT
extern void dasch_setmaxshift(int);
#endif /* LIMIT_SHIFT */

// Prototypes that libwcs doesn't provide in any headers. Usually this would
// mean that a function is private and shouldn't be used, but libwcs is super
// sloppy about this stuff.

// findstar.c:
extern int FindStars(
    char *header, // FITS header text
    char *image, // image pixel data
    double **xa, // X and Y coordinates of stars (returned arrays)
    double **ya,
    double **ba, // fluxes of stars in counts (returned array)
    int **pa, // peak counts of stars in counts (returned array)
    int verbose, // print each star's position
    int zap // if 1, "set star to background after reading"
);
extern char *getimcat(void);
extern void setparm(char *);
extern void setrotate(int);
extern void setmirror(int);
extern void setbmin(double);
extern void setstarsig(double);
extern void setimcat(char *);

// imgetwcs.c:
extern struct WorldCoor *GetFITSWCS(
    char *filename,
    char *header, // FITS header data
    int verbose,
    double *cra, // center coordinate in degrees (returned)
    double *cdec,
    double *dra, // ra/dec half-width in degrees (returned)
    double *ddec,
    double *secpix, // arcseconds per pixel (returned)
    int *wp, // image width/height in pixels (returned)
    int *hp,
    int *sysout, // coordinate system to use (0=use image's; in/out)
    double *eqout /// equinox to use (0=use image's; in/out)
);
extern void setdcenter(double, double);
extern void setrot(double);
extern void setsecpix(double);
extern void setsecpix2(double);
extern void setsys(int);
extern void getcenter(double *ra, double *dec);
extern void setcenter(char *, char *);
extern void getrefpix(double *x, double *y);
extern void setrefpix(double, double);

// imrotate.c:
extern char *RotFITS(
    char *pathname,
    char *header,
    char *image0,
    int xshift,
    int yshift,
    int rotate,
    int mirror,
    int bitpix2,
    int rotwcs,
    int verbose
);

// imsetwcs.c:
extern int SetWCSFITS(
    char *filename,
    char *header,
    char *image,
    char *refcatname,
    int verbose
);
extern double getsecpix(void);
extern void setnofit(void);

// platefit.c:
extern int FitPlate(
    struct WorldCoor *wcs,
    double *x, // image WCS coordinates
    double *y,
    double *x1, // image pixel coordinates
    double *y1,
    int np, // number of points to fit
    int ncoeff0, // order of polynomial terms in x and y
    int debug
);

#endif