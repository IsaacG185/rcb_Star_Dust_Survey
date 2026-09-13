// This is a customized version of WCSTools imwcs support code, derived from:
//
//   IMWCS WCSTools 3.9.7, 11 August 2022, Jessica Mink (jmink@cfa.harvard.edu)
//
// It inherits the license of the above: GPLv2 or later.

#include "dasch_imwcs.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <libwcs/lwcs.h>
#include <libwcs/wcscat.h>

#define NPAR 8
#define NPAR1 9

#define ABS(a) ((a) < 0 ? (-(a)) : (a))

static void wcs_amoeba(struct WorldCoor *wcs0);
static int dasch_ParamFit(int nbin);
static double wcs_chisqr(double *v, int iter);
static void dasch_amoeba(
    double **p,
    double y[],
    int ndim,
    double ftol,
    int itmax,
    double (*funk)(double *, int),
    int *nfunk
);
static double amotry (
    double **p,
    double *y,
    double *psum,
    int ndim,
    double (*funk)(double *, int),
    int ihi,
    int *nfunk,
    double fac
);

static struct WorldCoor *wcsf;

/* Statics used by the chisqr evaluator */
static double *sx_p;
static double *sy_p;
static double *gra_p;
static double *gdec_p;
static double xref_p, yref_p;
static double xrefpix, yrefpix;
static int nbin_p;
static int nfit; /* Number of parameters to fit */
static int pfit0 = 0; /* List of parameters to fit, 1 per digit */
static int cdfit = 0; /* 1 if CD matrix has been fit */
static int minbin = 2; /* Minimum number of coincidence hits needed */
static int minmatch0 = MINMATCH; /* matches to drop out of loop */
static int nitmax0 = NMAX;  /* max iterations to stop fit */
static int vfit[NPAR1]; /* Parameters being fit: index to value vector
                                1= RA,    2= Dec,
                                3= X plate scale, 4= Y plate scale
                                5= rotation,   6= second rotation (skew),
                                7= optical axis X,8= optical axis Y */

#ifdef LIMIT_SHIFT
static int maxshift0 = 0;
#endif /* LIMIT_SHIFT */

/* Find shift, scale, and rotation of image stars to best-match reference stars
 * Get best match by finding which offsets between pairs of s's and g's
 * work for the most other pairs of s's and g's
 * N.B. we assume rotation will be "small enough" so that initial guesses can
 *   be done using just shifts.
 * Return count of total coincidences found, else 0 if none or -1 if trouble.
 */

int
dasch_StarMatch (
    int ns,  /* Number of image stars */
    double *sx,  /* Image star X coordinates in pixels */
    double *sy,  /* Image star Y coordinates in pixels */
    int refcat,  /* Reference Catalog code */
    int ng,  /* Number of reference stars */
    double *gnum,  /* Reference star catalog numbers */
    double *gra,  /* Reference star right ascensions in degrees */
    double *gdec,  /* Reference star right ascensions in degrees */
    int *goff,  /* Reference star offscale flags */
    double *gx,  /* Reference star X coordinates in pixels */
    double *gy,  /* Reference star Y coordinates in pixels */
    double tol,  /* +/- this many pixels is a hit */
    struct WorldCoor *wcs, /* World coordinate structure (fit returned) */
    int debug
) {
    double dx, bestdx, dxi;
    double dy, bestdy, dyi;
    double dx2, dy2, dxy, dxys, dxs, dys, dxsum, dysum;
    double *mx, *my, *mxy;
    int nmatch;
    int s, g, si, gi, igs;
    int nbin;
    double *sbx, *sby; /* malloced array of s stars in best bin */
    double *gbra, *gbdec; /* malloced array of g stars in best bin */
    int peaks[NPEAKS+1]; /* history of bin counts */
    int dxpeaks[NPEAKS+1], dypeaks[NPEAKS+1]; /* history of dx/dy at peaks */
    int npeaks;  /* entries in use in peaks[] */
    int maxnbin, i, nmatchd;
    int minmatch;
    int *is, *ig, *ibs, *ibg;
    char rastr[32], decstr[32];
    double xref0, yref0, xinc0, yinc0, rot0, xrefpix0, yrefpix0, cd0[4];
    int pfit;  /* List of parameters to fit, 1 per digit */
    char vpar[16]; /* List of parameters to fit */
    char *vi;
    char vc;
    double tol2 = tol * tol;

    /* Set minimum number of matches between image and reference stars to fit */
    if (ns > ng) {
        minmatch = 0.5 * ng;

        if (minmatch > minmatch0)
            minmatch = 0.25 * ng;

        if (minmatch > minmatch0)
            minmatch = minmatch0;
    } else {
        minmatch = 0.5 * ns;

        if (minmatch > minmatch0)
            minmatch = 0.25 * ns;

        if (minmatch > minmatch0)
            minmatch = minmatch0;
    }

    /* Set maximum number of matches and allocate match indices */
    if (ng > ns)
        maxnbin = (int) ((double) ng * 1.25);
    else
        maxnbin = (int) ((double) ns * 1.25);

    if (debug)
        fprintf(
            stderr,
            "Match history: nim=%d nref=%d tol=%3.0f minbin=%d minmatch=%d):\n",
            ns,
            ng,
            tol,
            minbin,
            minmatch
        );

    /* Allocate arrays in which to save match information */
    is = (int *) calloc (maxnbin, sizeof(int));
    ig = (int *) calloc (maxnbin, sizeof(int));
    ibs = (int *) calloc (maxnbin, sizeof(int));
    ibg = (int *) calloc (maxnbin, sizeof(int));

    /* Try matching stars using the current WCS first */
    nmatch = 0;
    bestdx = 0.0;
    bestdy = 0.0;
    dxsum = 0.0;
    dysum = 0.0;
    dxs = 0.0;
    dys = 0.0;

    /* Vote for closest match */
    mx = (double *) calloc (maxnbin, sizeof(double));
    my = (double *) calloc (maxnbin, sizeof(double));
    mxy = (double *) calloc (maxnbin, sizeof(double));

    /* Loop through image stars */
    for (s = 0; s < ns; s++) {
        dxys = tol2;
        igs = -1;

        /* Loop through reference catalog stars */
        for (g = 0; g < ng; g++) {
            /* Try reference catalog star only if it is on the image */
            dx = gx[g] - sx[s];
            dy = gy[g] - sy[s];
            dx2 = dx * dx;
            dy2 = dy * dy;
            dxy = dx2 + dy2;

            /* Check offset less than tolerance or this star's closest match */
            if (dxy < dxys) {
                dxys = dxy;
                dxs = dx;
                dys = dy;
                igs = g;
                ibs[nmatch] = s;
                ibg[nmatch] = g;
            }
        }

        /* If a match was found */
        if (igs > -1) {
            int report = 0;

            /* if new match is closer than old match, replace it */
            if (mxy[igs] > 0.0) {
                if (dxy < mxy[igs]) {
                    dxsum = dxsum - mx[igs];
                    dysum = dysum - my[igs];
                    dxsum = dxsum + dxs;
                    dysum = dysum + dys;
                    report = 1;
                }
            } else {
                /* If not matched before, use new match */
                dxsum = dxsum + dxs;
                dysum = dysum + dys;
                nmatch++;
                mx[igs] = dxs;
                my[igs] = dys;
                mxy[igs] = dxy;
                report = 1;
            }

            if (report) {
#if 0
                if (debug) {
                    char numstr[32];

                    CatNum (refcat, 10, 0, gnum[ibg[nmatch]], numstr);
                    ra2str (rastr, 31, gra[ibg[nmatch]], 3);
                    dec2str (decstr, 31, gdec[ibg[nmatch]], 2);
                    fprintf(
                        stderr,
                        " %3d %s %s %s %7.2f %7.2f %7.2f %7.2f %5.2f %5.2f %5.2f\n",
                        nmatch,
                        numstr,
                        rastr,
                        decstr,
                        gx[ibg[nmatch]],
                        gy[ibg[nmatch]],
                        sx[ibs[nmatch]],
                        sy[ibs[nmatch]],
                        dxs,
                        dys,
                        sqrt(dxys)
                    );
                }
#endif
            }
        }
    }

    free (mxy);
    free (mx);
    free (my);

    /* If we found enough matches, we can proceed with this offset */
    if (nmatch >= minmatch) {
        bestdx = dxsum / (double) nmatch;
        bestdy = dysum / (double) nmatch;

        if (debug)
            fprintf(stderr, "%d matches found at mean offset %6.3f %6.3f\n", nmatch, bestdx, bestdy);
    }

    /* Otherwise, we will look for a coarse alignment assuming no additional rotation.
     * This will allow us to collect a set of stars that correspond and
     * establish an initial guess of the solution.
     */

    if (nmatch < minmatch) {
        if (debug)
            fprintf (stderr, "%d matches found  less than %d minimum\n", nmatch, minmatch);

        npeaks = 0;
        nmatch = 0;

        for (i = 0; i < NPEAKS; i++) {
            peaks[i] = 0;
            dxpeaks[i] = 0;
            dypeaks[i] = 0;
        }

        bestdx = 0.0;
        bestdy = 0.0;

        for (s = 0; s < ns; s++) {
            for (g = 0; g < ng; g++) {
                dx = gx[g] - sx[s];
                dy = gy[g] - sy[s];

#ifdef LIMIT_SHIFT
                if (maxshift0 > 0) {
                    /* Skip the loop if the offset is too large */
                    if (dx >= maxshift0 || dx <= -maxshift0 || dy >= maxshift0 || dy <= -maxshift0) {
                        continue;
                    }
                }
#endif /* LIMIT_SHIFT */

                nbin = 0;

                for (gi = 0; gi < ng; gi++) {
                    for (si = 0; si < ns; si++) {
                        dxi = gx[gi] - sx[si] - dx;
                        if (dxi < 0)
                            dxi = -dxi;

                        dyi = gy[gi] - sy[si] - dy;
                        if (dyi < 0)
                            dyi = -dyi;

                        if (dxi <= tol && dyi <= tol) {
                            is[nbin] = si;
                            ig[nbin] = gi;
                            nbin++;
                        }
                    }
                }

                if (nbin > 1 && nbin >= nmatch) {
                    int i;

                    nmatch = nbin;
                    bestdx = (double) dx;
                    bestdy = (double) dy;

                    for (i = 0; i < nbin; i++) {
                        ibs[i] = is[i];
                        ibg[i] = ig[i];
                    }

                    /* keep last NPEAKS nmatchs, dx and dy;
                     * put newest first in arrays */
                    if (npeaks > 0) {
                        for (i = npeaks; i > 0; i--) {
                            peaks[i] = peaks[i-1];
                            dxpeaks[i] = dxpeaks[i-1];
                            dypeaks[i] = dypeaks[i-1];
                        }
                    }

                    peaks[0] = nmatch;

                    if (bestdx > 0.0)
                        dxpeaks[0] = (int) (bestdx + 0.5);
                    else
                        dxpeaks[0] = (int) (bestdx - 0.5);

                    if (bestdy > 0)
                        dypeaks[0] = (int) (bestdy + 0.5);
                    else
                        dypeaks[0] = (int) (bestdy - 0.5);

                    if (npeaks < NPEAKS)
                        npeaks++;
#if 0
                    if (debug)
                        fprintf(
                            stderr,
                            "%d: %d/%d matches at image %d cat %d: dx= %d dy= %d\n",
                            npeaks,
                            nmatch,
                            minmatch,
                            s,
                            g,
                            dxpeaks[0],
                            dypeaks[0]
                        );
#endif
                }

                if (nmatch > minmatch)
                    break;
            }

            if (nmatch > minmatch)
                break;
        }

        /* peak is broad */
        if (npeaks < 2 || peaks[1] == peaks[0]) {
            if (debug)
                fprintf(
                    stderr,
                    "  Broad peak of %d bins at dx=%.0f dy=%.0f\n",
                    peaks[0],
                    bestdx,
                    bestdy
                );
        }
    }

    /* too few hits */
    if (nmatch < minbin) {
        fprintf(
            stderr,
            "StarMatch: not enough matches to attempt fit (got %d, need %d)\n",
            nmatch,
            minbin
        );
        return nmatch;
    }

    /* Get X and Y coordinates of matches from best binning */
    nmatchd = nmatch * sizeof (double);

    if (!(sbx = (double *) malloc (nmatchd)))
        fprintf (stderr," Could not allocate %d bytes for SBX\n", nmatchd);
    if (!(sby = (double *) malloc (nmatchd)))
        fprintf (stderr," Could not allocate %d bytes for SBY\n", nmatchd);
    if (!(gbra = (double *) malloc (nmatchd)))
        fprintf (stderr," Could not allocate %d bytes for GBRA\n", nmatchd);
    if (!(gbdec = (double *) malloc (nmatchd)))
        fprintf (stderr," Could not allocate %d bytes for GBDEC\n", nmatchd);

    for (i = 0; i < nmatch; i++) {
        sbx[i] = sx[ibs[i]];
        sby[i] = sy[ibs[i]];
        gbra[i] = gra[ibg[i]];
        gbdec[i] = gdec[ibg[i]];
    }

    /* Reset image center based on star matching */
    wcs->xref = wcs->xref + (bestdx * wcs->xinc);
    if (wcs->xref < 0.0)
        wcs->xref = 360.0 + wcs->xref;

    wcs->yref = wcs->yref + (bestdy * wcs->yinc);

#ifdef POLE_FIX
    if (wcs->yref > 90.0) {
        wcs->yref = 90.0;
    } else if (wcs->yref < -90.0) {
        wcs->yref = -90.0;
    }
#endif /* POLE_FIX */

    /* Fit WCS to matched stars */

    /* Provide non-parametric access to the star lists */
    sx_p = sbx;
    sy_p = sby;
    gra_p = gbra;
    gdec_p = gbdec;
    xref_p = wcs->xref;
    yref_p = wcs->yref;
    xrefpix = wcs->xrefpix;
    yrefpix = wcs->yrefpix;
    nbin_p = nmatch;

    /* Number of parameters to fit from command line or number of matches */
    pfit = dasch_ParamFit (nmatch);

    /* Get parameters to fit from digits of pfit */
    sprintf (vpar, "%d", pfit);
    nfit = 0;
    vfit[0] = -1;

    for (i = 1; i < NPAR1; i++) {
        vc = i + 48;
        vi = strchr (vpar, vc);

        if (vi != NULL) {
            vfit[i] = vi - vpar;
            nfit++;
        } else
            vfit[i] = -1;
    }

    /* Set initial guesses for parameters which are being fit */
    xref0 = wcs->xref;
    yref0 = wcs->yref;
    xinc0 = wcs->xinc;
    yinc0 = wcs->yinc;
    rot0 = wcs->rot;
    xrefpix0 = wcs->xrefpix;
    yrefpix0 = wcs->yrefpix;
    cd0[0] = wcs->cd[0];
    cd0[1] = wcs->cd[1];
    cd0[2] = wcs->cd[2];
    cd0[3] = wcs->cd[3];

    if (vfit[6] > -1)
        cdfit = 1;
    else
        cdfit = 0;

    /* Fit image star coordinates to reference star positions */
    wcs_amoeba (wcs);

    if (debug) {
        fprintf (stderr,"\nAmoeba fit:\n");
        ra2str (rastr, 31, xref0, 3);
        dec2str (decstr, 31, yref0, 2);
        fprintf (stderr,"   initial guess:\n");

        if (vfit[6] > -1)
            fprintf(
                stderr,
                " cra= %s cdec= %s cd = %9.7f,%9.7f,%9.7f,%9.7f ",
                rastr,
                decstr,
                cd0[0],
                cd0[1],
                cd0[2],
                cd0[3]
            );
        else
            fprintf(
                stderr,
                " cra= %s cdec= %s del=%7.4f,%7.4f rot=%7.4f ",
                rastr,
                decstr,
                xinc0 * 3600.0,
                yinc0 * 3600.0,
                rot0
            );

        fprintf (stderr,"(%8.2f,%8.2f\n", xrefpix0, yrefpix0);
        ra2str (rastr, 31, wcs->xref, 3);
        dec2str (decstr, 31, wcs->yref, 2);
        fprintf (stderr,"\nfirst solution:\n");

        if (vfit[6] > -1)
            fprintf(
                stderr,
                " cra= %s cdec= %s cd = %9.7f,%9.7f,%9.7f,%9.7f ",
                rastr,
                decstr,
                wcs->cd[0],
                wcs->cd[1],
                wcs->cd[2],
                wcs->cd[3]
            );
        else
            fprintf(
                stderr,
                " cra= %s cdec= %s del=%7.4f,%7.4f rot=%7.4f ",
                rastr,
                decstr,
                3600.0 * wcs->xinc,
                3600.0 * wcs->yinc,
                wcs->rot
            );

        fprintf (stderr,"(%8.2f,%8.2f)\n", wcs->xrefpix, wcs->yrefpix);
    }

    free (sbx);
    free (sby);
    free (gbra);
    free (gbdec);
    free (is);
    free (ig);
    free (ibs);
    free (ibg);
    return nmatch;
}


static int
dasch_ParamFit (int nbin)
{
    int pfit;

    if (pfit0 != 0) {
        if (pfit0 < 3)
            pfit = 12;
        else if (pfit0 == 3) /* Fit center and plate scale */
            pfit = 123;
        else if (pfit0 == 4) /* Fit center, plate scale, rotation */
            pfit = 1235;
        else if (pfit0 == 5) /* Fit center, x&y plate scales, rotation */
            pfit = 12345;
        else if (pfit0 == 6) /* Fit center, x&y plate scales, x&y rotations */
            pfit = 123456;
        else if (pfit0 == 7) /* Fit center, x&y plate scales, rotation, refpix */
            pfit = 1234578;
        else if (pfit0 == 8) /* Fit center, x&y plate scales, x&y rotation, refpix */
            pfit = 12345678;
        else
            pfit = pfit0;
    } else if (nbin < 4)
        pfit = 12;
    else if (nbin < 6)
        pfit = 123;
    else
        pfit = 12345;

    return pfit;
}


/* Set up the necessary temp arrays and call the amoeba() multivariate solver */
static void
wcs_amoeba (struct WorldCoor *wcs0)
{
    double *p[NPAR1];      /* used as p[NPAR1][NPAR] */
    double vguess[NPAR], vp[NPAR], vdiff[NPAR];
    double p0[NPAR], p1[NPAR], p2[NPAR], p3[NPAR], p4[NPAR],
           p5[NPAR], p6[NPAR], p7[NPAR], p8[NPAR]; /* used as px[0..NPAR-1] */
    double y[NPAR1];      /* used as y[1..NPAR] */
    double sumx, sumy, sumr;
    int iter;
    int i, j;
    int nfit1;
    char rastr[32],decstr[32];
    int nitmax;

    nitmax = nitmax0;
    if (nfit > NPAR)
        nfit = NPAR;
    nfit1 = nfit + 1;
    wcsf = wcs0;

    for (i = 0; i < NPAR; i++) {
        vguess[i] = 0.0;
        vdiff[i] = 0.0;
    }

    /* Optical axis center (RA and Dec degrees) */
    if (vfit[1] > -1) {
        vguess[vfit[1]] = 0.0;
        vdiff[vfit[1]] = 5.0 * wcsf->xinc;
    }

    if (vfit[2] > -1) {
        vguess[vfit[2]] = 0.0;
        vdiff[vfit[2]] = 5.0 * wcsf->yinc;
    }

    /* Second rotation about optical axis (degrees) -> CD matrix */
    if (vfit[6] > -1) {
        wcsf->rotmat = 1;
        vguess[vfit[3]] = wcsf->cd[0];
        vdiff[vfit[3]] = wcsf->xinc * 0.03;
        vguess[vfit[4]] = wcsf->cd[1];
        vdiff[vfit[4]] = wcsf->yinc * 0.03;
        vguess[vfit[5]] = wcsf->cd[2];
        vdiff[vfit[5]] = wcsf->xinc * 0.03;
        vguess[vfit[6]] = wcsf->cd[3];
        vdiff[vfit[6]] = wcsf->yinc * 0.03;
    } else {
        /* Plate scale at optical axis right ascension or both (degrees/pixel) */
        if (vfit[3] > -1) {
            vguess[vfit[3]] = wcsf->xinc;
            vdiff[vfit[3]] = wcsf->xinc * 0.03;
        }

        /* Plate scale in declination at optical axis (degrees/pixel) */
        if (vfit[4] > -1) {
            vguess[vfit[4]] = wcsf->yinc;
            vdiff[vfit[4]] = wcsf->yinc * 0.03;
        }

        /* Rotation about optical axis in degrees */
        if (vfit[5] > -1) {
            vguess[vfit[5]] = wcsf->rot;
            vdiff[vfit[5]] = 0.5;
        }
    }

    /* Reference pixel (optical axis) */
    if (vfit[7] > -1) {
        vguess[vfit[7]] = 0.0;
        vdiff[vfit[7]] = 10.0;
    }

    if (vfit[8] > -1) {
        vguess[vfit[8]] = 0.0;
        vdiff[vfit[8]] = 10.0;
    }

    /* Set up matrix of nfit+1 initial guesses.
    * The supplied guess, plus one for each parameter altered by a small amount
    */
    p[0] = p0;

    if (nfit > 0)
        p[1] = p1;

    if (nfit > 1)
        p[2] = p2;

    if (nfit > 2)
        p[3] = p3;

    if (nfit > 3)
        p[4] = p4;

    if (nfit > 4)
        p[5] = p5;

    if (nfit > 5)
        p[6] = p6;

    if (nfit > 6)
        p[7] = p7;

    if (nfit > 7)
        p[8] = p8;

    for (i = 0; i <= nfit; i++) {
        for (j = 0; j < nfit; j++)
            p[i][j] = vguess[j];

        if (i > 0 && i <= nfit)
            p[i][i-1] = vguess[i-1] + vdiff[i-1];

        y[i] = wcs_chisqr (p[i], -i);
    }

#if 0
    fprintf (stderr,"Before:\n");

    for (i = 0; i < nfit1; i++) {
        double xinc1, yinc1, rot, xrefpix1, yrefpix1;

        if (vfit[1] > -1)
            ra2str (rastr, 31, p[i][vfit[1]] + xref_p, 3);
        else
            ra2str (rastr, 31, wcsf->xref, 3);

        if (vfit[2] > -1)
            dec2str (decstr, 16, p[i][vfit[2]]+yref_p, 2);
        else
            dec2str (decstr, 16, wcsf->yref, 2);

        if (vfit[6] > -1) {
            double cd[4];

            cd[0] = p[i][vfit[3]];
            cd[1] = p[i][vfit[4]];
            cd[2] = p[i][vfit[5]];
            cd[3] = p[i][vfit[6]];
            fprintf(
                stderr,
                "%d: %s %s CD: %7.5f,%7.5f,%7.5f,%7.5f ",
                i,
                rastr,
                decstr,
                cd[0],
                cd[1],
                cd[2],
                cd[3]
            );
        } else {
            if (vfit[3] > -1)
                xinc1 = p[i][vfit[3]];
            else
                xinc1 = wcsf->xinc;

            if (vfit[4] > -1)
                yinc1 = p[i][vfit[4]];
            else if (vfit[3] > -1) {
                if (xinc1 < 0)
                    yinc1 = -xinc1;
                else
                    yinc1 = xinc1;
            } else
                yinc1 = wcsf->yinc;

            if (vfit[5] > -1)
                rot = p[i][vfit[5]];
            else
                rot = wcsf->rot;

            fprintf(
                stderr,
                "%d: %s %s del=%6.4f,%6.4f rot=%5.3f ",
                i,
                rastr,
                decstr,
                3600.0 * xinc1,
                3600.0 * yinc1,
                rot
            );
        }

        if (vfit[7] > -1)
            xrefpix1 = xrefpix + p[i][vfit[7]];
        else
            xrefpix1 = wcsf->xrefpix;

        if (vfit[8] > -1)
            yrefpix1 = yrefpix + p[i][vfit[8]];
        else
            yrefpix1 = wcsf->yrefpix;

        fprintf (stderr,"(%8.2f,%8.2f) y=%g\n", xrefpix1, yrefpix1, y[i]);
    }
#endif

    dasch_amoeba (p, y, nfit, FTOL, nitmax, wcs_chisqr, &iter);

#if 0
    fprintf (stderr,"\nAfter:\n");

    for (i = 0; i < nfit1; i++) {
        double xinc1, yinc1, rot, xrefpix1, yrefpix1;

        if (vfit[1] > -1)
            ra2str (rastr, 31, p[i][vfit[1]] + xref_p, 3);
        else
            ra2str (rastr, 31, wcsf->xref, 3);

        if (vfit[2] > -1)
            dec2str (decstr, 31, p[i][vfit[2]] + yref_p, 2);
        else
            dec2str (decstr, 31, wcsf->yref, 2);

        if (vfit[6] > -1) {
            double cd[4];

            cd[0] = p[i][vfit[3]];
            cd[1] = p[i][vfit[4]];
            cd[2] = p[i][vfit[5]];
            cd[3] = p[i][vfit[6]];
            fprintf (
                stderr,
                "%d: %s %s CD: %7.5f,%7.5f,%7.5f,%7.5f ",
                i,
                rastr,
                decstr,
                cd[0],
                cd[1],
                cd[2],
                cd[3]
            );
        } else {
            if (vfit[3] > -1)
                xinc1 = p[i][vfit[3]];
            else
                xinc1 = wcsf->xinc;

            if (vfit[4] > -1)
                yinc1 = p[i][vfit[4]];
            else if (vfit[3] > -1) {
                if (xinc1 < 0)
                    yinc1 = -xinc1;
                else
                    yinc1 = xinc1;
            } else
                yinc1 = wcsf->yinc;

            if (vfit[5] > -1)
                rot = p[i][vfit[5]];
            else
                rot = wcsf->rot;

            fprintf(
                stderr,
                "%d: %s %s del=%6.4f,%6.4f rot=%5.3f ",
                i,
                rastr,
                decstr,
                3600.0 * xinc1,
                3600.0 * yinc1,
                rot
            );
        }

        if (vfit[7] > -1)
            xrefpix1 = xrefpix + p[i][vfit[7]];
        else
            xrefpix1 = wcsf->xrefpix;

        if (vfit[8] > -1)
            yrefpix1 = yrefpix + p[i][vfit[8]];
        else
            yrefpix1 = wcsf->yrefpix;

        fprintf (stderr,"(%8.2f,%8.2f) y=%g\n", xrefpix1, yrefpix1, y[i]);
    }
#endif

    /* On return, all entries in p[1..NPAR] are within FTOL;
     * Return the average, though you could just pick the first one
     */
    for (j = 0; j < nfit; j++) {
        double sum = 0.0;

        for (i = 0; i < nfit1; i++)
            sum += p[i][j];

        vp[j] = sum / (double) nfit1;
    }

    if (vfit[1] > -1) {
        wcsf->xref = xref_p + vp[vfit[1]];
        if (wcsf->xref < 0.0)
            wcsf->xref = 360.0 + wcsf->xref;
    }

    if (vfit[2] > -1) {
        wcsf->yref = yref_p + vp[vfit[2]];
#ifdef POLE_FIX
        if (wcsf->yref > 90.0) {
            printf("matchstar.c line %d ERROR: adjusting yref %f\n", __LINE__, wcsf->yref);
            exit(1);
        } else if (wcsf->yref < -90.0) {
            printf("matchstar.c line %d ERROR: adjusting yref %f\n", __LINE__, wcsf->yref);
            exit(1);
        }
#endif /* POLE_FIX */
    }

    if (vfit[6] > -1) {
        wcsf->cd[0] = vp[vfit[3]];
        wcsf->cd[1] = vp[vfit[4]];
        wcsf->cd[2] = vp[vfit[5]];
        wcsf->cd[3] = vp[vfit[6]];
    } else {
        if (vfit[3] > -1)
            wcsf->xinc = vp[vfit[3]];

        if (vfit[4] > -1)
            wcsf->yinc = vp[vfit[4]];
        else if (vfit[3] > -1) {
            if (wcsf->xinc < 0)
                wcsf->yinc = -wcsf->xinc;
            else
                wcsf->yinc = wcsf->xinc;
        }

        if (vfit[5] > -1)
            wcsf->rot = vp[vfit[5]];
    }

    if (vfit[7] > -1)
        wcsf->xrefpix = xrefpix + vp[vfit[7]];

    if (vfit[8] > -1)
        wcsf->yrefpix = yrefpix + vp[vfit[8]];

    ra2str(rastr, 31, wcsf->xref, 3);
    dec2str(decstr, 31, wcsf->yref, 2);

    if (vfit[6] > -1)
        fprintf(
            stderr,
            "iter=%d cra=%s cdec=%s CD=%.7f,%.7f,%.7f,%.7f ",
            iter,
            rastr,
            decstr,
            wcsf->cd[0],
            wcsf->cd[1],
            wcsf->cd[2],
            wcsf->cd[3]
        );
    else
        fprintf(
            stderr,
            "iter=%d cra=%s cdec=%s del=%.4f,%.4f rot=%.4f ",
            iter,
            rastr,
            decstr,
            wcsf->xinc * 3600.0,
            wcsf->yinc * 3600.0,
            wcsf->rot
        );

    fprintf(stderr, "(%.2f, %.2f)\n", wcsf->xrefpix, wcsf->yrefpix);

    sumx = 0.0;
    sumy = 0.0;
    sumr = 0.0;

    for (i = 0; i < nbin_p; i++) {
        double mra, mdec, ex, ey, er, dra;

        pix2wcs (wcsf, sx_p[i], sy_p[i], &mra, &mdec);

        // Report reasonable values if we're encountering an RA wrap
        dra = mra - gra_p[i];

        while (dra > 180.)
            dra -= 360.;

        while (dra < -180.)
            dra += 360.;

        ex = 3600.0 * dra * cos(degrad(mdec));
        ey = 3600.0 * (mdec - gdec_p[i]);
        er = sqrt (ex * ex + ey * ey);
        sumx = sumx + ex;
        sumy = sumy + ey;
        sumr = sumr + er;

#if 0
        {
            char rastr[32], decstr[32];

            ra2str (rastr, 31, gra_p[i], 3);
            dec2str (decstr, 31, gdec_p[i], 2);
            fprintf (stderr,"%2d: c: %s %s ", i+1, rastr, decstr);
            ra2str (rastr, 31, mra, 3);
            dec2str (decstr, 31, mdec, 2);
            fprintf (
                stderr,
                "i: %s %s %6.3f %6.3f %6.3f\n",
                rastr,
                decstr,
                3600.0 * ex,
                3600.0 * ey,
                3600.0 * sqrt(ex*ex + ey*ey));
        }
#endif
    }

    sumx = sumx / (double) nbin_p;
    sumy = sumy / (double) nbin_p;
    sumr = sumr / (double) nbin_p;
    fprintf(stderr, "mean dra=%.3f ddec=%.3f dr=%.3f\n", sumx, sumy, sumr);
}


/* Compute the chisqr of the vector v, where
 * v[0]=cra, v[1]=cdec, v[2]=ra deg/pix, v[3]=dec deg/pix,
 * v[4]=rotation, v[5]=2nd rotation->CD matrix, v[6]=ref x, and v[7] = ref y
 * chisqr is in arcsec^2
 */

static double
wcs_chisqr (double *v, int iter)
{
    double chsq;
    double xmp, ymp, dx, dy, cd[4], *cdx;
    double crval1, crval2, cdelt1, cdelt2, crota, crpix1, crpix2;
    int i, offscale;

    /* Set WCS parameters from fit parameter vector */

    /* Sky coordinates at optical axis (degrees) */
    if (vfit[1] > -1)
        crval1 = xref_p + v[vfit[1]];
    else
        crval1 = wcsf->xref;

    if (vfit[2] > -1)
        crval2 = yref_p + v[vfit[2]];
    else
        crval2 = wcsf->yref;

    if (crval2 > 90.0 || crval2 < -90.0) {
        return 1000000000.0;
    }

    /* CD matrix */
    if (vfit[6] > -1) {
        cdelt1 = 0.0;
        cdelt2 = 0.0;
        crota = 0.0;
        cd[0] = v[vfit[3]];
        cd[1] = v[vfit[4]];
        cd[2] = v[vfit[5]];
        cd[3] = v[vfit[6]];
        cdx = cd;
    } else {
        /* Plate scale (degrees/pixel) */
        if (vfit[3] > -1)
            cdelt1 = v[vfit[3]];
        else
            cdelt1 = wcsf->xinc;

        if (vfit[4] > -1)
            cdelt2 = v[vfit[4]];
        else if (vfit[3] > -1) {
            if (cdelt1 < 0)
                cdelt2 = -cdelt1;
            else
                cdelt2 = cdelt1;
        } else
            cdelt2 = wcsf->yinc;

        /* Rotation angle (degrees) */
        if (vfit[5] > -1)
            crota = v[vfit[5]];
        else
            crota = wcsf->rot;

        cdx = NULL;
    }

    /* Optical axis pixel coordinates */
    if (vfit[7] > -1)
        crpix1 = xrefpix + v[vfit[7]];
    else
        crpix1 = wcsf->xrefpix;

    if (vfit[8] > -1)
        crpix2 = yrefpix + v[vfit[8]];
    else
        crpix2 = wcsf->yrefpix;

    if (wcsreset(wcsf, crpix1, crpix2, crval1, crval2, cdelt1, cdelt2, crota, cdx)) {
        fprintf (stderr,"CHISQR: Cannot reset WCS!\n");
#ifdef POLE_FIX
        return 1000000000.0;
#endif /* POLE_FIX */
    }

#ifdef POLE_FIX
    if (wcsf->offscl != 0) {
        return 1000000000.0;
    }
#endif /* POLE_FIX */

    /* Compute sum of squared residuals for these parameters */
    chsq = 0.0;

    for (i = 0; i < nbin_p; i++) {
        wcs2pix(wcsf, gra_p[i], gdec_p[i], &xmp, &ymp, &offscale);
        dx = xmp - sx_p[i];
        dy = ymp - sy_p[i];
        chsq += dx*dx + dy*dy;
    }

#if 0
    {
        char rastr[32], decstr[32];

        ra2str (rastr, 31, wcsf->xref, 3);
        dec2str (decstr, 31, wcsf->yref, 2);

        if (vfit[6] > -1)
            fprintf(
                stderr,
                "%4d: %s %s CD: %9.7f,%9.7f,%9.7f,%9.7f ",
                iter,
                rastr,
                decstr,
                wcsf->cd[0],
                wcsf->cd[1],
                wcsf->cd[2],
                wcsf->cd[3]
            );
        else
            fprintf(
                stderr,
                "%4d: %s %s %9.7f,%9.7f %8.5f ",
                iter,
                rastr,
                decstr,
                wcsf->xinc * 3600.0,
                wcsf->yinc * 3600.0,
                wcsf->rot
            );

        fprintf(stderr, "(%8.2f,%8.2f) -> %f\n", wcsf->xrefpix, wcsf->yrefpix, chsq);
    }
#endif

    return chsq;
}

/* The following subroutines are based on those in Numerical Recipes in C */

/* amoeba.c */

#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

static void
dasch_amoeba (
    double **p,
    double y[],
    int ndim,
    double ftol,
    int itmax,
    double (*funk)(double *, int),
    int *nfunk
) {
    int i, j, ilo, ihi, inhi, ndim1 = ndim + 1;
    double ytry, ysave, sum, rtol, *psum;

    psum = (double *) malloc ((unsigned)ndim * sizeof(double));
    *nfunk = 0;

    for (j = 0; j < ndim; j++) {
        for (i = 0, sum = 0.0; i < ndim1; i++)
            sum += p[i][j];

        psum[j] = sum;
    }

    for (;;) {
        ilo = 1;

        if (y[0] > y[1]) {
            inhi = 1;
            ihi = 0;
        } else {
            inhi = 0;
            ihi = 1;
        }

        for (i = 0; i < ndim1; i++) {
            if (y[i] < y[ilo])
                ilo = i;

            if (y[i] > y[ihi]) {
                inhi = ihi;
                ihi = i;
            } else if (y[i] > y[inhi])
                if (i != ihi)
                    inhi = i;
        }

        rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]));
        if (rtol < ftol)
            break;

        if (*nfunk >= itmax) {
            fprintf(stderr, "Too many iterations in amoeba fit %d > %d", *nfunk, itmax);
            return;
        }

        ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, -ALPHA);

        if (ytry <= y[ilo])
            ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, GAMMA);
        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, BETA);

            if (ytry >= ysave) {
                for (i = 0; i < ndim1; i++) {
                    if (i != ilo) {
                        for (j = 0; j < ndim; j++) {
                            psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                            p[i][j] = psum[j];
                        }

                        y[i] = (*funk)(psum, *nfunk);
                    }
                }

                *nfunk += ndim;

                for (j=0; j < ndim; j++) {
                    for (i = 0, sum = 0.0; i < ndim1; i++)
                        sum += p[i][j];

                    psum[j] = sum;
                }
            }
        }
    }

    free(psum);
}


static double
amotry (
    double **p,
    double *y,
    double *psum,
    int ndim,
    double (*funk)(double *, int),
    int ihi,
    int *nfunk,
    double fac
) {
    int j;
    double fac1, fac2, ytry, *ptry;

    ptry = (double *) malloc ((unsigned) ndim * sizeof(double));
    fac1 = (1.0 - fac) / ndim;
    fac2 = fac1 - fac;

    for (j = 0; j < ndim; j++)
        ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;

    ytry = (*funk)(ptry, *nfunk);
    ++(*nfunk);

    if (ytry < y[ihi]) {
        y[ihi] = ytry;

        for (j = 0; j < ndim; j++) {
            psum[j] +=  ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }

    free (ptry);
    return ytry;
}

void
dasch_setnfit (int nfit)
{
    if (nfit == 0)
        setnofit();
    else if (nfit < 0) {
        fprintf(stderr, "fatal error: negative nfit values not supported in dasch-imwcs\n");
        exit(1);
    } else {
        pfit0 = nfit;
    }
}

int
dasch_getnfit(void)
{
    return pfit0;
}

void
dasch_setminbin(int minbin1)
{
    minbin = minbin1;
}

#ifdef LIMIT_SHIFT
void
dasch_setmaxshift(int maxshift)
{
  maxshift0 = maxshift;
}
#endif /* LIMIT_SHIFT */
