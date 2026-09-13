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

#include <libwcs/wcs.h>
#include <libwcs/lwcs.h>
#include <libwcs/fitsfile.h>
#include <libwcs/fitswcs.h>
#include <libwcs/wcscat.h>

static double tolerance = PIXDIFF;	/* +/- this many pixels is a hit */
static double refmag1 = MAGLIM1;	/* reference catalog magnitude limit */
static double refmag2 = MAGLIM2;	/* reference catalog magnitude limit */
static double frac = 1.0;	/* Additional catalog/image stars */
static int nofit = 0;		/* if =1, do not fit WCS */
static int maxcat = MAXSTARS;	/* Maximum number of catalog stars to use */
static int fitwcs = 1;		/* If 1, fit WCS, else use current WCS */
static int fitplate = 0;	/* If 1, fit polynomial, else do not */
static int iterate0 = 0;	/* If 1, search field again */
static int toliterate0 = 0;	/* if 1, halve tolerances when iter */
static int nfiterate0 = 0;	/* if 1, add two parameters to fit */
static int recenter0 = 0;	/* If 1, search again with new center*/
static int irafout = 0;		/* if 1, write X Y RA Dec out */
static int magfit = 0;		/* If 1, write magnitude polynomial(s) */
static int sortmag = 1;		/* Magnitude by which to sort stars */
static int minstars0 = MINSTARS;	/* Number of star matches for fit */
static int nmagmax = MAXNMAG;	/* Maximum number of magnitudes (etc.) per entry */
static int nxydec = NXYDEC;	/* Number of decimal places in image coordinates */

static void PrintRes();
static void CompRes();

#ifdef CROP_BORDER
static int borderpix0 = 0;          /* Pixels to crop from the border */
static int squareborder0 = 0;       /* Make the cropping region square */
#endif /* CROP_BORDER */

static char *kwt = NULL;        /* Keyword returned by ctgread() */



int
dasch_SetWCSFITS (
    char *filename,
    char *header,
    char *image,
    char *refcatname,
    int verbose
) {
    double *gnum;	/* Reference star numbers */
    double *gra;	/* Reference star right ascensions in degrees */
    double *gdec;	/* Reference star declinations in degrees */
    double *gpra;	/* Reference star right ascension proper motions (deg)*/
    double *gpdec;	/* Reference star declination proper motions (deg) */
    double **gm;	/* Reference star magnitudes */
    double *gx;		/* Reference star image X-coordinates in pixels */
    double *gy;		/* Reference star image Y-coordinates in pixels */
    int *gc;		/* Reference object types */
    int *goff;		/* Reference star offscale flags */
    int ng;		/* Number of reference stars in image */
    int nbg;		/* Number of brightest reference stars from search */
    int nrg;		/* Number of brightest reference stars actually used */
    double *sx;		/* Image star image X-coordinates in pixels */
    double *sy;		/* Image star image X-coordinates in pixels */
    double *sm;		/* Image star instrumental magnitude */
    int *sp;		/* Image star peak fluxes in counts */
    int ns;		/* Number of image stars */
    int nbs;		/* Number of brightest image stars actually used */
    double cra, cdec;	/* Nominal center in degrees from RA/DEC FITS fields */
    double dra, ddec;	/* Image half-widths in degrees */
    double drac, ddecc; /* Cropped image half-width in degrees */
    double secpix;	/* Pixel size in arcseconds */
    int imw, imh;	/* Image size, pixels */
    int imsearch = 1;	/* Flag set if image should be searched for sources */
    int nmax;		/* Maximum number of matches possible (nrg or nbs) */
    double mag1,mag2;
    int refcat = 0;		/* reference catalog switch */
    int nmag, mprop;
    double dxys;
    int ngmax = 0;
    int nbin, nbytes;
    int iterate, toliterate, nfiterate;
    int imag, magsort;
    int niter = 0;
    int recenter = recenter0;
    int ret = 0;
    int is, ig, igs, i;
    char rstr[32], dstr[32];
    double refeq, refep;
    int refsys;
    char refcoor[8];
    char title[80];
    char *imcatname;	/* file name for image star catalog, if used */
    struct WorldCoor *wcs=0;	/* WCS structure */
    double *sx1, *sy1, *sm1, *gra1, *gdec1, *gnum1, *gm1;
    char **gobj, **gobj1;	/* Catalog star object names */
    int nmatch;
    double dx, dy, dx2, dy2, dxy;
    struct StarCat *starcat;

    iterate = iterate0;
    toliterate = toliterate0;
    nfiterate = nfiterate0;
    gnum = NULL;
    gra = NULL;
    gdec = NULL;
    gpra = NULL;
    gpdec = NULL;
    gm = NULL;
    gx = NULL;
    gy = NULL;
    gc = NULL;
    gobj = NULL;
    gobj1 = NULL;
    goff = NULL;
    sm = NULL;
    sx = NULL;
    sy = NULL;
    sp = NULL;
    starcat = NULL;
    imcatname = NULL;
    ns = 0;

    if (refmag1 == refmag2) {
        mag1 = 0.0;
        mag2 = 0.0;
    } else {
        mag1 = refmag1;
        mag2 = refmag2;
    }

    printf(
        "@* init iterate=%d toliterate=%d nfiterate=%d mag1=%.3f mag2=%.3f sortmag=%d\n",
        iterate, toliterate, nfiterate, mag1, mag2, sortmag
    );

    /* Set reference catalog coordinate system and epoch */
    if (nofit) {
        refsys = 0;
        refeq = 0.0;
        refcat = 0;
    } else {
        refcat = RefCat (refcatname,title,&refsys,&refeq,&refep,&mprop,&nmag);
        wcscstr (refcoor, refsys, refeq, refep);
    }

    /* get nominal position and scale */
getfield:
    wcs = GetFITSWCS(
        filename,
        header,
        0, // verbose
        &cra,
        &cdec,
        &dra,
        &ddec,
        &secpix,
        &imw,
        &imh,
        &refsys,
        &refeq
    );

    if (nowcs (wcs)) {
        ret = 0;
        goto out;
    }

    refep = wcs->epoch;

    printf(
        "@* getwcs cra=%.5f cdec=%.5f dra=%.5f ddec=%.5f secpix=%.5f imw=%d imh=%d refsys=%d refeq=%.2f refep=%.3f\n",
        cra, cdec, dra, ddec, secpix, imw, imh, refsys, refeq, refep
    );

    if (nofit) {
        SetFITSWCS (header, wcs);
        ret = 1;
        goto out;
    }

    if (fitwcs) {
        wcs->prjcode = WCS_TAN;
        wcseqset (wcs, refeq);
    }

    if (refcatname == NULL) {
        refcatname = CatName (refcat, refcatname);

        if (refcatname == NULL) {
            ret = 0;
            goto out;
        }
    }

    if (sortmag > 9)
        sortmag = CatMagNum (sortmag, refcat);

    // Allocate arrays for results of reference star search
    ngmax = maxcat;
    nbytes = ngmax * sizeof (double);

    if (!(gnum = (double *) calloc (ngmax, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gnum\n", nbytes);
    if (!(gra = (double *) calloc (ngmax, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gra\n", nbytes);
    if (!(gdec = (double *) calloc (ngmax, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gdec\n", nbytes);
    if (!(gpra = (double *) calloc (ngmax, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gpra\n", nbytes);
    if (!(gpdec = (double *) calloc (ngmax, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gpdec\n", nbytes);
    if (!(gm = (double **) calloc (nmagmax, sizeof(double *))))
        fprintf (stderr, "@! Could not calloc %zu bytes for gm\n", nmagmax*sizeof(double *));
    else {
        for (imag = 0; imag < nmagmax; imag++) {
            if (!(gm[imag] = (double *) calloc (ngmax, sizeof(double))))
                fprintf (stderr, "@! Could not calloc %d bytes for gm\n", nbytes);
        }
    }

    if (!(gc = (int *) calloc (ngmax, sizeof(int))))
        fprintf (stderr, "@! Could not calloc %zu bytes for gc\n", ngmax*sizeof(int));

    if (!(gobj = (char **) calloc (ngmax, sizeof(char *))))
        fprintf (stderr, "@! Could not calloc %zu bytes for obj\n", ngmax*sizeof(char *));
    else {
        for (i = 0; i < ngmax; i++)
            gobj[i] = NULL;
    }

    /* Find the nearby reference stars, in ra/dec */
getstars:
    printf(
        "\n@> SetWCSFITS iteration: iterate=%d toliterate=%d nfiterate=%d\n\n",
        iterate, toliterate, nfiterate
    );

    drac = dra;
    ddecc = ddec;

#ifdef CROP_BORDER
    {
        double cropdegrees;

        cropdegrees = getsecpix() * borderpix0 / 3600.0;
        if (cropdegrees < 0.0) {
            cropdegrees = -cropdegrees;
        }

        drac = dra - cropdegrees;
        ddecc = ddec - cropdegrees;

        if (squareborder0 && cdec < 89.0) {
            double temprac = drac * cos(degrad(cdec));

            if (ddecc < temprac) {
                drac = ddec/cos(degrad(cdec));
            } else {
                ddecc = temprac;
            }
        }
    }
#endif /* CROP_BORDER */

    ng = ctgread(
        refcatname, // refcat name
        refcat, // refcat ID
        0, // sort by distance from center?
        cra, // search center (degrees)
        cdec,
        drac, // search half-width (degrees)
        ddecc,
        0.0, // no limiting separation
        0.0, // no inner annulus size
        refsys, // search coordinate system
        refeq, // search coordinate equinox
        refep, // search proper motion epoch
        mag1, // limiting magnitude 1
        mag2, // limiting magnitude 2
        sortmag, // index of magnitude to limit and sort by
        ngmax, // maximum number of entries to return
        &starcat, // the catalog data structure
        gnum, // array of refcat ID numbers (filled in)
        gra, // array of coordinates (filled in)
        gdec,
        gpra, // array of proper motion values (filled in)
        gpdec,
        gm, // 2D array of magnitudes (filled in)
        gc, // array of fluxes (filled in)
        gobj, // array of object names (filled in)
        0 // "nlog"; verbosity parameter
    );

    printf(
        "@* ctgread ng=%d ngmax=%d drac=%.6f ddecc=%.6f\n",
        ng, ngmax, drac, ddecc
    );

    // ngmax: maximum number of of entries to return
    // ng: number of entries found in the catalog; may be larger or smaller than ngmax!
    // nrg = min(ng, ngmax): number of entries in our arrays

    if (ng > ngmax)
        nrg = ngmax;
    else
        nrg = ng;

    if (gobj[0] == NULL)
        gobj1 = NULL;
    else
        gobj1 = gobj;

    if (sortmag > 0 && sortmag <= nmag)
        magsort = sortmag - 1;
    else
        magsort = 0;

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (
        gnum,
        gra,
        gdec,
        gpra,
        gpdec,
        NULL,
        NULL,
        gm,
        gc,
        gobj1,
        nrg,
        nmagmax,
        sortmag
    );

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    if (!(gx = (double *) calloc (ngmax, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %zu bytes for gx\n", ngmax*sizeof(double));

    if (!(gy = (double *) calloc (ngmax, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %zu bytes for gy\n", ngmax*sizeof(double));

    if (!(goff = (int *) calloc (ngmax, sizeof(int))))
        fprintf (stderr, "@! Could not calloc %zu bytes for gy\n", ngmax*sizeof(double));

    if (!gx || !gy || !goff) {
        ret = 0;
        goto out;
    }

    /* use the nominal WCS info to find x/y on image */
    for (ig = 0; ig < nrg; ig++) {
        gx[ig] = 0.0;
        gy[ig] = 0.0;
        wcs2pix (wcs, gra[ig], gdec[ig], &gx[ig], &gy[ig], &goff[ig]);
    }

    // 2024 DASCH customization: filter out catalog sources off image. For
    // DASCH's large images, the catalog query can yield many sources that are
    // far off the image, which significantly decreases the number of actually
    // useful catalog sources due to the `ngmax` limitation that is applied.

    {
        int i, nrg_orig = nrg, next_rec = 0;

        for (i = 0; i < nrg; i++) {
            if (gx[i] >= 0 && gx[i] < imw && gy[i] >= 0 && gy[i] < imh) {
                if (next_rec != i) {
                    gnum[next_rec] = gnum[i];
                    gra[next_rec] = gra[i];
                    gdec[next_rec] = gdec[i];
                    gm[magsort][next_rec] = gm[magsort][i];
                    gpra[next_rec] = gpra[i];
                    gpdec[next_rec] = gpdec[i];
                    gx[next_rec] = gx[i];
                    gy[next_rec] = gy[i];
                    goff[next_rec] = goff[i];
                }

                next_rec++;
            }
        }

        nrg = next_rec;
        fprintf(stderr, "Spatial filter: %d => %d\n", nrg_orig, nrg);
    }

    // Log outcome of reference star selection

    if (verbose) {
        if (ng > nrg)
            printf("Using %d / %d reference stars brighter than %.1f\n", nrg, ng, gm[magsort][nrg-1]);
        else if (refmag1 > 0.0 && refmag2 > 0.0)
            printf("Using all %d reference stars from %.1f to %.1f mag.\n", ng, refmag1, refmag2);
        else if (refmag2 > 0.0)
            printf("Using all %d reference stars brighter than %.1f\n", ng, refmag2);
        else
            printf("Using all %d reference stars\n", ng);
    }

    printf(
        "@* maglim refmag1=%.3f refmag2=%.3f maglim=%.3f minstars0=%d\n",
        refmag1, refmag2, gm[magsort][nrg - 1], minstars0
    );

    // Report the catalog stars as CSV
    if (verbose) {
#if 0
        // XXX hardcoding nmag > 0
        printf("@> catalog %s stars (%d rows):\n\nx,y,r,d,m\n", refcatname, nrg);

        for (ig = 0; ig < nrg; ig++) {
            printf("%.1f,%.1f,%.6f,%.6f,%.2f\n", gx[ig], gy[ig], gra[ig], gdec[ig], gm[magsort][ig]);
        }

        printf("\n\n");
#endif
    }

    if (nrg < minstars0) {
        if (ng < 0)
            fprintf(stderr, "@! Error getting reference stars: %d\n", ng);
        else if (ng == 0)
            fprintf(stderr, "@! No reference stars found in image area\n");
        else if (fitwcs)
            fprintf(stderr, "@! Found only %d out of %d reference stars needed\n", nrg, minstars0);

        if (ng <= 0 || fitwcs) {
            ret = 0;
            goto out;
        }
    }

    // Discover star-like things in the image, in pixels. If setimcat() has been
    // called (which it is, in DASCH), this will read in an existing table file
    // rather than doing its own source detection.
    if (imsearch) {
        ns = FindStars(
            header,
            image,
            &sx,
            &sy,
            &sm,
            &sp,
            0, // verbose
            0 // zap
        );
        printf("@* findstars ns=%d\n", ns);

#ifdef CROP_BORDER
        if (borderpix0 > 0) {
            double minborderx = borderpix0;
            double minbordery = borderpix0;
            double maxborderx = imw - borderpix0;
            double maxbordery = imh - borderpix0;
            int borderindex = 0;

            if (squareborder0) {
                int deltaborderx = maxborderx - minborderx;
                int deltabordery = maxbordery - minbordery;

                if (deltaborderx < deltabordery) {
                    minbordery = (imh/2) - (deltaborderx/2);
                    maxbordery = (imh/2) + (deltaborderx/2);
                } else {
                    minborderx = (imw/2) - (deltabordery/2);
                    maxborderx = (imw/2) + (deltabordery/2);
                }
            }

            printf("Starting ns=%d borderpix=%d;", ns, borderpix0);

            while (borderindex < ns) {
                if ((sx[borderindex] < minborderx) ||
                    (sx[borderindex] > maxborderx) ||
                    (sy[borderindex] < minbordery) ||
                    (sy[borderindex] > maxbordery))
                {
                    /* Star is outside border, skip it */
                    sx[borderindex] = sx[ns-1];
                    sy[borderindex] = sy[ns-1];
                    sm[borderindex] = sm[ns-1];
                    sp[borderindex] = sp[ns-1];
                    ns--;
                } else {
                    borderindex++;
                }
            }

            printf(" ending ns=%d\n", ns);
        }
#endif /* CROP_BORDER */

        if (ns < minstars0) {
            if (ns < 0)
                fprintf(stderr, "@! Error getting image stars: %d\n", ns);
            else if (ns == 0)
                fprintf(stderr, "@! No stars found in image\n");
            else if (fitwcs)
                fprintf(stderr, "@! Need at least %d image stars but only found %d\n", minstars0, ns);

            if (ns <= 0 || fitwcs) {
                ret = 0;
                iterate = 0;
                recenter = 0;
                goto out;
            }
        }

        imsearch = 0;
    }

    /* Fit a world coordinate system if requested */
    if (fitwcs) {
        niter++;

        /* Sort star-like objects in image by brightness (magnitude) */
        MagSortStars(
            NULL,
            NULL,
            NULL,
            NULL,
            NULL,
            sx,
            sy,
            &sm,
            sp,
            NULL,
            ns,
            1,
            1
        );

        // If matching a catalog field the same size as the image field, use
        // only as many star-like objects as reference stars, scaling by `frac`
        // if specified.

        if (ns > nrg) {
            nbs = nrg * frac;
            if (nbs > ns)
                nbs = ns;
            nbg = nrg;
        } else {
            nbs = ns;
            nbg = nbs * frac;
            if (nbg > nrg)
                nbg = nrg;
        }

        printf("@* smagsort nbg=%d nbs=%d ng=%d frac=%.2lf\n", nbg, nbs, ng, frac);

        if (verbose) {
            if (nbg == ng)
                printf("Using all %d reference stars\n", ng);
            else
                printf("Using brightest %d / %d reference stars\n", nbg, ng);

            if (nbs == ns)
                printf("Using all %d image stars\n", ns);
            else
                printf("Using brightest %d / %d image stars\n", nbs, ns);
        }

        // Report image stars as CSV
#if 0
        if (verbose) {
            printf("@> image stars (%d rows):\n\nx,y,r,d,m\n", nbs);

            for (is = 0; is < nbs; is++) {
                double ra, dec;

                pix2wcs (wcs, sx[is], sy[is], &ra, &dec);
                printf("%.1f,%.1f,%.6f,%.6f,%.2f\n", sx[is], sy[is], ra, dec, sm[is]);
            }

            printf("\n\n");
        }
#endif

        // Match offsets between all pairs of image stars and reference stars
        // and fit WCS to matches

        nbin = StarMatch(
            nbs,
            sx,
            sy,
            refcat,
            nbg,
            gnum,
            gra,
            gdec,
            goff,
            gx,
            gy,
            tolerance,
            wcs,
            0 // verbose
        );

        if (nbin < 0) {
            fprintf(stderr, "@! Star registration failed.\n");
            ret = 0;
            goto done;
        } else if (nbin < minstars0) {
            fprintf(stderr, "@! Only %d matches, registration failed.\n", nbin);
            ret = 0;
            goto done;
        } else if (verbose)
            printf("%d / %d bin hits\n", nbin, nbg);

        printf(
            "@* starmatch nbin=%d minstars0=%d tolerance=%.3f\n",
            nbin, minstars0, tolerance
        );

        hputs (header, "WCSRFCAT", refcatname);

        imcatname = getimcat ();
        if (strlen (imcatname) == 0)
            hputs (header, "WCSIMCAT", filename);
        else
            hputs (header, "WCSIMCAT", imcatname);

        hputi4 (header, "WCSMATCH", nbin);

        if (ns < nbg)
            hputi4 (header, "WCSNREF", ns);
        else
            hputi4 (header, "WCSNREF", nbg);

        hputnr8 (header, "WCSTOL", 4, tolerance);
        SetFITSWCS (header, wcs);
    }

    /* Match reference and image stars */
    nmatch = 0;

    if (verbose || !fitwcs) {
        imcatname = getimcat();
        printf("Current WCS:\n");

        if (wcs->ncoeff1 > 0)
            printf ("  %d-term x, %d-term y polynomial fit\n", wcs->ncoeff1, wcs->ncoeff2);
        else
            printf ("  Arcsec/Pixel=%.6f %.6f  Rotation=%.6f degrees\n", 3600.0*wcs->xinc, 3600.0*wcs->yinc, wcs->rot);

        ra2str (rstr, 32, wcs->xref, 3);
        dec2str (dstr, 32, wcs->yref, 2);
        printf("  Optical axis=%s  %s %s x=%.2f y=%.2f\n", rstr, dstr, refcoor, wcs->xrefpix, wcs->yrefpix);
        printf(
            "@* wcs xref=%.6f yref=%.6f xinc=%.6f yinc=%.6f xrefpix=%.2f yrefpix=%.2f rot=%.3f\n",
            wcs->xref, wcs->yref, wcs->xinc, wcs->yinc, wcs->xrefpix, wcs->yrefpix, wcs->rot
        );
    }

    /* Find star matches for this offset and print them */

    /* Use the fit WCS info to find catalog star x/y on image */
    for (ig = 0; ig < nrg; ig++) {
        gx[ig] = 0.0;
        gy[ig] = 0.0;
        wcs2pix (wcs, gra[ig], gdec[ig], &gx[ig], &gy[ig], &goff[ig]);
    }

    /* Set maximum number of matches which are possible */
    if (nrg < ns)
        nmax = nrg;
    else
        nmax = ns;

    /* Find best catalog matches to stars in image */
    nmatch = 0;
    nbytes = ns * sizeof (double);

    if (!(gra1 = (double *) calloc (ns, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gra1\n", nbytes);
    if (!(gdec1 = (double *) calloc (ns, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gdec1\n", nbytes);
    if (!(gm1 = (double *) calloc (ns, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gm1\n", nbytes);
    if (!(gnum1 = (double *) calloc (ns, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for gnum1\n", nbytes);
    if (!(sx1 = (double *) calloc (ns, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for sx1\n", nbytes);
    if (!(sy1 = (double *) calloc (ns, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for sy1\n", nbytes);
    if (!(sm1 = (double *) calloc (ns, sizeof(double))))
        fprintf (stderr, "@! Could not calloc %d bytes for sm1\n", nbytes);

    for (is = 0; is < ns; is++) {
        // Find the smallest separation, if it is smaller than tol^2
        dxys = tolerance * tolerance;
        igs = -1;

        for (ig = 0; ig < nrg; ig++) {
            if (!goff[ig]) {
                dx = gx[ig] - sx[is];
                dy = gy[ig] - sy[is];
                dx2 = dx * dx;
                dy2 = dy * dy;
                dxy = dx2 + dy2;

                if (dxy < dxys) {
                    dxys = dxy;
                    igs = ig;
                }
            }
        }

        // Did we find a match?
        if (igs > -1) {
            gnum1[nmatch] = gnum[igs];
            if (gm != NULL && nmag > 0)
                gm1[nmatch] = gm[magsort][igs];
            else
                gm1[nmatch] = 0.0;

            gra1[nmatch] = gra[igs];
            gdec1[nmatch] = gdec[igs];
            sx1[nmatch] = sx[is];
            sy1[nmatch] = sy[is];
            sm1[nmatch] = sm[is];
            nmatch++;
        }
    }

    printf("@* match nmatch=%d nmax=%d tolerance=%.3f\n", nmatch, nmax, tolerance);

    /* If there were any matches found, print them */
    if (nmatch > 0) {
        int rprint = verbose || !fitwcs;

        hputi4 (header, "WCSMATCH", nmatch);
        hputi4 (header, "WCSNREF", nmax);
        hputnr8 (header, "WCSTOL", 4, tolerance);

        if (rprint) {
            PrintRes(
                header,
                wcs,
                nmatch,
                sx1,
                sy1,
                sm1,
                gra1,
                gdec1,
                gm1,
                gnum1,
                refcat,
                rprint
            );

            if (strlen(imcatname) == 0)
                printf("  nmatch=%d nstars=%d between %s and %s  niter=%d\n", nmatch, nmax, refcatname, filename, niter);
            else
                printf("  nmatch=%d nstars=%d between %s and %s  niter=%d\n", nmatch, nmax, refcatname, imcatname, niter);
        } else
            CompRes(header, wcs, nmatch, sx1, sy1, sm1, gra1, gdec1, gm1, gnum1);

        /* Fit the matched catalog and image stars with a polynomial */
        if (!iterate && !recenter && fitplate && refcatname != NULL) {
            if (verbose)
                printf("Fitting matched stars with a polynomial\n");

            /* Fit residuals */
            if (FitPlate (wcs, sx1, sy1, gra1, gdec1, nmatch, fitplate, verbose))
                fprintf(stderr, "FitPlate cannot fit matches\n");
            else if (rprint) {
                PrintRes (header,wcs,nmatch,sx1,sy1,sm1,gra1,gdec1,gm1,gnum1, refcat,verbose);

                if (strlen (imcatname) == 0)
                    printf("# nmatch=%d nstars=%d between %s and %s niter=%d\n", nmatch, nmax, refcatname, filename, niter);
                else
                    printf("# nmatch=%d nstars=%d between %s and %s niter=%d\n", nmatch, nmax, refcatname, imcatname, niter);

                SetFITSPlate(header, wcs);
            } else {
                CompRes (header,wcs,nmatch,sx1,sy1,sm1,gra1,gdec1,gm1,gnum1);
                SetFITSPlate (header, wcs);
            }
        }
    } else {
        if (strlen (imcatname) == 0)
            fprintf(stderr, "@! SetWCSFITS: No matches between %s and %s\n", refcatname, filename);
        else
            fprintf(stderr, "@! SetWCSFITS: No matches between %s and %s\n", refcatname, imcatname);

        hputi4 (header, "WCSMATCH", 0);
    }

    free ((char *)gra1);
    free ((char *)gdec1);
    free ((char *)gm1);
    free ((char *)gnum1);
    free ((char *)sx1);
    free ((char *)sy1);
    free ((char *)sm1);
    ret = 1;

out:

    if (iterate) {
        wcssize (wcs, &cra, &cdec, &dra, &ddec);

        if (cra < 0.0)
            cra = cra + 360.0;

        iterate--;
        goto getstars;
    }

    if (toliterate) {
        wcssize (wcs, &cra, &cdec, &dra, &ddec);
        tolerance = tolerance * 0.5;
        toliterate--;
        goto getstars;
    }

    if (recenter) {
        double ra, dec, x, y;

        x = 0.5 * wcs->nxpix;
        y = 0.5 * wcs->nypix;
        pix2wcs (wcs, x, y, &ra, &dec);
        setdcenter (ra, dec);
        setsys (wcs->syswcs);
        setrefpix (x, y);
        setsecpix (-3600.0 * wcs->xinc);
        setsecpix2 (3600.0 * wcs->yinc);
        setrot (wcs->rot);
        recenter = 0;

        if (wcs) {
            wcsfree (wcs);
            wcs = NULL;
        }

        goto getfield;
    }

    if (nfiterate) {
        int nfit = getnfit();

        wcssize (wcs, &cra, &cdec, &dra, &ddec);
        if (verbose)
            printf ("\n fitting %d instead of %d parameters\n", nfit+2, nfit);

        if (nfit < 7)
            nfit = nfit + 2;

        setnfit (nfit);
        nfiterate--;
        goto getstars;
    }

done:
    if (wcs) {
        wcsfree (wcs);
        wcs = NULL;
    }

    /* Free catalog source arrays */
    free ((char *)gra);
    free ((char *)gdec);
    free ((char *)gpra);
    free ((char *)gpdec);

    if (gm) {
        for (imag = 0; imag < nmagmax; imag++) {
            free ((char *)gm[imag]);
        }

        free ((char *)gm);
    }

    free ((char *)gnum);
    free ((char *)gx);
    free ((char *)gy);
    free ((char *)goff);
    free ((char *)gc);

    /* Free memory used for object names in reference catalog */
    if (gobj1 != NULL) {
        for (i = 0; i < ngmax; i++) {
            if (gobj[i] != NULL) {
                free (gobj[i]);
                gobj[i] = NULL;
            }
        }
    }

    if (gobj) {
        free ((char *) gobj);
        gobj = NULL;
    }

    /* Free image source arrays */
    free ((char *)sx);
    free ((char *)sy);
    free ((char *)sm);
    free ((char *)sp);
    printf("@* finish ret=%d\n", ret);
    return ret;
}


static void
PrintRes(
    char *header, // FITS header
    struct WorldCoor *wcs, // WCS
    int nmatch, // number of matches
    double *sx1, // image star pixel coordinates
    double *sy1,
    double *sm1, // image magnitudes
    double *gra1, // reference catalog sky coordinates
    double *gdec1,
    double *gm1, // reference catalog magnitudes
    double *gnum1, // reference catalog numbers
    int refcat, // refcat code
    int verbose
) {
    int i, goff;
    double gx, gy, dx, dy, dx2, dy2, dxy, mag0;
    double sep, sep2, rsep, rsep2, dsep, dsep2;
    double dmatch, dmatch1, sra, sdec;
    double sepsum = 0.0;
    double rsepsum = 0.0;
    double rsep2sum = 0.0;
    double dsepsum = 0.0;
    double dsep2sum = 0.0;
    double sep2sum = 0.0;
    double dxsum = 0.0;
    double dysum = 0.0;
    double dx2sum = 0.0;
    double dy2sum = 0.0;
    double dxysum = 0.0;
    double coeff[5];
    double msig;
    double maxnum;
    double cmax;
    int nnfld;
    int nxyfld;
    char rstr[32], dstr[32], numstr[32], xstr[32], ystr[32], mstr[8];

    maxnum = 0.0;
    for (i = 0; i < nmatch; i++) {
        if (i == 0)
            maxnum = gnum1[i];
        else if (gnum1[i] > maxnum)
            maxnum = gnum1[i];
    }

    nnfld = CatNumLen (refcat, maxnum, 0);
    CatMagName (sortmag, refcat, mstr);
    CatID (numstr, refcat);

#if 0
    if (irafout)
        printf ("#   x      y        ra2000   dec2000  %5s %s", mstr, numstr);
    else
        printf ("# %s ra2000       dec2000    %5s    X      Y     magi", mstr, numstr);

    printf ("    dra   ddec   sep\n");
#endif

    /* Find maximum image coordinates and set field size accordingly */
    cmax = 0.0;

    for (i = 0; i < nmatch; i++) {
        if (sx1[i] > cmax)
            cmax = sx1[i];
        if (sy1[i] > cmax)
            cmax = sy1[i];
    }

    if (cmax > 9999.0)
        nxyfld = 6 + nxydec;
    else if (cmax > 999.0)
        nxyfld = 5 + nxydec;
    else
        nxyfld = 4 + nxydec;

    for (i = 0; i < nmatch; i++) {
        wcs2pix (wcs, gra1[i], gdec1[i], &gx, &gy, &goff);

        dx = gx - sx1[i];
        dy = gy - sy1[i];
        dx2 = dx * dx;
        dy2 = dy * dy;
        dxy = dx2 + dy2;
        dxsum = dxsum + dx;
        dysum = dysum + dy;
        dx2sum = dx2sum + dx2;
        dy2sum = dy2sum + dy2;
        dxysum = dxysum + sqrt (dxy);
        pix2wcs (wcs, sx1[i], sy1[i], &sra, &sdec);

        sep = 3600.0 * wcsdist(gra1[i],gdec1[i],sra,sdec);
        rsep = 3600.0 * ((gra1[i]-sra) * cos(degrad(sdec)));

        if (rsep > sep)
            rsep = 3600.0 * ((gra1[i] - sra - 360.0) * cos(degrad(sdec)));

        rsep2 = rsep * rsep;
        dsep = 3600.0 * (gdec1[i] - sdec);
        dsep2 = dsep * dsep;
        sepsum = sepsum + sep;
        rsepsum = rsepsum + rsep;
        dsepsum = dsepsum + dsep;
        rsep2sum = rsep2sum + rsep2;
        dsep2sum = dsep2sum + dsep2;
        sep2sum = sep2sum + (sep*sep);
        ra2str (rstr, 32, gra1[i], 3);
        dec2str (dstr, 32, gdec1[i], 2);
        num2str (xstr, sx1[i], nxyfld, nxydec);
        num2str (ystr, sy1[i], nxyfld, nxydec);
        CatNum (refcat, -nnfld, 0, gnum1[i], numstr);

#if 0
        if (irafout)
            printf (" %s %s %s %s %5.2f %s", xstr, ystr, rstr, dstr, gm1[i], numstr);
        else
            printf ("%s %s %s %5.2f %s %s %6.2f ", numstr, rstr, dstr, gm1[i], xstr, ystr, sm1[i]);

        printf ("%6.2f %6.2f %6.2f\n", rsep, dsep, sep);
#endif
    }

    dmatch = (double) nmatch;
    dmatch1 = (double) (nmatch - 1);
    dx = dxsum / dmatch;
    dy = dysum / dmatch;
    dx2 = sqrt (dx2sum / dmatch1);
    dy2 = sqrt (dy2sum / dmatch1);
    dxy = dxysum / dmatch;
    rsep = rsepsum / dmatch;
    dsep = dsepsum / dmatch;
    rsep2 = sqrt (rsep2sum / dmatch1);
    dsep2 = sqrt (dsep2sum / dmatch1);
    sep = sepsum / dmatch;
    sep2 = sqrt (sep2sum / dmatch1);

    printf("Residuals:\n");
    printf("  Mean  dx=%.4f/%.4f  dy=%.4f/%.4f  dxy=%.4f\n", dx, dx2, dy, dy2, dxy);
    printf("  Mean dra=%.4f/%.4f  ddec=%.4f/%.4f sep=%.4f/%.4f\n", rsep, rsep2, dsep, dsep2, sep, sep2);

    printf(
        "@* mresid dx=%.4f dx2=%.4f dy=%.4f dy2=%.4f dxy=%.4f rsep=%.4f rsep2=%.4f dsep=%.4f dsep2=%.4f sep=%.4f sep2=%.4f\n",
        dx, dx2, dy, dy2, dxy, rsep, rsep2, dsep, dsep2, sep, sep2
    );

    /* Fit and save image to catalog magnitude calibration polynomial */
    if (magfit) {
        mag0 = sm1[0];
        coeff[0] = 0.0;
        coeff[1] = 0.0;
        coeff[2] = 0.0;
        coeff[3] = 0.0;
        coeff[4] = 0.0;
        polfit (sm1, gm1, mag0, nmatch, 4, coeff, &msig);
        printf ("# Plate to catalog mag: mag0=%.6f mcoeff0=%.6f mcoeff1=%.6f\n", mag0, coeff[0], coeff[1]);
        printf ("# Plate to catalog mag: mcoeff2=%.6f mcoeff3=%.6f sigma=%.3f\n", coeff[2], coeff[3], msig);
    }

    hputi4 (header, "WCSMATCH", nmatch);
    hputnr8 (header, "WCSSEP", 3, sep);
}


static void
CompRes (
    char *header, // image FITS header
    struct WorldCoor *wcs, // image WCS
    int nmatch, // number of matches
    double *sx1, // image star pixel coordinates
    double *sy1,
    double *sm1, // image magnitudes
    double *gra1, // reference catalog sky coordinates
    double *gdec1,
    double *gm1, // reference catalog magnitudes
    double *gnum1 // reference catalog numbers
) {
    int i, goff;
    double gx, gy, dx, dy, dx2, dy2, dxy;
    double sep, rsep, rsep2, dsep, dsep2;
    double dmatch, dmatch1, sra, sdec;
    double sepsum = 0.0;
    double rsepsum = 0.0;
    double rsep2sum = 0.0;
    double dsepsum = 0.0;
    double dsep2sum = 0.0;
    double sep2sum = 0.0;
    double dxsum = 0.0;
    double dysum = 0.0;
    double dx2sum = 0.0;
    double dy2sum = 0.0;
    double dxysum = 0.0;

    for (i = 0; i < nmatch; i++) {
        wcs2pix (wcs, gra1[i], gdec1[i], &gx, &gy, &goff);
        dx = gx - sx1[i];
        dy = gy - sy1[i];
        dx2 = dx * dx;
        dy2 = dy * dy;
        dxy = dx2 + dy2;
        dxsum = dxsum + dx;
        dysum = dysum + dy;
        dx2sum = dx2sum + dx2;
        dy2sum = dy2sum + dy2;
        dxysum = dxysum + sqrt (dxy);

        pix2wcs (wcs, sx1[i], sy1[i], &sra, &sdec);
        sep = 3600.0 * wcsdist(gra1[i],gdec1[i],sra,sdec);
        rsep = 3600.0 * ((gra1[i]-sra) * cos(degrad(sdec)));
        if (rsep > sep)
            rsep = 3600.0 * ((gra1[i] - sra - 360.0) * cos(degrad(sdec)));

        rsep2 = rsep * rsep;
        dsep = 3600.0 * (gdec1[i] - sdec);
        dsep2 = dsep * dsep;
        sepsum = sepsum + sep;
        rsepsum = rsepsum + rsep;
        dsepsum = dsepsum + dsep;
        rsep2sum = rsep2sum + rsep2;
        dsep2sum = dsep2sum + dsep2;
        sep2sum = sep2sum + (sep*sep);
    }

    dmatch = (double) nmatch;
    dmatch1 = (double) (nmatch - 1);
    dx = dxsum / dmatch;
    dy = dysum / dmatch;
    dx2 = sqrt (dx2sum / dmatch1);
    dy2 = sqrt (dy2sum / dmatch1);
    dxy = dxysum / dmatch;
    rsep = rsepsum / dmatch;
    dsep = dsepsum / dmatch;
    rsep2 = sqrt (rsep2sum / dmatch1);
    dsep2 = sqrt (dsep2sum / dmatch1);
    sep = sepsum / dmatch;

    hputi4 (header, "WCSMATCH", nmatch);
    hputnr8 (header, "WCSSEP", 3, sep);
}

/* Subroutines to initialize various parameters */

void
settolerance (tol)
double tol;
{ tolerance = tol; return; }


/* Number of decimal places in X and Y image coordinates of sources */
void
setnxydec (ndec)
int ndec;
{ nxydec = ndec; return; }

void
setirafout ()
{ irafout = 1; return; }

void
setreflim (lim1, lim2)
double lim1, lim2;
{ refmag2 = lim2;
  if (lim1 > -2.0) refmag1 = lim1;
  return; }

void
setfitwcs (wfit)
int wfit;
{ fitwcs = wfit; return; }

void
setfitplate (nc)
int nc;
{ fitplate = nc; return; }

void
setminstars (minstars)
int minstars;
{ minstars0 = minstars;
  setminbin (minstars);
  return; }

void
setnofit ()
{ nofit = 1; return; }

void
setfrac (frac0)
double frac0;
{ if (frac0 < 1.0) frac = 1.0 + frac0;
    else frac = frac0;
  return; }

void
setmaxcat (ncat)
int ncat;
{ if (ncat < 1) maxcat = 25;
  else maxcat = ncat;
  return; }

void
setiterate (iter)
int iter;
{ iterate0 = iterate0 + iter;
  return; }

void
setnfiterate (iter)
int iter;
{ nfiterate0 = nfiterate0 + iter;
  return; }

void
setiteratet (iter)
int iter;
{ toliterate0 = toliterate0 + iter;
  return; }

void
setrecenter (recenter)
int recenter;
{ recenter0 = recenter;
  return; }

void
setsortmag (imag)
int imag;
{ sortmag = imag;
  return; }

void
setmagfit ()
{magfit++; return;}

void settabkrw (keyword0)
char *keyword0;
{ kwt = keyword0; return; }

#ifdef CROP_BORDER
void
setcropborder(int borderpix) {
  borderpix0 = borderpix;
  return;
}
void
setsquareborder(int squareborder) {
  squareborder0 = squareborder;
  return;
}
#endif /* CROP_BORDER */
