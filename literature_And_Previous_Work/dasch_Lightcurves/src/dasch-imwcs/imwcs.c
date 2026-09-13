// This is a customized version of WCSTools imwcs, derived from:
//
//   IMWCS WCSTools 3.9.7, 11 August 2022, Jessica Mink (jmink@cfa.harvard.edu)
//
// It inherits the license of the above: GPLv2 or later.
//
// This version is customized to fix various convergence failures observed
// within the DASCH dataset.
//
// Arguments when called inside the DASCH astrometry pipeline. Numbers in brackets
// are changed between bin16 and bin01 passes.
//
//   dasch-imwcs
//     -c tycho2      -- reference catalog
//     -q it2         -- iterate; reduce tolerance by factor of 2 iteratively
//     -v             -- verbose
//     -w             -- write header
//     -h 2500 [25000] - maximum number of refstars to use
//     -n 8           -- number of parameters to fit
//     -t 20 [200]    -- offset tolerance in pixels
//     -B 15 [240]    -- CROP_BORDER pixels to crop from the border
//     -C 250         -- LIMIT_SHIFT pixel limit of correlation algorithm range
//     -o (output)    -- output file
//     -d (.sex file) -- get sources from this file instead of image search
//     (input)

#include "dasch_imwcs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>
#include <libwcs/wcscat.h>
#include <libwcs/lwcs.h>
#include <libwcs/fitswcs.h>

#define MAXFILES 1000
static int maxnfile = MAXFILES;

static void PrintUsage();
static void FitWCS();

static char *RevMsg = "DASCH imwcs";

static int verbose = 0;		/* verbose/debugging flag */
static int writeheader = 0;	/* write header fields; else read-only */
static int overwrite = 0;	/* allow overwriting of input image file */
static int rot = 0;		/* Angle to rotate image (multiple of 90 deg) */
static int mirror = 0;		/* If 1, flip image right-left before rotating*/
static int bitpix = 0;
static int fitsout = 0;		/* Output FITS file from IRAF input if 1 */
static int imsearch = 1;	/* set to 0 if image catalog provided */
static int erasewcs = 0;	/* Set to 1 to erase initial image WCS */
static int rotatewcs = 1;	/* If 1, rotate FITS WCS keywords in image */
static char outname[128];		/* Name for output image */
static char *refcatname;	/* Name of reference catalog to match */
static int version = 0;		/* If 1, print only program name and version */
static char *matchfile;		/* File of X Y RA Dec matches for initial fit */
static char *progname;		/* Name of program as executed */

int
main (int ac, char **av)
{
    char *str, *str1, c, c1, c2;
    double bmin, maglim1, maglim2, drot, arot;
    char rastr[32];
    char decstr[32];
    char errmsg[256];
    char **fn;
    double x, y;
    int i, imag;
    int ifile, nfile;

    outname[0] = 0;
    refcatname = NULL;
    matchfile = NULL;
    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    setrevmsg (RevMsg);

    // Some custom defaults
    dasch_setfrac(4.0);

    /* Check name used to execute programe and set catalog name accordingly */
    progname = ProgName (av[0]);
    refcatname = ProgCat (progname);

    if (ac == 1)
        PrintUsage (NULL);

    /* Loop through the arguments */
    for (av++; --ac > 0; av++) {
        str = *av;

        if (!strcmp (str, "--help"))
            PrintUsage (NULL);
        else if (!strcmp (str, "--version")) {
            version = 1;
            PrintUsage ("version");
        } else if (strchr(str, '='))
            setparm (str);
        else if (str[0] == '-') {
            while ((c = *++str) != '\0') {
                switch (c) {

                case 'a':       /* Initial rotation angle in degrees */
                    if (ac < 2)
                        PrintUsage (str);
                    drot = atof (*++av);
                    arot = fabs (drot);
                    if (arot != 90.0 && arot != 180.0 && arot != 270.0) {
                        setrot (drot);
                        rot = 0;
                        }
                    else
                        rot = atoi (*av);
                    ac--;
                    break;

                case 'b':   /* initial coordinates on command line in B1950 */
                    if (ac < 3)
                        PrintUsage (str);
                    setsys (WCS_B1950);
                    strcpy (rastr, *++av);
                    ac--;
                    strcpy (decstr, *++av);
                    ac--;
                    setcenter (rastr, decstr);
                    break;

                case 'c':       /* Set reference catalog */
                    if (ac < 2)
                        PrintUsage (str);
                    refcatname = *++av;
                    ac--;
                    break;

                case 'd':  /* Read image star positions from DAOFIND file */
                    if (ac < 2)
                        PrintUsage (str);
                    setimcat (*++av);
                    imsearch = 0;
                    ac--;
                    break;

                case 'e':	/* Erase WCS projection in image header */
                    erasewcs++;
                    break;

                    case 'f':	/* Write FITS file */
                    fitsout = 1;
                    break;

                case 'g':	/* Guide Star object class */
                    if (ac < 2)
                        PrintUsage (str);
                    setgsclass ((int) atof (*++av));
                    ac--;
                    break;

                case 'h':	/* Maximum number of reference stars */
                    if (ac < 2)
                        PrintUsage (str);
                    setmaxcat ((int) atof (*++av));
                    ac--;
                    break;

                case 'i':       /* Image star minimum peak value */
                    if (ac < 2)
                        PrintUsage (str);
                    bmin = atof (*++av);
                    if (bmin < 0)
                        setstarsig (-bmin);
                    else
                        setbmin (bmin);
                    ac--;
                    break;

                case 'j':  /* center coordinates on command line in J2000 */
                    if (ac < 3)
                        PrintUsage ("* Missing RA Dec or coordinate system");
                    setsys (WCS_J2000);
                    strcpy (rastr, *++av);
                    ac--;
                    strcpy (decstr, *++av);
                    ac--;
                    setcenter (rastr, decstr);
                    break;

                case 'k':  /* select magnitude to use from reference catalog */
                    if (ac < 2)
                        PrintUsage (str);
                    av++;
                    c1 = (*av)[0];
                    if (c1 > '9')
                        imag = (int) c1;
                    else
                        imag = (int) c1 - 48;
                    setsortmag (imag);
                    ac--;
                    break;

                case 'l':	/* Left-right reflection before rotating */
                    mirror = 1;
                    setmirror (mirror);
                    break;

                case 'm':	/* Limiting reference star magnitude */
                    if (ac < 2)
                        PrintUsage (str);
                    maglim1 = -99.0;
                    maglim2 = atof (*++av);
                    ac--;
                    if (ac > 1 && isnum (*(av+1))) {
                        maglim1 = maglim2;
                        maglim2 = atof (*++av);
                        ac--;
                    }
                    setreflim (maglim1, maglim2);
                    break;

                case 'n':	/* Number of parameters to fit */
                    if (ac < 2)
                        PrintUsage (str);
                    setnfit ((int) atof (*++av));
                    ac--;
                    break;

                case 'o':	/* Specifiy output image filename */
                    if (ac < 2)
                        PrintUsage (str);
                    if (*(av+1)[0] == '-' || *(str+1) != (char)0)
                        overwrite++;
                    else {
                        strcpy (outname, *(av+1));
                        overwrite = 0;
                        av++;
                        ac--;
                    }
                    writeheader++;
                    break;

                case 'p':  /* Initial plate scale in arcseconds per pixel */
                    if (ac < 2)
                        PrintUsage ("* Missing arcseconds per pixel");
                    setsecpix (atof (*++av));
                    ac--;
                    if (ac > 1 && isnum (*(av+1))) {
                        setsecpix2 (atof (*++av));
                        ac--;
                    }
                    break;

                case 'q':	/* Fit again */
                    if (ac < 2)
                        PrintUsage ("* Missing -q option");
                    str1 = *++av;
                    ac--;

                    while ((c1 = *str1) != 0) {
                        switch (c1) {

                        case 'b':	/* Bin star matches for speed */
                            sprintf (errmsg, "* The \"binarray\" setting (-q b) has been removed");
                            PrintUsage (errmsg);
                            break;

                        case 'i':	/* Iterate fit: new area */
                            c2 = *(str1+1);
                            if ((int)c2 > 47 && (int)c2 < 58) {
                                i = (int) c2 - 48;
                                str1++;
                            } else
                                i = 1;
                            setiterate (i);
                            break;

                        case 'n':	/* Increase number of parameters fit */
                            c2 = *(str1+1);
                            if ((int)c2 > 47 && (int)c2 < 58) {
                                i = (int) c2 - 48;
                                str1++;
                            } else
                                i = 1;
                            setnfiterate (i);
                            break;

                        case 'r':	/* Recenter fit and rerun */
                            setrecenter (1);
                            break;

                        case 's':	/* Use only matches within 2 sigma */
                            sprintf (errmsg, "* The \"resid_refine\" setting (-q s) has been removed");
                            PrintUsage (errmsg);
                            break;

                        case 't':	/* Iterate fit: tighten up */
                            c2 = *(str1+1);
                            if ((int)c2 > 47 && (int)c2 < 58) {
                                i = (int) c2 - 48;
                                str1++;
                            } else
                                i = 1;
                            setiteratet (i);
                            break;

                        case 'p':	/* Use polynomial WCS */
                            c2 = *(str1+1);
                            if ((int)c2 > 47 && (int)c2 < 58) {
                                i = (int) c2;
                                str1++;
                            } else
                                i = 6;
                            setfitplate (i);
                            break;

                        case '8':	/* Fit 8 polynomial parameters */
                            setfitplate (8);
                            break;

                        case 'w':	/* Do not rotate image WCS */
                            rotatewcs = 0;
                            break;

                        default:
                            sprintf (errmsg, "* Illegal q option -%s-", str1);
                            PrintUsage (errmsg);
                            break;
                        }

                        str1++;
                    }
                    break;

                case 'r':	/* Angle in degrees to rotate before fitting */
                    if (ac < 2)
                        PrintUsage (str);
                    rot = (int) atof (*++av);
                    setrotate (rot);
                    ac--;
                    break;

                case 's':   /* Fraction image stars over reference stars */
                    if (ac < 2)
                        PrintUsage (str);
                    setfrac (atof (*++av));
                    ac--;
                    break;

                case 't':	/* Tolerance in pixels for star match */
                    if (ac < 2)
                        PrintUsage (str);
                    settolerance (atof (*++av));
                    ac--;
                    break;

                case 'u':	/* File of prematched (x,y)/(ra,dec) */
                    sprintf (errmsg, "* The \"matchcat\" setting (-u) has been removed");
                    PrintUsage (errmsg);
                    break;

                case 'v':	/* More verbosity */
                    verbose++;
                    break;

                case 'w':	/* Update the fields in a new FITS file */
                    writeheader++;
                    break;

                case 'x':	/* X and Y coordinates of reference pixel */
                    if (ac < 3)
                        PrintUsage (str);
                    x = atof (*++av);
                    ac--;
                    y = atof (*++av);
                    ac--;
                    setrefpix (x, y);
                    break;

                case 'y':	/* Multiply dimensions of image by fraction */
                    sprintf (errmsg, "* The \"imfrac\" setting (-y) has been removed");
                    PrintUsage (errmsg);
                    break;

                case 'z':       /* Use AIPS classic WCS */
                    setdefwcs (WCS_ALT);
                    break;

#ifdef CROP_BORDER
                case 'B':	/* Limit shift in pixels of correlation algorithm */
                    if (ac < 2)
                        PrintUsage (str);
                    setcropborder((int) atof (*++av));
                    ac--;
                    break;

                case 'D':	/* Make the cropping region square */
                    setsquareborder(1);
                    break;
#endif /* CROP_BORDER */

#ifdef LIMIT_SHIFT
                case 'C':	/* Limit shift in pixels of correlation algorithm */
                    if (ac < 2)
                        PrintUsage (str);
                    setmaxshift ((int) atof (*++av));
                    ac--;
                    break;
#endif /* LIMIT_SHIFT */

                default:
                    sprintf (errmsg, "* Illegal command -%c-", c);
                    PrintUsage (errmsg);
                    break;
                }
            }
        } else if (isfits (str) || isiraf (str)) {
            if (nfile >= maxnfile) {
                maxnfile = maxnfile * 2;
                fn = (char **) realloc ((void *)fn, maxnfile);
            }

            fn[nfile] = str;
            nfile++;
        } else {
            sprintf (errmsg, "* %s is not a FITS or IRAF file.", str);
            PrintUsage (errmsg);
        }
    }

    /* If reference catalog is not set, exit with an error message */
    if (refcatname == NULL && matchfile == NULL) {
        PrintUsage ("* Must specifiy a reference catalog using -c or alias.");
    }

    if (!writeheader && !verbose) {
        PrintUsage ("* Must have either w or v argument");
    }

    if (nfile > 0) {
        for (ifile = 0; ifile < nfile; ifile++) {
            if ( verbose)
                    printf ("%s:\n", fn[ifile]);
            FitWCS (progname, fn[ifile]);
            if (verbose)
                printf ("\n");
            }
        }

    /* Print error message if no image files to process */
    else
        PrintUsage ("* No files to process.");

    return 0;
}

static void
PrintUsage (command)

char    *command;

{
    fprintf (stderr,"%s %s\n", progname, RevMsg);
    if (version)
        exit (-1);

    if (command != NULL) {
        if (command[0] == '*')
            fprintf (stderr, "%s\n", command);
        else
            fprintf (stderr, "* Missing argument for command %c\n", command[0]);
        exit (1);
        }

    fprintf(stderr,"DASCH\'s version of WCSTool imwcs\n");
    fprintf(stderr,"Usage: [-vwdfl][-o filename][-m mag][-n frac][-s mode][-g class]\n");
    fprintf(stderr,"       [-h maxref][-i peak][-c catalog][-p scale][-b ra dec][-j ra dec]\n");
    fprintf(stderr,"       [-r deg][-t tol][-u matchfile][-x x y][-y frac] FITS or IRAF file(s)\n");
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: initial center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -c: reference catalog (gsc, uac, usac, ujc, tab table file\n");
    fprintf(stderr,"  -d: Use following DAOFIND output catalog instead of search\n");
    fprintf(stderr,"  -e: Erase image WCS keywords\n");
    fprintf(stderr,"  -f: write FITS output no matter what input\n");
    fprintf(stderr,"  -g: Guide Star Catalog class (-1=all,0,3 (default -1)\n");
    fprintf(stderr,"  -h: maximum number of reference stars to use (10-200, default %d\n", MAXSTARS);
    fprintf(stderr,"  -i: minimum peak value for star in image (<0=-sigma)\n");
    fprintf(stderr,"  -j: initial center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -k: magnitude to use (1 to nmag)\n");
    fprintf(stderr,"  -l: reflect left<->right before rotating and fitting\n");
    fprintf(stderr,"  -m: reference catalog magnitude limit(s) (default none)\n");
    fprintf(stderr,"  -n: list of parameters to fit (12345678; negate for refinement)\n");
    fprintf(stderr,"  -o: name for output image, no argument to overwrite\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -q: <i>terate, <r>ecenter, <s>igma clip, <p>olynomial, <t>olerance reduce, <w>do not rotate WCS, <n>more params\n");
    fprintf(stderr,"  -r: rotation angle in degrees before fitting (default 0)\n");
    fprintf(stderr,"  -s: use this fraction extra stars (default 1.0)\n");
    fprintf(stderr,"  -t: offset tolerance in pixels (default %d)\n", PIXDIFF);
    fprintf(stderr,"  -u: [removed in DASCH]\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: write header (default is read-only)\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
    fprintf(stderr,"  -y: [removed in DASCH]\n");
    fprintf(stderr,"  -z: use AIPS classic projections instead of WCSLIB\n");
#ifdef CROP_BORDER
    fprintf(stderr,"  -B: pixels to crop from the border\n");
    fprintf(stderr,"  -D: make the cropping region square\n");
#endif /* CROP_BORDER */
#ifdef LIMIT_SHIFT
    fprintf(stderr,"  -C: pixel limit of correlation algorithm range\n");
#endif /* LIMIT_SHIFT */
    exit (1);
}


static void
FitWCS (progname, name)

char	*progname;	/* Name of program being executed */
char	*name;		/* FITS or IRAF image filename */

{
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    int bpix = 0;
    char *image;		/* Image */
    char *header;		/* FITS header */
    char *irafheader = NULL;	/* IRAF image header */
    char newname[256];		/* Name for revised image */
    char pixname[256];		/* Pixel file name for revised image */
    char temp[16];
    char *ext;
    char *fname;
    int lext, lname;
    int rename = 0;
    char *imext, *imext1;
    char *newimage;

    image = NULL;

    /* Open IRAF image if .imh extension is present */
    if (isiraf (name)) {
        iraffile = 1;
        if ((irafheader = irafrhead (name, &lhead)) != NULL) {
            header = iraf2fits (name, irafheader, lhead, &nbhead);
            if (header == NULL) {
                fprintf (stderr, "Cannot translate IRAF header %s/n",name);
                free (irafheader);
                return;
                }
            if (imsearch || writeheader || rot || mirror) {
                if ((image = irafrimage (header)) == NULL) {
                    hgetm (header,"PIXFIL", 255, pixname);
                    fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
                    free (irafheader);
                    free (header);
                    return;
                    }
                }
            }
        else {
            fprintf (stderr, "Cannot read IRAF header file %s\n", name);
            return;
            }
        }

    /* Open FITS file if .imh extension is not present */
    else {
        iraffile = 0;
        fitsout = 1;
        if ((header = fitsrhead (name, &lhead, &nbhead)) != NULL) {
            if (imsearch || rot || mirror) {
                if ((image = fitsrimage (name, nbhead, header)) == NULL) {
                    fprintf (stderr, "Cannot read FITS image %s\n", name);
                    free (header);
                    return;
                    }
                }
            }
        else {
            fprintf (stderr, "Cannot read FITS file %s\n", name);
            return;
            }
        }

    if (erasewcs) {
        if (strchr (name, ',') || strchr (name,'['))
            setheadshrink (0);

        DelWCSFITS (header, verbose);
    }

    /* Rotate and/or reflect image */
    if ((imsearch || writeheader) && (rot != 0 || mirror)) {
        if ((newimage = RotFITS (name,header,image,0,0,rot,mirror,bitpix,
                                 rotatewcs,verbose))
            == NULL) {
            fprintf (stderr,"Image %s could not be rotated\n", name);
            if (iraffile)
                free (irafheader);
            if (image != NULL)
                free (image);
            free (header);
            return;
            }
        else {
            if (image != NULL)
                free (image);
            image = newimage;
            }

        if (!overwrite)
            rename = 1;
        }

    /* Check for permission to overwrite */
    else if (overwrite)
        rename = 0;
    else
        rename = 1;

    /* Use output filename if it is set on the command line */
    if (outname[0] > 0)
        strcpy (newname, outname);

    /* Make up name for new FITS or IRAF output file */
    else if (rename) {

    /* Remove directory path and extension from file name */
        ext = strrchr (name, '.');
        fname = strrchr (name, '/');
        if (fname)
            fname = fname + 1;
        else
            fname = name;
        lname = strlen (fname);
        if (ext) {
            lext = strlen (ext);
            strncpy (newname, fname, sizeof(newname));
            *(newname + lname - lext) = 0;
            }
        else
            strcpy (newname, fname);

    /* Add image extension number or name to output file name */
        imext = strchr (fname, ',');
        imext1 = NULL;
        if (imext == NULL) {
            imext = strchr (fname, '[');
            if (imext != NULL) {
                imext1 = strchr (fname, ']');
                *imext1 = (char) 0;
                }
            }
        if (imext != NULL) {
            strcat (newname, "_");
            strcat (newname, imext+1);
            }

    /* Add rotation and reflection to image name */
        if (mirror)
            strcat (newname, "m");
        else if (rot != 0)
            strcat (newname, "r");
        if (rot < 10 && rot > -1)
            sprintf (temp,"%1d",rot);
        else if (rot < 100 && rot > -10)
            sprintf (temp,"%2d",rot);
        else if (rot < 1000 && rot > -100)
            sprintf (temp,"%3d",rot);
        else
            sprintf (temp,"%4d",rot);
        if (rot != 0)
            strcat (newname, temp);

    /* Add file extension preceded by a w */
        if (fitsout)
            strcat (newname, "w.fits");
        else {
            strcpy (pixname, "HDR$");
            strcat (pixname, newname);
            strcat (pixname, "w.pix");
            hputm (header, "PIXFIL", pixname);
            strcat (newname, "w.imh");
            }
        }
    else
        strcpy (newname, name);

    if (SetWCSFITS (name, header, image, refcatname, verbose)) {
        if (writeheader) {
            if (verbose)
                (void) PrintWCS (header, verbose);	/* print new WCS */

        /* Log WCS program version in the image header */
            hputs (header,"IMWCS",RevMsg);
            hgeti4 (header, "BITPIX", &bpix);

            if (fitsout) {
                if (bpix == 0) {
                    if (fitswhead (newname, header) > 0 && verbose) {
                        if (overwrite)
                            printf ("%s: rewritten successfully.\n", newname);
                        else
                            printf ("%s: written successfully.\n", newname);
                        }
                    else if (verbose)
                        printf ("%s could not be written.\n", newname);
                    }
                else if (image == NULL) {
                    if (fitscimage (newname, header, name) > 0 && verbose) {
                        if (overwrite)
                            printf ("%s: rewritten successfully.\n", newname);
                        else
                            printf ("%s: written successfully.\n", newname);
                        }
                    else if (verbose)
                        printf ("%s could not be written.\n", newname);
                    }
                else {
                    if (fitswimage (newname, header, image) > 0 && verbose) {
                        if (overwrite)
                            printf ("%s: rewritten successfully.\n", newname);
                        else
                            printf ("%s: written successfully.\n", newname);
                        }
                    else if (verbose)
                        printf ("%s could not be written.\n", newname);
                    }
                }
            else if (rename) {
                if (irafwimage (newname,lhead,irafheader,header,image) > 0 && verbose) {
                    if (overwrite)
                        printf ("%s: rewritten successfully.\n", newname);
                    else
                        printf ("%s: written successfully.\n", newname);
                    }
                else if (verbose)
                    printf ("%s could not be written.\n", newname);
                }
            else {
                if (irafwhead (newname,lhead,irafheader,header) > 0 && verbose) {
                    if (overwrite)
                        printf ("%s: rewritten successfully.\n", newname);
                    else
                        printf ("%s: written successfully.\n", newname);
                    }
                else if (verbose)
                    printf ("%s could not be written.\n", newname);
                }
            }
        else if (verbose)
            printf ("%s: file unchanged.\n", name);
        }
    else if (verbose)
        printf ("%s: file unchanged.\n", name);

    free (header);
    if (iraffile)
        free (irafheader);
    if (image != NULL)
        free (image);
    return;
}

// We need to define this function for libwcs matchstar.o, even though we don't
// use it here.
void 
CheckRunTimer(void) 
{
}