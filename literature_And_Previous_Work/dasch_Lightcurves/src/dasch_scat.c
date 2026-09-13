// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* dasch_scat.c
 *
 * This program runs Jessica Mink's the "scat" function on the DASCH MySQL database.
 * Adapted from WCSTools 3.8.1, 14 December 2009, Jessica Mink SAO
 *
 *  both fitsfile.c, hput.c, proj.c and wcslib.c need special changes.
 *
 * Apr  4, 2011 Edward J. Los - Add "pending" selector for unscanned mosaics.
 * Sep 13, 2011 Edward J. Los - Add "-B" for bin one center search
 * Mar 20, 2012 Edward J. Los - Support sequestered data
 * Mar 16, 2020 Edward J. Los - Correct number of arguments required for '-P'
 *                              Add "-L" so that the showtext web page can output a single column plate list without duplicates.
 *
 * cc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   dasch_scat.c pipelineutils.a -ltable -lutil  -lwcs -o dasch_scat -L/usr/lib/mysql -lmysqlclient
 *
 *
 * For M44 web access
 *  dasch_scat dateform=j -W showplate.php extractimage.php extracttarball.php FITS TAR -r 200 -s n -T wcsfit -j -n 2000 -c scanner 08:40:24.000 +19:40:58.80 J2000 -H localhost.localdomain -U scanweb -P arzga5R
 *
 * For M44 web text access
 *  dasch_scat -h -t dateform=j -r 200 -s n     -T wcsfit  -j -n 2000 -c scanner 08:40:24.000  19:40:58.80 J2000 -H localhost.localdomain -U scanweb -P arzga5R
 *
 * For an M44 plate list:
 *
 * dasch_scat dateform=f -r 200 -s n -T wcsfit -j -n 5000 -c scanner 12:29:06.700 +2:03:08.60 J2000 -H localhost.localdomain -U scanweb -M -P arzga5R
 *
 * For a complete plate list
 * dasch_scat dateform=f -r 200 -s n -T wcsfit -j -n 20000 -c scanner  -H localhost.localdomain -U scanweb -M -P arzga5R
 */
#include <math.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include "searchgsc.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "libwcs/wcs.h"
#include "libwcs/lwcs.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcscat.h"
#include "scanread.h"
#define SCANCAT -5

void setminpmqual (int n);
     static void SearchHead();
static int GetArea();
extern void setminid(int minid);

static int afile = 0;		/* True to append output file */
static double epoch2 = 0.0;	/* Latest epoch for catalog search */
static double epoch1 = 0.0;	/* Earliest epoch for catalog search */
static int classd = -1;		/* Guide Star Catalog object classes */
static int wfile = 0;		/* True to print output file */
static int tabout = 0;		/* 1 for tab table to standard output */
static int webtable = 0; /* Use web format for output */
static char *refcatname[5];	/* reference catalog names */
static int refcat = 0;		/* reference catalog switch */
static int ncat = 0;		/* Number of reference catalogs to search */
static int nmag = 0;	/* Number of magnitudes in reference catalog */
static int nmagmax = MAXNMAG;
static double ra0 = -99.0;	/* Initial center RA in degrees */
static double dec0 = -99.0;	/* Initial center Dec in degrees */
static int syscoor = 0;		/* Input search coordinate system */
static double eqcoor = 0.0;	/* Equinox of search center */
static double epoch0 = 0.0;	/* Epoch for coordinates */
static double eqout = 0.0;	/* Equinox for output coordinates */
static int idrun = 0;		/* If 1, 2MASS ID run from inside loop */
static char *ranges;		/* Catalog numbers to print */
static char *coorsys;		/* Coordinate system of search center */
static int match = 0;		/* If 1, match num exactly in BIN or ASC cats*/
static int readlist = 0;	/* If 1, search centers are from a list */
static char *listfile;		/* Name of catalog file with search centers */
static int debug = 0;		/* True for extra information */
static int webdump = 0;
static int verbose = 0;		/* Verbose/debugging flag */
static int printprog = 0;	/* 1 to print program name and version */
static int catsort = SORT_UNSET; /* Default to sort stars by magnitude */
static int closest;		/* 1 if printing only closest star */
static int sysout0 = 0;		/* Output coordinate system */
static char title[80];	/* Title of reference Catalog */
static int sysref = 0;	/* Coordinate system of reference catalog */
static double eqref;	/* Equinox of catalog to be searched */
static double epref;	/* Epoch of catalog to be searched */
static int mprop = 0;	/* 1 if proper motion in reference catalog */
static int degout = 0;	/* Set to 1 to print coordinates in degrees */
static int degout0 = 0;		/* 1 if degrees output instead of hms */
static int searchcenter = 0;	/* 1 to print simpler format */
static int printabhead = 1;	/* 1 to print tab table heading if no sources */
static int printhead = 0;	/* 1 to print table heading */
static int printobj = 0;	/* If 1, print object name instead of number */
static int lofld = 0;		/* Length of object name field in output */
static char *keyword;		/* Column to add to tab table output */
static int oneline = 0;		/* If 1, print center and closest on 1 line */
static int nstars = 0;		/* Number of brightest stars to list */
static int sortmag = 0;
static double maglim1 = MAGLIM1; /* Catalog bright magnitude limit */
static double maglim2 = MAGLIM2; /* Catalog faint magnitude limit */
static char *objname;		/* Object name for output */
static double rad0 = 0.0;	/* Search box radius */
static double rad1 = 0.0;	/* Inner search annulus radius */
static int rdra = 0;		/* If 1, dra is in ra units, not sky units */
static double dra0 = 0.0;	/* Search box width */
static double ddec0 = 0.0;	/* Search box height */
static int votab = 0;		/* If 1, print output as VOTable XML */
static int printxy = 0;		/* If 1, print X Y instead of object number */
static char *xstr, *ystr;	/* X and Y strings if printxy */
static char *progname;		/* Name of program as executed */
static char *RevMsg = "Adapted from WCSTools 3.8.1, 14 December 2009, Jessica Mink SAO";
static char voerror[80];	/* Error for Virtual Observatory */
static int nddeg = 7;		/* Number of decimal places in degree output */
static int ndra = 3;		/* Number of decimal places in RA seconds */
static int nddec = 2;		/* Number of decimal places in Dec seconds */
static int nalloc = 0;

static double **gm;		/* Catalog magnitudes */
static double *gra;		/* Catalog star right ascensions */
static double *gdec;		/* Catalog star declinations */
static double *gpra;		/* Catalog star RA proper motions */
static double *gpdec;		/* Catalog star declination proper motions */
static double *gnum;		/* Catalog star numbers */
static int *gc;			/* Catalog star object classes */
static double *gx;		/* Catalog star X positions on image */
static double *gy;		/* Catalog star Y positions on image */
static char **gobj;		/* Catalog star object names */
static char **gobj1;		/* Catalog star object names */
static int nohead = 1;		/* 1 to print table heading */
static int notprinted = 1;	/* If 1, print header */
static struct Star *srch;	/* Search center structure for catalog search */
static struct StarCat *starcat[5]; /* Star catalog data structure */
static int printepoch = 0;	/* 1 to print epoch of entry */
static struct StarCat *srchcat; /* Search catalog structure */
static int padspt = 0;		/* Set to one to pad out long spectral type */
static char cpname[16];		/* Name of program for error messages */
static int minpmqual = 3; /* Proper motion quality limit (0=bad, 9=good)*/
static int minid = 0; /* Minimum number of plate ID's (<0 excludes Tycho-2) */
static char *revmessage = NULL;	/* Version and date for calling program */
static char *revmsg0 = "";
static int http=0;		/* Set to one if http header needed on output */

extern void DS_setlimdeg (int degoutx); /* scanread.c */

/* From catutil.c */
/* DS_rgetn -- Return number of values from range structure */

int
DS_rgetn (range)

         struct Range *range;	/* Range structure */

{
    return (range->nvalues);
}
/* from catutil.c */

/* DS_isrange -- Return 1 if string is a range, else 0 */

int
DS_isrange (string)

         char *string;		/* String which might be a range of numbers */

{
    int i, lstr;

    /* If string is NULL or empty, return 0 */
    if (string == NULL || strlen (string) == 0)
        return (0);

    /* If range separators present, check to make sure string is range */
    else if (strchr (string+1, '-') || strchr (string+1, ',')) {
        lstr = strlen (string);
        for (i = 0; i < lstr; i++) {
        if (strchr ("0123456789-,.x", (int)string[i]) == NULL)
                return (0);
        }
        return (1);
    }
    else
        return (0);
}
/* from catutil.c (MODIFIED) */
/* Return code for reference catalog or its type */

int
DS_CatCode (refcatname)

         char	*refcatname;	/* Name of reference catalog */
{
    int refcat;

    refcat = 0;
    if (refcatname == NULL)
        refcat = 0;
    else if (strlen (refcatname) < 1)
        refcat = 0;
    else if (strncasecmp(refcatname,"scanner",7)==0 &&
                     strcsrch(refcatname, ".tab") == NULL)
        refcat = SCANCAT;
    else
        refcat = 0;
    return refcat;
}



/* from catutil.c */
/* DS_RangeInit -- Initialize range structure from string */

struct Range *
DS_RangeInit (string, ndef)

         char	*string;	/* String containing numbers separated by , and - */
         int	ndef;		/* Maximum allowable range value */

{
    struct Range *range;
    int ip, irange;
    char *slast;
    double first, last, step;

    if (!DS_isrange (string) && !isnum (string))
        return (NULL);
    ip = 0;
    range = (struct Range *)calloc (1, sizeof (struct Range));
    range->irange = -1;
    range->nvalues = 0;
    range->nranges = 0;

    for (irange = 0; irange < MAXRANGE; irange++) {

        /* Default to entire list */
        first = 1.0;
        last = ndef;
        step = 1.0;

        /* Skip delimiters to start of range */
        while (string[ip] == ' ' || string[ip] == '	' ||
                     string[ip] == ',')
        ip++;

        /* Get first limit
         * Must be a number, '-', 'x', or EOS.  If not return ERR */
        if (string[ip] == (char)0) {	/* end of list */
        if (irange == 0) {

                /* Null string defaults */
                range->ranges[0] = first;
                if (first < 1)
                    range->ranges[1] = first;
                else
                    range->ranges[1] = last;
                range->ranges[2] = step;
                range->nvalues = range->nvalues + 1 +
                    ((range->ranges[1]-range->ranges[0])/step);
                range->nranges++;
                return (range);
            }
        else
                return (range);
        }
        else if (string[ip] > (char)47 && string[ip] < 58) {
        first = strtod (string+ip, &slast);
        ip = slast - string;
        }
        else if (strchr ("-:x", string[ip]) == NULL) {
        free (range);
        return (NULL);
        }

        /* Skip delimiters */
        while (string[ip] == ' ' || string[ip] == '	' ||
                     string[ip] == ',')
        ip++;

        /* Get last limit
         * Must be '-', or 'x' otherwise last = first */
        if (string[ip] == '-' || string[ip] == ':') {
        ip++;
        while (string[ip] == ' ' || string[ip] == '	' ||
                         string[ip] == ',')
                ip++;
        if (string[ip] == (char)0)
                last = first + ndef;
        else if (string[ip] > (char)47 && string[ip] < 58) {
                last = strtod (string+ip, &slast);
                ip = slast - string;
            }
        else if (string[ip] != 'x')
                last = first + ndef;
        }
        else if (string[ip] != 'x')
        last = first;

        /* Skip delimiters */
        while (string[ip] == ' ' || string[ip] == '	' ||
                     string[ip] == ',')
        ip++;

        /* Get step
         * Must be 'x' or assume default step. */
        if (string[ip] == 'x') {
        ip++;
        while (string[ip] == ' ' || string[ip] == '	' ||
                         string[ip] == ',')
                ip++;
        if (string[ip] == (char)0)
                step = 1.0;
        else if (string[ip] > (char)47 && string[ip] < 58) {
                step = strtod (string+ip, &slast);
                ip = slast - string;
            }
        else if (string[ip] != '-' && string[ip] != ':')
                step = 1.0;
        }

        /* Output the range triple */
        range->ranges[irange*3] = first;
        range->ranges[irange*3 + 1] = last;
        range->ranges[irange*3 + 2] = step;
        range->nvalues = range->nvalues + ((last-first+(0.1*step)) / step + 1);
        range->nranges++;
    }

    return (range);
}
/* from catutil.c MODIFIED */

/* Return code for reference catalog or its type */

int
DS_RefCat (refcatname, title, syscat, eqcat, epcat, catprop, nmag)

         char	*refcatname;	/* Name of reference catalog */
         char	*title;		/* Description of catalog (returned) */
         int	*syscat;	/* Catalog coordinate system (returned) */
         double	*eqcat;		/* Equinox of catalog (returned) */
         double	*epcat;		/* Epoch of catalog (returned) */
         int	*catprop;	/* 1 if proper motion in catalog (returned) */
         int	*nmag;		/* Number of magnitudes in catalog (returned) */
{
    int refcat;

    *catprop = 0;

    refcat = DS_CatCode (refcatname);
    if (refcat == SCANCAT) {
        strcpy (title, refcatname);
        strcat (title, " Scanner MysQL database");
        *syscat = WCS_J2000;
        *eqcat = 2000.0;
        *epcat = 2000.0;
        *catprop = 0;
        *nmag = 1;
    }
    return refcat;
}


void PrintUsage (char *command) /* Command where error occurred or NULL */
{
    FILE *dev;	/* Output, stderr for command line, stdout for web */
    char *srcname;
    char catname[] = "scanner";

    dev = stderr;

    /* Print program name and version */
    fprintf (dev,"%s %s\n", progname, RevMsg);
    if (command != NULL && !strncasecmp (command, "ver", 3))
        exit (0);

    if (command != NULL) {
        if (command[0] == '*')
        fprintf (dev, "%s\n", command);
        else
        fprintf (dev, "* Missing argument for command: %c\n", command[0]);
        exit (1);
    }

    srcname = catname;
    fprintf (stderr,"List %s in a region on the sky\n", srcname);

    fprintf (dev,"Usage: %s [arguments] ra dec system (J2000, B1950, etc.)\n", progname);
    fprintf (dev,"  or : %s [arguments] list of catalog number ranges\n",
                     progname);
    fprintf (dev,"  or : %s [arguments] @file of either positions or numbers)\n",
                     progname);
    fprintf(dev,"  -a: List single closest catalog source\n");
    fprintf(dev,"  -b: Output B1950 (FK4) coordinates\n");
    if (!strcmp (progname, "scat"))
    fprintf(dev,"  -c name: Reference catalog (act, gsc, ua2, usa2, or local file\n");
    fprintf(dev,"  -d: Output RA and Dec in degrees instead of hms dms\n");
    fprintf(dev,"  -e: Output ecliptic coordinates\n");
    fprintf(dev,"  -f: Output search center for other programs\n");
    fprintf(dev,"  -g: Output galactic coordinates\n");
    fprintf(dev,"  -h: Print heading, else do not \n");
    fprintf(dev,"  -i [length]: Print catalog object name, not catalog number (length optional)\n");
    fprintf(dev,"  -j: Output J2000 (FK5) coordinates\n");
    fprintf(dev,"  -k kwd: Add this keyword to output from tab table search\n");
    fprintf(dev,"  -l: Print center and closest star on one line\n");
    fprintf(dev,"  -mx mag1[,mag2]: Magnitude #x limit(s) (only one set allowed, default none) \n");
    fprintf(dev,"  -n num: Number of brightest stars to print (-1=all as found)\n");
    fprintf(dev,"  -o name: Object name \n");
    fprintf(dev,"  -q year: Equinox of output positions in FITS date format or years\n");
    fprintf(dev,"  -r rad: Search radius (<0=-half-width) in arcsec\n");
    fprintf(dev,"  -r radi-rado: Inner and outer edges of search annulus in arcsec\n");
    fprintf(dev,"  -r dx,dy: Search halfwidths in ra,dec in great circle arcseconds\n");
    fprintf(dev,"  -rr dra,ddec: Search halfwidths in ra,dec in arcsec of RA and Dec\n");
    fprintf(dev,"  -s d|e|mx|n|p|r: Sort by r=RA d=Dec mx=Mag#x n=none p=distance e=merge\n");
    fprintf(dev,"  -t: Tab table to standard output as well as file\n");
    fprintf(dev,"  -u x y: Print x y instead of number in front of non-tab entry\n");
    fprintf(dev,"  -v: Verbose\n");
    fprintf(dev,"  -w: Write output file search[objname].[catalog]\n");
    if (!strcmp (progname, "scat") || !strcmp (progname, "sgsc"))
        fprintf(dev,"  -x type: GSC object type (0=stars 3=galaxies -1=all -2=bands)\n");
    fprintf(dev,"  -y year: Epoch of output positions in FITS date format or years\n");
    fprintf(dev,"     year,year: First and last acceptable catalog entry epochs\n");
    fprintf(dev,"  -z: Append to output file search[objname].[catalog]\n");
    fprintf(dev,"   x: Number of magnitude must be same for sort and limits\n");
    fprintf(dev,"      and x may be omitted from either or both -m and -s m\n");
    fprintf(dev,"  -B Search so that plate center is in bin 1\n");
    fprintf(dev,"  -C class: Plate Class for the MysQL catalog\n");
    fprintf(dev,"  -H hostname: Host name for the MysQL catalog\n");
    fprintf(dev,"  -L display first column (plate) only\n");
    fprintf(dev,"  -M display mosaic number instead of exposure number\n");
    fprintf(dev,"  -N number: Plate Number for the MysQL catalog\n");
    fprintf(dev,"  -O Search inner bins of unfitted plates\n");
    fprintf(dev,"  -P userpassword: User Password for the MysQL catalog\n");
    fprintf(dev,"  -S series: Plate Series for the MysQL catalog\n");
    fprintf(dev,"  -T [all|scanned|pending|wcsfit]: Search Type for the MysQL catalog\n");
    fprintf(dev,"  -U username: User name for the MysQL catalog\n");
    fprintf(dev,"  -W use web format output.\n");

    if (command != NULL)
        exit (1);
    else
        exit (0);
}

/* from catutil.c */
/*  RGETR8 -- Return next number from range structure as 8-byte f.p. number */

double
DS_rgetr8 (range)

         struct Range *range;	/* Range structure */

{
    int i;

    if (range == NULL)
        return (0.0);
    else if (range->irange < 0) {
        range->irange = 0;
        range->first = range->ranges[0];
        range->last = range->ranges[1];
        range->step = range->ranges[2];
        range->value = range->first;
    }
    else {
        range->value = range->value + range->step;
        if (range->value > (range->last + (range->step * 0.5))) {
        range->irange++;
        if (range->irange < range->nranges) {
                i = range->irange * 3;
                range->first = range->ranges[i];
                range->last = range->ranges[i+1];
                range->step = range->ranges[i+2];
                range->value = range->first;
            }
        else
                range->value = 0.0;
        }
    }
    return (range->value);
}

/* from catutil.c */
double
DS_CatRad (refcat)

         int	refcat;		/* Catalog code */
{
    if (refcat==GSC || refcat==GSCACT || refcat==UJC || refcat==USAC ||
            refcat==USA1 || refcat==USA2 ||
            refcat == UCAC1 || refcat == UCAC2 || refcat == UCAC3)
        return (900.0);
    else if (refcat==UAC  || refcat==UA1  || refcat==UA2)
        return (120.0);
    else if (refcat == UB1 || refcat==SDSS)
        return (120.0);
    else if (refcat==GSC2)
        return (120.0);
    else if (refcat==TMPSC || refcat==TMPSCE || refcat==TMIDR2)
        return (120.0);
    else if (refcat==TMXSC)
        return (900.0);
    else if (refcat==GSC2)
        return (120.0);
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat == SKY2K ||
                     refcat==SKYBOT)
        return (5000.0);
    else
        return (1800.0);
}
/* from catutil.c */

/* Return number of decimal places in a number */

int
DS_NumNdec (number)

         double number;	/* Floating point number */
{
    char nstring[16];
    char format[16];
    int fracpart;
    int ndec, ndmax;
    double shift;

    if (number < 10.0) {
        ndmax = 12;
        shift = 1000000000000.0;
    }
    else if (number < 100.0) {
        ndmax = 11;
        shift = 100000000000.0;
    }
    else if (number < 1000.0) {
        ndmax = 10;
        shift = 10000000000.0;
    }
    else if (number < 10000.0) {
        ndmax = 9;
        shift = 1000000000.0;
    }
    else if (number < 100000.0) {
        ndmax = 8;
        shift = 100000000.0;
    }
    else if (number < 1000000.0) {
        ndmax = 7;
        shift = 10000000.0;
    }
    else if (number < 10000000.0) {
        ndmax = 6;
        shift = 1000000.0;
    }
    else if (number < 100000000.0) {
        ndmax = 5;
        shift = 100000.0;
    }
    else if (number < 1000000000.0) {
        ndmax = 4;
        shift = 10000.0;
    }
    else if (number < 10000000000.0) {
        ndmax = 3;
        shift = 1000.0;
    }
    else if (number < 100000000000.0) {
        ndmax = 2;
        shift = 100.0;
    }
    else if (number < 1000000000000.0) {
        ndmax = 1;
        shift = 10.0;
    }
    else
        return (0);
    fracpart = (int) (((number - floor (number)) * shift) + 0.5);
    sprintf (format, "%%0%dd", ndmax);
    sprintf (nstring, format, fracpart);
    for (ndec = ndmax; ndec > 0; ndec--) {
        if (nstring[ndec-1] != '0')
        break;
    }
    return (ndec);
}


/* from catutil.c */

int
DS_CatNumLen (refcat, maxnum, nndec)

         int	refcat;		/* Catalog code */
         double	maxnum;		/* Maximum ID number */
         /* (Ignored for standard catalogs) */
         int	nndec;		/* Number of decimal places ( >= 0) */

{
    int ndp;		/* Number of characters for decimal point */

    /* USNO A1.0, A2.0, SA1.0, or SA2.0 Catalogs */
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
            refcat == UAC  || refcat == UA1  || refcat == UA2)
        return (13);

    /* USNO-B1.0 and YB6 */
    else if (refcat == UB1 || refcat == YB6)
        return (12);

    /* GSC II */
    else if (refcat == GSC2)
        return (13);

    /* 2MASS Point Source Catalog */
    else if (refcat == TMPSC || refcat == TMPSCE)
        return (11);
    else if (refcat == TMIDR2)
        return (10);

    /* 2MASS Extended Source Catalog */
    else if (refcat == TMXSC)
        return (11);

    /* UCAC1 Catalog */
    else if (refcat == UCAC1)
        return (10);

    /* UCAC2 Catalog */
    else if (refcat == UCAC2)
        return (10);

    /* UCAC3 Catalog */
    else if (refcat == UCAC3)
        return (10);

    /* USNO Plate Catalogs */
    else if (refcat == USNO)
        return (7);

    /* USNO UJ 1.0 Catalog */
    else if (refcat == UJC)
        return (12);

    /* SDSS Catalog */
    else if (refcat == SDSS)
        return (18);

    /* SkyBot Objects */
    else if (refcat == SKYBOT)
        return (6);

    /* HST Guide Star Catalog */
    else if (refcat == GSC || refcat == GSCACT)
        return (9);

    /* SAO, PPM, Hipparcos, or IRAS Point Source Catalogs (TDC binary format) */
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==BSC ||
                     refcat==HIP)
        return (6);

    /* SKY2000 Catalog (TDC binary format) */
    else if (refcat==SKY2K)
        return (7);

    /* Tycho, Tycho2, or ACT catalogs */
    else if (refcat == TYCHO || refcat == TYCHO2 ||
                     refcat == TYCHO2E || refcat == ACT)
        return (10);

    /* Starbase tab-separated, TDC binary, or TDC ASCII catalogs */
    else {
        if (nndec > 0)
        ndp = 1;
        else {
        if ((nndec = DS_NumNdec (maxnum)) > 0)
                ndp = 1;
        else
                ndp = 0;
        }
        if (maxnum < 10.0)
        return (1 + nndec + ndp);
        else if (maxnum < 100.0)
        return (2 + nndec + ndp);
        else if (maxnum < 1000.0)
        return (3 + nndec + ndp);
        else if (maxnum < 10000.0)
        return (4 + nndec + ndp);
        else if (maxnum < 100000.0)
        return (5 + nndec + ndp);
        else if (maxnum < 1000000.0)
        return (6 + nndec + ndp);
        else if (maxnum < 10000000.0)
        return (7 + nndec + ndp);
        else if (maxnum < 100000000.0)
        return (8 + nndec + ndp);
        else if (maxnum < 1000000000.0)
        return (9 + nndec + ndp);
        else if (maxnum < 10000000000.0)
        return (10 + nndec + ndp);
        else if (maxnum < 100000000000.0)
        return (11 + nndec + ndp);
        else if (maxnum < 1000000000000.0)
        return (12 + nndec + ndp);
        else if (maxnum < 10000000000000.0)
        return (13 + nndec + ndp);
        else
        return (14 + nndec + ndp);
    }
}

/* from catutil.c */
void
DS_CatNum (refcat, nnfld, nndec, dnum, numstr)

         int	refcat;		/* Catalog code */
         int	nnfld;		/* Number of characters in number (from DS_CatNumLen) */
         /* Print leading zeroes if negative */
         int	nndec;		/* Number of decimal places ( >= 0) */
         /* Omit leading spaces if negative */
         double	dnum;		/* Catalog number of source */
         char	*numstr;	/* Formatted number (returned) */

{
    char nform[128];	/* Format for star number */
    int lnum, i;

    /* USNO A1.0, A2.0, SA1.0, or SA2.0 Catalogs */
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
            refcat == UAC  || refcat == UA1  || refcat == UA2) {
        if (nnfld < 0)
        sprintf (numstr, "%013.8f", dnum);
        else
        sprintf (numstr, "%13.8f", dnum);
    }

    /* USNO-B1.0  and USNO-YB6 */
    else if (refcat == UB1 || refcat == YB6) {
        if (nnfld < 0)
        sprintf (numstr, "%012.7f", dnum);
        else
        sprintf (numstr, "%12.7f", dnum);
    }

    /* USNO-UCAC1 */
    else if (refcat == UCAC1) {
        if (nnfld < 0)
        sprintf (numstr, "%010.6f", dnum);
        else
        sprintf (numstr, "%10.6f", dnum);
    }

    /* USNO-UCAC2 */
    else if (refcat == UCAC2) {
        if (nnfld < 0)
        sprintf (numstr, "%010.6f", dnum);
        else
        sprintf (numstr, "%10.6f", dnum);
    }

    /* USNO-UCAC3 */
    else if (refcat == UCAC3) {
        if (nnfld < 0)
        sprintf (numstr, "%010.6f", dnum);
        else
        sprintf (numstr, "%10.6f", dnum);
    }

    /* SDSS */
    else if (refcat == SDSS) {
        sprintf (numstr, "582%015.0f", dnum);
    }

    /* GSC II */
    else if (refcat == GSC2) {
        if (nnfld < 0) {
        if (dnum > 0)
                sprintf (numstr, "N%.0f", (dnum+0.01));
        else
                sprintf (numstr, "S%.0f", (-dnum + 0.01));
        lnum = strlen (numstr);
        if (lnum < -nnfld) {
                for ( i = lnum; i < -nnfld; i++)
                    strcat (numstr, " ");
            }
        }
        else {
        if (dnum > 0)
                sprintf (numstr, "N%.0f", (dnum+0.5));
        else
                sprintf (numstr, "S%.0f", (-dnum + 0.5));
        }
    }

    /* 2MASS Point Source Catalogs */
    else if (refcat == TMPSC || refcat == TMPSCE) {
        if (nnfld < 0)
        sprintf (numstr, "%011.6f", dnum);
        else
        sprintf (numstr, "%11.6f", dnum);
    }

    /* 2MASS Extended Source Catalog */
    else if (refcat == TMXSC) {
        if (nnfld < 0)
        sprintf (numstr, "%011.6f", dnum);
        else
        sprintf (numstr, "%11.6f", dnum);
    }

    /* 2MASS Point Source Catalogs */
    else if (refcat == TMIDR2) {
        if (nnfld < 0)
        sprintf (numstr, "%010.7f", dnum);
        else
        sprintf (numstr, "%10.7f", dnum);
    }

    /* USNO Plate Catalog */
    else if (refcat == USNO) {
        if (nnfld < 0)
        sprintf (numstr, "%07d", (int)(dnum+0.5));
        else
        sprintf (numstr, "%7d", (int)(dnum+0.5));
    }

    /* USNO UJ 1.0 Catalog */
    else if (refcat == UJC) {
        if (nnfld < 0)
        sprintf (numstr, "%012.7f", dnum);
        else
        sprintf (numstr, "%12.7f", dnum);
    }

    /* HST Guide Star Catalog */
    else if (refcat == GSC || refcat == GSCACT) {
        if (nnfld < 0)
        sprintf (numstr, "%09.4f", dnum);
        else
        sprintf (numstr, "%9.4f", dnum);
    }

    /* SAO, PPM, or IRAS Point Source Catalogs (TDC binary format) */
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==BSC ||
                     refcat==HIP) {
        if (nnfld < 0)
        sprintf (numstr, "%06d", (int)(dnum+0.5));
        else
        sprintf (numstr, "%6d", (int)(dnum+0.5));
    }

    /* SKY2000 Catalog (TDC binary format) */
    else if (refcat==SKY2K) {
        if (nnfld < 0)
        sprintf (numstr, "%07d", (int)(dnum+0.5));
        else
        sprintf (numstr, "%7d", (int)(dnum+0.5));
    }


    /* Tycho or ACT catalogs */
    else if (refcat==TYCHO || refcat==TYCHO2 ||
                     refcat == TYCHO2E || refcat==ACT) {
        if (nnfld < 0)
        sprintf (numstr, "%010.5f", dnum);
        else
        sprintf (numstr, "%10.5f", dnum);
    }

    /* Starbase tab-separated, TDC binary, or TDC ASCII catalogs */
    else if (nndec > 0) {
        if (nnfld > 0)
        sprintf (nform,"%%%d.%df", nnfld, nndec);
        else if (nnfld < 0)
        sprintf (nform,"%%0%d.%df", -nnfld, nndec);
        else
        sprintf (nform,"%%%d.%df", nndec+5, nndec);
        sprintf (numstr, nform, dnum);
    }
    else if (nnfld > 10) {
        sprintf (nform,"%%%d.0f", nnfld);
        sprintf (numstr, nform, dnum+0.49);
    }
    else if (nnfld > 0) {
        sprintf (nform,"%%%dd", nnfld);
        sprintf (numstr, nform, (int)(dnum+0.49));
    }
    else if (nnfld < 0) {
        sprintf (nform,"%%0%dd", -nnfld);
        sprintf (numstr, nform, (int)(dnum+0.49));
    }
    else if (nndec < 0)
        sprintf (numstr, "%d", (int)(dnum+0.49));
    else
        sprintf (numstr, "%6d", (int)(dnum+0.49));

    return;
}

/* from catutil.c */
/* Return number of magnitude specified by int of letter */

int
DS_CatMagNum (imag, refcat)

         int	imag;		/* int of magnitude letter */
         int	refcat;		/* Catalog code */
{
    char cmag = (char) imag;	/* Letter name of magnitude */

    /* Make letter upper case */
    if (cmag > 96)
        cmag = cmag - 32;

    if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
            refcat == USAC || refcat == USA1 || refcat == USA2) {
        if (cmag == 'R')
        return (2);
        else
        return (1);	/* B */
    }
    if (refcat == UB1) {
        if (cmag == 'N')
        return (5);
        else if (cmag == 'R')
        return (4);
        else
        return (3);	/* B */
    }
    else if (refcat == YB6) {
        if (cmag == 'K')
        return (5);
        else if (cmag == 'H')
        return (4);
        else if (cmag == 'J')
        return (3);
        else if (cmag == 'R')
        return (2);
        else if (cmag == 'B')
        return (1);
        else
        return (3);	/* J */
    }
    else if (refcat == SKYBOT) {
        return (1);
    }
    else if (refcat == SDSS) {
        if (cmag == 'Z')
        return (5);
        else if (cmag == 'I')
        return (4);
        else if (cmag == 'R')
        return (3);
        else if (cmag == 'G')
        return (2);
        else if (cmag == 'B')
        return (1);
        else
        return (2);	/* G */
    }
    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT) {
        if (cmag == 'B')
        return (1);
        else
        return (2);	/* V */
    }
    else if (refcat==GSC2) {
        if (cmag == 'J')
        return (2);
        else if (cmag == 'N')
        return (3);
        else if (cmag == 'U')
        return (4);
        else if (cmag == 'B')
        return (5);
        else if (cmag == 'V')
        return (6);
        else if (cmag == 'R')
        return (7);
        else if (cmag == 'I')
        return (8);
        else
        return (1);	/* F */
    }
    else if (refcat==TMPSC || refcat == TMXSC) {
        if (cmag == 'J')
        return (1);
        else if (cmag == 'H')
        return (2);
        else
        return (3);	/* K */
    }
    else if (refcat==UCAC2) {
        if (cmag == 'J')
        return (1);
        else if (cmag == 'H')
        return (2);
        else if (cmag == 'K')
        return (3);
        else if (cmag == 'C')
        return (4);
        else
        return (3);	/* K */
    }
    else if (refcat==UCAC3) {
        if (cmag == 'R')
        return (2);
        else if (cmag == 'I')
        return (3);
        else if (cmag == 'J')
        return (4);
        else if (cmag == 'H')
        return (5);
        else if (cmag == 'K')
        return (6);
        else if (cmag == 'M')
        return (7);
        else if (cmag == 'A')
        return (8);
        else
        return (1);	/* B */
    }
    else
        return (1);
}
/* from DS_ctgread.c (MODIFIED) */

/* CTGREAD -- Read ASCII stars in specified region */

int
DS_ctgread (catfile, refcat, distsort, cra, cdec, dra, ddec, drad, dradi,
                 sysout, eqout, epout, mag1, mag2, sortmag, nsmax, starcat,
                 tnum, tra, tdec, tpra, tpdec, tmag, tc, tobj, nlog)

     char	*catfile;	/* Name of reference star catalog file */
     int	refcat;		/* Catalog code from wcscat.h */
     int	distsort;	/* 1 to sort stars by distance from center */
     double	cra;		/* Search center J2000 right ascension in degrees */
     double	cdec;		/* Search center J2000 declination in degrees */
     double	dra;		/* Search half width in right ascension in degrees */
     double	ddec;		/* Search half-width in declination in degrees */
     double	drad;		/* Limiting separation in degrees (ignore if 0) */
     double	dradi;		/* Inner edge of annulus in degrees (ignore if 0) */
     int	sysout;		/* Search coordinate system */
     double	eqout;		/* Search coordinate equinox */
     double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
     double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
     int	sortmag;	/* Number of magnitude by which to limit and sort */
     int	nsmax;		/* Maximum number of stars to be returned */
     struct StarCat **starcat; /* Catalog data structure */
     double	*tnum;		/* Array of ID numbers (returned) */
     double	*tra;		/* Array of right ascensions (returned) */
     double	*tdec;		/* Array of declinations (returned) */
     double	*tpra;		/* Array of right ascension proper motions (returned) */
     double	*tpdec;		/* Array of declination proper motions (returned) */
     double	**tmag;		/* 2-D array of magnitudes (returned) */
     int	*tc;		/* Array of fluxes (returned) */
     char	**tobj;		/* Array of object names (returned) */
     int	nlog;
{
  int nstar;

  nstar = 0;

  if (refcat != TXTCAT) {
    if (refcat == SCANCAT)
      nstar = scanread (catfile, cra,cdec,drad,sysout,eqout,epout,nsmax,starcat,
                                                tnum,tra,tdec,tmag,tobj,nlog);

  }
  return (nstar);
}

/* from catutil.c */
char *
DS_CatName (refcat, refcatname)

         int	refcat;		/* Catalog code */
         char	*refcatname;	/* Catalog file name */
{
    char *catname;

    if (refcat < 1 || refcat > NUMCAT)
        return (refcatname);

    /* Allocate string in which to return a catalog name */
    catname = (char *)calloc (16, 1);

    if (refcat ==  GSC)		/* HST Guide Star Catalog */
        strcpy (catname, "GSC");
    else if (refcat ==  GSCACT)	/* HST GSC revised with ACT */
        strcpy (catname, "GSC-ACT");
    else if (refcat ==  GSC2) {	/* GSC II */
        if (strsrch (refcatname, "22")) {
        strcpy (catname, "GSC 2.2");
        }
        else {
        strcpy (catname, "GSC 2.3");
        }
    }
    else if (refcat == YB6)	/* USNO YB6 Star Catalog */
        strcpy (catname, "USNO-YB6");
    else if (refcat ==  UJC)	/* USNO UJ Star Catalog */
        strcpy (catname, "UJC");
    else if (refcat ==  UAC)	/* USNO A Star Catalog */
        strcpy (catname, "USNO-A2.0");
    else if (refcat ==  USAC)	/* USNO SA Star Catalog */
        strcpy (catname, "USNO-SA2.0");
    else if (refcat ==  SAO)	/* SAO Star Catalog */
        strcpy (catname, "SAO");
    else if (refcat ==  IRAS)	/* IRAS Point Source Catalog */
        strcpy (catname, "IRAS PSC");
    else if (refcat ==  SDSS)	/* Sloan Digital Sky Survey */
        strcpy (catname, "SDSS");
    else if (refcat ==  PPM)	/* PPM Star Catalog */
        strcpy (catname, "PPM");
    else if (refcat ==  TYCHO)	/* Tycho Star Catalog */
        strcpy (catname, "TYCHO");
    else if (refcat ==  UA1)	/* USNO A-1.0 Star Catalog */
        strcpy (catname, "USNO-A1.0");
    else if (refcat ==  UB1)	/* USNO B-1.0 Star Catalog */
        strcpy (catname, "USNO-B1.0");
    else if (refcat ==  UCAC1)	/* USNO UCAC1 Star Catalog */
        strcpy (catname, "USNO-UCAC1");
    else if (refcat ==  UCAC2)	/* USNO UCAC2 Star Catalog */
        strcpy (catname, "USNO-UCAC2");
    else if (refcat ==  UCAC3)	/* USNO UCAC3 Star Catalog */
        strcpy (catname, "USNO-UCAC3");
    else if (refcat ==  UA2)	/* USNO A-2.0 Star Catalog */
        strcpy (catname, "USNO-A2.0");
    else if (refcat ==  USA1)	/* USNO SA-1.0 Star Catalog */
        strcpy (catname, "USNO-SA1.0");
    else if (refcat ==  USA2)	/* USNO SA-2.0 Star Catalog */
        strcpy (catname, "USNO-SA2.0");
    else if (refcat ==  HIP)	/* Hipparcos Star Catalog */
        strcpy (catname, "Hipparcos");
    else if (refcat ==  ACT)	/* USNO ACT Star Catalog */
        strcpy (catname, "ACT");
    else if (refcat ==  BSC)	/* Yale Bright Star Catalog */
        strcpy (catname, "BSC");
    else if (refcat ==  TYCHO2 ||
                     refcat == TYCHO2E)	/* Tycho-2 Star Catalog */
        strcpy (catname, "TYCHO-2");
    else if (refcat ==  TMPSC ||
                     refcat == TMPSCE)	/* 2MASS Point Source Catalog */
        strcpy (catname, "2MASS PSC");
    else if (refcat ==  TMXSC)	/* 2MASS Extended Source Catalog */
        strcpy (catname, "2MASS XSC");
    else if (refcat ==  TMIDR2)	/* 2MASS Point Source Catalog */
        strcpy (catname, "2MASS PSC IDR2");
    else if (refcat ==  SKY2K)	/* SKY2000 Master Catalog */
        strcpy (catname, "SKY2000");
    else if (refcat ==  SKYBOT)	/* SkyBot Solar System Objects */
        strcpy (catname, "SkyBot");
    return (catname);
}
/* from catutil.c */

void
DS_CatID (catid, refcat)

         char	*catid;		/* Catalog ID (returned) */
         int	refcat;		/* Catalog code */
{
    if (refcat == ACT)
        strcpy (catid, "act_id     ");
    else if (refcat == BSC)
        strcpy (catid, "bsc_id    ");
    else if (refcat == GSC || refcat == GSCACT)
        strcpy (catid, "gsc_id    ");
    else if (refcat == GSC2)
        strcpy (catid, "gsc2_id        ");
    else if (refcat == SDSS)
        strcpy (catid, "sdss_id            ");
    else if (refcat == USAC)
        strcpy (catid,"usac_id       ");
    else if (refcat == USA1)
        strcpy (catid,"usa1_id       ");
    else if (refcat == USA2)
        strcpy (catid,"usa2_id       ");
    else if (refcat == UAC)
        strcpy (catid,"usnoa_id      ");
    else if (refcat == UA1)
        strcpy (catid,"usnoa1_id     ");
    else if (refcat == UB1)
        strcpy (catid,"usnob1_id    ");
    else if (refcat == YB6)
        strcpy (catid,"usnoyb6_id   ");
    else if (refcat == UA2)
        strcpy (catid,"usnoa2_id     ");
    else if (refcat == UCAC1)
        strcpy (catid,"ucac1_id  ");
    else if (refcat == UCAC2)
        strcpy (catid,"ucac2_id  ");
    else if (refcat == UCAC3)
        strcpy (catid,"ucac3_id  ");
    else if (refcat == UJC)
        strcpy (catid,"usnoj_id     ");
    else if (refcat == TMPSC || refcat == TMPSCE)
        strcpy (catid,"2mass_id      ");
    else if (refcat == TMXSC)
        strcpy (catid,"2mx_id        ");
    else if (refcat == SAO)
        strcpy (catid,"sao_id ");
    else if (refcat == PPM)
        strcpy (catid,"ppm_id ");
    else if (refcat == IRAS)
        strcpy (catid,"iras_id");
    else if (refcat == TYCHO)
        strcpy (catid,"tycho_id  ");
    else if (refcat == TYCHO2 || refcat == TYCHO2E)
        strcpy (catid,"tycho2_id ");
    else if (refcat == HIP)
        strcpy (catid,"hip_id ");
    else if (refcat == SKY2K)
        strcpy (catid,"sky_id ");
    else if (refcat == SKYBOT)
        strcpy (catid,"skybot_id ");
    else
        strcpy (catid,"id              ");

    return;
}
/* from catutil.c */

char *
DS_ProgName (progpath0)

         char *progpath0;	/* Pathname by which program is invoked */
{
    char *progpath, *progname;
    int i, lpath;

    lpath = (strlen (progpath0) + 2) / 8;
    lpath = (lpath + 1) * 8;;
    progpath = (char *) calloc (lpath, 1);
    strcpy (progpath, progpath0);
    progname = progpath;
    for (i = strlen (progpath); i > -1; i--) {
        if (progpath[i] > 63 && progpath[i] < 90)
            progpath[i] = progpath[i] + 32;
        if (progpath[i] == '/') {
            progname = progpath + i + 1;
            break;
        }
    }
    return (progname);
}
/* from catutil.c */
void
DS_setrevmsg (revmsg)		/* Set version and date string*/
         char *revmsg;
{ revmessage = revmsg; return; }
char *
DS_getrevmsg ()			/* Return version and date string */
{ if (revmessage == NULL) return (revmsg0);
  else return (revmessage); }


/* from catutil.c */

// Our hacked wcstools adds this function, which is needed
// in scanread.c.
extern int getdateform(void);

/* DS_DateString-- Return string with epoch of position in desired format */
char *
DS_DateString (epoch, tabout)

         double	epoch;
         int	tabout;
{
    double year;
    char *temp, *temp1;
    int dateform = getdateform();

    temp = calloc (16, 1);

    if (dateform < 1)
        dateform = EP_MJD;

    if (dateform == EP_EP) {
        if (tabout)
        sprintf (temp, "	%9.4f", epoch);
        else
        sprintf (temp, " %9.4f", epoch);
    }
    else if (dateform == EP_JD) {
        if (epoch == 0.0)
        year = 0.0;
        else
        year = ep2jd (epoch);
        if (tabout)
        sprintf (temp, "	%13.5f", year);
        else
        sprintf (temp, " %13.5f", year);
    }
    else if (dateform == EP_MJD) {
        if (epoch == 0.0)
        year = 0.0;
        else
        year = ep2mjd (epoch);
        if (tabout)
        sprintf (temp, "	%11.5f", year);
        else
        sprintf (temp, " %11.5f", year);
    }
    else {
        if (epoch == 0.0) {
        if (tabout)
                sprintf (temp,"	0000-00-00");
        else
                sprintf (temp," 0000-00-00");
        if (dateform == EP_ISO)
                sprintf (temp,"T00:00");
        }
        else {
        temp1 = ep2fd (epoch);
        if (dateform == EP_FD && strlen (temp1) > 10)
                temp1[10] = (char) 0;
        if (dateform == EP_ISO && strlen (temp1) > 16)
                temp1[16] = (char) 0;
        if (tabout)
                sprintf (temp, "	%s", temp1);
        else
                sprintf (temp, " %s", temp1);
        free (temp1);
        }
    }
    return (temp);
}
/* from catutil.c */
/* Return name of specified magnitude */

void
DS_CatMagName (imag, refcat, magname)

         int	imag;		/* Sequence number of magnitude */
         int	refcat;		/* Catalog code */
         char	*magname;	/* Name of magnitude, returned */
{
    if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
            refcat == USAC || refcat == USA1 || refcat == USA2) {
        if (imag == 2)
        strcpy (magname, "MagR");
        else
        strcpy (magname, "MagB");
    }
    else if (refcat == UB1) {
        if (imag == 5)
        strcpy (magname, "MagN");
        else if (imag == 4)
        strcpy (magname, "MagR2");
        else if (imag == 3)
        strcpy (magname, "MagB2");
        else if (imag == 2)
        strcpy (magname, "MagR1");
        else
        strcpy (magname, "MagB1");
    }
    else if (refcat == YB6) {
        if (imag == 5)
        strcpy (magname, "MagK");
        else if (imag == 4)
        strcpy (magname, "MagH");
        else if (imag == 3)
        strcpy (magname, "MagJ");
        else if (imag == 2)
        strcpy (magname, "MagR");
        else
        strcpy (magname, "MagB");
    }
    else if (refcat == SDSS) {
        if (imag == 5)
        strcpy (magname, "Magz");
        else if (imag == 4)
        strcpy (magname, "Magi");
        else if (imag == 3)
        strcpy (magname, "Magr");
        else if (imag == 2)
        strcpy (magname, "Magg");
        else
        strcpy (magname, "Magu");
    }
    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT) {
        if (imag == 2)
        strcpy (magname, "MagV");
        else
        strcpy (magname, "MagB");
    }
    else if (refcat==TYCHO2E) {
        if (imag == 1)
        strcpy (magname, "MagB");
        else if (imag == 3)
        strcpy (magname, "MagBe");
        else if (imag == 4)
        strcpy (magname, "MagVe");
        else
        strcpy (magname, "MagV");
    }
    else if (refcat==GSC2) {
        if (imag == 2)
        strcpy (magname, "MagJ");
        else if (imag == 3)
        strcpy (magname, "MagN");
        else if (imag == 4)
        strcpy (magname, "MagU");
        else if (imag == 5)
        strcpy (magname, "MagB");
        else if (imag == 6)
        strcpy (magname, "MagV");
        else if (imag == 7)
        strcpy (magname, "MagR");
        else if (imag == 8)
        strcpy (magname, "MagI");
        else
        strcpy (magname, "MagF");
    }
    else if (refcat==SKY2K) {
        if (imag == 1)
        strcpy (magname, "MagB");
        else if (imag == 2)
        strcpy (magname, "MagV");
        else if (imag == 3)
        strcpy (magname, "MagP");
        else
        strcpy (magname, "MagPv");
    }
    else if (refcat==TMPSC || refcat == TMXSC) {
        if (imag == 1)
        strcpy (magname, "MagJ");
        else if (imag == 2)
        strcpy (magname, "MagH");
        else
        strcpy (magname, "MagK");
    }
    else if (refcat==TMPSCE) {
        if (imag == 1)
        strcpy (magname, "MagJ");
        else if (imag == 2)
        strcpy (magname, "MagH");
        else if (imag == 3)
        strcpy (magname, "MagK");
        else if (imag == 4)
        strcpy (magname, "MagJe");
        else if (imag == 5)
        strcpy (magname, "MagHe");
        else if (imag == 6)
        strcpy (magname, "MagKe");
    }
    else if (refcat==UCAC2) {
        if (imag == 1)
        strcpy (magname, "MagJ");
        else if (imag == 2)
        strcpy (magname, "MagH");
        else if (imag == 3)
        strcpy (magname, "MagK");
        else if (imag == 4)
        strcpy (magname, "MagC");
    }
    else if (refcat==UCAC3) {
        if (imag == 1)
        strcpy (magname, "MagB");
        else if (imag == 2)
        strcpy (magname, "MagR");
        else if (imag == 3)
        strcpy (magname, "MagI");
        else if (imag == 4)
        strcpy (magname, "MagJ");
        else if (imag == 5)
        strcpy (magname, "MagH");
        else if (imag == 6)
        strcpy (magname, "MagK");
        else if (imag == 7)
        strcpy (magname, "MagM");
        else if (imag == 8)
        strcpy (magname, "MagA");
    }
    else if (refcat==SKYBOT)
        strcpy (magname, "MagV");
    else
        strcpy (magname, "Mag");
    return;
}
/* from catutil.c */
/* Return number of decimal places in catalogued numbers, if known */
int
DS_CatNdec (refcat)

         int	refcat;		/* Catalog code */

{
    /* USNO A1.0, A2.0, SA1.0, or SA2.0 Catalogs */
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
            refcat == UAC  || refcat == UA1  || refcat == UA2)
        return (8);

    /* USNO B1.0 and USNO YB6 */
    else if (refcat == UB1 || refcat == YB6)
        return (7);

    /* GSC II */
    else if (refcat == GSC2)
        return (0);

    /* SDSS */
    else if (refcat == SDSS)
        return (0);

    /* SkyBot */
    else if (refcat == SKYBOT)
        return (0);

    /* 2MASS Point Source Catalog */
    else if (refcat == TMPSC || refcat == TMPSCE)
        return (6);

    /* 2MASS Extended Source Catalog */
    else if (refcat == TMXSC)
        return (6);

    /* 2MASS Point Source Catalog */
    else if (refcat == TMIDR2)
        return (7);

    /* USNO Plate Catalogs */
    else if (refcat == USNO)
        return (0);

    /* UCAC1 Catalog */
    else if (refcat == UCAC1)
        return (6);

    /* UCAC2 Catalog */
    else if (refcat == UCAC2)
        return (6);

    /* UCAC3 Catalog */
    else if (refcat == UCAC3)
        return (6);

    /* USNO UJ 1.0 Catalog */
    else if (refcat == UJC)
        return (7);

    /* HST Guide Star Catalog */
    else if (refcat == GSC || refcat == GSCACT)
        return (4);

    /* SAO, PPM, Hipparcos, or IRAS Point Source Catalogs (TDC binary format) */
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==BSC ||
                     refcat==HIP || refcat == SKY2K)
        return (0);

    /* Tycho, Tycho2, or ACT catalogs */
    else if (refcat == TYCHO || refcat == TYCHO2 ||
                     refcat == TYCHO2E || refcat == ACT)
        return (5);

    /* Starbase tab-separated, TDC binary, or TDC ASCII catalogs */
    else
        return (-1);
}
/* from catutil.c */

/* TMCID -- Return 1 if string is 2MASS ID, else 0 */

int
DS_tmcid (string, ra, dec)

         char	*string;	/* Character string to check */
         double	*ra;		/* Right ascension (returned) */
         double	*dec;		/* Declination (returned) */
{
    char *sdec;
    int idec, idm, ids, ira, irm, irs;

    /* Check first character */
    if (string[0] != 'J' && string[0] != 'j')
        return (0);

    /* Find declination sign */
    sdec = strsrch (string, "-");
    if (sdec == NULL)
        sdec = strsrch (string,"+");
    if (sdec == NULL)
        return (0);

    /* Parse right ascension */
    *sdec = (char) 0;
    ira = atoi (string+1);
    irs = ira % 10000;
    ira = ira / 10000;
    irm = ira % 100;
    ira = ira / 100;
    *ra = (double) ira + ((double) irm) / 60.0 + ((double) irs) / 360000.0;
    *ra = *ra * 15.0;

    /* Parse declination */
    idec = atoi (sdec+1);
    ids = idec % 1000;
    idec = idec / 1000;
    idm = idec % 100;
    idec = idec / 100;
    *dec = (double) idec + ((double) idm) / 60.0 + ((double) ids) / 36000.0;
    return (1);
}
/* from catutil.c */

/* Return number of decimal places in numeric string (-1 if not number) */

int
DS_StrNdec (string)

         char *string;	/* Numeric string */
{
    char *cdot;
    int lstr;

    if (notnum (string))
        return (-1);
    else {
        lstr = strlen (string);
        if ((cdot = strchr (string, '.')) == NULL)
        return (0);
        else
        return (lstr - (cdot - string));
    }
}
/* from catutil.c */

/* MOVEB -- Copy nbytes bytes from source+offs to dest+offd (any data type) */

void
DS_movebuff (source, dest, nbytes, offs, offd)

         char *source;	/* Pointer to source */
         char *dest;	/* Pointer to destination */
         int nbytes;	/* Number of bytes to move */
         int offs;	/* Offset in bytes in source from which to start copying */
         int offd;	/* Offset in bytes in destination to which to start copying */
{
    char *from, *last, *to;
    from = source + offs;
    to = dest + offd;
    last = from + nbytes;
    while (from < last) *(to++) = *(from++);
    return;
}
/* from catutil.c */


/* REFLIM-- Set limits in reference catalog coordinates given search coords */
void
DS_RefLim (cra, cdec, dra, ddec, sysc, sysr, eqc, eqr, epc, epr, secmarg,
                ramin, ramax, decmin, decmax, wrap, verbose)

         double	cra, cdec;	/* Center of search area  in degrees */
         double	dra, ddec;	/* Horizontal and vertical half-widths of area */
         int	sysc, sysr;	/* System of search, catalog coordinates */
         double	eqc, eqr;	/* Equinox of search, catalog coordinates in years */
         double	epc, epr;	/* Epoch of search, catalog coordinates in years */
         double	secmarg;	/* Margin in arcsec/century to catch moving stars */
         double	*ramin,*ramax;	/* Right ascension search limits in degrees (returned)*/
         double	*decmin,*decmax; /* Declination search limits in degrees (returned) */
         int	*wrap;		/* 1 if search passes through 0:00:00 RA */
         int	verbose;	/* 1 to print limits, else 0 */

{
    double ra, ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4;
    double dec, acdec, adec, adec1, adec2, dmarg, dist, dra1;
    int nrot;

    /* Deal with all or nearly all of the sky */
    if (ddec > 80.0 && dra > 150.0) {
        *ramin = 0.0;
        *ramax = 360.0;
        *decmin = -90.0;
        *decmax = 90.0;
        *wrap = 0;
        if (verbose)
        fprintf (stderr,"RefLim: RA: 0.0 - 360.0  Dec: -90.0 - 90.0\n");
        return;
    }

    /* Set declination limits for search */
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;

    /* dec1 is always the smallest declination */
    if (dec1 > dec2) {
        dec = dec1;
        dec1 = dec2;
        dec2 = dec;
    }
    dec3 = dec2;
    dec4 = dec1;

    /* Deal with south pole */
    if (dec1 < -90.0) {
        dec1 = 90.0 - (dec1 + 90.0);
        if (dec1 > dec2)
        dec2 = dec1;
        dec1 = -90.0;
        dra1 = 180.0;
    }

    /* Deal with north pole */
    if (dec2 > 90.0) {
        dec2 = 90.0 - (dec2 - 90.0);
        if (dec2 < dec1)
        dec1 = dec2;
        dec2 = 90.0;
        dra1 = 180.0;
    }

    /* Adjust width in right ascension to that at max absolute declination */
    adec1 = fabs (dec1);
    adec2 = fabs (dec2);
    if (adec1 > adec2)
        adec = adec1;
    else
        adec = adec2;
    acdec = fabs (cdec);
    if (adec < 90.0 && adec > acdec)
        dra1 = dra * (cos (degrad(acdec)) / cos (degrad(adec)));
    else if (adec == 90.0)
        dra1 = 180.0;

    /* Set right ascension limits for search */
    ra1 = cra - dra1;
    ra2 = cra + dra1;

    /* Keep right ascension limits between 0 and 360 degrees */
    if (ra1 < 0.0) {
        nrot = 1 - (int) (ra1 / 360.0);
        ra1 = ra1 + (360.0 * (double) nrot);
    }
    if (ra1 > 360.0) {
        nrot = (int) (ra1 / 360.0);
        ra1 = ra1 - (360.0 * (double) nrot);
    }
    if (ra2 < 0.0) {
        nrot = 1 - (int) (ra2 / 360.0);
        ra2 = ra2 + (360.0 * (double) nrot);
    }
    if (ra2 > 360.0) {
        nrot = (int) (ra2 / 360.0);
        ra2 = ra2 - (360.0 * (double) nrot);
    }

    if (ra1 > ra2)
        *wrap = 1;
    else
        *wrap = 0;

    ra3 = ra1;
    ra4 = ra2;

    /* Convert search corners to catalog coordinate system and equinox */
    ra = cra;
    dec = cdec;
    wcscon (sysc, sysr, eqc, eqr, &ra, &dec, epc);
    wcscon (sysc, sysr, eqc, eqr, &ra1, &dec1, epc);
    wcscon (sysc, sysr, eqc, eqr, &ra2, &dec2, epc);
    wcscon (sysc, sysr, eqc, eqr, &ra3, &dec3, epc);
    wcscon (sysc, sysr, eqc, eqr, &ra4, &dec4, epc);

    /* Find minimum and maximum right ascensions to search */
    *ramin = ra1;
    if (ra3 < *ramin)
        *ramin = ra3;
    *ramax = ra2;
    if (ra4 > *ramax)
        *ramax = ra4;

    /* Add margins to RA limits to get most stars which move */
    if (secmarg > 0.0 && epc != 0.0) {
        dmarg = (secmarg / 3600.0) * fabs (epc - epr);
        *ramin = *ramin - (dmarg * cos (degrad (cdec)));
        *ramax = *ramax + (dmarg * cos (degrad (cdec)));
    }
    else
        dmarg = 0.0;

    if (*wrap) {
        ra = *ramax;
        *ramax = *ramin;
        *ramin = ra;
    }

    /* Find minimum and maximum declinatons to search */
    *decmin = dec1;
    if (dec2 < *decmin)
        *decmin = dec2;
    if (dec3 < *decmin)
        *decmin = dec3;
    if (dec4 < *decmin)
        *decmin = dec4;
    *decmax = dec1;
    if (dec2 > *decmax)
        *decmax = dec2;
    if (dec3 > *decmax)
        *decmax = dec3;
    if (dec4 > *decmax)
        *decmax = dec4;

    /* Add margins to Dec limits to get most stars which move */
    if (dmarg > 0.0) {
        *decmin = *decmin - dmarg;
        *decmax = *decmax + dmarg;
    }

    /* Check for pole */
    dist = wcsdist (ra, dec, *ramax, *decmax);
    if (dec + dist > 90.0) {
        *ramin = 0.0;
        *ramax = 359.99999;
        *decmax = 90.0;
        *wrap = 0;
    }
    else if (dec - dist < -90.0) {
        *ramin = 0.0;
        *ramax = 359.99999;
        *decmin = -90.0;
        *wrap = 0;
    }


    /* Search zones which include the poles cover 360 degrees in RA */
    else if (*decmin < -90.0) {
        *decmin = -90.0;
        *ramin = 0.0;
        *ramax = 359.99999;
        *wrap = 0;
    }
    else if (*decmax > 90.0) {
        *decmax = 90.0;
        *ramin = 0.0;
        *ramax = 359.99999;
        *wrap = 0;
    }
    if (verbose) {
        char rstr1[16],rstr2[16],dstr1[16],dstr2[16];
        if (degout) {
        deg2str (rstr1, 16, *ramin, 6);
            deg2str (dstr1, 16, *decmin, 6);
        deg2str (rstr2, 16, *ramax, 6);
            deg2str (dstr2, 16, *decmax, 6);
        }
        else {
        ra2str (rstr1, 16, *ramin, 3);
            dec2str (dstr1, 16, *decmin, 2);
        ra2str (rstr2, 16, *ramax, 3);
            dec2str (dstr2, 16, *decmax, 2);
        }
        fprintf (stderr,"RefLim: RA: %s - %s  Dec: %s - %s",
                         rstr1,rstr2,dstr1,dstr2);
        if (*wrap)
        fprintf (stderr," wrap\n");
        else
        fprintf (stderr,"\n");
    }
    return;
}
/* from catutil.c */

/* SEARCHLIM-- Set RA and Dec limits for search given center and dimensions */
void
DS_SearchLim (cra, cdec, dra, ddec, syscoor, ra1, ra2, dec1, dec2, verbose)

         double	cra, cdec;	/* Center of search area  in degrees */
         double	dra, ddec;	/* Horizontal and vertical half-widths in degrees */
         int	syscoor;	/* Coordinate system */
         double	*ra1, *ra2;	/* Right ascension limits in degrees */
         double	*dec1, *dec2;	/* Declination limits in degrees */
         int	verbose;	/* 1 to print limits, else 0 */

{
    double dec;

    /* Set right ascension limits for search */
    *ra1 = cra - dra;
    *ra2 = cra + dra;

    /* Keep right ascension between 0 and 360 degrees */
    if (syscoor != WCS_XY) {
        if (*ra1 < 0.0)
        *ra1 = *ra1 + 360.0;
        if (*ra2 > 360.0)
        *ra2 = *ra2 - 360.0;
    }

    /* Set declination limits for search */
    *dec1 = cdec - ddec;
    *dec2 = cdec + ddec;

    /* dec1 is always the smallest declination */
    if (*dec1 > *dec2) {
        dec = *dec1;
        *dec1 = *dec2;
        *dec2 = dec;
    }

    /* Search zones which include the poles cover 360 degrees in RA */
    if (syscoor != WCS_XY) {
        if (*dec1 < -90.0) {
        *dec1 = -90.0;
        *ra1 = 0.0;
        *ra2 = 359.99999;
        }
        if (*dec2 > 90.0) {
        *dec2 = 90.0;
        *ra1 = 0.0;
        *ra2 = 359.99999;
        }
    }

    if (verbose) {
        char rstr1[16],rstr2[16],dstr1[16],dstr2[16];
        if (syscoor == WCS_XY) {
        num2str (rstr1, *ra1, 10, 5);
            num2str (dstr1, *dec1, 10, 5);
        num2str (rstr2, *ra2, 10, 5);
            num2str (dstr2, *dec2, 10, 5);
        }
        else if (degout) {
        deg2str (rstr1, 16, *ra1, 6);
            deg2str (dstr1, 16, *dec1, 6);
        deg2str (rstr2, 16, *ra2, 6);
            deg2str (dstr2, 16, *dec2, 6);
        }
        else {
        ra2str (rstr1, 16, *ra1, 3);
            dec2str (dstr1, 16, *dec1, 2);
        ra2str (rstr2, 16, *ra2, 3);
            dec2str (dstr2, 16, *dec2, 2);
        }
        fprintf (stderr,"SearchLim: RA: %s - %s  Dec: %s - %s\n",
                         rstr1,rstr2,dstr1,dstr2);
    }
    return;
}
/* from DS_ctgread.c (MODIFIED) */
/* CTGCLOSE -- Close ASCII catalog and free associated data structures */

void
DS_ctgclose (sc)

         struct	StarCat *sc;
{
  if (sc == NULL) {
    return;
  }
    if (sc->refcat == SCANCAT)
        scanclose (sc);
    else
        free (sc);

    sc = NULL;
    return;
}

/* Set parameter values from the command line as keyword=value
 * Return 1 if successful, else 0 */

int scatparm (char *parstring)
{
  char *parname;
  char *parvalue;
  char *parequal;
  char *temp;
  char *refcatn;
  int lcat, lrange;

  /* Check for request for command line help */
  if (!strncasecmp (parstring,"comhelp", 7)) {
    PrintUsage (NULL);
  }

  /* Check for scat version request */
  if (!strcasecmp (parstring, "version")) {
    PrintUsage ("version");
  }

  /* Separate parameter name and value */
  parname = parstring;
  if ((parequal = strchr (parname,'=')) == NULL)
    return (0);
  *parequal = (char) 0;
  parvalue = parequal + 1;

  /* Get closest source */
  if (!strcasecmp (parname, "closest")) {
    if (!strncasecmp (parvalue, "y", 1)) {
      catsort = SORT_DIST;
      nstars = 1;
      closest++;
    }
  }

  /* Set range of source numbers to print */
  else if (!strncasecmp (parname,"num",3)) {
    if (ranges) {
      temp = ranges;
      lrange = strlen(ranges) + strlen(parvalue) + 2;
      ranges = (char *) malloc (lrange);
      strcpy (ranges, temp);
      strcat (ranges, ",");
      strcat (ranges, parvalue);
      free (temp);
    }
    else {
      lrange = strlen (parvalue) + 2;
      ranges = (char *) malloc (lrange);
      if (strchr (parvalue,'.'))
                match = 1;
      strcpy (ranges, parvalue);
    }
  }

  /* Radius in arcseconds */
  else if (!strncasecmp (parname,"degree",6)) {
    if (strchr (parvalue, 'y'))
      degout0 = 1;
    else
      degout0 = 0;
  }

  /* Radius in arcseconds */
  else if (!strncasecmp (parname,"rad",3)) {
    if (strchr (parvalue,':'))
      rad0 = 3600.0 * str2dec (parvalue);
    else if (isnum (parvalue))
      rad0 = atof (parvalue);
    else {
      sprintf (voerror, "Search radius %s (seconds) is not a number", parvalue);
      return (-1);
    }
  }

  /* Radius in arcminutes */
  else if (!strncasecmp (parname,"mrad",4)) {
    if (strchr (parvalue,':'))	/* hours:minutes:seconds */
      rad0 = 3600.0 * str2dec (parvalue);
    else if (isnum (parvalue))
      rad0 = 60.0 * atof (parvalue);
    else {
      sprintf (voerror, "Search radius %s (minutes) is not a number", parvalue);
      return (-1);
    }
  }

  /* Inner Annulus radius in arcseconds */
  else if (!strncasecmp (parname,"inrad",5)) {
    if (strchr (parvalue,':'))	/* hours:minutes:seconds */
      rad1 = 3600.0 * str2dec (parvalue);
    else if (isnum (parvalue))
      rad1 = atof (parvalue);
    else {
      sprintf (voerror, "Search radius %s is not a number", parvalue);
      return (-1);
    }
  }

  /* Radius in degrees */
  else if (!strncasecmp (parname,"sr",2) ||
                     !strncasecmp (parname,"drad",4)) {
    if (strchr (parvalue,':'))
      rad0 = 3600.0 * str2dec (parvalue);
    else if (isnum (parvalue))
      rad0 = 3600.0 * atof (parvalue);
    else {
      sprintf (voerror, "Search radius %s (degrees) is not a number", parvalue);
      return (-1);
    }
    votab = 1;
    degout0 = 1;
  }

  /* Search center right ascension */
  else if (!strcasecmp (parname,"ra")) {
    if (!isnum (parvalue) && !strchr (parvalue,':')) {
      sprintf (voerror, "Right ascension %s is not a number", parvalue);
      return (-1);
    }
    ra0 = str2ra (parvalue);
  }

  /* Search center declination */
  else if (!strcasecmp (parname,"dec")) {
    if (!isnum (parvalue) && !strchr (parvalue,':')) {
      sprintf (voerror, "Declination %s is not a number", parvalue);
      return (-1);
    }
    dec0 = str2dec (parvalue);
  }

  /* Search center coordinate system */
  else if (!strncasecmp (parname,"sys",3)) {
    syscoor = wcscsys (parvalue);
    eqcoor = wcsceq (parvalue);
  }

  /* Output coordinate system */
  else if (!strcasecmp (parname, "outsys")) {

    /* B1950 (FK4) coordinates */
    if (!strcasecmp (parvalue, "B1950") ||
                !strcasecmp (parvalue, "FK4")) {
      sysout0 = WCS_B1950;
      eqout = 1950.0;
    }

    /* J2000 (FK5) coordinates */
    else if (!strcasecmp (parvalue, "J2000") ||
                         !strcasecmp (parvalue, "FK5")) {
      sysout0 = WCS_J2000;
      eqout = 2000.0;
    }

    /* Galactic coordinates */
    else if (!strncasecmp (parvalue, "GAL", 3))
      sysout0 = WCS_GALACTIC;

    /* Ecliptic coordinates */
    else if (!strncasecmp (parvalue, "ECL", 3))
      sysout0 = WCS_ECLIPTIC;
  }

  /* Set reference catalog */
  else if (!strcasecmp (parname, "catalog")) {
    lcat = strlen (parvalue) + 2;
    refcatn = (char *) malloc (lcat);
    strcpy (refcatn, parvalue);
    refcatname[ncat] = refcatn;
    ncat = ncat + 1;
  }

  /* Set output coordinate epoch */
  else if (!strcasecmp (parname, "epoch"))
    epoch0 = fd2ep (parvalue);

  /* Output equinox in years */
  else if (!strcasecmp (parname, "equinox"))
    eqout = fd2ep (parvalue);

  /* Output in degrees instead of sexagesimal */
  else if (!strcasecmp (parname, "cformat")) {
    if (!strncasecmp (parvalue, "deg", 3))
      degout0 = 1;
    else if (!strncasecmp (parvalue, "rad", 3))
      degout0 = 2;
    else
      degout0 = 0;
  }

  /* Number of decimal places in output positions */
  else if (!strcasecmp (parname, "ndec")) {
    if (isnum (parvalue)) {
      if (degout0) {
                nddeg = atoi (parvalue);
                if (nddeg < 0 || nddeg > 10)
                    nddeg = 7;
      }
      else {
                ndra = atoi (parvalue);
                if (ndra < 0 || ndra > 10)
                    ndra = 3;
                nddec = ndra - 1;
      }
    }
  }

  /* Minimum proper motion quality for USNO-B1.0 catalog */
  else if (!strcasecmp (parname, "minpmq")) {
    if (isnum (parvalue)) {
      minid = atoi (parvalue);
      setminpmqual (minid);
    }
  }

  /* Format for epoch of catalog source position */
  else if (!strcasecmp (parname, "dateform")) {
    if (parvalue[0] == 'j')
      setdateform (EP_JD);
    else if (parvalue[0] == 'm')
      setdateform (EP_MJD);
    else if (parvalue[0] == 'f')
      setdateform (EP_FD);
    else if (parvalue[0] == 'i')
      setdateform (EP_ISO);
    else
      setdateform (EP_EP);
  }

  /* Minimum number of plate ID's for USNO-B1.0 catalog */
  else if (!strcasecmp (parname, "minid")) {
    if (isnum (parvalue)) {
      minid = atoi (parvalue);
      setminid (minid);
    }
  }

  /* Output in VOTable XML instead of tab-separated table */
  else if (!strcasecmp (parname, "format")) {
    if (!strncasecmp (parvalue, "vot", 3)) {
      votab = 1;
      degout0 = 1;
    }
    else
      votab = 0;
  }

  /* Print center and closest star on one line */
  else if (!strcasecmp (parname, "oneline")) {
    if (parvalue[0] == 'y' || parvalue[0] == 'Y') {
      oneline++;
      catsort = SORT_DIST;
      closest++;
      nstars = 1;
    }
  }

  /* Magnitude limit */
  else if (!strncasecmp (parname,"mag",3)) {
    maglim2 = atof (parvalue);
    if (MAGLIM1 == MAGLIM2)
      maglim1 = -2.0;
  }
  else if (!strncasecmp (parname,"max",3))
    maglim2 = atof (parvalue);
  else if (!strncasecmp (parname,"min",3))
    maglim1 = atof (parvalue);

  /* Number of brightest stars to read */
  else if (!strncasecmp (parname,"nstar",5))
    nstars = atoi (parvalue);

  /* Object name */
  else if (!strcmp (parname, "object") || !strcmp (parname, "OBJECT")) {
    lcat = strlen (parvalue) + 2;
    objname = (char *) malloc (lcat);
    strcpy (objname, parvalue);
  }

  /* Magnitude by which to sort */
  else if (!strcasecmp (parname, "sortmag")) {
    if (isnum (parvalue+1))
      sortmag = atoi (parvalue);
  }

  /* Output sorting */
  else if (!strcasecmp (parname, "sort")) {

    /* Sort by distance from center */
    if (!strncasecmp (parvalue,"di",2))
      catsort = SORT_DIST;

    /* Sort by RA */
    else if (!strncasecmp (parvalue,"r",1))
      catsort = SORT_RA;

    /* Sort by Dec */
    else if (!strncasecmp (parvalue,"de",2))
      catsort = SORT_DEC;

    /* Sort by ID */
    else if (!strncasecmp (parvalue,"id",2))
      catsort = SORT_ID;

    /* Sort by magnitude */
    else if (!strncasecmp (parvalue,"m",1)) {
      catsort = SORT_MAG;
      if (strlen (parvalue) > 1) {
                if (isnum (parvalue+1))
                    sortmag = atoi (parvalue+1);
      }
    }

    /* No sort */
    else if (!strncasecmp (parvalue,"n",1))
      catsort = SORT_NONE;
    else
      catsort = SORT_NONE;
  }

  /* Search box half-width in RA */
  else if (!strcasecmp (parname,"dra")) {
    if (strchr (parvalue,':'))
      dra0 = 3600.0 * str2ra (parvalue);
    else
      dra0 = atof (parvalue);
    if (ddec0 <= 0.0)
      ddec0 = dra0;
  }

  /* Search box half-height in Dec */
  else if (!strcasecmp (parname,"ddec")) {
    if (strchr (parvalue,':'))
      ddec0 = 3600.0 * str2dec (parvalue);
    else
      ddec0 = atof (parvalue);
    if (dra0 <= 0.0)
      dra0 = ddec0;
  }

  /* Catalog to be searched */
  else if (!strncasecmp (parname,"cat",3)) {
    lcat = strlen (parvalue) + 2;
    refcatn = (char *) malloc (lcat);
    strcpy (refcatn, parvalue);
    refcatname[ncat] = refcatn;
    ncat = ncat + 1;
  }
  else {
    *parequal = '=';
    return (1);
  }
  *parequal = '=';
  return (0);
}



#define TABMAX 64

static int
ListCat (ranges, eqout)

     char	*ranges;	/* String with range of catalog numbers to list */
     double	eqout;		/* Equinox for output coordinates */

{
  double cra = 0.0;	/* Search center long or RA in degrees */
  double cdec = 0.0;	/* Search center lat or Dec in degrees */
  double crao, cdeco;	/* Output center long/lat or RA/Dec in degrees */
  double epout = 0.0;
  int ng = 0;		/* Number of catalog stars */
  int ns = 0;		/* Number of brightest catalog stars actually used */
  struct Range *range = NULL; /* Range of catalog numbers to list */
  int i, is, j, ngmax, nc;
  double das, dds, drs;
  int degout;
  double maxnum;
  FILE *fd = NULL;
  char rastr[32], decstr[32];	/* coordinate strings */
  char numstr[80];	/* Catalog number */
  char *catalog;
  double drad = 0.0;
  double dradi = 0.0;
  double dra, ddec, mag1, mag2;
  double gdist, da, dd, dec, gdmax;
  double epoch = 0;
  double date, time;
  double era = 0., edec = 0., epmr = 0., epmd = 0.;
  int nim = 0, nct = 0;
  int nlog;
  int magsort;
  int typecol = 0;
  int band = 0;
  int imag, nmagr = 0;
  int sysout = 0;
  char headline[160];
  char filename[80];
  char temp[80];
  char *dtemp;
  char isp[4];
  int ngsc = 0;
  int smag;
  int nid;
  int icat, nndec, nnfld = 0, nsfld;
  int gcset = 0;
  int ndist;
  int distsort;
  int lrv;
  int lobj;
  char tstr[32];
  double flux;
  double pra = 0.0;
  double pdec = 0.0;
  char magname[16];
  void ep2dt();
  void PrintNum();
  int LenNum();
  int printdist = 1; /* Show distance from center */

  /* Drop out if no catalog is specified */
  if (ncat < 1) {
    fprintf (stderr, "No catalog specified\n");
    PrintUsage(NULL);
    exit (-1);
  }

  /* Allocate space for returned catalog information */
  if (ranges != NULL) {
    int nfdef = 9;

    /* Allocate and fill list of numbers to read */
    range = DS_RangeInit (ranges, nfdef);
    ngmax = DS_rgetn (range) * 4;
  }
  else if (nstars != 0)
    ngmax = nstars;
  else
    ngmax = MAXCAT;

  if (ngmax > nalloc) {

    /* Free currently allocated buffers if more entries are needed */
    if (nalloc > 0) {
      if (gm) {
                for (imag = 0; imag < nmagmax; imag++)
                    if (gm[imag]) free ((char *) gm[imag]);
                free ((char *)gm);
      }
      if (gra) free ((char *)gra);
      if (gdec) free ((char *)gdec);
      if (gpra) free ((char *)gpra);
      if (gpdec) free ((char *)gpdec);
      if (gnum) free ((char *)gnum);
      if (gc) free ((char *)gc);
      if (gx) free ((char *)gx);
      if (gy) free ((char *)gy);
      if (gobj) {
                for (is = 1; is < nalloc; is++) {
                    if (gobj[is] != NULL)
                        free ((char *)gobj[is]);
                }
                free ((char *)gobj);
      }
    }
    gm = NULL;
    gra = NULL;
    gdec = NULL;
    gpra = NULL;
    gpdec = NULL;
    gnum = NULL;
    gc = NULL;
    gx = NULL;
    gy = NULL;
    gobj = NULL;

    if (!(gnum = (double *) calloc (ngmax, sizeof(double))))
      fprintf (stderr, "Could not calloc %lu bytes for gnum\n",
                             ngmax*sizeof(double));
    if (!(gra = (double *) calloc (ngmax, sizeof(double))))
      fprintf (stderr, "Could not calloc %lu bytes for gra\n",
                             ngmax*sizeof(double));
    if (!(gdec = (double *) calloc (ngmax, sizeof(double))))
      fprintf (stderr, "Could not calloc %lu bytes for gdec\n",
                             ngmax*sizeof(double));
    if (!(gm = (double **) calloc (nmagmax, sizeof(double *))))
      fprintf (stderr, "Could not calloc %lu bytes for gm\n",
                             nmagmax*sizeof(double *));
    else {
      for (imag = 0; imag < nmagmax; imag++) {
                if (!(gm[imag] = (double *) calloc (ngmax, sizeof(double))))
                    fprintf (stderr, "Could not calloc %lu bytes for gm\n",
                                     ngmax*sizeof(double));
      }
    }
    if (!(gc = (int *) calloc (ngmax, sizeof(int))))
      fprintf (stderr, "Could not calloc %lu bytes for gc\n",
                             ngmax*sizeof(int));
    if (!(gx = (double *) calloc (ngmax, sizeof(double))))
      fprintf (stderr, "Could not calloc %lu bytes for gx\n",
                             ngmax*sizeof(double));
    if (!(gy = (double *) calloc (ngmax, sizeof(double))))
      fprintf (stderr, "Could not calloc %lu bytes for gy\n",
                             ngmax*sizeof(double));
    if (!(gobj = (char **) calloc (ngmax, sizeof(char *))))
      fprintf (stderr, "Could not calloc %lu bytes for obj\n",
                             ngmax*sizeof(char *));
    else {
      for (i = 0; i < ngmax; i++)
                gobj[i] = NULL;
    }
    if (!(gpra = (double *) calloc (ngmax, sizeof(double))))
      fprintf (stderr, "Could not calloc %lu bytes for gpra\n",
                             ngmax*sizeof(double));
    if (!(gpdec = (double *) calloc (ngmax, sizeof(double))))
      fprintf (stderr, "Could not calloc %lu bytes for gpdec\n",
                             ngmax*sizeof(double));
    if (gnum==NULL || gra==NULL || gdec==NULL || gm==NULL || gc==NULL ||
                gx==NULL || gy==NULL || gobj==NULL || gpra == NULL || gpdec == NULL){
      if (gm) {
                for (imag = 0; imag < nmagmax; imag++)
                    if (gm[imag]) free ((char *) gm[imag]);
                free ((char *)gm);
                gm = NULL;
      }
      if (gra) free ((char *)gra);
      gra = NULL;
      if (gdec) free ((char *)gdec);
      gdec = NULL;
      if (gpra) free ((char *)gpra);
      gpra = NULL;
      if (gpdec) free ((char *)gpdec);
      gpdec = NULL;
      if (gnum) free ((char *)gnum);
      gnum = NULL;
      if (gc) free ((char *)gc);
      gc = NULL;
      if (gx) free ((char *)gx);
      gx = NULL;
      if (gy) free ((char *)gy);
      gy = NULL;
      if (gobj) {
                for (is = 1; is < nalloc; is++) {
                    if (gobj[is] != NULL)
                        free ((char *)gobj[is]);
                }
                free ((char *)gobj);
      }
      gobj = NULL;
      nalloc = 0;
      return (0);
    }

    /* Initialize catalog entry values */
    for (i = 0; i < ngmax; i++) {
      for (imag = 0; imag < nmagmax; imag++)
                gm[imag][i] = 99.0;
      gra[i] = 0.0;
      gdec[i] = 0.0;
      gpra[i] = 0.0;
      gpdec[i] = 0.0;
      gnum[i] = 0.0;
      gc[i] = 0;
      gx[i] = 0.0;
      gy[i] = 0.0;
    }
    nalloc = ngmax;
  }

  /* Start of per catalog loop */
  for (icat = 0; icat < ncat; icat++) {
    if (ncat > 1)
      nohead = 1;
    nndec = 0;
    isp[2] = (char) 0;
    isp[3] = (char) 0;

    /* Skip this catalog if no name is given */
    if (refcatname[icat] == NULL || strlen (refcatname[icat]) == 0) {
      fprintf (stderr, "Catalog %d not specified\n", icat);
      continue;
    }

    if (printprog && notprinted) {
      if (closest)
                printf ("\n%s %s  Find closest star\n", progname, RevMsg);
      else
                printf ("\n%s %s\n", progname, RevMsg);
    }

    /* Figure out which catalog we are searching */
    if (ncat > 1 || refcat == 0) {
      if (!(refcat = DS_RefCat (refcatname[icat],title,&sysref,&eqref,
                                                         &epref,&mprop,&nmag))) {
                fprintf (stderr,"ListCat: Catalog '%s' is missing\n", refcatname[icat]);
                return (0);
      }
    }

    if (webdump)
      nlog = -1;
    else if (debug)
      nlog = 1;
    else if (verbose) {
      if (refcat==UAC || refcat==UA1 || refcat==UA2 || refcat==UB1 ||
                    refcat==USAC || refcat==USA1 || refcat==USA2 || refcat==GSC ||
                    refcat==GSCACT || refcat==TMPSC || refcat==TMPSCE ||
                    refcat==TMIDR2 || refcat==TMXSC || refcat == YB6)
                nlog = 1000;
      else
                nlog = 100;
    }
    else
      nlog = 0;

    /* If more magnitudes are needed, allocate space for them */
    if (nmag > nmagmax) {
      if (gm) {
                for (imag = 0; imag < nmagmax; imag++)
                    free ((char *) gm[imag]);
                free ((char *)gm);
                gm = NULL;
      }
      nmagmax = nmag;
      if (!(gm = (double **) calloc (nmagmax, sizeof(double *))))
                fprintf (stderr, "Could not calloc %lu bytes for gm\n",
                                 nmagmax*sizeof(double *));
      else {
                for (imag = 0; imag < nmagmax; imag++) {
                    if (!(gm[imag] = (double *) calloc (ngmax, sizeof(double))))
                        fprintf (stderr, "Could not calloc %lu bytes for gm\n",
                                         ngmax*sizeof(double));
                }
      }
    }

    /* Set output coordinate system from command line or catalog */
    if (sysout0)
      sysout = sysout0;
    else if (srch!= NULL && srch->epoch != 0.0)
      sysout = srch->coorsys;
    if (!sysout)
      sysout = sysref;
    if (!sysout)
      sysout = WCS_J2000;

    /* Set equinox from command line, search catalog, or searched catalog */
    if (eqout == 0.0) {
      if (srch!= NULL && srch->equinox != 0.0)
                eqout = srch->equinox;
      else if (sysout0 == WCS_J2000)
                eqout = 2000.0;
      else if (sysout0 == WCS_B1950)
                eqout = 1950.0;
      else if (eqcoor != 0.0)
                eqout = eqcoor;
      else
                eqout = epref;
      if (eqout == 0.0)
                eqout = 2000.0;
    }

    /* Set epoch from command line, search catalog, or searched catalog */
    epout = epoch0;
    if (epout == 0.0) {
      if (srch!= NULL && srch->epoch != 0.0)
                epout = srch->epoch;
      else if (sysout0 == WCS_J2000)
                epout = 2000.0;
      else if (sysout0 == WCS_B1950)
                epout = 1950.0;
      else if (mprop != 1)
                epout = epref;
      else if (eqout != 0.0)
                epout = eqout;
      else {
                if (sysout == WCS_B1950)
                    epout = 1950.0;
                else
                    epout = 2000.0;
      }
    }

    /* Set degree flag for output */
    if (sysout == WCS_ECLIPTIC || sysout == WCS_GALACTIC)
      degout = 1;
    else
      degout = degout0;
    DS_setlimdeg (degout);

    /* Find stars specified by number */
    if ((ranges != NULL) && (refcat != SCANCAT))
      {
      }

    /* Find catalog entries specified by date */
    else
      if ((dec0 < -90.0) &&
                    (epoch1 > 0.0) &&
                    (epoch2 > 0.0) &&
                    (refcat != SCANCAT))

                {

                }

    /* Find stars specified by location */
      else {
                /* Default to sort stars by magnitude */
                if (catsort == SORT_UNSET)
                    catsort = SORT_MAG;

                /* Set search radius if finding closest star */
                if (rad0 == 0.0 && dra0 == 0.0) {
                    if (closest)
                        rad0 = DS_CatRad (refcat);
                    else
                        rad0 = 10.0;
                }

                /* Set limits from defaults and command line information */
                if (GetArea (verbose,syscoor,eqcoor,sysout,eqout,epout,
                                         &cra,&cdec,&dra,&ddec,&drad,&dradi,&crao,&cdeco))
                    return (0);

                if (srch != NULL) {
                    if (srchcat->stnum <= 0 && strlen (srch->objname) > 0) {
                        if (objname == NULL)
                            objname = (char *) malloc (32);
                        strcpy (objname, srch->objname);
                    }
                    else {
                        if (objname == NULL)
                            objname = (char *) calloc (1,32);
                        nsfld = DS_CatNumLen (TXTCAT, srch->num, srchcat->nndec);
                        DS_CatNum (TXTCAT, nsfld, srchcat->nndec, srch->num, objname);
                    }
                }
                nnfld = DS_CatNumLen (refcat, 0.0, nndec);

                /* Print search center and size in input and output coordinates */
                if (verbose || (printhead && !oneline)) {
                    if (sysout != syscoor || eqcoor != eqout)
                        SearchHead (icat,syscoor,eqcoor,epout,
                                                cra,cdec,dra,ddec,drad,dradi,nnfld,degout);
                    SearchHead (icat,sysout,eqout,epout,
                                            crao,cdeco,dra,ddec,drad,dradi,nnfld,degout);
                    if (!closest) {
                        if (sysref != syscoor && sysref != sysout) {
                            double cra2 = cra;
                            double cdec2 = cdec;
                            wcscon (syscoor, sysref, 0.0, 0.0, &cra2, &cdec2, epout);
                            SearchHead (icat,sysref,eqref, epref,cra2,cdec2,
                                                    dra,ddec,drad,dradi,nnfld,degout);
                        }
                    }
                }

                /* Set the magnitude limits for the catalog search */
                if (maglim2 == 0.0) {
                    mag1 = 0.0;
                    mag2 = 0.0;
                }
                else {
                    mag1 = maglim1;
                    mag2 = maglim2;
                }

                /* Find the nearby reference stars, in ra/dec */
                if (catsort == SORT_DIST)
                    distsort = 1;
                else
                    distsort = 0;
                if (sortmag > 9) {
                    sortmag = DS_CatMagNum (sortmag, refcat);
                }
                setSearchParameters(epoch1,epoch2);

                ng = DS_ctgread (refcatname[icat], refcat, distsort, crao, cdeco,
                                            dra,ddec,drad,dradi,sysout,eqout,epout,mag1,mag2,
                                            sortmag,ngmax,&starcat[icat],
                                            gnum,gra,gdec,gpra,gpdec,gm,gc,gobj,nlog);
                if ((starcat[icat] != NULL) &&
                        (starcat[icat]->refcat == SCANCAT) &&
                        ((crao <= -91.0) || (cdeco <= -91.0))) {
                    printdist = 0;
                }

                if (webtable == 0) {

                    if (ngmax < 1)
                        return (ng);

                    if ((verbose || printhead) && ncat == 1 && ng < 1) {
                        fprintf (stderr, "No stars found in %s\n",refcatname[icat]);
                        return (0);
                    }

                    if (gobj[0] == NULL)
                        gobj1 = NULL;
                    else
                        gobj1 = gobj;
                    if (ng > ngmax)
                        ns = ngmax;
                    else
                        ns = ng;

                    /* Set flag if any proper motions are non-zero */
                    if (mprop == 1 && !oneline) {
                        mprop = 0;
                        for (i = 0; i < ns; i++) {
                            if (gpra[i] != 0.0 || gpdec[i] != 0.0) {
                                mprop = 1;
                                break;
                            }
                        }
                    }

                    /* Convert coordinates to output coordinate system */
                    if (syscoor != sysout) {
                        if (mprop == 1) {
                            for (i = 0; i < ns; i++) {
                                wcsconp (syscoor,sysout,eqout,eqout,epout,epout,
                                                 &gra[i],&gdec[i],&gpra[i],&gpdec[i]);
                            }
                        }
                        else {
                            for (i = 0; i < ns; i++) {
                                wcscon (syscoor,sysout,eqout,eqout,
                                                &gra[i],&gdec[i],epout);
                            }
                        }
                    }

                    /* Set flag if radial velocity is included in catalog entry */
                    if (mprop == 0 && starcat[icat] != NULL && starcat[icat]->entrv > 0)
                        mprop = 2;

                    /* Check to see if epoch is contained in entries */
                    if (starcat[icat] != NULL && starcat[icat]->nepoch != 0)
                        printepoch = 1;
                    else
                        printepoch = 0;

                    /* Find largest catalog number to be printed */
                    if (starcat[icat] != NULL && starcat[icat]->nnfld != 0)
                        nnfld = starcat[icat]->nnfld;
                    else if (refcat == SKYBOT) {
                        nnfld = 6;
                        for (i = 0; i < ng; i++) {
                            lobj = strlen (gobj[i]);
                            if (lobj > nnfld)
                                nnfld = lobj;
                        }
                    }
                    else {
                        maxnum = 0.0;
                        for (i = 0; i < ns; i++ ) {
                            if (gnum[i] > maxnum)
                                maxnum = gnum[i];
                        }
                        nnfld = DS_CatNumLen (refcat, maxnum, nndec);
                    }

                    /* Check to see whether gc is set at all */
                    gcset = 0;
                    for (i = 0; i < ns; i++ ) {
                        if (gc[i] != 0) {
                            gcset = 1;
                            break;
                        }
                    }

                    /* Set flag for plate, class, or type column */
                    if (refcat == BINCAT || refcat == SAO  || refcat == PPM ||
                            refcat == BSC || refcat == SDSS || refcat == SKY2K)
                        typecol = 1;
                    else if ((refcat == GSC || refcat == GSCACT) && classd < -1)
                        typecol = 3;
                    else if (refcat == GSC || refcat == GSCACT ||
                                     refcat == UJC ||  refcat == UB1 ||
                                     refcat == USAC || refcat == USA1   || refcat == USA2 ||
                                     refcat == UAC  || refcat == UA1    || refcat == UA2 ||
                                     refcat == BSC  || (refcat == TABCAT&&gcset))
                        typecol = 2;
                    else if (starcat[icat] != NULL && starcat[icat]->sptype > 0)
                        typecol = 1;
                    else
                        typecol = 0;

                    /* Set number of magnitudes (n-1 if radial velocity present) */
                    if (mprop == 2)
                        nmagr = nmag - 1;
                    else
                        nmagr = nmag;
                    if (printepoch)
                        nmagr = nmagr - 1;

                    if (printdist)
                        /* Compute distance from star to search center */
                        for (i = 0; i < ns; i++ ) {
                            gx[i] = wcsdist (crao, cdeco, gra[i], gdec[i]);
                            gy[i] = 1.0;
                        }

                    /* Sort catalogued objects, if requested */
                    if (ns > 1) {

                        /* Sort found catalog objects from closest to furthest */
                        if (catsort == SORT_DIST)
                            XSortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,ns,
                                                    nmagmax);

                        /* Sort found catalog objects by right ascension */
                        else if (catsort == SORT_RA)
                            RASortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,ns,
                                                     nmagmax);

                        /* Sort found catalog objects by declination */
                        else if (catsort == SORT_DEC)
                            DecSortStars(gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,ns,
                                                     nmagmax);

                        /* Sort found catalog objects from brightest to faintest */
                        else if (catsort == SORT_MAG)
                            MagSortStars(gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,ns,
                                                     nmagmax,sortmag);

                        /* Sort found catalog objects by ID number */
                        else if (catsort == SORT_ID)
                            IDSortStars(gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,
                                                    ns,nmagmax);
                    }

                    /* Print one line with search center and found star */
                    if (oneline) {
                        if (ns > 0) {
                            das = dra0 * 3600.0;
                            dds = ddec0 * 3600.0;
                            drs = rad0;
                            if (drs < 0.0)
                                drs = -drs;
                            else if (drs == 0.0)
                                drs = sqrt (das*das + dds*dds);
                            if (das <= 0.0) {
                                das = drs;
                                dds = drs;
                            }
                            if (srchcat != NULL) {
                                smag = srchcat->nmag;
                                if (srchcat->entrv > 0)
                                    smag = smag - 1;
                                if (srchcat->nepoch)
                                    smag = smag - 1;
                            }
                            else
                                smag = 1;

                            if (nohead && tabout) {

                                /* Write tab table heading */
                                catalog = DS_CatName (refcat, refcatname[icat]);
                                printf ("catalog	%s\n", catalog);
                                if (listfile != NULL)
                                    printf ("search	%s\n", listfile);
                                if (sysout == WCS_GALACTIC)
                                    printf ("radecsys	galactic\n");
                                else if (sysout == WCS_ECLIPTIC)
                                    printf ("radecsys	ecliptic\n");
                                else if (sysout == WCS_B1950)
                                    printf ("radecsys	fk4\n");
                                else
                                    printf ("radecsys	fk5\n");
                                printf ("equinox	%.4f\n", eqout);
                                if (!printepoch)
                                    printf ("epoch	%.4f\n", epout);
                                if (minid != 0)		/* Min number of plate IDs for USNO-B1.0 */
                                    printf ("minid	%d\n", minid);
                                if (minpmqual > 0)	/* Min proper motion quality for USNO-B1.0 */
                                    printf ("minpmq	%d\n", minpmqual);
                                if (dra0 > 0.0 || rad0 < 0) {
                                    if (syscoor == WCS_GALACTIC) {
                                        printf ("glonsec	%.2f\n", das);
                                        printf ("glatsec	%.2f\n", dds);
                                    }
                                    else if (syscoor == WCS_ECLIPTIC) {
                                        printf ("elonsec	%.2f\n", das);
                                        printf ("elatsec	%.2f\n", dds);
                                    }
                                    else {
                                        printf ("drasec	%.2f\n", das);
                                        printf ("ddecsec	%.2f\n", dds);
                                    }
                                }
                                else if (rad0 > 0)
                                    printf ("radsec	%.2f\n", drs);
                                if (mprop==1 || (srchcat != NULL && srchcat->mprop==1)) {
                                    if (degout)
                                        printf ("pmunit	mas/yr\n");
                                    else {
                                        printf ("rpmunit	mas/yr\n");
                                        printf ("dpmunit	mas/yr\n");
                                    }
                                }
                                printf ("program	%s %s\n", progname, RevMsg);

                                /* Write column headings */
                                if (srch != NULL) {
                                    if (srchcat->keyid[0] > 0) {
                                        printf ("%s", srchcat->keyid);
                                        nc = strlen (srchcat->keyid);
                                    }
                                    else {
                                        printf ("srch_id");
                                        nc = 7;
                                    }
                                    if (srchcat->nnfld > nc) {
                                        for (i = nc; i < srchcat->nnfld; i++)
                                            printf (" ");
                                    }
                                    printf ("	");
                                }
                                printf ("srch_ra     	srch_dec    	");
                                if (srchcat != NULL) {
                                    if (smag > 0) {
                                        for (imag = 0; imag < smag; imag++) {
                                            if (strlen (srchcat->keymag[imag]) >0)
                                                printf ("s_%s	", srchcat->keymag[imag]);
                                            else if (smag > 1)
                                                printf ("s_mag%d	", imag+1);
                                            else
                                                printf ("s_mag	");
                                        }
                                    }
                                    if (srchcat->nepoch)
                                        printf ("epoch    	");
                                    if (srchcat->sptype) {
                                        if (strlen (srch->isp) > 2) {
                                            padspt = 1;
                                            if (strlen (srchcat->keytype) >0)
                                                printf ("s_%s 	",srchcat->keytype);
                                            else
                                                printf ("s_type  	");
                                        }
                                        else {
                                            padspt = 0;
                                            if (strlen (srchcat->keytype) >0)
                                                printf ("s_%s	",srchcat->keytype);
                                            else
                                                printf ("s_type	");
                                        }
                                    }
                                    if (srchcat->entrv > 0) {
                                        if ((lrv = strlen (srchcat->keyrv)) >0) {
                                            sprintf (tstr, "%s	            ",
                                                             srchcat->keyrv);
                                            tstr[9] = (char)0;
                                            printf ("%s	", tstr);
                                        }
                                        else
                                            printf ("s_rv    	");
                                    }
                                    if (srchcat->mprop > 0)
                                        printf ("s_pra 	s_pdec	");

                                }
                                if (refcat == TABCAT && starcat[icat]->keyid[0] >0) {
                                    strcpy (headline, starcat[icat]->keyid);
                                    strcat (headline, "                ");
                                }
                                else
                                    DS_CatID (headline, refcat);
                                headline[nnfld] = (char) 0;
                                printf ("%s", headline);
                                printf ("	ra          	dec        	");
                                if (refcat == GSC2)
                                    printf ("magf 	magj 	magn 	magv	");
                                else if (refcat == HIP)
                                    printf ("magb 	magv 	parlx	parer	");
                                else if (refcat == SKYBOT)
                                    printf ("magv 	gdist 	hdist	");
                                else if (refcat == IRAS)
                                    printf ("f10m 	f25m 	f60m 	f100m	");
                                else if (refcat == TMPSC || refcat == TMIDR2)
                                    printf ("magj   	magh    	magk   		");
                                else if (refcat == TMPSCE)
                                    printf ("magje 	maghe  	magke  	");
                                else if (refcat == TMXSC)
                                    printf ("magj  	magh   	magk   	size  	");
                                else if (refcat == UB1)
                                    printf ("magb1	magr1	magb2	magr2	magn 	");
                                else if (refcat == YB6)
                                    printf ("magb 	magr 	magj 	magh 	magk 	");
                                else if (refcat == SDSS)
                                    printf ("magu 	magg 	magr 	magi 	magz 	");
                                else if (refcat == UCAC2)
                                    printf ("raerr	decerr	magj 	magh 	magk 	magc 	");
                                else if (refcat == UCAC3)
                                    printf ("raerr	decerr	magb 	magr 	magi 	magj 	magh 	magk 	magm 	maga 	");
                                else if (nmagr > 0) {
                                    for (imag = 0; imag < nmagr; imag++) {
                                        if (starcat[icat] != NULL &&
                                                strlen (starcat[icat]->keymag[imag]) >0)
                                            printf ("%s	", starcat[icat]->keymag[imag]);
                                        else if (nmagr > 1)
                                            printf ("mag%d  ", imag+1);
                                        else
                                            printf ("mag    ");
                                    }
                                }
                                if (typecol == 1)
                                    printf ("spt   	");
                                if (mprop == 1)
                                    printf ("pmra  	pmdec 	");
                                if (refcat == UCAC2 || refcat == UCAC3)
                                    printf ("epmra	epmdec	ni	nc	");
                                if (refcat == UB1)
                                    printf ("pm	ni	sg	");
                                if (refcat == GSC2)
                                    printf ("class	");
                                if ((starcat[icat]!=NULL && starcat[icat]->entrv>0) &&
                                        mprop == 2)
                                    printf ("velocity	");
                                printf ("n 	dra");
                                for (i = 3; i < LenNum(das,2); i++)
                                    printf (" ");
                                printf ("	ddec");
                                for (i = 4; i < LenNum(dds,2); i++)
                                    printf (" ");
                                printf ("	drad");
                                for (i = 4; i < LenNum(drs,2); i++)
                                    printf (" ");
                                printf ("\n");
                                if (srch != NULL) {
                                    printf ("-------");
                                    if (srchcat->nnfld > 7) {
                                        for (i = 8; i < srchcat->nnfld; i++)
                                            printf ("-");
                                    }
                                    printf ("	");
                                }
                                printf ("------------	------------	");
                                if (srchcat != NULL) {
                                    for (imag = 0; imag < smag; imag++)
                                        printf ("-----	");
                                    if (srchcat->nepoch)
                                        printf ("---------	");
                                    if (srchcat->sptype != 0) {
                                        if (padspt)
                                            printf ("-------	");
                                        else
                                            printf ("----	");
                                    }
                                    if (srchcat->entrv > 0)
                                        printf ("---------	");
                                    if (srchcat->mprop == 1)
                                        printf ("------	------	");
                                }
                                strcpy (headline,"----------------------");
                                headline[nnfld] = (char) 0;
                                printf ("%s", headline);
                                printf ("	------------	------------	");
                                if (refcat == UCAC2 || refcat == UCAC3)
                                    printf ("-----	-----	");
                                if (refcat == GSC2)
                                    printf ("-----	-----	-----	-----	");
                                else if (refcat == HIP || refcat == IRAS)
                                    printf ("-----	-----	-----	-----	");
                                else if (refcat == TMPSC || refcat == TMIDR2)
                                    printf ("-------	-------	-------	");
                                else if (refcat == TMPSCE)
                                    printf ("-------	-------	-------	-------	-------	-------	");
                                else if (refcat == TMXSC)
                                    printf ("-------	-------	-------	------	");
                                else if (nmagr > 0) {
                                    for (imag = 0; imag < nmagr; imag++)
                                        printf ("-----	");
                                }
                                if (typecol == 1)
                                    printf ("---	");
                                if ((starcat[icat]!=NULL && starcat[icat]->entrv>0) &&
                                        mprop == 2)
                                    printf ("---------	");
                                if (mprop == 1)
                                    printf ("------	------	");
                                if (refcat == UCAC2 || refcat == UCAC3)
                                    printf ("-----	-----	--	--	");
                                if (refcat == UB1)
                                    printf ("--	--	--	");
                                printf ("--	");
                                for (i = 0; i < LenNum(das,2); i++)
                                    printf ("-");
                                printf ("	");
                                for (i = 0; i < LenNum(dds,2); i++)
                                    printf ("-");
                                printf ("	");
                                for (i = 0; i < LenNum(drs,2); i++)
                                    printf ("-");
                                printf ("\n");
                                nohead = 0;
                            }
                            if (srch != NULL) {
                                if (srchcat->keyid[0] > 0 && strlen (srch->objname))
                                    strcpy (numstr, srch->objname);
                                else if (srchcat->stnum <= 0 &&
                                                 strlen (srch->objname) > 0)
                                    strcpy (numstr, srch->objname);
                                else
                                    DS_CatNum (TXTCAT,-srchcat->nnfld,srchcat->nndec,
                                                    srch->num,numstr);
                                if (tabout)
                                    printf ("%s	", numstr);
                                else
                                    printf ("%s ", numstr);
                            }
                            if (degout) {
                                num2str (rastr, crao, 12, nddeg);
                                num2str (decstr, cdeco, 12, nddeg);
                            }
                            else {
                                ra2str (rastr, 32, crao, ndra);
                                dec2str (decstr, 32, cdeco, nddec);
                            }
                            if (tabout)
                                printf ("%s	%s", rastr, decstr);
                            else
                                printf ("%s %s", rastr, decstr);
                            if (srchcat != NULL && srch != NULL) {
                                if (smag > 0) {
                                    for (imag = 0; imag < smag; imag++) {
                                        if (tabout) {
                                            if (srch->xmag[imag] > 100.0)
                                                printf ("	%5.2fL", srch->xmag[imag]-100.0);
                                            else
                                                printf ("	%5.2f", srch->xmag[imag]);
                                        }
                                        else {
                                            if (srch->xmag[imag] > 100.0)
                                                printf (" %5.2fL", srch->xmag[imag]-100.0);
                                            else
                                                printf (" %5.2f", srch->xmag[imag]);
                                        }
                                    }
                                }
                                if (srchcat->nepoch) {
                                    ep2dt (srch->epoch, &date, &time);
                                    if (tabout)
                                        printf ("	%9.4f", date);
                                    else
                                        printf (" %9.4f", date);
                                }
                                if (srchcat->sptype != 0) {
                                    if (padspt) {
                                        for (i = 0; i < 7; i++) {
                                            if (srch->isp[i] == (char) 0)
                                                srch->isp[i] = ' ';
                                        }
                                        srch->isp[7] = (char) 0;
                                        if (tabout)
                                            printf ("	%7s", srch->isp);
                                        else
                                            printf ("  %7s ", srch->isp);
                                    }
                                    else if (tabout)
                                        printf ("	%s", srch->isp);
                                    else
                                        printf ("  %s ", srch->isp);
                                }
                                if (srchcat->entrv > 0) {
                                    if (tabout)
                                        printf ("	%9.2f", srch->radvel);
                                    else
                                        printf (" %9.2f", srch->radvel);
                                }
                                if (srchcat->mprop == 1) {
                                    pra = srch->rapm * 3600000.0 * cosdeg (srch->dec);
                                    pdec = srch->decpm * 3600000.0;
                                    if (tabout)
                                        printf ("	%6.1f	%6.1f", pra, pdec);
                                    else
                                        printf (" %6.1f %6.1f", pra, pdec);
                                }
                            }
                            if (mprop == 1) {
                                pra = gpra[0] * 3600000.0 * cosdeg (gdec[0]);
                                pdec = gpdec[0] * 3600000.0;
                            }

                            /* Set up object name or number to print */
                            if (refcat == SDSS || refcat == GSC2 || refcat == SKYBOT)
                                strcpy (numstr, gobj[0]);
                            else if (starcat[icat] != NULL) {
                                if (starcat[icat]->stnum < 0 && gobj1 != NULL) {
                                    strncpy (numstr, gobj1[0], 79);
                                    if (lofld > 0) {
                                        for (j = 0; j < lofld; j++) {
                                            if (!numstr[j])
                                                numstr[j] = ' ';
                                        }
                                    }
                                }
                                else
                                    DS_CatNum (refcat,-nnfld,starcat[icat]->nndec,gnum[0],numstr);
                            }
                            else
                                DS_CatNum (refcat, -nnfld, nndec, gnum[0], numstr);
                            /* if (gobj1 != NULL) {
                                 if (strlen (gobj1[0]) > 0)
                                 strcpy (numstr, gobj1[0]);
                                 }
                                 if (starcat[icat] != NULL)
                                 DS_CatNum (refcat,-nnfld,starcat[icat]->nndec,gnum[0],numstr);
                                 else
                                 DS_CatNum (refcat,-nnfld,nndec,gnum[0],numstr);  */
                            if (degout) {
                                num2str (rastr, gra[0], 12, nddeg);
                                num2str (decstr, gdec[0], 12, nddeg);
                            }
                            else {
                                ra2str (rastr, 32, gra[0], ndra);
                                dec2str (decstr, 32, gdec[0], nddec);
                            }
                            if (tabout)
                                printf ("	%s	%s	%s",
                                                numstr, rastr, decstr);
                            else
                                printf (" %s %s %s",
                                                numstr, rastr, decstr);
                            if (refcat == UCAC2 || refcat == UCAC3) {
                                era = gm[nmagr][0] * cosdeg (gdec[i]) * 3600.0;
                                edec = gm[nmagr+1][0] * 3600.0;
                                printf ("	%5.3f	%5.3f", era, edec);
                            }
                            if (refcat == GSC2) {
                                if (tabout)
                                    printf ("	%5.2f	%5.2f	%5.2f	%5.2f %5.2f",
                                                    gm[0][0],gm[1][0],gm[2][0],gm[4][0],gm[5][0]);
                                else
                                    printf (" %5.2f %5.2f %5.2f %5.2f %5.2f",
                                                    gm[0][0],gm[1][0],gm[2][0],gm[4][0],gm[5][0]);
                            }
                            else if (refcat == HIP) {
                                if (tabout)
                                    printf ("	%5.2f	%5.2f	%5.2f	%5.2f",
                                                    gm[0][0], gm[1][0], gm[2][0], gm[3][0]);
                                else
                                    printf (" %5.2f %5.2f %5.2f %5.2f",
                                                    gm[0][0], gm[1][0], gm[2][0], gm[3][0]);
                            }
                            else if (refcat == IRAS) {
                                for (imag = 0; imag < 4; imag++) {
                                    if (gm[imag][0] > 100.0) {
                                        flux = 1000.0 * pow (10.0,-(gm[imag][0]-100.0)/2.5);
                                        if (tabout)
                                            printf ("	%.2fL", flux);
                                        else
                                            printf (" %5.2fL", flux);
                                    }
                                    else {
                                        flux = 1000.0 * pow (10.0,-gm[imag][0]/2.5);
                                        if (tabout)
                                            printf ("	%.2f ", flux);
                                        else
                                            printf (" %5.2f ", flux);
                                    }
                                }
                            }
                            else if (refcat == TMPSC || refcat == TMIDR2 ||
                                             refcat == TMPSCE || refcat == TMXSC) {
                                for (imag = 0; imag < 3; imag++) {
                                    if (gm[imag][0] > 100.0) {
                                        if (tabout)
                                            printf ("	%6.3fL", gm[imag][0]-100.0);
                                        else
                                            printf (" %6.3fL", gm[imag][0]-100.0);
                                    }
                                    else {
                                        if (tabout)
                                            printf ("	%6.3f ", gm[imag][0]);
                                        else
                                            printf (" %6.3f ", gm[imag][0]);
                                    }
                                }
                                if (refcat == TMPSCE) {
                                    for (imag = 3; imag < 6; imag++) {
                                        if (tabout)
                                            printf ("	%6.3f ", gm[imag][0]);
                                        else
                                            printf (" %6.3f ", gm[imag][0]);
                                    }
                                }
                                if (refcat == TMXSC) {
                                    if (tabout)
                                        printf("	%6.1f",((double)gc[0])* 0.1);
                                    else
                                        printf (" %6.1f", ((double)gc[0])*0.1);
                                }
                            }
                            else if (nmagr > 0) {
                                for (imag = 0; imag < nmagr; imag++) {
                                    if (tabout)
                                        printf ("	%5.2f", gm[imag][0]);
                                    else
                                        printf (" %5.2f", gm[imag][0]);
                                }
                            }
                            if (typecol == 1) {
                                isp[0] = gc[0] / 1000;
                                isp[1] = gc[0] % 1000;
                                if (isp[0] == ' ' && isp[1] == ' ') {
                                    isp[0] = '_';
                                    isp[1] = '_';
                                }
                                if (tabout)
                                    printf ("	%2s ", isp);
                                else
                                    printf (" %2s ", isp);
                            }
                            if (starcat[icat] != NULL && starcat[icat]->entrv>0) {
                                if (tabout)
                                    printf ("	%9.2f", gm[nmagr][0]);
                                else
                                    printf (" %9.2f", gm[nmagr][0]);
                            }
                            if (starcat[icat] != NULL && starcat[icat]->nepoch>0) {
                                if (starcat[icat]->entrv>0)
                                    dtemp = DS_DateString (gm[nmagr+1][0], tabout);
                                else
                                    dtemp = DS_DateString (gm[nmagr][0], tabout);
                                printf ("%s", dtemp);
                                free (dtemp);
                            }
                            if (mprop == 1) {
                                if (tabout)
                                    printf ("	%6.1f	%6.1f", pra, pdec);
                                else
                                    printf (" %6.1f %6.1f", pra, pdec);
                            }
                            if (refcat == UCAC2 || refcat == UCAC3) {
                                epmr = gm[nmagr+2][0] * cosdeg (gdec[i]) * 3600000.0;
                                epmd = gm[nmagr+3][0] * 3600000.0;
                                nim = gc[0] / 1000;
                                nct = gc[0] % 1000;
                                if (tabout)
                                    printf ("	%5.1f	%5.1f	%2d	%2d",
                                                    epmr, epmd, nim, nct);
                                else
                                    printf (" %5.1f %5.1f %3d %3d",
                                                    epmr, epmd, nim, nct);
                            }
                            if (refcat == UB1) {
                                if (tabout)
                                    printf ("	%2d	%2d	%2d",
                                                    gc[0]%10000/100, gc[0]%100, gc[0]/10000);
                                else
                                    printf (" %2d %2d %2d",
                                                    gc[0]%10000/100, gc[0]%100, gc[0]/10000);
                            }
                            if (refcat == GSC2) {
                                if (tabout)
                                    printf ("	%d", gc[0]);
                                else
                                    printf ("  %2d  ", gc[0]);
                            }

                            /* Number of stars in search radius */
                            if (tabout)
                                printf ("	%d", ng);
                            else
                                printf (" %d", ng);
                            dec = (gdec[0] + cdeco) * 0.5;
                            if (degout) {
                                if ((gra[0] - crao) > 180.0)
                                    da = gra[0] - crao - 360.0;
                                else if ((gra[0] - crao) < -180.0)
                                    da = gra[0] - crao + 360.0;
                                else
                                    da = gra[0] - crao;
                                dd = gdec[0] - cdeco;
                                gdist = sqrt (da*da + dd*dd);
                                ndist = 5;
                            }
                            else {
                                if ((gra[0] - crao) > 180.0)
                                    da = 3600.0*(gra[0]-crao-360.0)*cos(degrad(dec));
                                else if ((gra[0] - crao) < -180.0)
                                    da = 3600.0*(gra[0]+360.0-crao)*cos(degrad(dec));
                                else
                                    da = 3600.0 * (gra[0] - crao) * cos (degrad (dec));
                                dd = 3600.0 * (gdec[0] - cdeco);
                                gdist = 3600.0 * gx[0];
                                ndist = 2;
                            }
                            if (tabout)
                                printf ("	");
                            else
                                printf (" ");
                            PrintNum (das, da, ndist);
                            if (tabout)
                                printf ("	");
                            else
                                printf (" ");
                            PrintNum (dds, dd, ndist);
                            if (tabout)
                                printf ("	");
                            else
                                printf (" ");
                            PrintNum (drs, gdist, ndist);
                            printf ("\n");
                        }
                        notprinted = 0;
                        continue;
                    }

                    /* List the brightest or closest MAXSTARS reference stars */
                    if (sortmag > 0 && sortmag <= nmag)
                        magsort = sortmag - 1;
                    else
                        magsort = 0;
                    DS_CatMagName (sortmag, refcat, magname);
                    if (ng > ngmax) {
                        if ((verbose || printhead) && !closest) {
                            if (distsort) {
                                if (ng > 1)
                                    printf ("Closest %d / %d %s (closer than %.2f arcsec)",
                                                    ns, ng, title, 3600.0*gx[ns-1]);
                                else
                                    printf ("Closest of %d %s",ng, title);
                            }
                            else if (maglim1 > 0.0) {
                                double magmin, magmax;
                                magmin = gm[magsort][0];
                                magmax = gm[magsort][0];
                                for (is = 0; is < ns; is++) {
                                    if (gm[magsort][is] > magmax)
                                        magmax = gm[magsort][is];
                                    if (gm[magsort][is] < magmin)
                                        magmax = gm[magsort][is];
                                }
                                printf ("%d / %d %s (%s between %.2f and %.2f)",
                                                ns, ng, title, magname, magmin, magmax);
                            }
                            else {
                                double magmax;
                                magmax = gm[magsort][0];
                                for (is = 0; is < ns; is++) {
                                    if (gm[magsort][is] > magmax)
                                        magmax = gm[magsort][is];
                                }
                                printf ("%d / %d %s (%s brighter than %.2f)",
                                                ns, ng, title, magname, magmax);
                            }
                            printf ("\n");
                        }
                    }
                    else {
                        if (verbose || printhead) {
                            if (maglim1 > 0.0)
                                printf ("%d %s, %s between %.2f and %.2f\n",
                                                ng, title, magname, maglim1, maglim2);
                            else if (maglim2 > 0.0)
                                printf ("%d %s, %s brighter than %.2f\n",
                                                ng, title, magname, maglim2);
                            else if (verbose)
                                printf ("%d %s\n", ng, title);
                        }
                    }

                    /* Open result catalog file */
                    if (wfile && icat == 0) {
                        if (objname)
                            strcpy (filename,objname);
                        else
                            strcpy (filename,"search");
                        for (i = 0; i < ncat; i++) {
                            strcat (filename,".");
                            strcat (filename,refcatname[icat]);
                        }
                        if (printxy)
                            strcat (filename,".match");

                        if (afile)
                            fd = fopen (filename, "a");
                        else
                            fd = fopen (filename, "w");

                        /* Free result arrays and return if cannot write file */
                        if (fd == NULL) {
                            if (afile)
                                fprintf (stderr, "%s:  cannot append to file %s\n",
                                                 cpname, filename);
                            else
                                fprintf (stderr, "%s:  cannot write file %s\n",
                                                 cpname, filename);
                            return (0);
                        }
                    }
                }
      }
    if (webtable == 0) {

      /* Write heading */
      if (votab) {
                printf("ERROR: votab\n");
                exit(-1);
      }

      else if (ng == 0 && (!tabout || (tabout && !printabhead))) {
                if (verbose)
                    printf ("No %s Stars Found\n", title);
      }

      else if (tabout && nohead) {
                catalog = DS_CatName (refcat, refcatname[icat]);
                sprintf (headline, "catalog	%s", catalog);
                if (wfile)
                    fprintf (fd, "%s\n", headline);
                else
                    printf ("%s\n", headline);

                if (!ranges) {
                    if (sysout == WCS_GALACTIC) {
                        num2str (rastr, crao, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "glon	%s\n", rastr);
                        else
                            printf ("glon	%s\n", rastr);
                        num2str (decstr, cdeco, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "glat	%s\n", decstr);
                        else
                            printf ("glat	%s\n", decstr);
                    }
                    else if (sysout == WCS_ECLIPTIC) {
                        num2str (rastr, crao, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "elon	%s\n", rastr);
                        else
                            printf ("elon	%s\n", rastr);
                        num2str (decstr, cdeco, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "elat	%s\n", decstr);
                        else
                            printf ("elat	%s\n", decstr);
                    }
                    else if (degout) {
                        num2str (rastr, crao, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "ra	%s\n", rastr);
                        else
                            printf ("ra	%s\n", rastr);
                        num2str (decstr, cdeco, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "dec	%s\n", decstr);
                        else
                            printf ("dec	%s\n", decstr);
                    }
                    else {
                        ra2str (rastr, 32, crao, ndra);
                        if (wfile)
                            fprintf (fd, "ra	%s\n", rastr);
                        else
                            printf ("ra	%s\n", rastr);
                        dec2str (decstr, 32, cdeco, nddec);
                        if (wfile)
                            fprintf (fd, "dec	%s\n", decstr);
                        else
                            printf ("dec	%s\n", decstr);
                    }
                    if (syscoor == WCS_GALACTIC && syscoor != sysout) {
                        num2str (rastr, cra, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "glon	%s\n", rastr);
                        else
                            printf ("glon	%s\n", rastr);
                        num2str (decstr, cdec, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "glat	%s\n", decstr);
                        else
                            printf ("glat	%s\n", decstr);
                    }
                    if (syscoor == WCS_ECLIPTIC && syscoor != sysout) {
                        num2str (rastr, cra, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "elon	%s\n", rastr);
                        else
                            printf ("elon	%s\n", rastr);
                        num2str (decstr, cdec, 12, nddeg);
                        if (wfile)
                            fprintf (fd, "elat	%s\n", decstr);
                        else
                            printf ("elat	%s\n", decstr);
                    }
                }

                /* Minimum number of plate IDs for USNO-B1.0 catalog */
                if (minid != 0) {
                    sprintf (headline, "minid	%d", minid);
                    if (wfile)
                        fprintf (fd, "%s\n", headline);
                    if (tabout)
                        printf ("%s\n", headline);
                }

                /* Minimum proper motion quality for USNO-B1.0 catalog */
                if (minpmqual > 0) {
                    sprintf (headline, "minpmq	%d", minpmqual);
                    if (wfile)
                        fprintf (fd, "%s\n", headline);
                    if (tabout)
                        printf ("%s\n", headline);
                }

                if (wfile) {
                    if (sysout == WCS_GALACTIC)
                        fprintf (fd,"radecsys	galactic\n");
                    else if (sysout == WCS_ECLIPTIC)
                        fprintf (fd,"radecsys	ecliptic\n");
                    else if (sysout == WCS_B1950)
                        fprintf (fd,"radecsys	fk4\n");
                    else
                        fprintf (fd,"radecsys	fk5\n");
                    fprintf (fd, "equinox	%.4f\n", eqout);
                    if (!printepoch)
                        fprintf (fd, "epoch	%.4f\n", epout);
                }
                else {
                    if (sysout == WCS_GALACTIC)
                        printf ("radecsys	galactic\n");
                    else if (sysout == WCS_ECLIPTIC)
                        printf ("radecsys	ecliptic\n");
                    else if (sysout == WCS_B1950)
                        printf ("radecsys	fk4\n");
                    else
                        printf ("radecsys	fk5\n");
                    printf ("equinox	%.4f\n", eqout);
                    if (!printepoch)
                        printf ("epoch	%.4f\n", epout);
                }

                if (mprop == 1) {
                    if (wfile) {
                        fprintf (fd, "rpmunit	mas/year\n");
                        fprintf (fd, "dpmunit	mas/year\n");
                    }
                    else {
                        printf ("rpmunit	mas/year\n");
                        printf ("dpmunit	mas/year\n");
                    }
                }

                das = dra0  / 3600.0;
                dds = ddec0 / 3600.0;
                drs = rad0;
                if (drs < 0.0)
                    drs = -drs;
                else if (drs == 0.0)
                    drs = sqrt (das*das + dds*dds);
                if (das <= 0.0) {
                    das = drs;
                    dds = drs;
                }
                if (dra0 > 0.0 || rad0 < 0.0) {
                    if (syscoor == WCS_GALACTIC || syscoor == WCS_ECLIPTIC) {
                        if (wfile) {
                            fprintf (fd, "dlonsec	%.2f\n", das);
                            fprintf (fd, "dlatsec	%.2f\n", dds);
                        }
                        else {
                            printf ("dlonsec	%.2f\n", das);
                            printf ("dlatsec	%.2f\n", dds);
                        }
                    }
                    else {
                        if (wfile) {
                            fprintf (fd, "drasec	%.2f\n", das);
                            fprintf (fd, "ddecsec	%.2f\n", dds);
                        }
                        else {
                            printf ("drasec	%.2f\n", das);
                            printf ("ddecsec	%.2f\n", dds);
                        }
                    }
                }
                else if (rad0 > 0) {
                    if (wfile)
                        fprintf (fd, "radsec	%.2f\n", drs);
                    else
                        printf ("radsec	%.2f\n", drs);
                }

                if (catsort > 0) {
                    switch (catsort) {
                    case SORT_DEC:
                        if (wfile)
                            fprintf (fd, "catsort	dec\n");
                        else
                            printf ("catsort	dec\n");
                        break;
                    case SORT_DIST:
                        if (wfile)
                            fprintf (fd, "catsort	dist\n");
                        else
                            printf ("catsort	dist\n");
                        break;
                    case SORT_MAG:
                        if (wfile)
                            fprintf (fd, "catsort	date\n");
                        else
                            printf ("catsort	date\n");
                        break;
                    case SORT_RA:
                        if (wfile)
                            fprintf (fd, "catsort	ra\n");
                        else
                            printf ("catsort	ra\n");
                        break;
                    case SORT_ID:
                        if (wfile)
                            fprintf (fd, "catsort	id\n");
                        else
                            printf ("catsort	id\n");
                        break;
                    default:
                        break;
                    }
                }

                if (wfile)
                    fprintf (fd, "program	%s %s\n", progname, RevMsg);
                else
                    printf ("program	%s %s\n", progname, RevMsg);

                /* Print column headings */
                if (refcat == TABCAT && strlen(starcat[icat]->keyid) > 0)
                    sprintf (headline,"%s          ", starcat[icat]->keyid);
                else if (refcat == SKYBOT) {
                    strcpy (headline, "object");
                    for (i = 6; i < nnfld; i++)
                        strcat (headline, " ");
                }
                else
                    DS_CatID (headline, refcat);
                headline[nnfld] = (char) 0;

                if (sysout == WCS_GALACTIC)
                    strcat (headline,"	long_gal   	lat_gal  ");
                else if (sysout == WCS_ECLIPTIC)
                    strcat (headline,"	long_ecl   	lat_ecl  ");
                else if (sysout == WCS_B1950)
                    strcat (headline,"	ra1950      	dec1950  ");
                else
                    strcat (headline,"	ra      	dec      ");
                if (refcat == UCAC2 || refcat == UCAC3)
                    strcat (headline,"	raerr	decerr");
                if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
                        refcat == UAC  || refcat == UA1  || refcat == UA2)
                    strcat (headline,"	magb	magr	plate");
                else if (refcat == TMPSC || refcat == TMIDR2)
                    strcat (headline,"	magj  	magh  	magk  ");
                else if (refcat == TMPSCE)
                    strcat (headline,"	magj  	magh  	magk 	magje 	maghe 	magke ");
                else if (refcat == TMXSC)
                    strcat (headline,"	magj  	magh  	magk  	size  ");
                else if (refcat == UB1)
                    strcat (headline, "	magb1	magr1	magb2	magr2	magn 	pm	ni	sg");
                else if (refcat == YB6)
                    strcat (headline, "	magb 	magr 	magj 	magh 	magk ");
                else if (refcat == SDSS)
                    strcat (headline, "	magu 	magg 	magr 	magi 	magz ");
                else if (refcat == GSC2)
                    strcat (headline,"	magf	magj	magn  	magv ");
                else if (refcat == UCAC2)
                    strcat (headline,"	magj	magh	magk  	magc ");
                else if (refcat == UCAC3)
                    strcat (headline,"	magb 	magr 	magi 	magj 	magh 	magk 	magm 	maga ");
                else if (refcat == SKYBOT)
                    strcat (headline,"	magv	gdist	hdist ");
                else if (refcat == IRAS)
                    strcat (headline,"	f10m  	f25m  	f60m   	f100m ");
                else if (refcat == HIP)
                    strcat (headline,"	magb	magv	parlx 	parer");
                else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
                    strcat (headline,"	magb	magv");
                else if (refcat==TYCHO2E)
                    strcat (headline,"	magb 	magv 	magbe	magve");
                else if (refcat==SKY2K)
                    strcat (headline,"	magb 	magv 	magph	magpv");
                else if (refcat==HIP)
                    strcat (headline,"	magb 	magv 	prllx	parer");
                else if (refcat==SKYBOT)
                    strcat (headline,"	magv 	gdist  	hdist");
                else if (refcat == GSC || refcat == GSCACT)
                    strcat (headline,"	mag	class	band	N");
                else if (refcat == UJC)
                    strcat (headline,"	mag	plate");
                else if (nmagr > 0) {
                    for (imag = 0; imag < nmagr; imag++) {
                        if (starcat[icat] != NULL &&
                                strlen (starcat[icat]->keymag[imag]) > 0)
                            sprintf (temp, "	%s ", starcat[icat]->keymag[imag]);
                        else if (nmagr > 1)
                            sprintf (temp, "	mag%d ", imag);
                        else
                            sprintf (temp, "	mag  ");
                        strcat (headline, temp);
                    }
                }
                if (typecol == 1)
                    strcat (headline,"	type");
                if (printepoch)
                    strcat (headline, "	epoch     ");
                if (mprop == 2)
                    strcat (headline," 	velocity");
                if (mprop == 1)
                    strcat (headline,"	pmra  	pmdec ");
                if (refcat == UCAC2 || refcat == UCAC3)
                    strcat (headline,"	epmra	epmdec	ni	nc");
                if (refcat == GSC2)
                    strcat (headline,"	class");
                if (ranges == NULL)
                    strcat (headline,"	arcsec");
                if (refcat == TABCAT && keyword != NULL) {
                    strcat (headline,"	");
                    strcat (headline, keyword);
                }
                if (catsort == SORT_MERGE)
                    strcat (headline,"	nmatch");
                if (gobj1 != NULL && refcat != GSC2 && refcat != SDSS && refcat != SKYBOT) {
                    if (starcat[icat] == NULL ||
                            (starcat[icat] != NULL && starcat[icat]->stnum > 0))
                        strcat (headline,"	object");
                }
                if (printxy)
                    strcat (headline, "	X      	Y      ");
                if (wfile)
                    fprintf (fd, "%s\n", headline);
                else
                    printf ("%s\n", headline);

                strcpy (headline, "-------------------------");
                headline[nnfld] = (char) 0;
                strcat (headline,"	------------	------------");
                if (refcat == UCAC2 || refcat == UCAC3)
                    strcat (headline,"	-----	-----");
                if (refcat == TMPSC || refcat == TMIDR2)
                    strcat (headline,"	-------	-------	-------");
                else if (refcat == TMPSCE)
                    strcat (headline,"	-------	-------	-------	------	------	------");
                else if (refcat == TMXSC)
                    strcat (headline,"	-------	-------	-------	------");
                else if (refcat == IRAS)
                    strcat (headline,"	-----	-----	-----	-----");
                else if (refcat == HIP)
                    strcat (headline,"	-----	-----	-----	-----");
                else if (refcat == GSC2)
                    strcat (headline,"	-----	-----	-----	-----");
                else if (nmagr > 0) {
                    for (imag = 0; imag < nmagr; imag++)
                        strcat (headline,"	-----");
                }
                if (refcat == GSC || refcat == GSCACT)
                    strcat (headline,"	-----	----	-");
                else if (refcat == UB1)
                    strcat (headline,"	--	--	--");
                else if (typecol == 1)
                    strcat (headline,"	----");
                else if (typecol == 2)
                    strcat (headline,"	-----");
                if (printepoch)
                    strcat (headline, "	----------");
                if (mprop == 2)
                    strcat (headline,"	--------");
                if (mprop == 1)
                    strcat (headline,"	------	------");
                if (refcat == UCAC2 || refcat == UCAC3)
                    strcat (headline,"	-----	-----	--	--");
                if (refcat == GSC2)
                    strcat (headline,"	-----");
                if (ranges == NULL)
                    strcat (headline, "	------");
                if (refcat == TABCAT && keyword != NULL)
                    strcat (headline,"	------");
                if (catsort == SORT_MERGE)
                    strcat (headline,"	------");
                if (gobj1 != NULL && refcat != GSC2 && refcat != SDSS && refcat != SKYBOT) {
                    if (starcat[icat] == NULL ||
                            starcat[icat]->stnum > 0)
                        strcat (headline,"	------");
                }
                if (printxy)
                    strcat (headline, "	-------	-------");
                if (wfile)
                    fprintf (fd, "%s\n", headline);
                else
                    printf ("%s\n", headline);
                nohead = 0;
      }

      else if (printhead && nohead) {
                if (printxy)
                    strcpy (headline, "  X     Y   ");
                else {
                    if (refcat == GSC || refcat == GSCACT)
                        strcpy (headline, "GSC_number ");
                    else if (refcat == GSC2)
                        strcpy (headline, "GSC2_id     ");
                    else if (refcat == USAC)
                        strcpy (headline, "USNO_SA_number ");
                    else if (refcat == USA1)
                        strcpy (headline, "USNO_SA1_number");
                    else if (refcat == USA2)
                        strcpy (headline, "USNO_SA2_number");
                    else if (refcat == UAC)
                        strcpy (headline, "USNO_A_number  ");
                    else if (refcat == UA1)
                        strcpy (headline, "USNO_A1_number ");
                    else if (refcat == UA2)
                        strcpy (headline, "USNO_A2_number ");
                    else if (refcat == UB1)
                        strcpy (headline, "USNO_B1_number ");
                    else if (refcat == YB6)
                        strcpy (headline, "USNO_YB6_number ");
                    else if (refcat == SDSS)
                        strcpy (headline, "SDSS_number          ");
                    else if (refcat == SKYBOT)
                        strcpy (headline, "Object          ");
                    else if (refcat == UCAC1)
                        strcpy (headline, "UCAC1_num    ");
                    else if (refcat == UCAC2)
                        strcpy (headline, "UCAC2_num    ");
                    else if (refcat == TMPSC || refcat == TMPSCE)
                        strcpy (headline, "2MASS_num.  ");
                    else if (refcat == TMXSC)
                        strcpy (headline, "2MASS_XSC   ");
                    else if (refcat == TMIDR2)
                        strcpy (headline, "2MIDR2_num.");
                    else if (refcat == UJC)
                        strcpy (headline, " UJ_number    ");
                    else if (refcat == SAO)
                        strcpy (headline, "SAO_number ");
                    else if (refcat == PPM)
                        strcpy (headline, "PPM_number ");
                    else if (refcat == SKY2K)
                        strcpy (headline, "SKY_id    ");
                    else if (refcat == BSC)
                        strcpy (headline, "BSC_number ");
                    else if (refcat == IRAS)
                        strcpy (headline, "IRASnum  ");
                    else if (refcat == TYCHO)
                        strcpy (headline, "Tycho_number ");
                    else if (refcat == TYCHO2 || refcat == TYCHO2E)
                        strcpy (headline, "Tycho2_num  ");
                    else if (refcat == HIP)
                        strcpy (headline, "Hip_num ");
                    else if (refcat == ACT)
                        strcpy (headline, "ACT_number  ");
                    else if (nnfld > 5) {
                        strcpy (headline, "Number      ");
                        headline[nnfld+3] = (char) 0;
                    }
                    else {
                        strcpy (headline, "ID       ");
                        headline[nnfld+3] = (char) 0;
                    }
                }
                if (sysout == WCS_B1950) {
                    if (degout) {
                        if (eqout == 1950.0)
                            strcat (headline, "   RA1950     Dec1950   ");
                        else {
                            sprintf (temp, " RAB%7.2f   DecB%7.2f  ", eqout, eqout);
                            strcat (headline, temp);
                        }
                    }
                    else {
                        if (eqout == 1950.0)
                            strcat (headline, "RA1950      Dec1950    ");
                        else {
                            sprintf (temp, "RAB%7.2f  DecB%7.2f ", eqout, eqout);
                            strcat (headline, temp);
                        }
                    }
                }
                else if (sysout == WCS_ECLIPTIC)
                    strcat (headline, "Ecl Lon    Ecl Lat  ");
                else if (sysout == WCS_GALACTIC)
                    strcat (headline, "Gal Lon    Gal Lat  ");
                else {
                    if (degout) {
                        if (eqout == 2000.0)
                            strcat (headline, "   RA2000     Dec2000 ");
                        else {
                            sprintf (temp," RAJ%7.2f    DecJ%7.2f ", eqout, eqout);
                            strcat (headline, temp);
                        }
                    }
                    else {
                        if (eqout == 2000.0)
                            strcat (headline, " RA2000       Dec2000   ");
                        else {
                            sprintf (temp,"RAJ%7.2f   DecJ%7.2f  ", eqout, eqout);
                            strcat (headline, temp);
                        }
                    }
                }
                if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
                        refcat == UAC  || refcat == UA1  || refcat == UA2)
                    strcat (headline, "  MagB  MagR");
                else if (refcat == UB1)
                    strcat (headline, "  MagB1  MagR1  MagB2  MagR2  MagN ");
                else if (refcat == YB6)
                    strcat (headline, "  MagB   MagR   MagJ   MagH   MagK ");
                else if (refcat == SDSS)
                    strcat (headline, "  Magu   Magg   Magr   Magi   Magz ");
                else if (refcat == UJC)
                    strcat (headline, "  Mag ");
                else if (refcat == GSC || refcat == GSCACT)
                    strcat (headline, "   Mag ");
                else if (refcat==SAO || refcat==PPM || refcat==BSC || refcat==UCAC1)
                    strcat (headline, "    Mag");
                else if (refcat==SKY2K)
                    strcat (headline, "   MagB   MagV  MagPh  MagPv");
                else if (refcat==UCAC2)
                    strcat (headline, "   MagJ   MagH   MagK   MagC");
                else if (refcat == UCAC3)
                    strcat (headline,"    MagB   MagR   MagI   MagJ   MagH   MagK   MagM   MagA");
                else if (refcat==SKYBOT)
                    strcat (headline, "     MagV   GDist  HDist");
                else if (refcat==TMPSC || refcat == TMIDR2)
                    strcat (headline, "   MagJ    MagH    MagK  ");
                else if (refcat==TMPSCE)
                    strcat (headline, "   MagJ    MagH    MagK   MagJe  MagHe  MagKe");
                else if (refcat==TMXSC)
                    strcat (headline, "   MagJ    MagH    MagK  ");
                else if (refcat==IRAS)
                    strcat (headline, "  f10m   f25m   f60m   f100m");
                else if (refcat==GSC2)
                    strcat (headline, "  MagF  MagJ  MagN  MagV");
                else if (refcat==HIP)
                    strcat (headline, "  MagB  MagV  Parlx Parer");
                else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
                    strcat (headline, "   MagB   MagV  ");
                else if (refcat==TYCHO2E)
                    strcat (headline, "   MagB   MagV MagBe MagVe");
                else if (nmagr > 0) {
                    for (imag = 0; imag < nmagr; imag++) {
                        if (strlen(starcat[icat]->keymag[imag]) > 0) {
                            strcat (headline," ");
                            sprintf (temp, " %s ", starcat[icat]->keymag[imag]);
                            strcat (headline, temp);
                        }
                        else if (nmagr > 1) {
                            sprintf (temp, "   Mag%d", imag+1);
                            strcat (headline, temp);
                        }
                        else
                            strcat (headline, "  Mag");
                    }
                }
                if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
                        refcat == UAC  || refcat == UA1  || refcat == UA2 ||
                        refcat == UJC)
                    strcat (headline, "  Plate");
                else if (refcat == UB1)
                    strcat (headline, " PM NI SG");
                else if (refcat==GSC2)
                    strcat (headline, " Class");
                else if (refcat == GSC || refcat == GSCACT)
                    strcat (headline, " Class Band N");
                else if (refcat==TMXSC)
                    strcat (headline, "   Size");
                else if (typecol == 1)
                    strcat (headline, " Type");
                else if (gcset && refcat != UCAC2 && refcat != UCAC3)
                    strcat (headline, "     Peak");
                if (printepoch)
                    strcat (headline, "   Epoch   ");
                /* if (mprop == 1)
                     strcat (headline, "   pmRA  pmDec"); */
                if (mprop == 2)
                    strcat (headline, " Velocity");
                if (refcat == UCAC2 || refcat == UCAC3)
                    strcat (headline, " nim ncat");
                if ((ranges == NULL) && (printdist == 1))
                    strcat (headline, "  Arcsec");
                if (gobj1 != NULL && refcat != GSC2 && refcat != SDSS && refcat != SKYBOT) {
                    if (starcat[icat] == NULL || starcat[icat]->stnum > 0) {
                        if (starcat[icat]->refcat == SCANCAT)
                            strcat (headline,"  Class");
                        else
                            strcat (headline,"  Object");
                    }
                }
                if (catsort == SORT_MERGE)
                    strcat (headline, "  Nmatch");
                if (wfile)
                    fprintf (fd, "%s\n", headline);
                else
                    printf ("%s\n", headline);
      }
      nohead = 0;

      /* Find maximum separation for formatting */
      gdmax = 0.0;
      for (i = 0; i < ns; i++) {
                if (gx[i] > 0.0) {
                    gdist = 3600.0 * gx[i];
                    if (gdist > gdmax)
                        gdmax = gdist;
                }
      }

      if (closest) ns = 1;
      for (i = 0; i < ns; i++) {

                /* Set source position epoch passed as a magnitude */
                if (printepoch) {
                    if (mprop == 2)
                        epoch = gm[nmagr+1][i];
                    else
                        epoch = gm[nmagr][i];
                    if (epoch == 99.0)
                        epoch = 0.0;
                }

                /* Check entry epoch to see if it is within specified range */
                if ((epoch1 > 0.0 && epoch2 > 0.0) &&
                        (epoch < epoch1 || epoch > epoch2))
                    continue;
                else if (epoch1 > 0.0 && epoch < epoch1)
                    continue;
                else if (epoch2 > 0.0 && epoch > epoch2)
                    continue;

                /* Set spectra type (or other 2-character code) */
                if (typecol == 1) {
                    isp[0] = gc[i] / 1000;
                    isp[1] = gc[i] % 1000;
                    if (isp[0] == ' ' && isp[1] == ' ') {
                        isp[0] = '_';
                        isp[1] = '_';
                    }
                }

                /* For HST Guide Star Catalog, set number of entries and band */
                if (refcat == GSC || refcat == GSCACT) {
                    ngsc = gc[i] / 10000;
                    gc[i] = gc[i] - (ngsc * 10000);
                    band = gc[i] / 100;
                    gc[i] = gc[i] - (band * 100);
                }

                /* For USNO-B1.0 catalog, drop sources on too few plates */
                if (refcat == UB1 && minid > 0) {
                    nid = gc[i]%100;
                    if (nid < minid)
                        continue;
                }

                /* For USNO UCAC2 and UCAC3 catalogs, set errors and number of epochs */
                if (refcat == UCAC2 || refcat == UCAC3) {
                    era = gm[nmagr][i] * cosdeg (gdec[i]) * 3600.0;
                    edec = gm[nmagr+1][i] * 3600.0;
                    epmr = gm[nmagr+2][i] * cosdeg (gdec[i]) * 3600000.0;
                    epmd = gm[nmagr+3][i] * 3600000.0;
                    nim = gc[i] / 1000;
                    nct = gc[i] % 1000;
                }

                if ((gy[i] > 0.0)  || (printdist == 0))

                    {
                        if (degout) {
                            deg2str (rastr, 32, gra[i], nddeg);
                            deg2str (decstr, 32, gdec[i], nddeg);
                        }
                        else {
#ifdef DEBUG_NOTAVAILABLE
                            ra2str (rastr, 32, gra[i], ndra);
                            dec2str (decstr, 32, gdec[i], nddec);
#else /* DEBUG_NOTAVAILABLE */
              if ((gra[i] < 999.0) && (gdec[i] < 99.0)) {
                ra2str (rastr, 32, gra[i], ndra);
                dec2str (decstr, 32, gdec[i], nddec);
              } else {
                strcpy(rastr,DASCH_NOTAVAILABLE);
                strcpy(decstr,DASCH_NOTAVAILABLE);
              }
#endif /* DEBUG_NOTAVAILABLE */

                        }
                        if (gx[i] > 0.0)
                            gdist = 3600.0 * gx[i];
                        else
                            gdist = 0.0;

                        /* Convert proper motion to milliarcsec/year from deg/year */
                        if (mprop == 1) {
                            pra = gpra[i] * 3600000.0 * cosdeg (gdec[i]);
                            pdec = gpdec[i] * 3600000.0;
                        }

                        /* Set up object name or number to print */
                        if (starcat[icat] != NULL) {
                            if (starcat[icat]->stnum < 0 && gobj1 != NULL) {
                                strncpy (numstr, gobj1[i], 32);
                                if (lofld > 0) {
                                    for (j = 0; j < lofld; j++) {
                                        if (!numstr[j])
                                            numstr[j] = ' ';
                                    }
                                }
                            }
                            else
                                if (starcat[icat]->refcat == SCANCAT) {
                                    ScanCatNum(starcat[icat],-nnfld,starcat[icat]->nndec,gnum[i],numstr);
                                } else {
                                    DS_CatNum (refcat,-nnfld,starcat[icat]->nndec,gnum[i],numstr);
                }
                        }

                        else if (refcat == SDSS || refcat == GSC2 || refcat == SKYBOT) {
                            strcpy (numstr, gobj[i]);
                            lobj = strlen (gobj[i]);
                            if (lobj < nnfld) {
                                for (j = lobj; j < nnfld; j++)
                                    strcat (numstr, " ");
                            }
                        }
                        else
                            DS_CatNum (refcat, -nnfld, nndec, gnum[i], numstr);

                        if (votab) {
                            sprintf (headline, "<tr>\n<td>%s</td><td>%s</td><td>%s</td>",
                                             numstr,rastr,decstr);
                            if (refcat == UCAC2 || refcat == UCAC3) {
                                sprintf (temp,"<td>%5.1f</td><td>%5.1f</td>", era, edec);
                                strcat (headline, temp);
                            }
                            for (imag = 0; imag < nmagr; imag++) {
                                sprintf (temp, "<td>%.2f</td>", gm[imag][i]);
                                strcat (headline, temp);
                            }
                            if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
                                    refcat == UAC  || refcat == UA1  || refcat == UA2) {
                                sprintf (temp, "<td>%d</td>", gc[i]);
                                strcat (headline, temp);
                            }
                            else if (refcat == GSC || refcat == GSCACT) {
                                sprintf (temp, "<td>%d</td><td>%d</td><td>%d</td>",
                                                 gc[i], band, ngsc);
                                strcat (headline, temp);
                            }
                            if (typecol == 1) {
                                sprintf (temp, "<td>%2s</td>", isp);
                                strcat (headline, temp);
                            }
                            if (mprop == 2) {
                                sprintf (temp, "<td>%8.2f</td>", gm[nmagr][i]);
                                strcat (headline, temp);
                            }
                            if (mprop == 1) {
                                sprintf (temp, "<td>%6.1f</td><td>%6.1f</td>", pra, pdec);
                                strcat (headline, temp);
                            }
                            if (refcat == UCAC2 || refcat == UCAC3) {
                                sprintf (temp,"<td>%5.1f</td><td>%5.1f</td><td>%2d</td><td>%2d</td>",
                                                 epmr, epmd, nim, nct);
                                strcat (headline, temp);
                            }
                            sprintf (temp, "<td>%.5f</td>\n</tr>", gdist / 3600.0);
                            strcat (headline, temp);

                            printf ("%s\n", headline);
                        }

                        /* Print or write tab-delimited output line for one star */
                        else if (tabout) {
                            if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
                                    refcat == UAC  || refcat == UA1  || refcat == UA2)
                                sprintf (headline, "%s	%s	%s	%.1f	%.1f	%d",
                                                 numstr, rastr, decstr, gm[0][i], gm[1][i], gc[i]);
                            else if (refcat == GSC || refcat == GSCACT)
                                sprintf (headline,
                                                 "%s	%s	%s	%.2f	%d	%d	%d",
                                                 numstr, rastr, decstr, gm[0][i], gc[i], band, ngsc);
                            else if (refcat==TMPSC || refcat==TMIDR2 ||
                                             refcat == TMPSCE || refcat==TMXSC) {
                                sprintf (headline, "%s	%s	%s", numstr, rastr, decstr);
                                for (imag = 0; imag < 3; imag++) {
                                    if (gm[imag][i] > 100.0)
                                        sprintf (temp, "	%6.3fL", gm[imag][i]-100.0);
                                    else
                                        sprintf (temp, "	%6.3f ", gm[imag][i]);
                                    strcat (headline, temp);
                                }
                                if (refcat == TMPSCE) {
                                    for (imag = 3; imag < 6; imag++) {
                                        sprintf (temp, "	%5.3f ", gm[imag][i]);
                                        strcat (headline, temp);
                                    }
                                }
                            }
                            else if (refcat == GSC2)
                                sprintf (headline, "%s	%s	%s	%.2f	%.2f	%.2f	%.2f",
                                                 numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i]);
                            else if (refcat == IRAS) {
                                sprintf (headline, "%s	%s	%s", numstr, rastr, decstr);
                                for (imag = 0; imag < 4; imag++) {
                                    if (gm[imag][i] > 100.0) {
                                        flux = 1000.0 * pow (10.0, -(gm[imag][i]-100.0) / 2.5);
                                        sprintf (temp, "	%.2fL", flux);
                                    }
                                    else {
                                        flux = 1000.0 * pow (10.0, -gm[imag][i] / 2.5);
                                        sprintf (temp, "	%.2f ", flux);
                                    }
                                    strcat (headline, temp);
                                }
                            }
                            else if (refcat == UB1)
                                sprintf (headline, "%s	%s	%s	%.2f	%.2f	%.2f	%.2f	%.2f	%2d	%2d	%2d",
                                                 numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i],gm[4][i],
                                                 gc[i]%10000/100, gc[i]%100, gc[i]/10000);
                            else if (refcat == YB6)
                                sprintf (headline, "%s	%s	%s	%.2f	%.2f	%.2f	%.2f	%.2f",
                                                 numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i],gm[4][i]);
                            else if (refcat == SDSS)
                                sprintf (headline, "%s	%s	%s	%.2f	%.2f	%.2f	%.2f	%.2f",
                                                 numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i],gm[4][i]);
                            else if (refcat == HIP)
                                sprintf (headline, "%s	%s	%s	%.2f	%.2f	%.2f	%.2f",
                                                 numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i]);
                            else if (refcat == UJC)
                                sprintf (headline, "%s	%s	%s	%.2f	%d",
                                                 numstr, rastr, decstr, gm[0][i], gc[i]);
                            else if (refcat==SAO || refcat==PPM || refcat == BSC)
                                sprintf (headline, "%s	%s	%s	%.2f",
                                                 numstr, rastr, decstr, gm[0][i]);
                            else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
                                sprintf (headline, "%s	%s	%s	%.2f	%.2f",
                                                 numstr, rastr, decstr, gm[0][i], gm[1][i]);
                            else if (refcat==TYCHO2E)
                                sprintf (headline, "%s	%s	%s	%.2f	%.2f	%.2f	%.2f",
                                                 numstr, rastr, decstr, gm[0][i], gm[1][i], gm[2][i], gm[3][i]);
                            else {
                                sprintf (headline, "%s	%s	%s", numstr, rastr, decstr);
                                if (refcat == UCAC2 || refcat == UCAC3) {
                                    sprintf (temp,"	%5.3f	%5.3f", era, edec);
                                    strcat (headline, temp);
                                }

                                for (imag = 0; imag < nmagr; imag++) {
                                    sprintf (temp, "	%.2f", gm[imag][i]);
                                    strcat (headline, temp);
                                }
                                if (typecol == 2) {
                                    sprintf (temp, "	%d", gc[i]);
                                    strcat (headline, temp);
                                }
                            }
                            if (typecol == 1) {
                                sprintf (temp, "	%2s", isp);
                                strcat (headline, temp);
                            }
                            if (mprop == 2) {
                                if (printepoch) {
                                    dtemp = DS_DateString (epoch, tabout);
                                    strcat (headline, dtemp);
                                    free (dtemp);
                                }
                                sprintf (temp, "	%8.2f", gm[nmagr][i]);
                                strcat (headline, temp);
                            }
                            else if (printepoch) {
                                dtemp = DS_DateString (epoch, tabout);
                                strcat (headline, dtemp);
                                free (dtemp);
                            }
                            if (mprop == 1) {
                                sprintf (temp, "	%6.1f	%6.1f", pra, pdec);
                                strcat (headline, temp);
                            }
                            if (refcat == TABCAT && gcset) {
                                sprintf (temp, "	%d", gc[i]);
                                strcat (headline, temp);
                            }
                            else if (refcat == TMXSC) {
                                sprintf (temp,"	%6.1f", ((double)gc[i])* 0.1);
                                strcat (headline, temp);
                            }
                            else if (refcat == GSC2) {
                                sprintf (temp, "	%d  ", gc[i]);
                                strcat (headline, temp);
                            }
                            else if (refcat == UCAC2 || refcat == UCAC3) {
                                sprintf (temp, "	%5.1f	%5.1f	%2d	%2d",
                                                 epmr, epmd, nim, nct);
                                strcat (headline, temp);
                            }

                            if (ranges == NULL) {
                                sprintf (temp, "	%.2f", gdist);
                                strcat (headline, temp);
                            }
                            if (refcat == TABCAT && keyword != NULL) {
                                strcat (headline, "	");
                                if (gobj[i] != NULL)
                                    strcat (headline, gobj[i]);
                                else
                                    strcat (headline, "___");
                            }
                            if (catsort == SORT_MERGE) {
                                sprintf (temp, "	%d", (int)gx[i]);
                                strcat (headline, temp);
                            }
                            if ((refcat == BINCAT || refcat == TXTCAT) &&
                                    gobj1 != NULL && gobj[i] != NULL) {
                                if (starcat[icat] == NULL || starcat[icat]->stnum > 0) {
                                    strcat (headline, "	");
                                    strcat (headline, gobj[i]);
                                }
                            }
                            if (printxy) {
                                strcat (headline, "	");
                                strcat (numstr, xstr);
                                strcat (headline, "	");
                                strcat (numstr, ystr);
                            }

                            if (wfile)
                                fprintf (fd, "%s\n", headline);
                            else
                                printf ("%s\n", headline);
                        }

                        /* Print or write space-delimited output line for one star */
                        else {
                            if (printxy) {
                                strcpy (numstr, xstr);
                                strcat (numstr, " ");
                                strcat (numstr, ystr);
                            }
                            if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
                                    refcat == UAC  || refcat == UA1  || refcat == UA2)
                                sprintf (headline,"%s %s %s %5.1f %5.1f ",
                                                 numstr,rastr,decstr,gm[0][i],gm[1][i]);
                            else if (refcat == TMPSC || refcat == TMIDR2 ||
                                             refcat == TMPSCE || refcat == TMXSC) {
                                sprintf (headline, "%s %s %s", numstr, rastr, decstr);
                                for (imag = 0; imag < 3; imag++) {
                                    if (gm[imag][i] > 100.0)
                                        sprintf (temp, " %6.3fL", gm[imag][i]-100.0);
                                    else
                                        sprintf (temp, " %6.3f ", gm[imag][i]);
                                    strcat (headline, temp);
                                }
                                if (refcat == TMPSCE) {
                                    for (imag = 3; imag < 6; imag++) {
                                        sprintf (temp, " %5.3f ", gm[imag][i]);
                                        strcat (headline, temp);
                                    }
                                }
                            }
                            else if (refcat == IRAS) {
                                sprintf (headline, "%s %s %s", numstr, rastr, decstr);
                                for (imag = 0; imag < 4; imag++) {
                                    if (gm[imag][i] > 100.0) {
                                        flux = 1000.0 * pow (10.0, -(gm[imag][i]-100.0) / 2.5);
                                        sprintf (temp, " %5.2fL", flux);
                                    }
                                    else {
                                        flux = 1000.0 * pow (10.0, -gm[imag][i] / 2.5);
                                        sprintf (temp, " %5.2f ", flux);
                                    }
                                    strcat (headline, temp);
                                }
                            }
                            else if (refcat == GSC2)
                                sprintf (headline, "%s %s %s %5.2f %5.2f %5.2f %5.2f",
                                                 numstr, rastr, decstr, gm[0][i], gm[1][i],
                                                 gm[2][i], gm[3][i]);
                            else if (refcat == HIP)
                                sprintf (headline, "%s %s %s %5.2f %5.2f %5.2f %5.2f",
                                                 numstr, rastr, decstr, gm[0][i], gm[1][i],
                                                 gm[2][i], gm[3][i]);
                            else if (refcat == GSC || refcat == GSCACT)
                                sprintf (headline, "%s %s %s %6.2f",
                                                 numstr, rastr, decstr, gm[0][i]);
                            else if (refcat == UJC)
                                sprintf (headline,"%s %s %s %6.2f",
                                                 numstr, rastr, decstr, gm[0][i]);
                            else if (refcat == UCAC1)
                                sprintf (headline,"%s  %s %s %6.2f",
                                                 numstr,rastr,decstr,gm[0][i]);
                            else if (refcat == UCAC2)
                                sprintf (headline,"%s  %s %s %6.2f %6.2f %6.2f %6.2f",
                                                 numstr,rastr,decstr,gm[0][i],
                                                 gm[1][i], gm[2][i], gm[3][i]);
                            else if (refcat==SAO || refcat==PPM || refcat == BSC)
                                sprintf (headline,"  %s  %s %s %6.2f",
                                                 numstr,rastr,decstr,gm[0][i]);
                            else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
                                sprintf (headline,"%s %s %s %6.2f %6.2f ",
                                                 numstr,rastr,decstr,gm[0][i],gm[1][i]);
                            else if (refcat==TYCHO2E)
                                sprintf (headline,"%s %s %s %6.2f %6.2f %5.2f %5.2f",
                                                 numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i]);
                            else {
                                sprintf (headline,"%s %s %s ", numstr, rastr, decstr);
                                for (imag = 0; imag < nmagr; imag++) {
                                    sprintf (temp, " %6.2f", gm[imag][i]);
                                    strcat (headline, temp);
                                }
                            }
                            if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
                                    refcat == UAC  || refcat == UA1  || refcat == UA2 ||
                                    refcat == UJC) {
                                sprintf (temp," %4d", gc[i]);
                                strcat (headline, temp);
                            }
                            else if (refcat == UB1) {
                                sprintf (temp," %2d %2d %2d",
                                                 gc[i]%10000/100, gc[i]%100, gc[i]/10000);
                                strcat (headline, temp);
                            }
                            else if (refcat == GSC || refcat == GSCACT) {
                                sprintf (temp, " %4d %4d %2d", gc[i], band, ngsc);
                                strcat (headline, temp);
                            }
                            else if (refcat == GSC2) {
                                sprintf (temp, " %2d  ", gc[i]);
                                strcat (headline, temp);
                            }
                            else if (typecol == 1) {
                                sprintf (temp, "  %2s", isp);
                                strcat (headline, temp);
                            }
                            else if (refcat == TMXSC) {
                                sprintf (temp," %6.1f", ((double)gc[i])* 0.1);
                                strcat (headline, temp);
                            }
                            else if (gcset && refcat != UCAC2 && refcat != UCAC3) {
                                sprintf (temp," %7d",gc[i]);
                                strcat (headline, temp);
                            }
                            if (mprop == 2) {
                                if (printepoch) {
                                    dtemp = DS_DateString (epoch, tabout);
                                    strcat (headline, dtemp);
                                    free (dtemp);
                                }
                                sprintf (temp, " %9.2f", gm[nmagr][i]);
                                strcat (headline, temp);
                            }
                            else if (printepoch) {
                                dtemp = DS_DateString (epoch, tabout);
                                strcat (headline, dtemp);
                                free (dtemp);
                            }
                            /* if (mprop == 1) {
                                 sprintf (temp, " %6.1f %6.1f", pra, pdec);
                                 strcat (headline, temp);
                                 } */

                            if (refcat == UCAC2 || refcat == UCAC3) {
                                sprintf (temp, "  %2d   %2d", nim, nct);
                                strcat (headline, temp);
                            }

                            /* Add distance from search center */
                            if ((ranges == NULL) && (printdist == 1))
                                {
                                    if (gdmax < 100.0)
                                        sprintf (temp, "  %5.2f", gdist);
                                    else if (gdmax < 1000.0)
                                        sprintf (temp, "  %6.2f", gdist);
                                    else if (gdmax < 10000.0)
                                        sprintf (temp, "  %7.2f", gdist);
                                    else if (gdmax < 100000.0)
                                        sprintf (temp, "  %8.2f", gdist);
                                    else
                                        sprintf (temp, "  %.2f", gdist);
                                    strcat (headline, temp);
                                }

                            /* Add specified keyword or object name */
                            if (refcat == TABCAT && keyword != NULL) {
                                if (gobj[i] != NULL)
                                    sprintf (temp, "  %s", gobj[i]);
                                else
                                    sprintf (temp, "  ___");
                                strcat (headline, temp);
                            }
                            else if ((refcat == BINCAT || refcat == TXTCAT || refcat == SCANCAT) &&
                                             gobj1 != NULL && gobj[i] != NULL)

                                {
                                    if (starcat[icat] == NULL || starcat[icat]->stnum > 0) {
                                        sprintf (temp, "  %s", gobj[i]);
                                        strcat (headline, temp);
                                    }
                                }

                            /* Add number of original entries in merged output */
                            if (catsort == SORT_MERGE) {
                                sprintf (temp, "  %2d", (int)gx[i]);
                                strcat (headline, temp);
                            }
              /* If necessary, truncate to the first column */
              if (getSingleColumnFlag()) {
                char *charPtr = headline;
                char charVal = *charPtr;
                while(charVal != 0) {
                  if (charVal == ' ') {
                    charPtr++;
                    charVal = *charPtr;
                  } else {
                    break;
                  }
                }
                while(charVal != 0) {
                  if (charVal != ' ') {
                    charPtr++;
                    charVal = *charPtr;
                  } else {
                    break;
                  }
                }
                *charPtr = 0;
              }

                            /* Write to file or standard output */
                            if (wfile) {
                                fprintf (fd, "%s\n", headline);
              } else {
                                printf ("%s\n", headline);
              }
                        }
                    }
      }
    }

    /* If searching more than one catalog, separate them with blank line */
    if (ncat > 0 && icat < ncat-1)
      printf ("\n");

    /* Free memory used for object names in current catalog */
    if (gobj1 != NULL) {
      for (i = 0; i < ns; i++) {
                if (gobj[i] != NULL) free (gobj[i]);
                gobj[i] = NULL;
      }
    }
  }

  /* Close output file */
  if (wfile)
    fclose (fd);

  return (ns);
}
/* Get a center and radius for a search area.  If the image center is not
 * given in the system of the reference catalog, convert it.
 * Return 0 if OK, else -1
 */

static int
GetArea (verbose,syscoor,eqcoor,sysout,eqout,epout,cra,cdec,dra,ddec,drad,dradi,
                 crao,cdeco)

         int	verbose;	/* Extra printing if =1 */
         int	syscoor;	/* Coordinate system of input search coordinates */
         double	eqcoor;		/* Equinox in years of input coordinates */
         int	sysout;		/* Coordinate system of output coordinates */
         double	eqout;		/* Equinox in years of output coordinates */
         double	epout;		/* Epoch in years of output coordinates (0=eqcoor */
         double	*cra;		/* Search center longitude/right ascension (degrees returned)*/
         double	*cdec;		/* Search center latitude/declination (degrees returned) */
         double	*dra;		/* Longitude/RA half-width (degrees returned) */
         double	*ddec;		/* Latitude/Declination half-width (degrees returned) */
         double	*drad;		/* Radius to search in degrees (0=box) (returned) */
         double	*dradi;		/* Radius of inner edge of search annulus (returned) */
         double	*crao;		/* Output search center longitude/right ascension */
         double	*cdeco;		/* Output search center latitude/declination */
{
    char rstr[32], dstr[32], cstr[32], dstri[32], dstro[32];

    *cra = ra0;
    *cdec = dec0;
    if (verbose) {
        if (syscoor == WCS_XY) {
        num2str (rstr, *cra, 10, 5);
            num2str (dstr, *cdec, 10, 5);
        }
        else if (syscoor==WCS_ECLIPTIC || syscoor==WCS_GALACTIC || degout0) {
        deg2str (rstr, 32, *cra, nddeg);
            deg2str (dstr, 32, *cdec, nddeg);
        }
        else {
        ra2str (rstr, 32, *cra, ndra);
            dec2str (dstr, 32, *cdec, nddec);
        }
        wcscstr (cstr, syscoor, 0.0, 0.0);
        fprintf (stderr,"Center:  %s   %s %s\n", rstr, dstr, cstr);
    }

    *crao = *cra;
    *cdeco = *cdec;
    if (syscoor != sysout || eqcoor != eqout) {
        wcscon (syscoor, sysout, eqcoor, eqout, crao, cdeco, epout);
        if (verbose) {
        if (syscoor == WCS_ECLIPTIC || syscoor == WCS_GALACTIC || degout0) {
                deg2str (rstr, 32, *crao, nddeg);
                deg2str (dstr, 32, *cdeco, nddeg);
            }
        else {
                ra2str (rstr, 32, *crao, ndra);
                dec2str (dstr, 32, *cdeco, nddec);
            }
        wcscstr (cstr, sysout, 0.0, 0.0);
        fprintf (stderr,"Center:  %s   %s %s\n", rstr, dstr, cstr);
        }
    }

    /* Set search box radius from command line, if it is there */
    if (dra0 > 0.0) {
        *drad = 0.0;
        if (syscoor == WCS_XY) {
        *dra = dra0;
        *ddec = ddec0;
        }
        else {
        *ddec = ddec0 / 3600.0;
        if (*cdec < 90.0 && *cdec > -90.0) {
                if (rdra)
                    *dra = dra0 / 3600.0;
                else
                    *dra = (dra0 / 3600.0) / cos (degrad (*cdec));
            }
        else
                *dra = 180.0;
        }
    }

    /* Search box */
    else if (rad0 < 0.0) {
        *drad = 0.0;
        if (syscoor == WCS_XY) {
        *dra = -rad0;
        *ddec = -rad0;
        }
        else {
        *ddec = -rad0 / 3600.0;
        *dra = *ddec / cos (degrad (*cdec));
        }
    }

    /* Search circle */
    else if (rad0 > 0.0) {
        if (syscoor == WCS_XY)
        *drad = rad0;
        else
        *drad = rad0 / 3600.0;
        *dra = *drad / cos (degrad (*cdec));
        *ddec = *drad;
        if (rad1 > 0.0)
        *dradi = rad1 / 3600.0;
        else
        *dradi = 0.0;
    }
    else {
        if (verbose)
        fprintf (stderr, "GetArea: Illegal radius, rad= %.5f\n",rad0);
        return (-1);
    }

    if (verbose) {
        if (syscoor == WCS_XY) {
        num2str (rstr, *dra, 10, 5);
        num2str (dstr, *ddec, 10, 5);
        num2str (dstro, *drad, 10, 5);
        num2str (dstri, *dradi, 10, 5);
        }
        else if (degout0) {
        deg2str (rstr, 32, *dra, 6);
        deg2str (dstr, 32, *ddec, 6);
        deg2str (dstro, 32, *drad, 6);
        deg2str (dstri, 32, *dradi, 6);
        }
        else {
        dec2str (rstr, 32, *dra, 2);
        dec2str (dstr, 32, *ddec, 2);
        dec2str (dstro, 32, *drad, 2);
        dec2str (dstri, 32, *dradi, 2);
        }
        if (*drad == 0.0)
        fprintf (stderr,"Area:    %s x %s\n", rstr, dstr);
        else if (*dradi > 0.0)
        fprintf (stderr, "Radius: %s-%s\n", dstri, dstro);
        else
        fprintf (stderr, "Radius: %s\n", dstro);
    }

    return (0);
}


static void
SearchHead (icat,sys,eq,ep,cra,cdec,dra,ddec,drad,dradi,nnfld,degout)

         int	icat;		/* Number of catalog in list */
         int	sys;		/* Coordinate system */
         double	eq;		/* Equinox */
         double	ep;		/* Epoch */
         double	cra, cdec;	/* Center coordinates of search box in degrees */
         double	dra, ddec;	/* Width and height of search box in degrees */
         double	drad;		/* Radius of search region in degrees */
         double	dradi;		/* Inner edge of annulus in degrees (ignore if 0) */
         int	nnfld;		/* Number of characters in ID field */
         int	degout;		/* Ouput in degrees if 1 */
{
    char rastr[32];
    char decstr[32];
    char cstr[16];
    char oform[16];
    char *catname;

    if (sys == WCS_XY) {
        num2str (rastr, cra, 10, 5);
        num2str (decstr, cdec, 10, 5);
    }
    else if (sys == WCS_ECLIPTIC || sys == WCS_GALACTIC || degout) {
        deg2str (rastr, 32, cra, nddeg);
        deg2str (decstr, 32, cdec, nddeg);
    }
    else {
        ra2str (rastr, 32, cra, ndra);
        dec2str (decstr, 32, cdec, nddec);
    }

    /* Set type of catalog being searched */
    if (!(icat == SCANCAT)) {
        fprintf (stderr,"ListCat: Catalog '%s' is missing\n",refcatname[icat]);
        return;
    }

    /* Label search center */
    sprintf (oform, "%%%ds", nnfld);
    if (objname)
        printf (oform, objname);
    else {
        catname = DS_CatName (icat, refcatname[icat]);
        printf (oform, catname);
    }
    wcscstr (cstr, sys, eq, ep);
    printf (" %s %s %s", rastr, decstr, cstr);
    if (degout) {
        if (drad != 0.0) {
        if (dradi > 0.0)
                printf (" r= %.2f - %.2f", dradi,drad);
        else
                printf (" r= %.2f", drad);
        }
        else
        printf (" +- %.2f %.2f", dra, ddec);
    }
    else {
        if (drad != 0.0) {
        dec2str (decstr, 32, drad, 1);
        if (dradi > 0.0) {
                dec2str (rastr, 32, dradi, 1);
                printf (" r= %s - %s", rastr, decstr);
            }
        else
                printf (" r= %s", decstr);
        }
        else {
        dec2str (rastr, 32, dra, 1);
        dec2str (decstr, 32, ddec, 1);
        printf (" +- %s %s", rastr, decstr);
        }
    }
    if (classd == 0)
        printf (" stars");
    else if (classd == 3)
        printf (" nonstars");
    if (ep != 0.0)
        printf (" at epoch %9.4f\n", ep);
    else
        printf ("\n");
    return;
}


void
PrintNum (maxnum, num, ndec)

         double	maxnum;	/* Maximum value of number to print */
         double	num;	/* Number to print */
         int	ndec;	/* Number of decimal places in output */
{
    char nform[120];
    int LenNum();

    if (ndec > 0)
        sprintf (nform, "%%%d.%df", LenNum (maxnum,ndec), ndec);
    else
        sprintf (nform, "%%%dd", LenNum (maxnum,ndec));
    printf (nform, num);
}


int
LenNum (maxnum, ndec)

         double	maxnum;	/* Maximum value of number to print */
         int	ndec;	/* Number of decimal places in output */
{
    if (ndec <= 0)
        ndec = -1;
    if (maxnum < 9.999)
        return (ndec + 3);
    else if (maxnum < 99.999)
        return (ndec + 4);
    else if (maxnum < 999.999)
        return (ndec + 5);
    else if (maxnum < 9999.999)
        return (ndec + 6);
    else
        return (ndec + 7);
}



int main(int argc,char *argv[])
{
  int ndcat;
  char rastr[32];
  char decstr[32];
  int lrange;
  char *newranges;
  int systemp = 0;		/* Input search coordinate system */
  char *rstr, *dstr, *astr, *cstr;
  char *str;
  char *str1;
  int lcat;
  char cs, cs1;
  char *refcatn;
  int nmag1;
  char  *ccom,*cep2;
  char coorout[32];
  char errmsg[256];
  int i;
  int imag;

  progname = DS_ProgName (argv[0]);

  ranges = NULL;
  keyword = NULL;
  objname = NULL;
  voerror[0] = (char) 0;
  for (i = 0; i < 5; i++)
    starcat[i] = NULL;

  /* Null out buffers before starting */
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
  DS_setrevmsg (RevMsg);
  coorout[0] = (char) 0;
  strcpy(cpname,"SCAT");
  refcatn = 0;
  ndcat = -1;


  /* crack arguments */
  for (argv++; --argc > 0; argv++) {



    /* Set parameters from keyword=value arguments */
    if (strchr (*argv, '=')) {
      if (scatparm (*argv))
                fprintf (stderr, "SCAT: %s is not a parameter.\n", *argv);
      refcat = SCANCAT;
      if (nmag > nmagmax)
                nmagmax = nmag;
      ndcat = DS_CatNdec (refcat);
    }

    /* Set search RA, Dec, and equinox if 2MASS ID */
    else if (DS_tmcid (*argv, &ra0, &dec0)) {
      syscoor = WCS_J2000;
      eqcoor = 2000.0;
      if (epoch0 == 0.0)
                epoch0 = 2000.0;
      if (eqout == 0.0)
                eqout = 2000.0;
      idrun = 0;
      ListCat (ranges, eqout);
      idrun = 1;
      ra0 = 0.0;
      dec0 = 0.0;
    }

    /* Set search RA, Dec, and equinox if colon in argument */
    else if (strsrch (*argv,":") != NULL) {
      if (argc < 2)
                PrintUsage (*argv);
      else {
                if (strlen (*argv) < 32)
                    strcpy (rastr, *argv);
                else {
                    strncpy (rastr, *argv, 31);
                    rastr[31] = (char) 0;
                }
                argc--;
                argv++;
                if (strlen (*argv) < 32)
                    strcpy (decstr, *argv);
                else {
                    strncpy (decstr, *argv, 31);
                    decstr[31] = (char) 0;
                }
                argc--;
                ra0 = str2ra (rastr);
                dec0 = str2dec (decstr);
                if (argc < 1) {
                    syscoor = WCS_J2000;
                    eqcoor = 2000.0;
                }
                else if ((syscoor = wcscsys (*(argv+1))) >= 0) {
                    coorsys = *++argv;
                    eqcoor = wcsceq (coorsys);
                }
                else {
                    syscoor = WCS_J2000;
                    eqcoor = 2000.0;
                }
      }
    }

    /* Set range and make a list of star numbers from it */
    else if (DS_isrange (*argv)) {
      if (ranges) {
                lrange = strlen(ranges) + strlen(*argv) + 16;
                newranges = (char *) calloc (lrange, 1);
                strcpy (newranges, ranges);
                strcat (newranges, ",");
                strcat (newranges, *argv);
                free (ranges);
                ranges = newranges;
                newranges = NULL;
      }
      else {
                lrange = strlen(*argv) + 16;
                ranges = (char *) calloc (lrange, 1);
                if (strchr (*argv,'.'))
                    match = 1;
                strcpy (ranges, *argv);
      }
      if (strchr (*argv, 'x') == NULL && ndcat > 0) {
                int n = ndcat;
                strcat (ranges, "x0.");
                while (n-- > 0)
                    strcat (ranges, "0");
                strcat (ranges, "1");
      }
    }

    /* Set decimal degree center or star number */
    else if (isnum (*argv)) {

      if (argc > 1 && isnum (*(argv+1))) {
                int ndec1, ndec2;

                /* Check for second number and coordinate system */
                if (argc > 2 && (systemp = wcscsys (*(argv + 2))) > 0) {
                    rstr = *argv++;
                    argc--;
                    dstr = *argv++;
                    argc--;
                    cstr = *argv;
                    eqcoor = wcsceq (cstr);
                }

                /* Check for two numbers which aren't catalog numbers */
                else {
                    ndec1 = DS_StrNdec (*argv);
                    ndec2 = DS_StrNdec (*(argv+1));
                    if (ndcat > -1 && ndec1 == ndcat && ndec2 == ndcat)
                        cstr = NULL;
                    else {
                        rstr = *argv++;
                        argc--;
                        dstr = *argv;
                        cstr = (char *) malloc (8);
                        strcpy (cstr, "J2000");
                        systemp = WCS_J2000;
                    }
                }
      }
      else
                cstr = NULL;

      /* Set decimal degree center */
      if (cstr != NULL) {
                ra0 = atof (rstr);
                dec0 = atof (dstr);
                syscoor = systemp;
                eqcoor = wcsceq (cstr);
      }

      /* Assume number to be star number if no coordinate system */
      else {
                if (strchr (*argv,'.'))
                    match = 1;
                if (ranges) {
                    lrange = strlen(ranges) + strlen(*argv) + 2;
                    newranges = (char *)calloc (lrange, 1);
                    strcpy (newranges, ranges);
                    strcat (newranges, ",");
                    strcat (newranges, *argv);
                    free (ranges);
                    ranges = newranges;
                    newranges = NULL;
                }
                else {
                    lrange = strlen(*argv) + 2;
                    ranges = (char *) calloc (lrange, 1);
                    strcpy (ranges, *argv);
                }
      }
    }

    else if (*(str = *argv) == '@') {
      readlist = 1;
      listfile = *argv + 1;
    }

    /* Otherwise, read command */
    else if ((*(str = *argv)) == '-') {
      char c;
      while ((c = *++str))
                switch (c) {

                case 'v':	/* more verbosity */
                    if (debug) {
                        webdump++;
                        debug = 0;
                        verbose = 0;
                    }
                    else if (verbose)
                        debug++;
                    else
                        verbose++;
                    printprog++;
                    break;

                case 'a':	/* Get closest source */
                    catsort = SORT_DIST;
                    closest++;
                    break;

                case 'b':	/* output coordinates in B1950 */
                    sysout0 = WCS_B1950;
                    eqout = 1950.0;
                    break;

                case 'c':       /* Set reference catalog */
                    if (argc < 2)
                        PrintUsage (str);
                    lcat = strlen (*++argv);
                    refcatn = (char *) calloc (1, lcat + 2);
                    strcpy (refcatn, *argv);
                    refcatname[ncat] = refcatn;
                    refcat = DS_RefCat (refcatn,title,&sysref,&eqref,&epref,&mprop,&nmag);
                    if (refcat == UCAC2 || refcat == UCAC3)
                        nmag1 = nmag + 4;
                    else if (refcat)
                        nmag1 = nmag;
                    else
                        PrintUsage ("Cannot find catalog");
                    if (nmag1 > nmagmax)
                        nmagmax = nmag1;
                    ndcat = DS_CatNdec (refcat);
                    ncat = ncat + 1;
                    argc--;
                    break;

                case 'd':	/* output in degrees instead of sexagesimal */
                    degout0++;
                    break;

                case 'e':	/* Set ecliptic coordinate output */
                    sysout0 = WCS_ECLIPTIC;
                    break;

                case 'f':	/* output in centering coordinates */
                    searchcenter++;
                    break;

                case 'g':	/* Set galactic coordinate output and optional center */
                    sysout0 = WCS_GALACTIC;
                    break;

                case 'h':	/* output descriptive header */
                    if (tabout)
                        printabhead = 0;
                    else {
                        printhead++;
                        printprog++;
                    }
                    break;

                case 'i':	/* ouput catalog object name instead of number */
                    printobj++;
                    if (!(*(str+1)) && argc > 1 && isnum (*(argv+1))) {
                        lofld = atoi (*++argv);
                        argc--;
                    }
                    break;

                case 'j':	/* center coordinates on command line in J2000 */
                    sysout0 = WCS_J2000;
                    eqout = 2000.0;
                    break;

                case 'k':	/* Keyword (column) to add to output from tab table */
                    if (argc < 2)
                        PrintUsage (str);
                    keyword = *++argv;
                    settabkey (keyword);
                    argc--;
                    break;

                case 'l':	/* Print center and closest star on one line */
                    oneline++;
                    catsort = SORT_DIST;
                    closest++;
                    nstars = 1;
                    break;

                case 'm':	/* Magnitude limit */
                    if (argc < 2)
                        PrintUsage (str);
                    cs1 = *(str+1);
                    if (cs1 != (char) 0) {
                        ++str;
                        if (cs1 > '9')
                            sortmag = (int) cs1;
                        else
                            sortmag = (int) cs1 - 48;
                    }
                    argv++;
                    argc--;
                    if ((ccom = strchr (*argv, ',')) != NULL) {
                        *ccom = (char) 0;
                        maglim1 = atof (*argv);
                        maglim2 = atof (ccom+1);
                    }
                    else {
                        maglim2 = atof (*argv);
                        if (argc > 1 && isnum (*(argv+1))) {
                            argv++;
                            argc--;
                            maglim1 = maglim2;
                            maglim2 = atof (*argv);
                        }
                        else if (MAGLIM1 == MAGLIM2)
                            maglim1 = -2.0;
                    }
                    break;

                case 'n':	/* Number of brightest stars to read */
                    if (argc < 2)
                        PrintUsage (str);
                    nstars = atoi (*++argv);
                    argc--;
                    break;

                case 'o':	/* Object name */
                    if (argc < 2)
                        PrintUsage (str);
                    objname = *++argv;
                    argc--;
                    break;

                case 'p':	/* Sort by distance from center */
                    catsort = SORT_DIST;
                    printprog++;
                    break;

                case 'q':	/* Output equinox in years */
                    if (argc < 2)
                        PrintUsage (str);
                    strcpy (coorout, *++argv);
                    if (coorout[0] == 'J' || coorout[0] == 'j')
                        sysout0 = WCS_J2000;
                    else if (coorout[0] == 'B' || coorout[0] == 'b')
                        sysout0 = WCS_B1950;
                    eqout = wcsceq (coorout);
                    argc--;
                    break;

                case 'r':	/* Search box or circle half-size in arcseconds */
                    if (argc < 2)
                        PrintUsage (str);
                    cs1 = *(str+1);
                    if (cs1 != (char) 0)
                        ++str;
                    argv++;

                    /* Separate 2 arguments for rectangles */
                    if ((dstr = strchr (*argv, ',')) != NULL) {
                        *dstr = (char) 0;
                        dstr++;
                    }

                    /* Separate 2 arguments for annulus */
                    if ((astr = strchr ((*argv)+1, '-')) != NULL) {
                        *astr = (char) 0;
                        astr++;
                    }

                    /* Convert radius or first argument to arcseconds */
                    if (strchr (*argv,':'))
                        rad0 = 3600.0 * str2dec (*argv);
                    else
                        rad0 = atof (*argv);

                    /* Convert outer radius of annulus to arcseconds */
                    if (astr != NULL) {
                        rad1 = rad0;
                        if (strchr (astr,':'))
                            rad0 = 3600.0 * str2dec (astr);
                        else
                            rad0 = atof (astr);
                    }

                    /* Convert second argument (=dec radius) to arcseconds */
                    if (dstr != NULL) {
                        if (cs1 == 'r')
                            rdra = 1;
                        else
                            rdra = 0;
                        dra0 = rad0;
                        rad0 = 0.0;
                        if (strchr (dstr, ':'))
                            ddec0 = 3600.0 * str2dec (dstr);
                        else
                            ddec0 = atof (dstr);
                        if (ddec0 <= 0.0)
                            ddec0 = dra0;
                    }
                    argc--;
                    break;

                case 's':	/* sort by RA, Dec, magnitude or nothing */
                    catsort = SORT_RA;
                    if (argc > 1) {
                        str1 = *(argv + 1);
                        cs = str1[0];
                        if (strchr ("ademinprs",(int)cs)) {
                            cs1 = str1[1];
                            argv++;
                            argc--;
                        }
                        else
                            cs = 'r';
                    }
                    else
                        cs = 'r';
                    switch (cs) {

                        /* Merge */
                    case 'e':
                        catsort = SORT_MERGE;
                        break;

                        /* Declination */
                    case 'd':
                        catsort = SORT_DEC;
                        break;

                        /* ID Number */
                    case 'i':
                        catsort = SORT_ID;
                        break;

                        /* Magnitude (brightest first) */
                    case 'm':
                        catsort = SORT_MAG;
                        if (cs1 != (char) 0) {
                            if (cs1 > '9')
                                sortmag = (int) cs1;
                            else
                                sortmag = (int) cs1 - 48;
                        }
                        break;

                        /* No sorting */
                    case 'n':
                        catsort = SORT_NONE;
                        break;

                        /* Distance from search center (closest first) */
                    case 'a':
                    case 'p':
                    case 's':
                        catsort = SORT_DIST;
                        break;

                        /* Right ascension */
                    case 'r':
                        catsort = SORT_RA;
                        break;
                    default:
                        catsort = SORT_RA;
                    }
                    break;

                case 't':	/* tab table to stdout */
                    if (tabout) {
                        tabout = 0;
                        votab = 1;
                        degout0 = 1;
                    }
                    else
                        tabout = 1;
                    if (printhead) {
                        printabhead = 0;
                        printhead = 0;
                        printprog = 0;
                    }
                    break;

                case 'u':       /* Print following 2 numbers at start of line */
                    if (argc > 2) {
                        printxy = 1;
                        xstr = *++argv;
                        argc--;
                        ystr = *++argv;
                        argc--;
                    }
                    break;

                case 'w':	/* write output file */
                    wfile++;
                    break;

                case 'x':       /* Guide Star object class */
                    if (argc < 2)
                        PrintUsage (str);
                    classd = (int) atof (*++argv);
                    setgsclass (classd);
                    argc--;
                    break;

                case 'y':	/* Set output coordinate epoch */
                    if (argc < 2)
                        PrintUsage (str);
                    if ((cep2 = strchr (*(argv+1), ','))) {
                        argv++;
                        *cep2 = (char) 0;
                        cep2 = cep2 + 1;
                        if (strchr (*argv, '.'))
                            epoch1 = atof (*argv);
                        else
                            epoch1 = fd2ep (*argv);
                        if (strchr (*argv, '.'))
                            epoch2 = atof (cep2);
                        else
                            epoch2 = fd2ep (cep2);
                        argc--;
                    }
                    else if (strchr (*(argv+1), '.')) {
                        argv++;
                        epoch0 = atof (*argv);
                        argc--;
                    }
                    else {
                        argv++;
                        epoch0 = fd2ep (*argv);
                        argc--;
                    }
                    break;

                case 'z':	/* Set append flag */
                    afile++;
                    wfile++;
                    break;

                case 'C':
                    if (argc < 2)
                        PrintUsage (str);
                    setPlateClass(*++argv);
                    argc--;
                    break;

                case 'H':
                    if (argc < 2)
                        PrintUsage (str);
                    // this no longer does anything: scanread() now looks at the
                    // standard environment variables directly.
                    argv++;
                    argc--;
                    break;

                case 'M':
                    setMosaicFlag(1);
                    break;

                case 'N':
                    if (argc < 2)
                        PrintUsage (str);
                    setPlateNumber(*++argv);
                    argc--;
                    break;

                case 'L':
                    setSingleColumnFlag(1);
                    break;

                case 'O':
                    setPessimisticFlag(1);
                    break;

                case 'B':
                    setPessimisticFlag(2);
                    break;


                case 'P':
                    if (argc < 1)
                        PrintUsage (str);
                    // this no longer does anything: scanread() now looks at the
                    // standard environment variables directly.
                    argv++;
                    argc--;
                    break;


                case 'S':
                    if (argc < 2)
                        PrintUsage (str);
                    setPlateSeries(*++argv);
                    argc--;
                    break;

                case 'T':
                    if (argc < 2)
                        PrintUsage (str);
                    setSearchType(*++argv);
                    argc--;
                    break;

                case 'U':
                    if (argc < 2)
                        PrintUsage (str);
                    // this no longer does anything: scanread() now looks at the
                    // standard environment variables directly.
                    argv++;
                    argc--;
                    break;


                case 'W':
                    if (argc < 7)
                        PrintUsage (str);
                    setWebTable(*++argv);
                    argc--;
                    setWebExtract(*++argv);
                    argc--;
                    setTarballExtract(*++argv);
                    argc--;
                    setWebFileType(*++argv);
                    argc--;
                    setWebArchiveType(*++argv);
                    argc--;
                    setRegionFlag(*++argv);
                    argc--;

                    webtable = 1;
                    break;



                default:
                    sprintf (errmsg, "* Illegal command -%s-", *argv);
                    PrintUsage (errmsg);
                    break;
                }
    }
    else {
      lcat = strlen (*argv);
      refcatn = (char *) calloc (1, lcat + 2);
      strcpy (refcatn, *argv);
      refcatname[ncat] = refcatn;
      ncat = ncat + 1;
    }
  }

  /* Set output equinox appropriately if output system is specified */
  if (eqout == 0.0) {
    if (eqcoor != 0.0)
      eqout = eqcoor;
    if (sysout0 != 0) {
      if (sysout0 == WCS_J2000)
                eqout = 2000.0;
      if (sysout0 == WCS_B1950)
                eqout = 1950.0;
    }
  }

  /* Set output epoch appropriately if output system is specified */
  if (sysout0 == 0) {
    if (strlen (coorout) > 0) {
      sysout0 = wcscsys (coorout);
      eqout = wcsceq (coorout);
    }
    else if (syscoor != 0)
      sysout0 = syscoor;
  }

  /* Set output equinox appropriately if output system is specified */
  if (eqout == 0.0) {
    if (eqcoor != 0.0)
      eqout = eqcoor;
    if (sysout0 != 0) {
      if (sysout0 == WCS_J2000)
                eqout = 2000.0;
      if (sysout0 == WCS_B1950)
                eqout = 1950.0;
    }
  }

  /* Set output epoch from output equinox if not otherwise set */
  if (epoch0 == 0.0)
    epoch0 = eqout;
  /* If http output, send header */
  if (votab) {
    printf ("Content-type: text/xml\n\n");
    printprog = 0;
  }
  else if (http) {
    printf ("Content-type: text/plain\n\n");
  }

  /* http is 0 */
  /* votab is 0 */
  /* readlist is 0 */

  if (sysout0 && !syscoor)
    syscoor = sysout0;
  if (syscoor) {
    if (eqout == 0.0) {
      if (eqcoor != 0.0)
                eqout = eqcoor;
      else if (syscoor == WCS_B1950)
                eqout = 1950.0;
      else
                eqout = 2000.0;
    }
    if (epoch0 == 0.0)
      epoch0 = eqout;
  }

  ListCat (ranges, eqout);

  for (i = 0; i < ncat; i++) {
    if (refcatname[i]) free (refcatname[i]);
    DS_ctgclose (starcat[i]);
  }

  /* Free memory used for search results and return */
  if (gx) {
    free ((char *)gx);
    gx = NULL;
  }
  if (gy) {
    free ((char *)gy);
    gy = NULL;
  }
  if (gm) {
    for (imag = 0; i < nmagmax; i++) {
      if (gm[imag]) {
                free ((char *)gm[imag]);
                gm[imag] = NULL;
      }
    }
    free ((char *)gm);
    gm = NULL;
  }
  if (gra) {
    free ((char *)gra);
    gra = NULL;
  }
  if (gdec) {
    free ((char *)gdec);
    gdec = NULL;
  }
  if (gnum) {
    free ((char *)gnum);
    gnum = NULL;
  }
  if (gc) {
    free ((char *)gc);
    gc = NULL;
  }
  if (gobj) {
    for (i = 0; i < nalloc; i++) {
      if (gobj[i] != NULL)
                free ((char *)gobj[i]);
    }
    free ((char *)gobj);
    gobj = NULL;
  }


  return(0);
}
