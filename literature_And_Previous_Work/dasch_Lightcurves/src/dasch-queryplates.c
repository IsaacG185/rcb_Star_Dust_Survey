// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// dasch-queryplates -- query the plate database
//
// Derived from dasch_scat.c, which is in turn derived from libwcs scat.c.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

// Starbase -- needed for File,Seek,etc.
#include <table.h>

// libwcs
#include <libwcs/wcs.h>
#include <libwcs/lwcs.h>
#include <libwcs/fitsfile.h>
#include <libwcs/wcscat.h>

// local
#include "pipelineutils.h"
#include "photometryutils.h"
#include "searchgsc.h"
#include "scanread.h"

static const int syscoor = WCS_J2000; // Input search coordinate system
static const double eqcoor = 2000.0; // Equinox of search center
static const double epoch0 = 2000.0; // Epoch for coordinates
static const int sysout0 = WCS_J2000; // Output coordinate system
static const double eqout = 2000.0; // Equinox for output coordinates

static double epoch2 = 0.0;	/* Latest epoch for catalog search */
static double epoch1 = 0.0;	/* Earliest epoch for catalog search */
static double ra0 = -99.0;	/* Initial center RA in degrees */
static double dec0 = -99.0;	/* Initial center Dec in degrees */
static double rad0 = 0.0;	/* Search box radius */
static struct StarCat *scancat = NULL; // scanner catalog data structure

extern void DS_setlimdeg (int degoutx); /* scanread.c */


static void
PrintUsage(char *command)
{
    if (command != NULL) {
        if (command[0] == '*')
            fprintf (stderr, "%s\n", command);
        else
            fprintf (stderr, "* Missing argument for command: %c\n", command[0]);
        exit(1);
    }

    fprintf(stderr, "Usage: dasch-queryplates [arguments]\n\n");
    fprintf(stderr, "  --help                           show this help\n");
    fprintf(stderr, "  -j <RA> <dec>                    position filter (J2000 coordinates)\n");
    fprintf(stderr, "  -r <arcsec>                      search radius\n");
    fprintf(stderr, "  -C class                         plate class filter\n");
    fprintf(stderr, "  -N number                        plate number filter\n");
    fprintf(stderr, "  -O                               \"pessimistic 1\": search inner bins of unfitted plates\n");
    fprintf(stderr, "  -B                               \"pessimistic 2\": search so that plate center is in bin 1\n");
    fprintf(stderr, "  -S series                        plate series filter\n");
    fprintf(stderr, "  -T [all|scanned|pending|wcsfit]  plate status filter (default is wcsfit)\n");

    if (command != NULL)
        exit(1);
    else
        exit(0);
}


static int
GetArea (
    int syscoor,   // Coordinate system of input search coordinates
    double eqcoor, // Equinox in years of input coordinates
    int sysout,    // Coordinate system of output coordinates
    double eqout,  // Equinox in years of output coordinates
    double epout,  // Epoch in years of output coordinates
    double *drad,  // Radius to search in degrees (returned)
    double *crao,  // Output search center longitude/right ascension (returned)
    double *cdeco  // Output search center latitude/declination (returned)
) {
    *crao = ra0;
    *cdeco = dec0;

    if (syscoor != sysout || eqcoor != eqout) {
        wcscon (syscoor, sysout, eqcoor, eqout, crao, cdeco, epout);
    }

    if (rad0 > 0.0) {
        *drad = rad0 / 3600.0;
    } else {
        return -1;
    }

    return 0;
}


static int
ListCat(double eqout)
{
    double crao, cdeco;	/* Output center long/lat or RA/Dec in degrees */
    double drad = 0.0;
    int ngmax = 800000; /* maximum number of results to print */

    if (rad0 == 0.0) {
        rad0 = 10.0;
    }

    if (
        GetArea(
            syscoor,
            eqcoor,
            sysout0,
            eqout,
            epoch0,
            &drad,
            &crao,
            &cdeco
        )
    )
        return 0;

    setSearchParameters(epoch1, epoch2);

    return scanread(
        "scanner",
        crao,
        cdeco,
        drad,
        sysout0,
        eqout,
        epoch0,
        ngmax,
        &scancat,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        0 // verbosity level
    );
}


int
main(int argc, char *argv[])
{
    char *str;
    char errmsg[256];

    for (argv++; --argc > 0; argv++) {
        if (strcmp(*argv, "--help") == 0) {
            PrintUsage(NULL);
        } else if ((*(str = *argv)) == '-') {
            /* Otherwise, read command */
            char c;

            while ((c = *++str)) {
                switch (c) {

                case 'r':
                    // Search radius
                    if (argc < 2)
                        PrintUsage(str);

                    argv++;
                    rad0 = atof(*argv);
                    argc--;
                    break;

                case 'j':
                    // position filter
                    if (argc < 3)
                        PrintUsage(*argv);

                    ra0 = str2ra(argv[1]);
                    dec0 = str2dec(argv[2]);
                    argv += 2;
                    argc -= 2;
                    break;

                case 'C':
                    if (argc < 2)
                        PrintUsage (str);
                    setPlateClass(*++argv);
                    argc--;
                    break;

                case 'N':
                    if (argc < 2)
                        PrintUsage (str);
                    setPlateNumber(*++argv);
                    argc--;
                    break;

                case 'O':
                    setPessimisticFlag(1);
                    break;

                case 'B':
                    setPessimisticFlag(2);
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

                case 'h':
                case 'H':
                case '?':
                    PrintUsage(NULL);
                    break;

                default:
                    sprintf (errmsg, "* unhandled option argument in \"-%s\"", *argv);
                    PrintUsage (errmsg);
                    break;
                }
            }
        } else {
            sprintf(errmsg, "* unhandled argument \"%s\"", *argv);
            PrintUsage(errmsg);
            break;
        }
    }

    // Other setup

    DS_setlimdeg(1);
    dasch_scanread_set_query_tool_mode(1);

    ListCat(eqout);
    scanclose(scancat);
    return 0;
}
