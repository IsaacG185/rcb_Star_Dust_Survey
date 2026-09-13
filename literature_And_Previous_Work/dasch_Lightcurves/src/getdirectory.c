// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 * Given a filename, generate the mosaic directory where that file will be found.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mysql.h>

#include "scandb.h"
#include "pipelineutils.h"

int main(int argc, char *argv[])
{
    int res;
    MYSQL my_connection;
    MYSQL_RES *res_ptr;
    MYSQL_ROW sqlrow;
    char filename[MAX_FILENAME_LEN];
    char series[MAX_SERIES_STRING];
    int plateNumber;
    int mosaicNumber;
    int rotation;
    int binning;
    int diskLocation;
    char queryString[MAX_QUERY_STRING];
    int totalCount = 0;
    int foundEntry = 0;
    int nvals;
    int errorFlag = 0;
    char *argstr;
    char cmdchar;
    int fullMosaic = 0;
    int noWarning = 0;
    char *raid_override = NULL;

    filename[0] = 0;

    /* Check qualifiers */
    for (argv++; --argc > 0; argv++) {
        argstr = *argv;

        /* Decode arguments */
        if (argstr[0] != '-') {
            /* This must be the list of plates */
            if (strlen(filename) > 0) {
                errorFlag = 1;
                printf("ERROR: list %s is being overwritten by %s\n", filename, argstr);
            } else {
                strncpy(filename, argstr, MAX_FILENAME_LEN-1);
                filename[MAX_FILENAME_LEN-1] = 0;
            }
        } else {
            while ((cmdchar = *++argstr) != 0) {
                switch(cmdchar) {
                case 'm': /* provide the full mosaic name */
                case 'M':
                    fullMosaic = 1;
                    break;

                case 'q': /* Do not print the warning message */
                case 'Q':
                    noWarning = 1;
                    break;

                default:
                    printf("* illegal command -%c-", cmdchar);
                    errorFlag = 1;
                }
            }
        }
    }

    if (filename[0] == 0) {
        fprintf(stderr, "ERROR: No filename was specified\n");
        errorFlag = 1;
    }

    if (errorFlag) {
        fprintf(stderr, "Usage: getdirectory filename [-q] [-m]\n");
        fprintf(stderr, "       where -q avoids printing a warning if the mosaicNumber and binning is not specified\n");
        fprintf(stderr, "             -m provides a full path to the mosaic if it has been WCS and SCAMP fitted\n");
        return -1;
    }

    dasch_init_scandb(&my_connection);

    if (ParseFilename(filename, series, &plateNumber, &mosaicNumber, &binning, &rotation) == 0) {
        if (ParseFilename2(filename, series, &plateNumber) == 0) {
            fprintf(stderr, "ERROR: Failed to parse filename %s as %s%05d_%02d_%02d rotation %d\n",
                    filename, series, plateNumber, mosaicNumber, binning, rotation);
            return -1;
        } else {
            /* Here we have the series and plateNumber but not the mosaic number.  We set the binning to 1 */
            if (SelectBestMosaic(&my_connection, series, plateNumber, &mosaicNumber, &rotation) == 0) {
                fprintf(stderr, "ERROR: Failed to find a candidate mosaic for plate filename %s, series %s plateNumber %d\n",
                        filename, series, plateNumber);
                return -1;
            } else {
                binning = 1;

                if (noWarning == 0) {
                    fprintf(stderr, "WARNING: Selecting mosaic for %s as %s%05d_%02d_%02d rotation %d\n",
                            filename, series, plateNumber, mosaicNumber, binning, rotation);
                }
            }
        }
    }

    if (binning != 1) {
        fprintf(stderr, "File binning must be 1\n");
        return -1;
    }

    sprintf(queryString, "SELECT diskLocation FROM mosaics where series = '%s' and plateNumber = %d and mosaicNumber = %d and solutionNumber = 0;", series, plateNumber, mosaicNumber);

    if (strlen(queryString) > MAX_QUERY_STRING) {
      fprintf(stderr, "ERROR: MAX_QUERY_STRING exceeded %zu\n", strlen(queryString));
      return -1;
    }

    res = mysql_query(&my_connection, queryString);
    if (!res) {
        res_ptr = mysql_store_result(&my_connection);
        if (res_ptr) {
            while ((sqlrow = mysql_fetch_row(res_ptr)) != NULL) {
                totalCount++;
                if (sqlrow[0] == NULL) {
                    fprintf(stderr, "No disk location found for %s\n", filename);
                    continue;
                }

                nvals = sscanf(sqlrow[0], "%d", &diskLocation);
                if (nvals != 1) {
                    fprintf(stderr, "ERROR: nvals is %d for diskLocation\n", nvals);
                    continue;
                }

                foundEntry++;
            }

            mysql_free_result(res_ptr);
        }
    } else {
        fprintf(stderr, "Select error %d: %s res %d\n", mysql_errno(&my_connection), mysql_error(&my_connection), res);
    }

    mysql_close(&my_connection);

    if (foundEntry != 1) {
        fprintf(stderr, "ERROR: %d database entries found\n", foundEntry);
        return -1;
    }

    raid_override = getenv("DASCH_RAID_OVERRIDE");
    if (raid_override == NULL) {
        printf("/dasch/raid%03d/ExposureData/Mosaics/%s/%05d_%02d", diskLocation, series, plateNumber, mosaicNumber);
    } else {
        printf("%s/ExposureData/Mosaics/%s/%05d_%02d", raid_override, series, plateNumber, mosaicNumber);
    }

    if (fullMosaic != 0) {
        if (rotation == 0) {
            printf("/%s%05d_%02d_01ww_tnx.fit", series, plateNumber, mosaicNumber);
        } else {
            printf("/%s%05d_%02d_01r%dww_tnx.fit", series, plateNumber, mosaicNumber, rotation);
        }
    }

    printf("\n");
    return EXIT_SUCCESS;
}

/*
 * Jan 15, 2008 Edward J. Los.  Change /mosaic/ to /raidnnn/ to handle multiple mosaic disks
 * Sep 22, 2009 Edward J. Los   Correct for multiple mosaic solutions.
 * Jun 27, 2011 Edward J. Los - Add SelectBestMosaic to work with only a series and plateNumber
 */
