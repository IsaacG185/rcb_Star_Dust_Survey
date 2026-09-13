// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*  parseidtable.c
 *
 *    Performs additional parsing on transient candidate tables
 *
     gcc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I/usr/include/mysql  -I /dasch/install/include   parseidtable.c -L /dasch/install/lib -lm pipelineutils.a -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread   -ltable -lutil -lwcs  plottransientstub.o  -o parseidtable

 * 
 *  
    parseidtable -v  -i /home/scanner/backup/2015_09_17/candidates_unmatched_apass_2015_09_17.db -o /home/scanner/junk/los.db -c /home/scanner/Pipeline/flarecandidates.txt  -d /home/scanner/junk/candidates 2015-09-17T15-28-24 2015-09-18T05-40-43 2015_09_17 -f /home/scanner/Pipeline/flarecandidates.list -p /home/scanner/Pipeline/lospossible.log

    parseidtable -v  -i /home/scanner/backup/2015_09_17/candidates_apass_2015_09_17.db -o /home/scanner/junk/los.db -c /home/scanner/Pipeline/flarecandidates.txt  -d /home/scanner/junk/candidates 2015-09-17T15-28-24 2015-09-18T05-40-43 2015_09_17 -f /home/scanner/Pipeline/flarecandidates.list -p /home/scanner/Pipeline/lospossible.log

    parseidtable -v  -i /dasch/Pipeline/candidates/candidates_unmatched_apass_2015_09_21.db -o /home/scanner/junk/los.db -c /home/scanner/Pipeline/flarecandidates.txt  -d /dasch/Pipeline/candidates 2015-09-21T19-56-49 2015-09-21T20-24-16 2015_09_21x -f /home/scanner/Pipeline/flarecandidates.list -p /home/scanner/Pipeline/lospossible.log

    parseidtable -v  -i /dasch/Pipeline/candidates/candidates_apass_2015_09_21.db -o /home/scanner/junk/los.db -c /home/scanner/Pipeline/flarecandidates.txt  -d /dasch/Pipeline/candidates 2015-09-21T19-56-49 2015-09-21T20-24-16 2015_09_21x -f /home/scanner/Pipeline/flarecandidates.list -p /home/scanner/Pipeline/lospossible.log



    
    column -i los.db  nearbyREFflag releaseField peakDays peakYear peakMag peakSlope peakRMS peakCount peakUpperCount peakNumber peakMaxMag peakMinMag  REF nearbyObjects 

 * 
 * Jun 16, 2015  Edward J. Los - Original Version
 * Jun 23, 2015  Edward J. Los - add gmag, gminusrmag, and nearestREFarcsec to parse the SDSS catalog
 * Aug  4, 2015  Edward J. Los - convert to starbase
 *                               Check new candidates against a list of existing candidates
 * Aug 17, 2015  Edward J. Los - Add npoints, ngood, peakCount, and peakExcess = npoints-ngood; fix bug in ALREADY CHECKED candidate comment printout
 * Aug 19, 2015  Edward J. Los - Read peakExcess from the input table (now npoints - points outside the suspected transient)
 * Aug 21, 2015  Edward J. Los - Add peakMaxMag for brightest magnitude
 * Sep  7, 2015  Edward J. Los - Correct buffer overrun errors
 * Sep 15, 2015  Edward J. Los - Rename peakExcess to peakExtra
 *                               Add peakCountNF, peakCountWF, peakOutside, peakOutsideNF, and peakOutsideWF
 *                               Add peakRA and peakDec with 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window
 *                               Add peakDefectCount for the number of defects in the selected transient window.
 * Sep 19, 2015  Edward J. Los - Parse nearbyObjects to separate out variables, NGC objects, nearby stars, and SDSS stars
 *                               Add peakNearbyDistance for the distance to nearby stars that are within 2 magnitudes of the brightest or average magnitude of the transient candidate. 
 * Sep 27, 2015  Edward J. Los - Add exclusion zones
 * Sep 29, 2015  Edward J. Los - List contents of the base directory
 * Oct  4, 2015  Edward J. Los - Add peakMultipleCount for the number of multiple exposure plates in the selected transient window
 *                             - add input thumbnail file showing the directories where candidates appear
 *                             - add output best candidate filename
 * Oct  9, 2015  Edward J. Los - Add '-m' to print file moving commands
 *                               Add '-w' to print warnings of missing possible flare candidates
 *                               Add '-g' to study location groupings
 *                               Change '-s' to print summary table
 *                               Add 'multiple33' class for transients with more 1/3 points coming from multiple exposures
 *                               Add 'multiple99' class for transients with more than 90% of points coming from multiple exposures
 *                               Add 'location' class for transients with at least one nearby flare candidate
 *                               Add 'unreleased' class for transients in regions not yet released
 *                               Add 'extra' class with for transients where 'peakExtra' is more than 10 points
 *                               Add 'sdss' class for transients that match SDSS objects within 5 arcsec and dimmer than magnitude 17.0
 *                               Prepend the peak count when copying the file names
 *                               Combine matched and unmatched into a single run
 *                               Add '-q' for a table of percentages of multiple exposure points in a flare
 * Oct 14, 2015  Edward J. Los - Correct behavior if previous candidate table contains no entries
 * Oct 19, 2015  Edward J. Los - Allow all entries into the possible table that matches the date of the run
 *                               Add "flag" column to the possible table.
 * Oct 24, 2015  Edward J. Los - Add gsc_bin_index for the current gsc bin of the transient
 *                             - Add peakDradRMS2 for the dradRMS2 values added in quadrature
 *                             - Add peakExtra to the png filenames:  Count_<peakCount>_<peakExtra>_<REF>_<timestamp>.png
 * Oct 29, 2015  Edward J. Los - Add peakDradRMS3 for the drad rms of points in the selected 90 day window relative to peakRA and peakDec
 *                             - Make rundate the matched and unmatched timestamps
 *                             - Remove "Count_" prefix
 *                             - The possible file will now have the following header: RA_TC  Dec_TC  Sig_TC  Npts_TC  Npts_ext name VSXname
 * Nov  5, 2015  Edward J. Los - add drad column for "nearby" distance in arcsec
 *                             - make "nearby" the least fatal.
 * Jan  3, 2015  Edward J. Los - Josh's memo of 12/31/2015 10:33 AM.  
 *                               1.  Make sdss the highest priority flag (did not work!)
 *                               2.  Change g > 17.0 to g > 17.5
 *                               3.  Offset from DASCH source <= 2.5*rms instead of simply 5 arcsec.  Add in the distance between the peak location and the DASCH object center.
 * May 18, 2015  Edward J. Los - Add REDUCE_LOCATION_PRIORITY
 * May 26, 2015  Edward J. Los - Add KD tree support to improve performance. (reduced execution time from 2006 seconds to 149 seconds
 * May 27, 2015  Edward J. Los - move "location" to just above "nearby"
 * May 30, 2015  Edward J. Los - Add a histogram of "location" distances.
 * Jun  6, 2015  Edward J. Los - Reduce the maximum SDSS distance to 12.0
 * Dec  7, 2016  Edward J. Los - Introduce MIN_TC_DATERANGE, the minimum allowable time coincidence between to separate candiates -- 'sametime' flag
 * Dec 16, 2016  Edward J. Los - Ignore multiple exposure precentages if at least MIN_TC_POINTS are not multiples. (PREFILTER_MULTDEFECT)
 * Dec 24, 2016  Edward J. Los - Add the "colorterm" condition for uncalibrated images on red plates and set reject flag if there are more than 33% bad colorterm points
 * Dec 26, 2016  Edward J. Los - group candidates in common subdirectories
 *                               Change "cp" to "mv"
 *                               Change "NOT A PREVIOUS CANDIDATE" to "NOT ON PREVIOUS CANDIDATE LIST"
 *                               Make a separate directory "lowdrad" for peakDradRMS3 < MAX_TC_PEAKDRADRMS3
 * Dec 27, 2016  Edward J. Los - Correct a bug in parsing the thumbnail filename.
 *                               Add RUNDATE_ONLY to enforce thumbnails with the correct date
 *                               Return to the original 33% multiple exposure filter.
 *                               Make GOOD_MASK more inclusive
 * Dec 28, 2016  Edward J. Los - Correct conditional error for the new 33% multile exposure filter.
 * Jan 24, 2017  Edward J. Los - Add "-b" qualifier to select only the "best" objects defined by 
 *                               Add peakDradRMS3 to possible file, replacing drad (Josh memorandum of Mon, 23 Jan 2017 12:25:57)
 *                               Add peakYear (Josh memorandum of Mon, 23 Jan 2017 12:25:57)
 *                               Add peakMaxMag (Josh memorandum of Mon, 23 Jan 2017 10:29:23)
 * Feb  9, 2017  Edward J. Los - Implement "LOWPOINTS" for objects with only three detections
 * Feb 17, 2017  Edward J. Los - Add the "-t" qualifier to list details of the "sametime" matches.
 *                               Back out the "LOWPOINTS" flag because parallel changes in limiting magnitude filter make this flag unecessary.
 * Mar  6, 2017  Edward J. Los - For the '-m' copy commands, add the subdirectory to the first part of the "cp" command.
 *                               Create a "best" subdirectory for "contamination", "knownvariable",and "nobadflag"
 *                               Add a '-a' file which will generate a previous candidates file (-c) for the next run.
 * Mar 17, 2017  Edward J. Los - Eliminate duplicates in the input file.
 *                               Add support for peakLimitingYears and peakLimitingPoints
 * Mar 27, 2017  Edward J. Los - Add SDSS entries that set the REJECT_FLAG_SDSS to nearbyVSX
 *                               Guarantee that filenames in the possible table are prefixed by "/n"
 *                               Add peakLimitingYears to the possible table.  Remove Sig_TC, reorder the table.
 *                               Change the text for sdss matches from "less than" to "dimmer than"
 *                               Print actual parameters in the log table for each flag.
 *                               Clear bestFlag for the "extra" flag, put these entries in a separate "extra" subdirectory
 *                               Move "previous" to the top of the list, combine with Josh Grindlay's TC-best-candidates-master.xlsx of Sat, 25 Mar 2017 22:26:39
 * Nov  7, 2017  Edward J. Los - Reverse the order of the summary table presentation: most likely first
 *                               Correct bug in which REJECT_FLAG_BEST was being included in the "best" directory
 *                               Allow multiple33 in the best category and put results in the "contaminated" subdirectory
 *                               Add the following codes immediately after the '#' in the candidate file name 
 *                                    "O" interesting    "interest1"
 *                                    "Y" less interesting "interest2"
 *                                    "N" not a TC,      "nottc"
 *                                    " " no judgement  "previous"
 * Dec 15, 2017 Edward J. Los -  Split "interest" into "orange" and "yellow"
 *                               Put all of the sdss matches in the sdss subdirectory
 * Feb  2, 2018 Edward J. Los   Add LARGE_TC_MAGNITUDE = 0.5 for selection of largest flares and REJECT_FLAG_LOWRANGE
 *                              Reorganize the selection table for complete rejection of unwanted images
 *                               'previous' directory no longer part of the cuts
 *                               'rejected' directory for many more entries.
 *                              Remove the unused '-u' qualifier
 *                              Remove the unused '-b' qualifier
 *                              Make the '-c' qualifier apply to the previous run's candidates
 *                              Add the '-j' qualifier for Josh's list of interest.
 *                              Change REJECT_FLAG_INTEREST1 to REJECT_FLAG_ORANGE
 *                              Change REJECT_FLAG_INTEREST2 to REJECT_FLAG_YELLOW
 *                              Remove the "possible" functionality in favor of Josh's list
 *                              Change MAX_SDSS_DISTANCE from 12 to 10
 *                              Get rid of unmatchedTimeStamp because it equals matchedTimeStamp
 *                              Output to "new"      - all not in previous runs, regardless of other conditions
 *                                        "interest" - all  yellow and orange colored entries in Josh's list regardless of other conditions
 *                                        "best"     - all unflagged, unreleased, and novae, and range > 0.5 mag even if drad > 10 arcsec.
 *                                        "secondbest" - all nearby and lowrange objects unless drad < 10 arcsec, then they become "best".  Any SDSS matches < 10 arcsec also become "best"
 *                                        "rejected" - all "NOT a TC" in Josh's list
 *                              Add the above strings to the possible_name flags column
 *                              Add nearbyObjectDist, the closest distance to any catalog object to the output table
 *                              Add peakDays, the duration of the flare, to the output table
 * Feb  5, 2018 Edward J. Los   correct the logic to eliminate duplicates from id_name
 * Feb  8, 2018 Edward J. Los   create "nova","ugem", and "variable" subdirectories
 * Mar  6, 2018 Edward J. Los   Change MAX_GALACTIC_LATITUDE from 15 to 0 degrees 
 * Mar 13, 2018 Edward J. Los   If a "##" appears, in the <candidates_name> or <interesting_name> files, then this number is the current Npts_TC or peakCount
 *                              Write the peakCount to the <all_name> file.  Do not put the runDate string in this file
 *                              If peakCount increases from the last run and we have a best or second best object, then put it in the "morepoints" directory.
 *                              If the peakNoDefectMagColumn exists, use it to set REJECT_FLAG_LOWRANGE
 * Apr  6, 2018 Edward J. Los   Support the "peakLongOutburst" column with the flag REJECT_FLAG_LONGFLARE
 *                              Remove support for older id file formats
 *                              Remove "INDEX_ALL_PLATES"
 * May  2, 2018 Edward J. Los   Split LONG_TC_DAYS into
 *                                    LONG_TC_FLARE_DAYS 365.0 days the flare is active
 *                                    LONG_TC_PRE_DAYS (5*365.0) days integration before the transient
 *                                    LONG_TC_SKIP_DAYS (2*365.0) days to skip after finding a transient
 * May 29, 2018 Edward J. Los   Support the GAIA calibration
 * Dec  4, 2020 Edward J. Los   Combine the flarecandidates.txt table with DASCH-Nov2018-DynBHLMXBs-DutyCycles_v2.xlsx 
 *                               Merge KV UMa 
 * Dec  9, 2020 Edward J. Los   Add InitFlareCandidateList
 * 
 */ 
#include <math.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h> /* needed for exit */
#include <errno.h>
#include <ctype.h>
#include <assert.h>
#include "kdtree.h"
double fmax(double x, double y);

#if 0
#define LOS_DEBUG1 1  /* Look for a particular object */
#define LOS_REF "DASCH_J054201.4+123050"
#define LOS_REF_NUMBER 305420141123050L
#endif
/* #define LOS_DEBUG2 1 */ /* Debug the flarecandidates.txt table */



#define MAX_BUFFER 2048
#define MAX_INPUT_NAME 1024
#define MAX_CANDIDATE_SIZE  700
#define MAX_COMMENT_SIZE 1000
#define MAX_SUBDIR_NAME 10
#define MAX_SDSS_DISTANCE_FACTOR 2.5 /* Multiply by Sig_TC to get radius */
#define MAX_SDSS_DISTANCE 10.0 /* Limit to 12 arcsec */
#define MAX_SDSS_MAGNITUDE 17.5 /* brightest magnitude */
/* #define SKIP_PREVIOUS 1 *//* skip the check of previous objects */ 
/* #define PREFILTER_MULTDEFECT 1 *//* if defined, ignore all other conditions for multiple90, multiple33, multdefect, and defect */
/* #define RUNDATE_ONLY *//* Use only thumbnails with the desired run date */
/* #define SKIP_MULTDEFECT 1 */ /* skip multdefect */
#define MAX_PEAK_EXTRA_COUNT 10 /* Maximum points outside the selected flare */
#define MAX_GALACTIC_LATITUDE  0 /* Maximum galactic latitude */
#define MAX_LOCATION_HISTOGRAM_ARCSEC (5*3600)
#define LOCATION_HISTOGRAM_BINSIZE 60
#define MAX_LOCATION_SEPARATION_ARCSEC 1000
#define MAX_DEFECT_PERCENTAGE 34
#define MAX_COLORTERM_PERCENTAGE 34
#define MAX_MULTIPLE_PERCENTAGE_33 34
#define MAX_MULTIPLE_PERCENTAGE_90 91
#define MAX_NEARBY_ARCSEC 50 /* parseidtable Limit to reject all catalog stars within 2 magnitudes of the transient */
#define MAX_DISTANCE (32) /* Report already seen flare candidate object if less than this distance */
#define MAX_FILE_DISTANCE (14)                        /* Eliminate duplicates if less than half this distance */
#define MIN_VARIABLE_DISTANCE 18.0                      /* display all variables within this distance or 3*Sig_TC, whichever is larger */




typedef struct _thumbnail {
  long long REFNumber;
  double ra;
  double dec;
  int refType;
  int haveCoords;
  int candidate_index;
  double candidateDistance;
  char filename[MAX_INPUT_NAME+1];
  char name[MAX_CANDIDATE_SIZE+1];
  char subdir[MAX_SUBDIR_NAME+1];
} THUMBNAIL,*PTHUMBNAIL;


#define MAX_NEARBY_STRING 60
#define MAX_NEARBY_COUNT 4
#define MAX_IDENTIFIER_STRING 30
#define MAX_NEARBY_TYPE_STRING  25
#define NEARBY_TYPE_UNKNOWN 0
#define NEARBY_TYPE_CATALOG 1
#define NEARBY_TYPE_SDSS    2
#define NEARBY_TYPE_NGC     3
#define NEARBY_TYPE_GALAXY  4
#define NEARBY_TYPE_IC      5
#define NEARBY_TYPE_VSX     6


typedef struct _nearbyimage {
  int nearbyType;
  int asteriskFlag;
  int sdssrejectflag; /* This is the one which set REJECT_FLAG_SDSS */
  double distance; /* Distance in arcsec */
  char nearbyObject[MAX_NEARBY_STRING+1];
  char identifier[MAX_IDENTIFIER_STRING+1];
  char type[MAX_NEARBY_TYPE_STRING+1];
  double catalogmag;
  double color;
} NEARBYIMAGE,*PNEARBYIMAGE;

#define REJECT_FLAG_NGC                    1  /*  0 */
#define REJECT_FLAG_NEARBY                 2  /*  1 */
#define REJECT_FLAG_DEFECT                 4  /*  2 */
#define REJECT_FLAG_UG                     8  /*  3 */
#define REJECT_FLAG_MIRA                  16  /*  4 */
#define REJECT_FLAG_RR                    32  /*  5 */
#define REJECT_FLAG_NOVA                  64  /*  6 */
#define REJECT_FLAG_VAR                  128  /*  7 */
#define REJECT_FLAG_PREVIOUS             256  /*  8 */
#define REJECT_FLAG_EXCLUSION            512  /*  9 */
#define REJECT_FLAG_ECLIPSING           1024  /* 10 */
#define REJECT_FLAG_POLAR               2048  /* 11 */
#define REJECT_FLAG_LOCATION            8192  /* 12 */
#define REJECT_FLAG_MULTIPLE33          4096  /* 13 */
#define REJECT_FLAG_UNRELEASED         16384  /* 14 */
#define REJECT_FLAG_EXTRA              32768  /* 15 */
#define REJECT_FLAG_SDSS               65536  /* 16 */
#define REJECT_FLAG_MULTIPLE90        131072  /* 17 */
#define REJECT_FLAG_MULTDEFECT        262144  /* 18 */
#define REJECT_FLAG_SAMETIME          524288  /* 19 */
#define REJECT_FLAG_COLORTERM        1048576  /* 20 */
#define REJECT_FLAG_LOWPOINTS        2097152  /* 21 */
#define REJECT_FLAG_LIM_YEARS        4194304  /* 22 */
#define REJECT_FLAG_LIM_POINTS       8388608  /* 23 */
#define REJECT_FLAG_NOTTC           16777216  /* 24 */
#define REJECT_FLAG_ORANGE          33554432  /* 25 */
#define REJECT_FLAG_YELLOW          67108864  /* 26 */
#define REJECT_FLAG_LOWRANGE       134217728  /* 27 */
#define REJECT_FLAG_MOREPOINTS     268435456  /* 28 */
#define REJECT_FLAG_LONGFLARE      536870912  /* 29 */
#define REJECT_TABLE_SIZE 30
#if 0
#define REJECT_FLAG_NEW          REJECT_FLAG_NOTTC         /*   rejectTableIndex > newIndex  */
#else
#define REJECT_FLAG_NEW          REJECT_FLAG_MIRA         /*   rejectTableIndex > newIndex  */
#endif
#define REJECT_FLAG_SECONDBEST   REJECT_FLAG_ORANGE       /* bestIndex >= rejectTableIndex > secondBestIndex  */
#define REJECT_FLAG_BEST         REJECT_FLAG_LOWRANGE     /* unusedIndex >= rejectTableIndex > bestIndex */
#define REJECT_FLAG_UNUSED       REJECT_FLAG_UNRELEASED   /* > unusedIndex */
#ifdef SKIP_MULTDEFECT
#ifdef SKIP_PREVIOUS
#else /* SKIP_PREVIOUS */
#endif /* SKIP_PREVIOUS */
#else /* SKIP_MULTDEFECT */
#ifdef SKIP_PREVIOUS
#else /* SKIP_PREVIOUS */
#endif /* SKIP_PREVIOUS */
#endif /* SKIP_MULTDEFECT */
#define LOCATION_MASK (~(REJECT_FLAG_UNRELEASED|REJECT_FLAG_SDSS|REJECT_FLAG_MULTIPLE33|REJECT_FLAG_MULTIPLE90|REJECT_FLAG_DEFECT|REJECT_FLAG_MULTDEFECT|REJECT_FLAG_PREVIOUS|REJECT_FLAG_LOCATION)) /* used to find nearby candidates */
#define MULTIPLE_MASK (~(REJECT_FLAG_UNRELEASED|REJECT_FLAG_SDSS|REJECT_FLAG_MULTIPLE33|REJECT_FLAG_MULTIPLE90|REJECT_FLAG_PREVIOUS)) /* used to find multiple exposure ratios */
#define GOOD_MASK     (~(REJECT_FLAG_NEARBY|REJECT_FLAG_UG|REJECT_FLAG_MIRA|REJECT_FLAG_RR|REJECT_FLAG_NOVA|REJECT_FLAG_VAR|REJECT_FLAG_PREVIOUS|REJECT_FLAG_ECLIPSING|REJECT_FLAG_POLAR|REJECT_FLAG_UNRELEASED|REJECT_FLAG_EXTRA|REJECT_FLAG_SDSS))
#define PREFILTER_MULTDEFECT_MASK (REJECT_FLAG_MULTIPLE33|REJECT_FLAG_MULTIPLE90|REJECT_FLAG_MULTDEFECT|REJECT_FLAG_DEFECT)
#define MAX_WEBNAME 200
typedef struct _candidatetype {
  int flagbit;
  int flagmask;
  char *flagname; /* C Language bit definition */
  char *subdir; /* Target subdirectory */
  char *webnametemplate; /* Web friendly name template */
  char webname[MAX_WEBNAME];  /* Web-Friendly name */
} CANDIDATETYPE,*PCANDIDATETYPE;


CANDIDATETYPE flareRejectTable[] = {
  {24,REJECT_FLAG_NOTTC      ,"nottc"    ,"rejected"     ,"Not a TC",0},
  {9,REJECT_FLAG_EXCLUSION  ,"exclusion" ,"rejected"     ,"Located in exclusion zone",0},
  {0,REJECT_FLAG_NGC        ,"ngc"       ,"rejected"     ,"Too close to an NGC object",0},  /* Fatal Josh's memo of Tue, 2 Jan 2018 18:28:27 */
  {2,REJECT_FLAG_DEFECT     ,"defect"    ,"rejected"     ,"Has more than %d percent defects",0}, /* Fatal Josh's memo of Mon, 1 Jan 2018 13:43:46 14:34:20 22:55:07 */
  {17,REJECT_FLAG_MULTIPLE90,"multiple90","rejected"     ,"Has more than %d percent multiple exposure plates",0},
  {20,REJECT_FLAG_COLORTERM ,"colorterm" ,"rejected"     ,"Missing colorterm correction",0},
  {13,REJECT_FLAG_LOCATION  ,"location"  ,"rejected"     ,"In a region of the sky with at least one nearby flare candidate",0},  /* Fatal Josh's memo of Mon, 1 Jan 2018 13:43:46 14:34:20 22:55:07 */
  {19,REJECT_FLAG_SAMETIME  ,"sametime"  ,"rejected"     ,"Within %.0f hour of a different candidate",0},  /* Fatal Josh's memo of Mon, 1 Jan 2018 13:43:46 14:34:20 22:55:07 */
  {23,REJECT_FLAG_LIM_POINTS,"coveragept","rejected"     ,"Has fewer than %d plates of coverage",0},
  {22,REJECT_FLAG_LIM_YEARS ,"coverageyr","rejected"     ,"Has less than %d years of coverage",0},
  {12,REJECT_FLAG_MULTIPLE33,"multiple33","rejected"     ,"Has more than %d percent multiple exposure plates",0},
  {4,REJECT_FLAG_MIRA       ,"mira"      ,"rejected"     ,"Mira type variable",0},
  /* Eliminated from "new" list above this line:  REJECT_FLAG_NEW  newIndex */
  {10,REJECT_FLAG_ECLIPSING ,"eclipsing" ,"rejected"     ,"Eclipsing binary",0},
  {11,REJECT_FLAG_POLAR     ,"polar"     ,"rejected"     ,"Polar",0},
  {3,REJECT_FLAG_UG         ,"ugem"      ,"rejected"     ,"U Gem type variable",0},
  {5,REJECT_FLAG_RR         ,"rr"        ,"rejected"     ,"RR Lyr type variable",0},
  {7,REJECT_FLAG_VAR        ,"variable"  ,"rejected"     ,"Other type variable",0},
  {15,REJECT_FLAG_EXTRA     ,"extra"     ,"rejected"     ,"Has more than %d points outside the selected flare region",0},
  {26,REJECT_FLAG_YELLOW    ,"yellow"    ,"interest"     ,"Objects of secondary interest",0},
  {25,REJECT_FLAG_ORANGE    ,"orange"    ,"interest"     ,"Objects of primary interest",0},
  /* COMPLETELY REJECTED ABOVE THIS LINE:     REJECT_FLAG_SECONDBEST, secondBestIndex */
  {1,REJECT_FLAG_NEARBY     ,"nearby"    ,"nearby"       ,"Within %d magnitudes and %d arcsec of a nearby catalog star",0},  /* Fatal Josh's memo of Mon, 1 Jan 2018 13:43:46 14:34:20 22:55:07, later rescinded */
  {27,REJECT_FLAG_LOWRANGE  ,"lowrange"  ,"lowrange"     ,"Flare magnitude range less than %.2f mag",0},
  /* GOOD POINTS BELOW THIS LINE:             REJECT_FLAG_BEST, bestIndex */
  {6,REJECT_FLAG_NOVA       ,"nova"      ,"nova"         ,"Nova type variable",0},
  {14,REJECT_FLAG_UNRELEASED,"unreleased","nobadflag"    ,"At a galactic latitude below %d degrees",0},
  /* UNUSED selection points below this line: REJECT_FLAG_UNUSED, unusedIndex */
  {16,REJECT_FLAG_SDSS      ,"sdss"      ,""             ,"Matches SDSS catalog if with min(%.1f,%.1f*peakDradRMS3) arcsec and dimmer than magnitude %.1f",0},
  {18,REJECT_FLAG_MULTDEFECT,"multdefect",""             ,"Has more than %d percent multiple exposure + defect plates",0},
  {8,REJECT_FLAG_PREVIOUS   ,"previous"  ,""             ,"Previously studied",0},
  {29,REJECT_FLAG_LONGFLARE,"longflare"  ,"rejected" ,"Change of >  %.1f mag for > %d points over %.0f days after %.0f days average; then skip %.0f days",0},
  {21,REJECT_FLAG_LOWPOINTS ,"lowpoints" ,""             ,"Has less than %d points in the selected flare region",0},
  {28,REJECT_FLAG_MOREPOINTS,"morepoints",""             ,"Has more points within the flare than previous run",0},
};
int flareRejectTableSize = sizeof(flareRejectTable)/sizeof(CANDIDATETYPE);
/* In the following tables, the last index is applies to candidates without any conditions */
int rejectTotalCounts[REJECT_TABLE_SIZE+1];     /* total number of candidates for which this condition applies */
int rejectInitialCounts[REJECT_TABLE_SIZE+1];   /* number of candidates placed in the subdirectory for the condition */
char rejectionString[REJECT_TABLE_SIZE*20];
char secondBestString[REJECT_TABLE_SIZE*20];
char bestString[REJECT_TABLE_SIZE*20];
char newString[REJECT_TABLE_SIZE*20];
char unusedString[REJECT_TABLE_SIZE*20];

/* Exclusion zones come from /dasch/data/scanner/backup/2015_09_24/candidatedistribution.png */
extern EXCLUSION exclusionTable[];
extern int exclusionTableSize;
extern GSCBIN gscBin04;
PGSCBIN pGscBin = &gscBin04;


typedef struct _inputimage {
  double peakSlope;
  double peakUpperCount;
  double ra;
  double declination;
  double lat; /* Galactic latitude */
  double minDistance;  /* Distance in arcsec to closest candidate */
  double Stdmag;
  double peakRA; /* 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window */
  double peakDec; /* 1/sqr(dradRMS2) weighted positions of the good points within the selected transient window */
  double peakYear; /* year of the transient peak */
  double peakOffset; /* Distance between (ra,dec) distance for transient and (peakRA,peakDec) in arcsec */
  double peakNearbyDistance; /* the distance to nearby stars that are within 2 magnitudes of the brightest or average magnitude of the transient candidate */
  double rejectNearbyDistance; /* The peakNearbyDistance that sets the REJECT_FLAG_NEARBY */
  int candidate_index; /* Closest candidate */
  int input_index;  /* Our index */
  int nearbyREFflag;
  int npoints;
  int ngood;
  int peakCount;
  int peakCountNF;  /* Number of good lightcurve points from narrow field telescopes */
  int peakCountWF;  /* Number of good lightcurve points from wide field telescopes */
  int peakOutside;  /* Number of good lightcurve points outside the detection window */
  int peakOutsideNF;  /* Number of good lightcurve points outside the detection window from narrow field telescopes */
  int peakOutsideWF;  /* Number of good lightcurve points outside the detection window from wide field telescopes */
  int peakExtra;
  int peakDefectCount; /* the number of lightcurve defects in the selected transient window */
  int peakMultipleCount; /* the number of multiple exposure plates in the selected transient window */
  int gsc_bin_index;
  int peakBadColorCount;
  int peakLimitingYears;
  int peakLimitingPoints;
  int rejectFlags;
  int previousFlag;
  int newDirectoryFlag;
  int otherLongDirectoryFlag;
  int bestLongDirectoryFlag;
  int badLongDirectoryFlag;
  int bestDirectoryFlag;
  int secondBestDirectoryFlag;
  int rejectedDirectoryFlag;
  int rejectFlag; /* Used to suppress printing */
  int subdirectoryFlag; /* This one was copied to a subdirectory */
  int peakLongOutburst;
  double peakMaxMag;
  double peakMinMag;
  double peakNoDefectMag;
  double peakDradRMS2;
  double peakDradRMS3;
  double peakDays;
  double nearbyObjectDist;
  double Sig_TC; /* Max of peakDradRMS2 and peakDradRMS3 */
  char REF[MAX_REF];
  int nearbyCount;
  NEARBYIMAGE nearbyImage[MAX_NEARBY_COUNT];
  char nearbyObjects[2*MAX_NEARBY_OBJECTS_STRING+1];

} INPUT, *PINPUT;



double str2dec(		        /* Return Dec in degrees from string */
	const char* in);	/* Character string (dd:mm:ss.sss or dd.dddd) */
double str2ra(		        /* Return RA in degrees from string */
	const char* in);	/* Character string (hh:mm:ss.sss or dd.dddd) */
double wcsdist(	/* Compute angular distance between 2 sky positions */
               double ra1,	/* First longitude/right ascension in degrees */
               double dec1,	/* First latitude/declination in degrees */
               double ra2,	/* Second longitude/right ascension in degrees */
               double dec2);	/* Second latitude/declination in degrees */


/* When comparing GSC2.3.2 references, be sure we have the entire string */
int compareStrings(char *fullString,char* compareString) {
  int charCount = strlen(compareString);
  int charIndex;
  char fullValue;
  char partialValue;
  for (charIndex = 0; charIndex < charCount; charIndex++) {
    fullValue = fullString[charIndex];
    partialValue = compareString[charIndex];
    if (partialValue == ' ') {
      return(0);
    }
    if (fullValue != partialValue) {
      return(-1);
    }
  }
  return(0);
}
void PrintPossibleTable(FILE* possible_handle,FILE * all_handle,char *runDate,PTHUMBNAIL pThumbnail,PFLARECANDIDATE pCandidate,   PINPUT pInput, int gotGscBinIndexColumn, int verbose, int *pSubdirectoryCount) 
{
  char filename[MAX_INPUT_NAME+3];
  char *pChar2 = NULL;
  char flagsString[MAX_BUFFER];
  char noThumbnail[] = "NO THUMBNAIL FOUND";
  char noNearbyVSX[] = "";
  char nearbyVSX[2*MAX_NEARBY_OBJECTS_STRING+1];
  int rejectTableIndex;
  char countString[10];
  PCANDIDATETYPE pRejectTable;
  int nearbyIndex;
  PNEARBYIMAGE pNearbyImage;
  countString[0] = 0;
  flagsString[0] = 0;
  filename[0] = 0;

  if (pThumbnail == NULL) {
    strcpy(filename,noThumbnail);
  } else {
#if 0
    if (strstr(pThumbnail->filename,"APASS_J063648.7+405400") != NULL) {
      printf("line %d At %s\n",__LINE__,pThumbnail->filename);
    }
#endif
    if (strstr(pThumbnail->filename,"/n") == pThumbnail->filename) {
      strcpy(filename,pThumbnail->filename);
    } else {
      strcpy(filename,"/n");
      strcat(filename,pThumbnail->filename);
    }
  }
  pChar2 = strstr(pCandidate->comment,"; ;");
  while (pChar2 != NULL) {
    pChar2[0] = ' ';
    pChar2 = strstr(pCandidate->comment,"; ;");
  }

  pChar2 = pCandidate->comment;
  while (*pChar2 == ' ') {
    pChar2++;
  }
  if ((pInput != NULL) && (pInput->peakCount > 0)) {
    sprintf(countString,"#%d",pInput->peakCount);
  } else if (pCandidate->newPeakCount > 0) {
    sprintf(countString,"#%d",pCandidate->newPeakCount);
  } else if (pCandidate->prevPeakCount > 0) {
    sprintf(countString,"#%d",pCandidate->prevPeakCount);
  }
  
  if (possible_handle == NULL) {
    if (pCandidate->interestingFlag == 1) {
      return;
    }
    if (pCandidate->flarecandidate_index < 0) {
      return;
    }
    if ((pCandidate->refType == REF_TYPE_DASCH) || (pCandidate->refType == REF_TYPE_APASS)) {
      fprintf(all_handle,"%s #%s %s\n",
              pCandidate->name,
              countString,
              pChar2);
      pCandidate->allNameOutputFlag = 1;
    } else if (pCandidate->unrecognized == 0) {
      fprintf(all_handle,"%f %f #%s %s\n",
              pCandidate->ra,
              pCandidate->dec,
              countString,
              pChar2);
      pCandidate->allNameOutputFlag = 1;
    } else {
      fprintf(all_handle,"%s #%s %s\n",
              pCandidate->name,
              countString,
              pChar2);
      pCandidate->allNameOutputFlag = 1;
    }      
    return;
  }
  

#if 0
  if (strstr(pCandidate->comment,"Possible strong cand. Not on DSS so B >20") != NULL) {
    printf("line %d comment %s\n",__LINE__,pCandidate->comment);
  }
#endif


  if (verbose) {
    printf("Record found for candidate %s ra %f dec %f Thumbnail %s Comment %s\n",pCandidate->name,pCandidate->ra,pCandidate->dec,filename,pChar2);
  }
  if (pInput != NULL) {
    flagsString[0] = 0;
    if (pInput->rejectFlags == 0) {
      strcpy(flagsString,"noflags");
      if ((gotGscBinIndexColumn != 0) && (pInput->peakDradRMS3 < MAX_TC_PEAKDRADRMS3)) {
        strcat(flagsString,",lowdrad");
      } 
    } else if ((gotGscBinIndexColumn != 0) && (pInput->peakDradRMS3 < MAX_TC_PEAKDRADRMS3)) {
      strcat(flagsString,"lowdrad");
    }

    for (rejectTableIndex = 0; rejectTableIndex < flareRejectTableSize; rejectTableIndex++) {
      pRejectTable = &flareRejectTable[rejectTableIndex];
      if ((pRejectTable->flagmask & pInput->rejectFlags) != 0) {
        if (flagsString[0] == 0) {
          strcpy(flagsString,pRejectTable->flagname);
        } else {
          strcat(flagsString,",");
          strcat(flagsString,pRejectTable->flagname);
        }
        if (strlen(flagsString) > (MAX_BUFFER/2)) {
          printf("ERROR: flagsString length %d exceeds %d\n",strlen(flagsString),MAX_BUFFER/2);
          exit(-1);
        }
      }
    }
    if (pInput->rejectedDirectoryFlag != 0) {
      strcat(flagsString,",rejected");
    }
    if (pInput->secondBestDirectoryFlag != 0) {
      strcat(flagsString,",secondbest");
    }
    if (pInput->bestDirectoryFlag != 0) {
      strcat(flagsString,",best");
    }
    if (pInput->newDirectoryFlag != 0) {
      strcat(flagsString,",new");
    }
    if (pInput->otherLongDirectoryFlag != 0) {
      strcat(flagsString,",longother");
    }
    if (pInput->bestLongDirectoryFlag != 0) {
      strcat(flagsString,",longbest");
    }
    if (pInput->badLongDirectoryFlag != 0) {
      strcat(flagsString,",longbad");
    }
    nearbyVSX[0] = 0;
    if ((pInput->rejectFlags & REJECT_FLAG_SDSS) != 0) {
      for (nearbyIndex = 0; nearbyIndex < pInput->nearbyCount; nearbyIndex++) {
        pNearbyImage = &pInput->nearbyImage[nearbyIndex];
        if (pNearbyImage->sdssrejectflag) {
          if (nearbyVSX[0] == 0) {
            strcpy(nearbyVSX,pNearbyImage->nearbyObject);
          } else {
            strcat(nearbyVSX,"; ");
            strcat(nearbyVSX,pNearbyImage->nearbyObject);
          }
        }
      }
    }
    for (nearbyIndex = 0; nearbyIndex < pInput->nearbyCount; nearbyIndex++) {
      pNearbyImage = &pInput->nearbyImage[nearbyIndex];
      if (pNearbyImage->nearbyType == NEARBY_TYPE_VSX)  {
        if ((pNearbyImage->distance < MIN_VARIABLE_DISTANCE) ||
            ((pInput->Sig_TC < 90) && (pNearbyImage->distance < 3*pInput->Sig_TC))) {
          if (nearbyVSX[0] == 0) {
            strcpy(nearbyVSX,pNearbyImage->nearbyObject);
          } else {
            strcat(nearbyVSX,"; ");
            strcat(nearbyVSX,pNearbyImage->nearbyObject);
          }
        }
      }
    }
    if ((runDate[0] != 0) &&
        (pInput != NULL) &&
        (pCandidate != NULL) &&
        (all_handle != NULL)) {
      if ((strstr(pInput->REF,"DASCH_") == pInput->REF) ||
          (strstr(pInput->REF,"APASS_") == pInput->REF)) {
        fprintf(all_handle,"%s #%s %s\n",pInput->REF,countString,pChar2);
        pCandidate->allNameOutputFlag = 1;
      } else {
        if (strstr(pChar2,pInput->REF) != NULL) {
          if ((pInput->peakRA < 990.0) && (pInput->peakDec < 99.00)) {
            fprintf(all_handle,"%.4f %.4f #%s %s\n",
                    pInput->peakRA,
                    pInput->peakDec,
                    countString,
                    pChar2);
            pCandidate->allNameOutputFlag = 1;

          } else {
            fprintf(all_handle,"%.4f %.4f #%s %s\n",
                    pInput->ra,
                    pInput->declination,
                    countString,
                    pChar2);
            pCandidate->allNameOutputFlag = 1;
          }

        } else {
          if ((pInput->peakRA < 990.0) && (pInput->peakDec < 99.00)) {
            fprintf(all_handle,"%.4f %.4f #%s %s %s\n",
                    pInput->peakRA,
                    pInput->peakDec,
                    countString,
                    pInput->REF,
                    pChar2);
            pCandidate->allNameOutputFlag = 1;

          } else {
            fprintf(all_handle,"%.4f %.4f #%s %s %s\n",
                    pInput->ra,
                    pInput->declination,
                    countString,
                    pInput->REF,
                    pChar2);
            pCandidate->allNameOutputFlag = 1;
          }
        }
      }
    }

    if (pInput->subdirectoryFlag != 0) {
      (*pSubdirectoryCount)++;
      if ((pInput->peakRA < 990.0) && (pInput->peakDec < 99.00)) {
        fprintf(possible_handle,"%.4f\t%.4f\t%.1f\t%d\t%d\t%d\t%.4f\t%.1f\t%.1f\t%.1f\t%s\t'%s'\t'%s'\t'%s'\n",
                pInput->peakRA,
                pInput->peakDec,
                pInput->peakDradRMS3,
                pInput->peakLimitingYears,
                pInput->peakCount,
                pInput->peakExtra,
                pInput->peakYear,
                pInput->peakMaxMag,
                pInput->nearbyObjectDist,
                pInput->peakDays,
                filename,
                nearbyVSX,
                flagsString,
                pChar2);
      } else {
        fprintf(possible_handle,"%.4f?\t%.4f?\t%.1f\t%d\t%d\t%d\t%.4f\t%.1f\t%.1f\t%s\t'%s'\t'%s'\t'%s'\n",
                pInput->ra,
                pInput->declination,
                pInput->peakDradRMS3,
                pInput->peakLimitingYears,
                pInput->peakCount,
                pInput->peakExtra,
                pInput->peakYear,
                pInput->peakMaxMag,
                pInput->nearbyObjectDist,
                pInput->peakDays,
                filename,
                nearbyVSX,
                flagsString,
                pChar2);
      }
    }
  } else {
    pInput = NULL;
    nearbyVSX[0] = 0;
    fprintf(possible_handle,"%.4f?\t%.4f?\t%.1f\t%d\t%d\t%d\t%.4f\t%.1f\t%.1f\t%.1f\t%s\t'%s'\t'%s'\t'%s'\n",
            pCandidate->ra,
            pCandidate->dec,
            99.0,
            0,
            0,
            0,
            999.,
            99.0,
            999.0,
            999.0,
            filename,
            nearbyVSX,
            flagsString,
            pChar2);
  }
  return;
} /* End of PrintPossibleTable() */
int REFCompare(const void *first, const void *second) 
{
  char* REFFirst = (char *)((PINPUT)first)->REF;
  char* REFSecond = (char *)((PINPUT)second)->REF;
  return(strcmp(REFSecond,REFFirst));
}
int main(int argc,char *argv[]) {
  char *argstr;
  char cmdchar;
	int tabCount;
  char *inBuffer;
  char candidate_name[MAX_INPUT_NAME];
  char interesting_name[MAX_INPUT_NAME];
  int iteration; /* 0 for candidate_name and 1 for interesting_name */
  char out_name[MAX_INPUT_NAME];
  char all_name[MAX_INPUT_NAME];
  char thumbnail_name[MAX_INPUT_NAME];
  char possible_name[MAX_INPUT_NAME];
  char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
	int lineCounter = 0;
  int goodLines = 0;
	int nvals;
	int errorFlag = 0;
  FILE * candidate_handle = NULL;
  FILE * interesting_handle = NULL;
  FILE * inHandle = NULL;
  FILE * out_handle = NULL;
  FILE * all_handle = NULL;
  FILE * thumbnail_handle = NULL;
  FILE * possible_handle = NULL;
	int lineLen;
	int tabIndex;
  int verbose = 0;
  int printLocationHistogram = 0;
  int locationHistogram[(MAX_LOCATION_HISTOGRAM_ARCSEC/LOCATION_HISTOGRAM_BINSIZE)];
  int bestIndex = 0;
  int newIndex = 0;
  int secondBestIndex = 0;
  int unusedIndex = 0;
  int promoteFlag = 0;
  int printMoveCommands = 0;
  int printWarnings = 0;
  int enforcePeakExcess = 0;
  int printSummaryTable = 0;
  int showsametimematches = 0;
  int printLocationTable = 0;
  int printMultiplePercentTable = 0;
  int nearbyREFflag;
  char *sdssPtr;
  double sdssGMag;
  double sdssColor;
  double sdssDist;
  double peakSlope;
  int peakCount;
  int peakUpperCount;
  int acceptFlag;
  int lval;
  int candidate_count = 0;
  int interesting_count = 0;
  int candidate_index = 0;
  int interesting_index = 0;
  int candidate2_index = 0;
  PFLARECANDIDATE candidate_table;
  PFLARECANDIDATE pCandidate;
  PFLARECANDIDATE pCandidate2;

  PFLARECANDIDATE unrecognized_table;
  PFLARECANDIDATE pUnrecognized;
  int unrecognized_count = 0;
  int unrecognized_index;
  FLARECANDIDATE dummy_candidate;


  PTHUMBNAIL thumbnail_table = NULL;
  PTHUMBNAIL pThumbnail;
  PTHUMBNAIL pBestThumbnail;
  char noThumbnail[] = "NO THUMBNAIL FOUND";
  int thumbnail_count = 0;
  int thumbnail_index = 0;
  int daynumber;
  int monthnumber;
  int yearnumber;

  File id_handle = NULL;
  char rootDirectory[MAX_BUFFER];
  char matchedTimestamp[MAX_BUFFER];
  char runDate[MAX_BUFFER];
  char id_name[MAX_BUFFER];



  TableHead input_header = NULL;
  PINPUT input_table = NULL;
  size_t input_nrecs = 0;
  int input_index;
  int input_index2;
  int input_index3;
  PINPUT pInput;
  PINPUT pInput2;
  PINPUT pInput3;
  PINPUT pInput4;
  double bestdistance;
  int intdistance;
  int cumdistance;
  char *pChar;
  char *pChar2;
  int gsc_bin_index;
  int *gsc_bin_count_single = NULL;
  int *gsc_bin_count_multiple = NULL;
  int gsc_single_candidate_count = 0;
  int gsc_bin_count_total = 0;
  int gsc_bin_count_maximum = 0;
  double tmpra;
  double tmpdec;
  int decBin;
  int raBin;
  int bin_count;
  int bin_index;
  int bin_list[MAX_ADJACENT_BINS+1];
  int ratio_count = 0;
  int sametimeCount = 0;
 

  char * pComment;
  double distance;
  char charVal;
  char *pSlash;
  char *pUnderscore;
  char *pSign;
  char runTimestamp[MAX_SUBDIR_NAME+1];
  int columnIndex;
  int gotPeakExtraColumn = 0;
  int gotPeakExcessColumn = 0;
  int gotPeakNearbyDistanceColumn = 0;
  int gotPeakMultipleCountColumn = 0;
  int gotGscBinIndexColumn = 0;
  int gotPeakBadColorCountColumn = 0;
  int gotPeakLimitingYearsColumn = 0;
  int gotPeakDaysColumn = 0;
  int gotPeakNoDefectMagColumn = 0;
  int gotPeakLongOutburstColumn = 0;
  char nearbyObjects[2*MAX_NEARBY_OBJECTS_STRING+1];
  PNEARBYIMAGE pNearbyImage;
  int nearbyIndex;
  int nearbyLength;
  int refType;
  long long REFNumber;
  char cmdLine[MAX_BUFFER];
  int rejectTableIndex;
  PCANDIDATETYPE pRejectTable;
  int exclusionIndex;
  PEXCLUSION pExclusion;
  int possible;
  int possibleCount = 0;
  int inputDupCount = 0;
  int possibleMissingCount = 0;
  time_t startTime;
  time_t curTime;
  double factor;
  double sdssDistance;
  void *candidate_ptree; /* kdtree for the previous candidates */
  void *input_ptree; /* kdtree for the input file */
  double x;
  double y;
  double z;
	double cosx;
	double pt[3];  
	double pos[3];
	double deltaDec;
	double deltaRA;
	struct kdres *presults;
	int resultSize;
  int foundMatchFlag = 0;
  int thumbnailCandidateNearbyCount = 0;
  int thumbnailCandidateBetterCount = 0;
  int thumbnailNoCoordsCount = 0;

  int promoteBothCount = 0;
  int promoteLowdradCount = 0;
  int promoteSDSSCount = 0;
  int bestCount = 0;
  int secondBestCount = 0;
  int morepointsCount = 0;
  int interestCount = 0;
  int newCount = 0;
  int otherLongFlareCount = 0;
  int bestLongFlareCount = 0;
  int badLongFlareCount = 0;
  int novaCount = 0;
  int ugemCount = 0;
  int variableCount = 0;
  int subdirectoryCount = 0;
  int workDone;
  char* tmpChar;
  int preallocFlag = 0;


  time(&startTime);
  assert(flareRejectTableSize <= REJECT_TABLE_SIZE);
  memset(rejectTotalCounts,0,sizeof(rejectTotalCounts));
  memset(rejectInitialCounts,0,sizeof(rejectInitialCounts));
  memset(locationHistogram,0,sizeof(locationHistogram));
 

  candidate_ptree = kd_create(3); /* Create kdtree for the candidates */
  input_ptree = kd_create(3);   /* Create kdtree for the input file */

#ifdef SKIP_PREVIOUS
  printf("ERROR: SKIP_PREVIOUS is defined\n");
#endif /* SKIP_PREVIOUS */
#ifdef PREFILTER_MULTDEFECT
  printf("ERROR: PREFILTER_MULTDEFECT is defined\n");
#endif /* PREFILTER_MULTDEFECT */


  for (exclusionIndex = 0; exclusionIndex < exclusionTableSize; exclusionIndex++) {
    pExclusion = &exclusionTable[exclusionIndex];
    if (pExclusion->ra2 < pExclusion->ra1) {
      printf("ERROR: exclusion table right ascensions are out of order for index %d\n",
             pExclusion->ra1,
             pExclusion->ra2,
             exclusionIndex);
      exit(-1);
    }
    if (pExclusion->dec2 < pExclusion->dec1) {
      printf("ERROR: exclusion table right ascensions are out of order for index %d\n",
             pExclusion->dec1,
             pExclusion->dec2,
             exclusionIndex);
      exit(-1);
    }
       
  }

	/* Loop through the arguments */
	id_name[0] = 0;
	out_name[0] = 0;
	all_name[0] = 0;
  candidate_name[0] = 0;
  interesting_name[0] = 0;
  thumbnail_name[0] = 0;
  possible_name[0] = 0;

  rootDirectory[0] = 0;
  matchedTimestamp[0] = 0;
  runDate[0] = 0;

	for (argv++; --argc > 0; argv++) {
		argstr = *argv;
		/* Decode arguments */
		if (argstr[0] != '-') {
			errorFlag = 1;
			printf("ERROR: stray argument %s\n",argstr);
		} else {
			while ((cmdchar = *++argstr) != 0) {
				switch(cmdchar) {
				case 'v':
				case 'V':
					verbose = 1;
					break;


				case 'm':
				case 'M':
					printMoveCommands = 1;
					break;

        case 'l':
        case 'L':
          printLocationHistogram = 1;
          break;

        case 'g':
        case 'G':
          printLocationTable = 1;
          break;


				case 'w':
				case 'W':
					printWarnings = 1;
					break;

				case 'e':
				case 'E':
					enforcePeakExcess = 1;
					break;

				case 's':
				case 'S':
					printSummaryTable = 1;
					break;

				case 't':
				case 'T':
					showsametimematches = 1;
					break;

				case 'q':
				case 'Q':
					printMultiplePercentTable = 1;
					break;

				case 'i': /* matched id file name */
				case 'I':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(id_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;


				case 'c': /* candidate file name */
				case 'C':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(candidate_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;

				case 'j': /* Josh's interesting object list */
				case 'J':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(interesting_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;


				case 'f': /* candidate thumbnail file list */
				case 'F':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(thumbnail_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;

				case 'p': /* best selection list */
				case 'P':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(possible_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
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
						strncpy(out_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;

				case 'a': /* all candidates files for the next input */
				case 'A':
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(all_name,*++argv,MAX_INPUT_NAME-2);
						if (strlen(*argv) >= MAX_INPUT_NAME-2) {
							fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
						}
					}
					break;

        case 'd':  /* Root directory and timestamp */
        case 'D':
          if (argc < 3) {
            printf("Error: Insufficient arguments for -%c",cmdchar);
            errorFlag = 1;
          } else {
            argc--;
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: rootDirectory length %d exceeds %d\n",strlen(*argv),MAX_BUFFER-1);
              errorFlag = 1;
            } else {
              strcpy(rootDirectory,*argv);
              argc--;
              ++argv;
              if (strlen(*argv) >= (MAX_BUFFER-1)) {
                printf("ERROR: runDate length %d exceeds %d\n",strlen(*argv),MAX_BUFFER-1);
                errorFlag = 1;
              } else {
                strcpy(runDate,*argv);
                strcpy(matchedTimestamp,runDate);
              }
            }
          }
          break;


				default:
					printf("* illegal command -%c-",cmdchar);
					errorFlag = 1;
					break;
				}
        
			}
		}
	}



  id_handle = Open(id_name,"r");
  if (id_handle == NULL) {
    printf("ERROR: No input file found  %s\n",id_name);
    errorFlag = 1;
  } else {
    if (verbose) {
      printf("Found plate input file %s\n",id_name);
    }
  }


	if ((thumbnail_name[0] != 0) && (candidate_name[0] == 0)) {
    printf("ERROR: thumbnail file specified without a candidate file\n");
    errorFlag = 1;
  } else {
    if (thumbnail_name[0] != 0) {
      thumbnail_handle = fopen(thumbnail_name,"rt");
      if (thumbnail_handle == NULL) {
        printf("ERROR: No thumbnail file found  %s\n",thumbnail_name);
        errorFlag = 1;
      } else {
        if (verbose) {
          printf("Found plate thumbnail file %s\n",thumbnail_name);
        }
      }
    }
  }

	if ((possible_name[0] != 0) && (candidate_name[0] == 0)) {
    printf("ERROR: possible file specified without a candidate file\n");
    errorFlag = 1;
  } else {
    if (possible_name[0] != 0) {
      possible_handle = fopen(possible_name,"wt");
      if (possible_handle == NULL) {
        printf("ERROR: Failed to open the possible file  %s\n",possible_name);
        errorFlag = 1;
      } else {
        if (verbose) {
          printf("Opened possible file %s\n",possible_name);
        }
      }
    }
  }


	if (out_name[0] == 0) {
		printf("ERROR: Out file name not specified\n");
		errorFlag = 1;
	} else {
		out_handle = fopen(out_name,"wt");
		if (out_handle == NULL) {
			errorFlag = 1;
			printf("Could not open out file %s\n",out_name);
		}
	}

	if (all_name[0] != 0) {
		all_handle = fopen(all_name,"wt");
		if (all_handle == NULL) {
			errorFlag = 1;
			printf("Could not open the next previous candidates file file %s\n",all_name);
		}
	}



	if (candidate_name[0] != 0) {
		candidate_handle = fopen(candidate_name,"rt");
		if (candidate_handle == NULL) {
			errorFlag = 1;
			printf("Could not open list file %s\n",candidate_name);
		}
	} else {
    errorFlag = 1;
    printf("ERROR: No flare candidate filename found\n");
  }

	if (interesting_name[0] != 0) {
		interesting_handle = fopen(interesting_name,"rt");
		if (interesting_handle == NULL) {
			errorFlag = 1;
			printf("Could not open list file %s\n",interesting_name);
		}
	} else {
    printf("WARNING: No objects of interest filename found\n");
  }

  for (rejectTableIndex = 0; rejectTableIndex < flareRejectTableSize; rejectTableIndex++) {
    pRejectTable = &flareRejectTable[rejectTableIndex];
    switch (pRejectTable->flagmask) {
    case REJECT_FLAG_SDSS:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MAX_SDSS_DISTANCE,MAX_SDSS_DISTANCE_FACTOR,MAX_SDSS_MAGNITUDE);
      break;
    case REJECT_FLAG_LOWPOINTS:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,LOW_TC_POINTS);
      break;

    case REJECT_FLAG_EXTRA:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MAX_PEAK_EXTRA_COUNT);
      break;

    case REJECT_FLAG_MULTIPLE33:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MAX_MULTIPLE_PERCENTAGE_33-1);
      break;
    case REJECT_FLAG_LONGFLARE:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,LONG_TC_MAGNITUDE,LONG_TC_POINTS,LONG_TC_FLARE_DAYS,LONG_TC_PRE_DAYS,LONG_TC_SKIP_DAYS);
      break;

    case REJECT_FLAG_MULTIPLE90:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MAX_MULTIPLE_PERCENTAGE_90-1);
      break;

    case REJECT_FLAG_MULTDEFECT:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MAX_DEFECT_PERCENTAGE-1);
      break;

    case REJECT_FLAG_DEFECT:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MAX_DEFECT_PERCENTAGE-1);
      break;

    case REJECT_FLAG_LIM_YEARS:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MIN_TC_LIMITING_YEARS);
      break;

    case REJECT_FLAG_LOWRANGE:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,LARGE_TC_MAGNITUDE);
      break;

    case REJECT_FLAG_LIM_POINTS:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MIN_TC_LIMITING_POINTS);
      break;

    case REJECT_FLAG_UNRELEASED:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MAX_GALACTIC_LATITUDE);
      break;

    case REJECT_FLAG_NEARBY:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,MAX_TC_NEARBY_MAG,MAX_TC_NEARBY_ARCSEC);
      break;

    case REJECT_FLAG_SAMETIME:
      sprintf(pRejectTable->webname,pRejectTable->webnametemplate,365.25*24.0*MIN_TC_DATERANGE);
      break;



    default: 
      strcpy(pRejectTable->webname,pRejectTable->webnametemplate);
      break;
    }
  }

  for (rejectTableIndex = 0; rejectTableIndex < flareRejectTableSize; rejectTableIndex++) {
    pRejectTable = &flareRejectTable[rejectTableIndex];
#if 0
    if (rejectTableIndex == 23) {
      printf("At rejectTableIndex %d\n",rejectTableIndex);
    }
#endif
    if (pRejectTable->flagmask == REJECT_FLAG_BEST) {
      bestIndex = rejectTableIndex;
    }
    if (pRejectTable->flagmask == REJECT_FLAG_SECONDBEST) {
      secondBestIndex = rejectTableIndex;
    }
    if (pRejectTable->flagmask == REJECT_FLAG_UNUSED) {
      unusedIndex = rejectTableIndex;
    }
    if (pRejectTable->flagmask == REJECT_FLAG_NEW) {
      newIndex = rejectTableIndex;
    }
  }

	if (errorFlag) {
		printf("Usage: parseidtable\n");
    printf(" (RT)  where -i <id_name> is the candidates table from search_none\n");
    printf(" (WT)  where -o <out_name> is the new candidates table (saved only when looking for new detections)\n");
    printf("       where -e means enforce peakExtra <= %d criterion\n",MAX_TC_EXCESS_POINTS);
    printf("       where -d <directory> <unmatched> <matched> <rundate> is the root directory, unmatched timestamp, matched timestamp, and run date\n");
    printf(" (RT)  where -f <thumbnail_name> is the list of candidate thumbnail files produced by plottransient.c (input only)\n");
    printf(" (RT)  where -c <candidate_name> is a list of previously reported candidates (input only)\n");
    printf(" (RT)  where -j <interesting_name> is the list of Josh's candidates of interest (input only)\n");
    printf(" (WT)  where -p <possible_name> is the list of best choice flares\n");
    printf(" (WT)  where -a <all_name> the generated new previous candidates file\n");
		printf("       where -v is the verbose flag\n");
    printf("       where -m prints mkdir, ls, and mv commands to move the files\n");
    printf("       where -w prints warnings for missing possible flare candidates\n");
    printf("       where -s means print summary table\n");
    printf("       where -q means print percentages of multiple exposure plates\n");
    printf("       where -l means to print the location histogram to %d arcsec\n");
    printf("       where -t means to print matches with %.2f hour or less difference in peak times\n",365.25*24.0*MIN_TC_DATERANGE);    

		return(-1);
	}

	printf("parseidtable of %s %s ID Filename %s\n Output Filename %s\n Thumbnail name %s\n Previous candidates %s\n Objects of interest %s\n Best choice name %s\n",
				 __DATE__,__TIME__,id_name,out_name,thumbnail_name,candidate_name,interesting_name,possible_name);

  printf("best if the '%s' condition or better in subdirectory '%s' (%s)\n",flareRejectTable[bestIndex].flagname,flareRejectTable[bestIndex].subdir,flareRejectTable[bestIndex].webname);
  printf("secondbest if the '%s' condition or better in subdirectory '%s' (%s)\n",flareRejectTable[secondBestIndex].flagname,flareRejectTable[secondBestIndex].subdir,flareRejectTable[secondBestIndex].webname);
  printf("unused if the '%s' condition or better in subdirectory '%s' (%s)\n",flareRejectTable[unusedIndex].flagname,flareRejectTable[unusedIndex].subdir,flareRejectTable[unusedIndex].webname);
  printf("new if the '%s' condition or better in subdirectory '%s' (%s)\n",flareRejectTable[newIndex].flagname,flareRejectTable[newIndex].subdir,flareRejectTable[newIndex].webname);
  if (rootDirectory[0] != 0) {
    for (rejectTableIndex = 0; rejectTableIndex < flareRejectTableSize; rejectTableIndex++) {
      pRejectTable = &flareRejectTable[rejectTableIndex];
      printf("index %2d flagbit %2d flagmask %10d Flag '%s': '%s' '%s'\n",
             rejectTableIndex,
             pRejectTable->flagbit,
             pRejectTable->flagmask,
             pRejectTable->flagname,
             pRejectTable->subdir,
             pRejectTable->webname);
      
      if ((1 << pRejectTable->flagbit) != pRejectTable->flagmask) {
        printf("ERROR: bit %d does not agree with mask 0x%x for %s %s\n",
               pRejectTable->flagbit,
               pRejectTable->flagmask,
               pRejectTable->flagname,
               pRejectTable->webname);
        exit(-1);
      } 
#if 0
      if (printMoveCommands) {
        if (rejectTableIndex > bestIndex) {
          printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,pRejectTable->subdir);
#if 0
          printf("ls %s/%s/%s/*.png | wc -l\n",rootDirectory,runDate,pRejectTable->subdir);
#endif
        }
      }
#endif
    }
    if (printMoveCommands) {
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"best");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"secondbest");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"interest");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"new");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"nova");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"ugem");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"variable");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"morepoints");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"bestlongflare");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"otherlongflare");
      printf("mkdir -p %s/%s/%s\n",rootDirectory,runDate,"badlongflare");
#if 0
      printf("ls %s/%s/%s/*.png | wc -l\n",rootDirectory,runDate,"best");
      printf("ls %s/%s/%s/*.png | wc -l\n",rootDirectory,runDate,"secondbest");
      printf("ls %s/%s/%s/*.png | wc -l\n",rootDirectory,runDate,"interest");
      printf("ls %s/%s/%s/*.png | wc -l\n",rootDirectory,runDate,"new");
#endif
    }
  }

  if (candidate_name[0] != 0) {
    /* First find out how many candidates we have */
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,candidate_handle);
      if (inBuffer == NULL) {
        break;
      }
      lineLen = strlen(inBuffer);
      /* Trim off the carriage return */
      if (inBuffer[lineLen-1] == 10) {
        inBuffer[lineLen-1] = 0;
        lineLen--;
      }
      if (strlen(inBuffer) >= (MAX_BUFFER-2)) {
        printf("ERROR parseid line %4d MAX_BUFFER exceeded\n",__LINE__);
      }
      strcpy(copyLine,inBuffer);
     
      lineCounter++;
    }
    fclose(candidate_handle);
    candidate_handle = NULL;
    candidate_count = lineCounter;
    candidate_index = 0;
  }
  if (interesting_name[0] != 0) {
    /* First find out how many objects of interest we have */
    lineCounter = 0;
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,interesting_handle);
      if (inBuffer == NULL) {
        break;
      }
      lineLen = strlen(inBuffer);
      /* Trim off the carriage return */
      if (inBuffer[lineLen-1] == 10) {
        inBuffer[lineLen-1] = 0;
        lineLen--;
      }
      if (strlen(inBuffer) >= (MAX_BUFFER-2)) {
        printf("ERROR parseid line %4d MAX_BUFFER exceeded\n",__LINE__);
      }
      strcpy(copyLine,inBuffer);
     
      lineCounter++;
    }
    fclose(interesting_handle);
    interesting_handle = NULL;
    interesting_count = lineCounter;
    interesting_index = 0;
  }
  if ((candidate_count == 0) &&
      (interesting_count == 0)) {
    printf("ERROR both candidate_count %d and intereting_count %d are zero\n");
    exit(-1);
  }

  candidate_table = (PFLARECANDIDATE)calloc(candidate_count+interesting_count+1,sizeof(FLARECANDIDATE));
  if (candidate_table == NULL) {
    printf("ERROR: failed to allocate the candidate table\n");
    exit(-1);
  }
  
  unrecognized_table = (PFLARECANDIDATE)calloc(candidate_count+interesting_count+1,sizeof(FLARECANDIDATE));
  if (unrecognized_table == NULL) {
    printf("ERROR: failed to allocate the unrecognized table\n");
    exit(-1);
  }


  preallocFlag = 1;
  for (candidate_index = 0; candidate_index < (candidate_count+interesting_count+1); candidate_index++)  {
    pCandidate = &candidate_table[candidate_index];
    pCandidate->gsc_bin_index = -1;
    pCandidate->ra = 999.;
    pCandidate->dec = 99.;
  }

  if (InitFlareCandidateList(pGscBin,interesting_name,&interesting_count,&candidate_table,candidate_ptree,&startTime,preallocFlag) != 0) {
    exit(-1);
  }

#if 1
 for (candidate_index = 0; candidate_index < interesting_count; candidate_index++)  {
    pCandidate = &candidate_table[candidate_index];
    printf("Line %d candidate_index %2d DEBUG ra %8.4f dec %8.4f gsc_bin_index %10d name %s REFNumber %lld\n",
           __LINE__,
           candidate_index,
           pCandidate->ra,
           pCandidate->dec,
           pCandidate->gsc_bin_index,
           pCandidate->name,
           pCandidate->REFNumber);
 }
#endif
  for (iteration = 0; iteration < 2; iteration++) {
    if (iteration == 0) {
      if (candidate_count == 0) {
        continue;
      }
      candidate_handle = fopen(candidate_name,"rt");
      if (candidate_handle == NULL) {
        errorFlag = 1;
        printf("Could not open list file %s\n",candidate_name);
        exit(-1);
      }
    } else {
      if (interesting_count == 0) {
        continue;
      }
      candidate_handle = fopen(candidate_name,"rt");
      if (candidate_handle == NULL) {
        errorFlag = 1;
        printf("Could not open list file %s\n",candidate_name);
        exit(-1);
      }
      interesting_handle = NULL;
    }



    lineCounter = 0;
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,candidate_handle);
      if (inBuffer == NULL) {
        break;
      }
      lineLen = strlen(inBuffer);
      /* Trim off the carriage return */
      if (inBuffer[lineLen-1] == 10) {
        inBuffer[lineLen-1] = 0;
        lineLen--;
      }
      /* Trim trailing spaces */
      while ((lineLen >= 0) && (inBuffer[lineLen-1] == ' ')) {
        inBuffer[lineLen-1] = 0;
        lineLen--;
      }
        
      lineCounter++;
#ifdef LOS_DEBUG1
      if (strstr(inBuffer,LOS_REF) != NULL) {
        printf("Line %d at REF %s\n",__LINE__,inBuffer);
      }
#endif
      strcpy(copyLine,inBuffer);
      local_strlwr(copyLine);
      /* Parse the line looking for strings of type yyyy_mm_dd */
      runTimestamp[0] = 0;
      if (iteration == 0) {
        while (1) {
          pChar = strrchr(copyLine,' ');
          if (pChar == NULL) {
            break;
          }
          *pChar = 0;
          pChar++;
          if (strlen(pChar) == MAX_SUBDIR_NAME) {
            nvals = sscanf(pChar,"%d_%d_%d",&yearnumber,&monthnumber,&daynumber);
            if (nvals == 3) {
              if (runTimestamp[0] == 0) {
                strcpy(runTimestamp,pChar);
              } else {
                if (strcmp(pChar,runTimestamp) > 0) {
                  strcpy(runTimestamp,pChar);
                }
              }
              break;
            }
          }
        }
      }
      strcpy(copyLine,inBuffer);
#if 0
      if (strstr(copyLine,"ngc5204") != NULL) {
        printf("At line %d\n",lineCounter);
      }
#endif

      pCandidate = &candidate_table[candidate_index];
#if 1
      printf("DEBUG: parseidtable line %4d candidate_index %5d \n",__LINE__,candidate_index);
#endif
      memset(pCandidate,0,sizeof(FLARECANDIDATE));
      pCandidate->flarecandidate_index = -1;
      pCandidate->input_index = -1;
      pCandidate->interestingFlag = iteration;
      pCandidate->ra = 999.;
      pCandidate->dec = 99.;
      pCandidate->gsc_bin_index = -1;
      if (iteration == 0) {
        pCandidate->previousFlag = 1;
      }
      strcpy(pCandidate->runTimestamp,runTimestamp);
      pComment = strstr(inBuffer,"#");
      if (pComment != NULL) {
#ifdef LOS_DEBUG2
        tmpChar = strstr(inBuffer,"#X");
        if (tmpChar != NULL) {
          printf("DEBUG: line %4d candidate_index %d copyLine '%s'\n",__LINE__,candidate_index,copyLine);
        }
#endif
        pCandidate->typeFlag = pComment[1];
        *pComment = 0;
        pComment++;
        if (pCandidate->typeFlag == '#') {
          /* We have a peakCount following */
          *pComment = 0;
          pComment++;
          nvals = sscanf(pComment,"%d",&pCandidate->prevPeakCount);
          if (nvals == 1) {
            while ((*pComment >= '0') && (*pComment <= '9')) {
              *pComment = 0;
              pComment++;
            }
          } else {
            printf("ERROR: Can not decode prevPeakCount from comment %s in line %d iteration %d\n",pComment,lineCounter,iteration);
            pCandidate->prevPeakCount = 0;
          }

        }

        if (strlen(pComment) >= MAX_COMMENT_SIZE) {
          printf("ERROR: comment %s, %d is too long in line %d iteration %d\n",pComment,strlen(pComment),lineCounter,iteration);
        } else {
          strcpy(pCandidate->comment,pComment);
#if 0
          if (strstr(pCandidate->comment,"Possible strong cand. Not on DSS so B >20") != NULL) {
            printf("line %d candidate_index %d comment %s\n",__LINE__,candidate_index,pCandidate->comment);
          }
#endif
        }
      }
      if (strlen(inBuffer) >= MAX_FLARECANDIDATE_SIZE) {
        printf("ERROR: candidate name %s, %d is too long in line %d iteration %d\n",inBuffer,strlen(inBuffer),lineCounter,iteration);
        strcpy(pCandidate->name,"ERROR: BAD SIZE");
      } else {
        strcpy(pCandidate->name,inBuffer);
      }
      /* Remove trailing spaces */
      lval = strlen(pCandidate->name);
      while (pCandidate->name[lval-1] == ' ') {
        pCandidate->name[lval-1] = 0;
        lval--;
      }

      if (strlen(pCandidate->name) == 0) {
#if 1
        printf("DEBUG: line %4d missing name copyLine '%s' \n",__LINE__,copyLine);
#endif
        continue;
      }
#if 0
      if (strstr(pCandidate->name,"DASCH") != NULL) {
        printf("DASCH ref type %s\n",pCandidate->name);
      } else if (strstr(pCandidate->name,"APASS") != NULL) {
        printf("APASS ref type %s\n",pCandidate->name);
      } else if (strstr(pCandidate->name,"257.060833") != NULL) {
        printf("Other %s\n",pCandidate->name);
      }
#endif    
        


      if (parse_coo(pCandidate->name,&pCandidate->ra,&pCandidate->dec) != 0) {
        if (GetREFNumber(pCandidate->name,&pCandidate->REFNumber,&pCandidate->refType,0,0) == 0) {
          switch(pCandidate->refType) {
          case REF_TYPE_NONE:
            pCandidate->unrecognized = 1;
            break;
          case REF_TYPE_GSC:
          case REF_TYPE_TYCHO2:
          case REF_TYPE_GAIA1:
          case REF_TYPE_GAIA2:
          case REF_TYPE_ATLAS2:
            pCandidate->unrecognized = 1;
            break;
          case REF_TYPE_KEPLER:
            pCandidate->unrecognized = 1;
            break;
          case REF_TYPE_DASCH:
            if (GetDASCHCoordinates(pCandidate->name,
                                    &pCandidate->ra,
                                    &pCandidate->dec,0,REF_TYPE_DASCH) != 0) {
              pCandidate->unrecognized = 1;
            }

            break;
          case REF_TYPE_APASS:
            if (GetDASCHCoordinates(pCandidate->name,
                                    &pCandidate->ra,
                                    &pCandidate->dec,0,REF_TYPE_APASS) != 0) {
              pCandidate->unrecognized = 1;
            }

            break;
          default:
            break;
          }

        } else {
          pCandidate->unrecognized = 1;
          pCandidate->REFNumber = 0;
        }
      } else {
#ifdef LOS_DEBUG2
        printf("DEBUG line %4d parse_coo succeeded\n",__LINE__);
#endif
      }
      if (pCandidate->unrecognized == 1) {
        printf("ERROR: line %4d lineCounter %3d unrecognized candidate: %s REFNumber %lld lineCounter %d iteration %d\n",__LINE__,lineCounter,pCandidate->name,pCandidate->REFNumber,lineCounter,iteration);
        if (pCandidate->REFNumber != 0) {
          printf("SELECT REFNumber,ra,declination from stars2 where REFNumber = %lld;\n",pCandidate->REFNumber);
        }
        /* Save this one in case we get the RA and DEC from the id table */
        pUnrecognized = &unrecognized_table[unrecognized_count];
        memcpy(pUnrecognized,pCandidate,sizeof(FLARECANDIDATE));
        unrecognized_count++;
      } else {
#if 0
        printf("found index %d Name: %s comment: %s ra: %f dec %f\n",
               candidate_index,
               pCandidate->name,
               pCandidate->comment,
               pCandidate->ra,
               pCandidate->dec);
#endif
        /* Query the kd tree for a nearby candidate */
        z = sin(DEGREES_TO_RAD*pCandidate->dec);
        cosx = cos(DEGREES_TO_RAD*pCandidate->dec);
        x = cosx*cos(DEGREES_TO_RAD*pCandidate->ra);
        y = cosx*sin(DEGREES_TO_RAD*pCandidate->ra);
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
        foundMatchFlag = 0;
        presults = kd_nearest_range(candidate_ptree, pt,(MAX_FILE_DISTANCE*DEGREES_TO_RAD)/3600.0);
        resultSize = kd_res_size(presults);
        while  (!kd_res_end( presults ) ) {
          pCandidate2 = (PFLARECANDIDATE)kd_res_item(presults,pos);
          factor = cos(DEGREES_TO_RAD* pCandidate2->dec);
          deltaDec = pCandidate2->dec - pCandidate->dec;
          if (deltaDec < 0) {
            deltaDec = - deltaDec;
          }
          deltaRA  = factor *(pCandidate2->ra  - pCandidate->ra);
          if (deltaRA < 0) {
            deltaRA = -deltaRA;
          }
          distance = sqrt(sqr(deltaDec)+sqr(deltaRA)) * 3600.;
          
          if (foundMatchFlag != 0) {
            printf("ERRORD: multiple matches %s to the new candidate %s in line %d lineCounter %d iteration %d\n",pCandidate2->name,pCandidate->name,__LINE__,lineCounter,iteration);
          } else {
            foundMatchFlag = 1;
#ifdef LOS_DEBUG1
            if (pCandidate2->REFNumber == LOS_REF_NUMBER) {
              printf("Line %d at REFNumber %lld\n",__LINE__,pCandidate->REFNumber);
            }
#endif
            printf("ERRORD: line %d candidates %s and %s are within %.1f arcsec lineCounter %d iteration %d\n",
                   __LINE__,pCandidate->name,pCandidate2->name,distance,lineCounter,iteration);
            if ((strlen(pCandidate->comment) + strlen(pCandidate2->comment) +2) > MAX_COMMENT_SIZE) {
              printf("ERRORD MAX_COMMENT_SIZE exceeded by %d in lineCounter %d iteration %d\n",(strlen(pCandidate->comment) + strlen(pCandidate2->comment) +1),lineCounter,iteration);
            } else {
              if (iteration == 0) {
                strcat(pCandidate2->comment,"; ");
                strcat(pCandidate2->comment,pCandidate->comment);
              } else {
                if (pCandidate->prevPeakCount > pCandidate2->prevPeakCount) {
                  pCandidate2->prevPeakCount = pCandidate->prevPeakCount; 
                }
                /* Throw away the old comment completely and set the interesting flag */
                strcpy(pCandidate2->comment,pCandidate->comment);
                pCandidate2->REFNumber = pCandidate->REFNumber;
                pCandidate2->ra = pCandidate->ra;
                pCandidate2->dec = pCandidate->dec;
                pCandidate2->refType = pCandidate->refType;
                pCandidate2->typeFlag = pCandidate->typeFlag;
                strcpy(pCandidate2->name,pCandidate->name);
                pCandidate2->interestingFlag = 1;
              }
            }
          }
          /* go to the next entry */
          kd_res_next( presults );
            
        }
        
        kd_res_free(presults);
   
        if (foundMatchFlag == 0) {
          /* Insert this one in the kdtree */
          z = sin(DEGREES_TO_RAD*pCandidate->dec);
          cosx = cos(DEGREES_TO_RAD*pCandidate->dec);
          x = cosx*cos(DEGREES_TO_RAD*pCandidate->ra);
          y = cosx*sin(DEGREES_TO_RAD*pCandidate->ra);
          if (kd_insert3(candidate_ptree,x,y,z,pCandidate) != 0) {
            printf("ERROR: fatal return from kd_insert3\n");
            exit(-1);
          }
        
          pCandidate->flarecandidate_index = candidate_index;
#ifdef LOS_DEBUG1
          if (pCandidate->REFNumber == LOS_REF_NUMBER) {
            printf("Line %d at REFNumber %lld\n",__LINE__,pCandidate->REFNumber);
          }
#endif

          candidate_index++;
          if (candidate_index > (candidate_count+interesting_count)) {
            printf("ERROR candidate_index too large\n");
            exit(-1);
          }
        }
      }
    }
    fclose(candidate_handle);
    candidate_handle = NULL;
    printf("Found %d candidates in %d lines %d seconds for iteration %d\n",candidate_index,lineCounter,curTime,iteration);

  } /* end of iteration loop */
  time(&curTime);
  curTime -= startTime;
  candidate_count = candidate_index;
  /* Search the candidates table for any reference to the unrecognized strings in the "comments" column */
  for (candidate_index = 0; candidate_index < candidate_count; candidate_index++) {
    pCandidate = &candidate_table[candidate_index];
#ifdef LOS_DEBUG1
    if (pCandidate->REFNumber == LOS_REF_NUMBER) {
      printf("Line %d at REFNumber %lld\n",__LINE__,pCandidate->REFNumber);
    }
#endif
    for (unrecognized_index = 0; unrecognized_index < unrecognized_count; unrecognized_index++) {
      pUnrecognized = &unrecognized_table[unrecognized_index];
      pChar = strstr(pCandidate->comment,pUnrecognized->name);
      if (pChar != NULL) {
        /* Merge these two entries */
        if (compareStrings(pUnrecognized->name,pChar) == 0) {
          printf("ERROR: unrecognized string %s is in %s # %s\n",pUnrecognized->name,pCandidate->name,pCandidate->comment);
          pUnrecognized->unrecognized = 0;
          pUnrecognized->ra = pCandidate->ra;
          pUnrecognized->dec = pCandidate->dec;
          if ((strlen(pCandidate->comment) + strlen(pUnrecognized->comment) +2) > MAX_COMMENT_SIZE) {
            printf("ERROR MAX_COMMENT_SIZE exceeded by %d for unrecognized name %s\n",(strlen(pCandidate->comment) + strlen(pUnrecognized->comment) +2),pUnrecognized->name);
          } else {
            strcat(pCandidate->comment,"; ");
            strcat(pCandidate->comment,pUnrecognized->comment);
          }
          break;
        }
      }
    }
  }  


         
  for (candidate_index = 0; candidate_index < candidate_count; candidate_index++) {
    pCandidate = &candidate_table[candidate_index];
#if 0
    fprintf(out_handle,"%s # %s\n",pCandidate->name,pCandidate->comment);
#else
    fprintf(out_handle,"%.4f %.4f %s # %s\n",pCandidate->ra,pCandidate->dec,pCandidate->name,pCandidate->comment);
#endif
  }      

  time(&curTime);
  curTime -= startTime;
  printf("Resolved unrecognized positions in %d seconds\n",curTime);

  if (thumbnail_name[0] != 0) {
    /* First find out how many thumbnails we have */
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,thumbnail_handle);
      if (inBuffer == NULL) {
        break;
      }
      lineLen = strlen(inBuffer);
      /* Trim off the carriage return */
      if (inBuffer[lineLen-1] == 10) {
        inBuffer[lineLen-1] = 0;
        lineLen--;
      }
      strcpy(copyLine,inBuffer);
     
      lineCounter++;
    }
    fclose(thumbnail_handle);
    thumbnail_count = lineCounter;
    thumbnail_index = 0;
    thumbnail_table = (PTHUMBNAIL)calloc(thumbnail_count+1,sizeof(THUMBNAIL));
    if (thumbnail_table == NULL) {
      printf("ERROR: failed to allocate the thumbnail table\n");
      exit(-1);
    }

    thumbnail_handle = fopen(thumbnail_name,"rt");
    if (thumbnail_handle == NULL) {
      errorFlag = 1;
      printf("Could not open list file %s\n",thumbnail_name);
    }
    lineCounter = 0;
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,thumbnail_handle);
      if (inBuffer == NULL) {
        break;
      }
      lineLen = strlen(inBuffer);
      /* Trim off the carriage return */
      if (inBuffer[lineLen-1] == 10) {
        inBuffer[lineLen-1] = 0;
        lineLen--;
      }
      /* Trim trailing spaces */
      while ((lineLen >= 0) && (inBuffer[lineLen-1] == ' ')) {
        inBuffer[lineLen-1] = 0;
        lineLen--;
      }
        
      lineCounter++;
#ifdef LOS_DEBUG1
      if (strstr(inBuffer,LOS_REF) != NULL) {
        printf("Line %d at REF %s\n",__LINE__,inBuffer);
      }
#endif
      strcpy(copyLine,inBuffer);
      pThumbnail = &thumbnail_table[thumbnail_index];
      memset(pThumbnail,0,sizeof(THUMBNAIL));
      pThumbnail->candidate_index = -1;
      if (strlen(inBuffer) >= MAX_INPUT_NAME) {
        printf("ERROR: thumbnail filename %s is too long %d\n",inBuffer,strlen(inBuffer));
        exit(-1);
      }
      strcpy(pThumbnail->filename,inBuffer);
      pSlash = strrchr(inBuffer,'/');
      if (pSlash == NULL) {
        printf("ERROR: no directory tree in %s\n",pThumbnail->filename);
        exit(-1);
      }
      *pSlash = 0;
      pSlash++;
      pUnderscore = strrchr(pSlash,'_');
      if (pUnderscore == NULL) {
        printf("ERROR: Thumbnail file %s has no underscore\n",pThumbnail->filename);
        exit(-1);
      }
      *pUnderscore = 0;
      pSign = strchr(pSlash,'+');
      if (pSign != NULL) {
        pSign += 7;
        *pSign = 0;
      } else {
        pSign = strchr(pSlash,'-');
        if (pSign != NULL) {
          pSign += 7;
          *pSign = 0;
        } else {
          /* must be a GSC or Tycho ID */
          pUnderscore = strchr(pSlash,'_');
          if (pUnderscore != NULL) {
            *pUnderscore = 0;
          }
        }
      } 
      if (strlen(pSlash) >= MAX_FLARECANDIDATE_SIZE) {
        printf("ERROR: Thumbnail filename %s is too long %d\n",pSlash,strlen(pSlash));
        exit(-1);
      }

      strcpy(pThumbnail->name,pSlash);
      /* Now look for a subdirectory of the form "yyyy_mm_dd" */
      while (1) {
        pSlash = strrchr(inBuffer,'/');
        if (pSlash == NULL) {
          break;
        }
        *pSlash = 0;
        pSlash++;
        if (strcmp(pSlash,"candidates") == 0) {
          break;
        }
        if (strlen(pSlash) == MAX_SUBDIR_NAME) {
          nvals = sscanf(pSlash,"%d_%d_%d",&yearnumber,&monthnumber,&daynumber);
          if (nvals == 3) {
            strcpy(pThumbnail->subdir,pSlash);
            break;
          }
        }        
      }
      /* Now see if we can find a position */
      if (GetREFNumber(pThumbnail->name,&pThumbnail->REFNumber,&pThumbnail->refType,0,0) == 0) {
        switch(pThumbnail->refType) {
        case REF_TYPE_NONE:
          break;
        case REF_TYPE_GSC:
        case REF_TYPE_TYCHO2:
        case REF_TYPE_GAIA1:
        case REF_TYPE_GAIA2:
        case REF_TYPE_ATLAS2:
          break;
        case REF_TYPE_KEPLER:
          break;
        case REF_TYPE_DASCH:
          if (GetDASCHCoordinates(pThumbnail->name,
                                  &pThumbnail->ra,
                                  &pThumbnail->dec,0,REF_TYPE_DASCH) == 0) {
            pThumbnail->haveCoords = 1;
          }

          break;
        case REF_TYPE_APASS:
          if (GetDASCHCoordinates(pThumbnail->name,
                                  &pThumbnail->ra,
                                  &pThumbnail->dec,0,REF_TYPE_APASS) == 0) {
            pThumbnail->haveCoords = 1;
          }

          break;
        default:
          break;
        }

      } else {
        pThumbnail->haveCoords = 0;
        pThumbnail->REFNumber = 0;
      }
      if (pThumbnail->haveCoords == 0) {
        /* Search the candidate table for the name string */
        for (candidate_index = 0; candidate_index < candidate_count; candidate_index++) {
          pCandidate = &candidate_table[candidate_index];
          if (pCandidate->unrecognized != 0) {
            continue;
          }
          pChar = strstr(pCandidate->comment,pThumbnail->name);
          if (pChar != 0) {
            charVal = pChar[strlen(pThumbnail->name)];
            if ((charVal < '0') || (charVal > '9')) {
              pThumbnail->haveCoords = 1;
              pThumbnail->ra = pCandidate->ra;
              pThumbnail->dec = pCandidate->dec;
              break;
            }
          }
              
        }
      }


      thumbnail_index++;
    }
    thumbnail_count = thumbnail_index;
    printf("found %d thumbnail filenames\n",thumbnail_count);

  }
  
  time(&curTime);
  curTime -= startTime;
  printf("Finished reading thumbnails in %d seconds\n",curTime);



  input_header = table_header(id_handle,TABLE_PARSE);
  if (input_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",id_name);
    return(-1);
  }



  /* Now look for "parseExcess" or its synonym "parseExtra" */
  for (columnIndex = 1; columnIndex <= input_header->header->ncol; columnIndex++) {
    if (strcmp(input_header->header->column[columnIndex],"peakExcess") == 0) {
      gotPeakExcessColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"peakExtra") == 0) {
      gotPeakExtraColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"peakNearbyDistance") == 0) {
      gotPeakNearbyDistanceColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"peakMultipleCount") == 0) {
      gotPeakMultipleCountColumn = 1;
    }
    
    if (strcmp(input_header->header->column[columnIndex],"gsc_bin_index") == 0) {
      gotGscBinIndexColumn = 1;
    }
      
    if (strcmp(input_header->header->column[columnIndex],"peakBadColorCount") == 0) {
      gotPeakBadColorCountColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"peakLimitingYears") == 0) {
      gotPeakLimitingYearsColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"peakDays") == 0) {
      gotPeakDaysColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"peakNoDefectMag") == 0) {
      gotPeakNoDefectMagColumn = 1;
    }
    if (strcmp(input_header->header->column[columnIndex],"peakLongOutburst") == 0) {
      gotPeakLongOutburstColumn = 1;
    }
  }
  if (gotPeakLongOutburstColumn) {
    input_table = table_loadva(id_handle,
                               &input_header,
                               NULL, /* hbase */
                               NULL, /* rows */
                               NULL,
                               sizeof(INPUT),
                               &input_nrecs,
                               TblInt,"nearbyREFflag",TblOff(PINPUT,nearbyREFflag),
                               TblInt,"npoints",TblOff(PINPUT,npoints),
                               TblInt,"ngood",TblOff(PINPUT,ngood),
                               TblInt,"peakCount",TblOff(PINPUT,peakCount),
                               TblInt,"peakCountNF",TblOff(PINPUT,peakCountNF),
                               TblInt,"peakCountWF",TblOff(PINPUT,peakCountWF),
                               TblInt,"peakOutside",TblOff(PINPUT,peakOutside),
                               TblInt,"peakOutsideNF",TblOff(PINPUT,peakOutsideNF),
                               TblInt,"peakOutsideWF",TblOff(PINPUT,peakOutsideWF),
                               TblInt,"peakExtra",TblOff(PINPUT,peakExtra),
                               TblInt,"peakDefectCount",TblOff(PINPUT,peakDefectCount),
                               TblInt,"peakMultipleCount",TblOff(PINPUT,peakMultipleCount),
                               TblInt,"gsc_bin_index",TblOff(PINPUT,gsc_bin_index),
                               TblInt,"peakBadColorCount",TblOff(PINPUT,peakBadColorCount),
                               TblInt,"peakLimitingYears",TblOff(PINPUT,peakLimitingYears),
                               TblInt,"peakLimitingPoints",TblOff(PINPUT,peakLimitingPoints),
                               TblInt,"peakLongOutburst",TblOff(PINPUT,peakLongOutburst),
                               TblDbl,"peakDradRMS2",TblOff(PINPUT,peakDradRMS2),
                               TblDbl,"peakDradRMS3",TblOff(PINPUT,peakDradRMS3),
                               TblDbl,"peakSlope",TblOff(PINPUT,peakSlope),
                               TblDbl,"peakUpperCount",TblOff(PINPUT,peakUpperCount),
                               TblDbl,"peakMaxMag",TblOff(PINPUT,peakMaxMag),
                               TblDbl,"peakMinMag",TblOff(PINPUT,peakMinMag),
                               TblDbl,"peakNoDefectMag",TblOff(PINPUT,peakNoDefectMag),
                               TblDbl,"Stdmag",TblOff(PINPUT,Stdmag),
                               TblDbl,"peakRA",TblOff(PINPUT,peakRA),
                               TblDbl,"peakDec",TblOff(PINPUT,peakDec),
                               TblDbl,"peakYear",TblOff(PINPUT,peakYear),
                               TblDbl,"peakNearbyDistance",TblOff(PINPUT,peakNearbyDistance),
                               TblDbl,"ra",TblOff(PINPUT,ra),
                               TblDbl,"declination",TblOff(PINPUT,declination),
                               TblDbl,"lat",TblOff(PINPUT,lat),
                               TblDbl,"peakDays",TblOff(PINPUT,peakDays),
                               TblBuf,"REF"    ,TblOff(PINPUT,REF),MAX_REF,
                               TblBuf,"nearbyObjects"    ,TblOff(PINPUT,nearbyObjects),MAX_NEARBY_OBJECTS_STRING,
                               0,"end",0);

  } else { 
    if (gotPeakNoDefectMagColumn) {
      input_table = table_loadva(id_handle,
                                 &input_header,
                                 NULL, /* hbase */
                                 NULL, /* rows */
                                 NULL,
                                 sizeof(INPUT),
                                 &input_nrecs,
                                 TblInt,"nearbyREFflag",TblOff(PINPUT,nearbyREFflag),
                                 TblInt,"npoints",TblOff(PINPUT,npoints),
                                 TblInt,"ngood",TblOff(PINPUT,ngood),
                                 TblInt,"peakCount",TblOff(PINPUT,peakCount),
                                 TblInt,"peakCountNF",TblOff(PINPUT,peakCountNF),
                                 TblInt,"peakCountWF",TblOff(PINPUT,peakCountWF),
                                 TblInt,"peakOutside",TblOff(PINPUT,peakOutside),
                                 TblInt,"peakOutsideNF",TblOff(PINPUT,peakOutsideNF),
                                 TblInt,"peakOutsideWF",TblOff(PINPUT,peakOutsideWF),
                                 TblInt,"peakExtra",TblOff(PINPUT,peakExtra),
                                 TblInt,"peakDefectCount",TblOff(PINPUT,peakDefectCount),
                                 TblInt,"peakMultipleCount",TblOff(PINPUT,peakMultipleCount),
                                 TblInt,"gsc_bin_index",TblOff(PINPUT,gsc_bin_index),
                                 TblInt,"peakBadColorCount",TblOff(PINPUT,peakBadColorCount),
                                 TblInt,"peakLimitingYears",TblOff(PINPUT,peakLimitingYears),
                                 TblInt,"peakLimitingPoints",TblOff(PINPUT,peakLimitingPoints),
                                 TblDbl,"peakDradRMS2",TblOff(PINPUT,peakDradRMS2),
                                 TblDbl,"peakDradRMS3",TblOff(PINPUT,peakDradRMS3),
                                 TblDbl,"peakSlope",TblOff(PINPUT,peakSlope),
                                 TblDbl,"peakUpperCount",TblOff(PINPUT,peakUpperCount),
                                 TblDbl,"peakMaxMag",TblOff(PINPUT,peakMaxMag),
                                 TblDbl,"peakMinMag",TblOff(PINPUT,peakMinMag),
                                 TblDbl,"peakNoDefectMag",TblOff(PINPUT,peakNoDefectMag),
                                 TblDbl,"Stdmag",TblOff(PINPUT,Stdmag),
                                 TblDbl,"peakRA",TblOff(PINPUT,peakRA),
                                 TblDbl,"peakDec",TblOff(PINPUT,peakDec),
                                 TblDbl,"peakYear",TblOff(PINPUT,peakYear),
                                 TblDbl,"peakNearbyDistance",TblOff(PINPUT,peakNearbyDistance),
                                 TblDbl,"ra",TblOff(PINPUT,ra),
                                 TblDbl,"declination",TblOff(PINPUT,declination),
                                 TblDbl,"lat",TblOff(PINPUT,lat),
                                 TblDbl,"peakDays",TblOff(PINPUT,peakDays),
                                 TblBuf,"REF"    ,TblOff(PINPUT,REF),MAX_REF,
                                 TblBuf,"nearbyObjects"    ,TblOff(PINPUT,nearbyObjects),MAX_NEARBY_OBJECTS_STRING,
                                 0,"end",0);
    } else {
      printf("ERROR: obsolete table without the peakNoDefectMagColumn or peakLongOutburst columns\n");
      exit(-1);
    }
  } 
  time(&curTime);
  curTime -= startTime;
  printf("read %d records for %s in %d seconds\n",input_nrecs,id_name,curTime);

  time(&curTime);
  curTime -= startTime;
  printf("Combined input tables %d records in %d seconds\n",input_nrecs,curTime);
  qsort((void *)input_table,input_nrecs,sizeof(INPUT),REFCompare);
  for (input_index = 1; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
    pInput2 = &input_table[input_index-1];
    workDone = 1;
    while (workDone != 0) {
      workDone = 0;
      if (strcmp(pInput->REF,pInput2->REF) == 0) {
        workDone = 1;
        for (input_index3 = input_index+1;input_index3 < input_nrecs; input_index3++) {
          pInput3 = &input_table[input_index3];
          pInput4 = &input_table[input_index3-1];
          memcpy(pInput4,pInput3,sizeof(INPUT));
        }
        input_nrecs--;
        inputDupCount++;
      }
    }
  }


  time(&curTime);
  curTime -= startTime;
  printf("Eliminate %d duplicates %d records in %d seconds\n",inputDupCount,input_nrecs,curTime);
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
#ifdef LOS_DEBUG1
    if (strcmp(pInput->REF,LOS_REF) == 0) {
      printf("Line %d at REF %s\n",__LINE__,pInput->REF);
    }
#endif
    pInput->minDistance = -1;
    pInput->candidate_index = -1;
    pInput->rejectFlags = 0;
    pInput->nearbyObjectDist = 999.;
    if (gotPeakLimitingYearsColumn == 0) {
      pInput->peakLimitingYears = 0;
      pInput->peakLimitingPoints = 0;
    }
    if (gotPeakBadColorCountColumn == 0) {
      pInput->peakBadColorCount = 0;
    }
    if (gotPeakNoDefectMagColumn == 0) {
      pInput->peakNoDefectMag = pInput->peakMaxMag;
    }
    if (gotPeakLongOutburstColumn == 0) {
      pInput->peakLongOutburst = 0;
    }

    if (gotPeakDaysColumn == 0) {
      pInput->peakDays = 999.;
    }
    if (gotGscBinIndexColumn == 0) {
      pInput->gsc_bin_index = -1;
      pInput->peakDradRMS2 = 99;
      pInput->peakDradRMS3 = 99;
      pInput->Sig_TC = 99.0;
    } else {
      pInput->Sig_TC = 99.0;
      if (pInput->peakDradRMS2 < 90.0) {
        if (pInput->peakDradRMS3 < 90.0) {
          pInput->Sig_TC = fmax(pInput->peakDradRMS2,pInput->peakDradRMS3);
        } else {
          pInput->Sig_TC = pInput->peakDradRMS2;
        }
      } else {
        if (pInput->peakDradRMS3 < 90.0) {
          pInput->Sig_TC = pInput->peakDradRMS3;
        }
      }
    }
    if (gotPeakMultipleCountColumn == 0) {
      pInput->peakMultipleCount = 0;
    }
    if (gotPeakNearbyDistanceColumn == 0) {
      pInput->peakNearbyDistance = 0;
    } else {
      pInput->nearbyObjectDist = pInput->peakNearbyDistance;
    }
    if (gotPeakExtraColumn == 0) {
      pInput->peakDec = 99.0;
      pInput->peakRA = 9999.0;
      pInput->peakOffset = 0.0;
      pInput->peakCountNF = -1;
      pInput->peakCountWF = -1;
      pInput->peakDefectCount = -1;
      pInput->peakOutside = -1;
      pInput->peakOutsideWF = -1;
      pInput->peakOutsideNF = -1;
    } else {
      factor = cos(((pInput->declination+pInput->peakDec)/2.0) * DEGREES_TO_RAD);
      pInput->peakOffset = 3600.0*sqrt(sqr(factor*(pInput->ra - pInput->peakRA)) + sqr(pInput->declination-pInput->peakDec));
    }
    if (pInput->peakCount < LOW_TC_POINTS) {
      pInput->rejectFlags |= REJECT_FLAG_LOWPOINTS;        
    }
    if ((gotPeakExcessColumn == 0) && (gotPeakExtraColumn == 0)) {
      pInput->ngood = 0;
      pInput->peakExtra = 0;
    } else {
      if (pInput->peakExtra > MAX_PEAK_EXTRA_COUNT) {
        pInput->rejectFlags |= REJECT_FLAG_EXTRA;        
      }
    }
    if (gotPeakMultipleCountColumn) {
#ifdef LOS_DEBUG1
      if (strcmp(pInput->REF,LOS_REF) == 0) {
        printf("Line %d at REF %s\n",__LINE__,pInput->REF);
      }
#endif /* LOS_DEBUG1 */
#if 1
      if ((100*pInput->peakMultipleCount) > (MAX_MULTIPLE_PERCENTAGE_33 * (pInput->peakCountWF+pInput->peakCountNF))) {
        pInput->rejectFlags |= REJECT_FLAG_MULTIPLE33;
      }

#else
#if 0
      if (pInput->peakMultipleCount > 1) {
        pInput->rejectFlags |= REJECT_FLAG_MULTIPLE33;
      }
#endif
      if (((pInput->peakCountWF+pInput->peakCountNF-pInput->peakMultipleCount) < MIN_TC_POINTS) && 
          ((100*pInput->peakMultipleCount) > (MAX_MULTIPLE_PERCENTAGE_33 * (pInput->peakCountWF+pInput->peakCountNF)))) {
        pInput->rejectFlags |= REJECT_FLAG_MULTIPLE33;
        
      }
#endif 
      if (((pInput->peakCountWF+pInput->peakCountNF-pInput->peakMultipleCount) < MIN_TC_POINTS) &&
          ((100*pInput->peakMultipleCount) > (MAX_MULTIPLE_PERCENTAGE_90 * (pInput->peakCountWF+pInput->peakCountNF)))) {
        pInput->rejectFlags |= REJECT_FLAG_MULTIPLE90;
      } 
      if (gotPeakNearbyDistanceColumn) {
        if ((100*(pInput->peakDefectCount+pInput->peakMultipleCount)) > (MAX_DEFECT_PERCENTAGE *(pInput->peakCountWF+pInput->peakCountNF))) {
          pInput->rejectFlags |= REJECT_FLAG_MULTDEFECT;
        }
      }


    }
    if ((pInput->peakMinMag - pInput->peakNoDefectMag) <  LARGE_TC_MAGNITUDE) {
      pInput->rejectFlags |= REJECT_FLAG_LOWRANGE;
    }
    if (gotPeakLimitingYearsColumn) {
      if (pInput->peakLimitingYears <  MIN_TC_LIMITING_YEARS) {
        pInput->rejectFlags |= REJECT_FLAG_LIM_YEARS;
      }
      if (pInput->peakLimitingPoints <  MIN_TC_LIMITING_POINTS) {
        pInput->rejectFlags |= REJECT_FLAG_LIM_POINTS;
      }
    }
    if (gotPeakBadColorCountColumn) {
      if ((100*(pInput->peakBadColorCount)) > (MAX_COLORTERM_PERCENTAGE *(pInput->peakCountWF+pInput->peakCountNF))) {
        pInput->rejectFlags |= REJECT_FLAG_COLORTERM;
      }
    }
    
    if (pInput->lat < MAX_GALACTIC_LATITUDE) {
      pInput->rejectFlags |= REJECT_FLAG_UNRELEASED;
    }

    for (exclusionIndex = 0; exclusionIndex < exclusionTableSize; exclusionIndex++) {
      pExclusion = &exclusionTable[exclusionIndex];
      if ((pInput->ra > pExclusion->ra1) &&
          (pInput->ra < pExclusion->ra2) &&
          (pInput->declination > pExclusion->dec1) &&
          (pInput->declination < pExclusion->dec2)) {
        pInput->rejectFlags |= REJECT_FLAG_EXCLUSION;
        break;
      }       
    }

    if (gotPeakNearbyDistanceColumn) {
      if ((100*pInput->peakDefectCount) > (MAX_DEFECT_PERCENTAGE *(pInput->peakCountWF+pInput->peakCountNF))) {
        pInput->rejectFlags |= REJECT_FLAG_DEFECT;
      }
    }
    if (pInput->peakLongOutburst > 0) {
      pInput->rejectFlags |= REJECT_FLAG_LONGFLARE;
    }
    /* Split the nearbyObjects field into individual records */
    pInput->nearbyCount = 0;
    memset(pInput->nearbyImage,0,sizeof(pInput->nearbyImage));
    pInput->rejectNearbyDistance = 999.0;
    strcpy(nearbyObjects,pInput->nearbyObjects);
    pChar = nearbyObjects;
    while (pChar != NULL) {
      pChar2 = strchr(pChar,';');
      if (pChar2 != NULL) {
        *pChar2 = 0;
        pChar2++;
      }
      while (pChar[0] == ' ') {
        pChar++;
      }
      if (pChar[0] == 0) {
        break;
      }
      if (strlen(pChar) > MAX_NEARBY_STRING) {
        printf("ERROR: String size %d exceeds MAX_NEARBY_STRING for '%s'\n",strlen(pChar),pChar);
        exit(-1);
      }
      if (pInput->nearbyCount >= MAX_NEARBY_COUNT) {
        printf("ERROR: MAX_NEARBY_COUNT exceeded for %s\n",pInput->nearbyObjects);
        exit(-1);
      }
      pNearbyImage = &pInput->nearbyImage[pInput->nearbyCount];
      strcpy(pNearbyImage->nearbyObject,pChar);
#if 0
      printf("NEARBY %d %s\n",pInput->nearbyCount,pNearbyImage->nearbyObject);
#endif
      pInput->nearbyCount++;
      pChar = pChar2;

    }
    /* Now parse the individual entries */
    for (nearbyIndex = 0; nearbyIndex < pInput->nearbyCount; nearbyIndex++) {
      pNearbyImage = &pInput->nearbyImage[nearbyIndex];
      strcpy(nearbyObjects,pNearbyImage->nearbyObject);
      nearbyLength = strlen(nearbyObjects);
      pChar = &nearbyObjects[nearbyLength-1];
      if (*pChar == '*') {
        pNearbyImage->asteriskFlag = 1;
        *pChar = 0;
      }
      pChar = strrchr(nearbyObjects,'@');
      if (pChar == NULL) {
        pNearbyImage->distance = -1;
      } else {
        *pChar = 0;
        pChar++;
        while (pChar[0] == ' ') {
          pChar++;
        }
        nvals = sscanf(pChar,"%lf",&pNearbyImage->distance);
        if (nvals != 1) {
          printf("ERROR: failed to decode distance in %s %s\n",pNearbyImage->nearbyObject,pInput->nearbyObjects);
          exit(-1);
        }
        if (pNearbyImage->distance < pInput->nearbyObjectDist) {
          pInput->nearbyObjectDist = pNearbyImage->distance;
        }
      }
      pChar = strrchr(nearbyObjects,'(');
      if (pChar == NULL) {
        pNearbyImage->catalogmag = 99.0;
      } else {
        *pChar = 0;
        pChar++;
        while (pChar[0] == ' ') {
          pChar++;
        }
        nvals = sscanf(pChar,"%lf",&pNearbyImage->catalogmag);
        if (nvals != 1) {
          printf("ERRORD: failed to decode catalogmag in %s %s\n",pNearbyImage->nearbyObject,pInput->nearbyObjects);
          pNearbyImage->catalogmag = 99;
        }
      }
      nearbyLength = strlen(nearbyObjects);
      nearbyLength--;
      pChar = &nearbyObjects[nearbyLength];
      while ((*pChar == ' ') &&
             nearbyLength >= 0) {
        *pChar = 0;
        pChar--;
        nearbyLength--;
      }
      /* CSS 090910:223418-035530  is a special exception */
      if ((strstr(nearbyObjects,"CSS ") == nearbyObjects) &&
          (nearbyObjects[10] == ':')) {
        pChar2 = &nearbyObjects[11];
      } else {
        pChar2 = nearbyObjects;
      }
      pChar = strchr(pChar2,':');
      if (pChar != NULL) {
        *pChar = 0;
        pChar++;
        while (*pChar == ' ') {
          *pChar = 0;
        }
        if (strlen(pChar) > MAX_NEARBY_TYPE_STRING) {
          printf("ERROR: string %d '%s' exceeds MAX_NEARBY_TYPE_STRING in '%s'\n",strlen(pChar),pChar,pNearbyImage->nearbyObject);
          exit(-1);
        }
        strcpy(pNearbyImage->type,pChar);
      }
      pChar = strrchr(pChar2,':');
      if (pChar != NULL) {
        printf("ERROR: Object name has a colon in '%s' '%s'\n",nearbyObjects,pNearbyImage->nearbyObject);
        exit(-1);
      }
      nearbyLength = strlen(nearbyObjects);
      nearbyLength--;
      pChar = &nearbyObjects[nearbyLength];
      while ((*pChar == ' ') &&
             nearbyLength >= 0) {
        *pChar = 0;
        pChar--;
        nearbyLength--;
      }
      if (strlen(nearbyObjects) > MAX_IDENTIFIER_STRING) {
        printf("ERROR: string %d '%s' exceeds MAX_IDENTIFIER_STRING in '%s'\n",strlen(nearbyObjects),nearbyObjects,pNearbyImage->nearbyObject);
        exit(-1);
      }
      strcpy(pNearbyImage->identifier,nearbyObjects);
#if 0
      printf("NEARBY '%25s' '%15s' (%4.1f) @ %6.1f *:%d for '%s'\n",
             pNearbyImage->identifier,
             pNearbyImage->type,
             pNearbyImage->catalogmag,
             pNearbyImage->distance,
             pNearbyImage->asteriskFlag,
             pNearbyImage->nearbyObject);
#endif
      /* Now figure out what we have */
#if 0
      if (strcmp(pInput->REF,"DASCH_J105230.7-180408") == 0) {
        printf("At ref %s\n",pInput->REF);
      }
#endif

      pNearbyImage->nearbyType = NEARBY_TYPE_UNKNOWN;
      if (GetREFNumber(pNearbyImage->identifier,&REFNumber,&refType,0,0) == 0) {
        pNearbyImage->nearbyType = NEARBY_TYPE_CATALOG;
        /* Is the brightest point of the transient within two magnitudes of the catalog mag */
        if (gotPeakNearbyDistanceColumn) {
          if ((pInput->peakNearbyDistance > 0) && (pInput->peakNearbyDistance < MAX_NEARBY_ARCSEC)) {
            pInput->rejectFlags |= REJECT_FLAG_NEARBY;
            pInput->rejectNearbyDistance = pInput->peakNearbyDistance;
          }
        } else {
          if ((fabs(pNearbyImage->catalogmag-pInput->peakMaxMag) < MAX_TC_NEARBY_MAG) &&
              (pNearbyImage->distance < MAX_NEARBY_ARCSEC)) {
            pInput->rejectFlags |= REJECT_FLAG_NEARBY;
            pInput->rejectNearbyDistance = pNearbyImage->distance;

          }
          /* Is the average brightness  of the transient within two magnitudes of the catalog mag */
          if ((fabs(pNearbyImage->catalogmag-((pInput->peakMaxMag+pInput->peakMinMag)/2.0)) < MAX_TC_NEARBY_MAG) &&
              (pNearbyImage->distance < MAX_NEARBY_ARCSEC)) {
            pInput->rejectFlags |= REJECT_FLAG_NEARBY;
            pInput->rejectNearbyDistance = pNearbyImage->distance;
          }
        }
      }
      if (pNearbyImage->nearbyType == NEARBY_TYPE_UNKNOWN) {
        if (strstr(pNearbyImage->type,"SDSSg") != NULL) {
          pNearbyImage->nearbyType = NEARBY_TYPE_SDSS;
          /* Now parse out the SDSS characteristics */
          sdssPtr = pNearbyImage->type+5;
          nvals = sscanf(sdssPtr,"%lfc%lf",&pNearbyImage->catalogmag,&pNearbyImage->color);
          if (nvals != 2) {
            printf("ERROR parsing SDSS type '%s' '%s'\n",pNearbyImage->type,pNearbyImage->nearbyObject);
            exit(-1);
          }
          sdssDistance = ((MAX_SDSS_DISTANCE_FACTOR*pInput->peakDradRMS3)+pInput->peakOffset);
          if (sdssDistance > MAX_SDSS_DISTANCE) {
            sdssDistance = MAX_SDSS_DISTANCE;
          }
          if ((pNearbyImage->catalogmag >= MAX_SDSS_MAGNITUDE) && 
              (pNearbyImage->distance < sdssDistance)) {
            pInput->rejectFlags |= REJECT_FLAG_SDSS;
            pNearbyImage->sdssrejectflag = 1;
          }
        }
      }
      if (pNearbyImage->nearbyType == NEARBY_TYPE_UNKNOWN) {
        if (strstr(pNearbyImage->identifier,"NGC") != NULL) {
          pNearbyImage->nearbyType = NEARBY_TYPE_NGC;
          pInput->rejectFlags |= REJECT_FLAG_NGC;
        }
      }
      if (pNearbyImage->nearbyType == NEARBY_TYPE_UNKNOWN) {
        if ((strstr(pNearbyImage->type,"UZC Galaxy") != NULL) ||
            (strstr(pNearbyImage->type,"PGC Galaxy") != NULL) ||
            (strstr(pNearbyImage->type,"SGC Galaxy") != NULL) ||
            (strstr(pNearbyImage->type,"UGC Galaxy") != NULL))
          
          {
            pNearbyImage->nearbyType = NEARBY_TYPE_GALAXY;
          }
      }
      if (pNearbyImage->nearbyType == NEARBY_TYPE_UNKNOWN) {
        if (strstr(pNearbyImage->identifier,"IC") == pNearbyImage->identifier) {
          pNearbyImage->nearbyType = NEARBY_TYPE_IC;
        }
      }

      if (pNearbyImage->nearbyType == NEARBY_TYPE_UNKNOWN) {
        pNearbyImage->nearbyType = NEARBY_TYPE_VSX;
        if (strstr(pNearbyImage->type,"UG") == pNearbyImage->type) {
          pInput->rejectFlags |= REJECT_FLAG_UG;
        } else if (pNearbyImage->type[0] == 'M') {
          pInput->rejectFlags |= REJECT_FLAG_MIRA;
        } else if (strstr(pNearbyImage->type,"RR") == pNearbyImage->type) {
          pInput->rejectFlags |= REJECT_FLAG_RR;
        } else if (pNearbyImage->type[0] == 'N') {
          pInput->rejectFlags |= REJECT_FLAG_NOVA;
        } else  if ((strstr(pNearbyImage->type,"XPRM") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"XM") == pNearbyImage->type)) {
          pInput->rejectFlags != REJECT_FLAG_POLAR;
        } else if (((strstr(pNearbyImage->type,"E") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"EA") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"EB") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"EW") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"GS") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"PN") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"RS") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"WD") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"WR") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"AR") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"D") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"DM") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"DS") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"DW") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"K") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"KE") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"KW") == pNearbyImage->type) ||
                    (strstr(pNearbyImage->type,"SD") == pNearbyImage->type)) &&
                   (strstr(pNearbyImage->type,"EP") == NULL) &&
                   (strstr(pNearbyImage->type,"EL") == NULL) &&
                   (strstr(pNearbyImage->type,"DCEPS") == NULL) &&
                   (strstr(pNearbyImage->type,"DCEP") == NULL) &&
                   (strstr(pNearbyImage->type,"DSCTC") == NULL) &&
                   (strstr(pNearbyImage->type,"DSCT") == NULL)) {
          pInput->rejectFlags |= REJECT_FLAG_ECLIPSING;
        } else { 
          pInput->rejectFlags |= REJECT_FLAG_VAR;
        }
#if 0
        printf("NEARBY '%-15s' %d '%-25s'  (%4.1f) @ %6.1f  *:%d for '%s'\n",
               pNearbyImage->type,
               pNearbyImage->nearbyType,
               pNearbyImage->identifier,
               pNearbyImage->catalogmag,
               pNearbyImage->distance,
               pNearbyImage->asteriskFlag,
               pNearbyImage->nearbyObject);
#endif
      }
#if 0
      printf("NEARBY %d '%-25s' '%-15s' (%4.1f) @ %6.1f *:%d for '%s'\n",
             pNearbyImage->nearbyType,
             pNearbyImage->identifier,
             pNearbyImage->type,
             pNearbyImage->catalogmag,
             pNearbyImage->distance,
             pNearbyImage->asteriskFlag,
             pNearbyImage->nearbyObject);
#endif
    }


    
  }
#ifdef PREFILTER_MULTDEFECT
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
#if 0
    if (strcmp(pInput->REF,"DASCH_J174033.6+414756") == 0) {
      printf("At %s\n",pInput->REF);
    }
#endif
    if ((pInput->rejectFlags & PREFILTER_MULTDEFECT_MASK) != 0) {    
      if ((pInput->rejectFlags & (~PREFILTER_MULTDEFECT_MASK)) != 0) {
        pInput->rejectFlags &= PREFILTER_MULTDEFECT_MASK;
      }
    }
  }

#endif /* PREFILTER_MULTDEFECT */
  /* Now do a search for similar dates */
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
#ifdef PREFILTER_MULTDEFECT
    if ((pInput->rejectFlags & PREFILTER_MULTDEFECT_MASK) != 0) {
      continue;
    }
#endif /* PREFILTER_MULTDEFECT */
    for (input_index2 = input_index+1; input_index2 < input_nrecs; input_index2++) {
      pInput2 = &input_table[input_index2];
#ifdef PREFILTER_MULTDEFECT
      if ((pInput2->rejectFlags & PREFILTER_MULTDEFECT_MASK) != 0) {
        continue;
      }
#endif /* PREFILTER_MULTDEFECT */
      if (fabs(pInput->peakYear-pInput2->peakYear) < MIN_TC_DATERANGE) {
        if ((pInput->rejectFlags &  REJECT_FLAG_SAMETIME) == 0) {
          if (showsametimematches) {
            printf("sametime match of %6.3f hours between %s and %s\n",
                   fabs(pInput->peakYear-pInput2->peakYear),
                   pInput->REF,
                   pInput2->REF);
          }
          sametimeCount++;
        }
        if ((pInput2->rejectFlags &  REJECT_FLAG_SAMETIME) == 0) {
          sametimeCount++;
        }
        pInput->rejectFlags |= REJECT_FLAG_SAMETIME;
        pInput2->rejectFlags |= REJECT_FLAG_SAMETIME;
      }
    }
  }



  time(&curTime);
  curTime -= startTime;
  printf("Parsed input tables in %d seconds\n",curTime);

  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
#ifdef LOS_DEBUG1
    if (strcmp(pInput->REF,LOS_REF) == 0) {
      printf("Line %d at REF %s\n",__LINE__,pInput->REF);
    }
#endif /* LOS_DEBUG1 */

#if 0
    if (strcmp("N02112319798",pInput->REF) == 0) {
      printf("At %s\n",pInput->REF);
    }
#endif
    /* Search the unrecognized candidate table for this one */
    for (unrecognized_index = 0; unrecognized_index < unrecognized_count; unrecognized_index++) {
      pUnrecognized = &unrecognized_table[unrecognized_index];
      
      if ((pUnrecognized->unrecognized > 0) && 
          (strcmp(pUnrecognized->name,pInput->REF) == 0)) {
        printf("RECOGNIZED %.4f %.4f # %s %s\n",pInput->ra,pInput->declination,pUnrecognized->name,pUnrecognized->comment);
        pUnrecognized->unrecognized = 0;
        pUnrecognized->ra = pInput->ra;
        pUnrecognized->dec = pInput->declination;
        break;
      }

    }

    
    if ((enforcePeakExcess > 0) && 
        ((gotPeakExtraColumn != 0) || (gotPeakExcessColumn != 0)) &&
        (pInput->peakExtra > MAX_TC_EXCESS_POINTS)) {
      continue;
    }
#if 0
    for (candidate_index = 0; candidate_index < candidate_count; candidate_index++) {
      pCandidate = &candidate_table[candidate_index];
      distance = 3600.*wcsdist(pCandidate->ra,pCandidate->dec,pInput->ra,pInput->declination);
      if (candidate_index == 0) {
        pInput->candidate_index = 0;
        pInput->previousFlag = pCandidate->previousFlag;
        pInput->minDistance = distance;
      } else {
        if (distance < pInput->minDistance) {
          pInput->candidate_index = candidate_index;
          pInput->previousFlag = pCandidate->previousFlag;
          pInput->minDistance = distance;
        }
      }
    }
#if 0
    if (strstr(pCandidate->comment,"Possible strong cand. Not on DSS so B >20") != NULL) {
      printf("line %d input_index %d ra %f dec %f candidate_index %d comment %s\n",__LINE__,input_index,pInput->ra,pInput->declination,pInput->candidate_index,pCandidate->comment);
    }
#endif
#else
#ifdef LOS_DEBUG1
    if (strcmp(pInput->REF,LOS_REF) == 0) {
      printf("Line %d at REF %s\n",__LINE__,pInput->REF);
    }
#endif /* LOS_DEBUG1 */
    /* Query the kd tree for a nearby candidate */
    z = sin(DEGREES_TO_RAD*pInput->declination);
    cosx = cos(DEGREES_TO_RAD*pInput->declination);
    x = cosx*cos(DEGREES_TO_RAD*pInput->ra);
    y = cosx*sin(DEGREES_TO_RAD*pInput->ra);
    pt[0] = x;
    pt[1] = y;
    pt[2] = z;
    presults = kd_nearest_range(candidate_ptree, pt,(MAX_DISTANCE*DEGREES_TO_RAD)/3600.);
    resultSize = kd_res_size(presults);
    foundMatchFlag = 0;
    while  (!kd_res_end( presults ) ) {
      pCandidate = (PFLARECANDIDATE)kd_res_item(presults,pos);
      factor = cos(DEGREES_TO_RAD* pCandidate->dec);
      deltaDec = pCandidate->dec - pInput->declination;
      if (deltaDec < 0) {
        deltaDec = - deltaDec;
      }
      deltaRA  = factor *(pCandidate->ra  - pInput->ra);
      if (deltaRA < 0) {
        deltaRA = -deltaRA;
      }
      distance = sqrt(sqr(deltaDec)+sqr(deltaRA)) * 3600.;
      if (foundMatchFlag == 0) {
        pInput->candidate_index = pCandidate->flarecandidate_index;
        pInput->previousFlag = pCandidate->previousFlag;
        pInput->minDistance = distance;
      } else {
        if (distance < pInput->minDistance) {
          pInput->candidate_index = pCandidate->flarecandidate_index;
          pInput->previousFlag = pCandidate->previousFlag;
          pInput->minDistance = distance;
        }
      }      
    
      foundMatchFlag = 1;

      /* go to the next entry */
      kd_res_next( presults );
            
    }
        
    kd_res_free(presults);

#endif


    if (printSummaryTable) {
      printf("REF %22s ra: %8.4f dec: %8.4f mag: %4.1f npoints: %2d ngood: %d peakCount %d NF %d WF %d peakOutside %d NF %d WF %d peakExtra %d %s",
             pInput->REF,
             pInput->ra,
             pInput->declination,
             pInput->peakMaxMag,
             pInput->npoints,
             pInput->ngood,
             pInput->peakCount,
             pInput->peakCountNF,
             pInput->peakCountWF,
             pInput->peakOutside,
             pInput->peakOutsideNF,
             pInput->peakOutsideWF,
             pInput->peakExtra,
             pInput->nearbyObjects);
    }

#if 0
    if (pInput->Stdmag > 90) {
      printf("ERROR Stdmag %.1f\n",pInput->Stdmag);
    }
#endif
    for (rejectTableIndex = 0; rejectTableIndex < flareRejectTableSize; rejectTableIndex++) {
      pRejectTable = &flareRejectTable[rejectTableIndex];
      if ((pRejectTable->flagmask & pInput->rejectFlags) != 0) {
        if (printSummaryTable) {
          printf(" %s",pRejectTable->flagname);
        }
      }
    }
    if (foundMatchFlag != 0) {
      pCandidate = &candidate_table[pInput->candidate_index];
      if ((pInput->minDistance >= 0) && (pInput->minDistance < MAX_DISTANCE)) {
        pCandidate->input_index = input_index;
        if (printSummaryTable) {
          printf(" ALREADY CHECKED: %f.1 arcsec from %s %s",
                 pInput->minDistance,
                 pCandidate->name,
                 pCandidate->comment);
        }
        if (pInput->peakCount > pCandidate->prevPeakCount) {
          if (pCandidate->prevPeakCount > 0) {
            pInput->rejectFlags |= REJECT_FLAG_MOREPOINTS;
          }
        }
        pCandidate->newPeakCount = pInput->peakCount;

        switch(pCandidate->typeFlag) {
        case 'O':
          pInput->rejectFlags |= REJECT_FLAG_ORANGE;
          break;
        case 'Y':
          pInput->rejectFlags |= REJECT_FLAG_YELLOW;
          break;
        case 'N':
          pInput->rejectFlags |= REJECT_FLAG_NOTTC;
          break;
        default:
          pInput->rejectFlags |= REJECT_FLAG_PREVIOUS;
          break;
        }
      }
    }
    if (printSummaryTable) {
      printf("\n");
    }
  }

  time(&curTime);
  curTime -= startTime;
  printf("Printed summary table in %d seconds\n",curTime);

  /* Find out the distribution of adjacent points */
  gsc_bin_count_single = (int *)calloc(pGscBin->total_gsc_bins,sizeof(int));
  if (gsc_bin_count_single == NULL) {
    printf("ERROR: failed to allocate gsc_bin_count_single of size %d\n",pGscBin->total_gsc_bins);
    exit(-1);
  }
  gsc_bin_count_multiple = (int *)calloc(pGscBin->total_gsc_bins,sizeof(int));
  if (gsc_bin_count_multiple == NULL) {
    printf("ERROR: failed to allocate gsc_bin_count_multiple of size %d\n",pGscBin->total_gsc_bins);
    exit(-1);
  }

  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
    pInput->input_index = input_index;
    /* Insert this one into the input ktree */
    z = sin(DEGREES_TO_RAD*pInput->declination);
    cosx = cos(DEGREES_TO_RAD*pInput->declination);
    x = cosx*cos(DEGREES_TO_RAD*pInput->ra);
    y = cosx*sin(DEGREES_TO_RAD*pInput->ra);
    if (kd_insert3(input_ptree,x,y,z,pInput) != 0) {
      printf("ERROR: fatal return from kd_insert3\n");
      exit(-1);
    }


    if ((pInput->rejectFlags & LOCATION_MASK) == 0) {
      if ((pInput->peakRA < 990.0) && (pInput->peakDec < 99.00)) {
        gsc_bin_index = GetGSCBin(pGscBin,pInput->peakRA,pInput->peakDec,&decBin,&raBin,"parseidtable");
      } else {
        gsc_bin_index = GetGSCBin(pGscBin,pInput->ra,pInput->declination,&decBin,&raBin,"parseidtable");
      }
      gsc_bin_count_single[gsc_bin_index]++;
      gsc_single_candidate_count++;
    }
  }
  /* Now create a distance histogram between all candidates */
  if (printLocationHistogram) {
    for (input_index = 0; input_index < input_nrecs; input_index++) {
      pInput = &input_table[input_index];
      /* Query the kd tree for a nearby candidate */
      if ((pInput->peakRA < 990.0) && (pInput->peakDec < 99.00)) {
        z = sin(DEGREES_TO_RAD*pInput->peakDec);
        cosx = cos(DEGREES_TO_RAD*pInput->peakDec);
        x = cosx*cos(DEGREES_TO_RAD*pInput->peakRA);
        y = cosx*sin(DEGREES_TO_RAD*pInput->peakRA);
      } else{
        z = sin(DEGREES_TO_RAD*pInput->declination);
        cosx = cos(DEGREES_TO_RAD*pInput->declination);
        x = cosx*cos(DEGREES_TO_RAD*pInput->ra);
        y = cosx*sin(DEGREES_TO_RAD*pInput->ra);
      }
      pt[0] = x;
      pt[1] = y;
      pt[2] = z;
      foundMatchFlag = 0;
      bestdistance = 1.0*MAX_LOCATION_HISTOGRAM_ARCSEC;
      presults = kd_nearest_range(input_ptree, pt,((1.0*MAX_LOCATION_HISTOGRAM_ARCSEC)*DEGREES_TO_RAD)/3600.0);
      resultSize = kd_res_size(presults);
      while  (!kd_res_end( presults ) ) {
        pInput2 = (PINPUT)kd_res_item(presults,pos);
        if (pInput2->input_index > pInput->input_index)  {
          /* consider only indices we have not searched yet to avoid duplicates */
          factor = cos(DEGREES_TO_RAD* pInput2->declination);
          deltaDec = pInput2->declination - pInput->declination;
          if (deltaDec < 0) {
            deltaDec = - deltaDec;
          }
          deltaRA  = factor *(pInput2->ra  - pInput->ra);
          if (deltaRA < 0) {
            deltaRA = -deltaRA;
          }
          distance = sqrt(sqr(deltaDec)+sqr(deltaRA)) * 3600.;
          if (distance < bestdistance) {
            foundMatchFlag = 1;
            bestdistance = distance;
          }
        }
        /* go to the next entry */
        kd_res_next( presults );
            
      }
        
      kd_res_free(presults);
      if (foundMatchFlag != 0) {
        intdistance = bestdistance/(1.0*LOCATION_HISTOGRAM_BINSIZE);
        if ((intdistance < 0) || (intdistance >= (MAX_LOCATION_HISTOGRAM_ARCSEC/LOCATION_HISTOGRAM_BINSIZE))) {
          printf("ERROR: parseidtable index error in line %d\n",__LINE__);
          exit(-1);
        }
        locationHistogram[intdistance]++;
      }
    }
    /* Now dump the histogram */
    cumdistance = 0;
    printf("HHHHdrad\tcount\tcount2\n");
    printf("HHHH----\t-----\t------\n");
    for (intdistance = 0; intdistance < (MAX_LOCATION_HISTOGRAM_ARCSEC/LOCATION_HISTOGRAM_BINSIZE); intdistance++) {
      if (locationHistogram[intdistance] > 0) {
        cumdistance += locationHistogram[intdistance];
        printf("HHHH%d\t\%d\t%d\n",intdistance*LOCATION_HISTOGRAM_BINSIZE,locationHistogram[intdistance],cumdistance);
      }
    }
  }
  


  /* now sum up everything in adjacent bins */
  for (gsc_bin_index = 0; gsc_bin_index < pGscBin->total_gsc_bins; gsc_bin_index++) {
    FindAdjacentBins(pGscBin,gsc_bin_index,bin_list,&bin_count);
    bin_list[bin_count] = gsc_bin_index;
    bin_count++;
    for (bin_index = 0; bin_index < bin_count; bin_index++) {
      gsc_bin_count_multiple[gsc_bin_index] += gsc_bin_count_single[bin_list[bin_index]];
    }
    if (gsc_bin_count_multiple[gsc_bin_index] > gsc_bin_count_maximum) {
      gsc_bin_count_maximum = gsc_bin_count_multiple[gsc_bin_index];
    }
    if (gsc_bin_count_multiple[gsc_bin_index] > 1) {
      gsc_bin_count_total++;
    }
  }
  /* Now set the rejection flag for anything with at least 1 adjacent candidate within a bin distance */
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
#ifdef PREFILTER_MULTDEFECT
    if ((pInput->rejectFlags & PREFILTER_MULTDEFECT_MASK) != 0) {
      continue;
    }
#endif /* PREFILTER_MULTDEFECT */
    if ((pInput->peakRA < 990.0) && (pInput->peakDec < 99.00)) {
      gsc_bin_index = GetGSCBin(pGscBin,pInput->peakRA,pInput->peakDec,&decBin,&raBin,"parseidtable");
    } else {
      gsc_bin_index = GetGSCBin(pGscBin,pInput->ra,pInput->declination,&decBin,&raBin,"parseidtable");
    }
    if (gsc_bin_count_multiple[gsc_bin_index] > 1) {
      pInput->rejectFlags |= REJECT_FLAG_LOCATION;
#if 0
      if ((pInput->peakRA < 990.0) && (pInput->peakDec < 99.00)) {
        printf("REJECT_FLAG_LOCATION ra %f dec %f gsc_bin_index %10d, gsc_bin_count_multiple %d\n",pInput->peakRA,pInput->peakDec,gsc_bin_index,gsc_bin_count_multiple[gsc_bin_index]);
      } else {
        printf("REJECT_FLAG_LOCATION ra %f dec %f gsc_bin_index %10d, gsc_bin_count_multiple %d\n",pInput->ra,pInput->declination,gsc_bin_index,gsc_bin_count_multiple[gsc_bin_index]);
      }
#endif

    }
  }
  time(&curTime);
  curTime -= startTime;
  printf("Searching with new location algorithm at %d seconds\n",curTime);
  /* New location rejection algorithm - reject any pair within MAX_LOCATION_SEPARATION_ARCSEC */
  for (input_index = 0; input_index < input_nrecs; input_index++) {
    pInput = &input_table[input_index];
#ifdef PREFILTER_MULTDEFECT
    if ((pInput->rejectFlags & PREFILTER_MULTDEFECT_MASK) != 0) {
      continue;
    }
#endif /* PREFILTER_MULTDEFECT */
    /* Query the kd tree for a nearby candidate */
    if ((pInput->peakRA < 990.0) && (pInput->peakDec < 99.00)) {
      z = sin(DEGREES_TO_RAD*pInput->peakDec);
      cosx = cos(DEGREES_TO_RAD*pInput->peakDec);
      x = cosx*cos(DEGREES_TO_RAD*pInput->peakRA);
      y = cosx*sin(DEGREES_TO_RAD*pInput->peakRA);
    } else {
      z = sin(DEGREES_TO_RAD*pInput->declination);
      cosx = cos(DEGREES_TO_RAD*pInput->declination);
      x = cosx*cos(DEGREES_TO_RAD*pInput->ra);
      y = cosx*sin(DEGREES_TO_RAD*pInput->ra);
    }
    pt[0] = x;
    pt[1] = y;
    pt[2] = z;
    foundMatchFlag = 0;
    presults = kd_nearest_range(input_ptree, pt,((1.0*MAX_LOCATION_SEPARATION_ARCSEC)*DEGREES_TO_RAD)/3600.0);
    resultSize = kd_res_size(presults);
    while  (!kd_res_end( presults ) ) {
      pInput2 = (PINPUT)kd_res_item(presults,pos);

#ifdef PREFILTER_MULTDEFECT
      if ((pInput2->rejectFlags & PREFILTER_MULTDEFECT_MASK) != 0) {
        /* go to the next entry */
        kd_res_next( presults );
        continue;
      }
#endif /* PREFILTER_MULTDEFECT */

#if 1
      if ((pInput->input_index == 53112) && (pInput2->input_index == 53154)) {
        printf("At indices %d and %d\n",pInput->input_index,pInput2->input_index);
      }
#endif
      if (pInput2->input_index > pInput->input_index)  {
        /* consider only indices we have not searched yet to avoid duplicates */
        factor = cos(DEGREES_TO_RAD* pInput2->declination);
        deltaDec = pInput2->declination - pInput->declination;
        if (deltaDec < 0) {
          deltaDec = - deltaDec;
        }
        deltaRA = pInput2->ra  - pInput->ra;
        if (deltaRA < 0) {
          deltaRA = -deltaRA;
        }
        if (deltaRA > 180.0) {
          deltaRA = 360.0-deltaRA;
        }
        deltaRA  = factor * deltaRA;
        distance = sqrt(sqr(deltaDec)+sqr(deltaRA)) * 3600.;
        if (distance > (1.02*MAX_LOCATION_SEPARATION_ARCSEC)) {
          printf("ERRORD: separation distance is %f instead of %f in index %d and %d line %d\n",
                 distance,
                 (1.02*MAX_LOCATION_SEPARATION_ARCSEC),
                 pInput->input_index,
                 pInput2->input_index,
                 __LINE__);
        }
        pInput->rejectFlags |= REJECT_FLAG_LOCATION;
        pInput2->rejectFlags |= REJECT_FLAG_LOCATION;
      }
      /* go to the next entry */
      kd_res_next( presults );
            
    }
        
    kd_res_free(presults);
  }

  time(&curTime);
  curTime -= startTime;
  printf("With %d candidates, found %d gsc bins out of %d total gsc bins. Maximum bin count is %d in %d seconds\n",
         gsc_single_candidate_count,
         gsc_bin_count_total,
         pGscBin->total_gsc_bins,
         gsc_bin_count_maximum,curTime);

  if (printLocationTable) {
    printf("ra\tdec\tcount\n");
    printf("--\t---\t-----\n");
    for (gsc_bin_index = 0; gsc_bin_index < pGscBin->total_gsc_bins; gsc_bin_index++) {
      if (gsc_bin_count_multiple[gsc_bin_index] > 1) {
        gsc_bin_count_total++;
        GetBinCenter(pGscBin,gsc_bin_index,&tmpra,&tmpdec,"parseidtable2");
        printf("%.4f\t%.4f\t%d\n",tmpra,tmpdec,gsc_bin_count_multiple[gsc_bin_index]);
      }
      
    }
    time(&curTime);
    curTime -= startTime;
    printf("Printed location table in %d seconds\n",curTime);

  }

  if (printMultiplePercentTable) {
    printf("ratio1\n");
    printf("------\n");
    for (input_index = 0; input_index < input_nrecs; input_index++) {
      pInput = &input_table[input_index];
      if ((pInput->rejectFlags & MULTIPLE_MASK) == 0) {
        printf("%f\n",(1.0*pInput->peakMultipleCount)/(1.0*(pInput->peakCountWF+pInput->peakCountNF)));
        ratio_count++;
      }
    }
    time(&curTime);
    curTime -= startTime;
    printf("Have multiple exposure plate ratios for %d candidates in %d seconds\n",ratio_count,curTime);
  }
  memset(rejectionString,0,sizeof(rejectionString));
  memset(secondBestString,0,sizeof(secondBestString));
  memset(bestString,0,sizeof(bestString));
  memset(newString,0,sizeof(newString));
  memset(unusedString,0,sizeof(unusedString));

  for (rejectTableIndex = 0; rejectTableIndex < flareRejectTableSize; rejectTableIndex++) {
    pRejectTable = &flareRejectTable[rejectTableIndex];
    
    /* The newIndex is treated separately because of an overlap */
    if ((rejectTableIndex > newIndex) && 
        (rejectTableIndex <= unusedIndex) &&
        (pRejectTable->flagmask != REJECT_FLAG_ORANGE) &&
        (pRejectTable->flagmask != REJECT_FLAG_YELLOW)) {

      if (newString[0] == 0) {
        strcpy(newString,pRejectTable->flagname);
      } else {
        strcat(newString,", ");
        strcat(newString,pRejectTable->flagname);
      }
      if (strlen(newString) > (sizeof(newString)-40)) {
        printf("ERROR: newString overrun\n");
        exit(-1);
      }
    } 

    if (rejectTableIndex > unusedIndex) {
      if (unusedString[0] == 0) {
        strcpy(unusedString,pRejectTable->flagname);
      } else {
        strcat(unusedString,", ");
        strcat(unusedString,pRejectTable->flagname);
      }
      if (strlen(unusedString) > (sizeof(unusedString)-40)) {
        printf("ERROR: unusedString overrun\n");
        exit(-1);
      }
    } 



    if (rejectTableIndex <= secondBestIndex) {
      if (rejectionString[0] == 0) {
        strcpy(rejectionString,pRejectTable->flagname);
      } else {
        strcat(rejectionString,", ");
        strcat(rejectionString,pRejectTable->flagname);
      }
      if (strlen(rejectionString) > (sizeof(rejectionString)-40)) {
        printf("ERROR: rejectionString overrun\n");
        exit(-1);
      }
    } else if ((rejectTableIndex <= bestIndex) && (rejectTableIndex > secondBestIndex)) {
      if (secondBestString[0] == 0) {
        strcpy(secondBestString,pRejectTable->flagname);
      } else {
        strcat(secondBestString,", ");
        strcat(secondBestString,pRejectTable->flagname);
      }
      if (strlen(secondBestString) > (sizeof(secondBestString)-40)) {
        printf("ERROR: secondBestString overrun\n");
        exit(-1);
      }
    } else if ((rejectTableIndex > bestIndex) && (rejectTableIndex <= unusedIndex)) {
      if (bestString[0] == 0) {
        strcpy(bestString,pRejectTable->flagname);
      } else {
        strcat(bestString,", ");
        strcat(bestString,pRejectTable->flagname);
      }
      if (strlen(bestString) > (sizeof(bestString)-40)) {
        printf("ERROR: bestString overrun\n");
        exit(-1);
      }
    }
  }


  if (rootDirectory[0] != 0) {
    for (input_index = 0; input_index < input_nrecs; input_index++) {
      pInput = &input_table[input_index];
      for (rejectTableIndex = 0; rejectTableIndex < flareRejectTableSize; rejectTableIndex++) {
        pRejectTable = &flareRejectTable[rejectTableIndex];
        if ((pRejectTable->flagmask & pInput->rejectFlags) != 0) {
          rejectTotalCounts[rejectTableIndex]++;
        }
      }
      for (rejectTableIndex = 0; rejectTableIndex < flareRejectTableSize; rejectTableIndex++) {
        pRejectTable = &flareRejectTable[rejectTableIndex];
        if ((pInput->rejectFlags & pRejectTable->flagmask) != 0) {
          if (rejectTableIndex <= unusedIndex) {
            rejectInitialCounts[rejectTableIndex]++;
          }
          if (printMoveCommands != 0) {
#if 0  /* No longer split into multiple directories */
            printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,pRejectTable->subdir,pInput->REF,matchedTimestamp);
#endif
            /* flag the rejected objects */
            if (rejectTableIndex <= secondBestIndex) {
              pInput->rejectedDirectoryFlag = 1;
            }

            /* Now consider second best candidates */
            if ((rejectTableIndex > secondBestIndex) &&
                (rejectTableIndex <= bestIndex)) {
              /* second best.  promote it if it has lowdrad or matches an SDSS object */
              promoteFlag = 0;
              if ((gotGscBinIndexColumn != 0) && (pInput->peakDradRMS3 < MAX_TC_PEAKDRADRMS3)) {
                promoteFlag = 1;
                if ((pInput->rejectFlags & REJECT_FLAG_SDSS) != 0) {
                  promoteBothCount++;
                } else {
                  promoteLowdradCount++;
                }
              } else  if ((pInput->rejectFlags & REJECT_FLAG_SDSS) != 0) {
                promoteFlag = 1;
                promoteSDSSCount++;
              }
              if (promoteFlag != 0) {
                printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"best",pInput->REF,matchedTimestamp);
                pInput->bestDirectoryFlag = 1;
                pInput->subdirectoryFlag = 1;
                bestCount++;
              } else {
                printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"secondbest",pInput->REF,matchedTimestamp);
                pInput->subdirectoryFlag = 1;
                pInput->secondBestDirectoryFlag = 1;
                secondBestCount++;
              }
              if ((pInput->rejectFlags & REJECT_FLAG_MOREPOINTS)) {
                printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"morepoints",pInput->REF,matchedTimestamp);
                morepointsCount++;
              }         
            }
            /* Now do the best candidates */
            if ((rejectTableIndex > bestIndex) &&
                (rejectTableIndex <= unusedIndex)) {
              printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"best",pInput->REF,matchedTimestamp);
              pInput->subdirectoryFlag = 1;
              pInput->bestDirectoryFlag = 1;
              bestCount++;
              if ((pInput->rejectFlags & REJECT_FLAG_MOREPOINTS)) {
                printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"morepoints",pInput->REF,matchedTimestamp);
                morepointsCount++;
              }         
            }
          }
          break;

        }
      }
    
      if (rejectTableIndex > unusedIndex) {
        /* Here no conditions are met */
        rejectTotalCounts[REJECT_TABLE_SIZE]++;
        rejectInitialCounts[REJECT_TABLE_SIZE]++;   
        if (printMoveCommands != 0) {
          printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"best",pInput->REF,matchedTimestamp);
          pInput->subdirectoryFlag = 1;
          pInput->bestDirectoryFlag = 1;
          bestCount++;
        }
      }
      /* Handle the "interest" and "new" tables. */
      if (printMoveCommands != 0) {
        if ((pInput->rejectFlags & REJECT_FLAG_NOVA)) {
          printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"nova",pInput->REF,matchedTimestamp);
          pInput->subdirectoryFlag = 1;
          novaCount++;
        }         
        if ((pInput->rejectFlags & REJECT_FLAG_UG)) {
          printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"ugem",pInput->REF,matchedTimestamp);
          pInput->subdirectoryFlag = 1;
          ugemCount++;
        }         
        if ((pInput->rejectFlags & REJECT_FLAG_VAR)) {
          printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"variable",pInput->REF,matchedTimestamp);
          pInput->subdirectoryFlag = 1;
          variableCount++;
        }         
        


        if ((pInput->rejectFlags & REJECT_FLAG_ORANGE) || (pInput->rejectFlags & REJECT_FLAG_YELLOW)) {
          printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"interest",pInput->REF,matchedTimestamp);
          pInput->subdirectoryFlag = 1;
          interestCount++;
        }         
        if ((pInput->previousFlag == 0) && (rejectTableIndex > newIndex)) {
          /* Handle the NEW list, which by defintion is not in the interest list */
          pInput->newDirectoryFlag = 1;
          pInput->subdirectoryFlag = 1;
          printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"new",pInput->REF,matchedTimestamp);
          newCount++;
        }

      }
      if ((pInput->rejectFlags & REJECT_FLAG_LONGFLARE) != 0) {
        if (pInput->subdirectoryFlag == 1) {
          pInput->bestLongDirectoryFlag = 1;
          printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"bestlongflare",pInput->REF,matchedTimestamp);
          bestLongFlareCount++;
        } else if (((pInput->rejectFlags & REJECT_FLAG_DEFECT) != 0) ||
                   ((pInput->rejectFlags & REJECT_FLAG_MULTIPLE90) != 0) ||
                   ((pInput->rejectFlags & REJECT_FLAG_MULTIPLE33) != 0) ||
                   ((pInput->rejectFlags & REJECT_FLAG_LOCATION) != 0) ||
                   ((pInput->rejectFlags & REJECT_FLAG_SAMETIME) != 0) ||
                   ((pInput->rejectFlags & REJECT_FLAG_MIRA) != 0) ||
                   ((pInput->rejectFlags & REJECT_FLAG_LIM_POINTS) != 0) ||
                   ((pInput->rejectFlags & REJECT_FLAG_LIM_YEARS) != 0) ||
                   ((pInput->rejectFlags & REJECT_FLAG_NOTTC) != 0)) {
          pInput->subdirectoryFlag = 1;
          pInput->badLongDirectoryFlag = 1;
          printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"badlongflare",pInput->REF,matchedTimestamp);
          badLongFlareCount++;
        }
      }
      if ((pInput->rejectFlags & REJECT_FLAG_LONGFLARE) && (pInput->subdirectoryFlag == 0)) {
        pInput->subdirectoryFlag = 1;
        pInput->otherLongDirectoryFlag = 1;
        printf("cp  %s/%s/%s_%s.png %s/%s/%s/%s_%s.png\n",rootDirectory,runDate,pInput->REF,matchedTimestamp,rootDirectory,runDate,"otherlongflare",pInput->REF,matchedTimestamp);
        otherLongFlareCount++;
      }
    }
    time(&curTime);
    curTime -= startTime;
#if 0
    printf("Total Count %5d for no rejection flags set in %d candidates in %d seconds\n",rejectTotalCounts[REJECT_TABLE_SIZE],input_nrecs,curTime);
#else
    printf("Total Count %5d Selection Count %5d for no rejection flags set in %d candidates in %d seconds\n",rejectTotalCounts[REJECT_TABLE_SIZE],rejectInitialCounts[REJECT_TABLE_SIZE],input_nrecs,curTime);
#endif
    for (rejectTableIndex = (flareRejectTableSize-1); rejectTableIndex >= 0; rejectTableIndex--) {
      pRejectTable = &flareRejectTable[rejectTableIndex];
#if 0
      printf("Total Count %5d for condition '%-10s' (%s)\n",rejectTotalCounts[rejectTableIndex],pRejectTable->flagname,pRejectTable->webname);
#else
      printf("Total Count %5d Selection Count %5d for condition '%-10s' in directory '%-14s' (%s)\n",rejectTotalCounts[rejectTableIndex],rejectInitialCounts[rejectTableIndex],pRejectTable->flagname,pRejectTable->subdir,pRejectTable->webname);
#endif
    }
    for (rejectTableIndex = (flareRejectTableSize-1); rejectTableIndex >= 0; rejectTableIndex--) {
      pRejectTable = &flareRejectTable[rejectTableIndex];
      printf("Condition '%-10s' (%s)\n",pRejectTable->flagname,pRejectTable->webname);
    }
  }
  if (printMoveCommands != 0) {
    printf("lowdrad candidates have flare RMS of less than %f arcsec\n",MAX_TC_PEAKDRADRMS3);
    printf("Start by rejecting all candidates with conditions %s\n",rejectionString);
    printf("Second best candidates have the conditions: %s\n",secondBestString);
    printf("Best candidates have no conditions, have lowdrad, or are close to an SDSS object or have the conditions: %s\n",bestString);
    printf("new candidates did not appear in previous runs have no conditions or have only the conditions: %s\n",newString);
    printf("The following conditions are currently not used unless otherwise mentioned: %s\n",unusedString);
    printf("%5d total best candidates in directory %s/%s/best\n",bestCount,rootDirectory,runDate);
    printf("    which includes %5d second best candidates with lowdrad\n",promoteLowdradCount);
    printf("    and            %5d second best candidates close to an SDSS object\n",promoteSDSSCount);
    printf("    and            %5d second best candidates with both lowdrad and near an SDSS object\n");
    printf("%5d second best candidates not included above in directory %s/%s/secondbest\n",secondBestCount,rootDirectory,runDate);
    printf("%5d best and second best candidates with more points than the previous search in directory  %s/%s/morepoints\n",morepointsCount,rootDirectory,runDate);
    printf("%5d orange and yellow candidates (ignoring any reject conditions) in directory %s/%s/interest\n",interestCount,rootDirectory,runDate);

    printf("%5d nova candidates (ignoring any reject conditions) in directory %s/%s/nova\n",novaCount,rootDirectory,runDate);
    printf("%5d ugem candidates (ignoring any reject conditions) in directory %s/%s/ugem\n",ugemCount,rootDirectory,runDate);
    printf("%5d other variable candidates (ignoring any reject conditions) in directory %s/%s/variable\n",variableCount,rootDirectory,runDate);
    printf("%5d candidates new to this run in directory %s/%s/new\n",newCount,rootDirectory,runDate);
    printf("%5d long flare candidates in above directories are also in directory %s/%s/bestlongflare\n",bestLongFlareCount,rootDirectory,runDate);
    printf("%5d long flare candidates with defect,multiple90,multiple33,location,sametime,mira,nottc,coveragept,coverageyr reject conditions in directory %s/%s/badlongflare\n",badLongFlareCount,rootDirectory,runDate);
    printf("%5d other long flare candidates with other reject conditions in directory %s/%s/otherlongflare\n",otherLongFlareCount,rootDirectory,runDate);


  }
  if ((candidate_name[0] != 0) && (possible_handle != NULL)) {
    fprintf(possible_handle,"RA_TC\tDec_TC\tpeakDradRMS3\tpeakLimitingYears\tNpts_TC\tNpts_ext\tpeakYear\tpeakMaxMag\tnearbyObjectDistance\tpeakDays\tname\tnearbyVSX\tflag\tcomments\n");
    fprintf(possible_handle,"-----\t------\t------------\t-----------------\t-------\t--------\t--------\t----------\t--------------------\t--------\t----\t---------\t----\t--------\n");
    /* Now for every thumbnail with coordinates, find the closest flare candidate */
    for (thumbnail_index = 0; thumbnail_index < thumbnail_count; thumbnail_index++) {
      pThumbnail = &thumbnail_table[thumbnail_index];
      if (pThumbnail->haveCoords == 0) {
        /* Search the input table for coordinates */
        for (input_index = 0; input_index < input_nrecs; input_index++) {
          pInput = &input_table[input_index];
          if (strcmp(pInput->REF,pThumbnail->name) == 0) {
            pThumbnail->haveCoords = 1;
            pThumbnail->ra = pInput->ra;
            pThumbnail->dec = pInput->declination;
            break;
          }
        }
#if 0
        if (input_index == input_nrecs) {
          printf("ERROR: Unrecognized thumbnail index %6d %s\n",thumbnail_index,pThumbnail->filename);
          exit(-1);
        }
#endif
      }
      if (pThumbnail->haveCoords == 1) {
#if 1 /* candidate_ptree */
        /* Query the kd tree for a nearby candidate */
        z = sin(DEGREES_TO_RAD*pThumbnail->dec);
        cosx = cos(DEGREES_TO_RAD*pThumbnail->dec);
        x = cosx*cos(DEGREES_TO_RAD*pThumbnail->ra);
        y = cosx*sin(DEGREES_TO_RAD*pThumbnail->ra);
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
        foundMatchFlag = 0;
        presults = kd_nearest_range(candidate_ptree, pt,(MAX_DISTANCE*DEGREES_TO_RAD)/3600.0);
        resultSize = kd_res_size(presults);
        while  (!kd_res_end( presults ) ) {
          pCandidate2 = (PFLARECANDIDATE)kd_res_item(presults,pos);
          factor = cos(DEGREES_TO_RAD* pCandidate2->dec);
          deltaDec = pCandidate2->dec - pThumbnail->dec;
          if (deltaDec < 0) {
            deltaDec = - deltaDec;
          }
          deltaRA  = factor *(pCandidate2->ra  - pThumbnail->ra);
          if (deltaRA < 0) {
            deltaRA = -deltaRA;
          }
          distance = sqrt(sqr(deltaDec)+sqr(deltaRA)) * 3600.;
          if (distance < MAX_DISTANCE) {
            if (pThumbnail->candidate_index < 0) {
              pThumbnail->candidate_index = candidate_index;
              pThumbnail->candidateDistance = distance;
#if 0
              if (strstr(pCandidate->comment,"Possible strong cand. Not on DSS so B >20") != NULL) {
                printf("line %d candidate_index %d comment %s\n",__LINE__,candidate_index,pCandidate->comment);
              }
#endif

            } else {
              if (pThumbnail->candidateDistance > distance) {
                pThumbnail->candidate_index = candidate_index;
                pThumbnail->candidateDistance = distance;
              }
            }
          } 
          /* go to the next entry */
          kd_res_next( presults );
            
        }
        
        kd_res_free(presults);
   

#else /* candidate_ptree */
        for (candidate_index = 0; candidate_index < candidate_count; candidate_index++) {
          pCandidate = &candidate_table[candidate_index];
          if (pCandidate->unrecognized != 0) {
            continue;
          }
          distance = 3600.*wcsdist(pCandidate->ra,pCandidate->dec,pThumbnail->ra,pThumbnail->dec);
          if (distance < MAX_DISTANCE) {
            thumbnailCandidateNearbyCount++;
            if (pThumbnail->candidate_index < 0) {
              pThumbnail->candidate_index = candidate_index;
              pThumbnail->candidateDistance = distance;
#if 0
              if (strstr(pCandidate->comment,"Possible strong cand. Not on DSS so B >20") != NULL) {
                printf("line %d candidate_index %d comment %s\n",__LINE__,candidate_index,pCandidate->comment);
              }
#endif

            } else {
              if (pThumbnail->candidateDistance > distance) {
                pThumbnail->candidate_index = candidate_index;
                pThumbnail->candidateDistance = distance;
                thumbnailCandidateBetterCount++;
              }
            }
          }
        }
#endif /* candidate_ptree */
      } else {
        thumbnailNoCoordsCount++;
      }
    }
    time(&curTime);
    curTime -= startTime;
    printf("At line %d in %d seconds\n",__LINE__,curTime);
    /* Now try to locate a thumbnail for every input record */
    
    
    for (candidate_index = 0; candidate_index < candidate_count; candidate_index++) {
      pCandidate = &candidate_table[candidate_index];
      possibleCount++;
      if (pCandidate->input_index < 0) {
        if (printWarnings) {
          printf("WARNING: no record found for candidate %s %s\n",pCandidate->name,pCandidate->comment);
        }
        possibleMissingCount++;
      } else {
        if (pCandidate->input_index >= input_nrecs) {
          printf("ERROR: Illegal input_index for candidate %s %s\n",pCandidate->name,pCandidate->comment);
          exit(-1);
        }
        pBestThumbnail = NULL;
        for (thumbnail_index = 0; thumbnail_index < thumbnail_count; thumbnail_index++) {
          pThumbnail = &thumbnail_table[thumbnail_index];
          if (pThumbnail->candidate_index  == candidate_index) {
            if (strcmp(pThumbnail->subdir,pCandidate->runTimestamp) == 0) {
              /* Perfect match with the timestamp.  We are done */
              pBestThumbnail = pThumbnail;
              break;
            }
#if 0
            if (pBestThumbnail == NULL) {
              pBestThumbnail = pThumbnail;
            } else {
              if (strcmp(pThumbnail->subdir,pBestThumbnail->subdir) > 0) {
                pBestThumbnail = pThumbnail;
              }
            }
#endif
          }
        }
      }
      
    }  


    time(&curTime);
    curTime -= startTime;
    printf("Total possible flares %d, flares found %d flares with missing records %d in %d seconds\n",possibleCount,possibleCount-possibleMissingCount,possibleMissingCount,curTime);
  }

  printf("Thumbnails %d  thumbnailNoCoordsCount %d, thumbnailCandidateNearbyCount %d, thumbnailCandidateBetterCount %d\n",
         thumbnail_count,
         thumbnailNoCoordsCount,
         thumbnailCandidateNearbyCount,
         thumbnailCandidateBetterCount);

  if ((candidate_name[0] != 0) && (possible_handle != NULL)) {
    for (thumbnail_index = 0; thumbnail_index < thumbnail_count; thumbnail_index++) {
      pThumbnail = &thumbnail_table[thumbnail_index];
#ifdef RUNDATE_ONLY
      if (strstr(pThumbnail->filename,runDate) == NULL) {
        continue; /* Not from our run */
      }
#endif /* RUNDATE_ONLY */
      pCandidate = NULL;
      pInput = NULL;
      for (input_index = 0; input_index < input_nrecs; input_index++) {
        pInput = &input_table[input_index];
        if (strstr(pThumbnail->filename,pInput->REF) != NULL) {
          if ((pInput->minDistance < MAX_DISTANCE) &&
              (pInput->candidate_index >= 0) &&
              (pInput->candidate_index < candidate_count)) {
            pCandidate =  &candidate_table[pInput->candidate_index];
#if 0
            if (strstr(pCandidate->comment,"Possible strong cand. Not on DSS so B >20") != NULL) {
              printf("line %d candidate_index %d comment %s\n",__LINE__,pInput->candidate_index,pCandidate->comment);
            }
#endif
          }
          break;
          
        }
      }
      if (input_index == input_nrecs) {
        pInput = NULL;
      }
      if (pCandidate == NULL) {
        pCandidate = &dummy_candidate;
        memset(pCandidate,0,sizeof(FLARECANDIDATE));
        strcpy(pCandidate->name,pThumbnail->name);
#if 1
        pCandidate->comment[0] = 0;
#else
        strcpy(pCandidate->comment,"NOT ON THE PREVIOUS FLARECANDIDATE LIST");
#endif
        if (pThumbnail->haveCoords) {
          pCandidate->ra = pThumbnail->ra;
          pCandidate->dec = pThumbnail->dec;
        } else {
          pCandidate->ra = 999.;
          pCandidate->dec = 99.;
        }
      }
      PrintPossibleTable(possible_handle,all_handle,runDate,pThumbnail,pCandidate,pInput,gotGscBinIndexColumn,verbose,&subdirectoryCount);
    } 
  }
  /* Now look for previous candidates not represented in the current run */
  for (candidate_index = 0; candidate_index < candidate_count; candidate_index++) {
    pCandidate = &candidate_table[candidate_index];
    if (pCandidate->allNameOutputFlag != 0) {
      continue;
    }
    PrintPossibleTable(NULL,all_handle,runDate,pThumbnail,pCandidate,NULL,gotGscBinIndexColumn,verbose,&subdirectoryCount);
  }
  printf("Total candidates %d are in %s\n",subdirectoryCount,possible_name);

  time(&curTime);
  curTime -= startTime;
  printf("Completed in %d seconds\n",curTime);
  printf("sametimeCount %d\n",sametimeCount);
  if (candidate_handle != NULL) {
    fclose(candidate_handle);
    candidate_handle = NULL;
  }
  if (interesting_handle != NULL) {
    fclose(interesting_handle);
    interesting_handle = NULL;
  }

  if (thumbnail_handle != NULL) {
    fclose(thumbnail_handle);
  }

  if (possible_handle != NULL) {
    fclose(possible_handle);
  }

  if (out_handle != NULL) {
    fclose(out_handle);
  }
  if (all_handle != NULL) {
    fclose(all_handle);
  }
  if (inHandle != NULL) {
    fclose(inHandle);
  }
  if (candidate_table != NULL) {
    free(candidate_table);
  }
  if (unrecognized_table != NULL) {
    free(unrecognized_table);
  }
  if (thumbnail_table != NULL) {
    free(thumbnail_table);
  }
  Close(id_handle);
  free(input_table);
  if (gsc_bin_count_single != NULL) {
    free(gsc_bin_count_single);
  }
  if (gsc_bin_count_multiple != NULL) {
    free(gsc_bin_count_multiple);
  }
  kd_free(candidate_ptree);
  kd_free(input_ptree);
}

