// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* filter_scratch.c
 *
 *  From correspondance of Mon 3/02/15 10:47 AM:
 *
 * Borrowing an idea from Sumin Tang, if there are N objects evenly distributed in an area A of the plate, then the mean 
 * separation R between these objects is approximately  PI * (R**2) * N ~ A, or  R ~ sqrt(A/(PI*N)).  In the attached figures, 
 * I set R = K*sqrt(A/(PI*N)) where experimentation shows that K = 0.3 appears to be the minimum necessary to form reasonably 
 * long chains from the scratches.  To find N without including lots of grain noise while recognizing the presence of multiple 
 * exposures, I selected a cutoff MAG_ISO such that N = 4 * (number of images matched to the GSC2.3.2 catalog).    I found that 
 * for deeper scratches, SExtractor already recognizes components of the scratch through the ELLIPTICITY, THETA_IMAGE, aLength, 
 * and bLength parameters.  I therefore plotted rectangles using these SExtractor parameters, but expanding aLength and bLength 
 * by R.  Results are given for the six examples listed below where the defect filter did not work.  Two of the images (i49523 
 * and ac17627) were too small for this algorithm to be effective, and might be addressed by enhancements in Sumin's defect algorithm.  
 * The other four images do show good chains that do not require straight line, vertical defects.  Needless to say, there are a lot of 
 * false positives in the images.  Additional experimental parameters will be needed to reduce the number of false positives.  These 
 * additional parameters could include the length of the chain and the percentage of area of the bounding rectangle occupied by the 
 * chain.  The additional parameters, however, could also raise the false negative rate for groups of parallel scratches.
 *
 * In non-production mode, the above algorith is executed only for a given radius around a (ra,dec) point on the plate
 * In production mode, the plate is divided into a 10 x 10 grid of "sections" and the radius is calculated for each of these sections
 * 
 *  Output file:
 *    series
 *    plateNumber
 *    mosaicNumber
 *    memberCount   - number of objects in a candidate scratch chain
 *    aveNeighborCt - ave number of neighbors of objects in a candidate scratch chain
 *    aveLengthRatio    - SExtractor aLength/bLength of objects in a candidate scratch chain
 * 
 *
 * gcc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include   filter_scratch.c pipelineutils.a -L /dasch/install/lib -lm  -L/usr/lib${lib64}/mysql -lmysqlclient -ldl -pthread  -ltable -lutil -lwcs     -o filter_scratch 

 *
 *
 * getlocation am06907_01_01ww -e 0 | awk '{print $4}'  gives 6.3551506
 *  createregion  -i /home/scanner/Pipeline/match/am06907_01_01ww_tnx.db -o /home/scanner/junk/scratch/los.reg -r 181.0353 -d 13.3566 -a 0.71 -l
 *
 *
 * /dasch/Pipeline/filter_scratch  -c 0.30 -v -g /dasch/raid022/ExposureData/Mosaics/i/52751_00/i52751_00_01ww_scratch.reg -f -p i52751_00_01ww -t 1.788 -w 17412 -h 22026 -m /home/scanner/Pipeline/match/match_i52751_00_01ww_tnx_u.db -s /home/scanner/Pipeline/match/i52751_00_01ww_tnx.db -o /home/scanner/Pipeline/match/i52751_00_01ww_scratch.db
 *     | grep "STARBASE"  | awk '{sub(/STARBASE/,"",$0);print $0}' | sort -u  > los.db
 *
 * Feb 28, 2015 Edward J. Los - Experimental initial version
 * May  1, 2015 Edward J. Los - Back out original limiting magnitude code
 * Dec 23, 2016 Edward J. Los - Fix query to the spatialbin table.
 * 
 */


#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "table.h"

#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "kdtree.h"

double round(double x);
/* #define LOS_DEBUG 1 */
/* #define LOS_DEBUG_NUMBER 185702 */
/* #define LOS_DEBUG_CHAIN 69096 */
#define MAX_BUFFER 200
#if 1
#define MAG_ISO_DIFF -0.625 /* setting that produces the best results */
#else
#define MAG_ISO_DIFF 0.00 /* Difference between the maximum (dimmest)  and the doubly clipped MAG_ISO_med */
#endif
/* #define DEBUG_SEXTRACTOR 1 */  /* If set, print all SExtractor regions and exit */



#define MAX_FLUX_LEVELS 20
#define MAX_MAG_LEVELS 30
#define MAX_LEVEL_VECTOR 200
#define NUM_SECTIONS 10  /* Divide the plate into a 10x10 grid for the purpose of getting average distances */
#define TOTAL_SECTIONS (NUM_SECTIONS*NUM_SECTIONS)

#define MINIMUM_CHAIN_COUNT 3 /* Minimum count for a chain to be selected */
#if 1
#define CYAN_CHAIN_COUNT 3 /* Minumum count for highlighted colors */
#else
#define CYAN_CHAIN_COUNT 10 /* Minumum count for highlighted colors */
#endif
#define CLIP_RMS_FACTOR 1.0

typedef struct _section {
  int section_index;
  int section_images;
  double maxaLength;
  double maxbLength;
  double MAG_ISO_radius_pixels;  /* expansion radius for this section in pixels */
  double search_radius_pixels; /* The search radius takes into account adjacent bins */\ 
} SECTION,*PSECTION;

/* This is an array of distances to the neighbors */
typedef struct _neighbor {
  int neighbor_index;  /* Index of this entry */
  int path_flag;       /* Set if this entry is on the selected path */
  double neighborDist; /* Distance to the neighbor in pixels */
  struct _starimage *pSextractor; /* Pointer to the neighbor */
  struct _neighbor *pNeighborFlink; /* Pointer to the next neigbor structure */
} NEIGHBOR,*PNEIGHBOR;


/* This is an array of matched chains, the same size as the starimage table */
typedef struct _chain {
  int selected;
  int chain_index;
  int memberCount; /* Count of objects in the chain */
  int pathCount; /* Count of objects on the shortest chain*/
  int matchCount; /* Count of objects that have catalog matches */
  double aveNeighborCt;
  double aveLengthRatio;
  struct _starimage *clink; /* link to first member of a chain */
  double max_separation; /* distance in pixels between the furthest two members, measured between connecting centers */
  int max_separation_count; /* The number of members along the max_separation route */
  int max_separation_flag; /* Flag the entry as being on the max_separation route (1 = in chain; 2 = chain endpoint*/
} CHAIN,*PCHAIN;

typedef struct _internal {
  int sextractor_index;
  int flag_match;      /* Object appears in the match table if nonzero */
  int flag_output;     /* Object candidate for the output plot */
  int flag_scratch;      /* Object believed to be a scratch */
  int flag_chain_length;  /* chain exceeds CYAN_CHAIN_COUNT */
  int flag_length_ratio; /* chain exceeds pct90LengthRatio */
  int flag_MAG_ISO;    /* Object makes the cut for MAG_ISO */
  int section_index; /* Plate section index = xsection + (NUM_SECTIONS + ysection) */  
  int neighbor_count; /* Number of matching neighbors */
  int spatial_bin;  
  struct _starimage *flink; /* forward link in a chain of matched images */
  int chain_id;  /* Chain identification NON-zero */
  int brightest_member; /* First member in the chain to be flagged (print text string) */
  PNEIGHBOR pNeighbor; /* Pointer to the neighbor chain */ 
  int path_flag;       /* Set if on the longest path */
  double testDist;     /* Distance from the selected image to this image - initial search */
  double testDist2;    /* Distance from the selected image to this image - path search */
} INTERNAL,*PINTERNAL;

/* Output of the limiting_iso column of the spatialbin table for the GSC2.3.2 calibration on April 6, 2015*/
typedef struct _testplate {
  char series[MAX_SERIES_STRING];
  int plateNumber;
  double mean_MAG_ISO;
  double median_MAG_ISO;
  double rms_MAG_ISO;
  double bin_MAG_ISO[MAX_SPATIAL_BINS+1];
} TESTPLATE, *PTESTPLATE;

TESTPLATE testPlateTable[] = {
  /* plate     mean   median rms           1      2      3      4      5      6      7      8      9  */
  {"ac",13076, -9.33, -9.29,0.30,  0.00, -9.12, -8.99, -9.12, -9.19, -9.28, -9.29, -9.41, -9.62, -9.96},
  {"ac",17076, -9.28, -9.28,0.21,  0.00, -9.00, -9.07, -9.12, -9.22, -9.28, -9.28, -9.38, -9.48, -9.68},
  {"ac",17627, -9.70, -9.72,0.13,  0.00, -9.64, -9.67, -9.69, -9.73, -9.93, -9.72, -9.52, -9.55, -9.86},
  {"ac",20795, -9.48, -9.51,0.14,  0.00, -9.30, -9.40, -9.36, -9.42, -9.52, -9.51, -9.50, -9.57, -9.76},
  {"ac",34690, -9.51, -9.53,0.35,  0.00, -9.68, -9.26, -9.22, -9.27, -9.21, -9.33, -9.53, -9.91,-10.18},
  {"ac",37705, -9.58, -9.58,0.26,  0.00, -9.79, -9.40, -9.34, -9.46, -9.38, -9.45, -9.58, -9.71,-10.14},
  {"ac",41119,-10.30,-10.35,0.32,  0.00,-10.18,-10.00, -9.97,-10.01,-10.14,-10.35,-10.49,-10.70,-10.86},
  {"ac",42227,-10.27,-10.31,0.44,  0.00, -9.88, -9.88, -9.93,-10.01,-10.08,-10.31,-10.52,-10.73,-11.14},
  {"am",26017,-10.05, -9.94,0.70,  0.00, -9.44, -9.56, -9.57, -9.64, -9.82, -9.94,-10.11,-10.81,-11.55},
  {"am", 6907, -8.94, -8.94,0.11,  0.00, -8.83, -8.86, -8.82, -8.94, -8.93, -8.97, -8.94, -9.20, -9.00},
  {"am", 9230, -9.08, -9.09,0.16,  0.00, -9.03, -9.08, -9.13, -9.09, -9.00, -8.85, -8.89, -9.30, -9.33},
  {"b" ,16412, -9.29, -9.32,0.26,  0.00, -8.92, -9.08, -9.26, -9.64, -9.47, -9.38,  0.00,  0.00,  0.00},
  {"b" ,46256, -9.22, -9.21,0.50,  0.00, -8.62, -9.28, -9.14, -9.83,  0.00,  0.00,  0.00,  0.00,  0.00},
  {"dnb",3393,-11.24,-11.35,0.56,  0.00,-10.45,-10.64,-11.06,-11.35,-11.68,-11.83,-12.12,-10.89,-11.19},
  {"i" ,31013, -9.84, -9.86,0.24,  0.00, -9.49, -9.64, -9.80, -9.86, -9.79,-10.15,-10.12,  0.00,  0.00},
  {"i" ,49523,-10.66,-11.06,0.67,  0.00, -9.88,-10.11,-10.38,-10.24,-11.06,-11.55,-11.40,  0.00,  0.00},
  {"i" ,52750,-10.42,-10.63,0.38,  0.00,-10.03,-10.03,-10.30,-10.19,-10.99,-10.80,-10.63,  0.00,  0.00},
  {"i" ,52751,-10.42,-10.50,0.37,  0.00,-10.05,-10.01,-10.22,-10.36,-10.50,-10.99,-10.80,  0.00,  0.00},
  {"mc",35417,-10.37,-10.37,0.14,  0.00,-10.17,-10.21,-10.37,-10.30,-10.34,-10.41,-10.34,-10.59,-10.57},
  {"mc",39048, -9.77, -9.78,0.19,  0.00, -9.83, -9.98,-10.11, -9.43, -9.76, -9.73, -9.72, -9.63, -9.78},
#if 0
  {"mc", 8621,  0.0 ,  0.0 ,0.0 ,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00},    /* No solution */
#endif
  {"mf",40561,-10.37,-10.65,0.36,  0.00,-10.65,-10.65,-10.65,-10.65, -9.93, -9.96, -9.89,-10.67,-10.29},
  {"rb",12499,-10.09,-10.10,0.64,  0.00, -9.52, -9.60, -9.63, -9.83, -9.86,-10.10,-10.17,-10.60,-11.54},
};

int testPlateTableSize = sizeof(testPlateTable)/sizeof(TESTPLATE);

/* Actual samples of scratches and of normal images.  Used for testing only! */
typedef struct _sample {
        char series[MAX_SERIES_STRING]; /* Plate Series */
        int plateNumber;                /* Plate Number */
        int mosaicNumber;               /* Mosaic Number */
        int NUMBER;                     /* SExtractor source number */
        int scratchFlag;                /* 1 if considered a scratch or important defect */
        char notes[MAX_NOTES_STRING];   /* Notes */
} SAMPLE,*PSAMPLE;


SAMPLE sampleTable[] = {
  {"am",9230 ,0,593818,1,"Meeting of August 7, 2014: S23101322531 (not on dell)"},
  {"i" ,49523,1,133103,1,"Mon 2/16/15 6:41 PM N12212025274 spreadsheet above noise level (not on dell)"},
  {"ac",17627,0,410385,1,"Mon 2/16/15 6:41 PM N133102013709,spreadsheet above noise level"},
  {"i" ,52751,0, 47323,0,"m:31 chain, individual stars"},
  {"i" ,52751,0, 66954,0,"m:11 chain, individual stars"},
  {"i" ,52751,0, 39475,1,"m:7 scratch"},
  {"i" ,52751,0, 46668,1,"m:36 parallel scratch"},
  {"i" ,52751,0, 31238,1,"m:46 parallel scratch"},
  {"i" ,52751,0,121885,1,"m:3 irregular scratch"},
  {"i" ,52751,0,229335,1,"m:9 straight line scratch"},
  {"mc", 8621,0,185702,1,"m:3 plate number"},
  {"b" ,46256,0,823155,1,"m:4 irregular scratch"},
  {"ac",17076,1, 83318,0,"m:3 Edge of calibration square"}, /* changed to no scratch on Apr 28, 2015 */
  {"i", 52750,0, 49165,0,"m:32 individual stars, scratch in background"},
  {"i", 52750,0, 37239,0,"m:3 single defect includes additional stars"},
  {"i", 52750,0,642573,0,"m:13 spot defects appear as scratch"},
  {"i", 52750,0, 28293,1,"m:34 broad scrape"},
  {"i", 31013,0,173825,1,"m:89 scratch complex"},
  {"i", 31013,0,201797,1,"m:88 scratch complex"},
  {"i", 31013,0,281970,1,"m:19 low brightness scuff"},
  {"i", 31013,0,765425,1,"m:27 straight line scratch"},
  {"i", 31013,0,696404,1,"m:127 parallel scratch region"},
  {"i", 31013,0,604894,0,"m:39 not scratch, may be defect"},
  {"i", 31013,0,542347,0,"m:17 not scratch"},
  {"ac",34690,1,380921,1,"MAG_ISO -12.8 on scratch"},
  {"ac",34690,1,383248,1,"MAG_ISO -12.1 short image, one SE hit"}, /* changed to scratch on Apr 28, 2015 */
  {"ac",34690,1,128999,1,"MAG_ISO -14.6 on scratch"},
  {"ac",34690,1,150417,1,"MAG_ISO -14.0 on scratch"},
  {"ac",34690,1,147116,1,"MAG_ISO -13.5 on scratch"},
  {"ac",34690,1, 89076,1,"MAG_ISO -16.4 on scratch"},
  {"am",26017,0,195157,0,"MAG_ISO -14.5 chain of stars"},
  {"am",26017,0,288151,1,"MAG_ISO -14.1 on scratch"},
  {"am",26017,0,129854,1,"MAG_ISO -15,6 on scratch"},
  {"b" ,16412,1,1141745,1,"M:14 straight line scratch"},
  {"b" ,16412,1, 985032,0,"m:27 individual stars"},
  {"b" ,16412,1,363130,1,"m:12 low brightness scratch"},
  {"b" ,16412,1,713091,1,"M:19 low brightness scratch"},
  {"i" ,52750,0,372494,1,"m:18 curved scratch"},
  {"i" ,52750,0,373049,1,"m:13 ring of joined"},
  {"i" ,31013,0,777840,1,"m:67 long curvy scratch"},
  {"i" ,31013,0, 46874,1,"m:17 irregular scratch"},
  {"i" ,31013,0,328240,0,"m:13 stars and galaxies?"},
  {"i" ,31013,0,454474,1,"m:15 parallel scratches"},
  {"i" ,31013,0,658256,0,"m:12 stars or defects"},
  {"i" ,31013,0,797608,0,"m:12 stars"},
  {"mc",39048,1,277073,1,"m:99 scratch"},
  {"mc",39048,1,566902,0,"m:12 individual stars"},
  {"mc",39048,1,445195,0,"m:10 individual stars near limiting mag?"},
  {"mc",39048,1,171680,1,"m:22 broad scuff?"},
  {"mc",39048,1,196512,1,"m:35 broad scuff?"},
  {"mc",39048,1,390335,0,"m:10 individual stars"},
  {"mc",39048,1,518271,0,"m:17 individual stars"},
  {"mc",39048,1,562769,0,"m:10 individual stars"},
  {"mc",39048,1,677740,1,"m:21 scrape"},
  {"mc",39048,1,658113,0,"m:15 individual stars"},
  {"mc",39048,1,422627,0,"m:11 individual stars"},
  {"am", 6907,1,431745,1,"Mon 2/16/15 6:41 PM N12021311447,spreadsheet"},
  {"am", 6907,1, 46545,1,"m:92 parallel scratches"},
  {"am", 6907,1,478369,0,"m:13 individual stars or defects?"},
  {"am", 6907,1,623592,0,"m:11 individual stars"},
  {"am", 6907,1,649072,1,"m:13 scratch"},
  {"am", 6907,1,555780,0,"m:10 individual stars or defects?"},
  {"am", 6907,1,352389,0,"m:4 individual stars"},
  {"am", 6907,1,164527,0,"m:15 individual stars or low level noise?"},
  {"am", 6907,1,187623,0,"m:40 background smudge?"},
  {"ac",20795,0,325294,1,"m:12 scratch"},
  {"ac",20795,0,146768,0,"m:10 individual stars or large galaxy?"},
  {"ac",20795,0,124998,1,"m:64 parallel scratches"},
  {"ac",20795,0,275749,0,"m:13 individual stars, galaxy?"},
  {"ac",20795,0,647183,0,"m:3 individual stars"},
  {"ac",20795,0,690819,0,"m:10 individual stars"},
  {"ac",20795,0,639622,0,"m:10 individual stars"},
  {"ac",20795,0,534397,1,"m:24 scratch"},
  {"ac",20795,0,537150,1,"m:38 scrape area"},
  {"ac",13076,1,174458,1,"Original example (c) of 12/06/13 11:50 AM N1312210861 "},
  {"ac",13076,1,157933,1,"Example b, single deep scratch"},
  {"ac",13076,1,104287,1,"Example d, perpendicular scratch"},
  {"ac",13076,1,103855,1,"Example d, single scratch"},
  {"ac",13076,1,240820,1,"m:59 Example h, broad streak"},
  {"ac",13076,1,240208,1,"m:13 Example h, broad streak"},
  {"ac",13076,1,148128,1,"m:23 Example e, curved scrape"},
  {"ac",13076,1,366679,1,"m:59 Example f, broad regional area"},
  {"ac",13076,1,375852,1,"m:70 Example f, broad regional area"},
  {"ac",13076,1,479615,1,"m:32 Example g, leftover ink mark?"},
  {"ac",13076,1,483480,1,"m:29 Example g, leftover ink mark?"},
  {"ac",13076,1,511190,1,"m:4  Example i, handwriting"},
  {"ac",13076,1,511364,1,"m:5  Example i, handwriting"},
  {"ac",13076,1,167570,1,"m:32 Example c, Less deep scratch, well connected"},
  {"ac",13076,1,221523,0,"m:11 individual stars"},
  {"ac",13076,1,172643,0,"m:4  individual stars"},
  {"ac",13076,1,460477,0,"m:11 individual stars"},
  {"ac",13076,1,458304,0,"m:13 individual stars"},
  {"mf",40561,1, 71543,1,"m:175 halation ring"},
  {"mf",40561,1, 68602,0,"m:15 dark area"},
  {"mf",40561,1,102205,0,"m:3 three stars"},
  {"mf",40561,1, 98386,0,"m:11 individual stars"},
  {"mf",40561,1,117714,0,"m:26 dark area"},
  {"mf",40561,1,143055,0,"m:4 individual stars, one large"},
  {"mf",40561,1,288463,0,"m:12 individual stars"},
  {"mf",40561,1,  5038,1,"m:27 scratch/defect"},
  {"mf",40561,1,  4082,1,"m:16 defect"},
  {"mf",40561,1,354967,1,"m:28 scratch"},
  {"mc",35417,1,257773,0,"m:25 multiple exposure"},
  {"mc",35417,1,415584,0,"m:85 multiple exposure"},
  {"mc",35417,1,476477,1,"MAG_ISO -15.3 on defect"},
  {"mc",35417,1,133622,1,"MAG_ISO -12.7 on hair"},
  {"ac",37705,1, 69236,1,"Mon 2/16/15 6:41 PM N12332305652,spreadsheet above noise level"},
  {"ac",37705,1,236128,0,"m:17 individual stars"},
  {"ac",37705,1,256575,0,"m:39 stars around two bright"},
  {"ac",37705,1,  2855,1,"m:53 scratch"},
  {"ac",37705,1,  4590,1,"m:5 bright star next to scratch"},
  {"ac",37705,1, 10745,1,"m:136 parallel scratch region"},
  {"ac",37705,1, 17087,1,"m:48 parallel scratch region"},
  {"rb",12499,0,306853,1,"m:156: background smudge?"},
  {"rb",12499,0,343764,0,"m:11: individual stars"},
  {"rb",12499,0,364441,1,"m:27: dual parallel scratches"},
  {"rb",12499,0,413332,0,"m:15: individual stars"},
  {"rb",12499,0,431778,1,"m:62: dual parallel scratches"},
  {"rb",12499,0,415370,0,"m:26: stars/smudge near limiting mag"},
  {"rb",12499,0,310718,1,"m:39: defect"},
  {"ac",41119,0, 45936,1,"Mon 2/16/15 6:41 PM N23231135063,spreadsheet kink "},
  {"ac",41119,0,319104,1,"m:60: scratch region"},
  {"ac",41119,0,284221,0,"m:4: stars at limiting mag"},
  {"ac",41119,0,271991,0,"m:13: stars near limiting mag"},
  {"ac",41119,0, 84056,1,"m:17: scratch complex"},
  {"ac",41119,0,102533,0,"m:13: stars near limiting mag"},
  {"ac",41119,0, 57972,0,"m:16: stars near limiting mag"},
  {"ac",42227,0,189134,1,"m:19: halation ring"},
  {"ac",42227,0, 18092,0,"m:16: individual stars"},
  {"ac",42227,0,  4370,1,"m:13: scratch"},
  {"ac",42227,0, 24964,1,"m:26: scratch"},
  {"ac",42227,0,482058,0,"m:13: stars near limiting magnitude"},
  {"dnb",3393,0,132795,1,"m:16: airplane trail"},
  {"dnb",3393,0,444085,0,"m:10: stars near limiting mag"},
  {"dnb",3393,0,433615,0,"m:11: stars near limiting mag"},
  {"dnb",3393,0,351076,0,"m:16: stars near limiting mag"},
  {"dnb",3393,0,350740,0,"m:53: stars near limiting mag"},
};

int sampleTableSize = sizeof(sampleTable)/sizeof(SAMPLE);

typedef struct _starimage {
  int NUMBER;         /* Sextractor reference number */
  double X_IMAGE;        /* Sextractor X location in pixels */
  double Y_IMAGE;        /* Sextractor Y location in pixels */
  double ra;
  double dec;
  double FLUX_MAX; 
  double FLUX_ISO;
  double MAG_ISO;     /* Sextractor isophotonic magnitude */
  double FWHM_WORLD;
  double THETA_J2000;
  double THETA_IMAGE;
  double ELLIPTICITY;
  double aLength;
  double bLength;
  double plate_dist;
  int ISO0;
  int AFLAGS;          /* flags word */
  int BFLAGS;          /* flags word */
  double drad;        /* Matching error */
  double Stdmag;
  /* Internal flags */
  INTERNAL internal;
} STARIMAGE,*PSTARIMAGE;

typedef struct _scratch_common {
  PSTARIMAGE sextractor_table;
  size_t sextractor_nrecs;

  PCHAIN chain_table;
  int total_chain_count;
  int max_chain_index;

  PNEIGHBOR neighbor_table;
  int max_neighbor_count;
  int neighbor_alloc;

  int iterations;  /* Recursion iterations */
  
} SCRATCHCOMMON,*PSCRATCHCOMMON;


/* Sort routine based on the Sextractor reference number */
int ImageCompare(const void *first, const void *second) 
{
  int numberFirst = ((PSTARIMAGE)first)->NUMBER;
  int numberSecond = ((PSTARIMAGE)second)->NUMBER;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}

void DumpChains(PSCRATCHCOMMON pScratchCommon,PCHAIN pStartChain,int verbose,int lineno)
{

  int chain_index;
  PCHAIN pChain;
  PSTARIMAGE pSextractor;
  PSTARIMAGE pSextractor2;
  int chain_count = 0;
  int chain_index_start = 0;
  int chain_index_end = pScratchCommon->max_chain_index;
  int memberCount;
  PNEIGHBOR pNeighbor;
  printf("Chain Dump max_chain_count %d max_chain_index %d from line %d \n",pScratchCommon->total_chain_count,pScratchCommon->max_chain_index,lineno);
  if (pStartChain != NULL) {
    chain_index_start = pStartChain->chain_index;
    chain_index_end = chain_index_start+1;
  }
  for (chain_index = chain_index_start; chain_index < chain_index_end; chain_index++) {
    pChain = &pScratchCommon->chain_table[chain_index];
    memberCount = 0;
    if (pChain->memberCount != 0) {
      printf("Chain %d memberCount %d aveNeighborCt %f aveLengthRatio %f\n",chain_index,pChain->memberCount,pChain->aveNeighborCt,pChain->aveLengthRatio);
      chain_count++;
      if (pChain->chain_index != chain_index) {
        printf("ERROR in chain %d, wrong chain index %d\n",chain_index,pChain->chain_index);
      }
      pSextractor = pChain->clink;
      while (pSextractor != NULL) {
        memberCount++;
        if (verbose != 0) {
          printf("chain %d member %d neighbors %d NUMBER %d chain %d testDistance %f flink 0x%X\n",
                 chain_index,
                 memberCount,
                 pSextractor->internal.neighbor_count,
                 pSextractor->NUMBER,
                 pSextractor->internal.chain_id,
                 pSextractor->internal.testDist,
                 pSextractor->internal.flink);
        }
        pNeighbor = pSextractor->internal.pNeighbor;
        while (pNeighbor != NULL) {
          pSextractor2 = pNeighbor->pSextractor;
          printf("            neighbor_index %d path_flag %d neighborDist %f NUMBER %d\n",
                 pNeighbor->neighbor_index,
                 pNeighbor->path_flag,
                 pNeighbor->neighborDist,
                 pSextractor2->NUMBER);
          pNeighbor = pNeighbor->pNeighborFlink;
        }

        pSextractor = pSextractor->internal.flink;
      }
      if (memberCount != pChain->memberCount) {
        printf("ERROR in chain %d, member count %d does not match %d\n",pChain->chain_index,memberCount,pChain->memberCount);
      }

    } else {
      if (pChain->clink != NULL) {
        printf("ERROR: clink for chain %d is not null 0x%x\n");
      }
    }
  }
  if (pStartChain == NULL) {
    if (pScratchCommon->total_chain_count != chain_count) {
      printf("ERROR: total chains %d does not equal total_chain_count %d \n",chain_count,pScratchCommon->total_chain_count);
    }
  }
}
/* Clear topology fields and recurse into the next entry */
void ClearFields(PSCRATCHCOMMON pScratchCommon,PSTARIMAGE pSextractor) {
  PNEIGHBOR pNeighbor;
  if (pSextractor->internal.testDist < 0) {
    /* We were here already, discontinue the chain */
    return;
  }
  pScratchCommon->iterations++;
  pSextractor->internal.testDist = -1;
  pSextractor->internal.testDist2 = -1;
  pSextractor->internal.path_flag = 0;
  pNeighbor = pSextractor->internal.pNeighbor;
  while (pNeighbor != NULL) {
    pNeighbor->path_flag = 0;
    ClearFields(pScratchCommon,pNeighbor->pSextractor);
    pNeighbor = pNeighbor->pNeighborFlink;
  }
  return;
}
/* Find the shortest distance to all members of the chain */
void FindDistance(PSCRATCHCOMMON pScratchCommon,PSTARIMAGE pSextractor) {
  PNEIGHBOR pNeighbor;
  PSTARIMAGE pSextractor2;
  pScratchCommon->iterations++;
  double curDistance;
  pNeighbor = pSextractor->internal.pNeighbor;
  while (pNeighbor != NULL) {
    curDistance = pSextractor->internal.testDist + pNeighbor->neighborDist;
    pSextractor2 = pNeighbor->pSextractor;
    if ((pSextractor2->internal.testDist < 0) ||
        (curDistance < pSextractor2->internal.testDist)) {
      pSextractor2->internal.testDist = curDistance;
      FindDistance(pScratchCommon,pSextractor2);
    }
    pNeighbor = pNeighbor->pNeighborFlink;
  }
  return;
}

/* We know the shortest path between two points, now find the path which produces that distance */
int FindPath(PSCRATCHCOMMON pScratchCommon,PSTARIMAGE pSextractor,double curMaxDistance) {
  PNEIGHBOR pNeighbor;
  PSTARIMAGE pSextractor2;
  pScratchCommon->iterations++;
  double curDistance;
  int path_flag = 0;
  pNeighbor = pSextractor->internal.pNeighbor;
  while (pNeighbor != NULL) {
    pSextractor2 = pNeighbor->pSextractor;
    curDistance = pSextractor->internal.testDist2 + pNeighbor->neighborDist;
    if (fabs(curDistance - curMaxDistance) < 0.001) {
      pSextractor2->internal.testDist2 = curDistance;
      pSextractor2->internal.path_flag = 1;
      pSextractor->internal.path_flag = 1;
      pNeighbor->path_flag = 1;
      return(1);
    }
    if ((pSextractor2->internal.testDist2 < 0) ||
        (curDistance < pSextractor2->internal.testDist2)) {
      pSextractor2->internal.testDist2 = curDistance;
      if (fabs(curDistance-pSextractor2->internal.testDist) < 0.001) {
        path_flag = FindPath(pScratchCommon,pSextractor2,curMaxDistance);
        if (path_flag == 1) {
          pSextractor2->internal.path_flag = 1;
          pSextractor->internal.path_flag = 1;
          pNeighbor->path_flag = 1;
          return(1);
        }
      }
    }
    pNeighbor = pNeighbor->pNeighborFlink;
  }
  return(path_flag);
}

/*
 *  This routine investigates the chain topology, ignoring any curves in the chain of images.
 *  For every line connecting any two chain members along a path through neighboring chain centers, find the shortest distance.  Find the two chain members whose shortest distance 
 *  is longer than that of any other member. For a straight line scratch, these two members will be at each end of the scratch
 */
void    MapChainTopology(PSCRATCHCOMMON pScratchCommon,PCHAIN pChain) {
  PSTARIMAGE pSextractor = NULL;
  PSTARIMAGE pRootSextractor = NULL;
  PSTARIMAGE pMaxSextractor;
  PNEIGHBOR pNeighbor;
  int neighbor_index;
  int NUMBER;
  double curMaxDistance;
  double bestMaxDistance = 0;
  double curDistance;
  int path_flag;
  /* Initialize chain variables */  
  pChain->max_separation = 0;
  pChain->max_separation_count = 0;
  pChain->max_separation_flag = 0;
  if (pChain->memberCount <= 2) {
    return;
  }
  /* Select the base of the chain as our starting test point */
  pRootSextractor = pChain->clink;
  /* Zero out the neighbor test fields */
  while (1) {
    ClearFields(pScratchCommon,pRootSextractor);
    pRootSextractor->internal.testDist = 0.0;
    pScratchCommon->iterations = 0;
    FindDistance(pScratchCommon,pRootSextractor);
#ifdef LOS_DEBUG_CHAIN
    if (pChain->chain_index == LOS_DEBUG_CHAIN) {
      DumpChains(pScratchCommon,pChain,1,__LINE__);
    }
#endif /* LOS_DEBUG_CHAIN */
    /* Now find the current maximum distance */
    curMaxDistance = 0;
    pMaxSextractor = pRootSextractor;
    pSextractor = pChain->clink;
    while (pSextractor != NULL) {
      if (pSextractor->internal.testDist > curMaxDistance) {
        curMaxDistance = pSextractor->internal.testDist;
        pMaxSextractor = pSextractor;
      }
      pSextractor = pSextractor->internal.flink;
    }
    if (curMaxDistance <= bestMaxDistance) {
      break;
    }
    /* pMaxSextractor is an endpoint.  See if we can go further from this endpoint */
    bestMaxDistance = curMaxDistance;
    pRootSextractor = pMaxSextractor;  
  }
  /* Now trace the best path, setting the path_flag on the way back */
  pRootSextractor->internal.testDist2 = 0.0;
  path_flag = FindPath(pScratchCommon,pRootSextractor,curMaxDistance);
  /* Count the members on the shortest path */
  pChain->pathCount = 0;
  pChain->matchCount = 0;
  pSextractor = pChain->clink;
  while (pSextractor != NULL) {
    if (pSextractor->internal.path_flag != 0) {
      pChain->pathCount++;
    }
    if (pSextractor->internal.flag_match != 0) {
      pChain->matchCount++;
    }
    pSextractor = pSextractor->internal.flink;
  }

#ifdef LOS_DEBUG_CHAIN
  if (pChain->chain_index == LOS_DEBUG_CHAIN) {
    printf("Max distance %f for chain %d\n",curMaxDistance,pChain->chain_index);
    printf("iterations %d\n",pScratchCommon->iterations);
    DumpChains(pScratchCommon,pChain,1,__LINE__);

    pSextractor = pRootSextractor;
    curDistance = 0;
    while (pSextractor != NULL) {
      pNeighbor = pSextractor->internal.pNeighbor;
      NUMBER = pSextractor->NUMBER;
      pSextractor = NULL;
      while (pNeighbor != NULL) {
        if (pNeighbor->path_flag != 0) {
          curDistance += pNeighbor->neighborDist;
          pSextractor = pNeighbor->pSextractor;
          printf("curDistance %f for NUMBER %d to NUMBER %d\n",curDistance,NUMBER,pSextractor->NUMBER);
          break;
        }
        pNeighbor = pNeighbor->pNeighborFlink;
      }
    }
    printf("path_flag %d\n",path_flag);
  }
#endif /* LOS_DEBUG_CHAIN */
  if (path_flag == 0) {
    printf("ERROR: no path_flag returned for chain %d\n",pChain->chain_index);
    exit(-1);
  }
}



int main(int argc,char *argv[])
{
  char *argstr;
  char cmdchar;
  time_t startTime;
  time_t curTime;
  int verbose = 0;
  int doPlots = 0;
  int errorFlag = 0;
  char fileroot[MAX_BUFFER];
  int nvals;
  FILE *summaryHandle = NULL;
  char summaryfile[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int *matchVector;
  int matchVectorCount = 0;
  char outfile[MAX_BUFFER];
  char regionfile[MAX_BUFFER];
  int outCount = 0;
  int candidateCount = 0; 
  int overlapCount = 0;
  FILE *outHandle = NULL;
  FILE *regionHandle = NULL;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int printFlag = 0;
  File match_handle = NULL;
  char match_name[MAX_BUFFER];
  TableHead match_header = NULL;
  PSTARIMAGE match_table = NULL;
  size_t match_nrecs = 0;
  int match_index;
  PSTARIMAGE pMatch;
  char textString[MAX_BUFFER];

  PNEIGHBOR pNeighbor;
  PNEIGHBOR pNeighbor2;

  File sextractor_handle = NULL;
  char sextractor_name[MAX_BUFFER];
  TableHead sextractor_header = NULL;
  int sextractor_index;
  int chain_index;
  PSTARIMAGE pSextractor;
  PSTARIMAGE pSextractor2;
  PSTARIMAGE pSextractor3;
  PSTARIMAGE pSextractor4;
  PSTARIMAGE pBrightestSextractor;
  int index;
  double *vector = NULL;
  int vector_size = 0;
  int vector_index;
  double *vector2[MAX_SPATIAL_BINS+1] = {NULL};
  double *vector3;
  int vector2_size[MAX_SPATIAL_BINS+1] ={0};
  int vector2_index;
  int curClipCount;
  double MAG_ISO_med[MAX_SPATIAL_BINS+1] = {0};
  double MAG_ISO_rms[MAX_SPATIAL_BINS+1] = {0};
  double edgeDist;
  double medLengthRatio = 0.0;
  double rmsLengthRatio = 0.0;
  double pct90LengthRatio = 0.0;
  PCHAIN pChain;
  PCHAIN pChain2;

  double MAG_ISO_radius;
  double magIsoDiff = MAG_ISO_DIFF;

#if 0
  double max_plate_dist = 0;
  double min_MAG_ISO = 99.0;
  double max_MAG_ISO = 99.0;
  int MAG_ISO_unmatched[MAG_ISO_LEVELS+1];
  int MAG_ISO_matched[MAG_ISO_LEVELS+1];
  int MAG_ISO_unmatched_total = 0;
  int MAG_ISO_matched_total = 0;
  double MAG_ISO_increment;
  double MAG_ISO_oldratio = 0;
  double MAG_ISO_oldlimit;
  double MAG_ISO_ratio;
  double MAG_ISO_reqRatio = 2.0; /* desired ratio of unmatched to matched */
#endif


  double MAG_ISO_limit[MAX_SPATIAL_BINS+1] = {0};
  double arcsecPerPixel = 0.0;
  double rightAscension = -99.0;
  double declination = -99.0;
  double matchRadius = -1.0; /* Degrees */
  double factor = 99.0;
  double objectSizeReduction = 0.5; /* Object size reduction factor */

  int mosaicWidth = 0;
  int mosaicHeight = 0;
  int fullMosaicMode = 0;  /* non-production mode or demonstration mode generates a DS9 region file centered around
                            * <right ascension> and <declination> with radius <match radius> 
                            */
  int xsection;
  int ysection;
  int xsection2;
  int ysection2;
  int section_index;  /* Plate section index = xsection + (NUM_SECTIONS + ysection) */
  int section_index2;
  SECTION section_table[TOTAL_SECTIONS];
  double section_area = 0.0;
  PSECTION pSection;
  PSECTION pSection2;
  void *ptree; /* kdtree for images */
  struct kdres *presults;
  int resultSize = 0;
  double pt[2];
  double pos[2]; /* position index into the kdtree */
  double search_radius_pixels;

  SCRATCHCOMMON scratch_common;
  PSCRATCHCOMMON pScratchCommon = &scratch_common;

  PSAMPLE pSample;
  int sample_index;
  char sample_text[40];
  int correctFlag;
  int defectFlag;
  int sample_mosaic_count;
  int sample_mosaic_found = 0;

  PTESTPLATE pTestPlate;
  int testPlateIndex;
  int spatial_bin;
  int catalogNumber = 0;
  char catalogString[MAX_BUFFER];
  time(&startTime);
  ptree = kd_create(2);
  match_name[0] = 0;
  sextractor_name[0] = 0;
  outfile[0] = 0;
  regionfile[0] = 0;
  fileroot[0] = 0;
  qualifier[0] = 0;
  summaryfile[0] = 0;
  memset(section_table,0,sizeof(section_table));
  memset(pScratchCommon,0,sizeof(SCRATCHCOMMON));
  pScratchCommon->max_chain_index = 1;
#ifdef DEBUG_SEXTRACTOR
  printf("sizeof(SAMPLE) %d  sizeof(sampleTable) %d sampleTableSize %d \n",sizeof(SAMPLE),sizeof(sampleTable),sampleTableSize);
  printf("ERROR: DEBUG_SEXTRACTOR is set\n");
#endif /* DEBUG_SEXTRACTOR */
#if 0
  for (sample_index = 0; sample_index < sampleTableSize; sample_index++) {
    pSample = &sampleTable[sample_index];
    printf("SELECT series,plateNumber,spatial_bin,limiting_iso FROM spatialbin%s INNER JOIN photseries USING (seriesId) WHERE series = '%s' AND plateNumber = %d and solutionNumber = 0;\n",catalogString,pSample->series,pSample->plateNumber);
    printf("Sample for %s%05d_%02d NUMBER %6d scratch %d %s\n",
           pSample->series,
           pSample->plateNumber,
           pSample->mosaicNumber,
           pSample->NUMBER,
           pSample->scratchFlag,
           pSample->notes);
  }
  exit(-1);
#endif
#if 0
  for (testPlateIndex = 0; testPlateIndex < testPlateTableSize; testPlateIndex++) {
    pTestPlate = &testPlateTable[testPlateIndex];
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      if (pTestPlate->bin_MAG_ISO[spatial_bin] != 0.0) {
        printf("%5s %5d %d %.2f\n",
               pTestPlate->series,
               pTestPlate->plateNumber,
               spatial_bin,
               pTestPlate->bin_MAG_ISO[spatial_bin]);
      }
    }
  }
  exit(-1);
#endif
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This is a stray argument */
      errorFlag = 1;
      printf("ERROR: argument %s does not have a qualifier\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'q': /* file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(qualifier,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: filter_scratch MAX_BUFFER too small for %s\n",*argv);
              exit(-1);
            }
            catalogNumber = GetCatalogNumber(qualifier);
            if (catalogNumber < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            }
            if (catalogNumber != 0) {
            sprintf(catalogString,"%d",catalogNumber);
            } else {
            catalogString[0] = 0;
            }
          }
          break;



        case 'm': /* match file name */
        case 'M':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(match_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 't': /* arcsec/pixel */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&arcsecPerPixel);
            if (nvals != 1) {
              printf("ERROR: Can not decode number of arcsec per pixel\n");
              errorFlag = 1;
            }
          }
          break;

        case 's': /* sextractor file name */
        case 'S':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(sextractor_name,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'p': /* file root */
        case 'P':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(fileroot,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        case 'w': /* mosaic width */
        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic width %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'h': /* mosaic height */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&mosaicHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode mosaic height %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        

        case 'a': /* matchRadius */
        case 'A':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&matchRadius);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the matchRadius %s\n",*argv);
              errorFlag = 1;
            } else if (strstr(*argv,":") || strstr(*argv,",") || strstr(*argv," ")) {
              printf("ERROR: matchRadius %s is not in decimal degrees\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 'b': /* magIsoDiff */
        case 'B':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&magIsoDiff);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the magIsoDiff %s\n",*argv);
              errorFlag = 1;
            } else if (strstr(*argv,":") || strstr(*argv,",") || strstr(*argv," ")) {
              printf("ERROR: magIsoDiff %s is not in decimal degrees\n",*argv);
              errorFlag = 1;
            }
          }
          break;



        case 'c': /* objectSizeReduction */
        case 'C':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&objectSizeReduction);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the objectSizeReduction %s\n",*argv);
              errorFlag = 1;
            } else if (strstr(*argv,":") || strstr(*argv,",") || strstr(*argv," ")) {
              printf("ERROR: objectSizeReduction %s is not in decimal degrees\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'r': /* rightAscension */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&rightAscension);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the rightAscension %s\n",*argv);
              errorFlag = 1;
            } else if (strstr(*argv,":") || strstr(*argv,",") || strstr(*argv," ")) {
              printf("ERROR: rightAscension %s is not in decimal degrees\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 'd': /* declination */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&declination);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the declination %s\n",*argv);
              errorFlag = 1;
            } else if (strstr(*argv,":") || strstr(*argv,",") || strstr(*argv," ")) {
              printf("ERROR: declination %s is not in decimal degrees\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 'f': /* verbose */
        case 'F':
          fullMosaicMode = 1;
          break;


        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;

        case 'g': /* do plots */
        case 'G':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(regionfile,*++argv,MAX_BUFFER-2);
            if (strlen(*argv) >= MAX_BUFFER-2) {
              printf("ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;

        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }
  /* Verify that we have everything */
  if (outfile[0] == 0) {
    printf("ERROR: No output filename was specified\n");
    errorFlag = 1;
  }
  if (match_name[0] == 0) {
    printf("ERROR: No match filename was specified\n");
    errorFlag = 1;
  }
  if (sextractor_name[0] == 0) {
    printf("ERROR: No sextractor filename was specified\n");
    errorFlag = 1;
  }
  if (fileroot[0] == 0) {
    printf("ERROR: No file root was specified\n");
    errorFlag = 1;
  } else {
    if (ParseFilename(fileroot,series,&plateNumber,&mosaicNumber,&binning,&rotation) == 0) {
      printf("ERROR: Failed to parse file root %s\n",fileroot);
      errorFlag = 1;
    }
  }



  if (arcsecPerPixel == 0) {
    printf("ERROR: no plate scale specified\n");
    errorFlag = 1;
  }
  if (fullMosaicMode == 0) {
    if (rightAscension < -90) {
      printf("ERROR: no right ascension specified\n");
      errorFlag = 1;
    }
    if (declination < -90) {
      printf("ERROR: no declination specified\n");
      errorFlag = 1;
    }
    if (matchRadius <= 0.0) {
      printf("ERROR: no match raidus specified\n");
      errorFlag = 1;
    }
  }
  if ((rightAscension >= -90) &&
      (declination >= -90) &&
      (matchRadius > 0.0)) {
    factor = cos(DEGREES_TO_RAD * declination);
    if (factor == 0) {
      printf("ERROR: declination can not be +/-90 degrees\n");
      errorFlag = 1;
    }
    if ((declination+matchRadius > 90.0) ||
        (declination-matchRadius < -90.0)) {
      printf("ERROR: range exceeds 90 degrees or -90 degrees declination\n");
      errorFlag = 1;
    }
  }
  if (fullMosaicMode != 0) {
    if (mosaicWidth == 0) {
      printf("ERROR: No mosaic width specified\n");
      errorFlag = 1;
    }
    if (mosaicHeight == 0) {
      printf("ERROR: No mosaic height specified\n");
      errorFlag = 1;
    }
    section_area = mosaicWidth*mosaicHeight/(1.0 * TOTAL_SECTIONS);
  }

  /* Open the sextractor file */
  sextractor_handle = Open(sextractor_name,"r");
  if (sextractor_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the sextractor file %s\n",sextractor_name);
  } else {
    if (verbose) {
      printf("Found sextractor file %s\n",sextractor_name);
    }
  }

  /* Open the match file */
  match_handle = Open(match_name,"r");
  if (match_handle == NULL) {
    errorFlag = 1;
    printf("ERROR: Failed to find the match file %s\n",match_name);
  } else {
    if (verbose) {
      printf("Found match file %s\n",match_name);
    }
  }

  outHandle = fopen(outfile,"wt");
  if (outHandle == NULL) {
    printf("ERROR: Failed to open the output file %s\n",outfile);
    errorFlag = 1;
  } else {
    if (verbose) {
      printf("Output file %s\n",outfile);
    }
  }

  fprintf(outHandle,"series\tplateNumber\tmosaicNumber\tmemberCount\taveNeighborCt\taveLengthRatio\n");
  fprintf(outHandle,"------\t-----------\t------------\t-----------\t-------------\t--------------\n");

  regionHandle = fopen(regionfile,"wt");
  if (regionHandle == NULL) {
    printf("ERROR: Failed to open the region file %s\n",regionfile);
    errorFlag = 1;
  } else {
    if (verbose) {
      printf("Region file %s\n",regionfile);
    }
  }



  if (errorFlag) {
    printf("Usage: filter_scratch [arglist]\n");
    printf("                  Non-production mode parameters:\n");
    printf("                      -r <right ascension (degrees)>\n");
    printf("                      -d <declination (degrees)>\n");
    printf("                      -a <match radius (degrees)>\n");
    printf("                      -b <magIsoDiff: adjustment to MAG_ISO cutoff\n");
    printf("                  Production Mode parameters:\n");
    printf("                      -w <mosaic width (pixels)> \n");
    printf("                      -h <mosaic height (pixels)> \n");
    printf("                  Other parameters:\n");
    printf("                      -c <object size factor>\n");
    printf("                      -g <region file>\n");
    printf("                      -m <match file> \n");
    printf("                      -o <output file>\n");
    printf("                      -p <file root> \n");
    printf("                      -q <catalog qualifier>\n");
    printf("                      -s <sextractor file> \n");
    printf("                      -t <plate scale (arcsec/pix)> \n");
    printf("                      -f full mosaic mode (production)\n");
    printf("                      -v  verbose\n");
    return(-1);
  }

  printf("filter_scratch of %s %s, magIsoDiff %f for %s\n",
         __DATE__,__TIME__,magIsoDiff,fileroot);
  /* Now read in the sextractor file */

  sextractor_header = table_header(sextractor_handle,TABLE_PARSE);
  if (sextractor_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",sextractor_name);
    return(-1);
  }

  pScratchCommon->sextractor_table = table_loadva(sextractor_handle,
                                                  &sextractor_header,
                                                  NULL, /* hbase */
                                                  NULL, /* rows */
                                                  NULL,
                                                  sizeof(STARIMAGE),
                                                  &pScratchCommon->sextractor_nrecs,
                                                  TblInt,"NUMBER",TblOff(PSTARIMAGE,NUMBER),
                                                  TblDbl,"X_IMAGE",TblOff(PSTARIMAGE,X_IMAGE),
                                                  TblDbl,"Y_IMAGE",TblOff(PSTARIMAGE,Y_IMAGE),
                                                  TblDbl,"ra",TblOff(PSTARIMAGE,ra),
                                                  TblDbl,"dec",TblOff(PSTARIMAGE,dec),
                                                  TblDbl,"FLUX_ISO",TblOff(PSTARIMAGE,FLUX_ISO),
                                                  TblDbl,"FLUX_MAX",TblOff(PSTARIMAGE,FLUX_MAX),
                                                  TblDbl,"MAG_ISO",TblOff(PSTARIMAGE,MAG_ISO),
                                                  TblDbl,"FWHM_WORLD",TblOff(PSTARIMAGE,FWHM_WORLD),
                                                  TblDbl,"THETA_J2000",TblOff(PSTARIMAGE,THETA_J2000),
                                                  TblDbl,"THETA_IMAGE",TblOff(PSTARIMAGE,THETA_IMAGE),
                                                  TblDbl,"ELLIPTICITY",TblOff(PSTARIMAGE,ELLIPTICITY),
                                                  TblDbl,"aLength",TblOff(PSTARIMAGE,aLength),
                                                  TblDbl,"bLength",TblOff(PSTARIMAGE,bLength),
                                                  TblDbl,"plate_dist",TblOff(PSTARIMAGE,plate_dist),
                                                  TblInt,"ISO0",TblOff(PSTARIMAGE,ISO0),

                                                  TblInt,"AFLAGS",TblOff(PSTARIMAGE,AFLAGS),
                                                  TblInt,"BFLAGS",TblOff(PSTARIMAGE,BFLAGS),

                                                  0,"end",0);
  if (pScratchCommon->sextractor_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",sextractor_name);
    return(-1);
  }
  if (verbose) {
    printf("read %d records for %s\n",pScratchCommon->sextractor_nrecs,sextractor_name);
  }
  pScratchCommon->chain_table = (PCHAIN)calloc(pScratchCommon->sextractor_nrecs,sizeof(CHAIN));
  if (pScratchCommon->chain_table == NULL) {
    printf("ERROR allocating chain_table of %d records\n",pScratchCommon->sextractor_nrecs);
    exit(-1);
  }
  vector = (double *)calloc(pScratchCommon->sextractor_nrecs,sizeof(double));
  if (vector == NULL) {
    printf("ERROR allocating vector of size %d\n",pScratchCommon->sextractor_nrecs);
  }
  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    vector2[spatial_bin] = (double *)calloc(pScratchCommon->sextractor_nrecs,sizeof(double));
    if (vector2[spatial_bin] == NULL) {
      printf("ERROR allocating vector2[%d] of size %d\n",spatial_bin,pScratchCommon->sextractor_nrecs);
    }
  }

  pScratchCommon->neighbor_table = (PNEIGHBOR)calloc(2*pScratchCommon->sextractor_nrecs,sizeof(NEIGHBOR));
  if (pScratchCommon->neighbor_table == NULL) {
    printf("ERROR allocating neighbor_table of %d records\n",pScratchCommon->sextractor_nrecs);
  }
  pScratchCommon->neighbor_alloc = 2*pScratchCommon->sextractor_nrecs;

  /* Sort the table according to sextractor NUMBER */
  qsort((void*)pScratchCommon->sextractor_table,pScratchCommon->sextractor_nrecs,sizeof(STARIMAGE),ImageCompare);
  for (sextractor_index = 0; sextractor_index < pScratchCommon->sextractor_nrecs; sextractor_index++) {
    pSextractor = &pScratchCommon->sextractor_table[sextractor_index];
#if 0
    if ((pSextractor->ra > 216.48601) &&
        (pSextractor->ra < 216.52824) &&
        (pSextractor->dec > 23.7919) &&
        (pSextractor->dec < 23.8277)) {
      printf("sextractor_index %d NUMBER %6d ra %f dec %f\n",
             sextractor_index,
             pSextractor->NUMBER,
             pSextractor->ra,
             pSextractor->dec);
    }
#endif


    if (pSextractor->ISO0 <= 0) {
      pSextractor->ISO0 = 1;
    }
    if (pSextractor->FWHM_WORLD <= 0) {
      pSextractor->FWHM_WORLD = 1.0/3600.;
    }

    pSextractor->drad = 999.;        /* Matching error */
    pSextractor->Stdmag = 99.;       /* Catalog magnitude */
    memset(&pSextractor->internal,0,sizeof(INTERNAL));
    pSextractor->internal.sextractor_index = sextractor_index;
    pSextractor->internal.section_index = -1;
    pSextractor->internal.spatial_bin = CalculateBin(mosaicWidth,mosaicHeight,pSextractor->X_IMAGE,pSextractor->Y_IMAGE,&edgeDist);
  }


#ifdef DEBUG_SEXTRACTOR
  if (regionHandle != NULL) {
    for (sextractor_index = 0; sextractor_index < pScratchCommon->sextractor_nrecs; sextractor_index++) {
      pSextractor = &pScratchCommon->sextractor_table[sextractor_index];
      sprintf(textString," text = {N: %d MAG_ISO %.1f}",pSextractor->NUMBER,pSextractor->MAG_ISO);
      fprintf(regionHandle,"CIRCLE(%f,%f,3) # color = yellow %s\n",
              pSextractor->X_IMAGE,
              pSextractor->Y_IMAGE,
              textString);
    }
  }
  exit(-1);
#endif /* DEBUG_SEXTRACTOR */

  matchVector = (int *)calloc(pScratchCommon->sextractor_nrecs,sizeof(int));
  if (matchVector == NULL) {
    printf("ERROR: failed to allocate matchVector for %d images\n",pScratchCommon->sextractor_nrecs);
    exit(-1);
  }



  /* Now read in the match file */

  match_header = table_header(match_handle,TABLE_PARSE);
  if (match_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",match_name);
    return(-1);
  }

  match_table = table_loadva(match_handle,
                             &match_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(STARIMAGE),
                             &match_nrecs,
                             TblInt,"NUMBER",TblOff(PSTARIMAGE,NUMBER),
                             TblDbl,"drad",TblOff(PSTARIMAGE,drad),
                             TblDbl,"Stdmag",TblOff(PSTARIMAGE,Stdmag),
                             TblInt,"AFLAGS",TblOff(PSTARIMAGE,AFLAGS),
                             TblInt,"BFLAGS",TblOff(PSTARIMAGE,BFLAGS),
                             0,"end",0);
  if (match_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",match_name);
    return(-1);
  }
  if (verbose) {
    printf("read %d records for %s\n",match_nrecs,match_name);
  }
  /* Sort the table according to sextractor NUMBER */
  qsort((void*)match_table,match_nrecs,sizeof(STARIMAGE),ImageCompare);

  /* Now we step through the match table and copy matches to SExtractor */
  match_index = 0;
  pMatch = &match_table[match_index];
  for (sextractor_index = 0; sextractor_index < pScratchCommon->sextractor_nrecs; sextractor_index++) {
    pSextractor = &pScratchCommon->sextractor_table[sextractor_index];
    vector3 = vector2[pSextractor->internal.spatial_bin];
    vector3[vector2_size[pSextractor->internal.spatial_bin]] = pSextractor->MAG_ISO;
    vector2_size[pSextractor->internal.spatial_bin]++;
    
    while (pMatch->NUMBER < pSextractor->NUMBER) {
      match_index++;
      if (match_index >= match_nrecs) {
        break;
      }
      pMatch = &match_table[match_index];
    }
    if (pSextractor->NUMBER == pMatch->NUMBER) {
      /* Copy match data to SExtractor table */
      pSextractor->drad = pMatch->drad;
      pSextractor->Stdmag = pMatch->Stdmag;
      pSextractor->AFLAGS = pMatch->AFLAGS;
      pSextractor->BFLAGS = pMatch->BFLAGS;
      pSextractor->internal.flag_match = 1;
      matchVector[matchVectorCount] = sextractor_index;
      matchVectorCount++;
      vector[vector_size] = pSextractor->aLength/pSextractor->bLength;
      vector_size++;
    }
  }
  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    curClipCount = vector2_size[spatial_bin];
    if (curClipCount == 0) {
      MAG_ISO_med[spatial_bin] = 0;
      MAG_ISO_rms[spatial_bin] = 0;
    } else {
      curClipCount = CalcMedianAndRMS(curClipCount,10,vector2[spatial_bin],&MAG_ISO_med[spatial_bin],&MAG_ISO_rms[spatial_bin],1,CLIP_RMS_FACTOR,0);
      if (curClipCount > 0) {
        curClipCount = CalcMedianAndRMS(curClipCount,10,vector2[spatial_bin],&MAG_ISO_med[spatial_bin],&MAG_ISO_rms[spatial_bin],1,CLIP_RMS_FACTOR,0);
      }
    }
    if (curClipCount == 0) {
      MAG_ISO_med[spatial_bin] = 0;
      MAG_ISO_rms[spatial_bin] = 0;
    } else {

      for (testPlateIndex = 0; testPlateIndex < testPlateTableSize; testPlateIndex++) {
        pTestPlate = &testPlateTable[testPlateIndex];
        if ((pTestPlate->plateNumber == plateNumber)  &&
            (strcmp(pTestPlate->series,series) == 0)) {
          break;
        }
      }
      if (testPlateIndex >= testPlateTableSize) {
        if (spatial_bin == 1) {
          printf("ERROR: no testPlate entry MAG_ISO_med %.2f MAG_ISO_rms %.2f for %s\n",MAG_ISO_med[spatial_bin],MAG_ISO_rms[spatial_bin],fileroot);
        }
      } else {
        if (pTestPlate->bin_MAG_ISO[spatial_bin] != 0.0) {

          printf("MAG_ISO_med %6.2f MAG_ISO_rms %6.2f limiting bin_MAG_ISO %6.2f rms_MAG_ISO %.2f factor %7.3f diff %.2f spatial_bin %d for %s\n",
                 MAG_ISO_med[spatial_bin],
                 MAG_ISO_rms[spatial_bin],
                 pTestPlate->bin_MAG_ISO[spatial_bin],
                 pTestPlate->rms_MAG_ISO,
                 (pTestPlate->bin_MAG_ISO[spatial_bin] - MAG_ISO_med[spatial_bin])/MAG_ISO_rms[spatial_bin],
                 (pTestPlate->bin_MAG_ISO[spatial_bin] - MAG_ISO_med[spatial_bin]),
                 spatial_bin,
                 fileroot);
        }
      }
    }
  }



  if (CalcMedianAndRMS(vector_size,0,vector,&medLengthRatio,&rmsLengthRatio,0,3.0,0) > 0) {
    pct90LengthRatio = vector[(vector_size*9)/10];
#if 0
    for (vector_index = 0; vector_index < vector_size; vector_index++) {
      printf("%d %.1f\n",vector_index,vector[vector_index]);
    }
#endif
  }
#if 0 /* Test of photometry pipeline MAG_ISO values */
  for (testPlateIndex = 0; testPlateIndex < testPlateTableSize; testPlateIndex++) {
    pTestPlate = &testPlateTable[testPlateIndex];
    if ((pTestPlate->plateNumber == plateNumber)  &&
        (strcmp(pTestPlate->series,series) == 0)) {
      break;
    }
  }
  if (testPlateIndex >= testPlateTableSize) {
    printf("ERROR: plate %s is not in the testPlateTable\n",fileroot);
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      MAG_ISO_limit[spatial_bin] = -999.00;
      printf("WARNING: Setting MAG_ISO_limit to %f from bin_MAG_ISO %f and magIsoDiff %f for spatial_bin %d of %s\n",MAG_ISO_limit[spatial_bin],pTestPlate->bin_MAG_ISO[spatial_bin],magIsoDiff,spatial_bin,fileroot);
    }
  } else {
    for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
      if (pTestPlate->bin_MAG_ISO[spatial_bin] != 0) {
        MAG_ISO_limit[spatial_bin] = pTestPlate->bin_MAG_ISO[spatial_bin] + magIsoDiff;
        printf("Setting MAG_ISO_limit to %f from bin_MAG_ISO %f and magIsoDiff %f for spatial_bin %d of %s\n",MAG_ISO_limit[spatial_bin],pTestPlate->bin_MAG_ISO[spatial_bin],magIsoDiff,spatial_bin,fileroot);
      } else {
        MAG_ISO_limit[spatial_bin] = -999.00;
        printf("WARNING: Setting MAG_ISO_limit to %f from bin_MAG_ISO %f and magIsoDiff %f for spatial_bin %d of %s\n",MAG_ISO_limit[spatial_bin],pTestPlate->bin_MAG_ISO[spatial_bin],magIsoDiff,spatial_bin,fileroot);
      }
    }
  }

#endif

#if 1
  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    MAG_ISO_limit[spatial_bin] = MAG_ISO_med[spatial_bin] + magIsoDiff;
    printf("Setting MAG_ISO_limit to %f from MAG_ISO_med %f and magIsoDiff %f for spatial_bin %d of %s\n",MAG_ISO_limit[spatial_bin],MAG_ISO_med[spatial_bin],magIsoDiff,spatial_bin,fileroot);
  }
#endif
  for (sextractor_index = 0; sextractor_index < pScratchCommon->sextractor_nrecs; sextractor_index++) {
    pSextractor = &pScratchCommon->sextractor_table[sextractor_index];
    if (pSextractor->MAG_ISO < MAG_ISO_limit[pSextractor->internal.spatial_bin]) {            
      if (fullMosaicMode == 0) {
        if ((pSextractor->ra < (rightAscension+(matchRadius/factor))) &&
            (pSextractor->ra > (rightAscension-(matchRadius/factor))) &&
            (pSextractor->dec < (declination+matchRadius)) &&
            (pSextractor->dec > (declination-matchRadius))) {
          pSextractor->internal.flag_output = 1;
          outCount++;
        }
      } else { /*  fullMosaicMode */ 
        pSextractor->internal.flag_MAG_ISO = 1;
        candidateCount++;
        xsection = ((1.0 * NUM_SECTIONS) *pSextractor->X_IMAGE)/(1.0 *mosaicWidth);
        if (xsection < 0) {
          xsection = 0;
        }
        if (xsection >= NUM_SECTIONS) {
          xsection = NUM_SECTIONS-1;
        }
        ysection = ((1.0 * NUM_SECTIONS) *pSextractor->Y_IMAGE) /(1.0 * mosaicHeight);
        if (xsection < 0) {
          xsection = 0;
        }
        if (xsection >= NUM_SECTIONS) {
          xsection = NUM_SECTIONS-1;
        }
        
        pSextractor->internal.section_index = xsection + (NUM_SECTIONS * ysection);
        pSection = &section_table[pSextractor->internal.section_index];
        pSection->section_images++;
        if (pSextractor->aLength > pSection->maxaLength) {
          pSection->maxaLength = pSextractor->aLength;
        }
        if (pSextractor->bLength > pSection->maxbLength) {
          pSection->maxbLength = pSextractor->bLength;
        }
        pos[0] = pSextractor->X_IMAGE;
        pos[1] = pSextractor->Y_IMAGE;
        kd_insert(ptree,pos,pSextractor);
      }
    }
  }
  
  if (fullMosaicMode != 0) {
    /* Calculate the section index for the bin of interest */
    for (section_index = 0; section_index < TOTAL_SECTIONS; section_index++) {
      pSection = &section_table[section_index];
      pSection->section_index = section_index;
      if (pSection->section_images > 0) {
        pSection->MAG_ISO_radius_pixels = objectSizeReduction * sqrt((section_area)/(3.14159*pSection->section_images));
      }
    }
    /* Now calculate the maximum search radius for the bin of interest */
    for (section_index = 0; section_index < TOTAL_SECTIONS; section_index++) {
      pSection = &section_table[section_index];
      xsection = section_index % NUM_SECTIONS;      
      ysection = section_index / NUM_SECTIONS;
      for (section_index2 = 0; section_index2 < TOTAL_SECTIONS; section_index2++) {
        if (section_index == section_index2) {
          continue;
        }
        pSection2 = &section_table[section_index2];
        xsection2 = section_index2 % NUM_SECTIONS;      
        ysection2 = section_index2 / NUM_SECTIONS;
        if ((abs(xsection-xsection2) <= 1) &&
            (abs(ysection-ysection2) <= 1)) {
          search_radius_pixels = sqrt(sqr(pSection->maxaLength+pSection->MAG_ISO_radius_pixels+pSection2->maxaLength+pSection2->MAG_ISO_radius_pixels) +
                                      sqr(pSection->maxbLength+pSection->MAG_ISO_radius_pixels+pSection2->maxbLength+pSection2->MAG_ISO_radius_pixels));
          if (search_radius_pixels > pSection->search_radius_pixels) {
            pSection->search_radius_pixels = search_radius_pixels;
          }
        }
      }
    }
  
#if 0
    if (verbose) {
      for (section_index = 0; section_index < TOTAL_SECTIONS; section_index++) {
        pSection = &section_table[section_index];
        printf("section_index %3d, images %5d maxalength %6.1f maxblength %6.1f MAG_ISO_radius_pixels %6.1f search_radius_pixels %6.1f\n",
               section_index,
               pSection->section_images,
               pSection->maxaLength,
               pSection->maxbLength,
               pSection->MAG_ISO_radius_pixels,
               pSection->search_radius_pixels);
   
      }
    }
#endif
    /* Now do the search for any overlapping images */
    for (sextractor_index = 0; sextractor_index < pScratchCommon->sextractor_nrecs; sextractor_index++) {
      pSextractor = &pScratchCommon->sextractor_table[sextractor_index];
      if ((verbose) && ((sextractor_index % 100000) == 0)) {
        time(&curTime);
        curTime -= startTime;
        printf("At sextractor_index %7d of %7d resultSize %d outputCount %d at %d seconds\n",sextractor_index,pScratchCommon->sextractor_nrecs,resultSize,overlapCount,curTime);
      }
      if (pSextractor->internal.flag_MAG_ISO == 0) {
        continue;
      }
      pSection = &section_table[pSextractor->internal.section_index];
      pt[0] = pSextractor->X_IMAGE;
      pt[1] = pSextractor->Y_IMAGE;
      presults = kd_nearest_range(ptree,pt,pSection->search_radius_pixels);
      resultSize = kd_res_size(presults);
      while (!kd_res_end(presults)) {
        pSextractor2 = (PSTARIMAGE)kd_res_item(presults,pos);
        pSection2 = &section_table[pSextractor2->internal.section_index];

        if (pSextractor2->internal.sextractor_index != sextractor_index) {
#if 0
          if (((pSextractor->NUMBER == 175932) &&
               (pSextractor2->NUMBER == 176006)) ||
              ((pSextractor->NUMBER == 176006) &&
               (pSextractor2->NUMBER == 175932))) {
            printf("At NUMBER %6d NUMBER %6d\n",pSextractor->NUMBER,pSextractor2->NUMBER);
          }
#endif

          if (CheckBlend(pSextractor->X_IMAGE,
                         pSextractor->Y_IMAGE,
                         pSextractor->THETA_IMAGE,
                         pSextractor->aLength+pSection->MAG_ISO_radius_pixels,
                         pSextractor->bLength+pSection->MAG_ISO_radius_pixels,
                         pSextractor2->X_IMAGE,
                         pSextractor2->Y_IMAGE,
                         pSextractor2->THETA_IMAGE,
                         pSextractor2->aLength+pSection2->MAG_ISO_radius_pixels,
                         pSextractor2->bLength+pSection2->MAG_ISO_radius_pixels)) {
            /* We have a match! Record it in the neighbor record */
            if (pScratchCommon->max_neighbor_count >= pScratchCommon->neighbor_alloc) {
              printf("ERROR filter_scratch neighbor count %d exceeds neighbor_alloc %d for %s\n",pScratchCommon->max_neighbor_count,pScratchCommon->neighbor_alloc,fileroot);
              exit(-1);
            }
            /* Be sure that this neighbor has not already been added */
            pNeighbor = pSextractor->internal.pNeighbor;
            while (pNeighbor != NULL) {
              if (pNeighbor->pSextractor == pSextractor2) {
                break;
              }
              pNeighbor = pNeighbor->pNeighborFlink;
            }
            /* Add the neighbor */
            if (pNeighbor == NULL) {
              pNeighbor = &pScratchCommon->neighbor_table[pScratchCommon->max_neighbor_count];
              memset(pNeighbor,0,sizeof(NEIGHBOR));
              pNeighbor->neighbor_index = pScratchCommon->max_neighbor_count;
              pNeighbor->neighborDist = sqrt(sqr(pSextractor->X_IMAGE-pSextractor2->X_IMAGE) + sqr(pSextractor->Y_IMAGE-pSextractor2->Y_IMAGE));
              pNeighbor->pSextractor = pSextractor2;
              pNeighbor->pNeighborFlink = pSextractor->internal.pNeighbor;
              pSextractor->internal.pNeighbor = pNeighbor;
              pScratchCommon->max_neighbor_count++;
            }
            /* Now repeat with the matching object */
            pNeighbor = pSextractor2->internal.pNeighbor;
            while (pNeighbor != NULL) {
              if (pNeighbor->pSextractor == pSextractor) {
                break;
              }
              pNeighbor = pNeighbor->pNeighborFlink;
            }
            /* Add the neighbor */
            if (pNeighbor == NULL) {
              pNeighbor = &pScratchCommon->neighbor_table[pScratchCommon->max_neighbor_count];
              memset(pNeighbor,0,sizeof(NEIGHBOR));
              pNeighbor->neighbor_index = pScratchCommon->max_neighbor_count;
              pNeighbor->neighborDist = sqrt(sqr(pSextractor->X_IMAGE-pSextractor2->X_IMAGE) + sqr(pSextractor->Y_IMAGE-pSextractor2->Y_IMAGE));
              pNeighbor->pSextractor = pSextractor;
              pNeighbor->pNeighborFlink = pSextractor2->internal.pNeighbor;
              pSextractor2->internal.pNeighbor = pNeighbor;
              pScratchCommon->max_neighbor_count++;
            }


            /* Either create a new chain, add an unmatched image to a chain, or combine chains */
            overlapCount++;
            pSextractor->internal.neighbor_count++;


            if ((pSextractor->internal.chain_id == 0) && 
                (pSextractor2->internal.chain_id == 0)) {
              /* Neither are chain members, so allocate a new chain and insert them */
              pChain = &pScratchCommon->chain_table[pScratchCommon->max_chain_index];
              pChain->chain_index = pScratchCommon->max_chain_index;
              pChain->clink = pSextractor;
              pSextractor->internal.flink = pSextractor2;
              pChain->memberCount = 2;
              pSextractor->internal.chain_id = pScratchCommon->max_chain_index;
              pSextractor2->internal.chain_id = pScratchCommon->max_chain_index;
              pScratchCommon->total_chain_count++;
              pScratchCommon->max_chain_index++;
            } else if ((pSextractor->internal.chain_id != 0) &&
                       (pSextractor2->internal.chain_id == 0)) {
              /* Add the second object to the tail of the first chain */
              pChain = &pScratchCommon->chain_table[pSextractor->internal.chain_id];
              pChain->memberCount++;
              pSextractor2->internal.chain_id = pSextractor->internal.chain_id;
              pSextractor2->internal.flink = pChain->clink;
              pChain->clink = pSextractor2;
            } else if ((pSextractor2->internal.chain_id != 0) &&
                       (pSextractor->internal.chain_id == 0)) {
              /* Add the first object to the tail of the second chain */
              pChain = &pScratchCommon->chain_table[pSextractor2->internal.chain_id];
              pChain->memberCount++;
              pSextractor->internal.chain_id = pSextractor2->internal.chain_id;
              pSextractor->internal.flink = pChain->clink;
              pChain->clink = pSextractor;
            } else if  (pSextractor2->internal.chain_id != pSextractor->internal.chain_id) {
              /* Members of two different chains.  Delete one chain and move all of its members to the other chain */
              pScratchCommon->total_chain_count--;
              pChain = &pScratchCommon->chain_table[pSextractor->internal.chain_id];
              pChain2 = &pScratchCommon->chain_table[pSextractor2->internal.chain_id];
              pSextractor3 = pChain2->clink;
              pChain2->clink = NULL;
              while (pSextractor3 != NULL) {
                pSextractor4 = pSextractor3->internal.flink;
                pChain2->memberCount--;
                pChain->memberCount++;
                pSextractor3->internal.flink = pChain->clink;
                pChain->clink = pSextractor3;
                pSextractor3->internal.chain_id = pSextractor->internal.chain_id;
                pSextractor3 = pSextractor4;
              }
              if ((pChain2->memberCount != 0) ||
                  (pChain2->clink != NULL)) {
                printf("ERROR: corrupted chain id %d members %d clink 0x%x\n",pChain2->chain_index,pChain2->memberCount,pChain2->clink);
              }
            }


            if ((factor > 90.0) || 
                ((pSextractor->ra < (rightAscension+(matchRadius/factor))) &&
                 (pSextractor->ra > (rightAscension-(matchRadius/factor))) &&
                 (pSextractor->dec < (declination+matchRadius)) &&
                 (pSextractor->dec > (declination-matchRadius)))) {
            }
          }
        }
        /* On to the next item */
        kd_res_next(presults);
      }
      kd_res_free(presults);
    }
  
    if (verbose != 0) {
      time(&curTime);
      curTime -= startTime;
      printf("At sextractor_index %7d of %7d resultSize %d outputCount %d at %d seconds\n",sextractor_index,pScratchCommon->sextractor_nrecs,resultSize,overlapCount,curTime);
    }


    /* At this point, calculate the mean number of neighbors for each chain */
    for (chain_index = 0; chain_index < pScratchCommon->max_chain_index; chain_index++) {
      pChain = &pScratchCommon->chain_table[chain_index];
      if (pChain->memberCount == 0) {
        continue;
      }
      pBrightestSextractor = NULL;
      pChain->aveNeighborCt = 0;
      pChain->aveLengthRatio = 0;
      pSextractor = pChain->clink;
      while (pSextractor != NULL) {
        pChain->aveNeighborCt += 1.0*pSextractor->internal.neighbor_count;
        pChain->aveLengthRatio += 1.0*pSextractor->aLength/pSextractor->bLength;
        pSextractor = pSextractor->internal.flink;
      }
      pChain->aveNeighborCt = pChain->aveNeighborCt/(1.0*pChain->memberCount);
      pChain->aveLengthRatio = pChain->aveLengthRatio/(1.0*pChain->memberCount);
      if (pChain->memberCount >= MINIMUM_CHAIN_COUNT) {
        pChain->selected = 1;
        pSextractor = pChain->clink;
        while (pSextractor != NULL) {
          if (factor < 90.0) {
            if ((pSextractor->ra < (rightAscension+(matchRadius/factor))) &&
                (pSextractor->ra > (rightAscension-(matchRadius/factor))) &&
                (pSextractor->dec < (declination+matchRadius)) &&
                (pSextractor->dec > (declination-matchRadius))) {
              pSextractor->internal.flag_output = 1;
              outCount++;
              if (pBrightestSextractor == NULL) {
                pBrightestSextractor = pSextractor;
              } else if (pSextractor->MAG_ISO < pBrightestSextractor->MAG_ISO) {
                pBrightestSextractor = pSextractor;
              }
            } 
          } else {
            pSextractor->internal.flag_output = 1;
            outCount++;
            if (pBrightestSextractor == NULL) {
              pBrightestSextractor = pSextractor;
            } else if (pSextractor->MAG_ISO < pBrightestSextractor->MAG_ISO) {
              pBrightestSextractor = pSextractor;
            }
          }
          pSextractor = pSextractor->internal.flink;
        }
      }
      if (pBrightestSextractor != NULL) {
        pBrightestSextractor->internal.brightest_member = 1;
      }
    }
    if (verbose != 0) {
      time(&curTime);
      curTime -= startTime;
      printf("Find chain neighbors at %d seconds\n",curTime);
    }

    /* Analyze the chain topology */
    for (chain_index = 0; chain_index < pScratchCommon->max_chain_index; chain_index++) {
      pChain = &pScratchCommon->chain_table[chain_index];
      MapChainTopology(pScratchCommon,pChain);
    }

    if (verbose != 0) {
      time(&curTime);
      curTime -= startTime;
      printf("MapChainTopology at %d seconds\n",curTime);
    }

    for (sextractor_index = 0; sextractor_index < pScratchCommon->sextractor_nrecs; sextractor_index++) {
      pSextractor = &pScratchCommon->sextractor_table[sextractor_index];
#ifdef LOS_DEBUG_NUMBER
      if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
        printf("At NUMBER %6d in line %d\n",pSextractor->NUMBER,__LINE__);
      }
#endif /* LOS_DEBUG_NUMBER */
      if (pSextractor->internal.flag_output == 0) {
        continue;
      }
#if 1
      textString[0] = 0;
      if (pSextractor->internal.section_index >= 0) {
        pSection = &section_table[pSextractor->internal.section_index];
        MAG_ISO_radius = (pSection->MAG_ISO_radius_pixels*arcsecPerPixel)/3600.;
        if (pSextractor->internal.chain_id >= 0) {
          pChain = &pScratchCommon->chain_table[pSextractor->internal.chain_id];
          fprintf(outHandle,"%s\t%d\t%d\t%d\t%.1f\t%.1f\n",
                  series,
                  plateNumber,
                  mosaicNumber,
                  pChain->memberCount,
                  pChain->aveNeighborCt,
                  pChain->aveLengthRatio);

          if (pSextractor->internal.brightest_member != 0) {
            sprintf(textString," text = {m:%d n %.1f a: %.1f N: %d}",pChain->memberCount,pChain->aveNeighborCt,pChain->aveLengthRatio,pSextractor->NUMBER);
          }
          if (
#if 0
              (pChain->pathCount > 18) || 
              (pChain->aveLengthRatio > pct90LengthRatio) ||
              (((1.0*pChain->matchCount)/(1.0*pChain->memberCount)) < 0.40) || 
              (((1.0*pChain->matchCount)/(1.0*pChain->pathCount)) < 0.20)
#else
              (pChain->memberCount >= CYAN_CHAIN_COUNT) || ( pChain->aveLengthRatio >= pct90LengthRatio)
#endif
              ) {
            pSextractor->internal.flag_scratch = 1;
            if (pChain->memberCount >= CYAN_CHAIN_COUNT) {
              pSextractor->internal.flag_chain_length = pChain->memberCount;
            }
            if (pChain->aveLengthRatio >= pct90LengthRatio) {
              pSextractor->internal.flag_length_ratio = 1;
            }
            if (regionHandle != NULL) {
              pNeighbor = pSextractor->internal.pNeighbor;
              while (pNeighbor != NULL) {
                if (pNeighbor->path_flag != 0) {
                  pSextractor2 = pNeighbor->pSextractor;
#ifdef LOS_DEBUG_NUMBER
                  if (pSextractor2->NUMBER == LOS_DEBUG_NUMBER) {
                    printf("At NUMBER %6d in line %d\n",pSextractor2->NUMBER,__LINE__);
                  }
#endif /* LOS_DEBUG_NUMBER */
                  fprintf(regionHandle,"J2000;LINE(%f,%f,%f,%f) # line = 1 1 color = yellow\n",
                          pSextractor->ra,
                          pSextractor->dec,
                          pSextractor2->ra,
                          pSextractor2->dec);
                }
                pNeighbor = pNeighbor->pNeighborFlink;
              }
            }
            if ( pChain->aveLengthRatio >= pct90LengthRatio) {
#ifdef LOS_DEBUG_NUMBER
              if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
                printf("At NUMBER %6d in line %d\n",pSextractor->NUMBER,__LINE__);
              }
#endif /* LOS_DEBUG_NUMBER */
              fprintf(regionHandle,"J2000;BOX(%f,%f,%f,%f,%f) # color = yellow %s\n",
                      pSextractor->ra,
                      pSextractor->dec,
                      2*((pSextractor->aLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                      2*((pSextractor->bLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                      pSextractor->THETA_IMAGE,
                      textString);

            } else if (pChain->memberCount >= CYAN_CHAIN_COUNT) {
              if (pSextractor->internal.path_flag != 0) {
#ifdef LOS_DEBUG_NUMBER
                if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
                  printf("At NUMBER %6d in line %d\n",pSextractor->NUMBER,__LINE__);
                }
#endif /* LOS_DEBUG_NUMBER */
                fprintf(regionHandle,"J2000;BOX(%f,%f,%f,%f,%f) # color = orange %s\n",
                        pSextractor->ra,
                        pSextractor->dec,
                        2*((pSextractor->aLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                        2*((pSextractor->bLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                        pSextractor->THETA_IMAGE,
                        textString);

              } else {
#ifdef LOS_DEBUG_NUMBER
                if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
                  printf("At NUMBER %6d in line %d\n",pSextractor->NUMBER,__LINE__);
                }
#endif /* LOS_DEBUG_NUMBER */
                fprintf(regionHandle,"J2000;BOX(%f,%f,%f,%f,%f) # color = cyan %s\n",
                        pSextractor->ra,
                        pSextractor->dec,
                        2*((pSextractor->aLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                        2*((pSextractor->bLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                        pSextractor->THETA_IMAGE,
                        textString);
              }
            } else {
#if 0
              fprintf(regionHandle,"J2000;BOX(%f,%f,%f,%f,%f) # color = yellow %s\n",
                      pSextractor->ra,
                      pSextractor->dec,
                      2*((pSextractor->aLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                      2*((pSextractor->bLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                      pSextractor->THETA_IMAGE,
                      textString);
#endif
  
            }
          }
        }
      } else {
        if (regionHandle != NULL) {
#ifdef LOS_DEBUG_NUMBER
          if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
            printf("At NUMBER %6d in line %d\n",pSextractor->NUMBER,__LINE__);
          }
#endif /* LOS_DEBUG_NUMBER */
        
          fprintf(regionHandle,"J2000;BOX(%f,%f,%f,%f,%f) # color = yellow %s\n",
                  pSextractor->ra,
                  pSextractor->dec,
                  2*((pSextractor->aLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                  2*((pSextractor->bLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                  pSextractor->THETA_IMAGE,
                  textString);
  
        }
      }
#endif


    } 

    for (sextractor_index = 0; sextractor_index < pScratchCommon->sextractor_nrecs; sextractor_index++) {
      pSextractor = &pScratchCommon->sextractor_table[sextractor_index];
#ifdef LOS_DEBUG_NUMBER
      if (pSextractor->NUMBER == LOS_DEBUG_NUMBER) {
        printf("At NUMBER %6d in line %d\n",pSextractor->NUMBER,__LINE__);
      }
#endif /* LOS_DEBUG_NUMBER */
      sample_mosaic_count = 0;
      for (sample_index = 0; sample_index < sampleTableSize; sample_index++) {
        pSample = &sampleTable[sample_index];
        if ((strcmp(pSample->series,series) == 0)  &&
            (pSample->plateNumber == plateNumber) &&
            (pSample->mosaicNumber == mosaicNumber)) {
          sample_mosaic_count++;
          if (pSextractor->NUMBER == pSample->NUMBER) {
            if (MAG_ISO_limit[spatial_bin] < -990.00) {
              printf("WARNING: sample %s NUMBER %6d is in a bad spatial bin\n,",fileroot,pSextractor->NUMBER);
              continue;
            }

            sample_mosaic_found++;
            if (pSample->scratchFlag != 0) {
              if (pSextractor->internal.flag_scratch != 0) {
                strcpy(sample_text,"Scratch      ");
                correctFlag = 1;
              } else {
                strcpy(sample_text,"FalseNegative");
                correctFlag = 0;
              }
            } else {
              if (pSextractor->internal.flag_scratch!= 0) {
                strcpy(sample_text,"FalsePositive");
                correctFlag = 0;
              } else {
                strcpy(sample_text,"NormalStar   ");
                correctFlag = 1;
              }
            }              
            if ((pSextractor->AFLAGS & (1 << FILTER_AFLAG_DEFECT)) != 0) {
              defectFlag = 1;
            } else {
              defectFlag = 0;
            }
            printf("Sample %s scratchActual %d defect %d scratchDeclared %d chainLength %d lengthRatio %d for %s%05d_%02d NUMBER %6d X %.0f Y %.0f %s\n",
                   sample_text,
                   pSample->scratchFlag,
                   defectFlag,
                   pSextractor->internal.flag_scratch,
                   pSextractor->internal.flag_chain_length,
                   pSextractor->internal.flag_length_ratio,
                   pSample->series,
                   pSample->plateNumber,
                   pSample->mosaicNumber,
                   pSample->NUMBER,
                   pSextractor->X_IMAGE,
                   pSextractor->Y_IMAGE,
                   pSample->notes);
            pChain = &pScratchCommon->chain_table[pSextractor->internal.chain_id];
            if (printFlag == 0) {
              printf("STARBASEflag_MAG_ISO\tdefectFlag\tscratchFlag\tpathCount\tmemberCount\tmatchCount\taveLengthRatio\tNUMBER\tX_IMAGE\tY_IMAGE\tPlate\tnotes\n");
              printf("STARBASE------------\t----------\t-----------\t---------\t-----------\t----------\t--------------\t------\t-------\t-------\t-----\t-----\n");
              printFlag = 1;
            }
            printf("STARBASE%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t%s%05d_%02d\t%s\n",
                   pSextractor->internal.flag_MAG_ISO,
                   defectFlag,
                   pSample->scratchFlag,
                   pChain->pathCount,
                   pChain->memberCount,
                   pChain->matchCount,
                   pChain->aveLengthRatio,
                   pSample->NUMBER,
                   pSextractor->X_IMAGE,
                   pSextractor->Y_IMAGE,
                   pSample->series,
                   pSample->plateNumber,
                   pSample->mosaicNumber,
                   pSample->notes);
                   

          }
        }    
      }
    }
  
    printf("Samples on plate: %d of %d entries, %d not found for %s\n",sample_mosaic_found,sample_mosaic_count,(sample_mosaic_count-sample_mosaic_found),fileroot);

  }


#if 0
  DumpChains(pScratchCommon,NULL,0,__LINE__);
#endif
  if ((fullMosaicMode == 0) && 
      (outCount > 0)) {
    MAG_ISO_radius = objectSizeReduction * sqrt((4*matchRadius*matchRadius)/(3.14159*outCount));
    if (regionHandle != NULL) {
#if 0
      fprintf(regionHandle,"NUMBER\tStdmag\tMAG_ISO\tX_IMAGE\tY_IMAGE\tra\tdec\tQ0\n");
      fprintf(regionHandle,"------\t------\t-------\t-------\t-------\t--\t---\t--\n");
#endif
      for (sextractor_index = 0; sextractor_index < pScratchCommon->sextractor_nrecs; sextractor_index++) {
        pSextractor = &pScratchCommon->sextractor_table[sextractor_index];
        if (pSextractor->internal.flag_output == 0) {
          continue;
        }

#if 0
        fprintf(regionHandle,"j2000;CIRCLE(%f,%f,%f) # color = green \n",pSextractor->ra,pSextractor->dec,MAG_ISO_radius);
#endif
#if 1
        fprintf(regionHandle,"J2000;BOX(%f,%f,%f,%f,%f) # color = red\n",
                pSextractor->ra,
                pSextractor->dec,
                2*((pSextractor->aLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                2*((pSextractor->bLength*arcsecPerPixel/3600.)+MAG_ISO_radius),
                pSextractor->THETA_IMAGE);
  

#endif
#if 0
        fprintf(regionHandle,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",
                pSextractor->NUMBER,
                pSextractor->Stdmag,
                pSextractor->MAG_ISO,
                pSextractor->X_IMAGE,
                pSextractor->Y_IMAGE,
                pSextractor->ra,
                pSextractor->dec,
                pSextractor->internal.flag_match);
#endif
      } 
    }
  }


  /* All done.  Clean up */
 

  if (pScratchCommon->sextractor_table != NULL) {
    Free(pScratchCommon->sextractor_table);
  }
  if (vector != NULL) {
    free(vector);
  }
  for (spatial_bin = 1; spatial_bin <= MAX_SPATIAL_BINS; spatial_bin++) {
    if (vector2[spatial_bin] != NULL) {
      free(vector2[spatial_bin]);
    }
  }
  if (pScratchCommon->chain_table != NULL) {
    free(pScratchCommon->chain_table);
  }
  if (pScratchCommon->neighbor_table != NULL) {
    free(pScratchCommon->neighbor_table);
  }
  if (sextractor_header != NULL) {
    table_hdrfree(sextractor_header);
  }
  if (sextractor_handle != NULL) {
    Close(sextractor_handle);
  }

  if (match_table != NULL) {
    Free(match_table);
  }
  if (match_header != NULL) {
    table_hdrfree(match_header);
  }
  if (match_handle != NULL) {
    Close(match_handle);
  }

  if (outHandle != NULL) {
    fclose(outHandle);
  }

  if (regionHandle != NULL) {
    fclose(regionHandle);
  }

  if (summaryHandle != NULL) {
    fclose(summaryHandle);
  }

  if (matchVector != NULL) {
    free(matchVector);
  }

  time(&curTime);
  curTime -= startTime;
  if (match_nrecs == 0) {
    match_nrecs = 1;
  }

  kd_free(ptree);
  printf("Right Ascension %f,  Declination %f, matchRadius %f objectSizeReduction %f for %s\n",
         rightAscension,
         declination,
         matchRadius,
         objectSizeReduction,
         fileroot);
  printf("Stars %d, matched %d max ratio %f candidates %d overlaps %d output %d  Execution time %d seconds for %s\n",
         pScratchCommon->sextractor_nrecs,
         match_nrecs,
         (1.0*(pScratchCommon->sextractor_nrecs-match_nrecs))/(1.0*match_nrecs),
         candidateCount,
         overlapCount,
         outCount,
         curTime,
         fileroot);
  printf("MAG_ISO_radius %f for %s\n",
         MAG_ISO_radius,
         fileroot);
  printf("medAlength %.1f rmsLengthRatio %.1f pct90LengthRatio %.1f for %s\n",medLengthRatio,rmsLengthRatio,pct90LengthRatio,fileroot);
  printf("max_neighbor_count %d neighbor_alloc %d for %s\n",pScratchCommon->max_neighbor_count,pScratchCommon->neighbor_alloc,fileroot);

  return(0);
}
