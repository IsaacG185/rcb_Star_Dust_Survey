// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatucac5.c
 *
 * See http://cdsarc.u-strasbg.fr/viz-bin/Cat?I/340
 *
 * Read the ucac5 input catalog and put it into the gsc232bin.db format
 *
 * NOTE: maga (UCAC aperture magnitude) is given for Stdmag, no color is given.
 *       the obj_type is copied into the gsc class field.  Only zero obj_types have good astrometry.
 *     
 * NOTE: the highest proper motions are RaPM -4.41 to +6.77 arcsec/year and DecPM -5.82 to +10.33 arcsec/yr. for UCAC4
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm formatucac5.c pipelineutils.a  -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -o formatucac5
 * 
 *  formatucac5 -v -o /dasch/Pipeline/catalogs/ucac5.db
 *
 *   source id field usage in UCAC5: see pipelineutils.h for bit definitions
 *  min_HEALPIX   0x0000000000000001 max_HEALPIX    0x000000000bfffff5
 *  min_PROCCTR   0x0000000000000000 max_PROCCTR    0x0000000000000000
 *  min_SPARE     0x0000000000000000 max_SPARE      0x0000000000000001
 *  min_RUNSEQ    0x0000000000000000 max_RUNSEQ     0x00000000002d9b0a
 *  min_COMPONENT 0x0000000000000000 max_COMPONENT  0x0000000000000000
 *  ALL_MASK      0x7ffffff91fffff80
 *  
 *
 * Feb  7, 2018  Adopted from formatucac4, but must read in two zones because the files are sorted for epoch 2015 instead of epoch 2000
 */   


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#define MAX_INPUT_NAME 512

#define MAGBINS  32  /* Number of magnitude bins */
#define MAXRECORDS (287000) /* Largest file is zone 306 with 286833 records, but we need to read in two zones at once*/
/* #define REFNUMBER_HACK 1 */ /* temporary until GetREFNumber is fixed for Gaia */

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  
/* #define HIGH_PROPER_MOTION_CUTOFF 200.0 */ /* Proper motion cutoff in mas/yr */
/* #define DUMP_PROPER_MOTION 1 */   /* Print proper motion values */
#define MAXSIGMAPM 650  /* Units of 0.1 mas/yr */

#define OBJ_TYPE_CLEAN        1
#define OBJ_TYPE_OVEREXPOSED  2
#define OBJ_TYPE_STREAK       3
#define OBJ_TYPE_HPM          4
#define OBJ_TYPE_HPM_EXTERNAL 5
#define OBJ_TYPE_POORMOTION   6
#define OBJ_TYPE_SUPPLEMENT   7
#define OBJ_TYPE_HPMNOPPMXL   8
#define OBJ_TYPE_HPMBADPPMXL  9
/* The following is from ucac5.h */
#pragma pack( 1)
typedef struct ucac5_star {
  long long srcid; /* Gaia source ID  */
  int32_t rag;   /* Gaia RA  at epoch 2015.0 (mas) */
  int32_t dcg;   /* Gaia Dec at epoch 2015.0 (mas) */
  int16_t erg;   /* Gaia DR1 position error RA  at epoch 2015.0 (0.1 mas) */
  int16_t edg;   /* Gaia DR1 position error Dec at epoch 2015.0 (0.1 mas) */     
  int8_t flg;    /* 1 = TGAS,  2 = other UCAC-Gaia star, 3 = other NOMAD */
  int8_t nu1;    /* number of images used for UCAC mean position */
  int16_t epu;   /* mean UCAC epoch (1/1000 yr after 1997.0) (myr) */
  int32_t  ira;  /* mean UCAC RA  at epu epoch on Gaia reference frame (mas) */
  int32_t idc;   /* mean UCAC Dec at epu epoch on Gaia reference frame (mas)  */
  int16_t pmur;  /* proper motion RA*cosDec (UCAC-Gaia) (0.1 mas/yr) */
  int16_t pmud;  /* proper motion Dec       (UCAC-Gaia) (0.1 mas/yr) */
  int16_t pmer;  /* formal error of UCAC-Gaia proper motion RA*cosDec (0.1 mas/yr) */
  int16_t pmed;  /* formal error of UCAC-Gaia proper motion Dec (0.1 mas/yr) */
  int16_t gmag;  /* Gaia DR1 G magnitude      (1/1000 mag) */
  int16_t umag;   /* mean UCAC model magnitude (1/1000 mag) */
  int16_t rmag;  /* NOMAD photographic R mag  {1/1000 mag) */
  int16_t jmag;  /* 2MASS J magnitude         (1/1000 mag) */
  int16_t hmag;  /* 2MASS H magnitude         (1/1000 mag) */
  int16_t kmag;  /* 2MASS K magnitude         (1/1000 mag) */ 
} UCAC5STAR,*PUCAC5STAR;
#pragma pack( )

typedef struct ucac5_convert {
  long long srcid; /* Gaia source ID  */
  double rag;   /* Gaia RA  at epoch 2015.0 (degrees) */
  double dcg;   /* Gaia Dec at epoch 2015.0 (degrees) */
  double erg;   /* Gaia DR1 position error RA  at epoch 2015.0 (arcsec) */
  double edg;   /* Gaia DR1 position error Dec at epoch 2015.0 (arcsec) */     
  int flg;    /* 1 = TGAS,  2 = other UCAC-Gaia star, 3 = other NOMAD */
  int nu1;    /* number of images used for UCAC mean position */
  double epu;   /* mean UCAC epoch  (yr) */
  double  ira;  /* mean UCAC RA  at epu epoch on Gaia reference frame (mas) */
  double idc;   /* mean UCAC Dec at epu epoch on Gaia reference frame (mas)  */
  double pmur;  /* proper motion RA*cosDec (UCAC-Gaia)  (mas/yr) */
  double pmud;  /* proper motion Dec       (UCAC-Gaia)  (mas/yr) */
  double pmer;  /* formal error of UCAC-Gaia proper motion RA*cosDec (mas/yr) */
  double pmed;  /* formal error of UCAC-Gaia proper motion Dec (mas/yr) */
  double gmag;  /* Gaia DR1 G magnitude      (1/1000 mag) */
  double umag;   /* mean UCAC model magnitude (1/1000 mag) */
  double rmag;  /* NOMAD photographic R mag  {1/1000 mag) */
  double jmag;  /* 2MASS J magnitude         (1/1000 mag) */
  double hmag;  /* 2MASS H magnitude         (1/1000 mag) */
  double kmag;  /* 2MASS K magnitude         (1/1000 mag) */ 
} UCAC5CONVERT,*PUCAC5CONVERT;

/*
  RAJ2000 deg 	DEJ2000 deg 	SrcIDgaia  	RAgaia deg 	e_ mas 	DEgaia deg 	e_ mas 	Org  	Nu  	EPucac yr 	pmRA mas/yr 	e_ (...) 	pmDE mas/yr 	e_ (...) 	Gmag mag 	f.mag mag 	Rmag mag 	Jmag mag 	Hmag mag 	Kmag mag
1	354.5066572	-89.9015083	6341070474361059328	354.5344572 	0.3	-89.9016417 	0.1	2 	0	      1998.745 	  10.6 	        30.8 	    -29.6	        30.8	    16.432	15.322	16.040	14.477	13.855	13.477
*pConvertRecord = {{srcid = 6341070474361059328, rag = 354.53445722222222, dcg = -89.901641666666663, erg = 0.00029999999999999997, edg = 0.0001, flg = 2, nu1 = 0, epu = 1998.7449999999999, ira = 354.5066572222222, idc = -89.901508333333339, pmur = 10.6, pmud = -29.600000000000001, pmer = 30.800000000000001, pmed = 30.800000000000001, gmag = 16.431999999999999, umag = 15.321999999999999, rmag = 16.039999999999999, jmag = 14.477, hmag = 13.855, kmag = 13.477}
 *pInput = {srcid = 6341070474361059328, rag = 1276324046, dcg = 3971321386, erg = 3, edg = 1, flg = 2 '\002', nu1 = 0 '\000', epu = 1745, ira = 1276223966, idc = 3971321866, pmur = 106, pmud = 65240, pmer = 308, pmed = 308, gmag = 16432, umag = 15322, rmag = 16040, jmag = 14477, hmag = 13855, kmag = 13477}

RAJ2000 deg  DEJ2000 deg 	SrcIDgaia  	        RAgaia deg 	  e_ mas 	DEgaia deg 	e_ mas 	Org  	Nu  	EPucac yr 	pmRA mas/yr 	e_DE (...) 	pmDE mas/yr 	e_RA (...) 	Gmag mag 	f.mag mag 	Rmag mag 	Jmag mag 	Hmag mag 	Kmag mag
263.0545906	-76.6292006	  5800984624723144704	263.0543697 	0.1	   -76.6292269 	0.2	    2 	  1	    1998.234 	 -11.0        	7.7 	      -5.6	        16.2	       15.972	    17.327	  15.740	    13.773	  13.093	  12.872
*pConvertRecord = {srcid = 5800984624723144704, rag = 263.05436972222225, dcg = -76.62922694444444, erg = 0.0001, edg = 0.00020000000000000001, flg = 2, nu1 = 1, epu = 1998.2339999999999, ira = 263.05459055555554, idc = -76.629200555555556, pmur = -11, pmud = -5.5999999999999996, pmer = 7.7000000000000002, pmed = 16.199999999999999, gmag = 15.972, umag = 17.327000000000002, rmag = 15.74, jmag = 13.773, hmag = 13.093, kmag = 12.872}
*pInput
$2 = {srcid = 5800984624723144704, rag = 946995731, dcg = 4019102079, erg = 1, edg = 2, flg = 2 '\002', nu1 = 1 '\001', epu = 1234, ira = 946996526, idc = 4019102174, pmur = 65426, pmud = 65480, pmer = 77, pmed = 162, gmag = 15972, umag = 17327, rmag = 15740, jmag = 13773, hmag = 13093, kmag = 12872}

idc = 1.0*deg + 1193.04665
deg = idc-1193.04665

min_totalPM 0.000000 max_totalPM 2867.153831 for srcid 540291183436392704 and REFNumber 74221024870596818
REFNumber           ra           dec      Stdmag  color             RaPM          DecPM    RaSigmaPM  DecSigmaPM
74221024870596818	359.999776	77.545924	16.931000	0.191000	0	0	0	-2867.000000	29.700001	18.500000	18.000000

RAJ2000 deg  DEJ2000 deg  SrcIDgaia           RAgaia deg    e_ mas  DEgaia deg  e_ mas  Org   Nu   EPucac yr  pmRA mas/yr  e_DE (...)  pmDE mas/yr  e_RA (...)  Gmag mag     f.mag mag  Rmag mag    Jmag mag    Hmag mag    Kmag mag
359.9997761  +77.5459244  540291183436392704  0.0002231     0.2     77.5460167  0.1      2     1   2003.805   -2867.0      18.5        29.7         18.0         16.931      16.581      16.740      15.499      15.080      15.048
359.999776    77.545924   540291183436392704  0.000223      0.2     77.546017   0.1      2     1   2003.805   -2867.0      18.5        29.7         18.0    gmag 16.931 umag 16.581 rmag 16.740 jmag 15.499 hmag 15.080 kmag 15.048

Barnard's star:
UCAC4: 64110809258	      117.563952	7.193521	16.339001	99.000000	4	0	0	207.000000	-1793.000000	8.000000	8.000000
UCAC5: 24569041783906881	117.563939	7.193560	16.233999	0.174000	0	0	0	209.600006	-1797.900024	14.000000	8.800000

RAJ2000 deg   DEJ2000 deg  SrcIDgaia               RAgaia deg    e_ mas  DEgaia deg  e_ mas  Org   Nu   EPucac yr  pmRA mas/yr  e_DE (...)  pmDE mas/yr  e_RA (...)  Gmag mag     f.mag mag  Rmag mag    Jmag mag    Hmag mag    Kmag mag
117.5639392  +07.1935600  3144837348340080768     117.5648064    1.0     7.1861778   0.5      2     1 2000.218     209.6         14.0      -1797.9  8.8               16.234      17.099      16.060      14.996      14.718      14.634
117.563939 idc 7.193560   3144837348340080768 rag 117.564806     1.0 dcg 7.186178    0.5      2     1 2000.218     209.6    pmer 14.0 pmud -1797.9  8.8          gmag 16.234 umag 17.099 rmag 16.060 jmag 14.996 hmag 14.718 kmag 14.634000




*/

#define UCAC5_ASCII_SIZE 280

         /* This flag suppresses stars that were matched with Tycho       */
         /* stars.  In some of my software,  stars are drawn from both    */
         /* UCAC-4 and Tycho.  By setting this flag in the output_format  */
         /* field,  I keep those stars from being plotted twice.          */
#define UCAC5_OMIT_TYCHO_STARS            0x1

         /* By default,  zero magnitudes and proper motions are written    */
         /* out as zeroes.  Setting this 'output_format' flag causes them  */
         /* to be written out as spaces.                                   */
#define UCAC5_WRITE_SPACES                0x2
#define UCAC5_FORTRAN_STYLE               0x4

         /* "Raw binary" means the data is written out in the same 78 byte */
         /* per star format in which it was read: no reformatting is done. */
#define UCAC5_RAW_BINARY                  0x8
/* Fields in the u4hpm.dat file */
#define MAX_U4HPM_RECORDS                 50
typedef struct u4hpm {
  int rnm;   /* unique star identifier (col. 51 main data) */
  int zn;    /* UCAC5 zone number for this star */
  int rnz;   /* unning record number along that zone for this star */
  int pmrc;  /* actual proper motion in RA*cos(dec) [0.1 mas/yr] */
  int pmd;   /* actual proper motion in Dec         [0.1 mas/yr] */
  int ra;    /* col. 1 of main data file (redundant) */
  int dec;   /* col. 2 of main data file (redundant) */
  int maga;  /* col. 4 of main data file (redundant) */
} U4HPM,*PU4HPM;

int UCAC5Compare(const void *first, const void *second) 
{
  double numberFirst = ((PUCAC5CONVERT)first)->idc;
  double numberSecond = ((PUCAC5CONVERT)second)->idc;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }
}

int main(int argc,char *argv[])
{
  char *argstr;
  char ucacpath[] = "/dasch/ucac5";
  char zonepath[MAX_INPUT_NAME];
  char hpmpath[] = "/dasch/ucac5/u4i/u4hpm.dat";
  FILE *fcat;
  FILE *hpmHandle;
  char inLine[MAX_INPUT_NAME];
  char *inBuffer;
  U4HPM u4hpm_table[MAX_U4HPM_RECORDS];
  PU4HPM pU4hpm = NULL;
  int hpmLines = 0;
  int hpmIndex;
  int hpmCount = 0;
	int hpmFlag;
  int lineLen;
  int nvals;
  char output_name[MAX_INPUT_NAME];
  File output_handle = NULL;
  int outputCount = 0;
  int zeromagCount = 0;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int linecounter = 0;
  char *charPtr;
  char cmdchar;
  int verbose = 0;
  char *dotPtr = NULL;
  int compareResult;
  int binNumber;
  int oldBinNumber = -1;
  int dec_bin;
  int dec_bin2;
  int oldDecBin = 0;
  int savedDecBin = -2;
  int ra_bin;
  int maxDecBin = 0;
  double maxRaError = 0;
  double maxDecError = 0;
  int maxWriteCount = 0;
  int totalBinCount = 0;

  File input_handle = NULL;
#if 0
  int input_nrecs = 0;
  int zone_nrecs = 0;
#endif 
  int input_index;
  int cur_index = 0;
  int cur_records = 0;
  int old_cur_records = 0;
  int zone_index;
  int temp_index;
  int total_nrecs = 0;
  int total_processed = 0;

  size_t nbr;
  int decint;
  double minOutputDec;
  double maxOutputDec;
  double minZoneDec;
  double maxZoneDec;

  
  PUCAC5STAR input_table;
  PUCAC5STAR pInput;
  GSCIMAGE gsc_record;
  PGSCIMAGE pGscImage = &gsc_record;
  PUCAC5CONVERT pConvertTable = NULL;
  PUCAC5CONVERT pConvertRecord = NULL;
  PUCAC5CONVERT pConvertRecord2 = NULL;

  UCAC5CONVERT max_convert_record;
  PUCAC5CONVERT pMaxConvertRecord = &max_convert_record;


 
  double tmpRaPM;
  double tmpDecPM;
  double tmpPM;
  double maxRaPM = 0;
  double maxDecPM = 0;
  double maxPM = 0;
  int countPM = 0;
  int tempmag;

  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int magcount[MAGBINS+1];
  int zone;
  int statResult;
  struct stat filestats;
  off_t filesize;
  int filerecords;
  long long min_id_number = -1;
  long long max_id_number = -1;
	int properMotionCount = 0;
	int pm_sigma;
	double PMSigma;
	FILE *properMotionHandle = NULL;

	int sigmaPMHist[MAXSIGMAPM+1];
	double totSigma;
	int sigmaPMIndex;
	int sigmaPMCumulative = 0;
	char REF[MAX_REF];
  long long max_totalPM_REFNumber;
	int refType;

  double min_epu;
  double max_epu;
  double min_ira;
  double max_ira;
  double min_idc;
  double max_idc;
  double min_pmur;
  double max_pmur;
  double min_pmud;
  double max_pmud;
  double min_pmer;
  double max_pmer;
  double min_pmed;
  double max_pmed;
  double min_gmag;
  double max_gmag;
  double min_umag;
  double max_umag;
  double min_rmag;
  double max_rmag;
  double min_jmag;
  double max_jmag;
  double min_hmag;
  double max_hmag;
  double min_kmag;
  double max_kmag;
  double totalPM;
  double min_totalPM;
  double max_totalPM;
  long long max_totalPM_srcid;
  double lastDec = -999.0;

  double plateepoch;
  double factor;
  
  long long REFNumber;
  long long HEALPIX_MASK;
  long long PROCCTR_MASK;
  long long SPARE_MASK;
  long long RUNSEQ_MASK;
  long long COMPONENT_MASK;
  long long ALL_MASK;

  long long max_HEALPIX = 0;
  long long min_HEALPIX = 0x7fffffffffffffffL;
  long long max_PROCCTR = 0;
  long long min_PROCCTR = 0x7fffffffffffffffL;
  long long max_SPARE = 0;
  long long min_SPARE = 0x7fffffffffffffffL;
  long long max_RUNSEQ = 0;
  long long min_RUNSEQ = 0x7fffffffffffffffL;
  long long max_COMPONENT = 0;
  long long min_COMPONENT = 0x7fffffffffffffffL;

	memset(sigmaPMHist,0,sizeof(sigmaPMHist));


#ifdef REFNUMBER_HACK
  printf("ERROR: REFNUMBER_HACK is set - fix GetREFNumber for Gaia\n");
  

#endif /* REFNUMBER_HACK */
#ifdef HIGH_PROPER_MOTION_CUTOFF
	printf("ERROR: HIGH_PROPER_MOTION_CUTOFF is set to %f\n",HIGH_PROPER_MOTION_CUTOFF);
#endif /* HIGH_PROPER_MOTION_CUTOFF */
#ifdef DUMP_PROPER_MOTION
	printf("ERROR: DUMP_PROPER_MOTION is set, writing to lospm.db\n");
	properMotionHandle = fopen("lospm.db","wt");
	if (properMotionHandle == NULL) {
		printf("ERROR writing lospm.db\n");
		exit(-1);
	}
	fprintf(properMotionHandle,"Stdmag\tclass\tRaSigmaPM\tDecSigmaPM\n");
	fprintf(properMotionHandle,"------\t-----\t---------\t----------\n");
#endif /* DUMP_PROPER_MOTION */

  /* Loop through the arguments */
  output_name[0] = 0;


  memset(magcount,0,sizeof(magcount));
  

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

        case 'o': /* output file name */
        case 'O':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(output_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
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



  if (output_name[0] == 0) {
    printf("ERROR: Output file name not specified\n");
    errorFlag = 1;
  } else {
    output_handle = fopen(output_name,"wt");
    if (output_handle == NULL) {
      errorFlag = 1;
      printf("Could not open output file %s\n",output_name);
    }
  }



  if (errorFlag) {
    printf("Usage: formatucac5 [-v][-s]\n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s provides statistics only and skips writing output files\n");

    return(-1);
  }

  printf("formatucac5 of %s %s Output Filename %s\n",
         __DATE__,__TIME__,output_name);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of UCAC5STAR is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(UCAC5STAR),sizeof(GSCIMAGE));
 
  HEALPIX_MASK = GAIA_HEALPIX_MASK << GAIA_HEALPIX_BIT;
  PROCCTR_MASK = GAIA_PROCCTR_MASK << GAIA_PROCCTR_BIT;
  SPARE_MASK = GAIA_SPARE_MASK << GAIA_SPARE_BIT;
  RUNSEQ_MASK = GAIA_RUNSEQ_MASK << GAIA_RUNSEQ_BIT;
  COMPONENT_MASK = GAIA_COMPONENT_MASK << GAIA_COMPONENT_BIT;
  ALL_MASK =  HEALPIX_MASK | PROCCTR_MASK |  SPARE_MASK | RUNSEQ_MASK |  COMPONENT_MASK;
  printf("HEALPIX_MASK:   0x%016llx\n",HEALPIX_MASK);
  printf("PROCCTR_MASK:   0x%016llx\n",PROCCTR_MASK);
  printf("SPARE_MASK:     0x%016llx\n",SPARE_MASK);
  printf("RUNSEQ_MASK:    0x%016llx\n",RUNSEQ_MASK);
  printf("COMPONENT_MASK: 0x%016llx\n",COMPONENT_MASK);
  printf("ALL_MASK:       0x%016llx\n",ALL_MASK);
  ALL_MASK = 0;
  if (HEALPIX_MASK != GAIA_HEALPIX) {
    printf("ERROR: HEALPIX_MASK is not GAIA_HEALPIX 0x%016llx\n",GAIA_HEALPIX);
  }
  if (PROCCTR_MASK != GAIA_PROCCTR) {
    printf("ERROR: PROCCTR_MASK is not GAIA_PROCCTR 0x%016llx\n",GAIA_PROCCTR);
  }
  if (SPARE_MASK != GAIA_SPARE) {
    printf("ERROR: SPARE_MASK is not GAIA_SPARE 0x%016llx\n",GAIA_SPARE);
  }
  if (RUNSEQ_MASK != GAIA_RUNSEQ) {
    printf("ERROR: RUNSEQ_MASK is not GAIA_RUNSEQ 0x%016llx\n",GAIA_RUNSEQ);
  }
  if (COMPONENT_MASK != GAIA_COMPONENT) {
    printf("ERROR: COMPONENT_MASK is not GAIA_COMPONENT 0x%016llx\n",GAIA_COMPONENT);
  }

 

  time(&startTime);
  input_table = (PUCAC5STAR)calloc(MAXRECORDS,sizeof(UCAC5STAR));
  if (input_table == NULL) {
    fprintf(stderr,"ERROR: failed to allocate %d bytes for %d input_table records\n",MAXRECORDS*sizeof(UCAC5STAR),MAXRECORDS);
  }
  pConvertTable = (PUCAC5CONVERT)calloc(2*MAXRECORDS,sizeof(UCAC5CONVERT));
  if (pConvertTable == NULL) {
    fprintf(stderr,"ERROR: failed to allocate %d bytes for %d pConvertTable records\n",2*MAXRECORDS*sizeof(UCAC5STAR),MAXRECORDS);
  }

  
  
  fprintf(output_handle,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tRaSigmaPM\tDecSigmaPM\n");
  fprintf(output_handle,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t---------\t----------\n");
  /* Now read in the input file */

  /* (zone = 1; zone <=  900; zone++) */
  for (zone = 1; zone <= 900; zone++) {
    /* Nominal zone boundaries (epoch 2015) */
    minZoneDec =  (0.2*zone)-90.2;
    maxZoneDec =  (0.2*(zone+1))-90.2;
    /* Overlapping zone boundaries for writing.  Put on a GSC bin index boundary */
    minOutputDec = (0.2*zone) -90.2 - (1./16.);
    maxOutputDec = (0.2*(zone+1)) -90.2 - (1./16.);
    if (zone == 900) {
      maxOutputDec = +99.9;
    }

    sprintf (zonepath, "%s/z%03d", ucacpath, zone);
    statResult = stat(zonepath,&filestats);
    if (statResult != 0) {
      /* No such file exists */
      fprintf (stderr,"UCACOPEN: failed to obtain file stats for UCAC file %s\n",zonepath);
      exit(-1);
    }
    filesize = filestats.st_size;
    filerecords = (filesize/sizeof(UCAC5STAR));
    if (filerecords > MAXRECORDS) {
      fprintf (stderr,"UCACOPEN: UCAC file %s has %d records while MAXRECORDS is %d\n",zonepath,filerecords,MAXRECORDS);
      exit(-1);
    }
    if (!(fcat = fopen (zonepath, "r"))) {
      fprintf (stderr,"UCACOPEN: UCAC file %s cannot be read\n",zonepath);
      exit(-1);
    }
    nbr = fread (input_table, sizeof(UCAC5STAR), filerecords, fcat);
    if (nbr != filerecords) {
      fprintf (stderr,"UCACOPEN: read %d of %d records for %s\n",nbr,filerecords,zonepath,filerecords);
      exit(-1);      
    }
    fclose(fcat);
    total_nrecs += filerecords;
 
    for (input_index = 0; input_index < filerecords; input_index++) {
			hpmFlag = 0;
      pInput = &input_table[input_index];
      pConvertRecord = &pConvertTable[cur_records];
      memset(pConvertRecord,0,sizeof(UCAC5CONVERT));
      total_processed++;
      if (verbose) {
        if (((total_processed+1) % 1000000) == 0) {
          time(&curTime);
          curTime -= startTime;
          printf("At record %7d, time %5d sec\n",total_processed,curTime);
        }
      }
      /* Convert the integer values to real */
      memset(pConvertRecord,0,sizeof(UCAC5CONVERT));

      ALL_MASK |= pInput->srcid;
      if (REFNumber < min_HEALPIX) {
        min_HEALPIX = REFNumber;
      } 
      if (REFNumber > max_HEALPIX) {
        max_HEALPIX = REFNumber;
      } 

      REFNumber = ((pInput->srcid >> GAIA_PROCCTR_BIT) & GAIA_PROCCTR_MASK);
      if (REFNumber < min_PROCCTR) {
        min_PROCCTR = REFNumber;
      } 
      if (REFNumber > max_PROCCTR) {
        max_PROCCTR = REFNumber;
      } 
      REFNumber = ((pInput->srcid >> GAIA_SPARE_BIT) & GAIA_SPARE_MASK);
      if (REFNumber < min_SPARE) {
        min_SPARE = REFNumber;
      } 
      if (REFNumber > max_SPARE) {
        max_SPARE = REFNumber;
      } 
      REFNumber = ((pInput->srcid >> GAIA_RUNSEQ_BIT) & GAIA_RUNSEQ_MASK);
      if (REFNumber < min_RUNSEQ) {
        min_RUNSEQ = REFNumber;
      } 
      if (REFNumber > max_RUNSEQ) {
        max_RUNSEQ = REFNumber;
      } 
      REFNumber = ((pInput->srcid >> GAIA_COMPONENT_BIT) & GAIA_COMPONENT_MASK);
      if (REFNumber < min_COMPONENT) {
        min_COMPONENT = REFNumber;
      } 
      if (REFNumber > max_COMPONENT) {
        max_COMPONENT = REFNumber;
      } 

#ifdef REFNUMBER_HACK
      REFNumber = ((pInput->srcid >> GAIA_HEALPIX_BIT) & GAIA_HEALPIX_MASK);
      
#else /* REFNUMBER_HACK */     


      if ((pInput->srcid % GAIA_MODULUS) != 0) {
        printf("ERROR: formatucac5 line %d invalid gaia modulus %d for source_id %lld\n",__LINE__,GAIA_MODULUS,pInput->srcid);
        exit(-1);
      }
#endif /* REFNUMBER_HACK */

      pConvertRecord->srcid = pInput->srcid;
      pConvertRecord->rag = pInput->rag/(1000.*3600.);
      pConvertRecord->dcg = pInput->dcg/(1000.*3600.);
      pConvertRecord->erg = pInput->erg/(10000.);
      pConvertRecord->edg = pInput->edg/(10000.);
      pConvertRecord->flg = pInput->flg;
      pConvertRecord->nu1 = pInput->nu1;
      pConvertRecord->epu = (pInput->epu/1000.)+1997.0;
      pConvertRecord->ira = pInput->ira/(1000.*3600.);
      pConvertRecord->idc = pInput->idc/(1000.*3600.);
      pConvertRecord->pmur = pInput->pmur/10.;
      pConvertRecord->pmud = pInput->pmud/10.;
      pConvertRecord->pmer = pInput->pmer/10.;
      pConvertRecord->pmed = pInput->pmed/10.;
      pConvertRecord->gmag = pInput->gmag/1000.;
      pConvertRecord->umag = pInput->umag/1000.;
      pConvertRecord->rmag = pInput->rmag/1000.;
      pConvertRecord->jmag = pInput->jmag/1000.;
      pConvertRecord->hmag = pInput->hmag/1000.;
      pConvertRecord->kmag = pInput->kmag/1000.;
   
      /* Now precess the ra and dec to J2000 */

      plateepoch = ((1.0*GSC_EQUINOX)-pConvertRecord->epu);
      pConvertRecord->idc =  pConvertRecord->idc + (pConvertRecord->pmud * plateepoch)/(3600.0*1000.0);
      factor = cos(DEGREES_TO_RAD*(pConvertRecord->idc));
      if (factor != 0) {
        pConvertRecord->ira = pConvertRecord->ira + ((pConvertRecord->pmur * plateepoch)/(3600.0*1000.0))/factor;
      }

      tempmag = pConvertRecord->gmag;
      if (tempmag < 0) {
        tempmag = 0;
      }
      if (tempmag >= MAGBINS) {
        tempmag = MAGBINS;
      }


      magcount[tempmag]++;
      cur_records++;
    }
    /* Now sort the records in increasing declination */
    qsort((void*)pConvertTable,cur_records,sizeof(UCAC5CONVERT),UCAC5Compare);
   
    for (cur_index = 0; cur_index < cur_records; cur_index++) {
      memset(pGscImage,0,sizeof(GSCIMAGE));
      pConvertRecord = &pConvertTable[cur_index];
      /* Now fill in our GSC record */

#ifdef REFNUMBER_HACK
      pGscImage->REFNumber = pInput->srcid;
#else /* REFNUMBER_HACK */     
      if ((pConvertRecord->srcid % GAIA_MODULUS) != 0) {
        printf("ERROR: formatucac5 line %d invalid gaia modulus %d for source_id %lld %s\n",__LINE__,GAIA_MODULUS,pConvertRecord->srcid,REF);
        exit(-1);
      }
      sprintf(REF,"GAIA%lld",pConvertRecord->srcid); /* UCAC5 used GAIA id's */
			if (GetREFNumber(REF,&pGscImage->REFNumber,&refType,0,0) != 0) {
        printf("ERROR: formatucac5: invalid REF %s line %d\n",REF,__LINE__);
				exit(-1);
			}
#endif /* REFNUMBER_HACK */
			pGscImage->ra =  pConvertRecord->ira;
      pGscImage->dec = pConvertRecord->idc;
      pGscImage->class = 0; /* Everything from UCAC5 is good for update */
      pGscImage->Stdmag = pConvertRecord->gmag;
      pGscImage->color = pConvertRecord->gmag-pConvertRecord->rmag;
      pGscImage->RaPM =  pConvertRecord->pmur;
      pGscImage->DecPM = pConvertRecord->pmud;
      pGscImage->RaSigmaPM = pConvertRecord->pmer;
      pGscImage->DecSigmaPM = pConvertRecord->pmed;
      /* Stop if we are going into the next zone */
      if (pGscImage->dec >= maxOutputDec) {
#if 0
        printf("line %d break at dec %f maxOutputDec %f\n",__LINE__,pGscImage->dec, maxOutputDec);
#endif
        break;
      }



#ifdef DUMP_PROPER_MOTION
      fprintf(properMotionHandle,"%.2f\t%d\t%.2f\t%.2f\n",pGscImage->Stdmag,pGscImage->class,pGscImage->RaSigmaPM,pGscImage->DecSigmaPM);
      if (++properMotionCount > 1000000) {
        printf("ERROR: properMotionCount exceeds %d\n",properMotionCount);
        exit(-1);
      }
#endif /* DUMP_PROPER_MOTION */

			totSigma = sqrt(sqr(pGscImage->DecSigmaPM)+sqr(pGscImage->RaSigmaPM));
			sigmaPMIndex = totSigma*10.0;
			if (sigmaPMIndex >= MAXSIGMAPM) {
				sigmaPMHist[MAXSIGMAPM]++;
			} else {
				sigmaPMHist[sigmaPMIndex]++;
			}

#if 0
			if (strcmp(REF,"U4000200137") == 0) {
				printf("%s\n",REF);
			}
#endif
      totalPM = sqrt(sqr(pConvertRecord->pmur)+sqr(pConvertRecord->pmud));
      if (outputCount == 0) {
        min_id_number = pConvertRecord->srcid;
        max_id_number = pConvertRecord->srcid;
      } else {
        if (pConvertRecord->srcid < min_id_number) {
          min_id_number = pConvertRecord->srcid;
        }
        if (pConvertRecord->srcid > max_id_number) {
          max_id_number = pConvertRecord->srcid;
        }
      }

      if (outputCount == 0) {
        min_epu = pConvertRecord->epu;
        max_epu = pConvertRecord->epu;
      } else {
        if (pConvertRecord->epu < min_epu) {
          min_epu = pConvertRecord->epu;
        }
        if (pConvertRecord->epu > max_epu) {
          max_epu = pConvertRecord->epu;
        }
      }        

      if (outputCount == 0) {
        min_ira = pConvertRecord->ira;
        max_ira = pConvertRecord->ira;
      } else {
        if (pConvertRecord->ira < min_ira) {
          min_ira = pConvertRecord->ira;
        }
        if (pConvertRecord->ira > max_ira) {
          max_ira = pConvertRecord->ira;
        }
      }        

      if (outputCount == 0) {
        min_idc = pConvertRecord->idc;
        max_idc = pConvertRecord->idc;
      } else {
        if (pConvertRecord->idc < min_idc) {
          min_idc = pConvertRecord->idc;
        }
        if (pConvertRecord->idc > max_idc) {
          max_idc = pConvertRecord->idc;
        }
      }        
      if (outputCount == 0) {
        min_pmur = pConvertRecord->pmur;
        max_pmur = pConvertRecord->pmur;
      } else {
        if (pConvertRecord->pmur < min_pmur) {
          min_pmur = pConvertRecord->pmur;
        }
        if (pConvertRecord->pmur > max_pmur) {
          max_pmur = pConvertRecord->pmur;
        }
      }        
      if (outputCount == 0) {
        min_pmud = pConvertRecord->pmud;
        max_pmud = pConvertRecord->pmud;
      } else {
        if (pConvertRecord->pmud < min_pmud) {
          min_pmud = pConvertRecord->pmud;
        }
        if (pConvertRecord->pmud > max_pmud) {
          max_pmud = pConvertRecord->pmud;
        }
      }        
      if (outputCount == 0) {
        min_pmer = pConvertRecord->pmer;
        max_pmer = pConvertRecord->pmer;
      } else {
        if (pConvertRecord->pmer < min_pmer) {
          min_pmer = pConvertRecord->pmer;
        }
        if (pConvertRecord->pmer > max_pmer) {
          max_pmer = pConvertRecord->pmer;
        }
      }        

      if (outputCount == 0) {
        min_pmed = pConvertRecord->pmed;
        max_pmed = pConvertRecord->pmed;
      } else {
        if (pConvertRecord->pmed < min_pmed) {
          min_pmed = pConvertRecord->pmed;
        }
        if (pConvertRecord->pmed > max_pmed) {
          max_pmed = pConvertRecord->pmed;
        }
      }        
      if (outputCount == 0) {
        min_gmag = pConvertRecord->gmag;
        max_gmag = pConvertRecord->gmag;
      } else {
        if (pConvertRecord->gmag < min_gmag) {
          min_gmag = pConvertRecord->gmag;
        }
        if (pConvertRecord->gmag > max_gmag) {
          max_gmag = pConvertRecord->gmag;
        }
      }        
      if (outputCount == 0) {
        min_umag = pConvertRecord->umag;
        max_umag = pConvertRecord->umag;
      } else {
        if (pConvertRecord->umag < min_umag) {
          min_umag = pConvertRecord->umag;
        }
        if (pConvertRecord->umag > max_umag) {
          max_umag = pConvertRecord->umag;
        }
      }        
      if (outputCount == 0) {
        min_rmag = pConvertRecord->rmag;
        max_rmag = pConvertRecord->rmag;
      } else {
        if (pConvertRecord->rmag < min_rmag) {
          min_rmag = pConvertRecord->rmag;
        }
        if (pConvertRecord->rmag > max_rmag) {
          max_rmag = pConvertRecord->rmag;
        }
      }        
      if (outputCount == 0) {
        min_jmag = pConvertRecord->jmag;
        max_jmag = pConvertRecord->jmag;
      } else {
        if (pConvertRecord->jmag < min_jmag) {
          min_jmag = pConvertRecord->jmag;
        }
        if (pConvertRecord->jmag > max_jmag) {
          max_jmag = pConvertRecord->jmag;
        }
      }        
      if (outputCount == 0) {
        min_hmag = pConvertRecord->hmag;
        max_hmag = pConvertRecord->hmag;
      } else {
        if (pConvertRecord->hmag < min_hmag) {
          min_hmag = pConvertRecord->hmag;
        }
        if (pConvertRecord->hmag > max_hmag) {
          max_hmag = pConvertRecord->hmag;
        }
      }        
      if (outputCount == 0) {
        min_kmag = pConvertRecord->kmag;
        max_kmag = pConvertRecord->kmag;
      } else {
        if (pConvertRecord->kmag < min_kmag) {
          min_kmag = pConvertRecord->kmag;
        }
        if (pConvertRecord->kmag > max_kmag) {
          max_kmag = pConvertRecord->kmag;
        }
      }
      if (outputCount == 0) {
        min_totalPM = totalPM;
        max_totalPM = totalPM;
        max_totalPM_srcid = pConvertRecord->srcid;
        max_totalPM_REFNumber = pGscImage->REFNumber;
      } else {
        if (totalPM < min_totalPM) {
          min_totalPM = totalPM;
        }
        if (totalPM > max_totalPM) {
          max_totalPM = totalPM;
          max_totalPM_srcid = pConvertRecord->srcid;
          max_totalPM_REFNumber = pGscImage->REFNumber;
        }
      }
      if ((pGscImage->DecPM > -1800.0) && (pGscImage->DecPM < -1790.0)) {
        memcpy(pMaxConvertRecord,pConvertRecord,sizeof(UCAC5CONVERT));
      }


			outputCount++;
      if (pGscImage->dec < lastDec) {
        printf("ERROR: dec is not sorted line %d\n",__LINE__);
        exit(-1);
      } 
#ifndef  REFNUMBER_HACK
			fprintf(output_handle,"%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",
							pGscImage->REFNumber,
							pGscImage->ra,
							pGscImage->dec,
							pGscImage->Stdmag,
							pGscImage->color,
							pGscImage->class,
							pGscImage->VFlag,
							pGscImage->MAGFlag,
							pGscImage->RaPM,
							pGscImage->DecPM,
							pGscImage->RaSigmaPM,
							pGscImage->DecSigmaPM);
#endif /* REFNUMBER_HACK */
		}
    /* If we have any left over records, copy them to the front of the table */
    for (temp_index = 0; temp_index < (cur_records-cur_index); temp_index++) {
      pConvertRecord = &pConvertTable[cur_index+temp_index];
      pConvertRecord2 = &pConvertTable[temp_index];
      memcpy(pConvertRecord2,pConvertRecord,sizeof(UCAC5CONVERT));
    }
    old_cur_records = cur_records;
    cur_records = cur_records-cur_index;
  

		if (verbose) {
			printf("Zone %3d has %8d records\n",zone,filerecords);
		}
	}
	if (input_table != NULL) {
		free(input_table);
	}
  if (pConvertTable != NULL) {
    free(pConvertTable);
  }

	Close(output_handle);
 
	time(&curTime);
	curTime -= startTime;

	printf("min_id_number %lld, max_id_number %lld hpmStars %d\n",min_id_number,max_id_number,hpmCount);
	printf("Execution Time: %d seconds for %d records %d processed  %d output\n",curTime,total_nrecs,total_processed,outputCount);

	for (tempmag = 0; tempmag <= MAGBINS; tempmag++) {
		printf("Mag %3d     ",tempmag);
		printf("%8d",magcount[tempmag]);
		printf("\n");
	}
	printf("Proper Motion Error (mas/yr),count,cumulative count:\n");
	for (sigmaPMIndex = 0; sigmaPMIndex <= MAXSIGMAPM; sigmaPMIndex++) {
		if (sigmaPMHist[sigmaPMIndex] > 0) {
			totSigma = (1.0*sigmaPMIndex)/10.0;
			sigmaPMCumulative += sigmaPMHist[sigmaPMIndex];
			printf("%4.1f %9d %9d\n",totSigma,sigmaPMHist[sigmaPMIndex],sigmaPMCumulative);
		}

	}
  printf("min_epu  %f max_epu  %f\n",min_epu,max_epu);
  printf("min_ira  %f max_ira  %f\n",min_ira,max_ira);
  printf("min_idc  %f max_idc  %f\n",min_idc,max_idc);
  printf("min_pmur %f max_pmur %f\n",min_pmur,max_pmur);
  printf("min_pmud %f max_pmud %f\n",min_pmud,max_pmud);
  printf("min_pmer %f max_pmer %f\n",min_pmer,max_pmer);
  printf("min_pmed %f max_pmed %f\n",min_pmed,max_pmed);
  printf("min_gmag %f max_gmag %f\n",min_gmag,max_gmag);
  printf("min_umag %f max_umag %f\n",min_umag,max_umag);
  printf("min_rmag %f max_rmag %f\n",min_rmag,max_rmag);
  printf("min_jmag %f max_jmag %f\n",min_jmag,max_jmag);
  printf("min_hmag %f max_hmag %f\n",min_hmag,max_hmag);
  printf("min_kmag %f max_kmag %f\n",min_kmag,max_kmag);
  printf("min_totalPM %f max_totalPM %f for srcid %lld and REFNumber %lld\n",min_totalPM,max_totalPM,max_totalPM_srcid,max_totalPM_REFNumber);

  printf("min_HEALPIX   0x%016llx max_HEALPIX    0x%016llx\n",min_HEALPIX,max_HEALPIX);
  printf("min_PROCCTR   0x%016llx max_PROCCTR    0x%016llx\n",min_PROCCTR,max_PROCCTR);
  printf("min_SPARE     0x%016llx max_SPARE      0x%016llx\n",min_SPARE,max_SPARE);
  printf("min_RUNSEQ    0x%016llx max_RUNSEQ     0x%016llx\n",min_RUNSEQ,max_RUNSEQ);
  printf("min_COMPONENT 0x%016llx max_COMPONENT  0x%016llx\n",min_COMPONENT,max_COMPONENT);
  printf("ALL_MASK      0x%016llx\n",ALL_MASK);  

  printf("srcid %lld rag %f dcg %f erg %f edg %f flg %d nu1 %d epu %f ira %f idc %f pmur %f pmud %f pmer %f pmed %f gmag %f umag %f rmag %f jmag %f hmag %f kmag %f\n",
         pMaxConvertRecord->srcid,
         pMaxConvertRecord->rag,
         pMaxConvertRecord->dcg,
         pMaxConvertRecord->erg,
         pMaxConvertRecord->edg,
         pMaxConvertRecord->flg,
         pMaxConvertRecord->nu1,
         pMaxConvertRecord->epu,
         pMaxConvertRecord->ira,
         pMaxConvertRecord->idc,
         pMaxConvertRecord->pmur,
         pMaxConvertRecord->pmud,
         pMaxConvertRecord->pmer,
         pMaxConvertRecord->pmed,
         pMaxConvertRecord->gmag,
         pMaxConvertRecord->umag,
         pMaxConvertRecord->rmag,
         pMaxConvertRecord->jmag,
         pMaxConvertRecord->hmag,
         pMaxConvertRecord->kmag);

	if (properMotionHandle != NULL) {
		fclose(properMotionHandle);
	}

	return(EXIT_SUCCESS);
}
