// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatgaia2.c
 *
 * Read and format the Gaia DR2 source list
 *
 * May 14, 2018 Edward J. Los - Adapted from formatgaia.c
 * May 18, 2018 Edward J. Los - Add tests of photometric accuracy and excess 
 * May 20, 2018 Edward J. Los - Accept stars with both phot_bp_mean_flux_over_error and phot_rp_mean_flux_over_error > 10
 *                              and phot_bp_rp_excess_factor < a + b ( phot_bp_mean_mag-phot_rp_mean_mag) **2 where a = 1.5 and b = 0.03
 * May 27, 2018 Edward J. Los - Add SEARCH_LOCATION_RA, SEARCH_LOCATION_DEC, and SEARCH_LOCATION_RADIUS
 * May 28, 2018 Edward J. Los - Correct proper motions for epoch 2000
 * May 30, 2018 Edward J. Los - Expand prefix to "GAIA2_"
 * Jun  2, 2018 Edward J. Los - Count each of the five contributing error factors.
 *                              Make the output file optional
 * Jul 19, 2018 Edward J. Los - Create a starbase file of gaia color relationships
 *
 * 
gcc -ggdb -O0   -D_FILE_OFFSET_BITS=64 -I/dasch/install/include  -I/usr/include/mysql -L/usr/lib${lib64}/mysql  formatgaia2.c -L/dasch/install/lib  pipelineutils.a  -lmysqlclient  -lwcs    -lcfitsio  -lm -o formatgaia2

*  On odyssey
gcc -ggdb -O0   -D_FILE_OFFSET_BITS=64 -I/dasch/install/include  -I/usr/include/mysql   formatgaia2.c  mysqlstub.o       -L/dasch/install/lib  -lcfitsio    -lm -o formatgaia2  pipelineutils.a  -lwcs

 * 
 *  

Testing:
find /home/scanner/Pipeline/catalogs/gaiadr2 -name "*.fits" > /home/scanner/junk/Gaia_test.list
formatgaia2 -v -i /home/scanner/junk/Gaia_test.list -o /dasch/Pipeline/catalogs/gaiadr2XXX.db -d /dasch/Pipeline/catalogs/gaiastarbase.db

head -n 3 /home/scanner/junk/Gaia_test.list > /home/scanner/junk/Gaia_testX.list 
formatgaia2 -v -i /home/scanner/junk/Gaia_testX.list


On Odyssey
setenv LD_LIBRARY_PATH "/n/home11/elos/dasch/install/lib"
/n/home11/elos/formatgaia2 -v -i /n/home11/elos/Gaia_DR2_source.list -o /n/home11/elos/gaiadr2.db

/n/home11/elos/formatgaia2 -v -i  /n/dasch8/scanner/backup/2018_05_03/Gaia_DR2_source.list -o /n/home11/elos/gaiadr2.db >& /n/home11/elos/formatgaia2.log

echo "formatgaia2 -v -i  /dasch/scanner/backup/2018_05_17/Gaia_DR2_source2.list -o /dasch/scanner/backup/2018_05_17/gaiadr2.db >& /dasch/scanner/backup/2018_05_17/formatgaia2_B.log" | at now

 *   source id field usage in Gaia DR2: see pipelineutils.h for bit definitions
min_HEALPIX   0x0000000000000000 max_HEALPIX    0x000000000bffffff
min_PROCCTR   0x0000000000000000 max_PROCCTR    0x0000000000000000
min_SPARE     0x0000000000000000 max_SPARE      0x0000000000000001
min_RUNSEQ    0x0000000000000000 max_RUNSEQ     0x00000000002ec6b0
min_COMPONENT 0x0000000000000000 max_COMPONENT  0x0000000000000000
ALL_MASK      0x7ffffff91fffff80
 *
 */   


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include "longnam.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#define MAX_INPUT_NAME 512

#define TblOff(p_type,field) \
        ((int) (((char *) (&(((p_type)NULL)->field))) - ((char *) NULL)))

#define FLUX_OVER_ERROR_MIN  20
#define MAGBINS  48  /* Number of magnitude bins */
#define MAGBINS_FACTOR 2
/* #define REFNUMBER_HACK 1 *//* temporary until GetREFNumber is fixed for Gaia */

#if 0
#define SEARCH_LOCATION_RA      226.258098
#define SEARCH_LOCATION_DEC   	42.075823	
#define SEARCH_LOCATION_RADIUS  0.01667  /* 1 arcmin expressed in degrees */
#endif
#define GAIADR2_EPOCH 2015.5
extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  
/* #define HIGH_PROPER_MOTION_CUTOFF 200.0 */ /* Proper motion cutoff in mas/yr */
/* #define DUMP_PROPER_MOTION 1  */ /* Print proper motion values */
/* #define DUMP_NEW_IDS 1 *//* Add Hipparchus and Tycho-2 ids. */
#define MAXSIGMAPM 650  /* Units of 0.1 mas/yr */
#define MAX_FILE_INDEX 16 /* Number of files in TGAS_SOURCE */

#define MAX_OUTPUT_TABLE 568200000  # expect  568134299 sources brighter than G: 19.0  x 56 bytes/GSCIMAGE = 1,750,749,672 bytes 


#define OBJ_TYPE_CLEAN        1
#define OBJ_TYPE_OVEREXPOSED  2
#define OBJ_TYPE_STREAK       3
#define OBJ_TYPE_HPM          4
#define OBJ_TYPE_HPM_EXTERNAL 5
#define OBJ_TYPE_POORMOTION   6
#define OBJ_TYPE_SUPPLEMENT   7
#define OBJ_TYPE_HPMNOPPMXL   8
#define OBJ_TYPE_HPMBADPPMXL  9
/* The following is from ucac4.h */
/* #pragma pack(1) */
typedef struct gaia_star {  /* Records are 412 bytes long */
 long long solution_id                     ; /* col  1 */
 char      designation[28]                 ; /* col  2 */
 long long source_id                       ; /* col  3 ** */
 long      random_index                    ; /* col  4 (long type not properly packed) */
 float     ref_epoch                       ; /* col  5 ** ! not J2015.5*/
 double    ra                              ; /* col  6 ** */
 float     ra_error                        ; /* col  7 */
 double    dec                             ; /* col  8 ** */
 float     dec_error                       ; /* col  9 */
 float     parallax                        ; /* col 10 */
 float     parallax_error                  ; /* col 11 */
 float     parallax_over_error             ; /* col 12 */
 float     pmra                            ; /* col 13 mas/year */
 float     pmra_error                      ; /* col 14 */
 float     pmdec                           ; /* col 15 mas/year */
 float     pmdec_error                     ; /* col 16 */
 float     ra_dec_corr                     ; /* col 17 */
 float     ra_parallax_corr                ; /* col 18 */
 float     ra_pmra_corr                    ; /* col 19 */
 float     ra_pmdec_corr                   ; /* col 20 */
 float     dec_parallax_corr               ; /* col 21 */
 float     dec_pmra_corr                   ; /* col 22 */
 float     dec_pmdec_corr                  ; /* col 23 */
 float     parallax_pmra_corr              ; /* col 24 */
 float     parallax_pmdec_corr             ; /* col 25 */
 float     pmra_pmdec_corr                 ; /* col 26 */
 short     astrometric_n_obs_al            ; /* col 27 */
 short     astrometric_n_obs_ac            ; /* col 28 */
 short     astrometric_n_good_obs_al       ; /* col 29 */
 short     astrometric_n_bad_obs_al        ; /* col 30 */
 float     astrometric_gof_al              ; /* col 31 */
 float     astrometric_chi2_al             ; /* col 32 */
 float     astrometric_excess_noise        ; /* col 33 */
 float     astrometric_excess_noise_sig    ; /* col 34 */
 short     astrometric_params_solved       ; /* col 35 */
 char      astrometric_primary_flag[5]     ; /* col 36 */
 float     astrometric_weight_al           ; /* col 37 */
 float     astrometric_pseudo_colour       ; /* col 38 */
 float     astrometric_pseudo_colour_error ; /* col 39 */
 float     mean_varpi_factor_al            ; /* col 40 */
 short     astrometric_matched_observations; /* col 41 */
 short     visibility_periods_used         ; /* col 42 */
 float     astrometric_sigma5d_max         ; /* col 43 */
 short     frame_rotator_object_type       ; /* col 44 */
 short     matched_observations            ; /* col 45 */
 char      duplicated_source[5]            ; /* col 46 */
 short     phot_g_n_obs                    ; /* col 47 */
 float     phot_g_mean_flux                ; /* col 48 */
 float     phot_g_mean_flux_error          ; /* col 49 */
 float     phot_g_mean_flux_over_error     ; /* col 50 */
 float     phot_g_mean_mag                 ; /* col 51 */
 short     phot_bp_n_obs                   ; /* col 52 */
 float     phot_bp_mean_flux               ; /* col 53 */
 float     phot_bp_mean_flux_error         ; /* col 54 */
 float     phot_bp_mean_flux_over_error    ; /* col 55 */
 float     phot_bp_mean_mag                ; /* col 56 */
 short     phot_rp_n_obs                   ; /* col 57 */
 float     phot_rp_mean_flux               ; /* col 58 */
 float     phot_rp_mean_flux_error         ; /* col 59 */
 float     phot_rp_mean_flux_over_error    ; /* col 60 */
 float     phot_rp_mean_mag                ; /* col 61 */
 float     phot_bp_rp_excess_factor        ; /* col 62 */
 short     phot_proc_mode                  ; /* col 63 */
 float     bp_rp                           ; /* col 64 */
 float     bp_g                            ; /* col 65 */
 float     g_rp                            ; /* col 66 */
 float     radial_velocity                 ; /* col 67 */
 float     radial_velocity_error           ; /* col 68 */
 short     rv_nb_transits                  ; /* col 69 */
 float     rv_template_teff                ; /* col 70 */
 float     rv_template_logg                ; /* col 71 */
 float     rv_template_fe_h                ; /* col 72 */
 char      phot_variable_flag[13]           ; /* col 73 */
 double    l                               ; /* col 74 */
 double    b                               ; /* col 75 */
 double    ecl_lon                         ; /* col 76 */
 double    ecl_lat                         ; /* col 77 */
 short     priam_flags                     ; /* col 78 */
 float     teff_val                        ; /* col 79 */
 float     teff_percentile_lower           ; /* col 80 */
 float     teff_percentile_upper           ; /* col 81 */
 float     a_g_val                         ; /* col 82 */
 float     a_g_percentile_lower            ; /* col 83 */
 float     a_g_percentile_upper            ; /* col 84 */
 float     e_bp_min_rp_val                 ; /* col 85 */
 float     e_bp_min_rp_percentile_lower    ; /* col 86 */
 float     e_bp_min_rp_percentile_upper    ; /* col 87 */
 short     flame_flags                     ; /* col 88 */
 float     radius_val                      ; /* col 89 */
 float     radius_percentile_lower         ; /* col 90 */
 float     radius_percentile_upper         ; /* col 91 */
 float     lum_val                         ; /* col 92 */
 float     lum_percentile_lower            ; /* col 93 */
 float     lum_percentile_upper            ; /* col 94 */
} GAIASTAR,*PGAIASTAR;
/* #pragma pack( ) */
 

#define GAIA_ASCII_SIZE 280

         /* This flag suppresses stars that were matched with Tycho       */
         /* stars.  In some of my software,  stars are drawn from both    */
         /* UCAC-4 and Tycho.  By setting this flag in the output_format  */
         /* field,  I keep those stars from being plotted twice.          */
#define GAIA_OMIT_TYCHO_STARS            0x1

         /* By default,  zero magnitudes and proper motions are written    */
         /* out as zeroes.  Setting this 'output_format' flag causes them  */
         /* to be written out as spaces.                                   */
#define GAIA_WRITE_SPACES                0x2
#define GAIA_FORTRAN_STYLE               0x4

         /* "Raw binary" means the data is written out in the same 78 byte */
         /* per star format in which it was read: no reformatting is done. */
#define GAIA_RAW_BINARY                  0x8
/* Fields in the u4hpm.dat file */
#define MAX_U4HPM_RECORDS                 50
typedef struct u4hpm {
  int rnm;   /* unique star identifier (col. 51 main data) */
  int zn;    /* GAIA zone number for this star */
  int rnz;   /* unning record number along that zone for this star */
  int pmrc;  /* actual proper motion in RA*cos(dec) [0.1 mas/yr] */
  int pmd;   /* actual proper motion in Dec         [0.1 mas/yr] */
  int ra;    /* col. 1 of main data file (redundant) */
  int dec;   /* col. 2 of main data file (redundant) */
  int maga;  /* col. 4 of main data file (redundant) */
} U4HPM,*PU4HPM;


typedef struct _datatype {
  int datatype;
  const char* vartype;
  const char* label;
} DATATYPE,*PDATATYPE;

DATATYPE datatype_table[] = {
  {0,  "unknown"  ,"UNKNOWN"},
  {1,  "unknown"  ,"TBIT"},
  {11, "unknown"  ,"TBYTE"},
  {12, "unknown"  ,"TSBYTE"},
  {14, "unknown"  ,"TLOGICAL"},
  {16, "char"     ,"TSTRING"},
  {20, "unknown"  ,"TUSHORT"},
  {21, "short"    ,"TSHORT"},
  {30, "unknown"  ,"TUINT"},
  {31, "unknown"  ,"TINT"},
  {40, "unknown"  ,"TULONG"},
  {41, "long"     ,"TLONG"},
  {41, "unknown"  ,"TINT32BIT"},
  {42, "float"    ,"TFLOAT"},
  {81, "long long","TLONGLONG"},
  {82, "double"   ,"TDOUBLE"},
  {83, "unknown"  ,"TCOMPLEX"},
  {163,"unknown"  ,"TDBLCOMPLEX"},
};
int datatype_size = sizeof(datatype_table)/sizeof(DATATYPE);

int GAIACompare(const void *first, const void *second) 
{
  double numberFirst = ((PGAIASTAR)first)->dec;
  double numberSecond = ((PGAIASTAR)second)->dec;
  if (numberFirst > numberSecond) {
    return(1);
  } else if (numberFirst < numberSecond) {
    return(-1);
  } else {
    return(0);
  }

}
void gaia_print_error(int columnnumber,int rownumber,char *fieldname,int status,int lineno) {
  printf("ERROR: read error in line %d status %d row %d column %d field %s\n",
         lineno,status,rownumber,columnnumber, fieldname);
}

int main(int argc,char *argv[])
{
  char *argstr;
  char zonepath[MAX_INPUT_NAME];
  FILE *fcat;
  char inLine[MAX_INPUT_NAME];
  char read_buffer[MAX_INPUT_NAME];
  char *inBuffer;
  U4HPM u4hpm_table[MAX_U4HPM_RECORDS];
  PU4HPM pU4hpm = NULL;
  int hpmLines = 0;
  int hpmIndex;
  int hpmCount = 0;
	int hpmFlag;
  int lineLen;
  char output_name[MAX_INPUT_NAME];
  FILE* output_handle = NULL;
  int outputCount = 0;
  char starbase_name[MAX_INPUT_NAME];
  FILE* starbase_handle = NULL;
  int starbaseCount = 0;
  char input_file_list[MAX_INPUT_NAME];
  FILE *input_file_handle = NULL;
  int input_file_count = 0;
  char input_file_index = 0;
  fitsfile *fptr;
  int inputCount = 0;
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
  int input_count;
  char input_name[MAX_INPUT_NAME];
  long input_nrecs = 0;
  int input_ncols = 0;
  int zone_nrecs = 0;
  char *arrayptr;
  char nullstr[1] = {0};
  int input_index;
  size_t nbr;
  int hdutype;
  int col_typecode;
  long col_repeat;
  long col_width;
  long total_width = 0;
  int index_ncols;
  int datatype_index;
  PDATATYPE pDatatype;
  GAIASTAR input_buffer;
  PGAIASTAR pInput = &input_buffer;
 
  GSCIMAGE gsc_record;
  PGSCIMAGE pGscImage = &gsc_record;

 
  double tmpRaPM;
  double tmpDecPM;
  double tmpPM;
  double maxRaPM = 0;
  double maxDecPM = 0;
  double maxPM = 0;
  double excess_limit;
  int countPM = 0;
  int tempmag;
  int status = 0;
  int skipOutput = 0;
  int magsBinNumber;
  int colorType;
  int g_magcount[MAGBINS+1];
  int bp_magcount[MAGBINS+1];
  int rp_magcount[MAGBINS+1];
  int good_magcount[MAGBINS+1];
  int bad_magcount[MAGBINS+1];
  int bpfluxerrorcount[MAGBINS+1];
  int rpfluxerrorcount[MAGBINS+1];
  int excessfactorcount[MAGBINS+1];
  int variablecount[MAGBINS+1];
  int duplicatecount[MAGBINS+1];

  int zone;
  int statResult;
  struct stat filestats;
  off_t filesize;
  int64_t min_id_number = -1;
  int64_t max_id_number = -1;
	int properMotionCount = 0;
	int pm_sigma;
	double PMSigma;
	FILE *properMotionHandle = NULL;
  char properm_name[MAX_INPUT_NAME];
  FILE *newIdHandle = NULL;
  char id_name[MAX_INPUT_NAME];
	int sigmaPMHist[MAXSIGMAPM+1];
	double totSigma;
	int sigmaPMIndex;
	int sigmaPMCumulative = 0;
	char REF[MAX_REF];
  char tempREF[2*MAX_REF];
	int refType;
  long nrows;
  long totnrows = 0;
  long firstchar;
  int ncols;
  int columnnumber = 0;
  char keystr[FLEN_KEYWORD];
  char keyname[FLEN_KEYWORD];
  
  char typevalue[80];
  char formvalue[80];
  char nullvalue[80];
  char ucdvalue[80];
  char unitvalue[80];
  char commentvalue[80];
  char comment[80];
  char formatname[80];
  char datatypename[80];
  char countname[80];
  int nvals;
  int anynull;
  int charsize;
  int selectedKeyIndex;
  int selectedKeys[] = {1,2,4,6,7,8,9,10,11,13,14,15,16,27,28,29,30,31,32,33,34,35,36,41,42,54,55,56,57};
  int numSelectedKeys = sizeof(selectedKeys)/sizeof(int);
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

  int variable_count = 0;
  int duplicated_source_count = 0;
  int astrometric_primary_count = 0;
  int g_count = 0;
  int bp_count = 0;
  int rp_count = 0;
  int g_bp_count = 0;
  int g_rp_count = 0;
  int bp_rp_count = 0;
  int file_bp_rp_count = 0;
  double min_dec = 9999;
  double max_dec = 9999;
  double min_ra = 9999;
  double max_ra = 9999;
  int gaia_modulus_error_count = 0;
  int calibrationStar = 0;
  int placeholderStar = 0;
  int output_nrecs = 0;
  int output_color_count = 0;
  int output_placeholder_count = 0;
  double factor;
  double dec2;
  double ra2;

#if 0
  {


    double ra2000; /* 45.570216  3h2.2"  where ra2008 = 45.675 or +0.013098 degrees/yr */
    double dec2000; /* 4.085550  4 deg 5 arcmin where dec2008 = 4.1167 or +0.0038937 degrees/yr */
    double ra2015;
    double dec2015;
    ra2000 = ra2015 = 45.03433035439128;   /* expect ra2000 = 44.841925 or  0.012827 degrees/yr */
    dec2000 = dec2015 =   0.23539164875137225; /* expect dec2000 = 0.176277 or 0.003941 degrees/yr */
    fk5prec(2015.0,2000.0,&ra2000,&dec2000);   
    printf("ra,dec 2000 %f %f ra,dec 2015 %f %f\n",ra2000,dec2000,ra2015,dec2015);
    exit(-1);
  }
#endif

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


#ifdef REFNUMBER_HACK
  printf("ERROR: REFNUMBER_HACK is set - fix GetREFNumber for Gaia\n");
#endif /* REFNUMBER_HACK */

	memset(sigmaPMHist,0,sizeof(sigmaPMHist));
#ifdef HIGH_PROPER_MOTION_CUTOFF
	printf("ERROR: HIGH_PROPER_MOTION_CUTOFF is set to %f\n",HIGH_PROPER_MOTION_CUTOFF);
#endif /* HIGH_PROPER_MOTION_CUTOFF */

  /* Loop through the arguments */
  output_name[0] = 0;
  starbase_name[0] = 0;
  input_file_list[0] = 0;

  memset(g_magcount,0,sizeof(g_magcount));
  memset(bp_magcount,0,sizeof(bp_magcount));
  memset(rp_magcount,0,sizeof(rp_magcount));
  memset(good_magcount,0,sizeof(good_magcount));
  memset(bad_magcount,0,sizeof(bad_magcount));
  memset(bpfluxerrorcount,0,sizeof(bpfluxerrorcount));
  memset(rpfluxerrorcount,0,sizeof(rpfluxerrorcount));
  memset(excessfactorcount,0,sizeof(excessfactorcount));
  memset(variablecount,0,sizeof(variablecount));
  memset(duplicatecount,0,sizeof(duplicatecount));

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

        case 'i': /* input file list name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(input_file_list,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-30) {
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
            strncpy(output_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 'd': /* color starbase file */
        case 'D':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(starbase_name,*++argv,MAX_INPUT_NAME-2);
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
  } else {
    output_handle = fopen(output_name,"wt");
    if (output_handle == NULL) {
      errorFlag = 1;
      printf("Could not open output file %s\n",output_name);
    }
  }

  if (starbase_name[0] != 0) {
    starbase_handle = fopen(starbase_name,"wt");
    if (starbase_handle == NULL) {
      errorFlag = 1;
      printf("Could not open starbase file %s\n",starbase_name);
    }
  }

  if (input_file_list[0] == 0) {
    printf("ERROR: Input directory not specified\n");
    errorFlag = 1;
  }


  if (errorFlag) {
    printf("Usage: formatgaia2 -i <input file list> -o <output file> [-v][-s]\n");
    printf("       where -v is the verbose flag\n");
    printf("             -s provides statistics only and skips writing output files\n");
    printf("             -d <database file> correlate all three GAIA colors\n");
    printf("             <input directory> input directory containing TgasSource_000-000-000.fits through TgasSource_000-000-015.fits\n");
    printf("             <output file> output starbase file name\n");

    return(-1);
  }

  printf("formatgaia2 of %s %s Output Filename %s\n",
         __DATE__,__TIME__,output_name);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of GAIASTAR is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(GAIASTAR),sizeof(GSCIMAGE));
 
  if (output_handle != NULL) {
    fprintf(output_handle,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tRaSigmaPM\tDecSigmaPM\n");
    fprintf(output_handle,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t---------\t----------\n");
  }
  if (starbase_handle != NULL) {
    fprintf(starbase_handle,"phot_g_mean_mag\tphot_bp_mean_mag\tphot_rp_mean_mag\n");
    fprintf(starbase_handle,"---------------\t----------------\t----------------\n");
  }

#ifdef DUMP_PROPER_MOTION
  strcpy(properm_name,output_name);
  charPtr = strrchr(properm_name,'/');
  if (charPtr == NULL) {
    properm_name[0] = 0;
  } else {
    charPtr += 1;
    *charPtr = 0;
  }
  strcat(properm_name,"lospm.db");
	printf("ERROR: DUMP_PROPER_MOTION is set, writing to %s\n",properm_name);
	properMotionHandle = fopen(properm_name,"wt");
	if (properMotionHandle == NULL) {
		printf("ERROR writing %s\n",properm_name);
		exit(-1);
	}
	fprintf(properMotionHandle,"Stdmag\tclass\tRaSigmaPM\tDecSigmaPM\n");
	fprintf(properMotionHandle,"------\t-----\t---------\t----------\n");
#endif /* DUMP_PROPER_MOTION */
#ifdef DUMP_NEW_IDS
  strcpy(id_name,output_name);
  charPtr = strrchr(id_name,'/');
  if (charPtr == NULL) {
    id_name[0] = 0;
  } else {
    charPtr += 1;
    *charPtr = 0;
  } 
  strcat(id_name,"losid.db");
	printf("ERROR: DUMP_NEW_IDS is set, writing to %s\n",id_name);
	newIdHandle = fopen(id_name,"wt");
	if (newIdHandle == NULL) {
		printf("ERROR writing %s\n",id_name);
		exit(-1);
	}
	fprintf(newIdHandle,"hip\ttycho2_id\tsource_id\tsource_idx\n");
	fprintf(newIdHandle,"---\t--------\t---------\t---------\n");
#endif /* DUMP_NEW_IDS */
  input_file_handle = fopen(input_file_list,"rt");
  if (input_file_handle == NULL) {
    printf("runpipeline failed to open file list %s line %d\n",input_file_list,__LINE__);
    return(0);
  }

  time(&startTime);
  while (1) {
    inBuffer = fgets(inLine,MAX_INPUT_NAME,input_file_handle);
    if (inBuffer == NULL) {
      break;
    }
    lineLen = strlen(inBuffer);
    /* Trim off the carriage return */
    if (inBuffer[lineLen-1] == 10) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    /* Trim off the line feed */
    if (inBuffer[lineLen-1] == 13) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    input_file_count++;
    strcpy(input_name,inBuffer);
    printf("File %5d is %s\n",input_file_count,input_name);
    fits_open_file(&fptr, input_name, READONLY, &status);
    if (status != 0) {
      printf("Error Opening %s\n",input_name);
      fits_report_error(stderr, status);
      errorFlag = 1;
    }

    fits_movabs_hdu (fptr, 2, &hdutype, &status);
    if ( status != 0 ) {
      fits_report_error(stderr, status);
      return(-1);
    }

    fits_get_num_rows (fptr,&nrows, &status);
    if ( status != 0 ) {
      fits_report_error(stderr, status);
      return(-1);
    }
    totnrows += nrows;
    fits_get_num_cols (fptr,&input_ncols, &status);
    if ( status != 0 ) {
      fits_report_error(stderr, status);
      return(-1);
    }


    printf("file %s has hdutype %d with %d records and %d columns\n",input_name,hdutype,nrows,input_ncols);
    if (input_file_count == 1) {
      total_width = 0;
      for (index_ncols = 1; index_ncols <= input_ncols; index_ncols++) {
        sprintf(keyname,"TTYPE%d",index_ncols);
        fits_read_key_str(fptr,keyname,keystr,NULL,&status);
        if ( status != 0 ) {
          fits_report_error(stderr, status);
          return(-1);
        }

        fits_get_coltype (fptr,index_ncols,&col_typecode, &col_repeat, &col_width,&status);
        if ( status != 0 ) {
          fits_report_error(stderr, status);
          return(-1);
        }
        for (datatype_index = 0; datatype_index < datatype_size; datatype_index++) {
          pDatatype = &datatype_table[datatype_index];
          if (pDatatype->datatype == col_typecode) {
            break;
          }
        }
        if (datatype_index == datatype_size) {
          datatype_index = 0;
          pDatatype = &datatype_table[datatype_index];
        }
          
        /* Column types
         *      TCOMM   ascii comment (optional)
         *      TFORM   format
         *        FORMAT '11A     '    ASCII string       
         *        FORMAT '13A     '           
         *        FORMAT 'D       '    double precision float
         *        FORMAT 'E       '    single precision float      
         *        FORMAT 'I       '    signed short      
         *        FORMAT 'J       '    signed 32 bit int   
         *        FORMAT 'K       '    64 bit long signed integer     
         *        FORMAT 'L       '    logicals (int for keywords and char for table cols)
         *      TNULL   blank value for integer
         *      TTYPE   name of the column
         *      TUCD    UCD for the column (optional)
         *       TUNIT   unit for the column (optional)
         */
#if 0
        printf("Column %2d has type %2d %10s repeat %2d width %2d offset %2d for %s = %s\n",index_ncols,col_typecode,pDatatype->label,col_repeat,col_width,total_width,keyname,keystr);
#endif
        total_width += col_width;

#if 0
        if (pDatatype->datatype == TSTRING) {
          sprintf(typevalue,"%s[%d]",keystr,col_width);
        } else {
          strcpy(typevalue,keystr);
        }
        local_strlwr(typevalue);
#if 1
        printf("fits_read_col(fptr,%s,%d,zone_nrecs+1,1,1,0,&pInput->%s,&anynull,&status);\n",pDatatype->label,index_ncols,typevalue);
        printf("if (status != 0) {\n");
        printf("  gaia_print_error(%d,zone_nrecs+1,\"%s\",status,__LINE__);\n",index_ncols,typevalue);
        printf("  exit(-1);\n");
        printf("}\n");
          
#endif

#if 0
        printf(" %-9s %-32s; /* col %2d */\n",pDatatype->vartype,typevalue,index_ncols);
#endif
#if 0
        /* The following becomes printf("field %-32s offset %2d\n","lum_percentile_upper",TblOff(PGAIASTAR,lum_percentile_upper)); */
        printf("XXX%sYYYYTblOff(PGAIASTAR,%s));\n",typevalue,typevalue);
#endif
#endif
      }
      printf("total_width %d sizeof GAIASTAR %d\n",total_width,sizeof(GAIASTAR));
    }
    min_dec = 9999;
    max_dec = 9999;
    min_ra = 9999;
    max_ra = 9999;
    file_bp_rp_count = 0;

    for (zone_nrecs = 0; zone_nrecs < nrows; zone_nrecs++) {
      memset(pInput,0,sizeof(GAIASTAR));
      memset(pGscImage,0,sizeof(GSCIMAGE));
      calibrationStar = 0;
      placeholderStar = 0;
      status = 0;
      firstchar = 1;

 
      fits_read_col(fptr,TLONGLONG,1,zone_nrecs+1,1,1,0,&pInput->solution_id,&anynull,&status);
      if (status != 0) {
        gaia_print_error(1,zone_nrecs+1,"solution_id",status,__LINE__);
        exit(-1);
      }
      arrayptr = pInput->designation;
      fits_read_col(fptr,TSTRING,2,zone_nrecs+1,1,1,0,&arrayptr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(2,zone_nrecs+1,"designation[28]",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TLONGLONG,3,zone_nrecs+1,1,1,0,&pInput->source_id,&anynull,&status);
      if (status != 0) {
        gaia_print_error(3,zone_nrecs+1,"source_id",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TLONG,4,zone_nrecs+1,1,1,0,&pInput->random_index,&anynull,&status);
      if (status != 0) {
        gaia_print_error(4,zone_nrecs+1,"random_index",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,5,zone_nrecs+1,1,1,0,&pInput->ref_epoch,&anynull,&status);
      if (status != 0) {
        gaia_print_error(5,zone_nrecs+1,"ref_epoch",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,6,zone_nrecs+1,1,1,0,&pInput->ra,&anynull,&status);
      if (status != 0) {
        gaia_print_error(6,zone_nrecs+1,"ra",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,7,zone_nrecs+1,1,1,0,&pInput->ra_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(7,zone_nrecs+1,"ra_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,8,zone_nrecs+1,1,1,0,&pInput->dec,&anynull,&status);
      if (status != 0) {
        gaia_print_error(8,zone_nrecs+1,"dec",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,9,zone_nrecs+1,1,1,0,&pInput->dec_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(9,zone_nrecs+1,"dec_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,10,zone_nrecs+1,1,1,0,&pInput->parallax,&anynull,&status);
      if (status != 0) {
        gaia_print_error(10,zone_nrecs+1,"parallax",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,11,zone_nrecs+1,1,1,0,&pInput->parallax_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(11,zone_nrecs+1,"parallax_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,12,zone_nrecs+1,1,1,0,&pInput->parallax_over_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(12,zone_nrecs+1,"parallax_over_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,13,zone_nrecs+1,1,1,0,&pInput->pmra,&anynull,&status);
      if (status != 0) {
        gaia_print_error(13,zone_nrecs+1,"pmra",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,14,zone_nrecs+1,1,1,0,&pInput->pmra_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(14,zone_nrecs+1,"pmra_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,15,zone_nrecs+1,1,1,0,&pInput->pmdec,&anynull,&status);
      if (status != 0) {
        gaia_print_error(15,zone_nrecs+1,"pmdec",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,16,zone_nrecs+1,1,1,0,&pInput->pmdec_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(16,zone_nrecs+1,"pmdec_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,17,zone_nrecs+1,1,1,0,&pInput->ra_dec_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(17,zone_nrecs+1,"ra_dec_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,18,zone_nrecs+1,1,1,0,&pInput->ra_parallax_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(18,zone_nrecs+1,"ra_parallax_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,19,zone_nrecs+1,1,1,0,&pInput->ra_pmra_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(19,zone_nrecs+1,"ra_pmra_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,20,zone_nrecs+1,1,1,0,&pInput->ra_pmdec_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(20,zone_nrecs+1,"ra_pmdec_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,21,zone_nrecs+1,1,1,0,&pInput->dec_parallax_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(21,zone_nrecs+1,"dec_parallax_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,22,zone_nrecs+1,1,1,0,&pInput->dec_pmra_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(22,zone_nrecs+1,"dec_pmra_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,23,zone_nrecs+1,1,1,0,&pInput->dec_pmdec_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(23,zone_nrecs+1,"dec_pmdec_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,24,zone_nrecs+1,1,1,0,&pInput->parallax_pmra_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(24,zone_nrecs+1,"parallax_pmra_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,25,zone_nrecs+1,1,1,0,&pInput->parallax_pmdec_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(25,zone_nrecs+1,"parallax_pmdec_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,26,zone_nrecs+1,1,1,0,&pInput->pmra_pmdec_corr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(26,zone_nrecs+1,"pmra_pmdec_corr",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,27,zone_nrecs+1,1,1,0,&pInput->astrometric_n_obs_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(27,zone_nrecs+1,"astrometric_n_obs_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,28,zone_nrecs+1,1,1,0,&pInput->astrometric_n_obs_ac,&anynull,&status);
      if (status != 0) {
        gaia_print_error(28,zone_nrecs+1,"astrometric_n_obs_ac",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,29,zone_nrecs+1,1,1,0,&pInput->astrometric_n_good_obs_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(29,zone_nrecs+1,"astrometric_n_good_obs_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,30,zone_nrecs+1,1,1,0,&pInput->astrometric_n_bad_obs_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(30,zone_nrecs+1,"astrometric_n_bad_obs_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,31,zone_nrecs+1,1,1,0,&pInput->astrometric_gof_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(31,zone_nrecs+1,"astrometric_gof_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,32,zone_nrecs+1,1,1,0,&pInput->astrometric_chi2_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(32,zone_nrecs+1,"astrometric_chi2_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,33,zone_nrecs+1,1,1,0,&pInput->astrometric_excess_noise,&anynull,&status);
      if (status != 0) {
        gaia_print_error(33,zone_nrecs+1,"astrometric_excess_noise",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,34,zone_nrecs+1,1,1,0,&pInput->astrometric_excess_noise_sig,&anynull,&status);
      if (status != 0) {
        gaia_print_error(34,zone_nrecs+1,"astrometric_excess_noise_sig",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,35,zone_nrecs+1,1,1,0,&pInput->astrometric_params_solved,&anynull,&status);
      if (status != 0) {
        gaia_print_error(35,zone_nrecs+1,"astrometric_params_solved",status,__LINE__);
        exit(-1);
      }
      arrayptr = pInput->astrometric_primary_flag;
      fits_read_col(fptr,TSTRING,36,zone_nrecs+1,1,1,0,&arrayptr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(36,zone_nrecs+1,"astrometric_primary_flag[5]",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,37,zone_nrecs+1,1,1,0,&pInput->astrometric_weight_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(37,zone_nrecs+1,"astrometric_weight_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,38,zone_nrecs+1,1,1,0,&pInput->astrometric_pseudo_colour,&anynull,&status);
      if (status != 0) {
        gaia_print_error(38,zone_nrecs+1,"astrometric_pseudo_colour",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,39,zone_nrecs+1,1,1,0,&pInput->astrometric_pseudo_colour_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(39,zone_nrecs+1,"astrometric_pseudo_colour_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,40,zone_nrecs+1,1,1,0,&pInput->mean_varpi_factor_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(40,zone_nrecs+1,"mean_varpi_factor_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,41,zone_nrecs+1,1,1,0,&pInput->astrometric_matched_observations,&anynull,&status);
      if (status != 0) {
        gaia_print_error(41,zone_nrecs+1,"astrometric_matched_observations",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,42,zone_nrecs+1,1,1,0,&pInput->visibility_periods_used,&anynull,&status);
      if (status != 0) {
        gaia_print_error(42,zone_nrecs+1,"visibility_periods_used",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,43,zone_nrecs+1,1,1,0,&pInput->astrometric_sigma5d_max,&anynull,&status);
      if (status != 0) {
        gaia_print_error(43,zone_nrecs+1,"astrometric_sigma5d_max",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,44,zone_nrecs+1,1,1,0,&pInput->frame_rotator_object_type,&anynull,&status);
      if (status != 0) {
        gaia_print_error(44,zone_nrecs+1,"frame_rotator_object_type",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,45,zone_nrecs+1,1,1,0,&pInput->matched_observations,&anynull,&status);
      if (status != 0) {
        gaia_print_error(45,zone_nrecs+1,"matched_observations",status,__LINE__);
        exit(-1);
      }
      arrayptr = pInput->duplicated_source;
      fits_read_col(fptr,TSTRING,46,zone_nrecs+1,1,1,0,&arrayptr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(46,zone_nrecs+1,"duplicated_source[5]",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,47,zone_nrecs+1,1,1,0,&pInput->phot_g_n_obs,&anynull,&status);
      if (status != 0) {
        gaia_print_error(47,zone_nrecs+1,"phot_g_n_obs",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,48,zone_nrecs+1,1,1,0,&pInput->phot_g_mean_flux,&anynull,&status);
      if (status != 0) {
        gaia_print_error(48,zone_nrecs+1,"phot_g_mean_flux",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,49,zone_nrecs+1,1,1,0,&pInput->phot_g_mean_flux_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(49,zone_nrecs+1,"phot_g_mean_flux_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,50,zone_nrecs+1,1,1,0,&pInput->phot_g_mean_flux_over_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(50,zone_nrecs+1,"phot_g_mean_flux_over_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,51,zone_nrecs+1,1,1,0,&pInput->phot_g_mean_mag,&anynull,&status);
      if (status != 0) {
        gaia_print_error(51,zone_nrecs+1,"phot_g_mean_mag",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,52,zone_nrecs+1,1,1,0,&pInput->phot_bp_n_obs,&anynull,&status);
      if (status != 0) {
        gaia_print_error(52,zone_nrecs+1,"phot_bp_n_obs",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,53,zone_nrecs+1,1,1,0,&pInput->phot_bp_mean_flux,&anynull,&status);
      if (status != 0) {
        gaia_print_error(53,zone_nrecs+1,"phot_bp_mean_flux",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,54,zone_nrecs+1,1,1,0,&pInput->phot_bp_mean_flux_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(54,zone_nrecs+1,"phot_bp_mean_flux_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,55,zone_nrecs+1,1,1,0,&pInput->phot_bp_mean_flux_over_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(55,zone_nrecs+1,"phot_bp_mean_flux_over_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,56,zone_nrecs+1,1,1,0,&pInput->phot_bp_mean_mag,&anynull,&status);
      if (status != 0) {
        gaia_print_error(56,zone_nrecs+1,"phot_bp_mean_mag",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,57,zone_nrecs+1,1,1,0,&pInput->phot_rp_n_obs,&anynull,&status);
      if (status != 0) {
        gaia_print_error(57,zone_nrecs+1,"phot_rp_n_obs",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,58,zone_nrecs+1,1,1,0,&pInput->phot_rp_mean_flux,&anynull,&status);
      if (status != 0) {
        gaia_print_error(58,zone_nrecs+1,"phot_rp_mean_flux",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,59,zone_nrecs+1,1,1,0,&pInput->phot_rp_mean_flux_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(59,zone_nrecs+1,"phot_rp_mean_flux_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,60,zone_nrecs+1,1,1,0,&pInput->phot_rp_mean_flux_over_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(60,zone_nrecs+1,"phot_rp_mean_flux_over_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,61,zone_nrecs+1,1,1,0,&pInput->phot_rp_mean_mag,&anynull,&status);
      if (status != 0) {
        gaia_print_error(61,zone_nrecs+1,"phot_rp_mean_mag",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,62,zone_nrecs+1,1,1,0,&pInput->phot_bp_rp_excess_factor,&anynull,&status);
      if (status != 0) {
        gaia_print_error(62,zone_nrecs+1,"phot_bp_rp_excess_factor",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,63,zone_nrecs+1,1,1,0,&pInput->phot_proc_mode,&anynull,&status);
      if (status != 0) {
        gaia_print_error(63,zone_nrecs+1,"phot_proc_mode",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,64,zone_nrecs+1,1,1,0,&pInput->bp_rp,&anynull,&status);
      if (status != 0) {
        gaia_print_error(64,zone_nrecs+1,"bp_rp",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,65,zone_nrecs+1,1,1,0,&pInput->bp_g,&anynull,&status);
      if (status != 0) {
        gaia_print_error(65,zone_nrecs+1,"bp_g",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,66,zone_nrecs+1,1,1,0,&pInput->g_rp,&anynull,&status);
      if (status != 0) {
        gaia_print_error(66,zone_nrecs+1,"g_rp",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,67,zone_nrecs+1,1,1,0,&pInput->radial_velocity,&anynull,&status);
      if (status != 0) {
        gaia_print_error(67,zone_nrecs+1,"radial_velocity",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,68,zone_nrecs+1,1,1,0,&pInput->radial_velocity_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(68,zone_nrecs+1,"radial_velocity_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,69,zone_nrecs+1,1,1,0,&pInput->rv_nb_transits,&anynull,&status);
      if (status != 0) {
        gaia_print_error(69,zone_nrecs+1,"rv_nb_transits",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,70,zone_nrecs+1,1,1,0,&pInput->rv_template_teff,&anynull,&status);
      if (status != 0) {
        gaia_print_error(70,zone_nrecs+1,"rv_template_teff",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,71,zone_nrecs+1,1,1,0,&pInput->rv_template_logg,&anynull,&status);
      if (status != 0) {
        gaia_print_error(71,zone_nrecs+1,"rv_template_logg",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,72,zone_nrecs+1,1,1,0,&pInput->rv_template_fe_h,&anynull,&status);
      if (status != 0) {
        gaia_print_error(72,zone_nrecs+1,"rv_template_fe_h",status,__LINE__);
        exit(-1);
      }
      arrayptr = pInput->phot_variable_flag;
      fits_read_col(fptr,TSTRING,73,zone_nrecs+1,1,1,0,&arrayptr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(73,zone_nrecs+1,"phot_variable_flag[13]",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,74,zone_nrecs+1,1,1,0,&pInput->l,&anynull,&status);
      if (status != 0) {
        gaia_print_error(74,zone_nrecs+1,"l",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,75,zone_nrecs+1,1,1,0,&pInput->b,&anynull,&status);
      if (status != 0) {
        gaia_print_error(75,zone_nrecs+1,"b",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,76,zone_nrecs+1,1,1,0,&pInput->ecl_lon,&anynull,&status);
      if (status != 0) {
        gaia_print_error(76,zone_nrecs+1,"ecl_lon",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,77,zone_nrecs+1,1,1,0,&pInput->ecl_lat,&anynull,&status);
      if (status != 0) {
        gaia_print_error(77,zone_nrecs+1,"ecl_lat",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,78,zone_nrecs+1,1,1,0,&pInput->priam_flags,&anynull,&status);
      if (status != 0) {
        gaia_print_error(78,zone_nrecs+1,"priam_flags",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,79,zone_nrecs+1,1,1,0,&pInput->teff_val,&anynull,&status);
      if (status != 0) {
        gaia_print_error(79,zone_nrecs+1,"teff_val",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,80,zone_nrecs+1,1,1,0,&pInput->teff_percentile_lower,&anynull,&status);
      if (status != 0) {
        gaia_print_error(80,zone_nrecs+1,"teff_percentile_lower",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,81,zone_nrecs+1,1,1,0,&pInput->teff_percentile_upper,&anynull,&status);
      if (status != 0) {
        gaia_print_error(81,zone_nrecs+1,"teff_percentile_upper",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,82,zone_nrecs+1,1,1,0,&pInput->a_g_val,&anynull,&status);
      if (status != 0) {
        gaia_print_error(82,zone_nrecs+1,"a_g_val",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,83,zone_nrecs+1,1,1,0,&pInput->a_g_percentile_lower,&anynull,&status);
      if (status != 0) {
        gaia_print_error(83,zone_nrecs+1,"a_g_percentile_lower",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,84,zone_nrecs+1,1,1,0,&pInput->a_g_percentile_upper,&anynull,&status);
      if (status != 0) {
        gaia_print_error(84,zone_nrecs+1,"a_g_percentile_upper",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,85,zone_nrecs+1,1,1,0,&pInput->e_bp_min_rp_val,&anynull,&status);
      if (status != 0) {
        gaia_print_error(85,zone_nrecs+1,"e_bp_min_rp_val",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,86,zone_nrecs+1,1,1,0,&pInput->e_bp_min_rp_percentile_lower,&anynull,&status);
      if (status != 0) {
        gaia_print_error(86,zone_nrecs+1,"e_bp_min_rp_percentile_lower",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,87,zone_nrecs+1,1,1,0,&pInput->e_bp_min_rp_percentile_upper,&anynull,&status);
      if (status != 0) {
        gaia_print_error(87,zone_nrecs+1,"e_bp_min_rp_percentile_upper",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,88,zone_nrecs+1,1,1,0,&pInput->flame_flags,&anynull,&status);
      if (status != 0) {
        gaia_print_error(88,zone_nrecs+1,"flame_flags",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,89,zone_nrecs+1,1,1,0,&pInput->radius_val,&anynull,&status);
      if (status != 0) {
        gaia_print_error(89,zone_nrecs+1,"radius_val",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,90,zone_nrecs+1,1,1,0,&pInput->radius_percentile_lower,&anynull,&status);
      if (status != 0) {
        gaia_print_error(90,zone_nrecs+1,"radius_percentile_lower",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,91,zone_nrecs+1,1,1,0,&pInput->radius_percentile_upper,&anynull,&status);
      if (status != 0) {
        gaia_print_error(91,zone_nrecs+1,"radius_percentile_upper",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,92,zone_nrecs+1,1,1,0,&pInput->lum_val,&anynull,&status);
      if (status != 0) {
        gaia_print_error(92,zone_nrecs+1,"lum_val",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,93,zone_nrecs+1,1,1,0,&pInput->lum_percentile_lower,&anynull,&status);
      if (status != 0) {
        gaia_print_error(93,zone_nrecs+1,"lum_percentile_lower",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,94,zone_nrecs+1,1,1,0,&pInput->lum_percentile_upper,&anynull,&status);
      if (status != 0) {
        gaia_print_error(94,zone_nrecs+1,"lum_percentile_upper",status,__LINE__);
        exit(-1);
      }

      ALL_MASK |= pInput->source_id;
      REFNumber =  ((pInput->source_id >> GAIA_HEALPIX_BIT) & GAIA_HEALPIX_MASK);
      if (REFNumber < min_HEALPIX) {
        min_HEALPIX = REFNumber;
      } 
      if (REFNumber > max_HEALPIX) {
        max_HEALPIX = REFNumber;
      } 

      REFNumber = ((pInput->source_id >> GAIA_PROCCTR_BIT) & GAIA_PROCCTR_MASK);
      if (REFNumber < min_PROCCTR) {
        min_PROCCTR = REFNumber;
      } 
      if (REFNumber > max_PROCCTR) {
        max_PROCCTR = REFNumber;
      } 
      REFNumber = ((pInput->source_id >> GAIA_SPARE_BIT) & GAIA_SPARE_MASK);
      if (REFNumber < min_SPARE) {
        min_SPARE = REFNumber;
      } 
      if (REFNumber > max_SPARE) {
        max_SPARE = REFNumber;
      } 
      REFNumber = ((pInput->source_id >> GAIA_RUNSEQ_BIT) & GAIA_RUNSEQ_MASK);
      if (REFNumber < min_RUNSEQ) {
        min_RUNSEQ = REFNumber;
      } 
      if (REFNumber > max_RUNSEQ) {
        max_RUNSEQ = REFNumber;
      } 
      REFNumber = ((pInput->source_id >> GAIA_COMPONENT_BIT) & GAIA_COMPONENT_MASK);
      if (REFNumber < min_COMPONENT) {
        min_COMPONENT = REFNumber;
      } 
      if (REFNumber > max_COMPONENT) {
        max_COMPONENT = REFNumber;
      } 

      if (pInput->phot_g_n_obs != 0) {
        g_count++;
        if ((pInput->phot_g_mean_mag > 0) &&
            (pInput->phot_g_mean_mag < MAX_STDMAG)) {
          placeholderStar = 1; /* Something will go into the output record */
          pGscImage->Stdmag = pInput->phot_g_mean_mag;
          pGscImage->color = 99.0;
        }
        if (pInput->phot_bp_n_obs != 0) {
          g_bp_count++;
        }
        if (pInput->phot_rp_n_obs != 0) {
          g_rp_count++;
        }
      }

      if ((pInput->phot_bp_mean_mag*MAGBINS_FACTOR) >= MAGBINS) {
        tempmag = MAGBINS;
      } else if (pInput->phot_bp_mean_mag < 0) {
        tempmag = 0;
      } else {
        tempmag = pInput->phot_bp_mean_mag*MAGBINS_FACTOR;
      }
      calibrationStar = 1; /* tentative calibration star*/
      if (strcmp(pInput->phot_variable_flag,"VARIABLE") == 0) {
        variable_count++;
        variablecount[tempmag]++;
        calibrationStar = 0;
      }
      if (strcmp(pInput->duplicated_source,"true") == 0) {
        duplicated_source_count++;
        duplicatecount[tempmag]++;
        calibrationStar = 0;
      }
      if (pInput->phot_bp_mean_flux_over_error <= FLUX_OVER_ERROR_MIN) {
        bpfluxerrorcount[tempmag]++;
      }
      if (pInput->phot_rp_mean_flux_over_error <= FLUX_OVER_ERROR_MIN) {
        rpfluxerrorcount[tempmag]++;
      }

      if (pInput->phot_bp_n_obs != 0) {
        bp_count++;
        if (pInput->phot_rp_n_obs != 0) {
          bp_rp_count++;
          if ((pInput->phot_bp_mean_mag > 0) &&
              (pInput->phot_bp_mean_mag < MAX_STDMAG)) {
            pGscImage->Stdmag = pInput->phot_bp_mean_mag;
            pGscImage->color = pInput->phot_bp_mean_mag - pInput->phot_rp_mean_mag;
          } else {
            calibrationStar = 0;
          }
          file_bp_rp_count++;
#if 0
          printf("phot_proc_mode %d phot_bp_rp_excess_factor %f\n",
                 pInput->phot_proc_mode,
                 pInput->phot_bp_rp_excess_factor);
#endif

        
          excess_limit = (1.5 + (0.03*sqr(pInput->phot_bp_mean_mag-pInput->phot_rp_mean_mag)));
          if (pInput->phot_bp_rp_excess_factor >= excess_limit) {
            excessfactorcount[tempmag]++;
          }

          if ((pInput->phot_bp_mean_flux_over_error > FLUX_OVER_ERROR_MIN) &&
              (pInput->phot_rp_mean_flux_over_error > FLUX_OVER_ERROR_MIN) &&
              (pInput->phot_bp_rp_excess_factor < excess_limit)) {
            good_magcount[tempmag]++;
          } else {
            bad_magcount[tempmag]++;
            calibrationStar = 0;
          }
        } else {
          calibrationStar = 0;
        }
      } else {
        calibrationStar = 0;
      }
      if (pInput->phot_rp_n_obs != 0) {
        rp_count++;
      }      
      if (strcmp(pInput->astrometric_primary_flag,"true") == 0) {
        astrometric_primary_count++;
      }
      

#if 0
      printf("phot_g_n_obs %3d phot_g_mean_mag %5.2f  phot_bp_n_obs %3d phot_bp_mean_mag %5.2f phot_rp_n_obs %3d phot_rp_mean_mag %5.2f\n",
             pInput->phot_g_n_obs,
             pInput->phot_g_mean_mag,
             pInput->phot_bp_n_obs,
             pInput->phot_bp_mean_mag,
             pInput->phot_rp_n_obs,
             pInput->phot_rp_mean_mag);
      printf("ref_epoch %f astrometric_primary_flag %5s duplicated_source %5s  phot_variable_flag %13s\n",
             pInput->astrometric_primary_flag,
             pInput->duplicated_source,
             pInput->phot_variable_flag);
      printf("ra %.4f dec %.4f,ra_error %f dec_error %f\n",
             pInput->ra,
             pInput->dec,
             pInput->ra_error,
             pInput->dec_error);
      printf("pmra %12.6f pmdec %12.6f pmra_error %8.6f pmdec_error %8.6f\n",
             pInput->pmra,
             pInput->pmdec,
             pInput->pmra_error,
             pInput->pmdec_error);
#endif

             
      if (min_dec == 9999) {
        min_dec = pInput->dec;
        max_dec = pInput->dec;
        min_ra = pInput->ra;
        max_ra = pInput->ra;
      } else {
        if (min_dec > pInput->dec) {
          min_dec = pInput->dec;
        } 
        if (max_dec < pInput->dec) {
          max_dec = pInput->dec;
        }
        if (min_ra > pInput->ra) {
          min_ra = pInput->ra;
        } 
        if (max_ra < pInput->ra) {
          max_ra = pInput->ra;
        }
         

      }


      if ((pInput->phot_g_mean_mag*MAGBINS_FACTOR) >= MAGBINS) {
        tempmag = MAGBINS;
      } else if (pInput->phot_g_mean_mag < 0) {
        tempmag = 0;
      } else {
        tempmag = pInput->phot_g_mean_mag*MAGBINS_FACTOR;
      }
      g_magcount[tempmag]++;

      if ((pInput->phot_bp_mean_mag*MAGBINS_FACTOR) >= MAGBINS) {
        tempmag = MAGBINS;
      } else if (pInput->phot_bp_mean_mag < 0) {
        tempmag = 0;
      } else {
        tempmag = pInput->phot_bp_mean_mag*MAGBINS_FACTOR;
      }
      bp_magcount[tempmag]++;

      if ((pInput->phot_rp_mean_mag*MAGBINS_FACTOR) >= MAGBINS) {
        tempmag = MAGBINS;
      } else if (pInput->phot_rp_mean_mag < 0) {
        tempmag = 0;
      } else {
        tempmag = pInput->phot_rp_mean_mag*MAGBINS_FACTOR;
      }
      rp_magcount[tempmag]++;



      if (min_id_number == -1) {
        min_id_number = pInput->source_id;
        max_id_number = pInput->source_id;
      } else {
        if (pInput->source_id < min_id_number) {
          min_id_number = pInput->source_id;
        }
        if (pInput->source_id > max_id_number) {
          max_id_number = pInput->source_id;
        }
      }
      
      if ((pInput->source_id % GAIA_MODULUS) != 0) {
        gaia_modulus_error_count++;
        calibrationStar = 0;
        placeholderStar = 0;
      }

      /* Reject any unwanted stars */
      if ((calibrationStar == 0) &&
          (placeholderStar == 0)) {
        continue;
      }
      output_nrecs++;
      if ((calibrationStar != 0) && (starbase_handle != NULL)) {
        starbaseCount++;
        fprintf(starbase_handle,"%f\t\%f\t\%f\n",
                pInput->phot_g_mean_mag,
                pInput->phot_bp_mean_mag,
                pInput->phot_rp_mean_mag);
      }

      /* Now put this in the GSC format */

      sprintf(REF,"GAIA2_%lld",pInput->source_id); /* A preceeding 'GAIA2_' flags gaia DR2 ids */
      if (strlen(REF) >= MAX_REF) {
        printf("ERROR: MAX_REF exceeded for pInput->source_id %lld\n",pInput->source_id);
        exit(-1);
      }
      if (GetREFNumber(REF,&pGscImage->REFNumber,&refType,0,0) != 0) {
        printf("ERROR: formatgaia: invalid REF %s\n",REF);
        exit(-1);
      }
      /* Now verify the REFNumber */
      if (GetREF(pGscImage->REFNumber,tempREF,0,0) != 0) {
        printf("ERROR: formatgaia: invalid REFNumber %lld for REF %s\n",pGscImage->REFNumber,REF);
        exit(-1);
      } else if (strcmp(REF,tempREF) != 0) {
        printf("ERROR: formatgaia: invalid REFNumber %lld for REF %s tempREF %s\n",pGscImage->REFNumber,REF,tempREF);
        exit(-1);
      }
      /* Precess to J2000 position */
      /* The GAIA Data Release 2 Documentation Release 1.0 states "At DR2 this reference epoch is always J2015.5 but in future releases this will
         be different and not necessarily the same for all sources".
         It was found that all RA and DEC values are actuall J2000, but with not corrected to J2000 */
      pGscImage->ra = pInput->ra;
      pGscImage->dec = pInput->dec;
#if 1
      if ((pGscImage->REFNumber == 810884079025127515L) ||
          (pGscImage->REFNumber == 812386570556416787L)) {
        char dstr[16];
        char rstr[16];
        dec2str(dstr,16,pGscImage->dec,2);
        ra2str(rstr,16,pGscImage->ra,2);
        printf("RA %f %s Dec %f %s for %lld in line %d\n",pGscImage->ra,rstr,pGscImage->dec,dstr,pGscImage->REFNumber,__LINE__);
        
      }
#endif          
#if 0
      fk5prec(GAIADR2_EPOCH,2000.0,&pGscImage->ra,&pGscImage->dec);
#endif
#if 1
      if ((pGscImage->REFNumber == 810884079025127515L) ||
          (pGscImage->REFNumber == 812386570556416787L)) {
        char dstr[16];
        char rstr[16];
        dec2str(dstr,16,pGscImage->dec,2);
        ra2str(rstr,16,pGscImage->ra,2);
        printf("RA %f %s Dec %f %s for %lld in line %d\n",pGscImage->ra,rstr,pGscImage->dec,dstr,pGscImage->REFNumber,__LINE__);
        
      }
#endif          

      pGscImage->class = 0;
      if (calibrationStar == 0) {
        pGscImage->VFlag = 1;
        output_placeholder_count++;
      } else {
        output_color_count++;
      }

      pGscImage->RaPM = pInput->pmra;
      pGscImage->RaSigmaPM = pInput->pmra_error;
      pGscImage->DecPM = pInput->pmdec;
      pGscImage->DecSigmaPM = pInput->pmdec_error;

      totSigma = sqrt(sqr(pGscImage->DecSigmaPM)+sqr(pGscImage->RaSigmaPM));
      sigmaPMIndex = totSigma*10.0;
      if (sigmaPMIndex >= MAXSIGMAPM) {
        sigmaPMHist[MAXSIGMAPM]++;
      } else {
        sigmaPMHist[sigmaPMIndex]++;
      }
      /* Now correct for proper motions to year 2000 */

#ifdef SEARCH_LOCATION_RA
      {
        double dist;
        factor = cos(DEGREES_TO_RAD*(pGscImage->dec));
        if ((fabs(pGscImage->ra - SEARCH_LOCATION_RA)/factor < SEARCH_LOCATION_RADIUS) &&
            (fabs(pGscImage->dec - SEARCH_LOCATION_DEC) < SEARCH_LOCATION_RADIUS)) {
          dist = sqrt(sqr((pGscImage->ra - SEARCH_LOCATION_RA)/factor) + sqr(pGscImage->dec - SEARCH_LOCATION_DEC));
          if (dist < SEARCH_LOCATION_RADIUS) {
            printf("Dist is (1) %f arcsec for %lld %f %f %f %f %d %d %d %f %f %f %f\n",
                   3600.*dist,
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
          }
        }
      }
#endif /* SEARCH_LOCATION_RA */

      dec2 = pGscImage->dec + (pGscImage->DecPM * (GSC_EQUINOX-GAIADR2_EPOCH))/(3600.0*1000.0);
      factor = cos(DEGREES_TO_RAD*(dec2));
      if (factor == 0) {
        ra2 = pGscImage->ra;
      } else {
        ra2  = pGscImage->ra + ((pGscImage->RaPM * (GSC_EQUINOX-GAIADR2_EPOCH))/(3600.0*1000.0*factor));
      }
      pGscImage->dec = dec2;
      pGscImage->ra = ra2;
#ifdef SEARCH_LOCATION_RA
      {
        double dist;
        factor = cos(DEGREES_TO_RAD*(pGscImage->dec));
        if ((fabs(pGscImage->ra - SEARCH_LOCATION_RA)/factor < SEARCH_LOCATION_RADIUS) &&
            (fabs(pGscImage->dec - SEARCH_LOCATION_DEC) < SEARCH_LOCATION_RADIUS)) {
          dist = sqrt(sqr((pGscImage->ra - SEARCH_LOCATION_RA)/factor) + sqr(pGscImage->dec - SEARCH_LOCATION_DEC));
          if (dist < SEARCH_LOCATION_RADIUS) {
            printf("Dist is (2) %f arcsec for %lld %f %f %f %f %d %d %d %f %f %f %f\n",
                   3600.*dist,
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
          }
        }
      }
#endif /* SEARCH_LOCATION_RA */

      if (output_handle != NULL) {
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

      }
    }



    status = 0;
    time(&curTime);
    curTime -= startTime;
    printf("File %5d file_bp_rp_count %8d min_dec %f max_dec %f min_ra %f max_ra %f\n",
           input_file_count,
           file_bp_rp_count,
           min_dec,
           max_dec,
           min_ra,
           max_ra);
    printf("Finished with record %d for file %s at %d seconds\n",input_nrecs,input_name,curTime);
    if ( fits_close_file(fptr, &status) ) {
      fits_report_error(stderr, status);
      return(-1);
    }
    

  }

  printf("totnrows %d\n",totnrows);
#if 0
  /* Now sort the records in increasing declination */
  qsort((void*)input_table,input_nrecs,sizeof(GAIASTAR),GAIACompare);
  for (input_count = 0; input_count < input_nrecs; input_count++) {
    pInput = &input_table[input_count];
 
    memset(pGscImage,0,sizeof(GSCIMAGE));
    if (verbose) {
      if (((input_count+1) % 1000000) == 0) {
        time(&curTime);
        curTime -= startTime;
        printf("At record %7d, time %5d sec\n",input_nrecs,curTime);
      }
    }
        


    pGscImage->class = 0;
    /* NOTE: Stdmag and color is handled above */
    


#ifdef DUMP_PROPER_MOTION
    fprintf(properMotionHandle,"%.2f\t%d\t%.2f\t%.2f\n",pGscImage->Stdmag,pGscImage->class,pGscImage->RaSigmaPM,pGscImage->DecSigmaPM);
#endif /* DUMP_PROPER_MOTION */
#ifdef DUMP_NEW_IDS
    fprintf(newIdHandle,"%d\t%s\t%lld\t0x%016llx\n",pInput->hip,pInput->tycho2_id,pInput->source_id,pInput->source_id);

#endif /* DUMP_NEW_IDS */

  
    totSigma = sqrt(sqr(pGscImage->DecSigmaPM)+sqr(pGscImage->RaSigmaPM));
    sigmaPMIndex = totSigma*10.0;
    if (sigmaPMIndex >= MAXSIGMAPM) {
      sigmaPMHist[MAXSIGMAPM]++;
    } else {
      sigmaPMHist[sigmaPMIndex]++;
    }

#ifdef HIGH_PROPER_MOTION_CUTOFF
#if 0
    if (hpmFlag == 0) {
      continue;
    }
#else
    if (((pGscImage->RaPM*pGscImage->RaPM) + (pGscImage->DecPM*pGscImage->DecPM)) < 
        ((HIGH_PROPER_MOTION_CUTOFF*HIGH_PROPER_MOTION_CUTOFF))) {
      continue;
    }
#endif
			
#endif /* HIGH_PROPER_MOTION_CUTOFF */
#if 0
    if (strcmp(REF,"G1000200137") == 0) {
      printf("%s\n",REF);
    }
#endif


    outputCount++;
    if (output_handle != NULL ) {
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

    }
  }
#endif

  if (output_handle != NULL) {
    fclose(output_handle);
  }
  if (starbase_handle != NULL) {
    fclose(starbase_handle);
  }
  
  time(&curTime);
  curTime -= startTime;

  printf("min_id_number %lld, max_id_number %lld hpmStars %d\n",min_id_number,max_id_number,hpmCount);
  printf("g, bp, and rp magnitude counts good and bad bad_bp_flux bad_rp_flux badexcessfactor variable duplicate\n");
  for (tempmag = 0; tempmag <= MAGBINS; tempmag++) {
    printf("Mag %4.1f     ",(1.0*tempmag)/(1.0*MAGBINS_FACTOR));
    printf("%11d",g_magcount[tempmag]);
    printf("%11d",bp_magcount[tempmag]);
    printf("%11d",rp_magcount[tempmag]);
    printf("%11d",good_magcount[tempmag]);
    printf("%11d",bad_magcount[tempmag]);
    printf("%11d",bpfluxerrorcount[tempmag]);
    printf("%11d",rpfluxerrorcount[tempmag]);
    printf("%11d",excessfactorcount[tempmag]);
    printf("%11d",variablecount[tempmag]);
    printf("%11d",duplicatecount[tempmag]);
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

  printf("min_HEALPIX   0x%016llx max_HEALPIX    0x%016llx\n",min_HEALPIX,max_HEALPIX);
  printf("min_PROCCTR   0x%016llx max_PROCCTR    0x%016llx\n",min_PROCCTR,max_PROCCTR);
  printf("min_SPARE     0x%016llx max_SPARE      0x%016llx\n",min_SPARE,max_SPARE);
  printf("min_RUNSEQ    0x%016llx max_RUNSEQ     0x%016llx\n",min_RUNSEQ,max_RUNSEQ);
  printf("min_COMPONENT 0x%016llx max_COMPONENT  0x%016llx\n",min_COMPONENT,max_COMPONENT);
  printf("ALL_MASK      0x%016llx\n",ALL_MASK);  
  printf("variable_count %10d duplicated_source_count %10d astrometric_primary_count %10d\n",
         variable_count,
         duplicated_source_count,
         astrometric_primary_count);

  printf("g_count %10d bp_count %10d rp_count %10d g_bp_count %10d g_rp_count %10d bp_rp_count %10d\n",
         g_count,
         bp_count,
         rp_count,
         g_bp_count,
         g_rp_count,
         bp_rp_count);
  printf("gaia_modulus_error_count %d\n",gaia_modulus_error_count);
  printf("output_nrecs %d output_color_count %d output_placeholder_count %d starbaseCount %d at %lld seconds\n",output_nrecs,output_color_count,output_placeholder_count,starbaseCount,curTime);

  if (properMotionHandle != NULL) {
    fclose(properMotionHandle);
  }
  if (newIdHandle != NULL) {
    fclose(newIdHandle);
  }
  if (input_file_handle != NULL) {
    fclose(input_file_handle);
  }
  return(EXIT_SUCCESS);
}
