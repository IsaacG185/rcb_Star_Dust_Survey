// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatgaia.c
 *
 * Read the Tycho-Gaia Astrometric Solution catalog and put it into the gsc232bin.db format
 *
 * NOTE: maga (UCAC aperture magnitude) is given for Stdmag, no color is given.
 *       the obj_type is copied into the gsc class field.  Only zero obj_types have good astrometry.
 *     
 * NOTE: the highest proper motions are RaPM -4.41 to +6.77 arcsec/year and DecPM -5.82 to +10.33 arcsec/yr.
 *
 * cc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm formatgaia.c pipelineutils.a  -ltable -lutil -lwcs /dasch/install/lib/libcfitsio.a -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -o formatgaia
 * 
 *  formatgaia -v -i /home/scanner/junk/gaia/TgasSource/ -o /dasch/Pipeline/catalogs/gaiatgasdr1.db

date
formatgaia -v -i /home/scanner/junk/gaia/TgasSource/ -o /dasch/Pipeline/catalogs/gaiatgasdr1.db
date
formatgsc  -p -q apass /dasch/Pipeline/catalogs/gaiatgasdr1.db
date
matchcatalogs  -c t -t 1.0 -m -w -f  /dasch/Pipeline/catalogs/gaiatgasdr1.dat -s   /dasch/Pipeline/catalogs/ucac4.dat
date

On Odyssey
formatgaia -v -i /n/dasch12/gaia/TgasSource/ -o /dasch/Pipeline/catalogs/gaiatgasdr1.db


 *
 * Nov 21, 2016 Edward J. Los - Initial version
 * Nov 25, 2016 Edward J. Los - Precess ra and dec to J2000.  (TODO: need also to precess the proper motion near the pole!).
 *                              Use GAIA as prefix to source_id
 * May 30, 2018 Edward J. Los - expand prefix to "GAIA_"
 */   


#include <math.h>
#include <time.h>
#include "fitsio.h"
#include "table.h"
#include "longnam.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#define MAX_INPUT_NAME 512

#define MAGBINS  21   /* Number of magnitude bins */
#define MAXRECORDS 2058000 /* Count for DR1 TgasSource */
/* #define REFNUMBER_HACK 1 *//* temporary until GetREFNumber is fixed for Gaia */


extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

/* #define LOS_DEBUG 1 */  
/* #define HIGH_PROPER_MOTION_CUTOFF 200.0 */ /* Proper motion cutoff in mas/yr */
/* #define DUMP_PROPER_MOTION 1  */ /* Print proper motion values */
/* #define DUMP_NEW_IDS 1 *//* Add Hipparchus and Tycho-2 ids. */
#define MAXSIGMAPM 650  /* Units of 0.1 mas/yr */
#define MAX_FILE_INDEX 16 /* Number of files in TGAS_SOURCE */

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
/* #pragma pack( 1) */
typedef struct gaia_star {
    int32_t                      hip     ; /* Key  1 -2147483648  'meta.id.cross'  */
       char                 tycho2_id[12]; /* Key  2   'meta.id.cross'  */
    int64_t                     source_id; /* Key  4   'meta.id;meta.main'  */
     double                     ref_epoch; /* Key  6  'Time[Julian Years]' 'meta.ref;time.epoch'  */
     double                      ra      ; /* Key  7  'Angle[deg]' 'pos.eq.ra;meta.main'  */
     double                      ra_error; /* Key  8  'Angle[mas]' 'stat.error;pos.eq.ra'  */
     double                      dec     ; /* Key  9  'Angle[deg]' 'pos.eq.dec;meta.main'  */
     double                     dec_error; /* Key 10  'Angle[mas]' 'stat.error;pos.eq.dec'  */
     double                      parallax; /* Key 11  'Angle[mas]' 'pos.parallax'  */
     double                      pmra    ; /* Key 13  'Angular Velocity[mas/year]' 'pos.pm;pos.eq.ra'  */
     double                    pmra_error; /* Key 14  'Angular Velocity[mas/year]' 'stat.error;pos.pm;pos.eq.ra'  */
     double                      pmdec   ; /* Key 15  'Angular Velocity[mas/year]' 'pos.pm;pos.eq.dec'  */
     double                   pmdec_error; /* Key 16  'Angular Velocity[mas/year]' 'stat.error;pos.pm;pos.eq.dec'  */
    int32_t          astrometric_n_obs_al; /* Key 27   'meta.number'  */
    int32_t          astrometric_n_obs_ac; /* Key 28   'meta.number'  */
    int32_t     astrometric_n_good_obs_al; /* Key 29   'meta.number'  */
    int32_t     astrometric_n_good_obs_ac; /* Key 30   'meta.number'  */
    int32_t      astrometric_n_bad_obs_al; /* Key 31   'meta.number'  */
    int32_t      astrometric_n_bad_obs_ac; /* Key 32   'meta.number'  */
      float           astrometric_delta_q; /* Key 33   'stat.value'  */
     double      astrometric_excess_noise; /* Key 34  'Angle[mas]' 'stat.value'  */
     double  astrometric_excess_noise_sig; /* Key 35   'stat.value'  */
       char      astrometric_primary_flag; /* Key 36   'meta.code'  */
    int16_t          matched_observations; /* Key 41   'meta.number'  */
       char             duplicated_source; /* Key 42     */
     double               phot_g_mean_mag; /* Key 54  'Magnitude[mag]' 'phot.mag;stat.mean;em.opt'  */
       char        phot_variable_flag[14]; /* Key 55  'Dimensionless[see description]' 'meta.code;src.var'  */
     double                      l       ; /* Key 56  'Angle[deg]' 'pos.galactic.lon'  */
     double                      b       ; /* Key 57  'Angle[deg]' 'pos.galactic.lat'  */
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
  char *inBuffer;
  U4HPM u4hpm_table[MAX_U4HPM_RECORDS];
  PU4HPM pU4hpm = NULL;
  int hpmLines = 0;
  int hpmIndex;
  int hpmCount = 0;
	int hpmFlag;
  int lineLen;
  char output_name[MAX_INPUT_NAME];
  File output_handle = NULL;
  int outputCount = 0;
  char input_directory[MAX_INPUT_NAME];
  char input_file_index = 0;
  char input_name[MAX_INPUT_NAME];
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
  int input_nrecs = 0;
  int zone_nrecs = 0;
  char *arrayptr;
  char nullstr[1] = {0};
  int input_index;
  size_t nbr;
  
  PGAIASTAR input_table;
  PGAIASTAR pInput;
  GSCIMAGE gsc_record;
  PGSCIMAGE pGscImage = &gsc_record;

 
  double tmpRaPM;
  double tmpDecPM;
  double tmpPM;
  double maxRaPM = 0;
  double maxDecPM = 0;
  double maxPM = 0;
  int countPM = 0;
  int tempmag;
  int status = 0;
  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int magcount[MAGBINS+1];
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
  int ncols;
  long nOptim=0;
  int columnnumber = 0;
  char keyname[80];
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
  int maxtype = 0;
  int maxform = 0;
  int maxnull = 0;
  int maxucd = 0;
  int maxunit = 0;
  int maxcomment = 0;
  int selectedKeyIndex;
  int selectedKeys[] = {1,2,4,6,7,8,9,10,11,13,14,15,16,27,28,29,30,31,32,33,34,35,36,41,42,54,55,56,57};
  int numSelectedKeys = sizeof(selectedKeys)/sizeof(int);

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


#ifdef REFNUMBER_HACK
  printf("ERROR: REFNUMBER_HACK is set - fix GetREFNumber for Gaia\n");
#endif /* REFNUMBER_HACK */

	memset(sigmaPMHist,0,sizeof(sigmaPMHist));
#ifdef HIGH_PROPER_MOTION_CUTOFF
	printf("ERROR: HIGH_PROPER_MOTION_CUTOFF is set to %f\n",HIGH_PROPER_MOTION_CUTOFF);
#endif /* HIGH_PROPER_MOTION_CUTOFF */

  /* Loop through the arguments */
  output_name[0] = 0;
  input_directory[0] = 0;

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

        case 'i': /* input directory name */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(input_directory,*++argv,MAX_INPUT_NAME-2);
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

  if (input_directory[0] == 0) {
    printf("ERROR: Input directory not specified\n");
    errorFlag = 1;
  }


  if (errorFlag) {
    printf("Usage: formatgaia -i <input directory> -o <output file> [-v][-s]\n");
    printf("       where -v is the verbose flag\n");
    printf("             -s provides statistics only and skips writing output files\n");
    printf("             <input directory> input directory containing TgasSource_000-000-000.fits through TgasSource_000-000-015.fits\n");
    printf("             <output file> output starbase file name\n");

    return(-1);
  }

  printf("formatgaia of %s %s Output Filename %s\n",
         __DATE__,__TIME__,output_name);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of GAIASTAR is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(GAIASTAR),sizeof(GSCIMAGE));
 

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


  time(&startTime);
  input_table = (PGAIASTAR)calloc(MAXRECORDS,sizeof(GAIASTAR));
  if (input_table == NULL) {
    fprintf(stderr,"ERROR: failed to allocate %d bytes for %d records\n",MAXRECORDS*sizeof(GAIASTAR),MAXRECORDS);
  }
    
  
  
  fprintf(output_handle,"REFNumber\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tRaSigmaPM\tDecSigmaPM\n");
  fprintf(output_handle,"---------\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t---------\t----------\n");
  /* Now read in the input file */
  for (input_file_index = 0; input_file_index < MAX_FILE_INDEX; input_file_index++) {

    sprintf(input_name,"%sTgasSource_000-000-%03d.fits",input_directory,input_file_index);
    fits_open_table(&fptr, input_name, READONLY, &status);
    if (status != 0) {
      errorFlag = 1;
      printf("Could not open input file %s\n",input_name);
      exit(-1);
    }
    fits_get_num_rows (fptr, &nrows, &status);
    if ( status != 0 ) {
      fits_report_error(stderr, status);
      return(-1);
    }
    totnrows += nrows;
    fits_get_num_cols (fptr, &ncols, &status);
    if ( status != 0 ) {
      fits_report_error(stderr, status);
      return(-1);
    }
    fits_get_rowsize(fptr, &nOptim, &status);
    if ( status != 0 ) {
      fits_report_error(stderr, status);
      return(-1);
    }
    printf("nrows %d ncols %d nOptim %d for %s\n",nrows,ncols,nOptim,input_name);


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
    /*
     * Key  1 name:                      'hip     ' format: 'J       ' null: -2147483648 ucd:                         'meta.id.cross' comment:                                        
     * Key  2 name:                     'tycho2_id' format: '11A     ' null:             ucd:                         'meta.id.cross' comment:                                        
     * Key  4 name:                     'source_id' format: 'K       ' null:             ucd:                     'meta.id;meta.main' comment:                                        
     * Key  6 name:                     'ref_epoch' format: 'D       ' null:             ucd:                   'meta.ref;time.epoch' comment:                                        
     * Key  7 name:                      'ra      ' format: 'D       ' null:             ucd:                   'pos.eq.ra;meta.main' comment:                                        
     * Key  8 name:                      'ra_error' format: 'D       ' null:             ucd:                  'stat.error;pos.eq.ra' comment:                                        
     * Key  9 name:                      'dec     ' format: 'D       ' null:             ucd:                  'pos.eq.dec;meta.main' comment:                                        
     * Key 10 name:                     'dec_error' format: 'D       ' null:             ucd:                 'stat.error;pos.eq.dec' comment:                                        
     * Key 11 name:                      'parallax' format: 'D       ' null:             ucd:                          'pos.parallax' comment:                                         unit:                     'Angle[mas]'
     * Key 13 name:                      'pmra    ' format: 'D       ' null:             ucd:                      'pos.pm;pos.eq.ra' comment:                                        
     * Key 14 name:                    'pmra_error' format: 'D       ' null:             ucd:           'stat.error;pos.pm;pos.eq.ra' comment:                                        
     * Key 15 name:                      'pmdec   ' format: 'D       ' null:             ucd:                     'pos.pm;pos.eq.dec' comment:                                        
     * Key 16 name:                   'pmdec_error' format: 'D       ' null:             ucd:          'stat.error;pos.pm;pos.eq.dec' comment:                                        
     * Key 27 name:          'astrometric_n_obs_al' format: 'J       ' null:             ucd:                           'meta.number' comment:                                        
     * Key 28 name:          'astrometric_n_obs_ac' format: 'J       ' null:             ucd:                           'meta.number' comment:                                        
     * Key 29 name:     'astrometric_n_good_obs_al' format: 'J       ' null:             ucd:                           'meta.number' comment:                                        
     * Key 30 name:     'astrometric_n_good_obs_ac' format: 'J       ' null:             ucd:                           'meta.number' comment:                                        
     * Key 31 name:      'astrometric_n_bad_obs_al' format: 'J       ' null:             ucd:                           'meta.number' comment:                                        
     * Key 32 name:      'astrometric_n_bad_obs_ac' format: 'J       ' null:             ucd:                           'meta.number' comment:                                        
     * Key 33 name:           'astrometric_delta_q' format: 'E       ' null:             ucd:                            'stat.value' comment:                                        
     * Key 34 name:      'astrometric_excess_noise' format: 'D       ' null:             ucd:                            'stat.value' comment:                                        
     * Key 35 name:  'astrometric_excess_noise_sig' format: 'D       ' null:             ucd:                            'stat.value' comment:                                        
     * Key 36 name:      'astrometric_primary_flag' format: 'L       ' null:             ucd:                             'meta.code' comment:                                        
     * Key 41 name:          'matched_observations' format: 'I       ' null:             ucd:                           'meta.number' comment:                                        
     * Key 42 name:             'duplicated_source' format: 'L       ' null:             ucd:                                         comment:                                        
     * Key 54 name:               'phot_g_mean_mag' format: 'D       ' null:             ucd:             'phot.mag;stat.mean;em.opt' comment:                                         unit:                 'Magnitude[mag]'
     * Key 55 name:            'phot_variable_flag' format: '13A     ' null:             ucd:                     'meta.code;src.var' comment:                                         unit: 'Dimensionless[see description]'
     * Key 56 name:                      'l       ' format: 'D       ' null:             ucd:                      'pos.galactic.lon' comment:                                        
     * Key 57 name:                      'b       ' format: 'D       ' null:             ucd:                      'pos.galactic.lat' comment:                                        
     */
    columnnumber = 1;
    while (1) {
      sprintf(keyname,"TTYPE%d",columnnumber);
      fits_read_keyword(fptr,keyname,typevalue,comment,&status);
      if (status != 0) {
        break;
      }
      sprintf(keyname,"TFORM%d",columnnumber);
      fits_read_keyword(fptr,keyname,formvalue,comment,&status);
      if (status != 0) {
        break;
      }
      sprintf(keyname,"TNULL%d",columnnumber);
      status = 0;
      fits_read_keyword(fptr,keyname,nullvalue,comment,&status);
      if (status != 0) {
        nullvalue[0] = 0;
      }
      sprintf(keyname,"TUCD%d",columnnumber);
      status = 0;
      fits_read_keyword(fptr,keyname,ucdvalue,comment,&status);
      if (status != 0) {
        ucdvalue[0] = 0;
      }
      sprintf(keyname,"TCOMM%d",columnnumber);
      status = 0;
      fits_read_keyword(fptr,keyname,commentvalue,comment,&status);
      if (status != 0) {
        commentvalue[0] = 0;
      }
      sprintf(keyname,"TUNIT%d",columnnumber);
      status = 0;
      fits_read_keyword(fptr,keyname,unitvalue,comment,&status);
      if (status != 0) {
        unitvalue[0] = 0;
      }
      status = 0;
      if (strlen(typevalue) > maxtype) {
        maxtype = strlen(typevalue);
      }
      if (strlen(formvalue) > maxform) {
        maxform = strlen(formvalue);
      }
      if (strlen(nullvalue) > maxnull) {
        maxnull = strlen(nullvalue);
      }
      if (strlen(ucdvalue) > maxucd) {
        maxucd = strlen(ucdvalue);
      }
      if (strlen(commentvalue) > maxcomment) {
        maxcomment = strlen(commentvalue);
      }
      if (strlen(unitvalue) > maxunit) {
        maxunit = strlen(unitvalue);
      }
      formatname[0] = 0;
      countname[0] = 0;
      charsize = 0;
#if 0
      printf("Key %2d name: %31s format: %10s null: %11s ucd: %39s comment: %39s unit: %32s\n",
             columnnumber,
             typevalue,
             formvalue,
             nullvalue,
             ucdvalue,
             commentvalue,
             unitvalue);
#endif
      switch (formvalue[1]) {
      case 'D':
        strcpy(formatname,"double");
        strcpy(datatypename,"TDOUBLE");
        break;
      case 'E':
        strcpy(formatname,"float");
        strcpy(datatypename,"TFLOAT");
        break;
      case 'I':
        strcpy(formatname,"int16_t");
        strcpy(datatypename,"TSHORT");
        break;
      case 'J':
        strcpy(formatname,"int32_t");
        strcpy(datatypename,"TINT32BIT");
        break;
      case 'K':
        strcpy(formatname,"int64_t");
        strcpy(datatypename,"TLONGLONG");
        break;
      case 'L':
        strcpy(formatname,"char");
        strcpy(datatypename,"TLOGICAL");
        break;
      default:
        nvals = sscanf(&formvalue[1],"%dA",&charsize);
        if (nvals == 1) {
          strcpy(formatname,"char");
          sprintf(countname,"[%d]",charsize+1);
        } else {
          sprintf(formatname,"UNKNOWN");
          break;
        }
      }
      charPtr = strstr(&typevalue[1],"'");
      if (charPtr != NULL) {
        *charPtr = 0;
      }
      for (selectedKeyIndex = 0; selectedKeyIndex < numSelectedKeys; selectedKeyIndex++) {
        if (columnnumber == selectedKeys[selectedKeyIndex]) {
          break;
        }
      } 
      if (selectedKeyIndex <  numSelectedKeys) {

#if 0
        printf(" %10s %29s%s; /* Key %2d %s %s %s %s */\n",
               formatname,
               &typevalue[1],
               countname,
               columnnumber,
               nullvalue,
               unitvalue,
               ucdvalue,
               commentvalue);
#endif               
#if 0
        if (strcmp(formatname,"char") != 0) {
          printf("fits_read_col(fptr,%s,%d,zone_nrecs+1,1,1,0,&pInput->%s,&anynull,&status);\n",
                 datatypename,
                 columnnumber,
                 &typevalue[1]);
        } else {
          if (charsize == 0) {
            printf("arrayptr = &pInput->%s;\n",&typevalue[1]);
          } else {
            printf("arrayptr = &pInput->%s[0];\n",&typevalue[1]);
          }
          printf("fits_read_col_str(fptr,%d,zone_nrecs+1,1,1,nullstr,&arrayptr,&anynull,&status);\n",
                 columnnumber);
        }

        printf("if (status != 0) {\n");
        printf("gaia_print_error(%d,zone_nrecs+1,\"%s\",status,__LINE__);\n",columnnumber, &typevalue[1]);
        printf("exit(-1);\n");
        printf("}\n");
#endif        
      }

      columnnumber++;
    }
    columnnumber--;
    /* Read the rows we want */
    for (zone_nrecs = 0; zone_nrecs < nrows; zone_nrecs++) {
      pInput = &input_table[input_nrecs];
      memset(pInput,0,sizeof(GAIASTAR));
      status = 0;
      fits_read_col(fptr,TINT32BIT,1,zone_nrecs+1,1,1,0,&pInput->hip     ,&anynull,&status);
      if (status != 0) {
        gaia_print_error(1,zone_nrecs+1,"hip     ",status,__LINE__);
        exit(-1);
      }
      arrayptr = &pInput->tycho2_id[0];

      fits_read_col(fptr,TINT32BIT,1,zone_nrecs+1,1,1,0,&pInput->hip     ,&anynull,&status);
      if (status != 0) {
        gaia_print_error(1,zone_nrecs+1,"hip     ",status,__LINE__);
        exit(-1);
      }
      arrayptr = &pInput->tycho2_id[0];
      fits_read_col_str(fptr,2,zone_nrecs+1,1,1,nullstr,&arrayptr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(2,zone_nrecs+1,"tycho2_id",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TLONGLONG,4,zone_nrecs+1,1,1,0,&pInput->source_id,&anynull,&status);
      if (status != 0) {
        gaia_print_error(4,zone_nrecs+1,"source_id",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,6,zone_nrecs+1,1,1,0,&pInput->ref_epoch,&anynull,&status);
      if (status != 0) {
        gaia_print_error(6,zone_nrecs+1,"ref_epoch",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,7,zone_nrecs+1,1,1,0,&pInput->ra      ,&anynull,&status);
      if (status != 0) {
        gaia_print_error(7,zone_nrecs+1,"ra      ",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,8,zone_nrecs+1,1,1,0,&pInput->ra_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(8,zone_nrecs+1,"ra_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,9,zone_nrecs+1,1,1,0,&pInput->dec     ,&anynull,&status);
      if (status != 0) {
        gaia_print_error(9,zone_nrecs+1,"dec     ",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,10,zone_nrecs+1,1,1,0,&pInput->dec_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(10,zone_nrecs+1,"dec_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,11,zone_nrecs+1,1,1,0,&pInput->parallax,&anynull,&status);
      if (status != 0) {
        gaia_print_error(11,zone_nrecs+1,"parallax",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,13,zone_nrecs+1,1,1,0,&pInput->pmra    ,&anynull,&status);
      if (status != 0) {
        gaia_print_error(13,zone_nrecs+1,"pmra    ",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,14,zone_nrecs+1,1,1,0,&pInput->pmra_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(14,zone_nrecs+1,"pmra_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,15,zone_nrecs+1,1,1,0,&pInput->pmdec   ,&anynull,&status);
      if (status != 0) {
        gaia_print_error(15,zone_nrecs+1,"pmdec   ",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,16,zone_nrecs+1,1,1,0,&pInput->pmdec_error,&anynull,&status);
      if (status != 0) {
        gaia_print_error(16,zone_nrecs+1,"pmdec_error",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TINT32BIT,27,zone_nrecs+1,1,1,0,&pInput->astrometric_n_obs_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(27,zone_nrecs+1,"astrometric_n_obs_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TINT32BIT,28,zone_nrecs+1,1,1,0,&pInput->astrometric_n_obs_ac,&anynull,&status);
      if (status != 0) {
        gaia_print_error(28,zone_nrecs+1,"astrometric_n_obs_ac",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TINT32BIT,29,zone_nrecs+1,1,1,0,&pInput->astrometric_n_good_obs_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(29,zone_nrecs+1,"astrometric_n_good_obs_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TINT32BIT,30,zone_nrecs+1,1,1,0,&pInput->astrometric_n_good_obs_ac,&anynull,&status);
      if (status != 0) {
        gaia_print_error(30,zone_nrecs+1,"astrometric_n_good_obs_ac",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TINT32BIT,31,zone_nrecs+1,1,1,0,&pInput->astrometric_n_bad_obs_al,&anynull,&status);
      if (status != 0) {
        gaia_print_error(31,zone_nrecs+1,"astrometric_n_bad_obs_al",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TINT32BIT,32,zone_nrecs+1,1,1,0,&pInput->astrometric_n_bad_obs_ac,&anynull,&status);
      if (status != 0) {
        gaia_print_error(32,zone_nrecs+1,"astrometric_n_bad_obs_ac",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TFLOAT,33,zone_nrecs+1,1,1,0,&pInput->astrometric_delta_q,&anynull,&status);
      if (status != 0) {
        gaia_print_error(33,zone_nrecs+1,"astrometric_delta_q",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,34,zone_nrecs+1,1,1,0,&pInput->astrometric_excess_noise,&anynull,&status);
      if (status != 0) {
        gaia_print_error(34,zone_nrecs+1,"astrometric_excess_noise",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,35,zone_nrecs+1,1,1,0,&pInput->astrometric_excess_noise_sig,&anynull,&status);
      if (status != 0) {
        gaia_print_error(35,zone_nrecs+1,"astrometric_excess_noise_sig",status,__LINE__);
        exit(-1);
      }
      arrayptr = &pInput->astrometric_primary_flag;
      fits_read_col_str(fptr,36,zone_nrecs+1,1,1,nullstr,&arrayptr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(36,zone_nrecs+1,"astrometric_primary_flag",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TSHORT,41,zone_nrecs+1,1,1,0,&pInput->matched_observations,&anynull,&status);
      if (status != 0) {
        gaia_print_error(41,zone_nrecs+1,"matched_observations",status,__LINE__);
        exit(-1);
      }
      arrayptr = &pInput->duplicated_source;
      fits_read_col_str(fptr,42,zone_nrecs+1,1,1,nullstr,&arrayptr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(42,zone_nrecs+1,"duplicated_source",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,54,zone_nrecs+1,1,1,0,&pInput->phot_g_mean_mag,&anynull,&status);
      if (status != 0) {
        gaia_print_error(54,zone_nrecs+1,"phot_g_mean_mag",status,__LINE__);
        exit(-1);
      }
      arrayptr = &pInput->phot_variable_flag[0];
      fits_read_col_str(fptr,55,zone_nrecs+1,1,1,nullstr,&arrayptr,&anynull,&status);
      if (status != 0) {
        gaia_print_error(55,zone_nrecs+1,"phot_variable_flag",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,56,zone_nrecs+1,1,1,0,&pInput->l       ,&anynull,&status);
      if (status != 0) {
        gaia_print_error(56,zone_nrecs+1,"l       ",status,__LINE__);
        exit(-1);
      }
      fits_read_col(fptr,TDOUBLE,57,zone_nrecs+1,1,1,0,&pInput->b       ,&anynull,&status);
      if (status != 0) {
        gaia_print_error(57,zone_nrecs+1,"b       ",status,__LINE__);
        exit(-1);
      }
#if 0
      /* precess to J2000 position */
      fk5prec(pInput->ref_epoch,2000.0,&pInput->ra,&pInput->dec);
#endif
      input_nrecs++;


    }

    status = 0;
    time(&curTime);
    curTime -= startTime;
    printf("Finished with record %d for file %s at %d seconds\n",input_nrecs,input_name,curTime);
    if ( fits_close_file(fptr, &status) ) {
      fits_report_error(stderr, status);
      return(-1);
    }
    

  }


  printf("totnrows %d total columns %d maxtype %d maxform %d maxnull %d maxucd %d maxcomment %d maxunit %d\n",totnrows,columnnumber,maxtype,maxform,maxnull,maxucd,maxcomment,maxunit);
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
    if (pInput->phot_g_mean_mag >= MAGBINS) {
      tempmag = MAGBINS;
    } else if (pInput->phot_g_mean_mag < 0) {
      tempmag = 0;
    } else {
      tempmag = pInput->phot_g_mean_mag;
    }
    magcount[tempmag]++;
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
        

    /* Now fill in our GSC record */
    sprintf(REF,"GAIA1_%lld",pInput->source_id); /* A preceeding 'GAIA1_' flags gaia DR1 ids */
#ifdef REFNUMBER_HACK
    pGscImage->REFNumber = pInput->source_id;
#else /* REFNUMBER_HACK */

#if 0
    if (strcmp(REF,"GAIA1_1000000534063130624") == 0) {
      printf("At REF %s\n",REF);
    }
#endif
    if ((pInput->source_id % GAIA_MODULUS) != 0) {
      printf("ERROR: formatgaia invalid modulus %d for source_id %lld %s\n",GAIA_MODULUS,pInput->source_id,REF);
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
#endif

    pGscImage->ra =  pInput->ra;
    pGscImage->dec = pInput->dec;

    pGscImage->class = 0;
    pGscImage->Stdmag = pInput->phot_g_mean_mag;
    pGscImage->color = 99.0;
    

    pGscImage->RaPM = pInput->pmra;
    pGscImage->RaSigmaPM = pInput->pmra_error;
    pGscImage->DecPM = pInput->pmdec;
    pGscImage->DecSigmaPM = pInput->pmdec_error;

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

  if (input_table != NULL) {
    free(input_table);
  }

  Close(output_handle);
 
  time(&curTime);
  curTime -= startTime;

  printf("min_id_number %lld, max_id_number %lld hpmStars %d\n",min_id_number,max_id_number,hpmCount);
  printf("Execution Time: %d seconds for %d records; output %d \n",curTime,input_nrecs,outputCount);

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

  if (properMotionHandle != NULL) {
    fclose(properMotionHandle);
  }
  if (newIdHandle != NULL) {
    fclose(newIdHandle);
  }

  return(EXIT_SUCCESS);
}
