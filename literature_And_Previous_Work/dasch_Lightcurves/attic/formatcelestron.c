// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* formatcelestron.c
 *
 *   Catalog comes from an older version of the Celestron catalog; and the latest version Celestron remote control 
 *   
 *   See ~/Pipeline/catalogs/celestron_st/README.txt and celestron_sc_v1.0_2021_11_15C.pdf
 *   Sources for celestron_sc_v1.0_2021_11_15C.pdf are in celestron_sc_v1.0_2021_11_15C.odt
 
     The "st" table was hand-transcribed from 

     Celestron Advanced Astro Master (#93900) catalog of 10,000 objects
      https://www.yumpu.com/en/document/view/11511012/celestron-advanced-astro-master-operating-manual 
      pages 35-37 "The Star Catalog" ST001 to ST591
    
     The "nt" table comes from the Celestron 6 SE and Celestron 8 SE telescope 
      NexRemote, NexTour and NexGpa version 1.5.0 applications documented in 
      NexStar_6_SE_and_8_SE_Manual_1152120348_1106811069manua.pdf
      Copyright 2006 Celestron (www.celestron.com)
      The NexTour.mdb table is a Windows JET4 database read by the open-source
      mdb-tables utility written by Brian Bruns (https://github.com/mdbtools/mdbtools)
      to produce the file 
      /home/scanner/Pipeline/catalogs/celestron_st/Master.csv containing 45927 entries with sequential "ID" numbers from 1 through 46036

     The combined table is assigned the following numbers
     CELESTRON_00001 NT id 00001
     CELESTRON_46999 NT id 46999
     CELESTRON_47001 ST id 00001
     CELESTRON_47591 ST id 00591
 * 
 *
 * Read the celestron input catalog and put it into the gsc232bin.db format
 *        Accept everything:               formatcelestron -s -v -n -o /dasch/Pipeline/catalogs/celestron.db
 *
 echo "formatcelestron -v -c 0.1 -o /dasch/Pipeline/catalogs/celestron_temp.db >& /dasch/raid020/junk/celestron.log" | at now

 formatgsc  -q celestron /dasch/Pipeline/catalogs/celestron_temp.db

 *
 *
 * here we do not use MySQL
 *

gcc -fPIC -ggdb -c -O2 -DNO_MYSQL  -I/dasch/install/include  pipelineutils.c  -o pipelineutils.o


echo "formatcelestron"
cc -ggdb -O2   -DNO_MYSQL -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  formatcelestron.c -L /dasch/install/lib -lm   -ltable -lutil -lwcs   -ldl -pthread -o pipelineutils.o -o formatcelestron


 * here is what is should be
echo "formatcelestron"
cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  formatcelestron.c pipelineutils.a -L /dasch/install/lib -lm   -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -o formatcelestron



 *  1.  Use this program to reformat the celestron catalog
 *        Accept nobs > 2, rms cutoff 0.1: formatcelestron -v -c 0.1 -o /dasch/Pipeline/catalogs/celestron_temp.db
 *        Execution Time: 529 seconds; rejected because of color 1390 stars written 7173325
 *  
 *
 *  2.  Sort the results into 0.01 degree bins:
 *      
 *      index -mb -n celestron_temp.db dec
 *      (takes 2m24s on dell)
 *
 *  4.  Run formatgsc to put this file in binary format and to index it
 *      
 *      formatgsc -r 1.0 -q celestron /dasch/Pipeline/catalogs/celestron_temp.db
 * 
 *  5.  build matchcatalogs.c with CELESTRON_CATALOG, INCLUDE_KEPLER_COLOR, WRITE_MATCH_STARBASE and  WRITE_MERGE_CATALOG defined (no longer necessary)
 *      
 *      /dasch/Pipeline/matchcatalogs -c a -t 2.0 -k -m -w -f /dasch/Pipeline/catalogs/ucac3.dat -s /dasch/Pipeline/catalogs/celestron_temp.dat
 *      cd /dasch/Pipeline/catalogs
 *      mv merge.idx celestron.idx
 *      mv merge.dat celestron.dat
 *      rm celestron_temp.dat
 *      rm celestron_temp.idx
 *      rm celestron_temp_zero.db
 *
 *
 *     matchcatalogs of Mar 23 2011 17:19:47, 
 *       matching /dasch/Pipeline/catalogs/ucac3.dat and /dasch/Pipeline/catalogs/celestron_temp.dat 
 *       output /dasch/Pipeline/catalogs/matchgsckepler.db total bins 168966386 TOLERANCE 0.000556 degrees
 *      Done  populatedBinCount 6121155 matches 6292508 for 7173325 kepler stars at 699 seconds
 *            median degDrad is 0.000085 deg or 0.3 arcsec max is 0.000556 or 2.0 arcsec.
 *      mv merge.dat celestron.dat
 *      mv merge.idx celestron.idx
 *      rm celestron_temp*
 *
 * Nov 23, 2021 Edward J. Los - Initial version, adapted from formatatlas.c
 *                                      pTab table adapted from filterdr1.c
 *                                      _galaxyrec from galaxyutils.h
 * Dec  6, 2021 Edward J. Los - Add -t <maxTabStorageSize> with default MAX_TAB_STORAGE_SIZE
 *                                      introduce gscbinutils.h from pipelineutils.c
 * Dec  9, 2021 Edward J. Los - Split the string table and string index
 *                              MAX_TAB_STORAGE_SIZE and STRING_TABLE_OVERHEAD go away
 *                              in favor of MAX_TAB_LENGTH.  maxTabStorageSize becomes maxTabLength
 *                              globalByteStorageSize becomes globalStaticStorageSize;
 *                              globalAllocatedStorageSize becomes globalAllocatedStorageSize
 *                              globalTotalStorageSize goes away
 * Dec 14, 2021 Edward J. Los - When a separation character is found, store the individual tabs in
 *                               char tabTable[MAX_TAB+1][MAX_BUFFER+1] and the indices in
 *                               int tabIndexTable[MAX_TAB+1];
 *                              The tabIndexTable doesn't work well.  For each new tab found, store only
 *                               the buffer AFTER the tab and search the previously stored buffer for the
 *                               separation character and zero that character.
 */ 
  


#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "gscbinutils.h"
#include "galaxyutils.h"
#include "pipelineutils.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_INPUT_NAME 512
#define MAX_BUFFER 1024
#define MAX_LINES 45928
#define MAGBINS  30  /* Number of magnitude bins */
#define CELESTRON_BP     0
#define CELESTRON_RP     1
#define CELESTRON_G      2
#define CELESTRON_R      3
#define CELESTRON_I      4
#define CELESTRON_Z      5
#define CELESTRON_MAXMAG 6
#define MAX_MODULUS  1 /* Base 10 modulus */
#define MODULUS_BASE 10
#define MAX_SKIP_BIN 1000 /* dec bins skipped between entries */
#define SKIP_BUFFER 4 /* extra lines to read before the skip point */
#define CELESTRON_MAXRMS 1.0
#define RMSBINS 100
#define RESIDENT_DEC_BINS 128
#define MAX_DEC_BIN 232860 /* dec_bin 3903 at Declination  -29.016 */
#define MAX_INPUT_LINES 992637834
#define GAIADR2_EPOCH 2015.5
#define MAX_TAB 60 /* No more than 60 tabs in a record */
#define MAX_TAB_LENGTH 400
#define TAB_BIN_SIZE 32
#if 1
#define STRING_TABLE_LENGTH 1000000 /* Length of the stringTable */
#else
#define STRING_TABLE_LENGTH 2048 /* Length of the stringTable */
#endif
/*

MAX_TAB_STORAGE_SIZE =  0 Byte:        0 Array:  6225811 Total:  6225811 Allocated  6225811
MAX_TAB_STORAGE_SIZE =  8 Byte:  4629691 Array:  1595832 Total:  6225523 Allocated  6222421
MAX_TAB_STORAGE_SIZE = 16 Byte:  5899597 Array:   325126 Total:  6224723 Allocated  6222421
MAX_TAB_STORAGE_SIZE = 32 Byte:  6021141 Array:   200446 Total:  6221587
MAX_TAB_STORAGE_SIZE = 64 Byte:  6026764 Array:   184423 Total:  6211187
MAX_TAB_STORAGE_SIZE =128 Byte:  6044683 Array:   125992 Total:  6170675
MAX_TAB_STORAGE_SIZE =400 Byte:  6078027 Array:        0 Total:  6078027

globalTotalTabs 885698  average size is 7.0 bytes
*/

#define NO_GETDECBIN 1
/* #define PLOT_RMS_GMAG */ /* starbase count of gmag vs rms */
#ifndef  NO_GETDECBIN
extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
#endif /* NO_GETDECBIN */

/* Temporary definitions */
#define MIN_CELESTRON_COLOR -99.9
#define MAX_CELESTRON_COLOR +99.9

double sqr(double a) {
  return(a*a);
}
/* 
   stringTable is a simple global storage array for very long strings of size >= MAX_TAB_STORAGE_SIZE 
   each entry consists of an absolute address of off_t bytes 
                          followed by a count of off_t bytes 
   and then the string itself aligned to an off_t byte boundry
*/
typedef struct _stringIndex {
  char *pStringEntry;
  off_t stringLength;
} STRINGINDEX,*PSTRINGINDEX;

STRINGINDEX stringTableIndex[STRING_TABLE_LENGTH];
char stringTable[STRING_TABLE_LENGTH];
char stringBuffer[2*MAX_BUFFER];
char tabTable[MAX_TAB+1][MAX_BUFFER+1];
int tabIndexTable[MAX_TAB+1];
int globalStringTableEntries = 0; /* Number of entries in the table */
off_t globalStaticStorageSize = 0;
off_t globalAllocatedStorageSize = 0;
off_t stringTableLast = 0;        /* Last entry in the table = total bytes allocated */
int DumpCounter = 0;

static void DumpStringTable(int lineno,char *fileName);
static void InitStringTable(int lineno,char *fileName) {
  if (globalStringTableEntries =  0) {
    return;
  }
  if (globalStringTableEntries != 0) {
    printf("ERROR line %4d InitStringTable called twice. fileName '%s'\n",__LINE__,fileName);
    exit(-1);
  }
  memset(stringTable,0,sizeof(stringTable));
  memset(stringTableIndex,0,sizeof(stringTableIndex));
  globalStaticStorageSize += sizeof(stringTable);
  globalStaticStorageSize += sizeof(stringTableIndex);
  globalStaticStorageSize += sizeof(tabTable);
  globalStaticStorageSize += sizeof(tabIndexTable);
  DumpStringTable(lineno,fileName);

  return;
}
static void AllocateStringTable(PSTRINGINDEX *pStringTableIndex,off_t stringBytes ,int lineno,char *fileName) {
  char* stringTableAddress;
  PSTRINGINDEX pStringIndex;

  DumpStringTable(lineno,fileName);
  if (globalStringTableEntries == 0) {
    InitStringTable(lineno,fileName);
  }
  if ((stringBytes +stringTableLast) > STRING_TABLE_LENGTH) {
    printf("ERROR line %4d: STRING_TABLE_LENGTH %lld cur %lld new %lld exceeded for file '%s'\n",lineno,STRING_TABLE_LENGTH,stringTableLast,stringBytes +stringTableLast,fileName);
    exit(-1);
  }
#if 0
  if (*pStringTableIndex != NULL) {
    printf("ERROR line %4d: Non-null pStringTableIndex for file '%s'\n",lineno,fileName);
    exit(-1);
  }
#endif
  pStringIndex = (PSTRINGINDEX)&stringTableIndex[globalStringTableEntries];
  
  globalStringTableEntries++;
  if (pStringIndex->pStringEntry != 0) {
    printf("ERROR: line %4d: Non-null pStringIndex->pStringEntry  for file '%s'\n",lineno,fileName);
    exit(-1);
  }
  if (pStringIndex->stringLength != 0) {
    printf("ERROR: line %4d: Non-null pStringIndex->stringLength  for file '%s'\n",lineno,fileName);
    exit(-1);
  }
  stringTableAddress = (char *)(&stringTable[stringTableLast]);
  pStringIndex->pStringEntry = stringTableAddress;
  pStringIndex->stringLength = stringBytes;
  stringTableLast += stringBytes;
  DumpStringTable(lineno,fileName);
  *pStringTableIndex = pStringIndex;
}
static void WriteStringTable(PSTRINGINDEX pStringIndex,char* pTab,off_t stringBytes ,int lineno,char *fileName) {
  if (stringBytes  >= (MAX_BUFFER-3)) {
    printf("ERROR line %4d: MAX_BUFFER  exceeded with %lld for file '%s'\n",lineno,stringBytes ,fileName);
    exit(-1);
  }
#if 0
  printf("At line %4d DumpCounter %d for file '%s'\n",__LINE__,DumpCounter,fileName);
  {
    int index5;
    int index6;
    PSTRINGINDEX pStringIndex5;
    for (index5 = 0; index5 < 20; index5++) {
      pStringIndex5 = &stringTableIndex[index5];
      printf("I: %3d E: 0x%08lld L: %04llx",index5,pStringIndex5->pStringEntry,pStringIndex5->stringLength);
      fflush(stdout);
      if (((index5+1) % 5) == 0) {
        printf("\n");
      } else {
        printf(" ");
      }
    }
    for (index6 = 0; index6 < 420; index6++) {
      if (((index6+1) % 40) == 1) {
        printf("I %03d ",index6);
      }
      if (stringTable[index6] == 0) {
        printf(". ");
      } else {
        printf(".%01c",stringTable[index6]);
      }
      if (((index6+1) % 40) == 0) {
        printf("\n",stringTable[index6]);
      } else {
        printf(" ",stringTable[index6]);
      }
      fflush(stdout);
    }
    printf("\n");
  }
#endif
#if 0
  strcpy(stringBuffer,"ABC");
#endif
  strcpy(pStringIndex->pStringEntry,pTab);
#if 0
  printf("At line %4d DumpCounter %d for file '%s'\n",__LINE__,DumpCounter,fileName);
#endif
}
static void GetStringTable(PSTRINGINDEX pStringTableIndex,void *stringPtr,off_t stringBytes ,int lineno,char *fileName) {
  char *pStringEntry = pStringTableIndex->pStringEntry;
  off_t stringLength = pStringTableIndex->stringLength;  

  if (stringBytes  >= (MAX_BUFFER-3)) {
    printf("ERROR line %4d: MAX_BUFFER  exceeded with %lld for file '%s'\n",lineno,stringBytes ,fileName);
    exit(-1);
  }
  if (stringBytes  >= (stringLength )) {
    printf("ERROR line %4d: stringBytes  %lld  exceeded string size %lld for file '%s'\n",lineno,stringBytes ,stringLength,fileName);
    exit(-1);
  }
}
static void DumpStringTable(int lineno,char *fileName) {
#if 0
  off_t stringIndex;
  PSTRINGINDEX pStringIndex;
  char* stringTableAddress;
  off_t stringTableBytes;
  DumpCounter++;
  if (DumpCounter == 3) {
    printf("At line %4d for file '%s'\n",__LINE__,fileName);
  }
  printf("\nDEBUG line %4d: DumpCounter %4d DumpStringTable for file '%s'\n",lineno,DumpCounter,fileName);
  printf("line %4d: globalStringTableEntries %lld\n",__LINE__,globalStringTableEntries);
  printf("line %4d: stringTable  start       0x%llx\n",__LINE__,stringTable);
  printf("line %4d: stringTableLast          0x%llx %lld\n",__LINE__,stringTableLast,stringTableLast);
  printf("line %4d: stringBuffer             0x%08llx '%s'\n",__LINE__,stringBuffer,stringBuffer);
  for (stringIndex = 0; stringIndex < globalStringTableEntries; stringIndex++) {
    pStringIndex = &stringTableIndex[stringIndex];
    stringTableAddress = pStringIndex->pStringEntry;
    stringTableBytes  = pStringIndex->stringLength;
    printf("line %4d: stringIndex %5d, stringTableAddress 0x%08llx stringTableBytes %08lld stringTableAddress '%s'\n",
           __LINE__,
           stringIndex,
           stringTableAddress,
           stringTableBytes,
           stringTableAddress);
  }

  printf("\n");
  return;
#endif
}


 /* Greek alphabet from https://en.wikipedia.org/wiki/Greek_alphabet */
#define GREEK_ALPHABET_MAX    23
char *greekAlphabet[GREEK_ALPHABET_MAX+1] = 
  {
    "alpha,",
    "beta,",
    "gamma,",
    "delta,",
    "epsilon,",
    "zeta,",
    "eta,",
    "theta,",
    "iota,",
    "kappa,",
    "lambda,",
    "mu,",
    "nu,",
    "xi,",
    "omicron,",
    "pi,",
    "rho,",
    "tau,",
    "upsilon,",
    "phi,",
    "chi,",
    "psi,",
    "omega,",
  };


/* #define LOS_DEBUG 1 */  
#define VARIABLE 1      /* see dupvar below */
#define NOT_AVAILABLE 2 /* see dupvar below */
#define DUPLICATE 4     /* see dupvar below */
typedef struct _celestronimage {
  long long objid; /* (none) Object ID           */
  double RA; /* (degrees) Right ascension from Gaia DR2, J2000, epoch 2015.5     */
  double Dec; /* (degrees) Declination from Gaia DR2, J2000, epoch 2015.5      */
  double plx; /* (mas) Parallax from Gaia DR2         */
  double dplx; /* (mas) Parallax uncertainty from Gaia DR2        */
  double pmra; /* (mas/yr) Proper motion in right ascension from Gaia DR2     */
  double dpmra; /* (mas/yr) Proper motion uncertainty in right ascension       */
  double pmdec; /* (mas/yr) Proper motion in declination from Gaia DR2      */
  double dpmdec; /* (mas/yr) Proper motion uncertainty in declination        */
  double Gaia; /* (mag) Gaia G magnitude          */
  double dGaia; /* (mag) Gaia G magnitude uncertainty         */
  double BP; /* (mag) Gaia G_bp magnitude          */
  double dBP; /* (mag) Gaia G_bp magnitude uncertainty         */
  double RP; /* (mag) Gaia G_rp magnitude          */
  double dRP; /* (mag) Gaia G_rp magnitude uncertainty         */
  int Teff; /* (K) Gaia stellar effective temperature         */
  double AGaia; /* (mag) Gaia estimate of G-band extinction for this star     */
  int dupvar; /* (none) Gaia variability and duplicate flags, 0/1/2 for "CONSTANT"/"VARIABLE"/"NOT AVAILABLE" + 4*DUPLICATE  */
  double Ag; /* (mag) SFD estimate of total g-band extinction       */
  double rp1; /* (arcsec) Radius where cummulative G flux exceeds 0.1 x this star   */
  double r1; /* (arcsec) Radius where cummulative G flux exceeds 1.0 x this star   */
  double r10; /* (arcsec) Radius where cummulative G flux exceeds 10.0 x this star   */
  double g; /* (mag) PanSTARRS g magnitude          */
  double dg; /* (mag) PanSTARRS g magnitude uncertainty         */
  double gchi; /* (none) chi^2 / DOF for contributors        */
  int gcontrib; /* (none) Bitmap of conributing catalogs to g       */
  double r; /* (mag) PanSTARRS r magnitude          */
  double dr; /* (mag) PanSTARRS r magnitude uncertainty         */
  double rchi; /* (none) chi^2 / DOF for contributors        */
  int rcontrib; /* (none) Bitmap of conributing catalogs to r       */
  double i; /* (mag) PanSTARRS i magnitude          */
  double di; /* (mag) PanSTARRS i magnitude uncertainty         */
  double ichi; /* (none) chi^2 / DOF for contributors        */
  int icontrib; /* (none) Bitmap of conributing catalogs to i       */
  double z; /* (mag) PanSTARRS z magnitude          */
  double dz; /* (mag) PanSTARRS z magnitude uncertainty         */
  double zchi; /* (none) chi^2 / DOF for contributors        */
  int zcontrib; /* (none) Bitmap of conributing catalogs to z       */
  int nstat; /* (none) Count of griz outliers rejected        */
  double J; /* (mag) 2MASS J magnitude          */
  double dJ; /* (mag) 2MASS J magnitude uncertainty         */
  double H; /* (mag) 2MASS H magnitude          */
  double dH; /* (mag) 2MASS H magnitude uncertainty         */
  double K; /* (mag) 2MASS K magnitude          */
  double dK; /* (mag) 2MASS K magnitude uncertainty         */  
} CELESTRONIMAGE,*PCELESTRONIMAGE;

typedef struct _gsccelestronimage {
	long long REFNumber; /* Encoded catalog reference number */
  long long objid; /* (none) Object ID           */
  int celestronNumber;     /* This is just the place in the .csv file set */
  double ra;           /* Right Ascension in degrees */
  double dec;          /* Declination in degrees */
  float Stdmag;        /* Blue magnitude */
  float color;         /* color (bmag - vmag) */
  float RaPM;          /* Right ascension proper motion in mas/yr */
  float DecPM;         /* Declination proper motion in mas/yr */
	float RaSigmaPM;     /* Error of the Right Ascension proper motion */
	float DecSigmaPM;    /* Error of the Declination proper motion */
  float rmag;          /* SDSS r color */
  float imag;          /* SDSS i color */
  float zmag;          /* SDSS z color */
  float gmagerr;       /* magnitude errors */
  float rmagerr;
  float imagerr;
  float zmagerr;
  char class;          /* class */
  char VFlag;          /* variable flag */
  char MAGFlag;        /* magnitude flag */
  char flag;           /* padding - used by check_photometry and formatgsc */
} GSCCELESTRONIMAGE,*PGSCCELESTRONIMAGE;
typedef struct _fileEntry {
  const char *fileName;
  char lastTab[MAX_BUFFER];
  PSTRINGINDEX pStringIndex;
  char *pTabTable; /* Allocation depends on the -t parameter */
  int fileIndex; /* Original index of this file */
  int maxLines; /* Maximum for this file */
  int numLines; /* Actual number of lines */
  int maxTabs;  /* Maximum tabs for this file */
  int tabCount; /* Actual number of tabs in this file */
  int maxLineLength;  /* Maximum line length for this file */
  int maxTabLength;   /* Maximum tab length for this file */
  int maxFieldLength; /* Maximum field length for this file */
  char errorFlag;      /* zero if no error, passCounter if an error */
  char separationChar; /* field separation character for this file */
} FILEENTRY,*PFILEENTRY;
#define FILE_CELESTRON_NT   0
#define FILE_CELESTRON_ST   1
#define FILE_FT70D_HOME     2
#define FILE_FT70D_MEMORIES 3
#define FILE_FT70D_PMS      4
#define FILE_FT70D_VFO      5
#define FILE_FT70D_SKIP     6
#define FILE_FT70D_WEATHER  7
#if 1   
#define MAX_FILE_TABLE_SIZE 1
#else
#define MAX_FILE_TABLE_SIZE 8
#endif
FILEENTRY fileTable[MAX_FILE_TABLE_SIZE] = {
#if 0
  {"/home/scanner/Pipeline/catalogs/celestron_st/Master.csv"                                          ,"",NULL,NULL,0,45928,0,17,0,0,0,0,0,','},
#endif
  {"/home/scanner/Pipeline/catalogs/celestron_st/celestron_st.txt"                                    ,"",NULL,NULL,1,  593,0, 7,0,0,0,0,0,';'}, 
#if 0
  {"/home/scanner/users/default/Misc/ham/yaesu_ft70-dr/FT-70D_saved_memories_Dec2021_12_02T082800.csv","",NULL,NULL,3,  900,0,52,0,0,0,0,0,','},
  {"/home/scanner/users/default/Misc/ham/yaesu_ft70-dr/FT-70D_saved_skip_Dec2021_12_02T082800.csv"    ,"",NULL,NULL,6,   99,0,26,0,0,0,0,0,','},
  {"/home/scanner/users/default/Misc/ham/yaesu_ft70-dr/FT-70D_saved_pms_Dec2021_12_02T082800.csv"     ,"",NULL,NULL,4,  100,0,27,0,0,0,0,0,','},
  {"/home/scanner/users/default/Misc/ham/yaesu_ft70-dr/FT-70D_saved_vfo_Dec2021_12_02T082800.csv"     ,"",NULL,NULL,5,    6,0,23,0,0,0,0,0,','},
  {"/home/scanner/users/default/Misc/ham/yaesu_ft70-dr/FT-70D_saved_home_Dec2021_12_02T082800.csv"    ,"",NULL,NULL,2,    6,0,25,0,0,0,0,0,','},
  {"/home/scanner/users/default/Misc/ham/yaesu_ft70-dr/FT-70D_saved_weather_Dec2021_12_02T082800.csv" ,"",NULL,NULL,7,   10,0,31,0,0,0,0,0,','},
#endif
};
int fileTableSize = sizeof(fileTable)/sizeof(FILEENTRY);
void ProcessCelestronLine(PFILEENTRY pFileEntry,char *output_name,FILE *output_handle);
void SaveOutputTableLine(int lineno,PFILEENTRY pFileEntry,char *output_name,FILE *output_handle);

/* Sort the records in decreasing declination */
int CELESTRONCompare(const void *first, const void *second) 
{
  double numberFirst = ((PGSCCELESTRONIMAGE )first)->dec;
  double numberSecond = ((PGSCCELESTRONIMAGE )second)->dec;
  if (numberFirst > numberSecond) {
    return(-1);
  } else if (numberFirst < numberSecond) {
    return(1);
  } else {
    return(0);
  }
}



void FlushOutputTable(FILE* output_handle,PGSCCELESTRONIMAGE pGscCelestronTable,int maxAllocOutput,int max_read_dec_bin,int* pCurOutputCount,int *pOutputWrittenCount,int *poutput_dec_bin,int *decBinCountsOutput,int finalFlag)
{
  int curOutputCount = *pCurOutputCount;
  int outputWrittenCount = *pOutputWrittenCount;
  int output_dec_bin = *poutput_dec_bin;
  int initial_dec_bin = output_dec_bin;
  int curIndex;
  int newOutputCount = 0;
  int dec_bin;
  printf("line %4d in FlushOutputTable\n",__LINE__);
  
  PGSCCELESTRONIMAGE pGscCelestronImage = NULL;  
  if (curOutputCount == 0) {
    return;
  }
  qsort((void*)pGscCelestronTable,curOutputCount,sizeof(GSCCELESTRONIMAGE),CELESTRONCompare);
  for (curIndex = (curOutputCount-1); curIndex >= 0; curIndex--) {
    pGscCelestronImage = &pGscCelestronTable[curIndex];
#ifndef NO_GETDECBIN
    dec_bin = GetDecBin(pGscBin,pGscCelestronImage->dec,"formatcelestron");
    if (dec_bin <= initial_dec_bin) {
      printf("ERROR: line %4d sort error dec_bin %d output_dec_bin %d, max_read_dec_bin %d\n",__LINE__,dec_bin,output_dec_bin,max_read_dec_bin);
    }
    if ((finalFlag == 0) && 
        (dec_bin > (max_read_dec_bin - RESIDENT_DEC_BINS))) { 
      break;
    }
    output_dec_bin = dec_bin;
    newOutputCount++;
    decBinCountsOutput[dec_bin]++;
#endif /* NO_GETDECBIN */
    fprintf(output_handle,"%lld\t%lld\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
            pGscCelestronImage->REFNumber,
            pGscCelestronImage->objid,
            pGscCelestronImage->ra,
            pGscCelestronImage->dec,
            pGscCelestronImage->Stdmag,
            pGscCelestronImage->color,
            pGscCelestronImage->class,
            pGscCelestronImage->VFlag,
            pGscCelestronImage->MAGFlag,
            pGscCelestronImage->RaPM,
            pGscCelestronImage->DecPM,
            pGscCelestronImage->RaSigmaPM,
            pGscCelestronImage->DecSigmaPM,
            pGscCelestronImage->rmag,
            pGscCelestronImage->imag,
            pGscCelestronImage->zmag,
            pGscCelestronImage->gmagerr,
            pGscCelestronImage->rmagerr,
            pGscCelestronImage->imagerr,
            pGscCelestronImage->zmagerr);
  }
  if (newOutputCount == 0) {
    printf("ERROR: line %4d failed to flush the buffer\n",__LINE__);
    exit(-1);
  }
  printf("line %4d FlushOutputTable wrote %d records from dec bin %d to %d curIndex %d\n",__LINE__,newOutputCount,initial_dec_bin,dec_bin-1,curIndex);
  *poutput_dec_bin = dec_bin -1;
  *pOutputWrittenCount = outputWrittenCount+newOutputCount;
  *pCurOutputCount = curIndex+1;
}

int main(int argc,char *argv[])
{
  char *argstr;
  char output_name[MAX_INPUT_NAME];
  char raw_name[MAX_INPUT_NAME];
  FILE* output_handle = NULL;
  FILE* raw_handle = NULL;
  int outputCount = 0;
  int colorReject = 0;
  int skipCount = 0;
  int rawSkipCount = 0; /* Count of raw file output skipped */
  int zeromagCount = 0;
  int rawSkipFlag; /* entry skipped for raw file output */
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int tabIndex;
  int tabIndex2;
  int tabCount2;
  char *pTab;
  int tabLength;
  int tabLengthIndex;
  int linecounter = 0;
  int index1;
  int index2;
  char *charPtr;
  char *charPtr2;
  char *charPtr1;
  char cmdchar;
  char tmpChar;
  char tmpChar2;
  int verbose = 0;
  int acceptNegative = 0;
  char *dotPtr = NULL;
  int compareResult;
  int result;
  int binNumber;
  int oldBinNumber = -1;
  int dec_bin;
  int old_dec_bin;
  double oldDec = -999;
  int skip_bin_index;
  int skip_bin_size;
  int oldDecBin = 0;
  int savedDecBin = -2;
  int ra_bin;
  int maxDecBin = 0;
  double maxRaError = 0;
  double maxDecError = 0;
  int maxWriteCount = 0;
  int totalBinCount = 0;
  FILE *input_handle;
  char *celestronMeasurementText[CELESTRON_MAXMAG] = {"bp","rp","g","r","i","z"};
  double newrms;
  double scalerms = 0;
  double maxrms2 = 0;
  int intrms;
  int intrms2;
  long long max_objid;
  long long min_objid;
  unsigned long long objidbits[MAX_MODULUS]; 
  unsigned long long modobjidbits[MAX_MODULUS];
  int modindex;
  unsigned long long temp_objid;
  unsigned long long modulus;
  int old_bin_index = -1;
  PFILEENTRY pFileEntry;
  double factor;
  double dec2;
  double ra2;
  int flipFlag; /* Require 12 hour telescope flip when crossing the pole */
  int dec_bin2;
  int globalTotalTabs = 0;
  int globalMaxTabLength = 0;
  off_t tempStorageSize = 0;
  int maxTabLength =  MAX_TAB_LENGTH;
  int tabLengthTable[MAX_TAB_LENGTH];

  memset(tabLengthTable,0,sizeof(tabLengthTable));
  memset(objidbits,0,sizeof(objidbits)); 
  memset(modobjidbits,0,sizeof(modobjidbits));




  int fileIndex;

  int input_nrecs = 0;
  CELESTRONIMAGE input_record;
  PCELESTRONIMAGE pInput = &input_record;
  PGSCCELESTRONIMAGE pGscCelestronImage = NULL;
  PGSCCELESTRONIMAGE pGscCelestronTable = NULL;
  int passCounter;

  int tempmag;


  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int celestronCount[CELESTRON_MAXMAG];
  int celestronDual[CELESTRON_MAXMAG][CELESTRON_MAXMAG];
  int magcount[CELESTRON_MAXMAG][MAGBINS+5];
  int rmscount[RMSBINS+5];
  int rmscountTotal[RMSBINS+5][MAGBINS+5];
  
  int celestronNoCQCount = 0;
  int celestronHeaderCount = 0;
  int maxcqLength = 0;
  char maxcq[MAX_REF];
  double deltacolor;
  double maxdeltacolor = 0.0;
  int maxcelestroncolorid;
  char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
  int lineLen;
  int maxLineLen = 0;
  char *inBuffer;
  int numLines;
  int numFileLines;
  int nvals;
  double rmsCutoff = CELESTRON_MAXRMS;
  double value1;
  double rms1;
  double value2;
  double rms2;
  long long REFNumber;
  int skipFlag; /* Skip this entry because of zero mag or too high rms */
  int skipFlagBase;  /* Applies to the first of a magnitude pair */
  int zeroBmagCount = 0;
  int zeroVmagCount = 0;
  int negZeroCount = 0;
  int posZeroCount = 0;
  double maxpmdec = 0;
  double minpmdec = 0;
  double mindec = 0;
  double maxdec = 0;
  double minra = 0;
  double maxra = 0;
  int *decBinCounts = NULL;
  int *decBinCountsInitial = NULL;
  int *decBinCountsFinal = NULL;
  int *decBinCountsOutput = NULL;
  int *skipBinCounts = NULL;
  int *skipBinCountsTotal = NULL;
  int maxDecBinCount = 0;
  int maxDecBinValue = -1;
  int expectedRecords = 0;
  int outputIndex;
  int curOutputCount = 0;
  int maxAllocOutput = 0;
  int outputWrittenCount = 0;
  long long maxAllocBytes;
  int maxSkipLines = 0;
  int skipLines;
  int allocOutputLines;
  int output_dec_bin = -1; /* Last bin to be written */
  int max_read_dec_bin = -1;  /* Highest bin read in pass 1 */ 

  memset(celestronCount,0,sizeof(celestronCount));
  memset(celestronDual,0,sizeof(celestronDual));
  memset(magcount,0,sizeof(magcount));
  memset(rmscount,0,sizeof(rmscount));
  memset(rmscountTotal,0,sizeof(rmscountTotal));
#ifndef NO_GETDECBIN
  decBinCounts = (int *)calloc(pGscBin->dec_bins,sizeof(int));
  decBinCountsInitial = (int *)calloc(pGscBin->dec_bins,sizeof(int));
  decBinCountsFinal = (int *)calloc(pGscBin->dec_bins,sizeof(int));
  decBinCountsOutput = (int *)calloc(pGscBin->dec_bins,sizeof(int));
  skipBinCounts = (int *)calloc((2*pGscBin->dec_bins)+2,sizeof(int));
  skipBinCountsTotal = (int *)calloc((2*pGscBin->dec_bins)+2,sizeof(int));
#endif /* NO_GETDECBIN */
  /* Loop through the arguments */
  output_name[0] = 0;
	raw_name[0] = 0;

  memset(magcount,0,sizeof(magcount));
  

  maxcq[0] = 0;
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


        case 'c': /* RMS cutoff */
        case 'C':
          argc--;
          nvals = sscanf(*++argv,"%lf",&rmsCutoff);
          if (nvals != 1) {
            fprintf(stderr,"ERROR: Can not decode the rmsCutoff\n");
            errorFlag = 1;
          }
          break;

        case 't': /* storage cutoff for preallocated per-tab storage */
        case 'T':
          argc--;
          nvals = sscanf(*++argv,"%d",&maxTabLength);
          if (nvals != 1) {
            fprintf(stderr,"ERROR: Can not decode the maxTabLength\n");
            errorFlag = 1;
          } else {
            if (maxTabLength < 0) {
              maxTabLength = 0;
              printf("WARNING: line %4d too low.  maxTabLength set to %d\n",__LINE__,maxTabLength); 
            } else if (maxTabLength > MAX_TAB_LENGTH) {
              maxTabLength = MAX_TAB_LENGTH;
              printf("WARNING: line %4d too high maxTabLength set to %d\n",__LINE__,maxTabLength); 
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

        case 'r': /* raw file name */
        case 'R':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(raw_name,*++argv,MAX_INPUT_NAME-2);
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

  /* Attempt to open the list */
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

  if (raw_name[0] != 0) {
		raw_handle = fopen(raw_name,"wt");
    if (raw_handle == NULL) {
      errorFlag = 1;
      printf("Could not open raw file %s\n",raw_name);
    } else {
			fprintf(raw_handle,"ra\tdec\tStdmag\trms\n");
			fprintf(raw_handle,"--\t---\t------\t---\n");

		}
  }

  scalerms = rmsCutoff * sqrt(2.0);
  maxrms2 = sqr(scalerms);

  if (errorFlag) {
    printf("Usage: formatcelestron -o {output name} [-r <raw file>] [-v][-n][-s][-c {cutoff}]\n");
		printf("       where -r <filename> is a raw table of RA and DEC for plotting \n");
    printf("       where -v is the verbose flag\n");
    printf("       where -c {cutoff} gives the RMS cutoff values\n");
    printf("       where -t {maxTabLength} uses a separate array to store large tab entries\n");

    return(-1);
  }

  printf("formatcelestron of %s %s Output Filename %s Raw filename %s, rmsCutoff %f, maxTabLength %d\n",
         __DATE__,__TIME__,output_name,raw_name,rmsCutoff,maxTabLength);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of CELESTRONIMAGE is %d Size of GSCCELESTRONIMAGE is %d Size of GSCIMAGE is %d\n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(CELESTRONIMAGE),sizeof(GSCCELESTRONIMAGE),sizeof(GSCIMAGE));
 

  time(&startTime);
  
#if 0  
  fprintf(output_handle,"REFNumber\tobjid\tra\tdec\tStdmag\tcolor\tclass\tVFlag\tMAGFlag\tRaPM\tDecPM\tRaSigmaPM\tDecSigmaPM\trmag\timag\tzmag\tgmagerr\trmagerr\timagerr\tzmagerr\n");
  fprintf(output_handle,"---------\t-----\t--\t---\t------\t-----\t-----\t-----\t-------\t----\t-----\t---------\t----------\t----\t----\t----\t-------\t-------\t-------\t-------\n");
#endif

  for (fileIndex = 0; fileIndex < fileTableSize; fileIndex++) {
    pFileEntry = &fileTable[fileIndex];
    pFileEntry->fileIndex = fileIndex;
    printf("line %4d: maxLines %d maxTabs %d separationChar '%c' for '%s'\n",
           __LINE__,
           pFileEntry->maxLines,
           pFileEntry->maxTabs,
           pFileEntry->separationChar,
           pFileEntry->fileName);
    fflush(stdout);
    AllocateStringTable(&pFileEntry->pStringIndex,(MAX_TAB*pFileEntry->maxTabs),__LINE__,(char *)pFileEntry->fileName);
    for (tabIndex = 0; tabIndex < MAX_TAB; tabIndex++) {
      pFileEntry->pTabTable = pFileEntry->pStringIndex->pStringEntry;
      pTab = (char*)&pFileEntry->pTabTable[tabIndex];
      memset(pTab,0,sizeof(stringBuffer));
#if 0
      strcpy(stringBuffer,"0BCDEFGHIJ1BCDEFGHIJ2BCDEFGHIJ3BCDEFGHIJ4BCDEFGHIJ5BCDEFGHIJ");
#endif
      WriteStringTable(pFileEntry->pStringIndex,pTab,MAX_TAB,__LINE__,(char *)pFileEntry->fileName);
    }
    
  }


  input_nrecs = 0;
  for (fileIndex = 0; fileIndex < fileTableSize; fileIndex++) {
  
    pFileEntry = &fileTable[fileIndex];
    
    printf("line %4d: Opening %s\n",__LINE__,pFileEntry->fileName);
    memset(&(tabTable[0][0]),0,sizeof(tabTable));
    for (tabIndex = 0; tabIndex <= MAX_TAB; tabIndex++) {
      tabIndexTable[tabIndex] = -1;
    }
    for (passCounter = 1; passCounter <= 3; passCounter++) {
      input_handle = fopen(pFileEntry->fileName,"rt");
      if (input_handle == NULL) {
        printf("ERROR: line %4d passCounter %d failed to open %s\n",__LINE__,passCounter,pFileEntry->fileName);
        exit(-1);
      }
      pFileEntry->numLines = 0;
      while(1) {
        lineLen = 0;
        tabIndex = 0;
        memset(pInput,0,sizeof(CELESTRONIMAGE));

        inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
        if (inBuffer == NULL) {
          break;
        }
        lineLen = strlen(inBuffer);
        if (lineLen >= MAX_BUFFER) {
          printf("ERROR line %4d: MAX_BUFFER exceeded by file '%s'\n",__LINE__,pFileEntry->fileName);
          exit(-1);
        }

        if (pFileEntry->numLines >= MAX_LINES) {
          printf("ERROR line %4d: MAX_LINES exceeded by file '%s'\n",__LINE__,pFileEntry->fileName);
          exit(-1);
        }
        if (pFileEntry->numLines >= pFileEntry->maxLines) {
          printf("ERROR line %4d: maxLines %d exceeded by file '%s'\n",__LINE__,pFileEntry->maxLines,pFileEntry->fileName);
        }
        pFileEntry->numLines++;
        lineLen = strlen(inBuffer);
        if (lineLen > pFileEntry->maxLineLength) {
          pFileEntry->maxLineLength = lineLen;
        }
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
        strcpy(copyLine,inLine);
#if 1
        printf("DEBUG: line %4d tabIndex %d lineLen %d\n",
               __LINE__,
               tabIndex,
               lineLen);
        printf("copyLine %s\n",copyLine);
#endif
#if 0
        printf("DEBUG: line %4d tabIndex %d lineLen %d pFileEntry->maxTabs %d pFileEntry->tabCount %d\n",
               __LINE__,
               tabIndex,
               lineLen,
               pFileEntry->maxTabs,
               pFileEntry->tabCount);
#endif
#if 0
        printf("DEBUG: line %4d tabIndex %d lineLen %d pFileEntry->maxTabs %d pFileEntry->tabCount %d separationChar '%c' for file '%s' \n",
               __LINE__,
               tabIndex,
               lineLen,
               pFileEntry->maxTabs,
               pFileEntry->tabCount,
               pFileEntry->separationChar,
               pFileEntry->fileName);
#endif
#if 0
        printf("DEBUG: line %4d tabIndex %d lineLen %d pFileEntry->maxTabs %d pFileEntry->tabCount %d separationChar '%c' for file '%s' \n",
               __LINE__,
               tabIndex,
               lineLen,
               pFileEntry->maxTabs,
               pFileEntry->tabCount,
               pFileEntry->separationChar,
               pFileEntry->fileName);
#endif
        if ((passCounter == 1) && (pFileEntry->errorFlag == 0)) {
          printf("DEBUG: line %4d passCounter %d for file '%s' \n",
                 __LINE__,
                 passCounter,
                 pFileEntry->fileName);
          /* Now count the tabs in each line and keep track of the maximum record size and the very last record parsed */
          pFileEntry->tabCount = 0;
          if (lineLen > 0) {
            globalTotalTabs++;
          }
          for (tabIndex = 0; tabIndex < lineLen; tabIndex++) {
            if (inBuffer[tabIndex] == pFileEntry->separationChar) {
              pTab = (char*)&pFileEntry->pTabTable[tabIndex];
              /* This is a tab */
#if 0
              printf("DEBUG: line %4d tabIndex %d lineLen %d pFileEntry->maxTabs %d pFileEntry->tabCount %d separationChar '%c' for file '%s' \n",
                     __LINE__,
                     tabIndex,
                     lineLen,
                     pFileEntry->maxTabs,
                     pFileEntry->tabCount,
                     pFileEntry->separationChar,
                     pFileEntry->fileName);
#endif
              inBuffer[tabIndex] = 0;
              tabIndexTable[pFileEntry->tabCount] = tabIndex;
              if (pFileEntry->tabCount == 0) {
                strcpy(&tabTable[pFileEntry->tabCount][0],&inBuffer[0]);
                strcpy(&tabTable[pFileEntry->tabCount+1][0],&inBuffer[tabIndex+1]); /* Assume that there is only one tab */
              } else {
                /* strcpy(&tabTable[pFileEntry->tabCount][0],&inBuffer[tabIndex+1]); */
                strcpy(&tabTable[pFileEntry->tabCount+1][0],&inBuffer[tabIndex+1]); /* Assume that this is the last tab. If not, it will be over-written */
              }
#if 0 /* Tue Dec 14 20:08:42 EST 2021 */
              for (tabCount2 = 0; tabCount2 < pFileEntry->tabCount; tabCount2++) {
                printf("DEBUG: line %4d tC2 %3d tabIndexTable[tC2] %3d tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tC2] %3d tabTable[tC2][0] '%s'\n",
                       __LINE__,tabCount2,tabIndexTable[tabCount2],tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tabCount2],&tabTable[tabCount2][0]);
              }
#endif
#if 0         /* Does not work */
              for (tabCount2 = 1; tabCount2 < pFileEntry->tabCount; tabCount2++) {
                tabTable[tabCount2][tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tabCount2]+1] = 0;
              }
#endif
              if (pFileEntry->tabCount > 1) {
                tabIndex2 = 0;
                while (1) {
                  tmpChar2 = tabTable[pFileEntry->tabCount][tabIndex2];
                  if (tmpChar2 == 0) {
                    break;
                  }
                  if (tmpChar2 ==  pFileEntry->separationChar) {
                    tabTable[pFileEntry->tabCount][tabIndex2] = 0;
                    break;
                  }
                  tabIndex2++;
                }
              }

#if 0 /* Tue Dec 14 20:08:42 EST 2021 */
              for (tabCount2 = 0; tabCount2 < pFileEntry->tabCount; tabCount2++) {
                printf("DEBUG: line %4d tC2 %3d tabIndexTable[tC2] %3d tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tC2] %3d tabTable[tC2][0] '%s'\n",
                       __LINE__,tabCount2,tabIndexTable[tabCount2],tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tabCount2],&tabTable[tabCount2][0]);
              }
#endif
#if 0 /* Tue Dec 14 20:08:42 EST 2021 */
              if (pFileEntry->numLines == 3) {
                printf("DEBUG: line %4d early exit\n",__LINE__); 
                exit(-1);
              } else {
                printf("\n");
              }
#endif
              pFileEntry->tabCount++;
              globalTotalTabs++;
              if (pFileEntry->tabCount >= MAX_TAB) {
#if 0
                printf("DEBUG: line %4d tabIndex %d lineLen %d pFileEntry->maxTabs %d pFileEntry->tabCount %d separationChar '%c' for file '%s' \n",
                       __LINE__,
                       tabIndex,
                       lineLen,
                       pFileEntry->maxTabs,
                       pFileEntry->tabCount,
                       pFileEntry->separationChar,
                       pFileEntry->fileName);
#endif
                printf("ERROR: line %4d too many tabs\n",__LINE__);
                pFileEntry->errorFlag = 1;

              } else {
                if (((pFileEntry->tabCount > pFileEntry->maxTabs) && (pFileEntry->fileIndex != FILE_CELESTRON_NT))) {
                  printf("ERROR: line %4d  too many tabs, got %d expect %d line %d in file '%s'\n",
                         __LINE__,
                         pFileEntry->tabCount,
                         pFileEntry->maxTabs,
                         pFileEntry->numLines,
                         pFileEntry->fileName);
#if 1
                  printf("ERROR: line %4d copyLine '%s'\n",__LINE__,copyLine);
#endif
                  pFileEntry->errorFlag = 1;
                }
#if 0
                printf("DEBUG: line %4d tabIndex %d lineLen %d pFileEntry->maxTabs %d pFileEntry->tabCount %d separationChar '%c' for file '%s' \n",
                       __LINE__,
                       tabIndex,
                       lineLen,
                       pFileEntry->maxTabs,
                       pFileEntry->tabCount,
                       pFileEntry->separationChar,
                       pFileEntry->fileName);
#endif
          
              

              }

            }

          }

#if 0
          printf("DEBUG: line %4d tabIndex %d lineLen %d pFileEntry->maxTabs %d pFileEntry->tabCount %d separationChar '%c' for file '%s' \n",
                 __LINE__,
                 tabIndex,
                 lineLen,
                 pFileEntry->maxTabs,
                 pFileEntry->tabCount,
                 pFileEntry->separationChar,
                 pFileEntry->fileName);
#endif

          if (pFileEntry->tabCount < pFileEntry->maxTabs) {
            printf("ERROR: line %4d  not enough tabs, got %d expect %d line %d in file '%s'\n",
                   __LINE__,
                   pFileEntry->tabCount,
                   pFileEntry->maxTabs,
                   pFileEntry->numLines,
                   pFileEntry->fileName);
#if 1
            printf("ERROR: line %4d copyLine '%s'\n",__LINE__,copyLine);
#endif
            pFileEntry->errorFlag = 1;
          } 
          
#if 0 /* Tue Dec 14 20:08:42 EST 2021 */
          for (tabCount2 = 0; tabCount2 < pFileEntry->tabCount; tabCount2++) {
            printf("DEBUG: line %4d tC2 %3d tabIndexTable[tC2] %3d tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tC2] %3d tabTable[tC2][0] '%s'\n",
                   __LINE__,tabCount2,tabIndexTable[tabCount2],tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tabCount2],&tabTable[tabCount2][0]);
          }
#endif
          printf("DEBUG: line %4d\n",__LINE__);
        } else if ((passCounter == 2)  && (pFileEntry->errorFlag == 0)){
          printf("DEBUG: line %4d passCounter %d for file '%s' \n",
                 __LINE__,
                 passCounter,
                 pFileEntry->fileName);
#if 1 /* Tue Dec 14 20:08:42 EST 2021 */
          for (tabCount2 = 0; tabCount2 < pFileEntry->tabCount; tabCount2++) {
            printf("DEBUG: line %4d tabCount %d  string: '%s'\n",
                   __LINE__,tabCount2,&(tabTable[tabCount2][0]));
            fflush(stdout);
          }
#endif
          for (tabIndex = 0; tabIndex < MAX_TAB; tabIndex++) {
            pTab = &tabTable[tabIndex][0];
            /* Trim off any leading space or quotation */
            tabLength = strlen(pTab);
            while (tabLength > 0) {
              tmpChar = *pTab;
              if ((tmpChar == 34) || (tmpChar == ' ')) {
                pTab++;
                tabLength--;
                tmpChar = *pTab;
              } else {
                break;
              }
            }
            /* Trim off trailing spaces or quotations */
            while(tabLength > 0) {
              tmpChar = pTab[tabLength-1];
              if ((tmpChar == 34) || (tmpChar == ' ')) {
                pTab--;
                tabLength--;
                tmpChar = *pTab;
              } else {
                break;
              }
            }
#if 1 /* Tue Dec 14 20:08:42 EST 2021 */
            for (tabCount2 = 0; tabCount2 < pFileEntry->tabCount; tabCount2++) {
              printf("DEBUG: line %4d tC2 %3d tabIndexTable[tC2] %3d tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tC2] %3d tabTable[tC2][0] '%s'\n",
                     __LINE__,tabCount2,tabIndexTable[tabCount2],tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tabCount2],&tabTable[tabCount2][0]);
            }
            printf("DEBUG: line %4d\n",__LINE__);
            SaveOutputTableLine(__LINE__,pFileEntry,output_name,output_handle);
#endif
            if (tabLength < MAX_TAB_LENGTH) {                   
              tabLengthTable[tabLength]++;
            } else {
              tabLengthTable[(MAX_TAB_LENGTH-1)]++;
            }

            if (tabLength > pFileEntry->maxTabLength) {
              pFileEntry->maxTabLength = tabLength;
            }
            if (globalMaxTabLength < tabLength) {
              globalMaxTabLength = tabLength;
            }
            if ((tabLength >= MAX_TAB_LENGTH) || (tabLength >= pFileEntry->maxTabLength)) {
              /* If the tab is too big, allocate more space from the pool */
#if 0 /* Tue Dec 14 20:08:42 EST 2021 */
              for (tabCount2 = 0; tabCount2 < pFileEntry->tabCount; tabCount2++) {
                printf("DEBUG: line %4d tC2 %3d tabIndexTable[tC2] %3d tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tC2] %3d tabTable[tC2][0] '%s'\n",
                       __LINE__,tabCount2,tabIndexTable[tabCount2],tabIndexTable[pFileEntry->tabCount]-tabIndexTable[tabCount2],&tabTable[tabCount2][0]);
              }
#endif
              AllocateStringTable(&pFileEntry->pStringIndex,tabLength,__LINE__,(char *)pFileEntry->fileName);

#if 1 /* Tue Dec 14 20:08:42 EST 2021 */
              for (tabCount2 = 0; tabCount2 < pFileEntry->tabCount; tabCount2++) {
                printf("DEBUG: line %4d tabCount %d  string: '%s'\n",
                       __LINE__,tabCount2,&(tabTable[tabCount2][0]));
                fflush(stdout);
              }
#endif
             
              WriteStringTable(pFileEntry->pStringIndex,pTab,tabLength,__LINE__,(char *)pFileEntry->fileName);
            } else {
              /* Here the tab fits, write, but perform one last sanity check */
#if 0
              DumpStringTable(__LINE__,(char *)pFileEntry->fileName);
#endif
#if 0
              GetStringTable(pFileEntry->pStringIndex,pTab,tabLength,__LINE__,(char *)pFileEntry->fileName);
#endif
              WriteStringTable(pFileEntry->pStringIndex,pTab,tabLength,__LINE__,(char *)pFileEntry->fileName);
                   
            }
          }

        } else if (pFileEntry->errorFlag == 0 ){
          printf("DEBUG: line %4d passCounter %d for file '%s' \n",
                 __LINE__,
                 passCounter,
                 pFileEntry->fileName);
          ProcessCelestronLine(pFileEntry,output_name,output_handle);
        }
      }
#if 0
      printf("DEBUG: line %4d tabIndex %d lineLen %d pFileEntry->maxTabs %d pFileEntry->tabCount %d separationChar '%c' for file '%s' \n",
             __LINE__,
             tabIndex,
             lineLen,
             pFileEntry->maxTabs,
             pFileEntry->tabCount,
             pFileEntry->separationChar,
             pFileEntry->fileName);
#endif
#if 0
      if (errorFlag != 0) {
        printf("line %4d errorFlag %d EARLY EXIT\n",__LINE__,errorFlag);
        exit(-1);
      }
#endif
      time(&curTime);
      curTime -= startTime;
      printf("line %4d: numLines   %5d  maxLineLength %3d for '%s'\n",__LINE__,pFileEntry->numLines,pFileEntry->maxLineLength,pFileEntry->fileName);
      printf("line %4d: maxTabLength %3d              lastTab '%s' for '%s'\n",__LINE__,pFileEntry->maxTabLength,pFileEntry->lastTab,pFileEntry->fileName);
#if 0
      for (dec_bin = 0; dec_bin < pGscBin->dec_bins; dec_bin++) {
        if (decBinCounts[dec_bin] > 0) {
          printf("dec_bin %d Declination %8.3f count %d\n",dec_bin,(pGscBin->bin_size*dec_bin)-90.0,decBinCounts[dec_bin]);
        }
      }
#endif
#ifndef NO_GETDECBIN
      memset(decBinCounts,0,pGscBin->dec_bins*sizeof(int));
#endif /* NO_GETDECBIN */
#if 0
      for (skip_bin_index = 0; skip_bin_index < (2*pGscBin->dec_bins)+2; skip_bin_index++) {
        if (skipBinCounts[skip_bin_index] > 0) {
          printf("skip_bin_index %d count %d\n",skip_bin_index-pGscBin->dec_bins,skipBinCounts[skip_bin_index]);
        }
      }
#endif
#ifndef NO_GETDECBIN
      memset(skipBinCounts,0,((2*pGscBin->dec_bins)+2)*sizeof(int));
#endif /* NO_GETDECBIN */
      printf("DEBUG: line %4d passCounter %d closing file '%s' \n",
             __LINE__,
             passCounter,
             pFileEntry->fileName);
      fclose(input_handle);
      input_handle = NULL;
    } /* passCounter loop */
    input_nrecs += numFileLines;
    printf("line %4d entering FlushOutputTable\n",__LINE__);
    FlushOutputTable(output_handle,pGscCelestronTable,maxAllocOutput,max_read_dec_bin,&curOutputCount,&outputWrittenCount,&output_dec_bin,decBinCountsOutput,0);
  }

  printf("line %4d: globalTotalTabs %d globalMaxTabLength %d for all files\n",__LINE__,globalTotalTabs,globalMaxTabLength);
  tempStorageSize = 0;
  for (tabLengthIndex = 0; tabLengthIndex < MAX_TAB_LENGTH; tabLengthIndex++) {
    if (tabLengthTable[tabLengthIndex] != 0) {
      tempStorageSize += tabLengthTable[tabLengthIndex] * tabLengthIndex; /* Count actual bytes */
      if (tabLengthIndex > maxTabLength) {
        tempStorageSize += maxTabLength; /* Need a new allocation because it does not fit in the pre-allocated table  */
      }
      if ((((tabLengthIndex+1) % TAB_BIN_SIZE) == 0) ||
          (tabLengthIndex == (MAX_TAB_LENGTH-1))) {
        /* Time to print a row */
        if (tabLengthIndex > maxTabLength) {
          globalAllocatedStorageSize += tempStorageSize;
        } 

        printf("line %4d tabLengthIndex %3d tabLengthTable[tabLengthIndex] %8d Allocates: %8d Static: %8d Percentage %5.1f\n",
               __LINE__,
               tabLengthIndex,
               tabLengthTable[tabLengthIndex],
               globalAllocatedStorageSize,
               globalStaticStorageSize,
               ((100.*globalAllocatedStorageSize)/(1.0*globalStaticStorageSize)));
        tempStorageSize = 0;
      }
    }
  }
#if 0
  printf("ERROR: line %4d ERROR: EARLY EXIT\n",__LINE__);
  exit(-1);
#endif

  printf("line %4d entering FlushOutputTable\n",__LINE__);
 
#ifndef NO_GETDECBIN
  FlushOutputTable(output_handle,pGscCelestronTable,maxAllocOutput,max_read_dec_bin,&curOutputCount,&outputWrittenCount,&output_dec_bin,decBinCountsOutput,1);
#endif /* NO_GETDECBIN */

  printf("min_pmdec %f max_pmdec %f mindec %f maxdec %f minra %f maxra %f\n",minpmdec,maxpmdec,mindec,maxdec,minra,maxra);
  
  printf("min_objid %lld, max_objid %lld\n",min_objid,max_objid);
  modulus = 1;
  for (modindex = 0; modindex < MAX_MODULUS; modindex++) {
    objidbits[modindex] |= temp_objid/modulus;
    modobjidbits[modindex] |= temp_objid % modulus;
    printf("modindex %2d modulus %20lld objidbits 0x%016llx modobjidbits 0x%016llx\n",modindex,modulus,objidbits[modindex],modobjidbits[modindex]);
    modulus = modulus * MODULUS_BASE;
  }



  if (raw_handle != NULL) {
    fclose(raw_handle);
  }
  fclose(output_handle);
 
  time(&curTime);
  curTime -= startTime;

  printf("Measurement:");
  for (index1 = 0; index1 < CELESTRON_MAXMAG; index1++) {
    printf("%10s",celestronMeasurementText[index1]);
  }
  printf("\n");

  printf("Count:      ");
  for (index1 = 0; index1 < CELESTRON_MAXMAG; index1++) {
    printf("%10d",celestronCount[index1]);
  }
  printf("\n");

 
  for (index1 = 0; index1 < CELESTRON_MAXMAG; index1++) {
    printf("%10s  ",celestronMeasurementText[index1]);
    for (index2 = 0; index2 < CELESTRON_MAXMAG; index2++) {
      printf("%10d",celestronDual[index1][index2]);
    }
    printf("\n");
  }





  for (tempmag = 0; tempmag <= MAGBINS+2; tempmag++) {
    printf("gmag %2d     ",tempmag);
    for (index1 = 0; index1 < CELESTRON_MAXMAG; index1++) {
      printf("%10d",magcount[index1][tempmag]);
    }
    printf("\n");
  }    

  for (intrms = 0; intrms <= RMSBINS+2; intrms++) {
    if (rmscount[intrms] > 0) {
      printf("g-r rms  (mag) %.2f, count %6d\n",(CELESTRON_MAXRMS*intrms)/(1.0*RMSBINS),rmscount[intrms]);
    }
  }
#ifdef PLOT_RMS_GMAG
  printf("gmag\trms\tcount\n");
  printf("----\t---\t-----\n");
  

  for (intrms = 0; intrms <= (RMSBINS+2); intrms++) {
    for (tempmag = 0; tempmag <= MAGBINS+2; tempmag++) {
      if (rmscountTotal[intrms][tempmag] != 0) {
        printf("%f\t%f\t%d\n",1.0*tempmag,(CELESTRON_MAXRMS*intrms)/(1.0*RMSBINS),rmscountTotal[intrms][tempmag]);
      }

    }
  }
#endif /* PLOT_RMS_GMAG */
#ifndef NO_GETDECBIN
  for (dec_bin = 0; dec_bin < pGscBin->dec_bins; dec_bin++) {
    if ((decBinCountsInitial[dec_bin] > 0) && (decBinCountsFinal[dec_bin] > 0) && (decBinCountsOutput[dec_bin] > 0)) {
      if (decBinCountsInitial[dec_bin] > maxDecBinCount) {
        maxDecBinCount = decBinCountsInitial[dec_bin];
        maxDecBinValue = dec_bin;
      }
#if 1
      printf("dec_bin %5d Declination %8.3f count initial: %6d final %6d output %6d\n",dec_bin,(pGscBin->bin_size*dec_bin)-90.0,decBinCountsInitial[dec_bin],decBinCountsFinal[dec_bin],decBinCountsOutput[dec_bin]);
#endif
    }
  }
  printf("maxDecBinValue %d maxDecBinCount %d total\n",maxDecBinValue,maxDecBinCount);
#if 0
  for (skip_bin_index = 0; skip_bin_index < (2*pGscBin->dec_bins)+2; skip_bin_index++) {
    if (skipBinCountsTotal[skip_bin_index] > 0) {
      printf("skip_bin_index %d count %d total\n",skip_bin_index-pGscBin->dec_bins,skipBinCountsTotal[skip_bin_index]);
    }
  }
#endif
  if (pGscCelestronTable != NULL) {
    free(pGscCelestronTable);
  }
#endif /* NO_GETDECBIN */

  time(&curTime);
  curTime -= startTime;

  printf("zeroBmagCount %d zeroVmagCount %d negZeroCount %d posZeroCount %d\n",zeroBmagCount,zeroVmagCount,negZeroCount,posZeroCount);
  printf("Execution Time: %d seconds; rejected %d rejected because of color %d stars written %d\n",curTime,skipCount,colorReject,outputCount);
  printf("outputWrittenCount %d\n",outputWrittenCount);
  printf("rawSkipCount %d\n",rawSkipCount);

  return(EXIT_SUCCESS);

}
void ProcessCelestronLine(PFILEENTRY pFileEntry,char *output_name,FILE* output_handle) {
#if 0
  nvals = sscanf(inBuffer,"%lld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%lf,%d,%lf,%lf,%lf,%d,%lf,%lf,%lf,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf",
                 &pInput->objid, /* (none) Object ID           */
                 &pInput->RA, /* (degrees) Right ascension from Gaia DR2, J2000, epoch 2015.5     */
                 &pInput->Dec, /* (degrees) Declination from Gaia DR2, J2000, epoch 2015.5      */
                 &pInput->plx, /* (mas) Parallax from Gaia DR2         */
                 &pInput->dplx, /* (mas) Parallax uncertainty from Gaia DR2        */
                 &pInput->pmra, /* (mas/yr) Proper motion in right ascension from Gaia DR2     */
                 &pInput->dpmra, /* (mas/yr) Proper motion uncertainty in right ascension       */
                 &pInput->pmdec, /* (mas/yr) Proper motion in declination from Gaia DR2      */
                 &pInput->dpmdec, /* (mas/yr) Proper motion uncertainty in declination        */
                 &pInput->Gaia, /* (mag) Gaia G magnitude          */
                 &pInput->dGaia, /* (mag) Gaia G magnitude uncertainty         */
                 &pInput->BP, /* (mag) Gaia G_bp magnitude          */
                 &pInput->dBP, /* (mag) Gaia G_bp magnitude uncertainty         */
                 &pInput->RP, /* (mag) Gaia G_rp magnitude          */
                 &pInput->dRP, /* (mag) Gaia G_rp magnitude uncertainty         */
                 &pInput->Teff, /* (K) Gaia stellar effective temperature         */
                 &pInput->AGaia, /* (mag) Gaia estimate of G-band extinction for this star     */
                 &pInput->dupvar, /* (none) Gaia variability and duplicate flags, 0/1/2 for "CONSTANT"/"VARIABLE"/"NOT AVAILABLE" + 4*DUPLICATE  */
                 &pInput->Ag, /* (mag) SFD estimate of total g-band extinction       */
                 &pInput->rp1, /* (arcsec) Radius where cummulative G flux exceeds 0.1 x this star   */
                 &pInput->r1, /* (arcsec) Radius where cummulative G flux exceeds 1.0 x this star   */
                 &pInput->r10, /* (arcsec) Radius where cummulative G flux exceeds 10.0 x this star   */
                 &pInput->g, /* (mag) PanSTARRS g magnitude          */
                 &pInput->dg, /* (mag) PanSTARRS g magnitude uncertainty         */
                 &pInput->gchi, /* (none) chi^2 / DOF for contributors        */
                 &pInput->gcontrib, /* (none) Bitmap of conributing catalogs to g       */
                 &pInput->r, /* (mag) PanSTARRS r magnitude          */
                 &pInput->dr, /* (mag) PanSTARRS r magnitude uncertainty         */
                 &pInput->rchi, /* (none) chi^2 / DOF for contributors        */
                 &pInput->rcontrib, /* (none) Bitmap of conributing catalogs to r       */
                 &pInput->i, /* (mag) PanSTARRS i magnitude          */
                 &pInput->di, /* (mag) PanSTARRS i magnitude uncertainty         */
                 &pInput->ichi, /* (none) chi^2 / DOF for contributors        */
                 &pInput->icontrib, /* (none) Bitmap of conributing catalogs to i       */
                 &pInput->z, /* (mag) PanSTARRS z magnitude          */
                 &pInput->dz, /* (mag) PanSTARRS z magnitude uncertainty         */
                 &pInput->zchi, /* (none) chi^2 / DOF for contributors        */
                 &pInput->zcontrib, /* (none) Bitmap of conributing catalogs to z       */
                 &pInput->nstat, /* (none) Count of griz outliers rejected        */
                 &pInput->J, /* (mag) 2MASS J magnitude          */
                 &pInput->dJ, /* (mag) 2MASS J magnitude uncertainty         */
                 &pInput->H, /* (mag) 2MASS H magnitude          */
                 &pInput->dH, /* (mag) 2MASS H magnitude uncertainty         */
                 &pInput->K, /* (mag) 2MASS K magnitude          */
                 &pInput->dK); /* (mag) 2MASS K magnitude uncertainty         */  
  if (nvals != 45) {
    printf("ERROR: line %4d nvals is %d in line %4d\n",__LINE__,nvals,numLines);
  } else {
    if (numFileLines == 1) {
      mindec = pInput->Dec;
      maxdec = pInput->Dec;
      minra = pInput->RA;
      maxra = pInput->RA;
    } else {
      if (pInput->Dec > maxdec) {
        maxdec = pInput->Dec;
      }
      if (pInput->Dec < mindec) {
        mindec = pInput->Dec;
      }
         
      if (pInput->RA > maxra) {
        maxra = pInput->RA;
      }
      if (pInput->RA < minra) {
        minra = pInput->RA;
      }
    }
#ifndef NO_GETDECBIN
    dec_bin = GetDecBin(pGscBin,pInput->Dec,"formatcelestron2");
#endif /* NO_GETDECBIN */
    if (passCounter != 0) {
      if (max_read_dec_bin < dec_bin) {
        max_read_dec_bin = dec_bin;
      }
    }
    decBinCounts[dec_bin]++;
    decBinCountsInitial[dec_bin]++;
    skip_bin_size = dec_bin - old_dec_bin;
    if ((skip_bin_size > MAX_SKIP_BIN) ||
        (skip_bin_size < -MAX_SKIP_BIN)) {
      printf("Skip of %5d dec bins %5d (%6.2f) -> %5d (%6.2f) at line %10d\n",skip_bin_size,old_dec_bin,oldDec,dec_bin,pInput->Dec,numLines);
    }

#ifndef NO_GETDECBIN
    skip_bin_index = skip_bin_size + pGscBin->dec_bins;
    old_dec_bin = dec_bin;
    oldDec = pInput->Dec;
    if (skip_bin_index < 0) {
      skip_bin_index = 0;
    }
    if (skip_bin_index > (2*pGscBin->dec_bins)+1) {
      skip_bin_index > (2*pGscBin->dec_bins)+1;
    }
#endif /* NO_GETDECBIN */
    skipBinCounts[skip_bin_index]++;
    skipBinCountsTotal[skip_bin_index]++;


    if ((input_nrecs+numFileLines) == 0) {
      max_objid = pInput->objid;
      min_objid = pInput->objid;
    } else {
      if (pInput->objid > max_objid) {
        max_objid = pInput->objid;
      }
      if (pInput->objid < min_objid) {
        min_objid = pInput->objid;
      }

    }

    if (pInput->pmdec > maxpmdec) {
      maxpmdec = pInput->pmdec;
    } 
    if (pInput->pmdec < minpmdec) {
      minpmdec = pInput->pmdec;
    }


  
    temp_objid = pInput->objid;
    modulus = 1;
    for (modindex = 0; modindex < MAX_MODULUS; modindex++) {
      objidbits[modindex] |= temp_objid/modulus;
      modobjidbits[modindex] |= temp_objid % modulus;
      modulus = modulus * MODULUS_BASE;
    }


        

    rawSkipFlag = 0;
    if (raw_handle != NULL) {
      if ((pInput->dg > 0) && 
          (pInput->dr > 0) &&
          (pInput->dg < rmsCutoff) &&
          (pInput->dr < rmsCutoff)) {
        newrms = sqr(pInput->dg) + sqr(pInput->dr);
        if (newrms < maxrms2) {
          newrms = sqrt(newrms);
        } else {
          newrms = 99.0;
          rawSkipCount++;
          rawSkipFlag = 1;
        }
      } else {
        newrms = 99.0;
        rawSkipCount++;
        rawSkipFlag = 1;
      }
      fprintf(raw_handle,"%f\t%f\t%.2f\t%2f\n",pInput->RA,pInput->Dec,pInput->g,newrms);
  
    }

    tempmag = pInput->g-0.001;
    if (tempmag >= MAGBINS) {
      tempmag = MAGBINS+1;
    }
    if (tempmag < 0) {
      tempmag = MAGBINS+2;
    }

    for (index1 = 0; index1 < CELESTRON_MAXMAG; index1++) {
      skipFlag = 0;
      skipFlagBase = 0;
      switch(index1) {
      case CELESTRON_BP:
        value1 = pInput->BP;
        rms1   = pInput->dBP;
        break;
      case CELESTRON_RP:
        value1 = pInput->RP;
        rms1   = pInput->dRP;
        break;
      case CELESTRON_G:
        value1 = pInput->g;
        rms1   = pInput->dg;
        break;
      case CELESTRON_R:
        value1 = pInput->r;
        rms1   = pInput->dr;
        break;
      case CELESTRON_I:
        value1 = pInput->i;
        rms1   = pInput->di;
        break;
      case CELESTRON_Z:
        value1 = pInput->z;
        rms1   = pInput->dz;
        break;
      default:
        printf("ERROR: line %4d\n",__LINE__);
        exit(-1);
      }
      if ((rms1 > rmsCutoff) || (rms1 == 0)) {
        skipFlagBase = 1;
      }
      if (skipFlagBase == 0) {
        magcount[index1][tempmag]++;
        celestronCount[index1]++;
      }
      for (index2 = 0; index2 < CELESTRON_MAXMAG; index2++) {
        skipFlag = skipFlagBase;
        switch(index2) {
        case CELESTRON_BP:
          value2 = pInput->BP;
          rms2   = pInput->dBP;
          break;
        case CELESTRON_RP:
          value2 = pInput->RP;
          rms2   = pInput->dRP;
          break;
        case CELESTRON_G:
          value2 = pInput->g;
          rms2   = pInput->dg;
          break;
        case CELESTRON_R:
          value2 = pInput->r;
          rms2   = pInput->dr;
          break;
        case CELESTRON_I:
          value2 = pInput->i;
          rms2   = pInput->di;
          break;
        case CELESTRON_Z:
          value2 = pInput->z;
          rms2   = pInput->dz;
          break;
        default:
          printf("ERROR: line %4d\n",__LINE__);
          exit(-2);
        }
        if ((rms2 > rmsCutoff) && (rms2 == 0)) {
          skipFlag = 1;
        }
        if ((index1 == CELESTRON_G) &&
            (index2 == CELESTRON_R)) { 
          newrms = sqrt(sqr(rms1) + sqr(rms2));

#if 0
          printf("line %4d rms1 %f rms2 %f sum of squares %f sqrt %f\n",__LINE__,rms1,rms2,sqr(rms1) + sqr(rms2),newrms);
#endif
          intrms = (newrms*RMSBINS)/CELESTRON_MAXRMS;
          if (intrms >= RMSBINS) {
            intrms = RMSBINS+1;
          }
          if (intrms < 0) {
            intrms = RMSBINS+2;
          }
          rmscount[intrms]++;
          rmscountTotal[intrms][tempmag]++;


        }
        if ((index1 == CELESTRON_G) &&
            (index2 == CELESTRON_R))  {
          if (curOutputCount < maxAllocOutput) {
          } else {
            if (passCounter != 0) {
              printf("line %4d entering FlushOutputTable\n",__LINE__);
#ifndef NO_GETDECBIN
              FlushOutputTable(output_handle,pGscCelestronTable,maxAllocOutput,max_read_dec_bin,&curOutputCount,&outputWrittenCount,&output_dec_bin,decBinCountsOutput,0);
#endif /* NO_GETDECBIN */
            } else {
              printf("ERROR: unable to flush the output table in pass %d curOutputCount %d maxAllocOutput %d\n",passCounter,curOutputCount,maxAllocOutput);
              exit(-1);
            }
          }

          pGscCelestronImage = &pGscCelestronTable[curOutputCount];
          curOutputCount++;
          memset(pGscCelestronImage,0,sizeof(GSCCELESTRONIMAGE));
          pGscCelestronImage->celestronNumber = input_nrecs+numLines;
          pGscCelestronImage->REFNumber = pGscCelestronImage->celestronNumber + 9000000000L;
          pGscCelestronImage->objid = pInput->objid;
          pGscCelestronImage->ra = pInput->RA;
          pGscCelestronImage->dec = pInput->Dec;
          pGscCelestronImage->Stdmag =  pInput->g;
          pGscCelestronImage->color = pInput->g - pInput->r;
          pGscCelestronImage->rmag = pInput->r;
          pGscCelestronImage->imag = pInput->i;
          pGscCelestronImage->zmag = pInput->z;
          pGscCelestronImage->gmagerr = pInput->dg;
          pGscCelestronImage->rmagerr = pInput->dr;
          pGscCelestronImage->imagerr = pInput->di;
          pGscCelestronImage->zmagerr = pInput->dz;
          pGscCelestronImage->RaPM = pInput->pmra;
          pGscCelestronImage->DecPM = pInput->pmdec;
          pGscCelestronImage->RaSigmaPM = pInput->dpmra;
          pGscCelestronImage->DecSigmaPM = pInput->dpmdec;
          
          /* Shift the positions to epoch J2000 (see formatgaiadr2.) */
          dec2 = pGscCelestronImage->dec + (pGscCelestronImage->DecPM * (GSC_EQUINOX-GAIADR2_EPOCH))/(3600.0*1000.0);
          factor = cos(DEGREES_TO_RAD*(dec2));
          if (factor == 0) {
            ra2 = pGscCelestronImage->ra;
          } else {
            ra2  = pGscCelestronImage->ra + ((pGscCelestronImage->RaPM * (GSC_EQUINOX-GAIADR2_EPOCH))/(3600.0*1000.0*factor));
          }
          while (ra2 < 0) {
            ra2 += 360.0;
          } 
          while (ra2 >= 360.0) {
            ra2 -= 360.0;
          }
          flipFlag = 0;
          if (dec2 == 90.0) {
            dec2 -= 0.0001;
          } else if (dec2 > 90.0) {
            dec2 = 180.0-dec2;
            flipFlag = 1;
          }
          if (dec2 == -90.0) {
            dec2 += 0.0001;
          } else if (dec2 < -90) {
            dec2 = -180-dec2;
            flipFlag = 1;
          }
          if (flipFlag) {
            ra2 += 180.0;
            while (ra2 >= 360.0) {
              ra2 -= 360.0;
            }
          }
#ifndef NO_GETDECBIN
          dec_bin2 = GetDecBin(pGscBin,dec2,"formatcelestron3");
#endif /* NO_GETDECBIN */
          decBinCountsFinal[dec_bin2]++;
          pGscCelestronImage->dec = dec2;
          pGscCelestronImage->ra = ra2;
        

          if ((pGscCelestronImage->color < MIN_CELESTRON_COLOR) ||
              (pGscCelestronImage->color > MAX_CELESTRON_COLOR)) {
            colorReject++;
            skipFlag = 1;
          }
          outputCount++;
          /* Covert the rms added in quadrature to an integer */
          newrms = (127.*(sqrt(sqr(rms1) + sqr(rms2))))/scalerms;
#if 0
          printf("line %4d rms1 %f rms2 %f sum of squares %f sqrt %f scalerms %f \n",__LINE__,rms1,rms2,sqrt(sqr(rms1) + sqr(rms2)),newrms,scalerms);
#endif
          if (newrms > 127) {
            newrms = 127;
        
          }
      
          intrms = newrms;

          if ((skipFlag) || ((pInput->dupvar & VARIABLE) != 0) || (sqrt(sqr(rms1) + sqr(rms2)) > rmsCutoff))  {
            pGscCelestronImage->VFlag = 1;
            skipCount++;
          }
#if 0
          if ((pGscCelestronImage->ra < 185.5) ||
              (pGscCelestronImage->ra > 186.0) ||
              (pGscCelestronImage->dec < 63.5) ||
              (pGscCelestronImage->dec > 64.5)) {
            continue;
          }
#endif
								
#if 0
          if ((rawSkipFlag == 1)  && (pGscCelestronImage->VFlag == 0)) {
            printf("X");
          }

#endif

        }
        if (skipFlag == 0) {
          celestronDual[index1][index2]++;
        }
      }

    }

    if (verbose && (((input_nrecs+numFileLines) % 10000000) == 0)) {
      time(&curTime);
      curTime -= startTime;
      printf("Reading record %d at mindec %f maxdec %f minra %f maxra %f %d seconds\n",(input_nrecs+numFileLines),mindec,maxdec,minra,maxra,curTime);
      mindec = pInput->Dec;
      maxdec = pInput->Dec;
      minra = pInput->RA;
      maxra = pInput->RA;
     
    }

        
  }

  printf("line %4d: Reading record %d at mindec %f maxdec %f minra %f maxra %f %d seconds\n",__LINE__,(input_nrecs+numFileLines),mindec,maxdec,minra,maxra,curTime);
  mindec = pInput->Dec;
  maxdec = pInput->Dec;
  minra = pInput->RA;
  maxra = pInput->RA;
#endif
}
void SaveOutputTableLine(int lineno,PFILEENTRY pFileEntry,char *output_name,FILE *output_handle) {
  int tabIndex;
  printf("line %4d lineno %4d SaveOutputTableLine called for line %3d file '%s'\n and output_name '%s'\n",__LINE__,lineno,pFileEntry->numLines,pFileEntry->fileName,output_name);
#if 1
  fprintf(output_handle,"line %4d lineno %4d SaveOutputTableLine called for line %3d file '%s'\n and output_name '%s'\n",__LINE__,lineno,pFileEntry->numLines,pFileEntry->fileName,output_name);
#endif
  for (tabIndex = 0; tabIndex < pFileEntry->tabCount; tabIndex++) {
    fprintf(output_handle,"%s",&tabTable[tabIndex][0]);
    if (tabIndex != pFileEntry->maxTabs-1) {
      fprintf(output_handle,"\t");
    }
  }
  fprintf(output_handle,"\n");
#if 0
  for (tabIndex = 0; tabIndex < pFileEntry->tabCount; tabIndex++) {
    printf("%s",&tabTable[tabIndex][0]);
    if (tabIndex != pFileEntry->maxTabs-1) {
      printf("\t");
    }
  }
  printf("\n");
#endif
  return;


} /* End of PrintOutputTable */
