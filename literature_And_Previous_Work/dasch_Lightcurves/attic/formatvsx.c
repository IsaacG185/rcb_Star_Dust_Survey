// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 * Read in the vsx and galaxy catalogs
 *
 *   VSX catalog comes from http://vizier.cfa.harvard.edu/viz-bin/Cat?B/vsx
 *   Variability types are in http://www.sai.msu.su/groups/cluster/gcvs/gcvs/iii/vartype.txt (/dasch/scanner/backup/2015_08_28/vartype.txt)
 *   UZC (Zwicki) catalog comes from http://tdc-www.harvard.edu/uzc/uzcJ2000.tab.gz  (Download the UZC as a 452-Kb gzipped Starbase tab-separated table file) 
 *   HYPERLEDA. I. Catalog of galaxies http://vizier.cfa.harvard.edu/viz-bin/Cat?VII/237
 *   http://ned.ipac.caltech.edu/samples/NEDmdb.html  http://ned.ipac.caltech.edu/level5/CATALOGS/pgc.html
 *   SGC catalog http://ned.ipac.caltech.edu/cgi-bin/ex_refcode?refcode=1985SGC...C...0000C Corwin, H. G. Jr., de Vaucouleurs, A., de Vaucouleurs, G. 
 *     "Southern Galaxy Catalogue. A Catalogue of 5481 Galaxies South of Declination -17 Degrees Found on 1.2m U.K. Schmidt IIIa-J Plates", 1985, Austin: U. Texas"
 *     ftp://cdsarc.u-strasbg.fr/cats/VII/116/
 *   UGC Uppsala General Catalog of Galaxies 1973, Acta Universitatis Upsalienis, Nova Regiae Societatis Upsaliensis, Series v: a Vol.	,Nilson, P.
 *   ftp://cdsarc.u-strasbg.fr/cats/VII/26D/
 *    ftp://cdsarc.u-strasbg.fr/cats/J/ApJS/125/409/
 *   NGC2000 catalog  VizieR Online Data Catalog: NGC 2000.0 (Sky Publishing, ed. Sinnott 1988) 04/1997
 *      http://adsabs.harvard.edu/abs/1997yCat.7118....0S
 *      http://vizier.cfa.harvard.edu/viz-bin/Cat?cat=VII%2F118&target=http&
 *
 *   OGLE collection of variable stars.  arXiv 1701.03105v  ftp://ftp.astrouw.edu.pl/ogle/ogle4/OCVS/blg/ecl/
 *   "The OGLE Collection of Variable Stars. Over 450 000 Eclipsing and Ellipsoidal Binary Systems Toward the Galactic Bulge" Soszyński, I. et al., Acta Astronomica (2016) 66, 405
 *   Note the OGLE Subtype (C - contact, NC - non-contact, CV - cataclysmic variable, ELL - ellipsoidal variable) is different from the VAX type and will be prefixed with OGLE:.

 * cc -ggdb -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  formatvsx.c pipelineutils.a -L /dasch/install/lib -lm   -ltable -lutil -lwcs  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread  -o formatvsx
 *
 *    Variable types are in  ftp://cdsarc.u-strasbg.fr/cats/B/gcvs/vartype.txt
 *
 *   
 *
 *
 *    formatvsx -o /dasch/Pipeline/catalogs/galaxy.db -b  /dasch/Pipeline/catalogs/galaxy.dat
 *    update_limiting -n -l /dasch/Pipeline/ingest/file_nearby_update_limiting.log -o /dasch/Pipeline/catalogs/nearbylimiting.dat
 *
      formatvsx -d -o /dasch/Pipeline/catalogs/ngconly.db -b  /dasch/Pipeline/catalogs/ngconly.dat
      rm ngconly.dat
      awk --field-separator='\t' 'NR>2{print "J2000;CIRCLE(" $2 "," $3 ",0.05) # text ={" $1 "}"}' < /dasch/Pipeline/catalogs/ngconly.db | grep "NGC" > /dasch/Pipeline/catalogs/ngc.reg
 *
 *
 * Dec 14, 2012 Edward J. Los - Initial version
 * Dec 21, 2013 Edward J. Los - Find candidate variables for folding
 * Jul 20, 2015 Edward J. Los - Handle new VSX catalog, which now has 334219 entries vs 214028 new total entries 1366534
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * Mar  1, 2017 Edward J. Los - Add the OGLE catalog of eclipsing binaries
 * Mar  4, 2018 Edward J. Los - Add "-d" for NGC catalog only
 * Dec 18, 2018 Edward J. Los - Handle new definition of the "V" flag for duplicates.
 */ 
  

#define _GNU_SOURCE
#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include <stddef.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <assert.h>
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"
#include "kdtree.h"
#include "daschunistd.h"
#define MAX_INPUT_NAME 512
#define MAX_BUFFER 512
#define MAX_INPUT_LINES (2030000+450600)
#define INPUT_FILES 1
#define MAX_FILTER_SIZE 100
#define GALAXY_DEFAULT_RADIUS 6.0 /* Default galaxy radius in arcsec */
#define VARIABLE_SEARCH_RADIUS (5.0/3600.0)
#define VARIABLE_TEST_RADIUS 0.25 /* The minimum distance between two OGLE eclipsing binaries is 0.517 arcsec */
/* #define DEBUG_TYPE 1 */
/* #define DEBUG_ZWICKI 1 */

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;
void *ptree = NULL;
extern char *releaseFieldText[RELEASE_FIELD_MAX+1];
/* #define LOS_DEBUG 1 */  


typedef struct _vsximage {
  char OID[7];   /* Internal identifier */
  char blank01;
  char Name[30]; /* Variable star identifier */
  char blank02; 
  char V;        /* Variability flag  0 = variable; 1 = suspected; 2 = not variable */
  char blank03;
  char RAdeg[9]; /* Right ascension (J2000) */
  char blank04;
  char DEdeg[9]; /* Declination (J2000) */
  char blank05;
  char Type[30]; /* Variability type, as in GCVS catalog */
  char blank06;
  char blank07;
  char l_max;    /* Limit flag  on max */
  char max[6];   /* ? Magnitude at maximum, or amplitude */
  char u_max;    /* Uncertainty flag on max */
  char n_max[3]; /* Passband on max magnitide */
  char blank08;
  char blank09;
  char blank10;
  char blank11;
  char f_min;    /* [(] '(' indicates an amplitude */
  char l_min;    /* Limit flag on min */
  char min[6];   /* ? Magnitude at minimum, or amplitude */
  char u_min;    /* Uncertainty flag on min */
  char n_min[3]; /* Passband on min magnitude */
  char blank12;
  char blank13;
  char blank14;
  char blank15;
  char Epoch[12]; /* ? Epoch of maximum or minimum (HJD) */
  char u_Epoch;   /* [:)] Uncertainty flag (:) on epoch */
  char blank16;
  char l_period;  /* [<>(] Limit flag on period */
  char Period[16]; /* Period of variable in days */
  char u_period[3];  /* [: ) * / N2 ] Uncertainty flag on Period */
} VSXIMAGE,*PVSXIMAGE;

#define MAX_PGCID 10
#define MAX_PGCANAME 22
typedef struct _pgcimage {
  char pgc[MAX_PGCID];     /* pcg identifier */
  char blank01;
  char j;           /* 'J' */
  char ra[8];       /* Right Ascension */
  char dec[7];      /* Declination */
  char blank02;
  char Otype[2];    /* Object type */
  char blank03;
  char MTYpe[5];    /* Morphological type */
  char logD25[5];   /* logD25 apparent diameter */
  char plusmin[3];  /* '+/-' */
  char e_logD25[4]; /* error of apparent diameter */
  char blank04;
  char blank05;
  char logR25[4];   /* log of axis ratio */
  char plusmin2[3]; /* '+/-' */
  char e_logR25[4]; /* error of axis ratio */
  char blank06;
  char blank07;
  char PA[4];       /* Position angle */
  char plusmin3[3]; /* '+/-' */
  char e_PA[4];     /* Error in position angle */
  char blank08;  
  char o_ANames[2]; /* Alternate names */
  char blank09;
  char ANames[263]; /* Alternate names, each 22 bytes long */
} PGCIMAGE,*PPGCIMAGE;



typedef struct _sgcimage {
  char RAh[2];        /* Right Ascension hours (B1950) */
  char blank01;
  char RAm[2];        /* Right Ascension minutes */
  char blank02; 
  char RAs[2];        /* Right Ascension seconds */
  char blank03;
  char DE;            /* Declination sign */
  char DEd[2];        /* Declination degrees */
  char blank04;
  char DEm[4];        /* Declination minutes */
  char blank05;
  char RA2000h[2];    /* Right Ascension (hours) (J2000) */
  char blank06;
  char RA2000m[2];    /* Right Ascension (minutes) (J2000) */
  char blank07;
  char RA2000s[2];    /* Right Ascension (seconds) (J2000) */
  char blank08;
  char DE2000;        /* Declination sign */
  char DE2000d[2];    /* Declination degrees */
  char blank09;
  char DE2000m[4];    /* Declination minutes */
  char blank10;
  char GLON[6];       /* Galactic longitude */
  char blank11;
  char GLAT[6];       /* Galactic latitude */
  char blank12;
  char RC2[13];       /* RC2 name  */
  char blank13;
  char EU[9];         /* ESO/Uppsala name */
  char n_EU;          /* 'r' denotes a remark */
  char blank14;
  char MType[7];      /* Morphological type */
  char blank15;
  char T[4];          /* Weighted mean numerical type */
  char blank16;
  char o_T;           /* Number of estimates of T */
  char blank17;   
  char w_T[4];        /* Total internal weight of T */
  char L[4];          /* Weighted mean luminosity class */
  char blank18;
  char o_L;           /* Number of estimates of L */
  char blank19;
  char w_L[3];        /* Total weight of L */
  char blank20;
  char logD[4];       /* log of outer diamer along the major axis */
  char blank21;
  char w_LogD[4];     /* Total weight of log(D); */
  char blank22;
  char logR[4];       /* Weighted mean log of the axis ratio */
  char blank23;
  char w_logR[4];     /* Weight of R */
  char blank24;
  char I1;            /* Number of independent estimates of R */
} SGCIMAGE,*PSGCIMAGE;

typedef struct _ugc1image {
  char head[3];        /* "UGC" */
  char blank01;
  char UGC[5];        /* Catalog number */
  char A;             /* "A" if the galaxy is from the addenda */
  char RAh[2];        /* Right Ascension hours (B1950) */
  char RAm[4];        /* Right Ascension minutes */
  char DE;            /* Declination sign */
  char DEd[2];        /* Declination degrees */
  char DEm[2];        /* Declination minutes */
  char MCG[13];       /* MCG designation */
  char POSS[4];       /* POSS field number */
  char MajAxis[6];    /* Major axis in arcmin */
  char MinAxis[5];    /* Minor axis in arcmin */
  char PA[3];         /* Position angle in degrees */
  char Hubble[7];     /* Hubble classification */
  char blank02;
  char mag[4];        /* Photographic magnitude */
  char RadVal[5];     /* Radial Velocity */
  char aR[6];         /* Major axis in red */
  char bR[5];         /* Minor axis in red */
  char I1;            /* Inclination */
} UGC1IMAGE,*PUGC1IMAGE;

typedef struct _ugc2image {
  char head[3];        /* "UGC" */
  char blank01;
  char UGC[5];        /* Catalog number */
  char m_UGC[2];         /* Multiplicity number */
  char blank03;
  char blank04;
  char blank05;
  char blank06;
  char blank07;  
  char RAh[2];        /* Right Ascension hours   (J2000) */
  char blank08;
  char RAm[2];        /* Right Ascension minutes (J2000) */
  char blank09;
  char RAs[5];        /* Right Ascension seconds (J2000) */
  char blank10;
  char e_RAs[4];      /* RMS uncertainty in RA */
  char blank11;
  char DE;            /* Declination sign */
  char DEd[2];        /* Declination degrees (J2000) */
  char blank13;
  char DEm[2];        /* Declination minutes (J2000) */
  char blank14;
  char DEs[4];        /* Declination in arcsec (J2000) */
  char blank15;
  char e_DEs[4];      /* rms uncertainty in declination */

  char blank16;  
  char RA1h[2];        /* Right Ascension hours   (B1950) */
  char blank17;
  char RA1m[2];        /* Right Ascension minutes (B1950) */
  char blank18;
  char RA1s[5];        /* Right Ascension seconds (B1950) */
  char blank19;
  char e_RA1s[4];      /* RMS uncertainty in RA */
  char blank12;
  char blank21;
  char blank20;
  char blank25;
  char blank26;
  char blank27;
  char blank28;
  char blank29;
  char DE1;            /* Declination sign */
  char DE1d[2];        /* Declination degrees (B1950) */
  char blank22;
  char DE1m[2];        /* Declination minutes (B1950) */
  char blank23;
  char DE1s[4];        /* Declination in arcsec (B1950) */
  char blank24;
  char e_DE1s[4];      /* rms uncertainty in declination */
} UGC2IMAGE,*PUGC2IMAGE;



typedef struct _ogleidentimage {
  char ogleid[19];    /* Star's ID */
  char blank01[2];   
  char subtype[3];    /* Subtype (C - contact, NC - non-contact, CV - cataclysmic variable, ELL - ellipsoidal variable) */
  char blank02;
  char RAh[2];        /* Right ascension, equinox J2000.0 (hours) */
  char blank03;
  char RAm[2];        /* Right ascension, equinox J2000.0 (minutes) */
  char blank04; 
  char RAs[5];        /* Right ascension, equinox J2000.0 (seconds) */
  char blank05;
  char DE;            /* Declination, equinox J2000.0 (sign) */
  char DEd[2];        /* Declination, equinox J2000.0 (degrees) */
  char blank06;
  char DEm[2];        /* Declination, equinox J2000.0 (arc minutes) */
  char blank07;
  char DEs[4];        /* Declination, equinox J2000.0 (arc seconds) */
  char blank08[2];
  char ogleIV[16];    /* OGLE-IV ID */
  char blank09;
  char ogleIII[15];   /* OGLE-III ID */
  char blank10;
  char ogleII[15];    /* OGLE-II ID */
  char other[2];     /* Other designation (from VSX) */
} OGLEIDENTIMAGE,*POGLEIDENTIMAGE;

typedef struct _ogledatimage {
  char ogleid[19];    /* Star's ID */
  char blank01[2];   
  char imaxmag[6];    /* I-band magnitude at the maximum light */
  char blank02;
  char vmaxmag[6];    /* V-band magnitude at the maximum light */
  char blank03;
  char period[12];    /* Orbital period (days) */
  char blank04[2]; 
  char epoch[9];      /* Epoch of of the primary eclipse (HJD-2450000) */
  char blank05[2];
  char pridepth[5];   /* Depth of the primary eclipse (mag) */
  char blank06;
  char secdepth[5];   /* Depth of the secondary eclipse (mag) */
  char blank07;
} OGLEDATIMAGE,*POGLEDATIMAGE;



/* Holding structure for the UGC1 data */
#define MAX_UGC_RECORDS 15000

typedef struct _outugc {
  int ugc;           /* UGC integer identifier */
  char suffix;       /* Identifier suffix */
  double catalogmag;
  double radius; /* Radius in arcsec */
} OUTUGC,*POUTUGC;

/* Holding structure for OGLE data */

#define MAX_OGLE_RECORDS 450600
    /* The order of the ident file should be first ECL records numbered 1-425193 next ELL records numbered 1-25405 */
#define OGLE_ELL_BASE 425192

typedef struct _outogle {
  char ogleid[20];    /* Star's ID */
  int ECLflag;        /* 1 if OGLE-BLG-ECL-; 0 if OGLE-BLG-ECL- */
  int idnumber;
  char subtype[4];   /* Subtype */
  double ra;
  double dec;
  double imaxmag;
  double vmaxmag;
  double period;
  double epoch;
  double pridepth;
  double secdepth;
  int matchFlag;      /* Successfully matched with ecl.dat or ell.dat */
} OUTOGLE,*POUTOGLE;




#define MAX_ZNAME 15
#define MAX_ONAME 15
#define MAX_COORD 10

typedef struct _zwicki {
  double Zmag;
  char ra2000[MAX_COORD];
  char dec2000[MAX_COORD];
  char Zname[MAX_ZNAME];
  char Oname[MAX_ONAME];
} ZWICKI,*PZWICKI;
/* NOTE: the hyperleda catalog has NGC designations */
#define NGC_MINIMUM_RADIUS 30 /* Because the catalog is accurate to 1 arcmin, use a minimum radius of 30 arcsec */
typedef struct _ngcimage {
  char Name[5];  /* NGC or IC designation (preceded by I) */
  char blank1;
  char Type[3];  /* Object Classification */
  char blank2;
  char RAh[2];        /* Right Ascension hours   (B2000) */
  char blank3;
  char RAm[4];        /* Right Ascension minutes (B2000) */
  char blank4;
  char blank5;
  char DE1;            /* Declination sign */
  char DEd[2];         /* Declination degrees (B1950) */
  char blank6;
  char DEm[2];         /* Declination minutes (B1950) */
  char blank7;
  char Source;         /* Source of entry */
  char blank8;
  char blank9;
  char Const[3];       /* Constellation */
  char l_size;         /* (<) Limit on size */
  char size[5];        /* ? Largest dimension in arcmin */
  char blank10;
  char blank11;
  char mag[4];         /* ? Integrated magnitude, visual or photographic */
  char n_mag;          /*  'p' if mag is photographic (blue) */
  char blank12;
  char Desc[50];       /* Description of the object */
} NGCIMAGE,*PNGCIMAGE;

#define SOURCE_VSX  1
#define SOURCE_UZC  2 /* Zwicki */
#define SOURCE_PGC  3 /* Hyperleda */
#define SOURCE_SGC  4 /* Southern Galaxy Catalogue */
#define SOURCE_UGC  5 /* Uppsala General Catalog of Galaxies */
#define SOURCE_NGC  6 /* NGC2000 Catalog */
#define SOURCE_OGLE 7 /* OGLE variable star catalog */
#define MAX_SOURCE  8 /* Number of sources */

typedef struct _outrec {
  struct _outrec *pOutRec; /* Pointer to primary copy */
  double ra;
  double dec;
  double catalogmag;
  double radius; /* Radius in arcsec */
  double distance; /* Distance from search source */
  
  char catalogname[MAX_GALAXY_NAME];
  char type[MAX_GALAXY_TYPE];
  char source;
  char duplicate;  /* Duplicate if 1 */
  char checked;    /* Checked if 1 */
} OUTREC,*POUTREC;


void trimblank(char *buffer) {
  int length = strlen(buffer);
  while ((length-1) > 0) {
    if (buffer[length-1] != ' ') {
      break;
    }
    buffer[length-1] = 0;
    length--;
  }

}
int formatRA(char *buffer,double *ra)
{
  int degrees;
  int minutes;
  double seconds;
  int nvals;
  nvals = sscanf(buffer,"%2d%2d%lf",&degrees,&minutes,&seconds);
  if (nvals == 3) {
    *ra = 15.0*((1.0*degrees)+((1.0*minutes)/60.0)+(seconds/3600.));
#if 0
    if (*ra > 25.0) {
      printf("buffer %s ra %f\n",buffer,*ra);
    }
#endif
    return(0);
  } else {
    return(-1);
  }



}
int formatdec(char *buffer,double *dec)
{
  int degrees;
  int minutes;
  double seconds;
  double sign = 1.0;
  int nvals;
  if (buffer[0] == '-') {
    nvals = sscanf(&buffer[1],"%2d%2d%lf",&degrees,&minutes,&seconds);
    sign = -1.0;
  } else {
    nvals = sscanf(buffer,"%3d%2d%lf",&degrees,&minutes,&seconds);
  }
  if (nvals == 3) {
    *dec = sign * ((1.0*degrees)+((1.0*minutes)/60.0)+(seconds/3600.));


    return(0);
  } else {
    return(-1);
  }

}
int OutRecCompare(const void *first, const void *second) 
{
  double distanceFirst = ((POUTREC)first)->distance;
  double distanceSecond = ((POUTREC)second)->distance;
  if (distanceFirst > distanceSecond) {
    return(1);
  } else if (distanceFirst < distanceSecond) {
    return(-1);
  } else {
    return(0);
  }
}

void CheckDuplicates(POUTREC pOutRec,POUTREC filter_table,int filterSize,int cursource,int *pRadiusCount,int *pCatalogCount,int *pDuplicateCount,int variableFlag)
{
  int index;
  POUTREC pChkRec;
  POUTREC pFilterRec;
  int source[MAX_SOURCE]; 
  int testradius;
  memset(source,0,sizeof(source));
  /* First we want to sort the output table */
  qsort(filter_table,filterSize,sizeof(OUTREC),OutRecCompare);
  for (index = 0; index < filterSize; index++) {
    pFilterRec = &filter_table[index];
#if 0
    if (cursource == SOURCE_SGC) {
      printf("%12.6f arcsec filter %3d %31s %31s is near %31s %31s\n",pFilterRec->distance,filterSize,pFilterRec->catalogname,pFilterRec->type,pOutRec->catalogname,pOutRec->type);
    }
#endif
    if (variableFlag) {
      if (pFilterRec->distance > VARIABLE_TEST_RADIUS) {
        (*pRadiusCount)++;
        continue;
      }

    } else {

      testradius = pOutRec->radius;
      if (testradius == 0) {
        testradius = GALAXY_DEFAULT_RADIUS;
      }
      if (pFilterRec->distance > testradius) {
        (*pRadiusCount)++;
        continue;
      }
      if (source[pFilterRec->source] != 0) {
        (*pCatalogCount)++;
        continue; /* Do not consider a catalog twice */
      }
    }
    /* This is a valid duplicate */
#if 0
    if (cursource != SOURCE_PGC) {
      printf("%12.6f arcsec filter %3d %31s %31s is near %31s %31s\n",pFilterRec->distance,filterSize,pFilterRec->catalogname,pFilterRec->type,pOutRec->catalogname,pOutRec->type);
    }
#endif
    if (variableFlag == 0) {
      if ((pOutRec->catalogmag > 90)  && (pFilterRec->catalogmag < 90)) {
        pOutRec->catalogmag = pFilterRec->catalogmag;
      }
    }
    pChkRec = pFilterRec->pOutRec;
    pChkRec->checked = 1;
    pChkRec->duplicate = 1;
    (*pDuplicateCount)++;
    source[pFilterRec->source]++;
   
  }



}
#ifdef DEBUG_ZWICKI
void CheckZwicki(POUTREC pOutRec,POUTREC filter_table,int filterSize,int cursource,int *pRadiusCount,int *pCatalogCount,int *pDuplicateCount)
{
  int index;
  POUTREC pChkRec;
  POUTREC pFilterRec;
  int source[MAX_SOURCE]; 
  memset(source,0,sizeof(source));
  /* First we want to sort the output table */
  qsort(filter_table,filterSize,sizeof(OUTREC),OutRecCompare);
  for (index = 0; index < filterSize; index++) {    
    pFilterRec = &filter_table[index];
    if (index == 0) {
      printf("Minimum distance %12.f arcsec\n",pFilterRec->distance);
    }
#if 0
    if (cursource == SOURCE_SGC) {
      printf("%12.6f arcsec filter %3d %31s %31s is near %31s %31s\n",pFilterRec->distance,filterSize,pFilterRec->catalogname,pFilterRec->type,pOutRec->catalogname,pOutRec->type);
    }
#endif
    if (source[pFilterRec->source] != 0) {
      (*pCatalogCount)++;
      continue; /* Do not consider a catalog twice */
    }
    /* This is a valid duplicate */
    printf("%12.6f arcsec filter %3d %31s %31s is near %31s %31s\n",pFilterRec->distance,filterSize,pFilterRec->catalogname,pFilterRec->type,pOutRec->catalogname,pOutRec->type);
#if 0
    if ((pOutRec->catalogmag > 90)  && (pFilterRec->catalogmag < 90)) {
      pOutRec->catalogmag = pFilterRec->catalogmag;
    }
    pChkRec = pFilterRec->pOutRec;
    pChkRec->checked = 1;
    pChkRec->duplicate = 1;
    (*pDuplicateCount)++;
#endif
    source[pFilterRec->source]++;
   
  }


}

#endif /* DEBUG_ZWICKI */

int main(int argc,char *argv[])
{
  char *argstr;
  char output_name[MAX_INPUT_NAME];
  char galaxy_name[MAX_INPUT_NAME];
	char fold_name[MAX_INPUT_NAME];
  FILE* output_handle = NULL;
  File galaxy_handle = NULL;
  FILE* raw_handle = NULL;
	FILE* fold_handle = NULL;
  int totalCount = 0;
  int finalCount = 0;
  int finalVariableCount = 0;
  int finalGalaxyCount = 0;
  int notCheckedCount = 0;
  int constantCount = 0;
  int zeromagCount = 0;
  int errorFlag = 0;
  time_t startTime;
  time_t curTime;
  int linecounter = 0;
  int index1;
  int index2;
  char *charPtr;
  char cmdchar;
  int verbose = 0;
  int acceptNegative = 0;
  int ngcOnlyFlag = 0;
  int doStats = 0;
  char *dotPtr = NULL;
  int compareResult;
  int result;
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
  FILE *input_handle;
  double newrms;
  double scalerms;
  int intrms;
  int degrees;
  int minutes;
  double realminutes;
  int intseconds;
  double seconds;
  int o_ANames;
  int max_o_ANames = 0;
  double logD25;
  double D25;
  int outputIndex;
  int ngcnumber;
  double RAold;
  double decold;


  int inputFile;


  VSXIMAGE vsximage;
  PVSXIMAGE pVsxImage = &vsximage;

  PGCIMAGE pgcimage;
  PPGCIMAGE pPgcImage = &pgcimage;

  SGCIMAGE sgcimage;
  PSGCIMAGE pSgcImage = &sgcimage;
  
  NGCIMAGE ngcimage;
  PNGCIMAGE pNgcImage = &ngcimage;


  UGC1IMAGE ugc1image;
  PUGC1IMAGE pUgc1Image = &ugc1image;

  UGC2IMAGE ugc2image;
  PUGC2IMAGE pUgc2Image = &ugc2image;
  int ugc;
  char ugcsuffix;
  char ugcsuffix2;

  OGLEIDENTIMAGE ogleidentimage;
  POGLEIDENTIMAGE pOgleIdentImage = &ogleidentimage;

  OGLEDATIMAGE ogledatimage;
  POGLEDATIMAGE pOgleDatImage = &ogledatimage;

  POUTREC output_table = NULL;
  POUTREC pOutRec;
  POUTREC pChkRec;
  OUTREC filter_table[MAX_FILTER_SIZE];
  POUTREC pFilterRec;
  
  PGALAXYREC galaxy_table = NULL;
  PGALAXYREC pGalaxyRec;
  int writeItems;

  POUTUGC ugc_table = NULL;
  POUTUGC pOutUgc;
  int ugcIndex;
  int ugcCount = 0;

  POUTOGLE ogle_table = NULL;
  POUTOGLE pOutOgle;
  int ogleIndex;
  int ogleCount = 0;
  int ogleMatch = 0;
  int ogleBadVMaxMagCount = 0;
  int ogleBadPeriodCount = 0;
  OUTOGLE outogletemp;
  POUTOGLE pOutOgleTemp = &outogletemp;


  GSCIMAGE gsc_record;
  PGSCIMAGE pGscImage = &gsc_record;

  int tempmag;


  int skipOutput = 0;
  int magBinNumber;
  int colorType;
  int apassNoCQCount = 0;
  int apassHeaderCount = 0;
  int maxcqLength = 0;
  char maxcq[MAX_REF];
  double deltacolor;
  double maxdeltacolor = 0.0;
  int maxapasscolorid;
  char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
  int lineLen;
  char *inBuffer;
  int numLines;
  int nvals;
  double *vector = NULL;
  double rmsCutoff = 9.0;
  double value1;
  double rms1;
  double value2;
  double rms2;
  long long REFNumber;
  int skipFlag; /* Skip this flag because of neg rms or too high rms */
  int skipFlagBase;


  char *inputnames[INPUT_FILES] = {
    "/dasch/Pipeline/catalogs/vsx/vsx.dat"};

  File zwicki_handle = NULL;
  char zwicki_name[] = "/dasch/Pipeline/catalogs/uzc/uzcJ2000.tab";

  char hyperleda_name[] = "/dasch/Pipeline/catalogs/hyperleda/pgc.dat";

  char sgc_name[] = "/dasch/Pipeline/catalogs/sgc/catalog.dat";

  char ugc1_name[] = "/dasch/Pipeline/catalogs/ugc/catalog.dat";
  char ugc2_name[] = "/dasch/Pipeline/catalogs/ugc/table1.dat";

  char ngc_name[] = "/dasch/Pipeline/catalogs/ngc2000/ngc2000.dat";
  char ogleident_name[] = "/dasch/Pipeline/catalogs/ogle/ident.dat";
  char ogleecl_name[] = "/dasch/Pipeline/catalogs/ogle/ecl.dat";
  char ogleell_name[] = "/dasch/Pipeline/catalogs/ogle/ell.dat";

  TableHead zwicki_header = NULL;
  size_t zwicki_nrecs;
  size_t zwicki_index;
  PZWICKI zwicki_table;
  PZWICKI pZwicki;

  double GLAT;
  double GLON;
  double x;
  double y;
  double z;
  double cosx;
  double pt[3];
  double pos[3];
  struct kdres *presults;
  double distance;
  double minDistance;
  int maxFilterSize = 0;
  int resultSize;
  int maxResultSize = 0;
  int rejectCount = 0;
  int filterSize = 0;
  int radiusCount = 0;
  int catalogCount = 0;
  int duplicateCount = 0;
  int duplicate2Count = 0;
  int iteration;
  int cursource;
  int maxCatalogName = 0;
  int maxType = 0;
  double radius;


  assert(MAX_GALAXY_NAME > MAX_ONAME);
  assert(MAX_BUFFER > sizeof(VSXIMAGE));
  assert(MAX_GALAXY_NAME > MAX_PGCID);
  assert(MAX_BUFFER > sizeof(SGCIMAGE));
  assert(MAX_BUFFER > sizeof(PGCIMAGE));
  assert(MAX_BUFFER > sizeof(ZWICKI));
  assert(MAX_BUFFER > sizeof(UGC1IMAGE));
  assert(MAX_BUFFER > sizeof(UGC2IMAGE));
  assert(MAX_BUFFER > sizeof(NGCIMAGE));
  assert(MAX_BUFFER > sizeof(OGLEIDENTIMAGE));
  assert(MAX_BUFFER > sizeof(OGLEDATIMAGE));

#ifdef DEBUG_TYPE
  printf("ERROR: DEBUG_TYPE is set\n");
#endif

#if 0
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,OID),"OID");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,Name),"Name");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,V),"V");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,RAdeg),"RAdeg");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,DEdeg),"DEdeg");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,blank05),"blank05");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,Type),"Type");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,l_max),"l_max");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,max),"max");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,u_max),"u_max");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,n_max),"n_max");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,f_min),"f_min");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,l_min),"l_min");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,mag),"mag");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,u_min),"u_min");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,n_min),"n_min");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,Epoch),"Epoch");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,u_Epoch),"u_Epoch");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,l_period),"l_period");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,Period),"Period");
  printf("Offset %3d for %s\n",1+offsetof(VSXIMAGE,u_period),"u_period");
  printf("Size of VSXIMAGE %d\n",sizeof(VSXIMAGE));
#endif
#if 0

  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,pgc),"pgc");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank01),"blank01");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,j),"j");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,ra),"ra");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,dec),"dec");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank02),"blank02");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,Otype),"Otype");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank03),"blank03");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,MTYpe),"MTYpe");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,logD25),"logD25");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,plusmin),"plusmin");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,e_logD25),"e_logD25");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank04),"blank04");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank05),"blank05");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,logR25),"logR25");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,plusmin2),"plusmin2");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,e_logR25),"e_logR25");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank06),"blank06");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank07),"blank07");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,PA),"PA");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,plusmin3),"plusmin3");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,e_PA),"e_PA");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank08),"blank08");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,o_ANames),"o_ANames");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,blank09),"blank09");
  printf("Offset %3d for %s\n",1+offsetof(PGCIMAGE,ANames),"ANames");
  printf("Size of PGCIMAGE %d\n",sizeof(PGCIMAGE));

#endif
#if 0
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,RAh),"RAh");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank01),"blank01");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,RAm),"RAm");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank02),"blank02");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,RAs),"RAs");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank03),"blank03");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,DE),"DE");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,DEd),"DEd");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank04),"blank04");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,DEm),"DEm");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank05),"blank05");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,RA2000h),"RA2000h");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank06),"blank06");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,RA2000m),"RA2000m");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank07),"blank07");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,RA2000s),"RA2000s");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,DE2000),"DE2000");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,DE2000d),"DE2000d");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank08),"blank08");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,DE2000m),"DE2000m");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank09),"blank09");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,GLON),"GLON");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank10),"blank10");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,GLAT),"GLAT");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank11),"blank11");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,RC2),"RC2");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank12),"blank12");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,EU),"EU");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,n_EU),"n_EU");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank13),"blank13");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,MType),"MType");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank14),"blank14");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,T),"T");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank15),"blank15");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,o_T),"o_T");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank16),"blank16");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,w_T),"w_T");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,L),"L");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank17),"blank17");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,o_L),"o_L");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank18),"blank18");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,w_L),"w_L");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank19),"blank19");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,logD),"logD");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank20),"blank20");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,w_LogD),"w_LogD");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank21),"blank21");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,logR),"logR");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank22),"blank22");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,w_logR),"w_logR");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,blank23),"blank23");
  printf("Offset %3d for %s\n",1+offsetof(SGCIMAGE,I1),"I1");
  printf("Size of SGCIMAGE %d\n",sizeof(SGCIMAGE));
#endif

#if 0
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,head),"head");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,blank01),"blank01");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,UGC),"UGC");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,A),"A");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,RAh),"RAh");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,RAm),"RAm");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,DE),"DE");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,DEd),"DEd");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,DEm),"DEm");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,MCG),"MCG");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,POSS),"POSS");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,MajAxis),"MajAxis");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,MinAxis),"MinAxis");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,PA),"PA");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,Hubble),"Hubble");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,mag),"mag");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,RadVal),"RadVal");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,aR),"aR");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,bR),"bR");
  printf("Offset %3d for %s\n",1+offsetof(UGC1IMAGE,I1),"I1");
  printf("Size of UGC1IMAGE %d\n\n",sizeof(UGC1IMAGE));

  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,head),"head");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank01),"blank01");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,UGC),"UGC");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,m_UGC),"m_UGC");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank03),"blank03");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank04),"blank04");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank05),"blank05");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank06),"blank06");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank07),"blank07");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,RAh),"RAh");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank08),"blank08");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,RAm),"RAm");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank09),"blank09");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,RAs),"RAs");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank10),"blank10");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,e_RAs),"e_RAs");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank11),"blank11");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,DE),"DE");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,DEd),"DEd");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank13),"blank13");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,DEm),"DEm");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank14),"blank14");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,DEs),"DEs");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank15),"blank15");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,e_DEs),"e_DEs");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank16),"blank16");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,RA1h),"RA1h");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank17),"blank17");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,RA1m),"RA1m");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank18),"blank18");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,RA1s),"RA1s");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank19),"blank19");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,e_RA1s),"e_RA1s");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank20),"blank20");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,DE1),"DE1");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank21),"blank21");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,DE1d),"DE1d");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank22),"blank22");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,DE1m),"DE1m");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank23),"blank23");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,DE1s),"DE1s");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,blank24),"blank24");
  printf("Offset %3d for %s\n",1+offsetof(UGC2IMAGE,e_DE1s),"e_DE1s");
  printf("Size of UGC2IMAGE %d\n",sizeof(UGC2IMAGE));
#endif

#if 0
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,Name),"Name");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank1),"blank1");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,Type),"Type");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,RAh),"RAh");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank2),"blank2");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,RAm),"RAm");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank3),"blank3");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank4),"blank4");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,DE1),"DE1");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,DEd),"DEd");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank5),"blank5");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,DEm),"DEm");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank6),"blank6");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,Source),"Source");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank8),"blank8");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank9),"blank9");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,Const),"Const");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,l_size),"l_size");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,size),"size");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank10),"blank10");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank11),"blank11");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,mag),"mag");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,n_mag),"n_mag");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,blank12),"blank12");
  printf("Offset %3d for %s\n",1+offsetof(NGCIMAGE,Desc),"Desc");
  printf("sizeof(NGCIMAGE) %d\n",sizeof(NGCIMAGE));

#endif

#if 0
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,ogleid),"ogleid");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank01),"blank01");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,subtype),"subtype");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank02),"blank02");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,RAh),"RAh");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank03),"blank03");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,RAm),"RAm");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank04),"blank04");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,RAs),"RAs");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank05),"blank05");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,DE),"DE");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,DEd),"DEd");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank06),"blank06");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,DEm),"DEm");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank07),"blank07");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,DEs),"DEs");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank08),"blank08");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,ogleIV),"ogleIV");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank09),"blank09");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,ogleIII),"ogleIII");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,blank10),"blank10");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,ogleII),"ogleII");
  printf("Offset %3d for %s\n",1+offsetof(OGLEIDENTIMAGE,other),"other");
#endif

#if 0
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,ogleid),"ogleid");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,blank01),"blank01");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,imaxmag),"imaxmag");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,blank02),"blank02");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,vmaxmag),"vmaxmag");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,blank03),"blank03");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,period),"period");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,blank04),"blank04");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,epoch),"epoch");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,blank05),"blank05");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,pridepth),"pridepth");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,blank06),"blank06");
  printf("Offset %3d for %s\n",1+offsetof(OGLEDATIMAGE,secdepth),"secdepth");


  exit(-1);
#endif 

  /* Loop through the arguments */
  output_name[0] = 0;
  galaxy_name[0] = 0;
	fold_name[0] = 0;

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

        case 'n':
        case 'N':
          acceptNegative = 1;
          break;

        case 'd':
        case 'D':
          ngcOnlyFlag = 1;
          break;

        case 's':
        case 'S':
          doStats = 1;
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

        case 'b': /* binary file name */
        case 'B':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(galaxy_name,*++argv,MAX_INPUT_NAME-2);
            if (strlen(*argv) >= MAX_INPUT_NAME-2) {
              fprintf(stderr,"ERROR: MAX_INPUT_NAME too small for %s\n",*argv);
            }
          }
          break;

        case 'f': /* fold file name */
        case 'F':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(fold_name,*++argv,MAX_INPUT_NAME-2);
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

	if (fold_name[0] != 0) {
    fold_handle = fopen(fold_name,"wt");
    if (fold_handle == NULL) {
      errorFlag = 1;
      printf("Could not open fold file %s\n",fold_name);
    }
		fprintf(fold_handle,"release field  amplitude    period     min or amplitude   maximum    ra        dec                             name     GCVS type\n");
		fprintf(fold_handle,"                 (mag)      (days)          (mag)          (mag)    (deg)     (deg)                                              \n");

	}

  scalerms = rmsCutoff * sqrt(2.0);


  zwicki_handle = Open(zwicki_name,"r");
  if (zwicki_handle == NULL) {
    printf("ERROR: Failed to find Zwicki filefile %s\n",zwicki_name);
    errorFlag = 1;
  } else {
    if (verbose) {
      printf("Found the Zwicki file file %s\n",zwicki_name);
    }
  }


  if (errorFlag) {
    printf("Usage: formatapass -o {output name} [-r <raw file>] [-v][-n][-s][-c {cutoff}]\n");
		printf("       where -b <filename> is the binary file name \n");
		printf("       where -r <filename> is a raw table of RA and DEC for plotting \n");
    printf("       where -v is the verbose flag\n");
    printf("       where -s computes general statistics\n");
    printf("       where -c {cutoff} gives the RMS cutoff values\n");
    printf("       where -n accepts negative RMS values indicating only one measurement\n");
		printf("       where -f <foldfile> provides a table of fold candidates\n");
    printf("       where -d formats the NGC catalog only\n");
    return(-1);
  }

  printf("formatapass of %s %s Output Filename %s Binary Filename %s Fold filename %d acceptNegative %d,rmsCutoff %f\n",
         __DATE__,__TIME__,output_name,galaxy_name,fold_name,acceptNegative,rmsCutoff);
  printf("Size of BININDEX is %d.  Size of STARINDEX is %d. Size of VSXIMAGE is %d Size of GALAXYREC is %d \n",sizeof(BININDEX),sizeof(STARINDEX),sizeof(VSXIMAGE),sizeof(GALAXYREC));
 

  time(&startTime);
  
#ifdef DEBUG_TYPE
  fprintf(output_handle,"catalogname\tra\tdec\tsequence\tcatalogmag\tradius\n");
#else /* DEBUG_TYPE */
  fprintf(output_handle,"catalogname\tra\tdec\ttype\tcatalogmag\tradius\n");
#endif /* DEBUG_TYPE */
  fprintf(output_handle,"-----------\t--\t---\t----\t----------\t------\n");


  output_table = (POUTREC)calloc(MAX_INPUT_LINES,sizeof(OUTREC));
  if (output_table == NULL) {
    printf("ERROR: failed to allocate input table of size %d\n",(MAX_INPUT_LINES*sizeof(OUTREC)));
    exit(-1);
  }

  ugc_table = (POUTUGC)calloc(MAX_UGC_RECORDS,sizeof(OUTUGC));
  if (ugc_table == NULL) {
    printf("ERROR: failed to allocate ugc table of size %d\n",(MAX_UGC_RECORDS*sizeof(OUTUGC)));
    exit(-1);
  }

  ogle_table = (POUTOGLE)calloc(MAX_OGLE_RECORDS,sizeof(OUTOGLE));
  if (ogle_table == NULL) {
    printf("ERROR: failed to allocate ogle table of size %d\n",(MAX_OGLE_RECORDS*sizeof(OUTOGLE)));
    exit(-1);
  }

  vector = (double *)calloc(MAX_INPUT_LINES,sizeof(double));
  if (vector == NULL) {
    printf("ERROR: failed to allocate vector of size %d\n",(MAX_INPUT_LINES*sizeof(VSXIMAGE)));
    exit(-1);
  }
  if (ngcOnlyFlag == 0) {
    for (inputFile = 0; inputFile < INPUT_FILES; inputFile++) {
      printf("Opening %s\n",inputnames[inputFile]);
      input_handle = fopen(inputnames[inputFile],"rt");
      if (input_handle == NULL) {
        printf("ERROR: failed to open %s\n",inputnames[inputFile]);
        exit(-1);
      }
      numLines = 0;
      while(1) {
        memset(pVsxImage,0,sizeof(VSXIMAGE));
        inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
        if (inBuffer == NULL) {
          break;
        }
        strcpy(copyLine,inBuffer);
        numLines++;
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
        if (lineLen < offsetof(VSXIMAGE,blank05)) {
          printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
          continue;
        }
        /* Pad the line with blanks */
        while (lineLen < sizeof(VSXIMAGE)) {
          inBuffer[lineLen] = ' ';
          inBuffer[lineLen+1] = 0;
          lineLen++;
        }
        memcpy(pVsxImage,inBuffer,sizeof(VSXIMAGE));
        if ((pVsxImage->blank02 != ' ') ||
            (pVsxImage->blank04 != ' ') ||
            (pVsxImage->blank05 != ' ') ||
            (pVsxImage->blank06 != ' ')) {
          printf("ERROR: line '%s' has an invalid record type\n",copyLine);
          continue;
        }
        pVsxImage->blank02 = 0;
        pVsxImage->blank05 = 0;
        pVsxImage->blank06 = 0;
        pOutRec = &output_table[totalCount];
        memset(pOutRec,0,sizeof(OUTREC));
        pOutRec->source = SOURCE_VSX;
        nvals = sscanf((char *)&pVsxImage->RAdeg,"%lf %lf",&pOutRec->ra,&pOutRec->dec);
        if (nvals != 2) {
          printf("ERROR: line '%s' has undecodable ra and dec\n",copyLine);
          continue;
        }
        switch (pVsxImage->V) {
        case '0':
        case '1':
          /* A variable or a suspected variable */
          break;
        case '2':
          constantCount++;
          continue;
        case '3':
          duplicate2Count++;
          continue;
        default:
          printf("ERROR: line '%s' has an illegal variability flag\n",copyLine);
          continue;
        }
        pVsxImage->u_max = 0;
        nvals = sscanf((char *)&pVsxImage->max,"%lf",&pOutRec->catalogmag);
        if (nvals != 1) {
          pOutRec->catalogmag = 99;
        }
#if 0
        if (pOutRec->catalogmag < 0) {
          printf("Bright: %s\n",copyLine);
        }
#endif
        strcpy(pOutRec->catalogname,pVsxImage->Name);
        trimblank(pOutRec->catalogname);
        strcpy(pOutRec->type,pVsxImage->Type);
        trimblank(pOutRec->type);
        if ((fold_handle != NULL) &&
            (pOutRec->catalogmag < 90.0)) {
          /* Here we are looking for fold candidates.
             1.  The period must be defined.
             2.  The magnitude must be less than 13.0
             3.  The amplitude must be greater than 0.5
          */
          int nvalsPeriod;
          int nvalsMinMagnitude;
          double period;
          double minMagnitude;
          int releaseField;
          double amplitude;
          pVsxImage->u_period[0] = 0;
          pVsxImage->u_min = 0;
          nvalsPeriod = sscanf((char *)&pVsxImage->Period,"%lf",&period);
          nvalsMinMagnitude = sscanf((char *)&pVsxImage->min,"%lf",&minMagnitude);

          if ((nvalsPeriod == 1) &&
              (nvalsMinMagnitude == 1) && 
              (pOutRec->catalogmag < 13.0)) {
            amplitude = 0.0;
            if ((minMagnitude > 0.5) && (minMagnitude < 5.0)) {
              amplitude = minMagnitude;
            } else if ((minMagnitude-pOutRec->catalogmag) > 0.5) {
              amplitude = minMagnitude-pOutRec->catalogmag;
            }
            if (amplitude > 0.5) {
              CheckAuthorization("guest",pOutRec->ra,pOutRec->dec,&releaseField);
              fprintf(fold_handle," %6s          %4.2f   %16s    %5.2f         %5.2f  %8.4f  %8.4f  %30s      %s \n",releaseFieldText[releaseField],amplitude,pVsxImage->Period,minMagnitude,pOutRec->catalogmag,pOutRec->ra,pOutRec->dec,pOutRec->catalogname,pOutRec->type);
            }
          }
			

        }



        totalCount++;
        if ((totalCount) >= MAX_INPUT_LINES) {
          printf("ERROR: MAX_INPUT_LINES %d exceeded line %d\n",MAX_INPUT_LINES,__LINE__);
          exit(-1);
        }
      }
      printf("Read %d lines from %s\n",numLines,inputnames[inputFile]);
      fclose(input_handle);

    }
    if (fold_handle != NULL) {
      fclose(fold_handle);
      fold_handle = NULL;
    }
		

    /* Now read in the Zwicki catalog */
    zwicki_header = table_header(zwicki_handle,TABLE_PARSE);
    if (zwicki_header == NULL) {
      printf("ERROR: Failed to read header for %s\n",zwicki_name);
    
      exit(-1);
    }
    zwicki_table = table_loadva(zwicki_handle,
                                &zwicki_header,
                                NULL, /* hbase */
                                NULL, /* rows */
                                NULL,
                                sizeof(ZWICKI),
                                &zwicki_nrecs,
                                TblDbl,"Zmag"   ,TblOff(PZWICKI,Zmag),
                                TblBuf,"ra2000" ,TblOff(PZWICKI,ra2000),MAX_COORD,
                                TblBuf,"dec2000",TblOff(PZWICKI,dec2000),MAX_COORD,
                                TblBuf,"Zname"  ,TblOff(PZWICKI,Zname),MAX_ZNAME,
                                TblBuf,"Oname"  ,TblOff(PZWICKI,Oname),MAX_ONAME,
                                0,"end",0);
    if (zwicki_table == NULL) {
      printf("ERROR: Failed to read table for %s\n",zwicki_name);
      exit(-1);
    }
    printf("Read %d records from %s\n",zwicki_nrecs,zwicki_name);
    if ((totalCount + zwicki_nrecs) >= MAX_INPUT_LINES) {
      printf("ERROR: MAX_INPUT_LINES %d needs to be %d\n",MAX_INPUT_LINES,totalCount + zwicki_nrecs+1);
      exit(-1);
    }
    for (zwicki_index = 0; zwicki_index < zwicki_nrecs; zwicki_index++) {
      pZwicki = &zwicki_table[zwicki_index];

      pOutRec = &output_table[totalCount];
      memset(pOutRec,0,sizeof(OUTREC));
      pOutRec->source = SOURCE_UZC;
      if (formatRA(pZwicki->ra2000,&pOutRec->ra) != 0) {
        printf("ERROR: failed to decode record %d ra %s\n",zwicki_index,pZwicki->ra2000);
        continue;
      }
      if (formatdec(pZwicki->dec2000,&pOutRec->dec) != 0) {
        printf("ERROR: failed to decode record %d dec %s\n",zwicki_index,pZwicki->dec2000);
        continue;
      }
      pOutRec->catalogmag = pZwicki->Zmag;
      strcpy(pOutRec->catalogname,pZwicki->Oname);
      trimblank(pOutRec->catalogname);
      strcpy(pOutRec->type,"UZC Galaxy");

      totalCount++;
      if ((totalCount) >= MAX_INPUT_LINES) {
        printf("ERROR: MAX_INPUT_LINES %d exceeded line %d\n",MAX_INPUT_LINES,__LINE__);
        exit(-1);
      }

    }
    if (zwicki_table != NULL) {
      Free(zwicki_table);
      zwicki_table = NULL;
    }

    if (zwicki_header != NULL) {
      table_hdrfree(zwicki_header);
      zwicki_header = NULL;
    }

    if (zwicki_handle != NULL) {
      Close(zwicki_handle);
      zwicki_handle = NULL;
    }
    /* Now read in the HYPERLEDA catalog */
    printf("Opening %s\n",hyperleda_name);
    input_handle = fopen(hyperleda_name,"rt");
    if (input_handle == NULL) {
      printf("ERROR: failed to open %s\n",hyperleda_name);
      exit(-1);
    }
    numLines = 0;
    while(1) {
      memset(pPgcImage,0,sizeof(PGCIMAGE));
      inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
      if (inBuffer == NULL) {
        break;
      }
      strcpy(copyLine,inBuffer);
      numLines++;
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
      if (lineLen < offsetof(PGCIMAGE,blank09)) {
        printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
        continue;
      }
      /* Pad the line with blanks */
      while (lineLen < sizeof(PGCIMAGE)) {
        inBuffer[lineLen] = ' ';
        inBuffer[lineLen+1] = 0;
        lineLen++;
      }
      memcpy(pPgcImage,inBuffer,sizeof(PGCIMAGE));
      if ((pPgcImage->pgc[0] != 'P') ||
          (pPgcImage->pgc[1] != 'G') ||
          (pPgcImage->pgc[2] != 'C')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      if ((pPgcImage->blank01 != ' ') ||
          (pPgcImage->blank02 != ' ') ||
          (pPgcImage->blank03 != ' ') ||
          (pPgcImage->blank04 != ' ') ||
          (pPgcImage->blank05 != ' ') ||
          (pPgcImage->blank06 != ' ') ||
          (pPgcImage->blank07 != ' ') ||
          (pPgcImage->blank08 != ' ') ||
          (pPgcImage->blank09 != ' ')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      if ((pPgcImage->plusmin[0] != '+') ||
          (pPgcImage->plusmin2[0] != '+') ||
          (pPgcImage->plusmin3[0] != '+') ||
          (pPgcImage->plusmin[1] != '/') ||
          (pPgcImage->plusmin2[1] != '/') ||
          (pPgcImage->plusmin3[1] != '/') ||
          (pPgcImage->plusmin[2] != '-') ||
          (pPgcImage->plusmin2[2] != '-') ||
          (pPgcImage->plusmin3[2] != '-')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      pPgcImage->blank09 = 0;
      nvals = sscanf(pPgcImage->o_ANames,"%d",&o_ANames);
      if (nvals != 1) {
        printf("ERROR: line '%s' has undecodable o_ANames\n",copyLine);
        continue;
      }
      if (o_ANames > max_o_ANames) {
        max_o_ANames = o_ANames;
      }
      pPgcImage->blank01 = 0;
      pOutRec = &output_table[totalCount];
      memset(pOutRec,0,sizeof(OUTREC));
      pOutRec->source = SOURCE_PGC;
      if (o_ANames > 0) {
        pPgcImage->ANames[MAX_PGCANAME] = 0;
        strcpy(pOutRec->catalogname,pPgcImage->ANames);
      } else {
        strcpy(pOutRec->catalogname,pPgcImage->pgc);
      }
      trimblank(pOutRec->catalogname);
      strcpy(pOutRec->type,"PGC Galaxy");
      pPgcImage->blank02 = 0;
      if (formatdec(pPgcImage->dec,&pOutRec->dec) != 0) {
        printf("ERROR: Failed to decode declination in %s\n",copyLine);
        continue;
      }
      pPgcImage->dec[0] = 0;
      if (formatRA(pPgcImage->ra,&pOutRec->ra) != 0) {
        printf("ERROR: Failed to decode right ascension in %s\n",copyLine);
        continue;
      }
      pOutRec->catalogmag = 99.0;
      pPgcImage->plusmin[0] = 0;
      nvals = sscanf(pPgcImage->logD25,"%lf",&logD25);
      if (nvals != 1) {
        printf("ERROR: Failed to decode logD25 in %s\n",copyLine);
        continue;
      }
      if ((logD25 <= 9.989) || (logD25 >= 10.00)) {
        D25 = 0.1*exp10(logD25); /* Diameter in arcmin */
        pOutRec->radius = 60.0*D25/2.0; /* Radius in arcsec */
     
      }
      totalCount++;
      if ((totalCount) >= MAX_INPUT_LINES) {
        printf("ERROR: MAX_INPUT_LINES %d exceeded line %d\n",MAX_INPUT_LINES,__LINE__);
        exit(-1);
      }
    }
    printf("Read %d lines max o_ANames is %d from %s\n",numLines,max_o_ANames,hyperleda_name);
    fclose(input_handle);



    /* Now read in the SGC catalog */
    printf("Opening %s\n",sgc_name);
    input_handle = fopen(sgc_name,"rt");
    if (input_handle == NULL) {
      printf("ERROR: failed to open %s\n",sgc_name);
      exit(-1);
    }
    numLines = 0;
    while(1) {
      memset(pSgcImage,0,sizeof(SGCIMAGE));
      inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
      if (inBuffer == NULL) {
        break;
      }
      strcpy(copyLine,inBuffer);
      numLines++;
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
      if (lineLen < offsetof(SGCIMAGE,I1)) {
        printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
        continue;
      }
      /* Pad the line with blanks */
      while (lineLen < sizeof(SGCIMAGE)) {
        inBuffer[lineLen] = ' ';
        inBuffer[lineLen+1] = 0;
        lineLen++;
      }
      memcpy(pSgcImage,inBuffer,sizeof(SGCIMAGE));
      if ((pSgcImage->blank01 != ' ') ||
          (pSgcImage->blank02 != ' ') ||
          (pSgcImage->blank03 != ' ') ||
          (pSgcImage->blank04 != ' ') ||
          (pSgcImage->blank05 != ' ') ||
          (pSgcImage->blank06 != ' ') ||
          (pSgcImage->blank07 != ' ') ||
          (pSgcImage->blank08 != ' ') ||
          (pSgcImage->blank09 != ' ') ||
          (pSgcImage->blank10 != ' ') ||
          (pSgcImage->blank11 != ' ') ||
          (pSgcImage->blank12 != ' ') ||
          (pSgcImage->blank13 != ' ') ||
          (pSgcImage->blank14 != ' ') ||
          (pSgcImage->blank15 != ' ') ||
          (pSgcImage->blank16 != ' ') ||
          (pSgcImage->blank17 != ' ') ||
          (pSgcImage->blank18 != ' ') ||
          (pSgcImage->blank19 != ' ') ||
          (pSgcImage->blank20 != ' ') ||
          (pSgcImage->blank21 != ' ') ||
          (pSgcImage->blank22 != ' ') ||
          (pSgcImage->blank23 != ' ')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      pOutRec = &output_table[totalCount];
      memset(pOutRec,0,sizeof(OUTREC));
      pOutRec->source = SOURCE_SGC;
      strcpy(pOutRec->type,"SGC Galaxy");
      pOutRec->catalogmag = 99.;
      if (pSgcImage->RC2[0] == ' ') {
        pSgcImage->n_EU = 0;
        strcpy(pOutRec->catalogname,pSgcImage->EU);
      } else {
        pSgcImage->blank13 = 0;
        strcpy(pOutRec->catalogname,pSgcImage->RC2);
      }
      trimblank(pOutRec->catalogname);
#if 0

      if (strcmp(pOutRec->catalogname,"N5890") == 0) {
        printf("At %s\n",pOutRec->catalogname);
      }
#endif

      pSgcImage->blank08 = 0;
      pSgcImage->blank10 = 0;
      nvals = sscanf(pSgcImage->RA2000h,"%2d %2d %2d",&degrees,&minutes,&intseconds);
      if (nvals == 3) {
        pOutRec->ra = 15.0*((1.0*degrees)+((1.0*minutes)/60.0)+((1.0*intseconds)/3600.));
      } else {
        printf("ERROR: failed to decode record %d ra %s\n",numLines,pSgcImage->RA2000h);
        continue;
      }
      nvals = sscanf(pSgcImage->DE2000d,"%2d %lf",&degrees,&realminutes);
      if (nvals == 2) {
        pOutRec->dec = ((1.0*degrees)+(realminutes/60.0));
        if (pSgcImage->DE2000 == '-') {
          pOutRec->dec = -pOutRec->dec;
        }
      } else {
        printf("ERROR: failed to decode record %d dec %s\n",numLines,pSgcImage->DE2000);
        continue;
      }


      pSgcImage->blank22 = 0;
      nvals = sscanf(pSgcImage->logD,"%lf",&logD25);
      if (nvals == 1) {
        D25 = 0.1*exp10(logD25); /* Diameter in arcmin */
        pOutRec->radius = 60.0*D25/2.0; /* Radius in arcsec */
      }
      totalCount++;
      if ((totalCount) >= MAX_INPUT_LINES) {
        printf("ERROR: MAX_INPUT_LINES %d exceeded line %d\n",MAX_INPUT_LINES,__LINE__);
        exit(-1);
      }
    }
    printf("Read %d lines from %s\n",numLines,sgc_name);
    fclose(input_handle);
  } /* ngcOnlyFlag */
  /* Now read in the ngc catalog */
  printf("Opening %s\n",ngc_name);
  input_handle = fopen(ngc_name,"rt");
  if (input_handle == NULL) {
    printf("ERROR: failed to open %s\n",ngc_name);
    exit(-1);
  }
  numLines = 0;
  while(1) {
    memset(pNgcImage,0,sizeof(NGCIMAGE));
    inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
    if (inBuffer == NULL) {
      break;
    }
    strcpy(copyLine,inBuffer);
    numLines++;
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
    if (lineLen < offsetof(NGCIMAGE,Desc)) {
      printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
      continue;
    }
    /* Pad the line with blanks */
    while (lineLen < sizeof(NGCIMAGE)) {
      inBuffer[lineLen] = ' ';
      inBuffer[lineLen+1] = 0;
      lineLen++;
    }
    memcpy(pNgcImage,inBuffer,sizeof(NGCIMAGE));
    if ((pNgcImage->blank1 != ' ') ||
        (pNgcImage->blank2 != ' ') ||
        (pNgcImage->blank3 != ' ') ||
        (pNgcImage->blank4 != ' ') ||
        (pNgcImage->blank5 != ' ') ||
        (pNgcImage->blank6 != ' ') ||
        (pNgcImage->blank7 != ' ') ||
        (pNgcImage->blank8 != ' ') ||
        (pNgcImage->blank9 != ' ') ||
        (pNgcImage->blank10 != ' ') ||
        (pNgcImage->blank11 != ' ') ||
        (pNgcImage->blank12 != ' ')) {
      printf("ERROR: line '%s' has an invalid record type\n",copyLine);
      continue;
    }
    pOutRec = &output_table[totalCount];
    memset(pOutRec,0,sizeof(OUTREC));
    pOutRec->source = SOURCE_NGC;
    pNgcImage->blank2 = 0;
#if 0
    strcpy(pOutRec->type,"NGC2000");
#else
    strcpy(pOutRec->type,pNgcImage->Type);
#endif
    pOutRec->catalogmag = 99.;
    pNgcImage->blank1 = 0;
    trimblank(pNgcImage->Name);
    if (pNgcImage->Name[0] == 'I') {
      nvals = sscanf(&pNgcImage->Name[1],"%d",&ngcnumber);
      if (nvals == 1) {
        sprintf(pOutRec->catalogname,"IC%d",ngcnumber);
      } else {
        printf("ERROR: failed to decode record %d NGC number %s\n",numLines,pNgcImage->Name);
        continue;
      }

    } else {
      nvals = sscanf(pNgcImage->Name,"%d",&ngcnumber);
      if (nvals == 1) {
        sprintf(pOutRec->catalogname,"NGC%d",ngcnumber);
      } else {
        printf("ERROR: failed to decode record %d NGC number %s\n",numLines,pNgcImage->Name);
        continue;
      }
    }
#if 0

    if (strcmp(pOutRec->catalogname,"N5890") == 0) {
      printf("At %s\n",pOutRec->catalogname);
    }
#endif

    pNgcImage->blank4 = 0;
    nvals = sscanf(pNgcImage->RAh,"%2d %lf",&degrees,&realminutes);
    if (nvals == 2) {
      pOutRec->ra = 15.0*((1.0*degrees)+((1.0*realminutes)/60.0));
    } else {
      printf("ERROR: failed to decode record %d ra %s\n",numLines,pNgcImage->RAh);
      continue;
    }
    pNgcImage->blank7 = 0;
    nvals = sscanf(pNgcImage->DEd,"%2d %d",&degrees,&minutes);
    if (nvals == 2) {
      pOutRec->dec = ((1.0*degrees)+((1.0*minutes)/60.0));
      if (pNgcImage->DE1 == '-') {
        pOutRec->dec = -pOutRec->dec;
      }
    } else {
      printf("ERROR: failed to decode record %d dec %s\n",numLines,pNgcImage->DEd);
      continue;
    }
#if 0
    RAold = pOutRec->ra;
    decold = pOutRec->dec;
    
    wcscon(WCS_B1950,WCS_J2000,2000.0,2000.0,&pOutRec->ra,&pOutRec->dec,2000.0);
    printf("B2000 ra %f dec %f  J2000 ra %f dec %f\n",RAold,decold,pOutRec->ra,pOutRec->dec);
#endif

    pNgcImage->blank12 = 0;
    nvals = sscanf(pNgcImage->mag,"%lf",&pOutRec->catalogmag);
    if (nvals != 1) {
      pOutRec->catalogmag = 99.0;
    }
    pNgcImage->blank10 = 0;
    nvals = sscanf(pNgcImage->size,"%lf",&radius);
    if (nvals == 1) {
      pOutRec->radius = 60.0*radius/2.0; /* Radius in arcsec */
    }
#if 0
    printf("line %d name %s type %s ra %f dec %f mag %f radius %f\n",
           __LINE__,pOutRec->catalogname,pOutRec->type,pOutRec->ra,pOutRec->dec,pOutRec->catalogmag,pOutRec->radius);
#endif           
    if (pOutRec->radius < NGC_MINIMUM_RADIUS) {
      pOutRec->radius = NGC_MINIMUM_RADIUS;
    }
    totalCount++;
#if 0
    printf("line %d %s name %s type %s ra %f dec %f mag %f radius %f\n",
           __LINE__,copyLine,pOutRec->catalogname,pOutRec->type,pOutRec->ra,pOutRec->dec,pOutRec->catalogmag,pOutRec->radius);
#endif           
    if ((totalCount) >= MAX_INPUT_LINES) {
      printf("ERROR: MAX_INPUT_LINES %d exceeded line %d\n",MAX_INPUT_LINES,__LINE__);
      exit(-1);
    }

  }
  printf("Read %d lines from %s\n",numLines,ngc_name);
  fclose(input_handle);

  if (ngcOnlyFlag == 0) {
    /* Now read in the UGC1 catalog.  Data for this catalog goes in a tmporary holding array until the associated record is read from the UGC2 catalog */
    printf("Opening %s\n",ugc1_name);
    input_handle = fopen(ugc1_name,"rt");
    if (input_handle == NULL) {
      printf("ERROR: failed to open %s\n",ugc1_name);
      exit(-1);
    }
    numLines = 0;
    while(1) {
      memset(pUgc1Image,0,sizeof(UGC1IMAGE));
      inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
      if (inBuffer == NULL) {
        break;
      }
      strcpy(copyLine,inBuffer);
      numLines++;
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
      if (lineLen < offsetof(UGC1IMAGE,Hubble)) {
        printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
        continue;
      }
      /* Pad the line with blanks */
      while (lineLen < sizeof(UGC1IMAGE)) {
        inBuffer[lineLen] = ' ';
        inBuffer[lineLen+1] = 0;
        lineLen++;
      }
      memcpy(pUgc1Image,inBuffer,sizeof(UGC1IMAGE));
      if ((pUgc1Image->blank01 != ' ') ||
          (pUgc1Image->blank02 != ' ')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      if ((pUgc1Image->head[0] != 'U') ||
          (pUgc1Image->head[1] != 'G') ||
          (pUgc1Image->head[2] != 'C')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      pOutUgc= &ugc_table[ugcCount];
      memset(pOutUgc,0,sizeof(OUTREC));
      pOutUgc->suffix = pUgc1Image->A;
      pUgc1Image->A = 0;
      nvals = sscanf(pUgc1Image->UGC,"%d",&pOutUgc->ugc);
      if (nvals != 1) {
        printf("ERROR: failed to decode UGC designator in %s\n",copyLine);
        continue;
      }
      pUgc1Image->RadVal[0] = 0;
      nvals = sscanf(pUgc1Image->mag,"%lf",&pOutUgc->catalogmag);
      if (nvals != 1) {
        pOutUgc->catalogmag = 99;
      }
      pUgc1Image->MinAxis[0] = 0;
      nvals = sscanf(pUgc1Image->MajAxis,"%lf",&pOutUgc->radius);
      if (nvals != 1) {
        pOutUgc->radius = 0;
      }
      pOutUgc->radius = (60.0*pOutUgc->radius)/2.0;
      ugcCount++;
      if ((ugcCount) >= MAX_UGC_RECORDS) {
        printf("ERROR: MAX_UGC_RECORDS %d exceeded\n",MAX_UGC_RECORDS);
        exit(-1);
      }
    }
    printf("Read %d lines from %s\n",numLines,ugc1_name);
    fclose(input_handle);

#if 0
    for (ugcIndex = 0; ugcIndex < ugcCount; ugcIndex++) {
      pOutUgc= &ugc_table[ugcIndex];
      printf("Index %5d designator %05d%c magnitude %4.1f diameter %6.2f (arcmin)\n",
             ugcIndex,
             pOutUgc->ugc,
             pOutUgc->suffix,
             pOutUgc->catalogmag,
             (2*(pOutUgc->radius/60.0)));
    }
#endif

    /* Now read in the UGC2 catalog */
    printf("Opening %s\n",ugc2_name);
    input_handle = fopen(ugc2_name,"rt");
    if (input_handle == NULL) {
      printf("ERROR: failed to open %s\n",ugc2_name);
      exit(-1);
    }
    numLines = 0;
    while(1) {
      memset(pUgc2Image,0,sizeof(UGC2IMAGE));
      inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
      if (inBuffer == NULL) {
        break;
      }
      strcpy(copyLine,inBuffer);
      numLines++;
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
      if (lineLen < offsetof(UGC2IMAGE,e_DE1s)) {
        printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
        continue;
      }
      /* Pad the line with blanks */
      while (lineLen < sizeof(UGC2IMAGE)) {
        inBuffer[lineLen] = ' ';
        inBuffer[lineLen+1] = 0;
        lineLen++;
      }
      memcpy(pUgc2Image,inBuffer,sizeof(UGC2IMAGE));
      if ((pUgc2Image->blank01 != ' ') ||
          (pUgc2Image->blank03 != ' ') ||
          (pUgc2Image->blank04 != ' ') ||
          (pUgc2Image->blank05 != ' ') ||
          (pUgc2Image->blank06 != ' ') ||
          (pUgc2Image->blank07 != ' ') ||
          (pUgc2Image->blank08 != ' ') ||
          (pUgc2Image->blank09 != ' ') ||
          (pUgc2Image->blank10 != ' ') ||
          (pUgc2Image->blank11 != ' ') ||
          (pUgc2Image->blank12 != ' ') ||
          (pUgc2Image->blank13 != ' ') ||
          (pUgc2Image->blank14 != ' ') ||
          (pUgc2Image->blank15 != ' ') ||
          (pUgc2Image->blank16 != ' ') ||
          (pUgc2Image->blank17 != ' ') ||
          (pUgc2Image->blank18 != ' ') ||
          (pUgc2Image->blank19 != ' ') ||
          (pUgc2Image->blank20 != ' ') ||
          (pUgc2Image->blank21 != ' ') ||
          (pUgc2Image->blank22 != ' ') ||
          (pUgc2Image->blank23 != ' ') ||
          (pUgc2Image->blank24 != ' ') ||
          (pUgc2Image->blank25 != ' ') ||
          (pUgc2Image->blank26 != ' ') ||
          (pUgc2Image->blank27 != ' ') ||
          (pUgc2Image->blank28 != ' ') ||
          (pUgc2Image->blank29 != ' ')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      if ((pUgc2Image->head[0] != 'U') ||
          (pUgc2Image->head[1] != 'G') ||
          (pUgc2Image->head[2] != 'C')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      pOutRec = &output_table[totalCount];
      memset(pOutRec,0,sizeof(OUTREC));
      pOutRec->source = SOURCE_UGC ;
      strcpy(pOutRec->type,"UGC Galaxy");
      pOutRec->catalogmag = 99.;
      ugcsuffix = pUgc2Image->m_UGC[0];
      ugcsuffix2 = pUgc2Image->m_UGC[1];
      pUgc2Image->m_UGC[0] = 0;
      nvals = sscanf(pUgc2Image->UGC,"%d",&ugc);
      if (nvals != 1) {
        printf("ERROR: failed to decode UGC designator in %s\n",copyLine);
        continue;
      }
      sprintf(pOutRec->catalogname,"UGC%05d%c%c",ugc,ugcsuffix,ugcsuffix2);
      trimblank(pOutRec->catalogname);
#if 0
      if (ugc == 4280) {
        printf("At %d%c%c\n",ugc,ugcsuffix,ugcsuffix2);
      }
#endif
      pUgc2Image->blank10 = 0;
      nvals = sscanf(pUgc2Image->RAh,"%2d %2d %lf",&degrees,&minutes,&seconds);
      if (nvals == 3) {
        pOutRec->ra = 15.0*((1.0*degrees)+((1.0*minutes)/60.0)+((1.0*seconds)/3600.));
      } else {
        printf("ERROR: failed to decode record %d ra %s\n",numLines,pUgc2Image->RAh);
        continue;
      }
      pUgc2Image->blank15 = 0;
      nvals = sscanf(pUgc2Image->DEd,"%2d %2d %lf",&degrees,&minutes,&seconds);
      if (nvals == 3) {
        pOutRec->dec = ((1.0*degrees)+((1.0*minutes)/60.0)) +(seconds/3600.);
        if (pUgc2Image->DE == '-') {
          pOutRec->dec = -pOutRec->dec;
        }
      } else {
        printf("ERROR: failed to decode record %d dec %s\n",numLines,pUgc2Image->DE);
        continue;
      }
      /* At this point, search the temporary list for the object and copy over the size and magnitude */
#if 0
      if (ugc == 7064) {
        printf("At %05d%c%c\n",ugc,ugcsuffix,ugcsuffix2);
      }
#endif

      for (ugcIndex = 0; ugcIndex < ugcCount; ugcIndex++) {
        pOutUgc= &ugc_table[ugcIndex];
        if (ugc == pOutUgc->ugc) {
          if (((ugcsuffix == 'A') && (ugcsuffix == pOutUgc->suffix)) ||
              ((ugcsuffix != 'A') && (pOutUgc->suffix == ' '))) {
#if 0
            if ((ugcsuffix != ' ') ||
                (pOutUgc->suffix != ' ')) {
              printf("Matched %5d%c%c and %5d%c\n",ugc,ugcsuffix,ugcsuffix2,pOutUgc->ugc,pOutUgc->suffix);
            }
#endif
            pOutRec->catalogmag = pOutUgc->catalogmag;
            pOutRec->radius = pOutUgc->radius;
            break;
          }


        }
        
      }



      totalCount++;
      if ((totalCount) >= MAX_INPUT_LINES) {
        printf("ERROR: MAX_INPUT_LINES %d exceeded line %d\n",MAX_INPUT_LINES,__LINE__);
        exit(-1);
      }
    }
    printf("Read %d lines from %s\n",numLines,ugc2_name);
    fclose(input_handle);

    /* Now read in the OGLE identification catalog.  Data for this catalog goes in a tmporary holding array until the associated record is read from the OGLE ecl or ell catalogs */
    printf("Opening %s\n",ogleident_name);
    input_handle = fopen(ogleident_name,"rt");
    if (input_handle == NULL) {
      printf("ERROR: failed to open %s\n",ogleident_name);
      exit(-1);
    }
    numLines = 0;
    while(1) {
      memset(pOgleIdentImage,0,sizeof(OGLEIDENTIMAGE));
      inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
      if (inBuffer == NULL) {
        break;
      }
      strcpy(copyLine,inBuffer);
      numLines++;
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
      if (lineLen < offsetof(OGLEIDENTIMAGE,other)) {
        printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
        continue;
      }
      /* Pad the line with blanks */
      while (lineLen < sizeof(OGLEIDENTIMAGE)) {
        inBuffer[lineLen] = ' ';
        inBuffer[lineLen+1] = 0;
        lineLen++;
      }
      memcpy(pOgleIdentImage,inBuffer,sizeof(OGLEIDENTIMAGE));
      if ((pOgleIdentImage->blank01[0] != ' ') ||
          (pOgleIdentImage->blank01[1] != ' ') ||
          (pOgleIdentImage->blank02 != ' ') ||
          (pOgleIdentImage->blank03 != ':') ||
          (pOgleIdentImage->blank04 != ':') ||
          (pOgleIdentImage->blank05 != ' ') ||
          (pOgleIdentImage->blank06 != ':') ||
          (pOgleIdentImage->blank07 != ':') ||
          (pOgleIdentImage->blank08[0] != ' ') ||
          (pOgleIdentImage->blank08[1] != ' ') ||
          (pOgleIdentImage->blank09 != ' ') ||
          (pOgleIdentImage->blank10 != ' ')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      pOutOgle= &ogle_table[ogleCount];
      memset(pOutOgle,0,sizeof(OUTREC));
      if (strstr(&pOgleIdentImage->ogleid[0],"OGLE-BLG-ECL-") == &pOgleIdentImage->ogleid[0]) {
        pOutOgle->ECLflag = 1;
      } else if (strstr(&pOgleIdentImage->ogleid[0],"OGLE-BLG-ELL-") == &pOgleIdentImage->ogleid[0]) {
        pOutOgle->ECLflag = 0;
      } else {
        printf("ERROR: line '%s' has an invalid OGLE id\n",copyLine);
        continue;
      }
      nvals = sscanf(&pOgleIdentImage->ogleid[13],"%d",&pOutOgle->idnumber);
      if (nvals != 1) {
        printf("ERROR line '%s' can not decode OGLE id number in %s\n",pOgleIdentImage->ogleid);
        continue;
      }
      /* The order of the ident file should be first ECL records numbered 1-425193 next ELL records numbered 1-25405 */
      if (pOutOgle->ECLflag == 1) {
        if (pOutOgle->idnumber != ogleCount+1) {
          printf("ERROR: ECL idnumber %d does not match record index %d\n",pOutOgle->idnumber,ogleCount);
        }
      } else {
        if ((pOutOgle->idnumber+OGLE_ELL_BASE) != ogleCount) {
          printf("ERROR: ECL idnumber %d (+OGLE_ELL_BASE = %d) does not match record index %d \n",pOutOgle->idnumber,pOutOgle->idnumber+OGLE_ELL_BASE,ogleCount);
          exit(-1);
        }
      }
      pOgleIdentImage->blank01[0] = 0;
      if (strlen(pOgleIdentImage->ogleid) > sizeof(pOutOgle->ogleid)) {
        printf("ERROR: line '%s' has an OGLE id that is too long\n",copyLine);
        exit(-1);
      }
      strcpy(pOutOgle->ogleid,pOgleIdentImage->ogleid);

      if (pOgleIdentImage->subtype[1] == ' ') {
        pOgleIdentImage->subtype[1] = 0;
      } else if (pOgleIdentImage->subtype[2] == ' ') {
        pOgleIdentImage->subtype[2] = 0;
      } else {
        pOgleIdentImage->blank02 = 0;
      }
      if ((strcmp(pOgleIdentImage->subtype,"C") != 0) &&
          (strcmp(pOgleIdentImage->subtype,"NC") != 0) &&
          (strcmp(pOgleIdentImage->subtype,"CV") != 0) &&
          (strcmp(pOgleIdentImage->subtype,"ELL") != 0)) {
        printf("ERROR: line '%s' has an invalid OGLE subtype %s\n",copyLine,pOgleIdentImage->subtype);
        continue;
      }
      strcpy(pOutOgle->subtype,pOgleIdentImage->subtype);
      pOgleIdentImage->blank05 = 0;
      nvals = sscanf(pOgleIdentImage->RAh,"%2d:%2d:%lf",&degrees,&minutes,&seconds);
      if (nvals == 3) {
        pOutOgle->ra = 15.0*((1.0*degrees)+((1.0*minutes)/60.0)+((1.0*seconds)/3600.));
      } else {
        printf("ERROR: failed to decode record %d ra %s\n",numLines,pOgleIdentImage->RAh);
        continue;
      }

      pOgleIdentImage->blank08[0] = 0;
      nvals = sscanf(pOgleIdentImage->DEd,"%2d:%2d:%lf",&degrees,&minutes,&seconds);
      if (nvals == 3) {
        pOutOgle->dec = ((1.0*degrees)+((1.0*minutes)/60.0)) +(seconds/3600.);
        if (pOgleIdentImage->DE == '-') {
          pOutOgle->dec = -pOutOgle->dec;
        }
      } else {
        printf("ERROR: failed to decode record %d dec %s\n",numLines,pOgleIdentImage->DE);
        continue;
      }
      ogleCount++;
      if ((ogleCount) >= MAX_OGLE_RECORDS) {
        printf("ERROR: MAX_OGLE_RECORDS %d exceeded\n",MAX_OGLE_RECORDS);
        exit(-1);
      }
    }

    printf("Read %d lines from %s\n",numLines,ogleident_name);
    fclose(input_handle);

    /* Now read in the OGLE ecl catalog */
    printf("Opening %s\n",ogleecl_name);
    input_handle = fopen(ogleecl_name,"rt");
    if (input_handle == NULL) {
      printf("ERROR: failed to open %s\n",ogleecl_name);
      exit(-1);
    }
    numLines = 0;
    while(1) {
      memset(pOgleDatImage,0,sizeof(OGLEDATIMAGE));
      inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
      if (inBuffer == NULL) {
        break;
      }
      strcpy(copyLine,inBuffer);
      numLines++;
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
      if (lineLen < offsetof(OGLEDATIMAGE,secdepth)) {
        printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
        continue;
      }
      /* Pad the line with blanks */
      while (lineLen < sizeof(OGLEDATIMAGE)) {
        inBuffer[lineLen] = ' ';
        inBuffer[lineLen+1] = 0;
        lineLen++;
      }
      memcpy(pOgleDatImage,inBuffer,sizeof(OGLEDATIMAGE));
      if ((pOgleDatImage->blank01[0] != ' ') ||
          (pOgleDatImage->blank01[1] != ' ') ||
          (pOgleDatImage->blank02 != ' ') ||
          (pOgleDatImage->blank03 != ' ') ||
          (pOgleDatImage->blank04[0] != ' ') ||
          (pOgleDatImage->blank04[1] != ' ') ||
          (pOgleDatImage->blank05[0] != ' ') ||
          (pOgleDatImage->blank05[1] != ' ') ||
          (pOgleDatImage->blank06 != ' ')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      pOgleDatImage->blank01[0] = 0;
      if (strstr(&pOgleDatImage->ogleid[0],"OGLE-BLG-ECL-") == &pOgleDatImage->ogleid[0]) {
        pOutOgleTemp->ECLflag = 1;
      } else if (strstr(&pOgleDatImage->ogleid[0],"OGLE-BLG-ELL-") == &pOgleDatImage->ogleid[0]) {
        pOutOgleTemp->ECLflag = 0;
      } else {
        printf("ERROR: line '%s' has an invalid OGLE id\n",copyLine);
        continue;
      }
      nvals = sscanf(&pOgleDatImage->ogleid[13],"%d",&pOutOgleTemp->idnumber);
      if (nvals != 1) {
        printf("ERROR line '%s' can not decode OGLE id number in %s\n",pOgleDatImage->ogleid);
        continue;
      }


      /* Now hunt for the ID in the ident table */
      if (pOutOgleTemp->ECLflag == 1) {
        ogleIndex = pOutOgleTemp->idnumber-1;
      } else {
        ogleIndex = pOutOgleTemp->idnumber+OGLE_ELL_BASE;
      }
      if ((ogleIndex < 0) || (ogleIndex >= ogleCount)) {
        printf("ERROR: line '%s' has an unknown ogle index %d\n",copyLine,ogleIndex);
        exit(-1);
      }


      pOutOgle= &ogle_table[ogleIndex];
      if ((pOutOgle->idnumber != pOutOgleTemp->idnumber) ||
          (pOutOgle->ECLflag != pOutOgleTemp->ECLflag)) {
        printf("ERROR bad OGLE index %d %d for ECLflag %d %d line %s\n",pOutOgle->idnumber,pOutOgleTemp->idnumber,pOutOgle->ECLflag,pOutOgleTemp->ECLflag,copyLine);
        exit(-1);
      }
      if (pOutOgle->matchFlag != 0) {
        printf("ERROR: matchFlag not zero for OGLE index %d in line %s\n",copyLine);
        exit(-1);
      }
      pOutOgle->matchFlag = 1;
      ogleMatch++;

      pOgleDatImage->blank02 = 0;
      nvals = sscanf(&pOgleDatImage->imaxmag[0],"%lf",&pOutOgle->imaxmag);
      if (nvals != 1) {
        printf("ERROR decoding imaxmag in %s\n",copyLine);
        continue;
      }

      pOgleDatImage->blank03 = 0;
      if (strcmp(pOgleDatImage->vmaxmag,"   -  ") == 0) {
        pOutOgle->vmaxmag = 99.0;
        ogleBadVMaxMagCount++;
      } else {
        nvals = sscanf(&pOgleDatImage->vmaxmag[0],"%lf",&pOutOgle->vmaxmag);
        if (nvals != 1) {
          printf("ERROR decoding vmaxmag in %s\n",copyLine);
          continue;
        }
      }
      pOgleDatImage->blank04[0] = 0;
      if (strcmp(pOgleDatImage->period,"     -      ") == 0) {
        pOutOgle->period = 0;
        ogleBadPeriodCount++;
      } else {
        nvals = sscanf(&pOgleDatImage->period[0],"%lf",&pOutOgle->period);
        if (nvals != 1) {
          printf("ERROR decoding period in %s\n",copyLine);
          continue;
        }
      }
      pOgleDatImage->blank05[0] = 0;
      nvals = sscanf(&pOgleDatImage->epoch[0],"%lf",&pOutOgle->epoch);
      if (nvals != 1) {
        printf("ERROR decoding epoch in %s\n",copyLine);
        continue;
      }

      pOgleDatImage->blank06 = 0;
      nvals = sscanf(&pOgleDatImage->pridepth[0],"%lf",&pOutOgle->pridepth);
      if (nvals != 1) {
        printf("ERROR decoding pridepth in %s\n",copyLine);
        continue;
      }

      pOgleDatImage->blank07 = 0;
      nvals = sscanf(&pOgleDatImage->secdepth[0],"%lf",&pOutOgle->secdepth);
      if (nvals != 1) {
        printf("ERROR decoding secdepth in %s\n",copyLine);
        continue;
      }


#if 0 /* This takes 593 sec for 425,193 records */
      if (strcmp(pOutOgle->ogleid,pOgleDatImage->ogleid) == 0) {
        break;
      }
#endif
    
    }
    printf("Read %d lines from %s bad vmaxmag %d bad period %d\n",numLines,ogleecl_name,ogleBadVMaxMagCount,ogleBadPeriodCount);
    fclose(input_handle);

    /* Now read in the OGLE ell catalog */
    printf("Opening %s\n",ogleell_name);
    input_handle = fopen(ogleell_name,"rt");
    if (input_handle == NULL) {
      printf("ERROR: failed to open %s\n",ogleell_name);
      exit(-1);
    }
    numLines = 0;
    while(1) {
      memset(pOgleDatImage,0,sizeof(OGLEDATIMAGE));
      inBuffer = fgets(inLine,MAX_BUFFER,input_handle);
      if (inBuffer == NULL) {
        break;
      }
      strcpy(copyLine,inBuffer);
      numLines++;
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
      if (lineLen < offsetof(OGLEDATIMAGE,secdepth)) {
        printf("ERROR: line '%s' is too small at length %d\n",copyLine,lineLen);
        continue;
      }
      /* Pad the line with blanks */
      while (lineLen < sizeof(OGLEDATIMAGE)) {
        inBuffer[lineLen] = ' ';
        inBuffer[lineLen+1] = 0;
        lineLen++;
      }
      memcpy(pOgleDatImage,inBuffer,sizeof(OGLEDATIMAGE));
      if ((pOgleDatImage->blank01[0] != ' ') ||
          (pOgleDatImage->blank01[1] != ' ') ||
          (pOgleDatImage->blank02 != ' ') ||
          (pOgleDatImage->blank03 != ' ') ||
          (pOgleDatImage->blank04[0] != ' ') ||
          (pOgleDatImage->blank04[1] != ' ') ||
          (pOgleDatImage->blank05[0] != ' ') ||
          (pOgleDatImage->blank05[1] != ' ') ||
          (pOgleDatImage->blank06 != ' ')) {
        printf("ERROR: line '%s' has an invalid record type\n",copyLine);
        continue;
      }
      pOgleDatImage->blank01[0] = 0;
      if (strstr(&pOgleDatImage->ogleid[0],"OGLE-BLG-ECL-") == &pOgleDatImage->ogleid[0]) {
        pOutOgleTemp->ECLflag = 1;
      } else if (strstr(&pOgleDatImage->ogleid[0],"OGLE-BLG-ELL-") == &pOgleDatImage->ogleid[0]) {
        pOutOgleTemp->ECLflag = 0;
      } else {
        printf("ERROR: line '%s' has an invalid OGLE id\n",copyLine);
        continue;
      }
      nvals = sscanf(&pOgleDatImage->ogleid[13],"%d",&pOutOgleTemp->idnumber);
      if (nvals != 1) {
        printf("ERROR line '%s' can not decode OGLE id number in %s\n",pOgleDatImage->ogleid);
        continue;
      }


      /* Now hunt for the ID in the ident table */
      if (pOutOgleTemp->ECLflag == 1) {
        ogleIndex = pOutOgleTemp->idnumber-1;
      } else {
        ogleIndex = pOutOgleTemp->idnumber+OGLE_ELL_BASE;
      }
      if ((ogleIndex < 0) || (ogleIndex >= ogleCount)) {
        printf("ERROR: line '%s' has an unknown ogle index %d\n",copyLine,ogleIndex);
        exit(-1);
      }


      pOutOgle= &ogle_table[ogleIndex];
      if ((pOutOgle->idnumber != pOutOgleTemp->idnumber) ||
          (pOutOgle->ECLflag != pOutOgleTemp->ECLflag)) {
        printf("ERROR bad OGLE index %d %d for ECLflag %d %d line %s\n",pOutOgle->idnumber,pOutOgleTemp->idnumber,pOutOgle->ECLflag,pOutOgleTemp->ECLflag,copyLine);
        exit(-1);
      }
      if (pOutOgle->matchFlag != 0) {
        printf("ERROR: matchFlag not zero for OGLE line %s\n",copyLine);
        exit(-1);
      }
      pOutOgle->matchFlag = 1;
      ogleMatch++;

      pOgleDatImage->blank02 = 0;
      nvals = sscanf(&pOgleDatImage->imaxmag[0],"%lf",&pOutOgle->imaxmag);
      if (nvals != 1) {
        printf("ERROR decoding imaxmag in %s\n",copyLine);
        continue;
      }

      pOgleDatImage->blank03 = 0;
      if (strcmp(pOgleDatImage->vmaxmag,"   -  ") == 0) {
        pOutOgle->vmaxmag = 99.0;
        ogleBadVMaxMagCount++;
      } else {
        nvals = sscanf(&pOgleDatImage->vmaxmag[0],"%lf",&pOutOgle->vmaxmag);
        if (nvals != 1) {
          printf("ERROR decoding vmaxmag in %s\n",copyLine);
          continue;
        }
      }
      pOgleDatImage->blank04[0] = 0;
      if (strcmp(pOgleDatImage->period,"     -      ") == 0) {
        pOutOgle->period = 0;
        ogleBadPeriodCount++;
      } else {
        nvals = sscanf(&pOgleDatImage->period[0],"%lf",&pOutOgle->period);
        if (nvals != 1) {
          printf("ERROR decoding period in %s\n",copyLine);
          continue;
        }
      }
      pOgleDatImage->blank05[0] = 0;
      nvals = sscanf(&pOgleDatImage->epoch[0],"%lf",&pOutOgle->epoch);
      if (nvals != 1) {
        printf("ERROR decoding epoch in %s\n",copyLine);
        continue;
      }

      pOgleDatImage->blank06 = 0;
      nvals = sscanf(&pOgleDatImage->pridepth[0],"%lf",&pOutOgle->pridepth);
      if (nvals != 1) {
        printf("ERROR decoding pridepth in %s\n",copyLine);
        continue;
      }

      pOgleDatImage->blank07 = 0;
      nvals = sscanf(&pOgleDatImage->secdepth[0],"%lf",&pOutOgle->secdepth);
      if (nvals != 1) {
        printf("ERROR decoding secdepth in %s\n",copyLine);
        continue;
      }

    
    }
    printf("Read %d lines from %s bad vmaxmag %d bad period %d\n",numLines,ogleell_name,ogleBadVMaxMagCount,ogleBadPeriodCount);
    fclose(input_handle);
    if (ogleCount != ogleMatch) {
      printf("ERROR: ogleCount %d does not equal ogleMatch %d\n",ogleCount,ogleMatch);
      exit(-1);
    }

    for (ogleIndex = 0; ogleIndex < ogleCount; ogleIndex++) {
      pOutOgle= &ogle_table[ogleIndex];
      pOutRec = &output_table[totalCount];
      memset(pOutRec,0,sizeof(OUTREC));
      pOutRec->source = SOURCE_OGLE;
      pOutRec->ra = pOutOgle->ra;
      pOutRec->dec = pOutOgle->dec;
      pOutRec->catalogmag = pOutOgle->vmaxmag;
      pOutRec->radius = 0;
      strcpy(pOutRec->catalogname,pOutOgle->ogleid);
      sprintf(pOutRec->type,"OGLE:%s",pOutOgle->subtype);
      totalCount++;
      if ((totalCount) >= MAX_INPUT_LINES) {
        printf("ERROR: MAX_INPUT_LINES %d exceeded line %d\n",MAX_INPUT_LINES,__LINE__);
        exit(-1);
      }
  
    }
#if 0
    printf("variableName\tvariableDescription\tra\tdec\timaxmag\tvmaxmag\tperiod\thJD\tpridepth\tsecdepth\n");
    printf("------------\t-------------------\t--\t---\t-------\t-------\t------\t---\t--------\t--------\n");
    for (ogleIndex = 0; ogleIndex < ogleCount; ogleIndex++) {
      pOutOgle= &ogle_table[ogleIndex];
      printf("%s\t%s\t%.5f\t%.5f\t%.3f\t%.3f\t%.7f\t%.4f\t%.3f\t%.3f\n",
             pOutOgle->ogleid,
             pOutOgle->subtype,
             pOutOgle->ra,
             pOutOgle->dec,
             pOutOgle->imaxmag,
             pOutOgle->vmaxmag,
             pOutOgle->period,
             pOutOgle->epoch+2450000.0,
             pOutOgle->pridepth,
             pOutOgle->secdepth);

    }
#endif
  } /* ngcOnlyFlag */
  printf("totalCount is %d\n",totalCount);

	time(&curTime);
	curTime -= startTime;
  printf("Populating vsx/ogle kdtree at %d seconds\n",curTime);
  /* create a k-d tree for 3-dimensional points */
  ptree = kd_create( 3 );
  for (outputIndex = 0; outputIndex < totalCount; outputIndex++) {
    pOutRec = &output_table[outputIndex];
    if ((pOutRec->source == SOURCE_VSX)  ||
        (pOutRec->source == SOURCE_OGLE)) {
      z = sin(DEGREES_TO_RAD*pOutRec->dec);
      cosx = cos(DEGREES_TO_RAD*pOutRec->dec);
      x = cosx*cos(DEGREES_TO_RAD*pOutRec->ra);
      y = cosx*sin(DEGREES_TO_RAD*pOutRec->ra);
      if (kd_insert3(ptree,x,y,z,pOutRec) != 0) {
        printf("ERROR: fatal return from kd_insert3\n");
        exit(-1);
      }
    }
  }
 
	time(&curTime);
	curTime -= startTime;
  printf("Checking for vsx/ogle duplicates at %d seconds\n",curTime);
  /* We start first with the variables from VSX and OGLE */
  minDistance = GALAXY_SEARCH_RADIUS*3600.;
  for (outputIndex = 0; outputIndex < totalCount; outputIndex++) {
    pOutRec = &output_table[outputIndex];
    /* Consider VSX or OGLE sources */
    if (pOutRec->source == SOURCE_VSX) {
      pOutRec->checked = 1;
      z = sin(DEGREES_TO_RAD*pOutRec->dec);
      cosx = cos(DEGREES_TO_RAD*pOutRec->dec);
      x = cosx*cos(DEGREES_TO_RAD*pOutRec->ra);
      y = cosx*sin(DEGREES_TO_RAD*pOutRec->ra);
      pt[0] = x;
      pt[1] = y;
      pt[2] = z;
      presults = kd_nearest_range( ptree, pt, VARIABLE_SEARCH_RADIUS*DEGREES_TO_RAD );
      resultSize = kd_res_size(presults);
      if (resultSize > maxResultSize) {
        maxResultSize = resultSize;
      }
      filterSize = 0;
      while( !kd_res_end( presults ) ) {
        /* get the data and position of the current result item */
        pChkRec = (POUTREC)kd_res_item( presults, pos );
        if ((pChkRec == pOutRec) ||
            (pChkRec->checked != 0) ||
            (pChkRec->duplicate != 0) ||
            (pChkRec->source == SOURCE_VSX)) {
          rejectCount++;
          /* go to the next entry */
          kd_res_next( presults );
          continue;
        }
        distance = 3600.0*RAD_TO_DEGREES*sqrt(sqr(pt[0]-pos[0])+sqr(pt[1]-pos[1])+sqr(pt[2]-pos[2]));
        if (distance < minDistance) {
          minDistance = distance;
        }
        pFilterRec = &filter_table[filterSize];
        memcpy(pFilterRec,pChkRec,sizeof(OUTREC));
        pFilterRec->pOutRec = pChkRec;
        pFilterRec->distance = distance;
        filterSize++;
        if (filterSize >= MAX_FILTER_SIZE) {
          printf("ERROR: MAX_FILTER_SIZE %d is too small\n",MAX_FILTER_SIZE);
          exit(-1);
        }
#if 0
        printf("%12.6f arcsec filter %3d %31s %31s is near %31s %31s\n",distance,filterSize,pChkRec->catalogname,pChkRec->type,pOutRec->catalogname,pOutRec->type);
#endif
#if 0
        printf("dx %9.6f y %9.6f dz %9.6f distance %9.6f %s %s\n",pt[0]-pos[0],pt[1]-pos[1],pt[2]-pos[2],distance,pChkRec->catalogname,pChkRec->type);
#endif
        /* go to the next entry */
        kd_res_next( presults );
      }
      CheckDuplicates(pOutRec,filter_table,filterSize,cursource,&radiusCount,&catalogCount,&duplicateCount,1);
      kd_res_free(presults);
      if (filterSize > maxFilterSize) {
        maxFilterSize = filterSize;
      }
      pOutRec->checked = 1;

      
    }
  }
  printf("OGLE variables: minDistance %f maxFilterSize %d\n",minDistance,maxFilterSize);
  printf("Done checking for duplicates rejectCount %d radiusReject %d catalogReject %d duplicates %d maxResultSize %d maxFilterSize %d at %d seconds\n",
         rejectCount,
         radiusCount,
         catalogCount,
         duplicateCount,
         maxResultSize,
         maxFilterSize,
         curTime);

  kd_free(ptree);


  ptree = kd_create(3);
	time(&curTime);
	curTime -= startTime;
  printf("Populating kdtree at %d seconds\n",curTime);
  /* create a k-d tree for 3-dimensional points */
  ptree = kd_create( 3 );
  for (outputIndex = 0; outputIndex < totalCount; outputIndex++) {
    pOutRec = &output_table[outputIndex];
    /* This time, ignore VSX and OGLE variables */
    if ((pOutRec->source != SOURCE_VSX)  &&
        (pOutRec->source != SOURCE_OGLE)) {
      z = sin(DEGREES_TO_RAD*pOutRec->dec);
      cosx = cos(DEGREES_TO_RAD*pOutRec->dec);
      x = cosx*cos(DEGREES_TO_RAD*pOutRec->ra);
      y = cosx*sin(DEGREES_TO_RAD*pOutRec->ra);
      if (kd_insert3(ptree,x,y,z,pOutRec) != 0) {
        printf("ERROR: fatal return from kd_insert3\n");
        exit(-1);
      }
    }
  }
  for (iteration = 0; iteration < 4; iteration++) {
    switch (iteration) {
    case 0:
      cursource = SOURCE_NGC;
    case 1:
      cursource = SOURCE_PGC;
      break;
    case 2:
      cursource = SOURCE_SGC;
      break;
    case 3:
      cursource = SOURCE_UGC;
      break;
    default:
      printf("ERROR: invalid iteration %d\n",iteration);
      exit(-1);
    }

    for (outputIndex = 0; outputIndex < totalCount; outputIndex++) {
      pOutRec = &output_table[outputIndex];
      if ((pOutRec->source == SOURCE_VSX) ||
          (pOutRec->source == SOURCE_OGLE)) {
        /* Do not consider VSX or OGLE sources */
        pOutRec->checked = 1;
        continue;
      }
      if (pOutRec->source == cursource) {
        z = sin(DEGREES_TO_RAD*pOutRec->dec);
        cosx = cos(DEGREES_TO_RAD*pOutRec->dec);
        x = cosx*cos(DEGREES_TO_RAD*pOutRec->ra);
        y = cosx*sin(DEGREES_TO_RAD*pOutRec->ra);
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
        presults = kd_nearest_range( ptree, pt, GALAXY_SEARCH_RADIUS*DEGREES_TO_RAD );
        resultSize = kd_res_size(presults);
        if (resultSize > maxResultSize) {
          maxResultSize = resultSize;
        }
        filterSize = 0;
        while( !kd_res_end( presults ) ) {
          /* get the data and position of the current result item */
          pChkRec = (POUTREC)kd_res_item( presults, pos );
          if ((pChkRec == pOutRec) ||
              (pChkRec->checked != 0) ||
              (pChkRec->duplicate != 0) ||
              (pChkRec->source == SOURCE_VSX) ||
              (pChkRec->source == SOURCE_OGLE) ||
              (pChkRec->source == cursource)) {
            rejectCount++;
            /* go to the next entry */
            kd_res_next( presults );
            continue;
          }
          distance = 3600.0*RAD_TO_DEGREES*sqrt(sqr(pt[0]-pos[0])+sqr(pt[1]-pos[1])+sqr(pt[2]-pos[2]));
          pFilterRec = &filter_table[filterSize];
          memcpy(pFilterRec,pChkRec,sizeof(OUTREC));
          pFilterRec->pOutRec = pChkRec;
          pFilterRec->distance = distance;
          filterSize++;
          if (filterSize >= MAX_FILTER_SIZE) {
            printf("ERROR: MAX_FILTER_SIZE %d is too small\n",MAX_FILTER_SIZE);
            exit(-1);
          }
#if 0
          printf("%12.6f arcsec filter %3d %31s %31s is near %31s %31s\n",distance,filterSize,pChkRec->catalogname,pChkRec->type,pOutRec->catalogname,pOutRec->type);
#endif
#if 0
          printf("dx %9.6f y %9.6f dz %9.6f distance %9.6f %s %s\n",pt[0]-pos[0],pt[1]-pos[1],pt[2]-pos[2],distance,pChkRec->catalogname,pChkRec->type);
#endif
          /* go to the next entry */
          kd_res_next( presults );
        }

        CheckDuplicates(pOutRec,filter_table,filterSize,cursource,&radiusCount,&catalogCount,&duplicateCount,0);
        kd_res_free(presults);
        if (filterSize > maxFilterSize) {
          maxFilterSize = filterSize;
        }
        pOutRec->checked = 1;

      }
    }
    time(&curTime);
    curTime -= startTime;
    printf("Done iteration %d checking for duplicates rejectCount %d radiusReject %d catalogReject %d duplicates %d maxResultSize %d maxFilterSize %d at %d seconds\n",
           iteration,
           rejectCount,
           radiusCount,
           catalogCount,
           duplicateCount,
           maxResultSize,
           maxFilterSize,
           curTime);
  }
#ifdef DEBUG_ZWICKI
  iteration = 4;
  cursource = SOURCE_UZC;
  for (outputIndex = 0; outputIndex < totalCount; outputIndex++) {
    pOutRec = &output_table[outputIndex];
    if ((pOutRec->source == SOURCE_VSX) ||
        (pOutRec->source == SOURCE_OGLE)) {
      /* Do not consider VSX sources */
      continue;
    }
    if ((pOutRec->source == cursource) &&
        (pOutRec->checked == 0)) {
      z = sin(DEGREES_TO_RAD*pOutRec->dec);
      cosx = cos(DEGREES_TO_RAD*pOutRec->dec);
      x = cosx*cos(DEGREES_TO_RAD*pOutRec->ra);
      y = cosx*sin(DEGREES_TO_RAD*pOutRec->ra);
      pt[0] = x;
      pt[1] = y;
      pt[2] = z;
      presults = kd_nearest_range( ptree, pt, GALAXY_SEARCH_RADIUS*DEGREES_TO_RAD );
      resultSize = kd_res_size(presults);
      if (resultSize > maxResultSize) {
        maxResultSize = resultSize;
      }
      filterSize = 0;
      while( !kd_res_end( presults ) ) {
        /* get the data and position of the current result item */
        pChkRec = (POUTREC)kd_res_item( presults, pos );
        if ((pChkRec == pOutRec) ||
            (pChkRec->duplicate != 0) ||
            (pChkRec->source == SOURCE_VSX) ||
            (pChkRec->source == SOURCE_OGLE) ||
            (pChkRec->source == cursource)) {
          /* go to the next entry */
          kd_res_next( presults );
          continue;
        }
        distance = 3600.0*RAD_TO_DEGREES*sqrt(sqr(pt[0]-pos[0])+sqr(pt[1]-pos[1])+sqr(pt[2]-pos[2]));
        pFilterRec = &filter_table[filterSize];
        memcpy(pFilterRec,pChkRec,sizeof(OUTREC));
        pFilterRec->pOutRec = pChkRec;
        pFilterRec->distance = distance;
        filterSize++;
        if (filterSize >= MAX_FILTER_SIZE) {
          printf("ERROR: MAX_FILTER_SIZE %d is too small\n",MAX_FILTER_SIZE);
          exit(-1);
        }
#if 0
        printf("%12.6f arcsec filter %3d %31s %31s is near %31s %31s\n",distance,filterSize,pChkRec->catalogname,pChkRec->type,pOutRec->catalogname,pOutRec->type);
#endif
#if 0
        printf("dx %9.6f y %9.6f dz %9.6f distance %9.6f %s %s\n",pt[0]-pos[0],pt[1]-pos[1],pt[2]-pos[2],distance,pChkRec->catalogname,pChkRec->type);
#endif
        /* go to the next entry */
        kd_res_next( presults );
      }

      CheckZwicki(pOutRec,filter_table,filterSize,cursource,&radiusCount,&catalogCount,&duplicateCount);
      kd_res_free(presults);
      if (filterSize > maxFilterSize) {
        maxFilterSize = filterSize;
      }
    }
  }
  time(&curTime);
  curTime -= startTime;
  printf("Done iteration %d checking for duplicates rejectCount %d radiusReject %d catalogReject %d duplicates %d maxResultSize %d maxFilterSize %d at %d seconds\n",
         iteration,
         rejectCount,
         radiusCount,
         catalogCount,
         duplicateCount,
         maxResultSize,
         maxFilterSize,
         curTime);
#endif /* DEBUG_ZWICKI */

	time(&curTime);
	curTime -= startTime;
  printf("Writing output table at %d seconds\n",curTime);

  if (galaxy_name[0] != 0) {
    galaxy_handle = Open(galaxy_name,"w");
    if (galaxy_handle == NULL) {
      errorFlag = 1;
      printf("Could not open binary file %s\n",galaxy_name);
    } else {
      galaxy_table = (PGALAXYREC)calloc(totalCount,sizeof(GALAXYREC));
      if (galaxy_table == NULL) {
        printf("ERROR: failed to allocate binary table of size %d\n",totalCount);
        exit(-1);
      }
    }

  }


  for (outputIndex = 0; outputIndex < totalCount; outputIndex++) {
    pOutRec = &output_table[outputIndex];
    if (pOutRec->checked == 0) {
      notCheckedCount++;
#if 0
      printf("Not checked %s %s\n",pOutRec->catalogname,pOutRec->type);
#endif
    }
    if (pOutRec->duplicate == 0) {
      if (strlen(pOutRec->catalogname) > maxCatalogName) {
        maxCatalogName = strlen(pOutRec->catalogname);
      }
      if (strlen(pOutRec->type) > maxType) {
        maxType = strlen(pOutRec->type);
      }
      
#ifdef DEBUG_TYPE
      if ((pOutRec->source == SOURCE_NGC) && (pOutRec->duplicate == 0)) {
        fprintf(output_handle,"%s\t%.5f\t%.5f\t%d\t%.2f\t%f\n",pOutRec->catalogname,pOutRec->ra,pOutRec->dec,pOutRec->source,pOutRec->catalogmag,pOutRec->radius);
      }
#else /* DEBUG_TYPE */
      fprintf(output_handle,"%s\t%.5f\t%.5f\t%s\t%.2f\t%f\n",pOutRec->catalogname,pOutRec->ra,pOutRec->dec,pOutRec->type,pOutRec->catalogmag,pOutRec->radius);
#endif /* DEBUG_TYPE */
      if (galaxy_table != NULL) {
        pGalaxyRec = &galaxy_table[finalCount];
        pGalaxyRec->galaxyflag       = GALAXY_FLAG;
        pGalaxyRec->galaxyversion    = GALAXY_VERSION;
        pGalaxyRec->ra         = pOutRec->ra;
        pGalaxyRec->dec        = pOutRec->dec;
        pGalaxyRec->catalogmag = pOutRec->catalogmag;
        pGalaxyRec->radius     = pOutRec->radius;
        strcpy(pGalaxyRec->catalogname,pOutRec->catalogname);
        strcpy(pGalaxyRec->galaxytype,pOutRec->type);
        if ((pOutRec->source == SOURCE_VSX) ||
            (pOutRec->source == SOURCE_OGLE)) {
          pGalaxyRec->variableFlag = 1;
          finalVariableCount++;
        } else {
          pGalaxyRec->variableFlag = 0;
          finalGalaxyCount++;
        }
      
      }
      finalCount++;
    }
  }
  if (galaxy_table != NULL) {
    writeItems = Write(galaxy_handle,galaxy_table,sizeof(GALAXYREC),finalCount);
    if (writeItems != finalCount) {
      printf("ERROR writing %s items %d\n",galaxy_name,writeItems);
      exit(-1);
    } else {
      printf("Wrote %d items of size %d to %s\n",writeItems,sizeof(GALAXYREC),galaxy_name);
    }

  }

  if (raw_handle != NULL) {
    fclose(raw_handle);
  }
  if (fold_handle != NULL) {
    fclose(fold_handle);
  }
	fclose(output_handle);
  if (galaxy_handle != NULL) {
    Close(galaxy_handle);
  }
	time(&curTime);
	curTime -= startTime;



 





	if (output_table != NULL) {
		free(output_table);
	}
	if (ugc_table != NULL) {
		free(ugc_table);
	}

	if (ogle_table != NULL) {
		free(ogle_table);
	}

	if (vector != NULL) {
		free(vector);
	}

	time(&curTime);
	curTime -= startTime;
  printf("maxCatalogName %d, maxType %d\n",maxCatalogName,maxType);
	printf("Execution Time: %d seconds; constant: %d duplicate: %d not checked %d total entries %d objects written: %d (%d variables %d galaxies) of %d \n",curTime,constantCount,duplicate2Count,notCheckedCount,totalCount,finalCount,finalVariableCount,finalGalaxyCount,MAX_INPUT_LINES);

	return(EXIT_SUCCESS);
}
