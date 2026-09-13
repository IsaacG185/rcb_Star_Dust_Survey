// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

// Originally derived from libwcs's scanread.h; major adaptation for DASCH.

#ifndef __DASCH_SCANREAD_H__
#define __DASCH_SCANREAD_H__

#define SCANCATSTARS 800000 /* Approximate number of plates in the collection */
#define SCANCAT -5

#define MAX_QUERY_STRING 3000 /* Maximum MySQL query string */
#define MAX_SERIES_STRING  20 /* Maximum series name */
#define MAX_RIGHTASCENSION_STRING 13
#define MAX_DECLINATION_STRING 16
#define MAX_DATE_STRING 25
#define MAX_CLASS_STRING 8
#define MAX_EXPOSURE_STRING 8
#define MAX_CTYPE_STRING 10
#define MAX_PLATE_NUMBER 100000
#define MAX_EXPOSURE_NUMBER 1000
#define MAX_ROTATION 10

#define ORIENTATION_NORTH_FORWARD 1
#define ORIENTATION_NORTH_FORWARD_MIRRORED 2
#define ORIENTATION_NORTH_LEFT 3
#define ORIENTATION_NORTH_LEFT_MIRRORED 4

/* FitWCS bitmask definitions */
#define FITWCS_SELECTED    1
#define FITWCS_MIKESHAWFIT 2
#define FITWCS_COMPLETED   4
#define FITWCS_COPIED      8

#define PIXELS_PER_MM 90.9090


typedef struct seriesTable {
  int seriesId;
  char series[MAX_SERIES_STRING]; /* Series */
  double plateScale;
  int orientation;
  int sequestered;
} SERIESTABLE, *PSERIESTABLE;


typedef struct sortTable {
  double julianDate; /* This item must be first for sorting */
  char *listString; /* Pointer to the list string in tclass */
} SORTTABLE, *PSORTTABLE;


void setMosaicFlag(int flag);
void setPessimisticFlag(int flag);
void setSingleColumnFlag(int flag);
int getSingleColumnFlag();
int getMosaicFlag();
void setPlateSeries(char *strPtr);
void setPlateNumber(char *strPtr);
void setPlateClass(char *strPtr);
void setSearchType(char *strPtr);
void setWebTable(char *strPtr);
void setWebExtract(char *strPtr);
void setWebFileType(char *strPtr);
void setWebArchiveType(char *strPtr);
void setRegionFlag(char *strPtr);
void setTarballExtract(char *strPtr);
void ScanCatNum(struct StarCat *sc,int nnfld,int nndec,double dnum,char * numstr);
void setSearchParameters(double epoch1,double epoch2);

int scanread(
  char	*bincat,	/* Name of reference star catalog file */
  double	cra,		/* Search center J2000 right ascension in degrees */
  double	cdec,		/* Search center J2000 declination in degrees */
  double drad,           /* Search half-width and half-height in degrees */
  int	sysout,		/* Search coordinate system */
  double	eqout,		/* Search coordinate equinox */
  double	epout,		/* Proper motion epoch (0.0 for no proper motion) */
  int	nstarmax,	/* Maximum number of sources to be returned */
  struct StarCat **starcat, /* Star catalog data structure */
  double	*tnum,		/* Array of catalog numbers (returned) */
  double	*tra,		/* Array of right ascensions (returned) */
  double	*tdec,		/* Array of declinations (returned) */
  double	**tmag,		/* 2-D Array of magnitudes (returned) */
  char	**tclass,	 /* Array of spectral classes (returned) */
  int	nlog
);

void scanclose();

void dasch_scanread_set_query_tool_mode(int mode);

#endif
