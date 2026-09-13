// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* computemag.c
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I /dasch/install/include  -L /dasch/install/lib -lm computemag.c pipelineutils.a  -ltable -lutil -lwcs -lmysqlclient -o computemag
 *
 *   computemag /dasch/Pipeline/ingest/master_all_good_scans.db -t -n 3 -c -b -r
 *  
 *  Given a master observation file, produce a calibration table with a magnitude estimate for
 *  each GSC2.3 ID.
 *
 *  Dec 17, 2007 Edward J. Los - Initial version
 *  Jan 19, 2008 Edward J. Los - Make the minimum good stars a parameter for debugging
 *  Feb 29, 2008 Edward J. Los - Add limiting_mag, spatial_bin, and magcal_local_error to the master file
 *  Mar  5, 2008 Edward J. Los - Add extinction
 *  Mar  8, 2008 Edward J. Los - Covert to starbase input
 *  Mar 25, 2008 Edward J. Los - Do not attempt to read in entire master database.
 *  Apr  1, 2008 Edward J. Los - Remove unmatched stars
 *  Apr 24, 2008 Edward J. Los - Change the difference between limiting_mag_local and magcal_local from 1.0 to 0.5
 *                               Change median_local to clip_median_local
 *                               Change rms_local to clip_rms_local
 *                               Change ngood to clip_ngood
 *  May 29, 2008 Edward J. Los - Remove high drad and Pickering wedge objects.
 *  Jun 10, 2008 Edward J. Los - Remove plate defects. cprrect drad amd wedge flags
 *  Jul  1, 2008 Edward J. Los - Add max_local2, the second brightest magnitude
 *                                   min_local2, the second dimmest magnitude
 *                                   range_local2 = max_local2-min_local2
 *                                   yrbegin, date of the earliest point
 *                                   yrend, date of the latest point
 *                                   Stdmag GSC magnitude
 *                               Reject all blended stars and unmatched blended stars when computing the magnitude
 *  Jul  8, 2008 Edward J. Los - Reject all stars with bad local bins (high zout).
 *  Jul 28, 2008 Edward J. Los - Reject all stars that are too bright.
 *  Aug 15, 2008 Edward J. Los - Add median_iso and range_iso, which are derived from the magcal_iso column.
 *                             - Restore blended stars to favor
 *  Sep 19, 2008 Edward J. Los - Return to rejecting blended stars
 *  Sep 25, 2008 Edward J. Los - Add a column showing the error_bar_factor defined in extract_lightcurves.c
 *  Dec 22, 2008 Edward J. Los - Check for low altitude (high extinction)
 *  Jan 28, 2009 Edward J. Los - Move the CAL_FLAGS (now in BFLAGS) to pipelineutils.h
 *                               Split FLAGS into AFLAGS and BFLAGS
 *  Feb 23, 2009 Edward J. Los - Implement common AFLAGS interpretation routine
 *  Mar  8, 2009 Edward J. Los - Add color
 *  Apr  4, 2009 Edward J. Los - Flag objects without colorterm correction in plots
 *  May 16, 2011 Edward J. Los - Change /media to /dasch
 */   


#include <math.h>
#include <time.h>
#include "table.h"
#include "pipelineutils.h"
#include "pipelineutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#define MAX_MASTER_NAME 256
#define MAX_BUFFER 1000
#define MAX_PLATE 2000


/* #define LOS_DEBUG 1 */  

typedef struct _starimage {
  char REF[MAX_REF];   /* GSC2.3.2 reference number */
  char Plate[MAX_REF];
  double magcal_local; /* Local magnitude calibration */
  double ra;          /* Right Ascension in degrees */
  double dec;         /* Declination in degrees */
  int spatial_bin;    /* Lowess magnitude spatial bin  number */
  int AFLAGS;          /* Flags keyword. See header of pipelineutils.h */
  int BFLAGS;          /* Flags keyword. See header of pipelineutils.h */
  double magcal_iso;  /* lowess magnitude estimate */
  double magcal_iso_rms; /* Overall error */
  double magcal_local_rms; /* Overall error */
  double limiting_mag_local; /* Limiting magnitude */
  double Stdmag;       /* GSC2.3.2 magnitude */
  double color;        /* GSC2.3.2 color */
  double Date;         /* Heliocentric Julian Date  */      
#if 0
  int NUMBER;         /* Sextractor reference number */
  double extinction;   /* extinction correction */
#endif
} STARIMAGE,*PSTARIMAGE;

int ReadMaster(PSTARIMAGE pMaster,
               File master_handle,
               TableHead master_header,
               TblDescriptor master_descriptor,
               TableRow* master_row)
{
  *master_row = table_rowget(master_handle,master_header,*master_row,NULL,NULL,0);
  if (*master_row == NULL) {
    return(0);
  }
  if (!table_loadrow(master_handle,master_header,*master_row,master_descriptor,(char *)pMaster)) {
    printf("ERROR: Read Master table_loadrow failed\n");
    return(0);
  }
  return(1);
}

/* Sorting routine borrowed from statstable.c by John Roll */

int main(int argc,char *argv[])
{
  char *argstr;
  char outputname[MAX_MASTER_NAME];
  FILE * outputHandle;
  int errorFlag = 0;
  char *inBuffer;
  char inLine[MAX_BUFFER];
  char copyLine[MAX_BUFFER];
  double vector1[MAX_PLATE]; /* Vector for locally corrected magnitudes (magcal_local) */
  double vector2[MAX_PLATE]; /* Vector for julian days */
  double vector3[MAX_PLATE]; /* Vector for lowess magnitudes (magcal_iso) */
  double vector4[MAX_PLATE]; /* Vector for magcal_local_rms */
  int nvals;
  int lineLen;
  time_t startTime;
  time_t curTime;
  int linecounter = 0;
  int tabCount;
  int index;
  char *charPtr;
  char cmdchar;
  int verbose = 0;
  char *rootPtr; 
  int rootOffset = 0;
  char *slashPtr = NULL;
  int compareResult;
  int outStars = 0;
  int npoints = 0;
  double min_local;
  double max_local;
  double min_local2;
  double max_local2;
  double min_date;
  double max_date;
  double min_date2;
  double max_date2;
  int clip_ngood;
  int full_ngood = 0;
  int max_full_ngood = 0;
  char max_REF[MAX_REF] = "";
  int minGoodStars = MIN_GOODSTARS;
  double full_med;
  double full_rms;
  double clip_med;
  double clip_rms;
  double median_iso;
  double median_iso_rms;
  double max_iso;
  double min_iso;
  double max_iso2;
  double min_iso2;
  int bin9count = 0;
  int pipelineCount = 0;
  int isoRmsCount = 0;
  int localRmsCount = 0;
  int limitingMagCount = 0;
  int doClip = 0;
  int filterBin = 0;
  int filterIsoRms = 0;
  int filterLocalRms = 0;
  int filterLimitingMag = 0;
  int numStars = 0;
  int highDradCount = 0;
  int wedgeCount = 0;
  int defectCount = 0;
  int hiZoutCount = 0;
  int tooBrightCount = 0;
  int lowAltitudeCount = 0;
  int blendCount = 0;

  File master_handle = NULL;
  char master_name[MAX_BUFFER];
  TableHead master_header = NULL;
  int master_nrecs = 0;
  int master_index;
  STARIMAGE old_master_record;
  PSTARIMAGE pOldMaster = &old_master_record;
  STARIMAGE master_record;
  PSTARIMAGE pMaster = &master_record;
  TblDescriptor master_descriptor = NULL;
  TableRow master_row = NULL;

  double error_bar_factor;
  double rawerrmed;
  double rawerrrms;
  double rawmed;
  double rawrms;


  double magcal_local;          /* blue magnitude */
  int result;
  int AFLAGSMASK = 0;
  int BFLAGSMASK = -1;
  double symbol;
  int rejectReason1;
  int rejectReason2;

  /* Loop through the arguments */
#ifdef LOS_DEBUG
#if 0
  strcpy(master_name,"/dasch/mosaic/ExposureData/Mosaics/rb/00174_00/match_rb00174_00_01r180ww_tnx.db");
#else
  strcpy(master_name,"/dasch/mosaic/ExposureData/Mosaics/rb/00174_00/match_duplicates.db");
#endif
  strcpy(outputname,"/dasch/mosaic/ExposureData/Mosaics/rb/00174_00/match_duplicates_u.db");
#else /* LOS_DEBUG */
  master_name[0] = 0;
  outputname[0] = 0;
#endif /* LOS_DEBUG */

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      if (strlen(master_name) > 0) {
        errorFlag = 1;
        printf("ERROR: list %s is being overwritten by %s\n",master_name,argstr);
      } else {
        strncpy(master_name,argstr,MAX_MASTER_NAME-1);
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'v':
        case 'V':
          verbose = 1;
          break;
        case 'c':
        case 'C':
          doClip = 1;
          break;
        case 'b':
        case 'B':
          filterBin = 1;
          break;
        case 'i':
        case 'I':
          filterIsoRms = 1;
          break;
        case 'r':
        case 'R':
          filterLocalRms = 1;
          break;
        case 't':
        case 'T':
          filterLimitingMag = 1;
          break;
        case 'n': /* Minimum number of good stars */
        case 'N':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&minGoodStars);
            if (nvals != 1) {
              fprintf(stderr,"ERROR: Unable to decode minGoodStars %s\n",*argv);
              errorFlag = 1;
            } else {
              if (minGoodStars < 1) {
                fprintf(stderr,"ERROR: illegal value for minGoodStars: %d  Using default %d\n",minGoodStars,MIN_GOODSTARS);
                minGoodStars = MIN_GOODSTARS;
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

  /* Attempt to open the list */
  if (master_name[0] == 0) {
    printf("Input file name not specified\n");
    errorFlag = 1;
  } else {

    master_handle = Open(master_name,"r");
    if (master_handle == NULL) {
      printf("Could not open file %s\n",master_name);
      errorFlag = 1;
    } 
  }
  /* Find the root */
  charPtr = master_name;
  while ((charPtr = strstr(charPtr,"/")) != NULL) {
    slashPtr = charPtr;
    charPtr++;
  }
  if (slashPtr != NULL) {
    rootOffset = slashPtr - master_name + 1;
  }

  rootPtr = strstr(&master_name[rootOffset],"master_");
  if ((rootPtr == NULL) || (rootPtr != &master_name[rootOffset])) {
    errorFlag = 1;
    printf("Error: expecting master_name %s to begin with 'master_'\n",master_name);
  } else {
    strcpy(outputname,master_name);
    outputname[rootOffset] = 0;
    strcat(outputname,"id_");
    strcat(outputname,&master_name[rootOffset+7]);
  }

  if (errorFlag) {
    printf("Usage: computemag <master_name> [-v] [-n <min good>\n");
    printf("       where -v is the verbose flag\n");
    printf("             -n is the minimum number of good stars (default 10)\n");
    printf("             -c perform three sigma clipping\n");
    printf("             -b filter out bin 9\n");
    printf("             -i filter if iso_rms > 1.0\n");
    printf("             -r filter if total rms > 1.0 \n");
    printf("             -t filter limiting_mag_local - magcal_local < 0.5\n");
    return(-1);
  }

  printf("computemag of %s %s \n Input Filename %s\n Output Filename %s\n",
         __DATE__,__TIME__,master_name,outputname);

 

  time(&startTime);
  
  outputHandle = fopen(outputname,"wt");
  if (outputHandle == NULL) {
    printf("Could not open file %s\n",outputname);
    exit(-1);
  } 
  /* Write the header of the output file */
  fprintf(outputHandle,"REF\tra\tdec\tStdmag\tcolor\tyrbegin\tyrend\tnpoints\tmin_local\tmax_local\trange_local\tmin_local2\tmax_local2\trange_local2\tmedian_local\trms_local\tngood\tclip_median_local\tclip_rms_local\tclip_ngood\tmedian_iso\trange_iso\terror_bar_factor\n");
  fprintf(outputHandle,"---\t--\t---\t------\t-----\t-------\t-----\t-------\t---------\t---------\t-----------\t----------\t----------\t------------\t------------\t---------\t-----\t-----------------\t--------------\t----------\t----------\t---------\t----------------\n");

  /* Now read in the master file */

  master_header = table_header(master_handle,TABLE_PARSE);
  if (master_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",master_name);
    return(-1);
  }

  master_descriptor = table_create_descrip(&master_nrecs,
                                           TblBuf,"REF"    ,TblOff(PSTARIMAGE,REF),MAX_REF,
                                           TblBuf,"Plate"    ,TblOff(PSTARIMAGE,Plate),MAX_REF,
                                           TblDbl,"magcal_local"  ,TblOff(PSTARIMAGE,magcal_local),
                                           TblDbl,"ra"     ,TblOff(PSTARIMAGE,ra),
                                           TblDbl,"dec"    ,TblOff(PSTARIMAGE,dec),
                                           TblDbl,"Stdmag"    ,TblOff(PSTARIMAGE,Stdmag),
                                           TblDbl,"color"    ,TblOff(PSTARIMAGE,color),
                                           TblDbl,"Date"    ,TblOff(PSTARIMAGE,Date),
                                           TblDbl,"magcal_local_rms"    ,TblOff(PSTARIMAGE,magcal_local_rms),
                                           TblDbl,"magcal_iso"    ,TblOff(PSTARIMAGE,magcal_iso),
                                           TblDbl,"magcal_iso_rms"    ,TblOff(PSTARIMAGE,magcal_iso_rms),
                                           TblDbl,"limiting_mag_local"    ,TblOff(PSTARIMAGE,limiting_mag_local),
                                           TblInt,"spatial_bin",TblOff(PSTARIMAGE,spatial_bin),
                                           TblInt,"AFLAGS",TblOff(PSTARIMAGE,AFLAGS),
                                           TblInt,"BFLAGS",TblOff(PSTARIMAGE,BFLAGS),
                                           0,"end",0);
  if (master_descriptor == NULL) {
    fprintf(stderr,"ERROR: Failed to allocate descriptor for %s\n",master_name);
    return(-1);
  }
  table_loadmap(master_header,master_descriptor);

  memset(pOldMaster,0,sizeof(STARIMAGE));
  memset(pMaster,0,sizeof(STARIMAGE));

  if (filterBin == 0) {
    AFLAGSMASK = (1 << FILTER_AFLAG_BIN9);
  }
  if (filterIsoRms == 0) {
    AFLAGSMASK = (1 << FILTER_AFLAG_ISO_RMS);
  }
  if (filterLocalRms == 0) {
    AFLAGSMASK = (1 << FILTER_AFLAG_LOCAL_RMS);
  }
  if (filterLimitingMag == 0) {
    AFLAGSMASK = (1 << FILTER_AFLAG_LIMITING_MAG);
  }



  while (1) {
    if(ReadMaster(pMaster,master_handle,master_header,master_descriptor,&master_row) == 0) {
      break;
    }
    master_nrecs++;

    if (strstr(pMaster->REF,"NONE") == NULL) {
      if ((pMaster->BFLAGS & (PIPELINE_FLAG << CAL_FLAG_SHIFT)) > 0) {
        pipelineCount++;
      }
      if (pMaster->magcal_iso_rms > MAX_ISO_RMS) {
        isoRmsCount++;
      }
      if (pMaster->magcal_local_rms > MAX_LOCAL_RMS) {
        localRmsCount++;
      }
      if (pMaster->limiting_mag_local - pMaster->magcal_local < MAX_LIMITING_MAG) {
        limitingMagCount++;
      }
      if (pMaster->spatial_bin == 9) {
        bin9count++;
      }
    }
    if (master_nrecs <= 1) {
      compareResult = 0;
    } else {
      compareResult = strcmp(pMaster->REF,pOldMaster->REF);
    }
    switch (compareResult) {

    case 1:
#if 0
      printf("At ref %s\n",pOldMaster->REF);
      if (strcmp(pOldMaster->REF,"S2202000100") == 0) {
        printf("At ref %s\n",pOldMaster->REF);
      }
#endif
      /* Here we have a new GSC2.3 ID, write out the old results */
      if (full_ngood > max_full_ngood) {
        max_full_ngood = full_ngood;
        strcpy(max_REF,pOldMaster->REF);
      }
      numStars++;
      if (CalcMedianAndRMS(full_ngood,minGoodStars,vector1,&full_med,&full_rms,0,3.0,0)) {
        ComputeRange(full_ngood,vector1,&min_local,&min_local2,&max_local2,&max_local);
        ComputeRange(full_ngood,vector2,&min_date,&min_date2,&max_date2,&max_date);

        CalcMedianAndRMS(full_ngood,minGoodStars,vector3,&median_iso,&median_iso_rms,0,3.0,0);
        ComputeRange(full_ngood,vector3,&min_iso,&min_iso2,&max_iso2,&max_iso);
      
        /* Find the zero-based clipped rms of our error bars */
        rawmed = full_med;
        rawrms = full_rms;
        if (CalcMedianAndRMS(full_ngood,minGoodStars,vector4,&rawerrmed,&rawerrrms,0,3.0,1) == 0) {
          rawerrmed = 0.0;
          rawerrrms = 99.0;
        }
        if ((rawrms == 0.0) || 
            (rawrms == 99.0) ||
            (rawerrrms == 0.0) || 
            (rawerrrms == 99.0)) {
          error_bar_factor = 1.0;
        } else {
          error_bar_factor = rawrms/rawerrrms;
        }
        if (error_bar_factor > 1.0) {
          error_bar_factor = 1.0;
        }



        if ((clip_ngood = CalcMedianAndRMS(full_ngood,minGoodStars,vector1,&clip_med,&clip_rms,doClip,3.0,0)) > 0) {
          outStars++;
          fprintf(outputHandle,"%s\t%f\t%f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.7f\t%d\t%.4f\t%.7f\t%d\t%.4f\t%.4f\t%.4f\n",
                  pOldMaster->REF,
                  pOldMaster->ra,
                  pOldMaster->dec,
                  pOldMaster->Stdmag,
                  pOldMaster->color,
                  jd2ep(min_date),
                  jd2ep(max_date),
                  npoints,
                  min_local,
                  max_local,
                  max_local-min_local,
                  min_local2,
                  max_local2,
                  max_local2-min_local2,
                  full_med,
                  full_rms,
                  full_ngood,
                  clip_med,
                  clip_rms,
                  clip_ngood,
                  median_iso,
                  max_iso-min_iso,
                  error_bar_factor);
        }
      }
      full_ngood = 0;
      npoints = 0;
 
      
      /* NOW FALL THROUGH TO TALLY THE LATEST LINE */

    case 0:
      if (strstr(pMaster->REF,"NONE") == NULL) {
#if 0
        if (strcmp(pMaster->REF,"S2202000100") == 0) {
          printf("Point %d for ref %s\n",npoints,pMaster->REF);
        }
#endif
#if 0
        if ((strcmp(pMaster->REF,"N120013339") == 0) &&
            (strcmp(pMaster->Plate,"rh12974") == 0)) {
          printf("At REF %s, plate %s\n",pMaster->REF,pMaster->Plate);
        }
#endif
        npoints++;


        result = DecodeAFLAGS(pMaster->AFLAGS,AFLAGSMASK,pMaster->BFLAGS,BFLAGSMASK,&rejectReason1,&rejectReason2,&symbol);
        if (result == 0) {
          if ((rejectReason1 & (1<<REJECT_REASON_DRAD)) != 0) {
            highDradCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_WEDGE)) != 0) {
            wedgeCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_DEFECT)) != 0) {
            defectCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_HIZOUT)) != 0) {
            hiZoutCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_TOO_BRIGHT)) != 0) {
            tooBrightCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_LOW_ALTITUDE)) != 0) {
            lowAltitudeCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_UNCERTAIN_DATE)) != 0) {
            lowAltitudeCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_BLEND)) != 0) {
            blendCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_MULTIPLE_BLEND)) != 0) {
            blendCount++;
          } else if ((rejectReason1 & (1<<REJECT_REASON_BLEND_NOMATCH)) != 0) {
            blendCount++;
          }
          
          break;
        }
       

        if (full_ngood >= MAX_PLATE) {
          fprintf(stderr,"ERROR: MAX_PLATE exceeded for %s\n",pMaster->REF);
          exit(-1);
        }
        vector1[full_ngood] = pMaster->magcal_local;
        vector2[full_ngood] = pMaster->Date;
        vector3[full_ngood] = pMaster->magcal_iso;
        vector4[full_ngood] = pMaster->magcal_local_rms;
#if 0
        printf("DEBUG: point %4d AFLAGS %08x plate %25s target %s\n",full_ngood,pMaster->AFLAGS,pMaster->Plate,pMaster->REF);
#endif
        full_ngood++;
      }
 
      break;
    default:
    
      fprintf(stderr,"ERROR: CATALOG NOT SORTED BY REF compare %d for REF: %s, oldREF %s\n",compareResult,pMaster->REF,pOldMaster->REF);
      exit(-1);
    }
    memcpy(pOldMaster,pMaster,sizeof(STARIMAGE));
    
  }
  if (full_ngood > max_full_ngood) {
    max_full_ngood = full_ngood;
    strcpy(max_REF,pOldMaster->REF);
  }
  numStars++;

  if (CalcMedianAndRMS(full_ngood,minGoodStars,vector1,&full_med,&full_rms,0,3.0,0)) {
    ComputeRange(full_ngood,vector1,&min_local,&min_local2,&max_local2,&max_local);
    ComputeRange(full_ngood,vector2,&min_date,&min_date2,&max_date2,&max_date);

    CalcMedianAndRMS(full_ngood,minGoodStars,vector3,&median_iso,&median_iso_rms,0,3.0,0);
    ComputeRange(full_ngood,vector3,&min_iso,&min_iso2,&max_iso2,&max_iso);

    /* Find the zero-based clipped rms of our error bars */
    rawmed = full_med;
    rawrms = full_rms;
    if (CalcMedianAndRMS(full_ngood,minGoodStars,vector4,&rawerrmed,&rawerrrms,0,3.0,1) == 0) {
      rawerrmed = 0.0;
      rawerrrms = 99.0;
    }
    if ((rawrms == 0.0) || 
        (rawrms == 99.0) ||
        (rawerrrms == 0.0) || 
        (rawerrrms == 99.0)) {
      error_bar_factor = 1.0;
    } else {
      error_bar_factor = rawrms/rawerrrms;
    }
    if (error_bar_factor > 1.0) {
      error_bar_factor = 1.0;
    }


 
    if ((clip_ngood = CalcMedianAndRMS(full_ngood,minGoodStars,vector1,&clip_med,&clip_rms,doClip,3.0,0)) > 0) {
      outStars++;
      fprintf(outputHandle,"%s\t%f\t%f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.7f\t%d\t%.4f\t%.7f\t%d\t%.4f\t%.4f\t%4f\n",
              pOldMaster->REF,
              pOldMaster->ra,
              pOldMaster->dec,
              pOldMaster->Stdmag,
              pOldMaster->color,
              jd2ep(min_date),
              jd2ep(max_date),
              npoints,
              min_local,
              max_local,
              max_local-min_local,
              min_local2,
              max_local2,
              max_local2-min_local2,
              full_med,
              full_rms,
              full_ngood,
              clip_med,
              clip_rms,
              clip_ngood,
              median_iso,
              max_iso-min_iso,
              error_bar_factor);
    }
  }


  if (master_header != NULL) {
    table_hdrfree(master_header);
  }
  if (master_handle != NULL) {
    Close(master_handle);
  }
  if (master_descriptor != NULL) {
    Free(master_descriptor);
  }
  if (master_row != NULL) {
    table_rowfree(master_row);
  }

  fclose(outputHandle);
 
  time(&curTime);
  curTime -= startTime;
  printf("Execution Time: %d seconds for Observations: %d Stars: %d Output: %d stars max points %d for %s bin9 %d isoRMS %d localRMS %d limitMag %d pipeline %d bin9 clip %d bin %d iso %d total %d limit %d high drad %d wedge %d defect %d hiZout %d tooBright %d too low %d blend %d\n",curTime,master_nrecs,numStars,outStars,max_full_ngood,max_REF,bin9count,isoRmsCount,localRmsCount,limitingMagCount,pipelineCount,doClip,filterBin,filterIsoRms,filterLocalRms,filterLimitingMag,highDradCount,wedgeCount,defectCount,hiZoutCount,tooBrightCount,lowAltitudeCount,blendCount);

  return(EXIT_SUCCESS);
}
