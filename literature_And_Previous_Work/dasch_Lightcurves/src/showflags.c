// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pipelineutils.h"

#define MAX_BUFFER 100

extern FITWCSBIT FitWCSMasks[];
extern QUALITYBIT qualityMasks[];
extern QUALITYBIT peakEvaluationMasks[];
extern int peakEvaluationTableSize;
extern FLAGSENTRY AflagsTable[];
extern int AflagsTableSize;
extern FLAGSENTRY BflagsTable[];
extern int BflagsTableSize;

typedef struct _classentry {
  int class;
  char *classname;
} CLASSENTRY, *PCLASSENTRY;

CLASSENTRY classTable[] = {
  { 0, "Star" },
  { 1, "Galaxy" },
  { 2, "Blend" },
  { 3, "Non-Star" },
  { 4, "Unclassified" },
  { 5, "Defect" },
  { 6, "Undefined" },
  { 7, "Duplicate Star" }
};

int classTableSize = sizeof(classTable) / sizeof(CLASSENTRY);

#ifdef NEW_COLOR_CALIBRATION

typedef struct _colorentry {
  int MAGFlag;
  char *colorname;
} COLORENTRY, *PCOLORENTRY;

COLORENTRY colorTable[] = {
  { 0, "FPGMAG and JPGMAG" },
  { 2, "DEPRECATED (old BMAG_FOR_JPGMAG)" },
  { 3, "COLOR_BMAG_FOR_JPGMAG" },
  { 4, "COLOR_VMAG_FOR_FPGMAG" },
  { 5, "COLOR_BMAG_VMAG_FOR_JPGMAG_FPGMAG" },
  { 6, "COLOR_TYCHO" },
  { 7, "COLOR_SKYMAP" },
  { 8, "DEPRECATED (old VMSG_FOR_FPGMAG)" },
  { 10, "DEPRECATED (old BMAG_FOR_JPGMAG|VMAG_FOR_FPGMAG) " },
  { 11, "Undefined" }
};

int colorTableSize = sizeof(colorTable) / sizeof(COLORENTRY);

#endif /* NEW_COLOR_CALIBRATION */


int
main(int argc, char *argv[])
{
  int errorFlag = 0;
  int nvals;
  int FLAGS;
  int tmpMAGFlag;
  int index;
  int class;
  PFLAGSENTRY pFlagEntry;
  PCLASSENTRY pClassEntry = NULL;
  PCOLORENTRY pColorEntry = NULL;
  char *argstr;
  char cmdchar;
  char numberString[MAX_BUFFER];
  int decodeRejectFlags = 0;
  int decodeBFLAGS = 0;
  int decodePeakEvaluationFlags = 0;
  int decodeQualityFlags = 0;
  int decodeFitWCSFLAGS = 0;
#ifdef NEW_COLOR_CALIBRATION
  int MAGFlag;
#endif /* NEW_COLOR_CALIBRATION */
  int webText = 0; /* Web friendly text */
  PFITWCSBIT pFitWCSBit = FitWCSMasks;
  PQUALITYBIT pQualityBit = qualityMasks;
  PQUALITYBIT pPeakEvaluationBit = peakEvaluationMasks;
  char *charPtr;
  char webname[128];

  numberString[0] = 0;

  /* Loop through the arguments */

  for (argv++; --argc > 0; argv++) {
    argstr = *argv;

    if (argstr[0] != '-') {
      /* This is the number to decode */
      if (strlen(numberString) > 0) {
        printf(
          "ERROR: number %s is being overwritten by %s\n",
          numberString,
          argstr
        );
        errorFlag = 1;
      } else {
        strcpy(numberString, argstr);
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'r': /* Decode reject flags */
        case 'R':
          decodeRejectFlags = 1;
          break;

        case 'b': /* Decode BFLAGS */
        case 'B':
          decodeBFLAGS = 1;
          break;

        case 'q': /* Decode quality flags */
        case 'Q':
          decodeQualityFlags = 1;
          break;

        case 't': /* Decode transient candidate peak evaluation flags */
        case 'T':
          decodePeakEvaluationFlags = 1;
          break;

        case 'f': /* Decode FitWCS */
        case 'F':
          decodeFitWCSFLAGS = 1;
          break;

        case 'w': /* Web Friendly Text */
        case 'W':
          webText = 1;
          break;

        default:
          fprintf(stderr,"ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;
        }
      }
    }
  }

  if (numberString[0] == 0) {
    printf("ERROR: Need at least one argument\n");
    errorFlag = 1;
  } else {
    nvals = sscanf(numberString, "%d", &FLAGS);

    if (nvals != 1) {
      printf("ERROR: can not decode argument\n");
      errorFlag = 1;
    }
  }

  if (errorFlag != 0) {
    printf("Usage: showflags [-r|-b|-f] <decimal FLAGS value>\n");
    printf("                     where -b decodes BFLAGS\n");
    printf("                     where -r decodes the reject_reason\n");
    printf("                     where -f decodes FitWCS flags\n");
    printf("                     where -q decodes QUALITY flags\n");
    printf("                     where -w provides web-friendly text\n");
    printf("                     where -t decodes rejection reasons for transient candiate flares\n");
    exit(-1);
  }

  if (
    decodeRejectFlags == 0 &&
    decodeBFLAGS == 0 &&
    decodeQualityFlags == 0 &&
    decodeFitWCSFLAGS == 0 &&
    decodePeakEvaluationFlags == 0
  ) {
    /* Here we decode the AFLAGS field */
    class = (FLAGS >> GSC_CLASS_BIT) & CLASS_MASK;

    for (index = 0; index < classTableSize; index++) {
      pClassEntry = & classTable[index];

      if (pClassEntry->class == class) {
        break;
      }
    }

    if (webText) {
      printf("AFLAGS: %10d(d) ", FLAGS);
      printf(" %s\n", pClassEntry->classname);
    } else {
      printf("AFLAGS value: %10d(d) %08x(x)", FLAGS, FLAGS);
      printf(" GSC class %d: %s\n", class, pClassEntry->classname);
    }

    for (index = 0; index < AflagsTableSize; index++) {
      pFlagEntry = &AflagsTable[index];
      pFlagEntry->flagmask = 1 << pFlagEntry->flagbit;
    }

    for (index = 0; index < AflagsTableSize; index++) {
      pFlagEntry = &AflagsTable[index];

      if (strlen(pFlagEntry->webname) + 1 > sizeof(webname)) {
        printf("ERROR: size of webname exceeded\n");
      }

      strcpy(webname, pFlagEntry->webname);

      if ((charPtr = strchr(webname, '(')) != NULL) {
        *charPtr = 0;
      }

      if ((FLAGS & pFlagEntry->flagmask) != 0) {
        if (webText) {
          printf(
            "Bit: %2d  %s\n",
            pFlagEntry->flagbit,
            webname
          );
        } else {
          printf(
            "Bit: %2d,  Mask: %08x, Name %s\n",
            pFlagEntry->flagbit,
            pFlagEntry->flagmask,
            pFlagEntry->flagname
          );
        }
      }
    }
  } else if (decodeBFLAGS) {
    if (webText) {
      printf("BFLAGS: %10d(d)", FLAGS);
    } else {
      printf("BFLAGS value: %10d(d) %08x(x)", FLAGS, FLAGS);
    }

#ifdef NEW_COLOR_CALIBRATION
    MAGFlag = (FLAGS >> GSC_MAGNITUDE_FLAG_BIT) & MAGNITUDE_MASK;

    for (index = 0; index < colorTableSize; index++) {
      pColorEntry = &colorTable[index];

      if (pColorEntry->MAGFlag == MAGFlag) {
        break;
      }
    }

    if ((FLAGS & (1 <<FILTER_BFLAG_KEPLER)) == 0) {
      if (webText) {
        printf(" Catalog: %s", pColorEntry->colorname);
      } else {
        printf(
          " MAGFlag: %2d(d) %02x(x) %s",
          MAGFlag,
          MAGFlag,
          pColorEntry->colorname
        );
      }
    }
#endif /* NEW_COLOR_CALIBRATION */

    /* Here we decode the BFLAGS field */
    for (index = 0; index < BflagsTableSize; index++) {
      pFlagEntry = &BflagsTable[index];
      pFlagEntry->flagmask = 1 << pFlagEntry->flagbit;
    }

    if ((FLAGS & (1 <<FILTER_BFLAG_KEPLER)) != 0) {
      tmpMAGFlag = FLAGS >> GSC_MAGNITUDE_FLAG_BIT;
      tmpMAGFlag &= MAGNITUDE_MASK;
      FLAGS &= ~(MAGNITUDE_MASK << GSC_MAGNITUDE_FLAG_BIT);

      if (tmpMAGFlag < 0 || tmpMAGFlag >= KEPLER_CQ_MAGFLAG_MAX) {
        printf(" Unrecognized Kepler Input Catalog source %d", tmpMAGFlag);
      } else {
        if (webText) {
          printf(" Catalog Source: %s", keplerSourceText[tmpMAGFlag]);
        } else {
          printf(" Kepler Input Catalog Source: %s", keplerSourceText[tmpMAGFlag]);
        }
      }
    }

    printf("\n");

    for (index = 0; index < BflagsTableSize; index++) {
      pFlagEntry = &BflagsTable[index];

      if ((FLAGS & pFlagEntry->flagmask) != 0) {
        if (webText) {
          printf(
            "Bit: %2d  %s\n",
            pFlagEntry->flagbit,
            pFlagEntry->webname
          );
        } else {
          printf(
            "Bit: %2d,  Mask: %08x, Name %s\n",
            pFlagEntry->flagbit,
            pFlagEntry->flagmask,
            pFlagEntry->flagname
          );
        }
      }
    }
  } else if (decodePeakEvaluationFlags) {
    if (webText) {
      printf("PEAK EVALUATION FLAGS: %10d(d)\n", FLAGS);
    } else {
      printf("PEAK EVALUATION FLAGS value: %10d(d) %08x(x)\n", FLAGS, FLAGS);
    }

    /* Here we decode the peakEvaluation FLAGS field */
    while (pPeakEvaluationBit->qualityMask != 0) {
      if ((pPeakEvaluationBit->qualityMask  & FLAGS) != 0) {
        printf(
          "Mask: %08x, Name %s\n",
          pPeakEvaluationBit->qualityMask,
          pPeakEvaluationBit->qualityDescr
        );
      }

      pPeakEvaluationBit++;
    }
  } else if (decodeFitWCSFLAGS) {
    if (webText) {
      printf("FitWCS FLAGS: %10d(d)\n", FLAGS);
    } else {
      printf("FitWCS FLAGS value: %10d(d) %08x(x)\n", FLAGS, FLAGS);
    }

    /* Here we decode the FitWCS FLAGS field */
    while (pFitWCSBit->FitWCSMask != 0) {
      if ((pFitWCSBit->FitWCSMask  & FLAGS) != 0) {
        printf(
          "Mask: %08x, Name %s\n",
          pFitWCSBit->FitWCSMask,
          pFitWCSBit->FitWCSDescr
        );
      }

      pFitWCSBit++;
    }
  } else if (decodeQualityFlags) {
    if (FLAGS == QUALITY_UNINITIALIZED) {
      printf("No valid quality flags available\n");
    } else {
      if (webText) {
        printf("QUALITY FLAGS: %10d(d)\n", FLAGS);
      } else {
        printf("QUALITY FLAGS value: %10d(d) %08x(x)\n", FLAGS, FLAGS);
      }

      /* Here we decode the FitWCS FLAGS field */
      while (pQualityBit->qualityMask != 0) {
        if ((pQualityBit->qualityMask & FLAGS) != 0) {
          printf(
            "Mask: %08x, Name %s\n",
            pQualityBit->qualityMask,
            pQualityBit->qualityDescr
          );
        }

        pQualityBit++;
      }
    }
  } else {
    /* Here we decode reject reasons */
    printf("reject_reason value: %10d(d) %08x(x)\n", FLAGS, FLAGS);

    for (index = 0; index < REJECT_REASON_MAX; index++) {
      if ((FLAGS & (1 << index)) != 0) {
        printf(
          "Bit: %2d,  Mask: %08x, Name %s\n",
          index,
          (1 << index),
          GetRejectReasonText(index)
        );
      }
    }
  }

  return 0;
}

/* Oct  7, 2008 Edward J. Los - Initial Version
 * Oct 20, 2008 Edward J. Los - Decode the GSC class field
 * Dec  8, 2008 Edward J. Los - Decode the reject reasons
 * Jan 28, 2009 Edward J. Los - split FLAGS into AFLAGS and BFLAGS
 * Feb  6, 2009 Edward J. Los - Convert FILTER_AFLAG_BLEND_COPY to FILTER_AFLAG_BLEND since these now have distinct meanings
 * Feb 17, 2009 Edward J. Los - Add FILTER_BFLAG_DRAD_ADJUST
 * Aug 20, 2009 Edward J. Los - Add kepler catalog support
 * Nov 16, 2009 Edward J. Los - Add FILTER_AFLAG_UNCERTAIN_DATE and FILTER_AFLAG_MULTIPLE_BLEND
 * Dec  7, 2009 Edward J. Los - Add FILTER_AFLAG_MULTIPLE_NONE
 * Dec 11, 2009 Edward J. Los - Add new color transformation support
 * Jan 22, 2010 Edward J. Los - Add plate quality support
 * Jun  7, 2010 Edward J. Los - Add web-friendly output
 * Oct 30, 2011 Edward J. Los - Add magnitude-dependent star correction support (FILTER_BFLAG_MAGDEP_MAGCOR)
 * Feb 25, 2013 Edward J. Los - Add FitWCS decoding
 * May 16, 2015 Edward J. Los - Truncate the bit number in webname.
 * Aug 30, 2015 Edward J. Los - Add support for the QUALITY bits
 * Feb 26, 2016 Edward J. Los - Do not display unitilized QUALITY bits field
 * Mar 28, 2017 Edward J. Los - Add reasons for rejection of transient candidate flares
 */
