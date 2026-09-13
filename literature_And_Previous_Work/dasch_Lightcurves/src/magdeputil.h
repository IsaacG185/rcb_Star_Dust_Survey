// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* magdeputil.h
 *
 *  Routines to load and unload magnitude-dependent calibration tables
 *
 * Nov  8, 2011 Edward J. Los - Initial version
 */


#include <math.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"

int PhotLoadMagdepCorrections(MYSQL *pPhotConnection,int seriesId,int plateNumber,int solutionNumber,char *catalogString,PMAGDEPLIMITS pMagdepLimits,PMAGDEPCORRECTION *ppMagdepTable,int verbose);
int LoadMagdepCorrections(char *magdep_name,PMAGDEPLIMITS pMagdepLimits,PMAGDEPCORRECTION *pMagdepTable,int verbose);
