// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* magdeputil.c
 *
 *  Routines to load and unload magnitude-dependent calibration tables
 *
 *  gcc -ggdb -c -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I/usr/include/mysql  -I /dasch/install/include  -L /dasch/install/lib  magdeputil.c -o magdeputil.o
 *
 * Nov  8, 2011 Edward J. Los - Initial version
 * Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
 * Nov 27, 2016 Edward J. Los - correct double free of magdep_table.
 * Jul  5, 2019 Edward J. Los - Improve error reporting for a table to fix double deallocation of magdep_table
 * Jul 15, 2019 Edward J. Los - Remove debugging statements and back out above change
 *
 */


#include <math.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "magdeputil.h"
int ValidateMagdepTable(PMAGDEPLIMITS pMagdepLimits,PMAGDEPCORRECTION magdep_table,size_t magdep_nrecs,char *magdep_name,PMAGDEPCORRECTION *pMagdepTable)
{
  PMAGDEPCORRECTION pMagdep;
  PMAGDEPCORRECTION pMagdep2;
  int magdep_index;
  int zeroFlag = 0;
  PMAGDEPCORRECTION magdep_table2 = NULL;
  int magdep_bin;


  if (*pMagdepTable != NULL) {
    printf("ERROR: LoadMagdepCorrections: pMagdepTable already allocated %x for %s\n",pMagdepTable,magdep_name);
    exit(-1);
  }

  for (magdep_index = 0; magdep_index < magdep_nrecs; magdep_index++) {
    pMagdep = &magdep_table[magdep_index];

    pMagdep->flag = 0;
    if (pMagdepLimits->kxb < pMagdep->ixb) {
      pMagdepLimits->kxb = pMagdep->ixb;
    }
    if (pMagdep->ixb == 0) {
      zeroFlag = 1;
    }
    if (pMagdepLimits->kyb < pMagdep->iyb) {
      pMagdepLimits->kyb = pMagdep->iyb;
    }
    if (pMagdep->iyb == 0) {
      zeroFlag = 1;
    }
    if (pMagdepLimits->kmagb < pMagdep->imagb) {
      pMagdepLimits->kmagb = pMagdep->imagb;
    }
    if (pMagdep->imagb == 0) {
      zeroFlag = 1;
    }
  }
  if (zeroFlag == 0) {
    printf("ERROR: ValidateMagdepTable zeroFlag is zero for %s\n",magdep_name);
    exit(-1);
  }
  pMagdepLimits->kxb++;
  pMagdepLimits->kyb++;
  pMagdepLimits->kmagb += 2; /* This array starts with magnitude 0 */
  if ((pMagdepLimits->kxb*pMagdepLimits->kyb*(pMagdepLimits->kmagb-1)) != magdep_nrecs) {
    printf("ERROR: ValidateMagdepTable kxb %d kyb %d kmagb %d not nrecs %d for %s\n",
           pMagdepLimits->kxb,
           pMagdepLimits->kyb,
           pMagdepLimits->kmagb,
           magdep_nrecs,
           magdep_name);
#if 0
    Free(magdep_table);
#endif
    magdep_table = NULL;
    memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
    return(-1);
  }

  /* Allocate linear arrays */
  pMagdepLimits->xcoord_magdep = (double *)calloc(pMagdepLimits->kxb,sizeof(double));
  pMagdepLimits->ycoord_magdep = (double *)calloc(pMagdepLimits->kyb,sizeof(double));
  pMagdepLimits->magdep_bin_edge  = (double *)calloc(pMagdepLimits->kmagb,sizeof(double));
  if ((pMagdepLimits->xcoord_magdep == NULL) ||
      (pMagdepLimits->ycoord_magdep == NULL) ||
      (pMagdepLimits->magdep_bin_edge  == NULL)) {
    printf("ERROR: failed to allocate xcoord_magdep %x, ycoord_magdep %x magdep_bin_edge %x for %s\n",
           pMagdepLimits->xcoord_magdep,
           pMagdepLimits->ycoord_magdep,
           pMagdepLimits->magdep_bin_edge,
           magdep_name);
#if 0
    Free(magdep_table);
#endif
    magdep_table = NULL;
    FreeMagdepSubarrays(pMagdepLimits);
    memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
    return(-1);
  }


  /* Now allocate our final magnitude-dependent correction table */
  magdep_table2 = (PMAGDEPCORRECTION)calloc(magdep_nrecs,sizeof(MAGDEPCORRECTION));
  if (magdep_table2 == NULL) {
    printf("ERROR: Failed to allocate magdep_table2 for %s\n",magdep_name);
#if 0
    Free(magdep_table);
#endif
    magdep_table = NULL;
    FreeMagdepSubarrays(pMagdepLimits);
    memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
    return(-1);
  }
  for (magdep_index = 0; magdep_index < magdep_nrecs; magdep_index++) {
    pMagdep = &magdep_table[magdep_index];

    magdep_bin = pMagdep->ixb + (pMagdep->iyb*pMagdepLimits->kxb) + (pMagdep->imagb * pMagdepLimits->kxb * pMagdepLimits->kyb);
    if ((magdep_bin < 0) || (magdep_bin >= magdep_nrecs)) {
      printf("ERROR: magdep_bin %d exceeds limits for kxb %d kyb %d kmagb %d for %s\n",
             magdep_bin,
             pMagdepLimits->kxb,
             pMagdepLimits->kyb,
             pMagdepLimits->kmagb,
             magdep_name);
#if 0
      Free(magdep_table);
#endif
      magdep_table = NULL;
      free(magdep_table2);
      FreeMagdepSubarrays(pMagdepLimits);
      memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
      return(-1);
    }
    pMagdep->magdep_bin = magdep_bin;
    pMagdep2 = &magdep_table2[magdep_bin];
    memcpy(pMagdep2,pMagdep,sizeof(MAGDEPCORRECTION));
    pMagdep2->flag++;
    if (pMagdep2->flag > 1) {
      printf("ERROR: magdep_bin %d appears twice for kxb %d kyb %d kmagb %d for %s\n",
             magdep_bin,
             pMagdepLimits->kxb,
             pMagdepLimits->kyb,
             pMagdepLimits->kmagb,
             magdep_name);
#if 0
      Free(magdep_table);
#endif
      magdep_table = NULL;
      free(magdep_table2);
      FreeMagdepSubarrays(pMagdepLimits);
      memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
      return(-1);
    }
    if (pMagdep->ixb >= pMagdepLimits->kxb) {
      printf("ERROR magdep_record %d ixb index %d exceeds %d for %s ",magdep_index,pMagdep->ixb,pMagdepLimits->kxb,magdep_name);
#if 0
      Free(magdep_table);
#endif
      magdep_table = NULL;
      free(magdep_table2);
      FreeMagdepSubarrays(pMagdepLimits);
      memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
      return(-1);
    }
    pMagdepLimits->xcoord_magdep[pMagdep->ixb] = pMagdep->xcoord_magdep;
    if (pMagdep->iyb >= pMagdepLimits->kyb) {
      printf("ERROR magdep_record %d iyb index %d exceeds %d for %s ",magdep_index,pMagdep->iyb,pMagdepLimits->kyb,magdep_name);
#if 0
      Free(magdep_table);
#endif
      magdep_table = NULL;
      free(magdep_table2);
      FreeMagdepSubarrays(pMagdepLimits);
      memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
      return(-1);
    }
    pMagdepLimits->ycoord_magdep[pMagdep->iyb] = pMagdep->ycoord_magdep;
    if ((pMagdep->imagb+1) >= pMagdepLimits->kmagb) {
      printf("ERROR magdep_record %d imagb index %d exceeds %d for %s ",magdep_index,pMagdep->imagb+1,pMagdepLimits->kmagb,magdep_name);
#if 0
      Free(magdep_table);
#endif
      magdep_table = NULL;
      free(magdep_table2);
      FreeMagdepSubarrays(pMagdepLimits);
      memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
      return(-1);
    }
    pMagdepLimits->magdep_bin_edge[pMagdep->imagb+1] = pMagdep->magdep_bin_edge;
  }
  if (pMagdepLimits->kxb > 1) {
    pMagdepLimits->dx = pMagdepLimits->xcoord_magdep[1] - pMagdepLimits->xcoord_magdep[0];
  }
  if (pMagdepLimits->kyb > 1) {
    pMagdepLimits->dy = pMagdepLimits->ycoord_magdep[1] - pMagdepLimits->ycoord_magdep[0];
  }



  for (magdep_index = 0; magdep_index < magdep_nrecs; magdep_index++) {
    pMagdep2 = &magdep_table2[magdep_index];
    if (pMagdep2->flag == 0) {
      printf("ERROR: magdep_index %d is not present for magdep_index %s\n",
             magdep_index,magdep_name);
      printf("ERROR: kxb %d kyb %d kmagb %d nrecs %d\n",
             pMagdepLimits->kxb,
             pMagdepLimits->kyb,
             pMagdepLimits->kmagb,
             pMagdepLimits->nrecs);
#if 0
      Free(magdep_table);
#endif
      magdep_table = NULL;
      free(magdep_table2);
      FreeMagdepSubarrays(pMagdepLimits);
      memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
      return(-1);

    }
  }
#if 0  /* Dump the linear arrays */
  index = 0;
  indexFlag = 1;
  printf("dx %f dy %f\n",pMagdepLimits->dx,pMagdepLimits->dy);
  while (indexFlag == 1) {
    indexFlag = 0;
    printf("index %d ",index);
    if (index < pMagdepLimits->kxb) {
      indexFlag = 1;
      printf("xcoord_magdep %f ",pMagdepLimits->xcoord_magdep[index]);
    }
    if (index < pMagdepLimits->kyb) {
      indexFlag = 1;
      printf("ycoord_magdep %f ",pMagdepLimits->ycoord_magdep[index]);
    }
    if (index < (pMagdepLimits->kmagb)) {
      indexFlag = 1;
      printf("magdep_bin_edge %f ",pMagdepLimits->magdep_bin_edge[index]);
    }
    printf("\n");
    if (indexFlag == 0) {
      break;
    }
    index++;
  }
#endif
  *pMagdepTable = magdep_table2;
  pMagdepLimits->nrecs = magdep_nrecs;

  return(0);
}

int LoadMagdepCorrections(char *magdep_name,PMAGDEPLIMITS pMagdepLimits,PMAGDEPCORRECTION *ppMagdepTable,int verbose)
{
  File magdep_handle = NULL;
  TableHead magdep_header = NULL;
  PMAGDEPCORRECTION magdep_table = NULL;
  size_t magdep_nrecs = 0;
  int returnStatus = 0;

  memset(pMagdepLimits,0,sizeof(MAGDEPLIMITS));
  /* Open the magnitude-dependent calibration file */
  magdep_handle = Open(magdep_name,"r");
  if (magdep_handle == NULL) {
    return -1;
  } else {
    if (verbose) {
      printf("Found magnitude-dependent calibration file %s\n",magdep_name);
    }
  }

  magdep_header = table_header(magdep_handle,TABLE_PARSE);
  if (magdep_header == NULL) {
    printf("ERROR: Failed to read header for %s\n",magdep_name);
    return(-1);
  }

  magdep_table = table_loadva(magdep_handle,
                             &magdep_header,
                             NULL, /* hbase */
                             NULL, /* rows */
                             NULL,
                             sizeof(MAGDEPCORRECTION),
                             &magdep_nrecs,
                             TblInt,"ixb" ,TblOff(PMAGDEPCORRECTION,ixb),
                             TblInt,"iyb" ,TblOff(PMAGDEPCORRECTION,iyb),
                             TblInt,"imagb" ,TblOff(PMAGDEPCORRECTION,imagb),
                             TblInt,"nstar_magdep" ,TblOff(PMAGDEPCORRECTION,nstar_magdep),
                             TblInt,"magdep_bin_size" ,TblOff(PMAGDEPCORRECTION,magdep_bin_size),
                             TblDbl,"xcoord_magdep",TblOff(PMAGDEPCORRECTION,xcoord_magdep),
                             TblDbl,"ycoord_magdep",TblOff(PMAGDEPCORRECTION,ycoord_magdep),
                             TblDbl,"magdep_bin_edge",TblOff(PMAGDEPCORRECTION,magdep_bin_edge),
                             TblDbl,"magdep_bin_median",TblOff(PMAGDEPCORRECTION,magdep_bin_median),
                             TblDbl,"magdep_bin_magcor",TblOff(PMAGDEPCORRECTION,magdep_bin_magcor),
                             TblDbl,"magcal_magdep_rms",TblOff(PMAGDEPCORRECTION,magcal_magdep_rms),
                             TblDbl,"magdep_bin_quality",TblOff(PMAGDEPCORRECTION,magdep_bin_quality),
                             0,"end",0);
  if (magdep_table == NULL) {
    printf("ERROR: Failed to read table for %s\n",magdep_name);
    return(-1);
  }
  if (verbose) {
    printf("read %d records for %s\n",magdep_nrecs,magdep_name);
  }
  returnStatus = ValidateMagdepTable(pMagdepLimits,magdep_table,magdep_nrecs,magdep_name,ppMagdepTable);

  if (magdep_handle != NULL) {
    Close(magdep_handle);
    magdep_handle = NULL;
  }
  if (magdep_header != NULL) {
    table_hdrfree(magdep_header);
    magdep_header = NULL;
  }
  if (magdep_table != NULL) {
    Free(magdep_table);
    magdep_table = NULL;
  }
  return(returnStatus);
}

void FreeMagdepSubarrays(PMAGDEPLIMITS pMagdepLimits)
{
  if (pMagdepLimits->xcoord_magdep != NULL) {
    free(pMagdepLimits->xcoord_magdep);
    pMagdepLimits->xcoord_magdep = NULL;
  }
  if (pMagdepLimits->ycoord_magdep != NULL) {
    free(pMagdepLimits->ycoord_magdep);
    pMagdepLimits->ycoord_magdep = NULL;
  }
  if (pMagdepLimits->magdep_bin_edge != NULL) {
    free(pMagdepLimits->magdep_bin_edge);
    pMagdepLimits->magdep_bin_edge = NULL;
  }
}

/* If rejectFlag is set, then set low quality bins to zero */
double GetMagdepBinMagcor(PMAGDEPLIMITS pMagdepLimits,PMAGDEPCORRECTION pMagdepTable,double X_IMAGE,double Y_IMAGE,double Stdmag,int *pmagdep_bin,double *pmagcal_magdep_rms,double *pmagdep_bin_magcor,int rejectFlag)
{
  int magdep_bin1;
  int magdep_bin2;
  int ix;
  int iy;
  int imag;
  PMAGDEPCORRECTION pMagdep1;
  PMAGDEPCORRECTION pMagdep2;
  double magcor_value;


  magdep_bin1 = GetMagdepBin(pMagdepLimits,X_IMAGE,Y_IMAGE,Stdmag,&ix,&iy,&imag);
  pMagdep1 = &pMagdepTable[magdep_bin1];
  *pmagdep_bin = magdep_bin1;
  *pmagcal_magdep_rms = 99.0;
  *pmagdep_bin_magcor = 0;

  if (rejectFlag) {
    if (fabs(pMagdep1->magdep_bin_quality) < MAGCAL_MAGDEP_QUALITY_LIMIT) {
      /* Bin has low quality, return no correction */
      return(0);
    }
  }
  *pmagcal_magdep_rms = pMagdep1->magcal_magdep_rms;
  if ((imag == 0) && (Stdmag < pMagdep1->magdep_bin_median)) {
    /* Brighter than the brightest median */
    *pmagdep_bin_magcor = pMagdep1->magdep_bin_magcor;
    return(pMagdep1->magdep_bin_magcor);
  } else if ((imag == (pMagdepLimits->kmagb-2)) && (Stdmag >= pMagdep1->magdep_bin_median)) {
    /* Dimmer than the dimmest median */
    *pmagdep_bin_magcor = pMagdep1->magdep_bin_magcor;
    return(pMagdep1->magdep_bin_magcor);
  } else {
    /* In between, must interpolate */
    if (Stdmag < pMagdep1->magdep_bin_median) {
      /* Interpolate with the next dimmest bin */
      magdep_bin2 = ix + (iy*pMagdepLimits->kxb) + ((imag-1)* pMagdepLimits->kxb * pMagdepLimits->kyb);
      pMagdep2 = &pMagdepTable[magdep_bin2];
    } else {
      /* Interpolate with the next brightest bin */
      magdep_bin2 = ix + (iy*pMagdepLimits->kxb) + ((imag+1)* pMagdepLimits->kxb * pMagdepLimits->kyb);
      pMagdep2 = &pMagdepTable[magdep_bin2];
    }
    if (pMagdep1->magdep_bin_median == pMagdep2->magdep_bin_median) {
      *pmagdep_bin_magcor = pMagdep1->magdep_bin_magcor;
      return(pMagdep1->magdep_bin_magcor);
    } else {
      magcor_value = pMagdep1->magdep_bin_magcor + ((pMagdep2->magdep_bin_magcor - pMagdep1->magdep_bin_magcor)*(Stdmag - pMagdep1->magdep_bin_median))/(pMagdep2->magdep_bin_median - pMagdep1->magdep_bin_median);
      *pmagdep_bin_magcor = magcor_value;
      return(magcor_value);
    }
  }

}


int GetMagdepBin(PMAGDEPLIMITS pMagdepLimits,double X_IMAGE,double Y_IMAGE,double Stdmag,int *pix,int *piy,int *pimag)
{
  int imag;
  int ix;
  int iy;
  int magdep_bin;
  *pix = 0;
  *piy = 0;
  *pimag = 0;



  for (imag = 0; imag < (pMagdepLimits->kmagb-1); imag++) {
    if ((Stdmag >= pMagdepLimits->magdep_bin_edge[imag]) &&
        (Stdmag < pMagdepLimits->magdep_bin_edge[imag+1])) {
      break;
    }
  }
  if (imag == pMagdepLimits->kmagb-1) {
    imag =  pMagdepLimits->kmagb-2;
  }
  if (Y_IMAGE <= pMagdepLimits->ycoord_magdep[0]) {
    iy = 0;
  } else if (Y_IMAGE >= pMagdepLimits->ycoord_magdep[pMagdepLimits->kyb-1]) {
     iy = pMagdepLimits->kyb-1;
  } else {
    for (iy = 0; iy < pMagdepLimits->kyb; iy++) { /*  spatial bin in the Y direction */
      if (Y_IMAGE < (pMagdepLimits->ycoord_magdep[iy]+(pMagdepLimits->dy/2)+0.001)) {
        break;
      }
    }
    if (iy ==  pMagdepLimits->kyb) {
      printf("ERROR: GetMagdepBin failure y\n");
      exit(-1);
    }
  }


  if (X_IMAGE <= pMagdepLimits->xcoord_magdep[0]) {
    ix = 0;
  } else if (X_IMAGE >= pMagdepLimits->xcoord_magdep[pMagdepLimits->kxb-1]) {
    ix = pMagdepLimits->kxb-1;
  } else {
#if 0
    for (ix = 0; ix < pMagdepLimits->kxb; ix++) { /*  spatial bin in the X direction */
      printf("ix %d xcoord_magdep %f X_IMAGE %f dx/2 %f\n",ix,pMagdepLimits->xcoord_magdep[ix],X_IMAGE,pMagdepLimits->dx/2);
    }
#endif
    for (ix = 0; ix < pMagdepLimits->kxb; ix++) { /*  spatial bin in the X direction */
      if (X_IMAGE <  (pMagdepLimits->xcoord_magdep[ix]+(pMagdepLimits->dx/2)+0.001)) {
        break;
      }
    }
    if (ix ==  pMagdepLimits->kxb) {
      printf("ERROR: GetMagdepBin failure x\n");
      exit(-1);
    }
  }
  magdep_bin = ix + (pMagdepLimits->kxb*iy) + (pMagdepLimits->kxb*pMagdepLimits->kyb*imag);
  if ((magdep_bin < 0) || (magdep_bin >= pMagdepLimits->nrecs)) {
    printf("ERROR: GetMagdepBin failure magdep_bin\n");
    exit(-1);

  }
  *pix = ix;
  *piy = iy;
  *pimag = imag;
  return(magdep_bin);
}
int PhotLoadMagdepCorrections(MYSQL *pPhotConnection,int seriesId,int plateNumber,int solutionNumber,char *catalogString,PMAGDEPLIMITS pMagdepLimits,PMAGDEPCORRECTION *ppMagdepTable,int verbose)
{
  int nvals;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int res;
  int numRows = 0;
  int curRow = 0;
  PMAGDEPCORRECTION pMagdepTable = NULL;
  PMAGDEPCORRECTION pMagdep;
  int versionId;
  int versionIdCurrent = 0;
  int returnStatus;

  char queryString[MAX_QUERY_STRING];
  sprintf(queryString,"SELECT versionId,magdep_bin,ixb,iyb,imagb,nstar_magdep,magdep_bin_size,xcoord_magdep,ycoord_magdep,magdep_bin_edge,magdep_bin_median,magdep_bin_magcor,magcal_magdep_rms,magdep_bin_quality from magdepcalibrate%s WHERE seriesId = %d and plateNumber = %d and solutionNumber = %d\n",catalogString,seriesId,plateNumber,solutionNumber);

  res = ExecuteQuery(pPhotConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pPhotConnection);
    if (res_ptr) {
      numRows = mysql_affected_rows(pPhotConnection);
      pMagdepTable = (PMAGDEPCORRECTION)calloc(numRows,sizeof(MAGDEPCORRECTION));
      if (pMagdepTable == NULL) {
        printf("ERROR: Failed to allocate pMagdepTable in PhotLoadMagdepCorrections\n");
        exit(-1);
      }
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        if (curRow >= numRows) {
          printf("ERROR: curRow %d exceeds numRows %d in PhotLoadMagdepCorrections\n");
          exit(-1);
        }
        pMagdep = &pMagdepTable[curRow];
        memset(pMagdep,0,sizeof(MAGDEPCORRECTION));
        if (sqlrow[0]) {
          nvals = sscanf(sqlrow[0],"%d",&versionId);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for versionId in PhotLoadMagdepCorrections\n");
            exit(-1);
          }
          if (versionId <= 0) {
            printf("ERROR: illegal versionId %d in PhotLoadMagdepCorrection\n");
            exit(-1);
          }
        } else {
          printf("ERROR: NULL versionId in PhotLoadMagdepCorrections\n");
          exit(-1);
        }
        if (versionId < versionIdCurrent) {
          continue;
        }
        if (versionId > versionIdCurrent) {
          if (versionIdCurrent > 0) {
            /* We have stale versions - get rid of them now! */
            curRow = 0;
            pMagdep = &pMagdepTable[curRow];
            memset(pMagdep,0,sizeof(MAGDEPCORRECTION));
          }
          versionIdCurrent = versionId;
        }
        if (sqlrow[1]) {

          nvals = sscanf(sqlrow[1],"%d",&pMagdep->magdep_bin);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magdep_bin\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[2]) {

          nvals = sscanf(sqlrow[2],"%d",&pMagdep->ixb);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ixb\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[3]) {

          nvals = sscanf(sqlrow[3],"%d",&pMagdep->iyb);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for iyb\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[4]) {

          nvals = sscanf(sqlrow[4],"%d",&pMagdep->imagb);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for imagb\n");
            continue;
          }
        } else {
          continue;
        }
#if 0
        if ((pMagdep->ixb == 2) &&
            (pMagdep->iyb == 17) &&
            (pMagdep->imagb == 1)) {
          printf("At ixb %d iyb %d imagb %d\n",pMagdep->ixb,pMagdep->iyb,pMagdep->imagb);
        }
#endif

        if (sqlrow[5]) {

          nvals = sscanf(sqlrow[5],"%d",&pMagdep->nstar_magdep);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for nstar_magdep\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[6]) {

          nvals = sscanf(sqlrow[6],"%d",&pMagdep->magdep_bin_size);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magdep_bin_size\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[7]) {

          nvals = sscanf(sqlrow[7],"%lf",&pMagdep->xcoord_magdep);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for xcoord_magdep\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[8]) {

          nvals = sscanf(sqlrow[8],"%lf",&pMagdep->ycoord_magdep);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for ycoord_magdep\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[9]) {

          nvals = sscanf(sqlrow[9],"%lf",&pMagdep->magdep_bin_edge);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magdep_bin_edge\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[10]) {

          nvals = sscanf(sqlrow[10],"%lf",&pMagdep->magdep_bin_median);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magdep_bin_median\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[11]) {

          nvals = sscanf(sqlrow[11],"%lf",&pMagdep->magdep_bin_magcor);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magdep_bin_magcor\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[12]) {

          nvals = sscanf(sqlrow[12],"%lf",&pMagdep->magcal_magdep_rms);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magcal_magdep_rms\n");
            continue;
          }
        } else {
          continue;
        }
        if (sqlrow[13]) {

          nvals = sscanf(sqlrow[13],"%lf",&pMagdep->magdep_bin_quality);
          if (nvals != 1) {
            printf("ERROR: nvals is %d for magdep_bin_quality\n");
            continue;
          }
        } else {
          continue;
        }


        curRow++;
      }
      mysql_free_result(res_ptr);
    } else {
      printf("ERROR: mysql_store_result failed in photometryutils.c line %d\n",__LINE__);
      return(-1);
    }
  } else {
    printf("ERROR: mysql query failed in magdeputil.c line %d\n",__LINE__);
    return(-1);
  }
  if (curRow > 0) {
    returnStatus = ValidateMagdepTable(pMagdepLimits,pMagdepTable,curRow,queryString,ppMagdepTable);
  } else {
    returnStatus = -1;
  }

  if (pMagdepTable != NULL) {
    free(pMagdepTable);
    pMagdepTable = NULL;
  }
  return(returnStatus);


}

