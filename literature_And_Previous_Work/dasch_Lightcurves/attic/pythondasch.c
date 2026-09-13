// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* pythondasch.c
 *
 * This routine provides a Python interface to the dasch photometry database
 *
 *   python setup.py install --prefix=/dasch/install
 *   python setup.py clean
 *   python test.py
 *
 *  Mar 25, 2014 V 1.0  Edward J. Los - Initial version
 *  Mar 26, 2014 V 1.1  Edward J. Los - Make database accesses read-only
 *                                      Correct documentation for GetFileSummaryMagnitudes
 *                                      Correct maximum index check for GetFileStarImage
 *
 *
 */
#define __USE_XOPEN2K8 1
#include <Python.h>
#include "pipelineutils.h"
#include "mysql.h"
#include "photometryutils.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

#define MAX_BUFFER 256

extern GSCBIN gscBin64;
PGSCBIN pGscBin = &gscBin64;

static  FILECOMMON fileCommon;
static  PFILECOMMON pFileCommon = &fileCommon;
static  MYSQL my_connection;
static  MYSQL *pConnection = &my_connection;
static  MYSQL my_phot_connection;
static  MYSQL *pPhotConnection = &my_phot_connection;
static  int catalogNumber = -1;
static  char catalogString[MAX_BUFFER];


static PyObject *
PyGetMaxGSCBin(PyObject* self,PyObject* next)
{
  int gsc_bin_index = pGscBin->total_gsc_bins-1;
  return Py_BuildValue("i",gsc_bin_index);
}

static char PyGetMaxGSCBin_docs[] =
    "GetMaxGSCBin(): Return the maximum gsc_bin_index";


static PyObject *
PyGetGSCBin(PyObject* self,PyObject* next)
{
  int gsc_bin_index = 0;
  double ra = 55.3;
  double dec = 18.5;
  int DecBin;
  int RaBin;
  if (!PyArg_ParseTuple(next,(char *)"dd", &ra,&dec)) return NULL;
  if ((dec < -90.0) || (dec >= 90.0)) {
    PyErr_SetString(PyExc_ValueError,"GetGSCBin: declination is out of bounds");
    return NULL;
  }

  gsc_bin_index = GetGSCBin(pGscBin,ra,dec,&DecBin,&RaBin,"Python DASCH");
  return Py_BuildValue("i",gsc_bin_index);
}

static char PyGetGSCBin_docs[] =
    "GetGSCBin(ra,dec): Return the gsc_bin_index";


static PyObject *
PyGetBinCenter(PyObject* self,PyObject* next)
{
  int gsc_bin_index;
  double ra;
  double dec;
  char errorString[MAX_BUFFER];

  if (!PyArg_ParseTuple(next,(char *)"i", &gsc_bin_index)) return NULL;
  if ((gsc_bin_index < 0) || (gsc_bin_index >= pGscBin->total_gsc_bins)) {
    sprintf(errorString,"GetGSCBin: gsc_bin_index is out of bounds: requres %d to %d",0,pGscBin->total_gsc_bins-1);
    PyErr_SetString(PyExc_ValueError,errorString);
    return NULL;
  }
  if(GetBinCenter(pGscBin,gsc_bin_index,&ra,&dec,"Python DASCH") != 0) {
    PyErr_SetString(PyExc_RuntimeError,errorString);
    return NULL;

  }

  return Py_BuildValue("(dd)",ra,dec);
}

static char PyGetBinCenter_docs[] =
    "GetBinCenter(gsc_bin_index): Return the ra and dec at the bin center as the tuple (ra,dec)";

static PyObject *
PyGetFileSummaryMagnitudes(PyObject* self,PyObject* next)
{
  int gsc_bin_index;
  char errorString[MAX_BUFFER];
  int curMagnitudes;

  if ((catalogNumber < 0) || (catalogNumber > MAX_CATALOG_NUMBER)) {
    PyErr_SetString(PyExc_RuntimeError,"GetFileSummaryMagnitudes: database access not initiated with dasch.InitDataBaseAccess(catalog)");
    return NULL;

  }


  if (!PyArg_ParseTuple(next,(char *)"i", &gsc_bin_index)) return NULL;
  if ((gsc_bin_index < 0) || (gsc_bin_index >= pGscBin->total_gsc_bins)) {
    sprintf(errorString,"GetFileSummaryMagnitudes: gsc_bin_index is out of bounds: requres %d to %d",0,pGscBin->total_gsc_bins-1);
    PyErr_SetString(PyExc_ValueError,errorString);
    return NULL;
  }


  curMagnitudes = GetFileSummaryMagnitudes(pGscBin,pFileCommon,gsc_bin_index,gsc_bin_index,0,0,0,0,NULL,catalogString,0,NULL,1,0);

  

  return Py_BuildValue("i",curMagnitudes);
}

static char PyGetFileSummaryMagnitudes_docs[] =
    "GetFileSummaryMagnitudes(gsc_bin_index): Reads the magnitude data for the gsc_bin_index and returns a count of magnitudes";

static PyObject *
PyGetFileStarImage(PyObject* self,PyObject* next)
{
  int star_index;
  char errorString[MAX_BUFFER];
  char REF[MAX_REF];
  PFILESTARIMAGE pFileStarImage = NULL;
  double short_year;
  char short_plate[MAX_BUFFER];
  double magcal_local_rms;
  double magcal_magdep;
  double limiting_mag_local;

  if ((catalogNumber < 0) || (catalogNumber > MAX_CATALOG_NUMBER)) {
    PyErr_SetString(PyExc_RuntimeError,"GetFileStarImage: database access not initiated with dasch.InitDataBaseAccess(catalog)");
    return NULL;

  }


  if (!PyArg_ParseTuple(next,(char *)"i", &star_index)) return NULL;
  if ((star_index < 0) || (star_index >= pFileCommon->curMagnitudes)) {
    sprintf(errorString,"GetFileStarImage: star_index is out of bounds: requres %d to %d",0,pFileCommon->curMagnitudes-1);
    PyErr_SetString(PyExc_ValueError,errorString);
    return NULL;
  }

  pFileStarImage = &pFileCommon->magnitudeBuffer[star_index];
  short_year = jd2ep(pFileStarImage->Date);
  GetREF(pFileStarImage->REFNumber,REF,1,1);
  if (pFileStarImage->solutionNumber == 0) {
    sprintf(short_plate,"%s%05d",GetSeriesString(pFileStarImage->seriesId,0),pFileStarImage->plateNumber);
  } else {
    sprintf(short_plate,"%s%05d_s%d",GetSeriesString(pFileStarImage->seriesId,0),pFileStarImage->plateNumber,pFileStarImage->solutionNumber);
  }

  magcal_local_rms = pFileStarImage->magcal_local_rms;
  magcal_magdep = pFileStarImage->magcal_magdep;
  limiting_mag_local = pFileStarImage->limiting_mag_local;

  return Py_BuildValue("{sssdsdsdsdsdsdsdsdsdsssisi}",
                       "REF",REF,
                       "Date",pFileStarImage->Date,
                       "magcal_magdep",magcal_magdep,
                       "magcal_local_rms",magcal_local_rms,
                       "limiting_mag_local",limiting_mag_local,
                       "ra",pFileStarImage->ra,
                       "dec",pFileStarImage->dec,
                       "THETA_J2000",pFileStarImage->THETA_J2000,
                       "ELLIPTICITY",pFileStarImage->ELLIPTICITY,
                       "year",short_year,
                       "Plate",short_plate,
                       "versionId",pFileStarImage->versionId,
                       "AFLAGS",pFileStarImage->AFLAGS);

}

static char PyGetFileStarImage_docs[] =
    "GetFileStarImage(index): Return a dictionary for the star index read from the last GetFileSummaryMagnitudes() call";


static PyObject *
PyInitDatabaseAccess(PyObject* self,PyObject* next)
{
  char *qualifier = NULL;
  char errorString[MAX_BUFFER];
  char *mysqlhost;
  char *username;
  char *password;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;


  /* Connect to the database */

  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    PyErr_SetString(PyExc_ValueError,"ERROR: DASCH_MYSQLHOST is not defined");
    return NULL;
  }
  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    PyErr_SetString(PyExc_ValueError,"ERROR: DASCH_USERNAME is not defined");
    return NULL;
  }
  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    PyErr_SetString(PyExc_ValueError,"ERROR: DASCH_PASSWORD is not defined");
    return NULL;
  }

  mysql_init(pConnection);


  if (!mysql_real_connect(pConnection,mysqlhost,username,password,"scanner",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      sprintf(errorString,"ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
    } else {
      strcpy(errorString,"Failed to connect to the scanner database");
    }    
    PyErr_SetString(PyExc_ValueError,"errorString");
    return NULL;
  }


  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    PyErr_SetString(PyExc_ValueError,"DASCH_PHOT_MYSQLHOST is not defined");
    return NULL;
  }

  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    PyErr_SetString(PyExc_ValueError,"DASCH_PHOT_USERNAME is not defined");
    return NULL;
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    PyErr_SetString(PyExc_ValueError,"DASCH_PHOT_PASSWORD is not defined");
    return NULL;
  }
  mysql_init(pPhotConnection);


  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      sprintf(errorString,"ERROR: MySQL error %d: %s\n",mysql_errno(pConnection),mysql_error(pConnection));
    } else {
      strcpy(errorString,"Failed to connect to the photometry database");
    }    
    PyErr_SetString(PyExc_ValueError,"errorString");
    return NULL;
  }



  if (!PyArg_ParseTuple(next,(char *)"s", &qualifier)) return NULL;
  catalogNumber = GetCatalogNumber(qualifier);
  if (qualifier == NULL) {
    sprintf(errorString,"InitDatabaseAccess: NULL catalog name");
    PyErr_SetString(PyExc_ValueError,errorString);
    return NULL;
  }    

  if (catalogNumber < 0) {
    sprintf(errorString,"InitDatabaseAccess: Illegal catalog name %s",qualifier);
    PyErr_SetString(PyExc_ValueError,errorString);
    return NULL;
  }    
  sprintf(catalogString,"%d",catalogNumber);

  InitSeriesTable(pConnection,pPhotConnection);
  InitFileCommon(stdout,pFileCommon,pConnection,pPhotConnection,catalogNumber,1);
  InitMaxPlateNumber(pConnection,pFileCommon->maxPlateNumber);
  InitCatalogAccess(pFileCommon);

  

  return Py_BuildValue("i",0);
}

static char PyInitDatabaseAccess_docs[] =
    "InitDatabaseAccess(catalog): Initializes access to database catalog 'gsc2.3.2', 'kepler', or 'apass'";


/* Can also use METH_KEYWORDS with a different parsing function */
static PyMethodDef dasch_funcs[] = {
  {"GetMaxGSCBin",(PyCFunction)PyGetMaxGSCBin,METH_NOARGS,PyGetMaxGSCBin_docs},
  {"GetGSCBin",(PyCFunction)PyGetGSCBin,METH_VARARGS,PyGetGSCBin_docs},
  {"GetBinCenter",(PyCFunction)PyGetBinCenter,METH_VARARGS,PyGetBinCenter_docs},
  {"GetFileSummaryMagnitudes",(PyCFunction)PyGetFileSummaryMagnitudes,METH_VARARGS,PyGetFileSummaryMagnitudes_docs},
  {"GetFileStarImage",(PyCFunction)PyGetFileStarImage,METH_VARARGS,PyGetFileStarImage_docs},
  {"InitDatabaseAccess",(PyCFunction)PyInitDatabaseAccess,METH_VARARGS,PyInitDatabaseAccess_docs},
   {NULL}
};

void
initdasch(void)
{
Py_InitModule3("dasch",dasch_funcs,
              "Python interface to the DASCH database");



}
