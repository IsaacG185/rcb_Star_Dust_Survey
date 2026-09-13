// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* votable.c
 *
 * Convert our starbase tables to VOTables using the "variables" table in the auxscanner database
 *
 * cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include  -L /dasch/install/lib -lm -lcfitsio   votable.c pipelineutils.a -ltable -lutil  -lwcs -o votable -L/usr/lib64/mysql -L/usr/lib/mysql -lmysqlclient
 *
 * STARBASE to VOTABLE:
 *    votable -i /dasch/Pipeline/ingest/i31090_00_01r180ww_allobjects.db -o /dasch/junk/i31090.xml
 *    votable -i /dasch/Pipeline/ingest/mc39048_01_01r270ww_allobjects.db  -o /dasch/junk/mc39048.xml
 *    votable -i /dasch/Pipeline/match/match_i31090_00_01r180ww_tnx_u.db -o /dasch/junk/match_i31090.xml
 *
 * VOTABLE VERIFICATON:
 *    java -jar ~/Download/stilts.jar votlint version=1.1 votable=i31090.xml
 * VOTABLE to FITS (does not covert PARAM fields)
 *    java -jar ~/Download/stilts.jar tcopy i31090.xml i31090.fits
 *
 * Mar 10, 2010 Edward J. Los - Original Version
 * Mar 11, 2010 Edward J. Los - make "plate" a synonym for "Plate"
 * Jan 17, 2014 Edward J. Los - change to the Odyssey website for variable references
 * Aug 19, 2015 Edward J. Los - Increment MAX_VARIABLE_NAME_SIZE to 25
 * May 29, 2018 Edward J. Los - Remove 'quality' parameter because of units conflict
 *
 */


#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>

#include <mysql.h>

#include <table.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "photometryutils.h"


#define MAX_VARIABLE_NAME_SIZE   25
#define DATATYPE_CHAR      1
#define DATATYPE_DOUBLE    2
#define DATATYPE_INT       3
#define DATATYPE_SHORT     4
#define DATATYPE_LONG      5
#define DATATYPE_ENUM      6
#define DATATYPE_SET       7
#define MAX_DATATYPE       7
#define MAX_SPACES 100
/* datatypes stored in the variables table */
char *datatypeString1[MAX_DATATYPE+1] = {"ERROR","char","double","int","short","long","enum","set"};
/* equivalent datatypes for VOTables */
char *datatypeString2[MAX_DATATYPE+1] = {"ERROR","char","double","int","short","long","char","char"};


#define MAX_UNIT_SIZE     16
#define MAX_UCD_SIZE     132
#define MAX_VARIABLE_DESCRIPTION_SIZE  512
#define MAX_BUFFER 512

typedef struct _variablesEntry {
  int datatype;          /* Data storage type */
  int arraysize;         /* Length of string */
  int width;             /* Field width */
  int precision;         /* Field precision */
  char variableName[MAX_VARIABLE_NAME_SIZE]; /* Column header */
  char unit[MAX_UNIT_SIZE];   /* physical unit */
  char ucd[MAX_UCD_SIZE];     /* Unified Content Descriptor */
  char variableDescription[MAX_VARIABLE_DESCRIPTION_SIZE]; /* Description */
} VARIABLESENTRY,*PVARIABLESENTRY;

#define TABLEVALUE_FIELD 1
#define TABLEVALUE_PARAM 2

typedef struct _columnEntry {
  PVARIABLESENTRY pVariablesEntry;
  int tableValue;
} COLUMNENTRY,*PCOLUMNENTRY;
char *spaces(int numSpaces) {
  static int curSpaces = -1;
  static char spacesBuffer[MAX_SPACES];
  int index;

  if (numSpaces != curSpaces) {
    if (numSpaces >= MAX_SPACES) {
      printf("ERROR: spaces %d exceeds MAX_SPACES %d\n",spaces,MAX_SPACES);
      exit(-1);
    }
    for (index = 0; index < numSpaces; index++) {
      spacesBuffer[index] = ' ';
    }
    spacesBuffer[numSpaces] = 0;
    curSpaces = numSpaces;
  }

  return(spacesBuffer);

}


void  LoadVariables(PVARIABLESENTRY *pvariables_table,size_t *pVariables_nrecs)
{
  PVARIABLESENTRY variables_table = *pvariables_table;
  PVARIABLESENTRY pVariablesEntry;
  size_t variables_nrecs = *pVariables_nrecs;
  char *mysqlhost;
  char *username;
  char *password;
  char *auxscanner;
  MYSQL aux_connection;
  MYSQL *pAuxConnection = &aux_connection;
  char queryString[MAX_QUERY_STRING];
  int nvals;
  int res;
  MYSQL_RES *res_ptr;
  MYSQL_ROW sqlrow;
  int result_size;
  int foundCount = 0;

  *pVariables_nrecs = 0;

  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    printf("ERROR: DASCH_MYSQLHOST is not defined\n");
    exit(-1);
  }

  username = getenv("DASCH_USERNAME");
  if (username == NULL) {
    printf("ERROR: DASCH_USERNAME is not defined\n");
    exit(-1);
  }

  password = getenv("DASCH_PASSWORD");
  if (password == NULL) {
    printf("ERROR: DASCH_PASSWORD is not defined\n");
    exit(-1);
  }

  auxscanner = getenv("DASCH_AUXSCANNER");
  if (auxscanner == NULL) {
    printf("ERROR: DASCH_AUXSCANNER is not defined\n");
    exit(-1);
  }

  mysql_init(pAuxConnection);

  if (!mysql_real_connect(pAuxConnection,mysqlhost,username,password,auxscanner,0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pAuxConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pAuxConnection),mysql_error(pAuxConnection));
    }
    exit(-1);
  }

  sprintf(queryString,"SELECT variableName,datatype+0,arraysize,width,decimalPlaces,unit,ucd,variableDescription from variables order by variableName;");
  res = ExecuteQuery(pAuxConnection,queryString);
  if (!res) {
    res_ptr = mysql_store_result(pAuxConnection);
    if (res_ptr) {
      result_size = (int)mysql_num_rows(res_ptr);
      if (result_size > 0) {
        if (variables_table != NULL) {
          free(variables_table);
        }
        variables_table = (PVARIABLESENTRY)calloc(result_size,sizeof(VARIABLESENTRY));
        *pvariables_table = variables_table;
        if (variables_table == NULL) {
          printf("ERROR: allocation failure for variables_table of size %d\n",result_size);
          exit(-1);
        }

      }
      while ((sqlrow = mysql_fetch_row(res_ptr))) {
        pVariablesEntry = &variables_table[foundCount];
	memset(pVariablesEntry,0,sizeof(VARIABLESENTRY));
	if (sqlrow[0] == NULL) {
	  printf("ERROR: variableName is NULL in variables table\n");
	  continue;
	}
	if (strlen(sqlrow[0]) >= MAX_VARIABLE_NAME_SIZE-1) {
	  printf("ERROR: MAX_VARIABLE_NAME_SIZE %d is smaller than %d %s\n",MAX_VARIABLE_NAME_SIZE,strlen(sqlrow[0]),sqlrow[0]);
	  exit(-1);
	}
	strcpy(pVariablesEntry->variableName,sqlrow[0]);
	if (sqlrow[1] != NULL) {
	  nvals = sscanf(sqlrow[1],"%d",&pVariablesEntry->datatype);
	  if (nvals != 1) {
	    printf("ERROR: nvals is %d for datatype in variables table for %s\n",nvals,sqlrow[0]);
	    continue;
	  }
	  if ((pVariablesEntry->datatype < 0) ||
	      (pVariablesEntry->datatype > MAX_DATATYPE)) {
	    printf("ERROR: unrecognized datatype %d for %s\n",pVariablesEntry->datatype,sqlrow[0]);
	    continue;
	  }
	  if (pVariablesEntry->datatype == 0) {
	    printf("ERROR: null datatype %d for %s\n",pVariablesEntry->datatype,sqlrow[0]);
	  }
	} else {
	  printf("ERROR: null datatype %d for %s\n",pVariablesEntry->datatype,sqlrow[0]);
	}
	if (sqlrow[2] != NULL) {
	  nvals = sscanf(sqlrow[2],"%d",&pVariablesEntry->arraysize);
	  if (nvals != 1) {
	    printf("ERROR: nvals is %d for arraysize in variables table for %s\n",nvals,sqlrow[0]);
	    continue;
	  }
	} else {
	  if ((pVariablesEntry->datatype == DATATYPE_CHAR) ||
	      (pVariablesEntry->datatype == DATATYPE_ENUM) ||
	      (pVariablesEntry->datatype == DATATYPE_SET)) {
	    printf("ERROR arraysize is null for %s of type %s\n",sqlrow[0],datatypeString1[pVariablesEntry->datatype]);
	  }

	}
	if (sqlrow[3] != NULL) {
	  nvals = sscanf(sqlrow[3],"%d",&pVariablesEntry->width);
	  if (nvals != 1) {
	    printf("ERROR: nvals is %d for width in variables table for %s\n",nvals,sqlrow[0]);
	    continue;
	  }
	}
	if (sqlrow[4] != NULL) {
	  nvals = sscanf(sqlrow[4],"%d",&pVariablesEntry->precision);
	  if (nvals != 1) {
	    printf("ERROR: nvals is %d for precision in variables table for %s\n",nvals,sqlrow[0]);
	    continue;
	  }
	}

	if (sqlrow[5] != NULL) {
	  if (strlen(sqlrow[5]) >= MAX_UNIT_SIZE-1) {
	    printf("ERROR: MAX_UNIT_SIZE %d is smaller than %d unit %s\n",MAX_UNIT_SIZE,strlen(sqlrow[5]),sqlrow[5]);
	    exit(-1);
	  }
	  strcpy(pVariablesEntry->unit,sqlrow[5]);
	}

	if (sqlrow[6] != NULL) {
	  if (strlen(sqlrow[6]) >= MAX_UCD_SIZE-1) {
	    printf("ERROR: MAX_UCD_SIZE %d is smaller than %d ucd %s\n",MAX_UCD_SIZE,strlen(sqlrow[6]),sqlrow[6]);
	    exit(-1);
	  }
	  strcpy(pVariablesEntry->ucd,sqlrow[6]);
	}

	if (sqlrow[7] != NULL) {
	  if (strlen(sqlrow[7]) >= MAX_VARIABLE_DESCRIPTION_SIZE-1) {
	    printf("ERROR: MAX_VARIABLE_DESCRIPTION_SIZE %d is smaller than %d variableDescription %s\n",MAX_VARIABLE_DESCRIPTION_SIZE,strlen(sqlrow[7]),sqlrow[7]);
	    exit(-1);
	  }
	  strcpy(pVariablesEntry->variableDescription,sqlrow[7]);
	} else {
	  printf("ERROR variableDescription is null for %s\n",sqlrow[0]);
	}
	foundCount++;
	if (foundCount > result_size) {
	  printf("ERROR: foundCount %d exceeds allocCount %d\n",foundCount,result_size);
	  exit(-1);
	}
      }
      mysql_free_result(res_ptr);
    }
  } else {
    printf("ERROR Select error %d: %s res %d\n",mysql_errno(pAuxConnection),mysql_error(pAuxConnection),res);
    exit(-1);
  }
  if (foundCount != result_size) {
    printf("ERROR: foundCount %d not equal to allocCount %d\n",foundCount,result_size);
    exit(-1);
  }

  mysql_close(pAuxConnection);

  variables_nrecs = foundCount;
  *pVariables_nrecs = variables_nrecs;
  *pvariables_table = variables_table;
  return;
}

PVARIABLESENTRY GetVariableEntry(PVARIABLESENTRY variables_table,size_t variables_nrecs,char *nameString)
{
  size_t index;
  PVARIABLESENTRY pVariablesEntry;
  char tmpString[MAX_VARIABLE_NAME_SIZE+2];
  if (strlen(nameString) > MAX_VARIABLE_NAME_SIZE) {
    printf("ERROR: column header %s  %d exceeds MAX_VARIABLE_NAME_SIZE %d\n",nameString,strlen(nameString),MAX_VARIABLE_NAME_SIZE);
    return(NULL);
  }

  if (strcmp(nameString,"plate") == 0) {
    strcpy(tmpString,"Plate");
  } else if (strcmp(nameString,"Date") == 0) {
    strcpy(tmpString,"ExposureDate");
  } else {
    strcpy(tmpString,nameString);
  }
  for (index = 0; index < variables_nrecs; index++) {
    pVariablesEntry = &variables_table[index];
    if (strcmp(pVariablesEntry->variableName,tmpString) == 0) {
      return(pVariablesEntry);
    }
  }
  printf("ERROR: column header %s is not in the variables table\n",nameString);
  return(NULL);

}

/* Note: if colPtr is NULL, then this routine writes a FIELD description.  Otherwise
 * it writes a PARAM description using colPtr as a pointer to the string containing
 * the PARAM value
 */
int WriteParamOrValue(int numSpaces,FILE *outHandle,char *colPtr,PVARIABLESENTRY pVariablesEntry)
{
  char *webPtr;
  char *termPtr;
  char newDescription[MAX_VARIABLE_DESCRIPTION_SIZE];
  if (colPtr == NULL) {
    fprintf(outHandle,"%s<FIELD ",spaces(numSpaces));
  } else {
    fprintf(outHandle,"%s<PARAM ",spaces(numSpaces));
  }

  fprintf(outHandle,"datatype=\"%s\" ",datatypeString2[pVariablesEntry->datatype]);
  if (pVariablesEntry->arraysize != 0) {
    fprintf(outHandle,"arraysize=\"%d\" ",pVariablesEntry->arraysize);
  }
  if (pVariablesEntry->width != 0) {
    fprintf(outHandle,"width=\"%d\" ",pVariablesEntry->width);
  }
  if (pVariablesEntry->precision != 0) {
    fprintf(outHandle,"precision=\"%d\" ",pVariablesEntry->precision);
  }
  if (pVariablesEntry->unit[0] != 0) {
    fprintf(outHandle,"unit=\"%s\" ",pVariablesEntry->unit);
  }
  if (pVariablesEntry->ucd[0] != 0) {
    fprintf(outHandle,"ucd=\"%s\" ",pVariablesEntry->ucd);
  }
  fprintf(outHandle,"name=\"%s\" ",pVariablesEntry->variableName);
  if (colPtr != NULL) {
    fprintf(outHandle,"value=\"%s\" ",colPtr);
  }

  if (pVariablesEntry->variableDescription[0] == 0) {
    fprintf(outHandle,"/>\n");
  } else {
    if ((webPtr = strstr(pVariablesEntry->variableDescription,"WEB{")) != NULL) {
      *webPtr = 0;
      webPtr += 4;
      termPtr = strstr(webPtr,"}");
      if (termPtr) {
	*termPtr = 0;
      } else {
	printf("ERROR: %s has no terminator on %s\n",pVariablesEntry->variableDescription,webPtr);
      }
#if 1
      sprintf(newDescription,"%s http://dasch.rc.fas.harvard.edu/%s",pVariablesEntry->variableDescription,webPtr);
#else
      sprintf(newDescription,"%s http://hea-www.harvard.edu/DASCH/%s",pVariablesEntry->variableDescription,webPtr);
#endif
    } else {
      strcpy(newDescription,pVariablesEntry->variableDescription);
    }

    fprintf(outHandle,">\n");
    fprintf(outHandle,"%s<DESCRIPTION>%s</DESCRIPTION>\n",spaces(numSpaces+1),newDescription);
    if (colPtr == NULL) {
      fprintf(outHandle,"%s</FIELD>\n",spaces(numSpaces));
    } else {
      fprintf(outHandle,"%s</PARAM>\n",spaces(numSpaces));
    }
  }
  return(0);
}



int main(int argc,char *argv[])
{
  char *argstr;
  FILE *outHandle = NULL;
  int errorFlag = 0;
  int verbose = 0;
  char outfile[MAX_BUFFER];
  char cmdchar;

  size_t variables_nrecs;
  PVARIABLESENTRY variables_table = NULL;
  PVARIABLESENTRY pVariablesEntry;


  File input_handle = NULL;
  char input_name[MAX_BUFFER];
  char resource_name[MAX_BUFFER];
  char base_name[MAX_BUFFER];
  TableHead input_header = NULL;
  size_t input_nrecs;
  size_t input_index;
  int columnIndex;
  int numColumns;
  PCOLUMNENTRY pColumnEntry;
  PCOLUMNENTRY pColumnTable = NULL;
  int havePlate = 0;
  int typeAllobjects = 0;
  char series[MAX_SERIES_STRING];
  int plateNumber;
  int mosaicNumber;
  int binning;
  int rotation;
  int firstDataRow = 1;
  TableRow input_row = NULL;
  char *colPtr;
  char *charPtr;
  char variable_name[MAX_BUFFER];
  char variable_buffer[MAX_BUFFER];
  SERIES_HEADER seriesHeader;
  PSERIES_HEADER pSeriesHeader = &seriesHeader;

  PLATE_ENTRY plateEntry;
  PPLATE_ENTRY pPlateEntry = &plateEntry;

  MYSQL connection;
  MYSQL *pConnection = &connection;
  double wLongitude;
  double latitude;
  double elevation;
  char name[MAX_LOCATION_NAME];
  int numSpaces = 0;

  dasch_init_scandb(pConnection);

  input_name[0] = 0;
  outfile[0] = 0;

  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d\n",argstr,argc);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'i': /* input file */
        case 'I':
          argc--;
          if (argc < 1) {
            fprintf(stderr,"ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            strncpy(input_name,*++argv,MAX_BUFFER-2);
            if (strlen(input_name) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
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
            strncpy(outfile,*++argv,MAX_BUFFER-2);
            if (strlen(outfile) >= MAX_BUFFER-3) {
              fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
            }
          }
          break;


        case 'v': /* verbose */
        case 'V':
          verbose = 1;
          break;


        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;

        }

      }

    }
  }



  if (outfile[0] == 0) {
    fprintf(stderr,"ERROR: No output filename was specified\n");
    errorFlag = 1;
  } else {


    outHandle = fopen(outfile,"wt");
    if (outHandle == NULL) {
      errorFlag = 1;
      fprintf(stderr,"ERROR: Failed to open the output file %s\n",outfile);
    } else {
      if (verbose) {
        fprintf(stderr,"Output file %s\n",outfile);
      }
    }
  }




  if (input_name[0] == 0) {
    fprintf(stderr,"ERROR: No input filename was specified\n");
    errorFlag = 1;
  } else {

    input_handle = Open(input_name,"r");
    if (input_handle == NULL) {
      fprintf(stderr,"ERROR: Failed to find the input name %s\n",input_name);
      errorFlag = 1;
    } else {
      if (verbose) {
        fprintf(stderr,"Found input file %s\n",input_name);
      }
    }
  }

  if (errorFlag) {
    printf("Usage: votable options\n");
    printf("  options: -v verbose\n");
    printf("           -i <input file>\n");
    printf("           -o <output file>\n");
    return(-1);
  }
  LoadVariables(&variables_table,&variables_nrecs);

  charPtr = strrchr(input_name,'/');
  if (charPtr != NULL) {
    charPtr++;
    strcpy(base_name,charPtr);
  } else {
    strcpy(base_name,input_name);
  }


  charPtr = strstr(base_name,"match_");
  if (charPtr) {
    charPtr += 6;
    strcpy(resource_name,charPtr);

  } else {
    strcpy(resource_name,base_name);
  }


  if (ParseFilename(resource_name,series,&plateNumber,&mosaicNumber,&binning,&rotation)) {
    if (GetPlateInfo(pConnection,series,plateNumber,pPlateEntry) == 1) {


      havePlate = 1;
      if (strstr(resource_name,"_allobjects.db")) {
        typeAllobjects = 1;
      }
      charPtr = strstr(resource_name,"ww");
      if (charPtr) {
        charPtr += 2;
        *charPtr = 0;

      } else {
        printf("ERROR: resource name %s does not contain 'ww'\n",resource_name);
        exit(-1);
      }

      if (GetSeriesInfo(pConnection,series,-1,pSeriesHeader) != 1) {
        printf("ERROR: Failed to get series information for %s\n",series);
        exit(-1);
      }
      if (GetLongitude(pConnection,pPlateEntry->locationId,&wLongitude,&latitude,&elevation,name,MAX_LOCATION_NAME) != 1) {
        printf("ERROR: Failed to get observatory name\n");
        exit(-1);
      }

    } else {
      printf("WARNING: Input filename has no recognizable plate number (%s%05d)\n",series,plateNumber);
    }
  }

  mysql_close(pConnection);
  fprintf(outHandle,"<?xml version=\"1.0\"?>\n");
  fprintf(outHandle,"<!DOCTYPE VOTABLE SYSTEM \"http://us-vo.org/xml/VOTable.dtd\">\n");
  fprintf(outHandle,"<VOTABLE version=\"1.1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns=\"http://www.ivoa.net/xml/VOTable/v1.1\" xsi:schemaLocation=\"http://www.ivoa.net/xml/VOTable/v1.1 http://www.ivoa.net/xml/VOTable/v1.1\" >\n");
  fprintf(outHandle," <DESCRIPTION>\n");
  fprintf(outHandle,"  DASCH: Digital Access to a Sky Century @ Harvard\n");
  fprintf(outHandle,"  Harvard College Observatory astronomical plate photometry data\n");
#if 1
  fprintf(outHandle,"  http://dasch.rc.fas.harvard.edu/\n");
#else
  fprintf(outHandle,"  http://hea-www.harvard.edu/DASCH/\n");
#endif
  fprintf(outHandle," </DESCRIPTION>\n");
  fprintf(outHandle,"<!-- VOTable description at http://cdsweb.u-strasbg.fr/doc-cds/VOTable/ -->\n");

  fprintf(outHandle," <DEFINITIONS>\n");
  fprintf(outHandle,"  <COOSYS ID=\"J2000\" system=\"eq_FK5\" equinox=\"J2000\"/>\n");
  fprintf(outHandle," </DEFINITIONS>\n");
  numSpaces = 1;

  /* Now read in the input */

  input_header = table_header(input_handle,TABLE_PARSE);
  if (input_header == NULL) {
    fprintf(stderr,"ERROR: Failed to read header for %s\n",input_name);
    return(-1);
  }
  numColumns = input_header->header->ncol;
  if (numColumns < 1) {
    printf("ERROR: Table has %d columns\n",numColumns);
    exit(-1);
  }
  pColumnTable = (PCOLUMNENTRY)calloc(numColumns,sizeof(COLUMNENTRY));
  if (pColumnTable == NULL) {
    printf("ERROR: failed to allocate pColumnTable\n");
    exit(-1);
  }

  if (havePlate) {
    fprintf(outHandle,"%s<RESOURCE name=\"%s\">\n",spaces(numSpaces),resource_name);
    numSpaces++;
    fprintf(outHandle,"%s<DESCRIPTION>Data for a single plate in the Harvard College Observatory collection</DESCRIPTION>\n",spaces(numSpaces));
    strcpy(variable_name,"series");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    WriteParamOrValue(numSpaces,outHandle,series,pVariablesEntry);


    strcpy(variable_name,"description");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    WriteParamOrValue(numSpaces,outHandle,pSeriesHeader->description,pVariablesEntry);


    strcpy(variable_name,"aperture");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    sprintf(variable_buffer,"%lf",pSeriesHeader->aperture);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);

    strcpy(variable_name,"nominalPlateScale");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    sprintf(variable_buffer,"%lf",pSeriesHeader->nominalPlateScale);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);


    strcpy(variable_name,"name");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    WriteParamOrValue(numSpaces,outHandle,name,pVariablesEntry);


    strcpy(variable_name,"wLongitude");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    sprintf(variable_buffer,"%lf",wLongitude);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);

    strcpy(variable_name,"latitude");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    sprintf(variable_buffer,"%lf",latitude);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);

    strcpy(variable_name,"elevation");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    sprintf(variable_buffer,"%lf",elevation);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);



    strcpy(variable_name,"plateNumber");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    sprintf(variable_buffer,"%d",plateNumber);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);


    strcpy(variable_name,"plateClass");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    WriteParamOrValue(numSpaces,outHandle,pPlateEntry->plateClass,pVariablesEntry);


#if 0 /* Removed on May 29, 2018.  pPlatEntry returns quality as a string while the varables table lists the integer form */
    strcpy(variable_name,"quality");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    WriteParamOrValue(numSpaces,outHandle,pPlateEntry->quality,pVariablesEntry);
#endif


    strcpy(variable_name,"mosaicNumber");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    sprintf(variable_buffer,"%d",mosaicNumber);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);

    strcpy(variable_name,"binning");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",binning);
      exit(-1);
    }
    sprintf(variable_buffer,"%d",binning);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);

    strcpy(variable_name,"rotation");
    if ((pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,variable_name)) == NULL) {
      printf("ERROR: failed to get variable entry for %s\n",variable_name);
      exit(-1);
    }
    sprintf(variable_buffer,"%d",rotation);
    WriteParamOrValue(numSpaces,outHandle,variable_buffer,pVariablesEntry);



  } else {

    fprintf(outHandle,"%s<RESOURCE >\n",spaces(numSpaces));
    numSpaces++;
  }


  fprintf(outHandle,"%s<TABLE name=\"%s\" ID=\"%s\" >\n",spaces(numSpaces),base_name,base_name);
  numSpaces++;
  if (typeAllobjects) {
    fprintf(outHandle,"%s<DESCRIPTION>Table of magnitude results from the plate photometry pipeline.</DESCRIPTION>\n",spaces(numSpaces));


  }

  for (columnIndex = 1; columnIndex <= numColumns; columnIndex++) {
    pColumnEntry = &pColumnTable[columnIndex-1];

    pVariablesEntry = GetVariableEntry(variables_table,variables_nrecs,input_header->header->column[columnIndex]);
    if (pVariablesEntry == NULL) {
      errorFlag = 1;
      printf("ERROR: unrecognized column %d %s\n",columnIndex,input_header->header->column[columnIndex]);
    } else {
      pColumnEntry->pVariablesEntry = pVariablesEntry;
      if ((typeAllobjects == 1) &&
          (strstr(pVariablesEntry->variableName,"THRESHOLD") != NULL)) {
        pColumnEntry->tableValue = TABLEVALUE_PARAM;
      } else {
        pColumnEntry->tableValue = TABLEVALUE_FIELD;
        WriteParamOrValue(numSpaces,outHandle,NULL,pVariablesEntry);
      }


#if 0
      printf("Name: %-24s  Datatype %s\n",pVariablesEntry->variableName,datatypeString1[pVariablesEntry->datatype]);
#endif
    }
  }
  if (errorFlag) {
    exit(-1);
  }
  /* Now read in our first row and print out any parameters */
  input_row = table_rowget(input_handle,input_header,input_row,NULL,NULL,0);
  if (input_row == NULL) {
    printf("ERROR: starbase file has only a header\n");
    exit(-1);
  }
  for (columnIndex = 1; columnIndex <= numColumns; columnIndex++) {
    pColumnEntry = &pColumnTable[columnIndex-1];
    if (pColumnEntry->tableValue == TABLEVALUE_PARAM) {
      pVariablesEntry = pColumnEntry->pVariablesEntry;
      colPtr = table_colval(input_row,columnIndex);
      WriteParamOrValue(numSpaces,outHandle,colPtr,pVariablesEntry);

    }
  }



  fprintf(outHandle,"%s<DATA>\n",spaces(numSpaces));
  numSpaces++;
  fprintf(outHandle,"%s<TABLEDATA>\n",spaces(numSpaces));
  numSpaces++;

  /* Now read in each row and write it out to the votable */
  while (1) {
    if (firstDataRow == 0) {
      input_row = table_rowget(input_handle,input_header,input_row,NULL,NULL,0);
      if (input_row == NULL) {
        break;
      }
    } else {
      firstDataRow = 0;
    }
    fprintf(outHandle,"%s<TR>\n",spaces(numSpaces));
    numSpaces++;

    for (columnIndex = 1; columnIndex <= numColumns; columnIndex++) {
      pColumnEntry = &pColumnTable[columnIndex-1];
      if (pColumnEntry->tableValue == TABLEVALUE_FIELD) {
        colPtr = table_colval(input_row,columnIndex);
        fprintf(outHandle,"%s<TD>%s</TD>\n",spaces(numSpaces),colPtr);
      }
    }
    numSpaces--;
    fprintf(outHandle,"%s</TR>\n",spaces(numSpaces));

  }


  if (input_row != NULL) {
    table_rowfree(input_row);
  }
  numSpaces--;
  fprintf(outHandle,"%s</TABLEDATA>\n",spaces(numSpaces));
  numSpaces--;
  fprintf(outHandle,"%s</DATA>\n",spaces(numSpaces));
  numSpaces--;
  fprintf(outHandle,"%s</TABLE>\n",spaces(numSpaces));
  numSpaces--;
  fprintf(outHandle,"%s</RESOURCE>\n",spaces(numSpaces));
  numSpaces--;
  fprintf(outHandle,"%s</VOTABLE>\n",spaces(numSpaces));

  fclose(outHandle);
  if (pColumnTable != NULL) {
    free(pColumnTable);
  }


  if (variables_table != NULL) {
    free(variables_table);
  }
  if (input_header != NULL) {
    table_hdrfree(input_header);
  }
  if (input_handle != NULL) {
    Close(input_handle);
  }
  return(EXIT_SUCCESS);
}

