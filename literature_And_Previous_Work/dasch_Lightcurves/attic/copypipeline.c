// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* copypipeline.c
 *
 * Given a list of files, synchronize these files between the local /dasch/Pipeline directory and
 * a directory specified by the environment variable DASCH_PIPELINE_DISK.  However, if the
 * local system is the system specified in DASCH_PIPELINE_SERVER, then no copy operaton occurs
 *
 * Format of the copyfiles.txt template file
 * Field 1: doNotDelete      0 = o.k. to delete the original, 1 = never delete the original
 * Field 2: solution number: 0 = current solution; 1 = next solution
 * Field 3: suffix type:     0:  _01ww, _01r90ww, _01r180ww, or _01r270ww  takes solutionString (_sn) and qualifiers (_kepler or _apass)
 *                           1:  _16ww  (only in astrometry and headers, takes solutionString and no rotation)
 *                           2:  _16 and _01, no suffix in astrometry directory only
 *
 *  gcc -ggdb -O0  -I/usr/include/mysql  -I/dasch/install/include - -L /dasch/install/lib   copypipeline.c  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread -lm -o copypipeline
 *
 *  copypipeline  /home/scanner/Pipeline/dell.list -s 0 -f -q apass
 *
 * Sep 12, 2011 Edward J. Los - Initial version
 * Sep 16, 2011 Edward J. Los - Add a flag to prevent deletion of Sextractor output files
 * Dec 23, 2013 Edward J. Los - correct the copy count and add more diagnostic information.
 * Jan 19, 2019 Edward J. Los - add diagnostic mode to check file list
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <dirent.h>
#include "mysql.h"
#include "pipelineutils.h"



#define MAX_FILENAME 512
#define MAX_BUFFER 512
#define MAX_ALTLOGFILES 3
#define MAX_LIST_STRING 50
#define MAX_FORMAT_STRING 50
#define MAX_TEMPLATES 200

typedef struct _filetemplate {
  int doNotDelete;   /* Do not delete if nonzero */
  int solutionNumber; /* Solution number */
  int suffixType;     /* Suffix type */
  char directory[MAX_FORMAT_STRING];
  char formatString[MAX_FORMAT_STRING];
} FILETEMPLATE,*PFILETEMPLATE;


int CopyFile(PFILETEMPLATE pTemplate,char *sourcename,char *targetname,int forwardMode,int *pCopyCount,int *pDeleteCount,int *pWorkDone,int verbose,off_t* pTotalBytes)
{
  int sourceresult;
  int targetresult;
  int copyresult;
  struct stat sourcestats;
  struct stat targetstats;
  struct stat copystats;
  int result;
  char cmdStr[(2*MAX_FILENAME)+5];
#if 0
  printf("Source: %s Target: %s\n",sourcename,targetname);
#endif
  sourceresult = stat(sourcename,&sourcestats);
  targetresult = stat(targetname,&targetstats);

  if (forwardMode) {
    /* Source does not exist.  Nothing to copy */
    if (sourceresult != 0) {
      return(0);
    }
    /* Target exists don't copy if it is the same file */
    if (targetresult == 0) {
      if ((sourcestats.st_dev == targetstats.st_dev) &&
          (sourcestats.st_ino == targetstats.st_ino)) {
        printf("ERROR: %s and %s are the same file\n",sourcename,targetname);
        exit(-1);
      }

    }
    /* It is ok to copy the file */
    sprintf(cmdStr,"cp %s %s\n",sourcename,targetname);
#if 0
		printf("%s",cmdStr);
#endif
    *pWorkDone = 1;
    *pTotalBytes += sourcestats.st_size;
    result = system(cmdStr);
    if (verbose) {
      printf("Result is %d 0x%x for %s\n",result,result,cmdStr);
    }
    copyresult = stat(targetname,&copystats);
    if (copyresult != 0) {
      printf("ERROR copying to %s\n",targetname);
      exit(-1);
    }
    (*pCopyCount)++;
    if (copystats.st_size != sourcestats.st_size) {
      printf("ERROR: wrong file sizes for %s (%d) and %s (%d)\n",sourcename,sourcestats.st_size,targetname,targetstats.st_size);
    }
    
  } else { /* reverse mode */
#if 0
    if (pTemplate->doNotDelete != 0) {
      printf(".");
    }
#endif
    if (sourceresult != 0) {
      if (targetresult == 0) {
        if (pTemplate->doNotDelete == 0) { 
          /* The source file no longer exists, but the target file does.  Delete the target file */
          *pWorkDone = 1;
          result = unlink(targetname);
          if (result != 0) {
            printf("ERROR: could not delete %s\n",targetname);
          }
          if (verbose) {
            printf("Delete file %s\n",targetname);
          }
          (*pDeleteCount)++;
          result = -1;
        }
      }
    } else {
      /* Target exists don't copy if it is the same file */
      if (targetresult == 0) {
        if ((sourcestats.st_dev == targetstats.st_dev) &&
            (sourcestats.st_ino == targetstats.st_ino)) {
          printf("ERROR: %s and %s are the same file\n",sourcename,targetname);
          exit(-1);
        }
      }
      /* It is ok to copy the file */
      sprintf(cmdStr,"cp %s %s\n",sourcename,targetname);
#if 0
			printf("%s",cmdStr);
#endif
      *pWorkDone = 1;
      *pTotalBytes += sourcestats.st_size;
      result = system(cmdStr);
      if (verbose) {
        printf("Result is %d 0x%x for %s\n",result,result,cmdStr);
      }
      copyresult = stat(targetname,&copystats);
      if (copyresult != 0) {
        printf("ERROR copying to %s\n",targetname);
        exit(-1);
      }
      (*pCopyCount)++;
      if (copystats.st_size != sourcestats.st_size) {
        printf("ERROR: wrong file sizes for %s (%d) and %s (%d)\n",sourcename,sourcestats.st_size,targetname,targetstats.st_size);
      } else {
        /* It is o.k. to delete the source file */
        *pWorkDone = 1;
        result = unlink(sourcename);
        if (result != 0) {
          printf("ERROR: could not delete %s\n",sourcename);
        }
      }
    }
  }

  return(0);
}


int main(int argc,char *argv[])
{
  int errorFlag = 0;
  char *argstr;
  char cmdchar;
  int nvals;
  char *argstr2;
  int result;
  time_t startTime = 0;
  time_t curTime;
  time_t staleTime;
  char listname[MAX_FILENAME]; /* List of plates */
  int solutionNumber = 0; /* current WCS solution */
  int verbose = 0;
  int forwardMode = 0;
  int reverseMode = 0;
  char *localhost;
  char *daschPipelineDisk;
  char *daschPipelineServer;
  int pipelineDisk = 0;
  char *catalogallDirectory;
  char *catalogbinDirectory;
  char *binoutputDirectory;
  char *ingestDirectory;
  char *astrometryDirectory;
  char *headersDirectory;
  char *matchDirectory;
  char *scampimagesDirectory;
  char *scriptDirectory;
  FILE * listFile;
  char inLine[MAX_BUFFER];
  char copyBuffer[MAX_BUFFER];
  char *inBuffer;
  int lineLen;
  int nlines = 0;
  int maxString = 0;
  char *listStrings = NULL;
  int curline = 0;
  int catalogNumber = 0;
  char source[MAX_BUFFER];
  char catalogString[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  PFILETEMPLATE templateTable = NULL;
  PFILETEMPLATE pTemplate;
  char templateName[MAX_BUFFER];
  FILE *templateFile;
  int templateCount = 0;
  char *charPtr;
  char *charPtr2;
  char filebasename[MAX_BUFFER];
  char timestampname[MAX_BUFFER];
  char lastfilename[MAX_BUFFER];
  int timestampresult; 
  struct stat timestampstats;
  FILE *timestamphandle;
  char rootname[MAX_BUFFER];
  char rootfilename[MAX_FILENAME];
  char sourcename[MAX_FILENAME];
  char targetname[MAX_FILENAME];
  char localDirectory[MAX_BUFFER];
  char saveDirectory[MAX_BUFFER];
  int fileIndex;
  int templateIndex;
  int iteration;
  char solutionString[MAX_BUFFER];
  int testcount = 0;
  int copycount = 0;
  int deletecount = 0;
  int workDone;
  off_t totalBytes = 0;
  int diagnosticMode = 0;
  int lineCounter = 0;

  source[0] = 0;
  catalogString[0] = 0;
  solutionString[0] = 0;
  qualifier[0] = 0;
  lastfilename[0] = 0;

  listname[0] = 0;
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      /* This must be the list of plates */
      if (strlen(listname) > 0) {
        errorFlag = 1;
        printf("ERROR: list %s is being overwritten by %s\n",listname,argstr);
      } else {
        strncpy(listname,argstr,MAX_FILENAME-1);
      }
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {

        case 'v': /* Verbose mode */
        case 'V':
          verbose = 1;
          break;

        case 'f': /* Verbose mode */
        case 'F':
          forwardMode = 1;
          break;

        case 'q': /* Catalog and file name qualifier */
        case 'Q':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            catalogNumber = GetCatalogNumber(*++argv);
            if (catalogNumber < 0) {
              printf("ERROR: Illegal catalog name %s\n",*argv);
              errorFlag = 1;
            } else {
              strcpy(source,*argv);
              if (catalogNumber > 0) {
                sprintf(catalogString,"%d",catalogNumber);
                strcpy(qualifier,"_");
                strcat(qualifier,*argv);
              }
            }
          }
          break;

        case 'r': /* Verbose mode */
        case 'R':
          reverseMode = 1;
          break;

        case 'd': /* Diagnostic mode */
        case 'D':
          diagnosticMode = 1;
          break;

        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }



  /* Validate arguments */

  if (listname[0]  == 0)  {
    printf("ERROR: the list of plates is not specified\n");
    errorFlag = 1;
  }
  if ((reverseMode ^ forwardMode) != 1) {
    printf("ERROR: either -f %d or -r %d must be specified\n",forwardMode,reverseMode);
    errorFlag = 1;
  }

  localhost = getenv("HOST");
  if (localhost == NULL) {
    printf("ERROR: HOST is not defined\n");
    errorFlag = 1;
  }
  
  daschPipelineDisk = getenv("DASCH_PIPELINE_DISK");
  if (daschPipelineDisk == NULL) {
    printf("ERROR: DASCH_PIPELINE_DISK is not defined\n");
    errorFlag = 1;
  } else {
    nvals = sscanf(daschPipelineDisk,"%d",&pipelineDisk);
    if (nvals != 1) {
      printf("ERROR: problem decoding %s as an integer\n",daschPipelineDisk);
      errorFlag = 1;
    }
  }
  
  daschPipelineServer = getenv("DASCH_PIPELINE_SERVER");
  if (daschPipelineServer == NULL) {
    printf("ERROR: DASCH_PIPELINE_SERVER is not defined\n");
    errorFlag = 1;
  }
  
  catalogallDirectory = getenv("DASCH_CATALOGALL");
  if (catalogallDirectory == NULL) {
    printf("ERROR: DASCH_CATALOGALL is not defined\n");
    errorFlag = 1;
  }
  
  catalogbinDirectory = getenv("DASCH_CATALOGBIN");
  if (catalogbinDirectory == NULL) {
    printf("ERROR: DASCH_CATALOGBIN is not defined\n");
    errorFlag = 1;
  }
  
  binoutputDirectory = getenv("DASCH_BINOUTPUT");
  if (binoutputDirectory == NULL) {
    printf("ERROR: DASCH_BINOUTPUT is not defined\n");
    errorFlag = 1;
  }
  
  ingestDirectory = getenv("DASCH_INGEST");
  if (ingestDirectory == NULL) {
    printf("ERROR: DASCH_INGEST is not defined\n");
    errorFlag = 1;
  }
  
  astrometryDirectory = getenv("DASCH_ASTROMETRY");
  if (astrometryDirectory == NULL) {
    printf("ERROR: DASCH_ASTROMETRY is not defined\n");
    errorFlag = 1;
  }
  
  headersDirectory = getenv("DASCH_HEADERS");
  if (headersDirectory == NULL) {
    printf("ERROR: DASCH_HEADERS is not defined\n");
    errorFlag = 1;
  }
  
  matchDirectory = getenv("DASCH_MATCH");
  if (matchDirectory == NULL) {
    printf("ERROR: DASCH_MATCH is not defined\n");
    errorFlag = 1;
  }
  
  scampimagesDirectory = getenv("DASCH_SCAMPIMAGES");
  if (scampimagesDirectory == NULL) {
    printf("ERROR: DASCH_SCAMPIMAGES is not defined\n");
    errorFlag = 1;
  }
  
  scriptDirectory = getenv("DASCH_SCRIPTS");
  if (scriptDirectory == NULL) {
    printf("ERROR: DASCH_SCRIPTS is not defined\n");
    errorFlag = 1;
  }
  
  




  if (errorFlag) {
    printf("Usage: copypipeline <plate list> -s  {-f|-r} \n");
    printf("                   where -d creates the plate list from the MySql database (first pass only)\n");
    printf("                         -v verbose mode\n");
    printf("                         -f forward mode: copy from primary storage to local system\n");
    printf("                         -r reverse mode: copy from local system to primary storage\n");
    printf("                         -q Catalog qualifier (apass or kepler)\n");
    exit(-1);
  }
  
  if (diagnosticMode != 0) {
    printf("ERROR: copypipeline diagnosticMode is enabled\n");
  } else {

    if (strstr(localhost,daschPipelineServer) != NULL) {
      printf("copypipeline exiting because %s matches %s\n",daschPipelineServer,localhost);
      return(0);
    }
  }
  /* Attempt to open the list */
  listFile = fopen(listname,"rt");
  if (listFile == NULL) {
    printf("Could not open file %s\n",listname);
    errorFlag = 1;
  } else {


    /* Count the number of records in the file */
    while (1) {
      inBuffer = fgets(inLine,MAX_BUFFER,listFile);
      if (inBuffer == NULL) {
        break;
      }
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
      if (lineLen > 5) {
        nlines++;
        if (maxString < lineLen) {
          maxString = lineLen;
        }
      }
  

    }
    fclose(listFile);
    if (nlines == 0) {
      printf("Error: no lines in the list file\n");
      errorFlag = 1;
    } else {
      maxString += 5;
      listStrings = calloc((nlines+1)*maxString,sizeof(char));
      if (listStrings == NULL) { 
        printf("Failed to allocate listStrings %x \n",listStrings);
        errorFlag = 1;
      } else {
        strcpy(&listStrings[nlines*maxString],"UNKNOWN");



        listFile = fopen(listname,"rt");
        if (listFile == NULL) {
          printf("Could not open file %s\n",listname);
          errorFlag = 1;
        } else {

          /* Read the file records */
          while (1) {
            inBuffer = fgets(inLine,MAX_BUFFER,listFile);
            if (inBuffer == NULL) {
              break;
            }
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
            if (lineLen > 5) {
              strcpy(&listStrings[curline*maxString],inBuffer);
              curline++;
            }
  
          }
          fclose(listFile);
        }
     
      }
    }
  }

  /* Now open the template file */
  templateTable = (PFILETEMPLATE)calloc(MAX_TEMPLATES,sizeof(FILETEMPLATE));
  if (templateTable == NULL) {
    printf("ERROR: failed to allocate the template table of size %d\n",MAX_TEMPLATES*sizeof(FILETEMPLATE));
    errorFlag = 1;
  } else {
    strcpy(templateName,scriptDirectory);
    strcat(templateName,"/copypipeline.txt");
    templateFile = fopen(templateName,"rt");
    if (templateFile == NULL) {
      printf("Could not open file %s\n",templateName);
      errorFlag = 1;
    } else {
      while (1) {
        pTemplate = &templateTable[templateCount];
        inBuffer = fgets(inLine,MAX_BUFFER,templateFile);
        if (inBuffer == NULL) {
          break;
        }
        lineCounter++;
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
        strcpy(copyBuffer,inBuffer);
        charPtr = strstr(inBuffer," ");
        if (charPtr == NULL) {
          printf("ERROR in template file at line %s\n",inBuffer);
          continue;
        }
        charPtr++;
        charPtr = strstr(charPtr," ");
        if (charPtr == NULL) {
          printf("ERROR in template file at line %s\n",inBuffer);
          continue;
        }
        charPtr++;
        charPtr = strstr(charPtr," ");
        if (charPtr == NULL) {
          printf("ERROR in template file at line %s\n",inBuffer);
          continue;
        }
        charPtr++;
        nvals = sscanf(inBuffer,"%d %d %d",&pTemplate->doNotDelete,&pTemplate->solutionNumber,&pTemplate->suffixType);
        if (nvals != 3) {
          printf("ERROR in template file at line %s\n",inBuffer);
          continue;
        }
        charPtr2 = strstr(charPtr," ");
        if (charPtr2 == NULL) {
          printf("ERROR in template file at line %s\n",inBuffer);
          continue;
        }
        *charPtr2 = 0;
        charPtr2++;

        if (strlen(charPtr) > (MAX_FORMAT_STRING-1)) {
          printf("ERROR: format string too long in %s\n",copyBuffer);
          continue;
        }
        strcpy(pTemplate->directory,charPtr);

        if (strlen(charPtr2) > (MAX_FORMAT_STRING-1)) {
          printf("ERROR: format string too long in %s\n",copyBuffer);
          continue;
        }
        strcpy(pTemplate->formatString,charPtr2);


        if (templateCount < MAX_TEMPLATES) {
          templateCount++;
        } else {
          printf("ERROR: template file too long\n");
        }

      }
      fclose(templateFile);
    }

  }
                   

  if (diagnosticMode != 0) {
    printf("templates %d from lines %d in %s\n",templateCount,lineCounter,templateName);
  }


  printf("copypipeline of %s %s, list file %s  entries %d, max string %d forward %d templates %d\n",
         __DATE__,__TIME__,listname,nlines,maxString,forwardMode,templateCount);


  time(&startTime);
  for (fileIndex = 0; fileIndex < nlines; fileIndex++) {
    strcpy(filebasename,&listStrings[fileIndex*maxString]);
    solutionNumber = 0;
    sprintf(timestampname,"%s/%s%s_copypipeline.txt",catalogbinDirectory,filebasename,qualifier);
    if (forwardMode) {
      /* Create a timestamp file to avoid two reverse copy operations without a forward operation */
      timestamphandle = fopen(timestampname,"wt");
      if (timestamphandle == NULL) {
        printf("ERROR: failed to open %s\n",timestampname);
        exit(-1);
      }
      fprintf(timestamphandle,"%s %d\n",timestampname,startTime);
      fclose(timestamphandle);
    } else {
      /* Look for a timestamp that says that this step was preceeded by a copy operation */
      timestampresult = stat(timestampname,&timestampstats);
      if (diagnosticMode == 0) {
        if (timestampresult != 0) {
          printf("ERROR: reverse copy without forward copy, failed to find %s\n",timestampname);
          exit(-1);
        }
        staleTime = startTime - timestampstats.st_mtim.tv_sec;
        if (staleTime > (3600*24)) {
          printf("ERROR: stale time by %d seconds in %s\n",staleTime,timestampname);
          exit(-1);
        }
        /* Safe to preceed.  Get rid of the file */
        timestampresult = unlink(timestampname);
        if (timestampresult != 0) {
          printf("ERROR: failed to remove %s\n",timestampname);
        }
      }
    }
    while (1) {  /* Solution Number Loop */
      workDone = 0;
      strcpy(lastfilename,filebasename);
      if (solutionNumber == 0) {
        solutionString[0] = 0;
      } else {
        sprintf(solutionString,"_s%d",solutionNumber);
      } 


      for (templateIndex = 0; templateIndex < templateCount; templateIndex++) {
        pTemplate = &templateTable[templateIndex];
        sprintf(saveDirectory,"/dasch/raid%03d/Pipeline/%s/",pipelineDisk,pTemplate->directory);
        if (strstr(pTemplate->directory,"astrometry") != NULL) {
          printf("ERROR: copypipeline astrometry has not been debugged\n");
          exit(-1);
          strcpy(localDirectory,astrometryDirectory);
        } else if (strstr(pTemplate->directory,"bin9") != NULL) {
          strcpy(localDirectory,binoutputDirectory);
        } else if (strstr(pTemplate->directory,"catalog9bin") != NULL) {
          strcpy(localDirectory,catalogbinDirectory);
        } else if (strstr(pTemplate->directory,"catalogall") != NULL) {
          strcpy(localDirectory,catalogallDirectory);
        } else if (strstr(pTemplate->directory,"headers") != NULL) {
          printf("ERROR: copypipeline headers has not been debugged\n");
          exit(-1);
          strcpy(localDirectory,headersDirectory);
        } else if (strstr(pTemplate->directory,"ingest") != NULL) {
          strcpy(localDirectory,ingestDirectory);
        } else if (strstr(pTemplate->directory,"match") != NULL) {
          strcpy(localDirectory,matchDirectory);
        } else if (strstr(pTemplate->directory,"scampimages") != NULL) {
          strcpy(localDirectory,scampimagesDirectory);
        } else {
          printf("ERROR: unrecognized directory %s\n",pTemplate->directory);
          exit(-1);
        }
        strcat(localDirectory,"/");
        if (pTemplate->solutionNumber == 0) {
          if (solutionNumber == 0) {
            solutionString[0] = 0;
          } else {
            sprintf(solutionString,"_s%d",solutionNumber);
          } 
        } else {
          sprintf(solutionString,"_s%d",solutionNumber+1);          
        }
        for (iteration = 0; iteration < 2; iteration++) {
          if ((iteration == 1) && (pTemplate->suffixType < 2)) {
            break;
          }
          switch (pTemplate->suffixType) {
          case 0:
            strcpy(rootname,filebasename);
            sprintf(rootfilename,pTemplate->formatString,rootname,solutionString,qualifier);
            break;
          case 1:
            printf("ERROR: copypipeline case 1 has not been debugged\n");
            exit(-1);

            strcpy(rootname,filebasename);
            charPtr = strstr(rootname,"_01ww");
            if (charPtr != NULL) {
              strcpy(charPtr,"_16ww");
            } else {
              charPtr = strstr(rootname,"_01r");
              if (charPtr != NULL) {
                strcpy(charPtr,"_16ww");
              } else {
                printf("ERROR: illegal file base name %s\n",rootname);
                break;
              }
            }
            sprintf(rootfilename,pTemplate->formatString,rootname,solutionString);
            break;
          case 2:
            printf("ERROR: copypipeline case 1 has not been debugged\n");
            exit(-1);
         
            strcpy(rootname,filebasename);
            charPtr = strstr(rootname,"_01ww");
            if (charPtr != NULL) {
              if (iteration == 0) {
                strcpy(charPtr,"_16");
              } else {
                strcpy(charPtr,"_01");
              }
            } else {
              charPtr = strstr(rootname,"_01r");
              if (charPtr != NULL) {
                if (iteration == 0) {
                  strcpy(charPtr,"_16");
                } else {
                  strcpy(charPtr,"_01");
                }
              } else {
                printf("ERROR: illegal file base name %s\n",rootname);
                break;
              }
            }
            sprintf(rootfilename,pTemplate->formatString,rootname);
 

            break;
          default:
            printf("ERROR: Illegal suffixType %d\n",pTemplate->suffixType);
            exit(-1);
          }
          if (forwardMode) {
            strcpy(sourcename,saveDirectory);
            strcat(sourcename,rootfilename);
            strcpy(targetname,localDirectory);
            strcat(targetname,rootfilename);
          } else {
            strcpy(sourcename,localDirectory);
            strcat(sourcename,rootfilename);
            strcpy(targetname,saveDirectory);
            strcat(targetname,rootfilename);
          } 
          if (diagnosticMode == 0) {
            CopyFile(pTemplate,sourcename,targetname,forwardMode,&copycount,&deletecount,&workDone,verbose,&totalBytes);
          } else {
            printf("DIAG: target %s\n",targetname);
          }
#if 0
          printf("Count %3d Check: %s and %s\n",copycount,sourcename,targetname);
#endif
          testcount++;

        }
      }
      if (diagnosticMode == 0) {
        if (workDone == 0) {
          /* No work done, so we can exit */
          break;
        }
      } else {
        if (solutionNumber >= 1) {
          break;
        }
      }
      solutionNumber++;
    }
  }



  time(&curTime);
  curTime -= startTime;
  printf("copypipeline completed in %d seconds %d forward %d checks %d copies %d deleted for %s s%d%s  %lld bytes copied\n",curTime,forwardMode,testcount,copycount,deletecount,lastfilename,solutionNumber,source,totalBytes);
  if (listStrings != NULL) {
    free(listStrings);
  }
  if (templateTable != NULL) {
    free(templateTable);
  }

  return(0);
}
