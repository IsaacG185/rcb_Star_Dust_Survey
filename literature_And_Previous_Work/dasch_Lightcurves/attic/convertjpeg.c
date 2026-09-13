// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* convertjpeg.c
 *
 * Convert a list of JPEG images to lower resolution for CameraPC verification
 *
 *
 *  gcc -ggdb -O0  -I/usr/include/mysql  -I/dasch/install/include  -L/dasch/install/lib   convertjpeg.c  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread -lm -o convertjpeg
 *
 *  
 *   The image file or input list should have the full path (raid023jpeg.list) the root directory replaces ExposureData (/dasch/raid020/junk/JPEG).  The "i" option is a single file, e.g. 
 *      /dasch/raid023/ExposureData/JacketJpegBackup/b/b31104_j.jpg
 *   The output file is reduced in size with the following commands:
 *     jpegtopnm <input file>  | ppmtopgm  | pamscale -quiet -ysize=2456 -nomix | cjpeg -quality 100   >  <output file>
 *  the exiftool is used to write DateTimeOriginal from the old file, the Copyright, and the original Location as /JacketJpegBackup/b/b31104_j.jpg

       convertjpeg  -i <input file>  -t target directory -l <image list>
       convertjpeg -l /home/scanner/junk/JPEG/nashua_data.list -t /home/scanner/junk/JPEG/junk
       convertjpeg -l /home/scanner/junk/JPEG/nashua_finished.txt -t /home/scanner/junk/JPEG/junk

 *
 * Jan 21, 2019 Edward J. Los - Initial version
 * Jan 24, 2019 Edward J. Los - Fix names for single character series.
 * Jan 28, 2019 Edward J. Los - Add a comment at the end of the list file to be copied to the final image
 * Jan 28, 2019 Edward J. Los - add '-b' option to generate backup sets
 * Jan 31, 2019 Edward J. Los - handle CameraPC JPG files
 * Feb  4, 2019 Edward J. Los - correct single file case
 * Feb  6, 2019 Edward J. Los - correct sky coordinates for fits files
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


#define MAX_BUFFER 1024
#define MAX_DATE_STR 30
#define MAX_CMD_STR ((2*MAX_BUFFER)+5)
#define MAX_ALTLOGFILES 3
#define MAX_LIST_STRING 50
#define MAX_FORMAT_STRING 50
#define MAX_TEMPLATES 200
#define MAX_TARGET 6
#define MAX_TARGET_COUNT 150
#define TARGET_ALLOC_INCREMENT ((2*MAX_BUFFER)+5)

int numTarget = 0;
int numTarball = 0;
char targetstring[MAX_TARGET][MAX_BUFFER];
char targetroot[MAX_TARGET][MAX_BUFFER];
int targetCount[MAX_TARGET];
int tarballAlloc[MAX_TARGET];
char *tarballString[MAX_TARGET];
char *tempString;


int invokeCmd(char *cmdStr,int *timeSpent) {
  int result;
  time_t beginTime;
  time_t curTime;
  time(&beginTime);
  *timeSpent = 0;
  result = system(cmdStr);
  time(&curTime);
  beginTime = curTime - beginTime;
  *timeSpent = (int)beginTime;
  if (WIFSIGNALED(result) &&
      (WTERMSIG(result) == SIGINT || WTERMSIG(result) == SIGQUIT)) {
    printf("Terminated with signal %d at %d seconds\n",WTERMSIG(result),(int)curTime);
    exit(-1);
  }
#if 0
  if (WIFEXITED(result)) {
    printf("Exited with status %d at %d sec delta %d sec s\n",WEXITSTATUS(result),(int)curTime,(int)beginTime);
  } else {
    printf("Unknown result %d at %d sec delta %d step %d\n",result,(int)curTime,(int)beginTime,processingStep);
  }
  printf("result: %d\n",result);
#endif

  return(result);
}


int main(int argc,char *argv[])
{
  int errorFlag = 0;
  char *argstr;
  char cmdchar;
  int fileIndex;
  time_t curTime;
  time_t startTime = 0;
  char listname[MAX_BUFFER]; /* List of images */
  char singlename[MAX_BUFFER]; /* Single image */
  char targetbasedirectory[MAX_BUFFER]; /* Target directory */
  char filebasename[MAX_BUFFER];
  char targetname[MAX_BUFFER];
  char datename[MAX_BUFFER];
  char datestring[MAX_DATE_STR];
  char locationstring[MAX_BUFFER];
  char tempbuffer[MAX_BUFFER];
  char rootname[MAX_BUFFER];
  char origrootname[MAX_BUFFER];
  int verbose = 0;
  int diagnosticMode = 0;
  int backupMode = 0;
  int JPGMode = 0;
  FILE * listHandle;
  FILE * dateHandle;
  char inLine[MAX_BUFFER];
  int linelen;
  int nlines = 0;
  int maxString = 0;
  char *listStrings = NULL;
  int lineCounter = 0;
  int curline = 0;
  char *inBuffer;
  struct stat filestats;
  int result;
  char *charPtr;
  char *charPtr2;
  char *charPtr3;
  char *charPtr4;
  char *charPtr5;
  char *charPtr6;
  char *charPtr7;
  char *charPtr8;
  char *charPtr9;
  char *commentPtr;
  char blankString[] = "BLANK";
  char cmdStr[MAX_CMD_STR];
  int timeSpent;
  off_t inputBytes = 0;
  off_t outputBytes = 0;
  int fitsFlag;
  int targetIndex;

  memset(targetCount,0,sizeof(targetCount));
  memset(tarballAlloc,0,sizeof(tarballAlloc));
  memset(tarballString,0,sizeof(tarballString));
  singlename[0] = 0;
  listname[0] = 0;
  targetbasedirectory[0] = 0;
  /* Loop through the arguments */
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: stray argument %s\n",argstr);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 'l':  /* list of file names */
        case 'L':          
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(listname,*++argv,MAX_BUFFER-2);
						if (strlen(*argv) >= MAX_BUFFER-2) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
						}
					}
					break;

        case 'i':  /* Single file entry */
        case 'I':          
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(singlename,*++argv,MAX_BUFFER-2);
						if (strlen(*argv) >= MAX_BUFFER-2) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
						}
					}
					break;

        case 't':  /* target directory */
        case 'T':          
					argc--;
					if (argc < 1) {
						fprintf(stderr,"ERROR: Insufficient arguments for -%c\n",cmdchar);
						errorFlag = 1;
					} else {
						strncpy(targetbasedirectory,*++argv,MAX_BUFFER-2);
						if (strlen(*argv) >= MAX_BUFFER-2) {
							fprintf(stderr,"ERROR: MAX_BUFFER too small for %s\n",*argv);
						}
            linelen = strlen(targetbasedirectory);
            if (targetbasedirectory[linelen-1] != '/') {
              strcat(targetbasedirectory,"/");
            }
					}
					break;

        case 'v': /* Verbose mode */
        case 'V':
          verbose = 1;
          break;


        case 'd': /* Diagnostic mode */
        case 'D':
          diagnosticMode = 1;
          break;

        case 'b': /* generate backup tarballs */
        case 'B':
          backupMode = 1;
          break;

        default:
          printf("ERROR: * illegal command -%c-",cmdchar);
          errorFlag = 1;

        }
        
      }

    }
  }



  /* Validate arguments */

  if ((listname[0]  == 0) && (singlename[0] == 0)) {
    printf("ERROR: the list of plates (-l) and a single name (-i) is not specified\n");
    errorFlag = 1;
  }
  if ((listname[0] != 0) && (singlename[0] != 0)) {
    printf("ERROR: both the list of plates (-l) and a single name (-i) is  specified\n");
    errorFlag = 1;
  }
  if (targetbasedirectory[0] == 0) {
    printf("ERROR: target directory name not specified\n");
    errorFlag = 1;
  }
 
  if (errorFlag) {
    printf("Usage: convertjpeg {-i <single file>| -l <image list>  \n");
    printf("                   where -i <jpeg image to be converted>\n");
    printf("                         -l <list of images to be converted>\n");
    printf("                         -t <root directory to copy the file to>\n");
    printf("                         -d disgnostic mode\n");
    printf("                         -v verbose mode\n");
    exit(-1);
  }

  
  if (diagnosticMode != 0) {
    printf("ERROR: convertjpeg diagnosticMode is enabled\n");
  }
  /* Attempt to open the list */
  if (listname[0] != 0) {
    listHandle = fopen(listname,"rt");
    if (listHandle == NULL) {
      printf("ERROR: Could not open file %s\n",listname);
      errorFlag = 1;
    } else {


      /* Count the number of records in the file */
      while (1) {
        inBuffer = fgets(inLine,MAX_BUFFER,listHandle);
        if (inBuffer == NULL) {
          break;
        }
        linelen = strlen(inBuffer);
        /* Trim off the carriage return */
        if (inBuffer[linelen-1] == 10) {
          inBuffer[linelen-1] = 0;
          linelen--;
        }
        /* Trim off the line feed */
        if (inBuffer[linelen-1] == 13) {
          inBuffer[linelen-1] = 0;
          linelen--;
        }
        if (linelen > 5) {
          nlines++;
          if (maxString < linelen) {
            maxString = linelen;
          }
        }
  

      }
      fclose(listHandle);
      if (maxString > (MAX_BUFFER/3)) {
        printf("ERROR maxString %d greater than MAX_BUFFER/3 %d\n",maxString,MAX_BUFFER/3);
        exit(-1);
      }
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



          listHandle = fopen(listname,"rt");
          if (listHandle == NULL) {
            printf("Could not open file %s\n",listname);
            errorFlag = 1;
          } else {

            /* Read the file records */
            while (1) {
              inBuffer = fgets(inLine,MAX_BUFFER,listHandle);
              if (inBuffer == NULL) {
                break;
              }
              lineCounter++;
              linelen = strlen(inBuffer);
              /* Trim off the carriage return */
              if (inBuffer[linelen-1] == 10) {
                inBuffer[linelen-1] = 0;
                linelen--;
              }
              /* Trim off the line feed */
              if (inBuffer[linelen-1] == 13) {
                inBuffer[linelen-1] = 0;
                linelen--;
              }
              if ((strstr(inBuffer,"/ExposureData/") == NULL) &&
                  (strstr(inBuffer,".JPG") == NULL)) { 
                printf("ERROR: '/ExposureData/' is not in line %d of %s\n",lineCounter,listname);
                continue;
                  
              }

              if (linelen > 5) {
                strcpy(&listStrings[curline*maxString],inBuffer);
                curline++;
                if (curline > nlines) {
                  printf("ERROR: nlines exceeded\n");
                  exit(-1);
                }
              }
  
            }
            fclose(listHandle);
            nlines = curline;
          }
     
        }
      }
    }
  } else { 
    /* Single file case */
    maxString = strlen(singlename)+5;
    nlines = 1;
    listStrings = calloc((nlines+1)*maxString,sizeof(char));
    if (listStrings == NULL) { 
      printf("Failed to allocate listStrings %x \n",listStrings);
      exit(-1);
    } else {
      strcpy(listStrings,singlename);
    }

  }



  printf("convertjpeg of %s %s, entries %d, max string %d\nlist file '%s'\n single name '%s'\n targetbasedirectory '%s'\n",
         __DATE__,__TIME__,nlines,maxString,listname,singlename,targetbasedirectory);


  time(&startTime);
  for (fileIndex = 0; fileIndex < nlines; fileIndex++) {
    strcpy(filebasename,&listStrings[fileIndex*maxString]);
    commentPtr = strchr(filebasename,' ');
    if (commentPtr != NULL) {
      *commentPtr = 0;
      commentPtr++;
      if (strlen(commentPtr) == 0) {
        commentPtr = NULL;
      }
        
    }

    if (strstr(filebasename,".fit") != NULL) {
      fitsFlag = 1;
    } else {
      fitsFlag = 0;
    }

    if (backupMode == 0) {
      result = stat(filebasename,&filestats);
      if (result != 0) {
        printf("ERROR: source file %s does not exist\n",filebasename);
        continue;
      }
      inputBytes += filestats.st_size;
    }
    strcpy(targetname,targetbasedirectory);
    strcpy(tempbuffer,filebasename);
    charPtr = strstr(tempbuffer,".JPG");
    if (charPtr != NULL) {
      JPGMode = 1;
      charPtr = strstr(tempbuffer,"CameraPC");
      strcpy(locationstring,charPtr);
    } else {
      JPGMode = 0;
      charPtr = strstr(tempbuffer,"/ExposureData/");
      if (charPtr == NULL) {
        printf("ERROR: source file %s does not contain '/ExposureData/'\n");
        continue;
      }
      charPtr += strlen("/ExposureData/"); /* charPtr now points to the subdirectory */
      locationstring[0] = '/';
      strcpy(&locationstring[1],charPtr);
    }
    /* Create a directory */
    strcat(targetname,charPtr);
    charPtr2 = strrchr(targetname,'/');
    if (charPtr2 == NULL) {
      printf("ERROR no '/' in target directory %s\n",targetname);
      continue;
    }
    *charPtr2 = 0;
    charPtr2++; /* charPtr2 now points to the final filename */
    strcpy(rootname,charPtr2);
    /* Now step back again to skip over series */
    charPtr3 = strrchr(targetname,'/');
    if (charPtr3 == NULL) {
      printf("ERROR no '/' in target directory %s\n",targetname);
      continue;
    }
    *charPtr3 = 0;
    if ((fitsFlag != 0) || (JPGMode != 0)) {
      /* For fits files, go back another level */
      charPtr3 = strrchr(targetname,'/');
      if (charPtr3 == NULL) {
        printf("ERROR no '/' in target directory %s\n",targetname);
        continue;
      }
      *charPtr3 = 0;
    }
    

    if (backupMode == 0) {
      result = stat(targetname,&filestats);
      if (result != 0) {
        /* Need to create the directory */
        sprintf(cmdStr,"mkdir -p %s\n",targetname);
        result = invokeCmd(cmdStr,&timeSpent);
        if (result != 0) {
          printf("ERROR: result %d for %s\n",result,cmdStr);
          continue;
        }
      }
    } else {
      for (targetIndex = 0; targetIndex < numTarget; targetIndex++) {
        if (strcmp(targetname,&(targetstring[targetIndex][0])) == 0) {
          break;
        }
      }
      if (targetIndex >= numTarget) {
        if (targetIndex >= MAX_TARGET) {
          printf("ERROR: line %d MAX_TARGET exceeded\n",__LINE__);
          exit(-1);
        }
        strcpy(targetstring[targetIndex],targetname);
        charPtr7 = &targetstring[targetIndex][0]+strlen(targetbasedirectory);
        targetroot[targetIndex][0] = '.';
        targetroot[targetIndex][1] = '/';
        strcat(&targetroot[targetIndex][2],charPtr7);
        numTarget++;
      }
        
    }
    *charPtr3 = '/';
    charPtr3++;
    *charPtr3 = 0;
    strcat(targetname,rootname);
    if (fitsFlag) {
      strcpy(origrootname,rootname);
      charPtr5 = strrchr(rootname,'_');
      if (charPtr5 == NULL) {
        printf("ERROR: failed to fine '_' in %s\n",targetname);
        continue;
      }
      *charPtr5 = 0;
      charPtr4 = strstr(targetname,".fit");
      if (charPtr4 == NULL) {
        printf("ERROR: failed to find '.fit' in %s\n",targetname);
        continue;
      }
      *charPtr4 = 0;
      strcat(targetname,".jpg");
    }
    if (JPGMode) {
      charPtr5 = strstr(rootname,".JPG");
      if (charPtr5 == NULL) {
        printf("ERROR: failed to fine 'JPG' in %s\n",targetname);
        continue;
      }
      *charPtr5 = 0;
      charPtr4 = strstr(targetname,".JPG");
      if (charPtr4 == NULL) {
        printf("ERROR: failed to find '.JPG' in %s\n",targetname);
        continue;
      }
      *charPtr4 = 0;
      strcat(targetname,".jpg");
    }



    printf("line %d source name %s targetname %s\n",fileIndex,filebasename,targetname);
    result = stat(targetname,&filestats);
    if (result == 0) {
      if (diagnosticMode == 0) {
        printf("ERROR: target file already exists %s \n",targetname);
        continue;
      }
    }
    if (backupMode != 0) {
      if (tarballString[targetIndex] == NULL) {
        tarballString[targetIndex] = (char*) calloc(TARGET_ALLOC_INCREMENT,sizeof(char));
        if (tarballString == NULL) {
          printf("ERROR: line %d allocation failure for tarballString\n",__LINE__);
          exit(-1);
        }
        tarballAlloc[targetIndex] = TARGET_ALLOC_INCREMENT;
      }
      if ((strlen(tarballString[targetIndex])+MAX_BUFFER+5) >= tarballAlloc[targetIndex]) {
        tarballAlloc[targetIndex] += TARGET_ALLOC_INCREMENT;
        tempString = (char *)calloc(tarballAlloc[targetIndex],sizeof(char));
        if (tempString == NULL) {
          printf("ERROR: line %d allocation failure for tempString\n",__LINE__);
          exit(-1);
        }
        strcpy(tempString,tarballString[targetIndex]);
        tarballString[targetIndex] = tempString;
      }
      if (targetCount[targetIndex] <= 0) {
        sprintf(tarballString[targetIndex],"tar -cvf %s_%03d.tar ",targetstring[targetIndex],numTarball);
        numTarball++;
      }
      charPtr8 = targetname+strlen(targetbasedirectory);
      if (charPtr8 == NULL) {
        printf("ERROR: line %d null charPtr8\n",__LINE__);
        exit(-1);
      }
      strcat(tarballString[targetIndex],"./");
      strcat(tarballString[targetIndex],charPtr8);
      strcat(tarballString[targetIndex]," ");
      targetCount[targetIndex]++;
      if (targetCount[targetIndex] >= MAX_TARGET_COUNT) {
        printf("%s\n",tarballString[targetIndex]);
        tarballString[targetIndex][0] = 0;
        targetCount[targetIndex] = 0;
      }
      continue;
    }

      

  


    /* Now get the date and time */
    strcpy(datename,targetname);
    strcat(datename,".date.txt");
    result = stat(datename,&filestats);
    if (result == 0) {
      if (diagnosticMode == 0) {
        printf("ERROR: date file %s already exists\n",datename);
        continue;
      }
    }
    if (fitsFlag == 0) {
      /* Get the date and time */
      sprintf(cmdStr,"exiftool -s3 -DateTimeOriginal %s > %s\n",filebasename,datename);
    } else { 
      /* get the plate center */
      charPtr9 = strrchr(origrootname,'_');
      if (charPtr9 == NULL) {
        printf("No full mosaic name for %s\n",origrootname);
        continue;
      }
      /* Convert from 16 bin to unbinned */
      charPtr9[1] = '0';
      charPtr9[2] = '1';
      sprintf(cmdStr,"getlocation %s > %s\n",origrootname,datename);
    }
    result = invokeCmd(cmdStr,&timeSpent);
    if (result != 0) {
      printf("ERROR: result %d for %s\n",result,cmdStr);
      strcpy(datestring,"NO LOCATION");
    } else {
      /* Now read the date from the file */
      dateHandle = fopen(datename,"rt");
      if (dateHandle == NULL) {
        printf("ERROR: Could not open file %s\n",datename);
        continue;
      }
      inBuffer = fgets(inLine,MAX_BUFFER,dateHandle);
      if (inBuffer == NULL) {
        printf("ERROR: failed to read from %s\n",datename);
        fclose(dateHandle);
        continue;
      }
      if (fitsFlag) {
        if (strstr(inBuffer,"WARNING") != NULL) {
          inBuffer = fgets(inLine,MAX_BUFFER,dateHandle);
          if (inBuffer == NULL) {
            printf("ERROR: failed to read from %s\n",datename);
            fclose(dateHandle);
            continue;
          }
        }
      }
      fclose(dateHandle);

      linelen = strlen(inBuffer);
      /* Trim off the carriage return */
      if (inBuffer[linelen-1] == 10) {
        inBuffer[linelen-1] = 0;
        linelen--;
      }
      /* Trim off the line feed */
      if (inBuffer[linelen-1] == 13) {
        inBuffer[linelen-1] = 0;
        linelen--;
      }
      if (fitsFlag) {
        /* look for the second space */
        charPtr6 = strchr(inBuffer,' ');
        if (charPtr6 == NULL) {
          printf("ERROR: failed to find first space in %s for %s\n",inBuffer,rootname);
          continue;
        }
        charPtr6++;
        charPtr6 = strchr(charPtr6,' ');
        if (charPtr6 == NULL) {
          printf("ERROR: failed to find second space in %s for %s\n",inBuffer,rootname);
          continue;
        }
        *charPtr6 = 0;
      }

      if (strlen(inBuffer) >= MAX_DATE_STR) {
        printf("ERROR: date string in %s %s is too long\n",datename,inBuffer);
        continue;
      }
      strcpy(datestring,inBuffer);
    }
    unlink(datename);
    if (fitsFlag == 0) {
      /* Now create the file */
      /* Originally 2456 */

      sprintf(cmdStr,"jpegtopnm %s | ppmtopgm  | pamscale -quiet -ysize=1228 -nomix | cjpeg -quality 100   > %s",filebasename,targetname);
    } else {
      /* fitsFlag = 1; Create the thumbnail */
      /* Originally ysize = 612 on the DASCH website */
      sprintf(cmdStr,"fitstopnm -quiet %s | pamscale -quiet -ysize=1224 -nomix | pnmnorm -quiet | pnmdepth -quiet 255 | pnminvert | pnmflip -tb  | pnmgamma 2.2  | cjpeg > %s ",filebasename,targetname);
    }
    if ((strlen(cmdStr) > MAX_CMD_STR)) {
      printf("ERROR: cmdStr is too long size %d for %s\n",strlen(cmdStr),cmdStr);
      exit(-1);
    }
    result = invokeCmd(cmdStr,&timeSpent);
    if (result != 0) {
      printf("ERROR: result %d for %s\n",result,cmdStr);
      continue;
    }
    result = stat(targetname,&filestats);
    if (result != 0) {
      printf("ERROR: failed to create target file %s\n",targetname);
      continue;
    }
    outputBytes += filestats.st_size;
    if (fitsFlag == 0) {
      sprintf(cmdStr,"exiftool -overwrite_original -DateTimeOriginal='%s' -Description='%s' %s",datestring,locationstring,targetname);
    } else {
      if (commentPtr == NULL) {
        commentPtr = blankString;
      }
      sprintf(cmdStr,"exiftool -overwrite_original  -Description='%s' -Location='%s' %s",datestring,commentPtr,targetname);
    }
    result = invokeCmd(cmdStr,&timeSpent);
    if (result != 0) {
      printf("ERROR: result %d for %s\n",result,cmdStr);
      unlink(targetname);
      continue;
    }
   


    if (diagnosticMode == 0) {
      /*  CopyFile(sourcename,targetname,forwardMode,&copycount,&deletecount,&workDone,verbose,&totalBytes); */
    } else {
      /* printf("DIAG: target %s\n",targetname); */
    }
  }
  if (backupMode != 0) {
    for (targetIndex = 0; targetIndex < numTarget; targetIndex++) {
      if (targetCount[targetIndex] > 0) {
        printf("%s\n",tarballString[targetIndex]);
        tarballString[targetIndex][0] = 0;
        targetCount[targetIndex] = 0;
      }
    }
  }



  time(&curTime);
  curTime -= startTime;
  /* printf("convertjpeg completed in %d seconds %d checks %d copies %d deleted for %s %s  %lld bytes copied\n",curTime,testcount,copycount,deletecount,lastfilename,source,totalBytes); */
  printf("convertjpeg completed in %lld seconds files %d maxString %d\n",curTime,nlines,maxString);
  if (outputBytes > 0) {
    printf("inputBytes %lld outputBytes %lld ratio %f\n",inputBytes,outputBytes,(1.0*inputBytes)/(1.0*outputBytes));
  }
  if (listStrings != NULL) {
    free(listStrings);
  }

  return(0);
}
