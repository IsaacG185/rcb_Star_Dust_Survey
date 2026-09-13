// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* filetest.c
 *  This program investigates the possibility of storing the magnitudes MySQL table as a collection of files.
 * 
 *  gcc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql  -I /dasch/install/include  -L /dasch/install/lib -lm -L/usr/lib64/mysql  -lmysqlclient filetest.c pipelineutils.a -ltable -lutil -lwcs -o filetest 
 *
 *  ./filetest
 *
 *         od -A d -t x4  mag000000000.dat
 *                                blocks             used         available
 *  Starting with a new ext3 file system: 984,296,796,160     209,477,632  934,087,901,184   1% /dasch/filetest
 *  2189 seconds to create 675878 dirs:   984,296,796,160   2,991,730,688  931,305,648,128   1% /dasch/filetest
 *   571 seconds to create 168967 dirs:   984,296,796,160     905,023,488  933,392,355,328   1% /dasch/filetest
 *
 *  filetest creating files               984,296,796,160   2,901,553,152  931,395.825.664   1% /dasch/filetest
 *  Exiting with result -1 errno 28 No space left on device 60,869,446 files  start 11:00 to 12:57
 *                                                  need   168,966,366 files can store 10 bins/file
 *                                        9822 sec for 60211312 removeFile
 *                                         499 sec for   168650 removedirectory
 *                                           4 sec for      169 removedirectory (top level)
 *
 *  Modulus 1024 for directories, 16 for files:
 *  2731 seconds with 165,007 directories and 10,560,400 file
 *                                        984,296,796.160     889,962,496  933,407,416,320
 *  Return to modulus 10 for directories.
 *                                       12095 sec for  168967 directories 21120800 files (408 bytes/file)
 *                                        984,296,796,160     44,161,802,240 890135576576   5% /dasch/filetest
 * 
 * Aug  28, 2009 Edward J. Los - Initial version
 */


#include <math.h>
#include <unistd.h>
#include "table.h"
#include "mysql.h"
#include "pipelineutils.h"
#include "photometryutils.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h> /* needed for exit */
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#define MAG_SUBDIR_MODULUS 1000
#define MAG_FILE_MODULUS 1024
#define OLD_FILE_MODULUS 16
/* #define ORIGINAL_CODE 1 */
#define MODULUS_REPAIR_CODE 1
#define OLD_MODULUS_RATIO (MAG_FILE_MODULUS/OLD_FILE_MODULUS)

#define MAX_BUFFER 512
#define MAG_FORMAT_DATA   0x64617461 /* 'data' */
#if 1
#define START_BIN 0
#define END_BIN TOTAL_GSC_BINS_64
#endif
#if 0
#define START_BIN 107778048
#define END_BIN   107779072
#endif
#if 0
#define START_BIN 108004352
#define END_BIN  (108004352+1024)
#endif
#if 0
#define START_BIN 57032704
#define END_BIN  (57032704+1024)
#endif


/* Disk resident file structure.  NOTE!!! observe the decreasing size of field so that this structure
 * is packed the same on 32 bit and 64 bit systems 
 * The normalized file size is 272 bytes, but for performance, fields designated with a "+"
 * are added for faster extraction of lightcurve data to bring the total to 308 bytes, a 13% increase
 */



int processingStep = 0;

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
    printf("Terminated with signal %d at %d seconds step %d\n",WTERMSIG(result),(int)curTime,processingStep);
    exit(-1);
  }
#if 0
  if (WIFEXITED(result)) {
    printf("Exited with status %d at %d sec delta %d sec step %d\n",WEXITSTATUS(result),(int)curTime,(int)beginTime,processingStep);
  } else {
    printf("Unknown result %d at %d sec delta %d step %d\n",result,(int)curTime,(int)beginTime,processingStep);
  }
  printf("result: %d\n",result);
#endif

  return(result);
}

#ifdef ORIGINAL_CODE
int main(int argc,char *argv[])
{
  time_t startTime;
  time_t curTime;
  int gsc_bin_number;
  int directoryCounter = 0;
  int fileCounter = 0;
  int removeCounter = 0;
  int dir1;
  int dir2;
  char cmdStr[MAX_BUFFER];
  char filename[MAX_BUFFER];
  int timeSpent;
  int result;
  int fileNumber;
  ssize_t writeBytes;
  FILESTARIMAGE starImage;
  PFILESTARIMAGE pStarImage = &starImage;
  struct stat filestats;
  memset(pStarImage,0,sizeof(FILESTARIMAGE));
  long long versionTag;                            /* Version tag */
  char *photfilebase;

  versionTag = MAG_FORMAT_VERSION;
  versionTag = versionTag << 32;
  versionTag = versionTag | MAG_FORMAT_DATA;

  photfilebase = getenv("DASCH_PHOT_MAGNITUDES");
  if (photfilebase == NULL) {
    fprintf(stderr,"DASCH_PHOT_MAGNITUDES is not defined\n");
    return(-1);
  }

  time(&startTime);

  
#if 0
  printf("filetest creating directories\n");
  for (gsc_bin_number = 0; gsc_bin_number < TOTAL_GSC_BINS_64; gsc_bin_number+= MAG_SUBDIR_MODULUS) {
    dir1 = gsc_bin_number / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
    dir2 = gsc_bin_number / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
    sprintf(cmdStr,"mkdir -p %s/mag%03d/mag%03d",photfilebase,dir1,dir2);
    result = invokeCmd(cmdStr,&timeSpent);
    if (result != 0) {
      printf("Exiting with result %d\n");
    }
    directoryCounter++;
#if 0
    if (directoryCounter > 20) {
      break;
    }
#endif
  }

#endif


#if 0
  printf("filetest creating files\n");
  for (gsc_bin_number = 0; gsc_bin_number < TOTAL_GSC_BINS_64; gsc_bin_number+=MAG_FILE_MODULUS) {
    dir1 = gsc_bin_number / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
    dir2 = gsc_bin_number / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
    sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,gsc_bin_number);
    fileNumber = open(filename,O_CREAT|O_RDWR|O_APPEND,S_IRWXU|S_IRGRP);
    if (fileNumber < 0) {
      printf("Exiting with result %d errno %d %s\n",fileNumber,errno,strerror(errno));
      exit(-1);
    }
    fileCounter++;
    close(fileNumber);
#if 0
    if (fileCounter > 10) {
      break;
    }
#endif
  }

#endif


#if 0
  printf("filetest writing blocks of size %d\n",sizeof(FILESTARIMAGE));
  for (gsc_bin_number = 0; gsc_bin_number < TOTAL_GSC_BINS_64; gsc_bin_number+=MAG_FILE_MODULUS) {
    dir1 = gsc_bin_number / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
    dir2 = gsc_bin_number / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
    sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,gsc_bin_number);
    fileNumber = open(filename,O_CREAT|O_RDWR|O_APPEND,S_IRWXU|S_IRGRP);
    if (fileNumber < 0) {
      printf("Exiting with result %d errno %d %s for %s\n",fileNumber,errno,strerror(errno),filename);
      exit(-1);
    }
    fileCounter++;
    writeBytes = write(fileNumber,pStarImage,sizeof(FILESTARIMAGE));
    if (writeBytes != sizeof(FILESTARIMAGE)) {
      printf("Wrote only %d bytes to %s\n",writeBytes,filename);
    }



    close(fileNumber);
#if 0
    if (fileCounter > 10) {
      break;
    }
#endif
  }

#endif

#if 0
  printf("filetest write format test with size %d\n",sizeof(FILESTARIMAGE));
  gsc_bin_number = 0;
  dir1 = gsc_bin_number / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
  dir2 = gsc_bin_number / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
  sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,gsc_bin_number);
  fileNumber = open(filename,O_CREAT|O_RDWR|O_APPEND,S_IRWXU|S_IRGRP);
  if (fileNumber < 0) {
    printf("Exiting with result %d errno %d %s for %s\n",fileNumber,errno,strerror(errno),filename);
    exit(-1);
  }
  fileCounter++;
  pStarImage->versionTag = versionTag;
  pStarImage->seriesId = 1;
  pStarImage->plateNumber = 2;
  pStarImage->exposureNumber = 3;
  pStarImage->local_bin_index = 4;
  pStarImage->REFNumber = 6;        /* Translated reference number */
  pStarImage->NUMBER = 7;                 /* Sextractor reference number */
  pStarImage->passBits = 8;               /* Pass identifier */
  pStarImage->X_IMAGE = 9;
  pStarImage->Y_IMAGE = 10;
  pStarImage->AFLAGS = 11;                 /* Flags keyword. See header of pipelineutils.h */
  pStarImage->BFLAGS = 12;                 /* Flags keyword. See header of pipelineutils.h */
  pStarImage->MAG_ISO = 14;
  pStarImage->ra = 15;                  /* Right Ascension in degrees */
  pStarImage->dec = 16;                 /* Declination in degrees */
  pStarImage->magcal_iso = 17;          /* lowess magnitude estimate */
  pStarImage->magcal_iso_rms = 18;      /* local  error */
  pStarImage->magcal_local = 19;        /* Local magnitude calibration */
  pStarImage->magcal_local_rms = 21;    /* Overall error */
  pStarImage->spatial_bin = 22;            /* Lowess magnitude spatial bin  number */
  pStarImage->Date = 23;                /* Heliocentric Julian Date  */      
  pStarImage->limiting_mag_local = 24;  /* Limiting magnitude */
  pStarImage->FLUX_ISO = 30;
  pStarImage->MAG_APER = 31;
  pStarImage->MAG_AUTO = 32;
  pStarImage->KRON_RADIUS = 33;
  pStarImage->BACKGROUND = 34;
  pStarImage->FLUX_MAX = 36;
  pStarImage->THETA_J2000 = 37;
  pStarImage->ELLIPTICITY = 38;
  pStarImage->ISOAREA_WORLD = 39;
  pStarImage->FWHM_IMAGE = 40;
  pStarImage->FWHM_WORLD = 41;
  pStarImage->ISO0 = 42;
  pStarImage->ISO1 = 43;
  pStarImage->ISO2 = 44;
  pStarImage->ISO3 = 45;
  pStarImage->ISO4 = 46;
  pStarImage->ISO5 = 47;
  pStarImage->ISO6 = 48;
  pStarImage->ISO7 = 49;
  pStarImage->plate_dist = 50;
  pStarImage->Blendedmag = 51;
  pStarImage->gsc_bin_index = 52;
  pStarImage->versionId = 53;
  pStarImage->dra = 54;
  pStarImage->ddec = 55;
 
  writeBytes = write(fileNumber,pStarImage,sizeof(FILESTARIMAGE));
  if (writeBytes != sizeof(FILESTARIMAGE)) {
    printf("Wrote only %d bytes to %s\n",writeBytes,filename);
  }
  result = stat(filename,&filestats);
  
  printf("File %s size is %d bytes status %d\n",filename,filestats.st_size,result);


  close(fileNumber);


#endif





#if 0
  printf("filetest removing files\n");
  for (gsc_bin_number = 0; gsc_bin_number < TOTAL_GSC_BINS_64; gsc_bin_number+=MAG_FILE_MODULUS) {
    dir1 = gsc_bin_number / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
    dir2 = gsc_bin_number / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
    sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,gsc_bin_number);
    fileNumber = remove(filename);
#if 0
    if (fileNumber < 0) {
      printf("Exiting with result %d errno %d %s\n",fileNumber,errno,strerror(errno));
      exit(-1);
    }
#endif
    if (fileNumber >= 0) {
      removeCounter++;
    }
#if 0
    if (removeCounter > 10) {
      break;
    }
#endif
  }

#endif


#if 1
  printf("filetest removing directories\n");
  for (gsc_bin_number = 0; gsc_bin_number < TOTAL_GSC_BINS_64; gsc_bin_number+= MAG_SUBDIR_MODULUS) {
    dir1 = gsc_bin_number / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
    dir2 = gsc_bin_number / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
    sprintf(filename,"%s/mag%03d/mag%03d",photfilebase,dir1,dir2);
    fileNumber = remove(filename);
#if 0
    if (fileNumber < 0) {
      printf("Exiting with result %d errno %d %s\n",fileNumber,errno,strerror(errno));
      exit(-1);
    }
#endif
    if (fileNumber >= 0) {
      removeCounter++;
    }


  }

#endif
#if 0
  printf("filetest removing root directories\n");
  for (gsc_bin_number = 0; gsc_bin_number < TOTAL_GSC_BINS_64; gsc_bin_number+= (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS)) {
    dir1 = gsc_bin_number / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
    dir2 = gsc_bin_number / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
    sprintf(filename,"%s/mag%03d",photfilebase,dir1);
    fileNumber = remove(filename);
#if 0
    if (fileNumber < 0) {
      printf("Exiting with result %d errno %d %s\n",fileNumber,errno,strerror(errno));
      exit(-1);
    }
#endif
    if (fileNumber >= 0) {
      removeCounter++;
    }


  }

#endif


  time(&curTime);
  curTime -= startTime;
  printf("filetest completed in %d seconds with %d directories and %d files and %d removeFile\n",
         curTime,directoryCounter,fileCounter,removeCounter);
}
#endif /* ORIGINAL_CODE */
#ifdef MODULUS_REPAIR_CODE
int main(int argc,char *argv[])
{
  time_t startTime;
  time_t curTime;

  int dir1;
  int dir2;
  int base_gsc_bin;
  int gsc_bin_index;
  int fileindex;
  off_t filesize[OLD_MODULUS_RATIO];
  char filename[MAX_BUFFER];
  char newfile[MAX_BUFFER];
  int statResult;
  struct stat statbuf;
  int subfilecount;
  char *mysqlphothost;
  char *photusername;
  char *photpassword;
  char *photfilebase;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  char queryString[MAX_QUERY_STRING];
  int res;
  int gotAnswer;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  int statCounter = 0;
  int fileCounter = 0;
  int cleanBinCount = 0;
  int binCounter = 0;
  char *filebuffer = NULL;
  int bufferAlloc = 0;
  int magnitudeCount = 0;
  off_t maxsize;
  int newFileNumber;
  int oldFileNumber;
  ssize_t writeBytes;
  ssize_t readBytes;
  int result;
  int errorFlag;
  char cmdStr[MAX_BUFFER];
  int timeSpent;
  int directoryCounter = 0;

  time(&startTime);


  printf("filetest modulus repair code with size %d\n",sizeof(FILESTARIMAGE));
  photusername = getenv("DASCH_PHOT_USERNAME");
  if (photusername == NULL) {
    fprintf(stderr,"DASCH_PHOT_USERNAME is not defined\n");
    return(-1);
  }
  photpassword = getenv("DASCH_PHOT_PASSWORD");
  if (photpassword == NULL) {
    fprintf(stderr,"DASCH_PHOT_PASSWORD is not defined\n");
    return(-1);
  }
  mysqlphothost = getenv("DASCH_PHOT_MYSQLHOST");
  if (mysqlphothost == NULL) {
    fprintf(stderr,"DASCH_PHOT_MYSQLHOST is not defined\n");
    return(-1);
  }
  photfilebase = getenv("DASCH_PHOT_MAGNITUDES");
  if (photfilebase == NULL) {
    fprintf(stderr,"DASCH_PHOT_MAGNITUDES is not defined\n");
    return(-1);
  }

  mysql_init(pPhotConnection);
  
  if (!mysql_real_connect(pPhotConnection,mysqlphothost,photusername,photpassword,"photometry",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pPhotConnection)) {
      fprintf(stderr,"ERROR: MySQL error %d: %s\n",mysql_errno(pPhotConnection),mysql_error(pPhotConnection));
    }
    return(-1);
  }


  sprintf(queryString,"update photglobal set keepalive = 'yes';");
  res = ExecuteQuery(pPhotConnection,queryString);

  for (base_gsc_bin = START_BIN; base_gsc_bin < END_BIN; base_gsc_bin+= MAG_FILE_MODULUS) {
    gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
    if ((gotAnswer != 1)  || (pPhotGlobal->keepalive != PHOT_KEEPALIVE_YES)) {
      printf("Terminating because of keepalive failure\n");
      break;
    }
    errorFlag = 0;
    subfilecount = 0;
    for (fileindex = 0; fileindex < OLD_MODULUS_RATIO; fileindex++) {
      filesize[fileindex] = 0;
      gsc_bin_index = base_gsc_bin + (OLD_FILE_MODULUS * fileindex);
      dir1 = gsc_bin_index / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
      dir2 = gsc_bin_index / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
      sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,gsc_bin_index);
#if 0
      printf("Checking %s\n",filename);
#endif
      statResult = stat(filename,&statbuf);
      statCounter++;
      if (statResult == 0) {
#if 0
        printf("Found %s\n",filename);
#endif
        fileCounter++;
        filesize[fileindex] = statbuf.st_size;
        magnitudeCount += statbuf.st_size/sizeof(FILESTARIMAGE);
        if ((filesize[fileindex] % sizeof(FILESTARIMAGE)) != 0) {
          printf("ERROR: File %s had a bad size %d\n",filename,filesize[fileindex]);
          errorFlag = 1;
          filesize[fileindex] = 0;
        } else {
          subfilecount++;
        }

      }
    }
#if 0
    if ((subfilecount > 1) && (filesize[0] > 0)) {
      printf("Base index %9d has %2d files\n",base_gsc_bin,subfilecount);
    }
#endif
    if (subfilecount == 0) {
      continue;
    }
    if ((subfilecount == 1) && (filesize[0] > 0)) {
      cleanBinCount++;
      continue;
    }

    binCounter++;
    maxsize = 0;
    for (fileindex = 0; fileindex < OLD_MODULUS_RATIO; fileindex++) {
      if (filesize[fileindex] > maxsize) {
        maxsize = filesize[fileindex];
      }
    }
    if (maxsize > bufferAlloc) {
      if (filebuffer != NULL) {
        free(filebuffer);
      }
      filebuffer = (char *)malloc(maxsize);
      if (filebuffer == NULL) {
        printf("ERROR: failed to allocate a buffer of size %d\n",maxsize);
        exit(-1);
      }
      bufferAlloc = maxsize;
    }
      
    gsc_bin_index = base_gsc_bin;
    dir1 = gsc_bin_index / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
    dir2 = gsc_bin_index / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
    sprintf(newfile,"%s/mag%03d/mag%03d",photfilebase,dir1,dir2);
    statResult = stat(newfile,&statbuf);
    if (statResult != 0) {
      sprintf(cmdStr,"mkdir -p %s/mag%03d/mag%03d",photfilebase,dir1,dir2);
      result = invokeCmd(cmdStr,&timeSpent);
      if (result != 0) {
        printf("ERROR: Exiting with result %d for %s\n",result,cmdStr);
        exit(-1);
      }
      directoryCounter++;
    

    }
    sprintf(newfile,"%s/mag%03d/mag%03d/new%09d.dat",photfilebase,dir1,dir2,gsc_bin_index);
    newFileNumber = open(newfile,O_CREAT|O_WRONLY|O_TRUNC,S_IRWXU|S_IRGRP);
    if (newFileNumber < 0) {
      printf("ERROR creating file %s result %d errno %d %s\n",newfile,newFileNumber,errno,strerror(errno));
      exit(-1);
    }
    for (fileindex = 0; fileindex < OLD_MODULUS_RATIO; fileindex++) {
      if (filesize[fileindex] > 0) {
        gsc_bin_index = base_gsc_bin + (OLD_FILE_MODULUS * fileindex);
        dir1 = gsc_bin_index / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
        dir2 = gsc_bin_index / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
        sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,gsc_bin_index);
        oldFileNumber = open(filename,O_RDONLY,S_IRWXU|S_IRGRP);
        if (oldFileNumber < 0) {
          printf("ERROR creating file %s result %d errno %d %s\n",filename,oldFileNumber,errno,strerror(errno));
          exit(-1);
        }
        readBytes = read(oldFileNumber,filebuffer,filesize[fileindex]);
        if (readBytes != filesize[fileindex]) {
          printf("ERROR: Read only %d bytes from %s\n",readBytes,filename);
          exit(-1);
        }
        writeBytes = write(newFileNumber,filebuffer,filesize[fileindex]);
        if (writeBytes != filesize[fileindex]) {
          printf("ERROR: Wrote only %d bytes to %s\n",writeBytes,newfile);
          exit(-1);
        }
        close(oldFileNumber);
      }
    }

    close(newFileNumber);
    if (errorFlag == 0) {
      gsc_bin_index = base_gsc_bin;
      dir1 = gsc_bin_index / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
      dir2 = gsc_bin_index / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
      sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,gsc_bin_index);
      result = rename(newfile,filename);
      if (result != 0) {
        printf("ERROR: Exiting with result %d rename %s %s\n",result,newfile,filename);
        exit(-1);
      }
      for (fileindex = 1; fileindex < OLD_MODULUS_RATIO; fileindex++) {
        if (filesize[fileindex] > 0) {
          gsc_bin_index = base_gsc_bin + (OLD_FILE_MODULUS * fileindex);
          dir1 = gsc_bin_index / (MAG_SUBDIR_MODULUS*MAG_SUBDIR_MODULUS);
          dir2 = gsc_bin_index / (MAG_SUBDIR_MODULUS) - (dir1 * MAG_SUBDIR_MODULUS);
          sprintf(filename,"%s/mag%03d/mag%03d/mag%09d.dat",photfilebase,dir1,dir2,gsc_bin_index);
          result = unlink(filename);
          if (result != 0) {
            printf("ERROR: Exiting with result %d for unlink %s\n",result,filename);
            exit(-1);
          }

        }

      }
    }

  }
  free(filebuffer);
  mysql_close(pPhotConnection);

  time(&curTime);
  curTime -= startTime;
  printf("filetest completed in %d seconds with %d stat checks and %d files and %d gsc bins %d clean bins %d max buffer %d magnitudes %d directories\n",
         curTime,statCounter,fileCounter,binCounter,cleanBinCount,bufferAlloc,magnitudeCount,directoryCounter);
  return(0);
}

#endif /* MODULUS_REPAIR_CODE */
