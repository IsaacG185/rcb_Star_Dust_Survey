// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*  parse_coo.c
 *
 *    Parse coordinates
 *

gcc -fPIC -ggdb -c -O2   -I/n/sw/plplot-5.9.9/include/plplot -I/usr/include/mysql -I/dasch/install/include -I/usr/include/plplot parse_coo.c  -D_FILE_OFFSET_BITS=64 -o parse_coo.o

gcc -fPIC -ggdb -c -O2  -DPARSECOOSTUB -I/n/sw/plplot-5.9.9/include/plplot -I/usr/include/mysql -I/dasch/install/include -I/usr/include/plplot parse_coo.c  -D_FILE_OFFSET_BITS=64 -o parse_coostub.o


 *
 * May  5, 2010 Edward J. Los - Adapted from find_lightcurves.c and asas_cat_input provided by
 *                              Grzegorz Pojmanski <gp@astrouw.edu.pl>  ASAS source author on Apr 26, 2010
 * Ded 11, 2020 Edward J. Los - Moved to a separate module
 * 
 */ 
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#if 0
#include "mysql.h"
#endif
#include <sys/types.h>
#include <grp.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <assert.h>

#include <fcntl.h>
#include <limits.h>
#include <assert.h>
#if 0
#include "pipelineutils.h"
#include "photometryutils.h"
#endif
		/* Return RA in degrees from string */
double str2ra(const char* in);	/* Character string (hh:mm:ss.sss or dd.dddd) */
		/* Return Dec in degrees from string */
    double str2dec(const char* in);	/* Character string (dd:mm:ss.sss or dd.dddd) */
int isalpha(int c);

/* Note parse_coo is the same as the routine in web_query.c and was originally written by Grzegorz Pojmanski */
int parse_coo(coi,a,d)
char *coi;
double *a, *d;
{
#ifdef PARSECOOSTUB
  printf("ERROR: Illegal call to PlotTransientCandidates\n");
  exit(-1);
#else /* PARSECOOSTUB */
  char *cp, *cp0, *cp1, coo[256];
  cp0 = coi;
  cp1 = coo;
  while (*cp0 == ' ') {
    cp0++;
  }
  while (*cp0 != 0) {
    *cp1++ = *cp0++;
  } 
  *cp1++ = *cp0++;
  cp1 = coo + strlen(coo);
  while(*cp1 == ' '){
    *cp1 = '\0';
    if(cp1 == coo)break;
  }
  if(isalpha(coo[0]))return(-1);
  
  
  {
    char buf[256];
#if 0
    double ra,dec;
#endif
    int cnt;

    strcpy(buf,coo);
    /* get tokens */
    cnt = 0;
    cp0 = buf;
    for(;;){
      cp = strtok(cp0," ");
      if(cp == NULL) break;
      cnt++;
      cp0 = NULL;
    }
    if(cnt <=2 ){
      strcpy(buf,coo);
      cp0 = strtok(buf,"+-,; \t");
      if(cp0 == NULL) return(-1);
      cp1 = strtok(NULL,",; \t");
      if(cp1 == NULL) return(-1);
#if 0
      ra = str2ra(cp0);
      dec = str2dec(cp1);
      *a = (double)ra;
      *d = (double)dec;
#else
      *a=(double)str2ra(cp0);
      *d=(double)str2dec(cp1);
#endif
      return(0);
    }else{
      char cpp0[256], cpp1[256];
      strcpy(buf,coo);
      cp1 = strchr(buf,'+');
      if( cp1 == buf){
        cp1 = strchr(buf+1,'+');
      }else if(cp1 == NULL){
        cp1 = strchr(buf,'-');
      }
      if(cp1 != NULL){
        strcpy(cpp1, cp1);
        strcpy(cpp0, buf);
        cpp0[cp1-buf]='\0';
        *a=(double)str2ra(cpp0);
        *d=(double)str2dec(cpp1);
        return(0);
      }else{
        return(-1);
      }
    }
  }
#endif /* PARSECOOSTUB */
}
