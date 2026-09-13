// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* filterdr1.c
 *
 *    Get release field star statistics.  Input is a starbase table with the headers  "ra" and "declination" extracted from an id table.
 *
 *  gcc -ggdb -O2  -I/usr/include/mysql -I /dasch/install/include  -L /dasch/install/lib filterdr1.c pipelineutils.a  -L/usr/lib64/mysql  -lmysqlclient -ldl -pthread  -lwcs -lm -o filterdr1
 * 
 * May 16, 2013 Edward J. Los - Original Version
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#define MAX_TAB 3
#define EXP_TAB 2
#define MAX_CATALOG_LINE 512
#include "pipelineutils.h"
extern char *releaseFieldText[RELEASE_FIELD_MAX+1];
int main(int argc,char *argv[]) {
  char *pTab[MAX_TAB];
	int tabCount;
  char *inBuffer;
  char inLine[MAX_CATALOG_LINE];
  char copyLine[MAX_CATALOG_LINE];
	int lineCounter = 0;
	int nvals;
	int errorFlag = 0;
  FILE * inHandle;
	int lineLen;
	int tabIndex;
	int releaseFieldPlates[RELEASE_FIELD_MAX+1];
	double ra;
	double declination;
	int releaseField;


	memset(releaseFieldPlates,0,sizeof(releaseFieldPlates));
	if (argc < 2) {
		printf("Usage filterdr1 infile\n");
		exit(-1);
	}
	
  inHandle = fopen(argv[1],"rt");
	if (inHandle == NULL) {
		printf("ERROR: failed to open %s\n",argv[1]);
		exit(-1);
	}

  while (1) {
    tabCount = 0;
    inBuffer = fgets(inLine,MAX_CATALOG_LINE,inHandle);
    if (inBuffer == NULL) {
      break;
    }
    lineLen = strlen(inBuffer);
    /* Trim off the carriage return */
    if (inBuffer[lineLen-1] == 10) {
      inBuffer[lineLen-1] = 0;
      lineLen--;
    }
    strcpy(copyLine,inBuffer);

    pTab[0] = inBuffer;
		lineCounter++;
		if ((lineCounter % 500000) == 0) {
			printf("At line %d\n",lineCounter);
		}
		
    for (tabIndex = 0; tabIndex < lineLen; tabIndex++) {
      if (inBuffer[tabIndex] == 9) {
        /* This is a tab */
        inBuffer[tabIndex] = 0;
        tabCount++;
        if (tabCount >= MAX_TAB) {
					printf("ERROR: too many tabs\n");
          errorFlag = 1;
        } else {
          pTab[tabCount] = &inBuffer[tabIndex+1];
        }
      }

    }
    if (tabCount < EXP_TAB) {
			printf("ERROR: not enough tabs, got %d\n",tabCount);
      errorFlag = 1;
    } 
		if (lineCounter == 1) {
			printf("tab 1 is %s\n",pTab[1]);
			printf("tab 2 is %s\n",pTab[2]);
		}
		if (lineCounter > 2) {
			nvals = sscanf(pTab[1],"%lf",&ra);
			if (nvals != 1) {
				continue;
			}
			nvals = sscanf(pTab[2],"%lf",&declination);
			if (nvals != 1) {
				continue;
			}
			CheckAuthorization(NULL,ra,declination,&releaseField);
			if ((releaseField >= 0) && (releaseField < RELEASE_FIELD_MAX)) {
				releaseFieldPlates[RELEASE_FIELD_MAX]++;
				releaseFieldPlates[releaseField]++;
				
			}
			
		}
	}
	printf("total lines %d\n",lineCounter++);
	for (releaseField = 0; releaseField <= RELEASE_FIELD_MAX; releaseField++) {
		printf("%6d lightcurves for field %s\n",releaseFieldPlates[releaseField],releaseFieldText[releaseField]);
	}

}
