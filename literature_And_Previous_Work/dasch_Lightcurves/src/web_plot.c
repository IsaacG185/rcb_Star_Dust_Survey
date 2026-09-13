// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include <table.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>
#include <libwcs/wcscat.h>

#include <gd.h>
#include <gdfontg.h>
#include <gdfonts.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include "photometryutils.h"
#include "searchgsc.h"
#include "daschunistd.h"

double fmax(double x, double y);
double fmin(double x, double y);

#define MAX_BUFFER 500
#define MAX_FILENAME 512
#define LIST_STRING_ALLOC 2000
#define LIST_SUBSTR_SIZE 28
#define MAX_MOSAICID_STRING 50
#define MAX_PERIOD_STRING 32 /* VSX has a period field with a maximum of 16 characters */
#define ERROR_BAR_CLIP 0.4

/* #define ACCEPT_LOW_ALTITUDE */

#define PLOT_AVE	0x01
#define PLOT_FLIP	0x02
#define PLOT_LINE	0x04
#define PLOT_EXACT	0x08
#define PLOT_UPPER_L	0x10
#define MAX_LABEL 128
#define UPPER_LIMIT 19.0

#define COLOR_BLACK  0    /* AFLAGS = 0 single-catalog display */
                          /* AFLAGS = 0 gsc2.3.2, combined display */
#define COLOR_BLUE   1    /* too bright  and limiting mangitude*/
#define COLOR_RED    2    /* AFLAGS != 0, single-catalog display*/
#define COLOR_GREEN  3    /* AFLAGS != 0 for apass, combined display */
#define COLOR_GRAY   5    /* undected */
#define COLOR_YELLOW 6    /* AFLAGS != 0 experimental, combined display */
#define COLOR_PURPLE 7    /* AFLAGS = 0 experimental, combined display */
#define COLOR_CYAN   8    /* AFLAGS = 0 apass combined display */
#define COLOR_MAX    9


extern FLAGSENTRY AflagsTable[];
extern int AflagsTableSize;
extern QUALITYBIT webQualityMasks[];
extern int webQualityTableSize;
extern char *catalogText[MAX_CATALOG_NUMBER];
extern GSCBIN gscBin64;


typedef struct sortTable {
  double julianDate; /* This item must be first for sorting */
  char listString[LIST_SUBSTR_SIZE+1]; /* Pointer to the list string in tclass */
} SORTTABLE,*PSORTTABLE;


static PGSCBIN pGscBin = &gscBin64;
static int ndec = 1;		/* Number of decimal places in non-angles */
static SERIESLIST seriesList[MAX_SERIES+1];
static int seriesCount = 0;
static int allSeriesMask = 0;


/* This routine was authored by  Grzegorz Pojmanski */
void dsortindx(int n,
         double *arrin,
         int *indx)
{
  int l,ir,i,j,indxt;
  double q;

  for ( i=0; i<n; i++ ) indx[i]=i;
  if ( n == 1 ) return;
  l=n/2;
  ir=n-1;
  while ( 10 ) {
    if ( l > 0 ) {
      l--;
      indxt=indx[l];
      q=arrin[indxt];
    } else {
      indxt=indx[ir];
      q=arrin[indxt];
      indx[ir]=indx[0];
      ir--;
      if ( ir == 0 ) {
        indx[0]=indxt;
        return;
      }
    }
    i=l;
    j=l+l+1;
    while ( j <= ir ) {
      if ( j < ir )
        if ( arrin[indx[j]] < arrin[indx[j+1]] ) j++;
      if ( q < arrin[indx[j]] ) {
        indx[i]=indx[j];
        i=j;
        j+=j+1;
      } else {
        j=ir+1;
      }
    }
    indx[i]=indxt;
  }
}

/* This routine was authored by  Grzegorz Pojmanski, but modified for the special needs of the DASCH project */

int PlotCurve(gdImagePtr im_out,
              int ix,
              int iy,
              int sx,
              int sy,
              int n,
              double *x,
              double *y,
              double *s,
              double *err,
              double *limiting,
              int *col,
              int *crossindex,
              PFILESTARIMAGEEXT db_table,
              double yymin,
              double yymax,
              int col1,
              int col2,
              int *colors,
              char *title,
              char *xlab,
              char *ylab,
              char *lab1,
              char *lab2,
              char *lab3,
              int style,
              char *store,
              double upper_limit,
              double xxmin,
              double xxmax,
              int grey,
              int verbose,
              int AFLAGSMASK,
              int qualitymask,
              int seriesmask,
              double reqstartdate,
              double reqenddate,
              double reqdimmag,
              double reqbrightmag,
              int reqXDown,
              int reqYDown,
              int reqXUp,
              int reqYUp,
              int plotLimiting,
              double averageRa,
              double averageDec,
              int foldflag,
              int enableRematch,
              double foldperiod,
              char *foldperiodString,
              double foldcenter)
{
  int i, x0, y0, x1, y1, xa, ya, xb, yb, len, yy;
  double xmax = -1.e32, xmin = 1.e32, ymax = -1.e32, ymin = 1.e32, dx, dy;
  double v, fx, fx0, fy, fy0, ym = 0, DX, DY, X0, Y0;
  int nt=20,NDX,NDY;
  char buf[32];
  int kkk;
  int colorCount[COLOR_MAX];

  memset (colorCount,0,sizeof(colorCount));

  x0=40;
  x1=sx-10;
  if (foldflag || enableRematch) {
    y0=45;
  } else {
    y0=35;
  }
  y1=sy-15;

  if(n > 0){
    int *indx, l;
    double ymed, sigma;
    indx = (int *)malloc(n*sizeof(int));
    dsortindx(n,x,indx);
    xmin = x[indx[0]];
    xmax = x[indx[n-1]];

    if ((xxmin == 0 && xxmax == 0) || (xxmin == -1.e32 && xxmax == 1.e32)) {
    } else {
      if (xxmin != -1.e32)
        xmin = xxmin;
      if (xxmax !=  1.e32)
        xmax = xxmax;
    }

    dsortindx(n,y,indx);
    l=n;
    if(style & PLOT_UPPER_L){
      for(l=n-1;l>=0;l--){
        if (y[indx[l]] < 29.99)break;
      }
    } else {
      l = n-1;
    }
    if(l < 0){
      ymin = upper_limit-0.2;
      ymax = upper_limit+0.2;
    }else
      if (style & PLOT_EXACT || l<20 ){
        for (kkk = 0; kkk < n; kkk++) {
          if (y[indx[kkk]] != 0) {
            break;
          }
        }
        if (kkk < n) {
          ymin = y[indx[kkk]];
          ymax = y[indx[l]];
          ym = y[indx[(int)(0.5*l)]];
        }
        ymin = reqbrightmag;
        ymax = reqdimmag;
        if ( l != n-1 ) upper_limit = ymax+0.1;

      }else{
        ymed = y[indx[(int)(0.5*l)]];
        sigma = (y[indx[(int)(0.83*l)]]-y[indx[(int)(0.17*l)]])/2.;
        if(verbose)printf("ymed, sigma %f %f\n",ymed,sigma);

        if(sigma == 0){
          ymin = y[indx[0]] - 0.2;
          ymax = y[indx[l]];
          if ( ymax < upper_limit ) ymax = upper_limit;
          upper_limit = ymax+0.1;
        }else{
          ymin = y[indx[(int)(0.00*l)]] - sigma;
          ymax = y[indx[(int)(0.97*l)]] + sigma;
          if(ymax < y[indx[l-1]] ) ymax = y[indx[l-1]] + 0.2;
          upper_limit = ymax+0.1;
          if(verbose)printf("0, 100 %f %f\n",y[indx[0]],y[indx[l-1]]);

        }
      }

    if(verbose)printf("ymin,ymax,upper_limit %f %f %f\n",ymin,ymax,upper_limit);
    free(indx);
  }else{
    if(verbose)printf("8\n");

    xmin=PIPELINE_MIN_DATE;
    xmax=PIPELINE_MAX_DATE;
    ymin=0; ymax=1.;
    ym=(ymax+ymin)/2.;
  }

  if(verbose)printf("b\n");
  dx=(xmax-xmin)*1.1;
  dy=(ymax-ymin)*1.1;
  if(dx == 0. ) dx=1.;
  if(dy < 0.16) {
    dy = 0.2;
  }else if (dy < .4){
    dy = 0.5;
  }else if (dy < 0.8){
    dy = 1.;
  }else if (dy < 2.){
    dy = 2.;
  }

  fx = (x1-x0)/dx;
  fx0 = (xmax+xmin)/2.-dx/2.;
  fy = (y1-y0)/dy;
  fy0 = (ymax+ymin)/2.-dy/2.;
  if((int)strlen(store) > 0){
    FILE *fp;
    if((fp=fopen(store,"w")) != NULL){
      fprintf(fp,"%d %d %d %lf %lf\n",sx,x0,x1,fx,fx0);
      fprintf(fp,"%d %d %d %lf %lf\n",sy,y0,y1,fy,fy0);
      fprintf(fp,"%lf %d %s %lf %d\n",upper_limit,foldflag,foldperiodString,foldcenter,enableRematch);
      fprintf(fp,"%d %d %d %f %f %f %f %f %f\n",AFLAGSMASK,qualitymask,seriesmask,reqstartdate,reqenddate,reqdimmag,reqbrightmag,averageRa,averageDec);
      fprintf(fp,"%d %d %d %d\n",reqXDown,reqYDown,reqXUp,reqYUp);
      fprintf(fp,"This file contains plot scaling coordinates\n");
      fprintf(fp,"Line 1: sx x0 x1 fx fx0\n");
      fprintf(fp,"Line 2: sy y0 y1 fy fy0\n");
      fprintf(fp,"Line 3: upper_limit foldflag foldperiod foldcenter enableRematch \n");
      fprintf(fp,"Line 4: AFLAGS qualitymask seriesmask reqstartdate reqenddate reqdimmag reqbrightmag averageRa averageDec\n");
      fprintf(fp,"Line 5: reqXDown reqYDown reqXUp reqYUp\n");
      fprintf(fp,"where sx,sy are the plot width and height in pixels\n");
      fprintf(fp,"      x0,x1 are the left and right graph borders \n");
      fprintf(fp,"      y0,y1 are the left and right graph borders \n");
      fprintf(fp,"      fx is pixels/(unit x value)\n");
      fprintf(fp,"      fy is pixels/(unit y value)\n");
      fprintf(fp,"      fx0 is the offset of xmin from the left border\n");
      fprintf(fp,"      fy0 is the offset of ymin from the y axis\n");
      for(i=0; i<n; i++){
        PFILESTARIMAGEEXT pFileStarImageExt;
        int magnitudeIndex;
        magnitudeIndex = crossindex[i];
        pFileStarImageExt = &db_table[magnitudeIndex];

        fprintf(fp,"%d x %f y %f %d %05d\n",i,x[i],y[i],pFileStarImageExt->filestarimage.seriesId,pFileStarImageExt->filestarimage.plateNumber);
      }
      fclose(fp);

    }
  }

  if(verbose)printf("c\n");
  gdImageRectangle(im_out,ix+x0,iy+y0,ix+x1,iy+y1,col1);
  if(style & PLOT_AVE){
    ya = iy+fy*(ym-fy0)+y0;
    gdImageLine(im_out,ix+x0,ya,ix+x1,ya,col1);
  }
  DX=pow(10.,(double)((int)(log10(dx/nt)+100)-100));
  if(dx / DX > 5*nt) DX *= 5;
  if(dx / DX > 2*nt) DX *= 2;
  DY=pow(10.,(double)((int)(log10(dy/nt)+100)-100));
  if(dy / DY > 5*nt) DY *= 5;
  if(dy / DY > 2*nt) DY *= 2;
  NDX = 1;
  if(dx / DX / NDX > 5) NDX *=2;
  if(dx / DX / NDX > 5) NDX *=2.5;
  if(dx / DX / NDX > 5) NDX *=2;
  NDY = 1;
  if(dy / DY / NDY > 6) NDY *=2;
  if(dy / DY / NDY > 6) NDY *=2.5;
  if(dy / DY / NDY > 6) NDY *=2;
  X0 = NDX*DX * floor(fx0/DX/NDX);
  Y0 = NDY*DY * floor(fy0/DY/NDY);
  for(v=X0,i=0; v < fx0+dx; v+=DX, i++){
    if(v < fx0)continue;
    xa = ix+fx*(v-fx0)+x0;
    if(i%NDX == 0){
      if(x1 - 5*(int)strlen(xlab)-4 > xa+5*(int)strlen(buf)/2){
        if(dx > 10) sprintf(buf,"%d",(int)v);
        else if(dx > 1) sprintf(buf,"%4.1f",v);
        else  sprintf(buf,"%4.3f",v);
        gdImageString(im_out,gdFontSmall,xa-5*(int)strlen(buf)/2,iy+y1+3, (unsigned char *) buf,col1);
      }
      len=6;
    }else len = 3;
    gdImageLine(im_out,xa,iy+y1,xa,iy+y1-len,col1);
    gdImageLine(im_out,xa,iy+y0,xa,iy+y0+len,col1);
  }
  for(v=Y0,i=0; v < fy0+dy; v+=DY,i++){
    if(v < fy0)continue;
    ya = fy*(v-fy0)+y0;

    if( ! (style & PLOT_FLIP) ) yy = iy+sy-ya;
    else yy=iy+ya;

    if(i%NDY == 0){
      sprintf(buf,"%4.2lf",v);
      gdImageString(im_out, gdFontSmall, ix+x0-6*strlen(buf)-4, yy-5, (unsigned char *) buf,col1);
      len=6;
    } else len = 3;
    gdImageLine(im_out,ix+x0,yy,ix+x0+len,yy,col1);
    gdImageLine(im_out,ix+x1,yy,ix+x1-len,yy,col1);
  }
  if(style & PLOT_LINE ){
    for(i=0; i<n-1; i++){
      xa = ix+fx*(x[i]-fx0)+x0;
      ya = iy+sy - (fy*(y[i]-fy0)+y0);
      xb = ix+fx*(x[i+1]-fx0)+x0;
      yb = iy+sy - (fy*(y[i+1]-fy0)+y0);
      if(col[i] == col[i+1])
        gdImageLine(im_out,xa,ya,xb,yb,colors[col[i]]);
    }
  } else /* if ((float)n/sx < 3.) */ {
    /* Plot all error bars first */
    for(i=0; i<n; i++){
      if(y[i] <= upper_limit || !(style & PLOT_UPPER_L)){
        if ((col[i] != COLOR_BLUE) &&
            (col[i] != COLOR_GRAY)) {
          /* Normal point */
          xa = ix+fx*(x[i]-fx0)+x0;
          ya = iy+(fy*(y[i]-fy0)+y0);
          {
            /* Now plot the error bars before the point itself */
            int errtop;
            int errbot;
            errtop = ya + (fy*err[i]);
            errbot = ya - (fy*err[i]);
            if (errtop > iy+y1) {
              errtop = iy+y1;
            }
            if (errbot < iy+y0) {
              errbot = iy+y0;
            }
            if (err[i] >= ERROR_BAR_CLIP) {
              gdImageDashedLine(im_out,xa,errtop,xa,errbot,colors[2]);
            } else {
              gdImageLine(im_out,xa,errtop,xa,errbot,colors[4]);
            }
          }
        }
        /* Now plot the limiting magnitudes */
        if (plotLimiting) {
          xa = ix+fx*(x[i]-fx0)+x0;
          ya = iy+(fy*(limiting[i]-fy0)+y0);
          if (col[i] != COLOR_GRAY) {
            /* Normal point - use a bar */
            if ((plotLimiting & QUALITY_LIMITING) != 0) {
              gdImageLine(im_out,xa,ya,xa+2,ya,grey);
              gdImageLine(im_out,xa,ya,xa-2,ya,grey);
            }
          } else {
            if ((plotLimiting & QUALITY_UNDETECTED) != 0) {
              colorCount[COLOR_BLUE]++;
              /* Limiting magnitude down arrow */
              gdImageLine(im_out,xa,ya,xa+2,ya-2,colors[COLOR_BLUE]);
              gdImageLine(im_out,xa,ya,xa-2,ya-2,colors[COLOR_BLUE]);
            }
          }
        }
      }else{
        xa = ix+fx*(x[i]-fx0)+x0;
        ya = iy+(fy*(upper_limit-fy0)+y0);
        //        gdImageLine(im_out,xa-2,ya-2,xa,ya,colors[col[i]]);
        //        gdImageLine(im_out,xa,ya,xa+2,ya-2,colors[col[i]]);
        gdImageLine(im_out,xa-2,ya-2,xa,ya,grey);
        gdImageLine(im_out,xa,ya,xa+2,ya-2,grey);
      }
    }


    /* Now plot symbols */
    for(i=0; i<n; i++){
      if (y[i] != 0.0) {
        if(y[i] <= upper_limit || !(style & PLOT_UPPER_L)){
          colorCount[col[i]]++;
          if (col[i] == COLOR_BLUE) {
            /* Too bright.  Use an up arrow */
            xa = ix+fx*(x[i]-fx0)+x0;
            ya = iy+(fy*(y[i]-fy0)+y0);
            gdImageLine(im_out,xa,ya,xa+2,ya+2,colors[col[i]]);
            gdImageLine(im_out,xa,ya,xa-2,ya+2,colors[col[i]]);
          } else {
            /* Normal point */
            xa = ix+fx*(x[i]-fx0)+x0;
            ya = iy+(fy*(y[i]-fy0)+y0);

            /* A box */
            gdImageLine(im_out,xa-1,ya-1,xa-1,ya+1,colors[col[i]]);
            gdImageLine(im_out,xa-1,ya+1,xa+1,ya+1,colors[col[i]]);
            gdImageLine(im_out,xa+1,ya+1,xa+1,ya-1,colors[col[i]]);
            gdImageLine(im_out,xa+1,ya-1,xa-1,ya-1,colors[col[i]]);
          }
        }else{
          xa = ix+fx*(x[i]-fx0)+x0;
          ya = iy+(fy*(upper_limit-fy0)+y0);
          //        gdImageLine(im_out,xa-2,ya-2,xa,ya,colors[col[i]]);
          //        gdImageLine(im_out,xa,ya,xa+2,ya-2,colors[col[i]]);
          gdImageLine(im_out,xa-2,ya-2,xa,ya,grey);
          gdImageLine(im_out,xa,ya,xa+2,ya-2,grey);
        }
      }
    }

  }

  if(verbose)printf("d\n");
  gdImageString(im_out, gdFontSmall, ix+sx/2-5*(int)strlen(title)/2,iy, (unsigned char *) title,col2);
  gdImageString(im_out, gdFontSmall, ix+sx/2-5*(int)strlen(lab1)/2, iy+10, (unsigned char *) lab1,col2);
  gdImageString(im_out, gdFontSmall, ix+sx/2-5*(int)strlen(lab2)/2, iy+20, (unsigned char *) lab2,col2);
  if (foldflag || enableRematch) {
    gdImageString(im_out, gdFontSmall, ix+sx/2-5*(int)strlen(lab3)/2, iy+30, (unsigned char *) lab3,col2);
  }
  gdImageString(im_out, gdFontSmall, ix+x0 - 5*strlen(ylab)-4, iy+y0-13, (unsigned char *) ylab,col2);
  gdImageString(im_out, gdFontSmall, ix+x1 - 5*strlen(xlab)-4, iy+y1+3, (unsigned char *) xlab,col2);

  return 0;
}

/*
 * AFLAGSMASK and qualitymask are the current settings
 * mode indicates the type of button
 *   mode = 0 leave unchanged
 *   mode = 1 set all flags to show everything
 *   mode = 2 clear AFLAGS and quality flags to show only good points
 *   mode = 3 toggle bit
 *   mode = 4 set the default
 * SETAFLAGS and setquality indicate the bit to be displayed
 *
 * The redisplay function invoked is in JavaScript in lightcurve_plot.php
 */

void AddCheckbox(char *REFurlencode,char *centerurlencode,double distance,char *source,char *rematchString,int AFLAGSMASK,int qualitymask,int seriesmask,int gsc_bin_index,char* tmpdir,int windowWidth,int windowHeight,int mode,int SETAFLAGS,int setquality,int setseries, char *reqWindow,char* foldperiodString,double foldcenter)
{
  int NEWAFLAGS = 0;
  int newquality = 0;
  int newseriesmask = 0;
  char setchecked[MAX_BUFFER];
  char outtext[MAX_BUFFER];
  char outtext2[MAX_BUFFER];
  char origtmpdir[MAX_BUFFER];
  int index;
  PFLAGSENTRY pFlagEntry;
  PQUALITYBIT pQualityEntry;
  char *charPtr;
  PSERIESLIST pSeriesList;
  int seriesIndex;
  if (mode == 1) {
    /* This is the set all mode */
    if ((AFLAGSMASK == 0x7fffffff) && (qualitymask == DEFAULT_QUALITY_MASK) && (seriesmask == allSeriesMask)) {
      /* Already at default: do nothing */
      return;
    } else {
      setchecked[0] = 0;
    }
    strcpy(outtext,"All points");
    NEWAFLAGS =  0x7fffffff;
    newquality = 0x7fffffff;
    newseriesmask = allSeriesMask;
  } else if (mode == 2) {
    /* This is the clear all mode */
    if ((AFLAGSMASK == 0) && ((qualitymask & (~(QUALITY_PATROL|QUALITY_NONPATROL))) == 0)) {
      /* Already at default: do nothing */
      return;
    } else {
      setchecked[0] = 0;
    }
    strcpy(outtext,"Unflagged points");
    NEWAFLAGS = 0;
    newquality = 0;
    newseriesmask = 0;
  } else if (mode == 4) {
    /* This is the default mode */
    if ((AFLAGSMASK == FILTER_AMASK_PLOT2) && (qualitymask == DEFAULT_QUALITY_MASK)) {
      /* Already at default: do nothing */
      return;
    }
    setchecked[0] = 0;
    strcpy(outtext,"Default points");
    NEWAFLAGS = FILTER_AMASK_PLOT2;
    newquality = DEFAULT_QUALITY_MASK;
    newseriesmask = allSeriesMask;
  } else if (mode == 3) {
    if (SETAFLAGS != 0) {
      if (SETAFLAGS == (GSC_CLASS_NONSTAR << GSC_CLASS_BIT)) {
        if ((SETAFLAGS & AFLAGSMASK) == SETAFLAGS) {
          strcpy(setchecked,"checked = \"checked\"");
        } else {
          setchecked[0] = 0;
        }
        NEWAFLAGS = SETAFLAGS;
        strcpy(outtext,"Non-star");
      } else {
        outtext[0] = 0;
        if ((SETAFLAGS & AFLAGSMASK) != 0) {
          strcpy(setchecked,"checked = \"checked\"");
          NEWAFLAGS = SETAFLAGS;
        } else {
          setchecked[0] = 0;
          NEWAFLAGS = SETAFLAGS;
        }
        for (index = 0; index < AflagsTableSize; index++) {
          pFlagEntry = &AflagsTable[index];
          if (SETAFLAGS == (1 << pFlagEntry->flagbit)) {
            strcpy(outtext,pFlagEntry->webname);
            break;
          }
        }
        if (outtext[0] == 0) {
          return;
        }
      }
    } else if (setquality != 0) {
      outtext[0] = 0;
      if ((setquality & qualitymask) != 0) {
        strcpy(setchecked,"checked = \"checked\"");
        newquality = setquality;
      } else {
        setchecked[0] = 0;
        newquality = setquality;
      }
      for (index = 0; index < webQualityTableSize; index++) {
        pQualityEntry = &webQualityMasks[index];
        if (setquality == pQualityEntry->qualityMask) {
          strcpy(outtext,pQualityEntry->qualityDescr);
          break;
        }
      }
      if (outtext[0] == 0) {
        return;
      }
    } else if (setseries != 0) {
      outtext[0] = 0;
      if ((setseries & seriesmask) != 0) {
        strcpy(setchecked,"checked = \"checked\"");
        newseriesmask = setseries;
      } else {
        setchecked[0] = 0;
        newseriesmask = setseries;
      }
      for (seriesIndex = 0; seriesIndex <= MAX_SERIES; seriesIndex++) {
        pSeriesList = &seriesList[seriesIndex];
        if (pSeriesList->bitMask == setseries) {
          strcpy(outtext,GetSeriesString(seriesIndex,1));
          strcat(outtext," series");
          break;
        }
      }
      if (outtext[0] == 0) {
        return;
      }


    } else {
      printf("ERROR: Neither SETAFLAGS nor setquality nor setseries is nonzero for mode %d\n",mode);
    }


  } else {
    printf("ERROR: Unrecognized mode %d\n",mode);
  }


  /* Escape the DASCH sign before we continue */

  charPtr = strstr(tmpdir,"tmp/");
  if (charPtr == NULL) {
    printf("ERROR: invalid format for tmpdir %s\n",tmpdir);
    exit(1);
  }
  charPtr += 4;
  strcpy(origtmpdir,charPtr);

  if (strlen(setchecked) != 0) {
    strcpy(outtext2,"Hide: ");
    strcpy(setchecked,"style = \"font-weight: bold;\"");
  } else {
    strcpy(outtext2,"Show: ");
  }
  strcat(outtext2,outtext);
  printf("<button  %s onClick=redisplay('%s','%s','%.0f','%s','%s','%d','%s','%d','%d','%d','%d','%d','%d','%d','%s','%s','%f'); return true;> %s</button>",
         setchecked,REFurlencode,centerurlencode,distance,source,rematchString,gsc_bin_index,origtmpdir,windowWidth,windowHeight,mode,NEWAFLAGS,newquality,newseriesmask,allSeriesMask,reqWindow,foldperiodString,foldcenter,outtext2);
  return;
}


int
PlateCompare(const void *first, const void *second)
{
  PFILESTARIMAGEEXT pPlateRecFirst = (PFILESTARIMAGEEXT)first;
  PFILESTARIMAGEEXT pPlateRecSecond = (PFILESTARIMAGEEXT)second;

  if (pPlateRecFirst->filestarimage.seriesId > pPlateRecSecond->filestarimage.seriesId) {
    return 1;
  } else if (pPlateRecFirst->filestarimage.seriesId < pPlateRecSecond->filestarimage.seriesId) {
    return -1;
  } else {

    if (pPlateRecFirst->filestarimage.plateNumber > pPlateRecSecond->filestarimage.plateNumber) {
      return 1;
    } else if (pPlateRecFirst->filestarimage.plateNumber < pPlateRecSecond->filestarimage.plateNumber) {
      return -1;
    } else {
      if (pPlateRecFirst->mosaicNumber > pPlateRecSecond->mosaicNumber) {
        return 1;
      } else if (pPlateRecFirst->mosaicNumber < pPlateRecSecond->mosaicNumber) {
        return -1;
      } else {
        if (pPlateRecFirst->filestarimage.solutionNumber > pPlateRecSecond->filestarimage.solutionNumber) {
          return 1;
        } else if (pPlateRecFirst->filestarimage.solutionNumber < pPlateRecSecond->filestarimage.solutionNumber) {
          return -1;
        } else {
          return 0;
        }
      }
    }
  }
}


int JulianDateCompare(const void *first, const void *second)
{
  double dateFirst = ((PFILESTARIMAGEEXT)first)->filestarimage.Date;
  double dateSecond = ((PFILESTARIMAGEEXT)second)->filestarimage.Date;
  if (dateFirst > dateSecond) {
  return 1;
  } else if (dateFirst < dateSecond) {
    return -1;
  } else {
    return 0;
  }
}


void AddLimitingMagnitudeEntry(PPLATELIMITINGREC pPlateLimitingRec,PFILESTARIMAGEEXT db_table,int *pNumMagnitudes,double centerRa,double centerDec,int RefType)
{
  int numMagnitudes = *pNumMagnitudes;
  PFILESTARIMAGEEXT pFileStarImageExt = &db_table[numMagnitudes];
  PFILESTARIMAGE pFileStarImage = &pFileStarImageExt->filestarimage;
  PPHOTPLATES pPhotPlates = &pFileStarImageExt->photplates;

  memset(pFileStarImageExt,0,sizeof(FILESTARIMAGEEXT));
  pFileStarImageExt->mosaicNumber = pPlateLimitingRec->mosaicNumber;
  pFileStarImageExt->quality = QUALITY_UNDETECTED;
  pFileStarImageExt->plateVersionId = pPlateLimitingRec->versionId;
  pFileStarImage->Date = pPlateLimitingRec->geoJulianDate;
  pFileStarImage->seriesId = pPlateLimitingRec->seriesId;
  pFileStarImage->plateNumber = pPlateLimitingRec->plateNumber;
  pFileStarImage->solutionNumber = pPlateLimitingRec->solutionNumber;
  pFileStarImage->limiting_mag_local = pPlateLimitingRec->limiting_mag_local;
  pFileStarImage->versionId = pPlateLimitingRec->versionId;
  pFileStarImage->timeAccuracy = -1;
  pPhotPlates->mosaicNumber = pPlateLimitingRec->mosaicNumber;
  pPhotPlates->quality = QUALITY_UNDETECTED;
  pPhotPlates->versionId = pPlateLimitingRec->versionId;

  if (RefType == REF_TYPE_NONE) {
    pFileStarImage->ra = centerRa;
    pFileStarImage->dec = centerDec;
  }

  numMagnitudes++;
  *pNumMagnitudes = numMagnitudes;

}


int main(int argc,char *argv[])
{
  int index;
  int nvals;
  char *argstr;
  char cmdchar;
  char REF[MAX_REF];          /* GSC2.3.2 reference number */
  char REFurlencode[4*MAX_REF];
  char *center = NULL;               /* Search center */
  char *centerurlencode = NULL;
  char tmpREF[MAX_REF];
  long long REFNumber;
  long long tmpREFNumber;
  int RefType;
  int tmpRefType;
  double ra;
  double dec;
  double radactual;
  double radminimum;
  long long REFNumberMinimum;
  int errorFlag = 0;
  char *tmpdir = NULL;
  int catalogNumber = 0;
  char source[MAX_BUFFER];
  char binaries[MAX_BUFFER];
  char lsbin[MAX_BUFFER];
  char catalogString[MAX_BUFFER];
  char qualifier[MAX_BUFFER];
  int AFLAGSMASK = -1;
  int qualitymask = -1;
  int seriesmask = -1;
  double distance = -1.0;
  int CURAFLAGSMASK = 0;
  int curqualitymask = 0;
  int gsc_bin_index = -1;
  int bitindex;
  int plotWidth = -1;
  int plotHeight = -1;
  int windowWidth = -1;
  int windowHeight = -1;
  char textfilename[MAX_BUFFER]; /* -t qualifier */
  char plotfilename[MAX_BUFFER]; /* -p qualifier */
  char storefilename[MAX_BUFFER]; /* -s qualifier */
  char fullplotname[MAX_BUFFER];
  char origtmpdir[MAX_BUFFER];
  char *charPtr;
  char db_name[MAX_FILENAME];
  char txt_name[MAX_FILENAME];
  char zipdb_name[MAX_FILENAME];
  char zipdb_name_gz[MAX_FILENAME];
  char gz_name[MAX_FILENAME];
  char vo_name[MAX_FILENAME];
  char plotshort_db_name[MAX_FILENAME];
  char plotshort_txt_name[MAX_FILENAME];
  char plotshort_zipdb_name[MAX_FILENAME];
  char plotshort_zipdb_name_gz[MAX_FILENAME];
  char plotshort_gz_name[MAX_FILENAME];

  char plotshort_vo_name[MAX_FILENAME];
  char plotshort_vo_name_gz[MAX_FILENAME];
  FILE *plotshort_dbhandle = NULL;
  FILE *plotshort_txthandle = NULL;
  int *plotshort_indices = NULL;
  int plotshort_index;
  int plotshort_ndata = 0;
  double plotshort_med;
  double plotshort_rms;
  int plotshort_statndata = 0; /* formerly statndata */

  char fullshort_db_name[MAX_FILENAME];
  char fullshort_txt_name[MAX_FILENAME];
  char fullshort_zipdb_name[MAX_FILENAME];
  char fullshort_zipdb_name_gz[MAX_FILENAME];
  char fullshort_gz_name[MAX_FILENAME];
  char fullshort_vo_name[MAX_FILENAME];
  char fullshort_vo_name_gz[MAX_FILENAME];
  FILE *fullshort_dbhandle = NULL;
  FILE *fullshort_txthandle = NULL;
  int *fullshort_indices = NULL;
  int fullshort_index;
  int fullshort_ndata = 0;

  double fullshort_med;
  double fullshort_rms;
  int fullshort_statndata = 0; /* formerly statndata */

  int *db_indices = NULL;
  int dbindex;
  int fullndata = 0;


  FILE *dbhandle = NULL;
  FILE *txthandle = NULL;
  File db_handle = NULL;
  TableHead db_header = NULL;
  PFILESTARIMAGEEXT db_table = NULL;
  PFILESTARIMAGEEXT pFileStarImageExt;
  size_t db_nrecs = 0;
  FILE *errorHandle;
  int printErrorFlag = 1;
  char error_name[MAX_FILENAME];

  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  PHOTGLOBAL basePhotGlobal;
  PPHOTGLOBAL pPhotGlobal = &basePhotGlobal;
  FILECOMMON fileCommon;
  PFILECOMMON pFileCommon0 = &fileCommon;
  FILECOMMON fileCommon1;
  PFILECOMMON pFileCommon1 = &fileCommon1;
  FILECOMMON fileCommon2;
  PFILECOMMON pFileCommon2 = &fileCommon2;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  int gotAnswer;
  int verbose = 0;
  int resetToOriginal = 0;
  int m44release = 0;
  int ndata = 0;
  int statndata2 = 0;
  double ramed = 0.0;
  double rarms;
  double decmed = 0.0;
  double decrms;
  double drad;
  double rawmed;
  double rawrms;
  char rstr[32], dstr[32];

  int *col = NULL;
  double *x = NULL;
  double *y = NULL;
  double *err = NULL;
  double *limiting = NULL;
  double *vector1 = NULL;
  double *vector2 = NULL;
  double *vector3 = NULL;
  double *vector1s = NULL;
  double *vector2s = NULL;
  double *vector3s = NULL;
  double *vector4s = NULL;
  double *vector1f = NULL;
  double *vector2f = NULL;
  double *vector3f = NULL;
  double *vector4f = NULL;
  int *crossindex = NULL;
  FILE *plotHandle = NULL;
  char *dotloc;
  int statResult;
  struct stat statbuf;
  int numMagnitudes;
  int numPlotted;
  int newNumMagnitudes;
  int magnitudeIndex;
  int curMagnitudeIndex;
  PFILESTARIMAGE pFileStarImage = NULL;
  PFILESTARIMAGE pCurFileStarImage = NULL;
  PPHOTSTARIMAGE pCurStarImage = NULL;
  int noneMagnitudeAlloc = 0;
  PFILESTARIMAGE pNoneMagnitudeTable = NULL;
  int histogramTable[MAX_NONE_HISTOGRAM];
  int groupNumber = 1;
  char cmdStr[MAX_BUFFER];
  int result;
  int imagesize = -1; /* Width and height of the extracted image in pixels */
  char *sizeunits = NULL;
  char defaultsizeunits[] = THUMBNAIL_PIXELS;
  PPHOTPLATES pPhotPlates = NULL;
  int dbCallocFlag = 0;
  char *listString = NULL;
  char *newListString = NULL;
  int listStringAlloc = 0;
  char tmpListString[MAX_MOSAICID_STRING];
  int listStringSize = 0;
  double actbrightmag = 99.0;
  double actdimmag = 0.0;
  double finalbrightmag = 99.0;
  double finaldimmag = 0.0;
  double finalbright;
  double finaldim;
  double reqbrightmag = 0;
  double reqdimmag = 0;
  double actstartdate = PIPELINE_MAX_DATE;
  double actenddate =   PIPELINE_MIN_DATE;
  double reqstartdate = 0;
  double reqenddate = 0;
  double dateyear;
  double tempmag;
  double templimiting;
  double actdimlimiting = 0.0;
  double actbrightlimiting = 99.0;
  int webquality; /* Augmented quality bits for web display */
  int reqXDown = -1; /* Upper left corner of new plot request */
  int reqYDown = -1;
  int reqXUp = -1;   /* Lower right corner of new plot request */
  int reqYUp = -1;
  int reqXValid = 0; /* Request is valid if nonzero - else rescale to max limits */
  int reqYValid = 0; /* Request is valid if nonzero - else rescale to max limits */
  int tempX;
  int tempY;
  char reqWindow[MAX_BUFFER];
  char reqWindow2[MAX_BUFFER];
  char starbase_title[MAX_BUFFER];
  char starbase_txt_title[MAX_BUFFER];
  PSORTTABLE pSortEntry;
  PSORTTABLE pSortTable = NULL;
  int sortTableSize = 0;
  GALAXYCOMMON galaxycommon;
  PGALAXYCOMMON pGalaxyCommon = &galaxycommon;
  int haveGalaxyTable = 1;
  int releaseField;

  /* Grzegorz Pomanski's plotting parameters */
  FILE *fp;
  int sx;
  int x0;
  int x1;
  double fx;
  double fx0;
  int sy;
  int y0;
  int y1;
  double fy;
  double fy0;
  char store[MAX_FILENAME] = "";
  char* regionFlag = NULL;
  double centerRa;
  double centerDec;
  double averageRa;
  double averageDec;
  double xxxd;
  int xxxi;
  PSERIESLIST pSeriesList;
  int seriesIndex;
  int curseriesmask = 0;
  int webseriesbit;
  PPHOTSTARIMAGE pMagnitudeTable = NULL;
  PFILESTARIMAGE pFileStarImageTable = NULL;
  PPHOTSTARIMAGE pMagnitudeTable1 = NULL;
  PPHOTSTARIMAGE pMagnitudeTable2 = NULL;

  int allocMagnitudes = 0;

  int goodPoints = 0;
  double error_bar_factor = 1.0;
  double rawerrrms = 0;
  double rawerrmed = 0;
  double flag;
  int reject_reason1;
  int reject_reason2;
  int quality = 0;
  int BFLAGSMASK = 0;
  int plateIndex;
  PPLATELIMITINGREC pPlateLimitingRec;
  int writeHeader = 1;
  int writetxtHeader = 1;
  double foldcenter = 0.0;
  double foldperiod = 0.0;
  char foldperiodString[MAX_PERIOD_STRING+1];
  double foldflag = 0;
  int foldfactor;
  int oldfoldflag = 0;
  int oldEnableRematch = 0;
  double oldfoldperiod = 0.0; /* Note: never use oldfoldperiod because it may be missing significant digits. */
  double oldfoldcenter = 0.0;
  double oldreqstartdate = 0.0;
  double oldreqenddate = 0.0;
  double oldreqbrightmag = 0.0;
  double oldreqdimmag = 0.0;
  int versionIdWarningFlag = 0;
  char AFLAGSBuffer[MAX_BITMAP_SIZE];
  char BFLAGSBuffer[MAX_BITMAP_SIZE];
  char qualityBuffer[MAX_BITMAP_SIZE];
  char AFLAGSbitsBuffer[MAX_BITMAP_SIZE];
  char BFLAGSbitsBuffer[MAX_BITMAP_SIZE];
  char qualitybitsBuffer[MAX_BITMAP_SIZE];
  int enableRematch = 0;
  char rematchString[10]; /* used for web redirection */

  memset(seriesList,0,sizeof(seriesList));
  memset(pGalaxyCommon,0,sizeof(GALAXYCOMMON));

  SetQueryCount(0);
  qualifier[0] = 0;
  source[0] = 0;
  catalogString[0] = 0;
  textfilename[0] = 0;
  plotfilename[0] = 0;
  storefilename[0] = 0;
  binaries[0] = 0;
  lsbin[0] = 0;

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


        case 'd': /* Temporary directory */
        case 'D':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            tmpdir = *++argv;
          }
          break;

        case 'e': /* Executable binary image directory */
        case 'E':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: binaries directory length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(binaries,*argv);
            }
          }
          break;

        case 'l': /* gzip binary image directory */
        case 'L':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: lsbin directory length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(lsbin,*argv);
            }
          }
          break;


        case 'u': /* Region Flag */
        case 'U':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            regionFlag = *++argv;
          }
          break;


        case 'f': /* AFLAGS */
        case 'F':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&AFLAGSMASK);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the AFLAGS %s\n",*argv);
              errorFlag = 1;
            }
            AFLAGSMASK |= (1 << GSC_VARIABLE_BIT);
          }
          break;

        case 'k': /* seriesmask */
        case 'K':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&seriesmask);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the seriesmask %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'y': /* quality */
        case 'Y':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&qualitymask);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the quality %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 's': /* Separation from the center source */
        case 'S':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%lf",&distance);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the distance %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;



        case 'g': /* gsc_bin_index */
        case 'G':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&gsc_bin_index);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the gsc_bin_index %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 'w': /* plot width */
        case 'W':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&plotWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the plotWidth %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;
        case 'h': /* plotHeight */
        case 'H':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&plotHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the plotHeight %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'a': /* window width */
        case 'A':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&windowWidth);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the windowWidth %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;
        case 'b': /* windowHeight */
        case 'B':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&windowHeight);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the windowHeight %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;


        case 'v': /* verbose */
          verbose += 1;
          break;

        case 'o': /* M44 Release */
          m44release += 1;
          break;

        case 'O':
          enableRematch = 1;
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
                strcpy(qualifier,*argv);
              }
            }
          }
          break;

        case 'r': /* Object designation */
        case 'R':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_REF-1)) {
              printf("ERROR: designation length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(REF,*argv);
              GetREFNumber(REF,&REFNumber,&RefType,1,1);
            }
          }
          break;


        case 'c': /* Search center */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            center = *argv;
          }
          break;

        case 'C': /* Fold center */
          argc--;
          if (argc < 1) {
            foldcenter = 0.0;
          } else {
            nvals = sscanf(*++argv,"%lf",&foldcenter);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the fold center %s\n",*argv);
              foldcenter = 0.0;
              foldflag = -1;
            }
            if ((foldcenter >= 1.0) &&
                (foldcenter < PIPELINE_MAX_DATE)) {
              foldcenter = ep2jd(foldcenter);
            }
          }
          break;

        case 'P': /* Fold period */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            nvals = sscanf(*++argv,"%lf",&foldperiod);
            strncpy(foldperiodString,*argv,MAX_PERIOD_STRING);
            foldperiodString[MAX_PERIOD_STRING-1] = 0;
            if (nvals != 1) {
              printf("ERROR: Unable to decode the fold period %s\n",*argv);
              foldperiod = 0.0;
              foldflag = -1;
            } else {
              if (foldperiod < 0) {
                printf("ERROR: Negative fold period %s\n",*argv);
                foldperiod = 0.0;
                foldflag = -1;
              } else {
                if (foldperiod == 0.0) {
                  foldflag = -1;
                } else {
                  if (foldflag >= 0) {
                    foldflag = 1;
                  }
                }
              }
            }
          }

          break;

        case 't': /* Text file name */
        case 'T':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: text file name length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(textfilename,*argv);
            }
          }
          break;

        case 'p': /* Plot file name */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: plot filename length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(plotfilename,*argv);
            }
          }
          break;

        case 'n': /* parameter store filename */
        case 'N':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) >= (MAX_BUFFER-1)) {
              printf("ERROR: parameter store filename length %zu for %s is too long\n",strlen(*argv),*argv);
              errorFlag = 1;
            } else {
              strcpy(storefilename,*argv);
            }
          }
          break;


        case 'i': /* Image width in pixels or arcsec */
        case 'I':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            nvals = sscanf(*++argv,"%d",&imagesize);
            if (nvals != 1) {
              printf("ERROR: Unable to decode the imagesize %s\n",*argv);
            }
          }
          break;

        case 'j': /* Units of image width: pixels or arcsec */
        case 'J':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            sizeunits = *++argv;
          }
          break;

        case 'm': /* Corners of new plot request */
        case 'M':
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for -%c\n",cmdchar);
          } else {
            nvals = sscanf(*++argv,"%d,%d,%d,%d",&reqXDown,&reqYDown,&reqXUp,&reqYUp);
            if (nvals != 4) {
              if (strcmp(*argv,"reset") == 0) {
                resetToOriginal = 1;
              } else {
                printf("ERROR: Unable to decode the requested plot dimensions %s\n",*argv);
              }
            }
          }
          break;


        default:
          printf("ERROR:  unknown command -%c\n",cmdchar);
          errorFlag = 1;

        }

      }

    }
  }

  /* Check validity of arguments */
  if (AFLAGSMASK < 0) {
    printf("ERROR: AFLAGS is not specified\n");
    errorFlag = 1;
  }
  if (seriesmask < 0) {
    printf("ERROR: seriesmask is not specified\n");
    errorFlag = 1;
  }
  if (qualitymask < 0) {
    printf("ERROR: quality is not specified\n");
    errorFlag = 1;
  }
  if (distance < 0) {
    printf("ERROR: distance is not specified\n");
    errorFlag = 1;
  }
  if (gsc_bin_index < 0) {
    printf("ERROR: gsc_bin_index is not specified\n");
    errorFlag = 1;
  }
  if (plotWidth < 0) {
    printf("ERROR: plotWidth is not specified\n");
    errorFlag = 1;
  }
  if (plotHeight < 0) {
    printf("ERROR: plotHeight is not specified\n");
    errorFlag = 1;
  }
  if (windowWidth < 0) {
    printf("ERROR: windowWidth is not specified\n");
    errorFlag = 1;
  }
  if (windowHeight < 0) {
    printf("ERROR: windowHeight is not specified\n");
    errorFlag = 1;
  }
  if (strlen(textfilename) == 0) {
    printf("ERROR: textfilename is not specified\n");
    errorFlag = 1;
  }
  if (strlen(plotfilename) == 0) {
    printf("ERROR: plotfilename is not specified\n");
    errorFlag = 1;
  }
  if (strlen(storefilename) == 0) {
    printf("ERROR: storefilename is not specified\n");
    errorFlag = 1;
  }
  if (center == NULL) {
    printf("ERROR: center is not specified\n");
    errorFlag = 1;
  }
  if (tmpdir == NULL) {
    printf("ERROR: plot directory is not specified\n");
    errorFlag = 1;
  }

  if (strlen(binaries) == 0) {
    printf("ERROR: binaries directory is not specified\n");
    errorFlag = 1;
  }

  if (strlen(lsbin) == 0) {
    printf("ERROR: lsbin directory is not specified\n");
    errorFlag = 1;
  }
  if (tmpdir != NULL) {
    if (strlen(textfilename)+strlen(tmpdir)+5 > MAX_BUFFER) {
      printf("ERROR: full name of the text file is too long \n");
      errorFlag = 1;
    }
    if (strlen(plotfilename)+strlen(tmpdir)+5 > MAX_BUFFER) {
      printf("ERROR: full name of the plot file is too long \n");
      errorFlag = 1;
    }
    if (strlen(storefilename)+strlen(tmpdir)+5 > MAX_BUFFER) {
      printf("ERROR: full name of the plot file is too long \n");
      errorFlag = 1;
    }
  }

  if (foldflag < 0) {
    foldperiod = 0.0;
    foldcenter = 0.0;
    foldflag = 0;
  }

  if (foldcenter < 1.0) {
    foldcenter = foldcenter*foldperiod;
  }

  if (imagesize < 0) {
    printf("ERROR: imagesize is not specified using %d\n",THUMBNAIL_IMAGESIZE);
    imagesize = THUMBNAIL_IMAGESIZE;
  }
  /* Validate requested graphing area */
  if (reqXDown < 0) {
    reqXDown = 0;
  }
  if (reqYDown < 0) {
    reqYDown = 0;
  }
  if (reqXUp < 0) {
    reqXUp = 0;
  }
  if (reqYUp < 0) {
    reqYUp = 0;
  }

  if (reqXDown > plotWidth) {
    reqXDown = plotWidth;
  }
  if (reqYDown > plotHeight) {
    reqYDown = plotHeight;
  }
  if (reqXUp > plotWidth) {
    reqXUp = plotWidth;
  }
  if (reqYUp > plotHeight) {
    reqYUp = plotHeight;
  }

  if (reqXDown > reqXUp) {
    tempX = reqXUp;
    reqXUp = reqXDown;
    reqXDown = tempX;
  }
  if (reqYDown > reqYUp) {
    tempY = reqYUp;
    reqYUp = reqYDown;
    reqYDown = tempY;
  }
  if (reqXDown == reqXUp) {
    if (reqXDown >= 1) {
      reqXDown--;
    } else {
      reqXUp++;
    }
  }
  if (reqYDown == reqYUp) {
    if (reqYDown >= 1) {
      reqYDown--;
    } else {
      reqYUp++;
    }
  }

  /* If we have drawn a graph, get the plot limits from the store file */
  if (tmpdir != NULL) {
    strcpy(store, tmpdir);
    strcat(store, "/");
    strcat(store, storefilename);

    if((fp = fopen(store,"r")) != NULL){
      fscanf(fp, "%d %d %d %lf %lf", &sx, &x0, &x1, &fx, &fx0);
      fscanf(fp, "%d %d %d %lf %lf", &sy, &y0, &y1, &fy, &fy0);
      fscanf(fp, "%lf %d %lf %lf %d", &xxxd, &oldfoldflag, &oldfoldperiod, &oldfoldcenter, &oldEnableRematch); /* NOTE: Never use oldfoldperiod, because it may be missing significant digits. */
      fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", &xxxi, &xxxi, &xxxi, &oldreqstartdate, &oldreqenddate, &oldreqdimmag, &oldreqbrightmag, &averageRa, &averageDec);
      fclose(fp);

      if (oldfoldflag == 0) {
        if ((reqXDown >= 0) && (reqXUp <= (plotWidth+1))) {
          reqXValid = 1;
          reqstartdate = ((reqXDown-x0)/fx) + fx0;
          reqenddate   = ((reqXUp  -x0)/fx) + fx0;
        }
        if ((reqYDown >= 0) && (reqYUp <= plotHeight)) {
          reqYValid = 1;
          reqbrightmag = ((reqYDown-y0)/fy) + fy0;
          reqdimmag  = ((reqYUp  -y0)/fy) + fy0;
        }
      }

      if (foldflag != 0) {
        /* For a lightcurve fold, ignore any swipes and use the previous selection */
        reqXValid = 1;
        reqstartdate = oldreqstartdate;
        reqenddate   = oldreqenddate;
        reqYValid = 1;
        reqbrightmag = oldreqbrightmag;
        reqdimmag  =   oldreqdimmag;
      }
    }
  }



  /* We are done with the requested plot size.  If the user checks on a box, then we need to rescale the Y-axis.  There is no need to rescale the X-axis, but make sure we
     display the full range */
  sprintf(reqWindow,"%d,%d,%d,%d",1,0,plotWidth-1,plotHeight);
  if (strlen(reqWindow) > (MAX_BUFFER-2)) {
    printf("ERROR: reqWindow buffer overflow\n");
    exit(1);
  }

  /* The following is for the reset function */
  strcpy(reqWindow2,"reset");

  if ((sizeunits == NULL) ||
      ((strcmp(sizeunits,THUMBNAIL_PIXELS) != 0) &&
       (strcmp(sizeunits,THUMBNAIL_ARCSEC) != 0))) {
    sizeunits = defaultsizeunits;
    printf("ERROR: sizeunits is not specified, using %s\n",sizeunits);
  }

  if (tmpdir != NULL) {
    strcpy(fullplotname,tmpdir);
    strcat(fullplotname,"/");
    strcat(fullplotname,plotfilename);
    plotHandle = fopen(fullplotname,"w");
    if (plotHandle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to open the output file %s\n",fullplotname);
    } else {
      if (verbose) {
        printf("Output file %s\n",fullplotname);
      }
    }
  }

  if (errorFlag == 1) {
    printf("Usage: web_plot -r <REF> -q <catalog> -g <gsc_bin_index> -d <directory> -t <text name> -p <plot name> [-v]\n");
    printf("       where -v is the verbose flag\n");
    printf("             -r is the object reference\n");
    printf("             -w is the plot width in pixels\n");
    printf("             -h is the plot height in pixels\n");
    printf("             -d is the plot and text directory\n");
    printf("             -t is the text file name\n");
    printf("             -p is the plot file name\n");
    printf("             -n is the parameter store filename\n");
    printf("             -a is the window width in pixels\n");
    printf("             -b is the window height in pixels\n");
    printf("             -e is the directory containing the executable binaries\n");
    printf("             -l is the directory containing the gzip binary\n");
    printf("             -s is the separation from the search center in arcsec\n");
    printf("             -c is the search center\n");
    printf("             -i is the width of the extracted image in units specified by -j below\n");
    printf("             -j units of width of the extracted images - pixels or arcsec\n");
    printf("             -m request plot size <reqXDown,reqYDown,reqXUp,reqYUp>\n");
    printf("             -y plate quality mask\n");
    printf("             -o M44 release\n");
    printf("             -O Optimize location of transients\n");
    printf("             -f AFLAGS\n");
    printf("             -k seriesmask\n");
    printf("             -u authorization region\n");
    printf("             -C fold center\n");
    printf("             -P foldperiod\n");

    exit(1);
  }
  /* Check our authorization */

  if ((GetBinCenter(pGscBin,gsc_bin_index,&centerRa,&centerDec,"web_plot") != 0) ||
      (CheckAuthorization(regionFlag,centerRa,centerDec,&releaseField) == 0)) {
    /* Authorization failure */
    printf("Sorry, these data are not available due to <a target=\"_blank\" href=\"https://dasch.cfa.harvard.edu/data-access/#restrictions\">temporary data access restrictions</a><br/>\n");
    exit(0);
  }

  if ((qualitymask & DATABASE_QUALITY_MASK) != 0) {
    AFLAGSMASK |= (1<<FILTER_AFLAG_QUALITY);
  } else {
    AFLAGSMASK &= ~(1<<FILTER_AFLAG_QUALITY);
  }

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pConnection);
  dasch_init_photdb(pPhotConnection);

  InitSeriesTable(pConnection,pPhotConnection);
  gotAnswer = GetPhotometryGlobal(pPhotConnection,pPhotGlobal);
  if (gotAnswer != 1)  {
    printf("ERROR: failed to get the global photometry table\n");
    exit(1);
  }

  if (pPhotGlobal->magnitudeFile != PHOT_MAGNITUDEFILE_YES) {
    printf("ERROR: database uses obsolute magnitude table\n");
    exit(1);
  }

  InitFileCommon(stdout,pFileCommon0,pConnection,pPhotConnection,catalogNumber,0);
  InitMaxPlateNumber(pConnection,pFileCommon0->maxPlateNumber);

  if (enableRematch == 0) {
    rematchString[0] = 0;
  } else {
    strcpy(rematchString,"checked");
    pFileCommon0->enableRematch = 1;
    pFileCommon1->enableRematch = 1;
    pFileCommon2->enableRematch = 1;
  }

  /* Now see if we have already written out a result file for this object */
  strcpy(db_name,tmpdir);
  strcat(db_name,"/");
  strcat(db_name,textfilename);

  strcpy(plotshort_db_name,tmpdir);
  strcat(plotshort_db_name,"/plotshort_");
  strcat(plotshort_db_name,textfilename);

  strcpy(fullshort_db_name,tmpdir);
  strcat(fullshort_db_name,"/fullshort_");
  strcat(fullshort_db_name,textfilename);

  dotloc = strrchr(db_name,'.');
  if (dotloc) {
    *dotloc = 0;
    strcpy(vo_name,db_name);
    strcpy(zipdb_name,db_name);
    strcpy(txt_name,db_name);
    strcat(db_name,".db");
    strcat(txt_name,".txt");
    strcat(vo_name,".xml");
    strcat(zipdb_name,"_zip.db");
    strcpy(zipdb_name_gz,zipdb_name);
    strcat(zipdb_name_gz,".gz");
    strcpy(gz_name,db_name);
    strcat(gz_name,".gz");
  } else {
    printf("ERROR: %s has no suffix\n",textfilename);
    exit(1);
  }

  dotloc = strrchr(plotshort_db_name,'.');
  if (dotloc) {
    *dotloc = 0;
    strcpy(plotshort_vo_name,plotshort_db_name);
    strcpy(plotshort_zipdb_name,plotshort_db_name);
    strcpy(plotshort_txt_name,plotshort_db_name);
    strcat(plotshort_db_name,".db");
    strcat(plotshort_txt_name,".txt");
    strcat(plotshort_vo_name,".xml");
    strcpy(plotshort_vo_name_gz,plotshort_vo_name);
    strcat(plotshort_vo_name_gz,".gz");
    strcat(plotshort_zipdb_name,"_zip.db");
    strcpy(plotshort_zipdb_name_gz,plotshort_zipdb_name);
    strcat(plotshort_zipdb_name_gz,".gz");
    strcpy(plotshort_gz_name,plotshort_db_name);
    strcat(plotshort_gz_name,".gz");
  } else {
    printf("ERROR: %s has no suffix\n",textfilename);
    exit(1);
  }
  /* Clear out files from a previous invocation */
  statResult = stat(plotshort_db_name,&statbuf);
  if (statResult == 0) {
    unlink(plotshort_db_name);
  }
  statResult = stat(plotshort_txt_name,&statbuf);
  if (statResult == 0) {
    unlink(plotshort_txt_name);
  }
  statResult = stat(plotshort_zipdb_name,&statbuf);
  if (statResult == 0) {
    unlink(plotshort_zipdb_name);
  }
  statResult = stat(plotshort_zipdb_name_gz,&statbuf);
  if (statResult == 0) {
    unlink(plotshort_zipdb_name_gz);
  }
  statResult = stat(plotshort_gz_name,&statbuf);
  if (statResult == 0) {
    unlink(plotshort_gz_name);
  }
  statResult = stat(plotshort_vo_name,&statbuf);
  if (statResult == 0) {
    unlink(plotshort_vo_name);
  }
  statResult = stat(plotshort_vo_name_gz,&statbuf);
  if (statResult == 0) {
    unlink(plotshort_vo_name_gz);
  }

  dotloc = strrchr(fullshort_db_name,'.');
  if (dotloc) {
    *dotloc = 0;
    strcpy(fullshort_vo_name,fullshort_db_name);
    strcpy(fullshort_zipdb_name,fullshort_db_name);
    strcpy(fullshort_txt_name,fullshort_db_name);
    strcat(fullshort_db_name,".db");
    strcat(fullshort_txt_name,".txt");
    strcat(fullshort_vo_name,".xml");
    strcpy(fullshort_vo_name_gz,fullshort_vo_name);
    strcat(fullshort_vo_name_gz,".gz");
    strcat(fullshort_zipdb_name,"_zip.db");
    strcpy(fullshort_zipdb_name_gz,fullshort_zipdb_name);
    strcat(fullshort_zipdb_name_gz,".gz");
    strcpy(fullshort_gz_name,fullshort_db_name);
    strcat(fullshort_gz_name,".gz");

  } else {
    printf("ERROR: %s has no suffix\n",textfilename);
    exit(1);
  }
  /* Clear out files from a previous invocation */
  statResult = stat(fullshort_db_name,&statbuf);
  if (statResult == 0) {
    unlink(fullshort_db_name);
  }
  statResult = stat(fullshort_txt_name,&statbuf);
  if (statResult == 0) {
    unlink(fullshort_txt_name);
  }
  statResult = stat(fullshort_zipdb_name,&statbuf);
  if (statResult == 0) {
    unlink(fullshort_zipdb_name);
  }
  statResult = stat(fullshort_zipdb_name_gz,&statbuf);
  if (statResult == 0) {
    unlink(fullshort_zipdb_name_gz);
  }
  statResult = stat(fullshort_gz_name,&statbuf);
  if (statResult == 0) {
    unlink(fullshort_gz_name);
  }
  statResult = stat(fullshort_vo_name,&statbuf);
  if (statResult == 0) {
    unlink(fullshort_vo_name);
  }
  statResult = stat(fullshort_vo_name_gz,&statbuf);
  if (statResult == 0) {
    unlink(fullshort_vo_name_gz);
  }
  if (RefType == REF_TYPE_NONE) {
    sprintf(starbase_title," Limiting Magnitude Data near %s using the %s calibration catalog",center,catalogText[catalogNumber]);
  } else {
    sprintf(starbase_title," Photometry data for %s from search centered on %s using the %s calibration catalog",REF,center,catalogText[catalogNumber]);
  }

  /* Hack here!  Adding an ASCII 0x7 = "BEL" = "\a" turns the file into a binary file in the eyes of firefox and chrome */
  strcpy(starbase_txt_title,starbase_title);
  strcat(starbase_txt_title,"\a");

  statResult = stat(db_name,&statbuf);
  if (statResult != 0) {
    if (OpenGalaxyFiles(pGalaxyCommon, qualifier, NULL, GetPhotFileBase(catalogString)) != 0) {
      haveGalaxyTable = 0;
    }

    if (haveGalaxyTable) {
      haveGalaxyTable = LoadGalaxyTable(pGalaxyCommon,centerRa,centerDec);
    }

    /* Find all of the images in this bin */
    if ((REFNumber == 0)  && (RefType == REF_TYPE_NONE)) {
      numMagnitudes = 0;
    } else {
      LocateNoneImages(pGscBin,pFileCommon0,gsc_bin_index,&pMagnitudeTable,&allocMagnitudes,&numMagnitudes,catalogString,1,0);
    }

    if (numMagnitudes > 0) {
      if (numMagnitudes >= noneMagnitudeAlloc) {
        noneMagnitudeAlloc = numMagnitudes + 1000;
        if (pNoneMagnitudeTable != NULL) {
          free(pNoneMagnitudeTable);
        }
        pNoneMagnitudeTable = (PFILESTARIMAGE)calloc(noneMagnitudeAlloc,sizeof(FILESTARIMAGE));
        if (pNoneMagnitudeTable == NULL) {
          printf("ERROR: failed to allocate pNoneMagnitudeTable of size %d\n",noneMagnitudeAlloc);
          exit(1);
        }
      }
      for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
        pCurStarImage = &pMagnitudeTable[magnitudeIndex];
        pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
        memcpy(pFileStarImage,pCurStarImage->pFileStarImage,sizeof(FILESTARIMAGE));
      }
    } else {
      if ((REFNumber != 0)  || (RefType != REF_TYPE_NONE)) {
        printf("ERROR: no magnitudes found at gsc bin %d\n",gsc_bin_index);
        exit(1);
      }
    }

    /* Now search through the results and extract everything that applies to our reference number */
    ProcessNoneImagesX(
      pGscBin,
      pFileCommon0,
      NULL,
      NULL,
      NULL,
      pNoneMagnitudeTable,
      numMagnitudes,
      verbose,
      histogramTable,
      &groupNumber,
      gsc_bin_index
    );

    if (RefType == REF_TYPE_DASCH) {
      /* This is a DASCH reference number.  We need to perform astrometric checks to get everything in the vicinity that overlaps */
      if (GetDASCHCoordinates(REF,&ra,&dec,1,REF_TYPE_DASCH) != 0) {
        exit(1);
      }
      /* Now assign DASCH numbers to everything in the magnitude table */
      radminimum = 90.0;
      REFNumberMinimum = 0;
      /* Now find the group which is closest to our desired location */
      for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
        pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];

        GetREF(pFileStarImage->REFNumber,tmpREF,1,1);
        GetREFNumber(tmpREF,&tmpREFNumber,&tmpRefType,1,1);

        if (pFileStarImage->REFNumber == REFNumber) {
          /* We found the one we were looking for, so just exit */
          REFNumberMinimum = REFNumber;
          break;
        }
        if (tmpRefType == REF_TYPE_DASCH) {
          radactual = wcsdist(ra,dec,pFileStarImage->ra,pFileStarImage->dec);
          if (radactual < radminimum) {
            radminimum = radactual;
            REFNumberMinimum = pFileStarImage->REFNumber;
          }
        }
      }
      if (REFNumberMinimum == 0) {
        printf("ERROR: No objects found for %s\n",REF);
        exit(1);
      }
    } else {
      /* Non DASCH object. Just use its reference number */
      REFNumberMinimum = REFNumber;
    }
    /* Now extract all objects with our desired reference number */
    curMagnitudeIndex = 0;
    db_table = calloc(numMagnitudes+pGalaxyCommon->plateCount,sizeof(FILESTARIMAGEEXT));
    dbCallocFlag = 1;
    if (db_table == NULL) {
      printf("ERROR: failed to allocate db_table of size %d\n",numMagnitudes);
      exit(1);
    }
    for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
      pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
      if (pFileStarImage->REFNumber == REFNumberMinimum) {
        if (curMagnitudeIndex != magnitudeIndex) {
          pCurFileStarImage = &pNoneMagnitudeTable[curMagnitudeIndex];
          memcpy(pCurFileStarImage,pFileStarImage,sizeof(FILESTARIMAGE));
        }
        curMagnitudeIndex++;

      }
    }
    numMagnitudes = curMagnitudeIndex;
    if ((REFNumber != 0)  || (RefType != REF_TYPE_NONE)) {
      if (numMagnitudes == 0) {
        printf("ERROR: No images found for %s\n",REF);
        exit(1);
      }
    }


    /* Write the image file out */
    dbhandle = fopen(db_name,"wt");
    if (dbhandle == NULL) {
      printf("ERROR: failed (1) to open %s for writing\n",db_name);
      exit(1);
    }
    txthandle = fopen(txt_name,"wt");
    if (txthandle == NULL) {
      printf("ERROR: failed (1) to open %s for writing\n",txt_name);
      exit(1);
    }


    for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
      pFileStarImage = &pNoneMagnitudeTable[magnitudeIndex];
      pFileStarImageExt = &db_table[magnitudeIndex];
      memcpy(&pFileStarImageExt->filestarimage,pFileStarImage,sizeof(FILESTARIMAGE));
      pPhotPlates = &pFileStarImageExt->photplates;

      if (GetPhotPlate(pPhotConnection,
                       GetSeriesString(pFileStarImage->seriesId,1),
                       pFileStarImage->plateNumber,
                       pPhotPlates,
                       catalogString,1) != 0) {
        printf("ERROR: Failed to get photplates record for %s%05d\n",GetSeriesString(pFileStarImage->seriesId,0),pFileStarImage->plateNumber);
        pPhotPlates->mosaicNumber=99;
        pPhotPlates->quality = QUALITY_UNINITIALIZED;
        pPhotPlates->versionId = 0;
      }

      if (pPhotPlates->versionId == 0) {
        pPhotPlates->mosaicNumber=99;
        pPhotPlates->quality = QUALITY_UNINITIALIZED;
        pPhotPlates->versionId = 0;

      }

      if (pPhotPlates->versionId != pFileStarImage->versionId) {
        printf("*");

        if (printErrorFlag) {
          printErrorFlag = 0;
          strcpy(error_name,tmpdir);
          charPtr = strstr(error_name,"tmp/");
          if (charPtr) {
            charPtr += 4;
            *charPtr = 0;

            strcat(error_name,"dasch.log");
            errorHandle = fopen(error_name,"at");
            if (errorHandle != NULL) {
              fprintf(errorHandle,"ERROR: versionId mismatch %d %d for %s%05d gsc_bin_index %d\n",pPhotPlates->versionId,pFileStarImage->versionId,GetSeriesString(pFileStarImage->seriesId,0),pFileStarImage->plateNumber,pFileStarImage->gsc_bin_index);
              fclose(errorHandle);
            }
          } else {
            printf("ERROR: illegal format for directory %s\n",tmpdir);
          }
        }


        pPhotPlates->mosaicNumber=99;
        pPhotPlates->quality = QUALITY_UNINITIALIZED;
        pPhotPlates->versionId = 0;
      }
      GetFullQuality(&pFileStarImageExt->filestarimage,pPhotPlates->quality,&quality);
      pPhotPlates->quality = quality;
      pFileStarImageExt->mosaicNumber = pPhotPlates->mosaicNumber;
      pFileStarImageExt->quality = pPhotPlates->quality;
      pFileStarImageExt->plateVersionId = pPhotPlates->versionId;

      averageRa += pFileStarImageExt->filestarimage.ra;
      averageDec += pFileStarImageExt->filestarimage.dec;
    }
    averageRa = averageRa/(1.0*numMagnitudes);
    averageDec = averageDec/(1.0*numMagnitudes);
    /* Sort the results by series and plateNumber so that we can add limiting magnitude records.  Limiting
       magnitude records are already sorted by seriesId and plateNumber */
    qsort(db_table,numMagnitudes,sizeof(FILESTARIMAGEEXT),PlateCompare);
    int magnitudeIndex = 0;
    pFileStarImageExt =  &db_table[magnitudeIndex];
    plateIndex = 0;
    newNumMagnitudes = numMagnitudes;
    while (plateIndex < pGalaxyCommon->plateCount) {

      pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];

      if (pPlateLimitingRec->seriesId == pFileStarImageExt->filestarimage.seriesId) {
        if (pPlateLimitingRec->plateNumber == pFileStarImageExt->filestarimage.plateNumber) {
          if (( versionIdWarningFlag == 0) && (pPlateLimitingRec->versionId != pFileStarImageExt->filestarimage.versionId)) {
            versionIdWarningFlag = 1;
            printf("WARNING: versionId mismatch %d %d for plate %s%05d mosaicNumber %d %d gsc_bin_index %d, line %d\n",pPlateLimitingRec->versionId, pFileStarImageExt->filestarimage.versionId,GetSeriesString(pFileStarImageExt->filestarimage.seriesId,0),pFileStarImageExt->filestarimage.plateNumber,pPlateLimitingRec->mosaicNumber, pFileStarImageExt->mosaicNumber,pFileStarImageExt->filestarimage.gsc_bin_index,__LINE__);
          }

          if (pPlateLimitingRec->mosaicNumber == pFileStarImageExt->mosaicNumber) {
            if (pPlateLimitingRec->solutionNumber == pFileStarImageExt->filestarimage.solutionNumber) {
              /* Have a match, so we already know the limiting magnitude */
              plateIndex++;
              if (plateIndex >= pGalaxyCommon->plateCount) {
                break;
              }
              pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
              magnitudeIndex++;
              if (magnitudeIndex >= numMagnitudes) {
                break;
              }
              pFileStarImageExt = &db_table[magnitudeIndex];
            } else  if (pPlateLimitingRec->solutionNumber < pFileStarImageExt->filestarimage.solutionNumber) {
              AddLimitingMagnitudeEntry(pPlateLimitingRec,db_table,&newNumMagnitudes,centerRa,centerDec,RefType);
              plateIndex++;
              if (plateIndex >= pGalaxyCommon->plateCount) {
                break;
              }
              pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
            } else {
              magnitudeIndex++;
              if (magnitudeIndex >= numMagnitudes) {
                break;
              }
              pFileStarImageExt = &db_table[magnitudeIndex];
            }
          } else if (pPlateLimitingRec->mosaicNumber < pFileStarImageExt->mosaicNumber) {
            AddLimitingMagnitudeEntry(pPlateLimitingRec,db_table,&newNumMagnitudes,centerRa,centerDec,RefType);
            plateIndex++;
            if (plateIndex >= pGalaxyCommon->plateCount) {
              break;
            }
            pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
          } else {
            magnitudeIndex++;
            if (magnitudeIndex >= numMagnitudes) {
              break;
            }
            pFileStarImageExt = &db_table[magnitudeIndex];
          }

        } else if (pPlateLimitingRec->plateNumber < pFileStarImageExt->filestarimage.plateNumber) {
          AddLimitingMagnitudeEntry(pPlateLimitingRec,db_table,&newNumMagnitudes,centerRa,centerDec,RefType);
          plateIndex++;
          if (plateIndex >= pGalaxyCommon->plateCount) {
            break;
          }
          pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
        } else {
          magnitudeIndex++;
          if (magnitudeIndex >= numMagnitudes) {
            break;
          }
          pFileStarImageExt = &db_table[magnitudeIndex];
        }

      } else if (pPlateLimitingRec->seriesId < pFileStarImageExt->filestarimage.seriesId) {
        AddLimitingMagnitudeEntry(pPlateLimitingRec,db_table,&newNumMagnitudes,centerRa,centerDec,RefType);
        plateIndex++;
        if (plateIndex >= pGalaxyCommon->plateCount) {
          break;
        }
        pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
      } else {
        magnitudeIndex++;
        if (magnitudeIndex >= numMagnitudes) {
          break;
        }
        pFileStarImageExt = &db_table[magnitudeIndex];
      }
    }
    while (plateIndex < pGalaxyCommon->plateCount) {
      pPlateLimitingRec = &pGalaxyCommon->plateLimitingBuffer[plateIndex];
      AddLimitingMagnitudeEntry(pPlateLimitingRec,db_table,&newNumMagnitudes,centerRa,centerDec,RefType);
      plateIndex++;
    }
    numMagnitudes = newNumMagnitudes;
    /* Sort our results by Julian Date */
    qsort(db_table,numMagnitudes,sizeof(FILESTARIMAGEEXT),JulianDateCompare);
    for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
      pFileStarImageExt = &db_table[magnitudeIndex];
      pFileStarImage = &pFileStarImageExt->filestarimage;
      pPhotPlates = &pFileStarImageExt->photplates;

      if (pFileStarImageExt->mosaicNumber != 99) {
        WriteStarbaseRecord(pFileStarImage,pPhotPlates,dbhandle,starbase_title,&writeHeader,catalogNumber);
        WriteStarbaseRecord(pFileStarImage,pPhotPlates,txthandle,starbase_txt_title,&writetxtHeader,catalogNumber);
      }
    }

    fclose(dbhandle);
    fclose(txthandle);

    sprintf(cmdStr,"%s/cp %s %s\n",lsbin,db_name,zipdb_name);
    result = system(cmdStr);
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,db_name);
    result = system(cmdStr);
    sprintf(cmdStr,"%s/mv %s %s\n",lsbin,zipdb_name,db_name);
    result = system(cmdStr);

    /* Now write a VOTable to vo_name */

    sprintf(cmdStr,"%s/votable -i %s -o %s\n",binaries,db_name,vo_name);

    result = system(cmdStr);
    if (result != 0) {
      printf("ERROR: result %d creating the votable\n",result);
    }
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,vo_name);
    result = system(cmdStr);
  } else {

    /* Here we have already prepared the data.  Read it in */
    db_handle = Open(db_name,"r");
    if (db_handle == NULL) {
      errorFlag = 1;
      printf("ERROR: Failed to find the db file %s\n",db_name);
    } else {
      if (verbose) {
        printf("Found db file %s\n",db_name);
      }
    }


    db_table = table_loadva(db_handle,
                            &db_header,
                            NULL, /* hbase */
                            NULL, /* rows */
                            NULL,
                            sizeof(FILESTARIMAGEEXT),
                            &db_nrecs,
                            TblDbl,"X_IMAGE",TblOff(PFILESTARIMAGEEXT,filestarimage.X_IMAGE),
                            TblDbl,"Y_IMAGE",TblOff(PFILESTARIMAGEEXT,filestarimage.Y_IMAGE),
                            TblDbl,"MAG_ISO",TblOff(PFILESTARIMAGEEXT,filestarimage.MAG_ISO),
                            TblDbl,"ra",TblOff(PFILESTARIMAGEEXT,filestarimage.ra),
                            TblDbl,"dec",TblOff(PFILESTARIMAGEEXT,filestarimage.dec),
                            TblDbl,"Date",TblOff(PFILESTARIMAGEEXT,filestarimage.Date),
                            TblDbl,"FLUX_ISO",TblOff(PFILESTARIMAGEEXT,filestarimage.FLUX_ISO),
                            TblDbl,"MAG_APER",TblOff(PFILESTARIMAGEEXT,filestarimage.MAG_APER),
                            TblDbl,"MAG_AUTO",TblOff(PFILESTARIMAGEEXT,filestarimage.MAG_AUTO),
                            TblDbl,"KRON_RADIUS",TblOff(PFILESTARIMAGEEXT,filestarimage.KRON_RADIUS),
                            TblDbl,"BACKGROUND",TblOff(PFILESTARIMAGEEXT,filestarimage.BACKGROUND),
                            TblDbl,"FLUX_MAX",TblOff(PFILESTARIMAGEEXT,filestarimage.FLUX_MAX),
                            TblDbl,"THETA_J2000",TblOff(PFILESTARIMAGEEXT,filestarimage.THETA_J2000),
                            TblDbl,"ELLIPTICITY",TblOff(PFILESTARIMAGEEXT,filestarimage.ELLIPTICITY),
                            TblDbl,"ISOAREA_WORLD",TblOff(PFILESTARIMAGEEXT,filestarimage.ISOAREA_WORLD),
                            TblDbl,"FWHM_IMAGE",TblOff(PFILESTARIMAGEEXT,filestarimage.FWHM_IMAGE),
                            TblDbl,"FWHM_WORLD",TblOff(PFILESTARIMAGEEXT,filestarimage.FWHM_WORLD),
                            TblDbl,"plate_dist",TblOff(PFILESTARIMAGEEXT,filestarimage.plate_dist),
                            TblDbl,"Blendedmag",TblOff(PFILESTARIMAGEEXT,filestarimage.Blendedmag),
                            TblDbl,"dradRMS2",TblOff(PFILESTARIMAGEEXT,filestarimage.dradRMS2),
                            TblDbl,"ra_2",TblOff(PFILESTARIMAGEEXT,filestarimage.ra_2),
                            TblDbl,"dec_2",TblOff(PFILESTARIMAGEEXT,filestarimage.dec_2),

                            TblFlt,"magcal_iso",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_iso),
                            TblFlt,"magcal_iso_rms",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_iso_rms),
                            TblFlt,"magcal_local",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_local),
                            TblFlt,"magcal_local_rms",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_local_rms),
                            TblFlt,"limiting_mag_local",TblOff(PFILESTARIMAGEEXT,filestarimage.limiting_mag_local),
                            TblFlt,"magcal_local_error",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_local_error),
                            TblFlt,"magcor_local",TblOff(PFILESTARIMAGEEXT,filestarimage.magcor_local),
                            TblFlt,"extinction",TblOff(PFILESTARIMAGEEXT,filestarimage.extinction),
                            TblFlt,"magcal_magdep",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_magdep),
                            TblFlt,"magcal_magdep_rms",TblOff(PFILESTARIMAGEEXT,filestarimage.magcal_magdep_rms),
                            TblFlt,"RaPM",TblOff(PFILESTARIMAGEEXT,filestarimage.RaPM),
                            TblFlt,"DecPM",TblOff(PFILESTARIMAGEEXT,filestarimage.DecPM),
                            TblInt,"A2FLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.A2FLAGS),
                            TblInt,"B2FLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.B2FLAGS),
                            TblFlt,"timeAccuracy",TblOff(PFILESTARIMAGEEXT,filestarimage.timeAccuracy),
                            TblInt,"maskIndex",TblOff(PFILESTARIMAGEEXT,filestarimage.maskIndex),

                            TblInt,"gsc_bin_index",TblOff(PFILESTARIMAGEEXT,filestarimage.gsc_bin_index),
                            TblInt,"plateNumber",TblOff(PFILESTARIMAGEEXT,filestarimage.plateNumber),
                            TblInt,"NUMBER",TblOff(PFILESTARIMAGEEXT,filestarimage.NUMBER),
                            TblInt,"versionId",TblOff(PFILESTARIMAGEEXT,filestarimage.versionId),
                            TblInt,"AFLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.AFLAGS),
                            TblInt,"BFLAGS",TblOff(PFILESTARIMAGEEXT,filestarimage.BFLAGS),
                            TblInt,"ISO0",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO0),
                            TblInt,"ISO1",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO1),
                            TblInt,"ISO2",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO2),
                            TblInt,"ISO3",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO3),
                            TblInt,"ISO4",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO4),
                            TblInt,"ISO5",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO5),
                            TblInt,"ISO6",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO6),
                            TblInt,"ISO7",TblOff(PFILESTARIMAGEEXT,filestarimage.ISO7),
                            TblInt,"npoints_local",TblOff(PFILESTARIMAGEEXT,filestarimage.npoints_local),
                            TblInt,"rejectFlag",TblOff(PFILESTARIMAGEEXT,filestarimage.rejectFlag),
                            TblSht,"local_bin_index",TblOff(PFILESTARIMAGEEXT,filestarimage.local_bin_index),
                            TblByt,"seriesId",TblOff(PFILESTARIMAGEEXT,filestarimage.seriesId),
                            TblByt,"exposureNumber",TblOff(PFILESTARIMAGEEXT,filestarimage.exposureNumber),
                            TblByt,"solutionNumber",TblOff(PFILESTARIMAGEEXT,filestarimage.solutionNumber),
                            TblByt,"spatial_bin",TblOff(PFILESTARIMAGEEXT,filestarimage.spatial_bin),
                            TblInt,"mosaicNumber",TblOff(PFILESTARIMAGEEXT,mosaicNumber),
                            TblInt,"quality",TblOff(PFILESTARIMAGEEXT,quality),
                            TblInt,"plateVersionId",TblOff(PFILESTARIMAGEEXT,plateVersionId),
                            TblInt,"catalogNumber",TblOff(PFILESTARIMAGEEXT,filestarimage.catalogNumber),
                            0,"end",0);
    if (db_table == NULL) {
      printf("ERROR: Failed to read table for %s\n",db_name);
      return -1;
    }
    numMagnitudes = db_nrecs;

    if (db_handle != NULL) {
      Close(db_handle);
    }
    if (db_header != NULL) {
      table_hdrfree(db_header);
    }
  }

  /* No decide how many series we have and initialize unused fields */

  for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
    pFileStarImageExt = &db_table[magnitudeIndex];
    pPhotPlates = &pFileStarImageExt->photplates;
    memset(pPhotPlates,0,sizeof(PHOTPLATES));
    pPhotPlates->mosaicNumber = pFileStarImageExt->mosaicNumber;
    /* Fix of Mar 22, 2016 - comply with original assumption that QUALITY_UNDETECTED stands alone */
    if ((pFileStarImageExt->quality & QUALITY_UNDETECTED) != 0) {
      pFileStarImageExt->quality = QUALITY_UNDETECTED;
    }
    pPhotPlates->quality = pFileStarImageExt->quality;
    pPhotPlates->versionId = pFileStarImageExt->plateVersionId;
    if ((pFileStarImageExt->filestarimage.seriesId > 0) &&
        (pFileStarImageExt->filestarimage.seriesId <= MAX_SERIES)) {
      pSeriesList = &seriesList[pFileStarImageExt->filestarimage.seriesId];
      pSeriesList->numPlates++;
    }
  }
  seriesCount = START_SERIES_BIT;
  for (seriesIndex = 0; seriesIndex <= MAX_SERIES; seriesIndex++) {
    pSeriesList = &seriesList[seriesIndex];
    if ((pSeriesList->numPlates > 0) &&
        (pSeriesList->bitMask == 0) &&
        (seriesCount <= MAX_SERIES_BIT)) {
      pSeriesList->bitPosition = seriesCount;
      pSeriesList->bitMask = 1 << seriesCount;
      allSeriesMask |= pSeriesList->bitMask;
      seriesCount++;

    } else {
      pSeriesList->numPlates = 0; /* Can not include this series! */
    }
  }

  /* Now create the plotting array */
  if (numMagnitudes > 0) {
    if (foldflag) {
      numPlotted = 2*numMagnitudes;
    } else {
      numPlotted = numMagnitudes;

    }

    col = (int *)calloc(numPlotted,sizeof(int));
    crossindex = (int *)calloc(numPlotted,sizeof(int));
    x = (double *)calloc(numPlotted,sizeof(double));
    y = (double *)calloc(numPlotted,sizeof(double));
    err = (double *)calloc(numPlotted,sizeof(double));
    limiting = (double *)calloc(numPlotted,sizeof(double));

    vector1 = (double *)calloc(numMagnitudes,sizeof(double));
    vector2 = (double *)calloc(numMagnitudes,sizeof(double));
    vector3 = (double *)calloc(numMagnitudes,sizeof(double));
    vector1s = (double *)calloc(numMagnitudes,sizeof(double));
    vector2s = (double *)calloc(numMagnitudes,sizeof(double));
    vector3s = (double *)calloc(numMagnitudes,sizeof(double));
    vector4s = (double *)calloc(numMagnitudes,sizeof(double));
    vector1f = (double *)calloc(numMagnitudes,sizeof(double));
    vector2f = (double *)calloc(numMagnitudes,sizeof(double));
    vector3f = (double *)calloc(numMagnitudes,sizeof(double));
    vector4f = (double *)calloc(numMagnitudes,sizeof(double));
    if ((vector1 == NULL) ||
        (vector2 == NULL) ||
        (vector3 == NULL) ||
        (vector1s == NULL) ||
        (vector2s == NULL) ||
        (vector3s == NULL) ||
        (vector4s == NULL) ||
        (vector1f == NULL) ||
        (vector2f == NULL) ||
        (vector3f == NULL) ||
        (vector4f == NULL)) {
      printf("ERROR: failed to allocate statistics vectors of size %d\n",numMagnitudes);
      exit(1);
    }

    if ((col == NULL) || (x == NULL) || (y == NULL)) {
      printf("ERROR allocating col %p, x %p, or y %p\n", col, x, y);
      exit(1);
    }
    CURAFLAGSMASK = 0;
    curqualitymask = QUALITY_UNDETECTED|QUALITY_LIMITING;
    curseriesmask = 0;
    if (listString == NULL) {
      listStringAlloc = LIST_STRING_ALLOC;
      listString = (char *)calloc(LIST_STRING_ALLOC,sizeof(char));
      if (listString == NULL) {
        printf("ERROR: Failed to allocate listString of size %d\n",LIST_STRING_ALLOC);
        exit(1);
      }
    }
    listString[0] = 0;
    listStringSize = 0;

    pSortTable = (PSORTTABLE)calloc(numMagnitudes,sizeof(SORTTABLE));
    if (pSortTable == NULL) {
      printf("ERROR: Failed to allocate pSortTable\n");
      exit(1);
    }
    /* Get the range in dates and magnitudes */
    for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
      pFileStarImageExt = &db_table[magnitudeIndex];

      dateyear = jd2ep(pFileStarImageExt->filestarimage.Date);
      tempmag = pFileStarImageExt->filestarimage.magcal_magdep;
      templimiting = pFileStarImageExt->filestarimage.limiting_mag_local;

      if (dateyear < actstartdate) {
        actstartdate = dateyear;
      }
      if (dateyear > actenddate) {
        actenddate = dateyear;
      }
      if (pFileStarImageExt->quality != QUALITY_UNDETECTED) {
        if ((tempmag < 90.0) && (tempmag > actdimmag)) {
          actdimmag = tempmag;
        }
        if (tempmag < actbrightmag) {
          actbrightmag = tempmag;
        }
      }
      if ((templimiting < 90.0) && (templimiting > actdimlimiting)) {
        actdimlimiting = templimiting;
      }
      if (templimiting < actbrightlimiting) {
        actbrightlimiting = templimiting;
      }
    }

    if (reqXValid == 0) {
      reqstartdate = actstartdate;
      reqenddate   = actenddate;
    }

    if ((reqYValid == 0) || (resetToOriginal != 0))
        {
          reqbrightmag = actbrightmag;
          reqdimmag = actdimmag;
          if ((qualitymask & (QUALITY_UNDETECTED|QUALITY_LIMITING)) != 0) {
            if (actbrightlimiting < reqbrightmag) {
              reqbrightmag = actbrightlimiting;
            }
            if (actdimlimiting > reqdimmag) {
              reqdimmag = actdimlimiting;
            }

          }

          if (resetToOriginal != 0)
              {
                reqstartdate = actstartdate;
                reqenddate = actenddate;

              }
        }

    /* Calculate the error_bar_factor */

    goodPoints = 0;
    for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
      pFileStarImageExt = &db_table[magnitudeIndex];
      if (pFileStarImageExt->quality == QUALITY_UNDETECTED) {
        continue;
      }

      DecodeAFLAGS( pFileStarImageExt->filestarimage.AFLAGS,AFLAGSMASK, pFileStarImageExt->filestarimage.BFLAGS,BFLAGSMASK,&reject_reason1,&reject_reason2,&flag);
      if ((flag == SYMBOL_GOOD) || (flag == SYMBOL_NOCOLOR)  || (flag == SYMBOL_NEIGHBORS) || (flag == SYMBOL_SATURATED)  || (flag == SYMBOL_NOMAGDEP)) {
        if ((pFileStarImageExt->filestarimage.magcal_magdep < 90.0) &&
            (pFileStarImageExt->filestarimage.magcal_local_rms >= 0.0) &&
            (pFileStarImageExt->filestarimage.magcal_local_rms < 90.0)) {
          /* This currently includes stars near the limiting magnitude */
          vector1[goodPoints] = pFileStarImageExt->filestarimage.magcal_magdep;
          vector2[goodPoints] = pFileStarImageExt->filestarimage.magcal_local_rms;
          goodPoints++;
        }
      }
    }
    if (CalcMedianAndRMS(goodPoints,2,vector1,&rawmed,&rawrms,0,3.0,0) == 0) {
      rawmed = 0.0;
      rawrms = 99.0;
    }
    /* Find the zero-based clipped rms of our error bars */
    if (CalcMedianAndRMS(goodPoints,2,vector2,&rawerrmed,&rawerrrms,0,3.0,1) == 0) {
      rawerrmed = 0.0;
      rawerrrms = 99.0;
    }
    if ((rawrms == 0.0) ||
        (rawrms == 99.0) ||
        (rawerrrms == 0.0) ||
        (rawerrrms == 99.0)) {
      error_bar_factor = 1.0;
    } else {
      error_bar_factor = rawrms/rawerrrms;
    }
    if (error_bar_factor > 1.0) {
      error_bar_factor = 1.0;
    }

    /* At this point decide if we are going to have a short file and prepare the header */

    db_indices = (int *)calloc(numMagnitudes,sizeof(int));
    if (db_indices == NULL) {
      printf("ERROR: failed to allocate db_indices of size %d\n",numMagnitudes);
    }
    plotshort_indices = (int *)calloc(numMagnitudes,sizeof(int));
    if (plotshort_indices == NULL) {
      printf("ERROR: failed to allocate plotshort_indices of size %d\n",numMagnitudes);
    }
    fullshort_indices = (int *)calloc(numMagnitudes,sizeof(int));
    if (fullshort_indices == NULL) {
      printf("ERROR: failed to allocate fullshort_indices of size %d\n",numMagnitudes);
    }

    for (magnitudeIndex = 0; magnitudeIndex < numMagnitudes; magnitudeIndex++) {
      pFileStarImageExt = &db_table[magnitudeIndex];

      CURAFLAGSMASK |= pFileStarImageExt->filestarimage.AFLAGS;
      curqualitymask |= pFileStarImageExt->quality;
      webquality = pFileStarImageExt->quality;
      pSeriesList = &seriesList[pFileStarImageExt->filestarimage.seriesId];
      curseriesmask |= pSeriesList->bitMask;


      /* Handle pseudo-quality bits */
      if ((pFileStarImageExt->filestarimage.BFLAGS & (1 << FILTER_BFLAG_PSFSATURATED)) != 0) {
        curqualitymask |= QUALITY_SATURATED;
        webquality |=  QUALITY_SATURATED;
      }
      if (((pFileStarImageExt->filestarimage.BFLAGS & (1 << FILTER_BFLAG_MAGDEP_MAGCOR)) == 0) &&
          (pFileStarImageExt->quality != QUALITY_UNDETECTED)) {
        curqualitymask |= QUALITY_NOMAGDEP;
        webquality |=  QUALITY_NOMAGDEP;
      }
      if (GetFittedPlateScale(pFileStarImageExt->filestarimage.seriesId,pFileStarImageExt->filestarimage.plateNumber) >= PATROL_PLATE_SCALE) {
        curqualitymask |= QUALITY_PATROL;
        webquality |= QUALITY_PATROL;
      } else {
        curqualitymask |= QUALITY_NONPATROL;
        webquality |= QUALITY_NONPATROL;
      }
      if (pFileStarImageExt->filestarimage.ELLIPTICITY > MAX_ELLIPTICITY) {
        webquality |= QUALITY_TRAILED;
        curqualitymask |= QUALITY_TRAILED;
      }

      if (pFileStarImageExt->filestarimage.versionId != pFileStarImageExt->plateVersionId) {
        /* Skip stale version IDs here */
        printf(".");
        continue;
      }

      /* Collect statistics for all valid points */
      if ((pFileStarImageExt->filestarimage.magcal_local_rms < 90) &&
          (pFileStarImageExt->quality != QUALITY_UNDETECTED)) {
        vector1f[fullshort_statndata] = pFileStarImageExt->filestarimage.ra;
        vector2f[fullshort_statndata] = pFileStarImageExt->filestarimage.dec;
        vector3f[fullshort_statndata] = pFileStarImageExt->filestarimage.magcal_magdep;
        vector4f[fullshort_statndata] = pFileStarImageExt->filestarimage.magcal_local_rms;
        fullshort_statndata++;
      }
      fullshort_indices[fullshort_ndata] = magnitudeIndex;
      fullshort_ndata++;
      webseriesbit = 0;
      if ((pFileStarImageExt->filestarimage.seriesId > 0) &&
          (pFileStarImageExt->filestarimage.seriesId <= MAX_SERIES)) {
        pSeriesList = &seriesList[pFileStarImageExt->filestarimage.seriesId];
        webseriesbit = pSeriesList->bitMask;
      }

      if (((pFileStarImageExt->filestarimage.AFLAGS & ~AFLAGSMASK) == 0) &&
          ((webquality & ~qualitymask) == 0) &&
          ((webseriesbit & ~seriesmask) == 0) &&
          (pFileStarImageExt->mosaicNumber != 99)) {
        db_indices[fullndata] = magnitudeIndex;
        fullndata++;
      }
    }

    for (dbindex = 0; dbindex < fullndata; dbindex++) {
      magnitudeIndex = db_indices[dbindex];
      pFileStarImageExt = &db_table[magnitudeIndex];

      CURAFLAGSMASK |= pFileStarImageExt->filestarimage.AFLAGS;
      curqualitymask |= pFileStarImageExt->quality;
      pSeriesList = &seriesList[pFileStarImageExt->filestarimage.seriesId];
      curseriesmask |= pSeriesList->bitMask;


      /* Handle pseudo-quality bits */
      if ((pFileStarImageExt->filestarimage.BFLAGS & (1 << FILTER_BFLAG_PSFSATURATED)) != 0) {
        curqualitymask |= QUALITY_SATURATED;
      }
      if (((pFileStarImageExt->filestarimage.BFLAGS & (1 << FILTER_BFLAG_MAGDEP_MAGCOR)) == 0) &&
          (pFileStarImageExt->quality != QUALITY_UNDETECTED)) {
        curqualitymask |= QUALITY_NOMAGDEP;
      }
      if (GetFittedPlateScale(pFileStarImageExt->filestarimage.seriesId,pFileStarImageExt->filestarimage.plateNumber) >= PATROL_PLATE_SCALE) {
        curqualitymask |= QUALITY_PATROL;
      } else {
        curqualitymask |= QUALITY_NONPATROL;
      }
      if (pFileStarImageExt->filestarimage.ELLIPTICITY > MAX_ELLIPTICITY) {
        curqualitymask |= QUALITY_TRAILED;
      }


      if (pFileStarImageExt->filestarimage.versionId != pFileStarImageExt->plateVersionId) {
        /* Skip stale version IDs here */
        continue;
      }

      webseriesbit = 0;
      if ((pFileStarImageExt->filestarimage.seriesId > 0) &&
          (pFileStarImageExt->filestarimage.seriesId <= MAX_SERIES)) {
        pSeriesList = &seriesList[pFileStarImageExt->filestarimage.seriesId];
        webseriesbit = pSeriesList->bitMask;
      }


      tempmag  = pFileStarImageExt->filestarimage.magcal_magdep;
      templimiting = pFileStarImageExt->filestarimage.limiting_mag_local;

      dateyear = jd2ep(pFileStarImageExt->filestarimage.Date);

      if (foldflag) {
        foldfactor = ((pFileStarImageExt->filestarimage.Date-foldcenter)/foldperiod);
        x[ndata] = (pFileStarImageExt->filestarimage.Date-foldcenter) - (foldperiod*foldfactor);
        x[ndata] = x[ndata]/foldperiod;
        while (x[ndata] < 0.0) {
          x[ndata] += 1.0;
        }
        while (x[ndata] > 1.0) {
          x[ndata] -= 1.0;
        }
      } else {
        x[ndata] = dateyear;
      }

      err[ndata] = pFileStarImageExt->filestarimage.magcal_local_rms * error_bar_factor;
      if ((err[ndata] > ERROR_BAR_CLIP) ||
          (err[ndata] < 0)) {
        err[ndata] = ERROR_BAR_CLIP;
      }

      if ((dateyear >= (reqstartdate-0.001)) &&
          (dateyear <= (reqenddate+0.001))) {
        if (((pFileStarImageExt->quality != QUALITY_UNDETECTED) &&
             ((tempmag >= (reqbrightmag-0.001)) &&
              (tempmag <= (reqdimmag+0.001)))) ||
            ((pFileStarImageExt->quality == QUALITY_UNDETECTED) &&
             ((templimiting >= (reqbrightmag-0.001)) &&
              (templimiting <= (reqdimmag+0.001))))) {
          y[ndata] = tempmag;
          limiting[ndata] = templimiting;
          if (pFileStarImageExt->quality != QUALITY_UNDETECTED) {
            vector1[statndata2] = pFileStarImageExt->filestarimage.ra;
            vector2[statndata2] = pFileStarImageExt->filestarimage.dec;
            vector3[statndata2] = tempmag;
            statndata2++;
          }
          if (pFileStarImageExt->quality == QUALITY_UNDETECTED) {
            finalbright = templimiting;
            finaldim = templimiting;
          } else {
            if ((qualitymask & QUALITY_LIMITING) == 0) {
              finalbright = tempmag;
              finaldim = tempmag;
            } else {
              finalbright = tempmag;
              finaldim = fmax(tempmag,templimiting);
            }
          }
          finalbrightmag = fmin(finalbrightmag,finalbright);
          finaldimmag    = fmax(finaldimmag,finaldim);
          crossindex[ndata] = magnitudeIndex;

          if (pFileStarImageExt->quality == QUALITY_UNDETECTED) {
            col[ndata] = COLOR_GRAY; /* Grey */
          } else {
            if ((pFileStarImageExt->filestarimage.AFLAGS & (1<<FILTER_AFLAG_TOO_BRIGHT)) != 0) {
              col[ndata] = COLOR_BLUE; /* Blue */
            } else if ((pFileStarImageExt->filestarimage.AFLAGS & (~(1 << GSC_VARIABLE_BIT))) != 0) {
              col[ndata] = COLOR_RED; /* Red */
            }
            /* Otherwise black CATALOG_GSC232 and CATALOG_GAIA and CATALOG_ATLAS */
          }

          sprintf(tmpListString,"%s %d %d %9.4f;",GetSeriesString(pFileStarImageExt->filestarimage.seriesId,1),pFileStarImageExt->filestarimage.plateNumber,pFileStarImageExt->mosaicNumber,jd2ep(pFileStarImageExt->filestarimage.Date));
          if (strlen(tmpListString) >= LIST_SUBSTR_SIZE) {
            printf("ERROR: string %s is larger than %d characters\n",tmpListString,LIST_SUBSTR_SIZE);
          } else {
            pSortEntry = &pSortTable[sortTableSize];
            pSortEntry->julianDate = pFileStarImageExt->filestarimage.Date;
            strcpy(pSortEntry->listString,tmpListString);
            sortTableSize++;
          }

          if ((strlen(tmpListString) + listStringSize) > (listStringAlloc-5)) {
            listStringAlloc += LIST_STRING_ALLOC;
            newListString = (char *)calloc(listStringAlloc,sizeof(char));
            if (newListString == NULL) {
              printf("ERROR: Failed to allocate listString of size %d\n",LIST_STRING_ALLOC);
            }
            strcpy(newListString,listString);
            free(listString);
            listString = newListString;
            newListString = NULL;
          }
          strcat(listString,tmpListString);
          listStringSize = strlen(listString);
          /* Collect statistics for those points plotted */
          if ((pFileStarImageExt->filestarimage.magcal_local_rms < 90) &&
              (pFileStarImageExt->quality != QUALITY_UNDETECTED)) {
            vector1s[plotshort_statndata] = pFileStarImageExt->filestarimage.ra;
            vector2s[plotshort_statndata] = pFileStarImageExt->filestarimage.dec;
            vector3s[plotshort_statndata] = pFileStarImageExt->filestarimage.magcal_magdep;
            vector4s[plotshort_statndata] = pFileStarImageExt->filestarimage.magcal_local_rms;
            plotshort_statndata++;
          }
          plotshort_indices[plotshort_ndata] = magnitudeIndex;
          plotshort_ndata++;

          ndata++;
        }
      }
    }
  }
  if (fullshort_ndata > 0) {
    /* write out the short table for all points now */

    fullshort_dbhandle = fopen(fullshort_db_name,"wt");
    if (fullshort_dbhandle == NULL) {
      printf("ERROR: failed (2) to open %s for writing\n",fullshort_db_name);
      exit(1);
    }
    fullshort_txthandle = fopen(fullshort_txt_name,"wt");
    if (fullshort_txthandle == NULL) {
      printf("ERROR: failed (2) to open %s for writing\n",fullshort_txt_name);
      exit(1);
    }
    fprintf(fullshort_dbhandle,"title %s\n",starbase_title);
    fprintf(fullshort_txthandle,"title %s\n",starbase_txt_title);
    if (fullshort_statndata > 0) {
      if (CalcMedianAndRMS(fullshort_statndata,2,vector1f,&fullshort_med,&fullshort_rms,0,3.0,0) != 0) {
        fprintf(fullshort_dbhandle,"Median of Right Ascension  %9.5f, rms %9.5f\n",fullshort_med,fullshort_rms);
        fprintf(fullshort_txthandle,"Median of Right Ascension  %9.5f, rms %9.5f\n",fullshort_med,fullshort_rms);
      }
      if (CalcMedianAndRMS(fullshort_statndata,2,vector2f,&fullshort_med,&fullshort_rms,0,3.0,0) != 0) {
        fprintf(fullshort_dbhandle,"Median of Declination      %9.5f, rms %9.5f\n",fullshort_med,fullshort_rms);
        fprintf(fullshort_txthandle,"Median of Declination      %9.5f, rms %9.5f\n",fullshort_med,fullshort_rms);
      }
      if (CalcMedianAndRMS(fullshort_statndata,2,vector3f,&fullshort_med,&fullshort_rms,0,3.0,0) != 0) {
        fprintf(fullshort_dbhandle,"Median of magcal_magdep    %9.2f, rms %9.2f\n",fullshort_med,fullshort_rms);
        fprintf(fullshort_txthandle,"Median of magcal_magdep    %9.2f, rms %9.2f\n",fullshort_med,fullshort_rms);
      }
      if (CalcMedianAndRMS(fullshort_statndata,2,vector4f,&fullshort_med,&fullshort_rms,0,3.0,0) != 0) {
        fprintf(fullshort_dbhandle,"Median of magcal_local_rms %9.2f, rms %9.2f\n",fullshort_med,fullshort_rms);
        fprintf(fullshort_txthandle,"Median of magcal_local_rms %9.2f, rms %9.2f\n",fullshort_med,fullshort_rms);
      }
    }

    fprintf(fullshort_dbhandle, "Date\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\n");
    fprintf(fullshort_txthandle,"Date\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\n");
    fprintf(fullshort_dbhandle, "----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t-------\t------------\t----------\t----------\t-----------\n");
    fprintf(fullshort_txthandle,"----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t-------\t------------\t----------\t----------\t-----------\n");

    for (fullshort_index = 0; fullshort_index < fullshort_ndata; fullshort_index++) {
      magnitudeIndex = fullshort_indices[fullshort_index];
      pFileStarImageExt = &db_table[magnitudeIndex];

      double fullshort_year;
      char fullshort_plate[MAX_FILENAME];

      if (pFileStarImageExt->filestarimage.solutionNumber == 0) {
        sprintf(fullshort_plate,"%s%05d",GetSeriesString(pFileStarImageExt->filestarimage.seriesId,0),pFileStarImageExt->filestarimage.plateNumber);
      } else {
        sprintf(fullshort_plate,"%s%05d_s%d",GetSeriesString(pFileStarImageExt->filestarimage.seriesId,0),pFileStarImageExt->filestarimage.plateNumber,pFileStarImageExt->filestarimage.solutionNumber);
      }
      fullshort_year = jd2ep(pFileStarImageExt->filestarimage.Date);

      FormatFlagsBits(pFileStarImageExt->filestarimage.AFLAGS,AFLAGSBuffer,AFLAGSbitsBuffer,sizeof(AFLAGSBuffer),0);
      FormatFlagsBits(pFileStarImageExt->quality,qualityBuffer,qualitybitsBuffer,sizeof(qualityBuffer),1);
      FormatFlagsBits(pFileStarImageExt->filestarimage.BFLAGS,BFLAGSBuffer,BFLAGSbitsBuffer,sizeof(BFLAGSBuffer),0);

      fprintf(fullshort_dbhandle,"%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\n",
              pFileStarImageExt->filestarimage.Date,
              fullshort_year,
              pFileStarImageExt->filestarimage.magcal_magdep,
              pFileStarImageExt->filestarimage.magcal_local_rms,
              pFileStarImageExt->filestarimage.limiting_mag_local,
              pFileStarImageExt->filestarimage.ra,
              pFileStarImageExt->filestarimage.dec,
              pFileStarImageExt->filestarimage.THETA_J2000,
              pFileStarImageExt->filestarimage.ELLIPTICITY,
              fullshort_plate,
              pFileStarImageExt->filestarimage.versionId,
              AFLAGSBuffer,
              BFLAGSBuffer,
              pFileStarImageExt->filestarimage.NUMBER,
              qualityBuffer,
              pFileStarImageExt->filestarimage.timeAccuracy,
              AFLAGSbitsBuffer,
              BFLAGSbitsBuffer,
              qualitybitsBuffer);
      fprintf(fullshort_txthandle,"%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\n",
              pFileStarImageExt->filestarimage.Date,
              fullshort_year,
              pFileStarImageExt->filestarimage.magcal_magdep,
              pFileStarImageExt->filestarimage.magcal_local_rms,
              pFileStarImageExt->filestarimage.limiting_mag_local,
              pFileStarImageExt->filestarimage.ra,
              pFileStarImageExt->filestarimage.dec,
              pFileStarImageExt->filestarimage.THETA_J2000,
              pFileStarImageExt->filestarimage.ELLIPTICITY,
              fullshort_plate,
              pFileStarImageExt->filestarimage.versionId,
              AFLAGSBuffer,
              BFLAGSBuffer,
              pFileStarImageExt->filestarimage.NUMBER,
              qualityBuffer,
              pFileStarImageExt->filestarimage.timeAccuracy,
              AFLAGSbitsBuffer,
              BFLAGSbitsBuffer,
              qualitybitsBuffer);
    }

    fclose(fullshort_dbhandle);
    fclose(fullshort_txthandle);
    sprintf(cmdStr,"%s/cp %s %s\n",lsbin,fullshort_db_name,fullshort_zipdb_name);
    result = system(cmdStr);
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,fullshort_db_name);
    result = system(cmdStr);
    sprintf(cmdStr,"%s/mv %s %s\n",lsbin,fullshort_zipdb_name,fullshort_db_name);
    result = system(cmdStr);


    sprintf(cmdStr,"%s/votable -i %s -o %s\n",binaries,fullshort_db_name,fullshort_vo_name);


    result = system(cmdStr);
    if (result != 0) {
      printf("ERROR: result %d creating the short votable\n",result);
    }
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,fullshort_vo_name);
    result = system(cmdStr);
  }

  if (plotshort_ndata > 0) {
    /* Write out the short table for plotted points now */

    plotshort_dbhandle = fopen(plotshort_db_name,"wt");
    if (plotshort_dbhandle == NULL) {
      printf("ERROR: failed (2) to open %s for writing\n",plotshort_db_name);
      exit(1);
    }
    plotshort_txthandle = fopen(plotshort_txt_name,"wt");
    if (plotshort_txthandle == NULL) {
      printf("ERROR: failed (2) to open %s for writing\n",plotshort_txt_name);
      exit(1);
    }
    fprintf(plotshort_dbhandle,"title %s\n",starbase_title);
    fprintf(plotshort_txthandle,"title %s\n",starbase_txt_title);
    if (plotshort_statndata > 0) {
      if (CalcMedianAndRMS(plotshort_statndata,2,vector1s,&plotshort_med,&plotshort_rms,0,3.0,0) != 0) {
        fprintf(plotshort_dbhandle,"Median of Right Ascension  %9.5f, rms %9.5f\n",plotshort_med,plotshort_rms);
        fprintf(plotshort_txthandle,"Median of Right Ascension  %9.5f, rms %9.5f\n",plotshort_med,plotshort_rms);
      }
      if (CalcMedianAndRMS(plotshort_statndata,2,vector2s,&plotshort_med,&plotshort_rms,0,3.0,0) != 0) {
        fprintf(plotshort_dbhandle,"Median of Declination      %9.5f, rms %9.5f\n",plotshort_med,plotshort_rms);
        fprintf(plotshort_txthandle,"Median of Declination      %9.5f, rms %9.5f\n",plotshort_med,plotshort_rms);
      }
      if (CalcMedianAndRMS(plotshort_statndata,2,vector3s,&plotshort_med,&plotshort_rms,0,3.0,0) != 0) {
        fprintf(plotshort_dbhandle,"Median of magcal_magdep    %9.2f, rms %9.2f\n",plotshort_med,plotshort_rms);
        fprintf(plotshort_txthandle,"Median of magcal_magdep    %9.2f, rms %9.2f\n",plotshort_med,plotshort_rms);
      }
      if (CalcMedianAndRMS(plotshort_statndata,2,vector4s,&plotshort_med,&plotshort_rms,0,3.0,0) != 0) {
        fprintf(plotshort_dbhandle,"Median of magcal_local_rms %9.2f, rms %9.2f\n",plotshort_med,plotshort_rms);
        fprintf(plotshort_txthandle,"Median of magcal_local_rms %9.2f, rms %9.2f\n",plotshort_med,plotshort_rms);
      }
    }

    fprintf(plotshort_dbhandle, "Date\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\n");
    fprintf(plotshort_txthandle,"Date\tyear\tmagcal_magdep\tmagcal_local_rms\tlimiting_mag_local\tra\tdec\tTHETA_J2000\tELLIPTICITY\tPlate\tversionId\tAFLAGS\tBFLAGS\tNUMBER\tquality\ttimeAccuracy\tAFLAGSBits\tBFLAGSBits\tqualitybits\n");
    fprintf(plotshort_dbhandle, "----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\n");
    fprintf(plotshort_txthandle,"----\t----\t-------------\t----------------\t------------------\t--\t---\t-----------\t-----------\t-----\t---------\t------\t------\t------\t-------\t------------\t----------\t----------\t-----------\n");

    for (plotshort_index = 0; plotshort_index < plotshort_ndata; plotshort_index++) {
      magnitudeIndex = plotshort_indices[plotshort_index];
      pFileStarImageExt = &db_table[magnitudeIndex];

      double plotshort_year;
      char plotshort_plate[MAX_FILENAME];
      if (pFileStarImageExt->filestarimage.solutionNumber == 0) {
        sprintf(plotshort_plate,"%s%05d",GetSeriesString(pFileStarImageExt->filestarimage.seriesId,0),pFileStarImageExt->filestarimage.plateNumber);
      } else {
        sprintf(plotshort_plate,"%s%05d_s%d",GetSeriesString(pFileStarImageExt->filestarimage.seriesId,0),pFileStarImageExt->filestarimage.plateNumber,pFileStarImageExt->filestarimage.solutionNumber);
      }
      plotshort_year = jd2ep(pFileStarImageExt->filestarimage.Date);

      FormatFlagsBits(pFileStarImageExt->filestarimage.AFLAGS,AFLAGSBuffer,AFLAGSbitsBuffer,sizeof(AFLAGSBuffer),0);
      FormatFlagsBits(pFileStarImageExt->quality,qualityBuffer,qualitybitsBuffer,sizeof(qualityBuffer),1);
      FormatFlagsBits(pFileStarImageExt->filestarimage.BFLAGS,BFLAGSBuffer,BFLAGSbitsBuffer,sizeof(BFLAGSBuffer),0);

      fprintf(plotshort_dbhandle,"%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\n",
              pFileStarImageExt->filestarimage.Date,
              plotshort_year,
              pFileStarImageExt->filestarimage.magcal_magdep,
              pFileStarImageExt->filestarimage.magcal_local_rms,
              pFileStarImageExt->filestarimage.limiting_mag_local,
              pFileStarImageExt->filestarimage.ra,
              pFileStarImageExt->filestarimage.dec,
              pFileStarImageExt->filestarimage.THETA_J2000,
              pFileStarImageExt->filestarimage.ELLIPTICITY,
              plotshort_plate,
              pFileStarImageExt->filestarimage.versionId,
              AFLAGSBuffer,
              BFLAGSBuffer,
              pFileStarImageExt->filestarimage.NUMBER,
              qualityBuffer,
              pFileStarImageExt->filestarimage.timeAccuracy,
              AFLAGSbitsBuffer,
              BFLAGSbitsBuffer,
              qualitybitsBuffer);

      fprintf(plotshort_txthandle,"%.6f\t%.6f\t%5.2f\t%5.2f\t%5.2f\t%9.5f\t%9.5f\t%5.1f\t%.3f\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%s\t%s\t%s\n",
              pFileStarImageExt->filestarimage.Date,
              plotshort_year,
              pFileStarImageExt->filestarimage.magcal_magdep,
              pFileStarImageExt->filestarimage.magcal_local_rms,
              pFileStarImageExt->filestarimage.limiting_mag_local,
              pFileStarImageExt->filestarimage.ra,
              pFileStarImageExt->filestarimage.dec,
              pFileStarImageExt->filestarimage.THETA_J2000,
              pFileStarImageExt->filestarimage.ELLIPTICITY,
              plotshort_plate,
              pFileStarImageExt->filestarimage.versionId,
              AFLAGSBuffer,
              BFLAGSBuffer,
              pFileStarImageExt->filestarimage.NUMBER,
              qualityBuffer,
              pFileStarImageExt->filestarimage.timeAccuracy,
              AFLAGSbitsBuffer,
              BFLAGSbitsBuffer,
              qualitybitsBuffer);
    }

    fclose(plotshort_dbhandle);
    fclose(plotshort_txthandle);
    sprintf(cmdStr,"%s/cp %s %s\n",lsbin,plotshort_db_name,plotshort_zipdb_name);
    result = system(cmdStr);
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,plotshort_db_name);
    result = system(cmdStr);
    sprintf(cmdStr,"%s/mv %s %s\n",lsbin,plotshort_zipdb_name,plotshort_db_name);
    result = system(cmdStr);


    sprintf(cmdStr,"%s/votable -i %s -o %s\n",binaries,plotshort_db_name,plotshort_vo_name);

    result = system(cmdStr);
    if (result != 0) {
      printf("ERROR: result %d creating the short votable\n",result);
    }
    sprintf(cmdStr,"%s/gzip %s\n",lsbin,plotshort_vo_name);
    result = system(cmdStr);
  }

  {
    /* Here begins Grzegorz Pojmanski's code */
    int blue;
    int black;
    int grey;
    int colors[COLOR_MAX];
    double *s = NULL;
    double min = 0;
    double max = 0;
    char title[MAX_BUFFER] = "";
    char tmpbuf[MAX_BUFFER];
    char xlabel[MAX_LABEL] = "year";
    char ylabel[MAX_LABEL] = "mag";
    char label1[MAX_BUFFER] = "";
    char label2[MAX_BUFFER] = "";
    char label3[MAX_LABEL] = "";
    int style = PLOT_EXACT|PLOT_AVE|PLOT_FLIP|PLOT_UPPER_L;
    double xmin = -1.e32;
    double xmax = 1.e32;
    double upper_limit= UPPER_LIMIT;
    gdImagePtr im_out;
    int curStars;
    STARENTRY curStarEntry;
    PSTARENTRY pCurStarEntry = &curStarEntry;


    strcpy(title,REF);
    if (foldflag) {
      strcpy(xlabel,"phase");
    } else {
      strcpy(xlabel,"year");
    }
    curStars = GetStarEntry2(pPhotConnection,pCurStarEntry,REF,0,1);
    if (curStars == 1) {
      if (pCurStarEntry->updateflag == UPDATEFLAG_NEED_UPDATE) {
        /* Here the ra and dec are unreliable */
        rstr[0] = 0;
        dstr[0] = 0;
      } else {
        ra2str (rstr, 16, pCurStarEntry->ra, ndec);
        dec2str (dstr, 16, pCurStarEntry->dec, ndec-1);
      }
      sprintf(tmpbuf," (%s %s Stdmag: %.2f color %.2f)",rstr,dstr,pCurStarEntry->Stdmag,pCurStarEntry->color);
      if ((strlen(title) + strlen(tmpbuf) + 2) < MAX_BUFFER) {
        strcat(title,tmpbuf);
      }
    }

    sprintf(tmpbuf," is %.0f arcsec from ",distance);
    if ((strlen(title) + strlen(tmpbuf) + 2) < MAX_BUFFER) {
      strcat(title,tmpbuf);
      if (strcmp(REF,center) == 0) {
        strcpy(tmpbuf,"search center.");
        if ((strlen(title) + strlen(tmpbuf) + 2) < MAX_BUFFER) {
          strcat(title,tmpbuf);
        }
      } else {

        if ((strlen(title) + strlen(center) + 2) < MAX_BUFFER) {
          strcat(title,center);
        }
        strcat(title,".");
      }
    }
    if (m44release) {
      if ((strlen(title) + strlen(" Preliminary lightcurve.") + 2) < MAX_BUFFER) {
        strcat(title," Preliminary lightcurve.");
      }
    }

    actstartdate -= 0.5;
    actenddate += 0.5;
    if (fullndata != ndata) {
      sprintf(tmpbuf," Years: %.0f to %.0f Mags: %.2f to %.2f",actstartdate,actenddate,actdimmag,actbrightmag);
      if ((strlen(title) + strlen(tmpbuf) + 2) < MAX_BUFFER) {
        strcat(title,tmpbuf);
      }
    }


    if (ndata <= 0) {
      sprintf(label1,"NO POINTS SELECTED for AFLAGS=%d",AFLAGSMASK);
      if (AFLAGSMASK != 0) {
        strcpy(label2,"USE BUTTONS BELOW TO ALLOW PLOTTING OF ADDITIONAL POINTS");
      }
      if (fullndata != ndata) {
        sprintf(tmpbuf," OR RESET TO ORIGINAL YEARS: %.0f to %.0f and ORIGINAL MAGS %.2f to %.2f",actstartdate,actenddate,actdimmag,actbrightmag);
        if ((strlen(label2) + strlen(tmpbuf) + 2) < MAX_BUFFER) {
          strcat(label2,tmpbuf);
        }
      }


    } else {
      if (CalcMedianAndRMS(statndata2,0,vector1,&ramed,&rarms,0,3.0,0) == 0) {
        ramed = 0.0;
        rarms = 99.0;
      }
      if (CalcMedianAndRMS(statndata2,0,vector2,&decmed,&decrms,0,3.0,0) == 0) {
        decmed = 0.0;
        decrms = 99.0;
      }
      if (CalcMedianAndRMS(statndata2,0,vector3,&rawmed,&rawrms,0,3.0,0) == 0) {
        rawmed = 0.0;
        rawrms = 99.0;
      }
      if ((rarms == 99.0) || (decrms == 99.0)) {
        drad = 99.0;
      } else {
        drad = 3600*sqrt(sqr(rarms)+sqr(decrms));
      }

      ra2str (rstr, 16, ramed, ndec);
      dec2str (dstr, 16, decmed, ndec-1);

      sprintf(label1,"%s points: %d AFLAGS: %d ra: %s dec: %s rms: %.1f arcsec",catalogText[catalogNumber],ndata,AFLAGSMASK,rstr,dstr,drad);
      if (strlen(label1) > MAX_BUFFER) {
        printf("ERROR: label 1 length is %zu\n",strlen(label1));
        exit(1);
      }
      sprintf(label2,"mag: %.2f mag rms: %.2f error_bar_factor %.2f. Errors >= %.1f mag are indeterminate.",rawmed,rawrms,error_bar_factor,ERROR_BAR_CLIP);
      if (strlen(label1) > MAX_BUFFER) {
        printf("ERROR: label 1 length is %zu\n",strlen(label1));
        exit(1);
      }
      if (foldflag) {
        if (foldcenter > 0.0) {
          if (foldcenter < foldperiod) {
            if (snprintf(label3, sizeof(label3), "Folded plot with a period of %s days and offset with %f phase",foldperiodString,foldcenter/foldperiod) < 0) {
              printf("ERROR: buffer overflow for label3\n");
              exit(1);
            }
          } else {
            dateyear = jd2ep(foldcenter);
            if (snprintf(label3, sizeof(label3), "Folded plot with a period of %s days and centered on year %f or Julian Day %f",foldperiodString,dateyear,foldcenter) < 0) {
              printf("ERROR: buffer overflow for label3\n");
              exit(1);
            }
          }
        } else {
          if (snprintf(label3, sizeof(label3), "Folded plot with a period of %s days.",foldperiodString) < 0) {
            printf("ERROR: buffer overflow for label3\n");
            exit(1);
          }
        }
      } else {
        label3[0] = 0;
      }
    }
    im_out = gdImageCreate(plotWidth,plotHeight);
    gdImageColorAllocate(im_out, 255, 255, 255);
    blue  = gdImageColorAllocate(im_out,   0,   0, 255);
    black = gdImageColorAllocate(im_out,   0,   0,   0);
    grey = gdImageColorAllocate(im_out,  150,  150,  150);
    colors[0] = gdImageColorAllocate(im_out,   0,   0,   0);
    colors[1] = gdImageColorAllocate(im_out,   0,   0, 255); /* Blue */
    colors[2] = gdImageColorAllocate(im_out, 255,   0,   0);
    colors[3] = gdImageColorAllocate(im_out,   0, 255,   0); /* Green */
    colors[4] = gdImageColorAllocate(im_out, 150, 150, 150); /* Grey for too bright */
    colors[5] = gdImageColorAllocate(im_out, 150, 150, 150); /* Grey for limiting magnitude only */
    colors[6] = gdImageColorAllocate(im_out, 255, 255,   0); /* Yellow */
    colors[7] = gdImageColorAllocate(im_out, 255,   0, 255); /* Purple or magenta*/
    colors[8] = gdImageColorAllocate(im_out,   0, 255, 255); /* Cyan */
    if (reqstartdate < actstartdate) {
      reqstartdate = actstartdate;
    }
    if (reqenddate > actenddate) {
      reqenddate = actenddate;
    }
    if ((qualitymask & (QUALITY_UNDETECTED|QUALITY_LIMITING)) == 0) {
      /* Bugfix of Jul  8, 2013 - reverse the sign of the comparison */
      if (reqdimmag > actdimmag) {
        reqdimmag = actdimmag;
      }
      if (reqbrightmag < actbrightmag) {
        reqbrightmag = actbrightmag;
      }
    } else {
      if (reqdimmag > fmax(actdimmag,actdimlimiting)) {
        reqdimmag = fmax(actdimmag,actdimlimiting);
      }
      if (reqbrightmag < fmin(actbrightmag,actbrightlimiting)) {
        reqbrightmag = fmin(actbrightmag,actbrightlimiting);
      }
    }
    reqbrightmag = fmax(reqbrightmag,finalbrightmag);
    reqdimmag = fmin(reqdimmag,finaldimmag);

    if (foldflag) {
      for (index = 0; index < ndata; index++) {
        x[index+ndata] = x[index] - 1.0;
        y[index+ndata] = y[index];
        crossindex[index+ndata] = crossindex[index];
        col[index+ndata] = col[index];
        err[index+ndata] = err[index];
        limiting[index+ndata] = limiting[index];

      }
      ndata = 2*ndata;
    }
    if (RefType == REF_TYPE_NONE) {
      sprintf(title,"Limiting Magnitudes near %s.",center);
      sprintf(label1,"%s points: %d AFLAGS: %d",catalogText[catalogNumber],ndata,AFLAGSMASK);
      label2[0] = 0;
    }

    PlotCurve(im_out,
              0,   /* ix */
              0,   /* iy */
              plotWidth, /* sx */
              plotHeight, /* sy */
              ndata,
              x,
              y,
              s,
              err,
              limiting,
              col,
              crossindex,
              db_table,
              min, /* yymin */
              max, /* yymax */
              blue, /* col1 */
              black, /* col2 */
              colors,
              title,
              xlabel,
              ylabel,
              label1,
              label2,
              label3,
              style,
              store,
              upper_limit,
              xmin, /* xxmin */
              xmax, /* xxmax */
              grey,
              verbose,
              AFLAGSMASK,
              qualitymask,
              seriesmask,
              reqstartdate,
              reqenddate,
              reqdimmag,
              reqbrightmag,
              reqXDown,
              reqYDown,
              reqXUp,
              reqYUp,
              (qualitymask & (QUALITY_UNDETECTED|QUALITY_LIMITING)),
              averageRa,
              averageDec,
              foldflag,
              enableRematch,
              foldperiod,
              foldperiodString,
              foldcenter);
    gdImageInterlace(im_out, 1);
    gdImageGif(im_out, plotHandle);
    gdImageDestroy(im_out);

  }

  fclose(plotHandle);

  /* Now print the plotting checkboxes */

  if (m44release == 0) {
    int tempIndex;
    listString[0] = 0;
    qsort(pSortTable,sortTableSize,sizeof(SORTTABLE),dcmp);
    for (tempIndex = 0; tempIndex < sortTableSize; tempIndex++) {
      pSortEntry = &pSortTable[tempIndex];
      strcat(listString,pSortEntry->listString);
      if (strlen(listString) > listStringSize) {
        printf("ERROR: list string size %zu is larger than allocation %d\n",strlen(listString),listStringSize);
      }
    }

    printf("<form action=\"extracttarball.php\" id=\"extracttarball\" method=POST target=\"tarball\">");
    /* New POST method to handle large tarballs */
    printf("<input type=\"hidden\" name=\"filetype\" value = \"FITS\" />");
    printf("<input type=\"hidden\" name=\"archivetype\" value = \"TAR\"/>");
    printf("<input type=\"hidden\" name=\"ra\" value = \"%f\"/>",ramed);
    printf("<input type=\"hidden\" name=\"dec\" value = \"%f\"/>",decmed);
    printf("<input type=\"hidden\" name=\"system\" value = \"J2000\"/>");
    printf("<input type=\"hidden\" name=\"imagesize\" value = \"%d\"/>",imagesize);
    printf("<input type=\"hidden\" name=\"sizeunits\" value = \"%s\"/>",sizeunits);
    printf("<input type=\"hidden\" name=\"selection\" value = \"%s\"/>",listString);
    printf("<INPUT TYPE=submit NAME=submit VALUE=\"Get Tarball of Extracted FITS images\" >");

    printf("</form>");
  }

  urlencode(REF,REFurlencode,4*MAX_REF);
  centerurlencode = (char *)calloc(4*strlen(center),sizeof(char));
  if (centerurlencode == NULL) {
    printf("ERROR: allocation failed on centerurlencode of size %zu\n",strlen(center));
    exit(1);
  }
  urlencode(center,centerurlencode,(4*strlen(center)));

  charPtr = strstr(tmpdir,"tmp/");
  if (charPtr == NULL) {
    printf("ERROR: invalid format for tmpdir %s\n",tmpdir);
    exit(1);
  }
  charPtr += 4;
  strcpy(origtmpdir,charPtr);

  printf("<button onClick= redisplay('%s','%s','%.0f','%s','%s','%d','%s','%d','%d','0','0','0','0','0','%s','%f','%f'); return true;>Reset plot to original size</button><br />\n",
         REFurlencode,centerurlencode,distance,source,rematchString,gsc_bin_index,origtmpdir,windowWidth,windowHeight,reqWindow2,0.0,0.0);

  printf(
    "Black points have no quality flags set, but note that these flags have "
    "varying probabilities of false positive and false negatives.  Red points "
    "have one or more issues. Error bars greater than +/- %.1f mag are "
    "truncated to this magnitude and shown as dashed.  Stars too bright are "
    "shown as blue up arrows and limiting magnitudes are shown as horizontal "
    "bars.  Undetected images are shown as down arrows at the average limiting "
    "magnitude of a square degree region. Buttons below control the display of "
    "red points: a button with bold type means that the point is plotted. "
    "Parenthesis in some buttons indicate "
    "<a target=\"_blank\" href=\"https://dasch.cfa.harvard.edu/data-guide/#AFLAGS_ext\">AFLAGS</a> "
    "bit numbers. <br /><b>A two button procedure is necessary because the web "
    "plotter must also support detailed examination of events in a limited "
    "date range of the plot: A 'Show' button click may need to be followed by "
    "the above 'Reset plot to original size' to display the entire lightcurve.</b><br />",
    ERROR_BAR_CLIP
  );

  AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,2,0,0,0,reqWindow,foldperiodString,foldcenter); /* Good points */
  AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,1,0,0,0,reqWindow,foldperiodString,foldcenter); /* All points */
  AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,4,0,0,0,reqWindow,foldperiodString,foldcenter); /* Default points */

  /* Change of November 6, 2013 - move more frequently used buttons to the top */
  printf("<br >");
  AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,3,0,QUALITY_LIMITING,0,reqWindow,foldperiodString,foldcenter);
  AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,3,0,QUALITY_UNDETECTED,0,reqWindow,foldperiodString,foldcenter);
  if ((CURAFLAGSMASK & (1 << FILTER_AFLAG_DRAD)) != 0) {
    AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,3,1 << FILTER_AFLAG_DRAD,0,0,reqWindow,foldperiodString,foldcenter);
  }
  if ((CURAFLAGSMASK & (1 << FILTER_AFLAG_DEFECT)) != 0) {
    AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,3,1 << FILTER_AFLAG_DEFECT,0,0,reqWindow,foldperiodString,foldcenter);
  }
  printf("<br >");

  if ((CURAFLAGSMASK & (GSC_CLASS_NONSTAR << GSC_CLASS_BIT)) == (GSC_CLASS_NONSTAR << GSC_CLASS_BIT)) {
    AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,3,GSC_CLASS_NONSTAR << GSC_CLASS_BIT,0,0,reqWindow,foldperiodString,foldcenter);
  }
  for (bitindex = 0; bitindex < 32; bitindex++) {
    if ((CURAFLAGSMASK & (1 << bitindex)) != 0) {
      if ((bitindex != FILTER_AFLAG_QUALITY) &&
          (bitindex != GSC_VARIABLE_BIT) &&
          (bitindex != FILTER_AFLAG_DRAD) &&
          (bitindex != FILTER_AFLAG_DEFECT)) {
        AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,3,1 << bitindex,0,0,reqWindow,foldperiodString,foldcenter);
      }
    }
  }
  for (bitindex = 0; bitindex < 32; bitindex++) {
    if ((curqualitymask & (1 << bitindex)) != 0) {
      if (((1 << bitindex) != QUALITY_LIMITING) &&
          ((1 << bitindex) != QUALITY_UNDETECTED)) {
        AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,3,0,1 << bitindex,0,reqWindow,foldperiodString,foldcenter);
      }
    }
  }
  for (bitindex = 0; bitindex < 32; bitindex++) {
    if ((curseriesmask & (1 << bitindex)) != 0) {
      AddCheckbox(REFurlencode,centerurlencode,distance,source,rematchString,AFLAGSMASK,qualitymask,seriesmask,gsc_bin_index,tmpdir,windowWidth,windowHeight,3,0,0,1 << bitindex,reqWindow,foldperiodString,foldcenter);
    }
  }

  FreeFileCommon(pFileCommon0,0);

  if (pMagnitudeTable != NULL) {
    free(pMagnitudeTable);
  }
  if (pFileStarImageTable != NULL) {
    free(pFileStarImageTable);
  }
  if (pMagnitudeTable1 != NULL) {
    free(pMagnitudeTable1);
  }
  if (pMagnitudeTable2 != NULL) {
    free(pMagnitudeTable2);
  }


  if (pNoneMagnitudeTable != NULL) {
    free(pNoneMagnitudeTable);
  }

  mysql_close(pPhotConnection);
  mysql_close(pConnection);

  free(centerurlencode);
  centerurlencode = NULL;

  if (col != NULL) {
    free(col);
  }
  if (crossindex != NULL) {
    free(crossindex);
  }
  if (x != NULL) {
    free(x);
  }
  if (y != NULL) {
    free(y);
  }
  if (err != NULL) {
    free(err);
  }
  if (limiting != NULL) {
    free(limiting);
  }
  if (vector1 != NULL) {
    free(vector1);
  }
  if (vector2 != NULL) {
    free(vector2);
  }
  if (vector3 != NULL) {
    free(vector3);
  }
  if (vector1s != NULL) {
    free(vector1s);
  }
  if (vector2s != NULL) {
    free(vector2s);
  }
  if (vector3s != NULL) {
    free(vector3s);
  }
  if (vector4s != NULL) {
    free(vector4s);
  }
  if (db_table != NULL) {
    if (dbCallocFlag) {
      free(db_table);
      db_table = NULL;
    } else {
      Free(db_table);
      db_table = NULL;
    }
  }
  if (db_indices != NULL) {
    free(db_indices);
  }
  if (fullshort_indices != NULL) {
    free(fullshort_indices);
  }
  if (plotshort_indices != NULL) {
    free(plotshort_indices);
  }

  if (listString != NULL) {
    free(listString);
  }
  if (pSortTable != NULL) {
    free(pSortTable);
  }

  return 0;
}

/* May 26, 2010 Edward J. Los - Adapted from asas_plot_gif.c provided by
 *                              Grzegorz Pojmanski <gp@astrouw.edu.pl>  ASAS source author on Apr 26, 2010
 * Jun  9, 2010 Edward J. Los - Correct AddCheckBox for the "Good points" case.
 * Jul  2, 2010 Edward J. Los - Save AFLAGS, quality, reqstartdate, and reqenddate for the "movie" feature.
 * Sep  8, 2010 Edward J. Los - Adapt for Solaris environment
 * Dec 21, 2010 Edward J. Los - Make sure that the AFLAGS quality bit is set or cleared in
 *                              conformance with the qualitymask bitfield
 * Jul 18, 2011 Edward J. Los - Change THUMBNAIL_SIZEUNITS to THUMBNAIL_PIXELS
 * Nov 11, 2011 Edward J. Los - Support magnitude-dependent calibration
 * Nov 29, 2011 Edward J. Los - Do not plot points from non-Harvard plates
 * Dec  9, 2011 Edward J. Los - Filter on magnitude-dependent correction available
 * Dec 13, 2011 Edward J. Los - Convert to drag selection of plot height and width
 * Mar 23, 2012 Edward J. Los - Support plot selection by series
 * Mar 27, 2012 Edward J. Los - Sort the tarball string by Julian Date
 * Jul 18, 2012 Edward J. Los - Make sure that all calls to GetRef do not return spaces in object names
 * Oct 30, 2012 Edward J. Los - Add error bars by plotting magcal_local_rms.  Adopt error_bar_factor from
 *                              find_lightcurves, with the exception that BFLAGSMASK = 0
 * Nov 30, 2012 Edward J. Los - Correct a bug in recovering the object found by web_query.
 * Dec  4, 2012 Edward J. Los - Clip error bars to 2.0 magnitudes.
 *                              Add a second label line
 *                              Sort the results by Julian Date.
 * Dec 10, 2012 Edward J. Los - Add a title to the starbase record.
 * Dec 19, 2012 Edward J. Los - Add ellipticity, position angle, ephemeris date to the short table
 *                              Restrict the short table to only those points that pass the criteria.
 *                              Change FILTER_AMASK_PLOT to FILTER_AMASK_PLOT2
 *                              Generate the sort table for every plot, using only the points which pass the plotting criteria
 * Dec 31, 2012 Edward J. Los - Add limiting magnitude display support
 *                              Correct magcal_magdep_rms when re-reading the lightcurve table.
 * Mar 11, 2013 Edward J. Los - Add RaPM, DecPM, ra_2 and dec_2 to the magnitude file (Version 5);
 * Mar 30, 2013 Edward J. Los - Support new data release mechanism
 * Apr  5, 2013 Edward J. Los - Add versionId to the short plot
 * Apr 19, 2013 Edward J. Los - Turn off versionId mismatch errors, usually caused by transient conditions during database updates.
 * May  1, 2013 Edward J. Los - Use blue, not bold for undetected and "too bright" up arrows
 * Jun 29, 2013 Edward J. Los - Correct a bug which prevents using the tail end of the limiting magnitude table
 * Jul  8, 2013 Edward J. Los - Correct a sign error when checking the rescaling limits just prior to calling PlotCurve
 * Jul  9, 2013 Edward J. Los - Always look for unmatched objects where CHECKBLEND_BUGFIX is set
 * Jul 23, 2013 Edward J. Los - Include only the plotted points in the "short" table
 * Jul 26, 2013 Edward J. Los - Do not use the GSC_VARIABLE_BIT as an issue flag.  Rename "Good Points" to "Unflagged Points"
 * Oct  4, 2013 Edward J. Los - Reduce ERROR_BAR_CLIP to 0.4
 * Oct  7, 2013 Edward J. Los - If the rms is negative, force it to ERROR_BAR_CLIP
 *                            - If a DASCH object and we have an exact match to the REFNumber, use the object instead of the closest point.
 * Oct 26, 2013 Edward J. Los - Always expand the magnitude scale to the limits of the plotted regions (finalbrightmag and finaldimmag);
 * Nov  6, 2013 Edward J. Los - Add a comment about saving the .gif image of the lightcurve.
 *                              Move the comment about dashed error bars to the descriptive text.
 *                              Display Limiting Magnitudes, Undetected, High astrometric error, and plate defects on a separate line for prominance.
 * Nov 15, 2013 Edward J. Los - Fix bug in previous change to not show the High astrometric error and plate defect button for transient candidates.
 * Mar  6, 2014 Edward J. Los - Allow specification of fold periods with the full complement of double precision significant digits
 * Sep 15, 2014 Edward J. Los - Reject negative rms values when calculating error_bar_factor
 * Jan 20, 2015 Edward J. Los - V6 data format: add A2FLAGS, B2FLAGS, timeAccuracy, and maskIndex
 * Feb 16, 2015 Edward J. Los - Add timeAccuracy to the short form lightcurve table
 * Feb 26, 2015 Edward J. Los - Add a table of all points using the short column format (fullshort_ vs plotshort_)
 * May 16, 2015 Edward J. Los - List the bits in AFLAGS for the short output formats.
 * Sep 15, 2015 Edward J. Los - Add daschunistd.h for table.h conflicts
 * Feb 16, 2016 Edward J. Los - Correct exceptions when tmpdir is not defined
 *                              Fold only the previously selected points (FOLD_SELECTED)
 * Feb 22, 2016 Edward J. Los - In formatAflags, split AFLAGS into two columns: AFLAGS and AFLAGSBits
 * Feb 24, 2016 Edward J. Los - Change "formatAflags" to "FormatFlagsBits" and add BFLAGS and quality bitmaps
 * Mar  8, 2016 Edward J. Los - Use two separate fields for FormatFlagsBits
 * Mar 22, 2016 Edward J. Los - Correct failure to display undetected plates as a result of change from Feb 24.
 * Mar 24, 2017 Edward J. Los - Change GetStarEntry to GetStarEntry2 to move from the stars table to the starcatalog table.
 * Jan 19, 2017 Edward J. Los - Display the catalog in the lightcurve plot.
 * Jan 23, 2018 Edward J. Los   Support the merged experimental table: add catalogNumber to WriteStarbaseRecord.
 * Mar 19, 2018 Edward J. Los   support direct read of the "experimental" binary database.
 * Apr 13, 2018 Edward J. Los   Use yellow instead of greenfor the experimental combined display for AFLAGS != 0
 *                              Correct comments on the color definitions
 * May  8, 2018 Edward J. Los   Do not set resetToOriginal to 1 when displaying limiting magnitudes only. Allows zooming
 *                              in on regions (See memo of Thu, 3 May 2018 17:30:00)
 * Aug 13, 2018 Edward J. Los   Define enableRematch to optimize location of transients (Redefine -O) to support this function)
 * Dec 10, 2018 Edward J. Los   Add NUMBER to the short plot for positive identification of the image
 * Jan  4, 2019 Edward J. Los   Add verbose to LocateNoneImages and LoadNoneImages to display the files searched.
 * Apr  9, 2020 Edward J. Los   Add a stronger warning about the use of the "Reset plot to original size" button.
 */
