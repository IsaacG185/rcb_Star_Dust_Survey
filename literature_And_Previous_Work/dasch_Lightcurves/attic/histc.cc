// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
*HISTC Histogram count.
*   N = HISTC(X,EDGES), for vector X, counts the number of values in X
*   that fall between the elements in the EDGES vector (which must contain
*   monotonically non-decreasing values).  N is a LENGTH(EDGES) vector
*   containing these counts.  
*
*   N(k) will count the value X(i) if EDGES(k) <= X(i) < EDGES(k+1).  The
*   last bin will count any values of X that match EDGES(end).  Values
*   outside the values in EDGES are not counted.  Use -inf and inf in
*   EDGES to include all non-NaN values.
*
*
*   [N,BIN] = HISTC(X,...) also returns an index matrix BIN.  If X is a
*   vector, N(K) = SUM(BIN==K).
*
*   Aug 27, 2010 Edward J. Los - Update for Octave 3.2.3
*
*/

#define __USE_XOPEN2K8 1
#include <octave/oct.h>
#include <stdarg.h>
#include <stdio.h>


DEFUN_DLD (histc,args,, 
           "histc.\n\
\n\
Reference:\n\
\n\
  histogram algorithm")
{
  
  octave_value  xx = args(0);
  octave_value  yy = args(1);
  octave_value_list retval;

  NDArray xxx = xx.array_value();
  NDArray yyy = yy.array_value();
  
  double *x;
  double *y;

  int nbins = yyy.numel();
  int xsize = xxx.numel();

  int index;
#if 0
  int maxindex;
#endif
  int *histCount;
  int *histBin;
#if 0
  int nargin = args.length();
  printf("nargin is %d\n",nargin);
  printf("x length %d, y length %d,output length %d\n",xxx.length(),yyy.length(),nbins);
#endif
  x = xxx.fortran_vec();
  y = yyy.fortran_vec();
#if 0
  maxindex = xsize;
  if (nbins > maxindex) {
    maxindex = nbins;
  }
  for (index = 0; index < maxindex; index++) {
    if ((index < xsize) && (index < nbins)) {
      printf("index: %5d, y %f x %f\n",index,y[index],x[index]);
    } else if (index < nbins) {
      printf("index: %5d, y %f\n",index,y[index]);
    } else {
      printf("index: %5d,           x %f\n",index,x[index]);

    }

  }
#endif

  histCount = (int *)calloc(nbins,sizeof(int));
  histBin   = (int *)calloc(xsize,sizeof(int));

  for (index = 0; index < xsize; index++) {
    int k = -2;
    /* Ignore values too small */
    if (x[index] >= y[0]) {
      if (x[index] == y[nbins-1]) {
        k = nbins-1;
      }
      else if (x[index] < y[nbins-1]) {
        /* Use a binary search */
        int k0 = 0;
        int k1 = nbins-1;
        k = (k0+k1)/2;
        while (k0 < (k1-1)) {
          if (x[index] >= y[k]) {
            k0 = k;
          } else {
            k1 = k;
          }
          k = (k0+k1)/2;
        }
        
      }
    }
    if (k >= 0) {
      histCount[k] += 1;
      histBin[index] = k+1;
    } else {
      histBin[index] = -1;
    }
  }
 
#if 0
  for (index = 0; index < maxindex; index++) {
    if ((index < xsize) && (index < nbins)) {
      printf("index: %5d, histCount %d histBin %d\n",index,histCount[index],histBin[index]);
    } else if (index < nbins) {
      printf("index: %5d, histCount %d\n",index,histCount[index]);
    } else {
      printf("index: %5d,           histBin %d\n",index,histBin[index]);

    }

  }
#endif
  dim_vector dimsy = yyy.dims();
  
  NDArray zzz(dimsy);
  
  for (index = 0; index < nbins; index++) {
    zzz(index) = histCount[index];
  }

  dim_vector dimsx = xxx.dims();
  NDArray qqq(dimsx);

  for (index = 0; index < xsize; index++) {
    qqq(index) = histBin[index];
  }


  retval(0) = zzz;
  retval(1) = qqq;
#if 0
  printf("%s Returning\n",__FUNCTION__);
#endif
  free(histCount);
  free(histBin);
  return (retval);

}

