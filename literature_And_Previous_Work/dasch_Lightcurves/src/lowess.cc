// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* lowess.cc - Test the matlab lowess smoothing algorithm
 *
 *  mkoctfile lowess.cc
 *
 *
 *    Oct  9, 2007 Edward J. Los - correct calcleastmeansquare1 degenerate case where all
 *                                 of y is zero.
 *                                 Return 0 if both y and x are zero (octave does the same
 *                                 for matrix inversion).
 *    Aug  2, 2010 Edward J. Los - Add loess support
 *    Aug 27, 2010 Edward J. Los - Update for octave 3.2.3
 *
 * Adapted from the MatLab smooth.m script.
 *     Copyright 2001-2002 The MathWorks, Inc.
 *     $Revision: 1.17 $  $Date: 2002/05/29 18:43:31 $
 *
 *  Sorting borrowed from the following:
 *
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1999-2006   Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 *
 */

#define __USE_XOPEN2K8 1

#include <octave/oct.h>
#include <octave/version.h>

#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define EPS 2.2204e-16

static double
fcube(double x)
{
    return x * x * x;
}


static double
fsquare(double x)
{
    return x * x;
}


static int
rcmp(double x, double y)
{
  if (x < y)
    return -1;
  if (x > y)
    return 1;

  return 0;
}


static void
rPsort2(double *x, int lo, int hi, int k)
{
  double v, w;
  int L, R, i, j;

  for (L = lo, R = hi; L < R; ) {
    v = x[k];

    for(i = L, j = R; i <= j;) {
      while (rcmp(x[i], v) < 0)
        i++;

      while (rcmp(v, x[j]) < 0)
        j--;

      if (i <= j) {
        w = x[i];
        x[i++] = x[j];
        x[j--] = w;
      }
    }

    if (j < k)
      L = i;

    if (k < i)
      R = j;
  }
}


static void
rPsort(double *x, int n, int k)
{
  rPsort2(x, 0, n - 1, k);
}


static void
calclimits(
  double *x,
  int span,
  int n,
  int loopIndex,
  double mx,
  int *pMinidx,
  int *pMaxidx,
  double *pDmax
) {
  int minidx;
  int maxidx;
  double dmax;

  *pMinidx = 0;
  *pMaxidx = 0;
  *pDmax = 0.0;

  minidx = loopIndex;
  maxidx = loopIndex;

  while((maxidx - minidx) < (span-1)) {
    if (minidx <= 0) {
      /* Can only bump upper loopIndex */
      if (maxidx >= (n-1)) {
        break; /* can't go any further */
      }
      maxidx++;
    } else if (maxidx >= (n-1)) {
      minidx--;
    } else {
      /* Can go either way.  Pick the smallest increment */
      if ((x[maxidx+1] - mx) < (mx - x[minidx-1])) {
        maxidx++;
      } else {
        minidx--;
      }
    }
  }

  dmax = x[maxidx] - mx;

  if ((mx - x[minidx]) > dmax) {
    dmax = (mx - x[minidx]);
  }

  while (minidx > 0) {
    if ((mx - x[minidx-1]) <= dmax) {
      minidx--;
    } else {
      break;
    }
  }

  while (maxidx < (n - 1)) {
    if ((x[maxidx+1] - mx) <= dmax) {
      maxidx++;
    } else {
      break;
    }
  }

  *pMinidx = minidx;
  *pMaxidx = maxidx;
  *pDmax = dmax;
}


/*
 *  Least mean square fit to minimize the sum of
 *    (A*x) + (B*y) - C;
 *  solve for x, the weighted "1" coefficient
 *
 */
static int
calcleastmeansquare1(
  double *A,
  double *B,
  double *C,
  double *x,
  double *y,
  int minidx,
  int maxidx,
  int loopIndex
) {
  int index;
  double sqrA = 0.0;
  double sqrB = 0.0;
  double AB = 0.0;
  double AC = 0.0;
  double BC = 0.0;
  double D,E,F,G,H,I;
  double numerator;
  double denominator;
  double XVAL;
  double YVAL;

  *x = 0.0;
  *y = 0.0;

  for (index = minidx; index <= maxidx; index++) {
    sqrA += fsquare(A[index]);
    sqrB += fsquare(B[index]);
    AB += A[index]*B[index];
    AC += A[index]*C[index];
    BC += B[index]*C[index];
  }

  /*  We have (A*A*x) + (A*B*y) = AC;   = Dx+Ey=F;
      (A*B*x) + (B*B*y) = BC;   = Gx+Hy=I;
  */

  D = sqrA;
  E = AB;
  F = AC;
  G = AB;
  H = sqrB;
  I = BC;

  numerator = (F*H) - (I*E);
  denominator = (H*D) - (E*G);

  if (denominator == 0.0) {
    if (sqrA == 0.0) {
      return 1; /* Octave also returns a zero for this matrix inversion */
    } else {
      XVAL = AC/sqrA;
      *x = XVAL;
      return 1;
    }
  }

  XVAL = numerator/denominator;

  YVAL = (F - (D*XVAL))/E;

  *x = XVAL;
  *y = YVAL;
  return 1;
}


/*
 *  Least mean square fit to minimize the sum of
 *    (A*x) + (B*y) + (C*z) - D;
 *  solve for x, the weighted "1" coefficient
 *
 */
static int
calcleastmeansquare2(
  double *A,
  double *B,
  double *C,
  double *D,
  double *x,
  double *y,
  double *z,
  int minidx,
  int maxidx,
  int loopIndex
) {
  int index;
  double sqrA = 0.0;
  double sqrB = 0.0;
  double sqrC = 0.0;
  double AB = 0.0;
  double AC = 0.0;
  double AD = 0.0;
  double BC = 0.0;
  double BD = 0.0;
  double CD = 0.0;
  double E,F,G,H,I;
  double J,K,L,M,P,Q,R,S,T,U,V,W;
  double D2;
  double numerator;
  double denominator;
  double XVAL;
  double YVAL;
  double ZVAL;

  *x = 0.0;
  *y = 0.0;
  *z = 0.0;

  for (index = minidx; index <= maxidx; index++) {
    sqrA += fsquare(A[index]);
    sqrB += fsquare(B[index]);
    sqrC += fsquare(C[index]);
    AB += A[index]*B[index];
    AC += A[index]*C[index];
    AD += A[index]*D[index];
    BC += B[index]*C[index];
    BD += B[index]*D[index];
    CD += C[index]*D[index];
  }

  /*  We have (A*A*x) + (A*B*y) + (A*C*z) = AD;   = Jx+Ky+Lz=M;
              (A*B*x) + (B*B*y) + (B*C*z) = BD;   = Px+Qy+Rz=S;
        (A*C*x) + (B*C*y) + (C*C*z) = CD;   = Tx+Uy+Vz=W;
  */

  J = sqrA;
  K = AB;
  L = AC;
  M = AD;
  P = AB;
  Q = sqrB;
  R = BC;
  S = BD;
  T = AC;
  U = BC;
  V = sqrC;
  W = CD;

  /* Eliminate y
     (Q*J-K*P)x + (Q*L-K*R)z = (Q*M-K*S) = D2x + Ez = F;
     (U*P-Q*T)x + (U*R-Q*V)z = (U*S-Q*W) =  Gx + Hz = I;
  */
  D2 = (Q*J) - (K*P);
  E  = (Q*L) - (K*R);
  F  = (Q*M) - (K*S);
  G  = (U*P) - (Q*T);
  H  = (U*R) - (Q*V);
  I  = (U*S) - (Q*W);

  numerator = (F*H) - (I*E);
  denominator = (H*D2) - (E*G);

  if (denominator == 0.0) {
    if (sqrA == 0.0) {
      return 1; /* Octave also returns a zero for this matrix inversion */
    } else {
      XVAL = AC/sqrA;
      *x = XVAL;
      return 1;
    }
  }

  XVAL = numerator/denominator;
  ZVAL = (F - (D2*XVAL))/E;
  YVAL = (S - (P*XVAL) - (R*ZVAL))/Q;

  *x = XVAL;
  *y = YVAL;
  *z = ZVAL;
  return 1;
}


static double
calcmedian(double *x, double * y, int n, int nmax)
{
  int index;
  double returnVal;
  int m1;
  int m2;

  m1 = n/2;
  for (index = 0; index < n; index++) {
    y[index] = x[index];
  }
  rPsort(y,n,m1);
  returnVal = y[m1];
  if(n % 2 == 0) {
    m2 = n-m1-1;
    rPsort(y, n, m2);
    returnVal = (returnVal+y[m2])/2.0;
  }

  return returnVal;
}


static void
clowess(
  double *x,
  double *y,
  double *c,
  int n,
  int robust,
  int span,
  int iter,
  int loessFlag
) {
  double seps = sqrt(EPS);
  int k;
  double mx;
  int loopIndex;
  int index;
  double dmax;
  int minidx;
  int maxidx;
  int weightFlag;
  double *weight = NULL;
  double *v1 = NULL;
  double *v2 = NULL;
  double *v3 = NULL;
  double *y1 = NULL;
  double *r = NULL;
  double *r1 = NULL;
  double *r2 = NULL;
  double mad;
  double minmad;
  double *rweight = NULL;
  int result = 1;
  double xval;
  double yval;
  double zval = 0;

  weight = (double *)calloc(n,sizeof(double));
  v1 = (double *)calloc(n,sizeof(double));
  v2 = (double *)calloc(n,sizeof(double));
  v3 = (double *)calloc(n,sizeof(double));
  y1 = (double *)calloc(n,sizeof(double));

  if (robust) {
    r = (double *)calloc(n,sizeof(double));
    r1 = (double *)calloc(n,sizeof(double));
    r2 = (double *)calloc(n,sizeof(double));
    rweight = (double *)calloc(n,sizeof(double));
  }

  for (loopIndex = 0; loopIndex < n; loopIndex++) {
    if (loopIndex >= 1) {
      if (x[loopIndex] == x[loopIndex-1]) {
        c[loopIndex] = c[loopIndex-1];
        continue;
      }
    }

    mx = x[loopIndex];
    calclimits(x, span, n, loopIndex, mx, &minidx, &maxidx, &dmax);

    /* Compute weights */
    weightFlag = 1;

    for (index = minidx; index <= maxidx; index++) {
      weight[index] = sqrt(fcube(1.0 - fcube((fabs(x[index]-mx)/dmax))));
      if (weight[index] > seps) {
        weightFlag = 0;
      }
    }

    if (weightFlag == 1) {
      /* All weights are zero, skip weighting */
      for (index = minidx; index <= maxidx; index++) {
        weight[index] = 1.0;
      }
    }

    for (index = minidx; index <= maxidx; index++) {
      v1[index] = 1.0 * weight[index];
      v2[index] = (x[index]-mx) * weight[index];
      if (loessFlag) {
        v3[index] =  (x[index]-mx) * (x[index]-mx) * weight[index];
      }
      y1[index] = y[index] * weight[index];
    }

    if (loessFlag == 0) {
      result = calcleastmeansquare1(v1, v2, y1, &xval, &yval, minidx, maxidx, loopIndex);
    } else {
      result = calcleastmeansquare2(v1, v2, v3, y1, &xval, &yval, &zval, minidx, maxidx,loopIndex);
    }

    if (result == 0) {
      printf("ERROR: result is 0 for loopIndex %d,minidx %d maxidx %d\n",loopIndex,minidx,maxidx);
    }

    c[loopIndex] = xval;
  } /* End of primary loopIndex loop */

  if (robust) {
    minmad = 0.0;

    for (index = 0; index < n; index++) {
      if (fabs(y[index]) > minmad) {
        minmad = fabs(y[index]);
      }
    }

    minmad = minmad * EPS;

    for (k = 1; k <= iter; k++) {
      for (index = 0; index < n; index++) {
        r[index] = y[index]-c[index];
      }

      for (loopIndex = 0; loopIndex < n; loopIndex++) {
        if (loopIndex >= 1) {
          if (x[loopIndex] == x[loopIndex-1]) {
            c[loopIndex] = c[loopIndex-1];
            continue;
          }
        }

        mx = x[loopIndex];
        calclimits(x, span, n, loopIndex, mx, &minidx, &maxidx, &dmax);

        for (index = minidx; index <= maxidx; index++) {
          y1[index] = y[index];
          r1[index] = r[index];
        }

        weightFlag = 1;

        for (index = minidx; index <= maxidx; index++) {
          weight[index] = sqrt(fcube(1.0 - fcube((fabs(x[index]-mx)/dmax))));
          if (weight[index] > seps) {
            weightFlag = 0;
          }
        }

        if (weightFlag == 1) {
          /* All weights are zero, skip weighting */
          for (index = minidx; index <= maxidx; index++) {
            weight[index] = 1.0;
          }
        }

        mad = calcmedian(&r1[minidx],r2,maxidx-minidx+1,n);

        for (index = minidx; index <= maxidx; index++) {
          r1[index] = fabs(r1[index]-mad);
        }

        mad = calcmedian(&r1[minidx],r2,maxidx-minidx+1,n);

        if (mad > minmad) {
          for (index = minidx; index <= maxidx; index++) {
            rweight[index] = r1[index]/(6.0 * mad);

            if (rweight[index] <= 1.0) {
              rweight[index] = 1.0 - fsquare(rweight[index]);
            } else {
              rweight[index] = 0.0;
            }

            weight[index] = weight[index] * rweight[index];
          }
        }

        for (index = minidx; index <= maxidx; index++) {
          v1[index] = 1.0 * weight[index];
          v2[index] = (x[index]-mx) * weight[index];

          if (loessFlag) {
            v3[index] = (x[index]-mx) * (x[index]-mx) * weight[index];
          }

          y1[index] = y[index] * weight[index];
        }

        if (loessFlag == 0) {
          result = calcleastmeansquare1(v1,v2,y1,&xval,&yval,minidx,maxidx,loopIndex);
        } else {
          result = calcleastmeansquare2(v1,v2,v3,y1,&xval,&yval,&zval,minidx,maxidx,loopIndex);
        }

        if (result == 0) {
          printf("ERROR: result (2) is 0 for loopIndex %d iter %d\n",loopIndex,iter);
        }

        c[loopIndex] = xval;
      } /* End of loopIndex */
    } /* End of iteration liip */
  } /* End of robust iterations */

  free(weight);
  free(v1);
  free(v2);
  free(v3);
  free(y1);

  if (robust) {
    free(r1);
    free(r2);
    free(rweight);
    free(r);
  }
}


/*
 *   x and y are input values
 *   n is the size of the input arrays
 *   f is the fraction of points used to compute each fitted value
 *   nsteps is the number of iterations
 *   delta
 *   ys output fitted values as a function of x
 *   rw output robustness weights
 *   res output residuals
 */

DEFUN_DLD (lowess,args,,
           "Test lowess.\n\
\n\
Reference:\n\
\n\
  Lowess algorithm from R")
{
  /* function c = lowess(x,y, span, method, robust,iter) */
  octave_value xx = args(0);
  octave_value yy = args(1);
  octave_value_list retval;
  octave_value xspan = args(2);
  octave_value xmethod = args(3);
  octave_value xrobust = args(4);
  octave_value xiter = args(5);

  NDArray xxx = xx.array_value();
  NDArray yyy = yy.array_value();

  double spand = xspan.double_value();
  int span = int(spand + EPS);
  int iter = xiter.int_value();
  int robust = xrobust.int_value();
  std::string xxmethod = xmethod.string_value();
  const char * method = xxmethod.data();
  double *x = NULL;
  double *y = NULL;
  double *c = NULL;
  int index;
  int loessFlag = 0;
  int n = xxx.numel();

  if (robust != 1) {
    if (robust != 0) {
      printf("ERROR, robust %d is not zero or 1\n",robust);
      return retval;
    } else {
      iter = 0;
    }
  }

  if (strcmp(method,"loess") == 0) {
    loessFlag = 1;
  } else if (strcmp(method,"lowess") != 0) {
    printf("ERROR, method %s is not lowess\n",method);
    return retval;
  }

  x = xxx.fortran_vec();
  y = yyy.fortran_vec();

  if (n < span) {
    span = n;
  }

  if (span < 0) {
    printf("ERROR: Span must be an integer between 1 and length of x\n");
    return retval;
  }

  if (xxx.numel() != yyy.numel()) {
    printf("ERROR: x length %ld does not equal y length %ld\n",xxx.numel(),yyy.numel());
    return retval;
  }

  for (index = 0; index < n; index++) {
    if (isnan(y[index])) {
      printf("ERROR: Y contains NaN\n");
      return retval;
    }
  }

  for (index = 1; index < n; index++) {
    if (x[index] < x[index-1]) {
      printf("ERROR: X does not increase monotonically at index %d x %f\n",index,x[index]);
      return retval;
    }
  }

  c = (double *) calloc(n, sizeof(double));

  dim_vector dims = xxx.dims();

  NDArray zzz(dims);

  if (span == 1) {
    for (index = 0; index < n; index++) {
      zzz(index) = y[index];
    }

    retval(0) = zzz;
    free(c);
    printf("%s Returning with span == 1 \n", __FUNCTION__);
    return retval;
  }

  clowess(x, y, c, n, robust, span, iter, loessFlag);

  for (index = 0; index < n; index++) {
    zzz(index) = c[index];
  }

  free(c);
  retval(0) = zzz;
  return retval;
}

