/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/

#include "ode_common.h"
#include "ode_matrix.h"

// misc defines
#define ALLOCA dALLOCA16

dReal dDot (const dReal *a, const dReal *b, int n)
{  
  dReal p0,q0,m0,p1,q1,m1,sum;
  sum = 0;
  n -= 2;
  while (n >= 0) {
    p0 = a[0]; q0 = b[0];
    m0 = p0 * q0;
    p1 = a[1]; q1 = b[1];
    m1 = p1 * q1;
    sum += m0;
    sum += m1;
    a += 2;
    b += 2;
    n -= 2;
  }
  n += 2;
  while (n > 0) {
    sum += (*a) * (*b);
    a++;
    b++;
    n--;
  }
  return sum;
}


void dSetZero (dReal *a, int n)
{
  dAASSERT (a && n >= 0);
  while (n > 0) {
    *(a++) = 0;
    n--;
  }
}


void dSetValue (dReal *a, int n, dReal value)
{
  dAASSERT (a && n >= 0);
  while (n > 0) {
    *(a++) = value;
    n--;
  }
}


void dMultiply0 (dReal *A, const dReal *B, const dReal *C, int p, int q, int r)
{
  int i,j,k,qskip,rskip,rpad;
  dAASSERT (A && B && C && p>0 && q>0 && r>0);
  qskip = dPAD(q);
  rskip = dPAD(r);
  rpad = rskip - r;
  dReal sum;
  const dReal *b,*c,*bb;
  bb = B;
  for (i=p; i; i--) {
    for (j=0 ; j<r; j++) {
      c = C + j;
      b = bb;
      sum = 0;
      for (k=q; k; k--, c+=rskip) sum += (*(b++))*(*c);
      *(A++) = sum; 
    }
    A += rpad;
    bb += qskip;
  }
}


void dMultiply1 (dReal *A, const dReal *B, const dReal *C, int p, int q, int r)
{
  int i,j,k,pskip,rskip;
  dReal sum;
  dAASSERT (A && B && C && p>0 && q>0 && r>0);
  pskip = dPAD(p);
  rskip = dPAD(r);
  for (i=0; i<p; i++) {
    for (j=0; j<r; j++) {
      sum = 0;
      for (k=0; k<q; k++) sum += B[i+k*pskip] * C[j+k*rskip];
      A[i*rskip+j] = sum;
    }
  }
}


void dMultiply2 (dReal *A, const dReal *B, const dReal *C, int p, int q, int r)
{
  int i,j,k,z,rpad,qskip;
  dReal sum;
  const dReal *bb,*cc;
  dAASSERT (A && B && C && p>0 && q>0 && r>0);
  rpad = dPAD(r) - r;
  qskip = dPAD(q);
  bb = B;
  for (i=p; i; i--) {
    cc = C;
    for (j=r; j; j--) {
      z = 0;
      sum = 0;
      for (k=q; k; k--,z++) sum += bb[z] * cc[z];
      *(A++) = sum; 
      cc += qskip;
    }
    A += rpad;
    bb += qskip;
  }
}


int dFactorCholesky (dReal *A, int n)
{
  int i,j,k,nskip;
  dReal sum,*a,*b,*aa,*bb,*cc,*recip;
  dAASSERT (n > 0 && A);
  nskip = dPAD (n);
  recip = (dReal*) ALLOCA (n * sizeof(dReal));
  aa = A;
  for (i=0; i<n; i++) {
    bb = A;
    cc = A + i*nskip;
    for (j=0; j<i; j++) {
      sum = *cc;
      a = aa;
      b = bb;
      for (k=j; k; k--) sum -= (*(a++))*(*(b++));
      *cc = sum * recip[j];
      bb += nskip;
      cc++;
    }
    sum = *cc;
    a = aa;
    for (k=i; k; k--, a++) sum -= (*a)*(*a);
    if (sum <= REAL(0.0)) return 0;
    *cc = dSqrt(sum);
    recip[i] = dRecip (*cc);
    aa += nskip;
  }
  return 1;
}


void dSolveCholesky (const dReal *L, dReal *b, int n)
{
  int i,k,nskip;
  dReal sum,*y;
  dAASSERT (n > 0 && L && b);
  nskip = dPAD (n);
  y = (dReal*) ALLOCA (n*sizeof(dReal));
  for (i=0; i<n; i++) {
    sum = 0;
    for (k=0; k < i; k++) sum += L[i*nskip+k]*y[k];
    y[i] = (b[i]-sum)/L[i*nskip+i];
  }
  for (i=n-1; i >= 0; i--) {
    sum = 0;
    for (k=i+1; k < n; k++) sum += L[k*nskip+i]*b[k];
    b[i] = (y[i]-sum)/L[i*nskip+i];
  }
}


int dInvertPDMatrix (const dReal *A, dReal *Ainv, int n)
{
  int i,j,nskip;
  dReal *L,*x;
  dAASSERT (n > 0 && A && Ainv);
  nskip = dPAD (n);
  L = (dReal*) ALLOCA (nskip*n*sizeof(dReal));
  memcpy (L,A,nskip*n*sizeof(dReal));
  x = (dReal*) ALLOCA (n*sizeof(dReal));
  if (dFactorCholesky (L,n)==0) return 0;
  dSetZero (Ainv,n*nskip);	// make sure all padding elements set to 0
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) x[j] = 0;
    x[i] = 1;
    dSolveCholesky (L,x,n);
    for (j=0; j<n; j++) Ainv[j*nskip+i] = x[j];
  }
  return 1;  
}


int dIsPositiveDefinite (const dReal *A, int n)
{
  dReal *Acopy;
  dAASSERT (n > 0 && A);
  int nskip = dPAD (n);
  Acopy = (dReal*) ALLOCA (nskip*n * sizeof(dReal));
  memcpy (Acopy,A,nskip*n * sizeof(dReal));
  return dFactorCholesky (Acopy,n);
}


/***** this has been replaced by a faster version
void dSolveL1T (const dReal *L, dReal *b, int n, int nskip)
{
  int i,j;
  dAASSERT (L && b && n >= 0 && nskip >= n);
  dReal sum;
  for (i=n-2; i>=0; i--) {
    sum = 0;
    for (j=i+1; j<n; j++) sum += L[j*nskip+i]*b[j];
    b[i] -= sum;
  }
}
*/


void dVectorScale (dReal *a, const dReal *d, int n)
{
  dAASSERT (a && d && n >= 0);
  for (int i=0; i<n; i++) a[i] *= d[i];
}


void dSolveLDLT (const dReal *L, const dReal *d, dReal *b, int n, int nskip)
{
  dAASSERT (L && d && b && n > 0 && nskip >= n);
  dSolveL1 (L,b,n,nskip);
  dVectorScale (b,d,n);
  dSolveL1T (L,b,n,nskip);
}


void dLDLTAddTL (dReal *L, dReal *d, const dReal *a, int n, int nskip)
{
  int j,p;
  dReal *W1,*W2,W11,W21,alpha1,alpha2,alphanew,gamma1,gamma2,k1,k2,Wp,ell,dee;
  dAASSERT (L && d && a && n > 0 && nskip >= n);

  if (n < 2) return;
  W1 = (dReal*) ALLOCA (n*sizeof(dReal));
  W2 = (dReal*) ALLOCA (n*sizeof(dReal));

  W1[0] = 0;
  W2[0] = 0;
  for (j=1; j<n; j++) W1[j] = W2[j] = a[j] * M_SQRT1_2;
  W11 = (REAL(0.5)*a[0]+1)*M_SQRT1_2;
  W21 = (REAL(0.5)*a[0]-1)*M_SQRT1_2;

  alpha1=1;
  alpha2=1;

  dee = d[0];
  alphanew = alpha1 + (W11*W11)*dee;
  dee /= alphanew;
  gamma1 = W11 * dee;
  dee *= alpha1;
  alpha1 = alphanew;
  alphanew = alpha2 - (W21*W21)*dee;
  dee /= alphanew;
  gamma2 = W21 * dee;
  alpha2 = alphanew;
  k1 = REAL(1.0) - W21*gamma1;
  k2 = W21*gamma1*W11 - W21;
  for (p=1; p<n; p++) {
    Wp = W1[p];
    ell = L[p*nskip];
    W1[p] =    Wp - W11*ell;
    W2[p] = k1*Wp +  k2*ell;
  }

  for (j=1; j<n; j++) {
    dee = d[j];
    alphanew = alpha1 + (W1[j]*W1[j])*dee;
    dee /= alphanew;
    gamma1 = W1[j] * dee;
    dee *= alpha1;
    alpha1 = alphanew;
    alphanew = alpha2 - (W2[j]*W2[j])*dee;
    dee /= alphanew;
    gamma2 = W2[j] * dee;
    dee *= alpha2;
    d[j] = dee;
    alpha2 = alphanew;

    k1 = W1[j];
    k2 = W2[j];
    for (p=j+1; p<n; p++) {
      ell = L[p*nskip+j];
      Wp = W1[p] - k1 * ell;
      ell += gamma1 * Wp;
      W1[p] = Wp;
      Wp = W2[p] - k2 * ell;
      ell -= gamma2 * Wp;
      W2[p] = Wp;
      L[p*nskip+j] = ell;
    }
  }
}


// macros for dLDLTRemove() for accessing A - either access the matrix
// directly or access it via row pointers. we are only supposed to reference
// the lower triangle of A (it is symmetric), but indexes i and j come from
// permutation vectors so they are not predictable. so do a test on the
// indexes - this should not slow things down too much, as we don't do this
// in an inner loop.

#define _GETA(i,j) (A[i][j])
//#define _GETA(i,j) (A[(i)*nskip+(j)])
#define GETA(i,j) ((i > j) ? _GETA(i,j) : _GETA(j,i))


void dLDLTRemove (dReal **A, const int *p, dReal *L, dReal *d,
		  int n1, int n2, int r, int nskip)
{
  int i;
  dAASSERT(A && p && L && d && n1 > 0 && n2 > 0 && r >= 0 && r < n2 &&
	   n1 >= n2 && nskip >= n1);
  #ifndef dNODEBUG
  for (i=0; i<n2; i++) dIASSERT(p[i] >= 0 && p[i] < n1);
  #endif

  if (r==n2-1) {
    return;		// deleting last row/col is easy
  }
  else if (r==0) {
    dReal *a = (dReal*) ALLOCA (n2 * sizeof(dReal));
    for (i=0; i<n2; i++) a[i] = -GETA(p[i],p[0]);
    a[0] += REAL(1.0);
    dLDLTAddTL (L,d,a,n2,nskip);
  }
  else {
    dReal *t = (dReal*) ALLOCA (r * sizeof(dReal));
    dReal *a = (dReal*) ALLOCA ((n2-r) * sizeof(dReal));
    for (i=0; i<r; i++) t[i] = L[r*nskip+i] / d[i];
    for (i=0; i<(n2-r); i++)
      a[i] = dDot(L+(r+i)*nskip,t,r) - GETA(p[r+i],p[r]);
    a[0] += REAL(1.0);
    dLDLTAddTL (L + r*nskip+r, d+r, a, n2-r, nskip);
  }

  // snip out row/column r from L and d
  dRemoveRowCol (L,n2,nskip,r);
  if (r < (n2-1)) memmove (d+r,d+r+1,(n2-r-1)*sizeof(dReal));
}


void dRemoveRowCol (dReal *A, int n, int nskip, int r)
{
  int i;
  dAASSERT(A && n > 0 && nskip >= n && r >= 0 && r < n);
  if (r >= n-1) return;
  if (r > 0) {
    for (i=0; i<r; i++)
      memmove (A+i*nskip+r,A+i*nskip+r+1,(n-r-1)*sizeof(dReal));
    for (i=r; i<(n-1); i++)
      memcpy (A+i*nskip,A+i*nskip+nskip,r*sizeof(dReal));
  }
  for (i=r; i<(n-1); i++)
    memcpy (A+i*nskip+r,A+i*nskip+nskip+r+1,(n-r-1)*sizeof(dReal));
}

/* solve L*X=B, with B containing 1 right hand sides.
 * L is an n*n lower triangular matrix with ones on the diagonal.
 * L is stored by rows and its leading dimension is lskip.
 * B is an n*1 matrix that contains the right hand sides.
 * B is stored by columns and its leading dimension is also lskip.
 * B is overwritten with X.
 * this processes blocks of 2*2.
 * if this is in the factorizer source file, n must be a multiple of 2.
 */

static void dSolveL1_1 (const dReal *L, dReal *B, int n, int lskip1)
{  
  /* declare variables - Z matrix, p and q vectors, etc */
  dReal Z11,m11,Z21,m21,p1,q1,p2,*ex;
  const dReal *ell;
  int i,j;
  /* compute all 2 x 1 blocks of X */
  for (i=0; i < n; i+=2) {
    /* compute all 2 x 1 block of X, from rows i..i+2-1 */
    /* set the Z matrix to 0 */
    Z11=0;
    Z21=0;
    ell = L + i*lskip1;
    ex = B;
    /* the inner loop that computes outer products and adds them to Z */
    for (j=i-2; j >= 0; j -= 2) {
      /* compute outer product and add it to the Z matrix */
      p1=ell[0];
      q1=ex[0];
      m11 = p1 * q1;
      p2=ell[lskip1];
      m21 = p2 * q1;
      Z11 += m11;
      Z21 += m21;
      /* compute outer product and add it to the Z matrix */
      p1=ell[1];
      q1=ex[1];
      m11 = p1 * q1;
      p2=ell[1+lskip1];
      m21 = p2 * q1;
      /* advance pointers */
      ell += 2;
      ex += 2;
      Z11 += m11;
      Z21 += m21;
      /* end of inner loop */
    }
    /* compute left-over iterations */
    j += 2;
    for (; j > 0; j--) {
      /* compute outer product and add it to the Z matrix */
      p1=ell[0];
      q1=ex[0];
      m11 = p1 * q1;
      p2=ell[lskip1];
      m21 = p2 * q1;
      /* advance pointers */
      ell += 1;
      ex += 1;
      Z11 += m11;
      Z21 += m21;
    }
    /* finish computing the X(i) block */
    Z11 = ex[0] - Z11;
    ex[0] = Z11;
    p1 = ell[lskip1];
    Z21 = ex[1] - Z21 - p1*Z11;
    ex[1] = Z21;
    /* end of outer loop */
  }
}

/* solve L*X=B, with B containing 2 right hand sides.
 * L is an n*n lower triangular matrix with ones on the diagonal.
 * L is stored by rows and its leading dimension is lskip.
 * B is an n*2 matrix that contains the right hand sides.
 * B is stored by columns and its leading dimension is also lskip.
 * B is overwritten with X.
 * this processes blocks of 2*2.
 * if this is in the factorizer source file, n must be a multiple of 2.
 */

static void dSolveL1_2 (const dReal *L, dReal *B, int n, int lskip1)
{  
  /* declare variables - Z matrix, p and q vectors, etc */
  dReal Z11,m11,Z12,m12,Z21,m21,Z22,m22,p1,q1,p2,q2,*ex;
  const dReal *ell;
  int i,j;
  /* compute all 2 x 2 blocks of X */
  for (i=0; i < n; i+=2) {
    /* compute all 2 x 2 block of X, from rows i..i+2-1 */
    /* set the Z matrix to 0 */
    Z11=0;
    Z12=0;
    Z21=0;
    Z22=0;
    ell = L + i*lskip1;
    ex = B;
    /* the inner loop that computes outer products and adds them to Z */
    for (j=i-2; j >= 0; j -= 2) {
      /* compute outer product and add it to the Z matrix */
      p1=ell[0];
      q1=ex[0];
      m11 = p1 * q1;
      q2=ex[lskip1];
      m12 = p1 * q2;
      p2=ell[lskip1];
      m21 = p2 * q1;
      m22 = p2 * q2;
      Z11 += m11;
      Z12 += m12;
      Z21 += m21;
      Z22 += m22;
      /* compute outer product and add it to the Z matrix */
      p1=ell[1];
      q1=ex[1];
      m11 = p1 * q1;
      q2=ex[1+lskip1];
      m12 = p1 * q2;
      p2=ell[1+lskip1];
      m21 = p2 * q1;
      m22 = p2 * q2;
      /* advance pointers */
      ell += 2;
      ex += 2;
      Z11 += m11;
      Z12 += m12;
      Z21 += m21;
      Z22 += m22;
      /* end of inner loop */
    }
    /* compute left-over iterations */
    j += 2;
    for (; j > 0; j--) {
      /* compute outer product and add it to the Z matrix */
      p1=ell[0];
      q1=ex[0];
      m11 = p1 * q1;
      q2=ex[lskip1];
      m12 = p1 * q2;
      p2=ell[lskip1];
      m21 = p2 * q1;
      m22 = p2 * q2;
      /* advance pointers */
      ell += 1;
      ex += 1;
      Z11 += m11;
      Z12 += m12;
      Z21 += m21;
      Z22 += m22;
    }
    /* finish computing the X(i) block */
    Z11 = ex[0] - Z11;
    ex[0] = Z11;
    Z12 = ex[lskip1] - Z12;
    ex[lskip1] = Z12;
    p1 = ell[lskip1];
    Z21 = ex[1] - Z21 - p1*Z11;
    ex[1] = Z21;
    Z22 = ex[1+lskip1] - Z22 - p1*Z12;
    ex[1+lskip1] = Z22;
    /* end of outer loop */
  }
}



void dFactorLDLT (dReal *A, dReal *d, int n, int nskip1)
{  
  int i,j;
  dReal sum,*ell,*dee,dd,p1,p2,q1,q2,Z11,m11,Z21,m21,Z22,m22;
  if (n < 1) return;
  
  for (i=0; i<=n-2; i += 2) {
    /* solve L*(D*l)=a, l is scaled elements in 2 x i block at A(i,0) */
    dSolveL1_2 (A,A+i*nskip1,i,nskip1);
    /* scale the elements in a 2 x i block at A(i,0), and also */
    /* compute Z = the outer product matrix that we'll need. */
    Z11 = 0;
    Z21 = 0;
    Z22 = 0;
    ell = A+i*nskip1;
    dee = d;
    for (j=i-6; j >= 0; j -= 6) {
      p1 = ell[0];
      p2 = ell[nskip1];
      dd = dee[0];
      q1 = p1*dd;
      q2 = p2*dd;
      ell[0] = q1;
      ell[nskip1] = q2;
      m11 = p1*q1;
      m21 = p2*q1;
      m22 = p2*q2;
      Z11 += m11;
      Z21 += m21;
      Z22 += m22;
      p1 = ell[1];
      p2 = ell[1+nskip1];
      dd = dee[1];
      q1 = p1*dd;
      q2 = p2*dd;
      ell[1] = q1;
      ell[1+nskip1] = q2;
      m11 = p1*q1;
      m21 = p2*q1;
      m22 = p2*q2;
      Z11 += m11;
      Z21 += m21;
      Z22 += m22;
      p1 = ell[2];
      p2 = ell[2+nskip1];
      dd = dee[2];
      q1 = p1*dd;
      q2 = p2*dd;
      ell[2] = q1;
      ell[2+nskip1] = q2;
      m11 = p1*q1;
      m21 = p2*q1;
      m22 = p2*q2;
      Z11 += m11;
      Z21 += m21;
      Z22 += m22;
      p1 = ell[3];
      p2 = ell[3+nskip1];
      dd = dee[3];
      q1 = p1*dd;
      q2 = p2*dd;
      ell[3] = q1;
      ell[3+nskip1] = q2;
      m11 = p1*q1;
      m21 = p2*q1;
      m22 = p2*q2;
      Z11 += m11;
      Z21 += m21;
      Z22 += m22;
      p1 = ell[4];
      p2 = ell[4+nskip1];
      dd = dee[4];
      q1 = p1*dd;
      q2 = p2*dd;
      ell[4] = q1;
      ell[4+nskip1] = q2;
      m11 = p1*q1;
      m21 = p2*q1;
      m22 = p2*q2;
      Z11 += m11;
      Z21 += m21;
      Z22 += m22;
      p1 = ell[5];
      p2 = ell[5+nskip1];
      dd = dee[5];
      q1 = p1*dd;
      q2 = p2*dd;
      ell[5] = q1;
      ell[5+nskip1] = q2;
      m11 = p1*q1;
      m21 = p2*q1;
      m22 = p2*q2;
      Z11 += m11;
      Z21 += m21;
      Z22 += m22;
      ell += 6;
      dee += 6;
    }
    /* compute left-over iterations */
    j += 6;
    for (; j > 0; j--) {
      p1 = ell[0];
      p2 = ell[nskip1];
      dd = dee[0];
      q1 = p1*dd;
      q2 = p2*dd;
      ell[0] = q1;
      ell[nskip1] = q2;
      m11 = p1*q1;
      m21 = p2*q1;
      m22 = p2*q2;
      Z11 += m11;
      Z21 += m21;
      Z22 += m22;
      ell++;
      dee++;
    }
    /* solve for diagonal 2 x 2 block at A(i,i) */
    Z11 = ell[0] - Z11;
    Z21 = ell[nskip1] - Z21;
    Z22 = ell[1+nskip1] - Z22;
    dee = d + i;
    /* factorize 2 x 2 block Z,dee */
    /* factorize row 1 */
    dee[0] = dRecip(Z11);
    /* factorize row 2 */
    sum = 0;
    q1 = Z21;
    q2 = q1 * dee[0];
    Z21 = q2;
    sum += q1*q2;
    dee[1] = dRecip(Z22 - sum);
    /* done factorizing 2 x 2 block */
    ell[nskip1] = Z21;
  }
  /* compute the (less than 2) rows at the bottom */
  switch (n-i) {
    case 0:
    break;
    
    case 1:
    dSolveL1_1 (A,A+i*nskip1,i,nskip1);
    /* scale the elements in a 1 x i block at A(i,0), and also */
    /* compute Z = the outer product matrix that we'll need. */
    Z11 = 0;
    ell = A+i*nskip1;
    dee = d;
    for (j=i-6; j >= 0; j -= 6) {
      p1 = ell[0];
      dd = dee[0];
      q1 = p1*dd;
      ell[0] = q1;
      m11 = p1*q1;
      Z11 += m11;
      p1 = ell[1];
      dd = dee[1];
      q1 = p1*dd;
      ell[1] = q1;
      m11 = p1*q1;
      Z11 += m11;
      p1 = ell[2];
      dd = dee[2];
      q1 = p1*dd;
      ell[2] = q1;
      m11 = p1*q1;
      Z11 += m11;
      p1 = ell[3];
      dd = dee[3];
      q1 = p1*dd;
      ell[3] = q1;
      m11 = p1*q1;
      Z11 += m11;
      p1 = ell[4];
      dd = dee[4];
      q1 = p1*dd;
      ell[4] = q1;
      m11 = p1*q1;
      Z11 += m11;
      p1 = ell[5];
      dd = dee[5];
      q1 = p1*dd;
      ell[5] = q1;
      m11 = p1*q1;
      Z11 += m11;
      ell += 6;
      dee += 6;
    }
    /* compute left-over iterations */
    j += 6;
    for (; j > 0; j--) {
      p1 = ell[0];
      dd = dee[0];
      q1 = p1*dd;
      ell[0] = q1;
      m11 = p1*q1;
      Z11 += m11;
      ell++;
      dee++;
    }
    /* solve for diagonal 1 x 1 block at A(i,i) */
    Z11 = ell[0] - Z11;
    dee = d + i;
    /* factorize 1 x 1 block Z,dee */
    /* factorize row 1 */
    dee[0] = dRecip(Z11);
    /* done factorizing 1 x 1 block */
    break;
    
    default: *((char*)0)=0;  /* this should never happen! */
  }
}

/* solve L*X=B, with B containing 1 right hand sides.
 * L is an n*n lower triangular matrix with ones on the diagonal.
 * L is stored by rows and its leading dimension is lskip.
 * B is an n*1 matrix that contains the right hand sides.
 * B is stored by columns and its leading dimension is also lskip.
 * B is overwritten with X.
 * this processes blocks of 4*4.
 * if this is in the factorizer source file, n must be a multiple of 4.
 */

void dSolveL1 (const dReal *L, dReal *B, int n, int lskip1)
{  
  /* declare variables - Z matrix, p and q vectors, etc */
  dReal Z11,Z21,Z31,Z41,p1,q1,p2,p3,p4,*ex;
  const dReal *ell;
  int lskip2,lskip3,i,j;
  /* compute lskip values */
  lskip2 = 2*lskip1;
  lskip3 = 3*lskip1;
  /* compute all 4 x 1 blocks of X */
  for (i=0; i <= n-4; i+=4) {
    /* compute all 4 x 1 block of X, from rows i..i+4-1 */
    /* set the Z matrix to 0 */
    Z11=0;
    Z21=0;
    Z31=0;
    Z41=0;
    ell = L + i*lskip1;
    ex = B;
    /* the inner loop that computes outer products and adds them to Z */
    for (j=i-12; j >= 0; j -= 12) {
      /* load p and q values */
      p1=ell[0];
      q1=ex[0];
      p2=ell[lskip1];
      p3=ell[lskip2];
      p4=ell[lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[1];
      q1=ex[1];
      p2=ell[1+lskip1];
      p3=ell[1+lskip2];
      p4=ell[1+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[2];
      q1=ex[2];
      p2=ell[2+lskip1];
      p3=ell[2+lskip2];
      p4=ell[2+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[3];
      q1=ex[3];
      p2=ell[3+lskip1];
      p3=ell[3+lskip2];
      p4=ell[3+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[4];
      q1=ex[4];
      p2=ell[4+lskip1];
      p3=ell[4+lskip2];
      p4=ell[4+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[5];
      q1=ex[5];
      p2=ell[5+lskip1];
      p3=ell[5+lskip2];
      p4=ell[5+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[6];
      q1=ex[6];
      p2=ell[6+lskip1];
      p3=ell[6+lskip2];
      p4=ell[6+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[7];
      q1=ex[7];
      p2=ell[7+lskip1];
      p3=ell[7+lskip2];
      p4=ell[7+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[8];
      q1=ex[8];
      p2=ell[8+lskip1];
      p3=ell[8+lskip2];
      p4=ell[8+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[9];
      q1=ex[9];
      p2=ell[9+lskip1];
      p3=ell[9+lskip2];
      p4=ell[9+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[10];
      q1=ex[10];
      p2=ell[10+lskip1];
      p3=ell[10+lskip2];
      p4=ell[10+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* load p and q values */
      p1=ell[11];
      q1=ex[11];
      p2=ell[11+lskip1];
      p3=ell[11+lskip2];
      p4=ell[11+lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* advance pointers */
      ell += 12;
      ex += 12;
      /* end of inner loop */
    }
    /* compute left-over iterations */
    j += 12;
    for (; j > 0; j--) {
      /* load p and q values */
      p1=ell[0];
      q1=ex[0];
      p2=ell[lskip1];
      p3=ell[lskip2];
      p4=ell[lskip3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      Z21 += p2 * q1;
      Z31 += p3 * q1;
      Z41 += p4 * q1;
      /* advance pointers */
      ell += 1;
      ex += 1;
    }
    /* finish computing the X(i) block */
    Z11 = ex[0] - Z11;
    ex[0] = Z11;
    p1 = ell[lskip1];
    Z21 = ex[1] - Z21 - p1*Z11;
    ex[1] = Z21;
    p1 = ell[lskip2];
    p2 = ell[1+lskip2];
    Z31 = ex[2] - Z31 - p1*Z11 - p2*Z21;
    ex[2] = Z31;
    p1 = ell[lskip3];
    p2 = ell[1+lskip3];
    p3 = ell[2+lskip3];
    Z41 = ex[3] - Z41 - p1*Z11 - p2*Z21 - p3*Z31;
    ex[3] = Z41;
    /* end of outer loop */
  }
  /* compute rows at end that are not a multiple of block size */
  for (; i < n; i++) {
    /* compute all 1 x 1 block of X, from rows i..i+1-1 */
    /* set the Z matrix to 0 */
    Z11=0;
    ell = L + i*lskip1;
    ex = B;
    /* the inner loop that computes outer products and adds them to Z */
    for (j=i-12; j >= 0; j -= 12) {
      /* load p and q values */
      p1=ell[0];
      q1=ex[0];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[1];
      q1=ex[1];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[2];
      q1=ex[2];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[3];
      q1=ex[3];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[4];
      q1=ex[4];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[5];
      q1=ex[5];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[6];
      q1=ex[6];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[7];
      q1=ex[7];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[8];
      q1=ex[8];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[9];
      q1=ex[9];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[10];
      q1=ex[10];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* load p and q values */
      p1=ell[11];
      q1=ex[11];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* advance pointers */
      ell += 12;
      ex += 12;
      /* end of inner loop */
    }
    /* compute left-over iterations */
    j += 12;
    for (; j > 0; j--) {
      /* load p and q values */
      p1=ell[0];
      q1=ex[0];
      /* compute outer product and add it to the Z matrix */
      Z11 += p1 * q1;
      /* advance pointers */
      ell += 1;
      ex += 1;
    }
    /* finish computing the X(i) block */
    Z11 = ex[0] - Z11;
    ex[0] = Z11;
  }
}

/* solve L^T * x=b, with b containing 1 right hand side.
 * L is an n*n lower triangular matrix with ones on the diagonal.
 * L is stored by rows and its leading dimension is lskip.
 * b is an n*1 matrix that contains the right hand side.
 * b is overwritten with x.
 * this processes blocks of 4.
 */

void dSolveL1T (const dReal *L, dReal *B, int n, int lskip1)
{  
  /* declare variables - Z matrix, p and q vectors, etc */
  dReal Z11,m11,Z21,m21,Z31,m31,Z41,m41,p1,q1,p2,p3,p4,*ex;
  const dReal *ell;
  int lskip2,lskip3,i,j;
  /* special handling for L and B because we're solving L1 *transpose* */
  L = L + (n-1)*(lskip1+1);
  B = B + n-1;
  lskip1 = -lskip1;
  /* compute lskip values */
  lskip2 = 2*lskip1;
  lskip3 = 3*lskip1;
  /* compute all 4 x 1 blocks of X */
  for (i=0; i <= n-4; i+=4) {
    /* compute all 4 x 1 block of X, from rows i..i+4-1 */
    /* set the Z matrix to 0 */
    Z11=0;
    Z21=0;
    Z31=0;
    Z41=0;
    ell = L - i;
    ex = B;
    /* the inner loop that computes outer products and adds them to Z */
    for (j=i-4; j >= 0; j -= 4) {
      /* load p and q values */
      p1=ell[0];
      q1=ex[0];
      p2=ell[-1];
      p3=ell[-2];
      p4=ell[-3];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      m21 = p2 * q1;
      m31 = p3 * q1;
      m41 = p4 * q1;
      ell += lskip1;
      Z11 += m11;
      Z21 += m21;
      Z31 += m31;
      Z41 += m41;
      /* load p and q values */
      p1=ell[0];
      q1=ex[-1];
      p2=ell[-1];
      p3=ell[-2];
      p4=ell[-3];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      m21 = p2 * q1;
      m31 = p3 * q1;
      m41 = p4 * q1;
      ell += lskip1;
      Z11 += m11;
      Z21 += m21;
      Z31 += m31;
      Z41 += m41;
      /* load p and q values */
      p1=ell[0];
      q1=ex[-2];
      p2=ell[-1];
      p3=ell[-2];
      p4=ell[-3];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      m21 = p2 * q1;
      m31 = p3 * q1;
      m41 = p4 * q1;
      ell += lskip1;
      Z11 += m11;
      Z21 += m21;
      Z31 += m31;
      Z41 += m41;
      /* load p and q values */
      p1=ell[0];
      q1=ex[-3];
      p2=ell[-1];
      p3=ell[-2];
      p4=ell[-3];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      m21 = p2 * q1;
      m31 = p3 * q1;
      m41 = p4 * q1;
      ell += lskip1;
      ex -= 4;
      Z11 += m11;
      Z21 += m21;
      Z31 += m31;
      Z41 += m41;
      /* end of inner loop */
    }
    /* compute left-over iterations */
    j += 4;
    for (; j > 0; j--) {
      /* load p and q values */
      p1=ell[0];
      q1=ex[0];
      p2=ell[-1];
      p3=ell[-2];
      p4=ell[-3];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      m21 = p2 * q1;
      m31 = p3 * q1;
      m41 = p4 * q1;
      ell += lskip1;
      ex -= 1;
      Z11 += m11;
      Z21 += m21;
      Z31 += m31;
      Z41 += m41;
    }
    /* finish computing the X(i) block */
    Z11 = ex[0] - Z11;
    ex[0] = Z11;
    p1 = ell[-1];
    Z21 = ex[-1] - Z21 - p1*Z11;
    ex[-1] = Z21;
    p1 = ell[-2];
    p2 = ell[-2+lskip1];
    Z31 = ex[-2] - Z31 - p1*Z11 - p2*Z21;
    ex[-2] = Z31;
    p1 = ell[-3];
    p2 = ell[-3+lskip1];
    p3 = ell[-3+lskip2];
    Z41 = ex[-3] - Z41 - p1*Z11 - p2*Z21 - p3*Z31;
    ex[-3] = Z41;
    /* end of outer loop */
  }
  /* compute rows at end that are not a multiple of block size */
  for (; i < n; i++) {
    /* compute all 1 x 1 block of X, from rows i..i+1-1 */
    /* set the Z matrix to 0 */
    Z11=0;
    ell = L - i;
    ex = B;
    /* the inner loop that computes outer products and adds them to Z */
    for (j=i-4; j >= 0; j -= 4) {
      /* load p and q values */
      p1=ell[0];
      q1=ex[0];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      ell += lskip1;
      Z11 += m11;
      /* load p and q values */
      p1=ell[0];
      q1=ex[-1];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      ell += lskip1;
      Z11 += m11;
      /* load p and q values */
      p1=ell[0];
      q1=ex[-2];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      ell += lskip1;
      Z11 += m11;
      /* load p and q values */
      p1=ell[0];
      q1=ex[-3];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      ell += lskip1;
      ex -= 4;
      Z11 += m11;
      /* end of inner loop */
    }
    /* compute left-over iterations */
    j += 4;
    for (; j > 0; j--) {
      /* load p and q values */
      p1=ell[0];
      q1=ex[0];
      /* compute outer product and add it to the Z matrix */
      m11 = p1 * q1;
      ell += lskip1;
      ex -= 1;
      Z11 += m11;
    }
    /* finish computing the X(i) block */
    Z11 = ex[0] - Z11;
    ex[0] = Z11;
  }
}
