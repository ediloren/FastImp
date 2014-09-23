/*
Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA.
All rights reserved.

This Agreement gives you, the LICENSEE, certain rights and obligations.
By using the software, you indicate that you have read, understood, and
will comply with the terms.

Permission to use, copy and modify for internal, noncommercial purposes
is hereby granted.  Any distribution of this program or any part thereof
is strictly prohibited without prior written consent of M.I.T.

Title to copyright to this software and to any associated documentation
shall at all times remain with M.I.T. and LICENSEE agrees to preserve
same.  LICENSEE agrees not to make any copies except for LICENSEE'S
internal noncommercial use, or to use separately any portion of this
software without prior written consent of M.I.T.  LICENSEE agrees to
place the appropriate copyright notice on any such copies.

Nothing in this Agreement shall be construed as conferring rights to use
in advertising, publicity or otherwise any trademark or the name of
"Massachusetts Institute of Technology" or "M.I.T."

M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By
way of example, but not limitation, M.I.T. MAKES NO REPRESENTATIONS OR
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR
THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS OR DOCUMENTATION WILL
NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
M.I.T. shall not be held liable for any liability nor for any direct,
indirect or consequential damages with respect to any claim by LICENSEE
or any third party on account of or arising from this Agreement or use
of this software.
*/
/* FILE: dense.c
* 
* Routines for definitions of and operations on/with dense
* matrices, e.g. LU-factorization and back-solves for direct
* solutions of problems.
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include "dense.h"
/* For interaction with cLapack routine: */
#include "dgesvd.h"

#define MAX(a,b) ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b) ( ((a)>(b)) ? (b) : (a) )
#define ABS(A) ( ( (A) > 0 ) ? (A) : (-(A)) )
#define FREE(p) { assert((p)!=NULL); free(p); (p)=NULL;}
#define SQDEX(i, j, siz) ((i)*(siz) + (j))

/*
* ==================================================================== 
* Make a (pseudo)inverse of a dense matrix using CLAPACK SVD.
*
* Note the different definitions of a matrix here (double **) and 
* in other routines (double *).
*
* The matrix A must be of size at least [max(rows,cols)][max(rows,cols)] 
* since the pseudoinverse (and not its transpose) is returned in A.
*
* Singular value decomposition is used to calculate the psuedo-inverse.
* If #rows=#cols then (an approximation of) the inverse is found. 
* Otherwise an approximation of the pseudoinverse is obtained.
*
* The return value of this routine is the estimated rank of the system.
* ==================================================================== */ 
int 
pseudo_inverse(
	       double **A, 
	       int rows, 
	       int cols) 
{
  char fctName[] = "pinv_new";
  int i,j,k;
  int nsv, rank;
  double tol;
  /* Variables needed for interaction with CLAPACK routines 
   * (FORTRAN style)                                                   */
  double *amat, *svals, *U,  *V, *work;
  long int lda,        ldu, ldvt, lwork, m,n,info;

  /* Short-hands for this routine only [undef'ed at end of routine]    */
#define AA(i,j)  amat[i + j * (int) lda]
#define UU(i,j)  U[i + j * (int) ldu]
#define VV(i,j)  V[i + j * (int) ldvt]

  /* Set up arrays for call to CLAPACK routine                         */
  /* Use leading dimensions large enough to calculate both under- 
   * and over-determined systems of equations. The way this is 
   * calculated, we need to store U*(Z^-1^T) in U. The array U must 
   * have at least cols columns to perform this operation, and at 
   * least rows columns to store the initial U. However, for an under-
   * determined system the last cols-rows columns of U*(Z^-1^T) will be 
   * zero. Exploiting this, the size of the array U is just mxm        */
  m = (long int) rows;
  n = (long int) cols;
  lda  = m;
  ldu  = m;
  ldvt = n;

  nsv = MIN(rows,cols);
  /* Allocate memory for:
   * amat  : The matrix to factor
   * svals : Vector of singular values 
   * U     : mxm unitary matrix
   * V     : nxn unitary matrix                                        */
  amat  = (double *) calloc((int) m*n+1, sizeof(double)); 
  svals = (double *) calloc((int) m+n,   sizeof(double)); 
  U     = (double *) calloc((int) ldu*m+1,sizeof(double));
  V     = (double *) calloc((int) ldvt*n+1, sizeof(double));
  /* In matrix form amat = U * diag(svals) * V^T                       */
  /* Leading dimensions [number of rows] of the above arrays           */
  /* Work storage */
  lwork = 100*(m+n);
  work  = (double *) calloc(lwork, sizeof(double));

  /* Copy A to temporary storage                                       */
  for (j=0; j<cols; j++) {
    for (i=0; i<rows; i++) {
      AA(i,j) = A[i][j];
    }
  }

  {
    char jobu,jobvt;
    jobu  = 'A'; /* Calculate all columns of matrix U                  */
    jobvt = 'A'; /* Calculate all columns of matrix V^T                */

    /* Note that dgesvd_ returns V^T rather than V                     */
    dgesvd_(&jobu, &jobvt, &m, &n, 
	    amat, &lda, svals, U, &ldu, V, &ldvt, 
	    work, &lwork, &info);

    /* Test info from dgesvd */
    if (info) {
      if (info<0) {
	printf("%s: ERROR: clapack routine 'dgesvd_' complained about\n"
	       "    illegal value of argument # %d\n",fctName,(int) -info);
	exit(1);
      }
      else{
	printf("%s: ERROR: clapack routine 'dgesvd_' complained that\n"
	       "    %d superdiagonals didn't converge.\n",
	       fctName,(int) info);
	exit(1);
      }
    }
  }

  /* Test the the singular values are returned in correct ordering 
   * (This is done because there were problems with this with a former 
   * implementation [using f2c'ed linpack] when high optimization 
   * was used)                                                         */
  if (svals[0]<0.0) {
    printf("%s: ERROR: First singular value returned by clapack \n"
	   "   is negative: %16.6e.\n",fctName,svals[0]);
    exit(1);
  }
  for (i=1; i<nsv; i++){
    if ( svals[i] > svals[i-1] ) {
      printf("%s: ERROR: Singular values returned by clapack \n"
	     "  are not approprately ordered!\n"
	     "  svals[%d] = %16.6e > svals[%d] = %16.6e",
	     fctName,i,svals[i],i-1,svals[i-1]);
      exit(1);
    }
  }

  /* Test rank of matrix by examining the singular values.             */
  /* The singular values of the matrix is sorted in decending order, 
   * so that svals[i] >= svals[i+1]. The first that is zero (to some 
   * precision) yields information on the rank (to some precision)     */
  rank = nsv;
  tol  = DBL_EPSILON * svals[0] * MAX(rows,cols); 
  for (i=0; i<nsv; i++){
    if ( svals[i] <= tol ){
      rank = i;
      break;
    }
  }

  /* Compute (pseudo-) inverse matrix using the computed SVD:
   *   A   = U*S*V'
   *   A^+ = V * (Zi) U'
   *       = V * (U * (Zi)')'
   *       = { (U * Zi') * V' }'
   * Here Zi is the "inverse" of the diagonal matrix S, i.e. Zi is 
   * diagonal and has the same size as S'. The non-zero entries in 
   * Zi is calculated from the non-zero entries in S as Zi_ii=1/S_ii
   * Note that Zi' is of the same size as S.        
   * The last line here is used in the present computation. This 
   * notation avoids any need to transpose the output from the CLAPACK
   * rotines (which deliver U, S and V).                               */

  /* Inverse of [non-zero part of] diagonal matrix                     */
  tol = 1.0e-10 * svals[0]; 
  for (i=0; i< nsv; i++) {
    if (svals[i] < tol)
      svals[i] = 0.0; 
    else 
      svals[i] = 1.0 / svals[i]; 
  }

  /* Calculate  UZ = U * Zi', ie. scale COLUMN j in U by 
   * the j'th singular value [nsv columns only - since the diagonal 
   * matrix in general is not square]. If rows>cols then the last 
   * columns in UZ wil be zero (no need to compute).                   */
  for (j=0; j<nsv; j++){
    for (i=0; i<rows; i++){
      UU(i,j) *= svals[j];
    }
  }
  /* U*Zi' is stored in array U. It has size (rows x nsv). 
   * If cols>rows, then it should be though of as the larger matrix 
   * (rows x cols) with zero columns added on the right.               */

  /* Zero out the full array A to avoid confusion upon return.
   * This is not abosolutely necessary, zeroin out A[j][i] could be
   * part of the next loop.                                            */
  for (i=0; i< MAX(cols,rows); i++) {
    for (j=0; j< MAX(cols,rows); j++) {
      A[i][j] = 0.0;
    }
  }

  /* Matrix-matrix multiply  (U*Zi') * V'. 
   * Only the first nsv columns in U*Zi' are non-zero, so the inner 
   * most loop will go to k=nsv (-1). 
   * The result will be the transpose of the (psuedo-) inverse of A, 
   * so store directly in A'=A[j][i]. 
   *  A[j][i] = sum_k (U*Zi)[i][k] * V'[k][j] :                        */
  for (i=0; i< rows; i++) {
    for (j=0; j< cols; j++) {
      /*A[j][i] = 0.0; */ /* Include if A is not zeroed out above      */
      for (k=0;k<nsv;k++){
	A[j][i] += UU(i,k)*VV(k,j);
      }
    }
  }

  /* Free memory allocated in this routine */
  FREE(amat);
  FREE(svals);
  FREE(U);
  FREE(V);
  FREE(work);

  return rank;

  /* Undefine macros for this routine */
#undef AA
#undef UU
#undef VV
} /* End of routine pinv */



