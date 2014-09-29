/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Zhenhai Zhu

  Description:

  interface for superLU c code

  Resources:

  See also:

  const static char cvsid[] = "$Id: superLU.cc,v 1.11 2002/12/31 21:30:39 zhzhu Exp $";

  ==========================================================================
*/

#include "superLU.h"
#include "service.h"
// Enrico
#include "zsp_defs.h"
#include "util.h"

using namespace surf;
using namespace std;
using namespace pfft;

/**********************************************************************
 * setup --
 * It apears that superLU only support LU of compressed column format.
 **********************************************************************/
void
SuperLU::setup (
		pfft::SpColMat<std::complex<double> >& systemMat)
{
  numRow = systemMat.numRow();
  numCol = systemMat.numCol();
  numNonZero = systemMat.numNonZero();

  if (!memoryAllocated) allocate();

  int elementCount = 0;
  for (size_t col = 0; col < numCol; col++) {
    colStart[col] = elementCount;
    for (size_t row = 0; row < systemMat[col].size(); row++) {
      mat[elementCount].r = systemMat[col].value(row).real();
      mat[elementCount].i = systemMat[col].value(row).imag();
      colIndex[elementCount] = systemMat[col].index(row);
      elementCount++;
    }
  }
  colStart[systemMat.numCol()] = elementCount;

  // Structure A has its internal pointers pointing to mat, colIndex and
  // colStart, so there is no redundant memory allocation.
  // This also means that I should only deallocate A.
  // Otherwie I will free dandling pointer. This will cause segmentation fault.
  zCreate_CompCol_Matrix(&A, numRow, numCol, numNonZero, mat, colIndex,
			 colStart, NC, _Z, GE);
}

/**********************************************************************
 * allocate --
 **********************************************************************/
void
SuperLU::allocate (
		   void)
{
  mat = (doublecomplex *) malloc(numNonZero * sizeof(doublecomplex));
  colIndex = (int *) malloc(numNonZero * sizeof(int));
  colStart = (int *) malloc((numCol+1) * sizeof(int));
  superLUrhs = (doublecomplex*)malloc(numRow * sizeof(doublecomplex));
  perm_r = (int*) malloc(numRow * sizeof(int));
  perm_c = (int*) malloc(numCol * sizeof(int));
  etree = (int*) malloc(numCol * sizeof(int));
  memoryAllocated = true;
}

/**********************************************************************
 * factor --
 **********************************************************************/
void
SuperLU::factor (
		 void)
{
  char refact[1];
  *refact = 'N';
  int relax      = sp_ienv(2);
  int panel_size = sp_ienv(1);
  StatInit(panel_size, relax);

  /*
   * GET_PERM_C obtains a permutation matrix Pc, by applying the multiple
   * minimum degree ordering code by Joseph Liu to matrix A'*A or A+A'.
   * or using approximate minimum degree column ordering by Davis et. al.
   * The LU factorization of A*Pc tends to have less fill than the LU
   * factorization of A.
   * permc_spec specifies the type of column ordering to reduce fill:
   *         = 1: minimum degree on the structure of A^T * A
   *         = 2: minimum degree on the structure of A^T + A
   *         = 3: approximate minimum degree for unsymmetric matrices
   *         If ispec == 0, the natural ordering (i.e., Pc = I) is returned.
   */
  // this option turns out to be much better for gspiral4.inp and wire with
  // triangular panels. So I assume it could be a better option for
  // other structures as well.
  int permc_spec = 2;
  get_perm_c(permc_spec, &A, perm_c);

  SuperMatrix AC;
  sp_preorder(refact, &A, perm_c, etree, &AC);

  int info;
  int lwork=0;
  void *work;
  double diag_pivot_thresh = 1.0; // this corresponding to partial pivoting
  //  double diag_pivot_thresh = 0.0; // no pivoting, turns out to be useless
  double drop_tol = 0.;

  zgstrf (refact, &AC, diag_pivot_thresh, drop_tol, relax, panel_size,
	  etree, work, lwork, perm_r, perm_c, &L, &U, &info);
  if (info) {

// enrico debug
fprintf(stderr, "error #%i, AC.ncol %i\n", info, AC.ncol);

    errorMessage("SuperLU::factor()",
		 "SuperLU fails to factor the sparse matrix");

  }
  Destroy_CompCol_Permuted(&AC);

  progressReport("\tThe sparse matrix has been LU factored");
  messageWithOneNumber("The number of nonzero's before LU :=", numNonZero);

  if ((L.nrow != L.ncol || L.Mtype != TRLU) ||
      (U.nrow != U.ncol || U.Mtype != TRU) ) {
    errorMessage("SuperLU::factor()",
		 "Run out of memory! LU of the sparse matrix terminated prematurely");
  } else {
    int numNonZeroInLU = ((SCformat *)L.Store)->nnz + ((NCformat *)U.Store)->nnz;
    int numFillin = numNonZeroInLU - numNonZero;
    messageWithOneNumber("The number of nonzero's after LU :=" , numNonZeroInLU);
    messageWithOneNumber("The number of fill-in's := ", numFillin);
  }
}

/**********************************************************************
 * estimateCondNum_OneNorm --
 **********************************************************************/
double
SuperLU::estimateCondNum_OneNorm (
				  void)
{
  char norm[1];
  *norm = '1';
  double anorm = 1;
  double rcond;
  int info;
  zgscon(norm, &L, &U, anorm, &rcond, &info);
  if (info) {
    errorMessage("SuperLU::estimateCondNum_OneNorm()",
		 "SuperLU fails to extimate the condition number of the sparse matrix");
  }

  return 1./rcond;
}

/**********************************************************************
 * estimateCondNum_InfNorm --
 **********************************************************************/
double
SuperLU::estimateCondNum_InfNorm (
				  void)
{
  char norm[1];
  *norm = 'I';
  double anorm = 2;
  double rcond;
  int info;
  zgscon(norm, &L, &U, anorm, &rcond, &info);

  return 1./rcond;
}

/**********************************************************************
 * ~SuperLU --
 **********************************************************************/
SuperLU::~SuperLU (
		   void)
{
  if (memoryAllocated) {
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    //    Destroy_Dense_Matrix(&X);
    // this is redundant because its store has been freed in solve(),
    // and the nzval will be freed since it is the same as superLUrhs
    free(superLUrhs);
    free(perm_c);
    free(perm_r);
    free(etree);
  }
}


