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

  interface for superLU c code, unknown is double

  Resources:

  See also:

  const static char cvsid[] = "$Id: realsuperLU.h,v 1.1 2003/02/11 02:47:36 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _REAL_SUPER_LU_H_
#define _REAL_SUPER_LU_H_

#include "dsp_defs.h" // interface to superLU
#include "util.h" // interface to superLU
#include "utils.h" // this is for pfft::errorMessage
#include "spColMat.h"

namespace superLU {

  class SuperLU {

  public:
    SuperLU(void) : memoryAllocated(false) {}
    ~SuperLU(void);
    void allocate (void);
    void setup(pfft::SpColMat<double>& preCondMat);
    void factor (void);
    template<class VecX, class VecRHS> void solve (const VecRHS&, VecX&);
    double estimateCondNum_OneNorm(void);
    double estimateCondNum_InfNorm(void);

  private:
    bool memoryAllocated;
    int numRow;
    int numCol;
    int numNonZero;
    SuperMatrix A, L, U, X;
    double *mat;
    int *colIndex, *colStart;
    int *perm_c, *perm_r, *etree;
    double *superLUrhs, *superLUsol;
  };

  /**********************************************************************
   * solve --
   **********************************************************************/
  template<class VecX, class VecRHS> 
  void
  SuperLU::solve (
		  const VecRHS& rhs,
		  VecX& sol)
  {
    for (size_t i = 0; i < numRow; i++) {
      superLUrhs[i] = rhs[i];
    }
    dCreate_Dense_Matrix(&X, numRow, 1, superLUrhs, numRow, DN, _D, GE);
    
    int info;
    char trant[1];
    *trant = 'N';
    dgstrs (trant, &L, &U, perm_r, perm_c, &X, &info);
    if (info) {
      pfft::errorMessage("SuperLU::solve()",
			 "SuperLU fails to solve the sparse matrix equation");
    }

    DNformat *Xstore = (DNformat *) X.Store;
    double* superLUsol = (double *) Xstore->nzval;
    for (size_t i = 0; i < numRow; i++) {
      sol[i] = superLUsol[i];
    }

    Destroy_SuperMatrix_Store(&X);
  }

} // namespace superLU

#endif
