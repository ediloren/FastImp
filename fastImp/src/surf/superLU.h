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

  const static char cvsid[] = "$Id: superLU.h,v 1.5 2002/12/31 21:30:39 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _SUPER_LU_H_
#define _SUPER_LU_H_

#include "zsp_defs.h" 
#include "util.h" 
#include "spColMat.h"

namespace surf {

  class SuperLU {

  public:
    SuperLU(void) : memoryAllocated(false) {}
    ~SuperLU(void);
    void allocate (void);
    void setup(pfft::SpColMat<std::complex<double> >& preCondMat);
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
    doublecomplex *mat;
    int *colIndex, *colStart;
    int *perm_c, *perm_r, *etree;
    doublecomplex *superLUrhs, *superLUsol;
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
      superLUrhs[i].r = rhs[i].real();
      superLUrhs[i].i = rhs[i].imag();
    }
    zCreate_Dense_Matrix(&X, numRow, 1, superLUrhs, numRow, DN, _Z, GE);
    
    int info;
    char trant[1];
    *trant = 'N';
    zgstrs (trant, &L, &U, perm_r, perm_c, &X, &info);
    if (info) {
      errorMessage("SuperLU::solve()",
		   "SuperLU fails to solve the sparse matrix equation");
    }

    DNformat *Xstore = (DNformat *) X.Store;
    doublecomplex* superLUsol = (doublecomplex *) Xstore->nzval;
    for (size_t i = 0; i < numRow; i++) {
      sol[i] = std::complex<double>(superLUsol[i].r, superLUsol[i].i);
    }

    Destroy_SuperMatrix_Store(&X);
  }

} // namespace surf

#endif
