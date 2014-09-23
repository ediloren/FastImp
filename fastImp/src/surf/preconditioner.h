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

  Resources:

  See also:

  const static char cvsid[] = "$Id: preconditioner.h,v 1.7 2002/08/11 15:06:42 zhzhu Exp $";

  ==========================================================================
*/


#ifndef _PRE_CONDITIONER_H_
#define _PRE_CONDITIONER_H_

#include <complex>
#include "surfConst.h"
#include "condInfo.h"
#include "mesh.h"
#include "element.h"
#include "calcpForEikrOverR.h"
#include "eikrOverR.h"
#include "fullwaveCollocation.h"
#include "spColMat.h"
#include "formulation.h"
#include "superLU.h"

namespace surf {

  class Preconditioner {

  public:
    Preconditioner(const Formulation& formulationIn, 
		   const PreconditionerType precondTypeIn,
		   const SurfConst&,
		   const bool useGlobalCoord);
    ~Preconditioner(void) {}
    void setup (const double freq);
    template<class VecX, class VecRHS> void solve (const VecRHS&, VecX&);
    const PreconditionerType type(void) const { return precondType; }
    const bool isLeftPreconditioner(void) const { 
      return precondType == LEFT_PRE_CONDITIONER; }
    const bool isRightPreconditioner(void) const { 
      return precondType == RIGHT_PRE_CONDITIONER; }

  private:
    const Formulation& fp;
    const PreconditionerType precondType;
    size_t totalNumE_;
    size_t totalNumRow_;
    pfft::SpColMat<std::complex<double> > preCondMat;
    SuperLU superLU;
    const SurfConst surfConst_;
    bool useGlobalCoord_;

    double K0;
    std::vector<std::complex<double> > KC;
    std::vector<pfft::DifferentialOperator> innerOperatorList1;
    std::vector<pfft::DifferentialOperator> outerOperatorList1;
    std::vector<pfft::DifferentialOperator> innerOperatorList2;
    std::vector<pfft::DifferentialOperator> outerOperatorList2;
    std::vector<Formulation::Calcp1> calcp1List;
    Formulation::Calcp2_FW calcp2_FW;
    Formulation::Calcp2_EMQS calcp2_EMQS;
    Formulation::Calcp3_FW calcp3_FW;
    Formulation::Calcp3_EMQS calcp3_EMQS;

    // These are for debugging purpose
    std::vector<int> rowBookKeeping;
    std::vector<int> colBookKeeping;
    std::vector<std::complex<double> > rowSum;
    std::vector<std::complex<double> > colSum;

    bool printPrecondSolveTime_;
    TNT::stopwatch timer;

    void copyNonIntegralEquaMat(void);
    void fillNonIntegralEquaMatInLocalCoord(void);
    void fillSelfTermOfIntegralEqua(void);
    void fillEqua1_old(void);
    void fillEqua1_new(void);
    void fillEqua2a(void);
    void fillEqua1_old_inLocalCoord(void);
    void fillEqua1_new_inLocalCoord(void);
    void fillEqua2a_inLocalCoord(void);
    void fillEqua1_old_inGlobalCoord(void);
    void fillEqua1_new_inGlobalCoord(void);
    void fillEqua2a_inGlobalCoord(void);
    void fillEqua3(void);
    void fillEqua4_freq_depend(double freq);
    void fillEqua6_inLocalCoord(const size_t, const size_t);
    void insertSpEle(const size_t, const size_t, const double);
    void insertSpEle(const size_t, const size_t, const std::complex<double>&);
    void updateRowBookKeeping(size_t row);
    void updateColBookKeeping(size_t col);
    void updateRowSum (size_t row, const std::complex<double>& value);
    void updateRowSum (size_t row, const double value);
    void updateColSum (size_t col, const std::complex<double>& value);
    void updateColSum (size_t col, const double value);
    void checkColBookKeeping(void);
    void checkRowBookKeeping(void);
    void checkRowSum (void);
    void checkColSum (void);
    void checkIndex(size_t index);

  };

  template<class VecX, class VecRHS> 
  void 
  Preconditioner::solve (
			 const VecRHS& rhs,
			 VecX& sol)
  {
    if (printPrecondSolveTime_) timer.start();

    if (precondType != NO_PRE_CONDITIONER) {
      superLU.solve(rhs, sol);      
    } else {
      if (&sol != &rhs) {
	for (size_t i = 0; i < totalNumRow_; i++) {
	  sol[i] = rhs[i];
	}
      }	
    }

    if (printPrecondSolveTime_) {
      timeReport("Time used for one pre-conditioner solve := ", timer.read());
      printPrecondSolveTime_ = false;
    }
  }

  
} // namespace surf

#endif
