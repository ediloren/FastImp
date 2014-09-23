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

  This is the shared part of the panel integration function calcpForEikROverR
  and calcpForOneOverR.

  Resources:

  See also:

  const static char cvsid[] = "$Id: kernelIntegration.h,v 1.2 2003/05/15 17:59:50 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _KERNEL_INTEGRATION_H_
#define _KERNEL_INTEGRATION_H_

#include <vector>
#include <utility> // for pair
#include "utils.h" // for errorMessage

namespace pfft {

  enum IntegrationScheme { COLLOCATION, GALERKIN };
  enum DifferentialOperator {NONE, D_DX, D_DY, D_DZ, D_DN};  

  class KernelIntegration {

  public:
    
    typedef std::pair<DifferentialOperator, DifferentialOperator> IntegralType;

    KernelIntegration (const std::vector<DifferentialOperator>&,
		       const std::vector<DifferentialOperator>&);
    KernelIntegration (void) {}

    const size_t numIntegral(void) const { return integralTypeList_.size(); }
    const IntegralType integralType(size_t i) const { return integralTypeList_[i]; }
    const size_t findIntegralIndex(const DifferentialOperator, 
				   const DifferentialOperator) const;
    const DifferentialOperator outerOperatorOfIntegral(size_t i) { 
      return integralTypeList_[i].first; }
    const DifferentialOperator innerOperatorOfIntegral(size_t i) { 
      return integralTypeList_[i].second; }

    const size_t numInnerOperator(void) const { return innerOperatorList_.size(); }
    const size_t findInnerOperatorIndex(const DifferentialOperator innerOperator) const;
    const DifferentialOperator innerOperator(size_t index) const { 
      return innerOperatorList_[index]; }

    const size_t numOuterOperator(void) const { return outerOperatorList_.size(); }
    const size_t findOuterOperatorIndex(const DifferentialOperator outerOperator) const;
    const DifferentialOperator outerOperator(size_t index) const { 
      return outerOperatorList_[index]; }

  private:
    std::vector<DifferentialOperator> innerOperatorList_;
    std::vector<DifferentialOperator> outerOperatorList_;
    std::vector<IntegralType> integralTypeList_;

  };

} //namespace pfft

#endif
