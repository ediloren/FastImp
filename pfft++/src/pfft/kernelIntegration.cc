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
 
 ==========================================================================
*/

const static char cvsid[] = "$Id: kernelIntegration.cc,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

#include <stdexcept> //for exception handling
#include <algorithm>
#include "kernelIntegration.h"

using namespace pfft;

/**********************************************************************
 * KernelIntegration --
 **********************************************************************/
KernelIntegration::
KernelIntegration (
		   const std::vector<DifferentialOperator>& outerOperatorList,
		   const std::vector<DifferentialOperator>& innerOperatorList)
  : innerOperatorList_(innerOperatorList), outerOperatorList_(outerOperatorList) 
{ 
  if (innerOperatorList_.size() != outerOperatorList_.size()) {
    errorMessage("kernelIntegration.h: ",
		 "inner and outer operator lists are not of same length");
    throw std::domain_error("kernelIntegration.h: inner and outer operator lists are not of same length");
  }
  
  for (size_t i = 0; i < innerOperatorList_.size(); i++) {
    integralTypeList_.push_back(IntegralType(outerOperatorList_[i], 
					     innerOperatorList_[i]));
  }

  sort(innerOperatorList_.begin(), innerOperatorList_.end());
  innerOperatorList_.erase(unique(innerOperatorList_.begin(), 
				  innerOperatorList_.end()),
			   innerOperatorList_.end());
  sort(outerOperatorList_.begin(), outerOperatorList_.end());
  outerOperatorList_.erase(unique(outerOperatorList_.begin(), 
				  outerOperatorList_.end()),
			   outerOperatorList_.end());
}

/**********************************************************************
 * findInnerOperatorIndex --
 **********************************************************************/
const size_t 
KernelIntegration::
findInnerOperatorIndex (
			const DifferentialOperator innerOperator) const 
{
  std::vector<DifferentialOperator>::const_iterator 
    pos = find(innerOperatorList_.begin(), innerOperatorList_.end(), innerOperator);
  if (pos == innerOperatorList_.end()) {
    errorMessage("kernelIntegration.h: ",
		 "Can not find innerOperatorIndex, there must be a bug");
    throw std::domain_error("kernelIntegration.h: can not find innerOperatorIndex, there must be a bug");
  }

  return  pos - innerOperatorList_.begin();
}

/**********************************************************************
 * findOuterOperatorIndex --
 **********************************************************************/
const size_t 
KernelIntegration::
findOuterOperatorIndex (
			const DifferentialOperator outerOperator) const 
{
  std::vector<DifferentialOperator>::const_iterator 
    pos = find(outerOperatorList_.begin(), outerOperatorList_.end(), outerOperator);
  if (pos == outerOperatorList_.end()) {
    errorMessage("kernelIntegration.h: ",
		 "Can not find outerOperatorIndex, there must be a bug");
    throw std::domain_error("kernelIntegration.h: can not find outerOperatorIndex, there must be a bug");
  }

  return  pos - outerOperatorList_.begin();
}

/**********************************************************************
* findIntegralIndex --
**********************************************************************/
const size_t 
KernelIntegration::
findIntegralIndex (
		   const DifferentialOperator outerOperator, 
		   const DifferentialOperator innerOperator) const
{
  IntegralType integralType(outerOperator, innerOperator);
  std::vector<KernelIntegration::IntegralType>::const_iterator 
    pos = find(integralTypeList_.begin(), integralTypeList_.end(), integralType);
  if (pos == integralTypeList_.end()) {
    errorMessage("kernelIntegration.h: ",
		 "Can not find integralTypeIndex, there must be a bug");
    throw std::domain_error("kernelIntegration.h: can not find outerOperatorIndex, there must be a bug");
  }

  return  pos - integralTypeList_.begin();
}
