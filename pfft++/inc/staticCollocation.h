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

  This is the interface between pfft and the c version fullwave collocation 
  panel integration function, surfCalcp.

  Resources:

  See also:

  const static char cvsid[] = "$Id: staticCollocation.h,v 1.9 2002/10/01 21:17:37 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _STATIC_COLLOCATION_H_
#define _STATIC_COLLOCATION_H_

#include <vector>
#include <algorithm>
#include <stdexcept> //for exception handling
#include <utility> // for pair
#include "vector3D.h"
#include "kernelIntegration.h"
#include "calcpForOneOverR.h"

namespace pfft {

  template <class Panel>

  class StaticCollocation : public KernelIntegration {

  public:
    
    StaticCollocation (const std::vector<DifferentialOperator>&,
		       const std::vector<DifferentialOperator>&,
		       const bool skipPanelIntegration = false);

    StaticCollocation (void) {};



    const IntegrationScheme itegrationScheme(void) const { return COLLOCATION; }
    void operator () (const Panel& srcPanel, const Panel& evalPanel);
    double result(const size_t integralIndex) const;

  private:
    
    CalcpForOneOverR<Panel> calcp;

    enum ResultType { SINGLE_LAYER, D_DN_SINGLE_LAYER, 
		      D_DX_SINGLE_LAYER, D_DY_SINGLE_LAYER, D_DZ_SINGLE_LAYER, 
		      GRADIENT_SINGLE_LAYER, GRADIENT_DOUBLE_LAYER,
		      DOUBLE_LAYER, D_DN_DOUBLE_LAYER};

    std::vector<ResultType> resultTypeList_;
    std::vector<ResultType> calcpResultTypeList_;

    double slp;
    double dlp;
    double dlp_dn;
    double slp_dn;
    vector3D<double> grad_slp;
    vector3D<double> grad_dlp;

    bool skipPanelIntegration_;

    void setupResultTypeList(void);
    ResultType mapDiffOperatorToResultType (const IntegralType integralType) const;
    
    bool IsWanted(ResultType rt);
    bool IsWantedFromCalcp(ResultType rt) ;
  };

  /**********************************************************************
   * StaticCollocation --
   **********************************************************************/
  template <class Panel>
  StaticCollocation<Panel>::
  StaticCollocation (
		     const std::vector<DifferentialOperator>& outerOperatorList,
		     const std::vector<DifferentialOperator>& innerOperatorList,
		     const bool skipPanelIntegration)
    : KernelIntegration(outerOperatorList, innerOperatorList), 
		       skipPanelIntegration_(skipPanelIntegration) 
  { 
    setupResultTypeList(); 

    if (! skipPanelIntegration_) {
      bool slpWanted = IsWantedFromCalcp(SINGLE_LAYER);
      bool dlpWanted = IsWantedFromCalcp(DOUBLE_LAYER);
      bool gradOfSlpWanted = IsWantedFromCalcp(GRADIENT_SINGLE_LAYER);
      bool gradOfDlpWanted = IsWantedFromCalcp(GRADIENT_DOUBLE_LAYER);
      calcp = CalcpForOneOverR<Panel>(slpWanted, dlpWanted, gradOfSlpWanted, 
				      gradOfDlpWanted);
    }
  }

  /**********************************************************************
   * result --
   **********************************************************************/
  template <class Panel>
  double 
  StaticCollocation<Panel>::
  result (
	  const size_t integralIndex) const
  {
    switch (resultTypeList_[integralIndex]) {
    case SINGLE_LAYER:
      return slp;
      break;
    case D_DN_SINGLE_LAYER:
      return slp_dn;
      break;
    case D_DX_SINGLE_LAYER:
      return grad_slp.x();
      break;
    case D_DY_SINGLE_LAYER:
      return grad_slp.y();
      break;
    case D_DZ_SINGLE_LAYER:
      return grad_slp.z();
      break;
    case DOUBLE_LAYER:
      return dlp;
      break;
    case D_DN_DOUBLE_LAYER:
      return dlp_dn;
      break;
    default:
      throw std::domain_error("staticCollocation.h: unknown resultType");
      break;
    }
  }

  /**********************************************************************
   * setupResultTypeList --
   **********************************************************************/
  template <class Panel>
  void
  StaticCollocation<Panel>::setupResultTypeList (
						 void)
  {
    for (size_t i = 0; i < numIntegral(); i++) {
      resultTypeList_.push_back(mapDiffOperatorToResultType(integralType(i)));
    }

    calcpResultTypeList_ = resultTypeList_;
    for (size_t i = 0; i < numIntegral(); i++) {
      switch (calcpResultTypeList_[i]) {
      case D_DN_SINGLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_SINGLE_LAYER;
	break;
      case D_DX_SINGLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_SINGLE_LAYER;
	break;
      case D_DY_SINGLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_SINGLE_LAYER;
	break;
      case D_DZ_SINGLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_SINGLE_LAYER;
	break;
      case D_DN_DOUBLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_DOUBLE_LAYER;
	break;
      }    
    }
    sort(calcpResultTypeList_.begin(), calcpResultTypeList_.end());
    calcpResultTypeList_.erase(unique(calcpResultTypeList_.begin(), calcpResultTypeList_.end()),
				calcpResultTypeList_.end());
  }

  /**********************************************************************
   * mapDiffOperatorToResultType --
   **********************************************************************/
  template <class Panel>
  typename StaticCollocation<Panel>::ResultType   
  StaticCollocation<Panel>::
  mapDiffOperatorToResultType (
			       const IntegralType integralType) const
  {
    ResultType resultType;

    switch (integralType.first) {
    case NONE:
      switch(integralType.second) {
      case NONE:
	resultType = SINGLE_LAYER;	
	break;
      case D_DN:
	resultType = DOUBLE_LAYER;	
	break;
      case D_DX:
      case D_DY:
      case D_DZ:
	throw std::domain_error("staticCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("staticCollocation.h: unknown inner deferential operator");
      }
      break;
    case D_DN:
      switch(integralType.second) {
      case NONE:
	resultType = D_DN_SINGLE_LAYER;	
	break;
      case D_DN:
	resultType = D_DN_DOUBLE_LAYER;	
	break;
      case D_DX:
      case D_DY:
      case D_DZ:
	//	resultType = HYPER_SINGULAR;	
	throw std::domain_error("staticCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("staticCollocation.h: unknown inner deferential operator");
      }
      break;
    case D_DX:
      switch(integralType.second) {
      case NONE:
	resultType = D_DX_SINGLE_LAYER;	
	break;
      case D_DN:
      case D_DX:
      case D_DY:
      case D_DZ:
	throw std::domain_error("staticCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("staticCollocation.h: unknown inner deferential operator");
      }
      break;
    case D_DY:
      switch(integralType.second) {
      case NONE:
	resultType = D_DY_SINGLE_LAYER;	
	break;
      case D_DN:
      case D_DX:
      case D_DY:
      case D_DZ:
	throw std::domain_error("staticCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
      }
      break;
    case D_DZ:
      switch(integralType.second) {
      case NONE:
	resultType = D_DZ_SINGLE_LAYER;	
	break;
      case D_DN:
      case D_DX:
      case D_DY:
      case D_DZ:
	throw std::domain_error("staticCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("staticCollocation.h: unknown inner deferential operator");
      }
      break;
    default:
      throw std::domain_error("staticCollocation.h: unknown OUTER deferential operator");
      break;
    }
    return resultType;
  }

  /**********************************************************************
   * operator () --
   **********************************************************************/
  template <class Panel>
  void 
  StaticCollocation<Panel>::operator () (
					 const Panel& srcPanel, 
					 const Panel& evalPanel)
  {
    if (skipPanelIntegration_) {
      // this is for debugging purpose. So I do not have to wait for the real calculation
      slp = 1.;
      dlp = 1;
      grad_slp = pfft::vector3D<double>(slp, slp, slp);
      grad_dlp = pfft::vector3D<double>(slp, slp, slp);
      dlp_dn = 1.;
      slp_dn = 1.;

    } else {
      calcp(srcPanel, evalPanel.centroid(), slp, dlp, grad_slp, grad_dlp);
      if (IsWanted(D_DN_SINGLE_LAYER)) {
	dotProd(slp_dn, evalPanel.normal(), grad_slp);
      }     
      if (IsWanted(D_DN_DOUBLE_LAYER)) {
	dotProd(dlp_dn, evalPanel.normal(), grad_dlp);
      }
    }
  }

  /**********************************************************************
   * IsWanted --
   **********************************************************************/
  template <class Panel>
  bool
  StaticCollocation<Panel>::IsWanted (
				      ResultType rt) 
  {
    return (find(resultTypeList_.begin(), resultTypeList_.end(), rt) 
	    != resultTypeList_.end());
  }


  /**********************************************************************
   * need --
   **********************************************************************/
  template <class Panel>
  bool
  StaticCollocation<Panel>::IsWantedFromCalcp (
					       ResultType rt) 
  {
    return (find(calcpResultTypeList_.begin(), calcpResultTypeList_.end(), 
		 rt) != calcpResultTypeList_.end());
  }


} //namespace pfft

#endif

