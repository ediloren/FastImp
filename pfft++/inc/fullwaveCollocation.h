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

  This is the interface for the fullwave collocation panel integration function.

  Resources:

  See also:

  const static char cvsid[] = "$Id: fullwaveCollocation.h,v 1.6 2003/01/08 16:08:37 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _FULLWAVE_COLLOCATION_H_
#define _FULLWAVE_COLLOCATION_H_

#include <vector>
#include <algorithm>
#include <stdexcept> //for exception handling
#include <utility> // for pair
#include "vector3D.h"
#include "kernelIntegration.h"
#include "calcpForEikrOverR.h"

namespace pfft {

  template <class KT, class Panel>

  class FullwaveCollocation : public KernelIntegration {

  public:
    
    FullwaveCollocation (void) {};
    FullwaveCollocation (const std::vector<DifferentialOperator>&,
			 const std::vector<DifferentialOperator>&,
			 const KT KIn,
			 const int numGaussQuadraturePoint = 24,
			 const CalcpEvalPntDirection evalPointDirection = 
			 FOLLOW_NORMAL_DIRECTION);

    const IntegrationScheme itegrationScheme(void) const { return COLLOCATION; }
    void operator () (const Panel& srcPanel, const Panel& evalPanel);
    std::complex<double> result(const size_t integralIndex) const;

  private:

    CalcpForEikrOverR<KT, Panel> calcp;
    KT K_;
    enum ResultType { SINGLE_LAYER, D_DN_SINGLE_LAYER, 
		      D_DX_SINGLE_LAYER, D_DY_SINGLE_LAYER, D_DZ_SINGLE_LAYER, 
		      GRADIENT_SINGLE_LAYER, 
		      DOUBLE_LAYER, D_DN_DOUBLE_LAYER};
    std::vector<ResultType> resultTypeList_;
    std::vector<ResultType> calcpResultTypeList_;

    std::complex<double> slp;
    std::complex<double> dlp;
    std::complex<double> dlp_dn;
    std::complex<double> slp_dn;
    pfft::vector3D<std::complex<double> > grad_slp;

    bool skipPanelIntegration_;

    void setupResultTypeList(void);
    ResultType mapDiffOperatorToResultType (const IntegralType integralType) const;
    std::complex<double> compHyperSingularIntegral (const Panel& srcPanel, 
						    const Panel& evalPanel) const;
    bool IsWanted(ResultType rt) ;
    bool IsWantedFromCalcp(ResultType rt) ;
  };

  /**********************************************************************
   * FullwaveCollocation --
   **********************************************************************/
  template <class KT, class Panel>
  FullwaveCollocation<KT, Panel>::
  FullwaveCollocation (
		       const std::vector<DifferentialOperator>& outerOperatorList,
		       const std::vector<DifferentialOperator>& innerOperatorList,
		       const KT KIn,
		       const int numGaussQuadraturePoint,
		       const CalcpEvalPntDirection evalPointDirection)
    : KernelIntegration(outerOperatorList, innerOperatorList), K_(KIn) 
  { 
    setupResultTypeList(); 

    if (numGaussQuadraturePoint == 0) {
      // I use this as a sign to skip panel integration. This is for debugging purpose.
      skipPanelIntegration_ = true;      

    } else {
      bool slpWanted = IsWantedFromCalcp(SINGLE_LAYER);
      bool dlpWanted = IsWantedFromCalcp(DOUBLE_LAYER);
      bool gradWanted = IsWantedFromCalcp(GRADIENT_SINGLE_LAYER);
      calcp = CalcpForEikrOverR<KT, Panel>(K_, evalPointDirection, slpWanted, 
					   dlpWanted, gradWanted, 
					   numGaussQuadraturePoint);
      skipPanelIntegration_ = false;
    }
  }

  /**********************************************************************
   * result --
   **********************************************************************/
  template <class KT, class Panel>
  std::complex<double> 
  FullwaveCollocation<KT, Panel>::result (
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
      throw std::domain_error("fullwaveCollocation.h: unknown resultType");
      break;
    }
  }

  /**********************************************************************
   * setupResultTypeList --
   **********************************************************************/
  template <class KT, class Panel>
  void
  FullwaveCollocation<KT, Panel>::setupResultTypeList (
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
      }    
    }
    sort(calcpResultTypeList_.begin(), calcpResultTypeList_.end());
    calcpResultTypeList_.erase(unique(calcpResultTypeList_.begin(), 
				      calcpResultTypeList_.end()),
			       calcpResultTypeList_.end());
  }

  /**********************************************************************
   * mapDiffOperatorToResultType --
   **********************************************************************/
  template <class KT, class Panel>
  typename FullwaveCollocation<KT, Panel>::ResultType   
  FullwaveCollocation<KT, Panel>::
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
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
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
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
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
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
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
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
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
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
      }
      break;
    default:
      throw std::domain_error("fullwaveCollocation.h: unknown OUTER deferential operator");
      break;
    }
    return resultType;
  }

  /**********************************************************************
   * operator () --
   **********************************************************************/
  template <class KT, class Panel>
  void 
  FullwaveCollocation<KT, Panel>::operator () (
					       const Panel& srcPanel, 
					       const Panel& evalPanel)
  {
    if (skipPanelIntegration_) {
      // this is for debugging purpose. So I do not have to wait for the real calculation
      slp = 1.;
      dlp = 1;
      grad_slp = pfft::vector3D<std::complex<double> >(slp, slp, slp);
      dlp_dn = 1.;
      slp_dn = 1.;

    } else {
      calcp(srcPanel, evalPanel.centroid(), slp, dlp, grad_slp);
      
      if (IsWanted(D_DN_SINGLE_LAYER)) {
	dotProd(slp_dn, evalPanel.normal(), grad_slp);
      }
      
      if (IsWanted(D_DN_DOUBLE_LAYER)) {
	dlp_dn = compHyperSingularIntegral(srcPanel, evalPanel);
      }
    }
  }

  /**********************************************************************
   * compHyperSingularIntegral --
   * use FD to approximate d_dn
   **********************************************************************/
  template <class KT, class Panel>
  std::complex<double> 
  FullwaveCollocation<KT, Panel>::
  compHyperSingularIntegral (
			     const Panel& srcPanel, 
			     const Panel& evalPanel) const
  {
    const double QLP_MUTUAL_FD_STEP = 0.02;
    const double QLP_SELF_FD_STEP = 0.01;

    double minEdgeLength = srcPanel.findMinEdgeLength();
    vector3D<double> evalPnt1, evalPnt2;
    double fdStep;
    if (srcPanel == evalPanel) {
      /* Since there is a jump caused by singularity here. I could only
	 choose evaluation point from the same side.
	 The evaluation point should be on the side from which the panel 
	 normal vector is leaving. */
      fdStep = QLP_SELF_FD_STEP * minEdgeLength;
      evalPnt1 = evalPanel.centroid() - 2*fdStep * evalPanel.normal();
      evalPnt2 = evalPanel.centroid() - fdStep * evalPanel.normal();
    } else {
      fdStep = QLP_MUTUAL_FD_STEP * minEdgeLength;
      evalPnt1 = evalPanel.centroid() - 0.5*fdStep * evalPanel.normal();
      evalPnt2 = evalPanel.centroid() + 0.5*fdStep * evalPanel.normal();
    }

    bool slpWanted = false;
    bool dlpWanted = true;
    bool gradWanted = false;
    CalcpEvalPntDirection evalPointDirection = FOLLOW_NORMAL_DIRECTION;
    CalcpForEikrOverR<KT, Panel> calcp(K_, evalPointDirection, slpWanted, 
				       dlpWanted, gradWanted);

    std::complex<double> fs, dfdn1, dfdn2;
    vector3D<std::complex<double> > grad_f;
    calcp(srcPanel, evalPnt1, fs, dfdn1, grad_f);
    calcp(srcPanel, evalPnt2, fs, dfdn2, grad_f);

    return (dfdn2 - dfdn1) / fdStep;
  }

  /**********************************************************************
   * IsWanted --
   **********************************************************************/
  template <class KT, class Panel>
  bool
  FullwaveCollocation<KT, Panel>::IsWanted (
					    ResultType rt) 
  {
    return (find(resultTypeList_.begin(), resultTypeList_.end(), rt) 
	    != resultTypeList_.end());
  }

  /**********************************************************************
   * IsWantedFromCalcp --
   **********************************************************************/
  template <class KT, class Panel>
  bool
  FullwaveCollocation<KT, Panel>::IsWantedFromCalcp (
						     ResultType rt) 
  {
    return (find(calcpResultTypeList_.begin(), calcpResultTypeList_.end(), rt) 
	    != calcpResultTypeList_.end());
  }

} //namespace pfft

#endif
