/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Zhenhai Zhu and Ben Song
  
  Description:

  This class computes single and double layer potential, 
  and the gradient of single layer potential
  
  slp = Int {exp(ikr)/r} 
  dlp = Int {d/dn' exp(ikr)/r}
  grad = grad { slp }
  
  The code also works if k=0, which means it could also handle 1/r kernel.

  The code could handle singularity case when the evaluation point is on the
  source panel. So you don't need to worry about the jump in the 
  double-layer potential.

  Resources:
 
  See also:  

  1) Hess and Smith's paper on decomposing integral over a polygon into 
  integrals over a set of triangles. This idea is used in the main function.
  2) Our paper published at ICCAD 2001 on piece-wise gauss quadrature scheme. 
  This idea is used in function compIntegralOnTriangle.

  const static char cvsid[] = "$Id: calcpForEikrOverR.h,v 1.9 2003/07/15 19:58:45 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _CALCP_FOR_EIKR_OVER_R_H_
#define _CALCP_FOR_EIKR_OVER_R_H_

#include <stdexcept>
#include <iostream>
#include <cfloat> // for DBL_EPSILON
#include "vector3D.h"
#include "element.h"
#include "utils.h"

namespace pfft{

  enum CalcpEvalPntDirection {
    AGAINST_NORMAL_DIRECTION = +1,
    FOLLOW_NORMAL_DIRECTION = -1
  };
  
  template<class KT, class Panel = element>
  class CalcpForEikrOverR{

    typedef vector3D<std::complex<double> > ComplexVector;
    
    static const size_t DEFAULT_NUM_GQ_POINT_ = 12;

  public:
    CalcpForEikrOverR() {};
    CalcpForEikrOverR(
		      const KT& K,
		      const CalcpEvalPntDirection direction,
		      const bool needSlp,
		      const bool needDlp,
		      const bool needGradient,
		      const size_t n = DEFAULT_NUM_GQ_POINT_);
    void operator() (Panel, vector3D<double>, std::complex<double>&, 
		     std::complex<double>&, vector3D<std::complex<double> >&);

  private:
    KT K_;
    CalcpEvalPntDirection evalPointDirection_; 
    CalcpEvalPntDirection evalPointDirectionBackup_; 
    size_t NumGQpoint_;
    std::vector<double> GQweight_;
    std::vector<double> GQpoint_;

    bool needSingleLayerPotential_;
    bool needDoubleLayerPotential_;
    bool needGradient_;
    bool needNothing_;
    bool projectPointOnEdge_;

    bool disablePiecewiseFlag;

    struct TriangleInfo_{
      vector3D<double> PH_unit;
      vector3D<double> AH;
      vector3D<double> HB;
      vector3D<double> X;
      vector3D<double> Y;
      double PH_len, AB_len, AH_len, area;
    };
    TriangleInfo_ triangleInfo_;

    void setupPanelLocalCoordSystem(Panel&, vector3D<double>&);
    void getEvalPointDirection(const vector3D<double>&);
    void setupCalcpTriangleInfo(size_t, const Panel&);
    double getContributionSign(void);

    void compIntegralOnTriangle (const vector3D<double>&, std::complex<double>&, 
				 ComplexVector&);
    void compOuterIntegral(const vector3D<double>&, double, double,
			   std::complex<double>&, std::complex<double>&, 
			   std::complex<double>&); 
    std::complex<double> compSingleLayerInnerIntegral(const std::complex<double>, 
						      const std::complex<double>,
						      const std::complex<double>, 
						      const double);
    std::complex<double> compDIDXInnerIntegral(const std::complex<double>, 
					       const std::complex<double>, 
					       const std::complex<double>,
					       const double, const double, const double, 
					       const double, const double);
    std::complex<double> compDIDYOnTriangle(const vector3D<double>&);
    std::complex<double> compDIDZInnerIntegral(const std::complex<double>, 
					       const std::complex<double>, const double, 
					       const double, const double, const double);

    ComplexVector transferToPanelLocalCoordSystem(const ComplexVector& grad1);
    void localDoubleLayerCorrection(std::complex<double>&);
    void flipSign(std::complex<double>&, ComplexVector&);

    void gaussQuadrature(const size_t, 
			 std::vector<double>&, std::vector<double>&);
  };

  /**********************************************************************
   * CalcpForEikrOverR --
   **********************************************************************/
  template<class KT, class Panel>
  CalcpForEikrOverR<KT, Panel>
  ::CalcpForEikrOverR (
		       const KT& K,
		       const CalcpEvalPntDirection direction,
		       const bool needSlp,
		       const bool needDlp,
		       const bool needGradient,
		       const size_t n)
    : K_(K), evalPointDirection_(direction), 
    evalPointDirectionBackup_(direction),
    needSingleLayerPotential_(needSlp),
    needDoubleLayerPotential_(needDlp), needGradient_(needGradient),
    projectPointOnEdge_(false), disablePiecewiseFlag(false)
  {    
    needNothing_ = (!needSlp) && (!needDlp) && (!needGradient);
    if (! needNothing_) {
      if (GQpoint_.size() != n) {
	// To avoid one of Gauss Quad point to be zero, please
	// use the even number of points;
	NumGQpoint_ = n;
	if(NumGQpoint_ % 2 != 0) {
	  NumGQpoint_ += 1;
	} 
	GQpoint_.resize(NumGQpoint_);
	GQweight_.resize(NumGQpoint_);
	gaussQuadrature(NumGQpoint_, GQweight_, GQpoint_);     
      }
#ifdef DISABLE_PIECEWISE_SCHEME
      if (!disablePiecewiseFlag) {
	// print out this message only once
	disablePiecewiseFlag = true;
	std::cout << "\n\n\t calcpForEikrOverR:" << endl
		  << "\t piecewise scheme has been disabled\n" << std::endl;
      }	
#endif
    }
  }

  /**********************************************************************
   * operator() --
   **********************************************************************/
  template<class KT, class Panel>
  void
  CalcpForEikrOverR<KT, Panel>::
  operator() (
	      Panel panel,
	      vector3D<double> evalPnt,
	      std::complex<double>& slp,
	      std::complex<double>& dlp,
	      vector3D<std::complex<double> >& grad)
  {
    if (needNothing_) return;

    setupPanelLocalCoordSystem(panel, evalPnt);
    if( needDoubleLayerPotential_ || needGradient_ ) {
      getEvalPointDirection(evalPnt);
    }
    
    slp = std::complex<double>(0., 0.);
    ComplexVector gradSumInPanelLocalCoordSys(ZERO_COMPLEX_VECTOR_3D);
    double totalArea = 0.;   
    for(size_t ii=0; ii<panel.shape(); ++ii){
      setupCalcpTriangleInfo(ii, panel);
      double contributionSign = getContributionSign();      
      totalArea += triangleInfo_.area * contributionSign;

      std::complex<double> slpLocal;
      ComplexVector gradInLocalTriangleCoordSys;
      compIntegralOnTriangle(evalPnt, slpLocal, gradInLocalTriangleCoordSys);

      if ( needSingleLayerPotential_ ) {
	slp += slpLocal * contributionSign;
      }

      if ( needDoubleLayerPotential_ || needGradient_ ) {
	ComplexVector gradInPanelLocalCoordSys;
	gradInPanelLocalCoordSys
	  = transferToPanelLocalCoordSystem(gradInLocalTriangleCoordSys);
	gradInPanelLocalCoordSys *= contributionSign;
	gradSumInPanelLocalCoordSys += gradInPanelLocalCoordSys;
      }
    }

    if (totalArea < 0) {
      flipSign(slp, gradSumInPanelLocalCoordSys);
    }
 
    if ( needDoubleLayerPotential_ || needGradient_ ) {
      if (evalPnt.z() == 0.)
	localDoubleLayerCorrection(gradSumInPanelLocalCoordSys.z());
    }

    // it is d/dn Int in grad. and Int d/dn'= - d/dn Int, hence the minus sign.
    dlp = - gradSumInPanelLocalCoordSys.z();

    if ( needGradient_ ) {
      // move back to global coordinate system
      grad = rotateCoord(gradSumInPanelLocalCoordSys, panel.tangent1(), 
			 panel.tangent2(), panel.normal());
    }

    evalPointDirection_ = evalPointDirectionBackup_;
  }

  /**********************************************************************
   * setupPanelLocalCoordSystem --
   **********************************************************************/
  template<class KT, class Panel>
  void
  CalcpForEikrOverR<KT, Panel>
  ::setupPanelLocalCoordSystem (
				Panel& panel,
				point3D& evalPnt)
  {
    // 1) find the projection of evaluation point on the panel
    point3D evalPntProjection(evalPnt);
    evalPntProjection.transferGlobalToLocalCoord(panel.centroid(),
						 panel.tangent1(), 
						 panel.tangent2(),
						 panel.normal());
    double h = evalPntProjection.z();
    evalPntProjection = point3D(evalPntProjection.x(), evalPntProjection.y(), 0.);

    // 2) shift back to the global coord system to be 
    // compatible with all verteices which are in the globl coord.
    evalPntProjection.transferLocalToGlobalCoord(panel.centroid(), 
						 panel.tangent1(), 
						 panel.tangent2(), 
						 panel.normal());

    // 3) use projection point as origin of the panel local coordinate system
    panel.shiftToLocalCoord(evalPntProjection);
    //    evalPnt = point3D(0., 0., h);
    evalPnt.transferGlobalToLocalCoord(evalPntProjection, panel.tangent1(), 
				       panel.tangent2(), panel.normal());
  }
  
  /**********************************************************************
   * getEvalPointDirection --
   **********************************************************************/
  template<class KT, class Panel>
  void
  CalcpForEikrOverR<KT, Panel>
  ::getEvalPointDirection (
			   const vector3D<double>& evalPnt)
  {
    if (abs(evalPnt.z()) == 0.) {
      if ((evalPointDirection_ != FOLLOW_NORMAL_DIRECTION) && 
	  (evalPointDirection_ != AGAINST_NORMAL_DIRECTION)) {
	pfft::errorMessage("surfCalcp.c : setupPanelLocalCoordSystem",
			   "illegal evalPointDirection");
      }
    } else {
      // presumably, evalPointDirection has not been specified by the user,
      // it is specified here.    
      if (evalPnt.z() > 0.) {
	evalPointDirection_ = AGAINST_NORMAL_DIRECTION;
      } else {
	evalPointDirection_ = FOLLOW_NORMAL_DIRECTION;
      }
    }
  }

  /**********************************************************************
   * setupCalcpTriangleInfo --
   **********************************************************************/
  template<class KT, class Panel>
  void
  CalcpForEikrOverR<KT, Panel>
  ::setupCalcpTriangleInfo (
			    size_t i,
			    const Panel& panel)
  {
    size_t j = (i+1) % panel.shape();
    vector3D<double> PA = panel.vertex(i);
    vector3D<double> PB = panel.vertex(j);
    
    vector3D<double> Y = PB - PA;
    triangleInfo_.AB_len = length(Y);
    Y.normalize();
    
    vector3D<double> Z(0., 0., 1.);

    triangleInfo_.X = crossProd(Y, Z);
    triangleInfo_.Y = Y;
    
    double project = - (Y * PA);
    triangleInfo_.AH = Y * project;
    triangleInfo_.AH_len = abs(project); 
  
    triangleInfo_.PH_unit = PA + triangleInfo_.AH;
    triangleInfo_.PH_len = length(triangleInfo_.PH_unit);
      if( abs(triangleInfo_.PH_len) < 
	(1e1*DBL_EPSILON * triangleInfo_.AB_len) ){
      triangleInfo_.PH_unit = vector3D<double>(0., 0., 0.);
      triangleInfo_.PH_len = 0.;
    } else {
      triangleInfo_.PH_unit /= triangleInfo_.PH_len;
    } 
    
    project = Y * PB;
    triangleInfo_.HB = Y * project;
    
    triangleInfo_.area = 0.5 * triangleInfo_.PH_len * triangleInfo_.AB_len;
    
    if (triangleInfo_.area == 0) {
      projectPointOnEdge_ = true;
    }
  }

  /**********************************************************************
   * getContributionSign --
   * 
   * It is assumed here that the left-hand rule is applied to the order 
   * of local vertices of the triangle together with the unit normal vector.
   * The system is expressed using the local coord system 
   * defined by the tangent and normal vectors of the panel.
   * The Y axis is always the vector AB, and Z axis is always the 
   * normal of the panel, which also serves as z axis in the panel local coord 
   * system. With Y znd Z fixed, the X axis is fixed. It is aligned
   * with PH, but not necessary in the same direction as PH.
   * Key points:
   * 1) Clock-wise or not is relative to the point P, so PH's direction
   *    with respect to x axis becomes the criteria.
   * 2) Clock-wise or not also depends on the relative position of the view
   *    point w.r.t. the panel normal, or the z axis in the panel local 
   *    coord system.
   *    Because of the left-hand rule in th order of panel vertices, 
   *    PH dot X < 0 ---> clock-wise.
   *    PH dot X > 0 ---> counter-clock-wise.
   * 3) When the projection point is right on the edge, there is no PH. And
   *    it does not make sense to say clock-wise or not at all. 
   *    However, it turned out that the contribution from this edge is not
   *    entirely zero. d/dX is non-zero. Using perturbation idea, we found
   *    that after perturbation, when P is a small distance away from the edge, 
   *    the sign of d/dX and the contribution sign are always the opposite if 
   *    the view point is AGAINST_NORMAL_VECTOR, or they are always the same
   *    if the view point is FOLLOW_NORMAL_VECTOR. Hence
   *    we simply let the signs of d/dX be 1 (see compDIDXInnerIntegral) 
   *    and let contribution sign in this
   *    function takes different view point into account. This is purely
   *    an artifact. Mathematically it is equivalent to other possibilities. 
   **********************************************************************/
  template<class KT, class Panel>
  double
  CalcpForEikrOverR<KT, Panel>
  ::getContributionSign (void)
  {
    double contributionSign = 1.;
    
    /* Assuming view point is AGAINST_NORMAL_VECTOR */
    if( triangleInfo_.area == 0. ) {
      contributionSign = -1.;
    } else {
      double sign = triangleInfo_.X * triangleInfo_.PH_unit;
      if (sign > 0.) {  /* counter-clock-wise */
	contributionSign = +1.;
      } else if (sign < 0.) {  /* clock-wise */
	contributionSign = -1.;
      }
    }
    
    /* flip the sign if the view point is not as assumed */
    if (evalPointDirection_ == FOLLOW_NORMAL_DIRECTION) {
      contributionSign *= -1.;
    }
    
    return contributionSign;
  }

  /**********************************************************************
   * compIntegralOnTriangle --
   * Local triangle coord. system is used in this function.
   **********************************************************************/
  template<class KT, class Panel>
  void
  CalcpForEikrOverR<KT, Panel>
  ::compIntegralOnTriangle (
			    const vector3D<double>& evalPnt, 
			    std::complex<double>& slp,
			    ComplexVector& grad)
  {
    // calculate slp, dI/dx & dI/dz. They are double integrals */
    double AH_len_withSign = triangleInfo_.AH * triangleInfo_.Y;
    if ( (AH_len_withSign <= 0.) || 
	 (triangleInfo_.AH_len >= triangleInfo_.AB_len) ) {
      // H is not on AB
      // use gauss quadrature on edge AB directly
      double length = triangleInfo_.AB_len;
      double yStart = -AH_len_withSign;
      compOuterIntegral(evalPnt, length, yStart, slp, grad.x(), grad.z());

    } else {
      // H is on AB 
#ifdef DISABLE_PIECEWISE_SCHEME
      // use gauss quadrature on edge AB directly
      // This is to show the difference made by piecewise scheme
      double length = triangleInfo_.AB_len;
      double yStart = -AH_len_withSign;
      compOuterIntegral(evalPnt, length, yStart, slp, grad.x(), grad.z());
#else
      // use gauss quadrature on segment HA and HB separately
      
      // first on HA
      std::complex<double> slpLocal1;
      std::complex<double> dIdxLocal1;
      std::complex<double> dIdzLocal1;
      double length = triangleInfo_.AH_len;
      double yStart = -length;
      compOuterIntegral(evalPnt, length, yStart, 
			slpLocal1, dIdxLocal1, dIdzLocal1);
      
      // then on HB
      // if HB is identical to HA, then I could save some CPU time here 
      std::complex<double> slpLocal2;
      std::complex<double> dIdxLocal2;
      std::complex<double> dIdzLocal2;
      double HB_len = triangleInfo_.AB_len - triangleInfo_.AH_len;
      double diff = fabs(triangleInfo_.AH_len - HB_len) / HB_len;
      if (diff <= 1.e3 * DBL_EPSILON) {
	slpLocal2 = slpLocal1;
	dIdxLocal2 = dIdxLocal1;
	dIdzLocal2 = dIdzLocal1;
      } else {
	length = HB_len;
	yStart = 0.;
	compOuterIntegral(evalPnt, length, yStart, 
			  slpLocal2, dIdxLocal2, dIdzLocal2);
      }
      slp = slpLocal1 + slpLocal2;
      grad.x() = dIdxLocal1 + dIdxLocal2;
      grad.z() = dIdzLocal1 + dIdzLocal2;
#endif
    }

    // compute dI/dY, its outer integral has been eliminated because of d/dY
    if( needGradient_ ) {
      if (triangleInfo_.area != 0.) {
	grad.y() = compDIDYOnTriangle(evalPnt);
      } else {
	grad.y() = std::complex<double>(0., 0.);
      }
    }

  }

  /**********************************************************************
   * compOuterIntegral --
   **********************************************************************/
  template<class KT, class Panel>
  void
  CalcpForEikrOverR<KT, Panel>
  ::compOuterIntegral (
		       const vector3D<double>& evalPnt, 
		       double length,
		       double yStart,
		       std::complex<double>& slp,
		       std::complex<double>& dIdx,
		       std::complex<double>& dIdz)
  {
    double h = evalPnt.z();
    double d = triangleInfo_.PH_len;
    std::complex<double> IK = std::complex<double>(0., 1.) * K_;
    double Ra = abs(h);
    
    slp = std::complex<double>(0., 0.);
    dIdx = std::complex<double>(0., 0.);
    dIdz = std::complex<double>(0., 0.);
    for (size_t jj=0; jj < NumGQpoint_; ++jj) {
      double w = length * GQweight_[jj]; 
      double y = length * GQpoint_[jj] + yStart;

      double r = sqrt(d*d + y*y);
      double Rb = sqrt(r*r + h*h);
      double R = Rb - Ra;
      std::complex<double> ZRb = exp(IK * Rb);
      std::complex<double> ZRa = exp(IK * Ra);

      if ( triangleInfo_.area != 0 ) {
	// this means the projection of the eval point is not on the edge, 
	// otherwise single layer and dIdx are all zero.

	if ( needSingleLayerPotential_ ) {
	  std::complex<double> slpLocal = compSingleLayerInnerIntegral(ZRa, ZRb, IK, R);
	  slp += slpLocal * w * d / (r*r);
	}
	    
	if ( needDoubleLayerPotential_ || needGradient_ ) {
	  std::complex<double> dIdzLocal = compDIDZInnerIntegral(ZRb, ZRa, Rb, r, d, h);
	  dIdz += dIdzLocal * w * d / (r*r);
	}
      }

      // always non-zero regardless if projection of eval point is on the edge
      if( needGradient_ ) {
	std::complex<double> tmp = compDIDXInnerIntegral(ZRb, ZRa, IK, Rb, y, R, r, d);
	dIdx += tmp * w;
      }
    }  
  }
  
  /**********************************************************************
   * compDIDXInnerIntegral --
   **********************************************************************/
  template<class KT, class Panel>
  std::complex<double>
  CalcpForEikrOverR<KT, Panel>
  ::compDIDXInnerIntegral (
			   const std::complex<double> ZRb, 
			   const std::complex<double> ZRa,
			   const std::complex<double> IK,
			   const double Rb, 
			   const double y,
			   const double R,
			   const double r, 
			   const double d)
  {
    std::complex<double> dIdx;
    
    if( d == 0. ) {
      dIdx = compSingleLayerInnerIntegral(ZRa, ZRb, IK, R);
      dIdx *= 1./(y*y);
    } else { 
      std::complex<double> tmp = compSingleLayerInnerIntegral(ZRa, ZRb, IK, R);
      tmp *= (y*y-d*d)/(r*r*r*r);
      dIdx = ZRb * (d*d)/(Rb*r*r) + tmp;
    }

    // partial d / dx term in my note
    double sign;
    if ( triangleInfo_.area != 0 ) {
      sign = triangleInfo_.PH_unit * triangleInfo_.X;
      if( sign < 0. ) {
	sign = 1.;
      } else if( sign > 0. ) {
	sign = -1.;
      }
    } else {
      sign = 1.;  // see comments of the function getContributionSign 
    }

    return dIdx * sign;
  }
 
  /**********************************************************************
   * compDIDZInnerIntegral --
   **********************************************************************/
  template<class KT, class Panel>
  std::complex<double>
  CalcpForEikrOverR<KT, Panel>
  ::compDIDZInnerIntegral (
			   const std::complex<double> ZRb, 
			   const std::complex<double> ZRa,
			   const double Rb, 
			   const double r, 
			   const double d, 
			   const double h)
  {
    double sign;
    if (h > 0.) {
      sign = 1.;
    } else if (h < 0.){
      sign = -1.;
    } else {
      if (evalPointDirection_ == AGAINST_NORMAL_DIRECTION) {
	sign = 1.;
      } else {
	sign = -1.;
      }
    }
    
    // grad = ZRb dRb/dn - ZRa dRa/dn = ZRb * h / Rb - ZRa * sign
    return ZRb * h / Rb - ZRa * sign;
  }

  /**********************************************************************
   * compDIDYOnTriangle --
   * dI/dy = -f(ThetaB,h) * d / (rB^2) + f(ThetaA,h) * d / (rA^2). 
   * The outer integral has disapeared because of d/dy.
   **********************************************************************/
  template<class KT, class Panel>
  std::complex<double>
  CalcpForEikrOverR<KT, Panel>
  ::compDIDYOnTriangle(
		       const vector3D<double>& evalPnt)
  {
    double h = evalPnt.z();
    double d = triangleInfo_.PH_len;
    double LA = triangleInfo_.AH_len;
    double LB = length(triangleInfo_.HB);
    
    double rB2 = LB*LB + d*d;
    double rA2 = LA*LA + d*d;
    double RB = sqrt(rB2 + h*h);
    double RA = sqrt(rA2 + h*h);
    double R0 = fabs(h);
    
    std::complex<double> IK = std::complex<double>(0., 1.) * K_;
    std::complex<double> ZRB = exp(IK * RB);
    std::complex<double> ZRA = exp(IK * RA);
    std::complex<double> ZR0 = exp(IK * R0);

    std::complex<double> fThetaA = compSingleLayerInnerIntegral(ZR0, ZRA, IK, RA-R0);
    std::complex<double> fThetaB = compSingleLayerInnerIntegral(ZR0, ZRB, IK, RB-R0);
    return -fThetaB * (d/rB2) + fThetaA * (d/rA2);
  }

  /**********************************************************************
   * compSingleLayerInnerIntegral --
   **********************************************************************/
  template<class KT, class Panel>
  std::complex<double>
  CalcpForEikrOverR<KT, Panel>
  ::compSingleLayerInnerIntegral (
				  const std::complex<double> ZRa, 
				  const std::complex<double> ZRb,
				  const std::complex<double> IK,
				  const double R)
  {
    std::complex<double> slp;

    std::complex<double> IKR = IK * R;
    if ( abs(IKR) > 1.e-3 ) {
      slp = (ZRb - ZRa) / IK;
    } else {
      std::complex<double> tmp = std::complex<double>(R, 0.);
      slp = tmp;
      for (int m=2; m<=4; m++) {
	tmp *= IKR / static_cast<double>(m);
	slp += tmp;
      }
      slp *= ZRa;
    }
    
    return slp;
  }

  /**********************************************************************
   * flipSign --
   * Even though getContributionSign has taken care of clocok or anti-clock
   * flag, this is dependent upon the local order of panel vertices.
   * If everything unchanged except the oreder of the local panel vertices,
   * I might get a negative single-layer potential. To get rid of this
   * dependence, I use the toalArea as indicator since it should always be
   * positive. This will take some tiny extra CPU time. But it is worthwhile.
   **********************************************************************/
  template<class KT, class Panel>
  void
  CalcpForEikrOverR<KT, Panel>
  ::flipSign (
	      std::complex<double>& slp,
	      ComplexVector& grad)
  {
    slp *= -1.;
    grad *= -1.;
  }
  
  /**********************************************************************
   * localDoubleLayerCorrection --
   * For the case when  evalPnt is on the plane where the source panel lies,
   * it is known that 
   * Int dG/dn' = (+ or -) 2*PAI, if evalPnt is inside the panel, 
   * Int dG/dn' == 0, if evalPnt is outside of the panel.
   * However the simple summation will not produce exact zero or 2*PAI,
   * to maintain high accuracy, I modify the results here since I already 
   * know them.
   **********************************************************************/
  template<class KT, class Panel>
  void
  CalcpForEikrOverR<KT, Panel>
  ::localDoubleLayerCorrection (
				std::complex<double>& dfdn)
  {
    double realPart = dfdn.real();
    if (realPart > 1.85 * PAI) {
      dfdn = std::complex<double>(2.*PAI, 0.);
    } else if (realPart < -1.85*PAI) {
      dfdn = std::complex<double>(-2.*PAI, 0.);
    } else if ((realPart > -0.1) && (realPart < 0.1) ) {
      dfdn = std::complex<double>(0., 0.);
    } else if (projectPointOnEdge_) {
      // doubleLayer should be PAI
      if (abs(realPart/PAI - 1.) < 5e-2) {
	dfdn = std::complex<double>(PAI, 0.);
      } else if (abs(realPart/(-PAI) - 1.) < 5e-2) {
	dfdn = std::complex<double>(-PAI, 0.);
      }
    } else {
      pfft::warningMessage("calcpForEikrOverR.h : localDoubleLayerCorrection",
			   "Inaccurate double-layer potential due to large panel aspect ratio. Automatically corrected");

      if (realPart > 0.1) {
	dfdn = std::complex<double>(2.*PAI, 0.);
      } else if (realPart < -0.1) {
	dfdn = std::complex<double>(-2.*PAI, 0.);
      } else {
	dfdn = std::complex<double>(0., 0.);
      } 

#ifdef DEBUG
      std::cout << dfdn << endl;
      throw domain_error("calcpForEikrOverR.h : localDoubleLayerCorrection");
#endif

    }
  }

  /**********************************************************************
   * transferToPanelLocalCoordSystem --
   * The grad1 is w.r.t. the local triangle coord system
   * The grad2 is w.r.t. the local [anel coord system
   * x2 = panel.tangent1, y2=panel.tangent2, and z2=panel.normal;
   * Please note z1 = z2.
   * 
   * df/dx2 = df/d(tan1) = df/dx1*x1.x2 + df/dy1*y1.x2 + df/dz1*z1.x2(=0)
   * df/dy2 = df/d(tan2) = df/dx1*x1.y2 + df/dy1*y1.y2 + df/dz1*z1.y2(=0)
   * df/dz2 = df/dz1
   **********************************************************************/
  template<class KT, class Panel>
  typename CalcpForEikrOverR<KT, Panel>::ComplexVector
  CalcpForEikrOverR<KT, Panel>
  ::transferToPanelLocalCoordSystem (
				     const ComplexVector& grad1)
  {
    return ComplexVector(grad1.x() * triangleInfo_.X.x() + 
  			 grad1.y() * triangleInfo_.Y.x(),
  			 grad1.x() * triangleInfo_.X.y() + 
  			 grad1.y() * triangleInfo_.Y.y(),
  			 grad1.z() );
  }			   

  /**********************************************************************
   * gaussQuadrature --
   **********************************************************************/
  template<class KT, class Panel>
  void 
  CalcpForEikrOverR<KT, Panel>
  ::gaussQuadrature (
		     const size_t n, // number of points
		     std::vector<double>& w,  // weights
		     std::vector<double>& x)   // points
  {  
    double x1 = 0;
    double x2 = 1;
    
    int m = (n+1)/2;
    double xm = .5 * (x2 + x1); 
    double xl = .5 * (x2 - x1);
    double pp;

    for (int i = 1; i <= m; i++) {

      // Initial value of the ith root 
      double z = cos(PAI * (i-.25) / (n+.5));
      
      // Get the ith root by Newton method 
      bool converge = false;

      do  {

	// Evaluate Leg polynomial. p1 is the value at z
	double p1 = 1.0; 
	double p2 = 0.0;
	for (int j=1; j <= n; j++) {
	  double p3 = p2; 
	  p2 = p1; 
	  p1 = ((2.0*j-1.0) * z * p2 - (j-1.0) * p3) / j;
	}

	// pp is the derivative 
	pp = n * (z*p1 - p2) / (z*z - 1.0);
	double z1 = z;
	z = z1 - p1/pp;
	
	if (abs(z-z1) < 1e1*DBL_EPSILON ) {
	  converge = true;
	} 

      } while( !converge);
      
      x[i-1] = xm - xl*z;
      x[n+1-i-1] = xm + xl*z;
      w[i-1] = 2. * xl / ((1 - z*z) * pp * pp);
      w[n+1-i-1] = w[i-1];
    }
  } 

} //pfft

#endif  
  


