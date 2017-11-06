/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Ben Song and Zhenhai Zhu
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: calcpForOneOverR.h,v 1.12 2003/03/25 15:14:13 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _CALCP_FOR_ONE_OVER_R_
#define _CALCP_FOR_ONE_OVER_R_

#include <cfloat> // for DBL_MIN, DBL_MAX, DBL_EPSILON
#include "utils.h" // for PAI
#include <assert.h>
// Enrico, for std::min and std::max
#include <algorithm>

namespace pfft {

  template<class Panel = pfft::element>
  class CalcpForOneOverR {

  public:
    CalcpForOneOverR(void) {};
    CalcpForOneOverR (const bool IsSlpNeeded,
		      const bool IsDlpNeeded,
		      const bool IsGradOfSlpNeeded,
		      const bool IsGradOfDlpNeeded) 
      : needSingleLayer_(IsSlpNeeded),
	needDoubleLayer_(IsDlpNeeded),
	needGradSingleLayer_(IsGradOfSlpNeeded),
	needGradDoubleLayer_(IsGradOfDlpNeeded),
	previousPanelCentroid_(point3D(DBL_MAX, DBL_MIN, 0.)) { };

    void operator() (const Panel&,
		     const vector3D<double>& evalPnt,
		     double& slp,
		     double& dlp,
		     vector3D<double>& gradOfSlp,
		     vector3D<double>& gradOfDlp); 

  private:
    bool needSingleLayer_;
    bool needDoubleLayer_;
    bool needGradSingleLayer_;
    bool needGradDoubleLayer_;
    point3D previousPanelCentroid_;

    enum PanelVertexOrder { LEFT_HAND_RULE, RIGHT_HAND_RULE };
    PanelVertexOrder panelVertexOrder_;

    typedef struct PanelInfo_ {
      vector3D<double> X;
      vector3D<double> Y;
      vector3D<double> Z;
      vector3D<double> centroid;
      
      double dtol;
      double minSideLength;
      double maxSideLength;
      //      std::vector<double> ct;
      //      std::vector<double> st;
      //      std::vector<double> edgeLength;
      double ct[4];
      double st[4];
      double edgeLength[4];
      std::vector<vector3D<double> > vertex;  
    } PanelInfo_;
    PanelInfo_ panelInfo_;

    bool isEvalPntInsideSrcPanel (const Panel&,
				  const vector3D<double>& evalPnt) const;
    void findPanelVertexOrder(const Panel&);
    void setupPanelInfo(const Panel&);
    bool isSamePanelAsPreviousOne(const point3D& currentSrcPanelCentroid) { 
      return (previousPanelCentroid_ == currentSrcPanelCentroid);
    }
  };

  /**********************************************************************
   * operator --
   * To setup panel info takes time. So the idea is to minize it.
   * I use the centroid of the src panel as its identity (since it is 
   * unique to each panel) and decide if the current srcPanel is different
   * from previous one.
   * It is obvious that the ideal way of calling this fucntion is to
   * loop through eval points in the inner loop and loop through srcPanel 
   * in the outer loop. This way we could "vectorize" the whole process
   * and significantly reduce the CPU time.
   *
   * The class Panel should at least have the following interface functions:
   * shape(), centroid(), vertex(i).
   **********************************************************************/
  template<class Panel>
  void
  CalcpForOneOverR<Panel> :: 
  operator() (
	      const Panel& srcPanel,
	      const vector3D<double>& evalPnt,
	      double& slp,
	      double& dlp,
	      vector3D<double>& gradOfSlp,
	      vector3D<double>& gradOfDlp)
  {
    // only deal with triangle and quadratic panel. 
    assert((srcPanel.shape() == 3) || (srcPanel.shape() == 4));

    if (! isSamePanelAsPreviousOne(srcPanel.centroid())) {
      previousPanelCentroid_ = srcPanel.centroid();
      findPanelVertexOrder(srcPanel);
      setupPanelInfo(srcPanel);
    }

    // First map the evaluation point into panel's coord system.
    vector3D<double> V1 = evalPnt - panelInfo_.centroid;
    vector3D<double> V0(panelInfo_.X*V1, panelInfo_.Y*V1, panelInfo_.Z*V1);

    double zn = V0.z();
    double znabs = fabs(zn);
    double evalDistance = length(V0);
    
    /* Once per vertex computation. */
    int OK, next;
    double v, fac, fd, fh1, fh2;
    double u1, u2, rr;
    double s1, c1, s2, c2, s12, c12, val;
    double intOfdGdX, intOfdGdY, intOfdGdZ;
    double dXIntOfDlp, dYIntOfDlp, dZIntOfDlp;
    double xmxv[4], ymyv[4], fe[4], r[4], xri[4], yri[4];

    OK = 1;
    size_t numOfVertex = srcPanel.shape();
    for(size_t i=0; i < numOfVertex; ++i) {
      double xc = V0.x() - panelInfo_.vertex[i].x();
      double yc = V0.y() - panelInfo_.vertex[i].y();
      double zc = V0.z() - panelInfo_.vertex[i].z();
      
      xmxv[i] = xc; 
      ymyv[i] = yc;
    
      fe[i] = xc*xc + zc*zc;
      r[i] = sqrt(yc*yc + fe[i]);

      if(r[i] < (1.005*znabs)) {
	OK = 0; 
      }
      if(needDoubleLayer_ || needGradSingleLayer_ || needGradDoubleLayer_ ) {
	xri[i] = xmxv[i]/r[i];
	yri[i] = ymyv[i]/r[i];
      }
    }
  
    /* Potential and dipole formed by summing contributions from each edge */
    slp = 0.0; 
    dlp = 0.0;
    if (needDoubleLayer_ || needGradSingleLayer_ || needGradDoubleLayer_ ) {
      intOfdGdX = 0.;
      intOfdGdY = 0.;
      dXIntOfDlp = 0.;
      dYIntOfDlp = 0.;
      dZIntOfDlp = 0.;
    }
    
    for(size_t i=0; i <numOfVertex; ++i) {
      if (i == (numOfVertex-1)) next = 0;
      else next = i + 1;
    
      /* v is the projection of the eval-i edge on the perpend to the side-i */
      /* Exploits the fact that corner points in panel coordinates.  */
      v = xmxv[i] * panelInfo_.st[i] - ymyv[i] * panelInfo_.ct[i];
    
      /* arg == zero if eval on next-i edge, but then v = 0. */
      double arg = (r[i]+r[next]-panelInfo_.edgeLength[i]) / 
	(r[i]+r[next]+panelInfo_.edgeLength[i]);
    
      double fln = -log(arg);
      if (arg > 0.0) {
	slp += v * fln;
	if (needDoubleLayer_ || needGradSingleLayer_ || needGradDoubleLayer_) {
	  fac = (r[i] + r[next] - panelInfo_.edgeLength[i])
	    * (r[i] + r[next] + panelInfo_.edgeLength[i]);
	  fac = v * (2*panelInfo_.edgeLength[i]) / fac;
	  intOfdGdX += fln * panelInfo_.st[i] - fac*(xri[i] + xri[next]);
	  intOfdGdY -= fln * panelInfo_.ct[i] + fac*(yri[i] + yri[next]);
	  dZIntOfDlp -= fac * (1.0/r[i] + 1.0/r[next]);
	}
      }
    
      /* OK means eval not near a vertex normal, use Hess-Smith: */
      if (OK == 1) {
	s1 = v * r[i];
	c1 = znabs * (xmxv[i] * panelInfo_.ct[i] + ymyv[i] * panelInfo_.st[i]);
	s2 = v * r[next];
	c2 = znabs * (xmxv[next] * panelInfo_.ct[i] + ymyv[next] * panelInfo_.st[i]);
      } else { /* Near a vertex normal, use Newman  */
	s1 = (fe[i] * panelInfo_.st[i])-(xmxv[i] * ymyv[i] * panelInfo_.ct[i]);
	c1 = znabs*r[i] * panelInfo_.ct[i];
	s2 = (fe[next]*panelInfo_.st[i])-(xmxv[next]*ymyv[next]*panelInfo_.ct[i]);
	c2 = znabs*r[next]*panelInfo_.ct[i];
      }
    
      s12 = (s1*c2)-(s2*c1);
      c12 = (c1*c2)+(s1*s2);
      val = atan2(s12, c12);
      dlp += val;
      if (needDoubleLayer_ || needGradSingleLayer_ || needGradDoubleLayer_) {
        u1 = xmxv[i] * panelInfo_.ct[i] + ymyv[i] * panelInfo_.st[i];
        u2 = xmxv[next] * panelInfo_.ct[i] + ymyv[next] * panelInfo_.st[i];
        if (OK == 0) {
          rr  = r[i] * r[i];
          fh1 = xmxv[i] * ymyv[i];
          fh2 = xmxv[next] * ymyv[next];
          fac = c1/((c1*c1+s1*s1)*rr );
          if(zn < 0.0) fac = -1.0 * fac;  /* Nick's 8/98 correpanelInfo_.ction. */
          dXIntOfDlp += (rr*v+fh1*u1)*fac;
          dYIntOfDlp -= fe[i]*u1*fac;
          rr  = r[next]*r[next];
          fac = c2/((c2*c2+s2*s2)*rr);
          if(zn < 0.0) fac = -1.0 * fac; /* Nick's 8/98 correpanelInfo_.ction. */
          dXIntOfDlp -= (rr*v+fh2*u2)*fac;
          dYIntOfDlp += fe[next]*u2*fac;
        } else {
	  if (zn <= 3*DBL_EPSILON) {
	    fac = 0.;
	  } else {
	    if (c1*c1+s1*s1 == 0.) {
	      pfft::errorMessage("calcpForOneOverR::operator()",
				 "Divided by zero, there must be bug here!! ");
	    } else {
	      fac = zn/(c1*c1+s1*s1);
	    }
	  }
          dXIntOfDlp += (u1*v*xri[i]+r[i]*ymyv[i])*fac;
          dYIntOfDlp += (u1*v*yri[i]-r[i]*xmxv[i])*fac;

	  if (zn <= 3*DBL_EPSILON) {
	    fac = 0.;
	  } else {
	    if (c2*c2+s2*s2 == 0.) {
	      pfft::errorMessage("calcpForOneOverR::operator()",
				 "Divided by zero, there must be bug here!! ");
	    } else {
	      fac = zn/(c2*c2+s2*s2);
	    }
	  }
          dXIntOfDlp -= (u2*v*xri[next]+r[next]*ymyv[next])*fac;
          dYIntOfDlp -= (u2*v*yri[next]-r[next]*xmxv[next])*fac;
        }
      }
    }

    if (dlp < 0.0) dlp += 2.0 * PAI; 
    if (zn < 0) dlp *= -1.0;
    slp -= zn*dlp;

    if (fabs(V0.z()) < panelInfo_.dtol &&
	isEvalPntInsideSrcPanel (srcPanel, evalPnt) ) {
      // if eval pnt is on the src panel, the  
      // double layer potential should be 2*PAI.
      dlp = -2.*PAI;
    }
  
    if (panelVertexOrder_ != LEFT_HAND_RULE) {
      dlp *= -1.0;
    }

    if(needDoubleLayer_ || needGradSingleLayer_ || needGradDoubleLayer_) {
      intOfdGdX -= zn*dXIntOfDlp;
      intOfdGdY -= zn*dYIntOfDlp;
    }
  
    /* intOfdGdX means the d(slp)/d(local X)
       intOfdGdY means the d(slp)/d(local Y)
       Similarly, 
       dXIntOfDlp means the d(dlp)/d(local X)
       dYIntOfDlp means the d(dlp)/d(local Y)
       so, I have to transfer them back to global coordinate system.
    */

    /* Save the computed potentials and derivatives. */
    if (needGradSingleLayer_) {
      gradOfSlp.x()
	= intOfdGdX * panelInfo_.X.x()
	+ intOfdGdY * panelInfo_.Y.x()
	- dlp * panelInfo_.Z.x();
      
      gradOfSlp.y() 
	= intOfdGdX * panelInfo_.X.y()
	+ intOfdGdY * panelInfo_.Y.y()
	- dlp * panelInfo_.Z.y();
      
      gradOfSlp.z()
	= intOfdGdX * panelInfo_.X.z()
	+ intOfdGdY * panelInfo_.Y.z()
	- dlp * panelInfo_.Z.z();
    }
    
    if (needGradDoubleLayer_) {
      gradOfDlp.x() 
	= dXIntOfDlp * panelInfo_.X.x()
	+ dYIntOfDlp * panelInfo_.Y.x()
	+ dZIntOfDlp * panelInfo_.Z.x();
    
      gradOfDlp.y() 
	= dXIntOfDlp * panelInfo_.X.y()
	+ dYIntOfDlp * panelInfo_.Y.y()
	+ dZIntOfDlp * panelInfo_.Z.y();

      gradOfDlp.z()
	= dXIntOfDlp * panelInfo_.X.z()
	+ dYIntOfDlp * panelInfo_.Y.z()
	+ dZIntOfDlp * panelInfo_.Z.z();
    }
  }
 
  /**********************************************************************
   * findPanelVertexOrder --
   **********************************************************************/
  template<class Panel>
  void
  CalcpForOneOverR<Panel>::findPanelVertexOrder(
						const Panel& anyPanel)
  {
    vector3D<double> v01 = anyPanel.vertex(1) - anyPanel.vertex(0);
    vector3D<double> v02 = anyPanel.vertex(2) - anyPanel.vertex(1);

    double tmp = crossProd(v01, v02) * anyPanel.normal();
    if (tmp < 0) {
      panelVertexOrder_ = LEFT_HAND_RULE; // LeftHand(0->1->2) = normal 
    } else {
      panelVertexOrder_ = RIGHT_HAND_RULE; // RightHand(0->1->2) = normal
    }
  }

  /**********************************************************************
   * isEvalPntInsideSrcPanel --
   **********************************************************************/
  template<class Panel>
  bool
  CalcpForOneOverR<Panel>::
  isEvalPntInsideSrcPanel (
			   const Panel& srcPanel,
			   const vector3D<double>& evalPnt) const
  {
    const size_t numOfVertex = srcPanel.shape();
    double sum = 0;
    for (size_t ii = 1; ii <= numOfVertex; ii ++) {
      vector3D<double> e1 = evalPnt - srcPanel.vertex((ii-1) % numOfVertex);
      vector3D<double> e2 = evalPnt - srcPanel.vertex((ii) % numOfVertex);
      double costheta;
      dotProd(costheta, e1, e2);
      if (costheta != 0.) {
	costheta /= length(e1) * length(e2);
      } 
      sum += acos(costheta); 
    }
  
    // Enrico, calling double std::abs(double) instead of casting to int abs(int)
//    if (abs(sum - 2*PAI) < 1.e-6) {
    if (std::abs(sum - 2*PAI) < 1.e-6) {
      return true;
    } 
    return false;
  }
 
  /**********************************************************************
   * setupPanelInfo --
   **********************************************************************/
  template<class Panel>
  void 
  CalcpForOneOverR<Panel> :: 
  setupPanelInfo (
		  const Panel& srcPanel)
  {
    size_t numOfVertex = srcPanel.shape();
    //    panelInfo_.edgeLength.resize(numOfVertex);
    //    panelInfo_.ct.resize(numOfVertex);
    //    panelInfo_.st.resize(numOfVertex);
    panelInfo_.vertex.resize(numOfVertex);

    panelInfo_.minSideLength = DBL_MAX;
    panelInfo_.maxSideLength = DBL_MIN;

    for(size_t ii=0; ii < numOfVertex; ii++) {
      size_t j = (ii+1) % numOfVertex;
      panelInfo_.edgeLength[ii] 
	= length(srcPanel.vertex(j) - srcPanel.vertex(ii));
	  // Enrico, changed to standard template function
      //panelInfo_.minSideLength = min(panelInfo_.minSideLength,
      panelInfo_.minSideLength = std::min(panelInfo_.minSideLength,
				     panelInfo_.edgeLength[ii]);
	  // Enrico, changed to standard template function
      //panelInfo_.maxSideLength = max(panelInfo_.maxSideLength,
      panelInfo_.maxSideLength = std::max(panelInfo_.maxSideLength,
				     panelInfo_.edgeLength[ii]);
    };
  
    panelInfo_.X = srcPanel.vertex(2) - srcPanel.vertex(0);
    double xnorm = length(panelInfo_.X);
  
    if (numOfVertex == 3) {
      panelInfo_.Y = srcPanel.vertex(1) - srcPanel.vertex(0);
    } else {
      panelInfo_.Y = srcPanel.vertex(1) - srcPanel.vertex(3);
    }
  
    double ynorm = length(panelInfo_.Y);
    double minDiagLength;
    if (numOfVertex == 3) {
      minDiagLength = panelInfo_.minSideLength;
    } else {
      minDiagLength = xnorm < ynorm ? xnorm : ynorm;
    }

    panelInfo_.dtol = 1e6 * DBL_EPSILON * minDiagLength;
  
    crossProd(panelInfo_.Z, panelInfo_.X, panelInfo_.Y);
    double znorm = length(panelInfo_.Z);

    double area = .5 * znorm;
    panelInfo_.Z /= znorm;
    panelInfo_.X /= xnorm;
    crossProd(panelInfo_.Y, panelInfo_.Z, panelInfo_.X);

    // put panel vertices in newly defined coordinate system.
    panelInfo_.centroid = srcPanel.centroid();
    for (size_t i=0; i<numOfVertex; ++i) {
      vector3D<double> V0 = srcPanel.vertex(i) - panelInfo_.centroid;
      panelInfo_.vertex[i].x() = panelInfo_.X * V0;
      panelInfo_.vertex[i].y() = panelInfo_.Y * V0;
      panelInfo_.vertex[i].z() = panelInfo_.Z * V0;
    }

    // Check that panel is in the x-y plane
    for(size_t i=0; i<numOfVertex; ++i) {
      if ( fabs(panelInfo_.vertex[i].z()) >= (1.0e-8 * xnorm) ) {
	std::cout << "xnorm := " << xnorm << std::endl;
	std::cout << "panelInfo_.vertex[i].z() := " << panelInfo_.vertex[i].z() << std::endl;
	errorMessage("CalcpForOneOverR::setupPanelInfo", "vertices are not on the same plane");
      }
    }

    // compute the contribution terms for each edge.
    for(size_t i = 0; i < numOfVertex; i++) {
      size_t j = (i+1) % numOfVertex;
      panelInfo_.ct[i] 
	= (panelInfo_.vertex[j].x() - panelInfo_.vertex[i].x()) / panelInfo_.edgeLength[i];

      panelInfo_.st[i] 
	= (panelInfo_.vertex[j].y() - panelInfo_.vertex[i].y()) / panelInfo_.edgeLength[i];
    }
  }

} // namespace pfft

#endif
