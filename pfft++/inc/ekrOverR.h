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

  const static char cvsid[] = "$Id: ekrOverR.h,v 1.2 2003/01/10 19:41:04 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _EKR_OVER_R_H_
#define _EKR_OVER_R_H_

#include "vector3D.h"

namespace pfft {

  class EkrOverR {

  public:
    EkrOverR(double KIn) : K_(KIn) {}
    EkrOverR(void) : K_(0) {}

    double operator () (const double x, const double y, const double z) const {
      vector3D<double> r(x, y, z); 
      if (length(r) == 0) {
	return 0.0;
      } else {
	return exp(K_ * length(r)) / length(r); 
      }
    }

    double operator () (const vector3D<double>& r) const {
      if (length(r) == 0) {
	return 0.0;
      } else {
	return exp(K_ * length(r)) / length(r); 
      }
    }

    const double K(void) const { return K_; }

  private:
    double K_; // k as in exp(kR)/R
  };

} //namespace pfft

#endif
