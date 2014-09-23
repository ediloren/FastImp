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

  const static char cvsid[] = "$Id: eikrOverR.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _EIKR_OVER_R_H_
#define _EIKR_OVER_R_H_

#include "vector3D.h"
#include <complex>

namespace pfft {

  template <class KT = std::complex<double> >

  class EikrOverR {

  public:
    EikrOverR(KT KIn) : K_(KIn) {}
    EikrOverR(void) : K_(0) {}

    std::complex<double> operator () (const double x, const double y, const double z) const {
      vector3D<double> r(x, y, z); 
      if (length(r) == 0) {
	return 0.0;
      } else {
	return exp(std::complex<double>(0., 1.) * K_ * length(r)) / length(r); 
      }
    }

    std::complex<double> operator () (const vector3D<double>& r) const {
      if (length(r) == 0) {
	return 0.0;
      } else {
	return exp(std::complex<double>(0., 1.) * K_ * length(r)) / length(r); 
      }
    }

    const KT K(void) const { return K_; }

  private:
    KT K_; // k as in exp(ikR)/R
  };

} //namespace pfft

#endif
