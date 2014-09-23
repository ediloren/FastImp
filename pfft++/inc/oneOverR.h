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

  const static char cvsid[] = "$Id: oneOverR.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _ONE_OVER_R_H_
#define _ONE_OVER_R_H_

#include "vector3D.h"

namespace pfft {

  class OneOverR {

  public:

    OneOverR(void) {}

    double operator () (const double x, const double y, const double z) const {
      vector3D<double> r(x, y, z);    
      if (length(r) == 0) {
	return 0.0;
      } else {
	return 1. / length(r);
      }
    }

    double operator () (const vector3D<double>& r) const { 
      if (length(r) == 0) {
	return 0.0;
      } else {
	return 1. / length(r);
      }
    }

  };

} //namespace pfft

#endif
