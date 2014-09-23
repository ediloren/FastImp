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

  const static char cvsid[] = "$Id: utils.h,v 1.4 2003/03/25 15:14:13 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _UTILS_H_
#define _UTILS_H_

#include <string>
#include <complex>
#include "vector3D.h"

namespace pfft {

  void timeReport (const std::string& message, double time);
  void errorMessage (char *fileName, char *message);
  void progressReport (const std::string& message);
  void warningMessage (char *fileName, char *message);
  template <class T> 
  void warningMessage (char *fileName, 
		       const std::string& message, 
		       const T& value) {
    std::cout << std::endl
	      << "\t Warning in " << fileName << std::endl
	      << "\t " << message << "  " << value << std::endl 
	      << std::endl;
  }
  template <class T> 
  void message (const std::string& message, const T& value) {
#ifdef PRINT_PROGRESS
    std::cout << "\t" << message << "  " << value << std::endl;
#endif
  }

  const double PAI = 3.1415926535897932385;

  const std::complex<double> J(0.,1.);
  const std::complex<double> COMPLEX_ZERO(0.,0.);
  const pfft::vector3D<double> ZERO_REAL_VECTOR_3D(0., 0., 0.);
  const pfft::vector3D<std::complex<double> > ZERO_COMPLEX_VECTOR_3D(COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO); 

}

#endif
