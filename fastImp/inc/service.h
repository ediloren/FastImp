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

#ifndef __SERVICE_H_
#define __SERVICE_H_

#include <string>

namespace surf {

  void errorMessage (const std::string& fileName, const std::string& message);
  void warningMessage (const std::string& fileName, const std::string& message);
  void progressReport (const std::string& message);
  void timeReport (const std::string& message, double time);

  template <class T> 
  void messageWithOneNumber (const std::string& message, const T& value) {
#ifdef PRINT_PROGRESS
    std::cout << "\t" << message << "  " << value << std::endl;
#endif
  }

  template <class T> 
  void message (const std::string& message, const T& v1, const T& v2) {
#ifdef PRINT_PROGRESS
    std::cout << "\t" << message << "  " << v1 << "  " << v2 << std::endl;
#endif
  }
} // namespace surf

#endif


