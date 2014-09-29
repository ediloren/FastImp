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

  const static char cvsid[] = "$Id: service.cc,v 1.6 2003/07/10 21:25:59 zhzhu Exp $";

  ==========================================================================
*/

#include "service.h"
#include <iostream>
// Enrico
#include <stdlib.h>

using namespace surf;
using namespace std;

/**********************************************************************
 * errorMessage --
 **********************************************************************/
void
surf::errorMessage (
		    const std::string& fileName,
		    const std::string& message)
{
  std::cout << std::endl
       << "\t Error in " << fileName << std::endl
       << "\t " << message << std::endl 
       << "\t Exit ... " <<  std::endl;
  exit(1);
}

/**********************************************************************
 * warningMessage --
 **********************************************************************/
void
surf::warningMessage (
		      const std::string& fileName,
		      const std::string& message)
{
  std::cout << std::endl
       << "\t Warning in " << fileName << std::endl
       << "\t " << message << std::endl;
}

/**********************************************************************
 * progressReport --
 **********************************************************************/
void
surf::progressReport (
		      const std::string& message)
{
#ifdef PRINT_PROGRESS
  std::cout << std::endl << message << std::endl;
#endif
}

/**********************************************************************
 * timeReport --
 **********************************************************************/
void 
surf::timeReport (
		  const std::string& message,
		  double time)
{
#ifdef PRINT_TIME
  if (time < 0) time = 0.;
  std::cout << std::endl << "\t" << message 
	    << time << " (seconds)" << std::endl;
#endif
}


