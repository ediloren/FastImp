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

  const static char cvsid[] = "$Id: utils.cc,v 1.3 2003/02/04 20:59:33 zhzhu Exp $";

  ==========================================================================
*/

#include "utils.h"
#include <iostream>

using namespace pfft;
using namespace std;

/**********************************************************************
 * timeReport --
 **********************************************************************/
void 
pfft::timeReport (
		  const std::string& message,
		  double time)
{
#ifdef TIME
  if (time < 0) time = 0.;
  cout << endl << message << time << " (seconds)" << endl;
#endif
}

/**********************************************************************
 * errorMessage --
 **********************************************************************/
void
pfft::errorMessage (
		    char *fileName,
		    char *message)
{
  cout << endl
       << "\t Error in " << fileName << endl
       << "\t " << message << endl 
       << "\t Exit ... " <<  endl;
  exit(1);
}
  
/**********************************************************************
 * warningMessage --
 **********************************************************************/
void
pfft::warningMessage (
		      char *fileName,
		      char *message)
{
  cout << endl
       << "\t Warning in " << fileName << endl
       << "\t " << message << endl << endl;
}

/**********************************************************************
 * progressReport --
 **********************************************************************/
void
pfft::progressReport (
		      const std::string& message)
{
#ifdef PRINT_PROGRESS
  cout << endl << message << endl;
#endif
}

