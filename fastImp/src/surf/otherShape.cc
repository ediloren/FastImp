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

  const static char cvsid[] = "$Id: otherShape.cc,v 1.2 2002/07/18 15:05:45 zhzhu Exp $";

  ==========================================================================
*/

#include "otherShape.h"
#include "fileIo.h"

using namespace std;
using namespace surf;

/**********************************************************************
 * input --
 **********************************************************************/
void
OtherShape::readInputFile (
			   void)
{
  Rd_Double(&conductivity_);
}

