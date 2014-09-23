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

  const static char cvsid[] = "$Id: wire.cc,v 1.3 2002/07/18 15:05:45 zhzhu Exp $";

  ==========================================================================
*/


#include "wire.h"
#include "fileIo.h"

using namespace std;
using namespace surf;

/**********************************************************************
 * input --
 **********************************************************************/
void
Wire::readInputFile (
		     void)
{
  leftEndPoint0_ = readPoint(unit_);
  leftEndPoint1_ = readPoint(unit_);
  leftEndPoint2_ = readPoint(unit_);
  rightEndPoint_ = readPoint(unit_);
  Rd_Int(&numPanelBetweenPoint01_);
  Rd_Int(&numPanelBetweenPoint02_);
  Rd_Int(&numPanelLength_);
  Rd_Double(&conductivity_);
  
  double p1 = length(leftEndPoint0_ - leftEndPoint1_) / numPanelBetweenPoint01_;
  double p2 = length(leftEndPoint0_ - leftEndPoint2_) / numPanelBetweenPoint02_;
  double p3 = length(leftEndPoint0_ - rightEndPoint_) / numPanelLength_;

  maxPanelSize_ = p1 > p2 ? p1 : p2;
  maxPanelSize_ = maxPanelSize_ > p3 ? maxPanelSize_ : p3;

}
