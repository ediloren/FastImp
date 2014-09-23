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

  const static char cvsid[] = "$Id: spiral.cc,v 1.3 2002/07/18 15:05:45 zhzhu Exp $";

  ==========================================================================
*/

#include "spiral.h"
#include "fileIo.h"

using namespace std;
using namespace surf;

/**********************************************************************
 * input --
 **********************************************************************/
void
Spiral::readInputFile (
		       void)
{
  origin1_ = readPoint(unit_);
  origin2_ = readPoint(unit_);
  edgePoint_ = readPoint(unit_);
  innerRadius_ = readSize(unit_);
  outerRadius_ = readSize(unit_);
  r3_ = readSize(unit_);
  Rd_Double(&numRound_);
  Rd_Int(&numPanelRadius_);
  Rd_Int(&numPanelThickness_);
  Rd_Int(&numPanelArc_);
  Rd_Double(&conductivity_);
  
  double p1 = length(origin1_ - origin2_) / numPanelThickness_;
  double p2 = (outerRadius_ - innerRadius_) / numPanelRadius_;
  double r2 = outerRadius_;
  double r3 = r3_;
  double length = 2. * 3.14 * r2 * numRound_ * (1. + (r3-r2)/r2 * numRound_ / 2);
  double p3 = length / numPanelArc_;

  maxPanelSize_ = p1 > p2 ? p1 : p2;
  maxPanelSize_ = maxPanelSize_ > p3 ? maxPanelSize_ : p3;

}
