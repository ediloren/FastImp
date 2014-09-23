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

  const static char cvsid[] = "$Id: g2gPoint.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _G2G_POINT_H_
#define _G2G_POINT_H_

#include "gridIndex.h"

namespace pfft {

  class G2GPoint {

  public:

    G2GPoint(const GridIndex& pos, int directStencilPointIndex, 
	     int interpStencilPointIndex) 
      : pos_(pos), directStencilPointIndex_(directStencilPointIndex),
	interpStencilPointIndex_(interpStencilPointIndex) {}
    G2GPoint(void) {}

    const GridIndex& pos(void) const { return pos_; }
    const int directStencilPointIndex(void) const { return directStencilPointIndex_; }
    const int interpStencilPointIndex(void) const { return interpStencilPointIndex_; }

  private:
    // each point in the g2gUnion could be unique decide by      
    // its location, the origin is the center of the directStencil
    GridIndex pos_;
    // each interpStencil is centered at one directStencil point
    int directStencilPointIndex_; 
    // the local index within each interpStencil
    int interpStencilPointIndex_; 
    // pos is used to computed the kernel value. The other two items are 
    // used to build matrix "interpPointToG2GRowIndexMap"
  };

  bool operator < (const G2GPoint& p1, const G2GPoint& p2);
  bool operator != (const G2GPoint& p1, const G2GPoint& p2);
  bool operator == (const G2GPoint& p1, const G2GPoint& p2);

} //namespace pfft

#endif
