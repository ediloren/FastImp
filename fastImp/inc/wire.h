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

  const static char cvsid[] = "$Id: wire.h,v 1.3 2002/07/18 15:05:44 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _WIRE_H_
#define _WIRE_H_

#include "vector3D.h"
#include "condInfo.h"

namespace surf {

  class Wire : public CondInfo {

  public:

    Wire(void) {}
    Wire(const CondType condType, const double unit) 
      : CondInfo(condType, unit) { }
    virtual void readInputFile(void);

    pfft::point3D leftEndPoint0(void) const { return leftEndPoint0_; }
    pfft::point3D leftEndPoint1(void) const { return leftEndPoint1_; }
    pfft::point3D leftEndPoint2(void) const { return leftEndPoint2_; }
    pfft::point3D rightEndPoint(void) const { return rightEndPoint_; }
    int numPanelBetweenPoint01(void) const { return numPanelBetweenPoint01_; }
    int numPanelBetweenPoint02(void) const { return numPanelBetweenPoint02_; }
    int numPanelLength(void) const { return numPanelLength_; }
    double maxPanelSize(void) const { return maxPanelSize_; }

  private:

    pfft::point3D leftEndPoint0_;
    pfft::point3D leftEndPoint1_;
    pfft::point3D leftEndPoint2_;
    pfft::point3D rightEndPoint_;
    int numPanelBetweenPoint01_;
    int numPanelBetweenPoint02_;
    int numPanelLength_;

  };
    
} // namespace surf

#endif
