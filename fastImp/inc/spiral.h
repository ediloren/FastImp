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

  const static char cvsid[] = "$Id: spiral.h,v 1.3 2002/07/18 15:05:44 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _SPIRAL_H_
#define _SPIRAL_H_

#include "vector3D.h"
#include "condInfo.h"

namespace surf {

  class Spiral : public CondInfo {


  public:

    Spiral(void) {}
    Spiral(CondType condType, const double unit) 
      : CondInfo(condType, unit) { }
    virtual void readInputFile(void);

    pfft::point3D origin1(void) const { return origin1_; }
    pfft::point3D origin2(void) const { return origin2_; }
    pfft::point3D edgePoint(void) const { return edgePoint_; }
    double outerRadius(void) const { return outerRadius_; }
    double innerRadius(void) const { return innerRadius_; }
    double r3(void) const { return r3_; }
    double numRound(void) const { return numRound_; }
    int numPanelRadius(void) const { return numPanelRadius_; }
    int numPanelArc(void ) const { return numPanelArc_; }
    int numPanelThickness(void) const { return numPanelThickness_; }

  private:

    pfft::point3D origin1_;
    pfft::point3D origin2_;
    pfft::point3D edgePoint_;
    double outerRadius_;
    double innerRadius_;
    double r3_;
    double numRound_;
    int numPanelRadius_;
    int numPanelArc_;
    int numPanelThickness_;
    
  };

} // namespace surf

#endif

