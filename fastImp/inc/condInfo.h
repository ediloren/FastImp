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

  const static char cvsid[] = "$Id: condInfo.h,v 1.9 2003/07/16 15:32:33 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _COND_INFO_H_
#define _COND_INFO_H_

#include "vector3D.h"
#include <vector>

namespace surf {

  class CondInfo {

  public:

    enum CondType { WIRE=1, RING, SPIRAL, GROUND,
		    OTHER_SHAPE, OTHER_SHAPE_GROUND, RECTSPIRAL};

    CondInfo(void) {}
    CondInfo(const CondType condType, const double unit) 
      : condType_(condType), unit_(unit) {}

    const double unit(void) const { return unit_; }
    const CondType condType(void) const { return condType_; }
    const double conductivity(void) const { 
#ifdef DISABLE_SCALING
      return conductivity_; // without scaling
#else
      return conductivity_ * unit_; // with scaling
#endif
    } 
    const double maxPanelSize(void) const { return maxPanelSize_; }

    virtual void readInputFile() = 0;

  protected:

    CondType condType_;
    double unit_;
    double conductivity_;
    double maxPanelSize_;

  };

  void readCondInfo(char* inputFileName, std::vector<CondInfo*>& condInfoList);
  double readSize(const double unit);
  pfft::point3D readPoint(const double unit);

} // namespace surf

#endif
