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

  const static char cvsid[] = "$Id: otherShape.h,v 1.2 2002/07/18 15:05:44 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _OTHER_SHAPE_H_
#define _OTHER_SHAPE_H_

#include "condInfo.h"

namespace surf {

  class OtherShape : public CondInfo {

  public:

    OtherShape(void) {}
    OtherShape(CondType condType, const double unit) 
      : CondInfo(condType, unit) { }
    virtual void readInputFile(void);

  };
    
} // namespace surf

#endif
