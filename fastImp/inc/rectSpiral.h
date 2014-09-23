/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: rectSpiral.h,v 1.1 2002/08/26 01:55:22 bsong Exp $";

  ==========================================================================
*/
#ifndef __RECT_SPIRAL_H_
#define __RECT_SPIRAL_H_

#include <iostream>
#include <string>
#include <vector>
#include "vector3D.h"
#include "condInfo.h"

namespace surf {

  class RectSpiral : public CondInfo {
  public:
    RectSpiral(void) {}
    RectSpiral(const CondType condType, const double unit)
      : CondInfo(condType, unit) {}
    virtual void readInputFile(void);

    void getTurns(size_t &f, size_t &h, size_t &q) const {
      f = fullTurns_;
      h = halfTurns_;
      q = quadTurns_;
    }

    size_t getNumPanelAlongWidth(void) const {return numPanelAlongWidth_;}
    size_t getNumPanelAlongThickness(void) const {return numPanelAlongThickness_;}
    size_t getNumPanelAlongShortestLength(void) const {return numPanelAlongShortestLength_;}

    double getWidth(void) const {return width_;}
    double getThickness(void) const {return thickness_;}
    double getSpacing(void) const {return spacing_;}
    double getStartLineLength(void) const {return startLineLength_;}

    pfft::point3D getHead(void) const {return head_;}
    pfft::point3D getPlaneNormal(void) const {return planeNormal_;}
    pfft::point3D getStartLineUnitY(void) const {return startLineUnitY_;}
    
  private:
    
    // -- geometric parameters --
    size_t fullTurns_;
    size_t halfTurns_;
    size_t quadTurns_;

    double width_;
    double spacing_;
    double thickness_;
    double startLineLength_;

    pfft::point3D head_;
    pfft::point3D startLineUnitY_;
    pfft::point3D planeNormal_;

    // -- discretization parameters --
    size_t numPanelAlongWidth_;
    size_t numPanelAlongThickness_;
    size_t numPanelAlongShortestLength_; // this is for the shortest one.
                                 // you need more for the larger one.
  };
}
    
#endif
