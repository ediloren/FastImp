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

  const static char cvsid[] = "$Id: corner.h,v 1.1 2002/08/26 01:55:22 bsong Exp $";

  ==========================================================================
*/
#ifndef __CORNER_H_
#define __CORNER_H_

#include <iostream>
#include <vector>
#include "vector3D.h"
#include "lineSegment.h"

namespace mesh {

  class Corner {
    friend class CornerList;
    
  private:
    size_t index_;

    size_t leftSegmentIdx_;
    size_t backSegmentIdx_;

    pfft::point3D leftSegCentroid_;
    pfft::point3D backSegCentroid_;
    
    pfft::point3D localUnitX_;
    pfft::point3D localUnitY_;
    pfft::point3D localUnitZ_;

    double length_;
    double width_;
    double thickness_;
    double x1, x2, y1, y2, z1, z2;

  public:
    // -- constructor and destructor --
    Corner(void) {};
    ~Corner() {};

    // -- interface funcs --
    double getLength(void) {return length_; }
    double getWidth(void) { return width_; }
    double getThickness(void) { return thickness_; }

    const pfft::point3D& getLeftSegCentroid(void) const {
      return leftSegCentroid_;
    }

    const pfft::point3D& getBackSegCentroid(void) const {
      return backSegCentroid_;
    }

    const pfft::point3D& getLocalUnitX(void) const {
      return localUnitX_;
    }

    const pfft::point3D& getLocalUnitY(void) const {
      return localUnitY_;
    }

    const pfft::point3D& getLocalUnitZ(void) const {
      return localUnitZ_;
    }
    
    // -- mesh generator -- 
    void genPanel (void);
    void genNode (void);
  };

  class CornerList {
    
  private:
    double width_;
    double thickness_;
    std::vector<Corner> cornerArray_;
    
  public:
    // -- constructor --
    CornerList(void) {}
    
    // -- setup --
    void setup(double, double, const LineSegmentList&);
    
    // -- interface -- 
    bool IsEmpty(void) const {return cornerArray_.empty(); }
    size_t getNumCorner(void) const {return cornerArray_.size(); }
    
    double getWidth(void) const {return width_;}
    double getThickness(void) const {return thickness_;}
    double getLength(void) const {return width_;}

    const pfft::point3D getLeftSegCentroid (const size_t idx) const {
      return cornerArray_[idx].getLeftSegCentroid();
    }

    const pfft::point3D getBackSegCentroid (const size_t idx) const {
      return cornerArray_[idx].getBackSegCentroid();
    }

    const pfft::point3D getLocalUnitX (const size_t idx) const {
      return cornerArray_[idx].getLocalUnitX();
    }
    
    const pfft::point3D getLocalUnitY (const size_t idx) const {
      return cornerArray_[idx].getLocalUnitY();
    }

    const pfft::point3D getLocalUnitZ (const size_t idx) const {
      return cornerArray_[idx].getLocalUnitZ();
    }

    // -- generate mesh --
    void genMesh(void);
    
    // -- I/O functions --
    void prnInfo(void);
  };


}

#endif
