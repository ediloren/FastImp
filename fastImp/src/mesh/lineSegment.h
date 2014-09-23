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

  const static char cvsid[] = "$Id: lineSegment.h,v 1.2 2002/08/27 15:06:00 bsong Exp $";

  ==========================================================================
*/
#ifndef __LINE_SEGMENT_H_
#define __LINE_SEGMENT_H_

#include <iostream>
#include <vector>
#include "wirePath.h"

namespace mesh {
  
  enum LineSegmentType {
    HEAD = 0,
    MIDDLE = 1,
    TAIL = 2,
  };

  class LineSegment {

    friend class LineSegmentList;

  public:
    
    // -- constructor and destructor -- 
    LineSegment(void) {}
    ~LineSegment() {}

    // --  interface funcs  --
    double getLength(void) const { return length_; }
    double getWidth(void) const { return width_; }
    double getThickness(void) const { return thickness_; }

    bool IsFirstSegment (void) const { return type_ == HEAD; }
    bool IsLastSegment (void) const { return type_ == TAIL; }
    
    void global2Local (pfft::point3D& pnt);
    void local2Global (pfft::point3D& pnt);

    const pfft::point3D& getLocalUnitX(void) const {
      return localUnitX_;
    }

    const pfft::point3D& getLocalUnitY(void) const {
      return localUnitY_;
    }
    
    const pfft::point3D& getLocalUnitZ(void) const {
      return localUnitZ_;
    }

    const pfft::point3D& getStartVertex(void) const {
      return startVertex_;
    } 

    const pfft::point3D& getEndVertex(void) const {
      return endVertex_;
    } 

    // -- mesh generator -- 
    void genPanel (void);
    void genNode (void);

  private:

    size_t index_;

    size_t startVtxIdx_; // serve as local origin.
    size_t endVtxIdx_;

    pfft::point3D startVertex_;
    pfft::point3D endVertex_;

    pfft::point3D localUnitX_;
    pfft::point3D localUnitY_;
    pfft::point3D localUnitZ_;
    
    double length_;
    double width_;
    double thickness_;
    double x1, x2, y1, y2, z1, z2; // all in local coord.

    LineSegmentType type_;

  };

  class LineSegmentList {

  public:
    // -- constructor --
    LineSegmentList(void) {}

    // -- setup -- 
    void setup(double w, double t, const WirePath& wirePath);

    // -- interfaces -- 
    bool IsEmpty(void) const {return lineSegArray_.empty(); }
    size_t getNumSegment (void) const { return lineSegArray_.size(); }
    double getWidth(void) const {return width_; }
    double getThickness(void) const {return thickness_; }
    double getLength(size_t idx) const {return lineSegArray_[idx].getLength(); }

    const pfft::point3D& getLocalUnitX(const size_t idx) const {
      return lineSegArray_[idx].getLocalUnitX();
    }
    
    const pfft::point3D& getLocalUnitY(const size_t idx) const {
      return lineSegArray_[idx].getLocalUnitY();
    }

    const pfft::point3D& getLocalUnitZ(const size_t idx) const {
      return lineSegArray_[idx].getLocalUnitZ();
    }

    const pfft::point3D& getStartVertex(const size_t idx) const {
      return lineSegArray_[idx].getStartVertex();
    }
    
    const pfft::point3D& getEndVertex(const size_t idx) const {
      return lineSegArray_[idx].getEndVertex();
    }

    // -- generate mesh --
    void genMesh (void);
 
    void prnInfo (void);
    
  private: 

    double width_;
    double thickness_;
    
    std::vector<LineSegment> lineSegArray_;

  };
      

}

#endif
