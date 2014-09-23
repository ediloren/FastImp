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

  const static char cvsid[] = "$Id: wirePath.h,v 1.2 2002/08/27 15:06:00 bsong Exp $";

  ==========================================================================
*/
#ifndef __WIRE_PATH_H_
#define __WIRE_PATH_H_

#include <iostream>
#include <exception>
#include <vector>
#include "vector3D.h"

namespace mesh{
  
  class WirePath {

  public:
   
    // constructors
    WirePath(void) {}
    
    // -- setup function --
    void setup(
	       double width,
	       double thickness,
	       const std::vector<pfft::point3D>& vertexArray,
	       const std::vector<pfft::point3D>& localUnitXArray,
	       const std::vector<pfft::point3D>& localUnitYArray,
	       const std::vector<double>& lengthArray) {
      this->width_ = width;
      this->thickness_ = thickness;
      this->vtxArray_ = vertexArray;
      this->localUnitXArray_ = localUnitXArray;
      this->localUnitYArray_ = localUnitYArray;
      this->lengthArray_ = lengthArray;
    }
							     
    // -- interface funcs --
    double getWidth (void) const { return width_; }
    double getThickness (void) const { return thickness_; }
    double getLength (size_t idx) const { return lengthArray_[idx]; }

    //get vertices related information
    bool IsEmpty (void) const { return vtxArray_.empty(); }
    size_t getNumVertex (void) const { return vtxArray_.size(); }
    size_t getNumSegment (void) const { return lengthArray_.size(); }
    const pfft::point3D& getVertex  (size_t vtxidx) const { 
      return vtxArray_[vtxidx]; // add exception handling later.
    } 
    const pfft::point3D& getLocalUnitX(size_t idx) const {
      return localUnitXArray_[idx];
    }
    const pfft::point3D& getLocalUnitY(size_t idx) const {
      return localUnitYArray_[idx];
    }
    // -- I/O function -- 
    // finally, I'm gonna implement it in overloaded << >>
    void prnInfo(void);

  private:

    // -- geometric information --
    // width -- the size along the local x direction.
    // thickness -- the size along the local z direction.
    // vtxArray -- the center of each line segement.
    // localUnitXArray -- the local unix x axis of each segment.
    // localUnitYArray -- the longitudinal direction of each seg.
    // lengthArray -- the length of each segment.
    double width_; 
    double thickness_; 

    std::vector<pfft::point3D> vtxArray_; 
    std::vector<pfft::point3D> localUnitXArray_;
    std::vector<pfft::point3D> localUnitYArray_;
    
    std::vector<double> lengthArray_;

  };
}    
      
#endif
