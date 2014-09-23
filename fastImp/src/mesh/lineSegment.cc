/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 * All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: lineSegment.cc,v 1.1 2002/08/26 01:55:22 bsong Exp $";

  ==========================================================================
*/
#include "lineSegment.h"
#include <fstream>

using namespace mesh;
using namespace pfft;
using namespace std;

/**********************************************************************
 * LineSegmentList::setup --
 **********************************************************************/
void LineSegmentList::setup (
			     double w,
			     double t,
			     const WirePath& wirePath)
{
  width_ = w;
  thickness_ = t;

  size_t numLineSegment = wirePath.getNumSegment();
  lineSegArray_.resize(numLineSegment);

  for (size_t i = 0; i < numLineSegment; i++) {
    lineSegArray_[i].index_ = i;

    lineSegArray_[i].startVtxIdx_ = i;
    lineSegArray_[i].endVtxIdx_ = i+1;
    
    lineSegArray_[i].width_ = w;
    lineSegArray_[i].thickness_ = t;
    
    lineSegArray_[i].localUnitX_ = wirePath.getLocalUnitX(i);
    lineSegArray_[i].localUnitY_ = wirePath.getLocalUnitY(i);
    lineSegArray_[i].localUnitZ_ 
      = crossProd(
		  lineSegArray_[i].localUnitX_,
		  lineSegArray_[i].localUnitY_);

    pfft::point3D pnt1
      = wirePath.getVertex(lineSegArray_[i].startVtxIdx_);
    pfft::point3D pnt2
      = wirePath.getVertex(lineSegArray_[i].endVtxIdx_);
  
    lineSegArray_[i].startVertex_ 
      = (i == 0) 
      ? pnt1
      : pnt1 + .5 * w * lineSegArray_[i].localUnitY_;

    lineSegArray_[i].endVertex_ 
      = (i == numLineSegment - 1) 
      ? pnt2
      : pnt2 - .5 * w * lineSegArray_[i].localUnitY_;

    lineSegArray_[i].length_ 
      = length(lineSegArray_[i].endVertex_ - lineSegArray_[i].startVertex_);

    lineSegArray_[i].type_
      = (i==0) ? HEAD 
      : ((i==numLineSegment-1) ? TAIL
	 : MIDDLE);

    lineSegArray_[i].x1 = -.5 * w;
    lineSegArray_[i].x2 = +.5 * w;
    lineSegArray_[i].z1 = -.5 * t;
    lineSegArray_[i].z2 = +.5 * t;
    lineSegArray_[i].y1 = 0.;
    lineSegArray_[i].y2 = lineSegArray_[i].length_;
  }
}

/**********************************************************************
 * LineSegmentList::prnInfo --
 **********************************************************************/
void LineSegmentList::prnInfo (void)
{
  std::ofstream fout("seginfo.dat");
  fout << " <-----  line segment info  -----> " << std::endl << std::endl;
  for (size_t i=0; i<lineSegArray_.size(); i++) {
    fout << " line " << i << " : " << std::endl;
    fout << " index => " << lineSegArray_[i].index_ << std::endl;
    fout << " start vertex index => " << lineSegArray_[i].startVtxIdx_ << std::endl;
    fout << " end vertex index => " << lineSegArray_[i].endVtxIdx_ << std::endl;
    fout << " start vertex  => }" << std::endl;
    fout << lineSegArray_[i].startVertex_ << "}" <<std::endl;
    fout << " end vertex => }" << std::endl;
    fout << lineSegArray_[i].endVertex_ << "}" <<std::endl;
    fout << " local unit X => {" << std::endl;
    fout << lineSegArray_[i].localUnitX_ << "}" << std::endl;
    fout << " local unit Y => {" << std::endl;
    fout << lineSegArray_[i].localUnitY_ << "}" << std::endl;
    fout << " local unit Z => {" << std::endl;
    fout << lineSegArray_[i].localUnitZ_ << "}" << std::endl;   
    fout << " x1 = " << lineSegArray_[i].x1 << " ; ";
    fout << " x2 = " << lineSegArray_[i].x2 << " ; " << std::endl;
    fout << " y1 = " << lineSegArray_[i].y1 << " ; ";
    fout << " y2 = " << lineSegArray_[i].y2 << " ; " << std::endl;
    fout << " z1 = " << lineSegArray_[i].z1 << " ; ";
    fout << " z2 = " << lineSegArray_[i].z2 << " ; " << std::endl;
  }
  fout << std::endl << std::endl;
  fout.close();
}

/**********************************************************************
 * LineSegment::genMesh --
 **********************************************************************/
/*
void LineSegment::genMesh (
			   size_t lnSegIdx,
			   std::vector<nodeLis)
{

}
*/

/**********************************************************************
 * LineSegment::global2Local --
 * from global coordinate system to line segment local coordinate
 * system. what's line segment local coordinate system ?
 * axes: local unix {X,Y,Z}; 
 * origin: wirePath.getVertex(startVtxIdx_);
 **********************************************************************/
/*
void LineSegment::global2Local (pfft::point3D& pnt)
{
  pfft::point3D offset = pnt - this->startVertex_;
  pnt = pfft::point3D(
		      dotProd(offset, this->localUnitX_),
		      dotProd(offset, this->localUnitY_),
		      dotProd(offset, this->localUnitZ_));
}
*/

/**********************************************************************
 * LineSegment::local2Global --
 * from line segment local coordinate system to global coord sys.
 **********************************************************************/
/*
void LineSegment::local2Global (pfft::point3D& pnt)
{
  pfft::point3D offset
    = pnt.x() * this->localUnitX_ 
    + pnt.y() * this->localUnitY_
    + pnt.z() * this->localUnitZ_;
  
  pnt = this->startVertex_ + offset;
}
*/

