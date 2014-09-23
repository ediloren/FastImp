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

  const static char cvsid[] = "$Id: corner.cc,v 1.1 2002/08/26 01:55:22 bsong Exp $";

  ==========================================================================
*/
#include "corner.h"
#include "lineSegment.h"
#include <fstream>

using namespace std;
using namespace pfft;
using namespace mesh;

/**********************************************************************
 * CornerList::setup --
 **********************************************************************/
void CornerList::setup (
			double w,
			double t,
			const LineSegmentList& lineSegmentList)
{
  width_ = w;
  thickness_ = t;
  size_t numCorner = lineSegmentList.getNumSegment() - 1;
  if (numCorner != 0) {
    cornerArray_.resize(numCorner);
    
    for (size_t i = 0; i < numCorner; i ++) {
      cornerArray_[i].index_ = i;
      
      cornerArray_[i].leftSegmentIdx_ = i;
      cornerArray_[i].backSegmentIdx_ = i+1;
      
      cornerArray_[i].width_  = lineSegmentList.getWidth();
      cornerArray_[i].length_ = cornerArray_[i].width_;
      cornerArray_[i].thickness_ = lineSegmentList.getThickness();
      
      cornerArray_[i].leftSegCentroid_
	= lineSegmentList.getEndVertex(
				       cornerArray_[i].leftSegmentIdx_);
      
      cornerArray_[i].backSegCentroid_
	= lineSegmentList.getStartVertex(
					 cornerArray_[i].backSegmentIdx_);
      
      cornerArray_[i].localUnitX_
	= lineSegmentList.getLocalUnitX(
					cornerArray_[i].leftSegmentIdx_);
      
      cornerArray_[i].localUnitY_
	= lineSegmentList.getLocalUnitY(
					cornerArray_[i].leftSegmentIdx_);
      
      cornerArray_[i].localUnitZ_ 
	= crossProd(
		    cornerArray_[i].localUnitX_,
		    cornerArray_[i].localUnitY_);
      
      cornerArray_[i].x1 = -.5 * cornerArray_[i].width_;
      cornerArray_[i].x2 = +.5 * cornerArray_[i].width_;
      cornerArray_[i].z1 = -.5 * cornerArray_[i].thickness_;
      cornerArray_[i].z2 = +.5 * cornerArray_[i].thickness_;
      cornerArray_[i].y1 = 0.;
      cornerArray_[i].y2 = cornerArray_[i].length_;
    }
  }
}

/**********************************************************************
 * CornerList::prnInfo --
 **********************************************************************/
void CornerList::prnInfo (void)
{
  std::ofstream fout("cornerinfo.dat");
  fout << " <-----   corner info  -----> " << std::endl << std::endl;
  for (size_t i=0; i<cornerArray_.size(); i++) {
    fout << " line " << i << " : " << std::endl;
    fout << " index => " << cornerArray_[i].index_ << std::endl;
    fout << " left segment index => " << cornerArray_[i].leftSegmentIdx_ << std::endl;
    fout << " back segment index => " << cornerArray_[i].backSegmentIdx_ << std::endl;
    fout << " left seg centroid  => }" << std::endl;
    fout << cornerArray_[i].leftSegCentroid_ << "}" <<std::endl;
    fout << " back seg centroid => }" << std::endl;
    fout << cornerArray_[i].backSegCentroid_ << "}" <<std::endl;
    fout << " local unit X => {" << std::endl;
    fout << cornerArray_[i].localUnitX_ << "}" << std::endl;
    fout << " local unit Y => {" << std::endl;
    fout << cornerArray_[i].localUnitY_ << "}" << std::endl;
    fout << " local unit Z => {" << std::endl;
    fout << cornerArray_[i].localUnitZ_ << "}" << std::endl;   
    fout << " x1 = " << cornerArray_[i].x1 << " ; ";
    fout << " x2 = " << cornerArray_[i].x2 << " ; " << std::endl;
    fout << " y1 = " << cornerArray_[i].y1 << " ; ";
    fout << " y2 = " << cornerArray_[i].y2 << " ; " << std::endl;
    fout << " z1 = " << cornerArray_[i].z1 << " ; ";
    fout << " z2 = " << cornerArray_[i].z2 << " ; " << std::endl;
  }
  fout << std::endl << std::endl;
  fout.close();
}

