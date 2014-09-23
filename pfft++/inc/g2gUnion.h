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

  const static char cvsid[] = "$Id: g2gUnion.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _G2G_UNION_H_
#define _G2G_UNION_H_

#include "cmat.h"
#include "gridIndex.h"
#include "stencil.h"
#include "g2gPoint.h"
#include <vector>
#include <algorithm>

namespace pfft {

  template <class GreenFunc, class KernelValueType>

  class G2GUnion {

  public:

    G2GUnion (const Stencil&, const Stencil&, const Stencil&, const GreenFunc&);
    G2GUnion(void) {};

    const TNT::Matrix<KernelValueType>& matrix(void) const { return g2gMat; }
    const std::vector<int>& pointToG2GRowMap(int directStencilPointIndex) const { 
      return pointToG2GRowIndexMap[directStencilPointIndex]; }

    void computeG2GMatrix (const double, const Stencil&, const Stencil&, 
			   const GreenFunc&, const GridIndex&, const GridIndex&,
			   TNT::Matrix<KernelValueType>&);

  private:
    std::vector<std::vector<int> > pointToG2GRowIndexMap;
    TNT::Matrix<KernelValueType> g2gMat;
    Stencil projectStencil;
    Stencil interpStencil; 
    Stencil directStencil;
    GreenFunc greenFunc;
    std::vector<G2GPoint> allUnionPoint;

    void findAllUnionPoint(void);
    void findPointToG2GRowMap(void);
    void fillG2GMat(void);
    size_t countNumPoint (void);
    void pointReport(void);

  };

  /**********************************************************************
   * Grid2Grid --
   **********************************************************************/
  template <class GreenFunc, class KernelValueType>
  G2GUnion<GreenFunc, KernelValueType>::G2GUnion (
						  const Stencil& projectStencilIn, 
						  const Stencil& interpStencilIn, 
						  const Stencil& directStencilIn,
						  const GreenFunc& greenFuncIn)
    : projectStencil(projectStencilIn), interpStencil(interpStencilIn)
    , directStencil(directStencilIn), greenFunc(greenFuncIn)
  {
    findAllUnionPoint();
    findPointToG2GRowMap();
    fillG2GMat();
  }

  /**********************************************************************
   * computeG2GMatrix --
   * Given the center of the projectStencil (hostPoint) and the center
   * of the interpStencil (guestPoint), compute the small g2g matrix
   * that contains the interaction between these two stencils
   * This function is used for irregular neighbors.
   *
   * Note: hostPoint and guestPoint are all in global coordinate system
   **********************************************************************/
  template <class GreenFunc, class KernelValueType>
  void
  computeG2GMatrix (
		    const double gridStep,
		    const Stencil& projectStencil,
		    const Stencil& interpStencil,
		    GreenFunc greenFunc,
		    const GridIndex& guestPoint,
		    const GridIndex& hostPoint, 
		    TNT::Matrix<KernelValueType>& g2g)
  {
    for (int interpStencilPointIndex = 0; 
	 interpStencilPointIndex < interpStencil.numPoint(); 
	 interpStencilPointIndex++) {
      // shift to the global coordinate system
      GridIndex interpPoint = interpStencil.point(interpStencilPointIndex) + guestPoint;

      for (int projectStencilPointIndex = 0; 
	   projectStencilPointIndex < projectStencil.numPoint(); 
	   projectStencilPointIndex++) {
	// shift to the global coordinate system
	GridIndex projectPoint = projectStencil.point(projectStencilPointIndex) + hostPoint;

	GridIndex pos = projectPoint - interpPoint;
	vector3D<double> coord(pos.x(), pos.y(), pos.z());
	g2g[interpStencilPointIndex][projectStencilPointIndex] = greenFunc(coord * gridStep);
      }
    }
  }

  /**********************************************************************
   * findAllUnionPoint --
   **********************************************************************/
  template <class GreenFunc, class KernelValueType>
  void
  G2GUnion<GreenFunc, KernelValueType>::findAllUnionPoint (
							   void)
  {
    for (int di = 0; di < directStencil.numPoint(); di++) {
      const GridIndex& directStencilPointPos = directStencil.point(di);
      for (int ii = 0; ii < interpStencil.numPoint(); ii++) {
	const GridIndex& interpStencilPointPos = interpStencil.point(ii);
	GridIndex pos = directStencilPointPos + interpStencilPointPos;
	allUnionPoint.push_back(G2GPoint(pos, di, ii));
      }
    }
    stable_sort(allUnionPoint.begin(), allUnionPoint.end());

  }

  /**********************************************************************
   * findPointToG2GRowMap --
   **********************************************************************/
  template <class GreenFunc, class KernelValueType>
  void
  G2GUnion<GreenFunc, KernelValueType>::findPointToG2GRowMap (
							      void)
  {
    pointToG2GRowIndexMap.resize(directStencil.numPoint());
    for (int row = 0; row < directStencil.numPoint(); row++) {
      pointToG2GRowIndexMap[row].resize(interpStencil.numPoint());
    }

    size_t sequenceIndex = 0;
    typename std::vector<G2GPoint>::size_type index;
    for (index = 0; index < allUnionPoint.size()-1; index++) {
      int rowIndex = allUnionPoint[index].directStencilPointIndex();
      int colIndex = allUnionPoint[index].interpStencilPointIndex();
      pointToG2GRowIndexMap[rowIndex][colIndex] = sequenceIndex;
      if ( allUnionPoint[index] != allUnionPoint[index+1] ) {
	sequenceIndex++;
      }
    }
    // I have to do this to prevent over-the-bound-read
    int rowIndex = allUnionPoint[index].directStencilPointIndex();
    int colIndex = allUnionPoint[index].interpStencilPointIndex();
    pointToG2GRowIndexMap[rowIndex][colIndex] = sequenceIndex;

  }

  /**********************************************************************
   * fillG2GMat --
   **********************************************************************/
  template <class GreenFunc, class KernelValueType>
  void
  G2GUnion<GreenFunc, KernelValueType>::fillG2GMat (
						    void)
  {
    allUnionPoint.erase(unique(allUnionPoint.begin(), allUnionPoint.end()), 
			allUnionPoint.end());
    //    pointReport();

    g2gMat = TNT::Matrix<KernelValueType>(allUnionPoint.size(), 
					  projectStencil.numPoint());
    int rowIndex = 0;
    for (std::vector<G2GPoint>::const_iterator interpPoint = allUnionPoint.begin();
	 interpPoint != allUnionPoint.end(); ++interpPoint, rowIndex++) {
      for (int projectPointIndex = 0; 
	   projectPointIndex < projectStencil.numPoint(); projectPointIndex++) {
	GridIndex projectPoint = projectStencil.point(projectPointIndex);
	double dx = (projectPoint.x() - interpPoint->pos().x()) * directStencil.step();
	double dy = (projectPoint.y() - interpPoint->pos().y()) * directStencil.step();
	double dz = (projectPoint.z() - interpPoint->pos().z()) * directStencil.step();
	g2gMat[rowIndex][projectPointIndex] = greenFunc(dx, dy, dz);
      }
    }	
  }     

  /**********************************************************************
   * countNumPoint --
   **********************************************************************/
  template <class GreenFunc, class KernelValueType>
  void 
  G2GUnion<GreenFunc, KernelValueType>::pointReport (
						     void)
  {
#ifdef DEBUG_ELEMENT
    std::cout << "total number of points in g2g union := " 
	      << allUnionPoint.size() 
	      << std::endl;
#endif
  }

  /**********************************************************************
   * countNumPoint --
   **********************************************************************/
  template <class GreenFunc, class KernelValueType>
  size_t 
  G2GUnion<GreenFunc, KernelValueType>::countNumPoint (
						       void)
  {
    std::vector<G2GPoint> pointList;
    
    for (int di = 0; di < directStencil.numPoint(); di++) {
      const GridIndex& directStencilPointPos = directStencil.point(di);
      for (int ii = 0; ii < interpStencil.numPoint(); ii++) {
	const GridIndex& interpStencilPointPos = interpStencil.point(ii);
	GridIndex pos = directStencilPointPos + interpStencilPointPos;
	pointList.push_back(G2GPoint(pos, di, ii));
      }
    }
    stable_sort(pointList.begin(), pointList.end());
    pointList.erase(unique(pointList.begin(), pointList.end()), 
		    pointList.end());

    return pointList.size();
  }

} //namespace pfft

#endif

