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

  const static char cvsid[] = "$Id: stencil.h,v 1.7 2003/02/11 03:06:25 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _STENCIL_H_
#define _STENCIL_H_

#include <vector>
#include "gridIndex.h"
#include "vector3D.h"
#include "element.h"
#include "cubature.h"
#include "kernelIntegration.h"
#include "cmat.h" // for TNT::Matrix

namespace pfft {

  class Grid;
  enum StencilShape { CUBE_STENCIL, SPHERE_STENCIL, COMPACT_STENCIL, OTHER_STENCIL_SHAPE };
  enum StencilType { DIRECT_STENCIL, PROJECT_STENCIL, 
		     INTERP_STENCIL, OTHER_STENCIL_TYPE };
  enum BasisType { CONSISTENT_POLY, CUBED_POLY, SPECIAL_POLY, OTHER_POLY, NO_POLY };

  class Stencil {

  public:

    Stencil(const StencilType type, const int size);
    Stencil(void) {};

    void setup(const double step);
    bool isInsideStencil(const int x, const int y, const int z) const;
    bool isInsideStencil(const GridIndex& index3D) const;
    StencilShape shape(void) const { return stencilShape_; }
    StencilType type(void) const { return stencilType_; }
    int size(void) const { return size_; }
    double step(void) const { return step_; }
    int numBasis(void) const { return numBasis_; }
    int numPoint(void) const { return totalNumPoint_; }
    bool setupComplete(void) const;
    void findGlobalIndexOfStencilPoint(const int, const int, const int, 
				       const Grid&, std::vector<int>&) const;
    void findGlobal3DIndexOfStencilPoint(const GridIndex& hostPoint, 
					 const Grid&, std::vector<GridIndex>&) const;
    void evaluateBasisFuncAtOnePoint(point3D&, const vector3D<double>&, 
				     TNT::Vector<double>&, 
				     const DifferentialOperator) const;
    void evaluateBasisFuncAtOnePoint(double, double, double, 
				     std::vector<double>&) const;
    void integrateBasisFuncOverElement(const element&, const point3D& origin, 
				       const vector3D<double>&, 
				       TNT::Vector<double>&, 
				       const DifferentialOperator); 
    void compCoef(const TNT::Vector<double>&, TNT::Vector<double>&) const;
    const GridIndex point(size_t pointIndex) const { return points_[pointIndex]; }
    const GridIndex point(int pointIndex) const { return points_[pointIndex]; }
    std::vector<GridIndex> points(void) const { return points_; }
    int transfer3DIndexTo1D (const GridIndex& gridIndex) const {
      if(stencilShape_ == CUBE_STENCIL) {
	return (gridIndex.z() + size_) + (2*size_+1)*((gridIndex.y()+size_) + 
		(2*size_+1) * (gridIndex.x() + size_));
      } else if (stencilShape_ == SPHERE_STENCIL){
	return local1DIndexArrayForSphereStencil_[gridIndex.x() + 
                size_][gridIndex.y() + size_][gridIndex.z() + size_];
      } else if (stencilShape_ == COMPACT_STENCIL) {
	return local1DIndexArrayForCompactStencil_[gridIndex.x() + 
                size_][gridIndex.y() + size_][gridIndex.z() + size_];
      } else {
	std::cerr<<" sorry! we only recommend to use cube or sphere stencil "
		 << std::endl;
	exit(1);
	return 0;
      }
    }
  
  private:

    // essential info
    StencilShape stencilShape_;
    StencilType stencilType_;
    // number of grid points along each direction, it is 1 for 3*3*3 cube stencil.
    int size_; 
    // compact stencil is union of many small spheres. This is the radius of 
    // the largest sphere
    int compactStencilSize_;
    // distance between adjacent points. It is assumed here that only 
    // uniform grid is used.
    double step_; 
    // The type of polynomial used for interploation or projection
    BasisType basisType_;

    // derived info
    int totalNumPoint_; // total number of points in this stencil
    int numBasis_; // number of basis functions (poly terms)
    std::vector<GridIndex> points_; // full list of stencil points, in terms of (i,j,k)

    // value of each basis function at each point
    TNT::Matrix<double> basisValue_; 
    // interpolation weight. It is the pseudo-inverse of basisValue
    TNT::Matrix<double> weight_; 
    bool equaScaleFlag_;
    std::vector<double> equaScale_;
    double polyScale_;
    std::vector<int> xMonoOrder_, yMonoOrder_, zMonoOrder_;
    int maxMonoOrder_;
    int maxPolyOrder_;
    Cubature cubature_;

    std::vector<TNT::Matrix<int> > local1DIndexArrayForSphereStencil_;
    std::vector<TNT::Matrix<int> > local1DIndexArrayForCompactStencil_;

    void findStencilShape(void);
    void findBasisType(void);
    void setupStencilPoint(void);
    void findSphereStencilPoints(void);
    void findCubeStencilPoints (void);
    void findCompactStencilPoints (void);
    void findNumBasisFunc(void);
    void setupEquaScale(void);
    void setupPolyOrder(void);
    void setupBasisValue(void);
    void setupWeight(void);
    void addIntegralOverTriangle (const std::vector<point3D>&, 
				  const vector3D<double>&, 
				  TNT::Vector<double>&, 
				  const DifferentialOperator) const;
    inline double compTriangleArea (const point3D&, const point3D&, 
				    const point3D&) const;
    inline double basisFuncValue (const point3D&, const int polyTermIndex) const;
    inline double basisFuncValue_dx (const point3D&, const int polyTermIndex) const;
    inline double basisFuncValue_dy (const point3D&, const int polyTermIndex) const;
    inline double basisFuncValue_dz (const point3D&, const int polyTermIndex) const;
    // Enrico, corrected to avoid gcc error
    inline double intpow(const double d, int e) const;
    //inline double Stencil::intpow(const double d, int e) const;

  };

  /*
   * The stencil type determines the shape of the stencil:
   * 1: "Cube." The number of points is (2n+1)^3 where n is integer.
   *    Make a cube with side length 2*stencilSize*"gridsize" centered 
   *    on the node of interest. All grid points inside the cube are
   *    neighbours. A small number will be added to stencilSize
   *    such that points on the boundary are counted as "inside".
   *
   *    Example (size 1 - 27 points in 3D)
   *       x: normal grid point 
   *       o: stencil point
   *
   *                  |   |   |   |   | 
   *                - x - x - x - x - x -
   *                  |   |   |   |   | 
   *                - x - o - o - o - x -
   *                  |   |   |   |   | 
   *                - x - o - o - o - x -
   *                  |   |   |   |   | 
   *                - x - o - o - o - x -
   *                  |   |   |   |   | 
   *                - x - x - x - x - x -
   *                  |   |   |   |   | 
   *
   * 2: "Sphere." Make a sphere of radius stencilSize*dx centered on 
   *    the node of interest. All grid points inside the sphere are
   *    neighbours. A small number will be added to stencilSize
   *    such that points on the boundary are counted as "inside".
   *    Note that the "sphere" never contains more stencil points 
   *    than does the cube of the same size.
   *    Example (size 1 - 7 points in 3D)
   *       x: normal grid point 
   *       o: stencil point
   *
   *                  |   |   |   |   | 
   *                - x - x - x - x - x -
   *                  |   |   |   |   | 
   *                - x - x - o - x - x -
   *                  |   |   |   |   | 
   *                - x - o - o - o - x -
   *                  |   |   |   |   | 
   *                - x - x - o - x - x -
   *                  |   |   |   |   | 
   *                - x - x - x - x - x -
   *                  |   |   |   |   | 
   *
   */

} //namespace pfft

#endif
