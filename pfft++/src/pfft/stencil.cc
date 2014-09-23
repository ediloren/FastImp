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

  ==========================================================================
*/

const static char cvsid[] = "$Id: stencil.cc,v 1.9 2003/02/11 03:13:48 zhzhu Exp $";

#include "dense.h" // for pseudo_inverse()
#include <cstdlib> // for new and delete
#include <stdexcept> //for exception handling
#include "stencil.h"
#include "grid.h"
#include "vector3D.h"
#include "utils.h" // for errorMessage()

using namespace std;
using namespace pfft;
using namespace TNT;

/**********************************************************************
 * stencil --
 **********************************************************************/
Stencil::Stencil(
		 const StencilType type,
		 const int size) 
  : stencilType_(type), size_(size)
{
  findBasisType();
  findStencilShape();
  setupStencilPoint();
  if ( (stencilType_ == PROJECT_STENCIL) || (stencilType_ == INTERP_STENCIL) ) {
    findNumBasisFunc();
  }
}

/**********************************************************************
 * findBasisType --  
 **********************************************************************/
void
Stencil::findBasisType (
			void)
{
  if ( (stencilType_ == PROJECT_STENCIL) || (stencilType_ == INTERP_STENCIL) ) {
#ifdef DISABLE_CONSISTENT_STENCIL
    basisType_ = CUBED_POLY;
#else
    switch (size_) {
    case 1:
      basisType_ = CUBED_POLY;
      break;
    case 2:
      basisType_ = CONSISTENT_POLY;
      compactStencilSize_ = 5;
      break;
    case 3:
      basisType_ = CUBED_POLY;
      // Empirical study on a sphere with 972 panels shows that 12 is the optimal value.
      // With this setting, the error decreases by 50% from size_=2. Not good enough.
      // I need to find the scaling factor, just like size_=2 case. Before this is done,
      // stick with cubed_poly because it guarantees good accuracy, with high cost though.
      //	basisType_ = CONSISTENT_POLY;
      //	compactStencilSize_ = 12; 
      break;
    default:
      errorMessage("stencil.cc::findBasisType()",
		   "The order of interpolation or projection is too high. Computational cost is prohibitive");
    }
#endif
  } else {
    basisType_ = NO_POLY;
  }
}

/**********************************************************************
 * findStencilShape --
 **********************************************************************/
void
Stencil::findStencilShape (
			   void)
{
  if ( (stencilType_ == PROJECT_STENCIL) || (stencilType_ == INTERP_STENCIL) ) {
    if (basisType_ == CUBED_POLY) {
      stencilShape_ = CUBE_STENCIL;
    } else if (basisType_ == CONSISTENT_POLY) {
      stencilShape_ = COMPACT_STENCIL;
    } else {
      errorMessage("Stencil::findStencilShape()", "unknown basis type!");
    }
  } else if (stencilType_ == DIRECT_STENCIL) {
    stencilShape_ = SPHERE_STENCIL;
  } else {
    stencilShape_ = OTHER_STENCIL_SHAPE;
  }
}

/**********************************************************************
 * findNumBasisFunc --
 **********************************************************************/
void
Stencil::findNumBasisFunc (
			   void)
{
  maxMonoOrder_ = 2*size_;
  switch (basisType_) {
  case CUBED_POLY:
    numBasis_ = totalNumPoint_;
    equaScaleFlag_ = false;
    break;
  case CONSISTENT_POLY:
    numBasis_ = (maxMonoOrder_ + 1) * (maxMonoOrder_ + 2) * (maxMonoOrder_ + 3) / 6;
    equaScaleFlag_ = true;
    break;
  case NO_POLY:
    break;
  default:
    errorMessage("Stencil::findNumBasisFunc()", "unknown basis type!");
    break;
  }
}

/**********************************************************************
 * setupStencilPoint --
 **********************************************************************/
void
Stencil::setupStencilPoint (
			    void)
{
  if (stencilShape_ == CUBE_STENCIL) {
    findCubeStencilPoints();
  } else if (stencilShape_ == SPHERE_STENCIL) {
    findSphereStencilPoints();
  } else if (stencilShape_ == COMPACT_STENCIL) {
    findCompactStencilPoints();
  } else {
    errorMessage("Stencil::setupStencilPoint()", "unknown stencil shape!");
  }
  totalNumPoint_ = points_.size();
}

/**********************************************************************
 * findSphereStencilPoints --
 * Note:  the loop sequence in this function should be consistent with
 * that in findGlobalIndexOfSphereStencilPoint
 **********************************************************************/
void
Stencil::findSphereStencilPoints (
				  void)
{
  int min = -size_;
  int max = size_;

  int NN = 2*size_ + 1;
  local1DIndexArrayForSphereStencil_ = 
    vector<TNT::Matrix<int> >(NN, Matrix<int>(NN, NN));
  
  points_.clear();

  int currentLocal1DIndex = 0;
  for (int x=min; x<=max; x++) {
    for (int y=min; y<=max; y++) {
      for (int z=min; z<=max; z++) {
	GridIndex grid(x, y, z);
	if (grid.length2() <= size_*size_) {
	  points_.push_back(grid);
	  local1DIndexArrayForSphereStencil_[x+size_][y+size_][z+size_]
	    = currentLocal1DIndex ++;
	} else {
	  local1DIndexArrayForSphereStencil_[x+size_][y+size_][z+size_]
	    = -1;
	}
      }
    }
  }
}

/**********************************************************************
 * findCubeStencilPoints --
 * Note:  the loop sequence in this function should be consistent with
 * that in findGlobalIndexOfCubeStencilPoint
 **********************************************************************/
void
Stencil::findCubeStencilPoints (
				void)
{
  int min = -size_;
  int max = size_;

  points_.clear();
  for (int x=min; x<=max; x++) {
    for (int y=min; y<=max; y++) {
      for (int z=min; z<=max; z++) {
	GridIndex grid(x, y, z);
	points_.push_back(grid);
      }
    }
  }
}

/**********************************************************************
 * findCompactStencilPoints --
 * This is only for 57 points in the s_{0,1,2,3,4,5} type stencil
 * Stencil s_i means all points whose grid.length2() == i.
 * So it is a union of many sphere stencils with radius being 0,1,2,3,4,5
 **********************************************************************/
void
Stencil::findCompactStencilPoints (
				   void)
{
  int min = -size_;
  int max = size_;

  int NN = 2*size_ + 1;
  local1DIndexArrayForCompactStencil_ = 
    vector<TNT::Matrix<int> >(NN, Matrix<int>(NN, NN));
  
  points_.clear();

  int currentLocal1DIndex = 0;
  for (int x=min; x<=max; x++) {
    for (int y=min; y<=max; y++) {
      for (int z=min; z<=max; z++) {
	GridIndex grid(x, y, z);
	if (grid.length2() <= compactStencilSize_) {
	  points_.push_back(grid);
	  local1DIndexArrayForCompactStencil_[x+size_][y+size_][z+size_]
	    = currentLocal1DIndex ++;
	} else {
	  local1DIndexArrayForCompactStencil_[x+size_][y+size_][z+size_]
	    = -1;
	}
      }
    }
  }
}

/**********************************************************************
 * findGlobalIndexOfStencilPoint --
 **********************************************************************/
void
Stencil::findGlobalIndexOfStencilPoint(
				       const int ihx, 
				       const int ihy, 
				       const int ihz,
				       const Grid& grid,
				       vector<int>& globalIndexList) const
{
  // remove any points outside the grid range
  int xmin = max(ihx - size_, 0);
  int xmax = min(ihx + size_, static_cast<int>(grid.numPointX()-1));
  int ymin = max(ihy - size_, 0);
  int ymax = min(ihy + size_, static_cast<int>(grid.numPointY()-1));
  int zmin = max(ihz - size_, 0);
  int zmax = min(ihz + size_, static_cast<int>(grid.numPointZ()-1));

  for (int igx = xmin; igx <= xmax; igx++) {
    for (int igy = ymin; igy <= ymax; igy++) {
      for (int igz = zmin; igz <= zmax; igz++) {
	if (stencilShape_ == CUBE_STENCIL) {
	  globalIndexList.push_back(grid.transfer3DIndexTo1D(igx, igy, igz));
	} else if ( (stencilShape_ == SPHERE_STENCIL) || 
		    (stencilShape_ == COMPACT_STENCIL) ) {
	  if(isInsideStencil(ihx-igx, ihy-igy, ihz-igz)) {
	    globalIndexList.push_back(grid.transfer3DIndexTo1D(igx, igy, igz));
	  }
	} else {
	  errorMessage("Stencil::findGlobalIndexOfStencilPoint()", "unknown stencil shape!");
	}
      }
    }
  }
}

/**********************************************************************
 * findGlobal3DIndexOfStencilPoint --
 **********************************************************************/
void
Stencil::findGlobal3DIndexOfStencilPoint(
					 const GridIndex& hostPoint, 
					 const Grid& grid,
					 vector<GridIndex>& globalIndexList) const
{
  // remove any points outside the grid range
  int xmin = max(hostPoint.x() - size_, 0);
  int xmax = min(hostPoint.x() + size_, static_cast<int>(grid.numPointX()-1));
  int ymin = max(hostPoint.y() - size_, 0);
  int ymax = min(hostPoint.y() + size_, static_cast<int>(grid.numPointY()-1));
  int zmin = max(hostPoint.z() - size_, 0);
  int zmax = min(hostPoint.z() + size_, static_cast<int>(grid.numPointZ()-1));

  for (int igx = xmin; igx <= xmax; igx++) {
    for (int igy = ymin; igy <= ymax; igy++) {
      for (int igz = zmin; igz <= zmax; igz++) {
	if (stencilShape_ == CUBE_STENCIL) {
	  globalIndexList.push_back(GridIndex(igx, igy, igz));
	} else if ( (stencilShape_ == SPHERE_STENCIL) || 
		    (stencilShape_ == COMPACT_STENCIL) ) {
	  if(isInsideStencil(hostPoint.x()-igx, hostPoint.y()-igy, hostPoint.z()-igz)) {
	    globalIndexList.push_back(GridIndex(igx, igy, igz));
	  }
	} else {
	  errorMessage("Stencil::findGlobal3DIndexOfStencilPoint()", "unknown stencil shape!");
	}
      }
    }
  }
}

/**********************************************************************
 * isInsideStencil --
 **********************************************************************/
bool
Stencil::isInsideStencil (
			  const int x, 
			  const int y,
			  const int z) const
{
  GridIndex index3D(x, y, z);
  return isInsideStencil(index3D);
}

/**********************************************************************
 * isInsideStencil --
 **********************************************************************/
bool 
Stencil::isInsideStencil (
			  const GridIndex& index3D) const
{
  if (stencilShape_ == CUBE_STENCIL) {
    return (abs(index3D.x()) <= size_) && 
      (abs(index3D.y()) <= size_) && 
      (abs(index3D.z()) <= size_);
  } else if (stencilShape_ == SPHERE_STENCIL) {
    return index3D.length2() <= size_*size_;
  } else if (stencilShape_ == COMPACT_STENCIL) {
    return index3D.length2() <= compactStencilSize_;
  } else {
    errorMessage("Stencil::isInsideStencil()", "unknown stencil shape!");
  }
}

/**********************************************************************
 * setupComplete --
 **********************************************************************/
bool
Stencil::setupComplete (
			void) const
{
  if ( (stencilType_ == PROJECT_STENCIL) || (stencilType_ == INTERP_STENCIL) ) {
    return (! basisValue_.empty()) && (! weight_.empty());
  } else {
    return (! points_.empty());
  }
}

/**********************************************************************
 * setup --
 **********************************************************************/
void
Stencil::setup(
	       const double step) 
{
  step_ = step;
  if ( (stencilType_ == PROJECT_STENCIL) || (stencilType_ == INTERP_STENCIL) ) {
    setupEquaScale();
    setupPolyOrder();
    setupBasisValue();
    setupWeight();
  }
}

/**********************************************************************
 * setupEquaScale --
 **********************************************************************/
void
Stencil::setupEquaScale (
			 void)
{
  double scalarForQuartic[6] = {1.0, 0.70, 0.43, 0.38, 0.005, 0.075};

  equaScale_ = vector<double>(totalNumPoint_, 1.);
  if (equaScaleFlag_) {
    int pointIndex = 0;
    for (vector<GridIndex>::const_iterator it = points_.begin(); 
	 it != points_.end(); ++it, pointIndex++) {
      int distanceToOrigin = it->length2();
      if (size_ == 2) {
	if (distanceToOrigin <= compactStencilSize_) {
	  equaScale_[pointIndex] = scalarForQuartic[distanceToOrigin];
	} else {
	  errorMessage("Stencil::setupEquaScale()", "grid point index error!");
	}
      } else {
	// For higher order, I don't have scaling factor. This needs further investigation
	equaScale_[pointIndex] = 1.;
      }
    }
  }
}

/**********************************************************************
 * setupPolyOrder --
 **********************************************************************/
void
Stencil::setupPolyOrder (
			 void)
{
  xMonoOrder_.reserve(numBasis_);
  yMonoOrder_.reserve(numBasis_);
  zMonoOrder_.reserve(numBasis_);

  switch (basisType_) {
  case CUBED_POLY:
    maxPolyOrder_ = 3*maxMonoOrder_;
    for (int x = 0; x <= maxMonoOrder_; x++) 
      for (int y = 0; y <= maxMonoOrder_; y++) 
	for (int z = 0; z <= maxMonoOrder_; z++) {
	  xMonoOrder_.push_back(x);
	  yMonoOrder_.push_back(y);
	  zMonoOrder_.push_back(z);
	}
    break;
  case CONSISTENT_POLY:
    maxPolyOrder_ = maxMonoOrder_;
    for (int x = 0; x <= maxMonoOrder_; x++) 
      for (int y = 0; y <= maxMonoOrder_ - x; y++) 
	for (int z = 0; z <= maxMonoOrder_ - x - y; z++) {
	  xMonoOrder_.push_back(x);
	  yMonoOrder_.push_back(y);
	  zMonoOrder_.push_back(z);
	}
    break;
  default:
    errorMessage("Stencil::setupPolyOrder()", "unknown basis type!");
    break;
  }

  if (static_cast<int>(xMonoOrder_.size()) != numBasis_)
    errorMessage("Stencil::setupPolyOrder()", "inconsistent numBasis and polyTermIndex!");
}

/**********************************************************************
 * setupBasisValue --
 **********************************************************************/
void
Stencil::setupBasisValue (
			  void)
{
  basisValue_ = TNT::Matrix<double>(totalNumPoint_, numBasis_);

  /* Testing in Matlab shows that the matrix 
   * starts to lose rank (numerically, not analytically) for dx<0.005.
   * At this point the result of the present computations would 
   * contain errors of RELATIVE magnitude one. The same problem occurs 
   * when dx becomes large. Obviously, the problem is a poor scaling 
   * of the columns of the "values" matrix. As a consequence a scaling 
   * of the values matrix is needed to ensure that the results are 
   * decent for all values of dx.
   */
  polyScale_ = 1./step_;
  vector<double> poly(numBasis_);
  for (int pointIndex = 0; pointIndex < totalNumPoint_; pointIndex++) {
    double x = points_[pointIndex].x() * step_;
    double y = points_[pointIndex].y() * step_;
    double z = points_[pointIndex].z() * step_;
    evaluateBasisFuncAtOnePoint(x, y, z, poly);
    for (int polyTermIndex = 0; polyTermIndex < numBasis_; polyTermIndex++) {
      basisValue_[pointIndex][polyTermIndex] = poly[polyTermIndex] * equaScale_[pointIndex];
    }
  }
}

/**********************************************************************
 * evaluateBasisFuncAtOnePoint --
 * This function is used by Stencil only.
 **********************************************************************/
void
Stencil::evaluateBasisFuncAtOnePoint (
				      double x,
				      double y,
				      double z,
				      vector<double>& funcValue) const
{
  x *= polyScale_;
  y *= polyScale_;
  z *= polyScale_;
  for (int polyTermIndex = 0; polyTermIndex < numBasis_; polyTermIndex++) {
    funcValue[polyTermIndex] = 
      intpow(x, xMonoOrder_[polyTermIndex]) *
      intpow(y, yMonoOrder_[polyTermIndex]) *
      intpow(z, zMonoOrder_[polyTermIndex]);
  }
}

/**********************************************************************
 * evaluateBasisFuncAtOnePoint --
 * This function is used by interpMat.
 **********************************************************************/
void
Stencil::evaluateBasisFuncAtOnePoint (
				      point3D& point,
				      const vector3D<double>& normal,
				      TNT::Vector<double>& funcValue,
				      const DifferentialOperator operatorType) const
{
  point *= polyScale_;
  for (int polyTermIndex = 0; polyTermIndex < numBasis_; polyTermIndex++) {
    switch (operatorType) {
    case NONE:
      funcValue[polyTermIndex] = basisFuncValue(point, polyTermIndex);
      break;
    case D_DN:
      funcValue[polyTermIndex] = 
	normal.x() * basisFuncValue_dx(point, polyTermIndex) + 
	normal.y() * basisFuncValue_dy(point, polyTermIndex) + 
	normal.z() * basisFuncValue_dz(point, polyTermIndex);
      break;
    case D_DX:
      funcValue[polyTermIndex] = basisFuncValue_dx(point, polyTermIndex);
      break;
    case D_DY:
      funcValue[polyTermIndex] = basisFuncValue_dy(point, polyTermIndex);
      break;
    case D_DZ:
      funcValue[polyTermIndex] = basisFuncValue_dz(point, polyTermIndex);
      break;
    default:
      errorMessage("Stencil::evaluateBasisFuncAtOnePoint()", 
		   "unknown stencil basis type!");
    }
  }
}

/**********************************************************************
 * compCoef --
 * This function is called by both interpMat and projectMat.
 * It computes a vector matrix product, vec * weight_
 * where vec contains either the value of many ploynomial terms evaluated
 * at collocation point (for the case of interpMat), or the value of
 * those polynomila terms integrated over the element (for the case of 
 * projectMat).
 * I don't want to reveal the internal data weight_. So I add this 
 * function to stencil class.
 **********************************************************************/
void
Stencil::compCoef (
		   const TNT::Vector<double>& vec,
		   TNT::Vector<double>& coef) const
{
  coef = vec * weight_;  
}


/**********************************************************************
 * setupWeight --
 **********************************************************************/
void
Stencil::setupWeight (
		      void)
{
  int size = max(numBasis_, totalNumPoint_);
  double **tmp;
  tmp = new double*[size];
  for (int i=0; i<size; i++) tmp[i] = new double[size];

  for (int pointIndex = 0; pointIndex < totalNumPoint_; pointIndex++) 
    for (int polyTermIndex = 0; polyTermIndex < numBasis_; polyTermIndex++) 
      tmp[pointIndex][polyTermIndex] = basisValue_[pointIndex][polyTermIndex];
  pseudo_inverse(tmp, totalNumPoint_, numBasis_);

  weight_ = TNT::Matrix<double>(numBasis_, totalNumPoint_);
  // rescale
  for (int polyTermIndex = 0; polyTermIndex < numBasis_; polyTermIndex++) 
    for (int pointIndex = 0; pointIndex < totalNumPoint_; pointIndex++) 
      weight_[polyTermIndex][pointIndex] = 
	tmp[polyTermIndex][pointIndex] * equaScale_[pointIndex];	

  for (int i=0; i<size; i++) delete [] tmp[i];
  delete [] tmp;
}

/**********************************************************************
 * integrateBasisFuncOverElement --
 * moments: integral of basis Functions, the monomials, over an element
 **********************************************************************/
void
Stencil::integrateBasisFuncOverElement (
					const element& ele,
					const point3D& origin,
					const vector3D<double>& normal,
					TNT::Vector<double>& moments,
					const DifferentialOperator operatorType)
{
  // I will add to this vector. So initializing to zero is necessary
  for (size_t i=0; i < moments.size(); i++)
    moments[i] = 0.;

  // find cubature points and weight in Berycentric coordinate system;
  if (! cubature_.setupComplete()) {
    int degree = maxPolyOrder_ + 1;
    int rule = 1;
    cubature_ = Cubature(degree, rule, Cubature::TRIANGLE);
  }

  // divide the element into many triangles because it is easier to use cubature
  vector<point3D> vertex(3);
  // adjust element verteices to the new origin
  // and scale the coordinate in accordance with similar precedure in 
  // computing weight_
  vertex[0] = (ele.vertex(0) - origin) * polyScale_;
  for (size_t i=1; i < ele.shape() - 1; i++) {
    vertex[1] = (ele.vertex(i) - origin) * polyScale_;
    vertex[2] = (ele.vertex(i+1) - origin) * polyScale_;
    addIntegralOverTriangle(vertex, normal, moments, operatorType);
  }
}

/**********************************************************************
 * addIntegralOverTriangle --
 **********************************************************************/
void
Stencil::addIntegralOverTriangle (
				  const vector<point3D>& vertex, 
				  const vector3D<double>& normal,
				  TNT::Vector<double>& moments,
				  const DifferentialOperator operatorType) const
{
  double area = compTriangleArea(vertex[0], vertex[1], vertex[2])
    / (0.5 * polyScale_ * polyScale_);

  for (int pointIndex = 0; pointIndex < cubature_.numPoint() ; pointIndex++) {
    point3D point = cubature_.point(vertex, pointIndex);
    double weight = cubature_.weight(area, pointIndex);
    for (int polyTermIndex = 0; polyTermIndex < numBasis_; polyTermIndex++) {

      switch (operatorType) {
      case NONE:
	moments[polyTermIndex] += basisFuncValue(point, polyTermIndex) * weight;
	break;
      case D_DN:
	moments[polyTermIndex] += 
	  (normal.x() * basisFuncValue_dx(point, polyTermIndex) + 
	   normal.y() * basisFuncValue_dy(point, polyTermIndex) + 
	   normal.z() * basisFuncValue_dz(point, polyTermIndex) ) 
	   * weight;
	break;
      case D_DX:
	moments[polyTermIndex] += basisFuncValue_dx(point, polyTermIndex) * weight;
	break;
      case D_DY:
	moments[polyTermIndex] += basisFuncValue_dy(point, polyTermIndex) * weight;
	break;
      case D_DZ:
	moments[polyTermIndex] += basisFuncValue_dz(point, polyTermIndex) * weight;
      break;
      default:
	errorMessage("Stencil::addIntegralOverTriangle()", 
		     "unknown stencil basis type!");
      }
    }
  }
}

/**********************************************************************
 * compTriangleArea --
 **********************************************************************/
inline double
Stencil::compTriangleArea (
			   const point3D& v1,
			   const point3D& v2,
			   const point3D& v3) const 
{
  vector3D<double> edge1 = v2 - v1;
  vector3D<double> edge2 = v2 - v3;
  vector3D<double> vec = crossProd(edge1, edge2);
  return 0.5 * length(vec);
}

/**********************************************************************
 * basisFuncValue --
 **********************************************************************/
inline double
Stencil::basisFuncValue (
			 const point3D& point,
			 const int polyTermIndex) const
{
  return 
    intpow(point.x(), xMonoOrder_[polyTermIndex]) *
    intpow(point.y(), yMonoOrder_[polyTermIndex]) *
    intpow(point.z(), zMonoOrder_[polyTermIndex]);

}

/**********************************************************************
 * basisFuncValue_dx --
 **********************************************************************/
inline double
Stencil::basisFuncValue_dx (
			    const point3D& point,
			    const int polyTermIndex) const
{
  if (xMonoOrder_[polyTermIndex] == 0) {
    return 0.;
  } else {
    return
      xMonoOrder_[polyTermIndex] * intpow(point.x(), xMonoOrder_[polyTermIndex]-1) *
      intpow(point.y(), yMonoOrder_[polyTermIndex]) *
      intpow(point.z(), zMonoOrder_[polyTermIndex]) * polyScale_;
  }
}

/**********************************************************************
 * basisFuncValue_dy --
 **********************************************************************/
inline double
Stencil::basisFuncValue_dy (
			    const point3D& point,
			    const int polyTermIndex) const
{
  if (yMonoOrder_[polyTermIndex] == 0) {
    return 0.;
  } else {
    return
      intpow(point.x(), xMonoOrder_[polyTermIndex]) *
      yMonoOrder_[polyTermIndex] * intpow(point.y(), yMonoOrder_[polyTermIndex]-1) *
      intpow(point.z(), zMonoOrder_[polyTermIndex]) * polyScale_;
  }
}

/**********************************************************************
 * basisFuncValue_dz --
 **********************************************************************/
inline double
Stencil::basisFuncValue_dz (
			    const point3D& point,
			    const int polyTermIndex) const
{
  if (zMonoOrder_[polyTermIndex] == 0) {
    return 0.;
  } else {
    return
      intpow(point.x(), xMonoOrder_[polyTermIndex]) *
      intpow(point.y(), yMonoOrder_[polyTermIndex]) *
      zMonoOrder_[polyTermIndex] * intpow(point.z(), zMonoOrder_[polyTermIndex]-1) * polyScale_;
  }
}

/**********************************************************************
 * intpow --
 * By Dennis Wu
 **********************************************************************/
inline double 
Stencil::intpow(
		const double d, 
		int e) const 
{
  double product = 1.;
  double increment = d;
  
  while(e>0){
    if ((e%2)!=0){
      product = product * increment;
      e=e-1;
    }else{
      increment=increment*increment;
      e=e/2;
    }
  }
  return product;
}


/*
    // for s012346
    for (vector<GridIndex>::const_iterator it = points_.begin(); 
	 it != points_.end(); ++it, pointIndex++) {
      switch (it->length2()) {
      case 0:
	equaScale_[pointIndex] = 1.00;
	break;
      case 1:
	equaScale_[pointIndex] = 0.75;
	break;
      case 2:
	equaScale_[pointIndex] = 0.56;
	break;
      case 3:
	equaScale_[pointIndex] = 0.65;
	break;
      case 4:
	equaScale_[pointIndex] = 0.10;
	break;
      case 5:
	equaScale_[pointIndex] = 0.0;
	break;
      case 6:
	equaScale_[pointIndex] = 0.062;
	break;
      default:
	errorMessage("Stencil::setupEquaScale()", "grid point index error!");
	break;
      }
    }
*/
