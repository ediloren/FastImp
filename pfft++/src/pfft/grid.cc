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

const static char cvsid[] = "$Id: grid.cc,v 1.6 2003/02/11 03:11:30 zhzhu Exp $";

#include "grid.h"
#include "gridElement.h"
#include "stencil.h"

using namespace pfft;
using namespace std;

/**********************************************************************
 * Grid --
 **********************************************************************/
Grid::Grid (
	    GridElement* srcGridElement,
	    GridElement* evalGridElement,
	    Stencil* projectStencil,
	    Stencil* interpStencil,
	    Stencil* directStencil,
	    double maxGridSize)
  : srcGridElementPtr(srcGridElement), evalGridElementPtr(evalGridElement), 
  projectStencilPtr(projectStencil), interpStencilPtr(interpStencil),
  directStencilPtr(directStencil), maxGridSize_(maxGridSize)
{
  averageElementSize_ = (srcGridElementPtr->averageElementSize() + 
			 evalGridElementPtr->averageElementSize() ) / 2.;
  minElementSize_ = min(srcGridElementPtr->minElementSize(),
			evalGridElementPtr->minElementSize());
  maxElementSize_ = max(srcGridElementPtr->maxElementSize(),
			evalGridElementPtr->maxElementSize());

  numPointX_ = numPointY_ = numPointZ_ = 4;
  bakNumPointX_ = bakNumPointY_ = bakNumPointZ_ = numPointX_;
  prevGridVolume_ = currentGridVolume_ = DBL_MAX;

  // find the range of problem domain, which includes src and eval gridElements
  domain_xmin_ = min(srcGridElementPtr->xmin(), evalGridElementPtr->xmin());
  domain_xmax_ = max(srcGridElementPtr->xmax(), evalGridElementPtr->xmax());
  domain_ymin_ = min(srcGridElementPtr->ymin(), evalGridElementPtr->ymin());
  domain_ymax_ = max(srcGridElementPtr->ymax(), evalGridElementPtr->ymax());
  domain_zmin_ = min(srcGridElementPtr->zmin(), evalGridElementPtr->zmin());
  domain_zmax_ = max(srcGridElementPtr->zmax(), evalGridElementPtr->zmax());

  // Add a (very) little slack to the grid, so that all elements are fully inside the
  // grid, i.e. no part of any element is on the border of the grid.
  domain_xmin_ -= DBL_EPSILON*abs(domain_xmin_);
  domain_xmax_ += DBL_EPSILON*abs(domain_xmax_);
  domain_ymin_ -= DBL_EPSILON*abs(domain_ymin_);
  domain_ymax_ += DBL_EPSILON*abs(domain_ymax_);
  domain_zmin_ -= DBL_EPSILON*abs(domain_zmin_);
  domain_zmax_ += DBL_EPSILON*abs(domain_zmax_);

  // Add extra layers to make sure interpolation and projection stencils are
  // still cube stencil.
  extraLayer_ = max(projectStencilPtr->size(), interpStencilPtr->size()) - 0.5;
}

/**********************************************************************
 * setupGridSize --
 **********************************************************************/
void 
Grid::setupGridSize (
		     void)
{
  // make sure the total number is bigger than number of extra layers on both sides
  while (numPointX_ <= static_cast<size_t>(ceil(2*extraLayer_))+1) numPointX_ *= 2;
  while (numPointY_ <= static_cast<size_t>(ceil(2*extraLayer_))+1) numPointY_ *= 2;
  while (numPointZ_ <= static_cast<size_t>(ceil(2*extraLayer_))+1) numPointZ_ *= 2;

  // extra layers are counted in the numPointX, but xmax and xmin should not include them
  // The extra layers are meant to wrap the computational domain, not to be part of it
  double dx = (domain_xmax_ - domain_xmin_) / (numPointX_ - 2*extraLayer_ - 1);
  double dy = (domain_ymax_ - domain_ymin_) / (numPointY_ - 2*extraLayer_ - 1);
  double dz = (domain_zmax_ - domain_zmin_) / (numPointZ_ - 2*extraLayer_ - 1);
  gridSize_ = max(dx, dy);
  gridSize_ = max(gridSize_, dz);

  // Find out along which direction the elements fit tightly with grids
  // This direction will be the potential direction along which the grid number
  // should be increased if necessary
  if (dx == gridSize_) tightFit_ = 0;
  else if (dy == gridSize_) tightFit_ = 1;
  else  tightFit_ = 2;

  // the uniform grid is used so that we could exploit the shift invariance of Green's 
  // function and use FFT
  dx_ = gridSize_;
  dy_ = gridSize_;
  dz_ = gridSize_;

  // add extra layer(s) to each side
  xmin_ = domain_xmin_ - extraLayer_ * dx_; 
  ymin_ = domain_ymin_ - extraLayer_ * dy_;
  zmin_ = domain_zmin_ - extraLayer_ * dz_;
  xmax_ = domain_xmin_ + (numPointX_ - 1) * dx_;
  ymax_ = domain_ymin_ + (numPointY_ - 1) * dy_;
  zmax_ = domain_zmin_ + (numPointZ_ - 1) * dz_;

  x_.resize(numPointX_);
  for (size_t i = 0; i < numPointX_; i++)
    x_[i] = xmin_ + i * dx_;

  y_.resize(numPointY_);
  for (size_t i = 0; i < numPointY_; i++)
    y_[i] = ymin_ + i * dy_;

  z_.resize(numPointZ_);
  for (size_t i = 0; i < numPointZ_; i++)
    z_[i] = zmin_ + i * dz_;

  gridReport();
}

/**********************************************************************
 * checkGridSize --
 * 
 **********************************************************************/
Grid::GridStatus 
Grid::checkGridSize (
		     void)
{
  // some hard-coded numbers used here
  //  const double minGridElementRatio = 0.5;
  //  const double maxGridElementRatio = 10;
  //  const double minInteractVolumeRatio = 0.5;
  //  const size_t maxNumElementToOneGrid = 110;

  // if a grid has too many elements mapped to it, it means the grid is too coarse
  // I could use rather generous threshold here because the memory usage check will 
  // be more precise
  size_t maxNumElementToOneGrid = 110;
  size_t maxNum = max(srcGridElementPtr->maxNumElementMappedToOneGrid(), 
		      evalGridElementPtr->maxNumElementMappedToOneGrid());
  if (maxNum > maxNumElementToOneGrid) {
    debugReport("Max number of elements mapped to one grid := ", maxNum);
    debugReport("Grid too coarse: too many elements associated with one grid point.");
    return TOO_COARSE;
  }

  // Make sure the domain occupied by current grid is no larger than previous one
  // Ideally, we want the grids barely cover the computational domain, i.e. we want to
  // minimize the grid domain.
  prevGridVolume_ = currentGridVolume_;
  currentGridVolume_ = (xmax_-xmin_) * (ymax_-ymin_) * (zmax_-zmin_);
  if (currentGridVolume_ > prevGridVolume_) {
    debugReport("Grid too coarse: current grid volume is larger than previous one");
    return TOO_COARSE;
  }

  // this check makes sure that the error introduced in the interpolation is no
  // worse than the error in discretization. I assume the piece-wise constant elements are
  // used. Therefore, this is only a rough bound.
  size_t minInterpOrder = 2 * min(projectStencilPtr->size(), interpStencilPtr->size());
  if (gridSize_ > minInterpOrder * maxElementSize_) {
    debugReport("Grid too coarse: interpolation error in pfft might surpass discretization error");
    return TOO_COARSE; 
  }

  // The ideal case is that no extrapolation is involved. This check ensures that 
  // an average element do not hang out of the intended interpolation grid by more than 
  // one grid size. Experiments show that use average is better than maxElement. 
  // The latter is more conservative and hence produces larger grid step. This compromises
  // the smaller element's accuracy.
  double minStencilSize = min(projectStencilPtr->size(), interpStencilPtr->size());
  if ((minStencilSize+1) * gridSize_ < averageElementSize_) {
    debugReport("Grid too fine: extrapolation might be used for some large elements");
    return TOO_FINE;
  }

  // make sure the interactive volume defined by directStencil is not too small
  // comparing to the interactive volume defined by element size.
  // This checck ensures that the number of irregular neighbors are not excessive.
  // The cost of pre-correction for the irregular neighbors is higher than the 
  // regular neighbors.
  const double minInteractVolumeRatio = 0.5;
  double stencilVolume = directStencilPtr->numPoint() * dx_ * dy_ * dz_;
  double interactVolumeRadius = averageElementSize_ * (1. + separationFactor());
  double interactVolume = 4./3*PAI*pow(interactVolumeRadius, 3);
  if (stencilVolume < interactVolume * minInteractVolumeRatio) {
    debugReport("Grid too fine: interactive volume defined by stencil is too small");
    return TOO_FINE;
  }

  // Make sure the gridSize is no larger than a user-specified upper bound
  // If user does not specify, this constrain has no effect
  // This is for oscillatory kernels, in which case the grid size should be no larger than
  // a tenth of a wavelength. Normally a good ratio between grid and element size should
  // guarantee this. But this parameter provides a double check.
  if (gridSize_ > maxGridSize_) {
    debugReport("Grid too coarse: grid size is larger than user specified threshold");
    return TOO_COARSE;
  }

  debugReport("Grid size OK");
  return SIZE_OK;
}

/**********************************************************************
 * doubleNumPoint --
 **********************************************************************/
void
Grid::doubleNumPoint (
		      void)
{
  if (tightFit_ == 0) 
    numPointX_ *= 2;
  else if (tightFit_ == 1) 
    numPointY_ *= 2;
  else  
    numPointZ_ *= 2;
}

/**********************************************************************
 * transfer1DIndexTo3D --
 **********************************************************************/
GridIndex
Grid::transfer1DIndexTo3D (
			   int index1D) const
{
  int x = static_cast<int>( floor(index1D / static_cast<double>(2*numPointZ_) / (2*numPointY_)) );
  index1D -= x * (2*numPointZ_) * (2*numPointY_);
  int y = static_cast<int>( floor(index1D / static_cast<double>((2*numPointZ_))) );
  int z = index1D - y * (2*numPointZ_);

  return GridIndex(x, y, z);
}

/**********************************************************************
 * gridPointDistance --
 **********************************************************************/
double
Grid::distance(
	       const int ix1, 
	       const int iy1, 
	       const int iz1,
	       const int ix2, 
	       const int iy2, 
	       const int iz2) const
{
  double x = dx_ * (ix2 - ix1);
  double y = dy_ * (iy2 - iy1);
  double z = dz_ * (iz2 - iz1);
  return sqrt(x*x + y*y + z*z);
}

/**********************************************************************
 * debugReport --
 **********************************************************************/
inline void 
Grid::debugReport (
		   string message)
{
#ifdef DEBUG_GRID
    cout << message << endl;
#endif
}

/**********************************************************************
 * debugReport --
 **********************************************************************/
template <class T>
inline void 
Grid::debugReport (
		   const string message,
		   const T number)
{
#ifdef DEBUG_GRID
  cout << message << " " << number << endl;
#endif
}

/**********************************************************************
 * gridReport --
 **********************************************************************/
void
Grid::gridReport (
		  void) const
{
#ifdef DEBUG_GRID
  vector3D<double> problemDomain_lu(domain_xmin_, domain_ymin_, domain_zmin_);
  vector3D<double> problemDomain_rl(domain_xmax_, domain_ymax_, domain_zmax_);
  vector3D<double> gridDomain_lu(xmin_, ymin_, zmin_);
  vector3D<double> gridDomain_rl(xmax_, ymax_, zmax_);
  cout << "number of grids := (" 
       << numPointX_ << "," << numPointY_ << "," << numPointZ_ << ")"
       << "  grid size := " << gridSize_ << endl
       << "problem domain := " << problemDomain_lu << "--" << problemDomain_rl << endl
       << "grid domain := " << gridDomain_lu << "--" << gridDomain_rl << endl;
#endif
}

