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

  const static char cvsid[] = "$Id: grid.h,v 1.8 2003/03/25 15:14:13 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _GRID_H_
#define _GRID_H_

#include "stencil.h"
#include "utils.h"
#include "gridIndex.h"
#include <vector>
#include <string>
#include <float.h>

namespace pfft {

  class GridElement;

  class Grid {
  public:
    enum GridStatus { SIZE_OK, TOO_COARSE, TOO_FINE };

    Grid (GridElement*, GridElement*, Stencil*, Stencil*, Stencil*,
	  double maxGridSize = DBL_MAX  // user could specify a constrain on grid size
	  );
    Grid (void) {};

    void setupGridSize();
    GridStatus checkGridSize(void);
    void doubleNumPoint(void);
    void backupNumPoint(void) {
      bakNumPointX_ = numPointX_; bakNumPointY_ = numPointY_; bakNumPointZ_ = numPointZ_; 
      bakGridSize_ = gridSize_; }
    void recoverNumPoint(void) {
      numPointX_ = bakNumPointX_; numPointY_ = bakNumPointY_; numPointZ_ = bakNumPointZ_; 
      gridSize_ = bakGridSize_; }
    
    // The error bound is:
    // regular neigbors: (1/d)^(n+1), where d is the directStencilSize
    // irregular neighbors: (r/(r+s))^(n+1), where r is panel size and 
    //                      s is the separation distance between src and eval elements.
    // To make these two on the same order, I pick s/r as the following
    //    double separationFactor(void) const { return (directStencilPtr->size() - 1); }
    double separationFactor(void) const { return 2.; }

    double xmin(void) const {return xmin_;}
    double ymin(void) const {return ymin_;}
    double zmin(void) const {return zmin_;}
    double xmax(void) const {return xmax_;}
    double ymax(void) const {return ymax_;}
    double zmax(void) const {return zmax_;}
    double dx(void) const {return dx_;}
    double dy(void) const {return dy_;}
    double dz(void) const {return dz_;}
    double gridStep(void) const {return gridSize_;}
    size_t numPointX(void) const {return numPointX_;}
    size_t numPointY(void) const {return numPointY_;}
    size_t numPointZ(void) const {return numPointZ_;}
    size_t totalNumPoint(void) const {return numPointX()*numPointY()*numPointZ();}
    size_t totalNumPointAfterPadded(void) const {
      return 4*numPointX()*numPointY()*numPointZ(); 
    }

    int transfer3DIndexTo1D (const int x, const int y, const int z) const {
      return z + (2*numPointZ_) * (y + (2*numPointY_) * x); }

    int transfer3DIndexTo1D (const GridIndex& gridIndex) const {
      return gridIndex.z() + (2*numPointZ_) * (gridIndex.y() + (2*numPointY_) * gridIndex.x()); }

    GridIndex transfer1DIndexTo3D (int index1D) const;
    double distance(const int ix1, const int iy1, const int iz1,
		    const int ix2, const int iy2, const int iz2) const;

    double x(size_t i) { return x_[i]; }
    double y(size_t i) { return y_[i]; }
    double z(size_t i) { return z_[i]; }
  
  private:
    // state variables
    GridElement* srcGridElementPtr;
    GridElement* evalGridElementPtr;
    Stencil* projectStencilPtr;
    Stencil* interpStencilPtr;
    Stencil* directStencilPtr;
    double maxGridSize_; // gridSize upper bound
    size_t numPointX_, numPointY_, numPointZ_;
    double dx_, dy_, dz_, gridSize_, bakGridSize_;
    double xmin_, xmax_, ymin_, ymax_, zmin_, zmax_; // grid domain
    std::vector<double> x_, y_, z_;
    double domain_xmin_, domain_xmax_;  // computational domain
    double domain_ymin_, domain_ymax_;
    double domain_zmin_, domain_zmax_;
    double extraLayer_;

    // working variables
    size_t bakNumPointX_, bakNumPointY_, bakNumPointZ_;
    double averageElementSize_, minElementSize_, maxElementSize_;
    int tightFit_; // along which direction hte grid fit elements tightly,  0:x, 1:y, 2:z
    double prevGridVolume_, currentGridVolume_;

    inline void debugReport(std::string message);
    template <class T> 
    inline void Grid::debugReport (const std::string message, const T number);
    inline void gridReport (void) const;
  };

} // namespace pfft

#endif
