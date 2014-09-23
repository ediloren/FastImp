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

  The idea of pre-correstion is to subtract the inaccurate near neighbor 
  intercation, which is computed through projection, grid-to-grid convolution 
  and interpolation, from the accurate near neighbor interaction, which is
  computed directly using panel integration.

  The convolution in this stage is different from the one in setting 
  up the system matrix. In the formal one I want to compute the contribution 
  to one eval panel from another individual source panel so that I could 
  subtract it from an entry in the direct matrix.
  This is a one-to-one panel interaction. Only a small number of grid points 
  are involved. So FFT will not speed up the process much.
  In setting up the  sytem matrix, I want to compute the contribution of 
  all source panels to one eval panel. This multiple-to-one interaction
  involves a large amount of grid points, hence FFT could speed up the
  process significantly.

  However, the invariance of the Green's function with respect to the position
  could still be taken advantage of in the pre-correction stage. And this is 
  the key idea that I could use to speed up this step.

  For a grid point in the direct stencil (denoted as A), a small matrix with the 
  grid-to-grid interactions between the projection stencil centered at the 
  center of the direct stencil and the interpolation stencil centered at 
  the grid point A is pre-computed here. One matrix is computed for each direct
  stencil point. And the size of each small matrix is numInterpStencilPoint by
  numProjectStencilPoint. These small matrices are independent of the source and 
  eval panels. Therefore, we could pre-compue them and use them whenever 
  appropriate. The interaction between two panels is simply interp*G2G*project.

  This idea works perfectly for small panels. However, for large source panels, who
  have some irregular neighboring panels that can not be reached "through the grid", 
  we have to handle them seperately. The idea is again to find interp, G2G and 
  project matrices and compute the individual contribution. But since the 
  interpolation grid points involved are dependent on eval panel, I could not 
  pre-compute G2G mtrices any more. So the G2g matrix has to be computed on 
  the fly.
  
  Resources:

  See also:

  const static char cvsid[] = "$Id: grid2grid.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _GRID_2_GRID_H_
#define _GRID_2_GRID_H_

#include "gridIndex.h"
#include "stencil.h"
#include <vector>

#include "cmat.h" // for TNT::Matrix

namespace pfft {

  template <class GreenFunc, class KT>

  class Grid2Grid {

  public:

    Grid2Grid(const Stencil&, const Stencil&, const Stencil&, const GreenFunc&);
    Grid2Grid(void) {};

    const TNT::Matrix<KT>& matrix(size_t index) const { return g2gMat[index]; }
    size_t numMat(void) const { return numMat_; }

  private:
    std::vector<TNT::Matrix<KT> > g2gMat;
    size_t numMat_;
  };

  /**********************************************************************
   * Grid2Grid --
   **********************************************************************/
  template <class GreenFunc, class KT>
  Grid2Grid<GreenFunc, KT>::
  Grid2Grid (
	     const Stencil& projectStencil, 
	     const Stencil& interpStencil, 
	     const Stencil& directStencil,
	     const GreenFunc& greenFunc)
    : numMat_(directStencil.numPoint())
  {
    g2gMat.reserve(directStencil.numPoint());

    for (int directStencilPointIndex = 0; 
	 directStencilPointIndex < directStencil.numPoint(); 
	 directStencilPointIndex++) {
      TNT::Matrix<KT> g2g(interpStencil.numPoint(), projectStencil.numPoint());

      // each directStencil point serves as the center for a new interpStencil
      GridIndex interpStencilCenter = directStencil.point(directStencilPointIndex);
      for (int interpStencilPointIndex = 0; 
	   interpStencilPointIndex < interpStencil.numPoint(); 
	   interpStencilPointIndex++) {
	GridIndex interpPoint = interpStencil.point(interpStencilPointIndex);
	// shift to the coordinate center at directStencil center
	interpPoint += interpStencilCenter;
      
	int rowIndex = interpStencilPointIndex;
	for (int projectStencilPointIndex = 0; 
	     projectStencilPointIndex < projectStencil.numPoint(); 
	     projectStencilPointIndex++) {
	  GridIndex projectPoint = projectStencil.point(projectStencilPointIndex);
	  double dx = (projectPoint.x() - interpPoint.x()) * directStencil.step();
	  double dy = (projectPoint.y() - interpPoint.y()) * directStencil.step();
	  double dz = (projectPoint.z() - interpPoint.z()) * directStencil.step();
	  int colIndex = projectStencilPointIndex;
	  g2g[rowIndex][colIndex] = greenFunc(dx, dy, dz);
	}
      }

      g2gMat.push_back(g2g);
    }
  }


} //namespace pfft

#endif
