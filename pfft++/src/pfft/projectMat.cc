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

const static char cvsid[] = "$Id: projectMat.cc,v 1.4 2003/02/11 03:13:48 zhzhu Exp $";

#include "projectMat.h"
#include <algorithm>

using namespace std;
using namespace pfft;

/**********************************************************************
* memoryEstimation --
**********************************************************************/
size_t
ProjectMat::memoryEstimation (
			      void) const
{
  int numEntry = srcGridElementPtr->numElement() * projectStencilPtr->numPoint();
  int numCol = srcGridElementPtr->numElement();
  return matList[0].memoryEstimate(numEntry, numCol) * numMat();
}

/**********************************************************************
 * setup --
 **********************************************************************/
void
ProjectMat::setup (
		   void)
{
  if (projectStencilPtr->setupComplete()) {
    for (size_t i = 0; i < numMat(); i++) {
      setupEachMatrix(integrationPtr->innerOperator(i), i);
    }
  }
}

/**********************************************************************
* setupEachMatrix --
**********************************************************************/
void
ProjectMat::setupEachMatrix (
			     const DifferentialOperator operatorType,
			     const size_t matIndex)
{
  TNT::Vector<double> integralOfBasisFunc(projectStencilPtr->numBasis());
  TNT::Vector<double> denseCol(projectStencilPtr->numPoint());
  matList[matIndex] = SpColMat<double>(srcGridElementPtr->numElement(), 
				       gridPtr->totalNumPointAfterPadded());

  size_t globalElementIndex;

  for (size_t ix=0; ix < gridPtr->numPointX(); ix++) {
    for (size_t iy=0; iy < gridPtr->numPointY(); iy++) {
      for (size_t iz=0; iz < gridPtr->numPointZ(); iz++) {
	int gridPointIndex = gridPtr->transfer3DIndexTo1D(ix, iy, iz);
	size_t numElement = 
	  srcGridElementPtr->numElementMappedToOneGrid(gridPointIndex);
	if (! numElement) continue;

	std::vector<int> globalIndexList;
	projectStencilPtr->findGlobalIndexOfStencilPoint(ix, iy, iz, *gridPtr, 
							 globalIndexList);

	point3D origin(gridPtr->x(ix), gridPtr->y(iy), gridPtr->z(iz));
	for (size_t localElementIndex = 0; localElementIndex < numElement; 
	     localElementIndex++) {
	  globalElementIndex = 
	    srcGridElementPtr->getElementMappedToOneGrid(gridPointIndex, 
							 localElementIndex);
	  const element& evalElement = 
	    srcGridElementPtr->getElement(globalElementIndex);

	  // Note that the coordinate origin used in project stencil is the 
	  // center point. Therefore, the element vertices should be 
	  // adjusted accordingly. This should be done in the following function.
	  projectStencilPtr->integrateBasisFuncOverElement(evalElement, origin, 
							   evalElement.normal(),
							   integralOfBasisFunc,
							   operatorType);
	  projectStencilPtr->compCoef(integralOfBasisFunc, denseCol);
	  copyFullVectorToSpVec(globalIndexList.begin(), globalIndexList.end(), 
				denseCol.begin(), matList[matIndex][globalElementIndex]);
	}
      }
    }
  }
}
