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

const static char cvsid[] = "$Id: interpMat.cc,v 1.2 2003/02/11 03:13:48 zhzhu Exp $";

#include "interpMat.h"
#include <algorithm>

using namespace std;
using namespace pfft;

/**********************************************************************
* memoryEstimation --
**********************************************************************/
size_t
InterpMat::memoryEstimation (
			     void) const
{
  int numEntry = evalGridElementPtr->numElement() * interpStencilPtr->numPoint();
  int numRow = evalGridElementPtr->numElement();
  return matList[0].memoryEstimate(numEntry, numRow) * numMat();
}

/**********************************************************************
 * setup --
 **********************************************************************/
void
InterpMat::setup (
		  void)
{
  if (interpStencilPtr->setupComplete()) {
    for (size_t i = 0; i < numMat(); i++) {
      setupEachMatrix(integrationPtr->outerOperator(i), i);
    }
  }
}

/**********************************************************************
* setupEachMatrix --
* For the potential at the centroid of each panel, a number of interpolation
* stencil points are to be used to compute it. The interpolation
* coefficients are pre-compted and stored in this interpMat matrix.
* The number of coulmns of this matrix should be the total number of grid 
* points. The number of rows should be the number of eval panels.
* Fortunately, each panel only interacts with a small number of grid 
* points, so this matrix is very sparse. I use compressed row format
* to store this matrix.
* In dense format, the ith row is v = b' * F, where b is a column vector 
* that stores the values of each interploation basis function at the 
* centroid of the ith panel. F is a pre-computed small matrix in class
* Stencil. This matrix is only dependent upon interpStencil. It is not
* related to the configuration of panels.
* 
* This function does the following:
* 1) find the centroid of the panel and compute b. 
* 2) compute v = b' * F, which is delegated to interpStencil.
* 3) scatter v to the sparse vector form and put it into interpMat
**********************************************************************/
void
InterpMat::setupEachMatrix (
			    const DifferentialOperator operatorType,
			    const size_t matIndex)
{
  TNT::Vector<double> basisFuncValueAtEvalPoint(interpStencilPtr->numBasis());
  TNT::Vector<double> denseRow(interpStencilPtr->numPoint());
  matList[matIndex] = SpRowMat<double>(evalGridElementPtr->numElement(), 
				       gridPtr->totalNumPointAfterPadded());
  size_t globalElementIndex;

  for (size_t ix=0; ix < gridPtr->numPointX(); ix++) {
    for (size_t iy=0; iy < gridPtr->numPointY(); iy++) {
      for (size_t iz=0; iz < gridPtr->numPointZ(); iz++) {
	int gridPointIndex = gridPtr->transfer3DIndexTo1D(ix, iy, iz);
	size_t numElement = 
	  evalGridElementPtr->numElementMappedToOneGrid(gridPointIndex);
	if (! numElement) continue;

	vector<int> globalIndexList;
	interpStencilPtr->findGlobalIndexOfStencilPoint(ix, iy, iz, *gridPtr, 
							globalIndexList);

	point3D origin(gridPtr->x(ix), gridPtr->y(iy), gridPtr->z(iz));
	for (size_t localElementIndex = 0; localElementIndex < numElement; 
	     localElementIndex++) {
	  globalElementIndex = 
	    evalGridElementPtr->getElementMappedToOneGrid(gridPointIndex, 
							  localElementIndex);
	  const element& evalElement = 
	    evalGridElementPtr->getElement(globalElementIndex);
	  // Note that the coordinate origin used in interp stencil is the 
	  // center point. Therefore, the collocation point should be adjusted 
	  // accordingly
	  point3D collocationPoint = evalElement.centroid() - origin;

	  interpStencilPtr->evaluateBasisFuncAtOnePoint(collocationPoint, 
							evalElement.normal(),
							basisFuncValueAtEvalPoint,
							operatorType);
	  interpStencilPtr->compCoef(basisFuncValueAtEvalPoint, denseRow);
	  copyFullVectorToSpVec(globalIndexList.begin(), globalIndexList.end(), 
				denseRow.begin(), 
				matList[matIndex][globalElementIndex]);
	}
      }
    }
  }
}
