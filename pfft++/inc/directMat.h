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

  const static char cvsid[] = "$Id: directMat.h,v 1.3 2002/07/18 14:32:19 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _DIRECT_MAT_H_

#define _DIRECT_MAT_H_

#include "spColMat.h"
#include "stencil.h"
#include "interpMat.h"
#include "projectMat.h"
#include "gridElement.h"
#include "grid.h"
#include "grid2grid.h"
#include "g2gUnion.h"

#include "stpwatch.h"
#include "utils.h"
#include "cmat.h"
#include "vec.h"
#include <map>

namespace pfft {

  template <class GreenFunc, class Integration, class KernelType>

  class DirectMat {

  public:

    DirectMat(void) {};
    DirectMat(
	      GridElement* const srcGridElement,
	      GridElement* const evalGridElement, 
	      Integration* const integration,
	      GreenFunc* const greenFunc,
	      Grid* const grid,
	      Stencil* const directStencil,
	      Stencil* const interpStencil,
	      Stencil* const projectStencil,
	      ProjectMat* const projectMat,
	      InterpMat* const interpMat)
      : srcGridElementPtr(srcGridElement), evalGridElementPtr(evalGridElement), 
      gridPtr(grid), directStencilPtr(directStencil), 
      interpStencilPtr(interpStencil), 
      projectStencilPtr(projectStencil), projectMatPtr(projectMat), 
      interpMatPtr(interpMat), integrationPtr(integration), greenFuncPtr(greenFunc)  
    {
      matList = 
	std::vector<SpColMat<KernelType> >(integrationPtr->numIntegral(), 
					   SpColMat<KernelType>(srcGridElementPtr->numElement(), 
								evalGridElementPtr->numElement()));
    }

    size_t memoryEstimation (void) const;
    void setup(void);
    const size_t numMat(void) const { return matList.size(); }
    const SpColMat<KernelType>& operator[] (size_t matIndex) const { 
      return matList[matIndex]; }  
    SpColMat<KernelType>& operator[] (size_t matIndex) { return matList[matIndex]; }  

    void debugReport (void);

  private:

    std::vector<SpColMat<KernelType> > matList;

    GridElement* srcGridElementPtr;
    GridElement* evalGridElementPtr;
    Grid* gridPtr;
    Stencil* directStencilPtr;
    Stencil* interpStencilPtr;
    Stencil* projectStencilPtr;
    ProjectMat* projectMatPtr;
    InterpMat* interpMatPtr;  
    Integration* integrationPtr;
    GreenFunc* greenFuncPtr;

    void fillMat(void);
    void precorrect (Grid2Grid<GreenFunc, KernelType>& g2g);
    void precorrect (G2GUnion<GreenFunc, KernelType>& g2gUnion);
    void precorrectRegularNeighbor (std::vector<size_t>& neighborIndex,
				    int numSrcElement,
				    int hostPoint1DIndex,
				    const GridIndex& hostPoint3DIndex,
				    G2GUnion<GreenFunc, KernelType>& g2gUnion);
    void precorrectRegularNeighbor (std::vector<size_t>& neighborIndex,
				    int numSrcElement,
				    int hostPoint1DIndex,
				    const GridIndex& hostPoint3DIndex,
				    Grid2Grid<GreenFunc, KernelType>& g2g);
    void precorrectIrregularNeighbor (std::vector<size_t>& neighborIndex,
				      int numSrcElement,
				      int hostPoint1DIndex,
				      const GridIndex& hostPoint3DIndex,
				      GreenFunc greenFunc,
				      std::map<GridIndex, TNT::Matrix<KernelType> >& g2gMatList);
    
  };

  /**********************************************************************
   * memoryEstimation --
   * Estimate memory usage before the matrix is actually filled.
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  size_t
  DirectMat<GreenFunc, Integration, KernelType>::memoryEstimation (
								   void) const
  {
    size_t numNonZero = srcGridElementPtr->numNeighbor();
    size_t numCol = srcGridElementPtr->numElement();
    return matList[0].memoryEstimate(numNonZero, numCol) * numMat();
  }

  /**********************************************************************
   * setup --
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  void
  DirectMat<GreenFunc, Integration, KernelType>::setup (
							void)
  {
    //    Grid2Grid<GreenFunc, KernelType> g2g(*projectStencilPtr, *interpStencilPtr, 
    //				 *directStencilPtr, *greenFuncPtr);
    G2GUnion<GreenFunc, KernelType> g2gUnion(*projectStencilPtr, 
					     *interpStencilPtr, 
					     *directStencilPtr, *greenFuncPtr);

    TNT::stopwatch time;  time.start();
    fillMat();
    pfft::timeReport("Time used for filling the direct matrix := ", time.read());

    time.reset();   time.start();
    precorrect(g2gUnion);
    pfft::timeReport("Time used for pre-correcting the direct matrix := ", time.read());
    
    //    debugReport();
  }

  /**********************************************************************
   * fillMat --
   * Compute the interaction between each element and its near neighbors.
   * This is the direct interaction without pre-correction
   *
   * Important note:
   * For the convenience of pre-correction, the neighbors are classified
   * as regular and irregular ones. 
   * For the sake of efficiency, the elements in SpVec is not ordered
   * according to their index as in <index, value> pair. These elements
   * are simply push_back into the SpVec. So they are accessed by the 
   * sequence index. So it is crucial that the order of regular and irregular 
   * neighbor in this function should be the same as in the pre-correction 
   * function. Otherwise, the pre-correction would be imposed upon wrong 
   * entry. 
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  void
  DirectMat<GreenFunc, Integration, KernelType>::fillMat (
							  void)
  {
    size_t numHostElement = srcGridElementPtr->numElement();
    for (size_t hostIndex = 0; hostIndex < numHostElement; hostIndex++) {
      const element& srcElement = srcGridElementPtr->getElement(hostIndex);

      size_t numIrregularNeighbor = 
	srcGridElementPtr->numIrregularNeighbor(hostIndex);
      size_t numRegularNeighbor = 
	srcGridElementPtr->numRegularNeighbor(hostIndex);
      size_t totalNumNeighbor = numIrregularNeighbor + numRegularNeighbor;

      // reserve memory for each column before push_back could save a lot of memory
      for (size_t i = 0; i < integrationPtr->numIntegral(); i++) {
	matList[i][hostIndex].reserve(totalNumNeighbor);
      }

      // first, irregular neighbors
      size_t neighborLocalIndex;
      for (neighborLocalIndex = 0; neighborLocalIndex < numIrregularNeighbor; 
	   neighborLocalIndex++) {
	size_t neighborGlobalIndex = 
	  srcGridElementPtr->getIrregularNeighborGlobalIndex(hostIndex, 
							     neighborLocalIndex);
	const element& evalElement = 
	  evalGridElementPtr->getElement(neighborGlobalIndex);

	(*integrationPtr)(srcElement, evalElement);
	for (size_t i = 0; i < integrationPtr->numIntegral(); i++) {
	  KernelType spElement = integrationPtr->result(i);
	  matList[i].insertElement(neighborGlobalIndex, hostIndex, spElement);
	}
      }

      // second, regular neighbors
      for (neighborLocalIndex = 0; neighborLocalIndex < numRegularNeighbor; 
	   neighborLocalIndex++) {
	int neighborGlobalIndex = 
	  srcGridElementPtr->getRegularNeighborGlobalIndex(hostIndex, 
							   neighborLocalIndex);
	const element& evalElement = 
	  evalGridElementPtr->getElement(neighborGlobalIndex);

	(*integrationPtr)(srcElement, evalElement);
	for (size_t i = 0; i < integrationPtr->numIntegral(); i++) {
	  KernelType spElement = integrationPtr->result(i);
	  matList[i].insertElement(neighborGlobalIndex, hostIndex, spElement);
	}
      }
    }
  }

  /**********************************************************************
   * precorrect --
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  void
  DirectMat<GreenFunc, Integration, KernelType>::
  precorrect (
	      Grid2Grid<GreenFunc, KernelType>& g2g)
  {
    std::map<GridIndex, TNT::Matrix<KernelType> > g2gMatList;

    for (int x = 0; x < static_cast<int>(gridPtr->numPointX()); x++) {
      for (int y = 0; y < static_cast<int>(gridPtr->numPointY()); y++) {
	for (int z = 0; z < static_cast<int>(gridPtr->numPointZ()); z++) {
	  int hostPoint1DIndex = gridPtr->transfer3DIndexTo1D(x, y, z);
	  size_t numSrcElement = 
	    srcGridElementPtr->numElementMappedToOneGrid(hostPoint1DIndex);

	  // if no src element mapped to this point, no need to proceed
	  if (numSrcElement == 0) continue;

	  std::vector<size_t> neighborIndex(numSrcElement, 0);
	  GridIndex hostPoint3DIndex(x, y, z);
	  precorrectIrregularNeighbor(neighborIndex, numSrcElement, hostPoint1DIndex, 
				      hostPoint3DIndex, *greenFuncPtr, g2gMatList);
	  precorrectRegularNeighbor(neighborIndex, numSrcElement, hostPoint1DIndex, 
				    hostPoint3DIndex, g2g);
	}
      }
    }
  }

  /**********************************************************************
   * precorrect --
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  void
  DirectMat<GreenFunc, Integration, KernelType>::
  precorrect (
	      G2GUnion<GreenFunc, KernelType>& g2gUnion)
  {
    std::map<GridIndex, TNT::Matrix<KernelType> > g2gMatList;

    for (int x = 0; x < static_cast<int>(gridPtr->numPointX()); x++) {
      for (int y = 0; y < static_cast<int>(gridPtr->numPointY()); y++) {
	for (int z = 0; z < static_cast<int>(gridPtr->numPointZ()); z++) {
	  int hostPoint1DIndex = gridPtr->transfer3DIndexTo1D(x, y, z);
	  size_t numSrcElement = 
	    srcGridElementPtr->numElementMappedToOneGrid(hostPoint1DIndex);

	  // if no src element mapped to this point, no need to proceed
	  if (numSrcElement == 0) continue;

	  std::vector<size_t> neighborIndex(numSrcElement, 0);
	  GridIndex hostPoint3DIndex(x, y, z);
	  precorrectIrregularNeighbor(neighborIndex, numSrcElement, hostPoint1DIndex, 
				      hostPoint3DIndex, *greenFuncPtr, g2gMatList);
	  precorrectRegularNeighbor(neighborIndex, numSrcElement, hostPoint1DIndex,
				    hostPoint3DIndex, g2gUnion);
	}
      }
    }
  }

  /**********************************************************************
   * precorrectIrregularNeighbor --
   * irregular neighbors that cannot be reached though grid
   *
   * I setup a small database for the g2gMat again here. I use map to store
   * all the small g2gMat and search in this list first everytime I need one.
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  void
  DirectMat<GreenFunc, Integration, KernelType>::
  precorrectIrregularNeighbor (
			       std::vector<size_t>& neighborIndex,
			       int numSrcElement,
			       int hostPoint1DIndex,
			       const GridIndex& hostPoint3DIndex,
			       GreenFunc greenFunc,
			       std::map<GridIndex, TNT::Matrix<KernelType> >& g2gMatList)
  {
    std::vector<TNT::Vector<KernelType> > gp(integrationPtr->numInnerOperator(), 
					     TNT::Vector<KernelType>(interpStencilPtr->numPoint()));
    TNT::Matrix<KernelType> g2gMat(interpStencilPtr->numPoint(), 
				   projectStencilPtr->numPoint());
 
    for (int srcElementLocalIndex = 0; srcElementLocalIndex < numSrcElement; 
	 srcElementLocalIndex++) {
      size_t srcElementGlobalIndex = 
	srcGridElementPtr->getElementMappedToOneGrid(hostPoint1DIndex, 
						     srcElementLocalIndex);

      int numIrregularNeighbor = 
	srcGridElementPtr->numIrregularNeighbor(srcElementGlobalIndex);
      if (numIrregularNeighbor == 0) continue;

      for (int localNeighborIndex = 0; localNeighborIndex < numIrregularNeighbor; 
	   localNeighborIndex++) {
	int evalElementGlobalIndex = 
	  srcGridElementPtr->getIrregularNeighborGlobalIndex(srcElementGlobalIndex, 
							     localNeighborIndex);
        int gridIndex = 
	  evalGridElementPtr->gridMappedToElement(evalElementGlobalIndex);
	GridIndex guestPoint3DIndex = gridPtr->transfer1DIndexTo3D(gridIndex);

	GridIndex relativePosition = hostPoint3DIndex - guestPoint3DIndex;
	if (g2gMatList.find(relativePosition) != g2gMatList.end()) {
	  for (size_t i = 0; i < integrationPtr->numInnerOperator(); i++) {
	    matMultVec(gp[i], g2gMatList[relativePosition], 
		       (*projectMatPtr)[i][srcElementGlobalIndex]);
	  }
	} else { 
	  computeG2GMatrix(directStencilPtr->step(), *projectStencilPtr, 
			   *interpStencilPtr, greenFunc, guestPoint3DIndex, 
			   hostPoint3DIndex, g2gMat);
	  for (size_t i = 0; i < integrationPtr->numInnerOperator(); i++) {
	    matMultVec(gp[i], g2gMat, (*projectMatPtr)[i][srcElementGlobalIndex]);
	  }
	  g2gMatList.insert(std::pair<GridIndex, TNT::Matrix<KernelType> >(relativePosition ,
									   g2gMat));
	}
	
	for (size_t i = 0; i < integrationPtr->numIntegral(); i++) {
	  size_t innerOperatorIndex = 
	    integrationPtr->findInnerOperatorIndex(integrationPtr->innerOperatorOfIntegral(i));
	  size_t outerOperatorIndex = 
	    integrationPtr->findOuterOperatorIndex(integrationPtr->outerOperatorOfIntegral(i));
	  KernelType corrrection;
	  compactDotProd(corrrection, (*interpMatPtr)[outerOperatorIndex][evalElementGlobalIndex],
			 gp[innerOperatorIndex]);
	  matList[i][srcElementGlobalIndex][neighborIndex[srcElementLocalIndex]] -= corrrection;
	}
	neighborIndex[srcElementLocalIndex] ++;
      }	
    }
  }

  /**********************************************************************
   * precorrectRegularNeighbor --
   * regular neighbors that can be reached though grid
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  void
  DirectMat<GreenFunc, Integration, KernelType>::
  precorrectRegularNeighbor (
			     std::vector<size_t>& neighborIndex,
			     int numSrcElement,
			     int hostPoint1DIndex,
			     const GridIndex& hostPoint3DIndex,
			     Grid2Grid<GreenFunc, KernelType>& g2g)
  {
    std::vector<TNT::Vector<KernelType> > gp(integrationPtr->numInnerOperator(), 
					     TNT::Vector<KernelType>(interpStencilPtr->numPoint()));
    std::vector<GridIndex> stencilPointGlobal3DIndex;
    directStencilPtr->findGlobal3DIndexOfStencilPoint(hostPoint3DIndex, *gridPtr, 
						      stencilPointGlobal3DIndex);
    for (size_t pointIndex = 0; pointIndex < stencilPointGlobal3DIndex.size(); 
	 pointIndex++) {
      int stencilPointGlobal1DIndex = 
	gridPtr->transfer3DIndexTo1D(stencilPointGlobal3DIndex[pointIndex]);

      // if no eval element is mapped to this point, no need to proceed
      size_t numEvalElement = 
	evalGridElementPtr->numElementMappedToOneGrid(stencilPointGlobal1DIndex);
      if (numEvalElement == 0) continue;
      
      GridIndex stencilPointLocal3DIndex = 
	stencilPointGlobal3DIndex[pointIndex] - hostPoint3DIndex;
      int stencilPointLocal1DIndex = 
	directStencilPtr->transfer3DIndexTo1D(stencilPointLocal3DIndex);
      const TNT::Matrix<KernelType>& g2gMatRef = g2g.matrix(stencilPointLocal1DIndex);

      for (size_t srcElementLocalIndex = 0; srcElementLocalIndex < numSrcElement; 
	   srcElementLocalIndex++) {
	size_t srcElementGlobalIndex = 
	  srcGridElementPtr->getElementMappedToOneGrid(hostPoint1DIndex, 
						       srcElementLocalIndex);

	for (size_t i = 0; i < integrationPtr->numInnerOperator(); i++) {
	  matMultVec(gp[i], g2gMatRef, (*projectMatPtr)[i][srcElementGlobalIndex]);
	}

	for (size_t evalElementLocalIndex = 0; 
	     evalElementLocalIndex < numEvalElement; evalElementLocalIndex++) {
	  size_t evalElementGlobalIndex = 
	    evalGridElementPtr->getElementMappedToOneGrid(stencilPointGlobal1DIndex, 
							  evalElementLocalIndex);

	  for (size_t i = 0; i < integrationPtr->numIntegral(); i++) {
	    size_t innerOperatorIndex = 
	      integrationPtr->findInnerOperatorIndex(integrationPtr->innerOperatorOfIntegral(i));
	    size_t outerOperatorIndex = 
	      integrationPtr->findOuterOperatorIndex(integrationPtr->outerOperatorOfIntegral(i));
	    KernelType corrrection;
	    compactDotProd(corrrection, (*interpMatPtr)[outerOperatorIndex][evalElementGlobalIndex],
			   gp[innerOperatorIndex]);
	    matList[i][srcElementGlobalIndex][neighborIndex[srcElementLocalIndex]] -= corrrection;
	  }

	  neighborIndex[srcElementLocalIndex] ++;
	}
      }
    }
  }

  /**********************************************************************
   * precorrectRegularNeighbor --
   * regular neighbors that can be reached though grid
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  void
  DirectMat<GreenFunc, Integration, KernelType>::
  precorrectRegularNeighbor (
			     std::vector<size_t>& neighborIndex,
			     int numSrcElement,
			     int hostPoint1DIndex,
			     const GridIndex& hostPoint3DIndex,
			     G2GUnion<GreenFunc, KernelType>& g2gUnion)
  {
    const TNT::Matrix<KernelType>& g2gMatRef = g2gUnion.matrix();
    std::vector<TNT::Vector<KernelType> > gp(integrationPtr->numInnerOperator(), 
					     TNT::Vector<KernelType>(g2gMatRef.num_rows()));
    std::vector<GridIndex> stencilPointGlobal3DIndex;
    directStencilPtr->findGlobal3DIndexOfStencilPoint(hostPoint3DIndex, *gridPtr, 
						      stencilPointGlobal3DIndex);

    for (int srcElementLocalIndex = 0; srcElementLocalIndex < numSrcElement; 
	 srcElementLocalIndex++) {
      size_t srcElementGlobalIndex = 
	srcGridElementPtr->getElementMappedToOneGrid(hostPoint1DIndex, 
						     srcElementLocalIndex);

      for (size_t i = 0; i < integrationPtr->numInnerOperator(); i++) {
	matMultVec(gp[i], g2gMatRef, (*projectMatPtr)[i][srcElementGlobalIndex]);
      }


      for (size_t pointIndex = 0; pointIndex < stencilPointGlobal3DIndex.size(); 
	   pointIndex++) {
	int stencilPointGlobal1DIndex = 
	  gridPtr->transfer3DIndexTo1D(stencilPointGlobal3DIndex[pointIndex]);
	
	// if no eval element is mapped to this point, no need to proceed
	size_t numEvalElement = 
	  evalGridElementPtr->numElementMappedToOneGrid(stencilPointGlobal1DIndex);
	if (numEvalElement == 0) continue;
	
	GridIndex stencilPointLocal3DIndex = 
	  stencilPointGlobal3DIndex[pointIndex] - hostPoint3DIndex;
	int stencilPointLocal1DIndex = 
	  directStencilPtr->transfer3DIndexTo1D(stencilPointLocal3DIndex);
	const std::vector<int>& pointToG2GRowIndexMap = 
	  g2gUnion.pointToG2GRowMap(stencilPointLocal1DIndex);

	for (size_t evalElementLocalIndex = 0; 
	     evalElementLocalIndex < numEvalElement; evalElementLocalIndex++) {
	  size_t evalElementGlobalIndex = 
	    evalGridElementPtr->getElementMappedToOneGrid(stencilPointGlobal1DIndex, 
							  evalElementLocalIndex);

	  for (size_t i = 0; i < integrationPtr->numIntegral(); i++) {
	    size_t innerOperatorIndex = 
	      integrationPtr->findInnerOperatorIndex(integrationPtr->innerOperatorOfIntegral(i));
	    size_t outerOperatorIndex = 
	      integrationPtr->findOuterOperatorIndex(integrationPtr->outerOperatorOfIntegral(i));
	    KernelType corrrection;	    
	    selectDotProd(corrrection, 
			  (*interpMatPtr)[outerOperatorIndex][evalElementGlobalIndex], 
			  gp[innerOperatorIndex], pointToG2GRowIndexMap);

	    matList[i][srcElementGlobalIndex][neighborIndex[srcElementLocalIndex]] -= corrrection;
	  }

	  neighborIndex[srcElementLocalIndex] ++;
	}
      }
    }
  }

  /**********************************************************************
   * debugReport --
   **********************************************************************/
  template <class GreenFunc, class Integration, class KernelType>
  void 
  DirectMat<GreenFunc, Integration, KernelType>::
  debugReport (
	       void)
  {
#ifdef DEBUG_ELEMENT
    int totalNonZeros = 0;
    for (size_t ii=0; ii<matList.size(); ++ii) {
      totalNonZeros += matList[ii].numNonZero();
    }
    std::cout << "\tTotal Non Zero of Direct Matrix: " 
	      << totalNonZeros<<std::endl;
#endif
  }
  
} // namespace pfft

#endif

