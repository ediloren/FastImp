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

const static char cvsid[] = "$Id: gridElement.cc,v 1.3 2003/02/11 03:12:00 zhzhu Exp $";

#include <fstream>
#include "gridElement.h"
#include <algorithm>
#include "stpwatch.h"

using namespace std;
using namespace pfft;

/**********************************************************************
 * GridElement --
 **********************************************************************/
GridElement::GridElement (
			  ElementList* elementListPtrIn) 
  : elementListPtr(elementListPtrIn)
{
  findScopeOfElementDomain();
  ElementSizeStatistics();
  maxNumNeighbor_ = 0;
}

/**********************************************************************
 * findScopeOfElementDomain --
 **********************************************************************/
void 
GridElement::findScopeOfElementDomain (
				       void)
{
  ElementListConstIterator it = elementListPtr->begin();
  xmin_ = it->centroid().x() - it->boundingSphereRadius();
  xmax_ = it->centroid().x() + it->boundingSphereRadius();
  ymin_ = it->centroid().y() - it->boundingSphereRadius();
  ymax_ = it->centroid().y() + it->boundingSphereRadius();
  zmin_ = it->centroid().z() - it->boundingSphereRadius();
  zmax_ = it->centroid().z() + it->boundingSphereRadius();
  for (; it != elementListPtr->end(); ++it) {
    xmin_ = min(xmin_, it->centroid().x() - it->boundingSphereRadius());
    xmax_ = max(xmax_, it->centroid().x() + it->boundingSphereRadius());
    ymin_ = min(ymin_, it->centroid().y() - it->boundingSphereRadius());
    ymax_ = max(ymax_, it->centroid().y() + it->boundingSphereRadius());
    zmin_ = min(zmin_, it->centroid().z() - it->boundingSphereRadius());
    zmax_ = max(zmax_, it->centroid().z() + it->boundingSphereRadius());
  }
}

/**********************************************************************
 * ElementSizeStatistics --
 **********************************************************************/
void 
GridElement::ElementSizeStatistics (
				    void)
{
  ElementListConstIterator it = elementListPtr->begin();
  double radiusSum = 0.;
  maxElementSize_ = 0.;
  minElementSize_ = it->boundingSphereRadius();
  for (; it != elementListPtr->end(); ++it) {
    radiusSum += it->boundingSphereRadius();
    maxElementSize_ = max(maxElementSize_, it->boundingSphereRadius());
    minElementSize_ = min(minElementSize_, it->boundingSphereRadius());
  }
  averageElementSize_ =  2. * ( radiusSum / elementListPtr->size() );
  maxElementSize_ *= 2.;
  minElementSize_ *= 2.;

}

/**********************************************************************
 * mapElementToGridPoint --
 **********************************************************************/
void 
GridElement::mapElementToGridPoint(
				   const Grid& grid)
{
  elementGridMap.clear();
  elementGridMap.resize(elementListPtr->size());
  gridElementMap.clear();
  gridElementMap.resize(grid.totalNumPointAfterPadded());

  for (size_t index = 0; index < elementListPtr->size(); index++) {
    // find a grid point nearest to the element centroid
    int x = static_cast<int>(floor(((*elementListPtr)[index].centroid().x() - grid.xmin()) / grid.dx() + 0.50));
    int y = static_cast<int>(floor(((*elementListPtr)[index].centroid().y() - grid.ymin()) / grid.dy() + 0.50));
    int z = static_cast<int>(floor(((*elementListPtr)[index].centroid().z() - grid.zmin()) / grid.dz() + 0.50));
    int gridPointIndex = grid.transfer3DIndexTo1D(x, y, z);
    elementGridMap[index] = gridPointIndex;    
    gridElementMap[gridPointIndex].push_back(index);
  }
}

/**********************************************************************
 * countNumNeighbor --
 * Only setup irregular neighbor list because I could afford the 
 * memory used by it. But regular neighbors might be many when grid
 * is coarse. So I can not set them up. But it is cheap to count them.
 **********************************************************************/
void
GridElement::countNumNeighbor (
			       GridElement& guestGridElement,
			       const Grid& grid, 
			       const Stencil& directStencil)
{
  setupIrregularNeighbor(guestGridElement, grid, directStencil);
  countNumIrregularNeighbor();
  countNumRegularNeighbor(guestGridElement, grid, directStencil);
  numNeighbor_ = totalNumIrregularNeighbor_ + totalNumRegularNeighbor_;

  neighborReport();
}

/**********************************************************************
 * setupIrregularNeighbor --
 * Important Note:
 * The irregularNeighbor is to be setup many times in the recusion 
 * of finding the optimal gridSize. So I clear up it every time
 * I want to update a new list.
 **********************************************************************/
void
GridElement::setupIrregularNeighbor (
				     GridElement& guestGridElement,
				     const Grid& grid, 
				     const Stencil& directStencil)
{
  irregularNeighbor.clear();
  if (*this != guestGridElement)
    guestGridElement.irregularNeighbor.clear();  

  findIrregularNeighbour(grid, directStencil, guestGridElement); 
  if (*this != guestGridElement) 
    guestGridElement.findIrregularNeighbour(grid, directStencil, *this); 
 
  removeDuplicateIrregularNeighbor();
  if (*this != guestGridElement) 
    guestGridElement.removeDuplicateIrregularNeighbor();
}

/**********************************************************************
 * findIrregularNeighbour --
 * Purely from mathematical point of view, the criteria to decide if an
 * element is the nearby neighbor to a host one is simple: 
 * the ratio of center-to-center distance over the radius of the host 
 * elment should be larger than a certain threshold. 
 * The error of the pfft algorithm largely depends on this ratio.
 * 
 * From algorithm efficiency point of view, we can not use physical distance
 * as the metric to find neighboring elements. A feasible solution is to
 * use the direct stencil. When element size is not too larger than the grid
 * size, it is safe to say that all neighboring elements of a host element 
 * could be reached through direct stencil using elementGridMap and gridElementMap. 
 *
 * However, when the host element is much larger than gridSize, using direct 
 * stencil alone is not sufficient. This might happen when mesh generator 
 * produces inhomogeneous panels or when the grid is too fine.
 * To remedy this situation, we could set up an "irregular" neighbor list
 * that contains only those neighbors that can not be reached through the
 * direct stencil. For these neighbors, physical distance is measured against
 * element radius to decide if they are eligible. 
 *
 *
 * Important note:
 *
 * This function should be called by both src and eval gridElement object.
 * Each time this function being called, the neighboring element list 
 * in both srcGridElement and evalGridElement will be updated. Some elements
 * will be double counted. That is why I have to remove duplicate later.
 * This conservative double counting strategy will guarantee that 
 * no element will be missed.
 **********************************************************************/
void 
GridElement::findIrregularNeighbour (
				     const Grid& grid,
				     const Stencil& directStencil,
				     GridElement& guest)
{
  if (irregularNeighbor.empty()) {
    irregularNeighbor.resize(elementListPtr->size());
  }
  if ( (*this != guest) && (guest.irregularNeighbor.empty()) ) {
    guest.irregularNeighbor.resize(guest.elementListPtr->size());
  }

  for (int hx = 0; hx < static_cast<int>(grid.numPointX()); hx++) 
    for (int hy = 0; hy < static_cast<int>(grid.numPointY()); hy++)
      for (int hz = 0; hz < static_cast<int>(grid.numPointZ()); hz++) {
	int hostGridPoint = grid.transfer3DIndexTo1D(hx, hy, hz);
	if (! gridElementMap[hostGridPoint].empty()) {
	  vector<int> nearbyNonStencilGridPointList = 
	    findNearbyNonStencilGridPoint(hostGridPoint, hx, hy, hz, grid, directStencil);
	  for (vector<int>::const_iterator 
		 nearbyGridPoint = nearbyNonStencilGridPointList.begin(); 
	       nearbyGridPoint != nearbyNonStencilGridPointList.end(); 
	       ++nearbyGridPoint) {
	    // each nearby grid is the representative to a list of 
	    // potential neighbouring elements
	    for (OneGridElementMapConstIterator 
		   potentialNeighbour = guest.gridElementMap[*nearbyGridPoint].begin();
	         potentialNeighbour != guest.gridElementMap[*nearbyGridPoint].end(); 
	         ++potentialNeighbour) {
	      for (OneGridElementMapConstIterator host = gridElementMap[hostGridPoint].begin();
		   host != gridElementMap[hostGridPoint].end(); ++host) {
		if(AreClose((*elementListPtr)[*host], 
			    (*guest.elementListPtr)[*potentialNeighbour], 
			    grid.separationFactor())) {    
		  // Update neighbouring element list in both host and guest gridElement
		  // This will guarantee that a small element could also "see" a large one,
		  // not just the other way.
		  irregularNeighbor[*host].push_back(*potentialNeighbour);
		  guest.irregularNeighbor[*potentialNeighbour].push_back(*host);
		}
	      }
	    }
	  }
	}
      }
}

/**********************************************************************
 * findNearbyNonStencilGridPoint --
 **********************************************************************/
vector<int>
GridElement::findNearbyNonStencilGridPoint (
					    int hostGridPointIndex,
					    int ihx,
					    int ihy,
					    int ihz,
					    const Grid& grid,
					    const Stencil& directStencil) const
{
  double maxRadius = findMaxBoundingSphereRadius(hostGridPointIndex);
  double gridCubeSize = sqrt(grid.dx()*grid.dx() + grid.dy()*grid.dy() + grid.dz()*grid.dz());
  double maxDistance = maxRadius * (grid.separationFactor() + 2) + gridCubeSize;
  int xDisp = static_cast<int>( ceil(maxDistance / grid.dx()) );
  int yDisp = static_cast<int>( ceil(maxDistance / grid.dy()) );
  int zDisp = static_cast<int>( ceil(maxDistance / grid.dz()) );

  // remove any points outside the grid range
  int xmin = max(ihx - xDisp, 0);
  int xmax = min(ihx + xDisp, static_cast<int>(grid.numPointX()-1));
  int ymin = max(ihy - yDisp, 0);
  int ymax = min(ihy + yDisp, static_cast<int>(grid.numPointY()-1));
  int zmin = max(ihz - zDisp, 0);
  int zmax = min(ihz + zDisp, static_cast<int>(grid.numPointZ()-1));

  vector<int> nearbyNonStencilGridPointList;
  nearbyNonStencilGridPointList.reserve((xmax-xmin) * (ymax-ymin) * (zmax-zmin));
  for (int igx = xmin; igx <= xmax; igx++) {
    for (int igy = ymin; igy <= ymax; igy++) {
      for (int igz = zmin; igz <= zmax; igz++) {
	// These two "if" statements could be switched. But I found that 
	// checking if the gridElementMap is empty first is more efficient.
	int guestGridPointIndex = grid.transfer3DIndexTo1D(igx, igy, igz);
	if (! gridElementMap[guestGridPointIndex].empty()) {
	  if ( (! directStencil.isInsideStencil(ihx-igx, ihy-igy, ihz-igz)) &&
	       (grid.distance(igx, igy, igz, ihx, ihy, ihz) < maxDistance) ) {
	    nearbyNonStencilGridPointList.push_back(guestGridPointIndex);
	  }
	}
      }
    }
  }

  return nearbyNonStencilGridPointList;
}

/**********************************************************************
 * removeDuplicateIrregularNeighbor --
 **********************************************************************/
void
GridElement::removeDuplicateIrregularNeighbor (
					       void)
{
  for (AllNeighborListIterator oneList = irregularNeighbor.begin();
       oneList != irregularNeighbor.end(); ++oneList) {
    sort(oneList->begin(), oneList->end());
    OneNeighborListIterator newEnd = unique(oneList->begin(), oneList->end());
    oneList->erase(newEnd, oneList->end());
  }
}

/**********************************************************************
 * countNumIrregularNeighbor --
 **********************************************************************/
void
GridElement::countNumIrregularNeighbor (
					void)
{
  totalNumIrregularNeighbor_ = 0;
  for (AllNeighborListIterator oneElementNeighborList = irregularNeighbor.begin();
       oneElementNeighborList != irregularNeighbor.end(); ++oneElementNeighborList) {
    totalNumIrregularNeighbor_ += oneElementNeighborList->size();
  }
}

/**********************************************************************
 * setupElementNeighbor --
 **********************************************************************/
void
GridElement::setupElementNeighbor (
				   GridElement& guestGridElement,
				   const Grid& grid, 
				   const Stencil& directStencil)
{
  setupIrregularNeighbor(guestGridElement, grid, directStencil);
  setupRegularNeighbor(guestGridElement, grid, directStencil);
}

/**********************************************************************
 * countNumRegularNeighbor --
 **********************************************************************/
void
GridElement::countNumRegularNeighbor (
				      GridElement& guestGridElement,
				      const Grid& grid, 
				      const Stencil& directStencil)
{
  totalNumRegularNeighbor_ = 0;
  for (int x = 0; x < static_cast<int>(grid.numPointX()); x++) {
    for (int y = 0; y < static_cast<int>(grid.numPointY()); y++) { 
      for (int z = 0; z < static_cast<int>(grid.numPointZ()); z++) {
	int hostGridPointIndex = grid.transfer3DIndexTo1D(x, y, z);
	if (! gridElementMap[hostGridPointIndex].empty()) {
	  size_t numNeighborToOneGrid = 
	    countNumRegularNeighborToOneGridPoint(x, y, z, directStencil, grid, guestGridElement);
	  size_t numSrcMappedToOneGrid = gridElementMap[hostGridPointIndex].size();
	  totalNumRegularNeighbor_ += numSrcMappedToOneGrid * numNeighborToOneGrid;
	}
      }
    }
  }
}

/**********************************************************************
 * countNumRegularNeighborToOneGridPoint --
 **********************************************************************/
size_t
GridElement::countNumRegularNeighborToOneGridPoint (
						    int hx,
						    int hy,
						    int hz,
						    const Stencil& directStencil,
						    const Grid& grid,
						    GridElement& guestGridElement)
{  
  vector<int> stencilPointGlobalIndex;
  directStencil.findGlobalIndexOfStencilPoint(hx, hy, hz, grid, stencilPointGlobalIndex);

  size_t numNeighbor = 0;
  for (vector<int>::const_iterator pointIndex = stencilPointGlobalIndex.begin();
       pointIndex != stencilPointGlobalIndex.end(); ++pointIndex) {
    numNeighbor += guestGridElement.gridElementMap[*pointIndex].size();
  }

  return numNeighbor;
}

/**********************************************************************
 * setupRegularNeighbor --
 **********************************************************************/
void
GridElement::setupRegularNeighbor (
				   GridElement& guestGridElement,
				   const Grid& grid, 
				   const Stencil& directStencil)
{
  if (! regularNeighbor.empty()) {
    regularNeighbor.clear();
  }
  regularNeighbor.resize(elementListPtr->size());

  for (int x = 0; x < static_cast<int>(grid.numPointX()); x++) 
    for (int y = 0; y < static_cast<int>(grid.numPointY()); y++) 
      for (int z = 0; z < static_cast<int>(grid.numPointZ()); z++) {
	int hostGridPointIndex = grid.transfer3DIndexTo1D(x, y, z);
	if (! gridElementMap[hostGridPointIndex].empty()) {
	  // find a list of neighbors to one host grid point through direct stencil
	  // and append this list to the neighbor list of each host element 
	  // associated with this host grid point
	  vector<int> neighborList;
	  findNeighborToOneGridPoint(x, y, z, directStencil, grid, guestGridElement, neighborList);

	  OneGridElementMap hostElementList = gridElementMap[hostGridPointIndex];
	  for (OneGridElementMapConstIterator hostElement = hostElementList.begin();
	       hostElement != hostElementList.end(); ++hostElement) {
	    copy(neighborList.begin(), neighborList.end(), back_inserter(regularNeighbor[*hostElement]) );
	  }
	}
      }
}

/**********************************************************************
 * findNeighborToOneGridPoint --
 **********************************************************************/
void
GridElement::findNeighborToOneGridPoint (
					 int hx,
					 int hy,
					 int hz,
					 const Stencil& directStencil,
					 const Grid& grid,
					 GridElement& guestGridElement,
					 vector<int>& neighborList)
{
  vector<int> stencilPointGlobalIndex;
  directStencil.findGlobalIndexOfStencilPoint(hx, hy, hz, grid, stencilPointGlobalIndex);

  for (vector<int>::const_iterator pointIndex = stencilPointGlobalIndex.begin();
       pointIndex != stencilPointGlobalIndex.end(); ++pointIndex) {
    OneGridElementMap 
      neighborMappedToOneStencilPoint = guestGridElement.gridElementMap[*pointIndex];
    copy(neighborMappedToOneStencilPoint.begin(), neighborMappedToOneStencilPoint.end(), 
	 back_inserter(neighborList));
  }
}

/**********************************************************************
 * findMaxBoundingSphereRadius --
 * find max boundingSphereRadius of all elements associated with one grid
 **********************************************************************/
double
GridElement::findMaxBoundingSphereRadius (
					  int hostGridPointIndex) const
{
  double maxRadius = 0.;
  for (OneGridElementMapConstIterator 
	 elementIndex = gridElementMap[hostGridPointIndex].begin();
       elementIndex != gridElementMap[hostGridPointIndex].end(); ++elementIndex) {
    //  for (ElementListConstIterator element = elementListPtr->begin(); 
    //       element != elementListPtr->end(); ++element) {
    //    maxRadius = max(element->boundingSphereRadius(), maxRadius);
    maxRadius = max((*elementListPtr)[*elementIndex].boundingSphereRadius(), 
		    maxRadius);
  }
  return maxRadius;
}

/**********************************************************************
 * AreClose --
 **********************************************************************/
bool
GridElement::AreClose(
		      element& host,
		      element& guest,
		      double seperationFactor) const
{
  vector3D<double> vec = host.centroid() - guest.centroid();
  double distance = length(vec);
  double threshold = (seperationFactor + 1) * host.boundingSphereRadius() + 
    guest.boundingSphereRadius();
  return (distance < threshold);  
}

/**********************************************************************
 * comparison --
 **********************************************************************/
bool
pfft::operator != (
		   const GridElement& ge1,
		   const GridElement& ge2)
{
  return (ge1.elementListPtr != ge2.elementListPtr);
}

/**********************************************************************
 * memoryUsage --
 * The memory usage by src and eval gridEleemnt is different. This 
 * function try to differentiate them by looking at the internal state
 * variable "numNeighbor_".
 * The irregular neighbor list of both src and eval gridElements have
 * to be established. However, to estimate the memory usage by direct 
 * matrix and set it up, only srcGridElement's regular neighbor list is 
 * needed. 
 * It is assumed before calling this function that the function 
 * "countNumNegihbor" has been called for the object "srcGridElement",
 * but not for the evalGridElement. Setting up regular neighbor is a
 * non-trivial task. It consumes a lot of memory for coarse grid. So I
 * do not want to do it until I have to. And I don't need regular neighbor 
 * list for evalGridElement.
 * So upon calling this function, if "numNeighbor_" is nonzero, it means
 * this might be srcGridElement calling. Otherwise, it is evalGridElement
 * that calls this function.
 **********************************************************************/
size_t
GridElement::memoryUsage(
			 const bool forSrcGridElement) 
{
  size_t memoryUsage = 0;

  // gridElementMap
  // Since there is no duplicate in grid to element map, 
  size_t totalNumMappedElement = elementListPtr->size();
  memoryUsage += sizeof(int) * totalNumMappedElement;

  // elementGridMap
  memoryUsage += sizeof(int) * elementGridMap.size();

  // elementNeighbors only exist for srcGridElement
  if (forSrcGridElement) {
    if (numNeighbor_ == 0) {
      countNumIrregularNeighbor();
      memoryUsage += sizeof(int) * totalNumIrregularNeighbor_;
    } else {
      memoryUsage += sizeof(int) * numNeighbor_;
    }
  }

  return memoryUsage;
}

/**********************************************************************
 * maxNumElementMappedToOneGrid --
 **********************************************************************/
size_t
GridElement::maxNumElementMappedToOneGrid (
					   void) const 
{
  size_t maxNumElement = 0;
  for (AllGridElementMapConstIterator oneGridToElements = gridElementMap.begin();
       oneGridToElements != gridElementMap.end(); ++oneGridToElements) {
    maxNumElement = max(maxNumElement, oneGridToElements->size());
  }
  return maxNumElement;
}

/**********************************************************************
 * countMinMaxNumNeighbor --
 **********************************************************************/
void 
GridElement::countMinMaxNumNeighbor (
				     void)
{
  maxNumNeighbor_ = 0;
  minNumNeighbor_ = elementListPtr->size();
  for (size_t hostIndex = 0; hostIndex < elementListPtr->size(); hostIndex++) {
    size_t numNeighbor = 
      irregularNeighbor[hostIndex].size() + regularNeighbor[hostIndex].size();
    maxNumNeighbor_ = max(maxNumNeighbor_, numNeighbor);
    minNumNeighbor_ = min(minNumNeighbor_, numNeighbor);
  }
}

/**********************************************************************
 * neighborReport --
 **********************************************************************/
void
GridElement::neighborReport (
			     void)
{
#ifdef DEBUG_ELEMENT
  cout << "Number of neighbors not to be reached through the direct stencil := " 
       << totalNumIrregularNeighbor_ << endl
       << "Number of neighbors reached through the direct stencil := " 
       << totalNumRegularNeighbor_ << endl;
#endif
}

/**********************************************************************
 * dumpIrregularNeighborInfo --
 **********************************************************************/
void 
GridElement::dumpIrregularNeighborInfo (void)
{
  std::ofstream fout("irregularInfo.dat");
  size_t hostElementIndex = 0;
  for(AllNeighborList::const_iterator hostElementIter = irregularNeighbor.begin();
      hostElementIter != irregularNeighbor.end(); 
      ++hostElementIter, ++hostElementIndex) {
    fout<<hostElementIndex<<" : ";
    for(OneNeighborList::const_iterator irregularNeighborIter = (*hostElementIter).begin();
	irregularNeighborIter != (*hostElementIter).end(); ++irregularNeighborIter) {
      fout<<*irregularNeighborIter<<" ";
    } 
    fout<<std::endl;
  }
  fout<<std::endl;

  fout.close();
}

/**********************************************************************
 * dumpIrregularNeighborInfo --
 **********************************************************************/
void 
GridElement::dumpRegularNeighborInfo (void)
{
  std::ofstream fout("regularInfo.dat");
  size_t hostElementIndex = 0;
  for(AllNeighborList::const_iterator hostElementIter = regularNeighbor.begin();
      hostElementIter != regularNeighbor.end(); 
      ++hostElementIter, ++hostElementIndex) {
    fout<<hostElementIndex<<" : ";
    for(OneNeighborList::const_iterator regularNeighborIter = (*hostElementIter).begin();
	regularNeighborIter != (*hostElementIter).end(); ++regularNeighborIter) {
      fout<<*regularNeighborIter<<" ";
    } 
    fout<<std::endl;
  }
  fout<<std::endl;

  fout.close();
}
