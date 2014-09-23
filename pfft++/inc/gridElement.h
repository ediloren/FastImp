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

  const static char cvsid[] = "$Id: gridElement.h,v 1.2 2003/02/11 02:58:03 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _GRID_ELEMENT_H_
#define _GRID_ELEMENT_H_

//#include <set>
#include <vector>
#include "element.h"
#include "grid.h"
#include "stencil.h"

namespace pfft {

class GridElement {

public:
  friend bool operator != (const GridElement&, const GridElement&);  

  typedef std::vector<element> ElementList;
  typedef ElementList::const_iterator ElementListConstIterator;
  typedef ElementList::iterator ElementListIterator;
  
  // A set of element indices of the opposite type to one host element
  typedef std::vector<int> OneNeighborList;
  typedef std::vector<OneNeighborList> AllNeighborList;
  typedef OneNeighborList::const_iterator OneNeighborListConstIterator;
  typedef OneNeighborList::iterator OneNeighborListIterator;
  typedef AllNeighborList::const_iterator AllNeighborListConstIterator;
  typedef AllNeighborList::iterator AllNeighborListIterator;

  // From one gridIndex to a list of elements 
  // I did not use stl class map (using class GridIndex as key) here. 
  // Storing data in map involves automatic sorting. And I need random access.
  // So vector seems to be better than map. 
  typedef std::vector<int> OneGridElementMap;
  typedef std::vector<OneGridElementMap> AllGridElementMap;
  typedef OneGridElementMap::const_iterator OneGridElementMapConstIterator;
  typedef OneGridElementMap::iterator OneGridElementMapIterator;
  typedef AllGridElementMap::const_iterator AllGridElementMapConstIterator;
  typedef AllGridElementMap::iterator AllGridElementMapIterator;

  // From one elementIndex to one gridIndex. I did not use map for the same
  // reason as above
  typedef std::vector<int> ElementGridMap;
  typedef ElementGridMap::const_iterator ElementGridMapConstIterator;
  typedef ElementGridMap::iterator ElementGridMapIterator;

  GridElement(ElementList*);
  GridElement(void) { };

  void mapElementToGridPoint(const Grid&); 
  void countNumNeighbor (GridElement&, const Grid&, const Stencil&);
  void countMinMaxNumNeighbor (void);
  size_t minNumNeighbor(void) {
    if (!maxNumNeighbor_) countMinMaxNumNeighbor();
    return minNumNeighbor_;
  }
  size_t maxNumNeighbor(void) {
    if (!maxNumNeighbor_) countMinMaxNumNeighbor();
    return maxNumNeighbor_;
  }
  void setupElementNeighbor(GridElement&, const Grid&, const Stencil&);
  size_t memoryUsage(const bool forSrcGridElement=true);
  double xmin(void) const { return xmin_; }
  double xmax(void) const { return xmax_; }
  double ymin(void) const { return ymin_; }
  double ymax(void) const { return ymax_; }
  double zmin(void) const { return zmin_; }
  double zmax(void) const { return zmax_; }
  double averageElementSize(void) const {return averageElementSize_; }
  double maxElementSize(void) const {return maxElementSize_; }
  double minElementSize(void) const {return minElementSize_; }
  size_t maxNumElementMappedToOneGrid (void) const ;

  const size_t numElement(void) const { return elementListPtr->size(); }
  const element& getElement(size_t index) const { return (*elementListPtr)[index]; }
  size_t numElementMappedToOneGrid(size_t gridIndex) { 
    return gridElementMap[gridIndex].size(); }
  size_t getElementMappedToOneGrid(size_t gridIndex, size_t localElementIndex) { 
    return gridElementMap[gridIndex][localElementIndex]; }
  size_t gridMappedToElement(size_t elementIndex) { 
    return elementGridMap[elementIndex]; }

  size_t numNeighbor(void) const { return numNeighbor_; }
  size_t numNeighbor(size_t hostIndex) const { 
    return numRegularNeighbor(hostIndex) + numIrregularNeighbor(hostIndex); }
  size_t numRegularNeighbor(size_t hostIndex) const { 
    return regularNeighbor[hostIndex].size(); }
  size_t numIrregularNeighbor(size_t hostIndex) const { 
    return irregularNeighbor[hostIndex].size(); }
  const size_t getIrregularNeighborGlobalIndex(size_t hostIndex, 
					       size_t neighborLocalIndex) const { 
    return irregularNeighbor[hostIndex][neighborLocalIndex]; }
  const size_t getRegularNeighborGlobalIndex(size_t hostIndex, 
					     size_t neighborLocalIndex) const { 
    return regularNeighbor[hostIndex][neighborLocalIndex]; }

private:
  ElementList* elementListPtr; 
  AllNeighborList irregularNeighbor; 
  AllNeighborList regularNeighbor; 
  AllGridElementMap gridElementMap;
  ElementGridMap elementGridMap; 

  double xmin_;
  double xmax_;
  double ymin_;
  double ymax_;
  double zmin_;
  double zmax_;
  double averageElementSize_;
  double maxElementSize_;
  double minElementSize_;
  size_t numNeighbor_;
  size_t totalNumIrregularNeighbor_;
  size_t totalNumRegularNeighbor_;
  size_t maxNumNeighbor_;
  size_t minNumNeighbor_;

  void findScopeOfElementDomain(void);
  void ElementSizeStatistics(void);
  double findMaxBoundingSphereRadius (int) const;
  void findIrregularNeighbour(const Grid&, const Stencil&, GridElement&);
  std::vector<int> findIrregularNeighborToOneGrid (int, int, int, const Grid&, 
					      const Stencil&, GridElement&) const;
  void removeDuplicateIrregularNeighbor(void);
  bool AreClose(element&, element&, double) const;
  std::vector<int> findNearbyNonStencilGridPoint (int, int, int, int, 
						  const Grid&, const Stencil&) const;
  void countNumIrregularNeighbor (void);
  void countNumRegularNeighbor (GridElement&, const Grid&, const Stencil&);
  void setupIrregularNeighbor(GridElement&, const Grid&, const Stencil&);
  void setupRegularNeighbor(GridElement&, const Grid&, const Stencil&);
  size_t countNumRegularNeighborToOneGridPoint(int, int, int, const Stencil&, 
					       const Grid&, GridElement&);
  void findNeighborToOneGridPoint(int, int, int, const Stencil&, const Grid&, 
				  GridElement&, std::vector<int>&);
  void neighborReport(void);

  void dumpIrregularNeighborInfo(void);
  void dumpRegularNeighborInfo(void);
};

} // namespace pfft


#endif


