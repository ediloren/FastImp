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

  const static char cvsid[] = "$Id: meshRing.cc,v 1.4 2002/12/01 15:50:06 zhzhu Exp $";

  ==========================================================================
*/

#include "meshRing.h"
#include <cmath>

using namespace mesh;
using namespace std;

/**********************************************************************
 * MeshRing --
 **********************************************************************/
MeshRing::MeshRing (
		    const surf::Ring* ringPtr, 
		    const bool useNonUniformMesh,
		    const bool useQuadPanel,
		    const bool aspectRatioWarning)
{
  useNonUniformMesh_ = useNonUniformMesh;
  useQuadPanel_ = useQuadPanel;
  aspectRatioWarning_ = aspectRatioWarning;
  origin = ringPtr->origin1();
  X = pfft::normalize(ringPtr->edgePoint() - ringPtr->origin1());
  Z = pfft::normalize(ringPtr->origin2() - ringPtr->origin1());
  Y = crossProd(Z, X);
  Y.normalize();

  mapToCube(ringPtr);
  meshCube(ringPtr);
  mapPolarToCard();
  transferToGlobalCoord();
}

/**********************************************************************
 * mapToCube --
 **********************************************************************/
void
MeshRing::mapToCube (
		     const surf::Ring* ringPtr)
{
  cube.x1 = ringPtr->innerRadius();
  cube.x2 = ringPtr->outerRadius();

  cube.y1 = 0;
  cube.y2 = (1 - ringPtr->arcGap()) * 2 * PI;

  cube.z1 = 0;
  cube.z2 = pfft::length(ringPtr->origin2() - ringPtr->origin1());

  cube.numPanelX = ringPtr->numPanelRadius();
  cube.numPanelY = ringPtr->numPanelArc();
  cube.numPanelZ = ringPtr->numPanelThickness();
}

/**********************************************************************
 * mapPolarToCard --
 **********************************************************************/
void
MeshRing::mapPolarToCard (
			  void)
{
  for (int i=0; i < numNode_; i++) {
    double r = nodeList[i].x();
    double theta = nodeList[i].y();
    nodeList[i].x() = r*cos(theta);
    nodeList[i].y() = r*sin(theta);
  }
}
