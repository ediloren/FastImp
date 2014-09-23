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

  const static char cvsid[] = "$Id: meshSpiral.cc,v 1.4 2002/12/01 15:50:06 zhzhu Exp $";

  ==========================================================================
*/

#include "meshSpiral.h"
#include <cmath>

using namespace mesh;
using namespace std;

/**********************************************************************
 * MeshSpiral --
 **********************************************************************/
MeshSpiral::MeshSpiral (
			const surf::Spiral* spiralPtr,
			const bool useNonUniformMesh,
			const bool useQuadPanel,
			const bool aspectRatioWarning)
{
  useNonUniformMesh_ = useNonUniformMesh;
  useQuadPanel_ = useQuadPanel;
  aspectRatioWarning_ = aspectRatioWarning;
  origin = spiralPtr->origin1();
  X = pfft::normalize(spiralPtr->edgePoint() - spiralPtr->origin1());
  Z = pfft::normalize(spiralPtr->origin2() - spiralPtr->origin1());
  Y = crossProd(Z, X);
  Y.normalize();

  mapToCube(spiralPtr);
  meshCube(spiralPtr);
  mapPolarToCard(spiralPtr);
  transferToGlobalCoord();
}

/**********************************************************************
 * mapToCube --
 **********************************************************************/
void
MeshSpiral::mapToCube (
		       const surf::Spiral* spiralPtr)
{
  cube.x1 = spiralPtr->innerRadius();
  cube.x2 = spiralPtr->outerRadius();

  cube.y1 = 0;
  cube.y2 = spiralPtr->numRound() * 2 * PI;

  cube.z1 = 0;
  cube.z2 = pfft::length(spiralPtr->origin2() - spiralPtr->origin1());

  cube.numPanelX = spiralPtr->numPanelRadius();
  cube.numPanelY = spiralPtr->numPanelArc();
  cube.numPanelZ = spiralPtr->numPanelThickness();
}

/**********************************************************************
 * mapPolarToCard --
 **********************************************************************/
void
MeshSpiral::mapPolarToCard (
			    const surf::Spiral* spiralPtr)
{
  double ramp = (spiralPtr->r3() - spiralPtr->innerRadius()) / (2 * PI);
  
  for (int i=0; i < numNode_; i++) {
    double theta = nodeList[i].y();
    double r = nodeList[i].x();
    r += theta * ramp;
    nodeList[i].x() = r*cos(theta);
    nodeList[i].y() = r*sin(theta);
  }
}
