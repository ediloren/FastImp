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

  const static char cvsid[] = "$Id: meshWire.cc,v 1.6 2003/03/14 02:07:37 zhzhu Exp $";

  ==========================================================================
*/

#include "meshWire.h"
#include <stdexcept>

using namespace mesh;
using namespace std;

/**********************************************************************
 * MeshWire --
 **********************************************************************/
MeshWire::MeshWire (
		    const surf::Wire* wirePtr, 
		    const bool useNonUniformMesh,
		    const bool useQuadPanel,
		    const bool aspectRatioWarning)
{
  useNonUniformMesh_ = useNonUniformMesh;
  useQuadPanel_ = useQuadPanel;
  aspectRatioWarning_ = aspectRatioWarning;
  origin = wirePtr->leftEndPoint0();
  X = normalize(wirePtr->leftEndPoint1() - wirePtr->leftEndPoint0());
  Z = normalize(wirePtr->leftEndPoint2() - wirePtr->leftEndPoint0());
  Y = crossProd(Z, X);
  Y.normalize();

  mapToCube(wirePtr);
  meshCube(wirePtr);
  transferToGlobalCoord();
}

/**********************************************************************
 * mapToCube --
 **********************************************************************/
void
MeshWire::mapToCube (
		     const surf::Wire* wirePtr)
{
  pfft::vector3D<double> y = 
    normalize(wirePtr->rightEndPoint() - wirePtr->leftEndPoint0());
  double projection = y*Y;
  if ( (projection != 1) && (projection != -1) ){
    cout << "\t meshWire.cc : MeshWire" 
	 << "\t The wire has unorthogonal sides" << endl;
    throw domain_error("Error in meshWire.cc");
  }

  if (projection == -1) {
     origin = wirePtr->rightEndPoint();
  }
  cube.x1 = 0;
  cube.x2 = length(wirePtr->leftEndPoint1() - wirePtr->leftEndPoint0());

  cube.z1 = 0;
  cube.z2 = length(wirePtr->leftEndPoint2() - wirePtr->leftEndPoint0());

  cube.y1 = 0;
  cube.y2 = fabs(length(wirePtr->rightEndPoint() - wirePtr->leftEndPoint0()));

  cube.numPanelX = wirePtr->numPanelBetweenPoint01();
  cube.numPanelY = wirePtr->numPanelLength();
  cube.numPanelZ = wirePtr->numPanelBetweenPoint02();
}


