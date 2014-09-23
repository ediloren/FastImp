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

  const static char cvsid[] = "$Id: meshCond.cc,v 1.17 2003/05/13 21:10:50 zhzhu Exp $";

  ==========================================================================
*/

#include "meshCond.h"
#include "condInfo.h"

using namespace std;
using namespace mesh;
using namespace surf;

/**********************************************************************
 * meshCube --
 **********************************************************************/
void
MeshCond::meshCube(
		   const CondInfo* condInfoPtr)
{
  numNode_ = 2 * (cube.numPanelX + cube.numPanelZ) * (cube.numPanelY + 1) + 
    2 * (cube.numPanelX - 1) * ( cube.numPanelZ - 1);
  numNonContactPanel_ = 0;

  if (useNonUniformMesh_) {
    findNonUniformNode();
  } else {
    findUniformNode();
  }
  findPanel(condInfoPtr);
}

/**********************************************************************
 * findPanel --
 **********************************************************************/
void
MeshCond::findPanel (
		     const CondInfo* condInfoPtr)
{
  int i, ix, iy, iz;
  vector<int> nodeIndex(4);
  int numNodeOnSides, numNodeOnEnd;
  int offset1, offset2;

  numPanel_ = 2 * (cube.numPanelX + cube.numPanelZ) * cube.numPanelY + 
    2 * cube.numPanelX * cube.numPanelZ;
  if ((! useQuadPanel_) && (cube.numPanelY > 4)) {
    numPanel_ += 2 * (cube.numPanelX + cube.numPanelZ) * (cube.numPanelY-4);
  }
  panelList.reserve(numPanel_);

  /* panels on the sides */
  int numPanelX = cube.numPanelX;
  int numPanelY = cube.numPanelY;
  int numPanelZ = cube.numPanelZ;
  int numPanelXZ = 2 * (numPanelX + numPanelZ);
  for (iy=0; iy<numPanelY; iy++) {
    for (i=0; i<numPanelXZ; i++) {
      nodeIndex[0] = iy*numPanelXZ + i;
      nodeIndex[1] = nodeIndex[0] + 1;
      if (i == numPanelXZ-1) {
	nodeIndex[1] = iy*numPanelXZ;
      }
      nodeIndex[2] = nodeIndex[1] + numPanelXZ;
      nodeIndex[3] = nodeIndex[0] + numPanelXZ;

      PanelType panelType;
      if ((iy < 2) || (iy > numPanelY-3) ) {
	panelType = BUFFER;
      } else {
	panelType = NON_CONTACT;
      }

      if ((! useQuadPanel_) && (panelType != BUFFER)) {
	// trianglation on non-buffer area
	panelList.push_back(Panel(nodeIndex[0], nodeIndex[1], nodeIndex[2], 
				  panelType));
	panelList.push_back(Panel(nodeIndex[0], nodeIndex[2], nodeIndex[3], 
				  panelType));
	numNonContactPanel_ += 2;
      } else {
	panelList.push_back(Panel(nodeIndex, panelType));
	numNonContactPanel_++;
      }
    }
  }

  /* panels on ends */
  numNodeOnSides =  numPanelXZ * (numPanelY + 1);
  numNodeOnEnd =  (numPanelX-1) * (numPanelZ-1);
  offset1 = numNodeOnSides - 1;
  offset2 =  numPanelXZ * numPanelY;

  /* right end */
  for (iz=0; iz < numPanelZ; iz++) {
    for (ix=0; ix < numPanelX; ix++) {
      if (iz == 0) {
	nodeIndex[0] = offset2 + ix;
	nodeIndex[1] = nodeIndex[0] + 1;
	if (numPanelZ == 1) {
	  nodeIndex[3] = offset1 - ix;
	  nodeIndex[2] = nodeIndex[3] - 1;
	} else {
	  nodeIndex[3] = offset1 + ix;
	  nodeIndex[2] = nodeIndex[3] + 1;
	}
	if (ix == numPanelX-1) {
	  nodeIndex[2] = offset2 + numPanelX + 1;
	}
      } else if (iz == (numPanelZ-1)) {
 	nodeIndex[0] = offset1 + (numPanelX-1)*(iz-1) + ix;
	nodeIndex[1] = nodeIndex[0] + 1;
	if (ix==0) nodeIndex[0] = offset1 - (iz - 1);
	if (ix == (numPanelX-1)) nodeIndex[1] = offset2 + numPanelX + iz;
	nodeIndex[2] = offset2 + numPanelX + numPanelZ + numPanelX-1-ix;
 	nodeIndex[3] = nodeIndex[2] + 1;
      } else {
	nodeIndex[0] = offset1 + (iz-1)*(numPanelX-1) + ix;
 	nodeIndex[1] = nodeIndex[0] + 1;
 	nodeIndex[2] = nodeIndex[1] + (numPanelX - 1);
 	nodeIndex[3] = nodeIndex[2] - 1;
	if (ix==0) {
	  nodeIndex[0] = offset1 - (iz - 1);
	  nodeIndex[3] = nodeIndex[0] - 1;
	}
	if (ix == numPanelX-1) {
	  nodeIndex[1] = offset2 + numPanelX + iz;
	  nodeIndex[2] = nodeIndex[1] + 1;
	}
      }

      PanelType type;
      if (condInfoPtr->condType() != CondInfo::GROUND) {
	type = RIGHT_CONTACT;
      } else {
	type = NON_CONTACT;
      }
      panelList.push_back(Panel(nodeIndex, type));
      numNonContactPanel_++;
    }
  }

  /* left end */
  offset1 = numNodeOnSides + numNodeOnEnd;
  for (iz=0; iz < numPanelZ; iz++) {
    for (ix=0; ix < numPanelX; ix++) {
      if (iz == 0) {
	nodeIndex[3] = ix;
	nodeIndex[2] = nodeIndex[3] + 1;
	if (numPanelZ == 1) {
	  nodeIndex[0] = numPanelXZ - 1 - ix;
	  nodeIndex[1] = nodeIndex[0] - 1;
	} else {
	  nodeIndex[0] = offset1 + ix - 1;
	  nodeIndex[1] = nodeIndex[0] + 1;
	  if (ix==0) {
	    nodeIndex[0] = numPanelXZ - 1;
	  }
	  if (ix == (numPanelX-1)) {
	    nodeIndex[1] = numPanelX + 1;
	  }
	}
      } else if (iz == (numPanelZ-1)) {
 	nodeIndex[0] = numPanelXZ - 1 - iz - ix;
	nodeIndex[1] = nodeIndex[0] - 1;
	nodeIndex[2] = offset1 + (iz-1)*(numPanelX-1) + ix;
 	nodeIndex[3] = nodeIndex[2] - 1;
	if (ix==0) {
	  nodeIndex[3] = nodeIndex[0] + 1;
	}
	if (ix == (numPanelX-1)) {
	  nodeIndex[2] = nodeIndex[1] - 1;
	}

      } else {
	nodeIndex[0] = offset1 + iz*(numPanelX-1) + ix - 1;
	nodeIndex[1] = nodeIndex[0] + 1;
	nodeIndex[2] = nodeIndex[1] - (numPanelX - 1);
	nodeIndex[3] = nodeIndex[2] - 1;
	if (ix==0) {
	  nodeIndex[0] = numPanelXZ - 1 - iz;
	  nodeIndex[3] = nodeIndex[0] + 1;
	}
	if (ix == numPanelX-1) {
	  nodeIndex[1] = numPanelX + iz + 1;
	  nodeIndex[2] = nodeIndex[1] - 1;
	}
      }

      PanelType type;
      if (condInfoPtr->condType() != CondInfo::GROUND) {
	type = LEFT_CONTACT;
      } else {
	type = NON_CONTACT;
      }
      panelList.push_back(Panel(nodeIndex, type));
      numNonContactPanel_++;
    }
  }

  //  assert(panelList.size() == numPanel_);
}

/**********************************************************************
 * findUniformNode --
 **********************************************************************/
void
MeshCond::findUniformNode (
			   void)
{
  int ix, iy, iz;
  int nodeIndex;
  double x, y, z;
  int offset1, offset2;
  int numNodeOnSides, numNodeOnEnd;

  nodeList.resize(numNode_);

  int numNodeX = cube.numPanelX + 1;
  int numNodeY = cube.numPanelY + 1;
  int numNodeZ = cube.numPanelZ + 1;
  int numNodeXZ = 2 * (numNodeX + numNodeZ - 2);
  double hx = (cube.x2-cube.x1) / cube.numPanelX;
  double hy = (cube.y2-cube.y1) / cube.numPanelY;
  double hz = (cube.z2-cube.z1) / cube.numPanelZ;
  minUniformStep = min(hx, min(hy, hz));

  /* nodes on the sides */
  for (iy=0; iy < numNodeY; iy++) {
    y = cube.y1 + iy*hy;

    /* bottom */
    z = 0;
    for (ix=0; ix<numNodeX-1; ix++) {
      x = cube.x1 + ix*hx;
      nodeIndex = iy*numNodeXZ + ix;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }

    /* back */
    x = cube.x2;
    for (iz=0; iz<numNodeZ-1; iz++) {
      z = cube.z1 + iz*hz;
      nodeIndex = iy*numNodeXZ + numNodeX-1 + iz;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }

    /* top */
    z = cube.z2;
    for (ix=0; ix<numNodeX-1; ix++) {
      x = cube.x2 - ix*hx;
      nodeIndex = iy*numNodeXZ + numNodeX-1 + numNodeZ-1 + ix;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }

    /* front */
    x = cube.x1;
    for (iz=0; iz<numNodeZ-1; iz++) {
      z = cube.z2 - iz*hz;
      nodeIndex = iy*numNodeXZ + 2*(numNodeX-1) + numNodeZ-1 + iz;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }
  }

  /* nodes on the ends */
  numNodeOnSides =  numNodeXZ * numNodeY;
  numNodeOnEnd =  (numNodeX-2) * (numNodeZ-2);
  offset1 = numNodeOnSides;
  offset2 = offset1 + numNodeOnEnd;

  /* right end */
  y = cube.y2;
  for (iz=0; iz<numNodeZ-2; iz++) {
    for (ix=0; ix<numNodeX-2; ix++) {
      nodeIndex = offset1 + iz*(numNodeX-2) + ix;
      x = cube.x1 + (ix+1)*hx;
      z = cube.z1 + (iz+1)*hz;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }
  }

  /* left end */
  y = cube.y1;
  for (iz=0; iz<numNodeZ-2; iz++) {
    for (ix=0; ix<numNodeX-2; ix++) {
      nodeIndex = offset2 + iz*(numNodeX-2) + ix;
      x = cube.x1 + (ix+1)*hx;
      z = cube.z1 + (iz+1)*hz;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }
  }
  assert(nodeList.size() == numNode_);
}

/**********************************************************************
 * findNonUniformNode --
 **********************************************************************/
void
MeshCond::findNonUniformNode (
			      void)
{
  int ix, iy, iz;
  int nodeIndex;
  double x, y, z;
  int offset1, offset2;
  int numNodeOnSides, numNodeOnEnd;

  nodeList.resize(numNode_);

  int numNodeX = cube.numPanelX + 1;
  int numNodeY = cube.numPanelY + 1;
  int numNodeZ = cube.numPanelZ + 1;
  int numNodeXZ = 2 * (numNodeX + numNodeZ - 2);
  findNonUniformStep(cube.x2-cube.x1, cube.numPanelX, hx);
  findNonUniformStep(cube.z2-cube.z1, cube.numPanelZ, hz);

  const int numBufferLayer = 2;

#ifdef USE_FIXED_BUFFER_SIZE
  double buffer_hy = findBuferStepSize();
  double bufferSize = 2 * numBufferLayer * buffer_hy;
  double hy = (cube.y2-cube.y1-bufferSize) /(cube.numPanelY -2*numBufferLayer);
  cout << endl << "\t buffer_hy = " << buffer_hy << "\t hy = " << hy << endl;
#else
  double hy = (cube.y2-cube.y1) / cube.numPanelY;
  double buffer_hy = hy;
#endif

  if (aspectRatioWarning_) {
    checkAspectRatio(hx, hz, hy);
  }

  /* nodes on the sides */
  y = cube.y1;
  for (iy=0; iy < numNodeY; iy++) {
    if (iy == 0) {
      y = y;
    } else if ( ((iy > 0) && (iy <= numBufferLayer)) || 
		(iy >= numNodeY-numBufferLayer) ) {
      y += buffer_hy;
    } else {
      y += hy;
    }

    /* bottom */
    z = cube.z1;
    for (ix=0; ix<numNodeX-1; ix++) {
      if (ix==0) {
	x = cube.x1;
      } else {
	x += hx[ix-1];
      }
      nodeIndex = iy*numNodeXZ + ix;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }

    /* back */
    x = cube.x2;
    for (iz=0; iz<numNodeZ-1; iz++) {
      if (iz==0) {
	z = cube.z1;
      } else {
	z += hz[iz-1];
      }
      nodeIndex = iy*numNodeXZ + numNodeX-1 + iz;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }

    /* top */
    z = cube.z2;
    for (ix=0; ix<numNodeX-1; ix++) {
      if (ix==0) {
	x = cube.x2;
      } else {
	x -= hx[ix-1];
      }
      nodeIndex = iy*numNodeXZ + numNodeX-1 + numNodeZ-1 + ix;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }

    /* front */
    x = cube.x1;
    for (iz=0; iz<numNodeZ-1; iz++) {
      if (iz==0) {
	z = cube.z2;
      } else {
	z -= hz[iz-1];
      }
      nodeIndex = iy*numNodeXZ + 2*(numNodeX-1) + numNodeZ-1 + iz;
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }
  }

  /* nodes on the ends */
  numNodeOnSides =  numNodeXZ * numNodeY;
  numNodeOnEnd =  (numNodeX-2) * (numNodeZ-2);
  offset1 = numNodeOnSides;
  offset2 = offset1 + numNodeOnEnd;

  /* right end */
  y = cube.y2;
  z = cube.z1;
  for (iz=1; iz < numNodeZ-1; iz++) {
    z += hz[iz-1];
    x = cube.x1;
    for (ix=1; ix < numNodeX-1; ix++) {
      x += hx[ix-1];
      nodeIndex = offset1 + (iz-1)*(numNodeX-2) + (ix-1);
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }
  }

  /* left end */
  y = cube.y1;
  z = cube.z1;
  for (iz=1; iz < numNodeZ-1; iz++) {
    z += hz[iz-1];
    x = cube.x1;
    for (ix=1; ix < numNodeX-1; ix++) {
      x += hx[ix-1];
      nodeIndex = offset2 + (iz-1)*(numNodeX-2) + (ix-1);
      nodeList[nodeIndex] = pfft::point3D(x, y, z);
    }
  }

  assert(nodeList.size() == numNode_);
}

/**********************************************************************
 * findBuferStepSize --
 **********************************************************************/
double
MeshCond::findBuferStepSize (
			     void)
{
  return min(hx.min(), hz.min());
  //  return max(hx.max(), hz.max());
}

/**********************************************************************
 * findNonUniformStep --
 * 1: 1
 * 2: 1 1
 * 3: 1 r1 1
 * 4: 1 r1 r1 1
 * 5: 1 r1 r1^2 r1 1
 * 6: 1 r1 r1^2 r1^2 r1 1
 * 7: 1 r1 r1^2 r1^2*r2 r1^2 r1 1
 * 8: 1 r1 r1^2 r1^2*r2 r1^2*r2 r1^2 r1 1
 * 9: 1 r1 r1^2 r1^2*r2 r1^2*r2*r3 r1^2*r2 r1^2 r1 1
 * 10:1 r1 r1^2 r1^2*r2 r1^2*r2*r3 r1^2*r2*r3 r1^2*r2 r1^2 r1 1
 * 11:1 r1 r1^2 r1^2*r2 r1^2*r2*r3 r1^2*r2*r3 r1^2*r2*r3 r1^2*r2 r1^2 r1 1
 *more:1 r1 r1^2 r1^2*r2 r1^2*r2*r3 ......... r1^2*r2*r3 r1^2*r2 r1^2 r1 1
 **********************************************************************/
void 
MeshCond::findNonUniformStep(
			     const double length,
			     const int numStep,
			     TNT::Vector<double>& h)
{
  TNT::Vector<double> stepSize(numStep);
#ifdef LARGE_STEP_RATIO
  const double r1 = 1.5;
  const double r2 = 3;
  const double r3 = 5;
#else
  const double r1 = 1.3;
  const double r2 = 2;
  const double r3 = 3;
#endif

  switch (numStep) {
  case 1:
  case 2:
    stepSize[0] = 1; 
    break;
  case 3:
  case 4:
    stepSize[0] = 1; 
    stepSize[1] = r1; 
    break;
  case 5:
  case 6:
    stepSize[0] = 1; 
    stepSize[1] = r1; 
    stepSize[2] = r1 * r1; 
    break;
  case 7:
  case 8:
    stepSize[0] = 1; 
    stepSize[1] = r1; 
    stepSize[2] = r1 * r1; 
    stepSize[3] = r1 * r1 * r2; 
    break;
  case 9:
  case 10:
    stepSize[0] = 1; 
    stepSize[1] = r1; 
    stepSize[2] = r1 * r1; 
    stepSize[3] = r1 * r1 * r2; 
    stepSize[4] = r1 * r1 * r2 * r3; 
    break;
  default:
    stepSize[0] = 1; 
    stepSize[1] = r1; 
    stepSize[2] = r1 * r1; 
    stepSize[3] = r1 * r1 * r2; 
    stepSize[4] = r1 * r1 * r2 * r3; 
    stepSize[5] = stepSize[4]; 
    break;
  }

  if (numStep%2 ==0) {
    // even numberStep

    if (numStep > 12) {
      for (int i = 6; i <= numStep/2-1; i++) {
	stepSize[i] = stepSize[5];
      }
    }
    
    // mirror the first half to the second halt
    for (int i = numStep/2; i < numStep; i++) {
      stepSize[i] = stepSize[numStep-i-1];
    }

  } else {
    //odd numStep
    int middlePoint = static_cast<int>(floor(numStep/2.));
    
    if (numStep > 12) {
      for (int i = 6; i <= middlePoint; i++) {
	stepSize[i] = stepSize[5];
      }
    }
    
    for (int i = middlePoint+1; i < numStep; i++) {
      stepSize[i] = stepSize[numStep-i-1];
    }
  }

  double resolution = length / sum(stepSize);
  h.newsize(numStep);
  h = stepSize * resolution;
}

/**********************************************************************
 * checkAspectRatio --
 **********************************************************************/
void 
MeshCond::checkAspectRatio(
			   const TNT::Vector<double>& hx, 
			   const TNT::Vector<double>& hz, 
			   const double hy)
{
  // these two hard-coded numbers intend to guarantee that the quad panel 
  // is close to long narrow shape. Experiments show that accuracy tend
  // to be good if the shape is long and narrow, particularly for the panels
  // close to the edges. This constrain is relaxed somewhat for the panels
  // far away from edges. 
  const double threshold1 = 10; 
  const double threshold2 = 2;
  double min_hx = hx.min();
  double min_hz = hz.min();
  double max_hx = hx.max();
  double max_hz = hz.max();

  if (hy / min_hx >= threshold1) {
    cout << endl << "\t Aspect ratio := " << hy / min_hx << endl;
    if (useQuadPanel_) {
      surf::warningMessage("MeshCond::checkAspectRatio",
			 "minimum step along width is too small");
    } else {
      surf::errorMessage("MeshCond::checkAspectRatio",
			 "minimum step along width is too small");
    }
  }

  if (hy / min_hz >= threshold1) {
    cout << endl << "\t Aspect ratio := " << hy / min_hz << endl;
    //    cout << endl << "\t hy := " << hy << "\t min_hz := " << min_hz << endl;  
    if (useQuadPanel_) {
      surf::warningMessage("MeshCond::checkAspectRatio",
			   "minimum step along thickness is too small");
    } else {
      surf::errorMessage("MeshCond::checkAspectRatio",
			 "minimum step along thickness is too small");
    }
  }

  if (max_hx / hy >= threshold2) {
    cout << endl << "\t Aspect ratio := " << max_hx / hy << endl;
    if (useQuadPanel_) {
      surf::errorMessage("MeshCond::checkAspectRatio",
			   "maximum step along width is too big");
    } else {
      surf::errorMessage("MeshCond::checkAspectRatio",
			 "maximum step along width is too big");
    }
  }

  if (max_hz / hy >= threshold2) {
    cout << endl << "\t Aspect ratio := " << max_hz / hy << endl;
    if (useQuadPanel_) {
      surf::errorMessage("MeshCond::checkAspectRatio",
			   "maximum step along thickness is too big");
    } else {
      surf::errorMessage("MeshCond::checkAspectRatio",
			 "maximum step along thickness is too big");
    }
  }
}

/**********************************************************************
 * transferToGlobalCoord --
 **********************************************************************/
void 
MeshCond::transferToGlobalCoord(
				void)
{
  for (int i=0; i < numNode_; i++) {
    nodeList[i].transferLocalToGlobalCoord(origin, X, Y, Z);
  }
}
