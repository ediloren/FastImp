/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: meshRectSpiral.cc,v 1.14 2002/12/31 21:32:10 zhzhu Exp $";


  ==========================================================================
*/
#include "meshRectSpiral.h"

#include <vector>
#include <iostream>
//#define OUTPUT_MESH_FOR_PLOT
//#define GEN_FASTHENRY_INP

using namespace std;
using namespace mesh;

/**********************************************************************
 * MeshRectSpiral::MeshRectSpiral --
 **********************************************************************/
MeshRectSpiral::MeshRectSpiral (
				const surf::RectSpiral* p,
				const bool useNonUniformMesh,
				const bool useQuadPanel,
				const bool aspectRatioWarning)
{
  p_ = p;
  useNonUniformMesh_ = useNonUniformMesh;
  useQuadPanel_ = useQuadPanel;
  aspectRatioWarning_ = aspectRatioWarning;

  setup();
  initMesh();
  genMesh();

#ifdef GEN_FASTHENRY_INP
  genFHOutput();
#endif

}

/**********************************************************************
 * MeshRectSpiral::setup --
 * setup those internal working objects.
 **********************************************************************/
void 
MeshRectSpiral::setup (void)
{
  setupWirePath();
  setupLineSegmentList();
  setupCornerList();

#ifdef DEBUG_MESH_RECT_SPIRAL
  wirePath_.prnInfo();
  lineSegmentList_.prnInfo();
  cornerList_.prnInfo();
#endif
}

/**********************************************************************
 * MeshRectSpiral::setupWirePath --
 **********************************************************************/
void 
MeshRectSpiral::setupWirePath (
			       void) 
{
  // fullTurns = 1; halfTurns = 0; quadTurns = 1; -> numVtx = 6;
  // fullTurns = 1; halfTurns = 1; quadTurns = 0; -> numVtx = 7;
  // fullTurns = 2; halfTurns = 1; quadTurns = 1; -> numVtx = 12;
  size_t fullTurns, halfTurns, quadTurns;
  p_->getTurns(fullTurns, halfTurns, quadTurns);

  size_t numVtx 
    = 4 * fullTurns 
    + 2 * halfTurns 
    + 1 * quadTurns + 1;

  size_t numLineSeg = numVtx - 1;

  vector<pfft::point3D> vtxArray(numVtx);
  vector<pfft::point3D> localUnitYArray(numLineSeg);
  vector<pfft::point3D> localUnitXArray(numLineSeg);
  vector<double> lengthArray(numLineSeg);

  vtxArray[0] = p_->getHead();
  localUnitYArray[0] = p_->getStartLineUnitY();
  localUnitXArray[0] = crossProd(p_->getStartLineUnitY(), p_->getPlaneNormal());
  lengthArray[0] = p_->getStartLineLength();
  
  int k = 1;
  double currentLength = p_->getStartLineLength();
  for (size_t i=1; i<numVtx; i++) { // the spiral is right rotary
    if(k == 2) {
      currentLength += p_->getSpacing() + p_->getWidth();
      k = 0;
    } k ++;
    
    if (i < numVtx - 1) {
      lengthArray[i] = currentLength;
    }
    vtxArray[i]= vtxArray[i-1] + localUnitYArray[i-1] * lengthArray[i-1];
    if (i < numVtx - 1) {
      localUnitYArray[i] = crossProd(localUnitYArray[i-1], p_->getPlaneNormal());
      localUnitXArray[i] = crossProd(localUnitYArray[i], p_->getPlaneNormal());
    }
  }
  
  wirePath_.setup(
		 p_->getWidth(),
		 p_->getThickness(),
		 vtxArray,
		 localUnitXArray,
		 localUnitYArray,
		 lengthArray);
}

/**********************************************************************
 * MeshRectSpiral::genFHOutput --
 **********************************************************************/
void 
MeshRectSpiral::genFHOutput (
			     void)
{
  std::ofstream fout("fh.inp");
  fout << "* A rectangular spiral" << std::endl;
  fout << ".units ";
  if (p_->unit() == 1e-6) {
    fout << "um";
  } else if (p_->unit() == 1e-3) {
    fout << "mm";
  } else if (p_->unit() == 1e-9) {
    fout << "nm";
  } else if (p_->unit() == 1) {
    fout <<"m";
  } else {
    std::cerr <<" unusal unit !";
    exit(1);
  }
  
  fout << std::endl;
  fout << ".default";
  fout << " nwinc=" << p_->getNumPanelAlongWidth();
  fout << " nhinc=" << p_->getNumPanelAlongThickness();
  fout << " sigma=" << p_->conductivity();
  fout << std::endl;
  fout <<".freq fmin=1e0 fmax=1e10 ndec=1" << std::endl;
  
  size_t numWirePath = wirePath_.getNumSegment();
  std::vector<pfft::point3D> nodeList;

  for (size_t i = 0; i < numWirePath; ++ i) {
    double len = wirePath_.getLength(i);
    double refLen = p_->getStartLineLength();
    double lenRatio = len/refLen;
    size_t numSegOnEachWirePath = 
      size_t(lenRatio * p_->getNumPanelAlongShortestLength());
    double pace = len / numSegOnEachWirePath;
    
    
    pfft::point3D v0 = wirePath_.getVertex(i);    
    for (size_t j = 0; j<numSegOnEachWirePath; ++j) {
      nodeList.push_back(v0 + j * pace * wirePath_.getLocalUnitY(i));
      fout << "N" << nodeList.size() << " "; 
      fout << "x=" << nodeList[nodeList.size() - 1].x() << "  ";
      fout << "y=" << nodeList[nodeList.size() - 1].y() << "  ";
      fout << "z=" << nodeList[nodeList.size() - 1].z() << std::endl;
    }
  }

  nodeList.push_back(wirePath_.getVertex(numWirePath));
  fout << "N" << nodeList.size() << " "; 
  fout << "x=" << nodeList[nodeList.size() - 1].x() << "  ";
  fout << "y=" << nodeList[nodeList.size() - 1].y() << "  ";
  fout << "z=" << nodeList[nodeList.size() - 1].z() << std::endl;

  for (size_t i = 0; i < nodeList.size()-1; ++ i) {
    fout << "E" << i+1 << " ";
    fout << "N" << i+1 << " ";
    fout << "N" << i+2 << " ";
    fout << "w=" << p_->getWidth() << " ";
    fout << "h=" << p_->getThickness() << std::endl;
  }
    
  fout << ".external" << " ";
  fout << "N1 " << "N" << nodeList.size() << std::endl;
  fout << ".end" << std::endl;
  fout.close();
}


/**********************************************************************
 * MeshRectSpiral::setupLineSegmentList --
 **********************************************************************/
void 
MeshRectSpiral::setupLineSegmentList (
				      void)
{
  lineSegmentList_.setup(p_->getWidth(), p_->getThickness(), wirePath_);
}

/**********************************************************************
 * setupCornerList --
 **********************************************************************/
void
MeshRectSpiral::setupCornerList (
				 void)
{
  cornerList_.setup(p_->getWidth(), p_->getThickness(), lineSegmentList_);
}

/**********************************************************************
 * MeshRectSpiral::initMesh --
 * Note:
 * Get 4 jobs done in this function:
 * 1. find the numPanel along the length of each segment.
 * 2. find the num Node & Panel on each line segment.
 * 3. find the num Node & Panel on each corner.
 * 4. allocate mem for vertexArray_ & panelList.
 **********************************************************************/
void MeshRectSpiral::initMesh (
			       void)
{
  size_t numLineSeg = lineSegmentList_.getNumSegment();
  size_t numCorner = cornerList_.getNumCorner();

  numNodeOnLineSegArray_.resize(numLineSeg);
  numPanelOnLineSegArray_.resize(numLineSeg);

  numNodeOnCornerArray_.resize(numLineSeg);
  numPanelOnCornerArray_.resize(numLineSeg);

#ifdef OUTPUT_MESH_FOR_PLOT
  numPanelYOnLineSegArray_.resize(numLineSeg, 1);
#else
  // setup numPanelYOnLineSegArray_; 
  numPanelYOnLineSegArray_.resize(numLineSeg);
  for (size_t i=0; i<numLineSeg; i++) {
    double length = lineSegmentList_.getLength(i);
    double ratio = length * p_->getNumPanelAlongShortestLength() / p_->getStartLineLength();
    
    /*     hard code !!
    if (ratio <= 1) {
      numPanelYOnLineSegArray_[i] = p_->getNumPanelAlongShortestLength();
    } else if (ratio < 2) {
      numPanelYOnLineSegArray_[i] = p_->getNumPanelAlongShortestLength() + 1;
    } else if (ratio < 3) {
      numPanelYOnLineSegArray_[i] = p_->getNumPanelAlongShortestLength() + 2;
    } else {
      numPanelYOnLineSegArray_[i] = p_->getNumPanelAlongShortestLength() + 3;
    }
	  
    if ((i == 0 || i == numLineSeg - 1) &&
	(numPanelYOnLineSegArray_[i] < 3)) {
      numPanelYOnLineSegArray_[i] = 3;
    } 
    if(numPanelYOnLineSegArray_[i] > 4) {
	numPanelYOnLineSegArray_[i] = 4;
    }
    */   

    numPanelYOnLineSegArray_[i] = int(floor(ratio));

#ifdef DEBUG_MESH_RECT_SPIRAL
    std::cout << " line segment =" << i << " ";
    std::cout << " length = " << length << " ";
    std::cout << " numPanelYOnLineSegArray_[i] = "<< numPanelYOnLineSegArray_[i] << " ";
    std::cout << " length / numPanelYOnLineSegArray_[i] = " << length/numPanelYOnLineSegArray_[i] << std::endl;
#endif
  }
#endif

  numNode_ = 0;
  numPanel_ = 0;
  numNonContactPanel_ = 0;

  // get num node on each conductor 
  for (size_t i=0; i<numLineSeg; i++) {
    calcNumNodePanelOnEachLineSeg(i);
  }
  for (size_t i=0; i<numCorner; i++) {
    calcNumNodePanelOnEachCorner(i);
  }
   
  nodeList.reserve(numNode_);
  panelList.reserve(numPanel_);
}

/**********************************************************************
 * MeshRectSpiral::calcNumNodePanelOnEachCorner --
 **********************************************************************/
void MeshRectSpiral::calcNumNodePanelOnEachCorner (
						   const size_t cornerIdx)
{
#ifdef OUTPUT_MESH_FOR_PLOT
  size_t numPanelX = 1;
  size_t numPanelZ = 1;
#else
  size_t numPanelX = p_->getNumPanelAlongWidth();
  size_t numPanelZ = p_->getNumPanelAlongThickness();
#endif

  size_t numPanelY = numPanelX;

  size_t numNodeY = numPanelY + 1; // if it's independant; but 
                                   // actually, the num of nodes 
                                   // along the first contour should
                                   // not be counted.

  size_t numNodeOnSides 
    = (2 * numPanelX + numPanelZ + 1) * (numPanelY - 1) 
    + 2 * (numPanelX + numPanelZ);

  size_t numNodeOnEnds = (numPanelX - 1) * (numPanelZ - 1);
    
  size_t numPanelOnSides = (2 * numPanelX + numPanelZ) * numPanelY;
  size_t numPanelOnEnds = (numPanelX * numPanelZ);

  numNodeOnCornerArray_[cornerIdx] = numNodeOnSides + numNodeOnEnds;
  numPanelOnCornerArray_[cornerIdx] = numPanelOnSides + numPanelOnEnds;
  
  numNode_ += numNodeOnCornerArray_[cornerIdx];
  numPanel_ += numPanelOnCornerArray_[cornerIdx];
}

/**********************************************************************
 * MeshRectSpiral::calcNumNodePanelOnEachLineSeg --
 **********************************************************************/
void 
MeshRectSpiral::calcNumNodePanelOnEachLineSeg (
					       const size_t lnSegIdx)
{
#ifdef OUTPUT_MESH_FOR_PLOT
  size_t numPanelX = 1;
  size_t numPanelZ = 1;
#else
  size_t numPanelX = p_->getNumPanelAlongWidth();
  size_t numPanelZ = p_->getNumPanelAlongThickness();
#endif

  size_t numPanelY = numPanelYOnLineSegArray_[lnSegIdx];

  // this is the real number along Y direction (the share cases considered)
  size_t numNodeY = (lnSegIdx == 0) ? numPanelY + 1 : numPanelY;

  size_t numNodeOnSides = 2 * (numPanelX + numPanelZ) * numNodeY;
  size_t numNodeOnEnds = 0;
  if (lnSegIdx == 0) {
    numNodeOnEnds += (numPanelX - 1) * (numPanelZ - 1);
  } 
  if (lnSegIdx == lineSegmentList_.getNumSegment()- 1) {
    numNodeOnEnds += (numPanelX - 1) * (numPanelZ - 1);
  }

  size_t numPanelOnSides = 2 * (numPanelX + numPanelZ) * numPanelY;
  size_t numPanelOnEnds = 0;
  if (lnSegIdx == 0) {
    numPanelOnEnds += (numPanelX * numPanelZ);
  } 
  if (lnSegIdx == lineSegmentList_.getNumSegment()- 1) {
    numPanelOnEnds += (numPanelX * numPanelZ);
  }
  
  numNodeOnLineSegArray_[lnSegIdx] = numNodeOnSides + numNodeOnEnds;
  numPanelOnLineSegArray_[lnSegIdx] = numPanelOnSides + numPanelOnEnds;

  numNode_ += numNodeOnLineSegArray_[lnSegIdx];
  numPanel_ += numPanelOnLineSegArray_[lnSegIdx];
}

/**********************************************************************
 * getMap4tartContour --
 **********************************************************************/
void 
MeshRectSpiral::getMap4StartContour (void)
{
  map4ContourList[0].bottom = map4LastContour.bottom;
  map4ContourList[0].back = map4LastContour.back;
  map4ContourList[0].top = map4LastContour.top;
  map4ContourList[0].front = map4LastContour.front;
}

/**********************************************************************
 * MeshRectSpiral::genMesh --
 **********************************************************************/
void 
MeshRectSpiral::genMesh (void)
{
  numNodeSofar_ = 0;
  size_t numLineSegment = lineSegmentList_.getNumSegment();

  for (size_t i=0; i<numLineSegment; i++) {
    size_t lineSegIdx = i;
    meshLineSegment(lineSegIdx);
    if (i != numLineSegment - 1) {
      size_t cornerIdx = i;
      meshCorner(cornerIdx);
    }
  }
}

/**********************************************************************
 * MeshRectSpiral::meshLineSegment --
 **********************************************************************/
void 
MeshRectSpiral::meshLineSegment (
				 const size_t lnSegIdx)
{
#ifdef OUTPUT_MESH_FOR_PLOT
  size_t numPanelX = 1;
  size_t numPanelZ = 1;
#else
  size_t numPanelX = p_->getNumPanelAlongWidth();
  size_t numPanelZ = p_->getNumPanelAlongThickness();
#endif

  size_t numPanelY = numPanelYOnLineSegArray_[lnSegIdx];

  size_t numNodeX = numPanelX + 1;
  size_t numNodeY = numPanelY + 1;
  size_t numNodeZ = numPanelZ + 1;

  map4ContourList.resize(numNodeY);
  for (size_t i=0; i<numNodeY; i++) {
    map4ContourList[i].bottom.resize(numPanelX);
    map4ContourList[i].back.resize(numPanelZ);
    map4ContourList[i].top.resize(numPanelX);
    map4ContourList[i].front.resize(numPanelZ);
  }
  
  if (lnSegIdx == 0) {
    map4TwoEnds.left.resize(numNodeX);
    for (size_t i=0; i<numNodeX; i++) {
      map4TwoEnds.left[i].resize(numNodeZ);
    }
  } 

  if (lnSegIdx == lineSegmentList_.getNumSegment() - 1) {
    map4TwoEnds.right.resize(numNodeX);
    for (size_t i=0; i<numNodeX; i++) {
      map4TwoEnds.right[i].resize(numNodeZ);
    }
  }
  
  if (lnSegIdx != 0) {
    getMap4StartContour(); // put map4LastContour into map4ContourList 
  }

  findNodeOnLineSegment(lnSegIdx); // now, map4LastContour has been updated !! 
  findPanelOnLineSegment(lnSegIdx);
}

/**********************************************************************
 * MeshRectSpiral::meshCorner --
 **********************************************************************/
void MeshRectSpiral::meshCorner (
				     const size_t cornerIdx)
{
#ifdef OUTPUT_MESH_FOR_PLOT
  size_t numPanelX = 1;
  size_t numPanelZ = 1;
  size_t numPanelY = numPanelX;
#else
  size_t numPanelX = p_->getNumPanelAlongWidth();
  size_t numPanelZ = p_->getNumPanelAlongThickness();
  size_t numPanelY = p_->getNumPanelAlongWidth();
#endif

  size_t numNodeX = numPanelX + 1;
  size_t numNodeY = numPanelY + 1;
  size_t numNodeZ = numPanelZ + 1;

  map4ContourList.resize(numNodeY);
  for (size_t i=0; i<numNodeY; i++) {
    map4ContourList[i].bottom.resize(numPanelX);
    map4ContourList[i].back.resize(numPanelZ);
    map4ContourList[i].top.resize(numPanelX);
    map4ContourList[i].front.resize(numPanelZ);
  }
  
  map4TwoEnds.right.resize(numNodeX);
  for (size_t i=0; i<map4TwoEnds.right.size(); i++) {
    map4TwoEnds.right[i].resize(numNodeZ);
  }
    
  getMap4StartContour(); // put map4LastContour into map4ContourList 

  findNodeOnCorner(cornerIdx); // now, map4LastContour has been updated !!
  findPanelOnCorner(cornerIdx);
}

/**********************************************************************
 * MeshRectSpiral::findNodeOnLineSegment --
 * note:
 * 1. get the coordinate of each node.
 * 2. get the map of local contour node index to global node index
 *    It'll be used in finding panel.
 * 3. put the global node index in map4LastContour, which will help
 *    to mesh the following corner.
 **********************************************************************/
void MeshRectSpiral::findNodeOnLineSegment (
					    const size_t lnSegIdx)
{
  double width = lineSegmentList_.getWidth();
  double thickness = lineSegmentList_.getThickness();
  double length = lineSegmentList_.getLength(lnSegIdx);

  const double x1 = -.5 * width;
  const double x2 = +.5 * width;
  const double y1 = 0.;
  const double y2 = length;
  const double z1 = -.5 * thickness;
  const double z2 = +.5 * thickness;

#ifdef OUTPUT_MESH_FOR_PLOT
  size_t numPanelX = 1;
  size_t numPanelZ = 1;
#else
  size_t numPanelX = p_->getNumPanelAlongWidth();
  size_t numPanelZ = p_->getNumPanelAlongThickness();
#endif

  size_t numPanelY = numPanelYOnLineSegArray_[lnSegIdx];

  size_t numNodeX = numPanelX + 1;
  size_t numNodeY = numPanelY + 1;
  size_t numNodeZ = numPanelZ + 1;

  // use for double check !!
  size_t numNodeOnSides = 0;
  size_t numNodeOnEnds = 0;
  
  pfft::point3D origin = lineSegmentList_.getStartVertex(lnSegIdx);
  pfft::point3D unitX = lineSegmentList_.getLocalUnitX(lnSegIdx);
  pfft::point3D unitY = lineSegmentList_.getLocalUnitY(lnSegIdx);
  pfft::point3D unitZ = lineSegmentList_.getLocalUnitZ(lnSegIdx);

  if (! useNonUniformMesh_) {
    
    double hx = width / numPanelX;
    double hz = thickness / numPanelZ;
    double hy = length / numPanelY;
    
    double x, y, z;
    for (size_t iy = (lnSegIdx == 0) ? 0 : 1; iy < numNodeY; iy ++) {
      y = y1 + iy * hy;
      
      // bottom
      z = z1;
      for (size_t ix=0; ix<numPanelX; ix++) {
	x = x1 + ix * hx;
	map4ContourList[iy].bottom[ix] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // back 
      x = x2;
      for (size_t iz=0; iz<numPanelZ; iz++) {
	z = z1 + iz * hz;
	map4ContourList[iy].back[iz] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // top
      z = z2;
      for (size_t ix=0; ix<numPanelX; ix++) {
	x = x2 - ix * hx;
	map4ContourList[iy].top[ix] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // front
      x = x1;
      for (size_t iz=0; iz<numPanelZ; iz++) {
	z = z2 - iz * hz;
	map4ContourList[iy].front[iz] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
    }
    
    if (lnSegIdx == lineSegmentList_.getNumSegment() - 1) {
      // setup right end !!
      y = y2;
      for (size_t iz=0; iz<numNodeZ; iz++) {
	z = z1 + iz * hz;
	for (size_t ix=0; ix<numNodeX; ix++) {
	  x = x1 + ix * hx;
	  if (iz == 0) {
	    if (ix == numNodeX - 1) {
	      map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].back[0];
	    } else {
	      map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].bottom[ix];
	    }
	  } else if (iz == numNodeZ - 1) { // iz != 0
	    if (ix == 0) {
	      map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].front[0];
	    } else {
	      map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].top[numPanelX - ix];
	    }
	  } else if (ix == 0) { // iz !=0 && iz != numNodeZ - 1
	    map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].front[numPanelZ - iz];
	  } else if (ix == numNodeX - 1) { // iz != 0 && iz != numNodeZ - 1
	    map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].back[iz];
	  } else {
	    map4TwoEnds.right[ix][iz] = numNodeSofar_ ++;
	    nodeList.push_back(
				   local2GlobalCoordSys(origin,
							unitX,
							unitY,
							unitZ,
							x, y, z));
	    numNodeOnEnds ++;
	  }
	}
      }
    } 

    if (lnSegIdx == 0) {
      // setup left end !!
      y = y1;
      for (size_t iz=0; iz<numNodeZ; iz++) {
	z = z1 + iz * hz;
	for (size_t ix=0; ix<numNodeX; ix++) {
	  x = x1 + ix * hx;
	  if (iz == 0) {
	    if (ix == numNodeX - 1) {
	      map4TwoEnds.left[ix][iz] = map4ContourList[0].back[0];
	    } else {
	      map4TwoEnds.left[ix][iz] = map4ContourList[0].bottom[ix];
	  }
	  } else if (iz == numNodeZ - 1) { // iz != 0
	    if (ix == 0) {
	      map4TwoEnds.left[ix][iz] = map4ContourList[0].front[0];
	    } else {
	      map4TwoEnds.left[ix][iz] = map4ContourList[0].top[numPanelX - ix];
	    }
	  } else if (ix == 0) { // iz !=0 && iz != numNodeZ - 1
	    map4TwoEnds.left[ix][iz] = map4ContourList[0].front[numPanelZ - iz];
	  } else if (ix == numNodeX - 1) { // iz != 0 && iz != numNodeZ - 1
	    map4TwoEnds.left[ix][iz] = map4ContourList[0].back[iz];
	  } else {
	    map4TwoEnds.left[ix][iz] = numNodeSofar_ ++;
	    nodeList.push_back(
				   local2GlobalCoordSys(origin,
							unitX,
							unitY,
							unitZ,
							x, y, z));
	    numNodeOnEnds ++;
	  }
	}
      }
    } 
  } else {
    
    findNonUniformStep(width, numPanelX, hx);
    findNonUniformStep(thickness, numPanelZ, hz);
    double hy = length / numPanelY;
    
    if (aspectRatioWarning_) {
      checkAspectRatio(hx, hz, hy);
    }

    // use for double check !!    
    double x, y, z;
    for (size_t iy = (lnSegIdx == 0) ? 0 : 1; iy < numNodeY; iy ++) {
      y = y1 + iy * hy;

      // bottom
      z = z1;
      for (size_t ix=0; ix<numPanelX; ix++) {
	x = (ix == 0) 
	  ? x1
	  : x + hx[ix-1];
	map4ContourList[iy].bottom[ix] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // back 
      x = x2;
      for (size_t iz=0; iz<numPanelZ; iz++) {
	z = (iz == 0)
	  ? z1
	  : z + hz[iz-1];
	map4ContourList[iy].back[iz] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // top
      z = z2;
      for (size_t ix=0; ix<numPanelX; ix++) {
	x = (ix == 0)
	  ? x2
	  : x - hx[ix-1];
	map4ContourList[iy].top[ix] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // front
      x = x1;
      for (size_t iz=0; iz<numPanelZ; iz++) {
	z = (iz == 0)
	  ? z2
	  : z - hz[iz-1];
	map4ContourList[iy].front[iz] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
    }
    
    if (lnSegIdx == lineSegmentList_.getNumSegment() - 1) {
      // setup right end !!
      y = y2;
      for (size_t iz=0; iz<numNodeZ; iz++) {
	z = (iz == 0)
	  ? z1 
	  : z + hz[iz-1];
	for (size_t ix=0; ix<numNodeX; ix++) {
	  x = (ix == 0) 
	    ? x1
	    : x + hx[ix-1];
	  if (iz == 0) {
	  if (ix == numNodeX - 1) {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].back[0];
	  } else {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].bottom[ix];
	  }
	  } else if (iz == numNodeZ - 1) { // iz != 0
	    if (ix == 0) {
	      map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].front[0];
	    } else {
	      map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].top[numPanelX - ix];
	    }
	  } else if (ix == 0) { // iz !=0 && iz != numNodeZ - 1
	    map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].front[numPanelZ - iz];
	  } else if (ix == numNodeX - 1) { // iz != 0 && iz != numNodeZ - 1
	    map4TwoEnds.right[ix][iz] = map4ContourList[numNodeY-1].back[iz];
	  } else {
	    map4TwoEnds.right[ix][iz] = numNodeSofar_ ++;
	    nodeList.push_back(
				   local2GlobalCoordSys(origin,
							unitX,
							unitY,
							unitZ,
							x, y, z));
	    numNodeOnEnds ++;
	  }
	}
      }
    } 
    
    if (lnSegIdx == 0) {
      // setup left end !!
      y = y1;
      for (size_t iz=0; iz<numNodeZ; iz++) {
	z = (iz == 0)
	  ? z1
	  : z + hz[iz-1];
	for (size_t ix=0; ix<numNodeX; ix++) {
	  x = (ix == 0) 
	    ? x1
	    : x + hx[ix-1];
	  if (iz == 0) {
	    if (ix == numNodeX - 1) {
	      map4TwoEnds.left[ix][iz] = map4ContourList[0].back[0];
	    } else {
	      map4TwoEnds.left[ix][iz] = map4ContourList[0].bottom[ix];
	    }
	  } else if (iz == numNodeZ - 1) { // iz != 0
	    if (ix == 0) {
	      map4TwoEnds.left[ix][iz] = map4ContourList[0].front[0];
	    } else {
	      map4TwoEnds.left[ix][iz] = map4ContourList[0].top[numPanelX - ix];
	    }
	  } else if (ix == 0) { // iz !=0 && iz != numNodeZ - 1
	    map4TwoEnds.left[ix][iz] = map4ContourList[0].front[numPanelZ - iz];
	  } else if (ix == numNodeX - 1) { // iz != 0 && iz != numNodeZ - 1
	    map4TwoEnds.left[ix][iz] = map4ContourList[0].back[iz];
	  } else {
	    map4TwoEnds.left[ix][iz] = numNodeSofar_ ++;
	    nodeList.push_back(
				   local2GlobalCoordSys(origin,
							unitX,
							unitY,
							unitZ,
							x, y, z));
	    numNodeOnEnds ++;
	  }
	}
      }
    }
  }
  
  // for debug use:
  if (numNodeOnEnds + numNodeOnSides != numNodeOnLineSegArray_[lnSegIdx]) {
    cerr << " Error! rectSpiral.cc: findNodeOnLineSegment() " << endl;
    cerr << " when meshing lineSeg " << lnSegIdx << endl;    
    cerr << " # nodes doesn't match ! (" << numNodeOnEnds + numNodeOnSides;
    cerr << " , " << numNodeOnLineSegArray_[lnSegIdx] << endl;
    cerr << " must be a bug!  ... exit " << endl;
    exit(1);
  }
  
  // fill in map4LastContour !
  if (lnSegIdx != lineSegmentList_.getNumSegment() - 1) {
    map4LastContour.bottom = map4ContourList[numNodeY - 1].bottom;
    map4LastContour.back = map4ContourList[numNodeY - 1].back;
    map4LastContour.top = map4ContourList[numNodeY - 1].top;
    map4LastContour.front = map4ContourList[numNodeY - 1].front;
  }
}

/**********************************************************************
 * MeshRectSpiral::findPanelOnLineSegment --
 **********************************************************************/
void MeshRectSpiral::findPanelOnLineSegment (
					     const size_t lnSegIdx)
{
  PanelType panelType;

  vector<int> nodeIndex(4);

  size_t numLineSegment = lineSegmentList_.getNumSegment();

#ifdef OUTPUT_MESH_FOR_PLOT
  size_t numPanelX = 1;
  size_t numPanelZ = 1;
#else
  size_t numPanelX = p_->getNumPanelAlongWidth();
  size_t numPanelZ = p_->getNumPanelAlongThickness();
#endif

  size_t numPanelY = numPanelYOnLineSegArray_[lnSegIdx];

  size_t numPanelOnSides = 0;
  size_t numPanelOnEnds = 0;

  for (size_t iy=0; iy<numPanelY; iy++) {
    // bottom !
    for (size_t ix=0; ix<numPanelX; ix++) {
      nodeIndex[0] = map4ContourList[iy].bottom[ix];

      nodeIndex[1] = (ix == numPanelX - 1) 
	   ? map4ContourList[iy].back[0]
	   : map4ContourList[iy].bottom[ix+1]; 

      nodeIndex[2] = (ix == numPanelX - 1)
	? map4ContourList[iy+1].back[0]
	: map4ContourList[iy+1].bottom[ix+1];

      nodeIndex[3] = map4ContourList[iy+1].bottom[ix];

      if (((lnSegIdx == 0) && (iy < 2))||
	  ((lnSegIdx == numLineSegment -1) && (iy > numPanelY - 3))){
	panelType = BUFFER;
      } else {
	panelType = NON_CONTACT;
      }

      panelList.push_back(Panel(nodeIndex, panelType));
      ++ numPanelOnSides;
      ++ numNonContactPanel_;
    }

    // back !
    for (size_t iz=0; iz<numPanelZ; iz++) {
      nodeIndex[0] = map4ContourList[iy].back[iz];

      nodeIndex[1] = (iz == numPanelZ - 1) 
	? map4ContourList[iy].top[0]
	: map4ContourList[iy].back[iz+1];

      nodeIndex[2] = (iz == numPanelZ - 1) 
	? map4ContourList[iy+1].top[0]
	: map4ContourList[iy+1].back[iz+1];
      
      nodeIndex[3] = map4ContourList[iy+1].back[iz];

      if (((lnSegIdx == 0) && (iy < 2))||
	  ((lnSegIdx == numLineSegment -1) && (iy > numPanelY - 3))){
	panelType = BUFFER;
      } else {
	panelType = NON_CONTACT;
      }
      
      panelList.push_back(Panel(nodeIndex, panelType));
      ++ numPanelOnSides;
      ++ numNonContactPanel_;
    } 
    
    // top !
    for (size_t ix=0; ix<numPanelX; ix++) {
      nodeIndex[0] = map4ContourList[iy].top[ix];

      nodeIndex[1] = (ix == numPanelX - 1) 
	? map4ContourList[iy].front[0]
	: map4ContourList[iy].top[ix+1];
      
      nodeIndex[2] = (ix == numPanelX - 1) 
	? map4ContourList[iy+1].front[0]
	: map4ContourList[iy+1].top[ix+1];
      
      nodeIndex[3] = map4ContourList[iy+1].top[ix];
      
      if (((lnSegIdx == 0) && (iy < 2))||
	  ((lnSegIdx == numLineSegment -1) && (iy > numPanelY - 3))){
	panelType = BUFFER;
      } else {
	panelType = NON_CONTACT;
      }
	    
      panelList.push_back(Panel(nodeIndex, panelType));
      ++ numPanelOnSides;
      ++ numNonContactPanel_;
    } 
  
    // front !
    for (size_t iz=0; iz<numPanelZ; iz++) {
      nodeIndex[0] = map4ContourList[iy].front[iz];

      nodeIndex[1] = (iz == numPanelZ - 1) 
	? map4ContourList[iy].bottom[0]
	: map4ContourList[iy].front[iz+1];

      nodeIndex[2] = (iz == numPanelZ - 1) 
	? map4ContourList[iy+1].bottom[0]
	: map4ContourList[iy+1].front[iz+1];

      nodeIndex[3] = map4ContourList[iy+1].front[iz];
      
      if (((lnSegIdx == 0) && (iy < 2))||
	  ((lnSegIdx == numLineSegment -1) && (iy > numPanelY - 3))){
	panelType = BUFFER;
      } else {
	panelType = NON_CONTACT;
      }

      panelList.push_back(Panel(nodeIndex, panelType));
      ++ numPanelOnSides;
      ++ numNonContactPanel_;
    }
  }

  // panels on ends 
  if (lnSegIdx == 0) {
    for (size_t iz=0; iz<numPanelZ; iz++) {
      for (size_t ix=0; ix<numPanelX; ix ++) {
	nodeIndex[0] = map4TwoEnds.left[ix][iz+1];
	nodeIndex[1] = map4TwoEnds.left[ix+1][iz+1];
	nodeIndex[3] = map4TwoEnds.left[ix][iz];
	nodeIndex[2] = map4TwoEnds.left[ix+1][iz];
	
	PanelType type = LEFT_CONTACT;
	panelList.push_back(Panel(nodeIndex,type));
	++ numPanelOnEnds;
      }
    }
  } 

  if (lnSegIdx == lineSegmentList_.getNumSegment() - 1) {
    for (size_t iz=0; iz<numPanelZ; iz++) {
      for (size_t ix=0; ix<numPanelX; ix++) {
	nodeIndex[0] = map4TwoEnds.right[ix][iz];
	nodeIndex[1] = map4TwoEnds.right[ix+1][iz];
	nodeIndex[2] = map4TwoEnds.right[ix+1][iz+1];
	nodeIndex[3] = map4TwoEnds.right[ix][iz+1];
	
	PanelType type = RIGHT_CONTACT;
	panelList.push_back(Panel(nodeIndex,type));
	++ numPanelOnEnds;
      }
    }
  } 

  if (numPanelOnSides + numPanelOnEnds != numPanelOnLineSegArray_[lnSegIdx]) {
    cerr << " Error! rectSpiral.cc: findPanelOnLineSegment() " << endl;
    cerr << " when meshing lineSeg " << lnSegIdx << endl;    
    cerr << " # panels doesn't match ! (" << numPanelOnEnds + numPanelOnSides;
    cerr << " , " << numPanelOnLineSegArray_[lnSegIdx] << endl;
    cerr << " must be a bug!  ... exit " << endl;
    exit(1);
  }
}

/**********************************************************************
 * findNodeOnCorner --
 **********************************************************************/
void  
MeshRectSpiral::findNodeOnCorner (
				  const size_t cornerIdx)
{
  double width = cornerList_.getWidth();
  double thickness = cornerList_.getThickness();
  double length = cornerList_.getLength();

  const double x1 = -.5 * width;
  const double x2 = +.5 * width;
  const double y1 = 0.;
  const double y2 = length;
  const double z1 = -.5 * thickness;
  const double z2 = +.5 * thickness;

#ifdef OUTPUT_MESH_FOR_PLOT
  size_t numPanelX = 1;
  size_t numPanelZ = 1;
#else
  size_t numPanelX = p_->getNumPanelAlongWidth();
  size_t numPanelZ = p_->getNumPanelAlongThickness();
#endif
  size_t numPanelY = numPanelX;

  size_t numNodeX = numPanelX + 1;
  size_t numNodeY = numPanelY + 1;
  size_t numNodeZ = numPanelZ + 1;

  pfft::point3D origin = cornerList_.getLeftSegCentroid(cornerIdx);
  pfft::point3D unitX = cornerList_.getLocalUnitX(cornerIdx);
  pfft::point3D unitY = cornerList_.getLocalUnitY(cornerIdx);
  pfft::point3D unitZ = cornerList_.getLocalUnitZ(cornerIdx);

  // use for double check !!
  size_t numNodeOnSides = 0;
  size_t numNodeOnEnds = 0;

  if (! useNonUniformMesh_) {
    double hx = width / numPanelX;
    double hz = thickness / numPanelZ;
    double hy = length / numPanelY;
    
    
    double x, y, z;
    for (size_t iy = 1; iy < numNodeY; iy ++) {
      y = y1 + iy * hy;
      
      // bottom
      z = z1;
      for (size_t ix=0; ix<numPanelX; ix++) {
	x = x1 + ix * hx;
	map4ContourList[iy].bottom[ix] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // back 
      x = x2;
      for (size_t iz=0; 
	   iz<((iy == numNodeY - 1) ? numPanelZ : 1);
	   iz ++) {
	z = z1 + iz * hz;
	map4ContourList[iy].back[iz] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // top
      z = z2;
      for (size_t ix=0; ix<numPanelX; ix++) {
	x = x2 - ix * hx;
	map4ContourList[iy].top[ix] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // front
      x = x1;
      for (size_t iz=0; iz<numPanelZ; iz++) {
	z = z2 - iz * hz;
	map4ContourList[iy].front[iz] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
    }
    
    // nodes on right end
    y = y2;
    for (size_t iz=0; iz<numNodeZ; iz++) {
      z = z1 + iz * hz;
      for (size_t ix=0; ix<numNodeX; ix++) {
	x = x1 + ix * hx;
	if (iz == 0) {
	  if (ix == numNodeX - 1) {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].back[0];
	  } else {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].bottom[ix];
	  }
	} else if (iz == numNodeZ - 1) { // iz != 0
	  if (ix == 0) {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].front[0];
	  } else {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].top[numPanelX - ix];
	  }
	} else if (ix == 0) { // iz !=0 && iz != numNodeZ - 1
	  map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].front[numPanelZ - iz];
	} else if (ix == numNodeX - 1) { // iz != 0 && iz != numNodeZ - 1
	  map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].back[iz];
	} else {
	  map4TwoEnds.right[ix][iz] = numNodeSofar_ ++;
	  nodeList.push_back(
				 local2GlobalCoordSys(
						      origin,
						      unitX,
						      unitY,
						      unitZ,
						      x, y, z));
	  
	  numNodeOnEnds ++;
	}
      }
    }
  } else { // non-uniform mesh !!

    TNT::Vector<double> hy;
    findNonUniformStep(width, numPanelX, hx);
    findNonUniformStep(thickness, numPanelZ, hz);
    findNonUniformStep(length, numPanelY, hy);
    
    double x, y, z;
    y = y1;
    for (size_t iy = 1; iy < numNodeY; iy ++) {
      y = y + hy[iy-1];
      
      // bottom
      z = z1;
      for (size_t ix=0; ix<numPanelX; ix++) {
	x = (ix == 0)
	  ? x1
	  : x + hx[ix-1];
	map4ContourList[iy].bottom[ix] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // back 
      x = x2;
      for (size_t iz=0; 
	   iz<((iy == numNodeY - 1) ? numPanelZ : 1);
	   iz ++) {
	z = (iz == 0)
	  ? z1
	  : z + hz[iz-1];
	map4ContourList[iy].back[iz] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // top
      z = z2;
      for (size_t ix=0; ix<numPanelX; ix++) {
	x = (ix == 0)
	  ? x2
	  : x - hx[ix-1];
	map4ContourList[iy].top[ix] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
      
      // front
      x = x1;
      for (size_t iz=0; iz<numPanelZ; iz++) {
	z = (iz == 0)
	  ? z2
	  : z - hz[iz-1];
	map4ContourList[iy].front[iz] = numNodeSofar_ ++;
	nodeList.push_back(
			       local2GlobalCoordSys(
						    origin,
						    unitX,
						    unitY,
						    unitZ,
						    x, y, z));
	numNodeOnSides ++;
      }
    }
    
    // nodes on right end
    y = y2;
    for (size_t iz=0; iz<numNodeZ; iz++) {
      z = (iz == 0) 
	? z1
	: z + hz[iz-1];
      for (size_t ix=0; ix<numNodeX; ix++) {
	x = (ix == 0)
	  ? x1
	  : x + hx[ix-1];
	if (iz == 0) {
	  if (ix == numNodeX - 1) {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].back[0];
	  } else {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].bottom[ix];
	  }
	} else if (iz == numNodeZ - 1) { // iz != 0
	  if (ix == 0) {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].front[0];
	  } else {
	    map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].top[numPanelX - ix];
	  }
	} else if (ix == 0) { // iz !=0 && iz != numNodeZ - 1
	  map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].front[numPanelZ - iz];
	} else if (ix == numNodeX - 1) { // iz != 0 && iz != numNodeZ - 1
	  map4TwoEnds.right[ix][iz] = map4ContourList[numPanelY].back[iz];
	} else {
	  map4TwoEnds.right[ix][iz] = numNodeSofar_ ++;
	  nodeList.push_back(
				 local2GlobalCoordSys(
						      origin,
						      unitX,
						      unitY,
						      unitZ,
						      x, y, z));
	  
	  numNodeOnEnds ++;
	}
      }
    }
  }

  // for debug use:
  if (numNodeOnEnds + numNodeOnSides != numNodeOnCornerArray_[cornerIdx]) {
    cerr << " Error! rectSpiral.cc: findNodeOnCorner() " << endl;
    cerr << " when meshing corner " << cornerIdx << endl;
    cerr << " # nodes doesn't match ! (" << numNodeOnEnds + numNodeOnSides;
    cerr << " , " << numNodeOnCornerArray_[cornerIdx] << endl;
    cerr << " must be a bug!  ... exit " << endl;
    exit(1);
  }
  
  // fill in map4LastContour !
  map4LastContour.bottom.resize(numPanelY);
  map4LastContour.back.resize(numPanelZ);
  map4LastContour.top.resize(numPanelY);
  map4LastContour.front.resize(numPanelZ);
  
  // bottom 
  for (size_t iy = 0; iy < numPanelY; iy ++) {
    map4LastContour.bottom[iy] = map4ContourList[numPanelY - iy].back[0];
  }
  // back
  for (size_t iz = 0; iz < numPanelZ; iz ++) {
    map4LastContour.back[iz] = map4ContourList[0].back[iz];
  }
  // top
  for (size_t iy = 0; iy < numPanelY; iy ++) {
    map4LastContour.top[iy] = map4ContourList[iy].top[0];
  }
  // front
  for (size_t iz = 0; iz < numPanelZ; iz ++) {
    map4LastContour.front[iz] = (iz == 0) 
      ? map4ContourList[numPanelY].top[0]
      : map4ContourList[numPanelY].back[numPanelZ - iz];
  }
}

/**********************************************************************
 * MeshRectSpiral::finPanelOnCorner --
 **********************************************************************/
void 
MeshRectSpiral::findPanelOnCorner (
				   const size_t cornerIdx)
{
  PanelType panelType;

  vector<int> nodeIndex(4);

  size_t numCorner = cornerList_.getNumCorner();

#ifdef OUTPUT_MESH_FOR_PLOT
  size_t numPanelX = 1;
  size_t numPanelZ = 1;
#else
  size_t numPanelX = p_->getNumPanelAlongWidth();
  size_t numPanelZ = p_->getNumPanelAlongThickness();
#endif
  size_t numPanelY = numPanelX;


  size_t numPanelOnSides = 0;
  size_t numPanelOnEnds = 0;

  for (size_t iy=0; iy<numPanelY; iy++) {
    // bottom !
    for (size_t ix=0; ix<numPanelX; ix++) {
      nodeIndex[0] = map4ContourList[iy].bottom[ix];

      nodeIndex[1] = (ix == numPanelX - 1) 
	   ? map4ContourList[iy].back[0]
	   : map4ContourList[iy].bottom[ix+1]; 

      nodeIndex[2] = (ix == numPanelX - 1)
	? map4ContourList[iy+1].back[0]
	: map4ContourList[iy+1].bottom[ix+1];
	
      nodeIndex[3] = map4ContourList[iy+1].bottom[ix];
      
      panelType = NON_CONTACT;
      
      panelList.push_back(Panel(nodeIndex, panelType));
      ++ numPanelOnSides;
      ++ numNonContactPanel_;
    }
    
    // top !
    for (size_t ix=0; ix<numPanelX; ix++) {
      nodeIndex[0] = map4ContourList[iy].top[ix];

      nodeIndex[1] = (ix == numPanelX - 1) 
	? map4ContourList[iy].front[0]
	: map4ContourList[iy].top[ix+1];
      
      nodeIndex[2] = (ix == numPanelX - 1) 
	? map4ContourList[iy+1].front[0]
	: map4ContourList[iy+1].top[ix+1];
      
      nodeIndex[3] = map4ContourList[iy+1].top[ix];
      
      panelType = NON_CONTACT;
      
      panelList.push_back(Panel(nodeIndex, panelType));
      ++ numPanelOnSides;
      ++ numNonContactPanel_;
    } 

    // front !
    for (size_t iz=0; iz<numPanelZ; iz++) {
      nodeIndex[0] = map4ContourList[iy].front[iz];

      nodeIndex[1] = (iz == numPanelZ - 1) 
	? map4ContourList[iy].bottom[0]
	: map4ContourList[iy].front[iz+1];

      nodeIndex[2] = (iz == numPanelZ - 1) 
	? map4ContourList[iy+1].bottom[0]
	: map4ContourList[iy+1].front[iz+1];

      nodeIndex[3] = map4ContourList[iy+1].front[iz];
      
      panelType = NON_CONTACT;
      
      panelList.push_back(Panel(nodeIndex, panelType));
      ++ numPanelOnSides;
      ++ numNonContactPanel_;
    }
  }

  // panels on ends 
  for (size_t iz=0; iz<numPanelZ; iz++) {
    for (size_t ix=0; ix<numPanelX; ix++) {
      nodeIndex[0] = map4TwoEnds.right[ix][iz];
      nodeIndex[1] = map4TwoEnds.right[ix+1][iz];
      nodeIndex[2] = map4TwoEnds.right[ix+1][iz+1];
      nodeIndex[3] = map4TwoEnds.right[ix][iz+1];
      
      PanelType type = NON_CONTACT;
      panelList.push_back(Panel(nodeIndex,type));
      ++ numPanelOnEnds;
      ++ numNonContactPanel_;
    }
  }
  
  if (numPanelOnSides + numPanelOnEnds != numPanelOnCornerArray_[cornerIdx]) {
    cerr << " Error! rectSpiral.cc: findPanelOnCorner() " << endl;
    cerr << " when meshing corner " << cornerIdx << endl;
    cerr << " # panels doesn't match ! (" << numPanelOnEnds + numPanelOnSides;
    cerr << " , " << numPanelOnCornerArray_[cornerIdx] << endl;
    cerr << " must be a bug!  ... exit " << endl;
    exit(1);
  }
}

/**********************************************************************
 * MeshRectSpiral::local2GlobalCoordSys --
 **********************************************************************/
pfft::point3D MeshRectSpiral::local2GlobalCoordSys (
						    const pfft::point3D& origin,
						    const pfft::point3D& X,
						    const pfft::point3D& Y,
						    const pfft::point3D& Z,
						    double x, 
						    double y, 
						    double z)
{
  pfft::point3D v = x * X + y * Y + z * Z + origin;
  return v;
}
