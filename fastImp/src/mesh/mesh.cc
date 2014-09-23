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

  const static char cvsid[] = "$Id: mesh.cc,v 1.12 2003/03/09 16:49:33 zhzhu Exp $";

  ==========================================================================
*/

#include "mesh.h"
#include "meshCond.h"
#include "meshWire.h"
#include "meshRing.h"
#include "meshSpiral.h"
#include "meshGround.h"
#include "meshRectSpiral.h"
#include "meshOtherShape.h"
#include "vertex.h"
#include <fstream> // for std::ofstream
#include <iomanip> // for std::setprecision()
#include <ctype.h> // for isspace

#include <sys/time.h>
#include <assert.h>
#include <string.h>

using namespace std;
using namespace mesh;
using namespace surf;
using namespace pfft;

/**********************************************************************
 * generateMesh --
 **********************************************************************/
void
Mesh::generateMesh (
		    const std::vector<CondInfo*>& condInfoPtrList,
		    const bool useNonUniformMesh,
		    const bool useQuadPanel,
		    const bool aspectRatioWarning)
{
  numCond_ = condInfoPtrList.size();
  meshCondPtrList.resize(numCond_, 0);
  unit_ = condInfoPtrList[0]->unit();
  for (size_t i = 0; i < numCond_; i++) {
    switch (condInfoPtrList[i]->condType()) {
    case CondInfo::WIRE: {
      meshCondPtrList[i] = 
	new MeshWire(dynamic_cast<Wire*>(condInfoPtrList[i]), 
		     useNonUniformMesh, useQuadPanel, aspectRatioWarning);
      break;
    }
    case CondInfo::RING: {
      meshCondPtrList[i] = 
	new MeshRing(dynamic_cast<Ring*>(condInfoPtrList[i]), 
		     useNonUniformMesh, useQuadPanel, aspectRatioWarning);
      break;
    }
    case CondInfo::SPIRAL: {
      meshCondPtrList[i] = 
	new MeshSpiral(dynamic_cast<Spiral*>(condInfoPtrList[i]), 
		       useNonUniformMesh, useQuadPanel, aspectRatioWarning);
      break;
    }
    case CondInfo::GROUND: {
      bool useQuad = true; // always use quad for ground
      meshCondPtrList[i] = 
	new MeshGround(dynamic_cast<Ground*>(condInfoPtrList[i]), 
		       useNonUniformMesh, useQuad, aspectRatioWarning);
      break;
    }
    case CondInfo::RECTSPIRAL: {
      bool useQuad = true; // always use quad for rectangular spiral
      meshCondPtrList[i] =
	new MeshRectSpiral(dynamic_cast<RectSpiral*>(condInfoPtrList[i]),
			   useNonUniformMesh, useQuad, aspectRatioWarning);
      break;
    }
    default: 
      surf::errorMessage ("mesh.cc::generateMesh", "Illegal condType in input file");
      break;
    }
  }

  countTotalNumPanelAndNode();
  setupVertexList();
  meshGenReport(useNonUniformMesh, useQuadPanel);
}

/**********************************************************************
 * ~Mesh --
 **********************************************************************/
Mesh::~Mesh (
	     void)
{
  for (size_t i = 0; i < meshCondPtrList.size(); i++) {
    delete meshCondPtrList[i];
  }
}

/**********************************************************************
 * meshGenReport --
 **********************************************************************/
void
Mesh::meshGenReport (
		     const bool useNonUniformMesh,
		     const bool useQuadPanel) const
{
  cout << endl << "\t Panels and nodes have been generated";
  cout << endl << "\t Number of conductors := " << numCond_;
  cout << endl 
       << "\t Number of panels := " << totalNumPanel_ 
       << ", number of nodes := " << totalNumNode_ 
       << endl;
  if (useNonUniformMesh) {
    cout << "\t Non-uniform mesh has been used" << endl;
  } else {
    cout << "\t Uniform mesh has been used" << endl;
  }
  if (useQuadPanel) {
    cout << "\t Only quadrilateral panels are used" << endl;
  } else {
    cout << "\t Both triangle panels and quadrilateral panels are used" << endl;
  }
}

/**********************************************************************
 * meshReadReport --
 **********************************************************************/
void
Mesh::meshReadReport (
		      const char* inputMeshFile) const
{
  cout << endl 
       << "\t Panels and nodes have been read from file  " << inputMeshFile << endl
       << "\t Number of conductors := " << numCond_ << endl 
       << "\t Number of panels := " << totalNumPanel_ 
       << ", number of nodes := " << totalNumNode_ 
       << endl;
}

/**********************************************************************
 * outputMesh --
 **********************************************************************/
void
Mesh::outputMesh (
		  const char* outputMeshFile,
		  const MeshFormat format)
{
  if (format == PATRAN) {
    outputMeshInPatranFormat(outputMeshFile);
  } else {
    outputMeshInFastCapFormat(outputMeshFile);
  }
}

/**********************************************************************
 * outputMeshInFastCapFormat --
 **********************************************************************/
void
Mesh::outputMeshInFastCapFormat (
				 const char* outputMeshFile) const
{
  std::ofstream fout(outputMeshFile);
  fout << "0 mesh in FastCap format output from FastImp" << endl;
  for (size_t condIndex = 0; condIndex < numCond_; condIndex++) {
    for (size_t panelIndex = 0; panelIndex < numPanel(condIndex); panelIndex++) {
      if (panelShape(condIndex, panelIndex) == QUAD) {
	//	fout << "Q " << condIndex << "  ";
	// It appears that in fastcap, condIndex starts at 1
	fout << "Q " << condIndex+1 << "  ";
      } else {
	fout << "T " << condIndex << "  ";
      }
      for (size_t vertexIndex = 0; 
	   vertexIndex < numVertex(condIndex, panelIndex); vertexIndex++) {
	point3D node = vertex(condIndex, panelIndex, vertexIndex);
	fout << node.x() * unit_ << " " 
	     << node.y() * unit_ << " " 
	     << node.z() * unit_ << "  ";
      }
      fout << endl;
    }
  }
  
  cout << endl << "\tmesh file in FastCap format " << outputMeshFile 
       << " has been generated" << endl;
}


/**********************************************************************
 * countTotalNumPanelAndNode --
 **********************************************************************/
void
Mesh::countTotalNumPanelAndNode (
				 void)
{
  totalNumPanel_ = 0;
  totalNumNode_ = 0;
  totalNumNonContactPanel_ = 0;
  for (size_t i = 0; i < numCond_; i++) {
    totalNumPanel_ += meshCondPtrList[i]->numPanel();
    totalNumNode_ += meshCondPtrList[i]->numNode();
    totalNumNonContactPanel_ += meshCondPtrList[i]->numNonContactPanel();
  }
}

/**********************************************************************
 * minStep --
 **********************************************************************/
double 
Mesh::minStep (
	       const size_t ci) const
{
  return meshCondPtrList[ci]->minStep(); 
}

/**********************************************************************
 * numNode --
 **********************************************************************/
int 
Mesh::numNode(
	      const size_t ci) const 
{ 
  return meshCondPtrList[ci]->numNode(); 
}

/**********************************************************************
 * numPanel --
 **********************************************************************/
int 
Mesh::numPanel (
		const size_t ci) const 
{
  return meshCondPtrList[ci]->numPanel(); 
}

/**********************************************************************
 * node --
 **********************************************************************/
pfft::point3D 
Mesh::node (
	    const size_t ci, 
	    const size_t i) const 
{ 
  return meshCondPtrList[ci]->node(i); 
}

/**********************************************************************
 * panelShape --
 **********************************************************************/
PanelShape 
Mesh::panelShape (
		  const size_t condIndex, 
		  const size_t i) const 
{ 
  return meshCondPtrList[condIndex]->panelShape(i); 
}

/**********************************************************************
 * panelType --
 **********************************************************************/
PanelType 
Mesh::panelType (
		 const size_t condIndex, 
		 const size_t i) const 
{ 
  return meshCondPtrList[condIndex]->panelType(i); 
}

/**********************************************************************
 * isNonContact --
 **********************************************************************/
bool 
Mesh::isNonContact(
		   const size_t condIndex, 
		   const size_t i) const 
{ 
  return meshCondPtrList[condIndex]->isNonContact(i); 
}

/**********************************************************************
 * isContact --
 **********************************************************************/
bool 
Mesh::isContact(
		   const size_t condIndex, 
		   const size_t i) const 
{ 
  return meshCondPtrList[condIndex]->isContact(i); 
}

/**********************************************************************
 * isNonBuffer --
 **********************************************************************/
bool 
Mesh::isNonBuffer(
		   const size_t condIndex, 
		   const size_t i) const 
{ 
  return meshCondPtrList[condIndex]->isNonBuffer(i); 
}

/**********************************************************************
 * isBuffer --
 **********************************************************************/
bool 
Mesh::isBuffer(
	       const size_t condIndex, 
	       const size_t i) const 
{ 
  return meshCondPtrList[condIndex]->isBuffer(i); 
}

/**********************************************************************
 * numVertex --
 **********************************************************************/
int Mesh::numVertex(
		    const size_t condIndex, 
		    const size_t i) const 
{ 
  return meshCondPtrList[condIndex]->numVertex(i); 
}

/**********************************************************************
 * globalNodeIndexOnOneCond --
 **********************************************************************/
int 
Mesh::globalNodeIndexOnOneCond(
			       const size_t condIndex, 
			       const size_t panelIndex, 
			       const size_t localIndex) const 
{ 
  return meshCondPtrList[condIndex]->globalNodeIndex(panelIndex, localIndex); 
}

/**********************************************************************
 * setupVertexList --
 * See vertex.h for comments
 **********************************************************************/
void
Mesh::setupVertexList (
		       void)
{
  vector<pfft::point3D> nodeList;
  nodeList.reserve(totalNumNode());
  for (size_t condIndex = 0; condIndex < numCond(); condIndex++) {
    for (size_t nodeIndex = 0; nodeIndex < numNode(condIndex);
	 nodeIndex++) {
      nodeList.push_back(node(condIndex, nodeIndex));
    }
  }

  vector<vector<size_t> > sharePanelIndexList(totalNumNode());
  vector<PanelType> typeList(totalNumNode(), NON_CONTACT);
  for (size_t condIndex = 0; condIndex < numCond(); condIndex++) {
    for (size_t panelIndex = 0; panelIndex < numPanel(condIndex);
	 panelIndex++) { 
      PanelType panType = panelType(condIndex, panelIndex);
      for (size_t li = 0; 
	   li < numVertex(condIndex, panelIndex); li++) {
	size_t gi = globalNodeIndex(condIndex, panelIndex, li);
	sharePanelIndexList[gi].push_back(panelIndex);
	if ((panType != NON_CONTACT) && (panType != BUFFER) )
	  typeList[gi] = panType;
      }
    }
  }
  vertexList.reserve(nodeList.size());
  for (size_t i = 0; i < nodeList.size(); i++) {
    vertexList.push_back(Vertex(nodeList[i], sharePanelIndexList[i], typeList[i]));
  }
}

/**********************************************************************
 * sharedPanelEdgeLength --
 **********************************************************************/
double 
Mesh::sharedPanelEdgeLength (
			     const size_t condIndex, 
			     const int panelIndex1, 
			     const int panelIndex2) const
{
  vector<pfft::point3D> sharedNodeList;
  for (int li1=0; li1 < numVertex(condIndex, panelIndex1) ; li1++) {
    int nodeIndex1 = globalNodeIndexOnOneCond(condIndex, panelIndex1, li1);
    for (int li2=0; li2 < numVertex(condIndex, panelIndex2) ; li2++) {
      int nodeIndex2 = globalNodeIndexOnOneCond(condIndex, panelIndex2, li2);
      if (nodeIndex1 == nodeIndex2) {
	sharedNodeList.push_back(node(condIndex, nodeIndex2));
      }
    }
  }

  if (sharedNodeList.size() != 2) {
    surf::errorMessage("mesh.cc : sharedPanelEdgeLength",
		       "Bug : Number of shared nodes between two adjacent panels is not 2");
  }

  return length(sharedNodeList[0] - sharedNodeList[1]);
}

/**********************************************************************
 * outputMeshInPatranFormat --
 **********************************************************************/
void
Mesh::outputMeshInPatranFormat (
				const char* outputMeshFile) const
{
  std::ofstream fout(outputMeshFile);
  fout << "# mesh in patran format output from fastimp" << endl;
  fout << "c " << numCond_ << endl;
  fout << "n " << totalNumNode_ << endl;
  fout << "p " << totalNumPanel_ << endl;
  fout << "u " << unit_ << endl;

  for (size_t condIndex = 0; condIndex < numCond_; condIndex++) {
    for (size_t nodeIndex = 0; nodeIndex < numNode(condIndex); nodeIndex++) {
      fout << "v " << condIndex << " " << nodeIndex << " " 
	   << std::setprecision(10)
	   << node(condIndex, nodeIndex).x() << " "
	   << node(condIndex, nodeIndex).y() << " "
	   << node(condIndex, nodeIndex).z()
	   << std::setprecision(6)
	   << endl;
    }
  }

  for (size_t condIndex = 0; condIndex < numCond_; condIndex++) {
    for (size_t panelIndex = 0; panelIndex < numPanel(condIndex); panelIndex++) {
      if (panelShape(condIndex, panelIndex) == QUAD) {
	fout << "Q " << condIndex << "  ";
      } else {
	fout << "T " << condIndex << "  ";
      }
      if (panelType(condIndex, panelIndex) == NON_CONTACT) {
	fout << "NC ";
      } else if (panelType(condIndex, panelIndex) == LEFT_CONTACT) {
	fout << "LC ";
      } else if (panelType(condIndex, panelIndex) == RIGHT_CONTACT) {
	fout << "RC ";
      } else if (panelType(condIndex, panelIndex) == CONTACT) {
	fout << "C ";
      } else if (panelType(condIndex, panelIndex) == BUFFER) {
	fout << "B ";
      }
      for (size_t vertexIndex = 0; 
	   vertexIndex < numVertex(condIndex, panelIndex); 
	   vertexIndex++) {
	fout << globalNodeIndexOnOneCond(condIndex, panelIndex, vertexIndex) 
	     << " ";
      }
      fout << endl;
    }
  }
  
  cout << endl << "\t patran mesh file " << outputMeshFile 
       << " has been generated" << endl;
}

/**********************************************************************
 * readMeshFile --
 * read the mesh file in patran format, as the following:
   # This is the mesh file for   <--- title line started with # or *
   * generated on 05/29/02        <--- comments, started with # or *
   c 1                           <--- number of conductors, started with c or C
   n 26                          <--- number of nodes, started with n or N
   p 24                          <--- number of panels, started with p or P
   u 1e-3                        <--- unit, mm, started with u or U
   v 1 11 x1 y1 z1               <--- condIndex, localNodeIndex, coordinates
   ...
   v 2 21 x2 y2 z2               <--- condIndex, localNodeIndex, coordinates
        localNodeIndex: the index on each conductor starts from zero
   Q 1 LC 11 19 21 22            <--- panel description
   Q 1 RC 11 19 21 22            <--- panel description
   T 2 NC 11 19 21               <--- panel description
   Q 5 B 11 19 21               <--- panel description
        Q, T :       quad or trangle
        1, 2 :       index of the conductor where the panel lies
        LC, RC, NC, B:  left contact, right contact, non-contact or buffer 
        11, 19, 21 : index of the vertices
 *
 * Note:
 * 1) function Mesh::outputMeshInPatranFormat() could output this format
 * 2) Because of the way I compute the current, it is mandatory that
 *    there should be two layers of quad panels around each contact. 
 *    These panels are also the non-contact ones. But they must be quadrilateral.
 **********************************************************************/
void
Mesh::readMeshFile (
		    const std::vector<CondInfo*>& condInfoPtrList,
		    char* inputMeshFile)
{
  FILE* fd = fopen(inputMeshFile, "r");
  if (!fd) {
    errorMessage("Mesh::readMeshFile()",
		 "Can not open the mesh file! Please double check its name.");
  }

  readMeshHead (condInfoPtrList, fd);
  int estimatedNumNodePerCond = totalNumNodeFromMeshFile_ / numCond_;
  int estimatedNumPanelPerCond = totalNumPanelFromMeshFile_ / numCond_;
  meshCondPtrList.reserve(numCond_); 
  for (size_t i=0; i<numCond_; i++) {  
    meshCondPtrList[i] = new MeshOtherShape(estimatedNumNodePerCond,
					    estimatedNumPanelPerCond);
  }
  numNonContactPanelFromMeshFile_.resize(numCond_);
  
  while (readNextLine(fd)) {
    switch (nextLine_[0]) {
    case 'v':
    case 'V':
      {
	point3D node;
	int condIndex, nodeIndex;
	readNode(&condIndex, &nodeIndex, node);
	//	meshCondPtrList[condIndex]->addNode(node);
	meshCondPtrList[condIndex]->addNode(node, nodeIndex);
      }
      break;
    case 'q':
    case 'Q':
      {
	Panel panel;
	int condIndex;
	readQuadPanel(&condIndex, panel);
	meshCondPtrList[condIndex]->addPanel(panel);
      }
      break;
    case 't':
    case 'T':
      {
	Panel panel;
	int condIndex;
	readTrianglePanel(&condIndex, panel);
	meshCondPtrList[condIndex]->addPanel(panel);
      }
      break;
    default:
      surf::errorMessage("mesh.cc::readMeshFile()",
			 "Illegal token in the mesh file");
      break;
    }
  }
  fclose(fd);

  for (size_t i=0; i<numCond_; i++) {
    meshCondPtrList[i]->setNumNonContactPanel(numNonContactPanelFromMeshFile_[i]);
    meshCondPtrList[i]->setNumNode();
    meshCondPtrList[i]->setNumPanel();
  }

  countTotalNumPanelAndNode();
  if (totalNumPanel_ != totalNumPanelFromMeshFile_) {
    surf::errorMessage("mesh.cc::readMeshFile()",
		       "Number of panels in the mesh file head is inconsistent with that in the mesh file body");
  }

  if (totalNumNode_ != totalNumNodeFromMeshFile_) {
    cout << endl << "totalNumNode :=" << totalNumNode_ << endl;
    cout << "totalNumNodeFromMeshFile :=" << totalNumNodeFromMeshFile_ << endl;
    surf::errorMessage("mesh.cc::readMeshFile()",
		       "Number of nodes in the mesh file head is inconsistent with that in the mesh file body");
  }

  setupVertexList();
  meshReadReport(inputMeshFile);
}

/**********************************************************************
 * readMeshHead --
 * read numCond, totalNumNode, totalNumPanel and unit
 **********************************************************************/
void
Mesh::readMeshHead (
		    const std::vector<CondInfo*>& condInfoPtrList,
		    FILE* fd)
{
  int counter = 0;
  int numItems = 4;

  while (counter < numItems) {
    readNextLine(fd);
    switch (nextLine_[0]) {
    case 'c':
    case 'C':
      counter++;
      readInt(&numCond_);
      break;
    case 'n':
    case 'N':
      counter++;
      readInt(&totalNumNodeFromMeshFile_);
      break;
    case 'p':
    case 'P':
      counter++;
      readInt(&totalNumPanelFromMeshFile_);
      break;
    case 'u':
    case 'U':
      counter++;
      readDouble(&unit_);
      break;
    default:
      surf::errorMessage("mesh.cc::readMeshHead()",
			 "Illegal token in the head of mesh file");
      break;
    } 
  }

  if (numCond_ != condInfoPtrList.size()) {
    surf::errorMessage("mesh.cc::readMeshHead()", 
		       "number of conductors in mesh file and structure file is inconsistent");
  }

  if (unit_ != condInfoPtrList[0]->unit()) {
    surf::errorMessage("mesh.cc::readMeshHead()", 
		       "unit in mesh file and structure file is inconsistent");
  }
}

/**********************************************************************
 * readNextLine --
 **********************************************************************/
bool
Mesh::readNextLine(
		   FILE* fd)
{
  static char lineBfr[300];

  do {
    /* Read line from file, return NULL if eof */
    nextLine_ = fgets(lineBfr, sizeof(lineBfr)-1, fd);
    if (nextLine_ == NULL) return false;
    assert(nextLine_ == lineBfr);

    /* Remove EOL and comments from the end of line buffer */
    char *p1;
    for(p1=nextLine_;
	(*p1)!= '\0' && (*p1) != '#' && (*p1) != '*' && (*p1) != '\n'; 
	p1++) { }
    *p1 = '\0';

    /* Skip leading spaces */
    while((*nextLine_) && isspace(*nextLine_)) nextLine_++;

    /* Remove trailing spaces */
    p1 = nextLine_ + strlen(nextLine_);
    while (p1 > nextLine_ && isspace(p1[-1])) p1--;
    *p1 = '\0';

  } while (strlen(nextLine_) == 0);

  return true;
}
	         
/**********************************************************************
 * readInt --
 **********************************************************************/
void
Mesh::readInt(
	      int* val) const 
{     
  char token[20];
  sscanf(nextLine_, "%s %d", token, val);
}

/**********************************************************************
 * readDouble --
 **********************************************************************/
void
Mesh::readDouble(
		 double* val) const 
{     
  char token[20];
  sscanf(nextLine_, "%s %lf", token, val);
}

/**********************************************************************
 * readNode --
 **********************************************************************/
void
Mesh::readNode (
		int* condIndex,
		int* nodeIndex,
		point3D& point) const 
{  
  char token[20];
  double x, y, z;
  sscanf(nextLine_, "%s %d %d %lf %lf %lf", 
	 token, condIndex, nodeIndex, &x, &y, &z);
  point = point3D(x, y, z);
}

/**********************************************************************
 * readQuadPanel --
 **********************************************************************/
void
Mesh::readQuadPanel (
		     int* condIndex,
		     Panel& panel)
{  
  char token[20], typeToken[3];
  vector<int> nodeIndex(4);
  PanelType type;
  sscanf(nextLine_, "%s %d %s %d %d %d %d", token, condIndex, typeToken, 
	 &nodeIndex[0], &nodeIndex[1], &nodeIndex[2], &nodeIndex[3]);
  if (!(strcmp(typeToken, "NC") && strcmp(typeToken, "nc"))) {
    type = NON_CONTACT;
    numNonContactPanelFromMeshFile_[*condIndex] ++;
  } else if (!(strcmp(typeToken, "RC") && strcmp(typeToken, "rc"))) {
    type = RIGHT_CONTACT;
  } else if (!(strcmp(typeToken, "LC") && strcmp(typeToken, "lc"))) {
    type = LEFT_CONTACT;
  } else if (!(strcmp(typeToken, "C") && strcmp(typeToken, "c"))) {
    type = CONTACT;
  } else if (!(strcmp(typeToken, "B") && strcmp(typeToken, "b"))) {
    type = BUFFER;
  } else {
    surf::errorMessage("mesh.cc::readQuadPanel()",
		       "Illegal panel type in mesh file");
  }
  panel = Panel(nodeIndex, type, QUAD);
}

/**********************************************************************
 * readTrianglePanel --
 **********************************************************************/
void
Mesh::readTrianglePanel (
			 int* condIndex,
			 Panel& panel)
{  
  char token[20], typeToken[3];
  vector<int> nodeIndex(3);
  PanelType type;
  sscanf(nextLine_, "%s %d %s %d %d %d", token, condIndex, typeToken, 
	 &nodeIndex[0], &nodeIndex[1], &nodeIndex[2]);
  if (!(strcmp(typeToken, "NC") && strcmp(typeToken, "nc"))) {
    type = NON_CONTACT;
    numNonContactPanelFromMeshFile_[*condIndex] ++;
  } else if (!(strcmp(typeToken, "RC") && strcmp(typeToken, "rc"))) {
    type = RIGHT_CONTACT;
  } else if (!(strcmp(typeToken, "LC") && strcmp(typeToken, "lc"))) {
    type = LEFT_CONTACT;
  } else if (!(strcmp(typeToken, "C") && strcmp(typeToken, "c"))) {
    type = CONTACT;
  } else if (!(strcmp(typeToken, "B") && strcmp(typeToken, "b"))) {
    type = BUFFER;
  } else {
    surf::errorMessage("mesh.cc::readQuadPanel()",
		       "Illegal panel type in mesh file");
  }
  panel = Panel(nodeIndex, type, TRIANGLE);
}

