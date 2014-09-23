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

  const static char cvsid[] = "$Id: mesh.h,v 1.7 2002/12/01 15:50:44 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _MESH_H_
#define _MESH_H_

#include <stdio.h> 
#include <vector>
#include "vector3D.h"
#include "condInfo.h"
#include "wire.h"
#include "spiral.h"
#include "ring.h"
#include "ground.h"
#include "meshConst.h"
#include "vertex.h"
#include "panel.h"
#include "service.h"

namespace mesh {

  class MeshCond;

  class Mesh {

  public:
    Mesh(void) : unit_(1.) {}

    // There is no constructor. I have a vector<pointer>
    // as the private member. Construction-assign-detruction causes
    // two pinters pointing to the same location and that location
    // has already be invalidated by deallocating one of the two pointers.
    // Using function rather than constructor, I have avoided the copying.
    void generateMesh(const std::vector<surf::CondInfo*>& condInfoPtrList, 
		      const bool useNonUniformMesh, const bool useQuadPanel,
		      const bool aspectRatioWarning);
    void readMeshFile(const std::vector<surf::CondInfo*>& condInfoPtrList, 
		      char* inputMeshFile);
    void outputMesh(const char* outputMeshFile, const MeshFormat format);
    ~Mesh(void);

    double unit(void) { return unit_; }
    int numNode(const size_t condIndex) const;
    int numPanel(const size_t condIndex) const;
    double minStep(const size_t condIndex) const;
    pfft::point3D node(const size_t condIndex, 
		       const size_t globalNodeIndexOnOneCond) const;
    std::vector<pfft::point3D> vertex(const size_t condIndex, 
				      const size_t panelIndex) const {
      std::vector<pfft::point3D> vertexList;
      for (int vi = 0; vi < numVertex(condIndex, panelIndex); vi++)
	vertexList.push_back(vertex(condIndex, panelIndex, vi));
      return vertexList;
    }
    pfft::point3D vertex(const size_t condIndex, 
			 const size_t panelIndex,
			 const size_t localVertexIndex) const {
      return node(condIndex, globalNodeIndexOnOneCond(condIndex, panelIndex, 
						      localVertexIndex));
    }
    PanelShape panelShape(const size_t condIndex, const size_t panelIndex) const;
    PanelType panelType(const size_t condIndex, const size_t panelIndex) const;
    bool isNonContact(const size_t condIndex, const size_t panelIndex) const;
    bool isContact(const size_t condIndex, const size_t panelIndex) const;
    bool isNonBuffer(const size_t condIndex, const size_t panelIndex) const;
    bool isBuffer(const size_t condIndex, const size_t panelIndex) const;
    int numVertex(const size_t condIndex, const size_t panelIndex) const;
    int globalNodeIndexOnOneCond(const size_t condIndex, 
				 const size_t panelIndex, 
				 const size_t localIndex) const;
    int globalNodeIndex(const size_t condIndex, 
			const size_t panelIndex, 
			const size_t localIndex) const { 
      int globalNodeIndex = 0;
      for (size_t i=0; i < condIndex; i++)
	globalNodeIndex += numNode(i);
      return globalNodeIndex + globalNodeIndexOnOneCond(condIndex, panelIndex, 
							localIndex); 
    }
    int globalNodeIndex(const size_t condIndex, 
			const size_t nodeIndexOnOneCond) const {
      int globalNodeIndex = 0;
      for (size_t i=0; i < condIndex; i++) 
	globalNodeIndex += numNode(i);
      return globalNodeIndex + nodeIndexOnOneCond;
    }
    int globalPanelIndex(const size_t condIndex, const size_t panelIndex) const {
      int globalPanelIndex = 0;
      for (size_t i = 0; i < condIndex; i++)
	globalPanelIndex += numPanel(i);
      return globalPanelIndex + panelIndex;
    }

    int numCond(void) const { return numCond_; }    
    int totalNumNode(void) const { return totalNumNode_; }    
    int totalNumPanel(void) const { return totalNumPanel_; }
    int totalNumNonContactPanel(void) const { return totalNumNonContactPanel_; }
    int totalNumContactPanel(void) const { 
      return totalNumPanel_ - totalNumNonContactPanel_; }

    PanelType vertexType (size_t index) const { 
      return vertexList[index].type(); }
    pfft::point3D vertexPos (size_t index) const {
      return vertexList[index].coordinate();  }
    size_t numSharedPanel (size_t index) const {
      return vertexList[index].numSharePanel(); }
    size_t sharedPanelIndex (size_t vertexIndex,size_t subPanelIndex) const {
      return vertexList[vertexIndex].sharePanelIndex(subPanelIndex); }
    double sharedPanelEdgeLength (const size_t condIndex, 
				  const int panelIndex1, 
				  const int panelIndex2) const;
    
  private:
    double unit_;
    std::vector<MeshCond*> meshCondPtrList;
    std::vector<Vertex> vertexList;
    int numCond_;
    int totalNumNode_;
    int totalNumPanel_;
    int totalNumNonContactPanel_;
    char* nextLine_;
    int totalNumNodeFromMeshFile_;
    int totalNumPanelFromMeshFile_;
    std::vector<int> numNonContactPanelFromMeshFile_;

    void countTotalNumPanelAndNode(void);
    void meshGenReport(const bool useNonUniformMesh, const bool useQuadPanel) const;
    void meshReadReport(const char* inputMeshFile) const;
    void outputMeshInFastCapFormat(const char* outputMeshFile) const;
    void outputMeshInPatranFormat(const char* outputMeshFile) const;
    bool readNextLine(FILE* fd);
    void readInt(int* val) const;
    void readDouble(double* val) const;
    void readNode (int* condIndex, int* nodeIndex, pfft::point3D& point) const;
    void readQuadPanel (int* condIndex, Panel& panel);
    void readTrianglePanel (int* condIndex, Panel& panel);
    void readMeshHead (const std::vector<surf::CondInfo*>& condInfoPtrList, 
		       FILE* fd);

    void setupVertexList(void);
  };

} //namespace mesh

#endif
