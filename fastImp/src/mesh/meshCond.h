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

  const static char cvsid[] = "$Id: meshCond.h,v 1.11 2003/05/13 21:10:50 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _MESH_COND_H_
#define _MESH_COND_H_

#include <vector>
#include "condInfo.h"
#include "vector3D.h"
#include "panel.h"
#include "meshConst.h"
#include "vec.h"
#include "service.h"
#include <cfloat> // for DBL_MAX

namespace mesh {

  const double PI = 3.14159265358979;

  class MeshCond {

  public:

    MeshCond(void) { }
    MeshCond(const surf::CondInfo* condInfoPtr, 
	     const bool useNonUniformMesh,
	     const bool useQuadPanel,
	     const bool aspectRatioWarning) {}
    MeshCond(const int numNode, const int numPanel) : actualNumNode_(0) { 
      //      nodeList.reserve(numNode);
      nodeList.resize(numNode);
      panelList.reserve(numPanel);
    }
    virtual ~MeshCond(void) {}

    double minStep(void) { 
      if (useNonUniformMesh_) {
	return std::min(hx.min(), hz.min()); 
      } else {
	return minUniformStep;
      }
    }

    int numNode(void) const { return nodeList.size(); }
    int numPanel(void) const { return panelList.size(); }
    int numNonContactPanel(void) const { return numNonContactPanel_; }
    pfft::point3D node(size_t i) const { return nodeList[i]; }
    PanelShape panelShape(size_t i) const { return panelList[i].shape(); }
    PanelType panelType(size_t i) const { return panelList[i].type(); }
    bool isNonContact(size_t i)  const { return panelList[i].isNonContact(); }
    bool isContact(size_t i)  const { return panelList[i].isContact(); }
    bool isNonBuffer(size_t i)  const { return panelList[i].isNonBuffer(); }
    bool isBuffer(size_t i)  const { return panelList[i].isBuffer(); }
    int numVertex(size_t i) const { return panelList[i].numVertex(); }
    int globalNodeIndex(size_t i, size_t localIndex) const { 
      return panelList[i].globalNodeIndex(localIndex); }

    void setNumNode(void) {nodeList.resize(actualNumNode_); numNode_ = actualNumNode_; }
    void setNumPanel(void) {numPanel_ = panelList.size();}
    void setNumNonContactPanel(const int numNonContactPanel) {
      numNonContactPanel_ = numNonContactPanel;}
    //    void addNode(const pfft::point3D& node) {nodeList.push_back(node);}
    // I can not just push_back, the nodeIndex here serves as the global index too.
    void addNode(const pfft::point3D& node, size_t nodeIndex) {
      if (nodeIndex >= nodeList.size()) {
	nodeList.resize(nodeIndex+1);
      }
      nodeList[nodeIndex] = node; 
      actualNumNode_ ++;
    }
    void addPanel(const Panel& panel) {panelList.push_back(panel);}

    void findNonUniformStep(const double length, const int numStep, 
			    TNT::Vector<double>& h);

  protected:
    int numNode_;
    int actualNumNode_;
    int numPanel_;
    int numNonContactPanel_;
    std::vector<pfft::point3D> nodeList;
    std::vector<Panel> panelList;
    TNT::Vector<double> hx;
    TNT::Vector<double> hz;
    double minUniformStep;
    bool useNonUniformMesh_;
    bool useQuadPanel_;
    bool aspectRatioWarning_;

    // global coordinate system
    pfft::vector3D<double> origin;
    pfft::vector3D<double> X;
    pfft::vector3D<double> Y;
    pfft::vector3D<double> Z;

    // every shape is transfered to a cube and mesh is generated on it
    typedef struct Cube {
      double x1;
      double x2;
      int numPanelX;
      double y1;
      double y2;
      int numPanelY;
      double z1;
      double z2;
      int numPanelZ;
    } Cube;
    Cube cube;

    void meshCube(const surf::CondInfo* condInfoPtr);
    void findPanel(const surf::CondInfo* condInfoPtr);
    void findNonUniformNode(void);
    void findUniformNode(void);
    double findBuferStepSize(void);
    void transferToGlobalCoord(void);
    void checkAspectRatio(const TNT::Vector<double>& hx, 
			  const TNT::Vector<double>& hz, 
			  const double hy);

  };

} //namespace mesh

#endif
