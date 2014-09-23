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

  const static char cvsid[] = "$Id: meshRectSpiral.h,v 1.3 2002/12/01 15:50:06 zhzhu Exp $";

  ==========================================================================
*/
#ifndef __MESH_RECT_SPIRAL_H_
#define __MESH_RECT_SPIRAL_H_

#include "rectSpiral.h"
#include "meshCond.h"
#include "wirePath.h"
#include "lineSegment.h"
#include "corner.h"

namespace mesh {

  class MeshRectSpiral : public MeshCond {

  public:
    MeshRectSpiral(void) {}
    MeshRectSpiral(
		   const surf::RectSpiral* p,
		   const bool useNonUniformMesh,
		   const bool useQuadPanel,
		   const bool aspectRatioWarning);
  
  private:

    const surf::RectSpiral* p_;

    // -- internal working objects --
    WirePath wirePath_;
    LineSegmentList lineSegmentList_;
    CornerList cornerList_;

    size_t numNodeSofar_;

    std::vector<size_t> numPanelYOnLineSegArray_;

    // used for back check 
    std::vector<size_t> numNodeOnLineSegArray_;
    std::vector<size_t> numNodeOnCornerArray_;
    std::vector<size_t> numPanelOnLineSegArray_;
    std::vector<size_t> numPanelOnCornerArray_;

    typedef struct {
      std::vector<size_t> bottom;
      std::vector<size_t> back;
      std::vector<size_t> top;
      std::vector<size_t> front;
    } Map4OneContour; // put the global node index

    typedef std::vector<Map4OneContour> Map4ContourList;

    typedef struct {
      std::vector<std::vector<size_t> > right;
      std::vector<std::vector<size_t> > left;
    } Map4TwoEnds;

    Map4ContourList map4ContourList;
    //std::vector<Map4OneContour> map4ContourList;
    Map4TwoEnds map4TwoEnds;
    Map4OneContour map4LastContour;

  public:

    // -- mesh generator --
    void initMesh (void);
    void genMesh(void);
    void getMap4StartContour(void);

    void calcNumNodePanelOnEachCorner (const size_t);
    void calcNumNodePanelOnEachLineSeg (const size_t);

    void meshLineSegment (const size_t);
    void meshCorner (const size_t);

    void findNodeOnLineSegment (const size_t);
    void findNodeOnCorner (const size_t);
    void findPanelOnLineSegment (const size_t);
    void findPanelOnCorner (const size_t);

    void setup(void); 
    void setupWirePath(void);
    void setupLineSegmentList(void);
    void setupCornerList(void);

    void genFHOutput (void);

    pfft::point3D local2GlobalCoordSys (
					const pfft::point3D& origin,
					const pfft::point3D& X,
					const pfft::point3D& Y,
					const pfft::point3D& Z,
					double x, 
					double y, 
					double z);
  };
}
    
#endif
