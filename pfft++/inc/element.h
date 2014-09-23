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

  const static char cvsid[] = "$Id: element.h,v 1.5 2002/09/24 16:01:07 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include <vector>
#include <string>
#include <iostream>
#include "vector3D.h"

namespace pfft {

class element {
  
  friend bool operator == (const element&, const element& );
  friend bool operator != (const element&, const element& );
  friend std::istream& operator>> (std::istream&, element&);

public:
  typedef std::vector<point3D>::size_type VertexIndex;

  element(void) : name_(0), boundaryIndex_(0), boundaryType_(1) {;}
  element(
	  const point3D& v1,
	  const point3D& v2,
	  const point3D& v3,
	  const int name = 0, 
	  const int boundaryIndex = 0,
	  const int boundaryType = 1);
  element(
	  const point3D& v1,
	  const point3D& v2,
	  const point3D& v3,
	  const point3D& v4,
	  const int name = 0, 
	  const int boundaryIndex = 0,
	  const int boundaryType = 1);
  element(
	  const std::vector<point3D>& vertex,
	  const int name = 0, 
	  const int boundaryIndex = 0,
	  const int boundaryType = 1);

  const vector3D<double> tangent1(void) const { return tangent_1; }
  const vector3D<double> tangent2(void) const { return tangent_2; }
  const vector3D<double> normal(void) const { return normal_; }
  const double tangent1(const size_t i) const { 
    if (i == 0) return tangent_1.x();
    else if (i == 1) return tangent_1.y();
    else return tangent_1.z();
  }
  const double tangent2(const size_t i) const { 
    if (i == 0) return tangent_2.x();
    else if (i == 1) return tangent_2.y();
    else return tangent_2.z();
  }
  const double normal(const size_t i) const { 
    if (i == 0) return normal_.x();
    else if (i == 1) return normal_.y();
    else return normal_.z();
  }
  const double linearShapeFuncCoe(const size_t vertexIndex, 
				  const size_t coeIndex) const { 
    if (coeIndex == 0) return linearShapeFuncList_[vertexIndex].a;
    else if (coeIndex == 1) return linearShapeFuncList_[vertexIndex].b;
    else return linearShapeFuncList_[vertexIndex].c;
  }

  void shiftToLocalCoord(const point3D&);

  const point3D centroid(void) const { return centroid_; }
  const point3D boundingSphereCenter(void) const { return centroid_; }
  const double boundingSphereRadius(void) const { return boundingSphereRadius_; }
  const double area(void) const { return area_; }
  // why this does not compile ??
  //  const point3D vertex(const size_t i) const { return vertex_.at(i); }  
  const point3D vertex(const VertexIndex i) const { return vertex_[i]; }
  const size_t name(void) const { return name_; }
  const size_t shape(void) const { return vertex_.size(); }
  const size_t numVertex(void) const { return vertex_.size(); }
  const size_t boundaryIndex(void) const { return boundaryIndex_; }
  const short int boundaryType(void) const { return boundaryType_; }
  const size_t index(void) const { return name_; }
  const double edgeMidPointDist1(void) const { return edgeMidPointDist1_; }
  const double edgeMidPointDist2(void) const { return edgeMidPointDist2_; }

  void setName(const size_t name) { name_ = name; } 
  void setBoundaryIndex(const size_t boundaryIndex) { boundaryIndex_ = boundaryIndex; }
  void setBoundaryType(const short int boundaryType) { boundaryType_ = boundaryType; }
  double findMinEdgeLength(void) const;

private:
  //
  // essential info
  //
  size_t name_;  // use the element index as its name
  size_t boundaryIndex_; // index of the conductor where the element lies on
  short int boundaryType_; // Dirichlet:1, Neumann:2, others:0
  std::vector<point3D> vertex_;

  //
  // derived info
  //

  // geometry
  vector3D<double> tangent_1; // served as x axis in local coordinate system 
  vector3D<double> tangent_2; // served as y axis in local coordinate system 
  vector3D<double> normal_;   // served as z axis in local coordinate system 
  double area_;
  point3D centroid_;
  double boundingSphereRadius_; // max dist between vertex and center
  double edgeMidPointDist1_; // only meaningful for quad panels.
  double edgeMidPointDist2_; // the distance between midpoint of opposite edges

  // basis function on th element
  typedef struct LinearShapeFunc {
    // y(t1, t2) = a*t1 + b*t2 + c, t1 and t2 are local coordinates
    // y(t11, t12) = 1, y(t21, t22) = 0, y(t31, t32) = 0,
    // where (ti1, ti2) are local coordinates of triangle vertices
    double a;
    double b;
    double c;    
  } LinearShapeFunc;
  std::vector<LinearShapeFunc> linearShapeFuncList_;

  void setupElement();
  void compCentroidAndArea();
  void compNormalAndTangent();
  void compBoundingSphereRadius();
  void compShapeFunction();
  double compTriangleArea(const point3D& v1, const point3D& v2, const point3D& v3);
};

std::ostream& operator << (std::ostream& os, const element& e);
bool readElementFile (const char* fileName, std::vector<element>& elementList, 
		      std::string& title, size_t&);
void readElementFile (const char* fileName, std::vector<element>& elementList);

} //namespace pfft 

#endif





