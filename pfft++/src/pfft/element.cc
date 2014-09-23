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

  ==========================================================================
*/

const static char cvsid[] = "$Id: element.cc,v 1.7 2003/02/11 03:10:26 zhzhu Exp $";

#include <fstream>
#include <string>
#include <stdexcept> //for exception handling
#include "element.h"
#include <cfloat> // for DBL_EPSILON
#include "utils.h" //for errorMessage()

using namespace std;
using namespace pfft;

/**********************************************************************
 * element --
 **********************************************************************/
element::element (
		  const point3D& v1,
		  const point3D& v2,
		  const point3D& v3,
		  const int name, 
		  const int boundaryIndex,
		  const int boundaryType)
  : name_(name), boundaryIndex_(boundaryIndex), boundaryType_(boundaryType)
{
  vertex_.push_back(v1);
  vertex_.push_back(v2);
  vertex_.push_back(v3);
  setupElement();
}

element::element (
		  const point3D& v1,
		  const point3D& v2,
		  const point3D& v3,
		  const point3D& v4,
		  const int name, 
		  const int boundaryIndex,
		  const int boundaryType)
  : name_(name), boundaryIndex_(boundaryIndex), boundaryType_(boundaryType)
{
  vertex_.push_back(v1);
  vertex_.push_back(v2);
  vertex_.push_back(v3);
  vertex_.push_back(v4);
  setupElement();
}

element::element (
		  const vector<point3D>& vertex,
		  const int name, 
		  const int boundaryIndex,
		  const int boundaryType)
  : name_(name), boundaryIndex_(boundaryIndex), boundaryType_(boundaryType)
  , vertex_(vertex)
{
  setupElement();
}

/**********************************************************************
 * setupElement --
 **********************************************************************/
void
element::setupElement (
		       void)
{
  compCentroidAndArea();
  compNormalAndTangent();
  compBoundingSphereRadius();
  compShapeFunction(); // this function call must be after others
}

/**********************************************************************
 * compCentroidAndArea --
 * The standard algorithm to calculate the area of a ploygon is to break it
 * up into triangles and sum up the area of each triangle.
 *
 * The standard algorithm to calculate the centroid, the center of gravity,
 * is to compute the centroid of each triangle and take the weighted average
 * of them. The weighting factor is the area of each triangle.
 *
 * For the sake of efficiency, these two functions are combined here.
 *
 * It is assumed here that the logically adjacent vertices are 
 * physically adjacent.
 * For example, vertex[0] is physically connected to vertex[1], 
 * and vertex[1] to vertex[2], etc
 **********************************************************************/
void
element::compCentroidAndArea (
			      void)
{
  double triangleArea;
  point3D triangleCentroid;
  area_ = 0;
  centroid_ = vector3D<double>(0., 0., 0.);
  for (size_t i=1; i < shape()-1; i++) {
    triangleArea = compTriangleArea(vertex_[0], vertex_[i], vertex_[i+1]);
    triangleCentroid = (vertex_[0] + vertex_[i] + vertex_[i+1]) / 3.;
    area_ += triangleArea;
    centroid_ += triangleCentroid * triangleArea;
  }
  if (area_==0) {
    pfft::errorMessage("element.cc:compCentroidAndArea",
		       "Panel area is zero! There must be a bug here!!");
  }    
  centroid_ /= area_;
}

/**********************************************************************
 * compTriangleArea --
 **********************************************************************/
double
element::compTriangleArea (
			   const point3D& v1,
			   const point3D& v2,
			   const point3D& v3)
{
  vector3D<double> edge1 = v2 - v1;
  vector3D<double> edge2 = v2 - v3;
  vector3D<double> vec = crossProd(edge1, edge2);
  return 0.5 * length(vec);
}

/**********************************************************************
 * compNormalAndTangent --
 * Order of verteices and normal vector follow the left-hand rule.
 * Vectors tangent1, tangent2 and normal follow the right-hand rule 
 * and could be used as the local cartesian coordinate system. 
 **********************************************************************/
void
element::compNormalAndTangent (
			       void)
{
  if (shape() == 4) {
    // special handling here is purely for my own application's sake
    point3D midPoint1 = (vertex_[0] + vertex_[1]) / 2.;
    point3D midPoint2 = (vertex_[2] + vertex_[3]) / 2.;
    tangent_1 = midPoint2 - midPoint1;
    edgeMidPointDist1_ = length(tangent_1);
    
    midPoint1 = (vertex_[0] + vertex_[3]) / 2.;
    midPoint2 = (vertex_[1] + vertex_[2]) / 2.;
    tangent_2 = midPoint2 - midPoint1;
    edgeMidPointDist2_ = length(tangent_2);

    if ((edgeMidPointDist1_ < 1e-16) || (edgeMidPointDist2_ < 1e-16) ) {
      std::cout << endl
		<< "\t element.cc:compNormalAndTangent" << endl
		<< "\t edgeMidPointDist is zero! There must be a bug here!!" 
		<< endl;
      throw domain_error("edgeMidPointDist is zero! There must be a bug here!!");      
    }

  } else {
    tangent_1 = vertex_[1] - vertex_[0];
    tangent_2 = vertex_[1] - vertex_[2];
  }
  tangent_1.normalize();
  normal_ = crossProd(tangent_1, tangent_2);
  normal_.normalize();
  tangent_2 = crossProd(normal_, tangent_1);
  tangent_2.normalize();
}

/**********************************************************************
 * compBoundingSphereRadius --
 **********************************************************************/
void
element::compBoundingSphereRadius (
				   void)
{
  boundingSphereRadius_ = 0.;
  for (size_t i=0; i<shape(); i++) {
    vector3D<double> vec = boundingSphereCenter() - vertex_[i];
    boundingSphereRadius_ = max(boundingSphereRadius_, length(vec));
  }
}

/**********************************************************************
 * compShapeFunction --
 * y(t1, t2) = a*t1 + b*t2 + c, t1 and t2 are local coordinates
 * y(t11, t12) = 1, y(t21, t22) = 0, y(t31, t32) = 0,
 * where (ti1, ti2) are local coordinates of triangle vertices.
 *
 * I could pick (t11, t12) as the origin of the local coordinates,
 * hence t11 = t12 = 0 and c=1
 * Solving
 * | t21  t22 | |a| = | -1 |
 * | t31  t32 | |b|   | -1 |
 * gives us one shape funtion associated with one vertex.
 *
 * Note:
 * For each shape function, the origin is different. But the coordinate axises
 * are always the same panel tangent and normal vectors.
 *
 * The other option is to fix the origin as the centroid. This way, all
 * three coefficients are to be computed. This might be more convenient
 * for panel integration. The difference between these two options is just the
 * coefficient c.
 **********************************************************************/
void
element::compShapeFunction (
			    void)
{
  if (shape() != 3) return; // shape function is only for triangle elements

  linearShapeFuncList_.resize(3);
  for (size_t shapeFuncIndex = 0; shapeFuncIndex < shape(); shapeFuncIndex++) {
    // assume local vertex order always being v0, v1 and v2, 
    // and v1 being the vertex with which the shape function is associated with
    size_t v1 = shapeFuncIndex;
    size_t v0 = (shapeFuncIndex - 1 + 3) % 3;
    size_t v2 = (shapeFuncIndex + 1) % 3;
    std::vector<point3D> node(2);
    node[0] = vertex_[v0];
    node[0].transferGlobalToLocalCoord(vertex_[v1], tangent_1, 
				       tangent_2, normal_);
    if (abs(node[0].z()) >= 10*DBL_EPSILON) {
      cout << "node[0].z = " << node[0].z() << endl;
      pfft::errorMessage("element::compShapeFunction",
			 "z should be zero, There must be a bug here!");
    }

    node[1] = vertex_[v2];
    node[1].transferGlobalToLocalCoord(vertex_[v1], tangent_1, 
				       tangent_2, normal_);
    if (abs(node[1].z()) >= 10*DBL_EPSILON) {
      cout << "node[1].z = " << node[1].z() << endl;
      pfft::errorMessage("element::compShapeFunction",
			 "z should be zero, There must be a bug here!");
    }

    double delta = node[0].x() * node[1].y() - node[0].y() * node[1].x();
    linearShapeFuncList_[shapeFuncIndex].a = (node[0].y() - node[1].y()) / delta;
    linearShapeFuncList_[shapeFuncIndex].b = (node[1].x() - node[0].x()) / delta;
    linearShapeFuncList_[shapeFuncIndex].c = 1;
  }
}

/**********************************************************************
 * comparison --
 **********************************************************************/
bool
pfft::operator == (
		   const element& e1, 
		   const element& e2) 
{
  return (e1.shape() == e2.shape()) && (e1.vertex_ == e2.vertex_);
}

bool
pfft::operator != (
		   const element& e1, 
		   const element& e2) 
{
  return (e1.shape() != e2.shape()) || (e1.vertex_ != e2.vertex_);
}

/**********************************************************************
 * output --
 **********************************************************************/
ostream& 
pfft::operator << (
		   ostream& os, 
		   const element& e) 
{
  os << endl 
     << "name = " << e.name() << endl
     << "boundaryIndex = " << e.boundaryIndex() << endl
     << "shape = " << e.shape() << endl;

  for (size_t i=0; i<e.shape(); i++) 
    os << "vertex[" << i << "] = " << e.vertex(i); 

  os << "area = " << e.area() << endl
     << "centroid = " << e.centroid()
     << "tangent1 = " << e.tangent1()
     << "tangent2 = " << e.tangent2()
     << "normal = " << e.normal()
     << "boundingSphereRadius = " << e.boundingSphereRadius()
     << endl;
  return os;
}

/**********************************************************************
 * input --
 **********************************************************************/
istream&
pfft::operator >> (
		   istream& is,
		   element& e)
{
  string shapeKey;
  int shape;

  if (is) {
    is >> shapeKey >> e.boundaryIndex_;
    if ( (shapeKey=="3") || (shapeKey=="t") || (shapeKey=="T") ) 
      shape = 3;
    else if ( (shapeKey=="4") || (shapeKey=="q") || (shapeKey=="Q") ) 
      shape = 4;
    else
      throw domain_error("corrupted element shape data in element line");

    e.vertex_.clear();
    e.vertex_.reserve(shape);
    for (int i=0; i<shape; i++) {
      vector3D<double> vt;
      if ( !(is >> vt) ) 
	throw domain_error("corrupted vertex data in element line");
      e.vertex_.push_back(vt);
    }
    e.setupElement();
  }
  return is;
}

/**********************************************************************
 * readElementFile --
 **********************************************************************/
void
pfft::readElementFile (
		       const char* fileName,
		       vector<element>& elementList)
{
  string title;
  size_t numOfConductor;
  readElementFile(fileName, elementList, title, numOfConductor);    
}

/**********************************************************************
 * readElementFile --
 **********************************************************************/
bool
pfft::readElementFile (
		       const char* fileName,
		       vector<element>& elementList,
		       string& title,
		       size_t& numOfConductor)
{
  numOfConductor = 0;

  element tmp;
  string s;

  ifstream dataSource(fileName);
  if (!dataSource) {
    errorMessage("element::readElementFile",
		 "fail to open the element file");
  }

  int elementIndex=0;
  char c;
  while (dataSource >> c) {
    switch (c) {
    case '0':
      getline(dataSource, title);
      break;
    case '3':
    case 't':
    case 'T':
    case '4':
    case 'q':
    case 'Q':
      dataSource.unget();
      dataSource >> tmp;
      tmp.setName(elementIndex++);
      numOfConductor = max(tmp.boundaryIndex(), numOfConductor);
      elementList.push_back(tmp);
      break;
    default:
      getline(dataSource, s);
    }
  }

  return true;
}

/**********************************************************************
 * findMinEdgeLength --
 **********************************************************************/
double
element::findMinEdgeLength (
			    void) const
{
  double minEdgeLength = 1e20;
  for (size_t i=0; i<vertex_.size(); i++) {
    size_t j = (i+1) % vertex_.size();
    double edgeLength = length(vertex_[i] - vertex_[j]);
    minEdgeLength = minEdgeLength < edgeLength ? minEdgeLength : edgeLength;
  }
  return minEdgeLength;
}

/**********************************************************************
 * shiftToLocalCoord --
 * Use tangent_1, tangent_2 and normal_ as X, Y and Z of the local coord.
 **********************************************************************/
void element::shiftToLocalCoord (
				 const point3D& localOrigin)
{
  for(size_t ii=0; ii<shape(); ii++){
    vertex_[ii].transferGlobalToLocalCoord(localOrigin, tangent_1, tangent_2, normal_);
  }
}

