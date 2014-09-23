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

  const static char cvsid[] = "$Id: cubature.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _CUBATURE_H_
#define _CUBATURE_H_

#include "vector3D.h"
#include <vector>

namespace pfft {

  class Cubature {

  public:
    enum Shape {TRIANGLE, SQUARE};

    Cubature(void) {};
    Cubature(const int degree, const int rule, const Shape shape);
    vector3D<double> point(const std::vector<vector3D<double> >& vertex, 
			   const size_t pointIndex) const;
    double weight(const double area, const size_t pointIndex) const;
    int numPoint(void) const { return numPoint_; }
    int degree(void) const { return degree_; }
    int rule(void) const { return rule_; }
    bool setupComplete(void) const { return (weight_.size() != 0); }

  private:
    std::vector<vector3D<double> > point_;
    std::vector<double> weight_;
    int degree_;
    int rule_;
    Shape shape_;
    int numPoint_;

    // state variable for setupCubature
    size_t currentPointIndex_;

    void setupSquareCubtature(void);
    void setupTriangleCubtature(void);
    void findNumPointInTriangleCubature (void);
    void findNumPointInSquareCubature (void);
    void findPointAndWeightInTriangle(void);
    void findPointAndWeightInSquare(void);
    void SIMPLEX2_FULLY1 (double, double, double);
    void SIMPLEX2_FULLY3 (double, double, double);
    void SIMPLEX2_FULLY6 (double, double, double);
    void SIMPLEX2_ROTATIONAL3 (double, double, double);

  };

} // namespace pfft

#endif




