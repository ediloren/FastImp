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

  const static char cvsid[] = "$Id: gridIndex.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _GRID_INDEX_H_
#define _GRID_INDEX_H_

#include <iostream>

namespace pfft {

  class GridIndex {  
  public:
    GridIndex(void):_x(0),_y(0),_z(0){ }
    GridIndex(const int x, const int y, const int z):_x(x),_y(y),_z(z) {}

    const int x(void) const { return _x; }
    const int y(void) const { return _y; }
    const int z(void) const { return _z; }
    GridIndex& operator -= (const GridIndex& g1);
    GridIndex& operator += (const GridIndex& g1);
    double length(void)  const;
    int length2(void)  const;

  private:
    int _x;
    int _y;
    int _z;
   };

  GridIndex operator - (const GridIndex& g1, const GridIndex& g2);
  GridIndex operator + (const GridIndex& g1, const GridIndex& g2);
  double distance(const GridIndex& g1, const GridIndex& g2);
  bool operator < (const GridIndex& g1, const GridIndex& g2);
  bool operator == (const GridIndex& g1, const GridIndex& g2);
  bool operator != (const GridIndex& g1, const GridIndex& g2);
  std::ostream& operator << (std::ostream& os, const GridIndex& p);

} //namespace pfft

#endif



