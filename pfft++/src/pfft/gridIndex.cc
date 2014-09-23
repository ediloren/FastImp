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
#include "gridIndex.h"
#include <cmath>

using namespace std;
using namespace pfft;

/**********************************************************************
 * length --
 **********************************************************************/
double 
GridIndex::length(
		  void) const
{
  return sqrt( static_cast<double>(_x * _x + _y * _y + _z * _z));
}

/**********************************************************************
 * length2 --
 **********************************************************************/
int 
GridIndex::length2(
		   void) const
{
  return _x * _x + _y * _y + _z * _z;
}

/**********************************************************************
 * operator -= --
 **********************************************************************/
GridIndex& GridIndex::operator -= (const GridIndex& g1){
  _x -= g1._x;
  _y -= g1._y;
  _z -= g1._z;
  return *this;
}

/**********************************************************************
 * operator += --
 **********************************************************************/

GridIndex& GridIndex::operator += (const GridIndex& g1){
  _x += g1._x;
  _y += g1._y;
  _z += g1._z;
  return *this;
}

/**********************************************************************
 * operator - --
 **********************************************************************/
GridIndex
pfft::operator - (
		  const GridIndex& g1,	    
		  const GridIndex& g2) 
{
  GridIndex gt = g1;
  gt -= g2;
  return gt;
}

/**********************************************************************
 * operator + --
 **********************************************************************/
GridIndex
pfft::operator + (
		  const GridIndex& g1,	    
		  const GridIndex& g2) 
{
  GridIndex gt = g1;
  gt += g2;
  return gt;
}

/**********************************************************************
 * distance --
 **********************************************************************/
double 
pfft::distance(
	       const GridIndex& g1, 
	       const GridIndex& g2)
{
  GridIndex gt = g1 - g2;
  return gt.length();
}

/**********************************************************************
 * operator == --
 **********************************************************************/
bool 
pfft::operator == (
		   const GridIndex& g1, 
		   const GridIndex& g2)
{
  return (g1.x()==g2.x()) &&  (g1.y()==g2.y()) && (g1.z()==g2.z());
}

/**********************************************************************
 * operator != --
 **********************************************************************/
bool 
pfft::operator != (
		   const GridIndex& g1, 
		   const GridIndex& g2)
{
  return (g1.x()!=g2.x()) ||  (g1.y()!=g2.y()) || (g1.z()!=g2.z());
}


/**********************************************************************
 * operator < --
 **********************************************************************/
bool 
pfft::operator < (
		  const GridIndex& g1, 
		  const GridIndex& g2)
{  
  if (g1.x() < g2.x()) {
    return true;
  } else if (g1.x() == g2.x()) {
    if (g1.y() < g2.y()) {
      return true;
    } else if (g1.y() == g2.y()) {
      if (g1.z() < g2.z()) {  
	return true;
      } else {
	return false;
      }
    } else {
      return false;
    }
  } else {
    return false;
  }
}

/**********************************************************************
 * output --
 **********************************************************************/
std::ostream& 
operator << (
	     std::ostream& os, 
	     const GridIndex& p)
{
  os <<  "(" << p.x() 
     << ", " << p.y()
     << ", " << p.z() 
     << ")";
  return os;
}



